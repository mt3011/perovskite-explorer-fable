"""Orquestração: para cada orgânico, consulta as fontes, grava um checkpoint e atualiza os relatórios.

Proteções contra interrupção:
  * cada resposta HTTP fica no cache SQLite (http.py);
  * cada orgânico concluído vira um arquivo JSON (escrita atômica) em <saida>/checkpoints;
  * os relatórios são regravados após cada orgânico (sempre há um resultado parcial utilizável);
  * falha de uma fonte não interrompe as outras e fica registrada (o orgânico é marcado como
    "busca incompleta" e refeito na próxima execução).
"""
from __future__ import annotations

import datetime as dt
import json
import logging
import os
import traceback
from pathlib import Path

from . import __version__, cod, csd, datacite, literature, pubchem, scoring, structures
from .chem import DEFAULT_HALOGENS, DEFAULT_METALS
from .http import SourceError
from .names import build_terms, find_terms_in_chemname, slugify

log = logging.getLogger("perovskite_explorer")

SCHEMA = 2   # muda quando o formato do checkpoint muda; checkpoints antigos são refeitos

DEFAULT_CONFIG = {
    "metais": DEFAULT_METALS,
    "halogenios": DEFAULT_HALOGENS,
    "fontes_literatura": ["openalex", "europepmc", "crossref"],
    "usar_cod": True,
    "usar_csd": True,          # só tem efeito se a CSD Python API estiver instalada
    "usar_datacite": True,     # CSD por nome e por artigo, via metadados públicos do DataCite
    "max_cifs": 300,           # limite de CIFs para candidatas só por fórmula com orgânico neutro (B2)
    "max_paginas": 2,          # páginas de 100 resultados por fonte de literatura
    "max_doi_datacite": 25,    # artigos estruturais cujos depósitos CCDC são listados
    "max_doi_crossref": 80,    # artigos de estruturas cujos metadados/citações são buscados
    "email": None,
    "openalex_api_key": None,
    "openalex_sem_chave": False,   # True = tenta o OpenAlex mesmo sem chave (cota anônima, instável)
    "semanticscholar_api_key": None,
}


def read_list(path) -> list:
    """Lê a lista de orgânicos (um por linha; '#' inicia comentário; duplicatas são removidas)."""
    names, seen = [], set()
    text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip().strip(";,")
        if line and line.lower() not in seen:
            seen.add(line.lower())
            names.append(line)
    return names


def atomic_write_json(path: Path, obj) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=1, default=_default)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _default(o):
    if isinstance(o, (set, frozenset)):
        return sorted(o)
    if hasattr(o, "item"):
        return o.item()
    return str(o)


def config_signature(cfg: dict) -> dict:
    """Parâmetros que mudam o resultado: se diferirem do checkpoint, o orgânico é refeito."""
    keys = ("metais", "halogenios", "fontes_literatura", "usar_cod", "usar_csd", "usar_datacite", "max_cifs")
    sig = {k: cfg.get(k) for k in keys}
    sig["openalex_ativo"] = bool(cfg.get("openalex_api_key") or cfg.get("openalex_sem_chave"))
    return sig


def process_organic(http, name: str, cfg: dict) -> dict:
    res = {"entrada": name, "versao": __version__, "schema": SCHEMA,
           "inicio": dt.datetime.now().isoformat(timespec="seconds"),
           "config": config_signature(cfg), "fontes": {}, "estruturas": [], "descartadas": [], "compostos": [], "literatura": [],
           "depositos_csd": []}

    # 1) identidade química -------------------------------------------------
    try:
        org = pubchem.resolve(http, name)
        res["fontes"]["PubChem"] = {"ok": True, "detalhe": org.get("metodo_resolucao", "")}
    except (SourceError, KeyError, IndexError, TypeError) as exc:
        org = {"entrada": name, "resolvido": False, "avisos": [f"PubChem falhou: {exc}"]}
        res["fontes"]["PubChem"] = {"ok": False, "erro": str(exc)}
    res["organico"] = org
    rel = org.get("formas_relacionadas") or {}
    terms = build_terms(name, org.get("titulo", ""), org.get("nome_iupac", ""), org.get("sinonimos"),
                        extra_cations=rel.get("cations"))
    res["termos"] = terms

    # 2) bancos cristalográficos -------------------------------------------
    all_structs = []
    res["verificacao"] = {"cifs_analisados": 0, "hibridos_nao_verificados": 0, "sais_omitidos": 0}

    def add_stats(stats):
        for k, v in stats.items():
            res["verificacao"][k] += v

    if cfg["usar_cod"]:
        try:
            entries, rep = cod.search(http, org, terms, cfg["metais"], cfg["halogenios"])
            cands = structures.prefilter(entries, org, terms, cfg["metais"], cfg["halogenios"])
            cands, stats = structures.verify(cands, org, lambda i: cod.fetch_cif(http, i), cfg["max_cifs"])
            add_stats(stats)
            all_structs += cands
            res["fontes"]["COD"] = {"ok": True, "consultas": rep["consultas"],
                                    "entradas_examinadas": len(entries), "candidatas": len(cands),
                                    "cifs_analisados": stats["cifs_analisados"],
                                    "falhas_parciais": rep["falhas"][:5]}
        except SourceError as exc:
            res["fontes"]["COD"] = {"ok": False, "erro": str(exc)}
    if cfg["usar_csd"] and csd.available():
        if org.get("heavy_elements"):
            try:
                entries, cifs = csd.search(org)
                cands = structures.prefilter(entries, org, terms, cfg["metais"], cfg["halogenios"])
                cands, stats = structures.verify(cands, org, cifs.get, max_cifs=10 ** 6, max_salts=10 ** 6)
                add_stats(stats)
                all_structs += cands
                res["fontes"]["CSD (API)"] = {"ok": True, "entradas_examinadas": len(entries), "candidatas": len(cands)}
            except Exception as exc:  # noqa: BLE001 - conector ainda não validado com licença real
                res["fontes"]["CSD (API)"] = {"ok": False, "erro": f"{type(exc).__name__}: {exc}"}
    if cfg["usar_datacite"]:
        try:
            found = datacite.search_names(http, terms, cfg["halogenios"])
            have = {s.get("ccdc") for s in all_structs if s.get("ccdc")}
            csd_structs = [x for x in (datacite.to_structure(r, terms, cfg["halogenios"], cfg["metais"])
                                       for r in found["registros"]) if x and x["ccdc"] not in have]
            all_structs += csd_structs
            res["csd_total_qualquer_composto"] = found["total_qualquer"]
            res["fontes"]["CSD (nomes via DataCite)"] = {
                "ok": True, "depositos_com_haleto_examinados": len(found["registros"]),
                "aceitos_pelo_nome": len(csd_structs), "depositos_citando_o_organico": found["total_qualquer"]}
            if found["truncado"]:
                res["fontes"]["CSD (nomes via DataCite)"]["falhas_parciais"] = [
                    f"lista truncada: {found['total_com_haleto']} depósitos, lidos {len(found['registros'])}"]
        except SourceError as exc:
            res["fontes"]["CSD (nomes via DataCite)"] = {"ok": False, "erro": str(exc)}
    for s in all_structs:
        (res["descartadas"] if s["confianca"] == structures.DESCARTADA else res["estruturas"]).append(s)
    res["estruturas"].sort(key=lambda s: (s["categoria"] != structures.HYBRID, structures.ORDER[s["confianca"]],
                                          -(s.get("ano") or 0)))
    res["compostos"] = structures.group_compounds(res["estruturas"])

    # 3) literatura ----------------------------------------------------------
    sources = list(cfg["fontes_literatura"])
    if "openalex" in sources and not cfg.get("openalex_api_key") and not cfg.get("openalex_sem_chave"):
        # sem chave, a cota anônima do OpenAlex é compartilhada por IP e acaba sem aviso: resultados
        # variariam de um dia para o outro. Preferimos não usar a fonte e dizer isso no relatório.
        sources.remove("openalex")
        res["fontes"]["OpenAlex"] = {"ok": True, "desativada": "sem chave de API (gratuita): defina "
                                     "OPENALEX_API_KEY para ativar, ou use --openalex-sem-chave"}
    papers, rep = literature.run(http, terms, cfg, sources)
    for label, info in rep.items():
        res["fontes"][label] = info if isinstance(info, dict) else {"ok": True, "detalhe": str(info)}

    # artigos das estruturas encontradas entram na lista (com metadados do Crossref)
    by_doi = {p["doi"]: p for p in papers if p.get("doi")}
    doi_structs: dict = {}
    for s in res["estruturas"]:
        if s.get("doi") and s["confianca"] in structures.FORTES and s["categoria"] == scoring.HYBRID:
            doi_structs.setdefault(s["doi"], []).append(f"{s['fonte']}:{s['id']}")
    lookup_fail = n_lookup = 0
    for doi, ids in doi_structs.items():
        p = by_doi.get(doi)
        if p is None:
            p = None
            if n_lookup < cfg["max_doi_crossref"]:    # limite para orgânicos com centenas de estruturas
                n_lookup += 1
                try:
                    p = literature.crossref_lookup(http, doi, cfg)
                except SourceError as exc:
                    lookup_fail += 1
                    log.warning("Crossref (DOI %s): %s", doi, exc)
            if p is None:
                s0 = next(s for s in res["estruturas"] if s.get("doi") == doi)
                p = {"doi": doi, "titulo": s0.get("titulo", ""), "resumo": "", "ano": s0.get("ano"),
                     "revista": s0.get("revista", ""), "autores": [a.strip() for a in s0.get("autores", "").split(";") if a.strip()],
                     "citacoes": None, "tipo": "article", "retratado": False, "fontes": [s0["fonte"]],
                     "url": f"https://doi.org/{doi}"}
            p = literature.classify(p, terms)
            papers.append(p)
            by_doi[doi] = p
        p["estruturas_ligadas"] = ids
        p["relevante"] = True
        if p["uso"] != literature.ESTRUTURAL:
            p["uso_pelo_texto"] = p["uso"]
            p["uso"] = literature.ESTRUTURAL
            p.setdefault("pistas_estruturais", []).append("estrutura cristalina depositada com este DOI")
    if lookup_fail:
        res["fontes"].setdefault("Crossref", {"ok": True})["falhas_parciais"] = [f"{lookup_fail} DOI(s) sem metadados"]

    # 4) depósitos CSD de cada artigo estrutural (via DataCite), com o nome químico de cada um
    if cfg["usar_datacite"]:
        targets = [p for p in papers if p.get("doi") and p.get("relevante")
                   and (p["uso"] == literature.ESTRUTURAL)]
        targets.sort(key=lambda p: (not p.get("estruturas_ligadas"), -(p.get("citacoes") or 0)))
        all_terms = terms.get("cations", []) + terms.get("neutros", [])
        n_fail, err = 0, ""
        for p in targets[: cfg["max_doi_datacite"]]:
            try:
                deps = datacite.depositions_for_doi(http, p["doi"])
            except SourceError as exc:
                n_fail += 1
                err = str(exc)
                if n_fail >= 3:
                    break
                continue
            for d in deps:
                d["nome_cita_o_organico"] = bool(find_terms_in_chemname(d["nome"], all_terms))
            p["depositos_csd"] = [d["ccdc"] for d in deps if d.get("ccdc")]
            p["depositos_csd_com_o_organico"] = [d["ccdc"] for d in deps if d["nome_cita_o_organico"] and d.get("ccdc")]
            res["depositos_csd"] += deps
        if targets:
            ok = n_fail == 0
            res["fontes"]["CSD (depósitos por artigo, DataCite)"] = {
                "ok": ok, "artigos_consultados": min(len(targets), cfg["max_doi_datacite"])}
            if not ok:
                res["fontes"]["CSD (depósitos por artigo, DataCite)"]["erro"] = err

    order = {literature.ESTRUTURAL: 0, literature.AMBIGUO: 1, literature.ADITIVO: 2, literature.MENCAO: 3,
             literature.NAO_VERIFICAVEL: 4}
    papers.sort(key=lambda p: (order.get(p.get("uso"), 9), -(p.get("citacoes") or 0)))
    res["literatura"] = papers
    res["resumo"] = scoring.summarize(res)
    res["completo"] = not any(not v.get("ok") or v.get("falhas_parciais") for v in res["fontes"].values())
    res["fim"] = dt.datetime.now().isoformat(timespec="seconds")
    return res


def load_results(out_dir: Path, names: list | None = None) -> list:
    """Carrega os checkpoints existentes (na ordem da lista, se fornecida)."""
    ck = Path(out_dir) / "checkpoints"
    found = {}
    for f in sorted(ck.glob("*.json")):
        try:
            r = json.loads(f.read_text(encoding="utf-8"))
            if r.get("schema") != SCHEMA:
                log.warning("checkpoint de versão antiga ignorado (será refeito pelo comando buscar): %s", f.name)
                continue
            # o resumo é sempre recalculado: assim `relatorio` reflete as regras da versão atual
            r["resumo"] = scoring.summarize(r)
            found[r["entrada"].lower()] = r
        except (ValueError, KeyError):
            log.warning("checkpoint ilegível ignorado: %s", f.name)
    if names is None:
        return list(found.values())
    return [found[n.lower()] for n in names if n.lower() in found]


def run(names: list, out_dir, http, cfg: dict, redo: bool = False, retry_incomplete: bool = True,
        on_progress=None) -> list:
    from . import report

    out_dir = Path(out_dir)
    ck = out_dir / "checkpoints"
    ck.mkdir(parents=True, exist_ok=True)
    state_file = out_dir / "estado.json"
    total = len(names)
    interrupted = False
    for i, name in enumerate(names, 1):
        path = ck / f"{slugify(name)}.json"
        if path.exists() and not redo:
            try:
                prev = json.loads(path.read_text(encoding="utf-8"))
                same_cfg = (prev.get("config") == json.loads(json.dumps(config_signature(cfg)))
                            and prev.get("schema") == SCHEMA)
                if same_cfg and (prev.get("completo") or not retry_incomplete):
                    log.info("[%d/%d] %s — já concluído (checkpoint), pulando", i, total, name)
                    continue
                log.info("[%d/%d] %s — %s", i, total, name,
                         "checkpoint incompleto, refazendo o que falhou" if same_cfg
                         else "parâmetros mudaram desde o checkpoint, refazendo")
            except ValueError:
                log.warning("[%d/%d] %s — checkpoint corrompido, refazendo", i, total, name)
        else:
            log.info("[%d/%d] %s — iniciando", i, total, name)
        atomic_write_json(state_file, {"total": total, "posicao": i, "atual": name,
                                       "atualizado": dt.datetime.now().isoformat(timespec="seconds")})
        try:
            res = process_organic(http, name, cfg)
        except KeyboardInterrupt:
            interrupted = True
            log.warning("Interrompido pelo usuário. Rode o mesmo comando para continuar de onde parou.")
            break
        except Exception as exc:  # noqa: BLE001 - um orgânico problemático não derruba a lista inteira
            log.error("[%d/%d] %s — erro interno: %s\n%s", i, total, name, exc, traceback.format_exc())
            res = {"entrada": name, "versao": __version__, "schema": SCHEMA, "completo": False, "erro_interno": f"{type(exc).__name__}: {exc}",
                   "organico": {"entrada": name, "resolvido": False, "avisos": []}, "termos": {},
                   "fontes": {"programa": {"ok": False, "erro": str(exc)}}, "estruturas": [], "descartadas": [],
                   "compostos": [], "literatura": [], "depositos_csd": []}
            res["resumo"] = scoring.summarize(res)
        atomic_write_json(path, res)
        s = res["resumo"]
        log.info("[%d/%d] %s — %s (índice %s)%s", i, total, name, s["nivel_exploracao"], s["indice_oportunidade"],
                 "" if s["busca_completa"] else " — BUSCA INCOMPLETA")
        try:
            report.write_all(load_results(out_dir, names), out_dir, cfg, total=total)
        except Exception as exc:  # noqa: BLE001 - relatório parcial é conveniência; o checkpoint já está salvo
            log.warning("não foi possível atualizar os relatórios parciais: %s", exc)
        if on_progress:
            on_progress(i, total, res)
    results = load_results(out_dir, names)
    report.write_all(results, out_dir, cfg, total=total)
    atomic_write_json(state_file, {"total": total, "concluidos": len(results), "interrompido": interrupted,
                                   "atualizado": dt.datetime.now().isoformat(timespec="seconds")})
    return results
