"""Avalia entradas cristalográficas candidatas (COD ou CSD) contra o orgânico procurado.

Três evidências independentes, em ordem crescente de força:
  fórmula  -> a composição comporta o orgânico (matching.match_formula)
  nome     -> algum nome/variante do orgânico aparece no nome químico ou no título
  grafo    -> o CIF contém uma molécula com a mesma conectividade (cif_analysis)
"""
from __future__ import annotations

import logging

from .chem import A_SITE_METALS, HALOGENS, METALS
from .cif_analysis import analyze_cif
from .matching import composition_key, match_formula
from .names import find_terms

log = logging.getLogger("perovskite_explorer")

CONFIRMADA = "confirmada (grafo no CIF)"
PROVAVEL = "provável (fórmula + nome)"
POSSIVEL = "possível (só fórmula; isômero não descartado)"
POSSIVEL_NOME = "possível (só nome)"
DESCARTADA = "descartada"
ORDER = {CONFIRMADA: 0, PROVAVEL: 1, POSSIVEL: 2, POSSIVEL_NOME: 3, DESCARTADA: 9}


def _category(formula: dict, organic_formula: dict, k: int) -> str:
    metals = [e for e in formula if e in METALS and e not in A_SITE_METALS]
    free_x = [e for e in formula if e in HALOGENS
              and formula[e] - k * organic_formula.get(e, 0) > 0.2]
    if metals and free_x:
        return "haleto metálico híbrido"
    if free_x:
        return "sal de haleto sem metal"
    return "outro"


def prefilter(entries: list, organic: dict, terms: dict) -> list:
    """Primeiro filtro (sem baixar CIF): fórmula e/ou nome."""
    all_terms = terms.get("cations", []) + terms.get("neutros", [])
    resolved = bool(organic.get("resolvido"))
    out = []
    for e in entries:
        if e.get("status") in ("retracted", "errors") or e.get("duplicata_de"):
            continue
        m = match_formula(e["formula"], organic) if resolved else None
        hits_name = find_terms(e.get("nome", ""), all_terms)
        hits_title = find_terms(e.get("titulo", ""), all_terms)
        if not m and not hits_name and not hits_title:
            continue
        k = m["k"] if m else 0
        cat = _category(e["formula"], organic.get("formula") or {}, k)
        if cat == "outro":
            continue  # sem haleto iônico: fora do escopo (ex.: cristal do orgânico puro)
        e = dict(e)
        e.update({"casamento_formula": m, "termos_no_nome": hits_name, "termos_no_titulo": hits_title,
                  "categoria": cat, "composto": composition_key(e["formula"]),
                  "metais": sorted(x for x in e["formula"] if x in METALS),
                  "haletos": sorted(x for x in e["formula"] if x in HALOGENS)})
        out.append(e)
    # prioridade para verificação por CIF: híbridos metálicos, fórmula exata, com nome
    out.sort(key=lambda e: (
        e["categoria"] != "haleto metálico híbrido",
        0 if (e["casamento_formula"] or {}).get("nivel") == "exata" else 1,
        not (e["termos_no_nome"] or e["termos_no_titulo"]),
        -(e.get("ano") or 0),
    ))
    return out


def verify(entries: list, organic: dict, get_cif, max_cifs: int = 60) -> list:
    """Baixa os CIFs (até `max_cifs`) e decide o nível de confiança de cada entrada."""
    query = None
    if organic.get("heavy_elements"):
        query = {"heavy_elements": organic["heavy_elements"], "heavy_edges": organic["heavy_edges"]}
    done = 0
    for e in entries:
        analysis = None
        if done < max_cifs:
            try:
                cif = get_cif(e["id"])
                analysis = analyze_cif(cif, query) if cif else {"ok": False, "erro": "CIF não disponível"}
            except Exception as exc:  # noqa: BLE001 - falha de rede em um CIF não derruba o resto
                analysis = {"ok": False, "erro": f"download falhou: {exc}"}
            done += 1
        else:
            analysis = {"ok": False, "erro": f"não verificado (limite de {max_cifs} CIFs por orgânico)"}
        e["analise"] = analysis
        e["confianca"], e["motivo_confianca"] = decide(e, analysis, has_graph=query is not None)
        if analysis.get("ok"):
            e["classificacao"] = analysis.get("classificacao")
            if not e.get("grupo_espacial"):
                e["grupo_espacial"] = analysis.get("grupo_espacial", "")
        else:
            e["classificacao"] = "não analisada (" + str(analysis.get("erro", "")) + ")"
    return entries


def decide(e: dict, analysis: dict, has_graph: bool = True):
    m = e.get("casamento_formula")
    named = bool(e.get("termos_no_nome"))
    titled = bool(e.get("termos_no_titulo"))
    status = (analysis.get("organico") or {}).get("status") if analysis.get("ok") else None
    if status == "confirmado":
        return CONFIRMADA, "molécula com a mesma conectividade encontrada no CIF"
    if status in ("isomero", "ausente"):
        why = (analysis.get("organico") or {}).get("motivo", "")
        if named:
            return POSSIVEL_NOME, ("o nome químico cita o orgânico, mas o grafo do CIF não confirmou "
                                   f"({why}) — revisar manualmente")
        return DESCARTADA, why
    why = (analysis.get("organico") or {}).get("motivo") or analysis.get("erro") or \
        ("orgânico sem grafo de referência" if not has_graph else "grafo não avaliado")
    if m and (named or titled):
        return PROVAVEL, f"fórmula compatível e nome citado; grafo não conclusivo ({why})"
    if m and m.get("nivel") == "exata":
        return POSSIVEL, f"fórmula compatível; grafo não conclusivo ({why})"
    if m:
        return POSSIVEL, f"fórmula compatível só admitindo solvente/co-cátion; grafo não conclusivo ({why})"
    if named:
        return POSSIVEL_NOME, f"apenas o nome químico cita o orgânico; grafo não conclusivo ({why})"
    return DESCARTADA, "orgânico citado apenas no título do artigo, sem fórmula compatível nem confirmação no CIF"


def group_compounds(structures: list) -> list:
    """Agrupa depósitos do mesmo composto: é aqui que se mede 'o quanto já foi reportado'."""
    groups: dict = {}
    for s in structures:
        if s["confianca"] == DESCARTADA:
            continue
        g = groups.setdefault(s["composto"], {
            "composto": s["composto"], "categoria": s["categoria"], "metais": s["metais"],
            "haletos": s["haletos"], "entradas": [], "dois": set(), "anos": [],
            "grupos_espaciais": set(), "classificacoes": set(), "melhor_confianca": s["confianca"],
        })
        g["entradas"].append(f"{s['fonte']}:{s['id']}")
        if s.get("doi"):
            g["dois"].add(s["doi"])
        if s.get("ano"):
            g["anos"].append(s["ano"])
        if s.get("grupo_espacial"):
            g["grupos_espaciais"].add(s["grupo_espacial"].replace(" ", ""))
        if s.get("classificacao") and (s.get("analise") or {}).get("ok"):
            g["classificacoes"].add(s["classificacao"])
        if ORDER[s["confianca"]] < ORDER[g["melhor_confianca"]]:
            g["melhor_confianca"] = s["confianca"]
    out = []
    for g in groups.values():
        g["n_depositos"] = len(g["entradas"])
        g["n_publicacoes"] = len(g["dois"])
        g["dois"] = sorted(g["dois"])
        g["primeiro_ano"] = min(g["anos"]) if g["anos"] else None
        g["ultimo_ano"] = max(g["anos"]) if g["anos"] else None
        g["grupos_espaciais"] = sorted(g["grupos_espaciais"])
        g["classificacoes"] = sorted(g["classificacoes"])
        del g["anos"]
        out.append(g)
    out.sort(key=lambda g: (g["categoria"] != "haleto metálico híbrido", ORDER[g["melhor_confianca"]],
                            -g["n_depositos"]))
    return out
