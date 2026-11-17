"""Avalia entradas cristalográficas candidatas (COD ou CSD com CIF) contra o orgânico procurado.

Três evidências independentes, em ordem crescente de força:
  fórmula  -> a composição comporta o orgânico (matching.match_formula)
  nome     -> algum nome/variante do orgânico aparece no nome químico ou no título
  grafo    -> o CIF contém uma molécula com a mesma conectividade (cif_analysis)

Como isômeros têm a mesma fórmula (ex.: fenazina e 1,10-fenantrolina), nenhuma entrada é
contada como evidência só pela fórmula: ela precisa do grafo (CIF) ou do nome.
"""
from __future__ import annotations

import logging

from .chem import A_SITE_METALS, HALOGENS, METALS
from .cif_analysis import analyze_cif
from .matching import composition_key, match_formula
from .names import find_terms, find_terms_in_chemname

log = logging.getLogger("perovskite_explorer")

CONFIRMADA = "confirmada (grafo no CIF)"
PROVAVEL = "provável (fórmula + nome)"
NOME_CSD = "provável (nome sistemático no CSD; sem CIF)"
POSSIVEL = "possível (só fórmula; isômero não descartado)"
POSSIVEL_NOME = "possível (só nome)"
DESCARTADA = "descartada"
ORDER = {CONFIRMADA: 0, PROVAVEL: 1, NOME_CSD: 1, POSSIVEL: 2, POSSIVEL_NOME: 3, DESCARTADA: 9}
FORTES = (CONFIRMADA, PROVAVEL, NOME_CSD)
HYBRID = "haleto metálico híbrido"
SALT = "sal de haleto sem metal"

CATION = "cátion (protonado)"
LIGAND = "ligante coordenado ao metal"
NEUTRAL = "neutro (co-cristal/solvato) ou H não localizado"
UNKNOWN = "indeterminado"


def _category(formula: dict, organic_formula: dict, k: int) -> str:
    metals = [e for e in formula if e in METALS and e not in A_SITE_METALS]
    free_x = [e for e in formula if e in HALOGENS
              and formula[e] - k * organic_formula.get(e, 0) > 0.2]
    if metals and free_x:
        return HYBRID
    if free_x:
        return SALT
    return "outro"


def prefilter(entries: list, organic: dict, terms: dict) -> list:
    """Primeiro filtro (sem baixar CIF). Atribui a cada candidata uma prioridade de verificação:

      A  nome do orgânico no nome químico ou no título do artigo
      B1 fórmula exata (só metais + haletos + água além do orgânico) com orgânico protonado
      B2 fórmula exata com orgânico neutro ou H não informado
    Fórmulas que só fecham admitindo solvente/co-cátion e sem nome citado são descartadas aqui.
    """
    all_terms = terms.get("cations", []) + terms.get("neutros", [])
    resolved = bool(organic.get("resolvido"))
    out = []
    for e in entries:
        if e.get("status") in ("retracted", "errors") or e.get("duplicata_de"):
            continue
        m = match_formula(e["formula"], organic) if resolved else None
        hits_name = find_terms_in_chemname(e.get("nome", ""), all_terms)
        hits_title = find_terms(e.get("titulo", ""), all_terms)
        named = bool(hits_name or hits_title)
        if not named and (not m or m["nivel"] != "exata"):
            continue
        k = m["k"] if m else 0
        cat = _category(e["formula"], organic.get("formula") or {}, k)
        if cat == "outro":
            continue  # sem haleto iônico: fora do escopo (ex.: cristal do orgânico puro)
        tier = "A" if named else ("B1" if m.get("cationico") else "B2")
        e = dict(e)
        e.update({"casamento_formula": m, "termos_no_nome": hits_name, "termos_no_titulo": hits_title,
                  "categoria": cat, "composto": composition_key(e["formula"]), "prioridade": tier,
                  "metais": sorted(x for x in e["formula"] if x in METALS),
                  "haletos": sorted(x for x in e["formula"] if x in HALOGENS)})
        out.append(e)
    rank = {(HYBRID, "A"): 0, (HYBRID, "B1"): 1, (SALT, "A"): 2, (HYBRID, "B2"): 3, (SALT, "B1"): 4, (SALT, "B2"): 5}
    out.sort(key=lambda e: (rank[(e["categoria"], e["prioridade"])], -(e.get("ano") or 0)))
    return out


def verify(entries: list, organic: dict, get_cif, max_cifs: int = 300, max_salts: int = 30):
    """Baixa os CIFs e decide o nível de confiança de cada entrada.

    Sempre verificadas: prioridade A e híbridos B1. Híbridos B2: até `max_cifs`. Sais sem metal
    identificados só por fórmula: até `max_salts`; os excedentes são omitidos (não são o foco).
    Retorna (entradas_avaliadas, estatísticas).
    """
    query = None
    if organic.get("heavy_elements"):
        query = {"heavy_elements": organic["heavy_elements"], "heavy_edges": organic["heavy_edges"]}
    stats = {"cifs_analisados": 0, "hibridos_nao_verificados": 0, "sais_omitidos": 0}
    n_b2 = n_salt = 0
    out = []
    for e in entries:
        tier, cat = e["prioridade"], e["categoria"]
        check = True
        if tier != "A":
            if cat == HYBRID and tier == "B2":
                n_b2 += 1
                check = n_b2 <= max_cifs
            elif cat == SALT:
                n_salt += 1
                if n_salt > max_salts:
                    stats["sais_omitidos"] += 1
                    continue
        if check:
            try:
                cif = get_cif(e["id"])
                analysis = analyze_cif(cif, query) if cif else {"ok": False, "erro": "CIF não disponível"}
            except Exception as exc:  # noqa: BLE001 - falha de rede em um CIF não derruba o resto
                analysis = {"ok": False, "erro": f"download falhou: {exc}"}
            stats["cifs_analisados"] += 1
        else:
            analysis = {"ok": False, "erro": f"não verificada: limite de {max_cifs} CIFs (aumente --max-cifs)"}
            stats["hibridos_nao_verificados"] += 1
        e["analise"] = analysis
        e["verificada"] = check
        e["confianca"], e["motivo_confianca"] = decide(e, analysis, has_graph=query is not None)
        e["papel_do_organico"] = role(e, analysis, organic)
        if analysis.get("ok"):
            e["classificacao"] = analysis.get("classificacao")
            if not e.get("grupo_espacial"):
                e["grupo_espacial"] = analysis.get("grupo_espacial", "")
        else:
            e["classificacao"] = "não analisada (" + str(analysis.get("erro", "")) + ")"
        out.append(e)
    return out, stats


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
    if m:
        return POSSIVEL, f"fórmula compatível; grafo não conclusivo ({why})"
    if named:
        return POSSIVEL_NOME, f"apenas o nome químico cita o orgânico; grafo não conclusivo ({why})"
    return DESCARTADA, "orgânico citado apenas no título do artigo, sem fórmula compatível nem confirmação no CIF"


def role(e: dict, analysis: dict, organic: dict) -> str:
    """Papel do orgânico no cristal: cátion, ligante coordenado ou neutro."""
    m = e.get("casamento_formula") or {}
    if analysis.get("ok"):
        if analysis.get("organico_alvo_coordenado_ao_metal"):
            return LIGAND
        found = (analysis.get("organico") or {}).get("h_por_molecula_no_cif") or []
        h_neutral = (organic.get("formula") or {}).get("H", 0)
        if found and max(found) > h_neutral:
            return CATION
    if organic.get("n_cationicos_permanentes"):
        return "cátion permanente"
    if m.get("cationico"):
        return CATION
    if m.get("h_informado") or analysis.get("ok"):
        return NEUTRAL
    return UNKNOWN


def unique_samples(strong: list) -> int:
    """Nº de amostras distintas somando COD e CSD sem contar duas vezes a mesma estrutura.

    Estruturas do mesmo artigo com os mesmos metais e haletos podem estar nos dois bancos;
    para cada (artigo, M, X) conta-se o maior dos dois totais. Sem DOI, cada entrada conta.
    """
    groups: dict = {}
    n = 0
    for s in strong:
        if not s.get("doi"):
            n += 1
            continue
        g = groups.setdefault((s["doi"], tuple(s["metais"]), tuple(s["haletos"])), {})
        g[s["fonte"]] = g.get(s["fonte"], 0) + 1
    return n + sum(max(g.values()) for g in groups.values())


def group_compounds(structures: list) -> list:
    """Agrupa depósitos do mesmo composto: é aqui que se mede 'o quanto já foi reportado'."""
    groups: dict = {}
    for s in structures:
        if s["confianca"] == DESCARTADA:
            continue
        g = groups.setdefault((s["fonte"], s["composto"]), {
            "composto": s["composto"], "fonte": s["fonte"], "categoria": s["categoria"], "metais": s["metais"],
            "haletos": s["haletos"], "entradas": [], "dois": set(), "anos": [], "papeis": set(),
            "grupos_espaciais": set(), "classificacoes": set(), "melhor_confianca": s["confianca"],
        })
        g["entradas"].append(f"{s['fonte']}:{s['id']}")
        if s.get("doi"):
            g["dois"].add(s["doi"])
        if s.get("ano"):
            g["anos"].append(s["ano"])
        if s.get("grupo_espacial"):
            g["grupos_espaciais"].add(s["grupo_espacial"].replace(" ", ""))
        if s.get("classificacao") and ((s.get("analise") or {}).get("ok") or s["fonte"] == "CSD"):
            g["classificacoes"].add(s["classificacao"])
        if s.get("papel_do_organico"):
            g["papeis"].add(s["papel_do_organico"])
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
        g["papeis"] = sorted(g["papeis"])
        del g["anos"]
        out.append(g)
    out.sort(key=lambda g: (g["categoria"] != HYBRID, ORDER[g["melhor_confianca"]], -g["n_depositos"]))
    return out
