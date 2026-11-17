"""Resumo por orgânico: quanto já foi explorado e índice de oportunidade (regras explícitas)."""
from __future__ import annotations

from collections import Counter

from . import literature as lit
from . import structures as st

INEXPLORADO = "INEXPLORADO"
SEM_USO = "SEM USO ESTRUTURAL CONHECIDO"
INDICIOS = "INDÍCIOS NÃO CONFIRMADOS"
LIGANTE = "USADO SÓ COMO LIGANTE OU NEUTRO"
POUCO = "POUCO EXPLORADO"
MODERADO = "MODERADAMENTE EXPLORADO"
ESTABELECIDO = "BEM ESTABELECIDO"
NOVELTY = {INEXPLORADO: 100, SEM_USO: 90, LIGANTE: 80, INDICIOS: 75, POUCO: 65, MODERADO: 35, ESTABELECIDO: 10}
HYBRID = st.HYBRID
SALT = st.SALT
CORE_METALS = ["Pb", "Sn", "Bi", "Sb", "Cu", "Ag"]


def viability(org: dict):
    """Fator 0–1: facilidade de o orgânico virar cátion (requisito usual para entrar na rede)."""
    if not org.get("resolvido"):
        return 0.5, "estrutura química não resolvida (fator neutro 0,5)"
    if org.get("n_cationicos_permanentes", 0) or org.get("n_basicos_fortes", 0):
        return 1.0, "tem amina alifática ou carga permanente (protona facilmente)"
    if org.get("n_basicos_medios", 0):
        return 0.85, "tem N sp2 básico (piridina/imidazol/amidina): protona em meio ácido HX"
    if org.get("n_basicos_fracos", 0):
        return 0.7, "apenas base fraca (tipo anilina): exige meio HX concentrado"
    return 0.3, "sem sítio básico: só como ligante neutro, co-cristal ou via coordenação"


def summarize(res: dict) -> dict:
    org = res.get("organico") or {}
    structs = [s for s in res.get("estruturas", []) if s["confianca"] != st.DESCARTADA]
    hybrid = [s for s in structs if s["categoria"] == HYBRID]
    salts = [s for s in structs if s["categoria"] == SALT]
    strong = [s for s in hybrid if s["confianca"] in st.FORTES]
    # "fracas": verificação tentada e inconclusiva, ou só nome. Candidatas ainda não verificadas
    # (limite de CIFs) são contadas à parte e não entram em nenhuma conclusão.
    weak = [s for s in hybrid if s["confianca"] not in st.FORTES and s.get("verificada", True)]
    pending = [s for s in hybrid if s["confianca"] not in st.FORTES and not s.get("verificada", True)]
    roles = Counter(s.get("papel_do_organico") or st.UNKNOWN for s in strong)

    papers = res.get("literatura", [])
    rel = [p for p in papers if p.get("relevante")]
    structural = [p for p in rel if p["uso"] == lit.ESTRUTURAL or p.get("estruturas_ligadas")]
    additive = [p for p in rel if p["uso"] == lit.ADITIVO and not p.get("estruturas_ligadas")]
    ambiguous = [p for p in rel if p["uso"] == lit.AMBIGUO and not p.get("estruturas_ligadas")]
    mentions = [p for p in rel if p["uso"] == lit.MENCAO and not p.get("estruturas_ligadas")]
    unverifiable = [p for p in papers if p.get("uso") == lit.NAO_VERIFICAVEL]

    strong_dois = {s["doi"] for s in strong if s.get("doi")}
    keys = {p.get("doi") or p.get("titulo") for p in structural} | strong_dois
    n_pubs = len(keys)
    if strong and n_pubs == 0:
        n_pubs = 1  # estrutura depositada sem DOI associado
    groups = lit.independent_groups(structural)
    if n_pubs and not groups:
        groups = 1
    cites = [p["citacoes"] for p in structural if p.get("citacoes") is not None]
    years = [p["ano"] for p in structural if p.get("ano")] + [s["ano"] for s in strong if s.get("ano")]

    # --- separação por papel: para perovskitas interessa o orgânico como CÁTION na rede ---
    # Papel indeterminado conta como cátion (lado conservador: não inflar a "novidade").
    def cationic(x):
        r = x.get("papel_do_organico") or st.UNKNOWN
        return r.startswith("cátion") or r == st.UNKNOWN

    cation_keys = {x.get("doi") or f"{x['fonte']}:{x['id']}" for x in strong if cationic(x)}
    # artigos estruturais sem estrutura localizada nos bancos também contam como uso catiônico
    lit_only = [p for p in structural if not p.get("estruturas_ligadas")]
    cation_keys |= {p.get("doi") or p.get("titulo") for p in lit_only}
    cation_papers = [p for p in structural if (p.get("doi") in cation_keys) or not p.get("estruturas_ligadas")]
    n_cat = len(cation_keys)
    groups_cat = lit.independent_groups(cation_papers) or (1 if n_cat else 0)
    n_lig = len({x.get("doi") or f"{x['fonte']}:{x['id']}" for x in strong if not cationic(x)} - cation_keys)

    if n_cat == 0:
        if n_pubs > 0:
            level = LIGANTE
        elif weak or ambiguous:
            level = INDICIOS
        elif additive or mentions:
            level = SEM_USO
        else:
            level = INEXPLORADO
    elif n_cat <= 2 or groups_cat <= 1:
        level = POUCO
    elif n_cat <= 9:
        level = MODERADO
    else:
        level = ESTABELECIDO

    combos = Counter()
    dims = Counter()
    seen_pairs = set()
    for s in strong:
        for m in s["metais"]:
            for x in s["haletos"]:
                key = (s.get("doi") or s["id"], m, x)
                if key not in seen_pairs:   # um artigo conta uma vez por par M–X
                    seen_pairs.add(key)
                    combos[f"{m}–{x}"] += 1
        a = s.get("analise") or {}
        if a.get("ok") and a.get("dimensionalidade") is not None:
            dims[f"{a['dimensionalidade']}D"] += 1
    reported_pairs = set(combos)
    gaps = [f"{m}–{x}" for m in CORE_METALS for x in ("Cl", "Br", "I") if f"{m}–{x}" not in reported_pairs]

    factor, why = viability(org)
    index = round(NOVELTY[level] * factor)
    failed = [k for k, v in (res.get("fontes") or {}).items() if not v.get("ok")]
    partial = [k for k, v in (res.get("fontes") or {}).items() if v.get("ok") and v.get("falhas_parciais")]

    summary = {
        "nivel_exploracao": level,
        "indice_oportunidade": index,
        "fator_viabilidade": factor,
        "viabilidade": why,
        "busca_completa": not failed and not partial,
        "fontes_com_falha": failed + [f"{p} (parcial)" for p in partial],
        "n_estruturas_hibridas_fortes": st.unique_samples(strong),
        "n_estruturas_hibridas_confirmadas": sum(1 for s in strong if s["confianca"] == st.CONFIRMADA),
        "n_estruturas_cod": sum(1 for s in strong if s["fonte"] == "COD"),
        "n_estruturas_csd": sum(1 for s in strong if s["fonte"] == "CSD"),
        "papeis_do_organico": dict(roles.most_common()),
        "n_como_cation": sum(v for k, v in roles.items() if k.startswith("cátion")),
        "n_estruturas_hibridas_possiveis": len(weak),
        "n_candidatas_nao_verificadas": len(pending),
        "n_compostos_hibridos": len({(s["fonte"], s["composto"]) for s in strong}),
        "n_sais_haleto_sem_metal": st.unique_samples([s for s in salts if s["confianca"] in st.FORTES]),
        "n_csd_qualquer_composto": res.get("csd_total_qualquer_composto"),
        "n_descartadas": len(res.get("descartadas", [])),
        "combinacoes_MX": dict(combos.most_common()),
        "dimensionalidades": dict(dims),
        "lacunas_MX": gaps if (strong or structural) else [],
        "n_publicacoes_como_cation": n_cat,
        "n_grupos_como_cation": groups_cat,
        "n_publicacoes_so_ligante_ou_neutro": n_lig,
        "n_publicacoes_estruturais": n_pubs,
        "n_grupos_independentes": groups,
        "citacoes_estruturais": sum(cites) if cites else 0,
        "primeiro_ano": min(years) if years else None,
        "ultimo_ano": max(years) if years else None,
        "n_artigos_aditivo": len(additive),
        "n_artigos_ambiguos": len(ambiguous),
        "n_mencoes": len(mentions),
        "n_nao_verificaveis": len(unverifiable),
        "n_depositos_csd_ligados": len({d.get("ccdc") or d.get("doi_deposito") for d in res.get("depositos_csd", [])}),
    }
    summary["parecer"] = _verdict(summary, org)
    return summary


def _verdict(s: dict, org: dict) -> str:
    lv = s["nivel_exploracao"]
    parts = []
    if lv == INEXPLORADO:
        parts.append("Nenhuma estrutura cristalina nem artigo sobre haletos metálicos/perovskitas foi encontrado "
                     "com este orgânico.")
    elif lv == SEM_USO:
        parts.append(f"Aparece na literatura de perovskitas apenas como aditivo/passivante ou menção "
                     f"({s['n_artigos_aditivo']} aditivo, {s['n_mencoes']} menções); nenhum uso como parte da "
                     "rede cristalina foi encontrado.")
    elif lv == INDICIOS:
        parts.append(f"Há indícios não confirmados ({s['n_estruturas_hibridas_possiveis']} estruturas com fórmula "
                     f"compatível, {s['n_artigos_ambiguos']} artigos ambíguos): verificar manualmente antes de "
                     "considerar inédito.")
    else:
        if lv == LIGANTE:
            parts.append("Já cristalizado com haletos metálicos, mas apenas como ligante coordenado ou molécula "
                         "neutra — nenhum caso como cátion na rede (o papel usual em perovskitas híbridas).")
        else:
            parts.append(f"Como cátion em haleto metálico: {s['n_publicacoes_como_cation']} publicação(ões) de "
                         f"~{s['n_grupos_como_cation']} grupo(s) independente(s).")
        parts.append(f"Uso estrutural total: {s['n_publicacoes_estruturais']} publicação(ões) de "
                     f"~{s['n_grupos_independentes']} grupo(s), {s['n_estruturas_hibridas_fortes']} estrutura(s)"
                     + (f", {s['primeiro_ano']}–{s['ultimo_ano']}" if s["primeiro_ano"] else "") + ".")
        if s["papeis_do_organico"]:
            parts.append("Papel do orgânico nas entradas: " + ", ".join(f"{k} ({v})" for k, v in s["papeis_do_organico"].items()) + ".")
        if s["combinacoes_MX"]:
            parts.append("Combinações M–X (nº de artigos): " + ", ".join(f"{k} ({v})" for k, v in s["combinacoes_MX"].items()) + ".")
        if s["dimensionalidades"]:
            parts.append("Dimensionalidades (COD): " + ", ".join(f"{k} ({v})" for k, v in sorted(s["dimensionalidades"].items())) + ".")
        if lv == POUCO:
            parts.append("Pouco reproduzido: vale reproduzir e explorar variações de metal/haleto.")
        if s["n_estruturas_hibridas_possiveis"]:
            parts.append(f"Além disso, {s['n_estruturas_hibridas_possiveis']} indício(s) inconclusivo(s) na aba Estruturas.")
    parts.append("Viabilidade: " + s["viabilidade"] + ".")
    if s["n_sais_haleto_sem_metal"]:
        parts.append(f"{s['n_sais_haleto_sem_metal']} estrutura(s) do sal de haleto sem metal (precursor caracterizado).")
    if org.get("grupos_acidos"):
        parts.append("Tem grupo ácido (" + ", ".join(sorted(set(org["grupos_acidos"]))) + "): pode coordenar ao metal "
                     "ou formar zwitterion.")
    if s["n_candidatas_nao_verificadas"]:
        parts.append(f"Pendente: {s['n_candidatas_nao_verificadas']} entradas do COD com fórmula compatível não tiveram "
                     "o CIF verificado (limite); use --max-cifs maior para concluir.")
    if not s["busca_completa"]:
        parts.append("ATENÇÃO: busca incompleta (falha em: " + ", ".join(s["fontes_com_falha"]) + "); rode novamente.")
    return " ".join(parts)
