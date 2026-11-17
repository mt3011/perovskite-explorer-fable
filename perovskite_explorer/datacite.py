"""CSD (Cambridge Structural Database) sem licença, através dos metadados públicos do DataCite.

O CCDC registra no DataCite um DOI (10.5517/ccdc.csd.*) para cada depósito do CSD. Esses
metadados são abertos e trazem: número CCDC, NOME QUÍMICO SISTEMÁTICO do composto (atribuído
pelos editores do CSD), ano e o artigo associado (autores, revista, DOI). Isso permite:
  * pesquisar o CSD pelo nome do orgânico (e variantes catiônicas);
  * listar todos os depósitos CSD de um artigo.
O que NÃO vem por aqui: fórmula, cela, grupo espacial e coordenadas. Para isso use o link do
depósito (visualização/CIF gratuitos no site do CCDC) + `analisar-cif`, ou a API do CSD (csd.py).
"""
from __future__ import annotations

import re

from .names import find_terms_in_chemname, normalize_text

API = "https://api.datacite.org/dois"
_CCDC = re.compile(r"(?:CCDC|CSD)\s*(\d{5,8})")
_GENERIC = {"crystal structure", "experimental 3d coordinates", "crystal system", "space group",
            "cell parameters", "crystallography"}

# nomes de elementos metálicos (como aparecem nos nomes do CSD) e radicais "-ato"
METAL_NAMES = {
    "lead": "Pb", "plumb": "Pb", "tin": "Sn", "stann": "Sn", "germanium": "Ge", "german": "Ge",
    "bismuth": "Bi", "antimony": "Sb", "antimon": "Sb", "indium": "In", "thallium": "Tl",
    "copper": "Cu", "cupr": "Cu", "silver": "Ag", "argent": "Ag", "gold": "Au", "aur": "Au",
    "cadmium": "Cd", "cadm": "Cd", "mercury": "Hg", "mercur": "Hg", "zinc": "Zn", "manganese": "Mn",
    "mangan": "Mn", "tellurium": "Te", "tellur": "Te", "palladium": "Pd", "pallad": "Pd",
    "platinum": "Pt", "platin": "Pt", "zirconium": "Zr", "zircon": "Zr", "hafnium": "Hf",
    "iron": "Fe", "ferr": "Fe", "cobalt": "Co", "nickel": "Ni", "chromium": "Cr", "chrom": "Cr",
    "vanadium": "V", "titanium": "Ti", "gallium": "Ga", "selenium": "Se", "molybdenum": "Mo",
    "molybd": "Mo", "tungsten": "W", "tungst": "W", "rhenium": "Re", "ruthenium": "Ru", "rhodium": "Rh",
    "iridium": "Ir", "osmium": "Os", "europium": "Eu", "ytterbium": "Yb", "aluminium": "Al",
    "aluminum": "Al", "lithium": "Li", "sodium": "Na", "potassium": "K", "rubidium": "Rb",
    "caesium": "Cs", "cesium": "Cs", "magnesium": "Mg", "calcium": "Ca", "strontium": "Sr",
    "barium": "Ba", "scandium": "Sc", "yttrium": "Y", "niobium": "Nb", "tantalum": "Ta",
    "lanthanum": "La", "cerium": "Ce", "praseodymium": "Pr", "neodymium": "Nd", "samarium": "Sm",
    "gadolinium": "Gd", "terbium": "Tb", "dysprosium": "Dy", "holmium": "Ho", "erbium": "Er",
    "thulium": "Tm", "lutetium": "Lu", "uranium": "U", "thorium": "Th", "arsenic": "As",
    "technetium": "Tc", "beryllium": "Be",
}
_ATE_STEMS = ("plumb", "stann", "german", "bismuth", "antimon", "cupr", "argent", "aur", "cadm", "mercur",
              "zinc", "mangan", "tellur", "pallad", "platin", "zircon", "hafn", "ferr", "cobalt", "nickel",
              "ind", "thall", "chrom", "molybd", "tungst", "rhen", "ruthen", "rhod", "irid", "osm")
_METAL_RX = re.compile(r"(?<![a-z])(" + "|".join(sorted(METAL_NAMES, key=len, reverse=True)) + r")(?![a-z])")
_ATE_RX = re.compile(r"(" + "|".join(_ATE_STEMS) + r")ates?(?![a-z])")
_HAL = {"fluor": "F", "chlor": "Cl", "brom": "Br", "iod": "I"}
_MULT = r"(?:di|tri|tetra|penta|hexa|hepta|octa|nona|deca|undeca|dodeca)?"
# haleto como ligante ("tetrabromo-", "(mu-chloro)"), fundido ao metalato ("hexachloroantimonate")
# ou como ânion livre ("chloride", "triiodide"). Halogênio de substituinte orgânico
# ("2-chloro-...", "4-bromophenyl") NÃO casa: exige ausência de localizador e fim de palavra.
_HAL_LIGAND = re.compile(r"(?<![a-z\d'])(?<![\d']-)" + _MULT + r"-?(fluor|chlor|brom|iod)o(?=[\-\)\]\s,]|$)(?!-\d)")
_HAL_FUSED = re.compile(r"(fluor|chlor|brom|iod)o(?:" + "|".join(_ATE_STEMS) + r")ates?(?![a-z])")
_HAL_ION = re.compile(r"(?<![a-z])" + _MULT + r"(fluor|chlor|brom|iod)ides?(?![a-z])")


def name_features(name: str) -> dict:
    """Extrai do nome sistemático do CSD: metais, haletos e se a rede é polimérica (catena)."""
    n = normalize_text(name).replace("μ", "mu")
    n = re.sub(r"mu\d*-", "mu-", n)
    metals = {METAL_NAMES[m] for m in _METAL_RX.findall(n)}
    metals |= {METAL_NAMES[m] for m in _ATE_RX.findall(n) if m in METAL_NAMES}
    n_lig = n.replace("mu-", " ")
    halides = {_HAL[h] for rx in (_HAL_LIGAND, _HAL_FUSED, _HAL_ION) for h in rx.findall(n_lig)}
    return {"metais": sorted(metals), "haletos": sorted(halides), "catena": "catena" in n}


def _record(item: dict) -> dict:
    att = item.get("attributes") or {}
    titles = " ".join(str(t.get("title", "")) for t in att.get("titles") or [])
    m = _CCDC.search(titles)
    names = [str(s.get("subject", "")).strip() for s in att.get("subjects") or []]
    names = [s for s in names if s and s.lower() not in _GENERIC]
    related = [str(r.get("relatedIdentifier", "")).lower() for r in att.get("relatedIdentifiers") or []
               if str(r.get("relatedIdentifierType", "")).upper() == "DOI"]
    authors, journal, year_art = "", "", None
    for d in att.get("descriptions") or []:
        text = str(d.get("description", ""))
        if text.startswith("Related Article:"):
            parts = text[len("Related Article:"):].split("|")
            authors = parts[0].strip()
            if len(parts) > 2:
                journal = parts[2].strip()
                try:
                    year_art = int(parts[1])
                except ValueError:
                    pass
            if not related:
                dm = re.search(r"doi:(\S+)", text)
                if dm:
                    related = [dm.group(1).lower()]
    ccdc = m.group(1) if m else None
    return {
        "ccdc": ccdc,
        "doi_deposito": att.get("doi"),
        "url": f"https://doi.org/{att.get('doi')}" if att.get("doi") else "",
        "nome": "; ".join(names),
        "ano": year_art or att.get("publicationYear"),
        "doi_artigo": related[0] if related else None,
        "autores": authors,
        "revista": journal,
        "titulo": titles,
    }


def _paged(http, query: str, max_records: int):
    out, total, page = [], 0, 1
    while len(out) < max_records:
        obj = http.get_json(API, params={"client-id": "ccdc.csd", "query": query, "page[size]": 500,
                                         "page[number]": page})
        data = (obj or {}).get("data") or []
        total = ((obj or {}).get("meta") or {}).get("total") or len(data)
        out += [_record(d) for d in data]
        if len(data) < 500 or len(out) >= total:
            break
        page += 1
    return out, total


def _phrase(term: str) -> str:
    return '"' + term.replace('"', " ").replace("\\", " ") + '"'


def search_names(http, terms: dict, halogens: list, max_records: int = 3000) -> dict:
    """Depósitos do CSD cujo nome químico cita o orgânico E contém haleto.

    Retorna {"registros": [...], "total_com_haleto": n, "total_qualquer": n, "truncado": bool}.
    """
    all_terms = terms.get("cations", []) + terms.get("neutros", [])
    if not all_terms:
        return {"registros": [], "total_com_haleto": 0, "total_qualquer": 0, "truncado": False}
    org_q = "subjects.subject:(" + " OR ".join(_phrase(t) for t in all_terms) + ")"
    stems = {"F": "fluor", "Cl": "chlor", "Br": "brom", "I": "iod"}
    hal = " OR ".join(f"*{stems[x]}o* OR *{stems[x]}ide*" for x in halogens if x in stems)
    records, total = _paged(http, f"{org_q} AND subjects.subject:({hal})", max_records)
    # contagem de todas as entradas do CSD que citam o orgânico (qualquer tipo de composto)
    obj = http.get_json(API, params={"client-id": "ccdc.csd", "query": org_q, "page[size]": 1})
    total_any = ((obj or {}).get("meta") or {}).get("total") or 0
    return {"registros": records, "total_com_haleto": total, "total_qualquer": total_any,
            "truncado": total > len(records)}


def depositions_for_doi(http, article_doi: str) -> list:
    """Todos os depósitos CSD associados a um artigo (com nome químico de cada um)."""
    article_doi = article_doi.lower().strip()
    recs, _ = _paged(http, f'relatedIdentifiers.relatedIdentifier:"{article_doi}"', 500)
    recs = [r for r in recs if r["doi_artigo"] == article_doi]
    recs.sort(key=lambda d: d["ccdc"] or "")
    return recs


def to_structure(rec: dict, terms: dict, halogens: list, metals: list):
    """Converte um depósito em 'estrutura' no formato comum, se o nome citar o orgânico de fato.

    Retorna None quando o nome não cita o orgânico como composto próprio, não tem haleto
    iônico/ligante, ou só contém metais fora da lista-alvo.
    """
    from . import structures as st

    cation_hits = find_terms_in_chemname(rec["nome"], terms.get("cations", []))
    neutral_hits = find_terms_in_chemname(rec["nome"], terms.get("neutros", []))
    if not cation_hits and not neutral_hits:
        return None
    feat = name_features(rec["nome"])
    hal = [x for x in feat["haletos"] if x in halogens]
    if not hal:
        return None
    target = [m for m in feat["metais"] if m in metals]
    if feat["metais"] and not target:
        return None  # ex.: complexos de Li/Al/U — fora do escopo de haletos semicondutores
    cat = st.HYBRID if target else st.SALT
    if cation_hits:
        role = st.CATION
    elif target:
        role = "forma neutra no nome (ligante coordenado ou co-cristal)"
    else:
        role = st.NEUTRAL
    if target:
        arch = ("rede polimérica (catena-): cadeias, camadas ou rede 3D" if feat["catena"]
                else "unidades discretas (sem catena-): provável 0D ou complexo molecular")
        arch = f"Pelo nome do CSD: {arch} [M = {', '.join(feat['metais'])}; X = {', '.join(hal)}]"
    else:
        arch = "Pelo nome do CSD: sal de haleto sem metal"
    return {
        "fonte": "CSD",
        "id": f"CCDC {rec['ccdc']}" if rec.get("ccdc") else (rec.get("doi_deposito") or "?"),
        "ccdc": rec.get("ccdc"),
        "url": rec.get("url", ""),
        "formula_txt": "",
        "formula": {},
        "nome": rec["nome"],
        "titulo": "",
        "autores": rec.get("autores", ""),
        "revista": rec.get("revista", ""),
        "ano": rec.get("ano"),
        "doi": rec.get("doi_artigo"),
        "grupo_espacial": "",
        "cela": [None] * 6,
        "origem": ["nome químico no CSD (DataCite)"],
        "casamento_formula": None,
        "termos_no_nome": cation_hits + neutral_hits,
        "termos_no_titulo": [],
        "categoria": cat,
        "composto": normalize_text(rec["nome"]),
        "metais": feat["metais"],
        "haletos": hal,
        "prioridade": "A",
        "verificada": False,
        "papel_do_organico": role,
        "analise": {"ok": False, "erro": "CIF não disponível por esta via (abra o link do depósito)"},
        "classificacao": arch,
        "confianca": st.NOME_CSD,
        "motivo_confianca": "nome sistemático atribuído pelo CSD cita o orgânico ("
                            + ", ".join(cation_hits + neutral_hits) + ")",
    }
