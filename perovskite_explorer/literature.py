"""Busca na literatura científica e classificação (por regras) do tipo de uso do orgânico.

Fontes: OpenAlex, Europe PMC e Crossref (abertas, sem chave obrigatória); Semantic Scholar
é opcional. Cada artigo recuperado é reavaliado localmente: só conta como evidência se o
nome do orgânico aparecer no título/resumo junto com contexto de haleto metálico/perovskita.
"""
from __future__ import annotations

import logging
import re

from .http import SourceError
from .names import find_terms, normalize_text, split_sentences, term_regex

log = logging.getLogger("perovskite_explorer")

# termos de contexto enviados às APIs (parte "E" da consulta booleana)
CONTEXT_TERMS = [
    "perovskite", "perovskites", "perovskitoid", "lead iodide", "lead bromide", "lead chloride",
    "lead halide", "tin iodide", "tin halide", "bismuth iodide", "bismuth halide", "antimony halide",
    "metal halide", "hybrid halide", "iodoplumbate", "bromoplumbate", "chloroplumbate",
    "iodobismuthate", "bromobismuthate", "iodoantimonate", "halometallate", "halometalate",
    "iodostannate", "iodocuprate", "iodoargentate", "organic-inorganic hybrid",
]
# contexto exigido localmente no título/resumo (haleto metálico / perovskita), deliberadamente
# estrito: palavras soltas como "hybrid" ou "iodide" não bastam.
_METAL = (r"(lead|tin|bismuth|antimony|copper|silver|gold|germanium|indium|thallium|cadmium|zinc|"
          r"manganese|mercury|tellurium|palladium|platinum|zirconium|hafnium)")
_HAL = r"(iodide|bromide|chloride|fluoride|halide)s?"
CONTEXT_RX = re.compile(
    r"perovskit|halometal|metal[- ]halide|"
    r"(plumb|bismuth|antimon|stann|cupr|argent|zinc|cadm|mercur|mangan|indat|tellur|german|aur|pallad|"
    r"platin|zircon|hafn)ates?\b|"
    + _METAL + r"(\s?\((ii|iii|iv|i)\))?[- ](based[- ])?" + _HAL + r"|"
    + _HAL + r"[- ](hybrid|perovskit|semiconductor)|"
    r"hybrid[- ](metal[- ])?" + _HAL + r"|"
    r"organic[- ]inorganic[- ](hybrid[- ])?(metal[- ])?(halide|iodide|bromide|chloride|perovskit)|"
    r"(pb|sn|ge|bi|sb|cu|ag|mn|zn|cd|hg|in|te)\s?\d?\s?(i|br|cl)\s?\d"
)

# pistas de uso ESTRUTURAL (o orgânico faz parte da rede cristalina): (rótulo, regex, peso)
STRUCT_CUES = [
    ("cristalografia", r"crystal structure|single[- ]crystal|x-ray diffraction|crystallograph|space group|"
                       r"crystalli[sz]es? in|isostructural|structurally characteri[sz]ed|structure determination", 3),
    ("dimensionalidade", r"\b[0-3]d\b|zero-dimensional|one-dimensional|two-dimensional|three-dimensional|"
                         r"low-dimensional|layered|ruddlesden|dion-jacobson|quantum[- ]well", 2),
    ("arquitetura", r"octahedra|\bchains?\b|inorganic layers?|corner-sharing|edge-sharing|face-sharing|"
                    r"templat|spacer|organic cations?|a-site|counter-?cation", 2),
    ("halometalato", r"plumbate|bismuthate|antimonate|stannate|cuprate|argentate|halometal|iodometal|"
                     r"bromometal|chlorometal|organic[- ]inorganic|"
                     r"hybrid (compound|material|halide|perovskite|iodide|bromide|chloride|salt|crystal)", 2),
    ("propriedades do material", r"band ?gap|photoluminescen|semiconduct|self-trapped exciton|"
                                 r"broadband emission|white[- ]light|ferroelectric|phase transition|"
                                 r"thermochrom|dielectric|scintillat", 1),
    ("síntese de novo composto", r"\bnovel\b|\bnew\b|synthesi[sz]ed|synthesis of|we report|"
                                 r"(was|were) (prepared|obtained|grown)", 1),
]
# pistas de uso NÃO estrutural (aditivo, passivação, camada de transporte...)
ADD_CUES = [
    ("aditivo/passivação", r"passivat|additive|defects?\b|trap[- ]stat|grain boundar|interfac|"
                           r"surface (treatment|modif)|post-treatment|modifier|dopant|doping|doped", 3),
    ("dispositivo fotovoltaico", r"power conversion efficienc|\bpce\b|solar cells?|photovoltaic performance|"
                                 r"open-circuit|fill factor|device (performance|stability)", 2),
    ("processamento/camadas", r"hole[- ]transport|electron[- ]transport|\bhtl\b|\betl\b|"
                              r"self-assembled monolayer|buffer layer|interlayer|antisolvent|anti-solvent|"
                              r"precursor solution|crystallization (kinetics|control)|lewis base|"
                              r"non-?radiative recombination", 2),
]
SKIP_TYPES = {"peer-review", "dataset", "component", "paratext", "erratum", "retraction", "grant",
              "reference-entry", "editorial", "letter"}

ESTRUTURAL = "estrutural"
ADITIVO = "aditivo/passivação/dispositivo"
AMBIGUO = "ambíguo (pistas estruturais e de aditivo)"
MENCAO = "menção sem contexto claro"
NAO_VERIFICAVEL = "não verificável (orgânico fora do título/resumo)"


def _or(terms) -> str:
    return " OR ".join(f'"{t}"' for t in terms)


def _query_terms(terms: dict, n_cation: int = 6, n_neutral: int = 6) -> list:
    return terms.get("cations", [])[:n_cation] + terms.get("neutros", [])[:n_neutral]


# --------------------------------------------------------------------------- #
# Fontes
# --------------------------------------------------------------------------- #

def search_openalex(http, terms: dict, cfg: dict) -> list:
    q = f"({_or(_query_terms(terms))}) AND ({_or(CONTEXT_TERMS)})"
    params = {"search": q, "per-page": 100}
    if cfg.get("email"):
        params["mailto"] = cfg["email"]
    if cfg.get("openalex_api_key"):
        params["api_key"] = cfg["openalex_api_key"]
    out = []
    for page in range(1, cfg.get("max_paginas", 2) + 1):
        obj = http.get_json("https://api.openalex.org/works", params={**params, "page": page})
        results = (obj or {}).get("results") or []
        for w in results:
            out.append(_from_openalex(w))
        if len(results) < 100:
            break
    return out


def _from_openalex(w: dict) -> dict:
    inv = w.get("abstract_inverted_index") or {}
    abstract = ""
    if inv:
        pos = sorted((p, word) for word, ps in inv.items() for p in ps)
        abstract = " ".join(word for _, word in pos)
    src = ((w.get("primary_location") or {}).get("source") or {})
    doi = (w.get("doi") or "").lower().replace("https://doi.org/", "") or None
    return {
        "doi": doi,
        "titulo": w.get("display_name") or w.get("title") or "",
        "resumo": abstract,
        "ano": w.get("publication_year"),
        "revista": src.get("display_name") or "",
        "autores": [((a.get("author") or {}).get("display_name") or "") for a in w.get("authorships") or []],
        "citacoes": w.get("cited_by_count"),
        "tipo": w.get("type") or "",
        "retratado": bool(w.get("is_retracted")),
        "fontes": ["OpenAlex"],
        "url": w.get("doi") or w.get("id") or "",
    }


def search_europepmc(http, terms: dict, cfg: dict) -> list:
    q = f"({_or(_query_terms(terms))}) AND ({_or(CONTEXT_TERMS)}) NOT SRC:PAT"
    out, cursor = [], "*"
    for _ in range(cfg.get("max_paginas", 2)):
        obj = http.get_json("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                            params={"query": q, "format": "json", "resultType": "core",
                                    "pageSize": 100, "cursorMark": cursor})
        results = ((obj or {}).get("resultList") or {}).get("result") or []
        for r in results:
            out.append(_from_epmc(r))
        nxt = (obj or {}).get("nextCursorMark")
        if len(results) < 100 or not nxt or nxt == cursor:
            break
        cursor = nxt
    return out


def _from_epmc(r: dict) -> dict:
    authors = [a.get("fullName") or "" for a in ((r.get("authorList") or {}).get("author") or [])]
    if not authors and r.get("authorString"):
        authors = [a.strip(" .") for a in r["authorString"].split(",")]
    journal = (((r.get("journalInfo") or {}).get("journal") or {}).get("title")) or r.get("journalTitle") or ""
    types = [str(t).lower() for t in ((r.get("pubTypeList") or {}).get("pubType") or [])]
    doi = (r.get("doi") or "").lower() or None
    return {
        "doi": doi,
        "titulo": r.get("title") or "",
        "resumo": r.get("abstractText") or "",
        "ano": _int(r.get("pubYear")),
        "revista": journal,
        "autores": authors,
        "citacoes": _int(r.get("citedByCount")),
        "tipo": "review" if any("review" in t for t in types) else "preprint" if r.get("source") == "PPR" else "article",
        "retratado": any("retract" in t for t in types),
        "fontes": ["Europe PMC"],
        "url": f"https://doi.org/{doi}" if doi else f"https://europepmc.org/article/{r.get('source')}/{r.get('id')}",
    }


def search_crossref(http, terms: dict, cfg: dict) -> list:
    out = []
    picks = terms.get("cations", [])[:2] + terms.get("neutros", [])[:2]
    for t in picks:
        params = {"query.bibliographic": f"{t} perovskite halide", "rows": 40,
                  "select": "DOI,title,abstract,author,issued,container-title,is-referenced-by-count,type"}
        if cfg.get("email"):
            params["mailto"] = cfg["email"]
        obj = http.get_json("https://api.crossref.org/works", params=params)
        for it in ((obj or {}).get("message") or {}).get("items") or []:
            out.append(_from_crossref(it))
    return out


def _from_crossref(it: dict) -> dict:
    doi = (it.get("DOI") or "").lower() or None
    year = None
    try:
        year = int(it["issued"]["date-parts"][0][0])
    except (KeyError, IndexError, TypeError, ValueError):
        pass
    authors = [" ".join(x for x in (a.get("given"), a.get("family")) if x) or a.get("name", "")
               for a in it.get("author") or []]
    return {
        "doi": doi,
        "titulo": " ".join(it.get("title") or []),
        "resumo": re.sub(r"<[^>]+>", " ", it.get("abstract") or ""),
        "ano": year,
        "revista": " ".join((it.get("container-title") or [])[:1]),
        "autores": authors,
        "citacoes": it.get("is-referenced-by-count"),
        "tipo": (it.get("type") or "").replace("journal-article", "article").replace("posted-content", "preprint"),
        "retratado": False,
        "fontes": ["Crossref"],
        "url": f"https://doi.org/{doi}" if doi else "",
    }


def crossref_lookup(http, doi: str, cfg: dict):
    """Metadados de um DOI específico (usado para os artigos ligados a estruturas)."""
    params = {"mailto": cfg["email"]} if cfg.get("email") else None
    from urllib.parse import quote

    obj = http.get_json(f"https://api.crossref.org/works/{quote(doi, safe='/')}", params=params)
    msg = (obj or {}).get("message")
    return _from_crossref(msg) if msg else None


def search_semanticscholar(http, terms: dict, cfg: dict) -> list:
    org = " | ".join(f'"{t}"' for t in _query_terms(terms, 4, 4))
    ctx = " | ".join(f'"{t}"' for t in CONTEXT_TERMS[:14])
    headers = {"x-api-key": cfg["semanticscholar_api_key"]} if cfg.get("semanticscholar_api_key") else None
    obj = http.get_json("https://api.semanticscholar.org/graph/v1/paper/search/bulk",
                        params={"query": f"({org}) + ({ctx})",
                                "fields": "title,abstract,year,venue,externalIds,citationCount,authors"},
                        headers=headers)
    out = []
    for p in (obj or {}).get("data") or []:
        doi = ((p.get("externalIds") or {}).get("DOI") or "").lower() or None
        out.append({
            "doi": doi, "titulo": p.get("title") or "", "resumo": p.get("abstract") or "",
            "ano": p.get("year"), "revista": p.get("venue") or "",
            "autores": [a.get("name") or "" for a in p.get("authors") or []],
            "citacoes": p.get("citationCount"), "tipo": "article", "retratado": False,
            "fontes": ["Semantic Scholar"],
            "url": f"https://doi.org/{doi}" if doi else f"https://www.semanticscholar.org/paper/{p.get('paperId')}",
        })
    return out


SOURCES = {
    "openalex": ("OpenAlex", search_openalex),
    "europepmc": ("Europe PMC", search_europepmc),
    "crossref": ("Crossref", search_crossref),
    "semanticscholar": ("Semantic Scholar", search_semanticscholar),
}


# --------------------------------------------------------------------------- #
# Fusão, relevância e classificação
# --------------------------------------------------------------------------- #

def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _paper_key(p: dict) -> str:
    if p.get("doi"):
        return "doi:" + p["doi"]
    return "t:" + re.sub(r"[^a-z0-9]", "", normalize_text(p.get("titulo", "")))[:80] + str(p.get("ano") or "")


def merge(papers: list) -> list:
    """Une registros do mesmo artigo vindos de fontes diferentes (por DOI ou título+ano)."""
    merged: dict = {}
    for p in papers:
        if not p.get("titulo") and not p.get("doi"):
            continue
        k = _paper_key(p)
        if k not in merged:
            merged[k] = dict(p)
            continue
        m = merged[k]
        m["fontes"] = sorted(set(m["fontes"]) | set(p["fontes"]))
        if len(p.get("resumo") or "") > len(m.get("resumo") or ""):
            m["resumo"] = p["resumo"]
        if len(p.get("autores") or []) > len(m.get("autores") or []):
            m["autores"] = p["autores"]
        cites = [c for c in (m.get("citacoes"), p.get("citacoes")) if c is not None]
        m["citacoes"] = max(cites) if cites else None
        for f in ("ano", "revista", "titulo", "url", "tipo"):
            m[f] = m.get(f) or p.get(f)
        m["retratado"] = m.get("retratado") or p.get("retratado")
    return list(merged.values())


def classify(p: dict, terms: dict) -> dict:
    """Decide relevância e tipo de uso. Guarda as pistas encontradas para auditoria."""
    all_terms = terms.get("cations", []) + terms.get("neutros", [])
    title, abstract = p.get("titulo") or "", p.get("resumo") or ""
    in_title = find_terms(title, all_terms)
    in_abs = find_terms(abstract, all_terms)
    text = normalize_text(title + " . " + abstract)
    p["termos_no_titulo"] = in_title
    p["termos_no_resumo"] = in_abs
    p["tem_resumo"] = bool(abstract.strip())
    p["trecho"] = ""
    p["pistas_estruturais"], p["pistas_aditivo"] = [], []
    if not in_title and not in_abs:
        p["uso"] = NAO_VERIFICAVEL
        p["relevante"] = False
        return p
    if not CONTEXT_RX.search(text):
        p["uso"] = MENCAO
        p["relevante"] = False
        p["trecho"] = "(orgânico citado, mas sem contexto de haleto metálico/perovskita no título/resumo)"
        return p
    p["relevante"] = True
    hits = in_title + in_abs
    # frases que citam o orgânico pesam o dobro do restante do resumo
    rxs = [term_regex(t) for t in set(hits)]
    focus = [s for s in split_sentences(abstract) if any(r.search(normalize_text(s)) for r in rxs)]
    focus_text = normalize_text(title + " . " + " ".join(focus))
    if focus:
        p["trecho"] = focus[0][:400]
    s_score = a_score = 0.0
    for label, rx, w in STRUCT_CUES:
        if re.search(rx, focus_text):
            s_score += w
            p["pistas_estruturais"].append(label)
        elif re.search(rx, text):
            s_score += w / 2
            p["pistas_estruturais"].append(label + " (fora da frase do orgânico)")
    for label, rx, w in ADD_CUES:
        if re.search(rx, focus_text):
            a_score += w
            p["pistas_aditivo"].append(label)
        elif re.search(rx, text):
            a_score += w / 2
            p["pistas_aditivo"].append(label + " (fora da frase do orgânico)")
    if any(t.endswith("ium") or "ammonio" in t for t in hits):
        s_score += 3
        p["pistas_estruturais"].append("orgânico citado na forma catiônica (-ium)")
    p["pontos_estrutural"], p["pontos_aditivo"] = s_score, a_score
    if s_score >= 4 and s_score >= a_score + 2:
        p["uso"] = ESTRUTURAL
    elif a_score >= 4 and a_score >= s_score + 2:
        p["uso"] = ADITIVO
    elif s_score >= 4 or a_score >= 4:
        p["uso"] = AMBIGUO
    else:
        p["uso"] = MENCAO
    return p


def author_key(name: str) -> str:
    """'John A. Smith', 'Smith JA' e 'Smith, John' -> 'smith j'."""
    n = normalize_text(name).replace(".", " ")
    if not n:
        return ""
    if "," in n:
        last, _, first = n.partition(",")
        parts_first = first.split()
        return f"{last.strip()} {parts_first[0][0] if parts_first else ''}".strip()
    parts = n.split()
    if len(parts) == 1:
        return parts[0]
    raw_last = str(name).replace(".", " ").split()[-1]
    if raw_last.isupper() and len(raw_last) <= 3:  # formato 'Smith JA'
        return f"{' '.join(parts[:-1])} {parts[-1][0]}"
    return f"{parts[-1]} {parts[0][0]}"


def independent_groups(papers: list) -> int:
    """Nº aproximado de grupos de pesquisa independentes (artigos ligados por coautoria).

    Dois artigos são do mesmo grupo se compartilham 2+ autores (ou 1, se algum tiver <= 2 autores).
    """
    sets = [set(filter(None, (author_key(a) for a in p.get("autores") or []))) for p in papers]
    sets = [s for s in sets if s]
    parent = list(range(len(sets)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            shared = len(sets[i] & sets[j])
            if shared >= 2 or (shared == 1 and min(len(sets[i]), len(sets[j])) <= 2):
                parent[find(i)] = find(j)
    return len({find(i) for i in range(len(sets))})


def run(http, terms: dict, cfg: dict, sources: list):
    """Executa as fontes pedidas. Retorna (artigos_classificados, relatório_por_fonte)."""
    report, raw = {}, []
    if not (terms.get("cations") or terms.get("neutros")):
        return [], {"aviso": "sem termos de busca"}
    for key in sources:
        if key not in SOURCES:
            continue
        label, fn = SOURCES[key]
        try:
            found = fn(http, terms, cfg)
            report[label] = {"ok": True, "registros": len(found)}
            raw += found
        except SourceError as exc:
            report[label] = {"ok": False, "erro": str(exc)}
            log.warning("Literatura: %s falhou: %s", label, exc)
        except Exception as exc:  # noqa: BLE001 - formato inesperado de uma fonte não derruba as outras
            report[label] = {"ok": False, "erro": f"resposta inesperada ({type(exc).__name__}: {exc})"}
            log.warning("Literatura: %s devolveu formato inesperado: %s", label, exc)
    papers = [classify(p, terms) for p in merge(raw)
              if (p.get("tipo") or "").lower() not in SKIP_TYPES]
    return papers, report
