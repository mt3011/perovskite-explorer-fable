"""Nomes: limpeza de sinônimos, geração de variantes catiônicas e busca de termos em texto.

Em perovskitas híbridas o orgânico quase sempre aparece protonado, com outro nome
(piridina -> piridínio). Aqui geramos essas variantes por regras fixas de nomenclatura.
"""
from __future__ import annotations

import re
import unicodedata

_REGISTRY = re.compile(
    r"^(\d{2,7}-\d{2}-\d$|unii|chebi|chembl|dtxsid|dtxcid|nsc[\s\-_]?\d|einecs|ec \d|mfcd|akos|schembl|"
    r"zinc\d|brn |ccris|hsdb|ai3-|bdbm|mls\d|smr\d|stk\d|brd-|cas-|ncgc|tox21|bmse|cs-\d|hy-|db\d{4}|"
    r"inchi|component|cid |sid |wln|bidd|gtpl|nci|bp-|sb\d|s\d{4}|z\d{6}|en300|f\d{4}-|pubchem|usepa|"
    r"epitope|spectrum|kbio|prestwick|lopac|sr-\d|bpbio|opera_|ids-|ccg-|cmap_|st\d{5}|as-\d|ft-\d|ks-)",
    re.I,
)


def normalize_text(text: str) -> str:
    """Minúsculas, sem acentos, hifens/apóstrofos unificados, espaços colapsados."""
    if not text:
        return ""
    text = unicodedata.normalize("NFKC", str(text))
    text = re.sub(r"<[^>]+>", " ", text)  # marcação JATS/HTML dos resumos
    for ch in "‐‑‒–—−":
        text = text.replace(ch, "-")
    for ch in "′’‘`´ʹ":
        text = text.replace(ch, "'")
    text = "".join(c for c in unicodedata.normalize("NFD", text) if unicodedata.category(c) != "Mn")
    text = re.sub(r"\s+", " ", text).strip().lower()
    return text


def clean_synonyms(synonyms, limit: int = 8) -> list:
    """Mantém só sinônimos que são nomes químicos legíveis (descarta CAS, códigos de catálogo etc.)."""
    out, seen = [], set()
    for s in synonyms or []:
        s = str(s).strip()
        low = s.lower()
        if not (5 <= len(s) <= 60) or low in seen:
            continue
        if _REGISTRY.match(low):
            continue
        letters = sum(c.isalpha() for c in s)
        digits = sum(c.isdigit() for c in s)
        if letters < 5 or digits > letters * 0.5:
            continue
        if s.isupper() and len(s) <= 6:  # siglas (DDS, NAM...) são ambíguas demais
            continue
        if not re.search(r"[a-z]{4}", low):
            continue
        if any(x in low for x in ("(van)", "[usan", "[inn", "[jan", "[mi]", "[who", "(tn)", "usp", "reference standard", "%", "analytical", "grade")):
            continue
        seen.add(low)
        out.append(s)
        if len(out) >= limit:
            break
    return out


_SUFFIX_RULES = [
    (re.compile(r"amine\b"), ["ammonium", "aminium"]),
    (re.compile(r"amide\b"), ["amidium"]),
    (re.compile(r"(?<!am)ine\b"), ["inium"]),
    (re.compile(r"ole\b"), ["olium"]),
    (re.compile(r"\bamino(?=[a-z])|(?<=[\-\d,'\s\(])amino(?=[a-z])|(?<=di)amino(?=[a-z])"), ["ammonio"]),
]
_SPECIAL_RULES = [
    (re.compile(r"^(\d)-aminobenzoic acid$"), [r"\1-carboxyanilinium", r"\1-carboxyphenylammonium"]),
    (re.compile(r"^pyridine-(\d)-carboxylic acid$"), [r"\1-carboxypyridinium", r"pyridin-1-ium-\1-carboxylic acid"]),
    (re.compile(r"^pyridine-(\d)-carboxamide$"), [r"\1-carbamoylpyridinium", r"\1-carbamoylpyridin-1-ium"]),
    (re.compile(r"^nicotinic acid$"), ["3-carboxypyridinium", "nicotinium"]),
    (re.compile(r"^isonicotinic acid$"), ["4-carboxypyridinium", "isonicotinium"]),
    (re.compile(r"^(nicotinamide|niacinamide)$"), ["3-carbamoylpyridinium", "nicotinamidium"]),
    (re.compile(r"^(.*)aniline$"), [r"\1anilinium"]),
]


def cation_variants(name: str) -> list:
    """Gera nomes prováveis da forma protonada (catiônica) a partir de um nome neutro."""
    n = normalize_text(name)
    if not n or n.endswith("ium"):
        return []
    out = []
    for rx, reps in _SPECIAL_RULES:
        if rx.search(n):
            out += [rx.sub(r, n) for r in reps]
    for rx, reps in _SUFFIX_RULES:
        if rx.search(n):
            out += [rx.sub(r, n) for r in reps]
    seen, res = {n}, []
    for v in out:
        if v not in seen and len(v) >= 6:
            seen.add(v)
            res.append(v)
    return res


def build_terms(input_name: str, title: str = "", iupac: str = "", synonyms=None,
                extra_cations=None, max_neutral: int = 8, max_cation: int = 10) -> dict:
    """Monta os termos de busca: nomes neutros e nomes catiônicos (sem duplicatas)."""
    neutral, seen = [], set()

    def add(lst, term, limit):
        t = normalize_text(term)
        if t and t not in seen and len(t) >= 5 and len(lst) < limit:
            seen.add(t)
            lst.append(t)

    for cand in [input_name, title, iupac] + clean_synonyms(synonyms, limit=max_neutral):
        if cand and not _REGISTRY.match(str(cand).lower()) and not re.fullmatch(r"[\d\-]+", str(cand)):
            add(neutral, cand, max_neutral)
    cations: list = []
    for cand in list(extra_cations or []):
        add(cations, cand, max_cation)
    for base in list(neutral):
        for v in cation_variants(base):
            add(cations, v, max_cation)
    return {"neutros": neutral, "cations": cations}


def term_regex(term: str) -> re.Pattern:
    """Regex para achar o termo como 'palavra inteira' (não casa dentro de nomes maiores)."""
    t = re.escape(normalize_text(term))
    t = t.replace(r"\ ", r"[\s\-]*").replace(r"\-", r"[\s\-]*")
    return re.compile(r"(?<![a-z0-9])" + t + r"(?![a-z])")


def find_terms(text: str, terms) -> list:
    """Lista os termos presentes em `text` (ambos são normalizados)."""
    norm = normalize_text(text)
    return [t for t in terms if t and term_regex(t).search(norm)]


def split_sentences(text: str) -> list:
    text = re.sub(r"<[^>]+>", " ", text or "")
    text = re.sub(r"\s+", " ", text).strip()
    return [s.strip() for s in re.split(r"(?<=[\.\!\?])\s+(?=[A-Z0-9\(\[])", text) if s.strip()]


def slugify(text: str) -> str:
    t = normalize_text(text)
    t = re.sub(r"[^a-z0-9]+", "_", t).strip("_")
    return t[:60] or "organico"
