"""Crystallography Open Database (COD) — banco cristalográfico aberto, com CIFs livres.

API REST: https://wiki.crystallography.net/RESTful_API/
Estratégia: (1) consultas por elementos (orgânico + metal + haleto), que não dependem de
nomes; (2) consultas por texto com os nomes/variantes catiônicas do orgânico.
"""
from __future__ import annotations

import logging

from .chem import parse_formula
from .http import SourceError

log = logging.getLogger("perovskite_explorer")
ENDPOINTS = ["https://www.crystallography.net/cod/result", "https://www.crystallography.net/cod/result.php"]
CIF_URL = "https://www.crystallography.net/cod/{id}.cif"
PAGE_URL = "https://www.crystallography.net/cod/{id}.html"


def _query(http, params: dict) -> list:
    params = {**params, "format": "json"}
    obj = None
    for n, url in enumerate(ENDPOINTS):
        obj = http.get_json(url, params=params)
        if obj is not None or n == len(ENDPOINTS) - 1:
            break
    if obj is None:
        return []
    if isinstance(obj, dict):
        obj = list(obj.values())
    return [e for e in obj if isinstance(e, dict)]


def _f(v):
    try:
        return round(float(v), 4)
    except (TypeError, ValueError):
        return None


def normalize_entry(raw: dict) -> dict:
    cod_id = str(raw.get("file") or raw.get("id") or "").strip()
    formula_txt = str(raw.get("formula") or raw.get("calcformula") or "").strip(" -")
    names = [str(raw.get(k)).strip() for k in ("chemname", "commonname", "mineral")
             if raw.get(k) and str(raw.get(k)).strip() not in ("?", "None")]
    year = raw.get("year")
    try:
        year = int(str(year)[:4])
    except (TypeError, ValueError):
        year = None
    return {
        "fonte": "COD",
        "id": cod_id,
        "url": PAGE_URL.format(id=cod_id),
        "formula_txt": formula_txt,
        "formula": parse_formula(formula_txt),
        "nome": "; ".join(dict.fromkeys(names)),
        "titulo": str(raw.get("title") or "").strip(),
        "autores": str(raw.get("authors") or "").strip(),
        "revista": str(raw.get("journal") or "").strip(),
        "ano": year,
        "doi": str(raw.get("doi") or "").strip().lower() or None,
        "grupo_espacial": str(raw.get("sg") or "").strip(),
        "cela": [_f(raw.get(k)) for k in ("a", "b", "c", "alpha", "beta", "gamma")],
        "volume": _f(raw.get("vol")),
        "Z": raw.get("Z"),
        "temperatura_K": _f(raw.get("celltemp")) or _f(raw.get("diffrtemp")),
        "R": _f(raw.get("Robs")) or _f(raw.get("Rall")),
        "metodo": str(raw.get("method") or "").strip(),
        "status": str(raw.get("status") or "").strip(),
        "duplicata_de": raw.get("duplicateof") or None,
    }


def search(http, organic: dict, terms: dict, metals: list, halogens: list, max_text_queries: int = 8):
    """Coleta entradas candidatas do COD. Retorna (entradas, relatório_das_consultas)."""
    entries: dict = {}
    report = {"consultas": 0, "falhas": [], "brutos": 0}

    def run(params, origin):
        report["consultas"] += 1
        try:
            raws = _query(http, params)
        except SourceError as exc:
            report["falhas"].append(f"{origin}: {exc}")
            return
        report["brutos"] += len(raws)
        for raw in raws:
            e = normalize_entry(raw)
            if not e["id"]:
                continue
            entries.setdefault(e["id"], e).setdefault("origem", set()).add(origin.split(":")[0])

    formula = organic.get("formula") or {}
    org_elements = sorted(e for e in formula if e != "H")
    if org_elements:
        for x in halogens:
            for m in metals:
                els = list(dict.fromkeys(org_elements + [m, x]))
                if len(els) > 8:  # a API aceita até 8 elementos obrigatórios
                    keep = [e for e in ("C", "N", m, x) if e in els]
                    els = keep + [e for e in els if e not in keep][: 8 - len(keep)]
                params = {f"el{i + 1}": e for i, e in enumerate(els)}
                run(params, f"elementos:{m}-{x}")
                if len(report["falhas"]) >= 6 and report["falhas"] and not entries:
                    # servidor claramente fora do ar: não insiste nas dezenas de consultas restantes
                    raise SourceError("COD indisponível: " + report["falhas"][-1])
            # sais de haleto sem metal (precursor orgânico·HX ou compostos livres de metal)
            els = list(dict.fromkeys(org_elements + [x]))[:8]
            n_el = len(set(formula) | {x})
            params = {f"el{i + 1}": e for i, e in enumerate(els)}
            params.update({"strictmin": n_el, "strictmax": n_el + 1})
            run(params, f"sem_metal:{x}")
    text_terms = (terms.get("cations", []) + terms.get("neutros", []))[:max_text_queries]
    for t in text_terms:
        run({"text": t}, f"texto:{t}")
    for e in entries.values():
        e["origem"] = sorted(e["origem"])
    if report["consultas"] and len(report["falhas"]) == report["consultas"]:
        raise SourceError("COD indisponível: " + report["falhas"][-1])
    return list(entries.values()), report


def fetch_cif(http, cod_id: str):
    return http.get_text(CIF_URL.format(id=cod_id), ttl_days=3650)
