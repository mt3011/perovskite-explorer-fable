"""Ponte para o CSD (Cambridge) sem licença: depósitos CCDC ligados a um artigo.

O CCDC registra um DOI (10.5517/ccdc.csd.*) no DataCite para cada depósito do CSD, com o DOI
do artigo associado. A API pública do DataCite permite, portanto, descobrir quantos
depósitos CSD um artigo tem e seus números CCDC. As coordenadas atômicas NÃO vêm por aqui:
para isso é preciso o site "Access Structures" do CCDC ou a API do CSD (ver csd.py).
"""
from __future__ import annotations

import re

API = "https://api.datacite.org/dois"
_CCDC = re.compile(r"(?:CCDC|CSD)\s*(\d{5,8})")


def depositions_for_doi(http, article_doi: str) -> list:
    article_doi = article_doi.lower().strip()
    obj = http.get_json(API, params={
        "query": f'relatedIdentifiers.relatedIdentifier:"{article_doi}"',
        "client-id": "ccdc.csd",
        "page[size]": 100,
    })
    out = []
    for item in (obj or {}).get("data") or []:
        att = item.get("attributes") or {}
        related = [str(r.get("relatedIdentifier", "")).lower() for r in att.get("relatedIdentifiers") or []]
        descs = " | ".join(str(d.get("description", "")) for d in att.get("descriptions") or [])
        if article_doi not in related and article_doi not in descs.lower():
            continue  # a busca do DataCite é textual; confirmamos o vínculo com o DOI
        titles = " ".join(str(t.get("title", "")) for t in att.get("titles") or [])
        m = _CCDC.search(titles)
        subjects = "; ".join(str(s.get("subject", "")) for s in att.get("subjects") or [])
        out.append({
            "ccdc": m.group(1) if m else None,
            "doi_deposito": att.get("doi"),
            "url": f"https://doi.org/{att.get('doi')}" if att.get("doi") else None,
            "titulo": titles,
            "ano": att.get("publicationYear"),
            "doi_artigo": article_doi,
            "detalhes": subjects,
        })
    out.sort(key=lambda d: d["ccdc"] or "")
    return out
