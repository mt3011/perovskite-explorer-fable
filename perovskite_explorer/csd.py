"""Conector para o Cambridge Structural Database (CSD) via CSD Python API.

ESTADO: pronto para uso, porém NÃO TESTADO — foi escrito seguindo a documentação pública da
CSD Python API (https://downloads.ccdc.cam.ac.uk/documentation/API/), sem acesso a uma
licença. Ele só é ativado quando o pacote `ccdc` está instalado (isto é, quando você tiver
o CSD/licença). Sem isso, o programa segue normalmente com as demais fontes.

Como funciona: busca por subestrutura usando o esqueleto de átomos pesados do orgânico
(SMARTS com ligações genéricas, para aceitar a forma protonada/aromática) + exigência de
um halogênio na estrutura. Cada acerto é devolvido no mesmo formato das entradas do COD,
junto com o CIF gerado pelo CSD, e passa pela MESMA verificação de grafo e análise
estrutural (cif_analysis). Assim os critérios de confiança são idênticos entre os bancos.

Se o seu acesso for por outra interface (ex.: uma API web), basta reescrever `search`
mantendo o formato de retorno: (lista_de_entradas, {id: texto_do_cif}).
"""
from __future__ import annotations

import logging

from .chem import ATOMIC_NUMBER, parse_formula

log = logging.getLogger("perovskite_explorer")


def available() -> bool:
    try:
        import ccdc  # noqa: F401
        return True
    except Exception:  # noqa: BLE001 - ausência de pacote ou de licença
        return False


def graph_to_smarts(heavy_elements: list, heavy_edges: list) -> str:
    """SMARTS do esqueleto: átomos por número atômico, ligações de qualquer ordem (~)."""
    n = len(heavy_elements)
    adj: dict = {i: [] for i in range(n)}
    for a, b in heavy_edges:
        adj[a].append(b)
        adj[b].append(a)
    visited, closures, counter = set(), {}, [0]
    ring_marks: dict = {i: [] for i in range(n)}

    # 1ª passada: descobre arestas de fechamento de anel
    def dfs(u, parent):
        visited.add(u)
        for v in sorted(adj[u]):
            if v == parent:
                continue
            if v in visited:
                key = (min(u, v), max(u, v))
                if key not in closures:
                    counter[0] += 1
                    closures[key] = counter[0]
                    ring_marks[u].append(counter[0])
                    ring_marks[v].append(counter[0])
            else:
                tree[u].append(v)
                dfs(v, u)

    tree: dict = {i: [] for i in range(n)}
    dfs(0, -1)

    def mark(k):
        return f"~{k}" if k < 10 else f"~%{k}"

    def write(u):
        s = f"[#{ATOMIC_NUMBER[heavy_elements[u]]}]" + "".join(mark(k) for k in ring_marks[u])
        kids = tree[u]
        for i, v in enumerate(kids):
            sub = "~" + write(v)
            s += sub if i == len(kids) - 1 else f"({sub})"
        return s

    return write(0)


def search(organic: dict, max_hits: int = 2000):
    """Retorna (entradas, cifs). Lança ImportError se a CSD Python API não estiver disponível."""
    from ccdc import search as csd_search

    smarts = graph_to_smarts(organic["heavy_elements"], organic["heavy_edges"])
    query = csd_search.SubstructureSearch()
    query.add_substructure(csd_search.SMARTSSubstructure(smarts))
    query.add_substructure(csd_search.SMARTSSubstructure("[F,Cl,Br,I]"))
    query.settings.max_hit_structures = max_hits
    entries, cifs, seen = [], {}, set()
    for hit in query.search():
        refcode = hit.identifier
        if refcode in seen:
            continue
        seen.add(refcode)
        try:
            entry = hit.entry
            crystal = entry.crystal
            pub = entry.publication
            lengths = list(crystal.cell_lengths) if crystal.cell_lengths else [None] * 3
            angles = list(crystal.cell_angles) if crystal.cell_angles else [None] * 3
            formula_txt = entry.formula or ""
            rec = {
                "fonte": "CSD",
                "id": refcode,
                "url": f"https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid={refcode}",
                "formula_txt": formula_txt,
                "formula": parse_formula(formula_txt),
                "nome": entry.chemical_name or "",
                "titulo": "",
                "autores": getattr(pub, "authors", "") or "",
                "revista": str(getattr(pub, "journal", "") or ""),
                "ano": getattr(pub, "year", None),
                "doi": (getattr(pub, "doi", None) or "").lower() or None,
                "grupo_espacial": crystal.spacegroup_symbol or "",
                "cela": [round(float(x), 4) if x is not None else None for x in lengths + angles],
                "volume": round(float(crystal.cell_volume), 2) if crystal.cell_volume else None,
                "Z": crystal.z_value,
                "temperatura_K": entry.temperature,
                "R": entry.r_factor,
                "metodo": "",
                "status": "",
                "duplicata_de": None,
                "ccdc": getattr(entry, "ccdc_number", None),
                "origem": ["subestrutura CSD"],
            }
            cifs[refcode] = crystal.to_string("cif")
            entries.append(rec)
        except Exception as exc:  # noqa: BLE001 - um registro problemático não interrompe a busca
            log.warning("CSD: entrada %s ignorada (%s)", refcode, exc)
    return entries, cifs
