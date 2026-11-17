"""Química básica: elementos, fórmulas e classificação de sítios do orgânico.

Tudo aqui é determinístico e funciona offline.
"""
from __future__ import annotations

import re
from collections import Counter

SYMBOLS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
]
ATOMIC_NUMBER = {s: i for i, s in enumerate(SYMBOLS)}

# Raios covalentes (Cordero et al., Dalton Trans. 2008), em Å.
COVALENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66,
    "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Ar": 1.06, "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.50, "Fe": 1.42, "Co": 1.38, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.20,
    "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16, "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40, "Cs": 2.44, "Ba": 2.15,
    "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96,
    "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87, "Lu": 1.87, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
    "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Po": 1.40, "At": 1.50, "Rn": 1.50, "Fr": 2.60, "Ra": 2.21,
    "Ac": 2.15, "Th": 2.06, "Pa": 2.00, "U": 1.96, "Np": 1.90, "Pu": 1.87, "Am": 1.80, "Cm": 1.69,
}

HALOGENS = frozenset({"F", "Cl", "Br", "I"})
# Elementos tratados como "parte orgânica" na análise estrutural.
ORGANIC_ELEMENTS = frozenset({"H", "B", "C", "N", "O", "Si", "P", "S"})
NOBLE = frozenset({"He", "Ne", "Ar", "Kr", "Xe", "Rn"})
# Qualquer outro elemento é tratado como centro metálico/metaloide (Pb, Sn, Bi, Sb, Te, Ge, Cu...).
METALS = frozenset(s for s in SYMBOLS[1:] if s not in ORGANIC_ELEMENTS | HALOGENS | NOBLE | {"At"})
# Cátions grandes que ocupam o sítio A (não fazem parte do arcabouço M-X).
A_SITE_METALS = frozenset({"Cs", "Rb", "Fr", "Ba", "Sr", "Ra"})
ALKALI_LIGHT = frozenset({"Li", "Na", "K", "Ca", "Mg"})

# Metais/metaloides pesquisados por padrão (formadores de haletos híbridos semicondutores).
DEFAULT_METALS = [
    "Pb", "Sn", "Ge", "Bi", "Sb", "In", "Tl", "Cu", "Ag", "Au",
    "Cd", "Hg", "Zn", "Mn", "Te", "Pd", "Pt", "Zr", "Hf",
]
EXTENDED_METALS = ["Fe", "Co", "Ni", "Cr", "V", "Ti", "Ga", "Se", "Mo", "W", "Re", "Ru", "Rh", "Ir", "Eu", "Yb"]
DEFAULT_HALOGENS = ["Cl", "Br", "I"]

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)\s*(\d+(?:\.\d+)?|\.\d+)?")


def parse_formula(text) -> dict:
    """Converte 'C12 H9 N2 Pb I3', '- C6 H7 N2 O -' ou 'C6H6N2O+' em {elemento: quantidade}.

    Deutério (D) e trítio (T) contam como H. Símbolos desconhecidos são ignorados.
    """
    out: Counter = Counter()
    if not text:
        return {}
    text = str(text).replace("·", " ").replace(",", " ")
    # remove cargas ("+", "2-") e parênteses simples (fórmulas-soma não os usam)
    text = re.sub(r"[\(\)\[\]]", " ", text)
    for sym, num in _FORMULA_TOKEN.findall(text):
        if sym in ("D", "T"):
            sym = "H"
        if sym not in ATOMIC_NUMBER or sym == "X":
            continue
        out[sym] += float(num) if num else 1.0
    return {k: (int(v) if abs(v - round(v)) < 1e-6 else v) for k, v in out.items() if v > 0}


def hill_formula(formula: dict) -> str:
    """Fórmula em notação de Hill (C, H, depois ordem alfabética)."""
    if not formula:
        return ""
    keys = sorted(formula)
    if "C" in formula:
        keys = ["C"] + (["H"] if "H" in formula else []) + [k for k in keys if k not in ("C", "H")]

    def fmt(v):
        if abs(v - round(v)) < 1e-6:
            v = int(round(v))
            return "" if v == 1 else str(v)
        return f"{v:g}"

    return "".join(f"{k}{fmt(formula[k])}" for k in keys)


def heavy(formula: dict) -> dict:
    return {k: v for k, v in formula.items() if k != "H"}


def metals_in(formula: dict) -> list:
    return sorted(e for e in formula if e in METALS)


def halogens_in(formula: dict) -> list:
    return sorted(e for e in formula if e in HALOGENS)


# --------------------------------------------------------------------------- #
# Grafo molecular (vindo do registro do PubChem) e classificação de sítios
# --------------------------------------------------------------------------- #

def largest_organic_unit(atoms: list, bonds: list):
    """Dado o grafo completo (com H), devolve o maior fragmento covalente que contém carbono.

    atoms: [{"id": int, "el": str, "charge": int}], bonds: [(id1, id2, ordem)].
    Útil quando o usuário lista um sal (ex.: 'cloridrato de X'): ficamos só com o orgânico.
    """
    parent = {a["id"]: a["id"] for a in atoms}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b, _ in bonds:
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    groups: dict = {}
    for a in atoms:
        groups.setdefault(find(a["id"]), []).append(a)

    def key(g):
        has_c = any(a["el"] == "C" for a in g)
        return (has_c, sum(1 for a in g if a["el"] != "H"))

    best = max(groups.values(), key=key)
    ids = {a["id"] for a in best}
    return best, [(a, b, o) for a, b, o in bonds if a in ids and b in ids], len(groups)


def describe_organic_graph(atoms: list, bonds: list) -> dict:
    """Extrai do grafo: fórmula, grafo de átomos pesados e sítios protonáveis/ácidos.

    A classificação de basicidade é uma heurística por grupos funcionais (documentada no guia).
    """
    el = {a["id"]: a["el"] for a in atoms}
    charge = {a["id"]: a.get("charge", 0) for a in atoms}
    nbrs: dict = {a["id"]: [] for a in atoms}
    for a, b, o in bonds:
        nbrs[a].append((b, o))
        nbrs[b].append((a, o))

    formula = dict(Counter(el.values()))
    heavy_ids = [i for i in el if el[i] != "H"]
    index = {i: n for n, i in enumerate(heavy_ids)}
    heavy_edges = sorted(
        (min(index[a], index[b]), max(index[a], index[b]))
        for a, b, _ in bonds if a in index and b in index
    )

    # átomos em anel: aresta (a,b) está em anel se a e b continuam conectados sem ela
    ring_atoms = _ring_atoms(heavy_ids, [(a, b) for a, b, _ in bonds if a in index and b in index])

    def n_h(i):
        return sum(1 for j, _ in nbrs[i] if el[j] == "H")

    def has_multiple(i, to=None):
        return any(o >= 2 and (to is None or el[j] in to) for j, o in nbrs[i])

    sites = []
    for i in heavy_ids:
        if el[i] != "N":
            continue
        heavy_n = [(j, o) for j, o in nbrs[i] if el[j] != "H"]
        orders = [o for _, o in nbrs[i]]
        nh = n_h(i)
        kind, basic = None, None
        if charge[i] > 0 and len(nbrs[i]) == 4 and nh == 0:
            kind, basic = "amônio quaternário (cátion permanente)", "permanente"
        elif charge[i] > 0 and any(el[j] == "O" for j, _ in heavy_n):
            kind, basic = "nitro/N-óxido", "não"
        elif charge[i] > 0:
            kind, basic = "N já protonado/catiônico", "permanente"
        elif 3 in orders:
            kind, basic = "nitrila", "não"
        elif 2 in orders:
            partner = [el[j] for j, o in heavy_n if o == 2][0]
            if partner == "O":
                kind, basic = "nitroso/nitro", "não"
            elif partner == "N":
                kind, basic = "azo/diazo", "fraco"
            elif i in ring_atoms:
                kind, basic = "N sp2 em anel (tipo piridina/imidazol)", "médio"
            else:
                kind, basic = "imina/amidina (N sp2)", "médio"
        else:
            # só ligações simples
            c_nb = [j for j, _ in heavy_n if el[j] == "C"]
            acyl = any(has_multiple(c, {"O", "S"}) for c in c_nb)
            sulfonyl = any(el[j] in ("S", "P") and has_multiple(j, {"O"}) for j, _ in heavy_n)
            amidine_like = any(has_multiple(c, {"N"}) for c in c_nb)
            unsat_c = any(has_multiple(c, {"C"}) for c in c_nb)
            if acyl:
                kind, basic = "amida/imida/ureia", "não"
            elif sulfonyl:
                kind, basic = "sulfonamida/fosforamida", "não"
            elif amidine_like:
                kind, basic = "NH de amidina/guanidina/aminopiridina (protona no N sp2 vizinho)", "não"
            elif unsat_c and i in ring_atoms:
                kind, basic = "N tipo pirrol/indol (par isolado no sistema π)", "não"
            elif unsat_c:
                kind, basic = "amina aromática (tipo anilina)", "fraco"
            elif any(el[j] == "N" for j, _ in heavy_n):
                kind, basic = "hidrazina", "médio"
            else:
                degree = {0: "amônia", 1: "amina primária", 2: "amina secundária", 3: "amina terciária"}
                kind, basic = degree.get(len(heavy_n), "amina") + " alifática", "forte"
        sites.append({"tipo": kind, "basicidade": basic, "n_H": nh})

    acids = []
    for i in heavy_ids:
        if el[i] == "O" and n_h(i) >= 1:
            anchor = [j for j, _ in nbrs[i] if el[j] != "H"]
            if anchor and el[anchor[0]] == "C" and has_multiple(anchor[0], {"O"}):
                acids.append("ácido carboxílico")
            elif anchor and el[anchor[0]] in ("S", "P"):
                acids.append("ácido sulfônico/fosfônico")
            elif anchor and el[anchor[0]] == "C" and has_multiple(anchor[0], {"C"}):
                acids.append("fenol/enol")
            elif anchor:
                acids.append("álcool")
        if el[i] == "S" and n_h(i) >= 1:
            acids.append("tiol")

    h_hetero = sum(n_h(i) for i in heavy_ids if el[i] in ("N", "O", "S", "P"))
    n_strong = sum(1 for s in sites if s["basicidade"] == "forte")
    n_medium = sum(1 for s in sites if s["basicidade"] == "médio")
    n_weak = sum(1 for s in sites if s["basicidade"] == "fraco")
    n_perm = sum(1 for s in sites if s["basicidade"] == "permanente")
    return {
        "formula": formula,
        "carga": sum(charge.values()),
        "heavy_elements": [el[i] for i in heavy_ids],
        "heavy_edges": [list(e) for e in heavy_edges],
        "sitios_N": sites,
        "grupos_acidos": acids,
        "h_em_heteroatomos": h_hetero,
        "n_basicos_fortes": n_strong,
        "n_basicos_medios": n_medium,
        "n_basicos_fracos": n_weak,
        "n_cationicos_permanentes": n_perm,
        "n_sitios_protonaveis": n_strong + n_medium + n_weak,
        "n_aneis_atomos": len(ring_atoms),
    }


def _ring_atoms(nodes, edges) -> set:
    """Átomos que pertencem a algum anel (via arestas que não são pontes)."""
    adj: dict = {n: set() for n in nodes}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    ring = set()
    for a, b in edges:
        # b alcançável a partir de a sem usar a aresta (a,b)?
        seen, stack = {a}, [a]
        found = False
        while stack and not found:
            u = stack.pop()
            for v in adj[u]:
                if u == a and v == b:
                    continue
                if v == b:
                    found = True
                    break
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        if found:
            ring.update((a, b))
    return ring
