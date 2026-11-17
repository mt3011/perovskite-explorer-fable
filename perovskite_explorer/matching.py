"""Casamento por fórmula: a fórmula de um cristal contém o orgânico procurado?

Ideia: fórmula_do_cristal = k x (orgânico ± prótons) + metais + haletos (+ água/solvente/co-cátion).
Isso independe do nome usado pelos autores. Isômeros têm a mesma fórmula, por isso este
teste é só o primeiro filtro; a confirmação vem do nome e, principalmente, do grafo no CIF.
"""
from __future__ import annotations

from .chem import HALOGENS, METALS, heavy

# (rótulo, átomos pesados, nº de H) de espécies pequenas que costumam acompanhar o orgânico.
EXTRAS = [
    ("H2O", {"O": 1}, 2),
    ("metilamônio", {"C": 1, "N": 1}, 6),
    ("formamidínio", {"C": 1, "N": 2}, 5),
    ("amônio", {"N": 1}, 4),
    ("dimetilamônio", {"C": 2, "N": 1}, 8),
    ("DMF", {"C": 3, "N": 1, "O": 1}, 7),
    ("DMSO", {"C": 2, "O": 1, "S": 1}, 6),
    ("metanol", {"C": 1, "O": 1}, 4),
    ("etanol", {"C": 2, "O": 1}, 6),
    ("acetonitrila", {"C": 2, "N": 1}, 3),
    ("acetona", {"C": 3, "O": 1}, 6),
    ("THF", {"C": 4, "O": 1}, 8),
    ("dioxano/GBL", {"C": 4, "O": 2}, 8),
    ("ácido acético/acetato", {"C": 2, "O": 2}, 4),
    ("H3PO2", {"P": 1, "O": 2}, 3),
]
_TOL = 0.26  # tolera fórmulas com ocupações fracionárias pequenas


def _sub(a: dict, b: dict, k: float = 1.0) -> dict:
    out = dict(a)
    for e, v in b.items():
        out[e] = out.get(e, 0) - k * v
    return out


def _is_zero(d: dict) -> bool:
    return all(abs(v) <= _TOL for v in d.values())


def _nonneg(d: dict) -> bool:
    return all(v >= -_TOL for v in d.values())


def _explain_rest(rest: dict):
    """Tenta escrever `rest` (átomos pesados não-metal/não-haleto) como combinação de EXTRAS.

    Retorna lista [(rótulo, n, nH)] ou None. Busca limitada: água + até 2 outras espécies.
    """
    rest = {e: v for e, v in rest.items() if abs(v) > _TOL}
    if not rest:
        return []
    if set(rest) == {"O"}:
        w = round(rest["O"])
        return [("H2O", w, 2)] if 1 <= w <= 12 else None
    others = EXTRAS[1:]
    for i, (lab1, h1, nh1) in enumerate(others):
        for n1 in range(1, 5):
            r1 = _sub(rest, h1, n1)
            if not _nonneg(r1):
                break
            r1c = {e: v for e, v in r1.items() if abs(v) > _TOL}
            if not r1c:
                return [(lab1, n1, nh1)]
            if set(r1c) == {"O"} and 1 <= round(r1c["O"]) <= 8:
                return [(lab1, n1, nh1), ("H2O", round(r1c["O"]), 2)]
            for lab2, h2, nh2 in others[i + 1:]:
                for n2 in range(1, 4):
                    r2 = _sub(r1, h2, n2)
                    if not _nonneg(r2):
                        break
                    r2c = {e: v for e, v in r2.items() if abs(v) > _TOL}
                    if not r2c:
                        return [(lab1, n1, nh1), (lab2, n2, nh2)]
                    if set(r2c) == {"O"} and 1 <= round(r2c["O"]) <= 8:
                        return [(lab1, n1, nh1), (lab2, n2, nh2), ("H2O", round(r2c["O"]), 2)]
    return None


def match_formula(entry_formula: dict, organic: dict):
    """Verifica se `entry_formula` pode conter `organic`.

    organic precisa de: formula (dict com H), h_em_heteroatomos, n_sitios_protonaveis,
    n_cationicos_permanentes.
    Retorna None ou um dict com k (nº de orgânicos por fórmula), extras, nível e balanço de H.
    """
    org = organic["formula"]
    org_heavy = heavy(org)
    if not org_heavy or not entry_formula:
        return None
    for e, v in org_heavy.items():
        if entry_formula.get(e, 0) < v - _TOL:
            return None
    ent_heavy = heavy(entry_formula)
    kmax = int(min((ent_heavy.get(e, 0) + _TOL) // v for e, v in org_heavy.items()))
    h_org = org.get("H", 0)
    h_ent = entry_formula.get("H", 0)
    best = None
    for k in range(kmax, 0, -1):
        res = _sub(ent_heavy, org_heavy, k)
        if not _nonneg(res):
            continue
        rest = {e: v for e, v in res.items() if e not in METALS and e not in HALOGENS}
        extras = _explain_rest(rest)
        if extras is None:
            continue
        cand = {"k": k, "extras": [{"especie": lab, "n": n} for lab, n, _ in extras]}
        n_water = sum(n for lab, n, _ in extras if lab == "H2O")
        h_extras = sum(n * nh for _, n, nh in extras)
        if h_ent > 0:
            p_total = h_ent - h_extras - k * h_org
            p = p_total / k
            lower = -organic.get("h_em_heteroatomos", 0) - 2 * n_water / k - 0.01 \
                - organic.get("n_cationicos_permanentes", 0)
            # o H da água muitas vezes não entra na fórmula declarada -> folga de 2 H por água
            n_sites = organic.get("n_sitios_protonaveis", 0)
            # no máximo um próton por sítio básico (se nenhum sítio foi identificado, tolera 1)
            upper = (n_sites if n_sites else 1) + n_water / k + 0.01
            cand["h_informado"] = True
            cand["protons_por_organico"] = round(p, 2)
            cand["h_coerente"] = lower <= p <= upper
            cand["cationico"] = p > 0.25
            if p > 0.25:
                cand["papel"] = f"cátion protonado (+{p:g} H por molécula)"
            elif p < -0.25:
                cand["papel"] = f"desprotonado/ânion ({p:g} H por molécula)"
            else:
                cand["papel"] = "neutro ou cátion permanente (sem H extra)"
        else:
            cand.update({"h_informado": False, "protons_por_organico": None, "h_coerente": None, "cationico": None,
                         "papel": "H não informado na fórmula"})
        only_water = all(lab == "H2O" for lab, _, _ in extras)
        cand["nivel"] = "exata" if only_water else "com_extras"
        if cand["h_coerente"] is False:
            continue
        # preferimos: sem extras > só água > outros; e, em empate, o maior k
        rank = (0 if not extras else 1 if only_water else 2, -k)
        if best is None or rank < best[0]:
            best = (rank, cand)
    return best[1] if best else None


def composition_key(formula: dict) -> str:
    """Chave de composição reduzida (sem H) para agrupar depósitos do mesmo composto."""
    from math import gcd
    from functools import reduce

    hv = {e: v for e, v in heavy(formula).items()}
    if not hv:
        return ""
    if all(abs(v - round(v)) < 1e-6 for v in hv.values()):
        g = reduce(gcd, (int(round(v)) for v in hv.values()))
        hv = {e: int(round(v)) // g for e, v in hv.items()}
    from .chem import hill_formula
    return hill_formula(hv)
