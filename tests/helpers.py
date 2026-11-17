"""Utilitários de teste: geração de CIFs sintéticos e um servidor HTTP falso.

ATENÇÃO: tudo aqui é dado SINTÉTICO, criado apenas para exercitar o código sem rede.
Não são resultados reais de nenhum banco de dados.
"""
from __future__ import annotations

import json
import math

from perovskite_explorer.http import SourceError

# grafos de átomos pesados (elementos, arestas) usados como "orgânico procurado"
PYRIDINE_RING = (["N", "C", "C", "C", "C", "C"], [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]])
NICOTINIC = (PYRIDINE_RING[0] + ["C", "O", "O"], PYRIDINE_RING[1] + [[2, 6], [6, 7], [6, 8]])
ISONICOTINIC = (PYRIDINE_RING[0] + ["C", "O", "O"], PYRIDINE_RING[1] + [[3, 6], [6, 7], [6, 8]])


def ring_molecule(center, attach_at=None, n_h_on_n=1, up=None):
    """Anel de 6 membros (N + 5 C) no plano xz, com H; opcionalmente um COOH no átomo `attach_at`.

    O átomo de índice `up` aponta para +z (padrão: o átomo que carrega o COOH), para que a
    molécula se estenda ao longo de z e não colida com a própria imagem periódica em x.
    Retorna lista de (elemento, x, y, z) em Å.
    """
    cx, cy, cz = center
    if up is None:
        up = attach_at or 0
    atoms = []
    for i in range(6):
        ang = math.radians(90 + 60 * (i - up))
        el = "N" if i == 0 else "C"
        atoms.append((el, cx + 1.39 * math.cos(ang), cy, cz + 1.39 * math.sin(ang)))
        if i == attach_at:
            ux, uz = math.cos(ang), math.sin(ang)
            ccx, ccz = cx + 2.89 * ux, cz + 2.89 * uz
            atoms.append(("C", ccx, cy, ccz))
            for sign in (1, -1):
                a2 = ang + sign * math.radians(60)
                atoms.append(("O", ccx + 1.25 * math.cos(a2), cy, ccz + 1.25 * math.sin(a2)))
        elif i != 0 or n_h_on_n:
            atoms.append(("H", cx + 2.45 * math.cos(ang), cy, cz + 2.45 * math.sin(ang)))
    return atoms


def hybrid_cif(attach_at=2, disorder_organic=False, coordinate_n=False):
    """Camada 2D tipo perovskita (PbI4) + anel orgânico entre camadas, em P1.

    attach_at=2 -> ácido nicotínico protonado; attach_at=3 -> isonicotínico (isômero).
    coordinate_n=True -> isonicotínico com o N ligado diretamente ao Pb (no lugar de um iodeto).
    """
    a, c = 6.30, 22.0
    atoms = [("Pb", 0, 0, 0), ("I", a / 2, 0, 0), ("I", 0, a / 2, 0), ("I", 0, 0, 3.15), ("I", 0, 0, -3.15)]
    if coordinate_n:
        atoms = atoms[:3] + [("I", 0, 0, -3.15)]
        center = (0.0, 0.0, 2.6 + 1.39)
        # N (índice 0) aponta para -z, em direção ao Pb; COOH (para) aponta para +z
        mol = ring_molecule(center, 3, n_h_on_n=0, up=3)
    else:
        mol = ring_molecule((a / 2, a / 2, 9.0), attach_at)
    lines = ["data_sintetico", "_chemical_formula_sum 'C6 H6 I4 N O2 Pb'",
             f"_cell_length_a {a}", f"_cell_length_b {a}", f"_cell_length_c {c}",
             "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
             "_symmetry_space_group_name_H-M 'P 1'", "loop_", "_symmetry_equiv_pos_as_xyz", "'x, y, z'",
             "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
             "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy", "_atom_site_disorder_group"]
    for n, (el, x, y, z) in enumerate(atoms + mol):
        organic = n >= len(atoms)
        occ, grp = ("0.5", "1") if (disorder_organic and organic and el != "H" and n % 2) else ("1", ".")
        if disorder_organic and organic and el != "H" and n % 4 == 0:
            continue  # átomos não localizados: o grafo fica incompleto
        lines.append(f"{el}{n} {el} {(x % a) / a:.5f} {(y % a) / a:.5f} {(z % c) / c:.5f} {occ} {grp}")
    return "\n".join(lines) + "\n"


def pyrazine_on_inversion_cif():
    """Pirazina (N2C4) sobre um centro de inversão em P-1: só meia molécula na unidade assimétrica."""
    a = 9.0
    half = []
    for i, el in ((0, "N"), (1, "C"), (2, "C")):
        ang = math.radians(90 + 60 * i)
        half.append((el, 1.39 * math.cos(ang), 0.0, 1.39 * math.sin(ang)))
    lines = ["data_pirazina", f"_cell_length_a {a}", f"_cell_length_b {a}", f"_cell_length_c {a}",
             "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
             "_symmetry_space_group_name_H-M 'P -1'", "loop_", "_symmetry_equiv_pos_as_xyz",
             "'x, y, z'", "'-x, -y, -z'", "loop_", "_atom_site_label", "_atom_site_type_symbol",
             "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
    for n, (el, x, y, z) in enumerate(half):
        lines.append(f"{el}{n} {el} {x / a:.5f} {y / a:.5f} {z / a:.5f}")
    lines.append("Cl1 Cl 0.5 0.5 0.5")
    return "\n".join(lines) + "\n"


class FakeHttp:
    """Substitui HttpClient nos testes: responde por regras (substring da URL -> resposta)."""

    def __init__(self, routes=None, fail=()):
        self.routes = list(routes or [])
        self.fail = tuple(fail)
        self.calls = []
        self.stats = {"rede": 0, "cache": 0, "falhas": 0}

    def _find(self, url, params, data):
        blob = url + " " + json.dumps(params or {}, sort_keys=True) + " " + json.dumps(data or {}, sort_keys=True)
        self.calls.append(blob)
        for bad in self.fail:
            if bad in blob:
                raise SourceError(f"falha simulada: {bad}")
        for key, value in self.routes:
            if key in blob:
                return value(url, params, data) if callable(value) else value
        return None

    def get_json(self, url, params=None, data=None, headers=None, ttl_days=None):
        return self._find(url, params, data)

    def get_text(self, url, params=None, ttl_days=None):
        return self._find(url, params, None)

    def close(self):
        pass
