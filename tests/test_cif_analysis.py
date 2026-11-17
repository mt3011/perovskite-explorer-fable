"""Análise de CIF: estruturas inorgânicas de referência + CIFs híbridos sintéticos."""
from pathlib import Path

import pytest

from helpers import ISONICOTINIC, NICOTINIC, hybrid_cif, pyrazine_on_inversion_cif
from perovskite_explorer.cif_analysis import analyze_cif

CIFS = Path(__file__).parent / "cifs"


def q(g):
    return {"heavy_elements": g[0], "heavy_edges": g[1]}


@pytest.mark.parametrize("name,dim,sharing,label", [
    ("CsPbBr3_cubica.cif", 3, {"vértices"}, "tipo perovskita"),      # perovskita 3D
    ("Cs2PbI2Cl2_2D.cif", 2, {"vértices"}, "tipo perovskita"),       # Ruddlesden-Popper n=1
    ("CsNiCl3_1D_faces.cif", 1, {"faces"}, "perovskitoide"),          # cadeias de octaedros por faces
    ("Cs2SnI6_0D.cif", 0, set(), "isolados"),                         # octaedros isolados
])
def test_reference_inorganic_structures(name, dim, sharing, label):
    r = analyze_cif((CIFS / name).read_text())
    assert r["ok"] and r["dimensionalidade"] == dim
    assert set(r["compartilhamento"]) == sharing
    assert r["esferas_de_coordenacao"] == ["MX6"]
    assert label in r["classificacao"]
    assert r["cations_inorganicos_sitio_A"] == ["Cs"]


def test_organic_confirmed_by_graph():
    r = analyze_cif(hybrid_cif(2), q(NICOTINIC))
    assert r["organico"]["status"] == "confirmado"
    assert r["dimensionalidade"] == 2 and "tipo perovskita" in r["classificacao"]
    assert not r["organico_alvo_coordenado_ao_metal"]


def test_isomer_is_rejected():
    r = analyze_cif(hybrid_cif(3), q(NICOTINIC))      # cristal contém o isômero isonicotínico
    assert r["organico"]["status"] == "isomero"


def test_disordered_organic_is_inconclusive():
    r = analyze_cif(hybrid_cif(2, disorder_organic=True), q(NICOTINIC))
    assert r["organico"]["status"] == "inconclusivo" and r["desordem"]


def test_organic_coordinated_to_metal_is_flagged():
    r = analyze_cif(hybrid_cif(coordinate_n=True), q(ISONICOTINIC))
    assert r["organico"]["status"] == "confirmado"
    assert r["organico_alvo_coordenado_ao_metal"]
    assert "COORDENADO" in r["classificacao"] and "tipo perovskita" not in r["classificacao"]


def test_molecule_on_special_position_is_completed_by_symmetry():
    pyrazine = (["N", "C", "C", "N", "C", "C"], [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]])
    r = analyze_cif(pyrazine_on_inversion_cif(), q(pyrazine))
    assert r["organico"]["status"] == "confirmado"
    assert "sem metal" in r["classificacao"]


def test_bad_cif_never_raises():
    assert analyze_cif("isto não é um CIF")["ok"] is False
    assert analyze_cif("data_x\n_cell_length_a 5\n")["ok"] is False
