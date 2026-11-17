"""Lógica das fontes (nomes do CSD, classificação da literatura) e cliente HTTP — sem rede."""
import pytest

from perovskite_explorer import datacite, literature
from perovskite_explorer import structures as st
from perovskite_explorer.http import HttpClient, SourceError

TERMS = {"cations": ["phenazinium", "phenazin-5-ium"], "neutros": ["phenazine"]}
METALS = ["Pb", "Sn", "Cd", "Hg", "Cu", "Sb"]


# ------------------------------------------------------------------ nomes do CSD
@pytest.mark.parametrize("name,metals,halides,catena", [
    ("bis(phenazinium) tetrabromo-cadmium(ii) dihydrate", ["Cd"], ["Br"], False),
    ("catena-[bis(3-carboxyanilinium) bis(mu-chloro)-bis(chloro)-cadmium]", ["Cd"], ["Cl"], True),
    ("catena-(bis(2-phenylethylammonium) tetrakis(mu2-iodo)-lead(ii))", ["Pb"], ["I"], True),
    ("4-chloroanilinium hexachloroantimonate dichloromethane solvate", ["Sb"], ["Cl"], False),
    ("3-carboxyanilinium 2-chloro-4-nitrobenzoate hydrate", [], [], False),          # Cl é substituinte
    ("2,2',3,3',5,5',6,6'-octafluoro-4,4'-diiodo-1,1'-biphenyl phenazine", [], [], False),
    ("3-Carboxyanilinium iodide", [], ["I"], False),
])
def test_csd_name_features(name, metals, halides, catena):
    f = datacite.name_features(name)
    assert (f["metais"], f["haletos"], f["catena"]) == (metals, halides, catena)


def rec(name):
    return {"ccdc": "1", "nome": name, "url": "u", "doi_artigo": "10.1/x", "ano": 2020, "autores": "", "revista": ""}


def test_csd_record_to_structure():
    s = datacite.to_structure(rec("bis(phenazinium) tetrabromo-cadmium(ii) dihydrate"), TERMS, ["Cl", "Br", "I"], METALS)
    assert s["categoria"] == st.HYBRID and s["confianca"] == st.NOME_CSD and s["papel_do_organico"] == st.CATION
    s = datacite.to_structure(rec("catena-(bis(mu2-chloro)-(phenazine)-mercury)"), TERMS, ["Cl", "Br", "I"], METALS)
    assert s["categoria"] == st.HYBRID and "neutra" in s["papel_do_organico"] and "catena" in s["classificacao"]
    s = datacite.to_structure(rec("phenazinium chloride"), TERMS, ["Cl", "Br", "I"], METALS)
    assert s["categoria"] == st.SALT
    # derivado, metal fora da lista-alvo e composto sem haleto são recusados
    assert datacite.to_structure(rec("5-methyl-phenazinium tetrachloro-copper"), TERMS, ["Cl"], METALS) is None
    assert datacite.to_structure(rec("bis(phenazinium) tetrachloro-dioxo-uranium(vi)"), TERMS, ["Cl"], METALS) is None
    assert datacite.to_structure(rec("phenazinium perchlorate"), TERMS, ["Cl"], METALS) is None


def test_unique_samples_do_not_double_count_cod_and_csd():
    def s(fonte, doi, m="Cd", x="Br"):
        return {"fonte": fonte, "doi": doi, "metais": [m], "haletos": [x]}
    strong = [s("COD", "d1"), s("CSD", "d1"), s("CSD", "d1"), s("COD", "d2"), s("CSD", None), s("COD", "d1", "Hg")]
    assert st.unique_samples(strong) == 2 + 1 + 1 + 1


# ------------------------------------------------------------------ literatura
def paper(title, abstract="", **kw):
    return dict({"titulo": title, "resumo": abstract, "doi": None, "autores": [], "fontes": ["x"]}, **kw)


def test_structural_paper():
    p = literature.classify(paper(
        "A Highly Optical Anisotropic Hybrid Perovskite",
        "Here, a new 1D hybrid halide perovskite PZPbCl3 (PZ = phenazine) is reported, which consists of "
        "parallel Pb-Cl chains and PZ cations. Single-crystal X-ray diffraction shows a large band gap."), TERMS)
    assert p["relevante"] and p["uso"] == literature.ESTRUTURAL
    assert "phenazine" in p["trecho"]


def test_additive_paper():
    p = literature.classify(paper(
        "Phenazine additive for efficient perovskite solar cells",
        "Phenazine is introduced as an additive to passivate defects at grain boundaries of the perovskite "
        "film, improving the power conversion efficiency and device stability."), TERMS)
    assert p["relevante"] and p["uso"] == literature.ADITIVO


def test_derivative_and_out_of_context_are_not_evidence():
    p = literature.classify(paper("Dibenzo[a,c]phenazine acceptors for perovskite solar cells"), TERMS)
    assert not p["relevante"] and p["uso"] == literature.NAO_VERIFICAVEL
    p = literature.classify(paper("Phenazine biosynthesis in Pseudomonas", "Phenazine antibiotics are hybrid molecules."), TERMS)
    assert not p["relevante"] and p["uso"] == literature.MENCAO
    p = literature.classify(paper("Safranin sensing", "adsorption of permanganate by phenazine dyes"), TERMS)
    assert not p["relevante"]


def test_merge_by_doi_and_by_title():
    a = paper("Crystallization modulation through optimized SnO2 ETL enables over 20% efficiency", doi="10.1/a", citacoes=3, fontes=["A"])
    b = paper("Crystallization Modulation Through Optimized SnO2 ETL Enables Over 20% Efficiency", doi="10.2/preprint", citacoes=0, fontes=["B"], resumo="texto")
    c = paper("Outro artigo", doi="10.1/a", fontes=["C"])
    m = literature.merge([a, b, c])
    assert len(m) == 1 and m[0]["doi"] == "10.1/a" and m[0]["fontes"] == ["A", "B", "C"] and m[0]["resumo"] == "texto"


def test_independent_groups_by_coauthorship():
    ps = [{"autores": ["John A. Smith", "Mary Jones", "Li Wang"]},
          {"autores": ["Smith JA", "Jones M"]},                 # mesmo grupo, outro formato de nome
          {"autores": ["Wang, Li", "Chen, Yu", "Zhao, Q."]},     # só 1 autor em comum com o 1º -> outro grupo
          {"autores": ["Ana Souza"]}]
    assert literature.independent_groups(ps) == 3


# ------------------------------------------------------------------ HTTP
class Resp:
    def __init__(self, status, text="{}", headers=None):
        self.status_code, self.text, self.headers = status, text, headers or {}


class Session:
    def __init__(self, script):
        self.script, self.n = list(script), 0
        self.headers = {}

    def request(self, *a, **k):
        self.n += 1
        item = self.script.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


@pytest.fixture
def no_sleep(monkeypatch):
    monkeypatch.setattr("perovskite_explorer.http.time.sleep", lambda s: None)


def test_cache_survives_restart(tmp_path, no_sleep):
    db = tmp_path / "c.sqlite"
    s = Session([Resp(200, '{"a": 1}')])
    h = HttpClient(db, session=s)
    assert h.get_json("https://x.org/q", params={"k": 1}) == {"a": 1}
    h.close()
    s2 = Session([])                                   # qualquer acesso à rede falharia
    h2 = HttpClient(db, session=s2)
    assert h2.get_json("https://x.org/q", params={"k": 1}) == {"a": 1}
    assert s2.n == 0 and h2.stats["cache"] == 1


def test_retry_then_success_and_404(tmp_path, no_sleep):
    s = Session([ConnectionError("queda"), Resp(503), Resp(200, "[1]"), Resp(404, "nada")])
    h = HttpClient(tmp_path / "c.sqlite", session=s)
    assert h.get_json("https://x.org/a") == [1] and s.n == 3
    assert h.get_json("https://x.org/b") is None
    assert h.get_json("https://x.org/b") is None and s.n == 4      # 404 também fica no cache


def test_quota_exhausted_fails_fast_and_opens_breaker(tmp_path, no_sleep):
    s = Session([Resp(429, "quota", {"Retry-After": "75000"})])
    h = HttpClient(tmp_path / "c.sqlite", session=s)
    with pytest.raises(SourceError, match="cota"):
        h.get_json("https://api.z.org/a")
    with pytest.raises(SourceError, match="indisponível"):
        h.get_json("https://api.z.org/b")
    assert s.n == 1


def test_offline_mode_uses_only_cache(tmp_path, no_sleep):
    h = HttpClient(tmp_path / "c.sqlite", session=Session([]), offline=True)
    with pytest.raises(SourceError, match="offline"):
        h.get_json("https://x.org/a")


# ------------------------------------------------------------------ conector CSD (API)
def test_csd_smarts_matches_protonated_form_and_rejects_isomer():
    Chem = pytest.importorskip("rdkit.Chem")       # opcional: só roda se o RDKit estiver instalado
    from helpers import NICOTINIC
    from perovskite_explorer.csd import graph_to_smarts

    patt = Chem.MolFromSmarts(graph_to_smarts(*NICOTINIC))
    assert patt is not None
    assert Chem.MolFromSmiles("OC(=O)c1cccnc1").HasSubstructMatch(patt)          # ácido nicotínico
    assert Chem.MolFromSmiles("OC(=O)c1ccc[nH+]c1").HasSubstructMatch(patt)      # forma protonada
    assert not Chem.MolFromSmiles("OC(=O)c1ccncc1").HasSubstructMatch(patt)      # isômero isonicotínico


def test_csd_connector_is_optional():
    from perovskite_explorer import csd

    assert csd.available() in (True, False)        # sem o pacote `ccdc`, apenas fica desativado


# ------------------------------------------------------------------ nível sensível ao papel
def _struct(role, doi, conf=st.CONFIRMADA):
    return {"fonte": "COD", "id": doi, "doi": doi, "confianca": conf, "categoria": st.HYBRID,
            "papel_do_organico": role, "metais": ["Cu"], "haletos": ["Cl"], "composto": doi, "ano": 2020,
            "verificada": True, "analise": {"ok": False}}


def _summary(structs):
    from perovskite_explorer import scoring
    res = {"organico": {"resolvido": True, "n_basicos_medios": 1}, "estruturas": structs, "literatura": [],
           "fontes": {"x": {"ok": True}}, "descartadas": [], "depositos_csd": []}
    return scoring.summarize(res)


def test_ligand_only_use_is_not_counted_as_cation_exploration():
    from perovskite_explorer import scoring
    s = _summary([_struct(st.LIGAND, f"10.1/l{i}") for i in range(12)])
    assert s["nivel_exploracao"] == scoring.LIGANTE and s["n_publicacoes_como_cation"] == 0
    assert s["n_publicacoes_estruturais"] == 12 and s["indice_oportunidade"] == round(80 * 0.85)


def test_cation_use_sets_the_level():
    from perovskite_explorer import scoring
    s = _summary([_struct(st.LIGAND, f"10.1/l{i}") for i in range(12)] + [_struct(st.CATION, "10.1/c1")])
    assert s["nivel_exploracao"] == scoring.POUCO and s["n_publicacoes_como_cation"] == 1
    assert s["n_publicacoes_so_ligante_ou_neutro"] == 12


def test_unverified_candidates_do_not_count_as_evidence():
    from perovskite_explorer import scoring
    pend = _struct(st.UNKNOWN, "10.1/p", conf=st.POSSIVEL)
    pend["verificada"] = False
    s = _summary([pend])
    assert s["nivel_exploracao"] == scoring.INEXPLORADO and s["n_candidatas_nao_verificadas"] == 1


def test_target_metal_and_halide_filter():
    org = {"C": 6, "H": 5, "N": 1, "O": 2}
    f = lambda formula: st._category(formula, org, 1, ["Pb", "Sn"], ["Cl", "Br", "I"])   # noqa: E731
    assert f({"C": 6, "H": 6, "N": 1, "O": 2, "Pb": 1, "I": 3}) == st.HYBRID
    assert f({"C": 6, "H": 6, "N": 1, "O": 2, "La": 1, "Cl": 3}) == "outro"      # metal fora da lista
    assert f({"C": 6, "H": 6, "N": 1, "O": 2, "Sn": 1, "F": 3}) == "outro"       # haleto não pedido
    assert f({"C": 6, "H": 6, "N": 1, "O": 2, "Cl": 1}) == st.SALT
