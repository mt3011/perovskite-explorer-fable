"""Pipeline completo com um servidor falso (dados SINTÉTICOS): resultado, checkpoints e retomada."""
import json

import pytest

from helpers import FakeHttp, hybrid_cif
from perovskite_explorer import literature, pipeline, scoring
from perovskite_explorer import structures as st

# ácido nicotínico: N1 piridínico, COOH no C3 (ordens de ligação kekulizadas, H explícito)
ELEMENTS = [7, 6, 6, 6, 6, 6, 6, 8, 8, 1, 1, 1, 1, 1]
BONDS = [(1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1), (5, 6, 2), (6, 1, 1), (3, 7, 1), (7, 8, 2), (7, 9, 1),
         (2, 10, 1), (4, 11, 1), (5, 12, 1), (6, 13, 1), (9, 14, 1)]
RECORD = {"PC_Compounds": [{"atoms": {"aid": list(range(1, 15)), "element": ELEMENTS},
                            "bonds": {"aid1": [b[0] for b in BONDS], "aid2": [b[1] for b in BONDS],
                                      "order": [b[2] for b in BONDS]}}]}
PROPS = {"PropertyTable": {"Properties": [{"CID": 938, "Title": "Nicotinic acid", "MolecularFormula": "C6H5NO2",
                                           "MolecularWeight": "123.11", "IUPACName": "pyridine-3-carboxylic acid",
                                           "SMILES": "C1=CC(=CN=C1)C(=O)O", "InChIKey": "PVNIIMVLHYAWGP-UHFFFAOYSA-N"}]}}


def cod_entry(i, formula, name="", doi="10.1000/estrutura"):
    return {"file": str(i), "formula": f"- {formula} -", "chemname": name, "title": "Synthetic entry", "year": "2021",
            "authors": "Silva, A.; Souza, B.", "journal": "J. Test", "doi": doi, "sg": "P 1", "a": "6.3"}


def cod_result(url, params, data):
    els = {v for k, v in (params or {}).items() if k.startswith("el")}
    if {"Pb", "I"} <= els:
        return [cod_entry(1, "C6 H6 I4 N O2 Pb"),            # contém o orgânico (grafo confirma)
                cod_entry(2, "C6 H6 I4 N O2 Pb"),            # mesmo fórmula, mas é o isômero
                cod_entry(3, "C10 H10 I4 N2 Pb")]            # fórmula incompatível
    return []


def datacite_item(ccdc, name, doi="10.1000/estrutura"):
    return {"attributes": {"doi": f"10.5517/ccdc.csd.{ccdc}", "publicationYear": 2021,
                           "titles": [{"title": f"CCDC {ccdc}: Experimental Crystal Structure Determination"}],
                           "subjects": [{"subject": "Crystal Structure"}, {"subject": name}],
                           "descriptions": [{"description": f"Related Article: A. Silva, B. Souza|2021|J.Test|1|1|doi:{doi}"}],
                           "relatedIdentifiers": [{"relatedIdentifier": doi, "relatedIdentifierType": "DOI"}]}}


def datacite(url, params, data):
    q = params["query"]
    if "relatedIdentifiers" in q:
        return {"data": [datacite_item(1000111, "catena-[bis(3-carboxypyridinium) tetrakis(mu-iodo)-lead(ii)]"),
                         datacite_item(1000112, "bis(acridinium) tetraiodo-lead")], "meta": {"total": 2}}
    if " AND " in q:
        return {"data": [datacite_item(1000111, "catena-[bis(3-carboxypyridinium) tetrakis(mu-iodo)-lead(ii)]"),
                         datacite_item(1000113, "3-carboxypyridinium 2-chloro-benzoate")], "meta": {"total": 2}}
    return {"data": [], "meta": {"total": 7}}


EPMC = {"resultList": {"result": [
    {"doi": "10.1000/ESTRUTURA", "title": "A layered lead iodide perovskite templated by 3-carboxypyridinium",
     "abstractText": "We report the crystal structure of a new 2D hybrid perovskite with 3-carboxypyridinium cations "
                     "between inorganic layers of corner-sharing octahedra.",
     "pubYear": "2021", "citedByCount": 4, "authorString": "Silva A, Souza B."},
    {"doi": "10.1000/aditivo", "title": "Nicotinic acid as additive in perovskite solar cells",
     "abstractText": "Nicotinic acid passivates defects at the perovskite interface and raises the power "
                     "conversion efficiency of the solar cells.",
     "pubYear": "2023", "citedByCount": 10, "authorString": "Lima C, Rocha D."}]}}

ROUTES = [
    ("compound/name/cids", {"IdentifierList": {"CID": [938]}}),
    ("cids_type", None),
    ("/property/", PROPS),
    ("/synonyms/", {"InformationList": {"Information": [{"CID": 938, "Synonym": ["Nicotinic acid", "59-67-6", "Niacin"]}]}}),
    ("/record/", RECORD),
    ("cod/result", cod_result),
    ("cod/1.cif", hybrid_cif(2)),
    ("cod/2.cif", hybrid_cif(3)),
    ("api.datacite.org", datacite),
    ("europepmc", EPMC),
    ("api.openalex.org", {"results": []}),
    ("api.crossref.org", {"message": {"items": []}}),
]
CFG = dict(pipeline.DEFAULT_CONFIG, metais=["Pb", "Sn"], halogenios=["I"], openalex_api_key="chave-de-teste")


def test_end_to_end_result(tmp_path):
    http = FakeHttp(ROUTES)
    res = pipeline.run(["Nicotinic acid"], tmp_path, http, CFG)[0]
    org = res["organico"]
    assert org["resolvido"] and org["formula_hill"] == "C6H5NO2" and org["n_sitios_protonaveis"] == 1
    assert "3-carboxypyridinium" in res["termos"]["cations"]

    by_id = {s["id"]: s for s in res["estruturas"]}
    assert by_id["1"]["confianca"] == st.CONFIRMADA and by_id["1"]["papel_do_organico"] == st.CATION
    assert "tipo perovskita" in by_id["1"]["classificacao"]
    assert by_id["CCDC 1000111"]["confianca"] == st.NOME_CSD
    assert [d["id"] for d in res["descartadas"]] == ["2"]           # o isômero foi rejeitado pelo grafo
    assert "3" not in by_id and "CCDC 1000113" not in by_id              # fórmula incompatível / Cl de substituinte

    s = res["resumo"]
    assert s["nivel_exploracao"] == scoring.POUCO and s["busca_completa"]
    assert s["n_estruturas_hibridas_fortes"] == 1                    # COD 1 e CCDC 111 são a mesma amostra
    assert s["n_estruturas_cod"] == 1 and s["n_estruturas_csd"] == 1
    assert s["n_publicacoes_estruturais"] == 1 and s["n_artigos_aditivo"] == 1
    assert s["combinacoes_MX"] == {"Pb–I": 1} and s["dimensionalidades"] == {"2D": 1}
    uses = {p["doi"]: p["uso"] for p in res["literatura"]}
    assert uses == {"10.1000/estrutura": literature.ESTRUTURAL, "10.1000/aditivo": literature.ADITIVO}
    deps = {d["ccdc"]: d["nome_cita_o_organico"] for d in res["depositos_csd"]}
    assert deps == {"1000111": True, "1000112": False}

    for f in ("resultados.xlsx", "relatorio.html", "resumo.csv", "estruturas.csv", "literatura.csv",
              "dados_completos.json", "estado.json", "checkpoints/nicotinic_acid.json"):
        assert (tmp_path / f).stat().st_size > 0, f
    assert "Nicotinic acid" in (tmp_path / "relatorio.html").read_text(encoding="utf-8")


def test_resume_skips_finished_and_retries_failed_sources(tmp_path):
    names = ["Nicotinic acid", "Pyridine-3-carboxylic acid"]
    broken = FakeHttp(ROUTES, fail=["europepmc"])
    first = pipeline.run(names, tmp_path, broken, CFG)
    assert all(not r["completo"] for r in first)
    assert all("Europe PMC" in r["resumo"]["fontes_com_falha"] for r in first)
    assert first[0]["resumo"]["n_estruturas_cod"] == 1               # as outras fontes seguiram funcionando

    fixed = FakeHttp(ROUTES)
    second = pipeline.run(names, tmp_path, fixed, CFG)
    assert all(r["completo"] for r in second) and fixed.calls           # incompletos foram refeitos

    idle = FakeHttp(ROUTES)
    third = pipeline.run(names, tmp_path, idle, CFG)
    assert idle.calls == [] and len(third) == 2                         # nada é refeito depois de completo

    changed = FakeHttp(ROUTES)
    pipeline.run(names, tmp_path, changed, dict(CFG, halogenios=["I", "Br"]))
    assert changed.calls                                                # parâmetros mudaram -> refaz


def test_interrupted_run_keeps_previous_organics(tmp_path, monkeypatch):
    real = pipeline.process_organic
    calls = []

    def flaky(http, name, cfg):
        calls.append(name)
        if name == "B":
            raise KeyboardInterrupt
        return real(http, name, cfg)

    monkeypatch.setattr(pipeline, "process_organic", flaky)
    done = pipeline.run(["Nicotinic acid", "B", "C"], tmp_path, FakeHttp(ROUTES), CFG)
    assert [r["entrada"] for r in done] == ["Nicotinic acid"] and calls == ["Nicotinic acid", "B"]
    assert json.loads((tmp_path / "estado.json").read_text())["interrompido"] is True
    assert "Execução parcial" in (tmp_path / "relatorio.html").read_text(encoding="utf-8")

    monkeypatch.setattr(pipeline, "process_organic", real)
    again = FakeHttp(ROUTES)
    done = pipeline.run(["Nicotinic acid", "B", "C"], tmp_path, again, CFG)
    assert len(done) == 3
    assert not any("Nicotinic" in c and "name/cids" in c for c in again.calls[:1])   # o 1º não foi refeito


def test_internal_error_in_one_organic_does_not_stop_the_list(tmp_path, monkeypatch):
    real = pipeline.pubchem.resolve

    def boom(http, name):
        if name == "ruim":
            raise RuntimeError("bug simulado")
        return real(http, name)

    monkeypatch.setattr(pipeline.pubchem, "resolve", boom)
    done = pipeline.run(["ruim", "Nicotinic acid"], tmp_path, FakeHttp(ROUTES), CFG)
    assert done[0]["erro_interno"].startswith("RuntimeError") and not done[0]["resumo"]["busca_completa"]
    assert done[1]["resumo"]["nivel_exploracao"] == scoring.POUCO


def test_unresolved_name_still_searches_by_name(tmp_path):
    routes = [("compound/name/cids", None)] + ROUTES[1:]
    res = pipeline.run(["Composto Inexistente XYZ"], tmp_path, FakeHttp(routes), CFG)[0]
    assert not res["organico"]["resolvido"]
    assert res["resumo"]["nivel_exploracao"] == scoring.INEXPLORADO
    assert res["resumo"]["fator_viabilidade"] == 0.5


def test_read_list_ignores_comments_and_duplicates(tmp_path):
    f = tmp_path / "l.txt"
    f.write_text("# comentário\nPhenazine\n\nphenazine\nNiacinamide  # nota\n", encoding="utf-8")
    assert pipeline.read_list(f) == ["Phenazine", "Niacinamide"]


def test_openalex_without_key_is_skipped_not_failed(tmp_path):
    http = FakeHttp(ROUTES)
    res = pipeline.run(["Nicotinic acid"], tmp_path, http, dict(CFG, openalex_api_key=None))[0]
    assert res["fontes"]["OpenAlex"]["desativada"] and res["completo"]
    assert not any("openalex" in c for c in http.calls)
