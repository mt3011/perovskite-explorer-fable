"""Química, nomes e casamento por fórmula (tudo offline)."""
from perovskite_explorer import chem, matching, names

# Registro real do PubChem para a nicotinamida (CID 936): átomos e ligações com H explícito.
NICOTINAMIDE_ELEMENTS = [8, 7, 7, 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1]
NICOTINAMIDE_BONDS = list(zip([1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 6, 7, 7, 9],
                              [8, 6, 9, 8, 14, 15, 5, 6, 8, 7, 10, 11, 9, 12, 13],
                              [2, 1, 2, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1]))

PHENAZINE = {"formula": {"C": 12, "H": 8, "N": 2}, "h_em_heteroatomos": 0, "n_sitios_protonaveis": 2,
             "n_cationicos_permanentes": 0}


def nicotinamide():
    atoms = [{"id": i + 1, "el": chem.SYMBOLS[z], "charge": 0} for i, z in enumerate(NICOTINAMIDE_ELEMENTS)]
    return chem.describe_organic_graph(atoms, NICOTINAMIDE_BONDS)


def test_parse_formula_variants():
    assert chem.parse_formula("- C12 H9 Cl4 N2 Sb -") == {"C": 12, "H": 9, "Cl": 4, "N": 2, "Sb": 1}
    assert chem.parse_formula("C6H6N2O") == {"C": 6, "H": 6, "N": 2, "O": 1}
    assert chem.parse_formula("C12H9N2+") == {"C": 12, "H": 9, "N": 2}
    assert chem.parse_formula("C2 D6 N")["H"] == 6
    assert chem.parse_formula("C12 H8.3 Br2")["H"] == 8.3
    assert chem.hill_formula({"N": 2, "H": 8, "C": 12}) == "C12H8N2"


def test_nitrogen_sites_of_nicotinamide():
    d = nicotinamide()
    assert d["formula"] == {"C": 6, "H": 6, "N": 2, "O": 1}
    kinds = sorted(s["basicidade"] for s in d["sitios_N"])
    assert kinds == ["médio", "não"]          # N piridínico protona; N da amida não
    assert d["n_sitios_protonaveis"] == 1
    assert len(d["heavy_elements"]) == 9 and len(d["heavy_edges"]) == 9


def test_salt_input_keeps_only_the_organic():
    atoms = [{"id": 1, "el": "N", "charge": 1}, {"id": 2, "el": "C", "charge": 0},
             {"id": 3, "el": "Cl", "charge": -1}]
    unit, bonds, n = chem.largest_organic_unit(atoms, [(1, 2, 1)])
    assert n == 2 and {a["el"] for a in unit} == {"N", "C"}


def test_cation_name_variants():
    assert "phenazinium" in names.cation_variants("Phenazine")
    assert "phenethylammonium" in names.cation_variants("phenethylamine")
    assert "nicotinamidium" in names.cation_variants("Nicotinamide")
    assert "3-carboxyanilinium" in names.cation_variants("3-Aminobenzoic acid")
    assert "3-carboxypyridinium" in names.cation_variants("pyridine-3-carboxylic acid")
    assert "imidazolium" in names.cation_variants("imidazole")
    assert names.cation_variants("pyridinium") == []


def test_synonym_cleaning_drops_registry_codes():
    syn = ["Dapsone", "80-08-0", "4,4'-Sulfonyldianiline", "DDS", "UNII-8W5C518302", "NSC 6091", "CHEBI:4325"]
    assert names.clean_synonyms(syn) == ["Dapsone", "4,4'-Sulfonyldianiline"]


def test_term_search_ignores_derivatives():
    t = ["phenazine", "phenazinium", "phenazin-5-ium"]
    assert names.find_terms("PZPbCl3 (PZ = phenazine) is reported", t) == ["phenazine"]
    assert names.find_terms("Dibenzo[a,c]phenazine acceptors", t) == []
    assert names.find_terms("N-methylphenazinium salts", t) == []
    assert names.find_terms("3,7-bis(dimethylamino)phenazin-5-ium chloride", t) == []
    assert names.find_terms("Phenazinium lead iodide", t) == ["phenazinium"]


def test_strict_chemical_name_matching():
    f = names.find_terms_in_chemname
    t = ["phenazinium", "phenazine", "pyridinium"]
    assert f("bis(phenazinium) tetrabromo-cadmium(ii) dihydrate", t) == ["phenazinium"]
    assert f("catena-[(mu-phenazine)-dichloro-copper]", t) == ["phenazine"]
    assert f("2-Amino-3-phenylamino-5-phenyl-5-phenazinium perchlorate", t) == []
    assert f("tris(mu-pyridinium-3-carboxylato)-diaqua-holmium", t) == []


def test_formula_match_cationic_hydrate():
    m = matching.match_formula(chem.parse_formula("C24 H22 Cl4 Hg N4 O2"), PHENAZINE)
    assert m["k"] == 2 and m["nivel"] == "exata" and m["cationico"] is True
    assert m["extras"] == [{"especie": "H2O", "n": 2}]
    assert m["protons_por_organico"] == 1.0


def test_formula_match_neutral_ligand_and_rejections():
    m = matching.match_formula(chem.parse_formula("C12 H8 Cd Cl2 N2"), PHENAZINE)
    assert m["k"] == 1 and m["cationico"] is False
    # H demais para ser fenazina (ex.: benzidina C12H12N2 + 2 H): balanço de H incoerente
    assert matching.match_formula(chem.parse_formula("C12 H14 Cl4 N2 Pb"), PHENAZINE) is None
    # carbono insuficiente
    assert matching.match_formula(chem.parse_formula("C10 H10 I4 N2 Pb"), PHENAZINE) is None


def test_formula_match_with_solvent_is_lower_tier():
    m = matching.match_formula(chem.parse_formula("C15 H16 Br3 N3 O Pb"), PHENAZINE)   # + DMF
    assert m["nivel"] == "com_extras" and m["extras"][0]["especie"] == "DMF"


def test_composition_key_reduces_and_drops_h():
    assert matching.composition_key({"C": 24, "H": 18, "N": 4, "Hg": 2, "I": 8}) == "C12Hg1I4N2".replace("Hg1", "Hg")
