"""PubChem (PUG REST): identifica o orgânico e coleta propriedades, sinônimos e grafo molecular.

Documentação da API: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
"""
from __future__ import annotations

import logging
import re

from .chem import HALOGENS, SYMBOLS, describe_organic_graph, heavy, hill_formula, largest_organic_unit, parse_formula
from .http import SourceError

log = logging.getLogger("perovskite_explorer")
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PROPS_NEW = ("Title,MolecularFormula,MolecularWeight,InChIKey,IUPACName,Charge,HBondDonorCount,"
             "HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,XLogP,TPSA,Complexity,"
             "CovalentUnitCount,SMILES,ConnectivitySMILES")
PROPS_OLD = PROPS_NEW.replace("SMILES,ConnectivitySMILES", "IsomericSMILES,CanonicalSMILES")
_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_SMILES_HINT = re.compile(r"[=#\(\)\[\]@]|^[A-Za-z0-9+\-\.\\/]+$")


def _cids_from(obj) -> list:
    try:
        return [int(c) for c in obj["IdentifierList"]["CID"] if int(c) > 0]
    except (KeyError, TypeError, ValueError):
        return []


def find_cids(http, query: str):
    """Resolve o texto do usuário em CIDs. Retorna (lista_de_cids, método_usado)."""
    q = query.strip()
    m = re.fullmatch(r"(?i)cid\s*[:=]?\s*(\d+)", q)
    if m:
        return [int(m.group(1))], "CID"
    if _INCHIKEY.match(q):
        return _cids_from(http.get_json(f"{BASE}/compound/inchikey/{q}/cids/JSON")), "InChIKey"
    m = re.fullmatch(r"(?i)smiles\s*[:=]\s*(\S+)", q)
    if m:
        return _cids_from(http.get_json(f"{BASE}/compound/smiles/cids/JSON", data={"smiles": m.group(1)})), "SMILES"
    # nome (inclui sinônimos e número CAS); POST evita problemas com / e ' no nome
    cids = _cids_from(http.get_json(f"{BASE}/compound/name/cids/JSON", data={"name": q}))
    if cids:
        return cids, "nome"
    if " " not in q and _SMILES_HINT.search(q):
        try:
            cids = _cids_from(http.get_json(f"{BASE}/compound/smiles/cids/JSON", data={"smiles": q}))
        except SourceError:
            cids = []
        if cids:
            return cids, "SMILES"
    return [], "não encontrado"


def _properties(http, cids, props_new=PROPS_NEW, props_old=PROPS_OLD) -> list:
    body = {"cid": ",".join(str(c) for c in cids)}
    try:
        obj = http.get_json(f"{BASE}/compound/cid/property/{props_new}/JSON", data=body)
    except SourceError as exc:
        if "HTTP 400" not in str(exc) or props_old == props_new:
            raise
        obj = http.get_json(f"{BASE}/compound/cid/property/{props_old}/JSON", data=body)
    try:
        return obj["PropertyTable"]["Properties"]
    except (KeyError, TypeError):
        return []


def _synonyms(http, cid) -> list:
    obj = http.get_json(f"{BASE}/compound/cid/{cid}/synonyms/JSON")
    try:
        return list(obj["InformationList"]["Information"][0].get("Synonym", []))
    except (KeyError, TypeError, IndexError):
        return []


def _record_graph(http, cid):
    """Átomos e ligações (com H explícito) do registro 2D do PubChem."""
    obj = http.get_json(f"{BASE}/compound/cid/{cid}/record/JSON")
    comp = obj["PC_Compounds"][0]
    aids = comp["atoms"]["aid"]
    elements = comp["atoms"]["element"]
    charges = {c["aid"]: c["value"] for c in comp["atoms"].get("charge", [])}
    atoms = [{"id": a, "el": SYMBOLS[z] if 0 < z < len(SYMBOLS) else "X", "charge": charges.get(a, 0)}
             for a, z in zip(aids, elements)]
    b = comp.get("bonds", {})
    bonds = list(zip(b.get("aid1", []), b.get("aid2", []), b.get("order", [])))
    return atoms, bonds


def _related_forms(http, cid, organic_formula: dict) -> dict:
    """Formas iônicas e sais do mesmo composto-pai no PubChem (melhor esforço).

    Daqui saem (a) nomes reais do cátion para as buscas e (b) sais de haleto conhecidos,
    que indicam que o precursor (orgânico·HX) existe/é comercial.
    """
    out = {"cations": [], "sais_de_haleto": []}
    obj = http.get_json(f"{BASE}/compound/cid/{cid}/cids/JSON",
                        params={"cids_type": "same_parent_connectivity"})
    cids = [c for c in _cids_from(obj) if c != cid][:300]
    if not cids:
        return out
    org_heavy = heavy(organic_formula)
    cation_cids = []
    for p in _properties(http, cids, "Title,MolecularFormula,Charge", "Title,MolecularFormula,Charge"):
        f = parse_formula(p.get("MolecularFormula", ""))
        hv = heavy(f)
        rest = {e: hv.get(e, 0) - org_heavy.get(e, 0) for e in set(hv) | set(org_heavy)}
        rest = {e: v for e, v in rest.items() if v}
        title = p.get("Title", "") or ""
        if not rest and int(p.get("Charge", 0) or 0) > 0:
            cation_cids.append(p["CID"])
            if title:
                out["cations"].append(title)
        elif rest and all(e in HALOGENS and v > 0 for e, v in rest.items()):
            out["sais_de_haleto"].append({"cid": p["CID"], "nome": title,
                                          "formula": hill_formula(f)})
            first = title.split(" ")[0] if title else ""
            if first.lower().endswith("ium"):
                out["cations"].append(first)
    for c in cation_cids[:2]:
        out["cations"] += [s for s in _synonyms(http, c)[:15] if "ium" in s.lower()]
    out["sais_de_haleto"] = out["sais_de_haleto"][:20]
    return out


def resolve(http, query: str) -> dict:
    """Monta a ficha do orgânico. Sempre retorna um dict; 'resolvido' indica sucesso."""
    org = {"entrada": query, "resolvido": False, "avisos": []}
    cids, how = find_cids(http, query)
    org["metodo_resolucao"] = how
    if not cids:
        org["avisos"].append("Não encontrado no PubChem: sem fórmula/estrutura, só a busca por nome na "
                             "literatura e no COD será feita.")
        return org
    cid = cids[0]
    if len(cids) > 1:
        org["avisos"].append(f"O nome corresponde a {len(cids)} registros do PubChem; usado o primeiro (CID {cid}).")
    props = (_properties(http, [cid]) or [{}])[0]
    org.update({
        "cid": cid,
        "url_pubchem": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
        "titulo": props.get("Title", ""),
        "nome_iupac": props.get("IUPACName", ""),
        "formula_pubchem": props.get("MolecularFormula", ""),
        "massa_molar": _num(props.get("MolecularWeight")),
        "inchikey": props.get("InChIKey", ""),
        "smiles": props.get("SMILES") or props.get("IsomericSMILES") or props.get("ConnectivitySMILES")
        or props.get("CanonicalSMILES") or "",
        "doadores_lig_H": props.get("HBondDonorCount"),
        "aceptores_lig_H": props.get("HBondAcceptorCount"),
        "ligacoes_rotaveis": props.get("RotatableBondCount"),
        "atomos_pesados": props.get("HeavyAtomCount"),
        "xlogp": props.get("XLogP"),
        "tpsa": props.get("TPSA"),
        "unidades_covalentes": props.get("CovalentUnitCount"),
    })
    org["sinonimos"] = _synonyms(http, cid)[:60]

    atoms, bonds = _record_graph(http, cid)
    unit_atoms, unit_bonds, n_units = largest_organic_unit(atoms, bonds)
    if n_units > 1:
        org["avisos"].append("O registro tem mais de um fragmento (sal/mistura); usado apenas o maior "
                             "fragmento orgânico.")
    org.update(describe_organic_graph(unit_atoms, unit_bonds))
    org["formula_hill"] = hill_formula(org["formula"])
    if "C" not in org["formula"]:
        org["avisos"].append("O composto não contém carbono; não é um orgânico.")
    if org["n_sitios_protonaveis"] + org["n_cationicos_permanentes"] == 0:
        org["avisos"].append("Nenhum sítio básico identificado: dificilmente formará cátion por protonação "
                             "(pode atuar como ligante neutro ou co-cristal).")
    try:
        org["formas_relacionadas"] = _related_forms(http, cid, org["formula"])
    except (SourceError, KeyError, TypeError, ValueError) as exc:
        org["formas_relacionadas"] = {"cations": [], "sais_de_haleto": []}
        org["avisos"].append(f"Formas iônicas relacionadas não consultadas ({exc}).")
    org["resolvido"] = True
    return org


def _num(v):
    try:
        return round(float(v), 3)
    except (TypeError, ValueError):
        return None
