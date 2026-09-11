"""sources.py — the PURE URL-builders and response parsers (no network in the gate). COD is
verified live end-to-end by hand (search → fetch → crystal_engine.load); these lock the parsing
against the REAL response shapes so a schema drift fails the gate. PubChem's host is DNS-blocked
in this sandbox, so only its offline logic is tested here."""
import json
from scripts import sources as src


# ---------------------------------------------------------------- COD
def test_cod_search_url_builds_formula_and_elements():
    url = src.cod_search_url(formula="C9 H8 O4", elements=["C", "O"], exclude=["N"])
    assert url.startswith("https://www.crystallography.net/cod/result?")
    assert "format=json" in url
    assert "formula=C9+H8+O4" in url          # urlencode → spaces as '+'
    assert "el1=C" in url and "el2=O" in url and "nel1=N" in url


def test_cod_cif_url():
    assert src.cod_cif_url(7247819) == "https://www.crystallography.net/cod/7247819.cif"


def test_cod_parse_results_real_shape():
    # the observed COD JSON: a flat list keyed by `file`; formula may carry '- ... -' decorators
    payload = json.dumps([
        {"file": "1515581", "formula": "- C9 H8 O4 -", "sg": "P 1 21/c 1",
         "a": "12.1016", "b": "6.4721", "c": "11.3344", "vol": "825.44", "doi": "10.1039/c1sc00430a"},
        {"file": "7247819", "formula": "C9 H8 O4", "sg": "P 1 21/c 1"},
        {"no_file_key": True},                # skipped (no COD id)
    ])
    hits = src.cod_parse_results(payload)
    assert len(hits) == 2
    assert hits[0]["id"] == "1515581" and hits[0]["formula"] == "C9 H8 O4"   # dashes stripped
    assert hits[0]["doi"] == "10.1039/c1sc00430a"
    assert hits[1]["id"] == "7247819"


# ---------------------------------------------------------------- COD identity (anti-collision guard)
# A minimal but realistic CIF fragment: the SAME formula (C5 H8 O4) can be glutaric acid OR an
# oxalate ester — identity lives in the name, never the formula (the real bug this guards against).
_CIF_ESTER = """data_2105040
_chemical_name_common            'Methyl ethyl oxalate'
_chemical_formula_moiety         'C5 H8 O4'
_symmetry_space_group_name_H-M   'C 1 2/c 1'
_journal_paper_doi               10.1107/S0108768111037487
"""
_CIF_ACID = """data_4116142
_chemical_name_systematic
;
 Glutaric acid
;
_chemical_formula_moiety         'C5 H8 O4'
_space_group_name_H-M_alt        'C 1 2/c 1'
"""


def test_cif_identity_parses_bare_quoted_and_text_block():
    ester = src.cif_identity(_CIF_ESTER)
    assert ester["name"] == "Methyl ethyl oxalate"
    assert ester["moiety"] == "C5 H8 O4" and ester["doi"] == "10.1107/S0108768111037487"
    acid = src.cif_identity(_CIF_ACID)             # ';'-delimited text-block name
    assert acid["name"] == "Glutaric acid"
    assert acid["spacegroup"] == "C 1 2/c 1"


def test_identity_matches_is_the_collision_guard():
    ester, acid = src.cif_identity(_CIF_ESTER), src.cif_identity(_CIF_ACID)
    # same formula, opposite identity: the guard must tell them apart
    assert src.identity_matches(acid, "glutaric") is True
    assert src.identity_matches(ester, "glutaric") is False
    assert src.identity_matches(ester, "oxalate") is True


def test_cod_search_formula_carries_a_verify_warning(monkeypatch):
    # a formula search must flag that a formula match is not an identity match
    monkeypatch.setattr(src, "_get", lambda *a, **k: json.dumps([{"file": "2105040", "formula": "C5 H8 O4"}]))
    res = src.cod_search(formula="C5 H8 O4")
    assert res["n"] == 1 and "warning" in res and "identity" in res["warning"].lower()
    # a text search makes no formula-collision claim → no warning
    assert "warning" not in src.cod_search(text="glutaric acid")


# ---------------------------------------------------------------- PubChem
def test_pubchem_property_url_default_and_namespace():
    url = src.pubchem_property_url("aspirin")
    assert url == ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/"
                   "property/MolecularFormula,MolecularWeight,XLogP,InChIKey/JSON")
    # smiles omitted from defaults (rename-safety); explicit request + non-name namespace
    u2 = src.pubchem_property_url("BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
                                  props=("ConnectivitySMILES",), namespace="inchikey")
    assert "/compound/inchikey/" in u2 and u2.endswith("property/ConnectivitySMILES/JSON")


def test_pubchem_parse_properties_real_shape():
    payload = json.dumps({"PropertyTable": {"Properties": [
        {"CID": 2244, "MolecularFormula": "C9H8O4", "MolecularWeight": "180.16",
         "XLogP": 1.2, "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"}]}})
    p = src.pubchem_parse_properties(payload)
    assert p["MolecularFormula"] == "C9H8O4" and p["InChIKey"].startswith("BSYNRYMUTXBXSQ")
    assert src.pubchem_parse_properties('{"PropertyTable":{"Properties":[]}}') == {}
