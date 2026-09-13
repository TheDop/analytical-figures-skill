"""One bond / H-bond predicate for the table AND the figures, and a PXRD cache that never
re-parses the CIF.

What these pin (ROADMAP 'From the pre-publication review', deferred items 1, 2 and 4):
- crystal_view used to re-implement the engine's bond and H-bond criteria and omitted the
  disorder-alternative exclusion, so a figure could dash a contact the validation table
  suppresses. Both now go through crystal_engine.is_bonded / is_hbond, `alt` included.
- The KD-tree searches must return exactly what the all-pairs loops returned.
- calc_pattern / reflection_list / peak_table on one structure parse the CIF once.
"""
import glob
import os

import numpy as np
import pytest

from scripts import crystal_engine as ce, crystal_view as cv, crystal_pxrd as cp
from scripts.config import Config

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIFS = os.path.join(ROOT, "skill_validation", "crystal", "cifs")
DISORDER = ["disorder-groups_flufenamic__COD4118081.cif", "zprime9-densflag_flufenamic__COD4118079.cif"]

gemmi = pytest.importorskip("gemmi")


def _load(name, **kw):
    cfg = Config(cif_path=os.path.join(CIFS, name), strict=False, **kw)
    return ce.load(cfg), cfg


# ------------------------------------------------------------ the predicate, in isolation
def _atom(label, sym, xyz, occ=1.0):
    return {"label": label, "sym": sym, "xyz": np.array(xyz, float), "occ": occ}


def test_view_suppresses_a_donor_own_alternative_acceptor_like_the_table():
    """A textbook O-H...O geometry whose acceptor is a disorder alternative of the donor: the
    engine reports it as 'disorder' (suppressed), and the view's H-bond search drops it once
    `alt` is passed — the exact contact the old figure path drew and the table did not list."""
    cfg = Config(strict=False)
    D = _atom("O1", "O", (0, 0, 0), occ=0.5)
    H = _atom("H1", "H", (0.97, 0, 0), occ=0.5)
    A = _atom("O1B", "O", (2.8, 0, 0), occ=0.5)
    cluster = [D, H, A]
    assert ce.hbond_geometry(D, H, A, cfg)["kind"] == "hbond"
    assert cv._hbond_pairs(cluster, cfg) == [(1, 2)]
    alt = {"O1": {"O1B"}, "O1B": {"O1"}, "H1": set()}
    assert ce.hbond_geometry(D, H, A, cfg, alt)["kind"] == "disorder"
    assert ce.is_hbond(D, H, A, cfg, alt) is False
    assert cv._hbond_pairs(cluster, cfg, alt) == []


def test_is_bonded_matches_the_engine_rules():
    assert ce.is_bonded("C", "C", 1.54)
    assert not ce.is_bonded("C", "C", 0.3)                       # coincident
    assert not ce.is_bonded("C", "C", 2.5)                       # beyond radii + tolerance
    assert not ce.is_bonded("O", "O", 0.6, 0.5, 0.5)             # near-coincident partial sites
    assert ce.is_bonded("O", "O", 0.6, 1.0, 0.5)                 # only one partial -> ordinary rule
    alt = {"C1": {"C1B"}}
    assert not ce.is_bonded("C", "C", 1.54, 1, 1, "C1", "C1B", alt)
    assert ce.is_bonded("C", "C", 1.54, 1, 1, "C1", "C2", alt)


# ------------------------------------------------------------ real CIFs: view <= table
@pytest.mark.parametrize("name", DISORDER)
def test_view_hbonds_are_a_subset_of_the_table_on_disorder_cifs(name):
    """Every dashed D-H...A in the H-bond-environment figure must be a row of the validation
    table (same labels), and none may be a pair the table suppressed."""
    s, cfg = _load(name)
    hb = ce.hbonds(s, cfg)
    table = {(r["D"], r["H"], r["A"]) for r in hb["bonds"]}
    suppressed = set(hb["suppressed"])
    alt = ce.disorder_alternatives(s)
    env = cv._hbond_env_atoms(s, Config(cif_path=cfg.cif_path, strict=False, hbond_neighbour="whole"))
    drawn = {(env[d]["label"], env[h]["label"], env[a]["label"]) for h, d, a in cv._hbond_records(env, cfg, alt)}
    assert drawn, "expected at least one H-bond in the environment figure"
    assert drawn <= table, f"figure asserts contacts the table does not list: {drawn - table}"
    assert not {(D, A) for D, _H, A in drawn} & suppressed


def test_kdtree_bond_search_equals_all_pairs():
    """The vectorised search returns exactly the pairs the O(N^2) loop did (one CIF, whole
    cluster; the full 9-CIF snapshot diff was run when this landed)."""
    s, cfg = _load("chiral_L-alanine__COD1574526.cif")
    atoms = cv.complete_molecules(s, cfg)
    alt = ce.disorder_alternatives(s)
    fast = set(cv._bond_pairs(atoms, alt))
    slow = set()
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            d = float(np.linalg.norm(atoms[i]["xyz"] - atoms[j]["xyz"]))
            if ce.is_bonded(atoms[i]["sym"], atoms[j]["sym"], d, atoms[i]["occ"], atoms[j]["occ"],
                            atoms[i]["label"], atoms[j]["label"], alt):
                slow.add((i, j))
    assert fast == slow and len(fast) > 5


# ------------------------------------------------------------ PXRD cache
def test_calc_pattern_parses_the_cif_once_per_structure(monkeypatch):
    dif = pytest.importorskip("Dans_Diffraction")
    s, cfg = _load("chiral_L-alanine__COD1574526.cif", pxrd_wavelength=1.5406, pxrd_crosscheck="bragg")
    calls = []
    real = dif.Crystal

    def counting(path):
        calls.append(path)
        return real(path)

    monkeypatch.setattr(dif, "Crystal", counting)
    p1 = cp.calc_pattern(s, cfg)
    p2 = cp.calc_pattern(s, cfg)
    refl = cp.reflection_list(s, cfg)
    rows = cp.peak_table(s, cfg)
    assert len(calls) == 1, f"CIF parsed {len(calls)} times"
    assert np.array_equal(p1["intensity"], p2["intensity"]) and p1["intensity"] is not p2["intensity"]
    assert refl and rows
    # a different wavelength is a different key: recomputed, not served from the cache
    p3 = cp.calc_pattern(s, Config(cif_path=cfg.cif_path, strict=False, pxrd_wavelength=0.71073,
                                   pxrd_crosscheck="bragg"))
    assert p3["peaks"][0] != p1["peaks"][0]
    assert len(calls) == 1
