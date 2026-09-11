"""cocrystal.py — the ΔpKa call, the NNLS sum-of-parents, and the predictor's PURE geometry/logic
(shape descriptors, complementarity). RDKit-dependent parts are smoke-tested only if RDKit is present;
the real predictor validation is skill_validation/cocrystal_predictor/validate.py."""
import numpy as np
import pytest
from scripts import cocrystal as CC


# ---------------------------------------------------------------- ΔpKa rule
def test_delta_pka_zones():
    assert CC.classify_ionisation(pka_acid=4.3, pka_base_conjugate=0.6)["zone"] == "cocrystal"  # Δ<0
    assert CC.classify_ionisation(pka_acid=1.6, pka_base_conjugate=8.0)["zone"] == "salt"        # Δ>3
    assert CC.classify_ionisation(pka_acid=1.6, pka_base_conjugate=3.5)["zone"] == "continuum"   # 0..3
    assert CC.classify_ionisation(1, 1)["confident"] is False


# ---------------------------------------------------------------- shape descriptors (PURE)
def test_shape_descriptors_are_invariant_and_bounded():
    rng = np.random.default_rng(0)
    pts = rng.standard_normal((40, 3)) * np.array([4., 2., 1.])
    base = CC.shape_descriptors(pts)
    Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    # rotation + translation invariant
    assert np.allclose(base, CC.shape_descriptors(pts @ Q + [9, -2, 3]), atol=1e-6)
    # scale invariant (ratios)
    assert np.allclose(base, CC.shape_descriptors(pts * 5.0), atol=1e-6)
    # ordered S/L <= M/L <= 1
    s_l, m_l = base
    assert 0 <= s_l <= m_l <= 1 + 1e-9


def test_shape_line_vs_sphere():
    line = np.c_[np.linspace(0, 10, 25), np.zeros(25), np.zeros(25)]
    i = np.arange(80); phi = np.arccos(1 - 2 * (i + .5) / 80); th = np.pi * (1 + 5 ** .5) * i
    sph = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    assert CC.shape_descriptors(line)[0] < 0.05        # elongated
    assert CC.shape_descriptors(sph)[0] > 0.9          # globular


# ---------------------------------------------------------------- complementarity logic (PURE)
def test_molecular_complementarity_thresholds():
    a = {"s_l": 0.5, "m_l": 0.7, "dipole_debye": 2.0}
    near = {"s_l": 0.55, "m_l": 0.72, "dipole_debye": 3.0}    # within thresholds → complementary
    far = {"s_l": 0.5, "m_l": 0.7, "dipole_debye": 9.0}       # |Δμ|=7 > 4 → not
    assert CC.molecular_complementarity(a, near)["complementary"] is True
    assert CC.molecular_complementarity(a, far)["complementary"] is False


# ---------------------------------------------------------------- RDKit-gated smoke test
rdkit = pytest.importorskip("rdkit", reason="RDKit optional")


def test_predictor_smoke_and_synthon():
    # theophylline + glutaric acid → predicted acid O–H···N(theophylline); no binary verdict emitted
    theo, glut = "CN1C2=C(C(=O)N(C1=O)C)NC=N2", "C(CC(=O)O)CC(=O)O"
    r = CC.predict_cocrystal(theo, glut, pka_acid=4.34, pka_base_conjugate=3.5, n_conf=3)
    assert "likely" not in r                              # honest API: no misleading go/no-go
    assert "acid" in r["synthon"].lower()                 # correct dominant donor
    assert r["ionisation"]["zone"] == "cocrystal"
    d = CC.molecule_descriptors(theo, n_conf=3)
    assert d["hbd"] == 1 and d["hba"] == 3                # theophylline, not caffeine
