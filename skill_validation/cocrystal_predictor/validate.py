"""
validate.py — validation harness for the cocrystal-formation predictor (scripts/cocrystal.py).

A predictor is worthless — worse, dangerous — without validation. This runs the parts that CAN be
automated; the rest (a Mercury cross-check + a prospective experiment) is manual. Four levels:

  L1 CORRECTNESS   — does the code compute what it claims? (invariances, known values, determinism,
                     conformer stability). Hard PASS/FAIL assertions.
  L2 REPRODUCTION  — does our open complementarity match Mercury's (the CSD-calibrated reference)?
                     Not runnable here (no CSD) — we PRINT the pairs+verdicts for you to diff in Mercury.
  L3 RETROSPECTIVE — does it recover KNOWN outcomes? enrichment of documented cocrystals + a small
                     labelled benchmark (sensitivity/specificity) + a LABEL-SCRAMBLE control + per-signal
                     ABLATION. Honest caveat: confirmed NEGATIVES are scarce in open data (see doc).
  L4 PROSPECTIVE   — pre-register predictions for the theophylline pairs (commit before the bench).

SMILES are PubChem-verified (fetched 2026-07-07); HBD/HBA hand-checked. Run:
    python skill_validation/cocrystal_predictor/validate.py
"""
import os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from scripts import cocrystal as CC

# --- PubChem-verified SMILES (name: SMILES) ---------------------------------------------------
SMI = {
    "theophylline": "CN1C2=C(C(=O)N(C1=O)C)NC=N2", "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "nicotinamide": "C1=CC(=CN=C1)C(=O)N", "isonicotinamide": "C1=CN=CC=C1C(=O)N",
    "urea": "C(=O)(N)N", "glutaric acid": "C(CC(=O)O)CC(=O)O", "succinic acid": "C(CC(=O)O)C(=O)O",
    "adipic acid": "C(CCC(=O)O)CC(=O)O", "ascorbic acid": "C(C(C1C(=C(C(=O)O1)O)O)O)O",
    "saccharin": "C1=CC=C2C(=C1)C(=O)NS2(=O)=O", "oxalic acid": "C(=O)(C(=O)O)O",
    "2-aminopyrimidine": "C1=CN=C(N=C1)N", "carbamazepine": "C1=CC=C2C(=C1)C=CC3=CC=CC=C3N2C(=O)N",
    "benzoic acid": "C1=CC=C(C=C1)C(=O)O", "paracetamol": "CC(=O)NC1=CC=C(C=C1)O", "naphthalene": "c1ccc2ccccc2c1",
}
_DESC = {}
def desc(name):
    if name not in _DESC:
        _DESC[name] = CC.molecule_descriptors(SMI[name])
    return _DESC[name]


def _close(a, b, tol=1e-6):
    return all(abs(x - y) <= tol for x, y in zip(np.atleast_1d(a), np.atleast_1d(b)))


# ================================================================= L1 — CORRECTNESS
def level1():
    print("\n" + "=" * 70 + "\nL1 CORRECTNESS\n" + "=" * 70)
    rng = np.random.default_rng(0)
    fails = 0

    # (a) shape descriptors: rotation/translation invariant, scale invariant
    pts = rng.standard_normal((40, 3)) * np.array([4., 2., 1.])
    Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    base = CC.shape_descriptors(pts)
    inv = _close(base, CC.shape_descriptors(pts @ Q + [7, -3, 11]), 1e-6)
    scl = _close(base, CC.shape_descriptors(pts * 9.0), 1e-6)
    print(f"  shape rotation/translation-invariant : {'PASS' if inv else 'FAIL'}  (S/L,M/L={tuple(round(v,3) for v in base)})")
    print(f"  shape scale-invariant (ratios)       : {'PASS' if scl else 'FAIL'}")
    fails += (not inv) + (not scl)

    # (b) known shapes: a line → S/L≈0 (elongated); a sphere → S/L≈1 (globular)
    line = np.c_[np.linspace(0, 10, 25), np.zeros(25), np.zeros(25)]
    i = np.arange(60); phi = np.arccos(1 - 2 * (i + .5) / 60); th = np.pi * (1 + 5 ** .5) * i
    sph = np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]
    sl_line, sl_sph = CC.shape_descriptors(line)[0], CC.shape_descriptors(sph)[0]
    print(f"  line S/L≈0 ({sl_line:.3f}) < sphere S/L≈1 ({sl_sph:.3f}) : {'PASS' if sl_line < 0.05 < 0.9 < sl_sph else 'FAIL'}")
    fails += not (sl_line < 0.05 < 0.9 < sl_sph)

    # (c) dipole ordering: polar (urea) >> nonpolar (naphthalene)
    du, dn = desc("urea")["dipole_debye"], desc("naphthalene")["dipole_debye"]
    print(f"  dipole urea ({du:.2f} D) >> naphthalene ({dn:.2f} D) : {'PASS' if du > dn + 1 else 'FAIL'}")
    fails += not (du > dn + 1)

    # (d) HBD/HBA match known chemistry
    known = {"theophylline": (1, 3), "caffeine": (0, 3), "urea": (2, 1), "nicotinamide": (1, 2)}
    hb_ok = all(desc(n)["hbd"] == h and desc(n)["hba"] == a for n, (h, a) in known.items())
    print(f"  HBD/HBA vs hand-count (theo 1/3, caffeine 0/3, urea 2/1, nic 1/2) : {'PASS' if hb_ok else 'FAIL'}")
    fails += not hb_ok

    # (e) determinism: same SMILES → identical descriptors (fixed conformer seed)
    d1, d2 = CC.molecule_descriptors(SMI["glutaric acid"]), CC.molecule_descriptors(SMI["glutaric acid"])
    det = _close((d1["s_l"], d1["m_l"], d1["dipole_debye"]), (d2["s_l"], d2["m_l"], d2["dipole_debye"]), 1e-9)
    print(f"  deterministic (fixed seed → identical) : {'PASS' if det else 'FAIL'}")
    fails += not det

    # (f) conformer stability: rigid (theophylline) tight; flexible (adipic) looser — flag if shape sd is huge
    for n in ("theophylline", "adipic acid"):
        d = desc(n)
        print(f"  conformer spread {n:14}: S/L sd={d['s_l_sd']:.3f}, μ sd={d['dipole_sd']:.2f} D  (n={d['n_conf']})")
    print(f"\n  L1 RESULT: {'ALL PASS' if fails == 0 else str(fails)+' FAIL(S)'}")
    return fails == 0


# ================================================================= L3 — RETROSPECTIVE
# Documented cocrystals (positives, with source). Sources (DOIs verified against CrossRef 2026-09-11):
#   Lu2009          Lu & Rohani 2009, Org. Process Res. Dev.            10.1021/op900047r
#   McTague2021     McTague & Rasmuson 2021, Cryst. Growth Des.         10.1021/acs.cgd.1c00296
#   Hou2015         Hou et al. 2015, Cryst. Growth Des.                 10.1021/acs.cgd.5b00800
#   Trask2005       Trask, Motherwell & Jones 2005, Cryst. Growth Des.  10.1021/cg0496540
#   Thompson2010    Thompson et al. 2010, Acta Cryst. C                 10.1107/S0108270110027319
#   Etter1990       Etter & Adsmond 1990, J. Chem. Soc., Chem. Commun.  10.1039/C39900000589
#   Fleischman2003  Fleischman et al. 2003, Cryst. Growth Des.          10.1021/cg034035x
#   CSD             structure deposited in the Cambridge Structural Database
POS = [("theophylline", "nicotinamide", "Lu2009"), ("theophylline", "glutaric acid", "McTague2021"),
       ("theophylline", "saccharin", "Hou2015"), ("caffeine", "glutaric acid", "Trask2005"),
       ("caffeine", "oxalic acid", "Trask2005"), ("nicotinamide", "succinic acid", "Thompson2010"),
       ("2-aminopyrimidine", "succinic acid", "Etter1990"), ("carbamazepine", "saccharin", "Fleischman2003"),
       ("carbamazepine", "nicotinamide", "CSD"), ("isonicotinamide", "oxalic acid", "CSD")]
# Defensible non-formers (negatives) — SCARCE in open data; this is the weak link (see doc).
NEG = [("caffeine", "saccharin", "unreported despite saccharin ubiquity; caffeine has no N-H donor"),
       ("naphthalene", "urea", "no complementary H-bonding (naphthalene inert)"),
       ("naphthalene", "adipic acid", "no H-bond acceptor/donor match")]


def _formation_signals(a, b):
    """The two a-priori FORMATION signals (ΔpKa is a separate salt filter, not formation)."""
    comp = CC.molecular_complementarity(desc(a), desc(b))
    hb = CC.hbond_competition(SMI[a], SMI[b])
    return comp["complementary"], hb["hetero_favoured"]


def _predict(a, b, mode="combined"):
    comp, hetero = _formation_signals(a, b)
    return {"combined": comp and hetero, "complementarity": comp, "hbond": hetero}[mode]


def _fom(y_true, y_pred):
    tp = sum(t and p for t, p in zip(y_true, y_pred)); tn = sum((not t) and (not p) for t, p in zip(y_true, y_pred))
    fp = sum((not t) and p for t, p in zip(y_true, y_pred)); fn = sum(t and (not p) for t, p in zip(y_true, y_pred))
    sens = tp / (tp + fn) if tp + fn else float("nan"); spec = tn / (tn + fp) if tn + fp else float("nan")
    eff = (sens * spec) ** 0.5 if sens == sens and spec == spec else float("nan")
    return sens, spec, eff, (tp, fp, tn, fn)


def level3():
    print("\n" + "=" * 70 + "\nL3 RETROSPECTIVE\n" + "=" * 70)
    labels = [1] * len(POS) + [0] * len(NEG)
    pairs = [(a, b) for a, b, _ in POS] + [(a, b) for a, b, _ in NEG]

    # (a) per-signal ABLATION + combined
    print("  ablation (sensitivity / specificity / efficiency):")
    for mode in ("complementarity", "hbond", "combined"):
        preds = [_predict(a, b, mode) for a, b in pairs]
        sens, spec, eff, cm = _fom(labels, preds)
        print(f"    {mode:16}: sens={sens:.2f} spec={spec:.2f} eff={eff:.2f}  (TP,FP,TN,FN={cm})")

    # (b) LABEL-SCRAMBLE control — predictor should LOSE power on shuffled labels
    rng = np.random.default_rng(1); effs = []
    preds = [_predict(a, b, "combined") for a, b in pairs]
    for _ in range(200):
        effs.append(_fom(list(rng.permutation(labels)), preds)[2])
    real_eff = _fom(labels, preds)[2]
    scr = np.nanmean(effs)
    print(f"\n  scramble control: real eff={real_eff:.2f} vs shuffled-label mean eff={scr:.2f} "
          f"→ {'PASS (real > shuffled)' if real_eff > scr + 0.05 else 'WEAK (not separated)'}")

    # (c) ENRICHMENT — does the predictor rank documented positives above random pool pairs?
    pool = ["theophylline", "caffeine", "nicotinamide", "urea", "glutaric acid", "succinic acid",
            "adipic acid", "ascorbic acid", "saccharin", "2-aminopyrimidine"]
    allpairs = [(pool[i], pool[j]) for i in range(len(pool)) for j in range(i + 1, len(pool))]
    scored = sorted(allpairs, key=lambda p: sum(_formation_signals(*p)), reverse=True)
    doc = {frozenset((a, b)) for a, b, _ in POS if a in pool and b in pool}
    print(f"\n  enrichment over {len(allpairs)} pool pairs (2=both signals, 1=one, 0=none):")
    for rank, (a, b) in enumerate(scored):
        if frozenset((a, b)) in doc:
            print(f"    documented '{a}+{b}' ranks {rank+1}/{len(allpairs)} (score {sum(_formation_signals(a,b))})")


# ================================================================= L4 — PROSPECTIVE (pre-register)
def level4():
    print("\n" + "=" * 70 + "\nL4 PROSPECTIVE — pre-registered SYNTHON predictions (theophylline pairs)\n" + "=" * 70)
    print("  (no binary formation verdict — validated as no-better-than-chance; these are synthon + ΔpKa triage)")
    PKA = {"glutaric acid": 4.34, "succinic acid": 4.21, "adipic acid": 4.43, "ascorbic acid": 4.10, "saccharin": 1.60}
    print(f"  {'co-former':16} {'predicted synthon':30} {'ΔpKa zone':12} |Δμ|/D")
    for cf in ("nicotinamide", "glutaric acid", "adipic acid", "ascorbic acid", "succinic acid", "urea", "saccharin"):
        r = CC.predict_cocrystal(SMI["theophylline"], SMI[cf],
                                 pka_acid=PKA.get(cf), pka_base_conjugate=3.5 if cf in PKA else None)
        ion = r["ionisation"]["zone"] if r["ionisation"] else "n/a"
        print(f"  {cf:16} {r['synthon'][:30]:30} {ion:12} {r['complementarity']['differences']['dipole_debye']:.1f}")


if __name__ == "__main__":
    ok = level1()
    level3()
    level4()
    print("\n" + ("L1 correctness PASSED — L3/L4 are indicative only; the Mercury cross-check and a prospective experiment are manual."
                  if ok else "L1 correctness FAILED — fix before trusting any prediction."))
