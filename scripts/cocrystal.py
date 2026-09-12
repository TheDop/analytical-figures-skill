"""
cocrystal.py — quantitative helpers for cocrystal IDENTIFICATION.

The DECISION layer is `references/cocrystal_id.md` (which evidence settles cocrystal vs salt
vs physical mixture vs polymorph); THIS module is the arithmetic that doc routes to. Built in
waves:

  wave 1:
    delta_pka / classify_ionisation      — the ΔpKa salt / cocrystal / continuum call [zero-dep]
    rwp / sum_of_parents / phase_report  — NNLS "new phase vs physical mixture"        [scipy]
  next:
    class_fom                             — one-class / classification figures of merit [numpy]

Nothing here decides the science. ΔpKa *ranks salt-risk*; it does not prove protonation — the
0–3 continuum is explicitly unpredictable, confirm by FTIR (the ionised-group band) or ssNMR
(see `cocrystal_id.md`). Every function returns numbers + a provenance caption, never a bare
label; `classify_ionisation` NEVER raises.

Validation: `skill_validation/cocrystal/` (closed-form + reference parity), bridged into
`python -m pytest tests/`.
"""
from __future__ import annotations
import numpy as np


# ============================================================= ΔpKa (salt vs cocrystal)
# Cruz-Cabeza 2012 (CrystEngComm 14:6362), the "pKa rule". ΔpKa is defined with the BASE's
# conjugate-acid pKa minus the ACID's pKa — both the pKa of the PROTONATED species:
#     ΔpKa = pKa(BH+) − pKa(AH)
# Simplified 3-zone rule (thresholds configurable): ΔpKa > 3 → salt, < 0 → cocrystal, else the
# salt–cocrystal continuum (proton transfer unpredictable — must be confirmed experimentally).
SALT_THRESHOLD = 3.0
COCRYSTAL_THRESHOLD = 0.0


def delta_pka(pka_acid, pka_base_conjugate):
    """ΔpKa for the pKa rule = pKa(base's conjugate acid, BH+) − pKa(acid, AH). Both inputs are
    the pKa of the PROTONATED species; returns the float. (The pKa values are look-ups, not
    computed here — fetch + cite them per the web-search discipline in SKILL.md.)"""
    return float(pka_base_conjugate) - float(pka_acid)


def classify_ionisation(pka_acid, pka_base_conjugate,
                        salt_threshold=SALT_THRESHOLD, cocrystal_threshold=COCRYSTAL_THRESHOLD):
    """Predict salt vs cocrystal from the ΔpKa rule [Cruz-Cabeza 2012]. Returns a dict:
    `delta_pka`, `zone` ('salt' | 'cocrystal' | 'continuum'), `confident` (bool), the
    thresholds, and a `caption`. This is TRIAGE, not proof: the continuum
    (cocrystal_threshold ≤ ΔpKa ≤ salt_threshold, default 0–3) is unpredictable — confirm with
    FTIR (ionised-group band) or ssNMR. Never raises."""
    d = delta_pka(pka_acid, pka_base_conjugate)
    if d > salt_threshold:
        zone = "salt"
    elif d < cocrystal_threshold:
        zone = "cocrystal"
    else:
        zone = "continuum"
    confident = zone != "continuum"
    tail = "" if confident else " — UNPREDICTABLE, confirm experimentally"
    caption = (f"ΔpKa = {d:.2f} [= pKa(BH+) {float(pka_base_conjugate):.2f} − "
               f"pKa(AH) {float(pka_acid):.2f}] → {zone}{tail}; rule: >{salt_threshold:g} salt, "
               f"<{cocrystal_threshold:g} cocrystal, else continuum [Cruz-Cabeza 2012 pKa rule]")
    return {"delta_pka": d, "zone": zone, "confident": confident,
            "salt_threshold": float(salt_threshold),
            "cocrystal_threshold": float(cocrystal_threshold), "caption": caption}


# ================================== NNLS sum-of-parents (new phase vs physical mixture)
# The #1 cocrystal-ID question: did a NEW phase form, or is the product just a physical
# mixture of the two starting materials? A physical mixture is a NON-NEGATIVE linear
# superposition of the parent patterns (you can't have a negative amount of a phase); a new
# phase (cocrystal / salt) has peaks the parents cannot build. So fit the product as
# sum-of-parents by NNLS and read the residual: low Rwp + no unexplained peaks ⇒ mixture;
# high Rwp + positive residual peaks ⇒ a new phase is plausible. Works on PXRD or FTIR
# profiles on a common grid. This is a DISCRIMINATION diagnostic, not quantitative phase
# analysis (no reference-intensity-ratio / Rietveld) — the verdict is heuristic; confirm by
# indexing / FTIR (cocrystal_id.md) / SCXRD.


def _resolve_weight(yo, weight):
    """'auto' → 'poisson' only when the profile carries a real background, min(y) > 1e-3·max(y)
    (counts-like data); a background-subtracted, floored or normalised profile gets 'unit'. A
    RELATIVE floor, not an exact zero, so the choice cannot flip on one channel and a profile
    floored at 1e-6 is not Poisson-weighted into a meaningless Rwp."""
    if weight == "auto":
        return "poisson" if (yo.size and float(np.min(yo)) > 1e-3 * float(np.max(yo))) else "unit"
    return weight


def rwp(y_obs, y_calc, weight="auto", eps=1e-12):
    """Weighted profile residual (Rietveld goodness-of-fit):
        Rwp = sqrt( Σ w_i (y_obs_i − y_calc_i)² / Σ w_i y_obs_i² )
    Returns the FRACTION (0 = perfect; ×100 for %). weight:
      'auto'    → (default) 'poisson' when every y_obs > 0, else 'unit' (see _resolve_weight).
      'poisson' → w_i = 1/max(y_obs_i, eps) — counting statistics, the crystallographic default
                  (assumes counts-like data with a non-zero background; on a normalised pattern
                  with true zeros prefer 'unit').
      'unit'    → w_i = 1 — an Rp-style unweighted profile factor, apt for normalised patterns.
    The NNLS fit itself (sum_of_parents) is UNWEIGHTED least squares; Rwp is a reported
    diagnostic on top of it."""
    yo = np.asarray(y_obs, float).ravel()
    yc = np.asarray(y_calc, float).ravel()
    weight = _resolve_weight(yo, weight)
    if weight == "unit":
        w = np.ones_like(yo)
    elif weight == "poisson":
        w = 1.0 / np.maximum(yo, eps)
    else:
        raise ValueError("weight must be 'auto', 'poisson' or 'unit'")
    den = float(np.sum(w * yo ** 2))
    if den <= 0:
        return float("nan")
    return float(np.sqrt(np.sum(w * (yo - yc) ** 2) / den))


def sum_of_parents(y_obs, parents, weight="auto", eps=1e-12):
    """Fit an observed pattern as a NON-NEGATIVE linear combination of the parent patterns
    (NNLS) — the physical-mixture model. A true mixture reconstructs well (low Rwp, no
    systematic unexplained peaks); a genuine new phase does NOT (the parents can't build its
    new peaks). Patterns must share a common x-grid.

    y_obs   : (n,) observed intensities.
    parents : (k, n) array or list of k (n,)-patterns on the SAME grid (a (n, k) array is
              accepted too).
    Returns dict: `coefficients` (k,), `fractions` (coefficients / Σ — RELATIVE scale, NOT
    quantitative phase % without an RIR), `y_calc` (n,), `residual` (n,), `rwp`, `resid_norm`,
    `weight`, `caption`. NNLS via scipy.optimize.nnls (lazy import)."""
    from scipy.optimize import nnls
    y = np.asarray(y_obs, float).ravel()
    P = np.asarray(parents, float)
    if P.ndim == 1:
        P = P[None, :]
    if P.shape[1] != y.size and P.shape[0] == y.size:
        P = P.T                                  # accept (n, k) too
    if P.shape[1] != y.size:
        raise ValueError(f"parents shape {np.asarray(parents).shape} incompatible with y (n={y.size})")
    A = P.T                                       # (n, k): columns are parents
    coef, rnorm = nnls(A, y)
    y_calc = A @ coef
    resid = y - y_calc
    total = float(coef.sum())
    frac = coef / total if total > 0 else np.full_like(coef, np.nan)
    weight = _resolve_weight(np.asarray(y, float).ravel(), weight)   # record the weighting actually used
    rw = rwp(y, y_calc, weight=weight, eps=eps)
    caption = (f"NNLS sum-of-parents: {coef.size} parents, relative scale "
               f"{np.array2string(frac, precision=3)} (NOT quantitative phase % — no RIR); "
               f"Rwp={rw*100:.1f}% [{weight} weights]. High Rwp / positive residual peaks ⇒ the "
               f"parents can't reconstruct the product ⇒ a new phase is plausible (confirm).")
    return {"coefficients": coef, "fractions": frac, "y_calc": y_calc, "residual": resid,
            "rwp": rw, "resid_norm": float(rnorm), "weight": weight, "caption": caption}


def unexplained_peaks(x, residual, reference=None, kind="new", rel_height=0.05, distance=None):
    """Peaks in the sum-of-parents residual. kind='new' → POSITIVE residual (intensity the
    parents can't explain — new-phase evidence); kind='lost' → NEGATIVE residual (parent
    intensity the product lacks). The height threshold is rel_height × max(|reference|) (or
    max(|residual|) if reference is None). Returns [(x, height), …] sorted by descending
    height. scipy.signal.find_peaks (lazy import)."""
    from scipy.signal import find_peaks
    x = np.asarray(x, float).ravel()
    r = np.asarray(residual, float).ravel()
    signal = r if kind == "new" else -r
    ref = float(np.max(np.abs(reference))) if reference is not None else float(np.max(np.abs(r)))
    thr = rel_height * ref if ref > 0 else 0.0
    idx, _ = find_peaks(signal, height=thr, distance=distance)
    peaks = [(float(x[i]), float(signal[i])) for i in idx]
    peaks.sort(key=lambda t: -t[1])
    return peaks


def phase_report(x, y_obs, parents, weight="auto", peak_rel_height=0.05):
    """Turn-key new-phase-vs-physical-mixture report: NNLS sum-of-parents fit + Rwp + the
    new/lost unexplained-peak lists + a HEURISTIC verdict. The verdict is a detection (are
    there unexplained peaks above `peak_rel_height`?), not a magic Rwp cutoff — read it WITH
    the Rwp and confirm structurally (indexing / FTIR ionised-group band / SCXRD). Returns the
    `sum_of_parents` dict plus `new_peaks`, `lost_peaks`, `verdict`."""
    fit = sum_of_parents(y_obs, parents, weight=weight)
    new = unexplained_peaks(x, fit["residual"], y_obs, kind="new", rel_height=peak_rel_height)
    lost = unexplained_peaks(x, fit["residual"], y_obs, kind="lost", rel_height=peak_rel_height)
    verdict = ("unexplained intensity present — inconsistent with a pure physical mixture; a "
               "new phase (cocrystal/salt) is plausible" if new else
               "no unexplained peaks above threshold — consistent with a physical mixture")
    caption = (f"{fit['caption']} {len(new)} new / {len(lost)} lost peak(s) at "
               f"rel_height {peak_rel_height:.0%}. Verdict (heuristic — confirm by "
               f"indexing/FTIR/SCXRD): {verdict}.")
    return {**fit, "new_peaks": new, "lost_peaks": lost, "verdict": verdict, "caption": caption}


# =============================================== COCRYSTAL-FORMATION PREDICTOR (a-priori screen)
# Transparent, OPEN reimplementation of the two established pre-screens — to be VALIDATED, never
# trusted blindly (see skill_validation/cocrystal_predictor/ + cocrystal_id.md):
#   • molecular complementarity  — shape + polarity similarity  [Fabián 2009, CGD 9:1436]
#   • H-bond synthon competition — best-donor→best-acceptor, hetero vs homo  [Etter 1990, Acc.Chem.Res. 23:120]
# combined with delta_pKa (Cruz-Cabeza 2012, above). This does NOT reproduce the CSD-TRAINED H-bond
# propensity (Galek 2007) — those coefficients come from CSD statistics we don't have; cross-check
# our complementarity against Mercury's Molecular Complementarity tool instead (the reference impl).
# Every function returns numbers + a provenance caption and a *per-signal* breakdown — the value is a
# legible prediction you can audit, not a black-box score. RDKit is a lazy import (optional dep).

_RDKIT_SEED = 0xC0CC  # fixed → the conformer ensemble is deterministic/reproducible


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
        return Chem, AllChem, rdMolDescriptors
    except ImportError as e:
        raise ImportError("the cocrystal predictor needs RDKit: python -m pip install rdkit") from e


def shape_descriptors(coords):
    """Two size-normalised shape descriptors (S/L, M/L) from the principal-axis bounding box of a
    point cloud (N×3): eigen-decompose the coordinate covariance, take the box EXTENT along each
    principal axis, sort S≤M≤L, return (S/L, M/L). PURE. Translation/rotation-invariant by
    construction and scale-invariant (ratios) — both are asserted in the validation harness."""
    X = np.asarray(coords, float)
    C = X - X.mean(0)
    _, V = np.linalg.eigh(C.T @ C)
    ext = np.ptp(C @ V, axis=0)
    S, M, L = np.sort(ext)
    return (float(S / L), float(M / L)) if L > 1e-9 else (1.0, 1.0)


def gasteiger_dipole_debye(mol, coords):
    """|dipole| in Debye from RDKit Gasteiger partial charges and a conformer's coordinates
    (Å): μ = |Σ qᵢ rᵢ| × 4.803. Approximate (Gasteiger, not QM) — a RELATIVE polarity descriptor."""
    from rdkit.Chem import AllChem
    AllChem.ComputeGasteigerCharges(mol)
    q = np.nan_to_num([a.GetDoubleProp("_GasteigerCharge") for a in mol.GetAtoms()])
    return float(np.linalg.norm((q[:, None] * np.asarray(coords, float)).sum(0)) * 4.803)


def molecule_descriptors(smiles, n_conf=12, seed=_RDKIT_SEED):
    """3-D descriptors for a SMILES via an RDKit conformer ENSEMBLE — flexible molecules (diacid
    chains!) change shape per conformer, so we report the MEDIAN and the spread (sd), which the
    validation uses as a stability check. Returns dict: s_l, m_l (+ *_sd), dipole_debye (+ sd),
    hbd, hba, formula, n_conf. HBD/HBA are graph-based (RDKit CalcNumHBD/HBA)."""
    Chem, AllChem, rd = _rdkit()
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(base)
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, randomSeed=seed))
    shapes, dips = [], []
    for cid in ids:
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid)
        except Exception:
            pass
        xyz = mol.GetConformer(cid).GetPositions()
        shapes.append(shape_descriptors(xyz))
        dips.append(gasteiger_dipole_debye(mol, xyz))
    shapes, dips = np.asarray(shapes), np.asarray(dips)
    return {"smiles": smiles, "formula": rd.CalcMolFormula(base),
            "s_l": float(np.median(shapes[:, 0])), "m_l": float(np.median(shapes[:, 1])),
            "s_l_sd": float(shapes[:, 0].std()), "m_l_sd": float(shapes[:, 1].std()),
            "dipole_debye": float(np.median(dips)), "dipole_sd": float(dips.std()),
            "hbd": int(rd.CalcNumHBD(base)), "hba": int(rd.CalcNumHBA(base)), "n_conf": len(ids)}


# Fabián-style thresholds: a pair is "complementary" (cocrystal-favoured) when its molecules are
# SIMILAR in shape and polarity. These cutoffs are APPROXIMATE and are CALIBRATED on the benchmark
# (validate.py) — not sacred constants; report the raw differences alongside the verdict.
MC_THRESHOLDS = {"s_l": 0.20, "m_l": 0.20, "dipole_debye": 4.0}


def molecular_complementarity(desc_a, desc_b, thresholds=None):
    """Fabián 2009 molecular-complementarity screen (approx): two molecules are cocrystal-favoured
    when their shape (S/L, M/L) and polarity (dipole) are SIMILAR. Returns per-descriptor absolute
    differences, which passed, an overall `complementary` bool, and a caption. TRIAGE — cross-check
    against Mercury's Molecular Complementarity tool (the CSD-calibrated reference)."""
    t = thresholds or MC_THRESHOLDS
    d = {"s_l": abs(desc_a["s_l"] - desc_b["s_l"]), "m_l": abs(desc_a["m_l"] - desc_b["m_l"]),
         "dipole_debye": abs(desc_a["dipole_debye"] - desc_b["dipole_debye"])}
    passes = {k: d[k] <= t[k] for k in d}
    ok = all(passes.values())
    caption = (f"molecular complementarity [Fabián 2009, open approx]: |Δ(S/L)|={d['s_l']:.2f}, "
               f"|Δ(M/L)|={d['m_l']:.2f}, |Δμ|={d['dipole_debye']:.1f} D → "
               f"{'complementary (cocrystal-favoured)' if ok else 'NOT complementary'} "
               f"(thresholds {t}; cross-check vs Mercury)")
    return {"complementary": ok, "differences": d, "passes": passes, "caption": caption}


# Functional groups (SMARTS) with ordinal donor/acceptor strengths (0 = n/a). Ranking follows the
# standard H-bond hierarchy behind Etter's rules (acid O–H > imide/amide N–H > O–H; carboxylate/
# sp2-N acceptors > carbonyl O > hydroxyl/ether O). Ordinal + heuristic — VALIDATE predicted
# synthons against known ones (validate.py), don't treat the numbers as energies.
_FG = [
    ("carboxylic acid O–H", "[CX3](=[OX1])[OX2H1]", 5, 4),
    ("imide N–H",           "[NX3H1]([CX3]=[OX1])[CX3]=[OX1]", 4, 3),
    ("amide N–H",           "[NX3;H1,H2][CX3]=[OX1]", 3, 3),
    ("aromatic N–H",        "[nH]", 4, 0),                 # imidazole/pyrrole (theophylline N7–H!)
    ("hydroxyl/enol O–H",   "[OX2H1]", 3, 2),
    ("pyridine/sp2 N",      "[nX2,$([NX2]=C)]", 0, 5),
    ("primary/sec amine",   "[NX3;H1,H2;!$(NC=O)]", 2, 4),
    ("carbonyl O",          "[CX3]=[OX1]", 0, 3),
    ("aromatic N (sub.)",   "[nX3]", 0, 2),
    ("sulfonyl O",          "[SX4](=[OX1])(=[OX1])", 0, 2),
]


def hbond_groups(smiles):
    """H-bond functional groups on a molecule → {groups, donors:[(name,strength)],
    acceptors:[(name,strength)]} via SMARTS. The input to the Etter synthon-competition call."""
    Chem, _, _ = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    groups, donors, acceptors = [], [], []
    for name, smarts, ds, as_ in _FG:
        n = len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))
        if n:
            groups.append((name, n))
            if ds:
                donors.append((name, ds))
            if as_:
                acceptors.append((name, as_))
    return {"groups": groups, "donors": donors, "acceptors": acceptors}


def hbond_competition(smiles_a, smiles_b):
    """Etter's rules as a synthon predictor: the best HETEROmeric donor→acceptor pair (A's donor to
    B's acceptor, or vice-versa) vs each partner's best HOMOmeric pair. If the heterosynthon score
    (donor+acceptor strengths) ≥ both homosynthons, a cocrystal is favoured. Returns the predicted
    `heterosynthon`, its score, the homo scores, a `hetero_favoured` bool, and a caption. Heuristic —
    validate the predicted synthon against the known one."""
    ga, gb = hbond_groups(smiles_a), hbond_groups(smiles_b)
    def pairs(donors, acceptors, tag_d, tag_a):
        return [(ds + as_, f"{tag_d}:{dn}···{tag_a}:{an}") for dn, ds in donors for an, as_ in acceptors]
    hetero = pairs(ga["donors"], gb["acceptors"], "A", "B") + pairs(gb["donors"], ga["acceptors"], "B", "A")
    homo_a = pairs(ga["donors"], ga["acceptors"], "A", "A")
    homo_b = pairs(gb["donors"], gb["acceptors"], "B", "B")
    best_h = max(hetero) if hetero else (0, "none")
    best_a = max(homo_a) if homo_a else (0, "none")
    best_b = max(homo_b) if homo_b else (0, "none")
    favoured = best_h[0] >= max(best_a[0], best_b[0]) and best_h[0] > 0
    caption = (f"H-bond synthon competition [Etter 1990]: best heterosynthon {best_h[1]} "
               f"(score {best_h[0]}) vs homo A {best_a[0]} / B {best_b[0]} → "
               f"{'heterosynthon competitive → cocrystal favoured' if favoured else 'homosynthons win → cocrystal disfavoured'}")
    return {"heterosynthon": best_h[1], "hetero_score": best_h[0],
            "homo_scores": (best_a[0], best_b[0]), "hetero_favoured": favoured, "caption": caption}


def predict_cocrystal(smiles_a, smiles_b, pka_acid=None, pka_base_conjugate=None,
                      thresholds=None, n_conf=12):
    """A-priori cocrystal SCREEN — three transparent, separately-reported signals (NOT a go/no-go
    oracle). The validation harness (skill_validation/cocrystal_predictor/) showed the combined
    binary does NOT beat a label-scramble control on a small diverse set, so this deliberately does
    NOT emit a single `likely` verdict. Use it for:
      • `synthon`      — the predicted heterosynthon [Etter 1990]: a genuinely useful FTIR-interpretation
                         aid ("expect THESE bands to shift"). Sensitive, but non-specific for yes/no.
      • `ionisation`   — the ΔpKa salt-vs-cocrystal filter [Cruz-Cabeza 2012] (if pKa given).
      • `complementarity` — shape/polarity [Fabián 2009]: report the raw Δ's to CROSS-CHECK against
                         Mercury's Molecular Complementarity; weak for rigid-API + flexible-coformer pairs.
    Go/no-go formation belongs to Mercury (CSD-calibrated) + the experiment. Pass pKa only for an acid/base pair."""
    da, db = molecule_descriptors(smiles_a, n_conf=n_conf), molecule_descriptors(smiles_b, n_conf=n_conf)
    comp = molecular_complementarity(da, db, thresholds)
    hb = hbond_competition(smiles_a, smiles_b)
    ion = classify_ionisation(pka_acid, pka_base_conjugate) if (pka_acid is not None and pka_base_conjugate is not None) else None
    caption = (f"cocrystal SCREEN(A,B): predicted synthon = {hb['heterosynthon']}"
               + (f"; ΔpKa → {ion['zone']}" if ion else "")
               + f"; complementarity |Δμ|={comp['differences']['dipole_debye']:.1f} D (cross-check Mercury). "
               f"NO binary verdict — validated as no-better-than-chance for formation; defer go/no-go to Mercury + experiment.")
    return {"descriptors": (da, db), "synthon": hb["heterosynthon"], "complementarity": comp,
            "hbond": hb, "ionisation": ion,
            "formation_confidence": "low — open screen did not beat chance in retrospective validation (skill_validation/cocrystal_predictor/validate.py, L3)",
            "caption": caption}
