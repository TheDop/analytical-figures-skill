"""NNLS sum-of-parents: closed-form Rwp + fraction recovery on a synthetic physical mixture +
new-phase discrimination + scipy.optimize.nnls parity + new/lost peak detection.

The synthetic sum_of_parents checks use `weight='unit'` (the patterns are normalised, with true
zeros between peaks where poisson weighting is ill-defined); the poisson Rwp FORMULA is validated
separately on zero-free hand cases.
"""
from __future__ import annotations
import numpy as np
from scripts import cocrystal as cc
from _harness import close, prop

A = "sum_of_parents"


def _gauss(x, c, w, h):
    return h * np.exp(-0.5 * ((x - c) / w) ** 2)


def _pattern(x, peaks, w=0.5):
    y = np.zeros_like(x)
    for c, h in peaks:
        y = y + _gauss(x, c, w, h)
    return y


def _grid():
    x = np.linspace(5.0, 40.0, 701)
    pa = _pattern(x, [(10, 1.0), (20, 0.6), (30, 0.4)])
    pb = _pattern(x, [(15, 0.8), (25, 1.0), (35, 0.5)])
    return x, pa / pa.sum(), pb / pb.sum()          # unit total → coefficients read as fractions


def collect():
    R = []

    # ---------- Rwp closed-form (T0) ----------
    # poisson w=1/yo → Rwp = sqrt( Σ(yo−yc)²/yo / Σyo );  [100,100] vs [110,90] → sqrt(2/200)=0.1
    R.append(close("rwp", "Rwp poisson hand == 0.1",
                   cc.rwp([100., 100.], [110., 90.], "poisson"), 0.1, "T0", "Σ(Δ²/yo)/Σyo"))
    R.append(close("rwp", "Rwp poisson asym == sqrt(1/125)",
                   cc.rwp([100., 25.], [100., 20.], "poisson"), float(np.sqrt(1 / 125)), "T0", "hand"))
    R.append(close("rwp", "Rwp unit asym == sqrt(25/10625)",
                   cc.rwp([100., 25.], [100., 20.], "unit"), float(np.sqrt(25 / 10625)), "T0", "hand"))
    R.append(prop("rwp", "Rwp == 0 for a perfect fit",
                  cc.rwp([1., 2., 3.], [1., 2., 3.], "unit") == 0.0, "T0", "identity"))

    # ---------- physical mixture: recover the known mixing fractions ----------
    x, pa, pb = _grid()
    y_mix = 0.7 * pa + 0.3 * pb
    fit = cc.sum_of_parents(y_mix, [pa, pb], weight="unit")
    R.append(close(A, "NNLS recovers mix fractions [0.7, 0.3]", fit["fractions"], [0.7, 0.3], "T1", "ground truth"))
    R.append(prop(A, "mixture Rwp ~ 0 (parents fully explain it)", fit["rwp"] < 1e-6, "T1", "consistent system"))
    R.append(prop(A, "NNLS coefficients are non-negative", bool(np.all(fit["coefficients"] >= 0)), "T3", "NNLS"))

    # parity: our coefficients == scipy.optimize.nnls on the same design (fit is weight-independent)
    from scipy.optimize import nnls
    ref_coef, _ = nnls(np.column_stack([pa, pb]), y_mix)
    R.append(close(A, "coefficients == scipy.optimize.nnls", fit["coefficients"], ref_coef, "T1", "scipy.optimize.nnls"))

    # ---------- new phase: parents CANNOT reconstruct new peaks ----------
    y_new = _pattern(x, [(12, 1.0), (22, 0.9), (32, 0.5)])
    y_new = y_new / y_new.sum()
    fit_new = cc.sum_of_parents(y_new, [pa, pb], weight="unit")
    R.append(prop(A, "new-phase Rwp >> mixture Rwp", fit_new["rwp"] > 100 * max(fit["rwp"], 1e-12), "T3", "discrimination"))
    R.append(prop(A, "new-phase Rwp is substantial (>10%)", fit_new["rwp"] > 0.10, "T3", "discrimination"))

    # ---------- new/lost peak detection ----------
    rep = cc.phase_report(x, y_new, [pa, pb], weight="unit")
    new_pos = [p[0] for p in rep["new_peaks"]]
    R.append(prop(A, "unexplained 'new' peak detected near x=12",
                  any(abs(p - 12.0) < 1.0 for p in new_pos), "T3", "find_peaks on +residual"))
    R.append(prop(A, "new-phase verdict flags a plausible new phase", "new phase" in rep["verdict"], "T3", "verdict"))
    rep_mix = cc.phase_report(x, y_mix, [pa, pb], weight="unit")
    R.append(prop(A, "pure mixture flags NO new peaks", len(rep_mix["new_peaks"]) == 0, "T3", "specificity"))

    return R
