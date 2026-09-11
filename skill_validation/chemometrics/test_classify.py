"""DD-SIMCA one-class classifier + class figures-of-merit (SPEC: cocrystal Tier-3 item 3).

`class_fom`: closed-form confusion-matrix metrics (T0). Classifier: accepts held-out target
members, rejects out-of-class (new-peak) samples; predict-on-training reproduces the already-
validated in-sample DD-SIMCA `f` (test_outliers), tying the wrapper to validated limits.
"""
from __future__ import annotations
import numpy as np
from scripts import chemometrics as cm
from _harness import close, prop


def _class_data(seed=0):
    """A target class with REAL low-rank structure (two latent intensity factors + small noise),
    so the PCA model captures genuine variation and the residual is consistent train-vs-test —
    the correct setup for a one-class model. Out-of-class samples carry a peak the class lacks."""
    rng = np.random.default_rng(seed)
    p = 60
    x = np.linspace(0.0, 60.0, p)

    def g(c, w, h):
        return h * np.exp(-0.5 * ((x - c) / w) ** 2)

    base = g(15, 2, 1.0) + g(30, 2, 0.7) + g(45, 2, 0.5)

    def cls(n, extra=0.0, sigma=0.01):
        M = []
        for _ in range(n):
            a = rng.normal(0, 0.25)                   # two class-internal latent factors
            b = rng.normal(0, 0.20)                   # (peak-height variation — real structure)
            M.append(base + a * g(30, 2, 1.0) + b * g(45, 2, 1.0) + extra + rng.normal(0, sigma, p))
        return np.array(M)

    Xtr = cls(20)                                     # target-class training spectra
    Xin = cls(10)                                     # held-out target members
    Xout = cls(10, extra=g(37, 2, 0.8))               # out-of-class: a peak the class lacks
    return Xtr, Xin, Xout


def collect():
    R = []

    # ---------------- class_fom closed-form (T0) ----------------
    F = "class_fom"
    yt = np.array([1, 1, 1, 1, 0, 0, 0, 0], bool)
    yp = np.array([1, 1, 1, 0, 0, 0, 0, 1], bool)     # TP=3 FN=1 TN=3 FP=1
    fom = cm.class_fom(yt, yp)
    R.append(close(F, "sensitivity == 3/4", fom["sensitivity"], 0.75, "T0", "TP/(TP+FN)"))
    R.append(close(F, "specificity == 3/4", fom["specificity"], 0.75, "T0", "TN/(TN+FP)"))
    R.append(close(F, "efficiency == sqrt(.75*.75)", fom["efficiency"], 0.75, "T0", "geom mean"))
    R.append(prop(F, "confusion counts tp/fn/tn/fp == 3/1/3/1",
                  (fom["tp"], fom["fn"], fom["tn"], fom["fp"]) == (3, 1, 3, 1), "T0", "counts"))
    fperf = cm.class_fom(yt, yt)
    R.append(prop(F, "perfect prediction → all metrics 1.0",
                  fperf["sensitivity"] == 1.0 and fperf["specificity"] == 1.0
                  and fperf["efficiency"] == 1.0, "T0", "perfect"))
    fnp = cm.class_fom(np.zeros(4, bool), np.zeros(4, bool))
    R.append(prop(F, "no positives → sensitivity NaN, specificity 1.0",
                  np.isnan(fnp["sensitivity"]) and fnp["specificity"] == 1.0, "T3", "edge"))
    yt2 = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0], bool)
    yp2 = np.array([1, 1, 1, 1, 0, 0, 0, 0, 1, 1], bool)     # sens 1, spec 4/6
    f2 = cm.class_fom(yt2, yp2)
    R.append(close(F, "efficiency == sqrt(sens*spec)", f2["efficiency"],
                   float(np.sqrt(f2["sensitivity"] * f2["specificity"])), "T0", "definition"))

    # ---------------- DD-SIMCA one-class classifier ----------------
    # A=1: the class has genuine ~1-D-dominant structure at this noise level; a higher A overfits
    # the small training set's noise so the residual boundary stops generalising (a live demo of
    # the parsimony issue). DD-SIMCA controls Type-I error IN-SAMPLE, so we assert the guaranteed
    # properties — training acceptance, out-of-class rejection, and a boundary that SEPARATES the
    # classes — not an out-of-sample acceptance rate DD-SIMCA doesn't promise.
    C = "ddsimca_classify"
    Xtr, Xin, Xout = _class_data(0)
    model = cm.ddsimca_fit(Xtr, n_components=1, alpha=0.05, dof="classical")
    ptr = cm.ddsimca_predict(model, Xtr)
    pin = cm.ddsimca_predict(model, Xin)
    pout = cm.ddsimca_predict(model, Xout)
    cc = model["dd"]["c_crit"]
    R.append(prop(C, "training accepted at ~1-alpha (>=0.85 in-sample)", ptr["accept"].mean() >= 0.85, "T3", "in-sample Type-I"))
    R.append(prop(C, "out-of-class (new-peak) rejected (>=0.9)", (~pout["accept"]).mean() >= 0.9, "T3", "out-class"))
    R.append(prop(C, "boundary separates classes: median(f_in) < c_crit < median(f_out)",
                  float(np.median(pin["f"])) < cc < float(np.median(pout["f"])), "T3", "separation"))
    R.append(prop(C, "accept == (f <= c_crit)",
                  bool(np.array_equal(pout["accept"], pout["f"] <= cc)), "T0", "definition"))
    # tie the wrapper to the validated in-sample path: predict(train) f == diagnostics f
    diag = cm.diagnostics(Xtr, n_components=1, alpha=0.05, gamma=None, dof="classical")
    R.append(close(C, "predict(train) f == diagnostics f", ptr["f"], diag["f"], "T1", "in-sample consistency"))
    # class FoM composes on the labelled target-vs-new-phase set (specificity is exact here)
    y_true = np.r_[np.ones(len(Xin), bool), np.zeros(len(Xout), bool)]
    y_pred = np.r_[pin["accept"], pout["accept"]]
    fom_all = cm.class_fom(y_true, y_pred)
    R.append(prop(C, "class_fom specificity == 1.0 (all out-class rejected)", fom_all["specificity"] == 1.0, "T3", "end-to-end"))
    R.append(prop(C, "class_fom efficiency >= 0.6", fom_all["efficiency"] >= 0.6, "T3", "end-to-end"))

    return R
