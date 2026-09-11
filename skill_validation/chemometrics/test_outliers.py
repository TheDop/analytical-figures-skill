"""Outlier / diagnostic gate: PCA-based T2 / Q / leverage + DD-SIMCA.

Primary anchor = the dependency-free hand-derived closed-form fixture (T0), each limit also
recomputed from scipy in-test (independent T0). Published cross-check (T2) = the R `mdatools`
people Hotelling T2 limits, which depend only on (n, A) so they need no dataset. Plus the
documented edge cases (SPEC.md §8 rows 11-14).
"""
from __future__ import annotations
import numpy as np
from scripts import config, chemometrics as cm
from _harness import close, prop, skip
from golden import closed_form as cf
from golden import manifest as gm


def collect():
    R = []
    om = cm.pca_outlier_model(cf.X, None, n_components=cf.A, pre="none", scale=False)
    d = cm.outlier_distances(om)

    # --- distances vs the hand-derived fixture. The literals in closed_form.py are printed to
    # ~7 sig figs, so compare at rtol=1e-5 (proves 5+ figure agreement with the independent
    # derivation); machine-precision T0 proof is the scipy recompute + invariants below.
    LIT = dict(rtol=1e-5)
    R.append(close("outliers", "eigenvalues vs closed-form", om["eig"][:3], cf.EIGENVALUES_DDOF1, "T0",
                   "hand-derived; SVD eig ddof=1", **LIT))
    R.append(close("outliers", "T2 per sample vs closed-form", d["t2"], cf.T2_PER_SAMPLE, "T0",
                   "Hotelling T2 = sum t_a^2/eig_a", **LIT))
    R.append(close("outliers", "Q per sample vs closed-form", d["q"], cf.Q_PER_SAMPLE, "T0",
                   "Q = ||Xc||^2 - sum t_A^2", **LIT))

    # --- invariants (T0/T3)
    R.append(prop("outliers", "T2 == (n-1)*leverage", bool(np.allclose(d["t2"], (om["n"] - 1) * d["leverage"])),
                  "T0", "identity"))
    R.append(close("outliers", "mean(leverage) == A/n", d["leverage"].mean(), cf.A / om["n"], "T0",
                   "DD-SIMCA SD normalization"))
    R.append(close("outliers", "sum(leverage) == A", d["leverage"].sum(), float(cf.A), "T0",
                   "hat trace (no mean term); matches chemometrics lib"))

    # --- T2 limits vs closed-form (T0) and an independent scipy recompute (T0)
    from scipy.stats import f as fdist, chi2
    n, A = om["n"], om["A"]
    forms = {"F_textbook": cf.T2_LIMITS["F_textbook"], "F_n2m1": cf.T2_LIMITS["F_n2m1"],
             "chi2": cf.T2_LIMITS["chi2"], "beta": cf.T2_LIMITS["beta_tracy_young"]}
    for form, want in forms.items():
        R.append(close("outliers", f"T2 limit {form} vs closed-form", cm.t2_limit(om, 0.05, form), want, "T0",
                       "Hotelling / Tracy-Young-Mason 1992", **LIT))
    R.append(close("outliers", "T2 limit F_textbook == scipy recompute", cm.t2_limit(om, 0.05, "F_textbook"),
                   A * (n - 1) / (n - A) * float(fdist.ppf(0.95, A, n - A)), "T0", "scipy.stats.f.ppf"))

    # --- Q limits vs closed-form (T0)
    R.append(close("outliers", "Q limit Jackson-Mudholkar vs closed-form",
                   cm.q_limit(om, 0.05, "jackson_mudholkar"), cf.Q_LIMITS["jackson_mudholkar"], "T0",
                   "Jackson & Mudholkar 1979 doi:10.1080/00401706.1979.10489779", **LIT))
    R.append(close("outliers", "Q limit Box/NM vs closed-form",
                   cm.q_limit(om, 0.05, "box", q_values=d["q"]), cf.Q_LIMITS["box_nomikos_macgregor"], "T0",
                   "Nomikos & MacGregor 1995 doi:10.1080/00401706.1995.10485888", **LIT))

    # --- DD-SIMCA limits vs locked fixture (T0) and scipy recompute (T0)
    dd = cm.ddsimca_limits(om, d["sd"], d["od"], 0.05, 0.01)
    R.append(prop("outliers", "DD-SIMCA Nh,Nq vs locked",
                  dd["Nh"] == cf.DDSIMCA["Nh"] and dd["Nq"] == cf.DDSIMCA["Nq"], "T0", "method-of-moments DoF"))
    R.append(close("outliers", "DD-SIMCA h0,q0 vs locked", [dd["h0"], dd["q0"]],
                   [cf.DDSIMCA["h0"], cf.DDSIMCA["q0"]], "T0", "Pomerantsev & Rodionova 2014"))
    R.append(close("outliers", "DD-SIMCA c_crit,c_out vs locked", [dd["c_crit"], dd["c_out"]],
                   [cf.DDSIMCA["Dcrit"], cf.DDSIMCA["c_out"]], "T0", "doi:10.1002/cem.2506"))
    R.append(close("outliers", "DD-SIMCA c_crit == scipy chi2", dd["c_crit"],
                   float(chi2.ppf(0.95, dd["Nh"] + dd["Nq"])), "T0", "scipy.stats.chi2.ppf"))

    # --- published cross-check: R mdatools `people` Hotelling T2 limits (depend only on n,A) [T2]
    for i, want in enumerate(gm.MDATOOLS_PEOPLE_T2LIM, start=1):
        R.append(close("outliers", f"T2 limit vs mdatools people (A={i}, n=32)",
                       cm.t2_limit({"A": i, "n": 32}, 0.05, "F_textbook"), want, "T2",
                       gm.MDATOOLS_SOURCE, rtol=1e-4,
                       note="Hotelling limit is a function of (n,A) only -> no dataset needed"))

    # --- edge cases (SPEC.md §8)
    om3 = cm.pca_outlier_model(cf.X, None, n_components=3, pre="none", scale=False)  # A==rank
    R.append(prop("outliers", "A>=rank -> JM Q limit = 0 (uninformative)",
                  cm.q_limit(om3, 0.05, "jackson_mudholkar") == 0.0, "T3", "no residual space"))
    R.append(close("outliers", "n-A<=0 -> T2 chi2 fallback", cm.t2_limit({"A": 5, "n": 5}, 0.05, "F_textbook"),
                   float(chi2.ppf(0.95, 5)), "T3", "F undefined -> chi2 guard"))

    # --- DoF estimators recover the true k from chi2(k)/k samples (self-consistent; T3)
    rng = np.random.default_rng(0)
    for k in (3, 8):
        u = rng.chisquare(k, size=6000) / k
        ncl = cm._dd_dof(u, "classical")[0]
        nrob = cm._dd_dof(u, "robust")[0]
        R.append(prop("outliers", f"classical DoF recovers k={k}", abs(ncl - k) <= 1, "T3", f"got {ncl}"))
        R.append(prop("outliers", f"robust DoF recovers k={k}", abs(nrob - k) <= 1, "T3", f"got {nrob}"))

    # --- end-to-end diagnostics structural check (classical DoF to match the locked dd above)
    diag = cm.diagnostics(cf.X, n_components=2, alpha=0.05, gamma=0.01, dof="classical")
    f_manual = dd["Nh"] * (d["sd"] / dd["h0"]) + dd["Nq"] * (d["od"] / dd["q0"])
    R.append(close("outliers", "diagnostics f == Nh*SD/h0 + Nq*OD/q0", diag["f"], f_manual, "T0", "combined distance"))
    R.append(prop("outliers", "flags valid & length n",
                  bool(len(diag["flags"]) == om["n"] and set(np.unique(diag["flags"])) <= {"regular", "extreme", "outlier"}),
                  "T3", "structural"))

    # --- SD-OD acceptance plot renders on Agg without error
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        cm.plot_influence(ax, diag, labels=[f"s{i}" for i in range(om["n"])])
        ok = len(ax.collections) >= 1
        plt.close(fig)
        fig2, ax2 = plt.subplots()
        cm.plot_influence(ax2, diag, axis_scale="sqrt")        # dynamic-range mode
        ok = ok and len(ax2.collections) >= 1
        plt.close(fig2)
    except Exception:
        ok = False
    R.append(prop("outliers", "plot_influence renders (linear + sqrt)", ok, "T3", "DD-SIMCA SD-OD acceptance plot"))

    return R
