"""Validation trio: permutation/y-scramble test, Nadeau-Bengio corrected paired t-test, RPD/RPIQ.
SPEC.md §6.9. The corrected-t golden is reconstructed THROUGH the function (build differences
with the exact mean/var that yield the published t), so it's a real reproduction, not a re-store.
"""
from __future__ import annotations
import numpy as np
from scripts import config, chemometrics as cm
from _harness import close, prop, skip
from golden import manifest as gm


def collect():
    R = []

    # ---------- RPD / RPIQ ----------
    y_sd10 = np.array([0.0, 10.0 * np.sqrt(2)])              # sample SD (ddof=1) == 10
    R.append(close("trio", "RPD unit example (SD=10, RMSEP=4 -> 2.50)", cm.rpd(y_sd10, 4.0),
                   gm.RPD_UNIT["rpd"], "T0", gm.RPD_SOURCE))
    rng = np.random.default_rng(1)
    y2 = rng.normal(50, 8, 400)
    R.append(close("trio", "RPD == SD(ddof1)/RMSE", cm.rpd(y2, 3.0), np.std(y2, ddof=1) / 3.0, "T0", "definition"))
    R.append(close("trio", "RPIQ == IQR/RMSE (R type-7)", cm.rpiq(y2, 3.0),
                   (np.percentile(y2, 75) - np.percentile(y2, 25)) / 3.0, "T0", "definition"))
    R.append(prop("trio", "RPIQ ~= 1.349*RPD for normal y",
                  abs(cm.rpiq(y2, 3.0) / cm.rpd(y2, 3.0) - 1.349) < 0.1, "T3", "Bellon-Maurel 2010"))
    R.append(prop("trio", "RPD -> inf when RMSE=0", np.isinf(cm.rpd(y2, 0.0)), "T3", "edge"))

    # ---------- Nadeau-Bengio corrected paired t ----------
    gold = gm.NADEAU_BENGIO_GOLDEN          # t=1.364576, p=0.205528, df=9, correction=0.35
    k, var = 10, 1.0
    n_train, n_test = 4, 1                   # 1/k + n_test/n_train = 0.1 + 0.25 = 0.35
    md = gold["t"] * np.sqrt(gold["correction"] * var)
    z = np.arange(k, dtype=float)
    d = md + (z - z.mean()) / z.std(ddof=1) * np.sqrt(var)   # mean=md, var=1.0 exactly
    res = cm.corrected_paired_t(d, np.zeros(k), n_train=n_train, n_test=n_test)
    R.append(close("trio", "corrected-t correction factor == 0.35", res["correction"], gold["correction"], "T0",
                   "1/k + n_test/n_train"))
    R.append(prop("trio", "corrected-t df == 9", res["df"] == gold["df"], "T0", ""))
    R.append(close("trio", "corrected-t t == golden 1.364576", res["t"], gold["t"], "T0",
                   gm.NADEAU_BENGIO_SOURCE, rtol=1e-5))
    R.append(close("trio", "corrected-t p == golden 0.205528", res["p_value"], gold["p"], "T2",
                   gm.NADEAU_BENGIO_SOURCE, rtol=1e-4))
    # edge: zero variance of differences, non-zero mean -> t = inf
    res0 = cm.corrected_paired_t(np.full(8, 0.5), np.zeros(8), n_train=4, n_test=1)
    R.append(prop("trio", "corrected-t inf when var(d)=0, mean!=0", np.isinf(res0["t"]), "T3", "edge"))

    # ---------- permutation / y-scramble ----------
    cfg = config.Config()
    x, X, y, g = _data_synth()
    perm = cm.permutation_test(X, y, 2, cfg, n_perm=99, groups=g, seed=0)
    expected_p = (int(np.sum(perm["null"] >= perm["observed"])) + 1) / (perm["n_perm"] + 1)
    R.append(close("trio", "permutation p == (#null>=obs + 1)/(n+1)", perm["p_value"], expected_p, "T0",
                   "sklearn permutation_test_score p-formula"))
    R.append(prop("trio", "permutation: observed beats 95th-pct null (real data)",
                  perm["observed"] > np.percentile(perm["null"], 95), "T3", "real relationship"))
    R.append(prop("trio", "permutation: p small for real relationship", perm["p_value"] < 0.05, "T3", ""))
    # between-level shuffle integrity
    rng2 = np.random.default_rng(0)
    yp = cm._permute_y(y, g, rng2)
    const_within = all(np.ptp(yp[g == gg]) == 0 for gg in np.unique(g))
    level_means = sorted(float(y[g == gg].mean()) for gg in np.unique(g))
    R.append(prop("trio", "between-level shuffle keeps reps' y constant per level", const_within, "T3", ""))
    R.append(prop("trio", "between-level shuffle permutes the level values",
                  np.allclose(sorted(np.unique(yp)), level_means), "T3", ""))
    return R


def _data_synth():
    import _data
    return _data.synth(seed=7)
