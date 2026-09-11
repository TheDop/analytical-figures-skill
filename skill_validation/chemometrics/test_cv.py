"""Cross-validation: RMSECV/bias definitions (T0), the SEP/RMSEP identity (T0),
the leakage-safe == hand-rolled loop guarantee (T1), and scheme/clamp properties (T3).
"""
from __future__ import annotations
import numpy as np
from scripts import config, chemometrics as cm
from _harness import close, prop, skip
import _data


def collect():
    R = []
    cfg = config.Config()
    x, X, y, g = _data.synth(seed=4)
    try:
        import sklearn  # noqa: F401
    except Exception as e:
        return [skip("cv", "all", f"sklearn unavailable ({type(e).__name__})")]

    cv = cm.cross_validate(X, y, 3, cfg, groups=g, report=False)   # leave-one-level-out
    pred = cv["pred"]
    resid = y - pred
    R.append(close("cv", "RMSECV == sqrt(mean(resid^2))", cv["rmsecv"],
                   float(np.sqrt(np.nanmean(resid ** 2))), "T0", "ddof=0; IUPAC Gold Book 10110"))
    R.append(close("cv", "bias == mean(pred-y)", cv["bias"], float(np.nanmean(pred - y)), "T0", ""))
    n = y.size
    sep = float(np.sqrt(np.nansum((resid - np.nanmean(resid)) ** 2) / (n - 1)))
    bias = cv["bias"]
    R.append(close("cv", "RMSEP^2 == SEP^2*(n-1)/n + bias^2", cv["rmsecv"] ** 2,
                   sep ** 2 * (n - 1) / n + bias ** 2, "T0", "SEP(ddof1)/RMSEP identity; ASTM E1655"))

    R.extend(_leakage(cfg, X, y, g))

    cv_loo = cm.cross_validate(X, y, 2, cfg, report=False)
    R.append(prop("cv", "LOO question names a SEEN level", "SEEN" in cv_loo["question"].upper(),
                  "T3", "optimistic-CV trap text"))
    cv_grp = cm.cross_validate(X, y, 2, cfg, groups=g, report=False)
    R.append(prop("cv", "groups -> leave-one-group-out scheme", "group" in cv_grp["scheme"].lower(), "T3", ""))

    cfg2 = config.Config(cv_folds=999)
    cvk = cm.cross_validate(X, y, 2, cfg2, scheme="kfold", report=False)
    R.append(prop("cv", "k-fold clamps k<=n", cvk["n_splits"] <= n, "T3", ""))
    return R


def _leakage(cfg, X, y, g):
    R = []
    pre = "msc+center"                       # stateful -> leaks if fit on the full set
    cv = cm.cross_validate(X, y, 3, cfg, pre=pre, groups=g, report=False)

    idx = np.arange(len(y))
    pred = np.full(len(y), np.nan)
    for grp in dict.fromkeys(g.tolist()):
        tr = idx[g != grp]
        te = idx[g == grp]
        pp = cm.Preprocessor(pre, cfg)
        Xtr = pp.fit_transform(X[tr])
        Xte = pp.transform(X[te])
        pred[te] = cm._fit_predict(Xtr, y[tr], Xte, 3, cfg, "pls")
    R.append(close("cv", "leakage-safe RMSECV == hand-rolled in-fold loop", cv["rmsecv"],
                   float(np.sqrt(np.nanmean((y - pred) ** 2))), "T1",
                   "preprocessing re-fit inside every fold"))

    # leaky: fit the stateful preprocessing on the WHOLE set, then CV the transformed matrix
    pp = cm.Preprocessor(pre, cfg)
    Xall = pp.fit_transform(X)
    leaky = cm.cross_validate(Xall, y, 3, cfg, pre="none", groups=g, report=False)
    R.append(prop("cv", "leaky preprocessing changes RMSECV (optimism)",
                  abs(leaky["rmsecv"] - cv["rmsecv"]) > 1e-6, "T3",
                  f"leaky={leaky['rmsecv']:.4g} vs safe={cv['rmsecv']:.4g}"))
    return R
