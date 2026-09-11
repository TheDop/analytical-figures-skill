"""PLS: engine identity vs a direct sklearn fit (T1), predict reconstruction (T1),
and the R `pls` gasoline RMSEP regression target (T2, skips if dataset absent).
"""
from __future__ import annotations
import numpy as np
from pathlib import Path
from scripts import config, chemometrics as cm
from _harness import close, prop, skip
import _data


def collect():
    R = []
    cfg = config.Config()                     # pls_preprocess='snv', pls_scale=True
    x, X, y, g = _data.synth(seed=2)
    try:
        from sklearn.cross_decomposition import PLSRegression
    except Exception as e:
        return [skip("pls", "all", f"sklearn unavailable ({type(e).__name__})")]

    model = cm.pls_fit(X, y, 3, cfg)
    nc = model["n_components"]
    Xp = cm.Preprocessor(cfg.pls_preprocess, cfg).fit_transform(X)
    m = PLSRegression(n_components=nc, scale=cfg.pls_scale).fit(Xp, y)

    coef_ref = np.ravel(m.coef_)
    if coef_ref.size != Xp.shape[1]:
        coef_ref = np.ravel(np.asarray(m.coef_).T)
    R.append(close("pls", "coef == direct PLSRegression", model["coef"], coef_ref, "T1",
                   "sklearn PLSRegression NIPALS (engine identity)"))
    R.append(close("pls", "yhat == direct PLSRegression", model["yhat"], np.ravel(m.predict(Xp)), "T1",
                   "sklearn PLSRegression"))

    # predict reconstruction identity (documents the formula used elsewhere)
    try:
        xmean = getattr(m, "x_mean_", getattr(m, "_x_mean", None))
        intercept = np.ravel(m.intercept_)
        recon = (Xp - xmean) @ np.atleast_2d(coef_ref).T
        recon = np.ravel(recon) + intercept
        R.append(close("pls", "yhat == (X-x_mean)@coef.T + intercept", np.ravel(m.predict(Xp)), recon, "T1",
                       "sklearn predict reconstruction; intercept_==y_mean"))
    except Exception as e:
        R.append(skip("pls", "predict reconstruction", f"attr layout differs ({type(e).__name__})"))

    R.append(prop("pls", "R2(cal) in (0,1] on separable synth", 0.0 < model["r2_cal"] <= 1.0, "T3", "sanity"))

    R.append(_gasoline())
    return R


def _gasoline():
    p = Path(__file__).resolve().parent / "datasets" / "gasoline.npz"
    if not p.exists():
        return skip("pls", "gasoline LOO RMSEP",
                    "datasets/gasoline.npz absent (run datasets/fetch.py)",
                    "R pls::gasoline, Kalivas 1997")
    import numpy as np
    d = np.load(p)
    X, y = d["X"], d["y"]
    cfg = config.Config(pls_scale=False)           # R pls default: scale=FALSE
    cv = cm.cross_validate(X[:50], y[:50], 2, cfg, pre="center", report=False)
    # CV-row RMSEP at 2 comps (LOO) printed in the pls vignette
    return close("pls", "gasoline LOO RMSEP @2c", cv["rmsecv"], 0.2966, "T2",
                 "R pls vignette summary(gas1) CV row, ncomp=2", rtol=5e-2,
                 note="cross-language (R kernel-PLS vs sklearn NIPALS); tolerance only")
