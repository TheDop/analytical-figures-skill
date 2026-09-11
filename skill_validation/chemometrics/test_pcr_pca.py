"""PCR / PCA: explained-variance vs sklearn (T1) and the PCR score->X back-map (T1)."""
from __future__ import annotations
import numpy as np
from scripts import config, chemometrics as cm
from _harness import close, skip
import _data


def collect():
    R = []
    cfg = config.Config()
    x, X, y, g = _data.synth(seed=3)
    try:
        from sklearn.decomposition import PCA
        from sklearn.linear_model import LinearRegression
    except Exception as e:
        return [skip("pcr_pca", "all", f"sklearn unavailable ({type(e).__name__})")]

    Xp = cm.Preprocessor(cfg.pls_preprocess, cfg).fit_transform(X)

    res = cm.pca_fit(X, cfg, n_components=4)
    p = PCA(n_components=4).fit(Xp)
    R.append(close("pcr_pca", "PCA explained_variance_ratio vs sklearn",
                   res["explained"], p.explained_variance_ratio_, "T1", "sklearn PCA"))

    res2 = cm.pcr_fit(X, y, 3, cfg)
    pca = PCA(n_components=3).fit(Xp)
    T = pca.transform(Xp)
    reg = LinearRegression().fit(T, y)
    coef = pca.components_.T @ np.ravel(reg.coef_)
    R.append(close("pcr_pca", "PCR coef back-map vs manual PCA+OLS", res2["coef"], coef, "T1",
                   "PCA-then-OLS, components_.T @ reg.coef_"))
    R.append(close("pcr_pca", "PCR yhat vs manual PCA+OLS", res2["yhat"], np.ravel(reg.predict(T)), "T1",
                   "PCA-then-OLS"))
    return R
