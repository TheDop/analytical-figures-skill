"""VIP: the Sum(VIP^2)==p normalization invariant (Chong & Jun 2005) and guards."""
from __future__ import annotations
import numpy as np
from scripts import config, chemometrics as cm
from _harness import close, prop, skip
import _data


def collect():
    R = []
    cfg = config.Config()
    x, X, y, g = _data.synth(seed=5)
    try:
        import sklearn  # noqa: F401
    except Exception as e:
        return [skip("vip", "all", f"sklearn unavailable ({type(e).__name__})")]

    model = cm.pls_fit(X, y, 4, cfg)
    v = cm.vip(model)
    p = int(np.asarray(model["coef"]).size)
    R.append(close("vip", "sum(VIP^2) == p (mean(VIP^2)=1)", float(np.sum(v ** 2)), float(p), "T0",
                   "Chong & Jun 2005 doi:10.1016/j.chemolab.2004.12.011; normalization"))
    R.append(prop("vip", "VIP finite & non-negative", bool(np.all(np.isfinite(v)) and np.all(v >= 0)), "T3", ""))

    fake = {"weights": model["weights"], "scores": model["scores"],
            "y_loadings": np.zeros_like(model["y_loadings"])}
    R.append(prop("vip", "VIP -> zeros when explained-Y == 0", bool(np.allclose(cm.vip(fake), 0.0)), "T3",
                  "guard: total explained Y SS == 0"))
    return R
