"""Edge cases for the existing surface (SPEC.md §8 rows 1-10). The new-feature edge
cases (11-19) land with their modules in step (2).
"""
from __future__ import annotations
import warnings
import numpy as np
from scripts import config, chemometrics as cm
from _harness import prop, skip
import _data


def collect():
    R = []
    cfg = config.Config()
    x, X, y, g = _data.synth(seed=6)

    # 1. SNV flat row (sd=0) -> finite, centred to ~0
    Xf = X.copy(); Xf[0] = 5.0
    snv = cm.Preprocessor(["snv"], cfg).fit_transform(Xf)
    R.append(prop("edge", "SNV flat row -> finite & ~0",
                  bool(np.all(np.isfinite(snv[0])) and np.allclose(snv[0], 0.0)), "T3", "sd=0 guard"))

    # 2. autoscale flat column (sd=0) -> finite
    Xc = X.copy(); Xc[:, 0] = 3.0
    au = cm.Preprocessor(["autoscale"], cfg).fit_transform(Xc)
    R.append(prop("edge", "autoscale flat column -> finite", bool(np.all(np.isfinite(au))), "T3", "sd=0 guard"))

    # 3+4. SG window <= polyorder -> passthrough unchanged AND warns (the SPEC.md §5 fix)
    Xsmall = X[:, :4]                          # n_features=4, poly=3 -> usable window 3 <= 3
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        out = cm.Preprocessor(["d1"], cfg).fit_transform(Xsmall)
    R.append(prop("edge", "SG too-small window -> passthrough unchanged",
                  bool(np.allclose(out, Xsmall)), "T3", "window<=poly"))
    R.append(prop("edge", "SG too-small window -> WARN emitted",
                  any("skipped" in str(w.message).lower() or "savitzky" in str(w.message).lower()
                      for w in wl), "T3", "SPEC.md §5 fix"))

    # 5. MSC flat reference -> no crash / no nan (polyfit is poorly conditioned here by design)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mscflat = cm.Preprocessor(["msc"], cfg).fit_transform(np.ones_like(X))
        ok = bool(np.all(np.isfinite(mscflat)))
    except Exception:
        ok = False
    R.append(prop("edge", "MSC flat reference -> no crash, finite", ok, "T3", "b=0 guard"))

    try:
        import sklearn  # noqa: F401
    except Exception as e:
        R.append(skip("edge", "model caps", f"sklearn unavailable ({type(e).__name__})"))
        return R

    # 6. PLS caps n_components at min(nc, n-1, p)
    m = cm.pls_fit(X, y, 999, cfg)
    R.append(prop("edge", "PLS caps n_components",
                  m["n_components"] <= min(X.shape[0] - 1, X.shape[1]), "T3", "cap min(nc,n-1,p)"))

    # 7. PCA caps n_components at min(n,p)
    pc = cm.pca_fit(X, cfg, n_components=999)
    R.append(prop("edge", "PCA caps n_components", pc["n_components"] <= min(X.shape), "T3", ""))

    # 9. R2 -> nan when total SS == 0 (constant y), via _fit_stats guard
    _, r2, _ = cm._fit_stats(np.full(10, 3.0), np.linspace(0, 1, 10))
    R.append(prop("edge", "R2 == nan when SST=0 (constant y)", bool(np.isnan(r2)), "T3", "guard"))

    return R
