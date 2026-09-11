"""Preprocessing transforms: SNV, MSC, Savitzky-Golay d1/d2, centre, autoscale.
Definition checks (T0) always run; chemotools parity (T1) skips if absent.
"""
from __future__ import annotations
import numpy as np
from scripts import config, chemometrics as cm
from _harness import close, prop, skip
import _data


def collect():
    R = []
    cfg = config.Config()
    x, X, y, g = _data.synth(seed=1)

    # --- SNV: definition (T0) and chemotools parity (T1)
    snv = cm.Preprocessor(["snv"], cfg).fit_transform(X)
    mu = X.mean(1, keepdims=True)
    sd = X.std(1, keepdims=True)              # ddof=0
    ref = (X - mu) / np.where(sd == 0, 1.0, sd)
    R.append(close("preprocessing", "SNV == (x-mean)/std(ddof0)", snv, ref, "T0",
                   "Barnes 1989 doi:10.1366/0003702894202201; ddof=0"))
    try:
        from chemotools.scatter import StandardNormalVariate
        ct = StandardNormalVariate().fit_transform(X)
        R.append(close("preprocessing", "SNV vs chemotools", snv, ct, "T1",
                       "chemotools.scatter.StandardNormalVariate"))
    except Exception as e:
        R.append(skip("preprocessing", "SNV vs chemotools", f"chemotools unavailable ({type(e).__name__})"))

    # --- MSC: definition (T0) and chemotools parity (T1)
    msc = cm.Preprocessor(["msc"], cfg).fit_transform(X)
    refspec = X.mean(0)
    man = np.empty_like(X)
    for i in range(X.shape[0]):
        b1, b0 = np.polyfit(refspec, X[i], 1)
        man[i] = (X[i] - b0) / (b1 if b1 != 0 else 1.0)
    R.append(close("preprocessing", "MSC == (x-b0)/b1 vs mean ref", msc, man, "T0",
                   "Geladi 1985 doi:10.1366/0003702854248656; mean ref, polyfit deg-1"))
    try:
        from chemotools.scatter import MultiplicativeScatterCorrection
        ctm = MultiplicativeScatterCorrection().fit_transform(X)
        R.append(close("preprocessing", "MSC vs chemotools", msc, ctm, "T1",
                       "chemotools mean-ref polyfit-1"))
    except Exception as e:
        R.append(skip("preprocessing", "MSC vs chemotools", f"chemotools unavailable ({type(e).__name__})"))

    # --- Savitzky-Golay derivatives vs scipy (T0: scipy IS the engine)
    try:
        from scipy.signal import savgol_filter
        win, poly = cfg.deriv_smooth
        w = win if win % 2 == 1 else win + 1
        for order, name in [(1, "d1"), (2, "d2")]:
            ours = cm.Preprocessor([name], cfg).fit_transform(X)
            scp = savgol_filter(X, w, poly, deriv=order, axis=1)
            R.append(close("preprocessing", f"{name} vs scipy.savgol_filter", ours, scp, "T0",
                           "scipy.signal.savgol_filter (engine); Savitzky-Golay 1964"))
    except Exception as e:
        R.append(skip("preprocessing", "SG vs scipy", f"scipy unavailable ({type(e).__name__})"))

    # --- centre (T0)
    cen = cm.Preprocessor(["center"], cfg).fit_transform(X)
    R.append(close("preprocessing", "center == X - col_mean", cen, X - X.mean(0), "T0",
                   "van den Berg 2006 doi:10.1186/1471-2164-7-142"))

    # --- autoscale: definition (T0) and StandardScaler parity (T1)
    auto = cm.Preprocessor(["autoscale"], cfg).fit_transform(X)
    sd0 = X.std(0)
    sd0 = np.where(sd0 == 0, 1.0, sd0)
    R.append(close("preprocessing", "autoscale == (X-mean)/std(ddof0)", auto, (X - X.mean(0)) / sd0, "T0",
                   "ddof=0"))
    try:
        from sklearn.preprocessing import StandardScaler
        ss = StandardScaler().fit_transform(X)
        R.append(close("preprocessing", "autoscale vs StandardScaler", auto, ss, "T1",
                       "sklearn.StandardScaler (ddof=0)"))
    except Exception as e:
        R.append(skip("preprocessing", "autoscale vs StandardScaler", f"sklearn unavailable ({type(e).__name__})"))

    return R
