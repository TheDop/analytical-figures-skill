"""
calibration.py  -  linear calibration done honestly, ICH Q2 flavour.

The correctness traps this guards against:
  * R^2 alone is NOT evidence of linearity -> a residual panel is MANDATORY.
  * Reading concentrations off the curve outside the calibrated range
    (extrapolation) -> predict() refuses it.
  * Quoting LOD/LOQ without saying how -> the method is named in the result.

    fit(x, y, cfg)                  -> dict(slope, intercept, r2, s_resid, n, ...)
    predict(model, signal, cfg)     -> concentration (raises on extrapolation)
    lod_loq(model, cfg)             -> dict(lod, loq, method)
    plot_calibration(x, y, cfg)     -> figure with curve + bands + MANDATORY residuals
"""
from __future__ import annotations
import numpy as np
from . import style, verify


_CAL_NO_SCIPY_WARNED = False


def _tmult(df, conf):
    try:
        from scipy import stats
        return float(stats.t.ppf(0.5 + conf / 2, df))
    except Exception:
        from statistics import NormalDist              # scipy absent: exact z, and say so (once)
        z = float(NormalDist().inv_cdf(min(max(0.5 + conf / 2, 1e-12), 1 - 1e-12)))
        global _CAL_NO_SCIPY_WARNED
        if not _CAL_NO_SCIPY_WARNED:
            _CAL_NO_SCIPY_WARNED = True
            print(f"  [WARN] calibration: scipy unavailable - normal quantile z={z:.3f} used instead of "
                  f"t at {df} dof (bands too narrow at small n); further calls stay silent")
        return z


def fit(x, y, cfg):
    """Ordinary least squares y = slope*x + intercept, with the statistics needed
    for honest bands and a slope CI."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    verify.check_trace(np.sort(x), y[np.argsort(x)], cfg, name="calibration",
                       min_points=3, require_monotonic=False)
    n = x.size
    xbar = x.mean()
    Sxx = float(np.sum((x - xbar) ** 2))
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    resid = y - yhat
    dof = n - 2
    s_resid = float(np.sqrt(np.sum(resid ** 2) / dof)) if dof > 0 else float("nan")
    se_slope = s_resid / np.sqrt(Sxx) if Sxx else float("nan")
    se_intercept = s_resid * np.sqrt(1.0 / n + xbar ** 2 / Sxx) if Sxx else float("nan")
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - np.sum(resid ** 2) / ss_tot if ss_tot else float("nan")
    t = _tmult(dof, cfg.conf_level)
    return {
        "slope": float(slope), "intercept": float(intercept), "r2": float(r2),
        "s_resid": s_resid, "se_slope": se_slope, "se_intercept": se_intercept,
        "n": int(n), "dof": int(dof),
        "xbar": xbar, "Sxx": Sxx, "t": t, "conf": cfg.conf_level,
        "x_range": (float(x.min()), float(x.max())),
        "x": x, "y": y, "resid": resid,
        "slope_ci": (float(slope - t * se_slope), float(slope + t * se_slope)),
        "intercept_ci": (float(intercept - t * se_intercept), float(intercept + t * se_intercept)),
    }


def predict(model, signal, cfg, x_for_signal=None):
    """Invert the curve: concentration from a measured signal. Refuses to read
    outside the calibrated x range when cfg.refuse_extrapolation (the default)."""
    conc = (signal - model["intercept"]) / model["slope"]
    lo, hi = model["x_range"]
    if cfg.refuse_extrapolation and not (lo <= conc <= hi):
        raise verify.GateError(
            f"predicted x={conc:.4g} is outside the calibrated range [{lo:.4g}, {hi:.4g}] "
            f"- extrapolation refused (extend the calibration instead)")
    return float(conc)


def lod_loq(model, cfg):
    """LOD = 3.3*sigma/slope, LOQ = 10*sigma/slope. With 'residual_sd' (default)
    sigma is the regression residual SD; the method is recorded in the result."""
    if cfg.lod_loq_method == "residual_sd":
        sigma = model["s_resid"]
    else:
        raise ValueError("blank_sd method needs blank replicates passed separately")
    m = abs(model["slope"])
    return {"lod": 3.3 * sigma / m, "loq": 10.0 * sigma / m,
            "method": f"{cfg.lod_loq_method} (3.3σ/m, 10σ/m)"}


def intercept_test(model, cfg):
    """Test whether the intercept differs significantly from zero — i.e. whether the
    method carries a CONSTANT BIAS. Two-sided t-test on b/SE(b) at (n-2) dof.

    Returns {intercept, se, t, p, ci, significant, caption}. A CI that INCLUDES 0
    (equivalently p > 1-conf) means no significant constant bias — the line is
    consistent with passing through the origin. A significant intercept flags a
    constant offset (baseline, matrix, or low-end lack-of-fit) to investigate."""
    b, se, dof = model["intercept"], model["se_intercept"], model["dof"]
    t_stat = b / se if se else float("nan")
    try:
        from scipy import stats
        p = float(2 * stats.t.sf(abs(t_stat), dof))
    except Exception:
        p = float("nan")
    lo, hi = model["intercept_ci"]
    significant = not (lo <= 0.0 <= hi)
    verdict = ("constant bias: intercept significantly != 0 (CI excludes 0)" if significant
               else "no significant constant bias (CI includes 0)")
    cap = (f"intercept {b:.4g} [{lo:.4g}, {hi:.4g}] ({int(model['conf']*100)}% CI), "
           f"t={t_stat:.2f}, p={p:.3f} -> {verdict}")
    return {"intercept": float(b), "se": float(se), "t": float(t_stat), "p": p,
            "ci": (lo, hi), "significant": bool(significant), "caption": cap}


def plot_calibration(x, y, cfg, model=None):
    """Calibration curve with confidence + prediction bands AND a mandatory
    residual panel underneath (R^2 alone does not certify linearity)."""
    if model is None:
        model = fit(x, y, cfg)
    style.apply_style(cfg)

    if cfg.force_residual_panel:
        fig, (ax, axr) = style.figure(cfg, nrows=2, height=None,
                                      height_ratios=[3, 1], sharex=True)
    else:
        fig, ax = style.figure(cfg); axr = None

    xs = np.linspace(*model["x_range"], 200)
    ys = model["slope"] * xs + model["intercept"]
    # CI of the mean response, and wider prediction band
    s, t, Sxx, xbar, n = (model["s_resid"], model["t"], model["Sxx"],
                          model["xbar"], model["n"])
    se_mean = s * np.sqrt(1.0 / n + (xs - xbar) ** 2 / Sxx)
    se_pred = s * np.sqrt(1.0 + 1.0 / n + (xs - xbar) ** 2 / Sxx)

    ax.scatter(model["x"], model["y"], zorder=3, label="data")
    ax.plot(xs, ys, label="fit")
    ax.fill_between(xs, ys - t * se_mean, ys + t * se_mean, alpha=0.25, lw=0,
                    label=f"{int(model['conf']*100)}% CI")
    ax.plot(xs, ys - t * se_pred, lw=0.6, ls="--", color="0.5")
    ax.plot(xs, ys + t * se_pred, lw=0.6, ls="--", color="0.5",
            label=f"{int(model['conf']*100)}% pred.")
    ax.set_ylabel("Response (y units)")
    ax.legend(loc="best")
    # report R^2 and slope CI in a corner annotation, not as a substitute for residuals
    lo, hi = model["slope_ci"]
    ax.annotate(f"$R^2$={model['r2']:.4f}\nslope {model['slope']:.3g} "
                f"[{lo:.3g}, {hi:.3g}]",
                xy=(0.03, 0.97), xycoords="axes fraction", va="top", fontsize="small")

    if axr is not None:
        axr.axhline(0, color="0.6", lw=0.6)
        axr.scatter(model["x"], model["resid"], zorder=3)
        axr.set_ylabel("Resid. (y units)")
        axr.set_xlabel("Concentration (x units)")
    else:
        ax.set_xlabel("Concentration (x units)")
    style.finalize_figure(fig)
    return fig, (ax, axr)
