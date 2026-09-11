#!/usr/bin/env python3
"""
Prompt: "Six pectin calibration standards as CSV (DM 3-70 %). Give me a publication figure of the
carbonyl region and the ester/(ester+carboxylate) band-area calibration, with the LOD."

Data: Wang 2023, ATR-FTIR raw spectra, Mendeley Data doi:10.17632/gkwbp3wc49.1 (CC BY 4.0).
"""
DESCRIPTION = ("Pectin DM by ATR-FTIR: ester/(ester+carboxylate) band-area ratio calibration on six "
               "standards (Wang 2023, Mendeley 10.17632/gkwbp3wc49.1)")
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))          # the skill root
sys.path.insert(0, ROOT)

import numpy as np
from scripts.config import Config
from scripts import style, verify, spectra, calibration, report

# ====================================================================== CONFIG
DATA = os.path.join(ROOT, "skill_validation", "ftir_integration", "data", "pectin", "calibration")
STANDARDS = {"DM3": 3.0, "DM20": 20.0, "DM37": 37.0, "DM55": 55.0, "DM62.8": 62.75, "DM70.5": 70.5}  # % DM, known
WINDOW = (1850.0, 1450.0)                                                # the carbonyl region shown
cfg = Config(
    output_dir=os.path.join(HERE, "out"),
    journal="nature", column="double", formats=("pdf", "png"),
    domain="ftir",
    integration_baseline="shared",            # the two bands overlap on one pedestal: ONE baseline, split at 1700
    integration_windows=[(1500.0, 1700.0, "carboxylate"), (1700.0, 1800.0, "ester")],
    conf_level=0.95, lod_loq_method="residual_sd",
    description=DESCRIPTION,
)


def main():
    os.makedirs(cfg.output_dir, exist_ok=True)

    # ---- ingest + gate every file, identically
    traces = {}
    for name in STANDARDS:
        path = os.path.join(DATA, f"{name}.csv")
        verify.check_ingest(path, cfg)
        x, y = spectra.load_xy(path)
        verify.check_trace(x, y, cfg, name=name)
        traces[name] = (np.asarray(x, float), np.asarray(y, float))

    # ---- integrate with the same windows and baseline for every standard
    ratio = []
    for name, (x, y) in traces.items():
        areas = {b["name"]: b["area"] for b in spectra.integrate_bands(x, y, cfg)}
        ratio.append(areas["ester"] / (areas["ester"] + areas["carboxylate"]))
    dm, ratio = np.array(list(STANDARDS.values())), np.array(ratio)

    # ---- calibration + figures of merit
    cal = calibration.fit(dm, ratio, cfg)
    ll = calibration.lod_loq(cal, cfg)
    bias = calibration.intercept_test(cal, cfg)
    print(f"  slope {cal['slope']:.5f} per % DM [{cal['slope_ci'][0]:.5f}, {cal['slope_ci'][1]:.5f}]  "
          f"R2 {cal['r2']:.3f}  n={cal['n']}")
    print(f"  LOD {ll['lod']:.1f} % DM, LOQ {ll['loq']:.1f} % DM  ({ll['method']})")
    print(f"  {bias['caption']}")

    # ---- figure: (a) calibration over its residual strip, (b) the carbonyl window, keyed at the right margin
    style.apply_style(cfg)
    fig, axes = style.figure(cfg, 2, 2, width_ratios=[1.0, 1.15], height_ratios=[3.2, 1.0])
    gs = axes[0, 1].get_gridspec()
    axes[0, 1].remove(); axes[1, 1].remove()
    b = fig.add_subplot(gs[:, 1])
    a, ar = axes[0, 0], axes[1, 0]
    ar.sharex(a)

    xx = np.linspace(*cal["x_range"], 100)
    yy = cal["slope"] * xx + cal["intercept"]
    se_mean = cal["s_resid"] * np.sqrt(1 / cal["n"] + (xx - cal["xbar"]) ** 2 / cal["Sxx"])
    se_pred = cal["s_resid"] * np.sqrt(1 + 1 / cal["n"] + (xx - cal["xbar"]) ** 2 / cal["Sxx"])
    a.fill_between(xx, yy - cal["t"] * se_mean, yy + cal["t"] * se_mean, alpha=0.25, lw=0, label="95 % CI")
    a.plot(xx, yy - cal["t"] * se_pred, lw=0.6, ls="--", color="0.5")
    a.plot(xx, yy + cal["t"] * se_pred, lw=0.6, ls="--", color="0.5", label="95 % prediction")
    a.plot(xx, yy, lw=1.0, label="fit")
    a.scatter(dm, ratio, zorder=3, s=14, label="standards")
    a.set_ylabel("band-area ratio I (dimensionless)")
    a.legend(loc="lower right", fontsize="x-small")
    a.annotate(f"$R^2$ = {cal['r2']:.3f}\nLOD = {ll['lod']:.1f} % DM", xy=(0.96, 0.96),
               xycoords="axes fraction", ha="right", va="top", fontsize="x-small")
    ar.axhline(0, color="0.6", lw=0.6)
    ar.scatter(dm, cal["resid"], zorder=3, s=14)
    ar.set_ylabel("resid. (ratio)"); ar.set_xlabel("DM / %")
    for lab in a.get_xticklabels():
        lab.set_visible(False)

    pal = style._palette(cfg)
    items = []
    for k, (name, (x, y)) in enumerate(traces.items()):
        w = (x <= WINDOW[0]) & (x >= WINDOW[1])
        b.plot(x[w], y[w], lw=0.8, color=pal[k % len(pal)])
        items.append((float(np.interp(WINDOW[1], x[np.argsort(x)], y[np.argsort(x)])),
                      f"DM {STANDARDS[name]:.0f} %", pal[k % len(pal)]))
    for (lo, hi, nm), col in zip(cfg.integration_windows, (pal[1], pal[0])):
        b.axvspan(lo, hi, alpha=0.12, color=col, lw=0)
    b.axvline(1700, color="0.3", lw=0.6, ls="--")
    spectra._apply_axis(b, cfg)
    b.set_xlim(*WINDOW)
    b.set_ylabel("Absorbance (a.u.)")
    spectra.edge_labels(b, items)

    style.finalize_figure(fig, wspace=0.12, hspace=0.06)
    for ax, L in ((a, "a"), (b, "b")):
        style.panel_letter(ax, L)
    verify.audit_layout(fig, cfg)
    style.save_fig(fig, os.path.join(cfg.output_dir, "pectin_calibration"), cfg)

    # ---- caption + methods text (numbers interpolated, never typed)
    caption = (
        f"Figure. (a) Band-area ratio I = A_ester / (A_ester + A_carboxylate) against the degree of "
        f"methyl-esterification for the six standards; line = least-squares fit, shaded = 95 % confidence "
        f"band of the mean, dashed = 95 % prediction band; residuals below. Slope {cal['slope']:.4f} per "
        f"% DM (95 % CI {cal['slope_ci'][0]:.4f} to {cal['slope_ci'][1]:.4f}), R² = {cal['r2']:.3f}, "
        f"n = {cal['n']}; LOD = {ll['lod']:.1f} % DM and LOQ = {ll['loq']:.1f} % DM from 3.3σ/m and "
        f"10σ/m with σ = residual SD. (b) The carbonyl region of the same standards; shaded windows are the "
        f"carboxylate (1500–1700 cm⁻¹) and ester (1700–1800 cm⁻¹) integration windows, integrated on ONE "
        f"shared linear baseline 1500–1800 cm⁻¹ with a vertical drop at 1700 cm⁻¹. Stack order top to "
        f"bottom follows the right-margin key.\n")
    open(os.path.join(HERE, "caption.txt"), "w", encoding="utf-8").write(caption)
    try:
        txt = report.methods_report(cfg, results={"slope": cal["slope"], "r2": cal["r2"], "n": cal["n"],
                                                   "lod": ll["lod"], "loq": ll["loq"]})
        open(os.path.join(HERE, "methods.txt"), "w", encoding="utf-8").write(txt if isinstance(txt, str) else str(txt))
    except Exception as e:                                   # methods text is a convenience, not a gate
        print(f"  [INFO] methods_report skipped: {type(e).__name__}: {e}")
    print("  wrote", os.path.join(cfg.output_dir, "pectin_calibration.pdf"), "+ .png, caption.txt")


if __name__ == "__main__":
    main()
