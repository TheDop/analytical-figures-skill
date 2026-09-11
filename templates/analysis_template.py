#!/usr/bin/env python3
"""
analysis_template.py  -  COPY THIS to start an analysis.

The shape to preserve: ONE Config block at the top is the only thing you edit
for routine work. Everything below it is the gated workflow:
    ingest -> gate the numbers -> process identically -> plot -> visual QA loop
    -> export vector master + raster proof.

Run from the skill root so `from scripts import ...` resolves, or adjust sys.path.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from scripts.config import Config
from scripts import style, verify, spectra, calibration, charts, report

# ====================================================================== CONFIG
# The one tunable block. If you want to change something not here, add it to
# Config first rather than burying it in the code below.
cfg = Config(
    output_dir="./out",
    journal="nature", column="single", palette="okabe_ito",
    formats=("pdf", "png"),
    domain="ftir",
    x_limits=(1820, 1650),                       # already inverted for FTIR
    baseline="arpls", arpls_lam=1e5,
    integration_windows=[(1745, 1700, "C=O")],   # identical for every spectrum
    normalize="max",
    conf_level=0.95,
)
os.makedirs(cfg.output_dir, exist_ok=True)


def demo():
    # --- 1. synth two FTIR reps so the template runs standalone -------------
    x = np.linspace(1650, 1820, 400)
    def rep(scale):
        peak = scale * np.exp(-0.5 * ((x - 1720) / 8) ** 2)
        drift = 0.03 * (x - x.min()) / np.ptp(x)          # baseline to remove
        return peak + drift + np.random.normal(0, 0.002, x.size)
    reps = [(x, rep(1.00)), (x, rep(0.97))]

    # --- 2. process IDENTICALLY + record it, then assert sameness -----------
    proc, records, disp = [], [], []
    for xr, yr in reps:
        verify.check_trace(xr, yr, cfg, name="ftir rep")
        yc, _ = spectra.correct_baseline(xr, yr, cfg)
        yn = spectra.normalize(xr, yc, cfg)
        disp.append((xr, yn))
        records.append({"baseline": cfg.baseline, "arpls_lam": cfg.arpls_lam,
                        "windows": cfg.integration_windows, "normalize": cfg.normalize})
        proc.append(spectra.integrate_bands(xr, yc, cfg))
    verify.same_processing(records, cfg)

    areas = [p[0]["area"] for p in proc]
    s = verify.summarize(areas, cfg)
    print(f"  C=O band area: {s['mean']:.4g}  ({s['caption']})")

    # --- 3. plot, then the VISUAL QA LOOP -----------------------------------
    fig, ax = spectra.plot_overlay({"rep 1": disp[0], "rep 2": disp[1]}, cfg)
    prev = verify.render_preview(fig, os.path.join(cfg.output_dir, "_preview_ftir.png"))
    verify.audit_layout(fig, cfg)
    # >>> here the model OPENS prev with the Read tool and walks
    #     verify.READ_IMAGE_CHECKLIST, fixes anything, and re-renders. <<<
    style.save_fig(fig, os.path.join(cfg.output_dir, "ftir_overlay"), cfg)

    # --- 4. a calibration curve (mandatory residual panel) ------------------
    xc = np.array([1., 2., 4., 8., 16.])
    yc = 0.42 * xc + 0.05 + np.random.normal(0, 0.05, xc.size)
    model = calibration.fit(xc, yc, cfg)
    print(f"  calib R^2={model['r2']:.4f}, slope CI={model['slope_ci']}")
    print(f"  LOD/LOQ: {calibration.lod_loq(model, cfg)}")
    figc, _ = calibration.plot_calibration(xc, yc, cfg, model)
    verify.render_preview(figc, os.path.join(cfg.output_dir, "_preview_calib.png"))
    verify.audit_layout(figc, cfg)
    style.save_fig(figc, os.path.join(cfg.output_dir, "calibration"), cfg)

    # --- 5. paste-ready methods text + SI table (form, not findings) --------
    lo, hi = model["slope_ci"]
    methods, si = report.methods_report(cfg, results={
        "C=O area (mean)": s["mean"], "C=O area %CI": s["ci_halfwidth"],
        "calib R^2": model["r2"], "slope": model["slope"],
        "LOD": calibration.lod_loq(model, cfg)["lod"],
    })
    print("\n--- METHODS (paste-ready) ---\n" + methods)
    print("\n--- SI TABLE ---\n" + si)

    print("\n  template OK ->", cfg.output_dir)


if __name__ == "__main__":
    demo()
