#!/usr/bin/env python3
"""
Prompt: "Doxorubicin peak areas from Bansal 2021, transmittance mode, 0.6-1.4 % w/w. Fit the
calibration, give me LOD and LOQ, and read off the concentration for peak areas of 3.5 and 6.0 mm2."

Data: Bansal, Singh & Kaur, BMC Chemistry 15:27 (2021), doi:10.1186/s13065-021-00752-3 — the
published per-point peak-area table (carbonyl band, baseline-corrected, transmittance mode).
"""
DESCRIPTION = ("Doxorubicin mid-IR calibration from the published carbonyl peak-area table "
               "(Bansal 2021, BMC Chem 15:27): fit, LOD/LOQ, and read-back of two unknowns")
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))          # the skill root
sys.path.insert(0, ROOT)

import numpy as np
from scripts.config import Config
from scripts import style, verify, calibration

# ====================================================================== CONFIG
CONC = [0.6, 0.8, 1.0, 1.2, 1.4]                 # % w/w      (Bansal 2021, Table: DOX transmittance)
AREA = [2.131, 2.850, 3.521, 3.913, 4.601]       # peak area / mm2
UNKNOWNS = [3.5, 6.0]                            # peak areas the user wants concentrations for
cfg = Config(
    output_dir=os.path.join(HERE, "out"),
    journal="nature", column="single", formats=("pdf", "png"),
    conf_level=0.95, lod_loq_method="residual_sd",
    refuse_extrapolation=True,                   # the default; shown here because the second unknown tests it
    description=DESCRIPTION,
)


def main():
    os.makedirs(cfg.output_dir, exist_ok=True)
    x, y = np.array(CONC), np.array(AREA)

    cal = calibration.fit(x, y, cfg)
    ll = calibration.lod_loq(cal, cfg)
    bias = calibration.intercept_test(cal, cfg)
    print(f"  y = {cal['slope']:.3f} x {cal['intercept']:+.3f}   R2 {cal['r2']:.4f}   n={cal['n']}")
    print(f"  LOD {ll['lod']:.3f} % w/w, LOQ {ll['loq']:.3f} % w/w  ({ll['method']})")
    print(f"  {bias['caption']}")

    # ---- read back the unknowns; the gate refuses anything outside the calibrated range
    readback = {}
    for area in UNKNOWNS:
        try:
            conc = calibration.predict(cal, area, cfg)
            readback[area] = conc
            print(f"  peak area {area:.2f} mm2 -> {conc:.3f} % w/w")
        except verify.GateError as e:
            readback[area] = None
            print(f"  peak area {area:.2f} mm2 -> REFUSED: {e}")

    # ---- the figure: curve + bands + the mandatory residual panel
    fig, (ax, axr) = calibration.plot_calibration(x, y, cfg, cal)
    ax.set_ylabel("peak area / mm²")
    axr.set_ylabel("resid. / mm²")
    axr.set_xlabel("doxorubicin / % w/w")
    verify.audit_layout(fig, cfg)
    style.save_fig(fig, os.path.join(cfg.output_dir, "bmc_dox_calibration"), cfg)

    ok = [f"{a:.2f} mm² → {c:.3f} % w/w" for a, c in readback.items() if c is not None]
    refused = [f"{a:.2f} mm²" for a, c in readback.items() if c is None]
    caption = (
        f"Figure. Doxorubicin calibration, carbonyl peak area (transmittance mode) against concentration, "
        f"from the published area table of Bansal et al. 2021. Line = least-squares fit "
        f"y = {cal['slope']:.3f}x {cal['intercept']:+.3f}, R² = {cal['r2']:.4f}, n = {cal['n']}; shaded = 95 % "
        f"confidence band, dashed = 95 % prediction band; residuals below. LOD = {ll['lod']:.3f} % w/w and "
        f"LOQ = {ll['loq']:.3f} % w/w (3.3σ/m and 10σ/m, σ = residual SD of the fit). Read-back: "
        f"{'; '.join(ok)}. Not reported: {', '.join(refused)} — outside the calibrated range "
        f"{cal['x_range'][0]:.1f}–{cal['x_range'][1]:.1f} % w/w, so no value is extrapolated.\n")
    open(os.path.join(HERE, "caption.txt"), "w", encoding="utf-8").write(caption)
    print("  wrote", os.path.join(cfg.output_dir, "bmc_dox_calibration.pdf"), "+ .png, caption.txt")


if __name__ == "__main__":
    main()
