#!/usr/bin/env python3
"""
validate_ftir_integration.py  -  VALIDATION of the analytical-figures skill against
real, published FTIR quantitative-analysis-by-peak-integration data.

It answers one question: does the skill, run on real literature data, reproduce
the authors' published figures of merit?

Two independent datasets:

  A. PECTIN degree-of-methyl-esterification (DM) by FTIR  -- Wang (2023), Mendeley
     Data doi:10.17632/gkwbp3wc49.1, CC BY 4.0. The ONLY open source that ships
     RAW spectra + the band areas + the calibration, so it validates the actual
     INTEGRATION step, not just the regression. DM is read from a band-area RATIO
        I = A_ester(~1745) / (A_ester + A_carboxylate(~1605))
     -- the internal-standard-ratio idea used for an API/excipient calibration. We
     re-integrate the raw spectra with the skill and re-fit I -> DM.

  B. DOXORUBICIN / ARTEROLANE by mid-IR  -- Bansal, Singh, Kaur, BMC Chemistry 15:27
     (2021), doi:10.1186/s13065-021-00752-3. Publishes per-point peak-AREA tables
     (carbonyl band, baseline-corrected) for transmission and reflectance modes, so
     each calibration is reproducible from the authors' own numbers. DOX's C=O band
     (1770-1680) is the analogue of the aspirin ester C=O.

The raw pectin spectra are NOT vendored (81 MB): download `Raw data.zip` from the Mendeley
record (doi:10.17632/gkwbp3wc49.1, CC BY 4.0) and unzip it so that this path exists:
    skill_validation/ftir_integration/pectin_mendeley/unzipped/Raw data/Calibration/
Part A2 (re-integration) runs only when it does; A1 and B need nothing external.

Run from anywhere (the sys.path shim resolves the skill), or bundle to a standalone with:
    python bundle.py skill_validation/ftir_integration/validate_ftir_integration.py -o standalone_validation.py
"""
import os, sys
# resolve `from scripts import ...` against the installed skill
_SKILL = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the skill root
sys.path.insert(0, _SKILL)

import numpy as np
from scripts.config import Config
from scripts import style, verify, spectra, calibration

# ====================================================================== CONFIG
HERE = os.path.dirname(os.path.abspath(__file__))
PECTIN_DIR = os.path.join(HERE, "pectin_mendeley", "unzipped", "Raw data", "Calibration")

cfg = Config(
    output_dir=os.path.join(HERE, "out"),
    journal="general", column="single", palette="okabe_ito",
    formats=("pdf", "png"),
    domain="ftir",
    # two OVERLAPPING carbonyl bands on a shared pedestal -> shared baseline:
    baseline=None,                         # raw trace; the shared band baseline does the work
    integration_baseline="shared",         # the mode this validation exercises
    integration_windows=[(1500, 1700, "carbox"), (1700, 1800, "ester")],
    conf_level=0.95,
    lod_loq_method="residual_sd",
    strict=False,                          # warn, don't raise, inside a batch validation
)
os.makedirs(cfg.output_dir, exist_ok=True)

# ---- published reference values (transcribed from the two papers) -----------
# Pectin: standard -> (known DM fraction, published A1745, published A1605)
PECTIN = {
    "DM3":    (0.030,  5.32,   14.61),
    "DM20":   (0.200,  4.46,    8.51),
    "DM37":   (0.370, 11.26,   16.44),
    "DM55":   (0.550,  3.2296,  3.7422),
    "DM62.8": (0.6275, 5.718,   6.133),
    "DM70.5": (0.705, 18.065,  17.179),
}
# a few of the authors' worked unknowns: measured I -> their reported DM
PECTIN_SAMPLES = [
    ("W-HWP-1", 0.47613390672147116, 0.591998448525743),
    ("CA-HWP-1", 0.304537934294493,  0.1106202670763411),
    ("E-SWP-1", 0.3810848400556328,  0.32535730180806666),
]
# BMC 2021: label -> (conc list, area list, published "slope,intercept,r,LOD,LOQ")
BMC = {
    "DOX transmittance": ([0.6, 0.8, 1.0, 1.2, 1.4], [2.131, 2.850, 3.521, 3.913, 4.601],
                          (3.0, 0.402, 0.992, 0.11, 0.30)),
    "DOX reflectance":   ([0.6, 0.8, 1.0, 1.2, 1.4], [0.614, 0.843, 1.03, 1.201, 1.382],
                          (0.95, 0.062, 0.996, 0.15, 0.48)),
    "ALM transmittance": ([0.2, 0.4, 0.6, 0.8, 1.0], [0.644, 0.811, 1.091, 1.40, 1.581],
                          (1.235, 0.363, 0.990, 0.07, 0.19)),
    "ALM reflectance":   ([0.2, 0.4, 0.6, 0.8, 1.0], [0.180, 0.486, 0.811, 1.060, 1.272],
                          (1.38, -0.068, 0.993, 0.06, 0.18)),
}


def _ratio_from_areas(a_ester, a_carbox):
    return a_ester / (a_ester + a_carbox)


def validate_pectin_published():
    """A1: reproduce the authors' calibration from their PUBLISHED band areas
    (isolates the regression / prediction math from any integration choice)."""
    print("\n" + "=" * 72)
    print("A1. PECTIN calibration from the authors' PUBLISHED areas (exact target)")
    print("=" * 72)
    DM = np.array([v[0] for v in PECTIN.values()])
    I = np.array([_ratio_from_areas(v[1], v[2]) for v in PECTIN.values()])

    # classical ICH orientation: x = known property (DM), y = measured response (I)
    model = calibration.fit(DM, I, cfg)
    print(f"  fit  I = {model['slope']:.4f}*DM + {model['intercept']:.4f}   "
          f"R^2 = {model['r2']:.5f}")
    print(f"  LOD/LOQ (in DM) = {calibration.lod_loq(model, cfg)['lod']:.3f} / "
          f"{calibration.lod_loq(model, cfg)['loq']:.3f}")

    # the authors actually read DM off an inverse fit DM = m*I + b; reproduce THAT
    inv_s, inv_b = np.polyfit(I, DM, 1)
    print(f"  authors' inverse fit  DM = {inv_s:.4f}*I + {inv_b:.4f}  (used to read unknowns)")
    print("  reproduce their worked unknowns from measured I:")
    worst = 0.0
    for name, I_meas, dm_pub in PECTIN_SAMPLES:
        dm_pred = inv_s * I_meas + inv_b
        worst = max(worst, abs(dm_pred - dm_pub))
        print(f"    {name:9s} I={I_meas:.5f}  DM_pred={dm_pred:.5f}  DM_paper={dm_pub:.5f}  "
              f"d={dm_pred - dm_pub:+.2e}")
    print(f"  -> max |DM_pred - DM_paper| over unknowns = {worst:.1e}  (rounding level)")

    fig, _ = calibration.plot_calibration(DM, I, cfg, model)
    fig.axes[0].set_ylabel("band-area ratio  I")
    fig.axes[-1].set_xlabel("DM / mole-ester fraction")
    verify.audit_layout(fig, cfg)
    style.save_fig(fig, os.path.join(cfg.output_dir, "pectin_calibration_published"), cfg)
    return model


def validate_pectin_integration():
    """A2: re-INTEGRATE the raw spectra with the skill and compare to the authors'
    areas/ratio, under both baseline modes -- the actual peak-integration test."""
    print("\n" + "=" * 72)
    print("A2. PECTIN re-integration of the RAW spectra with the skill")
    print("=" * 72)

    def run(mode):
        cfg.integration_baseline = mode
        DM, I_mine, I_pub = [], [], []
        for std, (dm, a1, a2) in PECTIN.items():
            path = os.path.join(PECTIN_DIR, std, f"{std}-raw data.CSV")
            x, y = spectra.load_xy(path)
            verify.check_trace(x, y, cfg, name=std)
            bands = {b["name"]: b["area"] for b in spectra.integrate_bands(x, y, cfg)}
            DM.append(dm)
            I_mine.append(_ratio_from_areas(bands["ester"], bands["carbox"]))
            I_pub.append(_ratio_from_areas(a1, a2))
        DM, I_mine, I_pub = map(np.array, (DM, I_mine, I_pub))
        m = calibration.fit(DM, I_mine, cfg)
        return m["r2"], float(np.mean(np.abs(I_mine - I_pub))), I_mine, I_pub

    r2_pw, mad_pw, *_ = run("per_window")
    r2_sh, mad_sh, Imine, Ipub = run("shared")
    print(f"  per_window baseline :  calib R^2 = {r2_pw:.4f}   mean|I_mine-I_pub| = {mad_pw:.4f}")
    print(f"  shared baseline     :  calib R^2 = {r2_sh:.4f}   mean|I_mine-I_pub| = {mad_sh:.4f}")
    print(f"  authors' published-area calibration R^2 = 0.9935 (target)")
    print("  per-standard I (shared mode) vs published:")
    for std, im, ip in zip(PECTIN, Imine, Ipub):
        print(f"    {std:7s} I_mine={im:.4f}  I_pub={ip:.4f}  d={im - ip:+.4f}")
    print("  NOTE: exact match needs the authors' (graphically-specified) baseline\n"
          "        anchors + valley split; shared mode is the correct method and\n"
          "        removes the per_window pedestal artefact (see DM3).")
    cfg.integration_baseline = "shared"


def validate_bmc():
    """B: reproduce the four BMC-2021 calibrations from the published area tables."""
    print("\n" + "=" * 72)
    print("B. BMC Chemistry 2021 (doxorubicin / arterolane) -- skill vs published")
    print("=" * 72)
    first_fig = None
    for name, (xs, ys, (m_pub, b_pub, r_pub, lod_pub, loq_pub)) in BMC.items():
        m = calibration.fit(xs, ys, cfg)
        ll = calibration.lod_loq(m, cfg)
        print(f"  {name}:")
        print(f"     skill     : y = {m['slope']:.3f}x {m['intercept']:+.3f}   "
              f"r = {np.sqrt(m['r2']):.4f}   LOD/LOQ = {ll['lod']:.3f}/{ll['loq']:.3f}")
        print(f"     published : y = {m_pub:.3f}x {b_pub:+.3f}   "
              f"r = {r_pub:.4f}   LOD/LOQ = {lod_pub:.3f}/{loq_pub:.3f}")
        if first_fig is None:
            fig, _ = calibration.plot_calibration(np.array(xs), np.array(ys), cfg, m)
            fig.axes[0].set_ylabel("peak area / mm$^2$")
            fig.axes[-1].set_xlabel("doxorubicin / % w/w")
            verify.audit_layout(fig, cfg)
            style.save_fig(fig, os.path.join(cfg.output_dir, "bmc_dox_transmittance"), cfg)
            first_fig = True
    print("  -> slopes/intercepts/linearity reproduce to the published rounding.\n"
          "     LOD/LOQ matches only for DOX-transmittance; the other three imply a\n"
          "     band-specific sigma the paper did not document (skill names its sigma).")


def main():
    validate_pectin_published()
    if os.path.isdir(PECTIN_DIR):
        validate_pectin_integration()
    else:
        print(f"\nA2 skipped: raw pectin spectra not found at {PECTIN_DIR}\n"
              "   (download doi:10.17632/gkwbp3wc49.1 'Raw data.zip' and unzip there)")
    validate_bmc()
    print("\nfigures written to", cfg.output_dir)


if __name__ == "__main__":
    main()
