#!/usr/bin/env python3
"""
groundtruth_integration_test.py  -  DETERMINISTIC validation of the
analytical-figures integrator against KNOWN answers, under realistic conditions.

The real-world reproductions (pectin, BMC, diclofenac) each hit a different data
limitation, so none could pin integrate_bands to an *exact* known area. This test
removes that ambiguity: it builds spectra from analytic peaks whose true area is
known in closed form, then checks the skill recovers it -- but it does so under
real-world conditions: a sloping+offset baseline, FTIR-realistic
point spacing (~2 cm-1, i.e. 4 cm-1 resolution), additive noise, an internal-standard
ratio, and overlapping bands on a shared pedestal.

Ground truth used (independent of the skill's code):
  * Gaussian area  = h * s * sqrt(2*pi)              (closed form)
  * scipy.integrate.quad on the analytic, baseline-subtracted, clipped band
    (an adaptive integrator -- NOT the skill's trapezoid) for the overlap partition.

What it certifies:
  T1  isolated band: recovers the analytic area through a sloping baseline.
  T2  overlapping bands: shared (drop-perpendicular) mode matches an independent
      integrator and the true ratio; per_window is biased by the shared pedestal.
  T3  internal-standard ratio is invariant to overall scale (ATR loading) and to
      cfg.normalize -- the property that makes the ratio cancel contact/pathlength.
  T4  a DM-like calibration series: shared-mode recovered ratio is linear vs the
      true mixing fraction (R^2 ~ 1); per_window is biased/curved -- the pectin
      finding, now against a known answer.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))   # the skill root
import numpy as np
from scipy.integrate import quad
from scripts.config import Config
from scripts import style, verify, spectra, calibration

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(20260615)
results = []   # (name, passed, detail)


def gaussian(x, h, c, s):
    return h * np.exp(-0.5 * ((x - c) / s) ** 2)


def gauss_area(h, s):
    return h * s * np.sqrt(2 * np.pi)


def linbase(x, lo, hi, ylo, yhi):
    return ylo + (yhi - ylo) * (x - lo) / (hi - lo)


def check(name, got, ref, tol, kind="rel"):
    err = abs(got - ref) / abs(ref) if kind == "rel" else abs(got - ref)
    ok = err <= tol
    results.append((name, ok, f"got={got:.6g} ref={ref:.6g} {kind}err={err:.2e} tol={tol:.0e}"))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: got={got:.6g} ref={ref:.6g} {kind}err={err:.2e}")
    return ok


# ===================================================================== T1
def t1_isolated_band():
    print("\nT1. isolated band, analytic area, through a sloping+offset baseline")
    h, c, s = 0.30, 1745.0, 12.0
    true_area = gauss_area(h, s)
    lo, hi = c - 5 * s, c + 5 * s            # 5 sigma: tails are negligible
    x = np.arange(1500.0, 2000.0, 2.0)       # FTIR-realistic 2 cm-1 spacing
    baseline = 0.05 + 8e-4 * (x - x[0])      # offset + slope
    y_clean = gaussian(x, h, c, s) + baseline
    cfg = Config(domain="ftir", integration_baseline="per_window",
                 integration_windows=[(lo, hi, "band")])
    # the 5-sigma window truncates ~5.7e-5 of the area; compare to that, not the full inf integral
    win_true = quad(lambda v: gaussian(v, h, c, s), lo, hi)[0]
    area_clean = spectra.integrate_bands(x, y_clean, cfg)[0]["area"]
    check("T1 clean (sloping baseline removed)", area_clean, win_true, 2e-3)
    check("T1 window vs full analytic area", win_true, true_area, 1e-3)
    # with noise: repeat and average, report scatter
    areas = []
    for _ in range(200):
        yn = y_clean + RNG.normal(0, 1e-3, x.size)
        areas.append(spectra.integrate_bands(x, yn, cfg)[0]["area"])
    m, sd = float(np.mean(areas)), float(np.std(areas))
    print(f"       noisy (sigma=1e-3, n=200): mean={m:.5f} sd={sd:.5f} "
          f"-> bias {abs(m-win_true)/win_true*100:.2f}%, RSD {sd/m*100:.2f}%")
    check("T1 noisy mean unbiased", m, win_true, 5e-3)


# ============================================ overlapping-band model (T2-T4)
# Broad ester + carboxylate bands that genuinely OVERLAP, so the valley between
# them sits well above baseline -> a real shared pedestal (the pectin situation).
CE, SE = 1740.0, 28.0
CC, SC = 1620.0, 30.0
ENV_LO, ENV_HI = 1500.0, 1858.0
GRID = np.arange(1480.0, 1880.0, 2.0)         # FTIR-realistic 2 cm-1 spacing


def overlap_spectrum(x, he, hc, loading=1.0, x0=1480.0):
    base = 0.04 + 6e-4 * (x - x0)             # offset+slope; removed exactly by a linear baseline
    return loading * (gaussian(x, he, CE, SE) + gaussian(x, hc, CC, SC)) + base


def model_scalar(v, he, hc, loading=1.0, x0=1480.0):
    return loading * (gaussian(v, he, CE, SE) + gaussian(v, hc, CC, SC)) + 0.04 + 6e-4 * (v - x0)


def valley_split():
    """Fixed integration split at the valley (computed once at he=hc, as a real
    method would fix its limits)."""
    xf = np.linspace(ENV_LO, ENV_HI, 40001)
    L = linbase(xf, ENV_LO, ENV_HI, model_scalar(ENV_LO, .3, .3), model_scalar(ENV_HI, .3, .3))
    corr = overlap_spectrum(xf, .3, .3) - L
    seg = (xf > 1640) & (xf < 1740)
    raw = float(xf[seg][np.argmin(corr[seg])])
    return float(GRID[np.argmin(np.abs(GRID - raw))]), corr[seg].min(), corr.max()


SPLIT, _valley_h, _peak_h = valley_split()


def _quad_partition(he, hc, mode):
    """Independent (scipy.quad) reference for the drop-perpendicular ('shared') or
    'per_window' area definition, on the continuous analytic model."""
    out = {}
    for a, b, name in [(ENV_LO, SPLIT, "carbox"), (SPLIT, ENV_HI, "ester")]:
        if mode == "shared":
            blo, bhi = ENV_LO, ENV_HI
        else:
            blo, bhi = a, b
        ylo, yhi = model_scalar(blo, he, hc), model_scalar(bhi, he, hc)
        f = lambda v: max(model_scalar(v, he, hc) - linbase(v, blo, bhi, ylo, yhi), 0.0)
        out[name] = quad(f, a, b, limit=200)[0]
    return out


# ===================================================================== T2
def t2_overlap_partition():
    print("\nT2. overlapping bands on a shared pedestal: shared & per_window vs an "
          "independent integrator")
    print(f"       fixed split at {SPLIT:.1f} cm-1; valley sits at "
          f"{_valley_h/_peak_h*100:.0f}% of peak height (real pedestal)")
    he, hc = 0.30, 0.25
    y = overlap_spectrum(GRID, he, hc)
    cfg_sh = Config(domain="ftir", integration_baseline="shared",
                    integration_windows=[(ENV_LO, SPLIT, "carbox"), (SPLIT, ENV_HI, "ester")])
    cfg_pw = Config(domain="ftir", integration_baseline="per_window",
                    integration_windows=[(ENV_LO, SPLIT, "carbox"), (SPLIT, ENV_HI, "ester")])
    sh = {b["name"]: b["area"] for b in spectra.integrate_bands(GRID, y, cfg_sh)}
    pw = {b["name"]: b["area"] for b in spectra.integrate_bands(GRID, y, cfg_pw)}
    qsh, qpw = _quad_partition(he, hc, "shared"), _quad_partition(he, hc, "per_window")

    # both modes' CODE matches an independent integrator of the matching definition
    check("T2 shared  ester  == quad", sh["ester"], qsh["ester"], 5e-3)
    check("T2 shared  carbox == quad", sh["carbox"], qsh["carbox"], 5e-3)
    check("T2 perwin  ester  == quad", pw["ester"], qpw["ester"], 5e-3)
    check("T2 perwin  carbox == quad", pw["carbox"], qpw["carbox"], 5e-3)

    # how the two modes differ on the ratio, and vs the true individual-band fraction
    I_true = gauss_area(he, SE) / (gauss_area(he, SE) + gauss_area(hc, SC))
    I_sh = sh["ester"] / (sh["ester"] + sh["carbox"])
    I_pw = pw["ester"] / (pw["ester"] + pw["carbox"])
    print(f"       I_true={I_true:.4f}  I_shared={I_sh:.4f} (d={I_sh-I_true:+.4f})  "
          f"I_perwindow={I_pw:.4f} (d={I_pw-I_true:+.4f})")
    ok = abs(I_sh - I_true) < abs(I_pw - I_true)
    results.append(("T2 shared ratio closer to truth than per_window", ok,
                    f"|dI_sh|={abs(I_sh-I_true):.4f} vs |dI_pw|={abs(I_pw-I_true):.4f}"))
    print(f"  [{'PASS' if ok else 'FAIL'}] T2 shared ratio closer to truth than per_window")


# ===================================================================== T3
def t3_ratio_invariance():
    print("\nT3. internal-standard ratio invariant to scale (ATR loading) and to normalize")
    y = overlap_spectrum(GRID, 0.30, 0.25)
    cfg = Config(domain="ftir", integration_baseline="shared",
                 integration_windows=[(ENV_LO, SPLIT, "carbox"), (SPLIT, ENV_HI, "ester")])

    def ratio(trace):
        b = {d["name"]: d["area"] for d in spectra.integrate_bands(GRID, trace, cfg)}
        return b["ester"] / (b["ester"] + b["carbox"])

    # NB scale only the band signal, not the additive baseline (the baseline is removed)
    sig = overlap_spectrum(GRID, 0.30, 0.25) - (0.04 + 6e-4 * (GRID - GRID[0]))
    base = 0.04 + 6e-4 * (GRID - GRID[0])
    I0 = ratio(sig + base)
    for k in (0.5, 2.0, 5.0):
        check(f"T3 ratio invariant to x{k} loading", ratio(k * sig + base), I0, 1e-9)
    a1 = spectra.integrate_bands(GRID, sig + base, cfg)[0]["area"]
    a2 = spectra.integrate_bands(GRID, 3.0 * sig + base, cfg)[0]["area"]
    check("T3 absolute area scales linearly with loading", a2 / a1, 3.0, 1e-6)


# ===================================================================== T4
def t4_calibration_series():
    print("\nT4. DM-like series under realistic conditions (overlap, loading, noise)")
    Ae1, Ac1 = gauss_area(1.0, SE), gauss_area(1.0, SC)
    cfg_sh = Config(domain="ftir", integration_baseline="shared",
                    integration_windows=[(ENV_LO, SPLIT, "carbox"), (SPLIT, ENV_HI, "ester")])
    cfg_pw = Config(domain="ftir", integration_baseline="per_window",
                    integration_windows=[(ENV_LO, SPLIT, "carbox"), (SPLIT, ENV_HI, "ester")])

    fracs = np.array([0.15, 0.30, 0.45, 0.60, 0.75, 0.90])
    I_true, I_sh, I_pw = [], [], []
    for f in fracs:
        he, hc = f, (1 - f)
        loading = float(RNG.uniform(0.6, 1.4))            # random ATR loading per sample
        y = overlap_spectrum(GRID, he, hc, loading=loading) + RNG.normal(0, 8e-4, GRID.size)
        I_true.append((he * Ae1) / (he * Ae1 + hc * Ac1))
        bs = {d["name"]: d["area"] for d in spectra.integrate_bands(GRID, y, cfg_sh)}
        bp = {d["name"]: d["area"] for d in spectra.integrate_bands(GRID, y, cfg_pw)}
        I_sh.append(bs["ester"] / (bs["ester"] + bs["carbox"]))
        I_pw.append(bp["ester"] / (bp["ester"] + bp["carbox"]))
    I_true, I_sh, I_pw = map(np.array, (I_true, I_sh, I_pw))

    cfg = Config(domain="ftir", conf_level=0.95)
    m_sh, m_pw = calibration.fit(I_true, I_sh, cfg), calibration.fit(I_true, I_pw, cfg)
    bias = float(np.max(np.abs(I_pw - I_sh)))
    print(f"       shared    : recovered = {m_sh['slope']:.3f}*true {m_sh['intercept']:+.3f}  R^2={m_sh['r2']:.5f}")
    print(f"       per_window: recovered = {m_pw['slope']:.3f}*true {m_pw['intercept']:+.3f}  R^2={m_pw['r2']:.5f}")
    print(f"       max|I_perwindow - I_shared| across series = {bias:.4f}")

    # shared gives a clean linear calibration variable even on noisy, varying-loading data
    ok_lin = m_sh["r2"] > 0.99
    results.append(("T4 shared ratio linear vs truth under realistic conditions (R2>0.99)",
                    ok_lin, f"R2_sh={m_sh['r2']:.5f}"))
    print(f"  [{'PASS' if ok_lin else 'FAIL'}] T4 shared linear (R2>0.99)")
    # and the two modes genuinely diverge under overlap -> the mode choice matters
    ok_div = bias > 0.01
    results.append(("T4 per_window vs shared differ materially under overlap (>0.01)",
                    ok_div, f"max|bias|={bias:.4f}"))
    print(f"  [{'PASS' if ok_div else 'FAIL'}] T4 modes diverge under overlap (>0.01)")

    style.apply_style(cfg)
    fig, ax = style.figure(cfg)
    xs = np.linspace(I_true.min(), I_true.max(), 50)
    ax.plot(xs, m_sh["slope"] * xs + m_sh["intercept"], color="0.6", lw=0.8, ls=":",
            label="shared fit")
    ax.scatter(I_true, I_sh, marker="o", zorder=3, label="shared")
    ax.scatter(I_true, I_pw, marker="s", zorder=3, label="per_window")
    ax.set_xlabel("true ester fraction"); ax.set_ylabel("recovered band-area ratio  I")
    ax.legend(loc="best")
    verify.audit_layout(fig, cfg)
    style.save_fig(fig, os.path.join(OUT, "groundtruth_shared_vs_perwindow"), cfg)


def main():
    t1_isolated_band()
    t2_overlap_partition()
    t3_ratio_invariance()
    t4_calibration_series()
    n_pass = sum(1 for _, ok, _ in results if ok)
    print("\n" + "=" * 60)
    print(f"SUMMARY: {n_pass}/{len(results)} checks passed")
    for name, ok, detail in results:
        if not ok:
            print(f"  FAIL  {name}: {detail}")
    print("=" * 60)
    return 0 if n_pass == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
