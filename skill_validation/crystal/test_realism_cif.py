"""Validate the Dans reflection-list adapter + realistic_pattern against REAL CIFs.

reflection_list (real structure factors) → pxrd_realism.simulate_pattern must reproduce the peak
POSITIONS of the plain calc_pattern (same reflections underneath), and the Kalpha2/Caglioti
realism must apply cleanly. Runs on aspirin and lactose (both COD). Exits
non-zero on any check failure; also renders a calc-vs-realistic overlay to eyeball.
"""
import sys, os
import numpy as np
SKILL = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the skill root
sys.path.insert(0, SKILL)
from dataclasses import replace
from scripts.config import Config
from scripts import crystal_engine as ce, crystal_pxrd as cp, style, spectra, verify
from scipy.signal import find_peaks

HERE = os.path.dirname(os.path.abspath(__file__))
CIF_DIR = os.path.join(HERE, "cifs")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
CIFS = ["monoclinic_aspirin__COD7247819.cif", "lactose__COD2206486.cif"]
CU_KA1 = 1.540598                       # simulate a Cu lab scan (both CIFs declare Mo); the
                                        # Kalpha2 default is Cu, so pin alpha1 to Cu for coherence

fails = 0


def check(name, ok, detail=""):
    global fails
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{('  ' + detail) if detail else ''}")
    if not ok:
        fails += 1


for cif in CIFS:
    print(f"\n== {cif} ==")
    cfg = Config(cif_path=os.path.join(CIF_DIR, cif), strict=False,
                 pxrd_two_theta_min=5, pxrd_two_theta_max=45, pxrd_wavelength=CU_KA1)
    s = ce.load(cfg)
    calc = cp.calc_pattern(s, cfg)
    calc_top = calc["peaks"][0][0]

    refl = cp.reflection_list(s, cfg)
    check("reflection_list non-empty", len(refl) > 5, f"{len(refl)} reflections")
    r_top = max(refl, key=lambda r: r[1])                       # strongest reflection
    check("strongest reflection ~ calc top peak", abs(r_top[0] - calc_top) < 0.20,
          f"refl {r_top[0]:.3f} (hkl {r_top[2]}) vs calc {calc_top:.3f}")

    rp = cp.realistic_pattern(s, cfg)
    tops = rp["two_theta"][find_peaks(rp["intensity"], height=5)[0]]
    rp_top = tops[np.argmin(np.abs(tops - calc_top))] if len(tops) else float("nan")
    check("realistic_pattern top peak ~ calc top peak", abs(rp_top - calc_top) < 0.20,
          f"realistic {rp_top:.3f} vs calc {calc_top:.3f}")
    check("realistic_pattern normalised to 100", abs(float(rp["intensity"].max()) - 100.0) < 1e-6)

    # Kalpha2 doublet on: a companion appears on the HIGH-2theta side of the strongest peak
    cfg2 = replace(cfg, pxrd_kalpha2=True)
    rp2 = cp.realistic_pattern(s, cfg2)
    exp_a2 = 2 * np.degrees(np.arcsin(np.sin(np.radians(r_top[0] / 2)) * 1.544426 / calc["wavelength"]))
    win = (rp2["two_theta"] > r_top[0] + 0.01) & (rp2["two_theta"] < exp_a2 + 0.10)
    check("Kalpha2 adds intensity between alpha1 and expected alpha2",
          bool(win.any() and rp2["intensity"][win].max() > 1.0),
          f"expected a2 at {exp_a2:.3f}")

# render aspirin calc vs realistic (Kalpha2 + Caglioti) to eyeball
cfg = Config(cif_path=os.path.join(CIF_DIR, CIFS[0]), strict=False,
             pxrd_two_theta_min=5, pxrd_two_theta_max=45, pxrd_kalpha2=True, pxrd_wavelength=CU_KA1, domain="pxrd")
s = ce.load(cfg)
calc = cp.calc_pattern(s, cfg)
rp = cp.realistic_pattern(s, cfg)
style.apply_style(cfg)
fig, ax = style.figure(cfg)
ax.plot(calc["two_theta"], calc["intensity"], lw=0.8, label="calc_pattern (single lambda)")
ax.plot(rp["two_theta"], rp["intensity"], lw=0.8, label="realistic (Cu Ka1/Ka2 + Caglioti)")
spectra._apply_axis(ax, cfg)
ax.legend(loc="best", fontsize=6)
out = os.path.join(OUT, "_realism_aspirin.png")
verify.render_preview(fig, out)
print(f"\nrendered {os.path.basename(out)}")

print(f"\n{'-'*50}\n{'ALL PASS' if fails == 0 else str(fails) + ' FAILED'}")
sys.exit(1 if fails else 0)
