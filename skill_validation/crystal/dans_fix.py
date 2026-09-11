"""Validate the Dans_Diffraction PXRD path against two API gotchas: the (000) reflection lands ON the
2theta grid unless min_twotheta is set, and peak_width is in Q (A^-1), not degrees. Cross-checks the
top-peak position against pymatgen when it is installed. Dans_Diffraction's mathematics is correct;
these are usage gotchas."""
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import Dans_Diffraction as dif
from scipy.signal import find_peaks

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
CIF = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "cifs", "monoclinic_aspirin__COD7247819.cif")
WL = 1.54178
TTMIN, TTMAX = 5.0, 50.0

xtl = dif.Crystal(CIF)
# setup FIRST: wavelength + restrict the grid to >=5 deg (keeps the (000) at 0 off-grid)
xtl.Scatter.setup_scatter(scattering_type="xray", wavelength_a=WL,
                          powder_units="twotheta",
                          min_twotheta=TTMIN, max_twotheta=TTMAX,
                          powder_lorentz=0.5, output=False)
# peak_width is in Q (A^-1) NOT degrees -- 0.1 was the over-broadening bug; 0.02 ~ 0.28 deg at low angle
tth, inten, refl = xtl.Scatter.powder("xray", units="tth", peak_width=0.02, lorentz_fraction=0.5)
tth, inten = np.asarray(tth, float), np.asarray(inten, float)
if inten.max() > 0:
    inten = inten / inten.max() * 100

pk, _ = find_peaks(inten, height=2.0, distance=5)
order = pk[np.argsort(inten[pk])[::-1]][:8]
print(f"Dans: {len(pk)} peaks >2% ; top-8 (2theta, rel-int):")
for i in order:
    print(f"   {tth[i]:7.3f}   {inten[i]:6.1f}")
dans_top = float(tth[order[0]]) if len(order) else float("nan")

fig, ax = plt.subplots(figsize=(6, 3))
ax.plot(tth, inten, lw=0.8)
ax.set_xlabel(r"2$\theta$ / degree"); ax.set_ylabel("Intensity (norm.)")
ax.set_title(f"{os.path.basename(CIF)} calc PXRD (Dans) Cu $\\lambda$={WL} A, (000) off-grid")
ax.set_xlim(TTMIN, TTMAX)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "probe_pxrd_fixed.png"), dpi=150)
print("saved out/probe_pxrd_fixed.png")

# pymatgen cross-check (now installed)
try:
    from pymatgen.core import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    st = Structure.from_file(CIF)
    pg = XRDCalculator(wavelength=WL).get_pattern(st, two_theta_range=(TTMIN, TTMAX))
    pg_top = float(pg.x[np.argmax(pg.y)])
    print(f"\nCROSS-CHECK top-peak 2theta:  Dans={dans_top:.3f}  pymatgen={pg_top:.3f}  "
          f"delta={abs(dans_top-pg_top):.3f}  (spec gate <0.05)")
except Exception as e:
    print("\npymatgen cross-check skipped:", repr(e)[:200])
