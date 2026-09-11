"""Smoke-test a HIGH-QUALITY static 3D render via PyVista (VTK), reusing the engine for
the data + orientation. Goal: gap-free, properly occluded, publication-grade still image."""
import sys, os
import numpy as np
SKILL = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the skill root
sys.path.insert(0, SKILL)
from scripts.config import Config
from scripts import crystal_engine as ce, crystal_view as cv
import pyvista as pv

HERE = os.path.dirname(os.path.abspath(__file__))
CIF_DIR = os.path.join(HERE, "cifs")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
cif = os.path.join(CIF_DIR, "monoclinic_aspirin__COD7247819.cif")
cfg = Config(cif_path=cif, view_hide_ch=True, strict=False)

s = ce.load(cif and cfg)
atoms = cv.complete_molecules(s, cfg)
if cfg.view_hide_ch:
    atoms = cv._hide_ch(atoms)
R, view, _ = cv.orient(s, cfg, atoms)
P = np.array([a["xyz"] for a in atoms])
ctr = P.mean(0)
bonds = cv._bond_pairs(atoms)
hbs = cv._hbond_pairs(atoms, cfg)


def crad(sym):
    try:
        r = ce._gemmi().Element(sym).covalent_r or 0.7
    except Exception:
        r = 0.7
    return 0.30 * (r / 0.7) + 0.08          # ball-and-stick: balls < vdW, H smaller


def rgb(hexstr):
    h = hexstr.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


pl = pv.Plotter(off_screen=True, window_size=[1800, 1800], lighting="light kit")
pl.set_background("white")

for a, p in zip(atoms, P):
    pl.add_mesh(pv.Sphere(radius=crad(a["sym"]), center=p, theta_resolution=48, phi_resolution=48),
                color=rgb(cv._color(a["sym"])), smooth_shading=True, specular=0.3, specular_power=15)

for i, j in bonds:
    a_, b_ = P[i], P[j]
    d = b_ - a_
    L = float(np.linalg.norm(d))
    ci, cj = cv._bond_color(atoms[i]["sym"]), cv._bond_color(atoms[j]["sym"])
    if ci == cj:
        pl.add_mesh(pv.Cylinder(center=(a_ + b_) / 2, direction=d, radius=0.11, height=L, resolution=28),
                    color=rgb(ci), smooth_shading=True, specular=0.3, specular_power=15)
    else:
        for c, q in ((ci, (a_ + (a_ + b_) / 2) / 2), (cj, (b_ + (a_ + b_) / 2) / 2)):
            pl.add_mesh(pv.Cylinder(center=q, direction=d, radius=0.11, height=L / 2, resolution=28),
                        color=rgb(c), smooth_shading=True, specular=0.3, specular_power=15)

for di, ai in hbs:                                   # H-bonds: thin dashed cylinders
    a_, b_ = P[di], P[ai]
    d = b_ - a_
    L = float(np.linalg.norm(d))
    n = max(3, int(L / 0.35))
    for t in range(0, n, 2):
        c0 = a_ + d * (t / n)
        c1 = a_ + d * (min(t + 1, n) / n)
        pl.add_mesh(pv.Cylinder(center=(c0 + c1) / 2, direction=d, radius=0.045,
                                height=float(np.linalg.norm(c1 - c0)), resolution=12),
                    color=(0.35, 0.35, 0.35))

# deterministic orthographic camera from the orientation engine (z = line of sight, y = up)
D = float(np.ptp(P, axis=0).max()) * 3 + 5
pl.enable_parallel_projection()
pl.camera_position = [tuple(ctr + R[2] * D), tuple(ctr), tuple(R[1])]
pl.reset_camera()
for fn in ("enable_anti_aliasing", "enable_ssao"):
    try:
        (pl.enable_anti_aliasing("ssaa") if fn == "enable_anti_aliasing" else pl.enable_ssao())
    except Exception as e:
        print(fn, "skipped:", repr(e)[:80])

out = os.path.join(OUT, "_pyvista_aspirin.png")
pl.screenshot(out)
print("saved", out)
