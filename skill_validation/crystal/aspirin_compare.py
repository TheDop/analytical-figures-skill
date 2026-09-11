"""Render aspirin (296 K -> strong anisotropy) ball-and-stick vs ORTEP ellipsoids, side by
side in out/, so the ellipsoid feature is visible on a room-temperature structure."""
import sys, os
from dataclasses import replace
SKILL = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the skill root
sys.path.insert(0, SKILL)
from scripts.config import Config
from scripts import crystal_engine as ce, crystal_view as cv

HERE = os.path.dirname(os.path.abspath(__file__))
CIF_DIR = os.path.join(HERE, "cifs")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
cif = os.path.join(CIF_DIR, "monoclinic_aspirin__COD7247819.cif")

base = Config(cif_path=cif, output_dir=OUT, view_hide_ch=True,
              view_label_atoms="hetero", strict=False)
s = ce.load(base)
cv.save(cv.render(s, base), os.path.join(OUT, "aspirin_ball_stick"), base)
ell = replace(base, view_style="ellipsoid", adp_probability=0.50)
cv.save(cv.render(s, ell), os.path.join(OUT, "aspirin_ellipsoid"), ell)
print("wrote aspirin_ball_stick.png and aspirin_ellipsoid.png ->", OUT)
