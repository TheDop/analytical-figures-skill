import sys, os
SKILL = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the skill root
sys.path.insert(0, SKILL)
from scripts.config import Config
from scripts import crystal_engine as ce, crystal_pxrd as cp, verify

HERE = os.path.dirname(os.path.abspath(__file__))
CIF_DIR = os.path.join(HERE, "cifs")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)


def _available(*names):
    """The named CIFs that are present in cifs/ (CCDC files are optional extras the user fetches; see README)."""
    return [n for n in names if os.path.exists(os.path.join(CIF_DIR, n))]

for cif in _available("monoclinic_aspirin__COD7247819.cif", "specialpos_urea__COD1008785.cif", "783049.cif"):
    cfg = Config(cif_path=os.path.join(CIF_DIR, cif), strict=False)
    s = ce.load(cfg)
    fig, ax, pat = cp.plot(s, cfg)
    out = os.path.join(OUT, "_pxrd_" + cif.split("__")[0].split(".")[0] + ".png")
    verify.render_preview(fig, out)
    print(f"{cif}: lambda={pat['wavelength']}  top peaks={pat['peaks'][:5]}  -> {os.path.basename(out)}")
