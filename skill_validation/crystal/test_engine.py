"""Smoke-test crystal_engine against the validation CIFs (run from anywhere)."""
import sys, os, glob
SKILL = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the skill root
sys.path.insert(0, SKILL)
from scripts.config import Config
from scripts import crystal_engine as ce

HERE = os.path.dirname(os.path.abspath(__file__))
CIF_DIR = os.path.join(HERE, "cifs")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)
for p in sorted(glob.glob(os.path.join(CIF_DIR, "*.cif"))):
    name = os.path.basename(p)
    cfg = Config(cif_path=p, strict=False)     # warn, don't raise, so we see every file
    try:
        s = ce.load(cfg)
        rep = ce.validate(s, cfg)
        print("\n" + "=" * 72)
        print(ce.table_text(rep))
    except Exception as e:
        import traceback
        print(f"\n### {name}: ERROR {type(e).__name__}: {e}")
        traceback.print_exc()
