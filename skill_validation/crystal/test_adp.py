"""Validate the ADP tensor transform: computed U_eq must match the CIF's U_iso_or_equiv."""
import sys, os
SKILL = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # the skill root
sys.path.insert(0, SKILL)
from scripts.config import Config
from scripts import crystal_engine as ce

HERE = os.path.dirname(os.path.abspath(__file__))
CIF_DIR = os.path.join(HERE, "cifs")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)


def _available(*names):
    """The named CIFs that are present in cifs/ (CCDC files are optional extras the user fetches; see README)."""
    return [n for n in names if os.path.exists(os.path.join(CIF_DIR, n))]


bad = False
for cif in _available("monoclinic_aspirin__COD7247819.cif", "disorder-groups_flufenamic__COD4118081.cif", "783049.cif"):
    cfg = Config(cif_path=os.path.join(CIF_DIR, cif), strict=False)
    s = ce.load(cfg)
    print(f"\n### {cif}  has_adp={ce.has_adp(s)}  n_aniso={len(s.aniso)}")
    print(f"  {'atom':6} {'U_eq(calc)':>11} {'U_iso(CIF)':>11} {'Δ':>8}")
    worst = 0.0
    for at in s.atoms:
        if at["label"] in s.aniso and at["u_iso"] is not None:
            ueq = ce.u_eq(s, s.aniso[at["label"]])
            d = abs(ueq - at["u_iso"])
            worst = max(worst, d)
            print(f"  {at['label']:6} {ueq:11.4f} {at['u_iso']:11.4f} {d:8.4f}")
    print(f"  --> worst |Δ| = {worst:.4f}  ({'OK' if worst < 0.002 else 'MISMATCH — transform wrong'})")
    bad = bad or worst >= 0.002

sys.exit(1 if bad else 0)
