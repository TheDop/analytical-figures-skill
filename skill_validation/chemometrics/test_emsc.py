"""EMSC: bit-for-bit parity against the biospectools / Kohler gold (T1), plus closed-form
coefficient-recovery and interferent-removal checks (T0). SPEC.md §6.8."""
from __future__ import annotations
import numpy as np
from pathlib import Path
from scripts import chemometrics as cm
from _harness import close, prop, skip


def collect():
    return _gold() + _synthetic()


def _gold():
    p = Path(__file__).resolve().parent / "datasets" / "emsc_testdata.xlsx"
    if not p.exists():
        return [skip("emsc", "biospectools gold", "datasets/emsc_testdata.xlsx absent (run datasets/fetch.py emsc)")]
    try:
        import openpyxl
        ws = openpyxl.load_workbook(p, data_only=True, read_only=True)["Sheet1"]
        rows = {r[0]: np.array(r[1:], float) for r in ws.iter_rows(min_row=1, values_only=True)}
    except Exception as e:
        return [skip("emsc", "biospectools gold", f"openpyxl/read failed ({type(e).__name__})")]

    wn = rows["wn"]
    X = np.vstack([rows["Raw spec 1"], rows["Raw spec 2"]])
    ref = X.mean(0)
    corr_gold = np.vstack([rows["EMSC corr quartic 1"], rows["EMSC corr quartic 2"]])
    res_gold = np.vstack([rows["EMSC res quartic 1"], rows["EMSC res quartic 2"]])
    corr, mdl = cm.emsc(X, wn, reference=ref, poly_order=4, return_model=True)
    return [
        close("emsc", "corrected vs biospectools gold (poly=4, ref=mean)", corr, corr_gold, "T1",
              "biospectools emsc_testdata.xlsx (vs Kohler MATLAB); Afseth & Kohler 2012", atol=1e-10),
        close("emsc", "residual vs biospectools gold", mdl["residual"], res_gold, "T1",
              "biospectools emsc_testdata.xlsx", atol=1e-10),
    ]


def _synthetic():
    R = []
    wn = np.linspace(1800.0, 1000.0, 200)
    g = lambda c, w: np.exp(-0.5 * ((wn - c) / w) ** 2)
    pure = g(1748, 12) + 0.5 * g(1300, 40)
    interf = g(1070, 30)

    # coefficient recovery: pure reference + poly_order=0 -> b is the true multiplicative factor
    b_true = 1.7
    corr, mdl = cm.emsc(b_true * pure, wn, reference=pure, poly_order=0, return_model=True)
    R.append(close("emsc", "b recovered (pure ref, poly=0)", mdl["coefs"][0, 0], b_true, "T0", "x=b*pure -> b"))
    R.append(close("emsc", "corrected == pure (pure ref, poly=0)", corr, pure, "T0", "(x-0)/b == pure"))

    # interferent removal: x = b*pure + d*interf -> corrected recovers pure
    x = 1.3 * pure + 0.8 * interf
    corr2 = cm.emsc(x, wn, reference=pure, poly_order=0, interferents=interf)
    R.append(close("emsc", "interferent removed -> corrected == pure", corr2, pure, "T0",
                   "subtract a known constituent (the pure-lactose use)"))

    # edge: constant spectrum -> no crash
    try:
        cm.emsc(np.ones(50), np.linspace(1800, 1000, 50), poly_order=2)
        ok = True
    except Exception:
        ok = False
    R.append(prop("emsc", "constant spectrum -> no crash", ok, "T3", "degenerate input"))
    return R
