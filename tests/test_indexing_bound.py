"""crystal_pxrd._all_reflections: the index bound must follow the cell, not a fixed cap.

Theophylline form II (COD 2006182, Pna2_1) has a = 24.61 A: at Cu K-alpha its (2 0 0) sits at
7.18 deg 2theta and the (h 0 0) series runs to h = 12 inside a 5-50 deg window. A fixed |h| <= 6
cap dropped every higher order along that axis, so the indexing aid could not index the very
peaks the form-II identification rests on."""
import math
import os
import sys
from types import SimpleNamespace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts import crystal_pxrd as cp  # noqa: E402

CU = 1.5406
CELL = SimpleNamespace(a=24.61, b=3.83, c=8.50, al=90.0, be=90.0, ga=90.0)   # theophylline form II


def _two_theta_h00(h):
    return 2 * math.degrees(math.asin(CU * h / (2 * CELL.a)))


def test_long_axis_high_orders_are_indexed():
    refl = cp._all_reflections(CELL, CU, 5.0, 50.0)
    hkls = {r[1] for r in refl}
    for h in (2, 8, 12):
        assert (h, 0, 0) in hkls, f"({h} 0 0) missing"
        t = next(r[0] for r in refl if r[1] == (h, 0, 0))
        assert abs(t - _two_theta_h00(h)) < 1e-6


def test_forced_cap_still_honoured():
    refl = cp._all_reflections(CELL, CU, 5.0, 50.0, hmax=6)
    assert (6, 0, 0) in {r[1] for r in refl} and (8, 0, 0) not in {r[1] for r in refl}


def test_nothing_beyond_the_window():
    refl = cp._all_reflections(CELL, CU, 5.0, 50.0)
    assert refl and all(5.0 <= r[0] <= 50.0 for r in refl)
