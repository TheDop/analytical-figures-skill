"""band_metric and integrate_bands are two callers of ONE primitive (band_area): on the same window
they must agree to machine precision, on either axis order, and the shared-envelope mode must equal
the hand-built construction."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts.config import Config  # noqa: E402
from scripts import spectra  # noqa: E402


def _trace(seed=0):
    rng = np.random.default_rng(seed)
    x = np.linspace(650, 4000, 1798)
    y = (0.02 + 1e-5 * (x - 650) + 0.8 * np.exp(-0.5 * ((x - 1748) / 8.0) ** 2)
         + 0.5 * np.exp(-0.5 * ((x - 1640) / 20.0) ** 2) + 0.004 * rng.normal(size=x.size))
    return x, y


def test_two_callers_agree_on_both_axis_orders():
    x, y = _trace()
    lo, hi = 1713.0, 1789.0
    cfg = Config(strict=False, integration_windows=[(lo, hi, "ester")])
    for xx, yy in ((x, y), (x[::-1], y[::-1])):
        a = spectra.band_metric(xx, yy, "area", anchors=(lo, hi))
        b = spectra.integrate_bands(xx, yy, cfg)[0]["area"]
        c = spectra.band_area(xx, yy, lo, hi)["area"]
        assert a > 0 and abs(a - b) < 1e-12 * a and abs(a - c) < 1e-12 * a


def test_shared_envelope_equals_hand_built_line():
    x, y = _trace(1)
    cfg = Config(strict=False, integration_baseline="shared",
                 integration_windows=[(1500.0, 1700.0, "carbox"), (1700.0, 1800.0, "ester")])
    got = {r["name"]: r["area"] for r in spectra.integrate_bands(x, y, cfg)}
    y_lo, y_hi = spectra._anchor_val(x, y, 1500.0), spectra._anchor_val(x, y, 1800.0)
    for a, b, name in cfg.integration_windows:
        line = tuple(np.interp([a, b], [1500.0, 1800.0], [y_lo, y_hi]))
        m = (x >= a) & (x <= b)
        corr = np.clip(y[m] - np.interp(x[m], [a, b], line), 0, None)
        assert abs(got[name] - np.trapezoid(corr, x[m])) < 1e-12 * got[name]


def test_signed_area_and_empty_window():
    x, y = _trace(2)
    r = spectra.band_area(x, y, 1713, 1789, clip=False)
    assert np.isfinite(r["area"]) and r["corr"].size == r["xs"].size
    assert np.isnan(spectra.band_area(x, y, 5000, 5001)["area"])
