"""integrate_bands / normalize('area') must give the same, positive, result on a descending axis
(the raw .spc export runs 4000 -> 650 cm-1) as on an ascending one."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts.config import Config  # noqa: E402
from scripts import spectra  # noqa: E402


def _band():
    x = np.linspace(650, 4000, 1798)
    y = 0.01 + 0.8 * np.exp(-0.5 * ((x - 1748) / 8.0) ** 2)
    return x, y


def test_area_sign_and_value_independent_of_axis_direction():
    x, y = _band()
    cfg = Config(strict=False, integration_windows=[(1713.0, 1789.0, "ester")])
    up = spectra.integrate_bands(x, y, cfg)[0]["area"]
    down = spectra.integrate_bands(x[::-1], y[::-1], cfg)[0]["area"]
    assert up > 0 and abs(up - down) < 1e-9 * up


def test_normalize_area_is_positive_on_a_descending_axis():
    x, y = _band()
    cfg = Config(strict=False, normalize="area")
    yn = spectra.normalize(x[::-1], y[::-1], cfg)
    assert np.all(yn >= 0) and abs(abs(np.trapezoid(yn, x[::-1])) - 1.0) < 1e-9
