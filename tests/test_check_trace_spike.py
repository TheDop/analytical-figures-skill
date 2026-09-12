"""check_trace's spike rule: flags an isolated single sample in any domain, passes sharp real peaks."""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts.config import Config  # noqa: E402
from scripts import verify  # noqa: E402


def _spiky(out):
    return any("spike" in m for _, m in out)


def test_sharp_ftir_band_passes():
    x = np.arange(4000.0, 650.0, -1.864)                  # the Cary 630 grid, descending
    y = 0.02 + np.exp(-0.5 * ((x - 1748) / 1.7) ** 2)     # 4 cm-1 FWHM: ~2 samples wide
    assert not _spiky(verify.check_trace(x, y, Config(strict=False, domain="ftir"), name="band"))


def test_bragg_peak_passes_and_zinger_is_caught():
    x = np.linspace(5, 80, 3859)
    y = 5 + 100 * np.exp(-0.5 * ((x - 25) / 0.042) ** 2)  # FWHM 0.1 deg ~ 5 samples
    cfg = Config(strict=False, domain="pxrd")
    assert not _spiky(verify.check_trace(x, y, cfg, name="bragg"))
    y[2000] += 500.0                                       # a one-channel zinger
    assert _spiky(verify.check_trace(x, y, cfg, name="zinger"))


def test_dropout_is_caught():
    x = np.linspace(0, 10, 200)
    y = 1 + np.sin(x)
    y[100] = -5.0
    assert _spiky(verify.check_trace(x, y, Config(strict=False), name="dropout"))
