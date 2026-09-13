"""simulate_pattern per-peak windows: `window_fwhm=None` is the exact full-grid sum; the default
window (±40 FWHM) agrees with it to within the far-tail bound; a descending grid gives the
reversed result of the ascending one. Reflection lists are the ones the pxrd validation suite
uses (skill_validation/pxrd/test_realism.py), no new data."""
import numpy as np
import pytest

from scripts import pxrd_realism as pr

REFL_POS = [(20.0, 100.0, (1, 0, 0)), (40.0, 50.0, (2, 0, 0)), (70.0, 40.0, (3, 0, 0))]
REFL_DIR = [(20.0, 100.0, (1, 0, 0)), (70.0, 100.0, (0, 0, 1))]
X = np.linspace(10.0, 80.0, 35001)
UVW = dict(U=0.03, V=-0.005, W=0.008)


def _full_sum(refl, x, eta=0.5, **uvw):
    y = np.zeros_like(x)
    for tth, I, _ in refl:
        y = y + I * pr.pseudo_voigt(x, tth, float(pr.caglioti_fwhm(tth, uvw["U"], uvw["V"], uvw["W"])), eta)
    return y


def test_window_none_is_the_exact_sum():
    y = pr.simulate_pattern(REFL_POS, X, eta=0.5, kalpha2=False, normalize=False, window_fwhm=None, **UVW)
    assert np.array_equal(y, _full_sum(REFL_POS, X, **UVW))


@pytest.mark.parametrize("refl", [REFL_POS, REFL_DIR])
def test_default_window_within_tail_bound(refl):
    yf = pr.simulate_pattern(refl, X, eta=0.5, kalpha2=False, window_fwhm=None, **UVW)
    yw = pr.simulate_pattern(refl, X, eta=0.5, kalpha2=False, **UVW)          # default window
    d = np.abs(yf - yw)
    assert d.max() < 1e-2                          # measured 6.5e-3 on a 0-100 scale at 40 FWHM
    # inside its own window each peak is exact: the difference at a peak top is only the other
    # peaks' far tails, bounded by their truncated area fraction eta*(1-2/pi*atan(80)) ~ 0.4 %
    tail = 0.5 * (1 - 2 / np.pi * np.arctan(2 * pr.DEFAULT_WINDOW_FWHM))
    top = np.argmax(yf)
    assert d[top] <= tail * yf.max()
    # the window never removes a peak: the same maxima survive
    assert np.argmax(yw) == top


def test_smaller_window_deviates_more():
    yf = pr.simulate_pattern(REFL_POS, X, eta=0.5, kalpha2=False, window_fwhm=None, **UVW)
    devs = [np.abs(yf - pr.simulate_pattern(REFL_POS, X, eta=0.5, kalpha2=False, window_fwhm=w, **UVW)).max()
            for w in (10, 20, 40)]
    assert devs[0] > devs[1] > devs[2]


def test_descending_grid_matches_reversed():
    ya = pr.simulate_pattern(REFL_POS, X, eta=0.5, kalpha2=False, **UVW)
    yd = pr.simulate_pattern(REFL_POS, X[::-1], eta=0.5, kalpha2=False, **UVW)
    assert np.array_equal(yd[::-1], ya)


def test_pick_fails_loudly_on_misspelt_field():
    from scripts import config

    class Cfg:  # a cfg object lacking the field: no silent default any more
        pass
    with pytest.raises(AttributeError):
        pr._pick(Cfg(), "pxrd_march_r", None)
    with pytest.raises(AttributeError):
        pr._pick(None, "pxrd_marchr", None)                 # misspelt: no literal to fall back on
    assert pr._pick(None, "pxrd_march_r", None) == config.Config().pxrd_march_r   # None -> the Config default
    assert pr._pick(None, "pxrd_march_r", 0.7) == 0.7
