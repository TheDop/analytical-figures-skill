"""cfg.view_tilt must actually move the camera on BOTH backends.

The bug this pins: orient() used to return the tilt only in its (elev, azim, roll) angles, but
_render_pyvista -- the DEFAULT renderer -- discards those angles and builds its camera from R
alone (`R, _, findings = orient(...)`). So view_tilt was a silent no-op on the default backend
while appearing to work in the matplotlib fallback. Hit for real: a 90 deg roll requested to fix a
panel's aspect ratio produced a byte-identical image.

The fix folds the tilt into R, so these tests assert on R, not on a rendered image -- no PyVista,
no VTK, and they run in milliseconds.
"""
import math

import numpy as np
import pytest

from scripts import crystal_view as CV
from scripts.config import Config


def _atoms():
    """A deliberately anisotropic little molecule: long in x, short in y, so a roll about the
    line of sight visibly swaps the projected width and height."""
    xyz = [(-4.0, 0.0, 0.0), (-2.0, 0.4, 0.0), (0.0, 0.0, 0.0),
           (2.0, -0.4, 0.0), (4.0, 0.0, 0.0), (0.0, 1.2, 0.0)]
    return [{"label": f"C{i}", "sym": "C", "xyz": np.array(p, float),
             "frac": np.zeros(3), "occ": 1.0} for i, p in enumerate(xyz)]


class _Struct:
    """orient() only reaches into the structure for the axis choice, which 'pca' does not use."""
    cell = (10.0, 10.0, 10.0, 90.0, 90.0, 90.0)
    spacegroup = "P 1"


def _R(roll=0.0, elev=0.0, azim=0.0):
    cfg = Config(strict=False, view_roll_objective="long_axis",
                 view_tilt=(elev, azim, roll))
    R, _view, _f = CV.orient(_Struct(), cfg, _atoms())
    return R


def _extent(R):
    P = np.array([a["xyz"] for a in _atoms()])
    Pr = (R @ (P - P.mean(0)).T).T
    return (float(Pr[:, 0].max() - Pr[:, 0].min()),
            float(Pr[:, 1].max() - Pr[:, 1].min()))


def test_zero_tilt_is_the_untilted_base():
    assert np.allclose(_R(), _R(roll=0.0, elev=0.0, azim=0.0))


def test_roll_changes_R():
    """The regression itself: a requested roll must reach R, not just the discarded angles."""
    assert not np.allclose(_R(), _R(roll=90.0)), "view_tilt roll did not reach the camera"


def test_roll_90_swaps_projected_width_and_height():
    w0, h0 = _extent(_R())
    w1, h1 = _extent(_R(roll=90.0))
    assert w1 == pytest.approx(h0, abs=1e-6)
    assert h1 == pytest.approx(w0, abs=1e-6)


@pytest.mark.parametrize("kw", [{"roll": 37.0}, {"elev": 25.0}, {"azim": -40.0},
                                {"roll": 15.0, "elev": 10.0, "azim": 20.0}])
def test_every_tilt_axis_is_live_and_stays_a_proper_rotation(kw):
    R = _R(**kw)
    assert not np.allclose(R, _R()), f"{kw} had no effect"
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-9), "a mirror would flip chirality"
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-9)


def test_angles_no_longer_double_count_the_tilt():
    """matplotlib applies R AND the returned angles, so the tilt must live in exactly one of
    them. It is in R, so the angles stay at the base."""
    cfg = Config(strict=False, view_roll_objective="long_axis", view_tilt=(11.0, 22.0, 33.0))
    _R_, view, _f = CV.orient(_Struct(), cfg, _atoms())
    assert view == (90.0, -90.0, 0.0)


def test_custom_orientation_still_takes_the_tilt_on_its_angles():
    """The 'custom' branch returns R=None and drives the camera from angles alone, so there the
    tilt must still be ADDED to them -- the opposite of the pca/axis path."""
    cfg = Config(strict=False, view_orientation="custom", view_angles=(10.0, 20.0, 30.0),
                 view_tilt=(1.0, 2.0, 3.0))
    R, view, _f = CV.orient(_Struct(), cfg, _atoms())
    assert R is None
    assert view == pytest.approx((11.0, 22.0, 33.0))
