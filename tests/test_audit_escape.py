"""audit_layout's data-escape check: a trace that leaves the panel must be caught.

The motivating failure was a waterfall normalised by subtracting the MEDIAN. On a sloping
baseline the window minimum then sits below zero, the lowest trace runs off the bottom of the
panel, and every other check in verify.py passes because no TEXT is clipped. These tests pin the
behaviour and, just as importantly, the cases that must stay silent -- a deliberate x-zoom,
axvline/axhline guides, and an explicit clip_on=False opt-out.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from scripts import verify
from scripts.config import Config

CFG = Config(domain="ftir", strict=False)


def _sloping_band():
    """A band on a baseline that slopes across the window -- the shape that triggers it."""
    x = np.linspace(650.0, 1800.0, 800)
    y = np.exp(-((x - 1700) / 25.0) ** 2) + 0.35 * np.exp(-((x - 1100) / 40.0) ** 2)
    return x, y - 0.30 * (x - 650.0) / 1150.0


def _escape_msgs(fig):
    return [m for sev, m in verify.audit_layout(fig, CFG) if "outside the y-limits" in m]


def test_median_normalised_waterfall_is_caught():
    x, y = _sloping_band()
    fig, ax = plt.subplots()
    for k in range(3):
        ax.plot(x, (y - np.median(y)) / np.abs(y).max() + k * 1.25)
    ax.set_ylim(-0.12, 3 * 1.25)
    msgs = _escape_msgs(fig)
    plt.close(fig)
    assert msgs, "a trace leaving the bottom of the panel must be reported"
    assert "MINIMUM" in msgs[0], "the message should name the fix, not just the fault"


def test_min_normalised_waterfall_is_clean():
    x, y = _sloping_band()
    fig, ax = plt.subplots()
    z = y - y.min()
    for k in range(3):
        ax.plot(x, z / z.max() + k * 1.25)
    ax.set_ylim(-0.06, 2 * 1.25 + 1.08)
    msgs = _escape_msgs(fig)
    plt.close(fig)
    assert not msgs


@pytest.mark.parametrize("build", [
    # a deliberate x-zoom clips in x by design and must not fire
    lambda ax, x, y: (ax.plot(x, y), ax.set_xlim(1650, 1780),
                      ax.set_ylim(y.min() - 0.05, y.max() + 0.05)),
    # axvline/axhline use a blended transform; their y data are axes fractions
    lambda ax, x, y: (ax.plot(x, y), ax.set_ylim(y.min() - 0.05, y.max() + 0.05),
                      ax.axvline(1700), ax.axhline(float(np.median(y)))),
    # clip_on=False is an explicit opt-out
    lambda ax, x, y: (ax.plot(x, y - 3.0, clip_on=False), ax.set_ylim(0, 1)),
])
def test_legitimate_cases_stay_silent(build):
    x, y = _sloping_band()
    fig, ax = plt.subplots()
    build(ax, x, y)
    msgs = _escape_msgs(fig)
    plt.close(fig)
    assert not msgs


def test_small_excursion_below_tolerance_is_ignored():
    """Sub-1 % overspill is invisible in print and must not become noise."""
    x = np.linspace(0.0, 10.0, 200)
    fig, ax = plt.subplots()
    ax.plot(x, np.sin(x))
    ax.set_ylim(-1.0 + 0.002, 1.0)          # clips by 0.1 % of the axis height
    msgs = _escape_msgs(fig)
    plt.close(fig)
    assert not msgs
