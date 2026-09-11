"""bundle._shadow_check: a body VARIABLE that shadows a skill function it calls is caught, not only a def."""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import bundle  # noqa: E402


def test_assigned_variable_shadowing_a_skill_call_is_flagged():
    body = (
        "from scripts import calibration\n"
        "def main(x, y, cfg):\n"
        "    fit = calibration.fit(x, y, cfg)\n"
        "    return fit['slope']\n"
    )
    assert bundle._shadow_check(body, None) == ["fit"]


def test_clean_body_is_not_flagged():
    body = (
        "from scripts import calibration\n"
        "def main(x, y, cfg):\n"
        "    cal = calibration.fit(x, y, cfg)\n"
        "    return cal['slope']\n"
    )
    assert bundle._shadow_check(body, None) == []


def test_def_shadowing_still_flagged():
    body = (
        "from scripts import style\n"
        "def figure(n):\n"
        "    return n\n"
        "fig, ax = style.figure(None)\n"
    )
    assert bundle._shadow_check(body, None) == ["figure"]
