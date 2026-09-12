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


def test_variable_named_like_a_skill_module_is_flagged():
    body = "charts = []\ncharts.append(1)\nreport = {}\nreport.update(a=1)\n"
    assert bundle._shadow_check(body, None) == ["charts", "report"]


def test_no_cross_module_top_level_collisions():
    assert bundle._cross_module_collisions() == {}


def test_bundling_works_without_utf8_mode(tmp_path):
    import subprocess
    body = tmp_path / "body.py"
    body.write_text("from scripts.config import Config\nfrom scripts import calibration\n"
                    "cfg = Config()\nprint(calibration.fit([1, 2, 3, 4], [2.0, 4.1, 5.9, 8.2], cfg)['slope'])\n",
                    encoding="utf-8")
    out = tmp_path / "standalone.py"
    r = subprocess.run([sys.executable, "-X", "utf8=0", os.path.join(ROOT, "bundle.py"), str(body), "-o", str(out)],
                       capture_output=True, text=True)
    assert r.returncode == 0 and out.exists(), r.stderr[-800:]
    r2 = subprocess.run([sys.executable, "-X", "utf8=0", str(out)], capture_output=True, text=True)
    assert r2.returncode == 0, r2.stderr[-800:]


def test_def_shadowing_still_flagged():
    body = (
        "from scripts import style\n"
        "def figure(n):\n"
        "    return n\n"
        "fig, ax = style.figure(None)\n"
    )
    assert bundle._shadow_check(body, None) == ["figure"]
