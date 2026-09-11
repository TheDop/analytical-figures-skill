"""cli.py --help must render: a bare % in an argparse help string makes argparse raise on Python 3.13."""
import os
import sys
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLI = os.path.join(ROOT, "cli.py")


def _help(*args):
    r = subprocess.run([sys.executable, CLI, *args, "--help"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    return r.stdout


def test_top_level_help():
    out = _help()
    assert "spc-rsd" in out and "calibrate" in out


def test_subcommand_help():
    assert "usage" in _help("spc-rsd").lower()
    assert "usage" in _help("calibrate").lower()
