"""Bridge the numeric validation suites into the pytest run, so `pytest` is one green gate for the
whole skill: the chemometrics (100 checks), cocrystal (26) and PXRD (23) suites, the FTIR ground-truth
integrator (14) and the doe self-check must all exit 0; the crystal set runs when gemmi + Dans_Diffraction are installed."""
import os
import sys
import subprocess

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run(path, *args):
    r = subprocess.run([sys.executable, path, *args], capture_output=True, text=True)
    return r


def test_chemometrics_validation_suite():
    r = _run(os.path.join(ROOT, "skill_validation", "chemometrics", "run.py"))
    assert r.returncode == 0, (r.stdout[-3000:] + "\n" + r.stderr[-1500:])


def test_cocrystal_validation_suite():
    r = _run(os.path.join(ROOT, "skill_validation", "cocrystal", "run.py"))
    assert r.returncode == 0, (r.stdout[-3000:] + "\n" + r.stderr[-1500:])


def test_pxrd_validation_suite():
    r = _run(os.path.join(ROOT, "skill_validation", "pxrd", "run.py"))
    assert r.returncode == 0, (r.stdout[-3000:] + "\n" + r.stderr[-1500:])


def test_ftir_groundtruth_integration():
    r = _run(os.path.join(ROOT, "skill_validation", "ftir_integration", "groundtruth_integration_test.py"))
    assert r.returncode == 0, (r.stdout[-3000:] + "\n" + r.stderr[-1500:])


def test_crystal_validation_set():
    pytest.importorskip("gemmi")
    pytest.importorskip("Dans_Diffraction")
    r = _run(os.path.join(ROOT, "skill_validation", "crystal", "run.py"))
    assert r.returncode == 0, (r.stdout[-3000:] + "\n" + r.stderr[-1500:])


def test_doe_selfcheck():
    r = _run(os.path.join(ROOT, "scripts", "doe.py"))
    assert r.returncode == 0, (r.stdout[-2000:] + "\n" + r.stderr[-1000:])
