"""Every module imports cleanly (catches drift as the skill grows). Heavy-dep crystal modules
skip if their optional deps are absent (they should lazy-import, so this rarely triggers)."""
import importlib
import pytest

CORE = ["config", "style", "verify", "spectra", "calibration", "chemometrics",
        "charts", "report", "doe", "cocrystal", "pxrd_realism", "sources"]
CRYSTAL = ["crystal_engine", "crystal_pxrd", "crystal_view"]


@pytest.mark.parametrize("mod", CORE)
def test_core_module_imports(mod):
    importlib.import_module(f"scripts.{mod}")


@pytest.mark.parametrize("mod", CRYSTAL)
def test_crystal_module_imports(mod):
    try:
        importlib.import_module(f"scripts.{mod}")
    except ImportError as e:
        pytest.skip(f"optional crystal dep missing: {e}")
