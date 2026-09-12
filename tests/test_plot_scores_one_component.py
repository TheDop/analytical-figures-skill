"""diagnostics_figure must build for a 1-component model (the parsimony pick is often 1 LV)."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts.config import Config  # noqa: E402
from scripts import chemometrics as cm  # noqa: E402


def test_one_lv_diagnostics_figure_builds(tmp_path):
    pytest.importorskip("sklearn")
    rng = np.random.default_rng(0)
    x = np.linspace(1800, 1000, 60)
    levels = np.repeat(np.linspace(5, 50, 6), 2)
    X = levels[:, None] * np.exp(-0.5 * ((x - 1748) / 12) ** 2)[None, :] + 0.01 * rng.normal(size=(12, 60))
    cfg = Config(output_dir=str(tmp_path), formats=("png",), strict=False)
    fig, axes, info = cm.diagnostics_figure(X, levels, x, cfg, viz_components=1)
    assert fig is not None and len(axes) == 4
