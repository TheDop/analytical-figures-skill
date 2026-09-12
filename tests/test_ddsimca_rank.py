"""DD-SIMCA with an empty residual space (components >= rank): finite limits, training rows accepted."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts import chemometrics as cm  # noqa: E402


def test_three_spectra_two_components_do_not_reject_themselves():
    pytest.importorskip("scipy")
    X = np.random.default_rng(1).normal(size=(3, 50))
    m = cm.ddsimca_fit(X, n_components=2)
    r = cm.ddsimca_predict(m, X)
    acc = np.asarray(r["accept"] if isinstance(r, dict) else r[0])
    assert acc.shape == (3,) and bool(acc.all())
