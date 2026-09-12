"""chemometrics._make_folds('kfold'): stratified by y means every fold spans the response range.
A contiguous split of the y-sorted order held out one y range per fold (end folds = extrapolation)."""
import os
import sys
from types import SimpleNamespace

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from scripts import chemometrics as cm  # noqa: E402


def test_every_fold_spans_the_y_range():
    y = np.arange(12, dtype=float)                       # 12 samples, y sorted
    cfg = SimpleNamespace(cv_folds=3, cv_scheme="kfold")
    folds, name, _ = cm._make_folds(len(y), y, None, "kfold", cfg)
    assert name.startswith("3-fold")
    tests = [te for _, te in folds]
    assert sorted(np.concatenate(tests).tolist()) == list(range(12))        # a partition
    for te in tests:
        assert y[te].min() < 4 and y[te].max() >= 8, f"fold {te} does not span the y range"


def test_train_and_test_are_disjoint_and_complete():
    y = np.random.default_rng(0).normal(size=17)
    cfg = SimpleNamespace(cv_folds=5, cv_scheme="kfold")
    folds, _, _ = cm._make_folds(len(y), y, None, "kfold", cfg)
    for tr, te in folds:
        assert not set(tr) & set(te)
        assert len(tr) + len(te) == 17
