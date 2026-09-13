"""component_scan fits ONE model per fold at the cap and truncates it to every smaller
component count (PLS1: leading block of the NIPALS rotation; PCR: leading OLS coefficients on
orthogonal scores). These tests pin that shortcut to the from-scratch refit it replaces, and
pin permutation_test's cached per-fold preprocessing to the uncached loop, on REAL data:
the six pectin calibration standards (Wang 2023, CC BY 4.0, vendored under
skill_validation/ftir_integration/data/pectin/) plus the validation suite's own fixture."""
import glob
import os
import sys
import warnings

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # conftest puts it on sys.path
sys.path.insert(0, os.path.join(ROOT, "skill_validation", "chemometrics"))

pytest.importorskip("sklearn")
from scripts import chemometrics as cm, config  # noqa: E402
import _data  # noqa: E402  (the suite's deterministic fixture recipe)

PECTIN = os.path.join(ROOT, "skill_validation", "ftir_integration", "data", "pectin", "calibration")


def _pectin():
    files = sorted(glob.glob(os.path.join(PECTIN, "DM*.csv")))
    X = np.array([np.loadtxt(f, delimiter=",")[:, 1] for f in files])
    y = np.array([float(os.path.basename(f)[2:-4]) for f in files])
    return X, y, None


def _synth():
    _x, X, y, g = _data.synth(seed=4)
    return X, y, g


def _cfg(**kw):
    return config.Config(**kw)


def _scan_from_scratch(X, y, cfg, pre_options, method, groups):
    """The refit-per-count reference component_scan replaced, built from the same primitives
    (same folds, same per-fold preprocessing, one _fit_predict per component count)."""
    _, X, y = cm._check_xy(X, y, cfg)
    folds, _, _ = cm._make_folds(len(y), y, groups, None, cfg)
    cap = max(int(min(cfg.pls_max_components, min(len(tr) for tr, _ in folds) - 1, X.shape[1])), 1)
    curves = {}
    for pre in pre_options:
        blocks = cm._fold_blocks(X, cfg, pre, folds)
        rms = []
        for a in range(1, cap + 1):
            pred = cm._cv_predict(X, y, a, cfg, pre, folds, method, blocks=blocks)
            rms.append(float(np.sqrt(np.nanmean((y - pred) ** 2))))
        curves[cm.Preprocessor(pre, cfg).name] = rms
    return {"curves": curves, "components": list(range(1, cap + 1))}


@pytest.mark.parametrize("data", [_pectin, _synth], ids=["pectin", "synth"])
@pytest.mark.parametrize("method", ["pls", "pcr"])
@pytest.mark.parametrize("pls_scale", [True, False])
def test_nested_scan_equals_from_scratch(data, method, pls_scale):
    X, y, groups = data()
    cfg = _cfg(pls_scale=pls_scale)
    kw = dict(pre_options=("snv", "snv+center", "d1", "msc"), method=method, groups=groups)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fast = cm.component_scan(X, y, cfg, **kw)
        slow = _scan_from_scratch(X, y, cfg, **kw)
    assert fast["components"] == slow["components"]
    for name in slow["curves"]:
        a, b = np.asarray(fast["curves"][name]), np.asarray(slow["curves"][name])
        assert np.allclose(a, b, rtol=0, atol=1e-10 * max(1.0, float(b.max()))), (name, a - b)


def test_truncation_identity_per_fold():
    """The per-fold primitive itself, column by column against fresh fits."""
    X, y, _ = _pectin()
    cfg = _cfg()
    tr, te = np.arange(1, len(y)), np.array([0])
    for method in ("pls", "pcr"):
        nested = cm._fit_predict_nested(X[tr], y[tr], X[te], 4, cfg, method)
        for a in range(1, 5):
            ref = cm._fit_predict(X[tr], y[tr], X[te], a, cfg, method)
            assert np.allclose(nested[:, a - 1], ref, rtol=0, atol=1e-10 * abs(ref).max())


def _null_uncached(X, y, nc, cfg, n_perm, groups, scheme, seed):
    """The loop permutation_test used before the per-fold cache: one cross_validate per draw."""
    rng = np.random.default_rng(seed)
    out = np.empty(n_perm)
    for i in range(n_perm):
        yp = cm._permute_y(y, groups, rng)
        out[i] = cm.cross_validate(X, yp, nc, cfg, groups=groups, scheme=scheme,
                                   report=False)["r2cv"]
    return out


@pytest.mark.parametrize("data, scheme", [(_pectin, "loo"), (_pectin, "kfold"), (_synth, None)],
                         ids=["pectin-loo", "pectin-kfold", "synth-groups"])
def test_permutation_null_bit_identical(data, scheme):
    X, y, groups = data()
    cfg = _cfg()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        perm = cm.permutation_test(X, y, 2, cfg, n_perm=25, groups=groups, scheme=scheme, seed=3)
        ref = _null_uncached(X, y, 2, cfg, 25, groups, scheme, seed=3)
    assert np.array_equal(perm["null"], ref)          # bit-identical, not merely close
    assert perm["p_value"] == (int(np.sum(ref >= perm["observed"])) + 1) / 26
