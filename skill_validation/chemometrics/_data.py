"""Deterministic synthetic spectra for the validation suite.

A binary aspirin/lactose-like mixture: two Gaussian pure components, a concentration
gradient over LEVELS x REPS, multiplicative loading variation (the project's real failure
mode), a broad baseline, and small noise. Fixed seed => reproducible. NOT a substitute for
the published golden datasets (gasoline etc.) — those are referenced where exactness matters.
"""
from __future__ import annotations
import numpy as np


def _gauss(x, c, w):
    return np.exp(-0.5 * ((x - c) / w) ** 2)


def synth(n_levels=6, reps=3, n_wave=80, seed=0, noise=0.002):
    rng = np.random.default_rng(seed)
    x = np.linspace(1800.0, 1000.0, n_wave)          # cm-1, descending (FTIR)
    aspirin = _gauss(x, 1748, 12) + 0.30 * _gauss(x, 1605, 15)
    lactose = _gauss(x, 1070, 30) + 0.20 * _gauss(x, 1340, 20)
    frac = np.repeat(np.linspace(0.05, 0.50, n_levels), reps)   # aspirin mass fraction
    loading = 1.0 + 0.15 * rng.standard_normal(frac.size)       # multiplicative loading swing
    X = frac[:, None] * aspirin + (1.0 - frac)[:, None] * lactose
    X = X * loading[:, None]
    X = X + 0.01 * _gauss(x, 1400, 200)[None, :]                # broad baseline
    X = X + noise * rng.standard_normal(X.shape)
    y = frac * 100.0                                            # % w/w
    groups = np.repeat(np.arange(n_levels), reps)               # concentration-level labels
    return x, X, y, groups
