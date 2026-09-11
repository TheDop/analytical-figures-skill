"""Tiny zero-dep test harness for the PXRD-realism validation suite.

Each test_*.py exposes `collect() -> list[Result]`. A Result records the tolerance TIER it was
judged at and the PROVENANCE of the reference number, so the run doubles as reproducibility
evidence. Mirror of skill_validation/chemometrics/_harness.py (kept local so the suite is
self-contained). No pytest dependency.
"""
from __future__ import annotations
import numpy as np

_TIER_TOL = {
    "T0": dict(rtol=1e-12, atol=1e-12),   # closed-form / hand-derived
    "T1": dict(rtol=1e-8,  atol=1e-10),   # same-language reference
    "T2": dict(rtol=2e-2,  atol=1e-6),    # numerical (quadrature) / published literal
    "T3": dict(rtol=0.0,   atol=0.0),     # property / invariant
}


class Result(dict):
    pass


def _mk(area, name, tier, status, got=None, want=None, tol=None, provenance="", note=""):
    return Result(area=area, name=name, tier=tier, status=status, got=got, want=want,
                  tol=tol, provenance=provenance, note=note)


def close(area, name, got, want, tier="T1", provenance="", rtol=None, atol=None, note=""):
    d = _TIER_TOL.get(tier, _TIER_TOL["T1"])
    rtol = d["rtol"] if rtol is None else rtol
    atol = d["atol"] if atol is None else atol
    g = np.asarray(got, float); w = np.asarray(want, float)
    if g.shape != w.shape:
        return _mk(area, name, tier, "FAIL", got=str(g.shape), want=str(w.shape),
                   tol="shape", provenance=provenance, note="shape mismatch")
    diff = np.abs(g - w)
    denom = np.maximum(np.abs(w), 1e-300)
    maxabs = float(np.max(diff)) if diff.size else 0.0
    maxrel = float(np.max(diff / denom)) if diff.size else 0.0
    ok = bool((maxabs <= atol) or (maxrel <= rtol))
    return _mk(area, name, tier, "PASS" if ok else "FAIL",
               got=f"abs={maxabs:.2e} rel={maxrel:.2e}",
               want=f"atol={atol:g} rtol={rtol:g}",
               tol=f"abs{maxabs:.1e}/rel{maxrel:.1e}", provenance=provenance, note=note)


def prop(area, name, ok, tier="T3", provenance="", note=""):
    return _mk(area, name, tier, "PASS" if bool(ok) else "FAIL",
               got=bool(ok), want=True, provenance=provenance, note=note)


def skip(area, name, reason, provenance=""):
    return _mk(area, name, "-", "SKIP", note=reason, provenance=provenance)
