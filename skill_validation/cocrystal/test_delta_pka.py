"""ΔpKa salt/cocrystal rule: closed-form arithmetic + zone logic + boundaries + regimes.

What's validated is the FORMULA and the Cruz-Cabeza ZONE thresholds — fully closed-form. The
pKa VALUES themselves are user look-ups (not computed or stored by the skill), so no external
fact is baked in here; the "regime" checks use qualitatively-unambiguous inputs, not stored
literature pKa.
"""
from __future__ import annotations
from scripts import cocrystal as cc
from _harness import close, prop

A = "delta_pka"
P = "Cruz-Cabeza 2012 pKa rule"


def collect():
    R = []

    # ---------- arithmetic (T0): ΔpKa = pKa(BH+) − pKa(AH) ----------
    R.append(close(A, "ΔpKa = 5.5 − 2.0 = 3.5", cc.delta_pka(2.0, 5.5), 3.5, "T0", "definition"))
    R.append(close(A, "ΔpKa = 3.0 − 4.0 = −1.0", cc.delta_pka(4.0, 3.0), -1.0, "T0", "definition"))
    R.append(close(A, "ΔpKa = 3.5 − 3.5 = 0.0", cc.delta_pka(3.5, 3.5), 0.0, "T0", "definition"))

    # ---------- zones (T0 logic): >3 salt, <0 cocrystal, else continuum ----------
    R.append(prop(A, "ΔpKa 3.5 → salt", cc.classify_ionisation(2.0, 5.5)["zone"] == "salt", "T0", P))
    R.append(prop(A, "ΔpKa −1.0 → cocrystal", cc.classify_ionisation(4.0, 3.0)["zone"] == "cocrystal", "T0", P))
    R.append(prop(A, "ΔpKa 1.5 → continuum", cc.classify_ionisation(3.0, 4.5)["zone"] == "continuum", "T0", P))

    # ---------- boundaries: continuum is INCLUSIVE (strict > and <) ----------
    R.append(prop(A, "ΔpKa == 3.0 exactly → continuum (not salt)",
                  cc.classify_ionisation(2.0, 5.0)["zone"] == "continuum", "T0", "boundary"))
    R.append(prop(A, "ΔpKa == 0.0 exactly → continuum (not cocrystal)",
                  cc.classify_ionisation(3.0, 3.0)["zone"] == "continuum", "T0", "boundary"))

    # ---------- confident flag + configurable thresholds ----------
    R.append(prop(A, "salt & cocrystal confident, continuum not",
                  cc.classify_ionisation(2.0, 5.5)["confident"]
                  and cc.classify_ionisation(4.0, 3.0)["confident"]
                  and not cc.classify_ionisation(3.0, 4.5)["confident"], "T0", "flag"))
    R.append(prop(A, "custom salt_threshold=4 moves ΔpKa 3.5 salt→continuum",
                  cc.classify_ionisation(2.0, 5.5, salt_threshold=4.0)["zone"] == "continuum", "T0", "config"))

    # ---------- illustrative REGIME checks (unambiguous inputs, not stored literature pKa) ----------
    R.append(prop(A, "strong base (pKaH 9) + mineral acid (pKa −6) → salt",
                  cc.classify_ionisation(-6.0, 9.0)["zone"] == "salt", "T3", "regime: amine·HCl"))
    R.append(prop(A, "weak base (pKaH 0.6) + weak acid (pKa 1.2) → cocrystal",
                  cc.classify_ionisation(1.2, 0.6)["zone"] == "cocrystal", "T3", "regime: caffeine-type"))

    # ---------- provenance carried in the caption ----------
    R.append(prop(A, "caption cites Cruz-Cabeza + shows the rule",
                  "Cruz-Cabeza" in cc.classify_ionisation(2.0, 5.5)["caption"], "T3", "provenance"))

    return R
