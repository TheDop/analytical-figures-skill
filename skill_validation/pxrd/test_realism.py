"""PXRD realism: March-Dollase (closed-form), reflection angles (cubic), Cu Kα2 doublet (Bragg
self-consistency + intensity ratio + splitting grows), Caglioti FWHM (formula + monotonic),
pseudo-Voigt (area + half-max for the Gaussian/Lorentzian limits), and simulate_pattern
composition (positions, normalisation, angle-dependent width, PO intensity)."""
from __future__ import annotations
import numpy as np
from scipy.signal import find_peaks
from scripts import pxrd_realism as pr
from _harness import close, prop


def _fwhm(x, y):
    y = np.asarray(y, float); x = np.asarray(x, float)
    half = y.max() / 2.0
    above = np.where(y >= half)[0]
    return float(x[above[-1]] - x[above[0]]) if above.size >= 2 else 0.0


def collect():
    R = []

    # ---------------- March-Dollase factor (closed-form) ----------------
    M = "march_dollase"
    R.append(prop(M, "r=1 → P=1 for all α",
                  bool(np.allclose(pr.march_dollase_factor([0, 0.5, 1.2, np.pi / 2], 1.0), 1.0)), "T0", "no PO"))
    R.append(close(M, "α=0, r=0.5 → r^-3 = 8", pr.march_dollase_factor(0.0, 0.5), 8.0, "T0", "(r^2)^-1.5"))
    R.append(close(M, "α=π/2, r=0.5 → r^1.5", pr.march_dollase_factor(np.pi / 2, 0.5), 0.5 ** 1.5, "T0", "(1/r)^-1.5"))
    R.append(close(M, "α=0, r=2 → 2^-3 = 0.125", pr.march_dollase_factor(0.0, 2.0), 0.125, "T0", "(r^2)^-1.5"))
    R.append(prop(M, "platy (r<1) enhances on-axis vs perpendicular",
                  pr.march_dollase_factor(0.0, 0.5) > pr.march_dollase_factor(np.pi / 2, 0.5), "T3", "habit"))

    # ---------------- reflection angles (cubic a=5) ----------------
    A = "reflection_angle"
    Gs = pr.reciprocal_metric(5, 5, 5, 90, 90, 90)
    R.append(close(A, "cubic ∠[(100),(010)] = 90°", np.degrees(pr.reflection_angle((1, 0, 0), (0, 1, 0), Gs)), 90.0, "T0", "orthogonal"))
    R.append(close(A, "cubic ∠[(100),(100)] = 0°", np.degrees(pr.reflection_angle((1, 0, 0), (1, 0, 0), Gs)), 0.0, "T0", "self"))
    R.append(close(A, "cubic ∠[(110),(100)] = 45°", np.degrees(pr.reflection_angle((1, 1, 0), (1, 0, 0), Gs)), 45.0, "T0", "diagonal"))

    # ---------------- Cu Kα2 doublet ----------------
    K = "kalpha2"
    refl = [(20.0, 100.0, (1, 0, 0)), (50.0, 60.0, (2, 0, 0))]
    dbl = pr.kalpha2_doublet(refl)
    R.append(prop(K, "doubles the reflection count", len(dbl) == 2 * len(refl), "T0", "α1+α2"))
    want20 = 2 * np.degrees(np.arcsin(np.sin(np.radians(10)) * pr.CU_KA2 / pr.CU_KA1))
    R.append(close(K, "α2(2θ=20°) matches Bragg (sinθ2=sinθ1·λ2/λ1)", dbl[1][0], want20, "T0", "Bragg λ2"))
    R.append(close(K, "α2 intensity = ratio·I", dbl[1][1], 100.0 * pr.CU_KA2_RATIO, "T0", "0.5·I"))
    R.append(prop(K, "splitting grows with angle (Δ2θ ∝ tanθ)",
                  (dbl[3][0] - 50.0) > (dbl[1][0] - 20.0), "T3", "tanθ"))

    # ---------------- Caglioti FWHM ----------------
    C = "caglioti"
    U, V, W = 0.01, -0.005, 0.008
    for tth in (30.0, 70.0):
        t = np.tan(np.radians(tth / 2))
        R.append(close(C, f"FWHM(2θ={tth:.0f}) = sqrt(U t²+V t+W)", pr.caglioti_fwhm(tth, U, V, W),
                       float(np.sqrt(U * t * t + V * t + W)), "T0", "formula"))
    R.append(prop(C, "FWHM grows from low to high angle",
                  pr.caglioti_fwhm(70, U, V, W) > pr.caglioti_fwhm(20, U, V, W), "T3", "instrumental"))

    # ---------------- pseudo-Voigt ----------------
    P = "pseudo_voigt"
    xx = np.linspace(-100.0, 100.0, 200001)
    R.append(close(P, "Gaussian (eta=0) area ≈ 1", float(np.trapezoid(pr.pseudo_voigt(xx, 0, 1.0, 0.0), xx)), 1.0, "T2", "∫=1", rtol=1e-3))
    R.append(close(P, "Lorentzian (eta=1) area ≈ 1", float(np.trapezoid(pr.pseudo_voigt(xx, 0, 1.0, 1.0), xx)), 1.0, "T2", "∫=1", rtol=2e-2))
    for eta, nm, tier, rt in ((0.0, "Gaussian", "T1", 1e-3), (1.0, "Lorentzian", "T0", 1e-9)):
        ratio = pr.pseudo_voigt(0.5, 0, 1.0, eta) / pr.pseudo_voigt(0.0, 0, 1.0, eta)
        R.append(close(P, f"{nm} half-max at ±FWHM/2", float(ratio), 0.5, tier, "FWHM def", rtol=rt))

    # ---------------- simulate_pattern composition ----------------
    S = "simulate_pattern"
    x = np.linspace(10.0, 80.0, 35001)                 # 0.002° grid — resolves the FWHM change
    Us, Vs, Ws = 0.03, -0.005, 0.008                   # bigger U → clearly-resolved high-angle broadening
    refl_pos = [(20.0, 100.0, (1, 0, 0)), (40.0, 50.0, (2, 0, 0)), (70.0, 40.0, (3, 0, 0))]
    y = pr.simulate_pattern(refl_pos, x, U=Us, V=Vs, W=Ws, eta=0.5, kalpha2=False)
    R.append(close(S, "normalised to 100", float(y.max()), 100.0, "T0", "norm"))
    tops = x[find_peaks(y, height=5)[0]]
    R.append(prop(S, "peaks at the reflection 2θ (20/40/70)",
                  all(np.any(np.abs(tops - t) < 0.3) for t in (20, 40, 70)), "T3", "positions"))
    fw_lo = _fwhm(x[(x > 18) & (x < 22)], y[(x > 18) & (x < 22)])
    fw_hi = _fwhm(x[(x > 68) & (x < 72)], y[(x > 68) & (x < 72)])
    R.append(prop(S, "high-angle peak broader than low (Caglioti)", fw_hi > fw_lo, "T3", "FWHM(θ)"))
    # March-Dollase: on-axis (100) enhanced vs perpendicular (001) when PO=[100], r<1
    refl_dir = [(20.0, 100.0, (1, 0, 0)), (70.0, 100.0, (0, 0, 1))]
    y0 = pr.simulate_pattern(refl_dir, x, kalpha2=False, normalize=False)
    ypo = pr.simulate_pattern(refl_dir, x, kalpha2=False, po_hkl=(1, 0, 0), march_r=0.4, Gs=Gs, normalize=False)
    def near(yy, c): m = (x > c - 1) & (x < c + 1); return float(yy[m].max())
    R.append(prop(S, "March-Dollase (r<1, PO=[100]) enhances (100) vs (001)",
                  (near(ypo, 20) / near(ypo, 70)) > (near(y0, 20) / near(y0, 70)), "T3", "PO intensity"))

    return R
