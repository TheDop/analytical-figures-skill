"""
pxrd_realism.py — make a calculated PXRD pattern look like a real lab Cu-Kα scan.

Three physical effects a bare stick pattern / single-wavelength profile misses, added as
NUMPY-ONLY functions that operate on a REFLECTION LIST — a list of
`(two_theta_deg, intensity, (h, k, l))` — plus the cell metric. Keeping the physics off the
Dans_Diffraction path makes it testable without heavy deps (Dans supplies the structure-factor
intensities that seed the list; see `crystal.md` for wiring):

  march_dollase   preferred-orientation INTENSITY correction (platy / needle habit)
  kalpha2_doublet Cu Kα1/Kα2 peak SPLITTING (α2 at ~half intensity; splitting grows with angle)
  caglioti_fwhm + pseudo_voigt   angle-dependent peak WIDTH and shape
  simulate_pattern  composes them onto a 2θ grid (defaults pulled from cfg); each peak is
                    evaluated only on a ±`window_fwhm` FWHM slice of the grid (searchsorted),
                    so a refinement loop calling it hundreds of times stays cheap

Nothing here computes structure factors or systematic absences — feed it real reflection
intensities. Preferred orientation is the per-reflection-orientation form; a fully rigorous PO
averages over symmetry-equivalents (pass those in the list if you have them). Validated
closed-form in skill_validation/pxrd/.
"""
from __future__ import annotations
import numpy as np

# Cu radiation (Å) and the Kα2/Kα1 integrated-intensity ratio (~0.5).
CU_KA1 = 1.540598
CU_KA2 = 1.544426
CU_KA2_RATIO = 0.5
# Per-peak evaluation half-width for simulate_pattern, in FWHM units (see its docstring).
DEFAULT_WINDOW_FWHM = 40.0


def reciprocal_metric(a, b, c, al, be, ga):
    """Reciprocal metric tensor G* (Å⁻²) from cell parameters (angles in degrees). For a
    reflection h, 1/d² = h·G*·h, and the angle between two reflections uses G* as the inner
    product."""
    ca, cb, cg = (np.cos(np.radians(t)) for t in (al, be, ga))
    G = np.array([[a * a, a * b * cg, a * c * cb],
                  [a * b * cg, b * b, b * c * ca],
                  [a * c * cb, b * c * ca, c * c]], float)
    return np.linalg.inv(G)


def reflection_angle(h1, h2, Gs):
    """Angle (radians) between two reciprocal-lattice vectors via the reciprocal metric:
    cos = (h1·G*·h2) / sqrt((h1·G*·h1)(h2·G*·h2))."""
    h1 = np.asarray(h1, float); h2 = np.asarray(h2, float)
    den = float(np.sqrt((h1 @ Gs @ h1) * (h2 @ Gs @ h2)))
    if den == 0:
        return 0.0
    return float(np.arccos(np.clip(float(h1 @ Gs @ h2) / den, -1.0, 1.0)))


def march_dollase_factor(alpha, r):
    """March-Dollase preferred-orientation factor for a reflection at angle `alpha` (rad) to the
    PO axis:  P = (r²·cos²α + sin²α/r)^(-3/2).  r=1 → 1 (no PO); r<1 platy (reflections along the
    PO axis enhanced); r>1 needle. Closed form: α=0 → r⁻³, α=π/2 → r^(3/2)."""
    a = np.asarray(alpha, float)
    return (r * r * np.cos(a) ** 2 + np.sin(a) ** 2 / r) ** (-1.5)


def apply_march_dollase(reflections, po_hkl, Gs, r):
    """Scale each reflection's intensity by its March-Dollase factor (angle between its hkl and
    `po_hkl` via G*). r=1 or po_hkl=None → intensities unchanged. Returns a NEW reflection list."""
    if po_hkl is None or r == 1.0:
        return [tuple(rf) for rf in reflections]
    out = []
    for tth, I, hkl in reflections:
        a = reflection_angle(hkl, po_hkl, Gs)
        out.append((tth, I * float(march_dollase_factor(a, r)), hkl))
    return out


def kalpha2_doublet(reflections, lam1=CU_KA1, lam2=CU_KA2, ratio=CU_KA2_RATIO):
    """Split each reflection into a Kα1 line (full I at its 2θ) plus a Kα2 line (`ratio`·I at the
    2θ for λ2 and the SAME d-spacing). θ₂ from sinθ₂ = sinθ₁·λ2/λ1, so the splitting Δ2θ grows
    with tanθ (the characteristic high-angle doublet). Returns the expanded list (≈2× entries)."""
    out = []
    for tth, I, hkl in reflections:
        out.append((tth, I, hkl))
        s2 = np.sin(np.radians(tth / 2.0)) * lam2 / lam1
        if s2 < 1.0:
            out.append((float(2 * np.degrees(np.arcsin(s2))), I * ratio, hkl))
    return out


def caglioti_fwhm(two_theta_deg, U, V, W):
    """Caglioti peak FWHM (deg) vs 2θ:  FWHM² = U·tan²θ + V·tanθ + W  (Caglioti, Paoletti &
    Ricci 1958), θ in radians, U/V/W in deg². Width grows with angle. Clamped to a positive
    floor so a bad (U,V,W) can't produce a non-finite width."""
    t = np.tan(np.radians(np.asarray(two_theta_deg, float) / 2.0))
    return np.sqrt(np.clip(U * t * t + V * t + W, 1e-6, None))


def pseudo_voigt(x, center, fwhm, eta):
    """AREA-normalized pseudo-Voigt (∫=1): eta·Lorentzian + (1−eta)·Gaussian of the SAME FWHM.
    eta=0 → Gaussian, eta=1 → Lorentzian. A reflection of intensity I contributes
    I·pseudo_voigt, so integrated area = I and the peak HEIGHT falls as the peak broadens
    (physically correct — the area is the structure-factor intensity)."""
    dx = np.asarray(x, float) - center
    hwhm = fwhm / 2.0
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    gauss = np.exp(-0.5 * (dx / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    lorentz = (hwhm / np.pi) / (dx * dx + hwhm ** 2)
    return eta * lorentz + (1 - eta) * gauss


def _pick(cfg, name, override):
    """Keyword override > declared cfg field. `cfg` is optional here (the validation suite calls
    with cfg=None) and then resolves to the skill defaults, so the ONLY copy of each default is
    the Config field — nothing here can drift from it; a misspelt name fails loudly."""
    if override is not None:
        return override
    from . import config
    return getattr(config.default_cfg(cfg), name)


def simulate_pattern(reflections, x_grid, cfg=None, *, U=None, V=None, W=None, eta=None,
                     kalpha2=None, lam1=CU_KA1, lam2=None, ratio=None,
                     po_hkl=None, march_r=None, Gs=None, normalize=True,
                     window_fwhm=DEFAULT_WINDOW_FWHM):
    """Build a realistic profile on `x_grid` (2θ deg) from a reflection list. Pipeline:
    March-Dollase (if po_hkl and r≠1 and Gs given) → Kα2 doublet (if kalpha2) → sum of
    area-normalized pseudo-Voigts with Caglioti(θ) widths. Unspecified parameters default from
    `cfg` (pxrd_caglioti → U,V,W; pxrd_lorentz_fraction → eta; pxrd_kalpha2 / pxrd_wavelength2 /
    pxrd_kalpha2_ratio; pxrd_po_axis → po_hkl; pxrd_march_r). Returns the intensity on `x_grid`
    (scaled to 100 at the max when `normalize`).

    `window_fwhm`: each reflection is evaluated only on the grid slice within ±window_fwhm·FWHM
    of its centre (np.searchsorted on the sorted grid), not on the whole grid. None = full grid
    (the exact sum). The truncated part is the far Lorentzian tail, area fraction
    η·(1 − (2/π)·arctan(2·window_fwhm)) per peak — 0.40 % at the default 40 FWHM with η=0.5
    (0.80 % for a pure Lorentzian). It is DROPPED, not renormalised: renormalising would raise
    every on-grid value of the peak by that fraction and change peak heights, whereas dropping
    leaves each peak bit-identical inside its window and only removes the smooth far-tail
    pedestal that other peaks contribute at its position — a background-like term. Measured on a
    real 328-line Cu Kα1/α2 aspirin list (3000-point 5–50° grid): max deviation from the full sum
    1.8e-4 of the pattern maximum (2.6e-3 relative on peaks >5 % of max), below the ~1e-3
    quantisation of a counted lab pattern; ±10 FWHM would already cost 2 % on weak peaks."""
    uvw = _pick(cfg, "pxrd_caglioti", None)
    U = uvw[0] if U is None else U
    V = uvw[1] if V is None else V
    W = uvw[2] if W is None else W
    eta = _pick(cfg, "pxrd_lorentz_fraction", eta)
    kalpha2 = _pick(cfg, "pxrd_kalpha2", kalpha2)
    lam2 = _pick(cfg, "pxrd_wavelength2", lam2) or CU_KA2
    ratio = _pick(cfg, "pxrd_kalpha2_ratio", ratio)
    po_hkl = _pick(cfg, "pxrd_po_axis", po_hkl)
    march_r = _pick(cfg, "pxrd_march_r", march_r)
    if lam1 == CU_KA1 and _pick(cfg, "pxrd_wavelength", None):
        lam1 = float(cfg.pxrd_wavelength)           # the α1 the reflection list was computed at
    if kalpha2 and abs(lam1 - CU_KA1) > 1e-3 and lam2 == CU_KA2:
        print(f"  [WARN] pxrd_realism: Kα2 doublet uses the Cu Kα2 line ({CU_KA2} Å) on an α1 of {lam1} Å - "
              f"set pxrd_wavelength2 for a non-Cu anode, or pin pxrd_wavelength to Cu")
    if po_hkl is not None and march_r != 1.0 and Gs is None:
        print("  [WARN] pxrd_realism: preferred orientation requested (pxrd_po_axis) but no reciprocal "
              "metric Gs was given - March-Dollase NOT applied")

    refl = list(reflections)
    if po_hkl is not None and march_r != 1.0 and Gs is not None:
        refl = apply_march_dollase(refl, po_hkl, Gs, march_r)
    if kalpha2:
        refl = kalpha2_doublet(refl, lam1, lam2, ratio)

    x = np.asarray(x_grid, float)
    y = np.zeros_like(x)
    if not refl or x.size == 0:
        return y
    w = np.inf if window_fwhm is None else float(window_fwhm)   # inf => slice (0, n) = the full sum
    order = np.argsort(x, kind="stable")          # searchsorted needs an ascending grid
    xs = x[order]
    ys = np.zeros_like(xs)
    tts = np.array([r[0] for r in refl], float)
    fws = caglioti_fwhm(tts, U, V, W)             # one vectorised call, not one per reflection
    i0s = np.searchsorted(xs, tts - w * fws)      # fw > 0 always (clamped), so inf*fw is never nan
    i1s = np.searchsorted(xs, tts + w * fws)
    for (tth, I, _hkl), fw, i0, i1 in zip(refl, fws, i0s, i1s):
        if i1 > i0:
            ys[i0:i1] += I * pseudo_voigt(xs[i0:i1], tth, float(fw), eta)
    y[order] = ys
    if normalize and y.max() > 0:
        y = y / y.max() * 100.0
    return y
