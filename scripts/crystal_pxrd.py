"""
crystal_pxrd.py  -  calculated PXRD from a CIF, handed to the spectra module (pxrd domain).

Calculated PXRD is a new *source*, not a new plotting family: this computes (2theta,
intensity) with Dans_Diffraction and plots it through the house pxrd conventions. The two
verified gotchas are handled here: the (000) direct beam is kept OFF-GRID (min_twotheta via
setup_scatter), and peak_width is in Q (A^-1), NOT degrees. Cross-checked against pymatgen
(or a Bragg-law fallback). Validated in skill_validation/dans_fix.py (ELAINM: 47 peaks,
Dans vs pymatgen top peak 0.011 deg).

    calc_pattern(struct, cfg)             -> {two_theta, intensity, wavelength, peaks}
    plot(struct, cfg, experimental=None)  -> (fig, ax, pattern)   via spectra pxrd axis
"""
from __future__ import annotations
import math
from dataclasses import replace
import numpy as np
from . import verify, spectra, style, crystal_engine


def _wavelength(struct, cfg):
    if cfg.pxrd_wavelength:
        return float(cfg.pxrd_wavelength), "cfg"
    wl = struct.declared.get("wavelength")
    if wl:
        return float(wl), "CIF"
    return 1.5406, "fallback Cu-Ka1 (CIF declared none)"


def _find_peaks(tt, inten, height=2.0, distance=5):
    try:
        from scipy.signal import find_peaks
        idx, _ = find_peaks(inten, height=height, distance=distance)
        order = idx[np.argsort(inten[idx])[::-1]]
        return [(round(float(tt[i]), 3), round(float(inten[i]), 1)) for i in order]
    except Exception:
        # crude fallback: just the global max
        i = int(np.argmax(inten))
        return [(round(float(tt[i]), 3), round(float(inten[i]), 1))]


def _bragg_nearest(struct, wl, tth, hmax=5):
    """Lightweight cross-check: 2theta of the reflection (|h,k,l|<=hmax) closest to `tth`,
    from the cell's reciprocal metric. A small gap confirms the wavelength/cell are
    self-consistent (the top peak is a real reflection), without a second full engine."""
    ca, cb, cg = (math.cos(math.radians(t)) for t in (struct.al, struct.be, struct.ga))
    a, b, c = struct.a, struct.b, struct.c
    G = np.array([[a * a, a * b * cg, a * c * cb],
                  [a * b * cg, b * b, b * c * ca],
                  [a * c * cb, b * c * ca, c * c]])
    Gs = np.linalg.inv(G)
    best = None
    for h in range(-hmax, hmax + 1):
        for k in range(-hmax, hmax + 1):
            for l in range(-hmax, hmax + 1):
                if (h, k, l) == (0, 0, 0):
                    continue
                hkl = np.array([h, k, l], float)
                inv_d2 = float(hkl @ Gs @ hkl)
                if inv_d2 <= 0:
                    continue
                s = wl * math.sqrt(inv_d2) / 2.0
                if s >= 1.0:
                    continue
                t2 = 2 * math.degrees(math.asin(s))
                if t2 < 1.0:
                    continue
                gap = abs(t2 - tth)
                if best is None or gap < best:
                    best = gap
    return best


def calc_pattern(struct, cfg):
    """Dans_Diffraction powder pattern, (000) off-grid, at the declared wavelength.
    Gates the calc/cross-check agreement per cfg.pxrd_crosscheck."""
    try:
        import Dans_Diffraction as dif
    except ImportError as e:
        raise ImportError("Required package not found: Dans-Diffraction. Install with:\n"
                          "    python -m pip install Dans-Diffraction") from e
    wl, src = _wavelength(struct, cfg)
    xtl = dif.Crystal(cfg.cif_path)
    xtl.Scatter.setup_scatter(scattering_type="xray", wavelength_a=wl,
                              powder_units="twotheta",
                              min_twotheta=cfg.pxrd_two_theta_min,
                              max_twotheta=cfg.pxrd_two_theta_max,
                              powder_lorentz=cfg.pxrd_lorentz_fraction, output=False)
    # peak_width is in Q (A^-1), NOT degrees; min_twotheta keeps the (000) at 0 off-grid
    tt, inten, _ = xtl.Scatter.powder("xray", units="tth",
                                      peak_width=cfg.pxrd_peak_width,
                                      lorentz_fraction=cfg.pxrd_lorentz_fraction)
    tt, inten = np.asarray(tt, float), np.asarray(inten, float)
    if inten.size and inten.max() > 0:
        inten = inten / inten.max() * 100.0
    peaks = _find_peaks(tt, inten)
    top = peaks[0][0] if peaks else float("nan")

    f = [("INFO", f"PXRD: lambda={wl:.5f} A ({src}); {len(peaks)} peaks >2%; "
                  f"strongest 2theta={top}")]
    if len(peaks) < 3:
        f.append(("WARN", "fewer than 3 peaks — check wavelength/cell or the (000) mask"))

    mode = cfg.pxrd_crosscheck
    if mode != "none" and not math.isnan(top):
        delta = method = other = None
        if mode in ("auto", "pymatgen"):
            try:
                from pymatgen.core import Structure as PMG
                from pymatgen.analysis.diffraction.xrd import XRDCalculator
                st = PMG.from_file(cfg.cif_path)
                pat = XRDCalculator(wavelength=wl).get_pattern(
                    st, two_theta_range=(cfg.pxrd_two_theta_min, cfg.pxrd_two_theta_max))
                other = float(pat.x[int(np.argmax(pat.y))])
                delta, method = abs(top - other), "pymatgen"
            except Exception:
                pass
        if delta is None and mode in ("auto", "bragg"):
            gap = _bragg_nearest(struct, wl, top)
            if gap is not None:
                delta, method, other = gap, "Bragg-law", None
        if delta is not None:
            sev = "WARN" if delta > 0.05 else "INFO"
            tail = f" vs {other:.3f}" if other is not None else " (nearest reflection)"
            f.append((sev, f"cross-check [{method}]: top peak {top:.3f}{tail}, "
                           f"delta {delta:.3f} deg (gate <0.05)"))
    verify._resolve(f, cfg, "pxrd")
    return {"two_theta": tt, "intensity": inten, "wavelength": wl, "peaks": peaks}


def reflection_list(struct, cfg):
    """The powder reflection list `[(two_theta_deg, intensity, (h,k,l)), ...]` from
    Dans_Diffraction — REAL structure-factor intensities with multiplicity + Lorentz-polarization
    already applied (Dans's `powder()` 3rd return: columns h,k,l,2theta,intensity, grouped by
    min_overlap). This is the correct input to `pxrd_realism.simulate_pattern` (feed it THESE, not
    the post-broadened peak list — that would double-broaden). Filtered to the cfg 2theta window."""
    try:
        import Dans_Diffraction as dif
    except ImportError as e:
        raise ImportError("Required package not found: Dans-Diffraction. Install with:\n"
                          "    python -m pip install Dans-Diffraction") from e
    wl, _ = _wavelength(struct, cfg)
    xtl = dif.Crystal(cfg.cif_path)
    xtl.Scatter.setup_scatter(scattering_type="xray", wavelength_a=wl, powder_units="twotheta",
                              min_twotheta=cfg.pxrd_two_theta_min, max_twotheta=cfg.pxrd_two_theta_max,
                              powder_lorentz=cfg.pxrd_lorentz_fraction, output=False)
    _, _, refl = xtl.Scatter.powder("xray", units="tth", peak_width=cfg.pxrd_peak_width,
                                    lorentz_fraction=cfg.pxrd_lorentz_fraction)
    refl = np.asarray(refl, float)
    lo, hi = cfg.pxrd_two_theta_min, cfg.pxrd_two_theta_max
    out = []
    for h, k, l, tth, I in refl:
        if lo <= tth <= hi and I > 0:
            out.append((float(tth), float(I), (int(round(h)), int(round(k)), int(round(l)))))
    return out


def realistic_pattern(struct, cfg, x_grid=None, npoints=3000):
    """Turn-key realistic Cu-Ka powder pattern from a CIF: real Dans reflections
    (`reflection_list`) fed to `pxrd_realism.simulate_pattern`, so the Kalpha2 doublet
    (`cfg.pxrd_kalpha2`), preferred orientation (`cfg.pxrd_po_axis`/`pxrd_march_r`) and Caglioti
    broadening (`cfg.pxrd_caglioti`) are applied ON TOP of physically-correct intensities. The
    alpha1 line is the CIF/cfg wavelength (not assumed Cu). Returns
    {two_theta, intensity, reflections, wavelength}. Use for overlaying calc vs a measured Cu-Ka
    scan; for positions-only phase ID, `calc_pattern` is enough.

    Coherence gotcha: the Kalpha2 default is **Cu** Ka2. To simulate a Cu lab scan from a
    Mo-refined CIF, set `cfg.pxrd_wavelength=1.540598` so alpha1 is Cu too (else you pair a Cu
    Ka2 onto a Mo alpha1); for a non-Cu anode set `cfg.pxrd_wavelength2` to its Ka2."""
    from . import pxrd_realism as pr
    wl, _ = _wavelength(struct, cfg)
    refl = reflection_list(struct, cfg)
    if x_grid is None:
        x_grid = np.linspace(cfg.pxrd_two_theta_min, cfg.pxrd_two_theta_max, npoints)
    Gs = pr.reciprocal_metric(struct.a, struct.b, struct.c, struct.al, struct.be, struct.ga)
    lam2 = cfg.pxrd_wavelength2 or pr.CU_KA2
    y = pr.simulate_pattern(refl, x_grid, cfg, lam1=wl, lam2=lam2, Gs=Gs)
    return {"two_theta": np.asarray(x_grid, float), "intensity": y, "reflections": refl,
            "wavelength": wl}


def plot(struct, cfg, experimental=None, pattern=None):
    """Plot the calculated pattern through the house pxrd conventions (normal 2theta axis,
    intensity). `experimental`=(2theta, I) overlays a measured trace for phase ID; that trace is
    an ingest and is gated with verify.check_trace (pxrd domain: structural checks, no FTIR
    spike heuristic). The calculated pattern is our own validated computation, not an ingest."""
    if pattern is None:
        pattern = calc_pattern(struct, cfg)
    pcfg = replace(cfg, domain="pxrd")
    style.apply_style(pcfg)
    fig, ax = style.figure(pcfg)
    ax.plot(pattern["two_theta"], pattern["intensity"], label="calculated")
    if experimental is not None:
        ex, ey = np.asarray(experimental[0], float), np.asarray(experimental[1], float)
        verify.check_trace(ex, ey, pcfg, name="experimental")
        if ey.size and ey.max() > 0:
            ey = ey / ey.max() * 100.0
        ax.plot(ex, ey, label="experimental")
    spectra._apply_axis(ax, pcfg)
    if experimental is not None:
        ax.legend(loc="best")
    style.finalize_figure(fig)
    return fig, ax, pattern


def _all_reflections(struct, wl, tth_min, tth_max, hmax=6):
    """(2theta, (h,k,l), d) for all |index|<=hmax reflections in the 2theta window, from the cell
    metric. For peak INDEXING only — no structure factors / systematic absences — an aid, not a
    full reflection list."""
    ca, cb, cg = (math.cos(math.radians(t)) for t in (struct.al, struct.be, struct.ga))
    a, b, c = struct.a, struct.b, struct.c
    G = np.array([[a * a, a * b * cg, a * c * cb],
                  [a * b * cg, b * b, b * c * ca],
                  [a * c * cb, b * c * ca, c * c]])
    Gs = np.linalg.inv(G)
    out = []
    for h in range(-hmax, hmax + 1):
        for k in range(-hmax, hmax + 1):
            for l in range(-hmax, hmax + 1):
                if (h, k, l) == (0, 0, 0):
                    continue
                hkl = np.array([h, k, l], float)
                inv_d2 = float(hkl @ Gs @ hkl)
                if inv_d2 <= 0:
                    continue
                d = 1.0 / math.sqrt(inv_d2)
                s = wl / (2 * d)
                if s >= 1.0:
                    continue
                t2 = 2 * math.degrees(math.asin(s))
                if tth_min <= t2 <= tth_max:
                    out.append((t2, (h, k, l), d))
    return out


def peak_table(struct, cfg, pattern=None, top=25):
    """(2theta, d, hkl, I_rel) for the strongest calculated peaks. hkl is the lowest-index
    cell-metric reflection nearest each peak (an indexing aid — no absences/structure factors).
    Returns a list of dicts, strongest first."""
    if pattern is None:
        pattern = calc_pattern(struct, cfg)
    wl = pattern["wavelength"]
    refl = _all_reflections(struct, wl, cfg.pxrd_two_theta_min, cfg.pxrd_two_theta_max)
    rows = []
    for tth, I in pattern["peaks"][:top]:
        th = math.radians(tth / 2.0)
        d = wl / (2 * math.sin(th)) if th > 0 else float("nan")
        near = [r for r in refl if abs(r[0] - tth) < 0.10]
        if near:
            h, k, l = min(near, key=lambda r: (sum(abs(x) for x in r[1]), abs(r[0] - tth)))[1]
            hkl = f"{h} {k} {l}"
        else:
            hkl = "?"
        rows.append({"two_theta": round(tth, 3), "d": round(d, 4), "hkl": hkl, "I": round(I, 1)})
    return rows


def write_peaks_csv(rows, path):
    """Write a peak list (2theta, d, hkl, I) to CSV — the SI companion to a PXRD phase-ID figure."""
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["2theta/deg", "d/A", "hkl", "I_rel"])
        for r in rows:
            w.writerow([r["two_theta"], r["d"], r["hkl"], r["I"]])
    return path


def plot_overlay(entries, cfg, experimental=None, offset=None):
    """Stacked calculated PXRD of several phases for phase ID (the cocrystal-vs-starting-materials
    figure), in the house waterfall idiom: **solid palette lines** (vertical position separates
    the traces) with the **right-margin per-trace key** (`spectra.edge_labels`), never labels over
    the data. `entries` = list of (label, cif_path); each pattern is computed at the house pxrd
    conventions and normalized to 100. `experimental` = (2theta, I) is drawn (grey) at the base.
    Returns (fig, ax, patterns) — patterns is a list of (label, pattern dict)."""
    pcfg = replace(cfg, domain="pxrd")
    style.apply_style(pcfg)
    fig, ax = style.figure(pcfg)
    pal = style._palette(pcfg)
    traces, patterns = [], []                       # (label, x, y, colour), bottom -> top
    if experimental is not None:
        ex, ey = np.asarray(experimental[0], float), np.asarray(experimental[1], float)
        if ey.size and ey.max() > 0:
            ey = ey / ey.max() * 100.0
        traces.append(("experimental", ex, ey, "0.3"))
    for i, (label, path) in enumerate(entries):
        ecfg = replace(pcfg, cif_path=path)
        s = crystal_engine.load(ecfg)
        pat = calc_pattern(s, ecfg)
        patterns.append((label, pat))
        traces.append((label, np.asarray(pat["two_theta"], float),
                       np.asarray(pat["intensity"], float), pal[i % len(pal)]))
    spans = [np.ptp(y) for _, _, y, _ in traces if len(y)]
    step = offset if offset is not None else (1.15 * max(spans) if spans else 100.0)
    edges = []                                      # (y_at_right_edge, label, colour) per trace
    for gi, (label, x, y, color) in enumerate(traces):
        ax.plot(x, y + gi * step, color=color, ls="-")   # solid always (house waterfall rule)
        if len(x):
            edges.append((y[int(np.argmax(x))] + gi * step, label, color))   # pxrd right edge = max 2theta
    spectra._apply_axis(ax, pcfg)
    ax.set_yticks([])                               # offsets are arbitrary; hide the y scale
    ax.set_ylabel("Intensity (normalized, offset)")
    spectra.edge_labels(ax, edges)                  # the house key: right margin, per trace
    return fig, ax, patterns
