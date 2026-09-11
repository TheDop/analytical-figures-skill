#!/usr/bin/env python3
# Self-contained, auto-bundled by analytical-figures/bundle.py.
# Depends only on the scientific stack; edit the CONFIG block to retune.
#
# --- provenance (auto-stamped by analytical-figures) ---
# skill:  analytical-figures v1.0
# env:    python 3.13.6  |  numpy 2.4.4  scipy 1.17.1  matplotlib 3.10.9  pandas 3.0.2  scikit-learn 1.9.0
# source: analysis.py
# shows:  Doxorubicin mid-IR calibration from the published carbonyl peak-area table (Bansal 2021, BMC Chem 15:27): fit, LOD/LOQ, and read-back of two unknowns
# reproduce: run this file as-is; edit the CONFIG block to retune.
# -------------------------------------------------------

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Sequence, List, Tuple
import os
import matplotlib
import matplotlib.pyplot as plt
import re
import numpy as np
import math
from dataclasses import replace
import os, sys


# ======================================================================
# from config.py
# ======================================================================
"""
py  -  THE one tunable block for the whole skill.

The discipline this enforces: every knob a user might reasonably want to turn
lives HERE, in one dataclass, with units stated per field. Nothing
reasonable-to-tune should live buried in the body of a script.

Field groups:  I/O  |  style  |  spectra  |  calibration  |  verification
Units are stated per field. Anything not here is, by definition, not meant to
be casually tuned -- if you find yourself wanting to, add it here first.
"""

#   Single source of truth for the skill version. Stamped into a bundled
#   deliverable's provenance header (bundle.py) and into a figure file's metadata
#   (save_fig), so a figure can always be traced back to the code + version
#   that produced it. Bump on a user-visible behaviour change.
SKILL_VERSION = "1.0"


@dataclass
class Config:
    # ---------------------------------------------------------------- I/O
    input_dir: str = "."                 # where raw data lives
    output_dir: str = "./out"            # where figures + tables are written
    formats: Sequence[str] = ("pdf", "png")   # vector first; png is the raster proof
    #   one-line PLAIN-LANGUAGE description of what the figure shows (species, band,
    #   design, n). Stamped into the figure file's metadata by save_fig so the
    #   file itself carries provenance — "figure X shows what, made by what". Left
    #   empty it falls back to the basename; set it for a report-grade deliverable.
    description: str = ""

    # ---------------------------------------------------------------- style
    journal: str = "nature"              # "nature" | "acs" | "ieee" | "general"  (sets widths/fonts)
    column: str = "single"               # "single" (~3.5in) | "double" (~7.2in)
    palette: str = "okabe_ito"           # "okabe_ito" | "colorblind"  (both colour-blind safe)
    #   grayscale-distinguishing channel for OVERLAID traces (colour alone isn't
    #   grayscale-safe). HOUSE RULE: spectra are SOLID by default ("none") — dashes
    #   fight high-frequency IR data and are opt-in only. Set "linestyle" (dashes)
    #   for a SMOOTH single-band overlay if you want it, or "marker" (sparse markers
    #   via markevery) for a choppy overlay. Waterfall/offset plots distinguish by
    #   POSITION, so they stay solid regardless of this setting.
    redundancy: str = "none"             # "none" (solid, house default) | "linestyle" | "marker"
    #   how a WATERFALL identifies its stacked traces. "edge" = direct labels at the
    #   right axis margin aligned to each trace (grayscale-safe, no floating text over
    #   data — the house default); "legend" = a legend box; "none" = caption only.
    waterfall_legend: str = "edge"       # "edge" | "legend" | "none"
    base_fontsize: float = 9.0           # pt at FINAL print size; floor enforced at 6 pt
    dpi_raster: int = 600                # for png/tiff; vector ignores this
    figsize: Optional[tuple] = None      # (w,h) inches; None => derived from journal/column

    # ---------------------------------------------------------------- spectra (FTIR + PXRD)
    domain: str = "ftir"                 # "ftir" | "pxrd"  -> flips axis/intensity conventions
    #   FTIR: x = wavenumber (cm-1), axis INVERTED, y = absorbance
    #   PXRD: x = 2theta (deg),      axis NORMAL,   y = intensity (counts / normalised)
    x_limits: Optional[tuple] = None     # (left,right) in x units; None => full range.
                                         #   FTIR example (1820,1650) is already inverted.
    baseline: Optional[str] = "arpls"    # "arpls" | "rubberband" | None
    arpls_lam: float = 1e5               # arPLS smoothness; bigger = stiffer baseline
    #   integration windows: list of (low_anchor, high_anchor, name) in x units.
    #   A local LINEAR baseline is drawn between the two anchors and applied
    #   IDENTICALLY to every spectrum (this is asserted by py).
    integration_windows: List[Tuple[float, float, str]] = field(default_factory=list)
    #   how the band baseline is drawn under those windows:
    #     "per_window" (default) - a separate local linear baseline per window, anchored
    #         at that window's own two edges. Correct for ISOLATED bands.
    #     "shared" - ONE linear baseline across the whole envelope (from the leftmost
    #         anchor to the rightmost anchor of all windows); each window then integrates
    #         its sub-interval above that single line, i.e. a vertical-drop / "drop
    #         perpendicular" partition at the window boundaries. Use for OVERLAPPING bands
    #         that sit on a shared pedestal (e.g. an analyte band adjacent to an
    #         internal-standard band) where a per-window baseline would carve out the
    #         valley and shrink the smaller band's area.
    integration_baseline: str = "per_window"   # "per_window" | "shared"
    #   --- algorithmic band anchors (removes the hand-picked window edges that make
    #   one analyst's %RSD differ from another's). With "auto", anchors are the
    #   flanking local minima found by find_anchors; detect them ONCE on a
    #   reference spectrum via auto_windows and lock them onto
    #   integration_windows (do NOT re-detect per replicate).
    integration_anchor: str = "manual"   # "manual" (use integration_windows) | "auto" (find_anchors)
    anchor_centers: List[Tuple[float, str]] = field(default_factory=list)  # (center_cm-1, name) for auto mode
    anchor_gap: float = 12.0             # cm-1 kept off the peak when searching for the flanking minimum
    anchor_maxhw: float = 60.0           # cm-1 max half-width of the anchor search, either side of center
    anchor_smooth: Tuple[int, int] = (13, 3)   # SG (window,poly) used to LOCATE anchors only (integrate on RAW)
    #   --- band quantification metric (what number feeds %RSD / the calibration)
    quant_metric: str = "area"           # "area" | "height" | "deriv2"  (deriv2 is anchor-free)
    #   SG (window,poly) for the 2nd-derivative ("deriv2") metric. savgol_filter(deriv=2)
    #   IS the smoothing+differentiation in one pass (don't smooth separately first).
    #   A 2nd derivative amplifies noise hard, so the window must be set to the band:
    #   ~1-2x its FWHM is a safe start; too narrow = noisy, too wide = the band is
    #   smeared and the derivative's resolving advantage is lost. poly 2 and 3 give
    #   IDENTICAL d2 filters (likewise 4 and 5). Tune per peak via the robustness panel.
    deriv_smooth: Tuple[int, int] = (17, 3)
    normalize: Optional[str] = None      # "max" | "area" | None  (display normalisation)
    ftir_yaxis: str = "absorbance"       # "absorbance" | "transmittance"  (FTIR only)
    smooth: Optional[str] = None         # None | "savgol"  (display smoothing; integrate on RAW)
    savgol_window: int = 11              # Savitzky-Golay window (odd); must exceed savgol_poly
    savgol_poly: int = 3                 # Savitzky-Golay polynomial order
    #   reference stick lines for PXRD phase ID: list of (2theta, label)
    pxrd_reference_lines: List[Tuple[float, str]] = field(default_factory=list)

    # ---------------------------------------------------------------- calibration
    conf_level: float = 0.95             # for CI / prediction bands and the slope CI
    force_residual_panel: bool = True    # a residual panel is MANDATORY; do not disable lightly
    refuse_extrapolation: bool = True    # predict() raises outside the calibrated range
    lod_loq_method: str = "residual_sd"  # "residual_sd" (3.3*s/m, 10*s/m) | "blank_sd"

    # ---------------------------------------------------------------- chemometrics (PLS / PCR / PCA)
    #   Multivariate calibration (the chemometrics family). scikit-learn installs
    #   lazily (only when a model is actually fit) — an FTIR/univariate job never
    #   imports it. The figures of merit are RMSECV (cross-validated) + R2(CV).
    pls_preprocess: str = "snv"          # "none"|"snv"|"msc"|"d1"|"d2"; chain with "+", e.g. "snv+center"
    #   d1/d2 use a Savitzky-Golay window from `deriv_smooth` (shared with the spectra family).
    #   STATEFUL steps (msc/center/autoscale) are fit on TRAINING rows inside each CV fold.
    pls_max_components: int = 10         # max latent variables / PCs scanned by component_scan
    pls_scale: bool = True               # sklearn PLS autoscale (divide each column by its SD)
    cv_scheme: str = "auto"              # "auto"/"loo" (per sample) | "kfold"; pass groups= for leave-one-LEVEL-out
    cv_folds: int = 5                    # k for the "kfold" scheme
    #   parsimony: choose the FEWEST components whose RMSECV is within this fraction of
    #   the global-min RMSECV (avoids over-fitting the component count to the CV noise).
    lv_parsimony_tol: float = 0.10       # fraction above the global-min RMSECV
    conc_unit: str = "a.u."              # response/concentration unit for axis labels (set e.g. "% w/w")

    # ---------------------------------------------------------------- verification
    min_n_for_mean_bar: int = 3          # below this, charts refuse a bare mean bar (show points)
    #   robust (median/MAD) outlier flag for a set of replicate band-metrics: a point
    #   is flagged when |x - median| > outlier_mad_n * scaled_MAD (scaled_MAD =
    #   1.4826*MAD ~ a robust SD). 3.5 is a sensible default; scverse single-cell QC
    #   uses ~5 to be permissive. Robust to the very outliers it detects, unlike mean±k·SD.
    outlier_mad_n: float = 3.5           # MAD multiplier for flag_outliers_mad
    tick_overlap_tol_px: float = 2.0     # audit_layout flags tick labels closer than this
    clip_tol_px: float = 2.0             # audit_layout flags non-tick text past the canvas edge by this
    strict: bool = True                  # True => any FAIL gate raises; False => warns only

    # ---------------------------------------------------------------- crystallography (CIF)
    #   The crystal family: CIF -> validation table + calculated PXRD + 3D viewer.
    #   gemmi / Dans_Diffraction / pymatgen install lazily at run time (heavy deps).
    cif_path: Optional[str] = None       # input CIF
    cif_block: Optional[str] = None      # data_ block name OR index; None => require a single-block file
    # -- validation / geometry
    xh_normalize: bool = True            # apply neutron X-H distances BEFORE H-bond geometry
    #   density triple-check decides its own severity: expansion-implicated => hard-fail
    #   (correctness, ignores strict); incomplete-list / inconsistent => loud banner.
    # -- H-bonds & contacts (donor and acceptor sets are SEPARATE)
    hbond_donors: Sequence[str] = ("N", "O")               # O-H / N-H; weak mode adds C
    hbond_acceptors: Sequence[str] = ("N", "O", "F", "S", "Cl")
    hbond_weak: bool = False             # opt-in C-H donors; emits an INFO nudge if the floor is still 120
    hbond_angle_min: float = 120.0       # deg; authoritative, never auto-mutated. ~90 for weak donors
    # -- PXRD from CIF
    pxrd_wavelength: Optional[float] = None  # A; None => the CIF's declared wavelength; numeric overrides
    pxrd_two_theta_min: float = 5.0      # deg; keeps the (000) at 0 deg OFF-GRID
    pxrd_two_theta_max: float = 50.0     # deg
    pxrd_peak_width: float = 0.02        # Q (A^-1) NOT deg! Dans FWHM; ~0.02 ~ 0.28 deg at low angle
    pxrd_lorentz_fraction: float = 0.5   # Dans pseudo-Voigt mix (0=Gauss, 1=Lorentz); reused as the
                                         #   pseudo-Voigt eta by simulate_pattern
    pxrd_crosscheck: str = "auto"        # "auto" (pymatgen if present else bragg) | "bragg" | "none"
    # -- PXRD realism (py): make a calc pattern match a real lab Cu-Ka scan.
    pxrd_kalpha2: bool = False           # add the Cu Ka2 line (doublet) — a real lab Cu tube emits Ka1+Ka2
    pxrd_wavelength2: Optional[float] = None  # Ka2 wavelength (A); None => Cu Ka2 1.544426
    pxrd_kalpha2_ratio: float = 0.5      # Ka2/Ka1 integrated-intensity ratio (~0.5 for Cu)
    pxrd_po_axis: Optional[Tuple[int, int, int]] = None  # March-Dollase preferred-orientation hkl (None => none)
    pxrd_march_r: float = 1.0            # March parameter (1 = no PO; <1 platy/plate, >1 needle)
    #   Caglioti instrumental broadening (U,V,W) in deg^2: FWHM^2 = U tan^2(theta) + V tan(theta) + W.
    #   Peak width GROWS with angle; refine against a standard (LaB6/Si) in practice. These are a
    #   sane lab-Cu-Ka starting point (~0.08-0.11 deg over 20-80 deg 2theta).
    pxrd_caglioti: Tuple[float, float, float] = (0.01, -0.005, 0.008)
    # -- 3D viewer
    view_renderer: str = "auto"          # "auto" (pyvista if installed, else matplotlib) | "pyvista" | "matplotlib"
    view_resolution: int = 1800          # px; pyvista raster size (structures are images, not vector)
    view_ssao: bool = True               # pyvista screen-space ambient occlusion (contact shadows -> depth)
    view_label_size: float = 1.0         # label-size multiplier (1.0 = the default size; raise to enlarge)
    view_label_offset: bool = False      # nudge labels off their atom; default places them AT the atom (halo keeps text legible)
    view_style: str = "ball_stick"       # "ball_stick" | "ellipsoid" (ORTEP ADP; needs _atom_site_aniso_U_*, else degrades)
    view_hide_ch: bool = True            # hide carbon-bound H (clean default; keep for H-bond figures)
    view_orientation: str = "pca"        # "pca" | "axis_a"|"axis_b"|"axis_c" | "vector" | "custom"
    view_vector: Optional[Tuple[float, float, float]] = None  # [u v w] direction to view down ("vector")
    view_angles: Optional[Tuple[float, float, float]] = None  # (elev,azim,roll) deg explicit camera ("custom")
    view_tilt: Tuple[float, float, float] = (0.0, 0.0, 0.0)   # (d_elev,d_azim,d_roll) nudge ON TOP of any base
    view_roll_objective: str = "hbond"   # "hbond" | "long_axis"  (auto bases only)
    view_label_atoms: str = "hetero"     # "hetero" (all N/O/S) | "all" | "none" | "hbond" (only synthon donor/acceptor atoms)
    view_interactive_html: bool = False  # also emit a py3Dmol HTML with the SAME camera
    # -- phase 2 (degrade gracefully until implemented)
    adp_probability: float = 0.50        # ellipsoid probability level (needs aniso U)
    pack_cells: Tuple[int, int, int] = (1, 1, 1)
    cell_fill: str = "molecule"          # unit-cell / packing: "molecule" (whole, by centroid) | "clip" (cell contents cut at the box)
    hbond_neighbour: str = "stub"        # render_hbond_environment neighbour extent: "whole" (full molecule) | "stub" (contact atom + 1 bonded shell) | "site" (contact atom only)
    color_by_component: bool = False      # cocrystal figures: keep the largest molecule in full element colour, desaturate the others (distinguish API vs coformer)


# ======================================================================
# from style.py
# ======================================================================
"""
py  -  the house-style LOCK.

Three jobs, all about consistency rather than "how to use matplotlib":
    apply_style(cfg)        install rcParams: journal widths, colour-blind palette,
                            font floor at 6 pt, vector-friendly defaults.
    figure(cfg, ...)        a Figure sized at FINAL print dimensions (never rescale later).
    save_fig(fig, base, cfg) export vector master + raster proof; refuses JPEG.

Degrades gracefully: if scienceplots is absent it falls back to a built-in
preset that encodes the same intent (no LaTeX requirement).
"""

# Okabe-Ito: colour-blind safe. Order chosen so the first few are maximally distinct.
OKABE_ITO = ["#000000", "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
# matplotlib's tableau-colorblind10 subset
COLORBLIND = ["#006BA4", "#FF800E", "#ABABAB", "#595959",
              "#5F9ED1", "#C85200", "#898989", "#A2C8EC"]

# Journal single/double column widths in inches.
_WIDTHS = {
    "nature":  {"single": 3.50, "double": 7.20},
    "ieee":    {"single": 3.50, "double": 7.16},
    # ACS (e.g. Analytical Chemistry): single up to 240 pt = 3.33 in;
    # double up to 504 pt = 7.0 in. Min font 4.5 pt / line 0.5 pt (we keep the
    # stricter 6 pt floor); B/W line art 1200 dpi, colour 300 dpi.
    "acs":     {"single": 3.33, "double": 7.00},
    "general": {"single": 3.50, "double": 7.00},
}
_FONT_FLOOR = 6.0  # pt at final size; nothing renders smaller than this


def _palette(cfg):
    return OKABE_ITO if cfg.palette == "okabe_ito" else COLORBLIND


# the hand-rolled preset used when scienceplots is absent (kept as a constant so
# export_mplstyle can write a COMPLETE standalone style).
_BASE_RC = {
    "axes.linewidth": 0.6, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "legend.frameon": False, "figure.autolayout": False,
}


def _main_rc(cfg):
    """The cfg-dependent house rcParams (font sizes, the colour-blind cycler with its redundant
    channel, vector-friendly savefig + font embedding). Single source of truth shared by
    apply_style and export_mplstyle."""
    from cycler import cycler
    fs = max(cfg.base_fontsize, _FONT_FLOOR)
    if cfg.base_fontsize < _FONT_FLOOR:
        import warnings
        warnings.warn(f"base_fontsize {cfg.base_fontsize} below {_FONT_FLOOR} pt floor; clamped.")
    pal = _palette(cfg)
    # redundant encoding so grayscale still separates overlaid traces. The CHANNEL depends on
    # cfg.redundancy: dashes for smooth data, sparse markers for choppy data, or none.
    red = getattr(cfg, "redundancy", "none")   # house default: spectra solid unless opted out
    line_extra = {}
    if red == "marker":
        marks = ["o", "s", "^", "D", "v", "P", "X", "*"]
        cyc = (cycler(color=pal) +
               cycler(marker=[marks[i % len(marks)] for i in range(len(pal))]))
        line_extra = {"lines.linestyle": "-", "lines.markevery": 0.08,
                      "lines.markersize": 3.0, "markers.fillstyle": "none"}
    elif red == "none":
        cyc = cycler(color=pal)
        line_extra = {"lines.linestyle": "-"}
    else:  # "linestyle" (default)
        styles = ["-", "--", "-.", ":"]
        cyc = (cycler(color=pal) +
               cycler(linestyle=[styles[i % len(styles)] for i in range(len(pal))]))
    rc = {
        "font.size": fs, "axes.titlesize": fs, "axes.labelsize": fs,
        "xtick.labelsize": fs - 1, "ytick.labelsize": fs - 1, "legend.fontsize": fs - 1,
        "font.family": "sans-serif",
        "axes.prop_cycle": cyc,
        "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42, "ps.fonttype": 42,   # embed real fonts, not paths (editable text)
        "svg.fonttype": "none",
        "lines.linewidth": 1.0, "lines.markersize": 3.5,
    }
    rc.update(line_extra)
    return rc


def apply_style(cfg):
    """Install the house  Call once before plotting."""
    # Prefer scienceplots' no-latex style; fall back to the hand-rolled preset.
    try:
        import scienceplots  # noqa: F401
        plt.style.use(["science", "no-latex"])
    except Exception:
        plt.rcParams.update(_BASE_RC)
    plt.rcParams.update(_main_rc(cfg))


def figure(cfg, nrows=1, ncols=1, height=None, height_ratios=None,
           width_ratios=None, **kw):
    """A Figure sized at FINAL print dimensions. Set the size ONCE here; never
    rescale in Word/LaTeX afterwards (that silently shrinks the pt fonts).

    For a grid the default height scales with the row:col ratio so panels stay
    roughly square instead of being squashed into a single row's height (pass
    `height` or cfg.figsize to override). `height_ratios`/`width_ratios` size the
    rows/columns relative to each other (e.g. a tall spectrum over a short residual
    strip, or a wide panel beside a narrow one)."""
    if cfg.figsize is not None:
        w, h = cfg.figsize
    else:
        w = _WIDTHS.get(cfg.journal, _WIDTHS["general"])[cfg.column]
        h = height if height is not None else w * 0.72 * nrows / ncols
    gridspec_kw = {}
    if height_ratios:
        gridspec_kw["height_ratios"] = height_ratios
    if width_ratios:
        gridspec_kw["width_ratios"] = width_ratios
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, h),
                             gridspec_kw=gridspec_kw or None, **kw)
    return fig, axes


def _figure_metadata(cfg, extra=None):
    """Provenance embedded IN the figure file, so the file itself records how it was
    made (Claude Science's 'every figure carries the exact code and environment that
    produced it', at the file level): the skill version, a plain-language description
    (cfg.description), and the producing software. Backends accept different metadata
    keys, so return a per-format dict; a format not listed gets None. Also OVERRIDES
    matplotlib's default 'Software'/date tEXt on PNG, which makes the raster proof
    byte-reproducible."""
    import sys
    desc = (getattr(cfg, "description", "") or "").strip()
    ver = f"analytical-figures v{SKILL_VERSION}"
    soft = f"{ver}; matplotlib {matplotlib.__version__}; python {sys.version.split()[0]}"
    pdf = {"Creator": ver, "Producer": soft}          # PDF/PS: fixed key set
    svg = {"Creator": ver}                            # SVG: Title/Description/Creator
    png = {"Software": soft}                           # PNG tEXt: arbitrary keys
    if desc:
        pdf["Title"] = pdf["Subject"] = desc
        svg["Title"] = svg["Description"] = desc
        png["Description"] = desc
    meta = {"pdf": pdf, "ps": pdf, "eps": pdf, "svg": svg, "svgz": svg, "png": png}
    if extra:
        meta = {k: {**v, **extra} for k, v in meta.items()}
    return meta


def save_fig(fig, basename, cfg, metadata=None):
    """Export every format in cfg.formats. Vector (pdf/svg) is the master;
    PNG/TIFF get dpi_raster; vector ignores it.
    Refuses JPEG (compression artefacts fail journal PDF checkers).

    Stamps provenance (skill version + cfg.description + producing software) into each
    file's metadata so the figure is traceable back to the code that made it; pass
    `metadata=` to add/override keys. A backend that rejects the metadata dict still
    gets the figure written (the stamp is best-effort, never fatal)."""
    os.makedirs(os.path.dirname(basename) or ".", exist_ok=True)
    written = []
    prov = _figure_metadata(cfg, metadata)
    for ext in cfg.formats:
        if ext.lower() in ("jpg", "jpeg"):
            raise ValueError("JPEG is not allowed for data figures (use pdf/svg/png/tiff)")
        path = f"{basename}.{ext}"
        dpi = cfg.dpi_raster if ext in ("png", "tiff") else None
        md = prov.get(ext.lower())
        try:
            fig.savefig(path, dpi=dpi, metadata=md)
        except (TypeError, ValueError):
            fig.savefig(path, dpi=dpi)   # backend rejected the metadata dict; write anyway
        written.append(path)
    return written


def check_figure_width(fig, cfg=None, journal=None, column=None, tol_mm=0.5):
    """Verify a figure's PHYSICAL width matches the journal/column spec, so it's provably
    submission-ready and won't be silently rescaled (rescaling shrinks the pt fonts). Returns
    (ok, message); journals quote widths in mm so the message is in mm. journal/column default
    from cfg. (ACS Anal. Chem. single 3.33 in / double 7.0 in; Nature 89 / 183 mm.)"""
    journal = journal or getattr(cfg, "journal", "general")
    column = column or getattr(cfg, "column", "single")
    spec = _WIDTHS.get(journal, _WIDTHS["general"])
    want_in = spec.get(column, spec["single"])
    got_in = float(fig.get_size_inches()[0])
    got_mm, want_mm = got_in * 25.4, want_in * 25.4
    ok = abs(got_mm - want_mm) <= tol_mm
    msg = (f"width {got_mm:.1f} mm {'==' if ok else '!='} {journal}/{column} spec "
           f"{want_mm:.1f} mm (delta {abs(got_mm-want_mm):.2f} mm, tol {tol_mm} mm)")
    return ok, msg


def export_mplstyle(cfg, path, include_figsize=True):
    """Write a standalone `.mplstyle` reproducing the house style WITHOUT importing this skill, so a
    teammate can `plt.style.use(path)` and match the figures. Encodes _BASE_RC + the cfg-dependent
    house rcParams (fonts, the colour-blind+redundant cycler, vector-friendly savefig/font
    embedding) and, optionally, the journal single/double figure width. Hex colours are written
    WITHOUT '#' (a '#' starts a comment in a .mplstyle)."""
    import os
    rc = dict(_BASE_RC)
    rc.update(_main_rc(cfg))
    lines = ["# analytical-figures house style (auto-exported by export_mplstyle;",
             "# no skill import needed -- plt.style.use(this_file))."]
    for k, v in rc.items():
        if k == "axes.prop_cycle":
            keyed = v.by_key()
            parts = ["cycler('color', [" + ", ".join(f"'{c.lstrip('#')}'" for c in keyed["color"]) + "])"]
            for ch in ("linestyle", "marker"):
                if ch in keyed:
                    parts.append(f"cycler('{ch}', [" + ", ".join(f"'{s}'" for s in keyed[ch]) + "])")
            lines.append("axes.prop_cycle: " + " + ".join(parts))
        else:
            lines.append(f"{k}: {v}")
    if include_figsize:
        w = _WIDTHS.get(cfg.journal, _WIDTHS["general"])[cfg.column]
        lines.append(f"figure.figsize: {w}, {round(w * 0.72, 3)}")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# ===================================================== composite / multi-panel
# Aligned panel letters (a/b/c) + a layout safety-net for composite figures.
# The alignment trick and the constrained->tight fallback are adapted from SciPilot
# (Haojae/scipilot-figure-skill, MIT): anchor every label at its axes' (0,1) corner
# and apply ONE points offset, so same-column labels share a figure-x and same-row
# labels a figure-y -> they line up both ways regardless of differing y-tick widths.
# Rewritten here to read the house cfg and to tag each label so audit_layout
# can confirm one per panel.

# journal convention for the panel letter; cfg.journal picks the default.
_PANEL_FMT = {
    "nature":  lambda s: s,          # bold lower-case  a b c   (Nature/Cell)
    "acs":     lambda s: s,          # bold lower-case  a b c   (ACS)
    "general": lambda s: s,          # a b c
    "ieee":    lambda s: f"({s})",   # (a) (b) (c)              (IEEE/Elsevier)
}


def _letter_sequence(n):
    import string
    L = string.ascii_lowercase
    return [L[i] if i < 26 else L[i // 26 - 1] + L[i % 26] for i in range(n)]


def _grid_axes(fig):
    """Real gridspec subplots only (drop colorbars / insets -- no subplotspec),
    sorted in reading order: top row first, left to right within a row."""
    axes = [ax for ax in fig.axes if ax.get_subplotspec() is not None]
    return sorted(axes, key=lambda ax: (-round(ax.get_position().y1, 3),
                                        round(ax.get_position().x0, 3)))


def finalize_figure(fig, prefer="constrained",
                    wspace=None, hspace=None, w_pad=None, h_pad=None, rect=None):
    """Layout safety-net: settle the margins so titles/labels aren't clipped, the
    legend doesn't sit on data, and panels don't overlap. Try constrained_layout,
    fall back to tight_layout, else leave it. CALL THIS BEFORE add_panel_labels --
    panel positions must be settled before the letters are anchored to them. Returns
    the engine actually used ('constrained' | 'tight' | 'none').

    Optional spacing controls WIDEN the gutters between panels -- needed when a panel
    carries right-margin edge labels (edge_labels) that would otherwise spill
    into the neighbouring panel. wspace/hspace are inter-panel gaps as a fraction of
    the panel size; w_pad/h_pad are outer padding. `rect=(left, bottom, width, height)`
    in figure fraction confines the axes to a sub-rectangle, RESERVING an outer margin
    -- use it to leave room for edge labels at the figure's right edge (the layout
    engine can't see clip_on=False text, so it won't reserve that space itself). The
    rect is given in constrained-layout convention and converted for the tight
    fallback (which uses left,bottom,right,top)."""
    if prefer == "constrained":
        try:
            fig.set_layout_engine("constrained")
            spacing = {k: v for k, v in (("wspace", wspace), ("hspace", hspace),
                                         ("w_pad", w_pad), ("h_pad", h_pad),
                                         ("rect", rect)) if v is not None}
            if spacing:
                fig.get_layout_engine().set(**spacing)
            fig.canvas.draw()            # force one layout pass so positions settle
            return "constrained"
        except Exception:
            pass
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tl = {k: v for k, v in (("w_pad", w_pad), ("h_pad", h_pad)) if v is not None}
            if rect is not None:
                l, b, w, h = rect
                tl["rect"] = (l, b, l + w, b + h)   # tight uses (left,bottom,right,top)
            fig.tight_layout(**tl)
        adj = {k: v for k, v in (("wspace", wspace), ("hspace", hspace)) if v is not None}
        if adj:
            fig.subplots_adjust(**adj)
        return "tight"
    except Exception:
        return "none"


def _auto_x_offset(fig, axs, pad_pt=4.0, fallback=-20.0):
    """ONE shared leftward offset (points) that clears the WIDEST left-hand furniture
    (y-tick labels + y-axis label) among the panels, so the letter never lands on a
    wide y-tick -- while staying a single shared value so the letters still line up.
    Clamped so the anchor can't run past the figure's left edge. Falls back to the
    fixed default if it can't measure."""
    try:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        dpi = fig.dpi
        worst_px = 0.0          # furthest any furniture sticks left of its own spine
        min_spine_px = float("inf")
        for ax in axs:
            spine_px = ax.get_window_extent(r).x0
            min_spine_px = min(min_spine_px, spine_px)
            lefts = [t.get_window_extent(r).x0 for t in ax.get_yticklabels() if t.get_text()]
            if ax.yaxis.label.get_text():
                lefts.append(ax.yaxis.label.get_window_extent(r).x0)
            if lefts:
                worst_px = max(worst_px, spine_px - min(lefts))
        need_pt = -(worst_px * 72.0 / dpi + pad_pt)
        # clamp: keep the anchor (label right edge) at least ~8 px inside the canvas
        floor_pt = (8.0 - min_spine_px) * 72.0 / dpi
        return max(need_pt, floor_pt)
    except Exception:
        return fallback


def add_panel_labels(fig, cfg=None, axes=None, labels=None, style=None,
                     x_offset_pt=-20.0, y_offset_pt=2.0, fontsize=None,
                     fontweight="bold", color="black"):
    """Place aligned a/b/c panel letters on a composite figure. Each letter is
    anchored at its panel's top-left (axes fraction (0,1)) and pushed by ONE shared
    points offset, so letters line up across panels whatever the tick-label widths.
    Placing the letter is form (the skill's job); WHICH panel says what stays in the
    caption (the model's). Letters are tagged gid='panel-label' so audit_layout can
    confirm one per panel. Returns the placed Text objects.

    x_offset_pt: a fixed points offset (default -20.0), or "auto" to derive one
    shared offset that clears the widest y-tick/label furniture across the panels
    (handy when some panels carry wide y-ticks like 0.0040 and others none).

    style defaults from cfg.journal ('nature'/'general' -> 'a', 'ieee' -> '(a)');
    pass labels=[...] to override the letters entirely. Run finalize_figure first."""
    axs = list(axes) if axes is not None else _grid_axes(fig)
    if not axs:
        return []
    if x_offset_pt == "auto":
        x_offset_pt = _auto_x_offset(fig, axs)
    if labels is None:
        if style is None:
            style = getattr(cfg, "journal", "general") if cfg else "general"
        fmt = _PANEL_FMT.get(style, _PANEL_FMT["general"])
        labels = [fmt(s) for s in _letter_sequence(len(axs))]
    elif len(labels) < len(axs):
        raise ValueError(f"{len(labels)} labels for {len(axs)} panels")
    if fontsize is None:
        fontsize = plt.rcParams.get("axes.labelsize", 9)
    placed = []
    for ax, lab in zip(axs, labels):
        t = ax.annotate(lab, xy=(0, 1), xycoords="axes fraction",
                        xytext=(x_offset_pt, y_offset_pt), textcoords="offset points",
                        ha="right", va="bottom", fontsize=fontsize,
                        fontweight=fontweight, color=color,
                        annotation_clip=False)     # let it sit outside the axes
        t.set_gid("panel-label")
        placed.append(t)
    return placed


def panel_letter(ax, label, loc="upper left", pad=0.04, fontsize=None,
                 fontweight="bold", color="black"):
    """Place ONE panel letter INSIDE an axes corner (default top-left). Use this for a
    SECTIONED or NESTED composite (subfigures / nested gridspecs) where
    add_panel_labels' single-flat-grid enumeration doesn't apply: loop the charts in
    reading order and call panel_letter(ax, s) on each. A coupled pair (e.g. a
    calibration over its residual strip) is ONE chart -> letter only its main axes.
    Tagged gid='panel-label' like add_panel_labels; `label` is placed verbatim (pass
    'a' or '(a)' to taste). `loc` is one of upper/lower x left/right.

    NB audit_layout's one-letter-per-panel check assumes a flat grid and will
    over-count a nested/sectioned figure -- confirm the lettering by reading the PNG."""
    if fontsize is None:
        fontsize = plt.rcParams.get("axes.labelsize", 9)
    va = "top" if "upper" in loc else "bottom"
    ha = "right" if "right" in loc else "left"
    x = (1 - pad) if ha == "right" else pad
    y = (1 - pad) if va == "top" else pad
    t = ax.text(x, y, label, transform=ax.transAxes, ha=ha, va=va, zorder=6,
                fontsize=fontsize, fontweight=fontweight, color=color)
    t.set_gid("panel-label")
    return t


# ======================================================================
# from verify.py
# ======================================================================
"""
py  -  the part no other figure skill has: verification in two tiers.

TIER 1  data-pipeline gates  (catch wrong NUMBERS before they reach a figure)
    check_trace(x, y, cfg)         finite / smooth / point-count / monotonic axis
    check_ingest(path, cfg)        raw-file size integrity
    same_processing(records)       assert identical baseline+window across a batch
    summarize(values, cfg)         honest stats: mean, SD, SEM, t-based CI, labelled

TIER 2  visual QA loop  (catch broken-LOOKING figures before final export)
    render_preview(fig, path)      rasterise to PNG so the model can *look* at it
    audit_layout(fig, cfg)         deterministic: glyphs / clipping / data escaping the axes /
                                   x+y tick overlap / panel letters
    READ_IMAGE_CHECKLIST           perceptual: the model opens the PNG and checks these

TIER 3  the critic  (catch numbers that don't trace to the code)
    check_number_provenance(src)   static: a figure-of-merit on the figure / in the
                                   caption must be INTERPOLATED from the computed value,
                                   never hand-typed (it silently desyncs on a refit)

Plus a robust univariate outlier flag that feeds Tier 1:
    flag_outliers_mad(values, cfg) median/MAD flag for replicate band-metrics (the
                                   univariate analog of diagnostics)

The intended loop:  draw -> render_preview -> audit_layout (fix any FAIL)
-> open the PNG with the Read tool, walk READ_IMAGE_CHECKLIST -> fix -> re-render
-> only then save_fig() the vector master (bundle.py runs the Tier-3 critic at handoff).
"""

SEV = {"INFO": 0, "WARN": 1, "FAIL": 2}


class GateError(RuntimeError):
    """Raised when a FAIL gate trips under cfg.strict."""


def _resolve(findings, cfg, context=""):
    """Print findings; raise on any FAIL if cfg.strict."""
    worst = max((SEV[s] for s, _ in findings), default=0)
    for sev, msg in findings:
        print(f"  [{sev}] {context}{': ' if context else ''}{msg}")
    if worst == SEV["FAIL"] and getattr(cfg, "strict", True):
        fails = "; ".join(m for s, m in findings if s == "FAIL")
        raise GateError(f"{context}: {fails}")
    return findings


# ============================================================ TIER 1: data gates
def check_trace(x, y, cfg, name="trace", min_points=8, require_monotonic=True):
    """Gate a single spectrum/trace. Returns list of (severity, msg); any FAIL
    means the decode/ingest is wrong and must not be plotted or integrated.
    min_points defaults to 8 (a spectrum with fewer is almost certainly a decode
    error); pass a smaller value for legitimately short series like calibration
    standards. require_monotonic defaults to True (a spectrum's x-axis must be
    monotonic); set False for calibration data, whose replicate standards
    legitimately repeat x-values."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    f = []
    if x.size != y.size:
        f.append(("FAIL", f"x/y length mismatch ({x.size} vs {y.size})"))
        return _resolve(f, cfg, name)
    if x.size < min_points:
        f.append(("FAIL", f"only {x.size} points (min {min_points}) - decode almost certainly wrong"))
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        f.append(("FAIL", "non-finite values present (NaN/inf)"))
    if require_monotonic:
        dx = np.diff(x)
        if not (np.all(dx > 0) or np.all(dx < 0)):
            f.append(("FAIL", "x axis is not monotonic (concatenation/parse error?)"))
    # crude noise/spike check: a single sample dominating the dynamic range
    if y.size and np.ptp(y) > 0:
        spikes = np.abs(np.diff(y, 2))
        if spikes.size and spikes.max() > 0.5 * np.ptp(y):
            f.append(("WARN", "large single-sample spike - check for a dead pixel/cosmic ray"))
    if not f:
        f.append(("INFO", f"{x.size} pts, monotonic, finite - OK"))
    return _resolve(f, cfg, name)


def check_ingest(path, cfg, min_bytes=64):
    """Raw-file integrity before parsing."""
    f = []
    if not os.path.exists(path):
        f.append(("FAIL", f"missing file: {path}"))
    elif os.path.getsize(path) < min_bytes:
        f.append(("FAIL", f"file too small ({os.path.getsize(path)} B) - truncated?"))
    else:
        f.append(("INFO", f"{os.path.basename(path)}: {os.path.getsize(path)} B"))
    return _resolve(f, cfg, "ingest")


def same_processing(records, cfg):
    """records: iterable of dicts each describing how a trace was processed, e.g.
    {"baseline":"arpls","arpls_lam":1e5,"windows":[...]}. Asserts every member of
    a batch was processed IDENTICALLY - rep-to-rep differences must be real, not
    an artefact of inconsistent baselines/windows."""
    records = list(records)
    f = []
    if len(records) <= 1:
        f.append(("INFO", "single trace - nothing to compare"))
        return _resolve(f, cfg, "batch")
    ref = records[0]
    keys = ("baseline", "arpls_lam", "windows", "normalize")
    for i, r in enumerate(records[1:], 1):
        diffs = [k for k in keys if r.get(k) != ref.get(k)]
        if diffs:
            f.append(("FAIL", f"trace {i} differs from trace 0 in {diffs} "
                              f"- batch not processed identically"))
    if not f:
        f.append(("INFO", f"{len(records)} traces share identical processing"))
    return _resolve(f, cfg, "batch")


def summarize(values, cfg):
    """Honest stats for an error bar. Returns a dict AND a caption string that
    NAMES the statistic, n, and (for CI) the t-multiplier. Never present an
    error bar without saying which statistic it is."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    n = v.size
    mean = float(np.mean(v)) if n else float("nan")
    sd = float(np.std(v, ddof=1)) if n > 1 else float("nan")
    sem = sd / np.sqrt(n) if n > 1 else float("nan")
    p = cfg.conf_level
    try:
        from scipy import stats
        t = float(stats.t.ppf(0.5 + p / 2, df=n - 1)) if n > 1 else float("nan")
    except Exception:
        # normal-approx fallback if scipy absent
        z = {0.90: 1.645, 0.95: 1.960, 0.99: 2.576}.get(round(p, 2), 1.960)
        t = z if n > 1 else float("nan")
    ci = t * sem if n > 1 else float("nan")
    caption = (f"mean ± {int(p*100)}% CI (t={t:.3f}·SEM), n={n}"
               if n > 1 else f"single value, n={n}")
    return {"mean": mean, "sd": sd, "sem": sem, "ci_halfwidth": ci,
            "n": n, "t": t, "caption": caption}


def flag_outliers_mad(values, cfg=None, n_mads=None, name="metric"):
    """Robust median/MAD outlier flag for a set of REPLICATE band-metrics (areas,
    heights, %-values). A point is flagged when its robust z, |x - median| /
    (1.4826*MAD), exceeds n_mads. This is the univariate, raw-spectrum analog of the
    multivariate `diagnostics` gate: it turns an ad-hoc "flag-don't-hide"
    outlier call (an under-loaded or a diluted replicate) into a principled, reportable
    rule. Robust BY CONSTRUCTION - a mean±k·SD rule is inflated by the very outlier it
    should catch (the outlier drags the mean and the SD); the median and MAD are not.
    n_mads defaults to cfg.outlier_mad_n (3.5).

    Small-n caveat: the MAD of 3 replicates is itself noisy, so this flags best over a
    POOLED set (a whole level, or all preps of one composition), not 3 lone reps - it
    WARNs when n<5. NEVER raises (flag-don't-hide: a flagged point is SURFACED for you
    to retain / re-measure / footnote, never silently dropped). Returns a dict with the
    boolean `mask` (same length as `values`, NaNs never flagged), the robust z, the
    accept band (lower, upper), and a caption naming the rule.

    Idiom for the figure: plot the raw replicate points, shade [lower, upper] as the
    accept band, and mark the flagged points - so the reader SEES why a point was
    called (this is `references/verification.md`'s before/after threshold overlay)."""
    v = np.asarray(values, float)
    finite = np.isfinite(v)
    vf = v[finite]
    n = vf.size
    if n_mads is None:
        n_mads = getattr(cfg, "outlier_mad_n", 3.5) if cfg is not None else 3.5
    med = float(np.median(vf)) if n else float("nan")
    abs_dev = np.abs(vf - med)
    mad = float(np.median(abs_dev)) if n else float("nan")
    scaled = 1.4826 * mad                       # MAD -> SD (normal consistency constant)
    note = ""
    if n and scaled == 0.0:
        # >=50% of the points are identical -> MAD collapses to 0. Fall back to the
        # (scaled) mean absolute deviation so a minority of outliers is still catchable.
        scaled = 1.2533 * float(np.mean(abs_dev))   # sqrt(pi/2)*MeanAD -> SD
        note = "; MAD=0, used mean-abs-dev fallback"
    mask = np.zeros(v.shape, bool)
    z_full = np.full(v.shape, np.nan)
    if n and scaled > 0:
        z_full[finite] = abs_dev / scaled
        mask[finite] = z_full[finite] > n_mads
    lower = med - n_mads * scaled if scaled > 0 else float("nan")
    upper = med + n_mads * scaled if scaled > 0 else float("nan")
    idx = [int(i) for i in np.nonzero(mask)[0]]
    n_flag = len(idx)
    findings = []
    if 0 < n < 5:
        findings.append(("WARN", f"{name}: n={n} < 5 - MAD is noisy at low n; flag over a "
                                 f"POOLED set and always read the flag WITH the raw points"))
    findings.append(("INFO", f"{name}: {n_flag}/{v.size} flagged at robust z>{n_mads} "
                             f"(median={med:.4g}, robust SD={scaled:.4g}, indices {idx}{note})"))
    for s, m in findings:                        # print-only; a flag is surfaced, never raised
        print(f"  [{s}] mad-outlier: {m}")
    caption = (f"outliers by robust z = |x-median|/(1.4826·MAD) > {n_mads}: "
               f"{n_flag} of {v.size} flagged (n={n})")
    return {"mask": mask, "robust_z": z_full, "median": med, "mad": mad,
            "scaled_mad": scaled, "lower": lower, "upper": upper,
            "n_flagged": n_flag, "flagged_index": idx, "n": n, "n_mads": n_mads,
            "caption": caption}


def adjust_pvalues(pvals, method="holm"):
    """Multiple-comparison correction for a FAMILY of p-values (returned in the input order).
        'bonferroni'  FWER control, p_adj = min(1, m*p)
        'holm'        step-down Holm-Bonferroni (uniformly more powerful than Bonferroni)
        'bh'/'fdr_bh' Benjamini-Hochberg FDR (step-up)
    Matches statsmodels.stats.multitest.multipletests (bonferroni / holm / fdr_bh)."""
    p = np.asarray(pvals, float)
    m = p.size
    if m == 0:
        return p.copy()
    order = np.argsort(p)
    ranked = p[order]
    if method == "bonferroni":
        adj_sorted = np.minimum(ranked * m, 1.0)
    elif method == "holm":
        terms = (m - np.arange(m)) * ranked
        adj_sorted = np.minimum(np.maximum.accumulate(terms), 1.0)
    elif method in ("bh", "fdr_bh"):
        ranks = np.arange(1, m + 1)
        terms = ranked * m / ranks
        adj_sorted = np.minimum(np.minimum.accumulate(terms[::-1])[::-1], 1.0)
    else:
        raise ValueError(f"unknown method {method!r} (bonferroni|holm|bh)")
    adj = np.empty(m)
    adj[order] = adj_sorted
    return adj


def multiplicity_check(pvals, cfg=None, alpha=0.05, method="holm", labels=None, context="multiplicity"):
    """Gate for annotating a FAMILY of p-values together (intercept-bias across levels, lack-of-fit
    across operators, several pairwise model comparisons). Returns the adjusted p-values + a caption,
    and WARNs that RAW p-values inflate significance when m>1 — annotate the ADJUSTED ones in the
    figure. (Encodes the 'visual claims must match the evidence' rule for multiple comparisons.)"""
    p = np.asarray(pvals, float)
    adj = adjust_pvalues(p, method)
    reject = adj <= alpha
    m = p.size
    f = []
    if m > 1:
        n_raw = int(np.sum(p <= alpha)); n_adj = int(np.sum(reject))
        if n_raw > n_adj:
            f.append(("WARN", f"{m} p-values: {n_raw} significant at RAW alpha but only {n_adj} after "
                              f"{method} correction — annotate ADJUSTED p (raw inflates significance)"))
        else:
            f.append(("INFO", f"{m} p-values, {method}-adjusted: {n_adj} significant at alpha={alpha}"))
    else:
        f.append(("INFO", "single p-value — no multiplicity correction needed"))
    if cfg is not None:
        _resolve(f, cfg, context)
    else:
        for s, msg in f:
            print(f"  [{s}] {context}: {msg}")
    return {"raw": p, "adjusted": adj, "reject": reject, "method": method, "alpha": alpha,
            "labels": list(labels) if labels is not None else None,
            "caption": f"{method}-adjusted p (family of {m}, alpha={alpha})"}


# ===================================================== TIER 1b: baseline robustness
def _locked_metric(cfg, center, reference):
    """Build a band_metric(x, y, method) closure that, for the anchor-based metrics
    (area/height), uses anchors detected ONCE on `reference` and locked - never
    re-detected per spectrum. deriv2 is anchor-free. If `reference` is None the
    area/height anchors fall back to per-spectrum detection (discouraged: drifts)."""
    locked = None
    if reference is not None:
        rx, ry = reference
        locked = find_anchors(rx, ry, center, gap=cfg.anchor_gap,
                                      maxhw=cfg.anchor_maxhw, smooth=cfg.anchor_smooth)
    def metric(x, y, mth):
        if mth == "deriv2":
            return band_metric(x, y, "deriv2", center=center,
                                       deriv_smooth=cfg.deriv_smooth)
        return band_metric(x, y, mth, anchors=locked, center=center,
                                   gap=cfg.anchor_gap, maxhw=cfg.anchor_maxhw,
                                   anchor_smooth=cfg.anchor_smooth)
    return metric, locked


def baseline_robustness(reps, center, cfg, methods=("area", "height", "deriv2"),
                        reference=None, context="band"):
    """ICH Q2 *robustness*, made concrete and automatic. Recompute one band's
    quantity under several baseline/metric methods for ONE replicate set, and
    report the per-method precision (%RSD). The point: show the conclusion
    ("precision is acceptable") does not hinge on the analyst's baseline choice -
    the exact failure mode where two analysts get different %RSD from the same
    

    reps      : iterable of (x, y) replicate spectra of the SAME sample.
    center    : band centre (cm-1), handed to band_metric.
    reference : (x, y) spectrum on which area/height anchors are detected ONCE and
                locked (recommended - pass the pure-analyte trace or a batch mean).
    Returns {method: {"values": [...], "summary": summarize(...), "rsd_pct": float}}.
    Findings WARN if the methods' %RSD disagree by more than ~3x (a sign the band
    model matters and you must state which metric you locked, and why)."""
    metric, _ = _locked_metric(cfg, center, reference)
    reps = [(np.asarray(x, float), np.asarray(y, float)) for x, y in reps]
    out, rsds, f = {}, [], []
    for mth in methods:
        vals = [metric(x, y, mth) for x, y in reps]
        s = summarize(vals, cfg)
        rsd = 100 * s["sd"] / s["mean"] if s["mean"] else float("nan")
        out[mth] = {"values": vals, "summary": s, "rsd_pct": rsd}
        f.append(("INFO", f"{mth}: %RSD={rsd:.2f} (mean={s['mean']:.4g}, n={s['n']})"))
        if np.isfinite(rsd):
            rsds.append(rsd)
    if len(rsds) >= 2 and min(rsds) > 0 and max(rsds) / min(rsds) > 3.0:
        f.append(("WARN", f"{context}: methods disagree on %RSD by >3x "
                          f"({min(rsds):.1f}-{max(rsds):.1f}%) - state and justify the locked metric"))
    _resolve(f, cfg, f"robustness[{context}]")
    return out


def calibration_robustness(levels, center, cfg, methods=("area", "height", "deriv2"),
                           check=None, reference=None):
    """The ICH Q2 robustness *panel* for 'which baseline?'. Run the WHOLE
    calibration under each band metric and tabulate the spread, so the locked
    primary metric is chosen from evidence, not taste.

    levels    : dict {concentration(float): [(x, y), ...replicate spectra]}.
    center    : band centre (cm-1).
    check     : optional (true_conc, [(x, y), ...]) independent sample -> % recovery.
    reference : (x, y) on which area/height anchors are detected ONCE and locked.
    Returns {method: {slope, intercept, r2, lod, loq, recovery, signal_rsd_mean}}.
    Prints a one-line-per-method comparison so the trade-offs are visible."""
    metric, _ = _locked_metric(cfg, center, reference)
    rows, f = {}, []
    for mth in methods:
        concs, sigs, per_level_rsd = [], [], []
        for c in sorted(levels):
            vals = [metric(x, y, mth) for x, y in levels[c]]
            for v in vals:
                concs.append(c); sigs.append(v)
            s = summarize(vals, cfg)
            if s["mean"] and len(vals) > 1:
                per_level_rsd.append(100 * s["sd"] / s["mean"])
        model = fit(np.array(concs, float), np.array(sigs, float), cfg)
        ll = lod_loq(model, cfg)
        rec = float("nan")
        if check is not None:
            true_c, creps = check
            csig = np.mean([metric(x, y, mth) for x, y in creps])
            pred = (csig - model["intercept"]) / model["slope"]
            rec = 100 * pred / true_c if true_c else float("nan")
        rows[mth] = {"slope": model["slope"], "intercept": model["intercept"],
                     "r2": model["r2"], "lod": ll["lod"], "loq": ll["loq"],
                     "recovery": rec, "signal_rsd_mean": float(np.mean(per_level_rsd)) if per_level_rsd else float("nan")}
        f.append(("INFO", f"{mth:7s} R2={model['r2']:.4f}  LOD={ll['lod']:.3g}  "
                          f"meanRSD={rows[mth]['signal_rsd_mean']:.1f}%  recovery={rec:.1f}%"))
    _resolve(f, cfg, "calibration-robustness")
    return rows


# ============================================================ TIER 2: visual QA
def render_preview(fig, path, dpi=200):
    """Rasterise to a quick PNG so the model can OPEN IT WITH THE READ TOOL and
    look. Low dpi on purpose - this is a proof, not the master."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=dpi)
    return path


def _glyph_render_warnings(fig):
    """Render once and harvest missing-glyph reports from BOTH matplotlib channels
    (older versions warn, newer ones log) - catches substituted tofu that leaves no
    U+FFFD in the text (minus signs, CJK, odd symbols). The U+FFFD scan alone misses
    those. Side effect: realises text extents for the measurements that follow."""
    import io, logging, warnings
    markers = ("missing from", "Glyph", "findfont")
    msgs = []

    class _H(logging.Handler):
        def emit(self, rec):
            m = rec.getMessage()
            if any(k in m for k in markers):
                msgs.append(m)

    lg = logging.getLogger("matplotlib")
    h = _H(); prev = lg.level
    lg.setLevel(logging.WARNING); lg.addHandler(h)
    try:
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter("always")
            buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=100); buf.close()
        for w in wl:
            s = str(w.message)
            if any(k in s for k in markers):
                msgs.append(s)
    finally:
        lg.removeHandler(h); lg.setLevel(prev)
    seen, out = set(), []
    for m in msgs:
        if m not in seen:
            seen.add(m); out.append(m)
    return out


def _ticks_overlap(labels, renderer, axis, tol):
    """Do adjacent tick labels collide? axis='x' tests horizontally, 'y' vertically."""
    boxes = []
    for t in labels:
        try:
            if t.get_text():
                boxes.append(t.get_window_extent(renderer))
        except Exception:
            pass
    if len(boxes) < 2:
        return False
    if axis == "x":
        boxes.sort(key=lambda b: b.x0)
        return any(a.x1 - b.x0 > tol for a, b in zip(boxes, boxes[1:]))
    boxes.sort(key=lambda b: b.y0)
    return any(a.y1 - b.y0 > tol for a, b in zip(boxes, boxes[1:]))


def _escaping_lines(ax, min_frac=0.01):
    """Lines whose DATA leaves the y-limits inside the visible x-range.

    Judged on y only, and only for points already on-screen in x: clipping in x is normally a
    deliberate zoom, whereas a trace dropping out of the bottom of a panel almost never is. The
    motivating failure was a waterfall normalised by subtracting the MEDIAN - on a sloping
    baseline the window minimum then sits below zero and the lowest trace runs off the panel,
    which every other check in this module passes silently because no TEXT is clipped.

    Returns (fraction of the axis height escaped, n points off-panel, n visible, label) per line.
    """
    (xlo, xhi) = sorted(ax.get_xlim())
    (ylo, yhi) = sorted(ax.get_ylim())
    span = yhi - ylo
    if not np.isfinite(span) or span <= 0:
        return []
    out = []
    for ln in ax.get_lines():
        # clip_on=False is an explicit opt-out; a blended transform means axvline/axhline, whose
        # y data are axes fractions and are meant to span the panel.
        if not ln.get_visible() or not ln.get_clip_on():
            continue
        if ln.get_transform() is not ax.transData:
            continue
        x = np.asarray(ln.get_xdata(orig=False), dtype=float)
        y = np.asarray(ln.get_ydata(orig=False), dtype=float)
        if x.size < 2 or x.size != y.size:
            continue
        m = np.isfinite(x) & np.isfinite(y) & (x >= xlo) & (x <= xhi)
        if int(m.sum()) < 2:
            continue
        yy = y[m]
        excess = max(ylo - yy.min() if yy.min() < ylo else 0.0,
                     yy.max() - yhi if yy.max() > yhi else 0.0)
        if excess / span > min_frac:
            out.append((excess / span, int(((yy < ylo) | (yy > yhi)).sum()), int(m.sum()),
                        ln.get_label()))
    return out


def audit_layout(fig, cfg=None):
    """Deterministic faults a program CAN catch (perceptual ones are for the Read
    tool). Single figures AND composites: missing-glyph tofu, off-canvas clipping of
    titles/labels/annotations, DATA running outside the y-limits, overlapping x AND y
    tick labels, empty axes, empty legends, unit-less numeric axes, and - for
    multi-panel figures - that every panel carries exactly one a/b/c letter (place
    them with add_panel_labels)."""
    import matplotlib.text as mtext
    f = []
    tol = getattr(cfg, "tick_overlap_tol_px", 2.0) if cfg else 2.0
    clip_tol = getattr(cfg, "clip_tol_px", 2.0) if cfg else 2.0
    # missing glyphs: intercept the render warning channels (catches tofu with no
    # U+FFFD); this render also realises the text extents used below.
    for g in _glyph_render_warnings(fig)[:2]:
        f.append(("FAIL", f"missing glyph in render - boxes/tofu will show: {g[:160]}"))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    W, H = float(fig.bbox.width), float(fig.bbox.height)
    data_axes = [ax for ax in fig.get_axes() if ax.get_subplotspec() is not None]

    for ax in fig.get_axes():
        # An axes with its frame switched off is a deliberate text/placeholder panel (a note,
        # a legend cell, an explanatory slot in a grid), not a data panel that failed to draw.
        if not ax.axison:
            continue
        if not (ax.lines or ax.collections or ax.patches or ax.images):
            f.append(("WARN", "an axes has no plotted data"))
        # missing-glyph tofu: matplotlib substitutes a box for unknown chars
        for t in ax.get_xticklabels() + ax.get_yticklabels() + [ax.xaxis.label, ax.yaxis.label]:
            if "\ufffd" in t.get_text():
                f.append(("FAIL", f"missing glyph in '{t.get_text()}' (font lacks the char)"))
        # data escaping the panel: report the worst line only, so a stack that all overflows
        # together gives one actionable message rather than one per trace
        esc = _escaping_lines(ax, getattr(cfg, "escape_tol_frac", 0.01) if cfg else 0.01)
        if esc:
            frac, n_out, n_vis, lab = max(esc)
            who = lab if lab and not str(lab).startswith("_") else "a trace"
            f.append(("WARN", f"data runs outside the y-limits: '{who}' has {n_out}/{n_vis} "
                              f"on-screen points off-panel, worst excursion "
                              f"{frac * 100:.0f} % of the axis height"
                              + (f" (+{len(esc) - 1} more line(s))" if len(esc) > 1 else "")
                              + " - widen ylim; for a waterfall normalise each trace by "
                                "subtracting the window MINIMUM, not the median"))
        if _ticks_overlap(ax.get_xticklabels(), renderer, "x", tol):
            f.append(("WARN", "x tick labels overlap - rotate, thin, or widen"))
        if _ticks_overlap(ax.get_yticklabels(), renderer, "y", tol):
            f.append(("WARN", "y tick labels overlap - fewer ticks or a taller panel"))
        # legend pointing at nothing
        leg = ax.get_legend()
        if leg is not None and len(leg.get_texts()) == 0:
            f.append(("WARN", "empty legend"))
        # units: an axis with numeric ticks must declare a unit in its label
        _UNIT_CUES = ("/", "%", "(", "cm$^{-1}$", "cm-1", "a.u", "ratio", "dimensionless",
                      "absorbance", "transmittance", "intensity", "counts", "$\\theta$",
                      "theta", "deg", "°", "r$^2$", "r²", "min", "hz", "ppm")
        for axis, getlab in ((ax.xaxis, ax.get_xlabel), (ax.yaxis, ax.get_ylabel)):
            ticks = [t.get_text() for t in axis.get_ticklabels()]
            has_numeric = any(any(ch.isdigit() for ch in t) for t in ticks)
            lab = getlab().strip()
            if has_numeric and lab and not any(c in lab.lower() for c in _UNIT_CUES):
                f.append(("WARN", f"axis '{lab}' has numeric ticks but no unit cue "
                                  f"(add a unit, or mark dimensionless/a.u./ratio)"))
    # off-canvas clipping of NON-tick text (title, axis labels, annotations, panel
    # letters). Ticks are excluded - tight/constrained layout handles those and the
    # overlap check covers them.
    tick_ids = set()
    for ax in fig.get_axes():
        for tl in (*ax.get_xticklabels(), *ax.get_xticklabels(minor=True),
                   *ax.get_yticklabels(), *ax.get_yticklabels(minor=True)):
            tick_ids.add(id(tl))
    clipped = []
    for t in fig.findobj(mtext.Text):
        if id(t) in tick_ids or not t.get_visible():
            continue
        txt = t.get_text().strip()
        if not txt:
            continue
        try:
            bb = t.get_window_extent(renderer)
        except Exception:
            continue
        if bb.x0 < -clip_tol or bb.y0 < -clip_tol or bb.x1 > W + clip_tol or bb.y1 > H + clip_tol:
            clipped.append(txt.replace("\n", " ")[:24])
    if clipped:
        f.append(("WARN", f"text may be clipped at the figure edge: "
                          f"{list(dict.fromkeys(clipped))[:5]} - run "
                          f"finalize_figure(fig) or export bbox_inches='tight'"))

    # composite: text from ONE panel spilling into a DIFFERENT panel -- manual
    # annotations / right-margin edge labels the layout engine doesn't account for
    # (constrained/tight only reserve space for ticks, axis labels, titles, legends).
    # Excludes those + the panel letters (handled by the one-per-panel check below).
    if len(data_axes) > 1:
        skip = set(tick_ids)
        for ax in fig.get_axes():
            skip.update((id(ax.xaxis.label), id(ax.yaxis.label), id(ax.title)))
            lg = ax.get_legend()
            if lg is not None:
                skip.update(id(tt) for tt in lg.get_texts())
        for lg in getattr(fig, "legends", []):
            skip.update(id(tt) for tt in lg.get_texts())
        boxes = [(ax, ax.get_window_extent(renderer)) for ax in data_axes]
        hits = []
        for t in fig.findobj(mtext.Text):
            if id(t) in skip or not t.get_visible() or t.get_gid() == "panel-label":
                continue
            if t.axes is None or not t.get_text().strip():
                continue
            try:
                bb = t.get_window_extent(renderer)
            except Exception:
                continue
            for ax, abox in boxes:
                if ax is t.axes:
                    continue
                ix = min(bb.x1, abox.x1) - max(bb.x0, abox.x0)
                iy = min(bb.y1, abox.y1) - max(bb.y0, abox.y0)
                if ix > clip_tol and iy > clip_tol:
                    hits.append(t.get_text().replace("\n", " ")[:20])
                    break
        if hits:
            f.append(("WARN", f"text overlaps a neighbouring panel: "
                              f"{list(dict.fromkeys(hits))[:5]} - widen the gutter "
                              f"(finalize_figure wspace=) or shorten the label"))

    # composite: every CHART needs exactly one aligned a/b/c letter. A residual strip (an
    # axes sharing x with a taller axes of the same figure) is part of that chart, not a
    # panel of its own: a calibration over its residual strip carries ONE letter.
    def _is_strip(ax):
        try:
            sib = [o for o in ax.get_shared_x_axes().get_siblings(ax) if o is not ax and o in data_axes]
        except Exception:
            return False
        h = ax.get_position().height
        return any(o.get_position().height > 2.0 * h for o in sib)
    charts = [ax for ax in data_axes if not _is_strip(ax)]
    if len(charts) > 1:
        n_lab = sum(1 for t in fig.findobj(mtext.Text) if t.get_gid() == "panel-label")
        if n_lab == 0:
            f.append(("WARN", f"{len(charts)}-panel figure has no panel letters - "
                              f"call add_panel_labels(fig, cfg)"))
        elif n_lab != len(charts):
            f.append(("WARN", f"panel letters ({n_lab}) != panels ({len(charts)}) - "
                              f"one per panel; check for a missing/duplicate label"))

    if not f:
        f.append(("INFO", "layout audit clean (data present, no tofu, ticks clear, "
                          "panels labelled)"))
    if cfg is not None:
        return _resolve(f, cfg, "layout")
    for s, m in f:
        print(f"  [{s}] layout: {m}")
    return f


READ_IMAGE_CHECKLIST = """\
PERCEPTUAL QA - open the preview PNG with the Read tool and verify each item.
A program cannot see these; you must actually look.

  1. Legend does not sit on top of any data or get clipped at an edge.
  2. Every axis has a label WITH UNITS (or dimensionless/a.u./ratio/%); every
     written number — annotation, residual axis, quoted value — carries its unit.
  3. Traces are distinguishable in GRAYSCALE (line style/marker differ, not just colour).
  4. FTIR: wavenumber axis runs HIGH -> LOW (inverted). PXRD: 2theta runs LOW -> HIGH.
  5. No text is clipped at the figure border; no overlapping annotations.
  6. Error bars (if any) are visible and the caption names the statistic + n.
  7. Calibration: residual panel present; residuals look random, not curved.
  8. Nothing is a default-matplotlib tell (rainbow jet colormap, untrimmed margins,
     box around legend, title doing the job of an axis label).
  9. The figure makes its intended point at FINAL print size, not just zoomed in.
 10. Composite: every panel has its a/b/c letter (one per panel) and the letters
     line up - same height across a row, same left edge down a column.
 11. Crystal structure: oriented FACE-ON (not edge-on/occluded) and the correct
     enantiomorph (not mirror-flipped); H-bonds visible + dashed, not foreshortened
     into the screen; C-bound H hidden when requested; heteroatom labels legible at
     final size; any orientation-degeneracy warning is surfaced/acknowledged.
 12. Crystal PXRD: NO monster peak at 2theta=0 (the (000) is removed); calc/exp
     overlay aligned on 2theta (a small offset is zero-point/cell, not a phase miss).
 13. Packing / unit cell: cell box drawn with a/b/c labelled; molecules whole (not
     chopped at the cell edge); the motif the figure exists to show is legible.
"""


# ============================================================ TIER 3: the critic
# A static "reviewer" pass over the analysis SCRIPT itself - not the data (Tier 1),
# not the pixels (Tier 2). It targets one silent failure mode: a figure-of-merit
# printed onto the figure or into the caption as a HAND-TYPED LITERAL that has drifted
# from what the code computed - write "R2 = 0.98" once, refit, and the code now
# computes 0.91 while the text still says 0.98. Adapted from the actor-critic reviewer
# in Anthropic's Claude Science (June 2026), which flags "untraceable numbers" and
# "figures that don't match their underlying code" - done our way: static,
# deterministic, and keyed on the skill's OWN figures-of-merit vocabulary, so it fires
# on a typed statistic but NOT on a fixed physical identifier ("1748 cm-1", a band
# centre - a constant, not a result). The rule is FORM: a computed FoM must be
# INTERPOLATED from the value the pipeline produced (an f-string), never pasted.

# FoM tokens that denote a COMPUTED statistic (so a number beside one must trace to code).
_FOM_PATTERNS = [
    r"R\s*\^?\s*2\b", r"R²", r"r-?squared",
    r"RMSE[CPV]{0,2}\b", r"\bLOD\b", r"\bLOQ\b", r"recover\w*",
    r"%?\s*RSD\b", r"\bCV\s*%", r"\bbias\b", r"\bslope\b", r"\bintercept\b",
    r"\bRPD\b", r"\bRPIQ\b", r"\bp\s*[=<>]", r"p-?value", r"\bn\s*=", r"\bLV\s*=",
]
_NUM = r"[-+]?\d*\.?\d+"                       # an int/float literal - the thing not to type
# matplotlib text sinks whose string args land ON the figure
_TEXT_SINKS = {"set_title", "set_xlabel", "set_ylabel", "set_zlabel", "suptitle",
               "text", "annotate", "set_label", "set_text", "figtext", "title",
               "xlabel", "ylabel", "set_suptitle"}
# variable-name hints that a string is caption/label prose bound for the document
_CAPTION_NAMES = ("caption", "footnote", "legend", "note", "annotation", "subtitle")


def _fom_literal_hits(s):
    """A FoM keyword FOLLOWED by a typed numeric literal within the same literal string
    (i.e. the number was pasted, not interpolated). Returns the offending snippets. The
    number must sit AFTER the keyword and OUTSIDE its own span - so the '2' in 'R2' is
    not mistaken for the value (that would flag the clean f-string segment 'R2 = '), and
    a unit that merely precedes a keyword ('5 mL, slope') doesn't trip it."""
    hits = []
    for pat in _FOM_PATTERNS:
        for m in re.finditer(pat, s, re.IGNORECASE):
            after = s[m.end(): m.end() + 12]              # a number just after the keyword
            if re.match(r"\s*[=:~]?\s*" + _NUM, after):
                hits.append(s[max(0, m.start() - 1): m.end() + 12].strip())
    for m in re.finditer(r"±\s*" + _NUM, s):             # a bare "± <value>" error bar
        hits.append(m.group(0))
    return list(dict.fromkeys(hits))


def check_number_provenance(source, cfg=None, context="critic"):
    """Static reviewer pass over an analysis SCRIPT (path or source text): flag any
    figure-of-merit written onto the figure / into a caption as a HAND-TYPED LITERAL
    rather than interpolated from the computed value. A FoM (R2 / RMSE / LOD / LOQ /
    recovery / %RSD / p / n / slope / bias / RPD / ±...) must come from an f-string
    ({value:.2f}) that traces to the pipeline; a fixed identifier in a label
    ("1748 cm-1", a band centre) is a constant and is NOT flagged. It scans only the
    strings that reach a reader: matplotlib text-sink arguments and assignments to a
    caption-ish variable name. WARN-only - it is a lint, not a data gate, and never
    raises; a hit is a prompt to interpolate the value. Returns the (severity, msg)
    findings.

    Example: `ax.set_title("R2 = 0.98")` -> WARN (typed);  `ax.set_title(f"R2 =
    {fit['r2']:.2f}")` -> clean (traces to the fit). Run it at handoff - bundle.py
    calls it automatically on the analysis body when it amalgamates the deliverable."""
    import ast
    if isinstance(source, str) and source.endswith(".py") and os.path.exists(source):
        src = open(source, encoding="utf-8").read()
    else:
        src = source
    findings = []
    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        findings.append(("WARN", f"could not parse the script for the critic pass ({e})"))
        return _critic_emit(findings, context)

    def scan(node, where):
        # a plain str constant, or an f-string's LITERAL segments only - the {...}
        # parts are interpolated (== fine), so they are never scanned.
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            segs = [node.value]
        elif isinstance(node, ast.JoinedStr):
            segs = [v.value for v in node.values
                    if isinstance(v, ast.Constant) and isinstance(v.value, str)]
        else:
            return
        for seg in segs:
            for snip in _fom_literal_hits(seg):
                findings.append(("WARN", f"{where} (line {getattr(node,'lineno','?')}): "
                                         f"hard-typed figure-of-merit '{snip}' - interpolate it "
                                         f"from the computed value (f\"...{{value:.2f}}...\") so it "
                                         f"cannot desync from the code"))

    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr in _TEXT_SINKS):
            for a in list(node.args) + [kw.value for kw in node.keywords]:
                scan(a, f"{node.func.attr}(...)")
        elif isinstance(node, ast.Assign):
            names = [t.id.lower() for t in node.targets if isinstance(t, ast.Name)]
            if any(any(c in nm for c in _CAPTION_NAMES) for nm in names):
                scan(node.value, f"{'/'.join(names)} =")
        elif (isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)
              and node.value is not None
              and any(c in node.target.id.lower() for c in _CAPTION_NAMES)):
            scan(node.value, f"{node.target.id} =")

    if not findings:
        findings.append(("INFO", "no hard-typed figures-of-merit - on-figure numbers trace to the code"))
    return _critic_emit(findings, context)


def _critic_emit(findings, context):
    for s, m in findings:
        print(f"  [{s}] {context}: {m}")
    return findings


# ======================================================================
# from spectra.py
# ======================================================================
"""
py  -  FTIR / IR / ATR  and  PXRD / XRD  figures and the analysis behind them.

Conventions are flipped by cfg.domain (set it and forget it):
    FTIR : x = wavenumber cm-1, axis INVERTED, y = absorbance
    PXRD : x = 2theta degrees,  axis NORMAL,   y = intensity (counts / normalised)

Functions take cfg and reuse style/ The matplotlib is assumed-known; what
this encodes is the right axis directions, the identical-processing discipline,
and the overlaid-reps + waterfall idiom.

    load_xy(path)                       parse a 2-column xy/csv spectrum
    correct_baseline(x, y, cfg)         arPLS or rubberband -> (y_corr, baseline)
    normalize(x, y, cfg)                max / area / none
    integrate_bands(x, y, cfg)          local LINEAR baseline per window; %RSD-ready areas
    plot_overlay(traces, cfg)           reps overlaid on one axis
    plot_waterfall(groups, cfg)         vertically offset groups (the house idiom)
"""

# np.trapz was removed in NumPy 2.0 in favour of np.trapezoid
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


# --------------------------------------------------------------------- ingest
def load_xy(path, delimiter=None):
    """Parse a 2-column ascii spectrum (xy / csv / dat). Returns (x, y) sorted
    ascending by x. Comment lines (#, ;, letters) are skipped."""
    xs, ys = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#;%":
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0])); ys.append(float(parts[1]))
            except ValueError:
                continue
    x = np.asarray(xs); y = np.asarray(ys)
    order = np.argsort(x)
    return x[order], y[order]


# ------------------------------------------------------------------ baselines
def _arpls(y, lam=1e5, ratio=1e-6, niter=50):
    """Asymmetrically reweighted penalised least squares baseline (Baek 2015).
    No external dep required; uses a sparse second-difference penalty."""
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    y = np.asarray(y, float)
    L = len(y)
    D = sparse.diags([1., -2., 1.], [0, -1, -2], shape=(L, L - 2))
    H = (lam * (D @ D.T)).tocsc()
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.diags(w, 0, shape=(L, L))
        z = spsolve((W + H).tocsc(), w * y)
        d = y - z
        dn = d[d < 0]
        if dn.size == 0:
            break
        m, s = np.mean(dn), np.std(dn)
        arg = np.clip(2 * (d - (2 * s - m)) / (s + 1e-12), -500, 500)  # avoid exp overflow
        wt = 1.0 / (1.0 + np.exp(arg))
        if np.linalg.norm(w - wt) / (np.linalg.norm(w) + 1e-12) < ratio:
            w = wt; break
        w = wt
    return z


def _rubberband(x, y):
    """Convex-hull 'rubber band' baseline - robust, dependency-free fallback."""
    try:
        from scipy.spatial import ConvexHull
        pts = np.column_stack([x, y])
        v = ConvexHull(pts).vertices
        v = np.roll(v, -v.argmin())
        v = v[: v.argmax() + 1]              # lower hull
        return np.interp(x, x[v], y[v])
    except Exception:
        return np.minimum.accumulate(y)      # crude monotone floor


def correct_baseline(x, y, cfg):
    """Return (y_corrected, baseline). Honours cfg.baseline; degrades from arPLS
    to rubberband if scipy is unavailable."""
    method = cfg.baseline
    if method is None:
        return np.asarray(y, float), np.zeros_like(y, float)
    if method == "arpls":
        try:
            b = _arpls(np.asarray(y, float), lam=cfg.arpls_lam)
        except Exception:
            b = _rubberband(np.asarray(x, float), np.asarray(y, float))
    elif method == "rubberband":
        b = _rubberband(np.asarray(x, float), np.asarray(y, float))
    else:
        raise ValueError(f"unknown baseline '{method}'")
    return np.asarray(y, float) - b, b


def normalize(x, y, cfg):
    y = np.asarray(y, float)
    if cfg.normalize in (None, "none"):
        return y
    if cfg.normalize == "max":
        m = np.max(np.abs(y))
        return y / m if m else y
    if cfg.normalize == "area":
        a = _trapz(np.abs(y), x)
        return y / a if a else y
    raise ValueError(f"unknown normalize '{cfg.normalize}'")


def apply_smooth(x, y, cfg):
    """Optional DISPLAY smoothing (Savitzky-Golay). Smooth for the figure only;
    always integrate/quantify on the RAW (baseline-corrected) trace, never on a
    smoothed one - smoothing changes peak areas. Returns y unchanged if
    cfg.smooth is None or scipy is unavailable."""
    if cfg.smooth in (None, "none"):
        return np.asarray(y, float)
    if cfg.smooth == "savgol":
        try:
            from scipy.signal import savgol_filter
            w = cfg.savgol_window if cfg.savgol_window % 2 == 1 else cfg.savgol_window + 1
            w = min(w, len(y) - (1 - len(y) % 2))       # window can't exceed length
            if w <= cfg.savgol_poly:
                return np.asarray(y, float)
            return savgol_filter(np.asarray(y, float), w, cfg.savgol_poly)
        except Exception:
            return np.asarray(y, float)
    raise ValueError(f"unknown smooth '{cfg.smooth}'")


def _sg(y, smooth, deriv=0, delta=1.0):
    """Savitzky-Golay smooth (deriv=0) or derivative of order `deriv`, on a COPY.
    `smooth` is (window, poly); window is forced odd and clipped to the series
    length. Degrades to the raw trace (deriv=0) or a finite-difference derivative
    if scipy is unavailable or the window is unworkable - never raises."""
    y = np.asarray(y, float)
    if not smooth:
        if deriv == 0:
            return y
        d = y
        for _ in range(deriv):
            d = np.gradient(d, delta)
        return d
    win, poly = smooth
    try:
        from scipy.signal import savgol_filter
        w = win if win % 2 == 1 else win + 1
        w = min(w, len(y) - (1 - len(y) % 2))          # window can't exceed length
        if w <= poly:
            raise ValueError
        return savgol_filter(y, w, poly, deriv=deriv, delta=delta)
    except Exception:
        if deriv == 0:
            return y
        d = y
        for _ in range(deriv):
            d = np.gradient(d, delta)
        return d


def _anchor_val(x, y, w, half=2.0):
    """Noise-robust absorbance AT an anchor: mean over w +/- half cm-1 (falls back
    to the single nearest point if the window is empty)."""
    m = (x >= w - half) & (x <= w + half)
    return float(y[m].mean()) if m.any() else float(y[np.argmin(np.abs(x - w))])


# ------------------------------------------------------------ automatic anchors
def find_anchors(x, y, center, gap=12.0, maxhw=60.0, smooth=(13, 3)):
    """Locate the flanking local minima (the 'continuum-return' points) either side
    of a band centred near `center`, and return them as (lo, hi) anchor wavenumbers.

    WHY: a local linear baseline needs two anchors; picking them by eye is the
    irreproducible step that makes one analyst's %RSD differ from another's. The
    valley on each side of a band IS the chemically-defensible anchor, and it can
    be found deterministically: the argmin of a lightly smoothed trace within a
    bounded side-window. For a band sitting next to a neighbour (e.g. an ester C=O
    on the wing of an acid C=O) the low-side minimum is the inter-band valley -
    exactly where a careful analyst would put it.

    The search runs on a Savitzky-Golay-smoothed COPY (location only - you still
    integrate on the raw trace). Each side is searched in [center +/- gap,
    center +/- maxhw]: `gap` keeps the search off the peak itself, `maxhw` bounds
    it so a distant band can't steal the anchor.

    DISCIPLINE: detect ONCE on a fixed reference spectrum (the pure-analyte trace,
    or a batch mean) and LOCK the result, then apply the same two anchors to every
    spectrum in the batch. Re-detecting per replicate re-introduces a small,
    avoidable variability and violates identical-processing (see `auto_windows`)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    o = np.argsort(x); x, y = x[o], y[o]
    ys = _sg(y, smooth)
    def side_min(a, b):
        m = (x >= min(a, b)) & (x <= max(a, b))
        if not m.any():
            raise ValueError(f"find_anchors: no points in search window [{a:.1f},{b:.1f}] cm-1")
        xs = x[m]
        return float(xs[int(np.argmin(ys[m]))])
    lo = side_min(center - maxhw, center - gap)
    hi = side_min(center + gap, center + maxhw)
    return (lo, hi)


def auto_windows(ref_x, ref_y, cfg):
    """Detect-once-and-lock helper. Runs find_anchors on a REFERENCE spectrum for
    every (center, name) in cfg.anchor_centers and returns a list of
    (lo, hi, name) tuples in the same shape as cfg.integration_windows. Assign the
    result back onto cfg ONCE, then process the whole batch with those locked
    anchors:

        cfg.integration_windows = auto_windows(ref_x, ref_y, cfg)

    This keeps the algorithmic, justified anchor choice while preserving the
    identical-processing guarantee that makes a rep-to-rep %RSD meaningful."""
    out = []
    for center, name in cfg.anchor_centers:
        lo, hi = find_anchors(ref_x, ref_y, center,
                              gap=cfg.anchor_gap, maxhw=cfg.anchor_maxhw,
                              smooth=cfg.anchor_smooth)
        out.append((lo, hi, name))
    return out


def band_metric(x, y, metric="area", anchors=None, center=None, window=None,
                gap=12.0, maxhw=60.0, anchor_smooth=(13, 3), deriv_smooth=(13, 3),
                anchor_half=2.0):
    """One scalar for one band, in three algorithmic flavours that trade off how
    much they depend on a baseline:

      "area"   - trapezoidal area above a local LINEAR baseline between `anchors`
                 (auto-found from `center` when anchors is None). Units: abs.cm-1.
      "height" - peak height above that same linear baseline, taken inside
                 `window` (defaults to the anchor span).
      "deriv2" - depth of the 2nd-derivative trough inside `window`, i.e.
                 max(-d2y/dx2). ANCHOR-FREE: a constant + linear background has a
                 zero 2nd derivative, so a sloping continuum (a neighbouring band's
                 wing) cancels and no baseline subtraction is needed. Pays for it
                 in noise -> always computed on an SG-smoothed derivative.

    x, y need not be sorted. Returns a float (NaN if the band window is empty)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    o = np.argsort(x); x, y = x[o], y[o]

    if metric == "deriv2":
        if window is None:
            if center is None:
                raise ValueError("band_metric(deriv2) needs `window` or `center`")
            window = (center - 8.0, center + 8.0)
        lo, hi = sorted(window)
        dx = float(np.mean(np.diff(x)))
        d2 = _sg(y, deriv_smooth, deriv=2, delta=dx)
        m = (x >= lo) & (x <= hi)
        return float((-d2[m]).max()) if m.any() else float("nan")

    if metric not in ("area", "height"):
        raise ValueError(f"unknown metric '{metric}'")
    if anchors is None:
        if center is None:
            raise ValueError(f"band_metric({metric}) needs `anchors` or `center`")
        anchors = find_anchors(x, y, center, gap=gap, maxhw=maxhw, smooth=anchor_smooth)
    a_lo, a_hi = sorted(anchors)
    m = (x >= a_lo) & (x <= a_hi)
    if m.sum() < 2:
        return float("nan")
    xs, ys = x[m], y[m]
    base = np.interp(xs, [a_lo, a_hi],
                     [_anchor_val(x, y, a_lo, anchor_half), _anchor_val(x, y, a_hi, anchor_half)])
    corr = ys - base
    if metric == "area":
        return float(_trapz(np.clip(corr, 0, None), xs))
    w = window if window else (a_lo, a_hi)
    wl, wh = sorted(w); mm = (xs >= wl) & (xs <= wh)
    return float(corr[mm].max()) if mm.any() else float("nan")


def as_transmittance(absorbance):
    """Convert absorbance to %T (T = 100 * 10**-A). FTIR is conventionally shown
    either way; pick one per figure via cfg.ftir_yaxis and stay consistent."""
    return 100.0 * np.power(10.0, -np.asarray(absorbance, float))


def label_peaks(ax, peaks, cfg, dy=0.02):
    """Annotate peaks the CALLER supplies - this styles labels consistently, it
    does NOT detect peaks or invent text (that judgement is yours). `peaks` is a
    list of (x, y, text). Labels are placed just above each point."""
    for px, py, text in peaks:
        ax.annotate(text, xy=(px, py), xytext=(px, py + dy),
                    ha="center", va="bottom", fontsize="small",
                    rotation=90 if cfg.domain == "ftir" else 0)
    return ax


def edge_labels(ax, items, dx=0.015, fontsize="small", fontweight="bold"):
    """The house waterfall KEY: direct labels at the RIGHT axis margin, each aligned
    to its own trace's row -> grayscale-safe (no colour-only legend box) and never
    floating over the data. `items` is an iterable of (y_data, text, colour). x is in
    AXES fraction just past the right spine (so it survives x-axis inversion); y is in
    DATA coords. Returns the placed Text objects so a composite can verify they don't
    collide with a neighbouring panel (audit_layout flags that). Used by
    plot_waterfall, and callable directly for hand-built waterfalls / composite panels."""
    placed = []
    for y, text, color in items:
        t = ax.text(1.0 + dx, y, text, transform=ax.get_yaxis_transform(),
                    va="center", ha="left", color=color, clip_on=False,
                    fontsize=fontsize, fontweight=fontweight)
        placed.append(t)
    return placed


# ----------------------------------------------------------------- integration
def integrate_bands(x, y, cfg):
    """For each (lo, hi, name) window in cfg.integration_windows, integrate the
    band above a linear baseline. The same anchors are used for every spectrum
    (caller should pass identical cfg) so rep-to-rep area differences are real,
    not a windowing artefact. Returns list of dicts: {name, area, lo, hi}.

    cfg.integration_baseline picks how the baseline is drawn:
      "per_window" (default): a separate local baseline per window, between that
          window's own two edges. Correct for ISOLATED bands.
      "shared": ONE baseline across the whole envelope (leftmost..rightmost anchor
          of all windows); each window integrates its sub-interval above that single
          line - a vertical-drop / "drop perpendicular" partition at the boundaries.
          Use for OVERLAPPING bands on a shared pedestal (e.g. analyte band next to
          an internal-standard band): a per-window baseline would follow the valley
          up and carve area off the smaller band.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    windows = cfg.integration_windows
    mode = getattr(cfg, "integration_baseline", "per_window")

    # shared mode draws the baseline once, across the union of all windows
    env = None
    if mode == "shared" and windows:
        env_lo = min(min(lo, hi) for lo, hi, _ in windows)
        env_hi = max(max(lo, hi) for lo, hi, _ in windows)
        em = (x >= env_lo) & (x <= env_hi)
        if em.sum() >= 2:
            xe = x[em]
            env = (env_lo, env_hi, y[em][0], y[em][-1])   # baseline through the envelope edges
    elif mode not in ("per_window", "shared"):
        raise ValueError(f"unknown integration_baseline '{mode}'")

    out = []
    for lo, hi, name in windows:
        a, b = sorted((lo, hi))
        m = (x >= a) & (x <= b)
        if m.sum() < 2:
            out.append({"name": name, "area": float("nan"), "lo": a, "hi": b})
            continue
        xs, ys = x[m], y[m]
        if env is not None:
            elo, ehi, eylo, eyhi = env
            base = np.interp(xs, [elo, ehi], [eylo, eyhi])   # the one shared line
        else:
            base = np.interp(xs, [xs[0], xs[-1]], [ys[0], ys[-1]])   # local per-window line
        area = float(_trapz(np.clip(ys - base, 0, None), xs))
        out.append({"name": name, "area": area, "lo": a, "hi": b})
    return out


# --------------------------------------------------------------------- plots
def _apply_axis(ax, cfg):
    if cfg.domain == "ftir":
        ax.set_xlabel(r"Wavenumber / cm$^{-1}$")
        if cfg.ftir_yaxis == "transmittance":
            ax.set_ylabel("Transmittance / %")
        else:
            ax.set_ylabel("Absorbance" if cfg.normalize in (None, "none") else "Absorbance (norm.)")
        lo, hi = (cfg.x_limits if cfg.x_limits else ax.get_xlim())
        ax.set_xlim(max(lo, hi), min(lo, hi))     # INVERTED
    elif cfg.domain == "pxrd":
        ax.set_xlabel(r"2$\theta$ / degree")
        ax.set_ylabel("Intensity" if cfg.normalize in (None, "none") else "Intensity (norm.)")
        if cfg.x_limits:
            ax.set_xlim(min(cfg.x_limits), max(cfg.x_limits))   # NORMAL
    else:
        raise ValueError(f"unknown domain '{cfg.domain}'")


def plot_overlay(traces, cfg, ax=None):
    """traces: dict label -> (x, y). Reps overlaid on a single axis."""
    apply_style(cfg)
    if ax is None:
        fig, ax = figure(cfg)
    else:
        fig = ax.figure
    for label, (x, y) in traces.items():
        check_trace(x, y, cfg, name=label)
        ax.plot(x, y, label=label)
    _apply_axis(ax, cfg)
    if len(traces) > 1:
        ax.legend(loc="best")
    return fig, ax


def plot_waterfall(groups, cfg, offset=None):
    """groups: dict label -> list of (x, y) reps. Each group is offset vertically
    (the house waterfall idiom); reps within a group are overlaid in one colour.
    PXRD reference stick lines from cfg.pxrd_reference_lines are drawn at the base."""
    apply_style(cfg)
    fig, ax = figure(cfg, height=None)
    # estimate a sensible offset from the data span if not given
    spans = [np.ptp(y) for trs in groups.values() for _, y in trs if len(y)]
    step = offset if offset is not None else (1.15 * max(spans) if spans else 1.0)
    pal = _palette(cfg)
    edges = []                                  # (y_at_right_edge, label, colour) per group
    for gi, (label, trs) in enumerate(groups.items()):
        color = pal[gi % len(pal)]
        for ri, (x, y) in enumerate(trs):
            check_trace(x, y, cfg, name=f"{label}#{ri}")
            # waterfalls distinguish groups by vertical POSITION -> solid lines
            # always (a dash period would fight a choppy fingerprint region)
            ax.plot(x, np.asarray(y, float) + gi * step,
                    color=color, ls="-", label=label if ri == 0 else None)
        x0, y0 = np.asarray(trs[0][0], float), np.asarray(trs[0][1], float)
        # the plot's RIGHT edge is min-x for ftir (inverted) and max-x for pxrd
        idx = int(np.argmin(x0)) if cfg.domain == "ftir" else int(np.argmax(x0))
        edges.append((y0[idx] + gi * step, label, color))
    if cfg.domain == "pxrd":
        for tt, lab in cfg.pxrd_reference_lines:
            ax.axvline(tt, color="0.6", lw=0.5, ls=":")
    _apply_axis(ax, cfg)
    ax.set_yticks([])                          # offsets are arbitrary; hide the y scale
    mode = getattr(cfg, "waterfall_legend", "edge")
    if len(groups) > 1 and mode == "legend":
        ax.legend(loc="best")
    elif len(groups) > 1 and mode == "edge":
        edge_labels(ax, edges)            # the house key: right-margin, per-trace
    return fig, ax


# ======================================================================
# from calibration.py
# ======================================================================
"""
py  -  linear calibration done honestly, ICH Q2 flavour.

The correctness traps this guards against:
  * R^2 alone is NOT evidence of linearity -> a residual panel is MANDATORY.
  * Reading concentrations off the curve outside the calibrated range
    (extrapolation) -> predict() refuses it.
  * Quoting LOD/LOQ without saying how -> the method is named in the result.

    fit(x, y, cfg)                  -> dict(slope, intercept, r2, s_resid, n, ...)
    predict(model, signal, cfg)     -> concentration (raises on extrapolation)
    lod_loq(model, cfg)             -> dict(lod, loq, method)
    plot_calibration(x, y, cfg)     -> figure with curve + bands + MANDATORY residuals
"""


def _tmult(df, conf):
    try:
        from scipy import stats
        return float(stats.t.ppf(0.5 + conf / 2, df))
    except Exception:
        return {0.90: 1.645, 0.95: 1.960, 0.99: 2.576}.get(round(conf, 2), 1.960)


def fit(x, y, cfg):
    """Ordinary least squares y = slope*x + intercept, with the statistics needed
    for honest bands and a slope CI."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    check_trace(np.sort(x), y[np.argsort(x)], cfg, name="calibration",
                       min_points=3, require_monotonic=False)
    n = x.size
    xbar = x.mean()
    Sxx = float(np.sum((x - xbar) ** 2))
    slope, intercept = np.polyfit(x, y, 1)
    yhat = slope * x + intercept
    resid = y - yhat
    dof = n - 2
    s_resid = float(np.sqrt(np.sum(resid ** 2) / dof)) if dof > 0 else float("nan")
    se_slope = s_resid / np.sqrt(Sxx) if Sxx else float("nan")
    se_intercept = s_resid * np.sqrt(1.0 / n + xbar ** 2 / Sxx) if Sxx else float("nan")
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - np.sum(resid ** 2) / ss_tot if ss_tot else float("nan")
    t = _tmult(dof, cfg.conf_level)
    return {
        "slope": float(slope), "intercept": float(intercept), "r2": float(r2),
        "s_resid": s_resid, "se_slope": se_slope, "se_intercept": se_intercept,
        "n": int(n), "dof": int(dof),
        "xbar": xbar, "Sxx": Sxx, "t": t, "conf": cfg.conf_level,
        "x_range": (float(x.min()), float(x.max())),
        "x": x, "y": y, "resid": resid,
        "slope_ci": (float(slope - t * se_slope), float(slope + t * se_slope)),
        "intercept_ci": (float(intercept - t * se_intercept), float(intercept + t * se_intercept)),
    }


def predict(model, signal, cfg, x_for_signal=None):
    """Invert the curve: concentration from a measured signal. Refuses to read
    outside the calibrated x range when cfg.refuse_extrapolation (the default)."""
    conc = (signal - model["intercept"]) / model["slope"]
    lo, hi = model["x_range"]
    if cfg.refuse_extrapolation and not (lo <= conc <= hi):
        raise GateError(
            f"predicted x={conc:.4g} is outside the calibrated range [{lo:.4g}, {hi:.4g}] "
            f"- extrapolation refused (extend the calibration instead)")
    return float(conc)


def lod_loq(model, cfg):
    """LOD = 3.3*sigma/slope, LOQ = 10*sigma/slope. With 'residual_sd' (default)
    sigma is the regression residual SD; the method is recorded in the result."""
    if cfg.lod_loq_method == "residual_sd":
        sigma = model["s_resid"]
    else:
        raise ValueError("blank_sd method needs blank replicates passed separately")
    m = abs(model["slope"])
    return {"lod": 3.3 * sigma / m, "loq": 10.0 * sigma / m,
            "method": f"{cfg.lod_loq_method} (3.3σ/m, 10σ/m)"}


def intercept_test(model, cfg):
    """Test whether the intercept differs significantly from zero — i.e. whether the
    method carries a CONSTANT BIAS. Two-sided t-test on b/SE(b) at (n-2) dof.

    Returns {intercept, se, t, p, ci, significant, caption}. A CI that INCLUDES 0
    (equivalently p > 1-conf) means no significant constant bias — the line is
    consistent with passing through the origin. A significant intercept flags a
    constant offset (baseline, matrix, or low-end lack-of-fit) to investigate."""
    b, se, dof = model["intercept"], model["se_intercept"], model["dof"]
    t_stat = b / se if se else float("nan")
    try:
        from scipy import stats
        p = float(2 * stats.t.sf(abs(t_stat), dof))
    except Exception:
        p = float("nan")
    lo, hi = model["intercept_ci"]
    significant = not (lo <= 0.0 <= hi)
    verdict = ("constant bias: intercept significantly != 0 (CI excludes 0)" if significant
               else "no significant constant bias (CI includes 0)")
    cap = (f"intercept {b:.4g} [{lo:.4g}, {hi:.4g}] ({int(model['conf']*100)}% CI), "
           f"t={t_stat:.2f}, p={p:.3f} -> {verdict}")
    return {"intercept": float(b), "se": float(se), "t": float(t_stat), "p": p,
            "ci": (lo, hi), "significant": bool(significant), "caption": cap}


def plot_calibration(x, y, cfg, model=None):
    """Calibration curve with confidence + prediction bands AND a mandatory
    residual panel underneath (R^2 alone does not certify linearity)."""
    if model is None:
        model = fit(x, y, cfg)
    apply_style(cfg)

    if cfg.force_residual_panel:
        fig, (ax, axr) = figure(cfg, nrows=2, height=None,
                                      height_ratios=[3, 1], sharex=True)
    else:
        fig, ax = figure(cfg); axr = None

    xs = np.linspace(*model["x_range"], 200)
    ys = model["slope"] * xs + model["intercept"]
    # CI of the mean response, and wider prediction band
    s, t, Sxx, xbar, n = (model["s_resid"], model["t"], model["Sxx"],
                          model["xbar"], model["n"])
    se_mean = s * np.sqrt(1.0 / n + (xs - xbar) ** 2 / Sxx)
    se_pred = s * np.sqrt(1.0 + 1.0 / n + (xs - xbar) ** 2 / Sxx)

    ax.scatter(model["x"], model["y"], zorder=3, label="data")
    ax.plot(xs, ys, label="fit")
    ax.fill_between(xs, ys - t * se_mean, ys + t * se_mean, alpha=0.25, lw=0,
                    label=f"{int(model['conf']*100)}% CI")
    ax.plot(xs, ys - t * se_pred, lw=0.6, ls="--", color="0.5")
    ax.plot(xs, ys + t * se_pred, lw=0.6, ls="--", color="0.5",
            label=f"{int(model['conf']*100)}% pred.")
    ax.set_ylabel("Response")
    ax.legend(loc="best")
    # report R^2 and slope CI in a corner annotation, not as a substitute for residuals
    lo, hi = model["slope_ci"]
    ax.annotate(f"$R^2$={model['r2']:.4f}\nslope {model['slope']:.3g} "
                f"[{lo:.3g}, {hi:.3g}]",
                xy=(0.03, 0.97), xycoords="axes fraction", va="top", fontsize="small")

    if axr is not None:
        axr.axhline(0, color="0.6", lw=0.6)
        axr.scatter(model["x"], model["resid"], zorder=3)
        axr.set_ylabel("Resid. (y units)")
        axr.set_xlabel("Concentration")
    else:
        ax.set_xlabel("Concentration")
    finalize_figure(fig)
    return fig, (ax, axr)


# ======================================================================
# from chemometrics.py
# ======================================================================
"""
py  -  multivariate calibration (PLS / PCR / PCA) for spectra, honestly.

This is the FOURTH family (after spectra, calibration, crystal). Univariate
`py` regresses ONE band area on concentration; this regresses the WHOLE
spectrum, which is what you reach for when no single band is selective (overlap,
matrix effects, scatter/loading variation). The figures are the standard PLS read-out:
RMSECV-vs-components, the regression-coefficient spectrum, scores, and loadings.

The correctness traps this guards against (the reason it lives in THIS skill and not
a notebook):

  * CV LEAKAGE. Preprocessing that uses cross-sample statistics - MSC's reference
    spectrum, mean-centring, autoscaling - must be fit on the TRAINING rows of each
    fold only. Fit it once on the whole set and the test rows have leaked into
    training, so RMSECV is optimistic. `cross_validate` re-fits a FRESH `Preprocessor`
    inside every fold. (SNV and derivatives are per-sample, so they don't leak - but
    the machinery is uniform so you can't get it wrong by accident.)

  * OVER-FITTING THE COMPONENT COUNT. The global-minimum RMSECV almost always sits at
    too many latent variables. `choose_n_components` applies a PARSIMONY rule: the
    fewest components whose RMSECV is within `cfg.lv_parsimony_tol` of the global min.

  * THE OPTIMISTIC-CV TRAP. Leave-one-out over replicate spectra that share a
    concentration LEVEL answers "predict a new sample at a level I've already seen",
    not "predict a new level". The two give different errors. Every CV result NAMES
    which question its scheme answers, and `cross_validate` WARNs if you LOO over data
    whose y replicates (telling you to pass `groups=` for leave-one-level-out).

scikit-learn is a LAZY dependency: it is imported only when a model is actually fit,
mirroring the crystal family's heavy deps. An FTIR or univariate-calibration job never
imports it.

    Preprocessor(steps, cfg)              leakage-aware SNV / MSC / d1 / d2 / centre / autoscale
    pls_fit(X, y, nc, cfg, pre=...)       -> model dict (coef, scores, loadings, weights, yhat...)
    pcr_fit / pca_fit                     PCR regression / PCA decomposition (same dict shape)
    cross_validate(X, y, nc, cfg, ...)    leakage-safe RMSECV + held-out predictions, scheme NAMED
    component_scan(X, y, cfg, ...)        RMSECV vs components, one curve per preprocessing
    choose_n_components(scan, cfg)        the parsimony pick
    vip(model)                            VIP scores (which wavenumbers drive the model)
    plot_rmsecv / plot_coefficients / plot_scores / plot_loadings   ax-level panels
    diagnostics_figure(X, y, x, cfg, ...) the gated + QA'd 2x2 composite deliverable
    methods_text(model, cv, choice, cfg)  paste-ready methods paragraph

(methods_text, not methods_report — the latter is py's name, and the bundler
flattens every module into one namespace, so module-level names must be unique.)
"""

_SKLEARN_HINT = ("Required package not found: scikit-learn. Install with:\n"
                 "    python -m pip install scikit-learn\n"
                 "(the chemometrics family needs sklearn for PLS / PCA / PCR).")


def _pls_cls():
    try:
        from sklearn.cross_decomposition import PLSRegression
        return PLSRegression
    except ImportError as e:
        raise ImportError(_SKLEARN_HINT) from e


def _pca_cls():
    try:
        from sklearn.decomposition import PCA
        return PCA
    except ImportError as e:
        raise ImportError(_SKLEARN_HINT) from e


def _linreg_cls():
    try:
        from sklearn.linear_model import LinearRegression
        return LinearRegression
    except ImportError as e:
        raise ImportError(_SKLEARN_HINT) from e


# ===================================================== leakage-aware preprocessing
_STATELESS = {"none", "snv", "d1", "d2"}
_STATEFUL = {"msc", "center", "autoscale"}


class Preprocessor:
    """A small preprocessing pipeline whose STATEFUL steps are fit on training rows
    only, so it is safe to drive cross-validation with (the whole reason it exists).

    steps: a list like ["snv"] or ["d1", "center"], or a "+"-joined string
    ("snv+center"). Order is applied left to right.
      stateless (identical train/test, no leakage):
        "none"  pass-through
        "snv"   standard normal variate - per-row centre & scale (scatter/path-length)
        "d1"    Savitzky-Golay 1st derivative (window from cfg.deriv_smooth)
        "d2"    Savitzky-Golay 2nd derivative (resolves overlap; amplifies noise)
      stateful (parameters estimated from TRAINING rows -> must go through fit):
        "msc"       multiplicative scatter correction against the training-mean spectrum
        "center"    subtract the training per-wavenumber mean
        "autoscale" centre and divide by the training per-wavenumber SD

    Use fit_transform on a training block and transform on the matching test block;
    cross_validate / component_scan do this for every fold automatically.

    Note on the scale convention: SNV and `autoscale` use the POPULATION standard
    deviation (ddof=0), matching sklearn's StandardScaler and the `chemotools` library.
    R/prospectr use the sample SD (ddof=1) -- a constant sqrt(n/(n-1)) factor. The ddof=0
    choice is intentional and is locked by the chemometrics validation suite
    (skill_validation/chemometrics/)."""

    def __init__(self, steps, cfg=None):
        if isinstance(steps, str):
            steps = [s.strip() for s in steps.split("+") if s.strip()]
        steps = list(steps) if steps is not None and len(steps) else ["none"]
        bad = [s for s in steps if s not in (_STATELESS | _STATEFUL)]
        if bad:
            raise ValueError(f"unknown preprocessing step(s): {bad} "
                             f"(allowed: {sorted(_STATELESS | _STATEFUL)})")
        self.steps = steps
        self.cfg = cfg
        self._state = {}        # per-step-index state for the stateful steps

    @property
    def name(self):
        return "+".join(self.steps)

    def _deriv(self, X, order):
        from scipy.signal import savgol_filter
        win, poly = (self.cfg.deriv_smooth if self.cfg is not None else (17, 3))
        w = win if win % 2 == 1 else win + 1
        w = min(w, X.shape[1] - (1 - X.shape[1] % 2))   # window can't exceed n_features
        if w <= poly:
            import warnings
            warnings.warn(
                f"Savitzky-Golay d{order} skipped: usable window ({w}) <= polyorder "
                f"({poly}) for {X.shape[1]} features; returning the input UNCHANGED "
                f"(no derivative applied).", stacklevel=2)
            return X
        return savgol_filter(X, w, poly, deriv=order, axis=1)

    def _run(self, X, fit):
        X = np.asarray(X, float)
        for i, s in enumerate(self.steps):
            st = self._state.get(i, {})
            if s == "none":
                pass
            elif s == "snv":
                mu = X.mean(axis=1, keepdims=True)
                sd = X.std(axis=1, keepdims=True)
                X = (X - mu) / np.where(sd == 0, 1.0, sd)
            elif s == "d1":
                X = self._deriv(X, 1)
            elif s == "d2":
                X = self._deriv(X, 2)
            elif s == "msc":
                if fit:
                    st["ref"] = X.mean(axis=0)
                ref = st["ref"]
                out = np.empty_like(X)
                for r in range(X.shape[0]):
                    b1, b0 = np.polyfit(ref, X[r], 1)        # X[r] ~ b1*ref + b0
                    out[r] = (X[r] - b0) / (b1 if b1 != 0 else 1.0)
                X = out
            elif s == "center":
                if fit:
                    st["mean"] = X.mean(axis=0)
                X = X - st["mean"]
            elif s == "autoscale":
                if fit:
                    st["mean"] = X.mean(axis=0)
                    sd = X.std(axis=0)
                    st["std"] = np.where(sd == 0, 1.0, sd)
                X = (X - st["mean"]) / st["std"]
            if fit:
                self._state[i] = st
        return X

    def fit(self, X):
        self._run(X, fit=True)
        return self

    def transform(self, X):
        return self._run(X, fit=False)

    def fit_transform(self, X):
        return self._run(X, fit=True)


# ===================================================== EMSC
def emsc(X, wavenumbers, reference=None, poly_order=2, interferents=None,
         weights=None, return_model=False):
    """Extended Multiplicative Signal Correction (Martens & Stark 1991; Afseth & Kohler 2012).

    Fits each spectrum    x ~= b*ref + sum_k a_k * z^k  (+ sum_j d_j * interferent_j),
    where z = wavenumbers affine-scaled to [-1, 1], then returns the corrected spectrum
        x_corr = (x - baseline - interferents) / b.
    `reference` defaults to the MEAN spectrum (the EMSC convention). Validated bit-for-bit
    (~5e-16) against the biospectools / Kohler MATLAB gold standard.

    The principled upgrade over SNV/MSC when the problem is multiplicative loading plus a
    baseline. `interferents` (rows of an array, e.g. a pure-excipient spectrum) models known
    constituents and subtracts them — the project's "model the matrix away" move.

    Caveat on b: with reference=mean AND a baseline term, b is NOT a physical dilution factor
    (the constant polynomial column absorbs part of the offset). To recover a known multiplicative
    factor, use a PURE-component reference with poly_order=0.

    Leakage note: a FIXED (pure/external) reference + interferents is leak-free — apply it before
    CV. The mean-reference variant is cross-sample; refit it per fold if used inside cross_validate.

    Returns the corrected array (same leading dim as X); also the coef/model dict if return_model.
    """
    Xin = np.asarray(X, float)
    X2 = np.atleast_2d(Xin)
    wn = np.asarray(wavenumbers, float)
    ref = X2.mean(0) if reference is None else np.asarray(reference, float)
    rng = wn.max() - wn.min()
    z = 2.0 * (wn - wn.min()) / (rng if rng else 1.0) - 1.0
    polys = np.vstack([z ** k for k in range(poly_order + 1)]).T        # (nwn, p+1)
    blocks = [ref[:, None], polys]
    I = None
    if interferents is not None:
        I = np.atleast_2d(np.asarray(interferents, float))             # (n_interf, nwn)
        blocks.append(I.T)
    M = np.hstack(blocks)                                              # [ ref | polys | interferents ]
    npoly = poly_order + 1
    w = None if weights is None else np.sqrt(np.asarray(weights, float))
    Mfit = M if w is None else M * w[:, None]

    corrected = np.empty_like(X2)
    coefs = np.empty((X2.shape[0], M.shape[1]))
    for i in range(X2.shape[0]):
        target = X2[i] if w is None else X2[i] * w
        c, *_ = np.linalg.lstsq(Mfit, target, rcond=None)
        b = c[0]
        baseline = polys @ c[1:1 + npoly]
        interf = (I.T @ c[1 + npoly:]) if I is not None else 0.0
        corrected[i] = (X2[i] - baseline - interf) / (b if b != 0 else 1.0)
        coefs[i] = c
    out = corrected if Xin.ndim == 2 else corrected[0]
    if return_model:
        return out, {"coefs": coefs, "ref": ref, "M": M, "poly_order": poly_order,
                     "n_interferents": 0 if I is None else I.shape[0],
                     "residual": X2 - coefs @ M.T}
    return out


# ===================================================== input gate
def _check_xy(X, y, cfg, n_components=None):
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    f = []
    if X.ndim != 2:
        f.append(("FAIL", f"X must be 2D (n_samples, n_features); got shape {X.shape}"))
        return _resolve(f, cfg, "chemometrics"), X, y
    if X.shape[0] != y.size:
        f.append(("FAIL", f"X has {X.shape[0]} rows but y has {y.size}"))
    if not np.isfinite(X).all():
        f.append(("FAIL", "X has non-finite values (NaN/inf) - check the ingest/baseline"))
    if not np.isfinite(y).all():
        f.append(("FAIL", "y has non-finite values"))
    if X.shape[0] < 4:
        f.append(("WARN", f"only {X.shape[0]} samples - CV and component estimates are unstable"))
    if n_components is not None and n_components > X.shape[0] - 1:
        f.append(("WARN", f"n_components={n_components} >= n_samples ({X.shape[0]}); over-fit risk"))
    if not f:
        f.append(("INFO", f"X {X.shape}, y n={y.size} - OK"))
    _resolve(f, cfg, "chemometrics")
    return f, X, y


# ===================================================== model fits (on the FULL set)
def _coef_vector(model, n_features):
    coef = np.ravel(np.asarray(model.coef_))
    if coef.size != n_features:                 # sklearn version differences in orientation
        coef = np.ravel(np.asarray(model.coef_).T)
    return coef


def _fit_stats(y, yhat):
    resid = y - yhat
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - float(np.sum(resid ** 2)) / ss_tot if ss_tot else float("nan")
    rmse = float(np.sqrt(np.mean(resid ** 2)))
    return resid, float(r2), rmse


def pls_fit(X, y, n_components, cfg, pre=None):
    """Fit PLS regression on the WHOLE set. Use this for the interpretable read-out
    (coefficients, scores, loadings) and the calibration fit; use cross_validate for
    the HONEST error. Returns a model dict."""
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, y = _check_xy(X, y, cfg, n_components)
    pp = Preprocessor(pre, cfg)
    Xp = pp.fit_transform(X)
    nc = min(n_components, X.shape[0] - 1, Xp.shape[1])
    PLS = _pls_cls()
    m = PLS(n_components=nc, scale=getattr(cfg, "pls_scale", True)).fit(Xp, y)
    coef = _coef_vector(m, Xp.shape[1])
    yhat = np.ravel(m.predict(Xp))
    resid, r2, rmsec = _fit_stats(y, yhat)
    return {"kind": "pls", "model": m, "pre": pp, "pre_name": pp.name,
            "n_components": int(nc), "coef": coef,
            "scores": np.asarray(m.x_scores_), "loadings": np.asarray(m.x_loadings_),
            "weights": np.asarray(m.x_weights_), "y_loadings": np.ravel(m.y_loadings_),
            "yhat": yhat, "resid": resid, "r2_cal": r2, "rmsec": rmsec, "y": y,
            "caption": f"PLS ({pp.name}, {nc} LV): R2(cal)={r2:.3f}, RMSEC={rmsec:.3g}"}


def pca_fit(X, cfg, n_components=None, pre=None):
    """PCA decomposition for exploratory scores/loadings (no y). sklearn PCA centres
    internally, so `pre` need not include 'center'."""
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, _y = _check_xy(X, np.zeros(X.shape[0]), cfg)
    pp = Preprocessor(pre, cfg)
    Xp = pp.fit_transform(X)
    nc = n_components or min(getattr(cfg, "pls_max_components", 10), *Xp.shape)
    nc = int(min(nc, *Xp.shape))            # cap at min(n_samples, n_features) (sklearn's own limit)
    PCA = _pca_cls()
    m = PCA(n_components=nc).fit(Xp)
    return {"kind": "pca", "model": m, "pre": pp, "pre_name": pp.name,
            "n_components": int(nc), "scores": np.asarray(m.transform(Xp)),
            "loadings": np.asarray(m.components_.T),
            "explained": np.asarray(m.explained_variance_ratio_),
            "caption": f"PCA ({pp.name}, {nc} PC): "
                       f"{100*float(np.sum(m.explained_variance_ratio_)):.1f}% variance"}


def pcr_fit(X, y, n_components, cfg, pre=None):
    """Principal-component regression: PCA, then OLS of y on the first n_components
    scores. Coefficients are mapped back to wavenumber space for the spectrum plot."""
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, y = _check_xy(X, y, cfg, n_components)
    pp = Preprocessor(pre, cfg)
    Xp = pp.fit_transform(X)
    nc = min(n_components, X.shape[0] - 1, Xp.shape[1])
    PCA = _pca_cls()
    pca = PCA(n_components=nc).fit(Xp)
    T = pca.transform(Xp)
    reg = _linreg_cls()().fit(T, y)
    yhat = np.ravel(reg.predict(T))
    coef = pca.components_.T @ np.ravel(reg.coef_)     # back to X space (p,)
    resid, r2, rmsec = _fit_stats(y, yhat)
    return {"kind": "pcr", "model": (pca, reg), "pre": pp, "pre_name": pp.name,
            "n_components": int(nc), "coef": coef, "scores": np.asarray(T),
            "loadings": np.asarray(pca.components_.T),
            "explained": np.asarray(pca.explained_variance_ratio_),
            "yhat": yhat, "resid": resid, "r2_cal": r2, "rmsec": rmsec, "y": y,
            "caption": f"PCR ({pp.name}, {nc} PC): R2(cal)={r2:.3f}, RMSEC={rmsec:.3g}"}


# ===================================================== cross-validation (leakage-safe)
def _make_folds(n, y, groups, scheme, cfg):
    """Return (folds, scheme_name, question). folds = list of (train_idx, test_idx)."""
    idx = np.arange(n)
    if groups is not None:
        groups = np.asarray(groups)
        uniq = list(dict.fromkeys(groups.tolist()))     # stable order
        folds = [(idx[groups != g], idx[groups == g]) for g in uniq]
        return (folds, f"leave-one-group-out ({len(uniq)} groups)",
                "predict a NEW group (e.g. an unseen concentration level)")
    sch = (scheme or getattr(cfg, "cv_scheme", "auto")).lower()
    if sch in ("auto", "loo"):
        folds = [(idx[idx != i], np.array([i])) for i in idx]
        return (folds, "leave-one-out (per sample)",
                "predict a new sample at a SEEN level "
                "(optimistic for a new level when samples replicate levels)")
    if sch == "kfold":
        k = int(min(getattr(cfg, "cv_folds", 5), n))
        order = np.argsort(y, kind="stable")            # spread each fold across the y range
        parts = [p for p in np.array_split(order, k) if len(p)]
        folds = [(np.setdiff1d(idx, p), p) for p in parts]
        return folds, f"{k}-fold (stratified by y)", "predict a held-out fold"
    raise ValueError(f"unknown cv_scheme: {sch!r} (use 'auto'/'loo'/'kfold' or pass groups=)")


def _fit_predict(Xtr, ytr, Xte, nc, cfg, method):
    nc = int(min(nc, len(Xtr) - 1, Xtr.shape[1]))
    if method == "pls":
        m = _pls_cls()(n_components=nc, scale=getattr(cfg, "pls_scale", True)).fit(Xtr, ytr)
        return np.ravel(m.predict(Xte))
    if method == "pcr":
        pca = _pca_cls()(n_components=nc).fit(Xtr)
        reg = _linreg_cls()().fit(pca.transform(Xtr), ytr)
        return np.ravel(reg.predict(pca.transform(Xte)))
    raise ValueError(f"unknown method: {method!r} (use 'pls' or 'pcr')")


def _cv_predict(X, y, nc, cfg, pre, folds, method):
    """Held-out predictions, preprocessing re-fit inside every fold (no leakage)."""
    pred = np.full(len(y), np.nan)
    for tr, te in folds:
        pp = Preprocessor(pre, cfg)                      # FRESH per fold
        Xtr = pp.fit_transform(X[tr])                    # fit on TRAIN rows only
        Xte = pp.transform(X[te])
        pred[te] = _fit_predict(Xtr, y[tr], Xte, nc, cfg, method)
    return pred


def _optimistic_loo(groups, scheme, y, cfg):
    sch = (scheme or getattr(cfg, "cv_scheme", "auto")).lower()
    return (groups is None and sch in ("auto", "loo")
            and np.unique(y).size < y.size)


def cross_validate(X, y, n_components, cfg, pre=None, groups=None, scheme=None,
                   method="pls", report=True):
    """Leakage-safe cross-validation. Re-fits the preprocessor inside each fold and
    NAMES which prediction question the scheme answers. Pass groups=level-labels for
    leave-one-level-out. Returns a dict with pred, rmsecv, r2cv, bias, scheme, question."""
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, y = _check_xy(X, y, cfg, n_components)
    folds, scheme_name, question = _make_folds(len(y), y, groups, scheme, cfg)
    pred = _cv_predict(X, y, n_components, cfg, pre, folds, method)
    resid = y - pred
    rmsecv = float(np.sqrt(np.nanmean(resid ** 2)))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2cv = 1 - float(np.nansum(resid ** 2)) / ss_tot if ss_tot else float("nan")
    bias = float(np.nanmean(pred - y))
    if report:
        f = []
        if _optimistic_loo(groups, scheme, y, cfg):
            f.append(("WARN", "LOO over data with replicated y answers 'predict a SEEN level'; "
                              "for a NEW level pass groups=level-labels (leave-one-level-out)"))
        f.append(("INFO", f"{scheme_name}: RMSECV={rmsecv:.3g} {getattr(cfg,'conc_unit','')}, "
                          f"R2cv={r2cv:.3f}, bias={bias:+.3g}, {len(folds)} folds"))
        _resolve(f, cfg, "cross-validate")
    return {"pred": pred, "resid": resid, "rmsecv": rmsecv, "r2cv": r2cv, "bias": bias,
            "n_components": int(n_components), "method": method,
            "pre_name": Preprocessor(pre, cfg).name, "scheme": scheme_name,
            "question": question, "n_splits": len(folds),
            "caption": f"RMSECV={rmsecv:.3g} {getattr(cfg,'conc_unit','')} "
                       f"({scheme_name}); answers: {question}"}


def component_scan(X, y, cfg, max_components=None, pre_options=None, groups=None,
                   scheme=None, method="pls"):
    """RMSECV against component count, ONE curve per preprocessing, all on the same
    folds (so the curves are comparable). Drives choose_n_components and plot_rmsecv."""
    _, X, y = _check_xy(X, y, cfg)
    n = len(y)
    folds, scheme_name, question = _make_folds(n, y, groups, scheme, cfg)
    min_train = min(len(tr) for tr, _ in folds)
    cap = int(min(max_components or getattr(cfg, "pls_max_components", 10),
                  min_train - 1, X.shape[1]))
    cap = max(cap, 1)
    if pre_options is None:
        pre_options = (cfg.pls_preprocess,)
    comps = list(range(1, cap + 1))
    curves = {}
    f = []
    if _optimistic_loo(groups, scheme, y, cfg):
        f.append(("WARN", "LOO over replicated y is optimistic for a NEW level; "
                          "pass groups=level-labels for leave-one-level-out"))
    for pre in pre_options:
        name = Preprocessor(pre, cfg).name
        rms = []
        for nc in comps:
            pred = _cv_predict(X, y, nc, cfg, pre, folds, method)
            rms.append(float(np.sqrt(np.nanmean((y - pred) ** 2))))
        curves[name] = rms
        f.append(("INFO", f"{name:14s} RMSECV " + " ".join(f"{r:.2f}" for r in rms)))
    _resolve(f, cfg, f"component-scan[{scheme_name}]")
    return {"curves": curves, "components": comps, "max": cap, "method": method,
            "scheme": scheme_name, "question": question}


def choose_n_components(scan, cfg, tol=None):
    """Parsimony pick: the FEWEST components whose RMSECV is within `tol` (fraction)
    of the global-min RMSECV across all scanned preprocessings. Ties broken by lower
    RMSECV. Avoids the over-fit you get from taking the global argmin."""
    tol = getattr(cfg, "lv_parsimony_tol", 0.10) if tol is None else tol
    flat = [(pre, nc, scan["curves"][pre][i])
            for pre in scan["curves"] for i, nc in enumerate(scan["components"])]
    gmin = min(v for _, _, v in flat)
    cands = [t for t in flat if t[2] <= gmin * (1 + tol)]
    pre, nc, rmsecv = min(cands, key=lambda t: (t[1], t[2]))
    return {"pre": pre, "n_components": int(nc), "rmsecv": float(rmsecv),
            "global_min": float(gmin), "tol": tol,
            "caption": f"{pre}, {nc} component(s) (RMSECV={rmsecv:.3g}); parsimony: "
                       f"fewest components within {int(tol*100)}% of the global min {gmin:.3g}"}


def vip(model):
    """VIP scores: per-wavenumber importance for a PLS model (>1 ~ above-average
    contribution). Needs a PLS model dict (weights, scores, y_loadings)."""
    W = np.asarray(model["weights"])                    # (p, A)
    T = np.asarray(model["scores"])                     # (n, A)
    q = np.ravel(model["y_loadings"])                   # (A,)
    p, A = W.shape
    Wn = W / np.where(np.linalg.norm(W, axis=0) == 0, 1.0, np.linalg.norm(W, axis=0))
    ssy = (q ** 2) * np.sum(T ** 2, axis=0)             # explained Y SS per component
    total = float(np.sum(ssy))
    if total == 0:
        return np.zeros(p)
    return np.sqrt(p * ((Wn ** 2) @ ssy) / total)


# ===================================================== ax-level plot helpers
_CYCLE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
_MARKERS = ["o", "s", "^", "D", "v", "P"]


def _xaxis(ax, x, cfg):
    """Honour FTIR inversion / cfg.x_limits for a wavenumber x-axis."""
    if cfg.x_limits:
        ax.set_xlim(*cfg.x_limits)
    elif getattr(cfg, "domain", "ftir") == "ftir":
        ax.set_xlim(float(np.max(x)), float(np.min(x)))


def plot_rmsecv(ax, scan, cfg, choice=None):
    """RMSECV vs component count, one line per preprocessing; rings the chosen point."""
    comps = scan["components"]
    for i, (name, rms) in enumerate(scan["curves"].items()):
        ax.plot(comps, rms, "-" + _MARKERS[i % len(_MARKERS)],
                color=_CYCLE[i % len(_CYCLE)], ms=4, lw=1, label=name)
    if choice is not None:
        ax.scatter([choice["n_components"]], [choice["rmsecv"]], s=90, zorder=5,
                   facecolor="none", edgecolor="k", linewidth=1.2)
    ax.set_xticks(comps)
    ax.set_xlabel("Components (latent variables)")
    ax.set_ylabel(f"RMSECV / {getattr(cfg, 'conc_unit', 'a.u.')}")
    if len(scan["curves"]) > 1:
        ax.legend(fontsize=6.5)
    return ax


def plot_coefficients(ax, model, x, cfg, ref=None, ref_label="reference", highlight=None):
    """Regression-coefficient spectrum (which wavenumbers drive the prediction), with
    an optional pure-component reference overlay (scaled) and a highlighted band."""
    x = np.asarray(x, float)
    coef = np.asarray(model["coef"], float)
    ax.axhline(0, color="0.7", lw=0.5)
    if highlight is not None:
        ax.axvspan(min(highlight), max(highlight), color="0.82", alpha=0.6, zorder=0)
    if ref is not None:
        ref = np.asarray(ref, float)
        scale = np.max(np.abs(coef)) / (np.max(np.abs(ref)) or 1.0)
        ax.fill_between(x, ref * scale, color="#56B4E9", alpha=0.35, lw=0, label=ref_label)
    ax.plot(x, coef, color="#222222", lw=0.8, zorder=3, label="coef.")
    _xaxis(ax, x, cfg)
    ax.set_xlabel(r"Wavenumber / cm$^{-1}$" if getattr(cfg, "domain", "ftir") == "ftir"
                  else r"2$\theta$ / deg")
    ax.set_ylabel("Regression coef. (a.u.)")
    ax.legend(fontsize=6, loc="best")
    return ax


def plot_scores(ax, model, cfg, color_by=None, comps=(0, 1), cbar=None, cbar_label=None):
    """Scores scatter (the calibration's sample structure). Colour by concentration to
    see the gradient. Returns the PathCollection so the caller can add a colorbar
    (do that AFTER finalize on a constrained-layout figure so it doesn't collide)."""
    T = np.asarray(model["scores"])
    i, j = comps
    kw = dict(s=26, edgecolor="white", linewidth=0.3)
    if color_by is not None:
        sc = ax.scatter(T[:, i], T[:, j], c=np.asarray(color_by, float), cmap="viridis", **kw)
    else:
        sc = ax.scatter(T[:, i], T[:, j], color=_CYCLE[0], **kw)
    tag = "LV" if model.get("kind") == "pls" else "PC"
    ax.set_xlabel(f"{tag}{i+1} score (a.u.)")
    ax.set_ylabel(f"{tag}{j+1} score (a.u.)")
    if cbar is not None and color_by is not None:
        cb = cbar.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(cbar_label or getattr(cfg, "conc_unit", "a.u."), fontsize=7)
    return sc


def plot_loadings(ax, model, x, cfg, comps=(0, 1), highlight=None):
    """Loadings (the latent variables themselves) vs wavenumber."""
    x = np.asarray(x, float)
    P = np.asarray(model["loadings"])
    tag = "LV" if model.get("kind") == "pls" else "PC"
    ax.axhline(0, color="0.7", lw=0.5)
    if highlight is not None:
        ax.axvspan(min(highlight), max(highlight), color="0.82", alpha=0.6, zorder=0)
    for k, c in enumerate(comps):
        ax.plot(x, P[:, c], color=_CYCLE[k % len(_CYCLE)], lw=0.8,
                alpha=1.0 if k == 0 else 0.75, label=f"{tag}{c+1}")
    _xaxis(ax, x, cfg)
    ax.set_xlabel(r"Wavenumber / cm$^{-1}$" if getattr(cfg, "domain", "ftir") == "ftir"
                  else r"2$\theta$ / deg")
    ax.set_ylabel(f"{tag} loading (a.u.)")
    ax.legend(fontsize=6.5)
    return ax


# ===================================================== the composite deliverable
def diagnostics_figure(X, y, x, cfg, groups=None, pre_options=None, ref=None,
                       ref_label="reference", highlight=None, viz_pre=None,
                       viz_components=None, basename="pls_diagnostics", method="pls"):
    """The 2x2 PLS diagnostic the family exists to produce, gated and QA'd:
       (a) RMSECV vs components, one curve per preprocessing, parsimony pick ringed
       (b) regression-coefficient spectrum (+ optional pure-component reference)
       (c) scores coloured by concentration
       (d) loadings of the first two components

    The honest error (a) is computed from a leakage-safe scan; (b)-(d) are drawn from
    a model fit on the whole set for interpretability. By default that model uses the
    parsimony pick, but pass viz_pre / viz_components to show a more band-interpretable
    model (e.g. raw, 2 LV) while still REPORTING the SNV/1-LV error in (a).

    Returns (fig, axes, info) where info has the scan, the parsimony choice, and the
    fitted viz model. Pass groups=level-labels to make (a) a leave-one-level-out error."""
    import matplotlib.pyplot as plt
    _, X, y = _check_xy(X, y, cfg)
    x = np.asarray(x, float)
    if pre_options is None:
        pre_options = tuple(dict.fromkeys((cfg.pls_preprocess, "none")))
    scan = component_scan(X, y, cfg, pre_options=pre_options, groups=groups, method=method)
    choice = choose_n_components(scan, cfg)
    fitter = {"pls": pls_fit, "pcr": pcr_fit}[method]
    model = fitter(X, y, viz_components or choice["n_components"], cfg,
                   pre=viz_pre or choice["pre"])

    apply_style(cfg)
    fig, ax = figure(cfg, 2, 2, layout="constrained")   # colorbar must be layout-aware
    (a, b), (c, d) = ax
    plot_rmsecv(a, scan, cfg, choice=choice)
    a.set_title("How many components?", fontsize=8.5)
    plot_coefficients(b, model, x, cfg, ref=ref, ref_label=ref_label, highlight=highlight)
    b.set_title(f"Coefficients ({model['pre_name']})", fontsize=8.5)
    plot_scores(c, model, cfg, color_by=y, cbar=fig)
    c.set_title("Scores", fontsize=8.5)
    comps = (0, 1) if model["scores"].shape[1] > 1 else (0,)
    plot_loadings(d, model, x, cfg, comps=comps, highlight=highlight)
    d.set_title("Loadings", fontsize=8.5)

    finalize_figure(fig, wspace=0.10, hspace=0.12)       # constrained-layout fractions
    # left column clears its own narrow y-furniture; right column gets a tight offset so
    # b/d sit at their corners instead of being dragged left by d's wide decimal ticks.
    add_panel_labels(fig, cfg, axes=[a, c], labels=["a", "c"], x_offset_pt="auto")
    add_panel_labels(fig, cfg, axes=[b, d], labels=["b", "d"], x_offset_pt=-10.0)
    import os
    render_preview(fig, os.path.join(cfg.output_dir, "_preview_" + basename + ".png"))
    audit_layout(fig, cfg)
    return fig, (a, b, c, d), {"scan": scan, "choice": choice, "model": model}


# ===================================================== methods text
def methods_text(model, cv, choice, cfg):
    """Paste-ready methods paragraph: matrix, preprocessing, model + component count
    (with the selection rule), CV scheme + the question it answers, and the figures of
    merit. Edit to taste; do not invent numbers it doesn't contain."""
    n = model["y"].size
    p = np.asarray(model["coef"]).size
    tag = "latent variables" if model["kind"] == "pls" else "principal components"
    lines = [
        f"Multivariate calibration was performed by {model['kind'].upper()} on the "
        f"baseline-corrected spectra ({n} samples x {p} wavenumbers).",
        f"Spectra were preprocessed by {model['pre_name']} before modelling.",
        f"The number of {tag} was chosen by cross-validation: {choice['caption']}.",
        f"Cross-validation used {cv['scheme']}, which answers '{cv['question']}'.",
        f"Figures of merit: RMSEC={model.get('rmsec', float('nan')):.3g} "
        f"{getattr(cfg, 'conc_unit', 'a.u.')}, R2(cal)={model.get('r2_cal', float('nan')):.3f}; "
        f"RMSECV={cv['rmsecv']:.3g} {getattr(cfg, 'conc_unit', 'a.u.')}, "
        f"R2(CV)={cv['r2cv']:.3f}, bias={cv['bias']:+.3g}.",
    ]
    return "\n".join(lines)


# ===================================================== outlier / diagnostic gate
# PCA-based per-sample diagnostics, the modern chemometric outlier suite:
#   * Hotelling T2  (score distance) -- "extreme but maybe valid" position in the model
#   * Q / SPE       (orthogonal distance) -- "doesn't fit the model" (the usual bad-load tell)
#   * leverage      (Sigma U^2; the unknown-sample extrapolation flag at 2A/n or 3A/n)
#   * DD-SIMCA      (Pomerantsev & Rodionova 2014) -- SD + OD combined into a chi^2(Nh+Nq)
#                   statistic with an extreme limit (alpha) and a Bonferroni outlier limit (gamma)
# numpy + scipy only (NO sklearn) so it runs as raw-spectrum QC BEFORE any  The
# limit forms are exactly those in the literature (Jackson & Mudholkar 1979 for Q; the F / chi^2
# / Tracy-Young-Mason 1992 forms for T2). Validated bit-for-bit against the hand-derived fixture
# in skill_validation/chemometrics/golden/closed_form.py, plus scipy F/chi2 closed forms.
#
# Why a separate PCA decomposition rather than reusing pca_fit's dict: this is intentionally
# self-contained (numpy SVD), stores the training mean/std so NEW unknowns can be scored, and
# matches the closed-form derivation exactly. T2 = (n-1)*leverage; mean(leverage) = A/n.

def pca_outlier_model(X, cfg=None, n_components=2, pre="none", scale=False):
    """Decompose X for outlier diagnostics: preprocess (`pre`), centre (store the mean),
    optionally autoscale (ddof=1), SVD, retain A = min(n_components, numerical rank).
    `pre` defaults to "none" (raw-spectrum QC); pass e.g. "snv" to QC the modelled space.
    Stores the training mean/std so `outlier_distances(om, Xnew)` can score new samples."""
    X = np.asarray(X, float)
    pp = Preprocessor(pre, cfg)
    Xp = pp.fit_transform(X)
    mean = Xp.mean(0)
    Xc = Xp - mean
    std = None
    if scale:
        sd = Xc.std(0, ddof=1)
        std = np.where(sd == 0, 1.0, sd)
        Xc = Xc / std
    n, p = Xc.shape
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    tol = (s.max() * max(n, p) * np.finfo(float).eps) if s.size else 0.0
    rank = int(np.sum(s > tol))
    A = int(max(1, min(n_components, rank)))
    eig = (s ** 2) / (n - 1)                        # eigenvalues (ddof=1)
    return {"kind": "pca_outlier", "pre": pp, "mean": mean, "std": std, "scale": scale,
            "s": s, "eig": eig, "P": Vt.T, "U": U, "scores": U * s,
            "A": A, "rank": rank, "n": n, "p": p}


def _outlier_centred(om, X):
    if X is None:
        return om["scores"] @ om["P"].T            # exact reconstruction of the training Xc
    Xp = om["pre"].transform(np.asarray(X, float))
    Xc = Xp - om["mean"]
    if om["scale"] and om["std"] is not None:
        Xc = Xc / om["std"]
    return Xc


def outlier_distances(om, X=None):
    """Per-sample T2, Q/SPE, leverage and the DD-SIMCA score/orthogonal distances, for the
    training set (X=None) or for new samples X. Q = ||Xc||^2 - sum(T_A^2) is general (it also
    captures the part of a new sample outside the model span). leverage = sum (T_a/s_a)^2."""
    A = om["A"]; n = om["n"]
    s = om["s"][:A]
    eig = om["eig"][:A]
    Xc = _outlier_centred(om, X)
    T = Xc @ om["P"]
    Ta = T[:, :A]
    t2 = np.sum(Ta ** 2 / np.where(eig == 0, np.inf, eig), axis=1)
    q = np.sum(Xc ** 2, axis=1) - np.sum(Ta ** 2, axis=1)
    q = np.clip(q, 0.0, None)                       # numerical floor (perfect fit -> 0)
    lev = np.sum((Ta / np.where(s == 0, np.inf, s)) ** 2, axis=1)   # = T2/(n-1); mean = A/n
    return {"t2": t2, "q": q, "leverage": lev, "sd": lev, "od": q, "scores": T}


def t2_limit(om, alpha=0.05, form="F_textbook"):
    """Hotelling T2 critical limit. forms: 'F_textbook' [A(n-1)/(n-A)]F (mdatools),
    'F_n2m1' [A(n^2-1)/(n(n-A))]F (chemometrics lib), 'chi2' (large-n), 'beta'
    (Tracy-Young-Mason, Phase-I). Falls back to chi2 when n-A<=0."""
    from scipy.stats import f as fdist, chi2, beta
    A = om["A"]; n = om["n"]
    if form == "chi2" or n - A <= 0:
        return float(chi2.ppf(1 - alpha, A))
    fq = float(fdist.ppf(1 - alpha, A, n - A))
    if form == "F_textbook":
        return A * (n - 1) / (n - A) * fq
    if form == "F_n2m1":
        return A * (n ** 2 - 1) / (n * (n - A)) * fq
    if form == "beta":
        if n - A - 1 <= 0:
            return float(chi2.ppf(1 - alpha, A))
        return (n - 1) ** 2 / n * float(beta.ppf(1 - alpha, A / 2.0, (n - A - 1) / 2.0))
    raise ValueError(f"unknown T2 limit form: {form!r}")


def q_limit(om, alpha=0.05, method="jackson_mudholkar", q_values=None):
    """Q/SPE critical limit. 'jackson_mudholkar' uses the discarded eigenvalues (theta1-3,h0);
    'box' fits g*chi2_h to the moments of the Q values (needs q_values). Returns 0.0 when there
    is no residual space (A>=rank) or a perfect fit -- the Q test is then uninformative."""
    from scipy.stats import norm, chi2
    if method == "jackson_mudholkar":
        res = om["eig"][om["A"]:om["rank"]]
        res = res[res > 0]
        if res.size == 0:
            return 0.0
        th1 = float(res.sum()); th2 = float((res ** 2).sum()); th3 = float((res ** 3).sum())
        if th2 == 0:
            return 0.0
        h0 = 1.0 - 2.0 * th1 * th3 / (3.0 * th2 ** 2)
        if h0 <= 0:
            h0 = 1e-6
        z = float(norm.ppf(1 - alpha))
        return float(th1 * (z * np.sqrt(2.0 * th2 * h0 ** 2) / th1 + 1.0
                            + th2 * h0 * (h0 - 1.0) / th1 ** 2) ** (1.0 / h0))
    if method == "box":
        if q_values is None:
            raise ValueError("box method needs q_values")
        qv = np.asarray(q_values, float)
        m = float(qv.mean()); v = float(qv.var(ddof=1))
        if m <= 0 or v <= 0:
            return 0.0
        g = v / (2.0 * m); h = 2.0 * m ** 2 / v
        return float(g * chi2.ppf(1 - alpha, h))
    raise ValueError(f"unknown Q limit method: {method!r}")


def _dd_dof(d, method="classical"):
    """DD-SIMCA data-driven degrees of freedom + scaling for a distance vector, under the
    parameterization u ~ (u0/N)*chi2(N) (so N*u/u0 ~ chi2(N)).

    method='classical' (mdatools ddmoments): N = round(2*mean^2/var), u0 = mean. EXACT but
        NON-ROBUST — a strong outlier inflates mean/var and masks itself. Right for clean data.
    method='robust' (mdatools ddrobust; Pomerantsev & Rodionova 2014): solve N from the
        IQR/median ratio of chi2(N) (which is scale-free), then u0 = N*median/median(chi2_N).
        Use this for DETECTION — the median/IQR are not dragged by the outliers being sought.
    Both clamped to [1,250]."""
    from scipy.stats import chi2
    d = np.asarray(d, float)
    if method == "robust":
        M = float(np.median(d))
        R = float(np.percentile(d, 75) - np.percentile(d, 25))
        if M <= 0 or R <= 0:
            return 1, max(M, 1e-12)
        ratio = R / M

        def g(nu):
            return (chi2.ppf(0.75, nu) - chi2.ppf(0.25, nu)) / chi2.ppf(0.5, nu) - ratio
        lo, hi = 1.0, 250.0
        glo, ghi = g(lo), g(hi)
        if glo * ghi < 0:
            try:
                from scipy.optimize import brentq
                nu = brentq(g, lo, hi, xtol=1e-4)
            except Exception:
                nu = lo
        else:
            nu = lo if abs(glo) <= abs(ghi) else hi      # ratio outside the chi2 range -> clamp
        N = max(1, min(250, int(round(nu))))
        u0 = N * M / float(chi2.ppf(0.5, N))
        return N, u0
    m = float(d.mean()); v = float(d.var(ddof=1))
    N = int(round(2.0 * m ** 2 / v)) if v > 0 else 1
    return max(1, min(250, N)), m


def ddsimca_limits(om, sd, od, alpha=0.05, gamma=0.01, dof="classical"):
    """DD-SIMCA combined-distance limits. f = Nh*(SD/h0) + Nq*(OD/q0) ~ chi2(Nh+Nq); extreme
    limit at alpha; per-sample Bonferroni outlier limit at gamma (None => no outlier line).
    `dof`: 'classical' (moments) or 'robust' (median/IQR — use for detection, see _dd_dof)."""
    from scipy.stats import chi2
    n = om["n"]
    Nh, h0 = _dd_dof(sd, dof)
    Nq, q0 = _dd_dof(od, dof)
    c_crit = float(chi2.ppf(1 - alpha, Nh + Nq))
    c_out = float(chi2.ppf((1 - gamma) ** (1.0 / n), Nh + Nq)) if gamma else None
    return {"Nh": Nh, "Nq": Nq, "h0": h0, "q0": q0, "c_crit": c_crit, "c_out": c_out, "dof": dof}


def diagnostics(X, cfg=None, n_components=2, pre="none", scale=False, alpha=0.05, gamma=0.01,
                t2_form="F_textbook", q_method="jackson_mudholkar", lev_factor=2, dof="robust"):
    """The full outlier/diagnostic read-out for a spectral matrix: builds the PCA model,
    computes T2/Q/leverage + their limits, the DD-SIMCA combined distance + flags
    (regular/extreme/outlier). Returns one dict; feed it to plot_influence. The DD-SIMCA
    path is the primary verdict; T2(F) and Q(Jackson-Mudholkar) are the classical cross-checks.
    `dof` defaults to 'robust' (median/IQR) so strong outliers don't mask themselves; pass
    'classical' for the method-of-moments DoF (right only when the set is already clean)."""
    om = pca_outlier_model(X, cfg, n_components, pre, scale)
    d = outlier_distances(om)
    dd = ddsimca_limits(om, d["sd"], d["od"], alpha, gamma, dof=dof)
    f = dd["Nh"] * (d["sd"] / dd["h0"]) + dd["Nq"] * (d["od"] / dd["q0"])
    if dd["c_out"] is not None:
        flags = np.where(f > dd["c_out"], "outlier",
                         np.where(f > dd["c_crit"], "extreme", "regular"))
    else:
        flags = np.where(f > dd["c_crit"], "extreme", "regular")
    return {"model": om, **d,
            "t2_limit": t2_limit(om, alpha, t2_form), "t2_form": t2_form,
            "q_limit": q_limit(om, alpha, q_method, q_values=d["q"]), "q_method": q_method,
            "leverage_warn": lev_factor * om["A"] / om["n"],
            "dd": dd, "f": f, "flags": flags, "alpha": alpha, "gamma": gamma,
            "caption": (f"PCA outlier diagnostics ({om['A']} PC, {om['pre'].name}): "
                        f"DD-SIMCA[{dd['dof']}] chi^2({dd['Nh']}+{dd['Nq']}), extreme alpha={alpha}"
                        + (f", outlier gamma={gamma}" if gamma else "")
                        + f"; {int(np.sum(flags!='regular'))} of {om['n']} flagged.")}


def plot_influence(ax, diag, cfg=None, labels=None, axis_scale="linear"):
    """DD-SIMCA acceptance ('influence') plot: orthogonal distance OD/q0 vs score distance
    SD/h0, with the chi^2 boundary drawn as a straight line (extreme dashed, outlier dotted).
    Points keyed regular/extreme/outlier by marker+colour (grayscale-safe). Annotate by
    passing labels (e.g. sample ids) — only flagged points get a tag.

    axis_scale: 'linear' (default) shows the true linear boundary — best when the distances are
    of comparable size. Use 'sqrt' when one or a few GROSS outliers dominate the range and squash
    the in-control cluster and the limit lines into the corner (the single-bad-spectrum case): a
    sqrt axis compresses the large values so the borderline/extreme points and the boundary stay
    legible. The boundary stays correct — it is drawn in data coordinates and the axis transform
    curves it; sqrt (not log) is the default transform because SD/OD are >=0 and often exactly 0."""
    dd = diag["dd"]
    xs = np.asarray(diag["leverage"]) / dd["h0"]
    ys = np.asarray(diag["q"]) / dd["q0"]
    flags = np.asarray(diag["flags"])
    sty = {"regular": ("#0072B2", "o"), "extreme": ("#E69F00", "s"), "outlier": ("#D55E00", "D")}
    for fl, (col, mk) in sty.items():
        m = flags == fl
        if np.any(m):
            ax.scatter(xs[m], ys[m], c=col, marker=mk, s=30, edgecolor="white",
                       linewidth=0.3, label=fl, zorder=3)
    xmax = max(float(xs.max()) * 1.1, dd["c_crit"] / dd["Nh"]) if xs.size else 1.0
    xl = np.linspace(0, xmax, 100)
    for c, ls, lab in [(dd["c_crit"], "--", f"extreme (α={diag['alpha']})"),
                       (dd["c_out"], ":", f"outlier (γ={diag['gamma']})" if dd["c_out"] else None)]:
        if c is None:
            continue
        yb = (c - dd["Nh"] * xl) / dd["Nq"]
        ok = yb >= 0
        ax.plot(xl[ok], yb[ok], ls, color="0.35", lw=0.9, label=lab, zorder=2)
    if labels is not None:
        for i, fl in enumerate(flags):
            if fl != "regular":
                ax.annotate(str(labels[i]), (xs[i], ys[i]), fontsize=6,
                            xytext=(2, 2), textcoords="offset points")
    if axis_scale == "sqrt":
        fwd = lambda a: np.sqrt(np.clip(np.asarray(a, float), 0, None))
        inv = lambda a: np.asarray(a, float) ** 2
        ax.set_xscale("function", functions=(fwd, inv))
        ax.set_yscale("function", functions=(fwd, inv))
    elif axis_scale not in ("linear", None):
        raise ValueError(f"axis_scale must be 'linear' or 'sqrt', got {axis_scale!r}")
    ax.set_xlabel("score distance / $h_0$ (ratio)")
    ax.set_ylabel("orthogonal distance / $q_0$ (ratio)")
    ax.legend(fontsize=6, loc="best")
    return ax


# =============================================== DD-SIMCA one-class classifier + class FoM
# The outlier gate above asks "is this sample weird for THIS set?"; a one-class CLASSIFIER asks
# "is this sample a member of the TARGET class?" — the right frame for cocrystal ID ("is this the
# target phase, or a physical mixture / a different polymorph?"). Fit the DD-SIMCA acceptance model
# on the target class's spectra, then accept/reject unknowns against its chi^2 boundary, and report
# the class figures-of-merit (sensitivity / specificity / efficiency) — NOT R2/LOD, the wrong
# vocabulary for a classifier. Reuses the validated pca_outlier_model / ddsimca_limits path, so its
# limits inherit test_outliers' validation. Routed to from references/cocrystal_id.md.

def ddsimca_fit(X, cfg=None, n_components=2, pre="none", scale=False,
                alpha=0.05, gamma=None, dof="classical"):
    """Fit a DD-SIMCA ONE-CLASS model on the TARGET class's spectra: a PCA model + the SD/OD
    scaling (h0, q0, Nh, Nq) + the acceptance limit c_crit at 1-alpha. Feed the model to
    `ddsimca_predict` to accept/reject unknowns ('is this the target phase?'). `dof` defaults to
    'classical' — a training class is assumed clean; use 'robust' only if it may itself contain
    outliers. gamma=None → no per-sample outlier line (a class model needs only the acceptance
    boundary)."""
    om = pca_outlier_model(X, cfg, n_components, pre, scale)
    d = outlier_distances(om)
    dd = ddsimca_limits(om, d["sd"], d["od"], alpha, gamma, dof=dof)
    return {"kind": "ddsimca", "om": om, "dd": dd, "alpha": alpha,
            "train_sd": d["sd"], "train_od": d["od"],
            "caption": (f"DD-SIMCA one-class model ({om['A']} PC, {om['pre'].name}): accept if the "
                        f"chi^2({dd['Nh']}+{dd['Nq']}) combined distance <= c_crit (alpha={alpha}).")}


def ddsimca_predict(model, X):
    """Accept/reject samples X against a `ddsimca_fit` model. Returns dict: `f` (combined
    distance), `accept` (bool — f <= c_crit → IN the target class), `sd`, `od`, `c_crit`. A
    rejected unknown is NOT the modelled class: a target-cocrystal model rejects a physical
    mixture or a different polymorph. X is projected through the STORED model (out-of-sample)."""
    om = model["om"]; dd = model["dd"]
    dn = outlier_distances(om, np.atleast_2d(np.asarray(X, float)))
    f = dd["Nh"] * (dn["sd"] / dd["h0"]) + dd["Nq"] * (dn["od"] / dd["q0"])
    return {"f": f, "accept": f <= dd["c_crit"], "sd": dn["sd"], "od": dn["od"],
            "c_crit": dd["c_crit"]}


def class_fom(y_true, y_pred):
    """Classification figures-of-merit — the CORRECT vocabulary for a class model (DD-SIMCA /
    SIMCA / PLS-DA), not R2/LOD. `y_true`, `y_pred` are boolean arrays (True = member of the
    target class / accepted):
      sensitivity (TPR) = TP/(TP+FN)  — target members correctly accepted,
      specificity (TNR) = TN/(TN+FP)  — non-members correctly rejected,
      efficiency        = sqrt(sensitivity·specificity)  — the SIMCA class efficiency.
    Returns those + the confusion counts + a caption. NaN-safe (a rate with no denominator is nan)."""
    yt = np.asarray(y_true, bool).ravel()
    yp = np.asarray(y_pred, bool).ravel()
    if yt.shape != yp.shape:
        raise ValueError("y_true and y_pred must have the same shape")
    tp = int(np.sum(yt & yp)); fn = int(np.sum(yt & ~yp))
    tn = int(np.sum(~yt & ~yp)); fp = int(np.sum(~yt & yp))
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    eff = float(np.sqrt(sens * spec)) if (np.isfinite(sens) and np.isfinite(spec)) else float("nan")
    caption = (f"class FoM — sensitivity={sens:.3f} (TPR {tp}/{tp+fn}), "
               f"specificity={spec:.3f} (TNR {tn}/{tn+fp}), efficiency={eff:.3f} "
               f"(geom. mean) [SIMCA class metrics, not R2/LOD]")
    return {"sensitivity": sens, "specificity": spec, "efficiency": eff,
            "tp": tp, "fn": fn, "tn": tn, "fp": fp, "caption": caption}


# ===================================================== validation trio
# Cheap, defensible figures-of-merit for a SMALL calibration: a permutation/y-scramble test
# (is the model better than chance?), a corrected paired t-test (is model A really better than
# B?), and RPD/RPIQ (is the calibration useful?).

def _permute_y(y, groups, rng):
    """Shuffle y for the null. With groups (level labels) permute BETWEEN levels (each level's
    reps keep ONE shared y) — replicate structure is preserved; sklearn's groups= permutes
    WITHIN groups, which is wrong for replicate  Without groups, a plain full shuffle."""
    y = np.asarray(y, float)
    if groups is None:
        return rng.permutation(y)
    groups = np.asarray(groups)
    uniq = list(dict.fromkeys(groups.tolist()))
    vals = rng.permutation(np.array([y[groups == g].mean() for g in uniq]))
    yp = np.empty_like(y)
    for g, v in zip(uniq, vals):
        yp[groups == g] = v
    return yp


def permutation_test(X, y, n_components, cfg, n_perm=199, groups=None, scheme=None,
                     method="pls", pre=None, seed=0):
    """Permutation / y-scrambling test: refit on shuffled y n_perm times to build the null for
    R2(CV), then p = (#{null >= observed} + 1) / (n_perm + 1). With groups, y is permuted between
    LEVELS (see _permute_y). Guards a small calibration against a chance correlation."""
    rng = np.random.default_rng(seed)
    obs = cross_validate(X, y, n_components, cfg, pre=pre, groups=groups, scheme=scheme,
                         method=method, report=False)["r2cv"]
    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = _permute_y(y, groups, rng)
        null[i] = cross_validate(X, yp, n_components, cfg, pre=pre, groups=groups, scheme=scheme,
                                 method=method, report=False)["r2cv"]
    p = (int(np.sum(null >= obs)) + 1) / (n_perm + 1)
    return {"observed": float(obs), "null": null, "p_value": float(p), "n_perm": n_perm,
            "caption": f"permutation test: R2(CV)={obs:.3f}, p={p:.3g} ({n_perm} permutations)"}


def corrected_paired_t(scores_a, scores_b, n_train, n_test):
    """Nadeau & Bengio (2003) CORRECTED resampled paired t-test comparing two models' per-resample
    scores/errors. The naive paired-t variance is inflated by (1/k + n_test/n_train) to account
    for the training sets overlapping across resamples. Returns t, two-sided p, df=k-1, correction.
    (mlxtend's paired_ttest_resampled is the UNcorrected Dietterich formula — not this.)"""
    from scipy.stats import t as tdist
    a = np.asarray(scores_a, float); b = np.asarray(scores_b, float)
    d = a - b
    k = d.size
    md = float(d.mean()); var = float(d.var(ddof=1))
    correction = 1.0 / k + float(n_test) / float(n_train)
    se = np.sqrt(correction * var)
    if se > 0:
        t = md / se
    else:
        t = 0.0 if md == 0 else float(np.inf) * np.sign(md)
    p = float(2 * tdist.sf(abs(t), k - 1)) if np.isfinite(t) else 0.0
    return {"t": float(t), "p_value": p, "df": k - 1, "correction": float(correction),
            "mean_diff": md,
            "caption": f"corrected paired t = {t:.3f}, p = {p:.3g} (df={k-1}, "
                       f"correction={correction:.3f}; Nadeau-Bengio)"}


def rpd(y_reference, rmse, ddof=1):
    """Ratio of Performance to Deviation = SD(reference)/RMSE. SD uses ddof=1 (sample). Rough
    interpretation (Chang 2001): >2 useful, >2.5 good, >3 excellent. Pass RMSEP or RMSECV."""
    sd = float(np.std(np.asarray(y_reference, float), ddof=ddof))
    return float(sd / rmse) if rmse > 0 else float("inf")


def rpiq(y_reference, rmse):
    """Ratio of Performance to IQR = IQR(reference)/RMSE (IQR via R-type-7 / numpy default
    linear interpolation). RPIQ ~= 1.349*RPD for normal y; more robust to skew."""
    y = np.asarray(y_reference, float)
    iqr = float(np.percentile(y, 75) - np.percentile(y, 25))
    return float(iqr / rmse) if rmse > 0 else float("inf")


# ======================================================================
# from charts.py
# ======================================================================
"""
py  -  generic line / bar for lab data, with the chart-choice guards that
matter. The value here is interception of bad chart choices, not the plotting API.

    line(series, cfg)        x/y line chart (e.g. timecourse, kinetics)
    bar(groups, cfg)         bar chart that REFUSES bare mean bars at low n -
                             it shows the individual points (and an error bar whose
                             statistic is named) instead. Mean bars at n<min hide
                             the distribution and overstate certainty.
"""


def line(series, cfg, xlabel="x", ylabel="y", ax=None):
    """series: dict label -> (x, y)."""
    apply_style(cfg)
    if ax is None:
        fig, ax = figure(cfg)
    else:
        fig = ax.figure
    for label, (x, y) in series.items():
        check_trace(np.asarray(x), np.asarray(y), cfg, name=label)
        ax.plot(x, y, marker="o", label=label)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if len(series) > 1:
        ax.legend(loc="best")
    return fig, ax


def bar(groups, cfg, ylabel="value"):
    """groups: dict label -> list of replicate values.
    If every group has n >= cfg.min_n_for_mean_bar, draw mean bars with a NAMED
    error bar (mean ± CI) and overlay the points. If any group is below the
    threshold, REFUSE the bare bar and show a points-only plot (with the mean as
    a tick) - a mean bar at n<3 hides the distribution and overstates certainty."""
    apply_style(cfg)
    fig, ax = figure(cfg)
    labels = list(groups.keys())
    stats = {k: summarize(v, cfg) for k, v in groups.items()}
    min_n = min(s["n"] for s in stats.values())
    pal = _palette(cfg)
    xpos = np.arange(len(labels))

    low_n = min_n < cfg.min_n_for_mean_bar
    if low_n:
        print(f"  [WARN] charts: min n={min_n} < {cfg.min_n_for_mean_bar} - "
              f"refusing mean bars, showing individual points instead")
        for i, k in enumerate(labels):
            pts = np.asarray(groups[k], float)
            ax.scatter(np.full_like(pts, xpos[i]), pts,
                       color=pal[i % len(pal)], zorder=3)
            ax.plot([xpos[i] - 0.2, xpos[i] + 0.2], [stats[k]["mean"]] * 2,
                    color="k", lw=1.2)            # mean tick
        cap = f"points + mean tick (n={min_n}); CIs omitted at low n"
    else:
        means = [stats[k]["mean"] for k in labels]
        errs = [stats[k]["ci_halfwidth"] for k in labels]
        ax.bar(xpos, means, yerr=errs, capsize=2,
               color=[pal[i % len(pal)] for i in range(len(labels))],
               edgecolor="k", linewidth=0.5)
        for i, k in enumerate(labels):          # redundant point overlay
            pts = np.asarray(groups[k], float)
            ax.scatter(np.full_like(pts, xpos[i]), pts, color="k", s=6, zorder=3)
        cap = list(stats.values())[0]["caption"]

    ax.set_xticks(xpos); ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.annotate(cap, xy=(0.02, 0.98), xycoords="axes fraction",
                va="top", fontsize="small")
    return fig, ax


# ======================================================================
# from report.py
# ======================================================================
"""
py  -  turn the CONFIG + results into prose you can paste into a 

    methods_report(cfg, results=None) -> (methods_text, si_table)

methods_text   a Methods paragraph describing the conventions the figure was made
               under (instrument axis convention, baseline, normalisation, the
               statistic behind any error bars, and the figure/software settings).
si_table       a plain-text table of whatever results dict you pass in.

This writes the FORM (how it was processed and rendered). It does NOT invent
findings or interpret the data - that judgement is yours.
"""


def _stat_sentence(cfg):
    return (f"Uncertainties are reported as {int(cfg.conf_level*100)}% confidence "
            f"intervals (Student's t, two-sided) unless stated otherwise; the "
            f"statistic and n are given with each value.")


def methods_report(cfg, results=None):
    parts = []

    if cfg.domain == "ftir":
        axis = ("Infrared spectra are presented with wavenumber on the abscissa "
                "(cm^-1), plotted high-to-low by convention")
        yax = ("transmittance (%)" if cfg.ftir_yaxis == "transmittance"
               else "absorbance")
        axis += f", and {yax} on the ordinate."
    elif cfg.domain == "pxrd":
        axis = ("Powder X-ray diffractograms are presented with diffraction angle "
                "2theta (degrees) on the abscissa, plotted low-to-high, and "
                "intensity on the ordinate.")
    else:
        axis = ""
    if axis:
        parts.append(axis)

    if cfg.baseline:
        bl = {"arpls": f"an asymmetrically reweighted penalised least-squares "
                       f"(arPLS, lambda={cfg.arpls_lam:g}) baseline",
              "rubberband": "a rubber-band (convex-hull) baseline"}.get(
                  cfg.baseline, f"a {cfg.baseline} baseline")
        parts.append(f"Each trace was baseline-corrected with {bl}, applied "
                     f"identically across every replicate in a set.")

    if cfg.integration_windows:
        wins = ", ".join(f"{n} ({min(a,b):g}-{max(a,b):g} cm^-1)"
                         for a, b, n in cfg.integration_windows)
        parts.append(f"Band areas were integrated above a local linear baseline "
                     f"between fixed window endpoints ({wins}), identical for all "
                     f"spectra so that replicate differences are not a windowing "
                     f"artefact.")

    if cfg.normalize and cfg.normalize != "none":
        how = {"max": "to unit maximum", "area": "to unit area"}.get(cfg.normalize, cfg.normalize)
        parts.append(f"Spectra were normalised {how} for display only; "
                     f"quantification used the un-normalised, baseline-corrected trace.")

    if cfg.smooth == "savgol":
        parts.append(f"A Savitzky-Golay filter (window {cfg.savgol_window}, order "
                     f"{cfg.savgol_poly}) was applied for display only; areas were "
                     f"measured on the unsmoothed trace.")

    parts.append(_stat_sentence(cfg))

    pal = {"okabe_ito": "the Okabe-Ito colour-blind-safe palette",
           "colorblind": "a colour-blind-safe palette"}.get(cfg.palette, cfg.palette)
    parts.append(f"Figures were prepared at final {cfg.journal} {cfg.column}-column "
                 f"dimensions using {pal} with redundant line-style encoding, and "
                 f"exported as vector graphics ({'/'.join(cfg.formats)}; raster "
                 f"proofs at {cfg.dpi_raster} dpi).")

    methods_text = " ".join(parts)

    # ---- SI table ----
    si_table = ""
    if results:
        rows = [(str(k), (f"{v:.4g}" if isinstance(v, (int, float)) else str(v)))
                for k, v in results.items()]
        w0 = max((len(k) for k, _ in rows), default=8)
        w1 = max((len(v) for _, v in rows), default=6)
        line = f"+{'-'*(w0+2)}+{'-'*(w1+2)}+"
        si_table = "\n".join(
            [line, f"| {'quantity'.ljust(w0)} | {'value'.ljust(w1)} |", line]
            + [f"| {k.ljust(w0)} | {v.ljust(w1)} |" for k, v in rows] + [line])

    return methods_text, si_table


# ======================================================================
# from crystal_engine.py
# ======================================================================
"""
py  -  CIF -> validated crystal data. The single source of truth for the
crystal family: parse, symmetry-expand, density triple-check, geometry, H-bonds. The
figures (crystal_pxrd, crystal_view) consume what this returns; they never recompute, so a
figure can never assert a contact the table doesn't list.

Pure computation + Tier-1 gates (no plotting). gemmi is imported lazily (heavy dep; a
non-CIF job never pulls it). Conventions/constants are pinned in the block below;
the logic was validated against the CIF edge-case set in skill_validation/crystal/ (ELAINM, urea, aspirin, HMT,
ferrocene, flufenamic, lactose, L-alanine, PTU-ellagic + synthetic broken-input CIFs).

    load(cfg)                 parse + block-select + symmetry + atoms (robust to CSD/SHELX quirks)
    expand(struct)            unique (frac,sym,occ) cell positions (global atom x symop dedup)
    densities(struct, cfg)    the signature triple-check + F(000) + element counts + branch gate
    special_positions(struct) atoms on a symmetry element (orbit < n_ops)
    bonds(struct, cfg)        covalent-radii bonds, disorder-aware
    hbonds(struct, cfg)       X-H-normalized D-H...A with symmetry-expanded search + geom_hbond xref
    validate(struct, cfg)     run all Tier-1 gates, assemble the validation-table rows
    table_text / write_csv    emit the validation table (the CORE deliverable)
"""

# ---- pinned constants (sources inline) ------------------------------------------
K_CODATA = 1.66054        # CODATA atomic-mass-constant factor -> assumption-free density (i)
K_CHECKCIF = 1.66042      # checkCIF's legacy DENSD01 factor -> formula density (ii)
NEUTRON_XH = {"C": 1.083, "N": 1.009, "O": 0.983}   # Allen & Bruno 2010 (X-H normalization)
# Bondi (1964) vdW radii (A); core set pinned, others fall back to gemmi.Element.vdw_r
VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
       "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "P": 1.80}
JEFFREY = ((2.2, 2.5, "strong"), (2.5, 3.2, "moderate"), (3.2, 4.0, "weak"))  # D...A ranges (A)
BOND_TOL = 0.40           # A added to covalent-radii sum for bond perception
DISORDER_MIN = 0.90       # A; sub-unity sites closer than this are mutually exclusive (no bond)
GEOM_OUTLIER = 0.25       # A; bond length off the covalent-radii sum by more than this -> flagged
GEOM_CONTACT_MIN = 90.0   # deg; contacts in [this, hbond_angle_min) are "geometric contacts" (listed, not asserted)
STD_WL = {"Cu-Ka": 1.5418, "Cu-Ka1": 1.5406, "Mo": 0.71073,   # standard Kα wavelengths (A)
          "Ag": 0.56086, "Ga": 1.3414, "In": 0.51359}


def _gemmi():
    try:
        import gemmi
        return gemmi
    except ImportError as e:
        raise ImportError(
            "Required package not found: gemmi. Install with:\n"
            "    python -m pip install gemmi\n"
            "(the crystal family needs gemmi for CIF parsing + symmetry).") from e


class Structure:
    """Parsed + cached crystal data. Build with load(cfg)."""
    def __init__(self, **kw):
        self.__dict__.update(kw)
        self._uniq = None
        self._alt = None

    def cart(self, frac):
        p = self.cell.orthogonalize(_gemmi().Fractional(float(frac[0]), float(frac[1]), float(frac[2])))
        return np.array([p.x, p.y, p.z])


# ----------------------------------------------------------------- parse helpers
def _num(s):
    try:
        return float(str(s).split("(")[0])
    except (ValueError, AttributeError, TypeError):
        return None


def _elem_from(type_symbol, label):
    """Element from _atom_site_type_symbol, stripping an oxidation-state suffix
    ('O2-'->O, 'Fe3+'->Fe); fall back to the label ('C1'->C) when type_symbol is
    absent/odd. '' if unrecognisable. (Both quirks seen in the validation set.)"""
    g = _gemmi()
    for src in (type_symbol, label):
        m = re.match(r"\s*([A-Za-z]{1,2})", src or "")
        if not m:
            continue
        a = m.group(1)
        cands = [a[0].upper() + a[1].lower(), a[0].upper()] if len(a) == 2 else [a[0].upper()]
        for c in cands:
            try:
                el = g.Element(c)
                if el.atomic_number > 0:
                    return el.name
            except Exception:
                pass
    return ""


def _parse_formula(s):
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)\s*(\d*\.?\d*)", s or ""):
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    return out


def _vdw(sym):
    if sym in VDW:
        return VDW[sym]
    try:
        r = _gemmi().Element(sym).vdw_r
        return r if r and r > 0 else 1.70
    except Exception:
        return 1.70


def _covsum(s1, s2):
    g = _gemmi()
    return g.Element(s1).covalent_r + g.Element(s2).covalent_r


def _angle(p, q, r):
    v1, v2 = p - q, r - q
    n = np.linalg.norm(v1) * np.linalg.norm(v2)
    if n == 0:
        return float("nan")
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2) / n)))))


# --------------------------------------------------------------------- load
def load(cfg):
    """Parse cfg.cif_path, select the data block, read cell / symmetry / atoms with the
    robustness the validation set demanded: oxidation-state type symbols, missing
    _atom_site_type_symbol or _atom_site_occupancy columns, symop loop OR a bare
    space-group symbol. Gates: file parses, ONE block selected (multi-block FAILs unless
    cfg.cif_block is set), cell + symmetry + atoms present."""
    g = _gemmi()
    if not cfg.cif_path:
        raise ValueError("cfg.cif_path is not set")
    doc = g.cif.read(cfg.cif_path)
    blocks = [b.name for b in doc]
    f = []

    # ---- block selection (spec §4 step 1)
    block = doc[0]
    if cfg.cif_block is not None:
        sel = str(cfg.cif_block)
        try:
            block = doc[int(sel)] if sel.lstrip("-").isdigit() else doc[sel]
        except Exception:
            f.append(("FAIL", f"cif_block '{cfg.cif_block}' not found; blocks={blocks}"))
    elif len(blocks) > 1:
        f.append(("FAIL", f"multi-block CIF ({len(blocks)} blocks: {blocks}); set cfg.cif_block"))

    def sval(tag):
        v = block.find_value(tag)
        return v.strip().strip("'\"") if v else None

    def fval(tag):
        return _num(sval(tag))

    a, b, c = fval("_cell_length_a"), fval("_cell_length_b"), fval("_cell_length_c")
    al, be, ga = fval("_cell_angle_alpha"), fval("_cell_angle_beta"), fval("_cell_angle_gamma")
    if None in (a, b, c, al, be, ga):
        f.append(("FAIL", "cell parameters missing/unparseable"))
        a, b, c = a or 1.0, b or 1.0, c or 1.0
        al, be, ga = al or 90.0, be or 90.0, ga or 90.0
    cell = g.UnitCell(a, b, c, al, be, ga)

    # ---- symmetry operators: explicit loop, else derive from H-M / Hall symbol
    raw = list(block.find_loop("_symmetry_equiv_pos_as_xyz")) or \
        list(block.find_loop("_space_group_symop_operation_xyz"))
    ops = [g.Op(str(o).strip().strip("'\"").replace(" ", "")) for o in raw]
    sg_hm = sval("_symmetry_space_group_name_H-M") or sval("_space_group_name_H-M_alt")
    if not ops:
        hall = sval("_symmetry_space_group_name_Hall") or sval("_space_group_name_Hall")
        sgo = None
        if sg_hm:
            try:
                sgo = g.find_spacegroup_by_name(sg_hm)
            except Exception:
                sgo = None
        if sgo is not None:
            ops = list(sgo.operations())
        elif hall:
            try:
                ops = list(g.symops_from_hall(hall))
            except Exception:
                ops = []
    if not ops:
        ops = [g.Op("x,y,z")]
        f.append(("WARN", "no symmetry operators found; assuming P1"))

    # ---- atoms: column-wise (type_symbol/occupancy optional)
    def col(tag):
        return [str(v) for v in block.find_loop(tag)]

    labels, occs = col("_atom_site_label"), col("_atom_site_occupancy")
    fx, fy, fz = col("_atom_site_fract_x"), col("_atom_site_fract_y"), col("_atom_site_fract_z")
    types, uiso = col("_atom_site_type_symbol"), col("_atom_site_U_iso_or_equiv")
    dis_asm, dis_grp = col("_atom_site_disorder_assembly"), col("_atom_site_disorder_group")
    atoms = []
    for i in range(min(len(labels), len(fx), len(fy), len(fz))):
        x, y, z = _num(fx[i]), _num(fy[i]), _num(fz[i])
        if None in (x, y, z):
            continue
        sym = _elem_from(types[i] if i < len(types) else "", labels[i])
        if not sym:
            continue
        occ = _num(occs[i]) if (i < len(occs) and occs[i] not in ("?", ".")) else 1.0
        atoms.append({"label": labels[i], "sym": sym,
                      "frac": np.array([x, y, z], float), "occ": occ if occ is not None else 1.0,
                      "u_iso": _num(uiso[i]) if i < len(uiso) else None,
                      "dis_asm": dis_asm[i] if i < len(dis_asm) else None,
                      "dis_grp": dis_grp[i] if i < len(dis_grp) else None})
    if not atoms:
        f.append(("FAIL", "no atom sites parsed"))

    # anisotropic displacement parameters (Phase 2): {label: 3x3 U tensor in the CIF basis}
    al_ = col("_atom_site_aniso_label")
    au = {k: col("_atom_site_aniso_U_" + k) for k in ("11", "22", "33", "23", "13", "12")}
    aniso = {}
    for i, lab in enumerate(al_):
        try:
            v = {k: _num(au[k][i]) for k in au}
            if any(x is None for x in v.values()):
                continue
            aniso[lab] = np.array([[v["11"], v["12"], v["13"]],
                                   [v["12"], v["22"], v["23"]],
                                   [v["13"], v["23"], v["33"]]], float)
        except Exception:
            continue

    _resolve(f, cfg, "load")     # raises under cfg.strict on any FAIL

    declared = {"V": fval("_cell_volume"), "Z": fval("_cell_formula_units_Z"),
                "MW": fval("_chemical_formula_weight"),
                "density": fval("_exptl_crystal_density_diffrn"),
                "F000": fval("_exptl_crystal_F_000"),
                "wavelength": fval("_diffrn_radiation_wavelength"),
                "rad_type": sval("_diffrn_radiation_type") or sval("_diffrn_radiation_probe"),
                "temperature": sval("_diffrn_ambient_temperature") or sval("_cell_measurement_temperature"),
                "formula_sum": sval("_chemical_formula_sum"),
                # echo-only refinement metadata (reported, never recomputed — keeps the CIF-only line honest)
                "R_gt": fval("_refine_ls_R_factor_gt"), "R_all": fval("_refine_ls_R_factor_all"),
                "gof": fval("_refine_ls_goodness_of_fit_ref"),
                "reflns_total": fval("_reflns_number_total"), "reflns_gt": fval("_reflns_number_gt"),
                "theta_max": fval("_diffrn_reflns_theta_max"),
                "size_max": fval("_exptl_crystal_size_max"), "size_mid": fval("_exptl_crystal_size_mid"),
                "size_min": fval("_exptl_crystal_size_min")}
    return Structure(name=block.name, blocks=blocks, block=block, cell=cell, ops=ops,
                     atoms=atoms, aniso=aniso, sg_hm=sg_hm, declared=declared,
                     a=a, b=b, c=c, al=al, be=be, ga=ga)


# --------------------------------------------------------------------- expansion
def expand(struct):
    """Every (atom x symop) image in the unit cell, deduplicated GLOBALLY by wrapped
    position+element. Robust to special positions AND symmetry-completed atom lists
    (CSD-issued CIFs list the whole molecule, not the minimal AU). Cached."""
    if struct._uniq is not None:
        return struct._uniq
    uniq = {}
    for at in struct.atoms:
        for op in struct.ops:
            gpos = np.array(op.apply_to_xyz(list(at["frac"]))) % 1.0
            key = (at["sym"], int(round(gpos[0] * 100)) % 100,
                   int(round(gpos[1] * 100)) % 100, int(round(gpos[2] * 100)) % 100)
            uniq.setdefault(key, {"sym": at["sym"], "frac": gpos,
                                  "occ": at["occ"], "label": at["label"]})
    struct._uniq = list(uniq.values())
    return struct._uniq


def special_positions(struct):
    """Atoms whose symmetry orbit is smaller than the number of operators (i.e. they sit
    on a symmetry element). Returns [(label, sym, 'orbit/nops')]."""
    nops = len(struct.ops)
    out = []
    for at in struct.atoms:
        seen = []
        for op in struct.ops:
            gp = np.array(op.apply_to_xyz(list(at["frac"]))) % 1.0
            if not any(np.all(np.abs((gp - s + 0.5) % 1.0 - 0.5) < 3e-3) for s in seen):
                seen.append(gp)
        if len(seen) < nops:
            out.append((at["label"], at["sym"], f"{len(seen)}/{nops}"))
    return out


# --------------------------------------------------------------- ADP / ellipsoids (Phase 2)
def has_adp(struct):
    return bool(getattr(struct, "aniso", None))


def _orth(struct):
    """Orthogonalization matrix M (frac -> cart); columns are the cell vectors a, b, c."""
    return np.column_stack([struct.cart([1, 0, 0]), struct.cart([0, 1, 0]), struct.cart([0, 0, 1])])


def u_cart(struct, U):
    """CIF anisotropic U (the a*_i a*_j basis) -> Cartesian mean-square-displacement tensor
    (Å²): U_cart = (M N) U (M N)^T, with M the orthogonalization matrix and N = diag(a*,b*,c*).
    Validated by U_eq = trace(U_cart)/3 matching the CIF's _atom_site_U_iso_or_equiv."""
    M = _orth(struct)
    rc = struct.cell.reciprocal()
    A = M @ np.diag([rc.a, rc.b, rc.c])
    return A @ np.asarray(U, float) @ A.T


def u_eq(struct, U):
    return float(np.trace(u_cart(struct, U)) / 3.0)


def _adp_scale(p):
    """Ellipsoid radius in σ units enclosing probability p of a 3-D Gaussian (chi-3 quantile).
    50% -> 1.5382 (the ORTEP default)."""
    try:
        from scipy.stats import chi2
        return float(chi2.ppf(p, 3)) ** 0.5
    except Exception:
        return {0.50: 1.5382, 0.30: 1.0972, 0.95: 2.7955, 0.99: 3.3682}.get(round(p, 2), 1.5382)


def adp_axes(struct, U, prob=0.5):
    """(semi_axes[3], V[3x3] eigenvectors as columns) of the displacement ellipsoid at `prob`.
    Returns (None, None) when U is not positive-definite (non-physical ADP -> caller degrades)."""
    w, V = np.linalg.eigh(u_cart(struct, U))
    if np.any(w <= 1e-9):
        return None, None
    return _adp_scale(prob) * np.sqrt(w), V


def op_rot_cart(struct, op):
    """Cartesian rotation of a symmetry operator — to rotate a symmetry image's ADP:
    U_image = R U_cart R^T."""
    M = _orth(struct)
    Rf = np.array(op.rot, float) / op.DEN
    return M @ Rf @ np.linalg.inv(M)


# --------------------------------------------------------------------- density
def densities(struct, cfg):
    """The signature density triple-check + F(000) + element-count cross-check, with the
    full branch tree (spec §5). The expansion-implicated case HARD-FAILS regardless of
    cfg.strict (drawing wrong atoms is a correctness error, not a preference)."""
    g = _gemmi()
    d = struct.declared
    V = d["V"] or struct.cell.volume
    uniq = expand(struct)
    mass = sum(g.Element(u["sym"]).weight * u["occ"] for u in uniq)
    elec = sum(g.Element(u["sym"]).atomic_number * u["occ"] for u in uniq)
    counts = {}
    for u in uniq:
        counts[u["sym"]] = counts.get(u["sym"], 0.0) + u["occ"]

    di = K_CODATA * mass / V if V else float("nan")
    dii = (K_CHECKCIF * d["MW"] * d["Z"] / V) if (d["MW"] and d["Z"] and V) else float("nan")
    diii = d["density"]
    fsum = _parse_formula(d["formula_sum"])
    mw_formula = sum(g.Element(e).weight * n for e, n in fsum.items()) if fsum else None

    def close(x, y, tol=0.03):
        return bool(x and y and not (isinstance(x, float) and math.isnan(x)) and abs(x - y) / y < tol)

    have_ii = not math.isnan(dii)
    i_ok, ii_ok = close(di, diii), close(dii, diii)
    mw_bad = bool(mw_formula and d["MW"] and abs(mw_formula - d["MW"]) / d["MW"] > 0.01)
    counts_match = bool(fsum and d["Z"]) and all(
        abs(counts.get(e, 0.0) - n * d["Z"]) <= max(0.05 * n * d["Z"], 0.15) for e, n in fsum.items())
    atoms_short = bool(fsum and d["Z"]) and \
        sum(counts.values()) < sum(n * d["Z"] for n in fsum.values()) - 0.5

    hard_fail = False
    if diii is None:
        branch, sev = "no declared density (cannot gate)", "INFO"
    elif i_ok and (ii_ok or not have_ii):
        branch, sev = "PASS (all consistent)", "INFO"
    elif have_ii and i_ok and not ii_ok:
        branch, sev = "A-ALERT: declared density/formula inconsistent (ii != iii)", "WARN"
    elif not i_ok and mw_bad and counts_match:
        branch, sev = (f"declared formula_weight wrong ({d['MW']} vs {mw_formula:.2f} "
                       f"from formula_sum); coordinates fine", "WARN")
    elif not i_ok and atoms_short:
        branch, sev = "incomplete atom list (i<iii: H-free / SQUEEZE'd solvent / disorder)", "WARN"
    elif not i_ok and counts_match and have_ii and ii_ok:
        branch, sev, hard_fail = "EXPANSION WRONG (atoms complete yet i != iii) — figure refused", "FAIL", True
    elif not i_ok:
        branch, sev = "(i) != (iii) — inspect (declared metadata or cell)", "WARN"
    else:
        branch, sev = "PASS", "INFO"

    rd = (diii / dii) if (diii and have_ii and dii) else None
    rd_level = None
    if rd is not None:
        rd_level = ("A" if not (0.90 <= rd <= 1.10) else
                    "B" if not (0.95 <= rd <= 1.05) else
                    "C" if not (0.99 <= rd <= 1.01) else None)
    out = {"i": di, "ii": dii, "iii": diii, "branch": branch, "rd": rd, "rd_level": rd_level,
           "mass": mass, "F000_calc": elec, "F000_decl": d["F000"],
           "counts": counts, "mw_formula": mw_formula, "mw_declared": d["MW"]}

    f = [(sev, f"density (i)={di:.4f} (ii)={dii:.4f} (iii)={diii} -> {branch}")]
    if rd_level:
        f.append(("WARN" if rd_level in ("A", "B") else "INFO",
                  f"density RD={rd:.4f} -> checkCIF {rd_level}-alert (DENSD01: declared/formula-DEN)"))
    if d["F000"] and abs(elec - d["F000"]) / d["F000"] > 0.05:
        f.append(("WARN", f"F(000) calc {elec:.0f} vs declared {d['F000']} (>5% off — corroborates the density flag)"))
    if hard_fail:
        for s, m in f:
            print(f"  [{s}] density: {m}")
        raise GateError(f"density: {branch}")     # always raises (correctness)
    _resolve(f, cfg, "density")
    return out


# --------------------------------------------------------------------- disorder
def disorder_alternatives(struct):
    """label -> set of mutually-exclusive disorder-alternative labels. Two sources, unioned:
    (1) NEAR-COINCIDENCE — sub-unity (occ<1) sites of the SAME element within DISORDER_MIN are
        alternatives (handles UNtagged disorder, e.g. ELAINM's water over 3 partial O sites);
    (2) DISORDER TAGS — same _atom_site_disorder_assembly, different _atom_site_disorder_group
        (handles A/B components farther apart than DISORDER_MIN, e.g. flufenamic's CF3).
    Transitively closed and cached, so a donor never H-bonds to its own alternative site."""
    if struct._alt is not None:
        return struct._alt
    alt = {a["label"]: set() for a in struct.atoms}

    def link(la, lb):
        if la != lb:
            alt[la].add(lb); alt[lb].add(la)

    sub = [a for a in struct.atoms if a["occ"] < 0.999]
    for i in range(len(sub)):
        for j in range(i + 1, len(sub)):
            if sub[i]["sym"] != sub[j]["sym"]:
                continue
            if float(np.linalg.norm(struct.cart(sub[i]["frac"]) - struct.cart(sub[j]["frac"]))) < DISORDER_MIN:
                link(sub[i]["label"], sub[j]["label"])

    tagged = [a for a in struct.atoms if a.get("dis_grp") not in (None, ".", "", "?")]
    for i in range(len(tagged)):
        for j in range(i + 1, len(tagged)):
            if tagged[i].get("dis_asm") == tagged[j].get("dis_asm") \
                    and tagged[i]["dis_grp"] != tagged[j]["dis_grp"]:
                link(tagged[i]["label"], tagged[j]["label"])

    changed = True                                        # transitive closure (chain -> full cluster)
    while changed:
        changed = False
        for k in alt:
            new = set().union(*(alt[m] for m in alt[k])) if alt[k] else set()
            new.discard(k)
            if not new <= alt[k]:
                alt[k] |= new; changed = True
    struct._alt = alt
    return alt


# --------------------------------------------------------------------- bonds
def bonds(struct, cfg):
    """Covalent bonds within the cell, DISORDER-AWARE: two sub-unity sites closer than
    DISORDER_MIN are alternative positions, not bonded (else naive perception invents
    0.5 A 'bonds'). Returns list of (i, j, dist) index pairs into expand(struct)."""
    atoms = expand(struct)
    alt = disorder_alternatives(struct)
    xyz = [struct.cart(u["frac"]) for u in atoms]
    out = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            dd = float(np.linalg.norm(xyz[i] - xyz[j]))
            if dd <= 0.4 or dd > _covsum(atoms[i]["sym"], atoms[j]["sym"]) + BOND_TOL:
                continue
            if atoms[i]["occ"] < 1 and atoms[j]["occ"] < 1 and dd < DISORDER_MIN:
                continue                                  # mutually-exclusive disorder (near-coincident)
            if atoms[j]["label"] in alt.get(atoms[i]["label"], set()):
                continue                                  # disorder alternative (incl. group-tagged)
            out.append((i, j, round(dd, 3)))
    return out


def _sym_code(struct, label, frac):
    """checkCIF-style symmetry code 'n_pqr' for an atom image at `frac` generated from the
    asym-unit atom `label`: n = 1-based symop index, pqr = 5 + lattice translation.
    '.' = identity in the reference cell; '?' if not recoverable."""
    asym = next((a for a in struct.atoms if a["label"] == label), None)
    if asym is None:
        return "?"
    fr = np.asarray(frac, float)
    for n, op in enumerate(struct.ops, 1):
        base = np.array(op.apply_to_xyz(list(asym["frac"])), float)
        t = np.round(fr - base)
        if np.allclose(base + t, fr, atol=2e-2):
            if n == 1 and np.allclose(t, 0):
                return "."
            return "%d_%d%d%d" % (n, 5 + int(t[0]), 5 + int(t[1]), 5 + int(t[2]))
    return "?"


# --------------------------------------------------------------------- geometry (Block B)
def geometry(struct, cfg):
    """Block B: selected bond lengths flagged against the covalent-radii sum (Cordero/gemmi
    radii ±0.25 A). This is a COARSE sanity bound, NOT a Mogul/CSD percentile (Mogul is
    licensed/unavailable — don't overclaim distributional rigour). Disorder-aware; bonds are
    deduped to unique (label-pair, length) types over the cell + nearest neighbours."""
    cell_atoms = [{**u, "xyz": struct.cart(u["frac"])} for u in expand(struct)]
    alt = disorder_alternatives(struct)
    origin = struct.cart(np.zeros(3))
    sup = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                shift = struct.cart(np.array([di, dj, dk], float)) - origin
                for a in cell_atoms:
                    sup.append({"sym": a["sym"], "label": a["label"], "occ": a["occ"],
                                "xyz": a["xyz"] + shift})
    rows, seen, n_out = [], set(), 0
    for a in cell_atoms:
        for b in sup:
            dd = float(np.linalg.norm(a["xyz"] - b["xyz"]))
            if dd <= 0.4:
                continue
            cs = _covsum(a["sym"], b["sym"])
            if dd > cs + BOND_TOL:
                continue
            if a["occ"] < 1 and b["occ"] < 1 and dd < DISORDER_MIN:
                continue                                          # mutually-exclusive disorder (near-coincident)
            if b["label"] in alt.get(a["label"], set()):
                continue                                          # disorder alternative (incl. group-tagged)
            key = tuple(sorted([a["label"], b["label"]]) + [round(dd, 2)])
            if key in seen:
                continue
            seen.add(key)
            outlier = abs(dd - cs) > GEOM_OUTLIER
            n_out += int(outlier)
            rows.append({"a": a["label"], "b": b["label"], "len": round(dd, 3),
                         "covsum": round(cs, 3), "outlier": outlier})
    f = [("INFO", f"Block B: {len(rows)} unique bonds, {n_out} outside the covalent envelope "
                  f"(±{GEOM_OUTLIER} A; coarse bound, not Mogul)")]
    for bd in rows:
        if bd["outlier"]:
            f.append(("WARN", f"bond {bd['a']}-{bd['b']} {bd['len']} A vs covalent sum "
                              f"{bd['covsum']} A — outlier"))
    _resolve(f, cfg, "geometry")
    return {"bonds": rows, "n_outliers": n_out}


def _wavelength_check(struct):
    """Flag the declared wavelength unless it is within 5e-4 A of a standard Kα, or the source
    is neutron/synchrotron/electron (where any λ is legitimate)."""
    wl = struct.declared["wavelength"]
    rt = (struct.declared["rad_type"] or "").lower()
    if wl is None:
        return ("INFO", "wavelength not declared")
    if any(k in rt for k in ("neutron", "synchrotron", "electron")):
        return ("INFO", f"wavelength {wl} A ({rt}) — non-Kα source, accepted")
    near = [k for k, v in STD_WL.items() if abs(wl - v) < 5e-4]
    if near:
        return ("INFO", f"wavelength {wl} A matches {near[0]}")
    return ("WARN", f"wavelength {wl} A not within 5e-4 of a standard Kα "
                    f"(type='{rt or '?'}') — verify")


def _sg_check(struct):
    """Cross-check the symmetry operators against the declared space-group name: the op count
    must equal the named group's order (catches an incomplete symop list / wrong setting), and
    where derivable, the group recovered FROM the ops must match the declared H-M symbol."""
    g = _gemmi()
    decl = struct.sg_hm
    named = None
    if decl:
        try:
            named = g.find_spacegroup_by_name(decl)
        except Exception:
            named = None
    if named is not None:
        order = len(list(named.operations()))
        if order != len(struct.ops):
            return ("WARN", f"space group '{decl}' has order {order} but {len(struct.ops)} "
                            f"symops are present — incomplete/extra symmetry?")
    try:
        from_ops = g.find_spacegroup_by_ops(g.GroupOps(list(struct.ops)))
        if from_ops is not None and decl and \
                from_ops.hm.replace(" ", "") != decl.replace(" ", ""):
            return ("INFO", f"space group from ops = '{from_ops.hm}' vs declared '{decl}' "
                            f"(setting/naming difference — verify)")
    except Exception:
        pass
    return ("INFO", f"space group {decl}: {len(struct.ops)} symops, consistent")


# --------------------------------------------------------------------- H-bonds
def hbonds(struct, cfg):
    """X-H-normalized D-H...A detection with a symmetry-expanded (3x3x3) acceptor search,
    classified by D...A (Jeffrey). Donor/acceptor sets are separate (cfg). Sub-floor
    contacts are demoted, not reported as H-bonds. Cross-checks the CIF's _geom_hbond
    loop on D...A (normalization-invariant) where present."""
    g = _gemmi()
    donors = set(cfg.hbond_donors) | ({"C"} if cfg.hbond_weak else set())
    acceptors = set(cfg.hbond_acceptors)
    floor = cfg.hbond_angle_min
    alt = disorder_alternatives(struct)
    dis_pairs = set()
    f = []
    if cfg.hbond_weak and floor >= 120.0:
        f.append(("INFO", "weak donors enabled but angle floor still 120 deg — most weak "
                          "H-bonds are bent; consider lowering hbond_angle_min toward 90"))

    cell_atoms = [{**u, "xyz": struct.cart(u["frac"])} for u in expand(struct)]
    Hs = [a for a in cell_atoms if a["sym"] == "H"]

    # 3x3x3 supercell of HEAVY atoms: the donor search must see molecules that straddle
    # the cell boundary (not just the [0,1) image, or O-H donors on a boundary-spanning
    # molecule are missed) + the acceptor subset for the contact search.
    sup_heavy, sup_acc = [], []
    origin = struct.cart(np.zeros(3))
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                shift = struct.cart(np.array([di, dj, dk], float)) - origin
                for a in cell_atoms:
                    if a["sym"] == "H":
                        continue
                    rec = {"sym": a["sym"], "label": a["label"], "cell": (di, dj, dk),
                           "frac": a["frac"] + np.array([di, dj, dk], float), "xyz": a["xyz"] + shift}
                    sup_heavy.append(rec)
                    if a["sym"] in acceptors:
                        sup_acc.append(rec)

    results, seen, contacts, seen_c = [], set(), [], set()
    for h in Hs:
        if not sup_heavy:
            break
        D = min(sup_heavy, key=lambda a: np.linalg.norm(a["xyz"] - h["xyz"]))
        dDH = float(np.linalg.norm(D["xyz"] - h["xyz"]))
        if dDH > _covsum(D["sym"], "H") + BOND_TOL or D["sym"] not in donors:
            continue
        hpos = h["xyz"]
        if cfg.xh_normalize and D["sym"] in NEUTRON_XH and dDH > 0:
            hpos = D["xyz"] + (h["xyz"] - D["xyz"]) / dDH * NEUTRON_XH[D["sym"]]
        altD = alt.get(D["label"], set())
        for A in sup_acc:
            DA = float(np.linalg.norm(A["xyz"] - D["xyz"]))
            if DA < 0.4 or DA > 4.0:
                continue
            if A["label"] in altD:                        # acceptor is a disorder-alternative of the donor
                dis_pairs.add((D["label"], A["label"]))
                continue
            HA = float(np.linalg.norm(A["xyz"] - hpos))
            if HA > _vdw("H") + _vdw(A["sym"]):
                continue
            ang = _angle(D["xyz"], hpos, A["xyz"])
            rec = {"D": D["label"], "Dsym": D["sym"], "H": h["label"],
                   "A": A["label"], "Asym": A["sym"], "DH": round(dDH, 3),
                   "HA": round(HA, 3), "DA": round(DA, 3), "angle": round(ang, 1),
                   "class": next((k for lo, hi, k in JEFFREY if lo <= DA < hi), "long"),
                   "sym": _sym_code(struct, A["label"], A["frac"]), "cell": A["cell"]}
            if ang < floor:                               # below the floor: geometric contact, not an asserted H-bond
                if ang >= GEOM_CONTACT_MIN:
                    ck = (D["label"], h["label"], A["label"], round(DA, 2))
                    if ck not in seen_c:
                        seen_c.add(ck); rec["class"] = "contact"; contacts.append(rec)
                continue
            key = (D["label"], h["label"], A["label"], round(DA, 2))
            if key in seen:
                continue
            seen.add(key)
            results.append(rec)

    # cross-check the CIF's own _geom_hbond loop (D...A is normalization-invariant)
    xref = {"rows": 0, "matched": 0}
    try:
        D_l = [str(v) for v in struct.block.find_loop("_geom_hbond_atom_site_label_D")]
        A_l = [str(v) for v in struct.block.find_loop("_geom_hbond_atom_site_label_A")]
        DA_l = [_num(v) for v in struct.block.find_loop("_geom_hbond_distance_DA")]
        found_DA = [r["DA"] for r in results]
        for i in range(min(len(D_l), len(A_l), len(DA_l))):
            xref["rows"] += 1
            if DA_l[i] is not None and any(abs(DA_l[i] - x) < 0.03 for x in found_DA):
                xref["matched"] += 1
        if xref["rows"]:
            f.append(("INFO", f"_geom_hbond cross-check: matched {xref['matched']}/{xref['rows']} "
                              f"deposited bonds on D...A (others may be below the {floor:.0f} deg floor)"))
    except Exception:
        pass

    if dis_pairs:
        f.append(("INFO", f"suppressed {len(dis_pairs)} donor/own-disorder-alternative acceptor "
                          f"pair(s): {sorted(dis_pairs)}"))
    f.insert(0, ("INFO", f"{len(results)} H-bonds + {len(contacts)} geometric contacts "
                         f"(donors={sorted(donors)}, acceptors={sorted(acceptors)}, floor={floor:.0f} deg)"))
    _resolve(f, cfg, "hbonds")
    return {"bonds": results, "contacts": contacts, "xref": xref}


# --------------------------------------------------------------------- validate + table
def validate(struct, cfg):
    """Run the Tier-1 gates and assemble the validation table. Returns a report dict
    {density, special, hbonds, volume, wavelength, ...} — the CORE deliverable that the
    figures consume."""
    g = _gemmi()
    # cell volume recompute (assumption-free) vs declared
    ca, cb, cg = (math.cos(math.radians(t)) for t in (struct.al, struct.be, struct.ga))
    Vcalc = struct.a * struct.b * struct.c * math.sqrt(
        max(0.0, 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg))
    Vdecl = struct.declared["V"]
    fV = []
    if Vdecl and abs(Vcalc - Vdecl) / Vdecl > 0.01:
        fV.append(("WARN", f"cell volume calc {Vcalc:.2f} vs declared {Vdecl} (>1% off)"))
    _resolve(fV, cfg, "volume")

    dens = densities(struct, cfg)
    spec = special_positions(struct)
    geom = geometry(struct, cfg)
    hb = hbonds(struct, cfg)
    wlchk = _wavelength_check(struct); _resolve([wlchk], cfg, "wavelength")
    sgchk = _sg_check(struct); _resolve([sgchk], cfg, "spacegroup")
    elems = sorted({a["sym"] for a in struct.atoms})
    missing = [e for e in elems if e not in VDW]
    fE = []
    if missing:
        fE.append(("INFO", f"elements outside the pinned vdW set {sorted(VDW)}: {missing} "
                           f"(using gemmi vdW radii as fallback)"))
    _resolve(fE, cfg, "elements")

    return {"name": struct.name, "blocks": struct.blocks, "sg_hm": struct.sg_hm,
            "nops": len(struct.ops), "elems": elems, "missing_vdw": missing,
            "V_calc": Vcalc, "V_decl": Vdecl, "density": dens, "special": spec,
            "geometry": geom, "wl_check": wlchk, "sg_check": sgchk, "floor": cfg.hbond_angle_min,
            "hbonds": hb["bonds"], "contacts": hb["contacts"], "hbond_xref": hb["xref"],
            "declared": struct.declared}


def table_text(report):
    """Human-readable validation table (the text deliverable)."""
    d, dn = report["density"], report["declared"]
    g = report["geometry"]
    rd = (f"  RD            {d['rd']:.4f}" +
          (f"  -> checkCIF {d['rd_level']}-alert" if d.get("rd_level") else "  (ok)")) if d.get("rd") else None
    size = "x".join(f"{dn[k]}" for k in ("size_max", "size_mid", "size_min")) if dn.get("size_max") else "?"
    L = [f"VALIDATION TABLE — {report['name']}  ({report['sg_hm']}, {report['nops']} ops)",
         "-" * 64,
         "Block A — crystal data",
         f"  cell volume   calc {report['V_calc']:.2f}   decl {report['V_decl']}",
         f"  density (i)   {d['i']:.4f}  (expanded-cell, assumption-free, CODATA)",
         f"  density (ii)  {d['ii']:.4f}  (checkCIF formula, 1.66042)",
         f"  density (iii) {d['iii']}  (declared)"]
    if rd:
        L.append(rd)
    L += [f"     -> {d['branch']}",
          f"  F(000)        calc {d['F000_calc']:.0f}   decl {d['F000_decl']}",
          f"  formula wt    from formula_sum {d['mw_formula']}   decl {d['mw_declared']}",
          f"  elements      {report['elems']}" + (f"  (outside vdW table: {report['missing_vdw']})"
                                                  if report["missing_vdw"] else ""),
          f"  wavelength    {dn['wavelength']}  ({dn['rad_type']})  [{report['wl_check'][1]}]",
          f"  space group   {report['sg_check'][1]}",
          f"  special pos   {report['special'] or 'none'}",
          f"  echo (not recomputed)  T={dn['temperature']}  R(gt)={dn.get('R_gt')}  GoF={dn.get('gof')}  "
          f"reflns={dn.get('reflns_total')}  th_max={dn.get('theta_max')}  size={size}",
          "Block B — geometry (bond lengths vs covalent-radii sum +/-0.25 A; coarse bound, not Mogul)",
          f"  {g['n_outliers']} outlier(s) of {len(g['bonds'])} unique bonds"]
    for bd in g["bonds"]:
        if bd["outlier"]:
            L.append(f"    {bd['a']}-{bd['b']}  {bd['len']:.3f} A  (covalent sum {bd['covsum']:.3f})  OUTLIER")
    L.append("Block C — H-bonds (D-H...A, X-H normalized; D...A class by Jeffrey; sym = operator on A)")
    if report["hbonds"]:
        L.append("  {:<6} {:<6} {:<6} {:>6} {:>6} {:>6} {:>6}  {:<9} {}".format(
            "D", "H", "A", "D-H", "H..A", "D..A", "ang", "class", "sym"))
        for b in report["hbonds"]:
            L.append("  {:<6} {:<6} {:<6} {:6.3f} {:6.3f} {:6.3f} {:6.1f}  {:<9} {}".format(
                b["D"], b["H"], b["A"], b["DH"], b["HA"], b["DA"], b["angle"], b["class"], b["sym"]))
        x = report["hbond_xref"]
        if x["rows"]:
            L.append(f"  (_geom_hbond cross-check: {x['matched']}/{x['rows']} matched on D..A)")
    else:
        L.append("  none detected")
    if report["contacts"]:
        L.append(f"  geometric contacts (below the {report['floor']:.0f} deg floor — listed, not asserted):")
        for c in report["contacts"]:
            L.append("    {:<6} {:<6} {:<6} D..A {:6.3f}  ang {:5.1f}  {}".format(
                c["D"], c["H"], c["A"], c["DA"], c["angle"], c["sym"]))
    return "\n".join(L)


def write_csv(report, path):
    """Write the H-bond table to CSV (the SI-ready artefact)."""
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["D", "H", "A", "D-H/A", "H..A/A", "D..A/A", "angle/deg", "class", "sym"])
        for b in report["hbonds"] + report["contacts"]:
            w.writerow([b["D"], b["H"], b["A"], b["DH"], b["HA"], b["DA"], b["angle"], b["class"], b["sym"]])
    return path


# ======================================================================
# from crystal_pxrd.py
# ======================================================================
"""
py  -  calculated PXRD from a CIF, handed to the spectra module (pxrd domain).

Calculated PXRD is a new *source*, not a new plotting family: this computes (2theta,
intensity) with Dans_Diffraction and plots it through the house pxrd conventions. The two
verified gotchas are handled here: the (000) direct beam is kept OFF-GRID (min_twotheta via
setup_scatter), and peak_width is in Q (A^-1), NOT degrees. Cross-checked against pymatgen
(or a Bragg-law fallback). Validated in skill_validation/dans_fix.py (ELAINM: 47 peaks,
Dans vs pymatgen top peak 0.011 deg).

    calc_pattern(struct, cfg)             -> {two_theta, intensity, wavelength, peaks}
    plot(struct, cfg, experimental=None)  -> (fig, ax, pattern)   via spectra pxrd axis
"""


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
    _resolve(f, cfg, "pxrd")
    return {"two_theta": tt, "intensity": inten, "wavelength": wl, "peaks": peaks}


def reflection_list(struct, cfg):
    """The powder reflection list `[(two_theta_deg, intensity, (h,k,l)), ...]` from
    Dans_Diffraction — REAL structure-factor intensities with multiplicity + Lorentz-polarization
    already applied (Dans's `powder()` 3rd return: columns h,k,l,2theta,intensity, grouped by
    min_overlap). This is the correct input to `simulate_pattern` (feed it THESE, not
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
    (`reflection_list`) fed to `simulate_pattern`, so the Kalpha2 doublet
    (`cfg.pxrd_kalpha2`), preferred orientation (`cfg.pxrd_po_axis`/`pxrd_march_r`) and Caglioti
    broadening (`cfg.pxrd_caglioti`) are applied ON TOP of physically-correct intensities. The
    alpha1 line is the CIF/cfg wavelength (not assumed Cu). Returns
    {two_theta, intensity, reflections, wavelength}. Use for overlaying calc vs a measured Cu-Ka
    scan; for positions-only phase ID, `calc_pattern` is enough.

    Coherence gotcha: the Kalpha2 default is **Cu** Ka2. To simulate a Cu lab scan from a
    Mo-refined CIF, set `cfg.pxrd_wavelength=1.540598` so alpha1 is Cu too (else you pair a Cu
    Ka2 onto a Mo alpha1); for a non-Cu anode set `cfg.pxrd_wavelength2` to its Ka2."""
    wl, _ = _wavelength(struct, cfg)
    refl = reflection_list(struct, cfg)
    if x_grid is None:
        x_grid = np.linspace(cfg.pxrd_two_theta_min, cfg.pxrd_two_theta_max, npoints)
    Gs = reciprocal_metric(struct.a, struct.b, struct.c, struct.al, struct.be, struct.ga)
    lam2 = cfg.pxrd_wavelength2 or CU_KA2
    y = simulate_pattern(refl, x_grid, cfg, lam1=wl, lam2=lam2, Gs=Gs)
    return {"two_theta": np.asarray(x_grid, float), "intensity": y, "reflections": refl,
            "wavelength": wl}


def plot(struct, cfg, experimental=None, pattern=None):
    """Plot the calculated pattern through the house pxrd conventions (normal 2theta axis,
    intensity). `experimental`=(2theta, I) overlays a measured trace for phase ID. We don't
    run check_trace here — the data is our own validated calc, not an ingest, and its
    sharp Bragg peaks would trip the FTIR dead-pixel spike heuristic."""
    if pattern is None:
        pattern = calc_pattern(struct, cfg)
    pcfg = replace(cfg, domain="pxrd")
    apply_style(pcfg)
    fig, ax = figure(pcfg)
    ax.plot(pattern["two_theta"], pattern["intensity"], label="calculated")
    if experimental is not None:
        ex, ey = np.asarray(experimental[0], float), np.asarray(experimental[1], float)
        if ey.size and ey.max() > 0:
            ey = ey / ey.max() * 100.0
        ax.plot(ex, ey, label="experimental")
    _apply_axis(ax, pcfg)
    if experimental is not None:
        ax.legend(loc="best")
    finalize_figure(fig)
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
    the traces) with the **right-margin per-trace key** (`edge_labels`), never labels over
    the data. `entries` = list of (label, cif_path); each pattern is computed at the house pxrd
    conventions and normalized to 100. `experimental` = (2theta, I) is drawn (grey) at the base.
    Returns (fig, ax, patterns) — patterns is a list of (label, pattern dict)."""
    pcfg = replace(cfg, domain="pxrd")
    apply_style(pcfg)
    fig, ax = figure(pcfg)
    pal = _palette(pcfg)
    traces, patterns = [], []                       # (label, x, y, colour), bottom -> top
    if experimental is not None:
        ex, ey = np.asarray(experimental[0], float), np.asarray(experimental[1], float)
        if ey.size and ey.max() > 0:
            ey = ey / ey.max() * 100.0
        traces.append(("experimental", ex, ey, "0.3"))
    for i, (label, path) in enumerate(entries):
        ecfg = replace(pcfg, cif_path=path)
        s = load(ecfg)
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
    _apply_axis(ax, pcfg)
    ax.set_yticks([])                               # offsets are arbitrary; hide the y scale
    ax.set_ylabel("Intensity (normalized, offset)")
    edge_labels(ax, edges)                  # the house key: right margin, per trace
    return fig, ax, patterns


# ======================================================================
# from crystal_view.py
# ======================================================================
"""
py  -  deterministic 3D structure render. A *view* of the engine's validated
computation, never an independent artifact.

The make-or-break is camera orientation; the invariant is REPRODUCIBILITY (the camera lives
in cfg as numbers), not PCA specifically. PCA face-on is the auto default (+ the degeneracy
gates); custom / vector / axis angles are first-class because they're equally reproducible.

    complete_molecules(struct, cfg)   grow whole molecules across symmetry (don't orient a fragment)
    orient(struct, cfg, atoms)        -> (R, (elev,azim,roll), findings)  the deterministic camera
    render(struct, cfg)               -> (fig, ax)   element-coloured, bonds + dashed H-bonds + labels
"""

# CPK / element colours (recognizable); labels carry identity so colour is redundant.
CPK = {"H": "#D0D0D0", "C": "#404040", "N": "#3050F8", "O": "#FF2010", "S": "#E0C020",
       "F": "#80D060", "Cl": "#20C020", "Br": "#A62929", "I": "#8000C0", "P": "#FF8000",
       "Fe": "#E06633"}


def _color(sym):
    return CPK.get(sym, "#B000B0")


def _bond_color(sym):
    # a bond half needs a VISIBLE colour: H's CPK near-white (#D0D0D0) disappears on a
    # white page and makes O-H/N-H bonds look "half broken" -> use a mid grey for H bonds
    # (the atom sphere stays light CPK with a black outline).
    return "#9A9A9A" if sym == "H" else _color(sym)


def _size(sym):
    try:
        r = _gemmi().Element(sym).covalent_r or 0.7
    except Exception:
        r = 0.7
    return float(40.0 * r * r)             # area ~ radius^2, readable at print size


# ----------------------------------------------------------- molecule completion
def complete_molecules(struct, cfg):
    """Grow the asymmetric-unit fragments into whole molecules by following covalent bonds
    across a 3x3x3 neighbourhood (so an inversion-centre half-molecule is completed, and a
    boundary-straddling molecule is made whole) — PCA must orient a real molecule, not an
    AU fragment. Disorder-aware (mutually-exclusive partial sites aren't bonded)."""
    cell_atoms = [{**u, "xyz": struct.cart(u["frac"])} for u in expand(struct)]
    alt = disorder_alternatives(struct)
    origin = struct.cart(np.zeros(3))
    sup = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                shift = struct.cart(np.array([di, dj, dk], float)) - origin
                for a in cell_atoms:
                    sup.append({"sym": a["sym"], "label": a["label"], "occ": a["occ"],
                                "xyz": a["xyz"] + shift})

    def k(p):
        return (int(round(p[0] * 50)), int(round(p[1] * 50)), int(round(p[2] * 50)))

    chosen, have = [], set()
    for at in struct.atoms:                       # seeds: one image of each AU atom
        rec = {"sym": at["sym"], "label": at["label"], "occ": at["occ"],
               "xyz": struct.cart(at["frac"])}
        if k(rec["xyz"]) not in have:
            have.add(k(rec["xyz"]))
            chosen.append(rec)
    frontier = list(chosen)
    while frontier:
        nxt = []
        for a in frontier:
            for s in sup:
                kk = k(s["xyz"])
                if kk in have:
                    continue
                d = float(np.linalg.norm(a["xyz"] - s["xyz"]))
                if not (0.4 < d <= _covsum(a["sym"], s["sym"]) + BOND_TOL):
                    continue
                if a["occ"] < 1 and s["occ"] < 1 and d < DISORDER_MIN:
                    continue
                if s["label"] in alt.get(a["label"], set()):       # disorder-alternative -> not bonded
                    continue
                have.add(kk)
                chosen.append(s)
                nxt.append(s)
        frontier = nxt
    return chosen


def _hide_ch(atoms):
    heavy = [a for a in atoms if a["sym"] != "H"]
    out = []
    for a in atoms:
        if a["sym"] == "H" and heavy:
            D = min(heavy, key=lambda h: np.linalg.norm(h["xyz"] - a["xyz"]))
            if D["sym"] == "C" and np.linalg.norm(D["xyz"] - a["xyz"]) <= _covsum("C", "H") + BOND_TOL:
                continue
        out.append(a)
    return out


def _bond_pairs(atoms, alt=None):
    xyz = [a["xyz"] for a in atoms]
    out = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if 0.4 < d <= _covsum(atoms[i]["sym"], atoms[j]["sym"]) + BOND_TOL:
                if atoms[i]["occ"] < 1 and atoms[j]["occ"] < 1 and d < DISORDER_MIN:
                    continue
                if alt and atoms[j]["label"] in alt.get(atoms[i]["label"], set()):
                    continue                                  # disorder alternative -> not bonded
                out.append((i, j))
    return out


def _hbond_pairs(atoms, cfg):
    """(donor_idx, acceptor_idx) for dashed lines, by the engine's criteria within the
    rendered cluster (same rule as hbonds -> the figure can't assert a bond
    the table wouldn't)."""
    donors = set(cfg.hbond_donors) | ({"C"} if cfg.hbond_weak else set())
    acc = set(cfg.hbond_acceptors)
    heavy = [a for a in atoms if a["sym"] != "H"]
    out = []
    for hi, h in enumerate(atoms):
        if h["sym"] != "H" or not heavy:
            continue
        D = min(heavy, key=lambda a: np.linalg.norm(a["xyz"] - h["xyz"]))
        dDH = float(np.linalg.norm(D["xyz"] - h["xyz"]))
        if dDH > _covsum(D["sym"], "H") + BOND_TOL or D["sym"] not in donors:
            continue
        hpos = h["xyz"]
        if cfg.xh_normalize and D["sym"] in NEUTRON_XH and dDH > 0:
            hpos = D["xyz"] + (h["xyz"] - D["xyz"]) / dDH * NEUTRON_XH[D["sym"]]
        for ai, a in enumerate(atoms):
            if a["sym"] not in acc:
                continue
            DA = float(np.linalg.norm(a["xyz"] - D["xyz"]))
            if DA < 0.4 or DA > 4.0:
                continue
            if np.linalg.norm(a["xyz"] - hpos) > _vdw("H") + _vdw(a["sym"]):
                continue
            if _angle(D["xyz"], hpos, a["xyz"]) < cfg.hbond_angle_min:
                continue
            out.append((hi, ai))     # dashed line starts at the HYDROGEN (H...A), not the donor O/N
    return out


def _hbond_label_set(atoms, cfg):
    """Atom indices to label under view_label_atoms='hbond': the DONOR heavy atom and the ACCEPTOR
    of each H-bond — i.e. only the atoms that make the synthon (the bridging H stays unlabelled).
    Empty set if the rendered cluster has no H-bonds."""
    idx = set()
    for hi, ai in _hbond_pairs(atoms, cfg):
        idx.add(ai)
        h = atoms[hi]["xyz"]
        donor = min((j for j, a in enumerate(atoms) if a["sym"] != "H"),
                    key=lambda j: float(np.linalg.norm(atoms[j]["xyz"] - h)), default=None)
        if donor is not None:
            idx.add(donor)
    return idx


# ------------------------------------------------------------------ orientation
def _view_axis(struct, cfg, P):
    """The line-of-sight axis (unit, cartesian) for the chosen mode, + PCA eigenvalues for
    the shape-degeneracy flag."""
    Pc = P - P.mean(0)
    evals, evecs = np.linalg.eigh(Pc.T @ Pc)          # ascending
    mode = cfg.view_orientation
    if mode == "pca":
        return evecs[:, 0], evals
    if mode in ("axis_a", "axis_b", "axis_c"):
        uvw = {"axis_a": (1, 0, 0), "axis_b": (0, 1, 0), "axis_c": (0, 0, 1)}[mode]
    elif mode == "vector":
        uvw = cfg.view_vector or (0, 0, 1)
    else:
        return evecs[:, 0], evals                      # custom handled separately
    d = struct.cart(np.array(uvw, float)) - struct.cart(np.zeros(3))
    n = np.linalg.norm(d)
    return (d / n if n else evecs[:, 0]), evals


def orient(struct, cfg, atoms):
    """Return (R, (elev,azim,roll), findings). For pca/axis/vector R rotates coords into a
    view frame (z = line of sight, x = largest in-plane spread) and the camera looks straight
    down z; the in-plane roll aligns the H-bond network (or long axis) horizontal. For
    'custom' R is None and the given angles are used. view_tilt nudges any base. det(R)=+1
    keeps a proper rotation (no mirror -> correct enantiomorph)."""
    f = []
    heavy = [a for a in atoms if a["sym"] != "H"]
    P = np.array([a["xyz"] for a in (heavy or atoms)])
    de, da, dr = cfg.view_tilt

    if cfg.view_orientation == "custom":
        ang = cfg.view_angles or (90.0, -90.0, 0.0)
        f.append(("INFO", f"custom camera elev/azim/roll={ang} (+tilt {cfg.view_tilt})"))
        return None, (ang[0] + de, ang[1] + da, ang[2] + dr), f

    z, evals = _view_axis(struct, cfg, P)
    z = z / (np.linalg.norm(z) or 1.0)
    # in-plane axes: PCA of coords projected off z (x = largest in-plane spread)
    Pc = P - P.mean(0)
    proj = Pc - np.outer(Pc @ z, z)
    ev2, evec2 = np.linalg.eigh(proj.T @ proj)
    x = evec2[:, -1] - (evec2[:, -1] @ z) * z
    x = x / (np.linalg.norm(x) or 1.0)
    y = np.cross(z, x)

    # roll: align the mean in-plane D->A H-bond vector horizontal (else long axis horizontal)
    theta = 0.0
    if cfg.view_roll_objective == "hbond":
        hb = _hbond_pairs(atoms, cfg)
        vecs = [atoms[a]["xyz"] - atoms[d]["xyz"] for d, a in hb]
        if vecs:
            inplane = np.array([[v @ x, v @ y] for v in vecs])
            mean = inplane.sum(0)
            if np.linalg.norm(mean) > 0.3 * np.abs(inplane).sum(0).max():
                theta = -math.atan2(mean[1], mean[0])
            else:
                f.append(("INFO", "H-bond roll under-determined (vectors cancel/out-of-plane) "
                                  "-> long-axis horizontal"))
    if theta:
        xr = math.cos(theta) * x + math.sin(theta) * y
        yr = -math.sin(theta) * x + math.cos(theta) * y
        x, y = xr, yr

    R = np.array([x, y, z])
    if np.linalg.det(R) < 0:                     # keep a proper rotation (no mirror)
        R = np.array([x, -y, z])

    # view_tilt is folded INTO R, not added to the returned angles.
    # Why: the PyVista renderer (the default) builds its camera from R alone and DISCARDS the
    # angles - `R, _, findings = orient(...)`. Adding the tilt to the angles therefore did
    # nothing at all on the default backend, silently, while appearing to work in the matplotlib
    # fallback. Folding it into R makes both backends honour it, and the angles below stay at the
    # base (90, -90, 0) so matplotlib does not apply the same nudge twice.
    if (de, da, dr) != (0.0, 0.0, 0.0):
        ce, se = math.cos(math.radians(de)), math.sin(math.radians(de))
        ca, sa = math.cos(math.radians(da)), math.sin(math.radians(da))
        cr, sr = math.cos(math.radians(dr)), math.sin(math.radians(dr))
        Rz = np.array([[cr, -sr, 0.0], [sr, cr, 0.0], [0.0, 0.0, 1.0]])   # roll, about the sight line
        Rx = np.array([[1.0, 0.0, 0.0], [0.0, ce, -se], [0.0, se, ce]])   # elevation
        Ry = np.array([[ca, 0.0, sa], [0.0, 1.0, 0.0], [-sa, 0.0, ca]])   # azimuth
        R = Rz @ Rx @ Ry @ R                     # view-frame nudge, applied after the base

    # shape-degeneracy flag (only meaningful for the auto PCA view)
    if cfg.view_orientation == "pca" and evals[2] > 0:
        w2w1 = evals[1] / evals[2]
        if w2w1 >= 0.85:
            f.append(("WARN", f"orientation degenerate (w2/w1={w2w1:.2f} >= 0.85): no unique "
                              f"plane; PCA view is a default, not a derived answer"))
    f.append(("INFO", f"orientation={cfg.view_orientation} (+tilt {cfg.view_tilt})"))
    return R, (90.0, -90.0, 0.0), f                # tilt is already in R (see above)


# ----------------------------------------------------------------------- render
def _render_matplotlib(struct, cfg, atoms=None, cell_box=False):
    """Zero-dependency FALLBACK renderer (matplotlib 3D). matplotlib has no depth buffer,
    so atom/bond occlusion at vertices is imperfect (small white wedges) — prefer the
    pyvista backend for publication output. Returns (fig, ax)."""
    if atoms is None:
        atoms = complete_molecules(struct, cfg)
    if cfg.view_hide_ch:
        atoms = _hide_ch(atoms)
    if len(atoms) < 2:
        raise ValueError("render: fewer than 2 atoms after completion")
    R, view, findings = orient(struct, cfg, atoms)
    _resolve(findings, cfg, "orientation")
    alt = disorder_alternatives(struct)
    desat = _desat_mask(atoms, alt) if getattr(cfg, "color_by_component", False) else {}

    P = np.array([a["xyz"] for a in atoms])
    Pc = P - P.mean(0)
    Pr = (R @ Pc.T).T if R is not None else Pc

    apply_style(cfg)
    fig, ax = figure(cfg, subplot_kw={"projection": "3d"})
    bonds = _bond_pairs(atoms, alt)
    for i, j in bonds:
        ci = _desat(_rgb(_bond_color(atoms[i]["sym"]))) if desat.get(id(atoms[i])) else _bond_color(atoms[i]["sym"])
        cj = _desat(_rgb(_bond_color(atoms[j]["sym"]))) if desat.get(id(atoms[j])) else _bond_color(atoms[j]["sym"])
        full = np.array([Pr[i], Pr[j]])
        # ONE continuous underlay -> no midpoint seam (the seam was showing the white page
        # through same-colour C-C bonds, splitting each aromatic bond in two).
        ax.plot(full[:, 0], full[:, 1], full[:, 2], color=cj, lw=2.0,
                solid_capstyle="round", zorder=1)
        if ci != cj:                                          # two-tone: overlay the i-half
            mid = (Pr[i] + Pr[j]) / 2
            seg = np.array([Pr[i], mid])
            ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color=ci, lw=2.0,
                    solid_capstyle="round", zorder=1)
    for hi, ai in _hbond_pairs(atoms, cfg):                  # dashed H...A (from the hydrogen)
        seg = np.array([Pr[hi], Pr[ai]])
        ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color="0.25", lw=0.9, ls=(0, (4, 3)), zorder=2)
    for a, p in zip(atoms, Pr):
        col = _desat(_rgb(_color(a["sym"]))) if desat.get(id(a)) else _color(a["sym"])
        ax.scatter(p[0], p[1], p[2], s=_size(a["sym"]), color=col,
                   edgecolors="k", linewidths=0.3, depthshade=True, zorder=3)
    if cfg.view_label_atoms != "none":
        hbset = _hbond_label_set(atoms, cfg) if cfg.view_label_atoms == "hbond" else None
        labelled = set()
        for i, (a, p) in enumerate(zip(atoms, Pr)):
            sup = a.get("sup")
            if a["occ"] < 0.99 or (a["label"], sup) in labelled:   # skip disordered + dedupe (label,sym-image)
                continue
            if (cfg.view_label_atoms == "all"
                    or (cfg.view_label_atoms == "hbond" and i in hbset)
                    or (cfg.view_label_atoms not in ("all", "hbond") and a["sym"] not in ("C", "H"))):
                from matplotlib import patheffects as _pe
                txt = f"$\\mathrm{{{a['label']}}}^{{\\mathrm{{{sup}}}}}$" if sup else a["label"]
                ax.text(p[0], p[1], p[2], txt, color="black", zorder=4,
                        fontsize=6.0 * getattr(cfg, "view_label_size", 1.0),
                        path_effects=[_pe.withStroke(linewidth=2.5, foreground="white")])
                labelled.add((a["label"], sup))

    if cell_box:
        c0 = P.mean(0)
        for e0, e1 in _cell_edges(struct):
            seg = np.array([e0, e1]) - c0
            seg = (R @ seg.T).T if R is not None else seg
            ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color="0.15", lw=1.0, zorder=2)

    # tight, undistorted framing: equal cubic limits (no axis stretched) + fill the frame
    ax.set_axis_off()
    mid = (Pr.max(0) + Pr.min(0)) / 2.0
    half = (float((Pr.max(0) - Pr.min(0)).max()) / 2.0 * 1.08) or 1.0
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)
    try:
        ax.set_box_aspect((1, 1, 1), zoom=1.35)       # zoom: matplotlib >= 3.8
    except TypeError:
        ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=view[0], azim=view[1], roll=view[2])
    return fig, ax


def _rgb(hexstr):
    h = hexstr.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _choose_backend(cfg):
    want = getattr(cfg, "view_renderer", "auto")
    if want == "matplotlib":
        return "matplotlib"
    try:
        import pyvista  # noqa: F401
        return "pyvista"
    except ImportError:
        if want == "pyvista":
            raise ImportError("Required package not found: pyvista. Install with:\n"
                              "    python -m pip install pyvista\n"
                              "(or set cfg.view_renderer='matplotlib' for the no-dep fallback).")
        return "matplotlib"


def _atom_op(struct, label, xyz_cart):
    """Recover the symmetry op that generated a rendered atom at xyz_cart (cartesian), so its
    ADP can be rotated to match (U_image = R U R^T). None if not found (treat as identity)."""
    asym = next((a for a in struct.atoms if a["label"] == label), None)
    if asym is None:
        return None
    g = _gemmi()
    fr = struct.cell.fractionalize(g.Position(float(xyz_cart[0]), float(xyz_cart[1]), float(xyz_cart[2])))
    target = np.array([fr.x, fr.y, fr.z]) % 1.0
    for op in struct.ops:
        gp = np.array(op.apply_to_xyz(list(asym["frac"]))) % 1.0
        if np.all(np.abs((gp - target + 0.5) % 1.0 - 0.5) < 5e-3):
            return op
    return None


def _overlay_labels(pl, img, sel, res, cfg):
    """Composite atom labels as HALO text (black glyphs + white outline) at each atom's
    projected pixel — legible over dark atoms WITHOUT a filled box that would hide the
    molecule (VTK's own labels can't draw a text halo). VTK world->display gives the pixel."""
    from PIL import Image, ImageDraw, ImageFont
    import matplotlib.font_manager as fm
    try:
        import vtk
        Coord = vtk.vtkCoordinate
    except Exception:
        from vtkmodules.vtkRenderingCore import vtkCoordinate as Coord
    coord = Coord(); coord.SetCoordinateSystemToWorld()
    ren = pl.renderer
    arr = img[..., :3].copy(); H = arr.shape[0]
    pim = Image.fromarray(arr); draw = ImageDraw.Draw(pim)
    fpx = max(12, int(res / 55 * getattr(cfg, "view_label_size", 1.0)))
    try:
        font = ImageFont.truetype(fm.findfont("DejaVu Sans"), fpx)
    except Exception:
        font = ImageFont.load_default()
    sw = max(2, fpx // 7)
    sub_fpx = max(9, int(fpx * 0.66))                    # symmetry-superscript size
    try:
        subfont = ImageFont.truetype(fm.findfont("DejaVu Sans"), sub_fpx)
    except Exception:
        subfont = font
    gap = 0.6 * fpx if getattr(cfg, "view_label_offset", False) else 0.0
    anchor = "lm" if gap else "mm"
    for ent in sel:
        lab, p = ent[0], ent[1]
        sup = ent[2] if len(ent) > 2 else None
        coord.SetValue(float(p[0]), float(p[1]), float(p[2]))
        dx, dy = coord.GetComputedDisplayValue(ren)
        x, y = dx + gap, H - dy
        draw.text((x, y), lab, font=font, fill=(0, 0, 0), anchor=anchor,
                  stroke_width=sw, stroke_fill=(255, 255, 255))
        if sup:                                          # raised superscript at the label's top-right
            bb = draw.textbbox((x, y), lab, font=font, anchor=anchor)
            draw.text((bb[2] + sub_fpx * 0.1, bb[1] - sub_fpx * 0.15), sup, font=subfont,
                      fill=(0, 0, 0), anchor="lm", stroke_width=max(1, sw // 2),
                      stroke_fill=(255, 255, 255))
    return np.asarray(pim)


def _desat(rgb, f=0.72):
    """Luminance-preserving desaturation (blend toward the colour's own grey) — mutes a
    cocrystal's secondary component without changing its brightness."""
    lum = 0.30 * rgb[0] + 0.59 * rgb[1] + 0.11 * rgb[2]
    return tuple(c * (1 - f) + lum * f for c in rgb)


def _desat_mask(atoms, alt):
    """{id(atom): True} for atoms NOT in the largest molecule — so color_by_component keeps the
    largest component (the 'main' molecule) in full element colour and mutes the rest (the
    coformer/solvent). Empty (no muting) when there's a single component."""
    comps = _components(atoms, alt)
    if len(comps) < 2:
        return {}
    mx = max(len(c) for c in comps)
    mask = {}
    for c in comps:
        muted = len(c) < mx
        for a in c:
            mask[id(a)] = muted
    return mask


def _render_pyvista(struct, cfg, atoms=None, cell_box=False):
    """HIGH-QUALITY offscreen 3D render (VTK): smooth-shaded spheres + cylinder bonds with a
    real depth buffer (correct occlusion -> NO vertex gaps), SSAA + SSAO, orthographic camera
    from the orientation engine. `atoms` overrides the default complete-molecule set; `cell_box`
    draws the unit-cell edges + a/b/c. Returns an RGB ndarray with halo labels composited."""
    import pyvista as pv
    if atoms is None:
        atoms = complete_molecules(struct, cfg)
    if cfg.view_hide_ch:
        atoms = _hide_ch(atoms)
    if len(atoms) < 2:
        raise ValueError("crystal_view: fewer than 2 atoms after completion")
    R, _, findings = orient(struct, cfg, atoms)
    P = np.array([a["xyz"] for a in atoms])
    ctr = P.mean(0)
    if R is None:                          # 'custom' angles target the matplotlib backend
        heavy = [a for a in atoms if a["sym"] != "H"]
        Q = np.array([a["xyz"] for a in (heavy or atoms)])
        Qc = Q - Q.mean(0)
        _, evec = np.linalg.eigh(Qc.T @ Qc)
        x = evec[:, 2]; z = evec[:, 0]; y = np.cross(z, x)
        R = np.array([x, y, z])
        if np.linalg.det(R) < 0:
            R = np.array([x, -y, z])
        findings.append(("INFO", "custom view_angles apply to the matplotlib backend; "
                                 "pyvista used a PCA camera (use view_orientation='vector' "
                                 "for an explicit pyvista view)"))
    _resolve(findings, cfg, "orientation")
    alt = disorder_alternatives(struct)
    desat = _desat_mask(atoms, alt) if getattr(cfg, "color_by_component", False) else {}

    res = int(getattr(cfg, "view_resolution", 1800))
    pl = pv.Plotter(off_screen=True, window_size=[res, res], lighting="light kit")
    pl.set_background("white")

    def crad(sym):
        try:
            r = _gemmi().Element(sym).covalent_r or 0.7
        except Exception:
            r = 0.7
        return 0.30 * (r / 0.7) + 0.08              # ball-and-stick: balls < vdW, H smaller

    sph = dict(smooth_shading=True, specular=0.3, specular_power=15)
    ellipsoid_mode = False
    if getattr(cfg, "view_style", "ball_stick") == "ellipsoid":
        if has_adp(struct):
            ellipsoid_mode = True
        else:
            _resolve([("WARN", "view_style='ellipsoid' but no anisotropic U in the CIF "
                                      "(_atom_site_aniso_U_*) — drawing ball-and-stick")], cfg, "adp")
    prob = getattr(cfg, "adp_probability", 0.5)
    for a, p in zip(atoms, P):
        col = _rgb(_color(a["sym"]))
        if desat.get(id(a)):
            col = _desat(col)
        if ellipsoid_mode and a["sym"] != "H" and a["label"] in struct.aniso:
            Uc = u_cart(struct, struct.aniso[a["label"]])
            op = _atom_op(struct, a["label"], p)
            if op is not None:                                    # rotate ADP for sym images
                Rc = op_rot_cart(struct, op)
                Uc = Rc @ Uc @ Rc.T
            w, Vv = np.linalg.eigh(Uc)
            if np.all(w > 1e-9):                                  # positive-definite -> ellipsoid
                semi = _adp_scale(prob) * np.sqrt(w)
                T = np.eye(4); T[:3, :3] = Vv @ np.diag(semi); T[:3, 3] = p
                m = pv.Sphere(radius=1.0, theta_resolution=36, phi_resolution=36)
                m.transform(T, inplace=True)
                pl.add_mesh(m, color=col, **sph)
                continue                                          # else fall through to a sphere
        pl.add_mesh(pv.Sphere(radius=crad(a["sym"]), center=p, theta_resolution=48,
                              phi_resolution=48), color=col, **sph)
    for i, j in _bond_pairs(atoms, alt):
        a_, b_ = P[i], P[j]; d = b_ - a_; L = float(np.linalg.norm(d)); mid = (a_ + b_) / 2
        ri, rj = _rgb(_bond_color(atoms[i]["sym"])), _rgb(_bond_color(atoms[j]["sym"]))
        if desat.get(id(atoms[i])):
            ri = _desat(ri)
        if desat.get(id(atoms[j])):
            rj = _desat(rj)
        if ri == rj:
            pl.add_mesh(pv.Cylinder(center=mid, direction=d, radius=0.11, height=L, resolution=28),
                        color=ri, **sph)
        else:
            for c, q in ((ri, (a_ + mid) / 2), (rj, (b_ + mid) / 2)):
                pl.add_mesh(pv.Cylinder(center=q, direction=d, radius=0.11, height=L / 2, resolution=28),
                            color=c, **sph)
    for hi, ai in _hbond_pairs(atoms, cfg):          # dashed H...A (originates at the hydrogen)
        a_, b_ = P[hi], P[ai]; d = b_ - a_; L = float(np.linalg.norm(d))
        n = max(3, int(L / 0.35))
        for t in range(0, n, 2):
            c0 = a_ + d * (t / n); c1 = a_ + d * (min(t + 1, n) / n)
            pl.add_mesh(pv.Cylinder(center=(c0 + c1) / 2, direction=d, radius=0.045,
                                    height=float(np.linalg.norm(c1 - c0)), resolution=12),
                        color=(0.35, 0.35, 0.35))
    # collect heteroatom labels (skip disordered + dedupe); composited as HALO text after
    # the render (below) so they sit at the atom without a box hiding the molecule.
    sel, seen = [], set()
    if cfg.view_label_atoms != "none":
        hbset = _hbond_label_set(atoms, cfg) if cfg.view_label_atoms == "hbond" else None
        for i, (a, p) in enumerate(zip(atoms, P)):
            sup = a.get("sup")
            dk = (a["label"], sup)
            if a["occ"] < 0.99 or dk in seen:
                continue
            if (cfg.view_label_atoms == "all"
                    or (cfg.view_label_atoms == "hbond" and i in hbset)
                    or (cfg.view_label_atoms not in ("all", "hbond") and a["sym"] not in ("C", "H"))):
                sel.append((a["label"], p, sup)); seen.add(dk)

    if cell_box:
        for p0, p1 in _cell_edges(struct):
            pl.add_mesh(pv.Line(p0, p1), color=(0.15, 0.15, 0.15), line_width=3)
        sel += [("a", struct.cart([1.08, 0, 0]), None), ("b", struct.cart([0, 1.08, 0]), None),
                ("c", struct.cart([0, 0, 1.08]), None)]

    D = float(np.ptp(P, axis=0).max()) * 3 + 5
    pl.enable_parallel_projection()
    pl.camera_position = [tuple(ctr + R[2] * D), tuple(ctr), tuple(R[1])]
    pl.reset_camera()
    try:
        pl.enable_anti_aliasing("ssaa")
    except Exception:
        pass
    if getattr(cfg, "view_ssao", True):
        try:
            pl.enable_ssao()
        except Exception:
            pass

    img = pl.screenshot(return_img=True)
    if sel:
        img = _overlay_labels(pl, img, sel, res, cfg)
    pl.close()
    return img


def render(struct, cfg, atoms=None, cell_box=False):
    """Render a structure view. Returns (backend, obj): ('pyvista', RGB ndarray) or
    ('matplotlib', (fig, ax)). `atoms` overrides the default complete-molecule set (the
    unit-cell / packing / H-bond-environment views pass their own); `cell_box` draws the
    unit-cell edges + a/b/c. PyVista default, matplotlib fallback. Write with save()."""
    backend = _choose_backend(cfg)
    if backend == "pyvista":
        return "pyvista", _render_pyvista(struct, cfg, atoms, cell_box)
    return "matplotlib", _render_matplotlib(struct, cfg, atoms, cell_box)


def save(rendered, basename, cfg):
    """Write the structure image. pyvista -> a high-res PNG (raster is the convention for
    structure images); matplotlib -> save_fig (vector + raster). Returns the paths
    written (open the PNG with the Read tool for the perceptual QA pass)."""
    backend, obj = rendered
    if backend == "pyvista":
        import os
        from PIL import Image
        png = basename + ".png"
        os.makedirs(os.path.dirname(png) or ".", exist_ok=True)
        Image.fromarray(obj).save(png)            # obj is the composited RGB ndarray
        return [png]
    fig, _ = obj
    return save_fig(fig, basename, cfg)


# ------------------------------------------------------------------- packing (Phase 2)
def _components(atoms, alt=None):
    """Connected components (whole molecules) of `atoms` by covalent bonds, disorder-aware."""
    from scipy.spatial import cKDTree
    n = len(atoms)
    if n == 0:
        return []
    xyz = np.array([a["xyz"] for a in atoms])
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    for i, j in cKDTree(xyz).query_pairs(r=2.6):
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        if d > _covsum(atoms[i]["sym"], atoms[j]["sym"]) + BOND_TOL:
            continue
        if atoms[i]["occ"] < 1 and atoms[j]["occ"] < 1 and d < DISORDER_MIN:
            continue
        if alt and atoms[j]["label"] in alt.get(atoms[i]["label"], set()):
            continue                                          # disorder alternative -> not bonded
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    comps = {}
    for idx in range(n):
        comps.setdefault(find(idx), []).append(atoms[idx])
    return list(comps.values())


def _pack_atoms(struct, cfg, cells=None):
    """WHOLE molecules whose CENTROID lies inside the (Nx,Ny,Nz) block — so a single cell shows
    exactly its Z formula units (each molecule ONCE), not every boundary fragment grown into a
    duplicate. Built over a -1..N+1 supercell, grouped into molecules, then centroid-filtered."""
    nx, ny, nz = cells or getattr(cfg, "pack_cells", (1, 1, 1))
    base = expand(struct)
    if getattr(cfg, "cell_fill", "molecule") == "clip":   # literal cell contents, molecules cut at the outer box
        out = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    t = np.array([i, j, k], float)
                    for u in base:
                        fr = u["frac"] + t
                        out.append({"sym": u["sym"], "label": u["label"], "occ": u["occ"],
                                    "frac": fr, "xyz": struct.cart(fr)})
        return out
    sup = []
    for i in range(-1, nx + 1):
        for j in range(-1, ny + 1):
            for k in range(-1, nz + 1):
                t = np.array([i, j, k], float)
                for u in base:
                    fr = u["frac"] + t
                    sup.append({"sym": u["sym"], "label": u["label"], "occ": u["occ"],
                                "frac": fr, "xyz": struct.cart(fr)})
    out = []
    for comp in _components(sup, disorder_alternatives(struct)):
        cen = np.mean([a["frac"] for a in comp], axis=0)
        if (-1e-4 <= cen[0] < nx) and (-1e-4 <= cen[1] < ny) and (-1e-4 <= cen[2] < nz):
            out.extend(comp)
    return out


def _cell_edges(struct):
    c = {(i, j, k): struct.cart([i, j, k]) for i in (0, 1) for j in (0, 1) for k in (0, 1)}
    edges = []
    for corner in c:
        for ax in range(3):
            if corner[ax] == 0:
                nb = list(corner); nb[ax] = 1
                edges.append((c[corner], c[tuple(nb)]))
    return edges


def _roman(n):
    out, vals = "", [(10, "x"), (9, "ix"), (5, "v"), (4, "iv"), (1, "i")]
    for v, s in vals:
        while n >= v:
            out += s; n -= v
    return out


def _sym_op_text(struct, code):
    """Human-readable symmetry op for an 'n_pqr' code (e.g. '-x+1, -y+1, -z+1') for the caption key."""
    g = _gemmi()
    try:
        op = struct.ops[int(code.split("_")[0]) - 1]
        new = g.Op(op.triplet())
        if "_" in code and len(code.split("_")[1]) == 3:
            t = code.split("_")[1]
            tr = list(new.tran)                          # gemmi returns a copy; assign the whole list back
            for i in range(3):
                tr[i] += (int(t[i]) - 5) * new.DEN
            new.tran = tr
        return new.triplet().replace(",", ", ")
    except Exception:
        return code


def _hbond_env_atoms(struct, cfg, central=None):
    """The asymmetric-unit molecule(s) + the symmetry neighbours they hydrogen-bond to, at the
    extent set by cfg.hbond_neighbour: 'whole' (full neighbour molecule), 'stub' (contact atom +
    one bonded shell), or 'site' (contact atom only). Symmetry-generated neighbour atoms are
    annotated with their 'n_pqr' code + a roman superscript, and the caption key is logged.
    `central` defaults to the whole asymmetric unit; pass a subset (e.g. one molecule) to CROP the
    figure to just that molecule's H-bond environment (useful for Z'>1 structures)."""
    if central is None:
        central = complete_molecules(struct, cfg)
    for a in central:
        a["neighbour"] = False
    alt = disorder_alternatives(struct)
    donors = set(cfg.hbond_donors) | ({"C"} if cfg.hbond_weak else set())
    acc = set(cfg.hbond_acceptors)
    floor = cfg.hbond_angle_min
    cell_atoms = [{**u, "xyz": struct.cart(u["frac"])} for u in expand(struct)]
    origin = struct.cart(np.zeros(3))
    sup = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                shift = struct.cart(np.array([di, dj, dk], float)) - origin
                for a in cell_atoms:
                    sup.append({"sym": a["sym"], "label": a["label"], "occ": a["occ"],
                                "frac": a["frac"] + np.array([di, dj, dk], float), "xyz": a["xyz"] + shift})

    def key(p):
        return (int(round(p[0] * 50)), int(round(p[1] * 50)), int(round(p[2] * 50)))

    central_keys = {key(a["xyz"]) for a in central}

    def hbond_ok(D, h, A):
        dDH = float(np.linalg.norm(D["xyz"] - h["xyz"]))
        if dDH > _covsum(D["sym"], "H") + BOND_TOL or D["sym"] not in donors:
            return False
        hpos = h["xyz"]
        if cfg.xh_normalize and D["sym"] in NEUTRON_XH and dDH > 0:
            hpos = D["xyz"] + (h["xyz"] - D["xyz"]) / dDH * NEUTRON_XH[D["sym"]]
        DA = float(np.linalg.norm(A["xyz"] - D["xyz"]))
        if not (0.4 < DA <= 4.0):
            return False
        if float(np.linalg.norm(A["xyz"] - hpos)) > _vdw("H") + _vdw(A["sym"]):
            return False
        return _angle(D["xyz"], hpos, A["xyz"]) >= floor

    cheavy = [a for a in central if a["sym"] != "H"]
    seeds = []
    for h in [a for a in central if a["sym"] == "H"]:               # central donor -> outside acceptor
        if not cheavy:
            break
        D = min(cheavy, key=lambda a: np.linalg.norm(a["xyz"] - h["xyz"]))
        for A in sup:
            if A["sym"] in acc and key(A["xyz"]) not in central_keys and hbond_ok(D, h, A):
                seeds.append(A)
    supheavy = [a for a in sup if a["sym"] != "H"]
    cacc = [a for a in central if a["sym"] in acc]
    for h in [a for a in sup if a["sym"] == "H" and key(a["xyz"]) not in central_keys]:  # outside donor -> central acceptor
        if not supheavy:
            break
        D = min(supheavy, key=lambda a: np.linalg.norm(a["xyz"] - h["xyz"]))
        if any(hbond_ok(D, h, A) for A in cacc):
            seeds.append(D)

    # ---- assemble neighbour atoms at the requested extent
    mode = getattr(cfg, "hbond_neighbour", "stub")
    seen_seed, seed_list = set(central_keys), []
    for s in seeds:                                   # unique contact atoms (the seeds)
        kk = key(s["xyz"])
        if kk not in seen_seed:
            seen_seed.add(kk); seed_list.append(s)
    have, neigh = set(central_keys), []
    for s in seed_list:
        have.add(key(s["xyz"])); neigh.append(s)
    if mode in ("stub", "whole"):                     # grow outward (1 shell for stub, fully for whole)
        depth = {id(s): 0 for s in seed_list}
        frontier = list(seed_list)
        while frontier:
            nf = []
            for a in frontier:
                if mode == "stub" and depth[id(a)] >= 1:
                    continue
                for s in sup:
                    kk = key(s["xyz"])
                    if kk in have or not (0.4 < float(np.linalg.norm(a["xyz"] - s["xyz"]))
                                          <= _covsum(a["sym"], s["sym"]) + BOND_TOL):
                        continue
                    if a["occ"] < 1 and s["occ"] < 1 and \
                            float(np.linalg.norm(a["xyz"] - s["xyz"])) < DISORDER_MIN:
                        continue
                    if s["label"] in alt.get(a["label"], set()):   # disorder-alternative -> not bonded
                        continue
                    have.add(kk); neigh.append(s); nf.append(s); depth[id(s)] = depth[id(a)] + 1
            frontier = nf

    # ---- annotate neighbours: symmetry code + roman superscript; log the caption key
    codes = []
    for a in neigh:
        a["neighbour"] = True
        a["symcode"] = _sym_code(struct, a["label"], a["frac"])
        if a["symcode"] and a["symcode"] != "." and a["symcode"] not in codes:
            codes.append(a["symcode"])
    roman = {c: _roman(i + 1) for i, c in enumerate(codes)}
    for a in neigh:
        a["sup"] = roman.get(a.get("symcode"))
    if codes:
        keytxt = "; ".join(f"({roman[c]}) {_sym_op_text(struct, c)}" for c in codes)
        _resolve([("INFO", f"symmetry key for the caption: {keytxt}")], cfg, "hbond_env")
    return central + neigh


def render_unit_cell(struct, cfg):
    """Contents of a SINGLE unit cell (whole molecules) + the unit-cell box and a/b/c axes."""
    return render(struct, cfg, atoms=_pack_atoms(struct, cfg, cells=(1, 1, 1)), cell_box=True)


def render_packing(struct, cfg):
    """Packing diagram: cfg.pack_cells unit cells of whole molecules with H-bonds dashed, plus
    the unit-cell box and a/b/c axes. Viewing down a cell axis (cfg.view_orientation='axis_a'|
    'axis_b'|'axis_c') usually reads best."""
    return render(struct, cfg, atoms=_pack_atoms(struct, cfg), cell_box=True)


def render_hbond_environment(struct, cfg):
    """Asymmetric-unit molecule(s) + the symmetry neighbours they hydrogen-bond to, with the
    intermolecular H-bonds dashed — the H-bonding environment."""
    return render(struct, cfg, atoms=_hbond_env_atoms(struct, cfg), cell_box=False)


# ======================================================================
# from pxrd_realism.py
# ======================================================================
"""
py — make a calculated PXRD pattern look like a real lab Cu-Kα scan.

Three physical effects a bare stick pattern / single-wavelength profile misses, added as
NUMPY-ONLY functions that operate on a REFLECTION LIST — a list of
`(two_theta_deg, intensity, (h, k, l))` — plus the cell metric. Keeping the physics off the
Dans_Diffraction path makes it testable without heavy deps (Dans supplies the structure-factor
intensities that seed the list; see `crystal.md` for wiring):

  march_dollase   preferred-orientation INTENSITY correction (platy / needle habit)
  kalpha2_doublet Cu Kα1/Kα2 peak SPLITTING (α2 at ~half intensity; splitting grows with angle)
  caglioti_fwhm + pseudo_voigt   angle-dependent peak WIDTH and shape
  simulate_pattern  composes them onto a 2θ grid (defaults pulled from cfg)

Nothing here computes structure factors or systematic absences — feed it real reflection
intensities. Preferred orientation is the per-reflection-orientation form; a fully rigorous PO
averages over symmetry-equivalents (pass those in the list if you have them). Validated
closed-form in skill_validation/pxrd/.
"""

# Cu radiation (Å) and the Kα2/Kα1 integrated-intensity ratio (~0.5).
CU_KA1 = 1.540598
CU_KA2 = 1.544426
CU_KA2_RATIO = 0.5


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
    x = np.asarray(x, float)
    hwhm = fwhm / 2.0
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    gauss = np.exp(-0.5 * ((x - center) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    lorentz = (hwhm / np.pi) / ((x - center) ** 2 + hwhm ** 2)
    return eta * lorentz + (1 - eta) * gauss


def _pick(cfg, name, default, override):
    if override is not None:
        return override
    return getattr(cfg, name, default) if cfg is not None else default


def simulate_pattern(reflections, x_grid, cfg=None, *, U=None, V=None, W=None, eta=None,
                     kalpha2=None, lam1=CU_KA1, lam2=None, ratio=None,
                     po_hkl=None, march_r=None, Gs=None, normalize=True):
    """Build a realistic profile on `x_grid` (2θ deg) from a reflection list. Pipeline:
    March-Dollase (if po_hkl and r≠1 and Gs given) → Kα2 doublet (if kalpha2) → sum of
    area-normalized pseudo-Voigts with Caglioti(θ) widths. Unspecified parameters default from
    `cfg` (pxrd_caglioti → U,V,W; pxrd_lorentz_fraction → eta; pxrd_kalpha2 / pxrd_wavelength2 /
    pxrd_kalpha2_ratio; pxrd_po_axis → po_hkl; pxrd_march_r). Returns the intensity on `x_grid`
    (scaled to 100 at the max when `normalize`)."""
    uvw = _pick(cfg, "pxrd_caglioti", (0.01, -0.005, 0.008), None)
    U = uvw[0] if U is None else U
    V = uvw[1] if V is None else V
    W = uvw[2] if W is None else W
    eta = _pick(cfg, "pxrd_lorentz_fraction", 0.5, eta)
    kalpha2 = _pick(cfg, "pxrd_kalpha2", False, kalpha2)
    lam2 = _pick(cfg, "pxrd_wavelength2", None, lam2) or CU_KA2
    ratio = _pick(cfg, "pxrd_kalpha2_ratio", CU_KA2_RATIO, ratio)
    po_hkl = _pick(cfg, "pxrd_po_axis", None, po_hkl)
    march_r = _pick(cfg, "pxrd_march_r", 1.0, march_r)

    refl = list(reflections)
    if po_hkl is not None and march_r != 1.0 and Gs is not None:
        refl = apply_march_dollase(refl, po_hkl, Gs, march_r)
    if kalpha2:
        refl = kalpha2_doublet(refl, lam1, lam2, ratio)

    x = np.asarray(x_grid, float)
    y = np.zeros_like(x)
    for tth, I, _hkl in refl:
        y = y + I * pseudo_voigt(x, tth, float(caglioti_fwhm(tth, U, V, W)), eta)
    if normalize and y.max() > 0:
        y = y / y.max() * 100.0
    return y


# ======================================================================
# from cocrystal.py
# ======================================================================
"""
py — quantitative helpers for cocrystal IDENTIFICATION.

The DECISION layer is `references/cocrystal_id.md` (which evidence settles cocrystal vs salt
vs physical mixture vs polymorph); THIS module is the arithmetic that doc routes to. Built in
waves:

  wave 1:
    delta_pka / classify_ionisation      — the ΔpKa salt / cocrystal / continuum call [zero-dep]
    rwp / sum_of_parents / phase_report  — NNLS "new phase vs physical mixture"        [scipy]
  next:
    class_fom                             — one-class / classification figures of merit [numpy]

Nothing here decides the science. ΔpKa *ranks salt-risk*; it does not prove protonation — the
0–3 continuum is explicitly unpredictable, confirm by FTIR (the ionised-group band) or ssNMR
(see `cocrystal_id.md`). Every function returns numbers + a provenance caption, never a bare
label; `classify_ionisation` NEVER raises.

Validation: `skill_validation/cocrystal/` (closed-form + reference parity), bridged into
`python -m pytest tests/`.
"""


# ============================================================= ΔpKa (salt vs cocrystal)
# Cruz-Cabeza 2012 (CrystEngComm 14:6362), the "pKa rule". ΔpKa is defined with the BASE's
# conjugate-acid pKa minus the ACID's pKa — both the pKa of the PROTONATED species:
#     ΔpKa = pKa(BH+) − pKa(AH)
# Simplified 3-zone rule (thresholds configurable): ΔpKa > 3 → salt, < 0 → cocrystal, else the
# salt–cocrystal continuum (proton transfer unpredictable — must be confirmed experimentally).
SALT_THRESHOLD = 3.0
COCRYSTAL_THRESHOLD = 0.0


def delta_pka(pka_acid, pka_base_conjugate):
    """ΔpKa for the pKa rule = pKa(base's conjugate acid, BH+) − pKa(acid, AH). Both inputs are
    the pKa of the PROTONATED species; returns the float. (The pKa values are look-ups, not
    computed here — fetch + cite them per the web-search discipline in SKILL.md.)"""
    return float(pka_base_conjugate) - float(pka_acid)


def classify_ionisation(pka_acid, pka_base_conjugate,
                        salt_threshold=SALT_THRESHOLD, cocrystal_threshold=COCRYSTAL_THRESHOLD):
    """Predict salt vs cocrystal from the ΔpKa rule [Cruz-Cabeza 2012]. Returns a dict:
    `delta_pka`, `zone` ('salt' | 'cocrystal' | 'continuum'), `confident` (bool), the
    thresholds, and a `caption`. This is TRIAGE, not proof: the continuum
    (cocrystal_threshold ≤ ΔpKa ≤ salt_threshold, default 0–3) is unpredictable — confirm with
    FTIR (ionised-group band) or ssNMR. Never raises."""
    d = delta_pka(pka_acid, pka_base_conjugate)
    if d > salt_threshold:
        zone = "salt"
    elif d < cocrystal_threshold:
        zone = "cocrystal"
    else:
        zone = "continuum"
    confident = zone != "continuum"
    tail = "" if confident else " — UNPREDICTABLE, confirm experimentally"
    caption = (f"ΔpKa = {d:.2f} [= pKa(BH+) {float(pka_base_conjugate):.2f} − "
               f"pKa(AH) {float(pka_acid):.2f}] → {zone}{tail}; rule: >{salt_threshold:g} salt, "
               f"<{cocrystal_threshold:g} cocrystal, else continuum [Cruz-Cabeza 2012 pKa rule]")
    return {"delta_pka": d, "zone": zone, "confident": confident,
            "salt_threshold": float(salt_threshold),
            "cocrystal_threshold": float(cocrystal_threshold), "caption": caption}


# ================================== NNLS sum-of-parents (new phase vs physical mixture)
# The #1 cocrystal-ID question: did a NEW phase form, or is the product just a physical
# mixture of the two starting materials? A physical mixture is a NON-NEGATIVE linear
# superposition of the parent patterns (you can't have a negative amount of a phase); a new
# phase (cocrystal / salt) has peaks the parents cannot build. So fit the product as
# sum-of-parents by NNLS and read the residual: low Rwp + no unexplained peaks ⇒ mixture;
# high Rwp + positive residual peaks ⇒ a new phase is plausible. Works on PXRD or FTIR
# profiles on a common grid. This is a DISCRIMINATION diagnostic, not quantitative phase
# analysis (no reference-intensity-ratio / Rietveld) — the verdict is heuristic; confirm by
# indexing / FTIR (cocrystal_id.md) / SCXRD.


def rwp(y_obs, y_calc, weight="poisson", eps=1e-12):
    """Weighted profile residual (Rietveld goodness-of-fit):
        Rwp = sqrt( Σ w_i (y_obs_i − y_calc_i)² / Σ w_i y_obs_i² )
    Returns the FRACTION (0 = perfect; ×100 for %). weight:
      'poisson' → w_i = 1/max(y_obs_i, eps) — counting statistics, the crystallographic default
                  (assumes counts-like data with a non-zero background; on a normalised pattern
                  with true zeros prefer 'unit').
      'unit'    → w_i = 1 — an Rp-style unweighted profile factor, apt for normalised patterns.
    The NNLS fit itself (sum_of_parents) is UNWEIGHTED least squares; Rwp is a reported
    diagnostic on top of it."""
    yo = np.asarray(y_obs, float).ravel()
    yc = np.asarray(y_calc, float).ravel()
    if weight == "unit":
        w = np.ones_like(yo)
    elif weight == "poisson":
        w = 1.0 / np.maximum(yo, eps)
    else:
        raise ValueError("weight must be 'poisson' or 'unit'")
    den = float(np.sum(w * yo ** 2))
    if den <= 0:
        return float("nan")
    return float(np.sqrt(np.sum(w * (yo - yc) ** 2) / den))


def sum_of_parents(y_obs, parents, weight="poisson", eps=1e-12):
    """Fit an observed pattern as a NON-NEGATIVE linear combination of the parent patterns
    (NNLS) — the physical-mixture model. A true mixture reconstructs well (low Rwp, no
    systematic unexplained peaks); a genuine new phase does NOT (the parents can't build its
    new peaks). Patterns must share a common x-grid.

    y_obs   : (n,) observed intensities.
    parents : (k, n) array or list of k (n,)-patterns on the SAME grid (a (n, k) array is
              accepted too).
    Returns dict: `coefficients` (k,), `fractions` (coefficients / Σ — RELATIVE scale, NOT
    quantitative phase % without an RIR), `y_calc` (n,), `residual` (n,), `rwp`, `resid_norm`,
    `weight`, `caption`. NNLS via scipy.optimize.nnls (lazy import)."""
    from scipy.optimize import nnls
    y = np.asarray(y_obs, float).ravel()
    P = np.asarray(parents, float)
    if P.ndim == 1:
        P = P[None, :]
    if P.shape[1] != y.size and P.shape[0] == y.size:
        P = P.T                                  # accept (n, k) too
    if P.shape[1] != y.size:
        raise ValueError(f"parents shape {np.asarray(parents).shape} incompatible with y (n={y.size})")
    A = P.T                                       # (n, k): columns are parents
    coef, rnorm = nnls(A, y)
    y_calc = A @ coef
    resid = y - y_calc
    total = float(coef.sum())
    frac = coef / total if total > 0 else np.full_like(coef, np.nan)
    rw = rwp(y, y_calc, weight=weight, eps=eps)
    caption = (f"NNLS sum-of-parents: {coef.size} parents, relative scale "
               f"{np.array2string(frac, precision=3)} (NOT quantitative phase % — no RIR); "
               f"Rwp={rw*100:.1f}% [{weight} weights]. High Rwp / positive residual peaks ⇒ the "
               f"parents can't reconstruct the product ⇒ a new phase is plausible (confirm).")
    return {"coefficients": coef, "fractions": frac, "y_calc": y_calc, "residual": resid,
            "rwp": rw, "resid_norm": float(rnorm), "weight": weight, "caption": caption}


def unexplained_peaks(x, residual, reference=None, kind="new", rel_height=0.05, distance=None):
    """Peaks in the sum-of-parents residual. kind='new' → POSITIVE residual (intensity the
    parents can't explain — new-phase evidence); kind='lost' → NEGATIVE residual (parent
    intensity the product lacks). The height threshold is rel_height × max(|reference|) (or
    max(|residual|) if reference is None). Returns [(x, height), …] sorted by descending
    height. scipy.signal.find_peaks (lazy import)."""
    from scipy.signal import find_peaks
    x = np.asarray(x, float).ravel()
    r = np.asarray(residual, float).ravel()
    signal = r if kind == "new" else -r
    ref = float(np.max(np.abs(reference))) if reference is not None else float(np.max(np.abs(r)))
    thr = rel_height * ref if ref > 0 else 0.0
    idx, _ = find_peaks(signal, height=thr, distance=distance)
    peaks = [(float(x[i]), float(signal[i])) for i in idx]
    peaks.sort(key=lambda t: -t[1])
    return peaks


def phase_report(x, y_obs, parents, weight="poisson", peak_rel_height=0.05):
    """Turn-key new-phase-vs-physical-mixture report: NNLS sum-of-parents fit + Rwp + the
    new/lost unexplained-peak lists + a HEURISTIC verdict. The verdict is a detection (are
    there unexplained peaks above `peak_rel_height`?), not a magic Rwp cutoff — read it WITH
    the Rwp and confirm structurally (indexing / FTIR ionised-group band / SCXRD). Returns the
    `sum_of_parents` dict plus `new_peaks`, `lost_peaks`, `verdict`."""
    fit = sum_of_parents(y_obs, parents, weight=weight)
    new = unexplained_peaks(x, fit["residual"], y_obs, kind="new", rel_height=peak_rel_height)
    lost = unexplained_peaks(x, fit["residual"], y_obs, kind="lost", rel_height=peak_rel_height)
    verdict = ("unexplained intensity present — inconsistent with a pure physical mixture; a "
               "new phase (cocrystal/salt) is plausible" if new else
               "no unexplained peaks above threshold — consistent with a physical mixture")
    caption = (f"{fit['caption']} {len(new)} new / {len(lost)} lost peak(s) at "
               f"rel_height {peak_rel_height:.0%}. Verdict (heuristic — confirm by "
               f"indexing/FTIR/SCXRD): {verdict}.")
    return {**fit, "new_peaks": new, "lost_peaks": lost, "verdict": verdict, "caption": caption}


# =============================================== COCRYSTAL-FORMATION PREDICTOR (a-priori screen)
# Transparent, OPEN reimplementation of the two established pre-screens — to be VALIDATED, never
# trusted blindly (see skill_validation/cocrystal_predictor/ + cocrystal_id.md):
#   • molecular complementarity  — shape + polarity similarity  [Fabián 2009, CGD 9:1436]
#   • H-bond synthon competition — best-donor→best-acceptor, hetero vs homo  [Etter 1990, Acc.Chem.Res. 23:120]
# combined with delta_pKa (Cruz-Cabeza 2012, above). This does NOT reproduce the CSD-TRAINED H-bond
# propensity (Galek 2007) — those coefficients come from CSD statistics we don't have; cross-check
# our complementarity against Mercury's Molecular Complementarity tool instead (the reference impl).
# Every function returns numbers + a provenance caption and a *per-signal* breakdown — the value is a
# legible prediction you can audit, not a black-box score. RDKit is a lazy import (optional dep).

_RDKIT_SEED = 0xC0CC  # fixed → the conformer ensemble is deterministic/reproducible


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
        return Chem, AllChem, rdMolDescriptors
    except ImportError as e:
        raise ImportError("the cocrystal predictor needs RDKit: python -m pip install rdkit") from e


def shape_descriptors(coords):
    """Two size-normalised shape descriptors (S/L, M/L) from the principal-axis bounding box of a
    point cloud (N×3): eigen-decompose the coordinate covariance, take the box EXTENT along each
    principal axis, sort S≤M≤L, return (S/L, M/L). PURE. Translation/rotation-invariant by
    construction and scale-invariant (ratios) — both are asserted in the validation harness."""
    X = np.asarray(coords, float)
    C = X - X.mean(0)
    _, V = np.linalg.eigh(C.T @ C)
    ext = np.ptp(C @ V, axis=0)
    S, M, L = np.sort(ext)
    return (float(S / L), float(M / L)) if L > 1e-9 else (1.0, 1.0)


def gasteiger_dipole_debye(mol, coords):
    """|dipole| in Debye from RDKit Gasteiger partial charges and a conformer's coordinates
    (Å): μ = |Σ qᵢ rᵢ| × 4.803. Approximate (Gasteiger, not QM) — a RELATIVE polarity descriptor."""
    from rdkit.Chem import AllChem
    AllChem.ComputeGasteigerCharges(mol)
    q = np.nan_to_num([a.GetDoubleProp("_GasteigerCharge") for a in mol.GetAtoms()])
    return float(np.linalg.norm((q[:, None] * np.asarray(coords, float)).sum(0)) * 4.803)


def molecule_descriptors(smiles, n_conf=12, seed=_RDKIT_SEED):
    """3-D descriptors for a SMILES via an RDKit conformer ENSEMBLE — flexible molecules (diacid
    chains!) change shape per conformer, so we report the MEDIAN and the spread (sd), which the
    validation uses as a stability check. Returns dict: s_l, m_l (+ *_sd), dipole_debye (+ sd),
    hbd, hba, formula, n_conf. HBD/HBA are graph-based (RDKit CalcNumHBD/HBA)."""
    Chem, AllChem, rd = _rdkit()
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(base)
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, randomSeed=seed))
    shapes, dips = [], []
    for cid in ids:
        try:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid)
        except Exception:
            pass
        xyz = mol.GetConformer(cid).GetPositions()
        shapes.append(shape_descriptors(xyz))
        dips.append(gasteiger_dipole_debye(mol, xyz))
    shapes, dips = np.asarray(shapes), np.asarray(dips)
    return {"smiles": smiles, "formula": rd.CalcMolFormula(base),
            "s_l": float(np.median(shapes[:, 0])), "m_l": float(np.median(shapes[:, 1])),
            "s_l_sd": float(shapes[:, 0].std()), "m_l_sd": float(shapes[:, 1].std()),
            "dipole_debye": float(np.median(dips)), "dipole_sd": float(dips.std()),
            "hbd": int(rd.CalcNumHBD(base)), "hba": int(rd.CalcNumHBA(base)), "n_conf": len(ids)}


# Fabián-style thresholds: a pair is "complementary" (cocrystal-favoured) when its molecules are
# SIMILAR in shape and polarity. These cutoffs are APPROXIMATE and are CALIBRATED on the benchmark
# (validate.py) — not sacred constants; report the raw differences alongside the verdict.
MC_THRESHOLDS = {"s_l": 0.20, "m_l": 0.20, "dipole_debye": 4.0}


def molecular_complementarity(desc_a, desc_b, thresholds=None):
    """Fabián 2009 molecular-complementarity screen (approx): two molecules are cocrystal-favoured
    when their shape (S/L, M/L) and polarity (dipole) are SIMILAR. Returns per-descriptor absolute
    differences, which passed, an overall `complementary` bool, and a caption. TRIAGE — cross-check
    against Mercury's Molecular Complementarity tool (the CSD-calibrated reference)."""
    t = thresholds or MC_THRESHOLDS
    d = {"s_l": abs(desc_a["s_l"] - desc_b["s_l"]), "m_l": abs(desc_a["m_l"] - desc_b["m_l"]),
         "dipole_debye": abs(desc_a["dipole_debye"] - desc_b["dipole_debye"])}
    passes = {k: d[k] <= t[k] for k in d}
    ok = all(passes.values())
    caption = (f"molecular complementarity [Fabián 2009, open approx]: |Δ(S/L)|={d['s_l']:.2f}, "
               f"|Δ(M/L)|={d['m_l']:.2f}, |Δμ|={d['dipole_debye']:.1f} D → "
               f"{'complementary (cocrystal-favoured)' if ok else 'NOT complementary'} "
               f"(thresholds {t}; cross-check vs Mercury)")
    return {"complementary": ok, "differences": d, "passes": passes, "caption": caption}


# Functional groups (SMARTS) with ordinal donor/acceptor strengths (0 = n/a). Ranking follows the
# standard H-bond hierarchy behind Etter's rules (acid O–H > imide/amide N–H > O–H; carboxylate/
# sp2-N acceptors > carbonyl O > hydroxyl/ether O). Ordinal + heuristic — VALIDATE predicted
# synthons against known ones (validate.py), don't treat the numbers as energies.
_FG = [
    ("carboxylic acid O–H", "[CX3](=[OX1])[OX2H1]", 5, 4),
    ("imide N–H",           "[NX3H1]([CX3]=[OX1])[CX3]=[OX1]", 4, 3),
    ("amide N–H",           "[NX3;H1,H2][CX3]=[OX1]", 3, 3),
    ("aromatic N–H",        "[nH]", 4, 0),                 # imidazole/pyrrole (theophylline N7–H!)
    ("hydroxyl/enol O–H",   "[OX2H1]", 3, 2),
    ("pyridine/sp2 N",      "[nX2,$([NX2]=C)]", 0, 5),
    ("primary/sec amine",   "[NX3;H1,H2;!$(NC=O)]", 2, 4),
    ("carbonyl O",          "[CX3]=[OX1]", 0, 3),
    ("aromatic N (sub.)",   "[nX3]", 0, 2),
    ("sulfonyl O",          "[SX4](=[OX1])(=[OX1])", 0, 2),
]


def hbond_groups(smiles):
    """H-bond functional groups on a molecule → {groups, donors:[(name,strength)],
    acceptors:[(name,strength)]} via SMARTS. The input to the Etter synthon-competition call."""
    Chem, _, _ = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    groups, donors, acceptors = [], [], []
    for name, smarts, ds, as_ in _FG:
        n = len(mol.GetSubstructMatches(Chem.MolFromSmarts(smarts)))
        if n:
            groups.append((name, n))
            if ds:
                donors.append((name, ds))
            if as_:
                acceptors.append((name, as_))
    return {"groups": groups, "donors": donors, "acceptors": acceptors}


def hbond_competition(smiles_a, smiles_b):
    """Etter's rules as a synthon predictor: the best HETEROmeric donor→acceptor pair (A's donor to
    B's acceptor, or vice-versa) vs each partner's best HOMOmeric pair. If the heterosynthon score
    (donor+acceptor strengths) ≥ both homosynthons, a cocrystal is favoured. Returns the predicted
    `heterosynthon`, its score, the homo scores, a `hetero_favoured` bool, and a caption. Heuristic —
    validate the predicted synthon against the known one."""
    ga, gb = hbond_groups(smiles_a), hbond_groups(smiles_b)
    def pairs(donors, acceptors, tag_d, tag_a):
        return [(ds + as_, f"{tag_d}:{dn}···{tag_a}:{an}") for dn, ds in donors for an, as_ in acceptors]
    hetero = pairs(ga["donors"], gb["acceptors"], "A", "B") + pairs(gb["donors"], ga["acceptors"], "B", "A")
    homo_a = pairs(ga["donors"], ga["acceptors"], "A", "A")
    homo_b = pairs(gb["donors"], gb["acceptors"], "B", "B")
    best_h = max(hetero) if hetero else (0, "none")
    best_a = max(homo_a) if homo_a else (0, "none")
    best_b = max(homo_b) if homo_b else (0, "none")
    favoured = best_h[0] >= max(best_a[0], best_b[0]) and best_h[0] > 0
    caption = (f"H-bond synthon competition [Etter 1990]: best heterosynthon {best_h[1]} "
               f"(score {best_h[0]}) vs homo A {best_a[0]} / B {best_b[0]} → "
               f"{'heterosynthon competitive → cocrystal favoured' if favoured else 'homosynthons win → cocrystal disfavoured'}")
    return {"heterosynthon": best_h[1], "hetero_score": best_h[0],
            "homo_scores": (best_a[0], best_b[0]), "hetero_favoured": favoured, "caption": caption}


def predict_cocrystal(smiles_a, smiles_b, pka_acid=None, pka_base_conjugate=None,
                      thresholds=None, n_conf=12):
    """A-priori cocrystal SCREEN — three transparent, separately-reported signals (NOT a go/no-go
    oracle). The validation harness (skill_validation/cocrystal_predictor/) showed the combined
    binary does NOT beat a label-scramble control on a small diverse set, so this deliberately does
    NOT emit a single `likely` verdict. Use it for:
      • `synthon`      — the predicted heterosynthon [Etter 1990]: a genuinely useful FTIR-interpretation
                         aid ("expect THESE bands to shift"). Sensitive, but non-specific for yes/no.
      • `ionisation`   — the ΔpKa salt-vs-cocrystal filter [Cruz-Cabeza 2012] (if pKa given).
      • `complementarity` — shape/polarity [Fabián 2009]: report the raw Δ's to CROSS-CHECK against
                         Mercury's Molecular Complementarity; weak for rigid-API + flexible-coformer pairs.
    Go/no-go formation belongs to Mercury (CSD-calibrated) + the experiment. Pass pKa only for an acid/base pair."""
    da, db = molecule_descriptors(smiles_a, n_conf=n_conf), molecule_descriptors(smiles_b, n_conf=n_conf)
    comp = molecular_complementarity(da, db, thresholds)
    hb = hbond_competition(smiles_a, smiles_b)
    ion = classify_ionisation(pka_acid, pka_base_conjugate) if (pka_acid is not None and pka_base_conjugate is not None) else None
    caption = (f"cocrystal SCREEN(A,B): predicted synthon = {hb['heterosynthon']}"
               + (f"; ΔpKa → {ion['zone']}" if ion else "")
               + f"; complementarity |Δμ|={comp['differences']['dipole_debye']:.1f} D (cross-check Mercury). "
               f"NO binary verdict — validated as no-better-than-chance for formation; defer go/no-go to Mercury + experiment.")
    return {"descriptors": (da, db), "synthon": hb["heterosynthon"], "complementarity": comp,
            "hbond": hb, "ionisation": ion,
            "formation_confidence": "low — open screen did not beat chance in retrospective validation (skill_validation/cocrystal_predictor/validate.py, L3)",
            "caption": caption}


# ======================================================================
# analysis body: analysis.py
# ======================================================================
#!/usr/bin/env python3
"""
Prompt: "Doxorubicin peak areas from Bansal 2021, transmittance mode, 0.6-1.4 % w/w. Fit the
calibration, give me LOD and LOQ, and read off the concentration for peak areas of 3.5 and 6.0 mm2."

Data: Bansal, Singh & Kaur, BMC Chemistry 15:27 (2021), doi:10.1186/s13065-021-00752-3 — the
published per-point peak-area table (carbonyl band, baseline-corrected, transmittance mode).
"""
DESCRIPTION = ("Doxorubicin mid-IR calibration from the published carbonyl peak-area table "
               "(Bansal 2021, BMC Chem 15:27): fit, LOD/LOQ, and read-back of two unknowns")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))          # the skill root


# ====================================================================== CONFIG
CONC = [0.6, 0.8, 1.0, 1.2, 1.4]                 # % w/w      (Bansal 2021, Table: DOX transmittance)
AREA = [2.131, 2.850, 3.521, 3.913, 4.601]       # peak area / mm2
UNKNOWNS = [3.5, 6.0]                            # peak areas the user wants concentrations for
cfg = Config(
    output_dir=os.path.join(HERE, "out"),
    journal="nature", column="single", formats=("pdf", "png"),
    conf_level=0.95, lod_loq_method="residual_sd",
    refuse_extrapolation=True,                   # the default; shown here because the second unknown tests it
    description=DESCRIPTION,
)


def main():
    os.makedirs(cfg.output_dir, exist_ok=True)
    x, y = np.array(CONC), np.array(AREA)

    cal = fit(x, y, cfg)
    ll = lod_loq(cal, cfg)
    bias = intercept_test(cal, cfg)
    print(f"  y = {cal['slope']:.3f} x {cal['intercept']:+.3f}   R2 {cal['r2']:.4f}   n={cal['n']}")
    print(f"  LOD {ll['lod']:.3f} % w/w, LOQ {ll['loq']:.3f} % w/w  ({ll['method']})")
    print(f"  {bias['caption']}")

    # ---- read back the unknowns; the gate refuses anything outside the calibrated range
    readback = {}
    for area in UNKNOWNS:
        try:
            conc = predict(cal, area, cfg)
            readback[area] = conc
            print(f"  peak area {area:.2f} mm2 -> {conc:.3f} % w/w")
        except GateError as e:
            readback[area] = None
            print(f"  peak area {area:.2f} mm2 -> REFUSED: {e}")

    # ---- the figure: curve + bands + the mandatory residual panel
    fig, (ax, axr) = plot_calibration(x, y, cfg, cal)
    ax.set_ylabel("peak area / mm²")
    axr.set_ylabel("resid. / mm²")
    axr.set_xlabel("doxorubicin / % w/w")
    audit_layout(fig, cfg)
    save_fig(fig, os.path.join(cfg.output_dir, "bmc_dox_calibration"), cfg)

    ok = [f"{a:.2f} mm² → {c:.3f} % w/w" for a, c in readback.items() if c is not None]
    refused = [f"{a:.2f} mm²" for a, c in readback.items() if c is None]
    caption = (
        f"Figure. Doxorubicin calibration, carbonyl peak area (transmittance mode) against concentration, "
        f"from the published area table of Bansal et al. 2021. Line = least-squares fit "
        f"y = {cal['slope']:.3f}x {cal['intercept']:+.3f}, R² = {cal['r2']:.4f}, n = {cal['n']}; shaded = 95 % "
        f"confidence band, dashed = 95 % prediction band; residuals below. LOD = {ll['lod']:.3f} % w/w and "
        f"LOQ = {ll['loq']:.3f} % w/w (3.3σ/m and 10σ/m, σ = residual SD of the fit). Read-back: "
        f"{'; '.join(ok)}. Not reported: {', '.join(refused)} — outside the calibrated range "
        f"{cal['x_range'][0]:.1f}–{cal['x_range'][1]:.1f} % w/w, so no value is extrapolated.\n")
    open(os.path.join(HERE, "caption.txt"), "w", encoding="utf-8").write(caption)
    print("  wrote", os.path.join(cfg.output_dir, "bmc_dox_calibration.pdf"), "+ .png, caption.txt")


if __name__ == "__main__":
    main()
