"""
config.py  -  THE one tunable block for the whole skill.

The discipline this enforces: every knob a user might reasonably want to turn
lives HERE, in one dataclass, with units stated per field. Nothing
reasonable-to-tune should live buried in the body of a script.

Field groups:  I/O  |  style  |  spectra  |  calibration  |  verification
Units are stated per field. Anything not here is, by definition, not meant to
be casually tuned -- if you find yourself wanting to, add it here first.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, Sequence, List, Tuple

#   Single source of truth for the skill version. Stamped into a bundled
#   deliverable's provenance header (bundle.py) and into a figure file's metadata
#   (style.save_fig), so a figure can always be traced back to the code + version
#   that produced it. Bump on a user-visible behaviour change.
SKILL_VERSION = "1.0"


@dataclass
class Config:
    # ---------------------------------------------------------------- I/O
    input_dir: str = "."                 # where raw data lives
    output_dir: str = "./out"            # where figures + tables are written
    formats: Sequence[str] = ("pdf", "png")   # vector first; png is the raster proof
    #   one-line PLAIN-LANGUAGE description of what the figure shows (species, band,
    #   design, n). Stamped into the figure file's metadata by style.save_fig so the
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
    #   IDENTICALLY to every spectrum (this is asserted by verify.py).
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
    #   flanking local minima found by spectra.find_anchors; detect them ONCE on a
    #   reference spectrum via spectra.auto_windows and lock them onto
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
    outlier_mad_n: float = 3.5           # MAD multiplier for verify.flag_outliers_mad
    tick_overlap_tol_px: float = 2.0     # audit_layout flags tick labels closer than this
    clip_tol_px: float = 2.0             # audit_layout flags non-tick text past the canvas edge by this
    escape_tol_frac: float = 0.01        # audit_layout flags data running past the y-limits by more than this fraction of the range
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
                                         #   pseudo-Voigt eta by pxrd_realism.simulate_pattern
    pxrd_crosscheck: str = "auto"        # "auto" (pymatgen if present else bragg) | "bragg" | "none"
    # -- PXRD realism (pxrd_realism.py): make a calc pattern match a real lab Cu-Ka scan.
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


def default_cfg(cfg):
    """Resolve an optional cfg: None -> ONE Config() of the skill defaults, so every threshold a
    `cfg=None` entry point uses is a declared field (a misspelt one fails loudly) instead of a
    per-call literal that can drift from Config. The one mechanism for every module."""
    return Config() if cfg is None else cfg
