"""
style.py  -  the house-style LOCK.

Three jobs, all about consistency rather than "how to use matplotlib":
    apply_style(cfg)        install rcParams: journal widths, colour-blind palette,
                            font floor at 6 pt, vector-friendly defaults.
    figure(cfg, ...)        a Figure sized at FINAL print dimensions (never rescale later).
    save_fig(fig, base, cfg) export vector master + raster proof; refuses JPEG.

Degrades gracefully: if scienceplots is absent it falls back to a built-in
preset that encodes the same intent (no LaTeX requirement).

pyplot is imported on FIRST USE (`_plt()`), not at module load: `spectra`, `verify`,
`calibration` and `charts` import this module, so a CSV-only standalone or `cli.py` would
otherwise pay the ~0.3 s matplotlib import for numbers that never reach a figure.
"""
from __future__ import annotations
import os
from .config import SKILL_VERSION, default_cfg   # dropped when bundled; then resolve to the inlined globals


def _plt():
    """matplotlib.pyplot, imported lazily (see the module docstring)."""
    import matplotlib.pyplot as plt
    return plt

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
    red = cfg.redundancy                        # house default "none": spectra solid unless opted out
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
    """Install the house style. Call once before plotting."""
    # Prefer scienceplots' no-latex style; fall back to the hand-rolled preset.
    plt = _plt()
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")          # scienceplots' own deprecation chatter
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
    fig, axes = _plt().subplots(nrows, ncols, figsize=(w, h),
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
    desc = (cfg.description or "").strip()
    ver = f"analytical-figures v{SKILL_VERSION}"
    import matplotlib
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
    cfg = default_cfg(cfg)
    journal = journal or cfg.journal
    column = column or cfg.column
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
    lines = ["# analytical-figures house style (auto-exported by style.export_mplstyle;",
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
# Rewritten here to read the house cfg and to tag each label so verify.audit_layout
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
    carries right-margin edge labels (spectra.edge_labels) that would otherwise spill
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
            style = default_cfg(cfg).journal
        fmt = _PANEL_FMT.get(style, _PANEL_FMT["general"])
        labels = [fmt(s) for s in _letter_sequence(len(axs))]
    elif len(labels) < len(axs):
        raise ValueError(f"{len(labels)} labels for {len(axs)} panels")
    if fontsize is None:
        fontsize = _plt().rcParams.get("axes.labelsize", 9)
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
        fontsize = _plt().rcParams.get("axes.labelsize", 9)
    va = "top" if "upper" in loc else "bottom"
    ha = "right" if "right" in loc else "left"
    x = (1 - pad) if ha == "right" else pad
    y = (1 - pad) if va == "top" else pad
    t = ax.text(x, y, label, transform=ax.transAxes, ha=ha, va=va, zorder=6,
                fontsize=fontsize, fontweight=fontweight, color=color)
    t.set_gid("panel-label")
    return t
