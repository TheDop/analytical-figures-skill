"""
spectra.py  -  FTIR / IR / ATR  and  PXRD / XRD  figures and the analysis behind them.

Conventions are flipped by cfg.domain (set it and forget it):
    FTIR : x = wavenumber cm-1, axis INVERTED, y = absorbance
    PXRD : x = 2theta degrees,  axis NORMAL,   y = intensity (counts / normalised)

Functions take cfg and reuse style/verify. The matplotlib is assumed-known; what
this encodes is the right axis directions, the identical-processing discipline,
and the overlaid-reps + waterfall idiom.

    load_xy(path)                       parse a 2-column xy/csv spectrum
    correct_baseline(x, y, cfg)         arPLS or rubberband -> (y_corr, baseline)
    normalize(x, y, cfg)                max / area / none
    integrate_bands(x, y, cfg)          local LINEAR baseline per window; %RSD-ready areas
    plot_overlay(traces, cfg)           reps overlaid on one axis
    plot_waterfall(groups, cfg)         vertically offset groups (the house idiom)
"""
from __future__ import annotations
import numpy as np
from . import style, verify

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
        except Exception as e:
            print(f"  [WARN] baseline: arPLS failed ({type(e).__name__}: {e}) - this trace fell back to "
                  f"rubberband; cfg.baseline no longer describes what ran for it")
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
        a = abs(_trapz(np.abs(y), x))                 # a descending axis must not flip the sign
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
    except Exception as e:
        print(f"  [WARN] smooth: Savitzky-Golay unavailable or window unworkable ({type(e).__name__}) - "
              f"{'raw trace' if deriv == 0 else 'finite-difference derivative'} used instead")
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


def band_area(x, y, lo, hi, anchor_half=2.0, line=None, clip=True):
    """THE band-area primitive: trapezoid area of y above the straight line between the two
    anchors, over lo <= x <= hi. Both `band_metric` (one band, locked anchors) and
    `integrate_bands` (declared windows, per-window or shared baseline) call this, so an
    axis-order or anchor-value bug cannot exist in one path and not the other.

    x, y need not be sorted (a descending instrument export is sorted here). The baseline
    value at each anchor is `_anchor_val`'s mean over +/- anchor_half x units, unless
    `line=(y_lo, y_hi)` supplies the values directly (the shared-envelope mode). `clip`
    drops the parts of the band that dip below the line (the house convention; pass False
    for a signed area). Returns dict(area, xs, corr): the area (NaN when the window holds
    fewer than 2 samples), the in-window x and the baseline-corrected trace."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    o = np.argsort(x, kind="stable"); x, y = x[o], y[o]
    lo, hi = sorted((float(lo), float(hi)))
    m = (x >= lo) & (x <= hi)
    if m.sum() < 2:
        return {"area": float("nan"), "xs": x[m], "corr": np.full(int(m.sum()), np.nan)}
    if line is None:
        line = (_anchor_val(x, y, lo, anchor_half), _anchor_val(x, y, hi, anchor_half))
    xs, ys = x[m], y[m]
    corr = ys - np.interp(xs, [lo, hi], [float(line[0]), float(line[1])])
    area = float(_trapz(np.clip(corr, 0, None) if clip else corr, xs))
    return {"area": area, "xs": xs, "corr": corr}


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

        cfg.integration_windows = spectra.auto_windows(ref_x, ref_y, cfg)

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
    r = band_area(x, y, a_lo, a_hi, anchor_half=anchor_half)          # the one primitive
    if r["xs"].size < 2:
        return float("nan")
    if metric == "area":
        return r["area"]
    xs, corr = r["xs"], r["corr"]
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
    collide with a neighbouring panel (verify.audit_layout flags that). Used by
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
    o = np.argsort(x, kind="stable")                 # a descending export (raw .spc: 4000 -> 650) would
    x, y = x[o], y[o]                                # flip the trapezoid sign and break np.interp
    half = 2.0 if cfg.domain == "ftir" else 0.0     # anchor averaging window (x units)
    windows = cfg.integration_windows
    mode = cfg.integration_baseline

    # shared mode draws the baseline once, across the union of all windows
    env = None
    if mode == "shared" and windows:
        env_lo = min(min(lo, hi) for lo, hi, _ in windows)
        env_hi = max(max(lo, hi) for lo, hi, _ in windows)
        em = (x >= env_lo) & (x <= env_hi)
        if em.sum() >= 2:
            env = (env_lo, env_hi, _anchor_val(x, y, env_lo, half), _anchor_val(x, y, env_hi, half))
    elif mode not in ("per_window", "shared"):
        raise ValueError(f"unknown integration_baseline '{mode}'")

    out = []
    for lo, hi, name in windows:
        a, b = sorted((lo, hi))
        if env is not None:                       # the one shared line, evaluated at this window's edges
            elo, ehi, eylo, eyhi = env
            line = tuple(np.interp([a, b], [elo, ehi], [eylo, eyhi]))
            r = band_area(x, y, a, b, line=line)
        else:                                     # local line between this window's own anchors
            r = band_area(x, y, a, b, anchor_half=half)
        out.append({"name": name, "area": r["area"], "lo": a, "hi": b})
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
    style.apply_style(cfg)
    if ax is None:
        fig, ax = style.figure(cfg)
    else:
        fig = ax.figure
    for label, (x, y) in traces.items():
        verify.check_trace(x, y, cfg, name=label)
        ax.plot(x, y, label=label)
    _apply_axis(ax, cfg)
    if len(traces) > 1:
        ax.legend(loc="best")
    return fig, ax


def plot_waterfall(groups, cfg, offset=None):
    """groups: dict label -> list of (x, y) reps. Each group is offset vertically
    (the house waterfall idiom); reps within a group are overlaid in one colour.
    PXRD reference stick lines from cfg.pxrd_reference_lines are drawn at the base."""
    style.apply_style(cfg)
    fig, ax = style.figure(cfg, height=None)
    # estimate a sensible offset from the data span if not given
    spans = [np.ptp(y) for trs in groups.values() for _, y in trs if len(y)]
    step = offset if offset is not None else (1.15 * max(spans) if spans else 1.0)
    pal = style._palette(cfg)
    edges = []                                  # (y_at_right_edge, label, colour) per group
    for gi, (label, trs) in enumerate(groups.items()):
        color = pal[gi % len(pal)]
        if not trs:
            print(f"  [WARN] waterfall: group '{label}' has no traces - skipped")
            continue
        for ri, (x, y) in enumerate(trs):
            verify.check_trace(x, y, cfg, name=f"{label}#{ri}")
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
    mode = cfg.waterfall_legend
    if len(groups) > 1 and mode == "legend":
        ax.legend(loc="best")
    elif len(groups) > 1 and mode == "edge":
        edge_labels(ax, edges)            # the house key: right-margin, per-trace
    return fig, ax
