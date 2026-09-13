"""
verify.py  -  the part no other figure skill has: verification in two tiers.

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
                                   univariate analog of chemometrics.diagnostics)

The intended loop:  draw -> render_preview -> audit_layout (fix any FAIL)
-> open the PNG with the Read tool, walk READ_IMAGE_CHECKLIST -> fix -> re-render
-> only then save_fig() the vector master (bundle.py runs the Tier-3 critic at handoff).
"""
from __future__ import annotations
import os
import re
import numpy as np

SEV = {"INFO": 0, "WARN": 1, "FAIL": 2}


class GateError(RuntimeError):
    """Raised when a FAIL gate trips under cfg.strict."""


def _resolve(findings, cfg, context=""):
    """Print findings; raise on any FAIL if cfg.strict (cfg=None -> the skill defaults)."""
    from . import config
    cfg = config.default_cfg(cfg)
    worst = max((SEV[s] for s, _ in findings), default=0)
    for sev, msg in findings:
        print(f"  [{sev}] {context}{': ' if context else ''}{msg}")
    if worst == SEV["FAIL"] and cfg.strict:
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
    # isolated single-sample spike (dead pixel / cosmic ray / zinger). For every sample, the local
    # baseline is the median of the outer ring (i-4..i-2, i+2..i+4) and the local noise its MAD.
    # A spike stands far above (or below) that baseline while BOTH neighbours stay within 25 %
    # of its excursion; a real band or Bragg peak >= ~1.3 samples wide lifts at least one
    # neighbour further than that. Domain-independent (the old second-difference rule
    # false-alarmed on a sharp FTIR band and was switched off for PXRD, where zingers occur).
    if y.size >= 9 and np.ptp(y) > 0:
        from numpy.lib.stride_tricks import sliding_window_view
        W = sliding_window_view(y, 9)
        ring = np.concatenate([W[:, 0:3], W[:, 6:9]], axis=1)
        base = np.median(ring, axis=1)
        noise = 1.4826 * np.median(np.abs(ring - base[:, None]), axis=1) + 1e-12
        centre, nb_hi, nb_lo = W[:, 4], np.maximum(W[:, 3], W[:, 5]), np.minimum(W[:, 3], W[:, 5])
        floor = 0.1 * np.ptp(y)
        up, dn = centre - base, base - centre
        spike = (((up > floor) & (up > 8 * noise) & (nb_hi - base < 0.25 * up))
                 | ((dn > floor) & (dn > 8 * noise) & (base - nb_lo < 0.25 * dn)))
        if spike.any():
            f.append(("WARN", f"isolated single-sample spike at {int(np.argmax(spike)) + 4} - "
                              f"check for a dead pixel/cosmic ray"))
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


_VER_NO_SCIPY_WARNED = False


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
        from statistics import NormalDist              # scipy absent: exact z, and say so (once)
        z = float(NormalDist().inv_cdf(min(max(0.5 + p / 2, 1e-12), 1 - 1e-12)))
        global _VER_NO_SCIPY_WARNED
        if not _VER_NO_SCIPY_WARNED:
            _VER_NO_SCIPY_WARNED = True
            print(f"  [WARN] summarize: scipy unavailable - normal quantile z={z:.3f} used instead of "
                  f"t at n-1={n - 1} dof (CI too narrow at small n); further calls stay silent")
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
    multivariate `chemometrics.diagnostics` gate: it turns an ad-hoc "flag-don't-hide"
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
        from . import config
        n_mads = config.default_cfg(cfg).outlier_mad_n
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
    _resolve(f, cfg, context)                   # only INFO/WARN here, so cfg=None prints the same
    return {"raw": p, "adjusted": adj, "reject": reject, "method": method, "alpha": alpha,
            "labels": list(labels) if labels is not None else None,
            "caption": f"{method}-adjusted p (family of {m}, alpha={alpha})"}


# ===================================================== TIER 1b: baseline robustness
def _locked_metric(cfg, center, reference):
    """Build a band_metric(x, y, method) closure that, for the anchor-based metrics
    (area/height), uses anchors detected ONCE on `reference` and locked - never
    re-detected per spectrum. deriv2 is anchor-free. If `reference` is None the
    area/height anchors fall back to per-spectrum detection (discouraged: drifts)."""
    from . import spectra                       # lazy: spectra imports verify
    locked = None
    if reference is not None:
        rx, ry = reference
        locked = spectra.find_anchors(rx, ry, center, gap=cfg.anchor_gap,
                                      maxhw=cfg.anchor_maxhw, smooth=cfg.anchor_smooth)
    def metric(x, y, mth):
        if mth == "deriv2":
            return spectra.band_metric(x, y, "deriv2", center=center,
                                       deriv_smooth=cfg.deriv_smooth)
        return spectra.band_metric(x, y, mth, anchors=locked, center=center,
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
    spectra.

    reps      : iterable of (x, y) replicate spectra of the SAME sample.
    center    : band centre (cm-1), handed to spectra.band_metric.
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
    from . import calibration                   # lazy: imports verify
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
        model = calibration.fit(np.array(concs, float), np.array(sigs, float), cfg)
        ll = calibration.lod_loq(model, cfg)
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
    them with style.add_panel_labels)."""
    import matplotlib.text as mtext
    from . import config
    cfg = config.default_cfg(cfg)
    f = []
    tol = cfg.tick_overlap_tol_px
    clip_tol = cfg.clip_tol_px
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
        esc = _escaping_lines(ax, cfg.escape_tol_frac)
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
            import re as _re
            low = lab.lower()
            def _cued(c):                       # word-aware for plain words ('ratio' must not match 'concentRATIOn')
                return (_re.search(r"(?<![a-z])" + _re.escape(c) + r"(?![a-z])", low) is not None
                        if c.isalpha() else c in low)
            if has_numeric and lab and not any(_cued(c) for c in _UNIT_CUES):
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
                          f"style.finalize_figure(fig) or export bbox_inches='tight'"))

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
                              f"call style.add_panel_labels(fig, cfg)"))
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
