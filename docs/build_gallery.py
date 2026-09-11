#!/usr/bin/env python3
"""
build_gallery.py — the README gallery, the GitHub Pages gallery and the social-preview card,
all from ONE manifest, using only figures the shipped code produces on public data and public
CIFs. Nothing here is drawn by hand and nothing is synthetic.

    python docs/build_gallery.py            # re-run the producers (~1-2 min), then compose
    python docs/build_gallery.py --no-run   # compose from the producers' existing out/ figures

Outputs (all under docs/):
    gallery/tiles/<id>.png    square thumbnails (the README table + the Pages grid)
    gallery/full/<id>.png     full-resolution copies (the click-through targets)
    gallery/social.png        2:1 mosaic, darkened, title centred (upload as the GitHub social preview)
    index.html                the Pages gallery: hover a tile for its placard (plain HTML + CSS)
    ../README.md              the tile table, written between <!-- gallery:start/end --> markers

Data behind the tiles:
    Wang 2023 pectin ATR-FTIR spectra, Mendeley Data 10.17632/gkwbp3wc49.1, CC BY 4.0
        (vendored: skill_validation/ftir_integration/data/pectin/)
    Bansal, Singh & Kaur 2021, BMC Chemistry 15:27 — the published peak-area tables
    Crystallography Open Database CIFs (aspirin 7247819, urea 1008785, lactose 2206486,
        flufenamic acid 4118081)

Rebuild only when a producer or a tile changes — the PNGs are committed, and every rebuild adds
to the repository history.
"""
import os
import sys
import glob
import html
import shutil
import subprocess
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                       # the skill root
GAL = os.path.join(HERE, "gallery")
TILES, FULL = os.path.join(GAL, "tiles"), os.path.join(GAL, "full")
INLINE = os.path.join(GAL, "_inline")               # scratch for the tiles rendered here
REPO_URL = "https://github.com/TheDop/analytical-figures-skill"
PAGES_URL = "https://thedop.github.io/analytical-figures-skill/"

V = "skill_validation"
CR, FT = f"{V}/crystal", f"{V}/ftir_integration"
PECTIN = os.path.join(ROOT, FT, "data", "pectin")
CIFS = os.path.join(ROOT, CR, "cifs")

WANG = "Wang 2023 pectin ATR-FTIR spectra, Mendeley Data 10.17632/gkwbp3wc49.1 (CC BY 4.0)"
BMC_SRC = "Bansal, Singh & Kaur 2021, BMC Chemistry 15:27, 10.1186/s13065-021-00752-3 (published peak-area tables)"
COD3 = "Aspirin COD 7247819, urea COD 1008785, lactose COD 2206486 (Crystallography Open Database)"
ASP = "Aspirin, COD 7247819"

# ---------------------------------------------------------------- producers (shipped scripts)
PRODUCERS = [
    f"{FT}/validate_ftir_integration.py",
    f"{FT}/groundtruth_integration_test.py",
    f"{CR}/crystal_demo.py",
    f"{CR}/test_realism_cif.py",
    f"{CR}/test_ellipsoid.py",
    f"{CR}/test_pyvista.py",          # optional: skipped when pyvista is absent
]

# ---------------------------------------------------------------- the manifest
# id, source PNG (relative to the skill root; INLINE:<name> for tiles rendered below),
# caption, data provenance, producing script
MANIFEST = [
    ("composite_method", "INLINE:composite_method",
     "Method-comparison page figure: four calibrations, each over its own residual strip, and the figures of merit with the direction of good in every title",
     BMC_SRC, "docs/build_gallery.py → scripts/calibration.fit + lod_loq, style.panel_letter"),
    ("composite_specificity", "INLINE:composite_specificity",
     "Specificity composite: the full spectra, the carbonyl window with both integration windows on one shared baseline, and the band-area ratio it yields against DM",
     WANG, "docs/build_gallery.py → scripts/spectra.integrate_bands (shared baseline), calibration.fit"),
    ("composite_pxrd", "INLINE:composite_pxrd",
     "Phase-ID composite: stacked calculated patterns over the full range and the low-angle window, one right-margin key serving both panels",
     COD3, "docs/build_gallery.py → scripts/crystal_pxrd.calc_pattern, spectra.edge_labels"),
    ("pls_diag", "INLINE:pls_diagnostics",
     "Leakage-safe PLS diagnostics on six real standards: RMSECV per preprocessing with the parsimony pick ringed, coefficients, scores, loadings",
     WANG + " — the six calibration standards, leave-one-out", "docs/build_gallery.py → scripts/chemometrics.diagnostics_figure"),
    ("calib_pectin", f"{FT}/out/pectin_calibration_published.png",
     "Calibration with confidence and prediction bands and the mandatory residual panel",
     WANG + " — the authors' published band areas", f"{FT}/validate_ftir_integration.py"),
    ("waterfall", "INLINE:waterfall",
     "Waterfall of the six calibration standards, right-margin keys instead of a legend",
     WANG, "docs/build_gallery.py → scripts/spectra.plot_waterfall"),
    ("ddsimca", "INLINE:ddsimca",
     "DD-SIMCA acceptance plot of eighteen replicate sample spectra: all inside the acceptance boundary, the extreme and outlier limits drawn",
     WANG + " — the eighteen sample spectra (6 preparations × 3 replicates), SNV, 2 PCs", "docs/build_gallery.py → scripts/chemometrics.diagnostics + plot_influence"),
    ("integration", "INLINE:integration",
     "Two overlapping bands on one shared baseline with the vertical drop at the window boundary — the drop-perpendicular method the validation showed is required",
     WANG + " — the DM 70.5 % standard", "docs/build_gallery.py → scripts/spectra.integrate_bands"),
    ("ortep", f"{CR}/out/crystal_demo/structure_ellipsoid.png",
     "ORTEP displacement ellipsoids at 50 % probability, deterministic PCA camera",
     ASP + " (296 K)", f"{CR}/crystal_demo.py"),
    ("calib_bmc", f"{FT}/out/bmc_dox_transmittance.png",
     "Doxorubicin carbonyl-area calibration; LOD/LOQ from the residual SD, and the sigma is named",
     BMC_SRC, f"{FT}/validate_ftir_integration.py"),
    ("groundtruth", f"{FT}/out/groundtruth_shared_vs_perwindow.png",
     "Shared vs per-window baseline against a closed-form answer: the per-window pedestal bias, measured",
     "Analytic Gaussian bands with known areas, checked against scipy.integrate.quad", f"{FT}/groundtruth_integration_test.py"),
    ("pxrd_calc", f"{CR}/out/crystal_demo/pxrd.png",
     "Calculated PXRD from a CIF at Cu Kα, top peak cross-checked against pymatgen",
     ASP, f"{CR}/crystal_demo.py"),
    ("pxrd_realism", f"{CR}/out/_realism_aspirin.png",
     "Lab-realistic pattern: Cu Kα₁/Kα₂ doublet and Caglioti broadening on the Dans reflection list",
     ASP, f"{CR}/test_realism_cif.py"),
    ("packing", f"{CR}/out/crystal_demo/packing.png",
     "Packing diagram down a with hydrogen bonds and the cell box", ASP, f"{CR}/crystal_demo.py"),
    ("hbond_env", f"{CR}/out/crystal_demo/hbond_environment.png",
     "Hydrogen-bond environment of the asymmetric unit, symmetry mates superscripted", ASP, f"{CR}/crystal_demo.py"),
    ("unit_cell", f"{CR}/out/crystal_demo/unit_cell.png",
     "Unit-cell contents with the cell box and a/b/c", ASP, f"{CR}/crystal_demo.py"),
    ("flufenamic", f"{CR}/out/_ellipsoid_disorder-groups_flufenamic.png",
     "Z′ = 3 with CF₃ disorder groups honoured: no bonds between mutually exclusive sites",
     "Flufenamic acid, COD 4118081", f"{CR}/test_ellipsoid.py"),
    ("pyvista", f"{CR}/out/_pyvista_aspirin.png",
     "PyVista/VTK still: real depth buffer, ambient occlusion, orthographic camera from the same orientation engine",
     ASP, f"{CR}/test_pyvista.py"),
]

TITLE = "analytical-figures"
TAGLINE = "Publication-grade figures and the gated analysis behind them, for analytical chemistry."
FOOT = "Claude Code skill  ·  MIT  ·  github.com/TheDop/analytical-figures-skill"


# ---------------------------------------------------------------- 1. run the producers
def run_producers():
    for rel in PRODUCERS:
        path = os.path.join(ROOT, rel)
        if rel.endswith("test_pyvista.py") and importlib.util.find_spec("pyvista") is None:
            print(f"  skip  {rel} (pyvista not installed)")
            continue
        r = subprocess.run([sys.executable, path], capture_output=True, text=True)
        print(f"  {'ok  ' if r.returncode == 0 else 'FAIL'}  {rel}")
        if r.returncode != 0:
            print(r.stdout[-1500:], r.stderr[-1000:])


# ---------------------------------------------------------------- 2. tiles rendered here
def _load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def render_inline():
    """Tiles drawn with the skill's own modules on the vendored pectin spectra, the BMC 2021
    area tables and the COD CIFs. Returns {INLINE name: PNG path}."""
    import numpy as np
    import matplotlib.pyplot as plt
    from dataclasses import replace
    sys.path.insert(0, ROOT)
    from scripts.config import Config
    from scripts import style, spectra, chemometrics, calibration

    # the transcribed published numbers live in ONE place: the validation script
    val = _load_module(os.path.join(ROOT, FT, "validate_ftir_integration.py"), "vfi")
    PECTIN_DM, BMC = val.PECTIN, val.BMC

    os.makedirs(INLINE, exist_ok=True)
    cfg = Config(output_dir=INLINE, formats=("png",), journal="nature", column="single",
                 domain="ftir", strict=False, description="gallery tile")
    cfg2 = replace(cfg, column="double")
    pal = style._palette(cfg)
    out = {}

    def save(fig, name, c=cfg):
        style.save_fig(fig, os.path.join(INLINE, name), c)
        out[name] = os.path.join(INLINE, name + ".png")

    # ---- the pectin standards (real ATR-FTIR, DM known) and samples (6 preparations × 3 reps)
    std = {}
    for s in PECTIN_DM:
        x, y = spectra.load_xy(os.path.join(PECTIN, "calibration", f"{s}.csv"))
        x, y = np.asarray(x, float), np.asarray(y, float)
        o = np.argsort(x)
        std[s] = (x[o], y[o])                                     # ascending
    dm_pct = {s: PECTIN_DM[s][0] * 100.0 for s in PECTIN_DM}
    x0 = std["DM3"][0]
    samples = {}
    for p in sorted(glob.glob(os.path.join(PECTIN, "samples", "*.csv"))):
        x, y = spectra.load_xy(p)
        x, y = np.asarray(x, float), np.asarray(y, float)
        o = np.argsort(x)
        samples[os.path.splitext(os.path.basename(p))[0]] = np.interp(x0, x[o], y[o])
    fp = (x0 >= 800) & (x0 <= 1900)                               # fingerprint + carbonyl window
    xd = x0[fp][::-1]                                             # descending, the FTIR convention

    # (a) waterfall of the six standards, right-margin keys
    wf = {f"DM {dm_pct[s]:.0f} %": [(xd, std[s][1][fp][::-1])] for s in PECTIN_DM}
    fig, ax = spectra.plot_waterfall(wf, cfg)
    save(fig, "waterfall")

    # (b) two overlapping bands on one SHARED baseline (vertical drop at 1700 cm-1)
    icfg = replace(cfg, integration_baseline="shared",
                   integration_windows=[(1500, 1700, "carbox"), (1700, 1800, "ester")])
    xs, ys = std["DM70.5"]
    lo, hi, split = 1500.0, 1800.0, 1700.0
    ylo, yhi = np.interp(lo, xs, ys), np.interp(hi, xs, ys)
    bands = {b["name"]: b["area"] for b in spectra.integrate_bands(xs, ys, icfg)}
    style.apply_style(cfg)
    fig, ax = style.figure(cfg)
    ax.plot(xs, ys, lw=0.9, color="k")
    for k, (a, b) in enumerate(((lo, split), (split, hi))):
        m = (xs >= a) & (xs <= b)
        base = np.interp(xs[m], [lo, hi], [ylo, yhi])
        ax.fill_between(xs[m], base, ys[m], alpha=0.35, lw=0, color=pal[k])
    ax.plot([lo, hi], [ylo, yhi], color="0.25", lw=0.8, ls="--")
    ax.plot([split, split], [np.interp(split, [lo, hi], [ylo, yhi]), np.interp(split, xs, ys)],
            color="0.25", lw=0.8, ls="--")
    ax.plot([lo, hi], [ylo, yhi], "o", color="0.25", ms=3)
    spectra._apply_axis(ax, cfg)
    ax.set_xlim(1850, 1450)
    ratio_705 = bands["ester"] / (bands["ester"] + bands["carbox"])
    ax.set_title(f"shared baseline {lo:.0f}–{hi:.0f} cm⁻¹, drop at {split:.0f}: "
                 f"ester/(ester+carboxylate) = {ratio_705:.3f}", fontsize="small")
    save(fig, "integration")

    # (c) PLS diagnostics on the six standards (leave-one-out)
    X = np.vstack([std[s][1][fp] for s in PECTIN_DM])
    y = np.array([dm_pct[s] for s in PECTIN_DM])
    fig, (pa, pb, pc, pd), info = chemometrics.diagnostics_figure(X, y, x0[fp], cfg2, basename="pls_diagnostics")
    cb = getattr(pc.collections[0], "colorbar", None)          # six labelled standards need no colour bar
    if cb is not None:
        cb.remove()                                           # Colorbar.remove also clears the axes bookkeeping
    for ax_ in [ax_ for ax_ in fig.axes if ax_ not in (pa, pb, pc, pd)]:
        ax_.remove()
    pc._colorbars = []                                        # belt and braces for the layout engine
    for (sx, sy), dm in zip(pc.collections[0].get_offsets(), y):
        pc.annotate(f"DM {dm:.0f} %", (sx, sy), xytext=(4, 4), textcoords="offset points", fontsize="x-small")
    style.finalize_figure(fig, wspace=0.08, hspace=0.16)
    save(fig, "pls_diagnostics", cfg2)

    # (d) DD-SIMCA on the eighteen replicate sample spectra
    names = list(samples)
    Xs = np.vstack([samples[n][fp] for n in names])
    diag = chemometrics.diagnostics(Xs, cfg, n_components=2, pre="snv")
    style.apply_style(cfg)
    fig, ax = style.figure(cfg)
    chemometrics.plot_influence(ax, diag, cfg, labels=names)
    save(fig, "ddsimca")

    # ---- composite A: method-comparison page figure (BMC 2021 — four calibrations + FoM row)
    style.apply_style(cfg2)
    fig = plt.figure(figsize=(7.0, 6.2), layout="constrained")
    fig.get_layout_engine().set(w_pad=0.04, h_pad=0.04)
    top, bot = fig.subfigures(2, 1, height_ratios=[3.0, 1.15], hspace=0.05)
    outer = top.add_gridspec(2, 2, hspace=0.20, wspace=0.16)
    models, charts = {}, []
    for i, (name, (xs_, ys_, _pub)) in enumerate(BMC.items()):
        r, c = divmod(i, 2)
        cell = outer[r, c].subgridspec(2, 1, height_ratios=[3.2, 1], hspace=0.05)
        a0 = top.add_subplot(cell[0]); a1 = top.add_subplot(cell[1], sharex=a0)
        m = calibration.fit(xs_, ys_, cfg2); models[name] = m
        xx = np.linspace(*m["x_range"], 100)
        se_mean = m["s_resid"] * np.sqrt(1.0 / m["n"] + (xx - m["xbar"]) ** 2 / m["Sxx"])
        se_pred = m["s_resid"] * np.sqrt(1.0 + 1.0 / m["n"] + (xx - m["xbar"]) ** 2 / m["Sxx"])
        yy = m["slope"] * xx + m["intercept"]
        a0.fill_between(xx, yy - m["t"] * se_mean, yy + m["t"] * se_mean, alpha=0.25, lw=0)
        a0.plot(xx, yy - m["t"] * se_pred, lw=0.6, ls="--", color="0.5")
        a0.plot(xx, yy + m["t"] * se_pred, lw=0.6, ls="--", color="0.5")
        a0.plot(xx, yy, lw=1.0)
        a0.scatter(m["x"], m["y"], zorder=3, s=12)
        a0.set_title(name, fontsize="small")
        a0.set_ylabel("peak area / mm²")
        a0.annotate(f"$R^2$={m['r2']:.4f}", xy=(0.96, 0.06), xycoords="axes fraction", ha="right", va="bottom", fontsize="x-small")
        plt.setp(a0.get_xticklabels(), visible=False)
        a1.axhline(0, color="0.6", lw=0.6)
        a1.scatter(m["x"], m["resid"], zorder=3, s=12)
        a1.set_ylabel("resid.")
        a1.set_xlabel("analyte / % w/w")
        charts.append(a0)
    axf = bot.subplots(1, 3)
    short = [n.replace(" transmittance", "\ntrans.").replace(" reflectance", "\nrefl.") for n in BMC]
    fom = {
        "Linearity $R^2$\n(higher = better)": [models[n]["r2"] for n in BMC],
        "LOD / % w/w\n(lower = better)": [calibration.lod_loq(models[n], cfg2)["lod"] for n in BMC],
        "LOQ / % w/w\n(lower = better)": [calibration.lod_loq(models[n], cfg2)["loq"] for n in BMC],
    }
    for k, (ttl, vals) in enumerate(fom.items()):
        axf[k].bar(range(len(vals)), vals, color=[pal[j % len(pal)] for j in range(len(vals))],
                   edgecolor="k", linewidth=0.5)
        axf[k].set_xticks(range(len(vals))); axf[k].set_xticklabels(short, fontsize="x-small")
        axf[k].set_title(ttl, fontsize="small")
        if "R^2" in ttl:
            axf[k].set_ylim(0.98, 1.0)
        else:
            axf[k].set_ylim(0, max(vals) * 1.28)
        charts.append(axf[k])
    fig.canvas.draw()                                             # settle the nested layout first
    style.add_panel_labels(fig, cfg2, axes=charts, x_offset_pt="auto")
    save(fig, "composite_method", cfg2)

    # ---- composite B: specificity (full spectra · carbonyl window · ratio vs DM)
    style.apply_style(cfg2)
    fig, axes = style.figure(cfg2, 2, 2, height_ratios=[1.0, 1.15])
    gs = axes[0, 0].get_gridspec()
    for ax_ in axes[0]:
        ax_.remove()
    a = fig.add_subplot(gs[0, :])
    b, c = axes[1, 0], axes[1, 1]
    for k, s in enumerate(PECTIN_DM):
        xs_, ys_ = std[s]
        a.plot(xs_, ys_, lw=0.6, color=pal[k % len(pal)])
        w = (xs_ >= 1450) & (xs_ <= 1850)
        b.plot(xs_[w], ys_[w], lw=0.8, color=pal[k % len(pal)])
    spectra._apply_axis(a, cfg2); a.set_ylabel("Absorbance")
    b.axvspan(1500, 1700, alpha=0.12, color=pal[0], lw=0)
    b.axvspan(1700, 1800, alpha=0.12, color=pal[1], lw=0)
    b.axvline(1700, color="0.3", lw=0.6, ls="--")
    spectra._apply_axis(b, cfg2); b.set_xlim(1850, 1450); b.set_ylabel("")
    b.set_title("carboxylate 1500–1700 · ester 1700–1800", fontsize="small")
    ratio = []
    for s in PECTIN_DM:
        xs_, ys_ = std[s]
        bb = {q["name"]: q["area"] for q in spectra.integrate_bands(xs_, ys_, icfg)}
        ratio.append(bb["ester"] / (bb["ester"] + bb["carbox"]))
    ratio = np.array(ratio)
    m = calibration.fit(y, ratio, cfg2)
    xx = np.linspace(*m["x_range"], 50)
    c.plot(xx, m["slope"] * xx + m["intercept"], lw=1.0, color="0.3")
    c.scatter(y, ratio, zorder=3, s=14, color=[pal[k % len(pal)] for k in range(len(y))])
    c.set_xlabel("DM / %"); c.set_ylabel("band-area ratio  I")
    c.annotate(f"$R^2$={m['r2']:.3f}", xy=(0.96, 0.06), xycoords="axes fraction", ha="right", va="bottom", fontsize="x-small")
    style.finalize_figure(fig, wspace=0.10, hspace=0.10)
    style.add_panel_labels(fig, cfg2, axes=[a, b, c], x_offset_pt="auto")
    save(fig, "composite_specificity", cfg2)

    # ---- composite C: PXRD phase ID (full range · low-angle zoom · one key) — needs gemmi + Dans
    try:
        from scripts import crystal_engine as ce, crystal_pxrd as cp
        pcfg = replace(cfg2, domain="pxrd", pxrd_two_theta_min=5, pxrd_two_theta_max=50, pxrd_wavelength=1.540598)
        pats = []
        for label, cif in (("aspirin", "monoclinic_aspirin__COD7247819.cif"),
                           ("urea", "specialpos_urea__COD1008785.cif"),
                           ("lactose", "lactose__COD2206486.cif")):
            ccfg = replace(pcfg, cif_path=os.path.join(CIFS, cif))
            pats.append((label, cp.calc_pattern(ce.load(ccfg), ccfg)))
        style.apply_style(pcfg)
        fig, (a, b) = style.figure(pcfg, 1, 2, width_ratios=[2.0, 1.0], sharey=True)
        items = []
        for k, (label, pat) in enumerate(pats):
            tt, ii = np.asarray(pat["two_theta"], float), np.asarray(pat["intensity"], float)
            if not np.isfinite(ii).all() or ii.max() <= 0:
                raise RuntimeError(f"empty calculated pattern for {label}")
            ii = ii / ii.max()
            off = k * 1.08
            col = pal[k % len(pal)]
            a.plot(tt, ii + off, lw=0.7, color=col)
            z = (tt >= 5) & (tt <= 20)
            b.plot(tt[z], ii[z] + off, lw=0.7, color=col)
            items.append((off + 0.5, label, col))
        spectra._apply_axis(a, pcfg); spectra._apply_axis(b, pcfg)
        a.set_ylabel("Intensity (normalised, offset)"); b.set_ylabel("")
        a.set_xlim(5, 50); b.set_xlim(5, 20)
        a.set_yticks([]); b.tick_params(labelleft=False)
        b.set_title("low-angle window", fontsize="small")
        spectra.edge_labels(b, items)
        style.finalize_figure(fig, wspace=0.10)
        style.add_panel_labels(fig, pcfg)
        save(fig, "composite_pxrd", pcfg)
    except ImportError as e:
        print(f"  composite_pxrd skipped ({e})")
    return out


# ---------------------------------------------------------------- 3. compose
def _font(name, size):
    from PIL import ImageFont
    import matplotlib
    p = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf", name)
    return ImageFont.truetype(p, size)


def _contain(im, size, pad_frac=0.03, bg="white"):
    from PIL import Image
    im = im.convert("RGB")
    pad = int(size * pad_frac)
    inner = size - 2 * pad
    scale = min(inner / im.width, inner / im.height)
    w, h = max(1, round(im.width * scale)), max(1, round(im.height * scale))
    im = im.resize((w, h), Image.LANCZOS)
    canvas = Image.new("RGB", (size, size), bg)
    canvas.paste(im, ((size - w) // 2, (size - h) // 2))
    return canvas


def compose(entries):
    from PIL import Image, ImageDraw
    for d in (TILES, FULL):
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
    tiles = {}
    for e in entries:
        im = Image.open(e["png"])
        cap = 1400                                    # long side capped so the repository stays small
        if max(im.size) > cap:
            s = cap / max(im.size)
            full = im.convert("RGB").resize((round(im.width * s), round(im.height * s)), Image.LANCZOS)
        else:
            full = im.convert("RGB")
        full.save(os.path.join(FULL, e["id"] + ".png"), optimize=True)
        t = _contain(im, 600)
        t.save(os.path.join(TILES, e["id"] + ".png"), optimize=True)
        tiles[e["id"]] = t

    # social card: full-bleed 6 x 3 mosaic, darkened, the title centred over it
    gap, tile, cols, rows = 12, 384, 6, 3
    W, H = cols * tile + (cols + 1) * gap, rows * tile + (rows + 1) * gap
    card = Image.new("RGB", (W, H), "#e6e9ee")
    slots = [(r, c) for r in range(rows) for c in range(cols)]
    for (r, c), e in zip(slots, entries):
        card.paste(tiles[e["id"]].resize((tile, tile), Image.LANCZOS),
                   (gap + c * (tile + gap), gap + r * (tile + gap)))
    card = Image.alpha_composite(card.convert("RGBA"),
                                 Image.new("RGBA", (W, H), (12, 14, 20, 150)))
    import numpy as np
    yy = np.arange(H)
    prof = np.clip(np.cos((yy - H / 2) / (0.32 * H) * np.pi / 2), 0, 1) ** 1.5   # 1 at centre -> 0 at ±32 %
    band = np.zeros((H, W, 4), dtype=np.uint8)
    band[..., :3] = (8, 10, 16)
    band[..., 3] = (prof * 120)[:, None].astype(np.uint8)
    card = Image.alpha_composite(card, Image.fromarray(band, "RGBA")).convert("RGB")
    d = ImageDraw.Draw(card)
    f_title, f_tag, f_foot = (_font("DejaVuSans-Bold.ttf", 118), _font("DejaVuSans.ttf", 36),
                              _font("DejaVuSans.ttf", 26))

    def centred(text, font, y, fill):
        x0, y0, x1, y1 = d.textbbox((0, 0), text, font=font)
        d.text(((W - (x1 - x0)) / 2 - x0, y), text, font=font, fill=fill)
        return y1 - y0

    y = H / 2 - 150
    h = centred(TITLE, f_title, y, "white")
    d.rectangle([W / 2 - 90, y + h + 34, W / 2 + 90, y + h + 40], fill="#4f8dff")
    centred(TAGLINE, f_tag, y + h + 74, "#e8ebf0")
    centred(FOOT, f_foot, y + h + 140, "#aab2bf")
    tmp = os.path.join(GAL, "_social_tmp.png")
    card.save(tmp, optimize=True)
    os.replace(tmp, os.path.join(GAL, "social.png"))
    print(f"  social.png {W}x{H}, {min(len(entries), len(slots))} tiles")


# ---------------------------------------------------------------- 4. README table + Pages
def write_readme_table(entries, cols=6):
    p = os.path.join(ROOT, "README.md")
    s = open(p, encoding="utf-8").read()
    a, b = "<!-- gallery:start -->", "<!-- gallery:end -->"
    if a not in s or b not in s:
        print("  README markers not found; table not written")
        return
    rows = []
    for k in range(0, len(entries), cols):
        cells = []
        for e in entries[k:k + cols]:
            tip = html.escape(f"{e['caption']} — data: {e['source']}", quote=True)
            cells.append(f'<td align="center"><a href="docs/gallery/full/{e["id"]}.png">'
                         f'<img src="docs/gallery/tiles/{e["id"]}.png" width="140" alt="{html.escape(e["caption"], quote=True)}" '
                         f'title="{tip}"></a></td>')
        rows.append("<tr>" + "".join(cells) + "</tr>")
    table = ("<table>\n" + "\n".join(rows) + "\n</table>\n\n"
             "<sub>Hover a tile for what it is and where the data came from; click for full size. Every tile is "
             f"produced by <code>docs/build_gallery.py</code> from public data and public CIFs — nothing is drawn "
             f"by hand and nothing is synthetic. Interactive version with placards: <a href=\"{PAGES_URL}\">{PAGES_URL}</a></sub>")
    s = s[: s.index(a) + len(a)] + "\n" + table + "\n" + s[s.index(b):]
    open(p, "w", encoding="utf-8", newline="\n").write(s)
    print(f"  README table: {len(entries)} tiles")


def write_pages(entries):
    cards = []
    for e in entries:
        cards.append(
            f'<a class="tile" href="gallery/full/{e["id"]}.png">'
            f'<img src="gallery/tiles/{e["id"]}.png" alt="{html.escape(e["caption"], quote=True)}" loading="lazy">'
            f'<div class="placard"><b>{html.escape(e["caption"])}</b>'
            f'<span class="src">Data: {html.escape(e["source"])}</span>'
            f'<span class="scr">Made by <code>{html.escape(e["script"])}</code></span></div></a>')
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{TITLE} — gallery</title>
<meta property="og:title" content="{TITLE}">
<meta property="og:description" content="{TAGLINE}">
<meta property="og:image" content="{PAGES_URL}gallery/social.png">
<style>
:root{{--bg:#f3f4f6;--card:#fff;--ink:#1b1f24;--muted:#57606a;--accent:#0b5fff}}
*{{box-sizing:border-box}}
body{{margin:0;font:16px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;background:var(--bg);color:var(--ink)}}
header,footer,main{{max-width:1240px;margin:0 auto;padding:0 20px}}
header{{padding-top:44px}}
h1{{font-size:2.1rem;margin:0;letter-spacing:-.01em}}
.tag{{color:var(--muted);margin:.3rem 0 0;max-width:60ch}}
.hint{{color:var(--muted);font-size:.9rem;margin:.6rem 0 0}}
main{{padding:22px 20px 8px;display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:16px}}
.tile{{position:relative;display:block;background:var(--card);border-radius:10px;overflow:hidden;aspect-ratio:1;
  box-shadow:0 1px 3px rgba(0,0,0,.08);transition:transform .15s ease,box-shadow .15s ease;color:inherit;text-decoration:none}}
.tile img{{width:100%;height:100%;object-fit:contain;display:block}}
.tile:hover,.tile:focus-visible{{transform:translateY(-3px);box-shadow:0 10px 26px rgba(0,0,0,.16);outline:none}}
.placard{{position:absolute;left:0;right:0;bottom:0;background:rgba(255,255,255,.97);border-top:3px solid var(--accent);
  padding:10px 12px 12px;transform:translateY(102%);transition:transform .18s ease;font-size:.84rem;line-height:1.35}}
.tile:hover .placard,.tile:focus-visible .placard{{transform:none}}
.placard b{{display:block;margin-bottom:4px}}
.placard .src,.placard .scr{{display:block;color:var(--muted);font-size:.76rem;margin-top:2px}}
.placard code{{font-size:.72rem;background:#eef0f3;padding:1px 4px;border-radius:4px}}
@media (hover:none){{.placard{{position:static;transform:none}}.tile{{aspect-ratio:auto}}.tile img{{aspect-ratio:1}}}}
footer{{padding:16px 20px 40px;color:var(--muted);font-size:.86rem}}
footer a{{color:var(--accent)}}
</style>
</head>
<body>
<header>
  <h1>{TITLE}</h1>
  <p class="tag">{TAGLINE} A <a href="{REPO_URL}">Claude Code skill</a>.</p>
  <p class="hint">Hover a tile for what it is and where the data came from; click for full size.</p>
</header>
<main>
{chr(10).join(cards)}
</main>
<footer>
  Every tile is produced by <code>docs/build_gallery.py</code> from public data and public CIFs — nothing is
  drawn by hand and nothing is synthetic. <a href="{REPO_URL}">Repository</a> · MIT licence.
</footer>
</body>
</html>
"""
    open(os.path.join(HERE, "index.html"), "w", encoding="utf-8", newline="\n").write(page)
    print("  index.html written")


# ---------------------------------------------------------------- main
def main():
    run = "--no-run" not in sys.argv
    os.makedirs(GAL, exist_ok=True)
    if run:
        print("producers:")
        run_producers()
    print("inline tiles:")
    inline = render_inline()
    entries = []
    for id_, src, caption, source, script in MANIFEST:
        png = inline.get(src[7:]) if src.startswith("INLINE:") else os.path.join(ROOT, src)
        if not png or not os.path.exists(png):
            print(f"  missing  {id_}: {src}")
            continue
        entries.append(dict(id=id_, png=png, caption=caption, source=source, script=script))
    print("compose:")
    compose(entries)
    write_readme_table(entries)
    write_pages(entries)
    shutil.rmtree(INLINE, ignore_errors=True)
    print(f"done: {len(entries)} tiles")


if __name__ == "__main__":
    main()
