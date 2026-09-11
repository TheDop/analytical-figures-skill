#!/usr/bin/env python3
"""
build_gallery.py — the README gallery, the GitHub Pages gallery and the social-preview card,
all from ONE manifest, using only figures the shipped code produces on public data, public
CIFs, or the synthetic validation set. Nothing here is drawn by hand.

    python docs/build_gallery.py            # re-run the producers (~1-2 min), then compose
    python docs/build_gallery.py --no-run   # compose from the producers' existing out/ figures

Outputs (all under docs/):
    gallery/tiles/<id>.png    square thumbnails (the README table + the Pages grid)
    gallery/full/<id>.png     full-resolution copies (the click-through targets)
    gallery/social.png        2:1 mosaic with the title card (upload as the GitHub social preview)
    index.html                the Pages gallery: hover a tile for its placard (plain HTML + CSS)
    ../README.md              the tile table, written between <!-- gallery:start/end --> markers

Rebuild only when a producer or a tile changes — the PNGs are committed, and every rebuild adds
to the repository history.
"""
import os
import sys
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

# ---------------------------------------------------------------- producers (shipped scripts)
PRODUCERS = [
    f"{FT}/validate_ftir_integration.py",
    f"{FT}/groundtruth_integration_test.py",
    f"{CR}/crystal_demo.py",
    f"{CR}/test_realism_cif.py",
    f"{CR}/test_ellipsoid.py",
    f"{CR}/test_pxrd.py",
    f"{CR}/test_pyvista.py",          # optional: skipped when pyvista is absent
]

# ---------------------------------------------------------------- the manifest
# id, source PNG (relative to the skill root; INLINE:<name> for tiles rendered below),
# caption, data provenance, producing script
MANIFEST = [
    ("calib_pectin", f"{FT}/out/pectin_calibration_published.png",
     "Calibration with confidence and prediction bands and the mandatory residual panel",
     "Pectin degree of methyl-esterification by FTIR — Wang 2023, Mendeley Data 10.17632/gkwbp3wc49.1 (CC BY 4.0)",
     f"{FT}/validate_ftir_integration.py"),
    ("waterfall", "INLINE:waterfall",
     "Replicate waterfall: reps overlaid per level, levels offset, right-margin keys instead of a legend",
     "Synthetic binary mixture, 6 levels × 3 reps with multiplicative loading variation (skill_validation/chemometrics/_data.py)",
     "docs/build_gallery.py → scripts/spectra.plot_waterfall"),
    ("pls_diag", "INLINE:pls_diagnostics",
     "Leakage-safe PLS diagnostics: RMSECV per preprocessing with the parsimony pick, coefficients, scores, loadings",
     "Synthetic binary mixture (skill_validation/chemometrics/_data.py), leave-one-level-out CV",
     "docs/build_gallery.py → scripts/chemometrics.diagnostics_figure"),
    ("pxrd_overlay", f"{CR}/out/crystal_demo/pxrd_overlay.png",
     "Stacked calculated PXRD patterns for phase discrimination",
     "Aspirin COD 7247819, urea COD 1008785, lactose COD 2206486 — Crystallography Open Database",
     f"{CR}/crystal_demo.py"),
    ("ortep", f"{CR}/out/crystal_demo/structure_ellipsoid.png",
     "ORTEP displacement ellipsoids at 50 % probability, deterministic PCA camera",
     "Aspirin, COD 7247819 (296 K)",
     f"{CR}/crystal_demo.py"),
    ("ddsimca", "INLINE:ddsimca",
     "DD-SIMCA acceptance plot: regular / extreme / outlier, with the chi-squared boundary",
     "Synthetic binary mixture plus one under-loaded and one baseline-tilted replicate injected",
     "docs/build_gallery.py → scripts/chemometrics.diagnostics + plot_influence"),
    ("integration", "INLINE:integration",
     "Band area on algorithmic anchors: the flanking minima found once and locked for the whole batch",
     "Synthetic binary mixture (skill_validation/chemometrics/_data.py)",
     "docs/build_gallery.py → scripts/spectra.find_anchors"),
    ("calib_bmc", f"{FT}/out/bmc_dox_transmittance.png",
     "Doxorubicin carbonyl-area calibration; LOD/LOQ from the residual SD, and the sigma is named",
     "Bansal, Singh & Kaur 2021, BMC Chemistry 15:27, 10.1186/s13065-021-00752-3 (published area table)",
     f"{FT}/validate_ftir_integration.py"),
    ("groundtruth", f"{FT}/out/groundtruth_shared_vs_perwindow.png",
     "Shared vs per-window baseline against a closed-form answer: the per-window pedestal bias, measured",
     "Analytic Gaussian bands with known areas, checked against scipy.integrate.quad",
     f"{FT}/groundtruth_integration_test.py"),
    ("pxrd_calc", f"{CR}/out/crystal_demo/pxrd.png",
     "Calculated PXRD from a CIF at Cu Kα, top peak cross-checked against pymatgen",
     "Aspirin, COD 7247819",
     f"{CR}/crystal_demo.py"),
    ("pxrd_realism", f"{CR}/out/_realism_aspirin.png",
     "Lab-realistic pattern: Cu Kα₁/Kα₂ doublet and Caglioti broadening on the Dans reflection list",
     "Aspirin, COD 7247819",
     f"{CR}/test_realism_cif.py"),
    ("packing", f"{CR}/out/crystal_demo/packing.png",
     "Packing diagram down a with hydrogen bonds and the cell box",
     "Aspirin, COD 7247819",
     f"{CR}/crystal_demo.py"),
    ("hbond_env", f"{CR}/out/crystal_demo/hbond_environment.png",
     "Hydrogen-bond environment of the asymmetric unit, symmetry mates superscripted",
     "Aspirin, COD 7247819",
     f"{CR}/crystal_demo.py"),
    ("unit_cell", f"{CR}/out/crystal_demo/unit_cell.png",
     "Unit-cell contents with the cell box and a/b/c",
     "Aspirin, COD 7247819",
     f"{CR}/crystal_demo.py"),
    ("flufenamic", f"{CR}/out/_ellipsoid_disorder-groups_flufenamic.png",
     "Z′ = 3 with CF₃ disorder groups honoured: no bonds between mutually exclusive sites",
     "Flufenamic acid, COD 4118081",
     f"{CR}/test_ellipsoid.py"),
    ("pxrd_urea", f"{CR}/out/_pxrd_specialpos_urea.png",
     "Calculated PXRD of a molecule on a special position: multiplicity handled in the expansion",
     "Urea, COD 1008785",
     f"{CR}/test_pxrd.py"),
    ("ballstick", f"{CR}/out/crystal_demo/structure.png",
     "Ball-and-stick, PCA face-on, C–H hidden, heteroatoms labelled",
     "Aspirin, COD 7247819",
     f"{CR}/crystal_demo.py"),
    ("pyvista", f"{CR}/out/_pyvista_aspirin.png",
     "PyVista/VTK still: real depth buffer, ambient occlusion, orthographic camera from the same orientation engine",
     "Aspirin, COD 7247819",
     f"{CR}/test_pyvista.py"),
]

TITLE = "analytical-figures"
TAGLINE = ["Publication-grade figures and the gated", "analysis behind them, for analytical chemistry."]
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
def render_inline():
    """Four tiles drawn with the skill's own modules on the synthetic validation set."""
    import numpy as np
    sys.path.insert(0, ROOT)
    from scripts.config import Config
    from scripts import style, spectra, chemometrics

    spec = importlib.util.spec_from_file_location("_data", os.path.join(ROOT, V, "chemometrics", "_data.py"))
    data = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(data)
    x, X, y, groups = data.synth(n_levels=6, reps=3, n_wave=400, seed=0, noise=0.002)
    os.makedirs(INLINE, exist_ok=True)
    cfg = Config(output_dir=INLINE, formats=("png",), journal="nature", column="single",
                 domain="ftir", strict=False, description="gallery tile — synthetic validation set")
    out = {}

    # (a) replicate waterfall with edge labels
    levels = sorted(set(groups))
    wf = {f"{y[groups == g][0]:.0f} % w/w": [(x, X[i]) for i in np.where(groups == g)[0]] for g in levels}
    fig, ax = spectra.plot_waterfall(wf, cfg)
    style.save_fig(fig, os.path.join(INLINE, "waterfall"), cfg)
    out["waterfall"] = os.path.join(INLINE, "waterfall.png")

    # (b) band area on locked algorithmic anchors
    i = int(np.where(groups == levels[3])[0][0])
    lo, hi = spectra.find_anchors(x, X[i], 1748.0, gap=12.0, maxhw=60.0)
    xs, ys = x[::-1], X[i][::-1]                                  # ascending for interp
    ylo, yhi = np.interp(lo, xs, ys), np.interp(hi, xs, ys)
    style.apply_style(cfg)
    fig, ax = style.figure(cfg)
    ax.plot(x, X[i], lw=0.9)
    m = (xs >= lo) & (xs <= hi)
    base = np.interp(xs[m], [lo, hi], [ylo, yhi])
    ax.fill_between(xs[m], base, ys[m], alpha=0.30, lw=0)
    ax.plot([lo, hi], [ylo, yhi], color="0.25", lw=0.8, ls="--")
    ax.plot([lo, hi], [ylo, yhi], "o", color="0.25", ms=3)
    spectra._apply_axis(ax, cfg)
    ax.set_xlim(1800, 1560)
    ax.set_title(f"ester C=O area, anchors locked at {lo:.0f} / {hi:.0f} cm⁻¹", fontsize="small")
    style.save_fig(fig, os.path.join(INLINE, "integration"), cfg)
    out["integration"] = os.path.join(INLINE, "integration.png")

    # (c) PLS diagnostics 2x2 (leave-one-LEVEL-out) — a 2x2 needs the double-column width;
    #     settle the constrained layout before export so no label is clipped
    from dataclasses import replace
    cfg2 = replace(cfg, column="double")
    fig, axes, info = chemometrics.diagnostics_figure(X, y, x, cfg2, groups=groups, basename="pls_diagnostics")
    style.finalize_figure(fig, wspace=0.12, hspace=0.22)
    style.save_fig(fig, os.path.join(INLINE, "pls_diagnostics"), cfg2)
    out["pls_diagnostics"] = os.path.join(INLINE, "pls_diagnostics.png")

    # (d) DD-SIMCA acceptance plot with two injected bad replicates
    X2 = X.copy()
    X2[4] = X2[4] * 0.45                                          # under-loaded
    X2[13] = X2[13] + 0.03 * np.linspace(0, 1, x.size)            # baseline tilt
    diag = chemometrics.diagnostics(X2, cfg, n_components=2, pre="none")
    style.apply_style(cfg)
    fig, ax = style.figure(cfg)
    chemometrics.plot_influence(ax, diag, cfg, labels=[f"s{k}" for k in range(X2.shape[0])])
    style.save_fig(fig, os.path.join(INLINE, "ddsimca"), cfg)
    out["ddsimca"] = os.path.join(INLINE, "ddsimca.png")
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
        # full-resolution copy, long side capped so the repository stays small
        cap = 1400
        if max(im.size) > cap:
            s = cap / max(im.size)
            full = im.convert("RGB").resize((round(im.width * s), round(im.height * s)), Image.LANCZOS)
        else:
            full = im.convert("RGB")
        full.save(os.path.join(FULL, e["id"] + ".png"), optimize=True)
        t = _contain(im, 600)
        t.save(os.path.join(TILES, e["id"] + ".png"), optimize=True)
        tiles[e["id"]] = t

    # social card: 6 x 3 grid, title card over the first two slots
    gap, tile, cols, rows = 12, 384, 6, 3
    W, H = cols * tile + (cols + 1) * gap, rows * tile + (rows + 1) * gap
    card = Image.new("RGB", (W, H), "#e6e9ee")
    slots = [(r, c) for r in range(rows) for c in range(cols)][2:]            # slot 0-1 = title card
    for (r, c), e in zip(slots, entries):
        x0, y0 = gap + c * (tile + gap), gap + r * (tile + gap)
        card.paste(tiles[e["id"]].resize((tile, tile), Image.LANCZOS), (x0, y0))
    x0, y0 = gap, gap
    box = Image.new("RGB", (2 * tile + gap, tile), "white")
    d = ImageDraw.Draw(box)
    d.text((34, 60), TITLE, font=_font("DejaVuSans-Bold.ttf", 66), fill="#1b1f24")
    for k, line in enumerate(TAGLINE):
        d.text((36, 170 + k * 38), line, font=_font("DejaVuSans.ttf", 27), fill="#3b4048")
    d.text((36, 312), FOOT, font=_font("DejaVuSans.ttf", 19), fill="#6b7280")
    d.rectangle([0, tile - 6, 2 * tile + gap, tile], fill="#0b5fff")
    card.paste(box, (x0, y0))
    card.save(os.path.join(GAL, "social.png"), optimize=True)
    print(f"  social.png {W}x{H}, {len(entries)} tiles")


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
             f"produced by <code>docs/build_gallery.py</code> from public data, public CIFs or the synthetic validation "
             f"set — nothing is drawn by hand. Interactive version with placards: <a href=\"{PAGES_URL}\">{PAGES_URL}</a></sub>")
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
<meta property="og:description" content="{' '.join(TAGLINE)}">
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
  <p class="tag">{' '.join(TAGLINE)} A <a href="{REPO_URL}">Claude Code skill</a>.</p>
  <p class="hint">Hover a tile for what it is and where the data came from; click for full size.</p>
</header>
<main>
{chr(10).join(cards)}
</main>
<footer>
  Every tile is produced by <code>docs/build_gallery.py</code> from public data, public CIFs or the synthetic
  validation set — nothing is drawn by hand. <a href="{REPO_URL}">Repository</a> · MIT licence.
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
