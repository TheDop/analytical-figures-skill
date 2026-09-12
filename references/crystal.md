# crystal.md — conventions for publication-grade crystal-structure figures

## Contents

- Gate vs judgment (where the hard line is)
- What to show
- Orientation (the judgment call, made reproducible)
- Rendering backend (static 3D)
- Displacement ellipsoids (ORTEP) and the ball-and-stick fallback
- View entry points (which function to call)
- Labelling
- Hydrogen bonds & contacts
- Colour — conventional element colours, reconciled with the house style
- Caption requirements (a crystal figure caption MUST carry)
- Calculated PXRD (reuses the `spectra` pxrd domain)
- Sources

The judgment layer for the `crystal` family. The *code* (`crystal_engine`/`_pxrd`/`_view`)
enforces correctness and reproducibility; *this doc* tells you what makes a crystal figure
**right** for the case in front of you. It is guidance, not matplotlib, and deliberately
not a decision tree — the point of doing this in Claude is the case-by-case call.

Conventions below track the journals that referee these figures: **IUCr Acta Cryst**
(Notes for Authors, Sections B/C/E/F), the **ORTEP** convention (Johnson 1965; Cruickshank
1956), and **ACS Crystal Growth & Design** / IUCr publication standards. Sources at the end.

## Gate vs judgment (where the hard line is)

The litmus test: *wrong → false/irreproducible* is a **gate** (the code refuses or fixes);
*wrong → a different-but-valid choice* is **judgment** (yours, per case).

| Enforced by the code (don't re-decide) | Your call (this doc guides it) |
|---|---|
| validation gates (density triple-check, disorder resolution, parse/symmetry, CheckCIF-style alerts) | which object: symmetry-independent unit vs whole molecule vs packing |
| camera **reproducibility** (lives in CONFIG as numbers) | which orientation tells the story (PCA default, or a custom view) |
| locked house style, size-at-print, vector+raster, the read-the-PNG QA loop | what to label; hide C–H or not; which contacts to draw |
| no mirror-flipped enantiomer; units/identity on every drawn quantity | the caption narrative — what the figure argues |

## What to show

- **Default: every symmetry-independent species, fully labelled** (IUCr requires this for a
  structure report). For a cocrystal that means *both* coformers (+ solvent if present).
- **Complete the molecules first.** Render symmetry-completed molecules, not the raw
  asymmetric unit — an AU is often a fragment (a molecule on an inversion centre gives half),
  and two coformers need not be coplanar. (The engine builds whole molecules before orienting.)
- **H atoms / solvent may be omitted** when they don't carry the message (IUCr allows it) —
  but say so in the caption. Hide **C-bound H** by default for a clean view; **keep H on
  donors/acceptors** when the figure is about H-bonding.
- **Extended/coordination structures:** show at least the chemically unique fragment and the
  full coordination environment of any metal.

## Orientation (the judgment call, made reproducible)

PCA face-on is the **automatic default** and what the degeneracy gates reason about. Reach
for a **custom angle** (`view_orientation = "vector"|"custom"`, or a `view_tilt` nudge) when
the science wants a specific viewpoint — it is equally valid because it is still written into
CONFIG (reproducible), which is the only invariant that matters.

- **Match the 2D scheme.** IUCr: the orientation of the species should correspond as closely
  as possible to the chemical drawing. If you have a reference orientation, use `custom`.
- **Single molecule / planar cocrystal → PCA face-on**, rolled so the H-bond network lies
  in-plane. The easy, good case.
- **Layered / π-stacked / channel structures → view *down* the stack or channel axis**
  (`vector = [u v w]`, or `axis_a/b/c`) so the layering/porosity reads. PCA face-on would hide it.
- **Packing diagrams → down a unit-cell axis**, with the **cell edges drawn and a/b/c
  labelled**; use an N×N×N block big enough to show the motif, no bigger.
- **Comparing polymorphs / before-after → put all panels in the *same* orientation** so the
  difference is the only thing that changes (a `custom` angle shared across panels).
- **Globular molecule → there is no "correct" view**; accept the PCA default and say so. The
  degeneracy warning will have fired.
- Whatever you pick, the exact camera goes in the caption (see below).

## Rendering backend (static 3D)

The static structure image is a **raster** (VTK rasterizes — the convention for structure
*images*; vector stays for the line-art figures, PXRD/calibration). Backend by `cfg.view_renderer`:

- **`pyvista` (default; "auto" picks it when installed)** — VTK with a real depth buffer, so
  occlusion is correct and there are **no vertex gaps**. Quality is maxed for a still image:
  smooth-shaded spheres + cylinder bonds, SSAA anti-aliasing, **SSAO** ambient occlusion
  (`cfg.view_ssao`), orthographic camera from the orientation engine, size via
  `cfg.view_resolution`. Heavy dep (`pyvista`/VTK), lazy-installed; needs a GL context (fine on
  a desktop; headless CI needs osmesa/xvfb).
- **`matplotlib` (zero-dependency fallback)** — runs anywhere, but has no depth buffer, so
  atom/bond joints show small white wedges; use only when VTK is unavailable.

## Displacement ellipsoids (ORTEP) and the ball-and-stick fallback

Set `view_style="ellipsoid"` (default `"ball_stick"`); `adp_probability` sets the level.
- **Draw ellipsoids at 50% probability** — the de-facto standard, and the probability level
  **must be stated in the caption** (ORTEP/IUCr). 50% unless you have a reason (very light/heavy
  data sometimes uses 30%). Ellipsoid size/shape *is* the thermal motion: a 100 K structure
  correctly shows tiny near-spheres, a room-temperature one shows large anisotropic ellipsoids —
  so pick a room-T CIF if you want the ADPs to read.
- **H atoms as small fixed-radius spheres**, not ellipsoids (their ADPs are usually riding/isotropic).
- **Isotropic-only CIF → degrade to ball-and-stick** and note it; an ellipsoid plot of riding
  isotropic atoms is meaningless. A clean ellipsoid plot is also a quality/disorder signal — it
  *should* reveal unusual displacements, don't hide them.

## View entry points (which function to call)

- `render(struct, cfg)` — the asymmetric unit's molecule(s), H-bonds dashed.
- `render_unit_cell(struct, cfg)` — one unit cell: whole molecules whose **centroid** lies in
  the cell, + the cell box and a/b/c. A molecule sitting ON a cell face/edge/corner is drawn on
  every face it touches (standard convention — e.g. an ellagic acid on an inversion centre on the
  b-face shows on both y=0 and y=1); the cell still holds Z molecules by count. `cfg.cell_fill`:
  `"molecule"` (default) or `"clip"` (literal cell contents cut at the outer box — for packing
  density / extended structures; whole-molecule reads better for molecular crystals).
- `render_packing(struct, cfg)` — `pack_cells=(Nx,Ny,Nz)` of whole molecules + box + a/b/c, with
  H-bonds dashed. Viewing down a cell axis (`view_orientation="axis_a"|"axis_b"|"axis_c"`) reads best.
- `render_hbond_environment(struct, cfg)` — the asymmetric unit + the symmetry neighbours it
  H-bonds to, intermolecular H-bonds dashed. `cfg.hbond_neighbour` sets the neighbour extent:
  **`"stub"`** (default — contact atom + one bonded shell; the focused structure-report figure),
  `"site"` (contact atom only), or `"whole"` (full neighbour molecules; the packing-shell view).
  Symmetry-generated neighbour atoms are auto-labelled with a roman superscript and the caption's
  symmetry-operation key is logged.

All return `(backend, obj)`; write with `save(rendered, basename, cfg)`. They share one renderer,
so ellipsoids / labels / H-bonds / quality settings apply uniformly.

## Labelling

- **Label the symmetry-independent heteroatoms** (and any atom the text refers to). Do **not**
  label every C/H — clutter buries the point; labels must not be obscured by bonds (IUCr).
- **Symmetry-generated atoms get a superscript** denoting the operation — `O3ⁱ`, `C5ⁱ` (IUCr
  prefers the roman-numeral superscript; `#`/`A` are tolerated) — and the **caption must define
  each symmetry operation** used to generate a labelled atom. *(`render_hbond_environment`
  auto-draws these superscripts on its symmetry-generated neighbour atoms and logs the caption
  key; the other views label without superscripts. The validation table's Block C also carries
  each acceptor's `n_pqr` operator for the caption.)*
- Bake only labels into the image; figure number and the symmetry-operation key live in the
  document caption (house rule: the letter/labels are in the image, the caption is not).
- By default labels sit **at the atom with a white halo** (black text + white outline) — the
  publication-standard way to keep text legible over dark atoms *without* a filled box that
  would hide the molecule. `cfg.view_label_offset` nudges them off the atom; `cfg.view_label_size`
  scales them. (pyvista composites the halo with Pillow since VTK labels can't outline text;
  matplotlib uses `path_effects`.) Disordered (partial-occupancy) sites and symmetry-image
  duplicates are skipped so labels don't pile up. H-bonds are dashed **from the hydrogen**
  (H···A), not the donor heavy atom.

## Hydrogen bonds & contacts

- Draw H-bonds as **dashed D···A lines**; don't write the distances on the plot.
- Put the **geometry in a caption line or SI table**: D–H, H···A, D···A, ∠D–H···A, and the
  **symmetry operator generating A** — reported on the engine's normalized geometry (the table
  is the source of truth; the figure is its view).
- Only assert bonds the engine classifies as H-bonds (≥ the angle floor); sub-floor contacts
  are "geometric contact", not laundered into the figure as H-bonds.

## Colour — conventional element colours, reconciled with the house style

Crystallography readers expect **CPK/element colours** (C grey, O red, N blue, H white/grey,
S yellow, Cl green, F yellow-green, Fe orange…), and recognizability is worth keeping. But the
house style demands colour never be the *only* channel and that grayscale still separates:

- **Labels carry identity**, so colour is redundant by construction — keep element colours.
- Ensure **luminance separation** (dark C, mid O/N, light H) so a grayscale print still reads;
  don't rely on a red/green distinction alone.
- **Distinguishing the two molecules of a cocrystal** is a *second* encoding problem: don't do
  it by element colour (both have C/O/N). Use a redundant channel — one coformer in full element
  colour and the other desaturated/outlined, or label each — and state it in the caption.
  `cfg.color_by_component=True` does exactly this: the **largest molecule** keeps full element
  colour, the rest are **luminance-desaturated** (muted but same brightness, so still legible).
- Never the matplotlib defaults (rainbow/jet); never JPEG (`save_fig` refuses it).

## Caption requirements (a crystal figure caption MUST carry)

The figure documents its own viewpoint and provenance. Include:
1. the **compound / refcode**;
2. **ellipsoid probability level** (e.g. "displacement ellipsoids at 50%") *or* "ball-and-stick"
   if isotropic;
3. **temperature** of the determination;
4. **what is omitted** ("H atoms omitted", "C-bound H omitted", "minor disorder component omitted");
5. the **symmetry-operation key** for every superscripted atom;
6. the **camera** — for a custom/`vector`/tilted view, state it ("viewed down [1 0 1]",
   "viewed down *b*", or the explicit angles) so the view is reproducible from the caption alone;
7. if the model is **reproduced/derived from a published structure**, a source citation
   (IUCr/ACS requirement).

## Calculated PXRD (reuses the `spectra` pxrd domain)

- Plot 2θ on a **normal (low→high) axis**; **state the wavelength** in the caption (use the
  CIF's declared λ; only fall back to Cu Kα₁ 1.5406 Å if undeclared, and say so).
- **Overlay calc vs experimental** for phase ID; align on 2θ; a small calc/exp offset is a
  zero-point or unit-cell/temperature effect, not a phase mismatch — flag, don't hand-shift.
- The (000) beam is removed before broadening (engine handles it). Calc PXRD is a *source*,
  not a new plotting family — house spectra conventions otherwise apply.
- **Cocrystal ID = stack the phases.** `crystal_pxrd.plot_overlay_patterns([(label, cif), …], cfg,
  experimental=…)` waterfalls calc(cocrystal) vs calc(API) vs calc(coformer) (+ the measured
  pattern, drawn at the bottom) so a genuine new cocrystal phase reads as *distinct* from a
  physical mix of its starting materials. Label each trace at the right edge; state λ in the
  caption. `peak_table` / `write_peaks_csv` emit the (2θ, d, hkl, I) list — the hkl is a
  cell-metric **indexing aid** (no structure factors / absences), not a full reflection list.
- **Matching a real lab Cu-Kα scan → `pxrd_realism`.** `crystal_pxrd`'s single-λ, fixed-width
  profile is idealised; a benchtop **Cu** diffractogram shows three things it misses, added by
  `scripts/pxrd_realism.py` on a reflection list `[(2θ, I, hkl)…]`: the **Cu Kα₁/Kα₂ doublet**
  (`kalpha2_doublet` — α2 at ~½ intensity, the split growing with tanθ), **preferred orientation**
  (`march_dollase` — platy/needle habit skews relative intensities; set `cfg.pxrd_po_axis` +
  `cfg.pxrd_march_r`), and **angle-dependent broadening** (`caglioti_fwhm` + `pseudo_voigt`, tuned
  by `cfg.pxrd_caglioti` + `cfg.pxrd_lorentz_fraction`). **`crystal_pxrd.realistic_pattern(struct,
  cfg)`** is the turn-key: it pulls the real Dans powder reflection list via
  `crystal_pxrd.reflection_list` (structure factors with multiplicity + Lorentz-polarization
  already applied — NOT the post-broadened peak list, which would double-broaden) and composes the
  realism with `simulate_pattern`. Use it when overlaying calc vs a measured Cu-Kα scan so peak
  *shapes/positions* line up; for pure phase-ID (positions only) the plain `calc_pattern` is enough.
  **Coherence gotcha:** the Kα2 default is *Cu* — to simulate a Cu scan from a Mo-refined CIF set
  `cfg.pxrd_wavelength=1.540598` (so α1 is Cu too). Validated on real CIFs (aspirin, lactose)
  in `skill_validation/crystal/test_realism_cif.py`.

## Sources

- IUCr, *Notes for Authors*, Acta Cryst. Sections B/C/E/F — labelled displacement-ellipsoid
  diagram of each symmetry-independent species; 50% probability stated in the caption;
  symmetry-related atoms superscripted with the operation defined in the caption; labels
  unobscured. journals.iucr.org/c/services/notesforauthors.html (and /b, /e, /f).
- ORTEP: C. K. Johnson (1965), *ORTEP*, ORNL-3794; D. W. J. Cruickshank (1956) — ellipsoid =
  probability contour; probability level always quoted (typically 50%).
- IUCr, *Publication standards for crystal structures* (2011); ACS *Crystal Growth & Design*
  Author Guidelines — fully labelled ORTEP-type figure; orientation to match the chemical
  scheme; CheckCIF validation; cite reproduced structures.
- Bondi (1964) vdW radii; Allen & Bruno (2010) X–H normalization; Jeffrey / Steiner H-bond
  ranges — the pinned numeric criteria are the constants block at the top of `scripts/crystal_engine.py`.
