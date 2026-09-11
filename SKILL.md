---
name: analytical-figures
description: >-
  Publication-grade figures and the gated analysis behind them for analytical
  chemistry / spectroscopy work. Use for ANY of: plotting FTIR / IR / ATR or
  PXRD / XRD spectra (overlays, waterfalls, baseline correction, band
  integration, %RSD/precision); building calibration curves (linear regression,
  confidence/prediction bands, residuals, LOD/LOQ, recovery, ICH Q2); generic
  line/bar charts for lab data; or turning raw instrument files (.spc, xy/csv)
  into journal-ready PDF+PNG. Trigger on:
  spectrum, spectra, FTIR, IR, ATR, PXRD, XRD, diffractogram, waterfall,
  baseline, band area, integration, RSD, calibration curve, linearity, LOD,
  LOQ, residuals, recovery, publication figure, journal figure; crystal structure,
  CIF, .cif, space group, hydrogen bond, displacement ellipsoid / ORTEP, packing
  diagram, calculated PXRD from a CIF, structure validation table.
---

# analytical-figures

A figure + analysis skill for spectroscopy and analytical QC. It assumes you
already know matplotlib. Its job is the two things that actually go wrong:
**the numbers are subtly wrong before they ever reach a plot**, and **the figure
looks "fine" but is broken in a way a quick glance misses**. So it forces
verification and it locks a house style. Nothing here is a plotting tutorial.

If a request is just "how do I do X in matplotlib", answer it directly; this
skill is for producing *deliverables* from analytical data.

## Output contract: a self-contained script, not a figure

The deliverable this skill produces is **one runnable, self-contained Python
script** — a CONFIG block plus the analysis body — that the user runs to generate
the figure. Do **not** render the figure in-session and hand over a PNG; hand over
the *script that produces it*, so it is reproducible, editable, and attachable to
a report. The modules in `scripts/` are the maintained source of truth; you write
the analysis against them, then `bundle.py` amalgamates the modules it uses + your
body into a single file that depends only on the scientific stack (numpy,
matplotlib, optionally scipy/scienceplots), not on this skill being installed:

```
python bundle.py my_analysis.py -o standalone.py    # hand over standalone.py
```

**Provenance is stamped, not narrated.** `bundle.py` writes a provenance header into the
deliverable — skill version, the Python + scientific-stack versions it was built against,
the source body, and a one-line plain-language description (`--describe "…"`, or a
`DESCRIPTION="…"` constant in the body) — so the "exact code + environment that produced
it" travels with the figure and is recoverable months later. `style.save_fig` embeds the
same version + `cfg.description` into each figure FILE's metadata (PDF/SVG/PNG). At handoff
`bundle.py` also runs a **critic pass** (below).

**Portability (NumPy 2.0+):** `np.trapz` was **removed** in NumPy 2.0 and renamed
**`np.trapezoid`**. A bundled standalone must run on current NumPy, so never write
`np.trapz` in the analysis body — use `np.trapezoid`, or the family shim
`spectra._trapz` (`= getattr(np, "trapezoid", getattr(np, "trapz", None))`, works on
both), and prefer `spectra.band_metric` for band integration.

## Form is the skill's job; judgment is the model's

The skill enforces **form** — the style, the gates, the script scaffold, the label
*typography*. The model supplies **judgment** — which bands to label and what each
label says, which chart answers the question, the narrative. Never bake a content
decision into the skill. (Concretely: `label_peaks` styles labels you pass in; it
does not detect peaks or write text. Guardrails that enforce *correctness* —
refusing low-n mean bars, mandating residual panels — are form, not content, and
stay.)

**Prefer captions to floating text.** Keep the plot area clean: avoid free-floating
annotations, value labels, target-line labels, and stat strings inside the axes —
put those details in the caption. Identify traces by position (waterfall stack
order; an x-axis that already encodes the variable) and state the order in the
caption, rather than by in-axes labels or a colour-only legend. Keep *silent*
visual pointers (shaded bands, position guides, a dashed target line) and explain
them in the caption. When you strip text, emit the details (band assignments, stack
order, numbers) for the caption — don't lose them.

## The one rule that holds it together: a single CONFIG block

Everything tunable lives in the `Config` dataclass (`scripts/config.py`) and is
set ONCE at the top of an analysis. Copy `templates/analysis_template.py` and
edit only that block. If you want to tune something that isn't in `Config`, add
it to `Config` first — do not bury knobs in script bodies.

## Workflow

0. **Choose the figure** — before touching code, settle WHAT to plot and WHY with
   `references/figure_selection.md` (data shape × argument × scale → figure type, with the
   citable rule behind each call). It standardises the *form* decision; it doesn't decide your
   science. A capable model can answer this ad hoc — the doc's value is a *consistent* answer.
1. **Configure** — copy the template, fill the CONFIG block (domain, journal,
   windows, conf level…).
2. **Ingest + gate** the raw numbers with `verify.check_ingest` / `check_trace`
   BEFORE plotting (finite, monotonic axis, sane point count). A decode error
   caught here saves a wrong figure later.
3. **Process identically across a batch.** Baseline/normalise/integrate every
   rep with the same cfg, record how, and assert it with `verify.same_processing`
   so rep-to-rep differences are real, not artefacts.
4. **Stats** via `verify.summarize` — returns mean, SD, SEM, t-based CI *and* a
   caption naming the error type and n. Never present an error bar without saying
   which statistic it is and what n is.
5. **Plot** with the family module (`spectra` / `calibration` / `charts`) using
   `style.figure`. For a **composite** (multi-panel) figure, after drawing all
   panels call `style.finalize_figure(fig)` **then** `style.add_panel_labels(fig,
   cfg)` — in that order (settle the layout, *then* anchor the letters).
6. **Visual QA loop** (`verify.render_preview` → `verify.audit_layout` [glyphs,
   clipping, data escaping the axes, x+y tick overlap, one a/b/c letter per panel] →
   open the PNG with the **Read tool** and walk `verify.READ_IMAGE_CHECKLIST` → fix →
   re-render). Export the vector master only after a clean pass.

**Normalise stacked traces by subtracting the window MINIMUM, not the median.** With median
subtraction a sloping baseline puts the window minimum below zero and the lowest trace runs off
the bottom of the panel — a fault no text-clipping check can see, because no text is clipped.
Subtracting the minimum maps each trace to exactly [0, 1]; size the panel from
`(n_traces - 1) * step + ~1.08`. `audit_layout` now catches the residual cases: it judges **y
only**, and only for points already on-screen in x, so a deliberate zoom, `axvline`/`axhline`
guides and an explicit `clip_on=False` are all left alone.
7. **Export** with `style.save_fig` (vector + raster proof, sized at final dims).

## Hard principles (these are gates, not preferences)

- **Size at final print dimensions.** Set `figsize`/journal+column once; never
  rescale in Word/LaTeX (pt fonts shrink with it).
- **Vector first + a raster proof.** PDF/SVG for data figures; PNG/TIFF only for
  images. **Never JPEG** (`save_fig` refuses it).
- **Colour-blind safe + redundant encoding.** Okabe-Ito/colorblind palette plus a
  second channel so grayscale still separates — line-style/marker for generic charts,
  but **spectra stay solid (continuous) by default** and take their redundancy from
  waterfall *position* (dashes are opt-in only; see `references/spectra.md`).
- **Fonts ≥ 6 pt at final size** (`style` enforces the floor).
- **Errors are always named** — SD vs SEM vs CI is a √n factor; state the
  statistic, n, and (for CI) the t-multiplier.
- **Numbers carry units** — every written quantity states its unit. Axis tick
  numbers get theirs from the axis label, so each axis label must name a unit (or
  mark it `dimensionless` / `a.u.` / `ratio` / `%` when it genuinely has none —
  e.g. absorbance, R²). The same applies to residual axes (same unit as the
  response) and to any value quoted in a caption. `audit_layout` WARNs on an axis
  whose label carries no unit cue.
- **On-figure numbers must trace to the code.** A figure-of-merit printed on the figure
  or written into a caption (R², RMSEP, LOD, recovery, %RSD, p, n, ±error…) is
  **interpolated from the computed value** (an f-string), never hand-typed — a pasted
  literal silently desyncs the moment you refit ("R2 = 0.98" in the text, `0.91` in the
  code). `verify.check_number_provenance` (the **Tier-3 critic**, run automatically at
  `bundle.py` handoff) flags a hard-typed statistic; a fixed physical identifier (a band
  centre like 1748 cm⁻¹) is a *constant*, not a result, and is exempt. Adapted from the
  actor-critic reviewer in Anthropic's Claude Science.
- **Calibration: a residual panel is mandatory** — R² alone may not stand as the
  linearity claim; refuse extrapolation outside the calibrated range.
- **Identical processing across a batch** — rep-to-rep differences must not be an
  artefact of inconsistent baselines/windows.
- **No positional jitter on a quantitative axis.** Replicates that share an x-value
  are drawn AT that value — the *other* axis separates them. Never nudge points along a
  measured/numeric axis to un-stack them: it fabricates coordinates that don't exist (a
  40 % w/w point drawn at 42 reads as 42 %). If markers overlap, separate them on the
  perpendicular axis, use open/transparent markers or small multiples, or make the axis
  deliberately CATEGORICAL (x carries no numeric meaning) — but a value axis shows true
  values only.

## Domain profiles (spectra)

| | FTIR / IR / ATR | PXRD / XRD |
|---|---|---|
| x axis | wavenumber cm⁻¹, **inverted** | 2θ degrees, **normal** |
| y | absorbance | intensity (counts / normalised) |
| baseline | arPLS / rubberband, then integrate bands | background subtract; watch amorphous halo |
| idiom | overlaid reps + waterfall; band area / %RSD | stacked patterns + reference stick lines; phase ID |

Set `cfg.domain` and the modules flip conventions automatically. Details and
recipes in `references/spectra.md`.

## Chemometrics (multivariate calibration: PLS / PCR / PCA)

When **no single band is selective** for the analyte (overlap, matrix effects, or
scatter/sample-loading variation that swings every band together), regress the WHOLE
spectrum instead of one integrated band. `scripts/chemometrics.py` is the family for
this — its value is the same as the rest of the skill: it gates the things that go
*silently* wrong in chemometrics, which are not the same as in univariate work. Three
gates, encoded so you can't bypass them by accident: **(1) CV leakage** — a fresh
`Preprocessor` is re-fit on the training rows INSIDE every fold (MSC reference,
mean-centring, autoscaling all leak if fit on the full set; SNV/derivatives don't, but
the path is uniform); **(2) component over-fit** — `choose_n_components` takes the
parsimony pick (fewest components within `cfg.lv_parsimony_tol` of the global-min
RMSECV), and the RMSECV-vs-components curve is shown not asserted; **(3) the optimistic-CV
trap** — leave-one-out over replicate spectra answers *"predict a SEEN level"*, so every
result NAMES its question and `cross_validate` WARNs you to pass `groups=level-labels`
for leave-one-LEVEL-out when predicting a new level. Report **RMSECV / R²(CV) / bias**,
never calibration R² as accuracy. `diagnostics_figure` is the turnkey 2×2 (RMSECV·
coefficients·scores·loadings) on a constrained-layout figure with the colourbar and
panel letters handled; `methods_text` writes the paste-ready paragraph. **scikit-learn
is a LAZY dependency** (imported only when a model is fit, like the crystal heavy deps) —
an FTIR/univariate job never pulls it, and a bundled deliverable keeps the import inside
the functions. A separate **outlier/diagnostic gate** (`diagnostics` / `plot_influence`) does
PCA-based per-sample QC — Hotelling T², Q/SPE, leverage and the **DD-SIMCA** SD–OD acceptance
plot (flag "bad load" vs "extreme but valid") — using **numpy+scipy only**, so it runs as
raw-spectrum QC *before* any calibration. **EMSC** (`emsc`) adds a jointly-fitted-baseline scatter
correction with an optional pure-component **interferent** term (model a known matrix component away);
a **validation trio** (`permutation_test`, `corrected_paired_t` = Nadeau–Bengio, `rpd`/`rpiq`)
gives defensible figures-of-merit for a small calibration. Conventions, the preprocessing menu, and
the coefficient-interpretability caveat (raw model to SEE bands, SNV model to REPORT error) are in
`references/chemometrics.md`.

## Crystal structures (CIF)

A third family, for **crystallography**: a CIF → (1) a **validation table** (density
triple-check, geometry, H-bonds, CheckCIF-style alerts), (2) **calculated PXRD**, and (3)
a **deterministic 3D structure view**. `crystal_engine` validates *before* anything is
drawn — a figure is a view of the validated computation, never independent. Heavy deps
(`gemmi`, `Dans_Diffraction`, optional `pymatgen`/`py3Dmol`) install lazily at run time, so
an FTIR/calibration job never pulls them. Presentation conventions (what to show,
orientation, ellipsoids, labels, captions — grounded in IUCr/ORTEP) live in
`references/crystal.md`. The hard floor is **correctness + reproducibility** (the camera is
numbers in CONFIG); *which* view and what to label is the model's case-by-case call —
PCA face-on is the auto default, custom/`vector`/axis cameras are first-class. The 3D view
renders via **PyVista/VTK** (real depth buffer → gap-free, SSAO+SSAA, offscreen) with a
matplotlib fallback; structure images are high-res raster, PXRD/table stay vector. Views:
`render` (molecule), `render_unit_cell`, `render_packing`, `render_hbond_environment`; ORTEP
displacement ellipsoids via `view_style='ellipsoid'` (`adp_probability`), packing via `pack_cells`.
`render_hbond_environment` neighbour extent = `hbond_neighbour` (`stub` default / `site` / `whole`,
symmetry mates auto-superscripted); unit-cell/packing fill = `cell_fill` (`molecule` / `clip`).
Cocrystal ID — **start with `references/cocrystal_id.md`** (the decision layer: cocrystal vs salt
vs physical mixture vs polymorph, and which evidence settles it — PXRD new-phase, the FTIR
proton-transfer/ΔpKa salt discriminator, the lab-capability gate — routing into the tools below).
`crystal_pxrd.plot_overlay` waterfalls calc(cocrystal) vs starting materials (+ experimental) for
phase discrimination; `peak_table`/`write_peaks_csv` emit a (2θ, d, hkl, I) list;
`cfg.color_by_component` mutes all but the largest molecule to distinguish API vs coformer.

## Composite / multi-panel figures

Most report figures are composites — an a/b panel pair, a spectrum over its
residual strip, three calibrations on shared axes. Build them as ONE figure, never
by pasting separate PNGs together (and panels that invite comparison must literally
`sharex=`/`sharey=` so the comparison is honest):

1. **Grid:** `fig, axes = style.figure(cfg, nrows, ncols, height_ratios=…,
   width_ratios=…)` — real subplots, sized once at final print dimensions. The
   default height scales with `nrows/ncols` so a grid isn't squashed into one
   row's height; `width_ratios` sizes columns (e.g. make one panel a bit narrower so
   its neighbour's right-margin edge labels have room).
2. **Draw** each panel with the family modules. A waterfall panel keys its traces
   with `spectra.edge_labels(ax, items)` — right-margin labels aligned per trace
   (grayscale-safe, beats a colour-only legend box). Prefer keying the RIGHTMOST panel
   so the labels fall in the figure's right margin (reserve it with `rect`, next step);
   if you must label an interior panel, its labels land in the gutter — widen that with
   `wspace`. When panels share a colour map + stack order, ONE such key serves them all.
3. **Finalise, THEN label — in that order.** `style.finalize_figure(fig,
   rect=…, wspace=…, hspace=…)` settles the margins (constrained layout, tight
   fallback) so nothing clips or overlaps. The layout engine can't see edge labels
   (clip_on=False), so `rect=(left,bottom,width,height)` RESERVES an outer margin for
   right-edge labels, and `wspace`/`hspace` widen inter-panel gutters for interior
   ones. Only then `style.add_panel_labels(fig, cfg)` places the a/b/c letters —
   anchored at each panel's top-left and pushed by ONE shared points offset, so they
   line up whatever the tick-label widths (pass `x_offset_pt="auto"` to derive that one
   offset from the widest y-furniture). Labelling *before* the layout settles misplaces them.
4. `audit_layout` confirms one letter per panel AND that no text spills from one
   panel into another (widen `wspace` or shorten the label if it does); the
   Read-the-PNG pass confirms the letters line up.

**House rule — keep gutters TIGHT (~0.1, never 0.2–0.4).** For a compact print figure
the inter-panel `wspace`/`hspace` belong **around 0.1** (as the chemometrics panel uses:
`wspace=0.10, hspace=0.12`), NOT the 0.2–0.4 that spreads panels apart and shrinks the
data. Start tight; only *widen one specific gutter* when an interior panel's right-margin
edge labels would otherwise spill — never inflate every gutter to make room for one.

**If panels still look far apart at `wspace≈0.1`, the gutter is NOT the cause — TEXT is.**
Constrained layout sizes each gutter from the neighbouring axes' *tight* bboxes, which
include y-labels, tick labels and any `clip_on=False` text. So a single long string silently
buys itself a wide gutter, and raising `wspace` to "fix" the look makes it worse. Diagnose by
printing `ax.get_position()` for every panel — if the axes occupy well under ~85 % of the
figure width, the space is going to text, not spacing. Two habitual offenders and their fixes:
- **A long y-axis label** (`"ΔRwp on deleting theophylline / pp"`) sits in the *inter-panel*
  gutter and pushes the columns apart. **Break it over two lines** (`\n`) — same information,
  half the width.
- **Long categorical tick labels** (`"1:1 mix (weighed)"`) eat the outer margin. Shorten the
  tick text and put the qualifier in the caption.

**Do not pass `rect` to reserve room for right-margin `edge_labels` under constrained layout** —
it already counts them in the axes' tight bbox, so `rect` pays for the same text twice and
squeezes every panel inwards. Reserve with `rect` only when the layout engine genuinely cannot
see the artist (a `fig.text` in figure coordinates), and confirm by reading the SAVED figure.

**The letter is in the IMAGE; the figure number and caption are not.** Bake only
the panel letter into the figure. "Figure 3", the caption, and cross-references
("see Figure 3") live in the document (Word/LaTeX caption feature), so the image
stays reusable and renumbering can't desync. Never render "Figure 3" into the plot.

### Sectioned composites & repeated-panel grids (the page-figure idioms)

A full-page figure that stacks logically distinct *sections* (e.g. "what we measure"
over "candidate calibrations" over "figures of merit") goes beyond one flat grid:

- **Sections via subfigures.** `fig = plt.figure(layout="constrained")`, then
  `subs = fig.subfigures(n, 1, height_ratios=…)`; give each section its own internal grid
  and, if it needs one, its own legend (`sf.legend(loc="outside lower center")`). Each
  section's spacing stays independent.
- **Guarantee an outer margin.** Constrained layout can pull axes flush to the edge.
  `fig.get_layout_engine().set(w_pad=…, h_pad=…)` reserves margin on ALL sides. **Judge the
  margin on the SAVED figure, not `render_preview`** — the preview tight-crops, so it always
  looks margin-less.
- **Repeated panel-units → a CSS-grid, NOT a tall grid with spacer rows.** For N identical
  units (a calibration over its residual strip, one per candidate), build an OUTER gridspec
  of the *cells* with ONE uniform `hspace`/`wspace`, and an INNER `subgridspec` per cell for
  the tightly-coupled sub-axes:
  ```
  outer = sf.add_gridspec(2, 2, hspace=0.22, wspace=0.34)        # one uniform gap between cells
  for i, k in enumerate(items):
      r, c = divmod(i, 2)
      cell = outer[r, c].subgridspec(2, 1, height_ratios=[3.4, 1], hspace=0.06)  # cal+resid, tight
      a0 = sf.add_subplot(cell[0]); a1 = sf.add_subplot(cell[1], sharex=a0)
  ```
  A single tall grid with spacer rows makes the gaps COMPOUND (spacer height + two hspaces)
  into an uneven "big gap in the middle"; the nested grid gives one even gap, like CSS `gap`.
- **Per-chart letters in reading order, OUTSIDE the corner.** Sectioning/nesting breaks
  `add_panel_labels`' automatic enumeration, so pass the charts explicitly, in reading order:
  `fig.canvas.draw(); style.add_panel_labels(fig, cfg, axes=[…], x_offset_pt="auto")` — the
  letters then sit outside the top-left corner, aligned, exactly as on a flat grid. A coupled
  cal+residual pair is ONE chart → letter only its main axes (`audit_layout` treats a residual
  strip as part of its chart). `style.panel_letter(ax, s)` puts a letter INSIDE the corner and is
  only for the rare panel with no margin to spare; never mix the two conventions in one figure.
- **Direction-of-good in comparison titles.** Any FoM/comparison panel where taller/larger
  isn't self-evidently "better" states the good direction IN its title — `"Precision\n(lower
  = better)"`, `"Accuracy\n(100 = ideal)"`, `"Linearity\n(higher = better)"` — so a bar chart
  can't be misread.
- **The title names the species, not just the band.** `"Aspirin ester C=O (1748 cm⁻¹)"`, not
  `"Ester C=O"` — so a metric named elsewhere ("ester area") resolves to a species. Band
  identifiers in titles/tick-labels carry their unit (`1066 cm⁻¹`, never bare `1066`).

## Files

- `scripts/config.py` — the `Config` dataclass (the one tunable block).
- `scripts/style.py` — `apply_style`, `figure` (takes `width_ratios`; grid-aware
  default height), `save_fig`, plus composite helpers `finalize_figure` (takes
  `wspace`/`hspace`/`w_pad`/`h_pad` to widen gutters, `rect` to reserve an outer margin
  for right-edge labels) + `add_panel_labels` (aligned a/b/c letters for a flat grid;
  `x_offset_pt="auto"` derives one shared offset clearing the y-furniture) + `panel_letter`
  (ONE inside-corner letter, for hand-lettering sectioned/nested composites in reading order).
  Plus `check_figure_width` (verify a figure is at the journal's physical width — no silent
  rescaling) and `export_mplstyle` (write a standalone `.mplstyle`; presets in `assets/styles/`).
  `save_fig` stamps provenance (skill version + `cfg.description` + producing software)
  into each figure FILE's metadata (best-effort; a backend that rejects it still writes the file).
- `scripts/verify.py` — data gates + `summarize` + the render/read-image QA loop
  (`audit_layout` also covers clipping, **data running outside the y-limits**, x+y tick overlap,
  composite panel letters, and text from one panel spilling into another). Plus `adjust_pvalues` / `multiplicity_check`
  (Holm/Bonferroni/BH correction for a FAMILY of p-values — annotate adjusted, not raw),
  `flag_outliers_mad` (robust median/MAD outlier flag for replicate band-metrics — the
  univariate analog of `chemometrics.diagnostics`; never raises, flag-don't-hide), and the
  **Tier-3 critic** `check_number_provenance` (static scan: a hand-typed figure-of-merit on
  the figure/in a caption that should be interpolated from the computed value — the desync-
  on-refit failure mode; keyed on the skill's FoM vocabulary so band identifiers are exempt).
- `scripts/spectra.py` — FTIR+PXRD ingest, baseline, integration, waterfall/overlay,
  optional Savitzky-Golay (`apply_smooth`), %T display (`as_transmittance` /
  `cfg.ftir_yaxis`), `edge_labels` (right-margin per-trace waterfall key, grayscale-safe;
  reusable for hand-built or composite panels), and `label_peaks` (styles the labels
  YOU supply — it does not pick peaks or write the text).
- `scripts/calibration.py` — regression, bands, residuals, LOD/LOQ, predict.
- `scripts/chemometrics.py` — multivariate calibration (PLS/PCR/PCA): leakage-safe
  `cross_validate`/`component_scan`, parsimony `choose_n_components`, `Preprocessor`
  (SNV/MSC/d1/d2/centre/autoscale), `vip`, ax-level panel plots, the `diagnostics_figure`
  composite, and `methods_text`. sklearn is a lazy dep. Plus the PCA-based **outlier gate**
  (`diagnostics`/`plot_influence`: T²/Q/leverage + DD-SIMCA; numpy+scipy only), a **one-class
  DD-SIMCA classifier** (`ddsimca_fit`/`ddsimca_predict`/`class_fom` — sensitivity/specificity/
  efficiency; the "is this the target phase?" call for cocrystal ID), **`emsc`** (EMSC
  scatter+baseline+interferents), and the **validation trio** (`permutation_test`,
  `corrected_paired_t`, `rpd`/`rpiq`). All validated by `skill_validation/chemometrics/` (100 checks).
- `scripts/charts.py` — generic line/bar with anti-pattern guards.
- `scripts/doe.py` — seeded, operator-balanced **run-order schedules** (operator orthogonal to
  concentration — fixes the operator/level confound) + CSV run sheet. Self-validates: `python
  scripts/doe.py` (determinism + balance + coverage).
- `scripts/report.py` — `methods_report(cfg, ...)`: paste-ready methods text + SI table.
- `scripts/crystal_engine.py` — CIF → validated crystal data (parse, symmetry-expand,
  density triple-check, geometry, H-bonds); the single source of truth for the crystal family.
- `scripts/crystal_pxrd.py` — calculated PXRD from a CIF (Dans_Diffraction + pymatgen/Bragg
  cross-check), plotted through the `spectra` pxrd domain. `reflection_list` extracts Dans's
  physically-correct powder reflections `[(2θ, I, hkl)…]`; `realistic_pattern` composes them with
  `pxrd_realism` (the turn-key Cu-Kα-scan simulator).
- `scripts/pxrd_realism.py` — make a calc pattern match a real lab Cu-Kα scan: `march_dollase`
  (preferred-orientation intensity), `kalpha2_doublet` (Cu Kα₁/Kα₂ splitting), `caglioti_fwhm` +
  `pseudo_voigt` (angle-dependent width/shape), composed by `simulate_pattern`. Numpy-only, so it
  operates on a reflection list `[(2θ, I, hkl)…]` (Dans supplies the structure factors) and is
  fully testable without heavy deps. Validated by `skill_validation/pxrd/` (physics, dep-light) and
  end-to-end on real CIFs (aspirin, lactose) via `crystal_pxrd.realistic_pattern`.
- `scripts/crystal_view.py` — deterministic 3D structure render (PCA or custom/vector/axis
  camera, bonds, dashed H-bonds, C–H hiding, heteroatom labels).
- `scripts/cocrystal.py` — cocrystal-ID arithmetic the `cocrystal_id.md` decision layer routes
  to. Wave 1: `delta_pka`/`classify_ionisation` (the ΔpKa salt/cocrystal/continuum call,
  Cruz-Cabeza 2012; zero-dep, never raises) + `rwp`/`sum_of_parents`/`phase_report` (NNLS
  new-phase-vs-physical-mixture discriminator with Rwp + new/lost-peak table; scipy, lazy).
  Validated by `skill_validation/cocrystal/` (26 checks).
- `templates/analysis_template.py` — copy this; shows the CONFIG-block shape.
- `bundle.py` — amalgamate the modules an analysis uses + its CONFIG/body into one
  self-contained deliverable script. Stamps a **provenance header** (skill version, the
  Python + scientific-stack versions, source, `--describe`/`DESCRIPTION` text) and runs the
  **Tier-3 critic** on the body at handoff (`--no-critic` to skip).
- `scripts/sources.py` — fetch reference data from the OPEN databases (stdlib-only, best-effort,
  cites the URL, retries transient DNS/timeouts): `cod_search`/`cod_fetch_cif` (Crystallography Open
  Database → CIF → `crystal_engine`) and `pubchem_properties` (PUG-REST identity/computed props) —
  both live-verified end-to-end. The full sourcing map —
  incl. the paywalled/manual ones — is `references/databases.md`. Pure URL/parse logic is gate-tested
  (`tests/test_sources.py`); the network wrappers aren't.
- `cli.py` — scripts-first QC shortcuts: `spc-rsd` (.spc replicates → band area/%RSD) and
  `calibrate` (conc,signal CSV → slope/R²/LOD/LOQ). Convenience for routine checks; the report
  deliverable is still a written analysis, not CLI output.
- `tests/` — pytest **skill self-tests** (module imports, SKILL.md structure, `bundle.py`
  round-trip, multiplicity vs statsmodels, `.mplstyle`/width, the CLI) that also bridge in the
  numeric suites. `python -m pytest tests/` is one green gate for the whole skill.
- `assets/styles/` — `.mplstyle` house-style presets per journal (nature/acs/ieee/general).
- `skill_validation/chemometrics/` — the numeric-parity validation suite (100 checks; see its SPEC.md).
- `skill_validation/cocrystal/` — the numeric-parity validation suite for the cocrystal-ID family
  (ΔpKa zones + NNLS sum-of-parents/Rwp, 26 checks), bridged into `pytest tests/`.
- `skill_validation/pxrd/` — the numeric-parity validation suite for PXRD realism (23 checks:
  March-Dollase, Kα₂ doublet, Caglioti, pseudo-Voigt, `simulate_pattern`), bridged into `pytest tests/`.
- `skill_validation/crystal/` — the CIF edge-case set (COD CIFs + synthetic broken inputs; CCDC extras
  optional) with the engine / independent-gate / ADP / realism checks (`run.py`; needs gemmi +
  Dans_Diffraction, pytest skips it otherwise) and `crystal_demo.py`, the worked example.
- `skill_validation/ftir_integration/` — the band integrator against published calibrations (pectin
  DM, BMC 2021) plus a deterministic closed-form ground-truth test (14 checks, bridged into `pytest tests/`).
- `references/` — `figure_selection.md` (the "what to plot" decision layer — start here),
  `cocrystal_id.md` (the "which evidence identifies a cocrystal" decision layer — start here
  for cocrystal ID), `databases.md` (where to source each external fact — open first,
  cite always), `verification.md`, `spectra.md`, `calibration.md`,
  `chemometrics.md`, `charts.md`, `crystal.md`, `adding_a_family.md`. View on demand, not all at once.

## When NOT to over-reach

Don't psychoanalyse the data or invent conclusions. Report what the gates and
the numbers show. If n is small, say the CIs are wide and lead with the trend.

**Build the figure that was asked for.** Produce the panels the request names; don't
annex an adjacent analysis, an extra panel, or a supporting figure nobody asked for.
Where a second figure looks genuinely warranted, say so in one line and let the user
decide — an unrequested panel costs a re-render and dilutes the one that answers the
question. The gates above are the exception: a residual panel, a named error bar and a
QA pass are part of the deliverable being *correct*, not scope creep.

**Length is calibrated to the artefact.** A caption is a caption; `methods_text` /
`methods_report` emit paste-ready prose, not a review. The standalone script carries the
analysis plus its provenance header — it does not need step-by-step narration in
comments, and the house comment rule applies (state a constraint the code can't show;
never restate what the next line does). Prose that isn't carrying information is prose
the reader skips.

**External facts get looked up, not recalled.** Any value that comes from *outside* the
data in front of you — a reagent/column price or availability, a literature reference
value, a pharmacopoeial limit, a band assignment you're quoting as authoritative — must be
fetched with a live web search (or read from a cited document) and reported **with its
source + date**, never stated from training memory. Prices and catalog availability drift;
a confidently recalled number that's stale silently corrupts a budget or a spec. Same
provenance discipline as the figure/number gates, applied to the world outside the dataset.
