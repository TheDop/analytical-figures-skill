# spectra.md — FTIR/IR/ATR and PXRD/XRD recipes

## Contents

- FTIR / IR / ATR
  - Line style: spectra are SOLID by default (house rule)
  - Choosing baseline anchors WITHOUT hand-picking them
  - Robustness: report it, don't argue it
- PXRD / XRD
- Gotchas

Set `cfg.domain` and the axis/intensity conventions flip automatically.

## FTIR / IR / ATR
- x = wavenumber cm-1, **axis inverted** (high -> low). `x_limits=(1820,1650)` is
  already in inverted order.
- y = absorbance (or normalised).
- Baseline: `arpls` (tune `arpls_lam`; bigger = stiffer) or `rubberband`.
- Band integration: put `(lo, hi, name)` windows in `cfg.integration_windows`.
  A **local linear baseline** is drawn between the two anchors and applied
  IDENTICALLY to every spectrum — that identical application is what makes a
  %RSD across reps meaningful. Compute areas with `integrate_bands`, then feed
  the per-rep areas to `verify.summarize` for mean +/- CI and %RSD.
- Baseline mode (`cfg.integration_baseline`): `"per_window"` (default) gives each
  window its own local baseline — correct for **isolated** bands. For **overlapping
  bands on a shared pedestal** (e.g. an analyte band right next to an
  internal-standard band, or two carbonyl bands that don't return to baseline
  between them) use `"shared"`: ONE baseline is drawn across the whole envelope
  (leftmost..rightmost anchor) and each window integrates its sub-interval above
  that single line — i.e. a vertical-drop / "drop perpendicular" partition at the
  window boundaries (set the boundary by where two adjacent windows meet). A
  per-window baseline on overlapping bands follows the valley *up* and carves area
  off the smaller band — which biases a band-area *ratio* and flattens its
  calibration. (Verified against the Wang 2023 pectin DM dataset: the ester/
  carboxylate area ratio reproduces only under the shared baseline.)
- Idiom: `plot_overlay` for reps on one axis; `plot_waterfall` for groups offset
  vertically (e.g. as-received vs ground vs heated).

### Line style: spectra are SOLID by default (house rule)
**Spectra use continuous (solid) lines unless a dashed/dotted style is explicitly
requested.** Dash/dot patterns fight the data's own high frequency (an IR fingerprint
is full of close peaks), read as clutter, and are NOT the default redundancy channel
for spectra. Get grayscale-safety another way:
- **Waterfall / offset (the default idiom):** traces are separated by POSITION, so
  solid lines are already grayscale-safe — `plot_waterfall` forces solid. For many
  busy spectra an offset waterfall is usually cleaner than any overlay. When a factor
  (e.g. concentration) subdivides the stack, offset the *groups* into tiers and
  overlay the within-group reps on top of one another (solid, colour = the sub-factor).
- **Overlay of a few traces:** rely on colour + a direct label/legend (and, where the
  bands sit in different regions, position itself identifies them). If you genuinely
  need a second channel on an overlay, prefer `cfg.redundancy="marker"` (solid lines +
  sparse `markevery` markers) over dashes.
- **Dashes are opt-in only** — `cfg.redundancy="linestyle"` remains available for a
  smooth single-band overlay *if the user asks for it*, but it is never the default.
- **Identifying waterfall traces** (`cfg.waterfall_legend`): default `"edge"` puts
  a direct label at the right margin aligned to each trace — grayscale-safe and no
  floating text over the data, the house default (preferred over a colour-only
  legend). `"legend"` draws a legend box; `"none"` leaves identity to the caption.

### Choosing baseline anchors WITHOUT hand-picking them
Typing window edges by eye is the irreproducible step: two analysts get different
%RSD from the *same* spectra purely from where they put the anchors (we saw a 2.5x
height difference and a 4% vs 2% precision split this way). Three fixes, in order
of how much they remove the judgement:

1. **`band_metric(x, y, "deriv2", center=...)`** — 2nd-derivative trough depth.
   ANCHOR-FREE: a constant+linear background has zero 2nd derivative, so a sloping
   continuum (e.g. an analyte band on a neighbour's wing) cancels. Pays in noise →
   computed on an SG derivative (`cfg.deriv_smooth`); weakest at trace levels.
   Notes that bite: `savgol_filter(deriv=2)` *is* Savitzky-Golay smoothing applied
   in the same pass as the derivative — do NOT SG-smooth and then finite-difference
   (worse). The SG **window** is the real knob (a 2nd derivative amplifies noise):
   set it to ~1-2x the band FWHM — too narrow is noisy, too wide smears the band and
   throws away the resolving advantage that justified the derivative. poly 2≡3 for
   d2 (and 4≡5). If the optimum window is several× the FWHM, the data is noise- not
   slope-limited → a linear-baseline area is the better-justified choice. (Example, the
   aspirin ester C=O: deriv2 needed ~80-95 cm⁻¹ to reach area's R²/LOD → area kept primary.)
2. **`find_anchors(x, y, center, gap, maxhw)`** → the flanking local minima (the
   continuum-return valleys) found as the argmin of an SG-smoothed copy in a
   bounded side-window. Reproduces a careful analyst's choice deterministically.
   Then `band_metric(..., metric="area"|"height", anchors=...)`.
3. **`baseline` = arpls/rubberband** (global) then integrate a fixed window — no
   local anchors at all.

**Detect-once-and-LOCK.** Re-detecting anchors per replicate re-introduces drift
and breaks identical-processing. Detect once on a reference (pure-analyte trace or
batch mean) and reuse:
```python
cfg.anchor_centers      = [(1748.0, "ester")]     # band(s) to find
cfg.integration_windows = spectra.auto_windows(ref_x, ref_y, cfg)   # lock here
# ...then process the whole batch with those locked windows
```

### Robustness: report it, don't argue it
There is no single "true" baseline — it's a *method parameter*, and ICH Q2 calls
its influence **robustness**. Quantify it instead of choosing in private:
- `verify.baseline_robustness(reps, center, cfg, reference=ref)` — per-method %RSD
  for one replicate set; WARNs if methods disagree >3x.
- `verify.calibration_robustness(levels, center, cfg, check=..., reference=ref)` —
  runs the WHOLE calibration under area/height/deriv2 and tabulates slope, R²,
  LOD/LOQ, recovery so the locked primary metric is chosen from evidence. The
  spread it prints IS the robustness result you report.
(Example, aspirin/lactose ATR: area on locked 1713/1789 anchors → best R²/LOD/precision,
recovery ~100%; all three metrics agree the grinding %RSD trend, which is the point.)

## PXRD / XRD
- x = 2theta degrees, **axis normal** (low -> high). y = intensity (counts / norm.).
- Background subtract; watch for the amorphous halo (don't integrate it as signal).
- Phase ID: put reference reflections in `cfg.pxrd_reference_lines` as
  `(2theta, label)`; `plot_waterfall` draws them as dotted stick lines.
- Idiom: stacked/offset patterns + reference sticks. Never invert the 2theta axis.

## Gotchas
- Don't mix baseline methods within a batch — `same_processing` will FAIL.
- Normalisation is for display; integrate on the baseline-corrected (not
  normalised) trace if you want absolute areas.
