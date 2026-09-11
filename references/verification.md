# verification.md — the two tiers, and how to actually use them

## Contents

- Tier 1 — data-pipeline gates (wrong NUMBERS)
- Tier 2 — visual QA loop (broken-LOOKING figures)
- Tier 3 — the critic (numbers that don't trace to the code)
- Provenance — the figure carries how it was made
- Composite figures

The skill's whole reason to exist. Other figure skills stop at styling; this one
gates the data and then makes you *look* at the render.

## Tier 1 — data-pipeline gates (wrong NUMBERS)

Run these BEFORE anything is plotted or integrated.

- `check_ingest(path, cfg)` — file exists and isn't truncated.
- `check_trace(x, y, cfg, name)` — finite, monotonic x, enough points, no lone
  spike. A FAIL here almost always means a decode/parse error (wrong column,
  concatenated files, mixed delimiters). Fix the ingest, don't plot around it.
- `same_processing(records, cfg)` — every rep in a batch was baseline/window/
  normalised identically. Build a small dict per trace recording how it was
  processed and pass the list. A FAIL means a rep-to-rep difference might be an
  artefact, not real.
- `summarize(values, cfg)` — returns mean, SD, SEM, t-based CI and a caption that
  NAMES the statistic + n + t-multiplier. Put that caption on the figure.
- `flag_outliers_mad(values, cfg)` — robust median/MAD outlier flag for a set of
  replicate band-metrics (a point is flagged when `|x − median| / (1.4826·MAD)` exceeds
  `cfg.outlier_mad_n`, default 3.5). The univariate, numpy-only analog of the
  multivariate `chemometrics.diagnostics` gate — use it to turn an ad-hoc outlier call
  (an under-loaded or a diluted replicate) into a principled, reportable rule. Robust *by
  construction*: a mean±k·SD rule is inflated by the very outlier it should catch; the
  median/MAD is not. It **never raises** (flag-don't-hide: the point is surfaced for you
  to retain / re-measure / footnote, never silently dropped) and WARNs when n<5 (the MAD
  of 3 reps is noisy — flag over a *pooled* set). **Idiom:** plot the raw replicate
  points, shade `[lower, upper]` as the accept band, mark the flagged points — so the
  reader *sees* why a point was called (the before/after threshold overlay).
- `adjust_pvalues(pvals, method)` / `multiplicity_check(pvals, cfg)` — when you annotate a
  FAMILY of p-values together (intercept-bias across levels, lack-of-fit across operators, several
  pairwise model comparisons), correct for multiplicity (Bonferroni / Holm / BH-FDR) and put the
  ADJUSTED p on the figure. Raw p-values across many comparisons inflate significance;
  `multiplicity_check` WARNs when the raw count exceeds the adjusted one. (Matches
  `statsmodels.stats.multitest.multipletests`.)

`cfg.strict=True` (default) turns any FAIL into a raised `GateError`. Set it
False only to triage; never to ship.

## Tier 2 — visual QA loop (broken-LOOKING figures)

1. `render_preview(fig, path)` -> a low-dpi PNG proof.
2. `audit_layout(fig, cfg)` -> deterministic faults a program *can* catch:
   missing-glyph tofu (both render warning channels, not just U+FFFD), off-canvas
   clipping of titles/labels/annotations, overlapping x AND y tick labels, empty
   axes, empty legends, unit-less numeric axes, and — on a composite — exactly one
   a/b/c letter per panel.
3. **Open the PNG with the Read tool** and walk `READ_IMAGE_CHECKLIST` — the
   perceptual faults no program sees: legend on top of data, indistinct-in-
   grayscale traces, FTIR axis not inverted, clipped text, curved residuals,
   default-matplotlib tells.
4. Fix, re-render, repeat until clean. ONLY THEN `save_fig` the vector master.

The point of step 3 is that the model critiques its own output by looking at it.
Skipping it is how a figure that "ran fine" ships with the legend over the peak.

## Tier 3 — the critic (numbers that don't trace to the code)

`check_number_provenance(source, cfg)` is a static reviewer pass over the analysis
**script** — not the data, not the pixels. It targets one silent failure: a
figure-of-merit written onto the figure or into a caption as a **hand-typed literal**
that has drifted from what the code computes. Write `"R2 = 0.98"` once, refit, and the
code now says `0.91` while the text still says `0.98`. The rule is *form*: a computed
FoM (R²/RMSE/LOD/LOQ/recovery/%RSD/p/n/slope/bias/RPD/±…) must be **interpolated** from
the value the pipeline produced —

```python
ax.set_title(f"R² = {fit['r2']:.2f}")     # traces to the code  ✓
ax.set_title("R² = 0.98")                  # hand-typed → WARN   ✗
```

It scans only the strings that reach a reader (matplotlib text-sink arguments and
assignments to a caption-ish variable) and keys on the skill's **own FoM vocabulary**,
so a fixed physical identifier — a band centre like `1748 cm⁻¹`, an axis label
`Recovery (%)` with no value — is a *constant*, not a result, and is left alone. It is
**WARN-only** (a lint, never a gate; it never raises). `bundle.py` runs it automatically
on the analysis body at handoff, so the deliverable is reviewed as it's produced. Adapted
from the **actor-critic reviewer** in Anthropic's Claude Science (June 2026), which flags
"untraceable numbers, and figures that don't match their underlying code."

## Provenance — the figure carries how it was made

Reproducibility is stamped, not narrated. `bundle.py` writes a **provenance header** into
the standalone deliverable (skill version · Python + numpy/scipy/matplotlib/… versions ·
source body · a one-line plain-language description from `--describe` or a `DESCRIPTION="…"`
constant), and `style.save_fig` embeds the skill version + `cfg.description` into each
figure FILE's metadata (PDF/SVG/PNG). So the "exact code and environment that produced it"
travels with the figure and is recoverable months later — ideal for the report appendix.

## Composite figures

A multi-panel figure is ONE figure: build it with `style.figure(cfg, nrows, ncols)`,
finish with `style.finalize_figure(fig)`, letter it with
`style.add_panel_labels(fig, cfg)` — in that order, so the layout is settled before
the a/b/c letters are anchored (each at its panel's corner, with one shared points
offset, so they line up across panels regardless of tick-label widths).
`audit_layout` then checks the letters are present and one-per-panel; the
Read-the-PNG pass checks they actually line up. The letter belongs in the image; the
figure number and caption belong in the document (Word/LaTeX cross-references), so
the image stays reusable and renumbering can't desync.
