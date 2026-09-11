# chemometrics.md — multivariate calibration (PLS / PCR / PCA)

## Contents

- When to reach for it (vs univariate `calibration.py`)
- The three correctness traps (this is why the family is gated)
- Preprocessing menu (`cfg.pls_preprocess`, chain with "+")
  - EMSC (`emsc`) — scatter + baseline + known interferents
- The four diagnostic panels (`diagnostics_figure`)
- What to report
- VIP
- Outlier / diagnostic gate (`diagnostics`, `plot_influence`)
- Validation trio — is it real / better / useful?
- Dependency

The conventions for the chemometrics family. The matplotlib is assumed known; what
this encodes is *what's right* — which error to report, what leaks, what to label.

## When to reach for it (vs univariate `calibration.py`)

Use `calibration.py` when ONE band is selective for the analyte: integrate it, regress
area on concentration, done. Reach for chemometrics when no single band carries the
quantity cleanly — overlapping bands, matrix interference, or scatter/sample-loading
variation that swings every band together. PLS regresses the *whole* spectrum and finds
the analyte direction in spite of those. It is the right tool exactly when the
univariate %RSD is loading-limited rather than detection-limited.

Don't reach for it to dress up a problem a single band already solves — a one-band
ratio that works is more transparent than a latent-variable model that works equally
well. Report both when you can.

## The three correctness traps (this is why the family is gated)

1. **CV leakage.** Preprocessing that uses cross-sample statistics — MSC's reference
   spectrum, mean-centring, autoscaling — must be fit on the TRAINING rows of each fold
   only. Fit it once on the whole set and the test rows have leaked into training, so
   RMSECV is optimistic. `cross_validate` / `component_scan` re-fit a fresh
   `Preprocessor` inside every fold; never pre-transform X and then cross-validate the
   transformed matrix. (SNV and derivatives are per-sample, so they can't leak — but
   the machinery is uniform so you can't get it wrong by accident.)

2. **Over-fitting the component count.** The global-minimum RMSECV almost always sits at
   too many latent variables; the extra ones fit CV noise. Take the parsimony pick
   (`choose_n_components`): the fewest components within `cfg.lv_parsimony_tol` of the
   global min. Show the RMSECV-vs-components curve so the choice is visible, not asserted.

3. **The optimistic-CV trap (name your question).** Leave-one-out over replicate spectra
   that share a concentration LEVEL answers *"predict a new sample at a level I've already
   measured"*. Predicting a *new* level is harder and is the relevant question for an
   unknown — use leave-one-LEVEL-out by passing `groups=level-labels`. The two RMSECVs
   differ (on the Group-C set, ~2.9 vs ~3.5 %w/w). Always state which you report;
   `cross_validate` WARNs if you LOO over replicated y without groups.

## Preprocessing menu (`cfg.pls_preprocess`, chain with "+")

| step | what it does | when |
|---|---|---|
| `snv` | per-row centre & scale | scatter / path-length / **sample-loading** variation (the usual win for ATR) |
| `msc` | regress each spectrum onto the training mean | same target as SNV; needs a reference (fit in-fold) |
| `d1` / `d2` | Savitzky-Golay derivative (window = `cfg.deriv_smooth`) | sloping baselines (`d1`); **resolving overlap** (`d2`) — but a 2nd derivative amplifies noise hard, so it often *loses* on precision-limited data |
| `center` / `autoscale` | subtract / standardise per-wavenumber (fit in-fold) | PCA/PCR; rarely needed for PLS (it centres internally) |

Choose the preprocessing the same way you choose the component count — from the RMSECV
scan, not taste. Run a few options through `component_scan(pre_options=(...))` and let the
curve decide.

### EMSC (`emsc`) — scatter + baseline + known interferents

`emsc(X, wavenumbers, reference=, poly_order=, interferents=)` is the principled upgrade over
SNV/MSC for multiplicative loading variation: it fits each spectrum as `b·ref + Σ a_k·z^k (+ Σ d_j·interferent_j)` (z =
axis scaled to [−1,1]) and returns `(x − baseline − interferents)/b`. Two powers:
(1) a fitted polynomial baseline *jointly* with the multiplicative term, and (2) **interferents** —
pass the **pure-excipient spectrum** (e.g. lactose) to model the matrix away. Caveat: with `reference=mean` + a
baseline term, `b` is **not** a physical dilution factor (the constant column eats part of the
offset); use a *pure* reference with `poly_order=0` to recover a true factor. Leakage: a fixed
(pure/external) reference is leak-free — apply before CV; the mean-reference variant is cross-sample,
so refit per fold if used inside `cross_validate`. Validated bit-for-bit against the biospectools/
Kohler gold (`skill_validation/chemometrics/test_emsc.py`).

## The four diagnostic panels (`diagnostics_figure`)

- **(a) RMSECV vs components** — one curve per preprocessing; ring the parsimony pick.
  Y-axis carries the concentration unit (`cfg.conc_unit`, e.g. `% w/w`) — RMSECV is a
  concentration, not dimensionless.
- **(b) regression-coefficient spectrum** — which wavenumbers drive the prediction.
  FTIR x-axis **inverted**. Overlay the pure analyte (scaled) so the coefficients can be
  read chemically. **Caveat:** SNV/autoscale couple channels, so their coefficients smear
  across the whole axis and aren't band-interpretable. To *see* the bands, draw a RAW
  (band-localised) model — pass `viz_pre="none"` — while still REPORTING the SNV error
  from (a). That split is intentional, not a contradiction; say so in the caption.
- **(c) scores** — the calibration's sample structure; colour by concentration to see the
  gradient. The colourbar carries the concentration unit. (Watch for operator/run
  structure hiding here — clustering by who-prepared-what is a confound, not signal.)
- **(d) loadings** — the latent variables themselves vs wavenumber. FTIR inverted.

`diagnostics_figure` builds this on a **constrained-layout** figure (so the scores
colourbar reserves its own space instead of colliding with panel d) and places the panel
letters per-column (left column cleared to the margin; right column tucked at its corner
so it isn't dragged out by panel d's wide decimal ticks). Both are baked in — don't
reintroduce a non-constrained figure here.

## What to report

RMSEC and calibration R² describe the FIT, not predictive ability — never quote them as
the method's accuracy. Lead with **RMSECV** (state the scheme), **R²(CV)**, and **bias**.
`methods_text` assembles a paste-ready paragraph with the scheme and the question named.

## VIP

`vip(model)` ranks wavenumbers by contribution (>1 ≈ above-average). Useful to confirm the
model leans on chemically sensible bands rather than an artefact region; report the top
few with their assignments, not the whole vector.

## Outlier / diagnostic gate (`diagnostics`, `plot_influence`)

When a spectrum is *wrong* (a bad press, an under/over-load, a contaminant) you want to KNOW,
and to separate "weird spectrum" from "extreme but valid concentration" — the exact call the
Group-C calibration needed for its 40 %#1 (under-loaded) and 50 %#3 (dilute-IS) points. This is
the **flag-don't-hide** tool, and ASTM E1655 expects these diagnostics for multivariate IR.

A PCA decomposition (numpy SVD, **no sklearn** — it's raw-spectrum QC you can run *before* any
calibration) gives four per-sample distances:

- **Hotelling T²** — score distance; "how extreme inside the model" (a true high/low concentration
  is legitimately extreme). `T² = (n−1)·leverage`; `mean(leverage) = A/n`.
- **Q / SPE** — orthogonal distance; "how badly it fits the model" — the usual **bad-load / new-feature** tell.
- **leverage** (Σ (t/s)²) — the unknown-sample **extrapolation** flag: a new sample above the
  calibration max (or the 2A/n · 3A/n line) is out of domain.
- **DD-SIMCA** (Pomerantsev & Rodionova 2014) — the modern verdict: SD + OD combined into a
  χ²(Nh+Nq) statistic with an **extreme** limit (α) and a per-sample **Bonferroni outlier** limit (γ).
  This is the primary classification; T²(F) and Q(Jackson–Mudholkar) are the classical cross-checks.

`plot_influence(ax, diag)` is the standard **acceptance plot** — OD/q0 vs SD/h0 with the χ² boundary
as a straight line (extreme dashed, outlier dotted), points keyed regular/extreme/outlier by
marker+colour (grayscale-safe). Pass `labels=` to tag only the flagged points. Use
`axis_scale="sqrt"` when one or a few **gross** outliers dominate the range and squash the
in-control cluster + limit lines into the corner (the single-bad-spectrum case): a √ axis keeps
the borderline/extreme points and the boundary legible. Keep the default `"linear"` for moderate
spread (true linear boundary, easiest to read).

Conventions (pinned, see `skill_validation/chemometrics/`): SD uses the singular-value normalization
(DDSimca/mdatools; `sum(leverage)=A`); DoF are method-of-moments (ddof=1); limits fall back to χ²
when `n−A≤0` and the Q test reports 0 ("uninformative") when there is no residual space (A≥rank).
Every limit form is validated bit-for-bit against a hand-derived fixture + scipy, and the Hotelling
limits against published R `mdatools` numbers.

## Validation trio — is it real / better / useful?

Cheap, defensible figures-of-merit for a SMALL calibration (e.g. 18 spectra):

- **`permutation_test`** — refit on **scrambled y** to build a null for R²(CV); `p = (#{null ≥
  observed}+1)/(n_perm+1)`. With `groups=` it permutes **between LEVELS** (each level's reps keep
  one shared y) — the right null for replicate spectra (sklearn's `groups=` permutes *within*, which
  is wrong here). Guards against a chance correlation.
- **`corrected_paired_t`** — the **Nadeau–Bengio** corrected resampled paired t-test: the naive
  paired-t variance is inflated by `(1/k + n_test/n_train)` for the overlapping training sets. This
  is what turns "SNV-PLS 2.9 % vs matrix-IS 5.8 %" into a *defensible* difference. (mlxtend's
  `paired_ttest_resampled` is the **un**corrected formula — don't use it as the reference.)
- **`rpd` / `rpiq`** — RPD = SD(ref)/RMSE (the one-number "is the calibration useful": >2 useful,
  >2.5 good, >3 excellent, Chang 2001); RPIQ = IQR/RMSE (robust to skew, ≈1.349·RPD for normal y).

All validated in `skill_validation/chemometrics/test_trio.py` (the Nadeau–Bengio golden is
reproduced *through* the function).

## Dependency

scikit-learn loads **lazily** (only when a model is fit) — an FTIR or univariate job never
imports it. A bundled deliverable keeps the import indented inside the functions, so the
file is self-contained but only pulls sklearn at run time. If it's missing the error names
the `pip install`. The **outlier gate**, **EMSC**, and the **validation trio** use only numpy +
scipy (no sklearn).
