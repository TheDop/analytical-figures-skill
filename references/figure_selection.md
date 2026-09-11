# figure_selection.md — decide WHAT to plot before HOW

## Contents

- How to use it (3 questions, then route)
- Axis 1 — data shape
- Axis 2 — the argument (intent → figure)
- Axis 3 — regime / scale
- Data shape → figure quick table (preferred · alternative · never)
- Same data, different argument → different figure
- When to split into multiple figures
- Figure-type semantic boundaries (the "never", with the reason)
- Calibration & chemometrics — the conventions that drive these figures
- Final form (output spec) — apply last
- What this doc does NOT do
- Sources

The skill's "start here" for any figure. Modelled on scipilot's `chart_selection.md`
(Haojae/scipilot-figure-skill, MIT), re-grounded for analytical chemistry and tied to
citable standards.

**Why this exists.** A capable model already chooses a sensible chart ad hoc — so the
value here is **consistency**: one repeatable, *grounded* answer instead of a slightly
different judgement each session. This is **form, not content** (cf. SKILL.md "Form is
the skill's job; judgment is the model's"): it gives you a shared reasoning framework and
the citable rule behind each call. It does **not** pick which band to integrate, write the
caption, or decide the scientific narrative — that stays yours.

**It routes, it doesn't re-implement.** Once you've chosen, the *how* lives elsewhere:
family modules (`spectra` / `calibration` / `chemometrics` / `charts` / crystal), the
house-style lock (`style.py`), and the gates in `verification.md`. Rules already enforced
in code (low-n bar guard in `charts.md`; palette + font floor + vector export in
`style.py`; residual panel + no-extrapolation in `calibration.py`) are pointed to, not
restated.

---

## How to use it (3 questions, then route)

1. **What SHAPE is the data?** (§ Axis 1) — replicate spectra? a calibration? a variance
   study? a whole-spectrum matrix? a CIF/pattern? generic lab quantities?
2. **What ARGUMENT is the figure making?** (§ Axis 2 — the axis everyone skips) — identity,
   precision, linearity, accuracy, detection, error-diagnosis, multivariate structure,
   phase ID, or method-ranking? *Same data + different argument → different figure.*
3. **What REGIME/SCALE?** (§ Axis 3) — n per group, dynamic range/saturation, is any single
   band selective (the univariate-vs-chemometrics fork), how many spectra for a model?

Then read off the **quick table**, confirm against the **semantic-boundary "never" list**,
build with the family module, and size/colour via `style.py`. Apply the **output spec**
(§ Final form) last.

---

## Axis 1 — data shape

| Shape | Typical in this project |
|---|---|
| One spectrum | qualitative ID / specificity check |
| A few replicate spectra of one sample | precision / homogeneity (%RSD) |
| A batch across a gradient (levels × reps) | calibration set, grinding-time series |
| One band-metric vs known concentration | univariate calibration |
| Whole-spectrum matrix vs a property | chemometrics (no single selective band) |
| Repeated measurements, multi-factor | error/variance study (operator × loading × instrument) |
| A CIF / calculated pattern vs references | crystal structure / phase (cocrystal) ID |
| Generic lab quantities | time series, group comparison, composition |

## Axis 2 — the argument (intent → figure)

This is the axis that changes the figure even when the data don't. Map intent → figure:

| The argument you're making | Figure |
|---|---|
| **Identity / specificity** — "the band is here; lactose is flat at the C=O; the ester position is invariant" | annotated overlay; mark the diagnostic window |
| **Precision / homogeneity** — "reps agree to X % RSD" | overlaid reps (zoom to band) + %RSD in caption; variance bar |
| **Linearity** — "response ∝ concentration" | calibration scatter + fit + CI/PI **+ residual panel** |
| **Accuracy / bias** — "recovery ≈ 100 %; is the intercept ≠ 0?" | recovery vs added conc with 100 % line; bias plot; intercept test |
| **Detection** — "can we see it; LOD/LOQ" | low-end calibration + blank SD; distinguish decision limit L_C from capability L_D |
| **Diagnosing the limiting error** — "it's loading/scatter, not particle size" | total-absorbance swing vs constant FWHM; IS-ratio flat vs raw swinging |
| **Multivariate structure** — "LV1 = the aspirin–lactose contrast" | RMSECV-vs-LV, scores, loadings, predicted-vs-reference |
| **Phase identification** — "new peaks ⇒ cocrystal, not a physical mixture" | calc-vs-starting-materials PXRD waterfall + new/lost-peak diff |
| **Ranking methods/strategies** — "no-IS < matrix-IS < SNV-PLS" | grouped FoM bars, direction-of-good in the title |

## Axis 3 — regime / scale

- **n per group → how much you may summarise.** Our calibration is **3 reps per level — small**, which lands at the boundary below: prefer showing every point.
  | n per group | Use | Avoid |
  |---|---|---|
  | n < 3 | plot every point | box / violin / mean bar (no meaningful summary) |
  | 3 ≤ n < 10 | strip / dot / scatter; box cautiously | **mean-only bar** (hides the distribution) |
  | 10 ≤ n < 30 | box / violin **+ overlaid points** | mean-only bar |
  | n ≥ 30 | box / violin / mean ± CI bar all OK | — |
  Basis: small samples should be shown in full, not as mean ± error [Weissgerber 2015]; at n = 3 "it is better to simply plot the individual data points" [Cumming 2007]. The `charts.bar` guard already refuses bare mean bars below `min_n_for_mean_bar`.
- **Dynamic range / saturation.** ATR absorbance **can saturate** (e.g. above ~50 % w/w for a strong band in a binary mixture) — don't fit or predict through the saturated region; never extrapolate past the calibrated range [ICH Q2(R2) §3.2.1]. Use a log axis only when the variable spans decades.
- **The univariate-vs-chemometrics fork (selectivity test).** If a band is **selective for the analyte and loading-corrected** → univariate calibration on its area. If **no single band is selective**, or scatter/loading swamps every band together → regress the whole spectrum (chemometrics), and report **RMSECV/RMSEP, not calibration R²** [ASTM E1655].
- **Spectra count for a model.** Small sets (~20 spectra) → **leave-one-LEVEL-out** CV (not leave-one-rep-out — that answers an easier question), parsimony LV count, and a permutation/y-scramble check; state #LVs and #calibration samples in the caption [ASTM E1655: ≥ 6 × (factors + 1)].

---

## Data shape → figure quick table (preferred · alternative · never)

| Data shape | Preferred | Alternative | Never |
|---|---|---|---|
| One spectrum, qualitative ID | annotated trace, key bands labelled | — | unlabelled trace; rainbow fill |
| Replicate spectra, show agreement | overlay zoomed to band + %RSD | offset waterfall | mean spectrum alone (hides spread) [Weissgerber 2015] |
| Batch across a gradient | waterfall ordered by level + per-level %RSD | small multiples | overlay-all (occludes); colour-only legend |
| Calibration (band vs conc) | scatter + OLS fit + CI/PI **+ residual panel** | standard-addition | **R² alone as the linearity claim** [Eurachem]; extrapolation [ICH §3.2.1] |
| Recovery / accuracy | % recovery vs added conc, reps, 100 % line | bias (pred − true) vs conc | a single level called "accuracy" |
| Whole spectrum, no selective band | chemometrics diagnostics (RMSECV/scores/loadings/pred-vs-ref) | VIP overlay | calibration R² as accuracy [ASTM E1655]; forcing a univariate band |
| Variance / error study | per-factor box/strip + variance-hierarchy bar | nested dot plot | mean-only bars at small n [Weissgerber 2015] |
| Method / strategy comparison (FoM) | grouped bars, direction-of-good in title, raw points | slope / dumbbell | bars where "taller" is ambiguous |
| Composition (parts of a whole) | stacked / 100 % bar | — | **pie / 3D pie** [Cleveland & McGill 1984; Wong 2010] |
| Two continuous vars (relationship) | scatter + fit + CI | hexbin / 2D-KDE (large n) | connecting line if x is unordered |
| Time / ordered sequence (trend) | line + error band | step | bars |
| Calc PXRD vs references (phase ID) | waterfall calc vs starting materials (+ exp) + diff | stick overlay | a lone pattern, no reference |
| Matrix (correlation / confusion) | heatmap, perceptually-uniform cmap | annotated table | rainbow/jet cmap |

---

## Same data, different argument → different figure

A **0–50 % w/w calibration: 6 levels × 3 reps = 18 spectra** can become four
unrelated figures depending on the claim — fix the argument first:

| Argument | Figure |
|---|---|
| "It's linear" | calibration + residual panel + lack-of-fit / Mandel |
| "It's **precision-limited, not non-linear**" | within-level %RSD overlay / variance bar |
| "The error is **loading, not particle size**" | total-absorbance swing (20–30 %) vs ~constant FWHM |
| "An IS / **SNV-PLS** fixes it" | IS-ratio flat vs raw swinging + RMSECV ranking |

One dataset, four figures. (Mirrors scipilot's drug-A/B example, in the analytical-chemistry domain.)

---

## When to split into multiple figures

Split if **any** holds: dimension combinations > ~12 panels; x-tick labels collide
(> ~8 categories needing rotation); legend > ~6 entries (beyond quick recall); the y-scale
spans orders of magnitude and can't go log; or you're trying to make **two arguments** in
one figure (→ two figures). Split by grouping dimension, by panel, by argument, or push
detail to SI. Build multi-panel as ONE figure (see SKILL.md "Composite / multi-panel"),
never by pasting PNGs.

---

## Figure-type semantic boundaries (the "never", with the reason)

- **No pie / 3D / bubble.** Quantities should be encoded by **position on a common scale**
  (or length); angle, area and volume rank lowest in perceptual accuracy — replot as a bar
  [Cleveland & McGill 1984; Wong 2010 "Design of data figures"]. (Already flagged in `charts.md`.)
- **No rainbow / jet colormaps** — not perceptually uniform, not CVD-safe; use viridis-family
  (sequential) or a diverging map centred at 0. (Already in `charts.md`.)
- **Don't connect unordered points with a line** — a line asserts continuity/trend the data
  may not have; connect only when x is ordered (time, concentration) [Weissgerber 2015;
  Rougier 2014 R7 "Do Not Mislead"].
- **Dual-y axes — avoid** (not an absolute ban): comparing two non-aligned scales is
  perceptually weak and easy to manipulate into a spurious correlation [Cleveland & McGill
  1984 (rank-2 penalty); Rougier 2014 R7]. Prefer two panels or normalised overlay.
- **No truncated / non-zero bar baselines** — exaggerates small differences. (In `charts.md`.)
- **No mean-only "dynamite" bars** for small/continuous data — show the points
  [Weissgerber 2015]. (Enforced by the `charts.bar` low-n guard.)
- **Never colour alone** — add a redundant channel (line style, marker, direct label) so the
  figure survives grayscale and colour-vision deficiency (~8 % of men) [Wong 2011 "Color
  blindness"]. (`style.py` cycles colour **with** linestyle/marker.)
- **No positional jitter on a quantitative axis** — a 40 % point drawn at 42 reads as 42 %;
  separate overlapping markers on the *perpendicular* axis or make the axis categorical. (SKILL.md hard principle.)
- **Always name the error bar** (SD vs SEM vs CI — a √n factor) and **n**; unnamed bars are
  meaningless [Cumming 2007 R1]. CI/SEM overlap is **not** a significance test — the threshold
  depends on n and bar type [Cumming 2007 R6–7]. (`verify.summarize` returns the named statistic.)
- **Every axis carries a unit** (or `a.u.` / `ratio` / `%` / `dimensionless`). (`audit_layout` WARNs.)

---

## Calibration & chemometrics — the conventions that drive these figures

- **A residual panel is mandatory; judge linearity by the residual pattern, not R².** Random
  scatter about zero = linear; a systematic trend = lack of fit / wrong model / changing
  variance [ICH Q2(R2) §3.2.2.1; Eurachem §6.3.4 + Quick Ref 5]. (`calibration.py` builds it.)
- **≥ 5 levels, evenly spaced, replicated** (2–3× per level → the error bars) [ICH §3.2.2.1;
  Eurachem §6.3.4].
- **Stay inside the validated range**; extrapolation needs explicit justification [ICH §3.2.1].
- **LOD = 3.3 σ/S, LOQ = 10 σ/S**, σ = blank SD or regression residual SD; state which
  [ICH §3.2.3.3]. Distinguish the **decision limit L_C** from the **detection capability
  L_D ≈ 3.29 σ** on any detection figure [IUPAC Currie 1995].
- **Chemometrics:** report **RMSECV / RMSEP** on held-out data (not calibration R²); show a
  **leverage / Mahalanobis vs spectral-residual** outlier-diagnostic plot (a sample above the
  calibration's max leverage is an out-of-domain extrapolation); **name the preprocessing and
  spectral window** in the caption [ASTM E1655].

---

## Final form (output spec) — apply last

`style.py` already sizes at final dimensions, holds the font floor, cycles a CVD-safe palette,
and exports vector + raster (refusing JPEG). Targets:

- **Cross-journal safe default:** single column **~85 mm (3.3 in)**, double **~170–178 mm**;
  max height **~230 mm**; **sans-serif (Helvetica/Arial)**; body text **≥ 7 pt** at final size
  (hard floor 5 pt Nature/Science, 4.5 pt ACS); **panel letters 8 pt bold**; **line weight
  ≥ 0.5 pt**; **vector PDF/EPS master + ≥ 300 dpi raster proof**; **RGB**.
- **Our most likely target is ACS *Analytical Chemistry*** (single 3.33 in / 240 pt, double up
  to 7 in / 504 pt; 4.5 pt / 0.5 pt floors; B/W line art 1200 dpi, colour 300 dpi). Set
  `cfg.journal = "acs"` — `style.py` carries `nature`/`acs`/`ieee`/`general` width+panel profiles
  (the 6 pt font floor stays, stricter than ACS's 4.5 pt). For B/W line art at 1200 dpi set
  `cfg.dpi_raster` accordingly, or rely on the vector PDF master.
- **Biggest per-journal divergence:** Nature wants **editable vector, no TIFF**; RSC wants
  **flattened TIFF at 600 dpi**. Default to vector PDF/EPS (passes Nature/Science/ACS); export a
  600 dpi TIFF only for RSC.

---

## What this doc does NOT do

It doesn't choose the band/window, write the caption, or decide the scientific story — that's
your judgement. It standardises the *form* decision and hands you the citable reason. If a
case isn't covered, reason from the three axes and the "never" list; don't force a fit.

---

## Sources

Graphical perception & integrity:
- Cleveland WS, McGill R. 1984. *Graphical Perception…* JASA 79(387):531–554. doi:10.2307/2288400
- Cumming G, Fidler F, Vaux DL. 2007. *Error bars in experimental biology.* J Cell Biol 177(1):7–11. doi:10.1083/jcb.200611141
- Krzywinski M, Altman N. 2013. *Error bars.* Nat Methods 10(10):921–922. doi:10.1038/nmeth.2659
- Weissgerber TL, et al. 2015. *Beyond Bar and Line Graphs.* PLoS Biol 13(4):e1002128. doi:10.1371/journal.pbio.1002128
- Wong B. 2010. *Design of data figures.* Nat Methods 7(9):665. doi:10.1038/nmeth0910-665
- Wong B. 2011. *Color blindness.* Nat Methods 8(6):441. doi:10.1038/nmeth.1618
- Rougier NP, Droettboom M, Bourne PE. 2014. *Ten Simple Rules for Better Figures.* PLoS Comput Biol 10(9):e1003833. doi:10.1371/journal.pcbi.1003833
- Tufte ER. 1983/2001. *The Visual Display of Quantitative Information.* Graphics Press.

Calibration / validation / chemometrics:
- ICH Q2(R2). 2023. *Validation of Analytical Procedures.* (linearity/residuals §3.2.2.1; range §3.2.1; LOD/LOQ §3.2.3.3; multivariate §3.2.2.3)
- Eurachem. 2014. *The Fitness for Purpose of Analytical Methods*, 2nd ed. (§6.2–6.3, Quick Ref 5)
- Currie LA (IUPAC). 1995. *Nomenclature… Detection and Quantification Capabilities.* Pure Appl Chem 67(10):1699. doi:10.1351/pac199567101699
- ASTM E1655-17. *Standard Practices for Infrared Multivariate Quantitative Analysis.* (verify clause text before quoting in a submission — paywalled)

Journal figure specs:
- ACS *Analytical Chemistry* author guidelines / "Preparing Graphics" · Nature "Guide to Preparing Final Artwork" · Science author figure-prep guide (2022) · RSC "Figures, graphics & images".

Scaffolding adapted from scipilot-figure-skill (Haojae/scipilot-figure-skill, MIT).
