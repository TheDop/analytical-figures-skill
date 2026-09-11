# analytical-figures — roadmap / ideas backlog

A planning doc, not operational skill context (lives at the skill root so it is NOT
auto-loaded; SKILL.md and `references/*` are the live docs). Captures ideas mined
from the wider Claude/Agent-Skills ecosystem + mature chemometrics & crystallography
tooling (research pass, 2026-06-23), prioritised for the project the skill was built in: an
FTIR API-quantification study (aspirin in lactose) followed by a cocrystal-identification
study. Items marked DONE were built and validated; the rest is backlog.

## Strategic finding (why this list is "features", not "redesign")

Across the public ecosystem — "Chemist Analyst", the big collections (K-Dense-AI
`scientific-agent-skills`, InternScience, Microck), scicraft, and the chemometrics
libraries — **nothing gates the data pipeline or ships a self-contained reproducible
script the way this skill does.** Others are persona/"interpret this data" skills,
library-usage guides, or a figure-style *menu* that emits a rendered PNG. So: borrow
*features and conventions*, keep our architecture (single CONFIG, verify gates,
bundle-to-standalone, locked house style).

The one genuinely architectural idea worth adopting is the **decision module** below.

---

## ★ Marquee item — a "What figure should this be?" decision module

**Source of the idea:** `scipilot-figure-skill` (https://github.com/Haojae/scipilot-figure-skill),
the skill this one descends from. Its `references/chart_selection.md` ("the brain":
*decide what to plot before deciding how*) + companion `references/data_profiling.md`
(reads an auto-profiler's facts → routes each to a plotting decision).

**Why it is wanted:** a capable model already answers "what chart
for this data?" well without a skill — so the value is **consistency**: one repeatable,
grounded source of that answer, instead of re-deriving it (slightly differently) each
session. This is *form*, not *content*: it's a shared reasoning framework the model
consults, not a lookup that overrides judgment — so it fits "Form is the skill's job;
judgment is the model's." It also consolidates rules we already enforce but that are
scattered across SKILL.md (no pie, no positional jitter, units on every axis, bar-vs-box
n-thresholds, direction-of-good titles, line-vs-scatter) into one "start here" entry point.

**STATUS: written 2026-06-23** — `references/figure_selection.md` exists, grounded in citable
sources (Weissgerber/Cumming/Cleveland-McGill/Wong/Rougier for viz; ICH Q2(R2)/Eurachem/IUPAC
Currie/ASTM E1655 for calibration; ACS/Nature/Science/RSC for output specs), and wired into
SKILL.md as Workflow step 0. The sketch below records the design.

**What to copy:** the *structure*, not scipilot's generic-bio rules. Adapt to analytical
chemistry. Proposed `references/figure_selection.md`:

- **Three decision axes** (adapted to our domain):
  1. *What data do you have* — single/few replicate spectra · a batch across a gradient
     (levels × reps) · one band-metric vs known conc (calibration) · whole-spectrum
     matrix vs property (chemometrics) · a nested error study (operator × loading ×
     instrument) · a CIF / calc pattern vs references (phase ID) · generic lab quantities.
  2. *What argument is the figure making* (the axis everyone skips) — **identity/
     specificity** · **precision/homogeneity** · **linearity** · **accuracy/bias** ·
     **detection (LOD/LOQ)** · **diagnosing the limiting error** (loading/scatter vs
     particle size) · **multivariate structure** · **phase identification** ·
     **ranking methods/strategies**.
  3. *Scale/regime* — n per level (ours is 3 reps → small: show all points, CIs wide,
     never a mean-only bar) · dynamic range & **saturation** (ATR saturates >50% → don't
     extrapolate) · **single-band selectivity test** → the univariate-vs-chemometrics fork
     · spectra count for chemometrics (→ leave-one-level-out, parsimony LVs, permutation).
- **A data-shape → figure quick table** (preferred / alternative / never), e.g.:
  - replicate spectra, show agreement → overlay (zoomed to band) + %RSD · *alt* waterfall ·
    *never* mean spectrum alone (hides rep spread)
  - calibration → scatter + fit + CI/PI + **mandatory residual panel** · *never* R² alone
    as the linearity claim; *never* extrapolate past the calibrated range
  - variance study → per-factor box/variance-hierarchy bar + raw points · *never* mean-only
    bars (n small)
  - method/strategy comparison (FoM) → grouped bars with **direction-of-good in the title** ·
    *never* bars where "taller" is ambiguous
  - whole spectrum, no selective band → chemometrics diagnostics (RMSECV/scores/loadings) ·
    *never* force a univariate band
  - phase ID → calc-vs-references waterfall + new/lost-peak diff · *never* a lone pattern
- **"Same data, different argument → different figure"** worked on a real dataset — an 18-spectrum
  0–50 % w/w binary-mixture set: "it's linear" (calibration+residual+lack-of-fit) vs "it's precision-limited, not
  nonlinear" (within-level RSD / variance bar) vs "the error is loading, not particle size"
  (total-absorbance swing vs constant FWHM) vs "an IS / SNV-PLS fixes it" (IS-ratio flat vs
  raw swinging + RMSECV ranking). One dataset, four figures — fix the argument first. (The module
  makes this the default instead of something re-derived each session.)
- **Figure-type semantic boundaries** — fold in the rules we already hold (line vs scatter;
  bar vs box; no pie; no positional jitter on a quantitative axis; no rainbow/jet colormaps;
  units on every axis) so they live in one place.

**Optional companion (secondary, medium effort):** an analytical *batch profiler* analog of
`profile_data.py` — given a batch (levels × reps), emit *facts* (n/level, within-level %RSD,
dynamic range, single-band selectivity score, saturation flags, operator/level confound check)
+ a suggested figure/method, which `figure_selection.md` then combines with intent. We already
have `verify.check_ingest`/`check_trace`/`summarize` to build on. Lower priority than the
reference doc itself.

**Effort:** the reference doc is a writing task (no code) — do this first. The profiler is a
medium code add.

---

## ★ Claude Science pass (2026-07-03)

Anthropic launched **Claude Science** on 2026-06-30 (a research workbench — "Claude Code for
scientists": 60+ curated skills/connectors, code execution + HPC, a reviewer agent). First-party
skills are open at **`anthropics/life-sciences`** (marketplace + `single-cell-rna-qc`,
`instrument-data-to-allotrope`, `nextflow-development`, `scvi-tools`, `scientific-problem-selection`).
Strategic finding **holds**: their skills are QC/pipeline + best-practices packages; none gates the
data pipeline or ships a bundle-to-standalone the way we do. So we borrowed *features*, not architecture.

**★ Built this session (all three, validated — `tests/test_critic.py`, 43-test suite green):**
- **Number-fidelity critic (Tier 3)** — `verify.check_number_provenance`: static AST scan flagging a
  figure-of-merit hand-typed onto the figure/caption instead of interpolated from the computed value
  (the desync-on-refit failure). Keyed on our FoM vocabulary so band identifiers are exempt; WARN-only;
  `bundle.py` runs it at handoff. **This is Claude Science's actor-critic reviewer, done our way** (static,
  deterministic). Also the guardrail that makes the EDA idea below safe.
- **Provenance stamp** — `bundle.py` header (skill version + Python/stack versions + source + `--describe`/
  `DESCRIPTION` text) and `style.save_fig` figure-file metadata (PDF/SVG/PNG). Realizes "the exact code +
  environment that produced it, recoverable months later." New: `config.SKILL_VERSION`, `cfg.description`.
- **MAD univariate outlier flag** — `verify.flag_outliers_mad`: robust median/MAD flag for replicate
  band-metrics (the numpy-only univariate analog of `chemometrics.diagnostics`; scverse-style, our
  flag-don't-hide framing). New: `cfg.outlier_mad_n`.

**★ Marquee for later — a "guarded exploratory data analysis" agent (idea logged 2026-07-03).**
Reading the launch, Claude Science reads like "send an agent on a data fishing trip." The version that
fits THIS skill is **guarded exploration**, not free fishing (an unconstrained EDA agent on 18 spectra /
3 reps is a multiple-comparisons machine — it will "find" noise; our gates exist to stop exactly that):
> profile the dataset → **propose** candidate arguments/findings → each candidate is adversarially checked
> by a skeptic *and* must clear the existing statistical gates (lack-of-fit, permutation null, multiplicity
> correction, "n small ⇒ CIs wide") → only survivors become verified figures.
That's Claude Science's actor-critic applied to EDA, and it composes with two things already here: the
`figure_selection.md` decision module ("same data, different argument → different figure") and the
roadmap's proposed **batch profiler** (emit *facts* about levels×reps). The fishing trip is the profiler's
big sibling with a critic attached — and the **Tier-3 number critic built above is its safety rail**. We
now have the Agent/Workflow tooling to implement the actor-critic with real subagents. Bigger project;
scope it as: (1) `profile_batch` facts emitter [medium code], (2) a candidate-finding + adversarial-verify
harness [larger], (3) survivors → the existing verified-figure pipeline. Keep the invariant: **exploration
proposes; the gates dispose.**

---

## ★ Method-development agent thought exercise (2026-07-03)

A follow-on to the Claude Science pass: role-played an ideal analytical-chemistry agent
through a method-development project turn-by-turn (an IND — quantify API + excipients +
identify all impurities; grounded on **acetaminophen 500 mg** + real excipients incl. sodium
metabisulphite; constraints FTIR+HPLC-UV unlimited / 1 NMR / no MS / $1000), mining each turn
for features that are **buildable in the skill**. Filter for "do now": *useful for cocrystal ID.*

**★ Built:**
- **`references/cocrystal_id.md`** — the cocrystal-ID decision module (cocrystal vs salt vs
  physical mixture vs polymorph; PXRD new-phase, the FTIR carboxylate/ΔpKa salt discriminator,
  the lab-capability gate). Sibling of `figure_selection.md`; routes into `crystal_pxrd`/
  `spectra` and names the roadmap Tier-3 items (ΔpKa, NNLS sum-of-parents, DD-SIMCA) it will
  route to once built. **Directly cocrystal-ID turn-one tooling.** Wired into SKILL.md.
- **Web-search-then-cite discipline** (SKILL.md "When NOT to over-reach") — external facts
  (prices, reference values, limits, band assignments) get looked up + cited with source/date,
  never recalled. Demonstrated live in the exercise (USP/Sigma/Phenomenex price pull).

**Parked as skill-features** (surfaced, buildable, but orphaned for now — revisit
if a stability/impurity job returns):
- **Capability→requirement matrix** (general form) — `cocrystal_id.md` is its first concrete
  instance; generalise if a second domain needs it.
- **Excipient→degradation-liability map** — cited lookup (Mg-stearate → hydrolysis+metal-ox;
  metabisulphite → antioxidant *but* sulphite-adduct + oxidant-scavenge + allergen; reducing
  sugars → Maillard). Feeds the next item.
- **Forced-degradation design generator** — sibling of `doe.py`: `{API, excipients+liabilities}`
  → stress panel + required control sample sets (API-alone / DP / placebo / provocation binary)
  + peak-origin truth table. Its value is forcing the control sets a novice forgets.
- **Cited procurement/budget ledger** — `{item, price, source, date, confidence, why}` = the
  report's projected-vs-actual budget appendix as a provenance-stamped object.

**Platform-vision primitives** (NOT skill-buildable — need a memory substrate/connectors/
dependency-graph/agent orchestration) → moved to **`VISION.md`**: problem-model+decomposition,
constraint admissibility gate, knowledge connection, and the ★ epistemic layer (three kinds of
not-knowing + value-of-information questioning), wrapped in a legible-traceable problem-state.

---

## ★ Open-database connectors (2026-07-04)

Prompted by Claude Science's 60-connector model — but analytical chemistry's reference data is mostly
closed, so only two corners are open + API-clean, and both are now wired: **`scripts/sources.py`** —
`cod_search`/`cod_fetch_cif` (Crystallography Open Database → CIF → `crystal_engine`; **live-verified
end-to-end**: search aspirin → 23 hits → fetch → load) and `pubchem_properties` (PUG-REST; SMILES
omitted from defaults after the 2025 rename, pKa flagged as not-in-PubChem). Pure URL/parse logic is
gate-tested (`tests/test_sources.py`, 5 checks); network wrappers are best-effort + cite their URL.
The full sourcing map (incl. paywalled CSD/ICDD/pharmacopoeia + free-but-manual DataWarrior-pKa/
AqSolDB/SDBS/NIST/SpectraBase) is **`references/databases.md`**. Key finding: **crystallography is
the open corner** (COD + free Rietveld = the open ICDD substitute, which is what the crystal family
already implements); spectral libraries + monographs stay closed. Both live-verified (PubChem aspirin →
C9H8O4/180.16/XLogP 1.2). The first PubChem attempts failed with `getaddrinfo failed` — NOT a block:
its hostname is Google-Cloud-fronted (34.107.x / 2600:1901::) via a CNAME chain whose DNS
transiently flaked from this sandbox's resolver (`www.ncbi`, a direct NCBI IP, never did), so `_get`
gained a retry — which resolved it.

## Feature backlog (from the research pass)

Effort key: **[fn]** small (function/param) · **[sub]** new sub-family · **[heavy]** needs a
heavy/uncommon dep. All cross-checked as NOT already in the skill.

### Tier 1 — strengthens FTIR calibration work

Targets the documented wall of a small ATR calibration: multiplicative loading/coverage variation, "flag-don't-hide"
outliers, and ranking no-IS < matrix-IS < SNV-PLS.

**Validation-first:** a full **chemometrics validation suite** lives at
`skill_validation/chemometrics/` (spec in `SPEC.md`; run `python …/run.py`) — reference numbers
+ edge cases for the WHOLE family, R-free at test time, tiered tolerances, provenance for the
report appendix. **Existing-surface harness BUILT & green (2026-06-23): 33 checks, 0 fail** —
already caught a `pca_fit` n_components-cap bug and locked the SG-passthrough WARN + ddof=0 note.
Remaining: implement each Tier-1 feature below with its validation module green (targets already
prepared in `golden/`).

- **Outlier/diagnostic gate sub-family [sub, highest payoff] — ★ DONE & validated (2026-06-23).**
  `chemometrics.diagnostics` / `plot_influence`: Hotelling-T² + Q/SPE (Jackson–Mudholkar + Box) +
  leverage + **DD-SIMCA** SD–OD acceptance plot, with critical limits (PCA-based, numpy+scipy only
  → runs as raw-spectrum QC before calibration). Separates "bad load/weird spectrum" (high
  orthogonal distance) from "extreme but valid conc" (high leverage) — turns ad-hoc outlier
  calls (an under-loaded or a dilute replicate) into a principled, reportable gate. Validated bit-exact
  vs the closed-form fixture + scipy + published R `mdatools` limits (`test_outliers.py`).
- **EMSC preprocessor with a pure-excipient interferent term — ★ DONE & validated (2026-06-23).**
  `chemometrics.emsc`: jointly-fitted baseline + multiplicative term + optional interferent (pass
  the pure-excipient spectrum to model the matrix away). Bit-exact vs the biospectools/Kohler gold.
- **Validation trio — ★ DONE & validated (2026-06-23).** `permutation_test` (between-level y-scramble
  null for R²(CV)); `corrected_paired_t` (Nadeau–Bengio — makes "SNV-PLS 2.9% vs matrix-IS 5.8%" a
  defensible claim; golden reproduced); `rpd`/`rpiq` (the "is it useful" FoM).
- **Seeded run-order schedule generator — ★ DONE & validated (2026-06-23).** `scripts/doe.py`:
  `run_order_schedule`/`write_schedule_csv`/`balance_report` — operator orthogonal to concentration
  (rotating round-robin) + seeded run order (reproducible/archivable). Self-validates (determinism +
  balance + coverage). **→ Tier 1 complete.**
- **Smaller [fn]:** Kennard–Stone / SPXY deterministic splits (defensible for 18 spectra);
  Hotelling-T² confidence ellipse on the scores plot.

### Tier 2 — packaging / discipline — ★ COMPLETE & green (2026-06-23)

`python -m pytest tests/` = one green gate (29 tests, bridges in the 85-check numeric suite + doe).

- **pytest skill self-tests — DONE.** `tests/` asserts module imports, SKILL.md structure,
  `bundle.py` round-trip (no skill imports survive, compiles), multiplicity vs statsmodels,
  `.mplstyle`/width, the CLI — and runs the numeric suites as subprocesses. Catches drift; citable.
- **Raw-vs-adjusted p-value (multiplicity) gate — DONE.** `verify.adjust_pvalues` /
  `multiplicity_check` (Bonferroni/Holm/BH) — validated bit-for-bit vs `statsmodels.multipletests`.
- **`.mplstyle` journal presets + `check_figure_width` — DONE.** `style.export_mplstyle` +
  `check_figure_width`; presets shipped in `assets/styles/` (nature/acs/ieee/general). House
  rcParams refactored to a single source (`_main_rc`/`_BASE_RC`).
- **Scripts-first CLI — DONE.** `cli.py`: `spc-rsd` (.spc → band area/%RSD) + `calibrate`
  (csv → slope/R²/LOD/LOQ). Framed as QC convenience; the report deliverable stays a written analysis.

### Tier 3 — cocrystal identification

The crystal family originally had no *identification/classification* tooling — which is what cocrystal ID is.

- **Classification sub-family [sub] — ★ DONE (one-class + FoM) & validated (2026-07-03).**
  `chemometrics.ddsimca_fit`/`ddsimca_predict` (one-class; *rejects* an unknown → the right frame for
  "is this the target cocrystal?", reusing the validated DD-SIMCA limits with out-of-sample projection)
  + `chemometrics.class_fom` (sensitivity/specificity/efficiency — the correct vocabulary, not R²/LOD).
  Validated in `skill_validation/chemometrics/test_classify.py` (13 checks: closed-form FoM, in-sample
  Type-I acceptance, out-of-class rejection, class-separating boundary, predict-vs-diagnostics
  consistency). **Still TODO: PLS-DA** multi-class (cocrystal vs physical mixture vs starting materials).
  **→ Wave 1 (ΔpKa · NNLS · one-class classify) COMPLETE.**
- **ΔpKa salt-vs-cocrystal calculator [fn, zero deps] — ★ DONE & validated (2026-07-03).**
  `cocrystal.classify_ionisation` / `delta_pka`: the Cruz-Cabeza 3-zone rule (>3 salt, <0 cocrystal,
  else continuum; configurable thresholds), never raises, cited caption. Orthogonal evidence to "new
  phase"; aspirin pKa≈3.5, so the coformer's basicity decides salt vs cocrystal. Validated closed-form
  in **`skill_validation/cocrystal/`** (13 checks: arithmetic, zones, boundaries, config, regimes),
  bridged into `pytest tests/`. **Wave-1 next: NNLS sum-of-parents → class-FoM (below).**
- **NNLS sum-of-parents fit + Rwp residual + new/lost-peak table [fn, scipy only] — ★ DONE & validated
  (2026-07-03).** `cocrystal.sum_of_parents`/`rwp`/`unexplained_peaks`/`phase_report`: fits the product
  as a non-negative combination of parent patterns (scipy NNLS), reports Rwp (poisson/unit) + new/lost
  peaks + a heuristic verdict — turns "new peaks appear" into a defensible number distinguishing
  cocrystal from physical mixture (PXRD *or* FTIR). Validated in `skill_validation/cocrystal/` (13
  checks: closed-form Rwp, fraction recovery on a synthetic mixture, new-phase discrimination,
  `scipy.optimize.nnls` parity, new-peak detection). **Wave-1 next: class-FoM + DD-SIMCA classify.**
- **chmpy Hirshfeld surfaces + decomposed fingerprints [sub, pure-Python pip]** — THE standard figure
  for how contacts change on cocrystallization (new O–H···N heterosynthon appearing). Rare in being
  dependency-light.
- **PXRD realism [fn] — ★ DONE & validated (2026-07-03; wave 2).** `scripts/pxrd_realism.py`:
  `march_dollase` (preferred orientation) + `caglioti_fwhm`/`pseudo_voigt` (angle-dependent
  broadening) + **Cu Kα1/Kα2 doublet** (`kalpha2_doublet` — no lib does the doublet automatically),
  composed by `simulate_pattern`. Numpy-only on a reflection list (Dans supplies structure factors),
  so the physics is testable without heavy deps. New cfg: `pxrd_kalpha2`/`pxrd_wavelength2`/
  `pxrd_kalpha2_ratio`, `pxrd_po_axis`/`pxrd_march_r`, `pxrd_caglioti`. Validated in
  `skill_validation/pxrd/` (23 checks). Doc: `crystal.md` calc-PXRD section. **Seam CLOSED (2026-07-04):**
  `crystal_pxrd.reflection_list` (Dans powder reflections — |F|²·multiplicity·LP) + `realistic_pattern`
  turn-key, validated end-to-end on real CIFs (aspirin, lactose) in
  `skill_validation/crystal/test_realism_cif.py` (calc positions reproduced; Kα2 companion appears; render
  eyeballed). Found + documented a coherence gotcha: the Kα2 default is Cu, so pin α1 to Cu (set
  `cfg.pxrd_wavelength`) when simulating a Cu scan from a Mo-refined CIF.
- **CIF validation depth [fn]** — ADP sanity panel (non-positive-definite ellipsoids, Hirshfeld
  rigid-bond) + **ADDSYM missed-symmetry sweep** via spglib (a wrong space group → wrong PXRD
  fingerprint).
- **RDKit conformer-overlay RMSD [sub-ish]** — API conformation in pure-API vs cocrystal; a change at
  aspirin's acetyl ester torsion is the structural correlate of the ester C=O FTIR shift (ties to the
  FTIR through-line).
- **Wildcard [sub]** — Indirect Hard Modeling / MCR-ALS (`chemometrics.IHM`, `McrAR`): fit overlapping
  bands as shifting/broadening pure-component shapes. Uniquely relevant — *peak shift IS the cocrystal
  signal*.

### Explicitly skip / cite-only

GA/CARS/MCUVE variable selection and full PDS/DS calibration transfer (disproportionate for a binary
mixture; slope-bias is enough for the "instrument transferability" extended task); **OPLS** (advertised
but its API isn't actually exposed — don't build on it); Rietveld/GSAS-II, PLATON SQUEEZE binary,
anything **CCDC** (COMPACK, Etter graph-sets, true BFDH, energy frameworks), **occpy** (no Windows
wheel). Cite these as references; don't embed.

---

## From the pre-publication review (2026-09-11)

Static review of `scripts/`, `bundle.py`, `cli.py` before the repository went public. Fixed at
once: `cli.py --help` crash (a bare `%` in an argparse help string), the chemometrics suite
crashing instead of skipping without scikit-learn, silent fallbacks in `correct_baseline` / `_sg` /
the t-multiplier (now they WARN and the z quantile is exact), the FTIR dead-pixel heuristic bypassing
the whole ingest gate for PXRD, `rwp`'s Poisson default on profiles with true zeros (now `auto`),
`simulate_pattern`'s hard-coded Cu α1 and two silent no-ops, an undeclared `escape_tol_frac` knob.
Deferred (real, larger):

- **One H-bond / bond predicate.** `crystal_view._hbond_pairs` / `_hbond_env_atoms` and the bond
  predicate copied across `crystal_engine.bonds` / `geometry` / `crystal_view.complete_molecules` /
  `_bond_pairs` / `_components` re-implement the engine's criteria and omit its disorder-alternative
  exclusion → a figure can draw a contact the table suppresses. Expose `crystal_engine.is_hbond` /
  `is_bonded` and call them everywhere (the module's stated invariant).
- **Vectorise the crystal geometry.** All-pairs Python loops with `np.linalg.norm` per pair and two
  gemmi `Element` constructions per pair; cart coords + the 27-image supercell rebuilt per helper.
  `cKDTree.query_pairs` / `query_ball_point` on a cached supercell, a covalent-radius dict, and
  `_hbond_pairs` computed once per render → seconds to tens of milliseconds.
- **`simulate_pattern` per-peak windows** (`searchsorted` ±10 FWHM) instead of a full-grid
  pseudo-Voigt per reflection; matters once cell refinement loops call it hundreds of times.
- **Cache the Dans structure / reflection list on `Structure`** so `calc_pattern`, `peak_table`,
  `realistic_pattern` and `plot_overlay` stop re-parsing the CIF and recomputing identical
  structure factors; run the pymatgen cross-check once per (cif, λ, window).
- **Nested-component PLS scan**: one fit at the cap per fold, predictions for a = 1..cap by
  truncation; cache per-fold preprocessed X across permutations in `permutation_test`.
- **Lazy pyplot** in `style` so `cli.py` and CSV-only standalones don't pay the matplotlib import.
- **Stop probing declared cfg fields with `getattr(cfg, name, default)`** — a dataclass field is
  safe to access directly and a misspelt field should fail loudly.

## Suggested sequence

1. Write `references/figure_selection.md` (the decision module) — pure writing, immediate consistency win.
2. Tier 1 chemometrics adds, in order: outlier/diagnostic gate → EMSC → validation trio → seeded schedule.
3. Tier 2 packaging when convenient (good before report submission).
4. Tier 3 for cocrystal ID; ΔpKa + NNLS-sum-of-parents are the cheapest first wins.

## Dependency notes

Tier-1 chemometrics adds ride on the existing lazy-sklearn path. Crystal Tier-3 deps (spglib, chmpy,
RDKit, pymatgen, PyXtal) are pip-installable, pure-Python-friendly, lazy-load alongside gemmi/
Dans_Diffraction. Genuinely blocked on this Windows setup: GSAS-II, the PLATON binary, occpy, anything CCDC.

## Sources

scipilot-figure-skill (https://github.com/Haojae/scipilot-figure-skill) · `chemometrics`
(readthedocs) · PyChemAuth (mahynski) · pyChemometrics (clicumu) · NIRPY Research · pybaselines ·
PyXtal · chmpy (peterspackman) · spglib · RDKit · K-Dense-AI scientific-agent-skills · scicraft
(jaechang-hits).
