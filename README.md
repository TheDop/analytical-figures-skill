# analytical-figures

A [Claude Code](https://claude.com/claude-code) skill for **publication-grade figures and the
gated analysis behind them** in analytical chemistry and spectroscopy: FTIR/ATR and PXRD spectra,
calibration curves with figures of merit (LOD/LOQ, recovery, ICH Q2), multivariate calibration
(PLS/PCR/PCA), crystal structures from CIFs (validation table, calculated PXRD, deterministic 3D
views) and cocrystal identification.

It is not a plotting tutorial. It does two things a plotting library does not:

- **It gates the numbers before they reach a figure.** Ingest checks, identical processing across
  a batch, named error statistics (SD vs SEM vs CI, with n), leakage-safe cross-validation,
  mandatory residual panels, no extrapolation past the calibrated range, and a static critic
  that flags a figure-of-merit hand-typed onto a figure instead of interpolated from the code.
- **It hands over a script, not a PNG.** The deliverable is one self-contained Python file
  (`bundle.py` amalgamates the modules an analysis uses) with a provenance header, so every
  figure is reproducible, editable and attachable to a report.

## Install

Clone into a Claude Code skills directory, user-level or per project:

```
git clone https://github.com/TheDop/analytical-figures-skill ~/.claude/skills/analytical-figures        # available in every project
git clone https://github.com/TheDop/analytical-figures-skill .claude/skills/analytical-figures           # this project only
python -m pip install -r ~/.claude/skills/analytical-figures/requirements.txt
```

`numpy` and `matplotlib` are required; `scipy` is strongly recommended. Everything else is
per-family and imported lazily (`scikit-learn` for chemometrics, `gemmi` + `Dans_Diffraction`
for crystal structures, `rdkit` for the cocrystal predictor), so an FTIR job never pulls them.
Claude Code loads the skill at session start. Then ask for what you need: a replicate overlay
with %RSD, a calibration curve with LOD/LOQ, a calculated PXRD pattern from a CIF, a
cocrystal-vs-mixture call.

## What is inside

| Path | Role |
|---|---|
| `SKILL.md` | The operating instructions Claude reads: workflow, hard principles, the single-CONFIG rule. |
| `scripts/` | The maintained modules: `config`, `style`, `verify`, `spectra`, `calibration`, `chemometrics`, `charts`, `doe`, `report`, `crystal_engine`, `crystal_pxrd`, `crystal_view`, `pxrd_realism`, `cocrystal`, `sources`. |
| `references/` | Decision docs: `figure_selection.md` (what to plot, and why), `cocrystal_id.md` (which evidence identifies a cocrystal), `databases.md` (where to source each external fact), plus per-family conventions. |
| `templates/analysis_template.py` | Copy this to start an analysis; `bundle.py` turns it into a standalone deliverable. |
| `cli.py` | QC shortcuts: `.spc` replicates to band area/%RSD, a conc/signal CSV to slope, R², LOD/LOQ. |
| `assets/styles/` | Journal `.mplstyle` presets (Nature, ACS, IEEE, general). |
| `tests/`, `skill_validation/` | The test gate and the numeric-parity suites, below. |
| `ROADMAP.md`, `VISION.md` | Backlog and build log; longer-range ideas. |

The modules are plain Python, so the skill also works without Claude: run from the skill root
and `from scripts import spectra, calibration` as the template does.

## Validation

`python -m pytest tests/` is one green gate for the whole skill (79 tests). It checks module
imports, the `SKILL.md` structure, the `bundle.py` round-trip, p-value correction against
`statsmodels`, the style presets and the CLI, and it bridges in three numeric-parity suites:

| Suite | Checks | Validated against |
|---|---|---|
| `skill_validation/chemometrics/` | 100 | closed-form fixtures, `scikit-learn`, published R `mdatools` critical limits, the biospectools EMSC gold |
| `skill_validation/cocrystal/` | 26 | closed-form Rwp, `scipy.optimize.nnls`, the Cruz-Cabeza ΔpKa zones |
| `skill_validation/pxrd/` | 23 | March–Dollase, Cu Kα₁/Kα₂ doublet, Caglioti, pseudo-Voigt |

Two external datasets are fetched on demand rather than vendored (`skill_validation/chemometrics/datasets/fetch.py`);
their checks skip cleanly when absent. The crystal family was additionally exercised on a CIF
edge-case set (disorder, special positions, Z′ > 1, non-tabulated elements, deliberately broken
inputs) and the FTIR integrator against published calibrations; see `references/crystal.md`
and `references/verification.md`.

## Provenance and licence

Built during an MSc analytical-chemistry project: quantifying an API in a binary excipient
mixture by ATR-FTIR, then identifying cocrystals by FTIR and PXRD. The worked examples in the
docs (aspirin in lactose, theophylline cocrystals) come from that work.

MIT licence (see `LICENSE`). The figure-selection scaffolding and the panel-label alignment
trick are adapted from [scipilot-figure-skill](https://github.com/Haojae/scipilot-figure-skill)
(MIT, Haojae); details in `THIRD_PARTY_NOTICES.md`.
