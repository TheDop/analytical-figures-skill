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

<!-- gallery:start -->
<table>
<tr><td align="center"><a href="docs/gallery/full/calib_pectin.png"><img src="docs/gallery/tiles/calib_pectin.png" width="140" alt="Calibration with confidence and prediction bands and the mandatory residual panel" title="Calibration with confidence and prediction bands and the mandatory residual panel — data: Pectin degree of methyl-esterification by FTIR — Wang 2023, Mendeley Data 10.17632/gkwbp3wc49.1 (CC BY 4.0)"></a></td><td align="center"><a href="docs/gallery/full/waterfall.png"><img src="docs/gallery/tiles/waterfall.png" width="140" alt="Replicate waterfall: reps overlaid per level, levels offset, right-margin keys instead of a legend" title="Replicate waterfall: reps overlaid per level, levels offset, right-margin keys instead of a legend — data: Synthetic binary mixture, 6 levels × 3 reps with multiplicative loading variation (skill_validation/chemometrics/_data.py)"></a></td><td align="center"><a href="docs/gallery/full/pls_diag.png"><img src="docs/gallery/tiles/pls_diag.png" width="140" alt="Leakage-safe PLS diagnostics: RMSECV per preprocessing with the parsimony pick, coefficients, scores, loadings" title="Leakage-safe PLS diagnostics: RMSECV per preprocessing with the parsimony pick, coefficients, scores, loadings — data: Synthetic binary mixture (skill_validation/chemometrics/_data.py), leave-one-level-out CV"></a></td><td align="center"><a href="docs/gallery/full/pxrd_overlay.png"><img src="docs/gallery/tiles/pxrd_overlay.png" width="140" alt="Stacked calculated PXRD patterns for phase discrimination" title="Stacked calculated PXRD patterns for phase discrimination — data: Aspirin COD 7247819, urea COD 1008785, lactose COD 2206486 — Crystallography Open Database"></a></td><td align="center"><a href="docs/gallery/full/ortep.png"><img src="docs/gallery/tiles/ortep.png" width="140" alt="ORTEP displacement ellipsoids at 50 % probability, deterministic PCA camera" title="ORTEP displacement ellipsoids at 50 % probability, deterministic PCA camera — data: Aspirin, COD 7247819 (296 K)"></a></td><td align="center"><a href="docs/gallery/full/ddsimca.png"><img src="docs/gallery/tiles/ddsimca.png" width="140" alt="DD-SIMCA acceptance plot: regular / extreme / outlier, with the chi-squared boundary" title="DD-SIMCA acceptance plot: regular / extreme / outlier, with the chi-squared boundary — data: Synthetic binary mixture plus one under-loaded and one baseline-tilted replicate injected"></a></td></tr>
<tr><td align="center"><a href="docs/gallery/full/integration.png"><img src="docs/gallery/tiles/integration.png" width="140" alt="Band area on algorithmic anchors: the flanking minima found once and locked for the whole batch" title="Band area on algorithmic anchors: the flanking minima found once and locked for the whole batch — data: Synthetic binary mixture (skill_validation/chemometrics/_data.py)"></a></td><td align="center"><a href="docs/gallery/full/calib_bmc.png"><img src="docs/gallery/tiles/calib_bmc.png" width="140" alt="Doxorubicin carbonyl-area calibration; LOD/LOQ from the residual SD, and the sigma is named" title="Doxorubicin carbonyl-area calibration; LOD/LOQ from the residual SD, and the sigma is named — data: Bansal, Singh &amp; Kaur 2021, BMC Chemistry 15:27, 10.1186/s13065-021-00752-3 (published area table)"></a></td><td align="center"><a href="docs/gallery/full/groundtruth.png"><img src="docs/gallery/tiles/groundtruth.png" width="140" alt="Shared vs per-window baseline against a closed-form answer: the per-window pedestal bias, measured" title="Shared vs per-window baseline against a closed-form answer: the per-window pedestal bias, measured — data: Analytic Gaussian bands with known areas, checked against scipy.integrate.quad"></a></td><td align="center"><a href="docs/gallery/full/pxrd_calc.png"><img src="docs/gallery/tiles/pxrd_calc.png" width="140" alt="Calculated PXRD from a CIF at Cu Kα, top peak cross-checked against pymatgen" title="Calculated PXRD from a CIF at Cu Kα, top peak cross-checked against pymatgen — data: Aspirin, COD 7247819"></a></td><td align="center"><a href="docs/gallery/full/pxrd_realism.png"><img src="docs/gallery/tiles/pxrd_realism.png" width="140" alt="Lab-realistic pattern: Cu Kα₁/Kα₂ doublet and Caglioti broadening on the Dans reflection list" title="Lab-realistic pattern: Cu Kα₁/Kα₂ doublet and Caglioti broadening on the Dans reflection list — data: Aspirin, COD 7247819"></a></td><td align="center"><a href="docs/gallery/full/packing.png"><img src="docs/gallery/tiles/packing.png" width="140" alt="Packing diagram down a with hydrogen bonds and the cell box" title="Packing diagram down a with hydrogen bonds and the cell box — data: Aspirin, COD 7247819"></a></td></tr>
<tr><td align="center"><a href="docs/gallery/full/hbond_env.png"><img src="docs/gallery/tiles/hbond_env.png" width="140" alt="Hydrogen-bond environment of the asymmetric unit, symmetry mates superscripted" title="Hydrogen-bond environment of the asymmetric unit, symmetry mates superscripted — data: Aspirin, COD 7247819"></a></td><td align="center"><a href="docs/gallery/full/unit_cell.png"><img src="docs/gallery/tiles/unit_cell.png" width="140" alt="Unit-cell contents with the cell box and a/b/c" title="Unit-cell contents with the cell box and a/b/c — data: Aspirin, COD 7247819"></a></td><td align="center"><a href="docs/gallery/full/flufenamic.png"><img src="docs/gallery/tiles/flufenamic.png" width="140" alt="Z′ = 3 with CF₃ disorder groups honoured: no bonds between mutually exclusive sites" title="Z′ = 3 with CF₃ disorder groups honoured: no bonds between mutually exclusive sites — data: Flufenamic acid, COD 4118081"></a></td><td align="center"><a href="docs/gallery/full/pxrd_urea.png"><img src="docs/gallery/tiles/pxrd_urea.png" width="140" alt="Calculated PXRD of a molecule on a special position: multiplicity handled in the expansion" title="Calculated PXRD of a molecule on a special position: multiplicity handled in the expansion — data: Urea, COD 1008785"></a></td><td align="center"><a href="docs/gallery/full/ballstick.png"><img src="docs/gallery/tiles/ballstick.png" width="140" alt="Ball-and-stick, PCA face-on, C–H hidden, heteroatoms labelled" title="Ball-and-stick, PCA face-on, C–H hidden, heteroatoms labelled — data: Aspirin, COD 7247819"></a></td><td align="center"><a href="docs/gallery/full/pyvista.png"><img src="docs/gallery/tiles/pyvista.png" width="140" alt="PyVista/VTK still: real depth buffer, ambient occlusion, orthographic camera from the same orientation engine" title="PyVista/VTK still: real depth buffer, ambient occlusion, orthographic camera from the same orientation engine — data: Aspirin, COD 7247819"></a></td></tr>
</table>

<sub>Hover a tile for what it is and where the data came from; click for full size. Every tile is produced by <code>docs/build_gallery.py</code> from public data, public CIFs or the synthetic validation set — nothing is drawn by hand. Interactive version with placards: <a href="https://thedop.github.io/analytical-figures-skill/">https://thedop.github.io/analytical-figures-skill/</a></sub>
<!-- gallery:end -->

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

`python -m pytest tests/` is one green gate for the whole skill (81 tests). It checks module
imports, the `SKILL.md` structure, the `bundle.py` round-trip, p-value correction against
`statsmodels`, the style presets and the CLI, and it bridges in three numeric-parity suites:

| Suite | Checks | Validated against |
|---|---|---|
| `skill_validation/chemometrics/` | 100 | closed-form fixtures, `scikit-learn`, published R `mdatools` critical limits, the biospectools EMSC gold |
| `skill_validation/cocrystal/` | 26 | closed-form Rwp, `scipy.optimize.nnls`, the Cruz-Cabeza ΔpKa zones |
| `skill_validation/pxrd/` | 23 | March–Dollase, Cu Kα₁/Kα₂ doublet, Caglioti, pseudo-Voigt |

| `skill_validation/ftir_integration/` | 14 | closed-form Gaussian areas and `scipy.integrate.quad` under sloping baselines, noise and band overlap; plus the integrator run against the published pectin-DM (Wang 2023) and BMC 2021 calibrations |
| `skill_validation/crystal/` | 9 CIFs | a CIF edge-case set (disorder groups, special positions, Z′ = 9, chiral, multi-block, synthetic broken inputs): the engine's table, an independent gemmi-only second opinion, ADP U_eq, the Dans reflection list vs `pymatgen`. Needs `gemmi` + `Dans_Diffraction`; pytest skips it otherwise |

External data is fetched on demand rather than vendored: the biospectools EMSC gold and the R
`pls` gasoline set (`skill_validation/chemometrics/datasets/fetch.py`), the Mendeley pectin
spectra, and the CCDC CIFs, which cannot be redistributed. Every check that needs one of them
skips cleanly when it is absent.

## Provenance and licence

Built during an MSc analytical-chemistry project: quantifying an API in a binary excipient
mixture by ATR-FTIR, then identifying cocrystals by FTIR and PXRD. The worked examples in the
docs (aspirin in lactose, theophylline cocrystals) come from that work.

MIT licence (see `LICENSE`). The figure-selection scaffolding and the panel-label alignment
trick are adapted from [scipilot-figure-skill](https://github.com/Haojae/scipilot-figure-skill)
(MIT, Haojae); details in `THIRD_PARTY_NOTICES.md`.
