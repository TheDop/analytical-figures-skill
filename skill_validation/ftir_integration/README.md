# analytical-figures — FTIR peak-integration validation

Goal: confirm the `analytical-figures` skill reproduces **real, published** FTIR
quantitative-analysis-by-peak-integration results.

Run: `python skill_validation/ftir_integration/validate_ftir_integration.py` (parts A1 and B
need nothing external; A2 needs the Mendeley raw spectra, see the table) and
`python skill_validation/ftir_integration/groundtruth_integration_test.py` (deterministic,
bridged into `pytest tests/`). Figures in `out/`.

## Datasets

| | Source | Why chosen |
|---|---|---|
| **A. Pectin DM** | Wang 2023, Mendeley Data `doi:10.17632/gkwbp3wc49.1`, CC BY 4.0 (**not vendored**: download `Raw data.zip` from the record and unzip into `pectin_mendeley/unzipped/`) | The **only open source with RAW spectra + areas + calibration**, so it validates the *integration* step. DM read from a band-area **ratio** `I = A_ester(~1745)/(A_ester+A_carboxylate(~1605))` — the internal-standard-ratio idea used for API/excipient calibrations. |
| **B. Doxorubicin / arterolane** | Bansal, Singh, Kaur, *BMC Chemistry* 15:27 (2021), `doi:10.1186/s13065-021-00752-3` | Pharma; publishes per-point peak-**area** tables (carbonyl band, baseline-corrected) for two modes. DOX C=O (1770–1680) is the analogue of an API ester C=O. |
| **C. Diclofenac sodium** (attempted) | Fahelelbom et al., *F1000Research* (2020) `doi:10.12688/f1000research.22274.2`; raw data Harvard Dataverse `doi:10.7910/DVN/6SJZ7W` | Procedure fully stated (first-derivative AUC, 1550–1605 cm⁻¹) **and** raw spectra deposited — structurally the "unicorn". **But the deposit is corrupted:** the 1.0% standard column is byte-identical to 0.8%, spectra are decimated to ~7.7 cm⁻¹ (~8 pts/band), and 0.6% is out of order. No AUC definition reproduces their Y=1.375X−0.014 / R²=0.9994 → **not reproducible**. |

## Results

**Calibration engine — reproduces exactly.**
- Pectin (their published areas → `I` → fit): **R² = 0.99352**; the regression
  reproduces all of the authors' worked-unknown DM values to **< 2×10⁻⁵**.
- BMC 2021, all four calibrations: slopes/intercepts match to the published
  rounding (e.g. DOX-trans **3.002x+0.402** vs **3x+0.402**; ALM-refl
  **1.379x−0.066** vs **1.38x−0.068**).

**LOD/LOQ — matches only when the paper used the same σ.**
DOX-transmittance matches (0.109/0.330 vs 0.11/0.30); the other three BMC bands do
not, because the paper used an undocumented band-specific σ. The skill's
`lod_loq` always *names* its σ basis (residual SD) — which is the point.

**Integration step — exposed a real limitation, now fixed.**
`integrate_bands` originally drew a **separate local baseline per window**. The
pectin ester/carboxylate bands overlap on a shared pedestal, so the per-window
baseline followed the valley *up* and carved area off the smaller band (DM3 ester
area came out ~5× too small) → calibration R² collapsed to **0.77**.
Added `cfg.integration_baseline = "shared"`: **one** baseline across the whole
envelope with a vertical-drop split at the window boundary (the standard
"drop-perpendicular" method). That removes the artefact and lifts the
re-integrated calibration to **R² 0.87** (DM3 ΔI −0.142 → −0.041).
Residual gap to their 0.9935 is the authors' exact baseline anchors + valley
split, which the dataset specifies only graphically — not a code defect.

## Deterministic ground-truth test (the close-out)

Three real-world sources, three *different* blockers (A: graphical-only limits; B:
no raw spectra; C: corrupted deposit) — so none pins the integrator to an *exact*
known area. `groundtruth_integration_test.py` removes that ambiguity: it builds
spectra from analytic peaks (known closed-form areas) under real-world conditions
(sloping+offset baseline, ~2 cm⁻¹ FTIR spacing, additive noise, random ATR loading,
overlapping bands) and checks against independent references (closed form +
`scipy.integrate.quad`, not the skill's trapezoid). **14/14 checks pass:**

- **T1** isolated band recovers its analytic Gaussian area through a sloping baseline (rel. err 2×10⁻⁵); noise unbiased (0.19% over n=200).
- **T2** both `shared` and `per_window` match an independent integrator of their respective definitions to ~10⁻⁵ (code is correct). With a real pedestal (valley 23% of peak), shared-mode ratio is off truth by −0.0009 vs per_window's +0.0136 → shared ~15× closer.
- **T3** the band-area ratio is invariant to overall scale/ATR loading (2×10⁻¹⁶) and absolute area scales exactly linearly — the internal-standard cancellation, proven.
- **T4** realistic DM-like series: `shared` gives a near-ideal linear calibration (slope 0.965, R²=0.99998); `per_window` is slope/intercept-biased (1.258, −0.126) and the modes diverge by up to 0.106 in ratio (`out/groundtruth_shared_vs_perwindow.png`). The pectin finding, reproduced against a known answer.

## Implication

The calibration/figures-of-merit engine is trustworthy. For any **band-area
ratio** (an API band vs an excipient reference band), if the two bands do not
return to baseline between them, use `integration_baseline="shared"` — a
per-window baseline biases the ratio and flattens the calibration. And quote LOD/
LOQ with the σ method named.
