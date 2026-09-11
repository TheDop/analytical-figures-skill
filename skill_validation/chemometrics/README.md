# chemometrics validation suite

Numeric-parity + edge-case regression tests for the `chemometrics` family. Reproduces
**exact numbers** from published worked examples / open reference implementations, and is
**R-free at test time** (R numbers are captured as literals in `golden/`). Design + the full
target list: **`SPEC.md`**.

## Run

```
python skill_validation/chemometrics/run.py
```

Exit 0 if all PASS (SKIP is not a failure); exit 1 on any FAIL/ERROR. Zero-dep runner — no
pytest. Each line shows the tolerance **tier** and the **provenance** of the reference number.

## Status (85 PASS, 0 FAIL, 1 SKIP)

- **Existing surface — green:** preprocessing (SNV/MSC/SG d1·d2/centre/autoscale, incl. chemotools
  parity), PLS (engine identity), PCR/PCA, leakage-safe CV, VIP, edge cases.
- **Outlier / diagnostic gate — green** (`test_outliers.py`): T²/Q/leverage + DD-SIMCA, validated
  bit-for-bit against the hand-derived `golden/closed_form.py` fixture + scipy recomputes, and
  against the published R `mdatools` Hotelling limits. SD–OD acceptance plot (linear + sqrt) renders.
- **EMSC — green** (`test_emsc.py`): bit-for-bit vs the biospectools/Kohler gold (`datasets/
  emsc_testdata.xlsx`), plus closed-form b-recovery + interferent removal.
- **Validation trio — green** (`test_trio.py`): permutation p-formula, Nadeau–Bengio corrected
  paired t (golden reproduced through the function), RPD/RPIQ.
- **1 SKIP:** the gasoline T2 check (needs `datasets/gasoline.npz`, exported from R `pls`).

## Tolerance tiers

| Tier | Meaning | Assert |
|---|---|---|
| T0 | closed-form / hand-derived | ≤ 1e-12 |
| T1 | same-language reference (sklearn, scipy, chemotools) | ≤ 1e-8 |
| T2 | cross-language / published literal (R `pls`, `mdatools`, papers) | documented (e.g. ≤ 0.5–5 % rel) + a named reason for the gap |
| T3 | property / invariant | exact logical assertion |

## Optional extras (upgrade SKIPs to checks)

- **`chemotools`** (`pip install chemotools`) → SNV/MSC same-language parity (T1).
- **Datasets** → `python skill_validation/chemometrics/datasets/fetch.py` (see its docstring):
  `gasoline.npz` enables the R `pls` RMSEP regression (T2); `emsc_testdata.xlsx` enables the
  biospectools EMSC parity (T1, step 2). Tests SKIP cleanly when absent.

## What it found / locked (first run, 2026-06-23)

- **Bug fixed:** `pca_fit` did not cap `n_components` when one was passed explicitly (only when
  `None`) → an over-large request raised an uncaught sklearn `ValueError`. Now capped at
  `min(n_samples, n_features)`, consistent with `pls_fit`/`pcr_fit`.
- **§5 fixes applied:** SG derivative now WARNs when it silently passes through (window ≤ poly);
  the **ddof=0** convention for SNV/autoscale is documented in the `Preprocessor` docstring.

## Key provenance (see `golden/manifest.py` for the literals)

| Check | Reference | Tier |
|---|---|---|
| SG d1/d2 | `scipy.signal.savgol_filter` (the engine) | T0 |
| SNV/MSC/centre/autoscale | closed form; `StandardScaler`; (chemotools) | T0/T1 |
| PLS coef & yhat | direct `sklearn.PLSRegression` (engine identity) | T1 |
| PCA variance / PCR back-map | direct sklearn | T1 |
| RMSECV / SEP-RMSEP identity | IUPAC Gold Book; ASTM E1655 | T0 |
| leakage-safe CV | hand-rolled in-fold loop | T1 |
| VIP | Σ VIP² = p invariant (Chong & Jun 2005) | T0 |
| gasoline RMSEP @2c = 0.2966 | R `pls` vignette (Kalivas 1997) | T2 (skips w/o dataset) |

## Files

`run.py` · `_harness.py` (tiered assert helpers) · `_data.py` (deterministic synthetic spectra)
· `test_*.py` (one per area) · `golden/manifest.py` + `golden/closed_form.py` · `datasets/fetch.py`.
