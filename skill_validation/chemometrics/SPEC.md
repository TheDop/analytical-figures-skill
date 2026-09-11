# Chemometrics validation suite — specification

A standing, R-free numeric-validation harness for the `chemometrics` family (existing
surface **and** the Tier-1 additions). Mirrors how the crystal family was validated:
reproduce **exact numbers** from published worked examples / open reference
implementations, and **pre-plan edge cases**. Every passing test records provenance so the
suite doubles as reproducibility evidence for the report appendix.

Status: SPEC for review. Nothing implemented yet. Research basis: 2026-06-23 fan-out
(four agents, cross-checked against live sklearn 1.9 / scipy 1.17 and the actual module).

---

## 1. Philosophy

- **Numeric parity, tiered by what's honestly achievable** (see §3). Bit-exact where the
  reference is same-language or closed-form; documented tolerance where it's cross-language
  (R↔Python, NIPALS↔SIMPLS) and the gap has a *named reason*.
- **R-free at test time.** R `pls`/`mdatools` numbers are captured as **literals** (from
  vignettes/books/committed package unit tests) in a golden manifest with provenance. The
  suite never shells out to R.
- **Optional Python references degrade gracefully.** scipy/sklearn are always present;
  `chemotools`/`biospectools`/`pychemauth` tests **skip with a clear message** if the dep is
  absent (like `scienceplots` in `style.py`). A dependency-free closed-form core always runs.
- **Edge cases are first-class** — their own module per area, enumerated in §6/§8.

## 2. Layout

```
skill_validation/chemometrics/
  SPEC.md                 # this file
  README.md               # how to run; the provenance table
  run.py                  # zero-dep runner: discovers test_*, prints PASS/FAIL + matched value + tol + provenance
  golden/
    manifest.py           # golden numbers as literals, each {value, source, version, settings, tier}
    closed_form.py        # the §7(a) self-contained fixture + its hand-derived expected values
  datasets/
    fetch.py              # optional: download gasoline / biospectools xlsx (cached); tests skip if offline
    gasoline.npz          # vendored if license permits, else fetched
  test_preprocessing.py   test_pls.py  test_pcr_pca.py  test_cv.py  test_vip.py
  test_outliers.py        test_emsc.py  test_trio.py
  test_edge_cases.py      # consolidated edge matrix (§8)
```

`run.py` is plain Python (no pytest required) to match the skill's zero-dep ethos; tests are
functions returning `(name, tier, passed, got, want, tol, provenance)`. A `--pytest` shim is
trivial to add later.

## 3. Tolerance tiers (every test is tagged)

| Tier | Meaning | Assert |
|---|---|---|
| **T0** | closed-form / hand-derived | ≤ 1e-12 |
| **T1** | same-language reference (sklearn, scipy, chemotools, biospectools) | ≤ 1e-8 |
| **T2** | cross-language / published literal (R `pls`, `mdatools`, papers) | documented (e.g. ≤ 0.5–2 % rel) **+ reason for the gap** |
| **T3** | property / invariant (no external number) | exact logical assertion |

## 4. Conventions we PIN (ratify these — they decide what "exact" means)

1. **Preprocessing ddof = 0** for SNV/autoscale (matches `chemotools` + `StandardScaler`;
   R uses ddof=1, a √(n/(n−1)) factor). **Keep ddof=0; document it in SKILL.md.**
2. **PLS engine = sklearn NIPALS.** Golden coef/yhat = a direct `PLSRegression` fit (T1,
   machine precision). Never golden-test scores/loadings cross-language (sign+scale free) —
   compare **yhat / original-unit coef only**.
3. **Outlier gate primary path = DD-SIMCA χ²(Nh+Nq) data-driven** (PyChemAuth conventions):
   best behaviour at our n=18, gives the SD–OD acceptance plot directly, has extreme (α) +
   outlier (γ, Bonferroni) limits. Expose **classical F-form T²** and **Jackson–Mudholkar Q**
   as documented cross-check alternatives. **Leverage 2p/n (warn) / 3p/n (large-n)** as the
   unknown-sample extrapolation flag.
4. **EMSC = biospectools conventions** (OLS; reference column ordering tracked; value-based
   polynomial axis scaled to [−1,1]; divide by b). **Do NOT use the `chemometrics` lib as the
   reference** (asym-LS default, index-based axis, broken on sklearn ≥1.6, no `/b` by default).
5. **Nadeau–Bengio corrected paired t — implement ourselves.** `mlxtend.paired_ttest_resampled`
   is the *uncorrected* Dietterich formula; it is NOT a valid reference. Validate against the
   RNG-free golden vector (§7) / R `correctR`.
6. **Permutation test — between-LEVEL shuffle, hand-rolled.** `sklearn`'s `groups=` permutes
   *within* groups, which is wrong for our 6×3 replicate structure. p = `(#{perm≥obs}+1)/(n_perm+1)`.
7. **RPD conventions:** SD with **ddof=1**; report both RMSEP- and RMSECV-based RPD; RPIQ uses
   numpy percentile `method='linear'` (= R type 7); cite Chang 2001 thresholds.

## 5. Two code fixes this surfaced (fold into implementation)

- **SG derivative silently passes through** when `window ≤ polyorder` (triggers at n_features ≤ 4):
  add a `WARN` so the user knows no derivative was applied.
- **Document the ddof=0 choice** for SNV/autoscale in SKILL.md / `chemometrics.md`.

---

## 6. Per-method validation plan

### 6.1 Preprocessing — `test_preprocessing.py`
| Method | Target | Tier | Pin / note |
|---|---|---|---|
| SNV | `chemotools.scatter.StandardNormalVariate` | T1 | ddof=0; formula `(x−mean)/sd` per row (Barnes 1989, doi:10.1366/0003702894202201) |
| MSC | `chemotools…MultiplicativeScatterCorrection` (mean ref) | T1 | `(x−c)/m`, polyfit deg-1 (Geladi 1985, doi:10.1366/0003702854248656); training-mean ref frozen in-fold |
| SG d1/d2 | `scipy.signal.savgol_filter` (it IS the engine) | T0/T1 | odd window, `delta` (default 1.0 → per-POINT units; note the 1.864 cm⁻¹ step in the report), `mode='interp'`, poly<window; output ∝ 1/delta^deriv |
| centre/autoscale | `sklearn.StandardScaler` (ddof=0) | T1 | column-wise (van den Berg 2006, doi:10.1186/1471-2164-7-142) |

### 6.2 PLS — `test_pls.py`
- **Engine identity (T1):** our `pls_fit` coef/yhat == a direct `PLSRegression(n_components, scale=…)` fit → max diff 0.
- **Predict reconstruction (T0):** `yhat == (X − x_mean_) @ coef_.T + intercept_`, `intercept_ == y_mean`.
- **Gasoline regression (T2, ~1–2 %):** see §7(b). Requires `scale=False` + LOO; cross-language so tolerance only.
- Pin sklearn ≥1.3 (coef orientation flipped permanently in 1.3); `_coef_vector` already handles both.

### 6.3 PCR / PCA — `test_pcr_pca.py`
- PCA `explained_variance_ratio_` on centred X vs sklearn (T1); cap `n_components ≤ min(n,p)`.
- PCR back-map `components_.T @ reg.coef_` reproduces a manual PCA-then-OLS (T0/T1).
- mdatools `people` X/Y cum-variance as a secondary **T2** sanity (scale=TRUE → tolerance only).

### 6.4 Cross-validation — `test_cv.py`
- **RMSECV/R²cv/bias definitions (T0):** divisor N (ddof=0); `bias=mean(pred−y)`; identity
  `RMSEP² = SEP²·(n−1)/n + bias²` with `SEP` = ddof=1 residual SD.
- **Leakage property (T3):** in-fold-preprocessed RMSECV (MSC/autoscale) ≥ the leaky
  (pre-transform-then-CV) RMSECV on a constructed case; and == a hand-rolled leakage-safe loop (T1).
- **Scheme questions (T3):** LOO-over-replicated-y WARNs and points to `groups=`; leave-one-level-out
  uses level labels; k-fold clamps k≤n.

### 6.5 VIP — `test_vip.py`
- **Invariant (T3/T0):** `Σ_j VIP_j² = p` (→ mean(VIP²)=1, the >1 rule). Verified = 40.0 for p=40.
- Formula = Chong & Jun 2005 (doi:10.1016/j.chemolab.2004.12.011), line-identical to `mdatools::vipscores`.
- (No public per-variable VIP literals exist; optional R-generated vector can be added to golden/ later.)

### 6.6 Component selection — `test_*` (property, T3)
`choose_n_components`: picked RMSECV ≤ global_min·(1+tol); no within-tol candidate has fewer
components; ties → lower RMSECV.

### 6.7 Outlier / diagnostics (NEW) — `test_outliers.py`
Definitions/citations: T² (Mahalanobis in PC space, λ ddof=1); Q/SPE = ‖x−x̂‖²; Jackson–Mudholkar
1979 (doi:10.1080/00401706.1979.10489779); Box/Nomikos–MacGregor g·χ²_h (doi:10.1080/00401706.1995.10485888);
DD-SIMCA Pomerantsev & Rodionova 2014 (doi:10.1002/cem.2506).
- **Closed-form core (T0, no deps):** §7(a) — eigenvalues, per-sample T² & Q, all four T² limit
  forms, JM & Box Q limits, DD-SIMCA SD/h0. This is the always-runs anchor.
- **mdatools `people` Hotelling limits (T2 literal):** per-component 4.15961510, 6.85271430, …,
  37.07020869 (comp 11); reproduced exactly via `[A·31/(32−A)]·F_.95(A,32−A)`.
- **chemometrics peach NIR (T1, installed lib):** `x_residual_std_=0.005448557285711451`,
  `crit_dmodx(.95)=1.0697876760`, `crit_dhypx(.95)=11.1869800486`, `sum(leverage)=4.0`.
- **PyChemAuth iris-setosa DD-SIMCA (T2 literal):** `n_components=3, scale_x=False, robust='semi'`
  → Nh=3, Nq=1, h0=3.2926643176398445, q0=0.005740692404315845, c_crit=χ²_.95,4=9.487729036781154.

### 6.8 EMSC (NEW) — `test_emsc.py`
- **biospectools gold (T1, bit-exact 6e-16):** `datasets/emsc_testdata.xlsx`, `poly_order=4`,
  `reference=mean(raw)`, axis = the `wn` row (4000→499.547, 3631 pts); reproduce the "corrected
  quartic" + "residual" rows. (Afseth & Kohler 2012; Martens & Stark 1991.)
- **Synthetic coefficient recovery (T0):** ref=pure spectrum, poly_order=0 → recover input
  `b=[1.7,0.5,1.0]`; the readable sanity check. (Note: with ref=mean + baseline term, `b` is NOT
  the dilution factor — document; don't interpret `b` quantitatively unless ref=pure & poly=0.)
- Track coefficient column order (ref-first vs poly-first) and odd-poly sign convention.

### 6.9 Validation trio (NEW) — `test_trio.py`
- **Permutation (T1 vs sklearn mechanics, T3 for grouping):** p=`(#{perm≥obs}+1)/(n_perm+1)`;
  golden recipe score≈0.829, p≈0.00498; assert the between-level shuffle keeps replicate blocks intact.
- **Nadeau–Bengio (T0 closed-form):** corrected variance factor `(1/k + n_test/n_train)`;
  golden vector → t=1.364576, p=0.205528, df=9, correction=0.35. (Nadeau & Bengio 2003; Bouckaert & Frank 2004.)
- **RPD/RPIQ/RER (T0):** RPD=SD(ddof=1)/RMSEP; unit test SD=10, RMSEP=4 → 2.50; RPIQ≈1.349·RPD;
  identity `RMSEP²=SEP²(n−1)/n+bias²`. (Williams; Bellon-Maurel 2010; Chang 2001 thresholds.)

---

## 7. Concrete golden fixtures (self-contained where possible)

**(a) Dependency-free closed-form (the always-runs anchor).**
`X = [[2,0,1],[0,2,0],[3,1,2],[1,3,1],[4,0,3],[2,2,2]]`, n=6, m=3, **centre only**, A=2, α=0.05:
- eigenvalues (ddof=1): `[3.7812886, 0.7632012, 0.0221769]`
- T² per sample: `[2.4705981, 2.2419840, 0.3937636, 1.9042242, 2.1356383, 0.8537918]`
- T² limits: F (n²−1 form) `20.2541264` · F (textbook (n−1)/(n−k)) `17.3606798` · χ² `5.9914645` · Tracy–Young Beta `3.6011630`
- Q per sample: `[0.00065204,0.00572200,0.04670982,0.01321247,0.00808042,0.03650782]`
- Q limits: Jackson–Mudholkar (degenerate single residual eigenvalue, h0=1/3) `0.08309166` · Box/NM (g=0.009407,h=1.964643) `0.05571389`
- DD-SIMCA: mean(SD)=A/n=`0.3333` (confirms singular-value normalization); moment DoF Nh/Nq + Dcrit
  to be **recomputed and locked at implementation** (noisy at n=6 — record exact values then).

**(b) R `pls` gasoline (T2 literal; needs `scale=False`+LOO).** 60 samples, octane + NIR 401 λ
(900→1700 nm @2 nm; Kalivas 1997). Split train[1:50]/test[51:60]. LOO RMSEP (CV row):
`Int 1.545, 1c 1.357, 2c 0.2966, 3c 0.2524, 4c 0.2476, 5c 0.2398, 6c 0.2319, 7c 0.2386, 8c 0.2316, 9c 0.2449, 10c 0.2673`.
Test predict @2c: `87.94,87.25,88.16,84.97,85.15,84.51,87.56,86.85,89.19,87.09`. (Chosen ncomp=2.)
MSEP divisor n vs n−1 undocumented in R man pages → tolerance-only; do not derive bit-exact MSEP.

**(c) EMSC** `emsc_testdata.xlsx` (biospectools, validated vs Kohler MATLAB 06.02.2020).
**(d) mdatools `people`, (e) chemometrics peach, (f) PyChemAuth iris** — values in §6.7.

## 8. Edge-case matrix — `test_edge_cases.py`
Each row is a test; current module behaviour was verified live (✓ = already guarded).
1. SNV flat row (sd=0) → centred, no div0 ✓
2. autoscale flat column (sd=0) → /1 ✓ (matches StandardScaler)
3. SG window > n_features → clamp ✓; **window ≤ poly → silent passthrough (ADD WARN)**
4. MSC flat reference (b=0) → /1 guard ✓
5. PLS/PCR n_components ≥ rank / ≥ n−1 → capped at `min(nc,n−1,p)` ✓
6. PCA n_components > min(n,p) → capped ✓
7. LOO over replicated y → WARNs, directs to `groups=` ✓
8. k-fold k>n → clamped ✓
9. R²(CV) with constant y (SST≈0) → nan ✓
10. VIP with explained-Y≈0 → zeros ✓
11. **A ≥ rank** (no residual space) → JM Q_lim undefined (0^(1/0)); switch to χ²/dd, flag "no residual space"
12. **Perfect fit, Q≈0** → JM degenerate / Box 0/0 → Q_lim→0⁺, flag "Q test uninformative"
13. **Negative/zero eigenvalues** (numerical) → clamp <tol·λ_max, exclude from λ⁻¹ & θ sums
14. **n−p ≤ 0** → F/Beta undefined → fall back to χ² limit + warning
15. EMSC: poly 0/1/2/6; ref=mean vs pure; ±interferent term; constant/zero spectrum
16. Permutation: ties, small n_perm, grouped 6×3 (between-level), zero-variance y
17. Paired-t: zero variance of fold diffs, k=2, df
18. RPD: ddof on SD, RMSEP=0, tiny n
19. Real shape: 18 samples × 1798 wavenumbers (n≪p; F/χ² limits valid but wide — prefer robust DoF)

---

## 9. Build order (recommendation)
1. **Harness + existing-surface suite first** — DONE & green (33 checks). §5 fixes applied;
   caught + fixed the `pca_fit` n_components cap.
2. **Implement Tier-1 features + their tests together**, in priority order:
   - outlier gate (6.7) — **DONE & green** (`diagnostics`/`plot_influence`; bit-exact vs the
     closed-form fixture + scipy + published mdatools limits). `plot_influence` gained `axis_scale="sqrt"`.
   - EMSC (6.8) — **DONE & green** (`emsc`; bit-exact vs the biospectools/Kohler gold + closed-form
     b-recovery/interferent-removal).
   - trio (6.9) — **DONE & green** (`permutation_test`, `corrected_paired_t`, `rpd`/`rpiq`;
     Nadeau–Bengio golden reproduced through the function).
3. §5 code fixes — DONE (SG passthrough WARN; ddof=0 documented).

Suite total: **85 PASS / 0 FAIL / 1 SKIP** (gasoline dataset). Remaining Tier-1: the seeded
run-order schedule (DoE — property-validated, lives outside this chemometrics suite).

## 10. Citations
Barnes 1989 (10.1366/0003702894202201) · Geladi 1985 (10.1366/0003702854248656) · Savitzky–Golay
1964 (10.1021/ac60214a047) · van den Berg 2006 (10.1186/1471-2164-7-142) · Chong & Jun 2005
(10.1016/j.chemolab.2004.12.011) · de Jong SIMPLS 1993 (10.1016/0169-7439(93)85002-X) · Jackson &
Mudholkar 1979 (10.1080/00401706.1979.10489779) · Nomikos & MacGregor 1995
(10.1080/00401706.1995.10485888) · Tracy/Young/Mason 1992 · De Maesschalck 2000
(10.1016/S0169-7439(99)00047-7) · Pomerantsev & Rodionova 2014 (10.1002/cem.2506) · Kucheryavskiy
mdatools (10.1016/j.chemolab.2020.103937; distances 10.1016/j.chemolab.2021.104304) · Afseth &
Kohler 2012; Martens & Stark 1991 · Nadeau & Bengio 2003; Bouckaert & Frank 2004 · Bellon-Maurel
2010; Chang 2001 · Kalivas 1997 (gasoline). Reference code: svkucheryavski/mdatools,
maruedt/chemometrics, mahynski/pychemauth, yzontov/dd-simca, BioSpecNorway/biospectools.
