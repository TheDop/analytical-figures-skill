# calibration.md — honest linear calibration (ICH Q2 flavour)

`fit(x, y, cfg)` returns slope, intercept, R^2, residual SD, SE(slope), the
t-multiplier, the slope CI, and the calibrated x-range.

## The mandatory residual panel
`plot_calibration` draws the curve with a 95% mean-response CI band and a wider
prediction band, AND a residual panel beneath it. This is deliberate: a high R^2
can hide curvature. If the residuals show a smile/frown, the linear model is
wrong no matter what R^2 says. Don't disable `force_residual_panel` to make a
deadline.

## Extrapolation is refused
`predict(model, signal, cfg)` inverts the curve to a concentration but raises a
`GateError` if the result falls outside the calibrated range. Extend the
calibration rather than reading off the extrapolated line.

## LOD / LOQ
`lod_loq(model, cfg)` -> LOD = 3.3*sigma/m, LOQ = 10*sigma/m, with sigma = residual
SD by default (`lod_loq_method="residual_sd"`). The method string is returned
with the numbers so the figure/caption can state how they were derived. For the
blank-SD approach, pass blank replicates and compute sigma from those instead.

## Slope CI
The slope CI is `slope +/- t*SE(slope)` at `cfg.conf_level`. Useful for asking
whether a slope differs significantly from a theoretical value (e.g. a Nernstian
-59 mV/decade): if the theoretical value is outside the CI, the difference is
significant.

## Intercept significance (constant-bias check)
`intercept_test(model, cfg)` tests whether the intercept differs significantly from
zero — a two-sided t-test on `b/SE(b)` at (n-2) dof. `fit` returns `se_intercept`
and `intercept_ci`. Read it as: **CI includes 0 (p > 1-conf) -> no significant
constant bias** (the line is consistent with the origin); a significant intercept
flags a constant offset (baseline, matrix effect, or low-end lack-of-fit) to chase
down. Report it alongside R² + the residual panel — a great R² can still sit on a
biased line. (Worked example, an aspirin/lactose ATR feasibility set: ester-area intercept
p=0.064 — not significant at 95 % but borderline, mirroring the low-end curvature.)

## Robustness across the band metric (ICH Q2)
`verify.calibration_robustness(levels, center, cfg, check=..., reference=ref)` runs
the entire calibration under each band metric (area / height / 2nd-derivative) and
tabulates slope, R², LOD/LOQ and recovery side by side. Use it to pick the locked
primary metric from evidence and to *report* baseline choice as a robustness result
rather than a private decision. Always pass `reference` so anchor-based metrics use
one locked anchor pair (see spectra.md "Detect-once-and-LOCK") — re-detecting per
standard inflates the calibration scatter (on an aspirin/lactose ATR set, locking alone
took the mean %RSD from 22.6 % to 7.5 %).
