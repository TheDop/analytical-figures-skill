# charts.md — generic line/bar, with the guards that matter

The value is interception of bad chart choices, not the plotting API.

## `bar(groups, cfg)` and the low-n guard
`groups` is `label -> list of replicate values`. If every group has
n >= `cfg.min_n_for_mean_bar` (default 3), it draws mean bars with a NAMED error
bar (mean +/- CI) and overlays the individual points. If ANY group is below the
threshold it REFUSES the bare bar and shows points + a mean tick instead — a
mean bar at n<3 hides the distribution and overstates certainty. The caption
always names the statistic and n.

## `line(series, cfg)`
`series` is `label -> (x, y)`. Each trace is gated with `check_trace` first.

## Anti-patterns to intercept (advise the user, don't just plot)
- Dual-Y axes that imply a correlation that isn't there.
- Pie charts for anything you'd compare quantitatively.
- Truncated y-axes that exaggerate small differences.
- Rainbow/jet colormaps (not perceptually uniform, not colour-blind safe).
- Bare mean bars at small n (the guard above handles this one programmatically).
