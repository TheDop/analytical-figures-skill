"""
charts.py  -  generic line / bar for lab data, with the chart-choice guards that
matter. The value here is interception of bad chart choices, not the plotting API.

    line(series, cfg)        x/y line chart (e.g. timecourse, kinetics)
    bar(groups, cfg)         bar chart that REFUSES bare mean bars at low n -
                             it shows the individual points (and an error bar whose
                             statistic is named) instead. Mean bars at n<min hide
                             the distribution and overstate certainty.
"""
from __future__ import annotations
import numpy as np
from . import style, verify


def line(series, cfg, xlabel="x", ylabel="y", ax=None):
    """series: dict label -> (x, y)."""
    style.apply_style(cfg)
    if ax is None:
        fig, ax = style.figure(cfg)
    else:
        fig = ax.figure
    for label, (x, y) in series.items():
        verify.check_trace(np.asarray(x), np.asarray(y), cfg, name=label)
        ax.plot(x, y, marker="o", label=label)
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if len(series) > 1:
        ax.legend(loc="best")
    return fig, ax


def bar(groups, cfg, ylabel="value"):
    """groups: dict label -> list of replicate values.
    If every group has n >= cfg.min_n_for_mean_bar, draw mean bars with a NAMED
    error bar (mean ± CI) and overlay the points. If any group is below the
    threshold, REFUSE the bare bar and show a points-only plot (with the mean as
    a tick) - a mean bar at n<3 hides the distribution and overstates certainty."""
    style.apply_style(cfg)
    fig, ax = style.figure(cfg)
    labels = list(groups.keys())
    stats = {k: verify.summarize(v, cfg) for k, v in groups.items()}
    min_n = min(s["n"] for s in stats.values())
    pal = style._palette(cfg)
    xpos = np.arange(len(labels))

    low_n = min_n < cfg.min_n_for_mean_bar
    if low_n:
        print(f"  [WARN] charts: min n={min_n} < {cfg.min_n_for_mean_bar} - "
              f"refusing mean bars, showing individual points instead")
        for i, k in enumerate(labels):
            pts = np.asarray(groups[k], float)
            ax.scatter(np.full_like(pts, xpos[i]), pts,
                       color=pal[i % len(pal)], zorder=3)
            ax.plot([xpos[i] - 0.2, xpos[i] + 0.2], [stats[k]["mean"]] * 2,
                    color="k", lw=1.2)            # mean tick
        cap = f"points + mean tick (n={min_n}); CIs omitted at low n"
    else:
        means = [stats[k]["mean"] for k in labels]
        errs = [stats[k]["ci_halfwidth"] for k in labels]
        ax.bar(xpos, means, yerr=errs, capsize=2,
               color=[pal[i % len(pal)] for i in range(len(labels))],
               edgecolor="k", linewidth=0.5)
        for i, k in enumerate(labels):          # redundant point overlay
            pts = np.asarray(groups[k], float)
            ax.scatter(np.full_like(pts, xpos[i]), pts, color="k", s=6, zorder=3)
        cap = list(stats.values())[0]["caption"]

    ax.set_xticks(xpos); ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel)
    ax.annotate(cap, xy=(0.02, 0.98), xycoords="axes fraction",
                va="top", fontsize="small")
    return fig, ax
