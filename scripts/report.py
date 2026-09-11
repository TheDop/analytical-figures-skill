"""
report.py  -  turn the CONFIG + results into prose you can paste into a report.

    methods_report(cfg, results=None) -> (methods_text, si_table)

methods_text   a Methods paragraph describing the conventions the figure was made
               under (instrument axis convention, baseline, normalisation, the
               statistic behind any error bars, and the figure/software settings).
si_table       a plain-text table of whatever results dict you pass in.

This writes the FORM (how it was processed and rendered). It does NOT invent
findings or interpret the data - that judgement is yours.
"""
from __future__ import annotations


def _stat_sentence(cfg):
    return (f"Uncertainties are reported as {int(cfg.conf_level*100)}% confidence "
            f"intervals (Student's t, two-sided) unless stated otherwise; the "
            f"statistic and n are given with each value.")


def methods_report(cfg, results=None):
    parts = []

    if cfg.domain == "ftir":
        axis = ("Infrared spectra are presented with wavenumber on the abscissa "
                "(cm^-1), plotted high-to-low by convention")
        yax = ("transmittance (%)" if cfg.ftir_yaxis == "transmittance"
               else "absorbance")
        axis += f", and {yax} on the ordinate."
    elif cfg.domain == "pxrd":
        axis = ("Powder X-ray diffractograms are presented with diffraction angle "
                "2theta (degrees) on the abscissa, plotted low-to-high, and "
                "intensity on the ordinate.")
    else:
        axis = ""
    if axis:
        parts.append(axis)

    if cfg.baseline:
        bl = {"arpls": f"an asymmetrically reweighted penalised least-squares "
                       f"(arPLS, lambda={cfg.arpls_lam:g}) baseline",
              "rubberband": "a rubber-band (convex-hull) baseline"}.get(
                  cfg.baseline, f"a {cfg.baseline} baseline")
        parts.append(f"Each trace was baseline-corrected with {bl}, applied "
                     f"identically across every replicate in a set.")

    if cfg.integration_windows:
        wins = ", ".join(f"{n} ({min(a,b):g}-{max(a,b):g} cm^-1)"
                         for a, b, n in cfg.integration_windows)
        parts.append(f"Band areas were integrated above a local linear baseline "
                     f"between fixed window endpoints ({wins}), identical for all "
                     f"spectra so that replicate differences are not a windowing "
                     f"artefact.")

    if cfg.normalize and cfg.normalize != "none":
        how = {"max": "to unit maximum", "area": "to unit area"}.get(cfg.normalize, cfg.normalize)
        parts.append(f"Spectra were normalised {how} for display only; "
                     f"quantification used the un-normalised, baseline-corrected trace.")

    if cfg.smooth == "savgol":
        parts.append(f"A Savitzky-Golay filter (window {cfg.savgol_window}, order "
                     f"{cfg.savgol_poly}) was applied for display only; areas were "
                     f"measured on the unsmoothed trace.")

    parts.append(_stat_sentence(cfg))

    pal = {"okabe_ito": "the Okabe-Ito colour-blind-safe palette",
           "colorblind": "a colour-blind-safe palette"}.get(cfg.palette, cfg.palette)
    parts.append(f"Figures were prepared at final {cfg.journal} {cfg.column}-column "
                 f"dimensions using {pal} with redundant line-style encoding, and "
                 f"exported as vector graphics ({'/'.join(cfg.formats)}; raster "
                 f"proofs at {cfg.dpi_raster} dpi).")

    methods_text = " ".join(parts)

    # ---- SI table ----
    si_table = ""
    if results:
        rows = [(str(k), (f"{v:.4g}" if isinstance(v, (int, float)) else str(v)))
                for k, v in results.items()]
        w0 = max((len(k) for k, _ in rows), default=8)
        w1 = max((len(v) for _, v in rows), default=6)
        line = f"+{'-'*(w0+2)}+{'-'*(w1+2)}+"
        si_table = "\n".join(
            [line, f"| {'quantity'.ljust(w0)} | {'value'.ljust(w1)} |", line]
            + [f"| {k.ljust(w0)} | {v.ljust(w1)} |" for k, v in rows] + [line])

    return methods_text, si_table
