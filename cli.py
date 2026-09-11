#!/usr/bin/env python3
"""
cli.py  -  scripts-first QC shortcuts for the recurring flows.

The skill's primary deliverable is still a self-contained, editable analysis script (see
bundle.py). This is a thin CONVENIENCE for routine QC — "drop new .spc files in, get the band
%RSD" / "point at a conc,signal CSV, get the calibration figures of merit" — without editing a
script each time. It wraps the maintained modules; for a report figure, write the analysis.

    python cli.py spc-rsd  DIR_OR_GLOB --center 1748 [--metric area|height|deriv2] [--reference ref.spc]
    python cli.py calibrate calib.csv [--conf 0.95]      # CSV: two columns conc,signal (header ok)
"""
from __future__ import annotations
import os
import sys
import glob
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)        # so `from scripts import ...` (the modules use relative imports)


def _load_spc(path):
    import numpy as np
    try:
        import spc_spectra
    except ImportError as e:
        raise SystemExit("spc-rsd needs the `spc-spectra` package: pip install spc-spectra") from e
    s = spc_spectra.File(path)
    return np.array(s.x, float), np.array(s.sub[0].y, float)


def _files(dir_or_glob, ext="*.spc"):
    if os.path.isdir(dir_or_glob):
        return sorted(glob.glob(os.path.join(dir_or_glob, ext)))
    return sorted(glob.glob(dir_or_glob))


def cmd_spc_rsd(args):
    import numpy as np
    from scripts import config, spectra, verify
    files = _files(args.dir)
    if not files:
        raise SystemExit(f"no .spc files found at {args.dir!r}")
    cfg = config.Config(strict=False)
    reps = []
    for fp in files:
        x, y = _load_spc(fp)
        verify.check_trace(x, y, cfg, name=os.path.basename(fp))
        reps.append((x, y, fp))

    if args.metric == "deriv2":
        anchors = None
        vals = [spectra.band_metric(x, y, "deriv2", center=args.center, deriv_smooth=cfg.deriv_smooth)
                for x, y, _ in reps]
    else:
        rx, ry = (_load_spc(args.reference) if args.reference else (reps[0][0], reps[0][1]))
        anchors = spectra.find_anchors(rx, ry, args.center, gap=cfg.anchor_gap,
                                       maxhw=cfg.anchor_maxhw, smooth=cfg.anchor_smooth)
        vals = [spectra.band_metric(x, y, args.metric, anchors=anchors, center=args.center,
                                    deriv_smooth=cfg.deriv_smooth) for x, y, _ in reps]

    s = verify.summarize(vals, cfg)
    rsd = 100 * s["sd"] / s["mean"] if s["mean"] else float("nan")
    print(f"\nband {args.metric} @ {args.center} cm-1"
          + (f"  (anchors locked at {anchors[0]:.1f}/{anchors[1]:.1f} cm-1 on "
             f"{os.path.basename(args.reference) if args.reference else os.path.basename(reps[0][2])})"
             if anchors else "  (anchor-free)"))
    for (x, y, fp), v in zip(reps, vals):
        print(f"  {os.path.basename(fp):32s} {v: .6g}")
    print(f"\n  n={s['n']}  mean={s['mean']:.6g}  SD={s['sd']:.3g}  %RSD={rsd:.2f}  ({s['caption']})")
    return 0


def cmd_calibrate(args):
    import numpy as np
    import csv
    from scripts import config, calibration
    conc, sig = [], []
    with open(args.csv, newline="") as fh:
        for row in csv.reader(fh):
            if len(row) < 2:
                continue
            try:
                c, y = float(row[0]), float(row[1])
            except ValueError:
                continue            # header / non-numeric line
            conc.append(c); sig.append(y)
    if len(conc) < 3:
        raise SystemExit("calibrate needs >=3 numeric (conc,signal) rows")
    cfg = config.Config(strict=False, conf_level=args.conf)
    model = calibration.fit(np.array(conc, float), np.array(sig, float), cfg)
    ll = calibration.lod_loq(model, cfg)
    print(f"\ncalibration ({len(conc)} points, {len(set(conc))} levels)")
    print(f"  slope     = {model['slope']:.6g}")
    print(f"  intercept = {model['intercept']:.6g}")
    print(f"  R2        = {model['r2']:.5f}")
    print(f"  LOD       = {ll['lod']:.4g}")
    print(f"  LOQ       = {ll['loq']:.4g}")
    print("  (figures of merit only; for a report figure write the analysis + residual panel)")
    return 0


def build_parser():
    ap = argparse.ArgumentParser(prog="cli.py", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("spc-rsd", help="band area/height/deriv2 + %%RSD over a set of .spc replicates")
    a.add_argument("dir", help="directory of .spc files, or a glob")
    a.add_argument("--center", type=float, required=True, help="band centre in cm-1 (e.g. 1748)")
    a.add_argument("--metric", choices=("area", "height", "deriv2"), default="area")
    a.add_argument("--reference", default=None, help="optional .spc to lock the anchors on (default: first)")
    a.set_defaults(func=cmd_spc_rsd)

    b = sub.add_parser("calibrate", help="fit a calibration from a conc,signal CSV; print slope/R2/LOD/LOQ")
    b.add_argument("csv", help="CSV with two columns: concentration, signal (header allowed)")
    b.add_argument("--conf", type=float, default=0.95, help="confidence level (default 0.95)")
    b.set_defaults(func=cmd_calibrate)
    return ap


if __name__ == "__main__":
    args = build_parser().parse_args()
    raise SystemExit(args.func(args))
