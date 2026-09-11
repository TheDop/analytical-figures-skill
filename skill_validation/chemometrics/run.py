"""Zero-dep runner for the chemometrics validation suite.

    python skill_validation/chemometrics/run.py

Discovers test_*.py modules, calls each module's collect(), prints a tiered
PASS/FAIL/SKIP report with the matched value and provenance, and exits non-zero if
anything FAILed or ERRORed (SKIP is not a failure). See SPEC.md for the design.
"""
from __future__ import annotations
import sys
import importlib
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                         # headless: scripts.style imports pyplot

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                         # .../analytical-figures (holds the `scripts` package)
sys.path.insert(0, str(ROOT))                 # for `from scripts import ...`
sys.path.insert(0, str(HERE))                 # for `_harness`, `_data`, test modules

MODULES = ["test_preprocessing", "test_pls", "test_pcr_pca", "test_cv", "test_vip",
           "test_outliers", "test_emsc", "test_trio", "test_edge_cases", "test_classify"]


def main():
    results = []
    for name in MODULES:
        try:
            mod = importlib.import_module(name)
            results.extend(mod.collect())
        except Exception as e:           # an import/collect blow-up is itself a failure
            results.append(dict(area=name, name="<import/collect>", tier="-",
                                status="ERROR", note=repr(e), provenance="", got=None, want=None))
            traceback.print_exc()

    area = None
    for r in sorted(results, key=lambda r: (r["area"], r["name"])):
        if r["area"] != area:
            area = r["area"]
            print(f"\n== {area} ==")
        mark = {"PASS": "PASS", "FAIL": "FAIL", "SKIP": "skip", "ERROR": "ERR "}[r["status"]]
        extra = r.get("note") or ""
        prov = f"  [{r['provenance']}]" if r.get("provenance") else ""
        got = f"  {r['got']}" if r.get("got") not in (None, "") and r["status"] == "FAIL" else ""
        print(f"  {mark}  {r['tier']:2}  {r['name']}{got}  {extra}{prov}")

    counts = {k: sum(1 for r in results if r["status"] == k) for k in ("PASS", "FAIL", "SKIP", "ERROR")}
    print(f"\n{'-'*60}")
    print(f"PASS {counts['PASS']}   FAIL {counts['FAIL']}   SKIP {counts['SKIP']}   ERROR {counts['ERROR']}"
          f"   (total {len(results)})")
    sys.exit(1 if (counts["FAIL"] or counts["ERROR"]) else 0)


if __name__ == "__main__":
    main()
