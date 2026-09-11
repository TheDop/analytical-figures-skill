"""doe.py  -  seeded, balanced run-order schedules for a factor study.

The lab counterpart to the analysis modules: it emits the actual run sheet the team follows,
built so a confounder can't sneak in. The Group-C calibration confounded OPERATOR with
CONCENTRATION (each person did a whole level) -- this generator assigns operators by a rotating
round-robin so operator is ORTHOGONAL to concentration (each operator does each level ~equally),
then shuffles the RUN ORDER with a fixed seed so the exact schedule is reproducible and archivable
in the plan/appendix. numpy + csv only.

    run_order_schedule(levels, reps, operators, seed=, instruments=)  -> list of run dicts
    write_schedule_csv(schedule, path)                                -> path
    balance_report(schedule)                                          -> (operator x level) count table
"""
from __future__ import annotations
import numpy as np


def run_order_schedule(levels, reps, operators, seed=0, instruments=None):
    """Build a randomized, operator-balanced run-order schedule.

    levels: the concentration levels (any labels). reps: preps per level. operators: the team.
    instruments: optional list to balance across spectrometers too. seed: fixes the run order.
    Operator for prep r of level i = operators[(r + i) % n_operators] -- a rotation that makes
    operator orthogonal to concentration (when reps == n_operators, each operator does every level
    exactly once). The run ORDER is then a seeded permutation. Returns a list of run dicts."""
    rng = np.random.default_rng(seed)
    nop = len(operators)
    nin = len(instruments) if instruments else 0
    runs = []
    for i, lev in enumerate(levels):
        for r in range(reps):
            run = {"level": lev, "rep": r + 1, "operator": operators[(r + i) % nop]}
            if nin:
                run["instrument"] = instruments[(r + i) % nin]
            runs.append(run)
    order = rng.permutation(len(runs))
    schedule = []
    for run_no, idx in enumerate(order, start=1):
        row = {"run": run_no}
        row.update(runs[idx])
        schedule.append(row)
    return schedule


def balance_report(schedule):
    """(operator, level) -> count. A balanced design has all counts within 1 of each other."""
    table = {}
    for row in schedule:
        table[(row["operator"], row["level"])] = table.get((row["operator"], row["level"]), 0) + 1
    return table


def write_schedule_csv(schedule, path):
    """Write the schedule to CSV (the run sheet). Column order = the keys of the first row."""
    import csv
    import os
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    keys = list(schedule[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(schedule)
    return path


if __name__ == "__main__":
    # self-validation (property-based: determinism, operator<->level balance, coverage)
    ok = True
    levels = [15, 25, 35, 45]
    operators = ["A", "B", "C"]
    s1 = run_order_schedule(levels, reps=3, operators=operators, seed=7, instruments=["FTIR-1", "FTIR-2"])
    s2 = run_order_schedule(levels, reps=3, operators=operators, seed=7, instruments=["FTIR-1", "FTIR-2"])

    det = s1 == s2
    different = run_order_schedule(levels, 3, operators, seed=8) != [
        {k: v for k, v in r.items() if k != "instrument"} for r in s1]
    counts = balance_report(s1)
    balanced = (max(counts.values()) - min(counts.values())) <= 1
    full = all((op, lev) in counts for op in operators for lev in levels)
    n_ok = len(s1) == len(levels) * 3
    runs_seq = [r["run"] for r in s1] == list(range(1, len(s1) + 1))

    for name, cond in [("determinism (same seed -> identical)", det),
                       ("different seed -> different order", different),
                       ("operator x level balanced (max-min <= 1)", balanced),
                       ("every operator covers every level", full),
                       ("run count == levels*reps", n_ok),
                       ("run numbers are 1..N", runs_seq)]:
        print(("PASS" if cond else "FAIL"), name)
        ok = ok and cond
    print("counts (op,level):", dict(sorted(counts.items())))
    import sys
    sys.exit(0 if ok else 1)
