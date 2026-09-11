"""Gate for the crystal validation set. Needs gemmi + Dans_Diffraction (pymatgen optional);
exits 0 with a SKIP message when they are absent so a pytest bridge can skip cleanly.

    python skill_validation/crystal/run.py
"""
import os, sys, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
STEPS = ["test_engine.py", "crystal_validate.py", "test_adp.py", "test_realism_cif.py", "dans_fix.py"]

for dep in ("gemmi", "Dans_Diffraction"):
    try:
        __import__(dep)
    except ImportError:
        print(f"SKIP: {dep} not installed -- the crystal validation set needs gemmi + Dans_Diffraction")
        sys.exit(0)

failed = []
for step in STEPS:
    t0 = time.time()
    r = subprocess.run([sys.executable, os.path.join(HERE, step)], capture_output=True, text=True)
    ok = r.returncode == 0
    print(f"[{'PASS' if ok else 'FAIL'}] {step:24s} {time.time() - t0:5.1f} s")
    if not ok:
        failed.append(step)
        print(r.stdout[-2500:]); print(r.stderr[-1500:])
print("-" * 50)
print("ALL PASS" if not failed else f"{len(failed)} FAILED: {failed}")
sys.exit(1 if failed else 0)
