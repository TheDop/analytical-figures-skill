"""pytest config for the skill self-tests: headless matplotlib + skill root on sys.path, and
`run_py` for the tests that must run in a fresh interpreter (import-time properties, the CLI)."""
import os
import subprocess
import sys
import matplotlib

matplotlib.use("Agg")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the skill root
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def run_py(*args, code=None, cwd=ROOT, check=True):
    """Run this interpreter from the skill root: `run_py(code=...)` for `-c`, `run_py("cli.py",
    "--help")` for a script. Returns the CompletedProcess (text, captured); asserts exit 0 unless
    check=False."""
    argv = [sys.executable] + (["-c", code] if code is not None else []) + list(args)
    r = subprocess.run(argv, cwd=cwd, capture_output=True, text=True)
    if check:
        assert r.returncode == 0, r.stderr
    return r
