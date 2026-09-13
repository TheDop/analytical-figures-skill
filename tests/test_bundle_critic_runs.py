"""The bundler's critic pass must RUN, not be skipped: its try/except turns an import error
(the `scripts` package not on sys.path when bundle.py lives inside scripts/, or a missing
`import sys`) into an '[INFO] critic: skipped' line that nothing else notices. Pin the
'critic pass' line in the output, in both layouts."""
import os
import shutil
import subprocess
import sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BODY = "import numpy as np\nfrom scripts import config\ncfg = config.Config()\nprint('LOD', 1.0)\n"


def _run_bundle(bundle_path, tmp_path):
    body = tmp_path / "analysis.py"
    body.write_text(BODY, encoding="utf-8")
    r = subprocess.run([sys.executable, str(bundle_path), str(body), "-o", str(tmp_path / "sa.py")],
                       capture_output=True, text=True, encoding="utf-8", errors="replace", cwd=ROOT)
    return r.stdout + r.stderr


def test_critic_runs_with_bundle_at_skill_root(tmp_path):
    out = _run_bundle(os.path.join(ROOT, "bundle.py"), tmp_path)
    assert "critic pass" in out, out
    assert "critic: skipped" not in out, out


def test_critic_runs_with_bundle_inside_scripts(tmp_path):
    """The layout used by skill collections: bundle.py beside the modules, not above them."""
    skill = tmp_path / "skill"
    shutil.copytree(os.path.join(ROOT, "scripts"), skill / "scripts",
                    ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copy(os.path.join(ROOT, "bundle.py"), skill / "scripts" / "bundle.py")
    src = (skill / "scripts" / "bundle.py").read_text(encoding="utf-8")
    assert 'SCRIPTS = os.path.join(HERE, "scripts")' in src
    (skill / "scripts" / "bundle.py").write_text(
        src.replace('SCRIPTS = os.path.join(HERE, "scripts")', "SCRIPTS = HERE"), encoding="utf-8")
    out = _run_bundle(skill / "scripts" / "bundle.py", tmp_path)
    assert "critic pass" in out, out
    assert "critic: skipped" not in out, out
