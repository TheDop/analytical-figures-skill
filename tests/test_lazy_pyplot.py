"""style imports pyplot on first use, so the numbers-only paths (module import, the CLI's
CSV calibrate) never pay the matplotlib import; and the figure path still works."""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _run(code):
    r = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def test_module_import_leaves_pyplot_out():
    out = _run("import sys\n"
               "from scripts import config, spectra, verify, calibration\n"
               "print('matplotlib.pyplot' in sys.modules, 'matplotlib' in sys.modules)")
    assert out.split() == ["False", "False"]


def test_cli_calibrate_leaves_pyplot_out(tmp_path):
    csv = tmp_path / "calib.csv"
    csv.write_text("conc,signal\n0,0.02\n10,0.55\n20,1.05\n30,1.52\n40,2.06\n50,2.49\n", encoding="utf-8")
    out = _run("import sys, cli\n"
               f"args = cli.build_parser().parse_args(['calibrate', {str(csv)!r}])\n"
               "assert args.func(args) == 0\n"
               "print('PYPLOT', 'matplotlib.pyplot' in sys.modules)")
    assert "PYPLOT False" in out


def test_figure_path_still_loads_pyplot():
    out = _run("import sys\n"
               "from scripts import config, style\n"
               "style.apply_style(config.Config())\n"
               "fig, _ = style.figure(config.Config())\n"
               "print('matplotlib.pyplot' in sys.modules, fig.get_size_inches()[0] > 0)")
    assert out.split() == ["True", "True"]
