"""The QC CLI builds its parser and the calibrate command runs end-to-end."""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_parser_builds():
    import cli
    p = cli.build_parser()
    assert p is not None


def test_calibrate_command(tmp_path):
    import cli
    csv = tmp_path / "calib.csv"
    csv.write_text("conc,signal\n0,0.02\n10,0.55\n20,1.05\n30,1.52\n40,2.06\n50,2.49\n", encoding="utf-8")
    args = cli.build_parser().parse_args(["calibrate", str(csv)])
    assert args.func(args) == 0
