"""House-style packaging: .mplstyle export round-trips, shipped assets load, check_figure_width."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts import config, style

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_export_mplstyle_roundtrip(tmp_path):
    cfg = config.Config(journal="acs")
    p = style.export_mplstyle(cfg, str(tmp_path / "s.mplstyle"))
    plt.style.use(str(p))
    assert plt.rcParams["font.size"] == cfg.base_fontsize
    assert plt.rcParams["pdf.fonttype"] == 42
    assert plt.rcParams["svg.fonttype"] == "none"
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    assert colors[0].lstrip("#").lower() == "000000"     # first Okabe-Ito entry survived


def test_shipped_assets_load():
    for j in ("nature", "acs", "ieee", "general"):
        plt.style.use(os.path.join(ROOT, "assets", "styles", f"analytical-{j}.mplstyle"))


def test_check_figure_width():
    cfg = config.Config(journal="acs", column="single")
    style.apply_style(cfg)
    fig, _ = style.figure(cfg)
    ok, msg = style.check_figure_width(fig, cfg)
    assert ok, msg
    fig2, _ = style.figure(config.Config(figsize=(5.0, 3.0)))
    ok2, _ = style.check_figure_width(fig2, cfg)
    assert not ok2
