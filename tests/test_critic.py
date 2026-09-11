"""The Claude-Science-inspired adds: the number-provenance critic (Tier 3), the
robust MAD outlier flag, and the provenance stamps (bundle header + figure metadata)."""
import os
import numpy as np
import pytest
from scripts import verify

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _warned(findings):
    return any(s == "WARN" for s, _ in findings)


# ------------------------------------------------------------- Tier 3: the critic
def test_critic_flags_hardtyped_fom():
    findings = verify.check_number_provenance("ax.set_title('R2 = 0.98')")
    assert _warned(findings)


def test_critic_passes_interpolated_fom():
    # the number comes from an f-string {..} -> traces to the code -> clean
    findings = verify.check_number_provenance("ax.set_title(f'R2 = {fit[\"r2\"]:.2f}')")
    assert not _warned(findings)


def test_critic_ignores_fixed_band_identifier():
    # a band centre is a physical CONSTANT, not a computed result -> not flagged
    findings = verify.check_number_provenance("ax.set_xlabel('Aspirin ester C=O (1748 cm-1)')")
    assert not _warned(findings)


def test_critic_flags_caption_variable():
    findings = verify.check_number_provenance("caption = 'RMSEP 3.5 %w/w, recovery 103%'")
    assert _warned(findings)


def test_critic_flags_hardtyped_n_in_fstring_literal_segment():
    # rmsep is interpolated (fine) but n=3 is typed in the literal segment (flagged)
    findings = verify.check_number_provenance("ax.text(0, 0, f'RMSEP {rmsep:.1f} at n=3')")
    assert _warned(findings)


def test_critic_survives_syntax_error():
    findings = verify.check_number_provenance("def (:oops")
    assert _warned(findings)  # a parse failure is itself surfaced as a WARN


# ------------------------------------------------------- robust MAD outlier flag
def test_mad_flags_the_obvious_outlier():
    vals = [10.0, 10.2, 9.8, 10.1, 9.9, 20.0]
    r = verify.flag_outliers_mad(vals, n_mads=3.5)
    assert r["flagged_index"] == [5]
    assert r["n_flagged"] == 1
    assert r["lower"] < r["median"] < r["upper"]


def test_mad_clean_set_flags_nothing():
    r = verify.flag_outliers_mad([5.0, 5.1, 4.9, 5.05, 4.95], n_mads=3.5)
    assert r["n_flagged"] == 0


def test_mad_zero_mad_fallback_no_crash():
    # >half identical -> MAD collapses to 0; the mean-abs-dev fallback still catches 7
    r = verify.flag_outliers_mad([5.0, 5.0, 5.0, 5.0, 7.0], n_mads=3.5)
    assert r["scaled_mad"] > 0
    assert r["mask"][4]


def test_mad_ignores_nan_and_keeps_length():
    r = verify.flag_outliers_mad([1.0, 2.0, 3.0, np.nan, 100.0], n_mads=3.5)
    assert r["mask"].shape == (5,)
    assert not r["mask"][3]        # the NaN slot is never flagged
    assert r["mask"][4]            # the real outlier is


def test_mad_uses_cfg_default():
    from scripts.config import Config
    cfg = Config(outlier_mad_n=100.0)     # absurdly lax -> nothing flagged
    r = verify.flag_outliers_mad([10.0, 10.2, 9.8, 20.0], cfg=cfg)
    assert r["n_mads"] == 100.0 and r["n_flagged"] == 0


# ----------------------------------------------------------- provenance stamps
def test_bundle_stamps_provenance():
    import bundle
    from scripts.config import SKILL_VERSION
    body = os.path.join(ROOT, "templates", "analysis_template.py")
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "standalone.py")
        bundle.bundle(body, out, description="unit-test figure", critic=False)
        src = open(out, encoding="utf-8").read()
    assert "provenance" in src
    assert f"analytical-figures v{SKILL_VERSION}" in src
    assert "shows:  unit-test figure" in src
    assert "python " in src        # env line present
    compile(src, out, "exec")      # still valid Python after the header


def test_bundle_extracts_body_description(tmp_path):
    import bundle
    body = tmp_path / "an.py"
    body.write_text('DESCRIPTION = "from the body"\nx = 1\n', encoding="utf-8")
    out = tmp_path / "out.py"
    bundle.bundle(str(body), str(out), critic=False)
    assert "shows:  from the body" in out.read_text(encoding="utf-8")


def test_savefig_embeds_description(tmp_path):
    Image = pytest.importorskip("PIL.Image")
    import matplotlib.pyplot as plt
    from scripts import style
    from scripts.config import Config
    cfg = Config(formats=("png",), description="Aspirin ester C=O calibration, n=3")
    style.apply_style(cfg)
    fig, ax = style.figure(cfg)
    ax.plot([0, 1, 2], [0, 1, 2])
    ax.set_xlabel("x (a.u.)"); ax.set_ylabel("y (a.u.)")
    base = str(tmp_path / "prov")
    written = style.save_fig(fig, base, cfg)
    plt.close(fig)
    assert base + ".png" in written
    img = Image.open(base + ".png")
    blob = " ".join(list(getattr(img, "text", {}).values())
                    + [str(v) for v in img.info.values()])
    assert "Aspirin ester" in blob
