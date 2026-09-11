"""SKILL.md keeps its required sections and every referenced doc/script exists."""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

REQUIRED_SECTIONS = ["## Output contract", "## Workflow", "## Hard principles",
                     "## Chemometrics", "## Files"]
REFERENCES = ["figure_selection", "cocrystal_id", "databases", "verification", "spectra",
              "calibration", "chemometrics", "charts", "crystal", "adding_a_family"]
SCRIPTS = ["config", "style", "verify", "spectra", "calibration", "chemometrics",
           "charts", "report", "doe", "crystal_engine", "crystal_pxrd", "crystal_view",
           "cocrystal", "pxrd_realism", "sources"]


def test_skill_md_sections():
    md = open(os.path.join(ROOT, "SKILL.md"), encoding="utf-8").read()
    missing = [h for h in REQUIRED_SECTIONS if h not in md]
    assert not missing, f"SKILL.md missing sections: {missing}"


def test_reference_docs_exist():
    missing = [f for f in REFERENCES if not os.path.exists(os.path.join(ROOT, "references", f + ".md"))]
    assert not missing, f"missing reference docs: {missing}"


def test_scripts_exist():
    missing = [m for m in SCRIPTS if not os.path.exists(os.path.join(ROOT, "scripts", m + ".py"))]
    assert not missing, f"missing scripts: {missing}"


def test_mplstyle_assets_exist():
    for j in ("nature", "acs", "ieee", "general"):
        assert os.path.exists(os.path.join(ROOT, "assets", "styles", f"analytical-{j}.mplstyle"))
