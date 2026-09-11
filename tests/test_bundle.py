"""bundle.py amalgamates the template into a standalone that (a) has no skill imports left,
(b) inlines the modules, and (c) is syntactically valid. We compile (not exec) — executing a
bundled analysis would try to load data / render figures."""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_bundle_roundtrip(tmp_path):
    import bundle
    body = os.path.join(ROOT, "templates", "analysis_template.py")
    out = tmp_path / "standalone.py"
    bundle.bundle(body, str(out))
    src = out.read_text(encoding="utf-8")

    # no skill imports survive the flattening
    assert "from scripts" not in src
    assert "from . import" not in src
    # modules were actually inlined (a couple of signature anchors)
    assert "def apply_style" in src
    assert "def save_fig" in src
    # and it parses
    compile(src, str(out), "exec")
