"""bundle.py must not mistake prose for an import.

A line-based scan cannot tell a real `from x import y` from a docstring line that merely begins
"from ...". Hoisting the latter tears the docstring in half and the standalone will not parse.
Hit for real by a docstring whose second line began "from the co-former, two things would
follow...". These tests pin the AST-based behaviour and the two things it must not break:
lazy skill imports inside functions must still be dropped, and multi-line parenthesised imports
must survive intact.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def _roundtrip(tmp_path, src, name="body.py"):
    import bundle
    p = tmp_path / name
    p.write_text(src, encoding="utf-8")
    out = tmp_path / "standalone.py"
    bundle.bundle(str(p), str(out), critic=False)
    return out.read_text(encoding="utf-8")


def test_docstring_line_starting_with_from_survives(tmp_path):
    src = ('"""Title line.\n\n'
           'A cocrystal and a salt differ by ONE PROTON. If the acceptor took the proton\n'
           'from the co-former, two things would follow and both are testable.\n'
           'import this line is prose too.\n"""\n'
           'import numpy as np\n\n'
           'DESCRIPTION = "x"\n\n'
           'def main():\n    return np.array([1.0])\n')
    out = _roundtrip(tmp_path, src)
    compile(out, "standalone.py", "exec")            # the actual regression: it must PARSE
    assert "from the co-former, two things" in out
    assert "import this line is prose too." in out


def test_lazy_skill_imports_are_still_dropped(tmp_path):
    src = ('"""Doc."""\n'
           'import numpy as np\n\n'
           'def f():\n'
           '    from . import config\n'
           '    from scripts import style\n'
           '    return np.array([1.0])\n')
    out = _roundtrip(tmp_path, src)
    compile(out, "standalone.py", "exec")
    assert "from . import config" not in out
    assert "from scripts import style" not in out


def test_multiline_import_is_not_split(tmp_path):
    src = ('"""Doc."""\n'
           'from math import (\n    sin,\n    cos,\n)\n\n'
           'def f():\n    return sin(0.0) + cos(0.0)\n')
    out = _roundtrip(tmp_path, src)
    compile(out, "standalone.py", "exec")            # a split parenthesised import would not
    assert "sin" in out and "cos" in out
