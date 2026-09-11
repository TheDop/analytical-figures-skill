#!/usr/bin/env python3
"""
bundle.py  -  amalgamate the skill's modules + an analysis body into ONE
self-contained Python script.

The skill's deliverable is not a PNG and not a script that imports this skill -
it is a single runnable file that depends only on the scientific stack
(numpy, matplotlib, optionally scipy/scienceplots). The maintained source of
truth is scripts/*.py; you write your analysis against those modules, then:

    python bundle.py my_analysis.py -o standalone.py     # hand over standalone.py

How it works: each module's external imports are hoisted and de-duplicated at the
top, intra-skill imports (from . import ..., from scripts ... ) are dropped, and
the `module.` qualifier is stripped from references so everything lives in one
flat namespace. The analysis body is appended last (with any sys.path shim
removed).
"""
import os
import re
import argparse

SKILL = ["config", "style", "verify", "spectra", "calibration", "chemometrics",
         "charts", "report", "crystal_engine", "crystal_pxrd", "crystal_view",
         "pxrd_realism", "cocrystal"]
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "scripts")

# matches a top-level import that pulls in one of THIS skill's own modules
_SKILL_IMPORT = re.compile(
    r'^\s*(?:'
    r'from\s+\.\S*\s+import\b'                                  # from . / from .mod import
    r'|from\s+scripts(?:\.\S+)?\s+import\b'                     # from scripts[.mod] import
    r'|from\s+(?:' + "|".join(SKILL) + r')\s+import\b'          # from <mod> import (flat layout)
    r'|import\s+(?:' + "|".join(SKILL) + r')\b'                 # import <mod>
    r')'
)


def _is_toplevel_import(line):
    if line[:1] in (" ", "\t"):          # indented -> lazy import, keep in body
        return False
    s = line.lstrip()
    return s.startswith("import ") or s.startswith("from ")


def _real_import_spans(src):
    """Line spans (1-based, inclusive) of GENUINE module-level import statements, via the AST.

    A line-based scan cannot do this job. Any unindented line inside a docstring that happens to
    begin 'from ' or 'import ' is indistinguishable from a real import, and hoisting it tears the
    docstring in half and produces a standalone that will not parse. Hit for real by a docstring
    whose second line began "from the co-former, ...". The AST also gets multi-line parenthesised
    imports right, which the line scan silently split into invalid fragments.

    Returns None if the source does not parse, so the caller can fall back to the line scan.
    """
    import ast
    try:
        tree = ast.parse(src)
    except SyntaxError:
        return None
    spans = []
    for node in tree.body:                       # module level only; lazy imports stay in place
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            spans.append((node.lineno, node.end_lineno or node.lineno))
    return spans


def _split(src):
    """Return (hoisted_external_imports, body_without_those_imports)."""
    lines = src.splitlines()
    spans = _real_import_spans(src)
    if spans is None:                            # unparseable: old line-based behaviour
        hoist, body = [], []
        for line in lines:
            if _SKILL_IMPORT.match(line):
                continue
            (hoist if _is_toplevel_import(line) else body).append(
                line.strip() if _is_toplevel_import(line) else line)
        return hoist, "\n".join(body)

    # The AST decides only what gets HOISTED. Dropping the skill's own imports still has to run
    # per line, because those also appear LAZILY inside functions (`from . import config`), which
    # never show up in tree.body and would otherwise survive into the standalone.
    consumed = {n for a, b in spans for n in range(a, b + 1)}
    hoist, body = [], []
    for a, b in spans:
        stmt = "\n".join(lines[a - 1:b])
        if _SKILL_IMPORT.match(stmt.splitlines()[0]):
            continue                             # drop the skill's own imports entirely
        hoist.append(stmt.strip() if a == b else stmt)
    for i, line in enumerate(lines, 1):
        if i in consumed:
            continue
        if _SKILL_IMPORT.match(line):
            continue                             # lazy/indented skill import
        body.append(line)
    return hoist, "\n".join(body)


def _debind(text, extra=()):
    """Strip the `module.` qualifier (style./verify./...) so calls resolve in the
    single flat namespace. `extra` adds import ALIASES that point at a skill module
    (e.g. `chemometrics as cm` -> also strip `cm.`). The negative lookbehind leaves
    things like `plt.style.use` untouched (the dot before 'style' blocks the match)."""
    names = list(SKILL) + [a for a in extra if a]
    return re.sub(r'(?<![\w.])(' + "|".join(names) + r')\.', '', text)


def _collect_aliases(src):
    """Map import aliases that point at one of THIS skill's modules, e.g.
    `from scripts import chemometrics as cm` -> {'cm': 'chemometrics'}, so `cm.` gets
    de-bound too. Scans the same skill-import lines that _split() drops."""
    out = {}
    for raw in src.splitlines():
        s = raw.strip()
        m = re.match(r'from\s+(?:\.\w*|scripts(?:\.\w+)?)\s+import\s+(.+)', s)
        if m:
            for part in m.group(1).split(','):
                am = re.match(r'(\w+)\s+as\s+(\w+)$', part.strip())
                if am and am.group(1) in SKILL:
                    out[am.group(2)] = am.group(1)
        m2 = re.match(r'import\s+scripts\.(\w+)\s+as\s+(\w+)$', s)
        if m2 and m2.group(1) in SKILL:
            out[m2.group(2)] = m2.group(1)
    return out


def _skill_version():
    try:
        from scripts.config import SKILL_VERSION
        return SKILL_VERSION
    except Exception:
        return "?"


def _env_versions(names=("numpy", "scipy", "matplotlib", "pandas", "scikit-learn")):
    """The installed versions of the scientific stack THIS bundle was produced against
    - the 'environment that produced it'. Only packages actually present are listed."""
    out = []
    try:
        from importlib.metadata import version, PackageNotFoundError  # noqa: F401
    except Exception:
        return out
    for n in names:
        try:
            out.append(f"{n} {version(n)}")
        except Exception:
            pass
    return out


def _extract_description(body_path):
    """A module-level DESCRIPTION / FIGURE_DESCRIPTION string constant in the analysis
    body, so the plain-language description lives WITH the analysis (CONFIG-block
    spirit) and is stamped without a separate flag."""
    import ast
    try:
        tree = ast.parse(open(body_path, encoding="utf-8").read())
    except Exception:
        return None
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant) \
                and isinstance(node.value.value, str):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id in ("DESCRIPTION", "FIGURE_DESCRIPTION"):
                    return node.value.value
    return None


def _provenance_header(body_path, description):
    """The reproducibility stamp for the deliverable (Claude Science: 'the exact code
    and environment that produced it, plus a plain-language description'). Version pins
    make the env that produced the figure recoverable months later; the description
    says what the figure shows. Comments only, so it sits fine above `from __future__`."""
    env = "  ".join(_env_versions())
    desc = description or _extract_description(body_path) \
        or '(add DESCRIPTION="..." to the analysis body, or pass --describe)'
    return "\n".join([
        "#",
        "# --- provenance (auto-stamped by analytical-figures) ---",
        f"# skill:  analytical-figures v{_skill_version()}",
        f"# env:    python {_pyver()}" + (f"  |  {env}" if env else ""),
        f"# source: {os.path.basename(body_path)}",
        f"# shows:  {desc}",
        "# reproduce: run this file as-is; edit the CONFIG block to retune.",
        "# -------------------------------------------------------",
    ])


def _pyver():
    import platform
    return platform.python_version()


def _run_critic(body_path):
    """Actor-critic at handoff: a static reviewer pass over the analysis body flagging
    figures-of-merit hard-typed onto the figure/caption instead of interpolated from
    the code (verify.check_number_provenance). Best-effort - never blocks the bundle."""
    try:
        from scripts import verify
        print("critic pass (number provenance):")
        verify.check_number_provenance(body_path)
    except Exception as e:
        print(f"  [INFO] critic: skipped ({e})")


def _shadow_check(body_src, skill_srcs):
    """Names the analysis body defines that will SHADOW a skill function after de-binding.

    _debind strips the `module.` qualifier, so `style.figure(cfg, 2, 3)` becomes `figure(cfg, 2, 3)`
    in the bundle. If the body also defines its own `def figure(panels)`, the two collide and the
    standalone dies with a TypeError -- while the ORIGINAL script keeps working, because there the
    call is still qualified. That makes it invisible until someone runs the deliverable, which is
    the one artefact that goes in a report appendix. Hit for real by a POM figure whose body
    defined figure() and called style.figure().
    """
    import ast
    try:
        body_tree = ast.parse(body_src)
    except SyntaxError:
        return []
    defined = {n.name for n in body_tree.body
               if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}
    # ...and any name the body ASSIGNS anywhere: `fit = calibration.fit(x, y, cfg)` de-binds to
    # `fit = fit(x, y, cfg)`, which is an UnboundLocalError inside a function and a silent
    # overwrite at module level. Hit for real by a worked example in docs/examples/.
    for node in ast.walk(body_tree):
        targets = []
        if isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign, ast.For, ast.comprehension)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        for t in targets:
            for leaf in ast.walk(t):
                if isinstance(leaf, ast.Name):
                    defined.add(leaf.id)
    if not defined:
        return []
    # what the body calls as `<skillmodule>.<name>` -- only those get de-bound into bare names
    used = set()
    for node in ast.walk(body_tree):
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) \
                and node.value.id in SKILL:
            used.add(node.attr)
    return sorted(defined & used)


def bundle(body_path, out_path, description=None, critic=True):
    hoist, chunks = [], []

    body = open(body_path).read()
    for name in _shadow_check(body, None):
        print(f"  [FAIL] '{name}' is defined or assigned in the analysis body AND used as a skill call "
              f"(module.{name}). After de-binding both become '{name}' and the STANDALONE will "
              f"break, though this script still runs. Rename the one in the body.")
    body = "\n".join(l for l in body.splitlines() if "sys.path" not in l)

    # collect skill-module aliases across ALL sources first, so _debind strips them too
    aliases = {}
    for mod in SKILL:
        p = os.path.join(SCRIPTS, mod + ".py")
        if os.path.exists(p):
            aliases.update(_collect_aliases(open(p).read()))
    aliases.update(_collect_aliases(body))
    extra = tuple(aliases)

    def add(src, title):
        h, b = _split(src)
        hoist.extend(h)
        chunks.append(f"# {'='*70}\n# {title}\n# {'='*70}\n" + _debind(b, extra).strip() + "\n")

    for mod in SKILL:                                  # dependency order
        p = os.path.join(SCRIPTS, mod + ".py")
        if os.path.exists(p):
            add(open(p).read(), f"from {mod}.py")
    add(body, f"analysis body: {os.path.basename(body_path)}")

    # de-dupe imports, keep __future__ first (it must precede other statements)
    seen, future, hdr = set(), [], []
    for i in hoist:
        if i in seen:
            continue
        seen.add(i)
        (future if i.startswith("from __future__") else hdr).append(i)

    out = ("#!/usr/bin/env python3\n"
           "# Self-contained, auto-bundled by analytical-figures/bundle.py.\n"
           "# Depends only on the scientific stack; edit the CONFIG block to retune.\n"
           + _provenance_header(body_path, description) + "\n\n"
           + "\n".join(future + hdr) + "\n\n\n" + "\n\n".join(chunks))
    with open(out_path, "w") as f:
        f.write(out)
    if critic:
        _run_critic(body_path)
    return out_path


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("body")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--describe", default=None,
                    help="one-line plain-language description stamped into the provenance "
                         "header (else a DESCRIPTION=... in the body, else a placeholder)")
    ap.add_argument("--no-critic", action="store_true",
                    help="skip the static number-provenance critic pass over the body")
    a = ap.parse_args()
    print("wrote", bundle(a.body, a.out, description=a.describe, critic=not a.no_critic))
