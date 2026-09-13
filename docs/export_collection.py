"""
export_collection.py — regenerate this skill in the layout K-Dense-AI/scientific-agent-skills
requires (skill dir = SKILL.md + references/ scripts/ assets/ only; bundle.py and cli.py under
scripts/, the template under assets/templates/, the third-party notice under references/;
upstream-only suite mentions reworded). Run after any upstream change, then re-verify in the
clone: `skills-ref validate skills/analytical-figures` and
`python tests/run_all.py --isolated analytical-figures`.

    python docs/export_collection.py <path to a clone of scientific-agent-skills>
"""
import os, re, shutil, sys
SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # this repo
DST = os.path.join(sys.argv[1], "skills", "analytical-figures")       # <clone of scientific-agent-skills>
if os.path.exists(DST): shutil.rmtree(DST)
ign = shutil.ignore_patterns("__pycache__", "*.pyc")
os.makedirs(DST)
shutil.copy(f"{SRC}/SKILL.md", f"{DST}/SKILL.md")
shutil.copytree(f"{SRC}/references", f"{DST}/references", ignore=ign)
shutil.copytree(f"{SRC}/scripts", f"{DST}/scripts", ignore=ign)
shutil.copytree(f"{SRC}/assets", f"{DST}/assets", ignore=ign)
os.makedirs(f"{DST}/assets/templates")
shutil.copy(f"{SRC}/templates/analysis_template.py", f"{DST}/assets/templates/analysis_template.py")
shutil.copy(f"{SRC}/bundle.py", f"{DST}/scripts/bundle.py")
shutil.copy(f"{SRC}/cli.py", f"{DST}/scripts/cli.py")
shutil.copy(f"{SRC}/THIRD_PARTY_NOTICES.md", f"{DST}/references/third_party_notices.md")

def rw(path, pairs, must=True):
    s = open(path, encoding="utf-8").read()
    for a, b in pairs:
        if must: assert a in s, (path, a)
        s = s.replace(a, b)
    open(path, "w", encoding="utf-8", newline="\n").write(s)

# bundle.py: modules live beside it now
rw(f"{DST}/scripts/bundle.py", [('SCRIPTS = os.path.join(HERE, "scripts")', 'SCRIPTS = HERE                       # the modules live beside this file (scripts/)')])
# cli.py: skill root is the parent of scripts/
rw(f"{DST}/scripts/cli.py", [('sys.path.insert(0, HERE)        # so `from scripts import ...` (the modules use relative imports)',
                             'sys.path.insert(0, os.path.dirname(HERE))   # the skill root, so `from scripts import ...` resolves')])
# template: assets/templates/x.py -> skill root is three levels up
rw(f"{DST}/assets/templates/analysis_template.py",
   [('sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))',
     'sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))   # the skill root')])

# Upstream-only mentions (validation suites / pytest tests are not shipped: skill dirs hold only what an agent loads)
rw(f"{DST}/SKILL.md", [
    ("All validated by `skill_validation/chemometrics/` (100 checks).",
     "All validated upstream against published reference values and reference implementations (100 checks)."),
    ("(`tests/test_sources.py`); the network wrappers aren't.",
     "(upstream test suite); the network wrappers aren't."),
])
rw(f"{DST}/references/chemometrics.md", [
    ("Kohler gold (`skill_validation/chemometrics/test_emsc.py`).", "Kohler gold (upstream validation suite)."),
    ("Conventions (pinned, see `skill_validation/chemometrics/`): SD uses", "Conventions (pinned by the upstream validation suite): SD uses"),
    ("All validated in `skill_validation/chemometrics/test_trio.py`", "All validated in the upstream validation suite"),
])
rw(f"{DST}/references/crystal.md", [("in `skill_validation/crystal/test_realism_cif.py`.", "in the upstream validation suite.")])
rw(f"{DST}/scripts/chemometrics.py", [("(skill_validation/chemometrics/).", "(the upstream validation suite).")])

# SKILL.md path rewrites
md = f"{DST}/SKILL.md"
rw(md, [("`templates/analysis_template.py`", "`assets/templates/analysis_template.py`"),
        ("python bundle.py ", "python scripts/bundle.py "),
        ("`bundle.py`", "`scripts/bundle.py`"),
        ("`cli.py`", "`scripts/cli.py`")])
s = open(md, encoding="utf-8").read()
# Files section: drop the tests/ and skill_validation/ bullets (not shipped here), point upstream
lines = s.split("\n")
out, skip = [], False
for ln in lines:
    if re.match(r"^- `(tests/|skill_validation/)", ln):
        skip = True; continue
    if skip and ln.startswith("  "):
        continue
    skip = False
    out.append(ln)
s = "\n".join(out)
s = s.replace("- `scripts/bundle.py` — amalgamate", "- `scripts/bundle.py` — amalgamate", 1)
marker = "- `scripts/cli.py` — "
assert marker in s
# append the upstream pointer right after the cli bullet (find end of that bullet)
i = s.index(marker); j = s.index("\n- ", i + 1) if "\n- " in s[i+1:] else s.index("\n\n", i)
# find the end of the cli bullet: next line that starts with "- " or blank
k = i
for m in re.finditer(r"\n(?=- |\n)", s[i+1:]):
    k = i + 1 + m.start(); break
upstream = ("\n- Validation suites (`skill_validation/…`, ~250 numeric-parity checks against published reference\n"
            "  values), the pytest self-tests, the worked examples and the gallery live in the upstream\n"
            "  repository: https://github.com/TheDop/analytical-figures-skill")
s = s[:k] + upstream + s[k:]
s = s.replace("Validated by `skill_validation/chemometrics/` (100 checks).", "Validated upstream (chemometrics suite, 100 checks).")
s = s.replace("Validated by `skill_validation/pxrd/` (physics, dep-light) and", "Validated upstream (PXRD physics suite) and")
s = s.replace("Validated by `skill_validation/cocrystal/` (26 checks).", "Validated upstream (cocrystal suite, 26 checks).")
open(md, "w", encoding="utf-8", newline="\n").write(s)
print("lines:", s.count("\n") + 1)
print("leftover skill_validation/tests mentions in SKILL.md:", [l for l in s.split("\n") if "skill_validation" in l or "`tests/" in l])
