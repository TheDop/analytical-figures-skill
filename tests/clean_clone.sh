#!/bin/bash
# Clean-clone release check: a fresh venv with numpy/matplotlib/scipy/pytest only, a clone from
# GitHub, pytest + cli + the three worked examples (stage A: the crystal example must fail with a
# clear "needs gemmi" message), then `pip install -r requirements.txt` + the crystal deps and the
# same again (stage B). Usage: bash tests/clean_clone.sh <scratch dir>   (needs `gh` logged in).
# Run it before every release; read the printed pass/skip counts against the last release.
set -u
R="$1"; rm -rf "$R"; mkdir -p "$R"; cd "$R"
python -m venv venv || exit 1
PY="$R/venv/Scripts/python"
"$PY" -m pip install -q --disable-pip-version-check numpy matplotlib scipy pytest 2>&1 | tail -2
gh repo clone TheDop/analytical-figures-skill skill -- -q || exit 1
cd skill
echo "== commit: $(git log --oneline | head -1)"
echo "== stage A: pytest with numpy/matplotlib/scipy/pytest only =="
"$PY" -m pytest tests/ -q -p no:cacheprovider 2>&1 | tail -6
echo "== stage A: cli --help =="
"$PY" cli.py --help 2>&1 | head -2
echo "== stage A: example standalones =="
for d in bmc_lod pectin_calibration cif_report; do
  "$PY" "docs/examples/$d/${d}_standalone.py" > "$R/${d}_A.log" 2>&1; echo "  $d exit $?  ($(grep -c 'FAIL\]' "$R/${d}_A.log") FAIL lines; tail: $(tail -1 "$R/${d}_A.log" | cut -c1-120))"
done
echo "== stage B: pip install -r requirements.txt =="
"$PY" -m pip install -q --disable-pip-version-check -r requirements.txt 2>&1 | tail -3
echo "== stage B: crystal deps =="
"$PY" -m pip install -q --disable-pip-version-check gemmi Dans-Diffraction 2>&1 | tail -2
"$PY" -m pip list 2>/dev/null | grep -i -E "^(numpy|matplotlib|scipy|scikit-learn|rdkit|gemmi|Dans|pybaselines|scienceplots|pymatgen|pyvista) " 
echo "== stage B: pytest with full deps =="
"$PY" -m pytest tests/ -q -p no:cacheprovider 2>&1 | tail -6
echo "== stage B: cif example =="
"$PY" docs/examples/cif_report/cif_report_standalone.py > "$R/cif_report_B.log" 2>&1; echo "  cif_report exit $?  tail: $(tail -1 "$R/cif_report_B.log" | cut -c1-120)"
echo "== DONE =="
