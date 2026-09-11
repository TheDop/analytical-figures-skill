# Third-party notices

## scipilot-figure-skill (MIT)

This skill descends from **scipilot-figure-skill** by Haojae
(<https://github.com/Haojae/scipilot-figure-skill>). The adapted portions are:

- the structure of `references/figure_selection.md` (after its `chart_selection.md`:
  *decide what to plot before deciding how*), re-grounded for analytical chemistry;
- in `scripts/style.py`, the panel-label alignment trick (anchor every letter at its
  axes' top-left corner, push by one shared points offset) and the constrained→tight
  layout fallback.

Its licence, reproduced as required:

```
MIT License

Copyright (c) 2026 Haojae

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Vendored data

- **Pectin ATR-FTIR raw spectra** — 王 裕鑫 (Wang), *FT-IR raw data*, Mendeley Data, 2023, V1,
  doi:10.17632/gkwbp3wc49.1, **CC BY 4.0**. The 24 raw-spectrum CSVs (6 calibration standards,
  18 samples) are vendored in `skill_validation/ftir_integration/data/pectin/`, content unchanged,
  file names normalised; see the README there. The record's other files are not included.

## Validation data fetched on demand (not redistributed)

- **biospectools EMSC test data** — `tests/data/emsc_testdata.xlsx` from
  BioSpecNorway/biospectools (<https://github.com/BioSpecNorway/biospectools>), the numeric
  gold for `skill_validation/chemometrics/test_emsc.py`. The repository declares no licence,
  so the file is **not** vendored here; `skill_validation/chemometrics/datasets/fetch.py emsc`
  downloads it, and the check skips cleanly when it is absent.

## Cited methods

Published methods the code implements (March–Dollase, Caglioti, Nadeau–Bengio corrected
t-test, DD-SIMCA, EMSC, Cruz-Cabeza ΔpKa rule, …) are cited in the docstring of the function
that implements each one and in `references/`.
