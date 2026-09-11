# Pectin ATR-FTIR raw spectra (vendored, CC BY 4.0)

**Source:** 王 裕鑫 (Wang), *FT-IR raw data*, Mendeley Data, 2023, V1,
doi:[10.17632/gkwbp3wc49.1](https://doi.org/10.17632/gkwbp3wc49.1). Licence: Creative Commons
Attribution 4.0 International.

**What is here:** the 24 raw-spectrum CSVs from the record (two columns: wavenumber / cm⁻¹,
absorbance; 3736 points, 400–4000 cm⁻¹), numeric content unchanged, file names normalised
(line endings are normalised to LF by the repository's `.gitattributes`):

- `calibration/DM<x>.csv` — the six calibration standards of known degree of
  methyl-esterification (DM 3, 20, 37, 55, 62.8, 70.5 %).
- `samples/<group>-<rep>.csv` — the eighteen sample spectra: six pectin preparations
  (CA/E/W extraction × HWP/SWP) × three replicates.

**Not vendored:** the record's `.SPA` instrument files, baseline-correction `.TIF` images,
the authors' calibration workbook and the SPSS output (~80 MB). They add nothing the
validation or the gallery uses.

**Used by:** `../../validate_ftir_integration.py` (re-integration of the standards against
the authors' published band areas) and `docs/build_gallery.py` (the spectral tiles).
