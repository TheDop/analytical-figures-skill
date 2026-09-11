# databases.md — where to source each external fact (open first, cite always)

The "where do I look this up" layer, the sourcing companion to the deciding layers
(`figure_selection.md`, `cocrystal_id.md`). It exists because analytical chemistry's reference
data is **mostly closed** — the genomics-style world of dozens of clean open APIs doesn't apply
here. So this doc records, per fact type: the best OPEN source, whether it has an API, and what
stays manual/paywalled. Ground rule (SKILL.md "When NOT to over-reach"): **external facts are
fetched + cited with source and date, never recalled** — and prefer an *experimental* value to a
computed/predicted one, saying which it is.

## Wired in code (`scripts/sources.py`, stdlib-only, best-effort, cites the URL)

- **COD — Crystallography Open Database** (`cod_search`, `cod_fetch_cif`). Free CIFs, no key.
  Search by Hill-notation formula (`"C9 H8 O4"`), `text=`, elements, cell, spacegroup, `doi`;
  fetch `.cif` by COD id → straight into `crystal_engine.load` → `crystal_pxrd`/`realistic_pattern`/
  `cocrystal.phase_report`. **Live-verified end-to-end.**
  - **HARD RULE — a formula match is NOT an identity match.** A COD `formula=` search returns *every*
    compound with that formula: polymorphs, redeterminations, isomers, **esters**, salts, and unrelated
    molecules (`C5 H8 O4` = glutaric acid **or** methyl-ethyl oxalate; `C13 H18 N4 O6` = a caffeine–glutaric
    cocrystal **or** a pyrazole-pyridine oxide). **Never compute, plot, or cite a hit until you have read
    its `_chemical_name` / `_publ_section_title`.** Use `cod_verify(id, expect)` / `cod_identity(id)`, or
    do search + name-filter in one call with **`cod_search_verified(name_contains=…, formula=…)`**. Every
    `cod_fetch_cif` return now carries an `identity` dict, and a formula `cod_search` carries a `warning`.
    *(This rule exists because it was violated — a whole caffeine-cocrystal reference library was built on
    formula collisions, 2026-07-04.)*
- **PubChem PUG-REST** (`pubchem_properties`). Free, ~5 req/s, no key. Identity + **computed**
  properties (MolecularFormula, MolecularWeight, XLogP, InChIKey…). *Two cautions:* (1) SMILES was
  renamed (`CanonicalSMILES` → `ConnectivitySMILES`/`SMILES`) and one bad property name 400s the
  whole call, so it's omitted from the defaults — request it explicitly with the current name.
  (2) **pKa is NOT reliably in PubChem** — for a ΔpKa input use an experimental source below.

## Free, but no clean API → dataset download or manual (cite the source + date)

- **Experimental pKa (for `cocrystal.classify_ionisation` ΔpKa):** **DataWarrior** bundles an open
  experimental-pKa set (~7,900 compounds, strongest acidic/basic in water) — the best free
  *experimental* pKa. Fallback **predicted**: **OPERA / EPA CompTox Chemicals Dashboard** (has an
  API) — say "predicted" when you use it.
- **Aqueous solubility / preformulation:** **AqSolDB** (curated, ~10k compounds, open dataset).
- **FTIR / Raman reference spectra (band assignments):** **SDBS** (AIST — IR/NMR/MS/Raman of
  organics; aspirin, salicylic acid, paracetamol all present), **NIST Chemistry WebBook** (evaluated
  IR), **SpectraBase** (Bio-Rad/Wiley, large free portal). All **web-only / no bulk API** and terms
  restrict scraping — use to *confirm a band assignment you quote*, manually.
- **Literature / DOIs (ties to the EndNote/.ris flow):** **Europe PMC** and **CrossRef** — free
  REST APIs (not wired here; the WebSearch/WebFetch tools already cover ad-hoc lookups).
- **Drug-product formulation / substances:** **DailyMed** and **FDA GSRS** — free APIs (excipients,
  substance records).

## Paywalled → cite only, do not try to "connect"

- **CSD / CCDC** — the definitive organic/**cocrystal** structure database; licensed (the CSD
  Python API needs a subscription). COD is the open substitute for what it has.
- **ICDD PDF (Powder Diffraction File)** — the standard reference powder-pattern DB for phase ID;
  licensed. **The open path is what this skill already does:** calc PXRD from COD CIFs
  (`realistic_pattern`) + `cocrystal.phase_report`/`ddsimca_fit` matching — plus free Rietveld tools
  (Profex/BGMN, FullProf, Jana) if you need full refinement.
- **Pharmacopoeias (USP, Ph. Eur., BP)** — monographs, impurity limits, official methods. No API;
  reference-standard *prices/availability* are on the USP store (fetch + cite, per the earlier
  procurement pass).

## Which source for which question (quick route)

| You need… | Source | Access |
|---|---|---|
| a CIF (parent / target / polymorph) | **COD** | `sources.cod_search`/`cod_fetch_cif` (wired) |
| molecular identity / MW / logP | **PubChem** | `sources.pubchem_properties` (wired; computed) |
| experimental **pKa** (ΔpKa) | **DataWarrior** set (or OPERA=predicted) | manual dataset |
| an FTIR band assignment | **SDBS / NIST / SpectraBase** | manual, cite |
| a reference **powder pattern** | COD CIF → `realistic_pattern` | wired (the open ICDD path) |
| a cocrystal crystal structure | **CSD/CCDC** (licensed) or COD if present | manual/paywalled |
| an impurity limit / official method | **USP / Ph. Eur.** monograph | paywalled, cite |

**The domain reality (worth stating in a report's method section):** crystallography is the open
corner — COD + free tools are a complete stack, built as the answer to the paywalled CSD/ICDD. The
*measurement* reference data (spectral libraries, monographs) stayed closed. So for
cocrystal PXRD work you are on fully open ground; for spectral-library confirmation you are not.
