# Crystal-family validation set — edge-case matrix

Purpose: a set of CIFs that exercises **every gate and convention** in the crystal
family (`scripts/crystal_engine.py`, `crystal_pxrd.py`, `crystal_view.py`), assembled
before the family was documented ("test before document"). Each row names the feature
under test, the spec section it came from, the **failure mode if the engine gets it
silently wrong**, a candidate real structure, and how it was supplied.

**Status legend**
- ✅ **covered** by the in-hand CIF `783049.cif` (ELAINM) — no fetch needed.
- 🔵 **FETCH** — a real CIF pulled from CSD (by refcode, via Access Structures) or COD (open).
- 🟠 **SYNTHESIZE** — made by `synth_broken.py` editing a known-good CIF (the honest way
  to test "broken input" branches; you can't reliably find published-wrong CIFs).
- 🔁 **re-run** — needs no new CIF, just a fixed probe on an existing one.

---

## What ELAINM (`783049.cif`, CCDC — optional, not vendored) covers ✅

Verified 2026-06-20:

- Triclinic **P-1**, inversion centre, Z=1 (low-symmetry expansion). Density (i)
  from symmetry-expanded atoms = 1.6441 vs declared 1.644 — expansion is correct.
- **Fractional formula** from partial occupancy (`C26 H19.32 N4 O11.68`) — element
  counts reproduced exactly; the occ-weighted comparison works.
- **Positional disorder WITHOUT disorder-group labels** (3 water O sites 0.56–0.77 Å
  apart) — drove the disorder gate; naive bonding invents 2 bogus O–O bonds.
- Density triple-check **all-consistent PASS** branch (solvent is modelled).
- **`_geom_hbond` oracle present** (5 bonds) incl. a **bifurcated** donor and a
  **sub-120° (112.8°)** contact that the 120° floor demotes.
- **X–H normalization** effect (X-ray H at 0.84/0.90 Å → neutron) — D···A invariant,
  H···A/angle change as expected.
- **Cu Kα weighted λ = 1.54178** (not Kα₁ 1.5406) — confirms use-declared-λ.
- Planar + H-bonded **good-case** for PCA orientation.

---

## Status — fetched / built / run (2026-06-20)

All validated by `crystal_validate.py` (7 real + 3 synthetic). Files marked **CCDC** are not
redistributable and are **not vendored**: fetch them from CCDC Access Structures by deposition
number into `cifs/` and every script picks them up automatically (they are git-ignored).

| file | role | rows | status |
|---|---|---|---|
| `783049.cif` **(CCDC 783049 — fetch)** | ELAINM (Cu) — baseline | — | ✅ PASS; disorder bonds suppressed; D···A 2/2 |
| `monoclinic_aspirin__COD7247819.cif` | P2₁/c, Mo | #1, #11 | ✅ PASS |
| `specialpos_urea__COD1008785.cif` | special pos + new symop tag + neutron + charged `O2-` | #2, #14, #24 | ✅ PASS (orbits match Wyckoff) |
| `chiral_L-alanine__COD1574526.cif` | chiral P2₁2₁2₁ + synchrotron λ | #7 | ✅ FLAGGED a real formula_weight error |
| `nohbond_naphthalene__COD2311088.cif` | π-stack / no H-bond, Mo | #17 | ✅ PASS, planar (RMS-oop 0.002 Å) |
| `lactose__COD2206486.cif` | chiral sugar, Mo | #22 | ✅ PASS; D···A 2/2 |
| `sulfur_ptu-ellagic__CCDC2018928.cif` **(CCDC 2018928 — fetch)** | the real PTU cocrystal — **S** | #8(partial) | ✅ PASS; D···A 2/2 |
| `1178952.cif` (HMT) **(CCDC 1178952 — fetch)** | **globular** + **I-centred** + special pos; CSD molecule-completed list | #1(globular), #3(I) | ✅ PASS (i=1.345); w2/w1=1.00 degeneracy flag fires; exposed the global-dedup bug |
| ~~`1154856.cif` (FEROCE)~~ | intended disorder+metal | #9, #8 | ❌ REJECTED — "no coordinates deposited" (1952 cell-only entry) |
| `element-Fe_ferrocene__CCDC131688.cif` **(CCDC 131688 — fetch)** | **Fe** (non-table element) + special pos + neutron; **no symop loop, no type_symbol col** | #8 | ✅ closes #8; forced symop-from-name + column-wise/element-from-label parser fixes; ⚠️ ring disorder has no occupancies → density unreliable (no decl density to gate), NOT the #9 fixture |
| `disorder-groups_flufenamic__COD4118081.cif` | **disorder WITH `_atom_site_disorder_group`** (CF₃ A/B) + F + Z′=3 + H-bonds | #9, #8 | ✅ **CLOSES #9**; tag=18 atoms, 5 mutually-exclusive F pairs; PASS (i=1.458); D···A 5/5 |
| `zprime9-densflag_flufenamic__COD4118079.cif` | extreme **Z′=9**; declared-density anomaly | bonus | ✅ FLAG: (i)=(ii)=1.486 vs declared 1.569 — a real metadata error caught |
| `multiblock__urea-aspirin.cif` | synthetic 2-`data_`-block file | block gate | ✅ flags "MULTI-BLOCK (2 blocks)" |
| `broken_wrongZ__ELAINM.cif` (needs ELAINM; `synth_broken.py`) | A-alert (ii≠iii) | #6 | ✅ fires |
| `broken_squeeze__ELAINM.cif` (needs ELAINM; `synth_broken.py`) | incomplete-atoms / SQUEEZE banner | #4 | ✅ fires |
| `broken_Hfree__aspirin.cif` | H-free banner | #5 | ✅ fires |

Wavelengths covered: **Cu / Mo / synchrotron / neutron**. PXRD path fixed (`dans_fix.py`):
47 clean Bragg peaks, Dans vs pymatgen 0.011°.

**Closed:** globular degeneracy + I-centred lattice ✅ (HMT 1178952); non-table element ✅
(Fe — ferrocene 131688; F — flufenamic); **disorder WITH `_atom_site_disorder_group` ✅
(flufenamic COD 4118081 — CF₃ A/B groups, also Z′=3)**; multi-block ✅ (synthetic
`multiblock__urea-aspirin.cif`). FEROCE 1154856 ❌ rejected (coordinate-less; deleted).

**Still open — optional only:**
- **#3 a C-centred lattice** (C2/c) — nice-to-have; HMT already covers I-centring. Everything
  else in P0/P1 is now exercised.

**Lesson:** `_atom_site_disorder_group` can't be searched directly, but COD *does* carry it
— **fetch-and-verify** (the harness reports the tag) beats trying to filter on it. Always
confirm an entry has **deposited coordinates** first (the FEROCE trap).

---

## P0 — must have (the gaps ELAINM leaves in the NEW spec logic)

| # | Feature / gate | Spec | Failure if unhandled | Candidate structure | Who |
|---|---|---|---|---|---|
| 1 | **Monoclinic P2₁/c** (commonest SG) | §4, §5 | wrong multiplicity in expansion → density/F000 off | **Aspirin form I** (refcode ACSALA), P2₁/c | 🔵 |
| 2 | **Atom ON a special position** (Z′<1) | §4, §5 | special-position occupancy not halved → density/F000 **double-count** | **Urea** (P-4̄2₁m, molecule on mm2) or **hexamethylenetetramine** (cubic, on special pos) | 🔵 |
| 3 | **Centred lattice** (C/I/F/R) | §4 | centering translations dropped → expansion misses half the cell | a **C2/c** or **I**-centred small molecule (e.g. paracetamol form I, monoclinic) | 🔵 |
| 4 | **SQUEEZE'd / omitted solvent** | §5 | should hit **"incomplete list → banner"**, must NOT hard-fail | a `_platon_squeeze`-flagged structure (common in cage/MOF papers) | 🔵 |
| 5 | **H-free refinement** (H not located) | §5 | density (i) underestimates by H mass → **banner branch** | an older / heavy-atom structure with no H sites | 🔵 |
| 6 | **Internally inconsistent CIF** (wrong Z / formula typo) | §5 | (ii)≠(iii) → **A-alert banner** | corrupt a good CIF (bump Z, wrong `_chemical_formula_weight`) | 🟠 |
| 7 | **Chiral / Sohncke** (P2₁2₁2₁ or P2₁) | §7 | eigenvector **sign** flips → wrong **enantiomorph** rendered | **L-alanine** (P2₁2₁2₁) or **sucrose** (P2₁) | 🔵 |
| 8 | **Element beyond the vdW table** (F/Br/I/P/metal) | §5 Block C, §8 | radius lookup KeyError or contact silently dropped | a halogenated drug (F/Cl/Br) or a phosphate | 🔵 |
| 9 | **Disorder WITH A/B group labels** | §4 | must honour `_atom_site_disorder_group` (ELAINM had none) | any structure with refined A/B occupancy disorder | 🔵 |
| 10 | **Globular molecule** (no plane) | §7 | w2/w1 ≥ 0.85 → must warn, view arbitrary | **adamantane** or a near-spherical cage | 🔵 |
| 11 | **Mo Kα data** (λ 0.71073) | §5, §6 | λ/keV handling assumes Cu | any Mo-collected structure (most modern small-mol) | 🔵 |
| 12 | **PXRD (000) leak + API** | §6 | FOUND broken: masked pattern had zero Bragg peaks | re-run fixed probe on ELAINM (+ a heavy-atom check) | 🔁 |

---

## P1 — should have

| # | Feature / gate | Spec | Failure if unhandled | Candidate structure | Who |
|---|---|---|---|---|---|
| 13 | **Z′ > 1** (≥2 independent molecules) | §4, §7 | per-molecule completion / orientation | a Z′=2 cocrystal | 🔵 |
| 14 | **New symop tag** `_space_group_symop_operation_xyz` | §4 | parser only reads the old `_symmetry_equiv_pos_as_xyz` | any post-~2010 CIF (most use the new tag) | 🔵 |
| 15 | **Very short strong O–H···O** (~2.4 Å) | §5 Block C | Jeffrey "strong" classification edge | an acid salt / co-crystal with a short H-bond | 🔵 |
| 16 | **Weak C–H···O dominated** (opt-in weak donors) | §5 Block C | weak-donor + angle-floor coupling | an aromatic cocrystal with only weak contacts | 🔵 |
| 17 | **No classical H-bonds / π-stacking** | §7 | long-axis roll fallback (not the H-bond objective) | **naphthalene** / **pyrene** | 🔵 |
| 18 | **Fully isotropic-only CIF** (no aniso U) | §1 Phase 2, §9 | ADP request must degrade to ball-and-stick + note | a powder-derived or old Uiso-only structure | 🔵 |
| 19 | **Roll-degenerate: planar but H-bonds ⊥ plane** | §7 | second degeneracy flag (objective curvature, not w2/w1) | curate from candidates (rarer) | 🔵 |
| 20 | **Substitutional/site-sharing disorder** (occ sum→1) | §4, §5 | element bookkeeping in density/formula | a solid solution / mixed-occupancy site | 🔵 |

---

## P2 — nice to have / cocrystal-ID relevance

| # | Feature / gate | Spec | Why | Candidate structure | Who |
|---|---|---|---|---|---|
| 21 | **Aspirin polymorph II** | §6 | PXRD polymorph discrimination | aspirin form II | 🔵 |
| 22 | **α- vs β-lactose** | §6 | excipient phase ID | lactose (α and β) — in COD | 🔵 |
| 23 | **Real aspirin/paracetamol cocrystal** | all | the cocrystal-ID use case | a published API cocrystal | 🔵 |
| 24 | **Neutron structure** (no Kα λ) | §5 | wavelength validator must accept declared neutron | any neutron determination | 🔵 |
| 25 | **Halogen-bonded cocrystal** (C–I···O) | §5 Block C | acceptor/contact handling beyond H-bonds | a halogen-bond cocrystal | 🔵 |
| 26 | **Non-positive-definite ADP** | Phase 2 | ellipsoid undrawable → must catch | a flagged poor-quality structure | 🔵 |

---

## Running

Everything here needs `gemmi`; the PXRD checks also need `Dans_Diffraction` (and `pymatgen` for
the cross-check). Outputs go to `out/` (git-ignored).

- `python run.py` — the gate, bridged into `pytest tests/` (skipped when the deps are absent):
  `test_engine.py` (the engine's validation table for every CIF in `cifs/`), `crystal_validate.py`
  (an independent gemmi-only re-implementation of the gates — a second opinion), `test_adp.py`
  (U_eq from the ADP tensor vs the CIF's U_iso), `test_realism_cif.py` (Dans reflection list →
  realistic Cu-Kα pattern; positions must match `calc_pattern`), `dans_fix.py` (the two
  Dans_Diffraction API gotchas + the pymatgen top-peak cross-check).
- Renders to eyeball: `test_view.py`, `test_ellipsoid.py`, `test_pxrd.py`, `aspirin_compare.py`,
  `test_pyvista.py` (needs `pyvista`), and `crystal_demo.py` — the full worked example
  (validation table + CSV, PXRD, six structure views, a three-phase overlay, peak list).
- `synth_broken.py` regenerates the synthetic broken-input CIFs from the real ones.

To extend the set, drop a CIF into `cifs/` named `<feature>__<COD-or-refcode>.cif`; every
script globs the folder.
