# cocrystal_id.md — decide WHAT evidence identifies a cocrystal, before running anything

## Contents

- The four outcomes it must tell apart
- How to use it (question → evidence → route)
- Evidence matrix (technique → what it *uniquely* answers → lab needs → route)
- Decision logic (the truth table)
- Same sample, different question → different measurement
- Pitfalls (the "never", with the reason)
- Lab-capability gate (the capability→requirement matrix, cocrystal edition)
- Routing into the skill (which function)
- What this doc does NOT do
- Sources

The "start here" for **cocrystal identification**. Sibling of
`figure_selection.md`: same job, one axis over. `figure_selection.md` answers *"what
should this figure be?"*; this answers *"which measurement actually settles the
question — is this a cocrystal?"* — so you commit an instrument, a coformer, and a week
of screening to the evidence that discriminates, not to the one that merely looks busy.

**Why this exists.** A capable model already reasons about cocrystal evidence ad hoc —
the value is **consistency**: one grounded, repeatable route instead of a slightly
different call each session. This is **form, not content** (SKILL.md "Form is the
skill's job; judgment is the model's"): it gives the shared framework + the citable rule
behind each call. It does **not** decide your coformer, read your spectrum, or write the
conclusion — that stays yours.

**It routes, it doesn't re-implement.** Once you've chosen, the *how* lives in the
crystal family (`crystal_engine` / `crystal_pxrd` / `crystal_view`), the FTIR family
(`spectra`), the cocrystal arithmetic (`cocrystal`: ΔpKa calculator, NNLS sum-of-parents),
and the gates in `verification.md`.

---

## The four outcomes it must tell apart

Mixing an API with a coformer and getting "a new solid" is **not** proof of a cocrystal.
Five things can happen; the method's whole job is discriminating them [Aitipamula 2012]:

| Outcome | What it is | The tell |
|---|---|---|
| **Physical mixture** | two crystalline phases coexisting, unreacted | PXRD = superposition of the two parents |
| **Cocrystal** | multi-component crystal, **all components neutral**, defined stoichiometry, all solid at RT | new PXRD phase **+** no proton transfer (neutral donor/acceptor bands only *shifted*; ΔpKa < 0) |
| **Salt** | proton transferred acid→base, ionic | new PXRD phase **+** proton transfer (a new *ionised-group* fingerprint — carboxylate / sulfonate / ammonium…; ΔpKa > ~3) |
| **Solvate / hydrate** | one component is a liquid (solvent/water) at RT | new phase + weight loss / solvent bands (TGA, ~1640/3400 for water) |
| **Polymorph** | **single** component, new packing | new PXRD phase but only ONE component present |

The **salt↔cocrystal boundary is a continuum**, not a switch [Childs 2007]. Report where
on it you land and the evidence, not a bare label.

---

## How to use it (question → evidence → route)

1. **Did a new phase form at all?** → PXRD: new peaks vs sum-of-parents. If it's just the
   two parents' peaks, stop — it's a physical mixture, nothing formed.
2. **Neutral (cocrystal) or ionised (salt)?** → the proton-transfer evidence: the FTIR
   ionised-group window (whichever group ionises) + the ΔpKa rule (+ ssNMR/SCXRD if it
   lands in the continuum).
3. **Is it the *target* phase, and single/pure?** → PXRD match to the calculated pattern
   of the intended structure; class/one-class discrimination.

Then read the **evidence matrix**, apply the **decision logic**, check the **pitfalls**,
and route into the family module.

---

## Evidence matrix (technique → what it *uniquely* answers → lab needs → route)

| Technique | What it uniquely settles | Needs | Route |
|---|---|---|---|
| **PXRD** | *New phase vs physical mixture* (positions), phase match to a target | a diffractometer **or** a supplied/experimental pattern | `crystal_pxrd.plot_overlay_patterns` (calc API vs coformer vs cocrystal + experimental); `peak_table`; **`cocrystal.phase_report`** (NNLS sum-of-parents + Rwp + new/lost-peak table) |
| **FTIR / Raman** | *Salt vs cocrystal* via proton transfer — the discriminator this skill is built for | an FTIR (benchtop ATR suffices) | `spectra` overlay/waterfall zoomed to **the window of whichever group ionises** (carbonyl/carboxylate, sulfonate, ammonium…) |
| **ΔpKa rule** | *Predicted* protonation state — orthogonal, zero-cost, run it first | just the two pKa values | `cocrystal.classify_ionisation(pka_acid, pka_base_conjugate)` → zone + caption (`ΔpKa = pKa(baseH⁺) − pKa(acid)`) |
| **DSC / hot-stage** | *Single new melt vs eutectic* — a fast new-phase screen | a DSC | `charts` (thermogram); not a crystal-family job |
| **ssNMR (¹⁵N/¹³C CP-MAS)** | *Definitive* proton location (salt vs cocrystal) | solid-state NMR (rare) | — (capability gap; name it) |
| **SCXRD** | *Gold standard*: proton, stoichiometry, packing | a single crystal + diffractometer | `crystal_engine` validation → `crystal_view` |
| **Solution NMR** | *Stoichiometry only* — the cocrystal **dissolves**, so it can't prove solid-state association | any NMR | integration ratio only; see pitfalls |

---

## Decision logic (the truth table)

- PXRD ≈ **sum of parents** → **physical mixture** (nothing formed). Stop.
- new phase **+** FTIR shows the **neutral** donor/acceptor bands only *shifted / broadened* (no new ionised-group band) **+** ΔpKa < 0 → **cocrystal**.
- new phase **+** FTIR shows a **new ionised-group fingerprint** (the protonation state changed — pick the band from the synthon table below) **+** ΔpKa > ~3 → **salt**.
- new phase **+** FTIR ambiguous / partial **+** 0 < ΔpKa < 3 → **salt–cocrystal continuum** → escalate to ssNMR / SCXRD; report the ambiguity, don't force a label [Childs 2007; Cruz-Cabeza 2012].
- new phase but only **one** component present → **polymorph** (wrong problem — that's a screen, not a cocrystal).

**FTIR is the primary, widely-accessible experimental discriminator** for salt vs cocrystal
because it reports the proton-transfer state directly — but *which bands are diagnostic
depends on which group changes protonation state*, so **identify the acid/base pair (the
synthon) first, then read the right window.** The universal rule:

- **Salt** → a **new ionised-group fingerprint appears** (protonation state changed).
- **Cocrystal** → the **neutral** donor/acceptor bands only **shift / broaden** (H-bonding), with *no* ionised-group band.

Which band, by the group that ionises (verify the exact window against a reference before quoting — [Socrates 2001]):

| Who ionises | Neutral (cocrystal: shifted) | Ionised (salt: new bands) |
|---|---|---|
| **Carboxylic acid** → carboxylate (the common case) | C=O ~1700–1730 | C=O lost; **COO⁻ asym ~1550–1650 + sym ~1300–1420** |
| **Sulfonic / sulfuric** → sulfonate/sulfate | S=O ~1350 / 1150 | shifted **SO₃⁻/SO₄ ~1200 / 1040** |
| **Basic amine or pyridine** (API is the base; acid = HCl, H₂SO₄, H₃PO₄, mesylate…) → ammonium / pyridinium | N–H sharp; ring bands | broad **N⁺–H ~2000–3000** + ammonium bends ~1600 / 1500; pyridinium ring shift |
| **Phenol / alcohol** (weak; rarely ionises) | O–H shift/broaden | phenolate C–O shift (uncommon) |

The third row is the trap the carboxyl assumption misses: a **basic API + a non-carboxylic
acid** (HCl, sulfate, phosphate, mesylate) is a very common salt with **no carboxyl at all** —
the tell is amine→ammonium, and a mineral-acid counterion is IR-silent. **FTIR-quiet synthons
exist too:** halogen-bond cocrystals (C–X···N/O, no O–H/N–H change) and salts whose only shift
hides in broad ammonium bands are weak by IR → lean on **ΔpKa + ssNMR/SCXRD** (and PXRD for the
new phase).

---

## Same sample, different question → different measurement

(The `figure_selection.md` "one dataset, four figures" idea, for cocrystal ID.)

| The question | The measurement |
|---|---|
| "Did anything form?" | PXRD: new peaks vs sum-of-parents (+ NNLS/Rwp) |
| "Cocrystal or salt?" | FTIR ionised-group window + ΔpKa (+ ssNMR in the continuum) |
| "What stoichiometry?" | solution NMR integration / SCXRD |
| "Is it the *target* cocrystal, and pure?" | PXRD vs **calc pattern of the target CIF** + one-class fit (`chemometrics.ddsimca_fit`/`ddsimca_predict` + `class_fom`) |
| "Does it survive / convert?" | variable-T or post-slurry PXRD |

Fix the question first; the same powder gives different figures.

---

## Pitfalls (the "never", with the reason)

- **Never call a cocrystal from ONE technique.** PXRD proves *a new phase*, not *cocrystal
  vs salt* — the proton isn't visible in a powder pattern. Pair PXRD (new phase) with FTIR
  + ΔpKa (protonation) at minimum [Aitipamula 2012].
- **Judge PXRD phase by peak POSITIONS, not intensities.** Preferred orientation and
  particle statistics distort intensities; the phase lives in the 2θ positions. (Intensity
  realism — March–Dollase, Kα₁/Kα₂ — is `pxrd_realism`, for *matching*, not ID.)
- **Solution NMR dissolves the association.** It reports stoichiometry, never "it's a
  cocrystal" — the solid-state contact is gone the moment it's in solvent.
- **Amorphous product → a PXRD halo, not sharp new peaks.** Don't read a halo as a new
  crystalline phase; it's disorder, and FTIR/DSC take over.
- **ΔpKa is predictive, not proof.** It ranks salt-risk; 0–3 is explicitly unpredictable
  [Cruz-Cabeza 2012]. Use it to *prioritise*, confirm by FTIR/ssNMR.
- **A physical mixture can hide a little cocrystal (and vice versa).** "New peak present"
  is qualitative; quantify the phase fraction (NNLS sum-of-parents + Rwp) before a purity
  claim.
- **Calc PXRD from a CIF needs the declared wavelength and no hand-shifting.** A small
  calc/exp 2θ offset is zero-point/cell/temperature, not a phase miss — flag it, don't
  slide the pattern to fake a match (`crystal.md`).

---

## Lab-capability gate (the capability→requirement matrix, cocrystal edition)

Before proposing the screen, check each discriminating question against what YOUR lab
actually has, and name the gaps as residual risk — don't promise a discrimination the
instruments can't deliver. This is a template to fill per project, not a fixed inventory:

| Question | Needs | Have it? | If not… |
|---|---|---|---|
| new phase? | PXRD | _(fill in)_ | lean on calc-PXRD vs an *experimental* pattern you can obtain; FTIR + DSC as the new-phase proxy |
| salt vs cocrystal? | FTIR + ΔpKa | _(fill in)_ | if no FTIR, ΔpKa + ssNMR/SCXRD carry it |
| definitive proton location | ssNMR / SCXRD | _(fill in)_ | report the continuum honestly; don't over-claim |
| stoichiometry | NMR / SCXRD | _(fill in)_ | solution NMR integration on the isolated phase |

Two things hold regardless of the lab: **(1)** FTIR is usually the most *accessible*
salt/cocrystal discriminator — a benchtop ATR reads the proton-transfer state directly — so
it's a strong default anchor when it's available. **(2)** The skill supplies **calculated**
PXRD from a CIF (`crystal_pxrd`) regardless of whether you own a diffractometer — so even
PXRD-limited, you can compare a *measured* pattern to calc(API)/calc(coformer)/calc(target
cocrystal) and to a sum-of-parents model.

---

## Routing into the skill (which function)

- **PXRD phase ID / cocrystal-vs-mixture:** `crystal_pxrd.plot_overlay_patterns([(label, cif), …],
  cfg, experimental=…)` — waterfall calc(cocrystal) vs calc(API) vs calc(coformer) with the
  measured pattern at the bottom; `peak_table` / `write_peaks_csv` for the (2θ, d, hkl, I)
  list. Set `cfg.color_by_component=True` on structure views to keep API vs coformer legible.
- **Salt/cocrystal FTIR:** `spectra` overlay/waterfall on the window of whichever group
  ionises (carbonyl/carboxylate, sulfonate, ammonium…); identical processing across samples
  (`verify.same_processing`); the diagnostic is the *presence/absence* of the ionised-group
  fingerprint, stated in the caption.
- **3D structure / H-bond heterosynthon:** `crystal_engine` (validate the CIF first) →
  `crystal_view.render_hbond_environment` for the new O–H···N contact.
- **ΔpKa salt/cocrystal call:** `cocrystal.classify_ionisation(pka_acid, pka_base_conjugate)`
  → `zone` (salt / cocrystal / continuum) + `confident` + a cited caption [Cruz-Cabeza 2012].
  Triage only — the 0–3 continuum is unpredictable, confirm by the FTIR ionised-group band.
- **New phase vs physical mixture:** `cocrystal.phase_report(x, y_obs, [parent_A, parent_B])`
  → NNLS sum-of-parents fit + Rwp + `new_peaks`/`lost_peaks` + a heuristic verdict. Low Rwp,
  no new peaks ⇒ mixture; high Rwp + positive residual peaks ⇒ a new phase is plausible
  (discrimination, not quantitative phase % — confirm structurally). Works on PXRD or FTIR.
- **Is it the target phase? (one-class):** `chemometrics.ddsimca_fit(target_spectra)` →
  `ddsimca_predict(model, unknown)` accepts/rejects an unknown against the target class's
  DD-SIMCA boundary; `chemometrics.class_fom(y_true, y_pred)` reports sensitivity / specificity
  / efficiency. A *rejected* unknown is NOT the target cocrystal (a mixture / other polymorph).
- **Still on the roadmap:** PLS-DA multi-class (cocrystal vs mixture vs starting materials);
  Hirshfeld surfaces; ADDSYM/ADP validation. (PXRD intensity realism — March-Dollase / Kα₁₂ /
  Caglioti — is built: `pxrd_realism`, via `crystal_pxrd.realistic_pattern`.)

---

## What this doc does NOT do

It doesn't pick the coformer, choose the crystallisation method, read your pattern, or
write the conclusion. It standardises the *identification-evidence* decision and hands you
the citable reason. If a case isn't covered, reason from the four outcomes and the decision
logic; don't force a label onto a continuum.

---

## Sources

- Aitipamula S, et al. 2012. *Polymorphs, Salts, and Cocrystals: What's in a Name?* Cryst. Growth Des. 12(5):2147–2152. doi:10.1021/cg3002948
- Childs SL, Stahly GP, Park A. 2007. *The Salt–Cocrystal Continuum: The Influence of Crystal Structure on Ionization State.* Mol. Pharmaceutics 4(3):323–338. doi:10.1021/mp0601345
- Cruz-Cabeza AJ. 2012. *Acid–base crystalline complexes and the pKa rule.* CrystEngComm 14:6362–6365. doi:10.1039/C2CE26055G
- FDA CDER. 2018. *Regulatory Classification of Pharmaceutical Co-Crystals — Guidance for Industry.* (EMA reflection paper on cocrystals, 2015, is the EU counterpart.)
- Socrates G. 2001. *Infrared and Raman Characteristic Group Frequencies*, 3rd ed., Wiley. (carboxylic acid C=O vs carboxylate asym/sym assignments — verify the exact window against your reference before quoting a submission.)
- Crystallographic conventions + calc-PXRD provenance: see `crystal.md` (IUCr Acta Cryst Notes for Authors; ORTEP; the `crystal_pxrd` calc-vs-experimental overlay rules).

Scaffolding shape adapted from `figure_selection.md` (itself after scipilot-figure-skill, MIT).
