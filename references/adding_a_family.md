# adding_a_family.md — extending the skill to a new figure type

To add a family (e.g. chromatography, DSC, UV-Vis):

1. **Add only the necessary knobs to `Config`** (units in comments). This
   preserves the single-CONFIG-block rule.

2. **Write a module `scripts/<family>.py`** exposing plain functions that take
   `cfg`. Reuse the shared infrastructure:
   - `style.figure(cfg)` / `style.save_fig(fig, base, cfg)` for sizing + export,
   - `verify.check_*` to gate inputs,
   - `verify.summarize` for any error bars,
   - `verify.render_preview` + `audit_layout` + `READ_IMAGE_CHECKLIST` for QA.

3. **Encode the conventions, not the matplotlib.** The reference doc should say
   what's *right* for this figure type (axis directions, what must be labelled,
   which chart is wrong for which data) — assume the plotting API is known.

4. **Add a row to the SKILL.md "Files" list and a `references/<family>.md`.**

5. **Add a gate if the family has a correctness trap.** Calibration refuses
   extrapolation and forces residuals; charts refuse low-n mean bars; spectra
   assert identical processing. A new family should ship with whatever the
   equivalent "easy to get silently wrong" guard is.

## Test before documenting

Run the new module against real (or realistic synthetic) data and *look at the
PNG* before writing the docs. Conventions specified against data you haven't
seen are how wrong assumptions get baked in — this is why PXRD profiles were
validated on synthetic patterns first and should be re-checked on real data.
