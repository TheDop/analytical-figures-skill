#!/usr/bin/env python3
"""
Prompt: "Here is aspirin's CIF (COD 7247819). Validate it, calculate its Cu K-alpha powder pattern,
and give me an ORTEP view with the hydrogen bonds."

Data: aspirin form I, Crystallography Open Database entry 7247819 (CrystEngComm), Mo K-alpha data.
Needs gemmi + Dans_Diffraction (pymatgen optional, for the cross-check).
"""
DESCRIPTION = ("Aspirin (COD 7247819): CIF validation table, calculated Cu K-alpha powder pattern with "
               "an independent cross-check, and an ORTEP view with hydrogen bonds")
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))          # the skill root
sys.path.insert(0, ROOT)

from dataclasses import replace
from scripts.config import Config
from scripts import style, verify, crystal_engine, crystal_pxrd, crystal_view

# ====================================================================== CONFIG
cfg = Config(
    cif_path=os.path.join(ROOT, "skill_validation", "crystal", "cifs", "monoclinic_aspirin__COD7247819.cif"),
    output_dir=os.path.join(HERE, "out"),
    journal="nature", column="single", formats=("pdf", "png"),
    pxrd_wavelength=1.5406,                      # the CIF is Mo data; the user asked for a Cu K-alpha pattern
    pxrd_two_theta_min=5, pxrd_two_theta_max=50,
    pxrd_crosscheck="auto",                      # pymatgen top-peak cross-check when it is installed
    view_style="ellipsoid", adp_probability=0.50, view_hide_ch=True, view_label_atoms="hetero",
    view_orientation="pca",
    strict=False,                                # warn, don't raise: the table reports every gate
    description=DESCRIPTION,
)


def main():
    os.makedirs(cfg.output_dir, exist_ok=True)

    # ---- 1. validate BEFORE anything is drawn
    s = crystal_engine.load(cfg)
    rep = crystal_engine.validate(s, cfg)
    print(crystal_engine.table_text(rep))
    crystal_engine.write_csv(rep, os.path.join(cfg.output_dir, "hbonds.csv"))

    # ---- 2. calculated powder pattern (Cu K-alpha) + peak list
    fig, ax, pat = crystal_pxrd.plot(s, cfg)
    verify.audit_layout(fig, cfg)
    style.save_fig(fig, os.path.join(cfg.output_dir, "pxrd"), cfg)
    rows = crystal_pxrd.peak_table(s, cfg)
    crystal_pxrd.write_peaks_csv(rows, os.path.join(cfg.output_dir, "peaks.csv"))
    top = rows[0] if rows else None
    print(f"  strongest reflection: {top}")

    # ---- 3. ORTEP view WITH the hydrogen bonds: aspirin's only H-bond is to a symmetry mate, so the
    #         H-bond-environment view (asymmetric unit + its H-bonded neighbours) is the one that shows it
    crystal_view.save(crystal_view.render_hbond_environment(s, cfg),
                      os.path.join(cfg.output_dir, "structure_ellipsoid"), cfg)

    caption = (
        f"Figure. Aspirin (COD 7247819). Left: powder pattern calculated at Cu Kα ({pat['wavelength']:.4f} Å) "
        f"from the deposited structure, 5–50° 2θ; the strongest reflection is at {top['two_theta']:.2f}° 2θ "
        f"(d = {top['d']:.3f} Å, hkl {top['hkl']}). Right: the asymmetric unit with its hydrogen-bonded neighbours, displacement "
        f"ellipsoids at the {int(cfg.adp_probability*100)} % probability level, C–H hydrogens omitted, heteroatoms "
        f"labelled; dashed lines are the hydrogen bonds listed in the validation table, symmetry mates "
        f"superscripted with their operator.\n")
    open(os.path.join(HERE, "caption.txt"), "w", encoding="utf-8").write(caption)
    print("  wrote pxrd.pdf/.png, structure_ellipsoid.png, hbonds.csv, peaks.csv, caption.txt")


if __name__ == "__main__":
    main()
