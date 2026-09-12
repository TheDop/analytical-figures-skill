#!/usr/bin/env python3
"""
crystal_demo.py — worked crystal-family analysis (copy + edit the CONFIG block).
Produces the three deliverables: validation table (+CSV), calculated PXRD, 3D structure.
Run against the skill (from the skill root), or bundle.py it into a standalone.
"""
import os, sys
from dataclasses import replace
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))          # the skill root

from scripts.config import Config
from scripts import style, verify, crystal_engine, crystal_pxrd, crystal_view

# ====================================================================== CONFIG
cfg = Config(
    cif_path=os.path.join(HERE, "cifs", "monoclinic_aspirin__COD7247819.cif"),
    output_dir=os.path.join(HERE, "out", "crystal_demo"),
    journal="nature", column="single", formats=("pdf", "png"),
    pxrd_crosscheck="auto",
    view_orientation="pca", view_hide_ch=True, view_label_atoms="hetero",
    strict=False,
)
os.makedirs(cfg.output_dir, exist_ok=True)


def main():
    s = crystal_engine.load(cfg)

    # --- 1. validation table (CORE) -> text + CSV
    rep = crystal_engine.validate(s, cfg)
    print(crystal_engine.table_text(rep))
    crystal_engine.write_csv(rep, os.path.join(cfg.output_dir, "hbonds.csv"))

    # --- 2. calculated PXRD -> visual QA loop -> vector master
    figp, axp, pat = crystal_pxrd.plot(s, cfg)
    verify.render_preview(figp, os.path.join(cfg.output_dir, "_pxrd_preview.png"))
    verify.audit_layout(figp, cfg)
    style.save_fig(figp, os.path.join(cfg.output_dir, "pxrd"), cfg)

    # --- 3a. 3D structure, ball-and-stick
    crystal_view.save(crystal_view.render(s, cfg), os.path.join(cfg.output_dir, "structure"), cfg)
    print("  structure ->", os.path.join(cfg.output_dir, "structure.png"))

    # --- 3b. ORTEP displacement ellipsoids (Phase 2; 50% probability)
    ecfg = replace(cfg, view_style="ellipsoid", adp_probability=0.50)
    crystal_view.save(crystal_view.render(s, ecfg), os.path.join(cfg.output_dir, "structure_ellipsoid"), ecfg)
    print("  ellipsoids ->", os.path.join(cfg.output_dir, "structure_ellipsoid.png"))

    # --- 3c. unit cell contents + cell box + a/b/c
    ucfg = replace(cfg, view_label_atoms="none")
    crystal_view.save(crystal_view.render_unit_cell(s, ucfg), os.path.join(cfg.output_dir, "unit_cell"), ucfg)
    print("  unit_cell ->", os.path.join(cfg.output_dir, "unit_cell.png"))

    # --- 3d. packing diagram with H-bonds (one layer down the a-axis, cell box + a/b/c)
    pcfg = replace(cfg, pack_cells=(1, 2, 2), view_orientation="axis_a", view_label_atoms="none")
    crystal_view.save(crystal_view.render_packing(s, pcfg), os.path.join(cfg.output_dir, "packing"), pcfg)
    print("  packing ->", os.path.join(cfg.output_dir, "packing.png"))

    # --- 3e. asymmetric unit + H-bonded symmetry neighbours (default 'stub' extent + symmetry superscripts)
    hcfg = replace(cfg, view_label_atoms="hetero")
    crystal_view.save(crystal_view.render_hbond_environment(s, hcfg), os.path.join(cfg.output_dir, "hbond_environment"), hcfg)
    print("  hbond_environment (stub) ->", os.path.join(cfg.output_dir, "hbond_environment.png"))

    # --- 3f. same view, full neighbour molecules (the packing-shell variant)
    hwcfg = replace(cfg, view_label_atoms="hetero", hbond_neighbour="whole")
    crystal_view.save(crystal_view.render_hbond_environment(s, hwcfg), os.path.join(cfg.output_dir, "hbond_environment_whole"), hwcfg)
    print("  hbond_environment (whole) ->", os.path.join(cfg.output_dir, "hbond_environment_whole.png"))

    # --- 3g. clipped unit cell (cell contents cut at the box, not whole molecules)
    clcfg = replace(cfg, view_label_atoms="none", cell_fill="clip")
    crystal_view.save(crystal_view.render_unit_cell(s, clcfg), os.path.join(cfg.output_dir, "unit_cell_clip"), clcfg)
    print("  unit_cell_clip ->", os.path.join(cfg.output_dir, "unit_cell_clip.png"))

    # --- 4. cocrystal-ID: stacked calculated PXRD (phase discrimination) + peak list
    here = os.path.dirname(cfg.cif_path)
    entries = [("aspirin", os.path.join(here, "monoclinic_aspirin__COD7247819.cif")),
               ("urea", os.path.join(here, "specialpos_urea__COD1008785.cif")),
               ("lactose", os.path.join(here, "lactose__COD2206486.cif"))]
    figo, axo, pats = crystal_pxrd.plot_overlay_patterns(entries, cfg)
    style.save_fig(figo, os.path.join(cfg.output_dir, "pxrd_overlay"), cfg)
    print("  pxrd_overlay ->", os.path.join(cfg.output_dir, "pxrd_overlay.pdf"))

    rows = crystal_pxrd.peak_table(s, cfg)
    crystal_pxrd.write_peaks_csv(rows, os.path.join(cfg.output_dir, "peaks.csv"))
    print("  peaks.csv (top peak) ->", rows[0] if rows else "none")

    # --- 5. components distinguished by colour (a no-op on a single-component structure; shown for the API)
    cc = replace(cfg, color_by_component=True, view_label_atoms="none")
    crystal_view.save(crystal_view.render(s, cc), os.path.join(cfg.output_dir, "structure_components"), cc)
    print("  structure_components ->", os.path.join(cfg.output_dir, "structure_components.png"))

    print("\n  crystal demo OK ->", cfg.output_dir)


if __name__ == "__main__":
    main()
