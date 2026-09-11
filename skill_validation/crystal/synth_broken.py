"""
synth_broken.py - make SYNTHETIC broken-input CIFs from known-good ones, to exercise
the density-triple-check branches that no clean fetched CIF hits. Edits are documented
in each output's header. These are NOT real structures.

Cases 1-2 derive from ELAINM (CCDC 783049), which is not vendored: they run only when
cifs/783049.cif has been fetched (see README). Cases 3-4 derive from COD files and always run.
"""
import os, re
HERE = os.path.dirname(os.path.abspath(__file__))
CIF_DIR = os.path.join(HERE, "cifs")

def read(p): return open(os.path.join(CIF_DIR, p), encoding="utf-8", errors="replace").read()

def write(p, s, note):
    hdr = (f"# ============================================================\n"
           f"# SYNTHETIC TEST FILE - {note}\n"
           f"# Edited from a real CIF by synth_broken.py; NOT a real structure.\n"
           f"# ============================================================\n")
    s = re.sub(r"(?m)^(data_\S+)", hdr + r"\1", s, count=1)
    open(os.path.join(CIF_DIR, p), "w", encoding="utf-8").write(s)
    print("wrote", p)

if os.path.exists(os.path.join(CIF_DIR, "783049.cif")):
    # 1) WRONG Z: ELAINM Z 1 -> 2. (ii)=checkCIF formula doubles; (i) from atoms and (iii)
    #    declared unchanged  => expect  A-ALERT (ii != iii).
    s = read("783049.cif")
    s2 = re.sub(r"(_cell_formula_units_Z\s+)1\b", r"\g<1>2", s, count=1)
    assert s2 != s, "Z substitution failed"
    write("broken_wrongZ__ELAINM.cif", s2, "Z falsified 1->2  (expect A-ALERT: ii != iii)")

    # 2) SQUEEZE-LIKE: drop the disordered water ATOM_SITE rows but keep _chemical_formula_sum
    #    and declared density (the SQUEEZE signature). (i) from atoms < (iii); (ii)~(iii)
    #    => expect BANNER (incomplete atom list).  Only atom_site rows removed (type_symbol
    #    is a bare element in col 2; aniso/geom rows for the same labels are left harmless).
    water = {"O21", "O22", "O23", "H231", "H232"}
    ELS = {"O", "H", "C", "N", "S"}
    out = []
    for ln in read("783049.cif").splitlines():
        t = ln.split()
        if t and t[0] in water and len(t) >= 6 and t[1] in ELS:   # an atom_site data row
            continue
        out.append(ln)
    write("broken_squeeze__ELAINM.cif", "\n".join(out),
          "disordered water atom_site rows removed, formula/density kept (expect BANNER i<iii)")
else:
    print("cifs/783049.cif (ELAINM, CCDC) absent -> skipping cases 1-2 (see README)")

# 3) H-FREE: drop all H atom_site rows from aspirin, keep formula C9H8O4. (i) underestimates
#    by the H mass => expect BANNER (i<iii).
out = []
for ln in read("monoclinic_aspirin__COD7247819.cif").splitlines():
    t = ln.split()
    if t and len(t) >= 6 and t[1] == "H" and t[0][:1] == "H":   # an atom_site H row
        continue
    out.append(ln)
write("broken_Hfree__aspirin.cif", "\n".join(out),
      "all H atom_site rows removed, formula kept (expect BANNER i<iii)")

# 4) MULTI-BLOCK: concatenate two single-block CIFs into one file with TWO data_ blocks
#    (tests the block-selection gate; a file-format case, not a broken structure).
ua = read("specialpos_urea__COD1008785.cif").rstrip()
ap = read("monoclinic_aspirin__COD7247819.cif").lstrip()
mb = ("#####################  SYNTHETIC MULTI-BLOCK CIF  #####################\n"
      "# Two independent structures (urea + aspirin) concatenated into one file.\n"
      "# An engine must DETECT both data_ blocks and refuse to silently use one.\n"
      "######################################################################\n\n"
      + ua + "\n\n" + ap + "\n")
open(os.path.join(CIF_DIR, "multiblock__urea-aspirin.cif"), "w", encoding="utf-8").write(mb)
print("wrote multiblock__urea-aspirin.cif")
