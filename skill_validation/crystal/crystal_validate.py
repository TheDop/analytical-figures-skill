"""
crystal_validate.py - independent multi-CIF gate harness (gemmi only; does NOT use crystal_engine).
Re-implements the structure-agnostic gates -- density triple-check, special positions, disorder,
_geom_hbond oracle, PCA orientation -- as a second opinion against crystal_engine. Runs over every
*.cif in cifs/ and prints a per-CIF report + a summary table; renders <stem>_pca.png into out/.
"""
import glob, os, math, re
import numpy as np
import gemmi
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

K_CODATA, K_CHECKCIF = 1.66054, 1.66042
VDW_TABLE = {"H", "O", "N", "C", "Cl", "S"}      # the original core vdW set (the engine now falls back to gemmi)
HERE = os.path.dirname(os.path.abspath(__file__))
CIF_DIR = os.path.join(HERE, "cifs")
OUT = os.path.join(HERE, "out"); os.makedirs(OUT, exist_ok=True)

def fnum(block, tag):
    raw = block.find_value(tag)
    if raw is None: return None
    raw = raw.strip().strip("'\"")
    if raw in ("?", ".", ""): return None
    try: return float(raw.split("(")[0])
    except ValueError: return None

def sval(block, tag):
    raw = block.find_value(tag)
    return raw.strip().strip("'\"") if raw else None

def elem_from(sym, label):
    """Element from _atom_site_type_symbol, stripping oxidation state ('O2-'->O,
    'Fe3+'->Fe); fall back to the label ('C1'->C) when type_symbol is absent/odd."""
    for src in (sym, label):
        m = re.match(r"\s*([A-Za-z]{1,2})", src or "")
        if not m:
            continue
        a = m.group(1)
        cands = [a[0].upper() + a[1].lower(), a[0].upper()] if len(a) == 2 else [a[0].upper()]
        for cand in cands:
            try:
                el = gemmi.Element(cand)
                if el.atomic_number > 0:
                    return el.name
            except Exception:
                pass
    return "X"

def parse_formula(s):
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)\s*(\d*\.?\d*)", s or ""):
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    return out

def wrapped_same(g, s, tol=3e-3):
    d = np.abs((g - s + 0.5) % 1.0 - 0.5)
    return np.all(d < tol)

def validate(path):
    name = os.path.basename(path)
    doc = gemmi.cif.read(path)
    blocks = [b.name for b in doc]
    block = doc[0]
    a, b, c = (fnum(block, f"_cell_length_{x}") for x in "abc")
    al, be, ga = (fnum(block, f"_cell_angle_{x}") for x in ("alpha", "beta", "gamma"))
    if None in (a, b, c, al, be, ga):
        print(f"\n### {name}: could not read cell - skipping"); return None
    cell = gemmi.UnitCell(a, b, c, al, be, ga)
    Vdecl = fnum(block, "_cell_volume")
    Z = fnum(block, "_cell_formula_units_Z")
    MW = fnum(block, "_chemical_formula_weight")
    dens_decl = fnum(block, "_exptl_crystal_density_diffrn")
    F000_decl = fnum(block, "_exptl_crystal_F_000")
    wl = fnum(block, "_diffrn_radiation_wavelength")
    rtype = sval(block, "_diffrn_radiation_type") or sval(block, "_diffrn_radiation_probe")
    sg = sval(block, "_symmetry_space_group_name_H-M") or sval(block, "_space_group_name_H-M_alt")

    # atoms
    # column-wise read: robust to a missing _atom_site_type_symbol (element from label)
    # and a missing _atom_site_occupancy (default 1) -- both happen in SHELX/older CIFs.
    def col(tag): return [str(v) for v in block.find_loop(tag)]
    def num(s):
        try: return float(s.split("(")[0])
        except Exception: return None
    labels = col("_atom_site_label")
    fxc, fyc, fzc = col("_atom_site_fract_x"), col("_atom_site_fract_y"), col("_atom_site_fract_z")
    types, occc = col("_atom_site_type_symbol"), col("_atom_site_occupancy")
    atoms = []
    for i in range(min(len(labels), len(fxc), len(fyc), len(fzc))):
        x, y, z = num(fxc[i]), num(fyc[i]), num(fzc[i])
        if None in (x, y, z):
            continue
        occ = num(occc[i]) if (i < len(occc) and occc[i] not in ("?", ".")) else 1.0
        atoms.append(dict(label=labels[i], sym=elem_from(types[i] if i < len(types) else "", labels[i]),
                          frac=np.array([x, y, z]), occ=occ if occ is not None else 1.0))

    raw_ops = list(block.find_loop("_symmetry_equiv_pos_as_xyz")) or \
              list(block.find_loop("_space_group_symop_operation_xyz"))
    ops = [gemmi.Op(str(o).strip().strip("'\"").replace(" ", "")) for o in raw_ops]
    if not ops:                       # no symop loop -> derive from the H-M / Hall symbol
        nm = sval(block, "_symmetry_space_group_name_H-M") or sval(block, "_space_group_name_H-M_alt")
        hall = sval(block, "_symmetry_space_group_name_Hall") or sval(block, "_space_group_name_Hall")
        sgobj = None
        if nm:
            try: sgobj = gemmi.find_spacegroup_by_name(nm)
            except Exception: sgobj = None
        if sgobj is not None:
            ops = list(sgobj.operations())
        elif hall:
            try: ops = list(gemmi.symops_from_hall(hall))
            except Exception: ops = []
    nops = len(ops)

    def orbit(frac):
        seen = []
        for op in ops:
            g = np.array(op.apply_to_xyz(list(frac))) % 1.0
            if not any(wrapped_same(g, s) for s in seen):
                seen.append(g)
        return len(seen)

    # density triple-check + F000 + special positions + elements
    Vcalc = a*b*c*math.sqrt(max(0.0, 1 - sum(math.cos(math.radians(t))**2 for t in (al,be,ga))
                            + 2*math.prod(math.cos(math.radians(t)) for t in (al,be,ga))))
    # GLOBAL dedup of every (atom x symop) position, then sum mass over uniques.
    # Robust to BOTH special positions AND symmetry-completed atom lists (CSD-issued
    # CIFs list the whole molecule, not the minimal AU -> per-atom orbit*mass overcounts).
    mass = elec = 0.0
    counts, specials, elems = {}, [], set()
    uniq = {}
    for at in atoms:
        el = gemmi.Element(at["sym"]); elems.add(at["sym"])
        m = orbit(at["frac"])
        if nops and m < nops:
            specials.append((at["label"], at["sym"], f"{m}/{nops}"))
        for op in ops:
            g = np.array(op.apply_to_xyz(list(at["frac"]))) % 1.0
            key = (at["sym"], int(round(g[0] * 100)) % 100,
                   int(round(g[1] * 100)) % 100, int(round(g[2] * 100)) % 100)
            uniq.setdefault(key, (el, at["occ"]))
    for el, occ in uniq.values():
        mass += el.weight * occ
        elec += el.atomic_number * occ
        counts[el.name] = counts.get(el.name, 0.0) + occ
    dens_i = K_CODATA * mass / Vdecl if Vdecl else float("nan")
    dens_ii = (K_CHECKCIF * MW * Z / Vdecl) if (MW and Z and Vdecl) else float("nan")
    fsum = parse_formula(sval(block, "_chemical_formula_sum"))
    mw_formula = sum(gemmi.Element(e).weight * n for e, n in fsum.items()) if fsum else None
    mw_bad = bool(mw_formula and MW and abs(mw_formula - MW) / MW > 0.01)
    counts_match = bool(fsum and Z) and all(
        abs(counts.get(e, 0.0) - n * Z) <= max(0.05 * n * Z, 0.15) for e, n in fsum.items())
    atoms_short = bool(fsum and Z) and sum(counts.values()) < sum(n * Z for n in fsum.values()) - 0.5

    def close(x, y, tol=0.03): return bool(x and y and not math.isnan(x) and abs(x - y) / y < tol)
    have_ii = not math.isnan(dens_ii)
    i_ok, ii_ok = close(dens_i, dens_decl), close(dens_ii, dens_decl)
    if dens_decl is None:
        branch = "n/a (no declared density)"
    elif i_ok and (ii_ok or not have_ii):
        branch = "PASS"
    elif have_ii and i_ok and not ii_ok:
        branch = "A-ALERT (ii!=iii: declared density/MW inconsistent)"
    elif not i_ok:
        if mw_bad and counts_match:
            branch = f"FLAG: formula_weight wrong (decl {MW} vs {mw_formula:.2f} from formula_sum)"
        elif atoms_short:
            branch = "BANNER (i<iii: incomplete atoms / H-free / squeeze)"
        else:
            branch = "FLAG: (i)!=(iii) - inspect"
    else:
        branch = "PASS"

    # disorder
    part = [at for at in atoms if at["occ"] < 0.999]
    mut = []
    for i in range(len(part)):
        for j in range(i + 1, len(part)):
            f1, f2 = part[i]["frac"], part[j]["frac"]
            d = np.linalg.norm(np.array([cell.orthogonalize(gemmi.Fractional(*f1)).x - cell.orthogonalize(gemmi.Fractional(*f2)).x,
                                         cell.orthogonalize(gemmi.Fractional(*f1)).y - cell.orthogonalize(gemmi.Fractional(*f2)).y,
                                         cell.orthogonalize(gemmi.Fractional(*f1)).z - cell.orthogonalize(gemmi.Fractional(*f2)).z]))
            if d < 0.9:
                mut.append((part[i]["label"], part[j]["label"], round(float(d), 3)))

    # _geom_hbond oracle (same-cell '.' rows only)
    hb_total = hb_checked = hb_ok = 0
    try:
        hcols = ["_geom_hbond_atom_site_label_D", "_geom_hbond_atom_site_label_A",
                 "_geom_hbond_distance_DA", "_geom_hbond_site_symmetry_A"]
        bylabel = {at["label"]: at for at in atoms}
        for row in block.find(hcols):
            hb_total += 1
            if row[3].strip() in (".", "1_555") and row[0] in bylabel and row[1] in bylabel:
                hb_checked += 1
                dd = np.linalg.norm(np.array([cell.orthogonalize(gemmi.Fractional(*bylabel[row[0]]["frac"])).x,
                                              cell.orthogonalize(gemmi.Fractional(*bylabel[row[0]]["frac"])).y,
                                              cell.orthogonalize(gemmi.Fractional(*bylabel[row[0]]["frac"])).z]) -
                                    np.array([cell.orthogonalize(gemmi.Fractional(*bylabel[row[1]]["frac"])).x,
                                              cell.orthogonalize(gemmi.Fractional(*bylabel[row[1]]["frac"])).y,
                                              cell.orthogonalize(gemmi.Fractional(*bylabel[row[1]]["frac"])).z]))
                if abs(dd - float(row[2].split("(")[0])) < 0.02:
                    hb_ok += 1
    except Exception:
        pass

    # PCA orientation
    heavy = [at for at in atoms if at["sym"] != "H" and at["occ"] > 0.5]
    w2w1 = rms_oop = float("nan")
    if len(heavy) >= 3:
        P = np.array([[cell.orthogonalize(gemmi.Fractional(*at["frac"])).x,
                       cell.orthogonalize(gemmi.Fractional(*at["frac"])).y,
                       cell.orthogonalize(gemmi.Fractional(*at["frac"])).z] for at in heavy])
        Pc = P - P.mean(0)
        evals, evecs = np.linalg.eigh(Pc.T @ Pc)
        if np.linalg.det(evecs) < 0: evecs[:, 0] = -evecs[:, 0]
        rms_oop = float(np.sqrt(((Pc @ evecs[:, 0])**2).mean()))
        w2w1 = float(evals[1] / evals[2]) if evals[2] else float("nan")
        XY = Pc @ np.column_stack([evecs[:, 2], evecs[:, 1]])
        col = {"C": "0.35", "N": "tab:blue", "O": "tab:red", "S": "y", "Cl": "g"}
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.scatter(XY[:, 0], XY[:, 1], s=70, c=[col.get(at["sym"], "k") for at in heavy],
                   edgecolors="k", linewidths=0.3)
        ax.set_aspect("equal"); ax.set_title(f"{name}\nw2/w1={w2w1:.2f} RMSoop={rms_oop:.2f}", fontsize=8)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, os.path.splitext(name)[0] + "_pca.png"), dpi=110)
        plt.close(fig)

    missing_vdw = elems - VDW_TABLE
    print(f"\n### {name}")
    if len(blocks) > 1:
        print(f"  !! MULTI-BLOCK ({len(blocks)} data_ blocks: {blocks}) — engine must FAIL "
              f"unless cif_block is set; harness uses [0]={blocks[0]}")
    print(f"  blocks={len(blocks)}  SG={sg}  nops={nops}  Z={Z}  asym_atoms={len(atoms)}  elems={sorted(elems)}")
    print(f"  V: calc {Vcalc:.2f} vs decl {Vdecl} (d {abs(Vcalc-(Vdecl or 0)):.3f})")
    print(f"  density (i)={dens_i:.4f} (ii)={dens_ii:.4f} (iii)={dens_decl}  -> {branch}")
    print(f"  formula_weight decl={MW} vs from-formula_sum={mw_formula if mw_formula is None else round(mw_formula,2)}  atom-counts-match-formula={counts_match}")
    if F000_decl: print(f"  F(000) calc {elec:.1f} vs decl {F000_decl}")
    print(f"  wavelength={wl} type={rtype}")
    print(f"  special positions: {specials or 'none'}")
    dg = [v for v in block.find_loop("_atom_site_disorder_group")]
    dg_tagged = sum(1 for v in dg if v.strip() not in (".", "?", ""))
    print(f"  disorder partial sites={len(part)}  mutually-exclusive pairs (would bogus-bond)={mut}")
    print(f"  _atom_site_disorder_group: {f'{dg_tagged} atoms tagged (A/B groups present)' if dg_tagged else 'tag ABSENT'}")
    print(f"  _geom_hbond rows={hb_total} same-cell-checked={hb_checked} D..A-match={hb_ok}")
    print(f"  PCA: w2/w1={w2w1:.3f} {'<-- ROLL/SHAPE DEGENERATE' if w2w1>=0.85 else ''}  RMSoop={rms_oop:.3f}")
    if missing_vdw: print(f"  !! elements outside spec vdW table: {sorted(missing_vdw)}")
    return dict(name=name, sg=sg, nops=nops, branch=branch, dens_i=dens_i, dens_iii=dens_decl,
                specials=len(specials), mut=len(mut), wl=wl, rtype=rtype, w2w1=w2w1,
                hb=f"{hb_ok}/{hb_checked}", missing=sorted(missing_vdw))

rows = []
for p in sorted(glob.glob(os.path.join(CIF_DIR, "*.cif"))):
    r = validate(p)
    if r: rows.append(r)

print("\n" + "=" * 96)
print(f"{'file':36} {'SG':12} {'nops':>4} {'dens_i/iii':>11} {'branch':28} {'sp':>2} {'mut':>3} {'w2/w1':>5}")
for r in rows:
    print(f"{r['name'][:36]:36} {str(r['sg'])[:12]:12} {r['nops']:>4} "
          f"{r['dens_i']:.3f}/{str(r['dens_iii']):>5} {r['branch'][:28]:28} "
          f"{r['specials']:>2} {r['mut']:>3} {r['w2w1']:>5.2f}")
