"""
crystal_engine.py  -  CIF -> validated crystal data. The single source of truth for the
crystal family: parse, symmetry-expand, density triple-check, geometry, H-bonds. The
figures (crystal_pxrd, crystal_view) consume what this returns; they never recompute, so a
figure can never assert a contact the table doesn't list.

Pure computation + Tier-1 gates (no plotting). gemmi is imported lazily (heavy dep; a
non-CIF job never pulls it). Conventions/constants are pinned in the block below;
the logic was validated against the CIF edge-case set in skill_validation/crystal/ (ELAINM, urea, aspirin, HMT,
ferrocene, flufenamic, lactose, L-alanine, PTU-ellagic + synthetic broken-input CIFs).

    load(cfg)                 parse + block-select + symmetry + atoms (robust to CSD/SHELX quirks)
    expand(struct)            unique (frac,sym,occ) cell positions (global atom x symop dedup)
    densities(struct, cfg)    the signature triple-check + F(000) + element counts + branch gate
    special_positions(struct) atoms on a symmetry element (orbit < n_ops)
    bonds(struct, cfg)        covalent-radii bonds, disorder-aware
    hbonds(struct, cfg)       X-H-normalized D-H...A with symmetry-expanded search + geom_hbond xref
    validate(struct, cfg)     run all Tier-1 gates, assemble the validation-table rows
    table_text / write_csv    emit the validation table (the CORE deliverable)
"""
from __future__ import annotations
import math
import re
import numpy as np
from . import verify

# ---- pinned constants (sources inline) ------------------------------------------
K_CODATA = 1.66054        # CODATA atomic-mass-constant factor -> assumption-free density (i)
K_CHECKCIF = 1.66042      # checkCIF's legacy DENSD01 factor -> formula density (ii)
NEUTRON_XH = {"C": 1.083, "N": 1.009, "O": 0.983}   # Allen & Bruno 2010 (X-H normalization)
# Bondi (1964) vdW radii (A); core set pinned, others fall back to gemmi.Element.vdw_r
VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
       "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "P": 1.80}
JEFFREY = ((2.2, 2.5, "strong"), (2.5, 3.2, "moderate"), (3.2, 4.0, "weak"))  # D...A ranges (A)
BOND_TOL = 0.40           # A added to covalent-radii sum for bond perception
DISORDER_MIN = 0.90       # A; sub-unity sites closer than this are mutually exclusive (no bond)
GEOM_OUTLIER = 0.25       # A; bond length off the covalent-radii sum by more than this -> flagged
GEOM_CONTACT_MIN = 90.0   # deg; contacts in [this, hbond_angle_min) are "geometric contacts" (listed, not asserted)
STD_WL = {"Cu-Ka": 1.5418, "Cu-Ka1": 1.5406, "Mo": 0.71073,   # standard Kα wavelengths (A)
          "Ag": 0.56086, "Ga": 1.3414, "In": 0.51359}


def _gemmi():
    try:
        import gemmi
        return gemmi
    except ImportError as e:
        raise ImportError(
            "Required package not found: gemmi. Install with:\n"
            "    python -m pip install gemmi\n"
            "(the crystal family needs gemmi for CIF parsing + symmetry).") from e


class Structure:
    """Parsed + cached crystal data. Build with crystal_engine.load(cfg)."""
    def __init__(self, **kw):
        self.__dict__.update(kw)
        self._uniq = None
        self._alt = None

    def cart(self, frac):
        p = self.cell.orthogonalize(_gemmi().Fractional(float(frac[0]), float(frac[1]), float(frac[2])))
        return np.array([p.x, p.y, p.z])


# ----------------------------------------------------------------- parse helpers
def _num(s):
    try:
        return float(str(s).split("(")[0])
    except (ValueError, AttributeError, TypeError):
        return None


def _elem_from(type_symbol, label):
    """Element from _atom_site_type_symbol, stripping an oxidation-state suffix
    ('O2-'->O, 'Fe3+'->Fe); fall back to the label ('C1'->C) when type_symbol is
    absent/odd. '' if unrecognisable. (Both quirks seen in the validation set.)"""
    g = _gemmi()
    for src in (type_symbol, label):
        m = re.match(r"\s*([A-Za-z]{1,2})", src or "")
        if not m:
            continue
        a = m.group(1)
        cands = [a[0].upper() + a[1].lower(), a[0].upper()] if len(a) == 2 else [a[0].upper()]
        for c in cands:
            try:
                el = g.Element(c)
                if el.atomic_number > 0:
                    return el.name
            except Exception:
                pass
    return ""


def _parse_formula(s):
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)\s*(\d*\.?\d*)", s or ""):
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    return out


def _vdw(sym):
    if sym in VDW:
        return VDW[sym]
    try:
        r = _gemmi().Element(sym).vdw_r
        return r if r and r > 0 else 1.70
    except Exception:
        return 1.70


def _covsum(s1, s2):
    g = _gemmi()
    return g.Element(s1).covalent_r + g.Element(s2).covalent_r


def _angle(p, q, r):
    v1, v2 = p - q, r - q
    n = np.linalg.norm(v1) * np.linalg.norm(v2)
    if n == 0:
        return float("nan")
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2) / n)))))


# --------------------------------------------------------------------- load
def load(cfg):
    """Parse cfg.cif_path, select the data block, read cell / symmetry / atoms with the
    robustness the validation set demanded: oxidation-state type symbols, missing
    _atom_site_type_symbol or _atom_site_occupancy columns, symop loop OR a bare
    space-group symbol. Gates: file parses, ONE block selected (multi-block FAILs unless
    cfg.cif_block is set), cell + symmetry + atoms present."""
    g = _gemmi()
    if not cfg.cif_path:
        raise ValueError("cfg.cif_path is not set")
    doc = g.cif.read(cfg.cif_path)
    blocks = [b.name for b in doc]
    f = []

    # ---- block selection (spec §4 step 1)
    block = doc[0]
    if cfg.cif_block is not None:
        sel = str(cfg.cif_block)
        try:
            block = doc[int(sel)] if sel.lstrip("-").isdigit() else doc[sel]
        except Exception:
            f.append(("FAIL", f"cif_block '{cfg.cif_block}' not found; blocks={blocks}"))
    elif len(blocks) > 1:
        f.append(("FAIL", f"multi-block CIF ({len(blocks)} blocks: {blocks}); set cfg.cif_block"))

    def sval(tag):
        v = block.find_value(tag)
        return v.strip().strip("'\"") if v else None

    def fval(tag):
        return _num(sval(tag))

    a, b, c = fval("_cell_length_a"), fval("_cell_length_b"), fval("_cell_length_c")
    al, be, ga = fval("_cell_angle_alpha"), fval("_cell_angle_beta"), fval("_cell_angle_gamma")
    if None in (a, b, c, al, be, ga):
        f.append(("FAIL", "cell parameters missing/unparseable"))
        a, b, c = a or 1.0, b or 1.0, c or 1.0
        al, be, ga = al or 90.0, be or 90.0, ga or 90.0
    cell = g.UnitCell(a, b, c, al, be, ga)

    # ---- symmetry operators: explicit loop, else derive from H-M / Hall symbol
    raw = list(block.find_loop("_symmetry_equiv_pos_as_xyz")) or \
        list(block.find_loop("_space_group_symop_operation_xyz"))
    ops = [g.Op(str(o).strip().strip("'\"").replace(" ", "")) for o in raw]
    sg_hm = sval("_symmetry_space_group_name_H-M") or sval("_space_group_name_H-M_alt")
    if not ops:
        hall = sval("_symmetry_space_group_name_Hall") or sval("_space_group_name_Hall")
        sgo = None
        if sg_hm:
            try:
                sgo = g.find_spacegroup_by_name(sg_hm)
            except Exception:
                sgo = None
        if sgo is not None:
            ops = list(sgo.operations())
        elif hall:
            try:
                ops = list(g.symops_from_hall(hall))
            except Exception:
                ops = []
    if not ops:
        ops = [g.Op("x,y,z")]
        f.append(("WARN", "no symmetry operators found; assuming P1"))

    # ---- atoms: column-wise (type_symbol/occupancy optional)
    def col(tag):
        return [str(v) for v in block.find_loop(tag)]

    labels, occs = col("_atom_site_label"), col("_atom_site_occupancy")
    fx, fy, fz = col("_atom_site_fract_x"), col("_atom_site_fract_y"), col("_atom_site_fract_z")
    types, uiso = col("_atom_site_type_symbol"), col("_atom_site_U_iso_or_equiv")
    dis_asm, dis_grp = col("_atom_site_disorder_assembly"), col("_atom_site_disorder_group")
    atoms = []
    for i in range(min(len(labels), len(fx), len(fy), len(fz))):
        x, y, z = _num(fx[i]), _num(fy[i]), _num(fz[i])
        if None in (x, y, z):
            continue
        sym = _elem_from(types[i] if i < len(types) else "", labels[i])
        if not sym:
            continue
        occ = _num(occs[i]) if (i < len(occs) and occs[i] not in ("?", ".")) else 1.0
        atoms.append({"label": labels[i], "sym": sym,
                      "frac": np.array([x, y, z], float), "occ": occ if occ is not None else 1.0,
                      "u_iso": _num(uiso[i]) if i < len(uiso) else None,
                      "dis_asm": dis_asm[i] if i < len(dis_asm) else None,
                      "dis_grp": dis_grp[i] if i < len(dis_grp) else None})
    if not atoms:
        f.append(("FAIL", "no atom sites parsed"))

    # anisotropic displacement parameters (Phase 2): {label: 3x3 U tensor in the CIF basis}
    al_ = col("_atom_site_aniso_label")
    au = {k: col("_atom_site_aniso_U_" + k) for k in ("11", "22", "33", "23", "13", "12")}
    aniso = {}
    for i, lab in enumerate(al_):
        try:
            v = {k: _num(au[k][i]) for k in au}
            if any(x is None for x in v.values()):
                continue
            aniso[lab] = np.array([[v["11"], v["12"], v["13"]],
                                   [v["12"], v["22"], v["23"]],
                                   [v["13"], v["23"], v["33"]]], float)
        except Exception:
            continue

    verify._resolve(f, cfg, "load")     # raises under cfg.strict on any FAIL

    declared = {"V": fval("_cell_volume"), "Z": fval("_cell_formula_units_Z"),
                "MW": fval("_chemical_formula_weight"),
                "density": fval("_exptl_crystal_density_diffrn"),
                "F000": fval("_exptl_crystal_F_000"),
                "wavelength": fval("_diffrn_radiation_wavelength"),
                "rad_type": sval("_diffrn_radiation_type") or sval("_diffrn_radiation_probe"),
                "temperature": sval("_diffrn_ambient_temperature") or sval("_cell_measurement_temperature"),
                "formula_sum": sval("_chemical_formula_sum"),
                # echo-only refinement metadata (reported, never recomputed — keeps the CIF-only line honest)
                "R_gt": fval("_refine_ls_R_factor_gt"), "R_all": fval("_refine_ls_R_factor_all"),
                "gof": fval("_refine_ls_goodness_of_fit_ref"),
                "reflns_total": fval("_reflns_number_total"), "reflns_gt": fval("_reflns_number_gt"),
                "theta_max": fval("_diffrn_reflns_theta_max"),
                "size_max": fval("_exptl_crystal_size_max"), "size_mid": fval("_exptl_crystal_size_mid"),
                "size_min": fval("_exptl_crystal_size_min")}
    return Structure(name=block.name, blocks=blocks, block=block, cell=cell, ops=ops,
                     atoms=atoms, aniso=aniso, sg_hm=sg_hm, declared=declared,
                     a=a, b=b, c=c, al=al, be=be, ga=ga)


# --------------------------------------------------------------------- expansion
def expand(struct):
    """Every (atom x symop) image in the unit cell, deduplicated GLOBALLY by wrapped
    position+element. Robust to special positions AND symmetry-completed atom lists
    (CSD-issued CIFs list the whole molecule, not the minimal AU). Cached."""
    if struct._uniq is not None:
        return struct._uniq
    uniq = {}
    for at in struct.atoms:
        for op in struct.ops:
            gpos = np.array(op.apply_to_xyz(list(at["frac"]))) % 1.0
            key = (at["sym"], int(round(gpos[0] * 100)) % 100,
                   int(round(gpos[1] * 100)) % 100, int(round(gpos[2] * 100)) % 100)
            uniq.setdefault(key, {"sym": at["sym"], "frac": gpos,
                                  "occ": at["occ"], "label": at["label"]})
    struct._uniq = list(uniq.values())
    return struct._uniq


def special_positions(struct):
    """Atoms whose symmetry orbit is smaller than the number of operators (i.e. they sit
    on a symmetry element). Returns [(label, sym, 'orbit/nops')]."""
    nops = len(struct.ops)
    out = []
    for at in struct.atoms:
        seen = []
        for op in struct.ops:
            gp = np.array(op.apply_to_xyz(list(at["frac"]))) % 1.0
            if not any(np.all(np.abs((gp - s + 0.5) % 1.0 - 0.5) < 3e-3) for s in seen):
                seen.append(gp)
        if len(seen) < nops:
            out.append((at["label"], at["sym"], f"{len(seen)}/{nops}"))
    return out


# --------------------------------------------------------------- ADP / ellipsoids (Phase 2)
def has_adp(struct):
    return bool(getattr(struct, "aniso", None))


def _orth(struct):
    """Orthogonalization matrix M (frac -> cart); columns are the cell vectors a, b, c."""
    return np.column_stack([struct.cart([1, 0, 0]), struct.cart([0, 1, 0]), struct.cart([0, 0, 1])])


def u_cart(struct, U):
    """CIF anisotropic U (the a*_i a*_j basis) -> Cartesian mean-square-displacement tensor
    (Å²): U_cart = (M N) U (M N)^T, with M the orthogonalization matrix and N = diag(a*,b*,c*).
    Validated by U_eq = trace(U_cart)/3 matching the CIF's _atom_site_U_iso_or_equiv."""
    M = _orth(struct)
    rc = struct.cell.reciprocal()
    A = M @ np.diag([rc.a, rc.b, rc.c])
    return A @ np.asarray(U, float) @ A.T


def u_eq(struct, U):
    return float(np.trace(u_cart(struct, U)) / 3.0)


def _adp_scale(p):
    """Ellipsoid radius in σ units enclosing probability p of a 3-D Gaussian (chi-3 quantile).
    50% -> 1.5382 (the ORTEP default)."""
    try:
        from scipy.stats import chi2
        return float(chi2.ppf(p, 3)) ** 0.5
    except Exception:
        return {0.50: 1.5382, 0.30: 1.0972, 0.95: 2.7955, 0.99: 3.3682}.get(round(p, 2), 1.5382)


def adp_axes(struct, U, prob=0.5):
    """(semi_axes[3], V[3x3] eigenvectors as columns) of the displacement ellipsoid at `prob`.
    Returns (None, None) when U is not positive-definite (non-physical ADP -> caller degrades)."""
    w, V = np.linalg.eigh(u_cart(struct, U))
    if np.any(w <= 1e-9):
        return None, None
    return _adp_scale(prob) * np.sqrt(w), V


def op_rot_cart(struct, op):
    """Cartesian rotation of a symmetry operator — to rotate a symmetry image's ADP:
    U_image = R U_cart R^T."""
    M = _orth(struct)
    Rf = np.array(op.rot, float) / op.DEN
    return M @ Rf @ np.linalg.inv(M)


# --------------------------------------------------------------------- density
def densities(struct, cfg):
    """The signature density triple-check + F(000) + element-count cross-check, with the
    full branch tree (spec §5). The expansion-implicated case HARD-FAILS regardless of
    cfg.strict (drawing wrong atoms is a correctness error, not a preference)."""
    g = _gemmi()
    d = struct.declared
    V = d["V"] or struct.cell.volume
    uniq = expand(struct)
    mass = sum(g.Element(u["sym"]).weight * u["occ"] for u in uniq)
    elec = sum(g.Element(u["sym"]).atomic_number * u["occ"] for u in uniq)
    counts = {}
    for u in uniq:
        counts[u["sym"]] = counts.get(u["sym"], 0.0) + u["occ"]

    di = K_CODATA * mass / V if V else float("nan")
    dii = (K_CHECKCIF * d["MW"] * d["Z"] / V) if (d["MW"] and d["Z"] and V) else float("nan")
    diii = d["density"]
    fsum = _parse_formula(d["formula_sum"])
    mw_formula = sum(g.Element(e).weight * n for e, n in fsum.items()) if fsum else None

    def close(x, y, tol=0.03):
        return bool(x and y and not (isinstance(x, float) and math.isnan(x)) and abs(x - y) / y < tol)

    have_ii = not math.isnan(dii)
    i_ok, ii_ok = close(di, diii), close(dii, diii)
    mw_bad = bool(mw_formula and d["MW"] and abs(mw_formula - d["MW"]) / d["MW"] > 0.01)
    counts_match = bool(fsum and d["Z"]) and all(
        abs(counts.get(e, 0.0) - n * d["Z"]) <= max(0.05 * n * d["Z"], 0.15) for e, n in fsum.items())
    atoms_short = bool(fsum and d["Z"]) and \
        sum(counts.values()) < sum(n * d["Z"] for n in fsum.values()) - 0.5

    hard_fail = False
    if diii is None:
        branch, sev = "no declared density (cannot gate)", "INFO"
    elif i_ok and (ii_ok or not have_ii):
        branch, sev = "PASS (all consistent)", "INFO"
    elif have_ii and i_ok and not ii_ok:
        branch, sev = "A-ALERT: declared density/formula inconsistent (ii != iii)", "WARN"
    elif not i_ok and mw_bad and counts_match:
        branch, sev = (f"declared formula_weight wrong ({d['MW']} vs {mw_formula:.2f} "
                       f"from formula_sum); coordinates fine", "WARN")
    elif not i_ok and atoms_short:
        branch, sev = "incomplete atom list (i<iii: H-free / SQUEEZE'd solvent / disorder)", "WARN"
    elif not i_ok and counts_match and have_ii and ii_ok:
        branch, sev, hard_fail = "EXPANSION WRONG (atoms complete yet i != iii) — figure refused", "FAIL", True
    elif not i_ok:
        branch, sev = "(i) != (iii) — inspect (declared metadata or cell)", "WARN"
    else:
        branch, sev = "PASS", "INFO"

    rd = (diii / dii) if (diii and have_ii and dii) else None
    rd_level = None
    if rd is not None:
        rd_level = ("A" if not (0.90 <= rd <= 1.10) else
                    "B" if not (0.95 <= rd <= 1.05) else
                    "C" if not (0.99 <= rd <= 1.01) else None)
    out = {"i": di, "ii": dii, "iii": diii, "branch": branch, "rd": rd, "rd_level": rd_level,
           "mass": mass, "F000_calc": elec, "F000_decl": d["F000"],
           "counts": counts, "mw_formula": mw_formula, "mw_declared": d["MW"]}

    f = [(sev, f"density (i)={di:.4f} (ii)={dii:.4f} (iii)={diii} -> {branch}")]
    if rd_level:
        f.append(("WARN" if rd_level in ("A", "B") else "INFO",
                  f"density RD={rd:.4f} -> checkCIF {rd_level}-alert (DENSD01: declared/formula-DEN)"))
    if d["F000"] and abs(elec - d["F000"]) / d["F000"] > 0.05:
        f.append(("WARN", f"F(000) calc {elec:.0f} vs declared {d['F000']} (>5% off — corroborates the density flag)"))
    if hard_fail:
        for s, m in f:
            print(f"  [{s}] density: {m}")
        raise verify.GateError(f"density: {branch}")     # always raises (correctness)
    verify._resolve(f, cfg, "density")
    return out


# --------------------------------------------------------------------- disorder
def disorder_alternatives(struct):
    """label -> set of mutually-exclusive disorder-alternative labels. Two sources, unioned:
    (1) NEAR-COINCIDENCE — sub-unity (occ<1) sites of the SAME element within DISORDER_MIN are
        alternatives (handles UNtagged disorder, e.g. ELAINM's water over 3 partial O sites);
    (2) DISORDER TAGS — same _atom_site_disorder_assembly, different _atom_site_disorder_group
        (handles A/B components farther apart than DISORDER_MIN, e.g. flufenamic's CF3).
    Transitively closed and cached, so a donor never H-bonds to its own alternative site."""
    if struct._alt is not None:
        return struct._alt
    alt = {a["label"]: set() for a in struct.atoms}

    def link(la, lb):
        if la != lb:
            alt[la].add(lb); alt[lb].add(la)

    sub = [a for a in struct.atoms if a["occ"] < 0.999]
    for i in range(len(sub)):
        for j in range(i + 1, len(sub)):
            if sub[i]["sym"] != sub[j]["sym"]:
                continue
            if float(np.linalg.norm(struct.cart(sub[i]["frac"]) - struct.cart(sub[j]["frac"]))) < DISORDER_MIN:
                link(sub[i]["label"], sub[j]["label"])

    tagged = [a for a in struct.atoms if a.get("dis_grp") not in (None, ".", "", "?")]
    for i in range(len(tagged)):
        for j in range(i + 1, len(tagged)):
            if tagged[i].get("dis_asm") == tagged[j].get("dis_asm") \
                    and tagged[i]["dis_grp"] != tagged[j]["dis_grp"]:
                link(tagged[i]["label"], tagged[j]["label"])

    changed = True                                        # transitive closure (chain -> full cluster)
    while changed:
        changed = False
        for k in alt:
            new = set().union(*(alt[m] for m in alt[k])) if alt[k] else set()
            new.discard(k)
            if not new <= alt[k]:
                alt[k] |= new; changed = True
    struct._alt = alt
    return alt


# --------------------------------------------------------------------- bonds
def bonds(struct, cfg):
    """Covalent bonds within the cell, DISORDER-AWARE: two sub-unity sites closer than
    DISORDER_MIN are alternative positions, not bonded (else naive perception invents
    0.5 A 'bonds'). Returns list of (i, j, dist) index pairs into expand(struct)."""
    atoms = expand(struct)
    alt = disorder_alternatives(struct)
    xyz = [struct.cart(u["frac"]) for u in atoms]
    out = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            dd = float(np.linalg.norm(xyz[i] - xyz[j]))
            if dd <= 0.4 or dd > _covsum(atoms[i]["sym"], atoms[j]["sym"]) + BOND_TOL:
                continue
            if atoms[i]["occ"] < 1 and atoms[j]["occ"] < 1 and dd < DISORDER_MIN:
                continue                                  # mutually-exclusive disorder (near-coincident)
            if atoms[j]["label"] in alt.get(atoms[i]["label"], set()):
                continue                                  # disorder alternative (incl. group-tagged)
            out.append((i, j, round(dd, 3)))
    return out


def _sym_code(struct, label, frac):
    """checkCIF-style symmetry code 'n_pqr' for an atom image at `frac` generated from the
    asym-unit atom `label`: n = 1-based symop index, pqr = 5 + lattice translation.
    '.' = identity in the reference cell; '?' if not recoverable."""
    asym = next((a for a in struct.atoms if a["label"] == label), None)
    if asym is None:
        return "?"
    fr = np.asarray(frac, float)
    for n, op in enumerate(struct.ops, 1):
        base = np.array(op.apply_to_xyz(list(asym["frac"])), float)
        t = np.round(fr - base)
        if np.allclose(base + t, fr, atol=2e-2):
            if n == 1 and np.allclose(t, 0):
                return "."
            return "%d_%d%d%d" % (n, 5 + int(t[0]), 5 + int(t[1]), 5 + int(t[2]))
    return "?"


# --------------------------------------------------------------------- geometry (Block B)
def geometry(struct, cfg):
    """Block B: selected bond lengths flagged against the covalent-radii sum (Cordero/gemmi
    radii ±0.25 A). This is a COARSE sanity bound, NOT a Mogul/CSD percentile (Mogul is
    licensed/unavailable — don't overclaim distributional rigour). Disorder-aware; bonds are
    deduped to unique (label-pair, length) types over the cell + nearest neighbours."""
    cell_atoms = [{**u, "xyz": struct.cart(u["frac"])} for u in expand(struct)]
    alt = disorder_alternatives(struct)
    origin = struct.cart(np.zeros(3))
    sup = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                shift = struct.cart(np.array([di, dj, dk], float)) - origin
                for a in cell_atoms:
                    sup.append({"sym": a["sym"], "label": a["label"], "occ": a["occ"],
                                "xyz": a["xyz"] + shift})
    rows, seen, n_out = [], set(), 0
    for a in cell_atoms:
        for b in sup:
            dd = float(np.linalg.norm(a["xyz"] - b["xyz"]))
            if dd <= 0.4:
                continue
            cs = _covsum(a["sym"], b["sym"])
            if dd > cs + BOND_TOL:
                continue
            if a["occ"] < 1 and b["occ"] < 1 and dd < DISORDER_MIN:
                continue                                          # mutually-exclusive disorder (near-coincident)
            if b["label"] in alt.get(a["label"], set()):
                continue                                          # disorder alternative (incl. group-tagged)
            key = tuple(sorted([a["label"], b["label"]]) + [round(dd, 2)])
            if key in seen:
                continue
            seen.add(key)
            outlier = abs(dd - cs) > GEOM_OUTLIER
            n_out += int(outlier)
            rows.append({"a": a["label"], "b": b["label"], "len": round(dd, 3),
                         "covsum": round(cs, 3), "outlier": outlier})
    f = [("INFO", f"Block B: {len(rows)} unique bonds, {n_out} outside the covalent envelope "
                  f"(±{GEOM_OUTLIER} A; coarse bound, not Mogul)")]
    for bd in rows:
        if bd["outlier"]:
            f.append(("WARN", f"bond {bd['a']}-{bd['b']} {bd['len']} A vs covalent sum "
                              f"{bd['covsum']} A — outlier"))
    verify._resolve(f, cfg, "geometry")
    return {"bonds": rows, "n_outliers": n_out}


def _wavelength_check(struct):
    """Flag the declared wavelength unless it is within 5e-4 A of a standard Kα, or the source
    is neutron/synchrotron/electron (where any λ is legitimate)."""
    wl = struct.declared["wavelength"]
    rt = (struct.declared["rad_type"] or "").lower()
    if wl is None:
        return ("INFO", "wavelength not declared")
    if any(k in rt for k in ("neutron", "synchrotron", "electron")):
        return ("INFO", f"wavelength {wl} A ({rt}) — non-Kα source, accepted")
    near = [k for k, v in STD_WL.items() if abs(wl - v) < 5e-4]
    if near:
        return ("INFO", f"wavelength {wl} A matches {near[0]}")
    return ("WARN", f"wavelength {wl} A not within 5e-4 of a standard Kα "
                    f"(type='{rt or '?'}') — verify")


def _sg_check(struct):
    """Cross-check the symmetry operators against the declared space-group name: the op count
    must equal the named group's order (catches an incomplete symop list / wrong setting), and
    where derivable, the group recovered FROM the ops must match the declared H-M symbol."""
    g = _gemmi()
    decl = struct.sg_hm
    named = None
    if decl:
        try:
            named = g.find_spacegroup_by_name(decl)
        except Exception:
            named = None
    if named is not None:
        order = len(list(named.operations()))
        if order != len(struct.ops):
            return ("WARN", f"space group '{decl}' has order {order} but {len(struct.ops)} "
                            f"symops are present — incomplete/extra symmetry?")
    try:
        from_ops = g.find_spacegroup_by_ops(g.GroupOps(list(struct.ops)))
        if from_ops is not None and decl and \
                from_ops.hm.replace(" ", "") != decl.replace(" ", ""):
            return ("INFO", f"space group from ops = '{from_ops.hm}' vs declared '{decl}' "
                            f"(setting/naming difference — verify)")
    except Exception:
        pass
    return ("INFO", f"space group {decl}: {len(struct.ops)} symops, consistent")


# --------------------------------------------------------------------- H-bonds
def hbonds(struct, cfg):
    """X-H-normalized D-H...A detection with a symmetry-expanded (3x3x3) acceptor search,
    classified by D...A (Jeffrey). Donor/acceptor sets are separate (cfg). Sub-floor
    contacts are demoted, not reported as H-bonds. Cross-checks the CIF's _geom_hbond
    loop on D...A (normalization-invariant) where present."""
    g = _gemmi()
    donors = set(cfg.hbond_donors) | ({"C"} if cfg.hbond_weak else set())
    acceptors = set(cfg.hbond_acceptors)
    floor = cfg.hbond_angle_min
    alt = disorder_alternatives(struct)
    dis_pairs = set()
    f = []
    if cfg.hbond_weak and floor >= 120.0:
        f.append(("INFO", "weak donors enabled but angle floor still 120 deg — most weak "
                          "H-bonds are bent; consider lowering hbond_angle_min toward 90"))

    cell_atoms = [{**u, "xyz": struct.cart(u["frac"])} for u in expand(struct)]
    Hs = [a for a in cell_atoms if a["sym"] == "H"]

    # 3x3x3 supercell of HEAVY atoms: the donor search must see molecules that straddle
    # the cell boundary (not just the [0,1) image, or O-H donors on a boundary-spanning
    # molecule are missed) + the acceptor subset for the contact search.
    sup_heavy, sup_acc = [], []
    origin = struct.cart(np.zeros(3))
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                shift = struct.cart(np.array([di, dj, dk], float)) - origin
                for a in cell_atoms:
                    if a["sym"] == "H":
                        continue
                    rec = {"sym": a["sym"], "label": a["label"], "cell": (di, dj, dk),
                           "frac": a["frac"] + np.array([di, dj, dk], float), "xyz": a["xyz"] + shift}
                    sup_heavy.append(rec)
                    if a["sym"] in acceptors:
                        sup_acc.append(rec)

    results, seen, contacts, seen_c = [], set(), [], set()
    for h in Hs:
        if not sup_heavy:
            break
        D = min(sup_heavy, key=lambda a: np.linalg.norm(a["xyz"] - h["xyz"]))
        dDH = float(np.linalg.norm(D["xyz"] - h["xyz"]))
        if dDH > _covsum(D["sym"], "H") + BOND_TOL or D["sym"] not in donors:
            continue
        hpos = h["xyz"]
        if cfg.xh_normalize and D["sym"] in NEUTRON_XH and dDH > 0:
            hpos = D["xyz"] + (h["xyz"] - D["xyz"]) / dDH * NEUTRON_XH[D["sym"]]
        altD = alt.get(D["label"], set())
        for A in sup_acc:
            DA = float(np.linalg.norm(A["xyz"] - D["xyz"]))
            if DA < 0.4 or DA > 4.0:
                continue
            if A["label"] in altD:                        # acceptor is a disorder-alternative of the donor
                dis_pairs.add((D["label"], A["label"]))
                continue
            HA = float(np.linalg.norm(A["xyz"] - hpos))
            if HA > _vdw("H") + _vdw(A["sym"]):
                continue
            ang = _angle(D["xyz"], hpos, A["xyz"])
            rec = {"D": D["label"], "Dsym": D["sym"], "H": h["label"],
                   "A": A["label"], "Asym": A["sym"], "DH": round(dDH, 3),
                   "HA": round(HA, 3), "DA": round(DA, 3), "angle": round(ang, 1),
                   "class": next((k for lo, hi, k in JEFFREY if lo <= DA < hi), "long"),
                   "sym": _sym_code(struct, A["label"], A["frac"]), "cell": A["cell"]}
            if ang < floor:                               # below the floor: geometric contact, not an asserted H-bond
                if ang >= GEOM_CONTACT_MIN:
                    ck = (D["label"], h["label"], A["label"], round(DA, 2))
                    if ck not in seen_c:
                        seen_c.add(ck); rec["class"] = "contact"; contacts.append(rec)
                continue
            key = (D["label"], h["label"], A["label"], round(DA, 2))
            if key in seen:
                continue
            seen.add(key)
            results.append(rec)

    # cross-check the CIF's own _geom_hbond loop (D...A is normalization-invariant)
    xref = {"rows": 0, "matched": 0}
    try:
        D_l = [str(v) for v in struct.block.find_loop("_geom_hbond_atom_site_label_D")]
        A_l = [str(v) for v in struct.block.find_loop("_geom_hbond_atom_site_label_A")]
        DA_l = [_num(v) for v in struct.block.find_loop("_geom_hbond_distance_DA")]
        found_DA = [r["DA"] for r in results]
        for i in range(min(len(D_l), len(A_l), len(DA_l))):
            xref["rows"] += 1
            if DA_l[i] is not None and any(abs(DA_l[i] - x) < 0.03 for x in found_DA):
                xref["matched"] += 1
        if xref["rows"]:
            f.append(("INFO", f"_geom_hbond cross-check: matched {xref['matched']}/{xref['rows']} "
                              f"deposited bonds on D...A (others may be below the {floor:.0f} deg floor)"))
    except Exception:
        pass

    if dis_pairs:
        f.append(("INFO", f"suppressed {len(dis_pairs)} donor/own-disorder-alternative acceptor "
                          f"pair(s): {sorted(dis_pairs)}"))
    f.insert(0, ("INFO", f"{len(results)} H-bonds + {len(contacts)} geometric contacts "
                         f"(donors={sorted(donors)}, acceptors={sorted(acceptors)}, floor={floor:.0f} deg)"))
    verify._resolve(f, cfg, "hbonds")
    return {"bonds": results, "contacts": contacts, "xref": xref}


# --------------------------------------------------------------------- validate + table
def validate(struct, cfg):
    """Run the Tier-1 gates and assemble the validation table. Returns a report dict
    {density, special, hbonds, volume, wavelength, ...} — the CORE deliverable that the
    figures consume."""
    g = _gemmi()
    # cell volume recompute (assumption-free) vs declared
    ca, cb, cg = (math.cos(math.radians(t)) for t in (struct.al, struct.be, struct.ga))
    Vcalc = struct.a * struct.b * struct.c * math.sqrt(
        max(0.0, 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg))
    Vdecl = struct.declared["V"]
    fV = []
    if Vdecl and abs(Vcalc - Vdecl) / Vdecl > 0.01:
        fV.append(("WARN", f"cell volume calc {Vcalc:.2f} vs declared {Vdecl} (>1% off)"))
    verify._resolve(fV, cfg, "volume")

    dens = densities(struct, cfg)
    spec = special_positions(struct)
    geom = geometry(struct, cfg)
    hb = hbonds(struct, cfg)
    wlchk = _wavelength_check(struct); verify._resolve([wlchk], cfg, "wavelength")
    sgchk = _sg_check(struct); verify._resolve([sgchk], cfg, "spacegroup")
    elems = sorted({a["sym"] for a in struct.atoms})
    missing = [e for e in elems if e not in VDW]
    fE = []
    if missing:
        fE.append(("INFO", f"elements outside the pinned vdW set {sorted(VDW)}: {missing} "
                           f"(using gemmi vdW radii as fallback)"))
    verify._resolve(fE, cfg, "elements")

    return {"name": struct.name, "blocks": struct.blocks, "sg_hm": struct.sg_hm,
            "nops": len(struct.ops), "elems": elems, "missing_vdw": missing,
            "V_calc": Vcalc, "V_decl": Vdecl, "density": dens, "special": spec,
            "geometry": geom, "wl_check": wlchk, "sg_check": sgchk, "floor": cfg.hbond_angle_min,
            "hbonds": hb["bonds"], "contacts": hb["contacts"], "hbond_xref": hb["xref"],
            "declared": struct.declared}


def table_text(report):
    """Human-readable validation table (the text deliverable)."""
    d, dn = report["density"], report["declared"]
    g = report["geometry"]
    rd = (f"  RD            {d['rd']:.4f}" +
          (f"  -> checkCIF {d['rd_level']}-alert" if d.get("rd_level") else "  (ok)")) if d.get("rd") else None
    size = "x".join(f"{dn[k]}" for k in ("size_max", "size_mid", "size_min")) if dn.get("size_max") else "?"
    L = [f"VALIDATION TABLE — {report['name']}  ({report['sg_hm']}, {report['nops']} ops)",
         "-" * 64,
         "Block A — crystal data",
         f"  cell volume   calc {report['V_calc']:.2f}   decl {report['V_decl']}",
         f"  density (i)   {d['i']:.4f}  (expanded-cell, assumption-free, CODATA)",
         f"  density (ii)  {d['ii']:.4f}  (checkCIF formula, 1.66042)",
         f"  density (iii) {d['iii']}  (declared)"]
    if rd:
        L.append(rd)
    L += [f"     -> {d['branch']}",
          f"  F(000)        calc {d['F000_calc']:.0f}   decl {d['F000_decl']}",
          f"  formula wt    from formula_sum {d['mw_formula']}   decl {d['mw_declared']}",
          f"  elements      {report['elems']}" + (f"  (outside vdW table: {report['missing_vdw']})"
                                                  if report["missing_vdw"] else ""),
          f"  wavelength    {dn['wavelength']}  ({dn['rad_type']})  [{report['wl_check'][1]}]",
          f"  space group   {report['sg_check'][1]}",
          f"  special pos   {report['special'] or 'none'}",
          f"  echo (not recomputed)  T={dn['temperature']}  R(gt)={dn.get('R_gt')}  GoF={dn.get('gof')}  "
          f"reflns={dn.get('reflns_total')}  th_max={dn.get('theta_max')}  size={size}",
          "Block B — geometry (bond lengths vs covalent-radii sum +/-0.25 A; coarse bound, not Mogul)",
          f"  {g['n_outliers']} outlier(s) of {len(g['bonds'])} unique bonds"]
    for bd in g["bonds"]:
        if bd["outlier"]:
            L.append(f"    {bd['a']}-{bd['b']}  {bd['len']:.3f} A  (covalent sum {bd['covsum']:.3f})  OUTLIER")
    L.append("Block C — H-bonds (D-H...A, X-H normalized; D...A class by Jeffrey; sym = operator on A)")
    if report["hbonds"]:
        L.append("  {:<6} {:<6} {:<6} {:>6} {:>6} {:>6} {:>6}  {:<9} {}".format(
            "D", "H", "A", "D-H", "H..A", "D..A", "ang", "class", "sym"))
        for b in report["hbonds"]:
            L.append("  {:<6} {:<6} {:<6} {:6.3f} {:6.3f} {:6.3f} {:6.1f}  {:<9} {}".format(
                b["D"], b["H"], b["A"], b["DH"], b["HA"], b["DA"], b["angle"], b["class"], b["sym"]))
        x = report["hbond_xref"]
        if x["rows"]:
            L.append(f"  (_geom_hbond cross-check: {x['matched']}/{x['rows']} matched on D..A)")
    else:
        L.append("  none detected")
    if report["contacts"]:
        L.append(f"  geometric contacts (below the {report['floor']:.0f} deg floor — listed, not asserted):")
        for c in report["contacts"]:
            L.append("    {:<6} {:<6} {:<6} D..A {:6.3f}  ang {:5.1f}  {}".format(
                c["D"], c["H"], c["A"], c["DA"], c["angle"], c["sym"]))
    return "\n".join(L)


def write_csv(report, path):
    """Write the H-bond table to CSV (the SI-ready artefact)."""
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["D", "H", "A", "D-H/A", "H..A/A", "D..A/A", "angle/deg", "class", "sym"])
        for b in report["hbonds"] + report["contacts"]:
            w.writerow([b["D"], b["H"], b["A"], b["DH"], b["HA"], b["DA"], b["angle"], b["class"], b["sym"]])
    return path
