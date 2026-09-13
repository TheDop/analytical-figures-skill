"""
crystal_view.py  -  deterministic 3D structure render. A *view* of the engine's validated
computation, never an independent artifact.

The make-or-break is camera orientation; the invariant is REPRODUCIBILITY (the camera lives
in cfg as numbers), not PCA specifically. PCA face-on is the auto default (+ the degeneracy
gates); custom / vector / axis angles are first-class because they're equally reproducible.

    complete_molecules(struct, cfg)   grow whole molecules across symmetry (don't orient a fragment)
    orient(struct, cfg, atoms)        -> (R, (elev,azim,roll), findings)  the deterministic camera
    render(struct, cfg)               -> (fig, ax)   element-coloured, bonds + dashed H-bonds + labels

Every bond and H-bond drawn here comes from crystal_engine.is_bonded / bond_pairs and
iter_hbond_candidates / hbond_geometry (the same predicates and traversal that build the
validation table, disorder exclusions included), applied through ONE crystal_engine.Supercell
per rendered cluster; nothing is re-derived in this module.
"""
from __future__ import annotations
import math
import numpy as np
from . import crystal_engine, style, verify
from .crystal_engine import Supercell, is_bonded, iter_hbond_candidates, hbond_sets

# CPK / element colours (recognizable); labels carry identity so colour is redundant.
CPK = {"H": "#D0D0D0", "C": "#404040", "N": "#3050F8", "O": "#FF2010", "S": "#E0C020",
       "F": "#80D060", "Cl": "#20C020", "Br": "#A62929", "I": "#8000C0", "P": "#FF8000",
       "Fe": "#E06633"}


def _color(sym):
    return CPK.get(sym, "#B000B0")


def _bond_color(sym):
    # a bond half needs a VISIBLE colour: H's CPK near-white (#D0D0D0) disappears on a
    # white page and makes O-H/N-H bonds look "half broken" -> use a mid grey for H bonds
    # (the atom sphere stays light CPK with a black outline).
    return "#9A9A9A" if sym == "H" else _color(sym)


def _size(sym):
    try:
        r = crystal_engine._cov_r(sym) or 0.7
    except Exception:
        r = 0.7
    return float(40.0 * r * r)             # area ~ radius^2, readable at print size


# ----------------------------------------------------------- molecule completion
def complete_molecules(struct, cfg):
    """Grow the asymmetric-unit fragments into whole molecules by following covalent bonds
    across a 3x3x3 neighbourhood (so an inversion-centre half-molecule is completed, and a
    boundary-straddling molecule is made whole) — PCA must orient a real molecule, not an
    AU fragment. Bonds by crystal_engine.is_bonded (disorder-aware) over the cached supercell."""
    sup = crystal_engine.supercell(struct)

    def k(p):
        return (int(round(p[0] * 50)), int(round(p[1] * 50)), int(round(p[2] * 50)))

    chosen, have = [], set()
    for at in struct.atoms:                       # seeds: one image of each AU atom
        rec = {"sym": at["sym"], "label": at["label"], "occ": at["occ"],
               "xyz": struct.cart(at["frac"])}
        if k(rec["xyz"]) not in have:
            have.add(k(rec["xyz"]))
            chosen.append(rec)
    frontier = list(chosen)
    while frontier:
        nxt = []
        for a in frontier:
            for j, _d in sup.bonded_to(a):
                kk = k(sup.xyz[j])
                if kk in have:
                    continue
                have.add(kk)
                rec = sup.rec(j)
                chosen.append(rec)
                nxt.append(rec)
        frontier = nxt
    return chosen


# ------------------------------------------------------- the rendered cluster as a block
def _cluster(atoms, alt=None):
    """The rendered cluster as ONE crystal_engine.Supercell (KD-trees, disorder map). `alt` is
    the disorder-alternative map — or an already-built Supercell of these atoms, which the
    renderers pass so every search in one render shares a single block."""
    return alt if isinstance(alt, Supercell) else Supercell.from_atoms(atoms, alt)


def _hide_ch(atoms, alt=None):
    """Drop hydrogens whose nearest heavy atom is a bonded carbon (view_hide_ch)."""
    blk = _cluster(atoms, alt)
    out = []
    for a in atoms:
        if a["sym"] == "H":
            k = blk.nearest_heavy(a["xyz"])
            if k is not None and blk.sym[k] == "C" and is_bonded(
                    "C", "H", float(np.linalg.norm(blk.xyz[k] - a["xyz"])), blk.occ[k], a["occ"],
                    blk.label[k], a["label"], blk.alt):
                continue
        out.append(a)
    return out


def _bond_pairs(atoms, alt=None):
    """(i, j) covalent bonds within the rendered cluster, by crystal_engine.bond_pairs."""
    return [(i, j) for i, j, _d in _cluster(atoms, alt).bond_pairs()]


def _hbond_records(atoms, cfg, alt=None):
    """(h_idx, donor_idx, acceptor_idx) for every asserted D-H...A within the rendered cluster,
    by crystal_engine.iter_hbond_candidates — the table's traversal and criteria, disorder
    exclusion included, so the figure can't assert a bond the table wouldn't. The renderers
    compute this ONCE and derive the dashed lines, the roll and the label set from it."""
    blk = _cluster(atoms, alt)
    idx = {id(a): i for i, a in enumerate(atoms)}
    return [(idx[id(h)], k, j) for h, k, j, geo in
            iter_hbond_candidates(blk, (a for a in atoms if a["sym"] == "H"), cfg)
            if geo["kind"] == "hbond"]


def _hbond_pairs(atoms, cfg, alt=None):
    """(hydrogen_idx, acceptor_idx) for the dashed H...A lines (see _hbond_records)."""
    return [(hi, ai) for hi, _di, ai in _hbond_records(atoms, cfg, alt)]


def _label_set(records):
    return {i for _hi, di, ai in records for i in (di, ai)}


# ------------------------------------------------------------------ orientation
def _view_axis(struct, cfg, P):
    """The line-of-sight axis (unit, cartesian) for the chosen mode, + PCA eigenvalues for
    the shape-degeneracy flag."""
    Pc = P - P.mean(0)
    evals, evecs = np.linalg.eigh(Pc.T @ Pc)          # ascending
    mode = cfg.view_orientation
    if mode == "pca":
        return evecs[:, 0], evals
    if mode in ("axis_a", "axis_b", "axis_c"):
        uvw = {"axis_a": (1, 0, 0), "axis_b": (0, 1, 0), "axis_c": (0, 0, 1)}[mode]
    elif mode == "vector":
        uvw = cfg.view_vector or (0, 0, 1)
    else:
        return evecs[:, 0], evals                      # custom handled separately
    d = struct.cart(np.array(uvw, float)) - struct.cart(np.zeros(3))
    n = np.linalg.norm(d)
    return (d / n if n else evecs[:, 0]), evals


def orient(struct, cfg, atoms, hbonds=None):
    """Return (R, (elev,azim,roll), findings). For pca/axis/vector R rotates coords into a
    view frame (z = line of sight, x = largest in-plane spread) and the camera looks straight
    down z; the in-plane roll aligns the H-bond network (or long axis) horizontal. For
    'custom' R is None and the given angles are used. view_tilt nudges any base. det(R)=+1
    keeps a proper rotation (no mirror -> correct enantiomorph). `hbonds` = the cluster's
    _hbond_pairs when the caller already has them (renders compute them once)."""
    f = []
    heavy = [a for a in atoms if a["sym"] != "H"]
    P = np.array([a["xyz"] for a in (heavy or atoms)])
    de, da, dr = cfg.view_tilt

    if cfg.view_orientation == "custom":
        ang = cfg.view_angles or (90.0, -90.0, 0.0)
        f.append(("INFO", f"custom camera elev/azim/roll={ang} (+tilt {cfg.view_tilt})"))
        return None, (ang[0] + de, ang[1] + da, ang[2] + dr), f

    z, evals = _view_axis(struct, cfg, P)
    z = z / (np.linalg.norm(z) or 1.0)
    # in-plane axes: PCA of coords projected off z (x = largest in-plane spread)
    Pc = P - P.mean(0)
    proj = Pc - np.outer(Pc @ z, z)
    ev2, evec2 = np.linalg.eigh(proj.T @ proj)
    x = evec2[:, -1] - (evec2[:, -1] @ z) * z
    x = x / (np.linalg.norm(x) or 1.0)
    y = np.cross(z, x)

    # roll: align the mean in-plane D->A H-bond vector horizontal (else long axis horizontal)
    theta = 0.0
    if cfg.view_roll_objective == "hbond":
        hb = hbonds if hbonds is not None else _hbond_pairs(atoms, cfg)
        vecs = [atoms[a]["xyz"] - atoms[d]["xyz"] for d, a in hb]
        if vecs:
            inplane = np.array([[v @ x, v @ y] for v in vecs])
            mean = inplane.sum(0)
            if np.linalg.norm(mean) > 0.3 * np.abs(inplane).sum(0).max():
                theta = -math.atan2(mean[1], mean[0])
            else:
                f.append(("INFO", "H-bond roll under-determined (vectors cancel/out-of-plane) "
                                  "-> long-axis horizontal"))
    if theta:
        xr = math.cos(theta) * x + math.sin(theta) * y
        yr = -math.sin(theta) * x + math.cos(theta) * y
        x, y = xr, yr

    R = np.array([x, y, z])
    if np.linalg.det(R) < 0:                     # keep a proper rotation (no mirror)
        R = np.array([x, -y, z])

    # view_tilt is folded INTO R, not added to the returned angles.
    # Why: the PyVista renderer (the default) builds its camera from R alone and DISCARDS the
    # angles - `R, _, findings = orient(...)`. Adding the tilt to the angles therefore did
    # nothing at all on the default backend, silently, while appearing to work in the matplotlib
    # fallback. Folding it into R makes both backends honour it, and the angles below stay at the
    # base (90, -90, 0) so matplotlib does not apply the same nudge twice.
    if (de, da, dr) != (0.0, 0.0, 0.0):
        ce, se = math.cos(math.radians(de)), math.sin(math.radians(de))
        ca, sa = math.cos(math.radians(da)), math.sin(math.radians(da))
        cr, sr = math.cos(math.radians(dr)), math.sin(math.radians(dr))
        Rz = np.array([[cr, -sr, 0.0], [sr, cr, 0.0], [0.0, 0.0, 1.0]])   # roll, about the sight line
        Rx = np.array([[1.0, 0.0, 0.0], [0.0, ce, -se], [0.0, se, ce]])   # elevation
        Ry = np.array([[ca, 0.0, sa], [0.0, 1.0, 0.0], [-sa, 0.0, ca]])   # azimuth
        R = Rz @ Rx @ Ry @ R                     # view-frame nudge, applied after the base

    # shape-degeneracy flag (only meaningful for the auto PCA view)
    if cfg.view_orientation == "pca" and evals[2] > 0:
        w2w1 = evals[1] / evals[2]
        if w2w1 >= 0.85:
            f.append(("WARN", f"orientation degenerate (w2/w1={w2w1:.2f} >= 0.85): no unique "
                              f"plane; PCA view is a default, not a derived answer"))
    f.append(("INFO", f"orientation={cfg.view_orientation} (+tilt {cfg.view_tilt})"))
    return R, (90.0, -90.0, 0.0), f                # tilt is already in R (see above)


# ----------------------------------------------------------------------- render
def _render_matplotlib(struct, cfg, atoms=None, cell_box=False):
    """Zero-dependency FALLBACK renderer (matplotlib 3D). matplotlib has no depth buffer,
    so atom/bond occlusion at vertices is imperfect (small white wedges) — prefer the
    pyvista backend for publication output. Returns (fig, ax)."""
    atoms, blk, hbonds, hbset = _prepare(struct, cfg, atoms)
    R, view, findings = orient(struct, cfg, atoms, hbonds)
    verify._resolve(findings, cfg, "orientation")
    desat = _desat_mask(atoms, blk) if cfg.color_by_component else {}

    P = np.array([a["xyz"] for a in atoms])
    Pc = P - P.mean(0)
    Pr = (R @ Pc.T).T if R is not None else Pc

    style.apply_style(cfg)
    fig, ax = style.figure(cfg, subplot_kw={"projection": "3d"})
    for i, j in _bond_pairs(atoms, blk):
        ci = _desat(_rgb(_bond_color(atoms[i]["sym"]))) if desat.get(id(atoms[i])) else _bond_color(atoms[i]["sym"])
        cj = _desat(_rgb(_bond_color(atoms[j]["sym"]))) if desat.get(id(atoms[j])) else _bond_color(atoms[j]["sym"])
        full = np.array([Pr[i], Pr[j]])
        # ONE continuous underlay -> no midpoint seam (the seam was showing the white page
        # through same-colour C-C bonds, splitting each aromatic bond in two).
        ax.plot(full[:, 0], full[:, 1], full[:, 2], color=cj, lw=2.0,
                solid_capstyle="round", zorder=1)
        if ci != cj:                                          # two-tone: overlay the i-half
            mid = (Pr[i] + Pr[j]) / 2
            seg = np.array([Pr[i], mid])
            ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color=ci, lw=2.0,
                    solid_capstyle="round", zorder=1)
    for hi, ai in hbonds:                                    # dashed H...A (from the hydrogen)
        seg = np.array([Pr[hi], Pr[ai]])
        ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color="0.25", lw=0.9, ls=(0, (4, 3)), zorder=2)
    for a, p in zip(atoms, Pr):
        col = _desat(_rgb(_color(a["sym"]))) if desat.get(id(a)) else _color(a["sym"])
        ax.scatter(p[0], p[1], p[2], s=_size(a["sym"]), color=col,
                   edgecolors="k", linewidths=0.3, depthshade=True, zorder=3)
    if cfg.view_label_atoms != "none":
        labelled = set()
        for i, (a, p) in enumerate(zip(atoms, Pr)):
            sup = a.get("sup")
            if a["occ"] < 0.99 or (a["label"], sup) in labelled:   # skip disordered + dedupe (label,sym-image)
                continue
            if (cfg.view_label_atoms == "all"
                    or (cfg.view_label_atoms == "hbond" and i in hbset)
                    or (cfg.view_label_atoms not in ("all", "hbond") and a["sym"] not in ("C", "H"))):
                from matplotlib import patheffects as _pe
                txt = f"$\\mathrm{{{a['label']}}}^{{\\mathrm{{{sup}}}}}$" if sup else a["label"]
                ax.text(p[0], p[1], p[2], txt, color="black", zorder=4,
                        fontsize=6.0 * cfg.view_label_size,
                        path_effects=[_pe.withStroke(linewidth=2.5, foreground="white")])
                labelled.add((a["label"], sup))

    if cell_box:
        c0 = P.mean(0)
        for e0, e1 in _cell_edges(struct):
            seg = np.array([e0, e1]) - c0
            seg = (R @ seg.T).T if R is not None else seg
            ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color="0.15", lw=1.0, zorder=2)

    # tight, undistorted framing: equal cubic limits (no axis stretched) + fill the frame
    ax.set_axis_off()
    mid = (Pr.max(0) + Pr.min(0)) / 2.0
    half = (float((Pr.max(0) - Pr.min(0)).max()) / 2.0 * 1.08) or 1.0
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)
    try:
        ax.set_box_aspect((1, 1, 1), zoom=1.35)       # zoom: matplotlib >= 3.8
    except TypeError:
        ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=view[0], azim=view[1], roll=view[2])
    return fig, ax


def _rgb(hexstr):
    h = hexstr.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _choose_backend(cfg):
    want = cfg.view_renderer
    if want == "matplotlib":
        return "matplotlib"
    try:
        import pyvista  # noqa: F401
        return "pyvista"
    except ImportError:
        if want == "pyvista":
            raise ImportError("Required package not found: pyvista. Install with:\n"
                              "    python -m pip install pyvista\n"
                              "(or set cfg.view_renderer='matplotlib' for the no-dep fallback).")
        return "matplotlib"


def _atom_op(struct, label, xyz_cart):
    """Recover the symmetry op that generated a rendered atom at xyz_cart (cartesian), so its
    ADP can be rotated to match (U_image = R U R^T). None if not found (treat as identity)."""
    asym = next((a for a in struct.atoms if a["label"] == label), None)
    if asym is None:
        return None
    g = crystal_engine._gemmi()
    fr = struct.cell.fractionalize(g.Position(float(xyz_cart[0]), float(xyz_cart[1]), float(xyz_cart[2])))
    target = np.array([fr.x, fr.y, fr.z]) % 1.0
    for op in struct.ops:
        gp = np.array(op.apply_to_xyz(list(asym["frac"]))) % 1.0
        if np.all(np.abs((gp - target + 0.5) % 1.0 - 0.5) < 5e-3):
            return op
    return None


def _overlay_labels(pl, img, sel, res, cfg):
    """Composite atom labels as HALO text (black glyphs + white outline) at each atom's
    projected pixel — legible over dark atoms WITHOUT a filled box that would hide the
    molecule (VTK's own labels can't draw a text halo). VTK world->display gives the pixel."""
    from PIL import Image, ImageDraw, ImageFont
    import matplotlib.font_manager as fm
    try:
        import vtk
        Coord = vtk.vtkCoordinate
    except Exception:
        from vtkmodules.vtkRenderingCore import vtkCoordinate as Coord
    coord = Coord(); coord.SetCoordinateSystemToWorld()
    ren = pl.renderer
    arr = img[..., :3].copy(); H = arr.shape[0]
    pim = Image.fromarray(arr); draw = ImageDraw.Draw(pim)
    fpx = max(12, int(res / 55 * cfg.view_label_size))
    try:
        font = ImageFont.truetype(fm.findfont("DejaVu Sans"), fpx)
    except Exception:
        font = ImageFont.load_default()
    sw = max(2, fpx // 7)
    sub_fpx = max(9, int(fpx * 0.66))                    # symmetry-superscript size
    try:
        subfont = ImageFont.truetype(fm.findfont("DejaVu Sans"), sub_fpx)
    except Exception:
        subfont = font
    gap = 0.6 * fpx if cfg.view_label_offset else 0.0
    anchor = "lm" if gap else "mm"
    for ent in sel:
        lab, p = ent[0], ent[1]
        sup = ent[2] if len(ent) > 2 else None
        coord.SetValue(float(p[0]), float(p[1]), float(p[2]))
        dx, dy = coord.GetComputedDisplayValue(ren)
        x, y = dx + gap, H - dy
        draw.text((x, y), lab, font=font, fill=(0, 0, 0), anchor=anchor,
                  stroke_width=sw, stroke_fill=(255, 255, 255))
        if sup:                                          # raised superscript at the label's top-right
            bb = draw.textbbox((x, y), lab, font=font, anchor=anchor)
            draw.text((bb[2] + sub_fpx * 0.1, bb[1] - sub_fpx * 0.15), sup, font=subfont,
                      fill=(0, 0, 0), anchor="lm", stroke_width=max(1, sw // 2),
                      stroke_fill=(255, 255, 255))
    return np.asarray(pim)


def _desat(rgb, f=0.72):
    """Luminance-preserving desaturation (blend toward the colour's own grey) — mutes a
    cocrystal's secondary component without changing its brightness."""
    lum = 0.30 * rgb[0] + 0.59 * rgb[1] + 0.11 * rgb[2]
    return tuple(c * (1 - f) + lum * f for c in rgb)


def _desat_mask(atoms, alt):
    """{id(atom): True} for atoms NOT in the largest molecule — so color_by_component keeps the
    largest component (the 'main' molecule) in full element colour and mutes the rest (the
    coformer/solvent). Empty (no muting) when there's a single component."""
    comps = _component_indices(_cluster(atoms, alt))
    if len(comps) < 2:
        return {}
    mx = max(len(c) for c in comps)
    return {id(atoms[i]): len(c) < mx for c in comps for i in c}


def _prepare(struct, cfg, atoms):
    """The shared front half of both renderers: the cluster (default: the completed molecules),
    C-H hiding, ONE Supercell over it, and from that block — once — the dashed H-bond pairs
    and the synthon label set. Returns (atoms, blk, hbonds, hbset)."""
    alt = crystal_engine.disorder_alternatives(struct)
    if atoms is None:
        atoms = complete_molecules(struct, cfg)
    if cfg.view_hide_ch:
        atoms = _hide_ch(atoms, alt)
    if len(atoms) < 2:
        raise ValueError("crystal_view.render: fewer than 2 atoms after completion")
    blk = Supercell.from_atoms(atoms, alt)
    records = _hbond_records(atoms, cfg, blk)
    hbonds = [(hi, ai) for hi, _di, ai in records]
    hbset = _label_set(records) if cfg.view_label_atoms == "hbond" else None
    return atoms, blk, hbonds, hbset


def _render_pyvista(struct, cfg, atoms=None, cell_box=False):
    """HIGH-QUALITY offscreen 3D render (VTK): smooth-shaded spheres + cylinder bonds with a
    real depth buffer (correct occlusion -> NO vertex gaps), SSAA + SSAO, orthographic camera
    from the orientation engine. `atoms` overrides the default complete-molecule set; `cell_box`
    draws the unit-cell edges + a/b/c. Returns an RGB ndarray with halo labels composited."""
    import pyvista as pv
    atoms, blk, hbonds, hbset = _prepare(struct, cfg, atoms)
    R, _, findings = orient(struct, cfg, atoms, hbonds)
    P = np.array([a["xyz"] for a in atoms])
    ctr = P.mean(0)
    if R is None:                          # 'custom' angles target the matplotlib backend
        heavy = [a for a in atoms if a["sym"] != "H"]
        Q = np.array([a["xyz"] for a in (heavy or atoms)])
        Qc = Q - Q.mean(0)
        _, evec = np.linalg.eigh(Qc.T @ Qc)
        x = evec[:, 2]; z = evec[:, 0]; y = np.cross(z, x)
        R = np.array([x, y, z])
        if np.linalg.det(R) < 0:
            R = np.array([x, -y, z])
        findings.append(("INFO", "custom view_angles apply to the matplotlib backend; "
                                 "pyvista used a PCA camera (use view_orientation='vector' "
                                 "for an explicit pyvista view)"))
    verify._resolve(findings, cfg, "orientation")
    desat = _desat_mask(atoms, blk) if cfg.color_by_component else {}

    res = int(cfg.view_resolution)
    pl = pv.Plotter(off_screen=True, window_size=[res, res], lighting="light kit")
    pl.set_background("white")

    def crad(sym):
        try:
            r = crystal_engine._cov_r(sym) or 0.7
        except Exception:
            r = 0.7
        return 0.30 * (r / 0.7) + 0.08              # ball-and-stick: balls < vdW, H smaller

    sph = dict(smooth_shading=True, specular=0.3, specular_power=15)
    ellipsoid_mode = False
    if cfg.view_style == "ellipsoid":
        if crystal_engine.has_adp(struct):
            ellipsoid_mode = True
        else:
            verify._resolve([("WARN", "view_style='ellipsoid' but no anisotropic U in the CIF "
                                      "(_atom_site_aniso_U_*) — drawing ball-and-stick")], cfg, "adp")
    prob = cfg.adp_probability
    for a, p in zip(atoms, P):
        col = _rgb(_color(a["sym"]))
        if desat.get(id(a)):
            col = _desat(col)
        if ellipsoid_mode and a["sym"] != "H" and a["label"] in struct.aniso:
            Uc = crystal_engine.u_cart(struct, struct.aniso[a["label"]])
            op = _atom_op(struct, a["label"], p)
            if op is not None:                                    # rotate ADP for sym images
                Rc = crystal_engine.op_rot_cart(struct, op)
                Uc = Rc @ Uc @ Rc.T
            w, Vv = np.linalg.eigh(Uc)
            if np.all(w > 1e-9):                                  # positive-definite -> ellipsoid
                semi = crystal_engine._adp_scale(prob) * np.sqrt(w)
                T = np.eye(4); T[:3, :3] = Vv @ np.diag(semi); T[:3, 3] = p
                m = pv.Sphere(radius=1.0, theta_resolution=36, phi_resolution=36)
                m.transform(T, inplace=True)
                pl.add_mesh(m, color=col, **sph)
                continue                                          # else fall through to a sphere
        pl.add_mesh(pv.Sphere(radius=crad(a["sym"]), center=p, theta_resolution=48,
                              phi_resolution=48), color=col, **sph)
    for i, j in _bond_pairs(atoms, blk):
        a_, b_ = P[i], P[j]; d = b_ - a_; L = float(np.linalg.norm(d)); mid = (a_ + b_) / 2
        ri, rj = _rgb(_bond_color(atoms[i]["sym"])), _rgb(_bond_color(atoms[j]["sym"]))
        if desat.get(id(atoms[i])):
            ri = _desat(ri)
        if desat.get(id(atoms[j])):
            rj = _desat(rj)
        if ri == rj:
            pl.add_mesh(pv.Cylinder(center=mid, direction=d, radius=0.11, height=L, resolution=28),
                        color=ri, **sph)
        else:
            for c, q in ((ri, (a_ + mid) / 2), (rj, (b_ + mid) / 2)):
                pl.add_mesh(pv.Cylinder(center=q, direction=d, radius=0.11, height=L / 2, resolution=28),
                            color=c, **sph)
    for hi, ai in hbonds:                            # dashed H...A (originates at the hydrogen)
        a_, b_ = P[hi], P[ai]; d = b_ - a_; L = float(np.linalg.norm(d))
        n = max(3, int(L / 0.35))
        for t in range(0, n, 2):
            c0 = a_ + d * (t / n); c1 = a_ + d * (min(t + 1, n) / n)
            pl.add_mesh(pv.Cylinder(center=(c0 + c1) / 2, direction=d, radius=0.045,
                                    height=float(np.linalg.norm(c1 - c0)), resolution=12),
                        color=(0.35, 0.35, 0.35))
    # collect heteroatom labels (skip disordered + dedupe); composited as HALO text after
    # the render (below) so they sit at the atom without a box hiding the molecule.
    sel, seen = [], set()
    if cfg.view_label_atoms != "none":
        for i, (a, p) in enumerate(zip(atoms, P)):
            sup = a.get("sup")
            dk = (a["label"], sup)
            if a["occ"] < 0.99 or dk in seen:
                continue
            if (cfg.view_label_atoms == "all"
                    or (cfg.view_label_atoms == "hbond" and i in hbset)
                    or (cfg.view_label_atoms not in ("all", "hbond") and a["sym"] not in ("C", "H"))):
                sel.append((a["label"], p, sup)); seen.add(dk)

    if cell_box:
        for p0, p1 in _cell_edges(struct):
            pl.add_mesh(pv.Line(p0, p1), color=(0.15, 0.15, 0.15), line_width=3)
        sel += [("a", struct.cart([1.08, 0, 0]), None), ("b", struct.cart([0, 1.08, 0]), None),
                ("c", struct.cart([0, 0, 1.08]), None)]

    D = float(np.ptp(P, axis=0).max()) * 3 + 5
    pl.enable_parallel_projection()
    pl.camera_position = [tuple(ctr + R[2] * D), tuple(ctr), tuple(R[1])]
    pl.reset_camera()
    try:
        pl.enable_anti_aliasing("ssaa")
    except Exception:
        pass
    if cfg.view_ssao:
        try:
            pl.enable_ssao()
        except Exception:
            pass

    img = pl.screenshot(return_img=True)
    if sel:
        img = _overlay_labels(pl, img, sel, res, cfg)
    pl.close()
    return img


def render(struct, cfg, atoms=None, cell_box=False):
    """Render a structure view. Returns (backend, obj): ('pyvista', RGB ndarray) or
    ('matplotlib', (fig, ax)). `atoms` overrides the default complete-molecule set (the
    unit-cell / packing / H-bond-environment views pass their own); `cell_box` draws the
    unit-cell edges + a/b/c. PyVista default, matplotlib fallback. Write with save()."""
    backend = _choose_backend(cfg)
    if backend == "pyvista":
        return "pyvista", _render_pyvista(struct, cfg, atoms, cell_box)
    return "matplotlib", _render_matplotlib(struct, cfg, atoms, cell_box)


def save(rendered, basename, cfg):
    """Write the structure image. pyvista -> a high-res PNG (raster is the convention for
    structure images); matplotlib -> style.save_fig (vector + raster). Returns the paths
    written (open the PNG with the Read tool for the perceptual QA pass)."""
    backend, obj = rendered
    if backend == "pyvista":
        import os
        from PIL import Image
        png = basename + ".png"
        os.makedirs(os.path.dirname(png) or ".", exist_ok=True)
        Image.fromarray(obj).save(png)            # obj is the composited RGB ndarray
        return [png]
    fig, _ = obj
    return style.save_fig(fig, basename, cfg)


# ------------------------------------------------------------------- packing (Phase 2)
def _component_indices(blk):
    """Connected components (whole molecules) of a block by its covalent bonds (crystal_engine
    .bond_pairs, disorder-aware): lists of block indices, each ascending, ordered by first atom."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    n = len(blk.xyz)
    if n == 0:
        return []
    pairs = blk.bond_pairs()
    i = np.array([p[0] for p in pairs], int)
    j = np.array([p[1] for p in pairs], int)
    adj = coo_matrix((np.ones(len(i), bool), (i, j)), shape=(n, n))
    _, lab = connected_components(adj, directed=False)
    order = np.argsort(lab, kind="stable")
    bounds = np.flatnonzero(np.diff(lab[order])) + 1
    return [c.tolist() for c in np.split(order, bounds)]


def _components(atoms, alt=None):
    """Connected components (whole molecules) of `atoms`, as lists of atom dicts."""
    return [[atoms[i] for i in c] for c in _component_indices(_cluster(atoms, alt))]


def _pack_atoms(struct, cfg, cells=None):
    """WHOLE molecules whose CENTROID lies inside the (Nx,Ny,Nz) block — so a single cell shows
    exactly its Z formula units (each molecule ONCE), not every boundary fragment grown into a
    duplicate. Built over a -1..N+1 supercell, grouped into molecules, then centroid-filtered."""
    nx, ny, nz = cells or cfg.pack_cells
    if cfg.cell_fill == "clip":   # literal cell contents, molecules cut at the outer box
        return crystal_engine.supercell(struct, (0, 0, 0), (nx - 1, ny - 1, nz - 1)).records()
    blk = crystal_engine.supercell(struct, (-1, -1, -1), (nx, ny, nz))
    out = []
    for comp in _component_indices(blk):
        cen = blk.frac[comp].mean(axis=0)
        if (-1e-4 <= cen[0] < nx) and (-1e-4 <= cen[1] < ny) and (-1e-4 <= cen[2] < nz):
            out.extend(blk.rec(i) for i in comp)
    return out


def _cell_edges(struct):
    c = {(i, j, k): struct.cart([i, j, k]) for i in (0, 1) for j in (0, 1) for k in (0, 1)}
    edges = []
    for corner in c:
        for ax in range(3):
            if corner[ax] == 0:
                nb = list(corner); nb[ax] = 1
                edges.append((c[corner], c[tuple(nb)]))
    return edges


def _roman(n):
    out, vals = "", [(10, "x"), (9, "ix"), (5, "v"), (4, "iv"), (1, "i")]
    for v, s in vals:
        while n >= v:
            out += s; n -= v
    return out


def _sym_op_text(struct, code):
    """Human-readable symmetry op for an 'n_pqr' code (e.g. '-x+1, -y+1, -z+1') for the caption key."""
    g = crystal_engine._gemmi()
    try:
        op = struct.ops[int(code.split("_")[0]) - 1]
        new = g.Op(op.triplet())
        if "_" in code and len(code.split("_")[1]) == 3:
            t = code.split("_")[1]
            tr = list(new.tran)                          # gemmi returns a copy; assign the whole list back
            for i in range(3):
                tr[i] += (int(t[i]) - 5) * new.DEN
            new.tran = tr
        return new.triplet().replace(",", ", ")
    except Exception:
        return code


def _hbond_env_atoms(struct, cfg, central=None):
    """The asymmetric-unit molecule(s) + the symmetry neighbours they hydrogen-bond to, at the
    extent set by cfg.hbond_neighbour: 'whole' (full neighbour molecule), 'stub' (contact atom +
    one bonded shell), or 'site' (contact atom only). Symmetry-generated neighbour atoms are
    annotated with their 'n_pqr' code + a roman superscript, and the caption key is logged.
    `central` defaults to the whole asymmetric unit; pass a subset (e.g. one molecule) to CROP the
    figure to just that molecule's H-bond environment (useful for Z'>1 structures)."""
    if central is None:
        central = complete_molecules(struct, cfg)
    for a in central:
        a["neighbour"] = False
    donors, acc = hbond_sets(cfg)
    sup = crystal_engine.supercell(struct)

    def key(p):
        return (int(round(p[0] * 50)), int(round(p[1] * 50)), int(round(p[2] * 50)))

    central_keys = {key(a["xyz"]) for a in central}
    cacc_keys = {key(a["xyz"]) for a in central if a["sym"] in acc}

    # Both directions run the engine's traversal over the 3x3x3 block: central hydrogens seed
    # OUTSIDE acceptors; the hydrogens of every donor atom within the 4.0 A D...A ceiling of a
    # central acceptor seed their OUTSIDE donor. (A central molecule's own atoms are excluded by
    # position key, so intramolecular contacts never become neighbours.)
    seeds = []
    for h, k, j, geo in iter_hbond_candidates(sup, (a for a in central if a["sym"] == "H"), cfg):
        if geo["kind"] == "hbond" and key(sup.xyz[j]) not in central_keys:
            seeds.append(sup.rec(j))
    cand = set()
    for a in central:
        if a["sym"] in acc:
            for k in sup.neighbours(a["xyz"], 4.0):
                if sup.sym[k] in donors and key(sup.xyz[k]) not in central_keys:
                    cand.update(j for j, _d in sup.bonded_to(sup.atom(k)) if sup.sym[j] == "H")
    for h, k, j, geo in iter_hbond_candidates(sup, (sup.atom(i) for i in sorted(cand)), cfg):
        if geo["kind"] == "hbond" and key(sup.xyz[j]) in cacc_keys:
            seeds.append(sup.rec(k))

    # ---- assemble neighbour atoms at the requested extent
    mode = cfg.hbond_neighbour
    seen_seed, seed_list = set(central_keys), []
    for s_ in seeds:                                  # unique contact atoms (the seeds)
        kk = key(s_["xyz"])
        if kk not in seen_seed:
            seen_seed.add(kk); seed_list.append(s_)
    have, neigh = set(central_keys), []
    for s_ in seed_list:
        have.add(key(s_["xyz"])); neigh.append(s_)
    if mode in ("stub", "whole"):                     # grow outward (1 shell for stub, fully for whole)
        depth = {id(s_): 0 for s_ in seed_list}
        frontier = list(seed_list)
        while frontier:
            nf = []
            for a in frontier:
                if mode == "stub" and depth[id(a)] >= 1:
                    continue
                for j, _d in sup.bonded_to(a):
                    kk = key(sup.xyz[j])
                    if kk in have:
                        continue
                    s_ = sup.rec(j)
                    have.add(kk); neigh.append(s_); nf.append(s_); depth[id(s_)] = depth[id(a)] + 1
            frontier = nf

    # ---- annotate neighbours: symmetry code + roman superscript; log the caption key
    codes = []
    for a in neigh:
        a["neighbour"] = True
        a["symcode"] = crystal_engine._sym_code(struct, a["label"], a["frac"])
        if a["symcode"] and a["symcode"] != "." and a["symcode"] not in codes:
            codes.append(a["symcode"])
    roman = {c: _roman(i + 1) for i, c in enumerate(codes)}
    for a in neigh:
        a["sup"] = roman.get(a.get("symcode"))
    if codes:
        keytxt = "; ".join(f"({roman[c]}) {_sym_op_text(struct, c)}" for c in codes)
        verify._resolve([("INFO", f"symmetry key for the caption: {keytxt}")], cfg, "hbond_env")
    return central + neigh


def render_unit_cell(struct, cfg):
    """Contents of a SINGLE unit cell (whole molecules) + the unit-cell box and a/b/c axes."""
    return render(struct, cfg, atoms=_pack_atoms(struct, cfg, cells=(1, 1, 1)), cell_box=True)


def render_packing(struct, cfg):
    """Packing diagram: cfg.pack_cells unit cells of whole molecules with H-bonds dashed, plus
    the unit-cell box and a/b/c axes. Viewing down a cell axis (cfg.view_orientation='axis_a'|
    'axis_b'|'axis_c') usually reads best."""
    return render(struct, cfg, atoms=_pack_atoms(struct, cfg), cell_box=True)


def render_hbond_environment(struct, cfg):
    """Asymmetric-unit molecule(s) + the symmetry neighbours they hydrogen-bond to, with the
    intermolecular H-bonds dashed — the H-bonding environment."""
    return render(struct, cfg, atoms=_hbond_env_atoms(struct, cfg), cell_box=False)
