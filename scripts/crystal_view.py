"""
crystal_view.py  -  deterministic 3D structure render. A *view* of the engine's validated
computation, never an independent artifact.

The make-or-break is camera orientation; the invariant is REPRODUCIBILITY (the camera lives
in cfg as numbers), not PCA specifically. PCA face-on is the auto default (+ the degeneracy
gates); custom / vector / axis angles are first-class because they're equally reproducible.

    complete_molecules(struct, cfg)   grow whole molecules across symmetry (don't orient a fragment)
    orient(struct, cfg, atoms)        -> (R, (elev,azim,roll), findings)  the deterministic camera
    render(struct, cfg)               -> (fig, ax)   element-coloured, bonds + dashed H-bonds + labels
"""
from __future__ import annotations
import math
import numpy as np
from . import crystal_engine, style, verify

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
        r = crystal_engine._gemmi().Element(sym).covalent_r or 0.7
    except Exception:
        r = 0.7
    return float(40.0 * r * r)             # area ~ radius^2, readable at print size


# ----------------------------------------------------------- molecule completion
def complete_molecules(struct, cfg):
    """Grow the asymmetric-unit fragments into whole molecules by following covalent bonds
    across a 3x3x3 neighbourhood (so an inversion-centre half-molecule is completed, and a
    boundary-straddling molecule is made whole) — PCA must orient a real molecule, not an
    AU fragment. Disorder-aware (mutually-exclusive partial sites aren't bonded)."""
    cell_atoms = [{**u, "xyz": struct.cart(u["frac"])} for u in crystal_engine.expand(struct)]
    alt = crystal_engine.disorder_alternatives(struct)
    origin = struct.cart(np.zeros(3))
    sup = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                shift = struct.cart(np.array([di, dj, dk], float)) - origin
                for a in cell_atoms:
                    sup.append({"sym": a["sym"], "label": a["label"], "occ": a["occ"],
                                "xyz": a["xyz"] + shift})

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
            for s in sup:
                kk = k(s["xyz"])
                if kk in have:
                    continue
                d = float(np.linalg.norm(a["xyz"] - s["xyz"]))
                if not (0.4 < d <= crystal_engine._covsum(a["sym"], s["sym"]) + crystal_engine.BOND_TOL):
                    continue
                if a["occ"] < 1 and s["occ"] < 1 and d < crystal_engine.DISORDER_MIN:
                    continue
                if s["label"] in alt.get(a["label"], set()):       # disorder-alternative -> not bonded
                    continue
                have.add(kk)
                chosen.append(s)
                nxt.append(s)
        frontier = nxt
    return chosen


def _hide_ch(atoms):
    heavy = [a for a in atoms if a["sym"] != "H"]
    out = []
    for a in atoms:
        if a["sym"] == "H" and heavy:
            D = min(heavy, key=lambda h: np.linalg.norm(h["xyz"] - a["xyz"]))
            if D["sym"] == "C" and np.linalg.norm(D["xyz"] - a["xyz"]) <= crystal_engine._covsum("C", "H") + crystal_engine.BOND_TOL:
                continue
        out.append(a)
    return out


def _bond_pairs(atoms, alt=None):
    xyz = [a["xyz"] for a in atoms]
    out = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if 0.4 < d <= crystal_engine._covsum(atoms[i]["sym"], atoms[j]["sym"]) + crystal_engine.BOND_TOL:
                if atoms[i]["occ"] < 1 and atoms[j]["occ"] < 1 and d < crystal_engine.DISORDER_MIN:
                    continue
                if alt and atoms[j]["label"] in alt.get(atoms[i]["label"], set()):
                    continue                                  # disorder alternative -> not bonded
                out.append((i, j))
    return out


def _hbond_pairs(atoms, cfg):
    """(donor_idx, acceptor_idx) for dashed lines, by the engine's criteria within the
    rendered cluster (same rule as crystal_engine.hbonds -> the figure can't assert a bond
    the table wouldn't)."""
    donors = set(cfg.hbond_donors) | ({"C"} if cfg.hbond_weak else set())
    acc = set(cfg.hbond_acceptors)
    heavy = [a for a in atoms if a["sym"] != "H"]
    out = []
    for hi, h in enumerate(atoms):
        if h["sym"] != "H" or not heavy:
            continue
        D = min(heavy, key=lambda a: np.linalg.norm(a["xyz"] - h["xyz"]))
        dDH = float(np.linalg.norm(D["xyz"] - h["xyz"]))
        if dDH > crystal_engine._covsum(D["sym"], "H") + crystal_engine.BOND_TOL or D["sym"] not in donors:
            continue
        hpos = h["xyz"]
        if cfg.xh_normalize and D["sym"] in crystal_engine.NEUTRON_XH and dDH > 0:
            hpos = D["xyz"] + (h["xyz"] - D["xyz"]) / dDH * crystal_engine.NEUTRON_XH[D["sym"]]
        for ai, a in enumerate(atoms):
            if a["sym"] not in acc:
                continue
            DA = float(np.linalg.norm(a["xyz"] - D["xyz"]))
            if DA < 0.4 or DA > 4.0:
                continue
            if np.linalg.norm(a["xyz"] - hpos) > crystal_engine._vdw("H") + crystal_engine._vdw(a["sym"]):
                continue
            if crystal_engine._angle(D["xyz"], hpos, a["xyz"]) < cfg.hbond_angle_min:
                continue
            out.append((hi, ai))     # dashed line starts at the HYDROGEN (H...A), not the donor O/N
    return out


def _hbond_label_set(atoms, cfg):
    """Atom indices to label under view_label_atoms='hbond': the DONOR heavy atom and the ACCEPTOR
    of each H-bond — i.e. only the atoms that make the synthon (the bridging H stays unlabelled).
    Empty set if the rendered cluster has no H-bonds."""
    idx = set()
    for hi, ai in _hbond_pairs(atoms, cfg):
        idx.add(ai)
        h = atoms[hi]["xyz"]
        donor = min((j for j, a in enumerate(atoms) if a["sym"] != "H"),
                    key=lambda j: float(np.linalg.norm(atoms[j]["xyz"] - h)), default=None)
        if donor is not None:
            idx.add(donor)
    return idx


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


def orient(struct, cfg, atoms):
    """Return (R, (elev,azim,roll), findings). For pca/axis/vector R rotates coords into a
    view frame (z = line of sight, x = largest in-plane spread) and the camera looks straight
    down z; the in-plane roll aligns the H-bond network (or long axis) horizontal. For
    'custom' R is None and the given angles are used. view_tilt nudges any base. det(R)=+1
    keeps a proper rotation (no mirror -> correct enantiomorph)."""
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
        hb = _hbond_pairs(atoms, cfg)
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
    if atoms is None:
        atoms = complete_molecules(struct, cfg)
    if cfg.view_hide_ch:
        atoms = _hide_ch(atoms)
    if len(atoms) < 2:
        raise ValueError("crystal_view.render: fewer than 2 atoms after completion")
    R, view, findings = orient(struct, cfg, atoms)
    verify._resolve(findings, cfg, "orientation")
    alt = crystal_engine.disorder_alternatives(struct)
    desat = _desat_mask(atoms, alt) if getattr(cfg, "color_by_component", False) else {}

    P = np.array([a["xyz"] for a in atoms])
    Pc = P - P.mean(0)
    Pr = (R @ Pc.T).T if R is not None else Pc

    style.apply_style(cfg)
    fig, ax = style.figure(cfg, subplot_kw={"projection": "3d"})
    bonds = _bond_pairs(atoms, alt)
    for i, j in bonds:
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
    for hi, ai in _hbond_pairs(atoms, cfg):                  # dashed H...A (from the hydrogen)
        seg = np.array([Pr[hi], Pr[ai]])
        ax.plot(seg[:, 0], seg[:, 1], seg[:, 2], color="0.25", lw=0.9, ls=(0, (4, 3)), zorder=2)
    for a, p in zip(atoms, Pr):
        col = _desat(_rgb(_color(a["sym"]))) if desat.get(id(a)) else _color(a["sym"])
        ax.scatter(p[0], p[1], p[2], s=_size(a["sym"]), color=col,
                   edgecolors="k", linewidths=0.3, depthshade=True, zorder=3)
    if cfg.view_label_atoms != "none":
        hbset = _hbond_label_set(atoms, cfg) if cfg.view_label_atoms == "hbond" else None
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
                        fontsize=6.0 * getattr(cfg, "view_label_size", 1.0),
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
    want = getattr(cfg, "view_renderer", "auto")
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
    fpx = max(12, int(res / 55 * getattr(cfg, "view_label_size", 1.0)))
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
    gap = 0.6 * fpx if getattr(cfg, "view_label_offset", False) else 0.0
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
    comps = _components(atoms, alt)
    if len(comps) < 2:
        return {}
    mx = max(len(c) for c in comps)
    mask = {}
    for c in comps:
        muted = len(c) < mx
        for a in c:
            mask[id(a)] = muted
    return mask


def _render_pyvista(struct, cfg, atoms=None, cell_box=False):
    """HIGH-QUALITY offscreen 3D render (VTK): smooth-shaded spheres + cylinder bonds with a
    real depth buffer (correct occlusion -> NO vertex gaps), SSAA + SSAO, orthographic camera
    from the orientation engine. `atoms` overrides the default complete-molecule set; `cell_box`
    draws the unit-cell edges + a/b/c. Returns an RGB ndarray with halo labels composited."""
    import pyvista as pv
    if atoms is None:
        atoms = complete_molecules(struct, cfg)
    if cfg.view_hide_ch:
        atoms = _hide_ch(atoms)
    if len(atoms) < 2:
        raise ValueError("crystal_view: fewer than 2 atoms after completion")
    R, _, findings = orient(struct, cfg, atoms)
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
    alt = crystal_engine.disorder_alternatives(struct)
    desat = _desat_mask(atoms, alt) if getattr(cfg, "color_by_component", False) else {}

    res = int(getattr(cfg, "view_resolution", 1800))
    pl = pv.Plotter(off_screen=True, window_size=[res, res], lighting="light kit")
    pl.set_background("white")

    def crad(sym):
        try:
            r = crystal_engine._gemmi().Element(sym).covalent_r or 0.7
        except Exception:
            r = 0.7
        return 0.30 * (r / 0.7) + 0.08              # ball-and-stick: balls < vdW, H smaller

    sph = dict(smooth_shading=True, specular=0.3, specular_power=15)
    ellipsoid_mode = False
    if getattr(cfg, "view_style", "ball_stick") == "ellipsoid":
        if crystal_engine.has_adp(struct):
            ellipsoid_mode = True
        else:
            verify._resolve([("WARN", "view_style='ellipsoid' but no anisotropic U in the CIF "
                                      "(_atom_site_aniso_U_*) — drawing ball-and-stick")], cfg, "adp")
    prob = getattr(cfg, "adp_probability", 0.5)
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
    for i, j in _bond_pairs(atoms, alt):
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
    for hi, ai in _hbond_pairs(atoms, cfg):          # dashed H...A (originates at the hydrogen)
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
        hbset = _hbond_label_set(atoms, cfg) if cfg.view_label_atoms == "hbond" else None
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
    if getattr(cfg, "view_ssao", True):
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
def _components(atoms, alt=None):
    """Connected components (whole molecules) of `atoms` by covalent bonds, disorder-aware."""
    from scipy.spatial import cKDTree
    n = len(atoms)
    if n == 0:
        return []
    xyz = np.array([a["xyz"] for a in atoms])
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    for i, j in cKDTree(xyz).query_pairs(r=2.6):
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        if d > crystal_engine._covsum(atoms[i]["sym"], atoms[j]["sym"]) + crystal_engine.BOND_TOL:
            continue
        if atoms[i]["occ"] < 1 and atoms[j]["occ"] < 1 and d < crystal_engine.DISORDER_MIN:
            continue
        if alt and atoms[j]["label"] in alt.get(atoms[i]["label"], set()):
            continue                                          # disorder alternative -> not bonded
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    comps = {}
    for idx in range(n):
        comps.setdefault(find(idx), []).append(atoms[idx])
    return list(comps.values())


def _pack_atoms(struct, cfg, cells=None):
    """WHOLE molecules whose CENTROID lies inside the (Nx,Ny,Nz) block — so a single cell shows
    exactly its Z formula units (each molecule ONCE), not every boundary fragment grown into a
    duplicate. Built over a -1..N+1 supercell, grouped into molecules, then centroid-filtered."""
    nx, ny, nz = cells or getattr(cfg, "pack_cells", (1, 1, 1))
    base = crystal_engine.expand(struct)
    if getattr(cfg, "cell_fill", "molecule") == "clip":   # literal cell contents, molecules cut at the outer box
        out = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    t = np.array([i, j, k], float)
                    for u in base:
                        fr = u["frac"] + t
                        out.append({"sym": u["sym"], "label": u["label"], "occ": u["occ"],
                                    "frac": fr, "xyz": struct.cart(fr)})
        return out
    sup = []
    for i in range(-1, nx + 1):
        for j in range(-1, ny + 1):
            for k in range(-1, nz + 1):
                t = np.array([i, j, k], float)
                for u in base:
                    fr = u["frac"] + t
                    sup.append({"sym": u["sym"], "label": u["label"], "occ": u["occ"],
                                "frac": fr, "xyz": struct.cart(fr)})
    out = []
    for comp in _components(sup, crystal_engine.disorder_alternatives(struct)):
        cen = np.mean([a["frac"] for a in comp], axis=0)
        if (-1e-4 <= cen[0] < nx) and (-1e-4 <= cen[1] < ny) and (-1e-4 <= cen[2] < nz):
            out.extend(comp)
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
    alt = crystal_engine.disorder_alternatives(struct)
    donors = set(cfg.hbond_donors) | ({"C"} if cfg.hbond_weak else set())
    acc = set(cfg.hbond_acceptors)
    floor = cfg.hbond_angle_min
    cell_atoms = [{**u, "xyz": struct.cart(u["frac"])} for u in crystal_engine.expand(struct)]
    origin = struct.cart(np.zeros(3))
    sup = []
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            for dk in (-1, 0, 1):
                shift = struct.cart(np.array([di, dj, dk], float)) - origin
                for a in cell_atoms:
                    sup.append({"sym": a["sym"], "label": a["label"], "occ": a["occ"],
                                "frac": a["frac"] + np.array([di, dj, dk], float), "xyz": a["xyz"] + shift})

    def key(p):
        return (int(round(p[0] * 50)), int(round(p[1] * 50)), int(round(p[2] * 50)))

    central_keys = {key(a["xyz"]) for a in central}

    def hbond_ok(D, h, A):
        dDH = float(np.linalg.norm(D["xyz"] - h["xyz"]))
        if dDH > crystal_engine._covsum(D["sym"], "H") + crystal_engine.BOND_TOL or D["sym"] not in donors:
            return False
        hpos = h["xyz"]
        if cfg.xh_normalize and D["sym"] in crystal_engine.NEUTRON_XH and dDH > 0:
            hpos = D["xyz"] + (h["xyz"] - D["xyz"]) / dDH * crystal_engine.NEUTRON_XH[D["sym"]]
        DA = float(np.linalg.norm(A["xyz"] - D["xyz"]))
        if not (0.4 < DA <= 4.0):
            return False
        if float(np.linalg.norm(A["xyz"] - hpos)) > crystal_engine._vdw("H") + crystal_engine._vdw(A["sym"]):
            return False
        return crystal_engine._angle(D["xyz"], hpos, A["xyz"]) >= floor

    cheavy = [a for a in central if a["sym"] != "H"]
    seeds = []
    for h in [a for a in central if a["sym"] == "H"]:               # central donor -> outside acceptor
        if not cheavy:
            break
        D = min(cheavy, key=lambda a: np.linalg.norm(a["xyz"] - h["xyz"]))
        for A in sup:
            if A["sym"] in acc and key(A["xyz"]) not in central_keys and hbond_ok(D, h, A):
                seeds.append(A)
    supheavy = [a for a in sup if a["sym"] != "H"]
    cacc = [a for a in central if a["sym"] in acc]
    for h in [a for a in sup if a["sym"] == "H" and key(a["xyz"]) not in central_keys]:  # outside donor -> central acceptor
        if not supheavy:
            break
        D = min(supheavy, key=lambda a: np.linalg.norm(a["xyz"] - h["xyz"]))
        if any(hbond_ok(D, h, A) for A in cacc):
            seeds.append(D)

    # ---- assemble neighbour atoms at the requested extent
    mode = getattr(cfg, "hbond_neighbour", "stub")
    seen_seed, seed_list = set(central_keys), []
    for s in seeds:                                   # unique contact atoms (the seeds)
        kk = key(s["xyz"])
        if kk not in seen_seed:
            seen_seed.add(kk); seed_list.append(s)
    have, neigh = set(central_keys), []
    for s in seed_list:
        have.add(key(s["xyz"])); neigh.append(s)
    if mode in ("stub", "whole"):                     # grow outward (1 shell for stub, fully for whole)
        depth = {id(s): 0 for s in seed_list}
        frontier = list(seed_list)
        while frontier:
            nf = []
            for a in frontier:
                if mode == "stub" and depth[id(a)] >= 1:
                    continue
                for s in sup:
                    kk = key(s["xyz"])
                    if kk in have or not (0.4 < float(np.linalg.norm(a["xyz"] - s["xyz"]))
                                          <= crystal_engine._covsum(a["sym"], s["sym"]) + crystal_engine.BOND_TOL):
                        continue
                    if a["occ"] < 1 and s["occ"] < 1 and \
                            float(np.linalg.norm(a["xyz"] - s["xyz"])) < crystal_engine.DISORDER_MIN:
                        continue
                    if s["label"] in alt.get(a["label"], set()):   # disorder-alternative -> not bonded
                        continue
                    have.add(kk); neigh.append(s); nf.append(s); depth[id(s)] = depth[id(a)] + 1
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
