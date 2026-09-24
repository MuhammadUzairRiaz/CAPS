"""Visualisation & analysis pipeline boards (OVITO-style workflow, CAPS design)."""
import math
from lib import *
from screen_app import card
from screen_studio2 import bar, spacer
from screen_row16b import gasteiger
from display import ps_cell, style_bar

CHAIN = ["#F0A83C", "#6CC4D8", "#DE775D", "#9B7AD5", "#7DC884", "#D6AC5C", "#E9ECEF", "#2271DB", "#C77DBA", "#8FB8A8"]
VIRIDIS = ["#440154", "#3B528B", "#21908C", "#5DC963", "#FDE725"]
DIVERGE = ["#2271DB", "#E9ECEF", "#E35049"]
GAFF = {("C", True): "ca", ("C", False): "c3", ("H", True): "ha", ("H", False): "hc"}
LTYPE = {"c3": 1, "ca": 2, "hc": 3, "ha": 4}
TYPECOL = {"c3": "#8B969E", "ca": "#D6AC5C", "hc": "#E9ECEF", "ha": "#6CC4D8"}


def ramp(stops, t):
    t = min(1.0, max(0.0, t)) * (len(stops) - 1)
    i = min(int(t), len(stops) - 2)
    return mix(stops[i], stops[i + 1], t - i)


def annotate(cell):
    """Per-atom records: id, mol, GAFF type, LAMMPS type, unit, charge, xyz — all from the built cell."""
    recs, aid = [], 0
    for m, (atoms, bonds, backbone) in enumerate(cell):
        q, _ = gasteiger(atoms, bonds)
        nb = [[] for _ in atoms]
        for a, b in bonds:
            nb[a].append(b); nb[b].append(a)
        nbk = len(backbone)
        unit = [0] * len(atoms)
        x = nbk
        for k in range(nbk):
            unit[k] = k // 2
            n_ = 12 if k % 2 == 1 else 2
            for j in range(x, x + n_):
                unit[j] = k // 2
            x += n_
        unit[x] = 0; unit[x + 1] = (nbk - 1) // 2
        for j, a in enumerate(atoms):
            if a["e"] == "H":
                arom = atoms[nb[j][0]].get("ar", False)
            else:
                arom = a.get("ar", False)
            g = GAFF[(a["e"], bool(arom))]
            aid += 1
            recs.append({"id": aid, "mol": m + 1, "gaff": g, "type": LTYPE[g], "unit": unit[j] + 1, "q": q[j], "p": a["p"], "e": a["e"], "chain": m, "j": j})
    return recs


def draw(sc, cell, recs, colour_of, h_scale=None):
    it = iter(recs)
    for atoms, bonds, backbone in cell:
        sel = []
        for a in atoms:
            r = next(it)
            sel.append(dict(a, c=colour_of(r)))
        sc.add_atoms(sel, bonds)


def _shell(title, crumb, content, status):
    body = (topbar(["PS_melt.lammpstrj"], "PS_melt.lammpstrj", crumb) + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Analyze")
            + f'<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">{content}</div></div>' + statusbar(*status))
    return page(title, body)


# ---------------------------------------------------------------- Pipeline
def pipeline():
    L = 33.0
    cell = ps_cell(10, 8, L, seed=21)
    recs = annotate(cell)
    n = len(recs)
    lp_w, rp_w = 300, 330
    cw = 1368 - lp_w - rp_w
    vh = 900 - 44 - 26 - 40 - 200 - 48 - 12
    sc = Scene("vp", cw - 6, vh, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.25)
    zcut = (L / 2 - 6, L / 2 + 6)
    kept = [r for r in recs]
    draw(sc, cell, recs, lambda r: CHAIN[r["chain"] % len(CHAIN)] if r["e"] != "H" else mix(CHAIN[r["chain"] % len(CHAIN)], "#FFFFFF", 0.55))
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    view = sc.svg()
    mods = [("Python script", "compute end-to-end vectors", False, "terminal"), ("Radial distribution", "C–C, 10 Å, 200 bins", True, "chart"),
            ("Cluster analysis", "by molecule · 10 clusters", True, "layers"), ("Expression selection", "Type == 2 && Position.Z > 13", False, "filter"),
            ("Slice", "normal (0,0,1) · width 12 Å", False, "scissors"), ("Colour coding", "Molecule Identifier · categorical", True, "eye")]
    items = []
    for i, (nm, sub_, on, ic) in enumerate(mods):
        sel = nm == "Colour coding"
        items.append(f'<div style="display: flex; align-items: center; gap: 8px; padding: 7px 10px; border-radius: 6px; background: {BG3 if sel else "transparent"}; border: 1px solid {ACC if sel else "transparent"}">'
                     f'<span style="width: 15px; height: 15px; border-radius: 3px; border: 1.5px solid {ACC if on else DIM}; background: {ACC if on else "transparent"}; flex-shrink: 0"></span>'
                     f'{icon(ic, 14, MUTED)}<span style="display: flex; flex-direction: column; min-width: 0"><span style="font-size: 12.5px; color: {TEXT if on else MUTED}">{nm}</span>'
                     f'<span style="font-size: 11px; color: {DIM}; white-space: nowrap; overflow: hidden; text-overflow: ellipsis">{sub_}</span></span></div>')
    src = (f'<div style="display: flex; flex-direction: column; gap: 3px; padding: 8px 10px; border-radius: 6px; background: {BG0}; border: 1px solid {LINE}">'
           f'<span style="font-size: 11px; color: {DIM}; letter-spacing: 0.06em; text-transform: uppercase">Data source</span><span style="font-size: 12.5px">LAMMPS dump · PS_melt.lammpstrj</span>'
           f'<span style="font-family: {MONO}; font-size: 11px; color: {MUTED}">{n:,} particles · frame 120 / 500 · + topology from PS_melt.data</span></div>'.replace(",", " "))
    vis = "".join(f'<div style="display: flex; align-items: center; gap: 8px; padding: 5px 10px; font-size: 12.5px">{icon("eye", 14, MUTED)}{v}</div>' for v in ("Particles", "Bonds", "Simulation cell", "Colour legend"))
    left = (f'<aside style="width: {lp_w}px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-right: 1px solid {LINE}">'
            + panel_head("Pipeline", row(tbtn("plus", "Add step"), tbtn("dots", "More"), gap=2))
            + f'<div style="padding: 8px; display: flex; flex-direction: column; gap: 3px"><span style="padding: 4px 10px; font-size: 11px; color: {DIM}; letter-spacing: 0.06em; text-transform: uppercase">Steps · run bottom to top</span>{"".join(items)}</div>'
            + f'<div style="padding: 4px 8px 8px">{src}</div>'
            + f'<div style="padding: 0 8px 8px; display: flex; flex-direction: column"><span style="padding: 4px 10px; font-size: 11px; color: {DIM}; letter-spacing: 0.06em; text-transform: uppercase">Visual elements</span>{vis}</div>'
            + '</aside>')
    legend = "".join(f'<div style="display: flex; align-items: center; gap: 6px; font-size: 11.5px"><span style="width: 12px; height: 12px; border-radius: 3px; background: {CHAIN[i]}"></span>mol {i + 1}</div>' for i in range(len(cell)))
    right = (f'<aside style="width: {rp_w}px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">'
             + panel_head("Colour coding", chip("step 1", MUTED, BG2))
             + section("Input", col(select("Property", "Molecule Identifier"), select("Operate on", "Particles"), seg(["Categorical", "Continuous"], "Categorical", full=True), gap=8))
             + section("Colours", col(select("Palette", "CAPS chains · colour-vision safe"), toggle("Lighten hydrogens", True), toggle("Only selected", False), gap=8))
             + section("Legend", f'<div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 6px">{legend}</div>')
             + '</aside>')
    rows_ = recs[:5] + recs[-1:]
    dt = table(["id", "mol", "type", "name", "q (e)", "x", "y", "z"], [[str(r["id"]), str(r["mol"]), str(r["type"]), r["gaff"], f"{r['q']:+.4f}", f"{r['p'][0]:.3f}", f"{r['p'][1]:.3f}", f"{r['p'][2]:.3f}"] for r in rows_],
               ["9%", "8%", "8%", "10%", "14%", "17%", "17%", "17%"], mono_cols=tuple(range(8)), align_right=(0, 1, 2, 4, 5, 6, 7), fs=11, rowh=18)
    data = (f'<div style="height: 200px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-top: 1px solid {LINE}">'
            f'{tabs(["Particles", "Bonds", "Global attributes", "Data tables"], "Particles", 12)}'
            f'<div style="padding: 4px 12px; overflow: hidden">{dt}</div></div>')
    timeline = (f'<div style="height: 40px; flex-shrink: 0; display: flex; align-items: center; gap: 10px; padding: 0 12px; background: {BG1}; border-top: 1px solid {LINE}">'
                f'{tbtn("play", "Play")}<span style="font-family: {MONO}; font-size: 11.5px; color: {MUTED}">frame 120</span>'
                f'<div style="flex-grow: 1; height: 6px; border-radius: 3px; background: {BG3}; position: relative"><div style="width: 24%; height: 6px; border-radius: 3px; background: {ACC}"></div></div>'
                f'<span style="font-family: {MONO}; font-size: 11.5px; color: {DIM}">500 frames · 1 ps each</span></div>')
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
              + bar(tbtn("rotate", "Orbit", True), tbtn("move", "Pan"), tbtn("search", "Zoom"), sep(), tbtn("cursor", "Pick"), tbtn("ruler", "Measure"), spacer(), chip("Perspective", MUTED, BG2))
              + f'<div style="position: relative; flex-grow: 1; background: {BG0}; overflow: hidden">{view}{style_bar("All atoms", "Colour: molecule · H lightened", "right: 12px; top: 10px")}</div>'
              + timeline + data + '</div>')
    content = f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{center}{right}</div>'
    return _shell("CAPS — visualize pipeline", f'{icon("chevr", 12, DIM)}<span>Analyze › Visualize</span>', content,
                  (f"<span>{n:,} particles · {len(cell)} molecules · 4 types</span>".replace(",", " "), "<span>pipeline cached · 6 steps</span>"))


# ---------------------------------------------------------------- Colour coding gallery
def colour_by():
    L = 25.0
    cell = ps_cell(6, 6, L, seed=8)
    recs = annotate(cell)
    zs = [r["p"][2] for r in recs]
    zmin, zmax = min(zs), max(zs)
    qmax = max(abs(r["q"]) for r in recs)
    umax = max(r["unit"] for r in recs)
    specs = [("Element", "C grey, H white", lambda r: ELEM[r["e"]][0], [("C", ELEM["C"][0]), ("H", ELEM["H"][0])]),
             ("Molecule (chain)", "Molecule Identifier", lambda r: CHAIN[r["chain"] % 10] if r["e"] != "H" else mix(CHAIN[r["chain"] % 10], "#FFFFFF", 0.55), [(f"mol {i + 1}", CHAIN[i]) for i in range(len(cell))]),
             ("Atom type", "force-field type (GAFF2)", lambda r: TYPECOL[r["gaff"]], [(k, v) for k, v in TYPECOL.items()]),
             ("Monomer index", f"unit 1 → {umax} along each chain", lambda r: ramp(VIRIDIS, (r["unit"] - 1) / max(1, umax - 1)), None),
             ("Position z", f"{zmin:.1f} → {zmax:.1f} Å", lambda r: ramp(VIRIDIS, (r["p"][2] - zmin) / (zmax - zmin)), None),
             ("Partial charge", f"Gasteiger · ±{qmax:.3f} e", lambda r: ramp(DIVERGE, 0.5 + r["q"] / (2 * qmax)), None)]
    pw, ph = 432, 250
    panels = []
    for i, (title, sub_, fn, cats) in enumerate(specs):
        sc = Scene(f"cb{i}", pw, ph, yaw=0.55, pitch=0.4, persp=0.25, fog=0.25, atom_k=1.3)
        draw(sc, cell, recs, fn)
        sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 0.8, "4 4", 0.5)
        sc._prep(); sc.scale *= 1.2
        if cats:
            leg = "".join(f'<span style="display: inline-flex; align-items: center; gap: 5px; font-size: 11px; color: {MUTED}"><span style="width: 10px; height: 10px; border-radius: 2px; background: {c}"></span>{k}</span>' for k, c in cats)
        else:
            stops = DIVERGE if "charge" in title else VIRIDIS
            lo, hi = (f"−{qmax:.3f}", f"+{qmax:.3f}") if "charge" in title else ((f"{zmin:.1f}", f"{zmax:.1f} Å") if "z" in title else ("1", str(umax)))
            leg = (f'<span style="font-family: {MONO}; font-size: 10.5px; color: {DIM}">{lo}</span><span style="flex-grow: 1; height: 8px; border-radius: 4px; background: linear-gradient(90deg, {", ".join(stops)})"></span>'
                   f'<span style="font-family: {MONO}; font-size: 10.5px; color: {DIM}">{hi}</span>')
        panels.append(f'<div style="display: flex; flex-direction: column; gap: 6px"><div style="display: flex; align-items: baseline; gap: 8px"><span style="font-size: 13px; font-weight: 600">{title}</span><span style="font-size: 11.5px; color: {DIM}">{sub_}</span></div>'
                      f'<div style="background: {BG0}; border: 1px solid {ACC if i == 1 else LINE}; border-radius: 8px; overflow: hidden">{sc.svg()}</div>'
                      f'<div style="display: flex; align-items: center; gap: 8px; flex-wrap: wrap; min-height: 16px">{leg}</div></div>')
    grid = f'<div style="display: grid; grid-template-columns: repeat(3, {pw}px); gap: 14px 16px">{"".join(panels)}</div>'
    n = len(recs)
    content = (f'<div style="flex-shrink: 0; display: flex; align-items: center; gap: 16px; padding: 16px 22px; background: {BG1}; border-bottom: 1px solid {LINE}">'
               f'<div style="display: flex; flex-direction: column; gap: 3px"><h1 style="margin: 0; font-size: 19px; font-weight: 600">Colour by</h1><div style="font-size: 12.5px; color: {MUTED}">One all-atom cell, six colourings · any per-atom property can drive the colour</div></div><div style="flex-grow: 1"></div>'
               f'{select("Property", "Molecule Identifier", 220)}{select("Map", "categorical · CAPS chains", 220)}</div>'
               f'<div style="flex-grow: 1; padding: 14px 22px; overflow: hidden">{grid}</div>')
    return _shell("CAPS — colour by property", f'{icon("chevr", 12, DIM)}<span>Analyze › Visualize › Colour by</span>', content,
                  (f"<span>{len(cell)} chains · {n} atoms · every value computed from this cell</span>", "<span>continuous maps: viridis, diverging blue–red</span>"))


# ---------------------------------------------------------------- Import: LAMMPS and other formats
def import_formats():
    L = 33.0
    cell = ps_cell(10, 8, L, seed=21)
    recs = annotate(cell)
    n = len(recs)
    lines = ["ITEM: TIMESTEP", "120000", "ITEM: NUMBER OF ATOMS", str(n), "ITEM: BOX BOUNDS pp pp pp", f"0.0000000000000000e+00 {L:.16e}", f"0.0000000000000000e+00 {L:.16e}", f"0.0000000000000000e+00 {L:.16e}",
             "ITEM: ATOMS id mol type q xu yu zu"]
    lines += [f"{r['id']} {r['mol']} {r['type']} {r['q']:.5f} {r['p'][0]:.4f} {r['p'][1]:.4f} {r['p'][2]:.4f}" for r in recs[:5]] + ["…"]
    pv = "".join(f'<div style="display: flex; gap: 12px"><span style="width: 24px; text-align: right; color: {DIM}">{i + 1}</span><span style="white-space: pre; color: {ACC if l.startswith("ITEM") else TEXT}">{esc(l)}</span></div>' for i, l in enumerate(lines))
    cols = [("id", "Particle Identifier", "int"), ("mol", "Molecule Identifier", "int"), ("type", "Particle Type", "int"), ("q", "Charge", "float"),
            ("xu", "Position.X (unwrapped)", "float"), ("yu", "Position.Y (unwrapped)", "float"), ("zu", "Position.Z (unwrapped)", "float")]
    ct = table(["Column", "Maps to", "Kind"], [[a, b, c] for a, b, c in cols], ["22%", "56%", "22%"], mono_cols=(0, 2), fs=12, rowh=26)
    tn = table(["Type", "Name", "Mass", "Element"], [[str(v), k, f"{12.011 if k[0] == 'c' else 1.008}", "C" if k[0] == "c" else "H"] for k, v in LTYPE.items()], ["16%", "28%", "28%", "28%"], mono_cols=(0, 1, 2), align_right=(2,), fs=12, rowh=26)
    formats = [("LAMMPS data", "atom_style full · molecular · charge · atomic; Masses, Bonds, Angles…"), ("LAMMPS dump", "text custom columns · gzip · multi-file wildcard · binary"),
               ("LAMMPS restart", "via LAMMPS library, when installed"), ("XYZ / extended XYZ", "per-frame lattice and properties"), ("PDB · mmCIF", "with CONECT and residues"),
               ("GROMACS", ".gro · .top · .xtc · .trr · .tpr"), ("DCD", "CHARMM/NAMD trajectories"), ("CIF · POSCAR · CONTCAR", "crystals, symmetry expanded"),
               ("mol2 · SDF", "bond orders and charges"), ("CAPS .caps", "full document with provenance")]
    fl = "".join(f'<div style="display: flex; gap: 10px; padding: 4px 0; border-bottom: 1px solid {BG2}; font-size: 11.5px"><span style="width: 150px; flex-shrink: 0; font-weight: 600">{a}</span><span style="color: {MUTED}">{b}</span></div>' for a, b in formats)
    left = (f'<div style="width: 560px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("PS_melt.lammpstrj · first lines", f'<div style="font-family: {MONO}; font-size: 11.5px; line-height: 1.6; overflow: hidden">{pv}</div>', chip("detected: LAMMPS dump", OK, BG2), 14)
            + card("Formats CAPS reads", f'<div>{fl}</div>', chip("write: same list", MUTED, BG2), 14) + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Column mapping", col(ct, f'<span style="font-size: 11.5px; color: {DIM}">Unknown columns (c_pe, v_temp…) import as custom per-atom properties and can drive colour, selection and analysis.</span>', gap=8), chip("7 of 7 mapped", OK, BG2), 14)
             + card("Types, masses, topology", col(tn, row(select("Bonds from", "PS_melt.data (topology file)"), select("Frames", "all 500 · load on demand"), gap=8), gap=10), "", 14)
             + '</div>')
    content = (f'<div style="flex-shrink: 0; display: flex; align-items: center; gap: 16px; padding: 16px 22px; background: {BG1}; border-bottom: 1px solid {LINE}">'
               f'<div style="display: flex; flex-direction: column; gap: 3px"><h1 style="margin: 0; font-size: 19px; font-weight: 600">Open · LAMMPS dump</h1><div style="font-size: 12.5px; color: {MUTED}">Format detected from content, columns mapped to named properties, topology joined from the data file</div></div>'
               f'<div style="flex-grow: 1"></div>{btn("Cancel")}{btn("Open 500 frames", True, "folder")}</div>'
               f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 22px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return _shell("CAPS — open LAMMPS and other formats", f'{icon("chevr", 12, DIM)}<span>Open file</span>', content,
                  (f"<span>preview rows are the real values of this cell · {n} atoms</span>", "<span>units: LAMMPS real</span>"))


# ---------------------------------------------------------------- Analysis library
def modifier_library():
    cats = [("Colour & style", [("Colour coding", "any property, categorical or continuous"), ("Assign colour", "to selection"), ("Colour by type / molecule", "one click"), ("Transparency", "by property"), ("Particle radius", "by type or property")]),
            ("Select", [("Expression selection", "Type == 2 && Position.Z > 13"), ("By type · molecule · chain", ""), ("Manual pick & lasso", ""), ("Expand selection", "by bonds or distance"), ("Invert · clear", "")]),
            ("Modify", [("Slice", "slab by normal and width"), ("Delete selected", ""), ("Replicate", "periodic images"), ("Affine transform", "shear, strain, rotate"), ("Wrap / unwrap", "image flags"), ("Compute property", "expressions per atom"), ("Combine datasets", "add a second file"), ("Freeze property", "keep frame-0 values"), ("Smooth trajectory", "time averaging")]),
            ("Structure", [("Common neighbour analysis", "fcc, bcc, hcp, icosahedral"), ("Polyhedral template matching", "orientation + strain"), ("Centrosymmetry", ""), ("Voronoi analysis", "volumes, faces"), ("Chain orientation", "P₂, Herman's f"), ("Crystallinity (polymer)", "local chord alignment")]),
            ("Measure", [("Coordination & RDF", "partial g(r) by type"), ("Cluster analysis", "by cutoff or bonds, sizes, Rg"), ("Displacements & MSD", "vs reference frame"), ("Wigner–Seitz defects", "vacancies, interstitials"), ("Spatial binning", "1-D/2-D/3-D profiles"), ("Histogram · scatter · time series", ""), ("Surface mesh", "free volume, pores")]),
            ("Visual", [("Create bonds", "cutoff or topology"), ("Trajectory lines", "particle paths"), ("Vectors", "forces, velocities, dipoles"), ("Polyhedra", "coordination shells"), ("Colour legend · text labels", "on renders")]),
            ("Automate", [("Python step", "caps.pipeline API in the same process"), ("Batch over files", "same pipeline, many trajectories"), ("Save pipeline", "YAML, reusable, versioned"), ("Render movie", "frames → MP4 / PNG series")])]
    cardsh = []
    for name, items in cats:
        li = "".join(f'<div style="display: flex; flex-direction: column; gap: 0; padding: 4px 0; border-bottom: 1px solid {BG2}"><span style="font-size: 12px">{a}</span>' + (f'<span style="font-size: 10.5px; color: {DIM}">{b}</span>' if b else "") + '</div>' for a, b in items)
        cardsh.append(card(name, f'<div>{li}</div>', chip(str(len(items)), MUTED, BG2, True), 12))
    total = sum(len(i) for _, i in cats)
    grid = f'<div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; align-items: start">{"".join(cardsh)}</div>'
    content = (f'<div style="flex-shrink: 0; display: flex; align-items: center; gap: 16px; padding: 10px 22px; background: {BG1}; border-bottom: 1px solid {LINE}">'
               f'<div style="display: flex; flex-direction: column; gap: 3px"><h1 style="margin: 0; font-size: 19px; font-weight: 600">Add pipeline step</h1><div style="font-size: 12.5px; color: {MUTED}">{total} steps in 7 groups · each is a command, so it also runs from the CLI and Python</div></div>'
               f'<div style="flex-grow: 1"></div><div style="width: 320px">{field("Search steps", "rdf", mono=False)}</div></div>'
               f'<div style="flex-grow: 1; padding: 10px 22px; overflow: hidden">{grid}</div>')
    return _shell("CAPS — pipeline steps", f'{icon("chevr", 12, DIM)}<span>Analyze › Visualize › Add step</span>', content,
                  ("<span>steps are non-destructive; the source file is never changed</span>", "<span>results cached per frame</span>"))


# ---------------------------------------------------------------- Import: GROMACS (+ format matrix)
def import_gromacs():
    L = 33.0
    cell = ps_cell(10, 8, L, seed=21)
    recs = annotate(cell)
    n = len(recs)
    counters = {}
    gro = ["PS melt, 10 chains, t= 120.00000", f"{n:5d}"]
    for r in recs[:10]:
        k = (r["mol"], r["e"]); counters[k] = counters.get(k, 0) + 1
        name = f"{r['e']}{counters[k]}"
        gro.append(f"{r['unit']:5d}{'STY':<5}{name:>5}{r['id']:5d}{r['p'][0] / 10:8.3f}{r['p'][1] / 10:8.3f}{r['p'][2] / 10:8.3f}")
    gro.append("  …")
    gro.append(f"{L / 10:10.5f}{L / 10:10.5f}{L / 10:10.5f}")
    pv = "".join(f'<div style="display: flex; gap: 12px"><span style="width: 24px; text-align: right; color: {DIM}">{i + 1}</span><span style="white-space: pre; color: {ACC if i in (0, 1) or i == len(gro) - 1 else TEXT}">{esc(l)}</span></div>' for i, l in enumerate(gro))
    top = ["[ molecules ]", "; name      count", f"PS_DP8       {len(cell)}"]
    tv = "".join(f'<div style="white-space: pre; color: {ACC if l.startswith("[") else (DIM if l.startswith(";") else TEXT)}">{esc(l)}</div>' for l in top)
    files = [("PS_melt.gro", "coordinates · box · frame 0", True), ("PS_melt.top + ps.itp", "bonds, types, charges", True), ("traj.xtc", "500 frames · compressed, 0.001 nm precision", True), ("PS_melt.tpr", "optional · run input", False)]
    fl = "".join(f'<div style="display: flex; align-items: center; gap: 10px; padding: 8px 10px; border-radius: 7px; background: {BG0}; border: 1px solid {LINE}">{icon("check", 15, OK) if on else icon("file", 15, DIM)}'
                 f'<span style="display: flex; flex-direction: column"><span style="font-family: {MONO}; font-size: 12px">{a}</span><span style="font-size: 11px; color: {DIM}">{b}</span></span></div>' for a, b, on in files)
    Y, N, P = f'<span style="color: {OK}">●</span>', f'<span style="color: {DIM}">—</span>', f'<span style="color: {ACC}">◐</span>'
    fm = [("LAMMPS data", Y, Y, Y, N, Y), ("LAMMPS dump", Y, Y, P, Y, Y), ("GROMACS .gro", Y, Y, N, Y, Y), ("GROMACS .top/.itp", Y, Y, Y, N, Y), ("GROMACS .xtc/.trr", Y, N, N, Y, Y),
          ("PDB · mmCIF", Y, Y, P, Y, Y), ("XYZ / extended XYZ", Y, P, N, Y, Y), ("DCD", Y, N, N, Y, Y), ("CIF · POSCAR", Y, Y, N, N, Y), ("mol2 · SDF", Y, P, Y, N, Y), ("AMBER prmtop · nc", Y, Y, Y, Y, P)]
    mt = table(["Format", "Coords", "Cell", "Bonds & charges", "Frames", "Write"], [list(r) for r in fm], ["32%", "12%", "10%", "20%", "12%", "14%"], fs=11.5, rowh=25)
    left = (f'<div style="width: 560px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("PS_melt.gro · first lines", f'<div style="font-family: {MONO}; font-size: 11.5px; line-height: 1.6; overflow: hidden">{pv}</div>', chip("fixed-width .gro · nm", OK, BG2), 14)
            + card("PS_melt.top", f'<div style="font-family: {MONO}; font-size: 11.5px; line-height: 1.7">{tv}</div>', chip("includes ps.itp", MUTED, BG2), 14) + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Files in this set", col(fl, f'<span style="font-size: 11.5px; color: {DIM}">Dropping any one file finds the others by name. Units convert on read: nm → Å (× 10), ps stays ps, kJ/mol shown in your unit setting.</span>', gap=6), chip("3 of 4 found", OK, BG2), 14)
             + card("What each format carries", col(mt, f'<div style="font-size: 11px; color: {DIM}">{Y} full · {P} partial · {N} not in the format. Missing bonds are perceived and flagged; missing charges come from Field.</div>', gap=6), "", 14)
             + '</div>')
    content = (f'<div style="flex-shrink: 0; display: flex; align-items: center; gap: 16px; padding: 16px 22px; background: {BG1}; border-bottom: 1px solid {LINE}">'
               f'<div style="display: flex; flex-direction: column; gap: 3px"><h1 style="margin: 0; font-size: 19px; font-weight: 600">Open · GROMACS set</h1><div style="font-size: 12.5px; color: {MUTED}">Coordinates, topology and trajectory opened together as one document</div></div>'
               f'<div style="flex-grow: 1"></div>{btn("Cancel")}{btn("Open set", True, "folder")}</div>'
               f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 22px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return _shell("CAPS — open GROMACS files", f'{icon("chevr", 12, DIM)}<span>Open file</span>', content,
                  (f"<span>preview rows are this cell in nm · {n} atoms</span>", "<span>GROMACS units: nm, ps, kJ/mol</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage16/project", exist_ok=True)
    for name, fn in (("VisPipeline", pipeline), ("ColourBy", colour_by), ("OpenLammps", import_formats), ("OpenGromacs", import_gromacs), ("PipelineSteps", modifier_library)):
        h = fn(); open(f"stage16/project/{name}.dc.html", "w").write(h); print(name, len(h))
