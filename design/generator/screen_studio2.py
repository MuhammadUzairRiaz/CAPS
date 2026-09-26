import math, random
from lib import *
from mols import *
from screen_app import shell, pagehead, card, progress
from screen_builders import VPW, VPH, footer, hud

KCAL_PER_K = 0.0019872036  # R in kcal/(mol K)


def studio_shell(doc, title, toolbar, center, right, status, right_w=380):
    pn = f'<aside aria-label="Panel" style="width: {right_w}px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{right}</aside>'
    body = (topbar([doc], doc) + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">' + toolbar
            + f'<div style="flex-grow: 1; display: flex; min-height: 0">{center}{pn}</div></div></div>' + statusbar(*status))
    return page(title, body)


def bar(*items):
    return (f'<div role="toolbar" style="height: 48px; flex-shrink: 0; display: flex; align-items: center; gap: 2px; padding: 0 10px; background: {BG1}; border-bottom: 1px solid {LINE}">'
            + "".join(items) + '</div>')


def spacer():
    return '<div style="flex-grow: 1"></div>'


def styrene3d():
    """Planar styrene: aromatic C-C 1.39, C(ar)-C 1.47, C=C 1.34, C-H 1.08 Angstrom, 120 degree angles."""
    def pol(o, r, deg):
        return (o[0] + r * math.cos(math.radians(deg)), o[1] + r * math.sin(math.radians(deg)), 0.0)
    atoms = [{"e": "C", "p": pol((0, 0, 0), 1.39, 90 + 60 * k)} for k in range(6)]
    atoms += [{"e": "H", "p": pol((0, 0, 0), 2.47, 90 + 60 * k)} for k in range(1, 6)]
    ca = pol(atoms[0]["p"], 1.47, 90)
    cb = pol(ca, 1.34, 30)
    atoms += [{"e": "C", "p": ca}, {"e": "C", "p": cb}, {"e": "H", "p": pol(ca, 1.08, 150)},
              {"e": "H", "p": pol(cb, 1.08, 90)}, {"e": "H", "p": pol(cb, 1.08, -30)}]
    return atoms


# ------------------------------------------------------------------ Start screen
def welcome():
    def thumb(sid, builder):
        sc = Scene(sid, 250, 140, yaw=0.5, pitch=0.4, persp=0.2, fog=0.5, outline=False, bond_w=0.22)
        builder(sc)
        return sc.svg()

    def t_ps(sc):
        a, _, _ = polystyrene(4, "atactic", seed=3, hydrogens=False, curve=False)
        sc.add_atoms(a)

    def t_pe(sc):
        sc.add_atoms(pe_crystal(2, 2, 3))

    def t_cell(sc):
        cols = [ACC, SEL, "#E07A5F", "#9B7BD6", OK]
        for k in range(5):
            sc.add_tube(random_chain((10 + k * 3, 10 + (k % 2) * 8, 10 + k * 2), 30, 300 + k, box=((0, 0, 0), (25, 25, 25)), persistence=0.6), cols[k], 0.6)
        sc.add_box((0, 0, 0), (25, 0, 0), (0, 25, 0), (0, 0, 25), MUTED, 1, None, 0.8)

    recents = [("PS_atactic_DP40.caps", "PS-tacticity study", "edited 2 h ago", t_ps, "ps"),
               ("PE_crystal.caps", "Crystallinity", "edited yesterday", t_pe, "pe"),
               ("PS_cell_20x40.caps", "PS-tacticity study", "grown Sep 21", t_cell, "cl")]
    rc = "".join(f'<a href="Main.dc.html" style="display: flex; flex-direction: column; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px; overflow: hidden; text-decoration: none; color: {TEXT}">'
                 f'<div style="background: {BG0}; border-bottom: 1px solid {LINE}">{thumb(s, b)}</div>'
                 f'<div style="padding: 10px 12px; display: flex; flex-direction: column; gap: 3px"><span style="font-size: 13px; font-weight: 500">{n}</span><span style="font-size: 11.5px; color: {MUTED}">{p} · {w}</span></div></a>' for n, p, w, b, s in recents)
    starts = [("hex", "Molecule", "From SMILES, InChI, a name or a 2D sketch"), ("grow", "Polymer", "Repeat units, sequence, tacticity"),
              ("cube", "Crystal", "Space group, lattice, CIF import"), ("pack", "Amorphous cell", "Grow and pack a periodic cell"),
              ("layers", "Surface or interface", "Cleave, stack, add vacuum"), ("flask", "Solvated system", "Box, solvent model, ions")]
    tile_link = {"Molecule": "Sketch.dc.html", "Polymer": "PolymerBuilder.dc.html", "Crystal": "CrystalBuilder.dc.html", "Amorphous cell": "Grow.dc.html",
                 "Surface or interface": "SurfaceBuilder.dc.html", "Solvated system": "SolvationBuilder.dc.html"}
    st = "".join(f'<a href="{tile_link[t]}" style="display: flex; gap: 12px; align-items: flex-start; padding: 14px; text-align: left; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px; text-decoration: none">'
                 f'<span style="width: 38px; height: 38px; border-radius: 8px; background: {BG3}; display: flex; align-items: center; justify-content: center; flex-shrink: 0">{icon(ic, 20, ACC)}</span>'
                 f'<span style="display: flex; flex-direction: column; gap: 3px"><span style="font-size: 13.5px; font-weight: 600; color: {TEXT}">{t}</span><span style="font-size: 12px; color: {MUTED}; line-height: 1.4">{d}</span></span></a>' for ic, t, d in starts)
    learn = [("Build and type a polymer chain", "8 min"), ("Grow an amorphous PS cell", "15 min"), ("Equilibrate with the 21-step protocol", "12 min"),
             ("Script a build in Python", "10 min"), ("Moving from LAMMPS scripts", "5 min")]
    lr = "".join(f'<a href="#" style="display: flex; align-items: center; gap: 10px; height: 36px; padding: 0 4px; border-bottom: 1px solid {BG2}; text-decoration: none; color: {TEXT}; font-size: 12.5px">{icon("play", 13, ACC)}<span style="flex-grow: 1">{t}</span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{m}</span></a>' for t, m in learn)
    quick = (f'<div style="display: flex; align-items: center; gap: 10px; height: 52px; padding: 0 16px; background: {BG1}; border: 1px solid {ACC}; border-radius: 10px">'
             f'{icon("search", 18, ACC)}<label for="qs" style="position: absolute; width: 1px; height: 1px; overflow: hidden">Quick start</label>'
             f'<input id="qs" value="C=Cc1ccccc1" style="flex-grow: 1; background: transparent; border: 0; outline: none; font-family: {MONO}; font-size: 15px; color: {TEXT}">'
             f'{chip("SMILES detected", OK, "#16261A")}{btn("Build 3D", True, "cube", small=True, href="Sketch.dc.html")}</div>')
    drop = (f'<div style="display: flex; align-items: center; justify-content: center; gap: 10px; height: 56px; border: 1.5px dashed {LINE}; border-radius: 10px; color: {MUTED}; font-size: 12.5px">'
            f'{icon("download", 16, DIM)}<span>Drop PDB, CIF, mol2, SDF, XYZ, LAMMPS data or GROMACS gro/top</span></div>')
    main = (f'<div style="flex-grow: 1; display: flex; gap: 28px; padding: 34px 40px; min-height: 0">'
            f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 22px">'
            f'<div style="display: flex; flex-direction: column; gap: 6px"><h1 style="margin: 0; font-size: 26px; font-weight: 600; letter-spacing: -0.015em">Start</h1>'
            f'<p style="margin: 0; color: {MUTED}; font-size: 13.5px">Type a structure, pick a builder, or open recent work.</p></div>'
            + quick
            + f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px">{st}</div>'
            + f'<div style="display: flex; flex-direction: column; gap: 10px"><h2 style="margin: 0; font-size: 13px; font-weight: 600">Recent</h2><div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px">{rc}</div></div>'
            + drop + '</div>'
            + f'<div style="width: 340px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
            + card("Learn", lr, chip("5 tutorials"))
            + card("From LAMMPS", col(f'<div style="font-size: 12.5px; color: {MUTED}; line-height: 1.5">Open a LAMMPS data file with its input script. Atoms, types and force-field coefficients carry over; the run steps become CAPS engine steps.</div>',
                                    btn("Import LAMMPS input", ic="download", small=True), gap=10))
            + card("This machine", col(kv("Threads", "10"), kv("GPU", "Apple M5 · CPU path"), kv("Remote hosts", "none"), gap=7), btn("Settings", ic="gear", small=True, href="Settings.dc.html"))
            + '</div></div>')
    body = (topbar(["Start"], "Start") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + f'<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">{main}</div></div>'
            + statusbar("<span>CAPS 0.1.0</span>", "<span>Deterministic mode available · Philox RNG</span>"))
    return page("CAPS — start", body)


# ------------------------------------------------------------------ 2D sketch → 3D
def sketch():
    cx, cy, R = 330, 330, 62
    pts = [(cx + R * math.cos(math.radians(-90 + 60 * k)), cy + R * math.sin(math.radians(-90 + 60 * k))) for k in range(6)]
    lines = []
    for k in range(6):
        a, b = pts[k], pts[(k + 1) % 6]
        lines.append(f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" stroke="{TEXT}" stroke-width="2.4" stroke-linecap="round"></line>')
        if k % 2 == 1:
            # inner double bond line
            ia = (a[0] + (cx - a[0]) * 0.2, a[1] + (cy - a[1]) * 0.2)
            ib = (b[0] + (cx - b[0]) * 0.2, b[1] + (cy - b[1]) * 0.2)
            sh = 0.12
            ia = (ia[0] + (ib[0] - ia[0]) * sh, ia[1] + (ib[1] - ia[1]) * sh)
            ib = (ib[0] + (ia[0] - ib[0]) * sh, ib[1] + (ia[1] - ib[1]) * sh)
            lines.append(f'<line x1="{ia[0]:.1f}" y1="{ia[1]:.1f}" x2="{ib[0]:.1f}" y2="{ib[1]:.1f}" stroke="{TEXT}" stroke-width="2.4" stroke-linecap="round"></line>')
    ip = pts[0]
    ca = (ip[0] + 62 * math.cos(math.radians(-30)), ip[1] + 62 * math.sin(math.radians(-30)))
    cb = (ca[0], ca[1] - 62)
    lines.append(f'<line x1="{ip[0]:.1f}" y1="{ip[1]:.1f}" x2="{ca[0]:.1f}" y2="{ca[1]:.1f}" stroke="{TEXT}" stroke-width="2.4" stroke-linecap="round"></line>')
    lines.append(f'<line x1="{ca[0]:.1f}" y1="{ca[1]:.1f}" x2="{cb[0]:.1f}" y2="{cb[1]:.1f}" stroke="{TEXT}" stroke-width="2.4" stroke-linecap="round"></line>')
    lines.append(f'<line x1="{ca[0] - 8:.1f}" y1="{ca[1] - 7:.1f}" x2="{cb[0] - 8:.1f}" y2="{cb[1] + 7:.1f}" stroke="{TEXT}" stroke-width="2.4" stroke-linecap="round"></line>')
    # hover highlight + cursor ghost bond
    lines.append(f'<circle cx="{cb[0]:.1f}" cy="{cb[1]:.1f}" r="13" fill="none" stroke="{SEL}" stroke-width="2"></circle>')
    gx, gy = cb[0] + 54, cb[1] - 30
    lines.append(f'<line x1="{cb[0]:.1f}" y1="{cb[1]:.1f}" x2="{gx:.1f}" y2="{gy:.1f}" stroke="{SEL}" stroke-width="2.4" stroke-dasharray="5 4" stroke-linecap="round"></line>')
    lines.append(f'<g transform="translate({gx + 6:.1f},{gy - 4:.1f})">{icon("cursor", 18, TEXT)}</g>')
    grid = "".join(f'<circle cx="{x}" cy="{y}" r="1" fill="{BG3}"></circle>' for x in range(20, 680, 24) for y in range(20, 700, 24))
    s2d = (f'<svg width="660" height="690" viewBox="0 0 660 690" role="img" aria-label="2D sketch of styrene" style="display: block">{grid}{"".join(lines)}'
           f'<text x="{cb[0] + 16:.1f}" y="{cb[1] - 14:.1f}" fill="{SEL}" font-size="11" font-family="IBM Plex Mono">CH2 · drag to grow</text></svg>')
    tools = [("cursor", "Select", False), ("bond", "Single bond", True), ("layers", "Double bond", False), ("ring", "Benzene ring", False), ("hex", "Cyclohexane", False),
             ("plus", "Charge +", False), ("tag", "Stereo wedge", False), ("scissors", "Erase", False)]
    tb = "".join(tbtn(i, l, on) for i, l, on in tools)
    els = "".join(f'<button aria-label="{e}" style="width: 34px; height: 30px; border-radius: 5px; border: 1px solid {ACC if e == "C" else LINE}; background: {BG3 if e == "C" else BG0}; font-family: {MONO}; font-size: 12.5px; color: {TEXT}; cursor: pointer">{e}</button>' for e in ["C", "N", "O", "S", "F", "Cl", "Br", "H", "…"])
    left = (f'<div style="width: 660px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG0}; border-right: 1px solid {LINE}">'
            f'<div style="height: 44px; display: flex; align-items: center; gap: 2px; padding: 0 10px; background: {BG1}; border-bottom: 1px solid {LINE}">{tb}{spacer()}<div style="display: flex; gap: 4px">{els}</div></div>'
            f'<div style="position: relative; flex-grow: 1">{s2d}<div style="position: absolute; left: 14px; top: 12px">{chip("2D sketch", TEXT, "#16191Ccc")}</div></div></div>')
    atoms = styrene3d()
    sc = Scene("z", 1440 - 72 - 660 - 380, 470, yaw=0.7, pitch=0.9, roll=0.3, persp=0.25, fog=0.4)
    sc.add_atoms(atoms)
    v3 = sc.svg()
    confs = table(["#", "E rel.", "Method"], [["1", "0.00", "ETKDG v3 + GAFF2 min."], ["2", "[run]", "ETKDG v3 + GAFF2 min."], ["3", "[run]", "ETKDG v3 + GAFF2 min."]],
                  ["12%", "28%", "60%"], mono_cols=(0, 1), align_right=(1,), hl={0})
    mid = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
           f'<div style="position: relative; background: {BG0}; border-bottom: 1px solid {LINE}">{v3}<div style="position: absolute; left: 14px; top: 12px; display: flex; gap: 6px">{chip("3D preview · live", TEXT, "#16191Ccc")}{chip("updates as you draw", MUTED, "#16191Ccc")}</div></div>'
           f'<div style="padding: 14px 16px; display: flex; flex-direction: column; gap: 10px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Conformers</h3>{confs}</div></div>')
    right = (panel_head("Structure")
             + section("Identifiers", col(field("SMILES", "C=Cc1ccccc1"), field("InChIKey", "PPBRXRYQALVLMV-UHFFFAOYSA-N"), gap=10))
             + section("Properties", col(kv("Formula", "C<sub>8</sub>H<sub>8</sub>"), kv("Mol. weight", "104.15 g/mol"), kv("Stereocentres", "0"), kv("Valence check", f'<span style="color: {OK}">ok</span>'), gap=7))
             + section("3D embedding", col(row(select("Method", "ETKDG v3"), field("Conformers", "10"), gap=10), row(select("Clean with", "CAPS Field · GAFF2"), select("Hydrogens", "Add, VSEPR"), gap=10),
                                            cite("Riniker &amp; Landrum, <i>J. Chem. Inf. Model.</i> 55, 2562 (2015)"), gap=10))
             + footer(btn("Use as repeat unit", ic="grow"), btn("Insert into document", True, "plus")))
    center = left + mid
    return studio_shell("styrene.caps", "CAPS Studio — sketch to 3D", "", center, right,
                        ("<span>16 atoms after hydrogens</span>", "<span>Sketch: click places, drag grows, click bond cycles order</span>"))


# ------------------------------------------------------------------ Visual styles & surfaces
def visual():
    atoms, backbone, stereo = polystyrene(7, "atactic", seed=5)
    vdw = {"C": 1.70, "H": 1.20}
    # charge-like colouring on first units (GAFF-style sign pattern, illustrative)
    ramp_neg, ramp_mid, ramp_pos = "#3B7DD8", "#D9DCDF", "#E8893A"
    n_sf = 0
    xs = sorted(a["p"][0] for a in atoms)
    cut = xs[int(len(xs) * 0.42)]
    styled = []
    for a in atoms:
        if a["p"][0] < cut:
            n_sf += 1
            if a["e"] == "H":
                c = mix(ramp_mid, ramp_pos, 0.55)
            elif a.get("ar"):
                c = mix(ramp_mid, ramp_neg, 0.45)
            else:
                c = mix(ramp_mid, ramp_neg, 0.8)
            styled.append({**a, "r": vdw[a["e"]] * 0.95, "c": c, "sf": True})
        else:
            styled.append(a)
    sc = Scene("w", 1440 - 72 - 380, 900 - 44 - 48 - 26, yaw=-0.28, pitch=0.42, roll=0.05, persp=0.25, fog=0.55, outline=True)
    bonds = perceive_bonds(atoms)
    bonds = [(i, j) for i, j in bonds if not (styled[i].get("sf") and styled[j].get("sf"))]
    sc.add_atoms(styled, bonds)
    sc._prep()
    # clip plane indicator
    ov = f'<line x1="{sc.w * 0.47:.0f}" y1="40" x2="{sc.w * 0.47:.0f}" y2="{sc.h - 60}" stroke="{SEL}" stroke-width="1" stroke-dasharray="6 5" stroke-opacity="0.7"></line>'
    ov += label_pill(sc.w * 0.47 + 8, 52, "selection boundary · 3 units", SEL, border=SEL)
    ov += gizmo(46, sc.h - 46, -0.28, 0.42)
    svg = sc.svg(overlay=ov)
    legend = (f'<div style="position: absolute; left: 14px; bottom: 14px; width: 260px; padding: 10px 12px; background: #16191Ce6; border: 1px solid {LINE}; border-radius: 8px; display: flex; flex-direction: column; gap: 6px">'
              f'<span style="font-size: 11.5px; color: {MUTED}">Partial charge (e) · colour-blind-safe diverging</span>'
              f'<svg width="236" height="12" aria-hidden="true"><defs><linearGradient id="lg" x1="0" x2="1"><stop offset="0" stop-color="{ramp_neg}"></stop><stop offset="0.5" stop-color="{ramp_mid}"></stop><stop offset="1" stop-color="{ramp_pos}"></stop></linearGradient></defs><rect width="236" height="12" rx="3" fill="url(#lg)"></rect></svg>'
              f'<div style="display: flex; justify-content: space-between; font-family: {MONO}; font-size: 10.5px; color: {DIM}"><span>−0.20</span><span>0</span><span>+0.20</span></div></div>')
    vp = f'<main aria-label="3D view" style="position: relative; flex-grow: 1; min-width: 0; background: {BG0}">{svg}{hud(chip("Mixed styles: 2 selections", TEXT, "#16191Ccc"), chip("Illustrative charges", MUTED, "#16191Ccc"))}{legend}</main>'
    styles = [("Ball &amp; stick", "atom"), ("Stick", "bond"), ("Wireframe", "relax"), ("Space filling", "pack"), ("Polyhedra", "cube"), ("Ribbon", "grow")]
    sg = "".join(f'<button style="height: 58px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 5px; border-radius: 6px; border: 1px solid {ACC if n == "Space filling" else LINE}; background: {BG3 if n == "Space filling" else BG0}; font-size: 11.5px; color: {TEXT if n == "Space filling" else MUTED}; cursor: pointer">{icon(i, 18, ACC if n == "Space filling" else MUTED, 1.4)}<span>{n}</span></button>' for n, i in styles)
    right = (panel_head("Appearance", chip("Selection: units 1–3"))
             + section("Style", f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 6px">{sg}</div>')
             + section("Colour", col(row(select("Colour by", "Partial charge"), select("Palette", "Blue–orange (CB-safe)"), gap=10), gap=8))
             + section("Labels", col(row(check("Element", False), check("R/S", True), check("Atom type", False), check("Charge", False), gap=12), gap=8))
             + section("Surfaces", col(row(select("Surface", "Solvent-accessible"), field("Probe", "1.4", "Å"), gap=10),
                                        row(select("Map onto surface", "Electrostatic potential"), select("Opacity", "60 %"), gap=10),
                                        cite("Connolly 1983 · marching cubes, Lorensen &amp; Cline 1987"), gap=10))
             + section("Rendering", col(toggle("Ambient occlusion", True), toggle("Edge outlines", True), toggle("Depth cueing", True), gap=9))
             + footer(btn("Export image", ic="download"), btn("Apply", True, "check")))
    body = (topbar(["PS_atactic_DP40.caps"], "PS_atactic_DP40.caps") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">'
            + bar(tbtn("cursor", "Select", True), tbtn("lasso", "Lasso"), sep(), tbtn("eye", "Style", True, "Appearance"), tbtn("layers", "Surfaces", False, "Surfaces"), tbtn("cube", "Clip planes", False, "Clip planes"), spacer(), tbtn("cube", "Projection", False, "Perspective"), tbtn("download", "Export", False, "Export"))
            + f'<div style="flex-grow: 1; display: flex; min-height: 0">{vp}<aside aria-label="Appearance" style="width: 380px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{right}</aside></div></div></div>'
            + statusbar(f"<span>{len(atoms)} atoms · {n_sf} space-filling</span>", "<span>Impostor spheres · AO · outlines</span>"))
    return page("CAPS Studio — appearance", body)


# ------------------------------------------------------------------ Trajectory player
def trajectory():
    Lb = 34.0
    sc = Scene("t", 1440 - 72 - 400, 560, yaw=0.55, pitch=0.4, persp=0.3, fog=0.55)
    cols = ["#F0A83C", "#6CC4D8", "#E07A5F", "#9B7BD6", "#7CC784", "#D6A45E", "#E9ECEF", "#4C7BD9"]
    for k in range(8):
        rnd = random.Random(40 + k)
        st = (rnd.uniform(4, Lb - 4), rnd.uniform(4, Lb - 4), rnd.uniform(4, Lb - 4))
        pts = random_chain(st, 55, 500 + k, box=((0.8, 0.8, 0.8), (Lb - 0.8,) * 3), persistence=0.62)
        sc.add_tube(pts, cols[k] if k == 2 else mix(cols[k], BG0, 0.45), 0.55 if k == 2 else 0.45)
        if k == 2:
            sc.add_atoms([{"e": "C", "p": pts[0], "r": 0.8, "c": "#E07A5F"}, {"e": "C", "p": pts[-1], "r": 0.8, "c": "#E07A5F"}], [])
            sc.add_line(pts[0], pts[-1], ACC, 1.6, "5 4")
            ends = (pts[0], pts[-1])
    sc.add_box((0, 0, 0), (Lb, 0, 0), (0, Lb, 0), (0, 0, Lb), MUTED, 1.2, None, 0.8)
    sc._prep()
    a, b = sc.P(ends[0]), sc.P(ends[1])
    ree = dist(*ends)
    ov = label_pill((a[0] + b[0]) / 2 + 10, (a[1] + b[1]) / 2, f"R_ee chain 3 · {ree:.1f} Å", ACC, border=ACC)
    view = sc.svg(overlay=ov)
    # timeline
    N = 200
    cur = 128
    ticks = "".join(f'<span style="position: absolute; left: {i / N * 100:.2f}%; top: 0; width: 1px; height: {10 if i % 50 == 0 else 5}px; background: {DIM}"></span>' for i in range(0, N + 1, 10))
    marks = "".join(f'<span title="checkpoint" style="position: absolute; left: {i / N * 100:.2f}%; top: 14px; width: 6px; height: 6px; margin-left: -3px; border-radius: 50%; background: {SEL}"></span>' for i in (50, 100, 150))
    timeline = (f'<div style="flex-shrink: 0; display: flex; flex-direction: column; gap: 10px; padding: 12px 18px; background: {BG1}; border-top: 1px solid {LINE}">'
                + row(tbtn("undo", "Previous frame"), f'<button aria-label="Pause" style="width: 40px; height: 40px; border-radius: 50%; background: {ACC}; border: 0; display: flex; align-items: center; justify-content: center; cursor: pointer">{icon("pause", 18, ACC_INK, 2.2)}</button>',
                      tbtn("redo", "Next frame"), f'<span style="font-family: {MONO}; font-size: 13px; margin-left: 10px">frame {cur} / {N}</span>',
                      f'<span style="font-family: {MONO}; font-size: 12px; color: {MUTED}">t = {cur * 25 / 1000:.2f} ns</span>', spacer(),
                      seg(["0.5×", "1×", "2×", "4×"], "1×"), toggle_small("Loop"), toggle_small("Smooth"), btn("Export movie", ic="download", small=True), gap=8)
                + f'<div style="position: relative; height: 26px"><div style="position: absolute; left: 0; right: 0; top: 4px; height: 4px; border-radius: 2px; background: {BG3}"></div>'
                  f'<div style="position: absolute; left: 0; width: {cur / N * 100:.1f}%; top: 4px; height: 4px; border-radius: 2px; background: {ACC}"></div>'
                  f'<span style="position: absolute; left: {cur / N * 100:.1f}%; top: -2px; width: 16px; height: 16px; margin-left: -8px; border-radius: 50%; background: {TEXT}; border: 3px solid {ACC}; box-sizing: border-box"></span>{marks}</div>'
                + '</div>')
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">'
              f'<div style="position: relative; flex-grow: 1">{view}{hud(chip("dyn-12 · NPT 298 K", TEXT, "#16191Ccc"), chip("wrapped · colour: molecule", MUTED, "#16191Ccc"))}</div>{timeline}</div>')
    xs = [i * 0.025 for i in range(0, 201)]
    rho = [(x, 1.035 + 0.006 * math.sin(x * 9) + 0.003 * math.sin(x * 23)) for x in xs]
    T = [(x, 298 + 4 * math.sin(x * 31) + 2 * math.sin(x * 7)) for x in xs]
    rg = [(x, 15.2 + 0.8 * math.sin(x * 2.2) + 0.3 * math.sin(x * 11)) for x in xs]
    tc = cur * 0.025

    def cplot(series, yr, yt, yl, c):
        p = plot(360, 150, [(series, c, 1.5, None), ([(tc, yr[0]), (tc, yr[1])], TEXT, 1, "3 3")], (0, 5), yr, [0, 1, 2, 3, 4, 5], yt, "t (ns)", yl, pad=(44, 10, 16, 28))
        return p
    right = (panel_head("Linked plots", chip("cursor follows frame"))
             + section("Density", cplot(rho, (1.02, 1.05), [1.02, 1.035, 1.05], "g/cm³", ACC), chip("1.041", TEXT, BG2, True), pad=10, gap=6)
             + section("Temperature", cplot(T, (290, 306), [290, 298, 306], "K", SEL), chip("297.6", TEXT, BG2, True), pad=10, gap=6)
             + section("R<sub>g</sub> chain 3", cplot(rg, (14, 16.5), [14, 15, 16], "Å", "#E07A5F"), chip("15.4", TEXT, BG2, True), pad=10, gap=6)
             + f'<div style="padding: 10px 14px; font-size: 11.5px; color: {DIM}">Plot values are illustrative.</div>')
    body = (topbar(["PS_cell_20x40.caps"], "PS_cell_20x40.caps") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">'
            + bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), tbtn("pin", "Pin monitor"), sep(), tbtn("cube", "Wrap", False, "Wrap molecules"), tbtn("layers", "Images", False, "Periodic images"), spacer(), select_inline("Trajectory", "dyn-12 · traj.dcd · 200 frames"))
            + f'<div style="flex-grow: 1; display: flex; min-height: 0">{center}<aside aria-label="Plots" style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{right}</aside></div></div></div>'
            + statusbar("<span>12 840 atoms · 200 frames · 25 ps/frame</span>", "<span>Streaming from disk · 2 frames cached ahead</span>"))
    return page("CAPS Studio — trajectory", body)


def toggle_small(label):
    return f'<button role="switch" aria-checked="true" style="height: 30px; padding: 0 10px; display: flex; align-items: center; gap: 6px; background: transparent; border: 1px solid {LINE}; border-radius: 6px; font-size: 12px; color: {MUTED}; cursor: pointer">{dot(ACC, 7)}{label}</button>'


def select_inline(label, value):
    return (f'<button aria-label="{label}" style="height: 34px; padding: 0 10px; display: flex; align-items: center; gap: 8px; background: {BG2}; border: 1px solid {LINE}; border-radius: 6px; font-size: 12px; color: {TEXT}; cursor: pointer">'
            f'<span style="color: {DIM}">{label}</span><span style="font-family: {MONO}">{value}</span>{icon("chev", 13, DIM)}</button>')


# ------------------------------------------------------------------ Torsion scan & conformers
def torsion():
    # n-butane, TraPPE-UA torsion (Martin & Siepmann 1998); energies in kcal/mol
    c1, c2, c3 = 355.03, -68.19, 791.32

    def U(phi_deg):
        p = math.radians(phi_deg)
        return (c1 * (1 + math.cos(p)) + c2 * (1 - math.cos(2 * p)) + c3 * (1 + math.cos(3 * p))) * KCAL_PER_K

    curve = [(d, U(d)) for d in range(-180, 181, 2)]
    # locate minima / maxima numerically
    fine = [(d / 10, U(d / 10)) for d in range(-1800, 1801)]
    mins = [fine[i] for i in range(1, len(fine) - 1) if fine[i][1] < fine[i - 1][1] and fine[i][1] < fine[i + 1][1]]
    maxs = [fine[i] for i in range(1, len(fine) - 1) if fine[i][1] > fine[i - 1][1] and fine[i][1] > fine[i + 1][1]]
    g = min((m for m in mins if m[0] > 0), key=lambda m: m[1])
    top = max(m[1] for m in maxs)
    ymax = math.ceil(top + 0.5)
    scan_pts = [(d, U(d)) for d in range(-180, 181, 15)]
    pl = plot(560, 300, [(curve, ACC, 2, None)], (-180, 180), (0, ymax), [-180, -120, -60, 0, 60, 120, 180], list(range(0, ymax + 1, 1)),
              "dihedral φ (degrees)", "E (kcal/mol)", pad=(44, 14, 20, 32), markers=[(a, b, SEL) for a, b in scan_pts])
    # butane UA 3D at gauche
    phi = math.radians(g[0])
    b1 = (0.0, 0.0, 0.0)
    b2 = (1.54, 0.0, 0.0)
    th = math.radians(114.0)
    b0 = (1.54 * math.cos(math.pi - th) * -1 + 0.0, 1.54 * math.sin(th), 0.0)
    b0 = (b1[0] - 1.54 * math.cos(th - math.pi / 2 - math.pi / 2), b1[1] + 1.54 * math.sin(th), 0.0)
    # place atoms with standard internal-coordinate construction
    def place(a, b, c, bond, ang, dih):
        bc = norm(sub(c, b))
        n = norm(cross(sub(b, a), bc))
        m = cross(n, bc)
        d2 = (-bond * math.cos(math.radians(ang)), bond * math.sin(math.radians(ang)) * math.cos(math.radians(dih)), bond * math.sin(math.radians(ang)) * math.sin(math.radians(dih)))
        return add(c, add(add(mul(bc, d2[0]), mul(m, d2[1])), mul(n, d2[2])))
    A = (0.0, 0.0, 0.0)
    B = (1.54, 0.0, 0.0)
    C = add(B, (1.54 * math.cos(math.radians(180 - 114)), 1.54 * math.sin(math.radians(180 - 114)), 0.0))
    D = place(A, B, C, 1.54, 114.0, g[0])
    measured = dihedral(A, B, C, D)
    beads = [{"e": "C", "p": p, "r": 0.75, "c": c} for p, c in ((A, "#9AA1A8"), (B, ACC), (C, ACC), (D, "#9AA1A8"))]
    sc = Scene("q", 1440 - 72 - 640, 470, yaw=0.4, pitch=0.3, persp=0.2, fog=0.3, bond_w=0.3)
    sc.add_atoms(beads, [(0, 1), (1, 2), (2, 3)])
    sc._prep()
    pB, pC = sc.P(B), sc.P(C)
    ov = label_pill((pB[0] + pC[0]) / 2 - 30, (pB[1] + pC[1]) / 2 + 30, f"φ = {measured:.1f}°", ACC, border=ACC)
    ov += f'<g>{label_pill(sc.P(D)[0] + 16, sc.P(D)[1], "rotating side", SEL, border=SEL)}</g>'
    view = sc.svg(overlay=ov)
    E_t = U(180.0)
    rows = [["1", "trans", "180.0°", f"{E_t:.2f}"], ["2", "gauche+", f"{g[0]:.1f}°", f"{g[1] - E_t:.2f}"], ["3", "gauche−", f"{-g[0]:.1f}°", f"{g[1] - E_t:.2f}"]]
    conf = table(["#", "State", "φ", "ΔE (kcal/mol)"], rows, ["10%", "30%", "28%", "32%"], mono_cols=(0, 2, 3), align_right=(2, 3), hl={0})
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">'
              f'<div style="position: relative">{view}{hud(chip("n-butane · united atom", TEXT, "#16191Ccc"), chip("TraPPE-UA", MUTED, "#16191Ccc"))}</div>'
              f'<div style="flex-grow: 1; padding: 16px 20px; background: {BG1}; border-top: 1px solid {LINE}; display: flex; flex-direction: column; gap: 10px">'
              f'<h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Conformers found</h3>{conf}'
              f'<div style="font-size: 11.5px; color: {DIM}">Energies computed from the TraPPE-UA torsion; the gauche/trans gap and barrier are values of that potential, not measurements.</div></div></div>')
    right = (f'<div style="width: 640px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}">'
             + panel_head("Torsion scan", chip(f"barrier {top:.2f} kcal/mol", ACC, "#3A2C14"))
             + f'<div style="padding: 14px 16px; display: flex; flex-direction: column; gap: 12px; border-bottom: 1px solid {LINE}">{pl}'
             + row(chip(f'{dot(ACC)} potential'), chip(f'{dot(SEL)} scan points · 15° step'), spacer(), chip(f"gauche at ±{g[0]:.1f}°", TEXT, BG3, True), gap=8) + '</div>'
             + section("Scan setup", col(row(field("Atoms", "C1 C2 C3 C4"), select("Moving side", "C4 side"), gap=10),
                                          row(field("From", "−180", "°"), field("To", "180", "°"), field("Step", "15", "°"), select("At each step", "Relax others · L-BFGS"), gap=8),
                                          cite("Torsion potential: Martin &amp; Siepmann, <i>J. Phys. Chem. B</i> 102, 2569 (1998)"), gap=10))
             + footer(btn("Export CSV", ic="download"), btn("Run scan", True, "play")) + '</div>')
    body = (topbar(["butane_scan.caps"], "butane_scan.caps") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">'
            + bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), tbtn("rotate", "Set dihedral", True, "Set dihedral"), sep(), tbtn("chart", "Scan", True, "Torsion scan"), tbtn("layers", "Conformers", False, "Conformer search"), spacer())
            + f'<div style="flex-grow: 1; display: flex; min-height: 0">{center}{right}</div></div></div>'
            + statusbar("<span>4 beads · 3 bonds</span>", "<span>Scan uses CAPS Field + CAPS Relax</span>"))
    return page("CAPS Studio — torsion scan", body)


# ------------------------------------------------------------------ Settings
def swatches(cs):
    return "".join(f'<span style="width: 18px; height: 18px; border-radius: 4px; background: {c}"></span>' for c in cs)


def settings():
    navs = [("General", "gear"), ("Appearance", "eye"), ("Input &amp; shortcuts", "terminal"), ("Rendering", "cube"), ("Force fields", "tag"),
            ("Compute &amp; remote", "server"), ("Accessibility", "check"), ("Python &amp; scripting", "copy")]
    nv = "".join(f'<a href="{ {"Appearance": "Settings.dc.html", "Compute &amp; remote": "RemoteCompute.dc.html"}.get(n, "#") }" aria-current="{"page" if n == "Appearance" else "false"}" style="display: flex; align-items: center; gap: 10px; height: 36px; padding: 0 12px; border-radius: 6px; text-decoration: none; font-size: 13px; background: {BG3 if n == "Appearance" else "transparent"}; color: {TEXT if n == "Appearance" else MUTED}">{icon(i, 16, ACC if n == "Appearance" else DIM)}{n}</a>' for n, i in navs)

    def theme_card(name, bg, panel, text, on):
        return (f'<button aria-pressed="{"true" if on else "false"}" style="display: flex; flex-direction: column; gap: 8px; padding: 10px; background: {BG0}; border: 1.5px solid {ACC if on else LINE}; border-radius: 10px; cursor: pointer; text-align: left">'
                f'<svg width="190" height="104" viewBox="0 0 190 104" aria-hidden="true"><rect width="190" height="104" rx="6" fill="{bg}"></rect><rect width="190" height="14" rx="3" fill="{panel}"></rect>'
                f'<rect x="0" y="14" width="22" height="90" fill="{panel}"></rect><rect x="140" y="14" width="50" height="90" fill="{panel}"></rect>'
                f'<circle cx="70" cy="56" r="11" fill="#8E959C"></circle><circle cx="92" cy="46" r="8" fill="#E5534B"></circle><circle cx="96" cy="68" r="9" fill="#4C7BD9"></circle>'
                f'<line x1="70" y1="56" x2="92" y2="46" stroke="{text}" stroke-width="2"></line><line x1="70" y1="56" x2="96" y2="68" stroke="{text}" stroke-width="2"></line>'
                f'<rect x="146" y="22" width="36" height="5" rx="2" fill="{text}" fill-opacity="0.5"></rect><rect x="146" y="32" width="28" height="5" rx="2" fill="{text}" fill-opacity="0.3"></rect><rect x="146" y="42" width="32" height="5" rx="2" fill="{ACC}"></rect></svg>'
                f'<span style="font-size: 12.5px; color: {TEXT}">{name}</span></button>')
    themes = row(theme_card("Graphite (dark)", "#0F1113", "#16191C", "#E8E6E1", True), theme_card("Paper (light)", "#F4F2EE", "#FFFFFF", "#1E2226", False), theme_card("Match system", "#0F1113", "#F4F2EE", "#A5ABB1", False), gap=12)
    pal = [("Element (CAPS default)", ["#8E959C", "#E9ECEF", "#E5534B", "#4C7BD9", "#E3C74A", "#57B26A"]),
           ("Okabe–Ito (colour-blind safe)", ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"]),
           ("Monochrome + shape", ["#F2F2F2", "#C8C8C8", "#9E9E9E", "#747474", "#4A4A4A", "#222222"])]
    pr = "".join(f'<button aria-pressed="{"true" if i == 1 else "false"}" style="display: flex; align-items: center; gap: 12px; height: 44px; padding: 0 12px; background: {BG3 if i == 1 else BG0}; border: 1px solid {ACC if i == 1 else LINE}; border-radius: 7px; cursor: pointer">'
                 f'<span style="display: flex; gap: 3px">{swatches(cs)}</span>'
                 f'<span style="font-size: 12.5px; color: {TEXT}">{n}</span></button>' for i, (n, cs) in enumerate(pal))
    content = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 22px; padding: 28px 36px">'
               f'<div style="display: flex; flex-direction: column; gap: 4px"><h1 style="margin: 0; font-size: 22px; font-weight: 600">Appearance</h1><p style="margin: 0; font-size: 13px; color: {MUTED}">Theme, scale and colour palettes. Changes apply live and are saved per user.</p></div>'
               + col(f'<h2 style="margin: 0; font-size: 13px; font-weight: 600">Theme</h2>', themes, gap=10)
               + col(f'<h2 style="margin: 0; font-size: 13px; font-weight: 600">Interface scale</h2>',
                     row(f'<div style="flex-grow: 1; max-width: 420px; position: relative; height: 18px"><div style="position: absolute; top: 7px; left: 0; right: 0; height: 4px; border-radius: 2px; background: {BG3}"></div><div style="position: absolute; top: 7px; left: 0; width: 40%; height: 4px; border-radius: 2px; background: {ACC}"></div><span style="position: absolute; left: 40%; top: 0; width: 18px; height: 18px; margin-left: -9px; border-radius: 50%; background: {TEXT}"></span></div>',
                         f'<span style="font-family: {MONO}; font-size: 13px">100 %</span>', chip("Ctrl + / Ctrl −"), gap=14), gap=10)
               + col(f'<h2 style="margin: 0; font-size: 13px; font-weight: 600">Element and category palette</h2>', f'<div style="display: flex; flex-direction: column; gap: 6px; max-width: 520px">{pr}</div>',
                     f'<div style="font-size: 12px; color: {MUTED}">Plots, chain colours and selection highlights follow the same palette. Categories also differ in lightness, so they read in greyscale.</div>', gap=10)
               + '</div>')
    side = (f'<div style="width: 340px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 28px 28px 28px 0">'
            + card("Accessibility", col(toggle("Reduce motion", False), toggle("High-contrast outlines", False), toggle("Always show focus rings", True), toggle("Announce job progress to screen readers", True), gap=10))
            + card("Shortcut conflicts", col(f'<div style="display: flex; gap: 10px; align-items: flex-start; font-size: 12.5px">{icon("alert", 16, WARN)}<span>⌘ ⇧ C is bound to <b style="font-weight: 600">Clean geometry</b> and a recorded macro.</span></div>',
                                             row(btn("Resolve", small=True), btn("Open shortcut editor", ic="terminal", small=True), gap=6), gap=10), chip("1", WARN, "#3A2C14"))
            + '</div>')
    body = (topbar(["Settings"], "Settings") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("")
            + f'<nav aria-label="Settings sections" style="width: 240px; flex-shrink: 0; display: flex; flex-direction: column; gap: 2px; padding: 20px 12px; background: {BG1}; border-right: 1px solid {LINE}">{nv}</nav>'
            + f'<div style="flex-grow: 1; display: flex; min-width: 0">{content}{side}</div></div>'
            + statusbar("<span>Settings saved to ~/.caps/settings.json</span>", "<span>Export · Import · Reset</span>"))
    return page("CAPS — settings", body)
