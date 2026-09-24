from lib import *
from mols import polystyrene


def build():
    VW, VH = 788, 596
    atoms, backbone, stereo = polystyrene(7, "atactic", seed=5)
    sc = Scene("m", VW, VH, yaw=-0.28, pitch=0.42, roll=0.05, persp=0.25, fog=0.6)
    sc.add_atoms(atoms)
    sel = [stereo[3], backbone[7], backbone[8], stereo[4]]
    # selected fragment: phenyl ring on stereo[3]
    ring_start = None
    for i, a in enumerate(atoms):
        if a.get("ar") and dist(a["p"], atoms[stereo[3]]["p"]) < 1.6:
            ring_start = i
    sc.halo = set([stereo[3]] + list(range(ring_start, ring_start + 6)))
    sc._prep()
    # dihedral monitor on backbone 5-6-7-8
    b = [atoms[backbone[k]]["p"] for k in (5, 6, 7, 8)]
    phi = dihedral(*b)
    d_ring = dist(atoms[stereo[3]]["p"], atoms[ring_start]["p"])
    ang = angle(atoms[backbone[5]]["p"], atoms[backbone[6]]["p"], atoms[backbone[7]]["p"])
    ov = []
    pts = [sc.P(p) for p in b]
    ov.append(f'<polyline points="{" ".join(f"{p[0]:.1f},{p[1]:.1f}" for p in pts)}" fill="none" stroke="{ACC}" stroke-width="1.6" stroke-dasharray="4 3"></polyline>')
    mx = sum(p[0] for p in pts[1:3]) / 2
    my = sum(p[1] for p in pts[1:3]) / 2
    ov.append(label_pill(mx + 14, my - 40, f"φ {phi:.1f}°", ACC, border=ACC))
    # stereo labels
    rs_names = ["R", "S", "S", "R", "S", "R", "R"]
    for n, k in enumerate(stereo[1:6]):
        x, y, _, _ = sc.screen(k)
        ov.append(label_pill(x - 10, y + 22, rs_names[n], SEL if k == stereo[3] else MUTED, border=SEL if k == stereo[3] else LINE, fs=10.5))
    ov.append(f'<g transform="translate(0,0)">{gizmo(46, VH - 46, -0.28, 0.42)}</g>')
    svg = sc.svg(overlay="".join(ov))

    toolbar = (f'<div role="toolbar" aria-label="Builder tools" style="height: 48px; flex-shrink: 0; display: flex; align-items: center; gap: 2px; padding: 0 10px; background: {BG1}; border-bottom: 1px solid {LINE}">'
               + tbtn("cursor", "Select", True) + tbtn("lasso", "Lasso select") + tbtn("move", "Translate") + tbtn("rotate", "Rotate") + sep()
               + tbtn("atom", "Place atom") + tbtn("bond", "Draw bond") + tbtn("hex", "Fragment") + tbtn("ring", "Fuse ring")
               + f'<button aria-label="Element: carbon" style="height: 36px; padding: 0 10px; margin-left: 4px; display: flex; align-items: center; gap: 6px; background: {BG2}; border: 1px solid {LINE}; border-radius: 6px; cursor: pointer">'
                 f'<span style="width: 22px; height: 22px; border-radius: 4px; background: #8E959C; color: #0F1113; font-weight: 700; font-size: 12px; display: flex; align-items: center; justify-content: center">C</span><span style="font-size: 12px; color: {MUTED}">sp³</span>{icon("chev", 13, DIM)}</button>'
               + sep() + tbtn("ruler", "Measure") + tbtn("pin", "Pin monitor") + tbtn("mirror", "Mirror / invert") + sep()
               + tbtn("wand", "Auto-clean", True, "Auto-clean")
               + tbtn("undo", "Undo") + tbtn("redo", "Redo")
               + '<div style="flex-grow: 1"></div>'
               + tbtn("eye", "Style", False, "Ball &amp; stick") + tbtn("layers", "Colour", False, "Colour: element") + tbtn("cube", "Projection", False, "Perspective")
               + '</div>')

    tree_items = [
        (0, "folder", "PS-tacticity study", MUTED, False),
        (1, "cube", "PS_atactic_DP40.caps", TEXT, True),
        (2, "hex", "Chain 1 · 642 atoms", MUTED, False),
        (2, "hex", "Chain 2 · 642 atoms", MUTED, False),
        (2, "tag", "Selection: ring-7", SEL, False),
        (1, "cube", "PS_isotactic_DP40.caps", MUTED, False),
        (1, "cube", "styrene_monomer.caps", MUTED, False),
        (0, "folder", "Jobs", MUTED, False),
        (1, "grow", "Grow · amorphous cell #3", MUTED, False),
    ]
    tree = []
    for lvl, ic, name, c, on in tree_items:
        bg = f"background: {BG3}; " if on else ""
        tree.append(f'<button style="{bg}height: 26px; width: 100%; display: flex; align-items: center; gap: 7px; padding: 0 8px 0 {10 + lvl * 14}px; border: 0; border-radius: 5px; background-color: {BG3 if on else "transparent"}; color: {c}; font-size: 12.5px; cursor: pointer; text-align: left">'
                    f'{icon("chev" if ic == "folder" else "chevr", 11, DIM) if lvl < 2 else "<span style=&quot;width: 11px&quot;></span>".replace("&quot;", chr(34))}{icon(ic, 14, ACC if on else DIM)}<span style="overflow: hidden; white-space: nowrap; text-overflow: ellipsis">{name}</span></button>')
    frags = [("Benzene", "ring"), ("Cyclohexane", "hex"), ("Ester", "link"), ("Amide", "link"), ("Methyl", "atom"), ("Vinyl", "bond"),
             ("Water", "flask"), ("Na⁺ / Cl⁻", "atom"), ("Styrene", "hex")]
    fr = []
    for n, ic in frags:
        fr.append(f'<button style="height: 62px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 5px; background: {BG2}; border: 1px solid {LINE}; border-radius: 6px; font-size: 11px; color: {MUTED}; cursor: pointer">{icon(ic, 20, TEXT, 1.4)}<span>{n}</span></button>')
    left = (f'<aside aria-label="Project" style="width: 260px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-right: 1px solid {LINE}">'
            + panel_head("Project", tbtn("plus", "Add to project"))
            + f'<div style="padding: 8px 6px; display: flex; flex-direction: column; gap: 1px">{"".join(tree)}</div>'
            + f'<div style="border-top: 1px solid {LINE}">' + panel_head("Fragments", chip("Monomers ▾")) + '</div>'
            + f'<div style="padding: 10px"><div style="display: flex; align-items: center; gap: 8px; height: 30px; padding: 0 8px; margin-bottom: 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 5px; color: {DIM}; font-size: 12px">{icon("search", 14, DIM)}<span>Search 1,240 fragments</span></div>'
            + f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 6px">{"".join(fr)}</div></div>'
            + '</aside>')

    hud_l = (f'<div style="position: absolute; left: 14px; top: 12px; display: flex; gap: 6px">'
             + chip(f'{dot(SEL, 7)} 7 atoms selected', TEXT, "#16191Ccc") + chip("Tacticity: atactic · 3 m / 3 r dyads", MUTED, "#16191Ccc") + '</div>')
    monitors = [("φ C5–C6–C7–C8", f"{phi:.1f}°", ACC), ("∠ C5–C6–C7", f"{ang:.1f}°", TEXT), ("d C(α)–C(ipso)", f"{d_ring:.3f} Å", TEXT)]
    mon = "".join(f'<div style="display: flex; justify-content: space-between; gap: 18px; font-size: 12px"><span style="color: {MUTED}">{a}</span><span style="font-family: {MONO}; color: {c}">{b}</span></div>' for a, b, c in monitors)
    hud_r = (f'<div style="position: absolute; right: 14px; top: 12px; width: 230px; padding: 10px 12px; display: flex; flex-direction: column; gap: 6px; background: #16191Ce6; border: 1px solid {LINE}; border-radius: 8px">'
             f'<div style="display: flex; align-items: center; gap: 6px; font-size: 11px; font-weight: 600; letter-spacing: 0.06em; text-transform: uppercase; color: {MUTED}">{icon("pin", 13, ACC)}Live monitors</div>{mon}</div>')
    clean = (f'<div style="position: absolute; left: 50%; bottom: 14px; transform: translateX(-50%); display: flex; align-items: center; gap: 10px; height: 34px; padding: 0 14px; background: #16191Ce6; border: 1px solid {LINE}; border-radius: 17px; font-size: 12px; color: {MUTED}">'
             f'{icon("wand", 14, ACC)}<span>Auto-clean on · GAFF2 · L-BFGS</span><span style="font-family: {MONO}; color: {OK}">converged · |F|<sub>max</sub> 0.08 kcal/mol/Å</span></div>')
    viewport = (f'<main aria-label="3D viewport" style="position: relative; flex-grow: 1; min-width: 0; background: {BG0}">{svg}{hud_l}{hud_r}{clean}</main>')

    rows_valid = [
        (WARN, "alert", "2 undefined stereocentres at chain ends", "Set R/S"),
        (WARN, "alert", "C31–H44 1.21 Å is longer than usual (1.09 Å)", "Clean"),
        (OK, "check", "Charge balanced · net 0.000 e", None),
    ]
    vr = []
    for c, ic, t, fix in rows_valid:
        f = f'<button style="margin-left: auto; height: 24px; padding: 0 8px; background: transparent; border: 1px solid {LINE}; border-radius: 5px; font-size: 11.5px; color: {ACC}; cursor: pointer; flex-shrink: 0">{fix}</button>' if fix else ""
        vr.append(f'<div style="display: flex; align-items: center; gap: 8px; font-size: 12px; line-height: 1.35">{icon(ic, 15, c)}<span style="color: {TEXT}">{t}</span>{f}</div>')
    right = (f'<aside aria-label="Inspector" style="width: 320px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">'
             + tabs(["Atom", "Molecule", "Force field"], "Atom")
             + section("Selected atom · C12", col(
                 row(field("Element", "C", w=70), select("Hybridisation", "sp³ tetrahedral"), gap=8),
                 row(field("x", "4.182", "Å"), field("y", "−1.037", "Å"), field("z", "0.566", "Å"), gap=6),
                 kv("Stereo (CIP)", f'<span style="color: {SEL}">R</span> · priority c1&gt;c7&gt;c11&gt;H'),
                 kv("Formal charge", "0"), gap=10), tbtn("dots", "More"))
             + section("Force-field type", col(
                 row(chip("GAFF2", ACC_INK, ACC), chip("c3", TEXT, BG3, True), chip("q −0.0960 e", MUTED, BG2, True), gap=6),
                 f'<div style="padding: 9px 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px; font-size: 12px; line-height: 1.5; color: {MUTED}">'
                 f'<span style="color: {TEXT}">Why c3?</span> Rule <span style="font-family: {MONO}; color: {TEXT}">gaff2.c3</span> matched <span style="font-family: {MONO}; color: {TEXT}">[CX4]</span> — sp³ carbon, 4 connections, no heteroatom neighbour.</div>',
                 gap=10))
             + section("Validation", col(*vr, gap=9), chip("2 warnings", WARN, "#3A2C14"))
             + '</aside>')

    hist = [("12", "fragment.attach  Benzene → C12 (open valence)", False), ("13", "stereo.set  C12 → R", False),
            ("14", "measure.dihedral  C5 C6 C7 C8  --pin", False), ("15", "clean.run  ff=gaff2 method=lbfgs", True)]
    hr = "".join(f'<div style="display: flex; gap: 10px; align-items: center; height: 24px; padding: 0 8px; border-radius: 4px; {"background: " + BG3 + "; " if on else ""}font-family: {MONO}; font-size: 11.5px; color: {TEXT if on else MUTED}"><span style="color: {DIM}; width: 18px">{n}</span><span>{esc(t)}</span></div>' for n, t, on in hist)
    code = [("#", "Recorded from Studio · replayable"), ("", "import caps"), ("", 'doc = caps.open("PS_atactic_DP40.caps")'),
            ("", 'doc.fragment.attach("benzene", at=doc.atom(12))'), ("", 'doc.stereo.set(doc.atom(12), "R")'),
            ("", 'doc.clean(ff="gaff2", method="lbfgs")'), (">>>", "")]
    cr = []
    for p, t in code:
        if p == "#":
            cr.append(f'<div style="color: {DIM}"># {t}</div>')
        elif p == ">>>":
            cr.append(f'<div><span style="color: {ACC}">&gt;&gt;&gt;</span> <span style="display: inline-block; width: 7px; height: 13px; background: {TEXT}; vertical-align: -2px"></span></div>')
        else:
            cr.append(f'<div style="color: {TEXT}">{esc(t)}</div>')
    dock = (f'<section aria-label="History and console" style="height: 190px; flex-shrink: 0; display: flex; background: {BG1}; border-top: 1px solid {LINE}">'
            f'<div style="width: 44%; display: flex; flex-direction: column; border-right: 1px solid {LINE}">'
            + tabs(["History", "Selection sets", "Log"], "History")
            + f'<div style="padding: 8px; display: flex; flex-direction: column; gap: 2px">{hr}</div></div>'
            f'<div style="flex-grow: 1; display: flex; flex-direction: column">'
            + tabs(["Python console", "Command line"], "Python console")
            + f'<div style="padding: 10px 14px; font-family: {MONO}; font-size: 12px; line-height: 1.6">{"".join(cr)}</div></div></section>')

    body = (topbar(["PS_atactic_DP40.caps", "PET_cell.caps", "PE_crystal.caps"], "PS_atactic_DP40.caps")
            + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">' + toolbar
            + '<div style="flex-grow: 1; display: flex; min-height: 0">' + left + viewport + right + '</div>' + dock + '</div></div>'
            + statusbar(f'<span>{len(atoms)} atoms · {len(sc.bonds)} bonds</span><span>Selection 7</span><span>Cell: none (non-periodic)</span>',
                        '<span>Mouse: select · drag rotates · ⇧ drag pans</span><span>Double precision</span><span>60 fps</span>'))
    return page("CAPS Studio — builder", body)
