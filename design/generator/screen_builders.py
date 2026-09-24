import math, random
from lib import *
from mols import *

BUILDERS = ["Polymer", "Crystal", "Surface", "Nanostructure", "Biomolecule", "Solvation"]
BUILDER_LINK = {"Polymer": "PolymerBuilder.dc.html", "Crystal": "CrystalBuilder.dc.html", "Surface": "SurfaceBuilder.dc.html",
                "Nanostructure": "NanoBuilder.dc.html", "Biomolecule": "BioBuilder.dc.html", "Solvation": "SolvationBuilder.dc.html"}


def builder_shell(doc, active, scene_svg, panel, vp_overlay="", status=("", "")):
    strip = []
    for b in BUILDERS:
        on = b == active
        st = f"background: {BG3}; color: {TEXT}; border: 1px solid {LINE}" if on else f"background: transparent; color: {MUTED}; border: 1px solid transparent"
        strip.append(f'<a href="{BUILDER_LINK[b]}" aria-current="{"page" if on else "false"}" style="height: 30px; padding: 0 12px; display: flex; align-items: center; border-radius: 6px; font-size: 12.5px; text-decoration: none; {st}">{b}</a>')
    head = (f'<div style="height: 48px; flex-shrink: 0; display: flex; align-items: center; gap: 4px; padding: 0 12px; background: {BG1}; border-bottom: 1px solid {LINE}">'
            f'<span style="font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {DIM}; margin-right: 8px">Builders</span>'
            + "".join(strip) + '<div style="flex-grow: 1"></div>'
            + tbtn("undo", "Undo") + tbtn("redo", "Redo") + sep() + tbtn("eye", "Style", False, "Ball &amp; stick") + tbtn("cube", "Projection", False, "Orthographic") + '</div>')
    vp = f'<main aria-label="3D preview" style="position: relative; flex-grow: 1; min-width: 0; background: {BG0}">{scene_svg}{vp_overlay}</main>'
    pn = f'<aside aria-label="{active} builder" style="width: 420px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{panel}</aside>'
    body = (topbar([doc], doc) + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">' + head
            + '<div style="flex-grow: 1; display: flex; min-height: 0">' + vp + pn + '</div></div></div>'
            + statusbar(status[0], status[1]))
    return page(f"CAPS Studio — {active} builder", body)


VPW, VPH = 1440 - 72 - 420, 900 - 44 - 48 - 26


def footer(*btns):
    return (f'<div style="margin-top: auto; display: flex; gap: 8px; justify-content: flex-end; padding: 12px 14px; border-top: 1px solid {LINE}; background: {BG1}">'
            + "".join(btns) + '</div>')


def hud(*chips, pos="left: 14px; top: 12px"):
    return f'<div style="position: absolute; {pos}; display: flex; gap: 6px; flex-wrap: wrap">' + "".join(chips) + '</div>'


# ---------------------------------------------------------------- Polymer
def polymer():
    atoms, backbone, stereo = polystyrene(9, "syndiotactic", curve=True)
    sc = Scene("p", VPW, VPH - 150, yaw=-0.12, pitch=0.5, persp=0.2, fog=0.55)
    sc.add_atoms(atoms)
    # head/tail markers
    sc._prep()
    h = sc.screen(backbone[0]); t = sc.screen(backbone[-1])
    ov = (label_pill(h[0] - 70, h[1] - 26, "head *", ACC, border=ACC) + label_pill(t[0] + 12, t[1] - 26, "tail *", ACC, border=ACC))
    svg = sc.svg(overlay=ov)
    # sequence strip under viewport
    cells = []
    seq = ["A"] * 20 + ["B"] * 20
    for i, s in enumerate(seq):
        c = ACC if s == "A" else SEL
        cells.append(f'<span title="{"styrene" if s == "A" else "MMA"} {i + 1}" style="flex-grow: 1; height: 26px; border-radius: 3px; background: {c}; opacity: {1 if i < 18 or i > 21 else 0.55}"></span>')
    seqbar = (f'<div style="position: absolute; left: 0; right: 0; bottom: 0; height: 150px; box-sizing: border-box; padding: 14px 18px; border-top: 1px solid {LINE}; background: {BG1}; display: flex; flex-direction: column; gap: 10px">'
              + row(f'<h3 style="margin: 0; font-size: 12px; font-weight: 600">Sequence preview</h3>', chip(f'{dot(ACC)} A · styrene × 20'), chip(f'{dot(SEL)} B · methyl methacrylate × 20'),
                    '<div style="flex-grow: 1"></div>', chip("Diblock A<sub>20</sub>-b-B<sub>20</sub>", TEXT, BG3), gap=8)
              + f'<div style="display: flex; gap: 3px">{"".join(cells)}</div>'
              + row(*[f'<span style="font-family: {MONO}; font-size: 10.5px; color: {DIM}; flex-grow: 1">{n}</span>' for n in ("1", "10", "20", "30", "40")], extra="justify-content: space-between")
              + row(chip("Tacticity along chain: r r r r r r r r · 100 % racemo", MUTED), chip("Linkage: head-to-tail", MUTED), chip("Caps: H / CH<sub>3</sub>", MUTED), gap=6)
              + '</div>')
    vp_ov = hud(chip("Preview · first 9 repeat units", TEXT, "#16191Ccc"), chip("syndiotactic", SEL, "#16191Ccc")) + seqbar

    ru = (f'<div style="display: flex; gap: 12px; padding: 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 8px">'
          f'<svg width="112" height="84" viewBox="0 0 112 84" role="img" aria-label="Styrene repeat unit sketch" style="flex-shrink: 0">'
          f'<path d="M10 58 L30 46 L50 58 L70 46" fill="none" stroke="{TEXT}" stroke-width="1.8"></path>'
          f'<text x="2" y="62" fill="{ACC}" font-size="13" font-family="IBM Plex Mono">*</text><text x="72" y="50" fill="{ACC}" font-size="13" font-family="IBM Plex Mono">*</text>'
          f'<path d="M30 46 L30 30" stroke="{TEXT}" stroke-width="1.8"></path>'
          f'<path d="M30 30 L44 22 L44 6 L30 -2 L16 6 L16 22 Z" transform="translate(0,6)" fill="none" stroke="{TEXT}" stroke-width="1.8"></path>'
          f'<circle cx="30" cy="20" r="6" fill="none" stroke="{MUTED}" stroke-width="1.2"></circle>'
          f'<text x="80" y="76" fill="{DIM}" font-size="10" font-family="IBM Plex Sans">sketch</text></svg>'
          f'<div style="display: flex; flex-direction: column; gap: 6px; min-width: 0">'
          f'<div style="font-weight: 600">Styrene</div>'
          f'<div style="font-family: {MONO}; font-size: 12px; color: {MUTED}">*CC(*)c1ccccc1</div>'
          f'<div style="display: flex; gap: 6px; flex-wrap: wrap">{chip("head C1", ACC, "#3A2C14", True)}{chip("tail C2", ACC, "#3A2C14", True)}{chip("104.15 g/mol", MUTED, BG2, True)}</div></div></div>')
    panel = (panel_head("Polymer builder", chip("Library: 124 polymers"))
             + '<div style="overflow: hidden; display: flex; flex-direction: column">'
             + section("Repeat units", col(ru, row(btn("From SMILES", ic="plus", small=True), btn("From library", ic="search", small=True), btn("Pick head/tail in 3D", ic="cursor", small=True), gap=6), gap=10), chip("A · B"))
             + section("Architecture", col(seg(["Linear", "Branched", "Star", "Comb", "Dendrimer", "Network"], "Linear", cols=3), gap=8))
             + section("Sequence", col(seg(["Homopolymer", "Alternating", "Block", "Random", "Gradient"], "Block", cols=3),
                                        row(field("Degree of polymerisation", "40"), field("Block lengths", "20 : 20"), gap=10), gap=10))
             + section("Stereo &amp; linkage", col(seg(["Isotactic", "Syndiotactic", "Atactic", "Custom"], "Syndiotactic", cols=4),
                                                   row(select("Linkage", "Head-to-tail"), select("Head cap", "H"), select("Tail cap", "CH<sub>3</sub>"), gap=8), gap=10))
             + '</div>'
             + footer(btn("Build in document", ic="cube"), btn("Send to CAPS Grow", True, "grow", href="Grow.dc.html")))
    return builder_shell("PS-b-PMMA_block.caps", "Polymer", svg, panel, vp_ov,
                         ('<span>Chain: 40 units · 2 blocks</span><span>Preview atoms: %d</span>' % len(atoms), '<span>RIS: Mattice &amp; Suter 1994</span>'))


# ---------------------------------------------------------------- Crystal
def crystal():
    a, b, c = 7.40, 4.93, 2.534
    atoms = pe_crystal(2, 3, 5, a, b, c)
    sc = Scene("c", VPW, VPH, yaw=0.62, pitch=0.5, persp=0.15, fog=0.5, atom_k=1.25)
    sc.add_atoms(atoms)
    sc.add_box((0, 0, 0), (2 * a, 0, 0), (0, 3 * b, 0), (0, 0, 5 * c / 2 * 1.0), MUTED, 1.0, "3 3", 0.7)
    sc.add_box((0, 0, 0), (a, 0, 0), (0, b, 0), (0, 0, c), ACC, 2.0, None, 1.0)
    svg = sc.svg()
    ov = hud(chip("Orthorhombic · Pnam", TEXT, "#16191Ccc"), chip(f'{dot(ACC)} unit cell', MUTED, "#16191Ccc"), chip("supercell 2 × 3 × 2", MUTED, "#16191Ccc"), chip("H hidden", MUTED, "#16191Ccc"))
    rows = [["C1", "C", "0.0380", "0.0650", "0.2500", "1.00"]]
    tbl = table(["Label", "El.", "x", "y", "z", "Occ."], rows, ["60px", "44px", "", "", "", "56px"], mono_cols=(2, 3, 4, 5), align_right=(2, 3, 4, 5))
    sg = []
    for name, num, on in (("Pnam", "62", True), ("Pnma", "62", False), ("P2₁/c", "14", False), ("Fm-3m", "225", False)):
        sg.append(f'<button style="display: flex; align-items: center; justify-content: space-between; height: 30px; padding: 0 10px; border-radius: 5px; border: 1px solid {ACC if on else LINE}; background: {BG3 if on else BG0}; color: {TEXT}; font-size: 12.5px; cursor: pointer"><span style="font-family: {MONO}">{name}</span><span style="color: {DIM}; font-family: {MONO}">No. {num}</span></button>')
    panel = (panel_head("Crystal builder", chip("ITA Vol. A · 230 groups"))
             + section("Space group", col(
                 f'<div style="display: flex; align-items: center; gap: 8px; height: 32px; padding: 0 8px; background: {BG0}; border: 1px solid {LINE}; border-radius: 5px; color: {TEXT}; font-size: 12.5px">{icon("search", 14, DIM)}<span style="font-family: {MONO}">Pnam</span></div>',
                 f'<div style="display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 6px">{"".join(sg)}</div>', gap=10), chip("Orthorhombic"))
             + section("Lattice", col(row(field("a", "7.40", "Å"), field("b", "4.93", "Å"), field("c", "2.534", "Å"), gap=8),
                                       row(field("α", "90", "°"), field("β", "90", "°"), field("γ", "90", "°"), gap=8),
                                       cite("Bunn, <i>Trans. Faraday Soc.</i> 35, 482 (1939) — orthorhombic polyethylene"), gap=10))
             + section("Asymmetric unit", col(tbl, row(btn("Add site", ic="plus", small=True), btn("Import CIF", ic="file", small=True), gap=6), gap=8), chip("1 site → 4 atoms / cell", MUTED))
             + section("Symmetry tools", col(row(btn("Apply symmetry", small=True), btn("Find symmetry", ic="search", small=True), btn("Primitive cell", small=True), gap=6),
                                              row(field("Supercell", "2 × 3 × 2"), select("Tolerance", "0.01 Å"), gap=10), gap=10))
             + footer(btn("Cancel"), btn("Build crystal", True, "cube")))
    return builder_shell("PE_crystal.caps", "Crystal", svg, panel, ov, (f"<span>{len(atoms)} C atoms (H hidden)</span><span>V cell 92.45 Å³</span>", "<span>Symmetry: spglib (Togo &amp; Tanaka)</span>"))


# ---------------------------------------------------------------- Surface / interface
def surface():
    slab = silica_slab(8, 8, 3)
    sc = Scene("s", VPW, VPH, yaw=0.55, pitch=0.62, persp=0.2, fog=0.5, outline=False)
    sc.add_atoms(slab)
    L = 8 * 2.6
    chains = []
    for k, seed in enumerate((11, 12, 13, 14, 15, 16)):
        st = (3 + (k % 3) * 6.5, 3 + (k // 3) * 10, 13.5 + (k % 2) * 3)
        pts = random_chain(st, 26, seed, box=((0.5, 0.5, 11.5), (L - 0.5, L - 0.5, 22.0)))
        chains.append([{"e": "C", "p": p, "c": "#B9BEC4"} for p in pts])
    for ch in chains:
        sc.add_atoms(ch, [(i, i + 1) for i in range(len(ch) - 1)])
    sc.add_box((-1.3, -1.3, -1.0), (L + 1.3, 0, 0), (0, L + 1.3, 0), (0, 0, 32.0), MUTED, 1.0, "3 3", 0.7)
    sc.add_box((-1.3, -1.3, 23.0), (L + 1.3, 0, 0), (0, L + 1.3, 0), (0, 0, 9.0), SEL, 1.0, "2 4", 0.6)
    sc._prep()
    p = sc.P((L + 1.5, -1.3, 27.5))
    ov_svg = label_pill(p[0] + 8, p[1], "vacuum 9.0 Å", SEL, border=SEL)
    p2 = sc.P((L + 1.5, -1.3, 17))
    ov_svg += label_pill(p2[0] + 8, p2[1], "PE film · 6 chains", TEXT)
    p3 = sc.P((L + 1.5, -1.3, 2))
    ov_svg += label_pill(p3[0] + 8, p3[1], "α-quartz (001) · 3 layers", "#D6A45E", border="#6B5634")
    svg = sc.svg(overlay=ov_svg)
    ov = hud(chip("Interface stack · 2 layers", TEXT, "#16191Ccc"), chip("slab + film + vacuum", MUTED, "#16191Ccc"))
    layers = [("3", "Vacuum", "9.0 Å", SEL), ("2", "PE film · amorphous (CAPS Grow)", "10.5 Å", "#B9BEC4"), ("1", "α-quartz (001) · O-terminated", "6.1 Å", "#D6A45E")]
    lr = "".join(f'<div style="display: flex; align-items: center; gap: 10px; height: 40px; padding: 0 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px">'
                 f'<span style="font-family: {MONO}; color: {DIM}; width: 12px">{n}</span><span style="width: 10px; height: 22px; border-radius: 2px; background: {c}"></span>'
                 f'<span style="flex-grow: 1; font-size: 12.5px">{t}</span><span style="font-family: {MONO}; font-size: 12px; color: {MUTED}">{h}</span>{icon("dots", 15, DIM)}</div>' for n, t, h, c in layers)
    panel = (panel_head("Surface &amp; interface builder")
             + section("Cleave", col(row(field("h", "0", w=60), field("k", "0", w=60), field("l", "1", w=60), select("Source", "α-quartz.cif"), gap=8),
                                      row(field("Slab thickness", "3", "layers"), field("Vacuum", "9.0", "Å"), gap=10),
                                      select("Termination", "O-terminated · 2 of 3 · non-polar"),
                                      check("Hydroxylate dangling bonds (Si–OH)", True), gap=10))
             + section("Layer stack", col(lr, row(btn("Add layer", ic="plus", small=True), btn("From CAPS Grow cell", ic="grow", small=True), gap=6), gap=6), chip("top → bottom"))
             + section("Lattice matching", col(
                 table(["Layer", "Supercell", "Strain a", "Strain b"],
                       [["quartz (001)", "4 × 4", "0.00 %", "0.00 %"], ["PE film", "fits cell", "—", "—"]],
                       mono_cols=(1, 2, 3), align_right=(2, 3)),
                 gap=10), chip("max strain 2.0 %")) 
             + footer(btn("Cancel"), btn("Build interface", True, "layers")))
    return builder_shell("quartz_PE_interface.caps", "Surface", svg, panel, ov, ("<span>Slab + film + vacuum</span>", "<span>Surface cleave from CIF</span>"))


# ---------------------------------------------------------------- Nanostructure
def nano():
    atoms, bonds, R = nanotube(10, 10, 24.0)
    sc = Scene("n", VPW, VPH, yaw=1.15, pitch=0.3, roll=1.2, persp=0.3, fog=0.62, atom_k=0.75, bond_w=0.16, outline=False)
    sc.add_atoms(atoms, bonds)
    svg = sc.svg()
    d = 2 * R
    ov = hud(chip(f"(10,10) armchair · d = {d:.2f} Å", TEXT, "#16191Ccc"), chip(f"{len(atoms)} C", MUTED, "#16191Ccc"), chip("periodic along z", SEL, "#16191Ccc"))
    kinds = [("Sheet", "hex"), ("Nanotube", "layers"), ("Particle", "atom"), ("Pore", "ring")]
    kb = "".join(f'<button style="height: 64px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 6px; border-radius: 6px; border: 1px solid {ACC if k == "Nanotube" else LINE}; background: {BG3 if k == "Nanotube" else BG0}; color: {TEXT if k == "Nanotube" else MUTED}; font-size: 12px; cursor: pointer">{icon(ic, 20, ACC if k == "Nanotube" else MUTED, 1.4)}<span>{k}</span></button>' for k, ic in kinds)
    panel = (panel_head("Nanostructure builder")
             + section("Type", f'<div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 6px">{kb}</div>')
             + section("Chirality", col(row(field("n", "10"), field("m", "10"), select("Material", "Carbon (graphene)"), gap=8),
                                         seg(["Armchair", "Zigzag", "Chiral"], "Armchair", True),
                                         kv("Chiral angle", "30.00°"), kv("Diameter", f"{d:.2f} Å"), kv("Translation |T|", "2.46 Å"),
                                         cite("d = a·√(n² + nm + m²)/π, a = 2.46 Å — Saito, Dresselhaus &amp; Dresselhaus (1998)"), gap=9))
             + section("Geometry", col(row(field("Length", "24.0", "Å"), field("C–C", "1.42", "Å"), gap=10),
                                        toggle("Periodic along tube axis", True), toggle("Cap ends with hydrogen", False), toggle("Multi-walled", False), gap=10))
             + footer(btn("Cancel"), btn("Build nanotube", True, "layers")))
    return builder_shell("CNT_10_10.caps", "Nanostructure", svg, panel, ov, (f"<span>{len(atoms)} atoms · {len(bonds)} bonds</span>", "<span>Bonds perceived (Cordero 2008 radii)</span>"))


# ---------------------------------------------------------------- Biomolecule
def bio():
    ca = helix_ca(22)
    # extend with a coil tail
    rnd = random.Random(4)
    tail = random_chain(add(ca[-1], (1.2, 1.0, 1.8)), 8, 9, bond=3.8, persistence=0.8)
    pts = ca + tail
    # smooth tube via catmull-rom
    sm = []
    for i in range(len(pts) - 1):
        p0 = pts[max(0, i - 1)]; p1 = pts[i]; p2 = pts[i + 1]; p3 = pts[min(len(pts) - 1, i + 2)]
        for s in range(6):
            t = s / 6
            t2, t3 = t * t, t * t * t
            sm.append(tuple(0.5 * ((2 * p1[k]) + (-p0[k] + p2[k]) * t + (2 * p0[k] - 5 * p1[k] + 4 * p2[k] - p3[k]) * t2 + (-p0[k] + 3 * p1[k] - 3 * p2[k] + p3[k]) * t3) for k in range(3)))
    sm.append(pts[-1])
    sc = Scene("b", VPW, VPH, yaw=0.3, pitch=0.2, roll=1.45, persp=0.25, fog=0.5, atom_k=1.0)
    nh = 22 * 6
    sc.add_tube(sm[:nh + 1], "#E07A5F", 0.95)
    sc.add_tube(sm[nh:], MUTED, 0.35)
    # CA atoms + CB stubs
    cas = [{"e": "C", "p": p, "r": 0.42, "c": "#E9ECEF"} for p in pts]
    sc.add_atoms(cas, [])
    for i, p in enumerate(ca):
        out = norm((p[0], p[1], 0))
        cb = add(p, mul(out, 1.53))
        side = [{"e": "C", "p": p}, {"e": "C", "p": cb}]
        if i % 3 == 0:
            side.append({"e": "O", "p": add(cb, mul(out, 1.43))})
        elif i % 4 == 1:
            side.append({"e": "N", "p": add(cb, mul(norm(add(out, (0, 0, 0.6))), 1.47))})
        sc.add_atoms(side, [(k, k + 1) for k in range(len(side) - 1)])
    svg = sc.svg()
    ov = hud(chip("Peptide · 30 residues", TEXT, "#16191Ccc"), chip(f'{dot("#E07A5F")} α-helix 1–22', MUTED, "#16191Ccc"), chip(f'{dot(MUTED)} coil 23–30', MUTED, "#16191Ccc"))
    seq = "AEAAAKEAAAKEAAAKAGGSPGSG"[:22] + "GSPGSGSG"
    ss = ["H"] * 22 + ["C"] * 8
    cells = []
    for i, (r, s) in enumerate(zip(seq, ss)):
        c = "#E07A5F" if s == "H" else BG3
        tc = "#1A0C08" if s == "H" else MUTED
        cells.append(f'<span style="height: 26px; display: flex; align-items: center; justify-content: center; border-radius: 3px; background: {c}; color: {tc}; font-family: {MONO}; font-size: 11.5px; font-weight: 500">{r}</span>')
    grid = f'<div style="display: grid; grid-template-columns: repeat(10, minmax(0, 1fr)); gap: 3px">{"".join(cells)}</div>'
    panel = (panel_head("Biomolecule builder")
             + tabs(["Peptide", "Nucleic acid"], "Peptide")
             + section("Sequence", col(
                 f'<div style="padding: 8px 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 5px; font-family: {MONO}; font-size: 12.5px; letter-spacing: 0.06em; word-break: break-all">{seq}</div>',
                 row(chip("30 residues", MUTED), chip("one-letter codes", MUTED), chip("Import FASTA", ACC), gap=6), gap=8))
             + section("Secondary structure", col(grid, seg(["α-helix", "β-strand", "PPII", "Coil"], "α-helix", cols=4),
                                                   row(field("φ", "−57", "°"), field("ψ", "−47", "°"), field("ω", "180", "°"), gap=8),
                                                   cite("Ideal α-helix φ/ψ — Pauling, Corey &amp; Branson, <i>PNAS</i> 37, 205 (1951)"), gap=10))
             + section("Termini &amp; state", col(row(select("N-terminus", "NH<sub>3</sub><sup>+</sup>"), select("C-terminus", "COO<sup>−</sup>"), gap=10),
                                                   row(select("Protonation", "pH 7.0 defaults"), select("Force field", "AMBER ff14SB"), gap=10), gap=10))
             + footer(btn("Cancel"), btn("Build peptide", True, "cube")))
    return builder_shell("peptide_EAAAK.caps", "Biomolecule", svg, panel, ov, ("<span>Ribbon + side-chain stubs</span>", "<span>Style: ribbon · CA spheres</span>"))


# ---------------------------------------------------------------- Solvation
def solvation():
    Lb = 30.0
    atoms, backbone, stereo = polystyrene(4, "atactic", seed=7, curve=False)
    cx = sum(a["p"][0] for a in atoms) / len(atoms)
    cy = sum(a["p"][1] for a in atoms) / len(atoms)
    cz = sum(a["p"][2] for a in atoms) / len(atoms)
    atoms = [{**a, "p": (a["p"][0] - cx + Lb / 2, a["p"][1] - cy + Lb / 2, a["p"][2] - cz + Lb / 2)} for a in atoms]
    sc = Scene("v", VPW, VPH, yaw=0.62, pitch=0.42, persp=0.25, fog=0.62, atom_k=1.0)
    sc.add_atoms(atoms)
    rnd = random.Random(21)
    placed = [a["p"] for a in atoms]
    nw = 0
    tries = 0
    while nw < 150 and tries < 20000:
        tries += 1
        o = (rnd.uniform(1, Lb - 1), rnd.uniform(1, Lb - 1), rnd.uniform(1, Lb - 1))
        if all(dist(o, q) > 3.0 for q in placed):
            w, b = water(o, rnd)
            sc.add_atoms(w, b)
            placed.append(o)
            nw += 1
    ions = []
    for e in ("Na", "Na", "Cl", "Cl"):
        while True:
            o = (rnd.uniform(3, Lb - 3), rnd.uniform(3, Lb - 3), rnd.uniform(3, Lb - 3))
            if all(dist(o, q) > 3.6 for q in placed):
                placed.append(o); ions.append({"e": e, "p": o}); break
    sc.add_atoms(ions, [])
    sc.add_box((0, 0, 0), (Lb, 0, 0), (0, Lb, 0), (0, 0, Lb), ACC, 1.6, None, 0.9)
    svg = sc.svg()
    ov = hud(chip("Preview · sparse sample of solvent shown", TEXT, "#16191Ccc"), chip(f'{dot("#9B7BD6")} Na⁺', MUTED, "#16191Ccc"), chip(f'{dot("#57B26A")} Cl⁻', MUTED, "#16191Ccc"))
    stages = [("Random sequential insertion", "done", OK), ("Overlap minimisation · L-BFGS", "done", OK), ("Soft-core compression MD", "running", ACC), ("Verify min. distance ≥ 2.0 Å", "pending", DIM)]
    st = "".join(f'<div style="display: flex; align-items: center; gap: 10px; font-size: 12.5px">{dot(c, 9)}<span style="flex-grow: 1">{n}</span><span style="font-family: {MONO}; font-size: 11.5px; color: {c}">{s}</span></div>' for n, s, c in stages)
    panel = (panel_head("Solvation &amp; ions", chip("uses CAPS Pack"))
             + section("Box", col(row(select("Shape", "Cubic"), field("Edge", "30.0", "Å"), field("Min. distance", "2.0", "Å"), gap=8), gap=10))
             + section("Solvent", col(row(select("Solvent", "Water"), select("Model", "TIP4P/2005"), gap=10),
                                       row(field("Density", "0.997", "g/cm³"), field("Molecules", "auto"), gap=10), gap=10))
             + section("Ions", col(seg(["Neutralise only", "Concentration", "Custom counts"], "Concentration", True),
                                    row(select("Salt", "NaCl"), field("Concentration", "0.15", "mol/L"), gap=10),
                                    kv("Counts (from free volume)", "computed on build", False, MUTED), gap=10))
             + section("Packing progress", col(st, f'<div style="height: 6px; border-radius: 3px; background: {BG3}"><div style="width: 68%; height: 6px; border-radius: 3px; background: {ACC}"></div></div>', gap=9), chip("stage 3 of 4", ACC, "#3A2C14"))
             + footer(btn("Cancel"), btn("Solvate", True, "flask")))
    return builder_shell("PS_in_water.caps", "Solvation", svg, panel, ov, ("<span>Solute + solvent + ions</span>", "<span>Packing: CAPS Pack (RSA → L-BFGS → compression)</span>"))
