import math, random
from lib import *
from mols import *
from screen_app import card, pagehead, shell, progress
from screen_builders import builder_shell, footer, hud, VPW, VPH
from screen_studio2 import spacer

ELEM.setdefault("Au", ("#E0B84C", 0.6))


def analyze_side(active):
    calcs = [("Structure", ["Density", "RDF", "S(q)", "X-ray / neutron"]), ("Chains", ["R<sub>g</sub>", "C<sub>∞</sub>", "Entanglements"]),
             ("Thermo", ["CED", "δ", "T<sub>g</sub>"]), ("Mechanics", ["C<sub>ij</sub> strain", "C<sub>ij</sub> fluct.", "Stress–strain"]),
             ("Dynamics", ["MSD", "D", "Relaxation"]), ("Free volume", ["Probe insertion", "Pore size"])]
    cl = []
    for grp, items in calcs:
        its = "".join(f'<button style="height: 26px; padding: 0 9px; border-radius: 5px; background: {BG3 if it in active else BG0}; border: 1px solid {ACC if it in active else LINE}; font-size: 12px; color: {TEXT}; cursor: pointer">{it}</button>' for it in items)
        cl.append(f'<div style="display: flex; flex-direction: column; gap: 6px"><span style="font-size: 11px; font-weight: 600; letter-spacing: 0.06em; text-transform: uppercase; color: {DIM}">{grp}</span><div style="display: flex; flex-wrap: wrap; gap: 5px">{its}</div></div>')
    return (f'<aside aria-label="Calculations" style="width: 280px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-right: 1px solid {LINE}">'
            + panel_head("Calculations") + f'<div style="padding: 14px; display: flex; flex-direction: column; gap: 16px">{"".join(cl)}</div>'
            + section("Source", col(select("Trajectory", "equil-7 · last 5 ns"), select("Frames", "every 10 ps · 500 frames"), gap=10)) + '</aside>')


def result_cells(items):
    return "".join(f'<div style="display: flex; flex-direction: column; gap: 5px; padding: 12px 14px; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px">'
                   f'<span style="font-size: 12px; color: {MUTED}">{a}</span><span style="font-family: {MONO}; font-size: 20px; color: {DIM}">{b}</span><span style="font-size: 11px; color: {DIM}">{c}</span></div>' for a, b, c in items)


# ------------------------------------------------------------------ Mechanics
def mechanics():
    eps = [i * 0.002 for i in range(0, 101)]
    E0 = 3.0
    ss = [(e * 100, (E0 * e if e < 0.03 else E0 * 0.03 + 0.03 * (1 - math.exp(-(e - 0.03) / 0.03))) + 0.002 * math.sin(e * 180)) for e in eps]
    ss = [(x, max(0, y)) for x, y in ss]
    fit = [(0, 0), (2.0, E0 * 0.02)]
    p = plot(640, 300, [(ss, ACC, 1.8, None), (fit, SEL, 2, "6 4")], (0, 20), (0, 0.14), [0, 5, 10, 15, 20], [0, 0.04, 0.08, 0.12], "engineering strain (%)", "true stress (GPa)", pad=(48, 14, 20, 32),
             bands=[(0, 2, SEL, 0.08)])
    p = p.replace('true stress', 'true stress')
    Cij = [["C<sub>11</sub>", "[result]", "C<sub>12</sub>", "[result]", "C<sub>13</sub>", "[result]"], ["", "", "C<sub>22</sub>", "[result]", "C<sub>23</sub>", "[result]"],
           ["", "", "", "", "C<sub>33</sub>", "[result]"], ["C<sub>44</sub>", "[result]", "C<sub>55</sub>", "[result]", "C<sub>66</sub>", "[result]"]]
    ct = table(["", "GPa", "", "GPa", "", "GPa"], Cij, ["12%", "21%", "12%", "21%", "13%", "21%"], mono_cols=(1, 3, 5), align_right=(1, 3, 5), fs=12)
    res = result_cells([("Young's modulus E", "[result]", "fit 0–2 % strain"), ("Poisson ratio ν", "[result]", "from lateral strain"),
                        ("Bulk modulus K", "[result]", "(C₁₁ + 2C₁₂)/3, Voigt"), ("Yield stress", "[result]", "0.2 % offset")])
    content = (pagehead("Analyze · mechanics", "Elastic constants by constant strain, and uniaxial stress–strain at a stated rate",
                        row(btn("Export CSV / LaTeX", ic="download"), btn("Run", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 18px 22px; min-height: 0">'
               + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px">'
               + f'<div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 10px">{res}</div>'
               + card("Uniaxial tension · z", col(p, row(chip(f'{dot(ACC)} stress–strain · illustrative'), chip(f'{dot(SEL)} linear fit window'), gap=6), gap=8), chip("3 directions averaged"), 12)
               + '</div>'
               + f'<div style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
               + card("Stiffness matrix", col(ct, f'<div style="font-size: 11.5px; color: {MUTED}; line-height: 1.5">± 0.5 % strain per component, stress from the virial, averaged over 3 relaxed replicas.</div>', gap=8))
               + card("Protocol", col(row(select("Method", "Constant strain"), field("Strain", "±0.5", "%"), gap=10),
                                      row(field("Tension rate", "1e8", "s⁻¹"), select("Lateral", "NPT · 1 atm"), gap=10),
                                      cite("Theodorou &amp; Suter, <i>Macromolecules</i> 19, 139 (1986) — static elastic constants of glassy polymers"), gap=10))
               + '</div></div>')
    return shell("Analyze", "CAPS Analyze — mechanics", f'{icon("chevr", 12, DIM)}<span>Analyze</span>', content,
                 ("<span>PS cell · 12 840 atoms · 3 replicas</span>", "<span>MD strain rates are orders of magnitude above experiment; reported with every value</span>"))


# ------------------------------------------------------------------ Scattering
def scattering():
    q = [i * 0.01 for i in range(20, 400)]

    def halo(x, c, w, a):
        return a * math.exp(-((x - c) ** 2) / (2 * w * w))
    xr = [(x, 0.35 + halo(x, 0.75, 0.12, 0.9) + halo(x, 1.4, 0.18, 1.6) + halo(x, 3.0, 0.35, 0.35)) for x in q]
    ns = [(x, 0.5 + halo(x, 0.75, 0.12, 0.5) + halo(x, 1.35, 0.2, 0.8) + halo(x, 2.9, 0.4, 0.5)) for x in q]
    p = plot(700, 330, [(xr, ACC, 2, None), (ns, SEL, 2, "6 4")], (0.2, 4), (0, 2.4), [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4], [0, 0.5, 1, 1.5, 2], "q (Å⁻¹)", "I(q) · normalised", pad=(44, 14, 20, 32))
    b = table(["Nucleus", "b<sub>coh</sub> (fm)"], [["¹H", "−3.739"], ["²H (D)", "6.671"], ["C", "6.646"]], ["55%", "45%"], mono_cols=(1,), align_right=(1,), fs=12)
    content = (pagehead("Analyze · scattering", "X-ray and neutron patterns from the simulated cell, ready to overlay on measured data",
                        row(btn("Load experimental curve", ic="download"), btn("Run", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 18px 22px; min-height: 0">'
               + card("I(q) · PS amorphous cell", col(p, row(chip(f'{dot(ACC)} X-ray (Cu Kα)'), chip(f'{dot(SEL)} neutron, deuterated backbone'), spacer(), chip("illustrative shape"), gap=6),
                                                      f'<div style="font-size: 12px; color: {MUTED}; line-height: 1.5">Glassy PS shows two broad maxima near q ≈ 0.75 Å⁻¹ (the \'polymerisation\' peak from phenyl packing) and q ≈ 1.4 Å⁻¹ (the amorphous halo). The simulated curve is compared at these positions.</div>', gap=10), extra="flex-grow: 1")
               + f'<div style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
               + card("Method", col(row(select("From", "Debye sum over pairs"), field("r<sub>max</sub>", "L/2"), gap=10), row(field("q range", "0.2 – 4.0", "Å⁻¹"), field("Δq", "0.01", "Å⁻¹"), gap=10), gap=10))
               + card("X-ray", col(select("Form factors", "Cromer–Mann, International Tables Vol. C"), select("Source", "Cu Kα · λ 1.5406 Å"), gap=10))
               + card("Neutron", col(b, select("Isotope pattern", "d-backbone, h-phenyl"), cite("Scattering lengths: Sears, <i>Neutron News</i> 3, 26 (1992)"), gap=10))
               + '</div></div>')
    return shell("Analyze", "CAPS Analyze — scattering", f'{icon("chevr", 12, DIM)}<span>Analyze</span>', content,
                 ("<span>500 frames · 12 840 atoms</span>", "<span>GPU pair sum · deterministic reduction</span>"))


# ------------------------------------------------------------------ Free volume
def free_volume():
    Lb = 30.0
    sc = Scene("fv", 700, 560, yaw=0.6, pitch=0.42, persp=0.3, fog=0.55)
    pts_all = []
    cols = ["#8E959C"] * 10
    for k in range(10):
        rnd = random.Random(700 + k)
        st = (rnd.uniform(3, Lb - 3), rnd.uniform(3, Lb - 3), rnd.uniform(3, Lb - 3))
        pts = random_chain(st, 70, 800 + k, box=((0.8, 0.8, 0.8), (Lb - 0.8,) * 3), persistence=0.55)
        pts_all += pts
        sc.add_tube(pts, mix("#8E959C", BG0, 0.25), 0.55)
    sc.add_box((0, 0, 0), (Lb, 0, 0), (0, Lb, 0), (0, 0, Lb), MUTED, 1.2, None, 0.8)
    rnd = random.Random(3)
    voids = []
    for _ in range(6000):
        p = (rnd.uniform(1.5, Lb - 1.5), rnd.uniform(1.5, Lb - 1.5), rnd.uniform(1.5, Lb - 1.5))
        d = min(dist(p, q) for q in pts_all)
        r = d - 1.7
        if r > 1.2 and all(dist(p, v) > v_r + r for v, v_r in voids):
            voids.append((p, r))
        if len(voids) >= 40:
            break
    sc._prep()
    ov = ""
    for p, r in sorted(voids, key=lambda v: sc.P(v[0])[2]):
        x, y, z, k = sc.P(p)
        ov += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r * sc.scale * k:.1f}" fill="{SEL}" fill-opacity="0.22" stroke="{SEL}" stroke-opacity="0.7" stroke-width="1"></circle>'
    view = sc.svg(overlay=ov)
    radii = [r for _, r in voids]
    hist = [0] * 8
    for r in radii:
        hist[min(7, int((r - 1.0) / 0.25))] += 1
    hp = [(1.0 + i * 0.25 + 0.125, h) for i, h in enumerate(hist)]
    bars = plot(380, 200, [], (1.0, 3.0), (0, max(hist) + 2), [1.0, 1.5, 2.0, 2.5, 3.0], list(range(0, max(hist) + 3, 4)), "largest-sphere radius (Å)", "voids (count)", pad=(40, 12, 18, 30))
    # draw bars into the plot svg
    L_, R_, T_, B_ = 40, 12, 18, 30
    pw, ph = 380 - L_ - R_, 200 - T_ - B_
    bx = ""
    for c, h in hp:
        x0 = L_ + (c - 0.125 - 1.0) / 2.0 * pw + 2
        w = 0.25 / 2.0 * pw - 4
        hh = h / (max(hist) + 2) * ph
        bx += f'<rect x="{x0:.1f}" y="{T_ + ph - hh:.1f}" width="{w:.1f}" height="{hh:.1f}" rx="2" fill="{SEL}"></rect>'
    bars = bars.replace("</svg>", bx + "</svg>")
    content = (pagehead("Analyze · free volume", "Probe insertion on a grid, with the distribution of void sizes",
                        row(btn("Export voids as PDB", ic="download"), btn("Run", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 18px 22px; min-height: 0">'
               + card("Voids in the cell", f'<div style="position: relative; margin: -12px; background: {BG0}">{view}{hud(chip(f"{len(voids)} largest voids shown", SEL, "#16191Ccc"), chip("chains as tubes", MUTED, "#16191Ccc"))}</div>', chip("preview"), 12, "flex-grow: 1")
               + f'<div style="width: 420px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
               + f'<div style="display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 10px">{result_cells([("Fractional free volume", "[result]", "probe 0 Å, Bondi radii"), ("Accessible to probe", "[result]", "probe 1.4 Å")])}</div>'
               + card("Void size distribution · preview", col(bars, f'<div style="font-size: 11.5px; color: {DIM}">Counts from the {len(voids)} voids drawn in this preview, not a converged analysis.</div>', gap=6), None or "", 12)
               + card("Method", col(row(field("Probe radius", "1.4", "Å"), field("Grid", "0.2", "Å"), gap=10), select("Atomic radii", "Bondi van der Waals"),
                                    cite("Bondi, <i>J. Phys. Chem.</i> 68, 441 (1964)"), gap=10))
               + '</div></div>')
    return shell("Analyze", "CAPS Analyze — free volume", f'{icon("chevr", 12, DIM)}<span>Analyze</span>', content,
                 ("<span>Grid 150³ points · 0.2 Å</span>", "<span>Void geometry computed from the drawn preview cell</span>"))


# ------------------------------------------------------------------ Nanoparticle
def nanoparticle():
    a = 4.078  # Au fcc lattice constant, Angstrom
    R = 12.0
    basis = [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)]
    n = int(R / a) + 2
    atoms = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                for bx, by, bz in basis:
                    p = ((i + bx) * a, (j + by) * a, (k + bz) * a)
                    if math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2) <= R:
                        atoms.append({"e": "Au", "p": p, "r": 1.30})
    sc = Scene("np", VPW, VPH, yaw=0.5, pitch=0.35, persp=0.25, fog=0.55, atom_k=1.0, outline=True)
    sc.add_atoms(atoms, [])
    svg = sc.svg()
    N = len(atoms)
    shapes = [("Sphere", True), ("Cuboctahedron", False), ("Truncated octahedron", False), ("Icosahedron", False)]
    sg = "".join(f'<button style="height: 34px; border-radius: 6px; border: 1px solid {ACC if on else LINE}; background: {BG3 if on else BG0}; font-size: 12px; color: {TEXT if on else MUTED}; cursor: pointer">{s}</button>' for s, on in shapes)
    kinds = [("Sheet", "hex"), ("Nanotube", "layers"), ("Particle", "atom"), ("Pore", "ring")]
    kb = "".join(f'<button style="height: 64px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 6px; border-radius: 6px; border: 1px solid {ACC if k == "Particle" else LINE}; background: {BG3 if k == "Particle" else BG0}; color: {TEXT if k == "Particle" else MUTED}; font-size: 12px; cursor: pointer">{icon(ic, 20, ACC if k == "Particle" else MUTED, 1.4)}<span>{k}</span></button>' for k, ic in kinds)
    panel = (panel_head("Nanostructure builder")
             + section("Type", f'<div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 6px">{kb}</div>')
             + section("Lattice", col(row(select("Material", "Au · fcc"), field("a", "4.078", "Å"), gap=10), cite("Lattice constant at room temperature: Davey, <i>Phys. Rev.</i> 25, 753 (1925)"), gap=10))
             + section("Shape", col(f'<div style="display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 6px">{sg}</div>', row(field("Radius", "12.0", "Å"), select("Centre", "On an atom"), gap=10), kv("Atoms", f"{N}"), gap=10))
             + section("Surface", col(toggle("Cap with thiolate ligands", False), toggle("Relax surface with force field", True), gap=9))
             + footer(btn("Cancel"), btn("Build particle", True, "atom")))
    ov = hud(chip(f"Au sphere · r = 12.0 Å · {N} atoms", TEXT, "#16191Ccc"), chip("space filling", MUTED, "#16191Ccc"))
    return builder_shell("Au_np_12A.caps", "Nanostructure", svg, panel, ov, (f"<span>{N} Au atoms</span>", "<span>Cut from an fcc lattice; count computed</span>"))


# ------------------------------------------------------------------ Slit pore
def pore():
    acc = 1.42
    a1 = (acc * math.sqrt(3), 0.0)
    a2 = (acc * math.sqrt(3) / 2, acc * 1.5)
    base = [(0.0, 0.0), (0.0, acc)]
    Lx, Ly = 26.0, 22.0
    H = 10.0  # carbon-centre to carbon-centre width
    sheet = []
    for i in range(-20, 30):
        for j in range(-5, 20):
            for bx, by in base:
                x = i * a1[0] + j * a2[0] + bx
                y = i * a1[1] + j * a2[1] + by
                if 0 <= x < Lx and 0 <= y < Ly:
                    sheet.append((x, y))
    sc = Scene("sp", VPW, VPH, yaw=0.5, pitch=0.3, persp=0.3, fog=0.55, outline=False, bond_w=0.16, atom_k=0.7)
    for z in (0.0, H):
        sc.add_atoms([{"e": "C", "p": (x, y, z)} for x, y in sheet])
    rnd = random.Random(9)
    fl = []
    while len(fl) < 16:
        p = (rnd.uniform(2, Lx - 2), rnd.uniform(2, Ly - 2), rnd.uniform(3.4, H - 3.4))
        if all(dist(p, q) > 3.8 for q in fl):
            fl.append(p)
    sc.add_atoms([{"e": "C", "p": p, "r": 1.86, "c": SEL} for p in fl], [])
    sc._prep()
    pa, pb = sc.P((Lx + 1.5, 0, 0)), sc.P((Lx + 1.5, 0, H))
    ov = (f'<line x1="{pa[0]:.1f}" y1="{pa[1]:.1f}" x2="{pb[0]:.1f}" y2="{pb[1]:.1f}" stroke="{ACC}" stroke-width="1.6"></line>'
          + label_pill(pa[0] + 10, (pa[1] + pb[1]) / 2, f"H = {H:.1f} Å", ACC, border=ACC))
    svg = sc.svg(overlay=ov)
    kinds = [("Sheet", "hex"), ("Nanotube", "layers"), ("Particle", "atom"), ("Pore", "ring")]
    kb = "".join(f'<button style="height: 64px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 6px; border-radius: 6px; border: 1px solid {ACC if k == "Pore" else LINE}; background: {BG3 if k == "Pore" else BG0}; color: {TEXT if k == "Pore" else MUTED}; font-size: 12px; cursor: pointer">{icon(ic, 20, ACC if k == "Pore" else MUTED, 1.4)}<span>{k}</span></button>' for k, ic in kinds)
    panel = (panel_head("Nanostructure builder")
             + section("Type", f'<div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 6px">{kb}</div>')
             + section("Pore", col(seg(["Slit", "Cylinder", "From CIF (zeolite, MOF)"], "Slit", True), row(field("Width H", "10.0", "Å", note="carbon centre to centre"), select("Walls", "1 graphene layer each"), gap=10), gap=10))
             + section("Periodic cell", col(row(field("x", f"{Lx:.1f}", "Å"), field("y", f"{Ly:.1f}", "Å"), field("z", "30.0", "Å"), gap=8), toggle("Vacuum above walls (non-periodic in z)", False), gap=10))
             + section("Fluid", col(row(select("Molecule", "Methane · TraPPE-UA"), field("Count", "16"), gap=10), select("Fill with", "CAPS Pack · inside pore region only"), gap=10))
             + footer(btn("Cancel"), btn("Build pore", True, "ring")))
    ov2 = hud(chip(f"slit pore · 2 × {len(sheet)} C", TEXT, "#16191Ccc"), chip(f'{dot(SEL)} 16 CH₄ (united atom)', MUTED, "#16191Ccc"))
    return builder_shell("graphite_slit_10A.caps", "Nanostructure", svg, panel, ov2, (f"<span>{2 * len(sheet) + 16} particles</span>", "<span>Walls fixed · fluid packed by CAPS Pack</span>"))


# ------------------------------------------------------------------ Failed job
RADIO = f'<span style="width: 8px; height: 8px; border-radius: 50%; background: {ACC}"></span>'


def failed_job():
    jl = [("equil-7", "Equilibrate · 21-step", "done", OK), ("dyn-12", "Dynamics · NPT 10 ns", "running", ACC), ("react-2", "React · epoxy 85 %", "failed", ERR), ("pack-5", "Pack · PS in water", "done", OK)]
    rows = "".join(f'<button style="display: flex; flex-direction: column; gap: 6px; padding: 10px 12px; text-align: left; background: {BG3 if j == "react-2" else "transparent"}; border: 1px solid {LINE if j == "react-2" else "transparent"}; border-radius: 7px; cursor: pointer">'
                   f'<span style="display: flex; align-items: center; gap: 8px; width: 100%">{dot(c)}<span style="font-size: 12.5px; font-weight: 500; color: {TEXT}; flex-grow: 1">{n}</span><span style="font-family: {MONO}; font-size: 11px; color: {c}">{s}</span></span>'
                   f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{j}</span></button>' for j, n, s, c in jl)
    joblist = (f'<aside aria-label="Jobs" style="width: 290px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-right: 1px solid {LINE}">'
               + panel_head("Jobs", btn("New", ic="plus", small=True)) + f'<div style="padding: 8px; display: flex; flex-direction: column; gap: 4px">{rows}</div></aside>')
    f = [(i, 20 + 8 * math.sin(i / 3) + (0 if i < 40 else (i - 40) ** 2.4 * 6)) for i in range(0, 46)]
    fp = plot(520, 200, [(f, ERR, 1.8, None), ([(0, 400), (46, 400)], WARN, 1.2, "5 4")], (0, 46), (0, 1600), [0, 10, 20, 30, 40], [0, 400, 800, 1200, 1600], "reaction cycle", "max |F| (kcal/mol/Å)", pad=(48, 12, 18, 30))
    steps = [("check", "Cycle 41 formed bond N(amine)–C(epoxide) between atoms 18 204 and 18 377", OK),
             ("alert", "Relaxation after the bond ran 5 ps; max force stayed above the 400 kcal/mol/Å cap", WARN),
             ("xcircle", "Step 184 220: atom 18 377 moved 9.4 Å in one step (more than half the cutoff); run stopped", ERR)]
    st = "".join(f'<div style="display: flex; gap: 10px; align-items: flex-start; font-size: 12.5px; line-height: 1.45">{icon(i, 16, c)}<span>{t}</span></div>' for i, t, c in steps)
    fixes = [("Restart from checkpoint 40 with a longer relaxation after each bond (5 → 20 ps)", True),
             ("Also cap the capture distance at 4.0 Å (now 4.5 Å)", False),
             ("Open the structure at cycle 41 in Studio to inspect the new bond", False)]
    fx = "".join(f'<div style="display: flex; gap: 10px; align-items: flex-start; padding: 10px 12px; background: {BG3 if on else BG0}; border: 1px solid {ACC if on else LINE}; border-radius: 7px; font-size: 12.5px">'
                 f'<span style="width: 16px; height: 16px; margin-top: 1px; border-radius: 50%; border: 2px solid {ACC if on else DIM}; flex-shrink: 0; display: flex; align-items: center; justify-content: center">{RADIO if on else ""}</span><span>{t}</span></div>' for t, on in fixes)
    main = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
            + pagehead("React · epoxy 85 % failed at cycle 41", "react-2 · DGEBA/IPDA · 38 % conversion reached · checkpoint 40 is intact",
                       row(chip(f'{dot(ERR)} failed', ERR, "#2A1414"), btn("Copy diagnostics", ic="copy"), btn("Restart with fix", True, "play", href="Jobs.dc.html"), gap=8))
            + '<div style="flex-grow: 1; display: flex; gap: 14px; padding: 18px 22px; min-height: 0">'
            + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px">'
            + card("What happened", col(st, gap=10), chip("from the structured log"))
            + card("Force before the failure", fp, chip(f'{dot(WARN)} force cap'), 12)
            + '</div>'
            + f'<div style="width: 460px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
            + card("Suggested recovery", col(fx, f'<div style="font-size: 12px; color: {MUTED}; line-height: 1.5">A newly formed bond starts far from its equilibrium length, so the relaxation after it needs enough time to bring forces down before normal dynamics resume.</div>', gap=8))
            + card("Nothing is lost", col(kv("Last checkpoint", "cycle 40 · 37.6 % conversion"), kv("Provenance", "failure recorded, reason and step"), kv("Other jobs", "unaffected"), gap=8))
            + '</div></div></div>')
    content = f'<div style="flex-grow: 1; display: flex; min-height: 0">{joblist}{main}</div>'
    return shell("Jobs", "CAPS Jobs — failed run", f'{icon("chevr", 12, DIM)}<span>Jobs</span>', content,
                 ("<span>1 failed · 1 running</span>", "<span>Failure reasons are written to the manifest</span>"))
