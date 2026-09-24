import math, random
from lib import *
from mols import *
from screen_app import shell, pagehead, card
from screen_studio2 import studio_shell, bar, spacer, styrene3d
from screen_builders import footer, hud
from screen_row16b import jacobi3

NA = 6.02214076e23
VW = 1440 - 72 - 380
VH = 900 - 44 - 48 - 26


def cell_matrix(a, b, c, al, be, ga):
    ca, cb, cg = (math.cos(math.radians(x)) for x in (al, be, ga))
    sg = math.sin(math.radians(ga))
    va = (a, 0.0, 0.0)
    vb = (b * cg, b * sg, 0.0)
    cx = c * cb
    cy = c * (ca - cb * cg) / sg
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return va, vb, (cx, cy, cz)


def volume(a, b, c, al, be, ga):
    ca, cb, cg = (math.cos(math.radians(x)) for x in (al, be, ga))
    return a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)


CELLS = [("Polyethylene · orthorhombic", (7.40, 4.93, 2.534, 90, 90, 90), 4, 14.027, "CH₂", "≈ 1.00", "Bunn 1939"),
         ("α-quartz · hexagonal", (4.913, 4.913, 5.405, 90, 90, 120), 3, 60.084, "SiO₂", "2.65", ""),
         ("PET · triclinic", (4.56, 5.94, 10.75, 98.5, 118.0, 112.0), 1, 192.17, "C₁₀H₈O₄", "1.455", "Daubeny, Bunn &amp; Brown 1954")]


# ---------------------------------------------------------------- Cell editor
def cell_editor():
    name, p, Z, M, unit, lit, src = CELLS[2]
    va, vb, vc = cell_matrix(*p)
    V = volume(*p)
    rho = Z * M / NA / (V * 1e-24)
    sc = Scene("ce", VW, VH, yaw=0.55, pitch=0.35, persp=0.25, fog=0.3)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                o = add(add(mul(va, i), mul(vb, j)), mul(vc, k))
                on = (i, j, k) == (0, 0, 0)
                sc.add_box(o, va, vb, vc, ACC if on else MUTED, 2.0 if on else 1.0, None if on else "4 4", 0.95 if on else 0.5)
    sc._prep()
    lbl = ""
    for v, t, col_ in ((va, "a", "#E35049"), (vb, "b", "#7DC884"), (vc, "c", "#2271DB")):
        sc.add_line((0, 0, 0), v, col_, 3.2, None, 1)
    sc.center = None
    sc._prep()
    for v, t in ((va, "a"), (vb, "b"), (vc, "c")):
        q = sc.P(mul(v, 1.06))
        lbl += label_pill(q[0] + 4, q[1] - 10, t, TEXT, border=LINE)
    view = sc.svg(overlay=lbl)
    rows = []
    for n, pp, z, m, u, l, s in CELLS:
        v = volume(*pp)
        rows.append([n.split(" · ")[0].replace("Polyethylene", "PE"), f"{v:.2f}", str(z), f"{z * m / NA / (v * 1e-24):.3f}", l])
    ct = table(["Cell", "V (Å³)", "Z", "ρ calc", "ρ lit."], rows, ["28%", "22%", "10%", "20%", "20%"], mono_cols=(1, 2, 3, 4), align_right=(1, 3, 4), fs=11, rowh=28, hl={2})
    mat = "".join(f'<div style="display: flex; gap: 10px"><span style="width: 14px; color: {c_}">{t}</span>' + "".join(f'<span style="width: 70px; text-align: right">{x:8.3f}</span>' for x in v) + '</div>'
                  for t, v, c_ in (("a", va, "#E35049"), ("b", vb, "#7DC884"), ("c", vc, "#2271DB")))
    right = (panel_head("Unit cell", chip("triclinic · P1̄", MUTED, BG2))
             + section("Parameters", col(row(field("a", f"{p[0]}", "Å"), field("b", f"{p[1]}", "Å"), field("c", f"{p[2]}", "Å"), gap=6),
                                         row(field("α", f"{p[3]}", "°"), field("β", f"{p[4]}", "°"), field("γ", f"{p[5]}", "°"), gap=6),
                                         row(field("Z", str(Z)), field("Formula unit", unit, mono=False), gap=6), gap=8), chip("PET", MUTED, BG2))
             + section("Derived", col(kv("Volume", f"{V:.2f} Å³"), kv("Density", f"{rho:.4f} g/cm³", vcol=OK), kv("Literature", f"{lit} g/cm³ · {src}"), gap=6))
             + section("Cell vectors (Å, a along x, b in xy)", f'<div style="font-family: {MONO}; font-size: 11.5px; display: flex; flex-direction: column; gap: 4px">{mat}</div>')
             + section("Presets", ct)
             + footer(btn("Supercell…"), btn("Apply", True)))
    tb = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("cube", "Cell", True, "Unit cell"), tbtn("layers", "Images", False, "2 × 2 × 2"), spacer())
    center = f'<div style="flex-grow: 1; min-width: 0; position: relative; background: {BG0}">{view}{hud(chip("PET · triclinic cell", TEXT, "#16191Ccc"), chip("2 × 2 × 2 cells · origin cell in amber", MUTED, "#16191Ccc"))}</div>'
    return studio_shell("PET_crystal.caps", "CAPS Studio — unit cell editor", tb, center, right,
                        (f"<span>V = {V:.2f} Å³ · ρ = {rho:.4f} g/cm³</span>", "<span>volume from a, b, c, α, β, γ</span>"))


# ---------------------------------------------------------------- SASA
BONDI = {"C": 1.70, "H": 1.20, "O": 1.52, "N": 1.55}


def sphere_pts(n):
    g = math.pi * (3 - math.sqrt(5))
    return [(math.cos(g * i) * math.sqrt(1 - (1 - 2 * (i + 0.5) / n) ** 2), math.sin(g * i) * math.sqrt(1 - (1 - 2 * (i + 0.5) / n) ** 2), 1 - 2 * (i + 0.5) / n) for i in range(n)]


def sasa(atoms, probe=1.4, n=200):
    pts = sphere_pts(n)
    R = [BONDI[a["e"]] + probe for a in atoms]
    P = [a["p"] for a in atoms]
    nb = [[j for j in range(len(atoms)) if j != i and dist(P[i], P[j]) < R[i] + R[j]] for i in range(len(atoms))]
    out = []
    for i in range(len(atoms)):
        free = 0
        for u in pts:
            q = add(P[i], mul(u, R[i]))
            if all(dist(q, P[j]) >= R[j] for j in nb[i]):
                free += 1
        out.append(4 * math.pi * R[i] ** 2 * free / n)
    return out


def sasa_board():
    atoms, backbone, stereo = polystyrene(6, "atactic", seed=11)
    area = sasa(atoms, 1.4, 200)
    total = sum(area)
    conv = [(n, sum(sasa(atoms, 1.4, n))) for n in (30, 60, 120, 200, 400)]
    ring = {i for i, a in enumerate(atoms) if a.get("ar")}
    # attribute H to its nearest heavy atom
    heavy = [i for i, a in enumerate(atoms) if a["e"] != "H"]
    owner = {i: (i if atoms[i]["e"] != "H" else min(heavy, key=lambda k: dist(atoms[i]["p"], atoms[k]["p"]))) for i in range(len(atoms))}
    ph = sum(area[i] for i in range(len(atoms)) if owner[i] in ring)
    bb = total - ph
    maxs = [4 * math.pi * (BONDI[a["e"]] + 1.4) ** 2 for a in atoms]
    frac = [area[i] / maxs[i] for i in range(len(atoms))]
    sc = Scene("sa", 560, 420, yaw=-0.15, pitch=0.45, persp=0.2, fog=0.3)
    sc.add_atoms([dict(a, c=mix("#39414A", SEL, frac[i] ** 0.8)) for i, a in enumerate(atoms)])
    view = sc.svg()
    cp = plot(420, 190, [(conv, ACC, 2, None)], (0, 420), (math.floor(min(v for _, v in conv) / 50) * 50, math.ceil(max(v for _, v in conv) / 50) * 50),
              [0, 100, 200, 300, 400], None or [math.floor(min(v for _, v in conv) / 50) * 50, math.ceil(max(v for _, v in conv) / 50) * 50], "points per atom", "SASA (Å²)", pad=(50, 10, 14, 30),
              markers=[(a, b, ACC) for a, b in conv])
    ct = table(["Points", "SASA (Å²)", "Δ vs 400"], [[str(n), f"{v:.1f}", f"{(v - conv[-1][1]) / conv[-1][1] * 100:+.2f} %"] for n, v in conv], ["30%", "36%", "34%"], mono_cols=(0, 1, 2), align_right=(1, 2), fs=12, rowh=25, hl={3})
    res = table(["Part", "Å²", "%"], [["Whole chain", f"{total:.1f}", "100"], ["Phenyl rings (+ their H)", f"{ph:.1f}", f"{ph / total * 100:.1f}"], ["Backbone (+ its H)", f"{bb:.1f}", f"{bb / total * 100:.1f}"]],
                ["54%", "26%", "20%"], mono_cols=(1, 2), align_right=(1, 2), fs=12, rowh=27, hl={0})
    left = (f'<div style="width: 600px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Exposure per atom", col(f'<div style="background: {BG0}; border-radius: 8px; overflow: hidden">{view}</div>',
                                            f'<div style="display: flex; align-items: center; gap: 10px; font-size: 11.5px; color: {DIM}"><span>buried</span><div style="flex-grow: 1; height: 8px; border-radius: 4px; background: linear-gradient(90deg, #39414A, {SEL})"></div><span>fully exposed</span></div>', gap=8), chip(f"{len(atoms)} atoms", MUTED, BG2, True), 14)
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Solvent-accessible surface area", col(res, row(field("Probe radius", "1.40", "Å"), select("Radii", "Bondi 1964"), select("Points", "200 per atom"), gap=8), gap=10), chip("computed", OK, BG2), 14)
             + card("Point-count convergence", f'<div style="display: flex; gap: 14px; align-items: flex-start"><div style="flex-shrink: 0">{cp}</div><div style="flex-grow: 1; min-width: 0">{ct}</div></div>', "", 14)
             + f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Shrake–Rupley (J. Mol. Biol. 79, 351, 1973) on the drawn PS hexamer, golden-spiral points. Radii: Bondi, J. Phys. Chem. 68, 441 (1964). A trajectory run reports the mean and spread per frame.</div>'
             + '</div>')
    content = (pagehead("Surface area", "Solvent-accessible area per atom and per group", row(btn("Export per-atom CSV", ic="download"), btn("Run on trajectory", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Analyze", "CAPS Analyze — surface area", "Analyze › Surface area", content, ("<span>PS · atactic · DP 6</span>", f"<span>SASA {total:.1f} Å² · probe 1.40 Å</span>"))


# ---------------------------------------------------------------- Molecule inspector
MASS = {"C": 12.011, "H": 1.008}
MONO_M = {"C": 12.0, "H": 1.00782503207}


def inspector():
    atoms = styrene3d()
    nC = sum(1 for a in atoms if a["e"] == "C"); nH = sum(1 for a in atoms if a["e"] == "H")
    mw = nC * MASS["C"] + nH * MASS["H"]
    mono = nC * MONO_M["C"] + nH * MONO_M["H"]
    dbe = nC - nH / 2 + 1
    m = [MASS[a["e"]] for a in atoms]
    com = mul(tuple(sum(m[i] * atoms[i]["p"][k] for i in range(len(atoms))) for k in range(3)), 1 / sum(m))
    I = [[0.0] * 3 for _ in range(3)]
    for i, a in enumerate(atoms):
        r = sub(a["p"], com)
        r2 = dot3(r, r)
        for x in range(3):
            for y in range(3):
                I[x][y] += m[i] * ((r2 if x == y else 0) - r[x] * r[y])
    Ia, Ib, Ic = jacobi3(I)
    rg = math.sqrt(sum(m[i] * dist(a["p"], com) ** 2 for i, a in enumerate(atoms)) / sum(m))
    planar = abs(Ic - (Ia + Ib))
    bonds = perceive_bonds(atoms)
    sc = Scene("mi", VW, VH - 40, yaw=0.3, pitch=0.9, persp=0.15, fog=0.2, atom_k=1.15, bond_w=0.24)
    sc.add_atoms(atoms, bonds)
    sc._prep()
    sc.scale *= 0.62
    view = sc.svg()
    formula = f"C<sub>{nC}</sub>H<sub>{nH}</sub>"
    right = (panel_head("Molecule", chip("styrene", MUTED, BG2))
             + section("Identity", col(kv("Formula", formula), kv("SMILES", "C=Cc1ccccc1"), kv("Atoms · bonds", f"{len(atoms)} · {len(bonds)}"), gap=6))
             + section("Mass", col(kv("Molecular weight", f"{mw:.3f} g/mol"), kv("Monoisotopic", f"{mono:.4f} u"), kv("Double-bond equivalents", f"{dbe:.0f}"), gap=6))
             + section("Shape (from this geometry)", col(kv("I<sub>A</sub> · I<sub>B</sub> · I<sub>C</sub>", f"{Ia:.1f} · {Ib:.1f} · {Ic:.1f}"), kv("units", "amu·Å²"), kv("I<sub>C</sub> − I<sub>A</sub> − I<sub>B</sub>", f"{abs(Ic - Ia - Ib):.3f} amu·Å²"),
                                                            kv("Radius of gyration (mass)", f"{rg:.3f} Å"), gap=6))
             + section("Needs a calculation", col(kv("Dipole moment", "[result] D"), kv("Polarisability", "[result] Å³"), kv("Conformers", "1 rotatable bond · scan →"), gap=6))
             + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Standard atomic weights (IUPAC); monoisotopic ¹²C and ¹H. The inertia defect is zero because the drawn molecule is planar.</div>'
             + footer(btn("Copy table"), btn("Add to library", True)))
    tb = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("atom", "Inspect", True, "Inspector"), spacer())
    center = f'<div style="flex-grow: 1; min-width: 0; position: relative; background: {BG0}">{view}{hud(chip("styrene · planar", TEXT, "#16191Ccc"), chip("C–C ar 1.39 · C=C 1.34 · C–H 1.08 Å", MUTED, "#16191Ccc", True))}</div>'
    return studio_shell("styrene.caps", "CAPS Studio — molecule inspector", tb, center, right,
                        (f"<span>{formula.replace('<sub>', '').replace('</sub>', '')} · {mw:.2f} g/mol</span>", "<span>values computed from the drawn geometry</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage11/project", exist_ok=True)
    for name, fn in (("CellEditor", cell_editor), ("SurfaceArea", sasa_board), ("MoleculeInspector", inspector)):
        open(f"stage11/project/{name}.dc.html", "w").write(fn())
        print(name, "ok")
