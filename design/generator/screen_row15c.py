import math
from lib import *
from mols import *
from screen_app import shell, pagehead, card
from display import style_bar

AMU_A3_TO_GCM3 = 1.66054
MASS = {"Si": 28.086, "O": 15.999, "CH2": 14.027}


def interface_cell():
    slab = silica_slab(8, 8, 3)
    Lx = Ly = 8 * 2.6
    zs = [a["p"][2] for a in slab]
    top = max(zs)
    chains = []
    for k in range(10):
        rnd = random.Random(500 + k)
        st = (rnd.uniform(1, Lx - 1), rnd.uniform(1, Ly - 1), rnd.uniform(top + 4, top + 26))
        chains.append(random_chain(st, 44, 600 + k, bond=1.54, box=((0.3, 0.3, top + 3.2), (Lx - 0.3, Ly - 0.3, top + 30)), persistence=0.62))
    return slab, chains, Lx, Ly, top


def interface():
    slab, chains, Lx, Ly, top = interface_cell()
    dz = 1.0
    zmax = top + 32
    nb = int(math.ceil((zmax + 2) / dz))
    z0 = -2.0
    prof = {g: [0.0] * nb for g in ("Si", "O", "CH2")}
    for a in slab:
        prof[a["e"]][int((a["p"][2] - z0) // dz)] += MASS[a["e"]]
    for ch in chains:
        for p in ch:
            prof["CH2"][int((p[2] - z0) // dz)] += MASS["CH2"]
    vol = Lx * Ly * dz
    rho = {g: [(z0 + (i + 0.5) * dz, v / vol * AMU_A3_TO_GCM3) for i, v in enumerate(vals)] for g, vals in prof.items()}
    poly = [r for _, r in rho["CH2"]]
    bulk_rng = [r for z, r in rho["CH2"] if top + 10 <= z <= top + 26]
    bulk = sum(bulk_rng) / len(bulk_rng)
    half = next(z for z, r in rho["CH2"] if z > top and r >= bulk / 2)
    ymax = max(max(r for _, r in s) for s in rho.values())
    ytop = math.ceil(ymax * 2) / 2
    pl = plot(600, 230, [(rho["Si"], ELEM["Si"][0], 1.8, None), (rho["O"], ELEM["O"][0], 1.8, None), (rho["CH2"], ACC, 2.2, None),
                         ([(z0, bulk), (zmax, bulk)], DIM, 1, "4 4")],
              (z0, zmax), (0, ytop), [0, 10, 20, 30, 40], [round(ytop * i / 4, 2) for i in range(5)], "z (Å) · normal to the surface", "ρ (g/cm³)",
              pad=(48, 12, 14, 32), bands=[(top, half, SEL, 0.10)])
    # 3D view
    sc = Scene("if", 520, 560, yaw=0.55, pitch=0.18, persp=0.25, fog=0.45, atom_k=0.9)
    sc.add_atoms(slab)
    for i, ch in enumerate(chains):
        sc.add_tube(ch, mix(ACC, "#8B969E", (i % 4) / 5), 0.5)
    sc.add_box((0, 0, -1.5), (Lx, 0, 0), (0, Ly, 0), (0, 0, zmax + 1.5), MUTED, 1.0, "4 4", 0.6)
    view = sc.svg()
    nbeads = sum(len(c) for c in chains)
    res = table(["Quantity", "This cell", "From trajectory"],
                [["Slab top (outer Si/O)", f"{top:.2f} Å", "—"], ["Polymer plateau ρ", f"{bulk:.2f} g/cm³", "[result]"],
                 ["Polymer reaches ½ plateau", f"{half:.1f} Å", "[result]"], ["Gap to surface", f"{half - top:.1f} Å", "[result]"],
                 ["Adsorbed layer thickness", "—", "[result]"], ["Work of adhesion", "—", "[result] mJ/m²"]],
                ["44%", "28%", "28%"], mono_cols=(1, 2), align_right=(1, 2), fs=12, rowh=25)
    lg = row(chip(f'{dot(ELEM["Si"][0])} Si'), chip(f'{dot(ELEM["O"][0])} O'), chip(f'{dot(ACC)} PE united atom (CH₂)'), chip(f'{dot(DIM)} plateau'), chip(f'{dot(SEL)} gap'), gap=6)
    setup = row(select("Axis", "z · surface normal"), field("Bin", f"{dz:.1f}", "Å"), select("Groups", "Si · O · polymer"), select("Frames", "this frame"), gap=8)
    left = (f'<div style="width: 540px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG0}; border: 1px solid {LINE}; border-radius: 10px; overflow: hidden; position: relative">{view}'
            + style_bar("Backbone", f"PE as backbone · H hidden · {len(chains) * (44 + 90):,} PE atoms in model".replace(",", " "), "right: 10px; bottom: 10px")
            + f'<div style="position: absolute; left: 12px; top: 10px; display: flex; gap: 6px">{chip("SiO₂ slab + PE melt", TEXT, "#16191Ccc")}{chip(f"{len(slab)} + {nbeads} sites", MUTED, "#16191Ccc", True)}</div></div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Density profile", col(setup, pl, lg, gap=10), chip("computed from the drawn cell", DIM, BG2), 14)
             + card("Interface", col(res, f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Middle column: one frame of this Grow snapshot ({len(chains)} chains), not an equilibrated system. After Equilibrate, CAPS averages over frames and fills the right-hand column.</div>', gap=8), "", 14)
             + '</div>')
    content = (pagehead("Interface analysis", "Density along the surface normal, layer thickness, adhesion", row(btn("Export CSV", ic="download"), btn("Run on trajectory", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Analyze", "CAPS Analyze — interface", "Analyze › Interface", content,
                 (f"<span>SiO2_PE_interface.caps · {Lx:.1f} × {Ly:.1f} Å</span>", "<span>Periodic in x, y · vacuum above the melt</span>"))


def diffusion():
    kB, T, eta, L, xi = 1.380649e-23, 298.15, 0.890e-3, 3.0e-9, 2.837297
    dD = kB * T * xi / (6 * math.pi * eta * L)
    # illustrative MSD shape: ballistic -> caged -> diffusive
    pts = []
    for i in range(0, 61):
        t = 10 ** (-3 + i * 0.1)  # ns
        m = 1 / (1 / (40 * t * t) + 1 / (0.02 + 2.0 * t)) if t > 0 else 0
        pts.append((math.log10(t), math.log10(m)))
    guide = [(0.0, math.log10(2.0)), (3.0, math.log10(2.0e3))]
    fitb = (-0.0, 2.3)
    pl = plot(600, 300, [(pts, ACC, 2.2, None), (guide, DIM, 1.2, "5 4")], (-3, 3), (-4, 4), [-3, -2, -1, 0, 1, 2, 3], [-4, -2, 0, 2, 4],
              "log₁₀ t (ns)", "log₁₀ MSD (Å²)", pad=(46, 12, 14, 32), bands=[(fitb[0], fitb[1], SEL, 0.10)])
    # local slope of the illustrative curve, to show the diagnostic
    slope = [((pts[i][0] + pts[i + 1][0]) / 2, (pts[i + 1][1] - pts[i][1]) / (pts[i + 1][0] - pts[i][0])) for i in range(len(pts) - 1)]
    sp = plot(600, 130, [(slope, SEL, 1.8, None), ([(-3, 1), (3, 1)], DIM, 1, "4 4")], (-3, 3), (0, 2.2), [-3, -2, -1, 0, 1, 2, 3], [0, 1, 2],
              "", "slope", pad=(46, 8, 14, 20), bands=[(fitb[0], fitb[1], SEL, 0.10)])
    setup = col(row(select("Molecules", "water · 900"), select("Position", "centre of mass"), gap=8),
                row(select("Time origins", "every 10 ps"), select("Remove drift", "system COM"), gap=8),
                row(field("Fit from", "1", "ns"), field("Fit to", "200", "ns"), select("Unwrap", "yes · image flags"), gap=8), gap=10)
    fs = col(row(field("Box edge L", f"{L * 1e9:.1f}", "nm"), field("Viscosity η", f"{eta * 1e3:.3f}", "mPa·s"), field("T", f"{T:.2f}", "K"), gap=8),
             f'<div style="font-family: {MONO}; font-size: 12.5px; padding: 10px 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px">D<sub>∞</sub> = D<sub>PBC</sub> + k<sub>B</sub>T ξ / (6π η L) &nbsp; · &nbsp; ξ = {xi}</div>',
             kv("Correction k<sub>B</sub>T ξ / 6πηL", f"{dD:.2e} m²/s".replace("e-", " × 10⁻").replace("10⁻10", "10⁻¹⁰").replace("10⁻09", "10⁻⁹")),
             f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Computed from the three inputs above (η of water at 298 K, 3 nm box). CAPS takes η from a Green–Kubo run or asks you for it. Yeh &amp; Hummer, <i>J. Phys. Chem. B</i> 108, 15873 (2004).</div>', gap=10)
    res = table(["Quantity", "Value"], [["D<sub>PBC</sub> (fit)", "[result] m²/s"], ["D<sub>∞</sub> (corrected)", "[result] m²/s"], ["Slope in fit window", "[result]"], ["Replicas · spread", "[result]"]],
                ["55%", "45%"], mono_cols=(1,), fs=12, rowh=25)
    left = (f'<div style="width: 640px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Mean-squared displacement", col(pl, sp, f'<div style="font-size: 11.5px; color: {DIM}">Illustrative curve. Fit only where the log–log slope is 1; CAPS warns if the window includes the caged or ballistic part.</div>', gap=6), chip("illustrative", ERR, BG2), 14)
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Setup", setup, "", 14)
             + card("Finite-size correction", fs, chip("computed", OK, BG2), 14)
             + card("Result", res, "", 14)
             + '</div>')
    content = (pagehead("Diffusion", "Einstein MSD, fit-window check and box-size correction", row(btn("Export CSV", ic="download"), btn("Run", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Analyze", "CAPS Analyze — diffusion", "Analyze › Diffusion", content,
                 ("<span>water_box.caps · trajectory prod-2 · 400 ns</span>", "<span>Unwrapped with image flags</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage6/project", exist_ok=True)
    for name, fn in (("Interface", interface), ("Diffusion", diffusion)):
        open(f"stage6/project/{name}.dc.html", "w").write(fn())
