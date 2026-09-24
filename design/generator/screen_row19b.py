import math
from lib import *
from mols import polystyrene
from screen_app import shell, pagehead, card
from screen_studio2 import studio_shell, bar, spacer
from screen_builders import footer, hud
from screen_row16b import ps_bonds
from display import ps_cell, render_cell, count, style_bar

VW = 1440 - 72 - 380
VH = 900 - 44 - 48 - 26
M = {"C": 12.011, "H": 1.008}


def pe_chain(n=20):
    """All-atom n-alkane: zigzag C–C 1.53 Å, C–H 1.09 Å, returns atoms, bonds, carbon indices."""
    dx, dy = 1.27, 0.43
    C = [((k - n / 2) * dx, dy * (1 if k % 2 == 0 else -1), 0.0) for k in range(n)]
    atoms = [{"e": "C", "p": p} for p in C]
    bonds = [(k, k + 1) for k in range(n - 1)]
    for k, c in enumerate(C):
        t = norm(sub(C[min(k + 1, n - 1)], C[max(k - 1, 0)]))
        up = (0.0, 1.0 if k % 2 == 0 else -1.0, 0.0)
        side = norm(cross(t, up))
        for s in (1, -1):
            atoms.append({"e": "H", "p": add(c, mul(norm(add(mul(up, 0.58), mul(side, s * 0.81))), 1.09))}); bonds.append((k, len(atoms) - 1))
        if k in (0, n - 1):
            d = mul(t, -1 if k == 0 else 1)
            atoms.append({"e": "H", "p": add(c, mul(norm(add(mul(d, 0.94), mul(up, -0.34))), 1.09))}); bonds.append((k, len(atoms) - 1))
    return atoms, bonds, list(range(n))


def model_resolution():
    n = 20
    atoms, bonds, C = pe_chain(n)
    nH = sum(1 for a in atoms if a["e"] == "H")
    m_aa = sum(M[a["e"]] for a in atoms)
    ua = [{"e": "C", "p": atoms[k]["p"], "r": 0.95, "c": "#8B969E"} for k in C]
    m_ua = sum(M["C"] + M["H"] * (3 if k in (0, n - 1) else 2) for k in C)
    per = 5
    beads = []
    for b in range(n // per):
        idx = C[b * per:(b + 1) * per]
        com = mul(tuple(map(sum, zip(*[atoms[k]["p"] for k in idx]))), 1 / per)
        mass = sum(M["C"] + M["H"] * (3 if k in (0, n - 1) else 2) for k in idx)
        beads.append((com, mass))
    m_cg = sum(m for _, m in beads)
    w = (VW - 66) // 3
    def pan(sid, title, sub_, fill):
        sc = Scene(sid, w, 300, yaw=0.25, pitch=0.9, persp=0.1, fog=0.15, atom_k=1.05)
        fill(sc)
        return (f'<div style="width: {w}px; display: flex; flex-direction: column; gap: 6px"><div style="display: flex; align-items: baseline; gap: 8px"><span style="font-size: 13px; font-weight: 600">{title}</span><span style="font-size: 11.5px; color: {DIM}">{sub_}</span></div>'
                f'<div style="background: {BG0}; border: 1px solid {LINE}; border-radius: 8px; overflow: hidden">{sc.svg()}</div></div>')
    p1 = pan("mr1", "All-atom", f"{len(atoms)} sites · {nH} H", lambda sc: sc.add_atoms(atoms, bonds))
    p2 = pan("mr2", "United-atom", f"{len(ua)} sites · H folded in", lambda sc: sc.add_atoms(ua, bonds[:n - 1]))
    def cg(sc):
        sc.add_atoms([{"e": "C", "p": a["p"], "r": 0.001, "c": BG0} for a in ua], [])
        sc.add_atoms([{"e": "C", "p": p, "r": 2.6, "c": ACC} for p, _ in beads], [(i, i + 1) for i in range(len(beads) - 1)])
    p3 = pan("mr3", "Coarse-grained", f"{len(beads)} beads · {per} CH<sub>x</sub> each", cg)
    tb = table(["", "All-atom", "United-atom", "Coarse-grained"],
               [["Sites", str(len(atoms)), str(len(ua)), str(len(beads))], ["Hydrogens", f"{nH} explicit", "inside CH₂/CH₃", "inside beads"],
                ["Mass (g/mol)", f"{m_aa:.3f}", f"{m_ua:.3f}", f"{m_cg:.3f}"], ["Example force field", "GAFF2 · OPLS-AA", "TraPPE-UA", "Kremer–Grest (generic)"],
                ["Charges", "on every atom incl. H", "on united sites", "usually none"], ["Use it for", "chemistry, H-bonds, spectra", "alkanes, melts, faster", "long times, entanglement"]],
               ["20%", "26%", "26%", "28%"], mono_cols=(1, 2, 3), fs=12, rowh=30)
    same = abs(m_aa - m_ua) < 1e-9 and abs(m_aa - m_cg) < 1e-9
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 14px 16px; background: {BG0}; overflow: hidden">'
              f'<div style="display: flex; gap: 16px">{p1}{p2}{p3}</div>'
              + card("Same chain, three resolutions · n-C₂₀H₄₂", tb, chip("mass conserved" if same else "mass differs", OK if same else ERR, BG2), 14)
              + '</div>')
    right = (panel_head("Model resolution", chip("per document", MUTED, BG2))
             + section("This document", col(seg(["All-atom", "United-atom", "Coarse-grained"], "All-atom", full=True), f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">All-atom is the default. Changing resolution is a conversion that creates a new document; the original stays.</span>', gap=8))
             + section("Converting down", col(kv("All-atom → united-atom", "merge H into carbons"), kv("→ coarse-grained", f"{per} carbons per bead, centre of mass"), gap=6))
             + section("Converting back up", col(kv("Coarse-grained → all-atom", "backmap, then relax"), kv("United-atom → all-atom", "add H at 1.09 Å, then relax"), f'<a href="CoarseGrained.dc.html" style="font-size: 12.5px">Open the backmap board →</a>', gap=6))
             + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Masses from IUPAC standard atomic weights; united-atom and bead masses include their hydrogens, so the chain mass is identical in all three.</div>')
    tbar = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("layers", "Resolution", True, "Model resolution"), spacer())
    return studio_shell("C20H42.caps", "CAPS Studio — model resolution", tbar, center, right, (f"<span>n-C₂₀H₄₂ · {len(atoms)} atoms · {m_aa:.3f} g/mol</span>", "<span>all-atom model</span>"))


def add_hydrogens():
    atoms, backbone, stereo = polystyrene(4, "atactic", seed=6)
    bonds = ps_bonds(atoms, backbone, stereo)
    heavy = [i for i, a in enumerate(atoms) if a["e"] != "H"]
    hmap = {i: n for n, i in enumerate(heavy)}
    hb = [(hmap[a], hmap[b]) for a, b in bonds if a in hmap and b in hmap]
    deg = [0] * len(heavy)
    for a, b in hb:
        deg[a] += 1; deg[b] += 1
    val = [3 if atoms[i].get("ar") else 4 for i in heavy]  # aromatic C: 3 sigma neighbours + pi
    need = [v - d for v, d in zip(val, deg)]
    nH_rule = sum(need)
    nH_true = sum(1 for a in atoms if a["e"] == "H")
    assert nH_rule == nH_true, (nH_rule, nH_true)
    groups = {}
    for i, k in enumerate(heavy):
        key = ("aromatic C" if atoms[k].get("ar") else "sp³ C") + f" with {deg[i]} C neighbour{'s' if deg[i] != 1 else ''}"
        groups.setdefault(key, [0, 0]); groups[key][0] += 1; groups[key][1] += need[i]
    gt = table(["Atom", "Count", "H added"], [[k, str(v[0]), str(v[1])] for k, v in sorted(groups.items())], ["58%", "20%", "22%"], mono_cols=(1, 2), align_right=(1, 2), fs=12, rowh=27)
    w = (VW - 36) // 2
    def pan(sid, title, sub_, ats, bds, halo=None):
        sc = Scene(sid, w, VH - 250, yaw=-0.1, pitch=0.5, persp=0.15, fog=0.25, atom_k=1.1)
        sc.add_atoms(ats, bds)
        if halo: sc.halo = halo
        return (f'<div style="width: {w}px; display: flex; flex-direction: column; gap: 6px"><div style="display: flex; align-items: baseline; gap: 8px"><span style="font-size: 13px; font-weight: 600">{title}</span><span style="font-size: 11.5px; color: {DIM}">{sub_}</span></div>'
                f'<div style="background: {BG0}; border: 1px solid {LINE}; border-radius: 8px; overflow: hidden">{sc.svg()}</div></div>')
    before = pan("ah1", "As imported", f"{len(heavy)} heavy atoms · 0 H", [atoms[i] for i in heavy], hb)
    after = pan("ah2", "After Add hydrogens", f"{len(atoms)} atoms · +{nH_true} H", atoms, bonds, {i for i, a in enumerate(atoms) if a["e"] == "H"})
    strip = (f'<div style="flex-grow: 1; display: flex; gap: 22px; padding: 14px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
             f'<div style="width: 460px; flex-shrink: 0">{gt}</div>'
             f'<div style="display: flex; flex-direction: column; gap: 8px; font-size: 12.5px; line-height: 1.5"><span>H per atom = valence − heavy neighbours (sp³ C: 4, aromatic C: 3).</span>'
             f'<span style="color: {MUTED}">Rule total {nH_rule} H = the {nH_true} H of the reference structure.</span><span style="color: {DIM}">New H sit at 1.09 Å with ideal angles, then only H are relaxed.</span></div></div>')
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">'
              f'<div style="display: flex; gap: 12px; padding: 14px 12px 12px">{before}{after}</div>{strip}</div>')
    right = (panel_head("Add hydrogens", chip("PS_frag.pdb", MUTED, BG2, True))
             + section("Source", col(kv("File", "PDB, heavy atoms only"), kv("Bond orders", "perceived · 6 aromatic per ring"), gap=6))
             + section("Options", col(toggle("Keep hydrogens already present", True), toggle("Relax new H only", True), toggle("Protonation by pH (biomolecules)", False), row(field("pH", "7.0"), select("Tool", "built-in rules"), gap=8), gap=8))
             + section("Result", col(kv("Heavy atoms", str(len(heavy))), kv("Hydrogens added", str(nH_true), vcol=OK), kv("Net charge", "0"), gap=6))
             + footer(btn("Undo"), btn("Add hydrogens", True)))
    tbar = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("plus", "Add H", True, "Add hydrogens"), spacer())
    return studio_shell("PS_frag.caps", "CAPS Studio — add hydrogens", tbar, center, right, (f"<span>{len(atoms)} atoms after · {nH_true} H</span>", "<span>computed from valence rules</span>"))


def lens_view():
    L = 52.0
    cell = ps_cell(34, 8, L, seed=31)
    cen = (L / 2, L / 2, L / 2)
    R = 10.0
    sc = Scene("lv", VW, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.45)
    cols = ["#F0A83C", "#6CC4D8", "#DE775D", "#9B7AD5", "#7DC884", "#D6AC5C", "#E9ECEF", "#2271DB"]
    render_cell(sc, cell, "Backbone", colours=[mix(c, BG0, 0.35) for c in cols], tube_w=0.45)
    inside, inH = [], 0
    for atoms, bonds, bb in cell:
        keep = [j for j, a in enumerate(atoms) if dist(a["p"], cen) <= R]
        idx = {j: n for n, j in enumerate(keep)}
        sc.add_atoms([atoms[j] for j in keep], [(idx[a], idx[b]) for a, b in bonds if a in idx and b in idx])
        inside += keep; inH += sum(1 for j in keep if atoms[j]["e"] == "H")
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, "5 4", 0.6)
    sc._prep()
    c2 = sc.P(cen)
    rpx = R * sc.scale * c2[3]
    ring = f'<circle cx="{c2[0]:.1f}" cy="{c2[1]:.1f}" r="{rpx:.1f}" fill="none" stroke="{SEL}" stroke-width="1.6" stroke-dasharray="6 4"></circle>'
    view = sc.svg(overlay=ring + label_pill(c2[0] + rpx * 0.72, c2[1] - rpx * 0.9, f"lens · {R:.0f} Å", SEL, border=SEL))
    n, h = count(cell)
    right = (panel_head("Lens", chip("hold L to move", MUTED, BG2, True))
             + section("Inside the lens", col(kv("Style", "All atoms"), kv("Atoms", f"{len(inside)}"), kv("of them H", f"{inH}"), kv("Radius", f"{R:.0f} Å"), gap=6))
             + section("Outside", col(kv("Style", "Backbone"), kv("Chains", f"{len(cell)}"), kv("Atoms in model", f"{n:,}".replace(",", " ")), gap=6))
             + section("Lens options", col(row(field("Radius", f"{R:.0f}", "Å"), select("Follows", "cursor"), gap=8), toggle("Dim outside", True), toggle("Measurements only inside", False), gap=8))
             + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">The lens is how large systems stay readable without losing atoms: detail where you look, tubes elsewhere, every atom still in the model.</div>')
    tbar = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("search", "Lens", True, "All-atom lens"), spacer())
    center = (f'<div style="flex-grow: 1; min-width: 0; position: relative; background: {BG0}">{view}'
              + style_bar("Backbone", f"Backbone outside · all atoms inside the lens · {n:,} atoms in model".replace(",", " ")) + '</div>')
    return studio_shell("PS_melt_34.caps", "CAPS Studio — all-atom lens", tbar, center, right, (f"<span>{len(cell)} chains · {n:,} atoms · {h} H</span>".replace(",", " "), f"<span>lens {len(inside)} atoms</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage15/project", exist_ok=True)
    for name, fn in (("ModelResolution", model_resolution), ("AddHydrogens", add_hydrogens), ("LensView", lens_view)):
        h = fn(); open(f"stage15/project/{name}.dc.html", "w").write(h); print(name, len(h))
