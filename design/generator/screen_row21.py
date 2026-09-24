"""More Visualize pipeline steps: clusters, RDF & coordination, slice & binning, expression selection."""
import math
from lib import *
from screen_app import card
from screen_studio2 import bar, spacer
from display import ps_cell, style_bar
from screen_row20 import annotate, draw, _shell, CHAIN, VIRIDIS, ramp

L = 33.0
MASS = {"C": 12.011, "H": 1.008}
NA = 6.02214076e23


def cell_recs():
    cell = ps_cell(10, 8, L, seed=21)
    return cell, annotate(cell)


def mi(a, b):
    d = [b[k] - a[k] for k in range(3)]
    return math.sqrt(sum((x - L * round(x / L)) ** 2 for x in d))


def inspector(title, body, w=330):
    return (f'<aside style="width: {w}px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">'
            + panel_head(title, chip("pipeline step", MUTED, BG2)) + body + '</aside>')


def viewport(sc_svg, note, style="All atoms"):
    return (f'<div style="position: relative; flex-grow: 1; min-height: 0; background: {BG0}; overflow: hidden">{sc_svg}'
            + style_bar(style, note, "right: 12px; top: 10px") + '</div>')


def frame(content_center, right, title, crumb, status):
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
              + bar(tbtn("rotate", "Orbit", True), tbtn("move", "Pan"), tbtn("search", "Zoom"), sep(), tbtn("cursor", "Pick"), tbtn("ruler", "Measure"), spacer(), chip("frame 120 / 500", MUTED, BG2, True))
              + content_center + '</div>')
    return _shell(title, f'{icon("chevr", 12, DIM)}<span>{crumb}</span>', f'<div style="flex-grow: 1; display: flex; min-height: 0">{center}{right}</div>', status)


CW = 1368 - 330
BOTTOM = 230
VH = 900 - 44 - 26 - 48 - BOTTOM - 28


# ---------------------------------------------------------------- Cluster analysis
def clusters():
    cell, recs = cell_recs()
    heavy = {}
    for r in recs:
        if r["e"] != "H": heavy.setdefault(r["mol"], []).append(r["p"])
    mols = sorted(heavy)
    def comps(cut):
        parent = {m: m for m in mols}
        def f(x):
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for i, a in enumerate(mols):
            for b in mols[i + 1:]:
                if f(a) != f(b) and any(mi(p, q) < cut for p in heavy[a][::1] for q in heavy[b]):
                    parent[f(a)] = f(b)
        groups = {}
        for m in mols: groups.setdefault(f(m), []).append(m)
        return sorted(groups.values(), key=lambda g: (-len(g), g))
    sweep = [(c, comps(c)) for c in (3.1, 3.3, 3.6, 4.0)]
    cut, groups = sweep[1]
    cid = {m: i for i, g in enumerate(groups) for m in g}
    cols = ["#F0A83C", "#2271DB", "#7DC884", "#DE775D", "#9B7AD5", "#6CC4D8", "#D6AC5C", "#E9ECEF", "#C77DBA", "#8FB8A8"]
    sc = Scene("cl", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.25)
    draw(sc, cell, recs, lambda r: cols[cid[r["mol"]] % 10] if r["e"] != "H" else mix(cols[cid[r["mol"]] % 10], "#FFFFFF", 0.55))
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    natoms = {m: sum(1 for r in recs if r["mol"] == m) for m in mols}
    ct = table(["Cluster", "Molecules", "Atoms"], [[f'<span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 10px; height: 10px; border-radius: 2px; background: {cols[i % 10]}"></span>{i + 1}</span>', ", ".join(map(str, g)), str(sum(natoms[m] for m in g))] for i, g in enumerate(groups)],
               ["22%", "52%", "26%"], mono_cols=(1, 2), align_right=(2,), fs=11.5, rowh=24)
    st = table(["Cutoff (Å)", "Clusters", "Largest"], [[f"{c:.1f}", str(len(g)), f"{len(g[0])} mol"] for c, g in sweep], ["36%", "30%", "34%"], mono_cols=(0, 1, 2), align_right=(1, 2), fs=11.5, rowh=24, hl={1})
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 12px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex: 1 1 0; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Clusters at {cut} Å</h3>{ct}</div>'
              f'<div style="width: 330px; flex-shrink: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Cutoff sweep</h3>{st}</div></div>')
    right = inspector("Cluster analysis", section("Neighbours", col(select("Mode", "cutoff between heavy atoms"), row(field("Cutoff", f"{cut}", "Å"), select("Unit", "whole molecules"), gap=8), toggle("Periodic (minimum image)", True), gap=8))
                      + section("Output", col(toggle("Colour by cluster", True), toggle("Sort by size", True), toggle("Compute Rg per cluster", False), toggle("Unwrap clusters", False), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Two molecules join a cluster when any heavy atoms are closer than the cutoff. Computed on this frame; a built, un-equilibrated cell.</div>')
    return frame(viewport(sc.svg(), f"Colour: cluster · {len(groups)} clusters at {cut} Å") + bottom, right, "CAPS — cluster analysis", "Analyze › Visualize › Cluster analysis",
                 (f"<span>{len(recs)} atoms · {len(mols)} molecules</span>", f"<span>{len(groups)} clusters</span>"))


# ---------------------------------------------------------------- RDF & coordination
def rdf_board():
    cell, recs = cell_recs()
    C = [r for r in recs if r["e"] == "C"]
    rmax, dr = 12.0, 0.2
    nb = int(rmax / dr)
    hist = [0] * nb
    coord = {r["id"]: 0 for r in C}
    for i, a in enumerate(C):
        for b in C[i + 1:]:
            if a["mol"] == b["mol"]: continue
            d = mi(a["p"], b["p"])
            if d < rmax:
                hist[int(d / dr)] += 2
            if d < 5.0:
                coord[a["id"]] += 1; coord[b["id"]] += 1
    n = len(C)
    rho = n / L ** 3
    g = []
    for k in range(nb):
        r0, r1 = k * dr, (k + 1) * dr
        shell_v = 4 / 3 * math.pi * (r1 ** 3 - r0 ** 3)
        g.append(((r0 + r1) / 2, hist[k] / (n * rho * shell_v)))
    gmax = max(v for _, v in g)
    pl = plot(620, 210, [(g, ACC, 2, None), ([(0, 1), (rmax, 1)], DIM, 1, "4 4")], (0, rmax), (0, math.ceil(gmax * 2) / 2), [0, 2, 4, 6, 8, 10, 12], [0, 0.5, 1, math.ceil(gmax * 2) / 2],
              "r (Å)", "g(r) · C–C, different molecules", pad=(48, 12, 16, 32))
    cvals = list(coord.values())
    cmean = sum(cvals) / len(cvals)
    cm = max(cvals)
    hb = [0] * (cm + 1)
    for v in cvals: hb[v] += 1
    hw, hh = 300, 150
    bars = "".join(f'<rect x="{30 + i * (hw - 40) / (cm + 1):.1f}" y="{hh - 22 - v / max(hb) * (hh - 36):.1f}" width="{(hw - 40) / (cm + 1) - 2:.1f}" height="{v / max(hb) * (hh - 36):.1f}" rx="1.5" fill="{SEL}"></rect>' for i, v in enumerate(hb))
    ticks = "".join(f'<text x="{30 + (i + 0.5) * (hw - 40) / (cm + 1):.1f}" y="{hh - 7}" text-anchor="middle" font-size="9.5" font-family="IBM Plex Mono" fill="{DIM}">{i}</text>' for i in range(0, cm + 1, max(1, (cm + 1) // 8)))
    hist_svg = f'<svg width="{hw}" height="{hh}" viewBox="0 0 {hw} {hh}" role="img" aria-label="Coordination number histogram" style="display: block">{bars}{ticks}<line x1="30" y1="{hh - 22}" x2="{hw - 10}" y2="{hh - 22}" stroke="{LINE}"></line></svg>'
    sc = Scene("rd", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.25)
    draw(sc, cell, recs, lambda r: ramp(VIRIDIS, coord.get(r["id"], 0) / cm) if r["e"] == "C" else "#39414A")
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-shrink: 0">{pl}</div>'
              f'<div style="display: flex; flex-direction: column; gap: 6px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Coordination · C within 5 Å, other molecules</h3>{hist_svg}'
              f'<span style="font-family: {MONO}; font-size: 12px">mean {cmean:.2f} · max {cm}</span></div></div>')
    right = inspector("Coordination & RDF", section("Pairs", col(row(select("A", "C (types 1, 2)"), select("B", "C (types 1, 2)"), gap=8), toggle("Only different molecules", True), gap=8))
                      + section("Range", col(row(field("r max", f"{rmax:.0f}", "Å"), field("Bin", f"{dr}", "Å"), gap=8), row(field("Coordination cutoff", "5.0", "Å"), gap=8), gap=8))
                      + section("Averaging", col(select("Frames", "this frame"), toggle("Average over frames", False), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Computed from this frame with minimum-image distances in the {L:.0f} Å box. A built cell is not equilibrated, so this g(r) shows packing gaps, not melt structure.</div>')
    return frame(viewport(sc.svg(), f"Colour: coordination number · 0 → {cm}") + bottom, right, "CAPS — coordination & RDF", "Analyze › Visualize › Coordination & RDF",
                 (f"<span>{n} carbons · ρ = {rho:.4f} Å⁻³</span>", f"<span>mean coordination {cmean:.2f}</span>"))


# ---------------------------------------------------------------- Slice & spatial binning
def slice_bin():
    cell, recs = cell_recs()
    zc, w = L / 2, 8.0
    inside = [r for r in recs if abs((r["p"][2] % L) - zc) <= w / 2]
    dz = 1.0
    nb = int(L / dz)
    m = [0.0] * nb
    for r in recs:
        m[int((r["p"][2] % L) / dz) % nb] += MASS[r["e"]]
    vol_cm3 = L * L * dz * 1e-24
    prof = [((k + 0.5) * dz, m[k] / NA / vol_cm3) for k in range(nb)]
    tot = sum(MASS[r["e"]] for r in recs) / NA / (L ** 3 * 1e-24)
    ymax = max(v for _, v in prof) * 1.15
    pl = plot(600, 200, [(prof, ACC, 2, None), ([(0, tot), (L, tot)], DIM, 1, "4 4")], (0, L), (0, ymax), [0, 5, 10, 15, 20, 25], [0, round(ymax / 2, 2), round(ymax, 2)], "z (Å)", "ρ (g/cm³)",
              pad=(48, 12, 16, 30), bands=[(zc - w / 2, zc + w / 2, SEL, 0.12)])
    ids = {r["id"] for r in inside}
    sc = Scene("sl", CW - 6, VH, yaw=0.3, pitch=1.15, persp=0.2, fog=0.2, atom_k=1.25)
    draw(sc, cell, [r for r in recs], lambda r: (CHAIN[r["chain"] % 10] if r["e"] != "H" else mix(CHAIN[r["chain"] % 10], "#FFFFFF", 0.55)) if r["id"] in ids else "#20252A")
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-shrink: 0">{pl}</div>'
              f'<div style="display: flex; flex-direction: column; gap: 8px; padding-top: 6px">{kv("Atoms in slab", str(len(inside)))}{kv("of all", f"{len(recs)}")}{kv("Cell density (mean)", f"{tot:.3f} g/cm³")}'
              f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5; max-width: 330px">Profile bins every atom by z with its mass (H included). The built cell is loosely packed, so its density is below a real melt.</span></div></div>')
    right = inspector("Slice", section("Plane", col(row(field("Normal", "0 0 1"), field("Distance", f"{zc:.1f}", "Å"), gap=8), row(field("Slab width", f"{w:.0f}", "Å"), select("Keep", "inside"), gap=8), toggle("Only selected", False), gap=8))
                      + panel_head("Spatial binning", "") + section("Profile", col(row(select("Property", "mass density"), select("Direction", "z"), gap=8), row(field("Bin", f"{dz:.0f}", "Å"), select("Reduce", "sum / volume"), gap=8), gap=8)))
    return frame(viewport(sc.svg(), f"Slab z = {zc - w / 2:.0f}–{zc + w / 2:.0f} Å · outside dimmed") + bottom, right, "CAPS — slice & binning", "Analyze › Visualize › Slice & binning",
                 (f"<span>{len(inside)} of {len(recs)} atoms in the slab</span>", f"<span>mean ρ {tot:.3f} g/cm³</span>"))


# ---------------------------------------------------------------- Expression selection
def expression_select():
    cell, recs = cell_recs()
    exprs = [("Type == 2", lambda r: r["type"] == 2), ("Type == 2 && Position.Z > 13", lambda r: r["type"] == 2 and r["p"][2] > 13),
             ("Charge < -0.05", lambda r: r["q"] < -0.05), ("MoleculeIdentifier == 3", lambda r: r["mol"] == 3),
             ("Element == \"H\" && Position.X < 10", lambda r: r["e"] == "H" and r["p"][0] < 10), ("Monomer <= 2 || Monomer >= 7", lambda r: r["unit"] <= 2 or r["unit"] >= 7)]
    counts = [(e, sum(1 for r in recs if f(r))) for e, f in exprs]
    active = exprs[1]
    sel = {r["id"] for r in recs if active[1](r)}
    sc = Scene("ex", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.25)
    draw(sc, cell, recs, lambda r: ERR if r["id"] in sel else "#2A3036")
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    et = table(["Expression", "Selected"], [[f'<span style="color: {SEL}">{esc(e)}</span>', f"{c}"] for e, c in counts], ["78%", "22%"], mono_cols=(0, 1), align_right=(1,), fs=11.5, rowh=26, hl={1})
    props = ["Position.X · .Y · .Z", "Type · Element", "Charge", "MoleculeIdentifier", "Monomer", "ParticleIdentifier", "Selection", "any imported column"]
    pl_ = "".join(f'<span style="font-family: {MONO}; font-size: 11px; padding: 3px 7px; border-radius: 4px; background: {BG0}; border: 1px solid {LINE}">{esc(p)}</span>' for p in props)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex: 1 1 0; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Saved expressions · counts on this frame</h3>{et}</div>'
              f'<div style="width: 320px; flex-shrink: 0; display: flex; flex-direction: column; gap: 8px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Properties you can use</h3>'
              f'<div style="display: flex; flex-wrap: wrap; gap: 5px">{pl_}</div><span style="font-size: 11.5px; color: {DIM}">Operators: == != &lt; &gt; &amp;&amp; || ! and math functions. Same language in Studio query select.</span></div></div>')
    edit = (f'<div style="display: flex; flex-direction: column; gap: 4px"><label style="font-size: 11.5px; color: {MUTED}">Expression</label>'
            f'<div style="min-height: 56px; padding: 8px 10px; background: {BG0}; border: 1px solid {ACC}; border-radius: 6px; font-family: {MONO}; font-size: 12.5px; line-height: 1.5">'
            f'<span style="color: {SEL}">Type</span> == <span style="color: {ACC}">2</span> <span style="color: {MUTED}">&amp;&amp;</span> <span style="color: {SEL}">Position.Z</span> &gt; <span style="color: {ACC}">13</span></div></div>')
    right = inspector("Expression selection", section("Select particles", col(edit, kv("Selected", f"{len(sel)} of {len(recs)}", vcol=OK), kv("Type 2 means", "ca · aromatic carbon"), gap=8))
                      + section("Then", col(select("Next step", "Assign colour"), toggle("Show selection count in legend", True), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Counts are evaluated on the 1 300 atoms of this frame. Expressions re-evaluate on every frame.</div>')
    return frame(viewport(sc.svg(), f"Selected in red · {len(sel)} atoms") + bottom, right, "CAPS — expression selection", "Analyze › Visualize › Expression selection",
                 (f"<span>{len(sel)} selected</span>", "<span>re-evaluated per frame</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage18/project", exist_ok=True)
    for name, fn in (("ClusterAnalysis", clusters), ("RdfCoordination", rdf_board), ("SliceBinning", slice_bin), ("ExpressionSelect", expression_select)):
        h = fn(); open(f"stage18/project/{name}.dc.html", "w").write(h); print(name, len(h))
