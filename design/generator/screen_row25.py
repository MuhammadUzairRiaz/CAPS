"""Row 25 — structure: voids, Voronoi volumes (grid estimate), local density field, compare cells, pipeline groups, figure bundle."""
import math, hashlib, json
from lib import *
from screen_app import card
from display import ps_cell
from screen_row20 import annotate, draw, _shell, CHAIN, VIRIDIS, ramp
from screen_row21 import cell_recs, mi, inspector, viewport, frame, CW, VH, BOTTOM, L, MASS, NA
from screen_row22 import mol_shape, chain_col, H3, FULL
from screen_row23 import heading, code, backbone_dihedrals, MONO_BLOCK
from screen_row24 import contacts, n_clusters

NG = 30                 # grid points per axis
GS = L / NG             # 1.1 Å
VDW = {"C": 1.70, "H": 1.20}   # Bondi
RS = 6.0                # search radius for nearest atom

_F = {}


def gfield(seed=21):
    """Per grid point: (distance to nearest vdW surface, index of nearest atom by centre distance)."""
    if seed in _F: return _F[seed]
    if seed == 21:
        cell, recs = cell_recs()
    else:
        cell = ps_cell(10, 8, L, seed=seed); recs = annotate(cell)
    P = [tuple(v % L for v in r["p"]) for r in recs]
    nc = int(L // 3.3); cs = L / nc
    grid = {}
    for i, p in enumerate(P):
        grid.setdefault(tuple(int(p[k] // cs) % nc for k in range(3)), []).append(i)
    reach = int(math.ceil(RS / cs))
    surf, near = [], []
    for gx in range(NG):
        for gy in range(NG):
            for gz in range(NG):
                q = ((gx + 0.5) * GS, (gy + 0.5) * GS, (gz + 0.5) * GS)
                kx, ky, kz = (int(q[k] // cs) % nc for k in range(3))
                best_s, best_d, bi = 99.0, 99.0, -1
                for dx in range(-reach, reach + 1):
                    for dy in range(-reach, reach + 1):
                        for dz in range(-reach, reach + 1):
                            for j in grid.get(((kx + dx) % nc, (ky + dy) % nc, (kz + dz) % nc), ()):
                                d2 = 0.0
                                for k in range(3):
                                    x = P[j][k] - q[k]; x -= L * round(x / L); d2 += x * x
                                d = math.sqrt(d2)
                                s = d - VDW[recs[j]["e"]]
                                if s < best_s: best_s = s
                                if d < best_d: best_d, bi = d, j
                surf.append(best_s); near.append(bi)
    _F[seed] = (cell, recs, P, surf, near)
    return _F[seed]


def gidx(x, y, z):
    return (x % NG) * NG * NG + (y % NG) * NG + (z % NG)


def void_components(surf, probe):
    on = [s > probe for s in surf]
    lab = [-1] * len(on)
    comps = []
    for start in range(len(on)):
        if not on[start] or lab[start] >= 0: continue
        cid = len(comps); stack = [start]; lab[start] = cid; n = 0; pts = []
        while stack:
            c = stack.pop(); n += 1; pts.append(c)
            x, y, z = c // (NG * NG), (c // NG) % NG, c % NG
            for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                j = gidx(x + dx, y + dy, z + dz)
                if on[j] and lab[j] < 0:
                    lab[j] = cid; stack.append(j)
        comps.append(pts)
    comps.sort(key=len, reverse=True)
    return comps


def gpos(c):
    return ((c // (NG * NG) + 0.5) * GS, ((c // NG) % NG + 0.5) * GS, (c % NG + 0.5) * GS)


# ---------------------------------------------------------------- Voids
def voids():
    cell, recs, P, surf, near = gfield()
    probe = 1.4
    comps = void_components(surf, probe)
    sweep = []
    for pr in (0.5, 1.0, 1.4, 2.0, 3.0):
        cc = void_components(surf, pr)
        sweep.append((pr, sum(len(c) for c in cc) / NG ** 3, len(cc), len(cc[0]) / max(1, sum(len(c) for c in cc)) if cc else 0))
    vox = GS ** 3
    sc = Scene("vd", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.0, center=(L / 2, L / 2, L / 2))
    for m, (atoms, bonds, bb) in enumerate(cell):
        sc.add_tube([atoms[j]["p"] for j in bb], mix(CHAIN[m % 10], BG0, 0.55), 0.35)
    cols = [SEL, ACC, "#C79BE8", OK, "#DE775D", "#E9ECEF"]
    dots = []
    for k, cmp in enumerate(comps[:6]):
        for c in cmp[::(5 if k == 0 else 1)]:
            dots.append({"e": "H", "p": gpos(c), "c": cols[k]})
    sc.add_atoms(dots, [])
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    st = table(["Probe r (Å)", "Accessible", "Voids", "Largest share"], [[f"{a:.1f}", f"{b * 100:.1f} %", str(c), f"{d * 100:.0f} %"] for a, b, c, d in sweep],
               ["25%", "27%", "20%", "28%"], mono_cols=(0, 1, 2, 3), align_right=(1, 2, 3), fs=11.5, rowh=24, hl={2})
    vt = table(["Void", "Volume (Å³)", "Points"], [[f'<span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 10px; height: 10px; border-radius: 50%; background: {cols[k]}"></span>{k + 1}</span>', f"{len(c) * vox:,.0f}".replace(",", " "), str(len(c))] for k, c in enumerate(comps[:5])],
               ["26%", "40%", "34%"], mono_cols=(1, 2), align_right=(1, 2), fs=11.5, rowh=24)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="width: 470px; flex-shrink: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Probe sweep</h3>{st}</div>'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Voids at {probe} Å · {len(comps)} found</h3>{vt}</div></div>')
    right = inspector("Voids & pores", section("Probe", col(row(field("Radius", f"{probe}", "Å"), select("Radii", "Bondi vdW"), gap=8), row(field("Grid", f"{GS:.2f}", "Å"), select("Connect", "faces · periodic"), gap=8), gap=8))
                      + section("Show", col(toggle("Void points, coloured by void", True), toggle("Surface mesh", False), toggle("Only voids > 50 Å³", False), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">A grid point is free when a sphere of the probe radius fits there without touching any atom\'s vdW sphere. This loosely packed, built cell is mostly one connected void; an equilibrated melt breaks it into many small ones.</div>')
    acc = sweep[2][1]
    return frame(viewport(sc.svg(), f"Void points · probe {probe} Å · largest void thinned to every 5th point", "Backbone") + bottom, right, "CAPS — voids and pores", "Analyze › Visualize › Voids",
                 (f"<span>{acc * 100:.1f} % of the cell accessible to a {probe} Å probe</span>", f"<span>grid {NG}³ · {GS:.2f} Å · periodic</span>"))


# ---------------------------------------------------------------- Voronoi volumes (grid estimate)
def voronoi():
    cell, recs, P, surf, near = gfield()
    vox = GS ** 3
    vol = [0.0] * len(recs)
    for j in near:
        vol[j] += vox
    by = {}
    for r, v in zip(recs, vol):
        by.setdefault(r["gaff"], []).append(v)
    hi = sorted(vol)[int(0.90 * len(vol))]
    sc = Scene("vo", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.2)
    it = iter(vol)
    draw(sc, cell, recs, lambda r: ramp(VIRIDIS, min(1.0, vol[r["id"] - 1] / hi)))
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    rows = []
    for g in ("c3", "ca", "hc", "ha"):
        v = by[g]
        rows.append([g, str(len(v)), f"{sum(v) / len(v):.1f}", f"{sorted(v)[len(v) // 2]:.1f}", f"{max(v):.1f}"])
    t = table(["Type", "Atoms", "mean (Å³)", "median (Å³)", "max (Å³)"], rows, ["16%", "18%", "22%", "22%", "22%"], mono_cols=(0, 1, 2, 3, 4), align_right=(1, 2, 3, 4), fs=12, rowh=26)
    tot = sum(vol)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Cell volume per atom, by type</h3>{t}</div>'
              f'<div style="width: 300px; flex-shrink: 0; display: flex; flex-direction: column; gap: 8px; padding-top: 18px">{kv("Sum of cells", f"{tot:,.0f} Å³".replace(",", " "))}{kv("Box volume", f"{L ** 3:,.0f} Å³".replace(",", " "))}'
              f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Estimate: each {GS:.2f} Å grid point goes to its nearest atom. Exact Voronoi faces (and radical Voronoi weighted by radius) are the next option down.</span></div></div>')
    right = inspector("Voronoi volumes", section("Method", col(seg(["Grid estimate", "Exact", "Radical"], "Grid estimate", full=True), row(field("Grid", f"{GS:.2f}", "Å"), select("Periodic", "x y z"), gap=8), gap=8))
                      + section("Output", col(kv("Property", "AtomicVolume"), toggle("Coordination (face count)", False), toggle("Colour by volume", True), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Atoms beside voids get large cells, so the colour shows where the packing is loose. The map is clipped at the 90th percentile ({hi:.0f} Å³).</div>')
    return frame(viewport(sc.svg(), f"Colour: cell volume · 0 → {hi:.0f} Å³") + bottom, right, "CAPS — Voronoi volumes", "Analyze › Visualize › Voronoi",
                 (f"<span>{len(recs)} cells · sum {tot:,.0f} Å³ = box</span>".replace(",", " "), f"<span>grid estimate · {GS:.2f} Å</span>"))


# ---------------------------------------------------------------- Local density field
def density_field():
    cell, recs, P, surf, near = gfield()
    sig = 1.5
    rho = [0.0] * NG ** 3
    norm_ = 1.0 / ((2 * math.pi) ** 1.5 * sig ** 3)
    reach = int(math.ceil(3 * sig / GS))
    for r, p in zip(recs, P):
        m = MASS[r["e"]]
        cx, cy, cz = (int(p[k] // GS) for k in range(3))
        for dx in range(-reach, reach + 1):
            for dy in range(-reach, reach + 1):
                for dz in range(-reach, reach + 1):
                    q = ((cx + dx + 0.5) * GS, (cy + dy + 0.5) * GS, (cz + dz + 0.5) * GS)
                    d2 = sum((q[k] - p[k]) ** 2 for k in range(3))
                    rho[gidx(cx + dx, cy + dy, cz + dz)] += m * norm_ * math.exp(-d2 / (2 * sig * sig))
    conv = 1.0 / NA * 1e24      # g/mol/Å³ → g/cm³
    rho = [v * conv for v in rho]
    mean = sum(rho) / len(rho)
    bulk = sum(MASS[r["e"]] for r in recs) / NA / (L * 1e-8) ** 3
    z = NG // 2
    sl = [[rho[gidx(x, y, z)] for x in range(NG)] for y in range(NG)]
    vmax = 1.2
    cellpx = 15
    rects = "".join(f'<rect x="{x * cellpx}" y="{(NG - 1 - y) * cellpx}" width="{cellpx}" height="{cellpx}" fill="{ramp(VIRIDIS, min(1.0, sl[y][x] / vmax))}"></rect>' for y in range(NG) for x in range(NG))
    heat = (f'<svg width="{NG * cellpx}" height="{NG * cellpx}" viewBox="0 0 {NG * cellpx} {NG * cellpx}" role="img" aria-label="Density slice at z = {(z + 0.5) * GS:.1f} Å" style="display: block; border: 1px solid {LINE}">{rects}</svg>')
    stops = "".join(f'<stop offset="{t / 8:.3f}" stop-color="{ramp(VIRIDIS, t / 8)}"></stop>' for t in range(9))
    legend = (f'<svg width="{NG * cellpx}" height="34" viewBox="0 0 {NG * cellpx} 34" role="img" aria-label="Colour scale 0 to {vmax} g/cm³" style="display: block"><defs><linearGradient id="dfl">{stops}</linearGradient></defs>'
              f'<rect x="0" y="0" width="{NG * cellpx}" height="10" fill="url(#dfl)"></rect><text x="0" y="26" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">0</text>'
              f'<text x="{NG * cellpx}" y="26" text-anchor="end" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">≥ {vmax} g/cm³</text></svg>')
    nb = 30
    hc = [0] * nb
    for v in rho: hc[min(nb - 1, int(v / 1.5 * nb))] += 1
    mx = max(hc)
    hw, hh = 460, 250
    bars = "".join(f'<rect x="{40 + i * (hw - 50) / nb:.1f}" y="{hh - 26 - c / mx * (hh - 44):.1f}" width="{(hw - 50) / nb - 1.5:.1f}" height="{c / mx * (hh - 44):.1f}" fill="{ramp(VIRIDIS, min(1.0, (i + 0.5) / nb * 1.5 / vmax))}"></rect>' for i, c in enumerate(hc) if c)
    xt = "".join(f'<text x="{40 + v / 1.5 * (hw - 50):.1f}" y="{hh - 10}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{v:g}</text>' for v in (0, 0.5, 1.0))
    mk = f'<line x1="{40 + mean / 1.5 * (hw - 50):.1f}" y1="14" x2="{40 + mean / 1.5 * (hw - 50):.1f}" y2="{hh - 26}" stroke="{ACC}" stroke-dasharray="4 3"></line><text x="{44 + mean / 1.5 * (hw - 50):.1f}" y="24" font-size="10" font-family="IBM Plex Mono" fill="{ACC}">mean {mean:.3f}</text>'
    hist = (f'<svg width="{hw}" height="{hh}" viewBox="0 0 {hw} {hh}" role="img" aria-label="Histogram of local density" style="display: block">{bars}{xt}{mk}'
            f'<line x1="40" y1="{hh - 26}" x2="{hw - 10}" y2="{hh - 26}" stroke="{LINE}"></line><text x="40" y="11" font-size="10" fill="{DIM}">grid points by local density (g/cm³) · last bar ≥ 1.45</text></svg>')
    empty = sum(1 for v in rho if v < 0.05) / len(rho)
    content = (f'<div style="flex-grow: 1; display: flex; gap: 22px; padding: 16px; min-height: 0; overflow: hidden">'
               f'<div style="flex-shrink: 0; display: flex; flex-direction: column; gap: 8px"><h3 style="{H3}">Slice z = {(z + 0.5) * GS:.1f} Å · x → right, y → up</h3>{heat}{legend}</div>'
               f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px"><h3 style="{H3}">Distribution</h3>{hist}'
               f'<div style="display: flex; flex-direction: column; gap: 6px; max-width: 420px">{kv("Mean of field", f"{mean:.3f} g/cm³")}{kv("Mass / box volume", f"{bulk:.3f} g/cm³")}{kv("Points below 0.05 g/cm³", f"{empty * 100:.0f} %")}</div>'
               f'<span style="font-size: 12px; color: {MUTED}; line-height: 1.6; max-width: 460px">The field mean equals the bulk density, which checks the smoothing conserves mass. The wide spread and the empty points are what a loosely packed, un-equilibrated build looks like.</span></div></div>')
    right = inspector("Density field", section("Field", col(select("Quantity", "mass density"), row(field("Grid", f"{GS:.2f}", "Å"), field("Smoothing σ", f"{sig}", "Å"), gap=8), gap=8))
                      + section("View", col(seg(["Slice", "Isosurface", "Volume"], "Slice", full=True), row(select("Normal", "z"), field("Position", f"{(z + 0.5) * GS:.1f}", "Å"), gap=8), gap=8))
                      + section("Export", col(select("As", "CUBE · VTK · NumPy"), gap=8)))
    return frame(content, right, "CAPS — density field", "Analyze › Visualize › Density field",
                 (f"<span>{NG}³ grid · σ = {sig} Å · mean {mean:.3f} g/cm³</span>", "<span>Gaussian-smoothed atomic masses, periodic</span>"))


# ---------------------------------------------------------------- Compare two cells
def rdf(recs, P, rmax=12.0, dr=0.4):
    C = [i for i, r in enumerate(recs) if r["e"] == "C"]
    nb = int(rmax / dr); h = [0] * nb
    for a in range(len(C)):
        for b in range(a + 1, len(C)):
            i, j = C[a], C[b]
            if recs[i]["mol"] == recs[j]["mol"]: continue
            d = mi(P[i], P[j])
            if d < rmax: h[int(d / dr)] += 2
    n = len(C); rho = n / L ** 3
    return [((k + 0.5) * dr, h[k] / (n * rho * 4 / 3 * math.pi * (((k + 1) * dr) ** 3 - (k * dr) ** 3))) for k in range(nb)]


def metrics(seed):
    cell, recs, P, surf, near = gfield(seed)
    sh = mol_shape(recs)
    mols = sorted(sh)
    pairs = contacts(P, recs, 4.5)
    ncl = n_clusters([p for p in pairs if mi(P[p[0]], P[p[1]]) < 3.3], recs, mols)
    dh = backbone_dihedrals(cell)
    comps = void_components(surf, 1.4)
    return dict(cell=cell, recs=recs, P=P, rg=sum(s["rg"] for s in sh.values()) / len(sh), k2=sum(s["k2"] for s in sh.values()) / len(sh),
                con=len(pairs), ncl=ncl, ft=sum(abs(p) > 120 for _, _, p in dh) / len(dh), acc=sum(len(c) for c in comps) / NG ** 3, nv=len(comps), g=rdf(recs, P))


def compare():
    A, B = metrics(21), metrics(23)
    hw = (CW - 6 - 6) // 2
    tiles = []
    for tag, M in (("A · seed 21", A), ("B · seed 23", B)):
        sc = Scene("cm" + tag[0], hw, 360, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.1, center=(L / 2, L / 2, L / 2))
        draw(sc, M["cell"], M["recs"], chain_col)
        sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
        tiles.append(f'<div style="position: relative; width: {hw}px; height: 360px; background: {BG0}; overflow: hidden">{sc.svg(overlay=label_pill(10, 18, tag, TEXT, fs=11))}</div>')
    split = f'<div style="flex-shrink: 0; display: flex; gap: 6px; background: {LINE}">{"".join(tiles)}</div>'
    gm = max(max(v for _, v in A["g"]), max(v for _, v in B["g"]))
    gm = math.ceil(gm * 2) / 2
    pl = plot(470, 250, [(A["g"], SEL, 2, None), (B["g"], ACC, 2, "5 3"), ([(0, 1), (12, 1)], DIM, 1, "2 4")], (0, 12), (0, gm), [0, 3, 6, 9, 12], [0, gm / 2, gm], "r (Å)", "g(r) · C–C, other molecules", pad=(40, 14, 12, 30))
    rows = []
    for lab, k, fmt in (("mean Rg (Å)", "rg", "{:.2f}"), ("mean κ²", "k2", "{:.3f}"), ("trans fraction", "ft", "{:.3f}"), ("contacts < 4.5 Å", "con", "{}"), ("clusters at 3.3 Å", "ncl", "{}"), ("accessible (1.4 Å probe)", "acc", "{:.1%}"), ("voids (1.4 Å probe)", "nv", "{}")):
        a, b = A[k], B[k]
        diff = b - a
        ds = (f"{diff:+.1%}" if k == "acc" else (f"{diff:+d}" if isinstance(diff, int) else f"{diff:+.3f}"))
        rows.append([lab, fmt.format(a), fmt.format(b), ds])
    t = table(["Metric", "A", "B", "B − A"], rows, ["40%", "20%", "20%", "20%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=12, rowh=25)
    bottom = (f'<div style="flex-grow: 1; min-height: 0; display: flex; gap: 18px; padding: 12px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-shrink: 0; display: flex; flex-direction: column; gap: 6px">{pl}<div style="display: flex; gap: 16px; font-size: 12px; color: {MUTED}">'
              f'<span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 16px; height: 2px; background: {SEL}"></span>A</span><span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 16px; height: 0; border-top: 2px dashed {ACC}"></span>B</span></div></div>'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Same pipeline, two inputs</h3>{t}</div></div>')
    right = inspector("Compare", section("Inputs", col(select("A", "cells/seed21/PS_melt.data"), select("B", "cells/seed23/PS_melt.data"), toggle("Link cameras", True), gap=8))
                      + section("Pipeline", col(select("Applied to both", "ps_structure_report"), toggle("Show B − A column", True), toggle("Overlay plots", True), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Two builds with the same recipe and different seeds. Differences here are seed-to-seed scatter, the size of error you should expect from one build.</div>')
    return frame(split + bottom, right, "CAPS — compare cells", "Analyze › Visualize › Compare",
                 ("<span>A and B: 10 × PS DP 8 in 33 Å, seeds 21 and 23</span>", "<span>same pipeline hash on both sides</span>"))


# ---------------------------------------------------------------- Pipeline groups & branches
def pipeline_groups():
    def node(title, sub, colour=MUTED, on=True, indent=0, sel=False):
        return (f'<div style="display: flex; align-items: center; gap: 10px; margin-left: {indent * 22}px; padding: 8px 10px; border-radius: 6px; background: {BG3 if sel else BG2}; border: 1px solid {ACC if sel else LINE}">'
                f'<span style="width: 4px; align-self: stretch; border-radius: 2px; background: {colour}"></span>'
                f'<div style="flex-grow: 1; display: flex; flex-direction: column"><span style="font-size: 12.5px{"; color: " + DIM if not on else ""}">{title}</span><span style="font-size: 11px; color: {DIM}; font-family: {MONO}">{sub}</span></div>'
                f'{icon("eye" if on else "xcircle", 14, MUTED if on else DIM)}</div>')
    def grp(title, n, colour, open_=True, indent=0):
        return (f'<div style="display: flex; align-items: center; gap: 8px; margin-left: {indent * 22}px; padding: 6px 4px; font-size: 12px; font-weight: 600; color: {colour}">'
                f'{icon("chev" if open_ else "chevr", 12, colour)}{icon("layers", 14, colour)}<span>{title}</span><span style="font-family: {MONO}; font-weight: 400; color: {DIM}">{n} steps</span></div>')
    tree = col(
        node("Source · PS_melt.lammpstrj", "500 frames · topology PS_melt.data", TEXT),
        grp("Prepare", 3, SEL), node("Unwrap", "image flags", SEL, indent=1), node("Create bonds", "topology + cutoff check", SEL, indent=1), node("Compute property", "DistanceToCOM", SEL, indent=1),
        grp("Structure", 3, ACC), node("Cluster analysis", "3.3 Å · molecules", ACC, indent=1), node("Coordination & RDF", "C–C · 12 Å", ACC, indent=1, sel=True), node("Voids", "probe 1.4 Å", ACC, on=False, indent=1),
        grp("Look", 2, "#C79BE8", open_=False), gap=6)
    branch = col(
        f'<div style="font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Branches after “Structure”</div>',
        node("Branch · Figure", "colour by DistanceToCOM → render 1920 × 1080", "#C79BE8"),
        node("Branch · Numbers", "tables → results.csv, rdf.svg", OK),
        node("Branch · Chain view", "backbone only, colour by molecule", "#DE775D"), gap=6)
    flow = (f'<svg width="500" height="223" viewBox="0 0 560 250" role="img" aria-label="Pipeline graph: source, three groups, three branches" style="display: block">'
            + "".join(f'<rect x="{x}" y="{y}" width="{w}" height="34" rx="6" fill="{BG2}" stroke="{c}"></rect><text x="{x + w / 2}" y="{y + 21}" text-anchor="middle" font-size="11.5" fill="{TEXT}" font-family="IBM Plex Sans">{t}</text>'
                      for x, y, w, c, t in ((10, 108, 90, TEXT, "Source"), (130, 108, 90, SEL, "Prepare"), (250, 108, 90, ACC, "Structure"), (400, 30, 150, "#C79BE8", "Figure"), (400, 108, 150, OK, "Numbers"), (400, 186, 150, "#DE775D", "Chain view")))
            + "".join(f'<path d="M{a} {b} C {a + 30} {b}, {c - 30} {d}, {c} {d}" fill="none" stroke="{DIM}" stroke-width="1.4"></path>' for a, b, c, d in ((100, 125, 130, 125), (220, 125, 250, 125), (340, 125, 400, 47), (340, 125, 400, 125), (340, 125, 400, 203)))
            + '</svg>')
    left = card("Pipeline", tree, row(btn("Group", small=True, ic="layers"), btn("Branch", small=True, ic="link"), gap=6), 14, "width: 460px; flex-shrink: 0")
    mid = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px">'
           + card("Graph", flow, chip("1 source · 3 groups · 3 branches", MUTED, BG2), 14)
           + card("Branches", branch, "", 14) + '</div>')
    right = (f'<div style="width: 300px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Group · Structure", col(field("Name", "Structure", mono=False), toggle("Enabled", True), toggle("Collapse in list", False), kv("Cached frames", "[measure]"), gap=8), "", 12)
             + card("Why groups and branches", f'<span style="font-size: 12px; color: {MUTED}; line-height: 1.6">Shared steps run once per frame; each branch reuses their cached output. Turning a group off skips all its steps; the saved YAML keeps the grouping.</span>', "", 12)
             + '</div>')
    content = (heading("Pipeline · groups & branches", "Organise long pipelines; fork one result into several views without recomputing", btn("Save pipeline", ic="save"))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 14px 22px; min-height: 0; overflow: hidden">{left}{mid}{right}</div>')
    return _shell("CAPS — pipeline groups", f'{icon("chevr", 12, DIM)}<span>Analyze › Visualize › Pipeline</span>', content,
                  ("<span>3 groups · 8 steps · 3 branches</span>", "<span>disabled steps stay in the file, marked off</span>"))


# ---------------------------------------------------------------- Figure bundle
def figure_bundle():
    A = metrics(21)
    csv = "r_A,g_r\n" + "".join(f"{r:.2f},{g:.4f}\n" for r, g in A["g"])
    prov = {"caps": "0.1.0", "input": {"file": "cells/seed21/PS_melt.data", "atoms": len(A["recs"])}, "box_A": [L, L, L],
            "build": {"polymer": "polystyrene", "dp": 8, "chains": 10, "tacticity": "atactic", "seed": 21},
            "pipeline": {"name": "ps_structure_report", "steps": 7}, "figure": {"size_px": [1920, 1080], "panels": ["render", "rdf"]}}
    pj = "{\n" + ",\n".join(f'  "{k}": {json.dumps(v)}' for k, v in prov.items()) + "\n}"
    h = lambda s: hashlib.sha256(s.encode()).hexdigest()
    files = [("figure.svg", "vector figure", "[measure]", "—"), ("figure.png", "1920 × 1080 raster", "[measure]", "—"),
             ("rdf.csv", f"{len(A['g'])} rows · real values", f"{len(csv.encode()) / 1024:.1f} KB", h(csv)[:12]),
             ("provenance.json", "inputs, build, pipeline", f"{len(pj.encode()) / 1024:.1f} KB", h(pj)[:12]),
             ("ps_structure_report.caps-pipeline.yaml", "the pipeline as run", "[measure]", "—"), ("README.txt", "how to reproduce", "[measure]", "—")]
    ft = table(["File", "Content", "Size", "sha256"], [[a, b, c, d] for a, b, c, d in files], ["36%", "30%", "14%", "20%"], mono_cols=(0, 2, 3), align_right=(2,), fs=12, rowh=30)
    pv = code(pj.split("\n"))
    csvp = code(csv.split("\n")[:8] + ["…"])
    left = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("PS_rdf_figure.caps-bundle.zip", col(ft, f'<span style="font-size: 11.5px; color: {DIM}">Hashes are shown for files whose content is already fixed on this screen; the rest are written on export.</span>', gap=8), chip("6 files", MUTED, BG2, True), 14, "flex-shrink: 0")
            + f'<div style="display: flex; gap: 12px; min-height: 0">' + card("provenance.json", pv, "", 12, "flex-grow: 1; min-width: 0") + card("rdf.csv", csvp, "", 12, "width: 220px; flex-shrink: 0") + '</div></div>')
    right = (f'<div style="width: 320px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Bundle", col(toggle("Include input structure", False), toggle("Include pipeline", True), toggle("Include data behind every plot", True), toggle("Include README", True), gap=8), "", 12)
             + card("Reproduce", f'<div style="{MONO_BLOCK}">caps reproduce PS_rdf_figure.caps-bundle.zip</div><span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Rebuilds from the recorded seed and pipeline, then compares hashes of the data files.</span>', "", 12)
             + '</div>')
    content = (heading("Export · figure bundle", "A figure that carries its data, pipeline and inputs", btn("Cancel") + btn("Export bundle", True, "download"))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 14px 22px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return _shell("CAPS — figure bundle", f'{icon("chevr", 12, DIM)}<span>Export › Figure bundle</span>', content,
                  (f"<span>rdf.csv sha256 {h(csv)[:16]}</span>", "<span>hashes computed from the text shown</span>"))


BOARDS = (("VoidAnalysis", voids), ("VoronoiVolumes", voronoi), ("DensityField", density_field),
          ("CompareCells", compare), ("PipelineGroups", pipeline_groups), ("FigureBundle", figure_bundle))

if __name__ == "__main__":
    import os, sys
    os.makedirs("stage25/project", exist_ok=True)
    only = sys.argv[1:]
    for name, fn in BOARDS:
        if only and name not in only: continue
        h = fn(); open(f"stage25/project/{name}.dc.html", "w").write(h); print(name, len(h), flush=True)
