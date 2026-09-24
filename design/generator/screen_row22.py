"""Visualize row 22: four viewports, wrap/unwrap, molecule shape, histogram & scatter, data inspector, render & overlays."""
import math
from lib import *
from screen_row20 import draw, CHAIN, VIRIDIS, ramp
from screen_row21 import cell_recs, mi, inspector, viewport, frame, CW, VH, BOTTOM, L, MASS, NA

FULL = VH + BOTTOM
H3 = f'margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}'


def chain_col(r):
    c = CHAIN[r["chain"] % 10]
    return c if r["e"] != "H" else mix(c, "#FFFFFF", 0.55)


def wrapped(p):
    return tuple(v % L for v in p)


def eig3(A):
    """Jacobi eigen-decomposition of a symmetric 3×3: ascending (value, vector) pairs."""
    A = [r[:] for r in A]
    V = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    for _ in range(80):
        p, q = max(((i, j) for i in range(3) for j in range(i + 1, 3)), key=lambda t: abs(A[t[0]][t[1]]))
        if abs(A[p][q]) < 1e-14:
            break
        th = 0.5 * math.atan2(2 * A[p][q], A[q][q] - A[p][p])
        c, s = math.cos(th), math.sin(th)
        for k in range(3):
            A[k][p], A[k][q] = c * A[k][p] - s * A[k][q], s * A[k][p] + c * A[k][q]
        for k in range(3):
            A[p][k], A[q][k] = c * A[p][k] - s * A[q][k], s * A[p][k] + c * A[q][k]
        for k in range(3):
            V[k][p], V[k][q] = c * V[k][p] - s * V[k][q], s * V[k][p] + c * V[k][q]
    return sorted((A[i][i], (V[0][i], V[1][i], V[2][i])) for i in range(3))


def mol_shape(recs):
    out = {}
    for m in sorted({r["mol"] for r in recs}):
        rs = [r for r in recs if r["mol"] == m]
        M = sum(MASS[r["e"]] for r in rs)
        com = tuple(sum(MASS[r["e"]] * r["p"][k] for r in rs) / M for k in range(3))
        G = [[sum(MASS[r["e"]] * (r["p"][i] - com[i]) * (r["p"][j] - com[j]) for r in rs) / M for j in range(3)] for i in range(3)]
        ev = eig3(G)
        l1, l2, l3 = (e[0] for e in ev)
        rg2 = l1 + l2 + l3
        bb = [r for r in rs if r["e"] == "C" and r["j"] < 16]
        ree = dist(bb[0]["p"], bb[-1]["p"])
        out[m] = dict(com=com, ev=ev, rg=math.sqrt(rg2), b=l3 - 0.5 * (l1 + l2), c=l2 - l1,
                      k2=1 - 3 * (l1 * l2 + l2 * l3 + l3 * l1) / rg2 ** 2, ree=ree, n=len(rs), M=M)
    return out


# ---------------------------------------------------------------- Four viewports
def four_viewports():
    cell, recs = cell_recs()
    gw, gh = (CW - 6 - 6) // 2, (FULL - 6) // 2
    views = [("Top", 0.0, math.pi / 2, 0.0, "x → right · z → down"), ("Front", 0.0, 0.0, 0.0, "x → right · y → up"),
             ("Left", math.pi / 2, 0.0, 0.0, "z → right · y → up"), ("Perspective", 0.55, 0.4, 0.3, "orbit camera")]
    tiles = []
    for i, (name, yaw, pitch, persp, axes) in enumerate(views):
        sc = Scene(f"v{i}", gw, gh, yaw=yaw, pitch=pitch, persp=persp, fog=0.3 if persp else 0.15, atom_k=1.1, center=(L / 2, L / 2, L / 2))
        draw(sc, cell, recs, chain_col)
        sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
        active = name == "Perspective"
        ov = label_pill(10, 18, name + (" · ortho" if not persp else ""), ACC if active else TEXT, fs=11) + gizmo(gw - 34, gh - 34, yaw, pitch, 20)
        tiles.append(f'<div style="position: relative; width: {gw}px; height: {gh}px; background: {BG0}; outline: {"2px solid " + ACC if active else "1px solid " + LINE}; outline-offset: -1px; overflow: hidden">'
                     f'{sc.svg(overlay=ov)}<span style="position: absolute; left: 10px; bottom: 8px; font-family: {MONO}; font-size: 10.5px; color: {DIM}">{axes}</span></div>')
    grid = f'<div style="flex-grow: 1; min-height: 0; display: grid; grid-template-columns: {gw}px {gw}px; gap: 6px; padding: 0 3px; background: {LINE}">{"".join(tiles)}</div>'
    right = inspector("Viewports", section("Layout", col(seg(["1", "2 × 1", "2 × 2", "1 + 3"], "2 × 2", full=True), toggle("Link zoom across ortho views", True), toggle("Link selection", True), gap=8))
                      + section("Active: Perspective", col(select("Projection", "perspective · 35° FOV"), select("Preview", "render frame 16:9"), toggle("Show cell", True), toggle("Show axis tripod", True), gap=8))
                      + section("Maximise", col(kv("Active view", "Enter"), kv("Swap view", "right-click title"), kv("Cycle", "Tab"), gap=6))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">All four views draw the same frame: {len(recs)} atoms, {L:.0f} Å cell. Ortho views keep lengths comparable, so measure there.</div>')
    return frame(grid, right, "CAPS — four viewports", "Analyze › Visualize › Viewports",
                 (f"<span>{len(recs)} atoms · 4 views of one frame</span>", "<span>ortho: 1 Å = same length in every axis</span>"))


# ---------------------------------------------------------------- Wrap / unwrap
def wrap_unwrap():
    cell, recs = cell_recs()
    outside = [r for r in recs if any(not (0 <= v < L) for v in r["p"])]
    img = lambda p: tuple(math.floor(v / L) for v in p)
    out_ids = {r["id"] for r in outside}
    cross_mol = sorted({r["mol"] for r in outside})
    nbond_cross = 0
    wcell = []
    for atoms, bonds, bb in cell:
        wa = [dict(a, p=wrapped(a["p"])) for a in atoms]
        keep = []
        for a, b in bonds:
            if img(atoms[a]["p"]) != img(atoms[b]["p"]):
                nbond_cross += 1
            else:
                keep.append((a, b))
        wcell.append((wa, keep, bb))
    wrecs = [dict(r, p=wrapped(r["p"])) for r in recs]
    hw = (CW - 6 - 6) // 2
    tiles = []
    for tag, c, rs, note in (("Unwrapped · chains whole", cell, recs, f"{len(outside)} atoms outside the cell"),
                             ("Wrapped · atoms in cell", wcell, wrecs, f"{nbond_cross} bonds cut at the faces")):
        sc = Scene("w" + tag[0], hw, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.1, center=(L / 2, L / 2, L / 2), scale=VH / (L * 2.35))
        draw(sc, c, rs, lambda r: (ERR if r["id"] in out_ids else chain_col(r)) if c is cell else chain_col(r))
        sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), ACC, 1.3, None, 0.8)
        tiles.append(f'<div style="position: relative; width: {hw}px; height: {VH}px; background: {BG0}; overflow: hidden">{sc.svg(overlay=label_pill(10, 18, tag, TEXT, fs=11) + label_pill(10, VH - 16, note, DIM, fs=10.5))}</div>')
    split = f'<div style="flex-grow: 1; min-height: 0; display: flex; gap: 6px; background: {LINE}">{"".join(tiles)}</div>'
    sample = sorted(outside, key=lambda r: r["id"])[:5]
    t = table(["id", "mol", "unwrapped x y z (Å)", "image ix iy iz", "wrapped x y z (Å)"],
              [[str(r["id"]), str(r["mol"]), " ".join(f"{v:7.2f}" for v in r["p"]), " ".join(f"{v:2d}" for v in img(r["p"])), " ".join(f"{v:6.2f}" for v in wrapped(r["p"]))] for r in sample],
              ["8%", "8%", "32%", "20%", "32%"], mono_cols=(0, 1, 2, 3, 4), align_right=(0, 1), fs=11.5, rowh=24)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Image flags · first atoms outside the cell</h3>{t}</div>'
              f'<div style="width: 250px; flex-shrink: 0; display: flex; flex-direction: column; gap: 8px; padding-top: 18px">{kv("Atoms outside", f"{len(outside)} of {len(recs)}")}{kv("Molecules crossing", f"{len(cross_mol)} of {len(cell)}")}{kv("Bonds across faces", str(nbond_cross))}'
              f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Rg, end-to-end and clusters always use unwrapped chains.</span></div></div>')
    right = inspector("Wrap & unwrap", section("Mode", col(seg(["Wrap", "Unwrap"], "Unwrap", full=True), select("Unwrap using", "image flags ix iy iz"), toggle("Keep molecules whole", True), gap=8))
                      + section("Fallback", col(select("No image flags", "follow bonds (minimum image)"), toggle("Warn when a bond > L/2", True), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Left: this cell as built, atoms outside the box in red. Right: the same atoms folded into [0, {L:.0f}) Å; bonds that would stretch across the box are hidden, not drawn as long sticks.</div>')
    return frame(split + bottom, right, "CAPS — wrap and unwrap", "Analyze › Visualize › Wrap & unwrap",
                 (f"<span>{len(outside)} atoms outside · {len(cross_mol)} molecules cross faces</span>", "<span>periodic in x, y, z</span>"))


# ---------------------------------------------------------------- Molecule shape (gyration tensor)
def molecule_shape():
    cell, recs = cell_recs()
    sh = mol_shape(recs)
    sc = Scene("ms", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.0)
    for m, (atoms, bonds, bb) in enumerate(cell):
        sc.add_tube([atoms[j]["p"] for j in bb], mix(CHAIN[m % 10], BG0, 0.5), 0.45)
    for m, s in sh.items():
        c = CHAIN[(m - 1) % 10]
        for lam, v in s["ev"]:
            h = math.sqrt(lam) * math.sqrt(3)
            sc.add_line(add(s["com"], mul(v, -h)), add(s["com"], mul(v, h)), c, 3.0 if lam == s["ev"][2][0] else 1.8, None, 1.0)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    rows = []
    for m, s in sh.items():
        sw = f'<span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 10px; height: 10px; border-radius: 2px; background: {CHAIN[(m - 1) % 10]}"></span>{m}</span>'
        rows.append([sw, f"{s['rg']:.2f}", f"{s['ree']:.2f}", f"{s['ree'] ** 2 / s['rg'] ** 2:.2f}", f"{s['b']:.2f}", f"{s['c']:.2f}", f"{s['k2']:.3f}"])
    t = table(["mol", "Rg (Å)", "Ree (Å)", "Ree²/Rg²", "b (Å²)", "c (Å²)", "κ²"], rows[:6], ["12%", "14%", "14%", "15%", "15%", "15%", "15%"], mono_cols=(1, 2, 3, 4, 5, 6), align_right=(1, 2, 3, 4, 5, 6), fs=11.5, rowh=22)
    k2m = sum(s["k2"] for s in sh.values()) / len(sh)
    rgm = sum(s["rg"] for s in sh.values()) / len(sh)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Per molecule · 6 of {len(sh)} shown</h3>{t}</div>'
              f'<div style="width: 250px; flex-shrink: 0; display: flex; flex-direction: column; gap: 8px; padding-top: 18px">{kv("Mean Rg", f"{rgm:.2f} Å")}{kv("Mean κ²", f"{k2m:.3f}")}'
              f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">κ² is 0 for a sphere and 1 for a rod. DP 8 chains are short, so shapes are far from Gaussian coils.</span></div></div>')
    right = inspector("Molecule shape", section("Group by", col(select("Unit", "molecule"), select("Atoms", "all (mass-weighted)"), toggle("Unwrap first", True), gap=8))
                      + section("Outputs", col(toggle("Rg, Ree", True), toggle("Principal moments λ₁ ≤ λ₂ ≤ λ₃", True), toggle("Asphericity b, acylindricity c, κ²", True), gap=8))
                      + section("Display", col(select("Glyph", "principal axes · ±√(3λ)"), toggle("Dim chains", True), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Eigen-decomposition of each molecule\'s gyration tensor on this frame. The thick line is the longest axis.</div>')
    return frame(viewport(sc.svg(), "Glyphs: principal axes per molecule", "Backbone") + bottom, right, "CAPS — molecule shape", "Analyze › Visualize › Molecule shape",
                 (f"<span>{len(sh)} molecules · gyration tensor</span>", f"<span>mean Rg {rgm:.2f} Å · mean κ² {k2m:.3f}</span>"))


# ---------------------------------------------------------------- Histogram & scatter
def hist_scatter():
    cell, recs = cell_recs()
    sh = mol_shape(recs)
    C = [r for r in recs if r["e"] == "C"]
    coord = {r["id"]: 0 for r in C}
    for i, a in enumerate(C):
        for b in C[i + 1:]:
            if a["mol"] != b["mol"] and mi(a["p"], b["p"]) < 5.0:
                coord[a["id"]] += 1; coord[b["id"]] += 1
    for r in C:
        r["dcom"] = dist(r["p"], sh[r["mol"]]["com"])
    x0, y0 = 4.0, 6
    brushed = {r["id"] for r in C if r["dcom"] >= x0 and coord[r["id"]] >= y0}
    xm = math.ceil(max(r["dcom"] for r in C))
    ym = max(coord.values())
    sw, shh = 470, BOTTOM - 20
    Lp, Tp, Rp, Bp = 40, 12, 10, 26
    X = lambda v: Lp + v / xm * (sw - Lp - Rp)
    Y = lambda v: Tp + (shh - Tp - Bp) * (1 - v / ym)
    pts = "".join(f'<circle cx="{X(r["dcom"]):.1f}" cy="{Y(coord[r["id"]]):.1f}" r="2.2" fill="{SEL if r["id"] in brushed else DIM}" fill-opacity="{0.9 if r["id"] in brushed else 0.55}"></circle>' for r in C)
    xt = "".join(f'<text x="{X(v):.1f}" y="{shh - 10}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{v}</text>' for v in range(0, xm + 1, 2))
    yt = "".join(f'<text x="{Lp - 6}" y="{Y(v) + 3.5:.1f}" text-anchor="end" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{v}</text>' for v in range(0, ym + 1, max(1, ym // 4)))
    scat = (f'<svg width="{sw}" height="{shh}" viewBox="0 0 {sw} {shh}" role="img" aria-label="Scatter: coordination versus distance to centre of mass" style="display: block">'
            f'<rect x="{X(x0):.1f}" y="{Tp}" width="{X(xm) - X(x0):.1f}" height="{Y(y0) - Tp:.1f}" fill="{SEL}" fill-opacity="0.08" stroke="{SEL}" stroke-dasharray="4 3"></rect>{pts}{xt}{yt}'
            f'<line x1="{Lp}" y1="{shh - Bp}" x2="{sw - Rp}" y2="{shh - Bp}" stroke="{LINE}"></line>'
            f'<text x="{sw - Rp}" y="{shh}" text-anchor="end" font-size="10" fill="{DIM}">distance to own COM (Å)</text><text x="{Lp}" y="9" font-size="10" fill="{DIM}">C neighbours within 5 Å, other molecules</text></svg>')
    types = [("c3", "#8E959C"), ("ca", ACC), ("hc", SEL), ("ha", "#9B7AD5")]
    qs = [r["q"] for r in recs]
    lo, hi = math.floor(min(qs) * 100) / 100, math.ceil(max(qs) * 100) / 100
    nb = 26
    hw, hh = 440, BOTTOM - 20
    cnt = {g: [0] * nb for g, _ in types}
    for r in recs:
        cnt[r["gaff"]][min(nb - 1, int((r["q"] - lo) / (hi - lo) * nb))] += 1
    mx = max(sum(cnt[g][k] for g, _ in types) for k in range(nb))
    bars = []
    for k in range(nb):
        yb = hh - 24
        for g, c in types:
            v = cnt[g][k]
            if not v: continue
            hgt = v / mx * (hh - 40)
            yb -= hgt
            bars.append(f'<rect x="{36 + k * (hw - 46) / nb:.1f}" y="{yb:.1f}" width="{(hw - 46) / nb - 2:.1f}" height="{hgt:.1f}" fill="{c}"></rect>')
    ticks = "".join(f'<text x="{36 + t * (hw - 46):.1f}" y="{hh - 8}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{lo + t * (hi - lo):+.3f}</text>' for t in (0, 0.5, 1))
    leg = "".join(f'<rect x="{hw - 190 + i * 46}" y="4" width="9" height="9" fill="{c}"></rect><text x="{hw - 178 + i * 46}" y="12" font-size="10" font-family="IBM Plex Mono" fill="{MUTED}">{g}</text>' for i, (g, c) in enumerate(types))
    hist = (f'<svg width="{hw}" height="{hh}" viewBox="0 0 {hw} {hh}" role="img" aria-label="Histogram of Gasteiger charge by atom type" style="display: block">{"".join(bars)}{ticks}{leg}'
            f'<line x1="36" y1="{hh - 24}" x2="{hw - 10}" y2="{hh - 24}" stroke="{LINE}"></line><text x="36" y="12" font-size="10" fill="{DIM}">Gasteiger q (e) · stacked by type</text></svg>')
    sc = Scene("hs", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.2)
    draw(sc, cell, recs, lambda r: SEL if r["id"] in brushed else ("#39414A" if r["e"] == "H" else "#5A636C"))
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 22px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-shrink: 0">{scat}</div><div style="flex-shrink: 0">{hist}</div></div>')
    right = inspector("Scatter plot", section("Axes", col(select("x", "DistanceToCOM"), select("y", "Coordination (C, 5 Å)"), select("Points", f"carbons · {len(C)}"), gap=8))
                      + section("Brush", col(kv("Region", f"x ≥ {x0:.0f} Å · y ≥ {y0}"), kv("Selected", f"{len(brushed)} atoms"), toggle("Brush selects in viewport", True), gap=8))
                      + section("Histogram", col(select("Property", "Charge"), row(field("Bins", str(nb)), select("Stack by", "GAFF type"), gap=8), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Dragging a box on the plot selects those atoms; the viewport shows them in cyan. Values from this frame.</div>')
    return frame(viewport(sc.svg(), f"Selected by brush · {len(brushed)} carbons") + bottom, right, "CAPS — histogram and scatter", "Analyze › Visualize › Histogram & scatter",
                 (f"<span>{len(C)} carbons plotted · {len(brushed)} brushed</span>", f"<span>q from {min(qs):+.3f} to {max(qs):+.3f} e</span>"))


# ---------------------------------------------------------------- Data inspector
def data_inspector():
    cell, recs = cell_recs()
    pick = recs[40]
    ntop = VH - 90
    sc = Scene("di", CW - 6, ntop, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.2)
    draw(sc, cell, recs, chain_col)
    idx = recs.index(pick)
    sc.halo.add(idx)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    mass = sum(MASS[r["e"]] for r in recs)
    rho = mass / NA / (L * 1e-8) ** 3
    nb = sum(len(b) for _, b, _ in cell)
    qtot = sum(r["q"] for r in recs)
    rows = []
    for r in recs[36:44]:
        rows.append([str(r["id"]), str(r["mol"]), str(r["type"]), r["gaff"], f"{r['q']:+.4f}", f"{r['p'][0]:.3f}", f"{r['p'][1]:.3f}", f"{r['p'][2]:.3f}", str(r["unit"])])
    t = table(["Identifier", "Molecule", "Type", "Name", "Charge", "Position.X", "Position.Y", "Position.Z", "Unit"], rows,
              ["13%", "10%", "7%", "8%", "12%", "13%", "13%", "13%", "11%"], mono_cols=tuple(range(9)), align_right=(0, 1, 2, 4, 5, 6, 7, 8), fs=11.5, rowh=23, hl={4})
    attrs = [("Timestep", "120000"), ("SourceFrame", "120"), ("Particles", str(len(recs))), ("Bonds", str(nb)), ("Molecules", str(len(cell))),
             ("Cell volume", f"{L ** 3:,.0f} Å³".replace(",", " ")), ("Mass", f"{mass:,.1f} g/mol".replace(",", " ")), ("Density", f"{rho:.3f} g/cm³"), ("Total charge", f"{qtot:+.2e} e")]
    at = "".join(f'<div style="display: flex; justify-content: space-between; gap: 10px; padding: 3px 0; border-bottom: 1px solid {BG2}; font-size: 11.5px"><span style="color: {MUTED}">{k}</span><span style="font-family: {MONO}">{v}</span></div>' for k, v in attrs)
    panel = (f'<div style="height: {BOTTOM + 90}px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-top: 1px solid {LINE}">'
             f'<div style="display: flex; align-items: center; gap: 12px; padding: 0 12px; border-bottom: 1px solid {LINE}">{tabs(["Particles", "Bonds", "Global attributes", "Data tables", "Dimensions"], "Particles")}'
             f'<div style="flex-grow: 1"></div><div style="width: 300px">{field("Filter", "Molecule == 1", mono=True)}</div></div>'
             f'<div style="flex-grow: 1; min-height: 0; display: flex; gap: 16px; padding: 10px 14px">'
             f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px">{t}<span style="font-size: 11.5px; color: {DIM}">Rows 37–44 of {sum(1 for r in recs if r["mol"] == 1)} matching the filter ({len(recs)} in frame). Row 41 is the picked atom (ringed in the viewport).</span></div>'
             f'<div style="width: 250px; flex-shrink: 0; display: flex; flex-direction: column; gap: 4px"><h3 style="{H3}">Global attributes</h3>{at}</div></div></div>')
    right = inspector("Data inspector", section("Picked atom", col(kv("Particle", str(pick["id"])), kv("Molecule", str(pick["mol"])), kv("Name", pick["gaff"]), kv("Charge", f"{pick['q']:+.4f} e"),
                                                                    kv("Position", " ".join(f"{v:.2f}" for v in pick["p"])), gap=6))
                      + section("Export table", col(select("As", "CSV · all columns"), toggle("Only filtered rows", True), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Density is mass ÷ cell volume of this built cell ({len(cell)} short chains, no solvent). It is low because the cell was packed loosely, not equilibrated.</div>')
    return frame(viewport(sc.svg(), f"Picked: particle {pick['id']}") + panel, right, "CAPS — data inspector", "Analyze › Visualize › Data inspector",
                 (f"<span>{len(recs)} particles · {nb} bonds · ρ = {rho:.3f} g/cm³</span>", "<span>every value here is a column you can colour, select or plot</span>"))


# ---------------------------------------------------------------- Render settings & overlays
def render_overlays():
    cell, recs = cell_recs()
    sh = mol_shape(recs)
    for r in recs:
        r["dcom"] = dist(r["p"], sh[r["mol"]]["com"])
    hi = math.ceil(max(r["dcom"] for r in recs))
    w, h = CW - 6, FULL
    sc = Scene("ro", w, h, yaw=0.55, pitch=0.4, persp=0.0, fog=0.25, atom_k=1.2, center=(L / 2, L / 2, L / 2))
    draw(sc, cell, recs, lambda r: ramp(VIRIDIS, r["dcom"] / hi))
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    sc._prep()
    fw, fh = int(w * 0.86), int(w * 0.86 * 9 / 16)
    fx, fy = (w - fw) // 2, (h - fh) // 2
    bar_px = 10 * sc.scale
    stops = "".join(f'<stop offset="{t / 8:.3f}" stop-color="{ramp(VIRIDIS, t / 8)}"></stop>' for t in range(9))
    mass = sum(MASS[r["e"]] for r in recs)
    rho = mass / NA / (L * 1e-8) ** 3
    ov = (f'<rect x="0" y="0" width="{w}" height="{fy}" fill="#000" fill-opacity="0.45"></rect><rect x="0" y="{fy + fh}" width="{w}" height="{h - fy - fh}" fill="#000" fill-opacity="0.45"></rect>'
          f'<rect x="0" y="{fy}" width="{fx}" height="{fh}" fill="#000" fill-opacity="0.45"></rect><rect x="{fx + fw}" y="{fy}" width="{w - fx - fw}" height="{fh}" fill="#000" fill-opacity="0.45"></rect>'
          f'<rect x="{fx}" y="{fy}" width="{fw}" height="{fh}" fill="none" stroke="{ACC}" stroke-dasharray="6 4"></rect>'
          f'<text x="{fx}" y="{fy - 8}" font-size="11" font-family="IBM Plex Mono" fill="{ACC}">render frame · 1920 × 1080</text>'
          f'<rect x="{fx + 10}" y="{fy + 12}" width="330" height="48" rx="6" fill="#0B0D0F" fill-opacity="0.8"></rect>'
          f'<text x="{fx + 18}" y="{fy + 30}" font-size="15" font-family="IBM Plex Sans" font-weight="600" fill="{TEXT}">PS melt · frame 120 · {len(recs)} atoms</text>'
          f'<text x="{fx + 18}" y="{fy + 50}" font-size="12" font-family="IBM Plex Mono" fill="{MUTED}">ρ = {rho:.3f} g/cm³ · L = {L:.0f} Å</text>'
          f'<defs><linearGradient id="roleg" x1="0" y1="1" x2="0" y2="0">{stops}</linearGradient></defs>'
          f'<rect x="{fx + fw - 46}" y="{fy + 24}" width="14" height="160" fill="url(#roleg)" stroke="{LINE}"></rect>'
          f'<text x="{fx + fw - 52}" y="{fy + 32}" text-anchor="end" font-size="11" font-family="IBM Plex Mono" fill="{TEXT}">{hi} Å</text>'
          f'<text x="{fx + fw - 52}" y="{fy + 184}" text-anchor="end" font-size="11" font-family="IBM Plex Mono" fill="{TEXT}">0 Å</text>'
          f'<text x="{fx + fw - 32}" y="{fy + 202}" text-anchor="end" font-size="11" font-family="IBM Plex Sans" fill="{MUTED}">DistanceToCOM</text>'
          f'<line x1="{fx + 20}" y1="{fy + fh - 24}" x2="{fx + 20 + bar_px:.1f}" y2="{fy + fh - 24}" stroke="{TEXT}" stroke-width="3"></line>'
          f'<text x="{fx + 20 + bar_px / 2:.1f}" y="{fy + fh - 32}" text-anchor="middle" font-size="11" font-family="IBM Plex Mono" fill="{TEXT}">10 Å</text>'
          + gizmo(fx + fw - 40, fy + fh - 40, 0.55, 0.4, 22))
    vp = f'<div style="position: relative; flex-grow: 1; min-height: 0; background: {BG0}; overflow: hidden">{sc.svg(overlay=ov)}</div>'
    layers = [("Text label", "[SourceFrame], [Particles], [Density]", True), ("Colour legend", "DistanceToCOM · viridis", True), ("Scale bar", "10 Å · orthographic only", True), ("Axis tripod", "bottom right", True), ("Python overlay", "draw your own", False)]
    ll = "".join(f'<div style="padding: 5px 0; border-bottom: 1px solid {BG2}">' + toggle(f'<span style="display: flex; flex-direction: column; text-align: left"><span style="font-size: 12px">{a}</span><span style="font-size: 10.5px; color: {DIM}; font-family: {MONO}">{b}</span></span>', on) + '</div>' for a, b, on in layers)
    right = inspector("Render", section("Output", col(row(field("Width", "1920", "px"), field("Height", "1080", "px"), gap=8), select("Renderer", "OpenGL · fast"), f'<div style="display: flex; flex-direction: column; gap: 4px"><span style="font-size: 11.5px; color: {MUTED}">Background</span>{seg(["Dark", "White", "Transparent"], "Dark", full=True)}</div>', gap=8))
                      + section("Quality", col(toggle("Ambient occlusion", True), toggle("Depth cue (fog)", True), row(select("Antialias", "4×"), select("Frames", "this frame"), gap=8), gap=8))
                      + section("Overlays", f'<div>{ll}</div>', pad=12)
                      + f'<div style="padding: 0 18px 12px; display: flex; gap: 8px">{btn("Render image", True, "eye")}{btn("Movie…")}</div>')
    return frame(vp, right, "CAPS — render and overlays", "Analyze › Visualize › Render",
                 ("<span>orthographic camera · scale bar is true 10 Å</span>", "<span>labels fill from live attributes</span>"))


BOARDS = (("FourViewports", four_viewports), ("WrapUnwrap", wrap_unwrap), ("MoleculeShape", molecule_shape),
          ("HistScatter", hist_scatter), ("DataInspector", data_inspector), ("RenderOverlays", render_overlays))

if __name__ == "__main__":
    import os
    os.makedirs("stage22/project", exist_ok=True)
    for name, fn in BOARDS:
        h = fn(); open(f"stage22/project/{name}.dc.html", "w").write(h); print(name, len(h))
