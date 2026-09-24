"""Row 24 — trajectories: displacements, MSD, time series, particle paths, smoothing, timeline.

All motion comes from a labelled TOY trajectory (rigid-chain Brownian motion + atom vibration) generated here,
so every curve is computed, but none of it is molecular dynamics."""
import math, random
from lib import *
from screen_studio2 import bar, spacer
from screen_row20 import draw, _shell, CHAIN, VIRIDIS, ramp
from screen_row21 import cell_recs, mi, inspector, viewport, CW, VH, BOTTOM, L, MASS
from screen_row22 import chain_col, H3, FULL
from screen_row23 import code

NF = 51          # frames 0..50
D_IN = 0.05      # Å² per τ, COM translational diffusion of the toy
SIG_V = 0.15     # Å, per-frame vibration
SIG_R = 0.02     # rad, per-frame rotation step
TOY = f"toy trajectory: rigid-chain Brownian motion, D = {D_IN} Å²/τ, vibration σ = {SIG_V} Å — not MD"

_T = None


def rotm(axis, th):
    x, y, z = axis
    c, s, C = math.cos(th), math.sin(th), 1 - math.cos(th)
    return [[c + x * x * C, x * y * C - z * s, x * z * C + y * s], [y * x * C + z * s, c + y * y * C, y * z * C - x * s], [z * x * C - y * s, z * y * C + x * s, c + z * z * C]]


def mm(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def mv(A, v):
    return tuple(sum(A[i][k] * v[k] for k in range(3)) for i in range(3))


def traj():
    """frames[t][i] = position of atom i (unwrapped); coms[t][m] = centre of chain m."""
    global _T
    if _T: return _T
    cell, recs = cell_recs()
    rnd = random.Random(2024)
    mols = sorted({r["mol"] for r in recs})
    idx = {m: [i for i, r in enumerate(recs) if r["mol"] == m] for m in mols}
    com0 = {}
    for m in mols:
        M = sum(MASS[recs[i]["e"]] for i in idx[m])
        com0[m] = tuple(sum(MASS[recs[i]["e"]] * recs[i]["p"][k] for i in idx[m]) / M for k in range(3))
    R = {m: [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)] for m in mols}
    com = dict(com0)
    step = math.sqrt(2 * D_IN)
    frames, coms = [], []
    for t in range(NF):
        if t:
            for m in mols:
                com[m] = tuple(com[m][k] + rnd.gauss(0, step) for k in range(3))
                ax = norm((rnd.gauss(0, 1), rnd.gauss(0, 1), rnd.gauss(0, 1)))
                R[m] = mm(rotm(ax, rnd.gauss(0, SIG_R)), R[m])
        pos = [None] * len(recs)
        for m in mols:
            for i in idx[m]:
                base = add(com[m], mv(R[m], sub(recs[i]["p"], com0[m])))
                pos[i] = tuple(base[k] + (rnd.gauss(0, SIG_V) if t else 0.0) for k in range(3))
        frames.append(pos)
        coms.append(dict(com))
    _T = (cell, recs, frames, coms, mols, idx)
    return _T


def cell_at(cell, pos):
    out, base = [], 0
    for atoms, bonds, bb in cell:
        out.append(([dict(a, p=pos[base + j]) for j, a in enumerate(atoms)], bonds, bb))
        base += len(atoms)
    return out


def tframe(content_center, right, title, crumb, status, fr=50):
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
              + bar(tbtn("rotate", "Orbit", True), tbtn("move", "Pan"), tbtn("search", "Zoom"), sep(), tbtn("cursor", "Pick"), tbtn("ruler", "Measure"), spacer(),
                    chip("toy trajectory", WARN, BG2), chip(f"frame {fr} / {NF - 1}", MUTED, BG2, True))
              + content_center + '</div>')
    return _shell(title, f'{icon("chevr", 12, DIM)}<span>{crumb}</span>', f'<div style="flex-grow: 1; display: flex; min-height: 0">{center}{right}</div>', status)


def toy_note():
    return f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Motion here is a {TOY}. It has no forces, so chains can drift into each other. Swap in a real dump and every number on this screen is recomputed.</div>'


def contacts(pos, recs, cut):
    heavy = [i for i, r in enumerate(recs) if r["e"] != "H"]
    nc = int(L // cut); cs = L / nc
    grid = {}
    for i in heavy:
        grid.setdefault(tuple(int((pos[i][k] % L) // cs) % nc for k in range(3)), []).append(i)
    pairs = set()
    for i in heavy:
        kx, ky, kz = (int((pos[i][k] % L) // cs) % nc for k in range(3))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get(((kx + dx) % nc, (ky + dy) % nc, (kz + dz) % nc), ()):
                        if j > i and recs[i]["mol"] != recs[j]["mol"] and mi(pos[i], pos[j]) < cut:
                            pairs.add((i, j))
    return pairs


def n_clusters(pairs, recs, mols):
    parent = {m: m for m in mols}
    def f(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, j in pairs:
        a, b = f(recs[i]["mol"]), f(recs[j]["mol"])
        if a != b: parent[a] = b
    return len({f(m) for m in mols})


def bars_hist(vals, lo, hi, nb, w, h, colour, label, xt):
    c = [0] * nb
    for v in vals:
        if lo <= v < hi: c[int((v - lo) / (hi - lo) * nb)] += 1
    mx = max(c) or 1
    bars = "".join(f'<rect x="{34 + i * (w - 44) / nb:.1f}" y="{h - 22 - v / mx * (h - 38):.1f}" width="{(w - 44) / nb - 1.5:.1f}" height="{v / mx * (h - 38):.1f}" fill="{colour(lo + (i + 0.5) * (hi - lo) / nb)}"></rect>' for i, v in enumerate(c) if v)
    ticks = "".join(f'<text x="{34 + (v - lo) / (hi - lo) * (w - 44):.1f}" y="{h - 7}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{v:g}</text>' for v in xt)
    return (f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-label="{label}" style="display: block">{bars}{ticks}'
            f'<line x1="34" y1="{h - 22}" x2="{w - 10}" y2="{h - 22}" stroke="{LINE}"></line><text x="34" y="11" font-size="10" fill="{DIM}">{label}</text></svg>')


# ---------------------------------------------------------------- Displacements
def displacements():
    cell, recs, frames, coms, mols, idx = traj()
    p0, p1 = frames[0], frames[-1]
    heavy = [i for i, r in enumerate(recs) if r["e"] != "H"]
    d = {i: dist(p0[i], p1[i]) for i in heavy}
    hi = math.ceil(max(d.values()))
    sc = Scene("dp", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.0)
    for m, (atoms, bonds, bb) in enumerate(cell):
        base = idx[m + 1][0]
        sc.add_tube([p1[base + j] for j in bb], mix(CHAIN[m % 10], BG0, 0.7), 0.3)
    for i in heavy[::2]:
        c = ramp(VIRIDIS, d[i] / hi)
        sc.add_line(p0[i], p1[i], c, 1.6, None, 0.95)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    hs = bars_hist(list(d.values()), 0, hi, 28, 470, BOTTOM - 20, lambda v: ramp(VIRIDIS, v / hi), "|displacement| frame 0 → 50 (Å), heavy atoms", list(range(0, hi + 1)))
    per = []
    for m in mols:
        ds = [d[i] for i in idx[m] if i in d]
        per.append([str(m), f"{sum(ds) / len(ds):.2f}", f"{max(ds):.2f}", f"{dist(coms[0][m], coms[-1][m]):.2f}"])
    t = table(["mol", "mean |d| (Å)", "max |d| (Å)", "COM shift (Å)"], per[:6], ["16%", "28%", "28%", "28%"], mono_cols=(0, 1, 2, 3), align_right=(0, 1, 2, 3), fs=11.5, rowh=22)
    mean = sum(d.values()) / len(d)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-shrink: 0">{hs}</div><div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Per molecule · 6 of {len(mols)}</h3>{t}</div></div>')
    right = inspector("Displacements", section("Reference", col(seg(["Frame 0", "Previous", "Relative"], "Frame 0", full=True), toggle("Use unwrapped positions", True), toggle("Subtract system drift", False), gap=8))
                      + section("Output", col(kv("Property", "Displacement, |Displacement|"), select("Show as", "arrows · every 2nd heavy atom"), select("Colour", "|d| · viridis"), gap=8))
                      + toy_note())
    return tframe(viewport(sc.svg(), f"Arrows: frame 0 → 50 · 0 → {hi} Å", "Backbone") + bottom, right, "CAPS — displacements", "Analyze › Visualize › Displacements",
                  (f"<span>{len(heavy)} heavy atoms · mean |d| = {mean:.2f} Å</span>", f"<span>{TOY}</span>"))


# ---------------------------------------------------------------- MSD
def msd_board():
    cell, recs, frames, coms, mols, idx = traj()
    heavy = [i for i, r in enumerate(recs) if r["e"] != "H"]
    lags = list(range(1, 26))
    msd_a, msd_c = [], []
    for lag in lags:
        sa = sc_ = 0.0; na = nc = 0
        for t0 in range(0, NF - lag):
            for i in heavy[::3]:
                sa += sum((frames[t0 + lag][i][k] - frames[t0][i][k]) ** 2 for k in range(3)); na += 1
            for m in mols:
                sc_ += sum((coms[t0 + lag][m][k] - coms[t0][m][k]) ** 2 for k in range(3)); nc += 1
        msd_a.append(sa / na); msd_c.append(sc_ / nc)
    fit = [(l, v) for l, v in zip(lags, msd_c) if 5 <= l <= 20]
    sx = sum(l for l, _ in fit); sy = sum(v for _, v in fit); sxx = sum(l * l for l, _ in fit); sxy = sum(l * v for l, v in fit); n = len(fit)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx); icpt = (sy - slope * sx) / n
    D_fit = slope / 6
    ym = math.ceil(max(msd_a) + 1)
    pl = plot(700, 420, [([(0, 0)] + list(zip(lags, msd_a)), SEL, 2, None), ([(0, 0)] + list(zip(lags, msd_c)), ACC, 2, None),
                         ([(0, 0), (25, 6 * D_IN * 25)], DIM, 1, "4 4"), ([(5, icpt + slope * 5), (20, icpt + slope * 20)], TEXT, 1.2, "2 3")],
              (0, 25), (0, ym), [0, 5, 10, 15, 20, 25], [0, ym / 4, ym / 2, 3 * ym / 4, ym], "lag (τ)", "MSD (Å²)", pad=(48, 16, 16, 30),
              markers=[(l, v, SEL) for l, v in zip(lags, msd_a)] + [(l, v, ACC) for l, v in zip(lags, msd_c)])
    leg = "".join(f'<div style="display: flex; align-items: center; gap: 8px; font-size: 12px"><span style="width: 18px; height: 3px; background: {c}; {"border-top: 1px dashed " + c + "; height: 0" if dsh else ""}"></span>{t}</div>'
                  for t, c, dsh in (("heavy atoms (every 3rd)", SEL, False), ("chain centres of mass", ACC, False), ("6·D·t for the toy's input D", DIM, True), ("fit, lag 5–20", TEXT, True)))
    tb = table(["", "D (Å²/τ)"], [["toy input", f"{D_IN:.4f}"], ["fitted from COM MSD", f"{D_fit:.4f}"], ["ratio", f"{D_fit / D_IN:.2f}"]], ["60%", "40%"], mono_cols=(1,), align_right=(1,), fs=12, rowh=28)
    content = (f'<div style="flex-grow: 1; display: flex; gap: 16px; padding: 16px; min-height: 0; overflow: hidden">'
               f'<div style="flex-shrink: 0; display: flex; flex-direction: column; gap: 10px"><h3 style="{H3}">Mean-squared displacement · time-origin averaged</h3>{pl}'
               f'<div style="display: flex; gap: 18px; flex-wrap: wrap">{leg}</div></div>'
               f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px"><h3 style="{H3}">Diffusion coefficient</h3>{tb}'
               f'<span style="font-size: 12px; color: {MUTED}; line-height: 1.6">Atom MSD starts high because of vibration (≈ 6σ² = {6 * SIG_V ** 2:.3f} Å² at any lag) and rotation. Chain centres only diffuse, so their slope gives D. '
               f'Ten chains and {NF} frames make the fit noisy: this is the check CAPS runs to show you the error bar you actually have.</span></div></div>')
    right = inspector("MSD", section("Groups", col(select("Atoms", "heavy · every 3rd"), toggle("Chain centres of mass", True), toggle("Per molecule curves", False), gap=8))
                      + section("Averaging", col(select("Time origins", "every frame"), row(field("Max lag", "25", "τ"), field("Fit", "5–20", "τ"), gap=8), toggle("Remove system drift", True), gap=8))
                      + toy_note())
    return tframe(content, right, "CAPS — MSD", "Analyze › Visualize › MSD",
                  (f"<span>D fit = {D_fit:.4f} Å²/τ vs input {D_IN}</span>", f"<span>{TOY}</span>"))


# ---------------------------------------------------------------- Time series & time averaging
_TS = None


def series():
    global _TS
    if _TS: return _TS
    cell, recs, frames, coms, mols, idx = traj()
    ncon, ncl = [], []
    for pos in frames:
        pairs = contacts(pos, recs, 4.5)
        ncon.append(len(pairs))
        ncl.append(n_clusters([p for p in pairs if mi(pos[p[0]], pos[p[1]]) < 3.3], recs, mols))
    _TS = (ncon, ncl)
    return _TS


def running(v, w):
    out = []
    for i in range(len(v)):
        a = max(0, i - w // 2); b = min(len(v), i + w // 2 + 1)
        out.append(sum(v[a:b]) / (b - a))
    return out


def time_series():
    cell, recs, frames, coms, mols, idx = traj()
    ncon, ncl = series()
    avg = running(ncon, 9)
    lo, hi = min(ncon), max(ncon)
    lo = (lo // 50) * 50; hi = (hi // 50 + 1) * 50
    pl = plot(CW - 40, 250, [(list(enumerate(ncon)), SEL, 1.4, None), (list(enumerate(avg)), ACC, 2.2, None)], (0, NF - 1), (lo, hi),
              [0, 10, 20, 30, 40, 50], [lo, (lo + hi) / 2, hi], "frame", "inter-chain heavy-atom contacts < 4.5 Å", pad=(48, 16, 16, 30))
    ymax = max(ncl) + 1
    pts = []
    for t, v in enumerate(ncl):
        pts += [(t, v), (t + 1, v)]
    pl2 = plot(CW - 40, 170, [(pts, "#C79BE8", 1.8, None)], (0, NF), (0, ymax), [0, 10, 20, 30, 40, 50], list(range(0, ymax + 1, max(1, ymax // 4))), "frame", "clusters at 3.3 Å", pad=(48, 16, 16, 30))
    ch = [t for t in range(1, NF) if ncl[t] != ncl[t - 1]]
    content = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 12px; padding: 14px 16px; min-height: 0; overflow: hidden">'
               f'<div style="display: flex; gap: 20px; font-size: 12px; color: {MUTED}"><span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 16px; height: 2px; background: {SEL}"></span>per frame</span>'
               f'<span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 16px; height: 3px; background: {ACC}"></span>running mean · 9 frames</span></div>{pl}{pl2}'
               f'<div style="display: flex; gap: 28px; font-family: {MONO}; font-size: 12px"><span>contacts {min(ncon)}–{max(ncon)} · mean {sum(ncon) / NF:.0f}</span><span>cluster count changes {len(ch)} times</span></div>'
               f'<span style="font-size: 12px; color: {MUTED}; line-height: 1.6; max-width: 900px">Contacts rise because the toy has no forces: chains drift into each other. A real trajectory would hold contacts near a plateau; this screen only shows how the series, the running mean and the markers are drawn.</span></div>')
    right = inspector("Time series", section("Quantities", col(select("Series 1", "Contacts (inter-chain, 4.5 Å)"), select("Series 2", "Cluster count (3.3 Å)"), toggle("Add series from any attribute", True), gap=8))
                      + section("Time averaging", col(row(field("Window", "9", "frames"), select("Kind", "centred mean"), gap=8), toggle("Also average per-atom properties", False), gap=8))
                      + section("Frames", col(kv("Evaluated", f"{NF} of {NF}"), kv("Cost", "[measure]"), gap=6))
                      + toy_note())
    return tframe(content, right, "CAPS — time series", "Analyze › Visualize › Time series",
                  (f"<span>{NF} frames evaluated through the pipeline</span>", f"<span>{TOY}</span>"))


# ---------------------------------------------------------------- Particle paths
def particle_paths():
    cell, recs, frames, coms, mols, idx = traj()
    sc = Scene("pp", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.25, atom_k=1.0, center=(L / 2, L / 2, L / 2))
    for m, (atoms, bonds, bb) in enumerate(cell):
        base = idx[m + 1][0]
        sc.add_tube([frames[-1][base + j] for j in bb], mix(CHAIN[m % 10], BG0, 0.6), 0.35)
    tot = []
    for m in mols:
        c = CHAIN[(m - 1) % 10]
        path = [coms[t][m] for t in range(NF)]
        for t in range(NF - 1):
            sc.add_line(path[t], path[t + 1], c, 2.4, None, 0.35 + 0.65 * t / (NF - 1))
        tot.append(sum(dist(path[t], path[t + 1]) for t in range(NF - 1)))
    for m in mols:
        base = idx[m][0]
        atoms = cell[m - 1][0]
        tip = atoms.index(next(a for a in atoms if a["e"] == "C" and a.get("ar")))
        for t in range(NF - 1):
            sc.add_line(frames[t][base + tip], frames[t + 1][base + tip], "#E9ECEF", 1.0, None, 0.25)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    rows = [[f'<span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 10px; height: 10px; border-radius: 2px; background: {CHAIN[(m - 1) % 10]}"></span>{m}</span>',
             f"{tot[k]:.2f}", f"{dist(coms[0][m], coms[-1][m]):.2f}", f"{dist(coms[0][m], coms[-1][m]) / tot[k]:.2f}"] for k, m in enumerate(mols)][:6]
    t = table(["mol", "path length (Å)", "net shift (Å)", "net / path"], rows, ["16%", "30%", "28%", "26%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=11.5, rowh=22)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Centre-of-mass paths · frames 0–50</h3>{t}</div>'
              f'<div style="width: 280px; flex-shrink: 0; font-size: 11.5px; color: {DIM}; line-height: 1.6; padding-top: 18px">Thick lines: chain centres, fading from old to new. Thin white lines: one ring carbon per chain. Net / path near 0 means the walk turned back on itself — expected for diffusion.</div></div>')
    right = inspector("Trajectory lines", section("Trace", col(select("Particles", "chain centres + 1 ring C per chain"), row(field("From", "0"), field("To", "50"), gap=8), toggle("Unwrapped", True), gap=8))
                      + section("Style", col(row(field("Width", "0.25", "Å"), select("Colour", "by molecule"), gap=8), toggle("Fade with age", True), toggle("Show only up to current frame", True), gap=8))
                      + toy_note())
    return tframe(viewport(sc.svg(), "Paths: centres of mass, frames 0 → 50", "Backbone") + bottom, right, "CAPS — trajectory lines", "Analyze › Visualize › Trajectory lines",
                  (f"<span>{len(mols)} centre paths · {NF} frames</span>", f"<span>{TOY}</span>"))


# ---------------------------------------------------------------- Smooth trajectory (time averaging of positions)
def smooth():
    cell, recs, frames, coms, mols, idx = traj()
    t0, w = 25, 9
    win = frames[t0 - w // 2: t0 + w // 2 + 1]
    avg = [tuple(sum(f[i][k] for f in win) / len(win) for k in range(3)) for i in range(len(recs))]
    inst = frames[t0]
    cc = []
    for m, (atoms, bonds, bb) in enumerate(cell):
        base = idx[m + 1][0]
        for a, b in bonds:
            if atoms[a]["e"] == atoms[b]["e"] == "C":
                cc.append((base + a, base + b))
    li = [dist(inst[a], inst[b]) for a, b in cc]
    la = [dist(avg[a], avg[b]) for a, b in cc]
    l0 = [dist(frames[0][a], frames[0][b]) for a, b in cc]
    rms = math.sqrt(sum(dist(inst[i], avg[i]) ** 2 for i in range(len(recs))) / len(recs))
    hw = (CW - 6 - 6) // 2
    tiles = []
    for tag, pos in ((f"Frame {t0} as stored", inst), (f"Frames {t0 - w // 2}–{t0 + w // 2} averaged", avg)):
        sc = Scene("sm" + tag[0:2].strip(), hw, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.1, center=(L / 2, L / 2, L / 2), scale=VH / (L * 2.1))
        draw(sc, cell_at(cell, pos), recs, chain_col)
        sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
        tiles.append(f'<div style="position: relative; width: {hw}px; height: {VH}px; background: {BG0}; overflow: hidden">{sc.svg(overlay=label_pill(10, 18, tag, TEXT, fs=11))}</div>')
    split = f'<div style="flex-grow: 1; min-height: 0; display: flex; gap: 6px; background: {LINE}">{"".join(tiles)}</div>'
    stat = lambda v: f"{sum(v) / len(v):.3f} ± {math.sqrt(sum((x - sum(v) / len(v)) ** 2 for x in v) / len(v)):.3f}"
    t = table(["C–C bond length", "mean ± s.d. (Å)"], [["built (frame 0)", stat(l0)], [f"frame {t0}, stored", stat(li)], [f"{w}-frame average", stat(la)]], ["55%", "45%"], mono_cols=(1,), align_right=(1,), fs=12, rowh=28, hl={2})
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="width: 440px; flex-shrink: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">What averaging does to geometry</h3>{t}</div>'
              f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 8px; padding-top: 18px">{kv("RMS shift stored → averaged", f"{rms:.3f} Å")}'
              f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.6">Averaging removes vibration, but chains that rotate during the window average to slightly shrunken shapes. CAPS labels smoothed frames and never feeds them to bond-length or angle analysis unless you ask.</span></div></div>')
    right = inspector("Smooth trajectory", section("Window", col(row(field("Frames", str(w)), select("Kind", "centred mean"), gap=8), toggle("Unwrap before averaging", True), gap=8))
                      + section("Apply to", col(toggle("Positions", True), toggle("Per-atom properties", False), toggle("Mark frames as smoothed", True), gap=8))
                      + toy_note())
    return tframe(split + bottom, right, "CAPS — smooth trajectory", "Analyze › Visualize › Smooth trajectory",
                  (f"<span>C–C {stat(li)} → {stat(la)} Å</span>", f"<span>{TOY}</span>"), fr=t0)


# ---------------------------------------------------------------- Timeline & live pipeline
def timeline():
    cell, recs, frames, coms, mols, idx = traj()
    ncon, ncl = series()
    t0 = 32
    sc = Scene("tl", CW - 6, VH + 20, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.15, center=(L / 2, L / 2, L / 2))
    draw(sc, cell_at(cell, frames[t0]), recs, chain_col)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    tw = CW - 32 - 60
    X = lambda t: 50 + t / (NF - 1) * tw
    mx, mn = max(ncon), min(ncon)
    spark = " ".join(f"{X(t):.1f},{60 - (v - mn) / max(1, mx - mn) * 44:.1f}" for t, v in enumerate(ncon))
    cached = t0
    cache = "".join(f'<rect x="{X(t) - tw / (NF - 1) / 2:.1f}" y="70" width="{tw / (NF - 1) - 1:.1f}" height="6" fill="{OK if t <= cached else BG3}"></rect>' for t in range(NF))
    ch = [t for t in range(1, NF) if ncl[t] != ncl[t - 1]]
    marks = "".join(f'<g><line x1="{X(t):.1f}" y1="10" x2="{X(t):.1f}" y2="80" stroke="#C79BE8" stroke-dasharray="2 3"></line><circle cx="{X(t):.1f}" cy="10" r="4" fill="#C79BE8"></circle></g>' for t in ch)
    ticks = "".join(f'<text x="{X(t):.1f}" y="96" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{t}</text>' for t in range(0, NF, 5))
    head = f'<line x1="{X(t0):.1f}" y1="4" x2="{X(t0):.1f}" y2="82" stroke="{ACC}" stroke-width="2"></line><rect x="{X(t0) - 16:.1f}" y="-2" width="32" height="14" rx="3" fill="{ACC}"></rect><text x="{X(t0):.1f}" y="9" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{ACC_INK}">{t0}</text>'
    tsvg = (f'<svg width="{tw + 60}" height="100" viewBox="0 -4 {tw + 60} 104" role="img" aria-label="Timeline with contacts sparkline, cache and cluster-change markers" style="display: block">'
            f'<text x="0" y="40" font-size="10" fill="{DIM}">contacts</text><polyline points="{spark}" fill="none" stroke="{SEL}" stroke-width="1.4"></polyline>'
            f'<text x="0" y="77" font-size="10" fill="{DIM}">cache</text>{cache}{marks}{ticks}{head}</svg>')
    transport = row(btn("", ic="undo", small=True), btn("", ic="play", small=True), btn("", ic="redo", small=True), gap=6).replace("<button", '<button aria-label="step"', 3)
    bottom = (f'<div style="height: {BOTTOM - 20}px; flex-shrink: 0; display: flex; flex-direction: column; gap: 8px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="display: flex; align-items: center; gap: 12px">{transport}<span style="font-family: {MONO}; font-size: 12px">frame {t0} / {NF - 1} · t = {t0} τ</span>'
              f'<span style="flex-grow: 1"></span><span style="display: inline-flex; align-items: center; gap: 6px; font-size: 11.5px; color: {MUTED}"><span style="width: 8px; height: 8px; border-radius: 50%; background: #C79BE8"></span>cluster count changes</span>'
              f'<span style="display: inline-flex; align-items: center; gap: 6px; font-size: 11.5px; color: {MUTED}"><span style="width: 10px; height: 6px; background: {OK}"></span>pipeline cached</span></div>{tsvg}</div>')
    steps = [("Unwrap", True), ("Create bonds", True), ("Coordination & RDF", True), ("Cluster analysis", True), ("Colour by molecule", True), ("Python · backbone conformation", False)]
    sl = "".join(f'<div style="display: flex; align-items: center; gap: 8px; padding: 5px 0; border-bottom: 1px solid {BG2}; font-size: 12px">{icon("check", 14, OK) if ok else icon("pause", 14, DIM)}<span style="flex-grow: 1">{s}</span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{"cached" if ok else "paused while scrubbing"}</span></div>' for s, ok in steps)
    right = inspector("Timeline", section("Pipeline at this frame", f'<div>{sl}</div>')
                      + section("Playback", col(row(field("Every", "1", "frame"), field("FPS", "24"), gap=8), toggle("Loop", True), toggle("Pause slow steps while scrubbing", True), gap=8))
                      + section("Markers", col(kv("Cluster changes", f"{len(ch)} in {NF} frames"), toggle("Show on timeline", True), gap=6))
                      + toy_note())
    vp = f'<div style="position: relative; flex-grow: 1; min-height: 0; background: {BG0}; overflow: hidden">{sc.svg()}</div>'
    return tframe(vp + bottom, right, "CAPS — timeline", "Analyze › Visualize › Timeline",
                  (f"<span>frames 0–{cached} cached · {len(ch)} cluster-change markers</span>", f"<span>{TOY}</span>"), fr=t0)


BOARDS = (("Displacements", displacements), ("MsdBoard", msd_board), ("TimeSeries", time_series),
          ("ParticlePaths", particle_paths), ("SmoothTrajectory", smooth), ("Timeline", timeline))

if __name__ == "__main__":
    import os, sys
    os.makedirs("stage24/project", exist_ok=True)
    only = sys.argv[1:]
    for name, fn in BOARDS:
        if only and name not in only: continue
        h = fn(); open(f"stage24/project/{name}.dc.html", "w").write(h); print(name, len(h))
