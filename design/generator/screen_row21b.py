import math
from lib import *
from screen_row20 import draw, CHAIN, VIRIDIS, ramp
from screen_row21 import cell_recs, mi, inspector, viewport, frame, CW, VH, BOTTOM, L


def hist_svg(vals, lo, hi, bins, w, h, colour_fn, label):
    hb = [0] * bins
    for v in vals:
        hb[min(bins - 1, max(0, int((v - lo) / (hi - lo) * bins)))] += 1
    mx = max(hb)
    bars = "".join(f'<rect x="{36 + i * (w - 46) / bins:.1f}" y="{h - 24 - c / mx * (h - 40):.1f}" width="{(w - 46) / bins - 2:.1f}" height="{c / mx * (h - 40):.1f}" rx="1.5" fill="{colour_fn((i + 0.5) / bins)}"></rect>' for i, c in enumerate(hb))
    ticks = "".join(f'<text x="{36 + t * (w - 46):.1f}" y="{h - 8}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{lo + t * (hi - lo):.1f}</text>' for t in (0, 0.25, 0.5, 0.75, 1))
    return (f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-label="{label}" style="display: block">{bars}{ticks}'
            f'<line x1="36" y1="{h - 24}" x2="{w - 10}" y2="{h - 24}" stroke="{LINE}"></line><text x="{w - 10}" y="12" text-anchor="end" font-size="10" fill="{DIM}">{label}</text></svg>')


def compute_property():
    cell, recs = cell_recs()
    com = {}
    for r in recs:
        c = com.setdefault(r["mol"], [0.0, 0.0, 0.0, 0.0])
        m = 12.011 if r["e"] == "C" else 1.008
        for k in range(3): c[k] += m * r["p"][k]
        c[3] += m
    com = {k: tuple(v[i] / v[3] for i in range(3)) for k, v in com.items()}
    for r in recs:
        r["dcom"] = dist(r["p"], com[r["mol"]])
    vals = [r["dcom"] for r in recs]
    lo, hi = 0.0, math.ceil(max(vals))
    sc = __import__("lib").Scene("cp", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.25)
    draw(sc, cell, recs, lambda r: ramp(VIRIDIS, r["dcom"] / hi))
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    hs = hist_svg(vals, lo, hi, 24, 560, 190, lambda t: ramp(VIRIDIS, t), "distance to own molecule's centre of mass (Å)")
    sample = table(["id", "mol", "name", "DistanceToCOM"], [[str(r["id"]), str(r["mol"]), r["gaff"], f"{r['dcom']:.3f}"] for r in recs[:5]], ["18%", "18%", "22%", "42%"], mono_cols=(0, 1, 2, 3), align_right=(0, 1, 3), fs=11.5, rowh=22)
    mean = sum(vals) / len(vals)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-shrink: 0">{hs}</div>'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 8px">{sample}'
              f'<span style="font-family: {MONO}; font-size: 12px">min {min(vals):.2f} · mean {mean:.2f} · max {max(vals):.2f} Å</span></div></div>')
    edit = (f'<div style="display: flex; flex-direction: column; gap: 4px"><label style="font-size: 11.5px; color: {MUTED}">Expression</label>'
            f'<div style="padding: 8px 10px; background: {BG0}; border: 1px solid {ACC}; border-radius: 6px; font-family: {MONO}; font-size: 12px; line-height: 1.6">'
            f'<span style="color: {SEL}">norm</span>(<span style="color: {SEL}">Position</span> − <span style="color: {SEL}">MoleculeCOM</span>(<span style="color: {SEL}">MoleculeIdentifier</span>))</div></div>')
    right = inspector("Compute property", section("Output", col(field("Property name", "DistanceToCOM", mono=True), select("Kind", "float · Å"), gap=8))
                      + section("Definition", col(edit, toggle("Only selected", False), toggle("Neighbour terms (sum over r < rc)", False), gap=8))
                      + section("Then", col(select("Colour by", "DistanceToCOM · viridis"), kv("Range", f"0 → {hi:.0f} Å"), gap=8))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Mass-weighted centres of each molecule (H included); values computed for all {len(recs)} atoms of this frame.</div>')
    return frame(viewport(sc.svg(), "Colour: DistanceToCOM") + bottom, right, "CAPS — compute property", "Analyze › Visualize › Compute property",
                 (f"<span>new property on {len(recs)} atoms</span>", f"<span>mean {mean:.2f} Å</span>"))


def create_bonds():
    cell, recs = cell_recs()
    cut = {("C", "C"): 1.70, ("C", "H"): 1.25, ("H", "H"): 0.0}
    topo = set()
    base = 0
    for atoms, bonds, _ in cell:
        for a, b in bonds:
            topo.add((base + a, base + b) if a < b else (base + b, base + a))
        base += len(atoms)
    P = [r["p"] for r in recs]
    E = [r["e"] for r in recs]
    cs = 1.8
    nc = int(L // cs)
    cs = L / nc
    grid = {}
    for i, p in enumerate(P):
        grid.setdefault(tuple(int((p[k] % L) // cs) % nc for k in range(3)), []).append(i)
    found = set()
    for i, p in enumerate(P):
        kx, ky, kz = (int((p[k] % L) // cs) % nc for k in range(3))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get(((kx + dx) % nc, (ky + dy) % nc, (kz + dz) % nc), ()):
                        if j <= i: continue
                        c = cut[tuple(sorted((E[i], E[j])))]
                        if c and mi(p, P[j]) < c:
                            found.add((i, j))
    both = topo & found
    extra = found - topo
    missing = topo - found
    lens = sorted(mi(P[a], P[b]) for a, b in topo)
    sc = __import__("lib").Scene("cb", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.1)
    draw(sc, cell, recs, lambda r: CHAIN[r["chain"] % 10] if r["e"] != "H" else mix(CHAIN[r["chain"] % 10], "#FFFFFF", 0.55))
    for a, b in list(extra)[:40]:
        sc.add_line(P[a], P[b], ERR, 3, None, 1)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    rows = [["In both", str(len(both)), "agree"], ["Cutoff only", str(len(extra)), "red · spurious"], ["Topology only", str(len(missing)), "missed by cutoff"]]
    ct = table(["Bonds", "Count", ""], rows, ["40%", "22%", "38%"], mono_cols=(1,), align_right=(1,), fs=12, rowh=26)
    pt = table(["Pair", "Cutoff (Å)"], [["C–C", "1.70"], ["C–H", "1.25"], ["H–H", "none"]], ["50%", "50%"], mono_cols=(0, 1), align_right=(1,), fs=12, rowh=26)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="width: 420px; flex-shrink: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Cutoff bonds vs file topology</h3>{ct}</div>'
              f'<div style="display: flex; flex-direction: column; gap: 8px; padding-top: 18px">{kv("Topology bonds", str(len(topo)))}{kv("Bond length range", f"{lens[0]:.2f} – {lens[-1]:.2f} Å")}'
              f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5; max-width: 380px">When a file has bonds, CAPS keeps them and uses cutoff bonds only as a check. Disagreements are listed, never silently merged.</span></div></div>')
    right = inspector("Create bonds", section("Mode", col(seg(["From topology", "Cutoff by pair"], "Cutoff by pair", full=True), pt, gap=8))
                      + section("Options", col(toggle("Periodic (minimum image)", True), toggle("Only between different molecules", False), toggle("Keep file bonds", True), gap=8)))
    return frame(viewport(sc.svg(), f"Cutoff bonds · {len(extra)} spurious shown in red") + bottom, right, "CAPS — create bonds", "Analyze › Visualize › Create bonds",
                 (f"<span>{len(found)} cutoff bonds · {len(topo)} in topology</span>", f"<span>{len(extra)} extra · {len(missing)} missing</span>"))


def replicate():
    cell, recs = cell_recs()
    n = len(recs)
    reps = (2, 2, 2)
    sc = __import__("lib").Scene("rp", CW - 6, VH - 8, yaw=0.55, pitch=0.35, persp=0.25, fog=0.35)
    for ix in range(reps[0]):
        for iy in range(reps[1]):
            for iz in range(reps[2]):
                o = (ix * L, iy * L, iz * L)
                orig = (ix, iy, iz) == (0, 0, 0)
                for m, (atoms, bonds, bb) in enumerate(cell):
                    pts = [add(atoms[j]["p"], o) for j in bb]
                    c = CHAIN[m % 10]
                    sc.add_tube(pts, c if orig else mix(c, BG0, 0.55), 0.5)
                sc.add_box(o, (L, 0, 0), (0, L, 0), (0, 0, L), ACC if orig else MUTED, 1.6 if orig else 0.8, None if orig else "4 4", 0.9 if orig else 0.5)
    total = n * reps[0] * reps[1] * reps[2]
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 28px; padding: 14px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="width: 360px; display: flex; flex-direction: column; gap: 8px">{kv("Images", " × ".join(map(str, reps)))}{kv("Atoms", f"{n:,} → {total:,}".replace(",", " "))}{kv("Molecules", f"{len(cell)} → {len(cell) * 8}")}{kv("Box", f"{L:.0f} Å → {L * 2:.0f} Å cubic")}</div>'
              f'<div style="flex-grow: 1; font-size: 12px; color: {MUTED}; line-height: 1.6">Replicas are views in the pipeline until you choose <b>Make real</b>, which writes a larger document with unique atom and molecule IDs. Original cell in amber; images dimmed. Display: backbone, because 8 × {n:,} atoms is the size where the lens takes over.</div></div>'.replace(",", " "))
    right = inspector("Replicate", section("Images", col(row(field("a", "2"), field("b", "2"), field("c", "2"), gap=6), toggle("Adjust simulation box", True), toggle("Unique IDs", True), gap=8))
                      + section("Operate on", col(toggle("Particles", True), toggle("Bonds", True), toggle("Surfaces and vectors", True), gap=8))
                      + f'<div style="padding: 0 18px 12px">{btn("Make real", ic="copy")}</div>')
    return frame(viewport(sc.svg(), f"Backbone · {total:,} atoms across 8 images".replace(",", " "), "Backbone") + bottom, right, "CAPS — replicate", "Analyze › Visualize › Replicate",
                 (f"<span>{total:,} atoms in view (8 images)</span>".replace(",", " "), "<span>images are virtual until made real</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage20/project", exist_ok=True)
    for name, fn in (("ComputeProperty", compute_property), ("CreateBonds", create_bonds), ("Replicate", replicate)):
        h = fn(); open(f"stage20/project/{name}.dc.html", "w").write(h); print(name, len(h))
