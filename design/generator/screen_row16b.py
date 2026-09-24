import math, random
from lib import *
from mols import *
from screen_app import shell, pagehead, card
from screen_studio2 import studio_shell, bar, spacer
from screen_builders import footer, hud

VW = 1440 - 72 - 380
VH = 900 - 44 - 48 - 26

# Gasteiger & Marsili, Tetrahedron 36, 3219 (1980): chi = a + b q + c q^2
GM = {"H": (7.17, 6.24, -0.56), "C3": (7.98, 9.18, 1.88), "C2": (8.79, 9.32, 1.51)}


def gasteiger(atoms, bonds, iters=6):
    typ = ["H" if a["e"] == "H" else ("C2" if a.get("ar") else "C3") for a in atoms]
    nb = [[] for _ in atoms]
    for i, j in bonds:
        nb[i].append(j); nb[j].append(i)
    q = [0.0] * len(atoms)
    for k in range(iters):
        chi = [GM[t][0] + GM[t][1] * q[i] + GM[t][2] * q[i] ** 2 for i, t in enumerate(typ)]
        dq = [0.0] * len(atoms)
        for i in range(len(atoms)):
            for j in nb[i]:
                if chi[j] > chi[i]:
                    t = typ[i]; plus = 20.02 if t == "H" else sum(GM[t])
                else:
                    t = typ[j]; plus = 20.02 if t == "H" else sum(GM[t])
                dq[i] += (chi[j] - chi[i]) / plus * 0.5 ** (k + 1)
        q = [q[i] + dq[i] for i in range(len(q))]
    return q, typ


def ps_bonds(atoms, backbone, stereo):
    """Bonds replayed from the construction order of polystyrene(hydrogens=True)."""
    nb = len(backbone)
    bonds = [(k, k + 1) for k in range(nb - 1)]
    x = nb
    for k in range(nb):
        if k % 2 == 1:
            r = list(range(x, x + 6))
            bonds += [(r[m], r[(m + 1) % 6]) for m in range(6)] + [(k, r[0])]
            bonds += [(r[m], x + 6 + m - 1) for m in range(1, 6)] + [(k, x + 11)]
            x += 12
        else:
            bonds += [(k, x), (k, x + 1)]
            x += 2
    bonds += [(0, x), (nb - 1, x + 1)]
    assert x + 2 == len(atoms)
    return bonds


def charges():
    atoms, backbone, stereo = polystyrene(8, "atactic", seed=5)
    bonds = ps_bonds(atoms, backbone, stereo)
    q, typ = gasteiger(atoms, bonds)
    net = sum(q)
    qmax = max(abs(x) for x in q)
    POS, NEG = "#2271DB", "#E35049"

    def qc(x):
        t = min(1.0, abs(x) / qmax)
        return mix("#E9ECEF", POS if x > 0 else NEG, t)
    sc = Scene("qq", VW, VH - 210, yaw=-0.1, pitch=0.5, persp=0.2, fog=0.4)
    sc.add_atoms([dict(a, c=qc(q[i])) for i, a in enumerate(atoms)], bonds)
    view = sc.svg()
    groups = {}
    for i, a in enumerate(atoms):
        if a["e"] == "H":
            nbr = next(j for b in bonds for j in b if i in b and j != i)
            key = "H on ring" if atoms[nbr].get("ar") else "H on backbone"
        elif a.get("ar"):
            key = "C aromatic"
        else:
            key = "C backbone"
        groups.setdefault(key, []).append(q[i])
    gt = table(["Group", "n", "mean q (e)", "range"], [[k, str(len(v)), f"{sum(v) / len(v):+.4f}", f"{min(v):+.3f} … {max(v):+.3f}"] for k, v in sorted(groups.items())],
               ["30%", "10%", "24%", "36%"], mono_cols=(1, 2, 3), align_right=(1, 2), fs=12, rowh=28)
    # histogram
    bins = 24
    lo, hi = -qmax, qmax
    h = [0] * bins
    for x in q:
        h[min(bins - 1, int((x - lo) / (hi - lo) * bins))] += 1
    hm = max(h)
    hw, hh = 560, 120
    hs = "".join(f'<rect x="{30 + i * (hw - 40) / bins:.1f}" y="{hh - 22 - v / hm * (hh - 34):.1f}" width="{(hw - 40) / bins - 2:.1f}" height="{v / hm * (hh - 34):.1f}" rx="1.5" fill="{qc(lo + (i + 0.5) * (hi - lo) / bins)}"></rect>' for i, v in enumerate(h))
    hist = (f'<svg width="{hw}" height="{hh}" viewBox="0 0 {hw} {hh}" role="img" aria-label="Histogram of partial charges" style="display: block">{hs}'
            f'<line x1="30" y1="{hh - 22}" x2="{hw - 10}" y2="{hh - 22}" stroke="{LINE}"></line>'
            + "".join(f'<text x="{30 + (v - lo) / (hi - lo) * (hw - 40):.1f}" y="{hh - 6}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{v:+.2f}</text>' for v in (lo, lo / 2, 0, hi / 2, hi))
            + '</svg>')
    strip = (f'<div style="height: 200px; flex-shrink: 0; display: flex; gap: 22px; padding: 14px 20px; background: {BG1}; border-top: 1px solid {LINE}">'
             f'<div style="display: flex; flex-direction: column; gap: 8px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Distribution · {len(q)} atoms</h3>{hist}</div>'
             f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 10px; justify-content: center">{kv("Net charge", f"{net:+.6f} e", vcol=OK)}{kv("Target", "0 (neutral chain)")}{kv("Largest |q|", f"{qmax:.4f} e")}'
             f'<div style="height: 8px; border-radius: 4px; background: linear-gradient(90deg, {NEG}, #E9ECEF, {POS})"></div><div style="display: flex; justify-content: space-between; font-family: {MONO}; font-size: 10.5px; color: {DIM}"><span>−{qmax:.3f}</span><span>0</span><span>+{qmax:.3f}</span></div></div></div>')
    methods = [("Gasteiger–Marsili", "instant · shown here", True), ("AM1-BCC", "needs AmberTools · [measure] per molecule", False), ("RESP · HF/6-31G*", "external QM · import .chg", False), ("From force field", "OPLS-AA / CGenFF library charges", False)]
    ml = "".join(f'<div style="display: flex; align-items: center; gap: 10px; padding: 8px 10px; border-radius: 7px; background: {BG3 if on else "transparent"}; border: 1px solid {ACC if on else LINE}">'
                 f'<span style="width: 14px; height: 14px; border-radius: 50%; border: 2px solid {ACC if on else DIM}; flex-shrink: 0"></span><span style="display: flex; flex-direction: column; gap: 2px"><span style="font-size: 12.5px">{a}</span><span style="font-size: 11px; color: {DIM}">{b}</span></span></div>' for a, b, on in methods)
    right = (panel_head("Partial charges", chip("computed", OK, BG2))
             + section("Method", col(ml, gap=6))
             + section("By group", gt)
             + f'<div style="padding: 0 18px 10px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Charges above are Gasteiger–Marsili (6 iterations) computed for the drawn chain; net charge is zero by construction. Gasteiger &amp; Marsili, <i>Tetrahedron</i> 36, 3219 (1980).</div>'
             + footer(btn("Export .mol2"), btn("Apply to topology", True)))
    tb = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("atom", "Charges", True, "Colour by charge"), tbtn("tag", "Labels", False, "Show q"), spacer())
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">'
              f'<div style="position: relative; flex-grow: 1">{view}{hud(chip("PS · atactic · DP 8", TEXT, "#16191Ccc"), chip("colour = partial charge", MUTED, "#16191Ccc"))}</div>{strip}</div>')
    return studio_shell("PS_atactic_DP8.caps", "CAPS Field — partial charges", tb, center, right,
                        (f"<span>{len(atoms)} atoms · net {net:+.1e} e</span>", "<span>Gasteiger–Marsili · 6 iterations</span>"))


# ---------------------------------------------------------------- Periodic images
def pbc_cell():
    L = 28.0
    chains = []
    for k in range(7):
        rnd = random.Random(900 + k)
        st = (rnd.uniform(3, L - 3), rnd.uniform(3, L - 3), rnd.uniform(3, L - 3))
        chains.append(random_chain(st, 70, 950 + k, bond=1.54, box=None, persistence=0.7))
    return L, chains


def wrap(p, L):
    return tuple(x % L for x in p)


def split_wrapped(pts, L):
    segs, cur = [], [wrap(pts[0], L)]
    for p in pts[1:]:
        w = wrap(p, L)
        if dist(w, cur[-1]) > L / 2:
            segs.append(cur); cur = [w]
        else:
            cur.append(w)
    segs.append(cur)
    return segs


def periodic():
    from display import style_bar
    L, chains = pbc_cell()
    cols = ["#F0A83C", "#6CC4D8", "#DE775D", "#9B7AD5", "#7DC884", "#D6AC5C", "#2271DB"]
    crossing = sum(1 for ch in chains if any(not 0 <= x < L for p in ch for x in p))
    pieces = sum(len(split_wrapped(ch, L)) for ch in chains)
    w2 = (VW - 36) // 2
    h2 = VH - 150

    def view(sid, mode):
        sc = Scene(sid, w2, h2 - 30, yaw=0.6, pitch=0.4, persp=0.25, fog=0.45)
        for k, ch in enumerate(chains):
            if mode == "wrapped":
                for s in split_wrapped(ch, L):
                    if len(s) > 1:
                        sc.add_tube(s, cols[k], 0.55)
            else:
                sc.add_tube(ch, cols[k], 0.55)
        sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.2, None, 0.9)
        if mode == "whole":
            sc.center = (L / 2, L / 2, L / 2)
        return sc.svg()
    # minimum-image example: first and last bead of the chain that crosses most
    ch = max(chains, key=lambda c: max(abs(x - L / 2) for p in c for x in p))
    a, b = wrap(ch[0], L), wrap(ch[-1], L)
    raw = dist(a, b)
    d = [b[k] - a[k] for k in range(3)]
    d = [x - L * round(x / L) for x in d]
    mi = math.sqrt(sum(x * x for x in d))
    true = dist(ch[0], ch[-1])
    panel = lambda title, v, note: (f'<div style="width: {w2}px; display: flex; flex-direction: column; gap: 6px"><div style="display: flex; align-items: baseline; gap: 10px"><span style="font-size: 13px; font-weight: 600">{title}</span><span style="font-size: 11.5px; color: {DIM}">{note}</span></div>'
                                    f'<div style="position: relative; background: {BG0}; border: 1px solid {LINE}; border-radius: 8px; overflow: hidden; height: {h2 - 30}px">{v}</div></div>')
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">'
              f'<div style="display: flex; gap: 12px; padding: 14px 12px 0">{panel("Wrapped · as stored", view("pw", "wrapped"), f"{pieces} pieces in the box")}{panel("Whole molecules", view("pu", "whole") + style_bar("Backbone", f"Backbone · H hidden · {len(chains) * 212:,} atoms in model (PE, C₇₀H₁₄₂)".replace(",", " "), "right: 10px; top: 10px"), "unwrapped with image flags")}</div>'
              f'<div style="flex-grow: 1; display: flex; gap: 24px; align-items: center; padding: 0 20px; margin-top: 12px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'{kv("Chains", str(len(chains)))}{kv("Cross a face", str(crossing))}{kv("Pieces when wrapped", str(pieces))}{kv("Box", f"{L:.1f} Å cubic")}</div></div>')
    ib = next(i for i in range(len(ch) - 1) if dist(wrap(ch[i], L), wrap(ch[i + 1], L)) > L / 2)
    pa, pb = wrap(ch[ib], L), wrap(ch[ib + 1], L)
    dv = [pb[k] - pa[k] for k in range(3)]
    bmi = math.sqrt(sum((x - L * round(x / L)) ** 2 for x in dv))
    cn = chains.index(ch) + 1
    mt = table(["Chain " + str(cn), "Wrapped", "Min. img", "True"],
               [[f"Bond {ib + 1}–{ib + 2} (crosses)", f"{dist(pa, pb):.2f}", f"{bmi:.2f}", f"{dist(ch[ib], ch[ib + 1]):.2f}"],
                ["End to end", f"{raw:.2f}", f"{mi:.2f}", f"{true:.2f}"]],
               ["37%", "20%", "23%", "20%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=11.5, rowh=28)
    right = (panel_head("Periodic box", chip(f"{L:.0f} Å", MUTED, BG2, True))
             + section("Show", col(seg(["Wrapped", "Whole", "Images"], "Whole", full=True), row(select("Images", "3 × 3 × 1"), select("Image style", "faded 30 %"), gap=8), gap=8))
             + section("Wrap on save", col(seg(["Atoms", "Molecules", "Off"], "Molecules", full=True), f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Image flags are kept either way, so unwrapping is always exact.</span>', gap=8))
             + section("Measuring across the boundary", col(mt, f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Distances in Å. The minimum image fixes the bond, but gives the wrong end-to-end distance for a chain longer than half the box. CAPS measures chains on the unwrapped copy.</span>', gap=8))
             + footer(btn("Make molecules whole"), btn("Centre on selection", True)))
    tb = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("cube", "Box", True, "Periodic box"), tbtn("layers", "Images", False, "Images"), spacer())
    return studio_shell("PE_melt_7chains.caps", "CAPS Studio — periodic box", tb, center, right,
                        (f"<span>{len(chains)} chains · {sum(map(len, chains))} beads</span>", "<span>periodic x y z</span>"))


# ---------------------------------------------------------------- Orientation order
def jacobi3(A):
    A = [r[:] for r in A]
    for _ in range(60):
        p, q = max(((i, j) for i in range(3) for j in range(i + 1, 3)), key=lambda t: abs(A[t[0]][t[1]]))
        if abs(A[p][q]) < 1e-12:
            break
        th = 0.5 * math.atan2(2 * A[p][q], A[q][q] - A[p][p])
        c, s = math.cos(th), math.sin(th)
        for k in range(3):
            A[k][p], A[k][q] = c * A[k][p] - s * A[k][q], s * A[k][p] + c * A[k][q]
        for k in range(3):
            A[p][k], A[q][k] = c * A[p][k] - s * A[q][k], s * A[p][k] + c * A[q][k]
    return sorted(A[i][i] for i in range(3))


def order_S(chains):
    Q = [[0.0] * 3 for _ in range(3)]
    n = 0
    for ch in chains:
        for i in range(len(ch) - 2):
            u = norm(sub(ch[i + 2], ch[i]))
            for a in range(3):
                for b in range(3):
                    Q[a][b] += 1.5 * u[a] * u[b] - (0.5 if a == b else 0)
            n += 1
    Q = [[x / n for x in r] for r in Q]
    return jacobi3(Q)[-1], n


def orientation():
    na, nb, nc = 3, 4, 8
    cr = pe_crystal(na, nb, nc)
    per = 2 * nc
    xchains = [[a["p"] for a in cr[i:i + per]] for i in range(0, len(cr), per)]
    S_x, n_x = order_S(xchains)
    L = 30.0
    am = [random_chain((random.Random(70 + k).uniform(3, L - 3), random.Random(80 + k).uniform(3, L - 3), random.Random(90 + k).uniform(3, L - 3)), 60, 1000 + k,
                       bond=1.54, box=((0.5, 0.5, 0.5), (L - 0.5, L - 0.5, L - 0.5)), persistence=0.62) for k in range(12)]
    S_a, n_a = order_S(am)
    sw = 330
    s1 = Scene("ox", sw, 250, yaw=0.5, pitch=0.25, persp=0.2, fog=0.4, atom_k=0.9)
    s1.add_atoms(cr)
    s2 = Scene("oa", sw, 250, yaw=0.6, pitch=0.4, persp=0.25, fog=0.45)
    for k, ch in enumerate(am):
        s2.add_tube(ch, mix(ACC, "#8B969E", (k % 4) / 5), 0.5)
    s2.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)

    def cell(title, svg, S, n, note):
        return (f'<div style="display: flex; flex-direction: column; gap: 8px; padding: 14px; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px">'
                f'<span style="font-size: 13px; font-weight: 600">{title}</span><div style="background: {BG0}; border-radius: 8px; overflow: hidden">{svg}</div>'
                f'<div style="display: flex; align-items: baseline; gap: 8px"><span style="font-family: {MONO}; font-size: 26px">{S:.3f}</span><span style="font-size: 12px; color: {MUTED}">S from {n} chord vectors</span></div>'
                f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">{note}</span></div>')
    # illustrative stretching curve
    strain = [i / 50 for i in range(0, 101)]
    curve = [(e, 0.02 + 0.55 * (1 - math.exp(-e / 0.7))) for e in strain]
    pl = plot(520, 220, [(curve, ACC, 2, None)], (0, 2), (0, 1), [0, 0.5, 1, 1.5, 2], [0, 0.5, 1], "true strain", "S (chord vectors)", pad=(44, 10, 14, 30))
    left = (f'<div style="display: flex; gap: 14px">'
            + cell("PE crystal · Bunn cell 3 × 4 × 8", s1.svg(), S_x, n_x, "Chords i → i+2 lie on the c axis, so S = 1 exactly: a check that the estimator is right.")
            + cell("Amorphous cell · 12 chains", s2.svg(), S_a, n_a, "Random chains give S near 0; the residual is finite-sample noise from so few chords.")
            + '</div>')
    defs = [("S = λ<sub>max</sub> of Q", "Q = ⟨3/2 u u − 1/2 I⟩ over chord vectors u"), ("Director", "eigenvector of λ<sub>max</sub>; shown as an arrow in Studio"), ("Local crystallinity", "fraction of chords with ≥ 8 neighbours aligned within 10° · [result]"), ("Herman's f", "⟨P₂(cos θ)⟩ against a chosen axis, e.g. draw direction")]
    dl = "".join(f'<div style="display: flex; gap: 12px; padding: 6px 0; border-bottom: 1px solid {BG2}; font-size: 12px"><span style="width: 150px; flex-shrink: 0">{a}</span><span style="color: {MUTED}">{b}</span></div>' for a, b in defs)
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("During uniaxial draw", col(pl, f'<div style="font-size: 11.5px; color: {DIM}">Illustrative shape. From a Mechanics run CAPS plots S and Herman\'s f per frame.</div>', gap=6), chip("illustrative", ERR, BG2), 14)
             + card("Definitions", dl, "", 14) + '</div>')
    content = (pagehead("Orientation & crystallinity", "Nematic order parameter from chord vectors, per frame or per region", row(btn("Export CSV", ic="download"), btn("Run on trajectory", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden"><div style="width: {2 * sw + 14 + 60}px; flex-shrink: 0">{left}</div>{right}</div>')
    return shell("Analyze", "CAPS Analyze — orientation", "Analyze › Orientation", content, ("<span>two cells · computed from the drawn chains</span>", "<span>chord vectors i → i+2</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage8/project", exist_ok=True)
    for name, fn in (("Charges", charges), ("PeriodicBox", periodic), ("Orientation", orientation)):
        open(f"stage8/project/{name}.dc.html", "w").write(fn())
        print(name, "ok")
