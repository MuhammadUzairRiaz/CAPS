import math
from lib import *
from mols import *
from screen_app import shell, pagehead, card
from screen_studio2 import studio_shell, bar, spacer
from screen_builders import footer, hud

VW = 1440 - 72 - 380


# ---------------------------------------------------------------- geometry: Horn quaternion superposition
def jacobi4(A):
    A = [r[:] for r in A]
    V = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    for _ in range(100):
        p, q = max(((i, j) for i in range(4) for j in range(i + 1, 4)), key=lambda t: abs(A[t[0]][t[1]]))
        if abs(A[p][q]) < 1e-12:
            break
        th = 0.5 * math.atan2(2 * A[p][q], A[q][q] - A[p][p])
        c, s = math.cos(th), math.sin(th)
        for k in range(4):
            akp, akq = A[k][p], A[k][q]
            A[k][p], A[k][q] = c * akp - s * akq, s * akp + c * akq
        for k in range(4):
            apk, aqk = A[p][k], A[q][k]
            A[p][k], A[q][k] = c * apk - s * aqk, s * apk + c * aqk
        for k in range(4):
            vkp, vkq = V[k][p], V[k][q]
            V[k][p], V[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    return [A[i][i] for i in range(4)], V


def superpose(P, Q, idx=None):
    """Rotate/translate Q onto P, fitting on atoms idx. Returns moved Q."""
    idx = idx or list(range(len(P)))
    cp = mul(tuple(map(sum, zip(*[P[i] for i in idx]))), 1 / len(idx))
    cq = mul(tuple(map(sum, zip(*[Q[i] for i in idx]))), 1 / len(idx))
    S = [[0.0] * 3 for _ in range(3)]
    for i in idx:
        a, b = sub(Q[i], cq), sub(P[i], cp)
        for r in range(3):
            for c in range(3):
                S[r][c] += a[r] * b[c]
        (Sxx, Sxy, Sxz), (Syx, Syy, Syz), (Szx, Szy, Szz) = S
    N = [[Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
         [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
         [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
         [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz]]
    w, V = jacobi4(N)
    k = max(range(4), key=lambda i: w[i])
    q0, qx, qy, qz = (V[r][k] for r in range(4))
    R = [[q0*q0+qx*qx-qy*qy-qz*qz, 2*(qx*qy-q0*qz), 2*(qx*qz+q0*qy)],
         [2*(qx*qy+q0*qz), q0*q0-qx*qx+qy*qy-qz*qz, 2*(qy*qz-q0*qx)],
         [2*(qx*qz-q0*qy), 2*(qy*qz+q0*qx), q0*q0-qx*qx-qy*qy+qz*qz]]
    out = []
    for p in Q:
        d = sub(p, cq)
        out.append(add(cp, tuple(sum(R[r][c] * d[c] for c in range(3)) for r in range(3))))
    return out


def rmsd(P, Q, idx=None):
    idx = idx or list(range(len(P)))
    return math.sqrt(sum(dist(P[i], Q[i]) ** 2 for i in idx) / len(idx))


# ---------------------------------------------------------------- History & snapshots
def history():
    n10, n20 = len(polystyrene(10)[0]), len(polystyrene(20)[0])
    atoms, backbone, stereo = polystyrene(10, "isotactic", seed=4)
    sc = Scene("hi", VW, 900 - 44 - 48 - 26, yaw=-0.15, pitch=0.45, persp=0.2, fog=0.5)
    sc.add_atoms(atoms)
    view = sc.svg()
    steps = [("file", "New document", "—", "", "done"),
             ("hex", "Polymer builder · PS, DP 10, atactic", f"{n10} atoms", "", "done"),
             ("mirror", "Set tacticity · isotactic", f"{n10} atoms", "", "done"),
             ("relax", "Clean geometry · L-BFGS, GAFF2", f"{n10} atoms", "", "current"),
             ("plus", "Extend chain · DP 10 → 20", f"{n20} atoms", "undone · Redo or ⌘⇧Z", "undone"),
             ("relax", "Clean geometry · L-BFGS, GAFF2", f"{n20} atoms", "undone", "undone")]
    it = []
    for i, (ic, t, n, note, st) in enumerate(steps):
        cur = st == "current"
        colr = TEXT if st != "undone" else DIM
        it.append(f'<div style="display: flex; gap: 10px; align-items: flex-start; padding: 8px 10px; border-radius: 7px; background: {BG3 if cur else "transparent"}; border: 1px solid {ACC if cur else "transparent"}">'
                  f'<span style="width: 18px; text-align: right; font-family: {MONO}; font-size: 11px; color: {DIM}; padding-top: 2px">{i + 1}</span>{icon(ic, 16, ACC if cur else (MUTED if st == "done" else DIM))}'
                  f'<div style="display: flex; flex-direction: column; gap: 2px; min-width: 0; flex-grow: 1"><span style="font-size: 12.5px; color: {colr}{"; text-decoration: line-through" if st == "undone" else ""}">{t}</span>'
                  f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{n}{" · " + note if note else ""}</span></div>'
                  f'{chip("you are here", ACC, "#3A2C14") if cur else ""}</div>')
    snaps = []
    for sid, (dp, tact, name, when) in enumerate(((10, "atactic", "as built", "step 2"), (10, "isotactic", "iso · cleaned", "step 4"), (20, "isotactic", "DP 20 draft", "step 6 · on undone branch"))):
        a, _, _ = polystyrene(dp, tact, seed=4)
        s = Scene(f"hs{sid}", 96, 64, yaw=-0.15, pitch=0.45, fog=0.4, outline=False, atom_k=1.1)
        s.add_atoms(a)
        snaps.append(f'<div style="display: flex; gap: 10px; align-items: center; padding: 8px; border: 1px solid {LINE}; border-radius: 8px; background: {BG0}">'
                     f'<div style="width: 96px; height: 64px; flex-shrink: 0; border-radius: 5px; overflow: hidden; background: {BG0}">{s.svg()}</div>'
                     f'<div style="display: flex; flex-direction: column; gap: 3px; flex-grow: 1; min-width: 0"><span style="font-size: 12.5px; font-weight: 600">{name}</span>'
                     f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{len(a)} atoms · {when}</span></div>'
                     f'{btn("Compare", small=True, href="Compare.dc.html")}</div>')
    right = (panel_head("History", row(tbtn("undo", "Undo"), tbtn("redo", "Redo"), gap=2))
             + f'<div style="padding: 8px 8px 4px; display: flex; flex-direction: column; gap: 2px">{"".join(it)}</div>'
             + f'<div style="padding: 4px 18px 10px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Editing after an undo keeps the undone steps as a branch; nothing is lost until you delete it.</div>'
             + section("Snapshots", col(*snaps, gap=8), btn("Take snapshot", ic="plus", small=True))
             )
    tb = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("undo", "Undo", False, "Undo"), tbtn("redo", "Redo", False, "Redo"), sep(), tbtn("history", "History", True, "History"), spacer())
    center = f'<div style="flex-grow: 1; min-width: 0; position: relative; background: {BG0}">{view}{hud(chip("PS · isotactic · DP 10", TEXT, "#16191Ccc"), chip("step 4 of 6", MUTED, "#16191Ccc"))}</div>'
    return studio_shell("PS_iso_DP10.caps", "CAPS Studio — history & snapshots", tb, center, right,
                        (f"<span>{len(atoms)} atoms</span><span>history · 6 steps · 1 branch</span>", "<span>⌘Z undo · ⌘⇧Z redo</span>"))


# ---------------------------------------------------------------- Compare two snapshots
def compare():
    A, bbA, _ = polystyrene(8, "atactic", seed=6, curve=False)
    B, bbB, _ = polystyrene(8, "atactic", seed=6, curve=True)
    assert len(A) == len(B)
    PA, PB = [a["p"] for a in A], [b["p"] for b in B]
    fit_all = superpose(PA, PB)
    r_all = rmsd(PA, fit_all)
    heavy = [i for i, a in enumerate(A) if a["e"] != "H"]
    fit_bb = superpose(PA, PB, bbA)
    r_bb = rmsd(PA, fit_bb, bbA)
    r_heavy = rmsd(PA, superpose(PA, PB, heavy), heavy)
    dev = [dist(PA[i], fit_all[i]) for i in range(len(A))]
    dmax = max(dev)
    stops = [(0.0, SEL), (0.5, ACC), (1.0, ERR)]

    def ramp(t):
        for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
            if t <= t1:
                return mix(c0, c1, (t - t0) / (t1 - t0))
        return stops[-1][1]
    bonds = perceive_bonds(A)
    sc = Scene("cmp", VW, 900 - 44 - 48 - 26 - 180, yaw=0.0, pitch=0.35, persp=0.2, fog=0.35)
    for i, j in bonds:
        sc.add_line(PA[i], PA[j], MUTED, 1.4, None, 0.55)
    sc.add_atoms([{"e": b["e"], "p": fit_all[i], "c": ramp(dev[i] / dmax)} for i, b in enumerate(B)], bonds)
    view = sc.svg()
    worst = sorted(range(len(A)), key=lambda i: -dev[i])
    wl = [[f"{A[i]['e']}{i + 1}", "backbone" if i in bbA else ("ring" if A[i].get("ar") else "H"), f"{dev[i]:.2f}"] for i in worst[:6]]
    wt = table(["Atom", "Where", "Shift (Å)"], wl, ["30%", "40%", "30%"], mono_cols=(0, 2), align_right=(2,), fs=12, rowh=27)
    legend = (f'<div style="display: flex; flex-direction: column; gap: 4px"><div style="height: 8px; border-radius: 4px; background: linear-gradient(90deg, {SEL}, {ACC}, {ERR})"></div>'
              f'<div style="display: flex; justify-content: space-between; font-family: {MONO}; font-size: 10.5px; color: {DIM}"><span>0 Å</span><span>{dmax / 2:.1f}</span><span>{dmax:.1f} Å</span></div></div>')
    right = (panel_head("Compare snapshots", chip("same topology · 1:1 atom map", OK, BG2))
             + section("Pair", col(row(select("Reference (grey)", "straight · as built"), gap=8), row(select("Moving (coloured)", "curved · after edit"), gap=8), gap=8))
             + section("Fit on", seg(["All atoms", "Heavy", "Backbone"], "All atoms", full=True))
             + section("RMSD after superposition", col(kv("All atoms", f"{r_all:.2f} Å"), kv("Heavy atoms", f"{r_heavy:.2f} Å"), kv("Backbone C", f"{r_bb:.2f} Å"), kv("Atoms", f"{len(A)}"), gap=6))
             + section("Largest shifts", wt)
             + footer(btn("Export CSV", ic="download"), btn("Open both in split view", True)))
    tb = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("layers", "Compare", True, "Compare"), spacer(), chip("Colour = shift after fit", MUTED, BG2))
    strip = (f'<div style="height: 170px; flex-shrink: 0; display: flex; gap: 24px; padding: 16px 20px; background: {BG1}; border-top: 1px solid {LINE}">'
             f'<div style="width: 260px; display: flex; flex-direction: column; gap: 10px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Shift per atom</h3>{legend}'
             f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Grey sticks: reference. Balls: moving snapshot after fit.</span></div>'
             f'<div style="flex-grow: 1; min-width: 0">{plot(620, 138, [([(i + 1, d) for i, d in enumerate(dev)], ACC, 1.6, None)], (1, len(A)), (0, math.ceil(dmax)), [1, 50, 100, len(A)], list(range(0, math.ceil(dmax) + 1)), "atom index", "Å", pad=(30, 8, 10, 26))}</div></div>')
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">'
              f'<div style="position: relative; flex-grow: 1">{view}{hud(chip("PS · atactic · DP 8", TEXT, "#16191Ccc"), chip(f"RMSD {r_all:.2f} Å", ACC, "#16191Ccc", True))}</div>{strip}</div>')
    return studio_shell("PS_atactic_DP8.caps", "CAPS Studio — compare snapshots", tb, center, right,
                        (f"<span>{len(A)} atoms · 2 snapshots · values computed from the drawn chains</span>", "<span>fit: all atoms · Horn 1987 quaternion superposition</span>"))


# ---------------------------------------------------------------- Glass transition protocol
def glass_transition():
    T0, T1, dT, per_ns, reps, heat_ns = 500, 200, 20, 1.0, 3, 2.0
    temps = list(range(T0, T1 - 1, -dT))
    rate = dT / (per_ns * 1e-9)
    total = (len(temps) * per_ns + heat_ns) * reps
    stair = [(0, T0), (heat_ns, T0)]
    t = heat_ns
    for T in temps:
        stair += [(t, T), (t + per_ns, T)]
        t += per_ns
    tl = plot(560, 170, [(stair, SEL, 2, None)], (0, t), (150, 550), [0, 4, 8, 12, 16, round(t)], [200, 300, 400, 500], "time (ns) · one replica", "T (K)", pad=(40, 10, 14, 28))
    # illustrative bilinear shape only
    Tg_ill = 370
    v = lambda T: 0.95 + (2.0e-4 if T < Tg_ill else 5.6e-4) * (T - Tg_ill)
    pts = [(T, v(T)) for T in range(200, 501, 20)]
    glass = [(200, v(200)), (Tg_ill + 40, 0.95 + 2.0e-4 * 40)]
    melt = [(Tg_ill - 50, 0.95 - 5.6e-4 * 50), (500, v(500))]
    pl = plot(600, 300, [(glass, DIM, 1.4, "5 4"), (melt, DIM, 1.4, "5 4")], (190, 510), (0.88, 1.04), [200, 300, 400, 500], [0.9, 0.95, 1.0], "temperature (K)", "specific volume (cm³/g)",
              pad=(52, 12, 14, 32), markers=[(a, b, ACC) for a, b in pts])
    setup = col(row(field("Start", str(T0), "K"), field("End", str(T1), "K"), field("Step", str(dT), "K"), gap=8),
                row(field("Time per step", f"{per_ns:.0f}", "ns"), field("Anneal first", f"{heat_ns:.0f}", "ns"), field("Replicas", str(reps)), gap=8),
                row(select("Ensemble", "NPT · 1 atm"), select("Barostat", "MTTK"), gap=8),
                row(select("Property", "Specific volume"), select("Fit", "Two lines · free break"), gap=8), gap=10)
    derived = col(kv("Temperatures", f"{len(temps)}"), kv("Cooling rate", f"{rate:.1e} K/s".replace("e+", " × 10^").replace("^10", "¹⁰")), kv("Total simulated", f"{total:.0f} ns"), kv("Data per step", "last 50 % averaged"), gap=6)
    res = table(["Quantity", "Value"], [["T<sub>g</sub>", "[result] ± [result] K"], ["α glass", "[result] K⁻¹"], ["α melt", "[result] K⁻¹"], ["Fit residual", "[result]"]], ["55%", "45%"], mono_cols=(1,), fs=12, rowh=27)
    left = (f'<div style="width: 600px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
            + card("Cooling protocol", setup, chip("PS atactic · 20 chains", MUTED, BG2))
            + card("Schedule", col(tl, derived, gap=10), chip("computed from the settings", DIM, BG2))
            + '</div>')
    right_ = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px">'
              + card("Volume against temperature", col(pl, f'<div style="font-size: 11.5px; color: {DIM}">Illustrative shape only. CAPS plots the replica means with error bars and fits two lines; the crossing is T<sub>g</sub>.</div>', gap=8), chip("illustrative", ERR, BG2))
              + card("Result", col(res, f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Simulated cooling rates are many orders of magnitude faster than DSC, so T<sub>g</sub> comes out higher than experiment. CAPS reports the rate beside every T<sub>g</sub>. Method: Soldera &amp; Metatla, <i>Phys. Rev. E</i> 74, 061803 (2006).</div>', gap=8))
              + '</div>')
    content = (pagehead("Glass transition", "Stepwise NPT cooling, volume–temperature fit", row(btn("Save as recipe", ic="save"), btn("Queue 3 replicas", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right_}</div>')
    return shell("Equilibrate", "CAPS — glass transition", "Equilibrate › Glass transition", content,
                 ("<span>PS_atactic_DP40.caps</span>", f"<span>{reps} × {len(temps)} steps queued on local GPU</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage5/project", exist_ok=True)
    for name, fn in (("History", history), ("Compare", compare), ("GlassTransition", glass_transition)):
        open(f"stage5/project/{name}.dc.html", "w").write(fn())
        print(name, "ok")
