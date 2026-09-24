import math, random
from lib import *
from mols import *
from screen_app import shell, pagehead, card
from screen_studio2 import studio_shell, bar, spacer
from screen_builders import footer, hud


# ---------------------------------------------------------------- Tacticity statistics (Bernoullian)
def tacticity():
    Pm, DP = 0.50, 200
    rnd = random.Random(99)
    dy = ["m" if rnd.random() < Pm else "r" for _ in range(DP - 1)]
    tri = {"mm": 0, "mr": 0, "rr": 0}
    for a, b in zip(dy, dy[1:]):
        k = a + b
        tri["mr" if k in ("mr", "rm") else k] += 1
    nt = sum(tri.values())
    model = {"mm": Pm * Pm, "mr": 2 * Pm * (1 - Pm), "rr": (1 - Pm) ** 2}
    pent = ["mmmm", "mmmr", "rmmr", "mmrr", "mmrm", "rmrm", "rmrr", "mrrm", "mrrr", "rrrr"]
    def pm(s):
        mcount = s.count("m"); return Pm ** mcount * (1 - Pm) ** (4 - mcount)
    mult = {"mmmm": 1, "mmmr": 2, "rmmr": 1, "mmrr": 2, "mmrm": 2, "rmrm": 2, "rmrr": 2, "mrrm": 1, "mrrr": 2, "rrrr": 1}
    pmodel = {p: mult[p] * pm(p) for p in pent}
    assert abs(sum(pmodel.values()) - 1) < 1e-12
    pc = {p: 0 for p in pent}
    canon = {}
    for p in pent:
        canon[p] = p; canon[p[::-1]] = p
    for i in range(len(dy) - 3):
        pc[canon["".join(dy[i:i + 4])]] += 1
    npn = sum(pc.values())
    W, H = 600, 230
    bw = (W - 70) / len(pent)
    ymax = max(max(pmodel.values()), max(v / npn for v in pc.values())) * 1.15
    bars = ""
    for i, p in enumerate(pent):
        x = 50 + i * bw
        hm = pmodel[p] / ymax * (H - 60); hs = pc[p] / npn / ymax * (H - 60)
        bars += (f'<rect x="{x + 4:.1f}" y="{H - 36 - hm:.1f}" width="{bw / 2 - 5:.1f}" height="{hm:.1f}" rx="2" fill="{ACC}" opacity="0.85"></rect>'
                 f'<rect x="{x + bw / 2:.1f}" y="{H - 36 - hs:.1f}" width="{bw / 2 - 5:.1f}" height="{hs:.1f}" rx="2" fill="{SEL}" opacity="0.85"></rect>'
                 f'<text x="{x + bw / 2:.1f}" y="{H - 20}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{p}</text>')
    for v in (0, round(ymax / 2, 2)):
        y = H - 36 - v / ymax * (H - 60)
        bars += f'<line x1="46" y1="{y:.1f}" x2="{W - 10}" y2="{y:.1f}" stroke="{BG3}"></line><text x="40" y="{y + 3.5:.1f}" text-anchor="end" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{v:g}</text>'
    chart = f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="Pentad fractions, model and sample" style="display: block">{bars}</svg>'
    strip = "".join(f'<span title="{d}" style="width: 4px; height: 28px; border-radius: 1px; background: {ACC if d == "m" else "#2271DB"}"></span>' for d in dy)
    tt = table(["Triad", "Model", f"Chain · {nt}"], [[k, f"{model[k]:.3f}", f"{tri[k] / nt:.3f}"] for k in ("mm", "mr", "rr")], ["34%", "33%", "33%"], mono_cols=(0, 1, 2), align_right=(1, 2), fs=12, rowh=27)
    atoms, backbone, stereo = polystyrene(8, "atactic", seed=99)
    sc = Scene("tc", 440, 250, yaw=-0.12, pitch=0.5, persp=0.2, fog=0.4)
    sc.add_atoms(atoms)
    left = (f'<div style="width: 470px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("First 8 units", f'<div style="background: {BG0}; border-radius: 8px; overflow: hidden">{sc.svg()}</div>', chip("atactic", MUTED, BG2), 14)
            + card("Triads", col(tt, row(field("P<sub>m</sub>", f"{Pm:.2f}"), select("Model", "Bernoulli"), field("DP", str(DP)), gap=8), gap=10), chip("computed", OK, BG2), 14) + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card(f"Dyad sequence · {len(dy)} dyads", col(f'<div style="display: flex; gap: 1px; flex-wrap: wrap">{strip}</div>', row(chip(f'{dot(ACC)} m (meso) · {dy.count("m")}'), chip(f'{dot("#2271DB")} r (racemo) · {dy.count("r")}'), gap=6), gap=10), chip("seed 99", MUTED, BG2, True), 14)
             + card("Pentads", col(chart, row(chip(f'{dot(ACC)} Bernoulli model'), chip(f'{dot(SEL)} this chain'), gap=6),
                                   f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Pentads are what ¹³C NMR resolves. Enter measured pentads and CAPS fits P<sub>m</sub>, or first-order Markov P<sub>m/r</sub> and P<sub>r/m</sub> if Bernoulli fails the test mm·rr = (mr/2)².</div>', gap=8), "", 14)
             + '</div>')
    content = (pagehead("Tacticity statistics", "Stereo sequence of the built chain against the statistical model", row(btn("Redraw"), btn("Apply to chains", True), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Studio", "CAPS — tacticity statistics", "Studio › Polymer builder › Tacticity", content, (f"<span>PS · DP {DP} · P<sub>m</sub> = {Pm}</span>", "<span>model values exact · chain values counted</span>"))


# ---------------------------------------------------------------- Blend phase diagram (Flory–Huggins)
def fh_parts(NA, NB):
    def f(p, chi): return p / NA * math.log(p) + (1 - p) / NB * math.log(1 - p) + chi * p * (1 - p)
    def d1(p, chi): return (math.log(p) + 1) / NA - (math.log(1 - p) + 1) / NB + chi * (1 - 2 * p)
    return f, d1


def bisect(fn, a, b, it=80):
    fa = fn(a)
    for _ in range(it):
        m = 0.5 * (a + b); fm = fn(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def spinodal_pts(NA, NB, chi):
    # 1/(NA p) + 1/(NB (1-p)) = 2 chi
    g = lambda p: 1 / (NA * p) + 1 / (NB * (1 - p)) - 2 * chi
    pc = math.sqrt(NB) / (math.sqrt(NA) + math.sqrt(NB))
    return bisect(g, 1e-9, pc), bisect(g, pc, 1 - 1e-9)


def binodal_pts(NA, NB, chi):
    f, d1 = fh_parts(NA, NB)
    s1, s2 = spinodal_pts(NA, NB, chi)
    def partner(p1):
        target = d1(p1, chi)
        return bisect(lambda q: d1(q, chi) - target, s2, 1 - 1e-12)
    def tangent(p1):
        p2 = partner(p1)
        return (f(p2, chi) - f(p1, chi)) - d1(p1, chi) * (p2 - p1)
    p1 = bisect(tangent, 1e-12, s1)
    return p1, partner(p1)


def blend_phase():
    NA, NB, A, B = 100, 200, -0.02, 15.0
    chic = 0.5 * (1 / math.sqrt(NA) + 1 / math.sqrt(NB)) ** 2
    phic = math.sqrt(NB) / (math.sqrt(NA) + math.sqrt(NB))
    Tc = B / (chic - A)
    T_of = lambda chi: B / (chi - A)
    spin, bino = [], []
    for i in range(1, 60):
        chi = chic * (1 + 0.02 * i * i / 10)
        T = T_of(chi)
        if T < 250: break
        s1, s2 = spinodal_pts(NA, NB, chi); b1, b2 = binodal_pts(NA, NB, chi)
        spin.append((T, s1, s2)); bino.append((T, b1, b2))
    sp_line = [(s1, T) for T, s1, s2 in reversed(spin)] + [(phic, Tc)] + [(s2, T) for T, s1, s2 in spin]
    bi_line = [(b1, T) for T, b1, b2 in reversed(bino)] + [(phic, Tc)] + [(b2, T) for T, b1, b2 in bino]
    tmin = 250
    tmax = math.ceil(Tc / 50) * 50 + 50
    pl = plot(640, 380, [(bi_line, ACC, 2.2, None), (sp_line, SEL, 1.6, "6 4")], (0, 1), (tmin, tmax), [0, 0.25, 0.5, 0.75, 1], list(range(tmin, tmax + 1, 50)),
              "φ<tspan>A</tspan> (volume fraction of A)", "T (K)", pad=(52, 12, 16, 34), markers=[(phic, Tc, TEXT)])
    pl = pl.replace("φ<tspan>A</tspan>", "φA")
    T300 = 300.0
    chi300 = A + B / T300
    b1, b2 = binodal_pts(NA, NB, chi300)
    s1, s2 = spinodal_pts(NA, NB, chi300)
    res = table(["Quantity", "Value"], [["χ<sub>c</sub> = ½(N<sub>A</sub><sup>−½</sup> + N<sub>B</sub><sup>−½</sup>)²", f"{chic:.5f}"], ["φ<sub>c</sub> = √N<sub>B</sub> / (√N<sub>A</sub> + √N<sub>B</sub>)", f"{phic:.4f}"],
                                        ["T<sub>c</sub> = B / (χ<sub>c</sub> − A)", f"{Tc:.1f} K"], ["χ at 300 K", f"{chi300:.4f}"],
                                        ["Coexisting at 300 K", f"{b1:.4f} · {b2:.4f}"], ["Spinodal at 300 K", f"{s1:.4f} · {s2:.4f}"]],
                ["62%", "38%"], mono_cols=(1,), align_right=(1,), fs=12, rowh=28)
    left = (f'<div style="width: 700px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Phase diagram", col(pl, row(chip(f'{dot(ACC)} binodal (coexistence)'), chip(f'{dot(SEL)} spinodal'), chip(f'{dot(TEXT)} critical point'), gap=6),
                                         f'<div style="font-size: 11.5px; color: {DIM}">Two phases below the binodal. χ(T) = A + B/T gives an upper critical solution temperature because B &gt; 0.</div>', gap=8), chip("computed", OK, BG2), 14) + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Inputs", col(row(field("N<sub>A</sub>", str(NA)), field("N<sub>B</sub>", str(NB)), gap=8), row(field("A", f"{A}"), field("B", f"{B}", "K"), gap=8),
                                  select("χ source", "entered · or fitted from MD at 3 temperatures"), gap=8), chip("example inputs", ACC, "#3A2C14"), 14)
             + card("Results", res, "", 14)
             + f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Binodal from equal chemical potentials (common tangent), solved numerically; spinodal from ∂²f/∂φ² = 0. A and B here are example inputs, not values for a real blend. Flory 1942; Huggins 1942.</div>'
             + '</div>')
    content = (pagehead("Blend phase diagram", "Flory–Huggins coexistence and stability for a binary polymer blend", row(btn("Export CSV", ic="download"), btn("Fit χ(T) from MD", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Analyze", "CAPS Analyze — blend phase diagram", "Analyze › Blend phase diagram", content, (f"<span>N<sub>A</sub> = {NA} · N<sub>B</sub> = {NB}</span>", "<span>lattice reference volume = 1 segment</span>"))


# ---------------------------------------------------------------- Electrostatics (Ewald / PME)
def good_fft(n):
    while True:
        m = n
        for p in (2, 3, 5, 7):
            while m % p == 0: m //= p
        if m == 1: return n
        n += 1


def electrostatics():
    L, rc, tol, sp = (45.3, 45.3, 45.3), 12.0, 1e-5, 1.2
    beta = bisect(lambda b: math.erfc(b * rc) - tol, 0.01, 2.0)
    grid = [good_fft(math.ceil(x / sp)) for x in L]
    real_sp = [x / n for x, n in zip(L, grid)]
    tols = [(t, bisect(lambda b: math.erfc(b * rc) - t, 0.01, 2.0)) for t in (1e-4, 1e-5, 1e-6)]
    curve = [(r / 10, math.erfc(beta * r / 10)) for r in range(5, 161)]
    curve2 = [(r / 10, 1.0) for r in range(5, 161)]
    pl = plot(600, 250, [(curve, ACC, 2.2, None), ([(rc, 0), (rc, 1)], DIM, 1.2, "4 4")], (0, 16), (0, 1), [0, 4, 8, 12, 16], [0, 0.5, 1],
              "r (Å)", "erfc(βr) · real-space share", pad=(46, 12, 16, 32), markers=[(rc, tol, ERR)])
    tt = table(["Tolerance erfc(βr<sub>c</sub>)", "β (Å⁻¹)", "βr<sub>c</sub>"], [["10⁻" + "⁰¹²³⁴⁵⁶⁷⁸⁹"[int(round(-math.log10(t)))], f"{b:.4f}", f"{b * rc:.3f}"] for t, b in tols],
               ["44%", "28%", "28%"], mono_cols=(0, 1, 2), align_right=(1, 2), fs=12, rowh=27, hl={1})
    gt = table(["Axis", "L (Å)", "grid", "spacing (Å)"], [[a, f"{x:.1f}", str(n), f"{s:.4f}"] for a, x, n, s in zip("xyz", L, grid, real_sp)],
               ["16%", "28%", "24%", "32%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=12, rowh=27)
    left = (f'<div style="width: 660px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Real-space part of the Coulomb sum", col(pl, f'<div style="font-size: 11.5px; color: {DIM}">At the cut-off the real-space term has fallen to the tolerance (red); the rest is carried by the reciprocal-space mesh.</div>', gap=6), chip("computed", OK, BG2), 14)
            + card("Tolerance vs splitting parameter", tt, "", 14) + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Settings", col(row(select("Method", "PME"), field("Order", "4"), gap=8), row(field("Cut-off r<sub>c</sub>", f"{rc}", "Å"), field("Tolerance", "1e-5"), gap=8),
                                    row(field("Max grid spacing", f"{sp}", "Å"), select("LJ", "cut-off + tail correction"), gap=8), gap=8), "", 14)
             + card("Derived", col(kv("β", f"{beta:.4f} Å⁻¹"), kv("Mesh", " × ".join(map(str, grid))), gt, gap=8), chip("FFT sizes: factors 2·3·5·7", MUTED, BG2), 14)
             + f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">β solves erfc(β r<sub>c</sub>) = tolerance, the convention GROMACS uses for ewald-rtol. Grid sizes round up to FFT-friendly numbers. Essmann et al., <i>J. Chem. Phys.</i> 103, 8577 (1995).</div>'
             + '</div>')
    content = (pagehead("Dynamics · electrostatics", "Ewald splitting, cut-off and PME mesh, derived from one tolerance", row(btn("Reset to defaults"), btn("Apply", True), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Dynamics", "CAPS Dynamics — electrostatics", "Dynamics › Electrostatics", content, (f"<span>box {L[0]} Å cubic</span>", f"<span>β = {beta:.4f} Å⁻¹ · mesh {grid[0]}³</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage13/project", exist_ok=True)
    for name, fn in (("Tacticity", tacticity), ("BlendPhase", blend_phase), ("Electrostatics", electrostatics)):
        open(f"stage13/project/{name}.dc.html", "w").write(fn())
        print(name, "ok")
