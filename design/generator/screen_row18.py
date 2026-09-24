import math, random
from lib import *
from mols import *
from screen_app import shell, pagehead, card
from screen_studio2 import studio_shell, bar, spacer
from screen_builders import footer, hud
from display import style_bar, zoom_inset

M0_PS = 104.15
R_J = 8.314462618


# ---------------------------------------------------------------- Polydispersity (Schulz–Zimm)
def polydispersity():
    Nn, D, nch = 40, 1.10, 20
    k = 1 / (D - 1)
    rnd = random.Random(2026)
    sample = [max(1, round(rnd.gammavariate(k, Nn / k))) for _ in range(nch)]
    Mn_s = sum(n * M0_PS for n in sample) / nch
    Mw_s = sum((n * M0_PS) ** 2 for n in sample) / sum(n * M0_PS for n in sample)
    D_s = Mw_s / Mn_s
    xs = list(range(1, 91))
    pn = [x ** (k - 1) * math.exp(-k * x / Nn) for x in xs]
    zn = sum(pn); pn = [v / zn for v in pn]
    pw = [x * v for x, v in zip(xs, pn)]
    zw = sum(pw); pw = [v / zw for v in pw]
    ymax = max(max(pn), max(pw)) * 1.15
    pl = plot(600, 260, [(list(zip(xs, pn)), SEL, 2, None), (list(zip(xs, pw)), ACC, 2, None)], (0, 90), (0, ymax), [0, 20, 40, 60, 80], [0, round(ymax / 2, 3), round(ymax, 3)],
              "degree of polymerisation N", "fraction", pad=(52, 12, 16, 32), markers=[(n, 0.0015, "#E9ECEF") for n in sample])
    sc = Scene("pd", 420, 330, yaw=0.55, pitch=0.4, persp=0.25, fog=0.4)
    L = 34.0
    order = sorted(range(nch), key=lambda i: sample[i])
    for r, i in enumerate(order):
        rr = random.Random(700 + i)
        st = (rr.uniform(3, L - 3), rr.uniform(3, L - 3), rr.uniform(3, L - 3))
        pts = random_chain(st, 2 * sample[i], 800 + i, bond=1.54, box=((0.5, 0.5, 0.5), (L - 0.5, L - 0.5, L - 0.5)), persistence=0.6)
        sc.add_tube(pts, mix(SEL, ACC, r / (nch - 1)), 0.45)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    lens = "".join(f'<span style="min-width: 30px; text-align: center; padding: 3px 0; border-radius: 4px; background: {BG0}; border: 1px solid {LINE}; font-family: {MONO}; font-size: 11.5px">{n}</span>' for n in sorted(sample))
    tgt = table(["", "Target", "This sample"], [["N<sub>n</sub>", f"{Nn}", f"{Mn_s / M0_PS:.1f}"], ["M<sub>n</sub> (g/mol)", f"{Nn * M0_PS:,.0f}".replace(",", " "), f"{Mn_s:,.0f}".replace(",", " ")],
                                              ["M<sub>w</sub> (g/mol)", f"{Nn * M0_PS * D:,.0f}".replace(",", " "), f"{Mw_s:,.0f}".replace(",", " ")], ["Đ = M<sub>w</sub>/M<sub>n</sub>", f"{D:.2f}", f"{D_s:.3f}"],
                                              ["Shortest · longest", "—", f"{min(sample)} · {max(sample)}"]], ["40%", "30%", "30%"], mono_cols=(1, 2), align_right=(1, 2), fs=12, rowh=27)
    left = (f'<div style="width: 460px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Cell · 20 chains, short → long", f'<div style="position: relative; background: {BG0}; border-radius: 8px; overflow: hidden">{sc.svg()}'
                   + style_bar("Backbone", f"Backbone · H hidden · {sum(sample) * 16 + 2 * nch:,} atoms in model".replace(",", " "), "right: 8px; top: 8px")
                   + zoom_inset("pdz", 170, 132, 2, 5, "left: 8px; bottom: 8px", "All atoms") + '</div>', chip("cyan short · amber long", DIM, BG2), 14)
            + card("Drawn lengths (N)", f'<div style="display: flex; flex-wrap: wrap; gap: 5px">{lens}</div>', chip("seed 2026", MUTED, BG2, True), 14) + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Schulz–Zimm distribution", col(pl, row(chip(f'{dot(SEL)} number fraction'), chip(f'{dot(ACC)} weight fraction'), chip(f'{dot("#E9ECEF")} the 20 drawn chains'), gap=6), gap=8),
                    chip(f"k = 1/(Đ − 1) = {k:.0f}", MUTED, BG2, True), 14)
             + '<div style="display: flex; gap: 12px">'
             + card("Target vs sample", tgt, chip("computed", OK, BG2), 14, "flex: 1 1 0; min-width: 0")
             + card("Options", col(select("Distribution", "Schulz–Zimm"), row(field("N<sub>n</sub>", str(Nn)), field("Đ", f"{D:.2f}"), gap=8), select("Match", "draw, then report the sample"), gap=8), "", 14, "flex: 1 1 0; min-width: 0")
             + '</div>'
             + f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">With only 20 chains the sample Đ differs from the target; CAPS reports the sample values in the provenance, never the target as if achieved.</div>'
             + '</div>')
    content = (pagehead("Grow · polydispersity", "Chain lengths drawn from a distribution, reported as drawn", row(btn("Redraw", ic="rotate"), btn("Use these lengths", True), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Grow", "CAPS Grow — polydispersity", "Grow › Polydispersity", content, (f"<span>PS · {nch} chains · Σ N = {sum(sample)}</span>", f"<span>M<sub>0</sub> = {M0_PS} g/mol</span>"))


# ---------------------------------------------------------------- Copolymer (Mayo–Lewis, terminal model)
def copolymer():
    r1, r2, f1 = 0.52, 0.46, 0.50
    f2 = 1 - f1
    F1 = (r1 * f1 * f1 + f1 * f2) / (r1 * f1 * f1 + 2 * f1 * f2 + r2 * f2 * f2)
    az = (1 - r2) / (2 - r1 - r2)
    p11 = r1 * f1 / (r1 * f1 + f2)
    p22 = r2 * f2 / (r2 * f2 + f1)
    rnd = random.Random(12)
    n = 80
    seq = [1 if rnd.random() < F1 else 2]
    for _ in range(n - 1):
        stay = p11 if seq[-1] == 1 else p22
        seq.append(seq[-1] if rnd.random() < stay else 3 - seq[-1])
    F1_s = seq.count(1) / n
    runs = {1: [], 2: []}
    cur, ln = seq[0], 1
    for s in seq[1:]:
        if s == cur:
            ln += 1
        else:
            runs[cur].append(ln); cur, ln = s, 1
    runs[cur].append(ln)
    mrl1_th, mrl2_th = 1 / (1 - p11), 1 / (1 - p22)
    curve = [(x / 100, (r1 * (x / 100) ** 2 + (x / 100) * (1 - x / 100)) / (r1 * (x / 100) ** 2 + 2 * (x / 100) * (1 - x / 100) + r2 * (1 - x / 100) ** 2)) for x in range(0, 101)]
    pl = plot(330, 230, [(curve, ACC, 2, None), ([(0, 0), (1, 1)], DIM, 1, "4 4")], (0, 1), (0, 1), [0, 0.25, 0.5, 0.75, 1], [0, 0.25, 0.5, 0.75, 1], "f₁ (styrene in feed)", "F₁ (styrene in chain)",
              pad=(46, 12, 14, 32), markers=[(f1, F1, SEL), (az, az, OK)])
    CA, CB = "#F0A83C", "#2271DB"
    strip = "".join(f'<span title="{"S" if s == 1 else "M"}" style="width: 11px; height: 26px; border-radius: 2px; background: {CA if s == 1 else CB}"></span>' for s in seq)
    pts = random_chain((0.0, 0.0, 0.0), n, 31, bond=2.6, box=None, persistence=0.7)
    sc = Scene("cp", VW, 380, yaw=0.3, pitch=0.35, persp=0.2, fog=0.35, atom_k=1.0)
    sc.add_atoms([{"e": "C", "p": p, "r": 1.05, "c": CA if s == 1 else CB} for p, s in zip(pts, seq)], [(i, i + 1) for i in range(n - 1)])
    sc._prep()
    sc.scale *= 1.3
    view = sc.svg()
    stats = table(["", "Model", "This chain"], [["F₁ (styrene)", f"{F1:.3f}", f"{F1_s:.3f}"], ["Mean run, S", f"{mrl1_th:.2f}", f"{sum(runs[1]) / len(runs[1]):.2f}"],
                                             ["Mean run, MMA", f"{mrl2_th:.2f}", f"{sum(runs[2]) / len(runs[2]):.2f}"], ["Longest run", "—", f"{max(max(runs[1]), max(runs[2]))}"]],
                  ["44%", "28%", "28%"], mono_cols=(1, 2), align_right=(1, 2), fs=12, rowh=26)
    right = (panel_head("Copolymer", chip("terminal model", MUTED, BG2))
             + section("Monomers", col(row(select("M₁", "styrene"), select("M₂", "methyl methacrylate"), gap=8), row(field("r₁", f"{r1}"), field("r₂", f"{r2}"), field("f₁ feed", f"{f1:.2f}"), gap=8),
                                       cite("Reactivity ratios are inputs; the values shown are commonly quoted for S/MMA free-radical copolymerisation. Check your source."), gap=8))
             + section("Composition", f'<div style="display: flex; justify-content: center">{pl}</div>')
             + footer(btn("Redraw"), btn("Build 20 chains", True)))
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">'
              f'<div style="position: relative">{view}{hud(chip("P(S-co-MMA) · 80 units", TEXT, "#16191Ccc"), chip(f"{dot(CA)} styrene · {dot(CB)} MMA", MUTED, "#16191Ccc"))}</div>'
              f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 10px; padding: 16px 20px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Sequence, head → tail</h3>'
              f'<div style="display: flex; gap: 1px">{strip}</div>'
              f'<div style="display: flex; gap: 20px; align-items: flex-start"><div style="width: 420px; flex-shrink: 0">{stats}</div>'
              f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Mayo–Lewis: F₁ = (r₁f₁² + f₁f₂) / (r₁f₁² + 2f₁f₂ + r₂f₂²) = {F1:.3f}. Azeotrope at f₁ = (1 − r₂)/(2 − r₁ − r₂) = {az:.3f}. '
              f'Low-conversion model: composition drift is not included. Mayo &amp; Lewis, <i>J. Am. Chem. Soc.</i> 66, 1594 (1944).</div></div></div></div>')
    tb = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("hex", "Copolymer", True, "Copolymer builder"), spacer())
    return studio_shell("PS_co_PMMA.caps", "CAPS Studio — copolymer builder", tb, center, right,
                        (f"<span>{n} units · {seq.count(1)} S · {seq.count(2)} MMA</span>", "<span>sequence from the terminal model, seed 12</span>"))


VW = 1440 - 72 - 380


# ---------------------------------------------------------------- Solvent screen (Hildebrand / Flory–Huggins)
SOLV = [("Toluene", 92.14, 0.867, 18.2, "good solvent"), ("Chloroform", 119.38, 1.489, 19.0, "good solvent"), ("THF", 72.11, 0.889, 19.4, "good solvent"),
        ("Cyclohexane", 84.16, 0.779, 16.8, "Θ solvent near 34.5 °C"), ("Acetone", 58.08, 0.791, 20.0, "non-solvent"), ("n-Hexane", 86.18, 0.655, 14.9, "non-solvent"),
        ("Ethanol", 46.07, 0.789, 26.5, "non-solvent"), ("Water", 18.015, 0.997, 47.8, "non-solvent")]


def solvent_screen():
    dP, T = 18.6, 298.15
    rows, bars = [], []
    for name, M, rho, d, known in SOLV:
        V = M / rho
        chi = V * (d - dP) ** 2 / (R_J * T) + 0.34
        pred = "borderline" if abs(chi - 0.5) < 0.05 else ("solvent" if chi < 0.5 else "non-solvent")
        agree = (pred != "non-solvent") == ("good" in known or "Θ" in known)
        rows.append((name, V, d, chi, pred, known, agree))
    rows.sort(key=lambda r: r[3])
    t = table(["Solvent", "V (cm³/mol)", "δ (MPa½)", "χ", "Predicted", "Known for PS", ""],
              [[n, f"{V:.1f}", f"{d:.1f}", f"{c:.2f}", p, k, icon("check", 14, OK) if a else icon("alert", 14, ACC)] for n, V, d, c, p, k, a in rows],
              ["16%", "14%", "12%", "10%", "15%", "27%", "6%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=12, rowh=30)
    W, H = 560, 300
    xmax = 3.0
    bh = (H - 50) / len(rows)
    bs = "".join(f'<text x="96" y="{20 + i * bh + bh / 2 + 4:.1f}" text-anchor="end" font-size="11" fill="{MUTED}">{n}</text>'
                 f'<rect x="104" y="{20 + i * bh + 4:.1f}" width="{min(c, xmax) / xmax * (W - 130):.1f}" height="{bh - 8:.1f}" rx="3" fill="{OK if c < 0.5 else ERR}" opacity="0.75"></rect>'
                 + (f'<text x="{110 + c / xmax * (W - 130):.1f}" y="{20 + i * bh + bh / 2 + 4:.1f}" font-size="10.5" font-family="IBM Plex Mono" fill="{TEXT}">{c:.2f}</text>' if c <= xmax else f'<text x="{98 + (W - 130):.1f}" y="{20 + i * bh + bh / 2 + 4:.1f}" text-anchor="end" font-size="10.5" font-family="IBM Plex Mono" fill="#0F1113">{c:.2f} →</text>')
                 for i, (n, V, d, c, p, k, a) in enumerate(rows))
    x05 = 104 + 0.5 / xmax * (W - 130)
    chart = (f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="Flory–Huggins chi per solvent" style="display: block">{bs}'
             f'<line x1="{x05:.1f}" y1="12" x2="{x05:.1f}" y2="{H - 26}" stroke="{TEXT}" stroke-dasharray="4 4"></line><text x="{x05 + 4:.1f}" y="{H - 12}" font-size="10.5" fill="{MUTED}">χ = 0.5</text></svg>')
    miss = [r[0] for r in rows if not r[6]]
    V_tol = 92.14 / 0.867
    left = (f'<div style="width: 620px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("χ against polystyrene", chart, chip(f"δ(PS) = {dP} MPa½ · {T} K", MUTED, BG2, True), 14)
            + card("How χ is estimated", f'<div style="font-family: {MONO}; font-size: 12.5px; line-height: 1.8">χ ≈ V<sub>s</sub>(δ<sub>s</sub> − δ<sub>p</sub>)² / RT + 0.34<br>toluene: {V_tol:.1f} × {(18.2 - dP) ** 2:.2f} / {R_J * T:.0f} + 0.34 = {V_tol * (18.2 - dP) ** 2 / (R_J * T) + 0.34:.3f}</div>', "", 14)
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Solvents ranked", t, chip("computed", OK, BG2), 14)
             + card("Where this estimate fails", f'<div style="font-size: 12.5px; line-height: 1.6">Hildebrand parameters ignore polarity and hydrogen bonding. Here the estimate calls <b>{", ".join(miss) or "none"}</b> a solvent, but it does not dissolve polystyrene. CAPS shows this check beside every screen and offers Hansen parameters or a direct χ from MD when polarity matters.</div>', chip(f"{len(miss)} mismatch", ACC, "#3A2C14"), 14)
             + f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Inputs: δ values from Barton, <i>CRC Handbook of Solubility Parameters</i> (1991); densities at 20–25 °C; V = M/ρ. δ(PS) is an input (literature 17.4–19.0).</div>'
             + '</div>')
    content = (pagehead("Solvent screen", "Flory–Huggins χ from solubility parameters, checked against known behaviour", row(select("Polymer", "polystyrene"), btn("Compute χ by MD", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Analyze", "CAPS Analyze — solvent screen", "Analyze › Solvent screen", content, ("<span>8 solvents · 298.15 K</span>", "<span>χ ≈ V(Δδ)²/RT + 0.34</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage12/project", exist_ok=True)
    for name, fn in (("Polydispersity", polydispersity), ("Copolymer", copolymer), ("SolventScreen", solvent_screen)):
        open(f"stage12/project/{name}.dc.html", "w").write(fn())
        print(name, "ok")
