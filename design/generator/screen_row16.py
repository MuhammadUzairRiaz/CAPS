import math
from lib import *
from mols import *
from screen_app import shell, pagehead, card
from screen_studio2 import studio_shell, bar, spacer, KCAL_PER_K
from screen_builders import footer, hud

VW = 1440 - 72 - 380
VH = 900 - 44 - 48 - 26

# Larsen, Lin, Hart & Colina, Macromolecules 44, 6944 (2011): 21-step compression/relaxation
# (kind, T: "max"/"final", P fraction of Pmax or None, ps)
LARSEN = [("NVT", "max", None, 50), ("NVT", "final", None, 50), ("NPT", "final", 0.02, 50),
          ("NVT", "max", None, 50), ("NVT", "final", None, 100), ("NPT", "final", 0.6, 50),
          ("NVT", "max", None, 50), ("NVT", "final", None, 100), ("NPT", "final", 1.0, 50),
          ("NVT", "max", None, 50), ("NVT", "final", None, 100), ("NPT", "final", 0.5, 5),
          ("NVT", "max", None, 5), ("NVT", "final", None, 10), ("NPT", "final", 0.1, 5),
          ("NVT", "max", None, 5), ("NVT", "final", None, 10), ("NPT", "final", 0.01, 5),
          ("NVT", "max", None, 5), ("NVT", "final", None, 10), ("NPT", "final", "final", 800)]


# ---------------------------------------------------------------- Recipe editor
def recipe():
    Tmax, Tfin, Pmax, Pfin = 600, 300, 50000, 1
    total = sum(s[3] for s in LARSEN)
    # schedule plot: T and log P over time
    Tser, Pser, t = [], [], 0
    for kind, Tk, Pf, ps in LARSEN:
        T = Tmax if Tk == "max" else Tfin
        Tser += [(t, T), (t + ps, T)]
        if kind == "NPT":
            P = Pfin if Pf == "final" else Pf * Pmax
            Pser += [(t, math.log10(P)), (t + ps, math.log10(P))]
        t += ps
    tp = plot(620, 150, [(Tser, ACC, 1.8, None)], (0, total), (200, 700), [0, 400, 800, 1200, total], [300, 600], "", "T (K)", pad=(44, 8, 14, 20))
    pp = plot(620, 150, [([Pser[i], Pser[i + 1]], SEL, 2.4, None) for i in range(0, len(Pser), 2)], (0, total), (-0.5, 5), [0, 400, 800, 1200, total], [0, 1, 2, 3, 4, 5],
              "time (ps)", "log₁₀ P (bar)", pad=(44, 8, 14, 28))
    rows = []
    for i, (kind, Tk, Pf, ps) in enumerate(LARSEN):
        T = Tmax if Tk == "max" else Tfin
        P = "—" if Pf is None else (f"{Pfin}" if Pf == "final" else f"{Pf * Pmax:,.0f}".replace(",", " "))
        rows.append([str(i + 1), kind, f"{T}", P, f"{ps}"])
    tbl = table(["#", "Ens.", "T (K)", "P (bar)", "ps"], rows, ["10%", "18%", "22%", "30%", "20%"], mono_cols=(0, 2, 3, 4), align_right=(2, 3, 4), fs=10.5, rowh=17, hl={20})
    nodes = [("grow", "Grow", "20 chains · DP 40", False), ("relax", "Relax", "L-BFGS · 10⁻⁴", False), ("equil", "Equilibrate", "21-step · Larsen", True),
             ("dyn", "Production", "NPT 300 K · 10 ns", False), ("chart", "Analyze", "ρ · RDF · T<sub>g</sub>", False)]
    nd = []
    for i, (ic, n, s, on) in enumerate(nodes):
        nd.append(f'<div style="flex: 1 1 0; min-width: 0; display: flex; flex-direction: column; gap: 6px; padding: 10px 12px; border-radius: 9px; background: {BG3 if on else BG1}; border: 1px solid {ACC if on else LINE}">'
                  f'{icon(ic, 16, ACC if on else MUTED)}<span style="font-size: 13px; font-weight: 600">{n}</span>'
                  f'<span style="font-family: {MONO}; font-size: 11px; line-height: 1.4; color: {MUTED}">{s}</span></div>')
        if i < len(nodes) - 1:
            nd.append(f'<div style="width: 14px; flex-shrink: 0; height: 2px; background: {LINE}; position: relative"><span style="position: absolute; right: -2px; top: -4px; width: 0; height: 0; border-left: 7px solid {LINE}; border-top: 5px solid transparent; border-bottom: 5px solid transparent"></span></div>')
    pipe = f'<div style="display: flex; align-items: center; gap: 4px">{"".join(nd)}</div>'
    left_list = "".join(f'<div style="padding: 9px 12px; border-radius: 7px; background: {BG3 if on else "transparent"}; display: flex; flex-direction: column; gap: 2px"><span style="font-size: 12.5px; color: {TEXT if on else MUTED}">{a}</span><span style="font-family: {MONO}; font-size: 10.5px; color: {DIM}">{b}</span></div>'
                        for a, b, on in (("Amorphous polymer · standard", "5 steps · shared", True), ("Glass transition", "3 steps", False), ("Crosslinked epoxy", "6 steps", False),
                                         ("Solvated protein", "4 steps", False), ("Interface + adhesion", "5 steps", False)))
    side = (f'<aside style="width: 240px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-right: 1px solid {LINE}">{panel_head("Recipes", tbtn("plus", "New recipe"))}'
            f'<div style="padding: 8px; display: flex; flex-direction: column; gap: 2px">{left_list}</div>'
            f'<div style="padding: 10px 16px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Recipes are plain YAML in the project (same format as the CLI); runs record the exact recipe hash.</div></aside>')
    insp = (f'<div style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; gap: 10px">'
            + card("Step 3 · Equilibrate", col(row(field("T max", str(Tmax), "K"), field("T final", str(Tfin), "K"), field("P max", f"{Pmax:,}".replace(",", " "), "bar"), field("P final", str(Pfin), "bar"), gap=6),
                                               row(select("Thermostat", "Nosé–Hoover chain"), select("Barostat", "MTTK"), gap=8),
                                               cite("Larsen, Lin, Hart &amp; Colina, <i>Macromolecules</i> 44, 6944 (2011). Step values transcribed; check against the paper before release."), gap=10), chip(f"{len(LARSEN)} steps · {total} ps", ACC, "#3A2C14", True), 14)
            + card("Schedule", tbl, "", 12)
            + '</div>')
    mid = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
           + card("Pipeline", pipe, row(chip("valid", OK, BG2), gap=6), 14)
           + card("Temperature and pressure, step 3", col(tp, pp, f'<div style="font-size: 11.5px; color: {DIM}">Drawn from the schedule table; P is shown only for NPT steps.</div>', gap=4), "", 14)
           + '</div>')
    content = (pagehead("Recipe · Amorphous polymer, standard", "A reusable chain of steps; every run of it is logged with the recipe hash",
                        row(btn("Duplicate"), btn("Run on PS_atactic_DP40", True, "play"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; min-height: 0">{side}<div style="flex-grow: 1; min-width: 0; display: flex; gap: 14px; padding: 16px 18px; overflow: hidden">{mid}{insp}</div></div>')
    return shell("Jobs", "CAPS — recipe editor", "Jobs › Recipes", content, ("<span>recipes/amorphous_standard.yaml</span>", "<span>sha256 · shown after save</span>"))


# ---------------------------------------------------------------- Smart selection
def smart_select():
    atoms, backbone, stereo = polystyrene(10, "atactic", seed=8)
    arom = [i for i, a in enumerate(atoms) if a.get("ar")]
    ch2 = [k for k in backbone if k % 2 == 0]
    rings = [arom[i:i + 6] for i in range(0, len(arom), 6)]
    r3 = rings[4]
    cen = mul(tuple(map(sum, zip(*[atoms[i]["p"] for i in r3]))), 1 / 6)
    near = [i for i, a in enumerate(atoms) if dist(a["p"], cen) <= 5.0 and i not in r3]
    near_heavy = [i for i in near if atoms[i]["e"] != "H"]
    nH = sum(1 for a in atoms if a["e"] == "H")
    sc = Scene("ss", VW, VH - 64, yaw=-0.1, pitch=0.5, persp=0.2, fog=0.5)
    sc.add_atoms(atoms)
    sc.halo = set(arom)
    view = sc.svg()
    q = (f'<div style="height: 64px; flex-shrink: 0; display: flex; align-items: center; gap: 10px; padding: 0 16px; background: {BG1}; border-bottom: 1px solid {LINE}">'
         f'{icon("search", 16, MUTED)}<div style="flex-grow: 1; height: 36px; display: flex; align-items: center; padding: 0 12px; background: {BG0}; border: 1px solid {ACC}; border-radius: 7px; font-family: {MONO}; font-size: 13.5px">'
         f'<span style="color: {SEL}">smarts</span>&nbsp;<span>"c1ccccc1"</span>&nbsp;<span style="color: {MUTED}">and</span>&nbsp;<span style="color: {SEL}">chain</span>&nbsp;<span>1</span></div>'
         f'{chip(f"{len(arom)} atoms · {len(rings)} rings", ACC, "#3A2C14", True)}</div>')
    presets = [("Aromatic rings", 'smarts "c1ccccc1"', f"{len(arom)}"), ("Backbone CH₂", 'smarts "[CH2;!R]"', f"{len(ch2)}"),
               ("Stereo centres", "stereo *", f"{len(stereo)}"), ("Hydrogens", "element H", f"{nH}"),
               ("Near ring 5", "within 5 of ring 5", f"{len(near)}"), ("… heavy only", "within 5 of ring 5 and not H", f"{len(near_heavy)}")]
    pt = table(["Saved query", "Expression", "Atoms"], [[a, b, c] for a, b, c in presets], ["26%", "60%", "14%"], mono_cols=(1, 2), align_right=(2,), fs=11, rowh=28, hl={0})
    grammar = [("element C N O", "by element"), ("smarts \"…\"", "RDKit SMARTS match"), ("within 5.0 of sel", "distance in Å"), ("chain 1-4", "chains or residues"),
               ("stereo R|S|*", "stereo centres"), ("and · or · not · ( )", "combine")]
    gr = "".join(f'<div style="display: flex; gap: 10px; padding: 5px 0; border-bottom: 1px solid {BG2}; font-size: 12px"><span style="width: 170px; flex-shrink: 0; font-family: {MONO}; color: {SEL}">{a}</span><span style="color: {MUTED}">{b}</span></div>' for a, b in grammar)
    right = (panel_head("Select by query", chip("⌘F", MUTED, BG2, True))
             + section("Saved queries", pt)
             + section("Grammar", f'<div>{gr}</div>')
             + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Counts are from the drawn chain ({len(atoms)} atoms). Same grammar in the command palette, macros and Python: <span style="font-family: {MONO}">caps.select("…")</span>.</div>'
             + footer(btn("Save as group", ic="tag"), btn("Select", True)))
    tb = bar(tbtn("cursor", "Select", True), tbtn("lasso", "Lasso"), tbtn("search", "Query", True, "Query"), sep(), tbtn("ruler", "Measure"), spacer())
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">{q}'
              f'<div style="position: relative; flex-grow: 1">{view}{hud(chip("PS · atactic · DP 10", TEXT, "#16191Ccc"), chip("selection halo", SEL, "#16191Ccc"))}</div></div>')
    return studio_shell("PS_atactic_DP10.caps", "CAPS Studio — select by query", tb, center, right,
                        (f"<span>{len(atoms)} atoms · {len(arom)} selected</span>", "<span>query parsed · 0.4 ms</span>"))


# ---------------------------------------------------------------- Figure composer
def figure():
    W_in, dpi = 7.0, 600
    Hin = 4.2
    px = (round(W_in * dpi), round(Hin * dpi))
    scale = 118  # screen px per inch on this board
    fw, fh = round(W_in * scale), round(Hin * scale)
    atoms, _, _ = polystyrene(6, "isotactic", seed=2)
    pw, ph = fw // 2 - 18, fh // 2 - 18
    sa = Scene("fg", pw, ph - 6, yaw=-0.1, pitch=0.45, fog=0.0, outline=True, bg="#FFFFFF")
    sa.add_atoms(atoms)
    c1, c2, c3 = 355.03, -68.19, 791.32
    U = lambda d: (c1 * (1 + math.cos(math.radians(d))) + c2 * (1 - math.cos(2 * math.radians(d))) + c3 * (1 + math.cos(3 * math.radians(d)))) * KCAL_PER_K
    curve = [(d, U(d)) for d in range(-180, 181, 3)]
    ink, grid = "#1B1E21", "#C9CDD1"

    def lightplot(series, xr, yr, xt, yt, xl, yl):
        s = plot(pw, ph - 6, series, xr, yr, xt, yt, xl, yl, pad=(34, 8, 8, 26))
        for a, b in ((TEXT, ink), (MUTED, "#4A5056"), (DIM, "#6A7178"), (BG3, grid), (LINE, grid), (BG2, "#E6E8EA")):
            s = s.replace(a, b)
        return s
    pb = lightplot([(curve, "#B26B00", 1.6, None)], (-180, 180), (0, 7), [-180, 0, 180], [0, 2, 4, 6], "φ (°)", "E (kcal/mol)")
    rs = [(r / 10, 0 if r < 22 else 1 + 1.4 * math.exp(-((r / 10 - 2.6) ** 2) / 0.03) + 0.3 * math.exp(-((r / 10 - 4.0) ** 2) / 0.2)) for r in range(0, 120)]
    pc = lightplot([(rs, "#2271DB", 1.6, None)], (0, 12), (0, 3), [0, 4, 8, 12], [0, 1, 2, 3], "r (Å)", "g(r)")
    bars = [("PS", "#F0A83C"), ("PE", "#2271DB"), ("PMMA", "#DE775D")]
    pd = lightplot([([(0.5 + i, 0), (0.5 + i, 0)], c, 1, None) for i, (_, c) in enumerate(bars)], (0, 3.5), (0, 1.2), [], [0, 0.5, 1.0], "", "ρ (g/cm³)")
    bw = (pw - 34 - 8) / 3.5 * 0.6
    bsvg = "".join(f'<rect x="{34 + (0.5 + i) / 3.5 * (pw - 42) - bw / 2:.1f}" y="{8 + (ph - 40) * 0.25:.1f}" width="{bw:.1f}" height="{(ph - 40) * 0.75 - 0:.1f}" fill="{c}" opacity="0.35"></rect>'
                   f'<text x="{34 + (0.5 + i) / 3.5 * (pw - 42):.1f}" y="{ph - 14}" text-anchor="middle" font-size="10" fill="#4A5056">{n}</text>' for i, (n, c) in enumerate(bars))
    pd = pd.replace("</svg>", bsvg + '<text x="50%" y="40%" text-anchor="middle" font-size="11" fill="#6A7178">[result]</text></svg>')
    panels = [("a", sa.svg()), ("b", pb), ("c", pc), ("d", pd)]
    pn = "".join(f'<div style="position: relative; width: {pw}px; height: {ph}px; {"outline: 2px solid " + ACC + "; outline-offset: 2px;" if l == "b" else ""}">'
                 f'<span style="position: absolute; left: 2px; top: 0; font-family: Arial, sans-serif; font-weight: 700; font-size: 13px; color: {ink}; z-index: 1">{l}</span>{s}</div>' for l, s in panels)
    sheet = (f'<div style="width: {fw}px; height: {fh}px; background: #FFFFFF; box-shadow: 0 8px 30px #0008; display: grid; grid-template-columns: repeat(2, {pw}px); gap: 12px; padding: 12px; box-sizing: border-box">{pn}</div>')
    ruler = f'<div style="width: {fw}px; display: flex; justify-content: space-between; font-family: {MONO}; font-size: 10.5px; color: {DIM}">' + "".join(f"<span>{i} in</span>" for i in range(8)) + '</div>'
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 8px; background: {BG0}">{ruler}{sheet}'
              f'<span style="font-size: 11.5px; color: {DIM}">Panels b–d plot curves already in the project; d shows [result] until the densities are computed.</span></div>')
    right = (panel_head("Figure", chip("4 panels", MUTED, BG2))
             + section("Page", col(select("Preset", "Double column · 7.0 in"), row(field("Width", f"{W_in}", "in"), field("Height", f"{Hin}", "in"), field("DPI", str(dpi)), gap=8),
                                   kv("Output size", f"{px[0]} × {px[1]} px"), gap=8))
             + section("Panel b · torsion", col(row(select("Source", "butane_scan · potential"), gap=8), row(field("Line", "1.0", "pt"), select("Colour", "amber · v2"), gap=8), gap=8))
             + section("Text", col(row(select("Font", "Arial / Helvetica"), field("Size", "7", "pt"), gap=8), check("Panel labels bold, lower case"), check("Same axis font in every panel"), gap=8))
             + section("Checks", col(row(icon("check", 14, OK), '<span style="font-size: 12px">Lines ≥ 0.5 pt at print size</span>', gap=8),
                                     row(icon("check", 14, OK), '<span style="font-size: 12px">Palette v2 · colour-vision safe</span>', gap=8),
                                     row(icon("alert", 14, ACC), '<span style="font-size: 12px">Panel d has no data yet</span>', gap=8), gap=6))
             + footer(btn("SVG"), btn("TIFF"), btn("Export PDF", True, "download")))
    tb = bar(tbtn("cursor", "Select", True), tbtn("move", "Move"), sep(), tbtn("layers", "Panels", False, "Layout 2 × 2"), tbtn("tag", "Labels", False, "a b c d"), spacer(), chip("100 % · 118 px per inch", MUTED, BG2, True))
    return studio_shell("PS_paper_fig2.caps", "CAPS Studio — figure composer", tb, center, right,
                        ("<span>Figure 2 · 4 panels</span>", f"<span>{px[0]} × {px[1]} px at {dpi} dpi</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage7/project", exist_ok=True)
    for name, fn in (("RecipeEditor", recipe), ("SmartSelect", smart_select), ("FigureComposer", figure)):
        open(f"stage7/project/{name}.dc.html", "w").write(fn())
        print(name, "ok")
    print("larsen total", sum(s[3] for s in LARSEN))
    a, b, st = polystyrene(10, "atactic", seed=8)
    print(len(a), sum(1 for x in a if x.get("ar")))
