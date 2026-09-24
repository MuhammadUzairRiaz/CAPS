import re, random, itertools
from lib import *
from mols import *
from screen_app import card
from screen_row14 import simulate, dE, lab
from palette_opt import cr

KINDS = ["Protanopia", "Deuteranopia", "Tritanopia"]
PAL = [
    ("Elements", ["C", "H", "O", "N", "Si", "S", "Na", "Cl"],
     ["#8E959C", "#E9ECEF", "#E5534B", "#4C7BD9", "#D6A45E", "#E3C74A", "#9B7BD6", "#57B26A"],
     ["#8B969E", "#E9ECEF", "#E35049", "#2271DB", "#D6A357", "#E4C84C", "#9B7AD5", "#58B573"]),
    ("Chains", ["1", "2", "3", "4", "5", "6", "7", "8"],
     ["#F0A83C", "#6CC4D8", "#E07A5F", "#9B7BD6", "#7CC784", "#D6A45E", "#E9ECEF", "#4C7BD9"],
     ["#F0A83C", "#6CC4D8", "#DE775D", "#9B7AD5", "#7DC884", "#D6AC5C", "#E9ECEF", "#2271DB"]),
    ("Status", ["ok", "run", "error", "select"],
     [OK, ACC, ERR, SEL],
     ["#83CC91", ACC, ERR, SEL]),
]


def worst(names, cols, kind):
    cs = [c if kind == "Normal" else simulate(c, kind) for c in cols]
    return min((dE(cs[i], cs[j]), names[i], names[j]) for i, j in itertools.combinations(range(len(cs)), 2))


def sim_svg(svg, kind):
    return re.sub(r"#[0-9A-Fa-f]{6}\b", lambda m: simulate(m.group(0).upper(), kind), svg)


def swatch_row(label, names, cols, ref=None):
    out = []
    for n, c in zip(names, cols):
        txt = "#0F1113" if lab(c)[0] > 55 else "#FFFFFF"
        d = dE(c, ref[names.index(n)]) if ref else 0
        tag = f'<span style="font-family: {MONO}; font-size: 10px; color: {ACC if d >= 0.5 else BG3}">{"Δ %.1f" % d if d >= 0.5 else "·"}</span>' if ref else ""
        out.append(f'<div style="display: flex; flex-direction: column; align-items: center; gap: 3px"><span style="width: 46px; height: 30px; border-radius: 5px; background: {c}; display: flex; align-items: center; justify-content: center; font-family: {MONO}; font-size: 11px; color: {txt}">{n}</span>{tag}</div>')
    return (f'<div style="display: flex; align-items: flex-start; gap: 10px"><span style="width: 26px; padding-top: 7px; flex-shrink: 0; font-family: {MONO}; font-size: 11px; color: {DIM}">{label}</span>'
            f'<div style="display: flex; gap: 5px">{"".join(out)}</div></div>')


def chain_scene(sid, cols, w, h):
    Lb = 30.0
    sc = Scene(sid, w, h, yaw=0.6, pitch=0.42, persp=0.3, fog=0.45, outline=True)
    for k in range(8):
        rnd = random.Random(300 + k)
        st = (rnd.uniform(4, Lb - 4), rnd.uniform(4, Lb - 4), rnd.uniform(4, Lb - 4))
        pts = random_chain(st, 46, 400 + k, bond=1.54, box=((0.8, 0.8, 0.8), (Lb - 0.8, Lb - 0.8, Lb - 0.8)), persistence=0.6)
        sc.add_tube(pts, cols[k], 0.7)
    sc.add_box((0, 0, 0), (Lb, 0, 0), (0, Lb, 0), (0, 0, Lb), MUTED, 1.0, None, 0.7)
    return sc.svg()


def ion_scene(sid, cols, w, h):
    em = dict(zip(PAL[0][1], cols))
    rnd = random.Random(7)
    L = 16.0
    atoms = []
    for e, n, r in (("N", 9, 0.62), ("Na", 9, 0.72), ("Cl", 7, 0.8), ("O", 6, 0.55)):
        for _ in range(n):
            for _t in range(200):
                p = (rnd.uniform(1, L - 1), rnd.uniform(1, L - 1), rnd.uniform(1, L - 1))
                if all(sum((p[k] - a["p"][k]) ** 2 for k in range(3)) > 9 for a in atoms):
                    break
            atoms.append({"e": e, "p": p, "r": r, "c": em[e]})
    sc = Scene(sid, w, h, yaw=0.5, pitch=0.38, persp=0.3, fog=0.35, outline=True)
    sc.add_atoms(atoms, [])
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    return sc.svg()


def palette_v2():
    cards = []
    summary = []
    for name, names, v1, v2 in PAL:
        w1 = min(worst(names, v1, k) for k in ["Normal"] + KINDS)
        w2 = min(worst(names, v2, k) for k in ["Normal"] + KINDS)
        changed = sum(1 for a, b in zip(v1, v2) if a != b)
        cards.append(card(name, col(swatch_row("v1", names, v1), swatch_row("v2", names, v2, v1), gap=8),
                          chip(f"min ΔE {w1[0]:.1f} → {w2[0]:.1f}", OK if w2[0] >= 12 else ERR, BG2, True), 14))
        for k in KINDS:
            a, b = worst(names, v1, k), worst(names, v2, k)
            summary.append([name, k, f"{a[1]}↔{a[2]}  {a[0]:.1f}", f"{b[1]}↔{b[2]}  {b[0]:.1f}"])
    tbl = table(["Palette", "Vision", "Worst · v1", "Worst · v2"], [r for r in summary if r[0] != "Status"], ["20%", "28%", "26%", "26%"], mono_cols=(2, 3), fs=11.5, rowh=24)
    def panels(make, v1, v2, tag):
        out = []
        for lbl, cols in (("v1", v1), ("v2", v2)):
            s = make(f"{tag}{lbl}", cols)
            for kind, svg in (("normal", s), ("protan", sim_svg(re.sub(r'(id="|url\(#|href="#)', r"\1p", s), "Protanopia"))):
                out.append(f'<div style="display: flex; flex-direction: column; gap: 5px"><div style="background: {BG0}; border: 1px solid {LINE}; border-radius: 8px; overflow: hidden; height: 150px">{svg}</div>'
                           f'<span style="font-size: 11.5px; color: {MUTED}"><span style="font-family: {MONO}; color: {TEXT}">{lbl}</span> · {kind}</span></div>')
        return f'<div style="display: grid; grid-template-columns: repeat(4, 200px); gap: 12px">{"".join(out)}</div>'
    chains = panels(lambda sid, c: chain_scene(sid, c, 200, 150), PAL[1][2], PAL[1][3], "pc")
    ions = panels(lambda sid, c: ion_scene(sid, c, 200, 150), PAL[0][2], PAL[0][3], "pi")
    minc = min(cr(c, BG0) for _, _, _, v2 in PAL for c in v2)
    notes = [("Method", "Machado 2009 matrices at severity 1.0, CIE76 ΔE*ab, target ≥ 12 in all four visions. Hues kept; mostly lightness and chroma moved."),
             ("Biggest move", "Nitrogen and chain 8 blue: #4C7BD9 → #2271DB (Δ 8.9), darker and more saturated so it separates from Na purple."),
             ("Shared hues", "Chain 4 reuses Na purple and chain 8 reuses N blue. Chain 6 cannot reuse Si tan: it collides with chain 1 under tritanopia (ΔE 10.1)."),
             ("Contrast", f"Every v2 colour keeps ≥ 3:1 against the viewport background {BG0} (lowest {minc:.1f}:1, WCAG 1.4.11)."),
             ("Not applied", "Other boards still use v1. Switching is one change in the element and chain token maps.")]
    nt = "".join(f'<div style="display: flex; gap: 10px; font-size: 12px; line-height: 1.5"><span style="width: 88px; flex-shrink: 0; color: {DIM}">{a}</span><span>{b}</span></div>' for a, b in notes)
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 16px; padding: 26px 36px; min-height: 0">'
            + f'<div style="display: flex; align-items: flex-end; gap: 16px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">Palette v2</h1>'
              f'<span style="font-size: 13px; color: {MUTED}">Retuned so no pair drops below ΔE 12 for protan, deutan or tritan vision</span><div style="flex-grow: 1"></div>{chip("Proposal", ACC, BG2)}</div>'
            + '<div style="display: flex; gap: 18px; min-height: 0">'
            + f'<div style="width: 470px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">{"".join(cards)}{card("Worst pair per vision", tbl, "", 12)}</div>'
            + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">{card("Ions and heteroatoms · N, Na, Cl, O", ions, chip("N ↔ Na", MUTED, BG2), 14)}'
              f'{card("Eight chains", chains, chip("chains 4 ↔ 8", MUTED, BG2), 14)}{card("Notes", col(nt, gap=7), "", 14)}</div>'
            + '</div></div>')
    return page("CAPS — palette v2", body)


if __name__ == "__main__":
    import os
    os.makedirs("stage4/project", exist_ok=True)
    open("stage4/project/PaletteV2.dc.html", "w").write(palette_v2())
    for _, n, v1, v2 in PAL:
        print(n, [round(min(worst(n, v, k)[0] for k in ["Normal"] + KINDS), 1) for v in (v1, v2)])
