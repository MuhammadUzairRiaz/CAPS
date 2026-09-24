import math
from lib import *
from screen_studio2 import studio_shell, bar, spacer
from screen_builders import footer, hud
from display import ps_cell, render_cell, count, STYLES

VW = 1440 - 72 - 380
VH = 900 - 44 - 48 - 26


def display_styles():
    L = 33.0
    cell = ps_cell(10, 8, L, seed=21)
    n, h = count(cell)
    heavy = n - h
    nbb = sum(len(b) for _, _, b in cell)
    pw, ph = (VW - 36) // 2, (VH - 36) // 2 - 26
    def panel(style, sub_):
        sc = Scene(f"ds{STYLES.index(style)}", pw, ph, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.3)
        render_cell(sc, cell, style, colours=None if style != "Backbone" else ["#F0A83C", "#6CC4D8", "#DE775D", "#9B7AD5", "#7DC884", "#D6AC5C", "#E9ECEF", "#2271DB"])
        sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, "5 4", 0.6)
        on = style == "All atoms"
        return (f'<div style="display: flex; flex-direction: column; gap: 6px"><div style="display: flex; align-items: baseline; gap: 10px"><span style="font-size: 13px; font-weight: 600">{style}</span><span style="font-size: 11.5px; color: {DIM}">{sub_}</span></div>'
                f'<div style="background: {BG0}; border: 1px solid {ACC if on else LINE}; border-radius: 8px; overflow: hidden; height: {ph}px">{sc.svg()}</div></div>')
    rules = [("≤ 20 k", "All atoms", "H shown"), ("20 k – 200 k", "No H", "H return on zoom"),
             ("200 k – 2 M", "Backbone", "atoms in a 30 Å lens"), ("> 2 M", "Backbone + LOD", "impostors")]
    rt = table(["Atoms in view", "Default", "Detail"], [[a, b, c] for a, b, c in rules], ["28%", "34%", "38%"], mono_cols=(0,), fs=11.5, rowh=34)
    auto = (f'<div style="display: flex; flex-direction: column; gap: 6px"><div style="display: flex; align-items: baseline; gap: 10px"><span style="font-size: 13px; font-weight: 600">Automatic by size</span><span style="font-size: 11.5px; color: {DIM}">thresholds are settings</span></div>'
            f'<div style="background: {BG1}; border: 1px solid {LINE}; border-radius: 8px; padding: 12px; height: {ph}px; box-sizing: border-box; display: flex; flex-direction: column; gap: 10px">{rt}'
            f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">The status bar always says which style is on and how many atoms the model holds, so a simplified view is never mistaken for a simplified model.</div></div></div>')
    grid = (f'<div style="display: grid; grid-template-columns: repeat(2, {pw}px); gap: 12px; padding: 12px">'
            + panel("All atoms", f"{n:,} atoms · {h} H".replace(",", " ")) + panel("No H", f"{heavy} heavy atoms drawn · {h} H hidden")
            + panel("Backbone", f"{len(cell)} tubes through {nbb} backbone C") + auto + '</div>')
    right = (panel_head("Display", chip("view only", MUTED, BG2))
             + section("Style", col(seg(STYLES, "All atoms", full=True), row(select("Colour by", "element"), select("Atoms", "ball-and-stick"), gap=8), gap=8))
             + section("Hydrogens", col(toggle("Show hydrogens", True), toggle("Show polar H only", False), toggle("Selections always show all atoms", True), gap=8))
             + section("What never changes", col(*[f'<div style="display: flex; gap: 8px; font-size: 12.5px; line-height: 1.5">{icon("check", 15, OK)}<span>{t}</span></div>' for t in (
                 "The model is all-atom with every H: typing, charges, energies, density and exports use all of them.",
                 "Hiding H or drawing tubes changes the picture only.",
                 "Measurements snap to real atoms, also in Backbone style.")], gap=8))
             + footer(btn("Reset"), btn("Save as default", True)))
    tb = bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("eye", "Display", True, "Display styles"), spacer(), chip(f"model: {n:,} atoms".replace(",", " "), MUTED, BG2, True))
    center = f'<div style="flex-grow: 1; min-width: 0; background: {BG0}; overflow: hidden">{grid}</div>'
    return studio_shell("PS_atactic_cell.caps", "CAPS Studio — display styles", tb, center, right,
                        (f"<span>{len(cell)} chains · DP 8 · {n:,} atoms · {h} H</span>".replace(",", " "), "<span>Display: All atoms</span>"))


def paper_of(fn_html, names=("Start", "Main", "PolymerBuilder", "Grow", "Jobs", "Analyze")):
    """Paper (light) variant, same recipe as the published Paper boards."""
    import re
    from screen_row7 import LIGHT
    M = dict(LIGHT); M["#0B0D0F"] = "#FFFFFF"
    h = re.sub(r"#[0-9A-Fa-f]{8}|#[0-9A-Fa-f]{6}", lambda m: M.get(m.group(0).upper(), m.group(0)), fn_html)
    for m in names:
        h = h.replace(f'href="{m}.dc.html"', f'href="Paper{m}.dc.html"')
    return h.replace("<title>", "<title>Paper · ", 1)


if __name__ == "__main__":
    import os, sys
    os.makedirs("stage14/project", exist_ok=True)
    import screen_app, screen_row18, screen_row16b, screen_row15c, screen_row17
    out = {"DisplayStyles": display_styles, "Grow": screen_app.grow, "GrowAllAtom": lambda: screen_app.grow("All atoms"),
           "Polydispersity": screen_row18.polydispersity, "PeriodicBox": screen_row16b.periodic, "Interface": screen_row15c.interface, "ChainStats": screen_row17.chain_stats}
    for n, f in out.items():
        h = f()
        if n == "GrowAllAtom":
            h = h.replace("<title>CAPS Grow — amorphous cell</title>", "<title>CAPS Grow — all-atom view</title>")
        open(f"stage14/project/{n}.dc.html", "w").write(h)
        print(n, len(h))
    # Paper Grow: Scene default background -> paper
    d = list(Scene.__init__.__defaults__); d[-1] = "#F4F2EE"; Scene.__init__.__defaults__ = tuple(d)
    h = paper_of(screen_app.grow())
    open("stage14/project/PaperGrow.dc.html", "w").write(h)
    print("PaperGrow", len(h), "leftover dark", {c: h.upper().count(c) for c in ["#0F1113", "#16191C", "#262B30"] if h.upper().count(c)})
