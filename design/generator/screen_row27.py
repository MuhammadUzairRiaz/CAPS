"""Row 27 — figure export backgrounds: dark, white, transparent."""
import math
from lib import *
from screen_row20 import draw, VIRIDIS, ramp
from screen_row21 import cell_recs, inspector, frame, CW, VH, BOTTOM, L, MASS
from screen_row22 import mol_shape, chain_col, H3, FULL

INK = {"Dark": ("#0F1113", TEXT, MUTED), "White": ("#FFFFFF", "#141413", "#5A6168"), "Transparent": (None, "#141413", "#5A6168")}


def checker(w, h, s=12):
    return (f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" aria-hidden="true" style="position: absolute; left: 0; top: 0"><defs><pattern id="chk" width="{2 * s}" height="{2 * s}" patternUnits="userSpaceOnUse">'
            f'<rect width="{2 * s}" height="{2 * s}" fill="#FFFFFF"></rect><rect width="{s}" height="{s}" fill="#D9DCDF"></rect><rect x="{s}" y="{s}" width="{s}" height="{s}" fill="#D9DCDF"></rect></pattern></defs>'
            f'<rect width="{w}" height="{h}" fill="url(#chk)"></rect></svg>')


def figure_background():
    cell, recs = cell_recs()
    sh = mol_shape(recs)
    for r in recs:
        r["dcom"] = dist(r["p"], sh[r["mol"]]["com"])
    hi = math.ceil(max(r["dcom"] for r in recs))
    pw, ph = (CW - 6 - 2 * 10 - 24) // 3, 400
    tiles = []
    for mode in ("Dark", "White", "Transparent"):
        bg, ink, sub = INK[mode]
        sc = Scene("fb" + mode[0], pw, ph - 34, yaw=0.55, pitch=0.4, persp=0.0, fog=0.25 if bg else 0.0, atom_k=1.2, center=(L / 2, L / 2, L / 2), bg=bg or "#FFFFFF")
        draw(sc, cell, recs, chain_col)
        sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), sub, 1.0, None, 0.8)
        sc._prep()
        bar = 10 * sc.scale
        ov = (f'<text x="12" y="20" font-size="12" font-family="IBM Plex Sans" font-weight="600" fill="{ink}">PS melt · {len(recs)} atoms</text>'
              f'<line x1="12" y1="{ph - 34 - 16}" x2="{12 + bar:.1f}" y2="{ph - 34 - 16}" stroke="{ink}" stroke-width="3"></line>'
              f'<text x="{12 + bar / 2:.1f}" y="{ph - 34 - 22}" text-anchor="middle" font-size="10.5" font-family="IBM Plex Mono" fill="{ink}">10 Å</text>')
        back = checker(pw, ph - 34) if bg is None else ""
        sel = mode == "Transparent"
        tiles.append(f'<div style="width: {pw}px; flex-shrink: 0; display: flex; flex-direction: column; border-radius: 8px; overflow: hidden; outline: {"2px solid " + ACC if sel else "1px solid " + LINE}; outline-offset: -1px">'
                     f'<div style="position: relative; height: {ph - 34}px; background: {bg or "transparent"}">{back}<div style="position: relative">{sc.svg(overlay=ov)}</div></div>'
                     f'<div style="height: 34px; display: flex; align-items: center; justify-content: space-between; padding: 0 12px; background: {BG1}; font-size: 12px"><span style="font-weight: 600{"; color: " + ACC if sel else ""}">{mode}</span>'
                     f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{ {"Dark": "screen, talks", "White": "papers, print", "Transparent": "slides, posters, overlays"}[mode] }</span></div></div>')
    row3 = f'<div style="display: flex; gap: 12px">{"".join(tiles)}</div>'
    sizes = [("Journal · single column", 85, 600), ("Journal · double column", 170, 600), ("Poster panel", 250, 300), ("Slide 16:9", None, None)]
    rows = []
    for name, mm, dpi in sizes:
        if mm:
            px = round(mm / 25.4 * dpi)
            rows.append([name, f"{mm} mm", str(dpi), f"{px} × {round(px * 9 / 16)}"])
        else:
            rows.append([name, "—", "—", "1920 × 1080"])
    st = table(["Preset", "Width", "dpi", "Pixels (16:9)"], rows, ["40%", "18%", "14%", "28%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=12, rowh=26, hl={0})
    rules = [("Transparent", "PNG and TIFF keep an alpha channel; SVG and PDF are written with no background shape"),
             ("Depth fog", "fades toward the background colour; with no background it is replaced by darkening far atoms"),
             ("Labels, scale bar, legend", "switch to dark ink on White and Transparent so they stay readable"),
             ("Edges", "antialiased against alpha, so no dark fringe when placed on a coloured slide")]
    rl = "".join(f'<div style="display: flex; gap: 10px; padding: 6px 0; border-bottom: 1px solid {BG2}; font-size: 12px"><span style="width: 150px; flex-shrink: 0; font-weight: 600">{a}</span><span style="color: {MUTED}; line-height: 1.5">{b}</span></div>' for a, b in rules)
    bottom = (f'<div style="flex-grow: 1; min-height: 0; display: flex; gap: 18px; padding: 12px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="width: 470px; flex-shrink: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Size presets</h3>{st}</div>'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 4px"><h3 style="{H3}">What changes with the background</h3>{rl}</div></div>')
    content = f'<div style="flex-shrink: 0; padding: 14px 12px; background: {BG0}">{row3}</div>' + bottom
    right = inspector("Export figure", section("Background", col(seg(["Dark", "White", "Transparent"], "Transparent", full=True), select("Custom colour", "none"), toggle("Remember per project", True), gap=8))
                      + section("File", col(row(select("Format", "PNG · RGBA"), select("Bit depth", "8-bit"), gap=8), row(field("Width", "85", "mm"), field("dpi", "600"), gap=8), kv("Pixels", "2008 × 1130"), gap=8))
                      + section("Content", col(toggle("Cell edges", True), toggle("Scale bar", True), toggle("Title", True), toggle("Colour legend", False), gap=8))
                      + f'<div style="padding: 0 18px 12px">{btn("Export figure", True, "download")}</div>')
    return frame(content, right, "CAPS — export figure background", "Export › Figure",
                 ("<span>same frame, three backgrounds · colour by molecule · orthographic, scale bar true 10 Å</span>", "<span>checkerboard marks transparent pixels</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage27/project", exist_ok=True)
    h = figure_background(); open("stage27/project/FigureBackground.dc.html", "w").write(h); print("FigureBackground", len(h))
    from screen_row22 import render_overlays
    h = render_overlays(); open("stage27/project/RenderOverlays.dc.html", "w").write(h); print("RenderOverlays", len(h))
