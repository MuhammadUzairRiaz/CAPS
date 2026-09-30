"""Row 29-30 — editing, selection & toolbars (A1–A15): each new tool drawn over the real CAPS Studio as it is today (a
screenshot of the Studio with an h-BN (001) slab under a cis-1,4-polyisoprene film, 1,426 atoms), so the design keeps
the app's layout and only adds to it. Toolbars can be collapsed, dragged off to float, docked on any edge and reset."""
import os
from lib import *

SHOT = "/_blob/2b05b6ef69520f6fa412a6de8ee11ce7"   # CAPS Studio, iface.data, 1440 × 900
H3 = f"margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}"
SHADOW = "box-shadow: 0 18px 48px #000000a0, 0 2px 6px #00000080"


def at(x, y, inner, w=None, h=None, extra=""):
    size = (f"width: {w}px; " if w else "") + (f"height: {h}px; " if h else "")
    return f'<div style="position: absolute; left: {x}px; top: {y}px; {size}{extra}">{inner}</div>'


def board(title, *layers, dim=None):
    """The Studio screenshot with layers over it; dim: [(x, y, w, h)] areas left bright while the rest is shaded."""
    shade = ""
    if dim is not None:
        shade = at(0, 0, "", 1440, 900, "background: #0B0D0F8c")
        for (x, y, w, h) in dim:
            shade += at(x, y, "", w, h, f"background: transparent; outline: 1px solid {ACC}66; box-shadow: 0 0 0 9999px transparent")
    body = (f'<div style="position: relative; width: 1440px; height: 900px; overflow: hidden">'
            f'<img src="{SHOT}" alt="CAPS Studio today: project tree, toolbars, 3D view of the h-BN slab under a polyisoprene film, analysis dock and inspector" '
            f'style="position: absolute; left: 0; top: 0; width: 1440px; height: 900px">'
            + shade + "".join(layers) + '</div>')
    return page(title, body)


def cover(x, y, w, h, color=BG1):
    """Paints over a part of the screenshot that the new design replaces."""
    return at(x, y, "", w, h, f"background: {color}")


def panel(inner, w, pad=14, gap=10, extra=""):
    return (f'<div style="width: {w}px; box-sizing: border-box; display: flex; flex-direction: column; gap: {gap}px; padding: {pad}px; background: {BG1}; '
            f'border: 1px solid {LINE}; border-radius: 10px; {SHADOW}; {extra}">{inner}</div>')


def callout(n, text, w=260):
    return (f'<div style="width: {w}px; display: flex; gap: 10px; align-items: flex-start; padding: 10px 12px; background: #1A1408; border: 1px solid {ACC}; '
            f'border-radius: 8px; {SHADOW}"><span style="width: 22px; height: 22px; flex-shrink: 0; border-radius: 50%; background: {ACC}; color: {ACC_INK}; '
            f'font-size: 12px; font-weight: 700; display: flex; align-items: center; justify-content: center">{n}</span>'
            f'<span style="font-size: 12px; line-height: 1.45; color: {TEXT}">{text}</span></div>')


def menu(items, w=260):
    """items: (label, shortcut, kind) kind: '' | 'sep' | 'head' | 'on' | 'off' | 'sub' | 'hl' | 'dis'."""
    out = []
    for label, sc, kind in items:
        if kind == "sep":
            out.append(f'<div style="height: 1px; background: {LINE}; margin: 4px 0"></div>'); continue
        if kind == "head":
            out.append(f'<div style="padding: 6px 12px 3px; font-size: 10.5px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {DIM}">{label}</div>'); continue
        mark = ""
        if kind in ("on", "off"):
            mark = icon("check", 14, ACC) if kind == "on" else '<span style="width: 14px; display: inline-block"></span>'
        bg = f"background: {BG3}; " if kind == "hl" else ""
        col = DIM if kind == "dis" else TEXT
        sub = icon("chevr", 14, MUTED) if kind == "sub" else ""
        out.append(f'<button style="width: 100%; height: 28px; display: flex; align-items: center; gap: 8px; padding: 0 12px; {bg}border: 0; border-radius: 5px; '
                   f'background-color: {BG3 if kind == "hl" else "transparent"}; cursor: pointer; text-align: left"><span style="width: 14px; display: flex">{mark}</span>'
                   f'<span style="flex-grow: 1; font-size: 12.5px; color: {col}">{label}</span>'
                   f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{sc}</span>{sub}</button>')
    return (f'<div role="menu" style="width: {w}px; box-sizing: border-box; padding: 5px; background: {BG2}; border: 1px solid {LINE}; border-radius: 8px; {SHADOW}">'
            + "".join(out) + '</div>')


def grip(vertical=False):
    dots = "".join(f'<span style="width: 3px; height: 3px; border-radius: 50%; background: {DIM}"></span>' for _ in range(6))
    lay = "grid-template-columns: repeat(3, 3px)" if vertical else "grid-template-columns: repeat(2, 3px)"
    return (f'<button aria-label="Drag the toolbar" title="Drag to move · double-click to float" style="width: 16px; height: 28px; padding: 0; background: transparent; '
            f'border: 0; cursor: grab; display: grid; {lay}; gap: 3px; align-content: center; justify-content: center">{dots}</button>')


def mini(label, w=None, active=False, mono=False):
    ff = f"font-family: {MONO}; " if mono else ""
    wd = f"width: {w}px; " if w else "padding: 0 7px; "
    return (f'<button style="{wd}{ff}height: 26px; display: inline-flex; align-items: center; justify-content: center; font-size: 12px; font-weight: 600; '
            f'color: {ACC if active else TEXT}; background: {BG3 if active else "transparent"}; border: 1px solid {LINE if active else "transparent"}; border-radius: 5px; cursor: pointer">{label}</button>')


def toolbar_modify(collapsed=False):
    """The MODIFY row as a movable toolbar: grip, name, tools, collapse chevron."""
    if collapsed:
        return (f'<div style="height: 30px; display: inline-flex; align-items: center; gap: 6px; padding: 0 6px 0 2px; background: {BG2}; border: 1px solid {LINE}; border-radius: 7px">'
                f'{grip()}<span style="font-size: 11px; font-weight: 600; letter-spacing: 0.08em; color: {MUTED}">MODIFY</span>'
                f'<button aria-label="Expand the Modify toolbar" style="width: 22px; height: 22px; padding: 0; border: 0; background: transparent; cursor: pointer; display: flex; align-items: center; justify-content: center">{icon("chevr", 14, MUTED)}</button></div>')
    tools = "".join(mini(e, 26) for e in ["H", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br"])
    return (f'<div style="height: 32px; display: flex; align-items: center; gap: 4px; padding: 0 4px 0 2px">{grip()}'
            f'<span style="font-size: 11px; font-weight: 600; letter-spacing: 0.08em; color: {MUTED}; margin-right: 4px">MODIFY</span>{tools}'
            f'{sep()}{mini("—", 28)}{mini("=", 28)}{mini("≡", 28)}{sep()}{mini("sp³ ⌄", 50)}'
            f'<div style="flex-grow: 1"></div><button aria-label="Collapse the Modify toolbar" title="Collapse" style="width: 24px; height: 24px; padding: 0; border: 0; background: transparent; '
            f'cursor: pointer; display: flex; align-items: center; justify-content: center">{icon("chev", 14, MUTED)}</button></div>')


def scope_chip(text):
    return (f'<span style="display: inline-flex; align-items: center; gap: 6px; height: 24px; padding: 0 10px; border-radius: 12px; background: #10262B; '
            f'border: 1px solid {SEL}; color: {SEL}; font-size: 11.5px; white-space: nowrap">{icon("lasso", 13, SEL)}{text}</span>')


def lasso_blue():
    """The lasso ring around the blue chain (molecule 4) in the screenshot's view."""
    return at(612, 452, f'<svg width="200" height="130" viewBox="0 0 200 130" aria-hidden="true"><path d="M30 18 C80 -2 160 6 186 40 C198 70 170 112 120 122 C70 128 18 108 8 74 C2 48 10 28 30 18 Z" '
                        f'fill="{SEL}1f" stroke="{SEL}" stroke-width="1.6" stroke-dasharray="6 4"></path></svg>')




# ================================================================== CAPS's own interaction vocabulary
def shelf(tools, name, glyph="layers", vertical=False, active=None):
    """A tool shelf: an amber tab (drag to tear off / snap) then its tools, in a capsule."""
    tab = (f'<button aria-label="Move the {name} shelf" title="{name} · drag to move, click to fold" style="{"width: 34px; height: 22px" if vertical else "width: 22px; height: 34px"}; '
           f'padding: 0; border: 0; border-radius: {"14px 14px 4px 4px" if vertical else "14px 4px 4px 14px"}; background: {ACC}; cursor: grab; display: flex; '
           f'align-items: center; justify-content: center">{icon(glyph, 13, ACC_INK, 2)}</button>')
    items = "".join(tbtn(i, n, n == active) for i, n in tools)
    lay = "flex-direction: column" if vertical else "flex-direction: row"
    return (f'<div role="toolbar" aria-label="{name}" style="display: inline-flex; {lay}; align-items: center; gap: 2px; padding: 2px; background: {BG1}; '
            f'border: 1px solid {LINE}; border-radius: 16px; {SHADOW}">{tab}{items}</div>')


def puck(glyph, name, n=None):
    """A folded shelf: a round puck parked on an edge; hover opens it in place."""
    badge = (f'<span style="position: absolute; right: -4px; top: -4px; min-width: 16px; height: 16px; padding: 0 4px; box-sizing: border-box; border-radius: 8px; '
             f'background: {BG3}; border: 1px solid {LINE}; font-size: 10px; color: {MUTED}; display: flex; align-items: center; justify-content: center">{n}</span>') if n else ""
    return (f'<button aria-label="Open the {name} shelf" title="{name}" style="position: relative; width: 34px; height: 34px; padding: 0; border-radius: 50%; '
            f'background: {BG1}; border: 1.5px solid {ACC}; cursor: pointer; display: flex; align-items: center; justify-content: center; {SHADOW}">{icon(glyph, 15, ACC)}{badge}</button>')


def selbar(count, what, how, actions=("eye", "tag", "copy", "atom", "wand", "xcircle", "dots")):
    names = {"eye": "Hide", "tag": "Tag", "copy": "Copy", "atom": "Probe", "wand": "Clean up", "xcircle": "Delete", "dots": "More", "ruler": "Measure",
             "sliders": "Edit", "layers": "Ghost"}
    acts = "".join(f'<button aria-label="{names.get(a, a)}" title="{names.get(a, a)}" style="width: 30px; height: 30px; padding: 0; border: 0; border-radius: 15px; '
                   f'background: transparent; cursor: pointer; display: flex; align-items: center; justify-content: center">{icon(a, 15, TEXT)}</button>' for a in actions)
    return (f'<div role="toolbar" aria-label="Selection" style="display: inline-flex; align-items: center; gap: 2px; height: 38px; padding: 0 4px 0 12px; '
            f'background: {BG2}; border: 1px solid {SEL}; border-radius: 19px; {SHADOW}">'
            f'<span style="font-family: {MONO}; font-size: 12.5px; color: {SEL}; font-weight: 500">{count}</span>'
            f'<span style="font-size: 12px; color: {MUTED}; margin: 0 8px 0 6px">{what} · {how}</span>'
            f'<span style="width: 1px; height: 20px; background: {LINE}; margin-right: 4px"></span>{acts}</div>')


def ring_menu(cx, cy, items, hot=0, r=86):
    """A radial menu around the cursor: wedges with an icon and a word; the one under the pointer lit."""
    import math as m
    n = len(items)
    svg = [f'<svg width="{2 * r + 80}" height="{2 * r + 80}" viewBox="0 0 {2 * r + 80} {2 * r + 80}" aria-hidden="true">']
    c = r + 40
    svg.append(f'<circle cx="{c}" cy="{c}" r="{r + 26}" fill="{BG2}f2" stroke="{LINE}"></circle>')
    for k in range(n):
        a0, a1 = m.radians(-90 + 360 * k / n - 180 / n), m.radians(-90 + 360 * k / n + 180 / n)
        R, r0 = r + 26, 30
        p = (f'M{c + r0 * m.cos(a0):.1f} {c + r0 * m.sin(a0):.1f} L{c + R * m.cos(a0):.1f} {c + R * m.sin(a0):.1f} '
             f'A{R} {R} 0 0 1 {c + R * m.cos(a1):.1f} {c + R * m.sin(a1):.1f} L{c + r0 * m.cos(a1):.1f} {c + r0 * m.sin(a1):.1f} A{r0} {r0} 0 0 0 {c + r0 * m.cos(a0):.1f} {c + r0 * m.sin(a0):.1f}Z')
        svg.append(f'<path d="{p}" fill="{ACC + "33" if k == hot else "transparent"}" stroke="{LINE}" stroke-width="1"></path>')
    svg.append(f'<circle cx="{c}" cy="{c}" r="26" fill="{BG1}" stroke="{SEL}"></circle></svg>')
    labels = ""
    for k, (ic, word) in enumerate(items):
        a = m.radians(-90 + 360 * k / n)
        x, y = c + (r - 4) * m.cos(a), c + (r - 4) * m.sin(a)
        col = ACC if k == hot else TEXT
        labels += at(x - 34, y - 18, f'<div style="width: 68px; display: flex; flex-direction: column; align-items: center; gap: 2px">{icon(ic, 16, col)}'
                                     f'<span style="font-size: 10.5px; color: {col}; white-space: nowrap">{word}</span></div>')
    centre = at(c - 26, c - 9, f'<div style="width: 52px; text-align: center; font-family: {MONO}; font-size: 11px; color: {SEL}">158</div>')
    return at(cx - c, cy - c, f'<div style="position: relative">{"".join(svg)}{labels}{centre}</div>')


def hist(w, h, bars, lo=None, hi=None, color=SEL, label="", ticks=()):
    """A small histogram with a brushed range [lo, hi] (bar indices)."""
    mx = max(bars) or 1
    bw = w / len(bars)
    s = [f'<svg width="{w}" height="{h + 16}" viewBox="0 0 {w} {h + 16}" aria-label="{label}">']
    if lo is not None:
        s.append(f'<rect x="{lo * bw:.1f}" y="0" width="{(hi - lo + 1) * bw:.1f}" height="{h}" fill="{color}1f" stroke="{color}" stroke-width="1"></rect>')
    for k, v in enumerate(bars):
        on = lo is not None and lo <= k <= hi
        bh = (h - 4) * v / mx
        s.append(f'<rect x="{k * bw + 1:.1f}" y="{h - bh:.1f}" width="{bw - 2:.1f}" height="{bh:.1f}" rx="1.5" fill="{color if on else "#3A4148"}"></rect>')
    for x, t in ticks:
        anchor = "start" if x < 12 else "end" if x > w - 12 else "middle"
        s.append(f'<text x="{x:.1f}" y="{h + 12}" fill="{DIM}" font-size="10" font-family="IBM Plex Mono" text-anchor="{anchor}">{t}</text>')
    s.append('</svg>')
    return "".join(s)


def spark(w, h, ys, color=ACC):
    mx, mn = max(ys), min(ys)
    pts = " ".join(f'{w * k / (len(ys) - 1):.1f},{h - (h - 2) * (y - mn) / ((mx - mn) or 1) - 1:.1f}' for k, y in enumerate(ys))
    return f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" aria-hidden="true"><polyline points="{pts}" fill="none" stroke="{color}" stroke-width="1.4"></polyline></svg>'


def workspace_switch(active="Sketch"):
    return (f'<div role="group" aria-label="Workspace" style="display: inline-flex; gap: 2px; padding: 2px; background: {BG0}; border: 1px solid {LINE}; border-radius: 15px">'
            + "".join(f'<button aria-pressed="{"true" if w == active else "false"}" style="height: 24px; padding: 0 10px; border-radius: 12px; font-size: 11.5px; cursor: pointer; '
                      f'{"background: " + ACC + "; color: " + ACC_INK + "; border: 0; font-weight: 600" if w == active else "background: transparent; color: " + MUTED + "; border: 0"}">{w}</button>'
                      for w in ("Sketch", "Build", "Analyse", "Present", "Mine"))
            + '</div>')


SKETCH_TOOLS = [("cursor", "Pick"), ("lasso", "Lasso"), ("bond", "Bond"), ("ring", "Ring"), ("plus", "Add H"), ("wand", "Clean up")]
SELECT_TOOLS = [("filter", "Brush select"), ("tag", "Tags"), ("eye", "Ghost / hide"), ("atom", "Probe"), ("ruler", "Measure")]


# ------------------------------------------------------------------ S1 · tool shelves and workspaces
def shelves():
    # the workspace switch replaces nothing: it sits in the top bar beside the compute chip
    top = cover(760, 40, 150, 28, BG1) + at(640, 40, workspace_switch("Sketch"))
    # floating Select shelf being dragged; snap guide on the left edge of the view
    drag = at(560, 300, shelf(SELECT_TOOLS, "Select", "filter", active="Brush select"))
    guide = (at(302, 214, "", 4, 427, f"background: {ACC}; border-radius: 2px; box-shadow: 0 0 14px {ACC}")
             + at(314, 400, chip("snap: left edge · vertical", ACC, "#2E2618")))
    parked = at(1060, 222, f'<div style="display: flex; flex-direction: column; gap: 8px">{puck("ruler", "Measure", 4)}{puck("layers", "Style", 6)}{puck("atom", "Probes", 3)}</div>')
    return board("CAPS — Tool shelves and workspaces",
                 top, drag, guide, parked,
                 at(1120, 230, callout(1, "Tools live on shelves: capsules with an amber tab. Drag the tab to tear a shelf off and float it anywhere, over the view or outside the window.", 290)),
                 at(1120, 340, callout(2, "Near an edge a light guide shows where it will snap; it turns vertical on the side edges.", 290)),
                 at(1120, 428, callout(3, "Click the tab to fold a shelf into a puck parked on the edge (the number: its tools). Hover opens it in place.", 290)),
                 at(1120, 530, callout(4, "Workspaces (Sketch, Build, Analyse, Present, Mine) swap which shelves are out and where. Mine is yours; every layout is saved.", 290)))


def shelf_editor():
    tools = [("cursor", "Pick"), ("lasso", "Lasso"), ("bond", "Bond"), ("ring", "Ring"), ("plus", "Add H"), ("wand", "Clean up"), ("filter", "Brush select"),
             ("tag", "Tags"), ("eye", "Ghost / hide"), ("atom", "Probe"), ("ruler", "Measure"), ("copy", "Copy"), ("mirror", "Mirror"), ("rotate", "Rotate"),
             ("move", "Move"), ("scissors", "Break bond"), ("history", "History"), ("layers", "Style")]
    grid = "".join(f'<div style="display: flex; flex-direction: column; align-items: center; gap: 4px; width: 66px">{tbtn(i, n)}'
                   f'<span style="font-size: 10.5px; color: {MUTED}; text-align: center">{n}</span></div>' for i, n in tools)
    mine = shelf([("filter", "Brush select"), ("tag", "Tags"), ("atom", "Probe"), ("wand", "Clean up")], "My interface shelf", "pin")
    body = (f'<div style="display: flex; align-items: center; justify-content: space-between"><h2 style="margin: 0; font-size: 15px; font-weight: 600">Make a shelf</h2>'
            f'{icon("close", 16, MUTED)}</div>'
            f'<span style="font-size: 12px; color: {MUTED}">Drag tools onto the shelf below; drag off to remove. Name it, pick its glyph.</span>'
            f'<div style="display: flex; flex-wrap: wrap; gap: 10px 4px">{grid}</div>'
            f'<div style="padding: 14px; border: 1px dashed {ACC}; border-radius: 10px; display: flex; align-items: center; gap: 12px">{mine}'
            f'<span style="font-size: 11.5px; color: {DIM}">drop here</span></div>'
            + row(field("Name", "My interface shelf", None, False), select("Size", "Medium · 36 px", 160), select("Labels", "on hover", 140), gap=10)
            + f'<div style="display: flex; gap: 8px; justify-content: flex-end">{btn("Reset workspace", False, None, True)}{btn("Done", True, None, True)}</div>')
    vert = at(304, 222, shelf(SKETCH_TOOLS, "Sketch", "bond", vertical=True, active="Pick"))
    return board("CAPS — Make your own shelf",
                 vert, at(0, 0, "", 1440, 900, "background: #0B0D0F80"), at(430, 170, panel(body, 600, pad=18, gap=12)),
                 at(1060, 230, callout(1, "Any tool can go on any shelf; a shelf can be yours alone. Here the Sketch shelf is snapped vertically to the view's left edge.", 300)))


# ------------------------------------------------------------------ S2 · A1 · A8 · A11 · the selection bar and the ring menu
def selection_bar():
    bar = at(476, 410, selbar("158", "molecule 4", "lasso"))
    ring = ring_menu(985, 322, [("eye", "Hide"), ("layers", "Ghost"), ("tag", "Tag"), ("atom", "Probe"), ("copy", "Copy"), ("wand", "Clean"), ("sliders", "Edit"), ("xcircle", "Delete")], hot=1)
    status = cover(0, 876, 1440, 24, BG1) + at(14, 880, f'<span style="font-family: {MONO}; font-size: 11px; color: {MUTED}">1,426 atoms · '
                                                         f'<span style="color: {SEL}">158 selected · molecule 4</span> · click pick · double-click molecule · ⌥ double-click same type · '
                                                         f'⌥ drag box · right-press ring · Esc clear</span>')
    return board("CAPS — Selection bar and ring menu",
                 lasso_blue(), bar, ring, status,
                 at(1120, 520, callout(1, "One selection, whatever made it (click, box, lasso, brush, tag). A bar floats above it with its count, what it is and how it was made, "
                                          "and the actions; every command in the shelves and menus acts on it.", 290)),
                 at(1120, 670, callout(2, "Press and hold the right button: a ring opens under the pointer; flick toward an action and release. "
                                          "Ghost keeps the atoms faint in the view, out of picking.", 290)))


# ------------------------------------------------------------------ S3 · A2 · A9 · layers with ghost, hide and lock
def layers():
    rows = [("h-BN slab", "molecule 1 · B N", "320", "#E0A060", "hidden", 0), ("polyisoprene film", "7 chains", "1,106", SEL, "shown", 0),
            ("chain 2", "C H", "158", SEL, "shown", 1), ("chain 3", "C H", "158", SEL, "shown", 1), ("chain 4", "C H", "158", SEL, "sel", 1),
            ("chain 5", "C H", "158", SEL, "ghost", 1), ("chain 6", "C H", "158", SEL, "shown", 1)]
    def state_icons(st):
        eye = icon("eye", 14, DIM if st == "hidden" else (MUTED if st != "ghost" else ACC))
        lock = icon("pin", 13, DIM)
        return (f'<button aria-label="Show, ghost or hide" style="width: 22px; height: 22px; padding: 0; border: 0; background: transparent; cursor: pointer">{eye}</button>'
                f'<button aria-label="Lock against picking and editing" style="width: 22px; height: 22px; padding: 0; border: 0; background: transparent; cursor: pointer">{lock}</button>')
    zs = [[0, 0, 9, 0, 0, 0, 0, 0], [0, 1, 3, 6, 8, 7, 5, 2], [0, 0, 2, 5, 7, 4, 1, 0]]
    items = ""
    for k, (n, sub, cnt, c, st, lvl) in enumerate(rows):
        bg = f"background: {BG3}; " if st == "sel" else ""
        op = "opacity: 0.5; " if st == "hidden" else ""
        tag = chip("ghost", ACC, "#2E2618") if st == "ghost" else (chip("hidden", DIM, BG2) if st == "hidden" else "")
        items += (f'<div style="{bg}{op}height: 30px; display: flex; align-items: center; gap: 7px; padding: 0 6px 0 {8 + 14 * lvl}px; border-radius: 6px">'
                  f'<span style="width: 10px; height: 10px; border-radius: 3px; background: {c}"></span>'
                  f'<span style="font-size: 12.5px; flex-grow: 1; white-space: nowrap; overflow: hidden; text-overflow: ellipsis">{n}</span>{tag}'
                  f'{spark(40, 14, zs[k % 3], c)}<span style="font-family: {MONO}; font-size: 11px; color: {DIM}; width: 38px; text-align: right">{cnt}</span>{state_icons(st)}</div>')
    lay = cover(0, 300, 298, 576) + at(0, 842, f'<div style="width: 298px; height: 34px; box-sizing: border-box; display: flex; align-items: center; gap: 8px; padding: 0 14px; '
                                                 f'border-top: 1px solid {LINE}">{icon("chevr", 13, MUTED)}<span style="font-size: 12.5px; color: {MUTED}">Properties · Fragments</span></div>') + at(6, 308, f'<div style="width: 288px; display: flex; flex-direction: column; gap: 1px">'
                                                 f'<div style="display: flex; align-items: center; justify-content: space-between; padding: 0 6px 4px"><h3 style="{H3}">Layers</h3>'
                                                 f'<span style="font-size: 11px; color: {DIM}">z profile · atoms · state</span></div>{items}</div>')
    ghost_view = at(575, 318, "", 120, 230, f"background: radial-gradient(closest-side, {BG0} 60%, transparent 100%); opacity: 0.85")
    hud = at(318, 262, chip(f'h-BN slab hidden · chain 5 ghosted · <span style="color: {ACC}">show everything</span>', TEXT, BG2))
    return board("CAPS — Layers: show, ghost, hide, lock",
                 lay, ghost_view, hud,
                 at(1120, 230, callout(1, "The structure as layers: what it is made of (the slab, the film and its chains), grouped by kind, each with a small profile along z, "
                                          "its atom count and its state.", 290)),
                 at(1120, 370, callout(2, "The eye cycles shown → ghost (faint, not pickable) → hidden. The pin locks a layer against picking and edits: "
                                          "the held slab, for instance. Nothing is removed from the structure, the exports or the calculations.", 290)),
                 at(1120, 530, callout(3, "Selecting in the view lights the layer rows; clicking a row selects it in the view.", 290)))


# ------------------------------------------------------------------ S4 · A3 · tags
def tags():
    tg = [("slab", "320", "#E0A060"), ("film", "1,106", SEL), ("interphase", "214", "#9B7AD5"), ("chain ends", "14", "#7CC784"), ("+ new tag", "", None)]
    pills = "".join(
        (f'<button style="height: 28px; padding: 0 10px; display: inline-flex; align-items: center; gap: 6px; border-radius: 14px; background: {BG2}; border: 1px solid {LINE}; '
         f'font-size: 12px; cursor: pointer">{dot(c, 9)}{n}<span style="font-family: {MONO}; font-size: 10.5px; color: {DIM}">{cnt}</span></button>') if c else
        f'<button style="height: 28px; padding: 0 10px; border-radius: 14px; background: transparent; border: 1px dashed {LINE}; font-size: 12px; color: {MUTED}; cursor: pointer">{n}</button>'
        for n, cnt, c in tg)
    tray = cover(312, 222, 170, 28, BG0) + at(318, 222, f'<div style="display: flex; gap: 6px; align-items: center"><span style="font-size: 11px; font-weight: 600; letter-spacing: 0.08em; color: {MUTED}">TAGS</span>{pills}</div>')
    bar = at(560, 404, selbar("214", "tagged interphase", "brush z 7–11 Å", ("tag", "eye", "atom", "copy", "dots")))
    band = at(560, 380, "", 290, 34, "background: #9B7AD533; border-top: 1px dashed #9B7AD5; border-bottom: 1px dashed #9B7AD5; transform: rotate(-24deg)")
    def urow(k, v, code=True):
        vv = (f'<code style="font-family: {MONO}; font-size: 11.5px; color: {TEXT}">{v}</code>' if code else f'<span style="font-size: 12px; color: {TEXT}">{v}</span>')
        return (f'<div style="display: grid; grid-template-columns: 110px 1fr; gap: 10px; align-items: baseline; padding: 5px 0; border-bottom: 1px solid {LINE}">'
                f'<span style="font-size: 12px; color: {MUTED}">{k}</span>{vv}</div>')
    use = panel(f'<h3 style="{H3}">Where a tag goes</h3>'
                + urow("LAMMPS", "group interphase id 322:341 418 …") + urow("GROMACS", "[ interphase ] in index.ndx")
                + urow("Analyze", "density profile, RDF, MSD, orientation of a tag", False)
                + urow("Python", "doc.tags[\"interphase\"]") + urow("Recipes", "analyze: { of: interphase }"), 480, gap=2)
    return board("CAPS — Tags",
                 tray, band, bar, at(318, 640, use),
                 at(1120, 230, callout(1, "Tag any selection with a name and colour. Tags sit in a strip over the view: click to select, ⇧ click to add, "
                                          "⌥ click to show only that tag.", 290)),
                 at(1120, 370, callout(2, "An atom can carry several tags (a chain end inside the interphase). Tags are saved in the structure and become groups everywhere.", 290)))


# ------------------------------------------------------------------ S5 · A4 · brush select
def brush_select():
    z = [9, 12, 30, 38, 22, 16, 18, 20, 21, 19, 20, 18, 17, 14, 8, 3]
    q = [2, 5, 18, 40, 60, 38, 12, 4, 2, 1]
    bonds = [14, 0, 30, 90, 0, 0]
    hy = [0, 22, 88]
    card = lambda title, inner, note="": (f'<div style="display: flex; flex-direction: column; gap: 6px; padding: 10px 12px; background: {BG2}; border: 1px solid {LINE}; border-radius: 8px">'
                                          f'<div style="display: flex; justify-content: space-between; font-size: 12px"><span>{title}</span><span style="color: {SEL}; font-family: {MONO}; font-size: 11px">{note}</span></div>{inner}</div>')
    body = (f'<div style="display: flex; align-items: center; justify-content: space-between"><h3 style="{H3}">Brush to select</h3>{chip("214 of 1,426", SEL, "#10262B", True)}</div>'
            + card("Height z (Å)", hist(272, 54, z, 3, 5, ticks=[(0, "0"), (68, "8"), (136, "16"), (204, "24"), (272, "31")]), "7.0 – 11.0")
            + card("Partial charge (e)", hist(272, 40, q, ticks=[(0, "−0.3"), (136, "0"), (272, "+0.3")]))
            + card("Bonds per atom", hist(272, 34, bonds, ticks=[(23, "1"), (68, "2"), (113, "3"), (159, "4"), (204, "5"), (249, "6")]))
            + card("Hybridisation", hist(272, 30, hy, ticks=[(45, "sp"), (136, "sp²"), (227, "sp³")]))
            + f'<code style="padding: 7px 9px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px; font-family: {MONO}; font-size: 11.5px; color: {MUTED}">'
              f'z 7..11 and not layer "h-BN slab"</code>'
            + f'<div style="display: flex; gap: 6px">{btn("Tag it", True, "tag", True)}{btn("Add more", False, None, True)}{btn("Distance to…", False, None, True)}</div>')
    band = at(560, 380, "", 290, 34, f"background: {SEL}26; border-top: 1px dashed {SEL}; border-bottom: 1px dashed {SEL}; transform: rotate(-24deg)")
    return board("CAPS — Brush to select",
                 cover(1100, 212, 340, 664), at(1110, 220, panel(body, 322, pad=12, gap=8, extra="box-shadow: none")), band,
                 at(318, 262, callout(1, "Each property of the atoms as a histogram: drag across one to select that range, across several to intersect. "
                                        "The view rings the match live; the line under them is the same selection as text, editable, reusable in recipes.", 320)),
                 at(318, 420, callout(2, "Distance to… adds a histogram of distance from a picked atom, a probe or a tag (e.g. rubber within 6 Å of the filler).", 320)))


# ------------------------------------------------------------------ S6 · A5 · clipboard tray, stamp, drop into
def stamp():
    ghost = at(470, 250, f'<svg width="170" height="120" viewBox="0 0 170 120" aria-hidden="true"><g fill="none" stroke="{ACC}" stroke-width="2">'
                         f'<path d="M20 90 L50 70 L80 88 L110 66 L140 84"></path><path d="M50 70 L52 40"></path><path d="M110 66 L112 36"></path></g>'
                         + "".join(f'<circle cx="{x}" cy="{y}" r="6" fill="{ACC}55" stroke="{ACC}"></circle>' for x, y in [(20, 90), (50, 70), (80, 88), (110, 66), (140, 84), (52, 40), (112, 36)])
                         + '</svg>')
    tip = at(440, 378, chip(f'stamp · click to place · scroll turns · ⌥ place many · 1.5 Å clear of the film', TEXT, BG2))
    clips = [("isoprene unit", "13 atoms · from iface.data"), ("TESPT", "46 atoms · TESPT.mol2"), ("chain end ×3", "39 atoms"), ("ZnO cluster", "68 atoms")]
    tray = "".join(f'<button style="width: 118px; height: 104px; box-sizing: border-box; flex-shrink: 0; display: flex; flex-direction: column; gap: 4px; padding: 8px; background: {BG2 if k else "#2E2618"}; '
                   f'border: 1px solid {ACC if k == 0 else LINE}; border-radius: 8px; cursor: pointer; text-align: left">'
                   f'<svg width="100" height="40" viewBox="0 0 100 40" aria-hidden="true"><path d="M8 30 L28 18 L48 30 L68 18 L90 28" fill="none" stroke="{MUTED}" stroke-width="2"></path></svg>'
                   f'<span style="font-size: 12px">{n}</span><span style="font-size: 10.5px; color: {DIM}">{d}</span></button>' for k, (n, d) in enumerate(clips))
    trayp = at(300, 690, f'<div style="width: 800px; box-sizing: border-box; display: flex; align-items: center; gap: 8px; padding: 10px 12px; background: {BG1}; border-top: 1px solid {LINE}">'
                         f'<span style="font-size: 11px; font-weight: 600; letter-spacing: 0.08em; color: {MUTED}; writing-mode: vertical-rl; transform: rotate(180deg)">CLIPBOARD</span>{tray}</div>')
    drop = (cover(1100, 212, 340, 664, "#0B0D0Fd9") + at(1110, 222, f'<div style="width: 318px; height: 150px; box-sizing: border-box; border: 2px dashed {ACC}; border-radius: 12px; background: #2E2618; display: flex; '
                          f'flex-direction: column; align-items: center; justify-content: center; gap: 6px">{icon("plus", 22, ACC)}<b style="font-size: 13px">Drop into iface.data</b>'
                          f'<span style="font-size: 11.5px; color: {MUTED}">as a stamp to place · undoable</span></div>')
            + at(1110, 386, f'<div style="width: 318px; height: 90px; box-sizing: border-box; border: 2px dashed {LINE}; border-radius: 12px; background: {BG1}; display: flex; '
                            f'align-items: center; justify-content: center; font-size: 12.5px; color: {MUTED}">Open as a new structure</div>'))
    return board("CAPS — Clipboard tray, stamps, drop into",
                 ghost, tip, trayp, drop,
                 at(1110, 500, callout(1, "Copy keeps a tray of the last pieces (atoms, bonds, types, charges). Pick one and it becomes a stamp under the pointer: click to place, "
                                          "⌥ to place many, pushed clear of what is there.", 318)),
                 at(1110, 650, callout(2, "Drag a file onto the window: two targets, into this structure (a stamp) or open it new.", 318)))


# ------------------------------------------------------------------ S7 · A6 · A15 · the Look popover
def look():
    def ring():
        import math as m
        pts = [(60 + 34 * m.cos(m.radians(90 + 60 * k)), 58 + 34 * m.sin(m.radians(90 + 60 * k))) for k in range(6)]
        s = ""
        for k in range(6):
            (x1, y1), (x2, y2) = pts[k], pts[(k + 1) % 6]
            s += f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#8E959C" stroke-width="5" stroke-linecap="round"></line>'
            s += (f'<line x1="{60 + (x1 - 60) * 0.74:.1f}" y1="{58 + (y1 - 58) * 0.74:.1f}" x2="{60 + (x2 - 60) * 0.74:.1f}" y2="{58 + (y2 - 58) * 0.74:.1f}" '
                  f'stroke="#8E959C" stroke-width="2.5" stroke-dasharray="4 3"></line>')
        return s + "".join(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="7" fill="#8E959C"></circle>' for x, y in pts)
    diene = ('<g transform="translate(130 22)"><line x1="0" y1="44" x2="34" y2="24" stroke="#8E959C" stroke-width="5"></line>'
             '<line x1="34" y1="24" x2="72" y2="44" stroke="#8E959C" stroke-width="4"></line><line x1="37" y1="34" x2="69" y2="53" stroke="#8E959C" stroke-width="4"></line>'
             '<line x1="72" y1="44" x2="106" y2="24" stroke="#8E959C" stroke-width="5"></line>'
             + "".join(f'<circle cx="{x}" cy="{y}" r="7" fill="#8E959C"></circle>' for x, y in [(0, 44), (34, 24), (72, 44), (106, 24)]) + '</g>')
    def slider(label, val, pct):
        return (f'<div style="display: flex; flex-direction: column; gap: 5px"><div style="display: flex; justify-content: space-between; font-size: 12px"><span>{label}</span>'
                f'<span style="font-family: {MONO}; color: {MUTED}">{val}</span></div><div style="height: 4px; border-radius: 2px; background: {BG3}; position: relative">'
                f'<div style="width: {pct}%; height: 4px; border-radius: 2px; background: {ACC}"></div>'
                f'<div style="position: absolute; left: calc({pct}% - 7px); top: -5px; width: 14px; height: 14px; border-radius: 50%; background: {TEXT}"></div></div></div>')
    body = (f'<div style="display: flex; align-items: center; justify-content: space-between"><h3 style="{H3}">Look · the selection</h3>{seg(["All", "Selection"], "Selection")}</div>'
            f'<div style="padding: 8px; background: {BG0}; border: 1px solid {LINE}; border-radius: 8px"><svg width="250" height="112" viewBox="0 0 250 112" aria-label="Live preview: an aromatic ring and the isoprene C=C at these sizes">{ring()}{diene}</svg></div>'
            + slider("Atoms", "0.60 × r cov", 40) + slider("Sticks", "0.16 Å", 32) + slider("Spheres (space filling)", "1.00 × r vdW", 50) + slider("Lines", "1.5 px", 30)
            + toggle("Bond orders: double, triple, aromatic", True) + toggle("Arcs on angle and torsion monitors", True)
            + f'<div style="display: flex; gap: 6px">{btn("Keep as my look", False, "save", True)}{btn("Reset", False, None, True)}</div>')
    pop = at(760, 262, panel(body, 290, pad=12, gap=10))
    trigger = at(1044, 222, f'<button aria-label="Look" style="width: 44px; height: 30px; border-radius: 15px; border: 1px solid {ACC}; background: #2E2618; color: {ACC}; '
                            f'font-size: 11.5px; cursor: pointer">Look</button>')
    arc = at(600, 470, f'<svg width="120" height="90" viewBox="0 0 120 90" aria-hidden="true"><path d="M20 70 A48 48 0 0 1 86 26" fill="{ACC}22" stroke="{ACC}" stroke-width="1.5"></path></svg>')
    return board("CAPS — Look: bond orders and sizes",
                 pop, trigger, arc, at(650, 452, chip("116.4°", ACC, BG2, True)),
                 at(1120, 280, callout(1, "A Look button sits on the view itself. It opens over the view with a live preview (a ring and the isoprene C=C) "
                                          "at the sizes you drag, for everything or only the selection.", 290)),
                 at(1120, 420, callout(2, "Bond orders drawn as parallel sticks and dashed aromatic bonds; arcs on angle monitors. Keep as my look saves it as your default.", 290)))


# ------------------------------------------------------------------ S8 · A7 · bond rules drawn over distance histograms
def bond_rules():
    def pair(name, bars, cut, never=False, note=""):
        n = len(bars)
        w = 360
        x = w * cut / n
        h = (f'<div style="position: relative">{hist(w, 46, bars, color=ERR if never else ACC)}'
             + ("" if never else f'<div style="position: absolute; left: {x:.0f}px; top: -4px; width: 2px; height: 54px; background: {TEXT}"></div>'
                                 f'<div style="position: absolute; left: {x - 8:.0f}px; top: -10px; width: 18px; height: 10px; border-radius: 3px; background: {TEXT}"></div>')
             + '</div>')
        right = (f'<button aria-pressed="{"true" if never else "false"}" style="height: 24px; padding: 0 9px; border-radius: 12px; font-size: 11px; cursor: pointer; '
                 f'{"background: #2E1E1E; color: " + ERR + "; border: 1px solid " + ERR if never else "background: transparent; color: " + MUTED + "; border: 1px solid " + LINE}">never bond</button>')
        return (f'<div style="display: grid; grid-template-columns: 64px 1fr 110px; gap: 10px; align-items: center; padding: 8px 0; border-bottom: 1px solid {LINE}">'
                f'<span style="font-family: {MONO}; font-size: 13px">{name}</span>{h}<div style="display: flex; flex-direction: column; gap: 4px; align-items: flex-end">{right}'
                f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{note}</span></div></div>')
    body = (f'<div style="display: flex; align-items: center; justify-content: space-between"><h2 style="margin: 0; font-size: 15px; font-weight: 600">Bond rules · iface.data</h2>{icon("close", 16, MUTED)}</div>'
            f'<span style="font-size: 12px; color: {MUTED}; line-height: 1.5">Every pair of elements in the structure, with how many atom pairs sit at each distance (0.8 – 3.2 Å). '
            f'Drag the white handle to set the bond cut-off; the count updates live.</span>'
            + pair("C – C", [0, 0, 0, 2, 30, 58, 12, 0, 0, 0, 0, 3, 8, 14, 20, 22], 8, note="≤ 1.68 Å · 832")
            + pair("C – H", [0, 40, 70, 8, 0, 0, 0, 0, 0, 0, 0, 0, 2, 6, 10, 12], 6, note="≤ 1.32 Å · 686")
            + pair("B – N", [0, 0, 0, 0, 12, 60, 8, 0, 0, 0, 0, 0, 0, 4, 12, 18], 9, note="≤ 1.76 Å · 480")
            + pair("Zn – O", [0, 0, 0, 0, 0, 0, 10, 30, 22, 4, 0, 0, 0, 3, 6, 9], 0, True, "ionic: none")
            + f'<div style="display: flex; align-items: center; gap: 10px; padding: 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 8px; font-size: 12px">'
              f'{icon("bond", 16, ACC)}<span>1,579 bonds → <b>1,612</b> · 33 B–N at the slab\'s edges · the new ones dashed in the view</span></div>'
            + row(select("Rule set", "Rubber + oxide fillers", 220), btn("Save rule set", False, "save", True), gap=8)
            + f'<div style="display: flex; gap: 8px; justify-content: flex-end">{btn("Close", False, None, True)}{btn("Bond it · undoable", True, None, True)}</div>')
    return board("CAPS — Bond rules",
                 at(0, 0, "", 1440, 900, "background: #0B0D0F99"), at(400, 110, panel(body, 620, pad=18, gap=8)),
                 at(1060, 230, callout(1, "Bonds from what the structure actually holds: per element pair, a histogram of distances; the cut-off sits in the gap between "
                                          "bonded and non-bonded pairs, where you can see it. Ionic pairs (a ZnO filler) set to never bond.", 300)))


# ------------------------------------------------------------------ S9 · A10 · probes
def probes():
    view = (at(680, 300, f'<svg width="260" height="230" viewBox="0 0 260 230" aria-hidden="true">'
                         f'<polygon points="10,150 170,90 250,150 90,210" fill="#F0A83C1f" stroke="{ACC}" stroke-width="1.4"></polygon>'
                         f'<line x1="130" y1="10" x2="130" y2="150" stroke="#9B7AD5" stroke-width="2" stroke-dasharray="7 4"></line>'
                         f'<circle cx="130" cy="60" r="7" fill="{SEL}" stroke="{BG0}" stroke-width="2"></circle>'
                         f'<line x1="130" y1="60" x2="130" y2="150" stroke="{TEXT}" stroke-width="1" stroke-dasharray="2 3"></line></svg>'))
    def pcard(title, value, c, ys, sub):
        return (f'<div style="width: 180px; display: flex; flex-direction: column; gap: 4px; padding: 9px 11px; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px; {SHADOW}">'
                f'<div style="display: flex; align-items: center; gap: 6px">{dot(c, 8)}<span style="font-size: 11.5px; color: {MUTED}; flex-grow: 1">{title}</span>{icon("pin", 12, DIM)}</div>'
                f'<span style="font-family: {MONO}; font-size: 18px">{value}</span>{spark(158, 26, ys, c)}'
                f'<span style="font-size: 10.5px; color: {DIM}">{sub}</span></div>')
    tray = at(318, 520, f'<div style="display: flex; gap: 10px">'
                        + pcard("chain 4 centre ↑ slab", "8.42 Å", SEL, [9.1, 8.9, 8.6, 8.7, 8.5, 8.4, 8.4, 8.3], "centroid → plane · per frame")
                        + pcard("chain 4 axis ∠ slab", "21.7°", "#9B7AD5", [35, 31, 28, 26, 24, 22, 22, 21], "axis ∠ plane")
                        + pcard("slab flatness", "0.00 Å", ACC, [0, 0, 0, 0, 0, 0, 0, 0], "plane rms · held") + '</div>')
    bar = at(560, 404, selbar("158", "molecule 4", "layer click", ("atom", "tag", "eye", "dots")))
    menu_p = at(885, 262, panel(f'<h3 style="{H3}">Probe from the selection</h3>'
                                + "".join(f'<button style="height: 28px; display: flex; align-items: center; gap: 8px; padding: 0 8px; border: 0; border-radius: 5px; background: {BG3 if k == 0 else "transparent"}; '
                                          f'cursor: pointer; font-size: 12.5px">{icon(i, 14, MUTED)}{t}</button>'
                                          for k, (i, t) in enumerate([("atom", "Centre (mass-weighted)"), ("layers", "Best-fit plane"), ("ruler", "Long axis"), ("ring", "Enclosing ellipsoid")])),
                                210, pad=8, gap=2))
    return board("CAPS — Probes",
                 view, bar, menu_p, tray, at(830, 346, chip("centre · chain 4", SEL, BG2)), at(880, 470, chip("plane · slab top", ACC, BG2)),
                 at(1120, 230, callout(1, "A probe is a point, plane, axis or ellipsoid made from the selection (the Probe action on its bar). Probes draw in the view and measure "
                                          "against each other or atoms.", 290)),
                 at(1120, 380, callout(2, "Each measurement pins a card with its value and a sparkline over the trajectory's frames; a card goes to Analyze as a curve with one click.", 290)))


# ------------------------------------------------------------------ S10 · A12 · A13 · A14 · the atom bubble, H autopilot, orient puck
def atom_bubble():
    els = ["C", "H", "N", "O", "S", "Si", "B", "Zn"]
    el_ring = "".join(f'<button style="width: 30px; height: 30px; border-radius: 15px; font-size: 12px; font-weight: 600; cursor: pointer; '
                      f'{"background: " + ACC + "; color: " + ACC_INK + "; border: 0" if e == "C" else "background: " + BG3 + "; color: " + TEXT + "; border: 1px solid " + LINE}">{e}</button>' for e in els)
    def ef(label, val, mixed=False):
        return (f'<label style="display: flex; align-items: center; gap: 8px; font-size: 12px"><span style="width: 66px; color: {MUTED}">{label}</span>'
                f'<input value="{val}" style="flex-grow: 1; height: 26px; box-sizing: border-box; padding: 0 8px; background: {BG0}; border: 1px solid {LINE}; border-radius: 5px; '
                f'font-family: {MONO}; font-size: 12px; color: {DIM if mixed else TEXT}"></label>')
    bubble = (f'<div style="width: 318px; display: flex; flex-direction: column; gap: 8px; padding: 12px; background: {BG1}; border: 1px solid {ACC}; border-radius: 12px; {SHADOW}">'
              f'<div style="display: flex; align-items: center; justify-content: space-between"><span style="font-size: 12.5px; font-weight: 600">3 atoms · chain 4</span>{chip("C2 CT", MUTED, BG2, True)}</div>'
              f'<div style="display: flex; flex-wrap: wrap; gap: 5px">{el_ring}</div>'
              + ef("Type", "mixed: CT, CM", True) + ef("Charge", "mixed", True) + ef("Formal", "0") + ef("Nudge", "0  0  +0.5 Å")
              + f'<span style="font-size: 11px; color: {DIM}">Enter applies · each change one undo step</span></div>')
    tail = at(716, 470, f'<svg width="20" height="30" viewBox="0 0 20 30" aria-hidden="true"><path d="M10 0 L10 30" stroke="{ACC}" stroke-width="1.5"></path></svg>')
    auto = at(318, 222, f'<div style="display: inline-flex; align-items: center; gap: 8px; height: 30px; padding: 0 12px; background: {BG1}; border: 1px solid {LINE}; border-radius: 15px; {SHADOW}">'
                        f'{icon("plus", 14, ACC)}<span style="font-size: 12px">H autopilot</span>'
                        f'<span style="width: 30px; height: 16px; border-radius: 8px; background: {ACC}; position: relative"><span style="position: absolute; right: 2px; top: 2px; width: 12px; height: 12px; border-radius: 6px; background: {ACC_INK}"></span></span>'
                        f'<span style="font-size: 11px; color: {DIM}">+3 H, −1 H after the last edit</span></div>')
    orient = at(1010, 540, f'<div style="display: flex; flex-direction: column; align-items: center; gap: 6px">'
                           f'<svg width="76" height="76" viewBox="0 0 76 76" aria-label="Orientation puck"><circle cx="38" cy="38" r="36" fill="{BG1}" stroke="{LINE}"></circle>'
                           f'<line x1="38" y1="38" x2="64" y2="46" stroke="#FF7B72" stroke-width="2.5"></line><line x1="38" y1="38" x2="30" y2="10" stroke="#7CC784" stroke-width="2.5"></line>'
                           f'<line x1="38" y1="38" x2="16" y2="54" stroke="#6CA8FF" stroke-width="2.5"></line><circle cx="38" cy="38" r="5" fill="{SEL}"></circle></svg>'
                           f'<div style="display: flex; gap: 4px">{chip("along", SEL, "#10262B")}{chip("onto", MUTED)}{chip("side", MUTED)}</div></div>')
    return board("CAPS — Atom bubble, H autopilot, orient puck",
                 lasso_blue(), at(560, 250, bubble), tail, auto, orient,
                 at(1120, 230, callout(1, "Double-click atoms to open a bubble at them: element ring, type, charge, formal charge and a nudge. "
                                          "Several atoms: values that differ read mixed and stay each atom's own.", 290)),
                 at(1120, 370, callout(2, "H autopilot: after each sketch edit hydrogens are added where missing and removed where surplus; the chip says what it did.", 290)),
                 at(1120, 480, callout(3, "The orient puck in the view's corner: drag to turn, or snap the camera along the selection's long axis, onto its plane or side-on.", 290)))


BOARDS = [("Shelves", shelves, "Tool shelves & workspaces"), ("ShelfEditor", shelf_editor, "Make your own shelf"),
          ("SelectionBar", selection_bar, "Selection bar & ring menu (A1 A8 A11)"), ("Layers", layers, "Layers: show, ghost, hide, lock (A2 A9)"),
          ("Tags", tags, "Tags (A3)"), ("BrushSelect", brush_select, "Brush to select (A4)"), ("Stamp", stamp, "Clipboard tray, stamps, drop into (A5)"),
          ("Look", look, "Look: bond orders & sizes (A6 A15)"), ("BondRules", bond_rules, "Bond rules (A7)"), ("Probes", probes, "Probes (A10)"),
          ("AtomBubble", atom_bubble, "Atom bubble, H autopilot, orient puck (A12 A13 A14)")]

if __name__ == "__main__":
    out = os.environ.get("OUT", "stage29/project")
    os.makedirs(out, exist_ok=True)
    for name, fn, _ in BOARDS:
        h = fn(); open(os.path.join(out, f"{name}.dc.html"), "w").write(h); print(name, len(h))
