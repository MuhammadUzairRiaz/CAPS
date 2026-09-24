"""Shared helpers for generating CAPS Studio artboards (.dc.html)."""
import math, random

W, H = 1440, 900

# ---- palette -------------------------------------------------------------
BG0 = "#0F1113"   # viewport / app ground
BG1 = "#16191C"   # panels
BG2 = "#1E2226"   # raised, inputs
BG3 = "#262B30"   # hover / selected rows
LINE = "#2C3137"
TEXT = "#E8E6E1"
MUTED = "#A5ABB1"
DIM = "#868D94"
ACC = "#F0A83C"     # instrument amber
ACC_INK = "#1A1206"
SEL = "#6CC4D8"     # selection cyan
OK = "#7CC784"
WARN = "#F0A83C"
ERR = "#FF7B72"
MONO = "'IBM Plex Mono',ui-monospace,monospace"

ELEM = {
    "C": ("#8E959C", 0.36), "H": ("#E9ECEF", 0.22), "O": ("#E5534B", 0.35),
    "N": ("#4C7BD9", 0.35), "Si": ("#D6A45E", 0.48), "Na": ("#9B7BD6", 0.55),
    "Cl": ("#57B26A", 0.52), "S": ("#E3C74A", 0.45),
}


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def page(title, body, w=W, h=H):
    props = '{"$preview":{"width":%d,"height":%d}}' % (w, h)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<script src="./support.js"></script>
</head>
<body>
<x-dc>
<helmet>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&amp;family=IBM+Plex+Sans:wght@400;500;600;700&amp;display=swap" rel="stylesheet">
<style>
body{{margin:0;background:{BG0};color:{TEXT};font-family:'IBM Plex Sans',system-ui,sans-serif}}
a{{color:{ACC}}}a:hover{{color:#FFC46B}}
button,input,select{{font:inherit;color:inherit}}
</style>
</helmet>
<div style="width: {w}px; height: {h}px; box-sizing: border-box; display: flex; flex-direction: column; background: {BG0}; color: {TEXT}; font-family: 'IBM Plex Sans', system-ui, sans-serif; font-size: 13px; overflow: hidden">
{body}
</div>
</x-dc>
<script type="text/x-dc" data-dc-script data-props='{props}'>
class Component extends DCLogic {{
renderVals() {{
return {{}};
}}
}}
</script>
</body>
</html>
"""

# ---- icons (24x24 stroke) ------------------------------------------------
ICONS = {
    "cursor": '<path d="M5 3l14 8-6 1.5L10 19z"></path>',
    "lasso": '<ellipse cx="12" cy="9" rx="8" ry="5"></ellipse><path d="M6 12c-1 3 0 5 2 6s3 2 2 3"></path>',
    "move": '<path d="M12 3v18M3 12h18M12 3l-3 3M12 3l3 3M12 21l-3-3M12 21l3-3M3 12l3-3M3 12l3 3M21 12l-3-3M21 12l-3 3"></path>',
    "rotate": '<path d="M20 12a8 8 0 1 1-2.3-5.7"></path><path d="M20 4v5h-5"></path>',
    "atom": '<circle cx="12" cy="12" r="3"></circle><circle cx="12" cy="12" r="8"></circle>',
    "bond": '<circle cx="6" cy="18" r="2.5"></circle><circle cx="18" cy="6" r="2.5"></circle><path d="M8 16l8-8"></path>',
    "hex": '<path d="M12 3l7.8 4.5v9L12 21l-7.8-4.5v-9z"></path>',
    "ring": '<path d="M12 3l7.8 4.5v9L12 21l-7.8-4.5v-9z"></path><circle cx="12" cy="12" r="3.5"></circle>',
    "ruler": '<path d="M3 17L17 3l4 4L7 21z"></path><path d="M7 13l2 2M10 10l2 2M13 7l2 2"></path>',
    "pin": '<path d="M12 17v5M8 3h8l-1 6 3 3H6l3-3z"></path>',
    "wand": '<path d="M4 20L16 8M14 4v3M18 8h3M17 5l2-2M19 11l1 1"></path>',
    "eye": '<path d="M2 12s4-7 10-7 10 7 10 7-4 7-10 7S2 12 2 12z"></path><circle cx="12" cy="12" r="3"></circle>',
    "layers": '<path d="M12 3l9 5-9 5-9-5z"></path><path d="M3 13l9 5 9-5"></path>',
    "undo": '<path d="M9 14L4 9l5-5"></path><path d="M4 9h10a6 6 0 0 1 0 12h-3"></path>',
    "redo": '<path d="M15 14l5-5-5-5"></path><path d="M20 9H10a6 6 0 0 0 0 12h3"></path>',
    "search": '<circle cx="11" cy="11" r="7"></circle><path d="M20 20l-4-4"></path>',
    "play": '<path d="M7 4l13 8-13 8z"></path>',
    "pause": '<path d="M8 4v16M16 4v16"></path>',
    "stop": '<rect x="5" y="5" width="14" height="14" rx="1"></rect>',
    "cube": '<path d="M12 2l9 5v10l-9 5-9-5V7z"></path><path d="M12 22V12M3 7l9 5 9-5"></path>',
    "flask": '<path d="M9 3h6M10 3v6L4 19a1.5 1.5 0 0 0 1.3 2h13.4A1.5 1.5 0 0 0 20 19l-6-10V3"></path><path d="M7 15h10"></path>',
    "chart": '<path d="M3 3v18h18"></path><path d="M7 15l4-5 3 3 5-7"></path>',
    "gauge": '<path d="M4 18a8 8 0 1 1 16 0"></path><path d="M12 18l4-6"></path>',
    "check": '<path d="M4 12l5 5L20 6"></path>',
    "alert": '<path d="M12 3l10 18H2z"></path><path d="M12 10v5M12 18v.5"></path>',
    "xcircle": '<circle cx="12" cy="12" r="9"></circle><path d="M9 9l6 6M15 9l-6 6"></path>',
    "chev": '<path d="M6 9l6 6 6-6"></path>',
    "chevr": '<path d="M9 6l6 6-6 6"></path>',
    "plus": '<path d="M12 5v14M5 12h14"></path>',
    "folder": '<path d="M3 6h7l2 2h9v11H3z"></path>',
    "file": '<path d="M6 3h9l4 4v14H6z"></path><path d="M15 3v4h4"></path>',
    "link": '<path d="M10 14a4 4 0 0 0 6 0l3-3a4 4 0 0 0-6-6l-1 1"></path><path d="M14 10a4 4 0 0 0-6 0l-3 3a4 4 0 0 0 6 6l1-1"></path>',
    "cpu": '<rect x="6" y="6" width="12" height="12" rx="1"></rect><path d="M9 2v4M15 2v4M9 18v4M15 18v4M2 9h4M2 15h4M18 9h4M18 15h4"></path>',
    "server": '<rect x="3" y="4" width="18" height="7" rx="1"></rect><rect x="3" y="13" width="18" height="7" rx="1"></rect><path d="M7 7.5h.01M7 16.5h.01"></path>',
    "download": '<path d="M12 3v12M7 10l5 5 5-5M4 21h16"></path>',
    "gear": '<circle cx="12" cy="12" r="3"></circle><path d="M12 2v3M12 19v3M2 12h3M19 12h3M4.9 4.9l2.1 2.1M17 17l2.1 2.1M4.9 19.1L7 17M17 7l2.1-2.1"></path>',
    "terminal": '<path d="M4 6l6 6-6 6M12 18h8"></path>',
    "history": '<path d="M3 12a9 9 0 1 0 3-6.7L3 8"></path><path d="M3 3v5h5M12 7v5l3 3"></path>',
    "filter": '<path d="M3 4h18l-7 9v6l-4 2v-8z"></path>',
    "sliders": '<path d="M4 6h10M18 6h2M4 12h4M12 12h8M4 18h12M20 18h0"></path><circle cx="16" cy="6" r="2"></circle><circle cx="10" cy="12" r="2"></circle><circle cx="18" cy="18" r="2"></circle>',
    "grow": '<path d="M4 20c4-1 6-4 6-8s3-7 10-8"></path><circle cx="4" cy="20" r="1.5"></circle><circle cx="20" cy="4" r="1.5"></circle>',
    "pack": '<rect x="3" y="3" width="18" height="18" rx="1"></rect><circle cx="8" cy="8" r="2"></circle><circle cx="15" cy="9" r="2.5"></circle><circle cx="9" cy="15" r="2.5"></circle><circle cx="16" cy="16" r="2"></circle>',
    "relax": '<path d="M3 6c3 0 3 12 6 12s3-9 6-9 3 4 6 4"></path>',
    "dyn": '<path d="M2 12h3l3-8 4 16 3-10 2 2h5"></path>',
    "equil": '<path d="M4 18c3-10 5-10 8 0s5 10 8 0"></path><path d="M3 12h18"></path>',
    "react": '<circle cx="6" cy="12" r="3"></circle><circle cx="18" cy="12" r="3"></circle><path d="M9 12h6"></path><path d="M12 5v3M12 16v3"></path>',
    "bench": '<path d="M4 20V10M10 20V4M16 20v-8M22 20H2"></path>',
    "jobs": '<rect x="3" y="4" width="18" height="16" rx="1"></rect><path d="M7 9h10M7 13h6M7 17h8"></path>',
    "tag": '<path d="M3 3h8l10 10-8 8L3 11z"></path><circle cx="7.5" cy="7.5" r="1.5"></circle>',
    "save": '<path d="M5 3h11l3 3v15H5z"></path><path d="M8 3v5h7V3M8 21v-7h8v7"></path>',
    "copy": '<rect x="8" y="8" width="13" height="13" rx="1"></rect><path d="M16 8V3H3v13h5"></path>',
    "mirror": '<path d="M12 3v18"></path><path d="M9 7L3 17h6zM15 7l6 10h-6z"></path>',
    "scissors": '<circle cx="6" cy="6" r="3"></circle><circle cx="6" cy="18" r="3"></circle><path d="M8.5 8L20 18M8.5 16L20 6"></path>',
    "dots": '<circle cx="5" cy="12" r="1"></circle><circle cx="12" cy="12" r="1"></circle><circle cx="19" cy="12" r="1"></circle>',
    "close": '<path d="M6 6l12 12M18 6L6 18"></path>',
}


def icon(name, size=16, color="currentColor", sw=1.7):
    return (f'<svg width="{size}" height="{size}" viewBox="0 0 24 24" fill="none" stroke="{color}" '
            f'stroke-width="{sw}" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">{ICONS[name]}</svg>')


# ---- small components -----------------------------------------------------
def tbtn(ic, label, active=False, text=None, w=None):
    bg = BG3 if active else "transparent"
    col = ACC if active else MUTED
    bd = f"1px solid {LINE}" if active else "1px solid transparent"
    inner = icon(ic, 18, col)
    if text:
        inner += f'<span style="font-size: 12px; color: {TEXT if active else MUTED}">{text}</span>'
    width = f"width: {w}px; " if w else ("" if text else "width: 36px; ")
    return (f'<button aria-label="{label}" title="{label}" style="{width}height: 36px; padding: 0 {8 if text else 0}px; display: flex; align-items: center; '
            f'justify-content: center; gap: 6px; background: {bg}; border: {bd}; border-radius: 6px; cursor: pointer">{inner}</button>')


def sep(v=True):
    if v:
        return f'<div style="width: 1px; height: 24px; background: {LINE}; margin: 0 6px"></div>'
    return f'<div style="height: 1px; background: {LINE}"></div>'


def chip(text, color=MUTED, bg=BG2, mono=False):
    ff = f"font-family: {MONO}; " if mono else ""
    return (f'<span style="{ff}display: inline-flex; align-items: center; gap: 5px; height: 22px; padding: 0 8px; border-radius: 11px; '
            f'background: {bg}; color: {color}; font-size: 11.5px; white-space: nowrap">{text}</span>')


def dot(color, s=8):
    return f'<span style="width: {s}px; height: {s}px; border-radius: 50%; background: {color}; display: inline-block; flex-shrink: 0"></span>'


def btn(text, primary=False, ic=None, small=False, href=None):
    h = 30 if small else 36
    if primary:
        st = f"background: {ACC}; color: {ACC_INK}; border: 1px solid {ACC}; font-weight: 600"
        icol = ACC_INK
    else:
        st = f"background: {BG2}; color: {TEXT}; border: 1px solid {LINE}; font-weight: 500"
        icol = MUTED
    i = icon(ic, 16, icol) if ic else ""
    tag = "a" if href else "button"
    hr = f' href="{href}"' if href else ""
    dec = "text-decoration: none; box-sizing: border-box; " if href else ""
    return (f'<{tag}{hr} style="{dec}height: {h}px; padding: 0 {12 if small else 14}px; display: inline-flex; align-items: center; justify-content: center; gap: 7px; '
            f'border-radius: 6px; font-size: {12 if small else 13}px; cursor: pointer; {st}">{i}<span>{text}</span></{tag}>')


def seg(options, active, full=False, cols=None):
    items = []
    for o in options:
        on = o == active
        st = (f"background: {BG3}; color: {TEXT}; border: 1px solid {ACC}" if on
              else f"background: transparent; color: {MUTED}; border: 1px solid transparent")
        grow = "flex-grow: 1; " if full else ""
        if cols:
            grow = "padding: 0 4px; "
        items.append(f'<button aria-pressed="{"true" if on else "false"}" style="{grow}height: 28px; padding: 0 10px; border-radius: 5px; font-size: 12px; cursor: pointer; {st}">{o}</button>')
    lay = f"display: grid; grid-template-columns: repeat({cols}, minmax(0, 1fr))" if cols else "display: flex"
    return (f'<div role="group" style="{lay}; gap: 2px; padding: 2px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px">'
            + "".join(items) + '</div>')


_fid = [0]


def field(label, value, unit=None, mono=True, w=None, note=None):
    _fid[0] += 1
    fid = f"f{_fid[0]}"
    ff = f"font-family: {MONO}; " if mono else ""
    wd = f"width: {w}px; " if w else "flex-grow: 1; min-width: 0; "
    u = f'<span style="font-size: 11.5px; color: {DIM}">{unit}</span>' if unit else ""
    n = f'<span style="font-size: 11px; color: {DIM}">{note}</span>' if note else ""
    return (f'<div style="{wd}display: flex; flex-direction: column; gap: 4px">'
            f'<label for="{fid}" style="font-size: 11.5px; color: {MUTED}">{label}</label>'
            f'<div style="display: flex; align-items: center; gap: 6px; height: 30px; padding: 0 8px; background: {BG0}; border: 1px solid {LINE}; border-radius: 5px">'
            f'<input id="{fid}" value="{value}" style="{ff}flex-grow: 1; min-width: 0; width: 20px; background: transparent; border: 0; outline: none; font-size: 12.5px; color: {TEXT}">{u}</div>{n}</div>')


def select(label, value, w=None):
    _fid[0] += 1
    fid = f"f{_fid[0]}"
    wd = f"width: {w}px; " if w else "flex-grow: 1; min-width: 0; "
    return (f'<div style="{wd}display: flex; flex-direction: column; gap: 4px">'
            f'<label for="{fid}" style="font-size: 11.5px; color: {MUTED}">{label}</label>'
            f'<button id="{fid}" style="height: 32px; padding: 0 8px; display: flex; align-items: center; justify-content: space-between; gap: 6px; background: {BG0}; border: 1px solid {LINE}; border-radius: 5px; font-size: 12.5px; color: {TEXT}; cursor: pointer; text-align: left">'
            f'<span style="overflow: hidden; white-space: nowrap; text-overflow: ellipsis">{value}</span>{icon("chev", 14, DIM)}</button></div>')


def check(label, on=True):
    box = (f'<span style="width: 16px; height: 16px; border-radius: 4px; display: inline-flex; align-items: center; justify-content: center; flex-shrink: 0; '
           + (f'background: {ACC}">{icon("check", 12, ACC_INK, 3)}' if on else f'border: 1.5px solid {DIM}">') + '</span>')
    return (f'<button role="checkbox" aria-checked="{"true" if on else "false"}" style="display: flex; align-items: center; gap: 8px; padding: 0; background: transparent; border: 0; cursor: pointer; font-size: 12.5px; color: {TEXT}; text-align: left">{box}<span>{label}</span></button>')


def toggle(label, on=True):
    track = ACC if on else BG3
    pos = "margin-left: 14px" if on else "margin-left: 2px"
    kn = ACC_INK if on else MUTED
    return (f'<button role="switch" aria-checked="{"true" if on else "false"}" style="display: flex; align-items: center; justify-content: space-between; gap: 10px; width: 100%; padding: 0; background: transparent; border: 0; cursor: pointer; font-size: 12.5px; color: {TEXT}">'
            f'<span>{label}</span><span style="width: 30px; height: 18px; border-radius: 9px; background: {track}; display: flex; align-items: center; flex-shrink: 0">'
            f'<span style="width: 14px; height: 14px; border-radius: 50%; background: {kn}; {pos}"></span></span></button>')


def section(title, body, right="", pad=14, gap=10):
    return (f'<div style="display: flex; flex-direction: column; gap: {gap}px; padding: {pad}px; border-bottom: 1px solid {LINE}">'
            f'<div style="display: flex; align-items: center; justify-content: space-between; gap: 8px">'
            f'<h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">{title}</h3>{right}</div>'
            f'{body}</div>')


def kv(k, v, mono=True, vcol=TEXT):
    ff = f"font-family: {MONO}; " if mono else ""
    return (f'<div style="display: flex; justify-content: space-between; gap: 10px; font-size: 12.5px">'
            f'<span style="color: {MUTED}">{k}</span><span style="{ff}color: {vcol}; text-align: right">{v}</span></div>')


def cite(text):
    return (f'<div style="display: flex; gap: 6px; align-items: flex-start; font-size: 11.5px; color: {MUTED}; line-height: 1.45">'
            f'{icon("link", 13, DIM)}<span>{text}</span></div>')


def row(*items, gap=8, align="center", extra=""):
    return f'<div style="display: flex; align-items: {align}; gap: {gap}px; {extra}">' + "".join(items) + '</div>'


def col(*items, gap=8, extra=""):
    return f'<div style="display: flex; flex-direction: column; gap: {gap}px; {extra}">' + "".join(items) + '</div>'


# ---- app chrome ------------------------------------------------------------
MODULE_LINK = {"Studio": "Main.dc.html", "Grow": "Grow.dc.html", "Pack": "Pack.dc.html", "Relax": "Relax.dc.html", "Dynamics": "Dynamics.dc.html",
               "Equilibrate": "Convergence.dc.html", "React": "React.dc.html", "Analyze": "Analyze.dc.html", "Field": "ForceField.dc.html",
               "Jobs": "Jobs.dc.html", "Bench": "Bench.dc.html"}
MODULES = [("Studio", "hex"), ("Grow", "grow"), ("Pack", "pack"), ("Relax", "relax"), ("Dynamics", "dyn"),
           ("Equilibrate", "equil"), ("React", "react"), ("Analyze", "chart"), ("Field", "tag"), ("Jobs", "jobs"),
           ("Bench", "bench")]


def topbar(tabs, active_tab, crumb=None):
    t = []
    for name in tabs:
        on = name == active_tab
        st = (f"background: {BG0}; color: {TEXT}; border: 1px solid {LINE}; border-bottom-color: {BG0}" if on
              else f"background: transparent; color: {MUTED}; border: 1px solid transparent")
        t.append(f'<button style="height: 34px; margin-top: 10px; padding: 0 12px; display: flex; align-items: center; gap: 8px; border-radius: 7px 7px 0 0; font-size: 12.5px; cursor: pointer; {st}">'
                 f'{icon("cube", 14, ACC if on else DIM)}<span>{name}</span>{icon("close", 12, DIM) if on else ""}</button>')
    logo = (f'<a href="Start.dc.html" aria-label="CAPS start" style="width: 72px; display: flex; align-items: center; justify-content: center; flex-shrink: 0">'
            f'<svg width="30" height="30" viewBox="0 0 30 30" aria-label="CAPS"><path d="M15 2l11.3 6.5v13L15 28 3.7 21.5v-13z" fill="none" stroke="{ACC}" stroke-width="2"></path>'
            f'<circle cx="15" cy="15" r="3.4" fill="{ACC}"></circle><circle cx="15" cy="7.5" r="2" fill="{TEXT}"></circle><circle cx="21.5" cy="18.8" r="2" fill="{TEXT}"></circle><circle cx="8.5" cy="18.8" r="2" fill="{TEXT}"></circle>'
            f'<path d="M15 9.5v2.2M19.8 17.8l-1.9-1.1M10.2 17.8l1.9-1.1" stroke="{TEXT}" stroke-width="1.4"></path></svg></a>')
    right = (f'<div style="margin-left: auto; display: flex; align-items: center; gap: 10px; padding-right: 14px">'
             f'<button style="width: 300px; height: 30px; display: flex; align-items: center; gap: 8px; padding: 0 10px; background: {BG2}; border: 1px solid {LINE}; border-radius: 6px; color: {DIM}; font-size: 12.5px; cursor: pointer">'
             f'{icon("search", 14, DIM)}<span style="flex-grow: 1; text-align: left">Run a command or find anything</span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}; border: 1px solid {LINE}; border-radius: 4px; padding: 1px 5px">⌘K</span></button>'
             f'<button style="height: 30px; display: flex; align-items: center; gap: 8px; padding: 0 10px; background: transparent; border: 1px solid {LINE}; border-radius: 6px; font-size: 12px; color: {MUTED}; cursor: pointer">{dot(OK)}<span>Local · 16 threads · RTX&#160;4070</span></button>'
             f'</div>')
    cr = ""
    if crumb:
        cr = f'<div style="display: flex; align-items: center; gap: 6px; margin-left: 4px; margin-right: 18px; font-size: 12.5px; color: {MUTED}">{crumb}</div>'
    return (f'<header style="height: 44px; flex-shrink: 0; display: flex; align-items: stretch; background: {BG1}; border-bottom: 1px solid {LINE}">'
            f'{logo}{cr}<nav aria-label="Documents" style="display: flex; gap: 2px; align-items: flex-end">{"".join(t)}'
            f'<button aria-label="New document" style="width: 30px; height: 30px; margin: 0 0 3px 4px; display: flex; align-items: center; justify-content: center; background: transparent; border: 0; cursor: pointer">{icon("plus", 15, DIM)}</button></nav>'
            f'{right}</header>')


def rail(active):
    items = []
    for name, ic in MODULES:
        on = name == active
        st = f"background: {BG3}; color: {TEXT}" if on else f"background: transparent; color: {DIM}"
        bar = f'<span style="position: absolute; left: 0; top: 10px; width: 3px; height: 30px; border-radius: 0 2px 2px 0; background: {ACC}"></span>' if on else ""
        items.append(f'<a href="{MODULE_LINK[name]}" aria-current="{"page" if on else "false"}" style="position: relative; width: 60px; height: 52px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 3px; border-radius: 8px; text-decoration: none; font-size: 10.5px; {st}">'
                     f'{bar}{icon(ic, 19, ACC if on else MUTED)}<span>{name}</span></a>')
    return (f'<nav aria-label="Modules" style="width: 72px; flex-shrink: 0; display: flex; flex-direction: column; align-items: center; gap: 2px; padding-top: 8px; background: {BG1}; border-right: 1px solid {LINE}">'
            + "".join(items) +
            f'<div style="flex-grow: 1"></div><a href="Settings.dc.html" aria-label="Settings" style="width: 60px; height: 44px; margin-bottom: 8px; display: flex; align-items: center; justify-content: center; border-radius: 8px">{icon("gear", 19, MUTED)}</a></nav>')


def statusbar(left, right):
    return (f'<footer style="height: 26px; flex-shrink: 0; display: flex; align-items: center; justify-content: space-between; padding: 0 14px; background: {BG1}; border-top: 1px solid {LINE}; font-family: {MONO}; font-size: 11px; color: {MUTED}">'
            f'<div style="display: flex; gap: 18px; align-items: center">{left}</div><div style="display: flex; gap: 18px; align-items: center">{right}</div></footer>')


def panel_head(title, right=""):
    return (f'<div style="height: 38px; flex-shrink: 0; display: flex; align-items: center; justify-content: space-between; padding: 0 14px; border-bottom: 1px solid {LINE}">'
            f'<h2 style="margin: 0; font-size: 13px; font-weight: 600">{title}</h2>{right}</div>')


def tabs(names, active, size=12.5):
    out = []
    for n in names:
        on = n == active
        st = f"color: {TEXT}; border-bottom: 2px solid {ACC}" if on else f"color: {MUTED}; border-bottom: 2px solid transparent"
        out.append(f'<button role="tab" aria-selected="{"true" if on else "false"}" style="height: 38px; padding: 0 12px; background: transparent; border: 0; font-size: {size}px; cursor: pointer; {st}">{n}</button>')
    return f'<div role="tablist" style="display: flex; gap: 2px; padding: 0 6px; border-bottom: 1px solid {LINE}; flex-shrink: 0">' + "".join(out) + '</div>'


def table(headers, rows, widths=None, mono_cols=(), align_right=(), fs=12, rowh=30, hl=None):
    th = []
    for i, h in enumerate(headers):
        w = f"width: {widths[i]}; " if widths else ""
        ta = "right" if i in align_right else "left"
        th.append(f'<th scope="col" style="{w}text-align: {ta}; padding: 0 10px; height: 30px; font-size: 11px; font-weight: 600; letter-spacing: 0.04em; text-transform: uppercase; color: {DIM}; border-bottom: 1px solid {LINE}; white-space: nowrap">{h}</th>')
    tr = []
    for ri, r in enumerate(rows):
        bg = f"background: {BG3}; " if hl is not None and ri in (hl if isinstance(hl, (set, list, tuple)) else {hl}) else ""
        tds = []
        for i, c in enumerate(r):
            ff = f"font-family: {MONO}; " if i in mono_cols else ""
            ta = "right" if i in align_right else "left"
            tds.append(f'<td style="{ff}text-align: {ta}; padding: 0 10px; height: {rowh}px; font-size: {fs}px; border-bottom: 1px solid {BG2}; white-space: nowrap">{c}</td>')
        tr.append(f'<tr style="{bg}">' + "".join(tds) + '</tr>')
    return (f'<table style="width: 100%; border-collapse: collapse; table-layout: fixed">'
            f'<thead><tr>{"".join(th)}</tr></thead><tbody>{"".join(tr)}</tbody></table>')


# ---- 3D scene rendering -----------------------------------------------------
def mix(c1, c2, t):
    a = [int(c1[i:i + 2], 16) for i in (1, 3, 5)]
    b = [int(c2[i:i + 2], 16) for i in (1, 3, 5)]
    return "#%02X%02X%02X" % tuple(round(a[k] * (1 - t) + b[k] * t) for k in range(3))


def rot(p, yaw, pitch, roll=0.0):
    x, y, z = p
    cy, sy = math.cos(yaw), math.sin(yaw)
    x, z = x * cy + z * sy, -x * sy + z * cy
    cp, sp = math.cos(pitch), math.sin(pitch)
    y, z = y * cp - z * sp, y * sp + z * cp
    cr, sr = math.cos(roll), math.sin(roll)
    x, y = x * cr - y * sr, x * sr + y * cr
    return (x, y, z)


def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def sub(a, b): return (a[0] - b[0], a[1] - b[1], a[2] - b[2])
def add(a, b): return (a[0] + b[0], a[1] + b[1], a[2] + b[2])
def mul(a, s): return (a[0] * s, a[1] * s, a[2] * s)
def dot3(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def cross(a, b): return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def norm(a):
    l = math.sqrt(dot3(a, a)) or 1.0
    return (a[0] / l, a[1] / l, a[2] / l)


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = sub(p0, p1), sub(p2, p1), sub(p3, p2)
    b1n = norm(b1)
    v = sub(b0, mul(b1n, dot3(b0, b1n)))
    w = sub(b2, mul(b1n, dot3(b2, b1n)))
    x = dot3(v, w)
    y = dot3(cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def angle(a, b, c):
    u, v = norm(sub(a, b)), norm(sub(c, b))
    return math.degrees(math.acos(max(-1, min(1, dot3(u, v)))))


def perceive_bonds(atoms, tol=0.45):
    cov = {"C": 0.76, "H": 0.31, "O": 0.66, "N": 0.71, "Si": 1.11, "Na": 1.66, "Cl": 1.02, "S": 1.05}
    bonds = []
    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            ei, ej = atoms[i]["e"], atoms[j]["e"]
            if ei in ("Na", "Cl") or ej in ("Na", "Cl"):
                continue
            if dist(atoms[i]["p"], atoms[j]["p"]) < cov[ei] + cov[ej] + tol - 0.25:
                bonds.append((i, j))
    return bonds


class Scene:
    """Painter's-algorithm ball-and-stick renderer producing inline SVG."""

    def __init__(self, sid, w, h, yaw=0.5, pitch=0.35, roll=0.0, scale=None, center=None, persp=0.0,
                 fog=0.55, outline=True, bond_w=0.2, atom_k=1.0, bg=BG0):
        self.sid, self.w, self.h = sid, w, h
        self.yaw, self.pitch, self.roll = yaw, pitch, roll
        self.scale, self.center, self.persp = scale, center, persp
        self.fog, self.outline, self.bond_w, self.atom_k, self.bg = fog, outline, bond_w, atom_k, bg
        self.atoms, self.bonds, self.lines, self.overlays, self.tubes = [], [], [], [], []
        self.halo = set()

    def add_atoms(self, atoms, bonds=None):
        off = len(self.atoms)
        self.atoms.extend(atoms)
        if bonds is None:
            bonds = perceive_bonds(atoms)
        self.bonds.extend((i + off, j + off) for i, j in bonds)
        return off

    def add_box(self, origin, a, b, c, color=MUTED, width=1.2, dash=None, opacity=0.9):
        o = origin
        pts = [o, add(o, a), add(o, b), add(o, c), add(add(o, a), b), add(add(o, a), c), add(add(o, b), c), add(add(add(o, a), b), c)]
        edges = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]
        for i, j in edges:
            self.lines.append((pts[i], pts[j], color, width, dash, opacity))

    def add_line(self, p, q, color, width=1.0, dash=None, opacity=1.0):
        self.lines.append((p, q, color, width, dash, opacity))

    def add_tube(self, pts, color, width):
        self.tubes.append((pts, color, width))

    def _prep(self):
        allp = [a["p"] for a in self.atoms] + [l[0] for l in self.lines] + [l[1] for l in self.lines] + [p for t in self.tubes for p in t[0]]
        if self.center is None:
            n = len(allp)
            self.center = (sum(p[0] for p in allp) / n, sum(p[1] for p in allp) / n, sum(p[2] for p in allp) / n)
        rp = [rot(sub(p, self.center), self.yaw, self.pitch, self.roll) for p in allp]
        if self.scale is None:
            ext = max(max(abs(p[0]) for p in rp) / (self.w * 0.44), max(abs(p[1]) for p in rp) / (self.h * 0.42))
            self.scale = 1.0 / ext
        zs = [p[2] for p in rp]
        self.zmin, self.zmax = min(zs), max(zs)

    def P(self, p):
        r = rot(sub(p, self.center), self.yaw, self.pitch, self.roll)
        k = 1.0 + self.persp * r[2] / max(1e-6, (self.zmax - self.zmin))
        return (self.w / 2 + r[0] * self.scale * k, self.h / 2 - r[1] * self.scale * k, r[2], k)

    def fogt(self, z):
        t = (z - self.zmin) / max(1e-6, self.zmax - self.zmin)
        return round((1 - t) * self.fog * 4) / 4  # quantised fog level 0..fog

    def svg(self, extra_defs="", overlay=""):
        self._prep()
        defs, used = [], set()
        items = []
        proj = [self.P(a["p"]) for a in self.atoms]

        def grad(col, f):
            gid = f"{self.sid}g{col[1:]}{int(f * 100)}"
            if gid not in used:
                used.add(gid)
                base = mix(col, self.bg, f)
                defs.append(f'<radialGradient id="{gid}" cx="0.36" cy="0.32" r="0.72">'
                            f'<stop offset="0" stop-color="{mix(base, "#FFFFFF", 0.55 * (1 - f))}"></stop>'
                            f'<stop offset="0.45" stop-color="{base}"></stop>'
                            f'<stop offset="1" stop-color="{mix(base, "#000000", 0.55)}"></stop></radialGradient>')
            return gid

        for i, j in self.bonds:
            pi, pj = proj[i], proj[j]
            mx, my = (pi[0] + pj[0]) / 2, (pi[1] + pj[1]) / 2
            bw = max(1.2, self.bond_w * self.scale * (pi[3] + pj[3]) / 2)
            for (pa, ai) in ((pi, i), (pj, j)):
                col = self.atoms[ai].get("c") or ELEM[self.atoms[ai]["e"]][0]
                f = self.fogt(pa[2])
                c2 = mix(col, self.bg, f)
                z = (pa[2] + (pi[2] + pj[2]) / 2) / 2 - 0.3
                s = ""
                if self.outline:
                    s += f'<line x1="{pa[0]:.1f}" y1="{pa[1]:.1f}" x2="{mx:.1f}" y2="{my:.1f}" stroke="#070809" stroke-width="{bw + 1.6:.1f}" stroke-linecap="round"></line>'
                s += f'<line x1="{pa[0]:.1f}" y1="{pa[1]:.1f}" x2="{mx:.1f}" y2="{my:.1f}" stroke="{c2}" stroke-width="{bw:.1f}" stroke-linecap="round"></line>'
                items.append((z, s))
        for idx, a in enumerate(self.atoms):
            px, py, z, k = proj[idx]
            col = a.get("c") or ELEM[a["e"]][0]
            r = a.get("r", ELEM[a["e"]][1]) * self.atom_k * self.scale * k
            f = self.fogt(z)
            s = f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{r:.1f}" fill="url(#{grad(col, f)})"' + (f' stroke="#070809" stroke-width="0.9"' if self.outline else "") + '></circle>'
            if idx in self.halo:
                s += f'<circle cx="{px:.1f}" cy="{py:.1f}" r="{r + 3.5:.1f}" fill="none" stroke="{SEL}" stroke-width="2"></circle>'
            items.append((z, s))
        for (p, q, colr, wd, dash, op) in self.lines:
            a, b = self.P(p), self.P(q)
            d = f' stroke-dasharray="{dash}"' if dash else ""
            items.append(((a[2] + b[2]) / 2 - 0.5, f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" stroke="{colr}" stroke-width="{wd}" stroke-opacity="{op}"{d}></line>'))
        for pts, colr, wd in self.tubes:
            pp = [self.P(p) for p in pts]
            for k in range(len(pp) - 1):
                a, b = pp[k], pp[k + 1]
                z = (a[2] + b[2]) / 2
                f = self.fogt(z)
                c2 = mix(colr, self.bg, f)
                ww = wd * self.scale * (a[3] + b[3]) / 2
                s = (f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" stroke="#070809" stroke-width="{ww + 2:.1f}" stroke-linecap="round"></line>'
                     f'<line x1="{a[0]:.1f}" y1="{a[1]:.1f}" x2="{b[0]:.1f}" y2="{b[1]:.1f}" stroke="{c2}" stroke-width="{ww:.1f}" stroke-linecap="round"></line>'
                     f'<line x1="{a[0]:.1f}" y1="{a[1] - ww * 0.22:.1f}" x2="{b[0]:.1f}" y2="{b[1] - ww * 0.22:.1f}" stroke="{mix(c2, "#FFFFFF", 0.35)}" stroke-width="{ww * 0.25:.1f}" stroke-linecap="round" stroke-opacity="0.7"></line>')
                items.append((z, s))
        items.sort(key=lambda t: t[0])
        return (f'<svg width="{self.w}" height="{self.h}" viewBox="0 0 {self.w} {self.h}" role="img" aria-label="3D structure view" style="display: block">'
                f'<defs>{"".join(defs)}{extra_defs}</defs>' + "".join(s for _, s in items) + overlay + '</svg>')

    def screen(self, i):
        return self.P(self.atoms[i]["p"])


def label_pill(x, y, text, color=TEXT, bg="#0B0D0Fcc", border=LINE, fs=11, mono=True):
    w = len(text) * fs * 0.62 + 12
    ff = MONO.replace("'", "") if mono else "IBM Plex Sans"
    return (f'<g><rect x="{x:.1f}" y="{y - 9:.1f}" width="{w:.1f}" height="18" rx="4" fill="#0B0D0F" fill-opacity="0.85" stroke="{border}"></rect>'
            f'<text x="{x + 6:.1f}" y="{y + 4:.1f}" fill="{color}" font-size="{fs}" font-family="{ff}">{esc(text)}</text></g>')


def gizmo(x, y, yaw, pitch, L=26):
    out = []
    axes = [((1, 0, 0), "#FF7B72", "x"), ((0, 1, 0), "#7CC784", "y"), ((0, 0, 1), "#6FA8FF", "z")]
    proj = []
    for v, c, n in axes:
        r = rot(v, yaw, pitch)
        proj.append((r[2], x + r[0] * L, y - r[1] * L, c, n))
    proj.sort()
    for z, ex, ey, c, n in proj:
        out.append(f'<line x1="{x}" y1="{y}" x2="{ex:.1f}" y2="{ey:.1f}" stroke="{c}" stroke-width="2" stroke-linecap="round"></line>'
                   f'<circle cx="{ex:.1f}" cy="{ey:.1f}" r="7" fill="{BG1}" stroke="{c}" stroke-width="1.2"></circle>'
                   f'<text x="{ex:.1f}" y="{ey + 3.5:.1f}" text-anchor="middle" fill="{c}" font-size="10" font-family="IBM Plex Mono">{n}</text>')
    return "".join(out)


# ---- plots ---------------------------------------------------------------------
def plot(w, h, series, xr, yr, xticks, yticks, xlabel="", ylabel="", pad=(40, 12, 12, 28), fmt=None, bands=None, markers=None):
    L, R, T, B = pad
    pw, ph = w - L - R, h - T - B

    def X(v): return L + (v - xr[0]) / (xr[1] - xr[0]) * pw
    def Y(v): return T + ph - (v - yr[0]) / (yr[1] - yr[0]) * ph
    out = [f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-label="{esc(ylabel)} versus {esc(xlabel)}" style="display: block">']
    if bands:
        for (x0, x1, c, op) in bands:
            out.append(f'<rect x="{X(x0):.1f}" y="{T}" width="{X(x1) - X(x0):.1f}" height="{ph}" fill="{c}" fill-opacity="{op}"></rect>')
    for t in yticks:
        out.append(f'<line x1="{L}" y1="{Y(t):.1f}" x2="{L + pw}" y2="{Y(t):.1f}" stroke="{BG3}" stroke-width="1"></line>'
                   f'<text x="{L - 6}" y="{Y(t) + 3.5:.1f}" text-anchor="end" fill="{DIM}" font-size="10" font-family="IBM Plex Mono">{t:g}</text>')
    for t in xticks:
        out.append(f'<text x="{X(t):.1f}" y="{T + ph + 15}" text-anchor="middle" fill="{DIM}" font-size="10" font-family="IBM Plex Mono">{t:g}</text>')
    out.append(f'<line x1="{L}" y1="{T + ph}" x2="{L + pw}" y2="{T + ph}" stroke="{LINE}"></line>')
    for pts, c, wd, dash in series:
        d = " ".join(f"{X(a):.1f},{Y(b):.1f}" for a, b in pts)
        da = f' stroke-dasharray="{dash}"' if dash else ""
        out.append(f'<polyline points="{d}" fill="none" stroke="{c}" stroke-width="{wd}" stroke-linejoin="round" stroke-linecap="round"{da}></polyline>')
    if markers:
        for (a, b, c) in markers:
            out.append(f'<circle cx="{X(a):.1f}" cy="{Y(b):.1f}" r="3" fill="{c}"></circle>')
    if xlabel:
        out.append(f'<text x="{L + pw}" y="{h - 2}" text-anchor="end" fill="{DIM}" font-size="10" font-family="IBM Plex Sans">{esc(xlabel)}</text>')
    if ylabel:
        out.append(f'<text x="{L}" y="{T - 2 if T > 10 else 9}" fill="{DIM}" font-size="10" font-family="IBM Plex Sans">{esc(ylabel)}</text>')
    out.append('</svg>')
    return "".join(out)
