import math
from lib import *
from screen_row7 import LIGHT
from screen_studio2 import spacer


def lum(hexc):
    r, g, b = (int(hexc[i:i + 2], 16) / 255 for i in (1, 3, 5))
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def contrast(a, b):
    la, lb = lum(a), lum(b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def heading(t, sub_=""):
    s = f'<p style="margin: 0; font-size: 13px; color: {MUTED}">{sub_}</p>' if sub_ else ""
    return f'<div style="display: flex; flex-direction: column; gap: 4px"><h2 style="margin: 0; font-size: 15px; font-weight: 600">{t}</h2>{s}</div>'


def logo_mark(size, fg, acc, bg=None, rx=None):
    s = size / 30
    back = f'<rect width="{size}" height="{size}" rx="{rx or size * 0.22}" fill="{bg}"></rect>' if bg else ""
    pad = size * 0.16 if bg else 0
    k = (size - 2 * pad) / 30
    g = (f'<g transform="translate({pad},{pad}) scale({k})"><path d="M15 2l11.3 6.5v13L15 28 3.7 21.5v-13z" fill="none" stroke="{acc}" stroke-width="2"></path>'
         f'<circle cx="15" cy="15" r="3.4" fill="{acc}"></circle><circle cx="15" cy="7.5" r="2" fill="{fg}"></circle><circle cx="21.5" cy="18.8" r="2" fill="{fg}"></circle><circle cx="8.5" cy="18.8" r="2" fill="{fg}"></circle>'
         f'<path d="M15 9.5v2.2M19.8 17.8l-1.9-1.1M10.2 17.8l1.9-1.1" stroke="{fg}" stroke-width="1.4"></path></g>')
    return f'<svg width="{size}" height="{size}" viewBox="0 0 {size} {size}" aria-hidden="true">{back}{g}</svg>'


# ------------------------------------------------------------------ Foundations
def foundations():
    lockups = (f'<div style="display: flex; gap: 14px">'
               f'<div style="flex-grow: 1; display: flex; align-items: center; gap: 16px; padding: 22px; background: {BG0}; border: 1px solid {LINE}; border-radius: 10px">{logo_mark(56, TEXT, ACC)}'
               f'<div style="display: flex; flex-direction: column; gap: 2px"><span style="font-size: 30px; font-weight: 700; letter-spacing: 0.04em">CAPS</span><span style="font-size: 12px; color: {MUTED}">Chain Assembly and Packing Suite</span></div></div>'
               f'<div style="flex-grow: 1; display: flex; align-items: center; gap: 16px; padding: 22px; background: #F4F2EE; border: 1px solid #D9D5CF; border-radius: 10px; color: #1E2226">{logo_mark(56, "#1E2226", "#A35C00")}'
               f'<div style="display: flex; flex-direction: column; gap: 2px"><span style="font-size: 30px; font-weight: 700; letter-spacing: 0.04em">CAPS</span><span style="font-size: 12px; color: #555B62">Chain Assembly and Packing Suite</span></div></div>'
               f'<div style="display: flex; align-items: flex-end; gap: 14px; padding: 18px 22px; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px">'
               f'{logo_mark(96, TEXT, ACC, "#1B1F23")}{logo_mark(48, TEXT, ACC, "#1B1F23")}{logo_mark(32, TEXT, ACC, "#1B1F23")}{logo_mark(16, TEXT, ACC, "#1B1F23", 3)}</div></div>')
    lock_note = f'<div style="font-size: 12px; color: {MUTED}">The mark is a hexagonal ring (a repeat unit) around a branch point with three chain ends. On light backgrounds the amber darkens to #A35C00 so it passes contrast as text.</div>'
    tokens = [("bg/base", BG0, "viewport, app ground"), ("bg/panel", BG1, "panels, cards"), ("bg/raised", BG2, "inputs, chips"), ("bg/selected", BG3, "selected rows"),
              ("line", LINE, "dividers, borders"), ("text", TEXT, "primary text"), ("text/muted", MUTED, "labels"), ("text/dim", DIM, "captions, axes"),
              ("accent", ACC, "primary action, focus"), ("select", SEL, "selection, links in 3D"), ("ok", OK, "success"), ("error", ERR, "errors, clashes")]
    rows = []
    for name, c, use in tokens:
        light = LIGHT.get(c.upper(), c)
        on_dark = contrast(c, BG1) if name.startswith(("text", "accent", "select", "ok", "error")) else None
        on_light = contrast(light, "#FFFFFF") if name.startswith(("text", "accent", "select", "ok", "error")) else None
        fmt = lambda v: "—" if v is None else f'<span style="color: {OK if v >= 4.5 else (WARN if v >= 3 else ERR)}">{v:.1f}:1</span>'
        rows.append([f'<span style="font-family: {MONO}">{name}</span>',
                     f'<span style="display: inline-flex; align-items: center; gap: 8px"><span style="width: 18px; height: 18px; border-radius: 4px; background: {c}; border: 1px solid {LINE}"></span><span style="font-family: {MONO}">{c}</span></span>',
                     f'<span style="display: inline-flex; align-items: center; gap: 8px"><span style="width: 18px; height: 18px; border-radius: 4px; background: {light}; border: 1px solid {LINE}"></span><span style="font-family: {MONO}">{light}</span></span>',
                     fmt(on_dark), fmt(on_light), f'<span style="color: {MUTED}">{use}</span>'])
    tt = table(["Token", "Graphite (dark)", "Paper (light)", "vs panel", "vs white", "Use"], rows, ["14%", "19%", "19%", "10%", "10%", "28%"], align_right=(3, 4), fs=12, rowh=30)
    scale = [("Display", 28, 600, "How CAPS is organised"), ("Title", 19, 600, "Amorphous cell · CAPS Grow"), ("Heading", 15, 600, "Crystal builder"),
             ("Body", 13, 400, "Grow chains into a periodic cell, then Pack."), ("Label", 11, 600, "SECTION LABEL"), ("Mono", 12.5, 400, "C=Cc1ccccc1 · 1.041 g/cm³")]
    ts = "".join(f'<div style="display: flex; align-items: baseline; gap: 16px; height: {max(30, s + 14)}px; border-bottom: 1px solid {BG2}">'
                 f'<span style="width: 70px; flex-shrink: 0; font-size: 11px; color: {DIM}">{n}</span><span style="width: 70px; flex-shrink: 0; font-family: {MONO}; font-size: 11px; color: {DIM}">{s}/{w}</span>'
                 f'<span style="font-size: {s}px; font-weight: {w}; {"font-family: " + MONO + "; " if n == "Mono" else ""}{"letter-spacing: 0.08em; " if n == "Label" else ""}white-space: nowrap; overflow: hidden; text-overflow: ellipsis">{t}</span></div>' for n, s, w, t in scale)
    sp = "".join(f'<div style="display: flex; flex-direction: column; align-items: center; gap: 6px"><span style="width: {v}px; height: {v}px; background: {ACC}; border-radius: 2px"></span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{v}</span></div>' for v in (2, 4, 6, 8, 10, 12, 14, 18, 22, 28))
    rad = "".join(f'<div style="display: flex; flex-direction: column; align-items: center; gap: 6px"><span style="width: 44px; height: 44px; border: 1.5px solid {MUTED}; border-radius: {v}px"></span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{v}</span></div>' for v in (4, 6, 8, 10, 12))
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 22px; padding: 32px 40px">'
            + f'<div style="display: flex; align-items: flex-end; gap: 20px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">Foundations</h1><span style="font-size: 13px; color: {MUTED}">IBM Plex Sans + IBM Plex Mono · two themes · WCAG 2.2 contrast computed</span></div>'
            + col(lockups, lock_note, gap=8)
            + '<div style="display: flex; gap: 28px; flex-grow: 1; min-height: 0">'
            + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 10px">{heading("Colour tokens")}{tt}</div>'
            + f'<div style="width: 470px; flex-shrink: 0; display: flex; flex-direction: column; gap: 18px">'
            + col(heading("Type scale"), f'<div style="display: flex; flex-direction: column">{ts}</div>', gap=8)
            + col(heading("Spacing (px)"), f'<div style="display: flex; align-items: flex-end; gap: 12px">{sp}</div>', gap=10)
            + col(heading("Radii (px)"), f'<div style="display: flex; gap: 14px">{rad}</div>', gap=10)
            + '</div></div></div>')
    return page("CAPS — foundations", body)


# ------------------------------------------------------------------ Components
def components():
    def state_row(label, *items):
        return (f'<div style="display: flex; align-items: center; gap: 14px; min-height: 44px; border-bottom: 1px solid {BG2}"><span style="width: 120px; flex-shrink: 0; font-size: 12px; color: {DIM}">{label}</span>'
                + "".join(items) + '</div>')
    hover = f'<button style="height: 36px; padding: 0 14px; display: inline-flex; align-items: center; gap: 7px; border-radius: 6px; font-size: 13px; font-weight: 600; background: #FFBC55; color: {ACC_INK}; border: 1px solid #FFBC55; cursor: pointer">{icon("play", 16, ACC_INK)}Hover</button>'
    focus = f'<button style="height: 36px; padding: 0 14px; display: inline-flex; align-items: center; gap: 7px; border-radius: 6px; font-size: 13px; font-weight: 600; background: {ACC}; color: {ACC_INK}; border: 1px solid {ACC}; outline: 2px solid {SEL}; outline-offset: 2px; cursor: pointer">{icon("play", 16, ACC_INK)}Focus</button>'
    disabled = f'<button disabled style="height: 36px; padding: 0 14px; display: inline-flex; align-items: center; gap: 7px; border-radius: 6px; font-size: 13px; font-weight: 600; background: {BG3}; color: {DIM}; border: 1px solid {LINE}">{icon("play", 16, DIM)}Disabled</button>'
    busy = f'<button style="height: 36px; padding: 0 14px; display: inline-flex; align-items: center; gap: 9px; border-radius: 6px; font-size: 13px; font-weight: 600; background: {ACC}; color: {ACC_INK}; border: 1px solid {ACC}; cursor: progress"><span style="width: 14px; height: 14px; border-radius: 50%; border: 2px solid {ACC_INK}; border-right-color: transparent"></span>Running</button>'
    buttons = col(state_row("Primary", btn("Default", True, "play"), hover, focus, busy, disabled),
                  state_row("Secondary", btn("Default", ic="download"), btn("Small", ic="copy", small=True), btn("Text only")),
                  state_row("Tool (icon)", tbtn("cursor", "Select", True), tbtn("rotate", "Rotate"), tbtn("ruler", "Measure"), tbtn("wand", "Auto-clean", True, "With label")),
                  gap=0)
    err_field = (f'<div style="width: 200px; display: flex; flex-direction: column; gap: 4px"><label for="ef" style="font-size: 11.5px; color: {MUTED}">Density</label>'
                 f'<div style="display: flex; align-items: center; gap: 6px; height: 30px; padding: 0 8px; background: {BG0}; border: 1px solid {ERR}; border-radius: 5px"><input id="ef" value="-1.04" style="font-family: {MONO}; flex-grow: 1; min-width: 0; width: 20px; background: transparent; border: 0; outline: none; font-size: 12.5px; color: {TEXT}"><span style="font-size: 11.5px; color: {DIM}">g/cm³</span></div>'
                 f'<span style="font-size: 11px; color: {ERR}">Must be greater than 0</span></div>')
    inputs = col(state_row("Field", field("Target density", "1.04", "g/cm³", w=200), field("Tolerance", "2.0", "Å", w=140), err_field),
                 state_row("Select", select("Force field", "GAFF2 2.11", w=200), select("Thermostat", "CSVR (Bussi 2007)", w=220)),
                 state_row("Choice", seg(["Iso", "Syndio", "Atactic"], "Atactic"), f'<div style="width: 180px">{toggle("Auto-clean", True)}</div>', check("Include hydrogens", True), check("Wrap", False)),
                 gap=0)
    fb = col(state_row("Chips", chip("neutral"), chip(f"{dot(OK)} done", OK, "#16261A"), chip(f"{dot(ACC)} running", ACC, "#3A2C14"), chip(f"{dot(ERR)} failed", ERR, "#2A1414"), chip("C=Cc1ccccc1", MUTED, BG2, True)),
             state_row("Progress", f'<div style="width: 320px; height: 6px; border-radius: 3px; background: {BG3}"><div style="width: 57%; height: 6px; border-radius: 3px; background: {ACC}"></div></div>', f'<span style="font-family: {MONO}; font-size: 12px">57 %</span>'),
             state_row("Inline issue", f'<div style="display: flex; align-items: center; gap: 8px; font-size: 12px">{icon("alert", 15, WARN)}<span>C31–H44 1.21 Å is longer than usual</span><button style="height: 24px; padding: 0 8px; background: transparent; border: 1px solid {LINE}; border-radius: 5px; font-size: 11.5px; color: {ACC}; cursor: pointer">Clean</button></div>'),
             state_row("Citation", f'<div style="width: 460px">{cite("Bussi, Donadio &amp; Parrinello, <i>J. Chem. Phys.</i> 126, 014101 (2007)")}</div>'),
             gap=0)
    toast = (f'<div style="width: 400px; box-sizing: border-box; display: flex; gap: 12px; align-items: flex-start; padding: 14px; background: {BG2}; border: 1px solid {LINE}; border-radius: 10px; box-shadow: 0 12px 32px #00000080">{icon("check", 18, OK)}'
             f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 4px"><span style="font-size: 13px; font-weight: 600">Growth finished</span><span style="font-size: 12px; color: {MUTED}">grow-3 · 20 chains · 0 close contacts</span>'
             f'<div style="display: flex; gap: 8px; margin-top: 4px">{btn("Open cell", True, small=True)}{btn("Queue Relax", small=True)}</div></div>{tbtn("close", "Dismiss")}</div>')
    empty = (f'<div style="width: 400px; box-sizing: border-box; display: flex; flex-direction: column; align-items: center; gap: 10px; padding: 26px; background: {BG1}; border: 1px dashed {LINE}; border-radius: 10px; text-align: center">'
             f'{icon("jobs", 28, DIM, 1.4)}<span style="font-size: 13.5px; font-weight: 600">No jobs yet</span><span style="font-size: 12px; color: {MUTED}">Long operations such as Grow, Pack and Dynamics appear here with progress, logs and checkpoints.</span>{btn("Start a Grow job", True, "grow", small=True)}</div>')
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 20px; padding: 32px 40px">'
            + f'<div style="display: flex; align-items: flex-end; gap: 20px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">Components</h1><span style="font-size: 13px; color: {MUTED}">Real &lt;button&gt;, &lt;input&gt; and &lt;label&gt; elements · focus ring 2 px select-cyan · 36 px default, 44 px touch</span></div>'
            + '<div style="display: flex; gap: 32px">'
            + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 22px">'
            + col(heading("Buttons"), buttons, gap=8) + col(heading("Inputs"), inputs, gap=8) + col(heading("Feedback"), fb, gap=8) + '</div>'
            + f'<div style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; gap: 22px">'
            + col(heading("Toast"), toast, gap=8) + col(heading("Empty state"), empty, gap=8)
            + col(heading("Panel section"), f'<div style="border: 1px solid {LINE}; border-radius: 10px; overflow: hidden; background: {BG1}">{section("Section label", col(kv("Key", "value"), kv("Formula", "C<sub>8</sub>H<sub>8</sub>"), gap=7), chip("chip"))}</div>', gap=8)
            + '</div></div></div>')
    return page("CAPS — components", body)


# ------------------------------------------------------------------ Icons
def icons():
    names = list(ICONS.keys())
    cells = "".join(f'<div style="display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 8px; height: 92px; background: {BG1}; border: 1px solid {LINE}; border-radius: 8px">'
                    f'{icon(n, 24, TEXT, 1.6)}<span style="font-family: {MONO}; font-size: 10.5px; color: {DIM}">{n}</span></div>' for n in names)
    sizes = "".join(f'<div style="display: flex; flex-direction: column; align-items: center; gap: 6px">{icon("hex", s, ACC, 1.7 if s > 16 else 2)}<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{s}</span></div>' for s in (14, 16, 18, 20, 24, 32))
    rules = [("Grid", "24 × 24, 2 px padding, drawn on whole and half pixels"), ("Stroke", "1.7 px at 16–20 px, 1.4 px at 24 px+, round caps and joins"),
             ("Colour", "one colour from the text or accent tokens; never multi-colour"), ("Meaning", "every icon-only button carries an aria-label and a tooltip"),
             ("Science icons", "ring = aromatic ring, hex = fragment, grow / pack / relax / dyn / equil = engine modules")]
    rl = "".join(f'<div style="display: flex; gap: 12px; padding: 8px 0; border-bottom: 1px solid {BG2}; font-size: 12.5px"><span style="width: 110px; flex-shrink: 0; color: {DIM}">{a}</span><span>{b}</span></div>' for a, b in rules)
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 20px; padding: 32px 40px">'
            + f'<div style="display: flex; align-items: flex-end; gap: 20px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">Icons</h1><span style="font-size: 13px; color: {MUTED}">{len(names)} original stroke icons used across every screen</span></div>'
            + '<div style="display: flex; gap: 28px">'
            + f'<div style="flex-grow: 1; min-width: 0; display: grid; grid-template-columns: repeat(10, minmax(0, 1fr)); gap: 8px; align-content: start">{cells}</div>'
            + f'<div style="width: 360px; flex-shrink: 0; display: flex; flex-direction: column; gap: 18px">'
            + col(heading("Sizes"), f'<div style="display: flex; align-items: flex-end; gap: 18px">{sizes}</div>', gap=10)
            + col(heading("Rules"), f'<div>{rl}</div>', gap=6) + '</div></div></div>')
    return page("CAPS — icons", body)


# ------------------------------------------------------------------ User flows (clickable)
def flows():
    def node(board, title, sub_, color=LINE):
        return (f'<a href="{board}" style="display: flex; flex-direction: column; gap: 3px; width: 138px; padding: 12px 12px; background: {BG1}; border: 1px solid {color}; border-radius: 9px; text-decoration: none; color: {TEXT}">'
                f'<span style="font-size: 13px; font-weight: 600">{title}</span><span style="font-size: 11.5px; color: {MUTED}; line-height: 1.4">{sub_}</span>'
                f'<span style="font-family: {MONO}; font-size: 10.5px; color: {SEL}">open ›</span></a>')
    arr = f'<span style="display: flex; align-items: center; justify-content: center; width: 28px; flex-shrink: 0">{icon("chevr", 18, DIM)}</span>'

    def flow(title, desc, steps):
        return (f'<section style="display: flex; flex-direction: column; gap: 10px; padding: 18px; background: {BG0}; border: 1px solid {LINE}; border-radius: 12px">'
                f'<div style="display: flex; align-items: baseline; gap: 12px"><h2 style="margin: 0; font-size: 15px; font-weight: 600">{title}</h2><span style="font-size: 12.5px; color: {MUTED}">{desc}</span></div>'
                f'<div style="display: flex; align-items: stretch; flex-wrap: wrap; row-gap: 10px">' + arr.join(node(*s) for s in steps) + '</div></section>')
    f1 = flow("1 · Amorphous cell to properties", "the core polymer workflow", [
        ("Start.dc.html", "Start", "type SMILES or pick Polymer", ACC), ("PolymerBuilder.dc.html", "Polymer builder", "repeat unit, sequence, tacticity"),
        ("Grow.dc.html", "Grow", "CBMC growth, real force field"), ("Pack.dc.html", "Pack", "regions, overlap removal"), ("Relax.dc.html", "Relax", "minimise, push-off"),
        ("Jobs.dc.html", "Equilibrate", "21-step protocol"), ("Convergence.dc.html", "Converged?", "automatic criteria"), ("Analyze.dc.html", "Analyze", "density, Tg, C∞, δ", OK)])
    f2 = flow("2 · Build, check and export", "hand-built structures to other engines", [
        ("Sketch.dc.html", "Sketch", "2D → 3D"), ("Main.dc.html", "Studio", "edit, measure, clean"), ("Interactions.dc.html", "Checks", "clashes, valences"),
        ("ForceField.dc.html", "Field", "types, missing parameters"), ("Dynamics.dc.html", "Dynamics", "LAMMPS / GROMACS deck, parity"), ("Provenance.dc.html", "Provenance", "what made this file", OK)])
    f3 = flow("3 · When something fails", "recovery without losing work", [
        ("Jobs.dc.html", "Jobs", "running and queued"), ("FailedJob.dc.html", "Failed run", "reason, force spike", ERR), ("React.dc.html", "Adjust settings", "longer relaxation"),
        ("Jobs.dc.html", "Restart", "from checkpoint 40", OK)])
    f4 = flow("4 · Crosslinked network", "epoxy curing to gel point", [
        ("Main.dc.html", "Studio", "resin + hardener"), ("SolvationBuilder.dc.html", "Pack mixture", "stoichiometric box"), ("React.dc.html", "React", "REACTER-style bonding"),
        ("Mechanics.dc.html", "Mechanics", "modulus of the network"), ("Provenance.dc.html", "Provenance", "cited and reproducible", OK)])
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 18px; padding: 32px 40px">'
            + f'<div style="display: flex; align-items: flex-end; gap: 20px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">User flows</h1><span style="font-size: 13px; color: {MUTED}">Press Play, then click a step to open that screen</span></div>'
            + f1 + f2 + f3 + f4 + '</div>')
    return page("CAPS — user flows", body)
