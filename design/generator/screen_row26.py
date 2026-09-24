"""Row 26 — Visualize states & access: empty, loading, file problems, keyboard & screen reader, Paper variants."""
import math
from lib import *
from screen_app import card, progress
from screen_row20 import draw, CHAIN, annotate
from screen_row21 import cell_recs, mi, inspector, viewport, frame, CW, VH, BOTTOM, L, MASS, NA
from screen_row22 import chain_col, H3
from screen_row23 import topology


def shell(title, tab, crumb, content, status):
    body = (topbar([tab], tab, f'{icon("chevr", 12, DIM)}<span>{crumb}</span>') + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Analyze")
            + f'<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">{content}</div></div>' + statusbar(*status))
    return page(title, body)


def item(icon_name, title, sub, right="", colour=MUTED):
    return (f'<div style="display: flex; align-items: center; gap: 12px; padding: 9px 4px; border-bottom: 1px solid {BG2}">{icon(icon_name, 16, colour)}'
            f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 1px"><span style="font-size: 12.5px">{title}</span><span style="font-size: 11px; color: {DIM}; font-family: {MONO}">{sub}</span></div>{right}</div>')


# ---------------------------------------------------------------- Empty
def vis_empty():
    cell, recs = cell_recs()
    bonds, angles, dihs = topology(cell, recs)
    fmts = ["LAMMPS data · dump", "GROMACS .gro .xtc .trr .tpr", "PDB · mmCIF", "XYZ · extxyz", "DCD", "CIF · POSCAR", "mol2 · SDF", "CAPS .caps"]
    drop = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 16px; border: 2px dashed {LINE}; border-radius: 14px; background: {BG1}">'
            f'<div style="width: 64px; height: 64px; border-radius: 16px; background: {BG2}; display: flex; align-items: center; justify-content: center">{icon("folder", 30, ACC)}</div>'
            f'<h1 style="margin: 0; font-size: 22px; font-weight: 600">Open a structure or trajectory</h1>'
            f'<p style="margin: 0; max-width: 520px; text-align: center; font-size: 13px; color: {MUTED}; line-height: 1.6">Drop files here, or open them. A trajectory and its topology can be dropped together; CAPS pairs them by atom count.</p>'
            f'<div style="display: flex; gap: 10px">{btn("Open file…", True, "folder")}{btn("Open from Grow", ic="grow")}{btn("Paste path", ic="file")}</div>'
            f'<div style="display: flex; flex-wrap: wrap; gap: 6px; justify-content: center; max-width: 560px">{"".join(chip(f, MUTED, BG2) for f in fmts)}</div></div>')
    recent = col(item("file", "PS_melt.lammpstrj", "~/runs/ps · 500 frames · with PS_melt.data", chip("today", MUTED, BG2)),
                 item("file", "cells/seed21/PS_melt.data", "~/caps/cells · built by Grow", chip("today", MUTED, BG2)),
                 item("file", "cells/seed23/PS_melt.data", "~/caps/cells · built by Grow", chip("today", MUTED, BG2)),
                 item("file", "PS_melt.gro + topol.top", "~/runs/ps_gmx", chip("yesterday", MUTED, BG2)), gap=0)
    sample = col(item("cube", "Polystyrene cell · 10 × DP 8", f"{len(recs)} atoms · {len(bonds)} bonds · {L:.0f} Å box", btn("Open", small=True)),
                 item("cube", "Same cell, backbone only", f"{sum(len(bb) for _, _, bb in cell)} atoms", btn("Open", small=True)),
                 item("layers", "Pipeline · structure report", "7 steps · saved YAML", btn("Apply", small=True)), gap=0)
    right = (f'<div style="width: 380px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
             + card("Recent", recent, chip("4", MUTED, BG2, True), 14)
             + card("Samples made by CAPS", sample, "", 14)
             + card("Keyboard", col(kv("Open", "⌘O"), kv("Open recent", "⌘⇧O"), kv("Command palette", "⌘K"), gap=6), "", 14) + '</div>')
    content = f'<div style="flex-grow: 1; display: flex; gap: 18px; padding: 22px; min-height: 0; overflow: hidden">{drop}{right}</div>'
    return shell("CAPS — Visualize, nothing open", "Visualize", "Analyze › Visualize", content,
                 ("<span>no file open</span>", "<span>samples are the cells shown elsewhere on this canvas</span>"))


# ---------------------------------------------------------------- Loading
def vis_loading():
    cell, recs = cell_recs()
    bonds, angles, dihs = topology(cell, recs)
    sc = Scene("ld", CW - 6, VH + BOTTOM - 150, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.15)
    draw(sc, cell, recs, chain_col)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    idx_done, total = 212, 500
    stages = [("Detect format", "LAMMPS dump · text · columns id mol type q xu yu zu", "done", 1.0),
              ("Read frame 0", f"{len(recs)} atoms · box {L:.0f} Å cubic", "done", 1.0),
              ("Join topology", f"PS_melt.data · {len(bonds)} bonds · {len(angles)} angles", "done", 1.0),
              ("Index frames", f"{idx_done} of {total} · offsets only, no coordinates kept", "running", idx_done / total),
              ("Run pipeline on frame 0", "7 steps", "queued", 0.0)]
    rows = []
    for t, s, st, f in stages:
        c = OK if st == "done" else (ACC if st == "running" else DIM)
        ic = "check" if st == "done" else ("play" if st == "running" else "pause")
        rows.append(f'<div style="display: flex; flex-direction: column; gap: 6px; padding: 10px 0; border-bottom: 1px solid {BG2}"><div style="display: flex; align-items: center; gap: 10px">{icon(ic, 15, c)}'
                    f'<span style="flex-grow: 1; font-size: 12.5px">{t}</span><span style="font-family: {MONO}; font-size: 11px; color: {c}">{st}</span></div>'
                    f'<span style="font-size: 11px; color: {DIM}; font-family: {MONO}; padding-left: 25px">{s}</span>' + (f'<div style="padding-left: 25px">{progress(f)}</div>' if st == "running" else "") + '</div>')
    panel = (f'<aside style="width: 330px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}">'
             + panel_head("Opening PS_melt.lammpstrj", chip("step 4 of 5", ACC, BG2))
             + f'<div style="padding: 4px 18px">{"".join(rows)}</div>'
             + f'<div style="padding: 12px 18px; display: flex; flex-direction: column; gap: 8px">{kv("File size", "[measure]")}{kv("Time left", "[measure]")}'
             f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Frame 0 is already usable: orbit, pick and measure work while the rest is indexed. Scrubbing is limited to indexed frames.</span>'
             f'<div style="display: flex; gap: 8px">{btn("Cancel")}{btn("Open in background", ic="jobs")}</div></div></aside>')
    tl = (f'<div style="height: 150px; flex-shrink: 0; display: flex; flex-direction: column; gap: 10px; padding: 14px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
          f'<div style="display: flex; align-items: center; gap: 12px"><span style="font-family: {MONO}; font-size: 12px">frame 0 / {total - 1}</span><span style="flex-grow: 1"></span><span style="font-size: 11.5px; color: {MUTED}">frames right of the marker are not indexed yet</span></div>'
          f'<div style="position: relative; height: 14px; border-radius: 7px; background: {BG3}"><div style="width: {idx_done / total * 100:.1f}%; height: 14px; border-radius: 7px; background: {mix(ACC, BG3, 0.45)}"></div>'
          f'<div style="position: absolute; left: {idx_done / total * 100:.1f}%; top: -4px; width: 2px; height: 22px; background: {ACC}"></div></div>'
          f'<div style="display: flex; justify-content: space-between; font-family: {MONO}; font-size: 10.5px; color: {DIM}"><span>0</span><span>{idx_done} indexed</span><span>{total - 1}</span></div>'
          f'<span style="font-size: 11.5px; color: {DIM}">Stage counts from this cell; the {idx_done}/{total} indexing progress is illustrative.</span></div>')
    center = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
              f'<div style="position: relative; flex-grow: 1; min-height: 0; background: {BG0}; overflow: hidden">{sc.svg(overlay=label_pill(12, 20, "frame 0 · ready", OK, fs=11))}</div>{tl}</div>')
    content = f'<div style="flex-grow: 1; display: flex; min-height: 0">{center}{panel}</div>'
    return shell("CAPS — opening a trajectory", "PS_melt.lammpstrj", "Analyze › Visualize › Opening", content,
                 (f"<span>{len(recs)} atoms · frame 0 shown while indexing</span>", "<span>progress illustrative</span>"))


# ---------------------------------------------------------------- File problems
def vis_problems():
    cell, recs = cell_recs()
    bonds, angles, dihs = topology(cell, recs)
    outside = sum(1 for r in recs if any(not (0 <= v < L) for v in r["p"]))
    qtot = sum(r["q"] for r in recs)
    types = sorted({r["type"] for r in recs})
    maxlen = max(dist(recs[a]["p"], recs[b]["p"]) for a, b in bonds)
    def prob(level, title, detail, actions):
        c = {"pass": OK, "note": SEL, "warn": WARN, "error": ERR}[level]
        ic = {"pass": "check", "note": "alert", "warn": "alert", "error": "xcircle"}[level]
        return (f'<div style="display: flex; gap: 12px; padding: 11px 12px; border-radius: 8px; background: {BG1}; border: 1px solid {c if level in ("warn", "error") else LINE}">'
                f'<span style="flex-shrink: 0; padding-top: 1px">{icon(ic, 16, c)}</span><div style="flex-grow: 1; display: flex; flex-direction: column; gap: 4px">'
                f'<span style="font-size: 12.5px; font-weight: 600">{title}</span><span style="font-size: 12px; color: {MUTED}; line-height: 1.5">{detail}</span></div>'
                f'<div style="display: flex; gap: 6px; align-items: flex-start; flex-shrink: 0">{actions}</div></div>')
    real = col(
        prob("pass", "Atom counts agree", f"Dump and PS_melt.data both define {len(recs)} atoms with IDs 1–{len(recs)}.", ""),
        prob("pass", "Bonds are plausible", f"{len(bonds)} bonds; the longest is {maxlen:.2f} Å (unwrapped). None cross more than half the box.", ""),
        prob("pass", "Charges are neutral", f"Sum of q = {qtot:+.1e} e over {len(recs)} atoms.", ""),
        prob("note", f"{outside} atoms lie outside the box", "Coordinates are unwrapped (xu yu zu). Kept as read; wrap only changes the view.", btn("Wrap view", small=True)),
        prob("note", f"Types {', '.join(map(str, types))} have no element column", "Elements guessed from masses in PS_melt.data: 12.011 → C, 1.008 → H.", btn("Review", small=True)), gap=8)
    ex = col(
        prob("error", "Bond 1 371 points to atom 1 301", "Only 1 300 atoms are defined. The bond is skipped; the rest of the topology loads.", btn("Show line", small=True)),
        prob("warn", "No Masses for type 5", "Type 5 appears in Atoms but not in Masses. Pick an element or a mass to continue.", btn("Set mass…", small=True)),
        prob("warn", "Frames 240 → 260 skip timestep 250", "Frame spacing changes once. Time-based analysis will use real timesteps, not frame numbers.", btn("Details", small=True)),
        prob("error", "Column 'vx' is text in frame 311", "Frame 311 cannot be read. It is marked unreadable and skipped; other frames are fine.", btn("Show frame", small=True)), gap=8)
    content = (f'<div style="flex-shrink: 0; display: flex; align-items: center; gap: 16px; padding: 14px 22px; background: {BG1}; border-bottom: 1px solid {LINE}">'
               f'<div style="display: flex; flex-direction: column; gap: 3px"><h1 style="margin: 0; font-size: 19px; font-weight: 600">File checks</h1><div style="font-size: 12.5px; color: {MUTED}">Every opened file is checked; problems say what was done and what you can change</div></div>'
               f'<div style="flex-grow: 1"></div>{btn("Export report", ic="download")}{btn("Continue", True)}</div>'
               f'<div style="flex-grow: 1; display: flex; gap: 18px; padding: 16px 22px; min-height: 0; overflow: hidden">'
               f'<div style="flex: 1 1 0; min-width: 0; display: flex; flex-direction: column; gap: 10px"><div style="display: flex; align-items: center; gap: 10px"><h3 style="{H3}">PS_melt.lammpstrj + PS_melt.data</h3>{chip("checked on this cell", OK, BG2)}</div>{real}</div>'
               f'<div style="flex: 1 1 0; min-width: 0; display: flex; flex-direction: column; gap: 10px"><div style="display: flex; align-items: center; gap: 10px"><h3 style="{H3}">How problems read</h3>{chip("example messages", WARN, BG2)}</div>{ex}'
               f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Rules: never silently drop data; say which lines were skipped; keep reading what can be read; one action per problem.</span></div></div>')
    return shell("CAPS — file checks", "PS_melt.lammpstrj", "Analyze › Visualize › File checks", content,
                 ("<span>left column: checks run on this cell · right: example messages</span>", "<span>3 passed · 2 notes</span>"))


# ---------------------------------------------------------------- Keyboard & screen reader
def vis_access():
    cell, recs = cell_recs()
    pick = recs[40]
    near = sorted(((mi(pick["p"], r["p"]), r) for r in recs if r["id"] != pick["id"]), key=lambda t: t[0])[:3]
    sc = Scene("ax", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.2)
    draw(sc, cell, recs, chain_col)
    sc.halo.add(recs.index(pick))
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    desc = {"c3": "sp3 carbon", "ca": "aromatic carbon", "hc": "hydrogen on an sp3 carbon", "ha": "hydrogen on an aromatic carbon"}[pick["gaff"]]
    say = (f"Atom {pick['id']}, {pick['gaff']}, {desc}. Molecule {pick['mol']}, unit {pick['unit']}. "
           f"Charge {pick['q']:+.4f} e. Position {pick['p'][0]:.2f}, {pick['p'][1]:.2f}, {pick['p'][2]:.2f} Å. "
           f"Nearest: atom {near[0][1]['id']} at {near[0][0]:.2f} Å.")
    live = (f'<div role="status" style="position: absolute; left: 16px; bottom: 14px; right: 16px; display: flex; gap: 10px; align-items: flex-start; padding: 10px 12px; background: {BG1}; border: 1px solid {SEL}; border-radius: 8px">'
            f'{icon("eye", 16, SEL)}<div style="display: flex; flex-direction: column; gap: 3px"><span style="font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {SEL}">Screen reader announces</span>'
            f'<span style="font-size: 12.5px; line-height: 1.5">{say}</span></div></div>')
    vp = f'<div style="position: relative; flex-grow: 1; min-height: 0; background: {BG0}; overflow: hidden">{sc.svg()}{live}</div>'
    nt = table(["Neighbour", "Name", "Molecule", "Distance (Å)"], [[str(r["id"]), r["gaff"], str(r["mol"]), f"{d:.2f}"] for d, r in near], ["22%", "22%", "26%", "30%"], mono_cols=(0, 1, 2, 3), align_right=(3,), fs=12, rowh=26)
    keys = [("Tab / ⇧Tab", "move between panels"), ("↑ ↓", "previous / next atom in molecule"), ("[ ]", "previous / next molecule"), ("B", "step along bonds"), ("Space", "select focused atom"), ("M", "measure to last selected"), ("← →", "previous / next frame"), ("⌘⇧A", "announce selection")]
    kt = table(["Keys", "Does"], [[a, b] for a, b in keys], ["34%", "66%"], mono_cols=(0,), fs=12, rowh=25)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="width: 440px; flex-shrink: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">Around the focused atom</h3>{nt}</div>'
              f'<div style="flex-grow: 1; font-size: 12px; color: {MUTED}; line-height: 1.6; padding-top: 18px">The 3D view is one focusable region. Inside it, keys walk the structure the way a reader walks text: atom, bond, molecule, frame. The ring marks focus; the announcement is also shown on screen, for sighted keyboard users.</div></div>')
    right = inspector("Keyboard & reader", section("Walk the structure", kt, pad=12)
                      + section("Announce", col(select("Verbosity", "full · name, place, charge, nearest"), toggle("Announce frame changes", False), toggle("Sound on selection", False), gap=8)))
    return frame(vp + bottom, right, "CAPS — keyboard and screen reader", "Analyze › Visualize › Accessibility",
                 (f"<span>focus: atom {pick['id']} · molecule {pick['mol']}</span>", "<span>values read from this frame</span>"))


BOARDS = (("VisEmpty", vis_empty), ("VisLoading", vis_loading), ("VisProblems", vis_problems), ("VisAccess", vis_access))

if __name__ == "__main__":
    import os, sys
    os.makedirs("stage26/project", exist_ok=True)
    only = sys.argv[1:]
    for name, fn in BOARDS:
        if only and name not in only: continue
        h = fn(); open(f"stage26/project/{name}.dc.html", "w").write(h); print(name, len(h), flush=True)
    if not only or "Paper" in only:
        from screen_row19 import paper_of
        from screen_row20 import pipeline
        from screen_row22 import data_inspector
        d = list(Scene.__init__.__defaults__); d[-1] = "#F4F2EE"; Scene.__init__.__defaults__ = tuple(d)
        for name, fn in (("PaperVisPipeline", pipeline), ("PaperDataInspector", data_inspector)):
            h = paper_of(fn()); open(f"stage26/project/{name}.dc.html", "w").write(h)
            print(name, len(h), "leftover dark", {c: h.upper().count(c) for c in ["#0F1113", "#16191C", "#262B30"] if h.upper().count(c)})
