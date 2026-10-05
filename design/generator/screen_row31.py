"""Row 31 — Projects: a CAPS project is a folder the user puts anywhere, with one NAME.capsproj file in it. CAPS keeps a
list of every project it has made or opened (wherever the folder is), opens one with all its structures and the history
of its sessions, and saves it as the user works. Drawn over the real CAPS Studio (Start, and the Studio with a structure)."""
import os
from lib import *

START = "/_blob/7e8eda7c5ae603f2acb7ed68e4083728"    # CAPS Start today, 1440 × 900
STUDIO = "/_blob/a61a8c76400d88abff03fc15d6d22a61"   # CAPS Studio with ps_melt.data, 1440 × 900
H3 = f"margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}"
SHADOW = "box-shadow: 0 18px 48px #000000a0, 0 2px 6px #00000080"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def at(x, y, inner, w=None, h=None, extra=""):
    size = (f"width: {w}px; " if w else "") + (f"height: {h}px; " if h else "")
    return f'<div style="position: absolute; left: {x}px; top: {y}px; {size}{extra}">{inner}</div>'


def board(title, shot, alt, *layers):
    body = (f'<div style="position: relative; width: 1440px; height: 900px; overflow: hidden">'
            f'<img src="{shot}" alt="{alt}" style="position: absolute; left: 0; top: 0; width: 1440px; height: 900px">'
            + "".join(layers) + '</div>')
    return page(title, body)


def cover(x, y, w, h, color=BG1):
    return at(x, y, "", w, h, f"background: {color}")


def panel(inner, w, pad=14, gap=10, extra=""):
    return (f'<div style="width: {w}px; box-sizing: border-box; display: flex; flex-direction: column; gap: {gap}px; padding: {pad}px; background: {BG1}; '
            f'border: 1px solid {LINE}; border-radius: 10px; {SHADOW}; {extra}">{inner}</div>')


def callout(n, text, w=280):
    return (f'<div style="width: {w}px; display: flex; gap: 10px; align-items: flex-start; padding: 10px 12px; background: #1A1408; border: 1px solid {ACC}; '
            f'border-radius: 8px; {SHADOW}"><span style="flex-shrink: 0; width: 20px; height: 20px; border-radius: 50%; background: {ACC}; color: {ACC_INK}; '
            f'font-size: 11.5px; font-weight: 700; display: inline-flex; align-items: center; justify-content: center">{n}</span>'
            f'<span style="font-size: 12px; line-height: 1.45; color: {TEXT}">{text}</span></div>')


def crop(shot, x, y, w, h, alt):
    """A part of a screenshot (a structure's thumbnail) shown at its own size."""
    return (f'<div style="position: relative; width: {w}px; height: {h}px; overflow: hidden; background: {BG0}">'
            f'<img src="{shot}" alt="{alt}" style="position: absolute; left: -{x}px; top: -{y}px; width: 1440px; height: 900px"></div>')


def mono(t, color=DIM, fs=11.5):
    return f'<span style="font-family: {MONO}; font-size: {fs}px; color: {color}; overflow: hidden; white-space: nowrap; text-overflow: ellipsis">{t}</span>'


def project_card(name, path, counts, when, thumb=None, missing=False):
    if missing:
        top = (f'<div style="height: 96px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 6px; '
               f'border-bottom: 1px dashed {LINE}; background: {BG0}">{icon("alert", 20, WARN)}'
               f'<span style="font-size: 12px; color: {WARN}">Folder not found</span></div>')
        foot = row(btn("Locate…", small=True), btn("Forget", small=True), gap=6)
    else:
        top = thumb
        foot = f'<span style="font-size: 11.5px; color: {DIM}">{when}</span>'
    return (f'<a href="#" style="text-decoration: none; color: {TEXT}; width: 176px; display: flex; flex-direction: column; background: {BG1}; '
            f'border: 1px solid {LINE}; border-radius: 10px; overflow: hidden">{top}'
            f'<div style="display: flex; flex-direction: column; gap: 4px; padding: 10px 12px 12px">'
            f'<span style="display: flex; align-items: center; gap: 7px; font-size: 13.5px; font-weight: 600">{icon("folder", 15, ACC)}{name}</span>'
            f'{mono(path)}<span style="font-size: 12px; color: {MUTED}">{counts}</span>{foot}</div></a>')


# ---------------------------------------------------------------- 1. Start: the projects CAPS knows

def projects_start():
    t1 = crop(START, 362, 550, 176, 96, "ENR–MAH crosslinking: its last structure, a crosslinked ENR cell")
    t2 = crop(START, 595, 550, 176, 96, "PE melt Tg study: its last structure, a polyethylene cell")
    cards = row(
        f'<button aria-label="New project" style="width: 120px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 10px; '
        f'background: transparent; border: 1px dashed {ACC}; border-radius: 10px; color: {ACC}; cursor: pointer">{icon("plus", 22, ACC)}'
        f'<span style="font-size: 13px; font-weight: 600">New project</span><span style="font-size: 11.5px; color: {MUTED}; padding: 0 14px; text-align: center">a folder anywhere</span></button>',
        project_card("ENR–MAH crosslinking", "~/Research/ENR-MAH", "4 structures · 6 sessions", "open now", t1),
        project_card("PE melt Tg study", "/Volumes/Lab/CAPS/PE-Tg", "9 structures · 12 sessions", "opened 2 Oct", t2),
        project_card("Fibre–rubber interface", "/Volumes/USB/fibre-rubber", "last seen 21 Sep", "", missing=True),
        gap=12, align="stretch")
    recent = row(*[f'<a href="#" style="display: inline-flex; align-items: center; gap: 6px; height: 26px; padding: 0 10px; border-radius: 13px; background: {BG2}; '
                   f'border: 1px solid {LINE}; color: {TEXT}; text-decoration: none; font-size: 12px">{icon("file", 13, DIM)}{n}</a>'
                   for n in ("ps_melt.data", "pe_xl.data", "enr50_alt.data")], gap=6)
    main = (f'<div style="width: 688px; display: flex; flex-direction: column; gap: 12px">'
            + row(f'<h2 style="margin: 0; font-size: 15px; font-weight: 600">Projects</h2>',
                  f'<span style="font-size: 12px; color: {DIM}">3 known · wherever they are saved</span>',
                  f'<span style="flex-grow: 1"></span>', btn("Open project file…", small=True, ic="folder"), gap=10)
            + cards
            + row(f'<span style="font-size: 12px; color: {MUTED}">Recent files</span>', recent, gap=12)
            + '</div>')
    left = (f'<div style="width: 299px; display: flex; flex-direction: column; gap: 10px; padding: 14px">'
            + row(f'<span style="font-size: 13.5px; font-weight: 600">Project</span>', f'<span style="flex-grow: 1"></span>', icon("plus", 15, MUTED))
            + f'<div style="display: flex; flex-direction: column; gap: 8px; padding: 12px; border: 1px dashed {LINE}; border-radius: 8px">'
              f'<span style="font-size: 12.5px">No project open</span>'
              f'<span style="font-size: 12px; color: {MUTED}; line-height: 1.45">What you build now is kept for this session only. Put it in a project to keep it with its history.</span>'
            + row(btn("New project", small=True, primary=True, ic="plus"), btn("Open…", small=True, ic="folder"), gap=6) + '</div></div>')
    return board("CAPS — Projects on Start", START, "CAPS Start today: builders, recent files, learn and recipes",
                 cover(340, 505, 700, 340, BG0), at(340, 510, main),
                 cover(0, 128, 299, 345), at(0, 128, left),
                 at(1064, 520, callout(1, "Every project CAPS has made or opened, wherever its folder is (a disk, a lab share, a USB stick). "
                                          "The list lives in CAPS's own settings, so CAPS knows them all however they were opened.", 330)),
                 at(1064, 640, callout(2, "A project whose folder moved or whose disk is not attached stays in the list: <b>Locate…</b> finds it again "
                                          "(the .capsproj file in its new place); <b>Forget</b> drops it from the list and leaves the files alone.", 330)),
                 at(1064, 772, callout(3, "Double-clicking a <span style=\"font-family: " + MONO + "\">.capsproj</span> file in Finder, Explorer or a Linux file manager "
                                          "opens CAPS on that project and adds it to the list.", 330)))


# ---------------------------------------------------------------- 2. New project

def new_project():
    tree_rows = [("folder", "ENR–MAH crosslinking/", "the project folder", TEXT),
                 ("file", "  ENR–MAH crosslinking.capsproj", "the project file: open this one", ACC),
                 ("folder", "  structures/", "every structure, saved as you work", TEXT),
                 ("folder", "  runs/", "Minimise, Equilibrate, Dynamics and React runs", TEXT),
                 ("folder", "  exports/", "LAMMPS and GROMACS files you export", TEXT)]
    tree = "".join(f'<div style="display: flex; align-items: center; gap: 8px; height: 24px">{icon(ic, 14, c if c == ACC else DIM)}'
                   f'<span style="font-family: {MONO}; font-size: 12px; color: {c}; white-space: pre; width: 250px">{n}</span>'
                   f'<span style="font-size: 11.5px; color: {DIM}">{d}</span></div>' for ic, n, d, c in tree_rows)
    sheet = panel(
        row(f'<h2 style="margin: 0; font-size: 17px; font-weight: 600">New project</h2>', f'<span style="flex-grow: 1"></span>',
            f'<button aria-label="Close" style="width: 28px; height: 28px; padding: 0; background: transparent; border: 0; cursor: pointer">{icon("close", 16, MUTED)}</button>')
        + field("Name", "ENR–MAH crosslinking", mono=False)
        + f'<div style="display: flex; flex-direction: column; gap: 4px"><span style="font-size: 11.5px; color: {MUTED}">Where</span>'
        + row(f'<div style="flex-grow: 1; height: 30px; display: flex; align-items: center; padding: 0 8px; background: {BG0}; border: 1px solid {LINE}; border-radius: 5px">'
              f'{mono("~/Research", TEXT, 12.5)}</div>', btn("Choose folder…", small=True, ic="folder"), gap=6) + '</div>'
        + f'<h3 style="{H3}; margin-top: 4px">CAPS will make</h3>'
        + f'<div style="padding: 8px 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px">{tree}</div>'
        + check("Bring the 3 structures open now into it (with their force fields and runs)", True)
        + check("Save as I work: after every change, and when CAPS closes", True)
        + row(f'<span style="font-size: 11.5px; color: {DIM}">The folder can be moved or copied whole; CAPS finds it again.</span>',
              f'<span style="flex-grow: 1"></span>', btn("Cancel"), btn("Create project", primary=True), gap=8),
        600, pad=20, gap=14)
    return board("CAPS — New project", START, "CAPS Start, dimmed behind the new project sheet",
                 at(0, 0, "", 1440, 900, "background: #0B0D0Fb8"), at(420, 170, sheet),
                 at(1060, 260, callout(1, "The name becomes the folder and the file: <span style=\"font-family: " + MONO + "\">NAME/NAME.capsproj</span>. "
                                          "Characters a disk cannot hold are replaced; an existing folder is never overwritten.", 330)),
                 at(1060, 400, callout(2, "<span style=\"font-family: " + MONO + "\">.capsproj</span> is CAPS's own extension: a small readable file (JSON) that names "
                                          "the project's structures, its sessions and what was done in each. The structures sit beside it in "
                                          "<span style=\"font-family: " + MONO + "\">structures/</span>, as LAMMPS data with their force field.", 330)),
                 at(1060, 570, callout(3, "Started without a project? Nothing is lost: <b>Bring the structures open now</b> moves this session's work into the new project.", 330)))


# ---------------------------------------------------------------- 3. A project open

def project_open():
    def item(name, sub, active=False):
        bg = f"background: {BG2}; border-left: 2px solid {ACC}; " if active else "border-left: 2px solid transparent; "
        return (f'<a href="#" style="{bg}display: flex; gap: 8px; padding: 6px 10px; text-decoration: none; color: {TEXT}; border-radius: 4px">{icon("file", 14, ACC if active else DIM)}'
                f'<span style="display: flex; flex-direction: column; gap: 1px"><span style="font-size: 12.5px; font-weight: {600 if active else 400}">{name}</span>'
                f'<span style="font-size: 11px; color: {DIM}">{sub}</span></span></a>')

    def session(when, span, what, now=False):
        d = dot(ACC if now else DIM, 8)
        return (f'<div style="display: flex; gap: 10px"><div style="display: flex; flex-direction: column; align-items: center; gap: 0; padding-top: 4px">{d}'
                f'<span style="width: 1px; flex-grow: 1; background: {LINE}"></span></div>'
                f'<div style="display: flex; flex-direction: column; gap: 2px; padding-bottom: 8px"><span style="font-size: 12px; font-weight: 600">{when} '
                f'<span style="font-weight: 400; color: {DIM}">{span}</span></span><span style="font-size: 11.5px; color: {MUTED}; line-height: 1.4">{what}</span></div></div>')

    left = (f'<div style="width: 299px; display: flex; flex-direction: column; gap: 8px; padding: 12px 12px 0">'
            + row(icon("folder", 16, ACC), f'<span style="font-size: 13.5px; font-weight: 600">ENR–MAH crosslinking</span>', f'<span style="flex-grow: 1"></span>', icon("dots", 15, MUTED), gap=7)
            + mono("~/Research/ENR-MAH/ENR–MAH crosslinking.capsproj", DIM, 11)
            + f'<h3 style="{H3}; margin-top: 6px">Structures · 4</h3>'
            + item("ps_melt.data", "1,300 atoms · 0.386 g/cm³ · opened from a file", True)
            + item("ENR-50 alternate 2×10", "298 atoms · OPLS-AA · grown")
            + item("ENR-50 + 2 maleic acid", "packed · reacted, 1 link")
            + item("ENR-50 + 2 maleic acid · relaxed", "minimised · NPT 50 ps")
            + f'<h3 style="{H3}; margin-top: 8px">Sessions · 6</h3>'
            + session("Today", "09:12 – now", "opened ps_melt.data for comparison", True)
            + session("4 Oct", "14:05 – 17:40", "packed 2 maleic acid · crosslinked (ENR + carboxylic acid, 1 link) · minimised · NPT 50 ps")
            + session("3 Oct", "10:20 – 11:02", "grew ENR-50 alternate, 2 chains × DP 10 · assigned OPLS-AA")
            + f'<a href="#" style="font-size: 12px; padding-left: 18px">3 earlier sessions</a></div>')

    menu_rows = []
    for nm, path, cur in (("ENR–MAH crosslinking", "~/Research/ENR-MAH", True), ("PE melt Tg study", "/Volumes/Lab/CAPS/PE-Tg", False)):
        menu_rows.append(f'<a href="#" style="display: flex; gap: 8px; align-items: flex-start; padding: 7px 10px; border-radius: 5px; text-decoration: none; color: {TEXT}; '
                         f'{"background: " + BG3 + "; " if cur else ""}">{icon("check" if cur else "folder", 14, ACC if cur else DIM)}'
                         f'<span style="display: flex; flex-direction: column; gap: 1px"><span style="font-size: 12.5px">{nm}</span>{mono(path, DIM, 11)}</span></a>')
    menu_rows.append(f'<div style="display: flex; gap: 8px; align-items: center; padding: 7px 10px; color: {DIM}">{icon("alert", 14, WARN)}'
                     f'<span style="font-size: 12.5px">Fibre–rubber interface</span><span style="flex-grow: 1"></span><a href="#" style="font-size: 12px">Locate…</a></div>')
    sep_ = f'<div style="height: 1px; background: {LINE}; margin: 4px 0"></div>'
    def mi(ic, t, key=""):
        k = f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{key}</span>' if key else ""
        return (f'<button style="display: flex; align-items: center; gap: 8px; width: 100%; height: 30px; padding: 0 10px; background: transparent; border: 0; '
                f'border-radius: 5px; cursor: pointer; font-size: 12.5px; color: {TEXT}; text-align: left">{icon(ic, 14, MUTED)}<span style="flex-grow: 1">{t}</span>{k}</button>')
    menu = (f'<div role="menu" style="width: 300px; box-sizing: border-box; padding: 5px; background: {BG2}; border: 1px solid {LINE}; border-radius: 9px; {SHADOW}">'
            + "".join(menu_rows) + sep_ + mi("save", "Save", "⌘S") + mi("copy", "Save a copy as…", "⇧⌘S") + mi("plus", "New project…")
            + mi("folder", "Open project file…", "⌘O") + mi("file", "Show in Finder") + sep_ + mi("close", "Close project") + '</div>')
    chip = (f'<button style="display: inline-flex; align-items: center; gap: 8px; height: 32px; padding: 0 10px 0 9px; background: {BG2}; border: 1px solid {ACC}; '
            f'border-radius: 7px; cursor: pointer; color: {TEXT}">{icon("folder", 15, ACC)}<span style="font-size: 13px; font-weight: 600">ENR–MAH crosslinking</span>'
            f'<span style="font-size: 11.5px; color: {DIM}">saved 10:42</span>{icon("chev", 13, MUTED)}</button>')
    tab = (f'<div style="display: inline-flex; align-items: center; gap: 8px; height: 34px; padding: 0 12px; background: {BG1}; border: 1px solid {LINE}; border-bottom: 0; '
           f'border-radius: 7px 7px 0 0">{icon("file", 14, ACC)}<span style="font-size: 13px">ps_melt.data</span>{icon("close", 12, DIM)}</div>')
    return board("CAPS — A project open", STUDIO, "CAPS Studio with ps_melt.data open",
                 cover(0, 126, 299, 748), at(0, 126, left),
                 cover(92, 36, 510, 40, BG0), at(96, 40, row(chip, tab, f'<span style="color: {MUTED}">{icon("plus", 15, MUTED)}</span>', gap=10)),
                 at(96, 78, menu),
                 at(1108, 540, callout(1, "The project's name and when it was last saved, always in view. CAPS saves after each change (a structure built, "
                                          "a run finished, a rename) and on closing; ⌘S saves at once.", 312)),
                 at(1108, 668, callout(2, "Each time the project is opened is a session: when, how long, and what was built or run in it, written from the "
                                          "structures' own history. Opening the project brings every structure back where it was left.", 312)),
                 at(1108, 800, callout(3, "Runs still going when CAPS closes stay with their structure and are picked up when the project opens again.", 312)))


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for fn, f in (("ProjectsStart.dc.html", projects_start), ("NewProject.dc.html", new_project), ("ProjectOpen.dc.html", project_open)):
        open(os.path.join(OUT, fn), "w").write(f())
        print("wrote", fn)
