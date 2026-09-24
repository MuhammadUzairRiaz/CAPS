import math, random
from lib import *
from mols import *
from screen_app import card, pagehead, shell
from screen_studio2 import spacer, bar
from screen_row7 import with_overlay
from screen_row10 import logo_mark
import screen_main


# ------------------------------------------------------------------ First-run tour
def tour():
    base = screen_main.build()
    # spotlight the toolbar row (x 72..1440, y 44..92) with a box-shadow cut-out
    spot = (f'<div style="position: absolute; left: 76px; top: 47px; width: 1360px; height: 42px; border-radius: 8px; border: 2px solid {ACC}; box-shadow: 0 0 0 9999px #07080AC0"></div>')
    dots = "".join(f'<span style="width: {18 if i == 1 else 7}px; height: 7px; border-radius: 4px; background: {ACC if i == 1 else BG3}"></span>' for i in range(5))
    callout = (f'<div role="dialog" aria-label="Tour step 2 of 5" style="position: absolute; left: 150px; top: 108px; width: 420px; padding: 18px; background: {BG1}; border: 1px solid {ACC}; border-radius: 12px; box-shadow: 0 20px 50px #000000A0; display: flex; flex-direction: column; gap: 12px">'
               f'<span style="position: absolute; left: 40px; top: -8px; width: 14px; height: 14px; background: {BG1}; border-left: 1px solid {ACC}; border-top: 1px solid {ACC}; transform: rotate(45deg)"></span>'
               f'<div style="display: flex; align-items: center; gap: 8px">{chip("2 of 5", ACC, "#3A2C14", True)}<span style="font-size: 15px; font-weight: 600">Tools act on what you select</span></div>'
               f'<p style="margin: 0; font-size: 13px; line-height: 1.55; color: {MUTED}">Pick <b style="color: {TEXT}; font-weight: 600">Select</b> or <b style="color: {TEXT}; font-weight: 600">Build</b>, then click in the 3D view. Every action here is a command you can undo, replay, or copy as Python.</p>'
               f'<div style="display: flex; gap: 8px; align-items: center; font-size: 12px; color: {DIM}"><kbd style="font-family: {MONO}; border: 1px solid {LINE}; border-radius: 4px; padding: 1px 6px; color: {MUTED}">S</kbd> select <kbd style="font-family: {MONO}; border: 1px solid {LINE}; border-radius: 4px; padding: 1px 6px; color: {MUTED}">B</kbd> build <kbd style="font-family: {MONO}; border: 1px solid {LINE}; border-radius: 4px; padding: 1px 6px; color: {MUTED}">⌘K</kbd> commands</div>'
               f'<div style="display: flex; align-items: center; gap: 8px"><div style="display: flex; gap: 4px">{dots}</div>{spacer()}{btn("Skip tour", small=True, href="Main.dc.html")}{btn("Back", small=True)}{btn("Next", True, small=True, href="Main.dc.html")}</div></div>')
    steps = [("1", "Modules", "Studio, then Grow, Pack, Relax … in workflow order", False), ("2", "Tools", "select, build, measure, clean", True),
             ("3", "3D view", "orbit, pan, zoom · live monitors", False), ("4", "Inspector", "atom, molecule, force-field type", False), ("5", "Console", "history and Python replay", False)]
    st = "".join(f'<div style="display: flex; gap: 10px; align-items: flex-start; opacity: {1 if on else 0.75}"><span style="width: 22px; height: 22px; border-radius: 50%; flex-shrink: 0; display: flex; align-items: center; justify-content: center; font-family: {MONO}; font-size: 11px; background: {ACC if on else BG3}; color: {ACC_INK if on else MUTED}">{n}</span>'
                 f'<div style="display: flex; flex-direction: column; gap: 2px"><span style="font-size: 12.5px; font-weight: 600; color: {TEXT}">{t}</span><span style="font-size: 11.5px; color: {MUTED}">{d}</span></div></div>' for n, t, d, on in steps)
    outline = (f'<div style="position: absolute; right: 24px; bottom: 60px; width: 290px; padding: 16px; background: {BG1}; border: 1px solid {LINE}; border-radius: 12px; display: flex; flex-direction: column; gap: 12px">'
               f'<span style="font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Studio tour · 2 min</span>{st}</div>')
    return with_overlay(base, spot + callout + outline).replace("<title>CAPS Studio — builder</title>", "<title>CAPS Studio — first-run tour</title>")


# ------------------------------------------------------------------ Compact 1280 x 800
def compact():
    W, H = 1280, 800
    atoms, backbone, stereo = polystyrene(7, "atactic", seed=5)
    VW = W - 52 - 40
    VH = H - 40 - 44 - 30 - 24
    sc = Scene("cp", VW, VH, yaw=-0.28, pitch=0.42, roll=0.05, persp=0.25, fog=0.6)
    sc.add_atoms(atoms)
    svg = sc.svg(overlay=gizmo(40, VH - 40, -0.28, 0.42))
    top = (f'<header style="height: 40px; flex-shrink: 0; display: flex; align-items: center; gap: 10px; padding: 0 10px; background: {BG1}; border-bottom: 1px solid {LINE}">'
           f'{logo_mark(24, TEXT, ACC)}<button style="height: 28px; padding: 0 10px; display: flex; align-items: center; gap: 6px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px; font-size: 12px; color: {TEXT}; cursor: pointer">{icon("cube", 13, ACC)}PS_atactic_DP40.caps{icon("chev", 12, DIM)}</button>'
           f'<span style="font-size: 11.5px; color: {DIM}">+2 open</span>{spacer()}'
           f'<button aria-label="Search commands" style="height: 28px; padding: 0 10px; display: flex; align-items: center; gap: 6px; background: {BG2}; border: 1px solid {LINE}; border-radius: 6px; color: {DIM}; font-size: 12px; cursor: pointer">{icon("search", 13, DIM)}⌘K</button>{dot(OK)}</header>')
    rail_items = "".join(f'<a href="{MODULE_LINK[n]}" aria-label="{n}" title="{n}" style="width: 40px; height: 40px; display: flex; align-items: center; justify-content: center; border-radius: 8px; background: {BG3 if n == "Studio" else "transparent"}">{icon(ic, 18, ACC if n == "Studio" else MUTED)}</a>' for n, ic in MODULES)
    rail_ = f'<nav aria-label="Modules" style="width: 52px; flex-shrink: 0; display: flex; flex-direction: column; align-items: center; gap: 2px; padding-top: 6px; background: {BG1}; border-right: 1px solid {LINE}">{rail_items}</nav>'
    tb = (f'<div role="toolbar" style="height: 44px; flex-shrink: 0; display: flex; align-items: center; gap: 2px; padding: 0 8px; background: {BG1}; border-bottom: 1px solid {LINE}">'
          + tbtn("cursor", "Select", True) + tbtn("move", "Translate") + tbtn("rotate", "Rotate") + sep() + tbtn("atom", "Place atom") + tbtn("bond", "Draw bond") + tbtn("hex", "Fragment")
          + sep() + tbtn("ruler", "Measure") + tbtn("wand", "Auto-clean", True) + tbtn("undo", "Undo") + tbtn("redo", "Redo") + spacer()
          + f'<button aria-label="More tools" style="height: 32px; padding: 0 10px; display: flex; align-items: center; gap: 6px; background: transparent; border: 1px solid {LINE}; border-radius: 6px; font-size: 12px; color: {MUTED}; cursor: pointer">{icon("dots", 15, MUTED)}More</button></div>')
    drawer = (f'<aside aria-label="Inspector" style="position: absolute; right: 0; top: 0; bottom: 0; width: 300px; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; box-shadow: -16px 0 32px #00000070">'
              + panel_head("Inspector", tbtn("close", "Close drawer"))
              + section("Selected atom · C12", col(kv("Element", "C · sp³"), kv("Type", "c3 · GAFF2"), kv("Charge", "−0.0960 e"), kv("Stereo", "R"), gap=7))
              + section("Validation", col(f'<div style="display: flex; gap: 8px; font-size: 12px">{icon("alert", 15, WARN)}<span>2 undefined stereocentres</span></div>', gap=6), chip("2", WARN, "#3A2C14"))
              + f'<div style="padding: 12px 14px; font-size: 11.5px; color: {DIM}; line-height: 1.5">Below 1440 px wide, side panels become drawers over the view, the project tree moves into the document menu, and the console collapses to a bar.</div></aside>')
    side_tabs = "".join(f'<button aria-label="{n}" style="width: 40px; height: 88px; display: flex; align-items: center; justify-content: center; background: {BG3 if n == "Inspector" else "transparent"}; border: 0; border-left: 2px solid {ACC if n == "Inspector" else "transparent"}; cursor: pointer">'
                        f'<span style="transform: rotate(90deg); white-space: nowrap; font-size: 11.5px; color: {TEXT if n == "Inspector" else MUTED}">{n}</span></button>' for n in ("Inspector", "Fragments", "Project", "Monitors"))
    tabs_ = f'<div style="width: 40px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}">{side_tabs}</div>'
    vp = f'<div style="position: relative; flex-grow: 1; min-width: 0; background: {BG0}">{svg}{drawer}</div>'
    dockbar = (f'<div style="height: 30px; flex-shrink: 0; display: flex; align-items: center; gap: 14px; padding: 0 12px; background: {BG1}; border-top: 1px solid {LINE}; font-size: 12px; color: {MUTED}">'
               f'{icon("history", 14, DIM)}<span>History 15</span>{icon("terminal", 14, DIM)}<span>Console</span>{spacer()}<span style="font-family: {MONO}; font-size: 11px">clean.run ff=gaff2 · converged</span>{icon("chev", 13, DIM)}</div>')
    body = (top + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail_
            + f'<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">{tb}<div style="flex-grow: 1; display: flex; min-height: 0">{vp}{tabs_}</div>{dockbar}</div></div>'
            + f'<footer style="height: 24px; flex-shrink: 0; display: flex; align-items: center; justify-content: space-between; padding: 0 12px; background: {BG1}; border-top: 1px solid {LINE}; font-family: {MONO}; font-size: 10.5px; color: {MUTED}"><span>{len(atoms)} atoms · 7 selected</span><span>compact layout · 1280 × 800</span></footer>')
    return page("CAPS Studio — compact layout", body, W, H)


# ------------------------------------------------------------------ Theory manual
def theory():
    tree = [("Integrators", ["Velocity Verlet", "r-RESPA"]), ("Thermostats", ["Nosé–Hoover chains", "CSVR (Bussi)", "Langevin · BAOAB"]),
            ("Barostats", ["Berendsen", "MTK", "Stochastic cell rescaling"]), ("Electrostatics", ["Ewald", "Smooth PME", "Damped shifted force"]),
            ("Growth", ["Theodorou–Suter", "Configurational bias"]), ("Packing", ["Overlap penalty"])]
    tl = ""
    for g, items in tree:
        tl += f'<div style="padding: 10px 12px 4px; font-size: 11px; font-weight: 600; letter-spacing: 0.06em; text-transform: uppercase; color: {DIM}">{g}</div>'
        tl += "".join(f'<a href="#" aria-current="{"page" if it == "CSVR (Bussi)" else "false"}" style="display: block; height: 30px; line-height: 30px; padding: 0 12px 0 20px; border-radius: 6px; text-decoration: none; font-size: 12.5px; background: {BG3 if it == "CSVR (Bussi)" else "transparent"}; color: {TEXT if it == "CSVR (Bussi)" else MUTED}">{it}</a>' for it in items)
    nav = (f'<nav aria-label="Theory manual" style="width: 260px; flex-shrink: 0; padding: 12px 8px; background: {BG1}; border-right: 1px solid {LINE}">'
           f'<div style="display: flex; align-items: center; gap: 8px; height: 32px; padding: 0 8px; margin: 0 4px 6px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px; color: {DIM}; font-size: 12px">{icon("search", 14, DIM)}<span>Search the manual</span></div>{tl}</nav>')
    eq = (f'<div style="display: flex; align-items: center; justify-content: space-between; padding: 18px 22px; background: {BG0}; border: 1px solid {LINE}; border-radius: 10px">'
          f'<div style="font-family: Georgia, \'Times New Roman\', serif; font-size: 22px; color: {TEXT}; letter-spacing: 0.01em"><i>dK</i> = (<i>K̄</i> − <i>K</i>) <i>dt</i>/<i>τ</i> + 2 √( <i>K K̄</i> / <i>N</i><sub>f</sub> ) <i>dW</i> / √<i>τ</i></div>'
          f'<span style="font-family: {MONO}; font-size: 12px; color: {DIM}">Bussi et al. 2007, eq. 7</span></div>')
    syms = table(["Symbol", "Meaning", "CAPS setting"], [["K", "instantaneous kinetic energy", "—"], ["K̄", "target kinetic energy, N<sub>f</sub> k<sub>B</sub>T / 2", "T"],
                                                          ["τ", "thermostat relaxation time", "τ<sub>T</sub> · default 0.1 ps"], ["N<sub>f</sub>", "degrees of freedom after constraints", "computed"],
                                                          ["dW", "Wiener noise", "Philox stream, keyed by step"]], ["16%", "50%", "34%"], mono_cols=(0,), fs=12.5)
    art = (f'<article style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 16px; padding: 30px 40px; max-width: 820px">'
           f'<div style="display: flex; flex-direction: column; gap: 6px"><span style="font-size: 12px; color: {DIM}">Thermostats</span><h1 style="margin: 0; font-size: 26px; font-weight: 600; letter-spacing: -0.015em">Stochastic velocity rescaling (CSVR)</h1></div>'
           f'<p style="margin: 0; font-size: 14px; line-height: 1.65; color: {MUTED}">CSVR rescales all velocities by one factor each step. The factor is drawn so that the kinetic energy follows the stochastic equation below, which samples the canonical ensemble exactly while disturbing the dynamics as little as a Berendsen thermostat does.</p>'
           + eq
           + syms
           + f'<h2 style="margin: 8px 0 0; font-size: 15px; font-weight: 600">When to use it</h2>'
           f'<p style="margin: 0; font-size: 14px; line-height: 1.65; color: {MUTED}">Default for NVT and NPT production in CAPS. Prefer Nosé–Hoover chains when you need a deterministic, time-reversible trajectory without noise; prefer Langevin (BAOAB) for implicit-solvent or coarse-grained systems.</p>'
           f'<h2 style="margin: 8px 0 0; font-size: 15px; font-weight: 600">Reference</h2>'
           + cite('G. Bussi, D. Donadio, M. Parrinello, "Canonical sampling through velocity rescaling", <i>J. Chem. Phys.</i> 126, 014101 (2007). doi:10.1063/1.2408420')
           + '</article>')
    side = (f'<aside style="width: 320px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 30px 28px 30px 0">'
            + card("Implementation", col(kv("Source", "dynamics/thermo/csvr.cpp"), kv("Tested by", "equipartition, NVT ⟨K⟩"), kv("Deviation from paper", "none", False), gap=8))
            + card("Used in this project", col(*[f'<a href="#" style="display: flex; align-items: center; gap: 8px; font-size: 12.5px; color: {TEXT}; text-decoration: none">{icon("jobs", 14, DIM)}{j}</a>' for j in ("dyn-12 · NPT 10 ns", "equil-7 · steps 2, 5, 8 …")], gap=8))
            + card("Cite", col(f'<span style="font-size: 12px; color: {MUTED}">CAPS adds this reference to every output that used CSVR.</span>', btn("Copy BibTeX", ic="copy", small=True), gap=8))
            + '</aside>')
    content = f'<div style="flex-grow: 1; display: flex; min-height: 0">{nav}{art}{side}</div>'
    return shell("Dynamics", "CAPS — theory manual", f'{icon("chevr", 12, DIM)}<span>Theory manual</span>', content,
                 ("<span>Every equation in CAPS has a page like this</span>", "<span>Manual version matches CAPS 0.1.0</span>"), ["Theory manual"])


# ------------------------------------------------------------------ Project home
def project_home():
    def thumb(sid, atoms, w=210, h=120):
        sc = Scene(sid, w, h, yaw=0.5, pitch=0.4, persp=0.2, fog=0.5, outline=False, bond_w=0.22)
        sc.add_atoms(atoms)
        return sc.svg()
    docs = []
    for i, (n, t, seed) in enumerate((("PS_isotactic_DP40", "isotactic", 1), ("PS_syndiotactic_DP40", "syndiotactic", 1), ("PS_atactic_DP40", "atactic", 5))):
        a, _, _ = polystyrene(4, t, seed=seed, hydrogens=False, curve=False)
        docs.append(f'<a href="#" style="display: flex; flex-direction: column; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px; overflow: hidden; text-decoration: none; color: {TEXT}">'
                    f'<div style="background: {BG0}; border-bottom: 1px solid {LINE}">{thumb("ph" + str(i), a)}</div>'
                    f'<span style="padding: 10px 12px; display: flex; flex-direction: column; gap: 3px"><span style="font-size: 12.5px; font-weight: 500">{n}.caps</span><span style="font-size: 11px; color: {DIM}">{t} · 642 atoms per chain</span></span></a>')
    res = table(["Cell", "Density 298 K (g/cm³)", "T<sub>g</sub> (K)", "C<sub>∞</sub>", "Status"],
                [["isotactic · 20 × 40", "[result]", "[result]", "[result]", f'<span style="color: {OK}">analysed</span>'],
                 ["syndiotactic · 20 × 40", "[result]", "[result]", "[result]", f'<span style="color: {ACC}">equilibrating</span>'],
                 ["atactic · 20 × 40", "[result]", "[result]", "[result]", f'<span style="color: {DIM}">queued</span>'],
                 ["experiment (atactic PS)", "1.04–1.065", "≈ 373", "≈ 9.5–10", f'<span style="color: {MUTED}">literature</span>']],
                ["28%", "20%", "14%", "14%", "24%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=12, hl={3})
    methods = ("Amorphous cells of polystyrene (20 chains, degree of polymerisation 40, 12 840 atoms per cell) were built with CAPS 0.1.0. "
               "Atom types and charges were assigned from GAFF2 [1]. Chains were grown by configurational-bias Monte Carlo with 16 trial directions per step, "
               "using the full force-field energy [2, 3], at an initial density of 0.50 g/cm³, and minimised with L-BFGS [4] after a soft-core push-off [5]. "
               "Cells were equilibrated with the 21-step compression–decompression protocol (P<sub>max</sub> = 50 000 bar) [6], followed by 10 ns of NPT dynamics "
               "at 298 K and 1 atm using stochastic velocity rescaling [7] and stochastic cell rescaling [8]. Electrostatics used smooth PME with a relative tolerance of 10⁻⁵ [9]; "
               "Lennard-Jones interactions were truncated at 12 Å with analytic tail corrections.")
    refs = ["Wang et al., J. Comput. Chem. 25, 1157 (2004)", "Siepmann &amp; Frenkel, Mol. Phys. 75, 59 (1992)", "Rosenbluth &amp; Rosenbluth, J. Chem. Phys. 23, 356 (1955)",
            "Liu &amp; Nocedal, Math. Program. 45, 503 (1989)", "Auhl et al., J. Chem. Phys. 119, 12718 (2003)", "Larsen et al., Macromolecules 44, 6944 (2011)",
            "Bussi et al., J. Chem. Phys. 126, 014101 (2007)", "Bernetti &amp; Bussi, J. Chem. Phys. 153, 114107 (2020)", "Essmann et al., J. Chem. Phys. 103, 8577 (1995)"]
    rl = "".join(f'<div style="font-size: 11.5px; color: {MUTED}; line-height: 1.5">[{i}] {r}</div>' for i, r in enumerate(refs, 1))
    content = (pagehead("PS-tacticity study", "Project · 3 documents · 7 jobs · last activity today",
                        row(btn("Share project", ic="link"), btn("New document", True, "plus"), gap=8))
               + '<div style="flex-grow: 1; display: flex; gap: 16px; padding: 18px 22px; min-height: 0">'
               + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px">'
               + f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px">{"".join(docs)}</div>'
               + card("Results", col(res, f'<div style="font-size: 11.5px; color: {DIM}">Experimental ranges from the Polymer Handbook; C∞ from Flory 1969 and Mattice &amp; Suter 1994.</div>', gap=8), btn("Export table", ic="download", small=True), 12)
               + '</div>'
               + f'<div style="width: 470px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
               + card("Methods section · generated from provenance", col(f'<p style="margin: 0; font-size: 12.5px; line-height: 1.6">{methods}</p>', rl, gap=10),
                      row(btn("Copy", ic="copy", small=True), btn("BibTeX", ic="download", small=True), gap=6), 14, "flex-grow: 1")
               + '</div></div>')
    return shell("Studio", "CAPS — project home", f'{icon("chevr", 12, DIM)}<span>Project</span>', content,
                 ("<span>~/research/ps-tacticity</span>", "<span>Methods text updates when provenance changes</span>"), ["PS-tacticity study"])


# ------------------------------------------------------------------ Parameter sweep
def sweep():
    tacts = ["isotactic", "syndiotactic", "atactic"]
    dps = [20, 40, 80]
    seeds = [1, 2, 3]
    rnd = random.Random(4)
    status_col = {"done": OK, "running": ACC, "queued": DIM, "failed": ERR}
    grid = []
    n = 0
    for t in tacts:
        cells = []
        for dp in dps:
            dots_ = []
            for s in seeds:
                n += 1
                st = "done" if n <= 11 else ("running" if n <= 15 else "queued")
                if n == 9:
                    st = "failed"
                dots_.append(f'<span title="{t} DP{dp} seed {s}: {st}" style="width: 30px; height: 30px; border-radius: 6px; background: {status_col[st]}; opacity: {1 if st != "queued" else 0.45}; display: flex; align-items: center; justify-content: center; font-family: {MONO}; font-size: 11px; color: {ACC_INK}">{s}</span>')
            cells.append(f'<td style="padding: 10px; border-bottom: 1px solid {BG2}"><div style="display: flex; gap: 6px">{"".join(dots_)}</div></td>')
        grid.append(f'<tr><th scope="row" style="text-align: left; padding: 10px; font-size: 12.5px; font-weight: 500; border-bottom: 1px solid {BG2}">{t}</th>{"".join(cells)}</tr>')
    gt = (f'<table style="border-collapse: collapse"><thead><tr><th></th>' + "".join(f'<th scope="col" style="text-align: left; padding: 8px 10px; font-size: 11px; font-weight: 600; letter-spacing: 0.06em; text-transform: uppercase; color: {DIM}">DP {d}</th>' for d in dps)
          + f'</tr></thead><tbody>{"".join(grid)}</tbody></table>')
    legend = row(*[chip(f'{dot(c)} {k}') for k, c in status_col.items()], gap=6)
    atoms_per = {20: 20 * 16 + 2, 40: 40 * 16 + 2, 80: 80 * 16 + 2}
    est = table(["DP", "Atoms per cell (20 chains)", "Runs"], [[str(d), f"{20 * atoms_per[d]:,}".replace(",", " "), "9"] for d in dps], ["20%", "55%", "25%"], mono_cols=(0, 1, 2), align_right=(1, 2), fs=12)
    content = (pagehead("Sweep · tacticity × chain length", "27 runs = 3 tacticities × 3 DP × 3 seeds · recipe build_ps.py → Grow → Equilibrate → Analyze",
                        row(chip("10 done · 1 failed · 4 running · 12 queued", TEXT, BG3), btn("Pause sweep", ic="pause"), btn("Add runs", True, "plus"), gap=8))
               + '<div style="flex-grow: 1; display: flex; gap: 16px; padding: 18px 22px; min-height: 0">'
               + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px">'
               + card("Runs", col(gt, legend, f'<div style="font-size: 12px; color: {MUTED}">Each square is one run; the number is its seed. Click a square to open its job, cell and provenance.</div>', gap=12))
               + card("Results by condition", col(table(["Condition", "Density (g/cm³)", "T<sub>g</sub> (K)", "Seeds done"],
                                                        [["isotactic · DP 20", "[mean ± sd]", "[mean ± sd]", "3 / 3"], ["isotactic · DP 40", "[mean ± sd]", "[mean ± sd]", "3 / 3"],
                                                         ["isotactic · DP 80", "[mean ± sd]", "[mean ± sd]", "2 / 3 · 1 failed"], ["syndiotactic · DP 20", "[mean ± sd]", "[mean ± sd]", "2 / 3"]],
                                                        mono_cols=(1, 2), align_right=(1, 2), fs=12), gap=8), chip("mean ± sd over seeds"), 12, "flex-grow: 1")
               + '</div>'
               + f'<div style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
               + card("Parameters", col(row(field("tacticity", "iso, syndio, atactic"), gap=8), row(field("dp", "20, 40, 80"), field("seed", "1, 2, 3"), gap=8), select("Combine as", "Full grid (27 runs)"), gap=10))
               + card("Size", col(est, f'<div style="font-size: 11.5px; color: {DIM}">Atom counts: 16 per styrene unit plus 2 end hydrogens per chain.</div>', gap=8))
               + card("Where", col(select("Host", "hpc-login2 · SLURM gpu"), select("Concurrency", "4 runs at a time"), gap=10))
               + '</div></div>')
    return shell("Jobs", "CAPS — parameter sweep", f'{icon("chevr", 12, DIM)}<span>Sweep</span>', content,
                 ("<span>sweep-2 · 27 runs</span>", "<span>Every run keeps its own provenance manifest</span>"), ["PS-tacticity study"])
