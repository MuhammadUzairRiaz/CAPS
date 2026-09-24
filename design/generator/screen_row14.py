import math, itertools
from lib import *
from mols import *
from screen_app import card, pagehead, shell
from screen_studio2 import spacer
from screen_row7 import with_overlay, modal
from screen_row10 import logo_mark, heading
import screen_main

# ---------------------------------------------------------------- colour science
MACHADO = {  # Machado, Oliveira & Fernandes, IEEE TVCG 15, 1291 (2009), severity 1.0, linear RGB
    "Protanopia": [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
    "Deuteranopia": [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
    "Tritanopia": [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
}


def to_lin(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def to_srgb(c):
    c = min(1.0, max(0.0, c))
    return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


def hex2rgb(h):
    return [int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)]


def rgb2hex(r):
    return "#%02X%02X%02X" % tuple(round(x * 255) for x in r)


def simulate(h, kind):
    lin = [to_lin(c) for c in hex2rgb(h)]
    m = MACHADO[kind]
    out = [sum(m[i][j] * lin[j] for j in range(3)) for i in range(3)]
    return rgb2hex([to_srgb(c) for c in out])


def lab(h):
    r, g, b = [to_lin(c) for c in hex2rgb(h)]
    X = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047
    Y = (0.2126 * r + 0.7152 * g + 0.0722 * b) / 1.0
    Z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883
    f = lambda t: t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116
    fx, fy, fz = f(X), f(Y), f(Z)
    return (116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz))


def dE(a, b):
    la, lb = lab(a), lab(b)
    return math.sqrt(sum((la[i] - lb[i]) ** 2 for i in range(3)))


# ---------------------------------------------------------------- Colour-blind check
def colourblind():
    palettes = [("Elements", [("C", ELEM["C"][0]), ("H", ELEM["H"][0]), ("O", ELEM["O"][0]), ("N", ELEM["N"][0]), ("Si", ELEM["Si"][0]), ("S", ELEM["S"][0]), ("Na", ELEM["Na"][0]), ("Cl", ELEM["Cl"][0])]),
                ("Chains (Grow)", [("1", "#F0A83C"), ("2", "#6CC4D8"), ("3", "#E07A5F"), ("4", "#9B7BD6"), ("5", "#7CC784"), ("6", "#D6A45E"), ("7", "#E9ECEF"), ("8", "#4C7BD9")]),
                ("Status", [("ok", OK), ("running", ACC), ("error", ERR), ("select", SEL)])]
    kinds = ["Normal", "Protanopia", "Deuteranopia", "Tritanopia"]
    THR = 12.0
    blocks = []
    flagged_all = []
    for pname, items in palettes:
        rows = []
        for k in kinds:
            sw = []
            cols = [(n, c if k == "Normal" else simulate(c, k)) for n, c in items]
            for n, c in cols:
                txt = "#0F1113" if lab(c)[0] > 55 else "#FFFFFF"
                sw.append(f'<span title="{n}" style="width: 44px; height: 34px; border-radius: 5px; background: {c}; display: flex; align-items: center; justify-content: center; font-family: {MONO}; font-size: 11px; color: {txt}">{n}</span>')
            if k != "Normal":
                for (n1, c1), (n2, c2) in itertools.combinations(cols, 2):
                    d = dE(c1, c2)
                    if d < THR:
                        flagged_all.append((pname, k, n1, n2, d))
            rows.append(f'<div style="display: flex; align-items: center; gap: 10px"><span style="width: 96px; flex-shrink: 0; font-size: 11.5px; color: {MUTED}">{k}</span><div style="display: flex; gap: 4px">{"".join(sw)}</div></div>')
        blocks.append(card(pname, col(*rows, gap=6), None or chip(f"{len(items)} colours"), 14))
    fl = sorted(flagged_all, key=lambda t: t[4])
    frows = [[p, k, f"{a} ↔ {b}", f"{d:.1f}"] for p, k, a, b, d in fl[:12]]
    ft = table(["Palette", "Vision", "Pair", "ΔE*ab"], frows or [["—", "—", "none below threshold", "—"]], ["30%", "26%", "26%", "18%"], mono_cols=(2, 3), align_right=(3,), fs=12)
    fixes = [("Chains", "When pairs collide, CAPS adds a dash pattern per chain and labels chain ends; colour is never the only cue."),
             ("Status", "Every status also has an icon and a word (done, running, failed)."),
             ("Elements", "Labels on hover and in the inspector; the Okabe–Ito palette is one click away in Settings.")]
    fx = "".join(f'<div style="display: flex; gap: 10px; font-size: 12.5px; line-height: 1.5"><span style="width: 70px; flex-shrink: 0; color: {DIM}">{a}</span><span>{b}</span></div>' for a, b in fixes)
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 18px; padding: 30px 40px">'
            + f'<div style="display: flex; align-items: flex-end; gap: 20px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">Colour-vision check</h1>'
              f'<span style="font-size: 13px; color: {MUTED}">Simulated with Machado, Oliveira &amp; Fernandes (2009), severity 1.0 · pairs flagged when ΔE*ab &lt; {THR:.0f}</span></div>'
            + '<div style="display: flex; gap: 20px; min-height: 0">'
            + f'<div style="width: 520px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">{"".join(blocks)}</div>'
            + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px">'
            + card(f"Pairs that become hard to tell apart · {len(fl)} found", col(ft, f'<div style="font-size: 11.5px; color: {DIM}">Lowest 12 shown. Computed from the palette hex values on this canvas.</div>', gap=8), None or "", 14)
            + card("How the UI compensates", col(fx, gap=10))
            + '</div></div></div>')
    return page("CAPS — colour-vision check", body), len(fl)


# ---------------------------------------------------------------- Motion
def cubic_bezier(p1x, p1y, p2x, p2y, n=60):
    pts = []
    for i in range(n + 1):
        t = i / n
        x = 3 * (1 - t) ** 2 * t * p1x + 3 * (1 - t) * t ** 2 * p2x + t ** 3
        y = 3 * (1 - t) ** 2 * t * p1y + 3 * (1 - t) * t ** 2 * p2y + t ** 3
        pts.append((x, y))
    return pts


def motion():
    curves = [("standard", (0.2, 0, 0, 1), ACC), ("enter", (0, 0, 0, 1), SEL), ("exit", (0.3, 0, 1, 1), "#E07A5F")]
    cp = "".join(f'<div style="display: flex; flex-direction: column; gap: 6px">{plot(144, 150, [(cubic_bezier(*v), c, 2.2, None), ([(0, 0), (1, 1)], BG3, 1, "3 3")], (0, 1), (0, 1), [0, 1], [0, 1], "time", "", pad=(26, 10, 10, 26))}'
                 f'<span style="font-size: 12.5px; font-weight: 600">{n}</span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">cubic-bezier{v}</span></div>' for n, v, c in curves)
    tokens = [("instant", "0 ms", "selection highlight, hover"), ("fast", "120 ms", "chips, toggles, tooltips"), ("base", "200 ms", "drawers, menus, toasts"),
              ("slow", "320 ms", "modal dialogs, page changes"), ("camera", "450 ms", "3D fly-to selection, reset view")]
    tt = table(["Token", "Duration", "Used for"], [[f'<span style="font-family: {MONO}">{a}</span>', b, c] for a, b, c in tokens], ["22%", "20%", "58%"], mono_cols=(1,), fs=12.5)
    atoms, backbone, stereo = polystyrene(5, "atactic", seed=5)
    frames = []
    for i, (sc_k, yaw, lbl) in enumerate(((1.0, -0.28, "0 ms"), (1.5, -0.1, "150 ms"), (2.1, 0.05, "300 ms"), (2.4, 0.12, "450 ms · settled"))):
        sc = Scene(f"mo{i}", 205, 140, yaw=yaw, pitch=0.42, persp=0.2, fog=0.5, outline=True)
        sc.add_atoms(atoms)
        sc._prep()
        sc.scale *= sc_k
        ring = [j for j, a in enumerate(atoms) if a.get("ar")][:6]
        cx = sum(atoms[j]["p"][0] for j in ring) / 6
        cy = sum(atoms[j]["p"][1] for j in ring) / 6
        cz = sum(atoms[j]["p"][2] for j in ring) / 6
        f = min(1.0, (sc_k - 1) / 1.4)
        sc.center = tuple(sc.center[k] * (1 - f) + (cx, cy, cz)[k] * f for k in range(3))
        sc.halo = set(ring)
        frames.append(f'<div style="display: flex; flex-direction: column; gap: 6px"><div style="background: {BG0}; border: 1px solid {LINE}; border-radius: 8px; overflow: hidden">{sc.svg()}</div><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{lbl}</span></div>')
    rules = [("Reduce motion", "camera moves become cuts; progress bars stay; nothing loops"), ("Never animate data", "plots, energies and trajectories are never tweened between values"),
             ("One thing at a time", "a drawer and a camera move never run together"), ("Long jobs", "progress is continuous, no spinners over 2 s; completion is a single toast")]
    rl = "".join(f'<div style="display: flex; gap: 12px; padding: 8px 0; border-bottom: 1px solid {BG2}; font-size: 12.5px"><span style="width: 150px; flex-shrink: 0; font-weight: 600">{a}</span><span style="color: {MUTED}">{b}</span></div>' for a, b in rules)
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 20px; padding: 30px 40px">'
            + f'<div style="display: flex; align-items: flex-end; gap: 20px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">Motion</h1><span style="font-size: 13px; color: {MUTED}">Short, purposeful, and off when the system asks for reduced motion</span></div>'
            + '<div style="display: flex; gap: 28px">'
            + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 18px">'
            + col(heading("Camera fly-to selection · 450 ms, standard easing"), f'<div style="display: flex; gap: 12px">{"".join(frames)}</div>', gap=10)
            + col(heading("Rules"), f'<div>{rl}</div>', gap=6) + '</div>'
            + f'<div style="width: 460px; flex-shrink: 0; display: flex; flex-direction: column; gap: 18px">'
            + col(heading("Durations"), tt, gap=8) + col(heading("Easing"), f'<div style="display: flex; gap: 14px">{cp}</div>', gap=8) + '</div>'
            + '</div></div>')
    return page("CAPS — motion", body)


# ---------------------------------------------------------------- Accessibility map
def a11y():
    base = screen_main.build()
    marks = [(1, 36, 120, "Module rail"), (2, 300, 66, "Document tabs"), (3, 520, 68, "Tools"), (4, 200, 200, "Project &amp; fragments"),
             (5, 700, 400, "3D view"), (6, 1280, 200, "Inspector"), (7, 700, 760, "History &amp; console")]
    ov = "".join(f'<div style="position: absolute; left: {x - 15}px; top: {y - 15}px; display: flex; align-items: center; gap: 6px"><span style="width: 30px; height: 30px; border-radius: 50%; background: {SEL}; color: #0F1113; font-weight: 700; font-size: 14px; display: flex; align-items: center; justify-content: center; box-shadow: 0 0 0 3px #0F1113">{n}</span>'
                 f'<span style="padding: 3px 8px; background: #0F1113E6; border: 1px solid {SEL}; border-radius: 5px; font-size: 12px; color: {TEXT}; white-space: nowrap">{t}</span></div>' for n, x, y, t in marks)
    focus = f'<div style="position: absolute; left: 332px; top: 92px; width: 788px; height: 596px; border: 2px solid {SEL}; border-radius: 4px; box-shadow: 0 0 0 4px #6CC4D840"></div>'
    atoms, backbone, stereo = polystyrene(7, "atactic", seed=5)
    ang = angle(atoms[backbone[5]]["p"], atoms[backbone[6]]["p"], atoms[backbone[7]]["p"])
    sr = [("Focus 3D view", f"3D view, polystyrene, {len(atoms)} atoms, 7 selected. Arrow keys walk bonds; Tab moves to the next molecule."),
          ("→ on C12", "C12, carbon, sp3, type c3, charge minus 0.096, stereo R. Bonded to C11, C13, C28 and H44."),
          ("M then 3 atoms", f"Angle C5, C6, C7: {ang:.1f} degrees. Pinned as monitor 2."),
          ("Job finishes", "Polite announcement: Growth finished, 20 chains, no close contacts.")]
    srl = "".join(f'<div style="display: flex; flex-direction: column; gap: 3px; padding: 8px 0; border-bottom: 1px solid {BG3}"><span style="font-family: {MONO}; font-size: 11px; color: {SEL}">{a}</span><span style="font-size: 12.5px; line-height: 1.45">“{b}”</span></div>' for a, b in sr)
    panel = (f'<div style="position: absolute; left: 350px; top: 360px; width: 460px; padding: 16px; background: {BG1}; border: 1px solid {SEL}; border-radius: 12px; box-shadow: 0 20px 50px #000000A0">'
             f'<div style="display: flex; align-items: center; gap: 8px; margin-bottom: 6px">{icon("eye", 16, SEL)}<span style="font-size: 13.5px; font-weight: 600">What a screen reader hears</span></div>{srl}'
             f'<div style="font-size: 11.5px; color: {DIM}; margin-top: 8px">Tab order follows the numbers. Every icon button has an aria-label; the 3D view is one focus stop with its own keyboard model.</div></div>')
    dim = '<div style="position: absolute; inset: 0; background: #07080A73"></div>'
    return with_overlay(base, dim + focus + ov + panel).replace("<title>CAPS Studio — builder</title>", "<title>CAPS Studio — accessibility map</title>")


# ---------------------------------------------------------------- System states
def states():
    def st(icn, colr, title, msg, actions, tag):
        acts = "".join(actions)
        return (f'<section style="display: flex; flex-direction: column; gap: 10px; padding: 18px; background: {BG1}; border: 1px solid {LINE}; border-radius: 12px">'
                f'<div style="display: flex; align-items: center; gap: 10px"><span style="width: 34px; height: 34px; border-radius: 8px; background: {BG3}; display: flex; align-items: center; justify-content: center">{icon(icn, 18, colr)}</span>'
                f'<h2 style="margin: 0; font-size: 14px; font-weight: 600; flex-grow: 1">{title}</h2>{chip(tag, colr)}</div>'
                f'<p style="margin: 0; font-size: 12.5px; line-height: 1.55; color: {MUTED}">{msg}</p><div style="display: flex; gap: 8px; margin-top: auto">{acts}</div></section>')
    cards = [
        st("cpu", MUTED, "No supported GPU", "CAPS runs every engine on the CPU. Packing and MD are slower; results are identical in deterministic mode.", [btn("Why?", small=True), btn("Continue", True, small=True)], "info"),
        st("server", WARN, "hpc-login2 unreachable", "The SSH tunnel timed out after 20 s. Running jobs continue on the cluster; CAPS will reconnect and fetch results.", [btn("Retry now", True, small=True), btn("Details", small=True)], "warning"),
        st("tag", ERR, "3 missing parameters", "Dynamics cannot start: 2 dihedrals and 1 improper have no GAFF2 parameters. CAPS never guesses them.", [btn("Open Field report", True, small=True, href="ForceField.dc.html"), btn("Import .frcmod", small=True)], "blocks run"),
        st("link", SEL, "Bonds were inferred", "PE_melt.xyz had no bond records. 2 998 bonds were perceived from covalent radii; review before typing.", [btn("Review bonds", True, small=True), btn("Undo import", small=True)], "check"),
        st("download", ACC, "Opening a large trajectory", "traj.dcd · 1 024 000 atoms · 4 000 frames. Frames stream from disk; the first frame is ready now.", [btn("Cancel", small=True)], "streaming"),
        st("file", MUTED, "Offline", "The theory manual and fragment library are cached locally. Remote hosts and update checks wait until you are back online.", [btn("OK", small=True)], "info"),
    ]
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 18px; padding: 30px 40px">'
            + f'<div style="display: flex; align-items: flex-end; gap: 20px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">System states</h1><span style="font-size: 13px; color: {MUTED}">Each message says what happened, what still works, and what to do next</span></div>'
            + f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 16px">{"".join(cards)}</div>'
            + f'<div style="display: flex; gap: 16px">{card("Writing rules", col(*[kv(a, b, False) for a, b in (("Lead with", "what happened, in plain words"), ("Then", "what still works"), ("Then", "one primary action"), ("Never", "blame the user, or say an error is unknown when a log exists"))], gap=8), extra="flex-grow: 1")}'
              f'{card("Severity", col(*[f"<div style=&#34;display: flex; gap: 10px; align-items: center; font-size: 12.5px&#34;>{dot(c)}<span style=&#34;width: 80px&#34;>{a}</span><span style=&#34;color: {MUTED}&#34;>{b}</span></div>".replace("&#34;", chr(34)) for a, b, c in (("info", "nothing to do", MUTED), ("check", "worth a look", SEL), ("warning", "may affect work", WARN), ("blocks run", "must be fixed first", ERR))], gap=8), extra="width: 460px; flex-shrink: 0")}</div>'
            + '</div>')
    return page("CAPS — system states", body)


# ---------------------------------------------------------------- Update dialog
def update_dialog():
    base = screen_main.build()
    head = (f'<div style="display: flex; align-items: center; gap: 12px; padding: 18px 22px; border-bottom: 1px solid {LINE}">{logo_mark(36, TEXT, ACC)}'
            f'<div style="display: flex; flex-direction: column; gap: 2px"><h2 style="margin: 0; font-size: 16px; font-weight: 600">CAPS 0.2.0 is available</h2><span style="font-size: 12px; color: {MUTED}">You have 0.1.0 · example release notes</span></div>{spacer()}{tbtn("close", "Close")}</div>')
    new = ["Force-field-aware growth for ring backbones (PET, PC)", "Cylinder and sphere regions in Pack", "Trajectory export to XTC"]
    fixed = ["Crystal builder kept the wrong setting for Pnam after undo", "Compact layout clipped the inspector on 1280 px screens"]
    alters = [("Smooth PME default tolerance", "1e−5 → 1e−6", "energies change in the 5th significant figure"), ("GAFF2 parameter database", "caps-ffdb 1.0 → 1.1", "3 dihedral terms corrected against the source")]

    def lst(items, colr):
        return "".join(f'<div style="display: flex; gap: 8px; font-size: 12.5px; line-height: 1.5">{dot(colr, 6)}<span>{t}</span></div>' for t in items)
    alt = "".join(f'<div style="display: flex; flex-direction: column; gap: 3px; padding: 10px 12px; background: #2A2314; border: 1px solid #5A4520; border-radius: 7px">'
                  f'<div style="display: flex; justify-content: space-between; gap: 10px"><span style="font-size: 12.5px; font-weight: 600">{a}</span><span style="font-family: {MONO}; font-size: 12px; color: {ACC}">{b}</span></div>'
                  f'<span style="font-size: 12px; color: {MUTED}">{c}</span></div>' for a, b, c in alters)
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 16px; padding: 18px 22px">'
            + col(f'<h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {WARN}">Changes that can alter results</h3>', alt, gap=8)
            + '<div style="display: flex; gap: 24px">'
            + col(f'<h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">New</h3>', lst(new, OK), gap=6, extra="flex-grow: 1")
            + col(f'<h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Fixed</h3>', lst(fixed, SEL), gap=6, extra="flex-grow: 1")
            + '</div>'
            + f'<div style="padding: 12px 14px; background: {BG0}; border: 1px solid {LINE}; border-radius: 8px; display: flex; flex-direction: column; gap: 10px">'
            + toggle("Keep existing projects on the 0.1.0 engines (recommended)", True)
            + f'<span style="font-size: 12px; color: {MUTED}; line-height: 1.5">Open projects keep reproducing their published numbers. New projects use 0.2.0. You can re-run a project on 0.2.0 and compare both manifests side by side.</span></div>'
            + '</div>')
    foot = f'<div style="display: flex; justify-content: flex-end; gap: 8px; padding: 12px 22px; border-top: 1px solid {LINE}">{btn("Later")}{btn("Full changelog", ic="file")}{btn("Install 0.2.0", True, "download")}</div>'
    return with_overlay(base, modal(head + body + foot, 720, 640)).replace("<title>CAPS Studio — builder</title>", "<title>CAPS — update available</title>")


# ---------------------------------------------------------------- Download page (website)
def download_page():
    W, H = 1440, 1700
    atoms, _, _ = polystyrene(8, "atactic", seed=5)
    sc = Scene("dl", 620, 380, yaw=-0.28, pitch=0.42, roll=0.05, persp=0.25, fog=0.55)
    sc.add_atoms(atoms)
    hero_img = sc.svg()
    nav = (f'<header style="height: 72px; flex-shrink: 0; display: flex; align-items: center; gap: 28px; padding: 0 64px; border-bottom: 1px solid {LINE}">'
           f'<a href="#" style="display: flex; align-items: center; gap: 10px; text-decoration: none; color: {TEXT}">{logo_mark(30, TEXT, ACC)}<span style="font-size: 18px; font-weight: 700; letter-spacing: 0.04em">CAPS</span></a>'
           + "".join(f'<a href="#" style="font-size: 14px; color: {MUTED}; text-decoration: none">{t}</a>' for t in ("Features", "Documentation", "Benchmarks", "Cite", "Source"))
           + f'<div style="flex-grow: 1"></div>{btn("Download", True, "download")}</header>')
    hero = (f'<section style="display: flex; align-items: center; gap: 48px; padding: 64px">'
            f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 20px">'
            f'<span style="font-family: {MONO}; font-size: 13px; color: {ACC}">Chain Assembly and Packing Suite · 0.1.0</span>'
            f'<h1 style="margin: 0; font-size: 50px; line-height: 1.08; font-weight: 600; letter-spacing: -0.02em">Build, pack and equilibrate polymer systems, with every step cited.</h1>'
            f'<p style="margin: 0; font-size: 17px; line-height: 1.6; color: {MUTED}; max-width: 560px">One self-contained program from monomer SMILES to an equilibrated cell. Its own engines for growth, packing, minimisation and dynamics, and input decks for LAMMPS and GROMACS that reproduce its energies.</p>'
            f'<div style="display: flex; gap: 12px">{btn("Download for macOS · Apple silicon", True, "download")}{btn("Other platforms", ic="chev")}</div>'
            f'<span style="font-size: 12.5px; color: {DIM}">BSD-3-Clause · Linux, macOS and Windows · [release date]</span></div>'
            f'<div style="width: 620px; flex-shrink: 0; background: {BG1}; border: 1px solid {LINE}; border-radius: 14px; overflow: hidden">{hero_img}</div></section>')
    plats = [("macOS", "Apple silicon and Intel · notarised .app"), ("Windows", "x64 · signed MSI installer"), ("Linux", "AppImage and Flatpak"), ("Python", "wheels for 3.12+ on all three")]
    pc = "".join(f'<div style="display: flex; flex-direction: column; gap: 6px; padding: 18px; background: {BG1}; border: 1px solid {LINE}; border-radius: 12px"><span style="font-size: 15px; font-weight: 600">{a}</span><span style="font-size: 13px; color: {MUTED}">{b}</span></div>' for a, b in plats)
    cmds = (f'<div style="display: flex; flex-direction: column; gap: 8px; padding: 18px 20px; background: #0A0B0C; border: 1px solid {LINE}; border-radius: 12px; font-family: {MONO}; font-size: 14px">'
            f'<div><span style="color: {DIM}"># conda-forge</span></div><div><span style="color: {ACC}">$</span> conda install -c conda-forge [package-name]</div>'
            f'<div style="margin-top: 6px"><span style="color: {DIM}"># PyPI</span></div><div><span style="color: {ACC}">$</span> pip install [package-name]</div></div>')
    feats = [("grow", "Force-field-aware growth", "Configurational-bias growth with the real force field, plus Theodorou–Suter with RIS statistics."),
             ("pack", "Packing without packmol", "Rigid-body overlap minimisation with region constraints and a guaranteed minimum distance."),
             ("dyn", "Its own MD engine", "Verlet and r-RESPA, CSVR and Nosé–Hoover, cell rescaling, SPME, LINCS."),
             ("tag", "Typing you can explain", "Every atom type comes with the rule that assigned it. Missing parameters stop the run."),
             ("jobs", "Provenance by default", "Every file carries its inputs, versions, seeds, approximations and citations."),
             ("bench", "Benchmarks you can rerun", "One command regenerates the validation tables, including where CAPS is slower.")]
    fc = "".join(f'<div style="display: flex; flex-direction: column; gap: 10px; padding: 22px; background: {BG1}; border: 1px solid {LINE}; border-radius: 12px">{icon(i, 24, ACC, 1.5)}<span style="font-size: 16px; font-weight: 600">{t}</span><span style="font-size: 13.5px; line-height: 1.55; color: {MUTED}">{d}</span></div>' for i, t, d in feats)
    cite_ = (f'<section style="display: flex; gap: 32px; padding: 0 64px 64px">'
             f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 12px"><h2 style="margin: 0; font-size: 26px; font-weight: 600">How to cite</h2>'
             f'<p style="margin: 0; font-size: 14px; color: {MUTED}; line-height: 1.6">If CAPS contributes to your work, cite the software paper and the methods it used. CAPS writes the method references for you into every provenance manifest.</p>'
             f'<div style="padding: 16px 18px; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px; font-family: {MONO}; font-size: 13px; line-height: 1.7; color: {TEXT}">[Authors]. CAPS: Chain Assembly and Packing Suite. [Journal] ([year]). doi:[pending]</div></div>'
             f'<div style="width: 420px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px"><h2 style="margin: 0; font-size: 26px; font-weight: 600">Coming from scripts?</h2>'
             f'<p style="margin: 0; font-size: 14px; color: {MUTED}; line-height: 1.6">CAPS reads YAML and JSON recipes. Steps that used packmol, LAMMPS or DL_FIELD now run on CAPS engines.</p>{btn("Migration guide", ic="file")}</div></section>')
    footer_ = (f'<footer style="margin-top: auto; display: flex; align-items: center; gap: 24px; padding: 24px 64px; border-top: 1px solid {LINE}; font-size: 13px; color: {DIM}">'
               f'<span>CAPS · BSD-3-Clause</span><a href="#" style="color: {MUTED}; text-decoration: none">Source code [repository]</a><a href="#" style="color: {MUTED}; text-decoration: none">Documentation</a><a href="#" style="color: {MUTED}; text-decoration: none">Theory manual</a></footer>')
    body = (nav + hero
            + f'<section style="display: flex; flex-direction: column; gap: 18px; padding: 0 64px 56px"><h2 style="margin: 0; font-size: 26px; font-weight: 600">Download</h2>'
              f'<div style="display: flex; gap: 20px"><div style="flex-grow: 1; display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 12px">{pc}</div><div style="width: 520px; flex-shrink: 0">{cmds}</div></div></section>'
            + f'<section style="display: flex; flex-direction: column; gap: 18px; padding: 0 64px 56px"><h2 style="margin: 0; font-size: 26px; font-weight: 600">What it does</h2><div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 14px">{fc}</div></section>'
            + cite_ + footer_)
    return page("CAPS — download", body, W, H)
