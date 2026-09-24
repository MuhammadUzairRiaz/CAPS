import math, random
from lib import *
from mols import *


def shell(module, title, crumb, content, status, tabs_=None):
    body = (topbar(tabs_ or ["PS_atactic_DP40.caps"], (tabs_ or ["PS_atactic_DP40.caps"])[0], crumb)
            + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail(module)
            + f'<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">{content}</div></div>'
            + statusbar(*status))
    return page(title, body)


def pagehead(title, sub, right=""):
    return (f'<div style="flex-shrink: 0; display: flex; align-items: center; gap: 16px; padding: 16px 22px; background: {BG1}; border-bottom: 1px solid {LINE}">'
            f'<div style="display: flex; flex-direction: column; gap: 3px"><h1 style="margin: 0; font-size: 19px; font-weight: 600; letter-spacing: -0.01em">{title}</h1>'
            f'<div style="font-size: 12.5px; color: {MUTED}">{sub}</div></div><div style="flex-grow: 1"></div>{right}</div>')


def card(title, body, right="", pad=16, extra=""):
    return (f'<section style="display: flex; flex-direction: column; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px; overflow: hidden; {extra}">'
            f'<div style="height: 42px; flex-shrink: 0; display: flex; align-items: center; justify-content: space-between; gap: 8px; padding: 0 {pad}px; border-bottom: 1px solid {LINE}">'
            f'<h2 style="margin: 0; font-size: 13px; font-weight: 600">{title}</h2>{right}</div>'
            f'<div style="padding: {pad}px; display: flex; flex-direction: column; gap: 12px; flex-grow: 1; min-height: 0">{body}</div></section>')


def progress(frac, color=ACC, h=6):
    return f'<div style="height: {h}px; border-radius: {h / 2}px; background: {BG3}"><div style="width: {frac * 100:.0f}%; height: {h}px; border-radius: {h / 2}px; background: {color}"></div></div>'


# ------------------------------------------------------------------ Grow
def grow(style="Backbone"):
    from display import ps_cell, render_cell, count, style_bar, zoom_inset
    Lb = 38.0
    sc = Scene("g", 740, 450, yaw=0.6, pitch=0.42, persp=0.3, fog=0.55, outline=True)
    cols = ["#F0A83C", "#6CC4D8", "#E07A5F", "#9B7BD6", "#7CC784", "#D6A45E", "#E9ECEF", "#4C7BD9", "#C77DBA", "#8FB8A8"]
    for k in range(12):
        rnd = random.Random(100 + k)
        st = (rnd.uniform(4, Lb - 4), rnd.uniform(4, Lb - 4), rnd.uniform(4, Lb - 4))
        n = 60 if k < 11 else 26
        pts = random_chain(st, n, 200 + k, bond=1.54, box=((0.8, 0.8, 0.8), (Lb - 0.8, Lb - 0.8, Lb - 0.8)), persistence=0.6)
        sc.add_tube(pts, cols[k % len(cols)], 0.55)
        if k == 11:
            sc.add_atoms([{"e": "C", "p": pts[-1], "r": 0.9, "c": ACC}], [])
            sc.halo = {len(sc.atoms) - 1}
    sc.add_box((0, 0, 0), (Lb, 0, 0), (0, Lb, 0), (0, 0, Lb), MUTED, 1.2, None, 0.8)
    if style == "All atoms":
        Lz = 34.0
        Lz = 34.0
        cell = ps_cell(14, 6, Lz, seed=12)
        sc = Scene("ga", 740, 450, yaw=0.6, pitch=0.42, persp=0.3, fog=0.4, outline=True, atom_k=1.3)
        render_cell(sc, cell, "All atoms")
        sc.add_box((0, 0, 0), (Lz, 0, 0), (0, Lz, 0), (0, 0, Lz), MUTED, 1.2, "5 4", 0.8)
        n_shown, n_h = count(cell)
        overlay = style_bar("All atoms", f"{Lz:.0f} Å region around chain 12 · {n_shown:,} of 12 840 atoms shown · {n_h} H".replace(",", " "))
    else:
        overlay = style_bar("Backbone", "Backbone tube · H hidden · 12 840 atoms in model") + zoom_inset("gz", 230, 170, 3, 7)
    view = sc.svg()

    comp = table(["Component", "Chains", "DP", "Dispersity", "Atoms"],
                 [[f'{dot(ACC)} PS atactic', "20", "40", "Schulz–Zimm · Đ 1.10", "12 840"],
                  [f'{dot(SEL)} Toluene (solvent)', "—", "—", "—", "fill to ρ"]],
                 ["34%", "14%", "11%", "26%", "15%"], mono_cols=(1, 2, 4), align_right=(1, 2, 4))
    method = []
    radio_dot = f'<span style="width: 8px; height: 8px; border-radius: 50%; background: {ACC}"></span>'
    for name, sub_, on in (("Theodorou–Suter + RIS", "Rosenbluth-weighted, soft-sphere excluded volume · Macromolecules 1985", False),
                           ("Configurational-bias, real force field", "Siepmann &amp; Frenkel 1992 · trial energies from CAPS Field (GAFF2)", True)):
        method.append(f'<button role="radio" aria-checked="{"true" if on else "false"}" style="display: flex; gap: 10px; align-items: flex-start; padding: 12px; text-align: left; background: {BG3 if on else BG0}; border: 1px solid {ACC if on else LINE}; border-radius: 8px; cursor: pointer">'
                      f'<span style="width: 16px; height: 16px; margin-top: 1px; border-radius: 50%; border: 2px solid {ACC if on else DIM}; display: flex; align-items: center; justify-content: center; flex-shrink: 0">{radio_dot if on else ""}</span>'
                      f'<span style="display: flex; flex-direction: column; gap: 3px"><span style="font-size: 13px; font-weight: 600; color: {TEXT}">{name}</span><span style="font-size: 12px; color: {MUTED}; line-height: 1.4">{sub_}</span></span></button>')
    left = (f'<div style="width: 560px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 0 18px 22px; overflow: hidden">'
            + card("Components", col(comp, row(btn("Add polymer", ic="plus", small=True), btn("Add solvent / gas", ic="flask", small=True), gap=6), gap=10), chip("Polydisperse"))
            + card("Cell", col(row(field("Target density", "1.04", "g/cm³"), field("Temperature", "298", "K"), select("Shape", "Cubic, periodic"), gap=10),
                                row(select("Orientation", "Isotropic"), field("Initial density", "0.50", "g/cm³", note="grow sparse, then compress"), gap=10), gap=10))
            + card("Growth method", col(*method, row(field("Trial directions k", "16"), field("Look-ahead", "2", "bonds"), select("Force field", "GAFF2 (typed)"), gap=10), gap=8))
            + '</div>')
    live = [("Chains grown", "11 / 20"), ("Current chain", "26 / 40 units"), ("Acceptance", "0.63"), ("Rosenbluth log W", "−412.7"),
            ("Density now", "0.50 g/cm³"), ("Close contacts &lt; 1.6 Å", "0")]
    lv = "".join(f'<div style="display: flex; flex-direction: column; gap: 3px; padding: 10px 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px"><span style="font-size: 11.5px; color: {MUTED}">{k}</span><span style="font-family: {MONO}; font-size: 15px">{v}</span></div>' for k, v in live)
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
             + card("Live cell", f'<div style="position: relative; margin: -16px; background: {BG0}">{view}'
                    f'<div style="position: absolute; left: 12px; top: 12px; display: flex; gap: 6px">{chip("Colour: molecule", TEXT, "#16191Ccc")}{chip(dot(ACC, 7) + " growing chain 12", ACC, "#16191Ccc")}</div>{overlay}</div>',
                    row(chip("8 domains · 16 threads", MUTED), tbtn("pause", "Pause"), tbtn("stop", "Cancel"), gap=4), pad=16)
             + f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 8px">{lv}</div>'
             + col(row(f'<span style="font-size: 12px; color: {MUTED}">Overall</span>', '<div style="flex-grow: 1"></div>', f'<span style="font-family: {MONO}; font-size: 12px">57 %</span>'), progress(0.57), gap=6)
             + '</div>')
    head = pagehead("Amorphous cell · CAPS Grow", "Grow chains into a periodic cell, then hand off to Pack, Relax and Equilibrate",
                    row(btn("Save recipe", ic="save"), btn("Copy as Python", ic="copy"), btn("Queue growth", True, "play", href="Jobs.dc.html"), gap=8))
    content = head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>'
    return shell("Grow", "CAPS Grow — amorphous cell", f'{icon("chevr", 12, DIM)}<span>Grow</span>', content,
                 ("<span>Job grow-3 · running</span><span>Seed 20260923 · Philox</span>", "<span>Provenance manifest written on finish</span>"))


# ------------------------------------------------------------------ Force field
def forcefield():
    atoms, backbone, stereo = polystyrene(2, "atactic", seed=2, curve=False)
    sc = Scene("f", 496, 330, yaw=-0.2, pitch=0.55, persp=0.15, fog=0.4)
    tcol = {"c3": "#8E959C", "ca": "#D6A45E", "hc": "#E9ECEF", "ha": "#6CC4D8"}
    typed = []
    for a in atoms:
        if a["e"] == "C":
            t = "ca" if a.get("ar") else "c3"
        else:
            # H on aromatic if near an aromatic C
            t = "ha" if any(b.get("ar") and dist(a["p"], b["p"]) < 1.2 for b in atoms) else "hc"
        typed.append({**a, "c": tcol[t], "t": t})
    sc.add_atoms(typed)
    sc._prep()
    ov = ""
    shown = set()
    for i, a in enumerate(typed):
        if a["t"] not in shown or i in (3, 9):
            shown.add(a["t"])
            x, y, _, _ = sc.screen(i)
            ov += label_pill(x + 8, y - 14, a["t"], tcol[a["t"]], fs=10.5)
    # missing param highlight on torsion c3-c3-ca-ca
    view = sc.svg(overlay=ov)

    rows = [
        ["C1", "C", "c3", "−0.0914", "[CX4;H2]", "gaff2.dat"],
        ["C2", "C", "c3", "−0.0960", "[CX4;H1]", "gaff2.dat"],
        ["C3", "C", "ca", "−0.0152", "[c;R1](c)(c)[CX4]", "gaff2.dat"],
        ["C4", "C", "ca", "−0.1180", "[c;H1]", "gaff2.dat"],
        ["H9", "H", "hc", "0.0612", "[H][CX4]", "gaff2.dat"],
        ["H14", "H", "ha", "0.1296", "[H]c", "gaff2.dat"],
        ["C17", "C", f'<span style="color: {ERR}">?</span>', "—", f'<span style="color: {ERR}">no rule matched</span>', "—"],
    ]
    tbl = table(["Atom", "El.", "Type", "q (e)", "Rule (SMARTS)", "Source"], rows, ["9%", "7%", "10%", "13%", "41%", "20%"], mono_cols=(0, 2, 3, 4, 5), align_right=(3,), hl={6})
    miss = [("Dihedral", "c3–c3–ca–ca", "Not in GAFF2 parameter set"), ("Atom type", "C17 (sp² C=N⁺)", "No typing rule matched")]
    mr = "".join(f'<div style="display: flex; gap: 10px; align-items: flex-start; padding: 10px 12px; background: #2A1414; border: 1px solid #5A2A2A; border-radius: 7px">{icon("xcircle", 16, ERR)}'
                 f'<div style="display: flex; flex-direction: column; gap: 3px; flex-grow: 1"><span style="font-size: 12.5px"><b style="font-weight: 600">{k}</b> <span style="font-family: {MONO}">{v}</span></span><span style="font-size: 12px; color: {MUTED}">{d}</span></div></div>' for k, v, d in miss)
    fix = col(mr,
              f'<div style="font-size: 12px; color: {MUTED}; line-height: 1.5">CAPS never guesses parameters. Choose a source; anything you enter by hand is flagged <span style="font-family: {MONO}; color: {WARN}">estimated</span> in the output and the provenance manifest.</div>',
              row(btn("Import .itp / .frcmod", ic="download", small=True), btn("Enter manually", ic="sliders", small=True), gap=6), gap=10)
    legend = "".join(f'<span style="display: flex; align-items: center; gap: 6px; font-family: {MONO}; font-size: 11.5px; color: {MUTED}">{dot(c)}{t}</span>' for t, c in tcol.items())
    head = pagehead("Force field · typing report", "Explain every type, override with provenance, and block runs with missing parameters",
                    row(select("Force field", "GAFF2 2.11", w=180), btn("Retype all", ic="rotate"), btn("Export .itp / .lt / data", ic="download"), gap=8))
    stats = row(chip(f'{dot(OK)} 1 204 typed', TEXT), chip(f'{dot(ERR)} 1 untyped', ERR, "#2A1414"), chip(f'{dot(ERR)} 1 missing torsion', ERR, "#2A1414"), chip(f'{dot(WARN)} 0 estimated', MUTED), chip("net charge 0.000 e", MUTED, BG2, True), gap=8)
    left = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 0 18px 22px">'
            + stats
            + card("Atoms", col(row(f'<div style="flex-grow: 1; display: flex; align-items: center; gap: 8px; height: 30px; padding: 0 8px; background: {BG0}; border: 1px solid {LINE}; border-radius: 5px; color: {DIM}; font-size: 12px">{icon("filter", 14, DIM)}<span>Filter: all · type, element, rule, charge</span></div>', seg(["All", "Untyped", "Overridden"], "All"), gap=8), tbl, gap=10), extra="flex-grow: 1")
            + card("Why C3 is ca", f'<div style="font-size: 12.5px; line-height: 1.6; color: {MUTED}">Rule <span style="font-family: {MONO}; color: {TEXT}">gaff2.ca</span> (priority 40) matched <span style="font-family: {MONO}; color: {TEXT}">[c;R1](c)(c)[CX4]</span>: aromatic carbon in a 6-ring with an sp³ substituent. Rules <span style="font-family: {MONO}">cp</span>, <span style="font-family: {MONO}">cq</span> (biphenyl bridge) were checked and rejected. Charges: AM1-BCC set shipped with the typed fragment library, neutral per molecule.</div>')
            + '</div>')
    right = (f'<div style="width: 540px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
             + card("Types in 3D", f'<div style="margin: -16px; background: {BG0}">{view}</div>', f'<div style="display: flex; gap: 12px">{legend}</div>')
             + card("Missing before run", fix, chip("Blocks Dynamics", ERR, "#2A1414"))
             + '</div>')
    content = head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>'
    return shell("Field", "CAPS Field — typing report", f'{icon("chevr", 12, DIM)}<span>Field</span>', content,
                 ("<span>Typer: SMARTS rule file gaff2.rules v3</span>", "<span>Parameter DB: caps-ffdb 1.0 · GAFF2 (Wang et al.)</span>"))


# ------------------------------------------------------------------ Jobs / Equilibrate monitor
def jobs():
    jl = [("equil-7", "Equilibrate · 21-step", "running", ACC, "Local · RTX 4070", 0.66),
          ("grow-3", "Grow · PS cell", "running", ACC, "Local · 16 threads", 0.57),
          ("dyn-12", "Dynamics · NPT 10 ns", "queued", MUTED, "hpc-login2 · gRPC", 0.0),
          ("pack-5", "Pack · PS in water", "done", OK, "Local", 1.0),
          ("react-2", "React · epoxy 85 %", "failed", ERR, "hpc-login2", 0.41)]
    rows = []
    for jid, name, st, c, host, fr in jl:
        on = jid == "equil-7"
        rows.append(f'<button style="display: flex; flex-direction: column; gap: 6px; padding: 10px 12px; text-align: left; background: {BG3 if on else "transparent"}; border: 1px solid {LINE if on else "transparent"}; border-radius: 7px; cursor: pointer">'
                    f'<span style="display: flex; align-items: center; gap: 8px; width: 100%">{dot(c)}<span style="font-size: 12.5px; font-weight: 500; color: {TEXT}; flex-grow: 1">{name}</span><span style="font-family: {MONO}; font-size: 11px; color: {c}">{st}</span></span>'
                    f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{jid} · {host}</span>'
                    + (progress(fr, c, 4) if 0 < fr < 1 else "") + '</button>')
    joblist = (f'<aside aria-label="Jobs" style="width: 290px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-right: 1px solid {LINE}">'
               + panel_head("Jobs", btn("New", ic="plus", small=True))
               + f'<div style="padding: 8px; display: flex; flex-direction: column; gap: 4px">{"".join(rows)}</div></aside>')
    # 21-step Larsen protocol: 7 cycles of 3 stages (NVT high T, NVT 300 K, NPT)
    steps = []
    for i in range(21):
        kind = ["NVT 600 K", "NVT 300 K", "NPT"][i % 3]
        c = OK if i < 13 else (ACC if i == 13 else BG3)
        steps.append(f'<span title="Step {i + 1}: {kind}" style="height: {34 if kind == "NPT" else 22}px; flex-grow: 1; border-radius: 3px; background: {c}; align-self: flex-end"></span>')
    stepbar = f'<div style="display: flex; gap: 3px; height: 36px; align-items: flex-end">{"".join(steps)}</div>'
    # illustrative density series across the compression/decompression
    pts = []
    t = 0.0
    rho = 0.50
    targets = [0.6, 0.62, 0.7, 0.74, 0.85, 0.9, 1.05, 1.1, 1.25, 1.2, 1.16, 1.12, 1.09, 1.06]
    for k, tg in enumerate(targets):
        for s in range(12):
            rho += (tg - rho) * 0.22
            pts.append((t, rho + 0.004 * math.sin(t * 7.3)))
            t += 0.05
    Tpts = []
    for k in range(14):
        T = [600, 300, 300][k % 3]
        for s in range(12):
            tt = k * 0.6 + s * 0.05
            Tpts.append((tt, T + 12 * math.sin(tt * 11) * (0.6 if T == 300 else 1)))
    dplot = plot(560, 210, [(pts, ACC, 2, None)], (0, 9), (0.4, 1.3), [0, 3, 6, 9], [0.4, 0.7, 1.0, 1.3], "time (ns)", "density (g/cm³)", pad=(44, 12, 18, 30))
    tplot = plot(400, 210, [(Tpts, SEL, 1.6, None)], (0, 9), (200, 700), [0, 3, 6, 9], [200, 450, 700], "time (ns)", "temperature (K)", pad=(44, 12, 18, 30))
    log = [("12:41:07", "step 13 NPT P=50000 bar done · ρ 1.12 g/cm³"),
           ("12:41:07", "checkpoint written equil-7/ckpt_013.caps"),
           ("12:41:08", "step 14 NVT 600 K start · Bussi thermostat τ=0.1 ps"),
           ("12:43:52", "E drift 3.1e−6 kT/ns/atom (NVE probe) · OK")]
    lg = "".join(f'<div style="display: flex; gap: 12px"><span style="color: {DIM}">{a}</span><span>{esc(b)}</span></div>' for a, b in log)
    manifest = [("caps", "0.1.0-dev+g8c1e2"), ("inputs", "sha256 9f2c…a41e"), ("ff", "GAFF2 · caps-ffdb 1.0"), ("electrostatics", "SPME · rc 12 Å · 1e−5"),
                ("vdW", "cut 12 Å + tail corr."), ("precision", "mixed (deterministic)"), ("seed", "20260923 · Philox4x32")]
    mf = "".join(kv(k, v) for k, v in manifest)
    main = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
            + pagehead("Equilibrate · Larsen 21-step protocol", "equil-7 · PS atactic 20 × DP40 · 12 840 atoms · started 11:02",
                       row(chip("step 14 of 21", ACC, "#3A2C14"), btn("Open cell in Studio", ic="cube", href="Main.dc.html"), btn("Pause", ic="pause"), btn("Cancel", ic="stop"), gap=8))
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px; min-height: 0">'
            + card("Protocol", col(stepbar, row(chip(f'{dot(OK)} done 13'), chip(f'{dot(ACC)} running 1'), chip(f'{dot(BG3)} remaining 7'), '<div style="flex-grow: 1"></div>',
                                                 f'<span style="font-size: 12px; color: {MUTED}">Cycles of NVT 600 K → NVT 300 K → NPT, P<sub>max</sub> 50 000 bar · Larsen et al., <i>Macromolecules</i> 44, 6944 (2011)</span>', gap=8), gap=10))
            + f'<div style="display: flex; gap: 14px">' + card("Density", dplot, chip("live"), 12) + card("Temperature", tplot, chip("live"), 12) + '</div>'
            + f'<div style="display: flex; gap: 14px; flex-grow: 1; min-height: 0">'
            + card("Log", f'<div style="font-family: {MONO}; font-size: 11.5px; line-height: 1.8; color: {TEXT}">{lg}</div>', chip("structured · JSON"), extra="flex-grow: 1")
            + card("Provenance manifest", col(mf, gap=6), btn("View JSON", small=True), extra="width: 380px; flex-shrink: 0")
            + '</div></div></div>')
    content = f'<div style="flex-grow: 1; display: flex; min-height: 0">{joblist}{main}</div>'
    return shell("Jobs", "CAPS Jobs — equilibration monitor", f'{icon("chevr", 12, DIM)}<span>Jobs</span>', content,
                 ("<span>2 running · 1 queued</span><span>Checkpoint every 50 ps</span>", "<span>Remote: hpc-login2 connected</span>"))


# ------------------------------------------------------------------ Analyze
def analyze():
    # RDF-like curve (illustrative shape)
    g = []
    for i in range(1, 300):
        r = i * 0.05
        v = 0 if r < 2.2 else 1 + 1.6 * math.exp(-((r - 2.55) ** 2) / 0.02) + 0.55 * math.exp(-((r - 3.9) ** 2) / 0.15) * 1 + 0.25 * math.sin(r * 4.2) * math.exp(-(r - 4) / 2.2) * (r > 4)
        if r < 2.6 and r >= 2.2:
            v = v * (1 - math.exp(-(r - 2.2) * 8))
        g.append((r, v))
    rdf = plot(470, 230, [(g, ACC, 1.8, None)], (0, 15), (0, 3), [0, 5, 10, 15], [0, 1, 2, 3], "r (Å)", "g(r) · C–C intermolecular", pad=(40, 12, 18, 30))
    # Tg cooling curve with bilinear fit (illustrative)
    Ts = list(range(200, 520, 20))
    sp = []
    for T in Ts:
        v = 0.930 + 0.00020 * (T - 200) if T < 370 else 0.930 + 0.00020 * 170 + 0.00058 * (T - 370)
        v += 0.0015 * math.sin(T * 0.7)
        sp.append((T, v))
    fit1 = [(200, 0.930), (400, 0.930 + 0.0002 * 200)]
    fit2 = [(320, 0.964 + 0.00058 * (-50)), (500, 0.964 + 0.00058 * 130)]
    tg = plot(470, 230, [(fit1, MUTED, 1.2, "4 3"), (fit2, MUTED, 1.2, "4 3")], (200, 500), (0.92, 1.05), [200, 300, 400, 500], [0.92, 0.96, 1.0, 1.04],
              "T (K)", "specific volume (cm³/g)", pad=(44, 12, 18, 30), markers=[(a, b, SEL) for a, b in sp])
    props = [("Density · 298 K", "[result]", "exp. 1.04–1.065 g/cm³", "Mark, Physical Properties of Polymers Handbook"),
             ("T<sub>g</sub> · bilinear fit", "[result]", "exp. ≈ 373 K", "Brandrup et al., Polymer Handbook"),
             ("C<sub>∞</sub>", "[result]", "lit. ≈ 9.5–10", "Flory 1969; Mattice &amp; Suter 1994"),
             ("δ Hildebrand", "[result]", "exp. 17.4–19.0 MPa<sup>½</sup>", "Polymer Handbook")]
    pc = "".join(f'<div style="display: flex; flex-direction: column; gap: 6px; padding: 14px; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px">'
                 f'<span style="font-size: 12px; color: {MUTED}">{a}</span><span style="font-family: {MONO}; font-size: 22px; color: {DIM}">{b}</span>'
                 f'<span style="font-size: 12px; color: {TEXT}">{c}</span><span style="font-size: 11px; color: {DIM}">{d}</span></div>' for a, b, c, d in props)
    calcs = [("Structure", ["Density", "RDF", "S(q)", "X-ray / neutron"]), ("Chains", ["R<sub>g</sub>", "R<sub>ee</sub>", "C<sub>n</sub>, C<sub>∞</sub>", "Persistence", "Entanglements"]),
             ("Thermo", ["CED", "δ", "T<sub>g</sub>"]), ("Mechanics", ["C<sub>ij</sub> strain", "C<sub>ij</sub> fluct.", "Stress–strain"]),
             ("Dynamics", ["MSD", "D", "Relaxation"]), ("Free volume", ["Probe insertion", "Pore size"])]
    cl = []
    for grp, items in calcs:
        its = "".join(f'<button style="height: 26px; padding: 0 9px; border-radius: 5px; background: {BG3 if it in ("RDF", "T<sub>g</sub>") else BG0}; border: 1px solid {ACC if it in ("RDF", "T<sub>g</sub>") else LINE}; font-size: 12px; color: {TEXT}; cursor: pointer">{it}</button>' for it in items)
        cl.append(f'<div style="display: flex; flex-direction: column; gap: 6px"><span style="font-size: 11px; font-weight: 600; letter-spacing: 0.06em; text-transform: uppercase; color: {DIM}">{grp}</span><div style="display: flex; flex-wrap: wrap; gap: 5px">{its}</div></div>')
    side = (f'<aside aria-label="Calculations" style="width: 300px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-right: 1px solid {LINE}">'
            + panel_head("Calculations")
            + f'<div style="padding: 14px; display: flex; flex-direction: column; gap: 16px">{"".join(cl)}</div>'
            + section("Source", col(select("Trajectory", "equil-7 · last 5 ns"), select("Frames", "every 10 ps · 500 frames"), row(select("Groups", "all PS"), gap=8), gap=10))
            + '</aside>')
    main = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
            + pagehead("Analyze · PS atactic cell", "Properties compared against experiment, with sources, ready for the paper tables",
                       row(btn("Export CSV / LaTeX", ic="download"), btn("Run selected", True, "play"), gap=8))
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
            + f'<div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px">{pc}</div>'
            + '<div style="display: flex; gap: 14px">'
            + card("Radial distribution", rdf, chip("illustrative"), 12, "flex-grow: 1")
            + card("T<sub>g</sub> from cooling", tg, chip("illustrative"), 12, "flex-grow: 1") + '</div>'
            + card("Cooling protocol", row(field("From", "500", "K"), field("To", "200", "K"), field("Step", "20", "K"), field("NPT per step", "2", "ns"), select("Fit", "Bilinear, free hinge"), gap=10))
            + '</div></div>')
    content = f'<div style="flex-grow: 1; display: flex; min-height: 0">{side}{main}</div>'
    return shell("Analyze", "CAPS Analyze — properties", f'{icon("chevr", 12, DIM)}<span>Analyze</span>', content,
                 ("<span>Trajectory equil-7 · 500 frames</span>", "<span>Results carry citations in provenance</span>"))


# ------------------------------------------------------------------ Bench
def bench():
    tables_ = [("1", "Energy &amp; force parity", "vs LAMMPS, GROMACS", "pass", 14, 14), ("2", "Packing vs packmol", "≥ 10× on ≥ 100k atoms", "running", 3, 8),
               ("3", "Amorphous cell builders", "MS, Polymatic, PySIMM, RadonPy", "not run", 0, 10), ("4", "MD throughput", "ns/day", "not run", 0, 6),
               ("5", "NVE energy conservation", "drift kT/ns/atom", "pass", 6, 6), ("6", "Properties vs experiment", "PE PP PS PMMA PET PC PA6 PDMS", "not run", 0, 8),
               ("7", "Chain statistics", "C∞ vs literature", "fail", 5, 8), ("8", "Parallel scaling", "cores / GPUs", "not run", 0, 5),
               ("9", "Cross-platform reproducibility", "bit-identical", "pass", 4, 4), ("10", "Builder feature coverage", "cited vendor docs", "draft", 0, 0),
               ("11", "Rendering &amp; editing performance", "fps, load, latency", "not run", 0, 4), ("12", "Builder accuracy", "CIP set, geometry, space groups", "not run", 0, 3)]
    stc = {"pass": OK, "running": ACC, "not run": DIM, "fail": ERR, "draft": MUTED}
    rows = []
    for n, name, sub_, st, a, b in tables_:
        prog = f'{a} / {b}' if b else "—"
        rows.append([f'<span style="font-family: {MONO}; color: {DIM}">T{n}</span>', f'<span style="color: {TEXT}">{name}</span>', f'<span style="color: {MUTED}">{sub_}</span>',
                     prog, f'<span style="display: inline-flex; align-items: center; gap: 6px; color: {stc[st]}">{dot(stc[st])}{st}</span>'])
    t = table(["", "Table", "Scope", "Rows", "Status"], rows, ["5%", "33%", "40%", "8%", "14%"], mono_cols=(3,), align_right=(3,), hl={1}, rowh=34, fs=11.5)
    t2rows = [["Water box", "33 333", "100 000", "1.00", "packmol", "1", "[measure]", "[measure]", "[5 runs]", "1.0×"],
              ["Water box", "33 333", "100 000", "1.00", "CAPS", "16", "[measure]", "[measure]", "[5 runs]", "—"],
              ["PS melt fragments", "2 000", "128 000", "1.04", "packmol", "1", "[measure]", "[measure]", "[5 runs]", "1.0×"],
              ["PS melt fragments", "2 000", "128 000", "1.04", "CAPS", "GPU", "[measure]", "[measure]", "[5 runs]", "—"]]
    t2 = table(["System", "Mol.", "Atoms", "ρ", "Tool", "Thr.", "Wall (s)", "d<sub>min</sub> (Å)", "Success", "Speed-up"], t2rows, ["17%", "8%", "9%", "6%", "10%", "6%", "11%", "11%", "11%", "11%"],
               mono_cols=(1, 2, 3, 5, 6, 7, 8, 9), align_right=(1, 2, 3, 6, 7, 8, 9), fs=11, hl={1, 3})
    env = [("CPU", "AMD Ryzen 9 7950X · 16 C"), ("GPU", "NVIDIA RTX 4070 · 12 GB"), ("Compiler", "clang 18.1 · -O3"), ("packmol", "20.14.4 (pinned)"),
           ("LAMMPS", "2Aug2023 (pinned)"), ("GROMACS", "2024.3 (pinned)"), ("Repeats", "5 · mean ± sd")]
    head = pagehead("Validation &amp; benchmarks", "One command regenerates Tables 1–12 as Markdown, CSV and LaTeX",
                    row(f'<code style="font-family: {MONO}; font-size: 12px; padding: 7px 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px; color: {TEXT}">caps bench run --all --repeats 5</code>',
                        btn("Export for paper", ic="download"), btn("Run suite", True, "play"), gap=8))
    left = (f'<div style="width: 620px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 0 18px 22px">'
            + card("Benchmark tables", t, row(chip(f'{dot(OK)} 3 pass'), chip(f'{dot(ERR)} 1 fail', ERR, "#2A1414"), gap=6), 12, "flex-grow: 1")
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
             + card("Table 2 · Packing performance against packmol", col(t2,
                    f'<div style="display: flex; gap: 10px; align-items: flex-start; padding: 10px 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px; font-size: 12px; color: {MUTED}; line-height: 1.5">{icon("alert", 16, WARN)}'
                    f'<span>No speed-up is claimed until all five repeats finish on pinned versions. Rows where CAPS is slower or less accurate are published as-is.</span></div>', gap=12), chip("running · 3 / 8", ACC, "#3A2C14"), 12)
             + '<div style="display: flex; gap: 14px; flex-grow: 1; min-height: 0">'
             + card("Table 7 · failing row", col(
                 kv("System", "PP isotactic · 413 K", False),
                 kv("Tolerance", "±10 % of literature C∞", False),
                 kv("Status", f'<span style="color: {ERR}">outside tolerance</span>', False),
                 f'<div style="font-size: 12px; color: {MUTED}; line-height: 1.5">Open the run to compare RIS weights against Suter &amp; Flory statistical weights before the phase-3 exit review.</div>',
                 btn("Open run", ic="chart", small=True), gap=8), extra="flex-grow: 1")
             + card("Pinned environment", col(*[kv(k, v, False) for k, v in env], gap=6), chip("hash 4b1e…"), extra="width: 330px; flex-shrink: 0")
             + '</div></div>')
    content = head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>'
    return shell("Bench", "CAPS Bench — validation suite", f'{icon("chevr", 12, DIM)}<span>Bench</span>', content,
                 ("<span>Suite v0.3 · 76 rows</span><span>Nightly on Linux, macOS arm64/x86, Windows</span>", "<span>Raw data archived with each run</span>"))
