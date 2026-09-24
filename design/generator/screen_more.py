import math, random
from lib import *
from mols import *
from screen_app import shell, pagehead, card, progress


def formula(html):
    return (f'<div style="padding: 12px 14px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px; font-family: {MONO}; font-size: 13px; line-height: 1.7; color: {TEXT}">{html}</div>')


def radio_list(items):
    out = []
    for name, sub_, on in items:
        dotc = f'<span style="width: 8px; height: 8px; border-radius: 50%; background: {ACC}"></span>' if on else ""
        out.append(f'<button role="radio" aria-checked="{"true" if on else "false"}" style="display: flex; gap: 10px; align-items: flex-start; padding: 10px 12px; text-align: left; background: {BG3 if on else BG0}; border: 1px solid {ACC if on else LINE}; border-radius: 7px; cursor: pointer">'
                   f'<span style="width: 16px; height: 16px; margin-top: 1px; border-radius: 50%; border: 2px solid {ACC if on else DIM}; display: flex; align-items: center; justify-content: center; flex-shrink: 0">{dotc}</span>'
                   f'<span style="display: flex; flex-direction: column; gap: 2px"><span style="font-size: 12.5px; font-weight: 600; color: {TEXT}">{name}</span><span style="font-size: 11.5px; color: {MUTED}; line-height: 1.4">{sub_}</span></span></button>')
    return col(*out, gap=6)


def stat_grid(items, cols=3):
    cells = "".join(f'<div style="display: flex; flex-direction: column; gap: 3px; padding: 10px 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px"><span style="font-size: 11.5px; color: {MUTED}">{k}</span><span style="font-family: {MONO}; font-size: 15px; color: {c}">{v}</span></div>' for k, v, c in items)
    return f'<div style="display: grid; grid-template-columns: repeat({cols}, minmax(0, 1fr)); gap: 8px">{cells}</div>'


# ------------------------------------------------------------------ Pack
def pack():
    Lb = 40.0
    sc = Scene("k", 700, 390, yaw=0.62, pitch=0.42, persp=0.3, fog=0.6)
    rnd = random.Random(8)
    placed = []
    # slab region of fixed silica-like beads at bottom, molecules above (region constraint demo)
    mols = 0
    while mols < 90:
        c = (rnd.uniform(2, Lb - 2), rnd.uniform(2, Lb - 2), rnd.uniform(9, Lb - 2))
        if all(dist(c, q) > 5.2 for q in placed):
            placed.append(c)
            d = norm((rnd.uniform(-1, 1), rnd.uniform(-1, 1), rnd.uniform(-1, 1)))
            a3 = [{"e": "C", "p": add(c, mul(d, k * 1.54 - 1.54)), "c": "#B9BEC4" if mols % 5 else SEL} for k in range(3)]
            sc.add_atoms(a3, [(0, 1), (1, 2)])
            mols += 1
    fixed = []
    for i in range(10):
        for j in range(10):
            fixed.append({"e": "Si", "p": (2 + i * 4.0, 2 + j * 4.0, 3.0 + (i + j) % 2 * 1.5)})
    sc.add_atoms(fixed, [])
    sc.add_box((0, 0, 0), (Lb, 0, 0), (0, Lb, 0), (0, 0, Lb), ACC, 1.6, None, 0.9)
    sc.add_box((0, 0, 0), (Lb, 0, 0), (0, Lb, 0), (0, 0, 7.0), "#D6A45E", 1.0, "3 3", 0.8)
    view = sc.svg()
    # penalty convergence (illustrative)
    it = [(i, 6.0 - 0.055 * i - 1.4 * math.exp(-i / 6) + (0.25 * math.sin(i) if i < 40 else 0)) for i in range(0, 120)]
    it = [(x, max(-3.2, y)) for x, y in it]
    pl = plot(470, 200, [(it, ACC, 1.8, None)], (0, 120), (-4, 6), [0, 40, 80, 120], [-4, -2, 0, 2, 4, 6], "iteration", "log₁₀ penalty f", pad=(40, 12, 18, 30),
              bands=[(0, 12, SEL, 0.08), (12, 90, ACC, 0.05), (90, 120, OK, 0.07)])
    regions = [("Box", "0 0 0 → 40 40 40 Å", "all molecules", ACC), ("Slab below z = 7 Å", "fixed", "silica surface · 100 sites", "#D6A45E"),
               ("Outside slab", "z &gt; 8 Å", "propane × 90", SEL)]
    rg = "".join(f'<div style="display: flex; align-items: center; gap: 10px; padding: 9px 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px">'
                 f'<span style="width: 4px; height: 26px; border-radius: 2px; background: {c}"></span><div style="display: flex; flex-direction: column; gap: 2px; flex-grow: 1"><span style="font-size: 12.5px">{a}</span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{b}</span></div>'
                 f'<span style="font-size: 11.5px; color: {MUTED}">{t}</span></div>' for a, b, t, c in regions)
    left = (f'<div style="width: 520px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 0 18px 22px">'
            + card("Molecules &amp; regions", col(rg, row(btn("Add molecule", ic="plus", small=True), btn("Add region", ic="cube", small=True), chip("box · slab · sphere · cylinder"), gap=6), gap=8))
            + card("Objective", col(formula('f = Σ<sub>i&lt;j</sub> max(0, d<sub>tol</sub>² − ‖x<sub>i</sub> − x<sub>j</sub>‖²)²<br>&#160;&#160;+ Σ<sub>i</sub> w · g<sub>region</sub>(x<sub>i</sub>)'),
                                    cite("Overlap penalty after Martínez et al., <i>J. Comput. Chem.</i> 30, 2157 (2009). Rigid bodies: centre of mass + unit quaternion, analytic gradients."),
                                    row(field("Tolerance d<sub>tol</sub>", "2.0", "Å"), select("Optimiser", "L-BFGS (m = 10)"), select("Neighbour list", "Linked cell · Morton"), gap=8), gap=10))
            + card("Stages", col(*[f'<div style="display: flex; gap: 10px; font-size: 12.5px">{dot(c, 9)}<span style="flex-grow: 1">{a}</span><span style="font-family: {MONO}; font-size: 11.5px; color: {c}">{s}</span></div>' for a, s, c in (("1 · Random sequential insertion", "done", OK), ("2 · Overlap minimisation (L-BFGS)", "running", ACC), ("3 · Soft-core compression, verify d<sub>min</sub>", "pending", DIM))], gap=10), chip("stage 2 of 3", ACC, "#3A2C14"))
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
             + card("Packed cell", f'<div style="position: relative; margin: -16px; background: {BG0}">{view}'
                    f'<div style="position: absolute; left: 12px; top: 12px; display: flex; gap: 6px">{chip("stage 2 · overlap minimisation", ACC, "#16191Ccc")}{chip("GPU kernel", MUTED, "#16191Ccc")}</div></div>',
                    row(tbtn("pause", "Pause"), tbtn("stop", "Cancel"), gap=4))
             + '<div style="display: flex; gap: 14px">'
             + card("Convergence", pl, chip("illustrative"), 12, "flex-grow: 1")
             + card("Guarantee", col(stat_grid([("d<sub>min</sub> now", "1.74 Å", WARN), ("target", "≥ 2.00 Å", TEXT), ("overlaps", "37", WARN), ("iteration", "58", TEXT)], 2),
                                      f'<div style="font-size: 11.5px; color: {MUTED}; line-height: 1.45">Run fails loudly if d<sub>min</sub> is not met; it never writes a cell below tolerance.</div>', gap=10), extra="width: 250px; flex-shrink: 0")
             + '</div></div>')
    head = pagehead("Pack · CAPS Pack", "Native packing with region constraints — no packmol subprocess",
                    row(btn("Copy as Python", ic="copy"), btn("Queue packing", True, "play"), gap=8))
    return shell("Pack", "CAPS Pack — packing", f'{icon("chevr", 12, DIM)}<span>Pack</span>', head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>',
                 ("<span>90 propane + fixed silica · 40 Å box</span>", "<span>Deterministic mode · Philox seed 7</span>"), ["propane_on_silica.caps"])


# ------------------------------------------------------------------ Relax
def relax():
    e = [(i, 3.2 * math.exp(-i / 18) + 0.4 * math.exp(-i / 90) - 0.02 * i / 200) for i in range(0, 400, 2)]
    f = [(i, 2.2 - 0.0105 * i + 0.12 * math.sin(i / 7) * math.exp(-i / 150)) for i in range(0, 400, 2)]
    pe = plot(430, 220, [(e, ACC, 1.8, None)], (0, 400), (0, 3.6), [0, 100, 200, 300, 400], [0, 1, 2, 3], "iteration", "E − E_final (10³ kcal/mol)", pad=(44, 12, 18, 30),
              bands=[(0, 60, SEL, 0.07)])
    pf = plot(320, 220, [(f, SEL, 1.8, None), ([(0, -1.7), (400, -1.7)], OK, 1.2, "5 4")], (0, 400), (-2.5, 2.5), [0, 100, 200, 300, 400], [-2, -1, 0, 1, 2], "iteration", "log₁₀ |F|max (kcal/mol/Å)", pad=(44, 12, 18, 30))
    methods = radio_list([("Steepest descent", "Robust first pass for bad contacts", False),
                          ("Polak–Ribière conjugate gradient", "Low memory, smooth landscapes", False),
                          ("L-BFGS", "Liu &amp; Nocedal, <i>Math. Program.</i> 45, 503 (1989) · m = 10", True),
                          ("FIRE", "Bitzek et al., <i>Phys. Rev. Lett.</i> 97, 170201 (2006)", False)])
    cons = [("Fixed atoms", "chain-end carbons · 4 atoms"), ("Distance restraint", "C12–C31 = 3.80 Å · k 50"), ("Dihedral restraint", "φ C5–C6–C7–C8 = 180° · k 10")]
    cr = "".join(f'<div style="display: flex; align-items: center; gap: 10px; height: 34px; padding: 0 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px">{icon("pin", 14, ACC)}<span style="font-size: 12.5px; flex-grow: 1">{a}</span><span style="font-family: {MONO}; font-size: 11.5px; color: {MUTED}">{b}</span></div>' for a, b in cons)
    left = (f'<div style="width: 480px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 0 18px 22px">'
            + card("Method", methods)
            + card("Push-off &amp; cell", col(toggle("Soft-core push-off before minimising", True),
                                              row(field("Force cap", "50", "kcal/mol/Å"), field("λ ramp", "0 → 1", note="over 20 ps"), gap=10),
                                              cite("Auhl et al., <i>J. Chem. Phys.</i> 119, 12718 (2003)"),
                                              toggle("Relax box shape and volume", True), row(select("Coupling", "Anisotropic"), field("Pressure", "1", "atm"), gap=10), gap=10))
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
             + '<div style="display: flex; gap: 14px">' + card("Energy", pe, chip("illustrative"), 12, "flex-grow: 1") + card("Max force", pf, chip("target 0.02"), 12, "flex-grow: 1") + '</div>'
             + '<div style="display: flex; gap: 14px; flex-grow: 1; min-height: 0">'
             + card("Constraints &amp; restraints", col(cr, row(btn("Add from selection", ic="plus", small=True), btn("Pick in Studio", ic="cursor", small=True), gap=6), gap=8), extra="flex-grow: 1")
             + card("Stop when any is met", col(kv("|F|<sub>max</sub>", "&lt; 0.02 kcal/mol/Å"), kv("ΔE per step", "&lt; 1e−8 relative"), kv("Iterations", "10 000"), kv("Box stress", "&lt; 10 bar"),
                                                 f'<div style="font-size: 11.5px; color: {MUTED}; line-height: 1.45">Criteria and the one that stopped the run go into the provenance manifest.</div>', gap=8), extra="width: 340px; flex-shrink: 0")
             + '</div></div>')
    head = pagehead("Relax · CAPS Relax", "Energy minimisation with the same force field as Dynamics",
                    row(btn("Copy as Python", ic="copy"), btn("Minimise", True, "play"), gap=8))
    return shell("Relax", "CAPS Relax — minimisation", f'{icon("chevr", 12, DIM)}<span>Relax</span>', head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>',
                 ("<span>PS_atactic_DP40 · 12 840 atoms</span>", "<span>Double precision</span>"))


# ------------------------------------------------------------------ Dynamics
def dynamics():
    ens = seg(["NVE", "NVT", "NPT", "NPH"], "NPT", cols=4)
    deck = [("#", "LAMMPS input written by CAPS 0.1.0 · must match CAPS step-0 energy"),
            ("", "units           real"), ("", "pair_style      lj/cut/coul/long 12.0"),
            ("", "pair_modify     mix arithmetic tail yes"), ("", "kspace_style    pppm 1.0e-5"), ("", "read_data       PS_atactic_DP40.data"),
            ("", "fix             1 all npt temp 298 298 100 iso 1 1 1000")]
    dk = "".join(f'<div style="color: {DIM if p == "#" else TEXT}">{("# " if p == "#" else "") + esc(t)}</div>' for p, t in deck)
    left = (f'<div style="width: 600px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 0 18px 22px">'
            + card("Ensemble &amp; integrator", col(ens, row(select("Integrator", "Velocity Verlet"), select("Multiple time step", "r-RESPA · off"), field("Δt", "1.0", "fs"), gap=8),
                                                    row(select("Thermostat", "Bussi–Donadio–Parrinello (CSVR)"), field("T", "298", "K"), field("τ<sub>T</sub>", "0.1", "ps"), gap=8),
                                                    row(select("Barostat", "Stochastic cell rescaling"), field("P", "1", "atm"), field("τ<sub>P</sub>", "2", "ps"), gap=8),
                                                    cite("Bussi et al., <i>J. Chem. Phys.</i> 126, 014101 (2007); Bernetti &amp; Bussi, <i>J. Chem. Phys.</i> 153, 114107 (2020)"), gap=10))
            + card("Interactions", col(row(select("Electrostatics", "Smooth PME"), field("r<sub>c</sub>", "12.0", "Å"), field("Ewald tol.", "1e−5"), gap=8),
                                       row(select("van der Waals", "Cut-off + tail correction"), select("Mixing", "Force-field default"), gap=8),
                                       row(select("Constraints", "LINCS · bonds to H"), select("Precision", "Mixed · deterministic"), gap=8), gap=10))
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
             + card("Pre-flight", col(*[f'<div style="display: flex; align-items: center; gap: 10px; font-size: 12.5px">{icon(i, 16, c)}<span style="flex-grow: 1">{t}</span><span style="font-family: {MONO}; font-size: 11.5px; color: {c}">{s}</span></div>'
                                        for i, t, s, c in (("check", "All atoms typed, no missing parameters", "ok", OK), ("check", "Net charge 0.000 e", "ok", OK),
                                                           ("check", "Box ≥ 2 r<sub>c</sub> in every direction", "ok", OK), ("alert", "Δt 1 fs with LINCS on H only", "check", WARN),
                                                           ("check", "NVE probe drift under threshold", "ok", OK))], gap=10), chip("4 / 5 ok"))
             + card("Run", row(field("Length", "10", "ns"), select("Run on", "hpc-login2 · 4 × A100 (gRPC)"), select("Checkpoint", "every 50 ps"), gap=8))
             + card("Export to other engines", col(row(seg(["LAMMPS", "GROMACS", "Both"], "LAMMPS"), '<div style="flex-grow: 1"></div>', chip(f'{dot(OK)} step-0 parity passed', OK, "#16261A"), gap=8),
                                                   f'<div style="padding: 12px 14px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px; font-family: {MONO}; font-size: 11.5px; line-height: 1.7">{dk}</div>',
                                                   row(btn("Save input deck", ic="download", small=True), btn("Compare energies", ic="chart", small=True), gap=6), gap=10), extra="flex-grow: 1")
             + '</div>')
    head = pagehead("Dynamics · CAPS Dynamics", "Molecular dynamics in CAPS, or a validated input deck for LAMMPS and GROMACS",
                    row(btn("Copy as Python", ic="copy"), btn("Queue run", True, "play"), gap=8))
    return shell("Dynamics", "CAPS Dynamics — MD setup", f'{icon("chevr", 12, DIM)}<span>Dynamics</span>', head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>',
                 ("<span>12 840 atoms · throughput [measure]</span>", "<span>Units: real (kcal/mol, Å, fs) at I/O · internal SI-based</span>"))


# ------------------------------------------------------------------ React
def react():
    # conversion vs time and largest-cluster fraction vs conversion (illustrative)
    conv = [(t, 0.9 * (1 - math.exp(-t / 1.8))) for t in [i * 0.05 for i in range(0, 160)]]
    pc = 1 / math.sqrt(1 * (2 - 1) * (4 - 1))
    clus = []
    for i in range(0, 91):
        a = i / 100
        v = 0.02 + 0.04 * a if a < pc - 0.02 else min(1.0, 0.1 + 2.6 * (a - pc + 0.02))
        clus.append((a, v))
    p1 = plot(336, 210, [(conv, ACC, 1.8, None), ([(0, 0.85), (8, 0.85)], OK, 1.2, "5 4")], (0, 8), (0, 1), [0, 2, 4, 6, 8], [0, 0.25, 0.5, 0.75, 1], "reaction cycle (×100)", "conversion α", pad=(44, 12, 18, 30))
    p2 = plot(336, 210, [(clus, SEL, 1.8, None), ([(pc, 0), (pc, 1)], ACC, 1.2, "4 3")], (0, 0.9), (0, 1), [0, 0.3, 0.6, 0.9], [0, 0.5, 1], "conversion α", "largest cluster fraction", pad=(44, 12, 18, 30))
    tmpl = (f'<div style="display: flex; align-items: center; gap: 12px">'
            f'<div style="flex-grow: 1; padding: 10px 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px; display: flex; flex-direction: column; gap: 4px">'
            f'<span style="font-size: 11px; color: {DIM}">PRE</span><span style="font-family: {MONO}; font-size: 12px">[CH2:1]1[O:2][CH:3]1R · H[N:4]HR′</span><span style="font-size: 11.5px; color: {MUTED}">epoxide + primary amine</span></div>'
            f'{icon("chevr", 20, ACC)}'
            f'<div style="flex-grow: 1; padding: 10px 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px; display: flex; flex-direction: column; gap: 4px">'
            f'<span style="font-size: 11px; color: {DIM}">POST</span><span style="font-family: {MONO}; font-size: 12px">R′HN:4–[CH2:1][CH:3](O:2H)R</span><span style="font-size: 11.5px; color: {MUTED}">β-hydroxy amine</span></div></div>')
    left = (f'<div style="width: 560px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 0 18px 22px">'
            + card("Reaction templates", col(tmpl, row(chip("atom-mapped · 4 reacting atoms"), chip("retype after bond"), chip("second step: secondary amine"), gap=6),
                                             row(btn("Edit template", ic="sliders", small=True), btn("Import template", ic="download", small=True), gap=6), gap=10), chip("2 templates"))
            + card("Protocol", col(radio_list([("REACTER-style during MD", "Gissinger, Jensen &amp; Wise, <i>Polymer</i> 128, 211 (2017)", True),
                                               ("Polymatic cycle", "Abbott, Hart &amp; Colina, <i>Theor. Chem. Acc.</i> 132, 1334 (2013)", False)]),
                                   row(field("Capture distance", "4.5", "Å"), field("Probability", "0.5"), field("Relax per step", "5", "ps"), gap=8),
                                   row(field("Target conversion", "85", "%"), select("Stop at", "Target or 800 cycles"), gap=8), gap=10))
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
             + '<div style="display: flex; gap: 14px">' + card("Conversion", p1, chip("illustrative"), 12, "flex-grow: 1") + card("Gel point", p2, chip("illustrative"), 12, "flex-grow: 1") + '</div>'
             + '<div style="display: flex; gap: 14px; flex-grow: 1; min-height: 0">'
             + card("Flory–Stockmayer check", col(formula(f'α<sub>c</sub> = 1 / √( r (f<sub>A</sub> − 1)(f<sub>B</sub> − 1) )<br>r = 1, f<sub>A</sub> = 2 (DGEBA), f<sub>B</sub> = 4 (IPDA) → α<sub>c</sub> = {pc:.3f}'),
                                                   kv("Simulated gel point", "[result] from cluster analysis", False, MUTED),
                                                   cite("Flory, <i>J. Am. Chem. Soc.</i> 63, 3083 (1941); Stockmayer, <i>J. Chem. Phys.</i> 11, 45 (1943)"), gap=10), extra="flex-grow: 1")
             + card("System", col(kv("Resin", "DGEBA × 200"), kv("Hardener", "IPDA × 100"), kv("Stoichiometry r", "1.00"), kv("Force field", "GAFF2 · retyped per bond"), gap=8), extra="width: 300px; flex-shrink: 0")
             + '</div></div>')
    head = pagehead("React · crosslinking", "Form bonds by distance and probability, retype, relax, and track the gel point",
                    row(btn("Copy as Python", ic="copy"), btn("Queue crosslinking", True, "play"), gap=8))
    return shell("React", "CAPS React — crosslinking", f'{icon("chevr", 12, DIM)}<span>React</span>', head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>',
                 ("<span>Epoxy network · 300 molecules</span>", "<span>Cluster analysis every 10 cycles</span>"), ["DGEBA_IPDA.caps"])


# ------------------------------------------------------------------ Information architecture
def ia():
    def node(title, sub_, color=LINE, w=None, accent=False):
        ww = f"width: {w}px; " if w else ""
        return (f'<div style="{ww}display: flex; flex-direction: column; gap: 4px; padding: 12px 14px; background: {BG1}; border: 1px solid {color}; border-radius: 9px; {"box-shadow: inset 3px 0 0 " + color + "; " if accent else ""}">'
                f'<span style="font-size: 13.5px; font-weight: 600">{title}</span><span style="font-size: 11.5px; color: {MUTED}; line-height: 1.45">{sub_}</span></div>')
    arrow = f'<div style="display: flex; align-items: center; justify-content: center; width: 26px; flex-shrink: 0">{icon("chevr", 18, DIM)}</div>'
    pipeline = [("Studio", "build · edit · type"), ("Grow", "chains into cells"), ("Pack", "molecules into regions"), ("Relax", "minimise"),
                ("Equilibrate", "published protocols"), ("Dynamics", "MD or export"), ("React", "crosslink"), ("Analyze", "properties")]
    pl = arrow.join(node(a, b, ACC if a == "Studio" else LINE, 132) for a, b in pipeline)
    cross = [("Field", "Typing, parameters and charges used by every engine. Missing parameters block runs."),
             ("Jobs", "Every long operation: progress, logs, checkpoints, local or remote."),
             ("Bench", "Validation suite that regenerates Tables 1–12."),
             ("Project", "Documents, recipes, selection sets, provenance manifests.")]
    cr = "".join(node(a, b, SEL, None, True) for a, b in cross)
    layers = [("Adapters", "Avalonia GUI · CLI · Python (nanobind) · gRPC daemon", SEL),
              ("Command bus", "Every user action is a named command with arguments — recorded, replayable, testable", ACC),
              ("Undo stack", "Commands produce inverse deltas; history is editable and saved with the document", ACC),
              ("Document model", "Atoms (SoA), bonds, cell, residues, types, selections, monitors", TEXT),
              ("Engines (C ABI)", "Chem · Field · Grow · Pack · Relax · Dynamics · Equilibrate · React · Analyze · IO", MUTED)]
    ly = "".join(f'<div style="display: flex; align-items: center; gap: 14px; padding: 12px 16px; background: {BG1}; border: 1px solid {LINE}; border-radius: 8px">'
                 f'<span style="width: 10px; height: 10px; border-radius: 2px; background: {c}; flex-shrink: 0"></span><span style="width: 170px; flex-shrink: 0; font-weight: 600">{a}</span>'
                 f'<span style="font-size: 12.5px; color: {MUTED}">{b}</span></div>' for a, b, c in layers)
    ex = [("GUI", "drag bond C12 → C31"), ("CLI", "caps edit bond.add 12 31 --order 1"), ("Python", "doc.bond.add(12, 31, order=1)")]
    exr = "".join(f'<div style="display: flex; gap: 10px; font-family: {MONO}; font-size: 12px"><span style="width: 56px; color: {DIM}">{a}</span><span>{esc(b)}</span></div>' for a, b in ex)
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 22px; padding: 36px 44px">'
            f'<div style="display: flex; flex-direction: column; gap: 6px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">How CAPS is organised</h1>'
            f'<p style="margin: 0; font-size: 14px; color: {MUTED}">The module rail follows the physical workflow. Four areas cut across every module. Every action runs through one command layer.</p></div>'
            + section("Workflow modules · left to right on the rail", f'<div style="display: flex; align-items: stretch">{pl}</div>', pad=0)
            + section("Cross-cutting areas", f'<div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px">{cr}</div>', pad=0)
            + '<div style="display: flex; gap: 28px; flex-grow: 1">'
            + f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 10px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Command architecture · top calls down</h3>{ly}</div>'
            + f'<div style="width: 420px; flex-shrink: 0; display: flex; flex-direction: column; gap: 10px"><h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">One action, three front ends</h3>'
              f'<div style="display: flex; flex-direction: column; gap: 10px; padding: 16px; background: {BG1}; border: 1px solid {LINE}; border-radius: 8px">{exr}'
              f'<div style="height: 1px; background: {LINE}"></div><div style="font-size: 12.5px; color: {MUTED}; line-height: 1.5">All three emit the same <span style="font-family: {MONO}; color: {TEXT}">bond.add</span> command. Undo reverts it the same way, and the history panel shows it once.</div></div></div>'
            + '</div></div>')
    return page("CAPS — information architecture", body)


# ------------------------------------------------------------------ Interaction map
def interactions():
    def key(k):
        return f'<kbd style="display: inline-flex; align-items: center; justify-content: center; min-width: 24px; height: 24px; padding: 0 7px; background: {BG2}; border: 1px solid {LINE}; border-bottom-width: 2px; border-radius: 5px; font-family: {MONO}; font-size: 11.5px; color: {TEXT}">{k}</kbd>'

    def combo(ks):
        return '<span style="display: inline-flex; gap: 4px; align-items: center">' + "".join(key(k) if k != "+" else f'<span style="color: {DIM}">+</span>' for k in ks) + '</span>'

    def group(title, rows, w=None):
        rr = "".join(f'<div style="display: flex; align-items: center; justify-content: space-between; gap: 12px; min-height: 30px; border-bottom: 1px solid {BG2}"><span style="font-size: 12.5px; color: {TEXT}">{a}</span>{b}</div>' for a, b in rows)
        return (f'<section style="display: flex; flex-direction: column; gap: 6px; padding: 16px; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px">'
                f'<h2 style="margin: 0 0 4px; font-size: 13px; font-weight: 600">{title}</h2>{rr}</section>')

    def mouse(label):
        return f'<span style="font-family: {MONO}; font-size: 11.5px; color: {SEL}">{label}</span>'

    view = group("View · any mode", [("Rotate", mouse("left drag on empty space")), ("Pan", mouse("middle drag · ⇧ + left drag")), ("Zoom", mouse("wheel · pinch")),
                                     ("Roll about view axis", mouse("⌥ + left drag")), ("Centre on atom", mouse("double-click atom")), ("Fit all", combo(["F"])),
                                     ("Axis views", combo(["1", "2", "3"])), ("Perspective / orthographic", combo(["5"]))])
    sel = group("Select mode · S", [("Select atom", mouse("click")), ("Add to selection", mouse("⇧ + click")), ("Toggle atom", mouse("⌘ + click")),
                                    ("Box select", mouse("left drag from empty space with ⇧")), ("Lasso", combo(["L"])), ("Molecule", mouse("triple-click")),
                                    ("Grow along bonds", combo(["]"])), ("Invert selection", combo(["⌘", "+", "I"]))])
    build = group("Build mode · B", [("Place atom", mouse("click empty space")), ("Draw bond", mouse("drag atom → atom")), ("Grow chain", mouse("drag atom → empty")),
                                     ("Cycle bond order", mouse("click bond")), ("Change element", combo(["C", "N", "O", "H"])), ("Periodic table", combo(["⇧", "+", "E"])),
                                     ("Attach fragment", mouse("drop fragment on open valence")), ("Add hydrogens", combo(["⌘", "+", "H"]))])
    edit = group("Edit &amp; measure", [("Move selection", combo(["G"])), ("Rotate selection", combo(["R"])), ("Constrain to axis", combo(["X", "Y", "Z"])),
                                        ("Measure", combo(["M"]) + ' ' + mouse(" then click 2–4 atoms")), ("Pin as monitor", combo(["P"])), ("Set exact value", mouse("double-click a measurement")),
                                        ("Invert stereocentre", combo(["⌥", "+", "I"])), ("Delete", combo(["⌫"]))])
    app = group("App", [("Command palette", combo(["⌘", "+", "K"])), ("Undo / redo", combo(["⌘", "+", "Z"]) + ' ' + combo(["⇧", "⌘", "+", "Z"])),
                        ("Clean geometry now", combo(["⌘", "+", "⇧", "+", "C"])), ("Toggle auto-clean", combo(["A"])), ("Record macro", combo(["⌘", "+", "⇧", "+", "R"])),
                        ("Switch module", combo(["⌘", "+", "1…9"])), ("Next document", combo(["⌃", "+", "Tab"])), ("Help for this panel", combo(["F1"]))])
    devices = group("Other devices", [("Trackpad rotate", mouse("two-finger drag")), ("Trackpad pan", mouse("⇧ + two-finger drag")), ("Trackpad zoom", mouse("pinch")),
                                      ("3D mouse", mouse("6-DOF: translate + rotate view")), ("Pen tablet", mouse("same as mouse; barrel = middle")),
                                      ("Keyboard only", mouse("Tab through atoms, arrows move focus")), ("Remap", mouse("Settings → Input · export as JSON")), ("Conflicts", mouse("flagged on save"))])
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 22px; padding: 36px 44px">'
            f'<div style="display: flex; align-items: flex-end; justify-content: space-between; gap: 20px"><div style="display: flex; flex-direction: column; gap: 6px"><h1 style="margin: 0; font-size: 28px; font-weight: 600; letter-spacing: -0.015em">Mouse &amp; keyboard map</h1>'
            f'<p style="margin: 0; font-size: 14px; color: {MUTED}">Two modes carry the editing: Select (S) and Build (B). View gestures work the same in every mode. Every shortcut can be remapped.</p></div>'
            f'<div style="display: flex; gap: 8px">{chip(mouse("cyan") + " = mouse or gesture")}{chip(key("K") + " = key")}</div></div>'
            f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 16px">{view}{sel}{build}{edit}{app}{devices}</div></div>')
    return page("CAPS — interaction map", body)
