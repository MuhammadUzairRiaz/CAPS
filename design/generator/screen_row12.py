import math, random
from lib import *
from mols import *
from screen_app import card, pagehead, shell
from screen_builders import builder_shell, footer, hud, VPW, VPH
from screen_studio2 import spacer, bar


def tube_cell(sid, W, H, groups, L, seed0=100, persistence=0.6, n=55, width=0.55):
    sc = Scene(sid, W, H, yaw=0.6, pitch=0.42, persp=0.3, fog=0.55)
    k = 0
    for color, count in groups:
        for _ in range(count):
            rnd = random.Random(seed0 + k)
            st = (rnd.uniform(3, L - 3), rnd.uniform(3, L - 3), rnd.uniform(3, L - 3))
            sc.add_tube(random_chain(st, n, seed0 + 50 + k, box=((0.8,) * 3, (L - 0.8,) * 3), persistence=persistence), color, width)
            k += 1
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.2, None, 0.8)
    return sc


# ------------------------------------------------------------------ Blend builder
def blend():
    sc = tube_cell("bl", VPW, VPH, [(ACC, 8), (SEL, 8)], 36.0, seed0=900)
    svg = sc.svg()
    comps = table(["Component", "wt %", "Chains", "DP"], [[f'{dot(ACC)} PS atactic', "50", "8", "40"], [f'{dot(SEL)} PMMA atactic', "50", "[from wt %]", "40"]],
                  ["44%", "16%", "22%", "18%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=12)
    morph = [("Mixed", "random, homogeneous start"), ("Two slabs", "pre-separated, for interfaces"), ("Droplet", "minor phase as a sphere")]
    mg = "".join(f'<button style="display: flex; flex-direction: column; gap: 3px; padding: 10px; text-align: left; border-radius: 7px; border: 1px solid {ACC if n == "Mixed" else LINE}; background: {BG3 if n == "Mixed" else BG0}; cursor: pointer"><span style="font-size: 12.5px; color: {TEXT}">{n}</span><span style="font-size: 11px; color: {DIM}">{d}</span></button>' for n, d in morph)
    panel = (panel_head("Blend builder", chip("components by weight"))
             + section("Components", col(comps, row(btn("Add component", ic="plus", small=True), btn("From library", ic="search", small=True), gap=6), gap=10))
             + section("Composition", col(seg(["Weight %", "Volume %", "Chain count"], "Weight %", True), f'<div style="font-size: 11.5px; color: {MUTED}">PMMA chain count is computed from the weight fraction and each chain\'s molar mass.</div>', gap=8))
             + section("Initial morphology", f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 6px">{mg}</div>')
             + section("Cell", col(row(field("Target density", "1.10", "g/cm³", note="from component densities"), select("Growth", "CBMC · real force field"), gap=10), gap=8))
             + footer(btn("Cancel"), btn("Build blend", True, "grow")))
    ov = hud(chip(f'{dot(ACC)} PS', TEXT, "#16191Ccc"), chip(f'{dot(SEL)} PMMA', TEXT, "#16191Ccc"), chip("mixed start · 16 chains shown", MUTED, "#16191Ccc"))
    return builder_shell("PS_PMMA_blend.caps", "Polymer", svg, panel, ov, ("<span>2 components · 50 / 50 wt %</span>", "<span>Miscibility (χ from CED) available in Analyze after equilibration</span>"))


# ------------------------------------------------------------------ Coarse-grained
def coarse():
    M, N, rho = 50, 100, 0.85
    L = (M * N / rho) ** (1 / 3)
    Ld = 14.0
    sc = Scene("cg", VPW, VPH, yaw=0.55, pitch=0.4, persp=0.3, fog=0.55, bond_w=0.18)
    cols = [ACC, SEL, "#E07A5F", "#9B7BD6", OK, "#D6A45E", "#E9ECEF", "#4C7BD9"]
    for k in range(10):
        rnd = random.Random(1200 + k)
        st = (rnd.uniform(1, Ld - 1), rnd.uniform(1, Ld - 1), rnd.uniform(1, Ld - 1))
        pts = random_chain(st, 45, 1300 + k, bond=0.97, box=((0.4,) * 3, (Ld - 0.4,) * 3), persistence=0.35)
        sc.add_atoms([{"e": "C", "p": p, "r": 0.42, "c": cols[k % len(cols)]} for p in pts], [(i, i + 1) for i in range(len(pts) - 1)])
    sc.add_box((0, 0, 0), (Ld, 0, 0), (0, Ld, 0), (0, 0, Ld), MUTED, 1.2, None, 0.8)
    svg = sc.svg()
    mp = table(["Bead", "Atoms", "Mass (amu)"], [["B1", "CH₂–CH(C₆H₅) · 1 styrene unit", "104.15"]], ["16%", "58%", "26%"], mono_cols=(2,), align_right=(2,), fs=12)
    panel = (panel_head("Coarse-grained model", chip("Kremer–Grest"))
             + section("Chains", col(row(field("Chains M", str(M)), field("Beads per chain N", str(N)), gap=10), row(field("Density ρσ³", f"{rho}"), field("Box L", f"{L:.2f}", "σ"), gap=10), gap=10))
             + section("Potential", col(kv("Bonds", "FENE · K = 30 ε/σ², R₀ = 1.5 σ"), kv("Non-bonded", "WCA · cut at 2<sup>1/6</sup> σ"), row(field("Bending k<sub>θ</sub>", "0", "ε"), select("Units", "Reduced (ε, σ, m)"), gap=10),
                                        cite("Kremer &amp; Grest, <i>J. Chem. Phys.</i> 92, 5057 (1990)"), gap=8))
             + section("Backmap to all atoms", col(mp, row(select("Insert fragments", "Rigid unit + CAPS Relax"), select("Force field", "GAFF2"), gap=10), gap=10))
             + footer(btn("Cancel"), btn("Build CG melt", True, "grow")))
    ov = hud(chip(f"{M} × {N} beads · preview shows 10 chains", TEXT, "#16191Ccc"), chip("bond ≈ 0.97 σ", MUTED, "#16191Ccc"))
    return builder_shell("KG_melt_50x100.caps", "Polymer", svg, panel, ov, (f"<span>{M * N} beads · L = {L:.2f} σ (from ρσ³ = {rho})</span>", "<span>Box length computed from M·N/ρ</span>"))


# ------------------------------------------------------------------ Reaction template editor
def reaction_editor():
    def atom(x, y, lab, mapn, color=TEXT, hl=False):
        ring = f'<circle cx="{x}" cy="{y}" r="17" fill="none" stroke="{SEL}" stroke-width="2"></circle>' if hl else ""
        return (f'{ring}<text x="{x}" y="{y + 6}" text-anchor="middle" fill="{color}" font-size="17" font-family="IBM Plex Sans" font-weight="600">{lab}</text>'
                f'<circle cx="{x + 16}" cy="{y - 15}" r="9" fill="{BG3}" stroke="{ACC}"></circle><text x="{x + 16}" y="{y - 11}" text-anchor="middle" fill="{ACC}" font-size="10" font-family="IBM Plex Mono">{mapn}</text>')

    def bond(a, b, dash=False, color=TEXT, shorten=13):
        dx, dy = b[0] - a[0], b[1] - a[1]
        l = math.hypot(dx, dy)
        ux, uy = dx / l, dy / l
        d = ' stroke-dasharray="5 4"' if dash else ""
        return f'<line x1="{a[0] + ux * shorten:.1f}" y1="{a[1] + uy * shorten:.1f}" x2="{b[0] - ux * shorten:.1f}" y2="{b[1] - uy * shorten:.1f}" stroke="{color}" stroke-width="2.2"{d}></line>'
    O = "#E5534B"
    N = "#6F9BF0"
    C1, O2, C3, N4, H5, H6 = (160, 210), (205, 140), (250, 210), (360, 300), (410, 250), (330, 350)
    pre = (bond(C1, O2) + bond(O2, C3) + bond(C1, C3) + bond(C3, (320, 190)) + bond(N4, H5) + bond(N4, H6) + bond(N4, (420, 340))
           + f'<path d="M {N4[0] - 20} {N4[1] - 10} Q 250 300 {C1[0] + 10} {C1[1] + 16}" fill="none" stroke="{ACC}" stroke-width="1.8" stroke-dasharray="6 4"></path>'
           + atom(*C1, "C", 1, hl=True) + atom(*O2, "O", 2, O, True) + atom(*C3, "C", 3) + atom(*N4, "N", 4, N, True) + atom(*H5, "H", 5, hl=True) + atom(*H6, "H", 6)
           + f'<text x="330" y="195" fill="{MUTED}" font-size="15" font-family="IBM Plex Sans">R</text><text x="428" y="352" fill="{MUTED}" font-size="15" font-family="IBM Plex Sans">R′</text>')
    pC1, pC3, pO2, pN4, pH5, pH6 = (210, 230), (290, 190), (290, 110), (130, 190), (360, 80), (80, 150)
    post = (bond(pN4, pC1, color=OK) + bond(pC1, pC3) + bond(pC3, pO2) + bond(pO2, pH5, color=OK) + bond(pC3, (360, 230)) + bond(pN4, pH6) + bond(pN4, (100, 260))
            + atom(*pC1, "C", 1, hl=True) + atom(*pO2, "O", 2, O, True) + atom(*pC3, "C", 3) + atom(*pN4, "N", 4, N, True) + atom(*pH5, "H", 5, hl=True) + atom(*pH6, "H", 6)
            + f'<text x="366" y="240" fill="{MUTED}" font-size="15" font-family="IBM Plex Sans">R</text><text x="72" y="276" fill="{MUTED}" font-size="15" font-family="IBM Plex Sans">R′</text>')

    def pane(title, svgbody, tag):
        return (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}; border: 1px solid {LINE}; border-radius: 10px; overflow: hidden">'
                f'<div style="height: 38px; display: flex; align-items: center; gap: 8px; padding: 0 12px; background: {BG1}; border-bottom: 1px solid {LINE}"><span style="font-size: 12.5px; font-weight: 600">{title}</span>{chip(tag)}</div>'
                f'<svg width="100%" height="420" viewBox="0 0 480 420" role="img" aria-label="{title} structure">{svgbody}</svg></div>')
    changes = table(["Change", "Atoms", ""], [["bond formed", "N4 – C1", f'<span style="color: {OK}">new</span>'], ["bond broken", "C1 – O2", f'<span style="color: {ERR}">ring opens</span>'],
                                             ["bond broken", "N4 – H5", f'<span style="color: {ERR}">proton leaves N</span>'], ["bond formed", "O2 – H5", f'<span style="color: {OK}">hydroxyl</span>']],
                    ["38%", "30%", "32%"], mono_cols=(1,), fs=12)
    checks = col(*[f'<div style="display: flex; gap: 8px; align-items: center; font-size: 12.5px">{icon("check", 15, OK)}<span>{t}</span></div>' for t in
                   ("Every mapped atom appears on both sides (6 of 6)", "Element of each mapped atom unchanged", "Net charge 0 → 0", "Edge atoms R, R′ keep their bonds")], gap=8)
    left = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 0 18px 22px">'
            + f'<div style="display: flex; gap: 14px; align-items: center">{pane("Pre-reaction", pre, "epoxide + primary amine")}{icon("chevr", 26, ACC)}{pane("Post-reaction", post, "β-hydroxy amine")}</div>'
            + row(chip(f'<span style="width: 14px; height: 14px; border-radius: 50%; border: 1px solid {ACC}; display: inline-block"></span> map number'), chip(f'<span style="width: 14px; height: 14px; border-radius: 50%; border: 2px solid {SEL}; display: inline-block"></span> reacting atom'), chip("drag between panes to map atoms"), gap=6)
            + '</div>')
    right = (f'<div style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 18px 22px">'
             + card("Bond changes", changes, chip("4"))
             + card("Checks", checks, chip("pass", OK, "#16261A"))
             + card("After reaction", col(kv("Retype", "atoms within 3 bonds"), kv("Charges", "from template, neutral"), kv("Relax", "5 ps constrained"), gap=7))
             + '</div>')
    content = (pagehead("Reaction template · epoxy–amine (primary)", "Atom-mapped pre/post templates used by CAPS React",
                        row(btn("Import template JSON", ic="download"), btn("Test on a pair", ic="play"), btn("Save template", True, "save"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>')
    return shell("React", "CAPS React — template editor", f'{icon("chevr", 12, DIM)}<span>React</span>', content,
                 ("<span>Template 1 of 2 · second step: secondary amine</span>", "<span>Saved as atom-mapped reaction SMARTS + JSON</span>"))


# ------------------------------------------------------------------ Large system LOD
def large():
    W, H = 1440 - 72 - 360, 900 - 44 - 48 - 26
    N = 1_024_000
    rnd = random.Random(2)
    L = 220.0
    pts = []
    for _ in range(5200):
        pts.append((rnd.uniform(0, L), rnd.uniform(0, L), rnd.uniform(0, L)))
    yaw, pitch = 0.6, 0.42
    c = (L / 2,) * 3
    proj = []
    for p in pts:
        r = rot(sub(p, c), yaw, pitch)
        proj.append(r)
    zmin = min(r[2] for r in proj)
    zmax = max(r[2] for r in proj)
    s = min(W, H) / (L * 1.9)
    dots = []
    for r in sorted(proj, key=lambda r: r[2]):
        t = (r[2] - zmin) / (zmax - zmin)
        k = 1 + 0.3 * (t - 0.5)
        col_ = mix("#8E959C", BG0, 0.75 - 0.6 * t)
        dots.append(f'<circle cx="{W / 2 + r[0] * s * k:.1f}" cy="{H / 2 - r[1] * s * k:.1f}" r="{0.9 + 0.9 * t:.2f}" fill="{col_}"></circle>')
    box = Scene("lg", W, H, yaw=yaw, pitch=pitch, persp=0.3)
    corners = [(0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), (L, L, 0), (L, 0, L), (0, L, L), (L, L, L)]
    edges = [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5), (2, 4), (2, 6), (3, 5), (3, 6), (4, 7), (5, 7), (6, 7)]
    lines = ""
    for i, j in edges:
        a = rot(sub(corners[i], c), yaw, pitch)
        b = rot(sub(corners[j], c), yaw, pitch)
        lines += f'<line x1="{W / 2 + a[0] * s:.1f}" y1="{H / 2 - a[1] * s:.1f}" x2="{W / 2 + b[0] * s:.1f}" y2="{H / 2 - b[1] * s:.1f}" stroke="{MUTED}" stroke-width="1.2" stroke-opacity="0.8"></line>'
    # near-field inset: impostor detail
    atoms, _, _ = polystyrene(3, "atactic", seed=4, curve=False)
    ins = Scene("li", 300, 190, yaw=0.3, pitch=0.5, persp=0.2, fog=0.3)
    ins.add_atoms(atoms)
    inset = (f'<div style="position: absolute; right: 16px; bottom: 16px; width: 300px; background: {BG1}; border: 1px solid {ACC}; border-radius: 10px; overflow: hidden">'
             f'<div style="height: 30px; display: flex; align-items: center; gap: 8px; padding: 0 10px; border-bottom: 1px solid {LINE}; font-size: 11.5px; color: {MUTED}">{icon("search", 13, ACC)}Near field · full impostors</div>{ins.svg()}</div>')
    svg = f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="Million-atom cell" style="display: block">{lines}{"".join(dots)}</svg>'
    perf = (f'<div style="position: absolute; left: 14px; top: 12px; width: 260px; padding: 10px 12px; display: flex; flex-direction: column; gap: 5px; background: #16191Ce6; border: 1px solid {LINE}; border-radius: 8px; font-family: {MONO}; font-size: 11.5px">'
            f'<div style="display: flex; justify-content: space-between"><span style="color: {MUTED}">atoms</span><span>1 024 000</span></div>'
            f'<div style="display: flex; justify-content: space-between"><span style="color: {MUTED}">frame time</span><span style="color: {DIM}">[measure]</span></div>'
            f'<div style="display: flex; justify-content: space-between"><span style="color: {MUTED}">target</span><span>60 fps · 16.7 ms</span></div>'
            f'<div style="display: flex; justify-content: space-between"><span style="color: {MUTED}">LOD</span><span style="color: {ACC}">points &gt; 40 Å</span></div></div>')
    vp = f'<main aria-label="3D view" style="position: relative; flex-grow: 1; min-width: 0; background: {BG0}">{svg}{perf}{inset}</main>'
    bytes_per = 12 + 4 + 4
    mem = N * bytes_per / 1e6
    right = (panel_head("Level of detail", chip("Table 11 target"))
             + section("Tiers", col(kv("Near (&lt; 40 Å)", "ray-cast sphere + cylinder impostors"), kv("Mid", "sphere impostors, no bonds"), kv("Far", "point sprites"), kv("Culling", "frustum + occlusion by cell"), gap=8))
             + section("GPU memory estimate", col(kv("Per atom", f"{bytes_per} B · xyz float32, colour, radius"), kv("1 024 000 atoms", f"{mem:.1f} MB"), f'<div style="font-size: 11.5px; color: {DIM}">Computed from the buffer layout; trajectories stream frames instead of loading them all.</div>', gap=8))
             + section("Quality", col(toggle("Ambient occlusion (near tier only)", True), toggle("Outlines", True), toggle("Adaptive: drop tiers to keep 60 fps", True), gap=9))
             + footer(btn("Benchmark this view", ic="bench"), btn("Apply", True, "check")))
    body = (topbar(["PE_melt_1M.caps"], "PE_melt_1M.caps") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">'
            + bar(tbtn("cursor", "Select", True), tbtn("rotate", "Rotate"), sep(), tbtn("eye", "Style", False, "Adaptive LOD"), tbtn("cube", "Cell", True, "Cell 220 Å"), spacer(), tbtn("gauge", "Performance HUD", True, "Performance HUD"))
            + f'<div style="flex-grow: 1; display: flex; min-height: 0">{vp}<aside aria-label="LOD" style="width: 360px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{right}</aside></div></div></div>'
            + statusbar("<span>1 024 000 atoms · 5 200 sample points drawn in this mockup</span>", "<span>Load target &lt; 5 s · edit target &lt; 16 ms at 10 k atoms</span>"))
    return page("CAPS Studio — million-atom view", body)


# ------------------------------------------------------------------ CLI
def cli():
    lines = [
        ("p", "caps --help"),
        ("o", "caps 0.1.0 — Chain Assembly and Packing Suite"),
        ("o", ""),
        ("o", "usage: caps <command> [options]"),
        ("o", ""),
        ("h", "build       make monomers, chains, crystals, surfaces, blends"),
        ("h", "type        assign a force field and report why each atom got its type"),
        ("h", "grow        grow chains into a periodic cell (Theodorou–Suter, CBMC)"),
        ("h", "pack        pack molecules into regions (no packmol)"),
        ("h", "relax       minimise energy (SD, CG, L-BFGS, FIRE)"),
        ("h", "md          run molecular dynamics"),
        ("h", "equilibrate run a published protocol (larsen21, pushoff, annealing)"),
        ("h", "react       crosslink with atom-mapped templates"),
        ("h", "analyze     density, RDF, S(q), Rg, C∞, Tg, elastic constants, MSD"),
        ("h", "export      write LAMMPS, GROMACS, PDB, mol2, CIF, xyz …"),
        ("h", "run         run a YAML/JSON recipe"),
        ("h", "bench       regenerate validation Tables 1–12"),
        ("o", ""),
        ("p", "caps run ps_cell.yaml --seed 20260923 --deterministic"),
        ("s", "[1/6] build        PS · DP 40 × 20 chains · 12 840 atoms                 done"),
        ("s", "[2/6] type         GAFF2 · 12 840 typed · 0 missing parameters            done"),
        ("s", "[3/6] grow         CBMC k=16 · chain 20/20 · close contacts 0              done"),
        ("s", "[4/6] relax        L-BFGS · |F|max 0.02 kcal/mol/Å                       done"),
        ("a", "[5/6] equilibrate  larsen21 · step 14/21  ████████████████░░░░░░░░  66 %"),
        ("d", "[6/6] export       lammps, gromacs                                        waiting"),
    ]
    out = []
    for kind, t in lines:
        if kind == "p":
            out.append(f'<div><span style="color: {OK}">~/ps-tacticity</span> <span style="color: {ACC}">$</span> <span style="color: {TEXT}">{esc(t)}</span></div>')
        elif kind == "h":
            name, _, desc = t.partition(" ")
            out.append(f'<div style="white-space: pre">  <span style="color: {SEL}">{name:<11}</span> <span style="color: {MUTED}">{esc(desc.strip())}</span></div>')
        elif kind == "s":
            out.append(f'<div style="white-space: pre; color: {TEXT}">{esc(t[:-4])}<span style="color: {OK}">done</span></div>')
        elif kind == "a":
            out.append(f'<div style="white-space: pre; color: {ACC}">{esc(t)}</div>')
        elif kind == "d":
            out.append(f'<div style="white-space: pre; color: {DIM}">{esc(t)}</div>')
        else:
            out.append(f'<div style="white-space: pre; color: {MUTED}">{esc(t) or "&#160;"}</div>')
    term = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: #0A0B0C; border: 1px solid {LINE}; border-radius: 10px; overflow: hidden">'
            f'<div style="height: 34px; display: flex; align-items: center; gap: 10px; padding: 0 14px; background: {BG1}; border-bottom: 1px solid {LINE}; font-family: {MONO}; font-size: 12px; color: {MUTED}">{icon("terminal", 14, ACC)}zsh · caps</div>'
            f'<div style="padding: 16px 18px; font-family: {MONO}; font-size: 13px; line-height: 1.7">{"".join(out)}</div></div>')
    recipe = ["recipe: 1", "name: ps_cell", "build:", "  polymer:", "    smiles: \"*CC(*)c1ccccc1\"", "    dp: 40", "    chains: 20", "    tacticity: atactic",
              "type: { forcefield: gaff2 }", "grow: { method: cbmc, trials: 16, density: 0.50 }", "relax: { method: lbfgs, fmax: 0.02 }",
              "equilibrate: { protocol: larsen21 }", "export: [lammps, gromacs]"]
    rc = "".join(f'<div style="white-space: pre; color: {DIM if l.startswith("recipe") else TEXT}">{esc(l)}</div>' for l in recipe)
    side = (f'<div style="width: 440px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
            + card("ps_cell.yaml", f'<div style="font-family: {MONO}; font-size: 12px; line-height: 1.75">{rc}</div>', chip("recipe"))
            + card("CLI rules", col(*[kv(a, b, False) for a, b in (("Binary", "one static executable"), ("Exit codes", "0 ok · 2 input · 3 missing params · 4 failed run"),
                                                                  ("Output", "human by default · --json for machines"), ("Same commands as", "Studio, Python and the daemon"))], gap=8))
            + '</div>')
    body = (f'<div style="flex-grow: 1; display: flex; gap: 20px; padding: 28px 32px">'
            + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px"><h1 style="margin: 0; font-size: 24px; font-weight: 600; letter-spacing: -0.01em">caps · command line</h1>{term}</div>'
            + side + '</div>')
    return page("CAPS — command line", body)


# ------------------------------------------------------------------ Jupyter
def notebook():
    def cellbox(n, code_lines, output=""):
        code = "".join(f'<div style="white-space: pre">{l}</div>' for l in code_lines)
        o = f'<div style="display: flex; gap: 12px"><span style="width: 52px; flex-shrink: 0; text-align: right; font-family: {MONO}; font-size: 12px; color: {ERR}">[{n}]:</span><div style="flex-grow: 1; min-width: 0">{output}</div></div>' if output else ""
        return (f'<div style="display: flex; flex-direction: column; gap: 8px">'
                f'<div style="display: flex; gap: 12px"><span style="width: 52px; flex-shrink: 0; text-align: right; font-family: {MONO}; font-size: 12px; color: {SEL}; padding-top: 10px">[{n}]:</span>'
                f'<div style="flex-grow: 1; min-width: 0; padding: 10px 14px; background: {BG1}; border: 1px solid {LINE}; border-radius: 6px; font-family: {MONO}; font-size: 12.5px; line-height: 1.7">{code}</div></div>{o}</div>')

    def kw(t):
        return f'<span style="color: {SEL}">{t}</span>'

    def st(t):
        return f'<span style="color: {OK}">{esc(t)}</span>'
    c1 = [f'{kw("import")} caps', f'doc = caps.polymer({st(chr(34) + "*CC(*)c1ccccc1" + chr(34))}, dp={ACC_S("40")}, tacticity={st(chr(34) + "atactic" + chr(34))}, seed={ACC_S("5")})', 'doc.view(style=' + st('"ball-and-stick"') + ')']
    atoms, _, _ = polystyrene(7, "atactic", seed=5)
    sc = Scene("nb", 760, 250, yaw=-0.28, pitch=0.42, roll=0.05, persp=0.25, fog=0.55)
    sc.add_atoms(atoms)
    widget = (f'<div style="border: 1px solid {LINE}; border-radius: 8px; overflow: hidden; background: {BG0}; width: 760px">{sc.svg()}'
              f'<div style="height: 30px; display: flex; align-items: center; gap: 10px; padding: 0 10px; border-top: 1px solid {LINE}; font-size: 11.5px; color: {MUTED}">{icon("rotate", 13, DIM)}drag to rotate · CAPS View widget · {len(atoms)} atoms shown</div></div>')
    c2 = ["cells = [caps.sweep.result(t) " + kw("for") + " t " + kw("in") + " (" + st('"isotactic"') + ", " + st('"syndiotactic"') + ", " + st('"atactic"') + ")]", "caps.table(cells, [" + st('"density"') + ", " + st('"tg"') + ", " + st('"c_inf"') + "])"]
    tb = table(["", "density (g/cm³)", "tg (K)", "c_inf"], [["isotactic", "[result]", "[result]", "[result]"], ["syndiotactic", "[result]", "[result]", "[result]"], ["atactic", "[result]", "[result]", "[result]"]],
               ["28%", "24%", "24%", "24%"], mono_cols=(1, 2, 3), align_right=(1, 2, 3), fs=12)
    c3 = ["doc.provenance.citations(fmt=" + st('"bibtex"') + ")[:" + ACC_S("2") + "]"]
    bib = (f'<div style="font-family: {MONO}; font-size: 12px; line-height: 1.6; color: {MUTED}; white-space: pre">[\'@article{{bussi2007, title={{Canonical sampling through velocity rescaling}}, …}}\',\n'
           f' \'@article{{larsen2011, journal={{Macromolecules}}, volume={{44}}, pages={{6944}}, …}}\']</div>')
    body = (f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 18px; padding: 24px 40px; max-width: 1100px">'
            + f'<div style="display: flex; align-items: center; gap: 12px"><h1 style="margin: 0; font-size: 22px; font-weight: 600">ps_tacticity.ipynb</h1>{chip("Python 3.12 · caps 0.1.0")}{chip(f"{dot(OK)} kernel idle")}</div>'
            + cellbox(1, c1, widget) + cellbox(2, c2, f'<div style="width: 760px">{tb}</div>') + cellbox(3, c3, bib) + '</div>')
    return page("CAPS — Jupyter notebook", body)


def ACC_S(t):
    return f'<span style="color: {ACC}">{t}</span>'
