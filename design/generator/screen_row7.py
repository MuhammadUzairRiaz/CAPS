import math, random, re
from lib import *
from mols import *
from screen_app import card
from screen_builders import footer, hud
from screen_studio2 import bar, spacer, select_inline
import screen_main

SYMS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
        "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr "
        "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()


def pt_pos(z):
    """(row, col) in an 18-column table; f-block on rows 9 and 10 (1-based)."""
    if z == 1: return (1, 1)
    if z == 2: return (1, 18)
    for start, period in ((3, 2), (11, 3)):
        if start <= z < start + 8:
            i = z - start
            return (period, i + 1 if i < 2 else i + 11)
    for start, period in ((19, 4), (37, 5)):
        if start <= z < start + 18:
            return (period, z - start + 1)
    for start, period in ((55, 6), (87, 7)):
        i = z - start
        if not 0 <= i < 32:
            continue
        if i < 2: return (period, i + 1)
        if 2 <= i < 17: return (period + 3, i + 1)       # La-Lu / Ac-Lr on f rows, cols 3..17
        return (period, i - 13)
    raise ValueError(z)


def with_overlay(html, overlay):
    html = html.replace('font-size: 13px; overflow: hidden">', 'font-size: 13px; overflow: hidden; position: relative">', 1)
    i = html.rindex('\n</div>\n</x-dc>')
    return html[:i] + overlay + html[i:]


def modal(inner, w, h, top=None):
    t = f"top: {top}px; " if top is not None else f"top: {(900 - h) // 2}px; "
    return (f'<div style="position: absolute; inset: 0; background: #07080AB3"></div>'
            f'<div role="dialog" aria-modal="true" style="position: absolute; left: {(1440 - w) // 2}px; {t}width: {w}px; height: {h}px; box-sizing: border-box; display: flex; flex-direction: column; '
            f'background: {BG1}; border: 1px solid {LINE}; border-radius: 12px; box-shadow: 0 24px 64px #000000A0; overflow: hidden">{inner}</div>')


# ------------------------------------------------------------------ Command palette
def palette():
    base = screen_main.build()
    q = (f'<div style="display: flex; align-items: center; gap: 10px; height: 56px; padding: 0 18px; border-bottom: 1px solid {LINE}">{icon("search", 18, ACC)}'
         f'<label for="cp" style="position: absolute; width: 1px; height: 1px; overflow: hidden">Command</label>'
         f'<input id="cp" value="tacti" style="flex-grow: 1; background: transparent; border: 0; outline: none; font-size: 16px; color: {TEXT}">'
         f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}; border: 1px solid {LINE}; border-radius: 4px; padding: 2px 6px">esc</span></div>')

    def item(ic, title, cmd, keys, on=False):
        k = "".join(f'<kbd style="font-family: {MONO}; font-size: 11px; color: {MUTED}; border: 1px solid {LINE}; border-bottom-width: 2px; border-radius: 4px; padding: 1px 6px">{x}</kbd>' for x in keys)
        return (f'<button style="display: flex; align-items: center; gap: 12px; height: 44px; padding: 0 14px; text-align: left; background: {BG3 if on else "transparent"}; border: 0; border-radius: 7px; cursor: pointer">'
                f'{icon(ic, 17, ACC if on else MUTED)}<span style="flex-grow: 1; display: flex; flex-direction: column; gap: 1px"><span style="font-size: 13.5px; color: {TEXT}">{title}</span>'
                f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{cmd}</span></span><span style="display: flex; gap: 3px">{k}</span></button>')

    def grp(t):
        return f'<div style="padding: 10px 14px 4px; font-size: 10.5px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {DIM}">{t}</div>'
    body = (grp("Commands")
            + item("tag", "Set tacticity of selected chain…", "stereo.tacticity.set", ["⌥", "T"], True)
            + item("check", "Check tacticity along chain", "stereo.tacticity.check", [])
            + item("cursor", "Select by tacticity (m / r dyads)", "select.by_tacticity", [])
            + item("grow", "Polymer builder: tacticity", "builder.polymer.open --section stereo", [])
            + grp("Documents &amp; help")
            + item("cube", "PS_isotactic_DP40.caps", "document.open", [])
            + item("file", "Help: tacticity and dyad statistics", "help.open tacticity", ["F1"]))
    foot = (f'<div style="margin-top: auto; display: flex; align-items: center; gap: 16px; height: 40px; padding: 0 16px; border-top: 1px solid {LINE}; font-size: 11.5px; color: {DIM}">'
            f'<span>↑ ↓ move</span><span>↵ run</span><span>⇥ fill arguments</span><span style="margin-left: auto; font-family: {MONO}">same command in Python: doc.stereo.tacticity.set(…)</span></div>')
    inner = q + f'<div style="padding: 6px; display: flex; flex-direction: column">{body}</div>' + foot
    return with_overlay(base, modal(inner, 680, 470, top=120)).replace("<title>CAPS Studio — builder</title>", "<title>CAPS Studio — command palette</title>")


# ------------------------------------------------------------------ Periodic table
def periodic():
    base = screen_main.build()
    cell = 38
    gap = 3
    common = {"H", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "Na", "K", "Ca", "Li", "I"}
    cells = []
    for z, s in enumerate(SYMS, start=1):
        r, c = pt_pos(z)
        x = (c - 1) * (cell + gap)
        y = (r - 1) * (cell + gap) + (10 if r >= 8 else 0)
        on = s == "N"
        col_ = ELEM.get(s, (None,))[0]
        bg = ACC if on else (BG3 if s in common else BG0)
        fg = ACC_INK if on else (TEXT if s in common else DIM)
        bd = ACC if on else (LINE if s in common else BG2)
        stripe = f'<span style="position: absolute; left: 0; top: 0; bottom: 0; width: 3px; background: {col_}"></span>' if col_ and not on else ""
        cells.append(f'<button aria-label="{s}, atomic number {z}" style="position: absolute; left: {x}px; top: {y}px; width: {cell}px; height: {cell}px; padding: 0; border-radius: 5px; border: 1px solid {bd}; background: {bg}; color: {fg}; cursor: pointer; overflow: hidden">'
                     f'{stripe}<span style="position: absolute; left: 5px; top: 2px; font-family: {MONO}; font-size: 8.5px; opacity: 0.8">{z}</span>'
                     f'<span style="position: absolute; left: 0; right: 0; top: 13px; text-align: center; font-size: 13.5px; font-weight: 600">{s}</span></button>')
    for r, lbl in ((6, "57–71"), (7, "89–103")):
        x = 2 * (cell + gap)
        y = (r - 1) * (cell + gap)
        cells.append(f'<span style="position: absolute; left: {x}px; top: {y}px; width: {cell}px; height: {cell}px; display: flex; align-items: center; justify-content: center; font-family: {MONO}; font-size: 8.5px; color: {DIM}; border: 1px dashed {LINE}; border-radius: 5px; box-sizing: border-box">{lbl}</span>')
    W = 18 * (cell + gap)
    H = 10 * (cell + gap) + 10
    grid = f'<div style="position: relative; width: {W}px; height: {H}px">{"".join(cells)}</div>'
    details = col(
        row(f'<span style="width: 64px; height: 64px; border-radius: 10px; background: {ACC}; color: {ACC_INK}; display: flex; align-items: center; justify-content: center; font-size: 30px; font-weight: 700">N</span>',
            col(f'<span style="font-size: 18px; font-weight: 600">Nitrogen</span>', f'<span style="font-family: {MONO}; font-size: 12px; color: {MUTED}">Z = 7 · 14.007 u</span>', gap=4), gap=14),
        section("Isotope", col(seg(["natural", "¹⁴N", "¹⁵N"], "natural", True), f'<span style="font-size: 11.5px; color: {DIM}">¹⁴N 99.636 % · ¹⁵N 0.364 %</span>', gap=8), pad=0, gap=8),
        section("Formal charge", seg(["−1", "0", "+1", "+2"], "0", True), pad=0, gap=8),
        section("Geometry on placement", seg(["sp³ · 3 bonds", "sp² · 2", "sp · 1", "N⁺ · 4"], "sp³ · 3 bonds", cols=2), pad=0, gap=8),
        kv("Covalent radius", "0.71 Å"), kv("Pauling electronegativity", "3.04"), cite("Radii: Cordero et al., <i>Dalton Trans.</i> 2832 (2008)"), gap=14)
    head = (f'<div style="height: 52px; display: flex; align-items: center; gap: 12px; padding: 0 20px; border-bottom: 1px solid {LINE}"><h2 style="margin: 0; font-size: 15px; font-weight: 600">Element</h2>'
            f'{chip("stripe = common in polymers", MUTED)}<div style="flex-grow: 1"></div><span style="font-size: 12px; color: {DIM}">type a symbol to jump · ⇧E opens this</span>{tbtn("close", "Close")}</div>')
    inner = (head + f'<div style="flex-grow: 1; display: flex; gap: 24px; padding: 20px">{grid}'
             f'<div style="flex-grow: 1; min-width: 0; padding-left: 20px; border-left: 1px solid {LINE}">{details}</div></div>'
             + f'<div style="display: flex; justify-content: flex-end; gap: 8px; padding: 12px 20px; border-top: 1px solid {LINE}">{btn("Replace selected atom", ic="rotate")}{btn("Place N", True, "atom")}</div>')
    return with_overlay(base, modal(inner, 1120, 640)).replace("<title>CAPS Studio — builder</title>", "<title>CAPS Studio — element picker</title>")


# ------------------------------------------------------------------ Selection & stereo
def selection():
    seed, units = 5, 7
    atoms, backbone, stereo = polystyrene(units, "atactic", seed=seed)
    rnd = random.Random(seed)
    sides = [rnd.choice([1, -1]) for _ in range(units)]
    dyads = ["m" if sides[i] == sides[i + 1] else "r" for i in range(units - 1)]
    triads = ["".join(dyads[i:i + 2]) for i in range(len(dyads) - 1)]
    tri_c = {k: sum(1 for t in triads if "".join(sorted(t)) == k) for k in ("mm", "mr", "rr")}
    ar = [i for i, a in enumerate(atoms) if a.get("ar")]
    W = 1440 - 72 - 400
    sc = Scene("s7", W, 900 - 44 - 48 - 26 - 150, yaw=-0.28, pitch=0.42, roll=0.05, persp=0.25, fog=0.6)
    sc.add_atoms(atoms)
    sc.halo = set(ar)
    sc._prep()
    ov = ""
    for n, k in enumerate(stereo):
        x, y, _, _ = sc.screen(k)
        ov += label_pill(x - 8, y + 24, f"{n + 1}", MUTED, fs=10.5)
    svg = sc.svg(overlay=ov)
    dy_cells = []
    for i, d in enumerate(dyads):
        c = ACC if d == "m" else SEL
        dy_cells.append(f'<span style="display: flex; flex-direction: column; align-items: center; gap: 4px"><span style="width: 44px; height: 30px; border-radius: 5px; background: {c}; color: {ACC_INK}; display: flex; align-items: center; justify-content: center; font-family: {MONO}; font-weight: 600">{d}</span>'
                        f'<span style="font-family: {MONO}; font-size: 10.5px; color: {DIM}">{i + 1}–{i + 2}</span></span>')
    strip = (f'<div style="height: 150px; box-sizing: border-box; padding: 14px 18px; background: {BG1}; border-top: 1px solid {LINE}; display: flex; flex-direction: column; gap: 10px">'
             + row(f'<h3 style="margin: 0; font-size: 12px; font-weight: 600">Tacticity along chain · dyads between stereocentres</h3>', spacer(), chip(f'{dot(ACC)} m (meso)'), chip(f'{dot(SEL)} r (racemo)'), gap=8)
             + f'<div style="display: flex; gap: 8px">{"".join(dy_cells)}</div>'
             + row(chip(f"m = {dyads.count('m')} · r = {dyads.count('r')}", TEXT, BG3, True), chip(f"triads mm {tri_c['mm']} · mr {tri_c['mr']} · rr {tri_c['rr']}", TEXT, BG3, True), chip("computed from the drawn chain", MUTED), gap=6)
             + '</div>')
    center = f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}"><div style="position: relative; flex-grow: 1">{svg}{hud(chip(f"{len(ar)} atoms match", SEL, "#16191Ccc"), chip("SMARTS c1ccccc1", MUTED, "#16191Ccc"))}</div>{strip}</div>'
    modes = [("SMARTS", True), ("Element", False), ("Atom type", False), ("Charge range", False), ("Within distance", False), ("Grow along bonds", False)]
    mg = "".join(f'<button style="height: 30px; border-radius: 5px; border: 1px solid {ACC if on else LINE}; background: {BG3 if on else BG0}; font-size: 12px; color: {TEXT if on else MUTED}; cursor: pointer">{m}</button>' for m, on in modes)
    sets = [("ring-7", "6 atoms", SEL), ("backbone", f"{len(backbone)} atoms", ACC), ("chain-ends", "2 atoms", MUTED)]
    sr = "".join(f'<div style="display: flex; align-items: center; gap: 10px; height: 32px; padding: 0 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px">{dot(c)}<span style="flex-grow: 1; font-size: 12.5px">{n}</span><span style="font-family: {MONO}; font-size: 11.5px; color: {DIM}">{k}</span>{icon("dots", 14, DIM)}</div>' for n, k, c in sets)
    right = (panel_head("Selection", chip(f"{len(ar)} selected", SEL))
             + section("Select by", col(f'<div style="display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 4px">{mg}</div>',
                                         field("SMARTS pattern", "c1ccccc1"), row(btn("Replace", small=True), btn("Add", small=True), btn("Subtract", small=True), btn("Invert", small=True), gap=6), gap=10))
             + section("Named sets", col(sr, btn("Save selection as set", ic="save", small=True), gap=6))
             + section("Stereo", col(kv("Stereocentres (CH)", f"{len(stereo)}"), kv("Tacticity", "atactic"),
                                      f'<div style="font-size: 12px; color: {MUTED}; line-height: 1.5">Along a long chain these centres are pseudo-asymmetric, so tacticity is reported as m/r dyads rather than R/S labels.</div>',
                                      row(btn("Make isotactic", small=True), btn("Make syndiotactic", small=True), gap=6), gap=8))
             + footer(btn("Invert centre", ic="mirror"), btn("Apply", True, "check")))
    body = (topbar(["PS_atactic_DP40.caps"], "PS_atactic_DP40.caps") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">'
            + bar(tbtn("cursor", "Select", True), tbtn("lasso", "Lasso"), tbtn("search", "Select by pattern", True, "Select by…"), sep(), tbtn("tag", "Stereo", False, "Stereo"), spacer(), tbtn("eye", "Style", False, "Ball &amp; stick"))
            + f'<div style="flex-grow: 1; display: flex; min-height: 0">{center}<aside aria-label="Selection" style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{right}</aside></div></div></div>'
            + statusbar(f"<span>{len(atoms)} atoms · {len(ar)} selected</span>", "<span>Selections are commands: select.smarts, select.save</span>"))
    return page("CAPS Studio — selection & stereo", body)


# ------------------------------------------------------------------ H-bonds, contacts, clashes
def contacts():
    rnd = random.Random(12)
    L = 15.0
    Os = []
    tries = 0
    while len(Os) < 70 and tries < 60000:
        tries += 1
        o = (rnd.uniform(0.8, L - 0.8), rnd.uniform(0.8, L - 0.8), rnd.uniform(0.8, L - 0.8))
        if all(dist(o, q) > 2.65 for q in Os):
            Os.append(o)
    # orient: first H toward the nearest neighbour O (typical donor geometry), second H at 104.52 degrees
    mols = []
    for i, o in enumerate(Os):
        nb = sorted((dist(o, q), j) for j, q in enumerate(Os) if j != i)
        d1 = norm(sub(Os[nb[0][1]], o))
        tmp = norm(sub(Os[nb[1][1]], o))
        perp = norm(sub(tmp, mul(d1, dot3(tmp, d1))))
        a = math.radians(104.52)
        h1 = add(o, mul(d1, 0.9572))
        h2 = add(o, mul(norm(add(mul(d1, math.cos(a)), mul(perp, math.sin(a)))), 0.9572))
        mols.append((o, h1, h2))
    # deliberately overlapping water (the clash to fix)
    co = add(Os[5], (1.25, 0.9, 0.3))
    mols.append((co, add(co, (0.9572, 0, 0)), add(co, (0.9572 * math.cos(math.radians(104.52)), 0.9572 * math.sin(math.radians(104.52)), 0))))
    sc = Scene("h", 1440 - 72 - 400, 900 - 44 - 48 - 26, yaw=0.55, pitch=0.4, persp=0.28, fog=0.55)
    for o, h1, h2 in mols:
        sc.add_atoms([{"e": "O", "p": o}, {"e": "H", "p": h1}, {"e": "H", "p": h2}], [(0, 1), (0, 2)])
    # H-bond detection: R(O...O) < 3.5 A and angle(H-Od...Oa) < 30 deg (Luzar & Chandler 1996)
    hb = []
    for i, (o, h1, h2) in enumerate(mols):
        for j, (oa, _, _) in enumerate(mols):
            if i == j:
                continue
            r = dist(o, oa)
            if r < 3.5:
                for h in (h1, h2):
                    if angle(h, o, oa) < 30.0:
                        hb.append((h, oa, r))
    clashes = []
    for i in range(len(mols)):
        for j in range(i + 1, len(mols)):
            r = dist(mols[i][0], mols[j][0])
            if r < 2.4:
                clashes.append((mols[i][0], mols[j][0], r))
    for h, oa, r in hb:
        sc.add_line(h, oa, SEL, 1.4, "3 3", 0.9)
    for a, b, r in clashes:
        sc.add_line(a, b, ERR, 2.4, None, 1.0)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, "4 4", 0.6)
    sc._prep()
    ov = ""
    for a, b, r in clashes:
        pa, pb = sc.P(a), sc.P(b)
        ov += (f'<circle cx="{(pa[0] + pb[0]) / 2:.1f}" cy="{(pa[1] + pb[1]) / 2:.1f}" r="26" fill="none" stroke="{ERR}" stroke-width="2" stroke-dasharray="5 4"></circle>'
               + label_pill((pa[0] + pb[0]) / 2 + 30, (pa[1] + pb[1]) / 2 - 20, f"clash O···O {r:.2f} Å", ERR, border=ERR))
    svg = sc.svg(overlay=ov)
    nh = len(hb)
    rows = [[f'{icon("xcircle", 15, ERR)}', f"Overlaps {len(clashes)} neighbours", f"O···O ≥ {min(c[2] for c in clashes):.2f} Å", f'<button style="height: 24px; padding: 0 8px; background: transparent; border: 1px solid {LINE}; border-radius: 5px; font-size: 11.5px; color: {ACC}; cursor: pointer">Relax</button>'],
            [f'{icon("alert", 15, WARN)}', "Box edge cuts 4 molecules", "wrap needed", f'<button style="height: 24px; padding: 0 8px; background: transparent; border: 1px solid {LINE}; border-radius: 5px; font-size: 11.5px; color: {ACC}; cursor: pointer">Wrap</button>'],
            [f'{icon("check", 15, OK)}', "Valences", "all O 2-bonded", ""],
            [f'{icon("check", 15, OK)}', "Net charge", "0.000 e", ""]]
    vt = table(["", "Issue", "Detail", ""], rows, ["7%", "45%", "32%", "16%"], mono_cols=(2,), fs=11.5)
    right = (panel_head("Interactions &amp; checks")
             + section("Hydrogen bonds", col(row(chip(f"{nh} found", SEL, "#12303A", True), chip(f"{len(mols)} waters", MUTED), gap=6),
                                              row(field("D···A max", "3.5", "Å"), field("H–D···A max", "30", "°"), gap=10),
                                              cite("Geometric criterion: Luzar &amp; Chandler, <i>Nature</i> 379, 55 (1996)"), gap=10))
             + section("Contacts", col(toggle("Show close contacts &lt; sum of vdW − 0.4 Å", True), toggle("Show clashes &lt; 0.75 × sum of vdW", True), toggle("Update while editing and during playback", True), gap=9))
             + section("Validation", vt, chip(f"{len(clashes)} errors", ERR, "#2A1414"))
             + footer(btn("Export list", ic="download"), btn("Fix all safe issues", True, "wand")))
    legend = hud(chip(f'<span style="width: 16px; border-top: 2px dashed {SEL}; display: inline-block"></span> H-bond', TEXT, "#16191Ccc"), chip(f'<span style="width: 16px; border-top: 2px solid {ERR}; display: inline-block"></span> clash', TEXT, "#16191Ccc"))
    vp = f'<main aria-label="3D view" style="position: relative; flex-grow: 1; min-width: 0; background: {BG0}">{svg}{legend}</main>'
    body = (topbar(["water_box.caps"], "water_box.caps") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">'
            + bar(tbtn("cursor", "Select", True), tbtn("ruler", "Measure"), sep(), tbtn("link", "H-bonds", True, "H-bonds"), tbtn("alert", "Clashes", True, "Clashes"), spacer(), tbtn("eye", "Style", False, "Ball &amp; stick"))
            + f'<div style="flex-grow: 1; display: flex; min-height: 0">{vp}<aside aria-label="Interactions" style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{right}</aside></div></div></div>'
            + statusbar(f"<span>{len(mols) * 3} atoms · {nh} H-bonds · {len(clashes)} clashes</span>", "<span>Counts computed from the drawn geometry</span>"))
    return page("CAPS Studio — interactions & checks", body)


# ------------------------------------------------------------------ Light theme split view
LIGHT = {"#0F1113": "#F4F2EE", "#16191C": "#FFFFFF", "#1E2226": "#F1EFEB", "#262B30": "#E7E4DF", "#2C3137": "#D9D5CF",
         "#E8E6E1": "#1E2226", "#A5ABB1": "#555B62", "#868D94": "#6B7178", "#F0A83C": "#A35C00", "#1A1206": "#FFFFFF",
         "#6CC4D8": "#0F7C93", "#7CC784": "#2E7D32", "#FF7B72": "#C62828", "#3A2C14": "#FBEBD2", "#16261A": "#E3F2E4",
         "#2A1414": "#FDE8E8", "#16191CCC": "#FFFFFFDD", "#16191CE6": "#FFFFFFEE", "#12303A": "#DDF1F5"}


def lightify(html):
    return re.sub(r"#[0-9A-Fa-f]{8}|#[0-9A-Fa-f]{6}", lambda m: LIGHT.get(m.group(0).upper(), m.group(0)), html)


def light_split():
    LBG = "#F4F2EE"
    Wv = (1440 - 72) // 2 - 1
    Hv = 900 - 44 - 48 - 34 - 222 - 26

    def scene(sid, tact, seed):
        atoms, backbone, stereo = polystyrene(7, tact, seed=seed)
        sc = Scene(sid, Wv, Hv, yaw=-0.28, pitch=0.42, roll=0.05, persp=0.25, fog=0.35, bg=LBG)
        sc.add_atoms(atoms)
        return sc.svg(), atoms, stereo
    s1, a1, st1 = scene("la", "atactic", 5)
    s2, a2, st2 = scene("lb", "isotactic", 1)

    def pane(title, tag, placeholder, active):
        head = (f'<div style="height: 34px; display: flex; align-items: center; gap: 8px; padding: 0 12px; background: {BG1}; border-bottom: 2px solid {ACC if active else LINE}">'
                f'{icon("cube", 14, ACC if active else DIM)}<span style="font-size: 12.5px; font-weight: 500">{title}</span>{chip(tag)}<div style="flex-grow: 1"></div>{tbtn("dots", "Pane menu")}</div>')
        return f'<div style="width: {Wv}px; display: flex; flex-direction: column; background: {BG0}">{head}<div style="position: relative; height: {Hv}px">{placeholder}</div></div>'
    rows = [["Stereocentres", "7", "7"], ["Tacticity", "atactic", "isotactic"], ["Atoms", str(len(a1)), str(len(a2))],
            ["Rg (drawn fragment)", f"{rg(a1):.2f} Å", f"{rg(a2):.2f} Å"]]
    cmp_ = table(["Property", "PS_atactic_DP40", "PS_isotactic_DP40"], rows, ["30%", "35%", "35%"], mono_cols=(1, 2), align_right=(1, 2), rowh=27)
    compare = (f'<div style="height: 222px; flex-shrink: 0; display: flex; gap: 20px; padding: 14px 20px; box-sizing: border-box; background: {BG1}; border-top: 1px solid {LINE}">'
               f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 8px"><h3 style="margin: 0; font-size: 12px; font-weight: 600">Compare documents</h3>{cmp_}</div>'
               f'<div style="width: 360px; display: flex; flex-direction: column; gap: 10px"><h3 style="margin: 0; font-size: 12px; font-weight: 600">Linked views</h3>'
               + toggle("Sync camera between panes", True) + toggle("Sync selection by atom index", False) + toggle("Sync trajectory frame", True) + '</div></div>')
    chrome = (topbar(["PS_atactic_DP40.caps", "PS_isotactic_DP40.caps"], "PS_atactic_DP40.caps") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
              + '<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">'
              + bar(tbtn("cursor", "Select", True), tbtn("rotate", "Rotate"), tbtn("ruler", "Measure"), sep(), tbtn("layers", "Split", True, "Split vertical"), tbtn("link", "Sync", True, "Camera synced"), spacer(), tbtn("eye", "Theme", False, "Theme: Paper (light)"))
              + f'<div style="display: flex; gap: 2px; background: {LINE}">{pane("PS_atactic_DP40.caps", "atactic", "%%S1%%", True)}{pane("PS_isotactic_DP40.caps", "isotactic", "%%S2%%", False)}</div>'
              + compare + '</div></div>'
              + statusbar("<span>2 documents · camera synced</span>", "<span>Theme: Paper (light) · follows system</span>"))
    html = lightify(page("CAPS Studio — light theme, split view", chrome))
    return html.replace("%%S1%%", s1).replace("%%S2%%", s2)


def rg(atoms):
    heavy = [a["p"] for a in atoms if a["e"] != "H"]
    c = tuple(sum(p[k] for p in heavy) / len(heavy) for k in range(3))
    return math.sqrt(sum(dist(p, c) ** 2 for p in heavy) / len(heavy))


# ------------------------------------------------------------------ Remote compute
def remote():
    navs = [("General", "gear"), ("Appearance", "eye"), ("Input &amp; shortcuts", "terminal"), ("Rendering", "cube"), ("Force fields", "tag"),
            ("Compute &amp; remote", "server"), ("Accessibility", "check"), ("Python &amp; scripting", "copy")]
    act = "Compute &amp; remote"
    nv = "".join(f'<a href="{ {"Appearance": "Settings.dc.html", "Compute &amp; remote": "RemoteCompute.dc.html"}.get(n, "#") }" aria-current="{"page" if n == act else "false"}" style="display: flex; align-items: center; gap: 10px; height: 36px; padding: 0 12px; border-radius: 6px; text-decoration: none; font-size: 13px; background: {BG3 if n == act else "transparent"}; color: {TEXT if n == act else MUTED}">{icon(i, 16, ACC if n == act else DIM)}{n}</a>' for n, i in navs)
    hosts = [("This Mac", "local · 10 threads · CPU path", "ready", OK), ("hpc-login2", "caps-daemon 0.1.0 over SSH tunnel · SLURM", "connected", OK), ("lab-gpu-01", "caps-daemon over SSH · no scheduler", "unreachable", ERR)]
    hr = "".join(f'<div style="display: flex; align-items: center; gap: 12px; padding: 12px 14px; background: {BG3 if n == "hpc-login2" else BG0}; border: 1px solid {ACC if n == "hpc-login2" else LINE}; border-radius: 8px">{icon("server", 18, MUTED)}'
                 f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 2px"><span style="font-size: 13px; font-weight: 500">{n}</span><span style="font-size: 11.5px; color: {MUTED}">{d}</span></div>'
                 f'<span style="display: flex; align-items: center; gap: 6px; font-family: {MONO}; font-size: 11.5px; color: {c}">{dot(c)}{s}</span></div>' for n, d, s, c in hosts)
    sbatch = ["#!/bin/bash", "#SBATCH --job-name=caps-{job}", "#SBATCH --partition=gpu", "#SBATCH --gres=gpu:1", "#SBATCH --cpus-per-task=8",
              "#SBATCH --time=24:00:00", "#SBATCH --output=caps-%j.log", "module load cuda/12.4", "caps run {recipe} --checkpoint-every 50ps"]
    sb = "".join(f'<div style="color: {DIM if l.startswith("#SBATCH") else TEXT}">{esc(l)}</div>' for l in sbatch)
    content = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 18px; padding: 28px 36px">'
               f'<div style="display: flex; flex-direction: column; gap: 4px"><h1 style="margin: 0; font-size: 22px; font-weight: 600">Compute &amp; remote</h1>'
               f'<p style="margin: 0; font-size: 13px; color: {MUTED}">Jobs run here or on a remote host through the headless caps-daemon. Files move with the job; results come back with their provenance manifest.</p></div>'
               + col(row(f'<h2 style="margin: 0; font-size: 13px; font-weight: 600">Hosts</h2>', spacer(), btn("Add host", ic="plus", small=True), gap=8), hr, gap=10)
               + col(f'<h2 style="margin: 0; font-size: 13px; font-weight: 600">hpc-login2</h2>',
                     row(field("Hostname", "hpc-login2.example.edu", mono=True), field("User", "username"), field("Port", "22", w=90), gap=10),
                     row(select("Authentication", "SSH agent (no passwords stored)"), select("Scheduler", "SLURM"), select("Default partition", "gpu"), gap=10),
                     row(field("Remote work directory", "/scratch/$USER/caps"), select("Daemon", "Start on demand · caps-daemon 0.1.0"), gap=10),
                     row(btn("Test connection", ic="link", small=True), chip(f'{dot(OK)} daemon reachable · version matches', OK, "#16261A"), gap=10), gap=12)
               + '</div>')
    side = (f'<div style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px; padding: 28px 28px 28px 0">'
            + card("Job template · SLURM", f'<div style="padding: 12px 14px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px; font-family: {MONO}; font-size: 11.5px; line-height: 1.75">{sb}</div>', chip("editable"))
            + card("Transfers", col(kv("Upload", "inputs + recipe, hashed"), kv("Download", "outputs + manifest"), kv("Large trajectories", "stay remote, stream on view", False), gap=8))
            + '</div>')
    body = (topbar(["Settings"], "Settings") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("")
            + f'<nav aria-label="Settings sections" style="width: 240px; flex-shrink: 0; display: flex; flex-direction: column; gap: 2px; padding: 20px 12px; background: {BG1}; border-right: 1px solid {LINE}">{nv}</nav>'
            + f'<div style="flex-grow: 1; display: flex; min-width: 0">{content}{side}</div></div>'
            + statusbar("<span>2 of 3 hosts reachable</span>", "<span>Credentials stay in your SSH agent</span>"))
    return page("CAPS — compute & remote", body)
