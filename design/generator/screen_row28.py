"""Row 28 — the pipeline in order: a step strip over every page, the force field chosen in Grow and assigned when growing
finishes, a Force field step that says which force fields fit the structure (and what is missing from the others), and
an Export center reachable from any step. Real data: a polychloroprene (CR) melt grown by CAPS (10 chains x 30 units,
3020 atoms), the coverage CAPS computed for it over the library, and the LAMMPS / GROMACS files CAPS wrote for it."""
import json, math, os, re
from lib import *
from screen_app import card, progress

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "boarddata")
H3 = f"margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}"

# ------------------------------------------------------------------ the structure
def read_data(path):
    lines = open(path).read().split("\n")
    L = None
    atoms, bonds, sec = [], [], None
    for l in lines:
        s = l.split("#")[0].strip()
        if s.endswith("xlo xhi"):
            lo, hi = map(float, s.split()[:2]); L = hi - lo
        if re.match(r"^[A-Z][A-Za-z ]+$", s):
            sec = s; continue
        w = s.split()
        if sec == "Atoms" and len(w) >= 7:
            ix = [int(v) for v in w[7:10]] if len(w) >= 10 else [0, 0, 0]
            atoms.append((int(w[0]), int(w[1]), int(w[2]), float(w[3]), [float(w[4 + k]) + ix[k] * L for k in range(3)]))
        if sec == "Bonds" and len(w) >= 4:
            bonds.append((int(w[2]), int(w[3])))
    atoms.sort()
    return L, atoms, bonds


L, ATOMS, BONDS = read_data(os.path.join(DATA, "cr.data"))
ELT = {1: "C", 2: "H", 3: "Cl"}
NB = {}
for i, j in BONDS:
    NB.setdefault(i, []).append(j); NB.setdefault(j, []).append(i)
# each chain shifted so its centre lies in the cell
POS = {}
for m in sorted({a[1] for a in ATOMS}):
    ids = [a for a in ATOMS if a[1] == m]
    c = [sum(a[4][k] for a in ids) / len(ids) for k in range(3)]
    sh = [-L * math.floor(c[k] / L) for k in range(3)]
    for a in ids:
        POS[a[0]] = tuple(a[4][k] + sh[k] for k in range(3))
EL = {a[0]: ELT[a[2]] for a in ATOMS}
MOL = {a[0]: a[1] for a in ATOMS}
CHAIN = ["#E8A33D", "#5FB3C8", "#8CC084", "#D98A6C", "#A993D6", "#D3B25C", "#6FA3D8", "#C98FB8", "#9CB27A", "#D07C7C"]


def backbone(m):
    """The chain's carbons in bond order (the backbone of a 1,4 unit: all four C)."""
    cs = [a[0] for a in ATOMS if a[1] == m and EL[a[0]] == "C"]
    cset = set(cs)
    ends = [c for c in cs if sum(1 for j in NB.get(c, []) if j in cset) == 1]
    path, prev, cur = [], None, ends[0]
    while cur is not None:
        path.append(cur)
        nxt = [j for j in NB[cur] if j in cset and j != prev]
        prev, cur = cur, (nxt[0] if nxt else None)
    return path


def cell_scene(sid, w, h, mode, bg=BG0, yaw=0.62, pitch=0.42):
    sc = Scene(sid, w, h, yaw=yaw, pitch=pitch, persp=0.0, fog=0.45, atom_k=1.0, center=(L / 2, L / 2, L / 2), bg=bg, bond_w=0.18)
    if mode == "chains":
        pass   # drawn by chains_overlay: one continuous line per chain, with its Cl atoms
    else:
        # heavy atoms; the atoms a force field could not treat drawn in its colour, the rest dimmed
        heavy = [a[0] for a in ATOMS if EL[a[0]] != "H"]
        idx = {a: k for k, a in enumerate(heavy)}
        hot = set(mode)
        recs = []
        for a in heavy:
            c = ("#57B26A" if EL[a] == "Cl" else ACC) if a in hot else mix(ELEM[EL[a]][0], bg, 0.55)
            recs.append(dict(p=POS[a], e=EL[a], c=c, r=0.34 if EL[a] == "C" else 0.5))
        bl = [(idx[i], idx[j]) for i, j in BONDS if i in idx and j in idx and dist(POS[i], POS[j]) < 3]
        sc.add_atoms(recs, bl)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), DIM, 1.0, None, 0.8)
    return sc


def chains_overlay(sc, bg):
    """Each chain as one polyline (dark rim, colour, highlight), its Cl atoms on it; chains painted back to front."""
    sc._prep()
    groups = []
    for m in range(1, 11):
        bb = backbone(m)
        pp = [sc.P(POS[i]) for i in bb]
        cls = [sc.P(POS[j]) for i in bb for j in NB[i] if EL[j] == "Cl"]
        z = sum(p[2] for p in pp) / len(pp)
        colr = CHAIN[(m - 1) % len(CHAIN)]
        w = 0.5 * sc.scale
        pts = " ".join(f"{p[0]:.1f},{p[1]:.1f}" for p in pp)
        s = (f'<polyline points="{pts}" fill="none" stroke="#070809" stroke-width="{w + 2.2:.1f}" stroke-linejoin="round" stroke-linecap="round"></polyline>'
             f'<polyline points="{pts}" fill="none" stroke="{colr}" stroke-width="{w:.1f}" stroke-linejoin="round" stroke-linecap="round"></polyline>'
             f'<polyline points="{pts}" fill="none" stroke="{mix(colr, "#FFFFFF", 0.4)}" stroke-width="{w * 0.28:.1f}" stroke-opacity="0.6" stroke-linejoin="round" stroke-linecap="round" transform="translate(0 {-w * 0.2:.1f})"></polyline>')
        r = 0.55 * sc.scale
        s += "".join(f'<circle cx="{c[0]:.1f}" cy="{c[1]:.1f}" r="{r:.1f}" fill="#57B26A" stroke="#070809" stroke-width="0.9"></circle>'
                     f'<circle cx="{c[0] - r * 0.3:.1f}" cy="{c[1] - r * 0.3:.1f}" r="{r * 0.35:.1f}" fill="#FFFFFF" fill-opacity="0.35"></circle>' for c in cls)
        groups.append((z, s))
    groups.sort(key=lambda g: g[0])
    return "".join(s for _, s in groups)


COV = json.load(open(os.path.join(DATA, "cov_cr.json")))["forcefields"]
TYPED = [f for f in COV if f["status"] != "no typing rules"]
NO_RULES = [f for f in COV if f["status"] == "no typing rules"]


TOTAL = f"{sum(os.path.getsize(os.path.join(DATA, 'cr_out', f)) for f in ('system.data', 'system.in', 'system.top', 'system.itp', 'system.gro', 'system.mdp')) / 1024:.0f} KB"


def fsize(p):
    n = os.path.getsize(os.path.join(DATA, "cr_out", p))
    return f"{n / 1024:.0f} KB" if n >= 10240 else f"{n / 1024:.1f} KB"


# ------------------------------------------------------------------ chrome
MODULES2 = [("Studio", "hex", "Main.dc.html"), ("Build", "flask", "PolymerBuilder.dc.html"), ("Grow", "grow", "PipelineGrow.dc.html"),
            ("Force field", "tag", "ForceFieldStep.dc.html"), ("Pack", "pack", "Pack.dc.html"), ("Relax", "relax", "Relax.dc.html"),
            ("Equilibrate", "equil", "Convergence.dc.html"), ("Dynamics", "dyn", "Dynamics.dc.html"), ("React", "react", "React.dc.html"),
            ("Analyze", "chart", "Analyze.dc.html"), ("Export", "download", "ExportCenter.dc.html"), ("Jobs", "jobs", "Jobs.dc.html"),
            ("Bench", "bench", "Bench.dc.html")]


def rail2(active):
    items = []
    for name, ic, href in MODULES2:
        on = name == active
        st = f"background: {BG3}; color: {TEXT}" if on else f"background: transparent; color: {DIM}"
        bar = f'<span style="position: absolute; left: 0; top: 10px; width: 3px; height: 28px; border-radius: 0 2px 2px 0; background: {ACC}"></span>' if on else ""
        items.append(f'<a href="{href}" aria-current="{"page" if on else "false"}" style="position: relative; width: 64px; height: 50px; display: flex; flex-direction: column; align-items: center; justify-content: center; gap: 3px; border-radius: 8px; text-decoration: none; font-size: 10px; text-align: center; line-height: 1.1; {st}">'
                     f'{bar}{icon(ic, 18, ACC if on else MUTED)}<span>{name}</span></a>')
        if name in ("Studio", "Analyze"):
            items.append(f'<div style="width: 36px; height: 1px; background: {LINE}; margin: 3px 0"></div>')
    return (f'<nav aria-label="Modules" style="width: 72px; flex-shrink: 0; display: flex; flex-direction: column; align-items: center; gap: 1px; padding-top: 6px; background: {BG1}; border-right: 1px solid {LINE}">'
            + "".join(items) +
            f'<div style="flex-grow: 1"></div><a href="Settings.dc.html" aria-label="Settings" style="width: 60px; height: 40px; margin-bottom: 6px; display: flex; align-items: center; justify-content: center; border-radius: 8px">{icon("gear", 18, MUTED)}</a></nav>')


# steps of the strip: (name, board, state, detail); state: done | current | auto | next | optional | todo
def strip(steps, next_label, next_href):
    cells = []
    for k, (name, href, state, detail) in enumerate(steps):
        ring = {"done": OK, "auto": OK, "current": ACC, "next": TEXT, "optional": DIM, "todo": DIM}[state]
        if state in ("done", "auto"):
            mark = f'<span style="width: 22px; height: 22px; border-radius: 50%; background: {OK}; display: flex; align-items: center; justify-content: center; flex-shrink: 0">{icon("check", 13, BG0, 3)}</span>'
        elif state == "current":
            mark = f'<span style="width: 22px; height: 22px; border-radius: 50%; border: 2px solid {ACC}; box-sizing: border-box; display: flex; align-items: center; justify-content: center; flex-shrink: 0; font-family: {MONO}; font-size: 11px; font-weight: 600; color: {ACC}">{k + 1}</span>'
        else:
            dash = "dashed" if state == "optional" else "solid"
            mark = f'<span style="width: 22px; height: 22px; border-radius: 50%; border: 1.5px {dash} {ring}; box-sizing: border-box; display: flex; align-items: center; justify-content: center; flex-shrink: 0; font-family: {MONO}; font-size: 11px; color: {MUTED}">{k + 1}</span>'
        lab_col = TEXT if state in ("current", "next", "done", "auto") else MUTED
        bg = f"background: {BG3}; box-shadow: inset 0 -2px 0 {ACC};" if state == "current" else ""
        cells.append(f'<a href="{href}" style="display: flex; align-items: center; gap: 7px; padding: 6px 8px; border-radius: 7px; text-decoration: none; flex-shrink: {0 if state in ("current", "done", "auto") else 1}; min-width: 0; {bg}">{mark}'
                     f'<span style="display: flex; flex-direction: column; min-width: 0"><span style="font-size: 12.5px; font-weight: 600; color: {lab_col}; white-space: nowrap">{name}</span>'
                     f'<span style="font-size: 11px; color: {OK if state == "auto" else DIM}; white-space: nowrap; overflow: hidden; text-overflow: ellipsis">{detail}</span></span></a>')
        if k < len(steps) - 1:
            cells.append(f'<span style="width: 10px; height: 1px; flex-shrink: 0; background: {LINE}"></span>')
    right = (f'<div style="margin-left: auto; display: flex; align-items: center; gap: 8px; padding-left: 12px; border-left: 1px solid {LINE}">'
             + btn("Export now", False, "download", True, "ExportCenter.dc.html").replace("height: 30px;", "height: 30px; white-space: nowrap;") + '</div>')
    return (f'<div role="navigation" aria-label="Pipeline" style="height: 56px; flex-shrink: 0; display: flex; align-items: center; gap: 2px; padding: 0 12px; background: {BG1}; border-bottom: 1px solid {LINE}">'
            + "".join(cells) + right + '</div>')


def steps_for(current):
    base = [("Build", "PolymerBuilder.dc.html", "done", "CR unit"),
            ("Grow", "PipelineGrow.dc.html", "done", "3020 atoms"),
            ("Force field", "ForceFieldStep.dc.html", "auto", "GAFF · auto"),
            ("Pack", "Pack.dc.html", "optional", "optional"),
            ("Relax", "Relax.dc.html", "next", "next"),
            ("Equilibrate", "Convergence.dc.html", "todo", ""),
            ("Dynamics", "Dynamics.dc.html", "todo", ""),
            ("Analyze", "Analyze.dc.html", "todo", "")]
    out = []
    for n, h, s, d in base:
        if n == current:
            s = "current"
        out.append((n, h, s, d))
    return out


def shell(active, title, crumb, strip_html, content, status):
    body = (topbar(["CR_melt.caps"], "CR_melt.caps", f'{icon("chevr", 12, DIM)}<span>{crumb}</span>')
            + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail2(active)
            + f'<div style="flex-grow: 1; display: flex; flex-direction: column; min-width: 0">{strip_html}{content}</div></div>' + statusbar(*status))
    return page(title, body)


def viewbar(*items):
    return (f'<div style="height: 44px; flex-shrink: 0; display: flex; align-items: center; gap: 4px; padding: 0 10px; background: {BG1}; border-bottom: 1px solid {LINE}">'
            + "".join(items) + '</div>')


def spacer():
    return '<div style="flex-grow: 1"></div>'


def okrow(text, state="ok", sub=""):
    ic, col = {"ok": ("check", OK), "warn": ("alert", WARN), "err": ("xcircle", ERR), "info": ("link", MUTED)}[state]
    s = f'<span style="display: block; font-size: 11.5px; color: {DIM}; margin-top: 1px">{sub}</span>' if sub else ""
    return (f'<div style="display: flex; gap: 8px; align-items: flex-start; font-size: 12.5px; line-height: 1.4">{icon(ic, 15, col, 2.2)}'
            f'<span>{text}{s}</span></div>')


# ------------------------------------------------------------------ board 1: Grow with the force field in it
def pipeline_grow(bg=BG0):
    cw = 1440 - 72 - 340
    vh = 900 - 44 - 26 - 56 - 44 - 176 - 29   # less the rules between the bars
    sc = cell_scene("pg", cw, vh, "chains", bg=bg)
    chains = chains_overlay(sc, bg)
    bar10 = 10 * sc.scale
    ov = chains + (f'<line x1="16" y1="{vh - 18}" x2="{16 + bar10:.1f}" y2="{vh - 18}" stroke="{MUTED}" stroke-width="2"></line>'
          f'<text x="{16 + bar10 / 2:.1f}" y="{vh - 24}" text-anchor="middle" font-size="10.5" font-family="IBM Plex Mono" fill="{MUTED}">10 Å</text>')
    legend = (f'<div style="position: absolute; left: 14px; top: 12px; display: flex; flex-direction: column; gap: 6px; padding: 10px 12px; background: {BG1}e6; border: 1px solid {LINE}; border-radius: 8px; font-size: 12px">'
              f'<span style="font-weight: 600">Polychloroprene melt · as grown</span>'
              f'<span style="font-family: {MONO}; font-size: 11px; color: {MUTED}">10 chains × 30 units · 3020 atoms · 41.90 Å · 0.600 g/cm³</span>'
              f'<span style="display: flex; gap: 10px; align-items: center; font-size: 11px; color: {MUTED}">{dot(CHAIN[0])}{dot(CHAIN[1])}{dot(CHAIN[2])}<span>backbone, one colour per chain</span>{dot("#57B26A")}<span>Cl</span></span></div>')
    view = (f'<div style="position: relative; flex-grow: 1; min-height: 0; background: {bg}; overflow: hidden">{sc.svg(overlay=ov)}{legend}</div>')
    vb = viewbar(tbtn("rotate", "Orbit", True), tbtn("move", "Pan"), tbtn("search", "Zoom"), sep(), tbtn("cursor", "Pick"), tbtn("ruler", "Measure"), spacer(),
                 seg(["Backbone", "Ball and stick", "Force-field types"], "Backbone"))
    # project state: what this project holds now
    def tile(label, value, sub, col=TEXT, ic="check", icol=OK):
        return (f'<div style="flex: 1 1 0; min-width: 0; display: flex; flex-direction: column; gap: 4px; padding: 12px 14px; background: {BG2}; border: 1px solid {LINE}; border-radius: 8px">'
                f'<span style="display: flex; align-items: center; gap: 6px; font-size: 11px; color: {MUTED}; text-transform: uppercase; letter-spacing: 0.06em; font-weight: 600">{icon(ic, 13, icol, 2.2)}{label}</span>'
                f'<span style="font-size: 14px; font-weight: 600; color: {col}; white-space: nowrap; overflow: hidden; text-overflow: ellipsis">{value}</span>'
                f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}; white-space: nowrap; overflow: hidden; text-overflow: ellipsis">{sub}</span></div>')
    proj = (f'<div style="height: 176px; flex-shrink: 0; display: flex; flex-direction: column; gap: 10px; padding: 12px 14px; background: {BG1}; border-top: 1px solid {LINE}">'
            f'<div style="display: flex; align-items: center; gap: 10px"><h3 style="{H3}">Project now</h3><span style="font-size: 12px; color: {DIM}">what the next step starts from · every item can be exported as it is</span>{spacer()}'
            f'<a href="ExportCenter.dc.html" style="font-size: 12px">Open the Export center</a></div>'
            f'<div style="display: flex; gap: 10px">'
            + tile("Structure", "3020 atoms", "10 chains × 30 units")
            + tile("Force field", "GAFF (AmberTools 25)", "auto after Grow · 5 types", OK)
            + tile("Charges", "Gasteiger", "net 0.000 e")
            + tile("Cell", "41.90 Å cubic", "0.600 g/cm³, loose", TEXT, "alert", WARN)
            + tile("Next", "Relax", "minimise, then NPT", ACC, "chevr", ACC)
            + '</div>'
            f'<div style="display: flex; align-items: center; gap: 8px; font-size: 11.5px; color: {DIM}; font-family: {MONO}">'
            f'<span>history</span><span style="color: {MUTED}">Build unit</span>{icon("chevr", 11, DIM)}<span style="color: {MUTED}">Grow 10 × 30 (seed 3)</span>{icon("chevr", 11, DIM)}'
            f'<span style="color: {OK}">Force field: GAFF (AmberTools 25), auto</span>{icon("chevr", 11, DIM)}<span>Relax …</span></div></div>')
    center = f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">{vb}{view}{proj}</div>'
    # the Grow inspector: chain settings, then the force field it is assigned with
    ffcard = (f'<div style="display: flex; flex-direction: column; gap: 10px; padding: 12px; background: {BG2}; border: 1px solid {LINE}; border-radius: 8px">'
              f'<div style="display: flex; align-items: center; gap: 8px">{icon("tag", 16, ACC)}<span style="font-size: 13px; font-weight: 600">Force field</span>{spacer()}{chip("assigned after growing", OK, BG3)}</div>'
              + select("Force field", "GAFF (AmberTools 25)") + select("Charges", "Automatic (force field, else Gasteiger)")
              + okrow("All 3020 atoms typed, every term found", "ok", "c3 · c2 · hc · ha · cl")
              + okrow("OPLS-AA 2024 would miss 3 terms", "warn", "the C–C(Cl)= angle and a CM–CT–CT–CM torsion")
              + f'<a href="ForceFieldStep.dc.html" style="font-size: 12px">Compare all force fields for this structure</a></div>')
    right = (f'<aside style="width: 340px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">'
             + panel_head("Grow", chip("step 2 of 8", MUTED, BG2))
             + section("Repeat unit", col(row(field("SMILES", "[*]C/C(Cl)=C/C[*]"), gap=8), row(kv("Unit", "C₄H₅Cl · 88.53 g/mol"), gap=8), gap=8))
             + section("Chains", col(row(field("Chains", "10"), field("Units per chain", "30"), gap=8), row(field("Density", "0.600", "g/cm³"), select("Linkage", "trans-1,4"), gap=8), gap=8))
             + f'<div style="padding: 12px 14px; border-bottom: 1px solid {LINE}">{ffcard}</div>'
             + f'<div style="padding: 12px 14px; display: flex; flex-direction: column; gap: 8px"><div style="display: flex; gap: 8px">{btn("Grow again", False, "grow")}<div style="flex-grow: 1; display: flex">{btn("Next: Relax", True, "chevr", False, "Relax.dc.html").replace("display: inline-flex", "display: flex; flex-grow: 1")}</div></div>'
             f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.4">Pack is optional: Grow made the cell. Use it for solvent, filler or a blend.</span></div></aside>')
    content = f'<div style="flex-grow: 1; display: flex; min-height: 0">{center}{right}</div>'
    return shell("Grow", "CAPS — pipeline: grow with the force field", "Grow", strip(steps_for("Grow"), "Next: Relax", "Relax.dc.html"), content,
                 ("<span>CR_melt · 3020 atoms · GAFF (AmberTools 25) · Gasteiger</span>", "<span>pipeline: Build ✓ · Grow ✓ · Force field ✓ · Relax next</span>"))


# ------------------------------------------------------------------ board 2: the Force field step
def short_type(t):
    m = re.match(r"^(\d+)_b([A-Za-z0-9]+)_a", t)
    return f"{m.group(1)} ({m.group(2)})" if m else t


def ff_step():
    by_id = {f["id"]: f for f in TYPED}
    sel = by_id["oplsaa2024-moltemplate"]
    order = {"complete": 0, "missing parameters": 1, "untyped atoms": 2}
    rows = sorted(TYPED, key=lambda f: (order.get(f["status"], 3), f["name"]))
    trs = []
    for f in rows:
        st = f["status"]
        on = f["id"] == sel["id"]
        if st == "complete":
            badge = chip(f'{icon("check", 12, OK, 2.6)}complete', OK, "#1F2B22")
        elif st == "missing parameters":
            badge = chip(f'{f["missing_count"]} terms missing', WARN, "#2E2618")
        else:
            badge = chip(f'{f["untyped"]} atoms untyped', ERR, "#2E1E1E")
        typed = 3020 - (f.get("untyped") or 0)
        ch = {"types": "force field", "gasteiger": "Gasteiger", "none": "none (UFF)"}.get(f.get("charges"), "—") if st != "untyped atoms" else "—"
        bgc = f"background: {BG3}; box-shadow: inset 3px 0 0 {ACC};" if on else ""
        trs.append(f'<tr style="height: 34px; border-bottom: 1px solid {BG2}; {bgc}"><td style="padding: 0 12px; font-weight: {600 if on else 400}">{f["name"]}</td>'
                   f'<td style="font-family: {MONO}; text-align: right; padding-right: 14px; color: {TEXT if typed == 3020 else ERR}">{typed}</td>'
                   f'<td style="padding-right: 10px; color: {MUTED}">{ch}</td><td style="padding-right: 12px; text-align: right">{badge}</td></tr>')
    ncomplete = sum(1 for f in TYPED if f["status"] == "complete")
    table_html = (f'<table style="width: 100%; border-collapse: collapse; font-size: 12.5px"><thead><tr style="height: 32px; color: {MUTED}; font-size: 11.5px; text-align: left; border-bottom: 1px solid {LINE}">'
                  f'<th style="padding: 0 12px; font-weight: 500">Force field</th><th style="font-weight: 500; text-align: right; padding-right: 14px">Atoms typed</th><th style="font-weight: 500">Charges</th><th style="font-weight: 500; text-align: right; padding-right: 12px">For this structure</th></tr></thead>'
                  f'<tbody>{"".join(trs)}</tbody></table>')
    left = (f'<div style="width: 600px; flex-shrink: 0; display: flex; flex-direction: column; border-right: 1px solid {LINE}; background: {BG1}">'
            f'<div style="padding: 14px 16px 10px; display: flex; flex-direction: column; gap: 4px; border-bottom: 1px solid {LINE}"><h2 style="margin: 0; font-size: 15px; font-weight: 600">Which force fields fit this structure</h2>'
            f'<span style="font-size: 12px; color: {MUTED}">Each force field\'s typing rules and parameters tried on all 3020 atoms · {ncomplete} of {len(TYPED)} complete · about 6 s</span></div>'
            f'<div style="flex-grow: 1; overflow: hidden">{table_html}</div>'
            f'<div style="padding: 10px 16px; border-top: 1px solid {LINE}; font-size: 11.5px; color: {DIM}; line-height: 1.45">{len(NO_RULES)} more force fields in the library have no CAPS typing rules (CHARMM36 proteins, TraPPE-UA, MARTINI, water models …): assign those by hand in the Types table, or type a matching structure.</div></div>')
    # viewport: the atoms the missing terms involve
    vw, vh = 1440 - 72 - 600 - 340 - 4, 360
    # one chain, close up: the atoms of the terms OPLS-AA 2024 lacks in colour, the rest dimmed
    ids = [a[0] for a in ATOMS if a[1] == 4]
    idx = {a: k for k, a in enumerate(ids)}
    recs = []
    for a in ids:
        if EL[a] == "Cl":
            c = "#57B26A"
        elif EL[a] == "C" and any(EL[j] == "Cl" for j in NB[a]):
            c = ACC
        else:
            c = mix(ELEM[EL[a]][0], BG0, 0.45)
        recs.append(dict(p=POS[a], e=EL[a], c=c))
    sc = Scene("ff", vw, vh, yaw=0.3, pitch=0.25, persp=0.0, fog=0.35, atom_k=1.5, bond_w=0.22)
    sc.add_atoms(recs, [(idx[i], idx[j]) for i, j in BONDS if i in idx and j in idx])
    view = (f'<div style="position: relative; height: {vh}px; flex-shrink: 0; background: {BG0}; overflow: hidden; border-bottom: 1px solid {LINE}">{sc.svg()}'
            f'<div style="position: absolute; left: 10px; top: 10px; display: flex; flex-direction: column; gap: 4px; padding: 8px 10px; background: {BG1}e6; border: 1px solid {LINE}; border-radius: 7px; font-size: 11.5px">'
            f'<span style="font-weight: 600">Where OPLS-AA 2024 falls short · chain 4</span><span style="display: flex; gap: 6px; align-items: center; color: {MUTED}">{dot(ACC)}C of C(Cl)={dot("#57B26A")}Cl<span style="color: {DIM}">· 300 of each in the cell</span></span></div></div>')
    how = (f'<div style="flex-grow: 1; min-height: 0; display: flex; flex-direction: column; gap: 10px; padding: 14px 16px; background: {BG0}">'
           f'<h3 style="{H3}">How CAPS decides</h3>'
           + okrow("Types from the chemical environment", "info", "element, bonded neighbours, bond orders, rings and aromaticity, matched against each force field's rules; the most specific rule wins")
           + okrow("Every bonded term looked up", "info", "bonds, angles, torsions and impropers by the types' classes and the file's wildcards; a term with no entry is reported, never guessed")
           + okrow("Charges must add up", "info", "per-type charges have to sum to the structure's formal charge; otherwise the assignment is incomplete")
           + '</div>')
    mid = f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">{view}{how}</div>'
    miss = []
    for m in sel["missing"]:
        k, *ts = m.split()
        miss.append(f'<div style="display: flex; flex-direction: column; gap: 2px; padding: 8px 10px; background: {BG2}; border: 1px solid {LINE}; border-radius: 6px">'
                    f'<span style="font-size: 12px; font-weight: 600">{k}</span><span style="font-family: {MONO}; font-size: 11.5px; color: {MUTED}">{" – ".join(short_type(t) for t in ts)}</span></div>')
    alts = [f for f in rows if f["status"] == "complete" and f["id"] in ("gaff-amber25", "cvff", "pcff", "dreiding-moltemplate")]
    alt_html = "".join(f'<div style="display: flex; align-items: center; gap: 8px; font-size: 12.5px">{icon("check", 14, OK, 2.4)}<span style="flex-grow: 1">{f["name"]}</span>{btn("Use", False, None, True)}</div>' for f in alts)
    right = (f'<aside style="width: 340px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">'
             + panel_head("OPLS-AA (2024 parameter file)", chip("3 terms missing", WARN, "#2E2618"))
             + section("Typing", col(okrow("3020 of 3020 atoms typed", "ok", "135 · 136 · 140 · 142 · 144 · 226 · 399"), okrow("Charges from the force field balance", "ok", "net 0.000 e: Cl −0.12 against C(Cl)= +0.12"), gap=8))
             + section("Not in this force field", col(*miss, f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.45">The 2024 table has no C–C(Cl)= angle and no CM–CT–CT–CM torsion (as in polybutadiene too). A run would leave these terms out, so CAPS stops here.</span>', gap=6))
             + section("Complete for this structure", col(alt_html, gap=8))
             + f'<div style="padding: 12px 14px; display: flex; flex-direction: column; gap: 8px">{btn("Add the 3 terms by hand", False, "plus")}'
             f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.45">Entered terms are saved with the project and marked as yours in every exported file.</span></div></aside>')
    content = f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{mid}{right}</div>'
    return shell("Force field", "CAPS — force field step: coverage", "Force field", strip(steps_for("Force field"), "Next: Relax", "Relax.dc.html"), content,
                 (f"<span>CR_melt · 3020 atoms · {ncomplete} complete · {sum(1 for f in TYPED if f['status'] == 'missing parameters')} missing terms · {sum(1 for f in TYPED if f['status'] == 'untyped atoms')} untyped</span>",
                  "<span>in use: GAFF (AmberTools 25)</span>"))


# ------------------------------------------------------------------ board 3: Export center
def export_center():
    pre_in = open(os.path.join(DATA, "cr_out", "system.in")).read().split("\n")
    shown = [l for l in pre_in if l.strip()][2:40]
    def hl(l):
        e = esc(l)
        if l.startswith("#"):
            return f'<span style="color: {DIM}">{e}</span>'
        m = re.match(r"^(\S+)(\s+)(.*)$", l)
        if not m:
            return e
        rest = re.sub(r"(#.*)$", lambda x: f'<span style="color: {DIM}">{x.group(1)}</span>', esc(m.group(3)))
        return f'<span style="color: {SEL}">{esc(m.group(1))}</span>{m.group(2)}{rest}'
    code = "".join(f'<div style="display: flex; gap: 12px"><span style="width: 22px; flex-shrink: 0; text-align: right; color: {DIM}">{k + 1}</span><span style="white-space: pre">{hl(l)}</span></div>' for k, l in enumerate(shown))
    ftabs = tabs(["system.in", "system.data", "system.top", "system.itp", "system.gro", "system.mdp"], "system.in", 12)
    preview = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px; overflow: hidden">'
               f'{ftabs}<div style="flex-grow: 1; min-height: 0; overflow: auto; padding: 10px 12px; background: {BG0}; font-family: {MONO}; font-size: 11.5px; line-height: 1.55">{code}</div></div>')

    def fileline(name, size, what):
        return (f'<div style="display: flex; align-items: center; gap: 8px; font-size: 12.5px">{icon("file", 14, MUTED)}<span style="font-family: {MONO}">{name}</span>'
                f'<span style="flex-grow: 1; color: {DIM}; font-size: 11.5px; overflow: hidden; white-space: nowrap; text-overflow: ellipsis">{what}</span><span style="font-family: {MONO}; font-size: 11.5px; color: {MUTED}">{size}</span></div>')
    lmp = card("LAMMPS", col(
        check("Write for LAMMPS", True),
        fileline("system.data", fsize("system.data"), "atoms, bonds, masses, bonded coefficients"),
        fileline("system.in", fsize("system.in"), "styles, all 15 pair_coeff lines, protocol"),
        f'<div style="height: 1px; background: {LINE}"></div>',
        seg(["Check", "Minimise", "NVT", "NPT"], "NPT", full=True),
        row(field("T", "300", "K"), field("P", "1", "atm"), field("Steps", "100000"), gap=8),
        toggle("Minimise first", True), gap=10), chip("lj/cut/coul/dsf · fourier", MUTED, BG2, True))
    gmx = card("GROMACS", col(
        check("Write for GROMACS", True),
        fileline("system.top", fsize("system.top"), "defaults, atom types, all pairs"),
        fileline("system.itp", fsize("system.itp"), "1 molecule type × 10 chains"),
        fileline("system.gro", fsize("system.gro"), "coordinates, nm, chains whole"),
        fileline("system.mdp", fsize("system.mdp"), "matching cut-offs and tail"),
        okrow("Coulomb: damped shifted force becomes PME", "warn", "GROMACS has no DSF; energies agree to the method difference"), gap=10), chip("top · itp · gro · mdp", MUTED, BG2, True))
    other = card("Structure only", col(
        check("PDB", False), check("MOL2 (with types and charges)", False), check("XYZ", False), check("CAPS project (.caps)", False),
        f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.45">Any step can be exported: the files hold the structure as it is now (grown, not yet relaxed).</span>', gap=9))
    top = f'<div style="display: grid; grid-template-columns: 1.15fr 1.15fr 0.9fr; gap: 12px">{lmp}{gmx}{other}</div>'
    checks = (f'<aside style="width: 330px; flex-shrink: 0; display: flex; flex-direction: column; gap: 10px; padding: 14px; background: {BG1}; border: 1px solid {LINE}; border-radius: 10px">'
              f'<h3 style="{H3}">Checked before writing</h3>'
              + okrow("3020 of 3020 atoms typed", "ok", "GAFF (AmberTools 25): c3 c2 hc ha cl")
              + okrow("Every bond, angle and torsion has parameters", "ok", "3010 bonds · 5400 angles · 5030 torsions · 600 impropers")
              + okrow("Net charge 0.000 e", "ok", "Gasteiger–Marsili")
              + okrow("All 15 type pairs written", "ok", "mixing applied by CAPS")
              + okrow("0.600 g/cm³, not yet relaxed", "warn", "the NPT protocol takes it to density")
              + f'<div style="flex-grow: 1"></div>'
              + row(field("Folder", "~/CAPS/CR_melt/export", None, True), field("Stem", "system", None, True, 76), gap=8)
              + btn("Write 6 files", True, "download")
              + '</aside>')
    head = (f'<div style="display: flex; align-items: flex-end; gap: 16px"><div style="display: flex; flex-direction: column; gap: 3px"><h1 style="margin: 0; font-size: 19px; font-weight: 600">Export</h1>'
            f'<span style="font-size: 12.5px; color: {MUTED}">CR_melt as it is now · after Grow and the force field · for LAMMPS, GROMACS or as a structure</span></div></div>')
    content = (f'<div style="flex-grow: 1; min-height: 0; display: flex; flex-direction: column; gap: 12px; padding: 16px 18px; background: {BG0}">{head}{top}'
               f'<div style="flex-grow: 1; min-height: 0; display: flex; gap: 12px">{preview}{checks}</div></div>')
    st = [(n, h, s if n != "Analyze" else s, d) for n, h, s, d in steps_for("")]
    return shell("Export", "CAPS — Export center", "Export", strip(st, "Back to Relax", "Relax.dc.html"), content,
                 (f"<span>CR_melt · GAFF (AmberTools 25) · 6 files · {TOTAL}</span>", "<span>export is open at every step</span>"))


# ------------------------------------------------------------------ Paper (light) variant
PAPER = {BG0: "#F4F2EE", BG1: "#FFFFFF", BG2: "#F1EFEB", BG3: "#E7E4DF", LINE: "#D9D5CF", TEXT: "#1E2226", MUTED: "#555B62", DIM: "#6B7178",
         ACC: "#A35C00", ACC_INK: "#FFFFFF", SEL: "#0F7C93", OK: "#2E7D32", ERR: "#C62828", "#1F2B22": "#E3F0E3", "#2E2618": "#F6EAD6", "#2E1E1E": "#F8E1DF"}


def paper(html):
    for a, b in PAPER.items():
        html = html.replace(a, b)
    return html


if __name__ == "__main__":
    os.makedirs("stage28/project", exist_ok=True)
    for name, fn in (("PipelineGrow", pipeline_grow), ("ForceFieldStep", ff_step), ("ExportCenter", export_center)):
        h = fn(); open(f"stage28/project/{name}.dc.html", "w").write(h); print(name, len(h))
    h = paper(pipeline_grow(bg="#F4F2EE")); open("stage28/project/PaperPipelineGrow.dc.html", "w").write(h); print("PaperPipelineGrow", len(h))
