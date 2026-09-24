import math, random
from lib import *
from mols import *
from screen_app import card, pagehead, shell
from screen_builders import footer, hud
from screen_studio2 import bar, spacer, styrene3d
from screen_row7 import with_overlay, modal
import screen_main


def small_scene(sid, atoms, w=150, h=100, bonds=None, yaw=0.6, pitch=0.5):
    sc = Scene(sid, w, h, yaw=yaw, pitch=pitch, persp=0.15, fog=0.3, outline=True)
    sc.add_atoms(atoms, bonds)
    return sc.svg()


def benzene():
    a = [{"e": "C", "p": (1.39 * math.cos(math.radians(60 * k)), 1.39 * math.sin(math.radians(60 * k)), 0)} for k in range(6)]
    a += [{"e": "H", "p": (2.47 * math.cos(math.radians(60 * k)), 2.47 * math.sin(math.radians(60 * k)), 0)} for k in range(6)]
    return a


def cyclohexane():
    R = math.sqrt(1.53 ** 2 - 0.25)  # chair with +-0.25 A puckering gives 1.53 A C-C
    return [{"e": "C", "p": (R * math.cos(math.radians(60 * k)), R * math.sin(math.radians(60 * k)), 0.25 if k % 2 else -0.25)} for k in range(6)]


def methane():
    d = 1.09 / math.sqrt(3)
    return [{"e": "C", "p": (0, 0, 0)}] + [{"e": "H", "p": (sx * d, sy * d, sx * sy * d)} for sx, sy in ((1, 1), (1, -1), (-1, 1), (-1, -1))]


def water1():
    a = math.radians(104.52 / 2)
    return [{"e": "O", "p": (0, 0, 0)}, {"e": "H", "p": (0.9572 * math.sin(a), 0.9572 * math.cos(a), 0)}, {"e": "H", "p": (-0.9572 * math.sin(a), 0.9572 * math.cos(a), 0)}]


def ester():
    # methyl acetate heavy atoms, planar approximate: C-C 1.51, C=O 1.20, C-O 1.34, O-CH3 1.44
    def pol(o, r, deg):
        return (o[0] + r * math.cos(math.radians(deg)), o[1] + r * math.sin(math.radians(deg)), 0.0)
    c1 = (0, 0, 0)
    c2 = pol(c1, 1.51, 0)
    o1 = pol(c2, 1.20, 60)
    o2 = pol(c2, 1.34, -60)
    c3 = pol(o2, 1.44, 0)
    return [{"e": "C", "p": c1}, {"e": "C", "p": c2}, {"e": "O", "p": o1}, {"e": "O", "p": o2}, {"e": "C", "p": c3}]


def amide():
    def pol(o, r, deg):
        return (o[0] + r * math.cos(math.radians(deg)), o[1] + r * math.sin(math.radians(deg)), 0.0)
    c1 = (0, 0, 0)
    c2 = pol(c1, 1.51, 0)
    o1 = pol(c2, 1.23, 60)
    n = pol(c2, 1.33, -60)
    c3 = pol(n, 1.45, 0)
    return [{"e": "C", "p": c1}, {"e": "C", "p": c2}, {"e": "O", "p": o1}, {"e": "N", "p": n}, {"e": "C", "p": c3}]


# ------------------------------------------------------------------ Fragment library
def fragments():
    items = [("Benzene", benzene(), "c1ccccc1", "1 attach · aromatic"), ("Styrene", styrene3d(), "C=Cc1ccccc1", "monomer · head/tail set"),
             ("Cyclohexane", cyclohexane(), "C1CCCCC1", "chair · H hidden"), ("Methyl", methane(), "[CH3]*", "1 attach"),
             ("Ester", ester(), "*C(=O)O*", "2 attach"), ("Amide", amide(), "*C(=O)N*", "2 attach · trans"),
             ("Water", water1(), "O", "solvent · TIP4P/2005"), ("Ethylene", [{"e": "C", "p": (0, 0, 0)}, {"e": "C", "p": (1.34, 0, 0)}], "C=C", "monomer")]
    tiles = []
    for i, (n, at, smi, meta) in enumerate(items):
        on = n == "Ester"
        tiles.append(f'<button style="display: flex; flex-direction: column; padding: 0; text-align: left; background: {BG1}; border: 1.5px solid {ACC if on else LINE}; border-radius: 10px; overflow: hidden; cursor: pointer">'
                     f'<div style="background: {BG0}; border-bottom: 1px solid {LINE}; display: flex; justify-content: center">{small_scene("fr" + str(i), at, 170, 110)}</div>'
                     f'<span style="padding: 10px 12px; display: flex; flex-direction: column; gap: 3px"><span style="font-size: 13px; font-weight: 600; color: {TEXT}">{n}</span>'
                     f'<span style="font-family: {MONO}; font-size: 11px; color: {MUTED}">{esc(smi)}</span><span style="font-size: 11px; color: {DIM}">{meta}</span></span></button>')
    cats = [("All", 1240), ("Rings", 86), ("Functional groups", 142), ("Monomers", 124), ("Amino acids", 20), ("Nucleotides", 8), ("Sugars", 24),
            ("Solvents", 48), ("Ions", 36), ("Ligands", 410), ("My fragments", 12)]
    cl = "".join(f'<a href="#" aria-current="{"page" if n == "Functional groups" else "false"}" style="display: flex; align-items: center; justify-content: space-between; height: 32px; padding: 0 12px; border-radius: 6px; text-decoration: none; font-size: 12.5px; background: {BG3 if n == "Functional groups" else "transparent"}; color: {TEXT if n == "Functional groups" else MUTED}"><span>{n}</span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{c}</span></a>' for n, c in cats)
    left = (f'<nav aria-label="Categories" style="width: 230px; flex-shrink: 0; display: flex; flex-direction: column; gap: 2px; padding: 14px 10px; background: {BG1}; border-right: 1px solid {LINE}">'
            f'<div style="display: flex; align-items: center; gap: 8px; height: 32px; padding: 0 8px; margin-bottom: 8px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px; color: {DIM}; font-size: 12px">{icon("search", 14, DIM)}<span>Search name or SMARTS</span></div>{cl}</nav>')
    grid = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 14px; padding: 20px 24px">'
            + row(f'<h1 style="margin: 0; font-size: 19px; font-weight: 600">Functional groups</h1>', chip("142"), spacer(), seg(["Grid", "List"], "Grid"), gap=10)
            + f'<div style="display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px">{"".join(tiles)}</div></div>')
    est = ester()
    sc = Scene("frd", 360, 220, yaw=0.3, pitch=0.3, persp=0.1, fog=0.2)
    sc.add_atoms(est)
    sc._prep()
    ov = ""
    for idx in (0, 4):
        x, y, _, _ = sc.screen(idx)
        ov += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="16" fill="none" stroke="{ACC}" stroke-width="2" stroke-dasharray="4 3"></circle>'
    ov += label_pill(sc.screen(0)[0] - 20, sc.screen(0)[1] + 30, "attach 1", ACC, border=ACC) + label_pill(sc.screen(4)[0] - 20, sc.screen(4)[1] + 30, "attach 2", ACC, border=ACC)
    detail = (panel_head("Ester", chip("2 attach points", ACC, "#3A2C14"))
              + f'<div style="background: {BG0}; border-bottom: 1px solid {LINE}">{sc.svg(overlay=ov)}</div>'
              + section("Definition", col(kv("SMILES", "*C(=O)O*"), kv("Attach 1", "carbonyl carbon"), kv("Attach 2", "ester oxygen"), kv("Geometry", "planar, s-cis ester"), gap=7))
              + section("Attach behaviour", col(toggle("Replace a hydrogen at the target", True), toggle("Align to open valence", True), toggle("Clean geometry after attaching", True), gap=9))
              + footer(btn("Save a copy to My fragments", ic="save"), btn("Attach to selection", True, "link")))
    body = (topbar(["PS_atactic_DP40.caps"], "PS_atactic_DP40.caps") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + f'<div style="flex-grow: 1; display: flex; min-width: 0">{left}{grid}<aside aria-label="Fragment" style="width: 360px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{detail}</aside></div></div>'
            + statusbar("<span>1 240 fragments · 12 yours</span>", "<span>Save any selection as a fragment: ⌘⇧S</span>"))
    return page("CAPS Studio — fragment library", body)


# ------------------------------------------------------------------ Macro recorder
def macro():
    code = [
        ("c", "# Recorded in Studio · 2026-09-24 · replayable"),
        ("k", "import caps"),
        ("", ""),
        ("k", "def build_ps(dp: int = 40, tacticity: str = \"atactic\", seed: int = 1):"),
        ("", "    doc = caps.Document.new(\"PS\")"),
        ("", "    unit = doc.chem.from_smiles(\"*CC(*)c1ccccc1\", head=1, tail=2)"),
        ("", "    chain = doc.polymer.build(unit, dp=dp, tacticity=tacticity, seed=seed)"),
        ("", "    doc.field.assign(\"gaff2\")"),
        ("", "    doc.relax.minimise(method=\"lbfgs\", fmax=0.02)"),
        ("", "    return doc"),
        ("", ""),
        ("k", "for t in (\"isotactic\", \"syndiotactic\", \"atactic\"):"),
        ("", "    build_ps(tacticity=t).save(f\"PS_{t}_DP40.caps\")"),
    ]
    lines = []
    for i, (kind, t) in enumerate(code, start=1):
        colr = DIM if kind == "c" else TEXT
        tt = esc(t)
        for kw in ("import ", "def ", "for ", "return ", " in "):
            tt = tt.replace(kw, f'<span style="color: {SEL}">{kw}</span>')
        tt = tt.replace("&quot;", '"')
        hl = f"background: {BG3}; " if i in (5, 6, 7) else ""
        lines.append(f'<div style="display: flex; {hl}"><span style="width: 36px; flex-shrink: 0; text-align: right; padding-right: 12px; color: {DIM}">{i}</span><span style="color: {colr}; white-space: pre">{tt}</span></div>')
    editor = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; background: {BG0}">'
              + tabs(["build_ps.py", "anneal.py"], "build_ps.py")
              + f'<div style="flex-grow: 1; padding: 14px 0; font-family: {MONO}; font-size: 13px; line-height: 1.75">{"".join(lines)}</div>'
              + f'<div style="height: 170px; flex-shrink: 0; border-top: 1px solid {LINE}; background: {BG1}; display: flex; flex-direction: column">'
              + tabs(["Output", "Problems"], "Output")
              + f'<div style="padding: 10px 14px; font-family: {MONO}; font-size: 12px; line-height: 1.7"><div style="color: {MUTED}">&gt;&gt;&gt; run build_ps.py</div>'
                f'<div>PS_isotactic_DP40.caps   642 atoms   typed   minimised</div><div>PS_syndiotactic_DP40.caps   642 atoms   typed   minimised</div>'
                f'<div style="color: {ACC}">PS_atactic_DP40.caps   running · relax 61 %</div></div></div></div>')
    rec = [("●", "Recording", ERR), ("1", "chem.from_smiles  *CC(*)c1ccccc1", MUTED), ("2", "polymer.build  dp=40  tacticity=atactic", MUTED), ("3", "field.assign  gaff2", MUTED), ("4", "relax.minimise  lbfgs", MUTED), ("5", "document.save", MUTED)]
    rl = "".join(f'<div style="display: flex; gap: 10px; align-items: center; height: 28px; font-family: {MONO}; font-size: 11.5px"><span style="width: 16px; color: {c}">{n}</span><span style="color: {TEXT if n != "●" else ERR}">{esc(t)}</span></div>' for n, t, c in rec)
    params = [("dp", "int", "40"), ("tacticity", "str", "atactic"), ("seed", "int", "1")]
    pr = "".join(f'<div style="display: flex; align-items: center; gap: 10px; height: 32px; padding: 0 10px; background: {BG0}; border: 1px solid {LINE}; border-radius: 6px; font-family: {MONO}; font-size: 12px"><span style="color: {ACC}">{n}</span><span style="color: {DIM}">{t}</span><span style="margin-left: auto">{v}</span></div>' for n, t, v in params)
    right = (panel_head("Recorder", row(f'<button aria-label="Stop recording" style="width: 28px; height: 28px; border-radius: 50%; border: 1px solid {ERR}; background: transparent; display: flex; align-items: center; justify-content: center; cursor: pointer"><span style="width: 10px; height: 10px; border-radius: 2px; background: {ERR}"></span></button>', gap=6))
             + section("Recorded commands", col(rl, gap=2))
             + section("Parameters", col(pr, f'<div style="font-size: 12px; color: {MUTED}">Select a literal in the editor and choose <b style="font-weight: 600; color: {TEXT}">Promote to parameter</b>.</div>', gap=6))
             + section("Run", col(row(select("Target", "3 new documents"), select("Where", "This Mac"), gap=10), toggle("Stop on first error", True), gap=10))
             + footer(btn("Save as tool", ic="save"), btn("Run", True, "play")))
    body = (topbar(["build_ps.py"], "build_ps.py") + '<div style="flex-grow: 1; display: flex; min-height: 0">' + rail("Studio")
            + f'<div style="flex-grow: 1; display: flex; min-width: 0">{editor}<aside aria-label="Recorder" style="width: 380px; flex-shrink: 0; display: flex; flex-direction: column; background: {BG1}; border-left: 1px solid {LINE}; overflow: hidden">{right}</aside></div></div>'
            + statusbar("<span>Python 3.12 · caps 0.1.0</span>", "<span>Every Studio action records as one line</span>"))
    return page("CAPS Studio — macro recorder", body)


# ------------------------------------------------------------------ Import dialog
def import_dialog():
    base = screen_main.build()
    xyz = ["3000", 'Lattice="36.2 0 0 0 36.2 0 0 0 36.2" Properties=species:S:1:pos:R:3', "C   1.2040   0.4312   0.0000", "C   2.4740  -0.4288   0.0000",
           "H   1.2040   1.0612   0.8900", "H   1.2040   1.0612  -0.8900", "C   3.7440   0.4312   0.0000", "…"]
    pv = "".join(f'<div style="display: flex; gap: 12px"><span style="width: 22px; text-align: right; color: {DIM}">{i}</span><span style="white-space: pre">{esc(l)}</span></div>' for i, l in enumerate(xyz, 1))
    # perceived bonds preview: short PE segment
    pts = []
    for k in range(10):
        pts.append({"e": "C", "p": (k * 1.27, 0.43 * (1 if k % 2 == 0 else -1), 0.0)})
    for k in range(10):
        c = pts[k]["p"]
        s = 1 if k % 2 == 0 else -1
        pts.append({"e": "H", "p": (c[0], c[1] + s * 0.63, 0.89)})
        pts.append({"e": "H", "p": (c[0], c[1] + s * 0.63, -0.89)})
    sc = Scene("imp", 420, 200, yaw=0.3, pitch=0.5, persp=0.1, fog=0.3)
    sc.add_atoms(pts)
    nb = len(sc.bonds)
    prev = sc.svg()
    head = (f'<div style="height: 52px; display: flex; align-items: center; gap: 12px; padding: 0 20px; border-bottom: 1px solid {LINE}"><h2 style="margin: 0; font-size: 15px; font-weight: 600">Import PE_melt.xyz</h2>'
            f'{chip("Extended XYZ detected", OK, "#16261A")}{chip("no bonds in file", WARN, "#3A2C14")}{spacer()}{tbtn("close", "Close")}</div>')
    left = (f'<div style="width: 460px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px; padding: 18px 20px; border-right: 1px solid {LINE}">'
            f'<h3 style="margin: 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">File preview</h3>'
            f'<div style="padding: 10px 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px; font-family: {MONO}; font-size: 11.5px; line-height: 1.7; overflow: hidden">{pv}</div>'
            + kv("Atoms", "3 000") + kv("Cell (from Lattice=)", "36.2 × 36.2 × 36.2 Å, orthogonal") + kv("Units", "Å (XYZ convention)", False)
            + f'<h3 style="margin: 6px 0 0; font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase; color: {MUTED}">Bond perception preview · first 10 carbons</h3>'
            + f'<div style="background: {BG0}; border: 1px solid {LINE}; border-radius: 7px">{prev}</div>'
            + f'<span style="font-size: 11.5px; color: {DIM}">{nb} bonds perceived in the preview fragment</span></div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column">'
             + section("Bonds", col(seg(["Perceive from distances", "Read from file", "None"], "Perceive from distances", True),
                                     row(select("Radii", "Covalent · Cordero 2008"), field("Tolerance", "0.45", "Å"), gap=10),
                                     toggle("Assign bond orders and aromaticity after perceiving", True), gap=10))
             + section("Molecules &amp; cell", col(toggle("Split into molecules by connectivity", True), toggle("Unwrap molecules across periodic boundaries", True), toggle("Use cell from file", True), gap=9))
             + section("After import", col(row(select("Force field", "Type with TraPPE-UA"), select("Charges", "From force field"), gap=10), toggle("Run validation checks", True), gap=10))
             + '</div>')
    inner = (head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>'
             + f'<div style="display: flex; justify-content: flex-end; gap: 8px; padding: 12px 20px; border-top: 1px solid {LINE}">{btn("Cancel")}{btn("Import as new document", True, "download")}</div>')
    return with_overlay(base, modal(inner, 1040, 720)).replace("<title>CAPS Studio — builder</title>", "<title>CAPS Studio — import</title>")


# ------------------------------------------------------------------ Export image / movie
def export_dialog():
    base = screen_main.build()
    atoms, backbone, stereo = polystyrene(7, "atactic", seed=5)
    sc = Scene("exp", 560, 315, yaw=-0.28, pitch=0.42, roll=0.05, persp=0.25, fog=0.5, bg="#FFFFFF")
    sc.add_atoms(atoms)
    prev = sc.svg()
    checker = f'background-color: #FFFFFF; background-image: linear-gradient(45deg, #E9E9E9 25%, transparent 25%, transparent 75%, #E9E9E9 75%), linear-gradient(45deg, #E9E9E9 25%, transparent 25%, transparent 75%, #E9E9E9 75%); background-size: 16px 16px; background-position: 0 0, 8px 8px'
    head = (f'<div style="height: 52px; display: flex; align-items: center; gap: 12px; padding: 0 20px; border-bottom: 1px solid {LINE}"><h2 style="margin: 0; font-size: 15px; font-weight: 600">Export</h2>'
            + tabs(["Image", "Movie"], "Image", 13) + f'{spacer()}{tbtn("close", "Close")}</div>')
    left = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 10px; padding: 18px 20px; border-right: 1px solid {LINE}">'
            f'<div style="border: 1px solid {LINE}; border-radius: 8px; overflow: hidden; {checker}">{prev}</div>'
            + row(chip("preview · white background"), chip("16 : 9"), spacer(), chip("3840 × 2160 px", TEXT, BG3, True), gap=6)
            + f'<div style="font-size: 12px; color: {MUTED}; line-height: 1.5">Exports re-render at full resolution with the same camera, styles and labels. Measurement pills and monitors are optional overlays.</div></div>')
    presets = [("Screen", "1920 × 1080"), ("4K", "3840 × 2160"), ("Journal figure", "3.5 in · 600 dpi"), ("Poster", "7680 × 4320")]
    pg = "".join(f'<button style="display: flex; flex-direction: column; gap: 2px; padding: 8px 10px; text-align: left; border-radius: 6px; border: 1px solid {ACC if n == "4K" else LINE}; background: {BG3 if n == "4K" else BG0}; cursor: pointer"><span style="font-size: 12.5px; color: {TEXT}">{n}</span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{d}</span></button>' for n, d in presets)
    right = (f'<div style="width: 400px; flex-shrink: 0; display: flex; flex-direction: column">'
             + section("Size", col(f'<div style="display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 6px">{pg}</div>', row(field("Width", "3840", "px"), field("Height", "2160", "px"), gap=10), gap=10))
             + section("Output", col(row(select("Format", "PNG · 16-bit"), select("Background", "White"), gap=10), row(select("Supersampling", "4×"), select("Colour profile", "sRGB"), gap=10),
                                      toggle("Include labels and measurements", True), toggle("Embed provenance manifest in metadata", True), gap=10))
             + '</div>')
    inner = (head + f'<div style="flex-grow: 1; display: flex; min-height: 0">{left}{right}</div>'
             + f'<div style="display: flex; justify-content: flex-end; gap: 8px; padding: 12px 20px; border-top: 1px solid {LINE}">{btn("Copy to clipboard", ic="copy")}{btn("Export PNG", True, "download")}</div>')
    return with_overlay(base, modal(inner, 1040, 640)).replace("<title>CAPS Studio — builder</title>", "<title>CAPS Studio — export</title>")


# ------------------------------------------------------------------ Equilibration convergence
def convergence():
    # illustrative series
    blocks = list(range(1, 21))
    rho = [(b, 1.02 + 0.03 * (1 - math.exp(-b / 4)) + 0.002 * math.sin(b * 1.7)) for b in blocks]
    E = [(b, -2.0 - 1.2 * (1 - math.exp(-b / 5)) + 0.02 * math.sin(b * 2.1)) for b in blocks]
    Rg = [(b, 14.2 + 1.3 * (1 - math.exp(-b / 7)) + 0.08 * math.sin(b * 1.3)) for b in blocks]
    n = list(range(1, 200, 3))
    ideal = [(k, 9.5 * (1 - 1.9 / (k + 1.9 * 1.0)) if k > 0 else 0) for k in n]
    now = [(k, 9.5 * (1 - 1.9 / (k + 1.9)) * (1 + 0.12 * math.exp(-k / 60)) - 0.25 * (k > 120) * (k - 120) / 80) for k in n]
    band = [(16, 20, OK, 0.10)]
    p1 = plot(420, 160, [(rho, ACC, 1.8, None)], (0, 20), (1.00, 1.06), [0, 5, 10, 15, 20], [1.00, 1.03, 1.06], "block (0.5 ns)", "density (g/cm³)", pad=(44, 10, 16, 28), bands=band)
    p2 = plot(420, 160, [(E, SEL, 1.8, None)], (0, 20), (-3.4, -1.8), [0, 5, 10, 15, 20], [-3.4, -2.6, -1.8], "block (0.5 ns)", "E pot (10⁴ kcal/mol)", pad=(44, 10, 16, 28), bands=band)
    p3 = plot(420, 160, [(Rg, "#E07A5F", 1.8, None)], (0, 20), (14, 16), [0, 5, 10, 15, 20], [14, 15, 16], "block (0.5 ns)", "⟨R_g⟩ (Å)", pad=(44, 10, 16, 28), bands=band)
    p4 = plot(620, 240, [(ideal, MUTED, 1.6, "6 4"), (now, ACC, 2.2, None)], (0, 200), (0, 11), [0, 50, 100, 150, 200], [0, 2, 4, 6, 8, 10], "n (bonds apart)", "⟨R²(n)⟩ / (n l²)", pad=(44, 14, 20, 32))
    crit = [("Density", "slope &lt; 0.1 % per block for 3 blocks", "pass", OK), ("Potential energy", "slope &lt; 0.05 % per block for 3 blocks", "pass", OK),
            ("⟨R<sub>g</sub>⟩", "slope &lt; 1 % per block for 3 blocks", "pass", OK), ("Internal distances", "within 5 % of target curve for n ≤ 150", "not yet", WARN)]
    cr = "".join(f'<div style="display: flex; align-items: center; gap: 10px; padding: 10px 12px; background: {BG0}; border: 1px solid {LINE}; border-radius: 7px">{icon("check" if s == "pass" else "alert", 16, c)}'
                 f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 2px"><span style="font-size: 12.5px; font-weight: 500">{a}</span><span style="font-size: 11.5px; color: {MUTED}">{b}</span></div>'
                 f'<span style="font-family: {MONO}; font-size: 11.5px; color: {c}">{s}</span></div>' for a, b, s, c in crit)
    content = (pagehead("Equilibrate · convergence", "equil-7 · decides automatically when the cell is equilibrated",
                        row(chip(f'{dot(WARN)} 3 of 4 criteria met', WARN, "#3A2C14"), btn("Extend 5 ns", ic="play"), btn("Accept now", ic="check", href="Analyze.dc.html"), gap=8))
               + '<div style="flex-grow: 1; display: flex; gap: 14px; padding: 18px 22px; min-height: 0">'
               + f'<div style="width: 452px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
               + card("Density", p1, chip("pass", OK, "#16261A"), 12) + card("Potential energy", p2, chip("pass", OK, "#16261A"), 12) + card("Radius of gyration", p3, chip("pass", OK, "#16261A"), 12) + '</div>'
               + f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
               + card("Mean-square internal distances", col(p4, row(chip(f'{dot(MUTED)} target from RIS chain statistics'), chip(f'{dot(ACC)} current cell'), spacer(), chip("illustrative"), gap=6),
                                                             f'<div style="font-size: 12px; color: {MUTED}; line-height: 1.5">Long-range deviations at large n mean chains are not yet relaxed at the scale of the whole chain, even though density and energy have levelled off.</div>', gap=10), extra="flex-grow: 1")
               + card("Criteria", col(cr, gap=6))
               + '</div></div>')
    return shell("Equilibrate", "CAPS Equilibrate — convergence", f'{icon("chevr", 12, DIM)}<span>Equilibrate</span>', content,
                 ("<span>Block averages · 0.5 ns blocks</span>", "<span>Criteria and the decision go into the provenance manifest</span>"))


# ------------------------------------------------------------------ Provenance viewer
def provenance():
    steps = [("chem.build", "polymer from *CC(*)c1ccccc1 · DP 40 × 20", "RDKit 2025.09 · ETKDG v3", "riniker2015"),
             ("field.assign", "GAFF2 · AM1-BCC library charges", "caps-ffdb 1.0", "wang2004"),
             ("grow.cbmc", "k = 16 · ρ₀ 0.50 g/cm³ · seed 20260923", "Philox4x32-10", "siepmann1992 · rosenbluth1955"),
             ("relax.lbfgs", "|F|max 0.02 kcal/mol/Å · push-off on", "", "liu1989 · auhl2003"),
             ("equilibrate.larsen21", "21 steps · Pmax 50 000 bar", "", "larsen2011"),
             ("dynamics.npt", "10 ns · CSVR · cell rescaling · SPME 1e−5", "", "bussi2007 · bernetti2020 · essmann1995")]
    tl = []
    for i, (e, p, v, c) in enumerate(steps):
        diff = e == "grow.cbmc"
        tl.append(f'<div style="display: flex; gap: 14px">'
                  f'<div style="display: flex; flex-direction: column; align-items: center; width: 14px"><span style="width: 12px; height: 12px; border-radius: 50%; background: {WARN if diff else ACC}; margin-top: 4px"></span>'
                  + (f'<span style="flex-grow: 1; width: 2px; background: {LINE}"></span>' if i < len(steps) - 1 else "") + '</div>'
                  f'<div style="flex-grow: 1; padding-bottom: 16px; display: flex; flex-direction: column; gap: 3px">'
                  f'<div style="display: flex; align-items: center; gap: 8px"><span style="font-family: {MONO}; font-size: 13px; font-weight: 500">{e}</span>{chip("differs from run B", WARN, "#3A2C14") if diff else ""}</div>'
                  f'<span style="font-size: 12.5px; color: {MUTED}">{p}</span>'
                  f'<span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{v + " · " if v else ""}cites {c}</span></div></div>')
    approx = table(["Approximation", "Value"], [["van der Waals", "cut-off 12 Å + tail correction"], ["Electrostatics", "SPME · relative tolerance 1e−5"], ["Constraints", "LINCS · bonds to H"],
                                                 ["Precision", "mixed, deterministic"], ["Estimated parameters", "none"]], ["45%", "55%"], mono_cols=(1,), fs=12)
    hashes = table(["Input", "SHA-256"], [["PS_atactic_DP40.caps", "9f2c…a41e"], ["gaff2.toml", "71be…09c2"], ["recipe.yaml", "c3d0…5e17"]], ["55%", "45%"], mono_cols=(0, 1), fs=12)
    content = (pagehead("Provenance · PS_cell_20x40.data", "Everything that produced this file, in order, with parameters, seeds, citations and approximations",
                        row(select_btn("Compare with", "run B · PS_cell_20x40_b.data"), btn("Export BibTeX", ic="download"), btn("View JSON", ic="file"), gap=8))
               + '<div style="flex-grow: 1; display: flex; gap: 14px; padding: 18px 22px; min-height: 0">'
               + card("Steps", f'<div style="display: flex; flex-direction: column">{"".join(tl)}</div>', chip("6 steps · 1 differs", WARN, "#3A2C14"), extra="flex-grow: 1")
               + f'<div style="width: 480px; flex-shrink: 0; display: flex; flex-direction: column; gap: 14px">'
               + card("Difference from run B", col(kv("grow.cbmc seed", f'20260923 <span style="color: {DIM}">vs</span> 20260924'), kv("Everything else", "identical", False),
                                                   f'<div style="font-size: 12px; color: {MUTED}; line-height: 1.5">Same inputs, versions and parameters; only the seed changed, so differences between the two cells are statistical.</div>', gap=8))
               + card("Approximations", approx) + card("Inputs", hashes) + '</div></div>')
    return shell("Jobs", "CAPS — provenance viewer", f'{icon("chevr", 12, DIM)}<span>Provenance</span>', content,
                 ("<span>caps-manifest/1.0 · deterministic run</span>", "<span>CAPS 0.1.0 · 8c1e2f4</span>"), ["PS_cell_20x40.data"])


def select_btn(label, value):
    return (f'<button style="height: 36px; padding: 0 12px; display: flex; align-items: center; gap: 8px; background: {BG2}; border: 1px solid {LINE}; border-radius: 6px; font-size: 12.5px; color: {TEXT}; cursor: pointer">'
            f'<span style="color: {DIM}">{label}</span><span>{value}</span>{icon("chev", 13, DIM)}</button>')
