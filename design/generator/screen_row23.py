"""Row 23 — export & automate: data-file export, Python step, saved pipelines, batch, topology distributions, vectors."""
import math, re, hashlib
from lib import *
from screen_app import card
from display import ps_cell
from screen_row20 import annotate, draw, _shell, CHAIN, VIRIDIS, ramp, LTYPE
from screen_row21 import cell_recs, mi, inspector, viewport, frame, CW, VH, BOTTOM, L, MASS, NA
from screen_row22 import mol_shape, chain_col, H3

MONO_BLOCK = f'font-family: {MONO}; font-size: 11.5px; line-height: 1.6'


def topology(cell, recs):
    """Global bonds, angles and proper dihedrals from the construction-order bonds."""
    bonds, angles, dihs = [], [], []
    base = 0
    for atoms, bl, _ in cell:
        nb = [[] for _ in atoms]
        for a, b in bl:
            nb[a].append(b); nb[b].append(a)
            bonds.append((base + a, base + b))
        for j in range(len(atoms)):
            for x in range(len(nb[j])):
                for y in range(x + 1, len(nb[j])):
                    angles.append((base + nb[j][x], base + j, base + nb[j][y]))
        for a, b in bl:
            for i in nb[a]:
                if i == b: continue
                for l in nb[b]:
                    if l == a or l == i: continue
                    dihs.append((base + i, base + a, base + b, base + l))
        base += len(atoms)
    return bonds, angles, dihs


def canon(t):
    return min(t, t[::-1])


def type_sets(recs, bonds, angles, dihs):
    g = [r["gaff"] for r in recs]
    bt = sorted({canon(tuple(g[i] for i in b)) for b in bonds})
    at = sorted({canon(tuple(g[i] for i in a)) for a in angles})
    dt = sorted({canon(tuple(g[i] for i in d)) for d in dihs})
    return bt, at, dt


def code(lines, start=1, hl_line=None):
    kw = r"\b(def|return|for|in|import|from|if|else|and|not|None|True|False|with|as|lambda)\b"
    out = []
    for i, ln in enumerate(lines):
        s = esc(ln)
        if "#" in s:
            k = s.index("#")
            body, com = s[:k], f'<span style="color: {DIM}">{s[k:]}</span>'
        else:
            body, com = s, ""
        body = re.sub(r"(&quot;[^&]*&quot;|\"[^\"]*\"|'[^']*')", lambda m: f'<span style="color: {OK}">{m.group(1)}</span>', body)
        body = re.sub(kw, lambda m: f'<span style="color: #C79BE8">{m.group(1)}</span>', body)
        bg = f"background: {BG3};" if hl_line == i + start else ""
        out.append(f'<div style="display: flex; gap: 14px; {bg}"><span style="width: 26px; flex-shrink: 0; text-align: right; color: {DIM}">{i + start}</span><span style="white-space: pre">{body}{com}</span></div>')
    return f'<div style="{MONO_BLOCK}; overflow: hidden">{"".join(out)}</div>'


def heading(title, sub, actions):
    return (f'<div style="flex-shrink: 0; display: flex; align-items: center; gap: 16px; padding: 14px 22px; background: {BG1}; border-bottom: 1px solid {LINE}">'
            f'<div style="display: flex; flex-direction: column; gap: 3px"><h1 style="margin: 0; font-size: 19px; font-weight: 600">{title}</h1><div style="font-size: 12.5px; color: {MUTED}">{sub}</div></div>'
            f'<div style="flex-grow: 1"></div>{actions}</div>')


# ---------------------------------------------------------------- Export data file
def export_data():
    cell, recs = cell_recs()
    bonds, angles, dihs = topology(cell, recs)
    bt, at, dt = type_sets(recs, bonds, angles, dihs)
    n = len(recs)
    img = lambda p: tuple(math.floor(v / L) for v in p)
    lines = ["CAPS 0.1 · PS_melt · frame 120 · atom_style full · units real", "",
             f"{n} atoms", f"{len(bonds)} bonds", f"{len(angles)} angles", f"{len(dihs)} dihedrals", "",
             f"4 atom types", f"{len(bt)} bond types", f"{len(at)} angle types", f"{len(dt)} dihedral types", "",
             f"0.0 {L:.1f} xlo xhi", f"0.0 {L:.1f} ylo yhi", f"0.0 {L:.1f} zlo zhi", "", "Masses", ""]
    lines += [f"{v} {12.011 if k[0] == 'c' else 1.008}  # {k}" for k, v in LTYPE.items()]
    lines += ["", "Atoms  # full", ""]
    for r in recs[:3]:
        w = [v % L for v in r["p"]]
        ix = img(r["p"])
        lines.append(f"{r['id']} {r['mol']} {r['type']} {r['q']:.6f} {w[0]:.5f} {w[1]:.5f} {w[2]:.5f} {ix[0]} {ix[1]} {ix[2]}")
    lines += ["…", "", "Bonds", ""]
    btid = {t: i + 1 for i, t in enumerate(bt)}
    g = [r["gaff"] for r in recs]
    for k, (a, b) in enumerate(bonds[:3]):
        lines.append(f"{k + 1} {btid[canon((g[a], g[b]))]} {a + 1} {b + 1}")
    lines.append("…")
    pv = "".join(f'<div style="display: flex; gap: 12px"><span style="width: 24px; flex-shrink: 0; text-align: right; color: {DIM}">{i + 1}</span><span style="white-space: pre; color: {ACC if l in ("Masses", "Atoms  # full", "Bonds") else (DIM if l.startswith("CAPS") or "#" in l else TEXT)}">{esc(l)}</span></div>' for i, l in enumerate(lines))
    est = sum(len(x) + 1 for x in lines[:-1])  # header part only; full size below is an estimate
    atom_line = len(f"1300 10 4 -0.061234 12.34567 12.34567 12.34567 0 0 0") + 1
    size_kb = (atom_line * n + 16 * len(bonds) + 20 * len(angles) + 24 * len(dihs) + 600) / 1024
    fmts = [("LAMMPS data", "atom_style full", True), ("LAMMPS dump", "custom columns", False), ("GROMACS .gro + .top", "nm, residues per unit", False),
            ("PDB · mmCIF", "CONECT, residues", False), ("XYZ · extended XYZ", "lattice in comment line", False), ("mol2 · SDF", "bond orders, charges", False), ("CAPS .caps", "everything + provenance", False)]
    fl = "".join(f'<div style="display: flex; align-items: center; gap: 10px; padding: 7px 10px; border-radius: 6px; {"background: " + BG3 + "; outline: 1px solid " + ACC if on else ""}"><span style="width: 150px; flex-shrink: 0; font-size: 12.5px; font-weight: 600">{a}</span><span style="font-size: 11.5px; color: {MUTED}">{b}</span></div>' for a, b, on in fmts)
    tt = table(["Section", "Count", "Types"], [["Atoms", str(n), "4"], ["Bonds", str(len(bonds)), str(len(bt))], ["Angles", str(len(angles)), str(len(at))], ["Dihedrals", str(len(dihs)), str(len(dt))]],
               ["44%", "30%", "26%"], mono_cols=(1, 2), align_right=(1, 2), fs=12, rowh=26)
    left = (f'<div style="width: 330px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Format", f'<div style="display: flex; flex-direction: column; gap: 2px">{fl}</div>', "", 12)
            + card("Options", col(select("Coordinates", "wrapped + image flags"), select("Charges", "Gasteiger (from Field)"), toggle("Write type labels (# c3)", True), toggle("Include Coeffs sections", False), gap=8), "", 12)
            + '</div>')
    mid = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
           + card("PS_melt.data · preview", f'<div style="{MONO_BLOCK}; line-height: 1.52; overflow: hidden">{pv}</div>', chip("real values of this frame", OK, BG2), 14, "flex-grow: 1; min-height: 0") + '</div>')
    right = (f'<div style="width: 300px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Topology written", col(tt, f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Angles and dihedrals are generated from the bonds. Types come from GAFF names, so c3–c3 and c3–ca are separate bond types.</span>', gap=8), "", 12)
             + card("File", col(kv("Name", "PS_melt.data"), kv("Size (est.)", f"{size_kb:.0f} KB"), kv("Coeffs", "none · use Field output"), gap=6), "", 12)
             + '</div>')
    content = (heading("Export · data file", "One frame, with topology, written for the engine you will run next", btn("Cancel") + btn("Export", True, "download"))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 14px 22px; min-height: 0; overflow: hidden">{left}{mid}{right}</div>')
    return _shell("CAPS — export data file", f'{icon("chevr", 12, DIM)}<span>Export › Data file</span>', content,
                  (f"<span>{n} atoms · {len(bonds)} bonds · {len(angles)} angles · {len(dihs)} dihedrals</span>", "<span>units: LAMMPS real</span>"))


# ---------------------------------------------------------------- Python step
def backbone_dihedrals(cell):
    out = []
    for m, (atoms, bl, bb) in enumerate(cell):
        for i in range(len(bb) - 3):
            out.append((m, i, dihedral(*(atoms[bb[i + k]]["p"] for k in range(4)))))
    return out


def python_step():
    cell, recs = cell_recs()
    dh = backbone_dihedrals(cell)
    per = {}
    for m, i, phi in dh:
        per.setdefault(m + 1, []).append(abs(phi) > 120)
    ree = {}
    for m, (atoms, bl, bb) in enumerate(cell):
        ree[m + 1] = dist(atoms[bb[0]]["p"], atoms[bb[-1]]["p"])
    ft = sum(sum(v) for v in per.values()) / sum(len(v) for v in per.values())
    src = ["import numpy as np", "import caps", "from caps.pipeline import step", "", "@step(name=\"Backbone conformation\")", "def modify(frame, data):",
           "    mol = data.particles[\"Molecule Identifier\"]", "    bb = data.particles[\"Backbone\"]           # 1 on backbone carbons", "    pos = data.particles.positions_unwrapped",
           "    ree, trans = {}, []", "    for m in np.unique(mol):", "        idx = np.flatnonzero((mol == m) & (bb == 1))", "        ree[m] = np.linalg.norm(pos[idx[-1]] - pos[idx[0]])",
           "        phi = caps.geometry.dihedrals(pos[idx])", "        trans += list(np.abs(phi) > 120)", "    data.attributes[\"TransFraction\"] = np.mean(trans)",
           "    data.tables[\"Ree\"] = ree                   # shows in Data inspector"]
    rows = [[str(m), f"{ree[m]:.2f}", f"{sum(per[m])}/{len(per[m])}"] for m in sorted(ree)]
    t = table(["mol", "Ree (Å)", "trans"], rows, ["26%", "38%", "36%"], mono_cols=(0, 1, 2), align_right=(0, 1, 2), fs=11.5, rowh=21)
    console = [">>> run step on frame 120", f"TransFraction = {ft:.3f}   ({sum(sum(v) for v in per.values())} of {len(dh)} backbone dihedrals)",
               f"Ree: {len(ree)} rows written to data.tables", "done · cached for this frame"]
    cons = "".join(f'<div style="white-space: pre; color: {DIM if l.startswith(">>>") else TEXT}">{esc(l)}</div>' for l in console)
    ed = card("backbone_conformation.py", code(src, hl_line=16), row(chip("Python 3.12 · numpy", MUTED, BG2), chip("in-process", OK, BG2), gap=6), 14, "flex-grow: 1; min-height: 0")
    co = card("Console", f'<div style="{MONO_BLOCK}">{cons}</div>', "", 12)
    left = f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">{ed}{co}</div>'
    right = (f'<div style="width: 330px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Output · Ree table", t, chip(f"{len(rows)} rows", MUTED, BG2, True), 12)
             + card("Step", col(kv("Inputs", "Molecule Identifier, Backbone"), kv("Writes", "TransFraction, Ree"), toggle("Re-run on frame change", True), gap=6), "", 12) + '</div>')
    content = (heading("Python pipeline step", "A function in the pipeline: reads properties, writes new ones, cached per frame", btn("Open in notebook") + btn("Run", True, "play"))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 14px 22px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return _shell("CAPS — Python step", f'{icon("chevr", 12, DIM)}<span>Analyze › Visualize › Python step</span>', content,
                  (f"<span>TransFraction = {ft:.3f} on this frame</span>", "<span>the same step runs from caps run and in notebooks</span>"))


# ---------------------------------------------------------------- Save & reuse pipeline
def save_pipeline():
    yaml = ["caps_pipeline: 1", "name: PS melt · structure report", "source:", "  file: PS_melt.lammpstrj", "  topology: PS_melt.data", "  columns: {id: Particle Identifier, mol: Molecule Identifier, type: Particle Type, q: Charge}",
            "steps:", "  - unwrap: {using: image_flags}", "  - create_bonds: {mode: topology, check_cutoffs: {C-C: 1.70, C-H: 1.25}}",
            "  - compute_property: {name: DistanceToCOM, expr: \"norm(Position - MoleculeCOM(MoleculeIdentifier))\"}",
            "  - cluster_analysis: {cutoff: 3.3, unit: molecule, periodic: true}", "  - rdf: {pairs: [C, C], rmax: 12.0, bin: 0.2, other_molecules: true}",
            "  - python: {file: backbone_conformation.py}", "  - colour_by: {property: DistanceToCOM, map: viridis}",
            "outputs:", "  - table: Ree -> ree.csv", "  - plot: rdf -> rdf.svg", "  - render: {size: [1920, 1080], overlays: [label, legend, scale_bar]}"]
    text = "\n".join(yaml) + "\n"
    h = hashlib.sha256(text.encode()).hexdigest()
    ed = card("ps_structure_report.caps-pipeline.yaml", code(yaml), chip(f"sha256 {h[:12]}…", MUTED, BG2, True), 14, "flex-grow: 1; min-height: 0")
    steps = [("Unwrap", "image flags"), ("Create bonds", "topology"), ("Compute property", "DistanceToCOM"), ("Cluster analysis", "3.3 Å"), ("Coordination & RDF", "C–C"), ("Python", "backbone_conformation.py"), ("Colour by", "viridis")]
    st = "".join(f'<div style="display: flex; align-items: center; gap: 10px; padding: 6px 0; border-bottom: 1px solid {BG2}"><span style="width: 20px; font-family: {MONO}; font-size: 11px; color: {DIM}">{i + 1}</span><span style="flex-grow: 1; font-size: 12.5px">{a}</span><span style="font-family: {MONO}; font-size: 11px; color: {MUTED}">{b}</span></div>' for i, (a, b) in enumerate(steps))
    right = (f'<div style="width: 380px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Steps in this pipeline", f'<div>{st}</div>', chip(f"{len(steps)}", MUTED, BG2, True), 12)
             + card("Save", col(field("Name", "PS melt · structure report", mono=False), select("Scope", "this project · shared library"), toggle("Store with every result it makes", True), toggle("Ask before running a Python file", True), gap=8), "", 12)
             + card("Run it elsewhere", f'<div style="{MONO_BLOCK}">caps run ps_structure_report.caps-pipeline.yaml \\<br>&nbsp;&nbsp;--input runs/*/PS_melt.lammpstrj</div>', "", 12)
             + '</div>')
    content = (heading("Save pipeline", "Plain text you can diff, review and version; results record its hash", btn("Cancel") + btn("Save", True, "save"))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 14px 22px; min-height: 0; overflow: hidden">{ed}{right}</div>')
    return _shell("CAPS — save pipeline", f'{icon("chevr", 12, DIM)}<span>Analyze › Visualize › Save pipeline</span>', content,
                  (f"<span>{len(text.encode())} bytes · sha256 {h[:16]}</span>", "<span>the hash is of the text shown</span>"))


# ---------------------------------------------------------------- Batch
def quick_metrics(seed):
    cell = ps_cell(10, 8, L, seed=seed)
    recs = annotate(cell)
    sh = mol_shape(recs)
    heavy = {}
    for r in recs:
        if r["e"] != "H": heavy.setdefault(r["mol"], []).append(r["p"])
    mols = sorted(heavy)
    parent = {m: m for m in mols}
    def f(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, a in enumerate(mols):
        for b in mols[i + 1:]:
            if f(a) != f(b) and any(mi(p, q) < 3.3 for p in heavy[a] for q in heavy[b]):
                parent[f(a)] = f(b)
    ncl = len({f(m) for m in mols})
    dh = backbone_dihedrals(cell)
    ft = sum(abs(p) > 120 for _, _, p in dh) / len(dh)
    out = sum(1 for r in recs if any(not (0 <= v < L) for v in r["p"]))
    return dict(seed=seed, n=len(recs), rg=sum(s["rg"] for s in sh.values()) / len(sh), k2=sum(s["k2"] for s in sh.values()) / len(sh), ncl=ncl, ft=ft, out=out)


def batch():
    seeds = [21, 22, 23, 25, 26, 27]
    res = []
    for sd in seeds:
        try:
            res.append(quick_metrics(sd))
        except RuntimeError as e:
            res.append(dict(seed=sd, err=str(e)))
    running = len(seeds) - 1
    rows = []
    for i, r in enumerate(res):
        if "err" in r:
            state, vals = chip("failed", ERR, BG2), ["—"] * 5
        elif i == running:
            state, vals = chip("running", ACC, BG2), ["—"] * 5
        else:
            state, vals = chip("done", OK, BG2), [f"{r['rg']:.2f}", f"{r['k2']:.3f}", f"{r['ft']:.3f}", str(r["ncl"]), str(r["out"])]
        rows.append([f"cells/seed{r['seed']}/PS_melt.data", state] + vals + ["[measure]" if state.find("done") > 0 else "—"])
    t = table(["Input", "State", "mean Rg (Å)", "mean κ²", "trans", "clusters", "atoms outside", "time"], rows,
              ["26%", "10%", "11%", "10%", "9%", "9%", "13%", "12%"], mono_cols=(0, 2, 3, 4, 5, 6, 7), align_right=(2, 3, 4, 5, 6), fs=12, rowh=30)
    fails = [r for r in res if "err" in r]
    errbox = "".join(f'<div style="display: flex; gap: 10px; align-items: flex-start; padding: 8px 10px; border-radius: 6px; background: {BG2}; border: 1px solid {ERR}"><span style="flex-shrink: 0">{icon("alert", 16, ERR)}</span><span style="font-size: 12px; line-height: 1.5"><b>seed{r["seed"]}</b> · {esc(r["err"])}. The batch skipped it and continued; the row stays in the results as failed.</span></div>' for r in fails)
    done = [r for i, r in enumerate(res) if "err" not in r and i != running]
    d = done
    mean = lambda k: sum(x[k] for x in d) / len(d)
    sd = lambda k: math.sqrt(sum((x[k] - mean(k)) ** 2 for x in d) / (len(d) - 1))
    summ = table(["Metric", "mean", f"s.d. (n = {len(d)})"], [["mean Rg (Å)", f"{mean('rg'):.2f}", f"{sd('rg'):.2f}"], ["mean κ²", f"{mean('k2'):.3f}", f"{sd('k2'):.3f}"], ["trans fraction", f"{mean('ft'):.3f}", f"{sd('ft'):.3f}"]],
                 ["44%", "28%", "28%"], mono_cols=(1, 2), align_right=(1, 2), fs=12, rowh=28)
    pts = [(seeds.index(x["seed"]) + 1, x["rg"]) for x in d]
    lo, hi = math.floor(min(p[1] for p in pts) * 2) / 2 - 0.5, math.ceil(max(p[1] for p in pts) * 2) / 2 + 0.5
    pl = plot(360, 170, [(pts, ACC, 1.5, "3 3")], (0.5, len(seeds) + 0.5), (lo, hi), list(range(1, len(seeds) + 1)), [lo, (lo + hi) / 2, hi], "input #", "mean Rg (Å)",
              markers=[(a, b, ACC) for a, b in pts])
    top = card("Batch · ps_structure_report", col(t, errbox, row(field("Inputs", "cells/seed*/PS_melt.data", mono=True), select("Workers", "4 local"), select("On error", "skip and continue"), gap=10), gap=12), row(chip(f"{len(d)} of {len(seeds)} done · {len(fails)} failed", ACC, BG2), gap=6), 14)
    bottom = (f'<div style="display: flex; gap: 14px">' + card("Across finished inputs", summ, "", 12, "flex-grow: 1")
              + card("mean Rg per input", pl, "", 12, "width: 420px; flex-shrink: 0")
              + card("Where results go", col(kv("Tables", "batch/results.csv"), kv("Per input", "batch/&lt;input&gt;/…"), kv("Pipeline", "sha256 recorded per row"), gap=6), "", 12, "width: 320px; flex-shrink: 0") + '</div>')
    content = (heading("Batch over files", "One saved pipeline, many inputs, one table out", btn("Pause") + btn("Open results", True, "folder"))
               + f'<div style="flex-grow: 1; display: flex; flex-direction: column; gap: 14px; padding: 14px 22px; min-height: 0; overflow: hidden">{top}{bottom}</div>')
    return _shell("CAPS — batch", f'{icon("chevr", 12, DIM)}<span>Analyze › Batch</span>', content,
                  ("<span>values: built cells (10 × PS DP 8, 33 Å), one per seed</span>", "<span>times are [measure] until CAPS runs</span>"))


# ---------------------------------------------------------------- Bond / angle / dihedral distributions
def bars_svg(counts, lo, hi, w, h, col_fn, label, xt):
    nb = len(counts)
    mx = max(counts) or 1
    bars = "".join(f'<rect x="{34 + i * (w - 44) / nb:.1f}" y="{h - 22 - c / mx * (h - 38):.1f}" width="{max(1.0, (w - 44) / nb - 1.5):.1f}" height="{c / mx * (h - 38):.1f}" fill="{col_fn(lo + (i + 0.5) * (hi - lo) / nb)}"></rect>' for i, c in enumerate(counts) if c)
    ticks = "".join(f'<text x="{34 + (v - lo) / (hi - lo) * (w - 44):.1f}" y="{h - 7}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{v:g}</text>' for v in xt)
    return (f'<svg width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-label="{label}" style="display: block">{bars}{ticks}'
            f'<line x1="34" y1="{h - 22}" x2="{w - 10}" y2="{h - 22}" stroke="{LINE}"></line><text x="34" y="11" font-size="10" fill="{DIM}">{label}</text></svg>')


def hist(vals, lo, hi, nb):
    c = [0] * nb
    for v in vals:
        if lo <= v < hi: c[int((v - lo) / (hi - lo) * nb)] += 1
    return c


def dihedral_state(phi):
    return "t" if abs(phi) > 120 else ("g+" if phi > 0 else "g−")


STATE_COL = {"t": ACC, "g+": SEL, "g−": "#C79BE8"}


def topo_dist():
    cell, recs = cell_recs()
    bonds, angles, dihs = topology(cell, recs)
    P = [r["p"] for r in recs]
    E = [r["e"] for r in recs]
    bl = [dist(P[a], P[b]) for a, b in bonds]
    cc = [x for x, (a, b) in zip(bl, bonds) if E[a] == E[b] == "C"]
    chh = [x for x, (a, b) in zip(bl, bonds) if E[a] != E[b]]
    ccc = [angle(P[a], P[b], P[c]) for a, b, c in angles if E[a] == E[b] == E[c] == "C"]
    dh = backbone_dihedrals(cell)
    phis = [p for _, _, p in dh]
    states = {s: sum(1 for p in phis if dihedral_state(p) == s) for s in ("t", "g+", "g−")}
    w3 = (CW - 6 - 32 - 36) // 3
    hb = bars_svg(hist(bl, 1.0, 1.7, 70), 1.0, 1.7, w3, BOTTOM - 44, lambda v: SEL if v < 1.2 else ACC, "bond length (Å) · all bonds", [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
    ha = bars_svg(hist(ccc, 100, 130, 60), 100, 130, w3, BOTTOM - 44, lambda v: ACC, "C–C–C angle (°)", [100, 110, 120, 130])
    hd = bars_svg(hist(phis, -180, 180, 36), -180, 180, w3, BOTTOM - 44, lambda v: STATE_COL[dihedral_state(v)], "backbone dihedral (°)", [-180, -120, -60, 0, 60, 120, 180])
    sc = Scene("td", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.0)
    for m, (atoms, bl_, bb) in enumerate(cell):
        for i in range(len(bb) - 3):
            s = dihedral_state(dihedral(*(atoms[bb[i + k]]["p"] for k in range(4))))
            sc.add_tube([atoms[bb[i + 1]]["p"], atoms[bb[i + 2]]["p"]], STATE_COL[s], 0.55)
        sc.add_tube([atoms[bb[0]]["p"], atoms[bb[1]]["p"]], "#4A525A", 0.4)
        sc.add_tube([atoms[bb[-2]]["p"], atoms[bb[-1]]["p"]], "#4A525A", 0.4)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    leg = "".join(f'<span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 12px; height: 4px; border-radius: 2px; background: {c}"></span>{k} {states[k]}</span>' for k, c in STATE_COL.items())
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; flex-direction: column; gap: 6px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="display: flex; gap: 18px; font-family: {MONO}; font-size: 11.5px; color: {MUTED}"><span>C–C {min(cc):.3f}–{max(cc):.3f} Å</span><span>C–H {min(chh):.3f}–{max(chh):.3f} Å</span><span>C–C–C {min(ccc):.1f}–{max(ccc):.1f}°</span><span style="flex-grow: 1"></span>{leg}</div>'
              f'<div style="display: flex; gap: 18px">{hb}{ha}{hd}</div></div>')
    right = inspector("Topology distributions", section("Terms", col(seg(["Bonds", "Angles", "Dihedrals"], "Dihedrals", full=True), select("Dihedrals", "backbone C–C–C–C"), select("States", "t |φ|>120°, g± otherwise"), gap=8))
                      + section("Counts", col(kv("Bonds", str(len(bonds))), kv("Angles", str(len(angles))), kv("Dihedrals (all)", str(len(dihs))), kv("Backbone dihedrals", str(len(phis))), gap=6))
                      + f'<div style="padding: 0 18px 12px; font-size: 11.5px; color: {DIM}; line-height: 1.5">A freshly built cell: bond lengths are the builder\'s exact values, so the bond histogram is a set of spikes. The dihedral mix is the builder\'s sampling, not an equilibrium result.</div>')
    return frame(viewport(sc.svg(), "Backbone bonds by dihedral state", "Backbone") + bottom, right, "CAPS — topology distributions", "Analyze › Visualize › Bonds, angles, dihedrals",
                 (f"<span>{len(phis)} backbone dihedrals · t {states['t']} · g+ {states['g+']} · g− {states['g−']}</span>", "<span>built cell, not equilibrated</span>"))


# ---------------------------------------------------------------- Vectors
def vectors():
    cell, recs = cell_recs()
    sc = Scene("vc", CW - 6, VH, yaw=0.55, pitch=0.4, persp=0.25, fog=0.3, atom_k=1.0)
    vecs = []
    for m, (atoms, bl, bb) in enumerate(cell):
        c = CHAIN[m % 10]
        sc.add_tube([atoms[j]["p"] for j in bb], mix(c, BG0, 0.6), 0.4)
        a, b = atoms[bb[0]]["p"], atoms[bb[-1]]["p"]
        v = sub(b, a)
        vecs.append((m + 1, v))
        sc.add_line(a, b, c, 3.2, None, 1.0)
        u = norm(v)
        side = norm(cross(u, (0, 0, 1) if abs(u[2]) < 0.9 else (1, 0, 0)))
        for s in (1, -1):
            sc.add_line(b, add(add(b, mul(u, -1.8)), mul(side, 0.9 * s)), c, 3.2, None, 1.0)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.0, None, 0.7)
    Q = [[0.0] * 3 for _ in range(3)]
    for _, v in vecs:
        u = norm(v)
        for i in range(3):
            for j in range(3):
                Q[i][j] += 1.5 * u[i] * u[j] - (0.5 if i == j else 0.0)
    Q = [[x / len(vecs) for x in r] for r in Q]
    from screen_row22 import eig3
    S, director = eig3(Q)[2]
    pz = sum(1.5 * norm(v)[2] ** 2 - 0.5 for _, v in vecs) / len(vecs)
    rows = [[f'<span style="display: inline-flex; align-items: center; gap: 6px"><span style="width: 10px; height: 10px; border-radius: 2px; background: {CHAIN[(m - 1) % 10]}"></span>{m}</span>', " ".join(f"{x:6.2f}" for x in v), f"{math.sqrt(dot3(v, v)):.2f}", f"{math.degrees(math.acos(abs(norm(v)[2]))):.0f}°"] for m, v in vecs[:6]]
    t = table(["mol", "R = r_last − r_first (Å)", "|R| (Å)", "angle to z"], rows, ["12%", "46%", "20%", "22%"], mono_cols=(1, 2, 3), align_right=(2, 3), fs=11.5, rowh=22)
    bottom = (f'<div style="height: {BOTTOM}px; flex-shrink: 0; display: flex; gap: 18px; padding: 10px 16px; background: {BG1}; border-top: 1px solid {LINE}">'
              f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 6px"><h3 style="{H3}">End-to-end vectors · 6 of {len(vecs)}</h3>{t}</div>'
              f'<div style="width: 260px; flex-shrink: 0; display: flex; flex-direction: column; gap: 8px; padding-top: 18px">{kv("Order S", f"{S:.3f}")}{kv("Director", " ".join(f"{x:.2f}" for x in director))}{kv("P₂ vs z", f"{pz:+.3f}")}'
              f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">From only {len(vecs)} vectors, S is noisy: even random directions give S well above 0 at this n.</span></div></div>')
    right = inspector("Vectors", section("Vector", col(select("Property", "end-to-end R (per molecule)"), select("Anchor", "first backbone carbon"), gap=8))
                      + section("Glyph", col(row(field("Width", "0.30", "Å"), field("Scale", "1.0"), gap=8), select("Colour", "by molecule"), toggle("Arrowheads", True), toggle("Flip direction", False), gap=8))
                      + section("Other sources", col(kv("Forces · velocities", "from dump columns"), kv("Dipoles", "from charges, per unit"), gap=6)))
    return frame(viewport(sc.svg(), "Arrows: end-to-end vectors", "Backbone") + bottom, right, "CAPS — vectors", "Analyze › Visualize › Vectors",
                 (f"<span>{len(vecs)} vectors · S = {S:.3f}</span>", "<span>order tensor from unit vectors</span>"))


BOARDS = (("ExportData", export_data), ("PythonStep", python_step), ("SavePipeline", save_pipeline),
          ("BatchRun", batch), ("TopologyDist", topo_dist), ("Vectors", vectors))

if __name__ == "__main__":
    import os, sys
    os.makedirs("stage23/project", exist_ok=True)
    only = sys.argv[1:]
    for name, fn in BOARDS:
        if only and name not in only: continue
        h = fn(); open(f"stage23/project/{name}.dc.html", "w").write(h); print(name, len(h))
