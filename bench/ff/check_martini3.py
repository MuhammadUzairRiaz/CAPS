#!/usr/bin/env python3
"""Martini 3 in CAPS against GROMACS: boxes of Martini 3's own molecules (built from the force field's molecule templates,
packed, typed and recognised again from the packed file), exported by CAPS (ff apply --gromacs) and evaluated by
GROMACS (grompp, mdrun -rerun): every energy term and every force compared.

Martini 3's run settings are the model's: reaction-field Coulomb (epsilon_r 15, epsilon_rf infinity), Lennard-Jones
shifted to zero at 1.1 nm, every pair from the table (no mixing), bonds 1-2 excluded only. The cases cover the ions,
solvents, lipids, the nucleobases (virtual sites: the weighted centre, exclusions), restricted-bending angles (TXE) and
harmonic impropers; the last case builds every molecule template once.

GROMACS is mixed precision: energies agree to ~1e-5 relative, forces to a few 1e-3 kcal/mol/A (plus ~1e-4 of the force
on beads held by constraints, which CAPS keeps as stiff bonds).
usage: check_martini3.py [--only substring] [--keep DIR]
"""
import json, math, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
GMX = os.environ.get("GMX", shutil.which("gmx") or "gmx")
FF = os.path.join(ROOT, "data", "forcefields", "martini3.json")
MOLS = json.load(open(os.path.join(ROOT, "data", "martini", "martini3-molecules.json")))["molecules"]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
only = arg("--only", "")
work = arg("--keep", "") or tempfile.mkdtemp()
os.makedirs(work, exist_ok=True)
KJ = 4.184

# (label, box edge in A, [(template, count)])
CASES = [
    ("Water and ions (W, NA, CL, CA, reaction field)", 40, [("W", 500), ("NA", 10), ("CL", 14), ("CA", 2)]),
    ("Lipids in water (POPC, POPS, DOPE, ions)", 55, [("POPC", 6), ("POPS", 4), ("DOPE", 4), ("NA", 4), ("W", 900)]),
    ("Nucleobases and TXE (virtual sites, exclusions, restricted bending, impropers)", 45,
     [("CYTO", 3), ("GUAN", 3), ("THYM", 3), ("URAC", 3), ("MIND", 3), ("TXE", 3), ("PCRE", 3), ("BIM", 3), ("W", 600)]),
    ("Every Martini 3 molecule template once", 90, [(k, 1) for k in sorted(MOLS)] + [("W", 1200)]),
]
# all-atom proteins mapped as martinize2 maps them (vermouth's Martini 3 tests; the files' own cells)
VREF = os.path.expanduser("~/vermouth-ref/tests-m3/tier-1")
PROTEINS = [("Ubiquitin 1UBQ, all-atom -> Martini 3 (restricted-bending backbone, side-chain fix dihedrals)", os.path.join(VREF, "1UBQ", "aa.pdb")),
            ("Lysozyme 3LZT, all-atom -> Martini 3 (tryptophan virtual sites, exclusions, impropers)", os.path.join(VREF, "lysozyme", "aa.pdb")),
            ("Lysozyme as martinize2 wrote it (cg.pdb + topol.top read by CAPS; elastic network 500)",
             (os.path.join(VREF, "lysozyme", "martinize2", "cg.pdb"), os.path.join(VREF, "lysozyme", "martinize2", "topol.top")))]


def gmx(d, stem):
    with open(os.path.join(d, stem + ".mdp"), "a") as f:
        f.write("nstfout = 1\n")
    run = lambda args, inp=None: subprocess.run([GMX, "-quiet"] + args, cwd=d, input=inp, capture_output=True, text=True)
    r = run(["grompp", "-f", stem + ".mdp", "-c", stem + ".gro", "-p", stem + ".top", "-o", stem + ".tpr", "-maxwarn", "10"])
    if r.returncode:
        raise RuntimeError("grompp: " + " | ".join([l for l in (r.stdout + r.stderr).splitlines() if l.strip()][-4:])[:300])
    r = run(["mdrun", "-s", stem + ".tpr", "-rerun", stem + ".gro", "-deffnm", stem, "-nt", "4"])
    if r.returncode:
        raise RuntimeError("mdrun: " + (r.stdout + r.stderr).strip().splitlines()[-1][:200])
    terms = ["Bond", "G96Angle", "Angle", "Restr.-Angles", "Proper-Dih.", "Improper-Dih.", "Per.-Imp.-Dih.", "LJ-(SR)",
             "Coulomb-(SR)", "Potential"]
    r = run(["energy", "-f", stem + ".edr", "-o", "e.xvg"], "\n".join(terms) + "\n\n")
    e = {}
    for l in r.stdout.split("Energy", 1)[-1].splitlines():
        m = re.match(r"^(\S.*?)\s+(-?[\d.]+(?:e[-+]?\d+)?)\s+--", l)
        if m:
            e[m.group(1).replace("-", " ")] = float(m.group(2))
    run(["traj", "-s", stem + ".tpr", "-f", stem + ".trr", "-of", "f.xvg", "-xvg", "none"], "0\n")
    vals = [float(x) for l in open(os.path.join(d, "f.xvg")) if l.strip() for x in l.split()[1:]]
    return e, [tuple(vals[3 * i:3 * i + 3]) for i in range(len(vals) // 3)]


def stiff_bonds(d, edge_nm):
    """Per atom of the exported system: the force uncertainty (kcal/mol/A) from its stiff bonds in single precision."""
    mols, cur, sec = {}, None, None
    for f in ("case.itp", "case.top"):
        if not os.path.exists(os.path.join(d, f)):
            continue
        for raw in open(os.path.join(d, f)):
            l = raw.split(";")[0].strip()
            m = re.match(r"\[\s*(\S+)\s*\]", l)
            if m:
                sec = m.group(1)
                continue
            w = l.split()
            if not w:
                continue
            if sec == "moleculetype":
                cur = mols.setdefault(w[0], {"n": 0, "k": {}})
            elif sec == "atoms" and cur is not None:
                cur["n"] += 1
                cur["q"] = cur.get("q", 0.0) + float(w[6])
            elif sec == "bonds" and cur is not None and len(w) >= 5 and float(w[4]) >= 1e5:
                for a in (int(w[0]), int(w[1])):
                    cur["k"][a] = cur["k"].get(a, 0) + float(w[4])
            elif sec == "molecules":
                for _ in range(int(w[1])):
                    order.append(w[0])
    out = []
    charge[0] = sum(mols[name].get("q", 0.0) for name in order)
    for name in order:
        mt = mols[name]
        out += [mt["k"].get(a, 0) * 6e-8 * edge_nm / 41.84 for a in range(1, mt["n"] + 1)]
    return out


order, charge = [], [0.0]
rows, fails = [], 0
for label, pdb in PROTEINS:
    if (only and only.lower() not in label.lower()) or not os.path.exists(pdb if isinstance(pdb, str) else pdb[0]):
        continue
    CASES.append((label, None, pdb))
for label, edge, content in CASES:
    if only and only.lower() not in label.lower():
        continue
    d = os.path.join(work, re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")[:40])
    os.makedirs(d, exist_ok=True)
    if edge is None:   # a protein: mapped and typed straight from the all-atom file
        content_file, content = content, []
    inp = [f"tolerance 4.3", "output box.data", f"pbc 0 0 0 {edge} {edge} {edge}", "seed 7"]
    bad = None
    for k, (t, n) in enumerate(content):
        r = subprocess.run([CAPS, "build", "--template", t, "--ff", FF, "-o", os.path.join(d, f"m{k}.data"), "--seed", "3"],
                           capture_output=True, text=True)
        if r.returncode:
            bad = f"{t}: {r.stderr.strip()[:200]}"
            break
        inp += [f"structure m{k}.data", f"  number {n}", f"  inside box 0 0 0 {edge} {edge} {edge}", "end structure"]
    if bad:
        rows.append((label, "CAPS failed to build " + bad, ""))
        fails += 1
        continue
    if edge is None:
        src = content_file
    else:
        open(os.path.join(d, "box.inp"), "w").write("\n".join(inp) + "\n")
        r = subprocess.run([CAPS, "pack", "box.inp", "-o", "box.data", "--quiet"], cwd=d, capture_output=True, text=True)
        if r.returncode:
            rows.append((label, "packing failed: " + r.stderr.strip()[:200], ""))
            fails += 1
            continue
        src = "box.data"
    topo = []
    if isinstance(src, tuple):   # coordinates and a GROMACS topology
        src, topo = src[0], ["--topology", src[1]]
    r = subprocess.run([CAPS, "ff", "apply", src, "--ff", FF, "--gromacs", "case", "--forces", "caps_f.txt"] + topo, cwd=d,
                       capture_output=True, text=True)
    m = re.search(r"energy \(kcal/mol\): bond (\S+)  angle (\S+)  dihedral (\S+)  improper (\S+)  vdW (\S+)  Coulomb (\S+)", r.stdout)
    if r.returncode or not m:
        rows.append((label, "CAPS failed: " + (r.stderr.strip() or r.stdout.strip())[-200:], ""))
        fails += 1
        continue
    recog = re.search(r"(\d+) molecules recognised", r.stdout) or re.search(r"onto (\d+) Martini 3 protein beads", r.stdout)
    ce = dict(zip(["bond", "angle", "dihedral", "improper", "vdw", "coulomb"], map(float, m.groups())))
    cf = [tuple(map(float, w[1:4])) for w in (l.split() for l in open(os.path.join(d, "caps_f.txt"))) if len(w) == 4]
    try:
        ge, gf = gmx(d, "case")
    except Exception as ex:
        rows.append((label, "GROMACS failed: " + str(ex), ""))
        fails += 1
        continue
    G = lambda *k: sum(ge.get(x, 0.0) for x in k) / KJ
    gm = {"bond": G("Bond"), "angle": G("Angle", "G96Angle", "Restr. Angles"), "dihedral": G("Proper Dih."),
          "improper": G("Improper Dih.", "Per. Imp. Dih."), "vdw": G("LJ (SR)"), "coulomb": G("Coulomb (SR)")}
    keys = list(gm)
    de = max(abs(ce[k] - gm[k]) / max(1.0, abs(ce[k])) for k in keys)
    # GROMACS writes the forces after spreading the virtual sites' onto their atoms; a site's own row is not a force on
    # anything: compare real atoms only (CAPS's site rows are zero)
    # Constraints are stiff bonds (1e6 kJ/mol/nm^2): a single-precision position (6e-8 relative, 5e-7 nm at 9 nm) moves
    # such a bond's force by ~0.5 kJ/mol/nm = 0.012 kcal/mol/A whatever its stretch, so a bead is allowed that per stiff bond
    order.clear()
    stiff = stiff_bonds(d, (edge or 80) / 10)
    pairs = [(a, tuple(x / 41.84 for x in b), stiff[i]) for i, (a, b) in enumerate(zip(cf, gf)) if any(a)] if len(cf) == len(gf) else []
    df = max((math.dist(a, b) for a, b, _ in pairs), default=float("inf"))
    # plus 1e-4 of the force itself: a mapped protein, not relaxed, has beads under 1e5 kcal/mol/A from stretched stiff
    # bonds and near-straight restricted-bending angles, which single precision carries to ~4e-5
    worst = max((math.dist(a, b) / (5e-3 + t + 1e-4 * math.hypot(*a)) for a, b, t in pairs), default=float("inf"))
    # the charges are the molecules' own (ions of the same bead type, Na+ / Cl-, told apart)
    want_q = sum(n * sum(a["charge"] for a in MOLS[t]["atoms"]) for t, n in content) if edge is not None else charge[0]
    dq = abs(charge[0] - want_q)
    ok = de < 5e-5 and worst < 1 and dq < 1e-6
    fails += 0 if ok else 1
    natoms = len(cf)
    rows.append((label, f"{'ok' if ok else 'DIFFERS'} · {natoms} beads, {(recog.group(1) + (' molecules recognised' if edge is not None else ' beads mapped')) if recog else 'topology read from ' + os.path.basename(topo[1])} · "
                 f"charge {charge[0]:+.4f} e (molecules {want_q:+.4f}) · energy terms {de:.1e} (relative) · forces {df:.1e} kcal/mol/A",
                 "CAPS/GROMACS (kcal/mol): " + " ".join(f"{k} {ce[k]:.4f}/{gm[k]:.4f}" for k in keys)))

for label, res, more in rows:
    print(f"{label}\n   {res}")
    if more:
        print(f"   {more}")
print(f"\n{sum(1 for r in rows if r[1].startswith('ok'))} of {len(rows)} Martini 3 cases agree with GROMACS"
      f"{'' if not fails else f'; {fails} failed'} (work: {work})")
sys.exit(1 if fails else 0)
