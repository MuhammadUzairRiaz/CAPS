#!/usr/bin/env python3
"""GROMACS files written by CAPS reproduce CAPS's energy and forces in GROMACS.

For each case: CAPS assigns the force field (caps ff apply --gromacs), writes the topology, coordinates and run parameters
and its own energy terms and forces; GROMACS evaluates the same configuration (gmx grompp, mdrun -rerun, gmx energy,
gmx traj -of) and every term is compared: bonds, angles, torsions (proper and improper), Lennard-Jones (pairs, 1-4 pairs
and the tail) and, for periodic cells with CAPS's PME, the Coulomb sum (real space, 1-4, reciprocal).

Two differences are by definition, not error, and are reported apart:
  · the dispersion correction: CAPS (as LAMMPS pair_modify tail) sums the tail over all ordered type pairs,
    E = (2π/V) Σ_ab N_a N_b 4ε_ab [σ¹²/(9rc⁹) − σ⁶/(3rc³)]; GROMACS (DispCorr AllEnerPres) takes the average C6 and C12
    over non-excluded atom pairs. The Lennard-Jones terms are compared without their tails; GROMACS's tail is checked
    against its own formula evaluated from the written topology;
  · without a cell, CAPS's damped shifted force Coulomb has no GROMACS counterpart (plain cut-off there): Coulomb is
    shown, not compared.
  · without tail corrections CAPS shifts every Lennard-Jones term to zero at the cut-off, 1-4 pairs included; GROMACS
    shifts pair-list terms only: the constant Σ over 1-4 pairs is added back (forces are the same).
GROMACS here is mixed precision: energies agree to ~1e-5 relative, forces to a few 1e-3 kcal/mol/Å (single-precision
positions).

usage: check_gromacs.py [--only substring] [--keep DIR]
Needs GMX (default gmx on PATH).
"""
import math, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import template_pdb

CAPS = os.path.join(ROOT, "build", "cli", "caps")
GMX = os.environ.get("GMX", shutil.which("gmx") or "gmx")
LIB = os.environ.get("DLFIELD", os.path.expanduser("~/project/dl_f_4.13")) + "/lib"
FF = os.path.join(ROOT, "data", "forcefields")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
only = arg("--only", "")
work = arg("--keep", "") or tempfile.mkdtemp()
os.makedirs(work, exist_ok=True)
PS = os.path.join(ROOT, "samples", "ps_melt.data")
PME = ["--pme", "--ewald-rtol", "1e-7", "--pme-spacing", "0.5", "--pme-order", "6"]

# (label, structure source, force field, charges, extra caps options)
CASES = [
    ("Polystyrene melt, GAFF2, PME (periodic, 1300 atoms)", ("file", PS), "gaff-amber25", "gasteiger", PME),
    ("Polystyrene melt, GAFF2, PME, shifted LJ (no tail)", ("file", PS), "gaff-amber25", "gasteiger", PME + ["--no-tail"]),
    ("Polystyrene melt, OPLS-AA, PME", ("file", PS), "opls2005", "gasteiger", PME),
    ("Polystyrene melt, OPLS-AA 2024 (its own charges), PME", ("file", PS), "oplsaa2024-moltemplate", "types", PME),
    ("Polystyrene melt, CVFF, PME (cvff impropers)", ("file", PS), "cvff", "types", PME),
    ("Polystyrene melt, GAFF2, PME, triclinic cell (xy 0.5, xz 0.3, yz 0.2 Å)", ("tilt", PS), "gaff-amber25", "gasteiger", PME),
    ("(6,6) nanotube, GAFF2, PME (bonds cross the cell: periodic-molecules)", ("nano", ["tube", "--n", "6", "--m", "6", "--length", "25"]),
     "gaff-amber25", "gasteiger", PME),
    ("GAFF toluene (vacuum)", ("template", "AMBER16_gaff", "toluene"), "gaff-amber16", "gasteiger", []),
    ("CVFF phenol (vacuum, cvff impropers)", ("template", "CVFF", "phenol"), "cvff", "types", []),
    ("OPLS-AA methyl vinyl ketone (vacuum)", ("template", "OPLS2005", "methyl_vinyl_ketone"), "opls2005", "gasteiger", []),
    ("CGenFF toluene (vacuum; separate 1-4 LJ, Urey-Bradley, harmonic impropers)", ("template", "CHARMM36_cgenff", "toluene"), "cgenff", "gasteiger", []),
    # GROMOS in its own settings: reaction field (ε_rf 61) at 1.4 nm, no dispersion correction, C6/C12 (comb-rule 1)
    ("SPC/E water box, GROMOS 54A7 (reaction field, periodic)", ("solvate", "40", "SPC/E", "1200"), "gromos-54a7", "keep", []),
    # ClayFF (LAMMPS's ClayFF test structure, its types and charges): the M-O-H bends are angles without an M-O bond, so
    # GROMACS keeps M and H interacting as CAPS does. IFF 1.5 (CVFF) montmorillonite from IFF's own model database. Both
    # cells are thinner than 2 x 10 A: compared at 8 and 9 A cut-offs.
    ("ClayFF pyrophyllite (M-O-H bends by contact; PME)", ("file", os.path.expanduser("~/lammps/tools/msi2lmp/test/PyAC_bulk-clayff.car")),
     "inorganic-clay", "keep", PME + ["--cutoff", "8", "--no-tail", "--names"]),
    ("IFF 1.5 (CVFF) Na-montmorillonite (bonded clay layers, cvff impropers; PME)",
     ("file", os.path.expanduser("~/iff-ref/INTERFACE_FF_1_5/MODEL_DATABASE/CLAY_MINERALS/mont0_333_Na_15_cell.car")),
     "iff-cvff", "keep", PME + ["--cutoff", "9", "--no-tail"]),
    ("PCFF polystyrene (class II: refused)", ("file", PS), "pcff", "types", []),
    ("UFF polystyrene (inversions: refused)", ("file", PS), "uff", "types", []),
]


def ff_file(fid):
    if fid == "uff":
        return "uff"
    p = os.path.join(FF, fid + ".json")
    return p if os.path.exists(p) else None


def structure(src, base):
    if src[0] == "file":
        return src[1] if os.path.exists(src[1]) else None
    if src[0] == "tilt":   # the same melt in a sheared cell (positions kept, so bonds are recomputed images)
        out = os.path.join(work, base + ".data")
        xy, xz, yz = 0.5, 0.3, 0.2
        lines, sec, lo, L = open(src[1]).read().split("\n"), "", {}, {}
        for k, l in enumerate(lines):
            m = re.match(r"\s*(\S+)\s+(\S+)\s+([xyz])lo [xyz]hi", l)
            if m:
                lo[m.group(3)], L[m.group(3)] = float(m.group(1)), float(m.group(2)) - float(m.group(1))
                if m.group(3) == "z":
                    lines[k] = l + f"\n{xy} {xz} {yz} xy xz yz"
                continue
            if re.match(r"^[A-Z][A-Za-z ]*$", l.split("#")[0].strip() or "-"):
                sec = l.split("#")[0].strip()
                continue
            w = l.split()
            if sec == "Atoms" and len(w) >= 7:   # full: id mol type q x y z [ix iy iz]; shear with the cell
                x, y, z = map(float, w[4:7])
                if len(w) >= 10:   # unwrapped positions, so the images shear with them
                    x += int(w[7]) * L["x"] + int(w[8]) * xy + int(w[9]) * xz
                    y += int(w[8]) * L["y"] + int(w[9]) * yz
                    z += int(w[9]) * L["z"]
                fy, fz = (y - lo["y"]) / L["y"], (z - lo["z"]) / L["z"]
                w[4:7] = [f"{x + xy * fy + xz * fz:.8f}", f"{y + yz * fz:.8f}", f"{z:.8f}"]
                lines[k] = " ".join(w[:7])
        open(out, "w").write("\n".join(lines))
        return out
    if src[0] == "solvate":   # a box of water from CAPS's solvent packing, the model's charges kept
        out = os.path.join(work, base + ".mol2")
        subprocess.run([CAPS, "solvate", "-o", out, "--edge", src[1], "--solvent", "water", "--model", src[2], "--no-ions",
                        "--molecules", src[3], "--tolerance", "2.6"],
                       capture_output=True, check=True)
        return out
    if src[0] == "nano":
        out = os.path.join(work, base + ".data")
        subprocess.run([CAPS, "nano"] + src[1] + ["-o", out], capture_output=True, check=True)
        return out
    sf = os.path.join(LIB, src[1] + ".sf")
    if not os.path.exists(sf):
        return None
    name = src[2]
    if name not in template_pdb.parse_sf(sf)[2]:
        names = template_pdb.candidates(sf, 30)
        if not names:
            return None
        name = names[0]
    pdb = os.path.join(work, base + ".pdb")
    template_pdb.build(sf, name, pdb)
    return pdb


def gmx(d, stem):
    """GROMACS energy terms (kJ/mol) and forces (kJ/mol/nm) of stem.gro."""
    with open(os.path.join(d, stem + ".mdp"), "a") as f:
        f.write("nstfout = 1\n")
    run = lambda args, inp=None: subprocess.run([GMX, "-quiet"] + args, cwd=d, input=inp, capture_output=True, text=True)
    r = run(["grompp", "-f", stem + ".mdp", "-c", stem + ".gro", "-p", stem + ".top", "-o", stem + ".tpr", "-maxwarn", "10"])
    if r.returncode:
        raise RuntimeError("grompp: " + [l for l in (r.stdout + r.stderr).splitlines() if l.strip()][-3][:200])
    r = run(["mdrun", "-s", stem + ".tpr", "-rerun", stem + ".gro", "-deffnm", stem, "-nt", "4"])
    if r.returncode:
        raise RuntimeError("mdrun: " + (r.stdout + r.stderr).strip().splitlines()[-1][:200])
    terms = ["Bond", "Morse", "Quartic-Bonds", "Angle", "G96Angle", "U-B", "Proper-Dih.", "Per.-Imp.-Dih.", "Improper-Dih.", "LJ-14",
             "Coulomb-14", "LJ-(SR)", "Disper.-corr.", "Coulomb-(SR)", "Coul.-recip.", "Potential"]
    # names absent from this run are reported by gmx energy and skipped; the table lists the rest
    r = run(["energy", "-f", stem + ".edr", "-o", os.path.join(d, "e.xvg")], "\n".join(terms) + "\n\n")
    e = {}
    for l in r.stdout.split("Energy", 1)[-1].splitlines():
        m = re.match(r"^(\S.*?)\s+(-?[\d.]+(?:e[-+]?\d+)?)\s+--", l)
        if m:
            e[m.group(1).replace("-", " ")] = float(m.group(2))
    r = run(["traj", "-s", stem + ".tpr", "-f", stem + ".trr", "-of", "f.xvg", "-xvg", "none"], "0\n")
    vals = [float(x) for l in open(os.path.join(d, "f.xvg")) if l.strip() for x in l.split()[1:]]
    f = [tuple(vals[3 * i:3 * i + 3]) for i in range(len(vals) // 3)]
    return e, f


def expand_topology(d, stem):
    """(pair table, per-atom type names, global exclusion pairs) of the whole system: the .top with its #include'd
    .itp, each molecule type repeated as [ molecules ] lists it."""
    text = []
    for l in open(os.path.join(d, stem + ".top")):
        m = re.match(r'\s*#include\s+"([^"]+)"', l)
        if m:
            text += open(os.path.join(d, m.group(1))).read().splitlines()
        else:
            text.append(l.rstrip("\n"))
    pairs, moltypes, cur, sec, molecules = {}, {}, None, "", []
    for raw in text:
        l = raw.split(";")[0].strip()
        if not l or l.startswith("#"):
            continue
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        w = l.split()
        if sec == "nonbond_params":
            pairs[(w[0], w[1])] = pairs[(w[1], w[0])] = (float(w[3]) * 10, float(w[4]) / 4.184)
        elif sec == "moleculetype":
            cur = w[0]
            moltypes[cur] = {"atoms": [], "excl": []}
        elif sec == "atoms" and cur:
            moltypes[cur]["atoms"].append(w[1])
        elif sec == "exclusions" and cur:
            moltypes[cur]["excl"] += [(int(w[0]) - 1, int(j) - 1) for j in w[1:]]
        elif sec == "molecules":
            molecules.append((w[0], int(w[1])))
    at, excl = [], []
    for name, count in molecules:
        mt = moltypes[name]
        for _ in range(count):
            off = len(at)
            at += mt["atoms"]
            excl += [(i + off, j + off) for i, j in mt["excl"]]
    return pairs, at, excl


def top_tail(d, stem, cutoff):
    """CAPS's tail term (kcal/mol) from the topology's pair coefficients and the box: (2π/V) Σ_ab N_a N_b 4ε[σ¹²/9rc⁹ − σ⁶/3rc³]."""
    pairs, at, _ = expand_topology(d, stem)
    counts = {}
    for a in at:
        counts[a] = counts.get(a, 0) + 1
    box = [float(x) * 10 for x in open(os.path.join(d, stem + ".gro")).read().strip().splitlines()[-1].split()]
    vol = box[0] * box[1] * box[2]
    rc = cutoff
    e = 0.0
    for a, na in counts.items():
        for b, nb in counts.items():
            s, eps = pairs[(a, b)]
            e += na * nb * 2 * math.pi * 4 * eps * (s ** 12 / (9 * rc ** 9) - s ** 6 / (3 * rc ** 3))
    return e / vol


def gmx_tail(d, stem, cutoff):
    """GROMACS's DispCorr AllEnerPres (kcal/mol): (2π N²/V)[⟨C12⟩/(9rc⁹) − ⟨C6⟩/(3rc³)], the averages over all atom
    pairs less the excluded ones."""
    P, at, excl = expand_topology(d, stem)
    c6 = lambda a, b: 4 * P[(a, b)][1] * P[(a, b)][0] ** 6
    c12 = lambda a, b: 4 * P[(a, b)][1] * P[(a, b)][0] ** 12
    n = len(at)
    cnt = {}
    for a in at:
        cnt[a] = cnt.get(a, 0) + 1
    s6 = s12 = 0.0
    for a, na in cnt.items():
        for b, nb in cnt.items():
            k = 0.5 * (na * nb if a != b else na * (na - 1))
            s6 += k * c6(a, b)
            s12 += k * c12(a, b)
    for i, j in excl:
        s6 -= c6(at[i], at[j])
        s12 -= c12(at[i], at[j])
    npair = n * (n - 1) / 2 - len(excl)
    box = [float(x) * 10 for x in open(os.path.join(d, stem + ".gro")).read().strip().splitlines()[-1].split()]
    vol = box[0] * box[1] * box[2]
    return 2 * math.pi * n * n / vol * (s12 / npair / (9 * cutoff ** 9) - s6 / npair / (3 * cutoff ** 3))


def pair14_shift(d, stem, cutoff):
    """Σ over [pairs] of 4ε'[(σ/rc)¹² − (σ/rc)⁶] (kcal/mol), each molecule type times its count: CAPS shifts its 1-4
    terms at the cut-off when there is no tail correction, GROMACS never shifts [pairs]; a constant, with no force."""
    text = []
    for l in open(os.path.join(d, stem + ".top")):
        m = re.match(r'\s*#include\s+"([^"]+)"', l)
        text += open(os.path.join(d, m.group(1))).read().splitlines() if m else [l]
    per, counts, cur, sec = {}, {}, None, ""
    for raw in text:
        l = raw.split(";")[0].strip()
        if not l or l.startswith("#"):
            continue
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        w = l.split()
        if sec == "moleculetype":
            cur = w[0]
            per[cur] = 0.0
        elif sec == "pairs" and cur:
            s, eps = float(w[3]) * 10, float(w[4]) / 4.184
            per[cur] += 4 * eps * ((s / cutoff) ** 12 - (s / cutoff) ** 6)
        elif sec == "molecules":
            counts[w[0]] = counts.get(w[0], 0) + int(w[1])
    return sum(per[k] * c for k, c in counts.items())

rows, fails = [], 0
KJ = 4.184
for label, src, fid, charges, extra in CASES:
    if only and only.lower() not in label.lower():
        continue
    base = re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")[:40]
    ffj, s = ff_file(fid), structure(src, base)
    if not ffj or not s:
        rows.append((label, "skipped (source not available)", ""))
        continue
    d = os.path.join(work, base)
    os.makedirs(d, exist_ok=True)
    cmd = [CAPS, "ff", "apply", s, "--ff", ffj, "--charges", charges, "--gromacs", os.path.join(d, "case"),
           "--forces", os.path.join(d, "caps_f.txt")] + extra
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode or "missing parameters" in r.stdout:
        why = (r.stderr.strip().splitlines() or r.stdout.strip().splitlines() or ["?"])[-1]
        expected = "refused" in label and "no GROMACS" in why
        rows.append((label, ("refused, as intended: " if expected else "CAPS failed: ") + why[:160], ""))
        fails += 0 if expected else 1
        continue
    m = re.search(r"energy \(kcal/mol\): bond (\S+)  angle (\S+)  dihedral (\S+)  improper (\S+)  vdW (\S+)  Coulomb (\S+)", r.stdout)
    ce = dict(zip(["bond", "angle", "dihedral", "improper", "vdw", "coulomb"], map(float, m.groups())))
    cf = [tuple(map(float, w[1:4])) for w in (l.split() for l in open(os.path.join(d, "caps_f.txt"))) if len(w) == 4]
    try:
        ge, gf = gmx(d, "case")
    except Exception as ex:
        rows.append((label, "GROMACS failed: " + str(ex), ""))
        fails += 1
        continue
    G = lambda *k: sum(ge.get(x, 0.0) for x in k) / KJ
    # a real reaction field (GROMOS's ε_rf 61), not the plain cut-off GROMACS writes as Reaction-Field with ε_rf 1 (vacuum)
    mdp = {l.split("=")[0].strip(): l.split("=")[1].split(";")[0].strip() for l in open(os.path.join(d, "case.mdp")) if "=" in l}
    rf = mdp.get("coulombtype") == "Reaction-Field" and float(mdp.get("epsilon-rf", "1")) != 1
    periodic = "Disper. corr." in ge or "Coul. recip." in ge or rf
    # the cut-off both sides used: the .mdp's rvdw (the force field's own, e.g. OPLS 12 Å)
    rc = next((float(l.split("=")[1].split(";")[0]) * 10 for l in open(os.path.join(d, "case.mdp")) if l.split("=")[0].strip() == "rvdw"), 10.0)
    tail_caps = top_tail(d, "case", rc) if ("--no-tail" not in extra and "Disper. corr." in ge) else 0.0
    gm = {"bond": G("Bond", "Morse", "Quartic Bonds"),
          "angle": G("Angle", "G96Angle", "U B"), "dihedral": G("Proper Dih."), "improper": G("Per. Imp. Dih.", "Improper Dih."),
          "vdw": G("LJ 14", "LJ (SR)"), "coulomb": G("Coulomb 14", "Coulomb (SR)", "Coul. recip.")}
    # GROMACS books a periodic improper (function 4) as "Per. Imp. Dih."; CAPS may count function-9 impropers as torsions
    cm = dict(ce)
    cm["vdw"] = ce["vdw"] - tail_caps
    shift14 = pair14_shift(d, "case", rc) if "--no-tail" in extra else 0.0
    cm["vdw"] += shift14
    keys = ["bond", "angle", "dihedral", "improper", "vdw"] + (["coulomb"] if periodic else [])
    de = max(abs(cm[k] - gm[k]) / max(1.0, abs(cm[k])) for k in keys)
    df = max(math.dist(a, tuple(x / 41.84 for x in b)) for a, b in zip(cf, gf)) if len(cf) == len(gf) else float("inf")
    # forces include the tail (none: it is a constant) and, without a cell, Coulomb of different methods: compare only with a cell
    # mixed precision: energies to 5e-5 relative (single-precision PME); positions are stored in single precision
    # (3.6e-7 nm at 3 nm), which on a C–H bond (2.8e5 kJ/mol/nm²) is 0.1 kJ/mol/nm = 2.4e-3 kcal/mol/Å of force
    dt = abs(G("Disper. corr.") - gmx_tail(d, "case", rc)) / max(1.0, abs(G("Disper. corr."))) if tail_caps else 0.0
    # the force tolerance from that single-precision position error on the stiffest bond of the case: 2 K δx, K the largest
    # harmonic constant (kJ/mol/nm², the .itp's function 1) and δx = 6e-8 of the box edge (SPC water's 4500 kcal/mol/Å²
    # bonds in a 4 nm box: 0.02 kcal/mol/Å)
    kmax = 0.0
    for l in open(os.path.join(d, "case.itp")) if os.path.exists(os.path.join(d, "case.itp")) else []:
        w = l.split(";")[0].split()
        if len(w) == 5 and w[2] == "1":
            try:
                kmax = max(kmax, float(w[4]))
            except ValueError:
                pass
    box = max(float(x) for x in open(os.path.join(d, "case.gro")).read().strip().splitlines()[-1].split()[:3])
    ftol = max(5e-3, 2 * kmax * 6e-8 * box / 41.84 * 2)
    ok = de < 5e-5 and dt < 5e-5 and (df < ftol or not periodic)
    fails += 0 if ok else 1
    extra_txt = []
    if periodic and tail_caps:
        extra_txt.append(f"tail CAPS {tail_caps:.4f} / GROMACS {G('Disper. corr.'):.4f} kcal/mol (definitions differ; GROMACS's own formula "
                         f"reproduced to {dt:.1e})")
    if shift14:
        extra_txt.append(f"CAPS's 1-4 cut-off shift {shift14:.4f} kcal/mol added back (GROMACS pairs are unshifted)")
    if not periodic:
        extra_txt.append(f"Coulomb CAPS DSF {ce['coulomb']:.4f} / GROMACS cut-off {gm['coulomb']:.4f} kcal/mol (methods differ; forces not compared)")
    rows.append((label, f"{'ok' if ok else 'DIFFERS'} · energy terms {de:.1e} (relative)" + (f" · forces {df:.1e} kcal/mol/Å" if periodic else ""),
                 " · ".join(extra_txt) + ("\n   CAPS/GROMACS: " + " ".join(f"{k} {cm[k]:.4f}/{gm[k]:.4f}" for k in keys) if not ok else "")))

for label, res, more in rows:
    print(f"{label}\n   {res}")
    if more:
        print(f"   {more}")
print(f"\n{sum(1 for r in rows if r[1].startswith('ok'))} of {len(rows)} cases agree with GROMACS; "
      f"{sum(1 for r in rows if 'as intended' in r[1])} refused as intended"
      f"{'' if not fails else f'; {fails} failed'} (work: {work})")
sys.exit(1 if fails else 0)
