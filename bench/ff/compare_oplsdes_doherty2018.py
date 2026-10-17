#!/usr/bin/env python3
"""OPLS-DES as CAPS built it from the paper (oplsdes-doherty2018) against a reference it was not built from:
Acevedo's own GROMACS files (github.com/orlandoacevedo/DES, GPL 3; a local copy). Nothing of them is copied into CAPS:
this script reads them where they are and reports agreement.

For each of the paper's molecules, on the reference's own coordinates nudged off their ideal geometry (0.08 A, fixed seed):
  Acevedo's files  CAPS types the molecule and writes a GROMACS topology; per atom the charge and Lennard-Jones terms,
                   per bond, bend and torsion (by atom indices) the parameters are compared with the reference .itp
                   (Fourier coefficients kcal/mol; impropers by their centre), and GROMACS evaluates both topologies
                   on the same coordinates (mdrun -rerun): every energy term.

usage: compare_oplsdes_doherty2018.py DES_DIR   (DES_DIR: the DES repository's DES folder, with ITP and PDBs)
"""
import math, os, random, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
FF = os.path.join(ROOT, "data", "forcefields", "oplsdes-doherty2018.json")
GMX = shutil.which("gmx") or "gmx"
DES = sys.argv[1]
work = tempfile.mkdtemp()
KJ = 4.184
MOLS = [  # name, reference .itp stem, reference .pdb
    ("choline", "choline", "choline"), ("chloride", "cl", "chloride"), ("ethylene glycol", "ethgly", "ethylene_glycol"),
    ("glycerol", "glycerol", "glycerol"), ("levulinic acid", "levulinic", "levulinic"), ("malonic acid", "malonic", "malonic"),
    ("oxalic acid", "oxalic", "oxalic_acid"), ("phenol", "phenol", "phenol"), ("urea", "urea", "urea")]


def sections(path):
    out, cur = {}, None
    for line in open(path):
        s = line.split(";")[0].strip()
        if not s or s.startswith("#"):
            continue
        m = re.fullmatch(r"\[\s*(\w+)\s*\]", s)
        if m:
            cur = m.group(1)
            out.setdefault(cur, [])
        elif cur:
            out[cur].append(s.split())
    return out


def read_pdb(path):
    atoms = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            atoms.append([l[12:16].strip(), float(l[30:38]), float(l[38:46]), float(l[46:54]), l[76:78].strip() or l[12:16].strip()[0]])
    return atoms


def write_pdb(path, atoms, resname="MOL"):
    with open(path, "w") as f:
        for i, (n, x, y, z, el) in enumerate(atoms, 1):
            f.write(f"HETATM{i:5d} {n:<4s} {resname:>3s} A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}\n")
        f.write("END\n")


def bonded_of(sec):
    """bonds, angles, torsions and impropers of a topology, keyed by atom indices (torsions as V1..V4 kcal/mol)"""
    bonds = {tuple(sorted(map(int, w[:2]))): (float(w[3]), float(w[4]) / 2 / KJ / 100) for w in sec.get("bonds", [])}
    angles = {}
    for w in sec.get("angles", []):
        i, j, k = map(int, w[:3])
        angles[(min(i, k), j, max(i, k))] = (float(w[4]), float(w[5]) / 2 / KJ)
    bset = set(bonds)
    tors, imps = {}, {}
    for w in sec.get("dihedrals", []):
        q = tuple(map(int, w[:4]))
        f = int(w[4])
        chain = all(tuple(sorted(q[a:a + 2])) in bset for a in range(3))
        if f == 5:   # Fourier: V1..V4 in kJ/mol
            v = [float(x) / KJ for x in w[5:9]]
        elif f in (9, 4, 1):   # k (1 + cos(n phi - phi_s))
            phs, k, n = float(w[5]), float(w[6]), int(w[7])
            v = [0.0] * 4
            sign = 1 if (n % 2 == 1) == (abs(phs) < 1) else -1
            v[n - 1] = sign * 2 * k / KJ
        else:
            raise ValueError(f"function {f}")
        if chain:
            key = q if q[1] < q[2] else q[::-1]
            old = tors.get(key, [0.0] * 4)
            tors[key] = [a + b for a, b in zip(old, v)]
        else:   # improper: centre = the atom bonded to the other three
            c = next(a for a in q if sum(tuple(sorted((a, b))) in bset for b in q if b != a) == 3)
            key = (c, tuple(sorted(a for a in q if a != c)))
            imps.setdefault(c, []).append((key, v))
    return bonds, angles, tors, imps


def gmx_energy(d, top, gro, mdp, stem):
    subprocess.run([GMX, "grompp", "-f", mdp, "-c", gro, "-p", top, "-o", stem + ".tpr", "-maxwarn", "10"], cwd=d, capture_output=True, check=True)
    subprocess.run([GMX, "mdrun", "-s", stem + ".tpr", "-rerun", gro, "-deffnm", stem, "-nt", "1"], cwd=d, capture_output=True, check=True)
    log = open(os.path.join(d, stem + ".log")).read()
    blk = log.split("Energies (kJ/mol)")[1].split("\n\n")[0].splitlines()[1:]
    out = {}
    for names, vals in zip(blk[0::2], blk[1::2]):
        ns = [names[i:i + 15].strip() for i in range(0, len(names), 15)]
        vs = [float(vals[i:i + 15]) for i in range(0, len(vals), 15)]
        out.update(zip(ns, vs))
    return out


def group(e):
    g = lambda *k: sum(e.get(x, 0.0) for x in k) / KJ
    return {"bond": g("Bond"), "angle": g("Angle"), "torsion+improper": g("Proper Dih.", "Ryckaert-Bell.", "Per. Imp. Dih.", "Improper Dih.", "Fourier Dih."),
            "LJ": g("LJ-14", "LJ (SR)"), "Coulomb": g("Coulomb-14", "Coulomb (SR)")}


def cmpd(label, c, r, tol, zero=None, who="ref"):
    """terms keyed by atom indices, CAPS's against the reference's; zero: a term one side leaves out counts as this
    (CAPS writes no all-zero torsion)"""
    keys = set(c) | set(r)
    if zero is not None:
        c = {k: c.get(k, zero) for k in keys}
        r = {k: r.get(k, zero) for k in keys}
    same = sum(1 for k in keys if k in c and k in r and all(abs(a - b) < tol for a, b in zip(c[k], r[k])))
    diff = [k for k in keys if k in c and k in r and not all(abs(a - b) < tol for a, b in zip(c[k], r[k]))]
    only_c, only_r = [k for k in keys if k not in r], [k for k in keys if k not in c]
    s = f"{label} {same}/{len(keys)}"
    if diff:
        s += " (differ: " + "; ".join(f"{k} CAPS {[round(x, 4) for x in c[k]]} {who} {[round(x, 4) for x in r[k]]}" for k in sorted(diff)[:3]) \
             + (" ..." if len(diff) > 3 else "") + ")"
    if only_c:
        s += f" (CAPS only: {len(only_c)})"
    if only_r:
        s += f" ({who} only: {len(only_r)})"
    return s


def compare_acevedo(name, stem, pdbname):
    d = os.path.join(work, stem)
    os.makedirs(d, exist_ok=True)
    atoms = read_pdb(os.path.join(DES, "PDBs", pdbname + ".pdb"))
    rnd = random.Random(7)
    nudged = [[n, x + rnd.uniform(-0.08, 0.08), y + rnd.uniform(-0.08, 0.08), z + rnd.uniform(-0.08, 0.08), el] for n, x, y, z, el in atoms]
    write_pdb(os.path.join(d, "mol.pdb"), nudged)
    r = subprocess.run([CAPS, "ff", "apply", os.path.join(d, "mol.pdb"), "--ff", FF, "--charges", "types", "--gromacs", os.path.join(d, "caps")],
                       capture_output=True, text=True)
    if r.returncode:
        return f"CAPS failed: {(r.stderr or r.stdout).strip().splitlines()[-1]}"
    cs, rs = sections(os.path.join(d, "caps.itp")), sections(os.path.join(DES, "ITP", stem + "_DES.itp"))
    ctop, rtyp = sections(os.path.join(d, "caps.top")), sections(os.path.join(DES, "ITP", stem + "_atomtypes_DES.itp"))
    cat = {w[0]: (float(w[5]) * 10, float(w[6]) / KJ) for w in ctop["atomtypes"]}
    rat = {w[0]: (float(w[4]) * 10, float(w[5]) / KJ) for w in rtyp["atomtypes"]}
    notes = []
    names, nq, nlj, ljd = [], 0, 0, set()
    for c, rr in zip(cs["atoms"], rs["atoms"]):
        names.append((c[1], rr[1]))
        nq += abs(float(c[6]) - float(rr[6])) < 1e-6
        same = all(abs(a - b) < 1e-5 for a, b in zip(cat[c[1]], rat[rr[1]]))
        nlj += same
        if not same:
            ljd.add(f"{rr[1]} sigma, eps CAPS {cat[c[1]][0]:.4f} {cat[c[1]][1]:.5f} ref {rat[rr[1]][0]:.4f} {rat[rr[1]][1]:.5f}")
    n = len(rs["atoms"])
    same_name = sum(c.split("_")[0] == r_ for c, r_ in names)
    mism = sorted({f"{c}~{r_}" for c, r_ in names if c.split("_")[0] != r_})
    out = [f"atoms {n}: type names {same_name}/{n}" + (f" ({', '.join(mism)})" if mism else "") + f", charges {nq}/{n}, LJ {nlj}/{n}" + (f" ({'; '.join(sorted(ljd))})" if ljd else "")]
    cb, ca, ct, ci = bonded_of(cs)
    rb, ra, rt, ri = bonded_of(rs)

    out.append(cmpd("bonds", cb, rb, 2e-3) + ", " + cmpd("bends", ca, ra, 2e-3) + ", " + cmpd("torsions", ct, rt, 2e-3, [0.0] * 4))
    fmt = lambda imp, at: {at[c - 1][1]: [round(v[1], 4) for _, v in t] for c, t in imp.items()}
    if ci or ri:
        out.append(f"impropers (centre: V2 of each, kcal/mol): CAPS {fmt(ci, cs['atoms'])} reference {fmt(ri, rs['atoms'])}")
    # energies: both topologies in GROMACS on CAPS's coordinates (the same atom order as the reference's .pdb)
    shutil.copy(os.path.join(DES, "ITP", stem + "_DES.itp"), os.path.join(d, "ref.itp"))
    shutil.copy(os.path.join(DES, "ITP", stem + "_atomtypes_DES.itp"), os.path.join(d, "ref_types.itp"))
    # the reference again with CAPS's 1-4 pair list added (its .itp has no [ pairs ]: gen-pairs from its own types, x 0.5)
    txt = open(os.path.join(d, "ref.itp")).read()
    pairs = "".join(f"{w[0]} {w[1]} 1\n" for w in cs.get("pairs", []))
    open(os.path.join(d, "refp.itp"), "w").write(re.sub(r"(\[\s*bonds\s*\])", "[ pairs ]\n" + pairs + "\n\\1", txt, count=1) if pairs else txt)
    molname = rs["moleculetype"][0][0]
    for top, itp in (("ref.top", "ref.itp"), ("refp.top", "refp.itp")):
        with open(os.path.join(d, top), "w") as f:
            f.write(f'[ defaults ]\n1 3 yes 0.5 0.5\n#include "ref_types.itp"\n#include "{itp}"\n[ system ]\nref\n[ molecules ]\n{molname} 1\n')
    # an 8 nm box (CAPS's own is the molecule plus 2 r_c, too small for GROMACS's pair list at r_c 1.6 nm for a small molecule)
    gro = open(os.path.join(d, "caps.gro")).read().splitlines()
    xyz = [tuple(map(float, l[20:].split()[:3])) for l in gro[2:-1]]
    cen = [sum(v[i] for v in xyz) / len(xyz) for i in range(3)]
    gro[2:-1] = [l[:20] + "".join(f"{v[i] - cen[i] + 4.0:10.5f}" for i in range(3)) for l, v in zip(gro[2:-1], xyz)]
    gro[-1] = "   8.00000   8.00000   8.00000"
    open(os.path.join(d, "box.gro"), "w").write("\n".join(gro) + "\n")
    if n == 1:   # a lone ion: nothing intramolecular (and grompp takes no lone charge in a cut-off box)
        return "\n      ".join(out)
    try:
        ec = group(gmx_energy(d, "caps.top", "box.gro", "caps.mdp", "ec"))
        er = group(gmx_energy(d, "ref.top", "box.gro", "caps.mdp", "er"))
        ep = group(gmx_energy(d, "refp.top", "box.gro", "caps.mdp", "ep"))
    except subprocess.CalledProcessError as e:
        out.append("GROMACS failed: " + (e.stderr or b"").decode()[-200:])
        return "\n      ".join(out)
    out.append("energy CAPS / reference / reference with 1-4 pairs (kcal/mol): "
               + "  ".join(f"{k} {ec[k]:.4f}/{er[k]:.4f}/{ep[k]:.4f}" for k in ec)
               + f"  total {sum(ec.values()):.4f}/{sum(er.values()):.4f}/{sum(ep.values()):.4f}")
    return "\n      ".join(out)


print(f"OPLS-DES (CAPS, from the paper) against Acevedo's GROMACS files in {DES}")
for name, stem, pdbname in MOLS:
    print(f"  {name}\n      " + compare_acevedo(name, stem, pdbname))
print(f"(work: {work})")
