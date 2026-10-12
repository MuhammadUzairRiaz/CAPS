#!/usr/bin/env python3
"""OPLS-AA atom types from CAPS (OPLS-AA 2024 rules) against hand-built GROMACS OPLS-AA topologies: the organic-liquid
set of Caleman et al., J. Chem. Theory Comput. 8 (2012) 61 (virtualchemistry.org; foyer ships it as opls_validation),
one .top and .gro per molecule with opls_NNN on every atom.

  compare_oplsaa_top.py prepare TOPDIR OUTDIR        (needs RDKit: bond orders from the connectivity and net charge)
  compare_oplsaa_top.py compare OUTDIR [--ff oplsaa2024-moltemplate] [--show 80]

prepare writes OUTDIR/NAME.sdf (bond orders from RDKit's DetermineBondOrders on the topology's bonds and charge) and
OUTDIR/NAME.ref (the opls number of each atom); compare types each sdf with CAPS and compares the numbers atom by atom.
"""
import os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


def read_top(path):
    atoms, bonds, sec = [], [], None
    for raw in open(path):
        l = raw.split(";")[0].strip()
        if not l or l.startswith("#"):
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", l)
        if m:
            sec = m.group(1); continue
        w = l.split()
        if sec == "atoms":
            atoms.append((w[1], w[4], float(w[6] if len(w) >= 8 else w[5])))   # the residue column can be empty
        elif sec in ("bonds", "constraints"):   # the nitriles hold C#N as a constraint
            bonds.append((int(w[0]) - 1, int(w[1]) - 1))
        elif sec == "molecules":
            break
    return atoms, bonds


def read_gro(path):
    L = open(path).read().splitlines()
    n = int(L[1])
    return [tuple(10 * float(L[2 + i][20 + 8 * k:28 + 8 * k]) for k in range(3)) for i in range(n)]


def atomic_numbers(itp):
    """opls_NNN -> atomic number, from the force field's ffnonbonded.itp"""
    z = {}
    for raw in open(itp):
        w = raw.split(";")[0].split()
        if len(w) > 3 and w[0].startswith("opls_") and w[2].isdigit():
            z[w[0]] = int(w[2])
    return z


def prepare(topdir, out):
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from rdkit.Geometry import Point3D
    os.makedirs(out, exist_ok=True)
    Z = atomic_numbers(os.path.join(topdir, "oplsaa.ff", "ffnonbonded.itp"))
    for name in sorted(os.listdir(topdir)):
        tp = os.path.join(topdir, name, name + ".top")
        gp = os.path.join(topdir, name, name + ".gro")
        if not (os.path.exists(tp) and os.path.exists(gp)):
            continue
        atoms, bonds = read_top(tp)
        xyz = read_gro(gp)
        if len(xyz) != len(atoms):
            print(f"{name}: {len(atoms)} atoms in top, {len(xyz)} in gro, skipped"); continue
        if any(t not in Z for t, *_ in atoms):   # massless virtual sites (the nitriles' ML): left out, with their bonds
            keep = [i for i, a in enumerate(atoms) if a[0] in Z]
            new = {i: k for k, i in enumerate(keep)}
            print(f"{name}: {len(atoms) - len(keep)} virtual sites left out")
            atoms = [atoms[i] for i in keep]
            xyz = [xyz[i] for i in keep]
            bonds = [(new[i], new[j]) for i, j in bonds if i in new and j in new]
        rw = Chem.RWMol()
        for t, nm, q in atoms:
            rw.AddAtom(Chem.Atom(Z[t]))
        for i, j in bonds:
            rw.AddBond(i, j, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, p in enumerate(xyz):
            conf.SetAtomPosition(i, Point3D(*p))
        mol.AddConformer(conf)
        charge = round(sum(a[2] for a in atoms))
        try:
            bonded = {i for b in bonds for i in b}
            if mol.GetNumAtoms() > 1 and len(bonded) < mol.GetNumAtoms():
                # atoms placed as virtual sites (the nitrile C and N) have no bonds in the topology: all bonds from the geometry
                mol = Chem.RWMol(mol)
                for b in list(mol.GetBonds()):
                    mol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                mol = mol.GetMol()
                rdDetermineBonds.DetermineConnectivity(mol)
            if mol.GetNumAtoms() == 1:   # a monatomic ion
                mol.GetAtomWithIdx(0).SetFormalCharge(charge); mol.GetAtomWithIdx(0).SetNoImplicit(True)
            else:
                rdDetermineBonds.DetermineBondOrders(mol, charge=charge)
            Chem.SanitizeMol(mol)
            Chem.Kekulize(mol, clearAromaticFlags=True)
        except Exception as e:
            print(f"{name}: no bond orders ({e})"); continue
        Chem.MolToMolFile(mol, os.path.join(out, name + ".sdf"), kekulize=True)
        open(os.path.join(out, name + ".ref"), "w").write("".join(t.replace("opls_", "") + "\n" for t, *_ in atoms))


# GROMACS numbers whose meaning differs from the 2024 table's, with the 2024 number of the same atom: 145B is the
# biphenyl bridge C (class C!, the table's 519); 925-931 are another alkyne set there, 927 the alkyne C carrying an alkyl
# (the table's 925) and 930 that alkyl's H (HC +0.06, the table's 140). In the 2024 table 927 and 930 are sp3 carbons.
GROMACS_TO_2024 = {"145B": "519", "927": "925", "930": "140"}
# the one deliberate difference: epoxides keep the table's epoxide types (three-ring angle terms); the topologies' ether
# types give methyloxirane's ring 109.5 deg angles, about 119 kcal/mol of strain against 44
EXPLAINED = {}
EXPLAINED_MOLECULE = {"methyloxirane": "epoxide: the table's epoxide types 1025-1029 with three-ring angle terms; the topologies' "
                                       "ether types put 109.5 deg angles on the ring (about 119 kcal/mol strain, 44 with these)"}


def compare(out, ff, show):
    caps = os.path.join(ROOT, "build", "cli", "caps")
    desc = {}
    import json
    for r in json.load(open(os.path.join(ROOT, "data", "typing", ff + ".typing.json")))["rules"]:
        desc.setdefault(r["type"], r.get("description", ""))
    ok = tot = untyped = full = n = 0
    diffs = []
    for f in sorted(os.listdir(out)):
        if not f.endswith(".sdf"):
            continue
        name = f[:-4]
        ref = [GROMACS_TO_2024.get(t, t) for t in open(os.path.join(out, name + ".ref")).read().split()]
        tf = os.path.join(out, name + ".caps")
        subprocess.run([caps, "ff", "type", os.path.join(out, f), "--ff", ff, "-o", tf], capture_output=True, text=True)
        n += 1
        if not os.path.exists(tf):
            diffs.append((name, "-", "CAPS failed", "")); continue
        got = [t.split("_")[0] for t in open(tf).read().split()]   # CAPS's OPLS-AA 2024 names start with the number
        good = True
        for i, (x, y) in enumerate(zip(ref, got)):
            tot += 1
            if y == "?":
                untyped += 1; good = False; diffs.append((name, i + 1, x, "untyped"))
            elif x == y:
                ok += 1
            else:
                good = False; diffs.append((name, i + 1, x, y))
        full += good
    why = lambda d: EXPLAINED_MOLECULE.get(d[0]) or EXPLAINED.get((d[2], d[3]))
    unexplained = [d for d in diffs if not why(d)]
    print(f"{ff} against the GROMACS OPLS-AA topologies: {ok}/{tot} atoms ({100 * ok / max(tot, 1):.1f} %), "
          f"{untyped} untyped, {full}/{n} molecules fully the same; {len(diffs) - len(unexplained)} differing atoms with a "
          f"stated reason, {len(unexplained)} without")
    for d in unexplained[:show]:
        print(f"   {d[0]:34s} atom {str(d[1]):4s} ref {d[2]:5s} CAPS {d[3]:8s} {desc.get(d[3], '')[:60]}")
    reasons = {}
    for d in diffs:
        if why(d):
            reasons.setdefault(why(d), set()).add(d[0])
    print("explained:")
    for r, ms in reasons.items():
        print(f"   {r} [{', '.join(sorted(ms))}]")


if sys.argv[1] == "prepare":
    prepare(sys.argv[2], sys.argv[3])
else:
    compare(sys.argv[2], arg("--ff", "oplsaa2024-moltemplate"), int(arg("--show", "80")))
