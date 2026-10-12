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


# differences with a stated reason, by (reference number, CAPS number). The topologies date from 2012 and were built by
# hand; where they and the 2024 table's own descriptions disagree, CAPS follows the table.
EXPLAINED = {
    ("140", "282"): "H alpha to an ester or aryl/cyclopropyl ketone C=O: the table gives 282 ('for R on C=O use 280-282'); "
                    "the topologies use 282 in acetone, 2-hexanone and acetic anhydride but 140 here",
    ("722", "975"): "alkyl Br: the topologies use the united-atom 722; the table's all-atom alkyl bromide is 975 (JPCB 2004)",
    ("152", "976"): "C of an alkyl bromide: the table's RCH2Br 976 (the topologies borrow RCH2Cl 152)",
    ("152", "978"): "C of a secondary alkyl bromide: the table's R2CHBr 978 (the topologies borrow RCH2Cl 152)",
    ("153", "977"): "H on an alkyl bromide C: the table's 977 (the topologies borrow 153 of RCH2Cl; same charge)",
    ("152", "1012"): "C of a secondary alkyl iodide: the table's iodoalkane 1012 (the topologies borrow RCH2Cl 152)",
    ("153", "140"): "H on an alkyl iodide C: the table has no number of its own; 140",
    ("732", "1014"): "alkyl I: the table's iodoalkane 1014 (the topologies borrow iodobenzene's 732)",
    ("469", "777"): "dialkyl carbonate: the table's carbonate types 771-779/789 (the topologies use the ester types)",
    ("467", "773"): "dialkyl carbonate: the table's carbonate types (the topologies use the ester types)",
    ("465", "772"): "dialkyl carbonate: the table's carbonate types (the topologies use the ester types)",
    ("466", "771"): "dialkyl carbonate: the table's carbonate types (the topologies use the ester types)",
    ("490", "774"): "dialkyl carbonate: the table's carbonate types (the topologies use the ester types)",
    ("468", "789"): "dimethyl carbonate: the table's 'methoxy C in carbonate' 789",
    ("180", "186"): "O of an acetal (O-CH2-O): the table's OPLS/2020 acetal O 186",
    ("154", "160"): "O of a tertiary alcohol: the table's OPLS/2020 tertiary alcohol O 160",
    ("185", "1029"): "epoxide: the table's epoxide types 1025-1029 (the topologies use the ether types)",
    ("180", "1025"): "epoxide: the table's epoxide types 1025-1029 (the topologies use the ether types)",
    ("182", "1026"): "epoxide: the table's epoxide types 1025-1029 (the topologies use the ether types)",
    ("183", "1027"): "epoxide: the table's epoxide types 1025-1029 (the topologies use the ether types)",
    ("226", "398"): "Cl2C=: the table's 398/399 for two Cl on the alkene C (the topologies use ClHC='s 226/227)",
    ("227", "399"): "Cl2C=: the table's 398/399 for two Cl on the alkene C (the topologies use ClHC='s 226/227)",
    ("472", "199"): "diphenyl ether: the table's OPLS/2020 aryl ether types 177/199 (the topologies use the phenyl ester's 472/473)",
    ("473", "177"): "diphenyl ether: the table's OPLS/2020 aryl ether types 177/199 (the topologies use the phenyl ester's 472/473)",
    ("135", "273"): "CH3 of acetate: the table's 'CH3, carboxylate ion' 273",
    ("145B", "519"): "biphenyl: 145B exists only in the GROMACS file; the table's biphenyl C1 519",
    ("135", "927"): "propyne: GROMACS numbers 925-930 differently (its 930 is an H, 927 the inner sp C); in the 2024 table 927 is the CH3",
    ("930", "140"): "propyne: GROMACS's 930 is an H; in the 2024 table 930 is a carbon (R3C-C#C)",
    ("260", "142"): "acrylonitrile: the topologies give the vinyl C carrying the nitrile the benzene ring C 260; CAPS the alkene C 142",
    ("140", "144"): "acrylonitrile: the topologies give the H on that vinyl C the sp3 H 140; CAPS the alkene H 144",
    ("927", "925"): "propyne: GROMACS's 927 is the inner sp C; in the 2024 table that C is 925 and 927 the CH3",
}
EXPLAINED_MOLECULE = {"isoquinoline": "the topologies give isoquinoline quinoline's numbers 603-619 position by position, "
                                      "though its N sits elsewhere; CAPS types it as a pyridine ring fused to a benzene ring"}


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
        ref = open(os.path.join(out, name + ".ref")).read().split()
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
