#!/usr/bin/env python3
"""OPLS-AA atom types from CAPS (OPLS-AA 2024 rules) against hand-built GROMACS OPLS-AA topologies: the organic-liquid
set of Caleman et al., J. Chem. Theory Comput. 8 (2012) 61 (virtualchemistry.org; foyer ships it as opls_validation),
one .top and .gro per molecule with opls_NNN on every atom.

  compare_oplsaa_top.py prepare TOPDIR OUTDIR        (needs RDKit: bond orders from the connectivity and net charge)
  compare_oplsaa_top.py compare OUTDIR [--ff oplsaa2024] [--show 80]

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
# (topology number, CAPS type): differences where the force field's own table, not the topology, decides. The topologies
# set every atom's charge in their [atoms] section, so they could type a group with a nearby number and fix the charge;
# CAPS takes each type's own charge, so a group must get the number (or CAPS supplement) whose charges balance
CL = "chlorinated C: the topologies put RCH2Cl's 152 / 153 on every chloroalkane carbon with charges of their own; CAPS's " \
     "CHCl2 / CH2Cl2 / CCl3 / CHCl3 carbons (supplement, +0.200 per Cl as the table's 152 and 973 have it) keep each molecule " \
     "neutral, and 972 (H in RCHCl) has 153's class and charge"
EXPLAINED = {
    **{(r, c): CL for r, c in (("152", "caps_chcl2"), ("152", "caps_ch2cl2"), ("152", "caps_ccl3"), ("152", "caps_chcl3"), ("153", "972"))},
    ("180", "186"): "acetal O: the 2024 table's 186 (acetal / hemiacetal ether O, OPLS/2020); the topologies predate it",
    ("186", "caps_os_acac"): "an O between two acetal C (paraldehyde, trioxane, polyacetal): CAPS supplement at -0.26 so each "
                             "acetal C group's +0.26 balances (186's -0.33 is for one acetal C)",
    ("472", "199"): "diphenyl ether: the 2024 table's 177 / 199 (diphenyl ether, anisole C, OPLS/2020); 472 / 473 are a phenyl ester's",
    ("473", "177"): "diphenyl ether: the 2024 table's 177 / 199 (diphenyl ether, anisole C, OPLS/2020); 472 / 473 are a phenyl ester's",
    ("145", "165"): "styrene's ring C carrying the vinyl group: the 2024 table's 165 (Cipso in styrene)",
    ("465", "caps_anh_c"): "anhydride: no OPLS-AA number; CAPS supplement with the ester C / O classes and charges that sum to zero",
    ("466", "caps_anh_o"): "anhydride: no OPLS-AA number; CAPS supplement with the ester C / O classes and charges that sum to zero",
    ("467", "caps_anh_os"): "anhydride: no OPLS-AA number; CAPS supplement with the ester C / O classes and charges that sum to zero",
    ("260", "caps_cm_cn"): "acrylonitrile's alkene C: the topology uses benzonitrile's ring C 260 (class CA); CAPS keeps the alkene "
                           "class CM (supplement, the files' CM-CZ terms)",
    ("140", "caps_hc_cn"): "the H on acrylonitrile's alkene C: alkene H 144's class and charge (the topology: alkane H 140)",
    ("142", "caps_cm_ves"): "vinyl ester C(O)=: CAPS supplement, alkene class CM charged as the table's phenyl ester ipso C",
    ("484", "136"): "sulfolane: the topology's sulfonamide alpha-C types (484 +0.39, H 485 -0.06) leave it +0.54 with the types' "
                    "charges; the table has no sulfone alpha-C, so CAPS takes the alkane numbers beside the neutral SO2",
    ("485", "140"): "sulfolane: the topology's sulfonamide alpha-C types (484 +0.39, H 485 -0.06) leave it +0.54 with the types' "
                    "charges; the table has no sulfone alpha-C, so CAPS takes the alkane numbers beside the neutral SO2",
    ("135", "927"): "propyne's CH3 (2024 table): H3C-C#C 927 (+0.02), whose group balances the terminal alkyne's 925 / 925 / 926 "
                    "(-0.20); the topology's alkane 135 is GROMACS's own alkyne set, where the inner sp C carries the balance",
    ("154", "160"): "tertiary alcohol O (2024 table): 160 (-0.65, OPLS/2020) beside its tertiary C 159 (0.232); with 154 the "
                    "alcohol is -0.033 e (GROMACS's 159 is 0.265, which balances 154)",
    ("135", "273"): "acetate's CH3: GROMACS's own carboxylate CH3 273 (-0.28), whose charges sum to -1 with 271 / 272; the topology "
                    "uses 135 with charges of its own",
}
BR = "alkyl bromide / iodide (2024 table): its all-atom alkyl bromide 975-979 and iodoalkane 1011-1015, which sum to zero " \
     "for RCH2X, R2CHX, R3CX; the topologies' 722 (a united-atom Br) / 732 (iodobenzene's I) with 152 / 153 need charges of their own"
CO3 = "open-chain carbonate (2024 table): its carbonate numbers 771-773 with the methoxy C 789 / 774 and H 777 sum to zero; " \
      "the topologies' ester types leave each carbonate at -0.08 e with the table's charges"
EXPLAINED.update({**{p: BR for p in (("152", "976"), ("152", "978"), ("153", "977"), ("722", "975"), ("152", "1012"),
                                      ("153", "1015"), ("732", "1014"))},
                  **{p: CO3 for p in (("469", "777"), ("467", "773"), ("490", "774"), ("468", "789"), ("466", "771"), ("465", "772"))},
                  ("227", "399"): "Cl2C= (2024 table): its own pair 399 (+0.12) / 398 (-0.06), neutral; 227 + 2 x 226 is -0.235 e",
                  ("226", "398"): "Cl2C= (2024 table): its own pair 399 (+0.12) / 398 (-0.06), neutral; 227 + 2 x 226 is -0.235 e",
                  ("186", "180"): "vinyl ether O (2024 table): 186 is there the OPLS/2020 acetal O (-0.33, +0.07 e left); the ether O 180 "
                                  "(-0.40) has the charge the topologies' 186 has"})
UNBAL = " has no numbers in the table whose charges balance; left untyped, CAPS reports it (with the topology's numbers and the " \
        "table's charges it would carry a net charge)"
# GROMACS's oplsaa.ff has no numbers whose own charges balance these groups (the topologies set charges of their own in
# [atoms]); the GROMACS-named library leaves them untyped and reports them, as the 2024 library does (check_type_charges.py)
GMX_UNBAL = {m: why for why, ms in {
    "alkyl bromide / iodide: Br 722 is a united-atom type (-0.22) and I 732 iodobenzene's (-0.10); with 152 / 153 the "
    "molecule is charged and the file has no all-atom alkyl bromide or iodoalkane": ("1-bromobutane", "1-bromopropane",
        "12-dibromoethane", "12-dibromopropane", "2-iodopropane", "bromoethane", "bromomethane", "dibromomethane", "iodomethane"),
    "polychlorinated C: 152 / 153 are RCH2Cl's; on CHCl2, CH2Cl2, CHCl3, CCl3 the molecule is charged (chloroform -0.50 e)":
        ("11-dichloroethane", "112-trichloroethane", "1122-tetrachloroethane", "dichloromethane", "pentachloroethane",
         "trichloromethane", "dichlorofluoromethane"),
    "Cl2C=: the file has no 398 / 399; 227 + 2 x 226 is -0.235 e": ("11-dichloroethene",),
    "acrylonitrile: no alkene C carrying a nitrile (260 is benzonitrile's ring C, +0.06 e on acrylonitrile)": ("2-propenenitrile",),
    "anhydride: the ester numbers leave C(=O)-O-C(=O) at -0.17 e": ("acetic-anhydride",),
    "biphenyl: 145B carries -0.115 with no H to balance it (-0.23 e)": ("biphenyl",),
    "dimethyl carbonate: the file has no methoxy C of a carbonate (789); the ester numbers leave -0.08 e": ("dimethylcarbonate",),
    "diaryl ether: no 177; the phenyl ester pair 472 / 473 leaves +0.055 e a link": ("diphenyl-ether",),
    "styrene: its 165 is not the 2024 table's styrene ipso C; benzene's 145 there leaves -0.115 e": ("styrene",),
    "vinyl ester: the phenyl ester O 473 with the alkene numbers leaves -0.135 e": ("vinyl-acetate",),
}.items() for m in ms}
EXPLAINED_MOLECULE = {"methyloxirane": "epoxide: the table's epoxide types 1025-1029 with three-ring angle terms; the topologies' "
                                       "ether types put 109.5 deg angles on the ring (about 119 kcal/mol strain, 44 with these)",
                      "triethyl-phosphate": "a phosphate triester: the table's 440-444 are the dimethyl phosphate anion's (-0.38 e on "
                                            "triethyl phosphate with their charges); left untyped in both libraries",
                      "thiophene": "thiophene (2024 table)" + UNBAL + " (furan's C / H with thiazole's S: +0.16 e)",
                      "bromomethane": "CH3Br (2024 table)" + UNBAL + " (976 with three H: +0.10 e)",
                      "dibromomethane": "CH2Br2 (2024 table)" + UNBAL,
                      "dichlorofluoromethane": "a C with F and Cl (2024 table)" + UNBAL,
                      "3-methylpyridine": "a 3- / 4-alkylpyridine (2024 table)" + UNBAL + " (672 / 674 are a pyridazine's and a pyrimidine's CH3)",
                      "4-methylpyridine": "a 3- / 4-alkylpyridine (2024 table)" + UNBAL + " (672 / 674 are a pyridazine's and a pyrimidine's CH3)",
                      "246-trimethylpyridine": "a 3- / 4-alkylpyridine (2024 table)" + UNBAL + " (672 / 674 are a pyridazine's and a pyrimidine's CH3)"}


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
        if ff != "oplsaa-gromacs":   # GROMACS's own numbers stay for its own library
            ref = [GROMACS_TO_2024.get(t, t) for t in ref]
        tf = os.path.join(out, name + ".caps")
        subprocess.run([caps, "ff", "type", os.path.join(out, f), "--ff", ff, "-o", tf], capture_output=True, text=True)
        n += 1
        if not os.path.exists(tf):
            diffs.append((name, "-", "CAPS failed", "")); continue
        got = [t[5:].lstrip("0") or "0" if re.fullmatch(r"opls_\d+[A-Z]?", t) else t for t in open(tf).read().split()]   # opls_NNN: the number
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
    why = lambda d: EXPLAINED_MOLECULE.get(d[0]) or EXPLAINED.get((d[2], d[3])) or \
        (GMX_UNBAL.get(d[0]) and ff == "oplsaa-gromacs" and GMX_UNBAL[d[0]] + (
            "; untyped, CAPS refuses the group" if d[3] == "untyped" else
            "; CAPS types the group with the OPLS-AA 2024 table's numbers (fallback types opls2024_NNN / caps_*)"))
    unexplained = [d for d in diffs if not why(d)]
    print(f"{ff} against the GROMACS OPLS-AA topologies: {ok}/{tot} atoms ({100 * ok / max(tot, 1):.1f} %), "
          f"{untyped} untyped, {full}/{n} molecules fully the same; {len(diffs) - len(unexplained)} differing atoms with a "
          f"stated reason, {len(unexplained)} without")
    for d in unexplained[:show]:
        print(f"   {d[0]:34s} atom {str(d[1]):4s} ref {d[2]:5s} CAPS {d[3]:8s} {desc.get('opls_' + d[3], desc.get(d[3], ''))[:60]}")
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
    compare(sys.argv[2], arg("--ff", "oplsaa2024"), int(arg("--show", "80")))
