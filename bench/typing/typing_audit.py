#!/usr/bin/env python3
"""Atom types of CAPS's organic force fields against independent typers, on one molecule set.

The set: bench/typing/gaff_reference_set.smi, bench/typing/rubber_set.smi (rubbers, curatives, accelerators,
antioxidants, coupling agents, fibres, solvents; repeat units with [*] are built as DP-3 oligomers), the solvents of
data/reference/solvents.json and every single-unit polymer of data/polymers/library.json as a DP-3 oligomer.

  typing_audit.py prepare OUT [--caps build/cli/caps]
      OUT/mol/NAME.mol2 (3D, explicit H, CAPS's bond orders, element symbols in the type column) and OUT/set.json
  typing_audit.py caps OUT [--ff ID,ID]            CAPS's types: OUT/caps/FF/NAME.types
  typing_audit.py antechamber OUT --amber BIN      AmberTools: OUT/ref/gaff/NAME.types, OUT/ref/gaff2/NAME.types
  typing_audit.py foyer OUT --python PY            foyer's OPLS-AA (MoSDeF): OUT/ref/foyer-oplsaa/NAME.types
  typing_audit.py pysimm OUT --python PY           pysimm's PCFF: OUT/ref/pysimm-pcff/NAME.types
  typing_audit.py compare OUT [--show 60] [--ff ID]

Each comparison counts the atoms the reference types with a type that has parameters in the CAPS library; a reference
that cannot type a molecule leaves that molecule out of its count (listed). Differences that are not CAPS errors are
listed apart with the reason (EXPLAINED), each checked against the force field's own definition.
"""
import glob, json, os, re, subprocess, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
CAPS = arg("--caps", os.path.join(ROOT, "build", "cli", "caps"))
TWO = {"Cl", "Br", "Si", "Na", "Mg", "Li", "Al", "Zn", "Ca", "Se", "Sn", "Ge", "As", "Ti", "Fe", "Cu", "Ni", "Co", "Mn",
       "Cr", "Zr", "Ag", "Au", "Pt", "Pd", "Ba", "Sr", "Cs", "Rb", "Xe", "Kr", "Ar", "Ne", "He", "Be", "Ga", "In", "Sb", "Te"}

# the libraries compared and their reference: (CAPS id, reference directory, how names are compared)
PAIRS = [
    ("gaff-1.81-ambertools", "gaff", "name"),
    ("gaff2-2.2.20-ambertools", "gaff2", "name"),
    ("gaff2-2.11-ambertools", "gaff2", "gaff2.11"),
    ("oplsaa2024", "foyer-oplsaa", "opls"),
    ("oplsaa-gromacs", "foyer-oplsaa", "opls"),
    ("pcff-frc", "pysimm-pcff", "pcff-equiv"),
    ("dreiding1990", "pysimm-dreiding", "name"),
]


def element(name):
    """the element of an atom name such as C12, Cl3, CL3, Si1 (two letters only when they make an element)"""
    s = re.match(r"([A-Za-z]+)", name).group(1)
    two = s[:2].capitalize()
    return two if len(s) >= 2 and two in TWO else s[0].upper()


def read_mol2(path):
    atoms, bonds, sec = [], [], None
    for line in open(path):
        if line.startswith("@<TRIPOS>"):
            sec = line.strip(); continue
        w = line.split()
        if sec == "@<TRIPOS>ATOM" and len(w) >= 6:
            atoms.append(w)
        elif sec == "@<TRIPOS>BOND" and len(w) >= 4:
            bonds.append((int(w[1]) - 1, int(w[2]) - 1, w[3]))
    return atoms, bonds


def write_element_mol2(src, dst, title):
    """the same molecule with element symbols in the atom-type column (what antechamber, parmed and pysimm read)"""
    atoms, bonds = read_mol2(src)
    els = [w[5] if re.fullmatch(r"[A-Z][a-z]?", w[5]) else element(w[1]) for w in atoms]
    with open(dst, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{title}\n{len(atoms)} {len(bonds)} 1 0 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
        cnt = Counter()
        for i, (w, el) in enumerate(zip(atoms, els)):
            cnt[el] += 1
            f.write("%7d %-6s %10.4f %10.4f %10.4f %-4s %3d %-4s %10.6f\n" % (i + 1, f"{el}{cnt[el]}", float(w[2]), float(w[3]), float(w[4]), el, 1, "MOL", 0.0))
        f.write("@<TRIPOS>BOND\n")
        for k, (i, j, o) in enumerate(bonds):
            f.write("%6d %5d %5d %s\n" % (k + 1, i + 1, j + 1, o))
    return els


def molecule_set():
    mols = {}

    def add(name, smi, src):
        name = re.sub(r"[^A-Za-z0-9_+-]", "_", name)
        if name not in mols:
            mols[name] = {"smiles": smi, "set": src, "polymer": "*" in smi}
    for f in ("gaff_reference_set.smi", "rubber_set.smi"):
        for line in open(os.path.join(HERE, f)):
            if line.strip() and not line.startswith("#"):
                w = line.split()
                add(w[1], w[0], f)
    for s in json.load(open(os.path.join(ROOT, "data", "reference", "solvents.json")))["solvents"]:
        add(s["name"].lower().replace(" ", "_"), s["smiles"], "solvents.json")
    for p in json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))["polymers"]:
        if isinstance(p.get("smiles"), str):
            add(p["id"], p["smiles"], "polymers/library.json")
    return mols


def prepare(out):
    os.makedirs(os.path.join(out, "mol"), exist_ok=True)
    os.makedirs(os.path.join(out, "raw"), exist_ok=True)
    mols = molecule_set()
    sys.path.insert(0, os.path.join(ROOT, "data", "python"))
    os.environ.setdefault("CAPS_LIB", os.path.join(ROOT, "build", "capi", "libcaps.dylib" if sys.platform == "darwin" else "libcaps.so"))
    import caps
    ok = {}
    for name, m in mols.items():
        raw = os.path.join(out, "raw", name + ".mol2")
        try:
            if m["polymer"]:
                caps.polymer(m["smiles"], dp=3).save(raw)
            else:
                r = subprocess.run([CAPS, "build", m["smiles"], "-o", raw, "--seed", "1"], capture_output=True, text=True)
                if r.returncode:
                    raise RuntimeError(r.stderr.strip()[-200:] or r.stdout.strip()[-200:])
            m["elements"] = write_element_mol2(raw, os.path.join(out, "mol", name + ".mol2"), name)
            ok[name] = m
        except Exception as e:
            print("build failed:", name, m["smiles"], str(e)[:200])
    json.dump(ok, open(os.path.join(out, "set.json"), "w"), indent=1)
    print(len(ok), "of", len(mols), "molecules built in", out)


def names(out):
    return list(json.load(open(os.path.join(out, "set.json"))))


def run_caps(out):
    ffs = arg("--ff", ",".join(dict.fromkeys(p[0] for p in PAIRS))).split(",")
    for ff in ffs:
        d = os.path.join(out, "caps", ff)
        os.makedirs(d, exist_ok=True)
        for n in names(out):
            subprocess.run([CAPS, "ff", "type", os.path.join(out, "mol", n + ".mol2"), "--ff", ff, "-o", os.path.join(d, n + ".types")],
                           capture_output=True, text=True)
        print(ff, len(glob.glob(os.path.join(d, "*.types"))), "molecules typed")


def run_antechamber(out):
    amber = arg("--amber")
    exe = os.path.join(amber, "antechamber") if amber else "antechamber"
    work = os.path.join(out, "ac")
    os.makedirs(work, exist_ok=True)
    for at in ("gaff", "gaff2"):
        d = os.path.join(out, "ref", at)
        os.makedirs(d, exist_ok=True)
        for n in names(out):
            o = os.path.join(work, f"{n}.{at}.mol2")
            if not os.path.exists(o):
                subprocess.run([exe, "-i", os.path.join(out, "mol", n + ".mol2"), "-fi", "mol2", "-o", o, "-fo", "mol2", "-at", at,
                                "-dr", "no", "-pf", "y"], capture_output=True, text=True, cwd=work)
            if os.path.exists(o):
                a, _ = read_mol2(o)
                open(os.path.join(d, n + ".types"), "w").write("\n".join(w[5] for w in a) + "\n")
        print(at, len(glob.glob(os.path.join(d, "*.types"))), "molecules typed by antechamber")


FOYER = r'''
import json, os, sys, warnings
warnings.filterwarnings("ignore")
import parmed, foyer
from foyer.atomtyper import find_atomtypes
ff = foyer.Forcefield(name=sys.argv[3])
out = {}
for n in json.load(open(sys.argv[1])):
    try:
        s = parmed.load_file(os.path.join(sys.argv[2], n + ".mol2"), structure=True)
        tm = find_atomtypes(s, forcefield=ff)
        out[n] = [tm[i]["atomtype"] for i in range(len(s.atoms))]
    except Exception as e:
        out[n] = "foyer: " + str(e).splitlines()[0][:300]
json.dump(out, open(sys.argv[4], "w"))
'''


def run_foyer(out):
    py = arg("--python", "python3")
    script = os.path.join(out, "foyer_types.py")
    open(script, "w").write(FOYER)
    res = os.path.join(out, "foyer-oplsaa.json")
    subprocess.run([py, script, os.path.join(out, "set.json"), os.path.join(out, "mol"), "oplsaa", res], check=True)
    d = os.path.join(out, "ref", "foyer-oplsaa")
    os.makedirs(d, exist_ok=True)
    bad = {}
    for n, t in json.load(open(res)).items():
        if isinstance(t, list):
            open(os.path.join(d, n + ".types"), "w").write("\n".join(t) + "\n")
        else:
            bad[n] = t
    json.dump(bad, open(os.path.join(out, "ref", "foyer-oplsaa.failed.json"), "w"), indent=1)
    print("foyer typed", len(glob.glob(os.path.join(d, "*.types"))), "molecules;", len(bad), "it cannot type")


PYSIMM = r'''
import contextlib, io, json, os, sys, warnings
warnings.filterwarnings("ignore")
from pysimm import system, forcefield
ff = forcefield.Pcff()
out = {}
ORD = {"1": 1, "2": 2, "3": 3, "ar": "A", "am": 1}
for n in json.load(open(sys.argv[1])):
    atoms, bonds, sec = [], [], None
    for line in open(os.path.join(sys.argv[2], n + ".mol2")):
        if line.startswith("@<TRIPOS>"):
            sec = line.strip(); continue
        w = line.split()
        if sec == "@<TRIPOS>ATOM" and len(w) >= 6: atoms.append(w)
        elif sec == "@<TRIPOS>BOND" and len(w) >= 4: bonds.append(w)
    # pysimm types aromatic C / N only on aromatic bond orders: antechamber's (same atoms, "ar" bonds) where it ran
    ac = os.path.join(sys.argv[4], n + ".gaff.mol2")
    if os.path.exists(ac):
        bonds, sec = [], None
        for line in open(ac):
            if line.startswith("@<TRIPOS>"):
                sec = line.strip(); continue
            w = line.split()
            if sec == "@<TRIPOS>BOND" and len(w) >= 4: bonds.append(w)
    s = system.System()
    for k, w in enumerate(atoms):   # the element from the type column (pysimm's read_mol2 takes the name's first letter)
        s.particles.add(system.Particle(tag=k + 1, elem=w[5], charge=0.0, molecule=1, x=float(w[2]), y=float(w[3]), z=float(w[4])))
    for k, w in enumerate(bonds):
        s.bonds.add(system.Bond(tag=k + 1, a=s.particles[int(w[1])], b=s.particles[int(w[2])], order=ORD.get(w[3], 1)))
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            r = ff.assign_ptypes(s)
        out[n] = [p.type.name if p.type else "?" for p in s.particles]
    except Exception as e:
        out[n] = "pysimm: " + type(e).__name__ + " " + str(e)[:200]
json.dump(out, open(sys.argv[3], "w"))
'''


def run_pysimm(out):
    """after antechamber (its mol2 files give pysimm the aromatic bond orders it needs)"""
    py = arg("--python", "python3")
    for cls, tag in (("Pcff", "pysimm-pcff"), ("Dreiding", "pysimm-dreiding")):
        script = os.path.join(out, "pysimm_types.py")
        open(script, "w").write(PYSIMM.replace("forcefield.Pcff()", f"forcefield.{cls}()"))
        res = os.path.join(out, tag + ".json")
        subprocess.run([py, script, os.path.join(out, "set.json"), os.path.join(out, "mol"), res, os.path.join(out, "ac")], check=True)
        d = os.path.join(out, "ref", tag)
        os.makedirs(d, exist_ok=True)
        bad = {}
        for n, t in json.load(open(res)).items():
            if isinstance(t, list) and "?" not in t:
                open(os.path.join(d, n + ".types"), "w").write("\n".join(t) + "\n")
            else:
                bad[n] = t if isinstance(t, str) else "pysimm leaves atoms untyped"
        json.dump(bad, open(os.path.join(out, "ref", tag + ".failed.json"), "w"), indent=1)
        print(f"pysimm ({cls}) typed", len(glob.glob(os.path.join(d, "*.types"))), "molecules;", len(bad), "it cannot type")


def types_file(path):
    return [w.split()[-1] if w.split() else "" for w in open(path) if not w.startswith("#")]


def library(ff):
    return json.load(open(os.path.join(ROOT, "data", "forcefields", ff + ".json")))


def comparer(ff, mode):
    """(same(ref, caps) -> bool, has_parameters(ref) -> bool)"""
    lib = library(ff)
    tnames = {t["name"] for t in lib["atom_types"]}
    if mode == "name":
        return (lambda r, c: r.lower() == c.lower()), (lambda r: r in tnames)
    if mode == "gaff2.11":   # GAFF2 2.11 has no c5 / c6 (sp3 C in five- and six-membered rings, added in 2.2): c3 there
        m = lambda r: "c3" if r in ("c5", "c6") else r
        return (lambda r, c: m(r) == c), (lambda r: m(r) in tnames)
    if mode == "opls":
        num = lambda t: (re.match(r"^opls_(\d+)", t or "") or [None, None])[1]
        return (lambda r, c: num(r) is not None and num(r) == num(c)), (lambda r: num(r) is not None and "opls_%03d" % int(num(r)) in tnames)
    if mode == "pcff-equiv":   # the same parameters: both names in one nonbond and one bond class through pcff.frc's equivalences
        eq = {t["name"]: (t.get("equivalence", {}).get("vdw", t["name"]), t.get("equivalence", {}).get("bond", t["name"])) for t in lib["atom_types"]}
        return (lambda r, c: r == c or (r in eq and c in eq and eq[r] == eq[c])), (lambda r: r in tnames)
    raise ValueError(mode)


GAFFS = ("gaff-1.81-ambertools", "gaff2-2.2.20-ambertools", "gaff2-2.11-ambertools")
OPLS = ("oplsaa2024", "oplsaa-gromacs")
# foyer's OPLS-AA SMARTS against OPLS-AA's own type table (both libraries' descriptions): where the table has a number for
# the environment, foyer often gives the generic one. (foyer number, CAPS number) -> reason
SPEC = "the table's own number for this environment, foyer the generic one"
FOYER_OPLS = {
    ("282", "140"): "H on the alkyl C of an ester's acyl R: 282 (H alpha to C=O) and 140 have one class, charge and LJ",
    ("465", "471"): "aryl ester C=O: 471 (Co in methyl benzoate, aryl ester): " + SPEC,
    ("472", "199"): "aryl ether C: 199 (anisole C) / 177 (diphenyl ether); foyer's 472 / 473 are a phenyl ester's",
    ("473", "177"): "aryl ether O: 177 (diphenyl ether, OPLS/2020); foyer's 473 is a phenyl ester's -O-",
    ("145", "488"): "ring C on a sulfone S: 488 (CA, 0.0) keeps the SO2 group and the ring neutral; 145 is a C-H carbon (-0.115)",
    ("142", "150"): "conjugated diene inner CH: 150 (diene =CH-CH=): " + SPEC,
    ("141", "178"): "conjugated diene inner CR: 178 (diene =CR-RC=): " + SPEC,
    ("137", "515"): "benzylic CH: 515 (isopropyl benzene CH): " + SPEC,
    ("139", "516"): "benzylic C: 516 (t-butyl benzene C): " + SPEC,
    ("203", "caps_s_poly"): "S between two S (trisulfide and longer): no OPLS number; CAPS supplement (203's class and LJ, 0 charge)",
    ("145", "519"): "biphenyl junction C: 519 (biphenyl C1, class C!): " + SPEC,
    ("145", "145B"): "biphenyl junction C: GROMACS's 145B ('Use #145B for biphenyl'), class C!",
    ("145", "1004"): "ring C on a diaryl-substituted sp3 C: 1004 (diphenylmethane Cipso): " + SPEC,
    ("238", "265"): "anilide N: 265 (N-phenylacetamide N): " + SPEC,
    ("145", "266"): "anilide ring C: 266 (ipso C in N-phenylacetamide): " + SPEC,
    ("183", "491"): "ester alkoxy CH: 491 (C(HOS) i-Pr ester); foyer's 183 is an ether's",
    ("185", "469"): "H on an ester alkoxy C: 469 (alkoxy H in esters); foyer's 185 is an ether's",
    ("145", "165"): "styrene ring C carrying the vinyl group: 165 (Cipso in styrene): " + SPEC,
    ("235", "234"): "benzamide C=O: 234 (C=O in benzamide): " + SPEC,
    ("136", "714"): "cyclobutane CH2: 714 (OPLS-2020 cyclobutane): " + SPEC,
    ("267", "470"): "benzoic acid C=O: 470 (Co in benzoic acid): " + SPEC,
    ("240", "250"): "urea: 247-250 (O=C(NH2)2): " + SPEC, ("237", "249"): "urea: 247-250 (O=C(NH2)2): " + SPEC,
    ("235", "247"): "urea: 247-250 (O=C(NH2)2): " + SPEC, ("236", "248"): "urea: 247-250 (O=C(NH2)2): " + SPEC,
    ("282", "144"): "H on the sp2 C of an enone / enal: alkene H 144; foyer's 282 is the H on an sp3 C alpha to C=O",
    ("233", "231"): "diaryl ketone C=O: 231 (benzophenone): " + SPEC,
    ("145", "221"): "ring C on a benzylic C-O: 221 (benzyl alcohol / ether ring C): " + SPEC,
    ("927", "939"): "2-butyne's sp C: 939 (RC#CR, MeCCMe); foyer gives the 2024 table's 927, an sp3 H3C-C#C carbon",
    ("927", "931"): "2-butyne's sp C: GROMACS's 931 (alkyne RC#CR); its 927 is the inner C of a terminal alkyne",
    ("145", "916"): "aniline ring C(N): 916 (C(NH2) aniline): " + SPEC,
    ("965", "726"): "CF2H end group: 966 / 967 / 726 (CF2H): " + SPEC, ("962", "966"): "CF2H end group: 966 / 967 / 726 (CF2H): " + SPEC,
    ("140", "967"): "CF2H end group: 966 / 967 / 726 (CF2H): " + SPEC,
    ("137", "caps_pyr2ch"): "2-alkylpyridine CH: CAPS supplement charged as 671 / 670 (2-ethyl / 2-methyl pyridine)",
    ("136", "671"): "2-alkylpyridine CH2: 671 (CH2, 2-ethyl pyridine): " + SPEC,
    ("185", "192"): "polyoxymethylene oligomer: the table's acetal / hemiacetal set 186-192: " + SPEC,
    ("180", "186"): "acetal O: 186 (acetal ether O): " + SPEC, ("157", "191"): "hemiacetal C: 191: " + SPEC,
    ("154", "187"): "hemiacetal OH: 187: " + SPEC, ("155", "188"): "hemiacetal H(O): 188: " + SPEC,
    ("180", "caps_os_acac"): "O between two acetal C: CAPS supplement at -0.26 so each acetal C group's +0.26 balances",
    ("183", "caps_bzo_ch"): "benzyl ether CH: CAPS supplement (183 + 0.055 for the ring, as 219 is to 158)",
    ("158", "219"): "benzylic alcohol CH: 219 (CHROH benzyl alcohols): " + SPEC,
    ("465", "caps_anh_c"): "anhydride: no OPLS number; CAPS supplement with the ester classes, charges summing to zero",
    ("466", "caps_anh_o"): "anhydride: no OPLS number; CAPS supplement with the ester classes, charges summing to zero",
    ("467", "caps_anh_os"): "anhydride: no OPLS number; CAPS supplement with the ester classes, charges summing to zero",
    ("135", "273"): "acetate CH3: GROMACS's carboxylate CH3 273, whose charges sum to -1 with 271 / 272",
    ("143", "1154"): "allene: the 2024 table's allene types 1153-1157", ("144", "1153"): "allene: the 2024 table's allene types 1153-1157",
    ("927", "1157"): "allene: the 2024 table's allene types 1153-1157 (foyer: an sp3 alkyne-substituent C for the central C)",
    **{p: "open-chain carbonate: the 2024 table's carbonate numbers 771-773 with the methoxy C 789 and H 777, which sum to "
          "zero (foyer's ester types leave -0.08 e)" for p in (("465", "772"), ("466", "771"), ("467", "773"), ("468", "789"), ("469", "777"))},
}
# molecules a library leaves untyped because no numbers of its table balance their charges (check_type_charges.py)
UNTYPED_BALANCE = {
    "trimethyl_phosphate": (OPLS, "a phosphate triester: the table's 440-444 are the dimethyl phosphate anion's (-0.47 e here)"),
    "thiophene": (("oplsaa2024",), "thiophene: no numbers in the 2024 table (furan's C / H with thiazole's S give +0.16 e)"),
    "lutidine": (("oplsaa2024",), "a 3-alkylpyridine: 672 is a pyridazine's CH3 and does not carry the ring H's charge (-0.16 e)"),
}
num = lambda t: re.sub(r"^opls_0*", "", t or "") if (t or "").startswith("opls_") else (t or "")
# pysimm's PCFF typer (pysimm/forcefield/pcff.py) is a short element-and-bond-order list: every non-aromatic sp2 C c=2,
# every 3-coordinate N n, every O with a non-sp3 neighbour o_2. pcff.frc's own descriptions and terms decide:
PYSIMM_PCFF = {
    ("c=2", "c="): "isolated C=C: c= is cff91's doubly bonded carbon; c=1 / c=2 were added for conjugated systems (pcff.frc "
                   "ref 3) and share c='s nonbond, angle and bond-increment classes",
    ("c=2", "c=1"): "conjugated diene carbon: c=1 (next to end doubly bonded carbon), c=-c=1 1.34 A in pcff.frc",
    ("c=2", "cp"): "aromatic five-ring C (furan, pyrrole, thiophene): cp (sp2 aromatic carbon); pysimm sees no aromatic bond",
    ("c_0", "cp"): "aromatic five-ring C next to O: cp (sp2 aromatic carbon); pysimm's c_0 is an aldehyde / ketone C",
    ("c=2", "cs"): "thiophene C next to S: cs (pcff.frc: aromatic C in a 5-membered ring next to S)",
    ("o_2", "oc"): "ether O on an aryl or vinyl C: oc (sp3 oxygen in ether or acetals); o_2 is the ester O",
    ("o_2", "o"): "peroxide O: o (generic sp3 oxygen); not an ester",
    ("o_2", "op"): "furan-type O: op (sp2 aromatic in 5 membered ring)",
    ("o_2", "oz"): "carbonate -O-: oz (ester oxygen in carbonate)", ("o_1", "oo"): "carbonate C=O: oo (carbonate only)",
    ("o", "o_2"): "carboxylic acid -O-: o_2, with ho2 (pcff.frc's o_2-ho2 bond, ver 2.1 ref 6)",
    ("h*", "ho2"): "carboxylic acid H: ho2 (hydroxyl hydrogen; the o_2-ho2 bond); h* is water's",
    ("h*", "hn2"): "urea / carbamate N-H: hn2 (amino hydrogen); h* is water's",
    ("n", "nb"): "N on an aromatic ring: nb (sp2 nitrogen in aromatic amines)",
    ("n", "na"): "amine N: na (sp3 nitrogen in amines); n is the amide N",
    ("nn", "np"): "two-coordinate ring N (pyridine, azoles): np (sp2 nitrogen in 5- or 6-membered ring, 2 connections)",
    ("nn", "nh"): "N-H in an aromatic five-ring (imidazole, pyrazole): nh (sp2 nitrogen in 5 or 6 membered ring, 3 connections)",
    ("n", "nh"): "N-H in an aromatic five-ring (pyrrole, indole, carbazole): nh",
    ("n", "n_2"): "urea / carbamate N: n_2 (nitrogen of urethane), beside c_2 (carbonyl carbon of carbamate, urea)",
    ("n", "n2"): "guanidinium NH2: n2", ("cr", "c+"): "guanidinium C: c+ (C in guanidinium group); cr is neutral arginine's",
    ("o_1", "o-"): "carboxylate O: o- (partial double oxygen)", ("c_1", "c-"): "carboxylate C: c- (C in charged carboxylate)",
    ("c=2", "c_1"): "formamide C=O (DMF): c_1 (carbonyl carbon of acid, ester, amide)",
    ("s'", "s"): "sulfate S: s with o, the pair pcff.frc's bond increments give (o-s, ver 2 ref 2); s' is the thioketone / sulfone S",
    ("o=", "o"): "sulfate / nitro O: o with s (o-s increments); pcff.frc has no nitro types (see the report)",
    ("o=", "o-"): "sulfonate O: o- (partial double oxygen) beside sf (S in sulfonate group)",
    ("s", "sp"): "thiophene-type S: sp (sulfur in an aromatic ring, e.g. thiophene)",
    ("n", "npc"): "nitro N / N-substituted carbazole N: npc (pcff.frc has no nitro N; for carbazole, sp2 N in a ring with a heavy atom)",
}
# pysimm's DREIDING typer (pysimm/forcefield/dreiding.py) never gives H_HB and reads aromaticity from bond orders only;
# the paper (Mayo, Olafson & Goddard, J. Phys. Chem. 94, 8897 (1990), section II.A) decides:
PYSIMM_DREIDING = {
    ("H_", "H__HB"): "H on N or O: H_HB, 'a hydrogen atom capable of forming hydrogen bonds' (section II.A, II.I)",
    ("N_2", "N_R"): "aromatic N, amide / amidine / guanidine N: R, 'an sp2 atom involved in a resonance situation'",
    ("N_3", "N_R"): "amide / anilide / guanidine N: R, a planar sp2 atom in resonance (not tetrahedral)",
    ("C_2", "C_R"): "aromatic five-ring C (furan, pyrrole, thiophene): R (resonant, aromatic)",
    ("O_3", "O_R"): "aromatic O (furan, oxazole): R (resonant, aromatic)",
    ("O_3", "O_2"): "one-bonded O of a carboxylate, nitro or sulfonate group: trigonal _2, not tetrahedral",
}
# (ff, test(mol, ref, caps, elements, index) -> bool, reason): differences that are not CAPS errors
EXPLAINED = [(ff, t, w) for ff in GAFFS for t, w in [
    (lambda m, r, c, e, i: m == "cyclooctatetraene" and {r, c} == {"cc", "cd"},
     "antechamber re-assigns the Kekule structure; both alternate cc / cd by their own double bonds (same parameters)"),
    (lambda m, r, c, e, i: m == "water" and (r, c) in (("oh", "ow"), ("ho", "hw")),
     "water: CAPS types it as GAFF's own water (ow / hw of gaff.dat), antechamber as an alcohol O-H (make_gaff_rules.py)"),
]] + [(ff, (lambda m, r, c, e, i: (num(r), num(c)) in FOYER_OPLS), "foyer") for ff in OPLS] + [
    ("oplsaa-gromacs", lambda m, r, c, e, i: m == "allene" and c in ("", "?"),
     "allene: GROMACS's oplsaa.ff has no allene types (the 2024 table's 1153-1157); CAPS leaves it untyped, foyer gives the "
     "central C 927, an alkyne C"),
    ("oplsaa-gromacs", lambda m, r, c, e, i: m == "P007" and num(r) == "962" and c in ("", "?"),
     "the CF2H end of the hydrogen-capped PVDF oligomer: GROMACS's file has no CF2H type (the 2024 table's 966 / 967); "
     "962 with its two F 965 and an H 140 would leave the group at +0.06"),
]


def compare(out):
    show = int(arg("--show", "60"))
    only = arg("--ff")
    mset = json.load(open(os.path.join(out, "set.json")))
    rows = []
    for ff, ref, mode in PAIRS:
        if only and ff != only:
            continue
        same, haspar = comparer(ff, mode)
        n = ok = nexp = nopar = mol_ok = mol_n = 0
        diffs, expl, missing = Counter(), defaultdict(set), []
        examples = defaultdict(list)
        for name in mset:
            rp = os.path.join(out, "ref", ref, name + ".types")
            cp = os.path.join(out, "caps", ff, name + ".types")
            if not os.path.exists(rp):
                missing.append(name); continue
            if not os.path.exists(cp):
                diffs[("", "(CAPS failed)", "")] += 1; examples[("", "(CAPS failed)", "")].append(name); continue
            r, c = types_file(rp), types_file(cp)
            if len(r) != len(c):
                diffs[("", "(atom count)", "")] += 1; continue
            els = mset[name].get("elements", [])
            mol_n += 1
            good = True
            for i, (x, y) in enumerate(zip(r, c)):
                if not haspar(x):
                    nopar += 1
                    continue
                n += 1
                if same(x, y):
                    ok += 1
                    continue
                why = next((w for f_, t_, w in EXPLAINED if f_ == ff and t_(name, x, y, els, i)), None)
                if why == "foyer":
                    why = FOYER_OPLS[(num(x), num(y))]
                if why is None and name in UNTYPED_BALANCE and ff in UNTYPED_BALANCE[name][0]:
                    why = UNTYPED_BALANCE[name][1] + "; left untyped (CAPS reports it), the rest of the molecule as the reference"
                if why is None and ff == "pcff-frc" and (x, y) in PYSIMM_PCFF:
                    why = PYSIMM_PCFF[(x, y)]
                if why is None and ff == "dreiding1990" and (x, y) in PYSIMM_DREIDING:
                    why = PYSIMM_DREIDING[(x, y)]
                if why:
                    expl[why].add(name); nexp += 1
                    continue
                good = False
                k = (els[i] if i < len(els) else "?", x, y or "(untyped)")
                diffs[k] += 1
                if len(examples[k]) < 6 and name not in examples[k]:
                    examples[k].append(name)
            mol_ok += good
        rows.append((ff, ref, n, ok, nexp, mol_n, mol_ok, nopar, len(missing)))
        print(f"\n== {ff} against {ref}: {ok}/{n} atoms the same type ({100 * ok / max(n, 1):.2f} %), {nexp} differ for a "
              f"stated reason (the force field's own definition decides), {n - ok - nexp} unexplained; {mol_ok}/{mol_n} molecules whole; "
              f"{nopar} atoms the reference types with a type the library has no parameters for, {len(missing)} molecules the reference cannot type")
        for w, ms in sorted(expl.items()):
            print(f"   explained: {w} ({len(ms)} molecules: {', '.join(sorted(ms)[:8])})")
        for k, v in diffs.most_common(show):
            print(f"   {k[0]:3s} ref {k[1]:10s} CAPS {k[2]:10s} {v:5d}  {', '.join(examples[k])}")
    print("\n| force field | reference | molecules | atoms | same type | reason stated | unexplained | CAPS right % |\n"
          "|---|---|---|---|---|---|---|---|")
    for ff, ref, n, ok, ne, mn, mo, _, _ in rows:
        print(f"| {ff} | {ref} | {mn} | {n} | {ok} | {ne} | {n - ok - ne} | {100 * (ok + ne) / max(n, 1):.2f} |")


if __name__ == "__main__":
    cmd, out = sys.argv[1], os.path.abspath(sys.argv[2])
    {"prepare": prepare, "caps": run_caps, "antechamber": run_antechamber, "foyer": run_foyer, "pysimm": run_pysimm,
     "compare": compare}[cmd](out)
