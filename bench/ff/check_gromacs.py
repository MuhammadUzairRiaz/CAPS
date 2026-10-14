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
import template_pdb, mix_box

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
    ("Polystyrene melt, GAFF2, PME (periodic, 1300 atoms)", ("file", PS), "gaff2-2.2.20-ambertools", "gasteiger", PME),
    ("Polystyrene melt, GAFF2, PME, shifted LJ (no tail)", ("file", PS), "gaff2-2.2.20-ambertools", "gasteiger", PME + ["--no-tail"]),
    ("Polystyrene melt, OPLS-AA, PME", ("file", PS), "opls2005", "gasteiger", PME),
    ("Polystyrene melt, OPLS-AA 2024 (its own charges), PME", ("file", PS), "oplsaa2024-moltemplate", "types", PME),
    ("Polystyrene melt, CVFF, PME (cvff impropers)", ("file", PS), "cvff", "types", PME),
    ("Polystyrene melt, GAFF2, PME, triclinic cell (xy 0.5, xz 0.3, yz 0.2 Å)", ("tilt", PS), "gaff2-2.2.20-ambertools", "gasteiger", PME),
    ("(6,6) nanotube, GAFF2, PME (bonds cross the cell: periodic-molecules)", ("nano", ["tube", "--n", "6", "--m", "6", "--length", "25"]),
     "gaff2-2.2.20-ambertools", "gasteiger", PME),
    ("GAFF toluene (vacuum)", ("template", "AMBER16_gaff", "toluene"), "gaff-1.81-ambertools", "gasteiger", []),
    ("CVFF phenol (vacuum, cvff impropers)", ("template", "CVFF", "phenol"), "cvff", "types", []),
    ("OPLS-AA methyl vinyl ketone (vacuum)", ("template", "OPLS2005", "methyl_vinyl_ketone"), "opls2005", "gasteiger", []),
    ("CGenFF toluene (vacuum; separate 1-4 LJ, Urey-Bradley, harmonic impropers)", ("template", "CHARMM36_cgenff", "toluene"), "cgenff", "gasteiger", []),
    # GROMOS in its own settings: reaction field (ε_rf 61) at 1.4 nm, no dispersion correction, C6/C12 (comb-rule 1)
    # ClayFF (Cygan 2004): LAMMPS's ClayFF test structure with its own types, and IFF's Na-montmorillonite typed by CAPS's
    # rules; the M-O-H bends are angles without an M-O bond, so GROMACS keeps M and H interacting as CAPS does. IFF 1.5
    # (CVFF) montmorillonite from IFF's own model database. The cells are thinner than 2 x 10 A: compared at 8 and 9 A.
    ("ClayFF (Cygan 2004) pyrophyllite (M-O-H bends by contact; PME)", ("file", os.path.expanduser("~/lammps/tools/msi2lmp/test/PyAC_bulk-clayff.car")),
     "clayff-cygan2004", "types", PME + ["--cutoff", "8", "--no-tail", "--names"]),
    ("ClayFF (Cygan 2004) Na-montmorillonite, typed by rule (PME)",
     ("file", os.path.expanduser("~/iff-ref/INTERFACE_FF_1_5/MODEL_DATABASE/CLAY_MINERALS/mont0_333_Na_15_cell.car")),
     "clayff-cygan2004", "types", PME + ["--cutoff", "9", "--no-tail"]),
    ("IFF 1.5 (CVFF) Na-montmorillonite (bonded clay layers, cvff impropers; PME)",
     ("file", os.path.expanduser("~/iff-ref/INTERFACE_FF_1_5/MODEL_DATABASE/CLAY_MINERALS/mont0_333_Na_15_cell.car")),
     "iff-cvff", "keep", PME + ["--cutoff", "9", "--no-tail"]),
    # the libraries built from their original sources on this machine (bench/ff/open), on molecules CAPS builds itself
    ("OPLS-AA (GROMACS oplsaa.ff) isoprene, ethylbenzene, ester (vacuum)", ("smiles", "CC=C(C)CC.CCc1ccccc1.CCOC(=O)C"), "oplsaa-gromacs", "types", []),
    ("OPLS-AA 2024 (Jorgensen) isoprene, ethylbenzene, ester, PVA unit (vacuum)", ("smiles", "CC=C(C)CC.CCc1ccccc1.CCOC(=O)C.CC(O)CC(O)C"), "oplsaa2024", "types", []),
    ("CL&P (il.ff) [C4C1im][NTf2] (vacuum)", ("smiles", "CCCCn1cc[n+](C)c1.FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F"), "clandp", "types", []),
    ("CL&P DES (des.ff) reline, 16 choline chloride + 32 urea, PME", ("mix-box", [("C[N+](C)(C)CCO", 16), ("[Cl-]", 16), ("NC(N)=O", 32)], 30,
     os.path.join(FF, "clandp-des.json")), "clandp-des", "types", PME),
    ("CL&P DES (des.ff) glyceline, 12 choline chloride + 24 glycerol, PME", ("mix-box", [("C[N+](C)(C)CCO", 12), ("[Cl-]", 12), ("OCC(O)CO", 24)], 30,
     os.path.join(FF, "clandp-des.json")), "clandp-des", "types", PME),
    ("CL&P DES (des.ff) betaine, P4444 Br, levulinic acid, ethylene glycol (vacuum)", ("smiles", "C[N+](C)(C)CC(=O)[O-].CCCC[P+](CCCC)(CCCC)CCCC.[Br-].CC(=O)CCC(=O)O.OCCO",
     os.path.join(FF, "clandp-des.json")), "clandp-des", "types", []),
    # 38 A boxes for the paper's 16 A cut-off (GROMACS's pair list needs more than twice 1.76 nm), reline at 0.37 g/cm3: a
    # sparser box leaves the van der Waals total small beside the single-precision pair sums it cancels from. Glyceline's
    # Coulomb total is a small rest of large 1-4, short-range and reciprocal sums; GROMACS's own rerun spread (mdrun -nt 1
    # to 4) on 36 + 72 molecules was 0.014 kcal/mol of 314, so it stays at 24 + 48
    ("OPLS-DES (Doherty & Acevedo 2018) reline, 48 choline chloride + 96 urea, PME", ("mix-box", [("C[N+](C)(C)CCO", 48), ("[Cl-]", 48), ("NC(N)=O", 96)], 38,
     os.path.join(FF, "oplsdes-doherty2018.json")), "oplsdes-doherty2018", "types", PME),
    ("OPLS-DES (Doherty & Acevedo 2018) glyceline, 24 choline chloride + 48 glycerol, PME", ("mix-box", [("C[N+](C)(C)CCO", 24), ("[Cl-]", 24), ("OCC(O)CO", 48)], 38,
     os.path.join(FF, "oplsdes-doherty2018.json")), "oplsdes-doherty2018", "types", PME),
    ("OPLS-DES (Doherty & Acevedo 2018) choline chloride, ethylene glycol, levulinic, malonic, oxalic acid, phenol (vacuum)",
     ("smiles", "C[N+](C)(C)CCO.[Cl-].OCCO.CC(=O)CCC(=O)O.OC(=O)CC(=O)O.OC(=O)C(=O)O.Oc1ccccc1", os.path.join(FF, "oplsdes-doherty2018.json")),
     "oplsdes-doherty2018", "types", []),
    ("GAFF 1.81 (gaff.dat) ester, epoxide, amide, thioether (vacuum)", ("smiles", "CCOC(=O)C.C1CO1.CC(=O)NC.CSC"), "gaff-1.81-ambertools", "gasteiger", []),
    ("GAFF2 2.11 (gaff211.dat) isoprene, styrene (vacuum)", ("smiles", "CC=C(C)CC.C=Cc1ccccc1"), "gaff2-2.11-ambertools", "gasteiger", []),
    ("TraPPE-UA (TraPPE database) hexane, DME, methyl methacrylate, ethanol, cis-2-butene (vacuum)", ("smiles", "CCCCCC.COCCOC.CC(=C)C(=O)OC.CCO.C/C=C\\C"),
     "trappe-ua-umn", "types", []),
    # the transferable hydrocarbon terms (Martin & Siepmann 1998, 1999; Wick, Martin & Siepmann 2000) on DP-4 chains, raw
    # from the builder (polyisobutylene, PVME and PVCP with seed 5: seed 3 leaves a contact single-precision GROMACS
    # follows to 5e-5 to 7e-5 only); the charged chains with PME
    ("TraPPE-UA (transferable terms) polypropylene DP 4 (one chain, periodic)", ("grow", "[*]CC(C)[*]", 4, 3), "trappe-ua-umn", "types", []),
    ("TraPPE-UA (transferable terms) cis-1,4-polybutadiene DP 4 (one chain, periodic)", ("grow", "[*]C/C=C\\C[*]", 4, 3), "trappe-ua-umn", "types", []),
    ("TraPPE-UA (transferable terms) trans-1,4-polybutadiene DP 4 (one chain, periodic)", ("grow", "[*]C/C=C/C[*]", 4, 3), "trappe-ua-umn", "types", []),
    ("TraPPE-UA (transferable terms) polystyrene DP 4 (one chain, periodic)", ("grow", "[*]CC([*])C1=CC=CC=C1", 4, 3), "trappe-ua-umn", "types", []),
    ("TraPPE-UA (transferable terms) polyisobutylene DP 4 (one chain, periodic)", ("grow", "*CC(*)(C)C", 4, 5), "trappe-ua-umn", "types", []),
    ("TraPPE-UA (transferable) PVA, poly(vinyl alcohol) DP 4 (one chain, periodic; PME)", ("grow", "[*]CC([*])O", 4, 3), "trappe-ua-umn", "types", PME),
    ("TraPPE-UA (transferable) PVME, poly(vinyl methyl ether) DP 4 (one chain, periodic; PME)", ("grow", "[*]CC([*])OC", 4, 5), "trappe-ua-umn", "types", PME),
    ("TraPPE-UA (transferable) PVCH, poly(vinyl cyclohexane) DP 4 (one chain, periodic)", ("grow", "[*]CC([*])C1CCCCC1", 4, 3), "trappe-ua-umn", "types", []),
    ("TraPPE-UA (transferable) PVCP, poly(vinyl cyclopentane) DP 4 (one chain, periodic)", ("grow", "[*]CC([*])C1CCCC1", 4, 5), "trappe-ua-umn", "types", []),
    ("TraPPE (TraPPE database) CO2, N2, NH3, H2S, O2, ethylene oxide, decane (vacuum; charge sites)", ("smiles", "O=C=O.N#N.N.S.O=O.C1CO1.CCCCCCCCCC"),
     "trappe-ua-umn", "types", []),
    ("TraPPE-EH (TraPPE database) octane, methylamine, acetamide, acetonitrile, nitroethane (vacuum)", ("smiles",
     "CCCCCCCC.CN.CC(N)=O.CC#N.CC[N+](=O)[O-]"), "trappe-eh-umn", "types", []),
    ("TraPPE-EH (TraPPE database) benzene, pyridine, phenol, 2-chlorothiophene, naphthalene (vacuum)", ("smiles",
     "c1ccccc1.c1ccncc1.Oc1ccccc1.Clc1cccs1.c1ccc2ccccc2c1"), "trappe-eh-umn", "types", []),
    ("O'Malley 2013 all-atom alkanes: octane, isobutane, neopentane (vacuum; A/B pairs)", ("smiles", "CCCCCCCC.CC(C)C.CC(C)(C)C"),
     "omalley2013-alkanes", "types", []),
    ("Plant 2007 methanol (vacuum; three methanols)", ("smiles", "CO.CO.CO"), "plant2007-methanol", "types", []),
    ("Peguin 2009 HFA134a FF no. 2 (vacuum; no intramolecular pairs)", ("smiles", "FCC(F)(F)F.FCC(F)(F)F"), "peguin2009-hfa134a", "types", []),
    ("TraPPE-UA (foyer) butanol, 2-propanol, 2,3-dimethylbutane (vacuum)", ("smiles", "CCCCO.CC(C)O.CC(C)C(C)C"), "trappe-ua-foyer", "types", []),
    # OPLS-UA: the 2024 file's united-atom types with the 1984 paper's terms where the file has none
    ("OPLS-UA 2024 + 1984 terms: polypropylene tetramer (vacuum)", ("smiles", "CC(C)CC(C)CC(C)CC(C)C"), "oplsua2024", "types", []),
    ("OPLS-UA 2024 + 1984 terms: 1,4-polybutadiene trimer, cis, trans, cis (vacuum)", ("smiles", "C/C=C\\CC/C=C/CC/C=C\\C"), "oplsua2024", "types", []),
    ("OPLS-UA 2024 + 1984 terms: cis-1,4-polyisoprene trimer (vacuum)", ("smiles", "C/C=C(C)\\CC/C=C(C)\\CC/C=C(C)\\C"), "oplsua2024", "types", []),
    ("OPLS-UA 2024 + 1986 terms: poly(vinyl alcohol) tetramer (vacuum)", ("smiles", "CC(O)CC(O)CC(O)CC(O)C"), "oplsua2024", "types", []),
    ("OPLS-UA 2024 + 1990 terms: polyoxytetramethylene trimer (vacuum)", ("smiles", "CCCCOCCCCOCCCCOCCCC"), "oplsua2024", "types", []),
    ("OPLS-UA 2024 + 1986 / 1990 terms: ethanol, 1-propanol, 2-propanol, tert-butanol, diethyl ether (vacuum)",
     ("smiles", "CCO.CCCO.CC(C)O.CC(C)(C)O.CCOCC"), "oplsua2024", "types", []),
    ("CGenFF (toppar) toluene (vacuum)", ("smiles", "Cc1ccccc1", os.path.join(FF, "cgenff-toppar.json")), "cgenff-toppar", "gasteiger", []),
    ("CGenFF (toppar) N-methylacetamide (vacuum)", ("smiles", "CNC(C)=O", os.path.join(FF, "cgenff-toppar.json")), "cgenff-toppar", "gasteiger", []),
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
    if src[0] == "mix-box":   # a periodic mixture: [(SMILES, count)], cube edge, the force field the molecules are cleaned up with
        return mix_box.build(CAPS, src[1], src[2], os.path.join(work, base + ".xyz"), src[3] if len(src) > 3 else None)
    if src[0] == "smiles":   # a molecule CAPS builds from SMILES
        out = os.path.join(work, base + ".mol2")
        # cleaned up with the force field under test (a raw embedding leaves steep 1-4 contacts that single-precision
        # GROMACS cannot follow to 1e-5)
        cmd = [CAPS, "build", src[1], "-o", out] + (["--ff", src[2]] if len(src) > 2 else [])
        subprocess.run(cmd, capture_output=True, check=True)
        return out
    if src[0] == "grow":   # one chain from CAPS's polymer builder: (grow, repeat-unit SMILES, degree of polymerisation, seed), in a
        # cell of 0.005 g/cm3 (over 2 r_c across, which GROMACS needs)
        out = os.path.join(work, base + ".xyz")
        subprocess.run([CAPS, "grow", "--units", src[1], "--chains", "1", "--dp", str(src[2]), "--density", "0.005", "--seed", str(src[3]),
                        "-o", out], capture_output=True, check=True)
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
    terms = ["Bond", "G96Bond", "Morse", "Quartic-Bonds", "Angle", "G96Angle", "U-B", "Proper-Dih.", "Per.-Imp.-Dih.", "Improper-Dih.", "LJ-14",
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

def force_switch_convention(d, stem, r1, rc):
    """CHARMM's force switch (Steinbach & Brooks, J. Comput. Chem. 15, 667 (1994): CAPS, CHARMM, LAMMPS lj/charmmfsw)
    minus GROMACS's vdw-modifier = Force-switch, in kcal/mol, from the files GROMACS ran. Two conventions differ:
    GROMACS switches each power with its own polynomial (another constant below r1: 0.535 against CHARMM's
    1/(r1³ rc³) = 0.579 nm⁻⁶ for r⁻⁶ at 1.0 / 1.2 nm), and its [pairs] are never switched where CHARMM shifts its 1-4
    pairs as every other pair. Both sums are over the same coefficients; the difference is a constant with small forces."""
    P, at, excl = expand_topology(d, stem)
    excl = {(min(i, j), max(i, j)) for i, j in excl}
    lines = open(os.path.join(d, stem + ".gro")).read().splitlines()
    n = int(lines[1])
    X = [tuple(float(x) * 10 for x in l[20:].split()[:3]) for l in lines[2:2 + n]]
    box = [float(x) * 10 for x in lines[2 + n].split()[:3]]
    def dist(i, j):
        return math.sqrt(sum(((X[i][k] - X[j][k] + box[k] / 2) % box[k] - box[k] / 2) ** 2 for k in range(3)))
    def charmm(r, s, e):
        if r >= rc:
            return 0.0
        if r <= r1:
            v12, v6 = r ** -12 - (r1 * rc) ** -6, r ** -6 - (r1 * rc) ** -3
        else:
            v12 = rc ** 6 / (rc ** 6 - r1 ** 6) * (r ** -6 - rc ** -6) ** 2
            v6 = rc ** 3 / (rc ** 3 - r1 ** 3) * (r ** -3 - rc ** -3) ** 2
        return 4 * e * (s ** 12 * v12 - s ** 6 * v6)
    def gromacs(r, s, e):
        if r >= rc:
            return 0.0
        out = 0.0
        for a, c in ((12, s ** 12), (6, -s ** 6)):
            A = -a * ((a + 4) * rc - (a + 1) * r1) / (rc ** (a + 2) * (rc - r1) ** 2)
            B = a * ((a + 3) * rc - (a + 1) * r1) / (rc ** (a + 2) * (rc - r1) ** 3)
            C = rc ** -a - A / 3 * (rc - r1) ** 3 - B / 4 * (rc - r1) ** 4
            out += 4 * e * c * (r ** -a - C - (A / 3 * (r - r1) ** 3 + B / 4 * (r - r1) ** 4 if r > r1 else 0.0))
        return out
    delta = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) not in excl:
                s, e = P[(at[i], at[j])]
                delta += charmm(dist(i, j), s, e) - gromacs(dist(i, j), s, e)
    # [pairs]: CHARMM shifts them, GROMACS does not
    text, sec, cur, per, counts = [], "", None, {}, []
    for l in open(os.path.join(d, stem + ".top")):
        m = re.match(r'\s*#include\s+"([^"]+)"', l)
        text += open(os.path.join(d, m.group(1))).read().splitlines() if m else [l]
    for raw in text:
        l = raw.split(";")[0].strip()
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        w = l.split()
        if not w or l.startswith("#"):
            continue
        if sec == "moleculetype":
            cur = w[0]
            per[cur] = []
        elif sec == "atoms" and cur:
            per.setdefault(cur + "#n", []).append(1)
        elif sec == "pairs" and cur:
            per[cur].append((int(w[0]) - 1, int(w[1]) - 1, float(w[3]) * 10, float(w[4]) / KJ))
        elif sec == "molecules":
            counts.append((w[0], int(w[1])))
    off = 0
    for name, c in counts:
        size = len(per.get(name + "#n", []))
        for _ in range(c):
            for i, j, s, e in per[name]:
                r = dist(i + off, j + off)
                delta += charmm(r, s, e) - 4 * e * ((s / r) ** 12 - (s / r) ** 6)
            off += size
    return delta


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
    gm = {"bond": G("Bond", "G96Bond", "Morse", "Quartic Bonds"),
          "angle": G("Angle", "G96Angle", "U B"), "dihedral": G("Proper Dih."), "improper": G("Per. Imp. Dih.", "Improper Dih."),
          "vdw": G("LJ 14", "LJ (SR)"), "coulomb": G("Coulomb 14", "Coulomb (SR)", "Coul. recip.")}
    # GROMACS books a periodic improper (function 4) as "Per. Imp. Dih."; CAPS may count function-9 impropers as torsions
    cm = dict(ce)
    cm["vdw"] = ce["vdw"] - tail_caps
    shift14 = pair14_shift(d, "case", rc) if "--no-tail" in extra else 0.0
    cm["vdw"] += shift14
    fsw = 0.0
    if mdp.get("vdw-modifier", "").lower() == "force-switch":
        fsw = force_switch_convention(d, "case", float(mdp["rvdw-switch"]) * 10, rc)
        cm["vdw"] -= fsw
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
    if fsw:
        extra_txt.append(f"CHARMM's force switch {fsw:+.4f} kcal/mol from GROMACS's (another polynomial; [pairs] unswitched), "
                         "taken off CAPS's van der Waals")
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
