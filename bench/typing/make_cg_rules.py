#!/usr/bin/env python3
"""CAPS typing rules and fixes for the library's coarse-grained force fields (MARTINI 2.0, Dry MARTINI, SDK).

Coarse-grained force fields type beads, not atoms. A bead structure (built by CAPS from the force field's bead
templates, or read from a file whose sites are named by bead type) is typed by site name: each rule gives a bead type to
the sites carrying that name. SDK also documents which atoms each bead stands for (Shinoda, DeVane & Klein 2007, 2010;
the moltemplate distribution's sdk_original_format/README.txt), so an all-atom structure is mapped onto SDK beads by
those fragments: every heavy atom covered by exactly one fragment (hydrogens join their heavy atom), each bead at its
fragment's centre of mass (the typing file's "beads"). MARTINI's sources give bead-level molecules only, no atomistic
mapping, so MARTINI types bead structures.

Force-field settings and fixes (each in the file's notes):
  MARTINI / Dry MARTINI  pair_settings from the source's In Init block: lj/gromacs from 9 to 12 Å, coul/gromacs from
                         1e-6 to 12 Å, relative permittivity 15; the 12 Å cut-off belongs to the model
  Dry MARTINI            the source's pair sigma values are in nm (0.47, 0.6 ...) while its epsilons and bonds are
                         converted: sigma x 10 (the file is marked "not tested" in the source)
  all three              beads have no element (the converter guessed Si, Ca, Al, K, Sc from bead masses)
  SDK                    bead charges from the source's parameter file (sdk_lipids.prm): NC and NH +1, PH and PHE -1
  bead templates         the MARTINI source's molecule templates (EMC bead SMILES) that use only this file's types

usage: make_cg_rules.py [MOLTEMPLATE_FORCE_FIELDS_DIR]
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/moltemplate/moltemplate/force_fields")
MARTINI_SETTINGS = {"lj_inner": 9.0, "coul_inner": 1e-6, "dielectric": 15.0, "model_cutoff": True}

# SDK fragments: sdk_original_format/sdk_lipids.prm gives each bead's chemistry in its MASS section (">" / "<" mark
# the atoms bonded to other beads; hydrogens follow from the bead masses), heavy atoms only as SMARTS
SDK_BEADS = [
    (["NC"], "[CH2X4][CH2X4][NX4+]([CH3X4])([CH3X4])[CH3X4]", ">CCN(C)(C)C, +1 (choline)"),
    (["NH"], "[CH2X4][CH2X4][N;$([NX4+;H3]),$([NX3;H2])]", ">CCN, +1 (ethanolamine)"),
    (["PHE"], "[OX2;$(O[CH2X4][CH2X4][N;$([NX4+;H3]),$([NX3;H2])])][PX4](~[OX1])(~[OX1])~[OX2]", ">OP(O)(O)O<, -1, of a PE lipid (README: PHE -PO4- (PE lipid))"),
    (["PH"], "[OX2][PX4](~[OX1])(~[OX1])~[OX2]", ">OP(O)(O)O<, -1"),
    (["GL"], "[CH2X4;$(*O)][CH1X4;$(*O)][CH2X4;$(*O)]", ">CC<C< (glycerol)"),
    (["EST1", "EST2"], "[CH2X4][CX3](=[OX1])[OX2]", ">CC(=O)O< (the two ester beads, told apart by their bond lengths to GL)"),
    (["CMD2"], "[CH1X3]=[CH1X3]", ">\\C=C/< (cis)"),
    (["CT"], "[CH3X4][CH2X4][CH2X4]", "CCC<"),
    (["CT2"], "[CH3X4][CH2X4]", "CC<"),
    (["EO"], "[CH2X4][OX2][CH2X4]", "-CH2OCH2- (LAMMPS SDK example peg-verlet: 44.05 g/mol)"),
    (["OA"], "[CH2X4][OX2H1]", "-CH2OH (peg-verlet: 31.03 g/mol)"),
    (["SO4"], "[OX2][SX4](~[OX1])(~[OX1])~[OX1]", "-OSO3-, -1 (sds-monolayer: 96.06 g/mol)"),
    (["CM"], "[CH2X4][CH2X4][CH2X4]", ">CCC<"),
]


def short(t):
    if t.get("aliases"):
        return t["aliases"][0]
    m = re.match(r"moltemplate @atom:(\S+)", t.get("description", ""))
    return m.group(1) if m else t["name"]


def note(ff, text):
    ff.setdefault("notes", [])
    if text not in ff["notes"]:
        ff["notes"].append(text)


def emc_templates(files):
    out = {}
    for f in files:
        txt = open(f).read()
        m = re.search(r"^ITEM\tTEMPLATES\n(.*?)^ITEM\tEND", txt, re.S | re.M)
        if not m:
            continue
        for l in m.group(1).splitlines():
            w = l.split()
            if len(w) == 2 and not l.startswith("#"):
                out[w[0]] = w[1]
    return out


def bead_types(smiles):
    return {re.sub(r"[+-]\d*$", "", b) for b in re.findall(r"\[([^\]]+)\]", smiles)}


def bead_rules(types, what):
    """a rule per bead type for sites named by its short name, and one for its full name (a data file CAPS wrote)"""
    out = []
    for t in sorted(types, key=short):
        n = short(t)
        out.append({"type": n, "smarts": "*", "atom_name": n, "priority": 0, "description": f"a {what} bead named {n}"})
        if t["name"] != n:
            out.append({"type": n, "smarts": "*", "atom_name": t["name"], "priority": 0, "description": f"a {what} bead named {t['name']}"})
    return out


def fix_elements(ff):
    bad = sorted({t.get("element") for t in ff["atom_types"] if t.get("element")})
    for t in ff["atom_types"]:
        t["element"] = ""
    if bad:
        note(ff, f"CAPS correction: beads have no element (the converter guessed {', '.join(bad)} from bead masses)")


LMP_EX = os.path.expanduser("~/lammps/examples/PACKAGES/cgspica")


def read_lammps_example(path, names):
    """{types, pairs, bonds, angles} of a LAMMPS SDK / SPICA example data file: names maps type index → bead name; bond and
    angle types are identified by the bead types they join in the example's own topology"""
    import gzip
    txt = gzip.open(path, "rt").read() if path.endswith(".gz") else open(path).read()
    secs, cur = {}, None
    for l in txt.splitlines():
        h = l.split("#")[0].strip()
        if h in ("Masses", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs", "Atoms", "Bonds", "Angles", "Velocities"):
            cur = h
            secs[cur] = []
            continue
        if cur and h:
            secs[cur].append(h.split())
    mass = {int(w[0]): float(w[1]) for w in secs["Masses"]}
    atype = {int(w[0]): int(w[2]) for w in secs["Atoms"]}
    bcoef = {int(w[0]): (float(w[1]), float(w[2])) for w in secs["Bond Coeffs"]}
    acoef = {int(w[0]): (float(w[1]), float(w[2])) for w in secs["Angle Coeffs"]}
    btypes, atypes = {}, {}
    for w in secs.get("Bonds", []):
        key = tuple(sorted((names[atype[int(w[2])]], names[atype[int(w[3])]])))
        btypes.setdefault(int(w[1]), set()).add(key)
    for w in secs.get("Angles", []):
        a, b, c = names[atype[int(w[2])]], names[atype[int(w[3])]], names[atype[int(w[4])]]
        atypes.setdefault(int(w[1]), set()).add(min((a, b, c), (c, b, a)))
    pairs = [(names[int(w[0])], names[int(w[1])], w[2], float(w[3]), float(w[4])) for w in secs.get("PairIJ Coeffs", [])]
    bonds = {k: (bcoef[t],) for t, ks in btypes.items() for k in ks}
    angles = {k: (acoef[t],) for t, ks in atypes.items() for k in ks}
    return {"mass": {names[t]: m for t, m in mass.items()}, "pairs": pairs, "bonds": bonds, "angles": angles}


def read_pair_coeffs(path, names):
    out = []
    for l in open(path):
        w = l.split("#")[0].split()
        if len(w) >= 5 and w[0] == "pair_coeff" and w[1].isdigit():
            out.append((names[int(w[1])], names[int(w[2])], w[3], float(w[4]), float(w[5])))
    return out


def sdk_examples(ff):
    """SDK beads for PEG (OA, EO) and SDS (SO4, SOD) from LAMMPS's SDK / SPICA examples (made by the model's developers):
    their pairs, bonds and angles; CM / CT / CT2 / W terms they share with the library are checked equal"""
    full = {short(t): t["name"] for t in ff["atom_types"]}
    fn = lambda n: full.get(n, f"{n}_b{n}_a{n}_d{n}_i{n}")
    peg = read_lammps_example(os.path.join(LMP_EX, "peg-verlet", "data.pegc12e8.gz"), {1: "OA", 2: "EO", 3: "CM", 4: "CT2", 5: "W"})
    sds_names = {1: "SO4", 2: "CM", 3: "CT", 4: "SOD", 5: "W"}   # the example's own pair_coeff comments
    sds = read_lammps_example(os.path.join(LMP_EX, "sds-monolayer", "data.sds.gz"), sds_names)
    sds["pairs"] = read_pair_coeffs(os.path.join(LMP_EX, "sds-monolayer", "in.sds-regular"), sds_names)
    what = {"OA": ("-CH2OH, the hydroxyl end of a PEG chain", 0.0), "EO": ("-CH2OCH2-, an ethylene-oxide unit", 0.0),
            "SO4": ("-OSO3-, the sulfate of SDS", -1.0), "SOD": ("Na+ with three waters", 1.0)}
    have_pair = {tuple(sorted((short_of(ff, p["match"][0]), short_of(ff, p["match"][-1])))): p for p in ff["pairs"]}
    mismatch, added = [], {"types": 0, "pairs": 0, "bonds": 0, "angles": 0}
    for ex, src in ((peg, "peg-verlet"), (sds, "sds-monolayer")):
        for n, m in ex["mass"].items():
            if n in full:
                continue
            ff["atom_types"].append({"name": fn(n), "element": "", "mass": m, "charge": what[n][1], "aliases": [n],
                                     "description": f"SDK {n}: {what[n][0]} (LAMMPS examples/PACKAGES/cgspica/{src})"})
            full[n] = fn(n)
            added["types"] += 1
        for a, b, style, eps, sig in ex["pairs"]:
            key = tuple(sorted((a, b)))
            if key in have_pair:
                p = have_pair[key]
                if abs(p["params"][0] - eps) > 1e-6 or abs(p["params"][1] - sig) > 1e-4 or p.get("style") != style:
                    mismatch.append(f"pair {a}-{b}")
                continue
            r = {"match": [fn(a)] if a == b else [fn(a), fn(b)], "style": style, "params": [eps, sig], "comment": f"{a}-{b} ({src})"}
            ff["pairs"].append(r)
            have_pair[key] = r
            added["pairs"] += 1
        for (a, b), ((k, r0),) in ex["bonds"].items():
            if any(sorted(short_bond(x) for x in bb["match"]) == sorted((a, b)) for bb in ff["bonds"]):
                old = next(bb for bb in ff["bonds"] if sorted(short_bond(x) for x in bb["match"]) == sorted((a, b)))
                if abs(old["params"][0] - k) > 1e-6 or abs(old["params"][1] - r0) > 1e-4:
                    mismatch.append(f"bond {a}-{b}")
                continue
            ff["bonds"].append({"name": f"{a}-{b}", "match": [f"*_b{a}_a*_d*_i*", f"*_b{b}_a*_d*_i*"], "style": "harmonic", "params": [k, r0],
                                "comment": f"{a}-{b} ({src})"})
            added["bonds"] += 1
        for (a, b, c), ((k, t0),) in ex["angles"].items():
            if any([short_angle(x) for x in aa["match"]] in ([a, b, c], [c, b, a]) for aa in ff["angles"]):
                continue
            ff["angles"].append({"name": f"{a}-{b}-{c}", "match": [f"*_b*_a{a}_d*_i*", f"*_b*_a{b}_d*_i*", f"*_b*_a{c}_d*_i*"], "style": "sdk",
                                 "params": [k, t0], "comment": f"{a}-{b}-{c} ({src})"})
            added["angles"] += 1
    note(ff, "CAPS extension: SDK beads for PEG (OA, EO) and SDS (SO4, SOD) with their pairs, bonds and angles from LAMMPS's SDK / SPICA "
             "examples (examples/PACKAGES/cgspica: peg-verlet, sds-monolayer); their CM, CT, CT2 and W terms equal this file's")
    if mismatch:
        raise SystemExit("the LAMMPS examples differ from the library on " + ", ".join(mismatch))
    return added


def short_of(ff, name):
    for t in ff["atom_types"]:
        if t["name"] == name:
            return short(t) if "aliases" not in t else t["aliases"][0]
    return name


def short_bond(glob):
    m = re.match(r"\*_b(\S+?)_a\*", glob)
    return m.group(1) if m else glob


def short_angle(glob):
    m = re.match(r"\*_b\*_a(\S+?)_d\*", glob)
    return m.group(1) if m else glob


if __name__ == "__main__":
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    entries = {e["id"]: e for e in cat["forcefields"]}
    mfmt = os.path.join(SRC, "martini_original_format")
    for fid, label in (("martini-moltemplate", "MARTINI"), ("drymartini-moltemplate", "Dry MARTINI")):
        p = os.path.join(DATA, "forcefields", fid + ".json")
        ff = json.load(open(p))
        fix_elements(ff)
        ff["pair_settings"] = dict(MARTINI_SETTINGS)
        ff["torsion_terms"] = "if_defined"   # the source: TORSION IGNORE (torsions only where a molecule's file gives them)
        ff["angle_terms"] = "if_defined"     # the source: ANGLE WARN (a lipid has no NC3-PO4-GL1 angle)
        ff["cutoff"] = 12.0
        note(ff, "pair settings from the source (In Init): lj/gromacs 9 to 12 A, coul/gromacs 1e-6 to 12 A, dielectric 15; special_bonds 0 1 1")
        if fid == "martini-moltemplate" and not any(b["match"] == ["*", "*"] for b in ff["bonds"]):
            # martini.prm's BOND_AUTO: any bonded pair the file does not list, 1250 kJ/mol/nm2 at 0.47 nm (first rule: the
            # listed pairs, later, win)
            ff["bonds"].insert(0, {"name": "auto", "match": ["*", "*"], "style": "harmonic", "params": [round(1250 / 2 / 4.184 / 100, 8), 4.7],
                                   "comment": "martini.prm BOND_AUTO: the default bond"})
            note(ff, "the default bond of the source (martini.prm BOND_AUTO: 1250 kJ/mol/nm2, 0.47 nm) for bonded pairs it does not list")
        if fid == "drymartini-moltemplate":
            n = 0
            for r in ff["pairs"]:
                if len(r.get("params", [])) >= 2 and 0 < r["params"][1] < 1.5:
                    r["params"][1] = round(r["params"][1] * 10, 6)
                    n += 1
            if n:
                note(ff, f"CAPS correction: {n} pair sigmas were in nm in the source (its epsilons and bonds are converted); multiplied by 10")
        names = sorted({short(t) for t in ff["atom_types"]})
        rules = bead_rules(ff["atom_types"], label)
        # the source's bead templates whose beads this file types
        tpl = {}
        if fid == "martini-moltemplate":
            # the files the library's MARTINI was converted from (martini.lt: martini, lipids, ions, cholesterol); the
            # other source files' molecules come with their own terms (bench/ff/convert_emc_martini.py overlays)
            all_t = emc_templates([os.path.join(mfmt, f) for f in ("martini.prm", "lipids.prm", "ions.prm", "cholesterol.prm")])
            tpl = {k: v for k, v in sorted(all_t.items()) if bead_types(v) <= set(names)}
            ff["bead_templates"] = tpl
        doc = {"format": "caps-typing", "version": 1, "forcefield": ff["name"], "coarse_grained": True, "unknown_types": "untyped",
               "description": f"CAPS rules for {label} (bench/typing/make_cg_rules.py): bead structures typed by site name "
                              "(bead templates, or files naming their sites by bead type). The sources give no atomistic mapping.",
               "rules": rules}
        out = os.path.join(DATA, "typing", fid + ".typing.json")
        json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
        ff["typing"] = f"../typing/{fid}.typing.json"
        json.dump(ff, open(p, "w"), ensure_ascii=False, indent=1)
        e = entries[fid]
        e["notes"] = f"{label}: lj/gromacs and coul/gromacs (dielectric 15, 12 A cut-off) evaluated by CAPS; bead structures typed by name" + \
                     (f"; {len(tpl)} bead templates" if tpl else "")
        e["typing"] = {"rules": f"typing/{fid}.typing.json", "evidence": f"CAPS rules (bench/typing/make_cg_rules.py): {len(rules)} bead names"}
        print(f"{out}: {len(rules)} rules; {len(tpl)} templates")

    # SDK
    fid = "sdk-moltemplate"
    p = os.path.join(DATA, "forcefields", fid + ".json")
    ff = json.load(open(p))
    fix_elements(ff)
    ff["torsion_terms"] = "if_defined"   # sdk_lipids.prm: TORSION IGNORE
    ff["angle_terms"] = "if_defined"     # sdk_lipids.prm: ANGLE WARN
    # the model's cut-off: 15 A in the LAMMPS SDK / SPICA examples (examples/PACKAGES/cgspica: lj/sdk 15.0); the moltemplate
    # file's 12 A leaves water at 0.85 g/cm3 (15 A: 1.03), the 12-4 water tail being long
    ff["cutoff"] = 15.0
    ff["pair_settings"] = {"model_cutoff": True}
    note(ff, "cut-off 15 A, the model's (LAMMPS SDK / SPICA examples: lj/sdk 15.0); the source file's 12 A leaves W water at 0.85 g/cm3")
    # SDK's electrostatics are screened: LAMMPS's SDS example gives its ions +-0.1118 e = 1/sqrt(80), i.e. unit charges with
    # a relative permittivity of 80
    ff["pair_settings"]["dielectric"] = 80.0
    note(ff, "relative permittivity 80 (LAMMPS's SDK SDS example: ionic beads +-0.1118 e = 1/sqrt(80))")
    added = sdk_examples(ff)
    print(f"SDK from the LAMMPS examples: {added}")
    # sdk_lipids.prm MASS charges (and the README); SO4 / SOD: the SDS example's +-0.1118 e times sqrt(80)
    q = {"NC": 1.0, "NH": 1.0, "PH": -1.0, "PHE": -1.0, "SO4": -1.0, "SOD": 1.0}
    for t in ff["atom_types"]:
        t["charge"] = q.get(short(t), 0.0)
    note(ff, "bead charges from the source's parameter file (sdk_lipids.prm: NC, NH +1; PH, PHE -1; others 0)")
    for t in ff["atom_types"]:
        if short(t) == "PH":
            t["equivalence"] = {"bond": "PHE_bPHE_aPHE_dPHE_iPHE", "angle": "PHE_bPHE_aPHE_dPHE_iPHE"}
    note(ff, "CAPS extension: PH takes PHE's bond and angle terms (the source gives the phosphate's bonded terms under PHE only, "
             "including the PC head's NC-PHE bond and GL-PHE-NC angle); PH keeps its own non-bonded pairs")
    names = sorted({short(t) for t in ff["atom_types"]})
    beads = [{("types" if len(t) > 1 else "type"): (t if len(t) > 1 else t[0]), "smarts": sm, "description": d} for t, sm, d in SDK_BEADS]
    water = {"type": "W", "molecule": "[OH2]", "count": 3, "description": "(H2O)3: three nearest waters"}
    doc = {"format": "caps-typing", "version": 1, "forcefield": ff["name"], "coarse_grained": True, "unknown_types": "untyped",
           "description": "CAPS rules for SDK (bench/typing/make_cg_rules.py): bead structures typed by site name; all-atom structures "
                          "mapped onto SDK beads by the source's fragments (every heavy atom in exactly one bead, hydrogens with their atom).",
           "rules": bead_rules(ff["atom_types"], "SDK"), "beads": beads, "bead_groups": [water]}
    out = os.path.join(DATA, "typing", fid + ".typing.json")
    json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
    ff["typing"] = f"../typing/{fid}.typing.json"
    json.dump(ff, open(p, "w"), ensure_ascii=False, indent=1)
    e = entries[fid]
    e["notes"] = "SDK: lj/sdk 9-6 / 12-4 and the sdk angle's 1-3 repulsion evaluated by CAPS; all-atom lipids, alkanes and water mapped onto beads"
    e["typing"] = {"rules": f"typing/{fid}.typing.json", "evidence": f"CAPS rules (bench/typing/make_cg_rules.py): {len(names)} bead names, {len(beads)} fragments, water triplets"}
    # Cooke-Deserno (cooke_deserno_lipid.lt and its supporting gen_potential-cooke.py): the tabulated pair is LAMMPS's
    # cosine/squared (WCA with b = 0.95 sigma between heads and heads / tails, WCA plus a cos^2 well of width w_c = 1.5
    # sigma between tails, the table the .lt names: TAIL_TAIL_Wc_1.5); FENE H-T and T-T bonds, a harmonic H-T bond
    # across the lipid. CAPS types bonds by type pair, so the lipid's end tail is TE (T's non-bonded terms)
    fid = "cooke-deserno-moltemplate"
    p = os.path.join(DATA, "forcefields", fid + ".json")
    ff = json.load(open(p))
    rm = 2 ** (1 / 6)
    ff["atom_types"] = [{"name": "H", "element": "", "mass": 1.0, "charge": 0.0, "description": "head bead"},
                        {"name": "T", "element": "", "mass": 1.0, "charge": 0.0, "description": "tail bead"},
                        {"name": "TE", "element": "", "mass": 1.0, "charge": 0.0, "description": "CAPS extension: the lipid's end tail bead (T's non-bonded terms)",
                         "equivalence": {"vdw": "T"}}]
    ff["styles"]["pair"] = "cosine/squared"
    ff["pairs"] = [{"name": "H-H", "match": ["H"], "style": "cosine/squared", "params": [1.0, round(0.95 * rm, 10), round(0.95 * rm, 10), 1], "comment": "WCA, b = 0.95"},
                   {"name": "H-T", "match": ["H", "T"], "style": "cosine/squared", "params": [1.0, round(0.95 * rm, 10), round(0.95 * rm, 10), 1], "comment": "WCA, b = 0.95"},
                   {"name": "T-T", "match": ["T"], "style": "cosine/squared", "params": [1.0, round(rm, 10), round(rm + 1.5, 10), 1],
                    "comment": "WCA, b = 1, plus the cos^2 attraction of width 1.5"}]
    ff["bonds"] = [{"name": "H-T", "match": ["H", "T"], "style": "fene", "params": [30.0, 1.5, 0.0, 0.0], "comment": "fene 30 1.5 0 0"},
                   {"name": "T-TE", "match": ["T", "TE"], "style": "fene", "params": [30.0, 1.5, 0.0, 0.0], "comment": "fene 30 1.5 0 0"},
                   {"name": "H-TE", "match": ["H", "TE"], "style": "harmonic", "params": [10.0, 4.0], "comment": "harmonic 10 4 (the lipid's straightening bond)"}]
    ff["angles"], ff["dihedrals"], ff["impropers"] = [], [], []
    ff["cutoff"] = 3.0
    ff["pair_settings"] = {"model_cutoff": True}
    ff["angle_terms"] = ff["torsion_terms"] = "if_defined"
    ff["special_lj"], ff["special_coul"] = [0, 1, 1], [0, 1, 1]   # the source: special_bonds lj 0 1 1 (no charges: Coulomb alike)
    ff["units"] = "lj"
    ff["bead_templates"] = {"lipid": "[H]1[T][TE]1"}
    ff["typing"] = f"../typing/{fid}.typing.json"
    ff["notes"] = [n for n in ff.get("notes", []) if "converted from moltemplate" in n] + [
        "pairs: the source's table (gen_potential-cooke.py, w_c = 1.5) as LAMMPS cosine/squared; bonds FENE 30 1.5 and harmonic 10 4",
        "CAPS extension: TE, the lipid's end tail (the source bonds the head to both tails, FENE to one and harmonic to the other)",
        "reduced units (sigma = epsilon = m = 1): CAPS reads the numbers as A, kcal/mol, g/mol; scale lengths, energies and times yourself"]
    json.dump(ff, open(p, "w"), ensure_ascii=False, indent=1)
    doc = {"format": "caps-typing", "version": 1, "forcefield": ff["name"], "coarse_grained": True, "unknown_types": "untyped",
           "description": "CAPS rules for the Cooke-Deserno lipid (bench/typing/make_cg_rules.py): beads typed by name (H, T, TE)",
           "rules": [{"type": n, "smarts": "*", "atom_name": n, "priority": 0, "description": f"a bead named {n}"} for n in ("H", "T", "TE")]}
    json.dump(doc, open(os.path.join(DATA, "typing", fid + ".typing.json"), "w"), ensure_ascii=False, indent=1)
    e = entries[fid]
    e["notes"] = "Cooke-Deserno lipid: cosine/squared pairs (w_c 1.5), FENE and harmonic bonds, evaluated by CAPS; reduced units"
    e["typing"] = {"rules": f"typing/{fid}.typing.json", "evidence": "CAPS rules (bench/typing/make_cg_rules.py): bead names H, T, TE"}
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    print(f"{out}: {len(names)} names, {len(beads)} fragments; Cooke-Deserno rewritten")
