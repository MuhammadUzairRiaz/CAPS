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
    (["CM"], "[CH2X4][CH2X4][CH2X4]", ">CCC<"),
]


def short(t):
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
    q = {"NC": 1.0, "NH": 1.0, "PH": -1.0, "PHE": -1.0}   # sdk_lipids.prm MASS charges (and the README)
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
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    print(f"{out}: {len(names)} names, {len(beads)} fragments")
