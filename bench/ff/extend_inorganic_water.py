#!/usr/bin/env python3
"""CAPS extensions to the library's inorganic and water force fields (documented in each file's notes):

  charges      each type's charge from the force field's compound templates (the converted type tables carry none;
               ionic Buckingham potentials are fitted with these charges: NaCl's ±0.988, MgO's shell-model split ...).
               LaNiO3's oxygen core / shell split (0.389 / -2.389) differs from the other LaBO3 perovskites' (0.24 / -2.24):
               it gets types O6cN / O6sN with O6c / O6s's parameters (vdw and bond equivalences)
  shell masses the shells are massless in the source (relaxed each step there); CAPS integrates them as adiabatic
               shells: 10 % of the ion's mass on the shell, 90 % on the core, as LAMMPS's core/shell examples
  halides      the source writes LiCl's Cl---Cl pair as "Cl3 Cl3" (RbCl's type); it is Cl4's, and RbCl keeps its own
  Hill-Sauer   the bond increments (INORGANIC_zeolite_Hill_Sauer.bci, BOND_INCREMENTS) that give that model its charges
  water        SPC/E, TIP3P (1983) and TIP3P (2004): charges, and the O-H bond and H-O-H angle terms of the moltemplate
               files (rigid models there, kept by SHAKE; CAPS has no constraints and runs them flexible, as the notes say)

usage: extend_inorganic_water.py [DL_LIB_DIR]   (the directory with the source .sf / .bci files)
"""
import collections, glob, json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FF = os.path.join(ROOT, "data", "forcefields")
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/project/dl_f_4.13/lib")
MASS = {"H": 1.008, "Li": 6.94, "Na": 22.99, "K": 39.098, "Rb": 85.468, "Cs": 132.905, "Mg": 24.305, "Ca": 40.078, "Sr": 87.62, "Ba": 137.327,
        "Al": 26.982, "Ga": 69.723, "Si": 28.085, "Ti": 47.867, "Zr": 91.224, "Fe": 55.845, "Cr": 51.996, "Mn": 54.938, "Co": 58.933, "Ni": 58.693,
        "Nb": 92.906, "Ta": 180.948, "Pb": 207.2, "U": 238.029, "La": 138.905, "Pr": 140.908, "Nd": 144.242, "Gd": 157.25, "Eu": 151.964,
        "Tb": 158.925, "Yb": 173.045, "Y": 88.906, "O": 15.999, "Bi": 208.98, "K_": 39.098}
FILES = {"binary_halides": "inorganic-binary-halides", "binary_misc": "inorganic-binary-misc", "binary_oxides": "inorganic-binary-oxides",
         "clay": "inorganic-clay", "glass": "inorganic-glass", "ternary_oxides": "inorganic-ternary-oxides", "zeolite": "inorganic-zeolite",
         "zeolite_Hill_Sauer": "inorganic-zeolite-hill-sauer"}


def note(ff, text):
    ff.setdefault("notes", [])
    if text not in ff["notes"]:
        ff["notes"].append(text)


def template_charges(sf):
    txt = open(sf).read()
    at = re.search(r"^ATOM_TYPE(.*?)^END ATOM_TYPE", txt, re.S | re.M).group(1)
    keys = {l.split()[1] for l in at.splitlines() if len(l.split()) >= 3 and not l.startswith("#")}
    ch = collections.defaultdict(list)
    for m in re.finditer(r"^MOLECULE (\S+).*?\n(.*?)^END MOLECULE", txt, re.S | re.M):
        for l in m.group(2).splitlines():
            w = l.split()
            if len(w) >= 3 and w[0] in keys:
                try:
                    q = float(w[2])
                except ValueError:
                    continue
                ch[w[0]].append((q, m.group(1)))
    return ch


for short, fid in FILES.items():
    p = os.path.join(FF, fid + ".json")
    ff = json.load(open(p))
    types = {t["name"]: t for t in ff["atom_types"]}
    ch = template_charges(os.path.join(SRC, f"INORGANIC_{short}.sf"))
    for name, vals in ch.items():
        if name in types and vals:
            types[name]["charge"] = vals[0][0]   # the first template's (O6c / O6s: the LaBO3 perovskites other than LaNiO3)
    if fid == "inorganic-ternary-oxides":
        for base, q in (("O6c", 0.389), ("O6s", -2.389)):
            if base in types and base + "N" not in types:
                t = dict(types[base], name=base + "N", charge=q, description=types[base].get("description", "") + " (LaNiO3's charge split)",
                         equivalence={"vdw": base, "bond": base, "angle": base})
                ff["atom_types"].append(t)
        note(ff, "CAPS extension: O6cN / O6sN are LaNiO3's oxygen core and shell (charges 0.389 / -2.389), with O6c / O6s's parameters")
    # core-shell pairs: a type Xc and Xs joined by a spring in the file (CLAYFF's obc / obs are two oxygens, not a pair)
    springs = {tuple(sorted(b["match"])) for b in ff.get("bonds", [])}
    is_core = lambda n: n.endswith("c") and n[:-1] + "s" in types and tuple(sorted((n, n[:-1] + "s"))) in springs
    # adiabatic shells: 10 % of the ion's mass on the shell
    shells = 0
    for name, t in list(types.items()):
        if name.endswith("s") and is_core(name[:-1] + "c") and not t.get("mass"):
            el = types[name[:-1] + "c"].get("element") or t.get("element")
            m = MASS.get(el)
            if m:
                t["mass"] = round(0.1 * m, 4)
                types[name[:-1] + "c"]["mass"] = round(0.9 * m, 4)
                shells += 1
    # cores carry charge only: the short-range potential acts on the shells (zero self pair where the file has none)
    zero = 0
    have_pair = {tuple(x["match"]) for x in ff.get("pairs", [])} | {(m,) for x in ff.get("pairs", []) for m in x["match"]}
    for name in list(types):
        if is_core(name) and (name,) not in have_pair:
            ff.setdefault("pairs", []).append({"name": name, "match": [name], "params": [0, 0],
                                               "comment": "CAPS extension: a core has no short-range potential (it acts on the shell)"})
            zero += 1
    if zero:
        note(ff, "CAPS extension: cores get a zero self pair (charge only; the short-range potential acts on the shell)")
    if shells:
        note(ff, "CAPS extension: shells carry 10 % of the ion's mass (adiabatic core-shell, as LAMMPS's core/shell examples); the source's shells are massless")
    if ch:
        note(ff, "CAPS extension: type charges from the force field's compound templates")
    if fid == "inorganic-binary-halides":
        cl3 = [k for k, x in enumerate(ff["pairs"]) if x["match"] == ["Cl3", "Cl3"]]
        if len(cl3) == 2 and not any(x["match"] == ["Cl4", "Cl4"] for x in ff["pairs"]):
            ff["pairs"][cl3[1]]["match"] = ["Cl4", "Cl4"]
            ff["pairs"][cl3[1]]["name"] = "Cl4-Cl4"
            ff["pairs"][cl3[1]]["comment"] = "CAPS correction: the source labels LiCl's Cl---Cl pair Cl3 (RbCl's type); LiCl's chlorine is Cl4"
            note(ff, "CAPS correction: LiCl's Cl---Cl Buckingham pair belongs to Cl4 (the source writes Cl3, RbCl's type, and RbCl's pair was overwritten)")
    if fid == "inorganic-zeolite-hill-sauer" and not ff.get("bond_increments"):
        bci = open(os.path.join(SRC, "INORGANIC_zeolite_Hill_Sauer.bci")).read()
        blk = re.search(r"^BOND_INCREMENTS\n(.*?)^END BOND_INCREMENTS", bci, re.S | re.M).group(1)
        ff["bond_increments"] = [{"name": f"{w[0]}-{w[1]}", "match": [w[0], w[1]], "params": [float(w[2]), float(w[3])]}
                                 for w in (l.split() for l in blk.splitlines()) if len(w) == 4]
        note(ff, "CAPS extension: bond increments (the force field's charges) from its .bci file")
    json.dump(ff, open(p, "w"), ensure_ascii=False, indent=1)
    print(f"{fid}: {sum(1 for t in ff['atom_types'] if 'charge' in t)} of {len(ff['atom_types'])} types with charges; {shells} shells given mass; {zero} zero core pairs")

WATER = {"spce-moltemplate": (-0.8476, 0.4238, 600.0, 1.0, 75.0, 109.47), "tip3p-1983-moltemplate": (-0.834, 0.417, 450.0, 0.9572, 55.0, 104.52),
         "tip3p-2004-moltemplate": (-0.830, 0.415, 450.0, 0.9572, 55.0, 104.52)}
for fid, (qo, qh, kb, r0, ka, a0) in WATER.items():
    p = os.path.join(FF, fid + ".json")
    ff = json.load(open(p))
    for t in ff["atom_types"]:
        t["charge"] = qo if t["name"] == "O" else qh
    if not ff.get("bonds"):
        ff["bonds"] = [{"name": "O-H", "match": ["O", "H"], "style": "harmonic", "params": [kb, r0], "comment": "moltemplate's OH bond (rigid there, by SHAKE)"}]
    if not ff.get("angles"):
        ff["angles"] = [{"name": "H-O-H", "match": ["H", "O", "H"], "style": "harmonic", "params": [ka, a0], "comment": "moltemplate's HOH angle (rigid there, by SHAKE)"}]
    note(ff, "CAPS extension: charges and the O-H / H-O-H terms of the moltemplate file; the model is rigid (SHAKE) there, CAPS has no constraints "
             "and runs it flexible with these terms")
    json.dump(ff, open(p, "w"), ensure_ascii=False, indent=1)
    print(f"{fid}: charges O {qo} H {qh}; bond {kb} {r0}; angle {ka} {a0}")
