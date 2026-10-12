#!/usr/bin/env python3
"""A CAPS force-field file from GROMACS's OPLS-AA force field (share/gromacs/top/oplsaa.ff: Jorgensen's OPLS-AA as
distributed with GROMACS under the LGPL).

  convert_gromacs_opls.py OPLSAA_FF_DIR OUT.json [--name NAME] [--typing ../typing/ID.typing.json]

Read:
  ffnonbonded.itp  [ atomtypes ]: opls_NNN, its bond class, element, mass, charge, sigma (nm), epsilon (kJ/mol)
  ffbonded.itp     [ bondtypes ] (b0 nm, kb kJ/mol/nm^2), [ angletypes ] (theta0, k kJ/mol/rad^2), [ dihedraltypes ]
                   function 3 (Ryckaert-Bellemans C0..C5, kJ/mol); the improper #defines (phase 180, k, n 2)
CAPS's form (real units):
  types      opls_NNN with its own charge; its bond class names it in bonds, angles, dihedrals and impropers
  pairs      epsilon, sigma; geometric mixing; 1-4 x 0.5 (GROMACS [ defaults ]: comb-rule 3, fudgeLJ 0.5, fudgeQQ 0.5)
  bonds      K = kb / 2, r0 (K (r - r0)^2)
  angles     K = k / 2, theta0
  dihedrals  opls V1..V4 from the Ryckaert-Bellemans form (psi = phi - 180): V1 = -2 C1 - 3/2 C3, V2 = -C2 - C4,
             V3 = -C3 / 2, V4 = -C4 / 4
  impropers  cvff K [1 - cos 2 phi] on the planar centres the force field names (carbonyl C, amide N, alkene C, ring C,
             guanidinium C, nitro N), K = k from its #defines; the centre third
The residue-specific torsions GROMACS defines for its amino-acid topologies (dih_SER_THR_chi1 ...) are applied by
residue name there, not by type, and are left out (listed in the notes).
"""
import json, os, re, sys

D, OUT = sys.argv[1], sys.argv[2]
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
KJ = 4.184
esc = lambda t: t.replace("*", "\\*").replace("?", "\\?")   # a literal * in a class name (C*, N*): CAPS's patterns are globs
ELEMENT = {1: "H", 3: "Li", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na", 12: "Mg", 14: "Si", 15: "P", 16: "S", 17: "Cl",
           19: "K", 20: "Ca", 26: "Fe", 29: "Cu", 30: "Zn", 35: "Br", 36: "Kr", 37: "Rb", 53: "I", 54: "Xe", 55: "Cs",
           2: "He", 10: "Ne", 18: "Ar", 38: "Sr", 56: "Ba", 0: "X"}


def sections(path):
    """(section, words, comment) for every data line, the #ifdef HEAVY_H alternatives left out (#else kept)"""
    sec, skip = None, False
    for line in open(path, errors="replace"):
        s = line.strip()
        if s.startswith("#ifdef"):
            skip = True
            continue
        if s.startswith("#else"):
            skip = False
            continue
        if s.startswith("#endif"):
            skip = False
            continue
        if skip or not s or s.startswith(";") or s.startswith("#"):
            continue
        m = re.match(r"\[\s*(\w+)\s*\]", s)
        if m:
            sec = m.group(1)
            continue
        data, _, comment = s.partition(";")
        yield sec, data.split(), comment.strip()


types, pairs, bonds, angles, dihedrals, notes = {}, [], {}, {}, {}, []
for sec, w, c in sections(os.path.join(D, "ffnonbonded.itp")):
    if sec != "atomtypes" or len(w) < 8:
        continue
    name, bt, z, mass, q, sig, eps = w[0], w[1], int(w[2]), float(w[3]), float(w[4]), float(w[6]), float(w[7])
    types[name] = {"name": name, "element": ELEMENT.get(z, "X"), "mass": mass, "charge": q,
                   "description": c, "equivalence": {"bond": bt, "angle": bt, "dihedral": bt, "improper": bt}}
    pairs.append({"name": name, "match": [esc(name)], "params": [eps / KJ, sig * 10]})
for sec, w, c in sections(os.path.join(D, "ffbonded.itp")):
    if sec == "bondtypes" and len(w) >= 5:
        k = tuple(sorted(w[:2]))
        bonds[k] = {"name": "-".join(w[:2]), "match": [esc(x) for x in w[:2]], "style": "harmonic", "params": [float(w[4]) / 2 / KJ / 100, float(w[3]) * 10],
                    **({"comment": c} if c else {})}
    elif sec == "angletypes" and len(w) >= 6:
        k = tuple(w[:3]) if w[0] <= w[2] else tuple(w[2::-1])
        angles[k] = {"name": "-".join(w[:3]), "match": [esc(x) for x in w[:3]], "style": "harmonic", "params": [float(w[5]) / 2 / KJ, float(w[4])],
                     **({"comment": c} if c else {})}
    elif sec == "dihedraltypes" and len(w) >= 11 and w[4] == "3":
        C = [float(x) / KJ for x in w[5:11]]
        if abs(C[5]) > 1e-9:
            notes.append(f"{'-'.join(w[:4])}: C5 = {C[5]} has no OPLS form, left out")
        V = [-2 * C[1] - 1.5 * C[3], -C[2] - C[4], -0.5 * C[3], -0.25 * C[4]]
        names = ["*" if x == "X" else esc(x) for x in w[:4]]
        k = tuple(names) if tuple(names) <= tuple(names[::-1]) else tuple(names[::-1])
        dihedrals[k] = {"name": "-".join(w[:4]), "match": names, "style": "opls", "params": [round(v, 10) for v in V],
                        **({"comment": c} if c else {})}
defines = {m.group(1): [float(x) for x in m.group(2).split()]
           for m in re.finditer(r"#define\s+(improper_\w+)\s+([-\d. ]+)", open(os.path.join(D, "ffbonded.itp")).read())}
residue_torsions = re.findall(r"#define\s+(dih_\w+)", open(os.path.join(D, "ffbonded.itp")).read())
if residue_torsions:
    notes.append(f"{len(residue_torsions)} residue-specific torsions of GROMACS's amino-acid topologies (dih_*) are applied there by "
                 "residue name, not by type: not in this file")

# impropers: the planar centres the force field names in its #defines, each centre third
RING_C = ["CA", "CB", "CN", "CV", "CW", "CR", "CK", "CQ", "CS", "C*"]
imp_rules = []
for define, centres, outer, what in [
        ("improper_Z_CA_X_Y", RING_C, None, "ring carbon (CA, CB, CN, CV, CW, CR, CK, CQ, CS, C*)"),
        ("improper_Z_CM_X_Y", ["CM", "C="], None, "alkene carbon (CM, C=)"),
        ("improper_Z_N_X_Y", ["N", "NA", "N2", "N*", "NC", "NB", "NT"], None, "planar N (amide N, ring N, N2)"),
        ("improper_O_C_X_Y", ["C", "C_2", "C_3"], ["O", "O2", "O_2", "O_3"], "carbonyl / carboxylate C"),
        ("improper_N2_X_N2_N2", ["CA", "C"], ["N2"], "guanidinium C bonded to N2"),
        ("improper_X_NO_ON_NO", ["NO"], ["ON"], "nitro N")]:
    if define not in defines:
        continue
    phase, k, n = defines[define]
    for cen in centres:
        for o in (outer or ["*"]):
            imp_rules.append({"name": f"{define}: {o}-*-{cen}-*" if outer else f"{define}: *-*-{cen}-*",
                              "match": ["*", "*", esc(cen), o] if outer else ["*", "*", esc(cen), "*"], "style": "cvff",
                              "params": [k / KJ, -1 if abs(phase - 180) < 1e-6 else 1, int(n)], "comment": what})
doc = {"format": "caps-forcefield", "format_version": 1, "name": arg("--name", "OPLS-AA (GROMACS oplsaa.ff)"),
       "source": "GROMACS share/top/oplsaa.ff (ffnonbonded.itp, ffbonded.itp)", "units": "real",
       "references": ["W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)",
                      "G. A. Kaminski, R. A. Friesner, J. Tirado-Rives, W. L. Jorgensen, J. Phys. Chem. B 105, 6474 (2001)"],
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"},
       "mixing": "geometric", "special_lj": [0, 0, 0.5], "special_coul": [0, 0, 0.5], "cutoff": 12,
       "improper_order": "center3_sorted", "equivalence": "replace",
       "atom_types": list(types.values()), "pairs": pairs, "bonds": list(bonds.values()), "angles": list(angles.values()),
       "dihedrals": list(dihedrals.values()), "impropers": imp_rules,
       "notes": ["converted by CAPS (bench/ff/convert_gromacs_opls.py) from GROMACS's oplsaa.ff (LGPL)"] + notes,
       "version": "OPLS-AA/L (GROMACS oplsaa.ff)"}
if arg("--typing"):
    doc["typing"] = arg("--typing")
json.dump(doc, open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"{OUT}: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} dihedrals, {len(imp_rules)} impropers")
