#!/usr/bin/env python3
"""La perovskite oxides for CAPS from their papers (shell models, formal charges):

  cherry1995  M. Cherry, M. S. Islam, C. R. A. Catlow, J. Solid State Chem. 118, 125 (1995), Table 1: LaMO3 (M = Cr, Mn, Fe,
              Co, all 3+); O shell Y -2.24, k 42.0
  read1999    M. S. D. Read, M. S. Islam, F. King, F. E. Hancock, J. Phys. Chem. B 103, 1558 (1999), Table 1:
              La2NiO4 (Ni2+ rigid) and its 2+ dopants (Mg, Ca, Sr, Ba, Mn, Fe, Co, Cu); O shell Y -2.389, k 42.0
Both: La-O 1545.21 exp(-r/0.3590), O-O 22764.3 exp(-r/0.1490) - 43.0/r^6 eV between shells; cation-cation Coulomb only;
springs 1/2 k r^2 (class II K2 = k/2); k >= 99999 a rigid ion. Each set typed by composition: a structure with Ni takes
Read's (its cations 2+), one without Cherry's. Short-range cut-off 12 A (Cherry: 3.0 lattice units).

usage: make_la_perovskites.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
EV = 23.060548
CHERRY = "M. Cherry, M. S. Islam, C. R. A. Catlow, J. Solid State Chem. 118, 125 (1995)"
READ = "M. S. D. Read, M. S. Islam, F. King, F. E. Hancock, J. Phys. Chem. B 103, 1558 (1999)"
FID = "la-perovskites"
MASS = {"La": 138.905, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Mg": 24.305, "Ca": 40.078, "Sr": 87.62,
        "Ba": 137.327, "Cu": 63.546, "O": 15.999}
Z = {"La": 57, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Mg": 12, "Ca": 20, "Sr": 38, "Ba": 56, "Cu": 29, "O": 8}
# set -> (reference, O shell Y, O k, [(element, charge, A, rho, Y, k)])
SETS = {"c": (CHERRY, -2.24, 42.0, [("La", 3, 1545.21, 0.3590, -0.25, 145.0), ("Cr", 3, 1690.90, 0.3010, 0.97, 67.0),
                                     ("Mn", 3, 1267.50, 0.3214, 3.00, 95.0), ("Fe", 3, 1156.36, 0.3299, 4.97, 304.7),
                                     ("Co", 3, 1329.82, 0.3087, 2.04, 196.3)]),
        "r": (READ, -2.389, 42.0, [("La", 3, 1545.21, 0.3590, -0.250, 145.0), ("Ni", 2, 641.20, 0.3372, 2.000, 99999.0),
                                   ("Mg", 2, 1428.5, 0.2945, 1.585, 361.6), ("Ca", 2, 1090.4, 0.3437, 3.135, 110.2),
                                   ("Sr", 2, 959.1, 0.3721, 3.251, 71.7), ("Ba", 2, 905.7, 0.3976, 9.203, 459.2),
                                   ("Mn", 2, 1007.4, 0.3262, 3.420, 95.0), ("Fe", 2, 1207.6, 0.3084, 2.997, 62.9),
                                   ("Co", 2, 1491.7, 0.2951, 3.503, 110.5), ("Cu", 2, 3799.3, 0.24273, 2.000, 999999.0)])}
types, pairs, bonds, rules, shells = [], [], [], [], {}
for tag, (ref, Yo, ko, ions) in SETS.items():
    paper = "Cherry 1995" if tag == "c" else "Read 1999"
    oc, os_ = f"O_{tag}c", f"O_{tag}s"
    ts = [{"name": oc, "element": "O", "mass": MASS["O"] * 0.9, "charge": round(-2 - Yo, 6), "description": f"O2- core ({paper}, Table 1)"},
          {"name": os_, "element": "O", "mass": MASS["O"] * 0.1, "charge": Yo, "description": f"O2- shell ({paper}, Table 1)"}]
    shells[oc] = os_
    bonds.append({"name": f"{oc}-{os_}", "match": [oc, os_], "style": "class2", "params": [0.0, ko / 2 * EV, 0.0, 0.0], "comment": f"k {ko} eV/A^2 ({paper})"})
    pairs.append({"name": f"{os_}-{os_}", "match": [os_, os_], "style": "buck", "params": [22764.3 * EV, 0.1490, 43.0 * EV],
                  "comment": f"O-O 22764.3 eV, 0.1490 A, C 43.0 eV A^6 ({paper}, Table 1)"})
    for el, q, A, rho, Y, k in ions:
        rigid = k >= 99999
        c = f"{el}_{tag}" if rigid else f"{el}_{tag}c"
        ts.append({"name": c, "element": el, "mass": MASS[el] * (1 if rigid else 0.9), "charge": float(q) if rigid else round(q - Y, 6),
                   "description": f"{el}{q}+ {'(rigid: k ' + str(k) + ')' if rigid else 'core'} ({paper}, Table 1)"})
        site = c
        if not rigid:
            ts.append({"name": f"{el}_{tag}s", "element": el, "mass": MASS[el] * 0.1, "charge": Y, "description": f"{el}{q}+ shell ({paper})"})
            shells[c] = site = f"{el}_{tag}s"
            bonds.append({"name": f"{c}-{site}", "match": [c, site], "style": "class2", "params": [0.0, k / 2 * EV, 0.0, 0.0], "comment": f"k {k} eV/A^2 ({paper})"})
        pairs.append({"name": f"{site}-{os_}", "match": [site, os_], "style": "buck", "params": [A * EV, rho, 0.0],
                      "comment": f"{el}{q}+-O2- A {A} eV, rho {rho} A ({paper}, Table 1)"})
        rules.append({"type": c, "smarts": f"[#{Z[el]}]", "priority": 0, "requires": ["O"] + (["Ni"] if tag == "r" and el != "Ni" else []),
                      "excludes": ["Ni"] if tag == "c" else [], "description": f"{el}{q}+ ({paper}: {'with' if tag == 'r' else 'no'} Ni)"})
    rules.append({"type": oc, "smarts": "[#8]", "priority": 0, "requires": ["La"] + (["Ni"] if tag == "r" else []), "excludes": ["Ni"] if tag == "c" else [],
                  "description": f"O2- ({paper})"})
    types += ts
    pairs += [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in ts]
for r in rules:
    if not r["excludes"]:
        del r["excludes"]
rules += [{"type": sh, "smarts": f"[#{Z[[t for t in types if t['name'] == sh][0]['element']]}]", "priority": 99, "atom_name": sh,
           "description": f"shell of {c}"} for c, sh in shells.items()]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "La perovskite oxides (Cherry 1995; Read 1999)",
       "source": f"{CHERRY}, Table 1; {READ}, Table 1", "references": [CHERRY, READ], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 12,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_la_perovskites.py) from the papers' Table 1",
                 "a structure with Ni takes Read's set (La2NiO4 and its 2+ dopants), one without Cherry's (LaMO3, M 3+)",
                 "short-range cut-off 12 A (Cherry: 3.0 lattice units)", "shells carry 10 % of the ion's mass (adiabatic core-shell)"],
       "version": "1995, 1999", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": shells,
           "description": "La perovskites by composition (bench/ff/make_la_perovskites.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(bonds)} springs")
