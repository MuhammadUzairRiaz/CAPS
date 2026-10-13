#!/usr/bin/env python3
"""MgO, CaO, SrO and BaO as rigid ions for CAPS from D. C. Sayle, G. W. Watson, Phys. Chem. Chem. Phys. 2, 5491 (2000),
Table 1.

  V(r) = A exp(-r/rho) - C r^-6, cut-off 10.0 A (Table 1 caption), Coulomb between formal charges M2+ and O2-
  M-O    A, rho of Table 1 (Lewis & Catlow's cation-O terms used "with the additional approximation of the rigid ion model")
  O-O    22764.0 exp(-r/0.1490) - 27.88/r^6 eV (Table 1)
  M-M    Coulomb only (Table 1 lists none)
One set for all four oxides and their interfaces (the paper's CaO / SrO / BaO films on MgO): one O type, one type per
cation, so a mixed system is typed as the paper modelled it.

usage: make_sayle2000.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "D. C. Sayle, G. W. Watson, Phys. Chem. Chem. Phys. 2, 5491 (2000)"
EV = 23.060548   # kcal/mol per eV
FID = "sayle2000-alkaline-earth-oxides"
# Table 1: A (eV), rho (A), C (eV A^6)
T1 = [("Mg", 1428.5, 0.2945, 0.0), ("Ca", 1203.2, 0.3375, 0.0), ("Sr", 1400.0, 0.3500, 0.0), ("Ba", 931.7, 0.3949, 0.0)]
OO = (22764.0, 0.1490, 27.88)
MASS = {"Mg": 24.305, "Ca": 40.078, "Sr": 87.62, "Ba": 137.327, "O": 15.999}
NUM = {"Mg": 12, "Ca": 20, "Sr": 38, "Ba": 56}
types = [{"name": f"{el}_sw", "element": el, "mass": MASS[el], "charge": 2.0, "description": f"{el}2+, rigid ion (formal charge)"} for el, *_ in T1]
types.append({"name": "O_sw", "element": "O", "mass": MASS["O"], "charge": -2.0, "description": "O2-, rigid ion (formal charge)"})
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
pairs += [{"name": f"{el}_sw-O_sw", "match": [f"{el}_sw", "O_sw"], "style": "buck", "params": [A * EV, rho, C * EV],
           "comment": f"{el}2+-O2- A {A} eV, rho {rho} A, C {C} eV A^6 (Table 1)"} for el, A, rho, C in T1]
pairs.append({"name": "O_sw-O_sw", "match": ["O_sw", "O_sw"], "style": "buck", "params": [OO[0] * EV, OO[1], OO[2] * EV],
              "comment": f"O2--O2- A {OO[0]} eV, rho {OO[1]} A, C {OO[2]} eV A^6 (Table 1)"})
rules = [{"type": f"{el}_sw", "smarts": f"[#{NUM[el]}]", "priority": 0, "requires": ["O"], "description": f"{el}2+ in its oxide"} for el, *_ in T1]
rules += [{"type": "O_sw", "smarts": "[#8]", "priority": 0, "requires": [el], "description": f"O2- in an alkaline-earth oxide (with {el})"} for el, *_ in T1]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "MgO, CaO, SrO, BaO rigid ions (Sayle & Watson 2000)", "source": REF + ", Table 1",
       "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10.0,
       "atom_types": types, "pairs": pairs, "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_sayle2000.py) from the paper's Table 1",
                 "rigid ions with formal charges (the paper's rigid-ion approximation of Lewis & Catlow's set); cation-cation Coulomb only",
                 "one set for the four oxides and their interfaces; short-range cut-off 10.0 A (Table 1 caption)"],
       "version": "2000 (Table 1)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined",
           "description": "Mg, Ca, Sr, Ba and O of the rock-salt oxides (bench/ff/make_sayle2000.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries")
