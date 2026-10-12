#!/usr/bin/env python3
"""TiO2 for CAPS from its paper: M. Matsui, M. Akaogi, Mol. Simul. 6, 239 (1991), Table 1.

  V_ij = q_i q_j / r - C_i C_j / r^6 + f (B_i + B_j) exp[(A_i + A_j - r) / (B_i + B_j)],  f = 4.184 kJ/(A mol)
as a Buckingham term: A_ij = f (B_i + B_j) exp[(A_i + A_j) / (B_i + B_j)], rho_ij = B_i + B_j, C_ij = C_i C_j
(C in A^3 (kJ/mol)^1/2). Rigid ions with the paper's partial charges (q_O = -q_Ti / 2). All three pairs (Ti-Ti,
Ti-O, O-O) are written: the model is the sum of the per-atom terms for every pair.

usage: make_matsui1991.py [DATA_DIR]
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "M. Matsui, M. Akaogi, Mol. Simul. 6, 239 (1991)"
KCAL = 4.184
FID = "matsui1991-tio2"
# Table 1: q (e), A (A), B (A), C (A^3 (kJ/mol)^1/2)
T1 = {"Ti": (2.196, 1.1823, 0.077, 22.5), "O": (-1.098, 1.6339, 0.117, 54.0)}
MASS = {"Ti": 47.867, "O": 15.999}
types = [{"name": f"{el}_ma", "element": el, "mass": MASS[el], "charge": q, "description": f"{el} in TiO2 (Table 1: q {q}, A {A}, B {B}, C {C})"}
         for el, (q, A, B, C) in T1.items()]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
for a, b in (("Ti", "Ti"), ("Ti", "O"), ("O", "O")):
    rho, A, C = T1[a][2] + T1[b][2], T1[a][1] + T1[b][1], T1[a][3] * T1[b][3]
    pairs.append({"name": f"{a}_ma-{b}_ma", "match": [f"{a}_ma", f"{b}_ma"], "style": "buck",
                  "params": [rho * math.exp(A / rho), rho, C / KCAL],   # f (B_i+B_j) exp(...) in kJ/mol, f/4.184 = 1 kcal/(A mol)
                  "comment": f"f (B_i+B_j) exp[(A_i+A_j)/(B_i+B_j)], B_i+B_j {rho:g} A, C_i C_j {C:g} A^6 kJ/mol (Table 1)"})
rules = [{"type": "Ti_ma", "smarts": "[#22]", "priority": 0, "requires": ["O"], "description": "Ti in TiO2"},
         {"type": "O_ma", "smarts": "[#8]", "priority": 0, "requires": ["Ti"], "description": "O in TiO2"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "TiO2 (Matsui & Akaogi 1991)", "source": REF + ", Table 1", "references": [REF],
       "units": "real", "styles": {"pair": "buck/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_matsui1991.py) from the paper's Table 1",
                 "rigid ions; the paper fitted rutile, anatase and brookite (and simulated TiO2 II)"],
       "version": "1991", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "rules": rules,
           "description": "Ti and O of titanium dioxide (bench/ff/make_matsui1991.py)"},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
for p in pairs[2:]:
    print(p["name"], "A %.2f eV rho %g C %.4f eV A^6" % (p["params"][0] / 23.060548, p["params"][1], p["params"][2] / 23.060548))
