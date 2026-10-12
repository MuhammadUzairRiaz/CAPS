#!/usr/bin/env python3
"""SiO2 for CAPS from its paper: M. J. Sanders, M. Leslie, C. R. A. Catlow, J. Chem. Soc., Chem. Commun. 1271 (1984), Table 1.

  Si-O   1283.9073 exp(-r/0.3205) - 10.6616/r^6 eV between Si and the O shell
  O-O    22764.0 exp(-r/0.149) - 27.88/r^6 eV between O shells
  Si-Si  Coulomb only
  charges  Si +4 (rigid), O core +0.8482, O shell -2.8482; O core-shell spring 74.9204 eV/A^2 (1/2 k r^2: class II K2 = k/2)
  O-Si-O E_B = K_B (theta - theta0)^2, K_B 2.097 eV/rad^2, theta0 109.47 deg, on the O shells (the paper uses the shell
         positions); every O shell pair within 1.9 A of a Si (CAPS's distance: the paper bends the O-Si-O "bonds" of each
         tetrahedron, Si-O 1.61 A in quartz) as a harmonic angle K = K_B: the paper's form has no 1/2
Si and O are not bonded: Si-O interacts through the Buckingham and Coulomb terms in full.

usage: make_sanders1984.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "M. J. Sanders, M. Leslie, C. R. A. Catlow, J. Chem. Soc., Chem. Commun. 1271 (1984)"
EV = 23.060548
FID = "sanders1984-sio2"
types = [{"name": "Si_s84", "element": "Si", "mass": 28.0855, "charge": 4.0, "description": "Si4+ (rigid ion, formal charge); Table 1"},
         {"name": "O_s84c", "element": "O", "mass": 15.999 * 0.9, "charge": 0.8482, "description": "O core (Table 1: +0.8482)"},
         {"name": "O_s84s", "element": "O", "mass": 15.999 * 0.1, "charge": -2.8482, "description": "O shell (Table 1: -2.8482)"}]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
pairs += [{"name": "Si_s84-O_s84s", "match": ["Si_s84", "O_s84s"], "style": "buck", "params": [1283.9073 * EV, 0.3205, 10.6616 * EV],
           "comment": "Si-O A 1283.9073 eV, B 0.3205 A, C 10.6616 eV A^6 (Table 1)"},
          {"name": "O_s84s-O_s84s", "match": ["O_s84s", "O_s84s"], "style": "buck", "params": [22764.0 * EV, 0.149, 27.88 * EV],
           "comment": "O-O A 22764.000 eV, B 0.149 A, C 27.88 eV A^6 (Table 1)"}]
bonds = [{"name": "O_s84c-O_s84s", "match": ["O_s84c", "O_s84s"], "style": "class2", "params": [0.0, 74.9204 / 2 * EV, 0.0, 0.0],
          "comment": "O core-shell spring 74.9204 eV/A^2 (Table 1)"}]
angles = [{"name": "O_s84s-Si_s84-O_s84s", "match": ["O_s84s", "Si_s84", "O_s84s"], "style": "harmonic", "params": [2.097 * EV, 109.47],
           "contact": [1.9], "comment": "E_B = K_B (theta - theta0)^2, K_B 2.097 eV/rad^2, theta0 109.47 (Table 1 and text); O shells within 1.9 A of Si"}]
rules = [{"type": "Si_s84", "smarts": "[#14]", "priority": 0, "requires": ["O"], "description": "Si in silica"},
         {"type": "O_s84c", "smarts": "[#8]", "priority": 0, "requires": ["Si"], "description": "O in silica"},
         {"type": "O_s84s", "smarts": "[#8]", "priority": 99, "atom_name": "O_s84s", "description": "shell of O_s84c"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "SiO2 (Sanders, Leslie & Catlow 1984)", "source": REF + ", Table 1",
       "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_sanders1984.py) from the paper's Table 1",
                 "O-Si-O bends between O shells within 1.9 A of the Si (the distance is CAPS's); E = K_B (theta - theta0)^2 as the paper writes it",
                 "shells carry 10 % of the oxygen's mass (adiabatic core-shell)"],
       "version": "1984 (Table 1)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined",
           "shells": {"O_s84c": "O_s84s"}, "description": "Si and O of silica (bench/ff/make_sanders1984.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types")
