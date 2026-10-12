#!/usr/bin/env python3
"""Silicalite with benzene for CAPS from its paper: T. R. Forester, W. Smith, J. Chem. Soc., Faraday Trans. 93, 3249 (1997),
Tables 1 and 2 and eq. 14 (Schroder and Sauer's shell model of silica, as the paper gives it).

  framework  Si 4.0 (rigid); O core 1.06237 (13.9994 u) / shell -3.06237 (2 u), spring 1/2 k_shell r^2, k_shell 10879.96 kJ/mol/A^2;
             A exp(-r/rho) between Si and the O shell only (A 149643.9 kJ/mol, rho 0.30017 A); every other pair Coulomb only;
             bends 1/2 k_theta (theta - theta0)^2 on O shell-Si-O shell, k_theta 17.7504 kJ/mol/rad^2, theta0 109.47, both shells
             in the first coordination sphere of the Si (within 2.0 A: CAPS's distance)
  benzene    C -0.075, H 0.075 (Table 2); A/r^12 - B/r^6 with the O shells (C-O: A 1.335e6 kJ/mol A^12, B 1.807e3 kJ/mol A^6;
             H-O: 1.599e5, 5.108e2), as Lennard-Jones sigma = (A/B)^(1/6), eps = B^2 / 4A; the paper's benzene is rigid (C-C
             1.4 A, C-H 1.09 A): no intramolecular terms here, run it as a rigid body
  short-range terms cut at 9.8 A as the paper

usage: make_forester1997.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "T. R. Forester, W. Smith, J. Chem. Soc., Faraday Trans. 93, 3249 (1997)"
KJ = 4.184
FID = "forester1997-silicalite"
types = [{"name": "Si_fs", "element": "Si", "mass": 28.0855, "charge": 4.0, "description": "Si (Table 1)"},
         {"name": "O_fsc", "element": "O", "mass": 13.9994, "charge": 1.06237, "description": "O core (Table 1; 13.9994 u as the paper)"},
         {"name": "O_fss", "element": "O", "mass": 2.0, "charge": -3.06237, "description": "O shell (Table 1; 2 u as the paper)"},
         {"name": "C_bz", "element": "C", "mass": 12.011, "charge": -0.075, "description": "benzene C (Table 2)"},
         {"name": "H_bz", "element": "H", "mass": 1.008, "charge": 0.075, "description": "benzene H (Table 2)"}]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
pairs.append({"name": "Si_fs-O_fss", "match": ["Si_fs", "O_fss"], "style": "buck", "params": [149643.9 / KJ, 0.30017, 0.0],
              "comment": "A 149643.9 kJ/mol, rho 0.30017 A (Table 1)"})
for a, A, Bc in (("C_bz", 1.335e6, 1.807e3), ("H_bz", 1.599e5, 5.108e2)):
    pairs.append({"name": f"{a}-O_fss", "match": [a, "O_fss"], "params": [Bc * Bc / (4 * A) / KJ, (A / Bc) ** (1 / 6)],
                  "comment": f"A {A:g} kJ/mol A^12, B {Bc:g} kJ/mol A^6 (Table 2)"})
names = [t["name"] for t in types]
listed = {("Si_fs", "O_fss"), ("C_bz", "O_fss"), ("H_bz", "O_fss")}
for i, a in enumerate(names):
    for b in names[i + 1:]:
        if (a, b) not in listed and (b, a) not in listed:
            pairs.append({"name": f"{a}-{b}", "match": [a, b], "params": [0.0, 0.0], "comment": "Coulomb only"})
bonds = [{"name": "O_fsc-O_fss", "match": ["O_fsc", "O_fss"], "style": "class2", "params": [0.0, 10879.96 / 2 / KJ, 0.0, 0.0],
          "comment": "1/2 k_shell r^2, k_shell 10879.96 kJ/mol/A^2 (Table 1, eq. 14)"}]
angles = [{"name": "O_fss-Si_fs-O_fss", "match": ["O_fss", "Si_fs", "O_fss"], "style": "harmonic", "params": [17.7504 / 2 / KJ, 109.47], "contact": [2.0],
           "comment": "1/2 k_theta (theta - theta0)^2, k_theta 17.7504 kJ/mol/rad^2 (Table 1, eq. 14); shells within 2.0 A of the Si"}]
rules = [{"type": "Si_fs", "smarts": "[#14]", "priority": 0, "requires": ["O"], "description": "Si of silicalite"},
         {"type": "O_fsc", "smarts": "[#8]", "priority": 0, "description": "O core"},
         {"type": "C_bz", "smarts": "[c;$(c1ccccc1)]", "priority": 0, "description": "benzene C"},
         {"type": "H_bz", "smarts": "[#1;$([#1]c)]", "priority": 0, "description": "benzene H"},
         {"type": "O_fss", "smarts": "[#8]", "priority": 99, "atom_name": "O_fss", "description": "shell of O_fsc"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Silicalite and benzene (Forester & Smith 1997)", "source": REF + ", Tables 1 and 2",
       "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 9.8,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_forester1997.py) from the paper's Tables 1 and 2 and eq. 14",
                 "the paper's benzene is rigid: no intramolecular terms (run it as a rigid body)",
                 "O-Si-O bends on shells within 2.0 A of the Si (the distance is CAPS's)", "short-range cut-off 9.8 A as the paper"],
       "version": "1997 (Tables 1, 2)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": {"O_fsc": "O_fss"},
           "description": "silicalite Si and O, benzene (bench/ff/make_forester1997.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries")
