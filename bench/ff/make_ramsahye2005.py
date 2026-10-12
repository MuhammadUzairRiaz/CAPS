#!/usr/bin/env python3
"""Zeolite NaY with chloroform for CAPS from its paper: N. A. Ramsahye, R. G. Bell, J. Phys. Chem. B 109, 4738 (2005), Table 1.

  framework  rigid ions with Kramer's charges (Si 2.4, Al 1.4, O -1.2, Na 1.0); Buckingham Si-O, Al-O, O-O, Na-O (every other
             pair Coulomb only); bends 1/2 k (theta - theta0)^2 on O-Si-O (k 12.1 eV/rad^2) and O-Al-O (2.2), each O within
             2.0 A of the T (CAPS's distance: first neighbours); the Si-O and Al-O bonds carry no term (kept for typing only)
  chloroform charges H 0.180, C -0.102, Cl -0.026; bonds 1/2 k (r - r0)^2 (C-Cl, C-H) and bends 1/2 k (theta - theta0)^2
             (Cl-C-Cl, H-C-Cl); Lennard-Jones 4 eps [(sigma/r)^12 - (sigma/r)^6] for the pairs Table 1 lists (chloroform-chloroform,
             chloroform-O, chloroform-Na); the rest Coulomb only
  cut-off 10 A (the paper gives none)

usage: make_ramsahye2005.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "N. A. Ramsahye, R. G. Bell, J. Phys. Chem. B 109, 4738 (2005)"
EV = 23.060548
FID = "ramsahye2005-nay"
types = [{"name": "Si_rb", "element": "Si", "mass": 28.0855, "charge": 2.4, "description": "framework Si (Table 1)"},
         {"name": "Al_rb", "element": "Al", "mass": 26.982, "charge": 1.4, "description": "framework Al (Table 1)"},
         {"name": "O_rb", "element": "O", "mass": 15.999, "charge": -1.2, "description": "framework O (Table 1)"},
         {"name": "Na_rb", "element": "Na", "mass": 22.98977, "charge": 1.0, "description": "extra-framework Na (Table 1)"},
         {"name": "C_cf", "element": "C", "mass": 12.011, "charge": -0.102, "description": "chloroform C (Table 1)"},
         {"name": "H_cf", "element": "H", "mass": 1.008, "charge": 0.180, "description": "chloroform H (Table 1)"},
         {"name": "Cl_cf", "element": "Cl", "mass": 35.453, "charge": -0.026, "description": "chloroform Cl (Table 1)"}]
BUCK = [("Si_rb", "O_rb", 30023.0, 0.1621, 12.840), ("Al_rb", "O_rb", 26998.0, 0.1622, 12.840), ("O_rb", "O_rb", 894.6, 0.3244, 0.0),
        ("Na_rb", "O_rb", 8200.0, 0.218, 11.8)]
LJ = [("C_cf", "O_rb", 0.750232e-2, 3.25), ("H_cf", "O_rb", 0.780134e-2, 2.70), ("Cl_cf", "O_rb", 0.142532e-1, 3.43),
      ("C_cf", "C_cf", 0.222836e-2, 3.75), ("H_cf", "C_cf", 0.230329e-2, 3.36), ("C_cf", "Cl_cf", 0.479530e-2, 3.79),
      ("H_cf", "Cl_cf", 0.495729e-2, 3.39), ("Na_rb", "Cl_cf", 0.183103e-1, 2.90), ("Cl_cf", "Cl_cf", 0.103230e-1, 3.82),
      ("H_cf", "Na_rb", 0.983246e-3, 3.10), ("C_cf", "Na_rb", 0.114095e-2, 3.69)]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types if t["name"] not in ("C_cf", "Cl_cf")]
pairs += [{"name": f"{a}-{b}", "match": [a, b], "style": "buck", "params": [A * EV, rho, C * EV], "comment": f"A {A} eV, rho {rho} A, C {C} eV A^6 (Table 1)"}
          for a, b, A, rho, C in BUCK]
pairs += [{"name": a if a == b else f"{a}-{b}", "match": [a] if a == b else [a, b], "params": [e * EV, s], "comment": f"eps {e:g} eV, sigma {s} A (Table 1)"}
          for a, b, e, s in LJ]
listed = {tuple(sorted(p["match"] * (2 if len(p["match"]) == 1 else 1))) for p in pairs}
names = [t["name"] for t in types]
for i, a in enumerate(names):
    for b in names[i + 1:]:
        if tuple(sorted((a, b))) not in listed:
            pairs.append({"name": f"{a}-{b}", "match": [a, b], "params": [0.0, 0.0], "comment": "not in Table 1: Coulomb only"})
bonds = [{"name": "C_cf-Cl_cf", "match": ["C_cf", "Cl_cf"], "style": "harmonic", "params": [27.2511 / 2 * EV, 1.761], "comment": "1/2 k (r - r0)^2, k 27.2511 eV/A^2 (Table 1)"},
         {"name": "C_cf-H_cf", "match": ["C_cf", "H_cf"], "style": "harmonic", "params": [29.5612 / 2 * EV, 1.105], "comment": "k 29.5612 eV/A^2 (Table 1)"}]
angles = [{"name": "Cl_cf-C_cf-Cl_cf", "match": ["Cl_cf", "C_cf", "Cl_cf"], "style": "harmonic", "params": [4.33652 / 2 * EV, 109.50], "comment": "k 4.33652 eV/rad^2 (Table 1)"},
          {"name": "H_cf-C_cf-Cl_cf", "match": ["H_cf", "C_cf", "Cl_cf"], "style": "harmonic", "params": [5.81095 / 2 * EV, 107.10], "comment": "k 5.81095 eV/rad^2 (Table 1)"},
          {"name": "O_rb-Si_rb-O_rb", "match": ["O_rb", "Si_rb", "O_rb"], "style": "harmonic", "params": [12.100 / 2 * EV, 109.47], "contact": [2.0],
           "comment": "1/2 k (theta - theta0)^2, k 12.100 eV/rad^2 (Table 1); O within 2.0 A of the Si"},
          {"name": "O_rb-Al_rb-O_rb", "match": ["O_rb", "Al_rb", "O_rb"], "style": "harmonic", "params": [2.200 / 2 * EV, 109.47], "contact": [2.0],
           "comment": "k 2.200 eV/rad^2 (Table 1); O within 2.0 A of the Al"}]
CHCL3 = "[#6;$([#6]([#1])(Cl)(Cl)Cl)]"
rules = [{"type": "Si_rb", "smarts": "[#14]", "priority": 0, "requires": ["O"], "description": "framework Si"},
         {"type": "Al_rb", "smarts": "[#13]", "priority": 0, "requires": ["O"], "description": "framework Al"},
         {"type": "O_rb", "smarts": "[#8]", "priority": 0, "description": "framework O"},
         {"type": "Na_rb", "smarts": "[#11]", "priority": 0, "description": "extra-framework Na"},
         {"type": "C_cf", "smarts": CHCL3, "priority": 0, "description": "chloroform C"},
         {"type": "H_cf", "smarts": f"[#1;$([#1]{CHCL3})]", "priority": 0, "description": "chloroform H"},
         {"type": "Cl_cf", "smarts": f"[Cl;$(Cl{CHCL3})]", "priority": 0, "description": "chloroform Cl"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Zeolite NaY and chloroform (Ramsahye & Bell 2005)", "source": REF + ", Table 1",
       "references": [REF], "units": "real",
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_ramsahye2005.py) from the paper's Table 1",
                 "framework: no T-O bond terms; O-T-O bends by distance (2.0 A, CAPS's)", "every pair Table 1 does not list: Coulomb only",
                 "cut-off 10 A (the paper gives none)"],
       "version": "2005 (Table 1)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined",
           "description": "NaY framework, Na and chloroform (bench/ff/make_ramsahye2005.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries")
