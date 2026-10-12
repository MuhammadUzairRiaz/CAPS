#!/usr/bin/env python3
"""Zeolites with bridging hydroxyls for CAPS from their paper: K.-P. Schroder, J. Sauer, M. Leslie, C. R. A. Catlow,
J. M. Thomas, Chem. Phys. Lett. 188, 320 (1992), Table 1 and eq. 1-2 (read from the page image).

  charges   Si 4.0, Al 3.0 (rigid); framework O core 0.86902 / shell -2.86902, spring 74.92 eV/A^2 (1/2 k r^2: class II
            K2 = k/2); hydroxyl O -1.426 and H 0.426 (rigid)
  pairs     Buckingham A exp(-r/rho) - C/r^6: Si-O, Al-O, O-O (O the framework shell), Si-OH, Al-OH (refitted in the paper),
            O-OH, O-H; every other pair Coulomb only
  O-H       Morse D [(1 - e^(-a(r - r0)))^2 - 1], D 7.0525 eV, a 2.1986 / A, r0 0.9485 A, as a pair term within 1.4 A (the cut-off
            of the published GULP transcriptions; the paper gives none) with the O-H Coulomb term kept, as in the paper's lattice code
  bends     k (theta - theta0)^2 (eq. 2: no 1/2), k 2.09724 eV/rad^2, theta0 109.47, on O-T-O (T Si or Al; O the framework shell or
            the hydroxyl O), each O within 2.2 A of the T (CAPS's distance: first neighbours only)
  short-range cut-off 7.5 A as the paper
The O-H bond of a hydroxyl carries no bond term: kept for typing only.

usage: make_schroder1992.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "K.-P. Schroder, J. Sauer, M. Leslie, C. R. A. Catlow, J. M. Thomas, Chem. Phys. Lett. 188, 320 (1992)"
EV = 23.060548
FID = "schroder1992-zeolite"
types = [{"name": "Si_sz", "element": "Si", "mass": 28.0855, "charge": 4.0, "description": "Si4+ (Table 1)"},
         {"name": "Al_sz", "element": "Al", "mass": 26.982, "charge": 3.0, "description": "Al3+ (Table 1)"},
         {"name": "O_szc", "element": "O", "mass": 15.999 * 0.9, "charge": 0.86902, "description": "framework O2- core (Table 1)"},
         {"name": "O_szs", "element": "O", "mass": 15.999 * 0.1, "charge": -2.86902, "description": "framework O2- shell (Table 1)"},
         {"name": "Oh_sz", "element": "O", "mass": 15.999, "charge": -1.426, "description": "hydroxyl O1.4- (Table 1)"},
         {"name": "H_sz", "element": "H", "mass": 1.008, "charge": 0.426, "description": "hydroxyl H0.4+ (Table 1)"}]
B = [("Si_sz", "O_szs", 1283.907, 0.32052, 10.66158), ("Al_sz", "O_szs", 1460.3, 0.29912, 0.0), ("O_szs", "O_szs", 22764.0, 0.149, 27.88),
     ("Si_sz", "Oh_sz", 983.5566, 0.32052, 10.66158), ("Al_sz", "Oh_sz", 1142.6775, 0.29912, 0.0), ("O_szs", "Oh_sz", 22764.0, 0.149, 27.88),
     ("O_szs", "H_sz", 311.97, 0.25, 0.0)]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
pairs += [{"name": f"{a}-{b}", "match": [a, b], "style": "buck", "params": [A * EV, rho, C * EV], "comment": f"A {A} eV, rho {rho} A, C {C} eV A^6 (Table 1)"}
          for a, b, A, rho, C in B]
pairs.append({"name": "Oh_sz-H_sz", "match": ["Oh_sz", "H_sz"], "style": "morse", "params": [7.0525 * EV, 2.1986, 0.9485, 1.4],
              "comment": "Morse D 7.0525 eV, a 2.1986 / A, r0 0.9485 A (Table 1), within 1.4 A"})
bonds = [{"name": "O_szc-O_szs", "match": ["O_szc", "O_szs"], "style": "class2", "params": [0.0, 74.92 / 2 * EV, 0.0, 0.0], "comment": "core-shell k 74.92 eV/A^2 (Table 1)"}]
angles = []
for T in ("Si_sz", "Al_sz"):
    for a, b in (("O_szs", "O_szs"), ("O_szs", "Oh_sz"), ("Oh_sz", "Oh_sz")):
        angles.append({"name": f"{a}-{T}-{b}", "match": [a, T, b], "style": "harmonic", "params": [2.09724 * EV, 109.47], "contact": [2.2],
                       "comment": "k (theta - theta0)^2, k 2.09724 eV/rad^2, theta0 109.47 (Table 1, eq. 2); O within 2.2 A of the T"})
OH = "[#8;$([#8][#1])]"
rules = [{"type": "Si_sz", "smarts": "[#14]", "priority": 0, "requires": ["O"], "description": "Si4+"},
         {"type": "Al_sz", "smarts": "[#13]", "priority": 0, "requires": ["O"], "description": "Al3+"},
         {"type": "O_szc", "smarts": "[#8]", "priority": 0, "description": "framework O2- core"},
         {"type": "Oh_sz", "smarts": OH, "priority": 1, "description": "hydroxyl O (bonded to an H)"},
         {"type": "H_sz", "smarts": "[#1;$([#1][#8])]", "priority": 0, "description": "hydroxyl H"},
         {"type": "O_szs", "smarts": "[#8]", "priority": 99, "atom_name": "O_szs", "description": "shell of O_szc"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Zeolites with bridging hydroxyls (Schroder et al. 1992)", "source": REF + ", Table 1",
       "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 7.5,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_schroder1992.py) from the paper's Table 1 and eq. 1-2",
                 "O-H Morse as a pair term within 1.4 A, the O-H Coulomb term kept (no bond, no exclusion)",
                 "O-T-O bends k (theta - theta0)^2 as eq. 2 writes it, O within 2.2 A of the T (the distance is CAPS's)",
                 "short-range cut-off 7.5 A as the paper", "shells carry 10 % of the oxygen's mass (adiabatic core-shell)"],
       "version": "1992 (Table 1)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": {"O_szc": "O_szs"},
           "description": "Si, Al, framework and hydroxyl O, hydroxyl H (bench/ff/make_schroder1992.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(angles)} bend rules")
