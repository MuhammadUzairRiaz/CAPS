#!/usr/bin/env python3
"""ZrO2 for CAPS from its paper: M. S. Khan, M. S. Islam, D. R. Bates, J. Mater. Chem. 8, 2299 (1998), Table 1.

  V(r) = -Z_a Z_b e^2 / (eps r) + A exp(-r/rho) - C/r^6   (eq. 1), formal charges Zr4+ and O2-
  Zr-O   985.869 exp(-r/0.3760) eV between the Zr shell and the O shell (Table 1(i))
  O-O    22764.300 exp(-r/0.1490) - 27.890/r^6 eV between O shells (Table 1(i))
  Zr-Zr  Coulomb only (Table 1 lists none)
  shells Table 1(ii): Zr Y 1.350, k 169.617; O Y -2.077, k 27.290 eV A^-2 ("Y and k refer to the shell charge and harmonic
         force constant"): core = Z - Y, spring 1/2 k r^2 (class II K2 = k/2)
Cut-off 10.254 A (Table 1 note). The paper's MD runs used the rigid-ion form of the same set (section "Simulation
methods": the shell model was not used in MD); the shell model of Table 1(ii) is what CAPS writes.

usage: make_khan1998_zro2.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "M. S. Khan, M. S. Islam, D. R. Bates, J. Mater. Chem. 8, 2299 (1998)"
EV = 23.060548   # kcal/mol per eV
FID = "khan1998-zro2"
Z = {"Zr": 4.0, "O": -2.0}
SHELL = {"Zr": (1.350, 169.617), "O": (-2.077, 27.290)}            # Table 1(ii): Y (e), k (eV A^-2)
SR = [("Zr", "O", 985.869, 0.3760, 0.000), ("O", "O", 22764.300, 0.1490, 27.890)]   # Table 1(i): A (eV), rho (A), C (eV A^6)
MASS = {"Zr": 91.224, "O": 15.999}
NUM = {"Zr": 40, "O": 8}
types, pairs, bonds, shells = [], [], [], {}
for el in ("Zr", "O"):
    Y, k = SHELL[el]
    c, s = f"{el}_k98c", f"{el}_k98s"
    types += [{"name": c, "element": el, "mass": MASS[el] * 0.9, "charge": round(Z[el] - Y, 6), "description": f"{el} core (Z {Z[el]:g} - Y {Y}; Table 1(ii))"},
              {"name": s, "element": el, "mass": MASS[el] * 0.1, "charge": Y, "description": f"{el} shell (Table 1(ii))"}]
    shells[c] = s
    bonds.append({"name": f"{c}-{s}", "match": [c, s], "style": "class2", "params": [0.0, k / 2 * EV, 0.0, 0.0],
                  "comment": f"{el} core-shell spring k {k} eV/A^2, 1/2 k r^2 (Table 1(ii))"})
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
for a, b, A, rho, C in SR:
    pairs.append({"name": f"{a}_k98s-{b}_k98s", "match": [f"{a}_k98s", f"{b}_k98s"], "style": "buck", "params": [A * EV, rho, C * EV],
                  "comment": f"{a}-{b} A {A} eV, rho {rho} A, C {C} eV A^6 (Table 1(i))"})
rules = [{"type": "Zr_k98c", "smarts": "[#40]", "priority": 0, "requires": ["O"], "description": "Zr in zirconia"},
         {"type": "O_k98c", "smarts": "[#8]", "priority": 0, "requires": ["Zr"], "description": "O in zirconia"}]
rules += [{"type": s, "smarts": f"[#{NUM[s.split('_')[0]]}]", "priority": 99, "atom_name": s, "description": f"shell of {c}"} for c, s in shells.items()]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "ZrO2 (Khan, Islam & Bates 1998)", "source": REF + ", Table 1",
       "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10.254,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_khan1998_zro2.py) from the paper's Table 1",
                 "formal charges (eq. 1); Zr and O shells (Table 1(ii)), the springs 1/2 k r^2 (k the harmonic force constant)",
                 "Zr-Zr is Coulomb only: Table 1 lists no Zr-Zr term",
                 "cut-off 10.254 A (Table 1 note)",
                 "the paper's MD used this set without shells; the shell model is the one Table 1 gives",
                 "shells carry 10 % of the ion's mass (adiabatic core-shell)"],
       "version": "1998 (Table 1)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": shells,
           "description": "Zr and O of zirconia (bench/ff/make_khan1998_zro2.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(bonds)} springs")
