#!/usr/bin/env python3
"""Rutile TiO2 for CAPS from C. R. A. Catlow, C. M. Freeman, R. L. Royle, Physica B+C 131, 1 (1985), Table IVa.

  Ti-O   A12 exp(-r/rho12) between the Ti shell and the O shell (Table IVa, footnote a)
  O-O    A22 exp(-r/rho22) - C22 r^-6 between O shells (footnote b)
  Ti-Ti  Coulomb only (Table IVa lists none)
  shells cation charge Z = 4.0, cation shell Y1 = 2.89, anion shell Y2 = -2.53 (anion charge -Z/2 = -2): Ti core 1.11,
         O core 0.53; springs K1 37.3 (Ti) and 86.4 (O) eV A^-2, written 1/2 K w^2
  quartic footnote c: "core-shell self-energy includes quartic as well as harmonic terms. For both Ti and O quartic
         constant, K2 = 1.1 x 10^5 eV/A^5". The paper prints neither the quartic energy expression nor a consistent unit
         (a w^4 constant is eV A^-4). CAPS reads it as K w^4 / 24 (4! in the denominator: the anharmonic core-shell
         form 1/2 k2 w^2 + 1/24 k4 w^4 of Catlow's group's lattice codes, which Sepliarsky et al. 2005 also print), with
         K in eV A^-4: class II bond K2 = K1/2, K4 = 1.1e5/24. Catlow & James, Proc. R. Soc. Lond. A 384, 157 (1982)
         does not settle it: its rutile potential (Table 1) is an earlier, harmonic-only shell model (K1 65974.0,
         K2 18.41 eV A^-2) with no quartic term.
The paper's alpha-quartz model (Table Ia) is Sanders, Leslie & Catlow's set, already the library sanders1984-sio2 (every
number equal), so it is not written again here.

usage: make_catlow1985_rutile.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "C. R. A. Catlow, C. M. Freeman, R. L. Royle, Physica B+C 131, 1 (1985)"
EV = 23.060548   # kcal/mol per eV
FID = "catlow1985-rutile"
# Table IVa
Z, A12, R12, A22, R22, C22, Y1, Y2, K1, K2 = 4.0, 754.2, 0.3879, 22764.3, 0.149, 27.88, 2.89, -2.53, 37.3, 86.4
K4 = 1.1e5   # footnote c, both ions
MASS = {"Ti": 47.867, "O": 15.999}
types = [{"name": "Ti_cfrc", "element": "Ti", "mass": MASS["Ti"] * 0.9, "charge": round(Z - Y1, 6), "description": f"Ti4+ core (Z {Z} - Y1 {Y1}; Table IVa)"},
         {"name": "Ti_cfrs", "element": "Ti", "mass": MASS["Ti"] * 0.1, "charge": Y1, "description": "Ti4+ shell (Table IVa)"},
         {"name": "O_cfrc", "element": "O", "mass": MASS["O"] * 0.9, "charge": round(-Z / 2 - Y2, 6), "description": f"O2- core (-Z/2 - Y2 {Y2}; Table IVa)"},
         {"name": "O_cfrs", "element": "O", "mass": MASS["O"] * 0.1, "charge": Y2, "description": "O2- shell (Table IVa)"}]
shells = {"Ti_cfrc": "Ti_cfrs", "O_cfrc": "O_cfrs"}
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
pairs += [{"name": "Ti_cfrs-O_cfrs", "match": ["Ti_cfrs", "O_cfrs"], "style": "buck", "params": [A12 * EV, R12, 0.0],
           "comment": f"Ti-O A12 {A12} eV, rho12 {R12} A (Table IVa, footnote a)"},
          {"name": "O_cfrs-O_cfrs", "match": ["O_cfrs", "O_cfrs"], "style": "buck", "params": [A22 * EV, R22, C22 * EV],
           "comment": f"O-O A22 {A22} eV, rho22 {R22} A, C22 {C22} eV A^6 (Table IVa, footnote b)"}]
bonds = [{"name": f"{c}-{s}", "match": [c, s], "style": "class2", "params": [0.0, K / 2 * EV, 0.0, K4 / 24 * EV],
          "comment": f"spring {K} eV/A^2 (1/2 K w^2) and quartic {K4:g} (read as K w^4/24) (Table IVa, footnote c)"}
         for (c, s), K in zip(shells.items(), (K1, K2))]
rules = [{"type": "Ti_cfrc", "smarts": "[#22]", "priority": 0, "requires": ["O"], "description": "Ti in rutile"},
         {"type": "O_cfrc", "smarts": "[#8]", "priority": 0, "requires": ["Ti"], "description": "O in rutile"},
         {"type": "Ti_cfrs", "smarts": "[#22]", "priority": 99, "atom_name": "Ti_cfrs", "description": "shell of Ti_cfrc"},
         {"type": "O_cfrs", "smarts": "[#8]", "priority": 99, "atom_name": "O_cfrs", "description": "shell of O_cfrc"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Rutile (Catlow, Freeman & Royle 1985)", "source": REF + ", Table IVa",
       "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_catlow1985_rutile.py) from the paper's Table IVa",
                 "Ti and O shells; Ti-Ti Coulomb only",
                 "the quartic core-shell constant 1.1e5 (footnote c, printed in eV/A^5) is read as K w^4/24 with K in eV/A^4; the paper prints no expression for it; Catlow & James 1982, the earlier rutile model, is harmonic only and does not settle it",
                 "the paper states no cut-off; 10 A is CAPS's choice",
                 "the paper's alpha-quartz set (Table Ia) is sanders1984-sio2",
                 "shells carry 10 % of the ion's mass (adiabatic core-shell)"],
       "version": "1985 (Table IVa)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": shells,
           "description": "Ti and O of rutile (bench/ff/make_catlow1985_rutile.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(bonds)} springs")
