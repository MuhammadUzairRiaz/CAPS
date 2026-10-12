#!/usr/bin/env python3
"""UO2 for CAPS from its paper: C. R. A. Catlow, Proc. R. Soc. Lond. A 353, 533 (1977), Table 3 (both potentials).

  U-O   A12 exp(-r/rho12) between the U shell and the O shell
  O-O   A22 exp(-r/rho22) - C22/r^6 between O shells
  U-U   Coulomb only
  shells  U: core Z - Y+, shell Y+, spring K+; O: core -Z/2 - Y-, shell Y-, spring K-; springs 1/2 K r^2 (class II K2 = K/2)
Potential 1 (Born-Mayer terms from Hartree-Fock, C fitted; for separations near the perfect lattice's) is the default;
potential 2 (C from Hartree-Fock, the ionic charge fitted; for short O-O separations) is the overlay
data/typing/sets/catlow1977-uo2/potential2.typing.json.

usage: make_catlow1977_uo2.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "C. R. A. Catlow, Proc. R. Soc. Lond. A 353, 533 (1977)"
EV = 23.060548
FID = "catlow1977-uo2"
# Table 3: Z, A12 (eV), rho12 (A), A22 (eV), rho22 (A), C22 (eV A^6), Y+, Y-, K+ (eV/A^2), K- (eV/A^2)
T3 = {1: (4.0, 1217.8, 0.3871, 22764.3, 0.1490, 112.2, 6.54, -4.40, 103.38, 292.98),
      2: (3.84, 1297.44, 0.3747, 22764.3, 0.1490, 20.37, 7.94, -3.06, 210.02, 80.21)}
MASS = {"U": 238.029, "O": 15.999}
types, pairs, bonds, rules, shells, overlay = [], [], [], [], {}, []
for n, (Zc, A12, r12, A22, r22, C22, Yp, Ym, Kp, Km) in T3.items():
    uc, us, oc, os_ = f"U{n}c", f"U{n}s", f"O{n}c", f"O{n}s"
    ts = [{"name": uc, "element": "U", "mass": MASS["U"] * 0.9, "charge": round(Zc - Yp, 6), "description": f"U core, potential {n} (Table 3: Z {Zc}, Y+ {Yp})"},
          {"name": us, "element": "U", "mass": MASS["U"] * 0.1, "charge": Yp, "description": f"U shell, potential {n}"},
          {"name": oc, "element": "O", "mass": MASS["O"] * 0.9, "charge": round(-Zc / 2 - Ym, 6), "description": f"O core, potential {n} (Table 3: Y- {Ym})"},
          {"name": os_, "element": "O", "mass": MASS["O"] * 0.1, "charge": Ym, "description": f"O shell, potential {n}"}]
    types += ts
    shells.update({uc: us, oc: os_})
    pairs += [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in ts]
    pairs += [{"name": f"{us}-{os_}", "match": [us, os_], "style": "buck", "params": [A12 * EV, r12, 0.0],
               "comment": f"U-O A12 {A12} eV, rho12 {r12} A (Table 3, potential {n})"},
              {"name": f"{os_}-{os_}", "match": [os_, os_], "style": "buck", "params": [A22 * EV, r22, C22 * EV],
               "comment": f"O-O A22 {A22} eV, rho22 {r22} A, C22 {C22} eV A^6 (Table 3, potential {n})"}]
    bonds += [{"name": f"{uc}-{us}", "match": [uc, us], "style": "class2", "params": [0.0, Kp / 2 * EV, 0.0, 0.0], "comment": f"K+ {Kp} eV/A^2 (Table 3, potential {n})"},
              {"name": f"{oc}-{os_}", "match": [oc, os_], "style": "class2", "params": [0.0, Km / 2 * EV, 0.0, 0.0], "comment": f"K- {Km} eV/A^2 (Table 3, potential {n})"}]
    rs = [{"type": uc, "smarts": "[#92]", "priority": 0, "requires": ["O"], "description": f"U in UO2 (potential {n})"},
          {"type": oc, "smarts": "[#8]", "priority": 0, "requires": ["U"], "description": f"O in UO2 (potential {n})"}]
    if n == 1:
        rules += rs
    else:
        overlay = rs
rules += [{"type": sh, "smarts": "[#92]" if sh.startswith("U") else "[#8]", "priority": 99, "atom_name": sh, "description": f"shell of {c}"}
          for c, sh in shells.items()]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "UO2 (Catlow 1977)", "source": REF + ", Table 3", "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_catlow1977_uo2.py) from the paper's Table 3",
                 "potential 1 by default; potential 2 (fitted charge 3.84) is the overlay 'potential2'",
                 "shells carry 10 % of the ion's mass (adiabatic core-shell)"],
       "version": "1977 (Table 3)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": shells,
           "description": "U and O of uranium dioxide (bench/ff/make_catlow1977_uo2.py); potential 2 is an overlay", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
os.makedirs(os.path.join(DATA, "typing", "sets", FID), exist_ok=True)
json.dump({"format": "caps-typing", "version": 1, "forcefield": f"{FID}: potential2", "unknown_types": "untyped",
           "description": "Overlay: Catlow's potential 2 (Table 3) in place of potential 1", "rules": [dict(r, priority=50) for r in overlay]},
          open(os.path.join(DATA, "typing", "sets", FID, "potential2.typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(bonds)} springs")
