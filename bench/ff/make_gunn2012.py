#!/usr/bin/env python3
"""Gd pyrochlores for CAPS from their paper: D. S. D. Gunn, N. L. Allan, H. Foxhall, J. H. Harding, J. A. Purton, W. Smith,
M. J. Stein, I. T. Todorov, K. P. Travis, J. Mater. Chem. 22, 4675 (2012), Table 1.

  rigid ions with formal charges (Gd +3, Ti/Zr +4, O -2); Buckingham A exp(-r/rho) - C/r^6 for O-O, Ti-O, Gd-O, Zr-O
  (cation-cation Coulomb only); short-range cut-off 12 A as the paper

usage: make_gunn2012.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = ("D. S. D. Gunn, N. L. Allan, H. Foxhall, J. H. Harding, J. A. Purton, W. Smith, M. J. Stein, I. T. Todorov, K. P. Travis, "
       "J. Mater. Chem. 22, 4675 (2012)")
EV = 23.060548
FID = "gunn2012-pyrochlore"
ION = {"Gd": (3.0, 157.25, 64), "Ti": (4.0, 47.867, 22), "Zr": (4.0, 91.224, 40), "O": (-2.0, 15.999, 8)}
T1 = [("O", 17428.9200, 0.1490, 27.89), ("Ti", 648.6791, 0.4065, 0.0), ("Gd", 1962.7412, 0.3250, 0.0), ("Zr", 806.0606, 0.4011, 0.0)]
types = [{"name": f"{el}_g", "element": el, "mass": m, "charge": q, "description": f"{el} ({q:+g}, formal; Table 1)"} for el, (q, m, z) in ION.items()]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
pairs += [{"name": f"{el}_g-O_g", "match": [f"{el}_g", "O_g"], "style": "buck", "params": [A * EV, rho, C * EV],
           "comment": f"{el}-O A {A} eV, rho {rho} A, C {C} eV A^6 (Table 1)"} for el, A, rho, C in T1]
rules = [{"type": f"{el}_g", "smarts": f"[#{z}]", "priority": 0, "requires": ["O"] if el != "O" else ["Gd"], "description": f"{el} in a Gd pyrochlore"}
         for el, (q, m, z) in ION.items()]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Gd pyrochlores (Gunn et al. 2012)", "source": REF + ", Table 1", "references": [REF],
       "units": "real", "styles": {"pair": "buck/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 12,
       "atom_types": types, "pairs": pairs, "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_gunn2012.py) from the paper's Table 1", "short-range cut-off 12 A as the paper"],
       "version": "2012 (Table 1)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "rules": rules,
           "description": "Gd, Ti, Zr and O of the Gd pyrochlores (bench/ff/make_gunn2012.py)"},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types")
