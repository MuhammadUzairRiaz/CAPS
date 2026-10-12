#!/usr/bin/env python3
"""Alkali halides for CAPS from their paper: C. R. A. Catlow, K. M. Diller, M. J. Norgett, J. Phys. C: Solid State
Phys. 10, 1395 (1977), Table 2 ("potential parameters - set 1"; the rigid-ion set with fitted charges).

  V+-(r) = A+- exp(-r/rho+-)            (per substance)
  V++(r) = A++ exp(-r/rho++) - C/r^6    (A++, rho++ per cation; C per substance)
  V--(r) = A-- exp(-r/rho--) - C/r^6    (A--, rho-- per anion; the same C)
  ion charges +-Z per substance

Each substance is its own set: its two ions typed by the structure's composition (the cation with that anion and no
other halide or alkali). The paper fitted V+- to nearest neighbours and V++, V-- to second neighbours; CAPS applies
them within the force field's cut-off (10 A).

usage: make_catlow1977.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "C. R. A. Catlow, K. M. Diller, M. J. Norgett, J. Phys. C: Solid State Phys. 10, 1395 (1977)"
EV = 23.060548
# Table 2(a), (b): A (eV), rho (A)
CATION = {"Li": (1153.8, 0.1364), "Na": (7895.4, 0.1709), "K": (3796.9, 0.2603), "Rb": (6808.1, 0.2704)}
ANION = {"F": (1127.7, 0.2753), "Cl": (1227.2, 0.3214), "Br": (2948.2, 0.3164), "I": (5502.5, 0.3066)}
# Table 2(c): substance -> Z, A+- (eV), rho+- (A), C (eV A^6)
SUBST = {("Li", "F"): (0.981, 485.6, 0.2614, 12.71), ("Li", "Cl"): (0.995, 1380.6, 0.2786, 0.34), ("Li", "Br"): (0.992, 1662.8, 0.2882, 0.01),
         ("Li", "I"): (0.979, 8894.7, 0.2617, 0.0), ("Na", "F"): (0.981, 1594.2, 0.2555, 11.68), ("Na", "Cl"): (0.988, 2314.7, 0.2903, 29.06),
         ("Na", "Br"): (0.988, 2593.1, 0.3014, 42.08), ("Na", "I"): (0.988, 3705.5, 0.3099, 27.38), ("K", "F"): (0.969, 2426.8, 0.2770, 44.60),
         ("K", "Cl"): (0.994, 4117.9, 0.3048, 124.9), ("K", "Br"): (0.988, 2944.0, 0.3295, 183.8), ("K", "I"): (0.994, 5842.4, 0.3264, 285.9),
         ("Rb", "F"): (0.960, 2209.3, 0.2945, 69.83), ("Rb", "Cl"): (0.991, 3559.9, 0.3233, 201.2), ("Rb", "Br"): (0.998, 4169.9, 0.3321, 261.4),
         ("Rb", "I"): (0.999, 5282.0, 0.3431, 400.6)}
MASS = {"Li": 6.94, "Na": 22.98977, "K": 39.0983, "Rb": 85.4678, "F": 18.998403, "Cl": 35.453, "Br": 79.904, "I": 126.90447}
Z = {"Li": 3, "Na": 11, "K": 19, "Rb": 37, "F": 9, "Cl": 17, "Br": 35, "I": 53}

types, pairs, rules = [], [], []
for (c, a), (q, A, rho, C) in SUBST.items():
    s = c + a
    tc, ta = f"{c}_{s}", f"{a}_{s}"
    types += [{"name": tc, "element": c, "mass": MASS[c], "charge": q, "description": f"{c}+ in {s} (Table 2c: Z {q})"},
              {"name": ta, "element": a, "mass": MASS[a], "charge": -q, "description": f"{a}- in {s} (Table 2c: Z {q})"}]
    pairs += [{"name": tc, "match": [tc], "params": [0.0, 0.0]}, {"name": ta, "match": [ta], "params": [0.0, 0.0]},
              {"name": f"{tc}-{tc}", "match": [tc, tc], "style": "buck", "params": [CATION[c][0] * EV, CATION[c][1], C * EV],
               "comment": f"{c}+-{c}+ A {CATION[c][0]} eV, rho {CATION[c][1]} A (Table 2a); C {C} eV A^6 ({s}, Table 2c)"},
              {"name": f"{ta}-{ta}", "match": [ta, ta], "style": "buck", "params": [ANION[a][0] * EV, ANION[a][1], C * EV],
               "comment": f"{a}--{a}- A {ANION[a][0]} eV, rho {ANION[a][1]} A (Table 2b); C {C} eV A^6 ({s}, Table 2c)"},
              {"name": f"{tc}-{ta}", "match": [tc, ta], "style": "buck", "params": [A * EV, rho, 0.0], "comment": f"{s}: A {A} eV, rho {rho} A (Table 2c)"}]
    others = [e for e in list(CATION) + list(ANION) if e not in (c, a)]
    rules += [{"type": tc, "smarts": f"[#{Z[c]}]", "priority": 0, "requires": [a], "excludes": others, "description": f"{c}+ in {s}"},
              {"type": ta, "smarts": f"[#{Z[a]}]", "priority": 0, "requires": [c], "excludes": others, "description": f"{a}- in {s}"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Alkali halides (Catlow, Diller & Norgett 1977)", "source": REF + ", Table 2 (set 1)",
       "references": [REF], "units": "real", "styles": {"pair": "buck/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_catlow1977.py) from the paper's Table 2 (set 1: rigid ions with fitted charges)",
                 "the paper fitted V+- to nearest and V++, V-- to second neighbours; CAPS applies them within 10 A"],
       "version": "1977 (set 1)", "typing": "../typing/catlow1977-halides.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", "catlow1977-halides.json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "rules": rules,
           "description": "each ion typed by the halide the structure holds (bench/ff/make_catlow1977.py)"},
          open(os.path.join(DATA, "typing", "catlow1977-halides.typing.json"), "w"), indent=1)
print(f"catlow1977-halides: {len(types)} types, {len(pairs)} pair terms, {len(rules)} rules")
