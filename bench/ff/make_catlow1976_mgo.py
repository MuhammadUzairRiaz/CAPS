#!/usr/bin/env python3
"""MgO for CAPS from its paper: C. R. A. Catlow, I. D. Faux, M. J. Norgett, J. Phys. C: Solid State Phys. 9, 419 (1976),
Table 1(b) (short-range potentials) and Table 1(c) (shell parameters).

  Mg-O   A+- exp(-r/rho+-) between the Mg ion and the O shell
  O-O    A-- exp(-r/rho--) - C--/r^6 between O shells
  Mg-Mg  Coulomb only ("We have neglected the cation-cation interaction", section 2.1)
  charges Z = +-2 (Table 1(b)); only the anion is polarisable (section 2.1): Mg a rigid ion +2, O core -2 - Y-, O shell
          Y-, spring 1/2 K- r^2 (class II K2 = K-/2), the harmonic shell-model spring as in Catlow's other papers
Model 1 (fitted to r0, C11 - C12 and C44) is the default; model 2 (fitted to r0, C11 and C12, many-body effects in an
effective pair potential) is the overlay data/typing/sets/catlow1976-mgo/model2.typing.json. The breathing shell model
(BSM, eq. 2.1.1) has no CAPS form and is not written; neither are Sangster's and Boswarva's columns, which the paper
quotes from other work.

usage: make_catlow1976_mgo.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "C. R. A. Catlow, I. D. Faux, M. J. Norgett, J. Phys. C: Solid State Phys. 9, 419 (1976)"
EV = 23.060548   # kcal/mol per eV
FID = "catlow1976-mgo"
# Table 1(b): A+- (eV), rho+- (A), A-- (eV), rho-- (A), C-- (eV A^6); Table 1(c): Y- (|e|), K- (eV A^-2)
T1 = {1: (1152.0, 0.3065, 22760.0, 0.149, 28.96, -2.62, 38.03),
      2: (553.9, 0.3664, 22760.0, 0.149, 255.2, -2.06, 18.30)}
MASS = {"Mg": 24.305, "O": 15.999}
types, pairs, bonds, rules, shells, overlay = [], [], [], [], {}, []
for n, (Apm, rpm, Amm, rmm, Cmm, Y, K) in T1.items():
    mg, oc, os_ = f"Mg_c76m{n}", f"O_c76m{n}c", f"O_c76m{n}s"
    ts = [{"name": mg, "element": "Mg", "mass": MASS["Mg"], "charge": 2.0, "description": f"Mg2+, rigid ion, model {n} (Table 1(b): Z = 2)"},
          {"name": oc, "element": "O", "mass": MASS["O"] * 0.9, "charge": round(-2.0 - Y, 6), "description": f"O2- core, model {n} (-2 - Y-, Table 1(c) Y- {Y})"},
          {"name": os_, "element": "O", "mass": MASS["O"] * 0.1, "charge": Y, "description": f"O2- shell, model {n} (Table 1(c))"}]
    types += ts
    shells[oc] = os_
    pairs += [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in ts]
    pairs += [{"name": f"{mg}-{os_}", "match": [mg, os_], "style": "buck", "params": [Apm * EV, rpm, 0.0],
               "comment": f"Mg-O A+- {Apm} eV, rho+- {rpm} A (Table 1(b), model {n})"},
              {"name": f"{os_}-{os_}", "match": [os_, os_], "style": "buck", "params": [Amm * EV, rmm, Cmm * EV],
               "comment": f"O-O A-- {Amm} eV, rho-- {rmm} A, C-- {Cmm} eV A^6 (Table 1(b), model {n})"}]
    bonds.append({"name": f"{oc}-{os_}", "match": [oc, os_], "style": "class2", "params": [0.0, K / 2 * EV, 0.0, 0.0],
                  "comment": f"O core-shell spring K- {K} eV/A^2, 1/2 K- r^2 (Table 1(c), model {n})"})
    rs = [{"type": mg, "smarts": "[#12]", "priority": 0, "requires": ["O"], "description": f"Mg in MgO (model {n})"},
          {"type": oc, "smarts": "[#8]", "priority": 0, "requires": ["Mg"], "description": f"O in MgO (model {n})"}]
    if n == 1:
        rules += rs
    else:
        overlay = rs
rules += [{"type": sh, "smarts": "[#8]", "priority": 99, "atom_name": sh, "description": f"shell of {c}"} for c, sh in shells.items()]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "MgO (Catlow, Faux & Norgett 1976)", "source": REF + ", Table 1(b), 1(c)",
       "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_catlow1976_mgo.py) from the paper's Table 1(b) and 1(c)",
                 "model 1 (fitted to r0, C11 - C12, C44) by default; model 2 (fitted to r0, C11, C12) is the overlay 'model2'",
                 "Mg rigid (+2), O a shell model: only the anion is polarisable (section 2.1); Mg-Mg Coulomb only",
                 "the spring is read as 1/2 K- r^2 (the paper gives K- in eV A^-2 without the energy expression)",
                 "the breathing shell model (BSM) is not written: CAPS has no breathing-shell term",
                 "the paper states no cut-off; 10 A is CAPS's choice",
                 "shells carry 10 % of the ion's mass (adiabatic core-shell)"],
       "version": "1976 (Table 1, models 1 and 2)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": shells,
           "description": "Mg and O of magnesium oxide (bench/ff/make_catlow1976_mgo.py); model 2 is an overlay", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
os.makedirs(os.path.join(DATA, "typing", "sets", FID), exist_ok=True)
json.dump({"format": "caps-typing", "version": 1, "forcefield": f"{FID}: model2", "unknown_types": "untyped",
           "description": "Overlay: Catlow, Faux & Norgett's model 2 (Table 1) in place of model 1", "rules": [dict(r, priority=50) for r in overlay]},
          open(os.path.join(DATA, "typing", "sets", FID, "model2.typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(bonds)} springs")
