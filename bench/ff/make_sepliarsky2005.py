#!/usr/bin/env python3
"""Ferroelectric perovskites for CAPS from M. Sepliarsky, A. Asthagiri, S. R. Phillpot, M. G. Stachiotti, R. L. Migoni,
Curr. Opin. Solid State Mater. Sci. 9, 107 (2005), Table 1 (eV, A, e).

Every ion a core and a shell; Coulomb between all cores and shells except an ion's own core and shell; short-range terms
between shells (section 2: "Short-ranged interactions are normally restricted to act between shells only").
  core-shell  V(w) = c2 w^2/2 + c4 w^4/24 (section 2): class II bond, r0 0, K2 = c2/2, K3 0, K4 = c4/24
  K(Nb,Ta)O3  Buckingham A exp(-r/rho) - C/r^6 for K-O, Nb-O, Ta-O, O-O (section 2); KNbO3, KTaO3 and their solid
              solutions share K, O and the O-O term (section 3)
  (Ba,Sr)TiO3 Buckingham for Ba-O, Sr-O, Ti-O, O-O; BaTiO3, SrTiO3 and their solid solutions share Ti, O
  PbTiO3      O-O Buckingham, every other short-range pair Rydberg (A + C r) exp(-r/rho): Pb-Ti, Pb-O, Ti-O (section 2)
Pairs Table 1 does not list (A-A, B-B, A-B in K(Nb,Ta)O3 and (Ba,Sr)TiO3) are Coulomb only. Each composition family is a
set typed by the elements present: K with Nb / Ta, Ba / Sr with Ti (no Pb), Pb with Ti.

PbTiO3 is transcribed but not written by default. With the form section 2 states, Table 1's Pb-Ti row
(0.09 + 12.57 r) exp(-r/2.4201) eV is +10.5 eV per Pb-Ti pair at 3.4 A and still +2.8 eV at 9 A: cubic PbTiO3 comes out at
+141 eV per formula unit (a = 3.97 A), i.e. unbound. The Pb-O and Ti-O rows read the same way are physical (without
Pb-Ti, -146.9 eV per formula unit with its minimum near a = 3.85-3.9 A, the LDA lattice constant the set was fitted to),
while reading all three rows as Buckingham terms leaves no minimum between 3.8 and 4.15 A. The Pb-Ti row (or its form)
is misprinted in the review; the parameterisation is in M. Sepliarsky, R. E. Cohen, AIP Conf. Proc. 626, 36 (2002) and
M. Sepliarsky, Z. Wu, A. Asthagiri, R. E. Cohen, Ferroelectrics 301, 55 (2004) (refs. 27, 28), not here.
SEPLIARSKY_SETS=knt,bst,pt writes it anyway (for checking the Rydberg form; not for use).

usage: [SEPLIARSKY_SETS=knt,bst] make_sepliarsky2005.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "M. Sepliarsky, A. Asthagiri, S. R. Phillpot, M. G. Stachiotti, R. L. Migoni, Curr. Opin. Solid State Mater. Sci. 9, 107 (2005)"
EV = 23.060548   # kcal/mol per eV
FID = "sepliarsky2005-perovskites"
MASS = {"K": 39.098, "Nb": 92.906, "Ta": 180.948, "O": 15.999, "Ba": 137.327, "Sr": 87.62, "Ti": 47.867, "Pb": 207.2}
NUM = {"K": 19, "Nb": 41, "Ta": 73, "O": 8, "Ba": 56, "Sr": 38, "Ti": 22, "Pb": 82}
# Table 1, per set: ions (element, core charge, shell charge, k2, k4) and short-range rows (a, b, form, A, rho, C)
SETS = {
    "knt": ("K(Nb,Ta)O3",
            [("K", 1.2394, -0.4189, 225.0, 0.0), ("Nb", -2.9878, 7.8265, 250.67, 400.0), ("Ta", -2.9878, 7.8265, 295.0, 400.0),
             ("O", 1.1236, -3.010, 75.0, 1500.0)],
            [("K", "O", "buck", 124872.44, 0.1950, 0.0), ("Nb", "O", "buck", 1036.63, 0.3900, 0.0), ("Ta", "O", "buck", 1046.63, 0.3900, 0.0),
             ("O", "O", "buck", 3600.26, 0.2834, 200.0)]),
    "bst": ("(Ba,Sr)TiO3",
            [("Ba", 5.62, -3.76, 251.8, 0.0), ("Sr", 5.62, -3.76, 623.3, 0.0), ("Ti", 4.76, -1.58, 322.0, 500.0), ("O", 0.91, -2.59, 31.0, 4000.0)],
            [("Ba", "O", "buck", 1061.30, 0.3740, 0.0), ("Sr", "O", "buck", 167.32, 0.4949, 0.0), ("Ti", "O", "buck", 3769.93, 0.2589, 0.0),
             ("O", "O", "buck", 4740.00, 0.2686, 160.0)]),
    "pt": ("PbTiO3",
           [("Pb", 4.958, -2.785, 119.5, 17968.0), ("Ti", 8.820, -5.158, 1428.0, 0.0), ("O", 0.563, -2.508, 23.3, 4514.0)],
           [("Pb", "Ti", "rydberg", 0.09, 2.4201, 12.57), ("Pb", "O", "rydberg", 6766.27, 0.2738, 127.78), ("Ti", "O", "rydberg", 1130.01, 0.3597, -160.84),
            ("O", "O", "buck", 3634.86, 0.3144, 331.60)]),
}
# which elements a set's ions need beside O, and which they must not see
SET_REQ = {"knt": (["K"], ["Ba", "Sr", "Ti", "Pb"]), "bst": (["Ti"], ["K", "Nb", "Ta", "Pb"]), "pt": (["Pb", "Ti"], ["K", "Nb", "Ta", "Ba", "Sr"])}
ONLY = [x for x in os.environ.get("SEPLIARSKY_SETS", "knt,bst").split(",") if x]
types, pairs, bonds, rules, shells = [], [], [], [], {}
for tag in ONLY:
    label, ions, sr = SETS[tag]
    need, avoid = SET_REQ[tag]
    for el, qc, qs, k2, k4 in ions:
        c, s = f"{el}_{tag}c", f"{el}_{tag}s"
        types += [{"name": c, "element": el, "mass": MASS[el] * 0.9, "charge": qc, "description": f"{el} core, {label} (Table 1)"},
                  {"name": s, "element": el, "mass": MASS[el] * 0.1, "charge": qs, "description": f"{el} shell, {label} (Table 1)"}]
        shells[c] = s
        bonds.append({"name": f"{c}-{s}", "match": [c, s], "style": "class2", "params": [0.0, k2 / 2 * EV, 0.0, k4 / 24 * EV],
                      "comment": f"c2 {k2} eV/A^2, c4 {k4} eV/A^4: c2 w^2/2 + c4 w^4/24 ({label}, Table 1)"})
        req = sorted(set(need + ["O"]) - {el}) if el != "O" else need
        rule = {"type": c, "smarts": f"[#{NUM[el]}]", "priority": 0, "requires": req, "excludes": avoid, "description": f"{el} in {label}"}
        rules.append(rule)
    for a, b, form, A, rho, C in sr:
        sa, sb = f"{a}_{tag}s", f"{b}_{tag}s"
        if form == "buck":
            pairs.append({"name": f"{sa}-{sb}", "match": [sa, sb], "style": "buck", "params": [A * EV, rho, C * EV],
                          "comment": f"{a}-{b} A {A} eV, rho {rho} A, C {C} eV A^6: A exp(-r/rho) - C/r^6 ({label}, Table 1)"})
        else:
            pairs.append({"name": f"{sa}-{sb}", "match": [sa, sb], "style": "rydberg", "params": [A * EV, rho, C * EV],
                          "comment": f"{a}-{b} A {A} eV, rho {rho} A, C {C} eV/A: (A + C r) exp(-r/rho) ({label}, Table 1)"})
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types] + pairs
rules += [{"type": s, "smarts": f"[#{NUM[s.split('_')[0]]}]", "priority": 99, "atom_name": s, "description": f"shell of {c}"} for c, s in shells.items()]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Ferroelectric perovskites (Sepliarsky et al. 2005)", "source": REF + ", Table 1",
       "references": [REF], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_sepliarsky2005.py) from the paper's Table 1 (signs as the rendered table prints them)",
                 "every ion a core and a shell; short-range terms between shells; pairs Table 1 does not list are Coulomb only",
                 "core-shell V(w) = c2 w^2/2 + c4 w^4/24 (section 2) as a class II bond: K2 = c2/2, K4 = c4/24",
                 "K(Nb,Ta)O3 and (Ba,Sr)TiO3: Buckingham pairs (section 2); each set typed by composition (K with Nb / Ta; Ba / Sr with Ti)",
                 "PbTiO3 is left out: Table 1's Pb-Ti row read as the stated Rydberg form (A + C r) exp(-r/rho) makes cubic PbTiO3 unbound (+141 eV per formula unit); the original parameter papers (Sepliarsky & Cohen 2002; Sepliarsky et al. 2004) are needed",
                 "the review states no cut-off; 10 A is CAPS's choice",
                 "Ta-O A 1046.63 eV and Nb-O 1036.63 eV as Table 1 prints them",
                 "shells carry 10 % of the ion's mass (adiabatic core-shell)"],
       "version": "2005 (Table 1)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": shells,
           "description": "ABO3 perovskites by composition (bench/ff/make_sepliarsky2005.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(bonds)} springs")
