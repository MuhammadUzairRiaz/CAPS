#!/usr/bin/env python3
"""Calcium carbonate for CAPS from its paper: S. Xiao, S. A. Edwards, F. Grater, J. Phys. Chem. C 115, 20067 (2011), Table 1
and the energy expression on p. 20069.

  charges   Ca +1.668, C +0.999, O -0.889
  pairs     C12/r^12 - C6/r^6 for the pairs Table 1 lists (in sigma/epsilon: sigma = (C12/C6)^(1/6), epsilon = C6^2 / 4 C12);
            every other pair Coulomb only. Ca-O has C6 = 0: written with C6 = 1e-12 kJ mol^-1 nm^6 (the sigma/epsilon form
            needs C6 > 0; less than 1e-8 kJ/mol at contact)
  bend      1/2 K_theta (theta - 120)^2 on O-C-O, K_theta 1852 kJ/mol/rad^2 (the paper prints (theta - theta0) without the
            square, a misprint of the harmonic bend)
  plane     K_psi (1 - cos 2 psi), K_psi 28.9 kJ/mol, psi the angle between the O1-C-O2 and O2-C-O3 planes: one cvff improper
            O1-C-O2-O3 per carbonate (the O sorted by index)
  The carbonate has no C-O bond term and no exclusions: C and O interact by Coulomb and their Lennard-Jones pair; the bend
  and the plane term are found by distance (1.6 A from the C).
  water     TIP3P (W. L. Jorgensen et al., J. Chem. Phys. 79, 926 (1983)): O -0.834, H +0.417, O-O A 582000 kcal A^12,
            C 595 kcal A^6; the paper's aragonite-water pairs (Table 1); the paper runs it rigid, CAPS flexible with
            harmonic O-H (450 kcal/mol/A^2) and H-O-H (55 kcal/mol/rad^2) terms about the TIP3P geometry
  Table 1's atom-wise C6, C12 (for mineral-peptide pairs by the geometric rule, eq. 3) are kept in the notes.

usage: make_xiao2011.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "S. Xiao, S. A. Edwards, F. Grater, J. Phys. Chem. C 115, 20067 (2011)"
TIP3P = "W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, M. L. Klein, J. Chem. Phys. 79, 926 (1983)"
KJ = 4.184
FID = "xiao2011-caco3"
types = [{"name": "Ca_x", "element": "Ca", "mass": 40.078, "charge": 1.668, "description": "Ca in CaCO3 (Table 1)"},
         {"name": "Cm_x", "element": "C", "mass": 12.011, "charge": 0.999, "description": "carbonate C (Table 1)"},
         {"name": "Om_x", "element": "O", "mass": 15.999, "charge": -0.889, "description": "carbonate O (Table 1)"},
         {"name": "Ow_x", "element": "O", "mass": 15.9994, "charge": -0.834, "description": "TIP3P water O (Jorgensen 1983)"},
         {"name": "Hw_x", "element": "H", "mass": 1.008, "charge": 0.417, "description": "TIP3P water H (Jorgensen 1983)"}]
# Table 1: C12 (kJ/mol nm^12), C6 (kJ/mol nm^6)
T1 = [("Ca_x", "Om_x", 9.49e-7, 0.0), ("Cm_x", "Cm_x", 4.61e-6, 1.43e-2), ("Cm_x", "Om_x", 9.04e-10, 3.08e-4), ("Om_x", "Om_x", 5.94e-7, 5.21e-5),
      ("Ca_x", "Ow_x", 8.85e-7, 2.00e-3), ("Cm_x", "Ow_x", 5.91e-6, 3.39e-3), ("Om_x", "Ow_x", 2.07e-6, 2.25e-3), ("Om_x", "Hw_x", 7.98e-9, 9.65e-7)]
C6_ZERO = 1e-12


def lj(c12, c6):
    c6e = c6 if c6 > 0 else C6_ZERO
    return [c6e * c6e / (4 * c12) / KJ, (c12 / c6e) ** (1 / 6) * 10]


pairs = [{"name": "Ca_x", "match": ["Ca_x"], "params": [0.0, 0.0], "comment": "Ca-Ca: Coulomb only (not in Table 1)"},
         {"name": "Hw_x", "match": ["Hw_x"], "params": [0.0, 0.0], "comment": "TIP3P H: no Lennard-Jones"},
         {"name": "Ow_x", "match": ["Ow_x"], "params": [(595.0 ** 2) / (4 * 582000.0), (582000.0 / 595.0) ** (1 / 6)],
          "comment": "TIP3P O-O A 582000 kcal A^12, C 595 kcal A^6 (Jorgensen 1983)"}]
for a, b, c12, c6 in T1:
    pairs.append({"name": a if a == b else f"{a}-{b}", "match": [a] if a == b else [a, b], "params": lj(c12, c6),
                  "comment": f"C12 {c12:g} kJ/mol nm^12, C6 {c6:g} kJ/mol nm^6 (Table 1)" + ("; C6 0 written as 1e-12" if c6 == 0 else "")})
explicit = {(p["match"][0], p["match"][-1]) for p in pairs}
names = [t["name"] for t in types]
for i, a in enumerate(names):   # every pair Table 1 does not list: Coulomb only
    for b in names[i + 1:]:
        if (a, b) not in explicit and (b, a) not in explicit:
            pairs.append({"name": f"{a}-{b}", "match": [a, b], "params": [0.0, 0.0], "comment": "not in Table 1: Coulomb only"})
bonds = [{"name": "Ow_x-Hw_x", "match": ["Ow_x", "Hw_x"], "style": "harmonic", "params": [450.0, 0.9572], "comment": "TIP3P geometry; flexible in CAPS (rigid in the paper)"}]
angles = [{"name": "Om_x-Cm_x-Om_x", "match": ["Om_x", "Cm_x", "Om_x"], "style": "harmonic", "params": [1852.0 / 2 / KJ, 120.0], "contact": [1.6],
           "comment": "1/2 K_theta (theta - theta0)^2, K_theta 1852 kJ/mol/rad^2, theta0 120 (Table 1); O within 1.6 A of the C"},
          {"name": "Hw_x-Ow_x-Hw_x", "match": ["Hw_x", "Ow_x", "Hw_x"], "style": "harmonic", "params": [55.0, 104.52], "comment": "TIP3P geometry; flexible in CAPS"}]
impropers = [{"name": "Om_x-Cm_x-Om_x-Om_x", "match": ["Om_x", "Cm_x", "Om_x", "Om_x"], "style": "cvff", "params": [28.9 / KJ, -1, 2], "contact": [1.6],
              "comment": "K_psi (1 - cos 2 psi), K_psi 28.9 kJ/mol (Table 1); O1-C-O2-O3, the three O within 1.6 A of the C"}]
CO3 = "[#6;$([#6](~[#8])(~[#8])~[#8])]"
rules = [{"type": "Ca_x", "smarts": "[#20]", "priority": 0, "requires": ["C", "O"], "description": "Ca of calcium carbonate"},
         {"type": "Cm_x", "smarts": CO3, "priority": 0, "description": "carbonate C (its C-O bonds carry no term: kept for typing only)"},
         {"type": "Om_x", "smarts": f"[#8;$([#8]~{CO3})]", "priority": 0, "description": "carbonate O"},
         {"type": "Ow_x", "smarts": "[OX2H2]", "priority": 1, "description": "water O"},
         {"type": "Hw_x", "smarts": "[H][OX2H2]", "priority": 1, "description": "water H"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "CaCO3 (Xiao, Edwards & Grater 2011)", "source": REF + ", Table 1",
       "references": [REF, TIP3P], "units": "real",
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "cvff"},
       "improper_order": "center2_sorted", "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": impropers,
       "notes": ["built by CAPS (bench/ff/make_xiao2011.py) from the paper's Table 1 and energy expression",
                 "carbonate: no C-O bond term and no exclusions; the O-C-O bend and the plane term by distance (1.6 A)",
                 "Ca-O C6 = 0 in the paper, written as 1e-12 kJ mol^-1 nm^6 (sigma/epsilon form)",
                 "Lennard-Jones cut at 1 nm as the paper",
                 "atom-wise C12, C6 for mineral-peptide pairs (Table 1, eq. 3, geometric): Ca 2.52e-7, 1.42e-3; C 1.44e-5, 4.61e-3; "
                 "O 1.77e-6, 2.03e-3 (kJ/mol nm^12, nm^6)"],
       "version": "2011 (Table 1)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined",
           "description": "calcium carbonate and TIP3P water (bench/ff/make_xiao2011.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries")
for p in pairs[3:11]:
    print(p["name"], [round(x, 6) for x in p["params"]])
