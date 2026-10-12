#!/usr/bin/env python3
"""Oxides for CAPS from G. V. Lewis, C. R. A. Catlow, J. Phys. C: Solid State Phys. 18, 1149 (1985), Table 3 (shell
model parameters, method (c)(i)), with the O-O potential of C. R. A. Catlow, Proc. R. Soc. Lond. A 353, 533 (1977),
Table 3: A22 and rho22 (the same in both of Catlow's potentials) with C22 = 20.37 eV A^6, potential 2's (potential 1 has
112.2). Lewis & Catlow name only "the potential derived by Catlow (1977)"; 20.37 is the value the published transcriptions
of this set use.

  cation-O   A exp(-r/rho) between the cation (its shell when it has one) and the oxide's O shell
  O-O        22764.3 exp(-r/0.1490) - 20.37/r^6 eV between O shells
  cation-cation Coulomb only (assumption (ii) of the paper)
  shells     cation shell Y+, spring k+ (k+ = 99999: unpolarisable, a rigid ion of formal charge); the oxide's own O
             shell Y-, k- (core +Z_O - Y-); springs 1/2 k r^2 (class II bond K2 = k/2)
Each oxide is a set typed by composition. A bracketed A (a cation in a tetrahedral oxygen environment) and the second
oxidation state of an element (Fe3+ beside Fe2+, Ti4+ beside Ti3+) are overlays: data/typing/sets/lewis-catlow1985-oxides/.
Tables 1, 2 (short-range terms only) and 4 (shell columns not self-consistent with alpha = Y^2/k) are not written.

usage: make_lewis_catlow1985.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "G. V. Lewis, C. R. A. Catlow, J. Phys. C: Solid State Phys. 18, 1149 (1985)"
REF_OO = "C. R. A. Catlow, Proc. R. Soc. Lond. A 353, 533 (1977)"
EV = 23.060548
FID = "lewis-catlow1985-oxides"
# Table 3: cation, charge Z, A (eV), A tetrahedral (or None), rho (A), Y+, k+ (eV/A^2), Y- (O shell), k- (eV/A^2)
T3 = [
    ("Mg", 2, 821.6, 710.5, 0.3242, 2.0, 99999, -2.00, 15.74),
    ("Ca", 2, 1228.9, None, 0.3372, 1.26, 34.00, -2.82, 27.45),
    ("Sr", 2, 1400.0, None, 0.3500, 1.33, 21.53, -3.51, 39.56),
    ("Ba", 2, 931.7, None, 0.3949, 1.46, 14.78, -4.84, 70.39),
    ("Mn", 2, 715.8, None, 0.3464, 3.00, 81.20, -2.30, 19.73),
    ("Fe", 2, 694.1, 599.4, 0.3399, 2.00, 10.92, -2.15, 17.75),
    ("Co", 2, 696.3, 600.9, 0.3362, 2.00, 10.74, -2.06, 16.52),
    ("Ni", 2, 683.5, 593.9, 0.3332, 2.00, 8.77, -2.00, 15.78),
    ("Zn", 2, 499.6, None, 0.3595, 2.05, 10.28, -2.00, 15.52),
    ("Al", 3, 1114.9, 1012.6, 0.3118, 3.00, 99999, -2.21, 27.29),
    ("Ti", 3, 1715.7, None, 0.3069, 1.23, 43.97, -2.21, 27.29),
    ("V", 3, 1790.2, None, 0.3061, 1.34, 36.27, -2.21, 27.29),
    ("Cr", 3, 1734.1, 1526.4, 0.3010, 0.97, 67.00, -2.21, 27.29),
    ("Fe", 3, 1102.4, 876.6, 0.3299, 4.97, 304.7, -2.21, 27.29),
    ("Ge", 4, 1035.5, None, 0.3464, 4.88, 317.2, -2.09, 61.5),
    ("Sn", 4, 938.7, None, 0.3813, 3.99, 191.5, -3.00, 117.8),
    ("Ti", 4, 754.2, None, 0.3879, 2.89, 37.3, -2.53, 86.4),
    ("Ce", 4, 1013.6, None, 0.3949, 5.85, 116.0, -2.67, 51.0),
    ("U", 4, 1055.5, None, 0.3949, 6.44, 129.0, -2.70, 51.6),
    ("Th", 4, 1147.7, None, 0.3949, 7.28, 193.1, -2.40, 39.0),
]
OO = (22764.3, 0.1490, 20.37)
MASS = {"Mg": 24.305, "Ca": 40.078, "Sr": 87.62, "Ba": 137.327, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Zn": 65.38,
        "Al": 26.982, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Ge": 72.63, "Sn": 118.71, "Ce": 140.116, "U": 238.029, "Th": 232.038, "O": 15.999}
Z = {"Mg": 12, "Ca": 20, "Sr": 38, "Ba": 56, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Zn": 30, "Al": 13, "Ti": 22, "V": 23, "Cr": 24,
     "Ge": 32, "Sn": 50, "Ce": 58, "U": 92, "Th": 90, "O": 8}
ALL = sorted({r[0] for r in T3})
types, pairs, bonds, shells = [], [], [], {}
default_rules, overlays = [], {}
seen_el = {}
for el, zc, A, Atet, rho, Yp, kp, Ym, km in T3:
    for tet in ((False, True) if Atet else (False,)):
        tag = f"{el}{zc}" + ("t" if tet else "")
        a = Atet if tet else A
        rigid = kp >= 99999
        c_core = f"{tag}c" if not rigid else tag
        o_core, o_shell = f"O_{tag}c", f"O_{tag}s"
        ctypes = [{"name": c_core, "element": el, "mass": MASS[el] * (1.0 if rigid else 0.9), "charge": float(zc) if rigid else zc - Yp,
                   "description": f"{el}{zc}+ {'(rigid: k+ 99999)' if rigid else 'core'}{' (tetrahedral)' if tet else ''}; {REF}, Table 3"}]
        if not rigid:
            ctypes.append({"name": f"{tag}s", "element": el, "mass": MASS[el] * 0.1, "charge": Yp, "description": f"{el}{zc}+ shell; Table 3"})
            shells[c_core] = f"{tag}s"
            bonds.append({"name": f"{c_core}-{tag}s", "match": [c_core, f"{tag}s"], "style": "class2", "params": [0.0, kp / 2 * EV, 0.0, 0.0],
                          "comment": f"cation core-shell spring k+ {kp} eV/A^2 (Table 3)"})
        ctypes += [{"name": o_core, "element": "O", "mass": MASS["O"] * 0.9, "charge": -2.0 - Ym, "description": f"O2- core in the {tag} oxide; Table 3"},
                   {"name": o_shell, "element": "O", "mass": MASS["O"] * 0.1, "charge": Ym, "description": f"O2- shell in the {tag} oxide; Table 3"}]
        shells[o_core] = o_shell
        bonds.append({"name": f"{o_core}-{o_shell}", "match": [o_core, o_shell], "style": "class2", "params": [0.0, km / 2 * EV, 0.0, 0.0],
                      "comment": f"O core-shell spring k- {km} eV/A^2 (Table 3)"})
        types += ctypes
        csite = c_core if rigid else f"{tag}s"
        pairs += [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in ctypes]
        pairs += [{"name": f"{csite}-{o_shell}", "match": [csite, o_shell], "style": "buck", "params": [a * EV, rho, 0.0],
                   "comment": f"{el}{zc}+-O2- A {a} eV, rho {rho} A (Table 3{', tetrahedral' if tet else ''})"},
                  {"name": f"{o_shell}-{o_shell}", "match": [o_shell, o_shell], "style": "buck", "params": [OO[0] * EV, OO[1], OO[2] * EV],
                   "comment": f"O-O {OO[0]} eV, {OO[1]} A, C {OO[2]} eV A^6 ({REF_OO}, Table 3, potential 2's C)"}]
        others = [e for e in ALL if e != el]
        rs = [{"type": c_core, "smarts": f"[#{Z[el]}]", "priority": 0, "requires": ["O"], "excludes": others, "description": f"{el}{zc}+ in its oxide"},
              {"type": o_core, "smarts": "[#8]", "priority": 0, "requires": [el], "excludes": others, "description": f"O2- in the {el}{zc}+ oxide"}]
        if tet or el in seen_el:
            overlays[tag] = rs
        else:
            default_rules += rs
            seen_el[el] = tag
rules = default_rules + [{"type": sh, "smarts": "[#8]" if sh.startswith("O_") else f"[#{Z[[t for t in types if t['name'] == sh][0]['element']]}]",
                          "priority": 99, "atom_name": sh, "description": f"shell of {c}"} for c, sh in shells.items()]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Oxides (Lewis & Catlow 1985)", "source": REF + ", Table 3; O-O: " + REF_OO + ", Table 3",
       "references": [REF, REF_OO], "units": "real",
       "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 6,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_lewis_catlow1985.py) from the paper's Table 3 (method (c)(i))",
                 "short-range terms truncated at 6 A as the paper's (Table 1 footnote)", "O-O: Catlow 1977 Table 3, A 22764.3 eV, rho 0.1490 A, C 20.37 eV A^6 (potential 2's C; Lewis & Catlow do not say which)",
                 "shells carry 10 % of the ion's mass (adiabatic core-shell)"],
       "version": "1985 (Table 3)", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
shell_rules = [r for r in rules if r.get("atom_name")]
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "shells": shells,
           "description": "each oxide typed by composition (bench/ff/make_lewis_catlow1985.py); tetrahedral sites and second oxidation "
                          "states are overlays", "rules": rules}, open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
os.makedirs(os.path.join(DATA, "typing", "sets", FID), exist_ok=True)
for name, rs in overlays.items():
    json.dump({"format": "caps-typing", "version": 1, "forcefield": f"{FID}: {name}", "unknown_types": "untyped",
               "description": f"Overlay: the {name} set (Table 3) in place of the element's first", "rules": [dict(r, priority=50) for r in rs]},
              open(os.path.join(DATA, "typing", "sets", FID, name + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(bonds)} springs, {len(default_rules)} default rules, {len(overlays)} overlays")
