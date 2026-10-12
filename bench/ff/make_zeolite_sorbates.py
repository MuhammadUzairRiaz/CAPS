#!/usr/bin/env python3
"""Sorbate models for zeolites, each from its paper, for CAPS:

  omalley2013-alkanes  A. J. O'Malley, C. R. A. Catlow, Phys. Chem. Chem. Phys. 15, 19024 (2013), Tables 2 and 3:
                       all-atom n-alkanes in silicalite. H +0.1, each C the charge that makes the alkane neutral (CH3 -0.3,
                       CH2 -0.2, CH -0.1, C 0, CH4 -0.4); bonds 1/2 K (r - r0)^2 and bends 1/2 K (theta - theta0)^2 as the
                       paper's DL_POLY input reads them (the table prints the bond K in kJ mol-1 rad-2: per A^2 here);
                       Lennard-Jones A / r^12 - B / r^6 for C-C, C-H, H-H and the framework O (as eps = B^2 / 4A, sigma =
                       (A / B)^(1/6)); the framework from Table 3 (formal charges, Buckingham Si-O and O-O within 10 A,
                       O-Si-O bend). The torsion is in Raj, Sastre & Catlow 1999 (ref. 18 of the paper), not given here:
                       the library holds a zero torsion, and 1-4 pairs count in full. The united-atom model of Table 1 is not included.
  plant2007-methanol   D. F. Plant, G. Maurin, R. G. Bell, J. Phys. Chem. B 111, 2836 (2007), Tables 1-3: methanol in NaY.
                       Framework and Na as Ramsahye & Bell (Table 1: Kramer charges, Buckingham, O-T-O bends); methanol's
                       charges, bonds, bends, H-C-O-H torsion k [1 + cos 3 phi] and Lennard-Jones pairs (Table 3); its pairs
                       with the framework O and with Na (Table 2: Lennard-Jones, Na-Om Buckingham). Pairs no table lists:
                       Coulomb only.
  peguin2009-hfa134a   R. P. S. Peguin, G. Kamath, J. J. Potoff, S. R. P. da Rocha, J. Phys. Chem. B 113, 178 (2009), FF no. 2:
                       HFA134a (CF3CH2F). Table 1 (read from the page image): OPLS-AA's fixed bond lengths (stiff harmonic
                       here, to be constrained), bends k_theta / 2 (theta - theta0)^2 in K / rad^2, F-C-C-F and F-C-C-H
                       torsions C0 + C1 (1 + cos phi) + C2 (1 - cos 2 phi) + C3 (1 + cos 3 phi) (eq. 3; the table's four
                       columns read as C0..C3); Table 2: charges and Lennard-Jones of FF no. 2, Lorentz-Berthelot; eq. 1 has no
                       intramolecular non-bonded energy, so every intramolecular pair is excluded. Cut-off 14 A (section 3.2).

usage: make_zeolite_sorbates.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
EV, KJ, KB = 23.060548, 1 / 4.184, 0.0019872043   # eV, kJ/mol, K -> kcal/mol
STYLES = {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"}


def fill_pairs(types, pairs):
    """every pair no table lists: Coulomb only (explicit zero, so no mixing rule invents one)"""
    listed = {tuple(sorted(p["match"] * (2 if len(p["match"]) == 1 else 1))) for p in pairs}
    names = [t["name"] for t in types]
    for i, a in enumerate(names):
        for b in names[i:]:
            if tuple(sorted((a, b))) not in listed:
                pairs.append({"name": a if a == b else f"{a}-{b}", "match": [a] if a == b else [a, b], "params": [0.0, 0.0],
                              "comment": "not in the paper's tables: Coulomb only"})
    return pairs


def write(fid, doc, rules, desc):
    doc = {"format": "caps-forcefield", "format_version": 1, **doc, "units": "real", "typing": f"../typing/{fid}.typing.json"}
    json.dump(doc, open(os.path.join(DATA, "forcefields", fid + ".json"), "w"), indent=1, ensure_ascii=False)
    json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined",
               "description": desc + " (bench/ff/make_zeolite_sorbates.py)", "rules": rules},
              open(os.path.join(DATA, "typing", fid + ".typing.json"), "w"), indent=1, ensure_ascii=False)
    print(f"{fid}: {len(doc['atom_types'])} types, {len(doc['pairs'])} pairs, {len(doc['bonds'])} bonds, {len(doc['angles'])} bends, "
          f"{len(doc['dihedrals'])} torsions")


# ------------------------------------------------------------------ O'Malley & Catlow 2013, all-atom alkanes in silicalite
REF = "A. J. O'Malley, C. R. A. Catlow, Phys. Chem. Chem. Phys. 15, 19024 (2013)"
FID = "omalley2013-alkanes"
C_KIND = [("C4_om", 4, -0.4, "methane C"), ("C3_om", 3, -0.3, "CH3 C"), ("C2_om", 2, -0.2, "CH2 C"), ("C1_om", 1, -0.1, "CH C"),
          ("C0_om", 0, 0.0, "quaternary C")]
types = [{"name": n, "element": "C", "mass": 12.011, "charge": q, "description": f"alkane {d}: H +0.1 each, the alkane neutral",
          "equivalence": {"vdw": "C_om", "bond": "C_om", "angle": "C_om", "dihedral": "C_om"}} for n, h, q, d in C_KIND]
types += [{"name": "C_om", "element": "C", "mass": 12.011, "charge": 0.0, "description": "alkane C (the class every alkane C's terms take)"},
          {"name": "H_om", "element": "H", "mass": 1.008, "charge": 0.1, "description": "alkane H (+0.1 a.u.)"},
          {"name": "Si_om", "element": "Si", "mass": 28.0855, "charge": 4.0, "description": "framework Si, formal charge (Table 3)"},
          {"name": "O_om", "element": "O", "mass": 15.999, "charge": -2.0, "description": "framework O, formal charge (Table 3)"}]


def ab(a, b, A, B):   # A / r^12 - B / r^6 (kJ/mol) as 4 eps [(s/r)^12 - (s/r)^6]
    eps, sig = B * B / (4 * A) * KJ, (A / B) ** (1 / 6)
    return {"name": a if a == b else f"{a}-{b}", "match": [a] if a == b else [a, b], "params": [eps, sig],
            "comment": f"A {A:g} kJ/mol A^12, B {B:g} kJ/mol A^6 (Table 2)"}


pairs = [ab("C_om", "H_om", 270160, 563.62), ab("H_om", "H_om", 37131.41, 191.69), ab("C_om", "C_om", 1899990, 1746.32),
         ab("H_om", "O_om", 150169.7, 537.5872), ab("C_om", "O_om", 1061340, 1703.35)]
pairs += [{"name": f"{a}-{b}", "match": [a, b], "style": "buck", "params": [A * KJ, rho, C * KJ, 10.0],
           "comment": f"A {A:g} kJ/mol, rho {rho} A, C {C} kJ/mol A^6 (Table 3); within 10 A"}
          for a, b, A, rho, C in (("Si_om", "O_om", 123878.16, 0.149, 2690.01), ("O_om", "O_om", 2196390, 0.32052, 1028.53))]
pairs = fill_pairs([t for t in types if not t.get("equivalence")], pairs)
bonds = [{"name": "C_om-H_om", "match": ["C_om", "H_om"], "params": [2990 / 2 * KJ, 1.054], "comment": "1/2 K (r - r0)^2, K 2990 kJ/mol/A^2 (Table 2)"},
         {"name": "C_om-C_om", "match": ["C_om", "C_om"], "params": [5790 / 2 * KJ, 1.405], "comment": "K 5790 kJ/mol/A^2, r0 1.405 A as Table 2 prints it"}]
angles = [{"name": "-".join(m), "match": list(m), "params": [327 / 2 * KJ, 109.47], "comment": "1/2 K (theta - theta0)^2, K 327 kJ/mol/rad^2 (Table 2)"}
          for m in (("C_om", "C_om", "C_om"), ("C_om", "C_om", "H_om"), ("H_om", "C_om", "H_om"))]
angles.append({"name": "O_om-Si_om-O_om", "match": ["O_om", "Si_om", "O_om"], "params": [202.353 / 2 * KJ, 109.47], "contact": [2.0],
               "comment": "1/2 K (theta - theta0)^2, K 202.353 kJ/mol/rad^2 (Table 3); each O within 2.0 A of the Si"})
rules = [{"type": n, "smarts": f"[CX4H{h}]" if h else "[CX4H0;!$(C~[!#6;!#1])]", "priority": 0, "description": f"alkane {d}"} for n, h, q, d in C_KIND]
rules += [{"type": "H_om", "smarts": "[#1][CX4]", "priority": 0, "description": "alkane H"},
          {"type": "Si_om", "smarts": "[#14]", "priority": 0, "requires": ["O"], "description": "framework Si"},
          {"type": "O_om", "smarts": "[#8]", "priority": 0, "description": "framework O"}]
write(FID, {"name": "All-atom alkanes in silicalite (O'Malley & Catlow 2013)", "source": REF + ", Tables 2 and 3", "references": [REF],
            "version": "2013 (Tables 2 and 3, all-atom model)", "styles": STYLES, "mixing": "geometric", "special_lj": [0, 0, 1],
            "special_coul": [0, 0, 1], "cutoff": 13, "equivalence": "replace", "atom_types": types, "pairs": pairs, "bonds": bonds,
            "angles": angles, "impropers": [],
            "dihedrals": [{"name": "X-C_om-C_om-X", "match": ["*", "C_om", "C_om", "*"], "style": "fourier", "params": [1, 0.0, 1, 0.0],
                           "comment": "zero: the paper's cosine torsion is in Raj, Sastre & Catlow 1999 (its ref. 18)"}],
            "notes": ["built by CAPS (bench/ff/make_zeolite_sorbates.py) from the paper's Tables 2 and 3",
                      "the paper's torsion is that of N. Raj, G. Sastre, C. R. A. Catlow, J. Phys. Chem. B 103, 11007 (1999) (its ref. 18), "
                      "not given in this paper: a zero torsion here; 1-2 and 1-3 pairs excluded, 1-4 pairs in full",
                      "bond K printed in kJ mol-1 rad-2 in Table 2: read as kJ mol-1 A-2 with DL_POLY's 1/2 K",
                      "framework: formal charges, Buckingham within 10 A, O-Si-O bends by distance (2.0 A, CAPS's); no Si-O bond term",
                      "Lennard-Jones cut-off 13 A (the paper's)"]},
      rules, "all-atom alkanes and the silicalite framework")

# ------------------------------------------------------------------ Plant, Maurin & Bell 2007, methanol in NaY
REF = "D. F. Plant, G. Maurin, R. G. Bell, J. Phys. Chem. B 111, 2836 (2007)"
FID = "plant2007-methanol"
types = [{"name": "Si_pm", "element": "Si", "mass": 28.0855, "charge": 2.4, "description": "framework Si (Table 1)"},
         {"name": "Al_pm", "element": "Al", "mass": 26.982, "charge": 1.4, "description": "framework Al (Table 1)"},
         {"name": "Oz_pm", "element": "O", "mass": 15.999, "charge": -1.2, "description": "framework O (Table 1)"},
         {"name": "Na_pm", "element": "Na", "mass": 22.98977, "charge": 1.0, "description": "extra-framework Na (Table 1)"},
         {"name": "C_pm", "element": "C", "mass": 12.011, "charge": -0.093, "description": "methanol C (Table 3)"},
         {"name": "H_pm", "element": "H", "mass": 1.008, "charge": 0.100, "description": "methanol methyl H (Table 3)"},
         {"name": "Om_pm", "element": "O", "mass": 15.999, "charge": -0.432, "description": "methanol O (Table 3)"},
         {"name": "Hm_pm", "element": "H", "mass": 1.008, "charge": 0.225, "description": "methanol hydroxyl H (Table 3)"}]
pairs = [{"name": f"{a}-{b}", "match": [a, b], "style": "buck", "params": [A * EV, rho, C * EV], "comment": f"A {A} eV, rho {rho} A, C {C} eV A^6 (Table {t})"}
         for a, b, A, rho, C, t in (("Si_pm", "Oz_pm", 30023.0, 0.1621, 12.84, 1), ("Al_pm", "Oz_pm", 26998.0, 0.1622, 12.84, 1),
                                    ("Na_pm", "Oz_pm", 8200.0, 0.2180, 11.80, 1), ("Oz_pm", "Oz_pm", 894.6, 0.3244, 0.0, 1),
                                    ("Na_pm", "Om_pm", 4521.60, 0.2177, 0.0, 2))]
pairs += [{"name": a if a == b else f"{a}-{b}", "match": [a] if a == b else [a, b], "params": [e * EV, s], "comment": f"eps {e} eV, sigma {s} A (Table {t})"}
          for a, b, e, s, t in (("Om_pm", "Oz_pm", 0.00762, 4.104, 2), ("Hm_pm", "Oz_pm", 0.00512, 3.374, 2), ("C_pm", "Oz_pm", 0.00591, 4.310, 2),
                                ("H_pm", "Oz_pm", 0.00387, 3.732, 2), ("Na_pm", "Hm_pm", 0.00878, 2.6285, 2), ("Na_pm", "H_pm", 0.00878, 2.6285, 2),
                                ("Na_pm", "C_pm", 0.01340, 3.0722, 2), ("H_pm", "H_pm", 0.00165, 2.450, 3), ("C_pm", "H_pm", 0.00338, 2.920, 3),
                                ("Om_pm", "H_pm", 0.00404, 2.650, 3), ("C_pm", "C_pm", 0.00694, 3.475, 3), ("Om_pm", "C_pm", 0.00828, 3.150, 3),
                                ("Om_pm", "Om_pm", 0.00988, 2.860, 3))]
pairs = fill_pairs(types, pairs)
bonds = [{"name": f"{a}-{b}", "match": [a, b], "params": [k / 2 * EV, r0], "comment": f"1/2 k (r - r0)^2, k {k} eV/A^2 (Table 3)"}
         for a, b, k, r0 in (("C_pm", "H_pm", 29.56, 1.105), ("C_pm", "Om_pm", 33.33, 1.420), ("Hm_pm", "Om_pm", 46.97, 0.945))]
angles = [{"name": f"{a}-{b}-{c}", "match": [a, b, c], "params": [k / 2 * EV, t0], "comment": f"1/2 k (theta - theta0)^2, k {k} eV/rad^2 (Table 3)"}
          for a, b, c, k, t0 in (("C_pm", "Om_pm", "Hm_pm", 5.60, 108.32), ("H_pm", "C_pm", "Om_pm", 5.50, 106.90), ("H_pm", "C_pm", "H_pm", 4.40, 108.38))]
angles += [{"name": f"Oz_pm-{t}-Oz_pm", "match": ["Oz_pm", t, "Oz_pm"], "params": [k / 2 * EV, 109.47], "contact": [2.0],
            "comment": f"1/2 k (theta - theta0)^2, k {k} eV/rad^2 (Table 1); each O within 2.0 A of the T"} for t, k in (("Si_pm", 12.10), ("Al_pm", 2.20))]
dihedrals = [{"name": "H_pm-C_pm-Om_pm-Hm_pm", "match": ["H_pm", "C_pm", "Om_pm", "Hm_pm"], "params": [1, 0.00762 * EV, 3, 0.0],
              "comment": "k [1 + cos 3 phi], k 0.00762 eV (Table 3)"}]
MEOH = "[CX4H3][OX2H1]"
rules = [{"type": "Si_pm", "smarts": "[#14]", "priority": 0, "requires": ["O"], "description": "framework Si"},
         {"type": "Al_pm", "smarts": "[#13]", "priority": 0, "requires": ["O"], "description": "framework Al"},
         {"type": "Oz_pm", "smarts": "[#8]", "priority": 0, "description": "framework O"},
         {"type": "Na_pm", "smarts": "[#11]", "priority": 0, "description": "extra-framework Na"},
         {"type": "C_pm", "smarts": f"[#6;$({MEOH})]", "priority": 1, "description": "methanol C"},
         {"type": "Om_pm", "smarts": f"[#8;$([OX2H1][CX4H3])]", "priority": 1, "description": "methanol O"},
         {"type": "H_pm", "smarts": f"[#1;$([#1][CX4H3][OX2H1])]", "priority": 1, "description": "methanol methyl H"},
         {"type": "Hm_pm", "smarts": f"[#1;$([#1][OX2][CX4H3])]", "priority": 1, "description": "methanol hydroxyl H"}]
write(FID, {"name": "Methanol in zeolite NaY (Plant, Maurin & Bell 2007)", "source": REF + ", Tables 1-3", "references": [REF],
            "version": "2007 (Tables 1-3)", "styles": STYLES, "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1],
            "cutoff": 10, "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": [],
            "notes": ["built by CAPS (bench/ff/make_zeolite_sorbates.py) from the paper's Tables 1-3",
                      "framework as Ramsahye & Bell 2005 (Table 1): no T-O bond terms, O-T-O bends by distance (2.0 A, CAPS's)",
                      "methanol's 1-4 pairs (H...Hm) in full: the paper names no scaling", "every pair no table lists: Coulomb only",
                      "cut-off 10 A (the paper gives none)"]},
      rules, "NaY framework, Na and methanol")

# ------------------------------------------------------------------ Peguin et al. 2009, HFA134a (FF no. 2)
REF = "R. P. S. Peguin, G. Kamath, J. J. Potoff, S. R. P. da Rocha, J. Phys. Chem. B 113, 178 (2009)"
FID = "peguin2009-hfa134a"
CF3 = "[CX4;$(C(F)(F)(F)[CX4H2]F)]"
CH2F = "[CX4;$([CX4H2](F)C(F)(F)F)]"
types = [{"name": "C3_hfa", "element": "C", "mass": 12.011, "charge": 0.5340, "description": "CF3 carbon (Table 2, FF no. 2)", "equivalence": {"vdw": "C_hfa", "bond": "C_hfa", "angle": "C_hfa", "dihedral": "C_hfa"}},
         {"name": "C2_hfa", "element": "C", "mass": 12.011, "charge": 0.0020, "description": "CH2F carbon (Table 2, FF no. 2)", "equivalence": {"vdw": "C_hfa", "bond": "C_hfa", "angle": "C_hfa", "dihedral": "C_hfa"}},
         {"name": "C_hfa", "element": "C", "mass": 12.011, "charge": 0.0, "description": "HFA134a carbon (the class both carbons' terms take)"},
         {"name": "F_hfa", "element": "F", "mass": 18.998, "charge": -0.1915, "description": "F (Table 2, FF no. 2)"},
         {"name": "H_hfa", "element": "H", "mass": 1.008, "charge": 0.1150, "description": "H (Table 2, FF no. 2)"}]
pairs = [{"name": n, "match": [n], "params": [e * KB, s], "comment": f"eps/kB {e} K, sigma {s} A (Table 2)"}
         for n, e, s in (("C_hfa", 47.0, 3.6), ("H_hfa", 10.0, 2.5), ("F_hfa", 24.5, 2.92))]
KBOND = 600.0
bonds = [{"name": f"C_hfa-{b}", "match": ["C_hfa", b], "params": [KBOND, r0], "comment": f"fixed length {r0} A (Table 1; constrain in MD)"}
         for b, r0 in (("C_hfa", 1.529), ("F_hfa", 1.360), ("H_hfa", 1.090))]
angles = [{"name": f"{a}-C_hfa-{c}", "match": [a, "C_hfa", c], "params": [k / 2 * KB, t0], "comment": f"k_theta/2 (theta - theta0)^2, k_theta {k} K/rad^2 (Table 1)"}
          for a, c, t0, k in (("C_hfa", "H_hfa", 110.7, 18889.4), ("C_hfa", "F_hfa", 109.5, 25185.8), ("H_hfa", "H_hfa", 107.8, 16622.6),
                              ("H_hfa", "F_hfa", 107.0, 20148.7), ("F_hfa", "F_hfa", 109.1, 38786.2))]


def tors(a, d, c0, c1, c2, c3):
    p = [(c1 * KB, 1, 0.0), (c2 * KB, 2, 180.0), (c3 * KB, 3, 0.0), (c0 * KB, 0, 90.0)]
    return {"name": f"{a}-C_hfa-C_hfa-{d}", "match": [a, "C_hfa", "C_hfa", d], "params": [len(p)] + [x for t in p for x in t],
            "comment": f"C0 + C1 (1 + cos phi) + C2 (1 - cos 2 phi) + C3 (1 + cos 3 phi); C0..C3 {c0}, {c1}, {c2}, {c3} K (Table 1)"}


dihedrals = [tors("F_hfa", "F_hfa", 11.35, 1.96, -24.50, 1070.18), tors("F_hfa", "H_hfa", -29.32, 1.59, 9.96, 1064.65)]
rules = [{"type": "C3_hfa", "smarts": CF3, "priority": 0, "description": "HFA134a CF3 carbon"},
         {"type": "C2_hfa", "smarts": CH2F, "priority": 0, "description": "HFA134a CH2F carbon"},
         {"type": "F_hfa", "smarts": f"[F;$(F{CF3}),$(F{CH2F})]", "priority": 0, "description": "HFA134a F"},
         {"type": "H_hfa", "smarts": f"[#1;$([#1]{CH2F})]", "priority": 0, "description": "HFA134a H"}]
write(FID, {"name": "HFA134a (Peguin et al. 2009, FF no. 2)", "source": REF + ", Tables 1 and 2 (FF no. 2)", "references": [REF],
            "version": "2009, FF no. 2", "styles": STYLES, "mixing": "arithmetic", "special_lj": [0, 0, 0], "special_coul": [0, 0, 0],
            "cutoff": 14, "equivalence": "replace", "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles,
            "dihedrals": dihedrals, "impropers": [],
            "notes": ["built by CAPS (bench/ff/make_zeolite_sorbates.py) from the paper's Tables 1 (read from the page image) and 2",
                      "bond lengths fixed in the paper: stiff harmonic bonds (600 kcal/mol/A^2), to be constrained in MD",
                      "Table 1's torsion columns read as C0..C3 of eq. 3", "no intramolecular non-bonded pairs (eq. 1)",
                      "charges belong to HFA134a only (CHELPG fit of that molecule)", "cut-off 14 A (section 3.2)"]},
      rules, "HFA134a")
