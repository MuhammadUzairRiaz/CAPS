#!/usr/bin/env python3
"""DREIDING for CAPS from its paper: S. L. Mayo, B. D. Olafson, W. A. Goddard III, J. Phys. Chem. 94, 8897 (1990).

Table I gives each atom type its bond radius R, bond angle theta, and van der Waals R0 and D0; the force field is
built from them by the paper's rules:
  bonds      E = 1/2 k (R - R_IJ)^2, R_IJ = R_I + R_J - 0.01 A, k = 700 n kcal/mol/A^2 for bond order n (resonant 1.5):
             CAPS's harmonic/order style, K = 350 per unit order times the structure's perceived order
  angles     E = 1/2 C (cos theta - cos theta_J)^2, C = K / sin^2 theta_J, K = 100 kcal/mol/rad^2 (cosine/squared with
             K = 50 / sin^2 theta_J); a linear centre E = K (1 + cos theta), K = 100
  torsions   by the hybridisation of the two central atoms, cases (a)-(j), V divided over the torsions about the bond:
             CAPS's built-in rules ("torsion_rules": "dreiding1990")
  inversions planar centres (X_2, X_R) E = K (1 - cos psi), K = 40 kcal/mol, averaged over the three permutations
  van der Waals  Lennard-Jones 12-6 with R0, D0; unlike pairs by geometric means of R0 and D0
  hydrogen bonds  E = D_hb [5 (R_hb/R)^12 - 6 (R_hb/R)^10] cos^4 theta, D_hb 4.0 kcal/mol, R_hb 2.75 A (with charges)
  1-2 and 1-3 interactions excluded, 1-4 kept in full; charges Gasteiger-Marsili (CAPS's)
The typing rules (written beside it): element and hybridisation as the paper defines its labels (_3 tetrahedral,
_R resonant, _2 trigonal, _1 linear), H___A the hydrogen of a hydrogen bond (on N or O).

usage: make_dreiding.py [DATA_DIR]
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "S. L. Mayo, B. D. Olafson, W. A. Goddard III, J. Phys. Chem. 94, 8897 (1990)"

# Table I: type, element, mass, R (A), theta (deg), R0 (A), D0 (kcal/mol)
TABLE = [
    ("H_", "H", 1.008, 0.330, 180.0, 3.195, 0.0152),
    ("H___A", "H", 1.008, 0.330, 180.0, 3.195, 0.0001),
    ("H___b", "H", 1.008, 0.510, 90.0, 3.195, 0.0152),
    ("B_3", "B", 10.81, 0.880, 109.471, 4.02, 0.095),
    ("B_2", "B", 10.81, 0.790, 120.0, 4.02, 0.095),
    ("C_3", "C", 12.011, 0.770, 109.471, 3.8983, 0.0951),
    ("C_R", "C", 12.011, 0.700, 120.0, 3.8983, 0.0951),
    ("C_2", "C", 12.011, 0.670, 120.0, 3.8983, 0.0951),
    ("C_1", "C", 12.011, 0.602, 180.0, 3.8983, 0.0951),
    ("N_3", "N", 14.007, 0.702, 106.7, 3.6621, 0.0774),
    ("N_R", "N", 14.007, 0.650, 120.0, 3.6621, 0.0774),
    ("N_2", "N", 14.007, 0.615, 120.0, 3.6621, 0.0774),
    ("N_1", "N", 14.007, 0.556, 180.0, 3.6621, 0.0774),
    ("O_3", "O", 15.999, 0.660, 104.51, 3.4046, 0.0957),
    ("O_R", "O", 15.999, 0.660, 120.0, 3.4046, 0.0957),
    ("O_2", "O", 15.999, 0.560, 120.0, 3.4046, 0.0957),
    ("O_1", "O", 15.999, 0.528, 180.0, 3.4046, 0.0957),
    ("F_", "F", 18.998, 0.611, 180.0, 3.4720, 0.0725),
    ("Al3", "Al", 26.982, 1.047, 109.471, 4.39, 0.31),
    ("Si3", "Si", 28.085, 0.937, 109.471, 4.27, 0.31),
    ("P_3", "P", 30.974, 0.890, 93.3, 4.15, 0.32),
    ("S_3", "S", 32.06, 1.040, 92.1, 4.03, 0.344),
    ("Cl", "Cl", 35.45, 0.997, 180.0, 3.9503, 0.2833),
    ("Ga3", "Ga", 69.723, 1.210, 109.471, 4.39, 0.40),
    ("Ge3", "Ge", 72.63, 1.210, 109.471, 4.27, 0.40),
    ("As3", "As", 74.922, 1.210, 92.1, 4.15, 0.41),
    ("Se3", "Se", 78.971, 1.210, 90.6, 4.03, 0.43),
    ("Br", "Br", 79.904, 1.167, 180.0, 3.95, 0.37),
    ("In3", "In", 114.818, 1.390, 109.471, 4.59, 0.55),
    ("Sn3", "Sn", 118.71, 1.373, 109.471, 4.47, 0.55),
    ("Sb3", "Sb", 121.76, 1.432, 91.6, 4.35, 0.55),
    ("Te3", "Te", 127.6, 1.280, 90.3, 4.23, 0.57),
    ("I_", "I", 126.904, 1.360, 180.0, 4.15, 0.51),
    ("Na", "Na", 22.99, 1.860, 90.0, 3.144, 0.5),
    ("Ca", "Ca", 40.078, 1.940, 90.0, 3.472, 0.05),
    ("Fe", "Fe", 55.845, 1.285, 90.0, 4.54, 0.055),
    ("Zn", "Zn", 65.38, 1.330, 109.471, 4.54, 0.055),
]
SIG = 2 ** (1 / 6)

types, pairs, bonds, angles = [], [], [], []
for t, el, m, R, th, R0, D0 in TABLE:
    types.append({"name": t, "element": el, "mass": m, "description": f"DREIDING {t}: R {R} A, theta {th}, R0 {R0} A, D0 {D0} kcal/mol"})
    pairs.append({"name": t, "match": [t], "params": [D0, R0 / SIG]})
    if abs(th - 180) < 1e-9:
        angles.append({"name": f"X-{t}-X", "match": ["*", t, "*"], "style": "cosine", "params": [100.0], "comment": "linear centre: K (1 + cos theta)"})
    else:
        s2 = math.sin(math.radians(th)) ** 2
        angles.append({"name": f"X-{t}-X", "match": ["*", t, "*"], "style": "cosine/squared", "params": [round(50.0 / s2, 8), th],
                       "comment": "1/2 C (cos theta - cos theta0)^2, C = 100 / sin^2 theta0"})
for a in range(len(TABLE)):
    for b in range(a, len(TABLE)):
        ta, tb = TABLE[a], TABLE[b]
        bonds.append({"name": f"{ta[0]}-{tb[0]}", "match": [ta[0], tb[0]], "style": "harmonic/order",
                      "params": [350.0, round(ta[3] + tb[3] - 0.01, 6)]})
impropers = [{"name": f"{c} planar centre", "match": [c, "*", "*", "*"], "style": "planar", "params": [40.0, 0.0],
              "comment": "DREIDING planar inversion, K (1 - cos psi)"} for c in ("*_2", "*_R")]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "DREIDING (Mayo 1990)", "source": REF + ", Table I and its rules",
       "references": [REF], "units": "real",
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "cosine/squared", "dihedral": "charmm", "improper": "umbrella",
                  "special": "dreiding"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10, "torsion_rules": "dreiding1990",
       "improper_order": "center1_sorted", "improper_max_neighbours": 3, "equivalence": "replace",
       "hbonds": {"style": "dreiding/lj", "power": 4, "inner": 6.0, "outer": 6.5, "angle": 90,
                  "source": REF + ": E_hb = D_hb [5 (R_hb/R)^12 - 6 (R_hb/R)^10] cos^4 theta, D_hb 4.0 kcal/mol, R_hb 2.75 A (with charges)",
                  "terms": [{"donor": d, "acceptor": a, "hydrogen": "H___A", "eps": 4.0, "sigma": 2.75, "n": 4}
                            for d in ("N_*", "O_*") for a in ("N_*", "O_*")]},
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": impropers,
       "notes": ["built by CAPS (bench/ff/make_dreiding.py) from the paper's Table I and rules"],
       "version": "1990", "typing": "../typing/dreiding1990.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", "dreiding1990.json"), "w"), indent=1)

R = []
def add(t, p, s, d):
    R.append({"type": t, "smarts": s, "priority": p, "description": d})
add("H_", 0, "[#1]", "hydrogen"); add("H___A", 2, "[H][#7,#8]", "hydrogen-bond hydrogen (on N or O)")
add("H___b", 3, "[H;D2]([#5])[#5]", "bridging hydrogen (diborane)")
add("B_3", 0, "[#5X4]", "tetrahedral B"); add("B_2", 0, "[#5X3]", "trigonal B")
add("C_3", 0, "[#6X4]", "tetrahedral C"); add("C_2", 0, "[#6X3]", "trigonal C"); add("C_1", 0, "[#6X2]", "linear C")
add("C_R", 1, "[c]", "resonant (aromatic) C")
add("N_3", 0, "[#7X4,#7X3]", "tetrahedral N"); add("N_2", 1, "[#7X2]", "trigonal N"); add("N_1", 1, "[#7X1]", "linear N")
add("N_R", 2, "[n]", "resonant (aromatic) N")
add("N_R", 2, "[#7X3;$(*-[#6X3]=[#8X1,#16X1]),$(*-c)]", "resonant N (amide, N on an aromatic ring)")
add("O_3", 0, "[#8X2]", "tetrahedral O"); add("O_2", 1, "[#8X1]", "trigonal O"); add("O_R", 2, "[o]", "resonant (aromatic) O")
add("F_", 0, "[F]", "fluorine"); add("Cl", 0, "[Cl]", "chlorine"); add("Br", 0, "[Br]", "bromine"); add("I_", 0, "[I]", "iodine")
for t, z in (("Al3", 13), ("Si3", 14), ("P_3", 15), ("S_3", 16), ("Ga3", 31), ("Ge3", 32), ("As3", 33), ("Se3", 34), ("In3", 49),
             ("Sn3", 50), ("Sb3", 51), ("Te3", 52), ("Na", 11), ("Ca", 20), ("Fe", 26), ("Zn", 30)):
    add(t, 0, f"[#{z}]", t)
json.dump({"format": "caps-typing", "version": 1, "forcefield": "DREIDING (Mayo 1990)", "unknown_types": "untyped",
           "description": "CAPS rules for DREIDING's labels as the paper defines them: element and hybridisation (_3 tetrahedral, "
                          "_R resonant, _2 trigonal, _1 linear), H___A the hydrogen of a hydrogen bond (bench/ff/make_dreiding.py)",
           "rules": R}, open(os.path.join(DATA, "typing", "dreiding1990.typing.json"), "w"), indent=1)
print(f"dreiding1990: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles; {len(R)} typing rules")
