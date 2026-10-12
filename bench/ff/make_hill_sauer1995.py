#!/usr/bin/env python3
"""Aluminosilicates for CAPS from their paper: J.-R. Hill, J. Sauer, J. Phys. Chem. 99, 9536 (1995), Table 8 (the cross terms
read from the page image), in the CFF form of the companion paper (J.-R. Hill, J. Sauer, J. Phys. Chem. 98, 1238 (1994), eq. 1-10):

  bonds       K2 (b - b0)^2 + K3 (b - b0)^3 + K4 (b - b0)^4                      (LAMMPS bond class2)
  angles      H2 (theta - theta0)^2 + H3 (...)^3 + H4 (...)^4                     (angle class2)
  torsions    V1 [1 - cos phi] + V2 [1 - cos 2phi] + V3 [1 - cos 3phi]            (dihedral class2)
  cross       bond-bond F_bb' (b - b0)(b' - b0'), bond-angle F_btheta, F_b'theta (a blank F_b'theta: the symmetric angle's own
              F_btheta), angle-angle F_theta theta' for i-j-k/l (the angles i-j-k and k-j-l: msi2lmp's angle-angle entry
              i j k l), angle-angle-torsion K_theta theta' phi (theta - theta0)(theta' - theta0') cos phi
  out-of-plane  K (chi - chi0)^2 on the bridging O (Si-Ob-Hb-Al)                 (improper class2, Wilson)
  nonbond     q_i q_j / r + A_i A_j / r^9 (the paper's B_i = 0); charges from the bond increments (half the ab initio charges);
              written as lj/class2 eps [2 (sigma/r)^9 - 3 (sigma/r)^6] with A_i A_j = 2 eps sigma^9 and B_ij = 3 eps sigma^6 =
              1e-10 kcal/mol A^6 (lj/class2 needs B > 0; < 1e-9 kcal/mol at 1 A)
The geometry-dependent charges of Table 9 (a Discover option) are not used; 1-2 and 1-3 pairs excluded, 1-4 in full (CFF).

usage: make_hill_sauer1995.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "J.-R. Hill, J. Sauer, J. Phys. Chem. 99, 9536 (1995)"
REF1 = "J.-R. Hill, J. Sauer, J. Phys. Chem. 98, 1238 (1994)"
FID = "hill-sauer1995"
B_TINY = 1e-10
# type: element, mass, description, A_i ((kcal/mol A^9)^1/2)
T = {"Si": ("Si", 28.0855, "Si", 432.3320), "Al": ("Al", 26.982, "Al", 143.1830),
     "Ob": ("O", 15.999, "bridging hydroxyl O, Si-O(H)-Al", 278.9970), "Hb": ("H", 1.008, "bridging hydroxyl H", 7.7660),
     "Oss": ("O", 15.999, "Si-O-Si O", 239.6090), "Oas": ("O", 15.999, "Si-O-Al O", 282.2060),
     "Osh": ("O", 15.999, "silanol O, Si-OH", 150.7170), "Hos": ("H", 1.008, "silanol H", 20.9470),
     "Oah": ("O", 15.999, "aluminol O, Al-OH", 168.1940), "Hoa": ("H", 1.008, "aluminol H", 8.4760)}
BONDS = [("Si", "Ob", 1.6581, 370.5744, -569.4988, 3102.2502), ("Si", "Oas", 1.6157, 494.1042, -36.7016, 2150.6836),
         ("Si", "Oss", 1.6104, 459.0786, -672.4445, 443.3651), ("Si", "Osh", 1.6238, 431.7825, -799.2202, 4446.0952),
         ("Al", "Ob", 1.9698, 86.8129, 26.3933, 758.4299), ("Al", "Oas", 1.7193, 328.6850, -341.0143, 2189.0349),
         ("Al", "Oah", 1.7096, 311.6053, -734.8390, 6022.3970), ("Ob", "Hb", 0.9540, 656.8465, -1627.4874, 3684.1929),
         ("Osh", "Hos", 0.9458, 707.2556, -1531.5564, 1603.1003), ("Oah", "Hoa", 0.9408, 730.0168, -1566.4841, 1703.5916)]
ANGLES = [("Ob", "Si", "Osh", 103.9653, 103.5467, -19.2943, 331.9679), ("Osh", "Si", "Osh", 111.3627, 58.2199, 96.8353, 527.5761),
          ("Ob", "Al", "Oas", 81.4988, 43.2829, -113.7790, 202.1476), ("Ob", "Al", "Oah", 97.9051, 27.6920, -41.0318, 35.4654),
          ("Oas", "Al", "Oah", 118.3370, 48.9048, 180.6775, 1557.5081), ("Oah", "Al", "Oah", 107.0368, 76.4087, -274.1490, 293.2546),
          ("Si", "Ob", "Al", 136.6645, 195.2950, 48.8575, 185.3481), ("Si", "Ob", "Hb", 119.2823, 43.0215, -15.2099, 13.1396),
          ("Al", "Ob", "Hb", 107.1200, 43.9725, -53.4364, 103.3838), ("Si", "Osh", "Hos", 118.9384, 20.0134, -33.2675, 42.4295),
          ("Al", "Oas", "Si", 162.4000, 10.3454, 12.6366, 9.0038), ("Al", "Oah", "Hoa", 131.3984, 7.8119, -1.6786, 53.6329),
          ("Oas", "Si", "Osh", 118.5582, 87.1068, 167.1975, 278.6751), ("Oas", "Si", "Oss", 112.4279, 88.0709, -57.0430, 92.4803),
          ("Oss", "Si", "Osh", 109.9424, 91.8675, 17.2882, 252.8231), ("Oss", "Si", "Oss", 112.0200, 81.9691, -36.5814, 116.9558),
          ("Oss", "Si", "Ob", 107.6066, 99.6715, -33.1512, 10.7204), ("Si", "Oss", "Si", 173.7651, 20.7015, 27.5506, 10.9930),
          ("Oas", "Al", "Oas", 113.4000, 300.0520, -32.5659, 67.2591), ("Oas", "Si", "Ob", 102.9438, 103.8270, 10.8990, 0.6674),
          ("Oas", "Si", "Oas", 110.6120, 154.1860, -68.6595, 23.6292)]
TORS = [("Osh", "Si", "Ob", "Al", 0.6228, 1.0917, -0.1696), ("Osh", "Si", "Ob", "Hb", 4.6176, -1.1913, -0.1141),
        ("Ob", "Si", "Osh", "Hos", 1.0464, -1.5833, -0.4804), ("Osh", "Si", "Osh", "Hos", 1.7740, -0.0772, 0.0326),
        ("Oas", "Al", "Ob", "Si", -4.1073, 3.5859, 0.2635), ("Oas", "Al", "Ob", "Hb", -3.6990, -3.0937, 0.8520),
        ("Oah", "Al", "Ob", "Si", -7.3022, 3.1717, 0.3184), ("Oah", "Al", "Ob", "Hb", -6.3428, -3.6734, -0.3907),
        ("Ob", "Al", "Oas", "Si", 2.1752, -0.9201, 1.1534), ("Oah", "Al", "Oas", "Si", 0.9847, -0.7285, -0.5777),
        ("Ob", "Al", "Oah", "Hoa", -0.9767, 0.7745, -0.9612), ("Oas", "Al", "Oah", "Hoa", -1.8010, 0.6647, 0.5364),
        ("Oah", "Al", "Oah", "Hoa", -1.1075, 0.8045, 0.3328), ("Al", "Oas", "Si", "Osh", 3.0650, 0.3404, 0.0146),
        ("Oas", "Si", "Osh", "Hos", 3.0272, -0.1573, -0.2321), ("Oss", "Si", "Oas", "Al", 2.1989, 0.7528, -0.4936),
        ("Oas", "Si", "Oss", "Si", 0.0397, 0.0133, -0.2453), ("Osh", "Si", "Oss", "Si", 0.1277, 0.0442, -0.0009),
        ("Oss", "Si", "Osh", "Hos", 3.2254, -0.2132, -0.4145), ("Oss", "Si", "Oss", "Si", 0.0306, -0.0105, 0.0804),
        ("Ob", "Si", "Oss", "Si", 0.5350, -0.4717, -0.0036), ("Oss", "Si", "Ob", "Al", 0.1641, 1.7914, 0.1970),
        ("Oss", "Si", "Ob", "Hb", 3.7091, -0.6560, 0.0856), ("Oas", "Al", "Oas", "Si", 6.0838, -0.5702, -1.4818),
        ("Ob", "Si", "Oas", "Al", 5.7178, 0.1922, -0.3374), ("Oas", "Si", "Ob", "Al", -1.8325, 2.2280, 0.2310),
        ("Oas", "Si", "Ob", "Hb", 5.7768, -2.1457, 0.1001), ("Oas", "Si", "Oas", "Al", 1.8597, -0.0231, -0.4345)]
# bond-bond F_bb' and bond-angle F_btheta, F_b'theta (None: the symmetric angle's own) per angle, from the page image
BB = {("Ob", "Si", "Osh"): 131.7627, ("Osh", "Si", "Osh"): 28.5745, ("Ob", "Al", "Oas"): 84.3372, ("Ob", "Al", "Oah"): 37.6888,
      ("Oas", "Al", "Oah"): 49.9526, ("Oah", "Al", "Oah"): 19.3587, ("Si", "Ob", "Al"): 130.8388, ("Si", "Ob", "Hb"): 18.0332,
      ("Al", "Ob", "Hb"): -18.4463, ("Si", "Osh", "Hos"): 9.6155, ("Al", "Oas", "Si"): 83.1658, ("Al", "Oah", "Hoa"): -2.6112,
      ("Oas", "Si", "Osh"): 72.4616, ("Oas", "Si", "Oss"): 45.7520, ("Oss", "Si", "Osh"): 70.6866, ("Oss", "Si", "Ob"): 63.3273,
      ("Si", "Oss", "Si"): 151.8742, ("Oas", "Al", "Oas"): -56.5050, ("Oas", "Si", "Ob"): 63.3273, ("Oas", "Si", "Oas"): 151.8742}
BA = {("Ob", "Si", "Osh"): (-0.9347, 96.7623), ("Osh", "Si", "Osh"): (43.7331, None), ("Ob", "Al", "Oas"): (12.9574, 20.0579),
      ("Ob", "Al", "Oah"): (25.6722, 20.7253), ("Oas", "Al", "Oah"): (87.4126, 15.6056), ("Oah", "Al", "Oah"): (-6.2892, None),
      ("Si", "Ob", "Al"): (23.0450, 3.6805), ("Si", "Ob", "Hb"): (36.5384, 53.8615), ("Al", "Ob", "Hb"): (-12.7238, 33.9457),
      ("Si", "Osh", "Hos"): (46.1485, 25.1227), ("Al", "Oas", "Si"): (9.1662, 13.7478), ("Al", "Oah", "Hoa"): (38.1658, 31.2114),
      ("Oas", "Si", "Osh"): (-54.0838, 63.3330), ("Oas", "Si", "Oss"): (20.1222, 77.9618), ("Oss", "Si", "Osh"): (41.1242, 75.6000),
      ("Oss", "Si", "Oss"): (78.1239, None), ("Oss", "Si", "Ob"): (52.7776, 12.2769), ("Si", "Oss", "Si"): (9.2390, None),
      ("Oas", "Al", "Oas"): (110.1523, None), ("Oas", "Si", "Ob"): (13.3600, -8.3197), ("Oas", "Si", "Oas"): (234.0716, None)}
# angle-angle i-j-k/l: (i, j, k, l)
AA = [("Osh", "Si", "Oss", "Osh", 52.6194), ("Oss", "Si", "Osh", "Osh", 60.1802), ("Osh", "Si", "Osh", "Osh", 38.4444),
      ("Oss", "Si", "Oss", "Osh", 52.3408), ("Oss", "Si", "Osh", "Oss", 86.3954), ("Oss", "Si", "Oss", "Oss", -6.3030),
      ("Osh", "Si", "Ob", "Osh", -13.6331), ("Ob", "Si", "Osh", "Osh", 17.4074), ("Oah", "Al", "Ob", "Oah", 0.0000),
      ("Ob", "Al", "Oah", "Oah", -20.0594), ("Oah", "Al", "Oah", "Oah", -51.8916), ("Al", "Ob", "Si", "Hb", 36.0652),
      ("Si", "Ob", "Al", "Hb", 58.8895), ("Si", "Ob", "Hb", "Al", 44.4559), ("Oss", "Si", "Ob", "Osh", 61.4849),
      ("Ob", "Si", "Oss", "Osh", 75.8868), ("Ob", "Si", "Osh", "Oss", 50.5839), ("Oas", "Al", "Ob", "Oah", 13.6760),
      ("Ob", "Al", "Oas", "Oah", 28.2428), ("Ob", "Al", "Oah", "Oas", -34.0834), ("Oah", "Al", "Oas", "Oah", -62.0442),
      ("Oas", "Al", "Oah", "Oah", -19.2127), ("Osh", "Si", "Oas", "Osh", 13.9083), ("Oas", "Si", "Osh", "Osh", 46.1582),
      ("Oas", "Si", "Oss", "Osh", 51.3039), ("Oss", "Si", "Oas", "Osh", 40.7833), ("Oss", "Si", "Osh", "Oas", 52.5469),
      ("Oas", "Al", "Ob", "Oas", -6.1120), ("Ob", "Al", "Oas", "Oas", 39.4951), ("Oas", "Al", "Oas", "Oah", 0.0000),
      ("Oas", "Al", "Oah", "Oas", 98.5159), ("Ob", "Si", "Oas", "Osh", 25.6772), ("Oas", "Si", "Ob", "Osh", -35.6378),
      ("Oas", "Si", "Osh", "Ob", 31.8901), ("Oas", "Si", "Oss", "Oss", 19.2372), ("Oss", "Si", "Oas", "Oss", 0.0000),
      ("Oss", "Si", "Ob", "Oss", 79.7785), ("Ob", "Si", "Oss", "Oss", 0.0000), ("Oas", "Si", "Osh", "Oas", 0.0000),
      ("Osh", "Si", "Oas", "Oas", -23.1336), ("Oas", "Si", "Ob", "Oas", 0.0000), ("Ob", "Si", "Oas", "Oas", 0.0000)]
AAT = [("Osh", "Si", "Ob", "Al", 213.6799), ("Osh", "Si", "Ob", "Hb", -21.4028), ("Ob", "Si", "Osh", "Hos", -9.3392),
       ("Osh", "Si", "Osh", "Hos", -2.9461), ("Oas", "Al", "Ob", "Si", 32.3177), ("Oas", "Al", "Ob", "Hb", 8.8676),
       ("Oah", "Al", "Ob", "Si", -5.9667), ("Oah", "Al", "Ob", "Hb", -22.3131), ("Ob", "Al", "Oas", "Si", 16.7797),
       ("Oah", "Al", "Oas", "Si", 7.8776), ("Ob", "Al", "Oah", "Hoa", -6.3680), ("Oas", "Al", "Oah", "Hoa", -52.1956),
       ("Oah", "Al", "Oah", "Hoa", 6.5018), ("Al", "Oas", "Si", "Osh", -14.1105), ("Oas", "Si", "Osh", "Hos", -4.6863),
       ("Oss", "Si", "Oas", "Al", 8.4655), ("Oas", "Si", "Oss", "Si", -10.4191), ("Osh", "Si", "Oss", "Si", 2.3216),
       ("Oss", "Si", "Osh", "Hos", 6.4184), ("Oss", "Si", "Oss", "Si", -4.5150), ("Ob", "Si", "Oss", "Si", 11.8037),
       ("Oss", "Si", "Ob", "Al", -71.4139), ("Oss", "Si", "Ob", "Hb", 15.1160), ("Oas", "Al", "Oas", "Si", -18.0184),
       ("Ob", "Si", "Oas", "Al", 17.2526), ("Oas", "Si", "Ob", "Al", 91.5202), ("Oas", "Si", "Ob", "Hb", 17.5674),
       ("Oas", "Si", "Oas", "Al", -10.9279)]
# bond increments delta(i) for the bond i-j (j gets -delta)
BI = [("Si", "Osh", 0.1303), ("Hos", "Osh", 0.0641), ("Al", "Oah", 0.3013), ("Hoa", "Oah", 0.0575), ("Si", "Ob", 0.1392),
      ("Al", "Ob", 0.0284), ("Hb", "Ob", 0.0839), ("Si", "Oss", 0.1309), ("Si", "Oas", 0.1265), ("Al", "Oas", 0.1694)]
names = list(T)
types = [{"name": n, "element": e, "mass": m, "description": f"{d} (Table 8)", "equivalence": {}} for n, (e, m, d, A) in T.items()]
pairs = []
for i, a in enumerate(names):
    for b in names[i:]:
        Aij = T[a][3] * T[b][3]
        sigma = (1.5 * Aij / B_TINY) ** (1 / 3)
        eps = B_TINY / (3 * sigma ** 6)
        pairs.append({"name": a if a == b else f"{a}-{b}", "match": [a] if a == b else [a, b], "params": [eps, sigma],
                      "comment": f"A_i A_j = {T[a][3]} x {T[b][3]} kcal/mol A^9 (Table 8), B 0 (written {B_TINY:g})"})
bonds = [{"name": f"{a}-{b}", "match": [a, b], "style": "class2", "params": [b0, k2, k3, k4], "comment": "Table 8"} for a, b, b0, k2, k3, k4 in BONDS]
angles = [{"name": f"{a}-{b}-{c}", "match": [a, b, c], "style": "class2", "params": [t0, h2, h3, h4], "comment": "Table 8"} for a, b, c, t0, h2, h3, h4 in ANGLES]
dihedrals = [{"name": f"{a}-{b}-{c}-{d}", "match": [a, b, c, d], "style": "class2", "params": [v1, 0.0, v2, 0.0, v3, 0.0], "comment": "Table 8"}
             for a, b, c, d, v1, v2, v3 in TORS]
impropers = [{"name": "Si-Ob-Hb-Al", "match": ["Si", "Ob", "Hb", "Al"], "style": "class2", "params": [0.8407, -0.2235], "comment": "out-of-plane (Table 8)"}]
cross = {"bb": [{"name": "-".join(k), "match": list(k), "params": [v], "comment": "Table 8"} for k, v in BB.items()],
         "ba": [{"name": "-".join(k), "match": list(k), "params": [v[0], v[0] if v[1] is None else v[1]], "comment": "Table 8"} for k, v in BA.items()],
         "aa": [{"name": "-".join(x[:4]), "match": list(x[:4]), "params": [x[4]], "comment": "Table 8 (i-j-k/l)"} for x in AA],
         "aat": [{"name": "-".join(x[:4]), "match": list(x[:4]), "params": [x[4]], "comment": "Table 8"} for x in AAT],
         "mbt": [], "ebt": [], "at": [], "bb13": []}
bi = [{"name": f"{a}-{b}", "match": [a, b], "params": [d, -d], "comment": "Table 8"} for a, b, d in BI]
SI, AL = "[Si]", "[Al]"
rules = [{"type": "Si", "smarts": "[Si]", "priority": 0, "description": "Si"},
         {"type": "Al", "smarts": "[Al]", "priority": 0, "description": "Al"},
         {"type": "Oss", "smarts": "[O;X2;$(O([Si])[Si])]", "priority": 1, "description": "Si-O-Si"},
         {"type": "Oas", "smarts": "[O;X2;$(O([Si])[Al])]", "priority": 1, "description": "Si-O-Al"},
         {"type": "Ob", "smarts": "[O;X3;$(O([Si])([Al])[#1])]", "priority": 2, "description": "bridging hydroxyl O"},
         {"type": "Osh", "smarts": "[O;X2;$(O([Si])[#1])]", "priority": 1, "description": "silanol O"},
         {"type": "Oah", "smarts": "[O;X2;$(O([Al])[#1])]", "priority": 1, "description": "aluminol O"},
         {"type": "Hb", "smarts": "[#1;$([#1]O([Si])[Al])]", "priority": 2, "description": "bridging hydroxyl H"},
         {"type": "Hos", "smarts": "[#1;$([#1][O;X2][Si])]", "priority": 1, "description": "silanol H"},
         {"type": "Hoa", "smarts": "[#1;$([#1][O;X2][Al])]", "priority": 1, "description": "aluminol H"}]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "Aluminosilicates (Hill & Sauer 1995)", "version": "1995 (Table 8)",
       "source": REF + ", Table 8", "references": [REF, REF1], "units": "real", "mixing": "sixthpower", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1],
       "cutoff": 10, "missing_bond_increment": "zero", "equivalence": "replace",
       "styles": {"pair": "lj/class2/coul/long", "bond": "class2", "angle": "class2", "dihedral": "class2", "improper": "class2"},
       "oop_scheme": "msi2lmp", "improper_order": "center2_sorted", "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles,
       "dihedrals": dihedrals, "impropers": impropers, "cross_rules": cross, "bond_increments": bi,
       "notes": ["built by CAPS (bench/ff/make_hill_sauer1995.py) from the paper's Table 8 (cross terms read from the page image)",
                 "nonbond A_i A_j / r^9 with B 0 as the paper, as lj/class2 with B 1e-10 kcal/mol A^6 (every pair explicit)",
                 "charges from the bond increments; the geometry-dependent charges of Table 9 are not used",
                 "cut-off 10 A (CAPS's; the paper's lattice sums are Ewald)"],
       "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped",
           "description": "Hill and Sauer's aluminosilicate O and H types by their neighbours (bench/ff/make_hill_sauer1995.py)", "rules": rules},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} torsions, "
      f"{sum(len(v) for v in cross.values())} cross terms")
