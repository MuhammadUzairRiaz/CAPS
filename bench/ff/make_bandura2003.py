#!/usr/bin/env python3
"""TiO2 surfaces and water for CAPS from A. V. Bandura, J. D. Kubicki, J. Phys. Chem. B 107, 11072 (2003), Tables 2, 3, 5,
6, 9 and 10 (kcal/mol, A, e).

  bulk     Matsui & Akaogi's rigid ions as Table 2 prints them: Buckingham A exp(-r/rho) - C/r^6 for Ti-O, Ti-Ti, O-O,
           q(Ti) 2.196, q(O) -1.098 (footnote a)
  surface  Table 9's charges: 5-fold Ti (TiV), 6-fold surface Ti (TiVI), Ti holding a terminal hydroxyl (TiVIoh) 2.156;
           bridging O (Ob) -1.018; 3-fold O in (Os) and under (Ou) the surface Ti layer -1.098; bridging hydroxyl O (Obh)
           -0.964 and H (Hbh) 0.455; terminal hydroxyl O (Ooh) -0.939 and H (Hoh) 0.430
           Table 10: TiVIoh-Ooh 315480.8 (Table 2's Ti-O^b), Ti(any)-Os 449223.3, O(any)-Os 570799.0, Os-Os 150150.9 (rho, C
           Table 2's); every other Ti-O, O-O and Ti-Ti pair Table 2's (Table 10 footnote a)
  hydroxyl (Ti)O-H Morse A {1 - exp[-alpha (r - r0)]}^2 - A (Table 5, from Fleming et al.), Ti-O-H 1/2 k2 (theta - theta0)^2
           (Table 6); no Coulomb between a hydroxyl's O and H; Ti-H Coulomb kept (section 3)
  water    SPC/E's charges and O-O term (Table 3, E = eps [(sigma/r)^12 - 2 (sigma/r)^6]: sigma is r_min), made flexible:
           (H)O-H Morse and H-O-H bend (Tables 5, 6); no Coulomb inside a molecule; Ti-O(H2) Buckingham (Table 2, row c);
           O(H2)-O(Ti) as SPC/E's O-O (section 4); H no van der Waals term
Ti-O bonds carry no term (all Ti-O interactions are pairs): typing reads them to count each Ti's and O's neighbours, then
they are dropped; the Ti-O-H bends take the Ti within 2.6 A of the hydroxyl O (CAPS's distance). The Morse constant -A is
left out (it shifts no force). Van der Waals cut-off 9.6 A (section 1).

Typing (CAPS's reading of Table 9's site names, by neighbours): TiV a Ti with five oxide O (or six counting an adsorbed
water O); TiVIoh a Ti holding a terminal hydroxyl; TiVI a six-fold Ti holding a bridging O (Ob or Obh); every other Ti
bulk; Ob an O between two Ti without H; Obh the same with H; Ooh an O on one Ti with H; Os a three-fold O bonded to a TiV
or TiVIoh (the in-plane O of the (110) surface layer); every other three-fold O bulk (Ou has the bulk O's charge and
terms); water O with two H.

usage: make_bandura2003.py [DATA_DIR] [--test-slab OUT.pdb]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "A. V. Bandura, J. D. Kubicki, J. Phys. Chem. B 107, 11072 (2003)"
FID = "bandura2003-tio2-water"
MASS = {"Ti": 47.867, "O": 15.999, "H": 1.008}
# Table 9 (and Table 2 footnote a, Table 3 footnote a): type, element, charge, description
TYPES = [("Ti_bk", "Ti", 2.196, "Ti, bulk (Table 9, Matsui & Akaogi's value)"),
         ("TiV_bk", "Ti", 2.156, "5-fold surface Ti (Table 9 TiV)"),
         ("TiVI_bk", "Ti", 2.156, "6-fold surface Ti (Table 9 TiVI)"),
         ("TiVIoh_bk", "Ti", 2.156, "6-fold Ti holding a terminal hydroxyl (Table 9 TiVIoh)"),
         ("O_bk", "O", -1.098, "O, bulk; also Ou, the 3-fold O under the surface Ti layer (Table 9)"),
         ("Os_bk", "O", -1.098, "3-fold O in the surface Ti layer (Table 9 Os)"),
         ("Ob_bk", "O", -1.018, "bridging O (Table 9 Ob)"),
         ("Obh_bk", "O", -0.964, "bridging hydroxyl O (Table 9 Obh)"),
         ("Ooh_bk", "O", -0.939, "terminal hydroxyl O (Table 9 Ooh)"),
         ("Hbh_bk", "H", 0.455, "bridging hydroxyl H (Table 9 Hbh)"),
         ("Hoh_bk", "H", 0.430, "terminal hydroxyl H (Table 9 Hoh)"),
         ("Ow_bk", "O", -0.8476, "water O (SPC/E, Table 3 footnote a)"),
         ("Hw_bk", "H", 0.4238, "water H (SPC/E, Table 3 footnote a)")]
TI = ["Ti_bk", "TiV_bk", "TiVI_bk", "TiVIoh_bk"]
OX = ["O_bk", "Os_bk", "Ob_bk", "Obh_bk", "Ooh_bk"]   # the oxide's O, every one of them "O" in Tables 2 and 10
# Table 2: A (kcal/mol), rho (A), C (kcal/mol A^6)
TIO, TIOB, TITI, OO, TIW = (391049.1, 0.194, 290.3317), (315480.8, 0.194, 290.3317), (717647.4, 0.154, 121.0676), (271716.3, 0.234, 696.8883), (28593.02, 0.265, 148.000)
# Table 10
TIC_OS, OD_OS, OS_OS = 449223.3, 570799.0, 150150.9
# Table 3: SPC/E O-O, eps (kcal/mol), r_min (A)
EPS_W, RMIN_W = 0.15539, 3.5532
types = [{"name": n, "element": el, "mass": MASS[el], "charge": q, "description": d} for n, el, q, d in TYPES]
pairs = [{"name": n, "match": [n], "params": [0.0, 0.0]} for n, *_ in TYPES]


def buck(a, b, A, rho, C, why):
    return {"name": f"{a}-{b}", "match": [a, b], "style": "buck", "params": [A, rho, C],
            "comment": f"A {A} kcal/mol, rho {rho} A, C {C} kcal/mol A^6 ({why})"}


for i, a in enumerate(TI):
    for b in TI[i:]:
        pairs.append(buck(a, b, *TITI, "Table 2 Ti-Ti"))
    for o in OX:
        if o == "Os_bk":
            pairs.append(buck(a, o, TIC_OS, TIO[1], TIO[2], "Table 10 Ti^c-Os, rho and C Table 2's"))
        elif a == "TiVIoh_bk" and o == "Ooh_bk":
            pairs.append(buck(a, o, *TIOB, "Table 10 TiVIoh-Ooh = Table 2 Ti-O^b"))
        else:
            pairs.append(buck(a, o, *TIO, "Table 2 Ti-O"))
    pairs.append(buck(a, "Ow_bk", *TIW, "Table 2 Ti-O(H2)"))
for i, a in enumerate(OX):
    for b in OX[i:]:
        if a == b == "Os_bk":
            pairs.append(buck(a, b, OS_OS, OO[1], OO[2], "Table 10 Os-Os, rho and C Table 2's"))
        elif "Os_bk" in (a, b):
            pairs.append(buck(a, b, OD_OS, OO[1], OO[2], "Table 10 O^d-Os, rho and C Table 2's"))
        else:
            pairs.append(buck(a, b, *OO, "Table 2 O-O"))
sig = RMIN_W / 2 ** (1 / 6)
for o in OX + ["Ow_bk"]:
    pairs.append({"name": f"{o}-Ow_bk", "match": [o, "Ow_bk"], "params": [EPS_W, sig],
                  "comment": f"SPC/E O-O: eps {EPS_W} kcal/mol, r_min {RMIN_W} A (sigma {sig:.6f} A) (Table 3; O(H2)-O(Ti) as H2O-OH2, section 4)"})
bonds = [{"name": "Ow_bk-Hw_bk", "match": ["Ow_bk", "Hw_bk"], "style": "morse", "params": [101.905, 2.3470, 1.00],
          "comment": "(H)O-H Morse A 101.905 kcal/mol, alpha 2.3470 1/A, r0 1.00 A (Table 5; the constant -A left out)"}]
bonds += [{"name": f"{o}-{h}", "match": [o, h], "style": "morse", "params": [125.094, 2.2682, 0.95],
           "comment": "(Ti)O-H Morse A 125.094 kcal/mol, alpha 2.2682 1/A, r0 0.95 A (Table 5, Fleming et al.'s; the constant -A left out)"}
          for o, h in (("Obh_bk", "Hbh_bk"), ("Ooh_bk", "Hoh_bk"))]
angles = [{"name": "Hw_bk-Ow_bk-Hw_bk", "match": ["Hw_bk", "Ow_bk", "Hw_bk"], "style": "harmonic", "params": [103.045 / 2, 109.47],
           "comment": "H-O-H 1/2 k2 (theta - theta0)^2, k2 103.045 kcal/(mol rad^2), theta0 109.47 (Table 6)"}]
angles += [{"name": f"{t}-{o}-{h}", "match": [t, o, h], "style": "harmonic", "params": [14.136 / 2, 90.85],
            "comment": "Ti-O-H 1/2 k2 (theta - theta0)^2, k2 14.136 kcal/(mol rad^2), theta0 90.85 (Table 6)"}
           for t in TI for o, h in (("Obh_bk", "Hbh_bk"), ("Ooh_bk", "Hoh_bk"))]
# a Ti's oxide neighbours counted by recursive SMARTS: bond perception also links Ti to Ti (2.96 A apart), so D and X
# count those too; a water O on a Ti is not counted (TiV keeps its name with water adsorbed)
OXO = "~[O;!H2]"
OX5 = "$(*(" + ")(".join([OXO] * 5) + "))"
OX6 = "$(*(" + ")(".join([OXO] * 6) + "))"
TIV = f"[Ti;{OX5};!{OX6}]"
TIOH = "[Ti;$(*~[O;X2;H1])]"
TIVI = "[Ti;$(*~[O;X2;H0]),$(*~[O;X3;H1])]"
rules = [
    {"type": "Ti_bk", "smarts": "[Ti]", "priority": 0, "description": "Ti, bulk"},
    {"type": "TiVI_bk", "smarts": TIVI, "priority": 1, "description": "6-fold surface Ti: holds a bridging O (Ob or Obh)"},
    {"type": "TiV_bk", "smarts": TIV, "priority": 2, "description": "Ti with five oxide O (TiV; a water O on its open site not counted)"},
    {"type": "TiVIoh_bk", "smarts": TIOH, "priority": 3, "description": "Ti holding a terminal hydroxyl (TiVIoh)"},
    {"type": "O_bk", "smarts": "[O;H0]", "priority": 0, "requires": ["Ti"], "description": "O, bulk (and Ou)"},
    {"type": "Os_bk", "smarts": f"[O;X3;H0;$(*~{TIV}),$(*~{TIOH});$(*~{TIVI})]", "priority": 2, "requires": ["Ti"],
     "description": "3-fold O between a TiV (or TiVIoh) and a TiVI: Os, the in-plane O of the surface Ti layer"},
    {"type": "Ob_bk", "smarts": "[O;X2;H0;$(*(~[Ti])~[Ti])]", "priority": 3, "description": "bridging O between two Ti (Ob)"},
    {"type": "Obh_bk", "smarts": "[O;X3;H1;$(*(~[Ti])~[Ti])]", "priority": 3, "description": "bridging hydroxyl O (Obh)"},
    {"type": "Ooh_bk", "smarts": "[O;X2;H1;$(*~[Ti])]", "priority": 3, "description": "terminal hydroxyl O (Ooh)"},
    {"type": "Hbh_bk", "smarts": "[H][O;X3;H1;$(*(~[Ti])~[Ti])]", "priority": 3, "description": "bridging hydroxyl H (Hbh)"},
    {"type": "Hoh_bk", "smarts": "[H][O;X2;H1;$(*~[Ti])]", "priority": 3, "description": "terminal hydroxyl H (Hoh)"},
    {"type": "Ow_bk", "smarts": "[O;H2]", "priority": 5, "description": "water O"},
    {"type": "Hw_bk", "smarts": "[H][O;H2]", "priority": 5, "description": "water H"},
]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "TiO2-H2O (Bandura & Kubicki 2003)",
       "source": REF + ", Tables 2, 3, 5, 6, 9, 10", "references": [REF], "units": "real",
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 9.6, "angle_contacts": 2.6,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_bandura2003.py) from the paper's Tables 2, 3, 5, 6, 9 and 10",
                 "Ti and O are rigid ions; every Ti-O, Ti-Ti and O-O interaction a Buckingham pair (Table 2, Table 10 for the surface)",
                 "Table 10's Ti^c (any Ti) and O^d (any oxide O) read as every Ti and every oxide O type; water O takes SPC/E's O-O with every O",
                 "surface sites by neighbours (CAPS's reading of Table 9's names): TiV five-fold, TiVI holds a bridging O, TiVIoh a terminal OH, Os a 3-fold O on a TiV / TiVIoh",
                 "Morse O-H without the constant -A of Table 5's expression; bends 1/2 k2 (theta - theta0)^2 written as K = k2/2",
                 "Ti-O-H bends on the Ti within 2.6 A of the hydroxyl O (no Ti-O bond is made; the distance is CAPS's)",
                 "1-2 and 1-3 pairs excluded: no Coulomb inside a hydroxyl or a water molecule; Ti-H across Ti-O-H kept (section 3)",
                 "van der Waals cut-off 9.6 A (section 1)"],
       "version": "2003", "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined",
           "description": "rutile surfaces, hydroxyls and water by neighbours (bench/ff/make_bandura2003.py); Ti-O bonds dropped after typing",
           "rules": rules}, open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1)
print(f"{FID}: {len(types)} types, {len(pairs)} pair entries, {len(bonds)} bonds, {len(angles)} bends")


def test_slab(path):
    """A rutile (110) slab for the tests (tests/data/rutile110_hydroxylated.pdb): three Ti layers, 2 x 3 surface cells
    (12.99 x 8.88 A) in a 40 A cell; on top one bridging O protonated, one TiV holding a terminal OH (the dissociated
    water of the paper's half-monolayer model) and one water on another TiV (Ti-O 2.2 A); the bottom left clean."""
    import math
    a, c, u = 4.594, 2.959, 0.3053   # rutile, the cell of data/crystals/rutile.cif
    s2 = math.sqrt(2)
    A, B, d = a * s2, c, a / s2      # [1-10] and [001] periods; (110) Ti-plane spacing
    base = [("Ti", 0, 0, 0), ("Ti", .5, .5, .5), ("O", u, u, 0), ("O", -u, -u, 0), ("O", .5 + u, .5 - u, .5), ("O", .5 - u, .5 + u, .5)]
    hb = d - 2 * u * d               # height of a bridging O above its Ti plane, (1 - 2u) a / sqrt(2)
    lo, hi = -hb - 0.05, 2 * d + hb + 0.05
    cell = {}
    for i in range(-6, 7):
        for j in range(-6, 7):
            for k in range(-2, 3):
                for el, x, y, z in base:
                    X, Y, Z = (x + i) * a, (y + j) * a, (z + k) * c
                    xp, yp, zp = (X - Y) / s2, Z, (X + Y) / s2
                    if lo <= zp <= hi:
                        fx, fy = xp % A, yp % B
                        fx = 0.0 if A - fx < 1e-6 else fx
                        fy = 0.0 if B - fy < 1e-6 else fy
                        cell[(el, round(fx, 2), round(fy, 2), round(zp, 2))] = (el, fx, fy, zp)
    atoms = [(el, x + p * A, y + q * B, z + 7.5) for p in range(2) for q in range(3) for el, x, y, z in cell.values()]
    top = 2 * d + 7.5
    ti_top = [t for t in atoms if t[0] == "Ti" and abs(t[3] - top) < 0.1]
    bridge = [t for t in atoms if t[0] == "O" and abs(t[3] - (top + hb)) < 0.1]
    # a TiV has no bridging O above its row: the Ti whose x differs from every bridging O's
    tiv = sorted([t for t in ti_top if all(abs(t[1] - o[1]) > 0.5 for o in bridge)], key=lambda t: (t[1], t[2]))
    ob = sorted(bridge, key=lambda t: (t[1], t[2]))[0]
    add = [("H", ob[1], ob[2], ob[3] + 0.95)]                                            # Obh-Hbh
    t1, t2 = tiv[0], tiv[-1]
    add += [("O", t1[1], t1[2], t1[3] + 1.85), ("H", t1[1] + 0.90, t1[2], t1[3] + 1.85 + 0.30)]   # Ooh-Hoh on a TiV
    add += [("O", t2[1], t2[2], t2[3] + 2.20), ("H", t2[1] + 0.757, t2[2], t2[3] + 2.20 + 0.586),
            ("H", t2[1] - 0.757, t2[2], t2[3] + 2.20 + 0.586)]                         # water on a TiV
    atoms += add
    with open(path, "w") as f:
        f.write("REMARK   rutile (110), 3 Ti layers, hydroxylated and with one water on top: bench/ff/make_bandura2003.py --test-slab\n")
        f.write("CRYST1%9.3f%9.3f%9.3f  90.00  90.00  90.00 P 1           1\n" % (2 * A, 3 * B, 40.0))
        for n, (el, x, y, z) in enumerate(atoms, 1):
            f.write("HETATM%5d %-4s MOL A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (n, el + str(n % 100), x, y, z, el))
        f.write("END\n")
    print(f"wrote {path}: {len(atoms)} atoms")


if "--test-slab" in sys.argv:
    test_slab(sys.argv[sys.argv.index("--test-slab") + 1])
