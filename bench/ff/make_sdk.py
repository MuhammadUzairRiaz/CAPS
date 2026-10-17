#!/usr/bin/env python3
"""SDK (Shinoda-DeVane-Klein, the model SPICA grew from) for CAPS, from the papers' tables (PDFs in
~/project/ff-sources/papers/cg, never in the repository):

  W. Shinoda, R. DeVane, M. L. Klein, Mol. Simul. 33, 27 (2007): Table 1 (segments W, CT, CM, CT2, EOT, EO, OA), Table 4
  (bonds, k (r - r0)^2 in kcal/mol/A^2), Table 5 (angles, k (theta - theta0)^2 in kcal/mol/rad^2), Table 6 (non-bonded:
  LJ12-4 for pairs with W, LJ9-6 otherwise, epsilon kcal/mol, sigma A, cut-off 15 A)
  W. Shinoda, R. DeVane, M. L. Klein, J. Phys. Chem. B 114, 6836 (2010): Table 2 (lipid segments NC, NH, PH, PHE, GL,
  EST1, EST2, CMD2; NC and NH +1, PH and PHE -1), Tables 3-5 (bonds, angles, non-bonded), eq 3 (the 1-3 correction:
  the pair's potential cut and shifted at its minimum, LAMMPS angle sdk), the Coulomb sum without truncation (Ewald)

Pairs separated by more than two bonds interact through the non-bonded terms (special 0 0 1). Every pair is given (no
mixing rule): a pair of types no table lists is reported missing. Masses are the segments' formulas (Table 1 / Table 2);
the papers screen Coulomb by a relative permittivity whose value they do not print: 80, water's, is CAPS's choice. The
cholesterol and cyclic-alkane beads of MacDermaid et al. 2015 are not included (their dihedral terms and the cholesterol
tables of that paper's Supporting Information are left for later).

usage: make_sdk.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
FID = "sdk-shinoda2007"
S2007 = "W. Shinoda, R. DeVane, M. L. Klein, Mol. Simul. 33, 27 (2007)"
S2010 = "W. Shinoda, R. DeVane, M. L. Klein, J. Phys. Chem. B 114, 6836 (2010)"
M = {"C": 12.011, "H": 1.008, "O": 15.9994, "N": 14.007, "P": 30.973762}
# segment: (formula counts, charge, description, paper)
SEG = {
    "W": ({"H": 6, "O": 3}, 0, "(H2O)3", "2007"), "CT": ({"C": 3, "H": 7}, 0, "CH3CH2CH2-", "2007"),
    "CM": ({"C": 3, "H": 6}, 0, "-CH2CH2CH2-", "2007"), "CT2": ({"C": 2, "H": 5}, 0, "CH3CH2-", "2007"),
    "EOT": ({"C": 2, "H": 5, "O": 1}, 0, "CH3-O-CH2-", "2007"), "EO": ({"C": 2, "H": 4, "O": 1}, 0, "-CH2-O-CH2-", "2007"),
    "OA": ({"C": 1, "H": 3, "O": 1}, 0, "HOCH2-", "2007"),
    "NC": ({"C": 5, "H": 13, "N": 1}, 1, "-CH2CH2-N-(CH3)3 (+1)", "2010"), "NH": ({"C": 2, "H": 7, "N": 1}, 1, "-CH2CH2-NH3 (+1)", "2010"),
    "PH": ({"P": 1, "O": 4}, -1, "-PO4- (-1)", "2010"), "PHE": ({"P": 1, "O": 4}, -1, "-PO4- (-1), PE head group", "2010"),
    "GL": ({"C": 3, "H": 5}, 0, "-CH2CH-CH2-", "2010"), "EST1": ({"C": 2, "H": 2, "O": 2}, 0, "-CH2CO2- (sn-2 chain)", "2010"),
    "EST2": ({"C": 2, "H": 2, "O": 2}, 0, "-CH2CO2- (sn-1 chain)", "2010"), "CMD2": ({"C": 2, "H": 2}, 0, "-HC=CH- (cis)", "2010"),
}
# 2007 Table 6 and 2010 Table 5: a b function epsilon sigma
PAIRS_2007 = """W W 12-4 0.8950 4.3710
CT CT 9-6 0.4690 4.5850
CT CM 9-6 0.4440 4.5455
CM CM 9-6 0.4200 4.5060
W CT 12-4 0.3600 4.4780
W CM 12-4 0.3400 4.4385
CT2 CT2 9-6 0.3120 4.2210
CT2 CM 9-6 0.3620 4.3635
W CT2 12-4 0.2900 4.2960
OA OA 9-6 0.4491 3.7130
EOT EOT 9-6 0.4370 4.2500
EO EO 9-6 0.4050 4.2500
EO EOT 9-6 0.4200 4.2500
EO OA 9-6 0.4400 3.8900
EO CT 9-6 0.4100 4.3400
EO CM 9-6 0.3770 4.2740
EO CT2 9-6 0.3700 4.1400
W OA 12-4 0.7000 3.9500
W EO 12-4 0.5700 4.3100
CT OA 9-6 0.4372 4.0330
CM OA 9-6 0.3650 3.9870
CT2 OA 9-6 0.3800 3.8400"""
PAIRS_2010 = """PH CM 9-6 0.300 4.9530
PH CT 9-6 0.320 4.9925
NC CM 9-6 0.400 5.1280
NC CT 9-6 0.420 5.1675
PH W 12-4 1.000 4.0300
NC W 12-4 0.900 4.6100
PH PH 9-6 1.400 5.4000
NC NC 9-6 0.700 5.7500
NC PH 9-6 1.150 4.2000
GL GL 9-6 0.420 4.5060
GL EST1 9-6 0.470 4.4030
GL EST2 9-6 0.470 4.4030
GL CM 9-6 0.420 4.5060
GL CT 9-6 0.444 4.5455
GL W 12-4 0.640 4.4385
GL CT2 9-6 0.362 4.3635
PH CT2 9-6 0.280 4.8105
NC CT2 9-6 0.320 4.9855
CMD2 GL 9-6 0.312 4.2555
CMD2 EST1 9-6 0.440 4.0050
CMD2 EST2 9-6 0.440 4.0050
PH CMD2 9-6 0.300 4.7025
NC CMD2 9-6 0.350 4.8775
GL NC 9-6 0.650 4.6200
GL PH 9-6 0.300 4.7500
EST1 NC 9-6 0.750 4.4750
EST2 NC 9-6 0.750 4.4750
EST1 PH 9-6 0.500 4.5500
EST2 PH 9-6 0.500 4.5500
PHE PHE 9-6 1.400 4.6000
PHE CM 9-6 0.300 4.9530
PHE CMD2 9-6 0.300 4.7025
PHE CT 9-6 0.320 4.9925
PHE CT2 9-6 0.280 4.8105
PHE GL 9-6 0.300 4.7500
PHE EST1 9-6 0.500 4.5500
PHE EST2 9-6 0.500 4.5500
PHE W 12-4 1.000 4.0300
NH CM 9-6 0.330 4.5530
NH CMD2 9-6 0.300 4.3025
NH CT2 9-6 0.320 4.4105
NH CT 9-6 0.340 4.5925
NH NH 9-6 1.100 4.6000
NH PHE 9-6 1.200 3.8000
NH GL 9-6 0.750 4.1900
NH EST1 9-6 0.850 4.1100
NH EST2 9-6 0.850 4.1100
NH W 12-4 0.800 3.9500
PH PHE 9-6 1.400 5.0000
NH PH 9-6 1.200 3.8000
NC PHE 9-6 1.150 4.2000
NC NH 9-6 0.880 5.1750
CT CT 9-6 0.469 4.5850
CT CM 9-6 0.444 4.5455
CM CM 9-6 0.420 4.5060
CMD2 CMD2 9-6 0.232 4.0050
CMD2 CT 9-6 0.330 4.2950
CMD2 CT2 9-6 0.269 4.1130
CMD2 CM 9-6 0.312 4.2555
CMD2 W 12-4 0.270 4.1880
EST1 CM 9-6 0.470 4.4030
EST2 CM 9-6 0.470 4.4030
EST1 CT 9-6 0.470 4.4425
EST2 CT 9-6 0.470 4.4425
EST1 CT2 9-6 0.390 4.2605
EST2 CT2 9-6 0.390 4.2605
EST1 EST1 9-6 0.495 4.3000
EST1 EST2 9-6 0.495 4.3000
EST2 EST2 9-6 0.495 4.3000
EST1 W 12-4 0.820 4.2900
EST2 W 12-4 0.820 4.2900
W CT 12-4 0.360 4.4780
W CM 12-4 0.340 4.4385
CT2 CT2 9-6 0.312 4.2210
CT2 CM 9-6 0.362 4.3635
W CT2 12-4 0.290 4.2960
CT2 CT 9-6 0.383 4.4030
W W 12-4 0.895 4.3710"""
BONDS_2007 = """CT CT 6.955 3.71
CM CM 6.160 3.64
CM CT 6.160 3.65
CT2 CM 9.000 3.13
OA OA 63.000 2.16
EOT EOT 5.500 3.33
EOT EO 5.400 3.34
EO EO 4.900 3.28
EO OA 15.000 2.79
CT EO 7.100 3.61
CM EO 7.100 3.56
CT2 EO 10.000 3.07
CT OA 6.955 3.05
CM OA 7.500 3.01
CT2 OA 14.000 2.53"""
BONDS_2010 = """NC PH 4.80 4.25
PH CM 12.00 3.69
PH GL 8.90 3.52
GL EST1 30.00 2.88
GL EST2 8.40 3.48
EST1 CM 4.70 3.55
EST2 CM 5.10 3.61
NH PHE 9.40 3.60
PHE GL 8.90 3.52
CM CM 6.16 3.64
CM CT 6.16 3.65
CT2 CM 9.00 3.13
CMD2 CM 8.00 3.03
CT CMD2 8.00 3.09"""
ANGLES_2007 = """CM CM CM 1.190 173.0
CM CM CT 1.190 175.0
CT2 CM CM 1.600 172.0
CT2 CM CT2 1.700 173.0
CT CM CT 1.093 175.5
EO EO EOT 3.500 135.0
EO EO EO 3.400 132.0
EOT EO EOT 3.500 134.0
OA EO OA 6.600 131.0
EO EO OA 3.000 131.0
CT EO EO 2.200 145.0
CT2 CM EO 1.300 178.0
CM EO CM 1.800 171.0
CT2 EO CT2 2.600 165.0
OA EO CT2 4.000 146.0
CM CM EO 1.500 172.0
CM EO EO 3.200 146.0
CT2 CM OA 1.800 166.0"""
ANGLES_2010 = """NC PH CM 3.300 112.0
PH CM CM 1.100 178.0
NC PH GL 3.100 112.0
PH GL EST1 1.400 124.0
PH GL EST2 2.000 138.0
GL EST1 CM 0.800 168.0
GL EST2 CM 0.800 172.0
EST1 GL EST2 1.000 95.0
EST1 CM CM 1.000 178.0
EST2 CM CM 1.000 178.0
NH PHE GL 4.000 102.0
PHE GL EST1 1.400 124.0
PHE GL EST2 2.000 138.0
CM CM CM 1.190 173.0
CM CM CT 1.190 175.0
CT2 CM CM 1.600 172.0
CT2 CM CT 1.600 172.0
CT2 CM CT2 1.700 173.0
CT CM CT 1.093 175.5
CT CMD2 CT 7.700 116.0
CT2 CMD2 CT2 12.000 110.0
CMD2 CM CM 1.900 161.0
CM CMD2 CM 6.000 110.0"""
# all-atom fragments onto beads (CAPS's SMARTS for the segments of 2007 Table 1 and 2010 Table 2), in the order tried;
# alkyl chains tile into CT / CT2 ends with CM between
BEADS = [
    {"type": "NC", "smarts": "[CH2X4][CH2X4][NX4+]([CH3X4])([CH3X4])[CH3X4]", "description": "-CH2CH2-N(CH3)3, +1 (choline)"},
    {"type": "NH", "smarts": "[CH2X4][CH2X4][N;$([NX4+;H3]),$([NX3;H2])]", "description": "-CH2CH2-NH3, +1 (ethanolamine)"},
    {"type": "PHE", "smarts": "[OX2;$(O[CH2X4][CH2X4][N;$([NX4+;H3]),$([NX3;H2])])][PX4](~[OX1])(~[OX1])~[OX2]",
     "description": "-PO4-, -1, of a PE head group"},
    {"type": "PH", "smarts": "[OX2][PX4](~[OX1])(~[OX1])~[OX2]", "description": "-PO4-, -1"},
    {"type": "GL", "smarts": "[CH2X4;$(*O)][CH1X4;$(*O)][CH2X4;$(*O)]", "description": "-CH2CH-CH2- (glycerol)"},
    {"types": ["EST1", "EST2"], "smarts": "[CH2X4][CX3](=[OX1])[OX2]",
     "description": "-CH2CO2- (the two ester beads, told apart by their bond lengths to GL, 2.88 and 3.48 A)"},
    {"type": "CMD2", "smarts": "[CH1X3]=[CH1X3]", "description": "-HC=CH- (cis)"},
    {"type": "CT", "smarts": "[CH3X4][CH2X4][CH2X4]", "description": "CH3CH2CH2-"},
    {"type": "CT2", "smarts": "[CH3X4][CH2X4]", "description": "CH3CH2-"},
    {"type": "EOT", "smarts": "[CH3X4][OX2][CH2X4]", "description": "CH3-O-CH2-"},
    {"type": "EO", "smarts": "[CH2X4][OX2][CH2X4]", "description": "-CH2-O-CH2-"},
    {"type": "OA", "smarts": "[CH2X4][OX2H1]", "description": "HOCH2-"},
    {"type": "CM", "smarts": "[CH2X4][CH2X4][CH2X4]", "description": "-CH2CH2CH2-"},
]


def rows(text):
    return [l.split() for l in text.strip().splitlines()]


def main():
    types = [{"name": t, "element": "", "mass": round(sum(M[e] * n for e, n in f.items()), 4), "charge": float(q),
              "description": f"{d} ({S2007 if p == '2007' else S2010})"} for t, (f, q, d, p) in SEG.items()]
    pairs, seen = [], {}
    for text, src in ((PAIRS_2007, "2007 Table 6"), (PAIRS_2010, "2010 Table 5")):
        for a, b, form, eps, sig in rows(text):
            key = tuple(sorted((a, b)))
            val = (form, float(eps), float(sig))
            if key in seen:
                assert seen[key] == val, (key, seen[key], val)   # the two papers agree where both list a pair
                continue
            seen[key] = val
            pairs.append({"name": f"{a}-{b}", "match": [a] if a == b else [a, b], "style": "lj9_6" if form == "9-6" else "lj12_4",
                          "params": [float(eps), float(sig)], "comment": f"{src}: LJ{form}"})
    bonds, bseen = [], {}
    for text, src in ((BONDS_2007, "2007 Table 4"), (BONDS_2010, "2010 Table 3")):
        for a, b, k, r0 in rows(text):
            key = tuple(sorted((a, b)))
            if key in bseen:
                assert bseen[key] == (float(k), float(r0)), (key, bseen[key], k, r0)
                continue
            bseen[key] = (float(k), float(r0))
            bonds.append({"name": f"{a}-{b}", "match": [a, b], "style": "harmonic", "params": [float(k), float(r0)], "comment": src})
    angles, aseen = [], {}
    for text, src in ((ANGLES_2007, "2007 Table 5"), (ANGLES_2010, "2010 Table 4")):
        for a, b, c, k, t0 in rows(text):
            key = (b, tuple(sorted((a, c))))
            if key in aseen:
                assert aseen[key] == (float(k), float(t0)), (key, aseen[key], k, t0)
                continue
            aseen[key] = (float(k), float(t0))
            angles.append({"name": f"{a}-{b}-{c}", "match": [a, b, c], "style": "sdk", "params": [float(k), float(t0)],
                           "comment": src + " (with the 1-3 correction of 2010 eq 3)"})
    doc = {"format": "caps-forcefield", "format_version": 1, "name": "SDK (Shinoda, DeVane & Klein 2007, 2010)",
           "source": S2007 + ", Tables 1, 4-6; " + S2010 + ", Tables 2-5", "references": [S2007, S2010], "units": "real",
           "styles": {"pair": "lj/sdk/coul/long", "bond": "harmonic", "angle": "sdk", "dihedral": "fourier", "improper": "cvff"},
           "mixing": "arithmetic", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 15.0, "timestep": 10,
           "pair_settings": {"model_cutoff": True, "dielectric": 80.0}, "angle_terms": "if_defined", "torsion_terms": "if_defined",
           "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
           "bead_templates": {"W": "[W]", "DMPC": "[NC+1][PH-1][GL]([EST1][CM][CM][CM][CT2])[EST2][CM][CM][CM][CT2]",
                              "DPPC": "[NC+1][PH-1][GL]([EST1][CM][CM][CM][CM][CT])[EST2][CM][CM][CM][CM][CT]",
                              "DMPE": "[NH+1][PHE-1][GL]([EST1][CM][CM][CM][CT2])[EST2][CM][CM][CM][CT2]",
                              "C12E8": "[OA]" + "[EO]" * 8 + "[CM][CM][CM][CT2]", "DOD": "[CT][CM][CM][CT]"},
           "notes": ["built by CAPS (bench/ff/make_sdk.py) from the 2007 and 2010 papers' tables",
                     "LJ12-4 for pairs with W, LJ9-6 for the others; every pair as listed (no mixing); cut-off 15 A",
                     "Coulomb without truncation (Ewald, as the 2010 paper); relative permittivity 80 is CAPS's choice (the papers do not print it)",
                     "masses from the segments' formulas"],
           "version": "2007 / 2010", "status": "validated", "typing": f"../typing/{FID}.typing.json"}
    json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1, ensure_ascii=False)
    typing = {"format": "caps-typing", "version": 1, "forcefield": doc["name"], "coarse_grained": True, "unknown_types": "untyped",
              "description": "CAPS rules for SDK (bench/ff/make_sdk.py): beads typed by name; all-atom structures mapped onto the papers' "
                             "segments (every heavy atom in one bead, hydrogens with their atom)",
              "beads": BEADS, "bead_groups": [{"type": "W", "molecule": "[OH2]", "count": 3, "description": "(H2O)3: three nearest waters"}],
              "rules": [{"type": t, "smarts": "*", "atom_name": t, "priority": 0, "description": f"a bead named {t}"} for t in SEG]}
    json.dump(typing, open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1, ensure_ascii=False)
    e = {"id": FID, "name": doc["name"], "version": doc["version"], "references": doc["references"], "origin": "built from source",
         "source": doc["source"] + "; bench/ff/make_sdk.py", "license": "published parameters of the papers", "status": "validated",
         "notes": f"{len(types)} segments (water, alkanes, PEG, phospholipids), {len(pairs)} pairs as tabulated; LAMMPS equals CAPS (bench/ff/check_data_lammps.py)",
         "file": FID + ".json", "typing": {"rules": f"typing/{FID}.typing.json", "evidence": "beads by name; all-atom lipids, alkanes and PEG mapped by segment"},
         "counts": {"atom_types": len(types), "pairs": len(pairs), "bonds": len(bonds), "angles": len(angles), "dihedrals": 0, "impropers": 0},
         "group": "cg", "aliases": ["sdk", "spica"], "list": {"label": "SDK (Shinoda-DeVane-Klein)", "key": "sdk", "order": 39}}
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    cat["forcefields"] = [x for x in cat["forcefields"] if x["id"] != FID] + [e]
    json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)
    print(f"{FID} written: {len(types)} types, {len(pairs)} pairs, {len(bonds)} bonds, {len(angles)} angles")


if __name__ == "__main__":
    main()
