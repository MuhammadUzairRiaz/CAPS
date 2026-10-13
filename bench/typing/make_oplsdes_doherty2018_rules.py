#!/usr/bin/env python3
"""CAPS typing rules for OPLS-DES (Doherty & Acevedo, J. Phys. Chem. B 122, 9982 (2018); bench/ff/make_opls_des_doherty2018.py),
following the SI's atom-type schemes: Scheme S1 (choline) and Scheme S2 (ethylene glycol, glycerol, levulinic acid,
malonic acid, oxalic acid, phenol, urea), plus chloride.

OPLS-DES is a set of whole molecules, not a transferable atom typing: each HBD's types carry its own charges and
Lennard-Jones terms. Every rule therefore holds only inside its molecule: `in_molecule` is a SMARTS that spells the
whole molecule with every heavy atom's connectivity and H count fixed, so it matches that molecule and nothing larger.
Anything else stays untyped.

  choline          NA, CA (methyls) with HA, CS (CH2 on N) with HS, CW (CH2 on O) with HW, OY, HY
  chloride         Cl
  ethylene glycol  OG (Scheme S2's OH), HO, CG, HG
  glycerol         CM (CH), HZ, OM on CM, HM; CB (CH2), HC, OH on CB, HO
  levulinic acid   CB (carboxyl C), OB (=O), OH, HO, CZ / HZ (CH2 on CB), CF / HF (CH2 on the ketone C), CD (ketone C),
                   OC, CT / HT (methyl)
  malonic acid     CD, OD, OH, HO, CT, HC
  oxalic acid      CD, OD, OH, HO
  phenol           CO (ring C holding the OH), CB, HB, OH, HO
  urea             C, O, N, HT; Scheme S2's HC (the H cis to O) has the same terms as HT and a graph cannot tell the two
                   apart, so every urea H is HT

usage: make_oplsdes_doherty2018_rules.py [OUT_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data", "typing")
R = []

# whole molecules, every heavy atom's degree and H count fixed
WHOLE = {
    "choline": "[CX4H3][#7X4]([CX4H3])([CX4H3])[CX4H2][CX4H2][OX2H1]",
    "ethylene glycol": "[OX2H1][CX4H2][CX4H2][OX2H1]",
    "glycerol": "[OX2H1][CX4H2][CX4H1]([OX2H1])[CX4H2][OX2H1]",
    "levulinic acid": "[CX4H3][CX3](=[OX1])[CX4H2][CX4H2][CX3](=[OX1])[OX2H1]",
    "malonic acid": "[OX2H1][CX3](=[OX1])[CX4H2][CX3](=[OX1])[OX2H1]",
    "oxalic acid": "[OX2H1][CX3](=[OX1])[CX3](=[OX1])[OX2H1]",
    "phenol": "[OX2H1][cX3]1[cX3H1][cX3H1][cX3H1][cX3H1][cX3H1]1",
    "urea": "[NX3H2][CX3](=[OX1])[NX3H2]",
}


def add(t, mol, s, d):
    R.append({"type": t, "smarts": s, "priority": 5, "description": f"{d} ({mol})", "in_molecule": [WHOLE[mol]]})


# choline (Scheme S1)
m = "choline"
add("NA", m, "[#7X4]", "N")
add("CA", m, "[CX4H3]", "methyl C")
add("CS", m, "[CX4H2]-[#7]", "CH2 on N")
add("CW", m, "[CX4H2]-[#8]", "CH2 on O")
add("OY", m, "[OX2]", "hydroxyl O")
add("HY", m, "[H][#8]", "hydroxyl H")
add("HA", m, "[H][CX4H3]", "H on a methyl")
add("HS", m, "[H][CX4H2]-[#7]", "H on the CH2 on N")
add("HW", m, "[H][CX4H2]-[#8]", "H on the CH2 on O")
# chloride (Table S1)
R.append({"type": "Cl", "smarts": "[Cl;X0]", "priority": 5, "description": "chloride"})
# ethylene glycol (Scheme S2)
m, x = "ethylene glycol", "_etg"
add("OG" + x, m, "[OX2]", "O (OH in Scheme S2)")
add("HO" + x, m, "[H][#8]", "hydroxyl H")
add("CG" + x, m, "[#6]", "C")
add("HG" + x, m, "[H][#6]", "H on C")
# glycerol
m, x = "glycerol", "_gly"
add("CM" + x, m, "[CX4H1]", "central CH")
add("CB" + x, m, "[CX4H2]", "terminal CH2")
add("OM" + x, m, "[OX2]-[CX4H1]", "O on the central C")
add("OH" + x, m, "[OX2]-[CX4H2]", "O on a terminal C")
add("HM" + x, m, "[H][#8]-[CX4H1]", "H on the central O")
add("HO" + x, m, "[H][#8]-[CX4H2]", "H on a terminal O")
add("HZ" + x, m, "[H][CX4H1]", "H on the central C")
add("HC" + x, m, "[H][CX4H2]", "H on a terminal C")
# levulinic acid
m, x = "levulinic acid", "_lev"
add("CB" + x, m, "[CX3]-[OX2]", "carboxyl C")
add("OB" + x, m, "[OX1]=[CX3]-[OX2]", "carboxyl =O")
add("OH" + x, m, "[OX2]", "carboxyl OH")
add("HO" + x, m, "[H][#8]", "carboxyl H")
add("CD" + x, m, "[CX3;!$(*-[OX2])]", "ketone C")
add("OC" + x, m, "[OX1]=[CX3;!$(*-[OX2])]", "ketone O")
add("CZ" + x, m, "[CX4H2]-[CX3]-[OX2]", "CH2 on the carboxyl C")
add("CF" + x, m, "[CX4H2]-[CX3;!$(*-[OX2])]", "CH2 on the ketone C")
add("CT" + x, m, "[CX4H3]", "methyl C")
add("HZ" + x, m, "[H][CX4H2]-[CX3]-[OX2]", "H on CZ")
add("HF" + x, m, "[H][CX4H2]-[CX3;!$(*-[OX2])]", "H on CF")
add("HT" + x, m, "[H][CX4H3]", "H on the methyl")
# malonic acid
m, x = "malonic acid", "_mal"
add("CD" + x, m, "[CX3]", "carboxyl C")
add("OD" + x, m, "[OX1]", "carboxyl =O")
add("OH" + x, m, "[OX2]", "carboxyl OH")
add("HO" + x, m, "[H][#8]", "carboxyl H")
add("CT" + x, m, "[CX4]", "CH2")
add("HC" + x, m, "[H][#6]", "H on the CH2")
# oxalic acid
m, x = "oxalic acid", "_ox"
add("CD" + x, m, "[CX3]", "C")
add("OD" + x, m, "[OX1]", "=O")
add("OH" + x, m, "[OX2]", "OH")
add("HO" + x, m, "[H][#8]", "H")
# phenol
m, x = "phenol", "_phe"
add("CO" + x, m, "[c]-[OX2]", "ring C holding the OH")
add("CB" + x, m, "[cH1]", "ring CH")
add("HB" + x, m, "[H][c]", "ring H")
add("OH" + x, m, "[OX2]", "O")
add("HO" + x, m, "[H][#8]", "hydroxyl H")
# urea
m, x = "urea", "_urea"
add("C" + x, m, "[CX3]", "C")
add("O" + x, m, "[OX1]", "O")
add("N" + x, m, "[NX3]", "N")
add("HT" + x, m, "[H][#7]", "H on N (Scheme S2's HT and HC, the same terms; a graph cannot tell cis from trans)")

json.dump({"format": "caps-typing", "version": 1, "forcefield": "OPLS-DES (Doherty & Acevedo 2018)", "unknown_types": "untyped",
           "description": "CAPS rules for OPLS-DES's molecules as the SI's Schemes S1 and S2 type them, each rule only inside its "
                          "own whole molecule (bench/typing/make_oplsdes_doherty2018_rules.py)",
           "rules": R}, open(os.path.join(OUT, "oplsdes-doherty2018.typing.json"), "w"), indent=1, ensure_ascii=False)
print(f"oplsdes-doherty2018: {len(R)} rules")
