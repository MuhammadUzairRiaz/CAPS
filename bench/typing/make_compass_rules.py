#!/usr/bin/env python3
"""CAPS typing rules for COMPASS (Sun, J. Phys. Chem. B 102, 7338 (1998)) in the library's two COMPASS files:
compass (diagonal terms, with extension types for alkenes, ketones, amines, halogens ...) and
compass-published-moltemplate (the published subset with every class II cross term).

The rules follow the published type definitions: c4 a generic sp3 carbon, c43 / c44 sp3 carbons with three / four
heavy atoms attached, c4o the carbon alpha to an oxygen, c41o / c43o the carbons of methanol / secondary alcohols,
c3a aromatic carbon, c3' a carbonyl carbon with one polar substituent (ester, acid, amide), o2e / o2s / o2h ether,
ester and hydroxyl oxygens, o1= a carbonyl (or NO2 / SO2) oxygen, o2z and si4c siloxanes, h1 / h1o non-polar and
polar hydrogens. Extension types of the diagonal file are ranked the same way. One rule set serves both files:
"unknown_types": "untyped", so a type the published subset lacks (ketone c3o, amide h1n, halogens) leaves its atoms
untyped instead of letting a more general rule take them. Charges come from each file's bond increments.

usage: make_compass_rules.py [DATA_DIR]
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data")
R = []


def add(t, prio, smarts, desc):
    R.append({"type": t, "smarts": smarts, "priority": prio, "description": desc})


# hydrogens
add("h1", 0, "[H][#6,#14]", "non-polar H (on C or Si)")
add("h1o", 1, "[H][O,F]", "strongly polar H (on O or F)")
add("h1n", 1, "[H][N,Cl]", "H on N")
add("h1h", 2, "[H][H]", "H in H2")
# carbon, sp3
add("c4", 0, "[CX4]", "generic sp3 C")
add("c43", 2, "[CX4H1]", "sp3 C with three heavy atoms attached")
add("c44", 2, "[CX4H0]", "sp3 C with four heavy atoms attached")
add("c4o", 4, "[CX4][OX2]", "sp3 C alpha to an oxygen")
add("c41o", 5, "[CX4H3][OX2H1]", "C of methanol")
add("c43o", 5, "[CX4H1][OX2H1]", "C of a secondary alcohol")
add("c4z", 5, "[CX4][NX2]=[NX2+]=[NX1-]", "sp3 C bonded to -N3")
add("c4x", 3, "[CX4][Cl]", "sp3 C bonded to Cl")
# carbon, sp2 / sp
add("c3a", 1, "c", "aromatic C")
add("c3=", 1, "[CX3]=[CX3]", "sp2 C of C=C")
add("c3'", 4, "[CX3](=[OX1])[OX2,NX3]", "carbonyl C with one polar substituent (ester, acid, amide)")
add("c3o", 3, "[CX3](=[OX1])([#6,#1])[#6,#1]", "carbonyl C of an aldehyde or ketone")
add("c3-", 5, "[CX3](=[OX1])[OX1-]", "carboxylate C")
add("c3n", 3, "[CX3]=[NX2]", "sp2 C bonded to N (imine)")
add("c2t", 1, "[CX2]#[C,N]", "sp C of a triple bond")
add("c2=", 5, "[CX2](=O)=O", "C in CO2")
add("c1o", 5, "[CX1-]#[OX1+]", "C in CO")
# oxygen
add("o2", 0, "[OX2]", "generic O with two bonds")
add("o1=", 1, "[OX1]=[#6,#7,#16]", "carbonyl O (and O in NO2, SO2)")
add("o2e", 2, "[OX2]([#6])[#6]", "ether O")
add("o2s", 3, "[OX2]([#6])[CX3]=[OX1]", "ester O")
add("o2h", 3, "[OX2H1]", "hydroxyl O")
add("o2c", 4, "[OX2H1][CX3]=[OX1]", "hydroxyl O of a carboxylic acid")
add("o1-", 5, "[OX1-][CX3]=[OX1]", "carboxylate O")
add("o2z", 3, "[OX2]([Si])[Si,Al]", "O in siloxanes and zeolites")
add("o12", 4, "[OX1]~[NX3](~[OX1])[#6]", "O in a nitro group")
add("o2n", 4, "[OX2][NX3](~[OX1])~[OX1]", "O in a nitrate")
# nitrogen
add("n3", 0, "[NX3]", "sp3 N in an amine")
add("n3h2", 1, "[NX3H2]c", "N of an aniline")
add("n3mh", 3, "[NX3;H1,H2][CX3]=[OX1]", "amide N with H")
add("n3m", 3, "[NX3H0][CX3]=[OX1]", "amide N without H")
add("n1t", 2, "[NX1]#[#6]", "sp N of a nitrile")
add("n3o", 3, "[NX3](~[OX1])~[OX1]", "N of a nitro group")
add("n2=", 1, "[NX2]=[#6,#7,#8]", "sp2 N (imine, nitroso)")
add("n3a", 2, "[nX3]", "aromatic N with three connections (pyrrole)")
add("n2a", 2, "[nX2]", "aromatic N with two connections (pyridine)")
add("n4+", 3, "[NX4+]", "ammonium N")
# silicon
add("si4", 0, "[SiX4]", "generic Si with four bonds")
add("si4c", 1, "[SiX4H0]", "Si with no H attached (siloxanes)")
# halogens (diagonal file's extension types)
add("f1", 0, "[F][CX4;!$(*(F)F)]", "F on a C with one halogen")
add("f12", 1, "[F][CX4;$(*(F)F);!$(*(F)(F)F)]", "F on a C with two halogens")
add("f13", 1, "[F][CX4;$(*(F)(F)F)]", "F on a C with three halogens")
add("cl1", 0, "[Cl][CX4;!$(*(Cl)[Cl,F,Br,I])]", "Cl on a C with one halogen")
add("cl12", 1, "[Cl][CX4;$(*(Cl)[Cl,F,Br,I]);!$(*(Cl)(Cl)(Cl)Cl)]", "Cl on a C with two or three halogens")
add("cl14", 2, "[Cl][CX4;$(*(Cl)(Cl)(Cl)Cl)]", "Cl on a C with four halogens")

# the published file spells two names differently
SPELL = {"compass-published-moltemplate": {"c3'": "c3prime", "o1=*": "o1=star"}}
LABEL = {"compass": "COMPASS (diagonal terms)", "compass-published-moltemplate": "COMPASS (published subset, full class II)"}

if __name__ == "__main__":
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for target in ("compass", "compass-published-moltemplate"):
        ff_p = os.path.join(DATA, "forcefields", target + ".json")
        ff = json.load(open(ff_p))
        have = {t["name"] for t in ff["atom_types"]} | {t["name"].split("~")[0] for t in ff["atom_types"]}
        sp = SPELL.get(target, {})
        rules = [dict(r, type=sp.get(r["type"], r["type"])) for r in R]
        absent = sorted({r["type"] for r in rules if r["type"] not in have})
        doc = {"format": "caps-typing", "version": 1, "forcefield": LABEL[target],
               "description": "CAPS rules for COMPASS from the published type definitions (Sun, J. Phys. Chem. B 102, 7338 (1998)); among "
                              "matching rules the highest priority wins. Types this file lacks leave their atoms untyped: " + (", ".join(absent) or "none") + ".",
               "unknown_types": "untyped", "rules": rules}
        out = os.path.join(DATA, "typing", target + ".typing.json")
        json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
        ff["typing"] = f"../typing/{target}.typing.json"
        json.dump(ff, open(ff_p, "w"), ensure_ascii=False, indent=1)
        for e in cat["forcefields"]:
            if e["id"] == target:
                e["typing"] = {"rules": f"typing/{target}.typing.json", "evidence": "CAPS COMPASS rules from the published type definitions (bench/typing/make_compass_rules.py)"}
        print(f"{out}: {len(rules)} rules; absent from {target}: {', '.join(absent) or 'none'}")
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
