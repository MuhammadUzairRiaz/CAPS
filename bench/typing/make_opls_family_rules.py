#!/usr/bin/env python3
"""CAPS typing rules for two OPLS-AA descendants in the library, from each one's type definitions:

  opls-clp   CL&P (Canongia Lopes & Pádua, J. Phys. Chem. B 108, 2038 (2004) and later): ionic-liquid ions by their
             structure (imidazolium, pyridinium, ammonium, phosphonium, guanidinium cations with their alkyl chains;
             NTf2, triflate / alkylsulfonate, alkyl sulfate, thiocyanate, dicyanamide, tricyanomethanide, PF6, BF4,
             carboxylate and halide anions)
  opls-des   OPLS-DES (Doherty & Acevedo, J. Phys. Chem. B 122, 9982 (2018)): the types belong to named components, so
             each component is recognised as a whole molecule (choline, ethylene glycol, glycerol, levulinic, malonic and
             oxalic acids, phenol, urea, chloride)

Neither file carries per-type charges: automatic charges are Gasteiger-Marsili (else QEq), as the Field
report says. Atoms no rule covers stay untyped and are reported.

(OPLS-AA/M, a protein force field, was removed on 2026-10-14: CAPS focuses on polymers, crystals and materials.)

usage: make_opls_family_rules.py [DATA_DIR]
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data")


class Rules(list):
    def add(self, t, prio, smarts, desc):
        self.append({"type": t, "smarts": smarts, "priority": prio, "description": desc})


def clp():
    R = Rules()
    IMID_N = "[#7;r5;$(*~[#6;r5]~[#7;r5])]"
    CATION_N = f"[$([#7+]),$({IMID_N})]"
    R.add("CT", 0, "[CX4]", "alkyl carbon"); R.add("HC", 0, "[H][CX4]", "alkyl H")
    R.add("CA", 0, "c", "aromatic carbon"); R.add("HA", 0, "[H]c", "aromatic H")
    # a chain carbon on the cationic N (not a ring carbon, not a benzyl carbon); phosphonium's carbons stay CT
    ONN = f"$(*{CATION_N});!R;!$(*c)"
    R.add("C1", 2, f"[CX4;{ONN}]", "chain carbon on the cationic nitrogen")
    R.add("H1", 2, f"[H][CX4;{ONN}]", "H on a chain carbon on the cationic nitrogen")
    # imidazolium
    R.add("NA", 3, IMID_N, "imidazolium nitrogen"); R.add("CR", 3, "[#6;r5](~[#7;r5])~[#7;r5]", "imidazolium C2")
    R.add("CW", 2, f"[#6;r5;$(*~{IMID_N});!$(*(~[#7;r5])~[#7;r5])]", "imidazolium C4 / C5")
    R.add("HCR", 3, "[H][#6;r5](~[#7;r5])~[#7;r5]", "H on imidazolium C2"); R.add("HCW", 2, f"[H][#6;r5;$(*~{IMID_N})]", "H on imidazolium C4 / C5")
    R.add("HNA", 3, "[H][#7;r5]", "H on an imidazolium nitrogen")
    # pyridinium
    R.add("NA", 3, "[#7+;r6]", "pyridinium nitrogen"); PYRC = "[#6X3;r6;$(*@[#7+;r6]),$(*@*@[#7+;r6]),$(*@*@*@[#7+;r6])]"
    R.add("CA", 2, PYRC, "pyridinium ring carbon")
    R.add("HA", 2, f"[H]{PYRC}", "H on a pyridinium ring carbon")
    # ammonium, phosphonium, guanidinium
    R.add("NT", 3, "[NX4]", "ammonium / pyrrolidinium nitrogen"); R.add("HN", 3, "[H][NX4]", "H on an ammonium nitrogen")
    R.add("PT", 3, "[PX4+]", "phosphonium phosphorus")
    R.add("CG", 4, "[CX3;$(C(~[#7])(~[#7])~[#7])]", "guanidinium carbon"); R.add("NG", 4, "[#7;$(*~[CX3](~[#7])~[#7])]", "guanidinium nitrogen")
    R.add("HG", 4, "[H][#7;$(*~[CX3](~[#7])~[#7])]", "H on a guanidinium nitrogen")
    # anions
    R.add("NB", 4, "[#7X2;$(*(~S)~S)]", "bis(sulfonyl)amide nitrogen (NTf2)"); R.add("SB", 4, "[SX4;$(S~[#7X2])]", "sulfonyl sulfur of NTf2")
    R.add("OB", 4, "[OX1]~[SX4]~[#7X2]", "sulfonyl oxygen of NTf2")
    R.add("CF", 4, "[CX4;$(*(F)F)]", "fluorinated carbon (CF2, CF3)")
    R.add("FB", 4, "[F][CX4]S", "F on a carbon on sulfur"); R.add("FB", 4, "[F][SX4]", "F on sulfur (FSI)"); R.add("F", 1, "[F][CX4]", "F on carbon")
    # triflate: the sulfonyl of NTf2 (SB / OB), unlike an alkylsulfonate's SO / OS
    R.add("SB", 5, "[SX4;$(S[CX4](F)(F)F);$(S(~[OX1])(~[OX1])~[OX1])]", "triflate sulfur")
    R.add("OB", 5, "[OX1]~[SX4;$(S[CX4](F)(F)F)](~[OX1])~[OX1]", "triflate oxygen")
    R.add("SO", 3, "[SX4;$(S(~[OX1])(~[OX1])~[OX1])]", "sulfonate / sulfate sulfur"); R.add("OS", 3, "[OX1]~[SX4](~[OX1])~[OX1]", "sulfonate / sulfate oxygen")
    R.add("CS3", 3, "[CX4;!$(*(F)(F)F)][SX4](~[OX1])(~[OX1])~[OX1]", "alkyl carbon on a sulfonate"); R.add("HS3", 3, "[H][CX4][SX4](~[OX1])(~[OX1])~[OX1]", "H on a sulfonate alkyl carbon")
    R.add("OC", 4, "[OX2]([#6])[SX4](~[OX1])(~[OX1])~[OX1]", "ester oxygen of an alkyl sulfate")
    R.add("CS4", 4, "[CX4][OX2][SX4](~[OX1])(~[OX1])~[OX1]", "alkyl carbon of an alkyl sulfate"); R.add("HS4", 4, "[H][CX4][OX2][SX4]", "H on an alkyl sulfate carbon")
    R.add("SK", 4, "[SX1]~[CX2]~[NX1]", "thiocyanate sulfur"); R.add("CK", 4, "[CX2](~[SX1])~[NX1]", "thiocyanate carbon"); R.add("NK", 4, "[NX1]~[CX2]~[SX1]", "thiocyanate nitrogen")
    R.add("N3", 4, "[NX2](~[CX2]~[NX1])~[CX2]~[NX1]", "dicyanamide central nitrogen"); R.add("CZ", 4, "[CX2](~[NX1])~[NX2]", "dicyanamide carbon")
    R.add("NZ", 4, "[NX1]~[CX2]~[NX2]", "dicyanamide terminal nitrogen")
    R.add("C3A", 4, "[CX3](C#N)(C#N)C#N", "tricyanomethanide central carbon")
    R.add("CN", 2, "[CX2]#[NX1]", "nitrile carbon"); R.add("NC", 2, "[NX1]#[CX2]", "nitrile nitrogen")
    R.add("B", 3, "[B]", "borate boron"); R.add("FP", 3, "[F][B]", "tetrafluoroborate fluorine"); R.add("P", 3, "[PX6]", "hexafluorophosphate phosphorus"); R.add("FP", 3, "[F][PX6]", "hexafluorophosphate fluorine")
    R.add("CO", 3, "[CX3](=[OX1])[OX1-]", "carboxylate carbon"); R.add("O2", 3, "[OX1]~[CX3]~[OX1]", "carboxylate oxygen")
    R.add("OH", 2, "[OX2H1]", "alcohol oxygen"); R.add("HO", 2, "[H][OX2]", "alcohol H")
    R.add("Li+", 0, "[Li]", "lithium ion"); R.add("Na+", 0, "[Na]", "sodium ion"); R.add("Cl-", 0, "[Cl-;X0]", "chloride ion"); R.add("Br-", 0, "[Br-;X0]", "bromide ion")
    R.add("OT3", 8, "[OX2H2]", "TIP3P water O"); R.add("HT3", 8, "[H][OX2H2]", "TIP3P water H")
    return R


def des():
    R = Rules()
    CHOL = "$([NX4+](C)(C)(C)CC[OX2H1])"                       # choline: the whole cation from its nitrogen
    R.add("NA1", 5, f"[{CHOL}]", "choline nitrogen")
    R.add("CA1", 5, "[CX4H3;$(*[NX4+](C)(C)CC[OX2H1])]", "choline N-methyl carbon"); R.add("HA1", 5, "[H][CX4H3;$(*[NX4+](C)(C)CC[OX2H1])]", "H on a choline N-methyl")
    R.add("CS1", 5, f"[CX4H2;$(*([NX4+](C)(C)C)C[OX2H1])]", "choline CH2 on nitrogen"); R.add("HS1", 5, f"[H][CX4H2;$(*([NX4+](C)(C)C)C[OX2H1])]", "H on the choline CH2 on nitrogen")
    R.add("CW1", 5, f"[CX4H2;$(*([OX2H1])C[NX4+](C)(C)C)]", "choline CH2 on oxygen"); R.add("HW1", 5, f"[H][CX4H2;$(*([OX2H1])C[NX4+](C)(C)C)]", "H on the choline CH2 on oxygen")
    R.add("OY1", 5, "[OX2H1;$(*CC[NX4+](C)(C)C)]", "choline hydroxyl oxygen"); R.add("HY1", 5, "[H][OX2;$(*CC[NX4+](C)(C)C)]", "choline hydroxyl H")
    EG = "$([CX4H2]([OX2H1])[CX4H2][OX2H1])"                   # ethylene glycol carbon
    R.add("CG2", 4, f"[{EG}]", "ethylene glycol carbon"); R.add("HG2", 4, f"[H][{EG}]", "H on an ethylene glycol carbon")
    R.add("OH2", 4, "[OX2H1;$(*[CX4H2][CX4H2][OX2H1])]", "ethylene glycol oxygen")
    R.add("HO2", 4, "[H][OX2;$(*[CX4H2][CX4H2][OX2H1])]", "ethylene glycol hydroxyl H")
    GLY_END = "$([CX4H2]([OX2H1])[CX4H1]([OX2H1])[CX4H2][OX2H1])"
    GLY_MID = "$([CX4H1]([OX2H1])([CX4H2][OX2H1])[CX4H2][OX2H1])"
    R.add("CB3", 5, f"[{GLY_END}]", "glycerol CH2"); R.add("CM3", 5, f"[{GLY_MID}]", "glycerol CH")
    R.add("HC3", 5, f"[H][{GLY_END},{GLY_MID}]", "H on a glycerol carbon")
    R.add("OH3", 5, f"[OX2H1;$(*[{GLY_END},{GLY_MID}])]", "glycerol oxygen"); R.add("HO3", 5, f"[H][OX2;$(*[{GLY_END},{GLY_MID}])]", "glycerol hydroxyl H")
    # levulinic acid CH3-C(=O)-CH2-CH2-C(=O)OH, atoms by position
    R.add("CT4", 5, f"[CX4H3;$(*C(=O)[CX4H2][CX4H2]C(=O)[OX2H1])]", "levulinic acid methyl carbon"); R.add("HT4", 5, f"[H][CX4H3;$(*C(=O)[CX4H2][CX4H2]C(=O)[OX2H1])]", "H on the levulinic methyl")
    R.add("CD4", 5, f"[CX3;$(*([CX4H3])(=O)[CX4H2][CX4H2]C(=O)[OX2H1])]", "levulinic ketone carbon"); R.add("OC4", 5, f"[OX1;$(*=C([CX4H3])[CX4H2][CX4H2]C(=O)[OX2H1])]", "levulinic ketone oxygen")
    R.add("CF4", 5, "[CX4H2;$(*(C(=O)[CX4H3])[CX4H2]C(=O)[OX2H1])]", "levulinic CH2 next to the ketone"); R.add("HF4", 5, "[H][CX4H2;$(*(C(=O)[CX4H3])[CX4H2]C(=O)[OX2H1])]", "H on it")
    R.add("CZ4", 5, "[CX4H2;$(*(C(=O)[OX2H1])[CX4H2]C(=O)[CX4H3])]", "levulinic CH2 next to the acid"); R.add("HZ4", 5, "[H][CX4H2;$(*(C(=O)[OX2H1])[CX4H2]C(=O)[CX4H3])]", "H on it")
    R.add("CB4", 5, f"[CX3;$(*(=O)([OX2H1])[CX4H2][CX4H2]C(=O)[CX4H3])]", "levulinic acid carbon"); R.add("OB4", 5, f"[OX1;$(*=C([OX2H1])[CX4H2][CX4H2]C(=O)[CX4H3])]", "levulinic acid C=O oxygen")
    R.add("OH4", 5, f"[OX2H1;$(*C(=O)[CX4H2][CX4H2]C(=O)[CX4H3])]", "levulinic acid OH oxygen"); R.add("HO4", 5, f"[H][OX2;$(*C(=O)[CX4H2][CX4H2]C(=O)[CX4H3])]", "levulinic acid H")
    MAL = "$([CX3](=O)([OX2H1])[CX4H2]C(=O)[OX2H1])"
    R.add("CD5", 5, f"[{MAL}]", "malonic acid carbon"); R.add("OD5", 5, "[OX1;$(*=C([OX2H1])[CX4H2]C(=O)[OX2H1])]", "malonic acid C=O oxygen")
    R.add("OH5", 5, "[OX2H1;$(*C(=O)[CX4H2]C(=O)[OX2H1])]", "malonic acid OH oxygen"); R.add("HO5", 5, "[H][OX2;$(*C(=O)[CX4H2]C(=O)[OX2H1])]", "malonic acid H")
    R.add("CT5", 5, "[CX4H2;$(*(C(=O)[OX2H1])C(=O)[OX2H1])]", "malonic CH2"); R.add("HC5", 5, "[H][CX4H2;$(*(C(=O)[OX2H1])C(=O)[OX2H1])]", "H on the malonic CH2")
    OX = "$([CX3](=O)([OX2H1])[CX3](=O)[OX2H1])"
    R.add("CD6", 5, f"[{OX}]", "oxalic acid carbon"); R.add("OD6", 5, "[OX1;$(*=C([OX2H1])C(=O)[OX2H1])]", "oxalic acid C=O oxygen")
    R.add("OH6", 5, "[OX2H1;$(*C(=O)C(=O)[OX2H1])]", "oxalic acid OH oxygen"); R.add("HO6", 5, "[H][OX2;$(*C(=O)C(=O)[OX2H1])]", "oxalic acid H")
    PHENOL = "$(c1ccccc1[OX2H1]),$(c1cccc(c1)[OX2H1]),$(c1ccc(cc1)[OX2H1])"
    R.add("CB7", 4, f"[c;{PHENOL}]", "phenol ring carbon"); R.add("HB7", 4, f"[H][c;{PHENOL}]", "phenol ring H")
    R.add("OH7", 5, "[OX2H1]c", "phenol oxygen"); R.add("HO7", 5, "[H][OX2]c", "phenol hydroxyl H")
    R.add("CB7", 5, "c[OX2H1]", "phenol ring carbon bearing OH")
    UREA = "$([CX3](=O)([NX3H2])[NX3H2])"
    R.add("C8", 5, f"[{UREA}]", "urea carbon"); R.add("O8", 5, "[OX1;$(*=C([NX3H2])[NX3H2])]", "urea oxygen")
    R.add("N8", 5, "[NX3H2;$(*C(=O)[NX3H2])]", "urea nitrogen"); R.add("H8", 5, "[H][NX3H2;$(*C(=O)[NX3H2])]", "urea H")
    R.add("CL", 0, "[Cl-;X0]", "chloride ion (the file's Lennard-Jones entry is CL)")
    return R


ANALOGIES = {
             # triflate and alkylsulfonate sulfonyls: the file lists S-CF3 only for NTf2's sulfonyl group
             "opls-clp": {"SO": ["SB"], "OS": ["OB"]}}

SETS = {"opls-clp": (clp(), "CL&P ionic liquids (OPLS-AA based)"), "opls-des": (des(), "OPLS-DES deep eutectic solvents")}

if __name__ == "__main__":
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for target, (rules, label) in SETS.items():
        ff_p = os.path.join(DATA, "forcefields", target + ".json")
        ff = json.load(open(ff_p))
        have = {t["name"] for t in ff["atom_types"]}
        absent = sorted({r["type"] for r in rules if r["type"] not in have})
        doc = {"format": "caps-typing", "version": 1, "forcefield": label,
               "description": f"CAPS rules for {label} from the force field's type definitions (bench/typing/make_opls_family_rules.py); "
                              "among matching rules the highest priority wins. Types this file lacks leave their atoms untyped: " + (", ".join(absent) or "none") + ".",
               "unknown_types": "untyped", "rules": list(rules)}
        if target in ANALOGIES:
            doc["analogies"] = ANALOGIES[target]
            doc["analogy_source"] = f"{label} analogues: the same element and hybridisation, nearest in hydrogen count or neighbourhood"
        out = os.path.join(DATA, "typing", target + ".typing.json")
        json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
        ff["typing"] = f"../typing/{target}.typing.json"
        json.dump(ff, open(ff_p, "w"), ensure_ascii=False, indent=1)
        for e in cat["forcefields"]:
            if e["id"] == target:
                e["typing"] = {"rules": f"typing/{target}.typing.json", "evidence": f"CAPS rules from {label}'s type definitions (bench/typing/make_opls_family_rules.py)"}
        print(f"{out}: {len(rules)} rules; absent from {target}: {', '.join(absent) or 'none'}")
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
