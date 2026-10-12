#!/usr/bin/env python3
"""CAPS typing rules for CL&P, the Canongia Lopes & Padua force field for ionic liquids, built from its own il.ff
(bench/ff/convert_fftool.py: clandp), in CL&P's type names, checked against CL&P's ion files (.zmat / .xyz, whose atoms
carry their types; bench/ff/validate_fftool_types.py).

CL&P types follow the ion family (each with its own charges): the atoms of a cation's head and the first carbons of
its side chains have their own types (C1 a methyl, C1A a CH2 on the charged N, C2 the next CH2, CE an ethyl's CH3, H1
the H on those first carbons), the rest of a chain is OPLS-AA's alkane (CS, CT, HC); every anion has its own set.

usage: make_clandp_rules.py [OUT_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data", "typing")
R = []


def add(t, p, s, d):
    R.append({"type": t, "smarts": s, "priority": p, "description": d})


IM = "[#6X3;r5;$(*1~[#7]~[#6]~[#6]~[#7]~1)]"                       # an imidazolium C2 (between the two N)
IMN = "[#7X3;r5;$(*1~[#6]~[#7]~[#6]~[#6]~1),$(*1~[#6]~[#6]~[#7]~[#6]~1)]"   # an imidazolium N
PYN = "[#7X3;r6;$(*1~[#6]~[#6]~[#6]~[#6]~[#6]~1)]"                  # a pyridinium N
AMN = "[#7X4]"                                                       # an ammonium N (also in a pyrrolidinium ring)
PHN = "[#15X4;!$(*~[#8,#7,#9])]"                                     # a phosphonium P (four C)
CAT = f"[$({IMN}),$({PYN}),$({AMN}),$({PHN})]"                       # a cation's charged N or P
# ---------------------------------------------------------------- alkyl chains (OPLS-AA), lowest priority
add("CT", 0, "[CX4H3]", "alkane CH3")
add("CS", 0, "[CX4H2]", "alkane CH2")
add("HC", 0, "[H][CX4]", "H on an alkane C")
add("CT3", 0, "[CX4H1;!$(*~[!#6;!#1])]", "alkane CH")
add("CT2", 1, "[CX4H2;R;!$(*~[!#6;!#1])]", "CH2 of a carbocycle (cyclohexane carboxylate)")
# ---------------------------------------------------------------- the first carbons on a cation's N or P
add("C1", 3, f"[CX4H3;$(*-{CAT})]", "methyl on the cation's N")
add("C1A", 3, f"[CX4H2;$(*-{CAT})]", "CH2 on the cation's N (side chain or pyrrolidinium ring)")
add("C2", 2, f"[CX4H2;$(*-[CX4H2]-{CAT})]", "CH2 next to the cation's first carbon")
add("CE", 3, f"[CX4H3;$(*-[CX4H2]-{CAT})]", "an ethyl's CH3 on the cation's N (imidazolium, ammonium, pyridinium)")
add("H1", 3, f"[H][CX4;$(*-{CAT})]", "H on a carbon bonded to the cation's N")
add("C1P", 4, f"[CX4;$(*-{PHN})]", "carbon on a phosphonium P")
add("P3", 4, PHN, "phosphonium P")
# ---------------------------------------------------------------- imidazolium ring
add("NA", 3, IMN, "imidazolium N")
add("CR", 3, IM, "imidazolium C2")
add("CW", 2, "[#6X3;r5;$(*1~[#6]~[#7]~[#6]~[#7]~1)]", "imidazolium C4 / C5")
add("HCR", 3, f"[H]{IM}", "H on the imidazolium C2")
add("HCW", 3, "[H][#6X3;r5;$(*1~[#6]~[#7]~[#6]~[#7]~1)]", "H on the imidazolium C4 / C5")
# 1-alkylimidazolium (protonated N3)
add("NAH", 5, "[#7X3H1;r5;$(*1~[#6]~[#7]~[#6]~[#6]~1),$(*1~[#6]~[#6]~[#7]~[#6]~1)]", "the N-H of an alkylimidazolium")
add("HNA", 5, "[H][#7X3H1;r5]", "H on the imidazolium N")
add("CRH", 5, f"[#6X3;r5;$(*(~[#7H1;r5])~[#7;r5])]", "C2 of an alkylimidazolium")
add("CWH", 5, f"[#6X3;r5;$(*~[#7H1;r5]);!$(*(~[#7])~[#7])]", "C4 / C5 next to the N-H")
# 1,2,3-trialkylimidazolium (C2 methylated)
add("CRM", 5, f"[#6X3;r5;$(*(~[#7;r5])~[#7;r5]);$(*-[CX4])]", "alkylated imidazolium C2")
add("NAM", 5, f"[#7X3;r5;$(*~[#6X3;r5;$(*(~[#7;r5])~[#7;r5]);$(*-[CX4])])]", "imidazolium N next to an alkylated C2")
add("CCR", 5, "[CX4H3;$(*-[#6X3;r5;$(*(~[#7;r5])~[#7;r5])])]", "methyl on the imidazolium C2")
add("H1", 5, "[H][CX4H3;$(*-[#6X3;r5;$(*(~[#7;r5])~[#7;r5])])]", "H on the C2 methyl")
# benzyl on the imidazolium N: the ring as toluene
BZ = f"[CX4H2;$(*-{IMN})]"
add("C1T", 6, f"[CX4H2;$(*-{IMN});$(*-c)]", "benzyl CH2 on the imidazolium N")
add("CAT", 6, f"[c;$(*-{BZ})]", "benzyl ipso C")
add("CAO", 6, f"[c;$(*:c-{BZ})]", "benzyl ortho C")
add("CAM", 6, f"[c;$(*:c:c-{BZ})]", "benzyl meta C")
add("CAP", 6, f"[c;$(*:c:c:c-{BZ})]", "benzyl para C")
add("HAT", 6, f"[H][c;$(*:c-{BZ}),$(*:c:c-{BZ}),$(*:c:c:c-{BZ})]", "H on the benzyl ring")
# fluoroalkyl side chain (C1H, CFH, CSF, CTF)
add("C1H", 6, f"[CX4H2;$(*-{IMN});$(*-[CX4H2]-[CX4](F)F)]", "CH2 on the imidazolium N of a fluoroalkyl chain")
add("CFH", 6, "[CX4;$(*(F)F);$(*-[CX4H2])]", "CF2 next to the CH2 of a fluoroalkyl chain")
add("CSF", 4, "[CX4H0;$(*(F)F);!$(*(F)(F)F)]", "perfluoro CF2")
add("CTF", 4, "[CX4H0;$(*(F)(F)F)]", "perfluoro CF3")
add("F", 4, "[F][CX4]", "F on a perfluoroalkyl C")
# hydroxyalkyl side chain
add("C2O", 5, f"[CX4H2;$(*-[OX2H1]);$(*-[CX4]-{CAT})]", "CH2-OH on a cation's side chain")
add("OH", 5, "[OX2H1;$(*-[CX4])]", "hydroxyl O")
add("HO", 5, "[H][OX2H1]", "hydroxyl H")
# ether side chains (OPLS-AA ethers): O OY, its CH2 CSO and CH3 C3O with H HY; a CH2 between the cation's first
# carbon and the O COM (its H HC); the phosphonium CH2 that leads to an ether C1PS
add("OY", 3, "[OX2;$(*(-[CX4])-[CX4])]", "ether O")
add("CSO", 3, "[CX4H2;$(*-[OX2]-[CX4])]", "ether CH2"); add("C3O", 4, "[CX4H3;$(*-[OX2]-[CX4])]", "methoxy CH3")
add("HY", 3, "[H][CX4;$(*-[OX2]-[CX4])]", "H on an ether C")
NCAT = f"[$({IMN}),$({PYN}),$({AMN})]"   # a cation's N (CL&P's COM is for ammonium and pyrrolidinium ether chains)
add("COM", 5, f"[CX4H2;$(*-[OX2]-[CX4]);$(*-[CX4]-{NCAT})]", "ether CH2 next to the cation's first carbon")
add("HC", 6, f"[H][CX4H2;$(*-[OX2]-[CX4]);$(*-[CX4]-{NCAT})]", "H on that CH2")
add("C1PS", 5, f"[CX4H2;$(*-{PHN});$(*-[CX4]-[OX2])]", "phosphonium CH2 leading to an ether")
# isoalkyl phosphonium: the CH2 on P CPI4, the CH CI4, its H HI4
add("CPI4", 5, f"[CX4H2;$(*-{PHN});$(*-[CX4H1])]", "phosphonium CH2 on a branch point")
add("CI4", 5, f"[CX4H1;$(*-[CX4H2]-{PHN})]", "branching CH of an isoalkyl phosphonium")
add("HI4", 5, f"[H][CX4H1;$(*-[CX4H2]-{PHN})]", "H on that CH")
# ---------------------------------------------------------------- pyridinium
add("NAP", 4, PYN, "pyridinium N")
add("CAPO", 4, f"[c;r6;$(*~{PYN})]", "pyridinium ortho C")
add("CAPM", 4, f"[c;r6;$(*~[c]~{PYN})]", "pyridinium meta C")
add("CAPP", 4, f"[c;r6;$(*~[c]~[c]~{PYN})]", "pyridinium para C")
add("HAP", 4, f"[H][c;r6;$(*~{PYN}),$(*~[c]~{PYN}),$(*~[c]~[c]~{PYN})]", "H on the pyridinium ring")
# isoquinolinium: the N NAQ, the ring-fusion C CAQ, the other ring C CA; H on the N ring HAP, on the benzo ring HA
ISQ = "[#7X3;r6;$(*1~[#6]~[#6;R2]~[#6;R2]~[#6]~[#6]~1),$(*1~[#6]~[#6]~[#6;R2]~[#6;R2]~[#6]~1)]"
add("NAQ", 6, ISQ, "isoquinolinium N")
add("CAQ", 6, f"[#6;R2;$(*@*@{ISQ}),$(*@*@*@{ISQ})]", "isoquinolinium ring-fusion C")
NRING = "$(*1~[#7]~*~*~*~*~1),$(*1~*~[#7]~*~*~*~1),$(*1~*~*~[#7]~*~*~1)"   # in the 6-ring that holds the N
R.append({"type": "CA", "smarts": "[c;R1]", "priority": 5, "description": "isoquinolinium ring C", "in_molecule": [ISQ]})
R.append({"type": "HA", "smarts": "[H][c;R1]", "priority": 5, "description": "H on the isoquinolinium benzo ring", "in_molecule": [ISQ]})
R.append({"type": "HAP", "smarts": f"[H][c;R1;{NRING}]", "priority": 6, "description": "H on the isoquinolinium N ring", "in_molecule": [ISQ]})
# ---------------------------------------------------------------- ammonium
add("N4", 3, "[#7X4;H0]", "quaternary ammonium N")
add("N3", 3, "[#7X4;H1]", "protonated tertiary amine N (R3NH+)")
add("H3", 3, "[H][#7X4;H1]", "H on R3NH+")
# ---------------------------------------------------------------- guanidinium
add("CG", 5, "[CX3;$(*(~[NX3])(~[NX3])~[NX3])]", "guanidinium C")
add("NG", 5, "[NX3;$(*~[CX3](~[NX3])~[NX3])]", "guanidinium N")
add("HG", 5, "[H][NX3;$(*~[CX3](~[NX3])~[NX3])]", "guanidinium H")
# ---------------------------------------------------------------- anions
add("B", 5, "[#5;$(*(F)(F)(F)F)]", "tetrafluoroborate B"); add("FB", 5, "[F][#5]", "tetrafluoroborate F")
add("P", 5, "[#15;$(*(F)(F)(F)(F)(F)F)]", "hexafluorophosphate P"); add("FP", 5, "[F][#15]", "hexafluorophosphate F")
add("Cl", 5, "[Cl;X0]", "chloride"); add("Br", 5, "[Br;X0]", "bromide")
add("Na", 5, "[Na]", "sodium cation"); add("Li", 5, "[Li]", "lithium cation")
# sulfonyl imides and triflate: S(=O)2 on N- or O-
SBT = "[#16X4;$(*(~[OX1])~[OX1]);$(*~[NX2]),$(*(~[OX1])(~[OX1])~[OX1])]"
add("SBT", 6, f"[#16X4;$(*~[NX2]~[#16X4]),$(*(~[OX1])(~[OX1])(~[OX1])-[CX4](F)(F)F)]", "S of a sulfonyl imide / triflate")
add("NBT", 6, "[NX2;$(*(~[#16X4])~[#16X4])]", "N of a sulfonyl imide")
add("OBT", 6, "[OX1;$(*~[#16X4]~[NX2])]", "O of a sulfonyl imide")
add("OTF", 6, "[OX1;$(*~[#16X4](~[OX1])(~[OX1])-[CX4](F)(F)F)]", "triflate O")
add("CBT", 6, "[CX4;$(*(F)(F)F);$(*-[#16X4])]", "CF3 on a sulfonyl S")
add("F1", 6, "[F][CX4;$(*-[#16X4])]", "F on a C bonded to the sulfonyl S")
add("FSI", 6, "[F][#16X4]", "F on S (bis(fluorosulfonyl)imide)")
add("C1F", 6, "[CX4;$(*(F)F);!$(*(F)(F)F);$(*-[#16X4])]", "CF2 on a sulfonyl S")
add("CEF", 6, "[CX4;$(*(F)(F)F);$(*-[CX4](F)(F)-[#16X4])]", "CF3 of a pentafluoroethyl sulfonyl")
# dicyanamide, thiocyanate, nitrate
add("N3A", 6, "[NX2;$(*(-[CX2]#[NX1])-[CX2]#[NX1])]", "dicyanamide central N")
add("CZA", 6, "[CX2;$(*#[NX1]);$(*-[NX2]-[CX2]#[NX1])]", "dicyanamide C"); add("NZA", 6, "[NX1;$(*#[CX2]-[NX2]-[CX2]#N)]", "dicyanamide terminal N")
add("SK", 6, "[SX1;$(*~[CX2]~[NX1])]", "thiocyanate S"); add("CK", 6, "[CX2;$(*~[SX1]);$(*~[NX1])]", "thiocyanate C")
add("NK", 6, "[NX1;$(*~[CX2]~[SX1])]", "thiocyanate N")
add("NO", 6, "[NX3;$(*(~[OX1])(~[OX1])~[OX1])]", "nitrate N"); add("ON", 6, "[OX1;$(*~[NX3](~[OX1])~[OX1])]", "nitrate O")
# tricyanomethanide
add("C3A", 6, "[CX3;$(*(-[CX2]#N)(-[CX2]#N)-[CX2]#N)]", "tricyanomethanide central C")
add("CN", 6, "[CX2;$(*#[NX1]);$(*-[CX3])]", "tricyanomethanide nitrile C"); add("NC", 6, "[NX1;$(*#[CX2]-[CX3])]", "tricyanomethanide N")
# carboxylates
add("CO2", 5, "[CX3;$(*(~[OX1])~[OX1])]", "carboxylate C")
add("O2", 5, "[OX1;$(*~[CX3]~[OX1])]", "carboxylate O")
add("CTA", 5, "[CX4H3;$(*-[CX3](~[OX1])~[OX1])]", "acetate CH3")
add("C2C", 5, "[CX4H2;$(*-[CX3](~[OX1])~[OX1])]", "CH2 on a carboxylate")
add("CPH", 6, "[c;$(*-[CX4H2]-[CX3](~[OX1])~[OX1])]", "phenylacetate ipso C")
PA = "[CX4H2]-[CX3](~[OX1])~[OX1]"
add("CAO", 6, f"[c;$(*:c-{PA})]", "phenylacetate ortho C"); add("CAM", 6, f"[c;$(*:c:c-{PA})]", "phenylacetate meta C")
add("CAP", 6, f"[c;$(*:c:c:c-{PA})]", "phenylacetate para C"); add("HA", 6, f"[H][c;$(*:c-{PA}),$(*:c:c-{PA}),$(*:c:c:c-{PA})]", "H on the phenylacetate ring")
add("C3C", 5, "[CX4H1;$(*-[CX3](~[OX1])~[OX1])]", "CH on a carboxylate")
BZ_O = "[CX3](~[OX1])~[OX1]"
add("CBI", 6, f"[c;$(*-{BZ_O})]", "benzoate ipso C"); add("CAO", 6, f"[c;$(*:c-{BZ_O})]", "benzoate ortho C")
add("CAM", 6, f"[c;$(*:c:c-{BZ_O})]", "benzoate meta C"); add("CBP", 6, f"[c;$(*:c:c:c-{BZ_O})]", "benzoate para C")
add("HA", 6, f"[H][c;$(*:c-{BZ_O}),$(*:c:c-{BZ_O})]", "H on the benzoate ring")
add("O2F", 7, "[OX1;$(*~[CX3](~[OX1])-[CX4](F)(F)F)]", "trifluoroacetate O")
add("CFA", 7, "[CX4;$(*(F)(F)F);$(*-[CX3](~[OX1])~[OX1])]", "trifluoroacetate CF3"); add("FFA", 7, "[F][CX4;$(*-[CX3](~[OX1])~[OX1])]", "trifluoroacetate F")
# alkylsulfonates and sulfates, tosylate
add("SO", 5, "[#16X4;$(*(~[OX1])(~[OX1])~[OX1]);!$(*-[CX4](F)(F)F)]", "sulfonate / sulfate S")
add("OS3", 5, "[OX1;$(*~[#16X4](~[OX1])(~[OX1])-[#6])]", "sulfonate O"); add("OS4", 6, "[OX1;$(*~[#16X4](~[OX1])(~[OX1])-[OX2])]", "sulfate terminal O")
add("OC4", 6, "[OX2;$(*-[#16X4]);$(*-[CX4])]", "sulfate ester O")
add("CS3", 5, "[CX4;$(*-[#16X4](~[OX1])(~[OX1])~[OX1])]", "C on a sulfonate S"); add("HS3", 5, "[H][CX4;$(*-[#16X4](~[OX1])(~[OX1])~[OX1])]", "H on that C")
add("CS4", 6, "[CX4;$(*-[OX2]-[#16X4])]", "C on a sulfate O"); add("HS4", 6, "[H][CX4;$(*-[OX2]-[#16X4])]", "H on that C")
TS = "[#16X4](~[OX1])(~[OX1])~[OX1]"
add("CAPS", 7, f"[c;$(*-{TS})]", "tosylate C on S"); add("CAMS", 7, f"[c;$(*:c-{TS})]", "tosylate C next to it")
add("CAOS", 7, f"[c;$(*:c:c-{TS})]", "tosylate C next to the methyl C"); add("CATS", 7, f"[c;$(*:c:c:c-{TS})]", "tosylate C bearing the methyl")
add("CTTS", 7, f"[CX4H3;$(*-c:c:c:c-{TS})]", "tosylate methyl"); add("HATS", 7, f"[H][c;$(*:c-{TS}),$(*:c:c-{TS})]", "H on the tosylate ring")
add("HTS", 7, f"[H][CX4H3;$(*-c:c:c:c-{TS})]", "H on the tosylate methyl")
json.dump({"format": "caps-typing", "version": 1, "forcefield": "CL&P (ionic liquids)", "unknown_types": "untyped",
           "description": "CAPS rules for CL&P's types by ion family (bench/typing/make_clandp_rules.py), checked against CL&P's own ion files",
           "rules": R}, open(os.path.join(OUT, "clandp.typing.json"), "w"), indent=1, ensure_ascii=False)
print(f"clandp: {len(R)} rules")
