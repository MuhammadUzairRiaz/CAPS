#!/usr/bin/env python3
"""Write CAPS typing rules for CGenFF (CHARMM General Force Field; data/typing/cgenff.typing.json).

The types follow CGenFF's own definitions (the comments of top_all36_cgenff.rtf: e.g. CG331 aliphatic C for methyl
group, HGA6 aliphatic H on fluorinated C, NG2S1 peptide nitrogen). Where DL_FIELD's templates deviate from them
(generic HA for alkane H in older templates, S_sulphate as SG301) the CGenFF definition is kept.

usage: make_cgenff_rules.py OUTDIR
"""
import json, os, sys

R = []


def add(t, p, *smarts, d=""):
    for s in smarts:
        R.append({"type": t, "smarts": s, **({"priority": p} if p else {}), **({"description": d} if d else {})})


def count(n, what, most=4):
    """exactly n neighbours matching `what` (any bond order)"""
    has = lambda k: "$(*" + "".join(f"(~{what})" for _ in range(k - 1)) + "~" + what + ")"
    return has(n) + (f";!{has(n + 1)}" if n < most else "")


# a neutral sp3 amine N (not amide, aniline, sulfonamide, enamine, imine)
AMN = "$([NX3;!$(N[#6X3]=[O,S,N]);!$(Na);!$(N[SX4]);!$(N[#6X3]=C);!$(NP);!$(N=*)])"
A5 = "$(*1aaaa1)"
A6 = "$(*1aaaaa1)"
# N in a cation: ammonium, iminium, amidinium, guanidinium (delocalised), not nitro
POSN = "[#7;!$(*~[OX1]);$([#7+]),$(*-[#6X3]=[#7X3+])]"
# a 6-ring aromatic C in a ring that carries a carbonyl or an N+ (nucleic bases, pyridinium)
CARB6 = ",".join(f"$(c1:{'a:' * k}[c;$(c=O)]:{'a:' * (3 - k)}a1)" for k in range(3)) + ",$(c1:[n+]:a:a:a:a1),$(c1:a:[n+]:a:a:a1),$(c1:a:a:[n+]:a:a1)"

# ---- hydrogen
add("HGA3", 0, "[H][CX4]", d="aliphatic H (methyl, methane)")
add("HGA1", 1, "[H][CX4H1]", d="aliphatic H (CH)")
add("HGA2", 1, "[H][CX4H2]", d="aliphatic H (CH2)")
add("HGA3", 1, "[H][CX4H3]", d="aliphatic H (CH3)")
add("HGA6", 2, "[H][CX4;" + count(1, "F") + "]", d="aliphatic H on monofluorinated C")
add("HGA7", 3, "[H][CX4](F)F", d="aliphatic H on difluorinated C")
for h, k in (("HGAAM0", 0), ("HGAAM1", 1), ("HGAAM2", 2)):
    add(h, 2, f"[H][CX4H3][NX3H{k};{AMN}]", d=f"H of a methyl on a neutral amine N with {k} H (methylamines)")
add("HGP5", 3, "[H][CX4][NX4+H0]", d="H on C next to a quaternary ammonium N")
add("HGA4", 1, "[H][CX3]=[CX3]", d="alkene H")
add("HGA5", 2, "[H][CX3H2]=[CX3]", d="terminal alkene H (=CH2)")
add("HGR61", 1, "[H]c", d="aromatic H")
add("HGR62", 2, f"[H][c;{A6};$(c:[#7,#8,#16]),$(c:c-F),$(c:c=O)]", f"[H][c;{A6};{CARB6}]", d="aromatic H, 6-ring C next to a heteroatom, C-F or in a nucleic-base ring")
add("HGR63", 3, f"[H][c;{A6};$(c:[n+]),$(c:a:[n+]),$(c:a:a:[n+])]", d="H on a pyridinium ring")
add("HGR51", 2, "[H][#6X3;r5]", d="H on sp2 C in a 5-ring")
add("HGR52", 3, "[H][#6X3;r5;$(*~[#7,#8,#16])]", "[H][CX3]=O", "[H][CX3]=[#7]", d="H on 5-ring C next to a heteroatom; aldehyde / formamide / imine H")
add("HGR53", 4, f"[H][#6X3;r5;$(*(~[#7])~[#7]);$(*~{POSN})]", d="H on C2 of imidazolium")
add("HGPAM1", 2, "[H]C#C", d="alkyne H")
add("HGP1", 0, "[H][#7]", "[H][#8]", d="polar H")
add("HGP2", 2, f"[H]{POSN}", f"[H][#7;$(*~[#6X3](~[#7X3])~[#7X3]),$(n1a[n+]aa1)]", d="polar H on N+ (ammonium, iminium, amidinium, guanidinium, imidazolium)")
add("HGP3", 1, "[H][#16]", d="thiol H")
add("HGP4", 2, "[H][NX3H2]c", d="H of an aromatic -NH2 (aniline, nucleic bases)")
for h, k in (("HGPAM1", 1), ("HGPAM2", 2), ("HGPAM3", 3)):
    add(h, 2, f"[H][NX3H{k};!R;{AMN},$(*-[CX3](=[NX2])[NX3])]", d=f"H on an acyclic neutral amine / guanidine N with {k} H")
add("HX", 3, "[H][OX1-]", d="hydroxide H")
add("H3P", 5, "[H][OX2H2]", d="water H (TIP3P)")

# ---- carbon, sp3
add("CG331", 1, "[CX4H3]", "[CX4H4]", d="methyl C / methane")
add("CG321", 1, "[CX4H2]", d="CH2")
add("CG311", 1, "[CX4H1]", d="CH")
add("CG301", 1, "[CX4H0]", d="quaternary C")
add("CG322", 3, "[CX4;" + count(1, "F") + "]", d="monofluorinated C")
add("CG312", 3, "[CX4;" + count(2, "F") + "]", d="difluorinated C")
add("CG302", 3, "[CX4](F)(F)F", d="trifluorinated C")
add("CG323", 3, "[CX4][SX1-]", d="C on a thiolate S")
add("CG314", 3, f"[CX4H1]{POSN}", d="CH next to a positive N")
add("CG324", 3, f"[CX4H2]{POSN}", d="CH2 next to a positive N")
add("CG334", 3, f"[CX4H3]{POSN}", d="CH3 next to a positive N")
for c, k in (("CG3AM0", 0), ("CG3AM1", 1), ("CG3AM2", 2)):
    add(c, 3, f"[CX4H3][NX3H{k};{AMN}]", d=f"methyl C on a neutral amine N with {k} H")
add("CG3C31", 4, "[CX4;r3]", d="cyclopropyl C")
add("CG3C41", 4, "[CX4;r4]", d="cyclobutyl C")
add("CG3C50", 4, "[CX4H0;r5]", d="5-ring quaternary sp3 C")
add("CG3C51", 4, "[CX4H1;r5]", d="5-ring sp3 CH")
add("CG3C52", 4, "[CX4H2;r5]", d="5-ring sp3 CH2")
add("CG3C53", 5, f"[CX4H1;r5]{POSN}", d="5-ring CH next to a positive N")
add("CG3C54", 5, f"[CX4H2;r5]{POSN}", d="5-ring CH2 next to a positive N")
add("CG3RC1", 6, "[CX4;x3;r3,r4,r5]", d="bridgehead of a bicyclic system with a ring of 5 or fewer")

# ---- carbon, sp2 / sp
add("CG2D1", 1, "[CX3]=[#6X3,#7X2,#7X3+]", d="alkene / imine C")
add("CG2D2", 2, "[CX3H2]=[CX3]", d="terminal alkene =CH2")
CONJ = "$(*-[#6X3]=[#6,#7,#8])"   # single bond to an sp2 C that is double-bonded (C=C, C=N, C=O)
add("CG2DC1", 3, f"[CX3;!H2;{CONJ}]=[#6X3,#7]", f"[CX3;!H2]=[#6X3;{CONJ}]", d="conjugated alkene C (with CG2DC2)")
add("CG2DC3", 4, f"[CX3H2]=[CX3;{CONJ}]", d="terminal =CH2 of a conjugated alkene")
add("CG2D1O", 3, "[CX3;!H2;$(*-[#7,#8,#16])]=[#6X3]", d="double-bond C next to a heteroatom (enol ether, enamine)")
add("CG2O1", 2, "[CX3](=O)[NX3]", d="amide C")
add("CG2O2", 2, "[CX3](=O)[OX2]", d="carboxylic acid / ester C")
add("CG2O5", 2, "[CX3](=O)([#6])[#6]", d="ketone C")
add("CG2O4", 3, "[CX3H1](=O)", "[CX3H2]=O", d="aldehyde C")
add("CG2O3", 3, "[CX3](=O)[OX1-]", d="carboxylate C")
add("CG2O6", 4, "[CX3](=[O,S])([#7,#8,#16])[#7,#8,#16]", d="urea / carbonate / carbamate C (and thio analogues)")
add("CG2O7", 3, "[CX2](=O)=O", d="CO2 C")
add("CG2N2", 4, "[CX3;$(*=[NX3+])](~[NX3])[#6,#1]", d="amidinium C")
add("CG2N1", 5, "[CX3](~[NX3])(~[NX3])~[NX3]", "[CX3](~[NX3])(~[NX3])=[NX2]", d="guanidine / guanidinium C")
add("CG2R61", 0, "[c]", d="aromatic C, 6-ring")
add("CG2R66", 2, "[c;r6]F", d="aromatic C bonded to F")
add("CG2R67", 2, "[c;r6](-!@[a;r6])", d="biphenyl / bipyridyl bridge C")
add("CG2R62", 2, f"[c;{A6};{CARB6}]", d="6-ring aromatic C in a ring with a carbonyl or N+ (nucleic bases, pyridinium)")
add("CG2R64", 3, "[c;r6](:n):n", "[c;r6](:n)-[NX3]", d="aromatic amidine C (between two N)")
add("CG2R63", 4, "[c;r6]=O", d="6-ring aromatic carbonyl C (nucleic bases)")
# five-membered rings: CGenFF's ring types win over the chain types (alkene, lactam carbonyl)
add("CG2R51", 6, "[#6X3;r5;!$(*=[O,S])]", d="sp2 C in a 5-ring")
add("CG2R52", 7, "[#6X3;r5]=[#7]", d="5-ring C double-bonded to N")
add("CG2R53", 8, "[#6X3;r5](=[#7,#8])~[#7,#8,#16]", f"[c;{A5}](~[#7,#8,#16])~[#7,#8,#16]", d="5-ring C double-bonded to N (or =O) next to another heteroatom")
add("CG2RC0", 9, "[#6X3;r5;r6;x3]", d="C at a 5/6 ring fusion")
add("CG1T1", 1, "[CX2]#C", d="alkyne C")
add("CG1T2", 2, "[CX2H1]#C", d="terminal alkyne C")
add("CG1N1", 2, "[CX2]#N", d="nitrile C")

# ---- nitrogen
add("NG301", 1, "[NX3H0]", d="tertiary amine N")
add("NG311", 1, "[NX3H1]", d="secondary amine N")
add("NG321", 1, "[NX3H2]", d="primary amine N")
add("NG331", 1, "[NX3H3]", d="ammonia")
add("NG2S0", 3, "[NX3H0][CX3]=O", d="tertiary amide N")
add("NG2S1", 3, "[NX3H1][CX3]=O", d="secondary amide N (peptide)")
add("NG2S2", 3, "[NX3H2][CX3]=O", d="primary amide N")
add("NG2S3", 2, "[NX3H2]c", "[NX3]P", d="aniline NH2, phosphoramide N")
add("NG2R60", 1, f"[nX2;{A6}]", d="pyridine N")
add("NG2R62", 2, f"[nX2;{A6};$(n:[n,o,s]),$(n:a:[n,o,s])]", d="6-ring N with a ring heteroatom ortho or meta (pyrimidine)")
add("NG2R61", 1, f"[nX3;{A6}]", d="6-ring aromatic N with a substituent (pyridone)")
add("NG2R67", 3, "[n;r6](-!@[a;r6])", d="6-ring N substituted with a 6-ring (N-phenyl pyridinone)")
add("NG2R50", 1, "[#7X2;r5]", d="5-ring N, double-bonded")
add("NG2R51", 1, "[#7X3;r5;{AR2},{AR3}]", d="5-ring N, single-bonded, planar ring (pyrrole type)")
add("NG3C51", 2, f"[NX3H1;r5;{AMN}]", d="secondary sp3 amine in a 5-ring")
add("NG3N1", 3, f"[NX3;{AMN}][NX3]", d="hydrazine N")
add("NG2R53", 4, "[NX3;r5][CX3;r5]=O", d="amide N in a non-planar 5-ring (lactam)")
add("NG2R52", 3, f"[#7X3;r5;$([#7+]),$(*-[#6X3;r5]=[#7X3+;r5])]", d="protonated imine / amidinium N in a 5-ring (imidazolium)")
add("NG2RC0", 3, "[nX3;R2]", d="N at a ring fusion (indolizine)")
add("NG2D1", 1, "[NX2]=[#6]", d="imine N")
add("NG2P1", 4, "[NX3]~[CX3](~[NX3])~[NX3]", "[NX3;$([NX3+]=[CX3]),$(*-[CX3]=[NX3+])]", d="guanidinium / amidinium / protonated imine N")
add("NG2O1", 3, "[NX3+](=O)[O-]", d="nitro N")
add("NG1T1", 2, "[NX1]#C", d="nitrile N")
for n_, k in (("NG3P0", 0), ("NG3P1", 1), ("NG3P2", 2), ("NG3P3", 3)):
    add(n_, 3, f"[NX4+H{k}]", d=f"ammonium N with {k} H")

# ---- oxygen
add("OG301", 0, "[OX2]", d="ether O")
add("OG311", 1, "[OX2H1]", d="hydroxyl O")
add("OG302", 2, "[OX2;!H1][CX3]=O", d="ester O")
add("OG303", 3, "[OX2;!H1][P,S]", d="phosphate / sulfate ester O")
add("OG304", 4, "[OX2](P)P", d="pyrophosphate bridging O")
add("OG3C31", 3, "[OX2;r3]", d="epoxide O")
add("OG3C51", 3, "[OX2;r5]", d="5-ring ether O")
add("OG3C61", 3, "[OX2;r6;H0]", d="6-ring ether O")
add("OG2R50", 3, "[o]", d="furan O")
add("OG2D1", 1, "[OX1]=[#6]", d="carbonyl O (amide, ester, acid, aldehyde)")
add("OG2D3", 2, "[OX1]=[CX3]([#6])[#6]", d="ketone O")
add("OG2D4", 2, "[OX1]=c", d="carbonyl O on an aromatic ring")
add("OG2D2", 3, "[OX1-][CX3]=O", "[OX1]=[CX3][OX1-]", d="carboxylate / carbonate O")
add("OG2D5", 3, "[OX1]=[CX2]=O", d="CO2 O")
add("OG2P1", 2, "[OX1]~[P,S]", d="=O / O- on P or S")
add("OG2N1", 2, "[OX1]~[N+]", d="nitro O")
add("OG312", 2, "[OX1-][#6]", d="alkoxide / phenoxide O")
add("OX", 3, "[OX1-][H]", d="hydroxide O")
add("OT", 5, "[OX2H2]", d="water O (TIP3P)")

# ---- sulfur, phosphorus, halogens, ions
add("SG311", 0, "[#16X2]", d="sulfide / thiol S")
add("SG301", 2, "[SX2]S", d="disulfide S")
add("SG302", 2, "[SX1-]", d="thiolate S")
add("SG2R50", 2, "[s]", d="thiophene S")
add("SG2D1", 2, "[SX1]=[#6]", d="thiocarbonyl S")
add("SG3O3", 1, "[SX3]", d="sulfoxide S")
add("SG3O2", 1, "[SX4]", d="sulfone / sulfonamide S")
add("SG3O1", 2, "[SX4](~[OX1])(~[OX1])~[OX1]", d="sulfate / sulfonate S")
add("PG0", 0, "[PX4]", d="neutral phosphate P")
add("PG1", 1, "[PX4;" + count(2, "[OX1]") + "]", d="phosphate P, -1")
add("PG2", 2, "[PX4;" + count(3, "[OX1]") + "]", d="phosphate P, -2")
for n in (1, 2, 3):
    add(f"FGA{n}", 1, f"[F][CX4;{count(n, 'F')}]", d=f"aliphatic F ({n} on the C)")
add("FGR1", 1, "[F]c", d="aromatic F")
add("CLGA1", 0, "[Cl]", d="aliphatic Cl")
add("CLGA3", 1, "[Cl][CX4](Cl)Cl", d="Cl of CCl3")
add("CLGR1", 1, "[Cl]c", d="aromatic Cl")
for n in (1, 2, 3):
    add(f"BRGA{n}", 1 if n == 1 else 2, f"[Br][CX4;{count(n, 'Br')}]", d=f"aliphatic Br ({n} on the C)")
add("BRGR1", 1, "[Br]c", d="aromatic Br")
add("IGR1", 0, "[I]", d="iodine (IGR1 for every I)")
add("SOD", 2, "[Na]"); add("POT", 2, "[K]"); add("CLA", 2, "[Cl-;X0]")

# Bonded terms CGenFF's file lacks, by analogy (Vanommeslaeghe & MacKerell, J. Chem. Inf. Model. 52, 3155 (2012) assign
# them by analogy with penalties): stand-ins of the same element and hybridisation, nearest first. Every term found this
# way is listed as estimated.
ANALOGIES = {
    # sp3 carbon by hydrogen count; ring, fluorinated, charged-neighbour and amine-side variants to the plain type
    "CG331": ["CG321", "CG311"], "CG321": ["CG311", "CG331", "CG301"], "CG311": ["CG321", "CG301", "CG331"], "CG301": ["CG311", "CG321"],
    "CG3C52": ["CG321"], "CG3C51": ["CG311"], "CG3C50": ["CG301"], "CG3C53": ["CG321"], "CG3C54": ["CG321"], "CG3C41": ["CG311", "CG321"],
    "CG3C31": ["CG311", "CG321"], "CG3RC1": ["CG311"], "CG322": ["CG321"], "CG312": ["CG311"], "CG302": ["CG301"],
    "CG334": ["CG331"], "CG324": ["CG321"], "CG314": ["CG311"], "CG323": ["CG321"], "CG3AM2": ["CG331"], "CG3AM1": ["CG331"], "CG3AM0": ["CG331"],
    # sp2 carbon: alkenes, conjugated alkenes, aromatic ring members, carbonyl carbons of the same group
    "CG2D1": ["CG2D2", "CG2DC1"], "CG2D2": ["CG2D1", "CG2DC3"], "CG2DC1": ["CG2D1"], "CG2DC2": ["CG2D1"], "CG2DC3": ["CG2D2"],
    "CG2D1O": ["CG2D1"], "CG2D2O": ["CG2D2"], "CG25C1": ["CG2DC1", "CG2D1"], "CG25C2": ["CG2DC3", "CG2D2"], "CG251O": ["CG2D1"], "CG252O": ["CG2D2"],
    "CG2R62": ["CG2R61"], "CG2R63": ["CG2R61"], "CG2R64": ["CG2R61"], "CG2R66": ["CG2R61"], "CG2R67": ["CG2R61"], "CG2R71": ["CG2R61"],
    "CG2RC0": ["CG2R61"], "CG2R52": ["CG2R51"], "CG2R53": ["CG2R51"], "CG2R57": ["CG2R51"],
    "CG2O4": ["CG2O5"], "CG2O5": ["CG2O4"], "CG2O6": ["CG2O2"],
    # oxygen: ring and other ethers to the plain ether; nitrogen: amines by hydrogen count, amides by hydrogen count
    "OG3C51": ["OG301"], "OG3C61": ["OG301"], "OG3R60": ["OG301"], "OG3C31": ["OG301"], "OG303": ["OG301"], "OG304": ["OG301"], "OG302": ["OG301"],
    "NG321": ["NG311", "NG331"], "NG311": ["NG321", "NG301"], "NG301": ["NG311"], "NG331": ["NG321"],
    "NG2S2": ["NG2S1"], "NG2S0": ["NG2S1"], "NG2S3": ["NG2S1"],
    # aliphatic and aromatic hydrogens by the carbon's hydrogen count
    "HGA1": ["HGA2", "HGA3"], "HGA2": ["HGA1", "HGA3"], "HGA3": ["HGA2", "HGA1"], "HGA6": ["HGA2"], "HGA7": ["HGA1"],
    "HGR62": ["HGR61"], "HGR63": ["HGR61"], "HGR52": ["HGR51"], "HGR53": ["HGR51"], "HGR71": ["HGR61"],
}

out = sys.argv[1]
json.dump({"format": "caps-typing", "version": 1, "forcefield": "CGenFF (CHARMM36)",
           "description": "CGenFF atom types from CGenFF's own definitions (top_all36_cgenff.rtf); conjugated alkenes pair "
                          "CG2DC1 / CG2DC2 (same type across a double bond); water is TIP3P (OT / H3P).",
           "pair_mode": "double_same", "pairs": [["CG2DC2", "CG2DC1"]], "rules": R,
           "analogies": ANALOGIES,
           "analogy_source": "CGenFF analogues: the same element and hybridisation, nearest in hydrogen count, ring or substituent (the CGenFF program "
                             "assigns missing bonded parameters by analogy with penalties; this table is CAPS's simpler version)"},
          open(os.path.join(out, "cgenff.typing.json"), "w"), indent=1)
print(len(R), "rules")
