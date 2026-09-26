#!/usr/bin/env python3
"""Write CAPS typing rules for OPLS-AA by type number (data/typing/oplsaa2024-moltemplate.typing.json and the 2008 set).

OPLS-AA (Jorgensen, Maxwell & Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996), and its later tables) gives every
chemical environment its own type number, each with its partial charge; the charges of a functional group and the
carbons and hydrogens around it add up to zero (a CH2 of an ether carries +0.20 against the oxygen's -0.40, an alkoxy
CH2 of an ester +0.25 against the ester's -0.25 ...). The bonded parameters come from each number's class (CT, CA, OS ...).
So the rules below name numbers, group by group, each with the carbons and hydrogens that go with it; nothing is
typed where the table has no number for the environment (a siloxane, an aryl carbonate, a CCl2 carbon): CAPS
then reports those atoms instead of guessing charges that would not balance.

Priorities: 0 the plain element and hybridisation (alkane, alkene, benzene C and H); higher values for functional groups
and the atoms next to them; among matching rules the highest wins.

Only for the 2024 table: the 2008 (BOSS 4.8) table numbers the same chemistry differently.

usage: make_oplsaa_rules.py OUTDIR
"""
import json, os, sys

R = []


def add(num, prio, smarts, desc):
    for s in ([smarts] if isinstance(smarts, str) else smarts):
        R.append({"type": str(num), "smarts": s, "priority": prio, "description": desc})


def by_h(nums, prio, env, desc):
    """sp3 carbon types by hydrogen count (H3, H2, H1, H0), each with the environment `env` as a recursive SMARTS."""
    for h, n in zip((3, 2, 1, 0), nums):
        if n:
            add(n, prio, f"[CX4H{h};{env}]", f"{desc} ({['C', 'CH', 'CH2', 'CH3'][h]})")


# ---------------------------------------------------------------- hydrocarbons
# the plain hydrocarbon types only where every neighbour is C or H: a carbon bearing any other atom must match a
# rule for its group, or stay untyped (its charge would not balance)
HC_ONLY = "!$(*~[!#6;!#1])"
add(135, 0, f"[CX4H3;{HC_ONLY}]", "alkane C (CH3)")
add(136, 0, f"[CX4H2;{HC_ONLY}]", "alkane C (CH2)")
add(137, 0, f"[CX4H1;{HC_ONLY}]", "alkane C (CH)")
add(139, 0, f"[CX4H0;{HC_ONLY}]", "alkane C (C)")
add(140, 0, "[H][CX4]", "H on sp3 C")
add(141, 1, f"[CX3H0;{HC_ONLY}]=[CX3]", "alkene C (R2C=)")
add(142, 1, f"[CX3H1;{HC_ONLY}]=[CX3]", "alkene C (RHC=)")
add(143, 1, f"[CX3H2;{HC_ONLY}]=[CX3]", "alkene C (H2C=)")
add(144, 1, "[H][CX3]=[CX3]", "alkene H")
add(145, 1, f"[c;{HC_ONLY}]", "aromatic C")
add(146, 1, "[H]c", "aromatic H")
add(147, 2, "[c;R2;$(c(:c)(:c):c)]", "fused aromatic C (naphthalene C9)")
add(519, 3, "[c;$(c-!@c)]", "biphenyl junction C")
by_h((148, 149, 515, 516), 2, f"$(*c);{HC_ONLY}", "benzylic C")

# ---------------------------------------------------------------- alcohols and phenols
ALC = "$(*[OX2H1])"
add(154, 5, "[OX2H1][CX4]", "alcohol O")
add(160, 6, "[OX2H1][CX4H0]", "tertiary alcohol O")
add(155, 5, "[H][OX2H1][CX4]", "alcohol H")
by_h((157, 157, 158, 159), 5, ALC, "alcohol C")
add(167, 5, "[OX2H1]c", "phenol O")
add(168, 5, "[H][OX2H1]c", "phenol H")
add(166, 5, "c[OX2H1]", "phenol C(OH)")

# ---------------------------------------------------------------- ethers (neither ester, carbonate nor acetal)
ETHO = "[OX2;!$(O[#1]);!$(O[#6]=[O,S,N]);!$(O[Si]);!$(O~[#6](~[OX2])~[OX2]);!$(O[CX4]([OX2])[OX2]);!r3]"
add(180, 5, f"{ETHO[:-1]};$(O([CX4])[CX4])]", "dialkyl ether O")
add(179, 5, f"{ETHO[:-1]};$(O(c)[CX4])]", "alkyl aryl ether O (anisole)")
add(177, 5, f"{ETHO[:-1]};$(O(c)c)]", "diaryl ether O")
by_h((181, 182, 183, 184), 5, f"$(*{ETHO[:-1]};$(O([CX4])[CX4,c])])", "ether alpha C")
add(185, 5, f"[H][CX4]{ETHO[:-1]};$(O([CX4])[CX4,c])]", "H on ether alpha C")
add(199, 5, f"c{ETHO[:-1]};$(O(c)[CX4,c])]", "aryl ether C(O)")

# ---------------------------------------------------------------- carbonyls: aldehydes, ketones
add(277, 5, "[CX3H1](=[OX1])[#6]", "aldehyde C")
add(278, 5, "[OX1]=[CX3H1][#6]", "aldehyde O")
add(279, 5, "[H][CX3]=[OX1]", "aldehyde H")
add(280, 5, "[CX3H0](=[OX1])([#6])[#6]", "ketone C")
add(281, 5, "[OX1]=[CX3H0]([#6])[#6]", "ketone O")
add(282, 4, "[H][CX4][CX3](=[OX1])[#6,#1]", "H on C alpha to a ketone / aldehyde")
# on an aromatic ring the ipso C keeps -0.115: the carbonyl C carries +0.115 more
add(232, 6, "[CX3H1](=[OX1])c", "aryl aldehyde C (benzaldehyde)")
add(233, 6, "[CX3H0](=[OX1])(c)[CX4]", "alkyl aryl ketone C (acetophenone)")
add(231, 7, "[CX3H0](=[OX1])(c)c", "diaryl ketone C (benzophenone)")

# ---------------------------------------------------------------- esters (alkyl or aryl acyl; alkyl or aryl alkoxy)
ESTC = "[CX3](=[OX1])[OX2][#6]"
add(465, 6, "[CX3;$(*(=[OX1])[OX2][#6]);$(*[#6,#1]);!$(*[OX2][#6][OX2]);!$(*(O)(O)=O)]", "ester C=O (alkyl acyl)")
add(471, 7, "[CX3;$(*(=[OX1])[OX2][#6]);$(*c)]", "ester C=O (aryl acyl, benzoate)")
add(466, 6, "[OX1]=[CX3;$(*[OX2][#6]);!$(*([OX2])[OX2])]", "ester =O")
add(467, 6, "[OX2;$(O[CX3](=[OX1])[#6,#1]);$(O[CX4])]", "ester -O- (alkoxy)")
add(473, 7, "[OX2;$(O[CX3](=[OX1])[#6,#1]);$(Oc)]", "ester -O- (phenyl ester)")
add(472, 7, "c[OX2][CX3](=[OX1])[#6,#1]", "phenyl ester C(ipso)")
by_h((468, 490, 491, 492), 6, "$(*[OX2][CX3](=[OX1])[#6,#1])", "ester alkoxy C")
add(469, 6, "[H][CX4][OX2][CX3](=[OX1])[#6,#1]", "H on ester alkoxy C")

# ---------------------------------------------------------------- carboxylic acids
add(267, 7, "[CX3](=[OX1])[OX2H1]", "carboxylic acid C")
add(269, 7, "[OX1]=[CX3][OX2H1]", "carboxylic acid =O")
add(268, 7, "[OX2H1][CX3]=[OX1]", "carboxylic acid -OH")
add(270, 7, "[H][OX2H1][CX3]=[OX1]", "carboxylic acid H")

# ---------------------------------------------------------------- carbonates (dialkyl; aryl carbonates have no OPLS-AA type)
add(772, 7, "[CX3](=[OX1])([OX2][CX4])[OX2][CX4]", "carbonate C=O")
add(771, 7, "[OX1]=[CX3]([OX2][CX4])[OX2][CX4]", "carbonate =O")
add(773, 7, "[OX2;$(O[CX3](=[OX1])[OX2][CX4]);$(O[CX4])]", "carbonate -O-")
by_h((789, 774, 775, None), 7, "$(*[OX2][CX3](=[OX1])[OX2][CX4])", "carbonate alkoxy C")
add(777, 7, "[H][CX4][OX2][CX3](=[OX1])[OX2][CX4]", "H on carbonate alkoxy C")

# ---------------------------------------------------------------- amides (urea, imide and carbamate not covered here)
AMC = "[CX3;$(*(=[OX1])[NX3]);!$(*([NX3])[NX3]);!$(*[OX2]);!$(*[NX3][CX3]=[OX1])]"
add(235, 6, AMC, "amide C=O")
add(236, 6, f"[OX1]={AMC}", "amide =O")
add(237, 6, f"[NX3H2]{AMC}", "primary amide N")
add(238, 6, f"[NX3H1;$(*{AMC})][#6]", "secondary amide N")
add(239, 6, f"[NX3H0;$(*{AMC})]([#6])[#6]", "tertiary amide N")
add(240, 6, f"[H][NX3H2]{AMC}", "H on primary amide N")
add(241, 6, f"[H][NX3H1]{AMC}", "H on secondary amide N")
add(234, 7, f"[CX3;$(*(=[OX1])[NX3]);$(*c);!$(*([NX3])[NX3]);!$(*[OX2])]", "aryl amide C=O (benzamide)")
# N-aryl amides: N and its ring C balance the amide H (N-phenylacetamide, N-methyl-N-phenylacetamide)
add(265, 7, f"[NX3H1;$(*{AMC})]c", "secondary N-aryl amide N")
add(266, 7, f"c[NX3H1]{AMC}", "ring C on a secondary amide N")
add(987, 7, f"[NX3H0;$(*{AMC})](c)[CX4]", "tertiary N-alkyl-N-aryl amide N")
add(988, 7, f"c[NX3H0;$(*{AMC})][CX4]", "ring C on a tertiary N-alkyl amide N")
by_h((242, 244, 229, 230), 6, f"$(*[NX3H1]{AMC})", "C on secondary amide N")
by_h((243, 245, 246, None), 6, f"$(*[NX3H0]{AMC})", "C on tertiary amide N")

# ---------------------------------------------------------------- amines and anilines
AMN = "!$(*[CX3]=[OX1,SX1,NX2]);!$(*[SX4]);!$(*[NX2,NX1])"
add(900, 5, f"[NX3H2;{AMN}][CX4,c]", "primary amine N")
add(901, 5, f"[NX3H1;{AMN}]([CX4])[CX4]", "secondary amine N")
add(902, 5, f"[NX3H0;{AMN}]([CX4])([CX4])[CX4]", "tertiary amine N")
add(909, 5, f"[H][NX3H2;{AMN}]", "H on primary amine N")
add(910, 5, f"[H][NX3H1;{AMN}]", "H on secondary amine N")
by_h((903, 906, 912, 913), 5, f"$(*[NX3H2;{AMN}])", "C on primary amine N")
by_h((904, 907, 914, None), 5, f"$(*[NX3H1;{AMN}])", "C on secondary amine N")
by_h((905, 908, 915, None), 5, f"$(*[NX3H0;{AMN}])", "C on tertiary amine N")
add(911, 5, f"[H][CX4][NX3;{AMN}]", "H on C bonded to an amine N")
add(916, 6, f"c[NX3H2;{AMN}]", "aniline C(NH2)")
add(917, 6, f"c[NX3H1;{AMN}]", "N-alkylaniline C(N)")
add(918, 6, f"c[NX3H0;{AMN}]", "N,N-dialkylaniline C(N)")

# ---------------------------------------------------------------- nitriles
add(753, 6, "[NX1]#[CX2][CX4]", "nitrile N")
add(754, 6, "[CX2](#[NX1])[CX4]", "nitrile C")
by_h((755, 756, 757, 758), 6, "$(*[CX2]#[NX1])", "C alpha to a nitrile")
add(759, 6, "[H][CX4][CX2]#[NX1]", "H alpha to a nitrile")
add(262, 6, "[NX1]#[CX2]c", "benzonitrile N")
add(261, 6, "[CX2](#[NX1])c", "benzonitrile C")
add(260, 6, "c[CX2]#[NX1]", "benzonitrile C(ipso)")

# ---------------------------------------------------------------- halogens (one halogen per carbon; CX2 / CX3 groups only when perfluorinated)
ONEX = lambda X: f"$(*[{X}]);!$(*([{X}])[{X}]);!$(*[F,Cl,Br,I;!{X}])"
add(151, 5, "[Cl][CX4;!$(*(Cl)[Cl,F,Br,I])]", "alkyl chloride Cl (one halogen on the C)")
by_h((None, 152, 973, 974), 5, ONEX("Cl"), "alkyl chloride C")
add(153, 5, "[H][CX4H2][Cl]", "H on RCH2Cl")
add(972, 5, "[H][CX4H1][Cl]", "H on R2CHCl")
add(975, 5, "[Br][CX4;!$(*(Br)[Cl,F,Br,I])]", "alkyl bromide Br (one halogen on the C)")
by_h((None, 976, 978, None), 5, ONEX("Br"), "alkyl bromide C")
add(956, 5, "[F][CX4;$(*[F]);!$(*(F)F)]", "monofluoroalkyl F")
by_h((None, 957, 959, 960), 5, ONEX("F"), "monofluoroalkyl C")
add(958, 5, "[H][CX4;$(*[F]);!$(*(F)F)]", "H on monofluoroalkyl C")
add(965, 6, "[F][CX4;$(*(F)F)]", "F of a CF2 / CF3 group")
add(962, 6, "[CX4;$(*(F)F);!$(*(F)(F)F);!$(*[H])]", "CF2")
add(961, 6, "[CX4;$(*(F)(F)F);!$(*(F)(F)(F)F)]", "CF3")
add(264, 6, "[Cl]c", "aryl chloride Cl")
add(263, 6, "c[Cl]", "aryl chloride C")
add(719, 6, "[F]c", "aryl fluoride F")
add(718, 6, "c[F]", "aryl fluoride C")
add(730, 6, "[Br]c", "aryl bromide Br")
add(729, 6, "c[Br]", "aryl bromide C")
add(226, 6, "[Cl][CX3;!$(*(Cl)Cl)]=[CX3]", "chloroalkene Cl, one Cl on the C (ClHC=, R-CCl=)")
add(227, 6, "[CX3H1]([Cl])=[CX3]", "chloroalkene C (ClHC=)")
# R-CCl= (chloroprene): an alkyl in place of 227's H. On an alkene C that swap adds +0.115 e (RHC= 86 against H2C= 87,
# with the H's +0.115 gone), so 227's +0.005 becomes +0.12, balancing the Cl's -0.12; the table's CM number with that
# charge is 399, so its number is used (same class, same bonded and LJ terms as 227)
add(399, 7, "[CX3H0;!$(*(Cl)Cl)]([Cl])([CX4])=[CX3]", "chloroalkene C with an alkyl group (R-CCl=, chloroprene): CM, +0.12 e")
add(398, 7, "[Cl][CX3;$(*(Cl)Cl)]=[CX3]", "chloroalkene Cl (Cl2C=)")
add(399, 7, "[CX3;$(*(Cl)Cl)]=[CX3]", "chloroalkene C (Cl2C=)")

# ---------------------------------------------------------------- sulfur
add(493, 6, "[SX4](=[OX1])(=[OX1])([#6])[#6]", "sulfone S")
add(494, 6, "[OX1]=[SX4](=[OX1])([#6])[#6]", "sulfone O")
add(488, 7, "c[SX4](=[OX1])(=[OX1])", "ring C on a sulfonyl S (charge 0: the SO2 group is neutral)")
add(202, 5, "[SX2]([CX4])[CX4]", "sulfide S")
add(203, 6, "[SX2]([#6])[SX2]", "disulfide S")
add(200, 5, "[SX2H1][CX4]", "thiol S")
add(204, 5, "[H][SX2H1]", "thiol H")
by_h((209, 210, 211, 212), 5, "$(*[SX2]([CX4])[CX4])", "C on a sulfide S")
by_h((213, 214, 215, 216), 6, "$(*[SX2][SX2])", "C on a disulfide S")
by_h((None, 206, 207, 208), 5, "$(*[SX2H1])", "C on a thiol S")

# ---------------------------------------------------------------- epoxides
add(1025, 6, "[OX2;r3]([CX4])[CX4]", "epoxide O")
by_h((None, 1026, 1027, 1028), 6, "$(*1[OX2][CX4]1)", "epoxide C")
add(1029, 6, "[H][CX4;$(*1[OX2][CX4]1)]", "H on epoxide C")


# ---------------------------------------------------------------- small molecules and special cases
add(138, 3, "[CX4H4]", "methane C")
add(156, 6, "[H][CX4H3][OX2H1]", "H on the C of methanol")
add(9999, 8, "[OX2H2]", "water O (TIP3P, the table's default water)")
add(9998, 8, "[H][OX2H2]", "water H (TIP3P)")
add(217, 6, "[CX4H3][SX2H1]", "C of methanethiol")
add(734, 6, "[SX2H1]c", "thiophenol S")
add(735, 6, "c[SX2H1]", "thiophenol C(S)")
# nitro groups: N and O, the C carrying them, and its H
add(760, 6, "[NX3+,NX3](=[OX1])(~[OX1])[CX4]", "nitro N (aliphatic)")
add(767, 7, "[NX3+,NX3](=[OX1])(~[OX1])c", "nitro N (aromatic)")
add(761, 6, "[OX1]~[NX3](~[OX1])[#6]", "nitro O")
by_h((762, 764, 765, 766), 6, "$(*[NX3](~[OX1])~[OX1])", "C carrying a nitro group")
add(763, 6, "[H][CX4][NX3](~[OX1])~[OX1]", "H alpha to a nitro group")
add(768, 7, "c[NX3](~[OX1])~[OX1]", "C(NO2) of a nitroarene")

# silicon: silanes, silanols and silyl ethers with at most one O on Si (siloxanes, Si(OR)2 / Si(OR)3 have no number)
SI0 = "!$(*~[O,N,S,Cl,F])"
for h, n in zip((0, 1, 2, 3), (1060, 1061, 1062, 1063)):
    add(n, 6, f"[SiX4H{h};{SI0}]", f"Si in a silane (SiH{h})" if h else "Si in a tetraalkylsilane")
for h, n in zip((0, 1, 2), (1070, 1071, 1072)):
    add(n, 6, f"[SiX4H{h};$(*[OX2H1]);!$(*(O)O)]", f"Si of a silanol (SiH{h})" if h else "Si of a silanol")
for h, n in zip((0, 1, 2), (1075, 1076, 1077)):
    add(n, 6, f"[SiX4H{h};$(*[OX2][#6]);!$(*(O)O)]", f"Si of a silyl ether (SiH{h})" if h else "Si of a silyl ether")
add(1064, 6, "[H][SiX4;!$(*(O)O)]", "H on Si")
add(1073, 6, "[OX2H1][SiX4;!$(*(O)O)]", "silanol O")
add(1074, 6, "[H][OX2H1][SiX4;!$(*(O)O)]", "silanol H")
add(1078, 6, "[OX2]([SiX4;!$(*(O)O)])[CX4]", "silyl ether O")
by_h((1065, 1066, 1067, 1068), 6, "$(*[SiX4;!$(*(O)O)])", "C on Si")
add(1069, 7, "c[SiX4;!$(*(O)O)]", "ring C on Si (phenylsilane)")
by_h((181, 182, 183, 184), 6, "$(*[OX2][SiX4;!$(*(O)O)])", "C on a silyl ether O")
add(185, 6, "[H][CX4][OX2][SiX4;!$(*(O)O)]", "H on C of a silyl ether O")

# five-membered heteroaromatics: OPLS charges of the parent ring; an alkyl at C2 takes a carbon type whose charge
# stands in for the H it replaces (2-methylfuran, 2-ethylpyrrole)
add(566, 7, "[o;r5;$(o1cccc1)]", "furan O")
add(567, 7, "[c;r5;$(c1occc1)]", "furan C2")
add(568, 7, "[c;r5;$(c1cocc1)]", "furan C3")
add(569, 7, "[H][c;r5;$(c1occc1)]", "furan H2")
add(570, 7, "[H][c;r5;$(c1cocc1)]", "furan H3")
add(680, 7, "[CX4H3][c;r5;$(c1occc1)]", "CH3 on furan C2")
add(681, 7, "[CX4H2][c;r5;$(c1occc1)]", "CH2 on furan C2")
add(542, 7, "[nH1;r5;$(n1cccc1)]", "pyrrole N")
add(543, 7, "[c;r5;$(c1[nH1]ccc1)]", "pyrrole C2")
add(544, 7, "[c;r5;$(c1c[nH1]cc1)]", "pyrrole C3")
add(545, 7, "[H][nH1;r5;$(n1cccc1)]", "pyrrole H1")
add(546, 7, "[H][c;r5;$(c1[nH1]ccc1)]", "pyrrole H2")
add(547, 7, "[H][c;r5;$(c1c[nH1]cc1)]", "pyrrole H3")
add(678, 7, "[CX4H3][c;r5;$(c1[nH1]ccc1)]", "CH3 on pyrrole C2")
add(679, 7, "[CX4H2][c;r5;$(c1[nH1]ccc1)]", "CH2 on pyrrole C2")


def rules_for(ff_json):
    """the rules whose numbers the force field has (the 2008 table lacks some of the 2024 numbers)"""
    ff = json.load(open(ff_json))
    have = {t["name"].split("_b")[0] for t in ff["atom_types"]}
    return [r for r in R if r["type"] in have], sorted({r["type"] for r in R if r["type"] not in have}, key=int)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "data/typing"
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    # the 2024 table's numbers only: BOSS 4.8's 2008 table numbers the same chemistry differently (its 135 is an acetal
    # carbon, 2024's an alkane CH3), so these rules would type it wrongly; it needs rules of its own
    for ff, label in (("oplsaa2024-moltemplate", "OPLS-AA (2024 parameter file)"),):
        rules, missing = rules_for(os.path.join(root, "data", "forcefields", ff + ".json"))
        doc = {"format": "caps-typing", "version": 1, "forcefield": label,
               "description": "CAPS rules giving OPLS-AA type numbers by chemical group, with the carbons and hydrogens that balance each group's "
                              "charge (Jorgensen, Maxwell & Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)). Environments the table has no "
                              "number for (siloxanes, aryl carbonates, CCl2 carbons, ureas, imides, carbamates) are left untyped and reported. "
                              "Among matching rules the highest priority wins.",
               "rules": rules}
        p = os.path.join(out, ff + ".typing.json")
        json.dump(doc, open(p, "w"), ensure_ascii=False, indent=1)
        print(f"{p}: {len(rules)} rules" + (f" (numbers not in this table: {', '.join(missing)})" if missing else ""))
