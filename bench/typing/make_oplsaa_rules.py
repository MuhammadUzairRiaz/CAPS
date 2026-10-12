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

For the 2024 table; the 2008 (BOSS 4.8) numbering gets these rules mapped by make_oplsaa2008_rules.py.

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
# small rings (OPLS-2020): cyclopropane and cyclobutane carbons have their own CY types; their H stay 140
by_h((None, 711, 712, 713), 1, f"r3;{HC_ONLY}", "cyclopropane C")
by_h((None, 714, 715, 716), 1, f"r4;{HC_ONLY}", "cyclobutane C")

# ---------------------------------------------------------------- alcohols and phenols
ALC = "$(*[OX2H1])"
add(154, 5, "[OX2H1][CX4]", "alcohol O")
add(155, 5, "[H][OX2H1][CX4]", "alcohol H")
by_h((157, 157, 158, 159), 5, ALC, "alcohol C")
# glycerol-like triols (three OH on consecutive carbons): the table's triol types
TRI = ["$(*([OX2H1])[CX4]([OX2H1])[CX4][OX2H1])", "$(*([OX2H1])([CX4][OX2H1])[CX4][OX2H1])"]
for env in TRI:
    by_h((None, 173, 174, 175), 8, env, "triol C")
    add(176, 8, f"[H][CX4;{env}]", "H on a triol C")
    add(171, 8, f"[OX2H1][CX4;{env}]", "triol O")
    add(172, 8, f"[H][OX2H1][CX4;{env}]", "triol H(O)")
by_h((None, 218, 219, 220), 7, "$(*c);$(*[OX2H1])", "benzyl alcohol C")
add(221, 7, "c[CX4][OX2H1]", "ring C of a benzyl alcohol")
add(167, 5, "[OX2H1]c", "phenol O")
add(168, 5, "[H][OX2H1]c", "phenol H")
add(166, 5, "c[OX2H1]", "phenol C(OH)")

# ---------------------------------------------------------------- ethers (neither ester, carbonate nor acetal)
ETHO = "[OX2;!$(O[#1]);!$(O[#6]=[O,S,N]);!$(O[Si]);!$(O~[#6](~[OX2])~[OX2]);!$(O[CX4]([OX2])[OX2]);!r3]"
add(180, 5, f"{ETHO[:-1]};$(O([CX4])[CX4])]", "dialkyl ether O")
add(179, 5, f"{ETHO[:-1]};$(O(c)[CX4])]", "alkyl aryl ether O (anisole)")
add(473, 5, f"{ETHO[:-1]};$(O(c)c)]", "diaryl ether O (the phenyl ester O, as the GROMACS OPLS-AA topologies type diphenyl ether)")
by_h((181, 182, 183, 184), 5, f"$(*{ETHO[:-1]};$(O([CX4])[CX4,c])])", "ether alpha C")
add(185, 5, f"[H][CX4]{ETHO[:-1]};$(O([CX4])[CX4,c])]", "H on ether alpha C")
add(199, 5, f"c{ETHO[:-1]};$(O(c)[CX4])]", "aryl ether C(O)")
add(472, 5, f"c{ETHO[:-1]};$(O(c)c)]", "diaryl ether C(O) (the phenyl ester C(ipso))")
# vinyl ethers (ethyl vinyl ether): the table has the vinyl C 517/518 but no O; the GROMACS OPLS-AA topologies give the O
# the acetal O 186, the vinyl CH 518 and the alkoxy C the ether types
add(186, 7, "[OX2;$(O[CX3]=[CX3]);$(O[CX4]);!$(O[CX3]=[O,S,N])]", "vinyl ether O")
add(518, 7, "[CX3;$(*[OX2][CX4]);!$(*=[O,S,N])]=[CX3]", "vinyl ether C(O)=")
by_h((181, 182, 183, 184), 7, "$(*[OX2][CX3]=[CX3])", "alkoxy C of a vinyl ether")
add(185, 7, "[H][CX4][OX2][CX3]=[CX3]", "H on the alkoxy C of a vinyl ether")

# ---------------------------------------------------------------- acetals and hemiacetals (OPLS/2020): an sp3 C bearing two O
# (OCH2O in poly(methylene oxide), a 1,3-dioxolane, a hemiacetal end) takes the acetal types, its O the acetal O 186;
# the carbons on the other side of that O stay ether alpha C (181-185), as in dimethoxymethane
ACO = "[OX2H0;!$(O[#6]=[O,S,N]);!$(O[Si]);!$(Oc)]"
ACET = f"$(*({ACO}){ACO})"
HEMI = f"$(*({ACO})[OX2H1])"
by_h((None, 189, 193, 197), 7, ACET, "acetal C")
add(190, 7, f"[H][CX4H2;{ACET}]", "H on acetal C (OCH2O)")
add(194, 7, f"[H][CX4H1;{ACET}]", "H on acetal C (OCHRO)")
by_h((None, 191, 195, 198), 7, HEMI, "hemiacetal C")
add(192, 7, f"[H][CX4H2;{HEMI}]", "H on hemiacetal C (OCH2OH)")
add(196, 7, f"[H][CX4H1;{HEMI}]", "H on hemiacetal C (OCHROH)")
add(186, 8, f"{ACO[:-1]};$(O([CX4][OX2])[CX4][OX2]);!$(O[CX4][OX2H1])]", "O between two acetal C (paraldehyde, 1,3,5-trioxane), as those topologies")
add(180, 7, f"{ACO[:-1]};$(O[CX4][OX2]);!$(O[CX4][OX2H1])]", "acetal ether O (the dialkyl ether O, as the GROMACS OPLS-AA topologies type dimethoxymethane and 1,3-dioxolane)")
add(186, 7, f"{ACO[:-1]};$(O[CX4][OX2H1])]", "hemiacetal ether O")
add(187, 7, f"[OX2H1;$(O[CX4]{ACO})]", "hemiacetal O(H)")
add(188, 7, f"[H][OX2H1;$(O[CX4]{ACO})]", "hemiacetal H(O)")

# ---------------------------------------------------------------- carbonyls: aldehydes, ketones
add(277, 5, "[CX3H1](=[OX1])[#6]", "aldehyde C")
add(278, 5, "[OX1]=[CX3H1][#6]", "aldehyde O")
add(279, 5, "[H][CX3]=[OX1]", "aldehyde H")
add(277, 5, "[CX3H2]=[OX1]", "formaldehyde C")
add(278, 5, "[OX1]=[CX3H2]", "formaldehyde O")
add(280, 5, "[CX3H0](=[OX1])([#6])[#6]", "ketone C")
add(281, 5, "[OX1]=[CX3H0]([#6])[#6]", "ketone O")
add(282, 4, "[H][CX4][CX3](=[OX1])[#6,#1]", "H on C alpha to a ketone / aldehyde")
# esters take the ketone's alpha parameters on the acyl side (the table: "=O: esters ketone params (see 280-282)"); an
# acid's alpha H stays 140
# as the GROMACS OPLS-AA topologies: alkyl ketones, aldehydes and anhydrides 282; the H alpha to an ester's acyl C, or to
# a ketone C that also carries an aryl or cyclopropyl group (acetophenone, cyclopropyl methyl ketone), stays 140
add(140, 5, "[H][CX4][CX3](=[OX1])[OX2][#6;!$([CX3]=[OX1])]", "H alpha to an ester carbonyl (acyl side)")
add(140, 5, "[H][CX4][CX3](=[OX1])[c,$([CX4;r3])]", "H alpha to an aryl or cyclopropyl ketone")
add(282, 5, "[H][CX4][CX3](=[OX1])[OX2][CX3]=[OX1]", "H alpha to an anhydride carbonyl")
add(140, 5, "[H][CX4;r3][CX3]=[OX1]", "H on a cyclopropyl C alpha to a ketone")
# diesters (glycol diacetates, diethyl malonate): the acyl-side H takes 282 there, as those topologies have it
for n in (2, 3, 4):
    add(282, 6, "[H][CX4][CX3](=[OX1])[OX2]" + "[CX4]" * n + "[OX2][CX3]=[OX1]", f"H alpha to an ester carbonyl of a diester (O-(C){n}-O)")
add(282, 6, "[H][CX4]([CX3](=[OX1])[OX2])[CX3](=[OX1])[OX2]", "H on the C between two ester carbonyls (malonate)")
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
# vinyl esters (vinyl acetate): the table has no vinyl ester O; the GROMACS OPLS-AA topologies take the phenyl ester O 473
# and the alkene types for the vinyl group
add(473, 7, "[OX2;$(O[CX3](=[OX1])[#6,#1]);$(O[CX3]=[CX3])]", "ester -O- (vinyl ester)")
add(142, 7, "[CX3H1;$(*[OX2][CX3]=[OX1])]=[CX3]", "vinyl ester C(O)= (RHC=)")
add(141, 7, "[CX3H0;$(*[OX2][CX3]=[OX1])]=[CX3]", "vinyl ester C(O)= (R2C=)")
add(467, 7, "[OX2;$(O([CX3]=[OX1])[CX3]=[OX1])]", "anhydride -O-")
add(472, 7, "c[OX2][CX3](=[OX1])[#6,#1]", "phenyl ester C(ipso)")
by_h((468, 490, 491, 492), 6, "$(*[OX2][CX3](=[OX1])[#6,#1])", "ester alkoxy C")
add(469, 6, "[H][CX4][OX2][CX3](=[OX1])[#6,#1]", "H on ester alkoxy C")

# ---------------------------------------------------------------- carboxylic acids
add(267, 7, "[CX3](=[OX1])[OX2H1]", "carboxylic acid C")
add(269, 7, "[OX1]=[CX3][OX2H1]", "carboxylic acid =O")
add(268, 7, "[OX2H1][CX3]=[OX1]", "carboxylic acid -OH")
add(270, 7, "[H][OX2H1][CX3]=[OX1]", "carboxylic acid H")
CXO = "$([CX3](~[OX1])~[OX1-])"
add(271, 7, f"[CX3;{CXO}]", "carboxylate C")
add(272, 7, f"[OX1;$(O~[CX3;{CXO}])]", "carboxylate O")
by_h((135, 136, 137, 139), 7, f"$(*[CX3;{CXO}])", "C alpha to a carboxylate (the alkane types, as the GROMACS OPLS-AA topologies type acetate)")

# ---------------------------------------------------------------- carbonates (dialkyl; aryl carbonates have no OPLS-AA type)
# cyclic carbonates (ethylene and propylene carbonate) 771-779; an open-chain dialkyl carbonate takes the ester types
# (465-469, 490), as the GROMACS OPLS-AA topologies type dimethyl and diethyl carbonate
add(772, 7, "[CX3;R](=[OX1])([OX2][CX4])[OX2][CX4]", "cyclic carbonate C=O")
add(771, 7, "[OX1]=[CX3;R]([OX2][CX4])[OX2][CX4]", "cyclic carbonate =O")
add(773, 7, "[OX2;R;$(O[CX3](=[OX1])[OX2][CX4]);$(O[CX4])]", "cyclic carbonate -O-")
by_h((789, 774, 775, None), 7, "$(*[OX2;R][CX3](=[OX1])[OX2][CX4])", "cyclic carbonate alkoxy C")
add(777, 7, "[H][CX4][OX2;R][CX3](=[OX1])[OX2][CX4]", "H on cyclic carbonate alkoxy C")
add(465, 7, "[CX3;!R](=[OX1])([OX2][CX4])[OX2][CX4]", "open-chain carbonate C=O (ester C)")
add(466, 7, "[OX1]=[CX3;!R]([OX2][CX4])[OX2][CX4]", "open-chain carbonate =O (ester =O)")
add(467, 7, "[OX2;!R;$(O[CX3](=[OX1])[OX2][CX4]);$(O[CX4])]", "open-chain carbonate -O- (ester -O-)")
by_h((468, 490, 491, 492), 7, "$(*[OX2;!R][CX3](=[OX1])[OX2][CX4])", "open-chain carbonate alkoxy C (ester alkoxy C)")
add(469, 7, "[H][CX4][OX2;!R][CX3](=[OX1])[OX2][CX4]", "H on open-chain carbonate alkoxy C")
add(778, 8, "[H][CX4H1][OX2;R][CX3](=[OX1])[OX2][CX4]", "H on carbonate alkoxy CH")
add(776, 7, "[CX4H3;$(*[CX4;R][OX2][CX3](=[OX1])[OX2][CX4])]", "CH3 on a cyclic carbonate ring C (propylene carbonate)")
add(779, 7, "[H][CX4H3;$(*[CX4;R][OX2][CX3](=[OX1])[OX2][CX4])]", "H of that CH3")

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
add(901, 5, f"[NX3H1;{AMN}]([CX4,c])[CX4,c]", "secondary amine N (alkyl or aryl: N-methylaniline, diarylamine)")
add(902, 5, f"[NX3H0;{AMN}]([CX4,c])([CX4,c])[CX4,c]", "tertiary amine N (alkyl or aryl: N,N-dimethylaniline)")
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
# a nitrile on an alkene C (acrylonitrile): the table has no numbers; the conjugated benzonitrile C and N, as the GROMACS
# OPLS-AA topologies (Caleman et al. 2012) take them
add(262, 6, "[NX1]#[CX2][CX3]=[CX3]", "conjugated nitrile N (acrylonitrile)")
add(261, 6, "[CX2](#[NX1])[CX3]=[CX3]", "conjugated nitrile C (acrylonitrile)")
add(260, 6, "[CX3;$(*[CX2]#[NX1])]=[CX3]", "alkene C carrying a nitrile (the benzonitrile ipso C, as the GROMACS OPLS-AA topologies type acrylonitrile)")
add(140, 7, "[H][CX3;$(*[CX2]#[NX1])]=[CX3]", "H on that C (140, as those topologies)")

# ---------------------------------------------------------------- halogens (one halogen per carbon; CX2 / CX3 groups only when perfluorinated)
ONEX = lambda X: f"$(*[{X}]);!$(*([{X}])[{X}]);!$(*[F,Cl,Br,I;!{X}])"
add(151, 5, "[Cl][CX4]", "alkyl chloride Cl")
by_h((None, 152, 973, 974), 5, ONEX("Cl"), "alkyl chloride C")
add(153, 5, "[H][CX4H2][Cl]", "H on RCH2Cl")
add(972, 5, "[H][CX4H1][Cl]", "H on R2CHCl")
# polychlorinated carbons (CH2Cl2, CHCl3, CHCl2-, CCl3-): the table has no numbers of its own; the GROMACS OPLS-AA
# topologies (Caleman et al. 2012) give every such carbon RCH2Cl's 152, its H 153 and the Cl 151
by_h((152, 152, 152, 152), 6, "$(*(Cl)Cl);!$(*[F,Br,I])", "polychlorinated C")
add(153, 6, "[H][CX4;$(*(Cl)Cl)]", "H on a polychlorinated C")
# a carbon with both F and Cl (CHFCl2): the fluoroalkyl C by H count, F 965, H 153, as those topologies have it
by_h((None, 957, 959, 960), 6, "$(*F);$(*Cl);!$(*(F)F)", "fluorochloro C")
add(965, 7, "[F][CX4;$(*Cl)]", "F on a C that also carries Cl")
add(153, 7, "[H][CX4;$(*F);$(*Cl)]", "H on a fluorochloro C")
# alkyl bromides and iodides as the GROMACS OPLS-AA topologies type them: Br 722, I 732, the carbon RCH2Cl's 152 and its
# H 153 whatever the number of H
add(722, 5, "[Br][CX4]", "alkyl bromide Br")
by_h((152, 152, 152, 152), 5, "$(*[Br,I]);!$(*[F,Cl])", "C of an alkyl bromide or iodide")
add(153, 5, "[H][CX4][Br,I]", "H on an alkyl bromide or iodide C")
add(732, 5, "[I][CX4]", "alkyl iodide I")
add(956, 5, "[F][CX4;$(*[F]);!$(*(F)F)]", "monofluoroalkyl F")
by_h((None, 957, 959, 960), 5, ONEX("F"), "monofluoroalkyl C")
add(958, 5, "[H][CX4;$(*[F]);!$(*(F)F)]", "H on monofluoroalkyl C")
add(965, 6, "[F][CX4;$(*(F)F)]", "F of a CF2 / CF3 group")
add(962, 6, "[CX4;$(*(F)F);!$(*(F)(F)F);!$(*[H])]", "CF2")
# a CF2H group (a fluoropolymer's H-capped end, HCF2-): the table's CF2H of difluoromethylbenzene, 966 / 967
add(966, 6, "[CX4H1;$(*(F)F);!$(*(F)(F)F)]", "CF2H carbon")
add(967, 6, "[H][CX4H1;$(*(F)F);!$(*(F)(F)F)]", "H on a CF2H carbon")
add(961, 6, "[CX4;$(*(F)(F)F);!$(*(F)(F)(F)F);$(*[#6])]", "CF3 on carbon (a perfluoroalkyl end)")
add(264, 6, "[Cl]c", "aryl chloride Cl")
# fluorobenzenes: a ring C-F with one C-F neighbour takes the difluorobenzene types, with two the hexafluorobenzene
# ones (1,2-difluoro-, 1,2,3,4- and 1,2,3,5-tetrafluorobenzene in the GROMACS OPLS-AA topologies)
add(727, 7, "[c;$(c[F]);$(c:c[F]);!$(c(:c[F]):c[F])]", "C(F) next to one C(F) (difluorobenzenes)")
add(728, 7, "[F][c;$(c:c[F]);!$(c(:c[F]):c[F])]", "F on a C(F) next to one C(F)")
add(720, 8, "[c;$(c[F]);$(c(:c[F]):c[F])]", "C(F) between two C(F) (hexafluorobenzene)")
add(721, 8, "[F][c;$(c(:c[F]):c[F])]", "F on a C(F) between two C(F)")
add(724, 8, "c[CX4](F)(F)F", "ring C of a trifluoromethylbenzene")
add(725, 8, "[CX4;$(*(F)(F)(F)c)]", "CF3 on a benzene ring")
add(726, 8, "[F][CX4;$(*(F)(F)(F)c)]", "F of a CF3 on a benzene ring")
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
add(226, 7, "[Cl][CX3;$(*(Cl)Cl)]=[CX3]", "chloroalkene Cl (Cl2C=: 226, as the GROMACS OPLS-AA topologies type 1,1-dichloroethene)")
add(227, 7, "[CX3;$(*(Cl)Cl)]=[CX3]", "chloroalkene C (Cl2C=: 227)")

# ---------------------------------------------------------------- sulfur
add(493, 6, "[SX4](=[OX1])(=[OX1])([#6])[#6]", "sulfone S")
add(494, 6, "[OX1]=[SX4](=[OX1])([#6])[#6]", "sulfone O")
add(488, 7, "c[SX4](=[OX1])(=[OX1])", "ring C on a sulfonyl S (charge 0: the SO2 group is neutral)")
add(202, 5, "[SX2]([CX4])[CX4]", "sulfide S")
add(203, 6, "[SX2]([#6])[SX2]", "disulfide S")
add(200, 5, "[SX2H1][CX4]", "thiol S")
add(204, 5, "[H][SX2H1]", "thiol H")
by_h((209, 210, 211, 212), 5, "$(*[SX2][CX4])", "C on a sulfide S")
add(496, 6, "[SX3](=[OX1])([CX4])[CX4]", "sulfoxide S")
add(497, 6, "[OX1]=[SX3]([CX4])[CX4]", "sulfoxide O")
add(495, 7, "[SX3](=[OX1])(c)[CX4]", "alkyl aryl sulfoxide S")
add(497, 7, "[OX1]=[SX3](c)[CX4]", "alkyl aryl sulfoxide O")
add(489, 7, "c[SX3](=[OX1])[CX4]", "ring C on a sulfoxide S")
by_h((498, 499, None, None), 6, "$(*[SX3]=[OX1])", "C on a sulfoxide S")
# the table has no number for the alkyl C of a sulfone; the GROMACS OPLS-AA topologies (Caleman et al. 2012, sulfolane)
# give it the sulfonamide alpha-C types, CH2 484 with H 485
add(484, 6, "[CX4H2;$(*[SX4](=[OX1])(=[OX1])[#6])]", "C alpha to a sulfone S (CH2)")
add(485, 6, "[H][CX4H2][SX4](=[OX1])(=[OX1])[#6]", "H on C alpha to a sulfone S")
by_h((213, 214, 215, 216), 6, "$(*[SX2][SX2])", "C on a disulfide S")
by_h((None, 206, 207, 208), 5, "$(*[SX2H1])", "C on a thiol S")

# ---------------------------------------------------------------- epoxides
# the table's epoxide types, with three-ring angle terms (CY-O^-CY 90 deg). The GROMACS OPLS-AA topologies type
# methyloxirane with the ether types instead; their 109.5 deg angles put about 119 kcal/mol of strain in each ring
# (44 with these), so CAPS keeps the epoxide types: the one place it differs from those topologies
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

# trialkyl phosphates (triethyl phosphate): the table's dimethyl phosphate numbers, as the GROMACS OPLS-AA topologies
PO4 = "[PX4;$(P(=[OX1])([OX2][#6])([OX2][#6])[OX2][#6])]"
add(440, 6, PO4, "phosphate triester P")
add(441, 6, f"[OX1]={PO4}", "phosphate =O")
add(442, 6, f"[OX2;$(O([#6]){PO4})]", "phosphate ester -O-")
by_h((443, 443, 443, None), 6, f"$(*[OX2]{PO4})", "C on a phosphate ester O")
add(444, 6, f"[H][CX4][OX2]{PO4}", "H on C of a phosphate ester O")

# alkynes: sp C 925 (939 when both ends carry carbon), terminal H 926, the sp3 C on an alkyne C 927-930 (their H 140)
add(925, 6, "[CX2]#[CX2]", "alkyne C")
add(939, 7, "[CX2;$(*([#6])#[CX2][#6])]", "alkyne C of RC#CR")
add(926, 6, "[H][CX2]#[CX2]", "alkyne terminal H")
# the sp3 C on an alkyne C keeps the alkane types (135-139), as the GROMACS OPLS-AA topologies type propyne

# monatomic ions
for z, n, d in (("[Cl-,Cl;X0]", 401, "Cl-"), ("[Br-,Br;X0]", 402, "Br-"), ("[Li+,Li;X0]", 406, "Li+"), ("[Na+,Na;X0]", 407, "Na+"), ("[K+,K;X0]", 408, "K+")):
    add(n, 8, z, d)

# pyridine (520-526), its methyl derivatives (670 2-methyl; 3- and 4-methyl as in the GROMACS OPLS-AA topologies: 672, 674),
# pyrimidine (530-536) and quinoline (603-619) by ring position
add(520, 4, "[nX2;$(n1ccccc1)]", "pyridine N")
add(521, 4, "[c;$(c1ncccc1)]", "pyridine C2")
add(522, 4, "[c;$(c1cnccc1)]", "pyridine C3")
add(523, 4, "[c;$(c1ccncc1)]", "pyridine C4")
add(524, 4, "[H][c;$(c1ncccc1)]", "pyridine H2")
add(525, 4, "[H][c;$(c1cnccc1)]", "pyridine H3")
add(526, 4, "[H][c;$(c1ccncc1)]", "pyridine H4")
# pyrazine (527-529) and neutral imidazole (557-565): the parent rings, every C bearing its H (a substituted ring has other
# numbers or none, and stays untyped); the charges of each ring sum to zero
add(527, 4, "[nX2;$(n1ccncc1)]", "pyrazine N")
add(528, 4, "[cH1;$(c1nccnc1)]", "pyrazine C")
add(529, 4, "[H][cH1;$(c1nccnc1)]", "pyrazine H")
add(557, 4, "[nX3H1;$([nH]1c[nX2]cc1)]", "imidazole N1 (NH)")
add(558, 4, "[cH1;$(c1[nX2]cc[nH]1)]", "imidazole C2")
add(559, 4, "[nX2;$(n1c[nH]cc1)]", "imidazole N3")
add(560, 4, "[cH1;$(c1[nX2]c[nH]c1)]", "imidazole C4")
add(561, 4, "[cH1;$(c1[nH]c[nX2]c1)]", "imidazole C5")
add(562, 4, "[H][nX3H1;$([nH]1c[nX2]cc1)]", "imidazole H1 (on N1)")
add(563, 4, "[H][cH1;$(c1[nX2]cc[nH]1)]", "imidazole H2")
add(564, 4, "[H][cH1;$(c1[nX2]c[nH]c1)]", "imidazole H4")
add(565, 4, "[H][cH1;$(c1[nH]c[nX2]c1)]", "imidazole H5")
# urea O=C(NH2)2 (247-250; Duffy et al., Isr. J. Chem. 33, 323 (1993)): above the amide rules
UREA = "[CX3;$(C(=[OX1])([NX3H2])[NX3H2])]"
add(247, 9, UREA, "urea C")
add(248, 9, f"[OX1]={UREA}", "urea O")
add(249, 9, f"[NX3H2]{UREA}", "urea N")
add(250, 9, f"[H][NX3H2]{UREA}", "urea H")
add(670, 7, "[CX4H3;$(*c1ncccc1)]", "CH3 on pyridine C2")
add(672, 7, "[CX4H3;$(*c1cnccc1)]", "CH3 on pyridine C3")
add(674, 7, "[CX4H3;$(*c1ccncc1),$(*c1ncncc1)]", "CH3 on pyridine C4 / pyrimidine C4")
add(530, 8, "[n;$(n1cnccc1)]", "pyrimidine N")
add(531, 8, "[c;$(c1ncccn1)]", "pyrimidine C2")
add(532, 8, "[c;$(c1ncncc1)]", "pyrimidine C4")
add(533, 8, "[c;$(c1cncnc1)]", "pyrimidine C5")
add(534, 8, "[H][c;$(c1ncccn1)]", "pyrimidine H2")
add(535, 8, "[H][c;$(c1ncncc1)]", "pyrimidine H4")
add(536, 8, "[H][c;$(c1cncnc1)]", "pyrimidine H5")
QUIN = {603: "n1cccc2ccccc12", 604: "c1nc2ccccc2cc1", 605: "c1cnc2ccccc2c1", 606: "c1ccnc2ccccc12", 607: "c1cccc2ncccc12",
        608: "c1ccc2ncccc2c1", 609: "c1ccc2cccnc2c1", 610: "c1cccc2cccnc12", 611: "c12ncccc1cccc2", 612: "c12ccccc1nccc2"}
for n, sm in QUIN.items():
    add(n, 8, f"[{sm[0]};$({sm})]", f"quinoline {'N1' if n == 603 else 'C' + str(n - 602)}")
for c, h in ((604, 613), (605, 614), (606, 615), (607, 616), (608, 617), (609, 618), (610, 619)):
    add(h, 8, f"[H][c;$({QUIN[c]})]", f"quinoline H{h - 611}")

# isoquinoline: quinoline's numbers position by position, as the GROMACS OPLS-AA topologies give them
ISOQ = {603: "n1cc2ccccc2cc1", 611: "c1nccc2ccccc12", 604: "c1ncc2ccccc2c1", 605: "c1cncc2ccccc12", 612: "c12ccccc1cncc2",
        606: "c1cccc2cnccc12", 607: "c1ccc2cnccc2c1", 608: "c1ccc2ccncc2c1", 609: "c1cccc2ccncc12", 610: "c12cnccc1cccc2"}
for n, sm in ISOQ.items():
    add(n, 8, f"[{sm[0]};$({sm})]", f"isoquinoline atom typed as quinoline {n}")
for c, h in ((611, 619), (604, 613), (605, 614), (606, 615), (607, 616), (608, 617), (609, 618)):
    add(h, 8, f"[H][c;$({ISOQ[c]})]", f"isoquinoline H typed as quinoline {h}")

# thiophene: no numbers of its own in the table; the GROMACS OPLS-AA topologies take furan's C and H and thiazole's S
add(633, 7, "[s;r5;$(s1cccc1)]", "thiophene S")
add(567, 7, "[c;r5;$(c1sccc1)]", "thiophene C2")
add(568, 7, "[c;r5;$(c1cscc1)]", "thiophene C3")
add(569, 7, "[H][c;r5;$(c1sccc1)]", "thiophene H2")
add(570, 7, "[H][c;r5;$(c1cscc1)]", "thiophene H3")

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
                              "number for (siloxanes, aryl carbonates, CCl2 carbons, substituted ureas, imides, carbamates) are left untyped and reported. "
                              "Among matching rules the highest priority wins.",
               "rules": rules}
        p = os.path.join(out, ff + ".typing.json")
        json.dump(doc, open(p, "w"), ensure_ascii=False, indent=1)
        print(f"{p}: {len(rules)} rules" + (f" (numbers not in this table: {', '.join(missing)})" if missing else ""))
