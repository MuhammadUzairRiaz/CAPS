#!/usr/bin/env python3
"""CAPS typing rules for the AMBER protein / nucleic-acid force fields built from AmberTools' own files
(bench/ff/open/build_amber.sh: amber-ff03-ambertools), in AMBER's own type names, checked against AmberTools' residue
libraries (bench/ff/validate_off_types.py: all_amino03.lib, all_aminont94.lib, all_aminoct94.lib, all_nucleic94.lib).

The types follow Cornell et al., J. Am. Chem. Soc. 117, 5179 (1995): CT for every sp3 carbon; its hydrogens HC, H1, H2,
H3 by the number of electron-withdrawing neighbours (N, O, S, halogen) of the carbon, HP next to an ammonium, H0 for
ff03's glycine; the aromatic and conjugated carbons by ring and neighbours (CA, CB, CN, C*, CW, CV, CC, CR, CK, CQ, CM),
their hydrogens HA / H4 / H5 by electron-withdrawing neighbours; C / O carbonyls, O2 carboxylate and phosphate oxygens,
OS ethers and esters, OH / HO hydroxyls; N amide, N2 conjugated amino, N3 sp3 amine, NA ring N-H, NB / NC ring N of five
and six rings, N* glycosidic N; S / SH, P; TIP3P water OW / HW and the ions of parm99.

usage: make_amber_rules.py [OUT_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_biomolecular_rules import Rules, glycam as glycam_rules

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data", "typing")
EW = "#7,#8,#16,#9,#17,#35,#53"   # AMBER's electron-withdrawing neighbours


def ew(n, exact=True):
    """a carbon with n electron-withdrawing neighbours (at least n when not exact)"""
    have = "$(*" + "".join(f"(~[{EW}])" for _ in range(n - 1)) + f"~[{EW}])"
    more = "$(*" + "".join(f"(~[{EW}])" for _ in range(n)) + f"~[{EW}])"
    return have + (f";!{more}" if exact else "")


N_IN_RING6 = "$(*1~[#7]~*~*~*~*~1),$(*1~*~[#7]~*~*~*~1),$(*1~*~*~[#7]~*~*~1)"   # a six-ring holding an N
BENZENE_FUSED = "$(*1:c:c:c:c:c:1)"


def amber():
    R = Rules()
    # carbon
    R.add("CT", 0, "[#6X4]", "sp3 C")
    R.add("CA", 1, "[#6X3]", "sp2 C (aromatic, guanidinium)")
    R.add("C", 3, "[#6X3]=[OX1]", "carbonyl C (amide, ester, ketone, ring C=O)")
    R.add("C", 3, "[#6X3](~[OX1])~[OX1]", "carboxylate C")
    R.add("C", 3, "[#6X3;R](-[OX2H1,OX1-])", "aromatic C bearing O (tyrosine CZ)")
    R.add("CB", 4, "[#6X3;r5;r6;!$(*~[#7;r6])]", "ring fusion C of indole away from N (Trp CD2), purine C5")
    R.add("CB", 4, "[#6X3;r5;r6;$(*~[#7;r6])]", "purine C4 / C5")
    R.add("CN", 5, f"[#6X3;r5;r6;{BENZENE_FUSED};$(*~[#7;r5])]", "indole ring fusion C bonded to N (Trp CE2)")
    R.add("C*", 4, "[#6X3;r5;!r6;$(*-[CX4]);$(*~[#6;r5;r6])]", "indole C3 (Trp CG)")
    R.add("CW", 4, "[#6X3;r5;!r6;$(*~[#7X3H1;r5]);!$(*(~[#7])~[#7])]", "five-ring C=C next to N-H (Trp CD1, His CD2 of HIE / HIP)")
    R.add("CV", 4, "[#6X3;r5;!r6;$(*~[#7X2;r5]);!$(*(~[#7])~[#7]);!$(*-[CX4])]", "five-ring C next to an N without H (His CD2 of HID)")
    R.add("CC", 5, "[#6X3;r5;!r6;$(*-[CX4]);$(*~[#7;r5])]", "imidazole C bearing the side chain (His CG)")
    R.add("CR", 5, "[#6X3;r5;!r6;$(*(~[#7;r5])~[#7;r5])]", "imidazole C between two N (His CE1)")
    R.add("CK", 6, "[#6X3;r5;!r6;$(*(~[#7;r5]~[#6;r5;r6])~[#7;r5]~[#6;r5;r6])]", "purine C8")
    R.add("CM", 2, f"[#6X3;r6;!r5;{N_IN_RING6}]", "pyrimidine C5 / C6")
    R.add("CQ", 6, "[#6X3H1;r6;$(*(~[#7;r6])~[#7;r6])]", "purine C2 between two ring N")
    R.add("CA", 6, "[#6X3;R;$(*-[NX3H2,NX3H1;!R])]", "ring C bearing an amino group (A C6, C C4, G C2)")
    R.add("C", 7, "[#6X3;R]=[OX1]", "ring carbonyl C")
    # hydrogen
    R.add("HC", 0, "[H][#6X4]", "H on sp3 C")
    R.add("H1", 1, f"[H][#6X4;{ew(1)}]", "H on sp3 C with one electron-withdrawing neighbour")
    R.add("H2", 1, f"[H][#6X4;{ew(2)}]", "H on sp3 C with two")
    R.add("H3", 1, f"[H][#6X4;{ew(3, False)}]", "H on sp3 C with three")
    R.add("HP", 2, "[H][#6X4][#7X4+]", "H on C next to an ammonium N")
    R.add("HP", 2, "[H][#6X4][#7X3;!$(*~[#6X3]);!$(*~[#6X2])]", "H on C next to an sp3 amine N (neutral Lys CE)")
    R.add("H0", 3, "[H][#6X4H2]([#7X3][#6X3]=[OX1])[#6X3](=[OX1])[#7X3]", "glycine H alpha inside a chain (ff03)")
    R.add("HA", 0, "[H][#6X3]", "H on sp2 C")
    R.add("H4", 1, f"[H][#6X3;{ew(1)}]", "H on sp2 C with one electron-withdrawing neighbour")
    R.add("H5", 1, f"[H][#6X3;{ew(2, False)}]", "H on sp2 C with two")
    R.add("H", 0, "[H][#7]", "H on N")
    R.add("HO", 0, "[H][#8]", "hydroxyl H")
    R.add("HS", 0, "[H][#16]", "thiol H")
    # nitrogen
    R.add("N3", 0, "[#7X4]", "sp3 N (ammonium)")
    R.add("N3", 0, "[#7X3;!$(*~[#6X3]);!$(*~[#6X2])]", "sp3 amine N")
    R.add("N2", 1, "[#7X3;!R;$(*~[#6X3])]", "sp2 amino N (guanidinium, base amines)")
    R.add("N", 2, "[#7X3;$(*-[#6X3]=[OX1])]", "amide N")
    R.add("NA", 4, "[#7X3H1;R]", "ring N with H")
    R.add("NB", 4, "[#7X2;r5]", "five-ring N without H")
    R.add("NC", 4, "[#7X2;r6]", "six-ring N without H")
    R.add("N*", 5, "[#7X3;R;$(*-[#6X4]);$(*@[#6X3])]", "glycosidic ring N")
    # oxygen
    R.add("O", 0, "[OX1]", "carbonyl O")
    R.add("O2", 1, "[OX1]~[#6X3]~[OX1]", "carboxylate O")
    R.add("O2", 1, "[OX1]~[#15]", "phosphate O"); R.add("O2", 1, "[OX1]~[#16X4]", "sulfate O")
    R.add("OS", 0, "[OX2]", "ether / ester O")
    R.add("OH", 1, "[OX2H1]", "hydroxyl O")
    # sulfur, phosphorus
    R.add("S", 0, "[#16]", "sulfide / disulfide S")
    R.add("SH", 1, "[#16X2H1]", "thiol S"); R.add("SH", 1, "[#16X1]", "thiolate S (CYM)")
    R.add("P", 0, "[#15]", "phosphate P")
    # water (TIP3P) and parm99's ions
    R.add("OW", 8, "[OX2H2]", "water O"); R.add("HW", 8, "[H][OX2H2]", "water H")
    for t, s in [("IP", "[Na+]"), ("IM", "[Cl-]"), ("Li", "[Li+]"), ("K", "[K+]"), ("Rb", "[Rb+]"), ("Cs", "[Cs+]"),
                 ("MG", "[Mg+2]"), ("C0", "[Ca+2]"), ("Zn", "[Zn+2]"), ("F", "[F-]"), ("Br", "[Br-]"), ("I", "[I-]")]:
        R.add(t, 9, s, "ion")
    return R


def glycam():
    """GLYCAM06 (Kirschner et al., J. Comput. Chem. 29, 622 (2008)) on AMBER's rules: a carbohydrate's sp3 C Cg (Cy at a
    sialic acid's C2), their H Hc / H1 / H2, hydroxyl Oh / Ho, ether and glycosidic O Os (Oy in a sialic acid ring),
    amide N Ng (bench/typing/make_biomolecular_rules.py glycam(), in GLYCAM's own names)"""
    R = amber()
    G = Rules()
    glycam_rules(G)
    for r in G:
        r["type"] = r["type"].removeprefix("G_")
    # the anomeric H of a 3-deoxy pyranose is not H1 as a class (GLYCAM types abequose's H2, tyvelose's H1: below)
    G = Rules(r for r in G if not r["description"].startswith("anomeric H of a 3-deoxy pyranose"))
    # GLYCAM types the anomeric H H1, not H2, in these sugars (and their enantiomers), told apart from their epimers by
    # the configuration of every ring carbon (@ / @@ from the structure's positions); written from GLYCAM_06j-1.prep's
    # residues: the ring from the anomeric C round to the ring O, each carbon's H and O (or CH3 / CH2) in that order
    for code, sugar, smarts in [
        ("?DD", "alpha-D-lyxofuranose", "[H][C@@]1([OX2])[C@@]([H])([OX2])[C@@]([H])([OX2])[C@@]([H])([CX4H2])[OX2]1"),
        ("?dD", "alpha-L-lyxofuranose", "[H][C@]1([OX2])[C@]([H])([OX2])[C@]([H])([OX2])[C@]([H])([CX4H2])[OX2]1"),
        ("?XD", "alpha-D-xylofuranose", "[H][C@@]1([OX2])[C@]([H])([OX2])[C@@]([H])([OX2])[C@@]([H])([CX4H2])[OX2]1"),
        ("?xD", "alpha-L-xylofuranose", "[H][C@]1([OX2])[C@@]([H])([OX2])[C@]([H])([OX2])[C@]([H])([CX4H2])[OX2]1"),
        ("?RU", "beta-D-ribofuranose", "[H][C@]1([OX2])[C@]([H])([OX2])[C@]([H])([OX2])[C@@]([H])([CX4H2])[OX2]1"),
        ("?rU", "beta-L-ribofuranose", "[H][C@@]1([OX2])[C@@]([H])([OX2])[C@@]([H])([OX2])[C@]([H])([CX4H2])[OX2]1"),
        ("?QA / ?QB", "D-quinovopyranose (6-deoxy-glucose), either anomer",
         "[H][CX4]1([OX2])[C@]([H])([OX2])[C@@]([H])([OX2])[C@]([H])([OX2])[C@@]([H])([CX4H3])[OX2]1"),
        ("?qA / ?qB", "L-quinovopyranose, either anomer",
         "[H][CX4]1([OX2])[C@@]([H])([OX2])[C@]([H])([OX2])[C@@]([H])([OX2])[C@]([H])([CX4H3])[OX2]1"),
        ("?TV / ?Tv", "D-tyvelose (3,6-dideoxy-arabino-hexopyranose), either anomer",
         "[H][CX4]1([OX2])[C@@]([H])([OX2])[CX4]([H])([H])[C@]([H])([OX2])[C@@]([H])([CX4H3])[OX2]1"),
        ("?tV / ?tv", "L-tyvelose, either anomer",
         "[H][CX4]1([OX2])[C@]([H])([OX2])[CX4]([H])([H])[C@@]([H])([OX2])[C@]([H])([CX4H3])[OX2]1")]:
        G.add("H1", 7, smarts, f"anomeric H of {sugar}: H1 as GLYCAM types it ({code})")
    return Rules(list(R) + list(G))


if __name__ == "__main__":
    for fid, label, make in [("amber-ff03-ambertools", "AMBER ff03", amber), ("glycam06j-ambertools", "GLYCAM06j", glycam)]:
        rules = make()
        doc = {"format": "caps-typing", "version": 1, "forcefield": label, "unknown_types": "untyped",
               "description": f"CAPS rules for {label} in AMBER's own type names, checked against AmberTools' residue "
                              "libraries (bench/typing/make_amber_rules.py)", "rules": list(rules)}
        json.dump(doc, open(os.path.join(OUT, fid + ".typing.json"), "w"), ensure_ascii=False, indent=1)
        print(f"{fid}: {len(rules)} rules")
