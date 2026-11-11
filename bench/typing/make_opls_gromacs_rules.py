#!/usr/bin/env python3
"""CAPS typing rules for OPLS-AA built from GROMACS's oplsaa.ff (bench/ff/convert_gromacs_opls.py: oplsaa-gromacs):
CAPS's OPLS-AA chemistry rules (the rules R of make_oplsaa_rules.py, which name Jorgensen's type numbers) written as
GROMACS's opls_NNN, keeping the rules whose type the GROMACS file has, plus the residue types of GROMACS's
own amino-acid topologies (checked with bench/ff/validate_rtp_types.py against aminoacids.rtp). Run it after
make_oplsaa_rules.py (which writes data/typing/oplsaa2024.typing.json from the same rules).

Groups GROMACS's numbers cannot type at an integer charge (a diaryl ether, an anhydride, chloroform, an epoxide, a
siloxane: the rules make_oplsaa_rules.py flags only2024, and the rules whose number the GROMACS file lacks) take the
OPLS-AA 2024 library's assignment, whose charges balance the group by construction: its types are added to
oplsaa-gromacs.json, a number the GROMACS file gives another atom or another charge renamed opls2024_NNN (GROMACS's
opls_165 is an H, the 2024 table's 165 styrene's ipso C), a CAPS supplement type caps_* under its own name, each with
"source" naming the fallback; a 2024 number GROMACS has with the same element, charge and bond class keeps GROMACS's
type. Bonded terms of classes only those types carry (the epoxide O$, the siloxane Si) come from the 2024 file with
"source" on each rule, which CAPS writes beside every term it gives (provenance). The GROMACS topology of such a
molecule is then GROMACS's numbers wherever they balance and the 2024 numbers in the fallback groups; never a residual
charge spread over atoms. A group neither table balances stays untyped and CAPS refuses the molecule.

usage: make_opls_gromacs_rules.py [OUT_DIR]
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, "typing")
FFPATH = os.path.join(DATA, "forcefields", "oplsaa-gromacs.json")
FALLBACK = "OPLS-AA 2024 fallback"   # the "source" of what this script adds to the GROMACS file (removed and rebuilt each run)
TAG = "@24"                          # a fallback type's first-choice classes: CA@24 (second choice CA)
GFF = json.load(open(FFPATH))
O24 = json.load(open(os.path.join(DATA, "forcefields", "oplsaa2024.json")))
NATIVE = {t["name"]: t for t in GFF["atom_types"] if not t.get("source", "").startswith(FALLBACK)}
T24 = {t["name"]: t for t in O24["atom_types"]}
have = set(NATIVE)


def opls(t):
    return "opls_%03d" % int(t) if re.fullmatch(r"\d+", t) else t


# every rule with its flags (make_oplsaa_rules.py's R)
_m = {}
exec(open(os.path.join(HERE, "make_oplsaa_rules.py")).read().split('\nif __name__ == "__main__":')[0], _m)
ALL = _m["R"]
WATER = {"9999", "9998"}   # TIP3P as the .par numbers it: GROMACS's own 111 / 112 (ADD below)


def same_atom(g, o):
    return g["element"] == o["element"] and abs(g["charge"] - o["charge"]) < 1e-9 and g["equivalence"]["bond"] == o["equivalence"]["bond"]


def gmx_name(num, fallback):
    """the GROMACS file's name for rule number num; fallback: the 2024 library's type, renamed where GROMACS's number is
    another atom (None: neither file has it)"""
    t = opls(num)
    if num in WATER:
        return None
    if not fallback and t in have:
        return t
    if t not in T24:
        return None
    if t in have and same_atom(NATIVE[t], T24[t]):
        return t
    return t.replace("opls_", "opls2024_") if t.startswith("opls_") else t


# 2024-only rules for chemistry GROMACS's own numbers already type at an integer charge, so GROMACS keeps its typing
# there (each found where the fallback changed a molecule GROMACS typed whole): conjugated diene central C (its 142 /
# 141), RC#CR (its 931), the alkyl C of R-C#CH (its alkyne set 925-930 with alkane C, propyne), cyclobutane (its alkane
# numbers), the tertiary alcohol O (its 154, t-butanol), the vinyl ether O (its 186, ethyl vinyl ether), the acetal /
# hemiacetal O (its 186 at -0.40 with its 180: glucose, cellobiose, polyacetals)
GMX_OWN = {"150", "178", "160", "180", "939", "927", "928", "929", "930", "714", "715", "716", "186", "caps_os_acac"}
rules, dropped, used24 = [], set(), {}
for r in ALL:
    if r.get("only2024") and r["type"] in GMX_OWN:
        continue
    fb = bool(r.get("only2024")) or opls(r["type"]) not in have
    t = gmx_name(r["type"], fb)
    if t is None:
        dropped.add(opls(r["type"]))
        continue
    if t not in have:
        used24[t] = opls(r["type"])
    r = {k: v for k, v in r.items() if k not in ("gmx", "only2024")}
    if fb:
        r["description"] += " (OPLS-AA 2024 fallback: GROMACS's numbers do not balance this group)"
    ov = [gmx_name(x, fb) for x in r.get("overrides", [])]
    for x, y in zip(r.get("overrides", []), ov):
        if y is not None and y not in have:
            used24[y] = opls(x)
    rules.append({**r, "type": t, **({"overrides": [x for x in ov if x]} if "overrides" in r else {})})
# GROMACS's OPLS-AA/L types beyond the general rules (type, priority, smarts, description), as its amino-acid topologies
# (aminoacids.rtp) and ions.itp use them
AA = "[CX4;$(*-[#7X3,#7X4]);$(*-[#6X3](~[OX1])~[OX1,#7X3,OX2H1])]"   # a peptide C alpha (a protonated C-terminus too)
GUAN_C = "[CX3;$(*(~[NX3])(~[NX3])~[NX3])]"                           # guanidinium C (three sp2 N)
HIS_A = f"[#6X3;r5;$(*-[CX4H2]-{AA})]"                                 # the ring C carrying His / Trp's CB
IMID_RING = ",".join("$(*1~" + "~".join(f"[#{z}]" for z in seq) + "~1)" for seq in
                     ((6, 7, 6, 6), (7, 6, 6, 7), (6, 6, 7, 6), (6, 7, 6, 7), (7, 6, 7, 6)))   # any atom of an imidazole
IMID = f"{IMID_RING};$({HIS_A}),$(*@{HIS_A}),$(*@*@{HIS_A})"   # an atom of His's imidazole
NH_R5, N_R5 = "[#7X3H1;r5]", "[#7X2;r5]"
TRP_CG = f"[#6X3;r5;!r6;$(*-[CX4H2]-{AA});$(*~[#6;r5;r6])]"           # Trp's CG (on the indole's ring fusion)
ADD = [
    # C alpha
    ("opls_224B", 9, f"[CX4H1;!R;{AA[5:-1]}]", "peptide C alpha (all but Gly, Pro, Aib)"),
    ("opls_223B", 9, f"[CX4H2;!R;{AA[5:-1]}]", "glycine C alpha"),
    ("opls_225B", 9, f"[CX4H0;!R;{AA[5:-1]}]", "Aib C alpha"),
    # carboxylate ion alpha carbons (Asp CB, Glu CG)
    *[(f"opls_{t}", 6, f"[CX4H{h};$(*-[CX3](~[OX1])~[OX1]);!$(*-[CX3]~[OX2H1])]", f"C alpha to a carboxylate ion ({lbl})")
      for t, h, lbl in ((273, 3, "CH3"), (274, 2, "CH2"), (275, 1, "CH"), (276, 0, "C"))],
    # ammonium (Lys NZ+)
    ("opls_287", 6, "[NX4H3+,NX4H3]", "N of RNH3+"), ("opls_290", 6, "[H][NX4H3]", "H of RNH3+"),
    *[(f"opls_{t}", 10, f"[CX4H{h};$(*-[NX4H3])]", f"C of {lbl}NH3+") for t, h, lbl in ((291, 3, "CH3"), (292, 2, "RCH2"), (293, 1, "R2CH"), (294, 0, "R3C"))],
    # the termini of a peptide (OPLS-AA's terminal C alpha types, as GROMACS's and OPLS-AA/M's topologies use them):
    # the C alpha of a C-terminal carboxylate 283 (Gly 284, Pro 285), of an N-terminal NH3+ 293 (Gly 292, opls_291-294
    # by H count), an N-terminal proline's NH2+ 309 / 310, its C alpha 295 and C delta 296
    ("opls_283", 10, "[CX4H1;!R;$(*-[CX3](~[OX1])~[OX1]);!$(*-[CX3]~[OX2H1]);$(*-[#7X3,#7X4])]", "C-terminal C alpha (carboxylate)"),
    ("opls_284", 10, "[CX4H2;!R;$(*-[CX3](~[OX1])~[OX1]);!$(*-[CX3]~[OX2H1]);$(*-[#7X3,#7X4])]", "C-terminal glycine C alpha"),
    ("opls_285", 10, "[CX4H1;R;$(*-[CX3](~[OX1])~[OX1]);!$(*-[CX3]~[OX2H1]);$(*-[#7X3,#7X4;R])]", "C-terminal proline C alpha"),
    ("opls_309", 10, "[NX4H2+,NX4H2;R]", "N-terminal proline N (R2NH2+)"), ("opls_310", 7, "[H][NX4H2;R]", "H on R2NH2+"),
    ("opls_295", 10, "[CX4H1;R;$(*-[NX4H2]);$(*-[CX3]=[OX1])]", "N-terminal proline C alpha"),
    ("opls_296", 10, "[CX4H2;R;$(*-[NX4H2])]", "N-terminal proline C delta"),
    # arginine, charged: guanidinium
    ("opls_302", 6, GUAN_C, "guanidinium C+"),
    ("opls_300", 6, f"[NX3H2;$(*~{GUAN_C})]", "guanidinium NH2"), ("opls_301", 6, f"[H][NX3H2;$(*~{GUAN_C})]", "H on guanidinium NH2"),
    ("opls_303", 6, f"[NX3H1;$(*-[CX4]);$(*~{GUAN_C})]", "guanidinium NHR"), ("opls_304", 6, f"[H][NX3H1;$(*-[CX4]);$(*~{GUAN_C})]", "H on guanidinium NHR"),
    ("opls_307", 6, f"[CX4H2;$(*-[NX3H1]~{GUAN_C})]", "Arg CD (CH2 on the guanidinium)"),
    ("opls_140", 7, f"[H][CX4H2;$(*-[NX3H1]~{GUAN_C})]", "H on Arg CD"),
    ("opls_308", 6, f"[CX4H2;$(*-[CX4H2]-[NX3H1]~{GUAN_C});$(*-[CX4H2]-{AA})]", "Arg CG"),
    ("opls_308", 7, f"[CX4H2;$(*-[CX4H2]-[NX3H1]-[CX3](=[NX2H1])-[NX3H2]);$(*-[CX4H2]-{AA})]", "CG of neutral Arg"),
    # arginine, neutral (HN=CZ)
    ("opls_752", 7, "[CX3;$(*(=[NX2H1])(-[NX3H2])-[NX3H1])]", "CZ of neutral Arg"),
    ("opls_750", 7, "[NX2H1;$(*=[CX3](-[NX3H2])-[NX3H1])]", "N1 of neutral Arg (HN=CZ)"),
    ("opls_751", 7, "[NX3H2;$(*-[CX3](=[NX2H1])-[NX3H1])]", "N2 of neutral Arg (H2N-CZ)"),
    ("opls_749", 7, "[NX3H1;$(*-[CX4]);$(*-[CX3](=[NX2H1])-[NX3H2])]", "NE of neutral Arg"),
    ("opls_748", 7, "[CX4H2;$(*-[NX3H1]-[CX3](=[NX2H1])-[NX3H2])]", "CD of neutral Arg"),
    ("opls_304", 7, "[H][NX3H1;$(*-[CX4]);$(*-[CX3](=[NX2H1])-[NX3H2])]", "HE of neutral Arg"),
    ("opls_301", 7, "[H][#7;$(*~[CX3](~[NX2H1,NX3H2])~[NX3H1,NX3H2]);!$(*-[CX4])]", "HH of neutral Arg"),
    # histidine: neutral (HID / HIE) and charged (HIP) imidazole
    ("opls_505", 7, f"[CX4H2;$(*-{HIS_A});$(*-{AA});$(*-[#6]1~[#6]~[#7]~[#6]~[#7]~1)]", "His CB"),
    ("opls_503", 7, f"[#7X3H1;r5;{IMID};$(*1~*~*~[#7X2]~*~1),$(*1~*~[#7X2]~*~*~1)]", "N-H of neutral His"),
    ("opls_504", 7, f"[H][#7X3H1;r5;{IMID};$(*1~*~*~[#7X2]~*~1),$(*1~*~[#7X2]~*~*~1)]", "H on the N of neutral His"),
    ("opls_511", 7, f"[#7X2;r5;{IMID}]", "N of neutral His without H"),
    ("opls_506", 7, f"[#6X3;r5;{IMID};$(*(~[#7X2])~[#7X3H1])]", "CE1 of neutral His"),
    ("opls_507", 7, f"[#6X3;r5;{IMID};!$(*(~[#7])~[#7]);$(*~[#7X2])]", "ring C next to the N without H (HID CD2, HIE CG)"),
    ("opls_508", 7, f"[#6X3;r5;{IMID};!$(*(~[#7])~[#7]);$(*~[#7X3H1]);!$(*1~*~[#7X3H1]~*~[#7X3H1]~1)]", "ring C next to the N-H (HID CG, HIE CD2)"),
    ("opls_512", 8, f"[#7X3H1;r5;{IMID};$(*1~*~[#7X3H1]~*~*~1),$(*1~*~*~[#7X3H1]~*~1)]", "N of charged His"),
    ("opls_513", 8, f"[H][#7X3H1;r5;{IMID};$(*1~*~[#7X3H1]~*~*~1),$(*1~*~*~[#7X3H1]~*~1)]", "H on the N of charged His"),
    ("opls_509", 8, f"[#6X3;r5;{IMID};$(*(~[#7X3H1])~[#7X3H1])]", "CE1 of charged His"),
    ("opls_510", 8, f"[#6X3;r5;{IMID};!$(*(~[#7])~[#7]);$(*1~*~[#7X3H1]~*~[#7X3H1]~1),$(*1~[#7X3H1]~*~[#7X3H1]~*~1)]", "CG / CD2 of charged His"),
    # tryptophan's indole
    ("opls_136", 8, f"[CX4H2;$(*-[#6X3;r5;!r6]~[#6;r5;r6]);$(*-{AA})]", "Trp CB"),
    ("opls_500", 8, f"[#6X3;r5;!r6;$(*-[CX4H2]-{AA});$(*~[#6;r5;r6])]", "Trp CG"),
    ("opls_514", 8, f"[#6X3;r5;!r6;$(*~{NH_R5});$(*~{TRP_CG})]", "Trp CD1"),
    ("opls_146", 8, f"[H][#6X3;r5;!r6;$(*~{NH_R5});$(*~{TRP_CG})]", "Trp HD1"),
    ("opls_501", 8, f"[#6X3;r5;r6;!$(*~[#7]);$(*~{TRP_CG})]", "Trp CD2"),
    ("opls_502", 8, f"[#6X3;r5;r6;$(*~{NH_R5}~[#6X3;r5]~[#6X3;r5;$(*-[CX4H2]-{AA})])]", "Trp CE2"),
    ("opls_503", 8, f"[{NH_R5[1:-1]};$(*~[#6X3;r5;!r6]~[#6X3;r5;!r6;$(*-[CX4H2]-{AA})]);$(*~[#6;r5;r6])]", "Trp NE1"),
    ("opls_504", 8, f"[H][{NH_R5[1:-1]};$(*~[#6X3;r5;!r6]~[#6X3;r5;!r6;$(*-[CX4H2]-{AA})]);$(*~[#6;r5;r6])]", "Trp HE1"),
    # cystine
    ("opls_203", 6, "[SX2;$(*-[SX2]);$(*-[#6])]", "disulfide S (on carbon: an S between two S has no number)"), ("opls_214", 6, "[CX4H2;$(*-[SX2]-[SX2])]", "CH2 on a disulfide S"),
    # pyroglutamate (a gamma-lactam closing on the C alpha)
    ("opls_238", 8, "[NX3H1;r5;$(*-[CX3;r5]=[OX1])]", "lactam N-H"), ("opls_241", 8, "[H][NX3H1;r5;$(*-[CX3;r5]=[OX1])]", "lactam N-H H"),
    ("opls_235", 8, "[CX3;r5;$(*=[OX1]);$(*-[NX3;r5])]", "lactam C=O"), ("opls_236", 8, "[OX1]=[CX3;r5;$(*-[NX3;r5])]", "lactam O"),
    ("opls_224B", 9, "[CX4H1;r5;$(*-[NX3H1;r5]-[CX3;r5]=[OX1]);$(*-[CX3;!R](~[OX1])~[OX1,#7X3])]", "pyroglutamate C alpha"),
    ("opls_274", 8, "[CX4H2;r5;$(*-[CX3;r5](=[OX1])-[NX3;r5]-[CX4H1;r5]-[CX3;!R](~[OX1])~[OX1,#7X3])]", "pyroglutamate CG"),
    # DMSO: OPLS-AA's all-atom sulfoxide types 496 / 497 / 498 (the general rules), as the GROMACS OPLS-AA topologies of
    # Caleman et al. 2012 type it; aminoacids.rtp's DMSO residue ("use at your own risk", 2006) takes the united-atom DMSO
    # S and O (124 / 125) with alkane C 139 and charges of its own, which with the types' charges leave DMSO at +0.04 e
    # (no 145B for biphenyl: GROMACS's file charges it -0.115 with no H to balance, -0.23 e on biphenyl; left untyped)
    # alkynes, GROMACS's own set (ffnonbonded.itp 925-931, other meanings than the 2024 table's): terminal C 925, its H 926,
    # the inner sp C of R-C#CH by R's H count (927 2-3 H, 928 1 H, 929 none or phenyl), the H on R's carbon 930 (that C:
    # 135-139), RC#CR 931; propyne 135 + 3 x 930 + 927 + 925 + 926 = 0, as the GROMACS topology of Caleman et al. types it
    ("opls_925", 2, "[CX2H1]#[CX2;$(*-[#6])]", "alkyne terminal C of R-C#CH (acetylene, -0.21 x 2 + 0.20 x 2, has no balance)"),
    ("opls_927", 2, "[CX2;$(*#[CX2H1]);$(*-[CX4H3,CX4H2])]", "alkyne inner C, R with 2-3 H"),
    ("opls_928", 2, "[CX2;$(*#[CX2H1]);$(*-[CX4H1])]", "alkyne inner C, R with 1 H"),
    ("opls_929", 2, "[CX2;$(*#[CX2H1]);$(*-[CX4H0,c])]", "alkyne inner C, R with no H or phenyl"),
    ("opls_930", 2, "[H][CX4;$(*-[CX2]#[CX2H1])]", "H on the alkyl C of R-C#CH"),
    ("opls_931", 2, "[CX2;$(*([#6])#[CX2][#6])]", "alkyne C of RC#CR"),
    # ions (ions.itp) and TIP3P water (111 / 112); other water models are overlays
    *[(t, 9, s, "ion") for t, s in (("opls_406", "[Li+]"), ("opls_407", "[Na+]"), ("opls_408", "[K+]"), ("opls_409", "[Rb+]"),
                                    ("opls_410", "[Cs+]"), ("opls_411", "[Mg+2]"), ("opls_412", "[Ca+2]"), ("opls_413", "[Sr+2]"),
                                    ("opls_414", "[Ba+2]"), ("opls_400", "[F-]"), ("opls_401", "[Cl-]"), ("opls_402", "[Br-]"),
                                    ("opls_403", "[I-]"))],
    ("opls_111", 9, "[OX2H2]", "TIP3P water O"), ("opls_112", 9, "[H][OX2H2]", "TIP3P water H"),
]
# GROMACS's own termini (aminoacids.n.tdb / c.tdb replace the C alpha): N-terminal NH3+ 293B (Gly 292B), a zwitterion
# 299 (Gly 298), an N-terminal proline's C alpha 246, a neutral NH2 terminus 912B (Gly 906B); C-terminal ones as ADD
GMX_TERMINI = [
    ("opls_293B", 11, "[CX4H1;!R;$(*-[NX4H3]);$(*-[CX3](=[OX1])-[NX3])]", "N-terminal C alpha (NH3+)"),
    ("opls_292B", 11, "[CX4H2;!R;$(*-[NX4H3]);$(*-[CX3](=[OX1])-[NX3])]", "N-terminal glycine C alpha (NH3+)"),
    ("opls_299", 11, "[CX4H1;!R;$(*-[NX4H3]);$(*-[CX3](~[OX1])~[OX1])]", "zwitterion C alpha"),
    ("opls_298", 11, "[CX4H2;!R;$(*-[NX4H3]);$(*-[CX3](~[OX1])~[OX1])]", "zwitterion glycine C alpha"),
    ("opls_246", 11, "[CX4H1;R;$(*-[NX4H2]);$(*-[CX3]=[OX1])]", "N-terminal proline C alpha (NH2+)"),
    ("opls_912B", 11, "[CX4H1;!R;$(*-[NX3H2]);$(*-[CX3](=[OX1])-[NX3])]", "N-terminal C alpha (neutral NH2)"),
    ("opls_906B", 11, "[CX4H2;!R;$(*-[NX3H2]);$(*-[CX3](=[OX1])-[NX3])]", "N-terminal glycine C alpha (neutral NH2)"),
]

def write_fallback(used24):
    """the 2024 types the fallback rules give, with their Lennard-Jones and the bonded terms of classes only they carry,
    into oplsaa-gromacs.json (what an earlier run added is removed first)"""
    import fnmatch
    src = (f"{FALLBACK}: Jorgensen et al., J. Phys. Chem. B 128, 250 (2024) (CAPS oplsaa2024); GROMACS's oplsaa.ff has no "
           "numbers that type this group at an integer charge")
    ff = GFF
    for key in ("atom_types", "pairs", "bonds", "angles", "dihedrals", "impropers"):
        ff[key] = [x for x in ff[key] if not x.get("source", "").startswith(FALLBACK)]
    # A fallback type's classes are its 2024 class with TAG (first choice) and the 2024 class itself (second choice,
    # bond2 / angle2 / dihedral2 / improper2). A term touching a fallback atom takes a TAG rule first: the 2024 term for a
    # class tuple GROMACS's file lacks, one rule per choice of which fallback-class positions carry TAG; otherwise the
    # plain names find GROMACS's own term. A molecule of GROMACS types only never sees a TAG name, so its terms, and its
    # refusal where GROMACS has none, stay as they were.
    pair24 = {p["match"][0]: p for p in O24["pairs"] if len(p["match"]) == 1}
    fbcls = set()
    for name, t24 in sorted(used24.items()):
        o = T24[t24]
        eq = {k: v + TAG for k, v in o["equivalence"].items()}
        eq.update({k + "2": v for k, v in o["equivalence"].items()})
        ff["atom_types"].append({**{k: v for k, v in o.items() if k != "source"}, "name": name, "equivalence": eq, "source": src})
        ff["pairs"].append({**pair24[t24], "name": name, "match": [name], "source": src})
        fbcls |= set(o["equivalence"].values())
    native_cls = {v for t in NATIVE.values() for v in t["equivalence"].values()}
    only = fbcls - native_cls   # classes no GROMACS type has: their 2024 terms as they are

    def matched(rs, m):
        return any(len(r["match"]) == len(m) and (all(fnmatch.fnmatchcase(a, b) for a, b in zip(m, r["match"])) or
                                                  all(fnmatch.fnmatchcase(a, b) for a, b in zip(m[::-1], r["match"])))
                   for r in rs)
    added = {}
    for key in ("bonds", "angles", "dihedrals"):
        native, seen = list(ff[key]), set()
        for r in O24[key]:
            m = r["match"]
            if matched(native, m):
                continue
            base = {k: v for k, v in r.items() if k != "source"}
            if any(c in only for c in m):
                if all(c == "*" or c in native_cls or c in fbcls for c in m):   # every class carried by some type here
                    ff[key].append({**base, "source": src})
                    added[key] = added.get(key, 0) + 1
                continue
            if "*" in m:
                continue
            pos = [i for i, c in enumerate(m) if c in fbcls]
            for bits in range(1, 1 << len(pos)):
                v = list(m)
                for b, i in enumerate(pos):
                    if bits >> b & 1:
                        v[i] += TAG
                k = min(tuple(v), tuple(v[::-1]))
                # a plain position matches only GROMACS's own types (a fallback atom's first name carries TAG)
                if k in seen or any(c != "*" and not c.endswith(TAG) and c not in native_cls for c in v):
                    continue
                seen.add(k)
                ff[key].append({**base, "name": "-".join(v), "match": v, "source": src})
                added[key] = added.get(key, 0) + 1
    if O24.get("ring3_torsions") == "zero" and "O$" in fbcls:
        # the 2024 library has no torsion inside a three-membered ring (ring3_torsions: zero); GROMACS's file has no such
        # setting, so the epoxide's ring torsions (X-C-O-C, every one in the ring) are written as zero terms
        ff["dihedrals"].append({"name": "X-CY-O$-CY (epoxide ring)", "match": ["*", "CY" + TAG, "O$" + TAG, "CY" + TAG], "style": "opls",
                                "params": [0.0, 0.0, 0.0, 0.0],
                                "source": src + "; zero, as the 2024 library leaves out every torsion inside a three-membered ring"})
        added["dihedrals"] = added.get("dihedrals", 0) + 1
    for r in O24["impropers"]:   # the 2024 file names the centre first, the GROMACS file third
        if r["match"][0] in only:
            m = r["match"]
            ff["impropers"].append({**{k: v for k, v in r.items() if k != "source"}, "match": [m[1], m[2], m[0], m[3]], "source": src})
            added["impropers"] = added.get("impropers", 0) + 1
    note = (f"{FALLBACK} (bench/typing/make_opls_gromacs_rules.py): {len(used24)} OPLS-AA 2024 types for groups GROMACS's "
            "numbers cannot type at an integer charge (opls2024_NNN where GROMACS's number is another atom, caps_* CAPS's "
            "supplement), with their Lennard-Jones; their classes are the 2024 class with @24 first (2024 terms for class tuples "
            "GROMACS's file lacks) and the plain class second (GROMACS's own terms); each carries its source")
    ff["notes"] = [n for n in ff["notes"] if not n.startswith(FALLBACK)] + [note]
    json.dump(ff, open(FFPATH, "w"), indent=1, ensure_ascii=False)
    print(f"oplsaa-gromacs.json: {len(used24)} fallback types, classes only they carry {sorted(only)}, terms added {added}")


if __name__ == "__main__":
    for t, p, sm, d in ADD + GMX_TERMINI:   # above the general rules' priorities (0-8)
        rules.append({"type": t, "priority": p + 5, "smarts": sm, "description": d})
    # SPC only: a TIP4P overlay (opls_113 O at 0, opls_114 H at +0.52) has no M site to carry the -1.04, so each water was
    # +1.04 e; TIP4P water comes from the tip4p-* water libraries, which add their M site
    for name, (o, h, label) in {"spc": ("opls_116", "opls_117", "SPC")}.items():
        json.dump({"format": "caps-typing", "version": 1, "forcefield": f"OPLS-AA (GROMACS) {label} water",
                   "description": f"{label} water for oplsaa-gromacs (an overlay on its rules)",
                   "rules": [{"type": o, "smarts": "[OX2H2]", "priority": 20, "description": f"{label} water O"},
                             {"type": h, "smarts": "[H][OX2H2]", "priority": 20, "description": f"{label} water H"}]},
                  open(os.path.join(OUT, f"oplsaa-gromacs-{name}.typing.json"), "w"), indent=1)
    doc = {"format": "caps-typing", "version": 1, "forcefield": "OPLS-AA (GROMACS oplsaa.ff)", "unknown_types": "untyped",
           "description": "CAPS's OPLS-AA rules in GROMACS's opls_NNN names (bench/typing/make_opls_gromacs_rules.py), with the "
                          "OPLS-AA 2024 library's types (opls2024_NNN, caps_*) for groups GROMACS's numbers cannot balance"
                          + ("; types neither file has (" + ", ".join(sorted(dropped)) + ") are left untyped" if dropped else ""),
           "rules": rules,
           "refused": [{"smarts": sm, "reason": why + " (neither GROMACS's oplsaa.ff nor the OPLS-AA 2024 table balances it)"}
                       for sm, why in _m["REFUSED"]]}
    json.dump(doc, open(os.path.join(OUT, "oplsaa-gromacs.typing.json"), "w"), indent=1, ensure_ascii=False)
    write_fallback(used24)
    print(f"oplsaa-gromacs: {len(rules)} rules, {len(dropped)} types the GROMACS file lacks")

