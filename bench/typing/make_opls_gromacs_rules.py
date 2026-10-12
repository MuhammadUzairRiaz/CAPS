#!/usr/bin/env python3
"""CAPS typing rules for OPLS-AA built from GROMACS's oplsaa.ff (bench/ff/convert_gromacs_opls.py: oplsaa-gromacs):
CAPS's OPLS-AA chemistry rules (data/typing/oplsaa2024-moltemplate.typing.json, which name Jorgensen's type numbers)
written as GROMACS's opls_NNN, keeping the rules whose type the GROMACS file has, plus the residue types of GROMACS's
own amino-acid topologies (checked with bench/ff/validate_rtp_types.py against aminoacids.rtp).

usage: make_opls_gromacs_rules.py [OUT_DIR]
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, "typing")
have = {t["name"] for t in json.load(open(os.path.join(DATA, "forcefields", "oplsaa-gromacs.json")))["atom_types"]}
src = json.load(open(os.path.join(DATA, "typing", "oplsaa2024-moltemplate.typing.json")))


def opls(t):
    return "opls_%03d" % int(t) if re.fullmatch(r"\d+", t) else t


# every rule with its flags (make_oplsaa_rules.py's R): the moltemplate file keeps only the numbers that library has
_m = {}
exec(open(os.path.join(HERE, "make_oplsaa_rules.py")).read().split('\nif __name__ == "__main__":')[0], _m)
ALL = _m["R"]
rules, dropped = [], set()
for r in [x for x in ALL if not x.get("only2024")]:
    t = opls(r["type"])
    if t not in have:
        dropped.add(t)
        continue
    rules.append({**r, "type": t, **({"overrides": [opls(x) for x in r["overrides"]]} if "overrides" in r else {})})
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
    ("opls_203", 6, "[SX2;$(*-[SX2])]", "disulfide S"), ("opls_214", 6, "[CX4H2;$(*-[SX2]-[SX2])]", "CH2 on a disulfide S"),
    # pyroglutamate (a gamma-lactam closing on the C alpha)
    ("opls_238", 8, "[NX3H1;r5;$(*-[CX3;r5]=[OX1])]", "lactam N-H"), ("opls_241", 8, "[H][NX3H1;r5;$(*-[CX3;r5]=[OX1])]", "lactam N-H H"),
    ("opls_235", 8, "[CX3;r5;$(*=[OX1]);$(*-[NX3;r5])]", "lactam C=O"), ("opls_236", 8, "[OX1]=[CX3;r5;$(*-[NX3;r5])]", "lactam O"),
    ("opls_224B", 9, "[CX4H1;r5;$(*-[NX3H1;r5]-[CX3;r5]=[OX1]);$(*-[CX3;!R](~[OX1])~[OX1,#7X3])]", "pyroglutamate C alpha"),
    ("opls_274", 8, "[CX4H2;r5;$(*-[CX3;r5](=[OX1])-[NX3;r5]-[CX4H1;r5]-[CX3;!R](~[OX1])~[OX1,#7X3])]", "pyroglutamate CG"),
    # DMSO as GROMACS's topology has it
    ("opls_124", 8, "[SX3](=[OX1])([CX4H3])[CX4H3]", "DMSO S"), ("opls_125", 8, "[OX1]=[SX3]([CX4H3])[CX4H3]", "DMSO O"),
    ("opls_139", 8, "[CX4H3][SX3](=[OX1])[CX4H3]", "DMSO C"),
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
if __name__ == "__main__":
    for t, p, sm, d in ADD + GMX_TERMINI:   # above the general rules' priorities (0-8)
        rules.append({"type": t, "priority": p + 5, "smarts": sm, "description": d})
    for name, (o, h, label) in {"spc": ("opls_116", "opls_117", "SPC"), "tip4p": ("opls_113", "opls_114", "TIP4P")}.items():
        json.dump({"format": "caps-typing", "version": 1, "forcefield": f"OPLS-AA (GROMACS) {label} water",
                   "description": f"{label} water for oplsaa-gromacs (an overlay on its rules)",
                   "rules": [{"type": o, "smarts": "[OX2H2]", "priority": 20, "description": f"{label} water O"},
                             {"type": h, "smarts": "[H][OX2H2]", "priority": 20, "description": f"{label} water H"}]},
                  open(os.path.join(OUT, f"oplsaa-gromacs-{name}.typing.json"), "w"), indent=1)
    doc = {"format": "caps-typing", "version": 1, "forcefield": "OPLS-AA (GROMACS oplsaa.ff)", "unknown_types": "untyped",
           "description": "CAPS's OPLS-AA rules in GROMACS's opls_NNN names (bench/typing/make_opls_gromacs_rules.py); types the "
                          "GROMACS file lacks (" + ", ".join(sorted(dropped)) + ") are left untyped",
           "rules": rules}
    json.dump(doc, open(os.path.join(OUT, "oplsaa-gromacs.typing.json"), "w"), indent=1, ensure_ascii=False)
    print(f"oplsaa-gromacs: {len(rules)} rules, {len(dropped)} types the GROMACS file lacks")

    # OPLS-AA 2024 from Jorgensen's own files (bench/ff/convert_boss_opls.py: oplsaa2024): the same chemistry rules in
    # opls_NNN names, water TIP3P as the .par numbers it (111 / 112)
    have24 = {t["name"] for t in json.load(open(os.path.join(DATA, "forcefields", "oplsaa2024.json")))["atom_types"]}
    WATER = {"9999": "opls_111", "9998": "opls_112"}
    rules24, lacking = [], set()
    for r in ALL:
        t = WATER.get(r["type"], opls(r["type"]))
        if t not in have24:
            lacking.add(t)
            continue
        r = {k: v for k, v in r.items() if k != "only2024"}
        rules24.append({**r, "type": t, **({"overrides": [WATER.get(x, opls(x)) for x in r["overrides"]]} if "overrides" in r else {})})
    json.dump({"format": "caps-typing", "version": 1, "forcefield": "OPLS-AA 2024 (Jorgensen)", "unknown_types": "untyped",
               "description": "CAPS's OPLS-AA rules in opls_NNN names for the library built from Jorgensen's 2024 BOSS files "
                              "(bench/typing/make_opls_gromacs_rules.py)", "rules": rules24},
              open(os.path.join(OUT, "oplsaa2024.typing.json"), "w"), indent=1, ensure_ascii=False)
    print(f"oplsaa2024: {len(rules24)} rules" + (f", types the file lacks: {sorted(lacking)}" if lacking else ""))
