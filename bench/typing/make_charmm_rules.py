#!/usr/bin/env python3
"""CAPS typing rules for the CHARMM libraries built from CHARMM's own files (bench/ff/open/build_charmm.sh: *-toppar ids),
checked against CHARMM's own topologies (top_all36_*.rtf, top_all22_prot.rtf; bench/ff/validate_rtf_types.py --cap-links).

The rules start from CAPS's chemistry rules for each family (make_biomolecular_rules.py) and add what CHARMM's residues
use: CT2A for the C beta of Asp / Glu / protonated His, proline's CP1 / CP2 / CP3, tryptophan's CAI, the protonated
histidine ring (NR3, HR1, HR2); water is TIP3P as toppar_water_ions.str names it (OT / HT).

usage: make_charmm_rules.py [OUT_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_biomolecular_rules as B

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data", "typing")
AA_CA = "[CX4]([NX3,NX4])[CX3]=[OX1]"          # an amino acid's C alpha
IMIDAZOLIUM_RING = "1~[#7H1]~[#6]~[#7H1]~[#6]~1"   # both ring N carrying H (HSP)


def protein(v36):
    R = B.charmm_prot(v36)
    if v36:
        # C beta of Asp, Glu and protonated His (CHARMM36's chi1 / chi2 refit)
        R.add("CT2A", 4, f"[CX4H2;$(*{AA_CA});$(*[CX3](~[OX1])~[OX1]),$(*[CX4H2][CX3](~[OX1])~[OX1]),$(*-[#6;r5]{IMIDAZOLIUM_RING})]",
              "C beta of Asp, Glu, protonated His")
        # tryptophan: the benzene CH next to the ring fusion
        R.add("CAI", 5, "[cH1;r6;$(c:[c;r5;r6])]", "indole C4 / C7 (Trp CE3, CZ2)")
    # proline ring: C alpha CP1, C beta / gamma CP2, C delta CP3
    PRO_N = "[NX3;r5;$(N[CX3]=O),$(N(C)(C)C)]"
    R.add("CP1", 6, f"[CX4H1;r5;$(*{PRO_N});$(*[CX3]=O)]", "proline C alpha")
    R.add("CP3", 6, f"[CX4H2;r5;$(*{PRO_N})]", "proline C delta")
    R.add("CP2", 6, f"[CX4H2;r5;!$(*{PRO_N});$(*1@*@*@{PRO_N}@*1),$(*1@*@{PRO_N}@*@*1)]", "proline C beta / gamma")
    # protonated histidine: both ring N NR3, the H between them HR2, the other ring H HR1
    R.add("NR3", 7, f"[#7H1;r5;$(*{IMIDAZOLIUM_RING}),$(*1~[#6]~[#7H1]~[#6]~[#6]~1)]", "imidazolium N (HSP)")
    R.add("HR2", 7, f"[H][#6;r5;$(*(~[#7H1])~[#7H1])]", "H between the two NH of an imidazolium (HSP HE1)")
    R.add("HR1", 7, f"[H][#6;r5;$(*1~[#7H1]~[#6]~[#7H1]~[#6]~1);!$(*(~[#7H1])~[#7H1])]", "imidazolium C4 / C5 H (HSP HD2)")
    return R


def lipid():
    """CHARMM36 lipids as CHARMM's own topology types them: H on sp3 C by the carbon's H count (HAL1 / 2 / 3), every choline
    methyl CTL5 and every H on a carbon bonded to the quaternary N HL (the library-template conventions of the old rules,
    all H HAL3 and CTL5 in DMPC only, are not CHARMM's)"""
    R = B.Rules([r for r in B.charmm_lipid() if r["type"] not in ("HAL3", "CTL5", "HL")])
    R.add("HAL1", 0, "[H][CX4H1]", "H on CH"); R.add("HAL2", 0, "[H][CX4H2]", "H on CH2"); R.add("HAL3", 0, "[H][CX4H3]", "H on CH3")
    R.add("CTL5", 2, "[CX4H3][NX4]", "methyl on a quaternary N (choline)")
    R.add("HL", 2, "[H][CX4][NX4H0]", "H on a carbon bonded to a quaternary N")
    return R


def carb():
    """CHARMM36 carbohydrates as CHARMM's own topology types them (top_all36_carb.rtf), beyond the ring rules"""
    R = B.Rules([r for r in B.charmm_carb() if r["type"] != "CC3263"])
    RO6 = "[OX2;r6]"
    R.add("CC3261", 1, f"[CX4H2;r6;!$(*@{RO6})]", "pyranose ring CH2 (deoxy C2-C4)")
    R.add("CC3263", 2, f"[CX4H2;r6;$(*@{RO6})]", "pyranose ring CH2 next to the ring O (pentose C5)")
    R.add("CC3062", 4, f"[CX4H0;r6;$(*@{RO6});$(*-[OX2;!R])]", "pyranose anomeric C without H (ketose, sialic acid)")
    R.add("CC3051", 4, "[CX4H0;r5;$(*@[OX2;r5]);$(*-[OX2;!R])]", "furanose anomeric C without H (fructose)")
    # open chains (alditols, open aldoses and ketoses, glycerol side chains)
    R.add("CC312", 2, "[CX4H1;!R;$(*[OX2])]", "open-chain CH bearing O")
    R.add("CC322", 2, "[CX4H2;!R;$(*[OX2]);!$(*[CX4;R])]", "open-chain CH2 bearing O (the exocyclic C6 on a ring stays CC321)")
    R.add("CC2O4", 4, "[CX3H1](=[OX1])[#6]", "aldehyde C (open aldose)"); R.add("OC2D4", 4, "[OX1]=[CX3H1][#6]", "aldehyde O")
    R.add("HCR1", 4, "[H][CX3]=[OX1]", "aldehyde H")
    R.add("CC2O3", 4, "[CX3H0](=[OX1])([#6])[#6]", "ketone C (open ketose)"); R.add("OC2D3", 4, "[OX1]=[CX3H0]([#6])[#6]", "ketone O")
    R.add("CC2O5", 4, "[CX3](=[OX1])[OX2][#6]", "ester carbonyl C (O-acetyl)"); R.add("OC2D1", 3, "[OX1]=[CX3][OX2][#6]", "ester carbonyl O")
    R.add("OC301", 3, "[OX2]([#6])[CX3]=[OX1]", "ester O")
    R.add("HCP1", 2, "[H][NX3][CX3]=[OX1]", "amide H (N-acetyl)")
    return R


def water(R):
    """TIP3P as CHARMM's own water stream names it (OT / HT)"""
    return B.Rules([r for r in R if r["smarts"] not in ("[OX2H2]", "[H][OX2H2]")] +
                   [{"type": "OT", "smarts": "[OX2H2]", "priority": 8, "description": "TIP3P water O"},
                    {"type": "HT", "smarts": "[H][OX2H2]", "priority": 8, "description": "TIP3P water H"}])


SETS = {
    "charmm36-prot-toppar": (water(protein(True)), "CHARMM36 proteins (toppar)"),
    "charmm22-prot-toppar": (water(protein(False)), "CHARMM22 proteins (toppar)"),
    "charmm36-nucl-toppar": (water(B.charmm_nucl()), "CHARMM36 nucleic acids (toppar)"),
    "charmm36-lipid-toppar": (water(lipid()), "CHARMM36 lipids (toppar)"),
    "charmm36-carb-toppar": (water(carb()), "CHARMM36 carbohydrates (toppar)"),
}

if __name__ == "__main__":
    for fid, (rules, label) in SETS.items():
        doc = {"format": "caps-typing", "version": 1, "forcefield": label, "unknown_types": "untyped",
               "description": f"CAPS rules for {label}, checked against CHARMM's own topology (bench/typing/make_charmm_rules.py).",
               "rules": list(rules)}
        if fid.replace("-toppar", "") in B.ANALOGIES:
            doc["analogies"] = B.ANALOGIES[fid.replace("-toppar", "")]
        json.dump(doc, open(os.path.join(OUT, fid + ".typing.json"), "w"), ensure_ascii=False, indent=1)
        print(f"{fid}: {len(rules)} rules")
