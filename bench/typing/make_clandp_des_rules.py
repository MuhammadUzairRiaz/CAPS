#!/usr/bin/env python3
"""CAPS typing rules for CL&P's deep-eutectic-solvent set, built from Padua's own des.ff (bench/ff/convert_fftool.py:
clandp-des), in des.ff's type names, checked against CL&P's molecule files (bench/ff/validate_fftool_types.py: the
cholinium file des.ff ships with, and every ion file of the CL&P repository whose types des.ff also defines).

des.ff gives types, not molecules: each block is one family of DES component, and the charges of a block add up to the
component's charge only in one way, which these rules follow.

  cations    cholinium (N4, C1, C1A, COL, H1, HC, OH, HO, as cholinium.xyz labels it), tetraalkylammonium (N4, C1, C1A,
             C2, CE, CS, CT, H1, HC), protonated tertiary amines (N3, H3), hydroxyethylammonium (COL on the CH2 holding
             the OH), tetraalkylphosphonium (P3, C1P, C2 and the alkane types), betaine (N4B: its trimethylammonium head
             is +1 like N1111's with N4B for N4, so its CH2 COO- is a carboxylate's CH2, C2C with HC)
  anions     chloride, bromide (OPLS's), carboxylates (CO2, O2, CTA, C2C, C3C)
  HBDs       urea (CU, OU, NU, HU); alcohols (OPLS-AA mono alcohols: CTO for CH2-OH, C2O for CH-OH, OH, HO, H on C HC);
             diols and triols (OPLS-AA's: OHG / HOG in a diol, OHT / HOT in a 1,2,3-triol, CTO on CH2-OH, CSO on CH-OH,
             H1O on those carbons); carboxylic acids (CO2 for the acid C, O3, OOH, HOO; des.ff's own CXT for an acetic
             CH3 and CX for a CH2 next to COOH, which make the COOH group's +0.18 neutral); ketones (CK, OK, the alpha
             CH3 / CH2 as alkane CT / CS: acetone, levulinic acid)

Molecules des.ff cannot make neutral are left untyped (the field reports them incomplete): methanol (CTO with three HC is
+0.06; des.ff's own H1O at 0.04 is commented out), malonic and oxalic acid (+0.18, +0.36), alpha-hydroxy and aromatic
acids, a CH next to COOH, tertiary alcohols, primary and secondary ammonium. OHC / HOC (an alternative cholinium
hydroxyl with the same sum as OH / HO) and CQ, C1B (no block uses them) are defined in the library and never assigned:
Padua's cholinium file labels its hydroxyl OH / HO.

usage: make_clandp_des_rules.py [OUT_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data", "typing")
R = []


def add(t, p, s, d, **kw):
    R.append({"type": t, "smarts": s, "priority": p, "description": d, **kw})


HET = "[!#6;!#1]"
AMN, PHN = "[#7X4]", "[#15X4;!$(*~[!#6])]"
CATN = f"[$({AMN}),$({PHN})]"                     # the ammonium N or phosphonium P
COO = "[CX3](~[OX1])~[OX1]"                       # a carboxylate
COOH = "[CX3](=[OX1])-[OX2H1]"                    # a carboxylic acid
ALC = "[OX2H1;$(*-[CX4])]"                        # an alcohol O
# an alcohol O with another one two to four carbons along the chain (a diol or a triol)
POLY = ("[OX2H1;$(*-[CX4]-[CX4]-[OX2H1]),$(*-[CX4]-[CX4]-[CX4]-[OX2H1]),$(*-[CX4]-[CX4]-[CX4]-[CX4]-[OX2H1])]")
TRIOL = "[CX4](-[OX2H1])-[CX4](-[OX2H1])-[CX4]-[OX2H1]"
NOT_CO = "!$(*-[CX3]=[OX1])"                      # not on a carbonyl or carboxyl C
# ---------------------------------------------------------------- alkyl chains (OPLS-AA), lowest priority
add("CT", 0, f"[CX4H3;!$(*~{HET});!$(*-{COOH})]", "alkane CH3 (also on a ketone C)")
add("CS", 0, f"[CX4H2;!$(*~{HET});!$(*-{COOH})]", "alkane CH2 (also next to a ketone C)")
add("HC", 0, "[H][CX4]", "H on an alkane C")
# ---------------------------------------------------------------- ammonium and phosphonium heads
add("N4", 3, "[#7X4H0]", "quaternary ammonium N")
add("N3", 3, "[#7X4H1]", "protonated tertiary amine N (R3NH+)")
add("H3", 3, "[H][#7X4H1]", "H on R3NH+")
add("C1", 3, f"[CX4H3;$(*-{AMN})]", "methyl on the ammonium N")
add("C1A", 3, f"[CX4H2;$(*-{AMN})]", "CH2 on the ammonium N")
add("C2", 2, f"[CX4H2;$(*-[CX4H2]-{CATN});!$(*~{HET})]", "CH2 next to the cation's first carbon")
add("CE", 3, f"[CX4H3;$(*-[CX4H2]-{AMN})]", "an ethyl's CH3 on the ammonium N")
add("H1", 3, f"[H][CX4;$(*-{CATN})]", "H on a carbon bonded to the ammonium N or phosphonium P")
add("P3", 4, PHN, "phosphonium P")
add("C1P", 4, f"[CX4;$(*-{PHN})]", "carbon on the phosphonium P")
# cholinium and hydroxyethylammonium: the CH2 holding the OH
add("COL", 5, f"[CX4H2;$(*-[OX2H1]);$(*-[CX4]-{AMN})]", "CH2-OH of cholinium (next to the ammonium's first carbon)")
# betaine: N4B, its CH2 a carboxylate's (C2C with HC)
add("N4B", 5, f"[#7X4H0;$(*-[CX4H2]-{COO})]", "betaine N")
add("HC", 6, f"[H][CX4;$(*-{COO})]", "H on a carbon bonded to a carboxylate C (betaine's CH2 included)")
# ---------------------------------------------------------------- halides
add("Cl", 5, "[Cl;X0]", "chloride")
add("Br", 5, "[Br;X0]", "bromide")
# ---------------------------------------------------------------- carboxylates
add("CO2", 5, f"[CX3;$(*(~[OX1])~[OX1])]", "carboxylate C")
add("O2", 5, "[OX1;$(*~[CX3]~[OX1])]", "carboxylate O")
add("CTA", 5, f"[CX4H3;$(*-{COO})]", "acetate CH3")
add("C2C", 6, f"[CX4H2;$(*-{COO})]", "CH2 on a carboxylate (betaine's CH2 included)")
add("C3C", 5, f"[CX4H1;$(*-{COO});!$(*~{HET})]", "CH on a carboxylate")
# ---------------------------------------------------------------- urea
add("CU", 5, f"[CX3;$(*(=[OX1])(-[NX3H2])-[NX3H2])]", "urea C")
add("OU", 5, "[OX1;$(*=[CX3](-[NX3H2])-[NX3H2])]", "urea O")
add("NU", 5, "[NX3H2;$(*-[CX3](=[OX1])-[NX3H2])]", "urea N")
add("HU", 5, "[H][NX3H2;$(*-[CX3](=[OX1])-[NX3H2])]", "urea H")
# ---------------------------------------------------------------- alcohols, diols, triols (OPLS-AA)
add("OH", 4, ALC, "alcohol O (mono alcohols, cholinium)")
add("HO", 4, f"[H]{ALC}", "alcohol H")
add("CTO", 4, f"[CX4H2;$(*-[OX2H1]);{NOT_CO};!$(*-[CX4]-{AMN});!$(*-{AMN})]", "CH2-OH")
add("C2O", 4, f"[CX4H1;$(*-[OX2H1]);{NOT_CO};!$(*~[!#6;!#1;!$([OX2H1])])]", "CH-OH of a mono alcohol")
add("OHG", 6, POLY, "diol O", not_in_molecule=[TRIOL])
add("HOG", 6, f"[H]{POLY}", "diol H", not_in_molecule=[TRIOL])
add("OHT", 7, ALC, "triol O", in_molecule=[TRIOL])
add("HOT", 7, f"[H]{ALC}", "triol H", in_molecule=[TRIOL])
add("CSO", 6, f"[CX4H1;$(*-{POLY});{NOT_CO};!$(*~[!#6;!#1;!$([OX2H1])])]", "CH-OH of a diol or triol")
add("H1O", 6, f"[H][CX4;$(*-{POLY});{NOT_CO}]", "H on a carbon holding a diol / triol OH")
# ---------------------------------------------------------------- carboxylic acids
add("CO2", 6, f"[CX3;$(*(=[OX1])(-[OX2H1])-[CX4])]", "carboxylic acid C (on an sp3 C)")
add("O3", 6, "[OX1;$(*=[CX3]-[OX2H1])]", "carboxylic acid =O")
add("OOH", 6, "[OX2H1;$(*-[CX3]=[OX1])]", "carboxylic acid OH")
add("HOO", 6, "[H][OX2H1;$(*-[CX3]=[OX1])]", "carboxylic acid H")
add("CXT", 6, f"[CX4H3;$(*-{COOH})]", "acetic acid CH3")
add("CX", 6, f"[CX4H2;$(*-{COOH});!$(*(-{COOH})-{COOH});!$(*~[!#6;!#1])]", "CH2 next to a COOH (not between two)")
# ---------------------------------------------------------------- ketones
add("CK", 5, "[CX3;$(*=[OX1]);$(*(-[CX4])-[CX4])]", "ketone C")
add("OK", 5, "[OX1;$(*=[CX3](-[CX4])-[CX4])]", "ketone O")
json.dump({"format": "caps-typing", "version": 1, "forcefield": "CL&P DES (deep-eutectic solvents)", "unknown_types": "untyped",
           "description": "CAPS rules for des.ff's types by DES component (bench/typing/make_clandp_des_rules.py), checked against "
                          "CL&P's own molecule files",
           "rules": R}, open(os.path.join(OUT, "clandp-des.typing.json"), "w"), indent=1, ensure_ascii=False)
print(f"clandp-des: {len(R)} rules")
