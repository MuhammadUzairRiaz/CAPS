#!/usr/bin/env python3
"""CAPS typing rules for the AMBER, CHARMM, GROMOS and TraPPE files of the library, written from each force field's own
type definitions (the meaning of every type name in its parameter and topology files):

  amber            Cornell et al., J. Am. Chem. Soc. 117, 5179 (1995): CT, HC / H1 / H2 by electron-withdrawing neighbours,
                   CA / HA, C / O, N / H, OS, OH / HO, N3 / N2, CY / NY, S, P; GLYCAM's alkene Ck / Ha
  charmm           the library's generic CHARMM set (CT1-3, CA, CE1/CE2, CD/CC amide / acid, OB, OS, ether CAE/OET ...)
  charmm22-prot,   MacKerell et al., J. Phys. Chem. B 102, 3586 (1998); Best et al., J. Chem. Theory Comput. 8, 3257 (2012):
  charmm36-prot    the protein types used chemically (CT1-3, HA*, C / O / NH1 / H amides, CD / OB / OS esters, CA / HP ...)
  charmm36-lipid   Klauda et al., J. Phys. Chem. B 114, 7830 (2010): CTL1-3 / HAL1-3, CEL1-2 / HEL1-2, CL / OBL / OSL esters
  charmm36-carb    Guvench et al., J. Chem. Theory Comput. 7, 3162 (2011): pyranose ring carbons by position, OC311 / HCP1 ...
  charmm19         Neria, Fischer & Karplus (1996), united atom: CH1E / CH2E / CH3E / CR1E, polar hydrogens explicit
  gromos-54a7      Schmid et al., Eur. Biophys. J. 40, 843 (2011), united atom: CH1 / CH2 / CH3 / CR1 / C, OA / OE / O / OM ...
  trappe-ua        Martin & Siepmann (1998) and later TraPPE-UA: CH4 / CH3 / CH2 / CH / C, ring CH2, alkenes, benzene CH / C,
                   alcohol, ether, ketone, thiol / sulfide, nitrile
  trappe1998-moltemplate   the same united-atom alkane rules (CH2, CH3, CH4)
  trappe-eh        Rai & Siepmann (2007), explicit-hydrogen aromatics: benzene, phenol, pyridine, pyrrole, furan, thiophene

United-atom files set "united_atom": the structure's hydrogens on carbon fold into their carbons before typing (CAPS's
prepare_for_forcefield); the rules then see CH2 and friends as carbons with implicit hydrogens. Where a file has no
type for an environment (an ether in CHARMM36 proteins, a quaternary carbon in the lipid set) its atoms are left
untyped and reported: "unknown_types": "untyped" also leaves atoms untyped when a rule names a type a file lacks.

usage: make_biomolecular_rules.py [DATA_DIR]
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data")

EWG = "N,O,S,F,Cl,Br,I"
ONE_EWG = f"$(*[{EWG}]);!$(*([{EWG}])[{EWG}])"
TWO_EWG = f"$(*([{EWG}])[{EWG}])"
ESTER_O = "[OX2]([#6])[CX3]=[OX1]"


class Rules(list):
    def add(self, t, prio, smarts, desc):
        self.append({"type": t, "smarts": smarts, "priority": prio, "description": desc})


def water_ions(R, o="OW", h="HW"):
    R.add(o, 8, "[OX2H2]", "water O")
    R.add(h, 8, "[H][OX2H2]", "water H")


# ---------------------------------------------------------------- AMBER (library keys)
def amber():
    R = Rules()
    R.add("C1", 0, "[CX4]", "sp3 carbon (CT)")
    R.add("CA", 1, "c", "aromatic carbon")
    R.add("G_Ck", 1, "[CX3]=[CX3]", "alkene sp2 carbon (GLYCAM Ck)")
    R.add("C", 2, "[CX3]=[OX1]", "carbonyl carbon (amide, ester, acid, ketone)")
    R.add("CY", 2, "[CX2]#[NX1]", "nitrile carbon")
    R.add("HC", 0, "[H][CX4]", "H on sp3 C without electron-withdrawing neighbours")
    R.add("H1", 1, f"[H][CX4;{ONE_EWG}]", "H on sp3 C with one electron-withdrawing neighbour (N, O, S, halogen)")
    R.add("H2", 2, f"[H][CX4;{TWO_EWG}]", "H on sp3 C with two electron-withdrawing neighbours")
    R.add("HP", 3, "[H][CX4][N+,n+]", "H on C next to a positively charged group")
    R.add("HA", 0, "[H]c", "aromatic H")
    R.add("G_Ha", 1, "[H][CX3]=[CX3]", "H on an alkene carbon (GLYCAM Ha)")
    R.add("HN", 0, "[H][N,n]", "H on N")
    R.add("HO", 0, "[H][OX2]", "hydroxyl H")
    R.add("HS", 0, "[H][SX2]", "thiol H")
    R.add("N1", 2, "[NX3][CX3]=[OX1]", "sp2 amide nitrogen")
    R.add("N4", 0, "[NX3;!$(N[CX3]=[OX1]);!$(Nc)]", "sp3 amino nitrogen")
    R.add("N2", 1, "[NX3;!$(N[CX3]=[OX1])]c", "sp2 nitrogen of an amino group on an aromatic ring")
    R.add("N3", 3, "[NX4+]", "sp3 ammonium nitrogen")
    R.add("NY", 2, "[NX1]#[CX2]", "nitrile nitrogen")
    R.add("N5", 2, "[nX2;r5]", "sp2 nitrogen in a five-membered ring with a lone pair (NB)")
    R.add("N5H", 2, "[nX3;r5]", "sp2 nitrogen in a five-membered ring with H or a substituent (NA)")
    R.add("O", 1, "[OX1]=[#6]", "carbonyl oxygen")
    R.add("OC", 3, "[OX1]~[CX3]~[OX1]", "carboxylate oxygen")
    R.add("OS", 0, "[OX2]([#6])[#6]", "ether or ester oxygen")
    R.add("OH", 1, "[OX2H1]", "hydroxyl oxygen")
    R.add("S1", 1, "[SX2][SX2]", "disulfide sulfur")
    R.add("S2", 1, "[SX2H1]", "thiol sulfur")
    R.add("P", 0, "[PX4]", "phosphate phosphorus")
    water_ions(R)
    R.add("Na+", 0, "[Na]", "sodium ion"); R.add("Cl-", 0, "[Cl-;X0]", "chloride ion"); R.add("Ca2+", 0, "[Ca]", "calcium ion")
    return R


# ---------------------------------------------------------------- CHARMM (library's generic keys)
def charmm():
    R = Rules()
    R.add("C1", 0, "[CX4H1]", "sp3 CH (CT1)"); R.add("C2", 0, "[CX4H2]", "sp3 CH2 (CT2)"); R.add("C3", 0, "[CX4H3]", "sp3 CH3 (CT3)")
    R.add("CAE1", 2, "[CX4H1][OX2;!$(O[CX3]=O);!$([OH1])]", "sp3 CH next to an ether oxygen")
    R.add("CAE2", 2, "[CX4H2][OX2;!$(O[CX3]=O);!$([OH1])]", "sp3 CH2 next to an ether oxygen")
    R.add("CAE3", 2, "[CX4H3][OX2;!$(O[CX3]=O);!$([OH1])]", "sp3 CH3 next to an ether oxygen")
    R.add("CTHF", 3, "[CX4H2;r5]@[OX2;r5]", "ring CH2 next to O in tetrahydrofuran")
    R.add("CTHP", 3, "[CX4H2;r6]@[OX2;r6]", "ring CH2 next to O in tetrahydropyran")
    R.add("CYCP", 1, "[CX4H2;r5;!$(*[O,N,S])]", "CH2 in a five-membered carbon ring")
    R.add("C4", 1, "c", "aromatic carbon")
    R.add("C15", 1, "[CX3H2]=[CX3]", "terminal alkene carbon H2C=")
    R.add("C16", 1, "[CX3;H1,H0]=[CX3]", "inner alkene carbon =CH- (and =CR-)")
    R.add("C6", 3, "[CX3](=[OX1])[NX3]", "amide carbonyl carbon")
    R.add("C5", 3, "[CX3](=[OX1])[OX2]", "carbonyl carbon of an acid or ester")
    R.add("C18", 4, "[CX3](=[OX1])([OX2])[OX2]", "carbonate carbon")
    R.add("H1", 0, "[H][CX4]", "aliphatic H")
    R.add("HAE1", 2, "[H][CX4H1][OX2;!$(O[CX3]=O);!$([OH1])]", "H on an ether CH")
    R.add("HAE2", 2, "[H][CX4H2][OX2;!$(O[CX3]=O);!$([OH1])]", "H on an ether CH2")
    R.add("HAE3", 2, "[H][CX4H3][OX2;!$(O[CX3]=O);!$([OH1])]", "H on an ether CH3")
    R.add("HTHF", 3, "[H][CX4H2;r5]", "H on a CH2 in a five-membered ring")
    R.add("H2", 0, "[H]c", "aromatic H")
    R.add("H7", 1, "[H][CX3H2]=[CX3]", "H on a terminal alkene carbon")
    R.add("H8", 1, "[H][CX3H1]=[CX3]", "H on an inner alkene carbon")
    R.add("HN", 0, "[H][N,n]", "H on N")
    R.add("HO", 0, "[H][OX2]", "hydroxyl H")
    R.add("HS", 0, "[H][SX2]", "thiol H")
    R.add("N1", 2, "[NX3H2][CX3]=[OX1]", "amide N with two H")
    R.add("N2", 2, "[NX3H1][CX3]=[OX1]", "amide N with one H")
    R.add("N5", 2, "[NX3H0][CX3]=[OX1]", "amide N without H (as in proline)")
    R.add("N3", 3, "[NX4+]", "ammonium nitrogen")
    R.add("OH", 1, "[OX2H1]", "hydroxyl oxygen")
    R.add("OET", 0, "[OX2]([#6])[#6]", "ether oxygen")
    R.add("OE", 2, ESTER_O, "ester oxygen")
    R.add("OC1", 2, "[OX1]=[CX3][OX2]", "carbonyl oxygen of an acid or ester")
    R.add("OC2", 2, "[OX1]=[CX3][NX3]", "carbonyl oxygen of an amide")
    R.add("OC3", 3, "[OX1]~[CX3]~[OX1]", "carboxylate oxygen")
    R.add("OCA", 3, "[OX1]=[CX3]([OX2])[OX2]", "carbonate carbonyl oxygen")
    R.add("S1", 1, "[SX2H1]", "thiol sulfur"); R.add("S2", 1, "[SX2][SX2]", "disulfide sulfur")
    R.add("S4", 1, "[SX3]=[OX1]", "sulfoxide sulfur"); R.add("S5", 1, "[SX4](=[OX1])=[OX1]", "sulfone / sulfate sulfur")
    R.add("P1", 0, "[PX4]", "phosphate phosphorus"); R.add("OP", 1, "[OX1]~P", "phosphate oxygen")
    water_ions(R, "OW", "HW")
    R.add("Na+", 0, "[Na]", "sodium ion"); R.add("Cl-", 0, "[Cl-;X0]", "chloride ion")
    return R


# ---------------------------------------------------------------- CHARMM proteins (standard names)
def charmm_prot(v36):
    R = Rules()
    R.add("CT1", 0, "[CX4H1]", "aliphatic sp3 CH"); R.add("CT2", 0, "[CX4H2]", "aliphatic sp3 CH2"); R.add("CT3", 0, "[CX4H3]", "aliphatic sp3 CH3")
    R.add("CT", 0, "[CX4H0]", "aliphatic sp3 C without H")
    if v36:
        R.add("HA1", 0, "[H][CX4H1]", "H on an aliphatic CH"); R.add("HA2", 0, "[H][CX4H2]", "H on an aliphatic CH2"); R.add("HA3", 0, "[H][CX4H3]", "H on an aliphatic CH3")
    else:
        R.add("HA", 0, "[H][CX4]", "aliphatic H")
    R.add("CA", 1, "c", "aromatic carbon"); R.add("HP", 0, "[H]c", "aromatic H")
    R.add("CE1", 1, "[CX3;H1,H0]=[CX3]", "alkene carbon RHC= (and R2C=)"); R.add("CE2", 1, "[CX3H2]=[CX3]", "alkene carbon H2C=")
    R.add("HE1", 1, "[H][CX3H1]=[CX3]", "H on an RHC= carbon"); R.add("HE2", 1, "[H][CX3H2]=[CX3]", "H on an H2C= carbon")
    R.add("C", 3, "[CX3](=[OX1])[NX3]", "amide (peptide) carbonyl carbon")
    R.add("CD", 3, "[CX3](=[OX1])[OX2]", "carbonyl carbon of an acid or ester")
    R.add("CC", 4, "[CX3](=[OX1])[OX1-]", "carboxylate carbon")
    R.add("CC", 4, "[CX3](=[OX1])[NX3H2]", "primary amide carbon (as the asparagine side chain)")
    R.add("O", 2, "[OX1]=[CX3][NX3]", "amide carbonyl oxygen")
    R.add("OB", 2, "[OX1]=[CX3][OX2]", "carbonyl oxygen of an acid or ester")
    R.add("OC", 3, "[OX1]~[CX3]~[OX1]", "carboxylate oxygen")
    R.add("OS", 2, ESTER_O, "ester oxygen")
    R.add("OH1", 1, "[OX2H1]", "hydroxyl oxygen")
    R.add("NH1", 2, "[NX3H1][CX3]=[OX1]", "amide nitrogen with one H"); R.add("NH2", 2, "[NX3H2][CX3]=[OX1]", "amide nitrogen with two H")
    R.add("N", 2, "[NX3H0][CX3]=[OX1]", "amide nitrogen without H (proline)"); R.add("NH3", 3, "[NX4+]", "ammonium nitrogen")
    R.add("H", 0, "[H][N,O]", "polar H"); R.add("HC", 1, "[H][NX4+]", "H on an ammonium nitrogen")
    R.add("S", 1, "[SX2]([#6])[#6]", "sulfide sulfur"); R.add("SM", 1, "[SX2][SX2]", "disulfide sulfur"); R.add("HS", 1, "[H][SX2]", "thiol H")
    water_ions(R, "OT", "HT")
    R.add("SOD", 0, "[Na]", "sodium ion"); R.add("CLA", 0, "[Cl-;X0]", "chloride ion"); R.add("POT", 0, "[K]", "potassium ion")
    return R


# ---------------------------------------------------------------- CHARMM36 lipids
def charmm_lipid():
    R = Rules()
    R.add("CTL1", 0, "[CX4H1]", "sp3 CH"); R.add("CTL2", 0, "[CX4H2]", "sp3 CH2"); R.add("CTL3", 0, "[CX4H3]", "sp3 CH3")
    R.add("CTL5", 2, "[CX4H3][NX4+]", "methyl on a quaternary ammonium (choline)")
    R.add("HAL1", 0, "[H][CX4H1]", "H on CH"); R.add("HAL2", 0, "[H][CX4H2]", "H on CH2"); R.add("HAL3", 0, "[H][CX4H3]", "H on CH3")
    R.add("HL", 2, "[H][CX4][NX4+]", "H on a carbon bonded to N+")
    R.add("CEL1", 1, "[CX3;H1,H0]=[CX3]", "alkene carbon =CH- (and =CR-)"); R.add("CEL2", 1, "[CX3H2]=[CX3]", "alkene carbon =CH2")
    R.add("HEL1", 1, "[H][CX3H1]=[CX3]", "H on =CH-"); R.add("HEL2", 1, "[H][CX3H2]=[CX3]", "H on =CH2")
    R.add("CL", 3, "[CX3](=[OX1])[OX2]", "ester / acid carbonyl carbon"); R.add("CCL", 4, "[CX3](=[OX1])[OX1-]", "carboxylate carbon")
    R.add("OBL", 2, "[OX1]=[CX3][OX2]", "ester / acid carbonyl oxygen"); R.add("OCL", 3, "[OX1]~[CX3]~[OX1]", "carboxylate oxygen")
    R.add("OSL", 2, ESTER_O, "ester oxygen"); R.add("OHL", 1, "[OX2H1]", "hydroxyl oxygen"); R.add("HOL", 1, "[H][OX2]", "hydroxyl H")
    R.add("NTL", 3, "[NX4+;H0]", "quaternary ammonium nitrogen"); R.add("NH3L", 3, "[NX4+;H3]", "ammonium nitrogen"); R.add("HCL", 3, "[H][NX4+]", "H on N+")
    R.add("PL", 0, "[PX4]", "phosphorus"); R.add("O2L", 1, "[OX1]~P", "phosphate non-bridging oxygen"); R.add("OSLP", 1, "[OX2](P)[#6]", "phosphate ester oxygen")
    water_ions(R, "OT", "HT")
    R.add("SOD", 0, "[Na]", "sodium ion"); R.add("CLA", 0, "[Cl-;X0]", "chloride ion"); R.add("POT", 0, "[K]", "potassium ion")
    return R


# ---------------------------------------------------------------- CHARMM36 carbohydrates (pyranoses, furanoses)
def charmm_carb():
    R = Rules()
    PYR = "[#6;r6;$(*@[OX2;r6])]"   # a pyranose ring carbon bonded to the ring oxygen
    R.add("CC331", 0, "[CX4H3]", "methyl carbon"); R.add("CC321", 0, "[CX4H2]", "exocyclic CH2 (C6)")
    R.add("CC3161", 1, "[CX4H1;r6]", "pyranose ring CH (C2-C4)")
    R.add("CC3263", 1, "[CX4H2;r6]", "pyranose ring CH2")
    R.add("CC3163", 2, "[CX4H1;r6;$(*@[OX2;r6]);!$(*[OX2;!r6])]", "pyranose C5 (ring oxygen and an exocyclic carbon)")
    R.add("CC3162", 3, "[CX4;r6;$(*@[OX2;r6]);$(*-[OX2;!R])]", "pyranose C1 (anomeric: ring oxygen and an exocyclic oxygen)")
    R.add("CC3151", 1, "[CX4H1;r5]", "furanose ring CH"); R.add("CC3251", 1, "[CX4H2;r5]", "furanose ring CH2")
    R.add("CC3152", 3, "[CX4;r5;$(*@[OX2;r5]);$(*-[OX2;!R])]", "furanose C1 (anomeric)")
    R.add("CC3153", 2, "[CX4H1;r5;$(*@[OX2;r5]);!$(*-[OX2;!R])]", "furanose C4 (ring oxygen and an exocyclic carbon)")
    R.add("OC3C61", 2, "[OX2;r6]", "pyranose ring oxygen"); R.add("OC3C51", 2, "[OX2;r5]", "furanose ring oxygen")
    R.add("OC311", 1, "[OX2H1]", "hydroxyl oxygen"); R.add("HCP1", 1, "[H][OX2]", "hydroxyl H")
    R.add("HCA1", 0, "[H][CX4H1]", "H on CH"); R.add("HCA2", 0, "[H][CX4H2]", "H on CH2"); R.add("HCA3", 0, "[H][CX4H3]", "H on CH3")
    R.add("CC2O1", 3, "[CX3](=[OX1])[NX3]", "amide carbon (N-acetyl)"); R.add("OC2D1", 2, "[OX1]=[CX3][NX3]", "amide oxygen"); R.add("NC2D1", 2, "[NX3][CX3]=[OX1]", "amide nitrogen")
    R.add("CC2O2", 3, "[CX3](=[OX1])[OX1-]", "carboxylate carbon (uronic acids)"); R.add("OC2D2", 3, "[OX1]~[CX3]~[OX1]", "carboxylate oxygen")
    water_ions(R, "OT", "HT")
    R.add("SOD", 0, "[Na]", "sodium ion"); R.add("CLA", 0, "[Cl-;X0]", "chloride ion"); R.add("POT", 0, "[K]", "potassium ion")
    return R


# ---------------------------------------------------------------- CHARMM36 nucleic acids
def charmm_nucl():
    """Types as CHARMM's nucleic-acid topology (top_all36_na.rtf) assigns them: sugar and phosphate by position, the bases
    atom by atom. Purine atoms lie within two bonds of a ring-fusion atom; pyrimidine atoms have none nearby."""
    R = Rules()
    PUR = "$(*@[R2]),$(*@*@[R2]),$([R2])"
    PYR = "r6;!$(*@[R2]);!$(*@*@[R2]);!R2"
    # phosphate
    R.add("P", 0, "[PX4]", "phosphodiester phosphorus"); R.add("P2", 1, "[PX4;$(P[OX2H1])]", "phosphorus of a terminal phosphate monoester")
    R.add("ON3", 1, "[OX1]~P", "non-bridging phosphate oxygen"); R.add("ON2", 1, "[OX2]([#6])P", "bridging (ester) phosphate oxygen")
    R.add("ON4", 2, "[OX2H1]P", "phosphate hydroxyl oxygen"); R.add("HN4", 2, "[H][OX2]P", "phosphate hydroxyl H")
    # sugar
    R.add("CN7", 0, "[CX4H1;r5]", "sugar ring CH (C3', C4')")
    R.add("CN7B", 2, "[CX4;r5;$(*[#7]);$(*@[OX2])]", "C1' (anomeric, bonded to the base)")
    R.add("CN7B", 2, "[CX4H1;r5;$(*[OX2H1]);$(*@[CX4]([#7]))]", "C2' of RNA (bonded to the 2'-OH)")
    R.add("CN8", 1, "[CX4H2;r5]", "C2' of DNA (CH2)")
    R.add("CN8B", 1, "[CX4H2;!R][CX4;r5]", "C5' (CH2 next to the ring)")
    R.add("ON6", 1, "[OX2;r5]", "sugar ring oxygen O4' (DNA)")
    R.add("ON6B", 2, "[OX2;r5;$(O1[CX4]([#7])[CX4]([OX2H1])[CX4][CX4]1)]", "sugar ring oxygen O4' (RNA)")
    R.add("ON5", 1, "[OX2H1][CX4]", "sugar hydroxyl oxygen (2'-OH, 5'/3' ends)"); R.add("HN5", 1, "[H][OX2][CX4]", "sugar hydroxyl H")
    R.add("HN7", 0, "[H][CX4H1]", "H on a sugar CH"); R.add("HN8", 0, "[H][CX4H2]", "H on a sugar CH2")
    # bases
    R.add("CN9", 2, "[CX4H3][#6;r6]", "thymine methyl carbon"); R.add("HN9", 2, "[H][CX4H3]", "thymine methyl H")
    R.add("CN1", 2, "[#6;R](=[OX1])(~[#7H1])~[#6]", "carbonyl carbon next to N-H and C (T/U C4, G C6)")
    R.add("CN1T", 3, f"[#6;{PYR}](=[OX1])(~[#7][CX4])~[#7H1]", "C2 of thymine / uracil")
    R.add("CN1", 3, f"[#6;{PYR}](=[OX1])(~[#7][CX4])~[#7X2]", "C2 of cytosine (CN1 in CHARMM36)")
    R.add("ON1", 1, "[OX1]=[#6;R]", "base carbonyl oxygen"); R.add("ON1C", 3, f"[OX1]=[#6;{PYR}](~[#7][CX4])~[#7X2]", "O2 of cytosine")
    R.add("CN2", 2, "[#6;R]~[NX3H2]", "ring carbon bearing an amino group (C C4, A C6, G C2)")
    R.add("NN1", 2, "[NX3H2][#6;R]", "amino nitrogen of a base"); R.add("HN1", 2, "[H][NX3H2]", "amino H of a base")
    R.add("CN3", 1, f"[#6;{PYR}]", "pyrimidine ring carbon (C5, C6)"); R.add("CN3T", 2, f"[#6;{PYR}][CX4H3]", "C5 of thymine")
    R.add("CN4", 1, f"[#6H1;!R2;{PUR}]", "purine CH (C8, adenine C2)")
    R.add("CN5", 1, "[#6;R2]", "purine ring-fusion carbon (C4, C5)"); R.add("CN5G", 2, "[#6;R2]~[#6]=[OX1]", "C5 of guanine")
    R.add("NN4", 1, "[#7X2;r5]", "purine N7")
    R.add("NN2", 1, "[#7X3;r5]([CX4])", "adenine N9 (glycosidic)"); R.add("NN2B", 2, "[#7X3;r5;$(*[CX4]);$(*~[#6]~[#6]~[#6]=[OX1])]", "guanine N9 (glycosidic)")
    R.add("NN2B", 1, f"[#7X3;{PYR}]([CX4])", "thymine / uracil N1 (glycosidic)")
    R.add("NN2", 3, f"[#7X3;{PYR};$(*([CX4])~[#6](=[OX1])~[#7X2])]", "cytosine N1 (glycosidic)")
    R.add("NN2U", 2, f"[#7H1;{PYR}]", "thymine / uracil N3 (N-H)"); R.add("NN2G", 2, f"[#7H1;r6;{PUR}]", "guanine N1 (N-H)")
    R.add("NN3", 1, f"[#7X2;{PYR}]", "cytosine N3"); R.add("NN3A", 1, f"[#7X2;r6;{PUR}]", "adenine N1 / N3")
    R.add("NN3G", 2, f"[#7X2;r6;{PUR};$(*~[#6](~[NX3H2])~[#7H1])]", "guanine N3")
    R.add("HN2", 2, "[H][#7;R]", "H on a base ring nitrogen"); R.add("HN3", 1, "[H][#6;R;X3]", "H on a base ring carbon")
    water_ions(R, "OT", "HT")
    R.add("SOD", 0, "[Na]", "sodium ion"); R.add("CLA", 0, "[Cl-;X0]", "chloride ion"); R.add("POT", 0, "[K]", "potassium ion")
    return R


# ---------------------------------------------------------------- united atom: CHARMM19, GROMOS 54A7, TraPPE-UA
def charmm19():
    R = Rules()
    R.add("CH1E", 0, "[CX4H1]", "united CH"); R.add("CH2E", 0, "[CX4H2]", "united CH2"); R.add("CH3E", 0, "[CX4H3]", "united CH3")
    R.add("CT", 0, "[CX4H0]", "aliphatic carbon without H"); R.add("CR1E", 1, "[cH1]", "united aromatic CH")
    R.add("C", 1, "[c,C;H0;X3]", "carbonyl or aromatic carbon without H")
    R.add("O", 1, "[OX1]=[#6]", "carbonyl oxygen"); R.add("OC", 3, "[OX1]~[CX3]~[OX1]", "carboxylate oxygen")
    R.add("OS", 0, "[OX2]([#6])[#6]", "ester (and ether) oxygen"); R.add("OH1", 1, "[OX2H1]", "hydroxyl oxygen")
    R.add("NH1", 1, "[NX3H1]", "nitrogen with one H"); R.add("NH2", 1, "[NX3H2]", "nitrogen with two H"); R.add("N", 0, "[NX3H0]", "nitrogen without H")
    R.add("NH3", 3, "[NX4+]", "ammonium nitrogen"); R.add("NR", 1, "[nX2]", "aromatic nitrogen without H"); R.add("NP", 1, "[nX3H1]", "pyrrole nitrogen")
    R.add("H", 0, "[H][N,O,n]", "polar H"); R.add("HC", 2, "[H][NX4+]", "H on a charged nitrogen")
    R.add("S", 0, "[SX2]", "sulfur")
    R.add("OT", 8, "[OX2H2]", "TIP3P water O"); R.add("H", 8, "[H][OX2H2]", "water H")
    return R


def gromos():
    R = Rules()
    R.add("CH1", 0, "[CX4H1]", "united CH"); R.add("CH2", 0, "[CX4H2]", "united CH2"); R.add("CH3", 0, "[CX4H3]", "united CH3")
    R.add("CH4", 0, "[CX4H4]", "united CH4"); R.add("CH0", 0, "[CX4H0]", "aliphatic carbon without H")
    R.add("CH3p", 2, "[CX4H3][NX4+]", "united CH3 on a charged nitrogen")
    R.add("CR1", 1, "[#6X3H1]", "united aromatic or alkene CH")
    R.add("C", 1, "[#6X3H0]", "bare sp2 carbon (carbonyl, aromatic without H)")
    R.add("O", 1, "[OX1]=[#6]", "carbonyl oxygen"); R.add("OM", 3, "[OX1]~[CX3,PX4]~[OX1]", "carboxylate / phosphate oxygen")
    R.add("OA", 0, "[OX2]", "hydroxyl or ether oxygen"); R.add("OE", 2, ESTER_O, "ester oxygen")
    R.add("N", 1, "[NX3][CX3]=[OX1]", "amide nitrogen"); R.add("NT", 0, "[NX3;!$(N[CX3]=[OX1])]", "amine nitrogen")
    R.add("NL", 3, "[NX4+]", "ammonium nitrogen")
    R.add("H", 0, "[H][N,O,S]", "polar H"); R.add("S", 0, "[SX2]", "sulfur"); R.add("P", 0, "[PX4]", "phosphorus")
    R.add("OWS", 8, "[OX2H2]", "SPC/E water O"); R.add("HWS", 8, "[H][OX2H2]", "SPC/E water H")
    return R


def trappe_ua(full=True):
    R = Rules()
    R.add("CH4", 0, "[CX4H4]", "methane"); R.add("CH3", 0, "[CX4H3]", "CH3"); R.add("CH2", 0, "[CX4H2]", "CH2")
    R.add("CH1", 0, "[CX4H1]", "CH"); R.add("CH0", 0, "[CX4H0]", "quaternary C")
    if not full:
        return R
    R.add("CYC6", 1, "[CX4H2;r6]", "CH2 in a six-membered ring"); R.add("CYC5", 1, "[CX4H2;r5]", "CH2 in a five-membered ring")
    R.add("CE2", 1, "[CX3H2]=[CX3]", "alkene CH2="); R.add("CE1", 1, "[CX3H1]=[CX3]", "alkene CH="); R.add("CE0", 1, "[CX3H0]=[CX3]", "alkene C=")
    R.add("CB", 1, "[cH1]", "benzene CH"); R.add("CB0", 1, "[cH0]", "substituted aromatic C")
    R.add("OA", 1, "[OX2H1]", "alcohol oxygen"); R.add("HA", 1, "[H][OX2]", "alcohol H"); R.add("OE", 0, "[OX2]([#6])[#6]", "ether oxygen")
    R.add("CK", 2, "[CX3](=[OX1])([#6])[#6]", "ketone carbon"); R.add("OK", 2, "[OX1]=[CX3]([#6])[#6]", "ketone oxygen")
    R.add("STI", 1, "[SX2H1]", "thiol sulfur"); R.add("HTI", 1, "[H][SX2]", "thiol H")
    R.add("SS", 1, "[SX2][SX2]", "disulfide sulfur"); R.add("S", 0, "[SX2]([#6])[#6]", "sulfide sulfur")
    R.add("CNI", 2, "[CX2]#[NX1]", "nitrile carbon"); R.add("NI", 2, "[NX1]#[CX2]", "nitrile nitrogen")
    return R


# ---------------------------------------------------------------- TraPPE-EH (explicit-hydrogen aromatics)
def trappe_eh():
    R = Rules()
    BENZ = "c;r6;!$(c~[!#6;!#1]);!$(c:[!#6])"
    R.add("C_1", 0, f"[{BENZ}]", "benzene carbon"); R.add("H_1", 0, "[H]c", "benzene H")
    R.add("N1_2", 2, "[nX2;r6;!$(n:c:n);!$(n:n);!$(n:c:c:n)]", "pyridine nitrogen")
    R.add("C2_2", 2, "[c;r6;$(c:[nX2]);!$(c:[n]:c:[n])]", "pyridine C2 (next to N)")
    R.add("C3_2", 2, "[c;r6;$(c:c:[nX2]);!$(c:[n])]", "pyridine C3")
    R.add("C4_2", 2, "[c;r6;$(c:c:c:[nX2]);!$(c:[n]);!$(c:c:[n])]", "pyridine C4")
    R.add("H_2", 2, "[H]c(:c):[nX2,$(c:c:n),$(c:c:c:n)]", "pyridine H")
    R.add("O_8", 3, "[OX2H1]c", "phenol oxygen"); R.add("HO_8", 3, "[H][OX2]c", "phenol hydroxyl H"); R.add("C1_8", 3, "c[OX2H1]", "phenol C1")
    R.add("N1_9", 3, "[nX3H1;r5]", "pyrrole nitrogen"); R.add("H6_9", 3, "[H][nX3;r5]", "pyrrole N-H")
    R.add("C2_9", 3, "[c;r5;$(c:[nH1])]", "pyrrole C2"); R.add("C3_9", 3, "[c;r5;$(c:c:[nH1]);!$(c:[nH1])]", "pyrrole C3")
    R.add("H_9", 3, "[H][c;r5;$(c:[nH1]),$(c:c:[nH1])]", "pyrrole C-H")
    R.add("O1_7", 3, "[o;r5]", "furan oxygen"); R.add("C2_7", 3, "[c;r5;$(c:o)]", "furan C2"); R.add("C3_7", 3, "[c;r5;$(c:c:o);!$(c:o)]", "furan C3")
    R.add("H_7", 3, "[H][c;r5;$(c:o),$(c:c:o)]", "furan H")
    R.add("S1_6", 3, "[s;r5]", "thiophene sulfur"); R.add("C2_6", 3, "[c;r5;$(c:s)]", "thiophene C2"); R.add("C3_6", 3, "[c;r5;$(c:c:s);!$(c:s)]", "thiophene C3")
    R.add("H_6", 3, "[H][c;r5;$(c:s),$(c:c:s)]", "thiophene H")
    return R


# Terms a file lacks, by analogy (as parmchk2 and the CGenFF program fill missing terms): the same element and
# hybridisation, nearest in hydrogen count or neighbourhood. Every term found this way is listed as estimated.
ALKYL = lambda a, b, c: {a: [b, c], b: [a, c], c: [b, a]}
ANALOGIES = {
    "amber": {"G_Ck": ["CA"], "G_Ha": ["HA"], "H1": ["HC"], "H2": ["H1", "HC"], "HP": ["HC"]},
    "charmm": {**ALKYL("C1", "C2", "C3"), "CAE1": ["C1"], "CAE2": ["C2"], "CAE3": ["C3"], "CTHF": ["CAE2", "C2"], "CTHP": ["CAE2", "C2"], "CYCP": ["C2"],
               "HAE1": ["H1"], "HAE2": ["H1"], "HAE3": ["H1"], "HTHF": ["H1"], "C15": ["C16"], "C16": ["C15"], "H7": ["H8"], "H8": ["H7"],
               "N1": ["N2"], "N2": ["N1"], "N5": ["N2"], "OET": ["OE"], "OE": ["OET"]},
    "charmm22-prot": {**ALKYL("CT1", "CT2", "CT3"), "CT": ["CT1"], "CE1": ["CE2"], "CE2": ["CE1"], "HE1": ["HE2"], "HE2": ["HE1"]},
    "charmm36-prot": {**ALKYL("CT1", "CT2", "CT3"), "CT": ["CT1"], "HA1": ["HA2", "HA3"], "HA2": ["HA3", "HA1"], "HA3": ["HA2", "HA1"],
                      "CE1": ["CE2"], "CE2": ["CE1"], "HE1": ["HE2"], "HE2": ["HE1"]},
    "charmm36-lipid": {**ALKYL("CTL1", "CTL2", "CTL3"), "HAL1": ["HAL2", "HAL3"], "HAL2": ["HAL3", "HAL1"], "HAL3": ["HAL2", "HAL1"]},
    "charmm36-nucl": {"P2": ["P"]},
    "charmm36-carb": {"CC331": ["CC321"], "CC321": ["CC331"], "HCA1": ["HCA2"], "HCA2": ["HCA1"], "HCA3": ["HCA2"]},
    "charmm19": {**ALKYL("CH1E", "CH2E", "CH3E"), "CT": ["CH1E"]},
    "gromos-54a7": {**ALKYL("CH1", "CH2", "CH3"), "CH0": ["CH1"], "CH3p": ["CH3"], "CR1": ["C"], "OE": ["OA"]},
    "trappe-ua": {**ALKYL("CH1", "CH2", "CH3"), "CH0": ["CH1"], "CYC6": ["CH2"], "CYC5": ["CH2"], "CE1": ["CE2", "CE0"], "CE2": ["CE1"], "CE0": ["CE1"]},
}
# Entries the converted files leave out although the force field defines them
EXTENSIONS = {
    "gromos-54a7": {"pairs": [{"name": "H", "match": ["H"], "params": [0, 0],
                               "comment": "CAPS extension: GROMOS polar hydrogens carry no Lennard-Jones (C6 = C12 = 0; Schmid et al. 2011)"}]},
}

SETS = {
    "amber": (amber(), "AMBER (Cornell et al. 1995 family)", False),
    "charmm": (charmm(), "CHARMM (library set)", False),
    "charmm22-prot": (charmm_prot(False), "CHARMM22 proteins", False),
    "charmm36-prot": (charmm_prot(True), "CHARMM36 proteins", False),
    "charmm36-lipid": (charmm_lipid(), "CHARMM36 lipids", False),
    "charmm36-carb": (charmm_carb(), "CHARMM36 carbohydrates", False),
    "charmm36-nucl": (charmm_nucl(), "CHARMM36 nucleic acids", False),
    "charmm19": (charmm19(), "CHARMM19 united atom", True),
    "gromos-54a7": (gromos(), "GROMOS 54A7 (united atom)", True),
    "trappe-ua": (trappe_ua(), "TraPPE-UA", True),
    "trappe1998-moltemplate": (trappe_ua(False), "TraPPE-UA alkanes", True),
    "trappe-eh": (trappe_eh(), "TraPPE-EH (explicit-hydrogen aromatics)", False),
}

if __name__ == "__main__":
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for target, (rules, label, ua) in SETS.items():
        ff_p = os.path.join(DATA, "forcefields", target + ".json")
        ff = json.load(open(ff_p))
        have = {t["name"] for t in ff["atom_types"]}
        absent = sorted({r["type"] for r in rules if r["type"] not in have})
        doc = {"format": "caps-typing", "version": 1, "forcefield": label,
               "description": f"CAPS rules for {label} from the force field's own type definitions (bench/typing/make_biomolecular_rules.py); "
                              "among matching rules the highest priority wins. " + ("United-atom: hydrogens on carbon fold into their carbons before typing. " if ua else "")
                              + "Types this file lacks leave their atoms untyped: " + (", ".join(absent) or "none") + ".",
               "unknown_types": "untyped", "rules": list(rules)}
        if ua:
            doc["united_atom"] = True
        if target in ANALOGIES:
            doc["analogies"] = ANALOGIES[target]
            doc["analogy_source"] = f"{label} analogues: the same element and hybridisation, nearest in hydrogen count"
        for kind, entries in EXTENSIONS.get(target, {}).items():
            have_e = {e["name"] for e in ff.get(kind, [])}
            ff.setdefault(kind, []).extend(e for e in entries if e["name"] not in have_e)
        out = os.path.join(DATA, "typing", target + ".typing.json")
        json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
        ff["typing"] = f"../typing/{target}.typing.json"
        json.dump(ff, open(ff_p, "w"), ensure_ascii=False, indent=1)
        for e in cat["forcefields"]:
            if e["id"] == target:
                e["typing"] = {"rules": f"typing/{target}.typing.json", "evidence": f"CAPS rules from {label}'s type definitions (bench/typing/make_biomolecular_rules.py)"
                               + ("; united atom" if ua else "")}
        print(f"{out}: {len(rules)} rules; absent from {target}: {', '.join(absent) or 'none'}")
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
