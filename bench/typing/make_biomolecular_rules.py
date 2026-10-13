#!/usr/bin/env python3
"""CAPS typing rules for the CHARMM (generic set) and TraPPE files of the library, written from each force field's own
type definitions (the meaning of every type name in its parameter and topology files):

  charmm           the library's generic CHARMM set (CT1-3, CA, CE1/CE2, CD/CC amide / acid, OB, OS, ether CAE/OET ...)
  trappe-ua        Martin & Siepmann (1998) and later TraPPE-UA: CH4 / CH3 / CH2 / CH / C, ring CH2, alkenes, benzene CH / C,
                   alcohol, ether, ketone, thiol / sulfide, nitrile
  trappe1998-moltemplate   the same united-atom alkane rules (CH2, CH3, CH4)
  trappe-eh        Rai & Siepmann (2007), explicit-hydrogen aromatics: benzene, phenol, pyridine, pyrrole, furan, thiophene

The rules for the biomolecular sets (AMBER, GLYCAM, the CHARMM protein / lipid / carbohydrate / nucleic-acid sets,
CHARMM19, GROMOS 54A7) were removed on 2026-10-14: CAPS focuses on polymers, crystals and materials.

United-atom files set "united_atom": the structure's hydrogens on carbon fold into their carbons before typing (CAPS's
prepare_for_forcefield); the rules then see CH2 and friends as carbons with implicit hydrogens. Where a file has no
type for an environment its atoms are left
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
    def add(self, t, prio, smarts, desc, **cond):
        """cond: in_molecule / not_in_molecule (SMARTS on the atom's own molecule)"""
        self.append({"type": t, "smarts": smarts, "priority": prio, "description": desc, **cond})


def sugar_ring_smarts():
    """an sp3 C of a carbohydrate ring: a five- or six-membered ring of carbons and one O, the C carrying an exocyclic O
    or N or next (in the ring) to a C that does (tetrahydrofuran and tetrahydropyran do not count)"""
    rings = []
    for n in (6, 5):
        seq = ["C"] * (n - 1) + ["O"]
        for k in range(n - 1):   # start at each carbon
            r = seq[k:] + seq[:k]
            rings.append("$(" + r[0].replace("C", "[CX4;R]") + "1" + "".join("@" + ("[CX4;R]" if x == "C" else "[OX2;R]") for x in r[1:]) + "1)")
    return "[CX4;" + ",".join(rings) + ";$(*!@[#8,#7]),$(*@[CX4]!@[#8,#7])]"


# ---------------------------------------------------------------- CHARMM (library's generic keys)
def charmm():
    """the library's generic CHARMM set, whose keys depend on the kind of molecule: amino acids and peptides, lipids
    (quaternary ammonium, phosphate and sulfate esters, sulfoxides, sterols: the _L keys), carbohydrates (CS / CF ...),
    ethers and saturated hydrocarbons (CAE / HAE / OET), and the rest; molecule classes by in_molecule"""
    R = Rules()
    STEROID = "[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]~[#6]~[#6]1~[#6]~2~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~2"
    LIPID = f"[$([NX4]([CH3])([CH3])[CH3]),$([PX4]),$([SX4]),$([SX3]=[OX1]),$({STEROID})]"
    SUGAR = sugar_ring_smarts()
    # ethers and saturated hydrocarbons: only C, H and ether O, nothing unsaturated
    ETHERS = {"not_in_molecule": ["[!#6;!#1;!#8]", "[OX2H1]", "[OX1]", "[#6;X3,X2]", "a", SUGAR]}
    LIP = {"in_molecule": [LIPID]}
    SUG = {"in_molecule": [SUGAR], "not_in_molecule": [LIPID]}
    OTHER = {"not_in_molecule": [LIPID, SUGAR]}
    # sp3 carbons and their H
    R.add("C1", 0, "[CX4;H1,H0]", "sp3 CH / C (CT1)"); R.add("C2", 0, "[CX4H2]", "sp3 CH2 (CT2)"); R.add("C3", 0, "[CX4H3]", "sp3 CH3 (CT3)")
    R.add("H1", 0, "[H][CX4]", "H on sp3 C")
    R.add("CAE1", 2, "[CX4;H1,H0]", "ether / alkane CH, C", **ETHERS); R.add("CAE2", 2, "[CX4H2]", "ether / alkane CH2", **ETHERS)
    R.add("CAE3", 2, "[CX4H3]", "ether / alkane CH3", **ETHERS)
    R.add("HAE1", 2, "[H][CX4;H1,H0]", "H on ether / alkane CH", **ETHERS); R.add("HAE2", 2, "[H][CX4H2]", "H on ether / alkane CH2", **ETHERS)
    R.add("HAE3", 2, "[H][CX4H3]", "H on ether / alkane CH3", **ETHERS)
    R.add("CTHF", 3, "[CX4;$(*1@*@*@*@[OX2]1),$(*1@*@*@[OX2]@*1)]", "tetrahydrofuran ring C", **ETHERS)
    R.add("HTHF", 3, "[H][CX4;$(*1@*@*@*@[OX2]1),$(*1@*@*@[OX2]@*1)]", "H on a tetrahydrofuran ring C", **ETHERS)
    R.add("CTHP", 3, "[CX4;$(*1@*@*@*@*@[OX2]1),$(*1@*@*@*@[OX2]@*1),$(*1@*@*@[OX2]@*@*1)]", "tetrahydropyran ring C", **ETHERS)
    R.add("CYCP", 3, "[CX4;$(*1@[CX4]@[CX4]@[CX4]@[CX4]1)]", "cyclopentane C", **ETHERS)
    R.add("HTHF", 3, "[H][CX4;$(*1@[CX4]@[CX4]@[CX4]@[CX4]1)]", "H on a cyclopentane C", **ETHERS)
    R.add("OET", 1, "[OX2]([CX4])[CX4]", "ether O")
    # lipids
    R.add("C_L1", 4, "[CX4;H1,H0]", "lipid CH / C", **LIP); R.add("C_L2", 4, "[CX4H2]", "lipid CH2", **LIP); R.add("C_L3", 4, "[CX4H3]", "lipid CH3", **LIP)
    R.add("C_L4", 5, "[CX4H3][NX4]", "methyl on a quaternary ammonium N", **LIP)
    R.add("H_L1", 4, "[H][CX4;H1,H0]", "H on lipid CH", **LIP); R.add("H_L2", 4, "[H][CX4H2]", "H on lipid CH2", **LIP)
    R.add("H_L3", 4, "[H][CX4H3]", "H on lipid CH3", **LIP)
    R.add("HL", 5, "[H][CX4][NX4H0]", "H on C of a quaternary ammonium", **LIP)
    R.add("NT", 5, "[NX4H0]", "quaternary ammonium N")
    R.add("N3L", 5, "[NX4;!H0]", "lipid ammonium N", **LIP)
    R.add("OSL", 5, "[OX2;$(O[CX4H2]);$(O[$([CX3]=[OX1]),$([SX4])])]", "lipid ester O (on a CH2)", **LIP)
    R.add("OE", 4, "[OX2]([#6])[$([CX3]=[OX1]),$([SX4]),$([PX4])]", "ester O (carboxylate, phosphate, sulfate)")
    R.add("OP", 5, "[OX1]~[P,SX4]", "phosphate / sulfate O"); R.add("P1", 5, "[PX4]", "phosphate P"); R.add("S5", 5, "[SX4]", "sulfate S")
    R.add("S4", 5, "[SX3]=[OX1]", "sulfoxide S"); R.add("OC3", 5, "[OX1]=[SX3]", "sulfoxide O")
    # carbohydrates
    PYR = "$(*1@*@*@*@*@[OX2]1)"; FUR = "$(*1@*@*@*@[OX2]1)"
    ADJ = "$(*@[OX2;R])"; EXOC = "$(*!@[#6])"; EXOO = "$(*!@[#8])"
    R.add("CS", 6, f"[CX4;R;{PYR},$(*1@*@*@*@[OX2]@*1),$(*1@*@*@[OX2]@*@*1)]", "pyranose ring C", **SUG)
    R.add("CS1", 7, f"[CX4;R;{ADJ};!{EXOC};$(*1@*@*@*@*@[OX2]1)]", "pyranose C1", **SUG)
    R.add("CS5", 7, f"[CX4;R;{ADJ};{EXOC};$(*1@*@*@*@*@[OX2]1)]", "pyranose C5", **SUG)
    R.add("CF", 6, f"[CX4;R;{FUR},$(*1@*@*@[OX2]@*1)]", "furanose ring C", **SUG)
    R.add("CFD", 7, f"[CX4H2;R;$(*1@*@*@[OX2]@*1)]", "furanose deoxy C2", **SUG)
    R.add("CF1", 7, f"[CX4;R;{ADJ};!{EXOC};{FUR}]", "furanose C1", **SUG)
    R.add("CF4", 7, f"[CX4;R;{ADJ};{EXOC};!{EXOO};{FUR}]", "furanose C4", **SUG)
    R.add("CFF2", 7, f"[CX4;R;{ADJ};{EXOC};{EXOO};{FUR}]", "fructofuranose C2", **SUG)
    R.add("CAS2", 6, "[CX4;!R;$(*[CX4;R])]", "exocyclic C of a sugar ring (C6)", **SUG)
    R.add("CAS3", 6, "[CX4H3;$(*[OX2][CX4;R])]", "O-methyl of a glycoside", **SUG)
    R.add("HAS1", 6, "[H][CX4H1]", "H on a sugar CH", **SUG); R.add("HAS2", 6, "[H][CX4H2]", "H on a sugar CH2", **SUG)
    R.add("HAS3", 6, "[H][CX4H3]", "H on a sugar CH3", **SUG)
    R.add("OHS", 6, "[OX2H1]", "sugar hydroxyl O", **SUG); R.add("HOS", 6, "[H][OX2H1]", "sugar hydroxyl H", **SUG)
    R.add("OS", 7, "[OX2;R;$(*1@*@*@*@*@*1)]([CX4])[CX4]", "pyranose ring O", **SUG)
    R.add("OF", 7, "[OX2;R;$(*1@*@*@*@*1)]([CX4])[CX4]", "furanose ring O", **SUG)
    # amino acids and the rest
    AAH = "[CX4]([NX3,NX4])[CX3]=[OX1]"
    R.add("H3", 5, f"[H]{AAH}", "backbone H alpha", **OTHER)
    R.add("C8", 6, "[CX4;r5]([NX4,NX3])[CX3]=[OX1]", "proline C alpha")
    R.add("C9", 6, "[CX4H2;r5;$(*1@*@*@[NX3,NX4]@[CX4]1([CX3]=[OX1])),$(*1@*@[NX3,NX4]@[CX4]([CX3]=[OX1])@*1)]", "proline C beta / gamma")
    R.add("C10", 6, "[CX4H2;r5;$(*1@[NX3,NX4]@[CX4]([CX3]=[OX1])@*@*1)]", "proline C delta")
    R.add("N4", 6, "[NX4;r5;$(*@[CX4][CX3]=[OX1])]", "proline N")
    R.add("C17", 5, "[CX4][SX1]", "thiolate C"); R.add("S3", 5, "[SX1]", "thiolate S")
    R.add("S1", 3, "[SX2;!$(*[SX2])]", "thiol / thioether S"); R.add("S2", 4, "[SX2]([#6])[SX2]", "disulfide S")
    R.add("HS", 3, "[H][SX2]", "thiol H")
    R.add("C4", 1, "c", "aromatic C"); R.add("H2", 1, "[H]c", "aromatic H")
    R.add("C15", 2, "[CX3H2]=[CX3]", "vinyl CH2="); R.add("C16", 2, "[CX3;H1,H0;!a]=[CX3;!a]", "alkene C")
    R.add("H7", 2, "[H][CX3H2]=[CX3]", "H on a vinyl CH2="); R.add("H8", 2, "[H][CX3;H1;!a]=[CX3;!a]", "H on an alkene CH")
    R.add("C5", 3, "[CX3](=[OX1])[OX2,#1,#6]", "acid / ester / aldehyde C=O")
    R.add("C6", 4, "[CX3](~[OX1])~[OX1]", "carboxylate C"); R.add("OC3", 4, "[OX1]~[CX3]~[OX1]", "carboxylate O")
    R.add("C6", 4, "[CX3](=[OX1])[NX3H2]", "primary amide C"); R.add("C7", 4, "[CX3](=[OX1])[NX3;!H2]", "peptide / secondary amide C")
    R.add("OC1", 3, "[OX1]=[CX3][OX2]", "acid / ester carbonyl O"); R.add("OC2", 3, "[OX1]=[CX3][#1,#6,#7]", "amide / aldehyde O")
    R.add("H4", 4, "[H][CX3]=[OX1]", "aldehyde H")
    R.add("OH", 2, "[OX2H1]", "hydroxyl O"); R.add("HO", 2, "[H][OX2]", "hydroxyl H")
    R.add("OC3", 5, "[OX1]~c", "phenoxide O")
    R.add("N3", 3, "[NX4]", "ammonium N"); R.add("HN", 2, "[H][N,n]", "H on N")
    R.add("N1", 3, "[NX3H2]", "primary amide / amine N")
    R.add("N2", 3, "[NX3H1,NX2;$(*[CX3]=[OX1]),$(*[CX4][CX3]=[OX1])]", "peptide N")
    R.add("HO", 4, "[H][NX3][CX3]=[OX1]", "amide H (the library's hydroxyl H)")
    # arginine, tryptophan, histidine
    R.add("C7", 5, "[CX3]([NX3])([NX3])=[NX3]", "guanidinium C"); R.add("N7", 5, "[NX3;$(*~[CX3](~[NX3])~[NX3])]", "guanidinium N")
    R.add("N6", 5, "[nX3;$(n1ccc2ccccc12)]", "indole N"); R.add("HO", 5, "[H][nX3;$(n1ccc2ccccc12)]", "indole H(N)")
    R.add("C11", 5, "[c;$(c1cnc2ccccc12);$(*[CX4])]", "indole C3"); R.add("C12", 5, "[c;R2;$(*1:*:*:[nX3]:*:1),$(*1:*:[nX3]:*:*:1),$(*1:[nX3]:*:*:*:1)]", "indole ring fusion C")
    IMI = "$(*1:*:*:*:*1)"
    R.add("C14", 5, "[c;r5;$(c(:n):n)]", "imidazole C2"); R.add("C13", 5, "[c;r5;$(c1:n:c:n:c1)]", "imidazole C4 / C5")
    R.add("N10", 6, "[n;r5;$(n1:c:[nH]:c:c1),$(n1:c:n:c:c1);$([nH]1:c:[nH]:c:c1),$([nH]1:c:[nH+]:c:c1),$([nH+]1:c:[nH]:c:c1)]", "imidazolium N")
    R.add("N8", 5, "[nX3H1;r5;$(n1:c:n:c:c1)]", "imidazole N-H"); R.add("N9", 5, "[nX2;r5;$(n1:c:n:c:c1)]", "imidazole N")
    R.add("HO", 5, "[H][nX3;r5;$(n1:c:n:c:c1)]", "imidazole H(N)")
    R.add("H4", 5, "[H][c;r5;$(c(:n):n)]", "imidazole H2")
    R.add("H6", 5, "[H][c;r5;$(c1:n:c:n:c1);!$(c(:n):n)]", "imidazole H4")
    R.add("H5", 6, "[H][c;r5;$(c(:[nH]):[nH])]", "imidazolium H2")
    R.add("H4", 6, "[H][c;r5;$(c1:[nH]:c:[nH]:c1);!$(c(:n):n)]", "imidazolium H4")
    # heme
    R.add("FE", 8, "[Fe]", "heme iron"); R.add("NPH", 8, "[n,N;$(*~[Fe])]", "heme pyrrole N")
    R.add("CPA", 8, "[#6;$(*~[n,N]~[Fe])]", "heme alpha C"); R.add("CPB", 8, "[#6;r5;$(*~[#6]~[n,N]~[Fe]);!$(*~[n,N]~[Fe])]", "heme beta C")
    R.add("CPM", 8, "[#6;!r5;$(*(~[#6]~[n,N])~[#6]~[n,N]~[Fe]);!$(*~[n,N])]", "heme meso C")
    R.add("H1", 8, "[H][#6;!r5;$(*(~[#6]~[n,N])~[#6]~[n,N]~[Fe])]", "heme meso H")
    # small molecules
    R.add("C19", 8, "[#6X1]~[#8X1]", "carbon monoxide C"); R.add("OM", 8, "[#8X1]~[#6X1]", "carbon monoxide O")
    R.add("CST", 8, "[CX2](=[OX1])=[OX1]", "carbon dioxide C"); R.add("OST", 8, "[OX1]=[CX2]=[OX1]", "carbon dioxide O")
    R.add("C18", 8, "[CX3](~[OX1])(~[OX1])~[OX1]", "carbonate C"); R.add("OCA", 8, "[OX1]~[CX3](~[OX1])~[OX1]", "carbonate O")
    R.add("CLAL", 6, "[Cl][CX4]", "alkyl chloride Cl")
    R.add("OT", 8, "[OX2H2]", "water O (TIP3P, CHARMM)"); R.add("HW", 8, "[H][OX2H2]", "water H (TIP3P, CHARMM)")
    for e, t in (("Na", "Na+"), ("Cl", "Cl-"), ("Cu", "Cu2+"), ("Zn", "Zn2+")):
        R.add(t, 8, f"[{e};X0]", f"{t} ion")
    return R


# ---------------------------------------------------------------- united atom: TraPPE-UA
def trappe_ua(full=True):
    R = Rules()
    R.add("CH4", 0, "[CX4H4]", "methane"); R.add("CH3", 0, "[CX4H3]", "CH3"); R.add("CH2", 0, "[CX4H2]", "CH2")
    R.add("CH1", 0, "[CX4H1]", "CH"); R.add("CH0", 0, "[CX4H0]", "quaternary C")
    if not full:
        return R
    R.add("CYC6", 1, "[CX4H2;r6]", "CH2 in a six-membered ring"); R.add("CYC5", 1, "[CX4H2;r5]", "CH2 in a five-membered ring")
    R.add("CYC", 1, "[CX4H2;R;!r3;!r4;!r5;!r6]", "CH2 in a ring of seven or more")
    # conjugated dienes (isoprene): the ends CEE2 / CEE1, the linked middle carbons CEL by H count, _1 on the more
    # substituted double bond, _2 on the other (butadiene's two middles are alike)
    DIENE = "[CX3]=[CX3]-[CX3]=[CX3]"
    R.add("CEE2", 2, "[CX3H2;$(*=[CX3]-[CX3]=[CX3])]", "terminal =CH2 of a conjugated diene")
    R.add("CEE1", 2, "[CX3H1;$(*=[CX3]-[CX3]=[CX3]);!$(*-[CX3]=[CX3])]", "terminal =CH- of a conjugated diene")
    R.add("CEL0_1", 3, "[CX3H0;$(*(=[CX3])-[CX3]=[CX3])]", "middle =C< of a conjugated diene")
    R.add("CEL1_2", 3, "[CX3H1;$(*(=[CX3])-[CX3H0]=[CX3])]", "middle =CH- of a diene whose other middle carbon is substituted")
    R.add("CEL1_1", 3, "[CX3H1;$(*(=[CX3])-[CX3H1]=[CX3])]", "middle =CH- of a diene (butadiene)")
    # O2 (TraPPE): two O and a massless centre site bonded to both
    R.add("O", 3, "[O;$(O(~[#0])~O)]", "O2 oxygen"); R.add("E1", 3, "[#0;$(*(~O)~O)]", "O2 centre site")
    R.add("CE2", 1, "[CX3H2]=[CX3]", "alkene CH2="); R.add("CE1", 1, "[CX3H1]=[CX3]", "alkene CH="); R.add("CE0", 1, "[CX3H0]=[CX3]", "alkene C=")
    R.add("CB", 1, "[cH1]", "benzene CH"); R.add("CB0", 1, "[cH0]", "substituted aromatic C")
    R.add("OA", 1, "[OX2H1]", "alcohol oxygen"); R.add("HA", 1, "[H][OX2]", "alcohol H"); R.add("OE", 0, "[OX2]([#6])[#6]", "ether oxygen")
    R.add("CK", 2, "[CX3](=[OX1])([#6])[#6]", "ketone carbon"); R.add("OK", 2, "[OX1]=[CX3]([#6])[#6]", "ketone oxygen")
    R.add("STI", 1, "[SX2H1]", "thiol sulfur"); R.add("HTI", 1, "[H][SX2]", "thiol H")
    R.add("SS", 1, "[SX2][SX2]", "disulfide sulfur"); R.add("S", 0, "[SX2]([#6])[#6]", "sulfide sulfur")
    R.add("CNI", 2, "[CX2]#[NX1]", "nitrile carbon"); R.add("NI", 2, "[NX1]#[CX2]", "nitrile nitrogen")
    # small gases (Potoff & Siepmann 2001): CO2, and N2 with its centre site (a virtual site CAPS adds, typed by name)
    R.add("CO2C", 3, "[CX2](=[OX1])=[OX1]", "carbon dioxide carbon"); R.add("CO2O", 3, "[OX1]=[CX2]=[OX1]", "carbon dioxide oxygen")
    R.add("N2N", 3, "[NX1]#[NX1]", "nitrogen (N2) atom"); R.add("N2M", 3, "[#0]", "nitrogen (N2) centre site", atom_name="N2M")
    return R


# ---------------------------------------------------------------- TraPPE-EH (explicit-hydrogen aromatics)
def trappe_eh():
    """TraPPE-EH (Rai & Siepmann 2007): each ring has its own keys by position (symmetry-equivalent positions share
    parameters; the library still numbers them apart, so the lower number is given)"""
    R = Rules()
    BENZ = "c;r6;!$(c~[!#6;!#1]);!$(c:[!#6])"
    R.add("C_1", 0, f"[{BENZ}]", "benzene carbon"); R.add("H_1", 0, "[H]c", "benzene H")
    # pyridine (one ring N)
    R.add("N1_2", 2, "[nX2;r6;!$(n:c:n);!$(n:n);!$(n:c:c:n)]", "pyridine N")
    R.add("C2_2", 2, "[c;r6;$(c:[nX2]);!$(c:[n]:c:[n]);!$(c:n:n);!$(c:n:c:c:n)]", "pyridine C2 / C6")
    R.add("C3_2", 2, "[c;r6;$(c:c:[nX2]);!$(c:[n]);!$(c:c:n:c:n);!$(c:c:n:n)]", "pyridine C3 / C5")
    R.add("C4_2", 2, "[c;r6;$(c:c:c:[nX2]);!$(c:[n]);!$(c:c:[n])]", "pyridine C4")
    R.add("HA", 2, "[H][c;r6;$(c:n),$(c:c:n),$(c:c:c:n);!$(c:n:n);!$(c:c:n:n)]", "H on a pyridine, pyrimidine or pyrazine C")
    # pyrimidine (1,3), pyrazine (1,4), pyridazine (1,2)
    R.add("N_3", 3, "[nX2;r6;$(n:c:[nX2])]", "pyrimidine N")
    R.add("C2_3", 3, "[c;r6;$(c(:[nX2]):[nX2])]", "pyrimidine C2")
    R.add("C4_3", 3, "[c;r6;$(c:[nX2]:c:[nX2]);!$(c(:n):n)]", "pyrimidine C4 / C6")
    R.add("C5_3", 3, "[c;r6;$(c:c:[nX2]:c:[nX2]);!$(c:n)]", "pyrimidine C5")
    R.add("N_4", 3, "[nX2;r6;$(n:c:c:[nX2])]", "pyrazine N")
    R.add("C2_4", 3, "[c;r6;$(c(:[nX2]):c:[nX2])]", "pyrazine C")
    R.add("N1_5", 3, "[nX2;r6;$(n:[nX2])]", "pyridazine N")
    R.add("C3_5", 3, "[c;r6;$(c:[nX2]:[nX2])]", "pyridazine C3 / C6")
    R.add("C4_5", 3, "[c;r6;$(c:c:[nX2]:[nX2]);!$(c:n)]", "pyridazine C4 / C5")
    R.add("H_5", 3, "[H][c;r6;$(c:[nX2]:[nX2]),$(c:c:[nX2]:[nX2])]", "pyridazine H")
    # phenol: its own ring keys
    R.add("O_8", 3, "[OX2H1]c", "phenol oxygen"); R.add("HO_8", 3, "[H][OX2]c", "phenol hydroxyl H"); R.add("C1_8", 3, "c[OX2H1]", "phenol C1")
    R.add("C2_8", 3, "[c;$(c:c[OX2H1])]", "phenol C2 / C6"); R.add("C3_8", 3, "[c;$(c:c:c[OX2H1])]", "phenol C3 / C5")
    R.add("C4_8", 3, "[c;$(c:c:c:c[OX2H1])]", "phenol C4")
    R.add("H_8", 3, "[H][c;$(c:c[OX2H1]),$(c:c:c[OX2H1]),$(c:c:c:c[OX2H1])]", "phenol ring H")
    # five-membered rings
    R.add("N1_9", 3, "[nX3H1;r5]", "pyrrole nitrogen"); R.add("H6_9", 3, "[H][nX3;r5]", "pyrrole N-H")
    R.add("C2_9", 3, "[c;r5;$(c:[nH1])]", "pyrrole C2 / C5"); R.add("C3_9", 3, "[c;r5;$(c:c:[nH1]);!$(c:[nH1])]", "pyrrole C3 / C4")
    R.add("H_9", 3, "[H][c;r5;$(c:[nH1]),$(c:c:[nH1])]", "pyrrole C-H")
    R.add("O1_7", 3, "[o;r5]", "furan oxygen"); R.add("C2_7", 3, "[c;r5;$(c:o)]", "furan C2 / C5"); R.add("C3_7", 3, "[c;r5;$(c:c:o);!$(c:o)]", "furan C3 / C4")
    R.add("H_7", 3, "[H][c;r5;$(c:o),$(c:c:o)]", "furan H")
    R.add("S1_6", 3, "[s;r5]", "thiophene sulfur"); R.add("C2_6", 3, "[c;r5;$(c:s)]", "thiophene C2 / C5"); R.add("C3_6", 3, "[c;r5;$(c:c:s);!$(c:s)]", "thiophene C3 / C4")
    R.add("H_6", 3, "[H][c;r5;$(c:s),$(c:c:s)]", "thiophene H")
    return R


# Terms a file lacks, by analogy (as parmchk2 and the CGenFF program fill missing terms): the same element and
# hybridisation, nearest in hydrogen count or neighbourhood. Every term found this way is listed as estimated.
ALKYL = lambda a, b, c: {a: [b, c], b: [a, c], c: [b, a]}
ANALOGIES = {
    "charmm": {**ALKYL("C1", "C2", "C3"), "CAE1": ["C1"], "CAE2": ["C2"], "CAE3": ["C3"], "CTHF": ["CAE2", "C2"], "CTHP": ["CAE2", "C2"], "CYCP": ["C2"],
               "HAE1": ["H1"], "HAE2": ["H1"], "HAE3": ["H1"], "HTHF": ["H1"], "C15": ["C16"], "C16": ["C15"], "H7": ["H8"], "H8": ["H7"],
               "N1": ["N2"], "N2": ["N1"], "N5": ["N2"], "OET": ["OE"], "OE": ["OET"]},
    "trappe-ua": {**ALKYL("CH1", "CH2", "CH3"), "CH0": ["CH1"], "CYC6": ["CH2"], "CYC5": ["CH2"], "CE1": ["CE2", "CE0"], "CE2": ["CE1"], "CE0": ["CE1"]},
}
# Entries the converted files leave out although the force field defines them
EXTENSIONS = {}

SETS = {
    "charmm": (charmm(), "CHARMM (library set)", False),
    "trappe-ua": (trappe_ua(), "TraPPE-UA", True),
    "trappe1998-moltemplate": (trappe_ua(False), "TraPPE-UA alkanes", True),
    "trappe-eh": (trappe_eh(), "TraPPE-EH (explicit-hydrogen aromatics)", False),
}

# water-model overlays: the same three atoms, the model the user picks (caps ff type --typing FF.typing.json,OVERLAY)
OVERLAYS = {
    "gaff-tip3p-c": ("GAFF / GAFF2 (library files): CHARMM's TIP3P water (OT / HT)", [("OT", "[OX2H2]", "TIP3P (CHARMM) water O"), ("HT", "[H][OX2H2]", "TIP3P (CHARMM) water H")]),
    "gaff-tip3p-pme": ("GAFF / GAFF2 (library files): TIP3P water for PME (OTP / HTP)", [("OTP", "[OX2H2]", "TIP3P (PME) water O"), ("HTP", "[H][OX2H2]", "TIP3P (PME) water H")]),
    "gaff-spc": ("GAFF / GAFF2 (library files): SPC / SPC/E water (OSPC / HSPC)", [("OSPC", "[OX2H2]", "SPC water O"), ("HSPC", "[H][OX2H2]", "SPC water H")]),
    "opls-tip3p-c": ("OPLS (2005, CL&P): CHARMM's TIP3P water (OTC / HTC)", [("OTC", "[OX2H2]", "TIP3P (CHARMM) water O"), ("HTC", "[H][OX2H2]", "TIP3P (CHARMM) water H")]),
    "opls-spc": ("OPLS (2005, CL&P): SPC / SPC/E water (OWS / HWS)", [("OWS", "[OX2H2]", "SPC water O"), ("HWS", "[H][OX2H2]", "SPC water H")]),
    "charmm-spc": ("CHARMM: SPC / SPC/E water", [("OSPC", "[OX2H2]", "SPC water O"), ("HSPC", "[H][OX2H2]", "SPC water H")]),
    "charmm-tip3p": ("CHARMM: original TIP3P water (H3P)", [("OT", "[OX2H2]", "TIP3P water O"), ("H3P", "[H][OX2H2]", "original TIP3P water H")]),
    "charmm-tip3p-pme": ("CHARMM: TIP3P water for PME (OTP / HTP)", [("OTP", "[OX2H2]", "TIP3P (PME) water O"), ("HTP", "[H][OX2H2]", "TIP3P (PME) water H")]),
}

if __name__ == "__main__":
    for name, (label, rules) in OVERLAYS.items():
        doc = {"format": "caps-typing", "version": 1, "forcefield": label,
               "description": f"{label}: an overlay on the force field's own rules that picks this water model (bench/typing/make_biomolecular_rules.py).",
               "rules": [{"type": t, "smarts": sm, "priority": 9, "description": d} for t, sm, d in rules]}
        json.dump(doc, open(os.path.join(DATA, "typing", name + ".typing.json"), "w"), ensure_ascii=False, indent=1)
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
