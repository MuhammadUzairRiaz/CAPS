#!/usr/bin/env python3
"""OPLS-DES for choline chloride deep eutectic solvents, from its paper and Supporting Information only:

  B. Doherty, O. Acevedo, "OPLS Force Field for Choline Chloride-Based Deep Eutectic Solvents",
  J. Phys. Chem. B 122, 9982-9993 (2018), doi 10.1021/acs.jpcb.8b06647, and its SI (jp8b06647_si_001).

What the paper says the force field is (pp. 9983-9984):
  - OPLS-AA functional form (eqs 1-4): harmonic bonds k (r - r0)^2 and bends k (theta - theta0)^2 (no 1/2), the OPLS
    Fourier torsion 1/2 V1 (1 + cos phi) + 1/2 V2 (1 - cos 2 phi) + 1/2 V3 (1 + cos 3 phi), Coulomb plus 12-6
    Lennard-Jones, geometric combining rules for sigma and epsilon, 1-2 and 1-3 pairs excluded, 1-4 pairs scaled by 0.5.
  - Charges built for a +-0.8 e scaling of the ions: the choline cation sums to +0.8, chloride is -0.8, every HBD neutral.
    The tables below already hold the scaled charges; nothing is scaled here.
  - Choline chloride's charges and Lennard-Jones terms refitted to the ND/EPSR radial distribution functions of reline;
    choline's geometry and torsions from OPLS-2009IL (Sambasivarao & Acevedo, J. Chem. Theory Comput. 5, 1038 (2009),
    the SI's ref. 1), refitted to LMP2 minima.
  - The HBDs' bonded terms are OPLS-AA's (Jorgensen, Maxwell & Tirado-Rives 1996, the SI's ref. 2), ethylene glycol's
    torsions Szefczyk & Cordeiro's (J. Phys. Chem. B 115, 3013 (2011), the SI's ref. 3); their non-bonded terms are tuned
    against choline chloride and are meant for DES mixtures only, not for the pure liquids.
  - GROMACS 5.0.7, PME, short-range electrostatic cut-off 16 A.

Every number below is transcribed from the SI's tables; each block names its table. The SI prints the HBDs' bonded terms
itself (Tables S5, S6), so none is taken from another file, with these exceptions and readings:
  - Ethylene glycol: Table S6 (footnote a, ref. 3) lists only OH-CG-CG-OH and CG-CG-OH-HO, where every other molecule's
    block lists all its torsions, zero ones included. Ref. 3 (Szefczyk & Cordeiro 2011) confirms these are the model's
    only torsions: its Table 2 gives OG-CG-CG-OG and CG-CG-OG-HO alone, and p. 3016 says "only three dihedral angles
    are present, one for the O-C-C-O angle and two for the C-C-O-H angles"; they were fitted to RHF/6-31G(d) scans with
    no other torsion in the force field. So the hydrogen torsions HG-CG-CG-HG, HG-CG-CG-OG and HG-CG-OG-HO are zero
    (OPLS-AA's 0.300, 0.468 and 0.352 kcal/mol, as the SI's glycerol rows hold them, are not added). The SI's two rows
    are not ref. 3's numbers (Table 2 in kJ/mol: OG-CG-CG-OG 12.521 -2.736 11.800, i.e. 2.993 -0.654 2.820 kcal/mol, for
    the SI's 3.887 -1.192 3.206; CG-CG-OG-HO 1.656 -3.139 4.300, i.e. 0.396 -0.750 1.028, for 0.413 -0.754 1.028); the
    SI's are kept, being the set the paper ran.
  - Ethylene glycol's rows use other names than Table S4 and Scheme S2: Table S5's bonds are in OPLS-AA class names
    (OG-HO, CT-OG, CT-CT, HC-CT) and Table S6 and Scheme S2 call the O "OH". Read as Table S4's OG, CG, HG, HO.
  - Levulinic acid: Table S5 prints the bend "HO-CB-CZ" (108.0, 70); HO is not bonded to CB. Read as OH-CB-CZ, OPLS-AA's
    C-C-OH bend (the same 108.0 / 70 as malonic acid's CT-CD-OH).
  - The rows X-NA-X-X, X-CD-X-X, X-C-X-X and X-N-X-X of Tables S3 and S6 have their centre in position 2 and are read
    as impropers on that atom, V2 / 2 (1 - cos 2 phi), over every triple of its neighbours (one on a three-connected
    atom, four on choline's N). X-CM-X-X and X-CT-X-X are all zero and are left out.
  - No improper is listed for phenol's ring or for levulinic acid's carbonyls; none is added.
  - Urea: Scheme S2 types the H cis to O as HC and the trans one as HT; Tables S4-S6 give both the same charge,
    Lennard-Jones terms, bonds, bends and torsions. CAPS's typing cannot tell them apart (graph rules) and assigns HT
    to all four urea H; HC is defined with its own terms, and the H-N-H bend (HC-N-HT) is matched as H?-N-H?.
Types carry a suffix per HBD (the SI reuses names such as CB, HO, OH, CD with different parameters); choline's and
chloride's keep the SI's names. Masses are the elements' standard masses (the paper gives none).

usage: make_opls_des_doherty2018.py [DATA_DIR]   (writes DATA_DIR/forcefields/oplsdes-doherty2018.json)
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
FID = "oplsdes-doherty2018"
REF = "B. Doherty, O. Acevedo, J. Phys. Chem. B 122, 9982 (2018)"
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "Cl": 35.453}

# ---------------------------------------------------------------------------------------------------------------------
# Molecules: suffix for its type names, name, SI tables. Choline and chloride keep the SI's names (suffix "").
# Non-bonded: (type, element, q (e), sigma (A), epsilon (kcal/mol), description)
MOLS = {}
# Table S1 (choline cation and chloride anion), Scheme S1
MOLS["choline"] = ("", "Table S1", [
    ("CS", "C", -0.131, 3.5, 0.066, "choline CH2 on N"),
    ("HS", "H", 0.068, 2.6, 0.03, "H on choline's CH2 on N"),
    ("NA", "N", 0.791, 3.25, 0.17, "choline N"),
    ("CA", "C", -0.100, 3.5, 0.066, "choline methyl C"),
    ("CW", "C", 0.132, 3.5, 0.066, "choline CH2 on O"),
    ("OY", "O", -0.468, 3.07, 0.17, "choline hydroxyl O"),
    ("HA", "H", 0.033, 2.5, 0.03, "H on choline's methyls"),
    ("HW", "H", 0.034, 2.2, 0.03, "H on choline's CH2 on O"),
    ("HY", "H", 0.275, 0.0, 0.0, "choline hydroxyl H")])
MOLS["chloride"] = ("", "Table S1", [("Cl", "Cl", -0.800, 3.77, 0.148, "chloride")])
# Table S4 (non-bonded parameters of the HBDs), Scheme S2
MOLS["ethylene glycol"] = ("_etg", "Table S4", [
    ("HO", "H", 0.348, 0.00, 0.000, "ethylene glycol hydroxyl H"),
    ("OG", "O", -0.560, 3.00, 0.2975, "ethylene glycol O (OH in Scheme S2)"),
    ("CG", "C", 0.116, 3.50, 0.1155, "ethylene glycol C"),
    ("HG", "H", 0.048, 2.50, 0.0525, "H on ethylene glycol's C")])
MOLS["glycerol"] = ("_gly", "Table S4", [
    ("CM", "C", 0.1420, 3.50, 0.1452, "glycerol central CH"),
    ("OM", "O", -0.5447, 3.07, 0.3740, "glycerol O on the central C"),
    ("CB", "C", 0.1600, 3.50, 0.1452, "glycerol terminal CH2"),
    ("HM", "H", 0.2912, 0.00, 0.0000, "H on glycerol's central O"),
    ("HZ", "H", 0.0221, 2.50, 0.0660, "H on glycerol's central C"),
    ("HC", "H", 0.0637, 2.50, 0.0660, "H on glycerol's terminal C"),
    ("OH", "O", -0.5470, 3.07, 0.3740, "glycerol O on a terminal C"),
    ("HO", "H", 0.3043, 0.00, 0.0000, "H on glycerol's terminal O")])
MOLS["levulinic acid"] = ("_lev", "Table S4", [
    ("CB", "C", 0.416, 3.75, 0.1260, "levulinic acid carboxyl C"),
    ("OH", "O", -0.424, 3.00, 0.2040, "levulinic acid carboxyl OH"),
    ("HO", "H", 0.360, 0.00, 0.0000, "levulinic acid carboxyl H"),
    ("OB", "O", -0.352, 2.96, 0.2520, "levulinic acid carboxyl =O"),
    ("CZ", "C", -0.096, 3.50, 0.0792, "levulinic acid CH2 on the carboxyl C"),
    ("HZ", "H", 0.048, 2.50, 0.0360, "H on CZ"),
    ("CF", "C", -0.096, 3.50, 0.0792, "levulinic acid CH2 on the ketone C"),
    ("HF", "H", 0.048, 2.50, 0.0360, "H on CF"),
    ("CD", "C", 0.376, 3.75, 0.1260, "levulinic acid ketone C"),
    ("CT", "C", -0.144, 3.50, 0.0792, "levulinic acid methyl C"),
    ("HT", "H", 0.048, 2.50, 0.0360, "H on the methyl"),
    ("OC", "O", -0.376, 2.96, 0.2520, "levulinic acid ketone O")])
MOLS["malonic acid"] = ("_mal", "Table S4", [
    ("CD", "C", 0.416, 3.75, 0.2625, "malonic acid carboxyl C"),
    ("CT", "C", -0.096, 3.50, 0.1650, "malonic acid CH2"),
    ("HC", "H", 0.048, 2.50, 0.0750, "H on malonic acid's CH2"),
    ("OD", "O", -0.352, 2.96, 0.5250, "malonic acid carboxyl =O"),
    ("OH", "O", -0.424, 3.00, 0.4250, "malonic acid carboxyl OH"),
    ("HO", "H", 0.360, 0.00, 0.0000, "malonic acid carboxyl H")])
MOLS["oxalic acid"] = ("_ox", "Table S4", [
    ("CD", "C", 0.416, 3.75, 0.1575, "oxalic acid C"),
    ("OH", "O", -0.394, 2.92, 0.2550, "oxalic acid OH"),
    ("HO", "H", 0.330, 0.00, 0.0000, "oxalic acid H"),
    ("OD", "O", -0.352, 2.96, 0.3150, "oxalic acid =O")])
MOLS["phenol"] = ("_phe", "Table S4", [
    ("CB", "C", -0.092, 3.55, 0.077, "phenol ring CH"),
    ("CO", "C", 0.120, 3.55, 0.077, "phenol ring C holding the OH"),
    ("HB", "H", 0.092, 2.42, 0.033, "phenol ring H"),
    ("OH", "O", -0.468, 3.07, 0.187, "phenol O"),
    ("HO", "H", 0.348, 0.00, 0.000, "phenol hydroxyl H")])
MOLS["urea"] = ("_urea", "Table S4", [
    ("C", "C", 0.124, 3.75, 0.1575, "urea C"),
    ("O", "O", -0.322, 2.96, 0.3150, "urea O"),
    ("N", "N", -0.453, 3.55, 0.2550, "urea N"),
    ("HT", "H", 0.276, 0.00, 0.0000, "urea H trans to O (Scheme S2); CAPS assigns it to every urea H"),
    ("HC", "H", 0.276, 0.00, 0.0000, "urea H cis to O (Scheme S2); the same terms as HT, not assigned by CAPS's rules")])

# ---------------------------------------------------------------------------------------------------------------------
# Bonds (r0 A, kr kcal/mol/A^2) and bends (theta0 deg, k kcal/mol/rad^2), in the SI's names (Tables S2, S5)
BONDS = {
    "choline": ("Table S2", [("HA-CA", 1.099, 340), ("CA-NA", 1.498, 490), ("NA-CS", 1.516, 490), ("CS-HS", 1.0805, 340),
                             ("CS-CW", 1.521, 317), ("CW-HW", 1.085, 340), ("CW-OY", 1.395, 450), ("OY-HY", 0.949, 553)]),
    # Table S5 prints these four in OPLS-AA class names: OG-HO, CT-OG, CT-CT, HC-CT
    "ethylene glycol": ("Table S5 (OG-HO, CT-OG, CT-CT, HC-CT read as Table S4's OG, CG, HG)",
                        [("OG-HO", 0.945, 553), ("CG-OG", 1.410, 320), ("CG-CG", 1.529, 268), ("HG-CG", 1.090, 340)]),
    "glycerol": ("Table S5", [("OM-CM", 1.410, 320), ("CB-CM", 1.529, 268), ("HM-OM", 0.945, 553), ("HZ-CM", 1.090, 340),
                              ("HC-CB", 1.090, 340), ("OH-CB", 1.410, 320), ("HO-OH", 0.945, 553)]),
    "levulinic acid": ("Table S5", [("OH-CB", 1.364, 450), ("OH-HO", 0.945, 553), ("OB-CB", 1.229, 570), ("CZ-CB", 1.522, 317),
                                    ("HZ-CZ", 1.090, 340), ("CF-CZ", 1.529, 268), ("HF-CF", 1.090, 340), ("CD-CF", 1.522, 317),
                                    ("CT-CD", 1.522, 317), ("HT-CT", 1.090, 340), ("OC-CD", 1.229, 570)]),
    "malonic acid": ("Table S5", [("OD-CD", 1.229, 570), ("OH-CD", 1.364, 450), ("CT-CD", 1.522, 317), ("HO-OH", 0.945, 553),
                                  ("HC-CT", 1.090, 340)]),
    "oxalic acid": ("Table S5", [("CD-CD", 1.510, 350), ("OD-CD", 1.229, 570), ("OH-CD", 1.364, 450), ("HO-OH", 0.945, 553)]),
    "phenol": ("Table S5", [("CB-CB", 1.400, 469), ("CO-CB", 1.400, 469), ("HB-CB", 1.080, 367), ("OH-CO", 1.365, 450),
                            ("HO-OH", 0.945, 553)]),
    "urea": ("Table S5", [("O-C", 1.229, 570), ("N-C", 1.335, 490), ("HT-N", 1.010, 434), ("HC-N", 1.010, 434)]),
}
ANGLES = {
    "choline": ("Table S2", [("HA-CA-HA", 110.01, 35), ("HA-CA-NA", 108.9, 35), ("CA-NA-CS", 110.2, 51.8),
                             ("NA-CS-HS", 106.4, 35), ("NA-CS-CW", 116.6, 70), ("CS-CW-HW", 108.3, 35),
                             ("CS-CW-OY", 109.6, 80), ("HS-CS-CW", 109.3, 35), ("HW-CW-OY", 111.6, 35),
                             ("CW-OY-HY", 110.9, 35), ("CA-NA-CA", 108.73, 55), ("HS-CS-HS", 108.6, 35),
                             ("HW-CW-HW", 107.4, 35)]),
    "ethylene glycol": ("Table S5", [("HO-OG-CG", 108.5, 55), ("OG-CG-CG", 108.0, 50), ("OG-CG-HG", 109.5, 35),
                                     ("CG-CG-HG", 110.7, 37.5), ("HG-CG-HG", 107.8, 33)]),
    "glycerol": ("Table S5", [("OM-CM-CB", 109.5, 50), ("CM-OM-HM", 108.5, 55), ("OM-CM-HZ", 109.5, 35),
                              ("CM-CB-HC", 110.7, 37.5), ("CM-CB-OH", 109.5, 50), ("CB-OH-HO", 108.5, 55),
                              ("CB-CM-CB", 112.7, 58.35), ("CB-CM-HZ", 110.7, 37.5), ("HC-CB-HC", 107.8, 33),
                              ("HC-CB-OH", 109.5, 35)]),
    # Table S5 prints "HO-CB-CZ 108.0 70"; HO is not bonded to CB: read as OH-CB-CZ
    "levulinic acid": ("Table S5 (its HO-CB-CZ read as OH-CB-CZ)",
                       [("CB-OH-HO", 113.0, 35), ("OH-CB-OB", 121.0, 80), ("OH-CB-CZ", 108.0, 70), ("CB-CZ-HZ", 109.5, 35),
                        ("CB-CZ-CF", 111.1, 63), ("CZ-CF-HF", 110.7, 37.5), ("CZ-CF-CD", 111.1, 63), ("CF-CD-CT", 116.0, 70),
                        ("CD-CT-HT", 109.5, 35), ("CF-CD-OC", 120.4, 80), ("OB-CB-CZ", 120.4, 80), ("HZ-CZ-HZ", 107.8, 33),
                        ("HZ-CZ-CF", 110.7, 37.5), ("HF-CF-HF", 107.8, 33), ("HF-CF-CD", 109.5, 35), ("CT-CD-OC", 120.4, 80),
                        ("HT-CT-HT", 107.8, 33)]),
    "malonic acid": ("Table S5", [("OD-CD-OH", 121.0, 80), ("OD-CD-CT", 120.4, 80), ("CD-OH-HO", 113.0, 35),
                                  ("CD-CT-CD", 111.1, 63), ("CD-CT-HC", 109.5, 35), ("CT-CD-OH", 108.0, 70),
                                  ("HC-CT-HC", 107.8, 33)]),
    "oxalic acid": ("Table S5", [("CD-CD-OD", 121.40, 80), ("CD-CD-OH", 118.03, 70.96), ("CD-OH-HO", 113.00, 35),
                                 ("OH-CD-OD", 121.00, 80)]),
    "phenol": ("Table S5", [("CB-CB-CB", 120.0, 63), ("CO-CB-CB", 120.0, 63), ("CB-CO-OH", 120.0, 70), ("CB-CB-HB", 120.0, 35),
                            ("CO-CB-HB", 120.0, 35), ("CO-OH-HO", 113.0, 35), ("CB-CO-CB", 120.0, 63)]),
    "urea": ("Table S5", [("O-C-N", 122.9, 80), ("C-N-HC", 119.8, 35), ("C-N-HT", 119.8, 35), ("N-C-N", 114.2, 70),
                          ("HC-N-HT", 120.0, 35)]),
}
# Torsions V1 V2 V3 (kcal/mol), Tables S3 and S6; rows "X-c-X-X" are the impropers (centre c)
TORSIONS = {
    "choline": ("Table S3", [("CW-CS-NA-CA", 0.100, 0.550, 0.650), ("CA-NA-CA-HA", 0.000, 0.000, 0.825),
                             ("CS-NA-CA-HA", 0.000, 0.000, 0.940), ("HS-CS-NA-CA", 0.000, 1.000, 0.700),
                             ("OY-CW-CS-NA", -6.000, -5.000, 3.200), ("OY-CW-CS-HS", -0.500, -2.500, 0.250),
                             ("HW-CW-CS-NA", -6.000, -7.000, 0.750), ("HW-CW-CS-HS", 6.000, -3.000, 2.000),
                             ("HY-OY-CW-CS", -0.356, -0.174, 0.350), ("HY-OY-CW-HW", -3.000, 1.000, -2.000),
                             ("X-NA-X-X", 0.000, 2.000, 0.000)]),
    # Table S6 (ref. 3) prints OH-CG-CG-OH and CG-CG-OH-HO: OH read as OG
    "ethylene glycol": ("Table S6 (ref. 3; its OH read as OG)", [("OG-CG-CG-OG", 3.887, -1.192, 3.206),
                                                                 ("CG-CG-OG-HO", 0.413, -0.754, 1.028)]),
    "glycerol": ("Table S6", [("HM-OM-CM-CB", -0.356, -0.174, 0.492), ("HZ-CM-OM-HM", 0.000, 0.000, 0.352),
                              ("HC-CB-CM-OM", 0.000, 0.000, 0.468), ("OH-CB-CM-OM", 12.234, 0.000, 0.000),
                              ("HO-OH-CB-CM", -0.356, -0.174, 0.492), ("HC-CB-CM-CB", 0.000, 0.000, 0.300),
                              ("OH-CB-CM-CB", -1.552, 0.000, 0.000), ("HC-CB-CM-HZ", 0.000, 0.000, 0.300),
                              ("OH-CB-CM-HZ", 0.000, 0.000, 0.468), ("HO-OH-CB-HC", 0.000, 0.000, 0.352)]),
    # X-CM-X-X 0 0 0 (Table S6): all zero, left out
    "levulinic acid": ("Table S6", [("OB-CB-OH-HO", 0.000, 5.500, 0.000), ("CZ-CB-OH-HO", 1.500, 5.500, 0.000),
                                    ("HZ-CZ-CB-OB", 0.000, 0.000, 0.000), ("CF-CZ-CB-OB", 0.000, 0.000, 0.000),
                                    ("HF-CF-CZ-CB", 0.000, 0.000, -0.100), ("CD-CF-CZ-CB", -0.550, 0.000, 1.000),
                                    ("CT-CD-CF-CZ", 1.454, -0.144, -0.775), ("HT-CT-CD-CF", 0.000, 0.000, 0.275),
                                    ("OC-CD-CF-CZ", 0.000, 0.000, 0.000), ("HZ-CZ-CB-OH", 0.000, 0.000, 0.000),
                                    ("CF-CZ-CB-OH", 1.000, 0.546, 0.450), ("HF-CF-CZ-HZ", 0.000, 0.000, 0.300),
                                    ("CD-CF-CZ-HZ", 0.000, 0.000, -0.100), ("CT-CD-CF-HF", 0.000, 0.000, 0.275),
                                    ("OC-CD-CF-HF", 0.000, 0.000, 0.000), ("OC-CD-CT-HT", 0.000, 0.000, 0.000)]),
    "malonic acid": ("Table S6", [("HO-OH-CD-CT", 1.500, 5.500, 0.000), ("HC-CT-CD-OD", 0.000, 0.000, 0.000),
                                  ("CD-CT-CD-OH", 1.000, 0.546, 0.450), ("HC-CT-CD-OH", 0.000, 0.000, 0.000),
                                  ("OD-CD-OH-HO", 0.000, 5.500, 0.000), ("OD-CD-CT-CD", 0.000, 0.000, 0.000),
                                  ("X-CD-X-X", 0.000, 21.000, 0.000)]),
    # X-CT-X-X 0 0 0 (Table S6): all zero, left out
    "oxalic acid": ("Table S6", [("HO-OH-CD-OD", 0.000, 5.500, 0.000), ("OD-CD-CD-OD", 1.600, 3.200, 0.000),
                                 ("OH-CD-CD-OH", 1.600, 3.200, 0.000), ("OD-CD-CD-OH", 1.600, 3.200, 0.000),
                                 ("HO-OH-CD-CD", 3.000, 5.500, 0.000), ("X-CD-X-X", 0.000, 21.000, 0.000)]),
    "phenol": ("Table S6", [("HO-OH-CO-CB", 0.000, 1.682, 0.000), ("CB-CB-CB-CB", 0.000, 7.250, 0.000),
                            ("CO-CB-CB-CB", 0.000, 7.250, 0.000), ("CB-CB-CO-CB", 0.000, 7.250, 0.000),
                            ("CB-CB-CB-HB", 0.000, 7.250, 0.000), ("CO-CB-CB-HB", 0.000, 7.250, 0.000),
                            ("CB-CO-CB-HB", 0.000, 7.250, 0.000), ("OH-CO-CB-CB", 0.000, 7.250, 0.000),
                            ("OH-CO-CB-HB", 0.000, 7.250, 0.000), ("HB-CB-CB-HB", 0.000, 7.250, 0.000)]),
    "urea": ("Table S6", [("HT-N-C-O", 0.000, 4.900, 0.000), ("HC-N-C-O", 0.000, 4.900, 0.000),
                          ("HT-N-C-N", 0.000, 4.900, 0.000), ("HC-N-C-N", 0.000, 4.900, 0.000),
                          ("X-C-X-X", 0.000, 21.000, 0.000), ("X-N-X-X", 0.000, 5.000, 0.000)]),
}
# Ethylene glycol's hydrogen torsions: Table S6 lists none, and ref. 3 (Szefczyk & Cordeiro 2011, Table 2 and p. 3016)
# has no torsion but O-C-C-O and C-C-O-H: zero
ZERO = {"ethylene glycol": ["HG-CG-CG-HG", "HG-CG-CG-OG", "HG-CG-OG-HO"]}


def main():
    types, pairs, bonds, angles, dihedrals, impropers = [], [], [], [], [], []
    names = {}
    for mol, (sfx, table, rows) in MOLS.items():
        local = {}
        for t, el, q, sig, eps, desc in rows:
            name = t + sfx
            assert name not in names, name
            names[name] = mol
            local[t] = name
            types.append({"name": name, "element": el, "mass": MASS[el], "charge": q,
                          "description": f"{desc} ({table}, Scheme {'S1' if mol in ('choline', 'chloride') else 'S2'})"})
            pairs.append({"name": name, "match": [name], "params": [eps, sig],
                          "comment": f"sigma {sig} A, epsilon {eps} kcal/mol ({table})"})
        MOLS[mol] = (sfx, table, rows, local)

    def tn(mol, t):
        return "*" if t == "X" else MOLS[mol][3][t]

    for mol, (table, rows) in BONDS.items():
        for lab, r0, k in rows:
            m = [tn(mol, t) for t in lab.split("-")]
            bonds.append({"name": "-".join(m), "match": m, "style": "harmonic", "params": [float(k), r0],
                          "comment": f"{mol} {lab}: kr {k} kcal/mol/A^2, r0 {r0} A ({table})"})
    for mol, (table, rows) in ANGLES.items():
        for lab, th, k in rows:
            m = [tn(mol, t) for t in lab.split("-")]
            note = ""
            if mol == "urea" and lab == "HC-N-HT":   # each N holds one HC and one HT: H-N-H is always this bend
                m, note = ["H?_urea", "N_urea", "H?_urea"], "; matched as H?-N-H? (CAPS types every urea H HT)"
            angles.append({"name": "-".join(m), "match": m, "style": "harmonic", "params": [float(k), th],
                           "comment": f"{mol} {lab}: k {k} kcal/mol/rad^2, theta0 {th} deg ({table}){note}"})
    for mol, (table, rows) in TORSIONS.items():
        for lab, v1, v2, v3 in rows:
            parts = lab.split("-")
            m = [tn(mol, t) for t in parts]
            if parts[0] == "X":   # an improper, centre in position 2: V2 / 2 (1 - cos 2 phi) = K [1 + d cos(n phi)]
                assert v1 == 0 and v3 == 0 and parts[2:] == ["X", "X"]
                impropers.append({"name": "-".join(m), "match": m, "style": "cvff", "params": [v2 / 2, -1, 2],
                                  "comment": f"{mol} {lab}: V2 {v2} kcal/mol, improper on {parts[1]} ({table})"})
            else:
                dihedrals.append({"name": "-".join(m), "match": m, "style": "opls", "params": [v1, v2, v3, 0.0],
                                  "comment": f"{mol} {lab}: V1 V2 V3 {v1} {v2} {v3} kcal/mol ({table})"})
    for mol, rows in ZERO.items():
        for lab in rows:
            m = [tn(mol, t) for t in lab.split("-")]
            dihedrals.append({"name": "-".join(m), "match": m, "style": "opls", "params": [0.0, 0.0, 0.0, 0.0],
                              "comment": f"{mol} {lab}: zero; Table S6 lists only the two heavy-atom torsions of ref. 3, "
                                         "and ref. 3 (Szefczyk & Cordeiro 2011, Table 2, p. 3016) has no other torsion"})
    doc = {
        "format": "caps-forcefield", "format_version": 1,
        "name": "OPLS-DES (Doherty & Acevedo 2018)",
        "source": f"{REF}, Supporting Information Tables S1-S6 and Schemes S1-S2",
        "references": [REF, "S. V. Sambasivarao, O. Acevedo, J. Chem. Theory Comput. 5, 1038 (2009)",
                       "W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)",
                       "B. Szefczyk, M. N. D. S. Cordeiro, J. Phys. Chem. B 115, 3013 (2011)"],
        "version": "2018 (SI Tables S1-S6)",
        "units": "real",
        "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"},
        "mixing": "geometric", "special_lj": [0, 0, 0.5], "special_coul": [0, 0, 0.5], "cutoff": 16, "ion_charge_scale": 0.8,
        "improper_order": "center2_sorted",
        "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": impropers,
        "notes": [
            "built by CAPS (bench/ff/make_opls_des_doherty2018.py) from the paper and its SI only; each term cites its table",
            "charges as the SI prints them, built for a +-0.8 e ion scaling: choline +0.8, chloride -0.8, every HBD neutral",
            "choline chloride's charges and Lennard-Jones terms refitted to reline's ND/EPSR structure; the HBDs' non-bonded "
            "terms tuned against choline chloride: for DES mixtures only, not the pure HBD liquids (paper p. 9984)",
            "HBD bonds, bends and torsions as the SI's Tables S5 and S6 print them (OPLS-AA's, ethylene glycol's torsions "
            "from Szefczyk & Cordeiro 2011); choline's from OPLS-2009IL (Table S2, S3)",
            "ethylene glycol: Table S6 lists only OH-CG-CG-OH and CG-CG-OH-HO; its hydrogen torsions H-C-C-H, H-C-C-O and "
            "H-C-O-H are zero, as in Szefczyk & Cordeiro 2011, whose model has only the O-C-C-O and C-C-O-H torsions (Table 2, "
            "p. 3016). Their Table 2 values (2.993 -0.654 2.820 and 0.396 -0.750 1.028 kcal/mol) differ from the SI's rows; "
            "the SI's are used",
            "rows X-c-X-X of Tables S3 and S6 read as impropers on c, V2 / 2 (1 - cos 2 phi), every triple of its neighbours "
            "(four on choline's N); the all-zero X-CM-X-X and X-CT-X-X left out; no improper for phenol's ring or "
            "levulinic acid's carbonyls (none listed)",
            "levulinic acid's bend HO-CB-CZ (Table S5) read as OH-CB-CZ; ethylene glycol's OPLS-AA class names in Table S5 "
            "(CT, HC) and the O called OH in Table S6 and Scheme S2 read as Table S4's CG, HG, OG",
            "urea's HC (cis to O) and HT (trans) have identical terms; CAPS's rules assign HT to every urea H, and the "
            "H-N-H bend (HC-N-HT, Table S5) is matched as H?-N-H?",
            "ion_charge_scale 0.8: an ion's types add up to its formal charge times 0.8 (the paper's scaling)",
            "geometric mixing, 1-4 pairs scaled by 0.5, cut-off 16 A (the paper's GROMACS set-up, PME)",
            "masses: standard atomic masses (the paper gives none)"],
        "typing": f"../typing/{FID}.typing.json",
    }
    path = os.path.join(DATA, "forcefields", FID + ".json")
    json.dump(doc, open(path, "w"), indent=1, ensure_ascii=False)
    print(f"{FID}: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} dihedrals, "
          f"{len(impropers)} impropers")


if __name__ == "__main__":
    main()
