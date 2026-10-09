// Element data for labels (generated; do not edit by hand):
//   names: IUPAC (CIAAW, Standard Atomic Weights, ciaaw.org/atomic-weights.htm, 2024).
//   mass number: the most abundant isotope in IUPAC's Isotopic Compositions of the Elements 2024 (CIAAW; the midpoint
//   of an interval) for the 84 elements it lists; for the others the longest-lived isotope: the bracketed value of
//   NIST's Atomic Weights and Isotopic Compositions (Tc, Pm, Po-Ac, Np, Pu), else the longest measured half-life of
//   NUBASE2020 (Kondev et al., Chin. Phys. C 45, 030001 (2021); Am-Og).
//   Pauling electronegativity and ground-state electron configuration: Blue Obelisk Data Repository via VTK 9.6.1 (MIT);
//   0 = not tabulated. Indexed by atomic number (0 unused).
struct ElementExtra { const char* name; double pauling; int mass_number; const char* configuration; };
static const ElementExtra kElementExtra[] = {
    {"", 0.0, 0, ""},
    {"Hydrogen", 2.20, 1, "1s1"},   // 1 H
    {"Helium", 0.00, 4, "1s2"},   // 2 He
    {"Lithium", 0.98, 7, "[He] 2s1"},   // 3 Li
    {"Beryllium", 1.57, 9, "[He] 2s2"},   // 4 Be
    {"Boron", 2.04, 11, "[He] 2s2 2p1"},   // 5 B
    {"Carbon", 2.55, 12, "[He] 2s2 2p2"},   // 6 C
    {"Nitrogen", 3.04, 14, "[He] 2s2 2p3"},   // 7 N
    {"Oxygen", 3.44, 16, "[He] 2s2 2p4"},   // 8 O
    {"Fluorine", 3.98, 19, "[He] 2s2 2p5"},   // 9 F
    {"Neon", 0.00, 20, "[He] 2s2 2p6"},   // 10 Ne
    {"Sodium", 0.93, 23, "[Ne] 3s1"},   // 11 Na
    {"Magnesium", 1.31, 24, "[Ne] 3s2"},   // 12 Mg
    {"Aluminium", 1.61, 27, "[Ne] 3s2 3p1"},   // 13 Al
    {"Silicon", 1.90, 28, "[Ne] 3s2 3p2"},   // 14 Si
    {"Phosphorus", 2.19, 31, "[Ne] 3s2 3p3"},   // 15 P
    {"Sulfur", 2.58, 32, "[Ne] 3s2 3p4"},   // 16 S
    {"Chlorine", 3.16, 35, "[Ne] 3s2 3p5"},   // 17 Cl
    {"Argon", 0.00, 40, "[Ne] 3s2 3p6"},   // 18 Ar
    {"Potassium", 0.82, 39, "[Ar] 4s1"},   // 19 K
    {"Calcium", 1.00, 40, "[Ar] 4s2"},   // 20 Ca
    {"Scandium", 1.36, 45, "[Ar] 3d1 4s2"},   // 21 Sc
    {"Titanium", 1.54, 48, "[Ar] 3d2 4s2"},   // 22 Ti
    {"Vanadium", 1.63, 51, "[Ar] 3d3 4s2"},   // 23 V
    {"Chromium", 1.66, 52, "[Ar] 3d5 4s1"},   // 24 Cr
    {"Manganese", 1.55, 55, "[Ar] 3d5 4s2"},   // 25 Mn
    {"Iron", 1.83, 56, "[Ar] 3d6 4s2"},   // 26 Fe
    {"Cobalt", 1.88, 59, "[Ar] 3d7 4s2"},   // 27 Co
    {"Nickel", 1.91, 58, "[Ar] 3d8 4s2"},   // 28 Ni
    {"Copper", 1.90, 63, "[Ar] 3d10 4s1"},   // 29 Cu
    {"Zinc", 1.65, 64, "[Ar] 3d10 4s2"},   // 30 Zn
    {"Gallium", 1.81, 69, "[Ar] 3d10 4s2 4p1"},   // 31 Ga
    {"Germanium", 2.01, 74, "[Ar] 3d10 4s2 4p2"},   // 32 Ge
    {"Arsenic", 2.18, 75, "[Ar] 3d10 4s2 4p3"},   // 33 As
    {"Selenium", 2.55, 80, "[Ar] 3d10 4s2 4p4"},   // 34 Se
    {"Bromine", 2.96, 79, "[Ar] 3d10 4s2 4p5"},   // 35 Br
    {"Krypton", 3.00, 84, "[Ar] 3d10 4s2 4p6"},   // 36 Kr
    {"Rubidium", 0.82, 85, "[Kr] 5s1"},   // 37 Rb
    {"Strontium", 0.95, 88, "[Kr] 5s2"},   // 38 Sr
    {"Yttrium", 1.22, 89, "[Kr] 4d1 5s2"},   // 39 Y
    {"Zirconium", 1.33, 90, "[Kr] 4d2 5s2"},   // 40 Zr
    {"Niobium", 1.60, 93, "[Kr] 4d4"},   // 41 Nb
    {"Molybdenum", 2.16, 98, "[Kr] 4d5 5s1"},   // 42 Mo
    {"Technetium", 1.90, 98, "[Kr] 4d6 5s1"},   // 43 Tc
    {"Ruthenium", 2.20, 102, "[Kr] 4d7 5s1"},   // 44 Ru
    {"Rhodium", 2.28, 103, "[Kr] 4d8 5s1"},   // 45 Rh
    {"Palladium", 2.20, 106, "[Kr] 4d10"},   // 46 Pd
    {"Silver", 1.93, 107, "[Kr] 4d10 5s1"},   // 47 Ag
    {"Cadmium", 1.69, 114, "[Kr] 4d10 5s2"},   // 48 Cd
    {"Indium", 1.78, 115, "[Kr] 4d10 5s2 5p1"},   // 49 In
    {"Tin", 1.96, 120, "[Kr] 4d10 5s2 5p2"},   // 50 Sn
    {"Antimony", 2.05, 121, "[Kr] 4d10 5s2 5p3"},   // 51 Sb
    {"Tellurium", 2.10, 130, "[Kr] 4d10 5s2 5p4"},   // 52 Te
    {"Iodine", 2.66, 127, "[Kr] 4d10 5s2 5p5"},   // 53 I
    {"Xenon", 2.60, 132, "[Kr] 4d10 5s2 5p6"},   // 54 Xe
    {"Caesium", 0.79, 133, "[Xe] 6s1"},   // 55 Cs
    {"Barium", 0.89, 138, "[Xe] 6s2"},   // 56 Ba
    {"Lanthanum", 1.10, 139, "[Xe] 5d1 6s2"},   // 57 La
    {"Cerium", 1.12, 140, "[Xe] 4f1 5d1 6s2"},   // 58 Ce
    {"Praseodymium", 1.13, 141, "[Xe] 4f3 6s2"},   // 59 Pr
    {"Neodymium", 1.14, 142, "[Xe] 4f4 6s2"},   // 60 Nd
    {"Promethium", 0.00, 145, "[Xe] 4f5 6s2"},   // 61 Pm
    {"Samarium", 1.17, 152, "[Xe] 4f6 6s2"},   // 62 Sm
    {"Europium", 0.00, 153, "[Xe] 4f7 6s2"},   // 63 Eu
    {"Gadolinium", 1.20, 158, "[Xe] 4f7 5d1 6s2"},   // 64 Gd
    {"Terbium", 0.00, 159, "[Xe] 4f9 6s2"},   // 65 Tb
    {"Dysprosium", 1.22, 164, "[Xe] 4f10 6s2"},   // 66 Dy
    {"Holmium", 1.23, 165, "[Xe] 4f11 6s2"},   // 67 Ho
    {"Erbium", 1.24, 166, "[Xe] 4f12 6s2"},   // 68 Er
    {"Thulium", 1.25, 169, "[Xe] 4f13 6s2"},   // 69 Tm
    {"Ytterbium", 0.00, 174, "[Xe] 4f14 6s2"},   // 70 Yb
    {"Lutetium", 1.27, 175, "[Xe] 4f14 5d1 6s2"},   // 71 Lu
    {"Hafnium", 1.30, 180, "[Xe] 4f14 5d2 6s2"},   // 72 Hf
    {"Tantalum", 1.50, 181, "[Xe] 4f14 5d3 6s2"},   // 73 Ta
    {"Tungsten", 2.36, 184, "[Xe] 4f14 5d4 6s2"},   // 74 W
    {"Rhenium", 1.90, 187, "[Xe] 4f14 5d5 6s2"},   // 75 Re
    {"Osmium", 2.20, 192, "[Xe] 4f14 5d6 6s2"},   // 76 Os
    {"Iridium", 2.20, 193, "[Xe] 4f14 5d7 6s2"},   // 77 Ir
    {"Platinum", 2.28, 195, "[Xe] 4f14 5d9 6s1"},   // 78 Pt
    {"Gold", 2.54, 197, "[Xe] 4f14 5d10 6s1"},   // 79 Au
    {"Mercury", 2.00, 202, "[Xe] 4f14 5d10 6s2"},   // 80 Hg
    {"Thallium", 1.62, 205, "[Xe] 4f14 5d10 6s2 6p1"},   // 81 Tl
    {"Lead", 2.33, 208, "[Xe] 4f14 5d10 6s2 6p2"},   // 82 Pb
    {"Bismuth", 2.02, 209, "[Xe] 4f14 5d10 6s2 6p3"},   // 83 Bi
    {"Polonium", 2.00, 209, "[Xe] 4f14 5d10 6s2 6p4"},   // 84 Po
    {"Astatine", 2.20, 210, "[Xe] 4f14 5d10 6s2 6p5"},   // 85 At
    {"Radon", 0.00, 222, "[Xe] 4f14 5d10 6s2 6p6"},   // 86 Rn
    {"Francium", 0.70, 223, "[Rn] 7s1"},   // 87 Fr
    {"Radium", 0.90, 226, "[Rn] 7s2"},   // 88 Ra
    {"Actinium", 1.10, 227, "[Rn] 6d1 7s2"},   // 89 Ac
    {"Thorium", 1.30, 232, "[Rn] 6d2 7s2"},   // 90 Th
    {"Protactinium", 1.50, 231, "[Rn] 5f2 6d1 7s2"},   // 91 Pa
    {"Uranium", 1.38, 238, "[Rn] 5f3 6d1 7s2"},   // 92 U
    {"Neptunium", 1.36, 237, "[Rn] 5f4 6d1 7s2"},   // 93 Np
    {"Plutonium", 1.28, 244, "[Rn] 5f6 7s2"},   // 94 Pu
    {"Americium", 1.30, 243, "[Rn] 5f7 7s2"},   // 95 Am
    {"Curium", 1.30, 247, "[Rn] 5f7 6d1 7s2"},   // 96 Cm
    {"Berkelium", 1.30, 247, "[Rn] 5f9 7s2"},   // 97 Bk
    {"Californium", 1.30, 251, "[Rn] 5f10 7s2"},   // 98 Cf
    {"Einsteinium", 1.30, 252, "[Rn] 5f11 7s2"},   // 99 Es
    {"Fermium", 1.30, 257, "[Rn] 5f12 7s2"},   // 100 Fm
    {"Mendelevium", 1.30, 258, "[Rn] 5f13 7s2"},   // 101 Md
    {"Nobelium", 1.30, 259, "[Rn] 5f14 7s2"},   // 102 No
    {"Lawrencium", 0.00, 266, "[Rn] 5f14 7s2 7p1"},   // 103 Lr
    {"Rutherfordium", 0.00, 267, "[Rn] 5f14 6d2 7s2"},   // 104 Rf
    {"Dubnium", 0.00, 268, "[Rn] 5f14 6d3 7s2"},   // 105 Db
    {"Seaborgium", 0.00, 269, ""},   // 106 Sg
    {"Bohrium", 0.00, 270, ""},   // 107 Bh
    {"Hassium", 0.00, 269, ""},   // 108 Hs
    {"Meitnerium", 0.00, 277, ""},   // 109 Mt
    {"Darmstadtium", 0.00, 282, ""},   // 110 Ds
    {"Roentgenium", 0.00, 282, ""},   // 111 Rg
    {"Copernicium", 0.00, 285, ""},   // 112 Cn
    {"Nihonium", 0.00, 286, ""},   // 113 Nh
    {"Flerovium", 0.00, 290, ""},   // 114 Fl
    {"Moscovium", 0.00, 290, ""},   // 115 Mc
    {"Livermorium", 0.00, 293, ""},   // 116 Lv
    {"Tennessine", 0.00, 294, ""},   // 117 Ts
    {"Oganesson", 0.00, 295, ""},   // 118 Og
};
