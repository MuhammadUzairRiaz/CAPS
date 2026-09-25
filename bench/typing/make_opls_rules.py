#!/usr/bin/env python3
"""Write CAPS typing rules for OPLS-AA as DL_FIELD's OPLS2005 library names its types (data/typing/opls2005).

Each type gets SMARTS for the chemistry its DL_FIELD atom types describe (the ATOM_TYPE table of OPLS2005.sf maps
e.g. C_ester, CR_ester, C_ketone to CO4). Generic element rules have priority 0; more specific environments higher
priorities. Types DL_FIELD reserves for single residues or molecules (amino-acid side-chain CB atoms, sugars, nucleic
bases) are not written as rules: those atoms get the generic type (CT, OAL, ...), and a user overlay can give the
residue-specific type.

usage: make_opls_rules.py OUTDIR
"""
import json, os, sys

X = "[F,Cl,Br,I]"
R = []


def add(t, p, *smarts, d=""):
    for s in smarts:
        R.append({"type": t, "smarts": s, **({"priority": p} if p else {}), **({"description": d} if d else {})})


def halo(n):
    """a carbon with exactly n halogens"""
    has = lambda k: "$(*" + "".join(f"({X})" for _ in range(k - 1)) + X + ")"
    return has(n) + (f";!{has(n + 1)}" if n < 4 else "")


# ---- hydrogen
add("HC", 0, "[H][CX4]", d="H on sp3 C")
add("HA", 1, "[H]c", d="H on aromatic C")
add("HE", 1, "[H][CX3]", d="H on alkene / imine C")
add("HAE", 2, "[H][CX3]=O", d="aldehyde H")
add("HALK", 3, "[H][CX4][OX1-]", d="H on alkoxide C")
add("HMET", 4, "[H][CX4H3][OX1-]", d="methoxide H")
add("HY", 2, "[H]C#C", d="alkyne H")
add("H", 0, "[H][#7]", d="H on N")
add("HAM", 2, "[H]N(C=O)C=O", d="imide H")
add("HNA", 3, "[H][NX3H3]", d="ammonia H")
add("HNP", 3, "[H][n+]", d="pyridinium H")
add("HO", 0, "[H][#8]", d="H on O")
add("HHY", 3, "[H][OX3+]", d="hydronium H")
add("H_OH", 3, "[H][OX1-]", d="hydroxide H")
add("HT3", 5, "[H][OX2H2]", d="water H (TIP3P)")
add("HS", 0, "[H][#16]")
add("HP", 0, "[H][#15]")
add("HSi", 0, "[H][Si]")

# ---- carbon
add("CT", 0, "[CX4]", d="sp3 C")
add("CYA", 1, "[CX4;R;!$(*~[#7;R]);!$(*~[#8;R;!r5])]", d="sp3 C in a carbocycle or tetrahydrofuran")
ring_n = ",".join(["$(*1~[#7]~*~*~*~*1)", "$(*1~*~[#7]~*~*~*1)", "$(*1~*~*~[#7]~*~*1)", "$(*1~[#7]~*~*~*1)", "$(*1~*~[#7]~*~*1)"])
add("CTNC", 2, f"[CX4;{ring_n}]", d="sp3 C in a saturated N-heterocycle (piperidine, pyrrolidine, morpholine)")
add("C3T", 3, "[CX4;r3]", d="sp3 C in a 3-ring")
add("C4T", 3, "[CX4;r4]", d="sp3 C in a 4-ring")
add("CO", 2, "[CX4](O)O", d="acetal / ketal C")
add("CT1", 3, "[CX4H1]([#7])C=O", d="amino-acid alpha C")
add("CT1G", 3, "[CX4H2]([NX4+])C(=O)[OX1-]", d="glycine alpha C (zwitterion)")
add("CALK", 3, "[CX4][OX1-]", d="alkoxide C")
add("CTHO", 3, f"[CX4;$(*{X});$(*Oc)]", d="halogenated C of an aryl ether")
add("CA", 0, "[c]", d="aromatic C")
add("CBP", 1, "[c;R1;!$(c:n)](-!@[c;!$(c:n)])", d="biphenyl bridge C (carbocyclic rings)")
add("CQ", 2, "[c;r6](:n):n", d="aromatic C between two ring N (pyrimidine C2, triazine)")
# five-membered aromatic rings
A5 = "$(c1aaaa1)"   # in an aromatic five-membered ring
add("C5B", 1, f"[c;{A5}]", d="C in a 5-ring aromatic (C4 of imidazole, C3 of thiophene)")
add("C5A", 2, f"[c;{A5}]~[#8,#16,$([nX3])]", d="C next to O, S or pyrrole-type N in a 5-ring aromatic")
add("C5BC", 3, "[c;$(c1cocc1),$(c1ccoc1)]", d="furan C3 / C4")
add("C5BE", 3, "[c;$(c1c[nX3]cc1),$(c1cc[nX3]c1)]", d="pyrrole C3 / C4")
add("C5BD", 3, "[c;$(c1cnnc1),$(c1cnoc1),$(c1conc1)]", d="pyrazole / isoxazole C4")
add("C5BB", 3, "[c;$(c1n[nX3,o]cc1)]", d="pyrazole / isoxazole C3")
add("CRA", 4, f"[c;{A5}](~[#7,#8,#16])~[#7,#8,#16]", d="C between two heteroatoms in a 5-ring aromatic")
add("C56B", 4, f"[c;R2;{A5};$(c1aaaaa1)]", d="fused 5/6 aromatic ring C (C3a)")
add("C56A", 5, f"[c;R2;{A5};$(c1aaaaa1);$(c~[#7,#8,#16;r5])]", d="fused 5/6 aromatic ring C next to the heteroatom (C7a)")
# sp2 chains
add("CM", 0, "[CX3]", d="sp2 C")
add("CME", 1, "[CX3]=[CX3]", d="alkene C")
add("CML", 2, "[CX3](=[CX3])-[CX3]=[CX3]", d="conjugated diene C (next to the end)")
add("CM", 2, "[CX3](=[CX3])[NX3]", d="enamine C")
add("CMEE", 2, "[CX3;!$(*[OX2H1])]=C[OX2H1]", d="enol C (end)")
add("CMLE", 2, "[CX3]([OX2H1])=C", d="enol C bearing OH")
add("CD", 1, "[CX3]=N", d="imine / oxime C")
add("CHZ", 2, "[CX3]=N[NX3]", d="hydrazone C")
add("CDX", 2, "[CX3](=N)[NX3]", d="amidine C")
add("CDXR", 3, "[CX3](=N)([NX3])[NX3]", d="guanidine C")
add("CG", 3, "[CX3]=[N+]", d="iminium C")
add("CG1", 4, "[CX3](~[NX3])(~[NX3])~[NX3]", d="guanidinium C")
add("C", 1, "[CX3]=O", d="carbonyl C (amide, acid, acyl halide, urea)")
add("CO4", 2, "[CX3](=O)[OX2;!H1]", "[CX3](=O)([#6,#1])[#6,#1]", "[CX3](=O)[SX2]", d="ester / ketone / aldehyde / thioester C")
add("CO3", 3, "[CX3](=O)[OX1-]", d="carboxylate C")
add("CDS", 2, "[CX3]=S", d="thiocarbonyl C")
add("CZ", 1, "[CX2]#N", d="nitrile C")
add("CZN", 1, "[CX2]#C", d="alkyne C")
add("CZA", 1, "[CX2](=*)=*", d="cumulated C (allene, isocyanate, carbodiimide)")

# ---- nitrogen
add("NT", 0, "[NX3]", d="amine N")
add("N3", 2, "[NX3H3]", d="ammonia")
add("N3T", 2, "[NX3;r3]", d="aziridine N")
add("N4T", 2, "[NX3;r4]", d="azetidine N")
add("NE", 1, "[NX3]c", d="aniline N")
add("NEA", 2, "[NX3][CX3]=[C,N]", d="enamine / amidine / guanidine N")
add("N", 3, "[NX3][CX3]=[O,S]", d="amide / urea / carbamate N")
add("NS", 3, "[NX3][SX4](=O)=O", d="sulfonamide N")
add("NE", 4, "[NX3](c)[SX4](=O)=O", d="aryl sulfonamide N")
add("NH", 2, "[NX3][NX3]", d="hydrazine N")
add("NHZ", 3, "[NX3]N=C", "[NX3]N=O", d="hydrazone / nitrosamine amine N")
add("NX", 3, "[NX3]P", d="phosphoramide N")
add("NP", 2, "[NX4+]", d="ammonium N")
add("NB", 1, "[nX2;r6]", d="pyridine-type N")
add("NBQ", 2, "[nX2;r6;$(n:[c;R2])]", d="quinoline-type N")
add("NAP", 2, "[n+;r6]", d="pyridinium N")
add("NA", 1, "[nX3;r5]", d="pyrrole-type N")
add("N5B", 1, "[nX2;r5]", d="pyridine-type N in a 5-ring")
add("N5A", 2, "[nX2;r5;$(n~[nX3,o])]", d="pyrazole / isoxazole N2")
add("NI", 0, "[NX2]=*", d="imine / isocyanate / oxime N")
add("NN", 1, "[NX2]=[NX2]", d="azo N")
add("NO", 2, "[NX2]=O", d="nitroso N")
add("NI", 3, "[NX2](=O)O", d="nitrite N")
add("NZ", 1, "[NX1]#C", d="nitrile N")
# azide in either resonance form (R-N=N+=N- or R-N(-)-N+#N)
add("NZA", 3, "[NX2;!+]~[NX2+]~[NX1]", d="azide N (on C)")
add("NZC", 3, "[NX2+](~[NX2])~[NX1]", d="azide central N")
add("NZT", 3, "[NX1]~[NX2+]~[NX2]", d="azide terminal N")
add("NO2", 2, "[NX3+](=O)[O-]", d="nitro N")
add("NO3", 3, "[NX3+](=O)([O-])O", d="nitrate ester N")
add("NOM", 4, "[NX3+](~[OX1])(~[OX1])~[OX1]", d="nitrate ion N")
add("NG", 4, "[NX3][CX3](~[NX3])~[NX3]", d="guanidinium N")

# ---- oxygen
add("OS", 0, "[OX2]", d="ether O")
add("OAL", 1, "[OX2H1][CX4]", d="alcohol O")
add("OHP", 1, "[OX2H1]c", d="phenol O")
add("OH", 2, "[OX2H1]C=O", "[OX2H1]O", d="carboxylic acid / hydroperoxide OH")
add("OHE", 2, "[OX2H1][CX3]=C", d="enol O")
add("OHX", 2, "[OX2H1]N=C", d="oxime O")
add("OHA", 3, "[OX2H1]NC=O", d="hydroxamic acid O")
add("OES", 2, "[OX2;!H1]C=O", d="ester single-bonded O")
add("OA", 2, "[o]", d="aromatic O")
add("O3T", 3, "[OX2;r3]", d="O in a 3-ring")
add("O4T", 3, "[OX2;r4]", d="O in a 4-ring")
add("OSP", 3, "[OX2]P", d="phosphate O (single)")
add("ONA", 3, "[OX2][N+](=O)[O-]", d="nitrate ester O")
add("O", 1, "[OX1]=C", d="carbonyl O")
add("O2Z", 2, "[OX1-][CX3]=O", "[OX1]=[CX3][OX1-]", d="carboxylate O")
add("ONI", 2, "[OX1]~[#7]", d="nitro / nitroso O")
add("ON", 3, "[OX1]~[NX3+](~[OX1])~[OX1]", d="nitrate ion O")
add("O2ZP", 2, "[OX1]~P", d="phosphate / phosphine oxide terminal O")
add("OY", 2, "[OX1]~[SX4]", d="sulfone / sulfate O")
add("OZ", 2, "[OX1]~[SX3]", d="sulfoxide O")
add("OM", 1, "[OX1-;$(*[CX4]),$(*[#1])]", d="alkoxide / hydroxide O")
add("OHY", 3, "[OX3+]", d="hydronium O")
add("OT3", 5, "[OX2H2]", d="water O (TIP3P)")

# ---- sulfur, phosphorus, halogens, others
add("S", 0, "[#16X2]", d="sulfide / disulfide / thioester S")
add("SH", 1, "[SX2H1]", d="thiol S")
add("SA", 2, "[s]", d="aromatic S")
add("SD", 1, "[SX1]=*", d="thiocarbonyl S")
add("SZ", 1, "[SX3]", d="sulfoxide S")
add("SY", 1, "[SX4]", d="sulfone / sulfate S")
add("SX6", 1, "[SX6]", d="SF6")
add("PR", 1, "[PX3]", d="phosphine P")
add("P1", 0, "[PX4]", d="neutral phosphate / phosphonate P")
add("P2", 1, "[PX4;$(*(~[OX1])~[OX1]);!$(*(~[OX1])(~[OX1])~[OX1])]", d="phosphate diester anion P")
add("P3", 2, "[PX4;$(*(~[OX1])(~[OX1])~[OX1]);!$(*(~[OX1])(~[OX1])(~[OX1])~[OX1])]", d="PO3 P")
add("P4", 3, "[PX4](~[OX1])(~[OX1])(~[OX1])~[OX1]", d="phosphate ion P")
add("F", 0, "[F]", d="F on aromatic C and others")
for n in (1, 2, 3, 4):
    add(f"FX{n}", 1, f"[F][CX4,CX3;{halo(n)}]", d=f"F on aliphatic C with {n} halogen(s)")
add("FG", 2, "[F][S]", d="F on S (SF6)")
for h in ("Cl", "Br", "I"):
    add(h, 0, f"[{h}]")
for ion in ("F-", "Cl-", "Br-", "I-", "Li+", "Na+", "K+", "Rb+", "Cs+", "Mg2+", "Ca2+", "Sr2+", "Ba2+", "Ga3+"):
    sym = ion.rstrip("+-23")
    add(ion, 2, f"[{sym};X0]", d="ion")
for g in ("He", "Ne", "Ar", "Kr", "Xe", "Rn", "Si", "B"):
    add(g, 0, f"[{g}]")

out = sys.argv[1]
json.dump({"format": "caps-typing", "version": 1, "forcefield": "OPLS-AA 2005 (DL_FIELD)",
           "description": "OPLS-AA atom types under DL_FIELD's OPLS2005 key names, from the chemistry each key's DL_FIELD atom types "
                          "describe. Residue-specific keys (amino-acid side-chain CB atoms, sugars, nucleic bases) are not "
                          "assigned automatically; water is TIP3P (OT3 / HT3).",
           "rules": R}, open(os.path.join(out, "opls2005.typing.json"), "w"), indent=1)
print(len(R), "rules")
