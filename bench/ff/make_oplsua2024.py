#!/usr/bin/env python3
"""OPLS-UA from the 2024 OPLS parameter file: data/forcefields/oplsua2024.json and data/typing/oplsua2024.typing.json.

The 2024 OPLS parameter file (Jorgensen et al., J. Phys. Chem. B 2024, supporting information; as converted in
oplsaa2024-moltemplate.json) keeps the united-atom types of the original OPLS (types 66-134: Jorgensen, Madura & Swenson,
J. Am. Chem. Soc. 106, 6638 (1984) for hydrocarbons; alcohols J. Phys. Chem. 90, 1276 (1986); sulfur J. Phys. Chem. 90,
6379 (1986); ethers J. Comput. Chem. 11, 958 (1990); nitriles Mol. Phys. 63, 547 (1988); chloromethanes J. Phys. Chem.
94, 1683 (1990); and 722-723 alkyl bromide), with their bonded terms by class (C2, C3, CH, C8, C9, CD ...) and UA
torsions. This overlay extends the OPLS-AA 2024 file (the same parameters) and gives each united-atom carbon site its
mass with its hydrogens (the converted file lists 12.011 for every carbon type); its typing rules fold hydrogens on
carbon into their carbons (united_atom) and type the sites by element, hydrogen count and neighbours.

A CH3 site takes its type by the number of carbons on the carbon it is bonded to, as the 1984 paper defines it (ethane
C1, n-alkanes C2, isobutane C3, neopentane C4). Chemistry the table has no united-atom type for (a substituted aromatic
carbon, an ester, an epoxide, an ether or sulfide CH) is left untyped and reported, not guessed.

usage: make_oplsua2024.py [DATA_DIR]
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data")
BASE = os.path.join(DATA, "forcefields", "oplsaa2024-moltemplate.json")
MC, MH = 12.011, 1.008

# united-atom carbon sites of the table and the hydrogens each carries (from the table's own descriptions)
UA_H = {66: 4, 67: 3, 68: 3, 69: 3, 70: 3, 71: 2, 72: 2, 73: 1, 74: 1, 75: 1, 76: 0, 77: 0, 80: 3, 81: 2, 88: 3, 89: 2,
        90: 3, 91: 2, 92: 3, 93: 2, 96: 3, 106: 1, 107: 0, 109: 3, 110: 2, 118: 2, 120: 1, 122: 0, 723: 2}

R = []


def add(num, prio, smarts, desc):
    R.append({"type": str(num), "smarts": smarts, "priority": prio, "description": desc})


HC = "!$(*~[!#6;!#1])"   # carbon with only carbon (and hydrogen) neighbours
# hydrocarbons (1984): CH4; CH3 by the carbons on its neighbour; CH2, CH, C; alkene and aromatic CH
add(66, 0, "[CX4H4]", "CH4")
for n, d in ((67, 1), (68, 2), (69, 3), (70, 4)):
    add(n, 0, f"[CX4H3;{HC};$(*[#6D{d}])]", f"CH3 on a carbon with {d} carbon(s) ({['ethane', 'n-alkane', 'isobutane', 'neopentane'][d - 1]})")
add(71, 0, f"[CX4H2;{HC}]", "CH2 (sp3)")
add(73, 0, f"[CX4H1;{HC}]", "CH (sp3)")
add(76, 0, f"[CX4H0;{HC}]", "C (sp3)")
add(72, 1, f"[CX3H2;{HC}]=[#6]", "CH2= (sp2, 1-alkene)")
add(74, 1, f"[CX3H1;{HC}]=[#6]", "CH= (sp2, 2-alkene)")
add(77, 1, f"[CX3H0;{HC}]=[#6]", "C= (sp2, isobutene)")
add(75, 1, f"[cH1;{HC}]", "aromatic CH (benzenoid)")
# alcohols (1986): the OH keeps its hydrogen (only hydrogens on carbon fold)
add(78, 5, "[OX2H1][CX4]", "alcohol O")
add(79, 5, "[H][OX2H1][CX4]", "alcohol H(O)")
for n, h, d in ((80, 3, "CH3"), (81, 2, "CH2"), (106, 1, "CH"), (107, 0, "C")):
    add(n, 5, f"[CX4H{h}][OX2H1]", f"{d} on an alcohol O")
# ethers (1990)
add(108, 5, "[OX2H0]([CX4])[CX4]", "ether O")
add(109, 5, "[CX4H3][OX2H0][CX4]", "CH3 on an ether O")
add(110, 5, "[CX4H2][OX2H0][CX4]", "CH2 on an ether O")
# sulfur (1986): thiol, sulfide, disulfide, H2S
add(82, 5, "[SX2H2]", "S in H2S")
add(86, 5, "[H][SX2H2]", "H in H2S")
add(83, 5, "[SX2H1][CX4]", "thiol S")
add(87, 5, "[H][SX2H1][CX4]", "thiol H(S)")
add(88, 5, "[CX4H3][SX2H1]", "CH3 on a thiol S")
add(89, 5, "[CX4H2][SX2H1]", "CH2 on a thiol S")
add(84, 5, "[SX2H0]([CX4])[CX4]", "sulfide S")
add(90, 5, "[CX4H3][SX2H0;!$(*[#16])][CX4]", "CH3 on a sulfide S")
add(91, 5, "[CX4H2][SX2H0;!$(*[#16])][CX4]", "CH2 on a sulfide S")
add(85, 6, "[SX2H0]([CX4])[SX2H0]", "disulfide S")
add(92, 6, "[CX4H3][SX2H0][SX2H0]", "CH3 on a disulfide S")
add(93, 6, "[CX4H2][SX2H0][SX2H0]", "CH2 on a disulfide S")
# nitriles (1988): acetonitrile's CH3 only
add(94, 5, "[NX1]#[CX2][CX4H3]", "nitrile N")
add(95, 5, "[CX2](#[NX1])[CX4H3]", "nitrile C")
add(96, 5, "[CX4H3][CX2]#[NX1]", "CH3 on a nitrile C")
# chloromethanes (1990) and alkyl bromides
add(118, 5, "[CX4H2](Cl)Cl", "CH2 in CH2Cl2")
add(119, 5, "Cl[CX4H2]Cl", "Cl in CH2Cl2")
add(120, 5, "[CX4H1](Cl)(Cl)Cl", "CH in CHCl3")
add(121, 5, "Cl[CX4H1](Cl)Cl", "Cl in CHCl3")
add(122, 5, "[CX4H0](Cl)(Cl)(Cl)Cl", "C in CCl4")
add(123, 5, "Cl[CX4H0](Cl)(Cl)Cl", "Cl in CCl4")
add(722, 5, "[Br][CX4H2]", "alkyl Br")
add(723, 5, "[CX4H2][Br]", "CH2 on an alkyl Br")

base = json.load(open(BASE))
names = {a["name"].split("_")[0]: a for a in base["atom_types"]}
types = []
for n, h in sorted(UA_H.items()):
    a = dict(names[str(n)])
    a["mass"] = round(MC + h * MH, 4)
    types.append(a)
rules = []
for r in R:
    t = names.get(r["type"])
    if t is None:
        raise SystemExit(f"type {r['type']} is not in {BASE}")
    rules.append(dict(r, type=t["name"]))

ff = {"format": "caps-forcefield", "format_version": 1, "name": "OPLS-UA (2024 parameter file)",
      "source": "the united-atom types of the 2024 OPLS parameter file",
      "references": ["W. L. Jorgensen et al., J. Phys. Chem. B (2024), doi:10.1021/acs.jpcb.3c06602 (the parameter file)",
                     "W. L. Jorgensen, J. D. Madura, C. J. Swenson, J. Am. Chem. Soc. 106, 6638 (1984)",
                     "W. L. Jorgensen, J. Phys. Chem. 90, 1276 (1986)", "W. L. Jorgensen, J. Phys. Chem. 90, 6379 (1986)",
                     "J. M. Briggs, T. Matsui, W. L. Jorgensen, J. Comput. Chem. 11, 958 (1990)"],
      "extends": "oplsaa2024-moltemplate.json", "atom_types": types, "version": "2024", "status": "converted",
      "notes": ["the united-atom carbon sites carry their hydrogens' mass (CH4 16.043, CH3 15.035, CH2 14.027, CH 13.019)"],
      "typing": "../typing/oplsua2024.typing.json"}
json.dump(ff, open(os.path.join(DATA, "forcefields", "oplsua2024.json"), "w"), ensure_ascii=False, indent=1)
doc = {"format": "caps-typing", "version": 1, "forcefield": "OPLS-UA (2024 parameter file)",
       "description": "CAPS rules for the united-atom types of the 2024 OPLS parameter file (bench/ff/make_oplsua2024.py): hydrogens on "
                      "carbon fold into their carbons, each site typed by element, hydrogen count and neighbours. Environments the table "
                      "has no united-atom type for (substituted aromatic C, esters, epoxides, ether or sulfide CH) are left untyped and "
                      "reported. Among matching rules the highest priority wins.",
       "unknown_types": "untyped", "united_atom": True, "rules": rules}
json.dump(doc, open(os.path.join(DATA, "typing", "oplsua2024.typing.json"), "w"), ensure_ascii=False, indent=1)
print(f"{len(types)} united-atom sites, {len(rules)} rules")
