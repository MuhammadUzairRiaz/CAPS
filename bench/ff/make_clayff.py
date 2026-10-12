#!/usr/bin/env python3
"""ClayFF for CAPS from its paper: R. T. Cygan, J.-J. Liang, A. G. Kalinichev, J. Phys. Chem. B 108, 1255 (2004),
Tables 1 (nonbond: charge, Do, Ro per species) and 2 (bond stretch, angle bend).

  E_VDW = Do [(Ro/r)^12 - 2 (Ro/r)^6]: Lennard-Jones epsilon = Do, sigma = Ro / 2^(1/6); Ro by the arithmetic mean,
          Do by the geometric mean (the paper's combination rules)
  bonds   k1 (r - ro)^2 for the O-H of water and hydroxyls; angles k2 (theta - thetao)^2 for water's H-O-H and every
          metal-O-H (metal-oxygen pairs are not bonded in ClayFF: CAPS finds the metal within 2.6 A of the hydroxyl O)
  1-2 and 1-3 pairs of the bonded O-H groups excluded; everything else non-bonded

usage: make_clayff.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "R. T. Cygan, J.-J. Liang, A. G. Kalinichev, J. Phys. Chem. B 108, 1255 (2004)"
SIX = 2 ** (1 / 6)
# Table 1: symbol, description, element, mass, charge, Do (kcal/mol), Ro (A), the source the paper names for the ion
T1 = [
    ("h*", "water hydrogen", "H", 1.00797, 0.4100, 0.0, 0.0, ""),
    ("ho", "hydroxyl hydrogen", "H", 1.00797, 0.4250, 0.0, 0.0, ""),
    ("o*", "water oxygen", "O", 15.9994, -0.8200, 0.1554, 3.5532, ""),
    ("oh", "hydroxyl oxygen", "O", 15.9994, -0.9500, 0.1554, 3.5532, ""),
    ("ob", "bridging oxygen", "O", 15.9994, -1.0500, 0.1554, 3.5532, ""),
    ("obos", "bridging oxygen with octahedral substitution", "O", 15.9994, -1.1808, 0.1554, 3.5532, ""),
    ("obts", "bridging oxygen with tetrahedral substitution", "O", 15.9994, -1.1688, 0.1554, 3.5532, ""),
    ("obss", "bridging oxygen with double substitution", "O", 15.9994, -1.2996, 0.1554, 3.5532, ""),
    ("ohs", "hydroxyl oxygen with substitution", "O", 15.9994, -1.0808, 0.1554, 3.5532, ""),
    ("st", "tetrahedral silicon", "Si", 28.0855, 2.1000, 1.8405e-6, 3.7064, ""),
    ("ao", "octahedral aluminum", "Al", 26.98154, 1.5750, 1.3298e-6, 4.7943, ""),
    ("at", "tetrahedral aluminum", "Al", 26.98154, 1.5750, 1.8405e-6, 3.7064, ""),
    ("mgo", "octahedral magnesium", "Mg", 24.305, 1.3600, 9.0298e-7, 5.9090, ""),
    ("mgh", "hydroxide magnesium", "Mg", 24.305, 1.0500, 9.0298e-7, 5.9090, ""),
    ("cao", "octahedral calcium", "Ca", 40.078, 1.3600, 5.0298e-6, 6.2484, ""),
    ("cah", "hydroxide calcium", "Ca", 40.078, 1.0500, 5.0298e-6, 6.2428, ""),
    ("feo", "octahedral iron", "Fe", 55.847, 1.5750, 9.0298e-6, 5.5070, ""),
    ("lio", "octahedral lithium", "Li", 6.941, 0.5250, 9.0298e-6, 4.7257, ""),
    ("Na", "aqueous sodium ion", "Na", 22.98977, 1.0, 0.1301, 2.6378, "Smith and Dang"),
    ("K", "aqueous potassium ion", "K", 39.0983, 1.0, 0.1000, 3.7423, "Koneshan et al."),
    ("Cs", "aqueous cesium ion", "Cs", 132.9054, 1.0, 0.1000, 4.3002, "Smith and Dang"),
    ("Ca", "aqueous calcium ion", "Ca", 40.078, 2.0, 0.1000, 3.2237, "Koneshan et al."),
    ("Ba", "aqueous barium ion", "Ba", 137.327, 2.0, 0.0470, 4.2840, "Aqvist"),
    ("Cl", "aqueous chloride ion", "Cl", 35.453, -1.0, 0.1001, 4.9388, "Smith and Dang"),
]
METALS = ["st", "ao", "at", "mgo", "mgh", "cao", "cah", "feo", "lio"]
esc = lambda t: t.replace("*", "\\*")
types = [{"name": s, "element": e, "mass": m, "charge": q, "description": f"{d} (Table 1" + (f"; {src}" if src else "") + ")"}
         for s, d, e, m, q, D, R, src in T1]
pairs = [{"name": s, "match": [esc(s)], "params": [D, R / SIX], "comment": f"Do {D:g} kcal/mol, Ro {R:g} A"} for s, d, e, m, q, D, R, src in T1]
bonds = [{"name": f"{a}-{b}", "match": [esc(a), esc(b)], "style": "harmonic", "params": [554.1349, 1.0], "comment": "Table 2"}
         for a, b in (("o*", "h*"), ("oh", "ho"), ("ohs", "ho"))]
angles = [{"name": "h*-o*-h*", "match": ["h\\*", "o\\*", "h\\*"], "style": "harmonic", "params": [45.7696, 109.47], "comment": "Table 2: water"}]
angles += [{"name": f"{m}-{o}-ho", "match": [m, o, "ho"], "style": "harmonic", "params": [30.0, 109.47], "comment": "Table 2: metal-O-H"}
           for o in ("oh", "ohs") for m in METALS]
doc = {"format": "caps-forcefield", "format_version": 1, "name": "ClayFF (Cygan 2004)", "source": REF + ", Tables 1 and 2", "references": [REF],
       "units": "real", "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "arithmetic", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10, "angle_contacts": 2.6,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_clayff.py) from the paper's Tables 1 and 2",
                 "metal-O-H bends: the metal within 2.6 A of the hydroxyl oxygen (ClayFF bonds no metal; the distance is CAPS's)"],
       "version": "2004", "typing": "../typing/clayff2004.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", "clayff2004.json"), "w"), indent=1, ensure_ascii=False)

AL_OCT, AL_TET, SUB = "[Al;!X4]", "[Al;X4]", "[Mg,Fe,Li]"
R = []
def rule(t, p, s, d):
    R.append({"type": t, "smarts": s, "priority": p, "description": d})
rule("st", 0, "[Si]", "tetrahedral silicon")
rule("ao", 0, "[Al]", "octahedral aluminium"); rule("at", 1, AL_TET, "tetrahedral aluminium")
rule("mgo", 0, "[Mg;!X0]", "octahedral magnesium"); rule("mgh", 1, "[Mg;!X0;!$(*~[O;H0])]", "magnesium of a hydroxide (every O an OH)")
rule("cao", 0, "[Ca;!X0]", "octahedral calcium"); rule("cah", 1, "[Ca;!X0;!$(*~[O;H0])]", "calcium of a hydroxide")
rule("feo", 0, "[Fe]", "octahedral iron"); rule("lio", 0, "[Li]", "octahedral lithium")
for ion, s in (("Na", "[Na;X0]"), ("K", "[K;X0]"), ("Cs", "[Cs;X0]"), ("Ca", "[Ca;X0]"), ("Ba", "[Ba;X0]"), ("Cl", "[Cl;X0]")):
    rule(ion, 2, s, f"aqueous {ion} ion")
rule("ob", 0, "[O;H0]", "bridging oxygen")
rule("obts", 2, f"[O;H0;$(*~{AL_TET})]", "bridging oxygen with tetrahedral substitution (Al for Si)")
rule("obos", 2, f"[O;H0;$(*~{SUB});$(*~[Si,Al])]", "bridging oxygen with octahedral substitution")
rule("obss", 3, f"[O;H0;$(*~{AL_TET});$(*~{SUB})]", "bridging oxygen with double substitution")
rule("oh", 1, "[O;H1]", "hydroxyl oxygen"); rule("ohs", 2, f"[O;H1;$(*~{SUB});$(*~[Al,Si])]", "hydroxyl oxygen with substitution")
rule("ho", 1, "[H][O;H1]", "hydroxyl hydrogen")
rule("o*", 5, "[OX2H2]", "SPC water oxygen"); rule("h*", 5, "[H][OX2H2]", "SPC water hydrogen")
json.dump({"format": "caps-typing", "version": 1, "forcefield": "ClayFF (Cygan 2004)", "unknown_types": "untyped", "bonds": "defined", "rules": R,
           "description": "CAPS rules for the paper's ClayFF species by coordination and neighbours (bench/ff/make_clayff.py); only the O-H "
                          "bonds the force field defines are kept"},
          open(os.path.join(DATA, "typing", "clayff2004.typing.json"), "w"), indent=1, ensure_ascii=False)
print(f"clayff2004: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(R)} rules")
