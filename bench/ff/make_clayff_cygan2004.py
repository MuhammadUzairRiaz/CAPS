#!/usr/bin/env python3
"""ClayFF for CAPS from its paper: R. T. Cygan, J.-J. Liang, A. G. Kalinichev, J. Phys. Chem. B 108, 1255 (2004).

Every number below is read from the paper (the page images, not only its text layer):
  Table 1 (p 1259)  charge, Do, Ro of each species
  Table 2 (p 1259)  k1, ro of the three O-H bonds; k2, theta_o of water's H-O-H and of every Metal-O-H bend
  Eq. 3   (p 1257)  E_VDW = Do,ij [(Ro,ij / r)^12 - 2 (Ro,ij / r)^6]
  Eqs. 4-5 (p 1257) Ro,ij = (Ro,i + Ro,j) / 2, Do,ij = sqrt(Do,i Do,j)
  Eq. 6   (p 1258)  E_bond = k1 (r - ro)^2;  Eq. 7 (p 1259) E_angle = k2 (theta - theta_o)^2
  p 1258            the O and H of a hydroxyl or a water molecule do not interact non-bonded (1-2 and 1-3 excluded);
                    hydrogens carry no Lennard-Jones term; Ewald sums for the Coulomb energy

CAPS's form: Lennard-Jones epsilon = Do, sigma = Ro / 2^(1/6), Lorentz-Berthelot mixing (the same as Eqs. 4-5, since sigma
is proportional to Ro). Metal-oxygen pairs are not bonded in ClayFF; CAPS finds the metal of a Metal-O-H bend within
2.6 A of the hydroxyl O ("angle_contacts"; the distance is CAPS's, not the paper's). The masses are standard atomic
weights (the paper lists none); the 10 A cut-off is CAPS's default.

Typing (data/typing/clayff-cygan2004.typing.json): CAPS's own rules by element, coordination and neighbours, following
the species the paper names (Table 1 and pp 1258-1259):
  st / at / ao      Si; Al with four O (tetrahedral) or more (octahedral)
  mgo, cao, feo, lio  octahedral metals; mgh, cah the metal of a hydroxide (every O around it an OH: brucite, portlandite)
  ob                bridging O (no H); obts next to a tetrahedral Al; obos next to one octahedral substitution (Mg, Fe, Li);
                    obss next to two octahedral substitutions or to a tetrahedral and an octahedral one (p 1259)
  oh / ohs          hydroxyl O; ohs next to a substitution in a silicate ("with substitution"; brucite's OH stays oh,
                    whose charges and mgh's make Mg(OH)2 neutral)
  ho, o*, h*        hydroxyl H; SPC water O and H
  Na K Cs Ca Ba Cl  the aqueous ions (Table 1, footnotes a-d); Ca is an ion unless it holds hydroxyl O (cao / cah)
Only the O-H bonds the force field defines are kept ("bonds": "defined").

usage: make_clayff_cygan2004.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
FID = "clayff-cygan2004"
REF = "R. T. Cygan, J.-J. Liang, A. G. Kalinichev, J. Phys. Chem. B 108, 1255 (2004)"
T1, T2 = "Table 1, p 1259", "Table 2, p 1259"
SIX = 2 ** (1 / 6)
# standard atomic weights (not in the paper)
MASS = {"H": 1.008, "O": 15.999, "Si": 28.085, "Al": 26.982, "Mg": 24.305, "Ca": 40.078, "Fe": 55.845, "Li": 6.94,
        "Na": 22.990, "K": 39.098, "Cs": 132.91, "Ba": 137.33, "Cl": 35.45}
# Table 1 as printed: species, symbol, element, charge (e), Do (kcal/mol), Ro (A), the source the footnote names
TABLE1 = [
    ("water hydrogen", "h*", "H", 0.4100, None, None, ""),
    ("hydroxyl hydrogen", "ho", "H", 0.4250, None, None, ""),
    ("water oxygen", "o*", "O", -0.8200, 0.1554, 3.5532, ""),
    ("hydroxyl oxygen", "oh", "O", -0.9500, 0.1554, 3.5532, ""),
    ("bridging oxygen", "ob", "O", -1.0500, 0.1554, 3.5532, ""),
    ("bridging oxygen with octahedral substitution", "obos", "O", -1.1808, 0.1554, 3.5532, ""),
    ("bridging oxygen with tetrahedral substitution", "obts", "O", -1.1688, 0.1554, 3.5532, ""),
    ("bridging oxygen with double substitution", "obss", "O", -1.2996, 0.1554, 3.5532, ""),
    ("hydroxyl oxygen with substitution", "ohs", "O", -1.0808, 0.1554, 3.5532, ""),
    ("tetrahedral silicon", "st", "Si", 2.1000, 1.8405e-6, 3.7064, ""),
    ("octahedral aluminum", "ao", "Al", 1.5750, 1.3298e-6, 4.7943, ""),
    ("tetrahedral aluminum", "at", "Al", 1.5750, 1.8405e-6, 3.7064, ""),
    ("octahedral magnesium", "mgo", "Mg", 1.3600, 9.0298e-7, 5.9090, ""),
    ("hydroxide magnesium", "mgh", "Mg", 1.0500, 9.0298e-7, 5.9090, ""),
    ("octahedral calcium", "cao", "Ca", 1.3600, 5.0298e-6, 6.2484, ""),
    ("hydroxide calcium", "cah", "Ca", 1.0500, 5.0298e-6, 6.2428, ""),   # 6.2428 as printed (cao 6.2484)
    ("octahedral iron", "feo", "Fe", 1.5750, 9.0298e-6, 5.5070, ""),
    ("octahedral lithium", "lio", "Li", 0.5250, 9.0298e-6, 4.7257, ""),
    ("aqueous sodium ion", "Na", "Na", 1.0, 0.1301, 2.6378, "footnote a: Smith and Dang, ref 32"),
    ("aqueous potassium ion", "K", "K", 1.0, 0.1000, 3.7423, "footnote b: Koneshan et al., ref 60"),
    ("aqueous cesium ion", "Cs", "Cs", 1.0, 0.1000, 4.3002, "footnote c: Smith and Dang, ref 59"),
    ("aqueous calcium ion", "Ca", "Ca", 2.0, 0.1000, 3.2237, "footnote b: Koneshan et al., ref 60"),
    ("aqueous barium ion", "Ba", "Ba", 2.0, 0.0470, 4.2840, "footnote d: Aqvist, ref 61"),
    ("aqueous chloride ion", "Cl", "Cl", -1.0, 0.1001, 4.9388, "footnote a: Smith and Dang, ref 32"),
]
# Table 2: bond stretch (i, j, k1 kcal/(mol A^2), ro A); angle bend (i, j, k, k2 kcal/(mol rad^2), theta_o deg)
BONDS = [("o*", "h*", 554.1349, 1.0000), ("oh", "ho", 554.1349, 1.0000), ("ohs", "ho", 554.1349, 1.0000)]
ANGLES = [("h*", "o*", "h*", 45.7696, 109.47), ("Metal", "oh", "ho", 30.0, 109.47), ("Metal", "ohs", "ho", 30.0, 109.47)]
METALS = ["st", "at", "ao", "mgo", "mgh", "cao", "cah", "feo", "lio"]   # the lattice metals of Table 1 ("Metal" of Table 2)

esc = lambda t: t.replace("*", "\\*")
types, pairs = [], []
for species, sym, el, q, D, R, src in TABLE1:
    cite = T1 + (f"; {src}" if src else "")
    types.append({"name": sym, "element": el, "mass": MASS[el], "charge": q, "description": f"{species} ({cite})"})
    if D is None:   # Table 1 gives the hydrogens no Do or Ro: Coulomb only (p 1258)
        pairs.append({"name": sym, "match": [esc(sym)], "params": [0.0, 0.0], "comment": f"no Lennard-Jones term ({cite}; p 1258)"})
    else:
        pairs.append({"name": sym, "match": [esc(sym)], "params": [D, R / SIX], "comment": f"Do {D:g} kcal/mol, Ro {R:.4f} A ({cite})"})
bonds = [{"name": f"{a}-{b}", "match": [esc(a), esc(b)], "style": "harmonic", "params": [k, r0],
          "comment": f"k1 {k} kcal/(mol A^2), ro {r0:.4f} A ({T2}; Eq. 6, p 1258: k1 (r - ro)^2)"} for a, b, k, r0 in BONDS]
angles = []
for i, j, k, k2, t0 in ANGLES:
    for m in (METALS if i == "Metal" else [i]):
        angles.append({"name": f"{m}-{j}-{k}", "match": [esc(m), esc(j), esc(k)], "style": "harmonic", "params": [k2, t0],
                       "comment": f"{i}-{j}-{k}: k2 {k2} kcal/(mol rad^2), theta_o {t0} deg ({T2}; Eq. 7, p 1259)"})
doc = {"format": "caps-forcefield", "format_version": 1, "name": "ClayFF (Cygan 2004)", "version": "2004 (Tables 1 and 2)",
       "source": REF + ", Tables 1 and 2 (p 1259), Eqs. 3-7 (pp 1257-1259)", "references": [REF],
       "units": "real", "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "arithmetic", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10, "angle_contacts": 2.6,
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_clayff_cygan2004.py) from the paper's Tables 1 and 2 (p 1259)",
                 "E = Do [(Ro/r)^12 - 2 (Ro/r)^6] (Eq. 3) as Lennard-Jones: epsilon = Do, sigma = Ro / 2^(1/6); Ro arithmetic, "
                 "Do geometric means (Eqs. 4-5, p 1257)",
                 "O-H of hydroxyls and water: 1-2 and 1-3 pairs excluded (p 1258); hydrogens have no Lennard-Jones term",
                 "Metal-O-H bends (Table 2) for the lattice metals st, at, ao, mgo, mgh, cao, cah, feo, lio, found within 2.6 A "
                 "of the hydroxyl O (no metal-oxygen bonds in ClayFF; the distance is CAPS's)",
                 "Table 1 prints Ro 6.2484 A for cao and 6.2428 A for cah: both kept as printed",
                 "masses: standard atomic weights (not in the paper); cut-off 10 A: CAPS's default"],
       "typing": f"../typing/{FID}.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1, ensure_ascii=False)

AL_TET, SUB, OH = "[Al;X4]", "[Mg,Fe,Li]", "[O;H1]"
R = []
def rule(t, p, s, d, **kw):
    R.append({"type": t, "smarts": s, "priority": p, "description": d, **kw})
rule("st", 0, "[Si]", "tetrahedral silicon")
rule("ao", 0, "[Al]", "octahedral aluminium"); rule("at", 1, AL_TET, "tetrahedral aluminium (four O)")
rule("mgo", 0, "[Mg;!X0]", "octahedral magnesium")
rule("mgh", 1, "[Mg;!X0;!$(*~[O;H0])]", "magnesium of a hydroxide (every O around it an OH)")
rule("Ca", 0, "[Ca]", "aqueous calcium ion (Ca holding no hydroxyl O: interlayer, in water)")
rule("cao", 1, f"[Ca;$(*~{OH})]", "octahedral calcium (holds a hydroxyl O)")
rule("cah", 2, f"[Ca;$(*~{OH});!$(*~[O;H0])]", "calcium of a hydroxide (every O around it an OH)")
rule("feo", 0, "[Fe]", "octahedral iron"); rule("lio", 0, "[Li]", "octahedral lithium")
for ion in ("Na", "K", "Cs", "Ba", "Cl"):   # whatever their contacts: ClayFF has no other species of these elements
    rule(ion, 0, f"[{ion}]", f"aqueous {ion} ion")
rule("ob", 0, "[O;H0]", "bridging oxygen")
rule("obts", 2, f"[O;H0;$(*~{AL_TET})]", "bridging oxygen next to a tetrahedral substitution (Al for Si)")
rule("obos", 2, f"[O;H0;$(*~{SUB});$(*~[Si,Al])]", "bridging oxygen next to an octahedral substitution")
rule("obss", 3, f"[O;H0;$(*(~{SUB})~{SUB});$(*~[Si,Al])]", "bridging oxygen next to two octahedral substitutions")
rule("obss", 3, f"[O;H0;$(*~{AL_TET});$(*~{SUB})]", "bridging oxygen next to a tetrahedral and an octahedral substitution")
rule("oh", 1, OH, "hydroxyl oxygen")
rule("ohs", 2, f"[O;H1;$(*~{SUB})]", "hydroxyl oxygen next to an octahedral substitution, in a silicate", requires=["Si"])
rule("ho", 1, f"[H]{OH}", "hydroxyl hydrogen")
rule("o*", 5, "[O;H2]", "SPC water oxygen"); rule("h*", 5, "[H][O;H2]", "SPC water hydrogen")
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined", "rules": R,
           "description": "CAPS rules for the ClayFF species of Cygan, Liang & Kalinichev 2004 (Table 1) by element, coordination and "
                          "neighbours (bench/ff/make_clayff_cygan2004.py); only the O-H bonds the force field defines are kept"},
          open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1, ensure_ascii=False)
print(f"{FID}: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(R)} rules")
