#!/usr/bin/env python3
"""Water models for CAPS, each its own library, every number from the model's paper or its authors' own distribution.

Each library has the types OW, HW (and MW for a four-site model, LP for TIP5P's two lone pairs), the model's charges, the
O-O Lennard-Jones (and CHARMM's hydrogen Lennard-Jones for its TIP3P), the O-H bond and H-O-H angle at the model's
geometry, and the massless sites as virtual sites: M = (1 - a) O + a/2 (H1 + H2) with a = d_OM / (r_OH cos(theta/2)) (LAMMPS's
tip4p styles and GROMACS's virtual_sites3 place it there), TIP5P's lone pairs O + b (r_OH1 + r_OH2) +- c (r_OH1 x r_OH2) with
b = -d_OL cos(phi/2) / (2 r_OH cos(theta/2)) and c = d_OL sin(phi/2) / (r_OH^2 sin(theta)) (GROMACS's 3out sites). mW (Molinero &
Moore 2009) is one Stillinger-Weber site per molecule.

The models are rigid (SHAKE / RATTLE in LAMMPS, SETTLE in GROMACS): their bond and angle constants serve a flexible run only
and are not the model's; each library names where they come from. Where a model's own paper could not be obtained the
numbers come from a paper by the same authors that tabulates them, and the library says so (SOURCES below).

Sources (PDFs and text in ~/project/ff-sources/papers/water, never in the repository):
  TIP3P, TIP4P   M. W. Mahoney, W. L. Jorgensen, J. Chem. Phys. 112, 8910 (2000), Table I (the TIP3P and TIP4P columns
                 restate Jorgensen et al. 1983 by its first author)
  TIP5P          the same paper, Table I
  TIP4P/2005     C. Vega, J. L. F. Abascal, M. M. Conde, J. L. Aragones, Faraday Discuss. 141, 251 (2009), Table I (the
                 model's authors; arXiv:0901.1803)
  OPC            S. Izadi, R. Anandakrishnan, A. V. Onufriev, J. Phys. Chem. Lett. 5, 3863 (2014), Table 2 (PMC4226301)
  OPC3           S. Izadi, A. V. Onufriev, J. Chem. Phys. 145, 074501 (2016), Table II (PMC4991989)
  SPC, SPC/E     the GROMOS distribution of the models' authors as GROMACS ships it (share/top/gromos43a1.ff: spc.itp,
                 spce.itp, ffnonbonded.itp OW C6 / C12; LGPL 2.1): charges, geometry (settles), Lennard-Jones and the
                 FLEXIBLE bond and angle constants
  TIP3P (CHARMM) CHARMM's toppar_water_ions.str (TIP3P with hydrogen Lennard-Jones, HT-OT 450 / 0.9572, HT-OT-HT 55 /
                 104.52), the model of MacKerell et al. 1998; checked against GROMACS's charmm27.ff (LGPL 2.1)
  mW             V. Molinero, E. B. Moore, J. Phys. Chem. B 113, 4008 (2009) (arXiv:0809.2811): lambda 23.15, eps 6.189
                 kcal/mol, sigma 2.3925 A; A 7.049556277, B 0.6022245584, p 4, q 0, gamma 1.2, a 1.8, cos theta0 -1/3

usage: make_water_models.py [DATA_DIR]
"""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
KJ = 4.184                      # kJ per kcal
KB = 0.0019872043               # kcal/mol/K (CODATA 2018 k_B N_A)
M_O, M_H = 15.9994, 1.008
RMIN = 2 ** (1 / 6)

JORGENSEN1983 = "W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, M. L. Klein, J. Chem. Phys. 79, 926 (1983)"
MAHONEY2000 = "M. W. Mahoney, W. L. Jorgensen, J. Chem. Phys. 112, 8910 (2000)"
VEGA2009 = "C. Vega, J. L. F. Abascal, M. M. Conde, J. L. Aragones, Faraday Discuss. 141, 251 (2009)"
GROMOS = "GROMACS share/top/gromos43a1.ff (spc.itp, spce.itp, ffnonbonded.itp; GNU LGPL 2.1)"
CHARMM_STR = "CHARMM toppar_water_ions.str (toppar_c36_jul17)"
# CHARMM's TIP3P bond and angle (HT-OT 450.0 kcal/mol/A^2 at 0.9572 A, HT-OT-HT 55.0 kcal/mol/rad^2 at 104.52 deg), LAMMPS's K
TIP_K = (450.0, 55.0, "CHARMM's TIP3P constants (toppar_water_ions.str: HT-OT 450.0, HT-OT-HT 55.0) at the model's geometry; "
                      "the model is rigid, they serve a flexible run only")
# GROMOS's flexible SPC: k_b 345000 kJ/mol/nm^2, k_theta 383 kJ/mol/rad^2 for k/2 (x - x0)^2, as LAMMPS's K = k/2
SPC_K = (345000 / KJ / 100 / 2, 383 / KJ / 2, "GROMOS's FLEXIBLE SPC constants (spc.itp: k_b 345000 kJ/mol/nm^2, k_theta 383 "
                                              "kJ/mol/rad^2) at the model's geometry; the model is rigid, they serve a flexible run only")


def c6c12(c6, c12):
    """sigma (A), epsilon (kcal/mol) of GROMACS C6 (kJ nm^6/mol), C12 (kJ nm^12/mol)"""
    return 10 * (c12 / c6) ** (1 / 6), c6 * c6 / (4 * c12) / KJ


SPC_SIG, SPC_EPS = c6c12(0.0026173456, 2.634129e-06)

# id: the model. q_h, r_oh (A), theta (deg), sigma (A), eps (kcal/mol), d_om (A) for 4 sites, d_ol / theta_l for 5 sites,
# eps_h / sigma_h (hydrogen LJ), k (bond K, angle K, where they come from), mixing, refs, source (what the numbers were read
# from), year, name
MODELS = {
    "spc-berendsen1981": dict(
        name="SPC water", short="SPC", q_h=0.41, r_oh=1.0, theta=109.47, sigma=SPC_SIG, eps=SPC_EPS, k=SPC_K, mixing="geometric", year="1981",
        refs=["H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, J. Hermans, in Intermolecular Forces, ed. B. Pullman "
              "(Reidel, Dordrecht, 1981), p. 331"],
        source=GROMOS + ": spc.itp (OW -0.82, H 0.41; settles d_OH 0.1 nm, d_HH 0.16330 nm), OW C6 0.0026173456, C12 2.634129e-06",
        values="q_H 0.41 e, r_OH 1.0 A, H-O-H 109.47 deg, O-O C6 0.0026173456 kJ nm^6/mol, C12 2.634129e-06 kJ nm^12/mol "
               "(sigma 3.165648 A, epsilon 0.650167 kJ/mol)"),
    "spce-berendsen1987": dict(
        name="SPC/E water", short="SPC/E", q_h=0.4238, r_oh=1.0, theta=109.47, sigma=SPC_SIG, eps=SPC_EPS, k=SPC_K, mixing="geometric", year="1987",
        refs=["H. J. C. Berendsen, J. R. Grigera, T. P. Straatsma, J. Phys. Chem. 91, 6269 (1987)"],
        source=GROMOS + ": spce.itp (OW -0.8476, H 0.4238; SPC's geometry), OW C6 0.0026173456, C12 2.634129e-06; "
               "the same charges in Izadi & Onufriev 2016, Table II",
        values="q_H 0.4238 e, r_OH 1.0 A, H-O-H 109.47 deg, SPC's Lennard-Jones"),
    "tip3p-jorgensen1983": dict(
        name="TIP3P water", short="TIP3P", q_h=0.417, r_oh=0.9572, theta=104.52, sigma=3.15061, eps=0.1521, k=TIP_K, mixing="geometric", year="1983",
        refs=[JORGENSEN1983], source=MAHONEY2000 + ", Table I (TIP3P column)",
        values="q_H 0.417 e, sigma 3.15061 A, epsilon 0.1521 kcal/mol, r_OH 0.9572 A, H-O-H 104.52 deg"),
    "tip3p-charmm-mackerell1998": dict(
        name="TIP3P water (CHARMM, hydrogen Lennard-Jones)", short="TIP3P (CHARMM)", q_h=0.417, r_oh=0.9572, theta=104.52,
        sigma=2 * 1.7682 / RMIN, eps=0.1521, eps_h=0.046, sigma_h=2 * 0.2245 / RMIN, k=(450.0, 55.0, "CHARMM's own: HT-OT 450.0 kcal/mol/A^2, "
        "HT-OT-HT 55.0 kcal/mol/rad^2 (toppar_water_ions.str); the model is rigid in CHARMM (SHAKE)"), mixing="arithmetic", year="1998",
        refs=["A. D. MacKerell Jr. et al., J. Phys. Chem. B 102, 3586 (1998)", JORGENSEN1983],
        source=CHARMM_STR + ": RESI TIP3 (OT -0.834, HT 0.417), OT epsilon -0.1521 Rmin/2 1.7682, HT epsilon -0.046 Rmin/2 0.2245; "
               "the same numbers in GROMACS share/top/charmm27.ff (LGPL 2.1)",
        values="q_H 0.417 e, O epsilon 0.1521 kcal/mol Rmin/2 1.7682 A, H epsilon 0.046 kcal/mol Rmin/2 0.2245 A, r_OH 0.9572 A, H-O-H 104.52 deg"),
    "tip4p-jorgensen1983": dict(
        name="TIP4P water", short="TIP4P", q_h=0.520, r_oh=0.9572, theta=104.52, sigma=3.15365, eps=0.1550, d_om=0.15, k=TIP_K, mixing="geometric",
        year="1983", refs=[JORGENSEN1983], source=MAHONEY2000 + ", Table I (TIP4P column)",
        values="q_H 0.520 e, sigma 3.15365 A, epsilon 0.1550 kcal/mol, r_OH 0.9572 A, H-O-H 104.52 deg, r_OM 0.15 A"),
    "tip4p2005-abascal2005": dict(
        name="TIP4P/2005 water", short="TIP4P/2005", q_h=0.5564, r_oh=0.9572, theta=104.52, sigma=3.1589, eps=93.2 * KB, d_om=0.1546, k=TIP_K,
        mixing="arithmetic", year="2005", refs=["J. L. F. Abascal, C. Vega, J. Chem. Phys. 123, 234505 (2005)"],
        source=VEGA2009 + ", Table I (the model's authors; arXiv:0901.1803): the original paper is not yet in hand",
        values="q_H 0.5564 e, sigma 3.1589 A, epsilon/k_B 93.2 K, r_OH 0.9572 A, H-O-H 104.52 deg, d_OM 0.1546 A"),
    "tip5p-mahoney2000": dict(
        name="TIP5P water", short="TIP5P", q_h=0.241, r_oh=0.9572, theta=104.52, sigma=3.12, eps=0.16, d_ol=0.70, theta_l=109.47, k=TIP_K,
        mixing="geometric", year="2000", refs=[MAHONEY2000], source=MAHONEY2000 + ", Table I (TIP5P column)",
        values="q_H 0.241 e (each lone pair -0.241 e), sigma 3.12 A, epsilon 0.16 kcal/mol, r_OH 0.9572 A, H-O-H 104.52 deg, r_OL 0.70 A, "
               "L-O-L 109.47 deg"),
    "opc-izadi2014": dict(
        name="OPC water", short="OPC", q_h=0.6791, r_oh=0.8724, theta=103.6, sigma=3.16655, eps=0.89036 / KJ, d_om=0.1594, k=TIP_K,
        mixing="arithmetic", year="2014", refs=["S. Izadi, R. Anandakrishnan, A. V. Onufriev, J. Phys. Chem. Lett. 5, 3863 (2014)"],
        source="the paper's Table 2 (PMC4226301)",
        values="q 0.6791 e, l 0.8724 A, z1 0.1594 A, theta 103.6 deg, sigma 3.16655 A, epsilon 0.89036 kJ/mol"),
    "opc3-izadi2016": dict(
        name="OPC3 water", short="OPC3", q_h=0.447585, r_oh=0.97888, theta=109.47, sigma=3.17427, eps=0.68369 / KJ, k=SPC_K,
        mixing="arithmetic", year="2016", refs=["S. Izadi, A. V. Onufriev, J. Chem. Phys. 145, 074501 (2016)"],
        source="the paper's Table II (PMC4991989)",
        values="q 0.447585 e, l 0.97888 A, theta 109.47 deg, sigma 3.17427 A, epsilon 0.68369 kJ/mol"),
}

MW = dict(name="mW water (monatomic, Stillinger-Weber)", year="2009",
          refs=["V. Molinero, E. B. Moore, J. Phys. Chem. B 113, 4008 (2009)"],
          source="the paper (arXiv:0809.2811, the accepted manuscript), section II",
          eps=6.189, sigma=2.3925, a=1.8, lam=23.15, gamma=1.2, cos0=-1 / 3, A=7.049556277, B=0.6022245584, p=4.0, q=0.0)

LICENSE = "published parameters of the paper"


def sites(m):
    return 5 if "d_ol" in m else 4 if "d_om" in m else 3


def m_alpha(m):
    return m["d_om"] / (m["r_oh"] * math.cos(math.radians(m["theta"]) / 2))


def lp_coefficients(m):
    ht, hl = math.radians(m["theta"]) / 2, math.radians(m["theta_l"]) / 2
    return (-m["d_ol"] * math.cos(hl) / (2 * m["r_oh"] * math.cos(ht)),
            m["d_ol"] * math.sin(hl) / (m["r_oh"] ** 2 * math.sin(2 * ht)))


def water_library(fid, m):
    n = sites(m)
    q_neg = -2 * m["q_h"]
    q_o = q_neg if n == 3 else 0.0
    kb, ka, ksrc = m["k"]
    types = [{"name": "OW", "element": "O", "mass": M_O, "charge": q_o, "description": f"{m['short']} oxygen"},
             {"name": "HW", "element": "H", "mass": M_H, "charge": m["q_h"], "description": f"{m['short']} hydrogen"}]
    pairs = [{"name": "OW", "match": ["OW"], "params": [m["eps"], m["sigma"]], "comment": "O-O Lennard-Jones: " + m["values"]},
             {"name": "HW", "match": ["HW"], "params": [m.get("eps_h", 0.0), m.get("sigma_h", 0.0)],
              "comment": "hydrogen Lennard-Jones (CHARMM's TIP3P)" if m.get("eps_h") else "no Lennard-Jones on H"}]
    vsites, notes = [], []
    if n == 4:
        a = m_alpha(m)
        types.append({"name": "MW", "element": "", "mass": 0.0, "charge": q_neg, "description": f"{m['short']} M site (massless, {m['d_om']} A from O on the bisector)"})
        pairs.append({"name": "MW", "match": ["MW"], "params": [0.0, 0.0], "comment": "charge only"})
        vsites.append({"smarts": "[OX2H2]([H])[H]", "type": "MW", "from": [0, 1, 2], "w": [1 - a, a / 2, a / 2],
                       "comment": f"M = (1 - a) O + a/2 (H1 + H2), a = d_OM / (r_OH cos(theta/2)) = {a:.9f}"})
    if n == 5:
        b, c = lp_coefficients(m)
        types.append({"name": "LP", "element": "", "mass": 0.0, "charge": -m["q_h"],
                      "description": f"{m['short']} lone pair (massless, {m['d_ol']} A from O, L-O-L {m['theta_l']} deg)"})
        pairs.append({"name": "LP", "match": ["LP"], "params": [0.0, 0.0], "comment": "charge only"})
        for sgn in (1, -1):
            vsites.append({"smarts": "[OX2H2]([H])[H]", "type": "LP", "from": [0, 1, 2], "w": [1 - 2 * b, b, b], "c": sgn * c,
                           "comment": "L = O + b (r_OH1 + r_OH2) %s c (r_OH1 x r_OH2), b = %.9f, c = %.9f 1/A" % ("+" if sgn > 0 else "-", b, c)})
        notes.append("five sites: GROMACS writes the lone pairs as virtual_sites3 (3out); LAMMPS has no style for them")
    doc = {"format": "caps-forcefield", "format_version": 1, "name": m["name"], "source": m["source"], "references": m["refs"],
           "units": "real",
           "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"},
           "mixing": m["mixing"], "special_lj": [0, 0, 0], "special_coul": [0, 0, 0], "cutoff": 10,
           "atom_types": types, "pairs": pairs,
           "bonds": [{"name": "OW-HW", "match": ["OW", "HW"], "style": "harmonic", "params": [kb, m["r_oh"]], "comment": ksrc}],
           "angles": [{"name": "HW-OW-HW", "match": ["HW", "OW", "HW"], "style": "harmonic", "params": [ka, m["theta"]], "comment": ksrc}],
           "dihedrals": [], "impropers": [],
           "notes": [f"built by CAPS (bench/ff/make_water_models.py) from {m['source']}",
                     "rigid model: run it with its bonds and angle constrained (SHAKE / RATTLE in LAMMPS, SETTLE in GROMACS)"] + notes,
           "version": m["year"], "status": "validated", "typing": f"../typing/{fid}.typing.json"}
    if vsites:
        doc["virtual_sites"] = vsites
    rules = [{"type": "OW", "smarts": "[OX2H2]", "priority": 0, "description": "water oxygen"},
             {"type": "HW", "smarts": "[H][OX2H2]", "priority": 0, "description": "water hydrogen"}]
    if n == 4:
        rules.append({"type": "MW", "smarts": "[#0]", "atom_name": "MW", "priority": 0, "description": "the M site (added by CAPS on the bisector)"})
    if n == 5:
        rules.append({"type": "LP", "smarts": "[#0]", "atom_name": "LP", "priority": 0, "description": "a lone pair (added by CAPS)"})
    typing = {"format": "caps-typing", "version": 1, "forcefield": m["name"], "unknown_types": "untyped",
              "description": f"CAPS rules for {m['short']} water (bench/ff/make_water_models.py): water's O and H; "
                             + ("the massless sites are added by the force field" if n > 3 else "three sites"),
              "rules": rules}
    return doc, typing


def mw_library():
    m = MW
    doc = {"format": "caps-forcefield", "format_version": 1, "name": m["name"], "source": m["source"], "references": m["refs"],
           "units": "real", "styles": {"pair": "sw", "bond": "harmonic", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"},
           "mixing": "arithmetic", "special_lj": [0, 0, 0], "special_coul": [0, 0, 0], "cutoff": 10,
           "atom_types": [{"name": "MW", "element": "O", "mass": M_O + 2 * M_H, "charge": 0.0,
                           "description": "mW water: one site per molecule at its oxygen, the molecule's mass"}],
           "pairs": [{"name": "MW", "match": ["MW"], "style": "sw",
                      "params": [m["eps"], m["sigma"], m["a"], m["lam"], m["gamma"], m["cos0"], m["A"], m["B"], m["p"], m["q"], 0.0],
                      "comment": "Stillinger-Weber (LAMMPS pair sw: epsilon sigma a lambda gamma cos(theta0) A B p q tol): epsilon 6.189 "
                                 "kcal/mol, sigma 2.3925 A, lambda 23.15 (the paper's section II); a 1.8, gamma 1.2, cos(theta0) -1/3, "
                                 "A 7.049556277, B 0.6022245584, p 4, q 0 (Stillinger-Weber's, as the paper gives them)"}],
           "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
           "notes": ["built by CAPS (bench/ff/make_water_models.py) from the paper",
                     "each water's hydrogens fold into its oxygen: one uncharged site per molecule, no bonds"],
           "version": m["year"], "status": "validated", "typing": "../typing/mw-molinero2009.typing.json"}
    typing = {"format": "caps-typing", "version": 1, "forcefield": m["name"], "united_atom": True, "united_atom_hosts": ["O"],
              "unknown_types": "untyped",
              "description": "CAPS rules for mW water (bench/ff/make_water_models.py): each water becomes one site at its oxygen",
              "rules": [{"type": "MW", "smarts": "[O;H2;X2]", "priority": 0, "description": "a water molecule (its oxygen, both hydrogens folded in)"},
                        {"type": "MW", "smarts": "[O;X0]", "atom_name": "MW", "priority": 0, "description": "an mW site read back (a site named MW)"}]}
    return doc, typing


def catalogue_entry(fid, doc, short, evidence, notes):
    return {"id": fid, "name": doc["name"], "version": doc["version"], "references": doc["references"], "origin": "built from source",
            "source": doc["source"], "license": LICENSE, "status": "validated", "notes": notes, "file": fid + ".json",
            "typing": {"rules": f"typing/{fid}.typing.json", "evidence": evidence},
            "counts": {"atom_types": len(doc["atom_types"]), **{k: len(doc.get(k, [])) for k in ("pairs", "bonds", "angles", "dihedrals", "impropers")}},
            "group": "water"}


def main():
    entries = []
    for fid, m in MODELS.items():
        doc, typing = water_library(fid, m)
        json.dump(doc, open(os.path.join(DATA, "forcefields", fid + ".json"), "w"), indent=1, ensure_ascii=False)
        json.dump(typing, open(os.path.join(DATA, "typing", fid + ".typing.json"), "w"), indent=1, ensure_ascii=False)
        n = sites(m)
        entries.append(catalogue_entry(fid, doc, m["short"], "water's O and H by SMARTS; massless sites added and placed by the force field",
                                       f"{m['short']}: {n}-site rigid water; geometry and charges reproduced exactly, LAMMPS and GROMACS "
                                       "energies equal CAPS's (bench/ff/check_water_libraries.py)"))
    doc, typing = mw_library()
    json.dump(doc, open(os.path.join(DATA, "forcefields", "mw-molinero2009.json"), "w"), indent=1, ensure_ascii=False)
    json.dump(typing, open(os.path.join(DATA, "typing", "mw-molinero2009.typing.json"), "w"), indent=1, ensure_ascii=False)
    entries.append(catalogue_entry("mw-molinero2009", doc, "mW", "each water one site at its oxygen (united atom, hydrogens folded in)",
                                   "mW: one Stillinger-Weber site per water; LAMMPS's pair sw equals CAPS's energy (bench/ff/check_water_libraries.py)"))
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    ids = {e["id"] for e in entries}
    cat["forcefields"] = [e for e in cat["forcefields"] if e["id"] not in ids] + entries
    json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)
    print(f"{len(entries)} water libraries written")


if __name__ == "__main__":
    main()
