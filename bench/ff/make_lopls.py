#!/usr/bin/env python3
"""L-OPLS for CAPS: the long-hydrocarbon refinement of OPLS-AA by S. W. I. Siu, K. Pluhackova, R. A. Boeckmann, J. Chem.
Theory Comput. 8, 1459 (2012), Table 2, as an overlay on CAPS's OPLS-AA 2024 (data/forcefields/oplsaa2024.json, "extends").
The paper's Supporting Information holds only the ab initio torsion scans (hexane, cis-3- and cis-2-hexene); every
parameter comes from Table 2 of the main text (PDF in ~/project/ff-sources/papers/lopls, never in the repository).

  torsions (Table 2, the Fourier columns, kJ/mol; GROMACS's Fourier form is LAMMPS's opls K1..K4 / 4.184; the table's V0
  is the constant the Ryckaert-Bellemans fit carries and has no force)
    CT-CT-CT-CT   2.697394  -0.896807  0.74567    (hexane)
    X-CM-CM-X     0          51.2551    0          (cis-3-hexene)
    CM-CM-CT-CT  -3.368171   1.34679   -0.4321105  (cis-2-hexene)
    CM-CT-CT-CT   2.017484   0.562197   0.74369    (cis-2-hexene)
  non-bonded (charge e, sigma nm, epsilon kJ/mol)
    CT in CH3  -0.222  0.35   0.276144      HC on CH3  0.074  0.25   0.125520
    CT in CH2  -0.148  0.35   0.276144      HC on CH2  0.074  0.25   0.110000
    CM in CH=  -0.160  0.355  0.317984      HC on CH=  0.160  0.242  0.125520

The paper's atom types (Table 1: CT = opls_135 / opls_136, CM = opls_142, HC = opls_140 / opls_144) become six L-OPLS
types (lopls_*) for the carbons and hydrogens of hydrocarbon chains: a CH3, CH2 or alkene CH whose neighbours are
hydrogens and carbons that carry no heteroatom and are not aromatic; every such group stays neutral. Everything else is
typed and parameterised by OPLS-AA 2024. The four torsions replace OPLS-AA 2024's for their OPLS torsion types, as the
paper's force field does (including those types elsewhere in a molecule).

usage: make_lopls.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
FID = "lopls-siu2012"
KJ = 4.184
REF = "S. W. I. Siu, K. Pluhackova, R. A. Boeckmann, J. Chem. Theory Comput. 8, 1459 (2012)"
# a carbon in a hydrocarbon chain: bonded to H and to carbons that are not aromatic and carry no heteroatom
CHAIN_NB = "!$(*~[!#6;!#1]);!$(*~[#6;$(*~[!#6;!#1])]);!$(*~a)"
EQ_CT = {"bond": "CT", "angle": "CT", "dihedral": "CT", "improper": "CT"}
EQ_CM = {"bond": "CM", "angle": "CM", "dihedral": "CM", "improper": "CM"}
EQ_HC = {"bond": "HC", "angle": "HC", "dihedral": "HC", "improper": "HC"}
# name: element, mass, charge, sigma (nm), epsilon (kJ/mol), equivalence, SMARTS, description
TYPES = [
    ("lopls_135", "C", 12.011, -0.222, 0.35, 0.276144, EQ_CT, f"[CX4H3;{CHAIN_NB}]", "L-OPLS CT in CH3 (Table 2, CTCH3)"),
    ("lopls_136", "C", 12.011, -0.148, 0.35, 0.276144, EQ_CT, f"[CX4H2;{CHAIN_NB}]", "L-OPLS CT in CH2 (Table 2, CTCH2)"),
    ("lopls_142", "C", 12.011, -0.160, 0.355, 0.317984, EQ_CM, f"[CX3H1;{CHAIN_NB}]=[CX3]", "L-OPLS CM in -CH= (Table 2, CMCH)"),
    ("lopls_140", "H", 1.008, 0.074, 0.25, 0.125520, EQ_HC, f"[H][CX4H3;{CHAIN_NB}]", "L-OPLS HC on CH3 (Table 2, HCCH3)"),
    ("lopls_140b", "H", 1.008, 0.074, 0.25, 0.110000, EQ_HC, f"[H][CX4H2;{CHAIN_NB}]", "L-OPLS HC on CH2 (Table 2, HCCH2)"),
    ("lopls_144", "H", 1.008, 0.160, 0.242, 0.125520, EQ_HC, f"[H][CX3H1;{CHAIN_NB}]=[CX3]", "L-OPLS HC on -CH= (Table 2, HCCH)"),
]
# Table 2, Fourier V1 V2 V3 (kJ/mol)
TORSIONS = [
    (["CT", "CT", "CT", "CT"], (2.697394, -0.896807, 0.74567), "hexane C5-C8-C11-C14"),
    (["??", "CM", "CM", "??"], (0.0, 51.2551, 0.0), "cis-3-hexene C5-C8-C10-C12 (X-CM-CM-X)"),
    (["CM", "CM", "CT", "CT"], (-3.368171, 1.34679, -0.4321105), "cis-2-hexene C5-C7-C9-C12"),
    (["CM", "CT", "CT", "CT"], (2.017484, 0.562197, 0.74369), "cis-2-hexene C7-C9-C12-C15"),
]
# OPLS-AA 2024 gives X-CM-CM-X again for particular ends; the hydrocarbon ones take L-OPLS's value too
CM_ENDS = [["CT", "CM", "CM", "CT"], ["CT", "CM", "CM", "HC"], ["HC", "CM", "CM", "HC"]]


def main():
    doc = {"format": "caps-forcefield", "format_version": 1, "name": "L-OPLS (Siu, Pluhackova & Boeckmann 2012) on OPLS-AA 2024",
           "source": REF + ", Table 2", "references": [REF], "extends": "oplsaa2024.json",
           "atom_types": [], "pairs": [], "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
           "notes": ["built by CAPS (bench/ff/make_lopls.py) from the paper's Table 2 as an overlay on oplsaa2024.json",
                     "L-OPLS types for CH3, CH2 and alkene CH of hydrocarbon chains (neutral groups); everything else is OPLS-AA 2024",
                     "the four torsions replace OPLS-AA 2024's for their OPLS torsion types"],
           "version": "2012", "status": "validated"}
    for name, el, mass, q, s_nm, e_kj, eq, smarts, desc in TYPES:
        doc["atom_types"].append({"name": name, "element": el, "mass": mass, "charge": q, "description": desc, "equivalence": eq,
                                  "smarts": smarts, "priority": 9})
        doc["pairs"].append({"name": name, "match": [name], "params": [e_kj / KJ, s_nm * 10],
                             "comment": f"Table 2: sigma {s_nm} nm, epsilon {e_kj:.6f} kJ/mol"})
    for match, (v1, v2, v3), what in TORSIONS:
        doc["dihedrals"].append({"name": "-".join(match), "match": match, "style": "opls", "params": [v1 / KJ, v2 / KJ, v3 / KJ, 0.0],
                                 "comment": f"Table 2 (Fourier, kJ/mol: {v1} {v2} {v3}), {what}"})
    v = TORSIONS[1][1]
    for m in CM_ENDS:
        doc["dihedrals"].append({"name": "-".join(m), "match": m, "style": "opls", "params": [v[0] / KJ, v[1] / KJ, v[2] / KJ, 0.0],
                                 "comment": "Table 2's X-CM-CM-X, for OPLS-AA 2024's hydrocarbon-specific entry"})
    json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1, ensure_ascii=False)
    entry = {"id": FID, "name": doc["name"], "version": "2012", "references": [REF], "origin": "built from source",
             "source": doc["source"] + "; overlay on oplsaa2024; bench/ff/make_lopls.py", "license": "published parameters of the paper",
             "status": "validated", "notes": "long hydrocarbon chains (alkanes, alkenes, lipid tails) with L-OPLS's charges, Lennard-Jones and torsions; "
             "the rest of a molecule as OPLS-AA 2024 (bench/ff/check_lopls.py)", "file": FID + ".json",
             "typing": {"rules": "typing/oplsaa2024.typing.json", "evidence": "OPLS-AA 2024's rules, with the six L-OPLS chain types ahead of them "
                        "(their SMARTS in the library)"},
             "counts": {"atom_types": len(doc["atom_types"]), "pairs": len(doc["pairs"]), "bonds": 0, "angles": 0, "dihedrals": len(doc["dihedrals"]),
                        "impropers": 0},
             "group": "organic", "aliases": ["lopls", "l-opls"], "list": {"label": "L-OPLS (long hydrocarbons)", "key": "lopls", "order": 7}}
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    cat["forcefields"] = [e for e in cat["forcefields"] if e["id"] != FID] + [entry]
    json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)
    print(f"{FID} written")


if __name__ == "__main__":
    main()
