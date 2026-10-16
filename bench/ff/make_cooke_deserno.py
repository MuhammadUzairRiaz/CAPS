#!/usr/bin/env python3
"""The Cooke-Kremer-Deserno solvent-free lipid model for CAPS, from the paper: I. R. Cooke, K. Kremer, M. Deserno,
Phys. Rev. E 72, 011506 (2005) (arXiv:cond-mat/0502418), equations 1-4, in its reduced units (sigma = epsilon = m = 1).

  lipid       one head bead H and two tail beads T, TE (TE: the end tail, which carries the straightening spring)
  Eq. 1       Weeks-Chandler-Andersen repulsion, b(head, head) = b(head, tail) = 0.95 sigma, b(tail, tail) = sigma
  Eq. 2       FENE bonds H-T and T-TE, k_bond = 30 epsilon/sigma^2, r_inf = 1.5 sigma
  Eq. 3       a harmonic spring H-TE, rest length 4 sigma, k_bend = 10 epsilon/sigma^2 for (1/2) k (r - 4 sigma)^2
              (LAMMPS harmonic K = k/2 = 5)
  Eq. 4       tail-tail attraction -epsilon for r < r_c = 2^(1/6) sigma, -epsilon cos^2(pi (r - r_c) / (2 w_c)) to r_c + w_c
              (LAMMPS cosine/squared epsilon, sigma = r_c, cutoff r_c + w_c, wca)

w_c is the model's tuning parameter; this library takes w_c = 1.6 sigma, where the paper measures the bending rigidity
(k_B T = 1.1 epsilon, a fluid bilayer); edit the T-T cutoff for another width. Bonded beads keep their repulsion (the
FENE bond carries the WCA term of its pair, b); every other pair inside a lipid interacts through the pair terms.

usage: make_cooke_deserno.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
FID = "cooke-deserno2005"
REF = "I. R. Cooke, K. Kremer, M. Deserno, Phys. Rev. E 72, 011506 (2005)"
R6 = 2 ** (1 / 6)
WC = 1.6


def main():
    hh, tt = 0.95 * R6, 1.0 * R6
    doc = {"format": "caps-forcefield", "format_version": 1, "name": "Cooke-Kremer-Deserno lipid model", "source": REF + ", equations 1-4 (arXiv:cond-mat/0502418)",
           "references": [REF], "units": "lj",
           "styles": {"pair": "cosine/squared", "bond": "fene", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"},
           "mixing": "arithmetic", "special_lj": [0, 1, 1], "special_coul": [0, 1, 1], "cutoff": tt + WC,
           "pair_settings": {"model_cutoff": True}, "angle_terms": "if_defined", "torsion_terms": "if_defined",
           "atom_types": [{"name": "H", "element": "", "mass": 1.0, "charge": 0.0, "description": "head bead"},
                          {"name": "T", "element": "", "mass": 1.0, "charge": 0.0, "description": "tail bead next to the head"},
                          {"name": "TE", "element": "", "mass": 1.0, "charge": 0.0, "description": "end tail bead (T's pair terms)", "equivalence": {"vdw": "T"}}],
           "pairs": [{"name": "H-H", "match": ["H"], "style": "cosine/squared", "params": [1.0, hh, hh, 1], "comment": "Eq. 1, WCA with b = 0.95"},
                     {"name": "H-T", "match": ["H", "T"], "style": "cosine/squared", "params": [1.0, hh, hh, 1], "comment": "Eq. 1, WCA with b = 0.95"},
                     {"name": "T-T", "match": ["T"], "style": "cosine/squared", "params": [1.0, tt, tt + WC, 1],
                      "comment": f"Eq. 1 (WCA, b = 1) and Eq. 4 (cos^2 attraction, w_c = {WC})"}],
           "bonds": [{"name": "H-T", "match": ["H", "T"], "style": "fene", "params": [30.0, 1.5, 1.0, 0.95], "comment": "Eq. 2: k 30, r_inf 1.5; the pair's WCA (b 0.95)"},
                     {"name": "T-TE", "match": ["T", "TE"], "style": "fene", "params": [30.0, 1.5, 1.0, 1.0], "comment": "Eq. 2: k 30, r_inf 1.5; the pair's WCA (b 1)"},
                     {"name": "H-TE", "match": ["H", "TE"], "style": "harmonic", "params": [5.0, 4.0], "comment": "Eq. 3: (1/2) 10 (r - 4)^2"}],
           "angles": [], "dihedrals": [], "impropers": [],
           "bead_templates": {"lipid": "[H]1[T][TE]1"},
           "notes": ["built by CAPS (bench/ff/make_cooke_deserno.py) from the paper's equations 1-4",
                     f"w_c = {WC} sigma (the paper's bending-rigidity state at k_B T = 1.1 epsilon); the T-T cutoff is r_c + w_c",
                     "reduced units (sigma = epsilon = m = 1): CAPS reads the numbers as A, kcal/mol, g/mol; scale lengths, energies and times yourself"],
           "version": "2005", "status": "validated", "typing": f"../typing/{FID}.typing.json"}
    json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1, ensure_ascii=False)
    typing = {"format": "caps-typing", "version": 1, "forcefield": doc["name"], "coarse_grained": True, "unknown_types": "untyped",
              "description": "CAPS rules for the Cooke-Kremer-Deserno lipid (bench/ff/make_cooke_deserno.py): beads typed by name (H, T, TE)",
              "rules": [{"type": t, "smarts": "*", "atom_name": t, "priority": 0, "description": f"a bead named {t}"} for t in ("H", "T", "TE")]}
    json.dump(typing, open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1, ensure_ascii=False)
    entry = {"id": FID, "name": doc["name"], "version": "2005", "references": [REF], "origin": "built from source",
             "source": doc["source"] + "; bench/ff/make_cooke_deserno.py", "license": "published model of the paper (its equations)",
             "status": "validated", "notes": f"three-bead lipids, implicit solvent, w_c = {WC}; LAMMPS equals CAPS's energy (bench/ff/check_data_lammps.py)",
             "file": FID + ".json", "typing": {"rules": f"typing/{FID}.typing.json", "evidence": "beads typed by name"},
             "counts": {"atom_types": 3, "pairs": 3, "bonds": 3, "angles": 0, "dihedrals": 0, "impropers": 0}}
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    cat["forcefields"] = [e for e in cat["forcefields"] if e["id"] != FID] + [entry]
    json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)
    print(f"{FID} written")


if __name__ == "__main__":
    main()
