#!/usr/bin/env python3
"""Noble gases for CAPS: one Lennard-Jones site per atom, Ne, Ar, Kr and Xe, from G. Rutkai, M. Thol, R. Span, J. Vrabec,
Mol. Phys. 115, 1104 (2017), Table 1 (the full Lennard-Jones column; the authors' copy, ~/project/ff-sources/papers/gases):

          epsilon/k_B (K)   sigma (A)
    Ne        33.921          2.801
    Ar       116.79           3.3952
    Kr       162.58           3.6274
    Xe       226.14           3.949

The paper fits each gas on its own (vapour pressure and saturated liquid density through the Lennard-Jones equation of
state); unlike pairs take the Lorentz-Berthelot rule, CAPS's choice for mixtures. Masses: IUPAC standard atomic weights.
Molecular gases (CO2, N2, O2, CH4, NH3, H2S) are in the TraPPE libraries (trappe-ua-umn).

usage: make_noble_gases.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
FID = "noble-gases-rutkai2017"
KB = 0.0019872043   # kcal/mol/K (CODATA 2018 k_B N_A)
REF = "G. Rutkai, M. Thol, R. Span, J. Vrabec, Mol. Phys. 115, 1104 (2017), doi:10.1080/00268976.2016.1246760"
# element: (mass, epsilon/k_B in K, sigma in A)
GASES = {"Ne": (20.1797, 33.921, 2.801), "Ar": (39.948, 116.79, 3.3952), "Kr": (83.798, 162.58, 3.6274), "Xe": (131.293, 226.14, 3.949)}
NAMES = {"Ne": "neon", "Ar": "argon", "Kr": "krypton", "Xe": "xenon"}


def main():
    doc = {"format": "caps-forcefield", "format_version": 1, "name": "Noble gases (Rutkai et al. 2017)", "source": REF + ", Table 1 (LJ column)",
           "references": [REF], "units": "real",
           "styles": {"pair": "lj/cut", "bond": "harmonic", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"},
           "mixing": "arithmetic", "special_lj": [0, 0, 0], "special_coul": [0, 0, 0], "cutoff": 12,
           "atom_types": [{"name": el, "element": el, "mass": m, "charge": 0.0, "description": f"{NAMES[el]} (Table 1)"} for el, (m, _, _) in GASES.items()],
           "pairs": [{"name": el, "match": [el], "params": [ek * KB, s], "comment": f"epsilon/k_B {ek} K, sigma {s} A (Table 1, LJ)"}
                     for el, (_, ek, s) in GASES.items()],
           "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
           "notes": ["built by CAPS (bench/ff/make_noble_gases.py) from the paper's Table 1 (the full Lennard-Jones potential, used with "
                     "long-range corrections in the paper; the truncated-shifted set is not this one)",
                     "unlike pairs: Lorentz-Berthelot, CAPS's choice (the paper fits each gas alone)"],
           "version": "2017", "status": "validated", "typing": f"../typing/{FID}.typing.json"}
    json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1, ensure_ascii=False)
    typing = {"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped",
              "description": "CAPS rules for the noble gases (bench/ff/make_noble_gases.py): one unbonded atom each",
              "rules": [{"type": el, "smarts": f"[{el};X0]", "priority": 0, "description": NAMES[el]} for el in GASES]}
    json.dump(typing, open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1, ensure_ascii=False)
    entry = {"id": FID, "name": doc["name"], "version": "2017", "references": [REF], "origin": "built from source",
             "source": doc["source"] + "; bench/ff/make_noble_gases.py", "license": "published parameters of the paper", "status": "validated",
             "notes": "Ne, Ar, Kr, Xe as single Lennard-Jones sites; LAMMPS equals CAPS's energy (bench/ff/check_water_libraries.py)",
             "file": FID + ".json", "typing": {"rules": f"typing/{FID}.typing.json", "evidence": "each noble-gas atom by its element"},
             "counts": {"atom_types": len(doc["atom_types"]), "pairs": len(doc["pairs"]), "bonds": 0, "angles": 0, "dihedrals": 0, "impropers": 0},
             "group": "gas", "aliases": ["noble-gases"], "list": {"label": "Noble gases (Ne, Ar, Kr, Xe)", "key": "noble_gases", "order": 33}}
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    cat["forcefields"] = [e for e in cat["forcefields"] if e["id"] != FID] + [entry]
    json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)
    print(f"{FID} written")


if __name__ == "__main__":
    main()
