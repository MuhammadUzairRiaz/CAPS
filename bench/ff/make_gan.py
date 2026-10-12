#!/usr/bin/env python3
"""GaN for CAPS from its paper: T. Wonglakhon, D. Zahn, J. Phys.: Condens. Matter 32, 205401 (2020) (CC BY), Table 1.

Formal charges Ga +3, N -3; Lennard-Jones for Ga-Ga and N-N (sigma, epsilon), Born-Mayer A exp(-r/rho) for Ga-N
(Buckingham with C = 0, no Lennard-Jones between Ga and N); Lorentz-Berthelot mixing towards other Lennard-Jones
species, as the paper proposes; short-range cut-off 12 A (the paper's).

usage: make_gan.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
REF = "T. Wonglakhon, D. Zahn, J. Phys.: Condens. Matter 32, 205401 (2020)"
EV = 23.060548
doc = {"format": "caps-forcefield", "format_version": 1, "name": "GaN (Wonglakhon & Zahn 2020)", "source": REF + ", Table 1", "references": [REF],
       "units": "real", "styles": {"pair": "buck/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
       "mixing": "arithmetic", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 12,
       "atom_types": [{"name": "Ga", "element": "Ga", "mass": 69.723, "charge": 3.0, "description": "Ga3+ (Table 1)"},
                      {"name": "N", "element": "N", "mass": 14.007, "charge": -3.0, "description": "N3- (Table 1)"}],
       "pairs": [{"name": "Ga", "match": ["Ga"], "params": [0.643 * EV, 2.390], "comment": "Lennard-Jones epsilon 0.643 eV, sigma 2.390 A"},
                 {"name": "N", "match": ["N"], "params": [1.474 * EV, 1.981], "comment": "Lennard-Jones epsilon 1.474 eV, sigma 1.981 A"},
                 {"name": "Ga-N", "match": ["Ga", "N"], "style": "buck", "params": [608.54 * EV, 0.435, 0.0], "comment": "Born-Mayer A 608.54 eV, rho 0.435 A"}],
       "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
       "notes": ["built by CAPS (bench/ff/make_gan.py) from the paper's Table 1"], "version": "2020", "typing": "../typing/gan2020.typing.json"}
json.dump(doc, open(os.path.join(DATA, "forcefields", "gan2020.json"), "w"), indent=1)
json.dump({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined",
           "description": "Ga and N ions (no bonds: the ions interact through the pair terms only)",
           "rules": [{"type": "Ga", "smarts": "[Ga]", "priority": 0, "description": "Ga3+"}, {"type": "N", "smarts": "[N]", "priority": 0, "description": "N3-"}]},
          open(os.path.join(DATA, "typing", "gan2020.typing.json"), "w"), indent=1)
print("gan2020 written")
