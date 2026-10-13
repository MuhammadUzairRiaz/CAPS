#!/usr/bin/env python3
"""The DES components of CL&P's des.ff (clandp-des) typed by CAPS from SMILES: each must type completely, carry its
integer charge with des.ff's own type charges and find every bonded term; the ones des.ff cannot make neutral must stay
untyped. Prints the types of each component.

usage: PYTHONPATH=data/python python3 bench/ff/check_clandp_des_components.py   (exit 1 on a failure)
"""
import collections, sys

import caps

TYPED = [  # name, SMILES, charge
    ("cholinium", "C[N+](C)(C)CCO", 1), ("tetramethylammonium", "C[N+](C)(C)C", 1),
    ("tetraethylammonium", "CC[N+](CC)(CC)CC", 1), ("tetrabutylammonium", "CCCC[N+](CCCC)(CCCC)CCCC", 1),
    ("triethylammonium", "CC[NH+](CC)CC", 1), ("N,N-diethylethanolammonium", "CC[NH+](CC)CCO", 1),
    ("N-ethyl-N,N-dimethylethanolammonium", "CC[N+](C)(C)CCO", 1),
    ("tetrabutylphosphonium", "CCCC[P+](CCCC)(CCCC)CCCC", 1), ("betaine", "C[N+](C)(C)CC(=O)[O-]", 0),
    ("chloride", "[Cl-]", -1), ("bromide", "[Br-]", -1), ("acetate", "CC(=O)[O-]", -1), ("propanoate", "CCC(=O)[O-]", -1),
    ("isobutyrate", "CC(C)C(=O)[O-]", -1),
    ("urea", "NC(N)=O", 0), ("ethanol", "CCO", 0), ("2-propanol", "CC(C)O", 0), ("ethylene glycol", "OCCO", 0),
    ("1,2-propanediol", "CC(O)CO", 0), ("1,4-butanediol", "OCCCCO", 0), ("glycerol", "OCC(O)CO", 0),
    ("acetic acid", "CC(=O)O", 0), ("propanoic acid", "CCC(=O)O", 0), ("decanoic acid", "CCCCCCCCCC(=O)O", 0),
    ("succinic acid", "OC(=O)CCC(=O)O", 0), ("glutaric acid", "OC(=O)CCCC(=O)O", 0),
    ("levulinic acid", "CC(=O)CCC(=O)O", 0), ("acetone", "CC(C)=O", 0), ("2-butanone", "CCC(C)=O", 0),
]
UNTYPED = [  # des.ff's charges cannot make these neutral: CAPS leaves an atom untyped
    ("methanol", "CO"), ("malonic acid", "OC(=O)CC(=O)O"), ("oxalic acid", "OC(=O)C(=O)O"), ("lactic acid", "CC(O)C(=O)O"),
    ("benzoic acid", "OC(=O)c1ccccc1"), ("tert-butanol", "CC(C)(C)O"), ("ethylammonium", "CC[NH3+]"),
]
bad = 0
for name, smi, q in TYPED:
    r = caps.build.smiles(smi, forcefield="uff").field.assign("clandp-des")
    types = collections.Counter(a["type"] for a in r["atoms"])
    ok = r["complete"] and r["untyped"] == 0 and not r["missing"] and abs(r["net_charge"] - q) < 1e-9
    bad += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {name:38s} q {r['net_charge']:+.3f}  " + " ".join(f"{t}×{n}" for t, n in sorted(types.items()))
          + ("" if ok else f"  untyped {r['untyped']} missing {r['missing'][:4]}"))
for name, smi in UNTYPED:
    r = caps.build.smiles(smi, forcefield="uff").field.assign("clandp-des")
    ok = r["untyped"] > 0
    bad += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {name:38s} left untyped ({r['untyped']} atoms)")
print(f"{len(TYPED)} components typed, neutral or at their ion charge, complete; {len(UNTYPED)} left untyped; {bad} failures")
sys.exit(1 if bad else 0)
