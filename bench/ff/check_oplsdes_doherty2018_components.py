#!/usr/bin/env python3
"""The molecules of OPLS-DES (Doherty & Acevedo 2018; oplsdes-doherty2018) typed by CAPS from SMILES: each must be
typed as the SI's Schemes S1 and S2 type it, carry its charge (choline +0.8, chloride -0.8, every HBD 0: the paper's
+-0.8 e scaling) and find every bonded term. Molecules outside the paper's set must stay untyped. Prints, per molecule,
how many atoms carry the scheme's type (urea's HC, the H cis to O, has HT's terms and is typed HT: counted by name and
by parameters).

usage: PYTHONPATH=data/python python3 bench/ff/check_oplsdes_doherty2018_components.py   (exit 1 on a failure)
"""
import collections, sys

import caps

FF = "oplsdes-doherty2018"
TYPED = [  # name, SMILES, charge, the scheme's types and counts
    ("choline", "C[N+](C)(C)CCO", 0.8, {"NA": 1, "CA": 3, "HA": 9, "CS": 1, "HS": 2, "CW": 1, "HW": 2, "OY": 1, "HY": 1}),
    ("chloride", "[Cl-]", -0.8, {"Cl": 1}),
    ("ethylene glycol", "OCCO", 0, {"OG_etg": 2, "HO_etg": 2, "CG_etg": 2, "HG_etg": 4}),
    ("glycerol", "OCC(O)CO", 0, {"CM_gly": 1, "OM_gly": 1, "HM_gly": 1, "HZ_gly": 1, "CB_gly": 2, "HC_gly": 4, "OH_gly": 2,
                                 "HO_gly": 2}),
    ("levulinic acid", "CC(=O)CCC(=O)O", 0, {"CB_lev": 1, "OB_lev": 1, "OH_lev": 1, "HO_lev": 1, "CZ_lev": 1, "HZ_lev": 2,
                                             "CF_lev": 1, "HF_lev": 2, "CD_lev": 1, "OC_lev": 1, "CT_lev": 1, "HT_lev": 3}),
    ("malonic acid", "OC(=O)CC(=O)O", 0, {"CD_mal": 2, "OD_mal": 2, "OH_mal": 2, "HO_mal": 2, "CT_mal": 1, "HC_mal": 2}),
    ("oxalic acid", "OC(=O)C(=O)O", 0, {"CD_ox": 2, "OD_ox": 2, "OH_ox": 2, "HO_ox": 2}),
    ("phenol", "Oc1ccccc1", 0, {"CO_phe": 1, "CB_phe": 5, "HB_phe": 5, "OH_phe": 1, "HO_phe": 1}),
    ("urea", "NC(N)=O", 0, {"C_urea": 1, "O_urea": 1, "N_urea": 2, "HT_urea": 2, "HC_urea": 2}),
]
SAME = {"HC_urea": "HT_urea"}   # the scheme's type -> the type with identical terms CAPS assigns
UNTYPED = [  # not in the paper's set: no rule applies
    ("1,2-propanediol", "CC(O)CO"), ("1,4-butanediol", "OCCCCO"), ("acetic acid", "CC(=O)O"), ("succinic acid", "OC(=O)CCC(=O)O"),
    ("methanol", "CO"), ("tetramethylammonium", "C[N+](C)(C)C"), ("betaine", "C[N+](C)(C)CC(=O)[O-]"), ("bromide", "[Br-]"),
    ("benzoic acid", "OC(=O)c1ccccc1"), ("p-cresol", "Cc1ccc(O)cc1"), ("acetamide", "CC(N)=O"), ("methylurea", "CNC(N)=O"),
]
bad = 0
tot_name = tot_par = tot = 0
for name, smi, q, want in TYPED:
    r = caps.build.smiles(smi, forcefield="uff").field.assign(FF)
    got = collections.Counter(a["type"] for a in r["atoms"])
    n = sum(want.values())
    by_name = sum(min(got[t], c) for t, c in want.items())
    par = collections.Counter()
    for t, c in want.items():
        par[SAME.get(t, t)] += c
    by_par = sum(min(got[t], c) for t, c in par.items())
    ok = r["complete"] and r["untyped"] == 0 and not r["missing"] and abs(r["net_charge"] - q) < 1e-9 and by_par == n == len(r["atoms"])
    bad += not ok
    tot += n
    tot_name += by_name
    tot_par += by_par
    print(f"{'ok  ' if ok else 'FAIL'} {name:16s} q {r['net_charge']:+.4f}  scheme types {by_name}/{n} by name, {by_par}/{n} by "
          f"parameters  " + " ".join(f"{t}x{c}" for t, c in sorted(got.items()))
          + ("" if ok else f"  untyped {r['untyped']} missing {r['missing'][:4]}"))
for name, smi in UNTYPED:
    r = caps.build.smiles(smi, forcefield="uff").field.assign(FF)
    ok = r["untyped"] == len(r["atoms"])
    bad += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {name:16s} left untyped ({r['untyped']} / {len(r['atoms'])} atoms)")
print(f"{len(TYPED)} molecules: {tot_name}/{tot} atoms typed as the schemes name them, {tot_par}/{tot} by parameters; "
      f"{len(UNTYPED)} others left untyped; {bad} failures")
sys.exit(1 if bad else 0)
