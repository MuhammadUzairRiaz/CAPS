#!/usr/bin/env python3
"""CAPS typing rules for DREIDING (Mayo, Olafson & Goddard, J. Phys. Chem. 94, 8897 (1990)) in both library files.

DREIDING types an atom by its element and hybridisation only: _3 sp3, _2 sp2, _1 sp, _R an sp2 atom in resonance
(aromatic rings; amide and aniline N, whose lone pair is conjugated). Its bond constant is 700 n kcal/mol/Å² for bond
order n (350 n as LAMMPS's harmonic K). Type names only carry the element and hybridisation, so the files add variants
that give a bond between two sp2 atoms its real order: moltemplate's C_2 / C_2_b1 / C_2_b2, C_R / C_R_b1 ... (a single bond
between two conjugated atoms is order 1) and the other file's C_2 / C_2S, C_R / C_RS (such a bond is order 1.5). CAPS
types by the paper's rules, then picks the variants that give every bond the constant of its order (bond_k_per_order,
conjugated_single_order; the core's refine_bond_order_variants).

moltemplate also encodes, as _dN, the number N of neighbours besides the bond partner that DREIDING's torsion barrier is
divided by (pyridine N_R_d1, amide N_R_d2; phosphine P_3_d2, phosphate P_3_d3): the rules read it from the atom.
moltemplate's hydrogen-bonding variants (H_HB, _ha / _hd) are for the explicit hbond/dreiding term and are not
assigned; the other file bonds O and N only to its H___A, which it gets.

usage: make_dreiding_rules.py [DATA_DIR]
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data")

# (canonical type, priority, SMARTS, description); canonical names are renamed per file below
BASE = [
    ("H", 0, "[#1]", "hydrogen"),
    ("C_3", 0, "[CX4]", "sp3 C"),
    ("C_2", 0, "[CX3]", "sp2 C"),
    ("C_R", 1, "[c]", "resonant C (aromatic)"),
    ("C_1", 0, "[CX2]", "sp C"),
    ("N_3", 0, "[NX4,NX3]", "sp3 N"),
    ("N_R/1", 1, "[n;X2]", "resonant N with two neighbours (pyridine)"),
    ("N_R/2", 1, "[n;X3]", "resonant N with three neighbours (pyrrole)"),
    ("N_R/2", 2, "[NX3;$(N[CX3]=[O,S,N]),$(Nc)]", "resonant N: amide, urea, aniline (conjugated lone pair)"),
    ("N_2/1", 1, "[NX2]=*", "sp2 N with two neighbours (imine, azo)"),
    ("N_2/2", 1, "[NX3]=*", "sp2 N with three neighbours (nitro)"),
    ("N_1", 1, "[NX1]#*", "sp N (nitrile)"),
    ("O_3", 0, "[OX2]", "sp3 O (ether, alcohol, ester)"),
    ("O_2", 1, "[OX1]=*", "sp2 O (carbonyl)"),
    ("O_2", 1, "[OX1]~[NX3](~[OX1])", "O of a nitro group"),
    ("O_R", 1, "[o]", "resonant O (furan)"),
    ("S_3", 0, "[SX2,SX4,SX3]", "S"),
    ("F", 0, "[F]", "F"),
    ("Cl", 0, "[Cl]", "Cl"),
    ("Br", 0, "[Br]", "Br"),
    ("I", 0, "[I]", "I"),
    ("Si_3", 0, "[Si]", "Si"),
    ("P_3/2", 0, "[PX3]", "P with three neighbours"),
    ("P_3/3", 0, "[PX4]", "P with four neighbours"),
    ("B_3", 0, "[BX4]", "sp3 B"),
    ("B_2/1", 0, "[BX2]", "sp2 B with two neighbours"),
    ("B_2/2", 0, "[BX3]", "sp2 B with three neighbours"),
    ("Na", 0, "[Na]", "Na"),
    ("Ca", 0, "[Ca]", "Ca"),
    ("Zn", 0, "[Zn]", "Zn"),
    ("Fe", 0, "[Fe]", "Fe"),
]

FILES = {
    "dreiding-moltemplate": {
        "label": "DREIDING (moltemplate)",
        "name": {"N_R/1": "N_R_d1", "N_R/2": "N_R_d2", "N_2/1": "N_2_d1", "N_2/2": "N_2_d2", "P_3/2": "P_3_d2", "P_3/3": "P_3_d3",
                 "B_2/1": "B_2_d1", "B_2/2": "B_2_d2"},
        "variants": {"C_2": ["C_2_b1", "C_2_b2"], "C_1": ["C_1_b1"], "C_R": ["C_R_b1"], "N_R_d2": ["N_R_b1_d2"],
                     "N_2_d1": ["N_2_b1_d1", "N_2_b2_d1"], "N_2_d2": ["N_2_b1_d2", "N_2_b2_d2"], "O_2": ["O_2_b1", "O_2_b2"],
                     "B_2_d1": ["B_2_b1_d1", "B_2_b2_d1"], "B_2_d2": ["B_2_b1_d2", "B_2_b2_d2"]},
        "conjugated_single_order": 1.0,
    },
    "dreiding": {
        "label": "DREIDING",
        "name": {"H": "H_", "N_R/1": "N_R", "N_R/2": "N_R", "N_2/1": "N_2", "N_2/2": "N_2", "F": "F_", "I": "I_", "Si_3": "Si3",
                 "P_3/2": "P_3", "P_3/3": "P_3", "B_2/1": "B_2", "B_2/2": "B_2"},
        "variants": {"C_2": ["C_2S"], "C_R": ["C_RS"]},
        "conjugated_single_order": 1.5,
        # this file bonds O and N only to the paper's hydrogen-bonding H (H___A, its H__HB)
        "extra": [("H__HB", 1, "[#1][O,N]", "H on O or N (DREIDING's H___A)")],
    },
}

if __name__ == "__main__":
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for target, spec in FILES.items():
        ff_p = os.path.join(DATA, "forcefields", target + ".json")
        ff = json.load(open(ff_p))
        have = {t["name"] for t in ff["atom_types"]}
        rules = [{"type": spec["name"].get(t, t), "smarts": s, "priority": pr, "description": d} for t, pr, s, d in BASE + spec.get("extra", [])]
        absent = sorted({r["type"] for r in rules if r["type"] not in have})
        doc = {"format": "caps-typing", "version": 1, "forcefield": spec["label"],
               "description": "CAPS rules for DREIDING (Mayo, Olafson & Goddard 1990): element and hybridisation, then the variants that give each "
                              "bond between sp2 atoms the constant of its order. Types this file lacks leave their atoms untyped: " + (", ".join(absent) or "none") + ".",
               "unknown_types": "untyped", "variants": spec["variants"], "bond_k_per_order": 350.0,
               "conjugated_single_order": spec["conjugated_single_order"], "rules": rules}
        out = os.path.join(DATA, "typing", target + ".typing.json")
        json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
        ff["typing"] = f"../typing/{target}.typing.json"
        json.dump(ff, open(ff_p, "w"), ensure_ascii=False, indent=1)
        for e in cat["forcefields"]:
            if e["id"] == target:
                e["typing"] = {"rules": f"typing/{target}.typing.json", "evidence": "CAPS DREIDING rules with bond-order variants (bench/typing/make_dreiding_rules.py)"
                               + ("" if target == "dreiding-moltemplate" else "; this file lists torsions for common neighbourhoods with the barrier already divided, "
                                  "so some sp3–sp2 torsions (e.g. next to Cl on an sp2 carbon) are missing: DREIDING (moltemplate) applies the paper's rule to every bond")}
        print(f"{out}: {len(rules)} rules; absent from {target}: {', '.join(absent) or 'none'}")
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
