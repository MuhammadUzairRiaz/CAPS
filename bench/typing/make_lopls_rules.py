#!/usr/bin/env python3
"""L-OPLS (Siu, Pluhackova & Böckmann, J. Chem. Theory Comput. 8, 1459 (2012); Pluhackova et al., J. Phys. Chem. B 119,
15287 (2015)) as a complete force field: the library's overlay file "extends" OPLS-AA 2024 (the base is loaded first,
the overlay's types and torsions on top), and these rules are OPLS-AA 2024's with the L-OPLS types where the papers use
them, in long hydrocarbon stretches:

  CH3 / CH2 (54L / 57L, H 60LCH3 / 60LCH2)   carbons whose heavy neighbours are carbons without heteroatoms, so the
                                             group stays neutral (-0.222 + 3 x 0.074, -0.148 + 2 x 0.074)
  RHC= / R2C= (142L with H 144L / 141L)      alkene carbons of hydrocarbons (neutral: -0.16 + 0.16, 0)

H2C= keeps OPLS-AA's 143 / 144: L-OPLS's 143L (-0.23) with two 144L hydrogens (+0.16) would leave +0.09 e. L-OPLS's
ester and alcohol types come with charges that balance only with neighbouring types the overlay does not define, so
those groups keep their OPLS-AA types. The same is written for the overlay on OPLS-AA 2008 (BOSS numbers 80L, 81L ...),
over the 2008 rules (bench/typing/make_oplsaa2008_rules.py).

usage: make_lopls_rules.py [DATA_DIR]
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data")

HYDROCARBON = "!$(*~[!#6;!#1]);!$(*~[#6]~[!#6;!#1])"   # no heteroatom on the carbon or on its carbon neighbours
# the L-OPLS types under each base table's numbers: CH3, CH2, H on CH3, H on CH2, RHC=, R2C=, H on RHC=
NUMBERS = {"2024": ("54L", "57L", "60LCH3", "60LCH2", "142L", "141L", "144L"), "2008": ("80L", "81L", "85LCH3", "85LCH2", "87L", "86L", "89L")}


def lopls(v):
    ch3, ch2, h3, h2, rhc, r2c, hc = NUMBERS[v]
    return [
        (ch3, 9, f"[CX4H3;{HYDROCARBON}]", "L-OPLS CH3 of a hydrocarbon chain"),
        (ch2, 9, f"[CX4H2;{HYDROCARBON}]", "L-OPLS CH2 of a hydrocarbon chain"),
        (h3, 9, f"[H][CX4H3;{HYDROCARBON}]", "L-OPLS H on a chain CH3"),
        (h2, 9, f"[H][CX4H2;{HYDROCARBON}]", "L-OPLS H on a chain CH2"),
        (rhc, 9, f"[CX3H1;{HYDROCARBON}]=[CX3;{HYDROCARBON}]", "L-OPLS alkene carbon RHC="),
        (r2c, 9, f"[CX3H0;{HYDROCARBON}]=[CX3;{HYDROCARBON}]", "L-OPLS alkene carbon R2C="),
        (hc, 9, f"[H][CX3H1;{HYDROCARBON}]=[CX3;{HYDROCARBON}]", "L-OPLS H on an RHC= carbon"),
    ]


if __name__ == "__main__":
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for v in ("2024", "2008"):
        base_id, target = f"oplsaa{v}-moltemplate", f"loplsaa{v}-moltemplate"
        base = json.load(open(os.path.join(DATA, "typing", base_id + ".typing.json")))
        L = lopls(v)
        rules = list(base["rules"]) + [{"type": t_, "smarts": s, "priority": p, "description": d} for t_, p, s, d in L]
        label = f"L-OPLS (on OPLS-AA {v})"
        doc = dict(base, rules=rules, forcefield=label, unknown_types="untyped",
                   description=f"OPLS-AA {v}'s rules with the L-OPLS types for hydrocarbon CH3 / CH2 and RHC= / R2C= carbons and their hydrogens "
                               "(bench/typing/make_lopls_rules.py); other groups keep their OPLS-AA types.")
        out = os.path.join(DATA, "typing", target + ".typing.json")
        json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
        ff_p = os.path.join(DATA, "forcefields", target + ".json")
        ff = json.load(open(ff_p))
        ff["extends"] = base_id + ".json"
        ff["typing"] = f"../typing/{target}.typing.json"
        ff["name"] = label
        ff["validation"] = f"overlay on OPLS-AA {v}, loaded as one force field (\"extends\")."
        json.dump(ff, open(ff_p, "w"), ensure_ascii=False, indent=1)
        for e in cat["forcefields"]:
            if e["id"] == target:
                e["name"] = label
                e["notes"] = f"L-OPLS types and torsions on OPLS-AA {v} (the base is loaded with it)"
                e["typing"] = {"rules": f"typing/{target}.typing.json", "evidence": f"OPLS-AA {v} rules with L-OPLS hydrocarbon types (bench/typing/make_lopls_rules.py)"}
        print(f"{out}: {len(rules)} rules ({len(L)} L-OPLS)")
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
