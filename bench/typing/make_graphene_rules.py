#!/usr/bin/env python3
"""CAPS typing rules for the library's graphene force field (moltemplate graphene.lt): one Lennard-Jones carbon, no
bonded terms (the source's sheet has no bonds: it is held in place). A carbon whose neighbours are all carbons takes the
type; "bonds": "defined" drops a built sheet's C-C bonds (the force field defines none) and "exclude_pairs" keeps its carbons
from interacting with each other (the source's neigh_modify exclude group Cgraphene Cgraphene): the sheet only interacts
with other molecules, and should be held in place (a held molecule), as the source uses it.

usage: make_graphene_rules.py
"""
import json, os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
FID = "graphene-moltemplate"

if __name__ == "__main__":
    doc = {"format": "caps-typing", "version": 1, "forcefield": "Graphene (LJ carbon)", "unknown_types": "untyped", "bonds": "defined",
           "exclude_pairs": [["C", "C"]],
           "description": "CAPS rules for the LJ graphene carbon (bench/typing/make_graphene_rules.py): carbons bonded only to carbons; "
                          "the sheet's bonds are dropped (the force field has none: hold the sheet in place)",
           "rules": [{"type": "C", "smarts": "[#6;!$(*~[!#6])]", "priority": 0, "description": "a carbon bonded only to carbons (a graphene sheet)"}]}
    json.dump(doc, open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), ensure_ascii=False, indent=1)
    p = os.path.join(DATA, "forcefields", FID + ".json")
    ff = json.load(open(p))
    ff["typing"] = f"../typing/{FID}.typing.json"
    note = ("no bonded terms (the source's sheet has no bonds): CAPS drops a built sheet's bonds and its carbons do not interact "
            "with each other (the source's neigh_modify exclude group Cgraphene Cgraphene); hold the sheet in place")
    if note not in ff.setdefault("notes", []):
        ff["notes"].append(note)
    json.dump(ff, open(p, "w"), ensure_ascii=False, indent=1)
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for e in cat["forcefields"]:
        if e["id"] == FID:
            e["notes"] = "LJ carbon for graphene sheets, non-bonded only (hold the sheet in place)"
            e["typing"] = {"rules": f"typing/{FID}.typing.json", "evidence": "CAPS rule (bench/typing/make_graphene_rules.py): carbons bonded only to carbons"}
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    print("graphene rules written")
