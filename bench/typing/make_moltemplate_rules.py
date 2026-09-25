#!/usr/bin/env python3
"""Typing rules for the moltemplate force fields whose type names are those of a force field CAPS already types:
GAFF (gaff.lt) and GAFF2 (gaff2.lt) take CAPS's antechamber-ordered GAFF / GAFF2 rules. Every rule is kept in its place
(antechamber's table is ordered: the first matching rule wins), with "unknown_types": "untyped", so a type the moltemplate
file lacks leaves its atoms untyped (reported by the Field page) instead of letting a later, more general rule take them.

usage: make_moltemplate_rules.py [DATA_DIR]
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data")
PAIRS = {"gaff-moltemplate": ("gaff-amber16", "GAFF (moltemplate gaff.lt)"),
         "gaff2-moltemplate": ("gaff-amber25", "GAFF2 (moltemplate gaff2.lt)")}
# types a newer version splits off from one the file has: GAFF2 from AmberTools 21 on types sp3 carbons in 5- and
# 6-membered rings c5 / c6; gaff2.lt predates that and types them c3, with c3's parameters
VERSION = {"gaff2-moltemplate": {"c5": "c3", "c6": "c3"}}

cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
cat = json.load(open(cat_p))
for target, (source, label) in PAIRS.items():
    src = json.load(open(os.path.join(DATA, "typing", source + ".typing.json")))
    ff_p = os.path.join(DATA, "forcefields", target + ".json")
    ff = json.load(open(ff_p))
    have = {t["name"] for t in ff["atom_types"]}
    vmap = VERSION.get(target, {})
    rules = [dict(r, type=vmap.get(r["type"], r["type"]), **({"description": (r.get("description", "") + f" ({r['type']} in newer versions)").strip()} if r["type"] in vmap else {}))
             for r in src["rules"]]
    pairs = [p for p in src.get("pairs", []) if not any(x in vmap for x in p)]
    absent = sorted({r["type"] for r in rules if r["type"] not in have})
    doc = dict(src)
    doc["rules"], doc["pairs"] = rules, pairs
    doc["forcefield"] = label
    doc["description"] = (src["description"] + f" Shared with {source}; types this file lacks ({', '.join(absent) or 'none'}) leave their atoms untyped.")
    doc["unknown_types"] = "untyped"
    out = os.path.join(DATA, "typing", target + ".typing.json")
    json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
    ff["typing"] = f"../typing/{target}.typing.json"
    json.dump(ff, open(ff_p, "w"), ensure_ascii=False, indent=1)
    for e in cat["forcefields"]:
        if e["id"] == target:
            e["typing"] = {"rules": f"typing/{target}.typing.json", "evidence": f"CAPS's {source} rules (antechamber order); types not in this file: {', '.join(absent) or 'none'}"}
    print(f"{out}: {len(doc['rules'])} rules; types absent from {target}: {', '.join(absent) or 'none'}")
json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
