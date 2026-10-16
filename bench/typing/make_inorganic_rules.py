#!/usr/bin/env python3
"""CAPS typing rules for the library's water force fields (SPC/E, TIP3P 1983 and TIP3P 2004): O and H of water.
ClayFF is not here: its library and rules come from the paper (bench/ff/make_clayff_cygan2004.py,
data/typing/clayff-cygan2004.typing.json).

usage: make_inorganic_rules.py
"""
import json, os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")


def rule(t, prio, smarts, desc, **kw):
    r = {"type": t, "smarts": smarts, "priority": prio, "description": desc}
    r.update({k: v for k, v in kw.items() if v})
    return r


def water(o, h):
    return [rule(o, 0, "[OX2H2]", "water oxygen"), rule(h, 0, "[H][OX2H2]", "water hydrogen")]


def write(fid, rules, label, extra):
    ff_p = os.path.join(DATA, "forcefields", fid + ".json")
    ff = json.load(open(ff_p))
    have = {t["name"] for t in ff["atom_types"]}
    absent = sorted({r["type"] for r in rules if r["type"] not in have})
    doc = {"format": "caps-typing", "version": 1, "forcefield": label, "unknown_types": "untyped", "rules": rules,
           "description": f"CAPS rules for {label} (bench/typing/make_inorganic_rules.py)." + (" " + extra if extra else "")
                          + (" Types this file lacks: " + ", ".join(absent) + "." if absent else "")}
    out = os.path.join(DATA, "typing", fid + ".typing.json")
    json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
    ff["typing"] = f"../typing/{fid}.typing.json"
    json.dump(ff, open(ff_p, "w"), ensure_ascii=False, indent=1)
    return out, doc


if __name__ == "__main__":
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    made = {}
    have = lambda fid: os.path.exists(os.path.join(DATA, "forcefields", fid + ".json"))
    for fid, rules, label in (("spce-moltemplate", water("O", "H"), "SPC/E water"), ("tip3p-1983-moltemplate", water("O", "H"), "TIP3P water (1983)"),
                              ("tip3p-2004-moltemplate", water("O", "H"), "TIP3P water (2004)")):
        if not have(fid):
            print(f"{fid}: no data/forcefields/{fid}.json here, skipped")
            continue
        write(fid, rules, label, "")
        made[fid] = f"{len(rules)} rules"
    for e in cat["forcefields"]:
        if e["id"] in made:
            e["typing"] = {"rules": f"typing/{e['id']}.typing.json", "evidence": f"CAPS rules (bench/typing/make_inorganic_rules.py): {made[e['id']]}"}
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    for k, v in made.items():
        print(f"{k}: {v}")
