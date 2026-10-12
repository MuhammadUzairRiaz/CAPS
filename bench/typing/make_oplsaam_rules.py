#!/usr/bin/env python3
"""CAPS typing rules for OPLS-AA/M built from the Jorgensen lab's own files (bench/ff/convert_charmm_prm.py on
top_opls_aam.inp / par_opls_aam.inp: oplsaam-jorgensen), whose types are the OPLS-AA numbers with the element in front
(C224 is opls_224): CAPS's OPLS-AA chemistry rules and the protein layer of make_opls_gromacs_rules.py, written in those
names, checked against the file's own residues (bench/ff/validate_rtf_types.py --cap-links).

usage: make_oplsaam_rules.py [OUT_DIR]
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_opls_gromacs_rules as G

DATA = os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, "typing")
lib = json.load(open(os.path.join(DATA, "forcefields", "oplsaam-jorgensen.json")))
have = {t["name"]: t["element"] for t in lib["atom_types"]}
by_number = {}
for name in have:   # C224 -> 224
    m = re.match(r"^([A-Z][a-z]?)(\d+)$", name)
    if m:
        by_number[m.group(2)] = name


def aam(t):
    """an OPLS number (135, opls_135, opls_224B) as OPLS-AA/M names it (C135), or None"""
    m = re.match(r"^(?:opls_)?(\d+)[A-Z]?$", t)
    return by_number.get(str(int(m.group(1)))) if m else None


rules, lacking = [], set()
for r in G.src["rules"]:
    t = aam(r["type"])
    if not t:
        lacking.add(r["type"])
        continue
    rules.append({**r, "type": t})
for t, p, sm, d in G.ADD:
    n = aam(t)
    if n:
        rules.append({"type": n, "priority": p + 5, "smarts": sm, "description": d})
# the ions and TIP3P water of the file (CHARMM names)
for t, s, d in (("SOD", "[Na;X0]", "sodium ion"), ("CLA", "[Cl;X0]", "chloride ion"), ("OT", "[OX2H2]", "TIP3P water O"), ("HT", "[H][OX2H2]", "TIP3P water H")):
    if t in have:
        rules.append({"type": t, "priority": 20, "smarts": s, "description": d})
json.dump({"format": "caps-typing", "version": 1, "forcefield": "OPLS-AA/M (Jorgensen lab)", "unknown_types": "untyped",
           "description": "CAPS's OPLS-AA rules in OPLS-AA/M's names (element + type number), protein types as the file's residues "
                          "use them (bench/typing/make_oplsaam_rules.py)", "rules": rules},
          open(os.path.join(OUT, "oplsaam-jorgensen.typing.json"), "w"), indent=1, ensure_ascii=False)
print(f"oplsaam-jorgensen: {len(rules)} rules")
