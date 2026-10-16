#!/usr/bin/env python3
"""CAPS extensions to the library's water force fields (documented in each file's notes):

  water        SPC/E, TIP3P (1983) and TIP3P (2004): charges, and the O-H bond and H-O-H angle terms of the moltemplate
               files (rigid models there, kept by SHAKE; CAPS has no constraints and runs them flexible, as the notes say)

usage: extend_inorganic_water.py
"""
import json, os

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FF = os.path.join(ROOT, "data", "forcefields")


def note(ff, text):
    ff.setdefault("notes", [])
    if text not in ff["notes"]:
        ff["notes"].append(text)


WATER = {"spce-moltemplate": (-0.8476, 0.4238, 600.0, 1.0, 75.0, 109.47), "tip3p-1983-moltemplate": (-0.834, 0.417, 450.0, 0.9572, 55.0, 104.52),
         "tip3p-2004-moltemplate": (-0.830, 0.415, 450.0, 0.9572, 55.0, 104.52)}
for fid, (qo, qh, kb, r0, ka, a0) in WATER.items():
    p = os.path.join(FF, fid + ".json")
    ff = json.load(open(p))
    for t in ff["atom_types"]:
        t["charge"] = qo if t["name"] == "O" else qh
    if not ff.get("bonds"):
        ff["bonds"] = [{"name": "O-H", "match": ["O", "H"], "style": "harmonic", "params": [kb, r0], "comment": "moltemplate's OH bond (rigid there, by SHAKE)"}]
    if not ff.get("angles"):
        ff["angles"] = [{"name": "H-O-H", "match": ["H", "O", "H"], "style": "harmonic", "params": [ka, a0], "comment": "moltemplate's HOH angle (rigid there, by SHAKE)"}]
    note(ff, "CAPS extension: charges and the O-H / H-O-H terms of the moltemplate file; the model is rigid (SHAKE) there, CAPS has no constraints "
             "and runs it flexible with these terms")
    json.dump(ff, open(p, "w"), ensure_ascii=False, indent=1)
    print(f"{fid}: charges O {qo} H {qh}; bond {kb} {r0}; angle {ka} {a0}")
