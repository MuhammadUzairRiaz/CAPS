#!/usr/bin/env python3
"""Write CAPS typing rules for the INTERFACE force field on PCFF and on CVFF (data/typing/iff-pcff.typing.json,
iff-cvff.typing.json): CAPS's own PCFF / CVFF rules for every organic atom (data/typing/pcff-frc.typing.json,
cvff.typing.json, as the catalogue states), then IFF's rules for what IFF adds: silica (sc4, oc23, oc24, hoy, Heinz et
al., Chem. Mater. 26, 2647 (2014)), the sodium ion and the metals (Heinz et al., J. Phys. Chem. C 112, 17281 (2008)),
which stay unbonded. IFF's rules carry priority 6, above every PCFF / CVFF rule (at most 5), so a silica Si is sc4 and
not PCFF's siloxane sio.

Run it after editing the PCFF or CVFF rules, so the IFF files never fall behind them.
usage: make_iff_rules.py [OUTDIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA, "typing")

SILICA_SI = "[SiX4;!$(*[OX2][#6])](~[#8])(~[#8])(~[#8])~[#8]"   # no O-C: an alkoxysilane (TEOS) keeps the PCFF / CVFF Si, IFF has no sc4-o term
UNBONDED = "!$(*~[#1,#5,#6,#7,#8,#9,#14,#15,#16,#17,#35,#53])"
METALS = ["Ag", "Al", "Au", "Cr", "Cu", "Fe", "Mo", "Ni", "Pb", "Pd", "Pt", "Sn", "W"]
OWN = [
    {"type": "sc4", "smarts": SILICA_SI, "priority": 6, "description": "IFF: Si in silica (four O, none of them on carbon)"},
    {"type": "oc23", "smarts": f"[OX2]([Si;$({SILICA_SI})])[Si;$({SILICA_SI})]", "priority": 6, "description": "IFF: O in bulk silica, Si-O-Si"},
    {"type": "oc24", "smarts": f"[OX2H1][Si;$({SILICA_SI})]", "priority": 6, "description": "IFF: silanol O on a silica surface"},
    {"type": "hoy", "smarts": f"[H][OX2][Si;$({SILICA_SI})]", "priority": 6, "description": "IFF: silanol H on a silica surface"},
    {"type": "na+", "smarts": "[Na]", "priority": 6, "description": "IFF: sodium ion (silica, clays, salts)"},
] + [{"type": m, "smarts": f"[{m};{UNBONDED}]", "priority": 6, "description": f"IFF: {m} metal"} for m in METALS]

for iff, base, ff, frc in (("iff-pcff", "pcff-frc", "INTERFACE (IFF 1.5, PCFF)", "pcff_interface_v1_5.frc"),
                           ("iff-cvff", "cvff", "INTERFACE (IFF 1.5, CVFF)", "cvff_interface_v1_5.frc")):
    B = json.load(open(os.path.join(DATA, "typing", base + ".typing.json")))
    lib = {t["name"] for t in json.load(open(os.path.join(DATA, "forcefields", iff + ".json")))["atom_types"]}
    assert max(r.get("priority", 0) for r in B["rules"]) < 6, "IFF's own rules must outrank the base rules"
    rules = [r for r in B["rules"] if r["type"] in lib]
    left = sorted({r["type"] for r in B["rules"]} - lib)
    for r in OWN:
        assert r["type"] in lib, (iff, r["type"])
    doc = {"format": "caps-typing", "version": 1, "forcefield": ff,
           "description": f"CAPS's {base} rules for organic atoms (bench/typing/make_iff_rules.py), then IFF's silica, sodium "
                          f"and metal rules (priority 6) for {frc}. Among matching rules the highest priority wins."
                          + (f" Base types {frc} does not define, left out: {', '.join(left)}." if left else ""),
           "rules": rules + OWN, "unbonded_types": METALS}
    for k in ("pairs", "analogies", "analogy_source", "ordered"):
        if k in B:
            doc[k] = B[k]
    json.dump(doc, open(os.path.join(OUT, iff + ".typing.json"), "w"), indent=1, ensure_ascii=False)
    print(f"{iff}: {len(rules)} rules from {base} + {len(OWN)} IFF rules")
