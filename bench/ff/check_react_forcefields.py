#!/usr/bin/env python3
"""The crosslinking and blend engines (unchanged) with each polymer force field CAPS builds from its original source:
every case is typed by the force field, crosslinked or blended, re-typed after every reaction cycle, and must stay
complete (no untyped atom, no missing parameter) with finite energies.

  sulfur cure of cis-1,4-polyisoprene (natural rubber; sulfur_allylic), epoxide-amine (epoxy_amine_primary),
  ENR + acid (enr_acid_ester), a PS / cis-PI blend

usage: check_react_forcefields.py [--ff pcff,compass ...]   (needs the caps Python package: PYTHONPATH, CAPS_LIB)
"""
import math, sys
import caps

FFS = sys.argv[sys.argv.index("--ff") + 1].split(",") if "--ff" in sys.argv else ["pcff", "cvff", "compass", "oplsaa2024", "gaff2", "dreiding"]
NR = "[*]C/C=C(C)\\C[*]"


def finite(doc):
    e = doc.energy()
    return all(math.isfinite(v) for v in e.values() if isinstance(v, (int, float))), e.get("total", e.get("potential", float("nan")))


def case(ff, label, make, templates, react_kw):
    try:
        doc = make()
        rep = doc.field.assign(ff)
        if not rep["complete"]:
            return f"{label:28s} typing incomplete before the reaction: {len([a for a in rep['atoms'] if not a['type']])} untyped, missing {rep['missing'][:3]}"
        if templates:
            doc.react(templates, **react_kw)
            s = doc.react_summary()
            rep2 = doc.field.report() if hasattr(doc.field, "report") else doc.field.assign(ff)
            ok, e = finite(doc)
            n = s.get("reactions", s.get("bonds_formed", "?"))
            if not rep2["complete"]:
                return f"{label:28s} INCOMPLETE after {n} reactions: missing {rep2['missing'][:4]}"
            return f"{label:28s} ok · {n} reactions · {doc.summary()['atoms']} atoms · energy {'finite' if ok else 'NOT finite'} ({e:.1f} kcal/mol)"
        ok, e = finite(doc)
        return f"{label:28s} ok · {doc.summary()['atoms']} atoms · energy {'finite' if ok else 'NOT finite'} ({e:.1f} kcal/mol)"
    except Exception as x:
        return f"{label:28s} FAILED: {str(x)[:200]}"


small = dict(cycles=10, per_cycle=3, capture=8.0, relax_iterations=200, seed=2)


def rubber_with_sulfur():
    d = caps.polymer(NR, dp=10, chains=4, density=0.6, seed=1)
    d.insert("SS", 8, seed=3)   # H-S-S-H donors, as React.SulfurCureOfNaturalRubber
    return d


for ff in FFS:
    print(f"== {ff}")
    print("  ", case(ff, "sulfur cure, natural rubber", rubber_with_sulfur, "sulfur_allylic", small))
    print("  ", case(ff, "epoxide-amine", lambda: caps.pack(molecules=[("C1CO1", 20), ("NCCN", 10)], box=16, seed=1), "epoxy_amine_primary", small))
    print("  ", case(ff, "ENR + acid", lambda: caps.pack(molecules=[("CC1(C)OC1CC", 10), ("CC(=O)O", 10)], box=16, seed=1), "enr_acid_ester", small))
    print("  ", case(ff, "PS / cis-PI blend", lambda: caps.blend([{"smiles": "*CC(*)c1ccccc1", "dp": 5, "chains": 2}, {"smiles": NR, "dp": 5, "chains": 2}],
                                                                 density=0.5, seed=1), None, None))
