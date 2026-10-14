#!/usr/bin/env python3
"""Does a force field cover what CAPS React makes? Every built-in reaction template on the reactants it is meant for:
a small cell of the reactants is cured with the template (to half and to full conversion, so partly reacted
crosslinkers are in it too), the product typed with the force field, and untyped atoms, missing terms and non-neutral
molecules are listed. A force field that is "complete" on the polymer library can still lack the terms a reaction
creates (an aryl glycidyl ether opened by an amine, a maleate half-ester): this finds them before a user does.

    react_coverage.py [--ff oplsaa2024] [--only TEMPLATE_OR_CASE] [--json out.json] [--save DIR]

--save DIR writes each product that is not complete (DIR/<case>.data, with its cell) for python -m caps.autofit
(bench/ff/autofit_sweep.py runs it over all of them).

Exit 0 only when every case is complete and neutral. The missing terms are listed once each with the cases that need
them; the torsions among them go to QM with python -m caps.torsionfit prepare (on a product fragment).
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
os.environ.setdefault("CAPS_LIB", os.path.join(ROOT, "build", "capi", "libcaps.dylib"))
sys.path.insert(0, os.path.join(ROOT, "data", "python"))
import caps  # noqa: E402

arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FF = arg("--ff", "oplsaa2024")
ONLY = arg("--only", "")

# oligomer stand-ins for the chains (a few repeat units, ends capped with H), the curatives and co-reactants
NR = "C/C(C)=C\\CC/C(C)=C\\CC/C(C)=C\\C"                       # cis-1,4-polyisoprene
BR = "C/C=C\\CC/C=C\\CC/C=C\\C"                                  # cis-1,4-polybutadiene
SBR = "CC/C=C\\CC(c1ccccc1)C/C=C\\C"                             # styrene-butadiene
EPDM = "CCCC(C)CCCC1CC2CC1C(=CC)C2"                              # ethylene-propylene with an ENB unit
PE = "CCCCCCCCCCCC"
PP = "CC(C)CC(C)CC(C)CC(C)C"
ENR = "CC/C=C(C)\\CCC1(C)OC1CC/C=C(C)\\CC"                       # epoxidised NR (one epoxide)
DGEBA = "CC(C)(c1ccc(OCC2CO2)cc1)c1ccc(OCC2CO2)cc1"
BGE = "CCCCOCC1CO1"                                              # butyl glycidyl ether (reactive diluent)
TGDDM = "C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1"

CASES = [
    # template, case name, [(smiles, count)], box (A)
    ("cc_crosslink", "PE", [(PE, 6)], 22),
    ("cc_crosslink", "NR", [(NR, 6)], 22),
    ("cc_crosslink", "EPDM", [(EPDM, 6)], 22),
    ("cc_crosslink_any", "PP", [(PP, 6)], 22),
    ("cc_crosslink_any", "PE + PP", [(PE, 3), (PP, 3)], 22),
    ("sulfur_allylic", "NR + H2S2", [(NR, 4), ("SS", 8)], 22),
    ("sulfur_allylic", "BR + H2S2", [(BR, 4), ("SS", 8)], 22),
    ("sulfur_allylic", "SBR + H2S2", [(SBR, 4), ("SS", 8)], 22),
    ("polysulfide_allylic", "NR + TESPT", [(NR, 4), ("CCO[Si](OCC)(OCC)CCCSSSSCCC[Si](OCC)(OCC)OCC", 2)], 24),
    ("polysulfide_allylic", "BR + TESPD", [(BR, 4), ("CCO[Si](OCC)(OCC)CCCSSCCC[Si](OCC)(OCC)OCC", 2)], 24),
    ("peroxide_allylic", "NR", [(NR, 6)], 22),
    ("peroxide_allylic", "BR", [(BR, 6)], 22),
    ("peroxide_allylic", "SBR", [(SBR, 6)], 22),
    ("enr_acid_ester", "ENR + acetic acid", [(ENR, 4), ("CC(=O)O", 6)], 22),
    ("enr_acid_ester", "ENR + PBS COOH end", [(ENR, 4), ("OCCCCOC(=O)CCC(=O)O", 4)], 24),
    ("enr_acid_ester", "ENR + maleic acid", [(ENR, 4), ("OC(=O)/C=C\\C(=O)O", 4)], 22),
    ("enr_acid_ester", "ENR + adipic acid", [(ENR, 4), ("OC(=O)CCCCC(=O)O", 4)], 22),
    ("enr_acid_ester", "ENR + stearic acid", [(ENR, 4), ("CCCCCCCCCCCCCCCCCC(=O)O", 3)], 26),
    ("enr_acid_ester", "ENR + benzoic acid", [(ENR, 4), ("OC(=O)c1ccccc1", 4)], 22),
    ("ester_condensation", "succinic acid + butanediol (PBS)", [("OC(=O)CCC(=O)O", 4), ("OCCCCO", 4)], 20),
    ("ester_condensation", "adipic acid + butanediol (PBA)", [("OC(=O)CCCCC(=O)O", 4), ("OCCCCO", 4)], 20),
    ("ester_condensation", "terephthalic acid + ethylene glycol (PET)", [("OC(=O)c1ccc(C(=O)O)cc1", 4), ("OCCO", 4)], 20),
    ("ester_condensation", "lactic acid (PLA)", [("CC(O)C(=O)O", 8)], 20),
    ("ester_condensation", "glycerol + adipic acid", [("OCC(O)CO", 3), ("OC(=O)CCCCC(=O)O", 4)], 20),
    ("anhydride_alcohol", "MAH + butanediol", [("O=C1C=CC(=O)O1", 4), ("OCCCCO", 4)], 20),
    ("anhydride_alcohol", "succinic anhydride + butanol", [("O=C1CCC(=O)O1", 4), ("CCCCO", 4)], 20),
    ("anhydride_alcohol", "phthalic anhydride + glycerol", [("O=C1OC(=O)c2ccccc21", 4), ("OCC(O)CO", 3)], 20),
    ("anhydride_alcohol", "MHHPA + butanediol", [("CC1CCC2C(=O)OC(=O)C2C1", 4), ("OCCCCO", 4)], 22),
    ("anhydride_alcohol", "MAH + ENR-OH (opened epoxide)", [("O=C1C=CC(=O)O1", 3), ("CC/C=C(C)\\CCC(C)(O)C(O)CC/C=C(C)\\CC", 3)], 22),
    ("epoxy_amine_primary", "DGEBA + DETA", [(DGEBA, 4), ("NCCNCCN", 3)], 24),
    ("epoxy_amine_primary", "DGEBA + IPDA", [(DGEBA, 4), ("CC1(C)CC(N)CC(C)(CN)C1", 3)], 24),
    ("epoxy_amine_primary", "DGEBA + Jeffamine D230 (n=2)", [(DGEBA, 4), ("CC(N)COCC(C)OCC(C)N", 3)], 24),
    ("epoxy_amine_primary", "DGEBA + DDM (aromatic amine)", [(DGEBA, 4), ("Nc1ccc(Cc2ccc(N)cc2)cc1", 3)], 24),
    ("epoxy_amine_primary", "DGEBA + DDS (aromatic sulfone amine)", [(DGEBA, 4), ("Nc1ccc(S(=O)(=O)c2ccc(N)cc2)cc1", 3)], 24),
    ("epoxy_amine_primary", "butyl glycidyl ether + ethylenediamine", [(BGE, 6), ("NCCN", 3)], 20),
    ("epoxy_amine_primary", "TGDDM + DDS", [(TGDDM, 2), ("Nc1ccc(S(=O)(=O)c2ccc(N)cc2)cc1", 3)], 24),
    ("epoxy_amine_secondary", "DGEBA + DETA (after the primary step)", [(DGEBA, 4), ("NCCNCCN", 3)], 24),
    ("epoxy_amine_secondary", "DGEBA + piperidine-type secondary amine", [(DGEBA, 4), ("C1CCNCC1", 4)], 24),
]
# the secondary-amine step needs its primary step first: those cases run both templates
TWO_STEP = {"epoxy_amine_secondary": ["epoxy_amine_primary", "epoxy_amine_secondary"]}


def cure(template, mols, box, target):
    d = caps.pack(mols, box=float(box))
    names = TWO_STEP.get(template, [template]) if "after the primary" in str(mols) or template in TWO_STEP else [template]
    d.react(names, cycles=80, per_cycle=8, target=target, relax=False, default_field=True, auto_capture=True,
            capture_max=12.0, capture_step=1.0, keep_byproducts=False, seed=3)
    return d


def check(d):
    r = d.field.assign(FF)
    return {"complete": bool(r["complete"]), "untyped": int(r.get("untyped", 0)), "missing": list(r["missing"]),
            "net_charge": float(r.get("net_charge", 0.0))}


def main():
    results, gaps = [], {}
    for template, name, mols, box in CASES:
        if ONLY and ONLY not in template and ONLY not in name:
            continue
        for target in (0.5, 1.0):
            label = f"{template}: {name} @ {int(target * 100)} %"
            try:
                d = cure(template, mols, box, target)
                summary = d.react_summary()
                res = check(d)
            except Exception as e:   # the cure or the assignment failed: reported as such
                res, summary = {"complete": False, "error": str(e)[:300], "missing": [], "untyped": -1, "net_charge": 0.0}, {}
            res.update({"case": label, "template": template, "crosslinks": summary.get("crosslinks")})
            ok = res["complete"] and abs(res["net_charge"]) < 1e-6 and res["untyped"] == 0
            if not ok and arg("--save") and "error" not in res:
                os.makedirs(arg("--save"), exist_ok=True)
                slug = "".join(c if c.isalnum() else "_" for c in label).strip("_")
                res["structure"] = os.path.join(os.path.abspath(arg("--save")), slug + ".data")
                d.save(res["structure"])
            print(f"{'ok  ' if ok else 'GAP '} {label}: untyped {res['untyped']}, missing {len(res['missing'])}, "
                  f"net charge {res['net_charge']:+.3f}" + (f", crosslinks {res['crosslinks']}" if res.get("crosslinks") is not None else "")
                  + (f"  [{res['error']}]" if "error" in res else ""))
            for m in res["missing"]:
                gaps.setdefault(m, []).append(f"{template}: {name}")
            results.append(res)
    if gaps:
        print("\nmissing terms (cases):")
        for m, where in sorted(gaps.items()):
            print(f"  {m}   [{'; '.join(sorted(set(where)))}]")
    if arg("--json"):
        json.dump({"forcefield": FF, "cases": results, "missing": {m: sorted(set(w)) for m, w in gaps.items()}},
                  open(arg("--json"), "w"), indent=1)
    bad = [r for r in results if not (r["complete"] and abs(r["net_charge"]) < 1e-6 and r["untyped"] == 0)]
    print(f"\n{len(results) - len(bad)} of {len(results)} cases complete and neutral with {FF}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
