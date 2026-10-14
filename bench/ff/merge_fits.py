#!/usr/bin/env python3
"""CAPS's QM fits into the OPLS-AA 2024 supplement (bench/ff/open/oplsaa2024_supplement.json), each term labelled with
where it comes from, then the library is rebuilt with bench/ff/convert_boss_opls.py --supplement.

    merge_fits.py SUPPLEMENT.json torsions FITTED_TORSIONS.json [...]      (caps.torsionfit fit output)
    merge_fits.py SUPPLEMENT.json angle FIT.json [...]                     (bench/ff/fit_angle_scan.py --json output)
    merge_fits.py SUPPLEMENT.json bond FIT.json [...]                      (bench/ff/fit_angle_scan.py --bond --json output)
    merge_fits.py UA_SUPPLEMENT.json ua-torsions FITTED_TORSIONS.json [...]  (united-atom fits into oplsua2024_supplement.json)

A term already in the supplement from a paper is kept (a paper's value wins over a CAPS fit); an earlier CAPS fit of
the same class is replaced. Refits of classes the force field has ("CAPS refit") are not merged: they stay overlays.
"""
import json, re, sys


def key(m):
    m = [x if isinstance(x, str) else "|".join(x) for x in m]   # a supplement match may list alternatives
    return tuple(m) if tuple(m) <= tuple(m[::-1]) else tuple(m[::-1])


def main():
    sup_path, kind, files = sys.argv[1], sys.argv[2], sys.argv[3:]
    S = json.load(open(sup_path))
    n = 0
    if kind == "torsions":
        lst = S.setdefault("dihedrals", [])
        for f in files:
            for r in json.load(open(f)).get("dihedrals", []):
                if r["comment"].startswith("CAPS refit"):
                    continue
                m = [x.replace("\\*", "*") for x in r["match"]]
                old = next((t for t in lst if key(t["match"]) == key(m)), None)
                if old and not old["source"].startswith("CAPS fit"):
                    continue
                if old:
                    lst.remove(old)
                lst.append({"match": m, "params": r["params"], "source": "CAPS fit (QM)",
                            "reason": r["comment"].replace("CAPS fit (QM): ", "") + f"; fragment {r['fit']['fragment']}"})
                n += 1
    elif kind == "angle":
        lst = S.setdefault("angles", [])
        for f in files:
            d = json.load(open(f))
            m = d["class"].split("-")
            km = tuple(m) if m[0] <= m[2] else tuple(m[::-1])
            old = next((t for t in lst if (tuple(t["match"]) if t["match"][0] <= t["match"][2] else tuple(t["match"][::-1])) == km), None)
            if old and not old["source"].startswith("CAPS fit"):
                continue
            if old:
                lst.remove(old)
            lst.append({"match": m, "params": [d["K"], d["theta0"]], "source": "CAPS fit (QM)",
                        "reason": f"{d['level']} relaxed scan of {d['molecule']} ({d['points']} points), "
                                  f"E = K (theta - theta0)^2, RMS {d['rms']:.2f} kcal/mol"})
            n += 1
    elif kind == "bond":
        lst = S.setdefault("bonds", [])
        for f in files:
            d = json.load(open(f))
            m = d["class"].split("-")
            km = tuple(sorted(m))
            old = next((t for t in lst if tuple(sorted(t["match"])) == km), None)
            if old and not old["source"].startswith("CAPS fit"):
                continue
            if old:
                lst.remove(old)
            lst.append({"match": m, "params": [d["K"], d["r0"]], "source": "CAPS fit (QM)",
                        "reason": f"{d['level']} relaxed scan of the bond in {d['molecule']} ({d['points']} points), "
                                  f"E = K (r - r0)^2, RMS {d['rms']:.2f} kcal/mol"})
            n += 1
    elif kind == "ua-torsions":
        # united-atom fits (caps.torsionfit fit on a plan whose ff is oplsua2024.json) into oplsua2024_supplement.json:
        # each site type read as its class (_b<class>_), then the same function for every carbon of its kind (the
        # supplement's class groups S and D) as an analogy listed after the exact classes, which come first and win
        group = {"C3": "S", "C2": "S", "CH": "S", "CT": "S", "C7": "D", "C8": "D", "C9": "D"}
        lst = S.setdefault("dihedrals", [])
        exact, generic = [], []
        for f in files:
            for r in json.load(open(f)).get("dihedrals", []):
                if r["comment"].startswith("CAPS refit") or r.get("fit", {}).get("zeroed"):
                    continue
                m = [re.search(r"_b([^_]+)_", x).group(1) if "_b" in x else x for x in r["match"]]
                why = r["comment"].replace("CAPS fit (QM): ", "") + f"; all-atom QM scan {r['fit']['fragment']}, MM with OPLS-UA"
                exact.append({"match": m, "params": r["params"], "source": "CAPS fit (QM)", "reason": why, "sulfur_rubber": True})
                g = [group.get(x, x) for x in m]
                if g != m:
                    generic.append({"match": g, "params": r["params"], "source": "CAPS fit (QM)",
                                    "reason": f"by_analogy: the {'-'.join(m)} fit for any carbon of the same kind; " + why,
                                    "sulfur_rubber": True})
        for e in exact + generic:
            if any(key(t["match"]) == key(e["match"]) and not t["source"].startswith("CAPS fit") for t in lst):
                continue   # a paper's value wins
            old = [t for t in lst if key(t["match"]) == key(e["match"])]
            if any(o["reason"].startswith("by_analogy") is False for o in old) and e["reason"].startswith("by_analogy"):
                continue   # an exact fit already there wins over an analogy
            for o in old:
                lst.remove(o)
            lst.append(e)
            n += 1
    else:
        raise SystemExit(__doc__)
    json.dump(S, open(sup_path, "w"), indent=1, ensure_ascii=False)
    print(f"{n} {kind} merged into {sup_path}")


if __name__ == "__main__":
    main()
