#!/usr/bin/env python3
"""CAPS's QM fits into the OPLS-AA 2024 supplement (bench/ff/open/oplsaa2024_supplement.json), each term labelled with
where it comes from, then the library is rebuilt with bench/ff/convert_boss_opls.py --supplement.

    merge_fits.py SUPPLEMENT.json torsions FITTED_TORSIONS.json [...]      (caps.torsionfit fit output)
    merge_fits.py SUPPLEMENT.json angle FIT.json [...]                     (bench/ff/fit_angle_scan.py --json output)

A term already in the supplement from a paper is kept (a paper's value wins over a CAPS fit); an earlier CAPS fit of
the same class is replaced. Refits of classes the force field has ("CAPS refit") are not merged: they stay overlays.
"""
import json, sys


def key(m):
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
    else:
        raise SystemExit(__doc__)
    json.dump(S, open(sup_path, "w"), indent=1, ensure_ascii=False)
    print(f"{n} {kind} merged into {sup_path}")


if __name__ == "__main__":
    main()
