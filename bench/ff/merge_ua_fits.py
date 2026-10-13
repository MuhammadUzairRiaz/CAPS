#!/usr/bin/env python3
"""CAPS's QM torsion fits for the united-atom force fields into the files their generators read.

    merge_ua_fits.py trappe FITTED.json [...]   -> bench/ff/open/trappe_ua_capsfits.json (exact TraPPE type names)
    merge_ua_fits.py oplsua FITTED.json [...]   -> bench/ff/open/oplsua2024_supplement.json (by dihedral class)

FITTED.json is caps.torsionfit fit output. A later fit of the same torsion replaces an earlier CAPS fit; refits of
terms the force field has ("CAPS refit") are not merged. Then regenerate the library (make_trappe.py SITE_DIR /
make_oplsua2024.py).
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def key(m):
    m = tuple(x if isinstance(x, str) else "|".join(x) for x in m)   # a supplement position may list several classes
    return m if m <= m[::-1] else m[::-1]


def ua_class(name):
    """the dihedral class of an OPLS-UA 2024 type name such as 108_bOS_aOS_dOS_iOS"""
    m = re.search(r"_d([^_]+)_i", name)
    return m.group(1) if m else name


def main():
    kind, files = sys.argv[1], sys.argv[2:]
    path = os.path.join(HERE, "open", "trappe_ua_capsfits.json" if kind == "trappe" else "oplsua2024_supplement.json")
    S = json.load(open(path))
    lst = S.setdefault("dihedrals", [])
    n = 0
    for f in files:
        for r in json.load(open(f)).get("dihedrals", []):
            if r["comment"].startswith("CAPS refit"):
                continue
            m = [x.replace("\\*", "*") for x in r["match"]]
            if kind == "oplsua":
                m = [ua_class(x) for x in m]
                entry = {"match": m, "params": r["params"], "source": "CAPS fit (QM)", "reason": r["comment"]}
                old = next((t for t in lst if key(t["match"]) == key(m)), None)
                if old and old.get("source") != "CAPS fit (QM)":
                    continue   # a paper's term wins
            else:
                entry = {"match": m, "params": r["params"], "comment": r["comment"]}
                old = next((t for t in lst if key(t["match"]) == key(m)), None)
            if old:
                lst.remove(old)
            lst.append(entry)
            n += 1
    json.dump(S, open(path, "w"), indent=1, ensure_ascii=False)
    print(f"{n} torsions merged into {path}")


if __name__ == "__main__":
    main()
