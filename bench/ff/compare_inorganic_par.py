#!/usr/bin/env python3
"""A per-paper inorganic library against DL_FIELD's .par / .sf files, term by term (DL_FIELD as an accuracy reference only:
nothing is copied from it). For each library type mapped to a DL_FIELD key (bench/typing/reference_type_maps.json): the
charge the reference's templates give that key, every short-range pair (Buckingham A, rho, C; a Rydberg pair against
whatever the reference writes for those keys), every core-shell spring (k2, k4) and Lennard-Jones pair, in eV. Prints the
terms that agree, those that differ (library value | reference value) and those one side lacks.

usage: compare_inorganic_par.py FF_ID [FF_ID ...]
Environment: DLFIELD (default ~/project/dl_f_4.13).
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
LIB = os.environ.get("DLFIELD", os.path.expanduser("~/project/dl_f_4.13")) + "/lib"
EV = 23.060548   # the libraries' factor (kcal/mol per eV)
MAPS = json.load(open(os.path.join(ROOT, "bench", "typing", "reference_type_maps.json")))


def read_reference(lib):
    """(charges by key, springs {(k1, k2): (k2, k4)}, pairs {(k1, k2): (form, a, b, c)}) from LIB.sf and LIB.par"""
    key_of, charges = {}, {}
    sec = None
    for line in open(os.path.join(LIB, lib + ".sf"), errors="replace"):
        w = line.split()
        if not w or w[0].startswith("#"):
            continue
        h = w[0].upper()
        if h == "ATOM_TYPE":
            sec = "types"; continue
        if h == "MOLECULE" and len(w) >= 2:
            sec = "mol"; continue
        if h == "END":
            sec = None; continue
        if sec == "types" and len(w) >= 3:
            key_of.setdefault(w[0], w[1])
        elif sec == "mol" and len(w) >= 3 and w[1] in key_of:
            try:
                charges.setdefault(key_of[w[1]], set()).add(float(w[2]))
            except ValueError:
                pass
    springs, pairs, sec = {}, {}, None
    for line in open(os.path.join(LIB, lib + ".par"), errors="replace"):
        w = line.split("!")[0].split()
        if not w or w[0].startswith("#"):
            continue
        h = w[0].upper()
        if h in ("SHELL", "VDW"):
            sec = h; continue
        if h == "END":
            sec = None; continue
        try:
            if sec == "SHELL" and len(w) >= 4:
                springs[frozenset((w[0], w[1]))] = (float(w[2]), float(w[3]))
            elif sec == "VDW" and len(w) >= 5:
                pairs[frozenset((w[0], w[1]))] = (w[2], float(w[3]), float(w[4]), float(w[5]) if len(w) > 5 and re.match(r"^-?[\d.]", w[5]) else 0.0)
        except ValueError:
            pass
    return charges, springs, pairs


def close(x, y):
    return abs(x - y) <= 1e-6 * max(1.0, abs(x), abs(y)) or round(x, 4) == round(y, 4)


def compare(fid):
    m = MAPS[fid]
    tmap, lib = m["map"], m["lib"]
    ff = json.load(open(os.path.join(ROOT, "data", "forcefields", fid + ".json")))
    q, springs, pairs = read_reference(lib)
    same, diff, missing = [], [], []
    for t in ff["atom_types"]:
        k = tmap.get(t["name"])
        if not k:
            continue
        ref = q.get(k)
        if ref is None:
            missing.append(f"charge {t['name']} ({k}): no reference template")
        elif any(close(t["charge"], r) for r in ref):
            same.append(f"charge {t['name']} {t['charge']}")
        else:
            diff.append(f"charge {t['name']} ({k}): {t['charge']} | {sorted(ref)}")
    for b in ff["bonds"]:
        a, c = b["match"]
        if a not in tmap or c not in tmap:
            continue
        lib2, lib4 = 2 * b["params"][1] / EV, 24 * b["params"][3] / EV
        r = springs.get(frozenset((tmap[a], tmap[c])))
        if r is None:
            missing.append(f"spring {a}-{c}: the reference has none ({tmap[a]}-{tmap[c]})")
        elif close(lib2, r[0]) and close(lib4, r[1]):
            same.append(f"spring {a}-{c} k2 {lib2:g} k4 {lib4:g}")
        else:
            diff.append(f"spring {a}-{c}: k2 {lib2:g} k4 {lib4:g} | k2 {r[0]:g} k4 {r[1]:g}")
    for p in ff["pairs"]:
        if len(p["match"]) != 2 or p["match"][0] not in tmap or p["match"][1] not in tmap:
            continue
        a, c = p["match"]
        st = p.get("style", "lj")
        if st == "lj" and not any(p["params"]):
            continue
        r = pairs.get(frozenset((tmap[a], tmap[c])))
        if st in ("buck", "rydberg"):
            mine = (p["params"][0] / EV, p["params"][1], p["params"][2] / EV)
        else:
            mine = (p["params"][0] / EV, p["params"][1])
        label = f"{st} {a}-{c}"
        if r is None:
            missing.append(f"{label} {tuple(round(x, 6) for x in mine)}: the reference has none ({tmap[a]}-{tmap[c]})")
            continue
        ref = (r[1], r[2], r[3]) if r[0] == "buckingham" else (r[1], r[2])
        if (st == "buck" and r[0] == "buckingham") or (st == "lj" and r[0] == "lj"):
            if all(close(x, y) for x, y in zip(mine, ref)):
                same.append(f"{label} {tuple(round(x, 6) for x in mine)}")
                continue
        diff.append(f"{label} {tuple(round(x, 6) for x in mine)} | {r[0]} {ref}")
    print(f"== {fid} against {lib}: {len(same)} terms agree, {len(diff)} differ, {len(missing)} on one side only")
    for x in diff:
        print("  differs:", x)
    for x in missing:
        print("  one side:", x)
    if "--all" in sys.argv:
        for x in same:
            print("  same:", x)


for fid in [a for a in sys.argv[1:] if not a.startswith("--")]:
    compare(fid)
