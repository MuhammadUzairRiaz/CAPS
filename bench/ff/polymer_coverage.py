#!/usr/bin/env python3
"""How much of the polymer library a force field covers: for every polymer of data/polymers/library.json, one DP-4 chain
grown by CAPS (caps grow --chains 1 --dp 4 --density 0.01 --seed 3; poly(p-phenylene vinylene) --seed 5) and the force
field applied (caps ff apply). Each polymer is complete (typed, every term found) and neutral, typed with missing terms
(the terms are counted by type class), or left untyped.

usage: polymer_coverage.py FF.json [--out DIR] [--only P001,P004] [--json result.json]
"""
import collections, json, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.environ.get("CAPS_BIN", os.path.join(ROOT, "build", "cli", "caps"))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
ff = sys.argv[1]
out = arg("--out") or tempfile.mkdtemp()
only = set(arg("--only", "").split(",")) - {""}
os.makedirs(out, exist_ok=True)


def charge(data):
    lines = open(data).read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.split("#")[0].strip() == "Atoms")
    q = 0.0
    for l in lines[i + 2:]:
        if not l.strip():
            break
        q += float(l.split()[3])
    return q


res, missing, who = {}, collections.Counter(), collections.defaultdict(set)
for p in json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))["polymers"]:
    if only and p["id"] not in only:
        continue
    xyz, data = os.path.join(out, p["id"] + ".xyz"), os.path.join(out, p["id"] + ".data")
    seed = "5" if "phenylene vinylene" in p["name"].lower() else "3"
    g = subprocess.run([CAPS, "grow", "--units", p["smiles"], "--chains", "1", "--dp", "4", "--density", "0.01", "--seed", seed,
                        "-o", xyz], capture_output=True, text=True)
    if g.returncode:
        res[p["id"]] = {"name": p["name"], "state": "not grown"}
        continue
    a = subprocess.run([CAPS, "ff", "apply", xyz, "--ff", ff, "-o", data], capture_output=True, text=True)
    txt = a.stdout + a.stderr
    if a.returncode == 0:
        q = charge(data)
        res[p["id"]] = {"name": p["name"], "state": "complete" if abs(q) < 1e-6 else "charged", "charge": q}
    elif "match no typing rule" in txt:
        res[p["id"]] = {"name": p["name"], "state": "untyped"}
    else:
        terms = []
        for kind, types in re.findall(r"(bond|angle|dihedral|improper) ((?:\S+_b\S+_a\S+_d\S+_i\S+ ?)+)", txt):
            cl = [re.search("_" + kind[0] + r"([^_]+)", t).group(1) for t in types.split()]
            key = kind + " " + "-".join(min(cl, cl[::-1]))
            terms.append(key)
            missing[key] += 1
            who[key].add(p["id"])
        res[p["id"]] = {"name": p["name"], "state": "missing terms", "missing": sorted(set(terms))}
states = collections.Counter(r["state"] for r in res.values())
for pid, r in res.items():
    print(f"{pid:5s} {r['state']:14s} {r['name']}")
print("\n" + ", ".join(f"{v} {k}" for k, v in states.most_common()) + f" of {len(res)}")
for k, v in missing.most_common():
    print(f"  {v:3d} {k}  ({', '.join(sorted(who[k]))})")
if arg("--json"):
    json.dump(res, open(arg("--json"), "w"), indent=1)
