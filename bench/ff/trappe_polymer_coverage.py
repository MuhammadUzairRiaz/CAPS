#!/usr/bin/env python3
"""TraPPE-UA on the polymer library: which repeat units CAPS can parameterise whole.

  trappe_polymer_coverage.py [--ff LIBRARY] [--json OUT] [--keep DIR] [--polymer SUBSTRING]

For each repeat unit of data/polymers/library.json CAPS grows one isolated chain (DP 4, hydrogen ends;
caps grow --units SMILES --chains 1 --dp 4 --density 0.01 --seed 3, seed 5 for poly(p-phenylene vinylene)) and
applies the library (default data/forcefields/trappe-ua-umn.json). A polymer counts as complete when every atom is
typed, no bonded term is missing and the net charge is zero; otherwise the first lines of what is missing are kept.
"""
import json, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FF = arg("--ff", os.path.join(ROOT, "data", "forcefields", "trappe-ua-umn.json"))
WORK = arg("--keep", "") or tempfile.mkdtemp(prefix="caps_trappe_poly_")
ONLY = arg("--polymer", "").lower()
os.makedirs(WORK, exist_ok=True)


def charge_of(path):
    """net charge of a LAMMPS data file (atom style full: id mol type q x y z)"""
    q, on = 0.0, False
    for l in open(path):
        s = l.split("#")[0].split()
        if l.startswith("Atoms"):
            on = True
            continue
        if on and l.strip() and not s:
            continue
        if on and s and not s[0].lstrip("-").isdigit():
            break
        if on and len(s) >= 7:
            q += float(s[3])
    return q


lib = json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))["polymers"]
out = []
for p in lib:
    if ONLY and ONLY not in p["name"].lower():
        continue
    stem = os.path.join(WORK, p["id"])
    seed = "5" if "phenylene vinylene" in p["name"].lower() else "3"
    g = subprocess.run([CAPS, "grow", "--units", p["smiles"], "--chains", "1", "--dp", "4", "--density", "0.01", "--seed", seed,
                        "-o", stem + ".xyz"], capture_output=True, text=True)
    if not os.path.exists(stem + ".xyz"):
        out.append({"id": p["id"], "name": p["name"], "complete": False, "why": "grow failed: " + (g.stderr or g.stdout).strip()[-160:]})
        continue
    r = subprocess.run([CAPS, "ff", "apply", stem + ".xyz", "--ff", FF, "-o", stem + ".data"], capture_output=True, text=True)
    text = r.stdout + r.stderr
    ok = r.returncode == 0 and os.path.exists(stem + ".data")
    q = charge_of(stem + ".data") if ok else None
    neutral = q is not None and abs(q) < 1e-6
    why = ""
    if not ok:
        lines = text.splitlines()
        miss = []
        for l in lines:
            m = re.search(r"(\d+ atoms match no typing rule: .*?) \(", l)
            if m:
                miss.append(m.group(1))
            m = re.match(r"\s+(bond|angle|dihedral|improper) (.*)", l)
            if m and l.strip() not in miss:
                miss.append(l.strip())
        why = "; ".join(dict.fromkeys(miss)) or text.strip()[-200:]
    elif not neutral:
        why = f"net charge {q:+.4f}"
    out.append({"id": p["id"], "name": p["name"], "complete": ok and neutral, "charge": q, "why": why})
done = [x for x in out if x["complete"]]
print(f"{len(done)} / {len(out)} polymers complete and neutral with {os.path.basename(FF)}")
for x in out:
    print(f"  {'ok  ' if x['complete'] else '--  '}{x['id']} {x['name']}" + ("" if x["complete"] else f": {x['why'][:300]}"))
if "--json" in sys.argv:
    json.dump(out, open(arg("--json", ""), "w"), indent=1)
