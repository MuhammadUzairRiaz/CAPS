#!/usr/bin/env python3
"""OPLS-AA by type number on the polymer library: every oligomer typed with data/typing/oplsaa2024-moltemplate.typing.json
and parameterised with the OPLS-AA 2024 table; each molecule must be fully typed, have every parameter and be neutral
(OPLS-AA's group charges balance). Where a reference directory holds another program's OPLS charges for the same xyz
(bench/ff/compare_polymers.py --keep DIR), the charges are compared atom by atom.

usage: check_oplsaa_numbers.py MOLDIR [--ref DIR] [--ff oplsaa2024-moltemplate]
"""
import glob, json, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
MOLDIR = sys.argv[1]
REF = arg("--ref", "")
FF = arg("--ff", "oplsaa2024-moltemplate")
ffj = os.path.join(ROOT, "data", "forcefields", FF + ".json")
rules = os.path.join(ROOT, "data", "typing", FF + ".typing.json")
names = {p["id"]: p["name"] for p in json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))["polymers"]}


def charges(path):
    q, sec = [], None
    for raw in open(path):
        l = raw.split("#")[0].strip()
        if re.match(r"^[A-Z][A-Za-z ]+$", l):
            sec = l
            continue
        w = l.split()
        if sec == "Atoms" and len(w) >= 7:
            q.append((int(w[0]), float(w[3])))
    return [c for _, c in sorted(q)]


ok = untyped = missing = charged = 0
rows = []
tmp = os.path.join(MOLDIR, "..", "opls_numbers")
os.makedirs(tmp, exist_ok=True)
for xyz in sorted(glob.glob(os.path.join(MOLDIR, "*.xyz"))):
    if xyz.endswith(".ext.xyz"):
        continue
    pid = os.path.basename(xyz)[:-4]
    out = os.path.join(tmp, pid + ".data")
    r = subprocess.run([CAPS, "ff", "apply", xyz, "--ff", ffj, "--typing", rules, "--charges", "types", "-o", out], capture_output=True, text=True)
    txt = r.stdout + r.stderr
    if r.returncode or not os.path.exists(out):
        m = re.search(r"(\d+) atoms? (?:match no|are untyped|untyped)[^\n]*", txt)
        why = (m.group(0) if m else (txt.strip().splitlines() or ["?"])[-1])[:160]
        if "match no typing rule" not in txt and ("missing" in txt or "no parameters" in txt or "dihedral" in why or "angle" in why or "bond" in why):
            missing += 1
            kind = "missing"
        else:
            untyped += 1
            kind = "untyped"
        rows.append((pid, kind, why))
        if os.path.exists(out):
            os.remove(out)
        continue
    q = charges(out)
    net = sum(q)
    ref = ""
    rp = os.path.join(REF, f"{pid}__opls2005", "dlfield", "lammps1.data") if REF else ""
    if rp and os.path.exists(rp):
        rq = charges(rp)
        if len(rq) == len(q):
            ref = f" · largest difference from the reference charges {max(abs(a - b) for a, b in zip(q, rq)):.3f} e"
    if abs(net) > 1e-3:
        charged += 1
        rows.append((pid, "charged", f"net charge {net:+.3f} e{ref}"))
    else:
        ok += 1
        rows.append((pid, "ok", f"neutral{ref}"))
for pid, kind, why in rows:
    print(f"{pid} {names.get(pid, '')[:32]:32s} {kind:8s} {why}")
print(f"\n{ok} neutral and complete · {charged} typed but not neutral · {untyped} with untyped atoms · {missing} with missing parameters (of {len(rows)})")
