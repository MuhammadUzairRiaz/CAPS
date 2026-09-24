#!/usr/bin/env python3
"""Check CAPS's automatic typing against the atom types in a CHARMM topology file (.rtf), e.g. CGenFF's
top_all36_cgenff.rtf, whose model-compound residues were typed by the force field's developers.

Each RESI becomes a mol2 file with its bonds but no bond orders (or, with --orders, the orders the RTF gives with
DOUBLE / TRIPLE); CAPS types it and the result is compared with the RTF types. Residues with lone pairs or dummy
atoms, and residues using types missing from the force-field file, are counted separately.

usage: validate_rtf_types.py TOP.rtf FF.json RULES.json [--orders] [--show N] [--all] [--only RESI,...] [--apply]
"""
import collections, json, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
rtf, ffjson, rules = sys.argv[1], sys.argv[2], sys.argv[3]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
show = int(arg("--show", "40"))
only = [x for x in arg("--only", "").split(",") if x]
known = {t["name"] for t in json.load(open(ffjson))["atom_types"]}

ELEM = [("CLG", "Cl"), ("BRG", "Br"), ("ALG", "Al"), ("IG", "I"), ("H", "H"), ("C", "C"), ("N", "N"), ("O", "O"),
        ("S", "S"), ("P", "P"), ("F", "F"), ("B", "B")]


def element(t):
    for p, e in ELEM:
        if t.upper().startswith(p):
            return e
    return None


def residues(path):
    res, cur = [], None
    for line in open(path, errors="replace"):
        s = line.split("!")[0].strip()
        w = s.split()
        if not w:
            continue
        k = w[0].upper()
        if k in ("RESI", "PRES"):
            cur = {"name": w[1], "patch": k == "PRES", "atoms": [], "bonds": []}
            res.append(cur)
        elif cur is None:
            continue
        elif k == "ATOM" and len(w) >= 3:
            cur["atoms"].append((w[1], w[2]))
        elif k in ("BOND", "DOUB", "DOUBLE", "TRIP", "TRIPLE", "AROM"):
            order = {"BOND": 1, "DOUB": 2, "DOUBLE": 2, "TRIP": 3, "TRIPLE": 3, "AROM": 4}[k]
            for a, b in zip(w[1::2], w[2::2]):
                cur["bonds"].append((a, b, order))
        elif k == "END":
            break
    return [r for r in res if not r["patch"]]


work = tempfile.mkdtemp()
total = right = 0
confusion, wrong, skipped, failed = collections.Counter(), [], collections.Counter(), []
applied, unparam = [0], []
for r in residues(rtf):
    if only and r["name"] not in only:
        continue
    names = [a for a, _ in r["atoms"]]
    types = [t for _, t in r["atoms"]]
    if any(element(t) is None or t.upper().startswith("LP") or t.upper().startswith("DUM") for t in types):
        skipped["lone pairs / dummy atoms"] += 1
        continue
    if any(t not in known for t in types):
        skipped["types not in the force-field file"] += 1
        continue
    idx = {a.upper(): i for i, a in enumerate(names)}
    bonds = {}
    ok = True
    for a, b, o in r["bonds"]:
        a, b = a.upper().lstrip("+-"), b.upper().lstrip("+-")
        if a not in idx or b not in idx:
            ok = False
            break
        key = (min(idx[a], idx[b]), max(idx[a], idx[b]))
        bonds[key] = max(bonds.get(key, 0), o)
    if not ok or not bonds and len(names) > 1:
        skipped["bonds to other residues"] += 1
        continue
    path = os.path.join(work, re.sub(r"[^A-Za-z0-9_]", "_", r["name"]) + ".mol2")
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{r['name']}\n{len(names)} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, (a, t) in enumerate(r["atoms"]):
            f.write(f"{i + 1} {a} {i * 1.5:.3f} 0 0 {element(t)} 1 {r['name'][:4]} 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for k, ((i, j), o) in enumerate(sorted(bonds.items())):
            code = {1: "1", 2: "2", 3: "3", 4: "ar"}[o] if "--orders" in sys.argv else "un"
            f.write(f"{k + 1} {i + 1} {j + 1} {code}\n")
    out = path + ".types"
    p = subprocess.run([CAPS, "ff", "type", path, "--ff", ffjson, "--typing", rules, "-o", out], capture_output=True, text=True)
    if not os.path.exists(out):
        failed.append((r["name"], (p.stderr or p.stdout).strip()[-120:]))
        continue
    got = [x.strip() for x in open(out)]
    if "--apply" in sys.argv:   # end to end: parameters for the automatically typed residue
        pa = subprocess.run([CAPS, "ff", "apply", path, "--ff", ffjson, "--typing", rules, "--charges", "gasteiger", "--allow-missing"],
                            capture_output=True, text=True)
        applied[0] += 1
        if pa.returncode or "missing parameters" in pa.stdout:
            miss = re.findall(r"^  (.+)$", pa.stdout.split("missing parameters", 1)[1], re.M) if "missing parameters" in pa.stdout else []
            unparam.append((r["name"], miss[:3] or [(pa.stderr.strip().splitlines() or ["?"])[-1][:80]]))
    bad = [(names[i], types[i], got[i]) for i in range(len(types)) if types[i] != got[i]]
    total += len(types)
    right += len(types) - len(bad)
    for _, e, g in bad:
        confusion[(e, g)] += 1
    if bad:
        wrong.append((r["name"], bad))

n = total and (right / total * 100)
tested = sum(1 for _ in residues(rtf)) - sum(skipped.values()) - len(failed)
print(f"{os.path.basename(rtf)}: {right}/{total} atoms ({n:.1f} %) as the RTF, "
      f"{tested - len(wrong)} of {tested} residues fully right" + ("" if not skipped else " · skipped: " +
      ", ".join(f"{k} {v}" for k, v in skipped.items())) + (f" · {len(failed)} failed" if failed else ""))
for (e, g), k in confusion.most_common(show):
    print(f"  {e:>8} → {g:<8} {k}")
for name, why in failed[:5]:
    print(f"  failed {name}: {why}")
if "--apply" in sys.argv:
    print(f"end to end: {applied[0] - len(unparam)} of {applied[0]} residues fully parameterised with automatic types")
    kinds = collections.Counter(m.split()[0] for _, ms in unparam for m in ms)
    print("  missing terms by kind (first 3 per residue): " + ", ".join(f"{k} {v}" for k, v in kinds.most_common()))
if "--all" in sys.argv:
    for name, bad in wrong:
        print(f"--- {name}: " + ", ".join(f"{a} {e}→{g}" for a, e, g in bad))
