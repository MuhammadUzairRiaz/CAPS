#!/usr/bin/env python3
"""Check CAPS's automatic typing against the molecule files of a force field in fftool's format (CL&P's ions: .zmat
z-matrices and .xyz files whose atom names are the force-field types).

A .zmat gives each atom its type and its first reference atom (the bond), and "connect i j" closes rings; a .xyz gives
positions (bonds from CAPS's perception). Elements come from the force field's types. CAPS types each molecule and the
result is compared atom by atom.

usage: validate_fftool_types.py FF.json RULES.json FILE.zmat|FILE.xyz ... [--all] [--show N] [--deviations FILE] [--list]

Each molecule's labels, with the library's type charges, must add up to an integer charge (reported otherwise); --list
names the molecules tested (a file whose labels the library does not all define is skipped).

--deviations: a JSON file {"molecules": {name: evidence}} of molecule files whose own labels are wrong, reported apart.
"""
import collections, json, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
ffjson, rules = sys.argv[1], sys.argv[2]
files = [a for a in sys.argv[3:] if a.endswith((".zmat", ".xyz"))]
deviations = json.load(open(sys.argv[sys.argv.index("--deviations") + 1]))["molecules"] if "--deviations" in sys.argv else {}
deviating = []
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
lib = json.load(open(ffjson))
element = {t["name"]: t["element"] for t in lib["atom_types"]}
# two names with the same class, charge and Lennard-Jones are the same atom to the force field (CL&P's CT2 and CS)
lj = {p_["name"]: tuple(p_["params"]) for p_ in lib["pairs"]}
by_name = {t["name"]: t for t in lib["atom_types"]}


def same(a, b):
    if a not in by_name or b not in by_name:
        return False
    ta, tb = by_name[a], by_name[b]
    return ta.get("charge") == tb.get("charge") and ta.get("equivalence") == tb.get("equivalence") and lj.get(a) == lj.get(b)


equivalent = collections.Counter()
work = os.environ.get("FFT_KEEP") or tempfile.mkdtemp()


def read_zmat(path):
    """(types, bonds) of a z-matrix: the first line names the molecule; an atom line is "[n] TYPE [ref dist [ref angle
    [ref dihedral]]]" (the number column is optional); "connect i j" adds a bond"""
    lines = [l.split("#")[0].strip() for l in open(path)]
    lines = [l for l in lines if l]
    atoms, bonds = [], []
    for l in lines[1:]:
        w = l.split()
        if w[0] == "connect" and len(w) >= 3:
            bonds.append((int(w[1]) - 1, int(w[2]) - 1))
            continue
        if w[0] in ("improper", "reconnect") or w[0].endswith((".ff", ".zmat")):
            continue
        if w[0].isdigit():
            w = w[1:]
        if not w or not re.match(r"^[A-Za-z]", w[0]) or (len(w) > 1 and not w[1].isdigit()):
            continue
        atoms.append(w[0])
        if len(w) >= 3:
            bonds.append((int(w[1]) - 1, len(atoms) - 1))
    return atoms, bonds


total = right = 0
confusion, wrong, failed, skipped = collections.Counter(), [], [], collections.Counter()
tested_names, fractional = [], []
for f in files:
    name = os.path.basename(f).rsplit(".", 1)[0]
    if f.endswith(".zmat"):
        types, bonds = read_zmat(f)
        if any(t not in element for t in types):
            skipped["types not in the force-field file"] += 1
            continue
        path = os.path.join(work, name + ".mol2")
        with open(path, "w") as o:
            o.write(f"@<TRIPOS>MOLECULE\n{name}\n{len(types)} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
            for i, t in enumerate(types):
                o.write(f"{i + 1} {t}{i + 1} {i * 1.5:.3f} 0 0 {element[t]} 1 MOL 0.0\n")
            o.write("@<TRIPOS>BOND\n")
            for k, (i, j) in enumerate(bonds):
                o.write(f"{k + 1} {i + 1} {j + 1} un\n")
    else:
        L = open(f).read().splitlines()
        n = int(L[0].split()[0])
        rows = [l.split() for l in L[2:2 + n]]
        types = [r[0] for r in rows]
        if any(t not in element for t in types):
            skipped["types not in the force-field file"] += 1
            continue
        path = os.path.join(work, name + ".xyz")
        open(path, "w").write(f"{n}\n{name}\n" + "".join(f"{element[r[0]]} {r[1]} {r[2]} {r[3]}\n" for r in rows))
    out = path + ".types"
    p = subprocess.run([CAPS, "ff", "type", path, "--ff", ffjson, "--typing", rules, "-o", out], capture_output=True, text=True)
    if not os.path.exists(out):
        failed.append((name, (p.stderr or p.stdout).strip()[-150:]))
        continue
    got = [x.strip() for x in open(out)][:len(types)]
    bad = [(f"{t}{i + 1}", t, g) for i, (t, g) in enumerate(zip(types, got)) if t != g and not same(t, g)]
    for t, g in zip(types, got):
        if t != g and same(t, g):
            equivalent[(t, g)] += 1
    if name in deviations:
        deviating.append((name, bad))
        continue
    total += len(types)
    right += len(types) - len(bad)
    tested_names.append(name)
    q = sum(by_name[t].get("charge", 0.0) for t in types)   # the file's own labels with the library's type charges
    if abs(q - round(q)) > 1e-6:
        fractional.append((name, q))
    for _, e, g in bad:
        confusion[(e, g)] += 1
    if bad:
        wrong.append((name, bad))
tested = len(files) - sum(skipped.values()) - len(failed) - len(deviating)
print(f"{right}/{total} atoms ({total and right / total * 100:.1f} %) as the files type them, {tested - len(wrong)} of {tested} molecules fully right"
      + ("" if not skipped else " · skipped: " + ", ".join(f"{k} {v}" for k, v in skipped.items())) + (f" · {len(failed)} failed" if failed else ""))
print(f"  {len(tested_names) - len(fractional)} of {len(tested_names)} molecules at an integer charge with the library's type charges"
      + "".join(f"; {n} {q:+.4f}" for n, q in fractional))
if "--list" in sys.argv:
    print("  tested: " + " ".join(tested_names))
for name, bad in deviating:
    print(f"  deviation {name}: {len(bad)} atoms differ; the file: {deviations[name]}")
for (e, g), k in equivalent.items():
    print(f"  equivalent names: {e} as {g} on {k} atoms (the same class, charge and Lennard-Jones)")
for (e, g), k in confusion.most_common(int(arg("--show", "40"))):
    print(f"  {e:>6} → {g:<6} {k}")
for name, why in failed[:8]:
    print(f"  failed {name}: {why}")
if "--all" in sys.argv:
    for name, bad in wrong:
        print(f"--- {name}: " + ", ".join(f"{a} {e}→{g}" for a, e, g in bad))
