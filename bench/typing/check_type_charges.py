#!/usr/bin/env python3
"""Net charge of molecules typed with a library's own type charges (OPLS-AA and the other libraries whose atom types
carry charges): a neutral molecule must come out at 0 and an ion at its formal charge. A group typed with numbers whose
charges do not add up (a carboxylate at -0.9, a phosphate triester with the anion's charges) shows here.

usage: check_type_charges.py --ff oplsaa2024[,oplsaa-gromacs] [--set OUT/set.json OUT/mol] [SMILES[:NAME] ...]
  --set: the molecules of a typing_audit.py set (its element mol2 files); otherwise the SMILES given (built by caps build)
Prints every molecule typed whole whose net charge is off by more than 0.005 e, and a count; exit status 1 if any.
"""
import json, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FFS = arg("--ff", "oplsaa2024").split(",")


def formal(smi):
    """the SMILES's formal charge: [..+], [..-], [..+2] ..."""
    q = 0
    for m in re.finditer(r"\[[^\]]*?([+-])(\d*)\]", smi):
        q += (1 if m.group(1) == "+" else -1) * int(m.group(2) or 1)
    return q


def molecules():
    if "--set" in sys.argv:
        i = sys.argv.index("--set")
        s, d = json.load(open(sys.argv[i + 1])), sys.argv[i + 2]
        return [(n, os.path.join(d, n + ".mol2"), formal(m["smiles"])) for n, m in s.items()]
    tmp = tempfile.mkdtemp()
    out = []
    args = [a for k, a in enumerate(sys.argv[1:], 1) if a != "--ff" and sys.argv[k - 1] != "--ff"]
    for a in args:
        smi, _, name = a.partition(":")
        name = name or re.sub(r"\W", "_", smi)
        p = os.path.join(tmp, name + ".mol2")
        r = subprocess.run([CAPS, "build", smi, "-o", p, "--seed", "1"], capture_output=True, text=True)
        if r.returncode == 0:
            out.append((name, p, formal(smi)))
        else:
            print("build failed:", smi)
    return out


bad = 0
mols = molecules()
for ff in FFS:
    q = {t["name"]: t.get("charge", 0.0) or 0.0 for t in json.load(open(os.path.join(ROOT, "data", "forcefields", ff + ".json")))["atom_types"]}
    whole = 0
    for name, path, z in mols:
        tf = path + "." + ff + ".types"
        subprocess.run([CAPS, "ff", "type", path, "--ff", ff, "-o", tf], capture_output=True, text=True)
        if not os.path.exists(tf):
            continue
        t = [w.strip() for w in open(tf) if w.strip()]
        if any(x not in q for x in t):
            continue   # untyped atoms: CAPS reports the molecule, no charge to check
        whole += 1
        net = sum(q[x] for x in t)
        if abs(net - z) > 0.005:
            bad += 1
            print(f"{ff:16s} {name:34s} net {net:+.4f} e (formal {z:+d})")
    print(f"{ff}: {whole} molecules typed whole")
print(bad, "molecules off their formal charge")
sys.exit(1 if bad else 0)
