#!/usr/bin/env python3
"""Check CAPS's automatic typing against the atom types in AMBER's own residue libraries (LEaP object files, .lib / .off:
all_amino03.lib, all_nucleic94.lib ...), whose residues the force field's developers typed.

Each residue becomes a mol2 file with its bonds but no bond orders. A residue bonded to its neighbours (the library's
head / tail atoms and any further connect atoms, e.g. a disulfide's S) gets a small cap on each link, chosen by the
partner the link reaches: an acetyl for an N (the previous residue's C), an N-methyl amide for a C, a methoxy for a P
(the previous nucleotide's O3'), a methyl phosphate for an O3', a methylthio for a disulfide S. Only the residue's own
atoms are compared.

usage: validate_off_types.py FF.json RULES.json LIB.lib [LIB.lib ...] [--show N] [--all] [--only RES,...]
"""
import collections, json, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
ffjson, rules = sys.argv[1], sys.argv[2]
libs = [a for a in sys.argv[3:] if not a.startswith("--") and (sys.argv[sys.argv.index(a) - 1] not in ("--show", "--only"))]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
show = int(arg("--show", "40"))
only = [x for x in arg("--only", "").split(",") if x]
known = {t["name"] for t in json.load(open(ffjson))["atom_types"]}
SYMBOL = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na", 12: "Mg", 15: "P", 16: "S", 17: "Cl", 19: "K", 20: "Ca",
          35: "Br", 37: "Rb", 53: "I", 55: "Cs", 3: "Li", 30: "Zn"}
METHYL = [("C", [1, 2, 3]), ("H", []), ("H", []), ("H", [])]
CAP = {   # by the residue's own linking atom: the cap standing in for the neighbour it bonds to
    "N": [("C", [1, 2]), ("O", []), ("C", [3, 4, 5]), ("H", []), ("H", []), ("H", [])],                       # acetyl
    "C": [("N", [1, 2]), ("H", []), ("C", [3, 4, 5]), ("H", []), ("H", []), ("H", [])],                       # N-methyl amide
    "P": [("O", [1])] + [(e, [j + 1 for j in js]) for e, js in METHYL],                                       # methoxy
    "O": [("P", [1, 2, 3]), ("O", []), ("O", []), ("O", [4])] + [(e, [j + 4 for j in js]) for e, js in METHYL],  # methyl phosphate
    "S": [("S", [1])] + [(e, [j + 1 for j in js]) for e, js in METHYL],                                       # methylthio
}


FULL = {"H": 1, "C": 4, "N": 3, "O": 2, "S": 2, "P": 4}   # neighbours of an atom with no free valence


def units(path):
    """{name: {"atoms": [(name, type, element)], "bonds": [(i, j)], "links": [atom index]}} from a LEaP object file"""
    out, cur, sec = {}, None, None
    for line in open(path, errors="replace"):
        if line.startswith("!"):
            m = re.match(r"!entry\.([^.]+)\.unit\.(\w+)", line)
            sec = None
            if m:
                cur = out.setdefault(m.group(1), {"atoms": [], "bonds": [], "links": set()})
                sec = m.group(2)
            continue
        w = line.split()
        if cur is None or not w:
            continue
        if sec == "atoms":
            cur["atoms"].append((w[0].strip('"'), w[1].strip('"'), int(w[6])))
        elif sec == "connectivity":
            cur["bonds"].append((int(w[0]) - 1, int(w[1]) - 1))
        elif sec in ("connect", "residueconnect"):
            cur["links"] |= {int(x) - 1 for x in w if int(x) > 0}
    return out


work = os.environ.get("OFF_KEEP") or tempfile.mkdtemp()
total = right = 0
confusion, wrong, skipped, failed, tested = collections.Counter(), [], collections.Counter(), [], 0
for lib in libs:
    for name, u in units(lib).items():
        if only and name not in only:
            continue
        names = [a for a, _, _ in u["atoms"]]
        types = [t for _, t, _ in u["atoms"]]
        elems = [SYMBOL.get(z) for _, _, z in u["atoms"]]
        if not names or any(e is None for e in elems) or any(t.startswith("EP") or t == "DU" for t in types):
            skipped["lone pairs / dummy atoms"] += 1
            continue
        if any(t not in known for t in types):
            skipped["types not in the force-field file"] += 1
            continue
        bonds, cap_atoms = {(min(i, j), max(i, j)) for i, j in u["bonds"]}, []
        degree = collections.Counter(x for b in bonds for x in b)
        for a in sorted(u["links"]):
            if degree[a] >= FULL.get(elems[a], 9):   # a cap residue's nominal head or tail (ACE's HH31, NME's CH3): no bond
                continue
            cap = CAP.get(elems[a])
            if cap is None:
                skipped["links the check cannot cap"] += 1
                break
            base = len(names) + len(cap_atoms)
            cap_atoms += [e for e, _ in cap]
            bonds.add((a, base))
            for k, (_, js) in enumerate(cap):
                bonds |= {(base + k, base + j) for j in js}
        else:
            path = os.path.join(work, re.sub(r"[^A-Za-z0-9_]", "_", name) + ".mol2")
            with open(path, "w") as f:
                f.write(f"@<TRIPOS>MOLECULE\n{name}\n{len(names) + len(cap_atoms)} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
                for i, a in enumerate(names):
                    f.write(f"{i + 1} {a} {i * 1.5:.3f} 0 0 {elems[i]} 1 {name[:4]} 0.0\n")
                for k, e in enumerate(cap_atoms):
                    i = len(names) + k
                    f.write(f"{i + 1} X{k} {i * 1.5:.3f} 0 0 {e} 1 CAP 0.0\n")
                f.write("@<TRIPOS>BOND\n")
                for k, (i, j) in enumerate(sorted(bonds)):
                    f.write(f"{k + 1} {i + 1} {j + 1} un\n")
            out = path + ".types"
            p = subprocess.run([CAPS, "ff", "type", path, "--ff", ffjson, "--typing", rules, "-o", out], capture_output=True, text=True)
            if not os.path.exists(out):
                failed.append((name, (p.stderr or p.stdout).strip()[-120:]))
                continue
            tested += 1
            got = [x.strip() for x in open(out)][:len(names)]
            bad = [(names[i], types[i], got[i]) for i in range(len(types)) if types[i] != got[i]]
            total += len(types)
            right += len(types) - len(bad)
            for _, e, g in bad:
                confusion[(e, g)] += 1
            if bad:
                wrong.append((name, bad))

print(f"{', '.join(os.path.basename(l) for l in libs)}: {right}/{total} atoms ({total and right / total * 100:.1f} %) as the library, "
      f"{tested - len(wrong)} of {tested} residues fully right" + ("" if not skipped else " · skipped: " +
      ", ".join(f"{k} {v}" for k, v in skipped.items())) + (f" · {len(failed)} failed" if failed else ""))
for (e, g), k in confusion.most_common(show):
    print(f"  {e:>8} → {g:<8} {k}")
for name, why in failed[:5]:
    print(f"  failed {name}: {why}")
if "--all" in sys.argv:
    for name, bad in wrong:
        print(f"--- {name}: " + ", ".join(f"{a} {e}→{g}" for a, e, g in bad))
