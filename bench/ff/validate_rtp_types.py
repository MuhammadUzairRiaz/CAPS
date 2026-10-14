#!/usr/bin/env python3
"""Check CAPS's automatic typing against the atom types in a GROMACS residue topology (.rtp: aminoacids.rtp of
oplsaa.ff ...), whose residues the force field's developers typed.

Each residue becomes a mol2 file with its bonds but no bond orders; a bond to the previous or next residue (-C, +N,
-O3', +P ...) is capped as validate_rtf_types.py caps it (an acetyl for -C, an N-methyl amide for +N, a methoxy for a
link O, a methyl phosphate for a link P). Only the residue's own atoms are compared. Elements come from the force
field's types.

usage: validate_rtp_types.py FILE.rtp FF.json RULES.json [--show N] [--all] [--only RES,...] [--exclude RES,...]
                             [--deviations FILE]

RULES may be "base,overlay" (a water model's overlay, as caps ff type takes it). --deviations: a JSON file
{"residues": {name: evidence}} of residues whose topology entry is itself wrong, reported apart.
"""
import collections, json, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
rtp, ffjson, rules = sys.argv[1], sys.argv[2], sys.argv[3]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
show = int(arg("--show", "40"))
only = [x for x in arg("--only", "").split(",") if x]
exclude = [x for x in arg("--exclude", "").split(",") if x]
deviations = json.load(open(arg("--deviations", "")))["residues"] if "--deviations" in sys.argv else {}
deviating = []
ff = json.load(open(ffjson))
element = {t["name"]: t["element"] for t in ff["atom_types"]}
METHYL = [("C", [1, 2, 3]), ("H", []), ("H", []), ("H", [])]
CAP = {"C": [("C", [1, 2]), ("O", []), ("C", [3, 4, 5]), ("H", []), ("H", []), ("H", [])],
       "N": [("N", [1, 2]), ("H", []), ("C", [3, 4, 5]), ("H", []), ("H", []), ("H", [])],
       "O": [("O", [1])] + [(e, [j + 1 for j in js]) for e, js in METHYL],
       "P": [("P", [1, 2, 3]), ("O", []), ("O", []), ("O", [4])] + [(e, [j + 4 for j in js]) for e, js in METHYL],
       "S": [("S", [1])] + [(e, [j + 1 for j in js]) for e, js in METHYL]}


def residues(path):
    out, cur, sec = [], None, None
    for line in open(path, errors="replace"):
        s = line.split(";")[0].strip()
        if not s:
            continue
        m = re.match(r"\[\s*(\S+)\s*\]", s)
        if m:
            name = m.group(1)
            if name in ("atoms", "bonds", "impropers", "dihedrals", "exclusions", "cmap", "angles"):
                sec = name
            elif name == "bondedtypes":
                cur, sec = None, None
            else:
                cur = {"name": name, "atoms": [], "bonds": []}
                out.append(cur)
                sec = None
            continue
        if cur is None:
            continue
        w = s.split()
        if sec == "atoms" and len(w) >= 2:
            cur["atoms"].append((w[0], w[1]))
        elif sec == "bonds" and len(w) >= 2:
            cur["bonds"].append((w[0], w[1]))
    return out


work = os.environ.get("RTP_KEEP") or tempfile.mkdtemp()
total = right = 0
confusion, wrong, skipped, failed, tested = collections.Counter(), [], collections.Counter(), [], 0
for r in residues(rtp):
    if (only and r["name"] not in only) or r["name"] in exclude:
        continue
    names = [a for a, _ in r["atoms"]]
    types = [t for _, t in r["atoms"]]
    if any(t not in element for t in types):
        skipped["types not in the force-field file"] += 1
        continue
    if any(element[t] in ("X", "") for t in types):
        skipped["virtual sites"] += 1
        continue
    idx = {a: i for i, a in enumerate(names)}
    bonds, cap_atoms, ok = set(), [], True
    links = list(r["bonds"])
    have = {(a, b) for a, b in links} | {(b, a) for a, b in links}
    for a, b in r["bonds"]:   # each peptide / phosphodiester bond is listed once: "+X" from atom A means this residue's
        for own, link in ((a, b), (b, a)):   # own X bonds to the previous residue's A; "-X" from A, that X to the next A
            if own[:1] in "+-" or link[:1] not in "+-" or link[1:] not in idx:
                continue
            twin = ("-" if link[:1] == "+" else "+") + own
            if (link[1:], twin) not in have:
                links.append((link[1:], twin))
    if not any(x[:1] in "+-" for b in links for x in b):   # a cap residue (ACE, NME): its C or N short of neighbours links on
        deg = collections.Counter(x for b in links for x in b)
        for a in ("C", "N"):
            if a in idx and element[types[idx[a]]] == a[0] and deg[a] < (3 if a == "C" else 2):
                links.append((a, "+N" if a == "C" else "-C"))
    deg = collections.Counter(x for b in links for x in b)
    for a in names:   # a disulfide S (its partner joined by GROMACS's specbond.dat): a methylthio
        if element[types[idx[a]]] == "S" and deg[a] == 1:
            links.append((a, "+S"))
    for a, b in links:
        if a[:1] in "+-" or b[:1] in "+-":
            link, own = (a, b) if a[:1] in "+-" else (b, a)
            cap = CAP.get(link[1:2].upper())
            if cap is None or own not in idx:
                ok = False
                break
            base = len(names) + len(cap_atoms)
            cap_atoms += [e for e, _ in cap]
            bonds.add((idx[own], base))
            for k, (_, js) in enumerate(cap):
                bonds |= {(base + k, base + j) for j in js}
            continue
        if a not in idx or b not in idx:
            ok = False
            break
        bonds.add((min(idx[a], idx[b]), max(idx[a], idx[b])))
    if not ok:
        skipped["bonds the check cannot cap"] += 1
        continue
    path = os.path.join(work, re.sub(r"[^A-Za-z0-9_]", "_", r["name"]) + ".mol2")
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{r['name']}\n{len(names) + len(cap_atoms)} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, a in enumerate(names):
            f.write(f"{i + 1} {a} {i * 1.5:.3f} 0 0 {element[types[i]]} 1 {r['name'][:4]} 0.0\n")
        for k, e in enumerate(cap_atoms):
            i = len(names) + k
            f.write(f"{i + 1} X{k} {i * 1.5:.3f} 0 0 {e} 1 CAP 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for k, (i, j) in enumerate(sorted(bonds)):
            f.write(f"{k + 1} {i + 1} {j + 1} un\n")
    out = path + ".types"
    p = subprocess.run([CAPS, "ff", "type", path, "--ff", ffjson, "--typing", rules, "-o", out], capture_output=True, text=True)
    if not os.path.exists(out):
        failed.append((r["name"], (p.stderr or p.stdout).strip()[-120:]))
        continue
    tested += 1
    got = [x.strip() for x in open(out)][:len(names)]
    bad = [(names[i], types[i], got[i]) for i in range(len(types)) if types[i] != got[i]]
    if r["name"] in deviations:   # the topology's own entry is wrong: reported apart
        deviating.append((r["name"], bad))
        tested -= 1
        continue
    total += len(types)
    right += len(types) - len(bad)
    for _, e, g in bad:
        confusion[(e, g)] += 1
    if bad:
        wrong.append((r["name"], bad))

print(f"{os.path.basename(rtp)}: {right}/{total} atoms ({total and right / total * 100:.1f} %) as the topology, "
      f"{tested - len(wrong)} of {tested} residues fully right" + ("" if not skipped else " · skipped: " +
      ", ".join(f"{k} {v}" for k, v in skipped.items())) + (f" · {len(failed)} failed" if failed else ""))
for name, bad in deviating:
    print(f"  deviation {name}: {len(bad)} atoms differ; the entry: {deviations[name]}")
for (e, g), k in confusion.most_common(show):
    print(f"  {e:>10} → {g:<10} {k}")
for name, why in failed[:5]:
    print(f"  failed {name}: {why}")
if "--all" in sys.argv:
    for name, bad in wrong:
        print(f"--- {name}: " + ", ".join(f"{a} {e}→{g}" for a, e, g in bad))
