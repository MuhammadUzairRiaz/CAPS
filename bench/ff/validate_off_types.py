#!/usr/bin/env python3
"""Check CAPS's automatic typing against the atom types in AMBER's own residue libraries (LEaP object files, .lib / .off:
all_amino03.lib, all_nucleic94.lib ...; prep files: GLYCAM_06j-1.prep), whose residues the force field's developers typed.

Each residue becomes a mol2 file with its bonds but no bond orders. A residue bonded to its neighbours (the library's
head / tail atoms and any further connect atoms, e.g. a disulfide's S) gets a small cap on each link, chosen by the
partner the link reaches: an acetyl for an N (the previous residue's C), an N-methyl amide for a C, a methoxy for a P
(the previous nucleotide's O3'), a methyl phosphate for an O3', a methylthio for a disulfide S. Only the residue's own
atoms are compared. A prep file names no links: an sp3 atom short of its neighbours is one, capped by what GLYCAM joins
there (a methoxy on a sugar C, an oxan-2-yl on a glycosidic O, a methyl on an N).

usage: validate_off_types.py FF.json RULES.json LIB.lib [LIB.lib ...] [--show N] [--all] [--only RES,...] [--deviations FILE]

--deviations: a JSON file {"residues": {name: evidence}} of residues whose library entry is itself wrong; their
mismatches are reported apart.
"""
import collections, json, math, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
ffjson, rules = sys.argv[1], sys.argv[2]
libs = [a for a in sys.argv[3:] if not a.startswith("--") and (sys.argv[sys.argv.index(a) - 1] not in ("--show", "--only", "--deviations"))]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
show = int(arg("--show", "40"))
only = [x for x in arg("--only", "").split(",") if x]
known = {t["name"] for t in json.load(open(ffjson))["atom_types"]}
deviations = json.load(open(arg("--deviations", "")))["residues"] if "--deviations" in sys.argv else {}
deviating = []
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
OXANYL = [("C", [1, 5, 6]), ("C", [2, 7, 8]), ("C", [3, 9, 10]), ("C", [4, 11, 12]), ("C", [5, 13, 14]), ("O", [])] + [("H", [])] * 9
OXANYLOXY = [("O", [1])] + [(e, [j + 1 for j in js]) for e, js in OXANYL]
# the partner GLYCAM joins there: a sugar's glycosidic O for a C (the anomeric C, an O-acetyl's or O-methyl's C) or an S
# (sulfate), a sugar C for an O, a methyl for an N
CAP_PREP = {"C": OXANYLOXY, "S": OXANYLOXY, "O": OXANYL, "N": METHYL}
# neighbours a complete atom of these prep types has (fewer: the rest are links)
PREP_FULL = {**{t: 4 for t in ("Cg", "Cy", "CT", "2C", "3C", "CX", "Cp")}, **{t: 2 for t in ("Os", "Oh", "Oy")},
             **{t: 3 for t in ("Ng", "NT", "C")}, "N3": 4, "S": 4}


def element_of(t):
    t = t.lstrip("0123456789")
    return {"C0": "Ca", "Na": "Na", "Cl": "Cl"}.get(t, t[:1].upper())


def zmat_place(pos, a, b, c, r, theta, phi):
    """an atom r from a, at angle theta to b, dihedral phi to c (degrees; the prep file's internal coordinates)"""
    if a not in pos:
        return (0.0, 0.0, 0.0)
    pa = pos[a]
    if b not in pos:
        return (pa[0] + r, pa[1], pa[2])
    pb = pos[b]
    pc = pos[c] if c in pos else (pb[0], pb[1] + 1.0, pb[2])
    sub = lambda u, v: [u[k] - v[k] for k in range(3)]
    crs = lambda u, v: [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
    unit = lambda u: [x / math.sqrt(sum(y * y for y in u)) for x in u]
    bc = unit(sub(pa, pb))
    n = crs(sub(pb, pc), bc)
    if sum(x * x for x in n) < 1e-12:   # collinear references: any perpendicular
        n = crs(bc, [1.0, 0.0, 0.0] if abs(bc[0]) < 0.9 else [0.0, 1.0, 0.0])
    n = unit(n)
    m = crs(n, bc)
    th, ph = math.radians(theta), math.radians(phi)
    d = [-r * math.cos(th), r * math.sin(th) * math.cos(ph), r * math.sin(th) * math.sin(ph)]
    return tuple(pa[k] + d[0] * bc[k] + d[1] * m[k] + d[2] * n[k] for k in range(3))


def prep_units(path):
    """the residues of an AMBER prep file: atoms with their types (dummy atoms left out), each bonded to its first
    reference atom, the LOOP closures, and as links the sp3 atoms short of neighbours"""
    L = open(path, errors="replace").read().splitlines()
    out, i = {}, 0
    while i < len(L):
        if i + 1 < len(L) and L[i + 1].split()[:1] in (["CORRECT"], ["CHANGE"]) and len(L[i].split()) >= 2 and L[i].split()[1] == "INT":
            name, atoms, bonds, idx, pos, xyz = L[i].split()[0], [], [], {}, {}, []
            i += 3
            while i < len(L) and L[i].strip():
                w = L[i].split()
                pos[int(w[0])] = zmat_place(pos, int(w[4]), int(w[5]), int(w[6]), float(w[7]), float(w[8]), float(w[9]))
                if w[2] != "DU":
                    xyz.append(pos[int(w[0])])
                    idx[int(w[0])] = len(atoms)
                    atoms.append((w[1], w[2], element_of(w[2])))
                    if int(w[4]) in idx:
                        bonds.append((idx[int(w[4])], idx[int(w[0])]))
                i += 1
            names = {a: k for k, (a, _, _) in enumerate(atoms)}
            while i < len(L) and L[i].strip() != "DONE":
                if L[i].strip() == "LOOP":
                    i += 1
                    while i < len(L) and L[i].strip():
                        a, b = L[i].split()[:2]
                        bonds.append((names[a], names[b]))
                        i += 1
                i += 1
            deg = collections.Counter(x for b in bonds for x in b)
            links = {k for k, (_, t, e) in enumerate(atoms) if deg[k] < PREP_FULL.get(t, 0)}
            out[name] = {"atoms": atoms, "bonds": bonds, "links": links, "prep": True, "xyz": xyz}
        i += 1
    return out


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
            cur["atoms"].append((w[0].strip('"'), w[1].strip('"'), SYMBOL.get(int(w[6]))))
        elif sec == "connectivity":
            cur["bonds"].append((int(w[0]) - 1, int(w[1]) - 1))
        elif sec in ("connect", "residueconnect"):
            cur["links"] |= {int(x) - 1 for x in w if int(x) > 0}
    return out


work = os.environ.get("OFF_KEEP") or tempfile.mkdtemp()
total = right = 0
confusion, wrong, skipped, failed, tested = collections.Counter(), [], collections.Counter(), [], 0
for lib in libs:
    for name, u in (prep_units(lib) if lib.endswith(".prep") else units(lib)).items():
        if only and name not in only:
            continue
        names = [a for a, _, _ in u["atoms"]]
        types = [t for _, t, _ in u["atoms"]]
        elems = [e for _, _, e in u["atoms"]]
        if not names or any(e is None for e in elems) or any(t.startswith("EP") or t == "DU" for t in types):
            skipped["lone pairs / dummy atoms"] += 1
            continue
        if any(t not in known for t in types):
            skipped["types not in the force-field file"] += 1
            continue
        bonds, cap_atoms, cap_xyz = {(min(i, j), max(i, j)) for i, j in u["bonds"]}, [], []
        xyz = u.get("xyz") or [(i * 1.5, 0.0, 0.0) for i in range(len(names))]
        degree = collections.Counter(x for b in bonds for x in b)
        for a in sorted(u["links"]):
            if not u.get("prep") and degree[a] >= FULL.get(elems[a], 9):   # a cap residue's nominal head or tail (ACE's HH31): no bond
                continue
            cap = (CAP_PREP if u.get("prep") else CAP).get(elems[a])
            if cap is None:
                skipped["links the check cannot cap"] += 1
                break
            base = len(names) + len(cap_atoms)
            cap_atoms += [e for e, _ in cap]
            # the cap's first atom 1.45 A out along the free valence (opposite the existing bonds), the rest beyond it
            nbv = [[xyz[o][k] - xyz[a][k] for k in range(3)] for b in bonds for o in b if a in b and o != a]
            out_dir = [-sum(v[k] / math.sqrt(sum(x * x for x in v)) for v in nbv) for k in range(3)] if nbv else [1.0, 0.0, 0.0]
            ln = math.sqrt(sum(x * x for x in out_dir)) or 1.0
            first = [xyz[a][k] + 1.45 * out_dir[k] / ln for k in range(3)]
            cap_xyz += [tuple(first[k] + (1.2 * j if k == 0 else 0.3 * j) * (out_dir[k] / ln if k == 0 else 1) for k in range(3)) for j in range(len(cap))]
            bonds.add((a, base))
            for k, (_, js) in enumerate(cap):
                bonds |= {(base + k, base + j) for j in js}
        else:
            path = os.path.join(work, re.sub(r"[^A-Za-z0-9_]", "_", name) + ".mol2")
            with open(path, "w") as f:
                f.write(f"@<TRIPOS>MOLECULE\n{name}\n{len(names) + len(cap_atoms)} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
                for i, a in enumerate(names):
                    f.write(f"{i + 1} {a} {xyz[i][0]:.4f} {xyz[i][1]:.4f} {xyz[i][2]:.4f} {elems[i]} 1 {name[:4]} 0.0\n")
                for k, e in enumerate(cap_atoms):
                    i = len(names) + k
                    f.write(f"{i + 1} X{k} {cap_xyz[k][0]:.4f} {cap_xyz[k][1]:.4f} {cap_xyz[k][2]:.4f} {e} 1 CAP 0.0\n")
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
            if name in deviations:   # the library's own entry is wrong: reported apart
                deviating.append((name, bad))
                tested -= 1
                continue
            total += len(types)
            right += len(types) - len(bad)
            for _, e, g in bad:
                confusion[(e, g)] += 1
            if bad:
                wrong.append((name, bad))

print(f"{', '.join(os.path.basename(l) for l in libs)}: {right}/{total} atoms ({total and right / total * 100:.1f} %) as the library, "
      f"{tested - len(wrong)} of {tested} residues fully right" + ("" if not skipped else " · skipped: " +
      ", ".join(f"{k} {v}" for k, v in skipped.items())) + (f" · {len(failed)} failed" if failed else ""))
for name, bad in deviating:
    print(f"  deviation {name}: {len(bad)} atoms differ; the entry: {deviations[name]}")
for (e, g), k in confusion.most_common(show):
    print(f"  {e:>8} → {g:<8} {k}")
for name, why in failed[:5]:
    print(f"  failed {name}: {why}")
if "--all" in sys.argv:
    for name, bad in wrong:
        print(f"--- {name}: " + ", ".join(f"{a} {e}→{g}" for a, e, g in bad))
