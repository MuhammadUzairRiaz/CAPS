#!/usr/bin/env python3
"""Validation structures from DL_FIELD's own molecule templates (.sf MOLECULE blocks): atoms, their DL_FIELD atom
types and CONNECT records, embedded in 3D with a small spring model (no toolkit needed), written as PDB with the
template's MOLECULE_TYPE key as residue name so DL_FIELD recognises it.

usage: template_pdb.py LIB.sf MOLECULE_NAME OUT.pdb [--seed N]
       template_pdb.py LIB.sf --list        (small, self-contained templates)
"""
import math, random, re, sys

COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
       "Si": 1.11, "B": 0.84, "Li": 1.28, "Na": 1.66, "K": 2.03, "Mg": 1.41, "Ca": 1.76, "Al": 1.21, "Zn": 1.22}


ATOM_KEY = {}   # DL_FIELD atom type → library key, from the ATOM_TYPE table of the last file parsed


def parse_sf(path):
    ATOM_KEY.clear()
    types, keys, mols = {}, {}, {}
    sec, cur = None, None
    for line in open(path, errors="replace"):
        s = re.split(r"(?:^|(?<=\s))#", line.split("!")[0])[0].rstrip("\n")   # '#' inside DL_F names is not a comment
        w = s.split()
        if not w:
            continue
        head = w[0].upper()
        if head == "ATOM_TYPE":
            sec = "types"; continue
        if head == "MOLECULE_TYPE":
            sec = "keys"; continue
        if head == "MOLECULE" and len(w) >= 3:
            cur = {"name": w[1], "atoms": [], "connect": {}, "remark": " ".join(w[4:])}
            mols[w[1]] = cur
            sec = "mol"; continue
        if head == "END":
            sec, cur = (None, None) if sec in ("mol", "types", "keys") else (sec, cur)
            continue
        if sec == "types" and len(w) >= 4:
            types.setdefault(w[0], w[2])
            ATOM_KEY.setdefault(w[0], w[1])
        elif sec == "keys" and len(w) >= 2:
            keys[w[0]] = w[1]
        elif sec == "mol":
            if head == "CONNECT" and len(w) >= 4 and w[2] == ">":
                cur["connect"][w[1]] = w[3:]
            elif len(w) >= 2 and w[1] in types:   # atom lines name a DL_FIELD atom type; other lines are directives
                cur["atoms"].append((w[0], w[1]))
    return types, keys, mols


def embed(atoms, bonds, elem, seed=1):
    n = len(atoms)
    rnd = random.Random(seed)
    nb = [[] for _ in range(n)]
    for i, j in bonds:
        nb[i].append(j); nb[j].append(i)
    r0 = {(min(i, j), max(i, j)): COV.get(elem[i], 0.8) + COV.get(elem[j], 0.8) for i, j in bonds}
    # 1-3 targets from the vertex degree
    d13 = {}
    for j in range(n):
        ang = math.radians(120.0 if len(nb[j]) == 3 else 180.0 if len(nb[j]) == 2 and elem[j] == "C" and False else 109.5)
        for a in range(len(nb[j])):
            for b in range(a + 1, len(nb[j])):
                i, k = nb[j][a], nb[j][b]
                ri, rk = r0[(min(i, j), max(i, j))], r0[(min(k, j), max(k, j))]
                d13[(min(i, k), max(i, k))] = math.sqrt(ri * ri + rk * rk - 2 * ri * rk * math.cos(ang))
    # grow outwards from atom 0 on random directions, then relax
    x = [None] * n
    order, seen = [], set()
    for s in range(n):
        if s in seen: continue
        x[s] = [rnd.uniform(-5, 5) for _ in range(3)]
        stack = [s]; seen.add(s)
        while stack:
            i = stack.pop(0); order.append(i)
            for j in nb[i]:
                if j not in seen:
                    v = [rnd.gauss(0, 1) for _ in range(3)]
                    l = math.sqrt(sum(c * c for c in v))
                    x[j] = [x[i][k] + v[k] / l * r0[(min(i, j), max(i, j))] for k in range(3)]
                    seen.add(j); stack.append(j)
    pairs13 = set(d13)
    bonded = set(r0)
    for it in range(3000):
        f = [[0.0, 0.0, 0.0] for _ in range(n)]
        def spring(i, j, d0, k):
            d = [x[j][c] - x[i][c] for c in range(3)]
            r = math.sqrt(sum(c * c for c in d)) or 1e-6
            g = k * (r - d0) / r
            for c in range(3):
                f[i][c] += g * d[c]; f[j][c] -= g * d[c]
        for (i, j), d0 in r0.items(): spring(i, j, d0, 1.0)
        for (i, j), d0 in d13.items(): spring(i, j, d0, 0.3)
        for i in range(n):
            for j in range(i + 1, n):
                if (i, j) in bonded or (i, j) in pairs13: continue
                d = [x[j][c] - x[i][c] for c in range(3)]
                r = math.sqrt(sum(c * c for c in d)) or 1e-6
                if r < 3.0:
                    g = -0.2 * (3.0 - r) / r
                    for c in range(3):
                        f[i][c] += g * d[c]; f[j][c] -= g * d[c]
        step = 0.2
        for i in range(n):
            for c in range(3): x[i][c] += step * f[i][c]
    return x


def write_pdb(path, name, key, atoms, x, elem, bonds):
    with open(path, "w") as fo:
        fo.write(f"REMARK   CAPS validation structure from DL_FIELD template {name}\n")
        for k, ((lab, ty), p) in enumerate(zip(atoms, x)):
            nm = lab if len(lab) == 4 else f" {lab:<3s}"
            # DL_FIELD's layout: 4-character residue in columns 18-21, no chain, sequence number in 22-26
            fo.write(f"HETATM{k + 1:5d} {nm:4s} {key:<4s}{1:5d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00  0.00          {elem[k]:>2s}\n")
        for i in range(len(atoms)):
            js = sorted({j for a, b in bonds for j in ((b,) if a == i else (a,) if b == i else ())})
            for c in range(0, len(js), 4):   # at most four partners per CONECT record
                fo.write("CONECT" + f"{i + 1:5d}" + "".join(f"{j + 1:5d}" for j in js[c:c + 4]) + "\n")
        fo.write("END\n")


def build(sf, molname, out, seed=1):
    types, keys, mols = parse_sf(sf)
    m = mols[molname]
    labels = [a for a, _ in m["atoms"]]
    idx = {a: i for i, a in enumerate(labels)}
    bonds = set()
    for a, js in m["connect"].items():
        for b in js:
            if a in idx and b in idx and a != b:
                bonds.add((min(idx[a], idx[b]), max(idx[a], idx[b])))
    sf_elem = [types.get(t, re.sub(r"[^A-Za-z]", "", a)[:1]) for a, t in m["atoms"]]
    # united atoms (CH2, CH3) and pseudo atoms: geometry from the leading symbol, but the PDB keeps the library's
    # element name (DL_FIELD matches united-atom templates by it)
    elem = [e if e in COV else (e[:2] if e[:2] in COV else e[:1]) for e in sf_elem]
    x = embed(labels, sorted(bonds), elem, seed)
    write_pdb(out, molname, keys.get(molname, molname[:4].upper()), m["atoms"], x, sf_elem, sorted(bonds))


def candidates(sf, max_atoms=40):
    types, keys, mols = parse_sf(sf)
    out = []
    for name, m in mols.items():
        text = (name + " " + m["remark"]).lower()
        if "water" in text or "link" in text or "(link" in text or name not in keys: continue
        if not (3 <= len(m["atoms"]) <= max_atoms) or not m["connect"]: continue
        if any(b not in [a for a, _ in m["atoms"]] for js in m["connect"].values() for b in js): continue   # links to neighbours
        out.append(name)
    return out


if __name__ == "__main__":
    if "--list" in sys.argv:
        print("\n".join(candidates(sys.argv[1])))
    else:
        seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 1
        build(sys.argv[1], sys.argv[2], sys.argv[3], seed)
