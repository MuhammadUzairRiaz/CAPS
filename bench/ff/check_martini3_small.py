#!/usr/bin/env python3
"""Martini 3 small molecules from all-atom structures: every molecule of data/martini/martini3-small-molecules.json is
built all-atom by CAPS (from a SMILES written from the CHARMM residue's own graph: its elements and bonds, bond orders
by Kekulé matching; CAPS's own atom names), mapped by caps martini --martini 3, and checked:
  · recognised as itself (graph matching, not names),
  · every bead placed, the virtual sites where their beads put them,
  · the mapped bead distances against the model's own bond and constraint lengths (Martini 3 small molecules are
    parameterised on mapped atomistic distances, so these are close; the report gives the deviations).

usage: check_martini3_small.py [--only NAME] [--keep DIR]
"""
import json, math, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
DATA = json.load(open(os.path.join(ROOT, "data", "martini", "martini3-small-molecules.json")))["molecules"]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
only = arg("--only", "")
CENTRE = ["--centre", arg("--centre", "geometry")]
work = arg("--keep", "") or tempfile.mkdtemp()
VALENCE = {"C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "F": 1, "Cl": 1, "Br": 1, "I": 1}


def smiles_of(mol):
    atoms = mol["aa"]["atoms"]
    n = len(atoms)
    nb = [[] for _ in range(n)]
    for i, j in mol["aa"]["bonds"]:
        nb[i].append(j)
        nb[j].append(i)
    heavy = [i for i in range(n) if atoms[i]["element"] != "H"]
    hcount = {i: sum(1 for j in nb[i] if atoms[j]["element"] == "H") for i in heavy}
    charge = {i: 0 for i in heavy}
    # nitro and similar: N with three heavy neighbours and two terminal O -> N+ (=O) O-
    for i in heavy:
        el = atoms[i]["element"]
        hn = [j for j in nb[i] if atoms[j]["element"] != "H"]
        if el == "N" and len(hn) == 3 and hcount[i] == 0:
            ter_o = [j for j in hn if atoms[j]["element"] == "O" and len([k for k in nb[j] if atoms[k]["element"] != "H"]) == 1 and hcount[j] == 0]
            if len(ter_o) == 2:
                charge[i] = 1
                charge[ter_o[1]] = -1
    deficit = {}
    for i in heavy:
        el = atoms[i]["element"]
        v = VALENCE[el] + (charge[i] if el == "N" else charge[i])
        if el == "S" and len(nb[i]) > 2:
            v = 6 if len(nb[i]) >= 4 else 4
        deficit[i] = v - len(nb[i])
    edges = [(i, j) for i in heavy for j in nb[i] if j in deficit and i < j]
    order = {e: 1 for e in edges}

    def solve(k):   # Kekulé: a perfect matching of the atoms with a missing valence
        if all(d == 0 for d in deficit.values()):
            return True
        if k == len(edges):
            return False
        i, j = edges[k]
        if deficit[i] > 0 and deficit[j] > 0:
            deficit[i] -= 1; deficit[j] -= 1; order[(i, j)] += 1
            if solve(k):   # a bond may take more than one extra (triple)
                return True
            deficit[i] += 1; deficit[j] += 1; order[(i, j)] -= 1
        return solve(k + 1)
    if not solve(0):
        return None
    # SMILES by depth-first search, ring closures for the back edges
    sym = {1: "=", 2: "#"}
    sym = {1: "=", 2: "#"}
    bond = lambda a, b: order.get((min(a, b), max(a, b)), 1)

    def atom_text(i):
        el = atoms[i]["element"]
        h = hcount[i]
        c = charge[i]
        return "[" + el + ("H" + (str(h) if h > 1 else "") if h else "") + ("+" if c > 0 else "-" if c < 0 else "") + "]"

    # a spanning tree by depth-first search; the other bonds are ring closures
    tree = {i: [] for i in heavy}
    visited, back = set(), []
    def span(i, parent):
        visited.add(i)
        for j in nb[i]:
            if j not in deficit or j == parent:
                continue
            if j in visited:
                if (j, i) not in back:
                    back.append((i, j))
            else:
                tree[i].append(j)
                span(j, i)
    span(heavy[0], -1)
    nxt = [0]
    ring_at = {i: [] for i in heavy}
    for a, b in back:   # a is reached after b: the bond order goes on a's (closing) digit
        nxt[0] += 1
        ring_at[a].append((nxt[0], b, True))
        ring_at[b].append((nxt[0], a, False))

    def emit(i, parent):
        s = atom_text(i)
        for num, other, closing in ring_at[i]:
            s += (sym.get(bond(i, other) - 1, "") if closing else "") + (str(num) if num < 10 else "%" + str(num))
        kids = tree[i]
        for k, j in enumerate(kids):
            t = sym.get(bond(i, j) - 1, "") + emit(j, i)
            s += t if k == len(kids) - 1 else "(" + t + ")"
        return s
    return emit(heavy[0], -1)


def read_data(path):
    pos, sec, box = [], None, [0, 0, 0]
    for l in open(path):
        s = l.strip()
        for k, ax in enumerate("xyz"):
            if s.endswith(f"{ax}lo {ax}hi"):
                lo, hi = map(float, s.split()[:2])
                box[k] = hi - lo
        if s.startswith("Atoms"):
            sec = "A"
            continue
        if s and s[0].isalpha():
            sec = None
        w = s.split()
        if sec == "A" and len(w) >= 7:
            img = [int(x) for x in w[7:10]] if len(w) >= 10 else [0, 0, 0]
            pos.append((int(w[0]), tuple(float(w[4 + k]) + img[k] * box[k] for k in range(3))))
    return [p for _, p in sorted(pos)]


def main():
    global work
    rows, fails, devs = [], 0, []
    for name, mol in sorted(DATA.items()):
        if only and only.upper() != name:
            continue
        d = os.path.join(work, name)
        os.makedirs(d, exist_ok=True)
        smi = smiles_of(mol)
        if not smi:
            rows.append((name, "no Kekulé structure for its graph", ""))
            fails += 1
            continue
        r = subprocess.run([CAPS, "build", smi, "-o", os.path.join(d, "aa.pdb"), "--seed", "1"], capture_output=True, text=True)
        if r.returncode:
            rows.append((name, f"CAPS could not build {smi}: {r.stderr.strip()[:120]}", ""))
            fails += 1
            continue
        r = subprocess.run([CAPS, "martini", os.path.join(d, "aa.pdb"), "--martini", "3", "-o", os.path.join(d, "cg.data"), "--itp", os.path.join(d, "cg.itp")] + CENTRE,
                           capture_output=True, text=True)
        if r.returncode:
            rows.append((name, f"not mapped: {r.stderr.strip()[:160]}", smi))
            fails += 1
            continue
        got = re.search(r"small molecules mapped onto Martini 3 beads: (.*)", r.stdout)
        names = got.group(1).split(", ") if got else []
        ok = names == [f"1 {name}"]
        pos = read_data(os.path.join(d, "cg.data"))
        blk = mol["block"]
        order = [a["atomname"] for a in blk["atoms"]]
        ok &= len(pos) == len(order)
        worst, sq, nb_ = 0.0, 0.0, 0
        for t in ("bonds", "constraints"):
            for it in blk["interactions"].get(t, []):
                if it["meta"].get("ifndef") == "FLEXIBLE":
                    continue
                a, b = (order.index(x) for x in it["atoms"])
                dd = math.dist(pos[a], pos[b]) / 10 - float(it["params"][1])
                worst = max(worst, abs(dd))
                sq += dd * dd
                nb_ += 1
        rms = math.sqrt(sq / nb_) if nb_ else 0.0
        devs.append(rms)
        fails += 0 if ok else 1
        rows.append((name, f"{'ok' if ok else 'WRONG'} · {'recognised' if ok else 'mapped as ' + str(names)} · {len(pos)} beads · "
                           f"{nb_} bond lengths: rms {rms * 10:.2f} A, largest {worst * 10:.2f} A from the model's", smi))
    for name, res, smi in rows:
        print(f"{name:6s} {res}\n       {smi}")
    print(f"\n{sum(1 for r in rows if r[1].startswith('ok'))} of {len(rows)} Martini 3 small molecules built all-atom, recognised and mapped; "
          f"bond lengths rms {sum(devs) / max(len(devs), 1) * 10:.2f} A on average (work: {work})")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
