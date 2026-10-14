#!/usr/bin/env python3
"""The pilot of CAPS's own training set: does DFT on the GPU (GPU4PySCF, Alpha) give what ORCA gives, and how long does
a census fragment take?

    make_pilot.py [--census bench/ff/census/oplsaa2024] [--n 20] [--out DIR]

picks N fragments from the census, spread over the parameter status (missing, CAPS fit, published), the elements
(S, N, Si, Cl, F next to C/H/O) and the fragment size, and writes three conformers of each: the fragment as cut, and
its rotatable bond turned by 120 and 240 degrees (60 / 180 / 300 when a turn puts two atoms closer than 1.5 A).
Out: DIR/xyz/<env>_<k>.xyz (comment: charge, multiplicity, env, angle) and DIR/pilot.json.

Both codes then compute the conformers at one level, wB97M-D3(BJ)/def2-TZVPPD: run_gpu4pyscf.py on Alpha (all of
them) and ORCA on Barnard (orca_inputs.py); compare.py lines them up.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "data", "python"))
from caps.torsionfit import neighbours, rotate, side  # noqa: E402

arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Br": 35, "I": 53}


def read_mol2(path):
    el, pos, bonds, sec = [], [], [], None
    for line in open(path):
        if line.startswith("@<TRIPOS>"):
            sec = line.strip()[9:]
            continue
        w = line.split()
        if sec == "ATOM" and len(w) >= 6:
            el.append(w[5])
            pos.append([float(x) for x in w[2:5]])
        elif sec == "BOND" and len(w) >= 3:
            bonds.append((int(w[1]) - 1, int(w[2]) - 1))
    return el, pos, bonds


def pick(envs, n):
    """a spread: every status with each rare element, then each status small to large, then the commonest
    environments; deterministic"""
    chosen, seen = [], set()

    def take(e):
        if e["env"] not in seen and len(chosen) < n:
            chosen.append(e)
            seen.add(e["env"])
    for el in ("S", "Si", "Cl", "F", "N"):
        for st in ("missing", "CAPS fit", "published"):
            for e in envs:
                if el in e["elements"] and e["status"] == st:
                    take(e)
                    break
    for st in ("missing", "CAPS fit", "published"):
        group = sorted((e for e in envs if e["status"] == st), key=lambda e: (e["fragment_atoms"], e["env"]))
        for k in range(4):   # small to large
            if group:
                take(group[min(len(group) - 1, k * len(group) // 3)])
    for e in sorted(envs, key=lambda e: (-e["count"], e["env"])):
        take(e)
    return chosen


def clash(el, pos, nb, tol=1.5):
    near = [set(nb[i]) | {k for j in nb[i] for k in nb[j]} for i in range(len(el))]
    for i in range(len(el)):
        for j in range(i + 1, len(el)):
            if j not in near[i] and sum((pos[i][k] - pos[j][k]) ** 2 for k in range(3)) < tol * tol:
                return True
    return False


def main():
    census = arg("--census", os.path.join(ROOT, "bench", "ff", "census", "oplsaa2024"))
    out = arg("--out", os.path.join(HERE, "pilot"))
    n = int(arg("--n", "20"))
    envs = json.load(open(os.path.join(census, "census.json")))["environments"]
    os.makedirs(os.path.join(out, "xyz"), exist_ok=True)
    rows = []
    for e in pick(envs, len(envs)):   # in order of preference; skipped ones make room for the next
        if len(rows) == n:
            break
        el, pos, bonds = read_mol2(os.path.join(census, "fragments", e["env"] + ".mol2"))
        if sum(Z[x] for x in el) % 2:   # an odd electron count: not a closed-shell fragment
            print(f"skip {e['env']}: odd electron count")
            continue
        nb = neighbours(len(el), bonds)
        b, c = (x - 1 for x in e["central_atoms"])
        moving = side(nb, c, b)
        confs = [(0, pos)]
        for want in ((120, 60), (240, 180, 300)):
            for a in want:
                p = rotate(pos, moving, b, c, a)
                if all(a != x for x, _ in confs) and not clash(el, p, nb):
                    confs.append((a, p))
                    break
        files = []
        for k, (a, p) in enumerate(confs):
            name = f"{e['env']}_{k}.xyz"
            with open(os.path.join(out, "xyz", name), "w") as f:
                f.write(f"{len(el)}\ncharge=0 mult=1 env={e['env']} angle={a}\n")
                for x, q in zip(el, p):
                    f.write(f"{x:2s} {q[0]:12.6f} {q[1]:12.6f} {q[2]:12.6f}\n")
            files.append(name)
        rows.append({"env": e["env"], "status": e["status"], "elements": e["elements"], "atoms": len(el),
                     "heavy": e["fragment_heavy"], "classes": e["classes"], "example": e["example"], "files": files})
        print(f"{e['env']}  {e['status']:9s} {len(el):3d} atoms  {'/'.join(e['elements']):10s} {len(files)} conformers  "
              f"{e['example'][:50]}")
    json.dump(rows, open(os.path.join(out, "pilot.json"), "w"), indent=1)
    print(f"{len(rows)} fragments, {sum(len(r['files']) for r in rows)} conformers -> {out}")


if __name__ == "__main__":
    main()
