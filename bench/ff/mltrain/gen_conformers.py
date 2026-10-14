#!/usr/bin/env python3
"""The training set, step 1: conformers of each census fragment, with AIMNet2 (MIT) on the GPU.

    gen_conformers.py TRAINSET ENV [ENV ...] [--step 15] [--rattle 0.05] [--model aimnet2]

For the fragment's rotatable bond (b-c, from the census) the dihedral a-b-c-d (a, d: the first heavy neighbours, else
H) is turned in STEP-degree steps round the full circle from its starting value. Each point starts from the previous
relaxed one turned by STEP and is optimised by AIMNet2 with the dihedral held (ASE FixInternals, BFGS to 0.05 eV/A with
steps of at most 0.1 A, at most 300 steps): a relaxed scan. A point where the constraint cannot be kept is skipped
and the scan goes on from the last good point. Each relaxed point also gets one distorted copy (every atom displaced by a normal
deviate of RATTLE A), so the training set has forces away from the minima. A point whose bonds changed in the
optimisation (a bond stretched past 1.4 x the covalent radii, or two atoms more than two bonds apart closer than their
covalent radii), or a copy with two atoms under 0.75 A apart, is left out.
Out: TRAINSET/confs/<env>.xyz (extended xyz; info: env, kind relaxed|rattled, angle, converged) with a last line
"# done" when complete, so a resubmitted task skips it.
"""
import os, sys, zlib

import numpy as np
from ase import Atoms
from ase.constraints import FixInternals
from ase.data import covalent_radii
from ase.io import write
from ase.optimize import BFGS

arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
STEP = float(arg("--step", "15"))
RATTLE = float(arg("--rattle", "0.05"))


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
    return el, np.array(pos), bonds


def changed(atoms, bonds, near):
    """did the optimisation break a bond (stretched past 1.4 x the covalent radii) or make one (two atoms more than
    two bonds apart closer than their covalent radii)? A distance criterion alone would count the close contacts of a
    crowded fragment as bonds from the start"""
    d = atoms.get_all_distances()
    r = covalent_radii[atoms.numbers]
    if any(d[i, j] > 1.4 * (r[i] + r[j]) for i, j in bonds):
        return True
    n = len(r)
    return any(d[i, j] < r[i] + r[j] for i in range(n) for j in range(i + 1, n) if j not in near[i])


def side(nb, start, block):
    seen, todo = {start}, [start]
    while todo:
        x = todo.pop()
        for y in nb[x]:
            if y != block and y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


def scan(trainset, env, central, calc):
    out = os.path.join(trainset, "confs", env + ".xyz")
    if os.path.exists(out) and open(out).read().rstrip().endswith("# done"):
        return 0
    el, pos, bonds = read_mol2(os.path.join(trainset, "fragments", env + ".mol2"))
    nb = [[] for _ in el]
    for i, j in bonds:
        nb[i].append(j)
        nb[j].append(i)
    b, c = central[0] - 1, central[1] - 1
    pick = lambda x, other: sorted((y for y in nb[x] if y != other), key=lambda y: (el[y] == "H", y))[0]
    a, d = pick(b, c), pick(c, b)
    moving = sorted(side(nb, c, b))
    atoms = Atoms(el, positions=pos)
    near = [set(nb[i]) | {k for j in nb[i] for k in nb[j]} for i in range(len(el))]   # 1-2 and 1-3 neighbours
    rng = np.random.default_rng(zlib.crc32(env.encode()))
    phi0 = atoms.get_dihedral(a, b, c, d)
    frames, current, skipped = [], atoms.copy(), 0
    for k in range(int(round(360 / STEP))):
        phi = (phi0 + k * STEP) % 360
        x = current.copy()
        x.set_dihedral(a, b, c, d, phi, indices=moving)
        x.set_constraint(FixInternals(dihedrals_deg=[[phi, [a, b, c, d]]]))
        x.calc = calc
        opt = BFGS(x, logfile=None, maxstep=0.1)   # short steps: the dihedral constraint follows them
        try:
            conv = bool(opt.run(fmax=0.05, steps=300))
        except Exception:   # the constraint could not be satisfied after a step: this point only is left out
            skipped += 1
            continue
        x.set_constraint()
        if changed(x, bonds, near):   # the optimisation broke or made a bond: not this fragment any more
            continue
        current = x
        r = Atoms(el, positions=x.positions)
        r.info.update({"env": env, "kind": "relaxed", "angle": round(phi, 2), "converged": conv})
        frames.append(r)
        y = Atoms(el, positions=x.positions + rng.normal(0, RATTLE, x.positions.shape))
        if y.get_all_distances()[np.triu_indices(len(el), 1)].min() >= 0.75:
            y.info.update({"env": env, "kind": "rattled", "angle": round(phi, 2), "converged": conv})
            frames.append(y)
    write(out, frames, format="extxyz")
    with open(out, "a") as f:
        f.write("# done\n")
    return len(frames)


def main():
    from aimnet.calculators import AIMNet2ASE
    import json
    trainset = sys.argv[1]
    envs = [x for x in sys.argv[2:] if not x.startswith("--") and sys.argv[sys.argv.index(x) - 1] not in
            ("--step", "--rattle", "--model")]
    manifest = {e["env"]: e for e in json.load(open(os.path.join(trainset, "manifest.json")))}
    os.makedirs(os.path.join(trainset, "confs"), exist_ok=True)
    calc = AIMNet2ASE(arg("--model", "aimnet2"))
    for env in envs:
        try:
            n = scan(trainset, env, manifest[env]["central_atoms"], calc)
            print(f"{env}: {n} conformers", flush=True)
        except Exception as ex:   # one fragment's failure does not stop the task
            print(f"{env}: FAILED {ex}", flush=True)


if __name__ == "__main__":
    main()
