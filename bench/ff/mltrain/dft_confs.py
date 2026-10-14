#!/usr/bin/env python3
"""The training set, step 2: DFT energies and forces of every conformer, with GPU4PySCF on one GPU, at the level the
pilot checked against ORCA: wB97M-D3(BJ)/def2-TZVPPD, density-fitted (def2-universal-jkfit), grid (99, 590).

    dft_confs.py TRAINSET ENV [ENV ...]

Reads TRAINSET/confs/<env>.xyz and appends one JSON line per conformer to TRAINSET/dft/<env>.jsonl (frame, energy in
Eh, forces in Eh/bohr, convergence, timings). Frames already in the file are skipped, so a task that ran out of time
is resubmitted as it is. An SCF that does not converge is written with converged false and left out of training.
"""
import io, json, os, sys, time

import pyscf
from ase.io import read
from gpu4pyscf.dft import rks

LEVEL = "wB97M-D3(BJ)/def2-TZVPPD, DF def2-universal-jkfit, grid (99,590), GPU4PySCF"


def run_env(trainset, env):
    src = os.path.join(trainset, "confs", env + ".xyz")
    if not os.path.exists(src) or not open(src).read().rstrip().endswith("# done"):
        print(f"{env}: no conformers yet", flush=True)
        return
    out = os.path.join(trainset, "dft", env + ".jsonl")
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try:
                done.add(json.loads(line)["frame"])
            except (ValueError, KeyError):   # a line cut short when the job was stopped
                pass
    text = open(src).read().replace("# done\n", "")   # the completion mark is not extended xyz
    frames = read(io.StringIO(text), index=":", format="extxyz") if text.strip() else []
    for k, at in enumerate(frames):
        if k in done:
            continue
        mol = pyscf.M(atom=[(s, tuple(p)) for s, p in zip(at.get_chemical_symbols(), at.positions)],
                      basis="def2-tzvppd", charge=0, spin=0, verbose=0, max_memory=90000)
        mf = rks.RKS(mol, xc="wb97m-d3bj").density_fit(auxbasis="def2-universal-jkfit")
        mf.grids.atom_grid = (99, 590)
        mf.conv_tol = 1e-10
        mf.max_cycle = 100
        t0 = time.time()
        e = float(mf.kernel())
        t1 = time.time()
        g = mf.Gradients().kernel()
        g = g.get() if hasattr(g, "get") else g
        t2 = time.time()
        with open(out, "a") as f:
            f.write(json.dumps({"frame": k, "kind": at.info.get("kind"), "angle": at.info.get("angle"),
                                "energy_Eh": e, "converged": bool(mf.converged),
                                "forces_Eh_per_bohr": [[-float(x) for x in row] for row in g],
                                "scf_s": round(t1 - t0, 2), "grad_s": round(t2 - t1, 2), "level": LEVEL}) + "\n")
    print(f"{env}: {len(frames)} conformers, {len(frames) - len(done)} computed now", flush=True)


if __name__ == "__main__":
    trainset = sys.argv[1]
    os.makedirs(os.path.join(trainset, "dft"), exist_ok=True)
    for env in sys.argv[2:]:
        try:
            run_env(trainset, env)
        except Exception as ex:
            print(f"{env}: FAILED {ex}", flush=True)
