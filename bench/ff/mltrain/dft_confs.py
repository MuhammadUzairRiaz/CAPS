#!/usr/bin/env python3
"""The training set, step 2: DFT energies and forces of every conformer, with GPU4PySCF on one GPU, at the level the
pilot checked against ORCA: wB97M-D3(BJ)/def2-TZVPPD, density-fitted (def2-universal-jkfit), grid (99, 590).

    dft_confs.py TRAINSET ENV [ENV ...]

Reads TRAINSET/confs/<env>.xyz and appends one JSON line per conformer to TRAINSET/dft/<env>.jsonl (frame, energy in
Eh, forces in Eh/bohr, convergence, timings). Frames already in the file are skipped, so a task that ran out of time
is resubmitted as it is. An SCF that does not converge is written with converged false and left out of training.
Fragments over DF_MAX_ATOMS atoms, and any conformer that runs out of GPU memory with density fitting, are computed
with the exact integrals (the pilot: 0.08 kcal/mol from the fitted energy, the same forces); one failed conformer does
not stop the fragment's others.
"""
import io, json, os, sys, time

import pyscf
from ase.io import read
from gpu4pyscf.dft import rks

LEVEL = "wB97M-D3(BJ)/def2-TZVPPD, DF def2-universal-jkfit, grid (99,590), GPU4PySCF"
LEVEL_NODF = "wB97M-D3(BJ)/def2-TZVPPD, no DF, grid (99,590), GPU4PySCF"
DF_MAX_ATOMS = 60   # 70-atom fragments ran the 94 GB H100 out of memory with density fitting


def point(at, df):
    """energy (Eh), gradient (Eh/bohr), convergence and timings of one conformer"""
    mol = pyscf.M(atom=[(s, tuple(p)) for s, p in zip(at.get_chemical_symbols(), at.positions)],
                  basis="def2-tzvppd", charge=0, spin=0, verbose=0, max_memory=90000)
    mf = rks.RKS(mol, xc="wb97m-d3bj")
    if df:
        mf = mf.density_fit(auxbasis="def2-universal-jkfit")
    mf.grids.atom_grid = (99, 590)
    mf.conv_tol = 1e-10
    mf.max_cycle = 100
    t0 = time.time()
    e = float(mf.kernel())
    t1 = time.time()
    g = mf.Gradients().kernel()
    g = g.get() if hasattr(g, "get") else g
    return e, g, bool(mf.converged), round(t1 - t0, 2), round(time.time() - t1, 2)


def free_gpu():
    """hand the GPU memory pool back, so a failed large conformer does not leave the next one short"""
    try:
        import cupy
        cupy.get_default_memory_pool().free_all_blocks()
    except Exception:
        pass


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
        # a large fragment's density-fitting tensor does not fit on the GPU: those are computed with the exact
        # integrals from the start (the whole fragment at one setting), and any conformer that still runs out of
        # memory is retried so; each line says which setting it used
        df = len(at) <= DF_MAX_ATOMS
        try:
            e, g, conv, t_scf, t_grad = point(at, df)
        except Exception as ex:
            if not df or "memory" not in str(ex).lower():
                print(f"{env} frame {k}: FAILED {ex}", flush=True)
                continue
            free_gpu()
            df = False
            try:
                e, g, conv, t_scf, t_grad = point(at, df)
            except Exception as ex2:
                print(f"{env} frame {k}: FAILED without DF too: {ex2}", flush=True)
                continue
        free_gpu()
        with open(out, "a") as f:
            f.write(json.dumps({"frame": k, "kind": at.info.get("kind"), "angle": at.info.get("angle"),
                                "energy_Eh": e, "converged": conv,
                                "forces_Eh_per_bohr": [[-float(x) for x in row] for row in g],
                                "scf_s": t_scf, "grad_s": t_grad, "level": LEVEL if df else LEVEL_NODF}) + "\n")
    print(f"{env}: {len(frames)} conformers, {len(frames) - len(done)} computed now", flush=True)


if __name__ == "__main__":
    trainset = sys.argv[1]
    os.makedirs(os.path.join(trainset, "dft"), exist_ok=True)
    for env in sys.argv[2:]:
        try:
            run_env(trainset, env)
        except Exception as ex:
            print(f"{env}: FAILED {ex}", flush=True)
