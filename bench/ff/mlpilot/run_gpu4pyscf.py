#!/usr/bin/env python3
"""Energy and forces of xyz files on one GPU with GPU4PySCF: wB97M-D3(BJ)/def2-TZVPPD, density-fitted
(def2-universal-jkfit), a (99, 590) grid, SCF converged to 1e-10 Eh.

    run_gpu4pyscf.py OUTDIR FILE.xyz [FILE.xyz ...]

The xyz comment line gives charge=... mult=... . Each file gives OUTDIR/<name>.json: energy (Eh), forces (Eh/bohr,
the negative gradient), the dispersion energy, SCF convergence and the wall time of the SCF and of the gradient.
A file whose json exists is skipped, so a job that ran out of time can be resubmitted as it is.
"""
import json, os, sys, time

import pyscf
from gpu4pyscf.dft import rks

LEVEL = "wB97M-D3(BJ)/def2-TZVPPD, DF def2-universal-jkfit, grid (99,590), GPU4PySCF"


def run(path, outdir):
    out = os.path.join(outdir, os.path.basename(path)[:-4] + ".json")
    if os.path.exists(out):
        return
    lines = open(path).read().splitlines()
    n = int(lines[0])
    meta = dict(w.split("=", 1) for w in lines[1].split() if "=" in w)
    mol = pyscf.M(atom="\n".join(lines[2:2 + n]), basis="def2-tzvppd", charge=int(meta.get("charge", 0)),
                  spin=int(meta.get("mult", 1)) - 1, verbose=0, max_memory=50000)
    mf = rks.RKS(mol, xc="wb97m-d3bj").density_fit(auxbasis="def2-universal-jkfit")
    mf.grids.atom_grid = (99, 590)
    mf.conv_tol = 1e-10
    t0 = time.time()
    e = float(mf.kernel())
    t1 = time.time()
    g = mf.Gradients().kernel()
    t2 = time.time()
    g = g.get() if hasattr(g, "get") else g   # a CuPy array
    json.dump({"file": os.path.basename(path), "level": LEVEL, "energy_Eh": e, "converged": bool(mf.converged),
               "dispersion_Eh": float(mf.scf_summary.get("dispersion", 0.0)),
               "forces_Eh_per_bohr": [[-float(x) for x in row] for row in g],
               "atoms": n, "nao": int(mol.nao), "scf_s": round(t1 - t0, 2), "grad_s": round(t2 - t1, 2)},
              open(out, "w"))
    print(f"{os.path.basename(path)}: E {e:.8f}  scf {t1 - t0:.1f} s  grad {t2 - t1:.1f} s  nao {mol.nao}", flush=True)


if __name__ == "__main__":
    os.makedirs(sys.argv[1], exist_ok=True)
    for p in sys.argv[2:]:
        try:
            run(p, sys.argv[1])
        except Exception as ex:   # one bad fragment does not stop the others
            print(f"{os.path.basename(p)}: FAILED {ex}", flush=True)
