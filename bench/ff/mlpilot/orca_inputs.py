#!/usr/bin/env python3
"""The ORCA half of the pilot: the same level as run_gpu4pyscf.py (wB97M-D3(BJ)/def2-TZVPPD, energy and gradient)
for the pilot's conformers, written for Barnard, and a Slurm array that runs them in the horse workspace.

    orca_inputs.py [--pilot pilot] [--all] [--nprocs 8]

By default only the first conformer of each fragment (one point per fragment is enough to check that both codes agree
on energies and forces; the GPU run covers the relative energies). --all writes every conformer.
Out: pilot/orca/<name>/point.inp + point.xyz, pilot/orca/list.txt and pilot/barnard_orca.sbatch.
ORCA's RIJCOSX (def2/J) and DefGrid3 differ from GPU4PySCF's density fitting and grid by far less than the differences
the pilot looks for (about 0.1 kcal/mol in a relative energy).
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
WS = "/data/horse/ws/muri912g-caps-ml"

SBATCH = """#!/bin/bash
# The ORCA half of the pilot on Barnard; inputs, scratch, results and logs all in the horse workspace.
#   cd {ws}/pilot && sbatch barnard_orca.sbatch
#SBATCH --job-name=caps-ml-orca
#SBATCH --nodes=1
#SBATCH --ntasks={nprocs}
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --array=1-{n}
#SBATCH --output={ws}/logs/orca-%A_%a.out

module purge >/dev/null 2>&1
module load release/24.04 GCC/12.3.0 OpenMPI/4.1.5 ORCA/6.0.1-avx2 >/dev/null 2>&1 || module load ORCA/6.0.1-avx2
d={ws}/pilot/orca/$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {ws}/pilot/orca/list.txt)
grep -q "ORCA TERMINATED NORMALLY" "$d/point.out" 2>/dev/null && exit 0
cd "$d" && /usr/bin/time -f "wall %e s" $(which orca) point.inp > point.out 2>&1
"""


def main():
    pilot = arg("--pilot", os.path.join(HERE, "pilot"))
    nprocs = int(arg("--nprocs", "8"))
    rows = json.load(open(os.path.join(pilot, "pilot.json")))
    names = []
    for r in rows:
        for f in (r["files"] if "--all" in sys.argv else r["files"][:1]):
            name = f[:-4]
            d = os.path.join(pilot, "orca", name)
            os.makedirs(d, exist_ok=True)
            lines = open(os.path.join(pilot, "xyz", f)).read().splitlines()
            meta = dict(w.split("=", 1) for w in lines[1].split() if "=" in w)
            open(os.path.join(d, "point.xyz"), "w").write("\n".join(lines) + "\n")
            open(os.path.join(d, "point.inp"), "w").write(
                f"! wB97M-D3BJ def2-TZVPPD def2/J RIJCOSX DefGrid3 TightSCF EnGrad\n"
                f"%pal nprocs {nprocs} end\n%maxcore 3500\n"
                f"* xyzfile {meta.get('charge', 0)} {meta.get('mult', 1)} point.xyz\n")
            names.append(name)
    open(os.path.join(pilot, "orca", "list.txt"), "w").write("\n".join(names) + "\n")
    open(os.path.join(pilot, "barnard_orca.sbatch"), "w").write(SBATCH.format(ws=WS, nprocs=nprocs, n=len(names)))
    print(f"{len(names)} ORCA points -> {pilot}/orca")


if __name__ == "__main__":
    main()
