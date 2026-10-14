#!/usr/bin/env python3
"""The training set of CAPS's own ML potential, step 0: the census fragments split into GPU tasks of similar cost.

    make_tasks.py [--census bench/ff/census/oplsaa2024] [--tasks 40] [--out DIR]

Every sound census fragment (no "problem") goes to one task; tasks are filled largest fragment first, each into the
task with the least estimated cost so far (cost ~ atoms^2.5, from the pilot's timings), so the array finishes together.
Out: DIR/fragments/<env>.mol2, DIR/manifest.json (env, central bond, atoms, status, example) and DIR/tasks/NN.txt.
On Alpha each task runs gen_conformers.py (AIMNet2 relaxed scan + distorted copies) and then dft_confs.py (GPU4PySCF).
"""
import json, os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


def main():
    census = arg("--census", os.path.join(ROOT, "bench", "ff", "census", "oplsaa2024"))
    out = arg("--out", os.path.join(HERE, "trainset"))
    ntask = int(arg("--tasks", "40"))
    envs = [e for e in json.load(open(os.path.join(census, "census.json")))["environments"] if not e.get("problem")]
    os.makedirs(os.path.join(out, "fragments"), exist_ok=True)
    os.makedirs(os.path.join(out, "tasks"), exist_ok=True)
    manifest = []
    for e in envs:
        shutil.copy(os.path.join(census, "fragments", e["env"] + ".mol2"), os.path.join(out, "fragments"))
        manifest.append({k: e[k] for k in ("env", "central_atoms", "fragment_atoms", "fragment_heavy", "status",
                                           "classes", "example", "elements")})
    json.dump(manifest, open(os.path.join(out, "manifest.json"), "w"), indent=1)
    load, tasks = [0.0] * ntask, [[] for _ in range(ntask)]
    for e in sorted(manifest, key=lambda e: -e["fragment_atoms"]):
        k = min(range(ntask), key=load.__getitem__)
        tasks[k].append(e["env"])
        load[k] += e["fragment_atoms"] ** 2.5
    for k, t in enumerate(tasks):
        open(os.path.join(out, "tasks", f"{k:02d}.txt"), "w").write("\n".join(t) + "\n")
    print(f"{len(manifest)} fragments in {ntask} tasks ({min(map(len, tasks))}-{max(map(len, tasks))} each), "
          f"largest {max(e['fragment_atoms'] for e in manifest)} atoms -> {out}")


if __name__ == "__main__":
    main()
