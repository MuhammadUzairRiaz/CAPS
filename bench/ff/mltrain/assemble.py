#!/usr/bin/env python3
"""The training set, step 3: conformers and their DFT energies and forces in one file for MACE.

    assemble.py TRAINSET [--test 0.1] [--seed 1] [--out DIR]

Reads TRAINSET/manifest.json, TRAINSET/exclude.txt, TRAINSET/confs/<env>.xyz and TRAINSET/dft/<env>.jsonl and writes
DIR/train.xyz and DIR/test.xyz (extended xyz: REF_energy in eV, REF_forces in eV/A, and per frame env, kind, angle,
level and config_type) plus DIR/summary.json. Left out: excluded environments, SCFs that did not converge and frames
with no DFT line, conformers with a force component over MAX_FORCE eV/A, and conformers of a double-bond rotor
twisted more than PLANAR degrees from planar. The test set holds out whole environments (a fraction --test of them, chosen with --seed), so it
measures the model on chemistry it has not seen, not on a neighbouring conformer of a training one.
"""
import io, json, os, random, sys

import numpy as np

from ase import Atoms
from ase.io import read, write

MAX_FORCE = 20.0          # eV/A: a conformer with a larger force component is an extreme distortion MD never visits
PLANAR = 30.0             # degrees: how far a double-bond rotor may be twisted from planar and stay in the set
DOUBLE = {("C", "C"): 1.40, ("C", "N"): 1.34, ("N", "N"): 1.30, ("C", "S"): 1.68, ("C", "O"): 1.28, ("N", "O"): 1.28}
HARTREE_EV = 27.211386245988
BOHR_A = 0.529177210903
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


def double_rotor(ts, env, central):
    """is the fragment's rotor a double bond (shorter than the single bond between those elements)? The census
    cannot tell bond orders, and a C=C turned through 90 degrees breaks its pi bond: no MD visits it, and a closed-shell
    DFT describes it badly"""
    el, pos, sec = [], [], None
    for line in open(os.path.join(ts, "fragments", env + ".mol2")):
        if line.startswith("@<TRIPOS>"):
            sec = line.strip()[9:]
            continue
        w = line.split()
        if sec == "ATOM" and len(w) >= 6:
            el.append(w[5])
            pos.append(np.array([float(x) for x in w[2:5]]))
    b, c = central[0] - 1, central[1] - 1
    key = tuple(sorted((el[b], el[c])))
    return key in DOUBLE and np.linalg.norm(pos[b] - pos[c]) < DOUBLE[key]


def main():
    ts = sys.argv[1]
    out = arg("--out", os.path.join(ts, "mace"))
    os.makedirs(out, exist_ok=True)
    manifest = json.load(open(os.path.join(ts, "manifest.json")))
    excl = {l.split()[0] for l in open(os.path.join(ts, "exclude.txt")) if l.strip() and not l.startswith("#")}
    envs = [e["env"] for e in manifest if e["env"] not in excl]
    central = {e["env"]: e["central_atoms"] for e in manifest}
    rng = random.Random(int(arg("--seed", "1")))
    test = set(rng.sample(sorted(envs), max(1, int(len(envs) * float(arg("--test", "0.1"))))))
    frames = {"train": [], "test": []}
    stats = {"environments": 0, "frames": 0, "unconverged": 0, "no_dft": 0, "no_df_level": 0, "missing_env": [],
             "large_force": 0, "twisted_double_bond": 0, "double_bond_rotors": []}
    for env in envs:
        src, dft = os.path.join(ts, "confs", env + ".xyz"), os.path.join(ts, "dft", env + ".jsonl")
        if not (os.path.exists(src) and os.path.exists(dft)):
            stats["missing_env"].append(env)
            continue
        text = open(src).read().replace("# done\n", "")
        confs = read(io.StringIO(text), index=":", format="extxyz") if text.strip() else []
        lines = {}
        for line in open(dft):
            try:
                d = json.loads(line)
                lines[d["frame"]] = d
            except (ValueError, KeyError):
                pass
        double = double_rotor(ts, env, central[env])
        if double:
            stats["double_bond_rotors"].append(env)
        n = 0
        for k, at in enumerate(confs):
            d = lines.get(k)
            if d is None:
                stats["no_dft"] += 1
                continue
            if not d["converged"]:
                stats["unconverged"] += 1
                continue
            if double:
                phi = float(at.info.get("angle", 0)) % 360
                if min(phi, abs(phi - 180), 360 - phi) > PLANAR:
                    stats["twisted_double_bond"] += 1
                    continue
            if np.abs(d["forces_Eh_per_bohr"]).max() * HARTREE_EV / BOHR_A > MAX_FORCE:
                stats["large_force"] += 1
                continue
            a = Atoms(at.get_chemical_symbols(), positions=at.positions)
            a.info.update({"REF_energy": d["energy_Eh"] * HARTREE_EV, "env": env, "kind": d.get("kind"),
                           "angle": d.get("angle"), "level": d["level"], "config_type": "Default"})
            a.arrays["REF_forces"] = np.array(d["forces_Eh_per_bohr"]) * (HARTREE_EV / BOHR_A)
            stats["no_df_level"] += "no DF" in d["level"]
            frames["test" if env in test else "train"].append(a)
            n += 1
        if n:
            stats["environments"] += 1
            stats["frames"] += n
    for part, fr in frames.items():
        write(os.path.join(out, f"{part}.xyz"), fr, format="extxyz")
    stats.update({"train_frames": len(frames["train"]), "test_frames": len(frames["test"]),
                  "test_environments": sorted(test), "excluded": len(excl)})
    json.dump(stats, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print({k: (len(v) if isinstance(v, list) else v) for k, v in stats.items()})


if __name__ == "__main__":
    main()
