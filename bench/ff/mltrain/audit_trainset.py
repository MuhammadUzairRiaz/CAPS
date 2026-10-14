#!/usr/bin/env python3
"""Audit of an assembled training set before it is trained on: are the data physically sensible?

    audit_trainset.py DIR      (DIR/train.xyz and DIR/test.xyz from assemble.py)

Checks every conformer and lists the worst cases of each by environment:
  barrier      the highest relaxed scan point of each environment above its lowest (kcal/mol); a torsion barrier is
               under about 25, far more means a scan through a broken bond or a clash
  composition  each conformer's energy against a least-squares fit on element counts (eV/atom): a wrong charge, spin
               or element stands out by far more than conformational energy
  relaxed      the largest force component on a relaxed scan point (only the held dihedral should carry force)
  geometry     the shortest interatomic distance, and bonds (within 1.25 x covalent radii) shorter or longer than
               any sound molecule has (X-H under 0.90 A, C-C under 1.15 A or over 1.70 A)
  duplicates   conformers with the same elements and coordinates
Writes DIR/audit.json and prints a summary; exits 1 when a check fails its limit.
"""
import collections, json, os, sys

import numpy as np
from ase.data import covalent_radii
from ase.io import iread

EV_KCAL = 23.0605478
LIMITS = {"barrier_kcal": 40.0, "composition_ev_per_atom": 0.5, "relaxed_force_ev_a": 5.0, "closest_a": 0.85}


def main():
    d = sys.argv[1]
    frames = []
    for part in ("train", "test"):
        for a in iread(os.path.join(d, f"{part}.xyz"), index=":", format="extxyz"):
            frames.append((part, a))
    els = sorted({s for _, a in frames for s in a.get_chemical_symbols()})
    X = np.array([[a.get_chemical_symbols().count(e) for e in els] for _, a in frames], float)
    E = np.array([a.info["REF_energy"] for _, a in frames])
    coef, *_ = np.linalg.lstsq(X, E, rcond=None)
    resid = (E - X @ coef) / X.sum(1)
    by_env = collections.defaultdict(list)
    worst = {"barrier": [], "composition": [], "relaxed": [], "geometry": [], "bonds": []}
    seen, dup = set(), 0
    for k, (part, a) in enumerate(frames):
        env = a.info["env"]
        by_env[env].append(k)
        f = np.abs(a.arrays["REF_forces"]).max()
        if a.info.get("kind") == "relaxed":
            worst["relaxed"].append((f, env, a.info.get("angle")))
        dm = a.get_all_distances()
        np.fill_diagonal(dm, 9)
        worst["geometry"].append((dm.min(), env, a.info.get("kind"), a.info.get("angle")))
        r = covalent_radii[a.numbers]
        sym = a.get_chemical_symbols()
        for i, j in zip(*np.where(np.triu(dm < 1.25 * (r[:, None] + r[None, :]), 1))):
            pair = "".join(sorted((sym[i], sym[j])))
            if ("H" in pair and dm[i, j] < 0.90) or (pair == "CC" and not 1.15 <= dm[i, j] <= 1.70):
                worst["bonds"].append((dm[i, j], pair, env, a.info.get("kind")))
        key = (tuple(sym), tuple(np.round(a.positions, 4).ravel()))
        dup += key in seen
        seen.add(key)
        worst["composition"].append((abs(resid[k]), env, a.info.get("kind")))
    for env, ks in by_env.items():
        rel = [E[k] for k in ks if frames[k][1].info.get("kind") == "relaxed"]
        if rel:
            worst["barrier"].append(((max(rel) - min(rel)) * EV_KCAL, env))
    for k in worst:
        worst[k].sort(key=lambda t: -t[0] if k != "geometry" else t[0])
    fails = {
        "barrier": sum(b > LIMITS["barrier_kcal"] for b, _ in worst["barrier"]),
        "composition": sum(c > LIMITS["composition_ev_per_atom"] for c, *_ in worst["composition"]),
        "relaxed": sum(f > LIMITS["relaxed_force_ev_a"] for f, *_ in worst["relaxed"]),
        "geometry": sum(g < LIMITS["closest_a"] for g, *_ in worst["geometry"]),
        "bonds": len(worst["bonds"]),
        "duplicates": dup,
    }
    fails = {k: int(v) for k, v in fails.items()}
    summary = {"frames": len(frames), "environments": len(by_env), "elements": els, "limits": LIMITS, "fails": fails,
               "worst": {k: [[round(float(t[0]), 3), *map(str, t[1:])] for t in v[:10]] for k, v in worst.items()},
               "barrier_percentiles_kcal": [round(float(np.percentile([b for b, _ in worst["barrier"]], p)), 1)
                                            for p in (50, 90, 99)],
               "composition_rmse_ev_per_atom": round(float(np.sqrt((resid ** 2).mean())), 4)}
    json.dump(summary, open(os.path.join(d, "audit.json"), "w"), indent=1)
    print(json.dumps({k: summary[k] for k in ("frames", "environments", "elements", "fails",
                                              "barrier_percentiles_kcal", "composition_rmse_ev_per_atom")}))
    for k, v in summary["worst"].items():
        print(f"worst {k}:", v[:5])
    sys.exit(1 if any(fails.values()) else 0)


if __name__ == "__main__":
    main()
