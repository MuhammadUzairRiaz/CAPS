#!/usr/bin/env python3
"""How well do ML potentials reproduce the DFT torsion profiles of fragments they never saw?

    eval_torsions.py TEST.xyz MODEL.model [--out DIR] [--aimnet] [--device cuda|cpu]

For each environment in TEST.xyz (held out of training whole), the relaxed scan points' energies relative to the
environment's lowest DFT point, against the same from the CAPS MACE model and (with --aimnet) AIMNet2. Per
environment: RMSE of the relative energies and the error in the barrier (highest minus lowest point), in kcal/mol.
Writes DIR/torsions.json (per environment and summary) and prints the summary: median and 90th percentile of the
per-environment RMSE, and how many environments are within 0.5 and 1 kcal/mol.
"""
import collections, json, os, sys

import numpy as np
from ase.io import read

EV_KCAL = 23.0605478
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


def profile(frames, calc):
    out = []
    for a in frames:
        b = a.copy()
        b.calc = calc
        out.append(b.get_potential_energy())
    return np.array(out)


def stats(per_env, key):
    r = np.array([v[key]["rmse"] for v in per_env.values()])
    b = np.array([abs(v[key]["barrier_error"]) for v in per_env.values()])
    return {"environments": len(r), "rmse_median": round(float(np.median(r)), 3),
            "rmse_p90": round(float(np.percentile(r, 90)), 3), "rmse_max": round(float(r.max()), 3),
            "within_0.5": int((r <= 0.5).sum()), "within_1.0": int((r <= 1.0).sum()),
            "barrier_error_median": round(float(np.median(b)), 3), "barrier_error_max": round(float(b.max()), 3)}


def main():
    test, model = sys.argv[1], sys.argv[2]
    out = arg("--out", os.path.dirname(os.path.abspath(test)))
    from mace.calculators import MACECalculator
    calcs = {"caps_mace": MACECalculator(model_paths=model, device=arg("--device", "cuda"), default_dtype="float64")}
    if "--aimnet" in sys.argv:
        from aimnet.calculators import AIMNet2ASE
        calcs["aimnet2"] = AIMNet2ASE("aimnet2")
    by_env = collections.defaultdict(list)
    for a in read(test, index=":", format="extxyz"):
        if a.info.get("kind") == "relaxed":
            by_env[str(a.info["env"])].append(a)
    per_env = {}
    for env, frames in by_env.items():
        if len(frames) < 3:
            continue
        frames.sort(key=lambda a: float(a.info.get("angle", 0)))
        ref = np.array([a.info["REF_energy"] for a in frames])
        ref = (ref - ref.min()) * EV_KCAL
        rec = {"points": len(frames), "dft_barrier": round(float(ref.max()), 3)}
        for name, calc in calcs.items():
            e = profile(frames, calc)
            e = (e - e[np.argmin(ref)]) * EV_KCAL    # zero at the DFT minimum
            rec[name] = {"rmse": float(np.sqrt(((e - ref) ** 2).mean())),
                         "barrier_error": float((e.max() - e.min()) - ref.max())}
        per_env[env] = rec
    summary = {name: stats(per_env, name) for name in calcs}
    json.dump({"summary": summary, "environments": per_env}, open(os.path.join(out, "torsions.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
