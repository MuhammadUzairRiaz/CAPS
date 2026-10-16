#!/usr/bin/env python3
"""The numbers behind the ML figures that are not already in a results file.

    collect.py parity TEST.xyz MODEL.model OUT.npz
    collect.py profiles DLPNO_WORKDIR MODEL.model OUT.json [CLASS ...]

parity    the model's energy and forces on every frame of the held-out test set, beside the DFT reference (REF_energy,
          REF_forces).
profiles  for each class, its DLPNO-CCSD(T) scan (the library workdir: plan.json with E_qm and the optimised points):
          the DLPNO energy, the model's energy at the same geometry, and CAPS's MM energy of OPLS-AA 2024 with the
          fitted torsion (oplsaa2024-dlpno.json overlay) and without it (the torsion missing, at zero), all in kcal/mol
          relative to each curve's lowest point.
MODEL is the uncompiled MACE model (the compiled one in data/ml is for the C++ engine); CAPS_BIN the caps binary.
"""
import json, os, re, subprocess, sys

import numpy as np
from ase import Atoms
from ase.io import read

EV_KCAL = 23.0605478
HARTREE_KCAL = 627.509474
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))


def calculator(model):
    from mace.calculators import MACECalculator
    return MACECalculator(model_paths=model, device="cpu", default_dtype="float64")


def parity(test, model, out):
    calc = calculator(model)
    env, n, e_ref, e_ml, f_ref, f_ml = [], [], [], [], [], []
    for a in read(test, index=":", format="extxyz"):
        env.append(str(a.info["env"]))
        n.append(len(a))
        e_ref.append(a.info["REF_energy"])
        f_ref.append(a.arrays["REF_forces"])
        a.calc = calc
        e_ml.append(a.get_potential_energy())
        f_ml.append(a.get_forces())
    np.savez_compressed(out, env=np.array(env), n=np.array(n), e_ref=np.array(e_ref), e_ml=np.array(e_ml),
                        f_ref=np.concatenate(f_ref), f_ml=np.concatenate(f_ml))
    print(f"{len(n)} frames -> {out}")


def mm_energy(mol2, overlay=None):
    cmd = [os.environ.get("CAPS_BIN", os.path.join(ROOT, "build-ml", "cli", "caps")), "ff", "apply", mol2, "--ff",
           os.path.join(ROOT, "data", "forcefields", "oplsaa2024.json"), "--charges", "types", "--allow-missing",
           "-o", mol2 + (".fit.data" if overlay else ".bare.data")]
    if overlay:
        cmd += ["--overlay", overlay]
    r = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return float(re.search(r"total (\S+)", r.stdout).group(1))


def profiles(work, model, out, wanted):
    calc = calculator(model)
    overlay = os.path.join(ROOT, "data", "forcefields", "oplsaa2024-dlpno.json")
    plan = json.load(open(os.path.join(work, "plan.json")))
    res = {}
    for s in plan["scans"]:
        if s["dir"] not in wanted:
            continue
        pts = sorted((p for p in s["points"] if "E_qm" in p), key=lambda p: p["phi"])
        e0 = min(p["E_qm"] for p in pts)
        pts = [p for p in pts if (p["E_qm"] - e0) * HARTREE_KCAL <= 30.0]   # as the fit: a broken point is left out
        phi, qm, ml, fit, bare = [], [], [], [], []
        try:
            for p in pts:
                a = Atoms(s["elements"], positions=p["positions"])
                a.calc = calc
                m2 = os.path.join(work, p["dir"], "opt.mol2")
                fit.append(mm_energy(m2, overlay))
                bare.append(mm_energy(m2))
                phi.append(p["phi"])
                qm.append(p["E_qm"] * HARTREE_KCAL)
                ml.append(a.get_potential_energy() * EV_KCAL)
        except subprocess.CalledProcessError as ex:   # a scan whose optimisation broke the fragment: CAPS cannot type it
            print(f"{s['dir']}: skipped, caps ff apply failed at {p['phi']:+.0f} ({ex.stderr.strip().splitlines()[-1]})")
            continue
        rel = lambda v: [round(x - min(v), 6) for x in v]
        res[s["dir"]] = {"phi": phi, "dlpno": rel(qm), "caps_ml": rel(ml), "opls_fitted": rel(fit),
                         "opls_without": rel(bare), "atoms": len(s["elements"])}
        print(f"{s['dir']}: {len(phi)} points")
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "parity":
        parity(*sys.argv[2:5])
    else:
        profiles(sys.argv[2], sys.argv[3], sys.argv[4], set(sys.argv[5:]))
