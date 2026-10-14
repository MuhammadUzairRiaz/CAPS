#!/usr/bin/env python3
"""AIMNet2 against CAPS's quantum-chemical torsion scans: can a machine-learned potential stand in for MP2 / DLPNO when
a torsion is missing (an automatic fit with no cluster run)?

    aimnet_torsion_benchmark.py PLAN.json [PLAN.json ...] [--model aimnet2] [--max-scans N] [--only CLASSES]
                                [--no-relax] [--json out.json]

PLAN.json: a collected caps.torsionfit plan (each point with its optimised positions and E_qm). For every scan:
  sp       AIMNet2's energy at each QM geometry (the energy surface alone)
  relaxed  AIMNet2's own relaxed scan: each QM geometry optimised by AIMNet2 with the scanned dihedral held (ASE
           FixInternals, BFGS to 0.01 eV/A), as an automatic fit would run it
and both profiles (relative to their lowest point; points more than 30 kcal/mol above the QM minimum left out, as
torsionfit's fit does) are compared with QM: RMS, largest deviation, and the barrier (max - min) error.
Runs in an environment with aimnet[ase] (MIT, isayevlab/aimnetcentral); the model downloads once to ~/.cache/aimnet.
"""
import json, math, sys, time

import numpy as np
from ase import Atoms
from ase.constraints import FixInternals
from ase.optimize import BFGS
from aimnet.calculators import AIMNet2ASE

HARTREE = 627.509474
EV = 23.0605478   # kcal/mol per eV
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
plans = [a for a in sys.argv[1:] if a.endswith(".json") and a != arg("--json")]
only = set(filter(None, (arg("--only", "") or "").split(",")))
relax = "--no-relax" not in sys.argv
calc = AIMNet2ASE(arg("--model", "aimnet2"))


def energy(at):
    at.calc = calc
    return at.get_potential_energy() * EV


def relaxed(at, quad, steps=300):
    a = at.copy()
    a.calc = calc
    phi = a.get_dihedral(*quad)
    a.set_constraint(FixInternals(dihedrals_deg=[[phi, list(quad)]]))
    opt = BFGS(a, logfile=None)
    conv = opt.run(fmax=0.01, steps=steps)
    return a.get_potential_energy() * EV, bool(conv)


def stats(ref, val):
    ref, val = np.array(ref), np.array(val)
    d = (val - val.min()) - (ref - ref.min())
    return {"rms": float(np.sqrt(np.mean(d ** 2))), "max": float(np.max(np.abs(d))),
            "barrier_err": float((val.max() - val.min()) - (ref.max() - ref.min())), "barrier_qm": float(ref.max() - ref.min())}


out = []
t0 = time.time()
for plan_path in plans:
    plan = json.load(open(plan_path))
    for s in plan["scans"]:
        key = "-".join(s["classes"])
        if only and key not in only:
            continue
        pts = sorted([p for p in s["points"] if "E_qm" in p and "positions" in p], key=lambda p: p["phi"])
        if len(pts) < 8:
            continue
        eqm = [p["E_qm"] * HARTREE for p in pts]
        keep = [i for i, e in enumerate(eqm) if e - min(eqm) <= 30.0]
        pts, eqm = [pts[i] for i in keep], [eqm[i] for i in keep]
        sp, rx, nconv = [], [], 0
        for p in pts:
            at = Atoms(symbols=s["elements"], positions=p["positions"])
            sp.append(energy(at))
            if relax:
                e, c = relaxed(at, s["quad"])
                rx.append(e)
                nconv += c
        rec = {"plan": plan_path, "scan": s["dir"], "classes": key, "level": plan.get("level", ""), "atoms": len(s["elements"]),
               "points": len(pts), "phi": [p["phi"] for p in pts], "E_qm": [e - min(eqm) for e in eqm],
               "E_sp": [e - min(sp) for e in sp], "sp": stats(eqm, sp)}
        if relax:
            rec.update({"E_relaxed": [e - min(rx) for e in rx], "relaxed": stats(eqm, rx), "converged": nconv})
        out.append(rec)
        print(f"{key:28s} {len(s['elements']):3d} atoms {len(pts):2d} pts  barrier {rec['sp']['barrier_qm']:5.2f}  "
              f"sp RMS {rec['sp']['rms']:.2f} max {rec['sp']['max']:.2f}"
              + (f"  relaxed RMS {rec['relaxed']['rms']:.2f} max {rec['relaxed']['max']:.2f} ({nconv}/{len(pts)} conv)" if relax else "")
              + f"  [{time.time() - t0:.0f} s]", flush=True)
        if arg("--max-scans") and len(out) >= int(arg("--max-scans")):
            break
    if arg("--max-scans") and len(out) >= int(arg("--max-scans")):
        break

if out:
    for k in (["sp", "relaxed"] if relax else ["sp"]):
        r = np.array([o[k]["rms"] for o in out])
        print(f"\n{k}: {len(out)} scans, RMS median {np.median(r):.2f}, mean {r.mean():.2f}, 90th pct {np.percentile(r, 90):.2f}, "
              f"worst {r.max():.2f} kcal/mol; scans under 0.5: {np.mean(r < 0.5) * 100:.0f} %, under 1.0: {np.mean(r < 1.0) * 100:.0f} %")
if arg("--json"):
    json.dump({"model": arg("--model", "aimnet2"), "scans": out}, open(arg("--json"), "w"), indent=1)
