#!/usr/bin/env python3
"""The pilot's verdict: GPU4PySCF (Alpha) against ORCA (Barnard) at wB97M-D3(BJ)/def2-TZVPPD, and what a conformer costs.

    compare.py [--pilot pilot]      (after copying back pilot/gpu/*.json and pilot/orca/*/point.out, point.engrad)

For every point both codes ran: the energy difference (kcal/mol; one code's absolute energy differs from the other's by
a near-constant offset from density fitting and grids, so the spread of the difference over the points is what
matters, and for the relative energies of a fragment's conformers only the spread does), and the force difference
(RMS and largest component, kcal/mol/A). Then the GPU's relative conformer energies per fragment and the time per point
on each code, and the cost of the full training set extrapolated from it. Writes pilot/compare.md.
"""
import glob, json, math, os, re, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
HARTREE = 627.509474           # kcal/mol
FORCE = HARTREE / 0.529177211  # Eh/bohr -> kcal/mol/A


def orca(d):
    out = os.path.join(d, "point.out")
    if not os.path.exists(out) or "ORCA TERMINATED NORMALLY" not in open(out).read():
        return None
    text = open(out).read()
    e = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(\S+)", text)[-1])
    wall = re.findall(r"wall (\S+) s", text)
    lines = open(os.path.join(d, "point.engrad")).read().splitlines()
    vals = [float(x) for x in lines if re.match(r"^\s*-?\d+\.\d+(E[-+]\d+)?\s*$", x, re.I)]
    n = int([x for x in lines if x.strip() and not x.startswith("#")][0])
    grad = vals[1:1 + 3 * n]   # after the energy
    return {"energy_Eh": e, "forces": [[-grad[3 * i + k] for k in range(3)] for i in range(n)],
            "wall_s": float(wall[-1]) if wall else None}


def main():
    pilot = arg("--pilot", os.path.join(HERE, "pilot"))
    rows = json.load(open(os.path.join(pilot, "pilot.json")))
    gpu = {os.path.basename(p)[:-5]: json.load(open(p)) for p in glob.glob(os.path.join(pilot, "gpu", "*.json"))}
    L = ["# Pilot: GPU4PySCF (Alpha, A100) vs ORCA 6.0.1 (Barnard, 8 cores)", "",
         "Level: wB97M-D3(BJ)/def2-TZVPPD, energy and gradient.", "",
         "| fragment | atoms | dE GPU-ORCA (kcal/mol) | force RMS diff | force max diff (kcal/mol/A) | GPU s | ORCA s |",
         "|---|---|---|---|---|---|---|"]
    diffs, frms, fmax, tg, to = [], [], [], [], []
    for name in sorted(gpu):
        o = orca(os.path.join(pilot, "orca", name))
        g = gpu[name]
        tg.append(g["scf_s"] + g["grad_s"])
        if not o:
            continue
        de = (g["energy_Eh"] - o["energy_Eh"]) * HARTREE
        d = [(a - b) * FORCE for ra, rb in zip(g["forces_Eh_per_bohr"], o["forces"]) for a, b in zip(ra, rb)]
        rms, mx = math.sqrt(sum(x * x for x in d) / len(d)), max(abs(x) for x in d)
        diffs.append(de), frms.append(rms), fmax.append(mx)
        if o["wall_s"]:
            to.append(o["wall_s"])
        L.append(f"| {name} | {g['atoms']} | {de:+.3f} | {rms:.3f} | {mx:.3f} | {tg[-1]:.0f} | {o['wall_s'] or '-'} |")
    if diffs:
        L += ["", f"Energy offset GPU - ORCA: mean {statistics.mean(diffs):+.3f}, spread (SD) "
                  f"{statistics.pstdev(diffs):.3f} kcal/mol over {len(diffs)} points.",
              f"Forces: RMS difference median {statistics.median(frms):.3f}, largest component {max(fmax):.3f} kcal/mol/A."]
    L += ["", "## Relative conformer energies (GPU, kcal/mol)", ""]
    for r in rows:
        es = [gpu[f[:-4]]["energy_Eh"] for f in r["files"] if f[:-4] in gpu]
        if es:
            L.append(f"- {r['env']} ({r['status']}, {r['example'][:40]}): "
                     + ", ".join(f"{(e - min(es)) * HARTREE:.2f}" for e in es))
    if tg:
        per = statistics.mean(tg)
        L += ["", "## Cost", "",
              f"- GPU: {per:.0f} s per conformer on average (SCF + gradient, one A100), {len(tg)} conformers.",
              f"- ORCA: {statistics.mean(to):.0f} s per point on 8 cores." if to else "- ORCA: no finished points yet.",
              f"- Full training set, 55,000 conformers: about {55000 * per / 3600:.0f} GPU-hours, "
              f"{55000 * per / 3600 / 32:.0f} h of wall time on 32 A100s."]
    open(os.path.join(pilot, "compare.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
