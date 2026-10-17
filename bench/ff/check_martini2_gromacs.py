#!/usr/bin/env python3
"""MARTINI 2 libraries in GROMACS: uncharged bead boxes built from the libraries' templates, CAPS's energy against GROMACS's
(grompp, then mdrun -rerun of the exported files; the shifted Lennard-Jones is GROMACS's own shift function). Charged MARTINI
2 systems are refused for GROMACS (the shifted Coulomb has no Verlet form) and are checked in LAMMPS
(bench/ff/check_data_lammps.py).

    PYTHONPATH=data/python CAPS_LIB=build/capi/libcaps.dylib python3 bench/ff/check_martini2_gromacs.py
"""
import os, re, shutil, subprocess, sys, tempfile

import caps

GMX = os.environ.get("GMX", shutil.which("gmx") or "gmx")
KJ = 4.184
CASES = [("martini", [("HD", 20), ("OCO", 20), ("CHEX", 20), ("W", 150), ("WF", 15)], 42.0),
         ("drymartini", [("HD", 30), ("OCO", 30), ("CHEX", 20)], 36.0),
         ("martini-polymers", [("PEO36", 2), ("PS10", 2), ("W", 120)], 44.0)]
work = tempfile.mkdtemp(prefix="caps_martini2_gmx_")
fails = 0
for ff, parts, edge in CASES:
    out = os.path.join(work, ff)
    os.makedirs(out, exist_ok=True)
    mols = []
    for name, n in parts:
        one = caps.build.beads(name, forcefield=ff)
        p = os.path.join(out, name.replace("+", "p").replace("-", "m") + ".data")
        one.save(p)
        mols.append((p, n))
    d = caps.pack(molecules=mols, box=edge, tolerance=3.0, seed=7)
    rep = d.field.assign(ff, charges="keep")
    d.relax(ftol=5.0, max_iterations=200)
    ce = d.energy()["total"]
    d.export_engines(out, stem="m", lammps=False, gromacs=True, run="check", cutoff=0)
    run = lambda a, i=None: subprocess.run([GMX, "-quiet"] + a, cwd=out, input=i, capture_output=True, text=True)
    r = run(["grompp", "-f", "m.mdp", "-c", "m.gro", "-p", "m.top", "-o", "m.tpr", "-maxwarn", "10"])
    if r.returncode:
        print(f"{ff}: grompp failed: {[l for l in (r.stdout + r.stderr).splitlines() if l.strip()][-3]}")
        fails += 1
        continue
    run(["mdrun", "-s", "m.tpr", "-rerun", "m.gro", "-deffnm", "m", "-nt", "4"])
    e = run(["energy", "-f", "m.edr", "-o", "e.xvg"], "Potential\n\n").stdout
    m = re.search(r"^Potential\s+(-?[\d.]+(?:e[-+]?\d+)?)", e, re.M)
    ge = float(m.group(1)) / KJ
    rel = abs(ce - ge) / max(1.0, abs(ge))
    bad = rel > 5e-5
    fails += bad
    print(f"{ff:18s} {d.atoms} beads: GROMACS {ge:.4f} / CAPS {ce:.4f} kcal/mol, {rel:.1e} relative{'  DIFFERS' if bad else ''}  ({rep['complete']})")
print(f"{'all agree' if not fails else f'{fails} failed'} (work: {work})")
sys.exit(1 if fails else 0)
