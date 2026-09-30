"""Cells with a force field per group of molecules, exported and checked against LAMMPS and GROMACS:

  A  toluene (GAFF) + cyclohexane (GAFF2) + water (TIP3P), packed together
  B  polystyrene (PCFF, class II 9-6) + water (TIP4P/2005): the cross pairs as 12-6 with the 9-6 site's ε and r_min
  C  toluene (OPLS-AA) + cyclohexane (GAFF): 1-4 scalings differ (0.5/0.5 against 0.5/0.8333) — refused by default

For each: CAPS's single-point energy (PME) against LAMMPS's run 0 with the exported files; GROMACS runs the exported
topology (grompp, mdrun -rerun) and its bonded terms and Coulomb sum are compared with CAPS's.

    PYTHONPATH=data/python CAPS_LIB=build/capi/libcaps.dylib python3 bench/ff/check_mixed_forcefields.py
"""
import ctypes as C
import os
import re
import shutil
import subprocess
import sys
import tempfile

import caps

L = caps.library()
L.caps_set_electrostatics.argtypes = [C.c_int32, C.c_double, C.c_double, C.c_int32]
L.caps_set_electrostatics(1, 1e-6, 0.6, 6)
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-caps/lmp"))
GMX = shutil.which("gmx")
samples = os.environ.get("CAPS_SAMPLES", "samples")
KJ = 4.184


def lammps_total(folder):
    subprocess.run([LMP, "-in", "s.in", "-log", "log.lammps"], cwd=folder, capture_output=True, text=True)
    log = open(os.path.join(folder, "log.lammps")).read()
    rows = re.findall(r"^\s*0\s+(-?[\d.]+)", log, re.M)
    if not rows:
        err = [l for l in log.splitlines() if "ERROR" in l]
        raise RuntimeError(err[-1] if err else "no thermo row")
    return float(rows[-1])


def gromacs_terms(folder):
    if not GMX:
        return None
    run = lambda a, inp=None: subprocess.run([GMX] + a, cwd=folder, capture_output=True, text=True, input=inp)
    r = run(["grompp", "-f", "s.mdp", "-c", "s.gro", "-p", "s.top", "-o", "s.tpr", "-maxwarn", "10"])
    if r.returncode:
        return "grompp: " + [l for l in (r.stdout + r.stderr).splitlines() if l.strip()][-2]
    r = run(["mdrun", "-s", "s.tpr", "-rerun", "s.gro", "-deffnm", "s", "-nt", "2"])
    if r.returncode:
        return "mdrun failed"
    r = run(["energy", "-f", "s.edr", "-o", "e.xvg"], "\n".join(str(i) for i in range(1, 30)) + "\n\n")
    vals = {}
    for line in r.stdout.splitlines():
        m = re.match(r"^(\S.*?)\s+(-?[\d.]+(?:e[-+]?\d+)?)\s+--.*\(kJ/mol\)", line)
        if m:
            vals[m.group(1).strip()] = float(m.group(2)) / KJ
    return vals


def check(label, d, groups, **kw):
    print(f"== {label}")
    try:
        d.field.assign_groups(groups, **kw)
    except caps.CapsError as e:
        print("   refused:", str(e)[:300])
        return
    e = d.energy()
    folder = tempfile.mkdtemp(prefix="caps_mixed_")
    r = d.export_engines(folder, stem="s", run="check", kspace_accuracy=1e-7, cutoff=0)
    styles = [l.strip() for l in open(os.path.join(folder, "s.in")) if l.startswith(("pair_style", "kspace_style", "bond_style", "angle_style", "dihedral_style", "special_bonds", "group "))]
    print("   LAMMPS:", " | ".join(styles))
    try:
        lt = lammps_total(folder)
        print(f"   CAPS {e['total']:.4f}  LAMMPS {lt:.4f}  Δ {e['total'] - lt:+.4f} kcal/mol")
    except RuntimeError as ex:
        print("   LAMMPS:", ex)
    g = gromacs_terms(folder) if r.get("gromacs_error") in (None, "") else "not written: " + r["gromacs_error"]
    if isinstance(g, dict):
        bonded = sum(v for k, v in g.items() if k in ("Bond", "Angle", "Proper Dih.", "Ryckaert-Bell.", "Improper Dih.", "Per. Imp. Dih.", "U-B"))
        coul = sum(v for k, v in g.items() if k in ("Coulomb (SR)", "Coul. recip.", "Coulomb-14"))
        caps_bonded = e["bond"] + e["angle"] + e["dihedral"] + e["improper"]
        print(f"   GROMACS bonded {bonded:.4f} (CAPS {caps_bonded:.4f}) · Coulomb {coul:.4f} (CAPS {e['coulomb']:.4f}) kcal/mol")
    elif g is not None:
        print("   GROMACS:", g)


# A: three solvents, three force fields (one 1-4 family: AMBER's)
a = caps.pack(molecules=[("Cc1ccccc1", 12), ("C1CCCCC1", 12), (os.path.join(samples, "water.pdb"), 80)], box=24.0, tolerance=2.0, seed=3)
a.edit(op="water_model", model="tip3p")
check("A toluene GAFF + cyclohexane GAFF2 + water TIP3P", a,
      [{"name": "toluene", "molecules": "1-12", "forcefield": "gaff"}, {"name": "cyclohexane", "molecules": "13-24", "forcefield": "gaff2"},
       {"name": "water", "molecules": "water", "water": "tip3p"}])

# B: a class II polymer with a four-site water
ps = caps.polymer("*CC(*)c1ccccc1", dp=6, chains=1, seed=2)
b = caps.pack(molecules=[(ps, 2), (os.path.join(samples, "water.pdb"), 60)], box=24.0, tolerance=2.0, seed=5)
b.edit(op="water_model", model="tip4p2005")
check("B polystyrene PCFF + water TIP4P/2005 (cross 9-6 → 12-6 r_min)", b,
      [{"name": "PS", "molecules": "1-2", "forcefield": "pcff"}, {"name": "water", "molecules": "water", "water": "tip4p2005"}], cross96="rmin")

# C: 1-4 scalings that differ
c = caps.pack(molecules=[("Cc1ccccc1", 12), ("C1CCCCC1", 12)], box=22.0, tolerance=2.0, seed=4)
check("C toluene OPLS-AA + cyclohexane GAFF (each its own 1-4 scaling)", c,
      [{"name": "toluene", "molecules": "1-12", "forcefield": "opls2005"}, {"name": "cyclohexane", "molecules": "13-24", "forcefield": "gaff"}])

# D: 1-4 scalings that differ, with a four-site water
dd = caps.pack(molecules=[("Cc1ccccc1", 8), ("C1CCCCC1", 8), (os.path.join(samples, "water.pdb"), 60)], box=22.0, tolerance=2.0, seed=6)
dd.edit(op="water_model", model="tip4p2005")
check("D toluene OPLS-AA + cyclohexane GAFF + water TIP4P/2005", dd,
      [{"name": "toluene", "molecules": "1-8", "forcefield": "opls2005"}, {"name": "cyclohexane", "molecules": "9-16", "forcefield": "gaff"},
       {"name": "water", "molecules": "water", "water": "tip4p2005"}])
