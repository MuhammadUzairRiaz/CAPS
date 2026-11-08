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
LMP = os.environ.get("LMP", shutil.which("lmp") or os.path.expanduser("~/lammps/build-caps/lmp"))
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
      [{"name": "toluene", "molecules": "1-12", "forcefield": "oplsaa2024"}, {"name": "cyclohexane", "molecules": "13-24", "forcefield": "gaff"}])

# D: 1-4 scalings that differ, with a four-site water
dd = caps.pack(molecules=[("Cc1ccccc1", 8), ("C1CCCCC1", 8), (os.path.join(samples, "water.pdb"), 60)], box=22.0, tolerance=2.0, seed=6)
dd.edit(op="water_model", model="tip4p2005")
check("D toluene OPLS-AA + cyclohexane GAFF + water TIP4P/2005", dd,
      [{"name": "toluene", "molecules": "1-8", "forcefield": "oplsaa2024"}, {"name": "cyclohexane", "molecules": "9-16", "forcefield": "gaff"},
       {"name": "water", "molecules": "water", "water": "tip4p2005"}])


# Force fields that change the atoms by group: massless sites on rigid molecules, written for LAMMPS as atoms of a rigid body
def gas_cell(gases):
    d = caps.polymer(r"[*]C/C=C\C[*]", dp=10, chains=2, density=0.6, seed=1)
    for smi, n in gases:
        d.insert(smi, n, tolerance=2.0, seed=1)
    return d


PB = {"name": "PB", "molecules": "1-2", "forcefield": "oplsaa2024"}
TR = "trappe-ua-umn"
check("E polybutadiene OPLS-AA + O2 TraPPE (a massless centre site)", gas_cell([("O=O", 3)]),
      [PB, {"name": "O2", "molecules": "3-5", "forcefield": TR}], eps_rule="geometric", sigma_rule="geometric")
check("F polybutadiene OPLS-AA + N2 TraPPE", gas_cell([("N#N", 3)]),
      [PB, {"name": "N2", "molecules": "3-5", "forcefield": TR}], eps_rule="geometric", sigma_rule="geometric")
check("G polybutadiene OPLS-AA + O2 + N2 + CO2, each TraPPE", gas_cell([("O=O", 2), ("N#N", 2), ("O=C=O", 2)]),
      [PB, {"name": "O2", "molecules": "3-4", "forcefield": TR}, {"name": "N2", "molecules": "5-6", "forcefield": TR},
       {"name": "CO2", "molecules": "7-8", "forcefield": TR}], eps_rule="geometric", sigma_rule="geometric")


def rigid_nve(label, d, groups, steps=4000, **kw):
    """The exported input with fix nve and fix rigid/small for fix nvt and fix rigid/nvt/small: etotal must stay put (the
    rigid thermostat's scalar repeats the bodies' kinetic energy, so it cannot be added to etotal), and the bodies rigid."""
    d.field.assign_groups(groups, **kw)
    folder = tempfile.mkdtemp(prefix="caps_rigid_nve_")
    d.export_engines(folder, stem="s", run="nvt", steps=steps, temperature=300, timestep=0.5, thermo_every=200, dump_every=500, kspace_accuracy=1e-5, cutoff=9, gromacs=False)
    text = open(os.path.join(folder, "s.in")).read()
    text = re.sub(r"rigid/nvt/small molecule temp \S+ \S+ \S+", "rigid/small molecule", text)
    text = re.sub(r"(fix\s+integrate mobile) nvt temp \S+ \S+ \S+", r"\1 nve", text)
    open(os.path.join(folder, "nve.in"), "w").write(text)
    subprocess.run([LMP, "-in", "nve.in", "-log", "nve.log"], cwd=folder, capture_output=True, text=True)
    log = open(os.path.join(folder, "nve.log")).read()
    if "ERROR" in log:
        print(f"== {label}\n   LAMMPS:", [l for l in log.splitlines() if "ERROR" in l][-1])
        return
    rows = [l.split() for l in log.splitlines() if re.match(r"^\s*\d+\s+-?\d", l) and len(l.split()) == 8]
    rows = rows[max(i for i, r in enumerate(rows) if r[0] == "0"):]
    e = [float(r[5]) for r in rows]
    print(f"== {label}\n   NVE etotal {e[0]:.3f} → {e[-1]:.3f}, range {max(e) - min(e):.3f} kcal/mol over {rows[-1][0]} steps")


rigid_nve("H polybutadiene + O2, NVE", gas_cell([("O=O", 3)]), [PB, {"name": "O2", "molecules": "3-5", "forcefield": TR}], eps_rule="geometric", sigma_rule="geometric")
rigid_nve("I polybutadiene + N2, NVE", gas_cell([("N#N", 3)]), [PB, {"name": "N2", "molecules": "3-5", "forcefield": TR}], eps_rule="geometric", sigma_rule="geometric")
