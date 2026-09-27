#!/usr/bin/env python3
"""CAPS's moltemplate export against its own LAMMPS files: a grown polystyrene cell typed with each force field, written as
LAMMPS data + input and as a moltemplate system (caps ff apply … --moltemplate); moltemplate.sh -overlay-all builds the
LAMMPS files from the .lt, and LAMMPS (LMP) evaluates both: every energy term must agree.

usage: python3 bench/ff/check_moltemplate.py [--ff id,id,…] [--keep DIR]
"""
import os, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
MT = os.environ.get("MOLTEMPLATE", os.path.expanduser("~/moltemplate/moltemplate/scripts/moltemplate.sh"))
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FFS = arg("--ff", "pcff-frc,compass-frc,cvff-frc,oplsaa2024-moltemplate,opls2005,gaff-amber25,dreiding,charmm36-prot").split(",")
WORK = arg("--keep", "") or tempfile.mkdtemp(prefix="caps_lt_")
THERMO = "thermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong\nthermo_modify format float %.8f\nrun 0\n"


def energies(workdir, infile):
    r = subprocess.run([LMP, "-in", infile, "-log", "none"], cwd=workdir, capture_output=True, text=True, timeout=600)
    lines = r.stdout.split("\n")
    for i, l in enumerate(lines):
        if l.split()[:1] == ["Step"]:
            return [float(x) for x in lines[i + 1].split()[1:]]
    raise RuntimeError("LAMMPS: " + " | ".join((r.stdout + r.stderr).strip().split("\n")[-3:]))


def main():
    os.makedirs(WORK, exist_ok=True)
    cell = os.path.join(WORK, "ps.data")
    subprocess.run([CAPS, "grow", "-o", cell, "--chains", "2", "--dp", "4", "--density", "0.3", "--seed", "2"], check=True, capture_output=True)
    ok = 0
    for ff in FFS:
        d = os.path.join(WORK, ff)
        os.makedirs(d, exist_ok=True)
        path = os.path.join(ROOT, "data", "forcefields", ff + ".json")
        r = subprocess.run([CAPS, "ff", "apply", cell, "--ff", path, "-o", "caps.data", "--lammps-input", "caps.in", "--moltemplate", "system.lt"],
                           cwd=d, capture_output=True, text=True)
        if not os.path.exists(os.path.join(d, "system.lt")):
            print(f"{ff:26s} not written: {(r.stdout + r.stderr).strip().splitlines()[-1:]}")
            continue
        m = subprocess.run([MT, "-overlay-all", "system.lt"], cwd=d, capture_output=True, text=True)
        if m.returncode != 0 or not os.path.exists(os.path.join(d, "system.data")):
            print(f"{ff:26s} moltemplate failed: {(m.stdout + m.stderr).strip().splitlines()[-2:]}")
            continue
        src = [l for l in open(os.path.join(d, "caps.in")) if not l.startswith(("thermo", "run "))]
        open(os.path.join(d, "c.in"), "w").write("".join(src) + THERMO)
        open(os.path.join(d, "m.in"), "w").write("include system.in.init\nread_data system.data\ninclude system.in.settings\n"
                                                 "neighbor 2 bin\nneigh_modify delay 0 every 1 check yes\n" + THERMO)
        try:
            a, b = energies(d, "c.in"), energies(d, "m.in")
        except Exception as e:
            print(f"{ff:26s} {e}")
            continue
        worst = max(abs(x - y) / max(1.0, abs(x)) for x, y in zip(a, b))
        good = worst < 1e-8
        ok += good
        print(f"{ff:26s} {'same energies' if good else 'differs'}  (largest relative difference {worst:.1e}; pe {a[0]:.6f} / {b[0]:.6f})")
    print(f"\n{ok} of {len(FFS)} force fields give the same energies through moltemplate; work {WORK}")


if __name__ == "__main__":
    main()
