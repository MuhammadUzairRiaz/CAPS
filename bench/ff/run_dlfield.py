#!/usr/bin/env python3
"""Run DL_FIELD (the oracle for converted DL_FIELD force fields) on one structure without touching the
user's DL_FIELD installation: the installation is copied once to build/dlfield_work and run there.

usage: run_dlfield.py FF STRUCTURE OUTDIR [--udff FILE]
  FF         DL_FIELD force-field keyword (pcff, cvff, compass, opls2005, amber, charmm36_cgenff, dreiding, ...)
  STRUCTURE  .xyz / .pdb / .mol2 (xyz: element typing through DL_F notation)
  OUTDIR     receives lammps.in, lammps1.data, dl_field.output, dlf_notation.output
Environment: DLFIELD (default ~/project/dl_f_4.13).
"""
import os, shutil, subprocess, sys

src = os.environ.get("DLFIELD", "~/project/dl_f_4.13")
here = os.path.dirname(os.path.abspath(__file__))
work = os.path.join(here, "..", "..", "build", "dlfield_work")
ff, structure, outdir = sys.argv[1], os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3])
udff = sys.argv[sys.argv.index("--udff") + 1] if "--udff" in sys.argv else "none"

if not os.path.exists(os.path.join(work, "dl_field")):
    shutil.copytree(src, work, ignore=shutil.ignore_patterns("output", "Examples", "source", "*.zip"))
os.makedirs(os.path.join(work, "out"), exist_ok=True)
for f in os.listdir(os.path.join(work, "out")):
    p = os.path.join(work, "out", f)
    shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
ext = os.path.splitext(structure)[1]
shutil.copy(structure, os.path.join(work, "caps_input" + ext))
with open(os.path.join(work, "dl_f_path"), "w") as f:
    f.write("library = lib/\nsolvent = solvent/\noutput = out/\ncontrol = caps.control\n")
template = open(os.path.join(src, "control_files", "ethanol.control")).read().splitlines()
lines = template[:]
lines[0] = "CAPS validation run"
# --dlpoly: families DL_FIELD cannot write for LAMMPS (GROMOS); the DL_POLY FIELD file is then the reference
dlpoly_only = "--dlpoly" in sys.argv
lines[2] = ("none" if dlpoly_only else "lammps") + "     * Seconday output files (gromacs, chemshell or none)."
lines[4] = f"{ff} * Type of force field require (see list below for choices)."
lines[9] = f"{udff}     * Include user-defined information. Put 'none' or a .udff filename"
lines[11] = f"caps_input{ext}  * Configuration file."
lines[14] = "0        * Optimise FIELD output size, if possible? 1=yes  0=no"
# DL_FIELD's own single-point energy by component ("MM calculation"), written to dl_field.output
for k, l in enumerate(lines):
    if "MM calculation" in l:
        lines[k] = "1          * MM calculation 1=on  0=off"
        break
# --keys: write DL_FIELD's library keys in the FIELD file ("Atom display 1") instead of the standard type names
if "--keys" in sys.argv:
    for k, l in enumerate(lines):
        if "Atom display" in l:
            lines[k] = "1        * Atom display: 1 = DL_FIELD format. 2 = Standard format"
            break
open(os.path.join(work, "caps.control"), "w").write("\n".join(lines) + "\n")
r = subprocess.run(["./dl_field"], cwd=work, capture_output=True, text=True)
os.makedirs(outdir, exist_ok=True)
open(os.path.join(outdir, "dl_field.stdout"), "w").write(r.stdout + r.stderr)
found = False
for root, dirs, files in os.walk(os.path.join(work, "out")):
    for f in files:
        if f in ("lammps.in", "lammps1.data", "dl_field.output", "dlf_notation.output", "dl_poly.FIELD", "dl_poly.CONFIG"):
            shutil.copy(os.path.join(root, f), os.path.join(outdir, f))
            found = found or f == ("dl_poly.FIELD" if dlpoly_only else "lammps1.data")
errors = [l for l in r.stdout.splitlines() if l.strip().startswith("Error")]
field = os.path.join(outdir, "dl_poly.FIELD")
if dlpoly_only and found and "close" not in open(field).read():
    errors.append("DL_FIELD stopped before finishing the FIELD file: " + " | ".join(r.stdout.strip().splitlines()[-3:]))
if not found or errors:
    print("\n".join(errors) if errors else r.stdout[-3000:])
    sys.exit(1)
print("DL_FIELD", ff, "->", outdir)
