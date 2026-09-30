"""Every CAPS water model against LAMMPS: 64 waters packed near 1 g/cm³, the model applied (M sites for the four-site
models), typed by the water group, CAPS's single-point energy (PME) against LAMMPS's run 0 with the exported files
(lj/cut/tip4p/long and pppm/tip4p for the four-site models). Needs lmp in ~/.local/bin with KSPACE and MOLECULE.

    PYTHONPATH=data/python CAPS_LIB=build/capi/libcaps.dylib python3 bench/ff/check_water_models.py [model …]
"""
import ctypes as C
import os
import re
import subprocess
import sys
import tempfile

import caps

L = caps.library()
L.caps_set_electrostatics.argtypes = [C.c_int32, C.c_double, C.c_double, C.c_int32]
L.caps_set_electrostatics(1, 1e-6, 0.6, 6)   # PME, tight
samples = os.environ.get("CAPS_SAMPLES", "samples")
worst = 0.0
for model in sys.argv[1:] or [m["id"] for m in caps.water_models()]:
    d = caps.pack(molecules=[(os.path.join(samples, "water.pdb"), 64)], box=12.42, tolerance=2.0, seed=4)
    d.edit(op="water_model", model=model)
    d.field.assign_groups([{"name": "water", "molecules": "water", "water": model}])
    caps_total = d.energy()["total"]
    out = tempfile.mkdtemp(prefix="caps_water_")
    d.export_engines(out, stem="w", gromacs=False, run="check", kspace_accuracy=1e-7, cutoff=0)
    subprocess.run([os.path.expanduser("~/.local/bin/lmp"), "-in", "w.in", "-log", "log.lammps"], cwd=out, capture_output=True, text=True)
    rows = re.findall(r"^\s*0\s+(-?[\d.]+)", open(os.path.join(out, "log.lammps")).read(), re.M)
    lmp_total = float(rows[-1])
    worst = max(worst, abs(caps_total - lmp_total))
    print(f"{model:14s} CAPS {caps_total:14.4f}  LAMMPS {lmp_total:14.4f}  Δ {caps_total - lmp_total:+.4f} kcal/mol")
print(f"largest difference {worst:.4f} kcal/mol")
sys.exit(0 if worst < 0.02 else 1)
