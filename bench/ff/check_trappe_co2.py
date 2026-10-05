"""TraPPE CO2 in CAPS (Potoff & Siepmann, AIChE J. 47, 1676 (2001)), checked three ways:

1. the library entry (trappe-ua: CO2C, CO2O) against gmso's TraPPE CO2 file (trappe_5site_co2.xml in the mosdef env):
   σ, ε and charges;
2. a box of 100 CO2 typed by CAPS: CAPS's single-point energy (PME) against LAMMPS's run 0 with the exported files;
3. sorption of CO2 into a TraPPE-UA polyethylene cell (Widom and two GCMC points), the sorbate taking TraPPE's rigid
   geometry (C=O 1.16 Å, linear) — reported, not compared (no amorphous-PE reference on disk).
The model's second virial coefficient and a pure-gas GCMC against the virial expansion are tests/test_sorption.cpp
Sorption.TrappeCo2VirialAndGcmc.

    PYTHONPATH=data/python CAPS_LIB=build/capi/libcaps.dylib python3 bench/ff/check_trappe_co2.py
"""
import ctypes as C
import glob
import json
import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import caps

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ok = True

# 1. the library against the reference file
lib = json.load(open(os.path.join(ROOT, "data", "forcefields", "trappe-ua.json")))
pairs = {p["name"]: p["params"] for p in lib["pairs"]}
charges = {t["name"]: t.get("charge") for t in lib["atom_types"]}
ref = glob.glob("/opt/homebrew/Caskroom/miniconda/base/envs/*/lib/python3*/site-packages/gmso/utils/files/gmso_xmls/test_ffstyles/trappe_5site_co2.xml")
if ref:
    x = ET.parse(ref[0]).getroot()
    for v in x.iter("VirtualSiteType"):
        name = "CO2C" if "Carbon" in v.get("name") else "CO2O"
        p = {q.get("name"): float(q.get("value")) for q in v.iter("Parameter")}
        eps_kcal, sig_A, q = p["epsilon"] / 4.184, p["sigma"] * 10, float(v.get("charge"))
        d_eps = abs(pairs[name][0] - eps_kcal) / eps_kcal
        good = d_eps < 2e-5 and abs(pairs[name][1] - sig_A) < 1e-9 and abs(charges[name] - q) < 1e-12
        ok &= good
        print(f"{name}: CAPS ε {pairs[name][0]:.7f} kcal/mol σ {pairs[name][1]} Å q {charges[name]:+.2f} · reference ε {eps_kcal:.7f} σ {sig_A:.3f} q {q:+.2f}"
              f" · {'same' if good else 'DIFFERENT'} (ε {d_eps:.1e} relative)")
else:
    print("gmso's TraPPE CO2 file not found: step 1 skipped")

# 2. CAPS against LAMMPS
L = caps.library()
L.caps_set_electrostatics.argtypes = [C.c_int32, C.c_double, C.c_double, C.c_int32]
L.caps_set_electrostatics(1, 1e-6, 0.6, 6)   # PME, tight
box = caps.pack(molecules=[("O=C=O", 100)], box=30.0, tolerance=2.5, seed=3)
rep = box.field.assign("trappe-ua")
e = box.energy()
out = tempfile.mkdtemp(prefix="caps_co2_")
box.export_engines(out, stem="co2", gromacs=False, run="check", kspace_accuracy=1e-7, cutoff=10)   # CAPS's own runs: 10 Å
subprocess.run([os.path.expanduser("~/.local/bin/lmp"), "-in", "co2.in", "-log", "log.lammps"], cwd=out, capture_output=True, text=True)
log = open(os.path.join(out, "log.lammps")).read()
rows = re.findall(r"^\s*0\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", log, re.M)
lmp_pe, lmp_vdw, lmp_coul, lmp_long = (float(rows[-1][i]) for i in (0, 5, 6, 7))
d = e["total"] - lmp_pe
ok &= abs(d) < 0.01
print(f"100 CO2 (PME): CAPS vdW {e['vdw']:.4f} Coulomb {e['coulomb']:.4f} total {e['total']:.4f} · LAMMPS vdW {lmp_vdw:.4f}"
      f" Coulomb {lmp_coul + lmp_long:.4f} total {lmp_pe:.4f} · Δ {d:+.5f} kcal/mol")

# 3. CO2 in TraPPE-UA polyethylene
pe = caps.polymer("*CC*", dp=40, chains=6, density=0.85, seed=2)
pe.field.assign("trappe-ua")   # automatic charges: CO2 its own types' (±0.70/0.35), the alkane sites Gasteiger summed (≈ 0)
s = pe.sorption("O=C=O", pressures_kpa=[100, 1000], temperature=300, insertions=40000, steps=40000, seed=1)
print(f"CO2 in TraPPE-UA PE (240 CH2/CH3 sites, 0.85 g/cm³, 300 K): S {s['solubility']:.3f} cm³(STP)/(cm³ atm),"
      f" μex {s['mu_ex']:.2f} kcal/mol · GCMC " + ", ".join(f"{p['pressure_kpa']:g} kPa {p['cm3stp_per_cm3']:.2f} cm³(STP)/cm³" for p in s["isotherm"]))
print("notes:", " · ".join(s["notes"]))
ok &= any("rigid sorbate" in n for n in s["notes"]) and any("electrostatics" in n for n in s["notes"])
print("all TraPPE CO2 checks passed" if ok else "TraPPE CO2 checks FAILED")
sys.exit(0 if ok else 1)
