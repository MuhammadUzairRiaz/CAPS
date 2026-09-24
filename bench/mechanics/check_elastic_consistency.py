#!/usr/bin/env python3
"""Two independent routes to the elastic constants of one amorphous polymer cell must meet at low temperature.

  static    caps elastic --method strain: minimise, strain ±1e-4, re-minimise (non-affine relaxation by minimisation)
  fluct     caps elastic --method fluct-run: NVT dynamics at 5 K from the same minimum, the virial stress sampled at
            every step and the Born term every 100 steps; Born minus stress fluctuations plus the kinetic term

In the harmonic limit (T → 0) the two are the same quantity. The fluctuation term cancels ~99 % of the bonded Born
term, so the fluctuation route needs very long runs (100 ps of 1300 atoms: block errors of ~4 GPa; errors fall as
1/√t): the check reports INCONCLUSIVE until the error is below a quarter of the constants. The conclusive
low-temperature agreement test is Mechanics.FluctuationConstantsOfColdCrystal (unit tests). (A LAMMPS reference for the static route is not practical on this model: LAMMPS's minimisers stop at
residual forces of ~0.2 kcal/mol/Å on the truncated LJ surface, so strained cells relax only partly and come out
stiffer. CAPS freezes the pair set at the minimum, which makes the surface smooth.)

usage: check_elastic_consistency.py [DATA] [--ps 100] [--keep DIR]
"""
import json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
args = [a for a in sys.argv[1:] if not a.startswith("--")]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
data = args[0] if args else os.path.join(ROOT, "out", "eq_full.data")
ps = float(arg("--ps", "100"))
work = arg("--keep", "") or tempfile.mkdtemp()
os.makedirs(work, exist_ok=True)
W = lambda f: os.path.join(work, f)


def run(*a):
    r = subprocess.run([CAPS, *map(str, a)], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr.strip())
    return r.stdout


def props(path):
    return {p["id"]: p for p in json.load(open(path))}


run("relax", data, "-o", W("min.data"), "--no-pushoff", "--quiet")
run("elastic", W("min.data"), "--method", "strain", "--strain", "1e-4", "--json", W("static.json"), "--quiet")
run("elastic", W("min.data"), "--method", "fluct-run", "--temp", 5, "--ps", ps, "--equilibrate", 5, "--new-velocities", "--seed", 11,
    "--json", W("fluct.json"), "--quiet")
s, f = props(W("static.json")), props(W("fluct.json"))
rows = [("C11", "cij", "C11 (GPa)"), ("C22", "cij", "C22 (GPa)"), ("C33", "cij", "C33 (GPa)"), ("C12", "cij", "C12 (GPa)"),
        ("C44", "cij", "C44 (GPa)"), ("C55", "cij", "C55 (GPa)"), ("C66", "cij", "C66 (GPa)")]
print(f"{'':6s}{'static':>10s}{'fluct 5 K':>12s}")
fails = 0
for name, pid, key in rows:
    a, b = s[pid]["extra"][key], f[pid + "_fluct"]["extra"][key]
    print(f"{name:6s}{a:10.3f}{b:12.3f}")
ferr = f["cij_fluct"]["extra"]
err11 = (ferr.get("C11 error (GPa)", 0) + ferr.get("C22 error (GPa)", 0) + ferr.get("C33 error (GPa)", 0)) / 3
print(f"block error of the fluctuation C11..C33: {err11:.3f} GPa")
conclusive = err11 < 0.25 * s["cij"]["value"]
for label, pid in [("bulk K", "bulk"), ("shear G", "shear"), ("Young E", "youngs")]:
    a, b = s[pid]["value"], f[pid + "_fluct"]["value"]
    ok = abs(a - b) <= 0.15 * abs(a) + 3 * err11
    fails += conclusive and not ok
    tag = ("ok  " if ok else "FAIL") if conclusive else "n/a "
    print(f"{tag} {label}: static {a:.3f} GPa, fluctuations at 5 K {b:.3f} GPa ({100 * (b / a - 1):+.1f} %)")
if not conclusive:
    print(f"INCONCLUSIVE: the fluctuation route's error ({err11:.2f} GPa) is too large to test the static constants; run longer (--ps)")
for n in f["cij_fluct"]["notes"]:
    print("   note:", n)
print(f"(work: {work})")
sys.exit(1 if fails else 0)
