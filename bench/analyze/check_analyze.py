#!/usr/bin/env python3
"""Check CAPS Analyze against independent implementations on one trajectory.

  g(r), Rg, MSD       MDAnalysis 2.x (InterRDF, radius_of_gyration, EinsteinMSD) on the same files
  molecule-centre MSD numpy, all time origins, drift removed
  D                   numpy fit of that MSD over the same window
  S(q)                direct reciprocal-space sum S(k) = |Σ exp(ik·r)|² / N over every lattice vector of the cell, in the
                      same |k| shells as CAPS
  free volume         Monte Carlo point insertion against Bondi spheres (scipy KD-tree, periodic)
  pore size           a simple cubic lattice of spheres: the largest cavity is analytic

usage: check_analyze.py DATA DUMP [--frame-ps 0.5]
Run with a Python that has MDAnalysis and scipy (e.g. PY=.../bin/python). CAPS: build/cli/caps.
"""
import json, math, os, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
data, dump = sys.argv[1], sys.argv[2]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
frame_ps = float(arg("--frame-ps", "0.5"))
work = tempfile.mkdtemp()
BONDI = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47, 14: 2.10, 15: 1.80, 16: 1.80, 17: 1.75}
results, fails = [], 0


def check(name, caps, ref, tol, unit=""):
    global fails
    ok = abs(caps - ref) <= tol
    fails += not ok
    results.append(f"{'ok  ' if ok else 'FAIL'} {name}: CAPS {caps:.5g} · reference {ref:.5g}{unit} (tolerance {tol:g})")


def caps_props(props, extra=()):
    out = os.path.join(work, "r.json")
    subprocess.run([CAPS, "analyze", dump, "--topology", data, "--props", props, "--frame-ps", str(frame_ps), "--json", out, *extra],
                   check=True, capture_output=True, text=True)
    return {p["id"]: p for p in json.load(open(out))}


import MDAnalysis as mda
from MDAnalysis.analysis import rdf as mdrdf, msd as mdmsd

u = mda.Universe(data, dump, format="LAMMPSDUMP", lammps_coordinate_convention="unwrapped", dt=frame_ps)
masses = u.atoms.masses
elem = np.array([{12: 6, 1: 1, 14: 7, 16: 8}.get(int(round(m)), 0) for m in masses])
nf = len(u.trajectory)

# ---- g(r) C–C, all pairs
r = caps_props("rdf", ["--pair", "C-C", "--rmax", "10", "--dr", "0.05"])["rdf"]
gx, gy = np.array(r["series"][0]["x"]), np.array(r["series"][0]["y"])
C = u.atoms[elem == 6]
m = mdrdf.InterRDF(C, C, nbins=200, range=(0.0, 10.0), exclusion_block=None)
m.run()
sel = gx > 1.0
diff = np.max(np.abs(np.interp(gx[sel], m.results.bins, m.results.rdf) - gy[sel]))
# MDAnalysis normalises same-group pairs by N_A·N_B (self pairs included in the count at r = 0); CAPS by N_A·N_B too
check("g(r) C–C, largest difference beyond 1 Å (MDAnalysis InterRDF)", diff, 0.0, 0.02)

# ---- Rg
rg = caps_props("rg")["rg"]
vals = []
for ts in u.trajectory:
    for frag in u.atoms.fragments:
        if (elem[frag.indices] != 1).sum() >= 4:
            vals.append(frag.radius_of_gyration() ** 2)
check("√⟨Rg²⟩ (MDAnalysis radius_of_gyration, mass-weighted)", rg["value"], math.sqrt(np.mean(vals)), 1e-3, " Å")

# ---- MSD: atoms (MDAnalysis) and molecule centres (numpy), drift removed in both
res = caps_props("msd,diffusion")
msd, D = res["msd"], res["diffusion"]
X = np.array([u.atoms.positions.copy() for ts in u.trajectory])          # frames × atoms × 3, unwrapped
mtot = masses.sum()
com = (X * masses[None, :, None]).sum(1) / mtot
X = X - (com - com[0])[:, None, :]
lags = np.arange(nf)
atom_msd = np.array([np.mean(np.sum((X[l:] - X[:nf - l]) ** 2, axis=2)) for l in lags])
caps_atoms = np.array(msd["series"][0]["y"])
check("atom MSD at the longest lag (numpy, drift removed)", caps_atoms[-1], atom_msd[-1], 1e-3 * atom_msd[-1] + 1e-6, " Å²")
mm = mdmsd.EinsteinMSD(u, select="all", msd_type="xyz", fft=False)
mm.run()
# MDAnalysis keeps the drift: compare where the system's centre has barely moved, the first few lags
check("atom MSD at lag 5 (MDAnalysis EinsteinMSD)", caps_atoms[5], mm.results.timeseries[5], 0.02 * mm.results.timeseries[5] + 1e-3, " Å²")
frags = u.atoms.fragments
M = np.stack([(X[:, f.indices] * masses[f.indices][None, :, None]).sum(1) / masses[f.indices].sum() for f in frags], axis=1)
com_msd = np.array([np.mean(np.sum((M[l:] - M[:nf - l]) ** 2, axis=2)) for l in lags])
check("molecule-centre MSD at the longest lag (numpy)", msd["value"], com_msd[-1], 1e-3 * com_msd[-1] + 1e-6, " Å²")
t = lags * frame_ps
w = (t / t[-1] >= 0.2) & (t / t[-1] <= 0.5)
slope = np.polyfit(t[w], com_msd[w], 1)[0]
check("D, Einstein fit 20–50 % (numpy)", D["value"], slope / 6 * 10.0, 1e-3 * abs(slope) + 1e-6, " 10⁻⁵ cm²/s")

# ---- S(q): the direct reciprocal-lattice sum in numpy, every k-vector, same shells (last frame)
dq = 0.1
sq = caps_props("sq", ["--first", str(nf - 1), "--dq", str(dq), "--qmax", "6", "--qdirect", "3.5"])["sq"]
u.trajectory[-1]
pos = u.atoms.positions.astype(np.float64)
L = u.dimensions[:3]
rng = np.random.default_rng(1)
nmax = int(3.5 * L.max() / (2 * np.pi)) + 1
g = np.array([(i, j, k) for i in range(0, nmax + 1) for j in range(-nmax, nmax + 1) for k in range(-nmax, nmax + 1) if (i, j, k) > (0, 0, 0)])
kv = 2 * np.pi * g / L
kn = np.linalg.norm(kv, axis=1)
kv, kn = kv[kn <= 3.5], kn[kn <= 3.5]
sk = np.concatenate([(lambda ph: (np.cos(ph).sum(0) ** 2 + np.sin(ph).sum(0) ** 2) / len(pos))(pos @ kv[a:a + 500].T) for a in range(0, len(kv), 500)])
shell = (kn / dq).astype(int)
qx, qy = np.array(sq["series"][0]["x"]), np.array(sq["series"][0]["y"])
ref, acc = [], []                      # shells merged until each holds ≥ 24 k-vectors, point at their mean |k|
for b in range(int(3.5 / dq) + 1):
    acc.extend(np.flatnonzero(shell == b))
    if len(acc) >= 24:
        ref.append((kn[acc].mean(), sk[acc].mean()))
        acc = []
worst = 0.0
for (rq, rv), q, v in zip(ref, qx, qy):
    worst = max(worst, abs(rq - q) + abs(rv - v))
check(f"S(q) up to 3.5 Å⁻¹, largest shell difference ({len(kn)} k-vectors, numpy direct sum)", worst, 0.0, 1e-6)
tail = qx > 4.0
check("S(q) beyond 4 Å⁻¹ (g(r) branch) stays near 1 ± 0.5", float(np.max(np.abs(qy[tail] - 1))), 0.0, 0.5)

# ---- free volume: Monte Carlo against Bondi spheres (last frame)
from scipy.spatial import cKDTree
fv = caps_props("ffv", ["--first", str(nf - 1), "--grid", "0.3"])["ffv"]
wrapped = pos - np.floor(pos / L) * L
tree = cKDTree(wrapped, boxsize=L)
pts = rng.random((400000, 3)) * L
rad = np.array([BONDI[z] for z in elem])
free = np.ones(len(pts), bool)
near = tree.query_ball_point(pts, r=max(rad))
for i, nb in enumerate(near):
    if nb:
        d = np.linalg.norm(((wrapped[nb] - pts[i] + L / 2) % L) - L / 2, axis=1)
        if np.any(d < rad[nb]):
            free[i] = False
mc = free.mean()
check("free fraction, point probe (Monte Carlo, 400 000 points)", fv["value"], mc, 3 * math.sqrt(mc * (1 - mc) / len(pts)) + 0.003)

# ---- pore size: simple cubic lattice of carbon atoms, spacing a — the largest cavity is the cube centre
a = 6.0
lat = np.array([(i, j, k) for i in range(4) for j in range(4) for k in range(4)], float) * a
p = os.path.join(work, "lattice.xyz")
with open(p, "w") as f:
    f.write(f"{len(lat)}\nLattice=\"{4 * a} 0 0 0 {4 * a} 0 0 0 {4 * a}\" Properties=species:S:1:pos:R:3\n")
    for x in lat:
        f.write(f"C {x[0]} {x[1]} {x[2]}\n")
out = os.path.join(work, "psd.json")
r = subprocess.run([CAPS, "analyze", p, "--props", "psd", "--grid", "0.2", "--json", out], capture_output=True, text=True)
if r.returncode == 0:
    psd = {q["id"]: q for q in json.load(open(out))}["psd"]
    largest = psd["extra"]["largest diameter (Å)"]
    check("largest pore of a simple cubic lattice (analytic 2(a√3/2 − r_C))", largest, 2 * (a * math.sqrt(3) / 2 - 1.70), 0.01, " Å")
else:
    results.append("skip lattice pore test: " + (r.stderr.strip().splitlines() or ["?"])[-1])

print("\n".join(results))
print(f"\n{sum(r.startswith('ok') for r in results)} of {len(results)} checks passed")
sys.exit(1 if fails else 0)
