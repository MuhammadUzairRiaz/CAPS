#!/usr/bin/env python3
"""Table 2 of the CAPS specification: CAPS Pack against packmol on identical packmol inputs.

Every output is checked independently with `caps contacts` (smallest distance between atoms of different molecules,
periodic minimum image). For one-species files the molecules are taken as consecutive atom blocks, because PDB files
over 99 999 atoms cannot number their CONECT records reliably.

usage: run_bench.py [case ...]   cases: water1000 water5000 water20000 water50000 ps400 slab
"""
import os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
CAPS = os.path.join(HERE, "..", "..", "build", "cli", "caps")
CASES = {  # name: (description, molecules, atoms, density, molecule size for the check or None)
    "water1000": ("water, periodic", 1000, 3000, "1.00", 3),
    "water5000": ("water, periodic", 5000, 15000, "1.00", 3),
    "water20000": ("water, periodic", 20000, 60000, "1.00", 3),
    "water50000": ("water, periodic", 50000, 150000, "1.00", 3),
    "ps400": ("polystyrene DP 8 chains, periodic", 400, 52000, "0.50", 130),
    "slab": ("water slab + PS chains (plane regions), periodic", 4060, 19800, "1.00 / 0.50", None),
}
names = sys.argv[1:] or list(CASES)
os.chdir(HERE)

LIMIT = 3600  # s per run; a run that does not finish in time is reported as such

def run(cmd, stdin=None):
    t = time.perf_counter()
    try:
        p = subprocess.run(cmd, stdin=open(stdin) if stdin else None, capture_output=True, text=True, timeout=LIMIT)
    except subprocess.TimeoutExpired:
        return time.perf_counter() - t, None
    return time.perf_counter() - t, p

def contacts(pdb, size):
    cmd = [CAPS, "contacts", pdb, "--tol", "2.0"] + (["--molecule-size", str(size)] if size else [])
    _, p = run(cmd)
    m = re.search(r"distance ([0-9.]+) Å · (\d+) pairs", p.stdout)
    return (float(m.group(1)), int(m.group(2))) if m else (float("nan"), -1)

rows = []
for name in names:
    desc, nmol, natoms, rho, size = CASES[name]
    for tool, threads, cmd, stdin in (("CAPS", 10, [CAPS, "pack", f"{name}_caps.inp", "--quiet"], None),
                                      ("CAPS", 1, [CAPS, "pack", f"{name}_caps.inp", "--quiet", "--threads", "1"], None),
                                      ("packmol 21.2.1", 1, ["packmol"], f"{name}_packmol.inp")):
        wall, p = run(cmd, stdin)
        out = f"{name}_caps.pdb" if tool == "CAPS" else f"{name}_packmol.pdb"
        finished = p is not None and p.returncode == 0 and (tool == "CAPS" or "Success!" in p.stdout)
        d, close = contacts(out, size) if finished else (float("nan"), -1)
        rows.append((name, tool, threads, wall, d, close, finished))
        print(f"{name:11s} {tool:15s} {threads:2d} thr {wall:8.2f} s  dmin {d:.3f} Å  close {close}  {'ok' if finished else 'FAILED'}", flush=True)

with open("results.md", "w") as f:
    f.write("| System | Molecules | Atoms | Target density (g/cm³) | Tool | Threads | Wall time (s) | Minimum distance (Å) | Success | Speed-up |\n")
    f.write("|---|---|---|---|---|---|---|---|---|---|\n")
    for name in names:
        desc, nmol, natoms, rho, _ = CASES[name]
        ref = next(r[3] for r in rows if r[0] == name and r[1] != "CAPS")
        for (_, tool, thr, wall, d, close, ok) in [r for r in rows if r[0] == name]:
            wt = f"{wall:.2f}" if ok else f"> {LIMIT} (stopped)" if wall >= LIMIT else f"{wall:.2f} (failed)"
            dm = f"{d:.3f}" if ok else "–"
            sp = f"{ref / wall:.1f}×" if ok else "–"
            f.write(f"| {desc} | {nmol:,} | {natoms:,} | {rho} | {tool} | {thr} | {wt} | {dm} | "
                    f"{'yes' if ok and close == 0 else 'no'} | {sp} |\n")
print(open("results.md").read())
