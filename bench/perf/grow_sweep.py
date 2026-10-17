#!/usr/bin/env python3
"""Robustness sweep of the chain grower: every polymer and copolymer of data/polymers/library.json, every growth method,
a grid of chain counts × DP and two target densities, each cell checked independently of CAPS.

    CAPS_LIB=build/capi/libcaps.dylib python3 bench/perf/grow_sweep.py [--jobs 8] [--timeout 60] [--polymer P0] \
        [--methods trials,rosenbluth,rosenbluth_lj] [--densities 0.3,1.1] [--json OUT.json]

Grid (chains × DP), scaled down per polymer so that no cell exceeds 15,000 atoms:
    1 × 1, 1 × 200, 50 × 5, 200 × 2, 20 × 40, and a mid cell of 10 chains near 14,000 atoms.
The Rosenbluth methods run on 1 × 1, 50 × 5 and 20 × 40 (and the mid cell with --full).
Cells grow as the recipe and Python's caps.polymer grow them by default: the contact scale on auto (lowered only where
a chain cannot be placed), seed 1, one thread per cell (the cells run in parallel).

Each cell runs in its own process with a time limit, then is checked:
  · it was grown (no error, no time-out);
  · molecules: as many as the chains asked for; residues per molecule: the DP asked for;
  · atoms: per molecule the sum of its units' atoms plus two end-cap hydrogens;
  · density: the cell's mass over its volume within 1 % of the target;
  · bonds: every bond 0.8 to 2.4 Å long in the unwrapped coordinates (none drawn across the cell);
  · contacts: no two atoms more than three bonds apart closer than 0.6 Å, the least the eased contact limits allow
    (minimum image); the closest contact of each cell is in the --json output.
Prints one row per polymer (a mark per cell: . pass, F check failed, E error, T time-out) and a summary."""
import argparse
import ctypes as C
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ["CAPS_SRC"]) if os.environ.get("CAPS_SRC") else HERE.parents[1]
GRID = [(1, 1), (1, 200), (50, 5), (200, 2), (20, 40), ("mid", None)]
GRID_ROSENBLUTH = [(1, 1), (50, 5), (20, 40)]
MAX_ATOMS = 15000
# closest non-bonded contact allowed: a chain stuck in a crowd may ease its limits down to 0.4 of the full ones (H–H
# 1.95 Å → 0.78 Å, pairs four bonds apart at 85 % of that, less the 0.05 Å accepted)
CONTACT = 0.6


def lib_path():
    return os.environ.get("CAPS_LIB") or str(ROOT / "build" / "capi" / ("libcaps.dylib" if sys.platform == "darwin" else "libcaps.so"))


def library_entries():
    lib = json.load(open(ROOT / "data" / "polymers" / "library.json"))
    by_id = {p["id"]: p for p in lib["polymers"]}
    out = [(p["id"], p["name"], {"units": [p["smiles"]], "sequence": "homopolymer"}) for p in lib["polymers"]]
    for c in lib.get("copolymers", []):
        spec = {"units": [by_id[u]["smiles"] if u in by_id else u for u in c["units"]], "sequence": c.get("sequence", "random")}
        for k in ("weights", "blocks", "pattern"):
            if k in c:
                spec[k] = c[k]
        out.append((c["id"], c["name"], spec))
    return out


# ---------------------------------------------------------------- one cell (child process)

class GrowOpts(C.Structure):
    _fields_ = [("chains", C.c_int32), ("dp", C.c_int32), ("tacticity", C.c_int32), ("seed", C.c_uint64), ("box", C.c_double),
                ("density", C.c_double), ("contact_scale", C.c_double), ("curve", C.c_int32)]


class Summary(C.Structure):
    _fields_ = [("atoms", C.c_int64), ("bonds", C.c_int64), ("molecules", C.c_int64), ("frames", C.c_int64),
                ("bonds_from_file", C.c_int32), ("has_charges", C.c_int32), ("cell_valid", C.c_int32), ("unwrapped", C.c_int32),
                ("cell_a", C.c_double), ("cell_b", C.c_double), ("cell_c", C.c_double), ("volume", C.c_double), ("density", C.c_double),
                ("total_mass", C.c_double), ("total_charge", C.c_double), ("format", C.c_char * 24)]


def load():
    L = C.CDLL(lib_path())
    L.caps_grow_chains.restype = C.c_void_p
    L.caps_grow_chains.argtypes = [C.c_char_p, C.POINTER(GrowOpts), C.c_void_p, C.c_void_p, C.c_char_p, C.c_int32]
    L.caps_close.argtypes = [C.c_void_p]
    L.caps_save.argtypes = [C.c_void_p, C.c_char_p]
    L.caps_summary_get.argtypes = [C.c_void_p, C.POINTER(Summary)]
    L.caps_last_error.restype = C.c_char_p
    L.caps_unit_info.argtypes = [C.c_char_p, C.c_char_p, C.c_int32]
    return L


def unit_info(L, smiles):
    buf = C.create_string_buffer(4096)
    L.caps_unit_info(smiles.encode(), buf, len(buf))
    j = json.loads(buf.value.decode())
    if not j.get("ok"):
        raise ValueError(j.get("error", "bad unit"))
    return j


def unit_atoms(L, smiles):
    return int(unit_info(L, smiles)["atoms"])


def read_data(path):
    """atoms {id: (mol, residue, (x, y, z) unwrapped)}, bonds, box edges of a LAMMPS data file written by CAPS."""
    box = [0.0, 0.0, 0.0]
    atoms, bonds, sec = {}, [], None
    for ln in Path(path).read_text().splitlines():
        body, _, comment = ln.partition("#")
        t = body.split()
        if not t:
            continue
        if len(t) >= 4 and t[2].endswith("lo") and t[3].endswith("hi"):
            box["xyz".index(t[2][0])] = float(t[1]) - float(t[0])
            continue
        if t[0] in ("Atoms", "Bonds", "Masses", "Velocities", "Angles", "Dihedrals", "Impropers") or t[0].endswith("Coeffs"):
            sec = t[0]
            continue
        if sec == "Atoms" and len(t) >= 7:
            x = [float(t[4]), float(t[5]), float(t[6])]
            if len(t) >= 10:
                x = [x[k] + int(t[7 + k]) * box[k] for k in range(3)]
            c = comment.split()
            res = (c[1], c[2] if len(c) > 2 else "") if len(c) >= 2 and c[0] == "res" else None
            atoms[int(t[0])] = (int(t[1]), res, x)
        elif sec == "Bonds" and len(t) >= 4:
            bonds.append((int(t[2]), int(t[3])))
    return atoms, bonds, box


def closest_contact(X, adj, box):
    n = len(X)
    rc = 3.0
    nc = [max(1, int(L // rc)) for L in box]
    cs = [box[k] / nc[k] for k in range(3)]
    bins, key = {}, []
    for i, p in enumerate(X):
        k = tuple(int((p[d] % box[d]) // cs[d]) % nc[d] for d in range(3))
        key.append(k)
        bins.setdefault(k, []).append(i)
    best = 1e9
    for i in range(n):
        near = {i}
        fr = [i]
        for _ in range(3):
            nx = []
            for u in fr:
                for w in adj[u]:
                    if w not in near:
                        near.add(w), nx.append(w)
            fr = nx
        kx, ky, kz = key[i]
        seen = set()
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    k = ((kx + dx) % nc[0], (ky + dy) % nc[1], (kz + dz) % nc[2])
                    if k in seen:
                        continue
                    seen.add(k)
                    for j in bins.get(k, ()):
                        if j <= i or j in near:
                            continue
                        d2 = 0.0
                        for d in range(3):
                            v = X[j][d] - X[i][d]
                            v -= box[d] * round(v / box[d])
                            d2 += v * v
                        if d2 < best * best:
                            best = math.sqrt(d2)
    # an atom against its own periodic images (a tiny cell): the shortest cell edge
    best = min(best, min(box))
    return best


def run_one(case):
    """Grows and checks one cell; returns a dict (printed as JSON by the child process)."""
    L = load()
    spec_in = case["spec"]
    names = [f"U{chr(65 + k)}" for k in range(len(spec_in["units"]))]
    spec = {"units": [{"name": n, "smiles": s} for n, s in zip(names, spec_in["units"])], "dp": case["dp"],
            "sequence": spec_in.get("sequence", "homopolymer"), "threads": 1}
    for k in ("weights", "blocks", "pattern"):
        if k in spec_in:
            spec[k] = spec_in[k]
    if case["method"] != "trials":
        spec["method"] = case["method"]
    n_unit = {n: unit_atoms(L, s) for n, s in zip(names, spec_in["units"])}
    g = GrowOpts(case["chains"], case["dp"], 0, case.get("seed", 1), 0.0, case["density"], -1.0, 1)
    rep = C.create_string_buffer(16384)
    t0 = time.perf_counter()
    h = L.caps_grow_chains(json.dumps(spec).encode(), C.byref(g), None, None, rep, len(rep))
    dt = time.perf_counter() - t0
    out = {"time": dt}
    if not h:
        out["status"] = "E"
        out["why"] = L.caps_last_error().decode()
        return out
    s = Summary()
    L.caps_summary_get(h, C.byref(s))
    notes = rep.value.decode()
    out["atoms"], out["density_got"] = s.atoms, s.density
    out["scale"] = next((ln.split("contact scale ")[1].split()[0] for ln in notes.splitlines() if "grown at contact scale" in ln), "1")
    with tempfile.TemporaryDirectory() as tmp:
        f = os.path.join(tmp, "c.data")
        L.caps_save(h, f.encode())
        L.caps_close(h)
        atoms, bonds, box = read_data(f)
    fails = []
    ids = sorted(atoms)
    idx = {a: k for k, a in enumerate(ids)}
    X = [atoms[a][2] for a in ids]
    adj = [[] for _ in ids]
    for a, b in bonds:
        adj[idx[a]].append(idx[b]), adj[idx[b]].append(idx[a])
    mols = {}
    for a in ids:
        mols.setdefault(atoms[a][0], []).append(a)
    if len(mols) != case["chains"]:
        fails.append(f"{len(mols)} molecules, {case['chains']} asked")
    bad_dp = bad_atoms = 0
    for m, members in mols.items():
        res = {}
        for a in members:
            r = atoms[a][1]
            if r:
                res.setdefault(r[0], r[1])
        if len(res) != case["dp"]:
            bad_dp += 1
        want = sum(n_unit.get(nm, 0) for nm in res.values()) + 2
        if len(members) != want:
            bad_atoms += 1
    if bad_dp:
        fails.append(f"{bad_dp} molecules without {case['dp']} residues")
    if bad_atoms:
        fails.append(f"{bad_atoms} molecules with the wrong atom count")
    if abs(s.density - case["density"]) > 0.01 * case["density"]:
        fails.append(f"density {s.density:.4f}")
    long_bonds = 0
    for a, b in bonds:
        d = math.dist(atoms[a][2], atoms[b][2])
        if d < 0.8 or d > 2.4:
            long_bonds += 1
    if long_bonds:
        fails.append(f"{long_bonds} bonds outside 0.8-2.4 Å")
    closest = closest_contact(X, adj, box)
    out["closest"] = closest
    if closest < CONTACT:
        fails.append(f"closest contact {closest:.2f} Å")
    out["status"] = "F" if fails else "."
    if fails:
        out["why"] = "; ".join(fails)
    return out


# ---------------------------------------------------------------- the sweep (parent)

def cases_for(L, pid, spec, methods, densities, full):
    n = max(unit_atoms(L, u) for u in spec["units"])
    out = []
    for method in methods:
        grid = GRID if method == "trials" or full else GRID_ROSENBLUTH
        for chains, dp in grid:
            if chains == "mid":
                chains, dp = 10, max(1, (14000 // 10 - 2) // n)
            # at most 15,000 atoms: fewer units per chain, then fewer chains
            while chains * (dp * n + 2) > MAX_ATOMS and dp > 1:
                dp = max(1, dp * MAX_ATOMS // (chains * (dp * n + 2)))
            while chains * (dp * n + 2) > MAX_ATOMS and chains > 1:
                chains -= 1
            for dens in densities:
                out.append({"id": pid, "spec": spec, "method": method, "chains": chains, "dp": dp, "density": dens})
    return out


def grid_label(r):
    """The grid point a cell stands for (its chains × DP before the 15,000-atom scaling)."""
    c, d = r["chains"], r["dp"]
    if (c, d) in ((1, 1), (50, 5), (200, 2), (20, 40)):
        return f"{c}x{d}"
    if c == 1:
        return "1xlong"
    if c == 10:
        return "mid"
    return "scaled"


def summary(results):
    """Pass / check failed / error / time-out per method, grid point and density, and the grow times."""
    table = {}
    for r in results:
        k = (r["method"], grid_label(r), round(r["density"], 2))
        t = table.setdefault(k, {".": 0, "F": 0, "E": 0, "T": 0, "times": []})
        t[r["status"]] += 1
        if r["status"] in ".F" and "time" in r:
            t["times"].append(r["time"])
    print(f"{'method':14s} {'grid':8s} {'g/cm3':>5s} {'pass':>5s} {'F':>4s} {'E':>4s} {'T':>4s} {'median s':>9s} {'max s':>7s}")
    for k in sorted(table, key=lambda k: (k[0] != "trials", k[0], k[1], k[2])):
        t = table[k]
        ts = sorted(t["times"])
        med = f"{ts[len(ts) // 2]:9.2f}" if ts else f"{'-':>9s}"
        mx = f"{ts[-1]:7.2f}" if ts else f"{'-':>7s}"
        print(f"{k[0]:14s} {k[1]:8s} {k[2]:5.1f} {t['.']:5d} {t['F']:4d} {t['E']:4d} {t['T']:4d} {med} {mx}")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=60)
    ap.add_argument("--polymer", default="", help="ids or names containing this text")
    ap.add_argument("--methods", default="trials,rosenbluth,rosenbluth_lj")
    ap.add_argument("--densities", default="0.3,1.1")
    ap.add_argument("--full", action="store_true", help="the whole grid for the Rosenbluth methods too")
    ap.add_argument("--json", default="")
    ap.add_argument("--case", default="", help=argparse.SUPPRESS)
    ap.add_argument("--summarise", default="", help="print the summary of an earlier --json result and stop")
    a = ap.parse_args()
    if a.summarise:
        results = json.load(open(a.summarise))
        for r in results:   # results of the first version kept the measured density under "density"
            if "density_got" not in r and r.get("status") in ".F":
                r["density"] = min((0.3, 1.1), key=lambda x: abs(x - r["density"]))
        summary(results)
        return
    if a.case:
        print(json.dumps(run_one(json.loads(a.case))))
        return
    L = load()
    entries = [e for e in library_entries() if not a.polymer or a.polymer in e[0] or a.polymer.lower() in e[1].lower()]
    methods = a.methods.split(",")
    dens = [float(x) for x in a.densities.split(",")]
    cases = []
    for pid, name, spec in entries:
        try:
            cases += cases_for(L, pid, spec, methods, dens, a.full)
        except ValueError as e:
            print(f"{pid} {name}: {e}", file=sys.stderr)

    def go(case):
        t0 = time.perf_counter()
        try:
            p = subprocess.run([sys.executable, __file__, "--case", json.dumps(case)], capture_output=True, text=True, timeout=a.timeout,
                               env={**os.environ, "CAPS_LIB": lib_path()})
            r = json.loads(p.stdout.strip().splitlines()[-1]) if p.stdout.strip() else {"status": "E", "why": p.stderr.strip()[-300:]}
        except subprocess.TimeoutExpired:
            r = {"status": "T", "why": f"over {a.timeout:.0f} s"}
        r["wall"] = time.perf_counter() - t0
        return {**case, **r}

    t0 = time.perf_counter()
    done = [0]

    def go_counted(case):
        r = go(case)
        done[0] += 1
        if done[0] % 100 == 0 or done[0] == len(cases):
            print(f"  {done[0]}/{len(cases)} cells ({time.perf_counter() - t0:.0f} s)", file=sys.stderr, flush=True)
        return r

    with ThreadPoolExecutor(a.jobs) as ex:
        results = list(ex.map(go_counted, cases))
    wall = time.perf_counter() - t0
    # one row per polymer: a mark per cell in case order
    print(f"{'polymer':8s} {'cells':>5s} {'pass':>4s}  marks (per method: grid × density)")
    by = {}
    for r in results:
        by.setdefault(r["id"], []).append(r)
    total = {".": 0, "F": 0, "E": 0, "T": 0}
    times = []
    for pid, rs in by.items():
        marks = "".join(r["status"] for r in rs)
        for r in rs:
            total[r["status"]] += 1
            if r["status"] in ".F":
                times.append(r["time"])
        print(f"{pid:8s} {len(rs):5d} {sum(r['status'] == '.' for r in rs):4d}  {marks}")
    print()
    summary(results)
    print(f"cells {len(results)}: pass {total['.']}, check failed {total['F']}, error {total['E']}, time-out {total['T']}  (wall {wall:.0f} s)")
    if times:
        times.sort()
        print(f"grow time of finished cells: median {times[len(times) // 2]:.2f} s, 95th percentile {times[int(0.95 * (len(times) - 1))]:.2f} s, max {times[-1]:.2f} s")
    print()
    print("failures:")
    for r in results:
        if r["status"] != ".":
            print(f"  {r['id']} {r['method']} {r['chains']}x{r['dp']} @{r['density']}: {r['status']} {r.get('why', '')[:200]}")
    if a.json:
        json.dump(results, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
