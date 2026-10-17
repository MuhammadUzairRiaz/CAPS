#!/usr/bin/env python3
"""Times every chain-growth path of CAPS on small, medium and large cells (all at most 15,000 atoms) and checks each
result: atom count, density, and the closest non-bonded contact (pairs more than three bonds apart, minimum image).

    CAPS_LIB=build/capi/libcaps.dylib python3 bench/perf/grow_bench.py [--only NAME] [--timeout 300] [--out FILE.tsv]

The cells go through the C library as CAPS Studio's Grow page does: caps_grow_chains (the repeat-unit grower, every
method, sequence, architecture and region) and caps_grow (the polystyrene grower behind `caps grow`); one case also
takes the live-view callback the Studio uses. Each case runs in its own process (a case over the time limit is reported
as a time-out); wall and CPU seconds of the growth call are given, and how many units were placed with eased limits.
Compare two builds by pointing CAPS_LIB at each."""
import argparse
import ctypes as C
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path


class GrowOpts(C.Structure):
    _fields_ = [("chains", C.c_int32), ("dp", C.c_int32), ("tacticity", C.c_int32), ("seed", C.c_uint64), ("box", C.c_double),
                ("density", C.c_double), ("contact_scale", C.c_double), ("curve", C.c_int32)]


class Summary(C.Structure):
    _fields_ = [("atoms", C.c_int64), ("bonds", C.c_int64), ("molecules", C.c_int64), ("frames", C.c_int64),
                ("bonds_from_file", C.c_int32), ("has_charges", C.c_int32), ("cell_valid", C.c_int32), ("unwrapped", C.c_int32),
                ("cell_a", C.c_double), ("cell_b", C.c_double), ("cell_c", C.c_double), ("volume", C.c_double), ("density", C.c_double),
                ("total_mass", C.c_double), ("total_charge", C.c_double), ("format", C.c_char * 24)]


LIVE = C.CFUNCTYPE(None, C.c_void_p, C.c_char_p, C.c_void_p)


def load():
    here = Path(__file__).resolve().parents[2]
    p = os.environ.get("CAPS_LIB") or str(here / "build" / "capi" / ("libcaps.dylib" if sys.platform == "darwin" else "libcaps.so"))
    L = C.CDLL(p)
    L.caps_grow_chains.restype = C.c_void_p
    L.caps_grow_chains.argtypes = [C.c_char_p, C.POINTER(GrowOpts), C.c_void_p, C.c_void_p, C.c_char_p, C.c_int32]
    L.caps_grow_chains_live.restype = C.c_void_p
    L.caps_grow_chains_live.argtypes = [C.c_char_p, C.POINTER(GrowOpts), C.c_void_p, LIVE, C.c_void_p, C.c_char_p, C.c_int32]
    L.caps_grow.restype = C.c_void_p
    L.caps_grow.argtypes = [C.POINTER(GrowOpts), C.c_void_p, C.c_void_p, C.c_char_p, C.c_int32]
    L.caps_close.argtypes = [C.c_void_p]
    L.caps_save.argtypes = [C.c_void_p, C.c_char_p]
    L.caps_summary_get.argtypes = [C.c_void_p, C.POINTER(Summary)]
    L.caps_last_error.restype = C.c_char_p
    return L


PE = "*CC*"
PS = "*CC(*)c1ccccc1"
PMMA = "*CC(C)(*)C(=O)OC"
CISPI = "[*]C/C=C(C)\\C[*]"
BD = "[*]C/C=C\\C[*]"
ENB = "[*][C@H]1[C@@H]2CC(=CC)[C@@H](C2)[C@H]1[*]"


def unit(s, name="u"):
    return {"name": name, "smiles": s}


# name, kind (chains | ps | live), spec / options, GrowOpts fields. contact_scale −1: automatic (the recipe's, Python's
# and the Studio's default: lowered only where a chain is stuck); 1: fixed full limits.
A = -1.0
CASES = [
    # homopolymers, best-of-k trials (the default method), melts at 0.8 g/cm³ unless noted
    ("PE 1x50 (302 at)", "chains", {"units": [unit(PE)], "dp": 50}, dict(chains=1, density=0.1, contact_scale=A)),
    ("PE 10x100 (6020 at)", "chains", {"units": [unit(PE)], "dp": 100}, dict(chains=10, density=0.8, contact_scale=A)),
    ("PE 24x100 (14448 at)", "chains", {"units": [unit(PE)], "dp": 100}, dict(chains=24, density=0.8, contact_scale=A)),
    ("PS 1x20 (322 at)", "chains", {"units": [unit(PS)], "dp": 20}, dict(chains=1, density=0.1, contact_scale=A)),
    ("PS 10x30 (4820 at)", "chains", {"units": [unit(PS)], "dp": 30}, dict(chains=10, density=0.8, contact_scale=A)),
    ("PS 15x60 (14430 at)", "chains", {"units": [unit(PS)], "dp": 60}, dict(chains=15, density=0.8, contact_scale=A)),
    ("PMMA 1x20 (302 at)", "chains", {"units": [unit(PMMA)], "dp": 20}, dict(chains=1, density=0.1, contact_scale=A)),
    ("PMMA 10x30 (4520 at)", "chains", {"units": [unit(PMMA)], "dp": 30}, dict(chains=10, density=0.8, contact_scale=A)),
    ("PMMA 16x60 (14432 at)", "chains", {"units": [unit(PMMA)], "dp": 60}, dict(chains=16, density=0.8, contact_scale=A)),
    ("cis-PI 1x30 (392 at)", "chains", {"units": [unit(CISPI)], "dp": 30}, dict(chains=1, density=0.1, contact_scale=A)),
    ("cis-PI 10x40 (5220 at)", "chains", {"units": [unit(CISPI)], "dp": 40}, dict(chains=10, density=0.8, contact_scale=A)),
    ("cis-PI 20x55 (14340 at)", "chains", {"units": [unit(CISPI)], "dp": 55}, dict(chains=20, density=0.8, contact_scale=A)),
    # copolymers: random SBR (styrene / butadiene) and a pattern
    ("SBR random 10x40", "chains", {"units": [unit(PS, "S"), unit(BD, "B")], "dp": 40, "sequence": "random", "weights": [0.25, 0.75]},
     dict(chains=10, density=0.8, contact_scale=A)),
    ("SBR random 20x60", "chains", {"units": [unit(PS, "S"), unit(BD, "B")], "dp": 60, "sequence": "random", "weights": [0.25, 0.75]},
     dict(chains=20, density=0.8, contact_scale=A)),
    ("PE/PS pattern AAB 10x30", "chains", {"units": [unit(PE, "E"), unit(PS, "S")], "dp": 30, "sequence": "pattern", "pattern": "AAB"},
     dict(chains=10, density=0.8, contact_scale=A)),
    # tacticity
    ("PS isotactic 10x30", "chains", {"units": [unit(PS)], "dp": 30, "tacticity": "isotactic"}, dict(chains=10, density=0.8, contact_scale=A)),
    ("PS syndiotactic 10x30", "chains", {"units": [unit(PS)], "dp": 30, "tacticity": "syndiotactic"}, dict(chains=10, density=0.8, contact_scale=A)),
    # Rosenbluth methods
    ("PE rosenbluth 10x100", "chains", {"units": [unit(PE)], "dp": 100, "method": "rosenbluth"}, dict(chains=10, density=0.8, contact_scale=A)),
    ("PE rosenbluth_lj 10x100", "chains", {"units": [unit(PE)], "dp": 100, "method": "rosenbluth_lj"}, dict(chains=10, density=0.8, contact_scale=A)),
    ("PS rosenbluth 10x30", "chains", {"units": [unit(PS)], "dp": 30, "method": "rosenbluth"}, dict(chains=10, density=0.8, contact_scale=A)),
    ("PS rosenbluth_lj 10x30", "chains", {"units": [unit(PS)], "dp": 30, "method": "rosenbluth_lj"}, dict(chains=10, density=0.8, contact_scale=A)),
    # architectures, regions, orientation, look-ahead, a ring backbone
    ("PE star 4 arms 6x25", "chains", {"units": [unit(PE)], "dp": 25, "architecture": "star", "arms": 4}, dict(chains=6, density=0.7, contact_scale=A)),
    ("PE comb 6x60 arm 6", "chains", {"units": [unit(PE)], "dp": 60, "architecture": "comb", "arm_dp": 6, "spacing": 6},
     dict(chains=6, density=0.7, contact_scale=A)),
    ("PE slab 10x60", "chains", {"units": [unit(PE)], "dp": 60, "region": {"shape": "slab", "thickness": 25, "vacuum": 20}},
     dict(chains=10, density=0.8, contact_scale=A)),
    ("PE oriented z 10x60", "chains", {"units": [unit(PE)], "dp": 60, "orientation": {"axis": "z", "strength": 4}}, dict(chains=10, density=0.8, contact_scale=A)),
    ("PS lookahead 3 10x30", "chains", {"units": [unit(PS)], "dp": 30, "lookahead": 3}, dict(chains=10, density=0.8, contact_scale=A)),
    ("ENB 12x10 0.3", "chains", {"units": [unit(ENB)], "dp": 10}, dict(chains=12, density=0.3, contact_scale=A, seed=7)),
    # fixed full contact limits (no easing): the engine itself
    ("PE fixed 10x100 @0.5", "chains", {"units": [unit(PE)], "dp": 100}, dict(chains=10, density=0.5)),
    ("PS fixed 10x30 @0.4", "chains", {"units": [unit(PS)], "dp": 30}, dict(chains=10, density=0.4)),
    # the Studio's path: live view on, auto contact scale
    ("PS Studio live 15x60", "live", {"units": [unit(PS)], "dp": 60}, dict(chains=15, density=0.8, contact_scale=A)),
    # the polystyrene grower (caps grow without --units)
    ("caps grow PS 1x20", "ps", None, dict(chains=1, dp=20, density=0.1, contact_scale=A)),
    ("caps grow PS 10x30 @0.5", "ps", None, dict(chains=10, dp=30, density=0.5, contact_scale=A)),
    ("caps grow PS 15x60 @0.5", "ps", None, dict(chains=15, dp=60, density=0.5, contact_scale=A)),
]


def read_lammps(path):
    """atoms (id → (x, y, z)), bonds, box edges from a LAMMPS data file."""
    lines = Path(path).read_text().splitlines()
    box = [0, 0, 0]
    pos, bonds, sec = {}, [], None
    for ln in lines:
        t = ln.split("#")[0].split()
        if not t:
            continue
        if len(t) >= 4 and t[2].endswith("lo") and t[3].endswith("hi"):
            box["xyz".index(t[2][0])] = float(t[1]) - float(t[0])
            continue
        if t[0] in ("Atoms", "Bonds", "Masses", "Velocities", "Angles", "Dihedrals", "Impropers") or t[0].endswith("Coeffs"):
            sec = t[0]
            continue
        if sec == "Atoms" and len(t) >= 7:
            pos[int(t[0])] = (float(t[4]), float(t[5]), float(t[6])) if len(t) >= 7 else None
        elif sec == "Bonds" and len(t) >= 4:
            bonds.append((int(t[2]), int(t[3])))
    return pos, bonds, box


def closest_contact(pos, bonds, box):
    """Smallest distance between atoms more than three bonds apart (minimum image), and the number of such pairs under 1.6 Å."""
    ids = sorted(pos)
    idx = {a: k for k, a in enumerate(ids)}
    n = len(ids)
    adj = [[] for _ in range(n)]
    for a, b in bonds:
        adj[idx[a]].append(idx[b]), adj[idx[b]].append(idx[a])
    X = [pos[a] for a in ids]
    rc = 3.0
    nc = [max(1, int(L // rc)) for L in box]
    cs = [box[k] / nc[k] for k in range(3)]
    bins = {}
    key = []
    for i, p in enumerate(X):
        k = tuple(int((p[d] % box[d]) // cs[d]) % nc[d] for d in range(3))
        key.append(k)
        bins.setdefault(k, []).append(i)
    best, close = 1e9, 0
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
                        if d2 < 1.6 * 1.6:
                            close += 1
    return best, close


def run_case(L, case, seed, check, tmp):
    name, kind, spec, opt = case
    g = GrowOpts(opt.get("chains", 1), opt.get("dp", 0), 0, opt.get("seed", seed), opt.get("box", 0.0), opt.get("density", 0.0),
                 opt.get("contact_scale", 1.0), 1)
    rep = C.create_string_buffer(8192)
    snaps = [0]

    @LIVE
    def live(doc, js, user):
        snaps[0] += 1
        L.caps_close(doc)

    t0, c0 = time.perf_counter(), time.process_time()
    if kind == "ps":
        h = L.caps_grow(C.byref(g), None, None, rep, len(rep))
    elif kind == "live":
        h = L.caps_grow_chains_live(json.dumps(spec).encode(), C.byref(g), None, live, None, rep, len(rep))
    else:
        h = L.caps_grow_chains(json.dumps(spec).encode(), C.byref(g), None, None, rep, len(rep))
    dt, dc = time.perf_counter() - t0, time.process_time() - c0
    if not h:
        return {"name": name, "time": dt, "cpu": dc, "error": L.caps_last_error().decode()}
    s = Summary()
    L.caps_summary_get(h, C.byref(s))
    notes = rep.value.decode()
    eased = next((ln.split(" units placed")[0] for ln in notes.splitlines() if "units placed with eased" in ln), "")
    out = {"name": name, "time": dt, "cpu": dc, "atoms": s.atoms, "bonds": s.bonds, "density": s.density, "eased": eased}
    if check:
        f = os.path.join(tmp, "g.data")
        L.caps_save(h, f.encode())
        pos, bonds, box = read_lammps(f)
        out["closest"], out["close_pairs"] = closest_contact(pos, bonds, box)
    L.caps_close(h)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--only", default="", help="cases whose name contains this text")
    ap.add_argument("--no-check", action="store_true", help="skip the closest-contact check")
    ap.add_argument("--timeout", type=float, default=300, help="seconds per case (each runs in its own process)")
    ap.add_argument("--out", default="", help="also write the table as TSV")
    ap.add_argument("--case", type=int, default=-1, help=argparse.SUPPRESS)
    a = ap.parse_args()
    if a.case >= 0:
        with tempfile.TemporaryDirectory() as tmp:
            print(json.dumps(run_case(load(), CASES[a.case], a.seed, not a.no_check, tmp)))
        return
    rows = []
    print(f"{'case':28s} {'wall s':>8s} {'cpu s':>8s} {'atoms':>6s} {'g/cm3':>6s} {'closest':>7s} {'<1.6':>5s}  eased units", flush=True)
    for i, case in enumerate(CASES):
        if a.only and a.only not in case[0]:
            continue
        cmd = [sys.executable, __file__, "--case", str(i), "--seed", str(a.seed)] + (["--no-check"] if a.no_check else [])
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=a.timeout)
            res = json.loads(p.stdout.strip().splitlines()[-1]) if p.stdout.strip() else {"name": case[0], "time": float("nan"), "cpu": float("nan"),
                                                                                           "error": p.stderr.strip()[-200:]}
        except subprocess.TimeoutExpired:
            res = {"name": case[0], "time": a.timeout, "cpu": float("nan"), "error": f"time-out (over {a.timeout:.0f} s)"}
        rows.append(res)
        if "error" in res:
            print(f"{res['name']:28s} {res['time']:8.2f} {res['cpu']:8.2f}  {res['error'][:90]}", flush=True)
        else:
            print(f"{res['name']:28s} {res['time']:8.2f} {res['cpu']:8.2f} {res['atoms']:6d} {res['density']:6.3f} {res.get('closest', float('nan')):7.3f} "
                  f"{res.get('close_pairs', -1):5d}  {res.get('eased', '')}", flush=True)
    if a.out:
        with open(a.out, "w") as f:
            f.write("case\twall_s\tcpu_s\tatoms\tdensity\tclosest\tclose_pairs\teased\terror\n")
            for r in rows:
                f.write(f"{r['name']}\t{r['time']:.3f}\t{r['cpu']:.3f}\t{r.get('atoms', '')}\t{r.get('density', '')}\t{r.get('closest', '')}\t"
                        f"{r.get('close_pairs', '')}\t{r.get('eased', '')}\t{r.get('error', '')}\n")


if __name__ == "__main__":
    main()
