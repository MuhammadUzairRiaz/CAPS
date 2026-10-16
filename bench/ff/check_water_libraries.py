#!/usr/bin/env python3
"""The water and gas libraries built from their sources (bench/ff/make_water_models.py, make_noble_gases.py) against
their papers, LAMMPS and GROMACS.

For each water library: a box of waters in the model's own geometry (r_OH, H-O-H) is typed by the library, which adds
the massless sites; then
  1. the model as published: every charge, the neutral molecule, r_OH, H-O-H, d_OM (four sites), d_OL and L-O-L (TIP5P)
     measured on the typed structure, and the Lennard-Jones the library assigns;
  2. LAMMPS (run 0 of the exported files: lj/cut/tip4p/long and pppm/tip4p for four sites) against CAPS's energy (PME);
  3. GROMACS (grompp, then mdrun -rerun of the exported files; M and the lone pairs as virtual sites) against CAPS's
     energy, the dispersion correction left out of both (their definitions differ);
  4. with --density: GROMACS NPT at 298 K and 1 bar, rigid (constraints = h-angles), 512 waters, 5000 steps of 2 fs, the mean density of
     the last half; a short check (10 ps), not a converged density.
mW (Stillinger-Weber) and the noble gases: LAMMPS against CAPS (no electrostatics).

    PYTHONPATH=data/python CAPS_LIB=build/capi/libcaps.dylib python3 bench/ff/check_water_libraries.py [--density] [--only id]
Needs lmp (~/.local/bin/lmp, KSPACE, MOLECULE, MANYBODY) and gmx.
"""
import ctypes as C
import json, math, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import make_water_models as W
import make_noble_gases as G
import caps

LMP = os.environ.get("LMP", os.path.expanduser("~/.local/bin/lmp"))
GMX = os.environ.get("GMX", shutil.which("gmx") or "gmx")
FF = os.path.join(ROOT, "data", "forcefields")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
only = arg("--only", "")
DENSITY = "--density" in sys.argv
work = arg("--keep", "") or tempfile.mkdtemp(prefix="caps_waterlib_")
L = caps.library()
L.caps_set_electrostatics.argtypes = [C.c_int32, C.c_double, C.c_double, C.c_int32]
L.caps_set_electrostatics(1, 1e-6, 0.6, 6)   # PME, tight
KJ = 4.184
fails = 0


def water_xyz(path, r_oh, theta, n, edge, seed=4):
    """n waters of the given geometry packed in a periodic cube of edge (A)"""
    h = math.radians(theta) / 2
    one = os.path.join(os.path.dirname(path), "one.xyz")
    with open(one, "w") as f:
        f.write(f"3\nwater\nO 0 0 0\nH {r_oh * math.sin(h):.10f} 0 {r_oh * math.cos(h):.10f}\nH {-r_oh * math.sin(h):.10f} 0 {r_oh * math.cos(h):.10f}\n")
    d = caps.pack(molecules=[(one, n)], box=edge, tolerance=2.0, seed=seed)
    return d


def lammps_total(d, out, cutoff=0):
    d.export_engines(out, stem="w", gromacs=False, run="check", kspace_accuracy=1e-7, cutoff=cutoff)
    r = subprocess.run([LMP, "-in", "w.in", "-log", "log.lammps"], cwd=out, capture_output=True, text=True)
    rows = re.findall(r"^\s*0\s+(-?[\d.]+)", open(os.path.join(out, "log.lammps")).read(), re.M) if os.path.exists(os.path.join(out, "log.lammps")) else []
    if not rows:
        raise RuntimeError("LAMMPS: " + (r.stdout + r.stderr).strip().splitlines()[-1][:200])
    return float(rows[-1])


def caps_tail(fid, d):
    """CAPS's (and LAMMPS's pair_modify tail) dispersion term, kcal/mol: (2 pi / V) sum_ab N_a N_b 4 eps_ab [sigma^12 / (9 rc^9)
    - sigma^6 / (3 rc^3)] over the library's types, unlike pairs by its mixing rule"""
    lib = json.load(open(os.path.join(FF, fid + ".json")))
    lj = {p["match"][0]: p["params"] for p in lib["pairs"] if len(p["match"]) == 1}
    rc = lib["cutoff"]
    n = {}
    for i in range(d.atoms):
        a = d.atom(i)
        t = {"O": "OW", "H": "HW"}.get(a["element"], None)
        if t:
            n[t] = n.get(t, 0) + 1
    vol = d.summary()["volume"] if "volume" in d.summary() else None
    if vol is None:
        c = d.summary()["cell"]
        vol = c[0] * c[1] * c[2]
    e = 0.0
    for a, na in n.items():
        for b, nb in n.items():
            ea, sa = lj[a]
            eb, sb = lj[b]
            eps = math.sqrt(ea * eb)
            sig = (sa + sb) / 2 if lib["mixing"] == "arithmetic" else math.sqrt(sa * sb)
            e += na * nb * 2 * math.pi * 4 * eps * (sig ** 12 / (9 * rc ** 9) - sig ** 6 / (3 * rc ** 3))
    return e / vol


def gmx_run(out, args, inp=None):
    return subprocess.run([GMX, "-quiet"] + args, cwd=out, input=inp, capture_output=True, text=True)


def gromacs_terms(d, out):
    d.export_engines(out, stem="w", lammps=False, gromacs=True, run="check", cutoff=0)
    mdp = os.path.join(out, "w.mdp")
    r = gmx_run(out, ["grompp", "-f", "w.mdp", "-c", "w.gro", "-p", "w.top", "-o", "w.tpr", "-maxwarn", "10"])
    if r.returncode:
        raise RuntimeError("grompp: " + [l for l in (r.stdout + r.stderr).splitlines() if l.strip()][-3][:300])
    r = gmx_run(out, ["mdrun", "-s", "w.tpr", "-rerun", "w.gro", "-deffnm", "w", "-nt", "4"])
    if r.returncode:
        raise RuntimeError("mdrun: " + (r.stdout + r.stderr).strip().splitlines()[-1][:300])
    r = gmx_run(out, ["energy", "-f", "w.edr", "-o", "e.xvg"], "Potential\nDisper.-corr.\n\n")
    e = {}
    for l in r.stdout.split("Energy", 1)[-1].splitlines():
        m = re.match(r"^(\S.*?)\s+(-?[\d.]+(?:e[-+]?\d+)?)\s+--", l)
        if m:
            e[m.group(1)] = float(m.group(2)) / KJ
    return e


def density(d, out):
    """GROMACS NPT, 298 K, 1 bar, SETTLE, 5000 x 2 fs: the mean density (kg/m3) of the last half"""
    d.export_engines(out, stem="w", lammps=False, gromacs=True, run="check", cutoff=0)
    mdp = open(os.path.join(out, "w.mdp")).read()
    keep = [l for l in mdp.splitlines() if l.split("=")[0].strip() in ("cutoff-scheme", "coulombtype", "rcoulomb", "rvdw", "vdwtype", "vdw-modifier",
                                                                       "pme-order", "DispCorr", "pbc")]
    run = keep + ["integrator = md", "dt = 0.002", "nsteps = 5000", "nstcalcenergy = 10", "nstenergy = 50", "tcoupl = v-rescale", "tc-grps = System",
                  "tau-t = 0.1", "ref-t = 298", "pcoupl = c-rescale", "tau-p = 1.0", "ref-p = 1.0", "compressibility = 4.5e-5", "constraints = h-angles",
                  "fourier-spacing = 0.12", "gen-vel = yes", "gen-temp = 298", "gen-seed = 7", "nstxout-compressed = 0"]
    open(os.path.join(out, "npt.mdp"), "w").write("\n".join(run) + "\n")
    r = gmx_run(out, ["grompp", "-f", "npt.mdp", "-c", "w.gro", "-p", "w.top", "-o", "npt.tpr", "-maxwarn", "10"])
    if r.returncode:
        raise RuntimeError("grompp npt: " + [l for l in (r.stdout + r.stderr).splitlines() if l.strip()][-3][:300])
    r = gmx_run(out, ["mdrun", "-s", "npt.tpr", "-deffnm", "npt", "-nt", "8"])
    if r.returncode:
        raise RuntimeError("mdrun npt: " + (r.stdout + r.stderr).strip().splitlines()[-1][:300])
    r = gmx_run(out, ["energy", "-f", "npt.edr", "-o", "rho.xvg", "-b", "5"], "Density\n\n")
    m = re.search(r"^Density\s+(-?[\d.]+)\s+(\S+)", r.stdout, re.M)
    return float(m.group(1)), float(m.group(2))


def check_model(fid, m):
    global fails
    n = W.sites(m)
    out = os.path.join(work, fid)
    os.makedirs(out, exist_ok=True)
    nw = 512   # a 24.8 A cube: GROMACS needs more than twice the 10 A cut-off
    edge = (nw * 18.0153 / 0.997 / 0.60221408) ** (1 / 3)
    d = water_xyz(os.path.join(out, "box.xyz"), m["r_oh"], m["theta"], nw, edge)
    rep = d.field.assign(os.path.join(FF, fid + ".json"), charges="forcefield")
    if not rep.get("complete"):
        print(f"{fid}: assignment incomplete: {rep.get('missing')}")
        fails += 1
        return
    at = [d.atom(i) for i in range(d.atoms)]
    pos = [a["position"] for a in at]
    # one molecule: O, H, H and its sites in order
    o = next(i for i, a in enumerate(at) if a["element"] == "O")
    mol = [i for i, a in enumerate(at) if a["molecule"] == at[o]["molecule"]]
    hs = [i for i in mol if at[i]["element"] == "H"]
    xs = [i for i in mol if at[i]["element"] not in ("O", "H")]
    P = lambda i: pos[i]
    sub = lambda a, b: [a[k] - b[k] for k in range(3)]
    nrm = lambda v: math.sqrt(sum(x * x for x in v))
    ang = lambda u, v: math.degrees(math.acos(sum(a * b for a, b in zip(u, v)) / nrm(u) / nrm(v)))
    q = sum(at[i]["charge"] for i in mol)
    r1, r2 = nrm(sub(P(hs[0]), P(o))), nrm(sub(P(hs[1]), P(o)))
    th = ang(sub(P(hs[0]), P(o)), sub(P(hs[1]), P(o)))
    errs = [abs(q), abs(r1 - m["r_oh"]), abs(r2 - m["r_oh"]), abs(th - m["theta"]) / 100, abs(at[hs[0]]["charge"] - m["q_h"])]
    geo = f"r_OH {r1:.6f} A, H-O-H {th:.4f} deg, q_H {at[hs[0]]['charge']:.6f}, molecule {q:+.1e} e"
    if n == 4:
        dm = nrm(sub(P(xs[0]), P(o)))
        errs += [abs(dm - m["d_om"]), abs(at[xs[0]]["charge"] + 2 * m["q_h"])]
        geo += f", d_OM {dm:.6f} A, q_M {at[xs[0]]['charge']:.6f}"
    if n == 5:
        d1, d2 = nrm(sub(P(xs[0]), P(o))), nrm(sub(P(xs[1]), P(o)))
        tl = ang(sub(P(xs[0]), P(o)), sub(P(xs[1]), P(o)))
        errs += [abs(d1 - m["d_ol"]), abs(d2 - m["d_ol"]), abs(tl - m["theta_l"]) / 100, abs(at[xs[0]]["charge"] + m["q_h"])]
        geo += f", d_OL {d1:.6f}/{d2:.6f} A, L-O-L {tl:.4f} deg"
    ok_geo = max(errs) < 1e-5 and len(xs) == n - 3
    line = f"{fid:28s} {'ok' if ok_geo else 'GEOMETRY DIFFERS'}: {geo}"
    if n < 5:   # the packed box's contacts eased (a few hundred steps), for single-precision GROMACS; a flexible TIP5P minimised
        # without constraints lets a hydrogen (no Lennard-Jones) fall onto another water's lone pair: compared as packed
        d.relax(ftol=5.0, max_iterations=300)
    caps_e = d.energy()
    if not ok_geo:
        fails += 1
    if n < 5:
        try:
            lt = lammps_total(d, os.path.join(out, "lammps"))
            dl = caps_e["total"] - lt
            line += f"\n{'':28s} LAMMPS {lt:.4f} / CAPS {caps_e['total']:.4f} kcal/mol: {dl:+.4f}"
            if abs(dl) > max(0.02, 1e-5 * abs(lt)):   # LAMMPS's pppm/tip4p at 1e-7: OPC's large M charge leaves 6e-6
                fails += 1
                line += "  DIFFERS"
        except Exception as ex:
            fails += 1
            line += f"\n{'':28s} LAMMPS failed: {ex}"
    try:
        g = gromacs_terms(d, os.path.join(out, "gromacs"))
        gt = g.get("Potential", 0) - g.get("Disper. corr.", 0)
        ct = caps_e["total"] - caps_tail(fid, d)
        dg = ct - gt
        line += f"\n{'':28s} GROMACS {gt:.4f} / CAPS {ct:.4f} kcal/mol (no dispersion correction): {dg:+.4f} ({abs(dg) / abs(gt):.1e} relative)"
        if abs(dg) / abs(gt) > 5e-5:
            fails += 1
            line += "  DIFFERS"
    except Exception as ex:
        fails += 1
        line += f"\n{'':28s} GROMACS failed: {ex}"
    if DENSITY:
        try:
            rho, err = density(d, os.path.join(out, "npt"))
            line += f"\n{'':28s} density, 512 waters, 298 K, 1 bar, 10 ps NPT (GROMACS, O-H and H-O-H constrained): {rho:.1f} +- {err:.1f} kg/m3 (last 5 ps)"
        except Exception as ex:
            line += f"\n{'':28s} density run failed: {ex}"
    print(line, flush=True)


def check_plain(fid, d, label):
    """a library without electrostatics: LAMMPS run 0 against CAPS"""
    global fails
    rep = d.field.assign(os.path.join(FF, fid + ".json"), charges="forcefield")
    if not rep.get("complete"):
        print(f"{fid}: assignment incomplete: {rep.get('missing')}")
        fails += 1
        return
    ce = d.energy()["total"]
    out = os.path.join(work, fid)
    os.makedirs(out, exist_ok=True)
    try:
        lt = lammps_total(d, out, cutoff=10)   # CAPS's single point here cuts at 10 A (the library's own 12 A in runs and exports)
        dl = ce - lt
        bad = abs(dl) > max(0.01, 2e-6 * abs(lt))
        fails += bad
        print(f"{fid:28s} {label}: LAMMPS {lt:.4f} / CAPS {ce:.4f} kcal/mol: {dl:+.4f}{'  DIFFERS' if bad else ''}", flush=True)
    except Exception as ex:
        fails += 1
        print(f"{fid:28s} LAMMPS failed: {ex}", flush=True)


for fid, m in W.MODELS.items():
    if only and only not in fid:
        continue
    check_model(fid, m)
if not only or "mw" in only:
    d = water_xyz(os.path.join(work, "mw.xyz"), 0.9572, 104.52, 343, 21.8)
    check_plain("mw-molinero2009", d, "343 mW sites (Stillinger-Weber)")
if not only or "noble" in only:
    gd = os.path.join(work, "gas")
    os.makedirs(gd, exist_ok=True)
    mols = []
    for el in G.GASES:
        p = os.path.join(gd, el + ".xyz")
        open(p, "w").write(f"1\n{el}\n{el} 0 0 0\n")
        mols.append((p, 100))
    d = caps.pack(molecules=mols, box=30.0, tolerance=3.0, seed=3)
    check_plain(G.FID, d, "100 each of Ne, Ar, Kr, Xe (Lorentz-Berthelot cross terms)")
print(f"{'all agree' if not fails else f'{fails} failed'} (work: {work})")
sys.exit(1 if fails else 0)
