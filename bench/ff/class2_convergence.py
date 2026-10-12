#!/usr/bin/env python3
"""CAPS against LAMMPS on the same data file, for PCFF and CVFF (from their .frc files) on four small melts: minimisation
(energy against iteration) and NPT at 300 K (density and energy against time).

Each melt (PE, PS, PMMA, cis-1,4-polyisoprene; under 3,000 atoms) is grown by CAPS at a low density, then for each force
field:

  · minimisation: caps ff apply writes the LAMMPS data file and input (lj/class2 or lj/cut with coul/dsf 0.2 Å⁻¹, 10 Å, tail
    correction, CAPS's own settings for relax and md); LAMMPS minimises it with min_style cg (etol 0, thermo every
    iteration), CAPS with caps relax --method cg and --method lbfgs (no push-off, the energy reported every 8 iterations);
    each final structure is evaluated again in LAMMPS;
  · NPT: caps relax --density packs the grown melt to the polymer's density; LAMMPS draws the velocities (300 K, Gaussian,
    zero momentum) and writes them into the data file both programs start from; LAMMPS fix npt (Nosé–Hoover chains, MTK,
    τ_T 100 fs, τ_P 1000 fs, iso 1 atm) and caps md --thermostat nose-hoover --barostat mtk with the same τ, 1 fs steps.

usage: class2_convergence.py --out DIR [--steps 5000] [--min-iter 2000] [--only PE] [--ff pcff-frc,cvff-frc] [-j 4] [--again]
(a case with its result.json is kept unless --again)
Needs LMP (default ~/lammps/build-class2/lmp) with CLASS2 and EXTRA-PAIR (coul/dsf). Local limits: ≤ 5000 MD steps.
"""
import concurrent.futures as cf, csv, json, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
OUT = os.path.abspath(arg("--out", "class2_convergence"))
STEPS = int(arg("--steps", "5000"))
MIN_IT = int(arg("--min-iter", "2000"))
FFS = arg("--ff", "pcff-frc,cvff-frc").split(",")
# name, repeat unit, chains, units per chain, grow density, packing density (g/cm³, the melt near 300 K)
MELTS = [("PE", "[*]CC[*]", 10, 48, 0.7, 0.85),
         ("PS", "[*]CC([*])C1=CC=CC=C1", 12, 15, 0.6, 1.04),
         ("PMMA", "[*]CC(C)(C(=O)OC)[*]", 10, 19, 0.6, 1.17),
         ("cisPI", "[*]C/C=C(C)\\C[*]", 10, 22, 0.6, 0.91)]
if "--only" in sys.argv:
    MELTS = [m for m in MELTS if m[0] in arg("--only", "").split(",")]
SEED = 4928459


def sh(cmd, cwd, log=None):
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if log:
        open(os.path.join(cwd, log), "w").write(r.stdout + "\n--- stderr ---\n" + r.stderr)
    if r.returncode not in (0, 3):   # caps relax returns 3 when it stops by iterations
        raise RuntimeError(f"{' '.join(cmd[:3])} failed in {cwd}: {(r.stderr or r.stdout).strip()[-400:]}")
    return r


def ff_json(ff):
    return os.path.join(ROOT, "data", "forcefields", ff + ".json")


def lammps_input(base_in, body, out):
    """CAPS's LAMMPS input with its run line replaced by `body`"""
    txt = open(base_in).read()
    txt = re.sub(r"(?m)^(thermo_style|thermo_modify|run|timestep) .*$", "", txt)
    open(out, "w").write(txt + "\n" + body + "\n")


def thermo_rows(log, first_col="Step"):
    """every thermo block of a LAMMPS log as dicts"""
    rows, head = [], None
    for line in open(log):
        w = line.split()
        if w and w[0] == first_col:
            head = w
            continue
        if head and len(w) == len(head):
            try:
                rows.append(dict(zip(head, map(float, w))))
                continue
            except ValueError:
                pass
        if head and (line.startswith("Loop time") or not w):
            head = None
    return rows


def lammps_point(d, base_in, data, ff):
    """the potential energy of a structure CAPS wrote, in LAMMPS with CAPS's settings: written again by caps ff apply in
    the force field's own styles (caps relax writes CVFF torsions as fourier), its charges kept"""
    stem = os.path.splitext(data)[0] + "_native"
    sh([CAPS, "ff", "apply", data, "--ff", ff_json(ff), "--charges", "keep", "--kspace", "dsf", "--lammps-cutoff", "10",
        "-o", stem + ".data", "--lammps-input", stem + ".in"], d)
    data = stem + ".data"
    lammps_input(base_in, "thermo_style custom step pe\nthermo_modify format float %.10f\nrun 0", os.path.join(d, "point.in"))
    txt = open(os.path.join(d, "point.in")).read()
    txt = re.sub(r"(?m)^read_data .*$", f"read_data {data}", txt)
    open(os.path.join(d, "point.in"), "w").write(txt)
    subprocess.run([LMP, "-in", "point.in", "-log", "point.log", "-screen", "none"], cwd=d)
    rows = thermo_rows(os.path.join(d, "point.log"))
    if not rows:
        raise RuntimeError(f"LAMMPS could not evaluate {data} in {d}")
    return rows[-1]["PotEng"]


def caps_relax_trace(err):
    """caps relax's progress (stderr): iteration and energy, every 8 iterations (the k-th report is after iteration k + 1)"""
    out = []
    for part in open(err).read().replace("\r", "\n").splitlines():
        m = re.search(r"it\s+(\d+)\s+E\s+(-?[\d.]+)\s+\|F\|max\s+([\d.]+)", part)
        if m:
            out.append((int(m.group(1)) + 1, float(m.group(2)), float(m.group(3))))
    return out


def one(melt, ff):
    name, smi, ch, dp, rho0, rho = melt
    d = os.path.join(OUT, f"{name}_{ff}")
    os.makedirs(d, exist_ok=True)
    grown = os.path.join(OUT, f"{name}.grown.data")
    # --- minimisation from the grown melt ---
    sh([CAPS, "ff", "apply", grown, "--ff", ff_json(ff), "--charges", "auto", "--kspace", "dsf", "--lammps-cutoff", "10",
        "-o", "start.data", "--lammps-input", "start.in"], d, "apply.log")
    lammps_input(os.path.join(d, "start.in"),
                 "thermo_style custom step pe fmax fnorm\nthermo_modify format float %.10f\nthermo 1\nmin_style cg\n"
                 f"minimize 0.0 1.0e-12 {MIN_IT} {100 * MIN_IT}", os.path.join(d, "min_lammps.in"))
    subprocess.run([LMP, "-in", "min_lammps.in", "-log", "min_lammps.log", "-screen", "none"], cwd=d)
    lam = thermo_rows(os.path.join(d, "min_lammps.log"))
    evals = re.search(r"Iterations, force evaluations\s*=\s*\d+\s+(\d+)", open(os.path.join(d, "min_lammps.log")).read())
    e0 = lam[0]["PotEng"]
    res = {"melt": name, "forcefield": ff, "atoms": None, "E0": e0, "lammps_cg": {"iterations": len(lam) - 1,
           "force_evaluations": int(evals.group(1)) if evals else None, "E_final": lam[-1]["PotEng"]}}
    with open(os.path.join(d, "min_trace.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["engine", "method", "iteration", "energy_kcal_mol", "fmax_kcal_mol_A"])
        for r in lam:
            w.writerow(["lammps", "cg", int(r["Step"]), f"{r['PotEng']:.6f}", f"{r['Fmax']:.6f}"])
        for method in ("cg", "lbfgs"):
            err = os.path.join(d, f"min_caps_{method}.err")
            r = subprocess.run([CAPS, "relax", "start.data", "--ff", ff_json(ff), "--method", method, "--iterations", str(MIN_IT),
                                "--ftol", "1e-6", "--no-pushoff", "--cutoff", "10", "-o", f"min_caps_{method}.data"],
                               cwd=d, capture_output=True, text=True)
            open(err, "w").write(r.stderr)
            open(os.path.join(d, f"min_caps_{method}.log"), "w").write(r.stdout)
            tr = caps_relax_trace(err)
            w.writerow(["caps", method, 0, f"{e0:.6f}", ""])
            for it, e, fm in tr:
                w.writerow(["caps", method, it, f"{e:.1f}", f"{fm:.3f}"])
            ef = lammps_point(d, os.path.join(d, "start.in"), f"min_caps_{method}.data", ff)
            ev = re.search(r"(\d+) iterations, (\d+) energy and force evaluations", r.stdout)
            res["caps_" + method] = {"iterations": int(ev.group(1)) if ev else (tr[-1][0] if tr else 0),
                                     "force_evaluations": int(ev.group(2)) if ev else None, "E_final": ef,
                                     "stopped_by": (re.findall(r"\((force|energy|iterations|line search)\)", r.stdout) or ["?"])[-1]}
    # --- NPT from the packed melt ---
    sh([CAPS, "relax", grown, "--ff", ff_json(ff), "--density", str(rho), "--cutoff", "10", "-o", "packed.data"], d, "pack.log")
    sh([CAPS, "ff", "apply", "packed.data", "--ff", ff_json(ff), "--charges", "keep", "--kspace", "dsf", "--lammps-cutoff", "10",
        "-o", "npt_start.data", "--lammps-input", "npt_start.in"], d, "apply_npt.log")
    lammps_input(os.path.join(d, "npt_start.in"), f"velocity all create 300.0 {SEED} dist gaussian mom yes rot no\n"
                 "write_data npt_start_v.data\nrun 0", os.path.join(d, "vel.in"))
    subprocess.run([LMP, "-in", "vel.in", "-log", "vel.log", "-screen", "none"], cwd=d)
    lammps_input(os.path.join(d, "npt_start.in"),
                 "timestep 1.0\nthermo_style custom step time temp pe ke etotal press density vol\nthermo_modify format float %.8f\n"
                 "thermo 50\nfix 1 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0\n" f"run {STEPS}", os.path.join(d, "npt_lammps.in"))
    txt = open(os.path.join(d, "npt_lammps.in")).read().replace("read_data       npt_start.data", "read_data npt_start_v.data")
    txt = re.sub(r"(?m)^read_data .*$", "read_data npt_start_v.data", txt)
    open(os.path.join(d, "npt_lammps.in"), "w").write(txt)
    subprocess.run([LMP, "-in", "npt_lammps.in", "-log", "npt_lammps.log", "-screen", "none"], cwd=d)
    lt = thermo_rows(os.path.join(d, "npt_lammps.log"))
    sh([CAPS, "md", "npt_start_v.data", "--ff", ff_json(ff), "--charges", "keep", "--steps", str(STEPS), "--dt", "1", "--temp", "300",
        "--thermostat", "nose-hoover", "--tau-t", "100", "--barostat", "mtk", "--pressure", "1", "--tau-p", "1000", "--cutoff", "10",
        "--precision", "double", "--thermo", "50", "--log", "npt_caps.csv", "-o", "npt_caps.data"], d, "npt_caps.log")
    with open(os.path.join(d, "npt_trace.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["engine", "step", "time_ps", "temperature_K", "potential_kcal_mol", "total_kcal_mol", "pressure_atm", "density_g_cm3"])
        for r in lt:
            w.writerow(["lammps", int(r["Step"]), r["Time"] / 1000.0, r["Temp"], r["PotEng"], r["TotEng"], r["Press"], r["Density"]])
        for r in csv.DictReader(open(os.path.join(d, "npt_caps.csv"))):
            w.writerow(["caps", int(r["step"]), float(r["time_ps"]), float(r["temperature_K"]), float(r["potential"]), float(r["total"]),
                        float(r["pressure_atm"]), float(r["density_g_cm3"])])
    res["atoms"] = int(re.search(r"(\d+) atoms", open(os.path.join(d, "start.data")).read()).group(1))
    res["velocities_from_file"] = "velocities drawn" not in open(os.path.join(d, "npt_caps.log")).read()
    json.dump(res, open(os.path.join(d, "result.json"), "w"), indent=1)
    return res


def main():
    os.makedirs(OUT, exist_ok=True)
    for name, smi, ch, dp, rho0, rho in MELTS:
        grown = os.path.join(OUT, f"{name}.grown.data")
        if not os.path.exists(grown):
            sh([CAPS, "grow", "--units", smi, "--chains", str(ch), "--dp", str(dp), "--density", str(rho0), "--seed", "7",
                "-o", grown], OUT, f"{name}.grow.log")
    jobs = [(m, ff) for m in MELTS for ff in FFS]

    def job(a):
        done = os.path.join(OUT, f"{a[0][0]}_{a[1]}", "result.json")
        if os.path.exists(done) and "--again" not in sys.argv:   # finished earlier: kept
            return json.load(open(done))
        try:
            return one(*a)
        except Exception as ex:   # one case failing leaves the others running
            return {"melt": a[0][0], "forcefield": a[1], "error": str(ex)}
    with cf.ThreadPoolExecutor(int(arg("-j", "4"))) as ex:
        results = list(ex.map(job, jobs))
    json.dump(results, open(os.path.join(OUT, "convergence.json"), "w"), indent=1)
    for r in results:
        print(json.dumps(r))


if __name__ == "__main__":
    main()
