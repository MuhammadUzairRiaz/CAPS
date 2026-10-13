#!/usr/bin/env python3
"""Density and solubility parameter of a library polymer at 298 K with CAPS's own engine, for the validation table of a
force field (the polymers' counterpart of liquid_validate.py).

    polymer_validate.py ID --ff FF.json --out DIR [--atoms 6000] [--dp 0] [--temp 298.15] [--scale 1] [--prod 1000000]
                        [--threads 8] [--seed 1] [--json result.json]

build    chains of the library's repeat unit (data/polymers/library.json) grown with the force field's typing at
         0.5 g/cm3, as many as fit --atoms (DP from --dp, else about 20 repeat units a chain, at least 4 chains)
equil    Larsen's 21-step compression/decompression protocol (caps equilibrate --protocol larsen21, --scale shortens
         every stage) ending at --temp and 1 atm
prod     NPT at --temp and 1 atm (--prod steps, 1 fs, PME), frames every 10 ps
result   density (mean and block error over the production run) and the cohesive energy density and Hildebrand
         solubility parameter delta = sqrt(CED) from the production frames (caps analyze --props density,ced,delta:
         each chain alone against the cell, the same force field)

Production runs belong on the cluster (no system above 15,000 atoms; on a laptop keep --scale and --prod small).
"""
import csv, json, math, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.environ.get("CAPS_BIN", os.path.join(ROOT, "build", "cli", "caps"))


def arg(k, d=None):
    return sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


def run(cmd, log):
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"failed: {' '.join(cmd)} (see {log})")


def unit_atoms(smiles, log, out):
    """atoms of one repeat unit, the two link points as H (a quick count from a capped build)"""
    capped = smiles.replace("[*]", "[H]").replace("*", "[H]")
    xyz = os.path.join(out, "unit.xyz")
    run([CAPS, "build", capped, "-o", xyz], log)
    return int(open(xyz).readline()) - 2


def block(v, nb=5):
    m = sum(v) / len(v)
    k = len(v) // nb
    if k < 1:
        return m, float("nan")
    b = [sum(v[i * k:(i + 1) * k]) / k for i in range(nb)]
    mb = sum(b) / nb
    return m, math.sqrt(sum((x - mb) ** 2 for x in b) / (nb * (nb - 1)))


def main():
    pid = sys.argv[1]
    lib = json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))
    lib = lib.get("polymers", lib) if isinstance(lib, dict) else lib
    P = next(p for p in lib if p["id"] == pid)
    ff, out = os.path.abspath(arg("--ff")), os.path.abspath(arg("--out"))
    T, threads, seed = arg("--temp", "298.15"), arg("--threads", "8"), arg("--seed", "1")
    os.makedirs(out, exist_ok=True)
    log = os.path.join(out, "run.log")
    na = unit_atoms(P["smiles"], log, out)
    target = int(arg("--atoms", "6000"))
    dp = int(arg("--dp", "0")) or max(5, min(40, round(20 * 12 / max(na, 1))) if na > 12 else 20)
    chains = max(4, target // (dp * na + 2))
    if chains * (dp * na + 2) > 15000:
        raise SystemExit("over 15,000 atoms")
    g = os.path.join(out, "grown.data")
    # grown loose (equilibration compresses it); a crowded cell of many short chains retries looser
    for rho0 in ("0.5", "0.35", "0.25"):
        try:
            run([CAPS, "grow", "-o", g, "--units", P["smiles"], "--chains", str(chains), "--dp", str(dp), "--density", rho0,
                 "--ff", ff, "--seed", seed], log)
            break
        except SystemExit:
            if rho0 == "0.25":
                raise
    eq = os.path.join(out, "equil.data")
    run([CAPS, "equilibrate", g, "-o", eq, "--protocol", "larsen21", "--tfinal", T, "--scale", arg("--scale", "1"),
         "--seed", seed, "--log", os.path.join(out, "equil.csv"), "--ff", ff], log)   # --ff: not taken from the data file
    prod = os.path.join(out, "prod.data")
    traj = os.path.join(out, "prod.lammpstrj")
    run([CAPS, "md", eq, "-o", prod, "--steps", arg("--prod", "1000000"), "--dt", "1", "--temp", T, "--thermostat", "bussi",
         "--barostat", "crescale", "--pressure", "1", "--pme", "--cutoff", "12", "--thermo", "1000", "--threads", threads,
         "--log", os.path.join(out, "prod.csv"), "--dump", traj, "--every", "10000", "--ff", ff], log)
    rows = list(csv.DictReader(open(os.path.join(out, "prod.csv"))))
    rho = [float(r["density_g_cm3"]) for r in rows]
    rho = rho[len(rho) // 5:]   # the first fifth relaxes from the equilibration
    rm, re_ = block(rho)
    aj = os.path.join(out, "analyze.json")
    run([CAPS, "analyze", traj, "--topology", eq, "--props", "density,ced,delta", "--ff", ff, "--json", aj], log)
    props = {p["id"]: p for p in json.load(open(aj))}   # caps analyze --json: a list of properties by id
    res = {"id": pid, "name": P["name"], "ff": os.path.basename(ff), "chains": chains, "dp": dp,
           "atoms": chains * (dp * na + 2), "T": float(T), "density": round(rm, 4), "density_err": round(re_, 4),
           "ced": props.get("ced", {}).get("value"), "ced_err": props.get("ced", {}).get("error"),
           "ced_unit": props.get("ced", {}).get("unit"), "delta": props.get("delta", {}).get("value"),
           "delta_err": props.get("delta", {}).get("error"), "delta_unit": props.get("delta", {}).get("unit")}
    print(json.dumps(res))
    if arg("--json"):
        json.dump(res, open(arg("--json"), "w"), indent=1)


if __name__ == "__main__":
    main()
