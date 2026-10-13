#!/usr/bin/env python3
"""Liquid density and heat of vaporisation of small molecules with CAPS's own engine: the test OPLS-AA's Lennard-Jones
terms were fitted to (Jorgensen, Maxwell & Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)).

    liquid_validate.py NAME SMILES --ff FF.json --out DIR [--n 300] [--temp 298.15] [--eq 5000] [--prod 5000]
                       [--gas 20000] [--rho 0.8] [--pressure 1] [--overlay OV.json] [--threads N] [--json result.json]

liquid   N molecules packed at 0.9 of --rho (a guess, or the experimental density: only the start), typed with FF.json (its own charges), minimised, then NPT at
         --temp and --pressure atm (Bussi thermostat, C-rescale barostat, PME, the force field's own cut-off (OPLS-AA 2024:
         12 A) with tail corrections), --eq steps
         to equilibrate and --prod steps averaged: density and potential energy per molecule. caps ff apply writes the
         force field (and --overlay) into the data files; relax and md take it, its charges, 1-4 weights and cut-off from them
gas      one molecule in a 100 A box, Langevin NVT, every intramolecular pair inside the 45 A cutoff and no periodic
         image (Coulomb as damped shifted force with alpha 0, the plain shifted Coulomb; for a molecule under 10 A
         the shift changes its energy by far less than the 0.5 kcal/mol tolerance)
result   dHvap = <E_gas> - <E_liq>/N + RT, in kcal/mol, with block errors (5 blocks) of both averages

Local checks stay short (<= 5000 steps); production runs go to the cluster. No system above 15,000 atoms.
"""
import csv, json, math, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.environ.get("CAPS_BIN", os.path.join(ROOT, "build", "cli", "caps"))
R = 0.0019872043   # kcal/mol/K
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Si": 28.086, "P": 30.974, "S": 32.06,
        "Cl": 35.45, "Br": 79.904, "I": 126.904}


def arg(k, d=None):
    return sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


def run(cmd, log):
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()   # before the command's own output
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode:
        raise SystemExit(f"failed: {' '.join(cmd)} (see {log})")


def thermo(path):
    rows = list(csv.DictReader(open(path)))
    return rows


def col(rows, *names):
    for n in names:
        if rows and n in rows[0]:
            return [float(r[n]) for r in rows]
    raise SystemExit(f"none of {names} in the thermo log (columns: {list(rows[0]) if rows else []})")


def block(v, nb=5):
    """mean and standard error from nb block averages"""
    m = sum(v) / len(v)
    k = len(v) // nb
    if k < 1:
        return m, float("nan")
    b = [sum(v[i * k:(i + 1) * k]) / k for i in range(nb)]
    mb = sum(b) / nb
    return m, math.sqrt(sum((x - mb) ** 2 for x in b) / (nb * (nb - 1)))


def main():
    name, smiles = sys.argv[1], sys.argv[2]
    ff, out = os.path.abspath(arg("--ff")), os.path.abspath(arg("--out"))
    n, T = int(arg("--n", "300")), float(arg("--temp", "298.15"))
    eq, prod, gas = int(arg("--eq", "5000")), int(arg("--prod", "5000")), int(arg("--gas", "20000"))
    threads = arg("--threads")
    th = ["--threads", threads] if threads else []
    ov = ["--overlay", os.path.abspath(arg("--overlay"))] if arg("--overlay") else []   # e.g. a Lennard-Jones refit on trial
    os.makedirs(out, exist_ok=True)
    log = os.path.join(out, "run.log")
    mol = os.path.join(out, "mol.mol2")   # with its bonds: nothing in the box is perceived from distances
    run([CAPS, "build", smiles, "-o", mol], log)
    el = []
    sec = None
    for l in open(mol):
        if l.startswith("@<TRIPOS>"):
            sec = l.strip()[9:]
        elif sec == "ATOM" and l.split():
            el.append("".join(c for c in l.split()[1] if c.isalpha()))
    if n * len(el) > 15000:
        raise SystemExit(f"{n} molecules of {len(el)} atoms exceed 15,000 atoms")
    mw = sum(MASS[e] for e in el)
    # packed at 0.9 x --rho (the barostat takes it from there), periodic (pbc: distances across the faces checked)
    edge = (n * mw / 0.6022 / (0.9 * float(arg("--rho", "0.8")))) ** (1 / 3)

    def pack(count, L, path):
        inp = path + ".inp"
        open(inp, "w").write(f"tolerance 2.0\npbc 0. 0. 0. {L:.4f} {L:.4f} {L:.4f}\noutput {path}\nstructure {mol}\n"
                             f"  number {count}\n  inside box 0. 0. 0. {L:.4f} {L:.4f} {L:.4f}\nend structure\n")
        run([CAPS, "pack", inp, "--quiet"], log)

    box = os.path.join(out, "box.data")
    pack(n, edge, box)
    liq = os.path.join(out, "liq.data")
    run([CAPS, "ff", "apply", box, "--ff", ff, "--charges", "types", "-o", liq] + ov, log)
    # relax and md run on the force field the data file holds (written by caps ff apply with the overlay), at its own
    # cut-off; a file whose force field CAPS cannot rebuild stops them
    run([CAPS, "relax", liq, "-o", os.path.join(out, "liq_min.data")], log)
    md = ["--temp", str(T), "--thermostat", "bussi", "--barostat", "crescale", "--pressure", arg("--pressure", "1"), "--pme",
          "--dt", "1", "--thermo", "100"] + th
    run([CAPS, "md", os.path.join(out, "liq_min.data"), "-o", os.path.join(out, "liq_eq.data"), "--steps", str(eq),
         "--log", os.path.join(out, "eq.csv")] + md, log)
    run([CAPS, "md", os.path.join(out, "liq_eq.data"), "-o", os.path.join(out, "liq_prod.data"), "--steps", str(prod),
         "--log", os.path.join(out, "prod.csv")] + md, log)
    # gas: one molecule in a 100 A box
    g = os.path.join(out, "gas_box.data")
    pack(1, 100.0, g)
    gd = os.path.join(out, "gas.data")
    run([CAPS, "ff", "apply", g, "--ff", ff, "--charges", "types", "-o", gd] + ov, log)
    run([CAPS, "md", gd, "-o", os.path.join(out, "gas_out.data"), "--steps", str(gas), "--temp", str(T),
         "--thermostat", "langevin", "--barostat", "none", "--cutoff", "45", "--dsf-alpha", "0", "--dt", "1",
         "--thermo", "50", "--log", os.path.join(out, "gas.csv")] + th, log)
    P = thermo(os.path.join(out, "prod.csv"))
    rho, rho_e = block(col(P, "density", "density_g_cm3", "rho"))
    El, El_e = block(col(P, "pe", "potential", "epot", "E_pot"))
    Gs = thermo(os.path.join(out, "gas.csv"))
    gpe = col(Gs, "pe", "potential", "epot", "E_pot")
    gpe = gpe[len(gpe) // 5:]   # the first fifth relaxes the built geometry
    Eg, Eg_e = block(gpe)
    dh = Eg - El / n + R * T
    dh_e = math.sqrt(Eg_e ** 2 + (El_e / n) ** 2)
    res = {"name": name, "smiles": smiles, "ff": os.path.basename(ff) + (" + " + os.path.basename(arg("--overlay")) if arg("--overlay") else ""), "n": n, "atoms": n * len(el), "T": T,
           "P_atm": float(arg("--pressure", "1")),
           "steps": {"eq": eq, "prod": prod, "gas": gas}, "density": round(rho, 4), "density_err": round(rho_e, 4),
           "E_liq_per_mol": round(El / n, 4), "E_gas": round(Eg, 4), "dHvap": round(dh, 3), "dHvap_err": round(dh_e, 3)}
    print(json.dumps(res))
    if arg("--json"):
        json.dump(res, open(arg("--json"), "w"), indent=1)


if __name__ == "__main__":
    main()
