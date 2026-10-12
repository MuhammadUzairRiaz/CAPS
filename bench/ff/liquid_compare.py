#!/usr/bin/env python3
"""Liquid results against experiment, and the job list that produces them.

    liquid_compare.py jobs SET.json FF.json OUTDIR [--threads 8] [--eq 100000] [--prod 200000] [--gas 1000000]
        one line per molecule of SET.json (bench/ff/open/liquids_oplsaa.json): the liquid_validate.py command, with
        its JSON result written to OUTDIR/NAME.json (run them as a job array)
    liquid_compare.py table OUTDIR REFERENCE.json [--csv out.csv]
        density and heat of vaporisation of every finished molecule against the reference values (NIST WebBook /
        ThermoML only), the errors, and a flag where the density is off by more than 2 % or dHvap by more than
        0.5 kcal/mol: the tolerances OPLS-AA's own liquid fits worked to (Jorgensen, Maxwell & Tirado-Rives 1996)

A reference flagged as weak in REFERENCE.json (a dHvap away from 298 K or from a compilation) is reported but does
not raise a flag.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOL_RHO, TOL_H = 0.02, 0.5


def arg(k, d=None):
    return sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


def slug(name):
    return "".join(c if c.isalnum() else "_" for c in name).strip("_")


def jobs():
    S, ff, out = json.load(open(sys.argv[2])), os.path.abspath(sys.argv[3]), os.path.abspath(sys.argv[4])
    th, eq, prod, gas = arg("--threads", "8"), arg("--eq", "100000"), arg("--prod", "200000"), arg("--gas", "1000000")
    script = os.path.join(HERE, "liquid_validate.py")
    for m in S["molecules"]:
        s = slug(m["name"])
        print(f"python3 {script} {s} '{m['smiles']}' --ff {ff} --out {out}/{s} --n {m['n']} --temp {m['T']} "
              f"--pressure {m['P_atm']} --rho {m['rho_start']} --eq {eq} --prod {prod} --gas {gas} --threads {th} "
              f"--json {out}/{s}.json")


def weak(r):
    """a reference value that is reported but not used to flag: dHvap away from 298 K, or from a compilation"""
    T = r.get("dHvap_T_K") or 298.15
    note = (r.get("note") or "").lower()
    return abs(T - 298.15) > 5 or "compilation" in note or "too low" in note


def table():
    out, ref = sys.argv[2], {slug(r["name"]): r for r in json.load(open(sys.argv[3]))}
    rows = []
    for k, r in sorted(ref.items()):
        p = os.path.join(out, k + ".json")
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        rho_x, h_x = r.get("density_g_cm3"), r.get("dHvap_kcal_mol")
        e_rho = (d["density"] - rho_x) / rho_x if rho_x else None
        e_h = d["dHvap"] - h_x if h_x else None
        flags = []
        if e_rho is not None and abs(e_rho) > TOL_RHO:
            flags.append("density")
        if e_h is not None and abs(e_h) > TOL_H and not weak(r):
            flags.append("dHvap")
        rows.append({"molecule": r["name"], "rho_calc": d["density"], "rho_err_calc": d["density_err"], "rho_exp": rho_x,
                     "rho_dev_pct": None if e_rho is None else round(100 * e_rho, 2),
                     "dHvap_calc": d["dHvap"], "dHvap_err_calc": d["dHvap_err"], "dHvap_exp": h_x,
                     "dHvap_dev": None if e_h is None else round(e_h, 3), "weak_reference": weak(r), "flag": " ".join(flags)})
    f = lambda v, fmt: "" if v is None else format(v, fmt)
    print(f"{'molecule':22s} {'rho':>7} {'exp':>7} {'dev%':>6}  {'dHvap':>6} {'exp':>6} {'dev':>6}  flag")
    for x in rows:
        print(f"{x['molecule'][:22]:22s} {f(x['rho_calc'], '.4f'):>7} {f(x['rho_exp'], '.4f'):>7} {f(x['rho_dev_pct'], '+.2f'):>6}  "
              f"{f(x['dHvap_calc'], '.2f'):>6} {f(x['dHvap_exp'], '.2f'):>6} {f(x['dHvap_dev'], '+.2f'):>6}  "
              f"{x['flag']}{' (weak ref)' if x['weak_reference'] else ''}")
    dr = [abs(x["rho_dev_pct"]) for x in rows if x["rho_dev_pct"] is not None]
    dh = [abs(x["dHvap_dev"]) for x in rows if x["dHvap_dev"] is not None and not x["weak_reference"]]
    if dr:
        print(f"density: mean |dev| {sum(dr) / len(dr):.2f} % over {len(dr)} liquids")
    if dh:
        print(f"dHvap: mean |dev| {sum(dh) / len(dh):.2f} kcal/mol over {len(dh)} liquids (weak references left out)")
    if arg("--csv"):
        with open(arg("--csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["molecule"])
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    {"jobs": jobs, "table": table}[sys.argv[1]]()
