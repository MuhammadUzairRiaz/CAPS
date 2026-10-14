#!/usr/bin/env python3
"""Polymer cells against experiment: density and solubility parameter of every finished polymer_validate.py run.

    polymer_compare.py OUTDIR REFERENCE.json [--other OUTDIR2 --labels A,B] [--csv out.csv]

REFERENCE.json: a list of {"id", "name", "density_g_cm3", "density_state", "solubility_parameter_MPa05", ...} (one per
library polymer; the sources are kept with the paper's run data). A density more than 2 % off or a solubility parameter
more than 1 MPa^1/2 off is flagged. The cells are amorphous, so a reference marked semicrystalline or crystalline
is listed but not compared (nor counted in the mean). With --other, a second set of runs (another force field, or longer chains) is shown
beside the first on the polymers both finished, each with its mean unsigned deviation over that common set.
"""
import csv, json, os, sys

TOL_RHO, TOL_SP = 0.02, 1.0


def arg(k, d=None):
    return sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


def runs(d):
    out = {}
    for f in sorted(os.listdir(d)):
        if f.endswith(".json"):
            r = json.load(open(os.path.join(d, f)))
            out[r["id"]] = r
    return out


def main():
    A, ref = runs(sys.argv[1]), {r["id"]: r for r in json.load(open(sys.argv[2]))}
    B = runs(arg("--other")) if arg("--other") else None
    la, lb = ((arg("--labels") or "A,B") + ",B").split(",")[:2]
    ids = sorted(k for k in A if k in ref and (B is None or k in B))
    rows = []
    for k in ids:
        r = ref[k]
        row = {"id": k, "polymer": r["name"], "rho_exp": r.get("density_g_cm3"), "rho_state": r.get("density_state"),
               "sp_exp": r.get("solubility_parameter_MPa05")}
        for lab, S in ((la, A), (lb, B)) if B else ((la, A),):
            x = S[k]
            row[f"chains_{lab}"], row[f"dp_{lab}"] = x.get("chains"), x.get("dp")
            row[f"rho_{lab}"], row[f"sp_{lab}"] = x["density"], x.get("delta")
            # a simulated cell is amorphous: a semicrystalline or crystalline reference is shown, not compared
            comparable = row["rho_exp"] and (row["rho_state"] or "amorphous") == "amorphous"
            row[f"rho_dev_pct_{lab}"] = round(100 * (x["density"] - row["rho_exp"]) / row["rho_exp"], 2) if comparable else None
            row[f"sp_dev_{lab}"] = round(x["delta"] - row["sp_exp"], 2) if row["sp_exp"] and x.get("delta") is not None else None
        rows.append(row)
    labs = (la, lb) if B else (la,)
    f = lambda v, fmt: "" if v is None else format(v, fmt)
    head = f"{'polymer':30s} {'rho exp':>7}" + "".join(f" {'rho ' + l:>9} {'dev%':>6}" for l in labs) + f"  {'d exp':>5}" + \
           "".join(f" {'d ' + l:>7} {'dev':>5}" for l in labs)
    print(head)
    for x in rows:
        line = f"{x['polymer'][:30]:30s} {f(x['rho_exp'], '.3f'):>7}"
        for l in labs:
            line += f" {f(x['rho_' + l], '.4f'):>9} {f(x['rho_dev_pct_' + l], '+.1f'):>6}"
        line += f"  {f(x['sp_exp'], '.1f'):>5}"
        for l in labs:
            line += f" {f(x['sp_' + l], '.2f'):>7} {f(x['sp_dev_' + l], '+.1f'):>5}"
        flags = [l for l in labs if (x[f"rho_dev_pct_{l}"] is not None and abs(x[f"rho_dev_pct_{l}"]) > 100 * TOL_RHO)]
        note = "" if (x["rho_state"] or "amorphous") == "amorphous" else f"  (reference {x['rho_state']}: not compared)"
        print(line + ("  density off: " + ",".join(flags) if flags else "") + note)
    for l in labs:
        dr = [abs(x[f"rho_dev_pct_{l}"]) for x in rows if x[f"rho_dev_pct_{l}"] is not None]
        ds = [abs(x[f"sp_dev_{l}"]) for x in rows if x[f"sp_dev_{l}"] is not None]
        if dr:
            print(f"{l}: density mean |dev| {sum(dr) / len(dr):.2f} % over {len(dr)} polymers"
                  + (f"; solubility parameter mean |dev| {sum(ds) / len(ds):.2f} MPa^1/2 over {len(ds)}" if ds else ""))
    if arg("--csv") and rows:
        with open(arg("--csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
