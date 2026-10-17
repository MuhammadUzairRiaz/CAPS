#!/usr/bin/env python3
"""Rubber validation bench (rubber request F7): the diene rubbers, their sulfur and peroxide networks and the other
common rubbers, each built and typed with OPLS-UA 2024 and the OPLS-AA files, checked for untyped atoms, missing terms
and neutrality; every complete one is exported to LAMMPS and run for 0 steps there (the file must read and evaluate).

    rubber_bench.py [--ff oplsua2024,oplsaa2024,...] [--only NAME] [--lmp PATH] [--json out.json] [--keep DIR]

Networks: 4 chains of 10 repeat units at 0.6 g/cm3, 6 sulfur donors (H-Sn-H), sulfur_allylic between chains to 4
links, unreacted donors removed (react(remove_unreacted=True)); the peroxide network has 3 C-C links. The table lists
per system and force field: untyped, missing, net charge, the sulfur bridges by rank, and the LAMMPS check; then every
CAPS-added atom type the complete systems use, with its source, for a methods section.
"""
import json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
os.environ.setdefault("CAPS_LIB", os.path.join(ROOT, "build", "capi", "libcaps.dylib"))
sys.path.insert(0, os.path.join(ROOT, "data", "python"))
import caps  # noqa: E402

arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FFS = (arg("--ff") or "oplsua2024,oplsaa2024").split(",")
ONLY = arg("--only", "")
LMP = arg("--lmp") or shutil.which("lmp") or os.path.expanduser("~/.local/bin/lmp")
work = arg("--keep") or tempfile.mkdtemp(prefix="rubber_bench_")

PB, PI = "[*]C/C=C\\C[*]", "[*]C/C(C)=C\\C[*]"
CHAINS = {"cis-PB": PB, "trans-PB": "[*]C/C=C/C[*]", "cis-PI": PI, "1,2-PB": "[*]CC([*])C=C", "IIR": "[*]CC([*])(C)C",
          "EPDM (ethylene unit)": "[*]CC[*]", "CR": "[*]C/C(Cl)=C\\C[*]", "NBR (acrylonitrile unit)": "[*]CC([*])C#N",
          "SBR (styrene unit)": "[*]CC([*])c1ccccc1", "PDMS": "[*]O[Si](C)(C)[*]"}
NETWORKS = {f"{p} + S{n}": (u, "S" * n, "sulfur_allylic", 4) for p, u in (("PB", PB), ("PI", PI)) for n in (1, 2, 3, 8)}
NETWORKS["PB peroxide"] = (PB, None, "peroxide_allylic", 3)


def make(name):
    if name in CHAINS:
        return caps.polymer(CHAINS[name], dp=10, chains=4, density=0.6, seed=3)
    u, donor, tmpl, links = NETWORKS[name]
    d = caps.polymer(u, dp=10, chains=4, density=0.6, seed=3)
    if donor:
        d.insert(donor, 6)
    d.react(tmpl, cycles=60, per_cycle=2, between_chains=True, auto_capture=True, capture_max=9, relax=False,
            crosslinks=("crosslinks", links), remove_unreacted=bool(donor))
    return d


def lammps_check(d, name, ff):
    folder = os.path.join(work, "".join(c if c.isalnum() else "_" for c in f"{name}_{ff}"))
    try:
        r = d.export_engines(folder, stem="sys", lammps=True, gromacs=False, run="check")
    except Exception as e:
        return f"export refused: {str(e)[:80]}", ""
    inp = os.path.join(folder, "sys.in")
    if not os.path.exists(inp) or not os.path.exists(LMP or ""):
        return "exported (no LAMMPS to run)", ""
    style = next((l.split()[1] for l in open(inp) if l.startswith("pair_style")), "")
    special = next((" ".join(l.split()[1:]) for l in open(inp) if l.startswith("special_bonds")), "")
    p = subprocess.run([LMP, "-in", "sys.in", "-log", "none", "-echo", "none"], cwd=folder, capture_output=True, text=True, timeout=600)
    ok = p.returncode == 0 and "ERROR" not in p.stdout
    return ("LAMMPS run 0 ok" if ok else "LAMMPS failed: " + next((l for l in p.stdout.splitlines() if "ERROR" in l), p.stderr.strip()[-80:])), f"{style}; special_bonds {special}"


def main():
    rows, caps_types = [], {}
    for name in list(CHAINS) + list(NETWORKS):
        if ONLY and ONLY not in name:
            continue
        for ff in FFS:
            try:
                d = make(name)   # a fresh document per force field (a united-atom assignment folds the hydrogens)
                sb = d.sulfur_bridges() if name in NETWORKS else {}   # topology, before the force field folds anything
                r = d.field.assign(ff)
            except Exception as e:
                rows.append({"system": name, "ff": ff, "error": str(e)[:120]})
                print(f"ERR  {name:24s} {ff:24s} {str(e)[:100]}")
                continue
            complete = r["complete"] and r.get("untyped", 0) == 0 and abs(r.get("net_charge", 0)) < 1e-6
            lmp, styles = lammps_check(d, name, ff) if complete else ("", "")
            row = {"system": name, "ff": ff, "untyped": r.get("untyped", 0), "missing": r["missing"],
                   "net_charge": round(r.get("net_charge", 0), 6), "complete": complete, "out_of_domain": r.get("out_of_domain", ""),
                   "bridges": sb.get("bridges", {}), "removed_donors": getattr(d, "removed_unreacted", {}).get("molecules", 0),
                   "lammps": lmp, "styles": styles}
            rows.append(row)
            if complete:
                for a in r.get("atoms", []):
                    t = a.get("type", "")
                    if t.startswith("caps_"):
                        caps_types.setdefault((ff, t), a.get("desc", ""))
            note = (f"bridges {row['bridges']} " if row["bridges"] else "") + (lmp or (row["out_of_domain"][:90] if row["out_of_domain"] else f"missing {r['missing'][:2]}"))
            print(f"{'ok ' if complete else 'GAP'}  {name:24s} {ff:24s} untyped {row['untyped']:3d} missing {len(row['missing']):2d} "
                  f"net {row['net_charge']:+.4f}  {note}")
    if caps_types:
        print("\nCAPS-added atom types in the complete systems (each with its source in the force-field file):")
        lib = {}
        for ff, t in sorted(caps_types):
            if ff not in lib:
                path = os.path.join(ROOT, "data", "forcefields", ff + ".json")
                lib[ff] = {x["name"]: x for x in json.load(open(path))["atom_types"]} if os.path.exists(path) else {}
            desc = lib[ff].get(t, {}).get("description", caps_types[(ff, t)])
            print(f"  {ff}: {t}: {desc[:220]}")
    n_ok = sum(1 for r in rows if r.get("complete"))
    print(f"\n{n_ok} of {len(rows)} system / force-field pairs complete and neutral (files in {work})")
    if arg("--json"):
        json.dump({"rows": rows, "caps_types": [{"ff": f, "type": t} for f, t in sorted(caps_types)]}, open(arg("--json"), "w"), indent=1)


if __name__ == "__main__":
    main()
