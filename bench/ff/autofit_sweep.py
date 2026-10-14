#!/usr/bin/env python3
"""Every gap bench/ff/react_coverage.py found, fitted by python -m caps.autofit (AIMNet2).

    react_coverage.py --save DIR/products --json DIR/coverage.json      # first: the incomplete products
    autofit_sweep.py DIR/coverage.json --ff oplsaa2024 --out DIR/autofit [--hold-neighbours]

The products are taken in turn, each against the force field plus everything fitted for the ones before, so a term is
fitted once (in the first product that has it). Products with untyped atoms or non-neutral charges are not fitted: they
are listed (typing rules and charges are not what auto-fit does). At the end every product is checked again against the
force field plus DIR/autofit/autofit_all.json, and the terms flagged needs_qm are listed for the cluster.
Run with a Python that has aimnet[ase] (e.g. the research conda environment).
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "data", "python"))
from caps import autofit as A   # noqa: E402

arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
cov = json.load(open(sys.argv[1]))
ff_id = arg("--ff", "oplsaa2024")
base = ff_id if ff_id.endswith(".json") else os.path.join(ROOT, "data", "forcefields", ff_id + ".json")
out = os.path.abspath(arg("--out", "autofit"))
os.makedirs(out, exist_ok=True)
extra = {"bonds": [], "angles": [], "dihedrals": []}
current = base
skipped, runs = [], []
for case in cov["cases"]:
    st = case.get("structure")
    if not st:
        continue
    if case["untyped"]:
        skipped.append((case["case"], f"{case['untyped']} untyped atoms: a typing rule is missing"))
        continue
    miss = A.missing_terms(st, current)
    if not any(miss[k] for k in ("bond", "angle", "dihedral")):
        if miss["charge"]:
            skipped.append((case["case"], "; ".join(miss["charge"])))
        continue
    slug = os.path.splitext(os.path.basename(st))[0]
    work = os.path.join(out, slug)
    cmd = [sys.executable, "-m", "caps.autofit", st, "--ff", current, "--out", work] + (["--hold-neighbours"] if "--hold-neighbours" in sys.argv else [])
    r = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "PYTHONPATH": os.path.join(ROOT, "data", "python")})
    rep = json.load(open(os.path.join(work, "report.json"))) if os.path.exists(os.path.join(work, "report.json")) else {}
    ov = json.load(open(os.path.join(work, "autofit.json"))) if os.path.exists(os.path.join(work, "autofit.json")) else {}
    for k in extra:
        extra[k] += ov.get(k, [])
    current = A.combined_ff(base, extra, os.path.join(out, "ff_cumulative.json"))
    fitted = sum(len(ov.get(k, [])) for k in extra)
    runs.append({"case": case["case"], "structure": st, "fitted": fitted, "complete": rep.get("complete"),
                 "not_fitted": rep.get("not_fitted", []), "needs_qm": rep.get("needs_qm", []), "returncode": r.returncode,
                 "stderr": (r.stderr or "")[-500:] if r.returncode and not rep else ""})
    print(f"{'ok ' if rep.get('complete') else 'GAP'} {case['case']}: {fitted} terms fitted"
          + (f", needs QM {rep['needs_qm']}" if rep.get("needs_qm") else "")
          + (f", not fitted {[n['term'] for n in rep['not_fitted']]}" if rep.get("not_fitted") else "")
          + (f"  [{runs[-1]['stderr'].strip().splitlines()[-1] if runs[-1]['stderr'] else ''}]" if r.returncode and not rep else ""),
          flush=True)

overlay = {"format": "caps-forcefield", "format_version": 1, "name": f"CAPS auto-fit (AIMNet2) for React products, {os.path.basename(base)}",
           "atom_types": [], **{k: v for k, v in extra.items() if v},
           "notes": ["bonded terms the force field lacks in CAPS React products (bench/ff/react_coverage.py), fitted by "
                     "python -m caps.autofit to AIMNet2 relaxed scans; entries marked needs_qm want a QM check"]}
json.dump(overlay, open(os.path.join(out, "autofit_all.json"), "w"), indent=1)
# the check: every saved product against the force field plus the overlay
final = []
for case in cov["cases"]:
    st = case.get("structure")
    if not st:
        continue
    m = A.missing_terms(st, current)
    ok = not any(m[k] for k in ("bond", "angle", "dihedral", "untyped", "charge")) and not case["untyped"]
    final.append((case["case"], ok, m))
n_ok = sum(1 for f in final if f[1])
needs = sorted({r["name"] for r in extra["dihedrals"] if r.get("needs_qm")})
json.dump({"runs": runs, "skipped": skipped, "needs_qm": needs, "terms": {k: len(v) for k, v in extra.items()},
           "after": [{"case": c, "complete": ok, "missing": {k: v for k, v in m.items() if v}} for c, ok, m in final]},
          open(os.path.join(out, "sweep_report.json"), "w"), indent=1)
print(f"\nfitted {sum(len(v) for v in extra.values())} terms ({len(extra['bonds'])} bonds, {len(extra['angles'])} bends, "
      f"{len(extra['dihedrals'])} torsions) -> {os.path.join(out, 'autofit_all.json')}")
print(f"previously incomplete products now complete and neutral: {n_ok} of {len(final)}")
for c, why in skipped:
    print(f"  not auto-fitted: {c}: {why}")
if needs:
    print(f"  needs QM confirmation: {', '.join(needs)}")
