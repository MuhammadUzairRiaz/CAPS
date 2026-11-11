#!/usr/bin/env python3
"""PCFF atom types from CAPS against pysimm's own PCFF typer (pysimm 1.1.0, github.com/polysimtools/pysimm, MIT), an
open, independent implementation. The molecules are the 122 of bench/typing/gaff_reference_set.smi: pysimm reads
NAME.ac_gaff.mol2 (antechamber's bond orders), CAPS types NAME.caps.mol2 (bare connectivity), same atom order.
pysimm uses a reduced PCFF type set (generic c for every sp3 carbon, h for every C-H), so two scores are given:
  types       the same type name
  parameters  the same parameters: both types map to the same nonbond and bond classes through pcff.frc's own
              equivalence table (c3 -> c, hc -> h, ...), so the two programs would build the same force field there
usage: compare_pysimm_pcff.py MOLDIR [--caps build/cli/caps] [--show 40]
"""
import argparse, glob, json, os, subprocess, tempfile, warnings, contextlib, io
warnings.filterwarnings("ignore")
from pysimm import system, forcefield

ap = argparse.ArgumentParser(); ap.add_argument("moldir"); ap.add_argument("--caps", default="build/cli/caps")
ap.add_argument("--show", type=int, default=40); a = ap.parse_args()
lib = json.load(open("data/forcefields/pcff-frc.json"))
eq = {t["name"]: (t.get("equivalence", {}).get("vdw", t["name"]), t.get("equivalence", {}).get("bond", t["name"]))
      for t in lib["atom_types"]} if "atom_types" in lib else {}
if not eq:   # the library keeps its types under another key
    for k, v in lib.items():
        if isinstance(v, list) and v and isinstance(v[0], dict) and "equivalence" in v[0] and "element" in v[0]:
            eq = {t["name"]: (t["equivalence"].get("vdw", t["name"]), t["equivalence"].get("bond", t["name"])) for t in v}
ff = forcefield.Pcff(); tmp = tempfile.mkdtemp()
n_mol = n_ok = at = same = same_par = full = full_par = 0; diffs = {}
for f in sorted(glob.glob(os.path.join(a.moldir, "*.ac_gaff.mol2"))):
    name = os.path.basename(f).split(".")[0]; n_mol += 1
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            s = system.read_mol2(f); ff.assign_ptypes(s)
        ps = [p.type.name if p.type else "?" for p in s.particles]
    except Exception as e:
        diffs.setdefault(("pysimm failed", type(e).__name__), []).append(name); continue
    out = os.path.join(tmp, name + ".types")
    subprocess.run([a.caps, "ff", "type", os.path.join(a.moldir, name + ".caps.mol2"), "--ff", "pcff-frc", "-o", out], capture_output=True)
    if not os.path.exists(out):
        diffs.setdefault(("caps failed", ""), []).append(name); continue
    cs = [w.split()[-1] for w in open(out) if w.strip() and not w.startswith("#")]
    if len(cs) != len(ps):
        diffs.setdefault(("atom count differs", ""), []).append(name); continue
    n_ok += 1; m = mp = 0
    for p, c in zip(ps, cs):
        at += 1
        if p == c:
            m += 1; mp += 1
        elif p in eq and c in eq and eq[p] == eq[c]:
            mp += 1
        else:
            diffs.setdefault((p, c), []).append(name)
    same += m; same_par += mp; full += m == len(ps); full_par += mp == len(ps)
print(f"molecules {n_mol}, typed by both {n_ok}; atoms {at}")
print(f"same type        {same}/{at} ({100*same/at:.1f} %), molecules {full}/{n_ok}")
print(f"same parameters  {same_par}/{at} ({100*same_par/at:.1f} %), molecules {full_par}/{n_ok}")
for k, v in sorted(diffs.items(), key=lambda kv: -len(kv[1]))[:a.show]:
    print(f"  pysimm {k[0]:>14}  CAPS {k[1]:<8} {len(v):4d}  {', '.join(sorted(set(v))[:6])}")
