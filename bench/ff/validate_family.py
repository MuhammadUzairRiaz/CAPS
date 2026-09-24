#!/usr/bin/env python3
"""Validate a converted DL_FIELD library on its own molecule templates: for each template, build a PDB
(template_pdb.py), let DL_FIELD write its DL_POLY FIELD (run_dlfield.py --dlpoly), and compare CAPS's assignment term
by term (compare_field.py).

usage: validate_family.py DLFIELD_KEYWORD LIBRARY FF.json [--max N] [--only a,b] [--work DIR]
"""
import os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.environ.get("DLFIELD", "~/project/dl_f_4.13") + "/lib"
kw, lib, ffjson = sys.argv[1], sys.argv[2], sys.argv[3]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
nmax = int(arg("--max", "12"))
work = arg("--work", tempfile.mkdtemp())
os.makedirs(work, exist_ok=True)
sf = os.path.join(LIB, lib + ".sf")
builder = "ionic_pdb.py" if "--ionic" in sys.argv else "template_pdb.py"   # --ionic: lattices of ionic formula units
names = arg("--only", "")
names = names.split(",") if names else subprocess.run([sys.executable, os.path.join(HERE, builder), sf, "--list"],
                                                      capture_output=True, text=True).stdout.split()
passed, failed, skipped = [], [], []
for name in names:
    if len(passed) + len(failed) >= nmax: break
    pdb = os.path.join(work, name + ".pdb")
    out = os.path.join(work, name)
    r = subprocess.run([sys.executable, os.path.join(HERE, builder), sf, name, pdb], capture_output=True, text=True)
    if r.returncode: skipped.append((name, "template: " + r.stderr.strip()[-80:])); continue
    r = subprocess.run([sys.executable, os.path.join(HERE, "run_dlfield.py"), kw, pdb, out, "--dlpoly"], capture_output=True, text=True)
    if r.returncode or not os.path.exists(os.path.join(out, "dl_poly.FIELD")):
        skipped.append((name, "DL_FIELD: " + (r.stdout.strip().splitlines() or ["?"])[-1][:100])); continue
    # the template's DL_FIELD atom types → the library's keys
    sys.path.insert(0, HERE)
    import template_pdb
    types, keys, mols = template_pdb.parse_sf(sf)
    key_of = dict(template_pdb.ATOM_KEY)
    # DL_FIELD re-types some atoms from their surroundings ("Atom 2 - HC_alkane changed to HCp_alkane"): apply that,
    # then the ATOM_TYPE table gives the key the parameters are looked up with
    import re as _re
    dltypes = [t for _, t in mols[name]["atoms"]]
    for line in open(os.path.join(out, "dl_field.output"), errors="replace"):
        m = _re.match(r"\s*Atom\s+(\d+)\s+-\s+(\S+)\s+changed to\s+(\S+)", line)
        if m and 1 <= int(m.group(1)) <= len(dltypes): dltypes[int(m.group(1)) - 1] = m.group(3)
    tf = os.path.join(out, "keys.txt")
    open(tf, "w").write("\n".join(key_of.get(t, t) for t in dltypes) + "\n")
    extra = ["--from-field"] if builder == "ionic_pdb.py" else ["--types", tf]   # ionic: compare per-type terms on DL_FIELD's own atom list
    c = subprocess.run([sys.executable, os.path.join(HERE, "compare_field.py"), out, pdb, ffjson] + extra, capture_output=True, text=True)
    last = (c.stdout.strip().splitlines() or ["?"])[-1]
    if c.returncode == 0: passed.append(name)
    else:
        failed.append(name)
        print(f"--- {name}: {last}")
        for l in c.stdout.splitlines():
            if "MISMATCH" in l or "EXTRA" in l or "Trace" in l or "caps ff" in l: print("   " + l.strip()[:220])
print(f"\n{kw}: {len(passed)} passed, {len(failed)} failed, {len(skipped)} skipped (DL_FIELD could not build them)")
for n, why in skipped[:8]: print(f"   skipped {n}: {why}")
print("passed:", " ".join(passed))
