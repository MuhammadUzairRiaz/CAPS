#!/usr/bin/env python3
"""OPLS-AA atom types of the polymer library against an independent typer: every library polymer as a DP-3 oligomer
(hydrogen ends) typed by CAPS with OPLS-AA 2024 (data/typing/oplsaa2024-moltemplate.typing.json) and by foyer's OPLS-AA
(MoSDeF, its own SMARTS rules), compared atom by atom on the OPLS type number.

Two interpreters: CAPS's (this one, with data/python on the path) builds and types; foyer's (--foyer-python, a MoSDeF
environment) types the same mol2 files. A polymer foyer cannot type is reported as such, not counted.

usage: compare_foyer_oplsaa.py OUTDIR [--foyer-python PATH] [--dp 3] [--only P001,P002]
"""
import json, os, re, subprocess, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "data", "python"))
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
OUT = os.path.abspath(sys.argv[1])
FOYER_PY = arg("--foyer-python", "python3")
DP = int(arg("--dp", "3"))
ONLY = set(x for x in arg("--only", "").split(",") if x)
os.makedirs(OUT, exist_ok=True)

FOYER = r'''
import json, sys, parmed, foyer
from foyer.atomtyper import find_atomtypes
ff = foyer.Forcefield(name="oplsaa")
out = {}
for pid, path in json.load(open(sys.argv[1])).items():
    try:
        s = parmed.load_file(path, structure=True)
        tm = find_atomtypes(s, forcefield=ff)
        out[pid] = [tm[i]["atomtype"] for i in range(len(s.atoms))]
    except Exception as e:
        out[pid] = str(e).splitlines()[0][:300]
json.dump(out, open(sys.argv[2], "w"))
'''


def number(t):
    m = re.match(r"^(?:opls_)?(\d+)", t or "")
    return int(m.group(1)) if m else None


def main():
    import caps
    lib = json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))["polymers"]
    caps_types, files, notes = {}, {}, {}
    for p in lib:
        if ONLY and p["id"] not in ONLY:
            continue
        smi = p.get("smiles")
        if not isinstance(smi, str):
            notes[p["id"]] = "not a single repeat unit"
            continue
        try:
            d = caps.polymer(smi, dp=DP)
            r = d.field.assign("oplsaa2024-moltemplate", "forcefield")
            path = os.path.join(OUT, p["id"] + ".mol2")
            d.save(path)
            files[p["id"]] = path
            caps_types[p["id"]] = {"types": [a.get("type", "") for a in r["atoms"]], "el": [a["el"] for a in r["atoms"]],
                                   "desc": [a.get("desc", "") for a in r["atoms"]]}
        except Exception as e:
            notes[p["id"]] = "CAPS: " + str(e)[:200]
    jf, ff = os.path.join(OUT, "files.json"), os.path.join(OUT, "foyer.json")
    json.dump(files, open(jf, "w"))
    script = os.path.join(OUT, "foyer_types.py")
    open(script, "w").write(FOYER)
    subprocess.run([FOYER_PY, script, jf, ff], check=True)
    foy = json.load(open(ff))
    names = {p["id"]: p["name"] for p in lib}
    total = agree = compared = 0
    pairs = Counter()
    rows = []
    for pid in files:
        c = caps_types[pid]
        f = foy.get(pid)
        if not isinstance(f, list):
            rows.append((pid, names[pid], "foyer cannot type it: " + str(f)[:120]))
            continue
        if len(f) != len(c["types"]):
            rows.append((pid, names[pid], "atom counts differ"))
            continue
        untyped = sum(1 for t in c["types"] if not t)
        diff = []
        for i, (ct, ft) in enumerate(zip(c["types"], f)):
            a, b = number(ct), number(ft)
            if a is None:
                continue
            compared += 1
            if a == b:
                agree += 1
            else:
                diff.append(i)
                pairs[(c["el"][i], a, b, c["desc"][i])] += 1
        total += 1
        rows.append((pid, names[pid], "same" if not diff and not untyped else
                     f"{len(diff)} of {len(c['types'])} atoms differ" + (f", {untyped} untyped by CAPS" if untyped else "")))
    with open(os.path.join(OUT, "report.md"), "w") as o:
        o.write(f"# OPLS-AA types: CAPS (OPLS-AA 2024 rules) against foyer, DP-{DP} oligomers\n\n")
        o.write(f"{total} polymers compared, {agree} of {compared} atoms with the same OPLS number ({100 * agree / max(compared, 1):.1f} %)\n\n")
        o.write("| id | polymer | result |\n|---|---|---|\n")
        for r in rows:
            o.write(f"| {r[0]} | {r[1]} | {r[2]} |\n")
        o.write("\n## differing assignments (element, CAPS number, foyer number, CAPS rule) × atoms\n\n")
        for (el, a, b, desc), n in pairs.most_common():
            o.write(f"- {el} CAPS {a} ({desc}) · foyer {b} · {n} atoms\n")
        for pid, n in notes.items():
            o.write(f"- {pid}: {n}\n")
    print(open(os.path.join(OUT, "report.md")).read())


if __name__ == "__main__":
    main()
