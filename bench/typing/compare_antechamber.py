#!/usr/bin/env python3
"""GAFF and GAFF2 atom types from CAPS against AmberTools' antechamber, molecule by molecule and atom by atom.

  compare_antechamber.py prepare SET.smi OUTDIR [--caps build/cli/caps]
      builds each molecule with `caps build` (3D, explicit hydrogens), writes OUTDIR/NAME.mol2 with element symbols in the
      atom-type column and CAPS's bond orders (what antechamber reads), and CAPS's own types: NAME.gaff.types (GAFF,
      gaff-amber16) and NAME.gaff2.types (GAFF2, gaff-amber25). Also OUTDIR/run_antechamber.sh.
  bash OUTDIR/run_antechamber.sh          (where AmberTools is: antechamber -at gaff / gaff2, atom order kept)
  compare_antechamber.py compare OUTDIR [--show 50]
      reads NAME.ac_gaff.mol2 and NAME.ac_gaff2.mol2 and lists every atom whose types differ.
"""
import argparse, os, re, subprocess, sys

TAGS = (("gaff", "gaff-amber16"), ("gaff2", "gaff-amber25"))
REF = {"gaff": "gaff", "gaff2": "gaff2"}   # antechamber's -at for each CAPS entry
SYM = {"CL": "Cl", "BR": "Br", "SI": "Si", "NA": "Na", "MG": "Mg"}


def element(name):
    m = re.match(r"([A-Za-z]+)", name)
    s = m.group(1)
    for two in ("Cl", "Br", "Si", "Na", "Mg", "Li", "Al", "Zn", "Ca"):
        if s.lower().startswith(two.lower()):
            return two
    return s[0].upper()


def rewrite_mol2(src, dst):
    """the same mol2, element symbols as atom types and no charges column games"""
    out, sec = [], None
    for line in open(src):
        if line.startswith("@<TRIPOS>"):
            sec = line.strip()
            out.append(line)
            continue
        if sec == "@<TRIPOS>ATOM" and line.strip():
            w = line.split()
            el = element(w[1])
            out.append("%7s %-6s %10s %10s %10s %-6s %3s %-4s %10s\n" % (w[0], w[1], w[2], w[3], w[4], el, w[6], w[7], "0.000000"))
        else:
            out.append(line)
    open(dst, "w").write("".join(out))


def types_of_mol2(path):
    t, sec = [], None
    for line in open(path):
        if line.startswith("@<TRIPOS>"):
            sec = line.strip()
            continue
        if sec == "@<TRIPOS>ATOM" and line.strip():
            t.append(line.split()[5])
    return t


def types_file(path):
    return [w.split()[-1] for w in open(path) if w.strip() and not w.startswith("#")]


def prepare(a):
    os.makedirs(a.outdir, exist_ok=True)
    names = []
    for line in open(a.smi):
        if not line.strip() or line.startswith("#"):
            continue
        smi, name = line.split()[:2]
        raw = os.path.join(a.outdir, name + ".caps.mol2")
        r = subprocess.run([a.caps, "build", smi, "-o", raw, "--seed", "1"], capture_output=True, text=True)
        if r.returncode or not os.path.exists(raw):
            print("build failed:", name, smi, r.stderr.strip()[-200:])
            continue
        rewrite_mol2(raw, os.path.join(a.outdir, name + ".mol2"))
        ok = True
        for tag, ff in TAGS:
            r = subprocess.run([a.caps, "ff", "type", raw, "--ff", ff, "-o", os.path.join(a.outdir, f"{name}.{tag}.types")],
                               capture_output=True, text=True)
            if r.returncode:
                print(f"CAPS {tag} typing failed:", name, (r.stdout + r.stderr).strip()[-300:])
                ok = False
        names.append(name)
    with open(os.path.join(a.outdir, "run_antechamber.sh"), "w") as f:
        f.write("#!/bin/bash\n# antechamber's GAFF and GAFF2 types for every molecule here (atom order kept, no charges)\n")
        f.write('cd "$(dirname "$0")"\nantechamber -L > /dev/null 2>&1; antechamber 2>&1 | head -3 > antechamber_version.txt\n')
        for n in names:
            for at in ("gaff", "gaff2"):
                f.write(f"antechamber -i {n}.mol2 -fi mol2 -o {n}.ac_{at}.mol2 -fo mol2 -at {at} -dr no -pf y > {n}.ac_{at}.log 2>&1\n")
        f.write("ls *.ac_gaff*.mol2 | wc -l\n")
    print(len(names), "molecules in", a.outdir)


# differences that are not typing errors, each with its reason (listed apart from the count of disagreements)
EXPLAINED = [
    (lambda tag, mol, x, y: x == "Si" and y in ("?", "(untyped)"), "GAFF / GAFF2 have no Si type: untyped by both"),
    (lambda tag, mol, x, y: mol == "cyclooctatetraene" and {x, y} == {"cc", "cd"},
     "antechamber re-assigns the Kekule structure; both alternate cc/cd by their own bond orders (CAPS keeps the file's)"),
    (lambda tag, mol, x, y: tag == "gaff" and x in ("op", "oq") and y == "os", "AmberTools 16 GAFF parameters have no op / oq"),
]


def compare(a):
    tot = {t: [0, 0, 0] for t, _ in TAGS}
    explained = {}
    rows = []
    for f in sorted(os.listdir(a.outdir)):
        if not f.endswith(".mol2") or ".caps." in f or ".ac_" in f:
            continue
        name = f[:-5]
        for tag, _ in TAGS:
            ref = os.path.join(a.outdir, f"{name}.ac_{REF[tag]}.mol2")
            mine = os.path.join(a.outdir, f"{name}.{tag}.types")
            if not (os.path.exists(ref) and os.path.exists(mine)):
                rows.append((tag, name, "missing", "", ""))
                continue
            r, c = types_of_mol2(ref), types_file(mine)
            if len(r) != len(c):
                rows.append((tag, name, "atom count", len(r), len(c)))
                continue
            names = [l.split()[1] for l in open(os.path.join(a.outdir, f)) if re.match(r"^\s*\d+\s+\S+\s+-?\d", l)]
            tot[tag][1] += len(r)
            for i, (x, y) in enumerate(zip(r, c)):
                why = next((w for f_, w in EXPLAINED if f_(tag, name, x, y)), None)
                if x.lower() == y.lower():
                    tot[tag][0] += 1
                elif why:
                    tot[tag][2] += 1
                    explained.setdefault((tag, why), set()).add(name)
                else:
                    rows.append((tag, name, i + 1, x, y))
    for tag, (ok, n, ex) in tot.items():
        print(f"{tag.upper()} ({dict(TAGS)[tag]}): {ok}/{n} atoms as antechamber ({100 * ok / max(n, 1):.1f} %); "
              f"{ex} more differ for a stated reason; {n - ok - ex} unexplained")
    for (tag, why), mols in sorted(explained.items()):
        print(f"  explained, {tag}: {why} ({', '.join(sorted(mols))})")
    shown = 0
    for r in rows:
        if shown >= a.show:
            break
        print(f"  {r[0]:7s} {r[1]:28s} atom {str(r[2]):4s} antechamber {r[3]!s:5s} CAPS {r[4]}")
        shown += 1


ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
p = sub.add_parser("prepare"); p.add_argument("smi"); p.add_argument("outdir"); p.add_argument("--caps", default="build/cli/caps")
p = sub.add_parser("compare"); p.add_argument("outdir"); p.add_argument("--show", type=int, default=50)
a = ap.parse_args()
prepare(a) if a.cmd == "prepare" else compare(a)
