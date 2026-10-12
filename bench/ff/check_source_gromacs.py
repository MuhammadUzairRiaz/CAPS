#!/usr/bin/env python3
"""CAPS's parameter assignment against the force field's own program: GROMACS with its own force-field directory.

For each structure CAPS types it and writes a GROMACS topology with every parameter it assigned (caps ff apply --gromacs).
A second topology keeps CAPS's atom types, charges and term lists but no parameters: it includes the force field's own
directory (e.g. oplsaa.ff/forcefield.itp) and GROMACS looks every bond, angle, torsion, 1-4 pair and Lennard-Jones
parameter up itself (impropers keep CAPS's, as GROMACS's own come from residue topologies). Both run in GROMACS
(grompp, mdrun -rerun) and every energy term is compared: a difference is a parameter CAPS assigned differently from
the force field's own program, or a term GROMACS has no parameter for (grompp names it).

usage: check_source_gromacs.py FF_ID GMX_FF_DIR STRUCTURE|SMILES... [--charges types] [--keep DIR]
       e.g. check_source_gromacs.py oplsaa-gromacs oplsaa.ff "CCOC(=O)C" peptide:GHWRK
Needs GMX (default gmx on PATH); GMX_FF_DIR is looked up in GMXLIB / GROMACS's share/top.
"""
import os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
GMX = os.environ.get("GMX", shutil.which("gmx") or "gmx")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
fid, ffdir = sys.argv[1], sys.argv[2]
items = [a for i, a in enumerate(sys.argv[3:], 3) if not a.startswith("--") and sys.argv[i - 1] not in ("--charges", "--keep")]
charges = arg("--charges", "types")
work = arg("--keep", "") or tempfile.mkdtemp()
os.makedirs(work, exist_ok=True)
TERMS = ["Bond", "Angle", "Proper Dih.", "Ryckaert-Bell.", "Improper Dih.", "Per. Imp. Dih.", "LJ-14", "Coulomb-14", "LJ (SR)", "Coulomb (SR)"]


def energies(d, stem):
    run = lambda a, inp=None: subprocess.run([GMX, "-quiet"] + a, cwd=d, input=inp, capture_output=True, text=True)
    r = run(["grompp", "-f", "case.mdp", "-c", "case.gro", "-p", stem + ".top", "-o", stem + ".tpr", "-maxwarn", "10"])
    if r.returncode:
        return None, " ".join(l.strip() for l in (r.stdout + r.stderr).splitlines() if "No default" in l or "ERROR" in l)[:300] or "grompp failed"
    run(["mdrun", "-s", stem + ".tpr", "-rerun", "case.gro", "-deffnm", stem, "-nt", "1"])
    names = run(["energy", "-f", stem + ".edr", "-o", stem + "_e.xvg"], "0\n").stdout + run(["energy", "-f", stem + ".edr", "-o", stem + "_e.xvg"], "0\n").stderr
    pick = [t for t in TERMS if re.search(r"\b\d+\s+" + re.escape(t.replace(" ", "-")) + r"(\s|$)", names)]   # gmx energy writes spaces as -
    run(["energy", "-f", stem + ".edr", "-o", stem + "_e.xvg"], "\n".join(t.replace(" ", "-") for t in pick) + "\n\n")
    vals = [l.split() for l in open(os.path.join(d, stem + "_e.xvg")) if l.strip() and l[0] not in "#@"]
    return dict(zip(pick, map(float, vals[0][1:]))), ""


def native_top(d):
    """case.top with the force field's own directory and no parameters on bonds, pairs, angles and proper torsions"""
    itp = open(os.path.join(d, "case.itp")).read().splitlines()
    out, sec, seen = [], None, set()
    for l in itp:
        m = re.match(r"\s*\[\s*(\w+)\s*\]", l)
        if m:
            sec = m.group(1)
            out.append(l)
            continue
        w = l.split(";")[0].split()
        if not w:
            out.append(l)
            continue
        if sec == "bonds":
            out.append(f"{w[0]} {w[1]} 1")
        elif sec == "pairs":
            out.append(f"{w[0]} {w[1]} 1")
        elif sec == "angles":
            out.append(f"{w[0]} {w[1]} {w[2]} 1")
        elif sec == "dihedrals" and w[4] in ("9", "3", "1", "5"):   # a proper torsion: once, the force field's own form
            if tuple(w[:4]) not in seen:
                seen.add(tuple(w[:4]))
                out.append(f"{w[0]} {w[1]} {w[2]} {w[3]} 3")
        else:
            out.append(l)
    open(os.path.join(d, "native.itp"), "w").write("\n".join(out) + "\n")
    top = open(os.path.join(d, "case.top")).read()
    body = top[top.index("[ system ]"):] if "[ system ]" in top else ""
    head = f'#include "{ffdir}/forcefield.itp"\n#include "native.itp"\n\n'
    open(os.path.join(d, "native.top"), "w").write(head + body)


fails = 0
for k, it in enumerate(items):
    d = os.path.join(work, f"case{k}")
    os.makedirs(d, exist_ok=True)
    if it.startswith("peptide:"):
        s = os.path.join(d, "in.pdb")
        subprocess.run([CAPS, "peptide", it.split(":", 1)[1], "-o", s], capture_output=True, check=True)
    elif os.path.exists(it):
        s = it
    else:
        s = os.path.join(d, "in.mol2")
        subprocess.run([CAPS, "build", it, "--ff", "uff", "-o", s], capture_output=True, check=True)
    r = subprocess.run([CAPS, "ff", "apply", s, "--ff", os.path.join(ROOT, "data", "forcefields", fid + ".json"), "--charges", charges,
                        "--gromacs", os.path.join(d, "case")], capture_output=True, text=True)
    if r.returncode or "missing parameters" in r.stdout:
        print(f"{it}: CAPS failed: {(r.stderr or r.stdout).strip().splitlines()[-1][:200]}")
        fails += 1
        continue
    if os.path.exists(os.path.join(d, "case.itp")) and "#include" in open(os.path.join(d, "case.top")).read():
        pass
    native_top(d)
    ec, why_c = energies(d, "case")
    en, why_n = energies(d, "native")
    if ec is None or en is None:
        print(f"{it}: {'CAPS topology: ' + why_c if ec is None else 'GROMACS with ' + ffdir + ': ' + why_n}")
        fails += 1
        continue
    # torsions: CAPS writes multi-term periodic (Proper Dih.), the force field Ryckaert-Bellemans: compared as one sum
    tors = lambda e: e.get("Proper Dih.", 0.0) + e.get("Ryckaert-Bell.", 0.0)
    cmp = {"bonds": (ec.get("Bond", 0), en.get("Bond", 0)), "angles": (ec.get("Angle", 0), en.get("Angle", 0)),
           "torsions": (tors(ec), tors(en)), "LJ-14": (ec.get("LJ-14", 0), en.get("LJ-14", 0)),
           "Coulomb-14": (ec.get("Coulomb-14", 0), en.get("Coulomb-14", 0)), "LJ": (ec.get("LJ (SR)", 0), en.get("LJ (SR)", 0)),
           "Coulomb": (ec.get("Coulomb (SR)", 0), en.get("Coulomb (SR)", 0))}
    worst = max(abs(a - b) / max(1.0, abs(b)) for a, b in cmp.values())
    ok = worst < 1e-4
    fails += 0 if ok else 1
    print(f"{it}: {'ok' if ok else 'DIFFERS'} · largest difference {worst:.1e} (relative) · " +
          " ".join(f"{n} {a:.4f}/{b:.4f}" for n, (a, b) in cmp.items()) + " kJ/mol (CAPS / " + ffdir + ")")
print(f"{len(items) - fails} of {len(items)} agree with GROMACS's own {ffdir} (work: {work})")
sys.exit(1 if fails else 0)
