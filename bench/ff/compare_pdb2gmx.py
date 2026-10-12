#!/usr/bin/env python3
"""CAPS's building-block topology against pdb2gmx's, term by term (GROMOS 54A7).

  compare_pdb2gmx.py STRUCTURE FF_DIR CAPS_FF.json [--keep DIR]

pdb2gmx (-ignh, -water none, the default termini) builds the reference topology and its coordinates; CAPS applies its
building blocks to pdb2gmx's own .gro and writes its GROMACS topology. Both are read with the #define'd codes resolved:
atom types and charges, bonds, angles, proper dihedrals, impropers (each with its parameters) and 1-4 pairs must be the
same sets; then grompp + mdrun -rerun on both give each energy term (vacuum, the same .mdp).
"""
import math, os, re, subprocess, sys, tempfile

struct, ffdir, capsff = sys.argv[1:4]
work = sys.argv[sys.argv.index("--keep") + 1] if "--keep" in sys.argv else tempfile.mkdtemp()
os.makedirs(work, exist_ok=True)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
GMX = os.environ.get("GMX", "gmx")
ffname = os.path.basename(ffdir.rstrip("/")).replace(".ff", "")
gmxtop = ffdir if os.path.isdir(ffdir) else os.path.join(os.path.dirname(subprocess.run(["which", GMX], capture_output=True, text=True).stdout.strip()), "..", "share", "gromacs", "top", ffdir)

defines = {}
for l in open(os.path.join(ffdir, "ffbonded.itp")):
    m = re.match(r"#define\s+(\S+)\s+(.*)", l)
    if m:
        defines[m.group(1)] = m.group(2).split(";")[0].split()


def read_top(path):
    """{section: [rows]} of the molecule, #include'd .itp files of the work dir followed, macros resolved"""
    out, sec = {}, None
    lines = []
    for l in open(path):
        m = re.match(r'\s*#include\s+"([^"]+)"', l)
        if m and os.path.exists(os.path.join(os.path.dirname(path), m.group(1))):
            lines += open(os.path.join(os.path.dirname(path), m.group(1))).read().splitlines()
        else:
            lines.append(l.rstrip("\n"))
    for raw in lines:
        l = raw.split(";")[0].strip()
        if not l or l.startswith("#"):
            continue
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        w = l.split()
        w = sum(([*defines[x]] if x in defines else [x] for x in w), [])
        out.setdefault(sec, []).append(w)
    return out


def terms(top, sec, n):
    got = {}
    for w in top.get(sec, []):
        at = tuple(int(x) for x in w[:n])
        func = int(w[n])
        func = 1 if (sec == "dihedrals" and func == 9) else func   # 9: function 1 with several terms on one quadruple
        par = tuple(round(float(x), 4) for x in w[n + 1:])
        key = at if at <= at[::-1] else at[::-1]   # an improper reversed (l k j i) is the same angle
        got.setdefault((key, func), []).append(par)
    return {k: sorted(v) for k, v in got.items()}


ref_gro = os.path.join(work, "ref.gro")
p = subprocess.run([GMX, "-quiet", "pdb2gmx", "-f", struct, "-o", ref_gro, "-p", os.path.join(work, "ref.top"), "-ff", ffname,
                    "-water", "none", "-ignh"], capture_output=True, text=True, cwd=work)
if p.returncode:
    raise SystemExit("pdb2gmx failed: " + p.stderr[-400:])
# the coordinates in a box 3 nm wider than the molecule (pdb2gmx keeps the input's cell, which may be smaller than the molecule)
subprocess.run([GMX, "-quiet", "editconf", "-f", ref_gro, "-o", os.path.join(work, "box.gro"), "-d", "3"], capture_output=True, cwd=work)
# CAPS reads pdb2gmx's atoms with pdb2gmx's bonds (a mol2): pdb2gmx places some hydrogens too close to other atoms
# (a proline N-terminus's H1 0.8 A from CD) for bonds perceived from distances
R0 = read_top(os.path.join(work, "ref.top"))
gl = open(os.path.join(work, "box.gro")).read().splitlines()
na = int(gl[1])
EL = {"C": "C", "N": "N", "O": "O", "H": "H", "S": "S", "P": "P"}
with open(os.path.join(work, "ref.mol2"), "w") as f:
    f.write(f"@<TRIPOS>MOLECULE\nref\n{na} {len(R0['bonds'])} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
    for i, l in enumerate(gl[2:2 + na]):
        x, y, z = (float(v) * 10 for v in l[20:].split()[:3])
        name, res, rid = l[10:15].strip(), l[5:10].strip(), int(l[:5])
        f.write(f"{i + 1} {name} {x:.4f} {y:.4f} {z:.4f} {EL.get(name[0], name[0])} {rid} {res} 0.0\n")
    f.write("@<TRIPOS>BOND\n")
    for k, w in enumerate(R0["bonds"]):
        f.write(f"{k + 1} {w[0]} {w[1]} 1\n")
    f.write("@<TRIPOS>CRYSIN\n" + " ".join(f"{float(v) * 10:.4f}" for v in gl[2 + na].split()[:3]) + " 90 90 90 1 1\n")
c = subprocess.run([CAPS, "ff", "apply", os.path.join(work, "ref.mol2"), "--ff", capsff, "--gromacs", os.path.join(work, "caps")], capture_output=True, text=True)
if c.returncode:
    raise SystemExit("CAPS failed: " + (c.stderr or c.stdout)[-400:])
R, C = read_top(os.path.join(work, "ref.top")), read_top(os.path.join(work, "caps.top"))
bad = 0
ra = [(w[1], round(float(w[6]), 4)) for w in R["atoms"]]
ca = [(w[1], round(float(w[6]), 4)) for w in C["atoms"]]
diff = [i + 1 for i, (x, y) in enumerate(zip(ra, ca)) if x[0].upper() != y[0].upper() or abs(x[1] - y[1]) > 1e-4]
print(f"atoms: {len(ra)} / {len(ca)}, {len(diff)} differ in type or charge" + (f" (first {diff[:5]})" if diff else ""))
bad += len(diff) + (len(ra) != len(ca))
# pdb2gmx writes impropers as [ dihedrals ] function 2; split them out of both
for T in (R, C):
    T["impropers"] = [w for w in T.get("dihedrals", []) if w[4] == "2"]
    T["dihedrals"] = [w for w in T.get("dihedrals", []) if w[4] != "2"]
for sec, n in (("bonds", 2), ("angles", 3), ("dihedrals", 4), ("impropers", 4), ("pairs", 2)):
    r, c_ = terms(R, sec, n), terms(C, sec, n)
    if sec == "pairs":   # pdb2gmx writes function 1 with the pairtypes; CAPS the pair's own sigma and epsilon
        r = {k[0]: 1 for k in r}
        c_ = {k[0]: 1 for k in c_}
    only_r = [k for k in r if k not in c_]
    only_c = [k for k in c_ if k not in r]
    if sec in ("dihedrals",):   # several multiplicities on one quadruple: compare the sets of (k, n, phi)
        pass
    pdiff = [k for k in r if k in c_ and sec != "pairs" and not all(
        len(a) == len(b) and all(math.isclose(x, y, rel_tol=2e-4, abs_tol=2e-4) for x, y in zip(a, b)) for a, b in zip(r[k], c_[k]))]
    print(f"{sec}: pdb2gmx {len(r)}, CAPS {len(c_)}; only pdb2gmx {len(only_r)}, only CAPS {len(only_c)}, parameters differ {len(pdiff)}")
    for k in (only_r[:3] + only_c[:3]):
        print("    ", "pdb2gmx" if k in r else "CAPS", k, r.get(k) or c_.get(k))
    for k in pdiff[:3]:
        print("     differ", k, r[k], c_[k])
    bad += len(only_r) + len(only_c) + len(pdiff)
# energies: both topologies on pdb2gmx's coordinates, CAPS's .mdp (plain cut-off in a large box)
MDP = os.path.join(work, "vac.mdp")
open(MDP, "w").write("integrator = md\nnsteps = 0\ncutoff-scheme = Verlet\npbc = xyz\nrvdw = 1.4\nrcoulomb = 1.4\nvdw-modifier = None\n"
                     "coulombtype = Reaction-Field\nepsilon-rf = 1\nDispCorr = no\nconstraints = none\nnstcalcenergy = 1\nnstenergy = 1\n")


def energies(top):
    stem = os.path.join(work, os.path.basename(top)[:-4])
    g = subprocess.run([GMX, "-quiet", "grompp", "-f", MDP, "-c", os.path.join(work, "box.gro"), "-p", top,
                        "-o", stem + ".tpr", "-maxwarn", "10"], capture_output=True, text=True, cwd=work)
    if g.returncode:
        raise SystemExit("grompp failed on " + top + ": " + g.stderr[-300:])
    subprocess.run([GMX, "-quiet", "mdrun", "-s", stem + ".tpr", "-rerun", os.path.join(work, "box.gro"), "-e", stem + ".edr", "-g", stem + ".log",
                    "-nt", "1"], capture_output=True, text=True, cwd=work)
    names = "G96Bond\nG96Angle\nProper-Dih.\nImproper-Dih.\nLJ-14\nCoulomb-14\nLJ-(SR)\nCoulomb-(SR)\n\n"
    o = subprocess.run([GMX, "-quiet", "energy", "-f", stem + ".edr", "-o", stem + "_e.xvg"], input=names, capture_output=True, text=True, cwd=work)
    e = {}
    for l in o.stdout.splitlines():
        m = re.match(r"^(\S.*?)\s{2,}(-?[\d.e+-]+)\s", l)
        if m and m.group(1) in ("G96Bond", "G96Angle", "Proper Dih.", "Improper Dih.", "LJ-14", "Coulomb-14", "LJ (SR)", "Coulomb (SR)"):
            e[m.group(1)] = float(m.group(2))
    return e
er, ec = energies(os.path.join(work, "ref.top")), energies(os.path.join(work, "caps.top"))
worst = max(abs(er[k] - ec.get(k, float("nan"))) / max(1.0, abs(er[k])) for k in er) if er else float("nan")
print("energies (kJ/mol, pdb2gmx / CAPS): " + "  ".join(f"{k} {er[k]:.4f}/{ec.get(k, float('nan')):.4f}" for k in er) + f" · largest relative difference {worst:.1e}")
bad += 0 if worst < 1e-5 else 1
print("same topology and energies" if not bad else f"{bad} differences")
sys.exit(1 if bad else 0)
