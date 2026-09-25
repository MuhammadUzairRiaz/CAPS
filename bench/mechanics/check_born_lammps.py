#!/usr/bin/env python3
"""The Born term of CAPS's fluctuation elastic constants against LAMMPS compute born/matrix, on one configuration.

  bonded terms (bonds, angles, Fourier torsions, cvff impropers): LAMMPS numdiff (finite differences of the virial with
      its stress add-on terms) with the pair style switched off. Bonded energies are smooth, so the two finite
      differences must agree closely.
  bonded + Lennard-Jones pairs (charges off): LAMMPS numdiff with a small step (1e-7), where pairs crossing the
      truncated cut-off are rare. CAPS takes finite differences with the pair set frozen (the analytic definition).
  argon crystal: LAMMPS numdiff against the lattice sum Σ(φ'' − φ'/r)X⁴/r², and LAMMPS's analytic pair Born term,
      which differs from both (reported).

usage: check_born_lammps.py [DATA] [--ff ID] [--keep DIR]   (default: out/ps500.data, gaff-amber25)
Needs LMP (default ~/.local/bin/lmp) with EXTRA-COMPUTE, MOLECULE, EXTRA-MOLECULE.
"""
import os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/.local/bin/lmp"))
args = [a for a in sys.argv[1:] if not a.startswith("--")]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
data = args[0] if args else os.path.join(ROOT, "out", "ps500.data")
ff = os.path.join(ROOT, "data", "forcefields", arg("--ff", "gaff-amber25") + ".json")
work = arg("--keep", "") or tempfile.mkdtemp()
os.makedirs(work, exist_ok=True)
NAMES = ["C11", "C22", "C33", "C44", "C55", "C66", "C12", "C13", "C14", "C15", "C16", "C23", "C24", "C25", "C26", "C34", "C35", "C36",
         "C45", "C46", "C56"]


def caps_born(which):
    r = subprocess.run([CAPS, "ff", "apply", data, "--ff", ff, "--charges", "gasteiger", "-o", os.path.join(work, "case.data"),
                        "--lammps-input", os.path.join(work, "case.in"), "--born", which, "--no-tail"], capture_output=True, text=True, check=True)
    m = re.search(r"born matrix \(kcal/mol, \w+\): (.*)", r.stdout)
    return list(map(float, m.group(1).split()))


def lammps_born(which, delta=None):
    txt = open(os.path.join(work, "case.in")).read()
    txt = re.sub(r"thermo_style.*\n|thermo_modify.*\n|run 0\n", "", txt)
    txt = txt.replace("pair_modify shift yes\n", "")
    if which == "bonded":
        txt = re.sub(r"pair_style .*", "pair_style zero 10.0 nocoeff", txt)
        txt += "compute vir all pressure NULL virial\ncompute born all born/matrix numdiff 1.0e-6 vir\n"
    elif delta:
        txt = re.sub(r"pair_style .*", "pair_style lj/cut 10.0", txt)
        txt += f"compute vir all pressure NULL virial\ncompute born all born/matrix numdiff {delta} vir\n"
    else:
        txt = re.sub(r"pair_style .*", "pair_style lj/cut 10.0", txt)
        txt += "compute born all born/matrix pair\n"
    txt += "thermo_style custom step " + " ".join(f"c_born[{k}]" for k in range(1, 22)) + "\nthermo_modify format float %.12g\nrun 0\n"
    inp = os.path.join(work, f"in.{which}")
    open(inp, "w").write(txt)
    r = subprocess.run([LMP, "-in", inp, "-log", "none"], capture_output=True, text=True, cwd=work)
    if r.returncode:
        raise RuntimeError((r.stdout + r.stderr).strip().splitlines()[-1])
    lines = r.stdout.splitlines()
    k = next(i for i, l in enumerate(lines) if l.split()[:1] == ["Step"])
    return list(map(float, lines[k + 1].split()[1:]))


ARGON = """units real
atom_style atomic
lattice fcc 5.40
region box block 0 4 0 4 0 4
create_box 1 box
create_atoms 1 box
mass 1 39.948
pair_style lj/cut 7.9
pair_coeff 1 1 0.238 3.405
compute vir all pressure NULL virial
compute born all {mode}
thermo_style custom step c_born[1] c_born[7] c_born[4]
thermo_modify format float %.12g
run 0
"""


def argon(mode):
    inp = os.path.join(work, "in.argon")
    open(inp, "w").write(ARGON.format(mode=mode))
    out = subprocess.run([LMP, "-in", inp, "-log", "none"], capture_output=True, text=True, cwd=work).stdout.splitlines()
    k = next(i for i, l in enumerate(out) if l.split()[:1] == ["Step"])
    return list(map(float, out[k + 1].split()[1:]))


def argon_lattice_sum():
    """(φ'' − φ'/r) X⁴/r² over the pairs of the fcc crystal: C11, C12, C44 (kcal/mol)."""
    import math
    a, n, eps, sig, rc = 5.40, 4, 0.238, 3.405, 7.9
    L = a * n
    basis = [(0, 0, 0), (.5, .5, 0), (.5, 0, .5), (0, .5, .5)]
    P = [((i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a) for i in range(n) for j in range(n) for k in range(n) for b in basis]
    c = [0.0, 0.0, 0.0]
    for p in P[1:]:
        d = [p[k] - P[0][k] for k in range(3)]
        d = [v - L * round(v / L) for v in d]
        r = math.sqrt(sum(v * v for v in d))
        if r >= rc:
            continue
        s6 = (sig / r) ** 6
        d1, d2 = 4 * eps * (-12 * s6 * s6 + 6 * s6) / r, 4 * eps * (156 * s6 * s6 - 42 * s6) / r ** 2
        w = (d2 - d1 / r) / r ** 2 * len(P) / 2
        c[0] += w * d[0] ** 4
        c[1] += w * d[0] ** 2 * d[1] ** 2
        c[2] += w * d[1] ** 2 * d[2] ** 2
    return c


fails = 0
def report(ok, text):
    global fails
    fails += not ok
    print(f"{'ok  ' if ok else 'FAIL'} {text}")

bonded_c, bonded_l = caps_born("bonded"), lammps_born("bonded")
scale = max(map(abs, bonded_l))
worst = max(abs(a - b) for a, b in zip(bonded_c, bonded_l)) / scale
report(worst < 1e-6, f"bonded terms (bonds, angles, Fourier torsions, cvff impropers) vs LAMMPS numdiff 1e-6: largest difference {worst:.1e} of C11")
lj_c = caps_born("lj")
total_c = [a + b for a, b in zip(bonded_c, lj_c)]
total_l = lammps_born("lj", 1e-7)   # numdiff differentiates the whole virial: bonded + LJ
worst = max(abs(a - b) for a, b in zip(total_c, total_l)) / max(map(abs, total_l))
report(worst < 1e-5, f"bonded + Lennard-Jones vs LAMMPS numdiff 1e-7: largest difference {worst:.1e} "
       f"(C11 {total_c[0]:.2f} vs {total_l[0]:.2f} kcal/mol; larger steps let pairs cross the truncated cut-off)")
an = lammps_born("lj")
print(f"     LAMMPS analytic pair term for the same LJ pairs: C11 {an[0]:.2f}, CAPS {lj_c[0]:.2f} kcal/mol ({100 * (an[0] / lj_c[0] - 1):+.1f} %)")
ls, nd, pa = argon_lattice_sum(), argon("born/matrix numdiff 1.0e-6 vir"), argon("born/matrix")
report(abs(nd[0] - ls[0]) < 1e-6 * ls[0] and abs(nd[1] - ls[1]) < 1e-6 * ls[1],
       f"argon fcc (a = 5.40 Å): LAMMPS numdiff C11 {nd[0]:.4f}, C12 {nd[1]:.4f} = lattice sum (φ'' − φ'/r) {ls[0]:.4f}, {ls[1]:.4f}")
print(f"     LAMMPS analytic born/matrix gives C11 {pa[0]:.4f}, C12 {pa[1]:.4f}: its pair prefactor ½φ'' − φ'/r is summed over a full "
      f"neighbour list, i.e. φ'' − 2φ'/r per pair; its own numdiff disagrees. CAPS follows the numdiff and the lattice sum.")
print(f"(work: {work})")
sys.exit(1 if fails else 0)
