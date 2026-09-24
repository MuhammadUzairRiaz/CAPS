#!/usr/bin/env python3
"""Validate a CAPS force field converted from moltemplate against moltemplate itself and LAMMPS.

usage: compare_moltemplate.py DIR NAME FF.json
  DIR contains NAME.data, NAME.in.init and NAME.in.settings written by moltemplate.

Checks, for the structure moltemplate built:
  1. every bonded interaction moltemplate generated has the same atoms and the same coefficients in CAPS
     (dihedral and improper terms compared as sets of non-zero Fourier terms);
  2. bond, angle, dihedral, improper and van der Waals energies agree with LAMMPS (pair style swapped for
     lj/cut/coul/dsf with tail corrections so both codes evaluate the same functional form);
  3. per-atom forces agree with LAMMPS.
Exit status 0 only when all three pass.
"""
import math, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CAPS = os.path.join(HERE, "..", "..", "build", "cli", "caps")
LMP = os.environ.get("LMP", "lmp")
d, name, ffjson = sys.argv[1], sys.argv[2], sys.argv[3]
data, init, settings = (os.path.join(d, name + x) for x in (".data", ".in.init", ".in.settings"))

def sections(path):
    out, cur = {}, None
    for line in open(path):
        s = line.split("#")[0].strip()
        if not s:
            continue
        if re.fullmatch(r"[A-Z][A-Za-z ]+", s):
            cur = s; out[cur] = []; continue
        if cur:
            out[cur].append(s.split())
    return out

sec = sections(data)
# default styles from the init file (used when coefficient lines carry no style keyword)
default = {}
hybrid_pair = False
for line in open(init):
    w = line.split("#")[0].split()
    if len(w) >= 2 and w[0] in ("bond_style", "angle_style", "dihedral_style", "improper_style"):
        default[w[0].split("_")[0]] = w[2] if w[1] == "hybrid" and len(w) > 2 else w[1]
    if len(w) >= 2 and w[0] == "pair_style" and w[1] == "hybrid":
        hybrid_pair = True
coef = {"bond": {}, "angle": {}, "dihedral": {}, "improper": {}}
for line in open(settings):
    w = line.split("#")[0].split()
    if len(w) >= 3 and w[0] in ("bond_coeff", "angle_coeff", "dihedral_coeff", "improper_coeff"):
        kind = w[0].split("_")[0]
        rest = w[2:]
        style = None
        if rest and re.match(r"[a-z]", rest[0]):
            style, rest = rest[0], rest[1:]
        coef[kind][w[1]] = (style or default.get(kind), [float(x) for x in rest])

def fourier_terms(style, p):
    """Non-zero (K, n, delta) terms of a dihedral / improper."""
    t = []
    if style == "fourier":
        m = int(p[0])
        t = [(p[1 + 3 * q], int(p[2 + 3 * q]), p[3 + 3 * q] % 360) for q in range(m)]
    elif style == "opls":
        t = [(0.5 * p[0], 1, 0), (0.5 * p[1], 2, 180), (0.5 * p[2], 3, 0), (0.5 * p[3], 4, 180)]
    elif style == "cvff":
        t = [(p[0], int(p[2]), 0 if p[1] > 0 else 180)]
    elif style == "harmonic" and len(p) == 3:
        t = [(p[0], int(p[2]), 0 if p[1] > 0 else 180)]
    return sorted((round(k, 6), n, round(dd, 3)) for k, n, dd in t if abs(k) > 1e-12)

# moltemplate's interactions
ref = {"bond": {}, "angle": {}, "dihedral": {}, "improper": {}}
for kind, key in (("bond", "Bonds"), ("angle", "Angles"), ("dihedral", "Dihedrals"), ("improper", "Impropers")):
    for w in sec.get(key, []):
        atoms = tuple(int(x) for x in w[2:])
        style, p = coef[kind][w[1]]
        ref[kind][atoms] = (style, p)

# CAPS interactions
out = subprocess.run([CAPS, "ff", "apply", data, "--ff", ffjson, "--charges", "keep", "--list", "--forces", os.path.join(d, "caps_f.txt"),
                      "--cutoff", "10"], capture_output=True, text=True)
if out.returncode not in (0, 3):
    print(out.stdout, out.stderr); sys.exit(2)
caps = {"bond": {}, "angle": {}, "dihedral": {}, "improper": {}}
for line in out.stdout.splitlines():
    w = line.split()
    if not w or w[0] not in caps:
        continue
    kind = w[0]
    na = {"bond": 2, "angle": 3, "dihedral": 4, "improper": 4}[kind]
    atoms = tuple(int(x) for x in w[1:1 + na])
    vals = w[1 + na:]
    if kind == "dihedral":
        caps[kind].setdefault(atoms, []).append((round(float(vals[0]), 6), int(vals[1]), round(float(vals[2]) % 360, 3)))
    elif kind == "improper":
        if vals[0] == "cvff":
            caps[kind].setdefault(atoms, []).append((round(float(vals[1]), 6), int(vals[2]), round(float(vals[3]) % 360, 3)))
        else:
            caps[kind][atoms] = ("harmonic", float(vals[1]), float(vals[2]))
    else:
        caps[kind][atoms] = [float(x) for x in vals]
energy_line = next(l for l in out.stdout.splitlines() if l.startswith("energy"))

def canon(kind, a):
    if kind in ("bond", "angle", "dihedral"):
        return min(a, a[::-1])
    return a

problems = 0
for kind in ("bond", "angle", "dihedral", "improper"):
    R = {canon(kind, a): v for a, v in ref[kind].items()}
    C = {canon(kind, a): v for a, v in caps[kind].items()}
    checked = 0
    for a, (style, p) in R.items():
        if kind in ("bond", "angle"):
            c = C.get(a)
            ok = c is not None and all(abs(x - y) < 1e-6 for x, y in zip(c, p[:2]))
        elif kind == "dihedral":
            want = fourier_terms(style, p)
            ok = sorted(C.get(a, [])) == want
        else:
            if style == "harmonic" and len(p) == 2:
                c = C.get(a)
                ok = c is not None and c[0] == "harmonic" and abs(c[1] - p[0]) < 1e-6 and abs(c[2] - p[1]) < 1e-6
            else:
                ok = sorted(C.get(a, [])) == fourier_terms(style, p)
        checked += 1
        if not ok:
            problems += 1
            print(f"  MISMATCH {kind} {a}: moltemplate {style} {p} · CAPS {C.get(a)}")
    extra = [a for a in C if a not in R and (kind not in ("dihedral", "improper") or C[a])]
    for a in extra:
        problems += 1
        print(f"  EXTRA in CAPS {kind} {a}: {C[a]}")
    print(f"{kind:9s} moltemplate {len(R):4d} · checked {checked:4d} · CAPS {len(C):4d}")

# LAMMPS energies and forces with the same pair form as CAPS
tmp = tempfile.mkdtemp()
init_txt = open(init).read()
init_txt = re.sub(r"(?m)^\s*pair_style.*$", "pair_style " + ("hybrid " if hybrid_pair else "") + "lj/cut/coul/dsf 0.2 10.0", init_txt)
init_txt = re.sub(r"(?m)^\s*kspace_style.*$", "", init_txt)
init_txt = re.sub(r"(?m)^\s*pair_modify\s+mix\s+(\w+).*$", r"pair_modify mix \1 tail yes", init_txt)
if "pair_modify" not in init_txt:
    init_txt += "\npair_modify tail yes\n"
set_txt = re.sub(r"(pair_coeff\s+\S+\s+\S+\s+)lj/\S+", r"\1lj/cut/coul/dsf", open(settings).read())
open(os.path.join(tmp, "in.init"), "w").write(init_txt)
open(os.path.join(tmp, "in.settings"), "w").write(set_txt)
open(os.path.join(tmp, "in.run"), "w").write(f"""include {tmp}/in.init
read_data {os.path.abspath(data)}
include {tmp}/in.settings
thermo_style custom ebond eangle edihed eimp evdwl ecoul
thermo_modify format float %.8f
run 0
write_dump all custom {tmp}/f.dump id fx fy fz modify sort id format float %.10f
""")
p = subprocess.run([LMP, "-in", os.path.join(tmp, "in.run"), "-log", os.path.join(tmp, "log"), "-screen", "none"], capture_output=True, text=True)
log = open(os.path.join(tmp, "log")).read().splitlines()
i = next(k for k, l in enumerate(log) if l.split()[:1] == ["E_bond"])
lv = [float(x) for x in log[i + 1].split()]
cv = [float(x) for x in re.findall(r"(?:bond|angle|dihedral|improper|vdW|Coulomb) (-?[\d.]+)", energy_line)]
names = ["bond", "angle", "dihedral", "improper", "vdW", "Coulomb"]
print(f"{'term':9s} {'CAPS':>16s} {'LAMMPS':>16s}")
for k, n in enumerate(names):
    flag = "" if n == "Coulomb" or abs(cv[k] - lv[k]) < 1e-4 * max(1, abs(lv[k])) else "   MISMATCH"
    if flag:
        problems += 1
    print(f"{n:9s} {cv[k]:16.6f} {lv[k]:16.6f}{flag}")
cf = {int(l.split()[0]): [float(x) for x in l.split()[1:]] for l in open(os.path.join(d, "caps_f.txt"))}
worst = 0
for l in open(os.path.join(tmp, "f.dump")).read().splitlines()[9:]:
    w = l.split()
    worst = max(worst, math.dist([float(x) for x in w[1:4]], cf[int(w[0])]))
print(f"largest per-atom force difference {worst:.2e} kcal/mol/Å (Coulomb energy differs by the DSF constant only)")
if worst > 1e-3:
    problems += 1
print("PASS" if problems == 0 else f"FAIL ({problems} problems)")
sys.exit(0 if problems == 0 else 1)
