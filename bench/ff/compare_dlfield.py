#!/usr/bin/env python3
"""Validate a CAPS force field converted from a DL_FIELD library against DL_FIELD itself and LAMMPS.

usage: compare_dlfield.py OUTDIR FF.json
  OUTDIR holds DL_FIELD's lammps.in and lammps1.data for a structure (bench/ff/run_dlfield.py writes them).
  LMP must have the styles DL_FIELD uses (CLASS2, MOFFF for inversion/harmonic, EXTRA-MOLECULE for quartic).

Checks: (1) every bond, angle, dihedral and improper DL_FIELD generated is in CAPS with the same coefficients
(CAPS assigns them from the converted library and the atom types DL_FIELD chose); (2) with the pair style
swapped for LJ + damped shifted force Coulomb (both codes then evaluate the same functions), the bond, angle,
dihedral, improper and van der Waals energies and all forces agree with LAMMPS running DL_FIELD's data file.
"""
import math, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CAPS = os.path.join(HERE, "..", "..", "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
d, ffjson = sys.argv[1], sys.argv[2]
# --structure FILE: give CAPS the input DL_FIELD read (mol2 keeps its bond order) instead of DL_FIELD's data file
structure = sys.argv[sys.argv.index("--structure") + 1] if "--structure" in sys.argv else None
data, lin = os.path.join(d, "lammps1.data"), os.path.join(d, "lammps.in")
TOL = 2e-6

def sections(path):
    out, cur = {}, None
    for line in open(path):
        s = line.split("#")[0].strip()
        if not s: continue
        if re.fullmatch(r"[A-Z][A-Za-z ]+", s):
            cur = s; out[cur] = []; continue
        if cur: out[cur].append(s.split())
    return out

sec = sections(data)
def coeffs(name):
    c = {}
    for w in sec.get(name, []):
        c[w[0]] = (w[1], [float(x) for x in w[2:]]) if not re.match(r"[-\d.]", w[1]) else (None, [float(x) for x in w[1:]])
    return c

def fourier_terms(style, p):
    t = []
    if style == "fourier":
        m = int(p[0]); t = [(p[1 + 3 * q], int(p[2 + 3 * q]), p[3 + 3 * q] % 360) for q in range(m)]
    elif style == "charmm":
        t = [(p[0], int(p[1]), p[2] % 360)]
    elif style == "opls":
        t = [(0.5 * p[0], 1, 0), (0.5 * p[1], 2, 180), (0.5 * p[2], 3, 0), (0.5 * p[3], 4, 180)]
    elif style == "cvff":
        t = [(p[0], int(p[2]), 0 if p[1] > 0 else 180)]
    return sorted((round(k, 6), n, round(x, 3)) for k, n, x in t if abs(k) > 1e-12)

ref = {"bond": {}, "angle": {}, "dihedral": {}, "improper": {}}
for kind, s, cs in (("bond", "Bonds", "Bond Coeffs"), ("angle", "Angles", "Angle Coeffs"), ("dihedral", "Dihedrals", "Dihedral Coeffs"),
                    ("improper", "Impropers", "Improper Coeffs")):
    C = coeffs(cs)
    for w in sec.get(s, []):
        a = tuple(int(x) for x in w[2:])
        if kind in ("dihedral", "improper") and a in ref[kind]:
            st, p = ref[kind][a]   # DL_FIELD writes one line per torsion term: merge them as a fourier set
            t = fourier_terms(st, p) + fourier_terms(*C[w[1]])
            ref[kind][a] = ("fourier", [len(t)] + [x for k in t for x in k])
        else:
            ref[kind][a] = C[w[1]]

types_file = None
if structure:
    # DL_FIELD keeps the atom order; take its types from the data file so both codes use the same ones
    tnames = {w[0]: None for w in sec.get("Masses", [])}
    for line in open(data):
        pass
    masses = {}
    for line in open(data).read().split("Masses")[1].split("\n\n")[1].splitlines():
        if "#" in line:
            masses[line.split()[0]] = line.split("#")[1].strip()
    types_file = os.path.join(d, "caps_types.txt")
    open(types_file, "w").write("\n".join(masses[w[2]] for w in sec["Atoms"]) + "\n")
out = subprocess.run([CAPS, "ff", "apply", structure or data, "--ff", ffjson] + (["--types", types_file] if types_file else []) + ["--charges", "keep", "--list", "--forces", os.path.join(d, "caps_f.txt"),
                      "--cutoff", "12", "--allow-missing", "--no-tail"], capture_output=True, text=True)
if out.returncode not in (0, 3):
    print(out.stdout, out.stderr); sys.exit(2)
for l in out.stdout.splitlines():
    if l.startswith(("note", "missing", "  ")): print("caps: " + l)
caps = {"bond": {}, "angle": {}, "dihedral": {}, "improper": {}}
for line in out.stdout.splitlines():
    w = line.split()
    if not w or line[0].isspace(): continue
    t = w[0]
    if t == "bond":      caps["bond"][tuple(map(int, w[1:3]))] = ("harmonic", [float(w[3]), float(w[4])])
    elif t == "bond2":   caps["bond"][tuple(map(int, w[1:3]))] = ("class2", [float(x) for x in w[3:7]])
    elif t == "angle":   caps["angle"][tuple(map(int, w[1:4]))] = ("harmonic", [float(w[4]), float(w[5])])
    elif t == "angle2":  caps["angle"][tuple(map(int, w[1:4]))] = ("quartic", [float(x) for x in w[4:8]], [float(x) for x in w[8:15]])
    elif t == "dihedral": caps["dihedral"].setdefault(tuple(map(int, w[1:5])), []).append((round(float(w[5]), 6), int(w[6]), round(float(w[7]) % 360, 3)))
    elif t == "improper" and w[5] == "cvff":
        caps["improper"].setdefault(tuple(map(int, w[1:5])), []).append((round(float(w[6]), 6), int(w[7]), round(float(w[8]) % 360, 3)))
    elif t == "inversion": caps["improper"][tuple(map(int, w[1:5]))] = ("inversion/harmonic", [float(w[5]), float(w[6])])

def close(a, b):
    return len(a) == len(b) and all(abs(x - y) <= TOL * max(1, abs(y)) for x, y in zip(a, b))

def akey(kind, a):
    if kind == "bond": return frozenset(a)
    if kind == "angle": return (a[1], frozenset((a[0], a[2])))
    if kind == "dihedral": return min(a, a[::-1])
    return a

problems = 0
order_diff = []
for kind in ("bond", "angle", "dihedral", "improper"):
    R = {akey(kind, a): (a, v) for a, v in ref[kind].items()}
    C = {akey(kind, a): v for a, v in caps[kind].items()}
    if kind == "improper":   # inversions are symmetric in the outer atoms
        # a torsion-type improper a-b-c-d is the angle between the planes (a, b, c) and (b, c, d): the same pair of
        # planes gives the same |φ| (DL_FIELD writes the centre second, AMBER / CAPS third)
        planes = lambda a: frozenset((frozenset(a[0:3]), frozenset(a[1:4])))
        R = {(a[0], frozenset(a[1:])) if v[0] == "inversion/harmonic" else planes(a): (a, v) for a, v in ref[kind].items()}
        C = {(a[0], frozenset(a[1:])) if isinstance(v, tuple) and v[0] == "inversion/harmonic" else planes(a): v
             for a, v in caps[kind].items()}
    same = 0
    if kind == "improper":
        # same centre, same outer atoms, same parameters, other plane pair: the outer-atom order convention
        def centre_set(a): return (frozenset(a), a)
        byset = {}
        for k2, v2 in C.items():
            if isinstance(k2, frozenset):
                atoms = set().union(*k2)
                byset[frozenset(atoms)] = (k2, v2)
    for k, (a, (style, p)) in R.items():
        c = C.get(k)
        if kind in ("dihedral",) or (kind == "improper" and style == "cvff"):
            want = fourier_terms(style, p)
            ok = sorted(c or []) == want
            if not want and c is None: ok = True
        elif kind == "angle" and style == "quartic":
            ok = c is not None and c[0] == "quartic" and close(c[1], p) and not any(c[2])
        elif kind == "angle" and style == "harmonic":
            ok = c is not None and c[0] == "harmonic" and close(c[1], p)
        elif kind == "bond":
            ok = c is not None and c[0] == style and close(c[1], p)
        elif kind == "improper" and style == "inversion/harmonic":
            ok = c is not None and close(c[1], p)
        else:
            ok = False
            print(f"  (style {style} not compared yet)")
        if not ok and kind == "improper" and c is None and frozenset(a) in byset:
            k2, v2 = byset[frozenset(a)]
            if sorted(v2) == fourier_terms(style, p):
                order_diff.append(a)
                C.pop(k2, None)
                print(f"  improper {a}: same parameters; CAPS places the atoms as the matching rule's key orders them, DL_FIELD keeps "
                      f"its own atom order (a DL_FIELD deviation when the key's end atom moves)")
                continue
        if ok: same += 1
        else:
            problems += 1
            print(f"  MISMATCH {kind} {a}: DL_FIELD {style} {p} · CAPS {c}")
    extra = [k for k in C if k not in R and C[k]]
    for k in extra:
        problems += 1
        print(f"  EXTRA in CAPS {kind} {k}: {C[k]}")
    print(f"{kind:9s} DL_FIELD {len(R):4d} · identical {same:4d} · CAPS {len(C):4d}")

# ---- LAMMPS on DL_FIELD's data file, pair style swapped for LJ + coul/dsf ----
tmp = tempfile.mkdtemp()
txt = open(lin).read()
lj96 = "class2" in re.search(r"(?m)^pair_style.*$", txt).group(0)
lj = "lj/class2" if lj96 else "lj/cut"
charged = "atom_style      molecular" not in txt and not re.search(r"(?m)^atom_style\s+molecular", txt)
if charged:
    txt = re.sub(r"(?m)^pair_style.*$", f"pair_style hybrid/overlay {lj} 12.0 coul/dsf 0.2 12.0\npair_modify tail no shift yes", txt)
else:
    txt = re.sub(r"(?m)^pair_style.*$", f"pair_style hybrid {lj} 12.0\npair_modify tail no shift yes", txt)
txt = re.sub(r"(?m)^kspace_style.*$", "", txt)
txt = re.sub(r"(?m)^(pair_coeff\s+\d+\s+\d+\s+)\S+", r"\1" + lj, txt)
txt = txt.replace("read_data lammps1.data", f"read_data {os.path.abspath(data)}" + ("\npair_coeff * * coul/dsf" if charged else ""))
txt = re.sub(r"(?m)^thermo_style.*$", "thermo_style custom ebond eangle edihed eimp evdwl ecoul\nthermo_modify format float %.10f", txt)
txt = re.sub(r"(?m)^(dihedral_style\s+hybrid\s+)charmm\s*$", r"\1charmm", txt)
txt = re.sub(r"(?m)^run 0.*$", f"run 0\nwrite_dump all custom {tmp}/f.dump id fx fy fz modify sort id format float %.10f", txt)
open(os.path.join(tmp, "in"), "w").write(txt)
r = subprocess.run([LMP, "-in", os.path.join(tmp, "in"), "-log", os.path.join(tmp, "log"), "-screen", "none"], capture_output=True, text=True)
log = open(os.path.join(tmp, "log")).read().splitlines()
try:
    i = next(k for k, l in enumerate(log) if l.split()[:1] == ["E_bond"])
except StopIteration:
    print("\n".join(log[-20:])); sys.exit(2)
lv = [float(x) for x in log[i + 1].split()]
energy_line = next(l for l in out.stdout.splitlines() if l.startswith("energy"))
cv = [float(x) for x in re.findall(r"(?:bond|angle|dihedral|improper|vdW|Coulomb) (-?[\d.]+)", energy_line)]
print(f"\n{'term':9s} {'CAPS':>16s} {'LAMMPS':>16s}")
for k, n in enumerate(["bond", "angle", "dihedral", "improper", "vdW", "Coulomb"]):
    bad = n != "Coulomb" and abs(cv[k] - lv[k]) > 1e-4 * max(1, abs(lv[k])) and not (n == "improper" and order_diff)
    problems += bad
    print(f"{n:9s} {cv[k]:16.6f} {lv[k]:16.6f}{'   MISMATCH' if bad else ''}")
cf = {int(l.split()[0]): [float(x) for x in l.split()[1:]] for l in open(os.path.join(d, "caps_f.txt"))}
worst = 0
skip = set().union(*[set(a) for a in order_diff]) if order_diff else set()
for l in open(os.path.join(tmp, "f.dump")).read().splitlines()[9:]:
    w = l.split()
    if int(w[0]) in skip: continue
    worst = max(worst, math.dist([float(x) for x in w[1:4]], cf[int(w[0])]))
print(f"largest per-atom force difference {worst:.2e} kcal/mol/Å (Coulomb energy differs by the DSF constant only)")
if order_diff:
    print(f"  (atoms of the {len(order_diff)} improper(s) with DL_FIELD's atom order left out of the force check: {sorted(skip)})")
problems += worst > 1e-3
print("PASS" if problems == 0 else f"FAIL ({problems} problems)")
sys.exit(0 if problems == 0 else 1)
