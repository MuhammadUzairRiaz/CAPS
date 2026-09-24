#!/usr/bin/env python3
"""Validate a class II (COMPASS / PCFF) CAPS force field against moltemplate and LAMMPS.

usage: compare_class2.py DIR NAME FF.json
  DIR contains NAME.data, NAME.in.init, NAME.in.settings (and NAME.in.charges) written by moltemplate.
  LMP (environment) must be a LAMMPS binary with the CLASS2 package.

Checks:
  1. every bond, angle, dihedral and improper moltemplate generated is in CAPS with the same coefficients,
     including every cross term (bb, ba, mbt, ebt, at, aat, bb13, aa). Cross terms are compared by the atoms
     they belong to, so the check does not depend on the direction an interaction is listed in.
     moltemplate stores each interaction in a canonical (sorted) atom order after matching; where a
     direction-dependent rule matched the other way round, its end-specific coefficients land on the wrong
     atoms. Those interactions are reported as "moltemplate orientation" and are not CAPS errors.
  2. bond-increment charges agree with moltemplate's (moltemplate doubles them when the molecule is written
     with "Data Bond List"; that is detected and reported).
  3. CAPS energies and forces equal LAMMPS's for the data file CAPS writes (class2 styles, lj/class2 +
     coul/dsf), and, when no orientation differences were found, for moltemplate's own data file.
Exit status 0 when 1 (apart from moltemplate orientation), 2 and 3 pass.
"""
import math, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CAPS = os.path.join(HERE, "..", "..", "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
d, name, ffjson = sys.argv[1], sys.argv[2], sys.argv[3]
data, init, settings, charges_in = (os.path.join(d, name + x) for x in (".data", ".in.init", ".in.settings", ".in.charges"))
TOL = 1e-6

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

# ---- moltemplate coefficients: main line plus cross-term lines per type ----
coef = {k: {} for k in ("bond", "angle", "dihedral", "improper")}
for line in open(settings):
    w = line.split("#")[0].split()
    if len(w) >= 3 and w[0] in ("bond_coeff", "angle_coeff", "dihedral_coeff", "improper_coeff"):
        kind, t, rest = w[0].split("_")[0], w[1], w[2:]
        group = "main"
        if re.match(r"[a-z]", rest[0]):
            group, rest = rest[0], rest[1:]
        coef[kind].setdefault(t, {})[group] = [float(x) for x in rest]

Z = lambda n: [0.0] * n
def full_angle(c): return c.get("main", Z(4)), c.get("bb", Z(3)), c.get("ba", Z(4))
def full_dihedral(c):
    return c.get("main", Z(6)), c.get("mbt", Z(4)), c.get("ebt", Z(8)), c.get("at", Z(8)), c.get("aat", Z(3)), c.get("bb13", Z(3))
def full_improper(c): return c.get("main", Z(2)), c.get("aa", Z(6))

# ---- orientation-independent forms: every end- or angle-specific coefficient keyed by its atoms ----
def canon_angle(a, v):
    i, j, k = a
    m, bb, ba = v
    return {"main": tuple(m), "bbM": bb[0], ("bbr", i): bb[1], ("bbr", k): bb[2],
            ("baN", i): ba[0], ("baN", k): ba[1], ("bar", i): ba[2], ("bar", k): ba[3]}

def canon_dihedral(a, v):
    i, j, k, l = a
    m, mbt, ebt, at, aat, bb13 = v
    return {"main": tuple(m), "mbt": tuple(mbt),
            ("ebt", i): tuple(ebt[0:3]), ("ebt", l): tuple(ebt[3:6]), ("ebtr", i): ebt[6], ("ebtr", l): ebt[7],
            ("at", j): tuple(at[0:3]), ("at", k): tuple(at[3:6]), ("att", j): at[6], ("att", k): at[7],
            "aatM": aat[0], ("aatt", j): aat[1], ("aatt", k): aat[2],
            "bb13N": bb13[0], ("bb13r", i): bb13[1], ("bb13r", l): bb13[2]}

def canon_improper(a, v):
    A, B, C, D = a
    m, aa = v
    ac, ad, cd = frozenset((A, C)), frozenset((A, D)), frozenset((C, D))
    # K (mean Wilson angle − χ0)²: invariant under the atom order when χ0 = 0 (always so in COMPASS / PCFF)
    return {"main": (m[0], m[1]), ("th", ac): aa[3], ("th", ad): aa[4], ("th", cd): aa[5],
            ("M", frozenset((ac, cd))): aa[0], ("M", frozenset((ac, ad))): aa[1], ("M", frozenset((ad, cd))): aa[2]}

def same(x, y):
    if x.keys() != y.keys():
        return False
    for k in x:
        a, b = x[k], y[k]
        a = a if isinstance(a, tuple) else (a,)
        b = b if isinstance(b, tuple) else (b,)
        if len(a) != len(b) or any(abs(p - q) > TOL * max(1, abs(q)) for p, q in zip(a, b)):
            return False
    return True

def key(kind, a):
    if kind == "bond": return frozenset(a)
    if kind == "angle": return (a[1], frozenset((a[0], a[2])))
    if kind == "dihedral": return min(tuple(a), tuple(a[::-1]))
    return (a[1], frozenset((a[0], a[2], a[3])))

def canon(kind, a, v):
    if kind == "bond": return {"main": tuple(v)}
    return {"angle": canon_angle, "dihedral": canon_dihedral, "improper": canon_improper}[kind](a, v)

def other_orders(kind, a):
    if kind in ("angle", "dihedral"): return [a[::-1]]
    if kind == "improper":
        A, B, C, D = a
        return [(p, B, q, r) for p, q, r in ((A, D, C), (C, A, D), (C, D, A), (D, A, C), (D, C, A))]
    return []

sec = sections(data)
ref = {k: {} for k in coef}
for kind, secname in (("bond", "Bonds"), ("angle", "Angles"), ("dihedral", "Dihedrals"), ("improper", "Impropers")):
    for w in sec.get(secname, []):
        a = tuple(int(x) for x in w[2:])
        c = coef[kind][w[1]]
        v = {"bond": lambda c: c["main"], "angle": full_angle, "dihedral": full_dihedral, "improper": full_improper}[kind](c)
        ref[kind][key(kind, a)] = (a, v)

# ---- CAPS ----
ffout = os.path.join(d, "caps_" + name + ".data")
cmd = [CAPS, "ff", "apply", data, "--ff", ffjson, "--charges", "types", "--list", "--forces", os.path.join(d, "caps_f.txt"),
       "--cutoff", "10", "-o", ffout, "--allow-missing"]
out = subprocess.run(cmd, capture_output=True, text=True)
if out.returncode not in (0, 3):
    print(out.stdout, out.stderr); sys.exit(2)
for l in out.stdout.splitlines():
    if l.startswith(("note", "missing", "  ")) or "types," in l:
        print("caps: " + l)
caps = {k: {} for k in coef}
qcaps = {}
for line in out.stdout.splitlines():
    w = line.split()
    if not w:
        continue
    if w[0] == "charge":
        qcaps[int(w[1])] = float(w[2]); continue
    kinds = {"bond2": ("bond", 2), "angle2": ("angle", 3), "dihedral2": ("dihedral", 4), "improper2": ("improper", 4)}
    if w[0] not in kinds:
        continue
    kind, na = kinds[w[0]]
    a = tuple(int(x) for x in w[1:1 + na])
    p = [float(x) for x in w[1 + na:]]
    v = {"bond": lambda p: p, "angle": lambda p: (p[0:4], p[4:7], p[7:11]),
         "dihedral": lambda p: (p[0:6], p[6:10], p[10:18], p[18:26], p[26:29], p[29:32]),
         "improper": lambda p: (p[0:2], p[2:8])}[kind](p)
    caps[kind][key(kind, a)] = (a, v)

problems, orientation = 0, 0
for kind in ("bond", "angle", "dihedral", "improper"):
    R, C = ref[kind], caps[kind]
    ok_n = orient_n = 0
    for k, (a, v) in R.items():
        if k not in C:
            problems += 1
            print(f"  MISSING in CAPS {kind} {a}")
            continue
        ca, cv = C[k]
        if same(canon(kind, ca, cv), canon(kind, a, v)):
            ok_n += 1
        elif any(same(canon(kind, ca, cv), canon(kind, o, v)) for o in other_orders(kind, a)):
            orient_n += 1
            print(f"  moltemplate orientation {kind}: moltemplate lists {a}, the rule matched as CAPS lists {ca}")
        else:
            problems += 1
            print(f"  MISMATCH {kind} {a}: moltemplate {v} · CAPS {ca} {cv}")
    for k in C:
        if k not in R:
            problems += 1
            print(f"  EXTRA in CAPS {kind} {C[k][0]}")
    orientation += orient_n
    print(f"{kind:9s} moltemplate {len(R):4d} · identical {ok_n:4d} · moltemplate orientation {orient_n:3d} · CAPS {len(C):4d}")

# ---- physical consistency: cross-term reference lengths against the r0 of the bond they refer to ----
def r0_consistency(table, bonds):
    good = bad = 0
    def chk(a, b, ref):
        nonlocal good, bad
        r0 = bonds.get(frozenset((a, b)))
        if r0 is None or ref == 0: return
        if abs(r0 - ref) < 1e-3: good += 1
        else: bad += 1
    for a, v in table["angle"].values():
        i, j, k = a
        m, bb, ba = v
        chk(i, j, bb[1]); chk(j, k, bb[2]); chk(i, j, ba[2]); chk(j, k, ba[3])
    for a, v in table["dihedral"].values():
        i, j, k, l = a
        m, mbt, ebt, at, aat, bb13 = v
        chk(j, k, mbt[3]); chk(i, j, ebt[6]); chk(k, l, ebt[7]); chk(i, j, bb13[1]); chk(k, l, bb13[2])
    return good, bad
def theta_consistency(table, angles):
    good = bad = 0
    def chk(i, j, k, ref):
        nonlocal good, bad
        t0 = angles.get((j, frozenset((i, k))))
        if t0 is None or ref == 0: return
        if abs(t0 - ref) < 1e-3: good += 1
        else: bad += 1
    for a, v in table["dihedral"].values():
        i, j, k, l = a
        m, mbt, ebt, at, aat, bb13 = v
        chk(i, j, k, at[6]); chk(j, k, l, at[7]); chk(i, j, k, aat[1]); chk(j, k, l, aat[2])
    for a, v in table["improper"].values():
        A, B, C, D = a
        m, aa = v
        chk(A, B, C, aa[3]); chk(A, B, D, aa[4]); chk(C, B, D, aa[5])
    return good, bad
angles_caps = {k: v[1][0][0] for k, v in caps["angle"].items()}
t1, u1 = theta_consistency(ref, angles_caps)
t2, u2 = theta_consistency(caps, angles_caps)
print(f"reference angles equal to the angle's theta0: moltemplate {t1}/{t1 + u1} · CAPS {t2}/{t2 + u2}")
if u2 > u1:
    problems += 1
    print("  CAPS has more cross-term reference angles on the wrong angle than moltemplate")
bonds_caps = {k: v[1][0] for k, v in caps["bond"].items()}
g1, b1 = r0_consistency(ref, bonds_caps)
g2, b2 = r0_consistency(caps, bonds_caps)
print(f"reference lengths equal to the bond's r0: moltemplate {g1}/{g1 + b1} · CAPS {g2}/{g2 + b2}")
if b2 > b1:
    problems += 1
    print("  CAPS has more cross-term reference lengths on the wrong bond than moltemplate")

# ---- charges ----
qmt = {}
if os.path.exists(charges_in):
    for l in open(charges_in):
        w = l.split()
        if len(w) >= 5 and w[0] == "set" and w[1] == "atom":
            qmt[int(w[2])] = float(w[4])
if qmt:
    worst = max(abs(qcaps[i] - q) for i, q in qmt.items())
    worst2 = max(abs(2 * qcaps[i] - q) for i, q in qmt.items())
    if worst < 1e-6:
        print("charges   bond increments identical to moltemplate")
    elif worst2 < 1e-6:
        print("charges   moltemplate's are exactly twice CAPS's: moltemplate applied each bond increment twice\n"
              "          (charge_by_bond.py is given both 'Data Bonds' and 'Data Bond List'); CAPS applies each once")
    else:
        problems += 1
        print(f"charges   MISMATCH: largest difference {worst:.4f} e")
print(f"charges   total {sum(qcaps.values()):+.6f} e")

# ---- LAMMPS ----
def lammps(datafile, header, extra=""):
    tmp = tempfile.mkdtemp()
    # the pair coefficients go to the input (hybrid/overlay needs the sub-style name on each)
    lines, pij, skip = [], [], False
    for l in open(datafile):
        s = l.split("#")[0].strip()
        if l.startswith("PairIJ Coeffs") or l.startswith("Pair Coeffs"):
            skip = l.startswith("PairIJ") and "ij" or "ii"; continue
        if skip:
            if s and re.fullmatch(r"[A-Z][A-Za-z ]+", s):
                skip = False
            elif s:
                w = s.split()
                pij.append((w[0], w[1], w[2], w[3]) if skip == "ij" else (w[0], w[0], w[1], w[2]))
                continue
            else:
                continue
        lines.append(l)
    open(os.path.join(tmp, "d.data"), "w").writelines(lines)
    pc = "\n".join(f"pair_coeff {a} {b} lj/class2 {e} {s}" for a, b, e, s in pij)
    open(os.path.join(tmp, "in"), "w").write(f"""units real
atom_style full
bond_style class2
angle_style class2
dihedral_style class2
improper_style class2
pair_style hybrid/overlay lj/class2 10.0 coul/dsf 0.2 10.0
pair_modify mix sixthpower tail yes
special_bonds lj/coul 0.0 0.0 1.0
read_data {tmp}/d.data
{header}
{pc}
pair_coeff * * coul/dsf
{extra}
thermo_style custom ebond eangle edihed eimp evdwl ecoul
thermo_modify format float %.10f
run 0
write_dump all custom {tmp}/f.dump id fx fy fz modify sort id format float %.10f
""")
    r = subprocess.run([LMP, "-in", os.path.join(tmp, "in"), "-log", os.path.join(tmp, "log"), "-screen", "none"], capture_output=True, text=True)
    log = open(os.path.join(tmp, "log")).read().splitlines()
    try:
        i = next(k for k, l in enumerate(log) if l.split()[:1] == ["E_bond"])
    except StopIteration:
        print("\n".join(log[-15:])); raise
    ev = [float(x) for x in log[i + 1].split()]
    f = {}
    for l in open(os.path.join(tmp, "f.dump")).read().splitlines()[9:]:
        w = l.split()
        f[int(w[0])] = [float(x) for x in w[1:4]]
    return ev, f

energy_line = next(l for l in out.stdout.splitlines() if l.startswith("energy"))
cv = [float(x) for x in re.findall(r"(?:bond|angle|dihedral|improper|vdW|Coulomb) (-?[\d.]+)", energy_line)]
cf = {int(l.split()[0]): [float(x) for x in l.split()[1:]] for l in open(os.path.join(d, "caps_f.txt"))}
names = ["bond", "angle", "dihedral", "improper", "vdW", "Coulomb"]

def report(title, ev, f, strict):
    global problems
    print(f"\n{title}\n{'term':9s} {'CAPS':>16s} {'LAMMPS':>16s}")
    for k, n in enumerate(names):
        bad = n != "Coulomb" and abs(cv[k] - ev[k]) > 1e-4 * max(1, abs(ev[k]))
        if bad and strict:
            problems += 1
        print(f"{n:9s} {cv[k]:16.6f} {ev[k]:16.6f}{'   MISMATCH' if bad else ''}")
    worst = max(math.dist(f[i], cf[i]) for i in f)
    if worst > 1e-3 and strict:
        problems += 1
    print(f"largest per-atom force difference {worst:.2e} kcal/mol/Å" + ("   MISMATCH" if worst > 1e-3 else "") +
          "\n(Coulomb energies differ by the DSF self-energy constant only; forces include Coulomb)")

ev, f = lammps(ffout, "")
report("LAMMPS on the data file CAPS wrote (CAPS's assignment, LAMMPS's class2 styles):", ev, f, True)
qset = "\n".join(f"set atom {i} charge {q:.10f}" for i, q in sorted(qcaps.items()))
# moltemplate's data with its own coefficients; pair lines get the sub-style name for hybrid/overlay
tmpset = tempfile.mktemp()
# With hybrid/overlay, "pair_coeff * * coul/dsf" marks every i-j pair as set, so LAMMPS would never mix the
# lj/class2 cross pairs: write the whole sixth-power matrix from moltemplate's i-i coefficients.
self_lj, mt_set = {}, []
for l in open(settings).read().splitlines():
    w = l.split("#")[0].split()
    if len(w) >= 5 and w[0] == "pair_coeff" and w[1] == w[2]:
        self_lj[int(w[1])] = (float(w[3]), float(w[4]))
    elif not (w and w[0] == "pair_coeff"):
        mt_set.append(l)
for i in sorted(self_lj):
    for j in sorted(self_lj):
        if j < i: continue
        (ei, si), (ej, sj) = self_lj[i], self_lj[j]
        s6 = si ** 6 + sj ** 6
        e = 2 * math.sqrt(ei * ej) * si ** 3 * sj ** 3 / s6 if s6 > 0 else 0.0
        mt_set.append(f"pair_coeff {i} {j} lj/class2 {e:.12g} {(0.5 * s6) ** (1 / 6):.12g}")
open(tmpset, "w").write("\n".join(mt_set) + "\n" + qset + "\n")
ev2, f2 = lammps(data, "", f"include {tmpset}")
print("\n(moltemplate's data file, CAPS charges, moltemplate coefficients" +
      (", some interactions in moltemplate orientation: energies expected to differ)" if orientation else ")"))
report("LAMMPS on moltemplate's data file:", ev2, f2, orientation == 0)
print("\nPASS" if problems == 0 else f"\nFAIL ({problems} problems)")
sys.exit(0 if problems == 0 else 1)
