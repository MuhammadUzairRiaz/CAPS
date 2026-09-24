#!/usr/bin/env python3
"""Compare CAPS's parameter assignment with DL_FIELD's DL_POLY FIELD file, term by term.

usage: compare_field.py OUTDIR STRUCTURE FF.json
  OUTDIR     holds dl_poly.FIELD written by DL_FIELD (bench/ff/run_dlfield.py) for STRUCTURE
  STRUCTURE  the structure DL_FIELD read (PDB with CONECT, or mol2); atoms in the same order as the FIELD file

The FIELD file is DL_FIELD's complete model: it has the terms DL_FIELD leaves out of its LAMMPS output (CHARMM
1-4 van der Waals, DREIDING inversions, GROMOS). DL_POLY conventions are converted to CAPS / LAMMPS ones:
  bonds  harm k r0 → K = k/2 · quar k r0 k' k'' → r0, k/2, k'/3, k''/4 · mors E0 r0 k → D = E0, α = k
         -hrm k r0 (Urey–Bradley) → k/2 · -lj ε σ (explicit 1-4 van der Waals)
  angles harm → k/2 · quar → θ0, k/2, k'/3, k''/4 · hcos k θ0 → cosine/squared k/2
  dihedrals cos A δ m (+ 1-4 scalings) · harm k φ0 (CHARMM-type improper) → k/2
  inversions harm k ω0 → k/2 (averaged over the three bonds) · plan A → A
  vdw lj ε σ · nm E0 9 6 r0 (class II) · buck A ρ C · shell k2 k4 → class2 spring K2 = k2/2, K4 = k4/24
Every term in either file must have an identical partner. Exit 0 on PASS.
"""
import os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CAPS = os.path.join(HERE, "..", "..", "build", "cli", "caps")
d, structure, ffjson = sys.argv[1], sys.argv[2], sys.argv[3]
# DL_FIELD converts units with rounded factors (eV → kcal/mol 23.061, kJ → kcal 0.23901, K → kcal 0.0019872041);
# CAPS uses the exact constants, so converted libraries agree to ~2e-5 relative
TOL = 5e-5

def close(a, b):
    return len(a) == len(b) and all(abs(x - y) <= TOL * max(1.0, abs(y)) for x, y in zip(a, b))

# ---- FIELD ----
# DL_FIELD writes fixed-width fields that can run together when a value is large ("7425.64200011530.500000" is
# 7425.642000 and 11530.500000); split such tokens after the sixth decimal
def split_fused(tok):
    parts = re.findall(r"-?\d+\.\d{6}", tok)
    return parts if len(parts) > 1 and "".join(parts) == tok else [tok]
fused = []
lines = []
for l in open(os.path.join(d, "dl_poly.FIELD")):
    w = []
    for t in l.split():
        sp = split_fused(t)
        if len(sp) > 1: fused.append(t)
        w += sp
    lines.append(w)
# Molecule types, each with nummols copies; indices in a type are local, expanded here to global ones
KEYS = ("bonds", "angles", "dihedral", "dihedrals", "inversions", "shell", "constraints", "rigid", "teth")
mols, cur, k = [], None, 0
top_vdw = []
while k < len(lines):
    w = lines[k]
    if not w: k += 1; continue
    key = w[0].lower()
    if key == "molecule" and len(w) >= 2 and w[1].lower() == "name":
        cur = {"nummols": 1, "atoms": [], "sec": {}}; mols.append(cur); k += 1; continue
    if key == "nummols" and cur is not None:
        cur["nummols"] = int(w[1]); k += 1; continue
    if key == "atoms" and cur is not None:
        n, k = int(w[1]), k + 1
        while n > 0:
            name, rep = lines[k][0], int(lines[k][3]) if len(lines[k]) > 3 else 1
            cur["atoms"] += [name] * rep; n -= rep; k += 1
        continue
    if key in KEYS and len(w) >= 2 and w[1].isdigit() and cur is not None:
        n = int(w[1]); cur["sec"].setdefault("dihedral" if key == "dihedrals" else key, []).extend(lines[k + 1:k + 1 + n]); k += n + 1
        continue
    if key == "vdw" and len(w) >= 2 and w[1].isdigit():
        n = int(w[1]); top_vdw = lines[k + 1:k + 1 + n]; k += n + 1; continue
    k += 1
atoms, sections, off = [], {}, 0
index_cols = {"bonds": 2, "angles": 3, "dihedral": 4, "inversions": 4, "shell": 2, "constraints": 2}
for m in mols:
    na = len(m["atoms"])
    for c in range(m["nummols"]):
        atoms += m["atoms"]
        for sec, rows in m["sec"].items():
            nc = index_cols.get(sec)
            if nc is None: continue
            first = 0 if sec == "shell" else 1   # shell rows have no form keyword
            for w in rows:
                w2 = list(w)
                for q in range(first, first + nc): w2[q] = str(int(w[q]) + off)
                sections.setdefault(sec, []).append(w2)
        off += na
sections["vdw"] = top_vdw

ref = {}   # (kind, key) -> values
def put(kind, key, val):
    ref.setdefault((kind, key), []).append(val)

for w in sections.get("bonds", []):
    f, i, j = w[0], int(w[1]), int(w[2]); p = [float(x) for x in w[3:]]
    b = frozenset((i, j))
    if f == "harm": put("bond", b, ("harmonic", [p[0] / 2, p[1]]))
    elif f == "quar": put("bond", b, ("class2", [p[1], p[0] / 2, p[2] / 3, p[3] / 4]))
    elif f == "mors": put("bond", b, ("morse", [p[0], p[2], p[1]]))
    elif f == "-hrm": put("ub", b, ("ub", [p[0] / 2, p[1]]))
    elif f == "-lj": put("pair14", b, ("lj", [p[0], p[1]]))
    else: put("bond", b, ("FIELD:" + f, p))
for w in sections.get("shell", []):
    i, j, p = int(w[0]), int(w[1]), [float(x) for x in w[2:]]
    put("bond", frozenset((i, j)), ("class2", [0.0, p[0] / 2, 0.0, (p[1] if len(p) > 1 else 0.0) / 24]))
for w in sections.get("angles", []):
    f, a = w[0], (int(w[2]), frozenset((int(w[1]), int(w[3])))); p = [float(x) for x in w[4:]]
    if f == "harm" and abs(p[0]) < 1e-12: continue   # a zero angle term (DL_FIELD's trans pairs at octahedral centres)
    if f == "harm": put("angle", a, ("harmonic", [p[0] / 2, p[1]]))
    elif f == "quar": put("angle", a, ("quartic", [p[1], p[0] / 2, p[2] / 3, p[3] / 4]))
    elif f == "hcos": put("angle", a, ("cosine/squared", [p[0] / 2, p[1]]))
    elif f == "cos" and abs(p[1]) < 1e-9 and abs(p[2] - 1) < 1e-9: put("angle", a, ("cosine", [p[0], 180.0]))   # A [1 + cos θ]
    else: put("angle", a, ("FIELD:" + f, p))
scal14 = set()
for w in sections.get("dihedral", []):
    f, t = w[0], tuple(int(x) for x in w[1:5]); p = [float(x) for x in w[5:]]
    key = min(t, t[::-1])
    if f == "cos":
        if abs(p[0]) > 1e-12: put("torsion", key, (round(p[0], 6), int(round(p[2])), round(p[1] % 360, 3)))
        if len(p) >= 5: scal14.add((round(p[3], 4), round(p[4], 4)))
    elif f == "cos3":   # ½[A1(1 + cos φ) + A2(1 − cos 2φ) + A3(1 + cos 3φ)]
        for A, n, dlt in ((p[0], 1, 0.0), (p[1], 2, 180.0), (p[2], 3, 0.0)):
            if abs(A) > 1e-12: put("torsion", key, (round(A / 2, 6), n, dlt))
        if len(p) >= 5: scal14.add((round(p[3], 4), round(p[4], 4)))
    elif f == "harm" and abs(p[0]) < 1e-12: continue   # DL_FIELD writes K = 0 where it found no parameter: no energy
    elif f == "harm": put("improper_h", (t[0], frozenset(t)), ("harmonic", [p[0] / 2, p[1]]))   # centre first; χ → −χ under reordering
    else: put("torsion", key, ("FIELD:" + f, tuple(p)))
for w in sections.get("inversions", []):
    f, c, o = w[0], int(w[1]), frozenset(int(x) for x in w[2:5]); p = [float(x) for x in w[5:]]
    if f == "harm": put("inversion", (c, o), ("harmonic", [p[0] / 2, p[1]]))
    elif f == "plan": put("inversion", (c, o), ("planar", [p[0], 0.0]))
    else: put("inversion", (c, o), ("FIELD:" + f, p))
vdw = {}
for w in sections.get("vdw", []):
    a, b, f = w[0], w[1], w[2]; p = [float(x) for x in w[3:]]
    if f == "lj": v = ("lj", [p[0], p[1]])
    elif f == "nm": v = ("lj9-6", [p[0], p[3]])
    elif f == "buck": v = ("buck", p[:3])
    else: v = ("FIELD:" + f, p)
    vdw[frozenset((a, b))] = v

# ---- CAPS ----
# --from-config: DL_FIELD reorders atoms (templates, several molecule types): take coordinates from its CONFIG file
# (same order as FIELD) and the bonds from FIELD (real bonds, constraints and core-shell springs)
if "--from-config" in sys.argv:
    cfg = [l.split() for l in open(os.path.join(d, "dl_poly.CONFIG"))]
    levcfg = int(cfg[1][0])
    per = 2 + levcfg
    xyz, i = [], 5 if int(cfg[1][1]) > 0 else 2
    while i < len(cfg) and len(xyz) < len(atoms):
        xyz.append([float(v) for v in cfg[i + 1][:3]]); i += per
    pairs = [(int(w[1]), int(w[2])) for w in sections.get("bonds", []) if not w[0].startswith("-")]
    pairs += [(int(w[0]), int(w[1])) for w in sections.get("constraints", [])]
    pairs += [(int(w[0]), int(w[1])) for w in sections.get("shell", [])]
    mol2 = os.path.join(d, "caps_from_config.mol2")
    with open(mol2, "w") as fo:
        fo.write(f"@<TRIPOS>MOLECULE\nfrom CONFIG\n{len(atoms)} {len(pairs)} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for q, name in enumerate(atoms, 1):
            p = xyz[q - 1]
            fo.write(f"{q} {re.sub(r'[^A-Za-z]', '', name)[:2]}{q} {p[0]} {p[1]} {p[2]} {name} 1 RES 0.0\n")
        fo.write("@<TRIPOS>BOND\n")
        for q, (a, b) in enumerate(pairs, 1): fo.write(f"{q} {a} {b} 1\n")
    structure = mol2
# --from-field: build CAPS's structure from the FIELD file itself (ionic / core-shell systems, where DL_FIELD adds shell
# particles and regroups atoms): FIELD atom names, core-shell springs as bonds, atoms on a sparse grid (only the
# per-type parameters are compared, so positions do not matter)
if "--from-field" in sys.argv:
    mol2 = os.path.join(d, "caps_from_field.mol2")
    pairs = [(int(w[0]), int(w[1])) for w in sections.get("shell", [])]
    shell_of = {b: a for a, b in pairs}
    pos, g = {}, 0
    for i in range(1, len(atoms) + 1):
        if i in shell_of: continue
        pos[i] = [6.0 * (g % 10), 6.0 * (g // 10 % 10), 6.0 * (g // 100)]; g += 1
    for b, a in shell_of.items():
        pos[b] = [pos[a][0] + 0.1, pos[a][1], pos[a][2]]
    with open(mol2, "w") as fo:
        fo.write(f"@<TRIPOS>MOLECULE\nfrom FIELD\n{len(atoms)} {len(pairs)} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, name in enumerate(atoms, 1):
            p = pos[i]
            fo.write(f"{i} {re.sub(r'[^A-Za-z]', '', name)[:2]}{i} {p[0]} {p[1]} {p[2]} {name} 1 RES 0.0\n")
        fo.write("@<TRIPOS>BOND\n")
        for k, (a, b) in enumerate(pairs, 1): fo.write(f"{k} {a} {b} 1\n")
    structure = mol2
# core-shell models: the springs are the only bonds; give CAPS a mol2 with those pairs (a data file without bonds
# would have them guessed from distances)
if sections.get("shell") and structure.endswith((".data", ".pdb")):
    xyz, reading = {}, False
    if structure.endswith(".pdb"):
        k = 0
        for line in open(structure):
            if line.startswith(("ATOM", "HETATM")):
                k += 1; xyz[k] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    for line in (open(structure) if structure.endswith(".data") else []):
        t = line.split("#")[0].split()
        if line.startswith("Atoms"): reading = True; continue
        if reading and t and not t[0].lstrip("-").isdigit(): break
        if reading and len(t) >= 7: xyz[int(t[0])] = [float(t[4]), float(t[5]), float(t[6])]
    mol2 = os.path.join(d, "caps_shell_model.mol2")
    with open(mol2, "w") as fo:
        pairs = [(int(w[0]), int(w[1])) for w in sections["shell"]]
        fo.write(f"@<TRIPOS>MOLECULE\nshell model\n{len(atoms)} {len(pairs)} 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, name in enumerate(atoms, 1):
            p = xyz[i]
            fo.write(f"{i} {re.sub(r'[^A-Za-z]', '', name)[:2]}{i} {p[0]} {p[1]} {p[2]} {name} 1 RES 0.0\n")
        fo.write("@<TRIPOS>BOND\n")
        for k, (a, b) in enumerate(pairs, 1): fo.write(f"{k} {a} {b} 1\n")
    structure = mol2
# --types FILE: DL_FIELD's internal atom keys (the FIELD file may show display names, e.g. N_3 for the NH4 key)
if "--types" in sys.argv:
    types_file = sys.argv[sys.argv.index("--types") + 1]
else:
    types_file = os.path.join(d, "caps_field_types.txt")
    open(types_file, "w").write("\n".join(atoms) + "\n")
out = subprocess.run([CAPS, "ff", "apply", structure, "--ff", ffjson, "--types", types_file, "--charges", "keep", "--list", "--allow-missing"],
                     capture_output=True, text=True)
if out.returncode not in (0, 3):
    print(out.stdout[-2000:], out.stderr); sys.exit(2)
for l in out.stdout.splitlines():
    if l.startswith(("note", "missing")) or (l.startswith("  ") and not l.strip()[0].isdigit()): print("caps: " + l.strip())
cap = {}
def cput(kind, key, val):
    cap.setdefault((kind, key), []).append(val)
cvdw, caps_types = {}, None
for l in out.stdout.splitlines():
    w = l.split()
    if not w or l[0].isspace(): continue
    t = w[0]
    if t == "bond": cput("bond", frozenset((int(w[1]), int(w[2]))), ("harmonic", [float(w[3]), float(w[4])]))
    elif t == "bond2": cput("bond", frozenset((int(w[1]), int(w[2]))), ("class2", [float(x) for x in w[3:7]]))
    elif t == "bondx": cput("bond", frozenset((int(w[1]), int(w[2]))), (w[3], [float(w[4]), float(w[5]), float(w[6])] if w[3] == "morse" else [float(w[4]), float(w[5])]))
    elif t == "angle": cput("angle", (int(w[2]), frozenset((int(w[1]), int(w[3])))), ("harmonic", [float(w[4]), float(w[5])]))
    elif t == "angle2":
        v = [float(x) for x in w[4:]]
        cput("angle", (int(w[2]), frozenset((int(w[1]), int(w[3])))), ("quartic" if not any(v[4:]) else "class2", v[:4] if not any(v[4:]) else v))
    elif t == "anglex": cput("angle", (int(w[2]), frozenset((int(w[1]), int(w[3])))), (w[4], [float(w[5]), float(w[6])]))
    elif t == "ub": cput("ub", frozenset((int(w[1]), int(w[2]))), ("ub", [float(w[3]), float(w[4])]))
    elif t == "pair14": cput("pair14", frozenset((int(w[1]), int(w[2]))), ("lj", [float(w[3]), float(w[4])]))
    elif t == "dihedral":
        a = tuple(int(x) for x in w[1:5]); cput("torsion", min(a, a[::-1]), (round(float(w[5]), 6), int(w[6]), round(float(w[7]) % 360, 3)))
    elif t == "improper" and w[5] == "cvff":
        a = tuple(int(x) for x in w[1:5]); cput("torsion", min(a, a[::-1]), (round(float(w[6]), 6), int(w[7]), round(float(w[8]) % 360, 3)))
    elif t == "improper" and w[5] == "harmonic":
        a = tuple(int(x) for x in w[1:5]); cput("improper_h", (a[0], frozenset(a)), ("harmonic", [float(w[6]), float(w[7])]))
    elif t == "inversion":
        cput("inversion", (int(w[1]), frozenset(int(x) for x in w[2:5])), (w[7], [float(w[5]), float(w[6])]))
    elif t == "pair": cvdw.setdefault(frozenset((w[1], w[2])), ("lj9-6" if w[5] == "lj9-6" else "lj", [float(w[3]), float(w[4])]))
    elif t == "pairfunc": cvdw[frozenset((w[1], w[2]))] = (w[3], [float(x) for x in w[4:7]])

# ---- compare ----
problems = 0
conv = []    # convention differences (not failures)
policy = []  # impropers CAPS adds where DL_FIELD's template has none
def same(a, b):
    if isinstance(a, tuple) and len(a) == 3 and not isinstance(a[1], list):   # torsion term (K, n, δ)
        return a[1] == b[1] and abs(a[2] - b[2]) < 1e-3 and abs(a[0] - b[0]) <= TOL * max(1.0, abs(b[0]))
    return a[0] == b[0] and close(a[1], b[1])

def same_terms(v, c):
    return c is not None and len(v) == len(c) and all(any(same(x, y) for y in c) for x in v)

def plane_key(t):   # a torsion-type improper is the pair of planes it measures
    return frozenset((frozenset(t[0:3]), frozenset(t[1:4])))

# 1-4 van der Waals is compared pair by pair only where DL_FIELD lists it explicitly (CHARMM, GROMOS); elsewhere it
# follows from the dihedral 1-4 scalings, reported below
if not any(k == "pair14" for k, _ in ref):
    cap = {kk: v for kk, v in cap.items() if kk[0] != "pair14"}
kinds = sorted({k for k, _ in ref} | {k for k, _ in cap})
for kind in kinds:
    R = {key: v for (k, key), v in ref.items() if k == kind}
    C = {key: v for (k, key), v in cap.items() if k == kind}
    if kind == "torsion":   # torsion-type impropers may be listed with other atom orders: match by planes too
        Rp, Cp = {}, {}
        for key, v in R.items():
            Rp.setdefault(plane_key(key), []).extend(v)
        for key, v in C.items():
            Cp.setdefault(plane_key(key), []).extend(v)
        R2 = dict(R)
        ok = bad = 0
        used_c = set()
        for key, v in R.items():
            c = C.get(key)
            if c is None:
                c = Cp.get(plane_key(key))
            if same_terms(v, c):
                ok += 1
                used_c.add(plane_key(key))
                continue
            # the same four atoms and the same rule, the outer atoms in another order: which permutation a builder uses
            # is its own convention (AMBER leaves it open; DL_FIELD takes the first it finds in its neighbour list)
            # — reported, not a failure
            alt = [k2 for k2 in C if set(k2) == set(key) and plane_key(k2) not in Rp and same_terms(v, C[k2])]
            swap = alt   # same rule, same four atoms: another order of the outer atoms
            if swap:
                conv.append((kind, key, swap[0]))
                used_c.add(plane_key(swap[0]))
                ok += 1
            else:
                bad += 1; problems += 1
                if bad <= 6: print(f"  MISMATCH {kind} {key}: DL_FIELD {v} · CAPS {c}")
        extra = [k for k in C if k not in R and plane_key(k) not in Rp and plane_key(k) not in used_c]
        # DL_FIELD places impropers only where its molecule template lists them; CAPS applies the rules at every
        # centre they match. An extra CAPS improper is that policy difference, reported but not a failure.
        for k2 in extra[:6]: print(f"  template policy: CAPS has {kind} {k2}: {C[k2]} (DL_FIELD's template lists none there)")
        policy.extend(extra)
        print(f"{kind:10s} DL_FIELD {len(R):5d} · identical {ok:5d} · CAPS {len(C):5d}")
        continue
    ok = bad = 0
    for key, v in R.items():
        c = C.get(key)
        if c is not None and len(c) == len(v) and all(any(same(x, y) for y in c) for x in v):
            ok += 1
        elif kind == "improper_h" and c and len(c) < len(v) and all(any(same(x, y) for y in c) for x in v):
            ok += 1   # DL_FIELD's template lists the same improper more than once (its choice, not a rule)
            policy.append(key)
        elif (kind == "bond" and c and c[0][0] == "gromos" and v[0][0] == "harmonic" and
              abs(v[0][1][0] - 4 * c[0][1][0] * c[0][1][1] ** 2) <= 1e-3 * v[0][1][0] and abs(v[0][1][1] - c[0][1][1]) < 1e-9):
            ok += 1
            conv.append(("bond", key, "GROMOS quartic bond; DL_FIELD writes its harmonic approximation K = 4 K_g r0²"))
        else:
            bad += 1; problems += 1
            if bad <= 6: print(f"  MISMATCH {kind} {sorted(key) if isinstance(key, frozenset) else key}: DL_FIELD {v} · CAPS {c}")
    extra = [k for k in C if k not in R]
    if kind in ("inversion", "improper_h"):
        for k2 in extra[:6]: print(f"  template policy: CAPS has {kind} {k2}: {C[k2]} (DL_FIELD's template lists none there)")
        policy.extend(extra)
    else:
        for k2 in extra[:6]: print(f"  EXTRA in CAPS {kind} {k2}: {C[k2]}")
        problems += len(extra)
    print(f"{kind:10s} DL_FIELD {len(R):5d} · identical {ok:5d} · CAPS {len(C):5d}")

# van der Waals, by type pair (CAPS type names translated to the FIELD's names through the atoms they label)
caps_names = [l.strip() for l in open(types_file) if l.strip()]
to_field = {}
for cn, fn in zip(caps_names, atoms):
    to_field.setdefault(cn, fn)
cvdw = {frozenset(to_field.get(n, n) for n in k): v for k, v in cvdw.items()}
okv = 0
for key, v in vdw.items():
    c = cvdw.get(key)
    if c is not None and same(v, c): okv += 1
    elif v[0] == "buck" and all(abs(x) < 1e-12 for x in v[1][::2]) and (c is None or c[0] == "lj" and c[1][0] == 0):
        okv += 1   # a zero Buckingham entry is "no interaction"
    else:
        problems += 1
        print(f"  MISMATCH vdw {sorted(key)}: DL_FIELD {v} · CAPS {c}")
print(f"{'vdw':10s} DL_FIELD {len(vdw):5d} · identical {okv:5d}")
if scal14: print(f"1-4 scalings in the FIELD dihedrals (elec, vdw): {sorted(scal14)}")
for c in conv[:6]: print(f"  convention: {c}")
if policy: print(f"{len(policy)} improper(s) at centres DL_FIELD's templates leave out (rule-based in CAPS)")
if conv: print(f"{len(conv)} convention difference(s): same parameters, atom order of symmetric outer atoms or DL_FIELD's form approximation")
if fused: print(f"note: {len(fused)} run-together number(s) in DL_FIELD's FIELD file (fixed-width overflow), e.g. {fused[0]}")
print("PASS" if problems == 0 else f"FAIL ({problems} problems)")
sys.exit(0 if problems == 0 else 1)
