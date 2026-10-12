"""Torsion parameters a force field lacks, fitted to quantum-chemical rotational profiles.

    python -m caps.torsionfit prepare STRUCTURE[,STRUCTURE...] FF.json WORKDIR [--engine orca|psi4] [--step 15] [--radius 2]
                                    [--level oplsaa-l|dlpno]
    python -m caps.torsionfit collect WORKDIR
    python -m caps.torsionfit fit WORKDIR [--out overlay.json]

prepare  finds every proper torsion of STRUCTURE that FF.json has no parameters for (CAPS's own assignment), groups them
         by the classes the force field looks torsions up with, and for each class cuts a fragment round one example:
         the atoms within RADIUS bonds of its central bond, each cut bond capped by a hydrogen. The fragment is kept only
         when CAPS types the torsion's atoms and their neighbours as in the structure. A scan of the central bond
         follows (360 / STEP points): each point a constrained optimisation with the torsion frozen, then a single-point
         energy, written as inputs for ORCA or Psi4 and a job list (jobs.txt: one shell command per point) to run where
         the program is installed. The levels follow OPLS-AA/L (Kaminski, Friesner, Tirado-Rives & Jorgensen,
         J. Phys. Chem. B 105, 6474 (2001)): geometries HF/6-31G**, energies MP2/cc-pVTZ (RI / density fitted).
         --level dlpno: geometries wB97X-D3/def2-TZVP, energies DLPNO-CCSD(T)/cc-pVTZ with TightPNO (ORCA), for
         checking a force field's torsions against a higher level or refitting them as an overlay.
collect  reads every point's optimised geometry and its MP2 energy.
fit      evaluates the fragment at each optimised geometry with CAPS (the missing torsions at zero, everything else the
         force field's own) and fits, by least squares, the OPLS Fourier terms V1-V4 of every missing class round the
         central bond to E_QM - E_MM (a constant included; points over 10 kcal/mol above the minimum weighted 0.2).
         The result is a force-field overlay (caps ff apply --overlay) whose torsions say where they come from:
         "CAPS fit (QM)", the level, the fragment, the scan and its RMS deviation.

CAPS never fills a torsion without saying so: these are CAPS's fits, not the force field's published parameters.
"""
import json, math, os, re, shutil, subprocess, sys

# --level: the OPLS-AA/L level (the default, as the scans CAPS fitted first), or a higher one for checking or refitting
# (geometries from a range-separated hybrid with dispersion, energies from local coupled cluster). The functional's own
# D3 is part of its name, so nothing of it carries into the single point after $new_job
LEVELS = {
    "oplsaa-l": ("! HF 6-31G** Opt TightSCF", "! RI-MP2 cc-pVTZ cc-pVTZ/C TightSCF",
                 "MP2/cc-pVTZ//HF/6-31G** (as OPLS-AA/L)"),
    "dlpno": ("! wB97X-D3 def2-TZVP def2/J RIJCOSX Opt TightSCF", "! DLPNO-CCSD(T) cc-pVTZ cc-pVTZ/C TightPNO TightSCF",
              "DLPNO-CCSD(T)/cc-pVTZ (TightPNO)//wB97X-D3/def2-TZVP"),
}
ORCA_OPT, ORCA_SP, LEVEL = LEVELS["oplsaa-l"]
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Br": 35, "I": 53}
CAP_BOND = {"C": 1.09, "N": 1.01, "O": 0.96, "S": 1.34, "Si": 1.48, "P": 1.42}


def caps_binary():
    for c in (os.environ.get("CAPS_BIN"), shutil.which("caps"),
              os.path.join(os.path.dirname(__file__), "..", "..", "..", "build", "cli", "caps")):
        if c and os.path.exists(c):
            return os.path.abspath(c)
    raise SystemExit("caps binary not found: set CAPS_BIN")


def arg(k, d=None):
    return sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


# ------------------------------------------------------------------ structures
def read_data(path):
    """atoms (element, type label, position), bonds and dihedral quadruples of a LAMMPS data file CAPS wrote"""
    label, mass, atoms, bonds, dihs, sec = {}, {}, {}, [], set(), None
    for raw in open(path):
        line = raw.split("#")[0].strip()
        if re.match(r"^[A-Z][A-Za-z ]+$", line):
            sec = line
            continue
        if not line:
            continue
        w = line.split()
        if sec == "Masses":
            label[int(w[0])] = raw.split("#")[1].split()[0] if "#" in raw else w[0]
            mass[int(w[0])] = float(w[1])
        elif sec == "Atoms" and len(w) >= 7:
            atoms[int(w[0])] = (int(w[2]), tuple(map(float, w[4:7])))
        elif sec == "Bonds":
            bonds.append((int(w[2]) - 1, int(w[3]) - 1))
        elif sec == "Dihedrals":
            q = tuple(int(x) - 1 for x in w[2:6])
            dihs.add(min(q, q[::-1]))
    ids = sorted(atoms)
    types = [label[atoms[i][0]] for i in ids]
    masses = [mass[atoms[i][0]] for i in ids]
    pos = [atoms[i][1] for i in ids]
    return types, masses, pos, bonds, dihs


def element_of_mass(m):
    best = min(((abs(m - mm), e) for e, mm in (("H", 1.008), ("C", 12.011), ("N", 14.007), ("O", 15.999), ("F", 18.998),
                                               ("Si", 28.086), ("P", 30.974), ("S", 32.06), ("Cl", 35.45), ("Br", 79.904),
                                               ("I", 126.904))))
    return best[1]


def neighbours(n, bonds):
    nb = [[] for _ in range(n)]
    for i, j in bonds:
        nb[i].append(j)
        nb[j].append(i)
    return nb


def propers(nb):
    out = set()
    for b in range(len(nb)):
        for c in nb[b]:
            if b < c:
                for a in nb[b]:
                    for d in nb[c]:
                        if a != c and d != b and a != d:
                            out.add(min((a, b, c, d), (d, c, b, a)))
    return out


def classes(ffjson):
    ff = json.load(open(ffjson))
    return {t["name"]: t.get("equivalence", {}).get("dihedral", t["name"]) for t in ff["atom_types"]}


def write_mol2(path, el, pos, bonds, name="fragment"):
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{name}\n{len(el)} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, (e, p) in enumerate(zip(el, pos)):
            f.write(f"{i + 1} {e}{i + 1} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f} {e} 1 FRG 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for k, (i, j) in enumerate(bonds):
            f.write(f"{k + 1} {i + 1} {j + 1} un\n")


def caps_apply(structure, ffjson, out, overlay=None):
    cmd = [caps_binary(), "ff", "apply", structure, "--ff", ffjson, "--charges", "types", "--allow-missing", "-o", out]
    if overlay:
        cmd += ["--overlay", overlay]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError((r.stderr or r.stdout).strip().splitlines()[-1])
    m = re.search(r"total (\S+)", r.stdout)
    return float(m.group(1)) if m else None


def listed_missing(structure, ffjson, cl):
    """the torsion classes CAPS reports it has no parameters for. A data file leaves out torsions whose terms are all
    zero (HC-CT-C-O in OPLS-AA ...), so a quadruple missing from it is not by itself a missing torsion"""
    tmp = structure + ".missing.data"
    r = subprocess.run([caps_binary(), "ff", "apply", structure, "--ff", ffjson, "--charges", "types", "-o", tmp],
                       capture_output=True, text=True)
    if os.path.exists(tmp):
        os.remove(tmp)
    out = set()
    for line in (r.stdout + r.stderr).splitlines():
        m = re.match(r"^\s+dihedral\s+(.+)$", line)
        if m:
            key = tuple(cl.get(t, t) for t in re.findall(r"(\S+?)\(", m.group(1)))
            out.add(min(key, key[::-1]))
    return out


# ------------------------------------------------------------------ geometry
def sub(a, b): return [a[k] - b[k] for k in range(3)]
def dot(a, b): return sum(a[k] * b[k] for k in range(3))
def cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
def norm(a): return math.sqrt(dot(a, a))


def dihedral(p, a, b, c, d):
    b0, b1, b2 = sub(p[a], p[b]), sub(p[c], p[b]), sub(p[d], p[c])
    n1 = norm(b1)
    v = [b0[k] - dot(b0, b1) / n1 ** 2 * b1[k] for k in range(3)]
    w = [b2[k] - dot(b2, b1) / n1 ** 2 * b1[k] for k in range(3)]
    x, y = dot(v, w), dot(cross([b1[k] / n1 for k in range(3)], v), w)
    return math.degrees(math.atan2(y, x))


def rotate(p, moving, b, c, angle):
    """positions with the atoms `moving` turned by `angle` degrees about the axis b -> c"""
    axis = sub(p[c], p[b])
    n = norm(axis)
    u = [x / n for x in axis]
    t = math.radians(angle)
    ct, st = math.cos(t), math.sin(t)
    out = [list(x) for x in p]
    for i in moving:
        v = sub(p[i], p[c])
        r = [v[k] * ct + cross(u, v)[k] * st + u[k] * dot(u, v) * (1 - ct) for k in range(3)]
        out[i] = [p[c][k] + r[k] for k in range(3)]
    return out


def in_ring(nb, b, c):
    """is the bond b-c in a ring (c reachable from b without that bond)?"""
    seen, todo = {b}, [b]
    while todo:
        x = todo.pop()
        for y in nb[x]:
            if (x, y) in ((b, c), (c, b)):
                continue
            if y == c:
                return True
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return False


def side(nb, start, block):
    seen, todo = {start}, [start]
    while todo:
        x = todo.pop()
        for y in nb[x]:
            if y != block and y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


# ------------------------------------------------------------------ prepare
def fragment(el, pos, nb, b, c, radius):
    """atoms within `radius` bonds of b or c; cut bonds capped by H along the bond"""
    dist = {b: 0, c: 0}
    todo = [b, c]
    while todo:
        x = todo.pop(0)
        for y in nb[x]:
            if y not in dist and dist[x] < radius:
                dist[y] = dist[x] + 1
                todo.append(y)
    keep = sorted(dist)
    idx = {a: k for k, a in enumerate(keep)}
    fel, fpos, fbonds = [el[a] for a in keep], [list(pos[a]) for a in keep], []
    for a in keep:
        for y in nb[a]:
            if y in idx and idx[a] < idx[y]:
                fbonds.append((idx[a], idx[y]))
            elif y not in idx:   # a cut bond: H on the kept atom, along the bond
                v = sub(pos[y], pos[a])
                L = CAP_BOND.get(el[a], 1.09)
                fel.append("H")
                fpos.append([pos[a][k] + v[k] / norm(v) * L for k in range(3)])
                fbonds.append((idx[a], len(fel) - 1))
    return keep, idx, fel, fpos, fbonds


def orca_input(el, pos, quad, charge=0, mult=1, level="oplsaa-l"):
    a, b, c, d = quad
    opt, sp, _ = LEVELS[level]
    xyz = "\n".join(f"{e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}" for e, p in zip(el, pos))
    mem = "" if level == "oplsaa-l" else "%maxcore 6000\n"   # local coupled cluster wants more memory per core
    return (f"{opt}\n{mem}%geom Constraints {{ D {a} {b} {c} {d} C }} end end\n* xyz {charge} {mult}\n{xyz}\n*\n\n"
            f"$new_job\n{sp}\n{mem}* xyzfile {charge} {mult} point.xyz\n")


def psi4_input(el, pos, quad, charge=0, mult=1):
    a, b, c, d = (x + 1 for x in quad)
    xyz = "\n".join(f"  {e} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}" for e, p in zip(el, pos))
    return (f"memory 4 GB\nmolecule {{\n{charge} {mult}\n{xyz}\nno_com\nno_reorient\n}}\n"
            f"set {{ basis 6-31G**\n  scf_type pk\n  optking {{ frozen_dihedral = (\"{a} {b} {c} {d}\") }} }}\n"
            f"optimize('scf')\nset {{ basis cc-pvtz\n  scf_type df\n  mp2_type df }}\n"
            f"e = energy('mp2')\nwith open('point.result', 'w') as f:\n"
            f"    f.write(repr(e) + '\\n')\n    f.write(psi4.core.get_active_molecule().save_string_xyz())\n")


def prepare():
    structures, ffjson, work = sys.argv[2].split(","), os.path.abspath(sys.argv[3]), os.path.abspath(sys.argv[4])
    engine, step, radius = arg("--engine", "orca"), int(arg("--step", "15")), int(arg("--radius", "1"))
    level = arg("--level", "oplsaa-l")
    if level not in LEVELS:
        raise SystemExit(f"--level is one of {', '.join(LEVELS)}")
    if level != "oplsaa-l" and engine != "orca":
        raise SystemExit("--level dlpno needs --engine orca")
    os.makedirs(work, exist_ok=True)
    cl = classes(ffjson)
    missing, mols = {}, []
    for n, structure in enumerate(structures):   # each missing class from the first structure that has it
        data = os.path.join(work, f"structure{n}.data")
        caps_apply(structure, ffjson, data)
        types, masses, pos, bonds, have = read_data(data)
        el = [element_of_mass(m) for m in masses]
        nb = neighbours(len(el), bonds)
        mols.append((types, pos, el, nb))
        listed = listed_missing(structure, ffjson, cl)
        for q in propers(nb) - have:
            key = tuple(cl.get(types[x], types[x]) for x in q)
            key = min(key, key[::-1])
            if key not in listed:   # an all-zero torsion of the force field, left out of the data file
                continue
            missing.setdefault(key, (n, []))[1].append(q) if missing.get(key, (n,))[0] == n else None
    plan = {"ff": ffjson, "engine": engine, "level": LEVELS[level][2], "level_key": level, "step": step, "scans": [],
            "not_scanned": []}
    jobs = []
    for key, (n, qs) in sorted(missing.items()):
        types, pos, el, nb = mols[n]
        a, b, c, d = qs[0]
        if in_ring(nb, b, c):   # the central bond in a ring: it cannot turn freely
            plan["not_scanned"].append({"classes": key, "why": "the central bond is in a ring"})
            continue
        for rad in range(radius, radius + 3):   # the smallest fragment that keeps the types
            keep, idx, fel, fpos, fbonds = fragment(el, pos, nb, b, c, rad)
            fdir = os.path.join(work, "-".join(key).replace("*", "x"))
            os.makedirs(fdir, exist_ok=True)
            m2 = os.path.join(fdir, "fragment.mol2")
            write_mol2(m2, fel, fpos, fbonds)
            fdata = os.path.join(fdir, "fragment.data")
            try:
                caps_apply(m2, ffjson, fdata)
            except RuntimeError:   # the cut left atoms CAPS cannot type: a larger fragment
                continue
            ftypes = read_data(fdata)[0]
            # the torsion's atoms and their neighbours keep the structure's types
            check = {x for x in (a, b, c, d)} | set(nb[b]) | set(nb[c])
            if all(x in idx and ftypes[idx[x]] == types[x] for x in check):
                break
        else:
            plan["not_scanned"].append({"classes": key, "why": "no fragment keeps the structure's types"})
            continue
        quad = [idx[a], idx[b], idx[c], idx[d]]
        fnb = neighbours(len(fel), fbonds)
        moving = side(fnb, quad[2], quad[1])
        if len(moving) > len(fel) / 2:   # turn the smaller half (the same dihedral change, fewer clashes)
            moving = side(fnb, quad[1], quad[2])
            sign = -1
        else:
            sign = 1
        moving = sorted(moving)
        phi0 = dihedral(fpos, *quad)
        skipped = 0
        scan = {"classes": key, "dir": os.path.relpath(fdir, work), "quad": quad, "elements": fel, "bonds": fbonds, "points": []}
        for k in range(0, 360, step):
            target = -180 + k
            p = rotate(fpos, moving, quad[1], quad[2], sign * (target - phi0))
            near = [set(fnb[i]) | {z for y in fnb[i] for z in fnb[y]} for i in range(len(p))]   # 1-2 and 1-3 partners
            close = min(norm(sub(p[i], p[j])) for i in range(len(p)) for j in range(i) if j not in near[i])
            if close < 1.0:   # atoms on top of each other: no point, no input, no job (the optimiser would start broken)
                skipped += 1
                continue
            pdir = os.path.join(fdir, f"p{target:+04d}")
            os.makedirs(pdir, exist_ok=True)
            rel = os.path.relpath(pdir, work)   # jobs.txt and plan.json name points relative to WORKDIR: copy it anywhere
            if engine == "orca":
                open(os.path.join(pdir, "point.inp"), "w").write(orca_input(fel, p, quad, level=level))
                jobs.append(f"cd {rel} && \"${{ORCA:-orca}}\" point.inp > point.out")
            else:
                open(os.path.join(pdir, "point.in"), "w").write(psi4_input(fel, p, quad))
                jobs.append(f"cd {rel} && ${{PSI4:-psi4}} point.in point.out")
            scan["points"].append({"phi": target, "dir": rel, "closest_nonbonded": round(close, 3)})
        scan["skipped_points"] = skipped
        plan["scans"].append(scan)
    json.dump(plan, open(os.path.join(work, "plan.json"), "w"), indent=1)
    # one line per planned point, each ended by a newline: the array is 1-N with N the number of lines
    assert len(jobs) == sum(len(s["points"]) for s in plan["scans"])
    open(os.path.join(work, "jobs.txt"), "w").write("".join(j + "\n" for j in jobs))
    print(f"{len(missing)} missing torsion classes: {len(plan['scans'])} scans, {len(jobs)} QM points (jobs.txt; array 1-{len(jobs)}); "
          f"not scanned: {len(plan['not_scanned'])}")
    for s in plan["scans"]:
        print(f"  {'-'.join(s['classes'])}: fragment of {len(s['elements'])} atoms, closest non-bonded pair over the scan "
              f"{min(p['closest_nonbonded'] for p in s['points']):.2f} A" + (f", {s['skipped_points']} points clash (not run)" if s["skipped_points"] else ""))
    for s in plan["not_scanned"]:
        print(f"  {'-'.join(s['classes'])}: not scanned ({s['why']})")


# ------------------------------------------------------------------ collect, fit
def read_point(pdir, engine):
    """(energy in hartree, element list, positions) of one finished point, or None"""
    if engine == "orca":
        out, xyz = os.path.join(pdir, "point.out"), os.path.join(pdir, "point.xyz")
        if not (os.path.exists(out) and os.path.exists(xyz)):
            return None
        text = open(out, errors="replace").read()
        # ORCA also ends normally after an optimisation that ran out of cycles: such a point is not used, nor one
        # whose correlated single point is missing
        if ("ORCA TERMINATED NORMALLY" not in text or "OPTIMIZATION RUN DONE" not in text
                or "did not converge" in text or not re.search(r"MP2 CORRELATION ENERGY|MP2 TOTAL ENERGY|E\(CCSD\(T\)\)", text)):
            return None
        e = float(re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)[-1])
        lines = open(xyz).read().splitlines()[2:]
    else:
        res = os.path.join(pdir, "point.result")
        if not os.path.exists(res):
            return None
        lines = open(res).read().splitlines()
        e = float(lines[0])
        lines = lines[2:]
    w = [l.split() for l in lines if l.strip()]
    return e, [x[0] for x in w], [tuple(map(float, x[1:4])) for x in w]


def collect():
    work = os.path.abspath(sys.argv[2])
    plan = json.load(open(os.path.join(work, "plan.json")))
    done = total = 0
    for s in plan["scans"]:
        for p in s["points"]:
            total += 1
            r = read_point(os.path.join(work, p["dir"]), plan["engine"])
            if r:
                p["E_qm"], p["positions"] = r[0], r[2]
                done += 1
    json.dump(plan, open(os.path.join(work, "plan.json"), "w"), indent=1)
    print(f"{done} of {total} QM points finished")


def lstsq(A, y, w):
    """weighted least squares by the normal equations (a tiny ridge keeps an unused column at zero)"""
    n = len(A[0])
    M = [[sum(w[r] ** 2 * A[r][i] * A[r][j] for r in range(len(A))) + (1e-9 if i == j else 0.0) for j in range(n)] for i in range(n)]
    b = [sum(w[r] ** 2 * A[r][i] * y[r] for r in range(len(A))) for i in range(n)]
    for c in range(n):   # Gauss-Jordan with partial pivoting
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p], b[c], b[p] = M[p], M[c], b[p], b[c]
        for r in range(n):
            if r != c and M[c][c]:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(n)]
                b[r] -= f * b[c]
    return [b[i] / M[i][i] if M[i][i] else 0.0 for i in range(n)]


def fit():
    work = os.path.abspath(sys.argv[2])
    plan = json.load(open(os.path.join(work, "plan.json")))
    ffjson, cl = plan["ff"], classes(plan["ff"])
    rules = []
    for s in plan["scans"]:
        pts = [p for p in s["points"] if "E_qm" in p]
        if len(pts) < 8:
            print(f"{'-'.join(s['classes'])}: {len(pts)} points finished, not fitted")
            continue
        fnb = neighbours(len(s["elements"]), [tuple(b) for b in s["bonds"]])
        b, c = s["quad"][1], s["quad"][2]
        rows, emm, eqm = [], [], []
        groups = None
        for p in pts:
            m2 = os.path.join(work, p["dir"], "opt.mol2")
            write_mol2(m2, s["elements"], p["positions"], [tuple(x) for x in s["bonds"]])
            data = os.path.join(work, p["dir"], "opt.data")
            emm.append(caps_apply(m2, ffjson, data))
            eqm.append(p["E_qm"] * 627.509474)
            types, _, pos, _, have = read_data(data)
            if groups is None:   # the missing torsions round the central bond, by class
                groups = {}
                listed = listed_missing(m2, ffjson, cl)
                for q in propers(fnb) - have:
                    if {q[1], q[2]} == {b, c}:
                        key = tuple(cl.get(types[x], types[x]) for x in q)
                        if min(key, key[::-1]) not in listed:
                            continue
                        groups.setdefault(min(key, key[::-1]), []).append(q)
                keys = sorted(groups)
            row = [1.0]
            for k in keys:
                phis = [math.radians(dihedral(pos, *q)) for q in groups[k]]
                row += [sum(0.5 * (1 + math.cos(f)) for f in phis), sum(0.5 * (1 - math.cos(2 * f)) for f in phis),
                        sum(0.5 * (1 + math.cos(3 * f)) for f in phis), sum(0.5 * (1 - math.cos(4 * f)) for f in phis)]
            rows.append(row)
        y = [q - m for q, m in zip(eqm, emm)]
        wts = [0.2 if q - min(eqm) > 10 else 1.0 for q in eqm]
        sol = lstsq(rows, y, wts)
        fitted = [sum(a * x for a, x in zip(r, sol)) for r in rows]
        rms = math.sqrt(sum((f - v) ** 2 for f, v in zip(fitted, y)) / len(y))
        # the profiles themselves, for the record: QM, MM without the fitted torsions, MM with them
        scan_rows = [(p["phi"], q - min(eqm), m, f) for p, q, m, f in zip(pts, eqm, emm, fitted)]
        for n, k in enumerate(keys):
            V = [round(float(v), 4) for v in sol[1 + 4 * n: 5 + 4 * n]]
            rules.append({"name": "-".join(k), "match": [x.replace("*", "\\*") for x in k], "style": "opls", "params": V,
                          "comment": f"CAPS fit (QM): {plan['level']}, {len(pts)}-point scan of a {len(s['elements'])}-atom "
                                     f"fragment, RMS {rms:.2f} kcal/mol",
                          "fit": {"fragment": s["dir"], "scan": [{"phi": a, "E_qm": round(q, 4), "E_mm_without": round(m, 4),
                                                                  "E_fit_correction": round(f, 4)} for a, q, m, f in scan_rows]}})
            print(f"{'-'.join(k)}: V = {V} (RMS {rms:.2f} kcal/mol over {len(pts)} points)")
    out = arg("--out", os.path.join(work, "fitted_torsions.json"))
    json.dump({"format": "caps-forcefield", "format_version": 1, "name": "CAPS-fitted torsions", "atom_types": [],
               "dihedrals": rules, "notes": [f"torsions CAPS fitted to {plan['level']} rotational profiles "
                                             "(python -m caps.torsionfit); not parameters of the force field's authors"]},
              open(out, "w"), indent=1)
    print(f"wrote {out}: {len(rules)} torsions (caps ff apply --overlay {out})")


if __name__ == "__main__":
    {"prepare": prepare, "collect": collect, "fit": fit}[sys.argv[1]]()
