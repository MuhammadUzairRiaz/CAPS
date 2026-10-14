"""Automatic parameters for the bonded terms a force field lacks, from AIMNet2 instead of a quantum-chemistry run.

    python -m caps.autofit STRUCTURE --ff FF.json --out DIR [--model aimnet2] [--step 15] [--confirm-qm]

STRUCTURE is anything CAPS types (mol2, data file, a React product). CAPS assigns the force field and lists what has no
parameters; every missing term is then fitted, in order, so each stage sees the terms fitted before it (a torsion fitted
while a bend is still missing takes up that bend's error):

  bonds      a 7-point scan of the bond (r_e +- 0.06 A), each point relaxed by AIMNet2 with the bond held; E = K (r - r0)^2
  bends      a 9-point scan of the angle (theta_e +- 12 deg), relaxed with the angle held; E = K (theta - theta0)^2
  torsions   caps.torsionfit's own fragments, scan points, fill-in of clashing points and joint Fourier fit (OPLS V1-V4), the
             QM energies replaced by AIMNet2's, each point relaxed with the dihedral (and any held neighbours) fixed

Every term is fitted against the MM energy of the same geometries with that term absent (its classes at zero), so the
fit holds exactly what the force field lacks. The result is DIR/autofit.json, a caps-forcefield overlay
(caps ff apply ... --overlay DIR/autofit.json), each entry marked "CAPS auto-fit (AIMNet2 ...)" with its scan and RMS,
and DIR/report.json.

How far to trust it: against CAPS's MP2 and DLPNO-CCSD(T) torsion scans (bench/ff/aimnet_torsion_benchmark.py), AIMNet2
reproduces carbon-backbone torsions to 0.1-0.7 kcal/mol, but not torsions about a C-O bond whose carbon carries a second
O or N (acetals, anomeric centres: up to 2.5 kcal/mol) nor an ester's C(=O)-O rotation (its barrier ~10 % low). Those
are fitted all the same and flagged "needs QM"; --confirm-qm writes their DLPNO-CCSD(T) scan set (caps.torsionfit
prepare --level dlpno) into DIR/qm_confirm for the cluster.

Charges are not fitted: a group CAPS cannot make neutral is reported, not filled.

Needs aimnet[ase] (MIT, isayevlab/aimnetcentral) and numpy in the Python that runs it; CAPS's own caps binary for the
MM energies (CAPS_BIN, else the build tree).
"""
import json, math, os, re, shutil, subprocess, sys

from . import torsionfit as T

HARTREE_EV = 27.211386245988
EV_KCAL = 23.0605478


def arg(k, d=None):
    return sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d


# ------------------------------------------------------------------ AIMNet2
_CALC = None


def calculator(model):
    global _CALC
    if _CALC is None:
        from aimnet.calculators import AIMNet2ASE
        _CALC = AIMNet2ASE(model)
    return _CALC


def relax(el, pos, model, dihedrals=(), angles=(), bonds=(), fmax=0.01, steps=400):
    """AIMNet2 energy (kcal/mol) and positions after a BFGS optimisation with the given internal coordinates held at
    their present values (index tuples)"""
    from ase import Atoms
    from ase.constraints import FixInternals
    from ase.optimize import BFGS
    a = Atoms(symbols=el, positions=pos)
    a.calc = calculator(model)
    cons = {}
    if dihedrals:
        cons["dihedrals_deg"] = [[a.get_dihedral(*q), list(q)] for q in dihedrals]
    if angles:
        cons["angles_deg"] = [[a.get_angle(*q), list(q)] for q in angles]
    if bonds:
        cons["bonds"] = [[a.get_distance(*q), list(q)] for q in bonds]
    if cons:
        a.set_constraint(FixInternals(**cons))
    BFGS(a, logfile=None).run(fmax=fmax, steps=steps)
    return a.get_potential_energy() * EV_KCAL, [tuple(map(float, p)) for p in a.get_positions()]


# ------------------------------------------------------------------ what is missing
def missing_terms(structure, ffjson):
    """{"bond": [names...], "angle": [...], "dihedral": [...]} from caps ff apply, each a list of type names"""
    tmp = structure + ".autofit.data"
    r = subprocess.run([T.caps_binary(), "ff", "apply", structure, "--ff", ffjson, "--charges", "types", "-o", tmp],
                       capture_output=True, text=True)
    if os.path.exists(tmp):
        os.remove(tmp)
    text = r.stdout + r.stderr
    out = {"bond": [], "angle": [], "dihedral": [], "untyped": [], "charge": []}
    # a failure that is not a report of missing terms (typing rules not found, unreadable structure) must never read as
    # "nothing missing"
    if r.returncode and not re.search(r"no parameters for|match no typing rule|charges add up", text):
        raise RuntimeError("caps ff apply failed: " + (text.strip().splitlines() or ["(no output)"])[-1])
    for line in text.splitlines():
        m = re.match(r"^\s+(bond|angle|dihedral)\s+(.+)$", line)
        if m:
            names = re.findall(r"(\S+?)\(", m.group(2)) or m.group(2).split()
            if names not in out[m.group(1)] and names[::-1] not in out[m.group(1)]:
                out[m.group(1)].append(names)
        if "match no typing rule" in line:
            out["untyped"].append(line.strip())
        if "charges add up" in line:
            out["charge"].append(line.strip())
    return out


def needs_qm(el, nb, quad):
    """the chemistries AIMNet2 does not reproduce (see the module text): a C-O central bond whose carbon has a second O
    or N, or an ester C(=O)-O bond"""
    b, c = quad[1], quad[2]
    for x, y in ((b, c), (c, b)):
        if el[x] == "C" and el[y] == "O":
            others = [z for z in nb[x] if z != y and el[z] in ("O", "N")]
            if others:
                return "a C-O bond whose carbon carries another O / N (acetal, ester, carbamate)"
    return None


# ------------------------------------------------------------------ a fragment that keeps the types
def fragment_around(structure, ffjson, work, atoms, label):
    """the smallest fragment round the bond atoms[0]-atoms[1] that keeps the types of `atoms` and their neighbours:
    elements, positions, bonds and the index map"""
    data = os.path.join(work, label + ".structure.data")
    T.caps_apply(structure, ffjson, data)
    types, masses, pos, bonds, _ = T.read_data(data)
    box = T.read_box(data)
    el = [T.element_of_mass(m) for m in masses]
    nb = T.neighbours(len(el), bonds)
    check = set(atoms) | {y for x in atoms for y in nb[x]}
    for rad in range(1, max(4, len(el))):
        keep, idx, fel, fpos, fbonds = T.fragment(el, pos, nb, atoms[0], atoms[1], rad, box)
        m2 = os.path.join(work, label + ".fragment.mol2")
        T.write_mol2(m2, fel, fpos, fbonds)
        try:
            T.caps_apply(m2, ffjson, os.path.join(work, label + ".fragment.data"))
        except RuntimeError:
            continue
        ft = T.read_data(os.path.join(work, label + ".fragment.data"))[0]
        if all(x in idx and ft[idx[x]] == types[x] for x in check):
            return fel, fpos, fbonds, idx, types
    raise RuntimeError(f"{label}: no fragment keeps the types")


def find_occurrence(structure, ffjson, work, names, kind):
    """atom indices of one bond / angle with these type names in the structure"""
    data = os.path.join(work, "occ.data")
    T.caps_apply(structure, ffjson, data)
    types, masses, pos, bonds, _ = T.read_data(data)
    nb = T.neighbours(len(types), bonds)
    want = [names, names[::-1]]
    if kind == "bond":
        for i, j in bonds:
            if [types[i], types[j]] in want:
                return (i, j)
    else:
        for j in range(len(types)):
            for i in nb[j]:
                for k in nb[j]:
                    if i < k and [types[i], types[j], types[k]] in want:
                        return (i, j, k)
    raise RuntimeError(f"{kind} {'-'.join(names)} not found")


# ------------------------------------------------------------------ bonds and bends
def class_of(ffjson, kind):
    d = json.load(open(ffjson))
    return {t["name"]: (t.get("equivalence") or {}).get(kind, t["name"]) for t in d["atom_types"]}


def fit_stretch_or_bend(structure, ffjson, work, names, kind, model):
    """K and the equilibrium value of a bond (kcal/mol/A^2, A) or bend (kcal/mol/rad^2, deg) by an AIMNet2 scan"""
    label = f"{kind}_{'-'.join(names)}"
    occ = find_occurrence(structure, ffjson, work, names, kind)
    fel, fpos, fbonds, idx, _ = fragment_around(structure, ffjson, work, list(occ), label)
    q = [idx[x] for x in occ]
    e0, p0 = relax(fel, fpos, model)
    import numpy as np
    from ase import Atoms
    a0 = Atoms(symbols=fel, positions=p0)
    # every occurrence of the class in the fragment (a class met twice is fitted as the sum of its terms)
    cl = class_of(ffjson, kind)
    fdata = os.path.join(work, label + ".fragment.data")
    ftypes = T.read_data(fdata)[0]
    fnb = T.neighbours(len(fel), fbonds)
    key = [cl.get(n, n) for n in names]
    def same(tup):
        k = [cl.get(ftypes[x], ftypes[x]) for x in tup]
        return k == key or k == key[::-1]
    if kind == "bond":
        occs = [tuple(b) for b in fbonds if same(tuple(b))]
        x0 = a0.get_distance(*q)
        grid = [x0 + d for d in (-0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06)]
    else:
        occs = [(i, j, k) for j in range(len(fel)) for i in fnb[j] for k in fnb[j] if i < k and same((i, j, k))]
        x0 = a0.get_angle(*q)
        grid = [x0 + d for d in range(-12, 13, 3)]
    rows, eai, emm, done = [], [], [], []
    for target in grid:
        a = a0.copy()
        try:
            if kind == "bond":
                a.set_distance(q[0], q[1], target, fix=0.5)
                e, p = relax(fel, a.get_positions(), model, bonds=[tuple(q)])
            else:
                a.set_angle(q[0], q[1], q[2], target)
                e, p = relax(fel, a.get_positions(), model, angles=[tuple(q)])
        except Exception:   # the constraint solver gave up on this point: left out
            continue
        m2 = os.path.join(work, f"{label}_{len(done)}.mol2")
        T.write_mol2(m2, fel, p, [tuple(b) for b in fbonds])
        emm.append(T.caps_apply(m2, ffjson, os.path.join(work, f"{label}_{len(done)}.data")))
        eai.append(e)
        b = Atoms(symbols=fel, positions=p)
        vals = [b.get_distance(*o) if kind == "bond" else math.radians(b.get_angle(*o)) for o in occs]
        rows.append([sum(v * v for v in vals), sum(vals), 1.0])
        done.append({"value": round(b.get_distance(*q) if kind == "bond" else b.get_angle(*q), 4)})
    y = np.array(eai) - np.array(emm)
    A = np.array(rows)
    # E = K sum (x - x0)^2 = K sum x^2 - 2 K x0 sum x + const: K and x0 from the first two coefficients
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    K = sol[0]
    x0f = -sol[1] / (2 * K) if K != 0 else float("nan")
    fit = A @ sol
    rms = float(np.sqrt(np.mean((fit - y) ** 2)))
    if kind == "angle":
        x0f = math.degrees(x0f)
    for d, a_, m, f in zip(done, eai, emm, fit):
        d.update({"E_ai": round(a_ - min(eai), 4), "E_mm_without": round(m - min(emm), 4)})
    return {"name": "-".join(key), "match": key, "style": "harmonic",
            "params": [round(float(K), 4), round(float(x0f), 4 if kind == "bond" else 3)],
            "comment": f"CAPS auto-fit (AIMNet2 {MODEL_LABEL}): {len(grid)}-point relaxed scan of the {kind} in a "
                       f"{len(fel)}-atom fragment, RMS {rms:.2f} kcal/mol",
            "fit": {"fragment_atoms": len(fel), "occurrences": len(occs), "scan": done, "rms": round(rms, 3)}}, rms


# ------------------------------------------------------------------ torsions through caps.torsionfit
def parse_inp(path):
    """elements and positions of an ORCA input's * xyz block"""
    lines = open(path).read().splitlines()
    i = next(n for n, l in enumerate(lines) if l.startswith("* xyz"))
    el, pos = [], []
    for l in lines[i + 1:]:
        if l.strip() == "*":
            break
        w = l.split()
        el.append(w[0])
        pos.append(tuple(map(float, w[1:4])))
    return el, pos


def held(plan_scan):
    return [tuple(h) for h in plan_scan.get("hold", [])]


def run_points(work, model, only_new=False):
    plan_path = os.path.join(work, "plan.json")
    plan = json.load(open(plan_path))
    n = 0
    for s in plan["scans"]:
        for p in s["points"]:
            if "E_qm" in p or "failed" in p:
                continue
            el, pos = parse_inp(os.path.join(work, p["dir"], "point.inp"))
            try:
                e, newpos = relax(el, pos, model, dihedrals=[tuple(s["quad"])] + held(s))
            except Exception as err:   # ASE's constraint solver gave up on this start: the point is left out, as a QM point that failed
                p["failed"] = str(err)[:200]
                continue
            p["E_qm"] = e / EV_KCAL / HARTREE_EV   # hartree, as torsionfit reads a QM point
            p["positions"] = newpos
            n += 1
    plan["level"] = f"AIMNet2 {MODEL_LABEL} (relaxed, dihedral held)"
    json.dump(plan, open(plan_path, "w"), indent=1)
    return n


def torsionfit(*args):
    r = subprocess.run([sys.executable, "-m", "caps.torsionfit", *args], capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": os.path.dirname(os.path.dirname(os.path.abspath(__file__)))})
    if r.returncode:
        raise RuntimeError(f"torsionfit {args[0]}: {(r.stderr or r.stdout).strip()[-400:]}")
    return r.stdout


def fit_torsions(structure, ffjson, work, model, step, hold):
    tw = os.path.join(work, "torsions")
    if os.path.exists(tw):
        shutil.rmtree(tw)
    torsionfit("prepare", structure, ffjson, tw, "--step", str(step), *(["--hold-neighbours"] if hold else []))
    plan = json.load(open(os.path.join(tw, "plan.json")))
    if not plan["scans"]:
        return [], plan.get("not_scanned", []), ""

    run_points(tw, model)
    for _ in range(3):   # clashing angles: rebuilt from the nearest finished geometry, as torsionfit fill does for QM
        out = torsionfit("fill", tw)
        if out.startswith("0 new"):
            break
        run_points(tw, model)
    out = os.path.join(tw, "fitted.json")
    log = torsionfit("fit", tw, "--joint", "--out", out)
    rules = json.load(open(out))["dihedrals"]
    # flags: the chemistries AIMNet2 does not reproduce
    plan = json.load(open(os.path.join(tw, "plan.json")))
    flags = {}
    for s in plan["scans"]:
        el = s["elements"]
        nb = T.neighbours(len(el), [tuple(b) for b in s["bonds"]])
        why = needs_qm(el, nb, s["quad"])
        if why:
            flags["-".join(s["classes"])] = why
    for r in rules:
        r["comment"] = r["comment"].replace("CAPS fit (QM)", "CAPS auto-fit (AIMNet2)").replace("CAPS refit (QM)", "CAPS auto-refit (AIMNet2)")
        key = "-".join(x.replace("\\*", "*") for x in r["match"])
        why = flags.get(key) or flags.get("-".join(key.split("-")[::-1]))
        if why:
            r["comment"] += f"; needs QM: {why}"
            r["needs_qm"] = True
    return rules, plan.get("not_scanned", []), log


# ------------------------------------------------------------------ the whole run
MODEL_LABEL = "wB97M-D3"


def combined_ff(ffjson, extra, path):
    d = json.load(open(ffjson))
    # a library names its typing rules relative to itself (../typing/...): the copy keeps the absolute path
    if isinstance(d.get("typing"), str) and not os.path.isabs(d["typing"]):
        d["typing"] = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(ffjson)), d["typing"]))
    for kind in ("bonds", "angles", "dihedrals"):
        d.setdefault(kind, [])
        d[kind] = [r for r in d[kind] if r.get("match") not in [e["match"] for e in extra.get(kind, [])]] + extra.get(kind, [])
    json.dump(d, open(path, "w"))
    return path


def main():
    structure, ffjson = os.path.abspath(sys.argv[1]), os.path.abspath(arg("--ff"))
    work = os.path.abspath(arg("--out", "autofit"))
    model = arg("--model", "aimnet2")
    step = int(arg("--step", "15"))
    os.makedirs(work, exist_ok=True)
    report = {"structure": structure, "forcefield": ffjson, "model": model, "terms": {}, "not_fitted": []}
    miss = missing_terms(structure, ffjson)
    report["missing_before"] = {k: ["-".join(n) for n in v] if k in ("bond", "angle", "dihedral") else v for k, v in miss.items()}
    if miss["untyped"]:
        raise SystemExit("untyped atoms: give the types first (caps ff type --explain); " + "; ".join(miss["untyped"]))
    chk = os.path.join(work, "geometry.data")
    T.caps_apply(structure, ffjson, chk)
    _, masses, pos, bonds, _ = T.read_data(chk)
    box = T.read_box(chk)
    el = [T.element_of_mass(m) for m in masses]
    blen = lambda a, b: T.norm(T.min_image(T.sub(pos[a], pos[b]), box))   # periodic: a network loops through the cell
    long = [(a, b) for a, b in bonds if blen(a, b) > 1.25 * (T.COV.get(el[a], 0.8) + T.COV.get(el[b], 0.8))]
    if long:
        raise SystemExit(f"{len(long)} bond(s) far longer than covalent (e.g. atoms {long[0][0] + 1}-{long[0][1] + 1}, "
                         f"{blen(*long[0]):.2f} A): relax the structure first (an unrelaxed React product)")
    extra = {"bonds": [], "angles": [], "dihedrals": []}
    current = ffjson
    for kind, key in (("bond", "bonds"), ("angle", "angles")):
        for names in miss[kind]:
            try:
                rule, rms = fit_stretch_or_bend(structure, current, work, names, kind, model)
            except Exception as e:   # reported, never filled with a guess
                report["not_fitted"].append({"term": f"{kind} {'-'.join(names)}", "why": str(e)[:300]})
                continue
            extra[key].append(rule)
            print(f"{kind:8s} {rule['name']:24s} K {rule['params'][0]:9.3f}  x0 {rule['params'][1]:9.4f}  RMS {rms:.2f}")
        if extra[key]:
            current = combined_ff(ffjson, extra, os.path.join(work, "ff_with_autofit.json"))
    if miss["dihedral"]:
        rules, not_scanned, log = fit_torsions(structure, current, work, model, step, "--hold-neighbours" in sys.argv)
        for line in log.splitlines():
            if " V = " in line or "zero" in line:
                print("torsion  " + line)
        extra["dihedrals"] = rules
        for n in not_scanned:
            report["not_fitted"].append({"term": "dihedral " + "-".join(n["classes"]), "why": n["why"]})
    overlay = {"format": "caps-forcefield", "format_version": 1, "name": f"CAPS auto-fit (AIMNet2) for {os.path.basename(structure)}",
               "atom_types": [], **{k: v for k, v in extra.items() if v},
               "notes": [f"bonded terms {os.path.basename(ffjson)} lacks, fitted by python -m caps.autofit to AIMNet2 ({model}, "
                         f"{MODEL_LABEL}) relaxed scans against CAPS's MM energy without the term; not parameters of the force "
                         "field's authors. Entries marked needs_qm are in a chemistry AIMNet2 does not reproduce well."]}
    out = os.path.join(work, "autofit.json")
    json.dump(overlay, open(out, "w"), indent=1)
    # the check: everything found now, and charges neutral
    chk = os.path.join(work, "check.data")
    r = subprocess.run([T.caps_binary(), "ff", "apply", structure, "--ff", ffjson, "--overlay", out, "--charges", "types", "-o", chk],
                       capture_output=True, text=True)
    still = missing_terms(structure, combined_ff(ffjson, extra, os.path.join(work, "ff_with_autofit.json")))
    report["missing_after"] = {k: ["-".join(n) for n in v] if k in ("bond", "angle", "dihedral") else v for k, v in still.items()}
    report["complete"] = r.returncode == 0 and not any(still[k] for k in ("bond", "angle", "dihedral", "untyped"))
    report["needs_qm"] = [r_["name"] for r_ in extra["dihedrals"] if r_.get("needs_qm")]
    report["overlay"] = out
    json.dump(report, open(os.path.join(work, "report.json"), "w"), indent=1)
    if "--confirm-qm" in sys.argv and report["needs_qm"]:
        torsionfit("prepare", structure, ffjson, os.path.join(work, "qm_confirm"), "--level", "dlpno", "--nprocs", "8",
                   "--classes", ",".join(report["needs_qm"]), "--only")
    print(f"wrote {out}: {sum(len(v) for v in extra.values())} terms; complete {report['complete']}"
          + (f"; needs QM: {', '.join(report['needs_qm'])}" if report["needs_qm"] else "")
          + (f"; not fitted: {len(report['not_fitted'])}" if report["not_fitted"] else "")
          + (f"; charges: {still['charge']}" if still["charge"] else ""))
    return 0 if report["complete"] else 1


if __name__ == "__main__":
    sys.exit(main())
