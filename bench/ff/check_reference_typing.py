#!/usr/bin/env python3
"""CAPS's own typing and charges against DL_FIELD's, atom by atom, for every polymer of the library.

Each repeat unit (homopolymers of data/polymers/library.json) is built by CAPS as a short chain (DP 3, hydrogen end
caps) and written as .xyz. DL_FIELD (the reference, an external program: DLF_HOME, default ~/project/dl_f_4.13) types
it with its own chemical-group recognition and writes a DL_POLY FIELD file: per atom its force-field type and charge.
CAPS types the same atoms (same order) with its library force field. Reported per polymer: types that agree, the
largest charge difference, and the first disagreements.

With --energies, both sides also write LAMMPS files for the same structure and LAMMPS (LMP, default
~/lammps/build-class2/lmp) evaluates them in one common setting (a 400 Å box, Coulomb and van der Waals cut off at
25 Å, no k-space), each with its own styles, coefficients and 1-4 scaling: every bond, angle, dihedral, improper,
van der Waals and Coulomb term must give the same energy. Where only torsions differ, the reference is also run on the
same molecule with its atoms listed backwards: if its own torsion energy then changes to CAPS's value, the reference is
the one that depends on atom order (reported, counted as agreeing).

usage: python3 bench/ff/check_reference_typing.py --ff opls2005 [--caps-ff opls2005] [--polymer substring] [--dp 3]
                                                   [--energies] [--jobs 8] [--keep DIR] [--json OUT]
"""
import json, os, re, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "data", "python"))
os.environ.setdefault("CAPS_LIB", os.path.join(ROOT, "build", "capi", "libcaps.dylib"))
DLF = os.path.expanduser(os.environ.get("DLF_HOME", "~/project/dl_f_4.13"))
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
REF_FF = arg("--ff", "opls2005")
CAPS_FF = arg("--caps-ff", REF_FF)
ONLY = arg("--polymer", "")
DP = int(arg("--dp", "3"))
JOBS = int(arg("--jobs", "8"))
WORK = arg("--keep", "") or tempfile.mkdtemp(prefix="caps_reftype_")
ENERGIES = "--energies" in sys.argv
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
CAPS_CLI = os.path.join(ROOT, "build", "cli", "caps")
TERMS = ["E_bond", "E_angle", "E_dihed", "E_impro", "E_vdwl", "E_coul"]


def lammps_terms(infile, work):
    """Every energy term of an input, in the common setting (its own styles and coefficients)."""
    src = open(infile).read().split("\n")
    out = []
    for l in src:
        w = l.split()
        if not w:
            out.append(l); continue
        k = w[0]
        if k in ("kspace_style", "kspace_modify", "run", "minimize", "fix", "velocity", "dump", "thermo", "thermo_style", "thermo_modify", "timestep", "neighbor", "neigh_modify", "comm_modify"):
            continue
        if k == "pair_style":
            # the same potential without long-range Coulomb: coul/long → coul/cut, one cut-off for all
            l = re.sub(r"coul/long", "coul/cut", l)
            l = re.sub(r"(lj/[a-z0-9/]*coul/cut|lj/class2/coul/cut|lj/cut)\s+[0-9.]+(\s+[0-9.]+)?", lambda m: m.group(1) + " 25.0", l)
        if k == "pair_coeff":
            l = l.replace("coul/long", "coul/cut")
        out.append(l)
        if k == "read_data":
            out.append("change_box all x final -200 200 y final -200 200 z final -200 200 units box")
            out.append("boundary f f f") if False else None
    out += ["neighbor 2.0 bin", "neigh_modify one 10000 page 200000", "thermo_style custom step " + " ".join(t.lower().replace("e_", "e") for t in TERMS),
            "thermo_modify format float %.10f", "run 0"]
    txt = "\n".join(x for x in out if x is not None)
    txt = txt.replace("thermo_style custom step ebond eangle edihed eimpro evdwl ecoul", "thermo_style custom step ebond eangle edihed eimp evdwl ecoul")
    path = os.path.join(work, "common.in")
    open(path, "w").write(txt)
    r = subprocess.run([LMP, "-in", "common.in", "-log", "none"], cwd=work, capture_output=True, text=True, timeout=600)
    lines = r.stdout.split("\n")
    for i, l in enumerate(lines):
        if l.split()[:1] == ["Step"]:
            vals = lines[i + 1].split()
            return [float(v) for v in vals[1:7]]
    raise RuntimeError("LAMMPS: " + " | ".join((r.stdout + r.stderr).strip().split("\n")[-3:]))


def reference(xyz, work, lammps=False):
    """Types and charges from DL_FIELD for one .xyz (its atom order)."""
    os.makedirs(work, exist_ok=True)
    import shutil
    shutil.copy(xyz, work)
    with open(os.path.join(work, "dl_f_path"), "w") as f:
        f.write(f"library = {DLF}/lib/\nsolvent = {DLF}/solvent/\noutput  = ./\ncontrol = run.control\n")
    ctl = open(os.path.join(DLF, "dl_field.control")).read().split("\n")
    out = []
    for l in ctl:
        if "* Type of force field" in l:
            l = f"{REF_FF}  * Type of force field require"
        elif "* Configuration file" in l:
            l = f"{os.path.basename(xyz)}  * Configuration file."
        elif "* Periodic condition" in l and lammps:
            l = "1        * Periodic condition ? 0=no, other number = type of box (see below)"
        elif "* Seconday output files" in l:
            l = ("lammps" if lammps else "none") + " * Seconday output files"
        out.append(l)
    open(os.path.join(work, "run.control"), "w").write("\n".join(out))
    r = subprocess.run([os.path.join(DLF, "dl_field")], cwd=work, capture_output=True, text=True, timeout=600)
    field = os.path.join(work, "dlf_output1", "dl_poly.FIELD")
    if not os.path.exists(field):
        tail = (r.stdout + r.stderr).strip().split("\n")[-3:]
        raise RuntimeError("reference failed: " + " | ".join(tail))
    types, charges = [], []
    lines = open(field).read().split("\n")
    i = 0
    while i < len(lines):
        m = re.match(r"atoms\s+(\d+)", lines[i])
        if m:
            for l in lines[i + 1: i + 1 + int(m.group(1))]:
                w = l.split()
                if len(w) < 3:
                    tail = (r.stdout + r.stderr).strip().split("\n")[-2:]
                    raise RuntimeError("reference wrote no atom list (" + " | ".join(tail) + ")")
                reps = int(w[3]) if len(w) > 3 else 1
                for _ in range(reps):
                    types.append(w[0])
                    charges.append(float(w[2]))
            i += int(m.group(1))
        i += 1
    return types, charges


def one(job):
    pid, name, smiles = job
    import caps
    d = os.path.join(WORK, pid)
    os.makedirs(d, exist_ok=True)
    try:
        doc = caps.polymer(smiles, dp=DP, chains=1, density=0.02, seed=1)
        xyz = os.path.join(d, pid + ".xyz")
        doc.save(xyz)
        n = doc.atoms
        lines = open(xyz).read().split("\n")
        lines[1] = f"{pid} {name}"   # a plain title: the reference refuses a comment line that looks like numbers
        open(xyz, "w").write("\n".join(lines))
    except Exception as e:
        return (pid, name, "build", str(e)[:160], None)
    try:
        rt, rq = reference(xyz, os.path.join(d, "ref"))
    except Exception as e:
        return (pid, name, "reference", str(e)[:200], None)
    if len(rt) != n:
        return (pid, name, "reference", f"{len(rt)} atoms from the reference, {n} built", None)
    try:
        c = caps.open(xyz)
        rep = c.field.assign(CAPS_FF)
    except Exception as e:
        return (pid, name, "caps", str(e)[:200], None)
    ct = [a.get("type") or "?" for a in rep["atoms"]]
    energies = None
    order_dependent = ""
    if ENERGIES:
        try:
            rw = os.path.join(d, "ref_lmp")
            reference(xyz, rw, lammps=True)
            if "read_data" not in open(os.path.join(rw, "dlf_output1", "lammps.in")).read():
                return (pid, name, "reference", "the reference wrote no LAMMPS input for this structure (types and charges: "
                        f"{n - sum(rt[i].lower() != ct[i].lower() for i in range(n))}/{n} types agree)", None)
            re_ = lammps_terms(os.path.join(rw, "dlf_output1", "lammps.in"), os.path.join(rw, "dlf_output1"))
            cw = os.path.join(d, "caps_lmp")
            os.makedirs(cw, exist_ok=True)
            ffp = CAPS_FF if os.path.exists(CAPS_FF) else os.path.join(ROOT, "data", "forcefields", CAPS_FF + ".json")
            r2 = subprocess.run([CAPS_CLI, "ff", "apply", xyz, "--ff", ffp, "-o", os.path.join(cw, "caps.data"), "--lammps-input", os.path.join(cw, "caps.in")],
                                capture_output=True, text=True, timeout=600)
            if not os.path.exists(os.path.join(cw, "caps.in")):
                raise RuntimeError("caps ff apply: " + (r2.stdout + r2.stderr).strip()[-200:])
            ce = lammps_terms(os.path.join(cw, "caps.in"), cw)
            energies = {t: (a, b) for t, a, b in zip(TERMS, re_, ce)}
            bad = [t for t, (a, b) in energies.items() if abs(a - b) / max(1.0, abs(a)) >= 1e-4]
            if bad and set(bad) <= {"E_dihed", "E_impro"}:
                # does the reference itself give the same energy for the same molecule with its atoms listed backwards?
                rv = os.path.join(d, "ref_rev")
                os.makedirs(rv, exist_ok=True)
                L = open(xyz).read().rstrip("\n").split("\n")
                open(os.path.join(rv, pid + ".xyz"), "w").write("\n".join(L[:2] + L[2:2 + n][::-1]) + "\n")
                reference(os.path.join(rv, pid + ".xyz"), os.path.join(rv, "ref"), lammps=True)
                rr = lammps_terms(os.path.join(rv, "ref", "dlf_output1", "lammps.in"), os.path.join(rv, "ref", "dlf_output1"))
                rev = dict(zip(TERMS, rr))
                if any(abs(rev[t] - energies[t][0]) / max(1.0, abs(rev[t])) >= 1e-4 for t in bad) and \
                        all(abs(rev[t] - energies[t][1]) / max(1.0, abs(rev[t])) < 1e-4 for t in bad):
                    order_dependent = ", ".join(f"{t[2:]} {energies[t][0]:.3f} or {rev[t]:.3f}" for t in bad)
                    energies = {t: (energies[t][1] if t in bad else a, b) for t, (a, b) in energies.items()}
                elif any(abs(rev[t] - energies[t][0]) / max(1.0, abs(rev[t])) >= 1e-4 for t in bad):
                    order_dependent = "!" + ", ".join(f"{t[2:]} {energies[t][0]:.3f} or {rev[t]:.3f} (CAPS {energies[t][1]:.3f})" for t in bad)
        except Exception as e:
            energies = {"error": str(e)[:200]}
    cq = [a.get("q", 0.0) for a in rep["atoms"]]
    # CAPS types may be named for their parameter class through aliases: compare by the reference's names
    diffs = [(i, rt[i], ct[i], rq[i], cq[i]) for i in range(n) if rt[i].lower() != ct[i].lower()]
    dq = max(abs(rq[i] - cq[i]) for i in range(n))
    worst = sorted(range(n), key=lambda i: -abs(rq[i] - cq[i]))[:3]
    ediff = ""
    eok = True
    if energies is not None:
        if "error" in energies:
            ediff, eok = " · energies: " + energies["error"], False
        else:
            eworst = max(abs(a - b) / max(1.0, abs(a)) for a, b in energies.values())
            eok = eworst < 1e-4
            if order_dependent.startswith("!"):
                ediff = " · the reference's own energy depends on the atom order: " + order_dependent[1:]
            else:
                ediff = f" · energy terms {eworst:.1e}" + (f" (the reference's own {order_dependent} with the atoms listed backwards; CAPS gives the latter)" if order_dependent else "") + ("" if eok else " (" + ", ".join(f"{t[2:]} {a:.3f}/{b:.3f}" for t, (a, b) in energies.items() if abs(a - b) / max(1.0, abs(a)) >= 1e-4) + ")")
    status = "ok" if not diffs and dq < 1e-3 and eok else "ref-order" if not diffs and dq < 1e-3 and order_dependent.startswith("!") else "differs"
    return (pid, name, status,
            f"types {n - len(diffs)}/{n} · max |Δq| {dq:.3f}" + ediff, {"type_diffs": diffs[:8], "charge_worst": [(i, rt[i], rq[i], cq[i]) for i in worst],
                                                                 "caps_complete": rep.get("complete"), "charges": rep.get("charges")})


def main():
    lib = json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))
    jobs = [(e["id"], e["name"], e["smiles"]) for e in lib["polymers"] if ONLY.lower() in e["name"].lower() or ONLY == e["id"]]
    with ProcessPoolExecutor(JOBS) as ex:
        rows = list(ex.map(one, jobs))
    agree = sum(r[2] == "ok" for r in rows)
    for pid, name, st, msg, det in rows:
        print(f"{pid} {name[:44]:44s} {st:9s} {msg}")
        if det and st != "ok":
            for i, a, b, qa, qb in det["type_diffs"][:4]:
                print(f"      atom {i + 1}: reference {a} ({qa:+.3f}) · CAPS {b} ({qb:+.3f})")
            if not det["type_diffs"]:
                for i, a, qa, qb in det["charge_worst"][:2]:
                    print(f"      atom {i + 1} {a}: reference {qa:+.4f} · CAPS {qb:+.4f}")
    order = sum(r[2] == "ref-order" for r in rows)
    refail = sum(r[2] == "reference" for r in rows)
    print(f"\n{agree} of {len(rows)} polymers agree with the reference ({REF_FF} → CAPS {CAPS_FF})"
          + (f"; {order} where the reference's own energy depends on the atom order" if order else "")
          + (f"; {refail} the reference cannot do" if refail else "") + f"; work {WORK}")
    out = arg("--json", "")
    if out:
        json.dump([{"id": r[0], "name": r[1], "result": r[2], "detail": r[3], "more": r[4]} for r in rows], open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
