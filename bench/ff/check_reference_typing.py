#!/usr/bin/env python3
"""CAPS's own typing and charges against DL_FIELD's, atom by atom, for every polymer of the library.

Each repeat unit (homopolymers of data/polymers/library.json) is built by CAPS as a short chain (DP 3, hydrogen end
caps) and written as .xyz. DL_FIELD (the reference, an external program: DLF_HOME, default ~/project/dl_f_4.13) types
it with its own chemical-group recognition and writes a DL_POLY FIELD file: per atom its force-field type and charge.
CAPS types the same atoms (same order) with its library force field. Reported per polymer: types that agree, the
largest charge difference, and the first disagreements.

usage: python3 bench/ff/check_reference_typing.py --ff opls2005 [--caps-ff opls2005] [--polymer substring] [--dp 3]
                                                   [--jobs 8] [--keep DIR] [--json OUT]
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


def reference(xyz, work):
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
        elif "* Seconday output files" in l:
            l = "none * Seconday output files"
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
    cq = [a.get("q", 0.0) for a in rep["atoms"]]
    # CAPS types may be named for their parameter class through aliases: compare by the reference's names
    diffs = [(i, rt[i], ct[i], rq[i], cq[i]) for i in range(n) if rt[i].lower() != ct[i].lower()]
    dq = max(abs(rq[i] - cq[i]) for i in range(n))
    worst = sorted(range(n), key=lambda i: -abs(rq[i] - cq[i]))[:3]
    return (pid, name, "ok" if not diffs and dq < 1e-3 else "differs",
            f"types {n - len(diffs)}/{n} · max |Δq| {dq:.3f}", {"type_diffs": diffs[:8], "charge_worst": [(i, rt[i], rq[i], cq[i]) for i in worst],
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
    print(f"\n{agree} of {len(rows)} polymers agree with the reference ({REF_FF} → CAPS {CAPS_FF}); work {WORK}")
    out = arg("--json", "")
    if out:
        json.dump([{"id": r[0], "name": r[1], "result": r[2], "detail": r[3], "more": r[4]} for r in rows], open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
