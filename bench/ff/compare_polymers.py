#!/usr/bin/env python3
"""CAPS against DL_FIELD on the polymer library: the same oligomer typed and parameterised by both, compared atom by
atom and in LAMMPS.

For each repeat unit of data/polymers/library.json CAPS grows one isolated chain (DP 4, H ends) and writes it as plain
xyz. DL_FIELD reads that xyz (its automatic typing) and writes LAMMPS files; CAPS applies the matching force field of
its library to the same xyz and writes its own. Compared:

  · types: the force-field type of every atom (DL_FIELD's Atoms comments against CAPS's Masses labels);
  · charges: the largest per-atom difference (DL_FIELD takes OPLS / PCFF / CVFF / COMPASS charges from bond increments,
    CAPS the force field's own or Gasteiger–Marsili, as its automatic charge mode does);
  · energies: both LAMMPS files run in LAMMPS (single point), bonded terms compared (bond, angle, dihedral, improper),
    with the counts of each kind of term; non-bonded terms depend on the cut-off and electrostatics each program writes
    and are shown, not compared.

usage: compare_polymers.py [--ff opls2005,pcff,...] [--polymer substring] [--keep DIR] [--json OUT]
Needs DL_FIELD (DLFIELD, default ~/project/dl_f_4.13) and LAMMPS with CLASS2 (LMP).
"""
import json, math, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
WORK = arg("--keep", "") or tempfile.mkdtemp(prefix="caps_vs_dlf_")
ONLY_P = arg("--polymer", "")
# DL_FIELD keyword → CAPS library id
FFS = {"opls2005": "opls2005", "pcff": "pcff", "cvff": "cvff", "compass": "compass", "amber16_gaff": "gaff-amber16",
       "amber25_gaff": "gaff-amber25", "charmm36_cgenff": "cgenff", "dreiding": "dreiding"}
if "--ff" in sys.argv:
    FFS = {k: v for k, v in FFS.items() if k in arg("--ff", "").split(",") or v in arg("--ff", "").split(",")}


def oligomer(pid, smiles):
    d = os.path.join(WORK, "mol")
    os.makedirs(d, exist_ok=True)
    ext, xyz = os.path.join(d, pid + ".ext.xyz"), os.path.join(d, pid + ".xyz")
    if not os.path.exists(xyz):
        r = subprocess.run([CAPS, "grow", "--units", smiles, "--chains", "1", "--dp", "4", "--density", "0.01", "--seed", "3", "-o", ext],
                           capture_output=True, text=True)
        if r.returncode:
            return None
        lines = open(ext).read().splitlines()
        # plain xyz: DL_FIELD takes an extended-xyz title for coordinates
        open(xyz, "w").write("\n".join([lines[0], f"{pid} oligomer"] + lines[2:]) + "\n")
    return xyz


def atom_types_from_data(path):
    """per-atom type names (by atom id order) and charges from a LAMMPS data file with Masses labels."""
    label, types, charges, sec = {}, {}, {}, None
    for raw in open(path):
        l = raw.split("#")[0].strip()
        if not l:
            continue
        if re.match(r"^[A-Z][A-Za-z ]+$", l):
            sec = l.strip()
            continue
        if sec == "Masses" and "#" in raw:
            w = l.split()
            if len(w) >= 2 and w[0].isdigit():
                label[int(w[0])] = raw.split("#")[1].strip().split()[0]
            continue
        w = l.split()
        if sec == "Atoms" and len(w) >= 7:
            types[int(w[0])] = int(w[2])
            charges[int(w[0])] = float(w[3])
    ids = sorted(types)
    return [label.get(types[i], str(types[i])) for i in ids], [charges[i] for i in ids]


def lammps_terms(d, script):
    txt = open(os.path.join(d, script)).read()
    txt = re.sub(r"(?m)^thermo_style.*$", "", txt)
    txt = re.sub(r"(?m)^run .*$", "", txt)
    txt += "\nthermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong\nthermo_modify format float %.8f\nrun 0\n"
    open(os.path.join(d, "check.in"), "w").write(txt)
    r = subprocess.run([LMP, "-in", "check.in", "-log", "none"], cwd=d, capture_output=True, text=True, timeout=300)
    if r.returncode:
        err = [l for l in (r.stdout + r.stderr).splitlines() if "ERROR" in l]
        return None, (err[0] if err else "LAMMPS failed")[:150]
    lines = r.stdout.splitlines()
    k = next((i for i, l in enumerate(lines) if l.split()[:2] == ["Step", "PotEng"]), None)
    if k is None:
        return None, "no thermo"
    head, vals = lines[k].split(), [float(x) for x in lines[k + 1].split()]
    return dict(zip(head, vals)), ""


def counts(path):
    c = {}
    for l in open(path):
        m = re.match(r"\s*(\d+)\s+(bonds|angles|dihedrals|impropers)\s*$", l)
        if m:
            c[m.group(2)] = int(m.group(1))
    return c


def case(pid, pname, xyz, kw, ffid):
    d = os.path.join(WORK, f"{pid}__{ffid}")
    os.makedirs(d, exist_ok=True)
    out = {"id": pid, "polymer": pname, "forcefield": ffid}
    # DL_FIELD
    ddir = os.path.join(d, "dlfield")
    r = subprocess.run([sys.executable, os.path.join(HERE, "run_dlfield.py"), kw, xyz, ddir], capture_output=True, text=True, timeout=600)
    if r.returncode or not os.path.exists(os.path.join(ddir, "lammps1.data")):
        out["dlfield"] = "failed: " + (r.stdout.strip().splitlines() or ["?"])[-1][:150]
    # CAPS
    cdir = os.path.join(d, "caps")
    os.makedirs(cdir, exist_ok=True)
    rc = subprocess.run([CAPS, "ff", "apply", xyz, "--ff", os.path.join(ROOT, "data", "forcefields", ffid + ".json"), "--charges", "auto",
                         "-o", os.path.join(cdir, "system.data"), "--lammps-input", os.path.join(cdir, "system.in")], capture_output=True, text=True)
    if rc.returncode:
        why = (rc.stderr.strip().splitlines() or rc.stdout.strip().splitlines() or ["?"])[-1]
        out["caps"] = "failed: " + why[:200]
    if "dlfield" in out or "caps" in out:
        return out
    dt, dq = atom_types_from_data(os.path.join(ddir, "lammps1.data"))
    ct, cq = atom_types_from_data(os.path.join(cdir, "system.data"))
    if len(dt) != len(ct):
        out["result"] = f"atom counts differ: DL_FIELD {len(dt)}, CAPS {len(ct)}"
        return out
    diff = [(i + 1, a, b) for i, (a, b) in enumerate(zip(dt, ct)) if a != b]
    out["atoms"] = len(dt)
    out["types_same"] = len(dt) - len(diff)
    out["type_diffs"] = sorted({f"{a}→{b}" for _, a, b in diff})
    out["charge_max_diff"] = max(abs(a - b) for a, b in zip(dq, cq))
    out["counts"] = {"dlfield": counts(os.path.join(ddir, "lammps1.data")), "caps": counts(os.path.join(cdir, "system.data"))}
    ed, errd = lammps_terms(ddir, "lammps.in")
    ec, errc = lammps_terms(cdir, "system.in")
    if ed is None or ec is None:
        out["result"] = "LAMMPS: " + (errd or errc)
        return out
    terms = {k: (ed[k], ec[k]) for k in ("E_bond", "E_angle", "E_dihed", "E_impro", "E_vdwl", "E_coul")}
    out["terms"] = terms
    out["bonded_max_rel"] = max(abs(a - b) / max(1.0, abs(a)) for k, (a, b) in terms.items() if k in ("E_bond", "E_angle", "E_dihed", "E_impro"))
    return out


if __name__ == "__main__":
    lib = json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))
    polys = [(p["id"], p["name"], p["smiles"]) for p in lib["polymers"] if ONLY_P.lower() in (p["id"] + " " + p["name"]).lower()]
    res = []
    for pid, pname, smi in polys:
        xyz = oligomer(pid, smi)
        if not xyz:
            res.append({"id": pid, "polymer": pname, "forcefield": "-", "result": "CAPS could not grow the oligomer"})
            continue
        for kw, ffid in FFS.items():
            res.append(case(pid, pname, xyz, kw, ffid))
            r = res[-1]
            tag = ("DL_FIELD " + r["dlfield"]) if "dlfield" in r else ("CAPS " + r["caps"]) if "caps" in r else r.get("result") or \
                  f"types {r['types_same']}/{r['atoms']} · bonded {r['bonded_max_rel']:.1e} · q {r['charge_max_diff']:.3f}"
            print(f"{pid} {pname[:28]:28s} {ffid:13s} {tag[:170]}", flush=True)
    print()
    for kw, ffid in FFS.items():
        mine = [r for r in res if r["forcefield"] == ffid]
        both = [r for r in mine if "types_same" in r]
        atoms = sum(r["atoms"] for r in both)
        same = sum(r["types_same"] for r in both)
        full = sum(1 for r in both if r["types_same"] == r["atoms"])
        bon = sum(1 for r in both if r.get("bonded_max_rel", 1) < 1e-4)
        dlf_fail = sum(1 for r in mine if "dlfield" in r)
        caps_fail = sum(1 for r in mine if "caps" in r)
        diffs = {}
        for r in both:
            for x in r["type_diffs"]:
                diffs[x] = diffs.get(x, 0) + 1
        top = ", ".join(f"{k} ×{v}" for k, v in sorted(diffs.items(), key=lambda kv: -kv[1])[:6])
        print(f"{ffid:13s} both typed {len(both):3d} · atoms same {same}/{atoms} ({100 * same / max(1, atoms):.1f} %) · polymers identical {full} · "
              f"bonded energy equal {bon} · DL_FIELD failed {dlf_fail} · CAPS failed {caps_fail}" + (f"\n              most frequent differences: {top}" if top else ""))
    if "--json" in sys.argv:
        json.dump(res, open(arg("--json", "compare.json"), "w"), indent=1)
