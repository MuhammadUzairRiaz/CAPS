#!/usr/bin/env python3
"""Every polymer of the library against every force field that can describe a polymer: does CAPS build a correct cell?

For each repeat unit (and copolymer) of data/polymers/library.json and each force field of data/forcefields/catalogue.json
(plus the built-in UFF), a recipe builds 3 chains of DP 10, types them (charges automatic: the force field's own, else
Gasteiger–Marsili, else QEq), grows them at 0.3 g/cm³, relaxes them (L-BFGS to |F|max 1 kcal/mol/Å) and writes LAMMPS files.
Then the cell is checked, independently of CAPS:

  · bonds: every bond within 0.25 Å of its force field's equilibrium length r0 (from the Bond Coeffs of the data file);
  · contacts: no two atoms that are not bonded or 1-3 closer than 1.0 Å (minimum image);
  · LAMMPS reads the data file and the input script and returns a finite energy (run 0), with at most +30 kcal/mol
    per atom (a relaxed cell has a small or negative potential energy per atom).

Force fields that cannot describe an organic polymer by design (water models, inorganic crystals, coarse-grained
bead models), library entries without a parameter file (aliases, templates) and force fields CAPS has no typing rules
for yet are listed and skipped with the reason. Failures print CAPS's own reason (untyped atoms, missing
parameters) so each can be traced to a typing rule or a parameter the library lacks.

usage: python3 bench/polymers/check_matrix.py [--ff substring] [--polymer substring] [--jobs 8] [--keep DIR] [--json OUT]
"""
import json, math, os, re, subprocess, sys, tempfile
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
ONLY_FF, ONLY_P = arg("--ff", ""), arg("--polymer", "")
JOBS = int(arg("--jobs", "8"))
WORK = arg("--keep", "") or tempfile.mkdtemp(prefix="caps_matrix_")

# not meant for organic polymers: the reason is printed
NOT_POLYMER = [(r"^(spce|tip3p|tip5p|mw)", "water model"), (r"^inorganic", "inorganic crystals"),
               (r"^(martini|drymartini|sdk|cooke-deserno)", "coarse-grained beads"), (r"^graphene", "carbon sheets only (no H)"),
               (r"^(spc_|tip3p_|tip5p_)", "water model on OPLS-AA")]


def forcefields():
    cat = json.load(open(os.path.join(ROOT, "data", "forcefields", "catalogue.json")))
    out = [("uff", "UFF (built in)", None)]
    for e in cat["forcefields"]:
        skip = next((why for pat, why in NOT_POLYMER if re.search(pat, e["id"])), None)
        if not skip and not isinstance(e.get("file"), str):
            skip = f"no parameter file ({e.get('status', '')}: {e.get('notes', '')})"
        if not skip and not isinstance(e.get("typing"), dict):
            skip = "no CAPS typing rules yet (types by hand only)"
        out.append((e["id"], e.get("name", e["id"]), skip))
    return out


def polymers():
    lib = json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))
    by_id = {p["id"]: p for p in lib["polymers"]}
    out = [(p["id"], p["name"], {"smiles": p["smiles"]}) for p in lib["polymers"]]
    for c in lib.get("copolymers", []):
        units = [by_id[u]["smiles"] if u in by_id else u for u in c["units"]]
        spec = {"units": units, "sequence": c.get("sequence", "random")}
        if c.get("weights"):
            spec["weights"] = c["weights"]
        out.append((c["id"], c["name"], spec))
    return out


def recipe(pid, spec, ff):
    lines = ["recipe: 1", f"name: {pid}", "build:", "  polymer:"]
    if "smiles" in spec:
        lines.append(f"    smiles: {json.dumps(spec['smiles'])}")
    else:
        lines.append("    units: [" + ", ".join(json.dumps(u) for u in spec["units"]) + "]")
        lines.append(f"    sequence: {spec['sequence']}")
        if "weights" in spec:
            lines.append("    weights: [" + ", ".join(str(w) for w in spec["weights"]) + "]")
    lines += ["    dp: 10", "    chains: 3", f"type: {{ forcefield: {ff}, charges: auto }}", "grow: { density: 0.3, seed: 7 }",
              "relax: { method: lbfgs, fmax: 1.0 }", "export: [lammps]"]
    return "\n".join(lines) + "\n"


def read_data(path):
    """box lengths, atoms {id: (type, x, y, z)}, bonds [(type, i, j)], bond coeffs {type: tokens}."""
    L = [0, 0, 0]
    atoms, bonds, coeffs = {}, [], {}
    sec = None
    for raw in open(path):
        l = raw.split("#")[0].strip()
        m = re.match(r"(\S+)\s+(\S+)\s+([xyz])lo\s+[xyz]hi", l)
        if m:
            L["xyz".index(m.group(3))] = float(m.group(2)) - float(m.group(1))
            continue
        if not l:
            continue
        if re.match(r"^[A-Z][A-Za-z ]+$", l):
            sec = l
            continue
        w = l.split()
        if sec == "Atoms" and len(w) >= 7:
            atoms[int(w[0])] = (int(w[2]), float(w[4]), float(w[5]), float(w[6]))
        elif sec == "Bonds" and len(w) >= 4:
            bonds.append((int(w[1]), int(w[2]), int(w[3])))
        elif sec == "Bond Coeffs" and len(w) >= 2:
            coeffs[int(w[0])] = w[1:]
    return L, atoms, bonds, coeffs


def r0_of(tokens, style):
    st = style.split()
    if st and st[0] == "hybrid":   # hybrid: the sub-style is the first token of the coefficients
        sub, tokens = tokens[0], tokens[1:]
    else:
        sub = st[0] if st else "harmonic"
    v = [float(x) for x in tokens if re.match(r"^-?[\d.]+(e[-+]?\d+)?$", x, re.I)]
    if sub in ("harmonic", "gromos"):
        return v[1]
    if sub == "class2":
        return v[0]
    if sub == "morse":
        return v[2]
    return None   # zero or other: no equilibrium length to check


def geometry(data, script):
    L, atoms, bonds, coeffs = read_data(data)
    style = ""
    for l in open(script):
        m = re.match(r"bond_style\s+(.*)", l)
        if m:
            style = m.group(1).strip()
    mi = lambda d, k: d - L[k] * round(d / L[k]) if L[k] > 0 else d
    dist = lambda a, b: math.sqrt(sum(mi(atoms[a][k + 1] - atoms[b][k + 1], k) ** 2 for k in range(3)))
    worst_bond, worst_what = 0.0, ""
    nb = {i: set() for i in atoms}
    for t, i, j in bonds:
        nb[i].add(j), nb[j].add(i)
        r0 = r0_of(coeffs.get(t, []), style) if t in coeffs else None
        if r0 is None:
            continue
        d = abs(dist(i, j) - r0)
        if d > worst_bond:
            worst_bond, worst_what = d, f"bond {i}-{j} {dist(i, j):.2f} Å vs r0 {r0:.2f}"
    # contacts: cell list on 2.5 Å cells, pairs that are neither bonded nor 1-3
    near = {i: nb[i] | {k for j in nb[i] for k in nb[j]} for i in atoms}
    cs = 2.5
    n = [max(1, int(L[k] // cs)) for k in range(3)]
    cells = {}
    for i, a in atoms.items():
        key = tuple(int((a[k + 1] % L[k]) / L[k] * n[k]) % n[k] if L[k] > 0 else 0 for k in range(3))
        cells.setdefault(key, []).append(i)
    closest, pair = 9e9, None
    for key, members in cells.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    other = cells.get(((key[0] + dx) % n[0], (key[1] + dy) % n[1], (key[2] + dz) % n[2]), [])
                    for i in members:
                        for j in other:
                            if j <= i or j in near[i]:
                                continue
                            d = dist(i, j)
                            if d < closest:
                                closest, pair = d, (i, j)
    return worst_bond, worst_what, closest, pair, len(atoms)


def lammps_energy(d, stem):
    r = subprocess.run([LMP, "-in", stem + ".in", "-log", "none"], cwd=d, capture_output=True, text=True, timeout=300)
    if r.returncode:
        err = [l for l in (r.stdout + r.stderr).splitlines() if "ERROR" in l]
        return None, (err[0] if err else (r.stdout + r.stderr).strip().splitlines()[-1])[:160]
    lines = r.stdout.splitlines()
    k = next((i for i, l in enumerate(lines) if l.split()[:2] == ["Step", "PotEng"]), None)
    if k is None:
        return None, "no thermo output"
    return float(lines[k + 1].split()[1]), ""


def case(job):
    pid, pname, spec, ffid = job
    d = os.path.join(WORK, f"{pid}__{ffid}")
    os.makedirs(d, exist_ok=True)
    rp = os.path.join(d, "r.yaml")
    open(rp, "w").write(recipe(pid, spec, ffid))
    r = subprocess.run([CAPS, "run", rp, "--out", d, "--threads", "1"], capture_output=True, text=True, timeout=600)
    out = r.stdout + r.stderr
    if r.returncode:
        why = next((l for l in out.splitlines() if "failed" in l), out.strip().splitlines()[-1] if out.strip() else "?")
        why = re.sub(r"\s+", " ", why).strip()
        stage = "type" if r.returncode == 3 or "type " in why[:12] else "build" if "build" in why[:12] else "run"
        return pid, pname, ffid, "fail-" + stage, why[:220]
    data, script = os.path.join(d, pid + ".data"), os.path.join(d, pid + ".in")
    try:
        wb, what, closest, pair, n = geometry(data, script)
    except Exception as e:
        return pid, pname, ffid, "fail-check", f"could not read the data file: {e}"
    e, err = lammps_energy(d, pid)
    probs = []
    if wb > 0.25:
        probs.append(f"stretched {what}")
    if closest < 1.0:
        probs.append(f"overlap {closest:.2f} Å between atoms {pair[0]} and {pair[1]}")
    if e is None:
        probs.append("LAMMPS: " + err)
    elif not math.isfinite(e) or e / n > 30:
        probs.append(f"LAMMPS energy {e:.1f} kcal/mol ({e / n:.1f} per atom)")
    detail = f"{n} atoms · worst bond {wb:.3f} Å · closest contact {closest:.2f} Å · LAMMPS {e if e is None else round(e, 1)} kcal/mol"
    return pid, pname, ffid, ("ok" if not probs else "bad-structure"), (detail if not probs else "; ".join(probs))


if __name__ == "__main__":
    ffs = [f for f in forcefields() if ONLY_FF.lower() in f[0].lower()]
    polys = [p for p in polymers() if ONLY_P.lower() in (p[0] + " " + p[1]).lower()]
    skipped = [(fid, why) for fid, _, why in ffs if why]
    jobs = [(pid, pname, spec, fid) for fid, _, why in ffs if not why for pid, pname, spec in polys]
    print(f"{len(polys)} polymers × {len(ffs) - len(skipped)} force fields = {len(jobs)} cells ({len(skipped)} force fields skipped by design); work {WORK}", flush=True)
    res = []
    with ProcessPoolExecutor(JOBS) as ex:
        for k, r in enumerate(ex.map(case, jobs, chunksize=4)):
            res.append(r)
            if (k + 1) % 100 == 0:
                print(f"  {k + 1} / {len(jobs)}", flush=True)
    # per force field: counts and the first reasons
    print()
    for fid, name, why in ffs:
        if why:
            print(f"{fid:34s} skipped: {why}")
            continue
        mine = [r for r in res if r[2] == fid]
        ok = sum(1 for r in mine if r[3] == "ok")
        kinds = {}
        for r in mine:
            if r[3] != "ok":
                kinds.setdefault(r[3], []).append(r)
        print(f"{fid:34s} {ok:3d} / {len(mine)} correct" + "".join(f" · {len(v)} {k}" for k, v in sorted(kinds.items())))
        for k, v in sorted(kinds.items()):
            for r in v[:3]:
                print(f"      {k}: {r[1]} — {r[4]}")
    tot = sum(1 for r in res if r[3] == "ok")
    print(f"\n{tot} of {len(res)} cells correct")
    if "--json" in sys.argv:
        json.dump([dict(zip(["id", "polymer", "forcefield", "result", "detail"], r)) for r in res], open(arg("--json", "matrix.json"), "w"), indent=1)
    sys.exit(0 if tot == len(res) else 1)
