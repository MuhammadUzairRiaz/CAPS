"""CAPS's DL_POLY FIELD against the reference's, term by term, for every polymer of the library.

Each repeat unit is built by CAPS as a short chain (DP 3) and written as .xyz (as check_reference_typing.py does). The
reference (DL_FIELD, DLF_HOME, default ~/project/dl_f_4.13) writes its DL_POLY FIELD for it; CAPS types the same atoms
with its library force field and writes its own (caps ff apply --dlpoly). Compared, atom by atom and term by term:
names, masses and charges; bonds, angles, dihedrals (Fourier terms normalised: cos3 as three cos terms, zero amplitudes
dropped, a quadruple read either way), inversions (centre first, the others in any order), the 1-4 scale factors per
quadruple, and every van der Waals pair. Masses that differ alone are listed but not counted (the reference's
templates give chlorine 35.065 in places; CAPS writes its standard atomic weight).

usage: python3 bench/ff/check_dlpoly.py --ff opls2005 [--caps-ff opls2005] [--polymer substring] [--jobs 8] [--keep DIR]
"""
import os, re, subprocess, sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_reference_typing as R   # the reference runner and the polymer builder settings (same arguments)


def parse_field(path):
    """atoms [(name, mass, q)], and per term kind a dict key -> normalised terms (global 1-based indices)."""
    L = [l.rstrip("\n") for l in open(path)]
    atoms, T = [], {"bond": {}, "angle": {}, "dihed": {}, "inv": {}, "vdw": {}, "s14": {}}
    i = 0
    while i < len(L):
        w = L[i].split()
        if w[:1] and w[0].lower() == "nummols":
            nmol = int(w[1])
            i += 1
            m = re.match(r"atoms\s+(\d+)", L[i], re.I)
            na = int(m.group(1))
            mol_atoms = []
            for l in L[i + 1:i + 1 + na]:
                x = l.split()
                for _ in range(int(x[3]) if len(x) > 3 else 1):
                    mol_atoms.append((x[0].lower(), float(x[1]), float(x[2])))
            i += 1 + na
            sections = []
            while not L[i].lower().startswith("finish"):
                h = L[i].split()
                cnt = int(h[1])
                sections.append((h[0].lower(), L[i + 1:i + 1 + cnt]))
                i += 1 + cnt
            for c in range(nmol):
                off = len(atoms)
                atoms += mol_atoms
                for kind, rows in sections:
                    for r in rows:
                        x = r.split()
                        key = x[0].lower()
                        if kind == "bonds":
                            a, b = off + int(x[1]), off + int(x[2])
                            T["bond"].setdefault((min(a, b), max(a, b)), []).append((key,) + tuple(float(v) for v in x[3:]))
                        elif kind == "angles":
                            a, j, b = off + int(x[1]), off + int(x[2]), off + int(x[3])
                            T["angle"].setdefault((min(a, b), j, max(a, b)), []).append((key,) + tuple(float(v) for v in x[4:]))
                        elif kind in ("dihedrals", "dihedral"):
                            q = tuple(off + int(v) for v in x[1:5])
                            q = min(q, q[::-1])
                            p = [float(v) for v in x[5:]]
                            terms = T["dihed"].setdefault(q, [])
                            if key == "cos3":
                                for m, (A, d) in enumerate(((p[0] / 2, 0.0), (p[1] / 2, 180.0), (p[2] / 2, 0.0)), 1):
                                    if abs(A) > 1e-9:
                                        terms.append(("cos", round(A, 5), d, m))
                                s = (p[3], p[4])
                            elif key == "cos":
                                if abs(p[0]) > 1e-9:
                                    terms.append(("cos", round(p[0], 5), round(p[1] % 360.0, 4), int(round(p[2]))))
                                s = (p[3], p[4]) if len(p) > 4 else (0, 0)
                            else:
                                terms.append((key,) + tuple(round(v, 5) for v in p[:2]))
                                s = (p[2], p[3]) if len(p) > 3 else (0, 0)
                            e = T["s14"].get(q, (0.0, 0.0))
                            T["s14"][q] = (max(e[0], s[0]), max(e[1], s[1]))
                        elif kind == "inversions":
                            if all(abs(float(v)) < 1e-12 for v in x[5:6]):
                                continue   # a zero force constant: no energy (the reference writes some, CAPS leaves them out)
                            c, others = off + int(x[1]), tuple(sorted(off + int(v) for v in x[2:5]))
                            T["inv"].setdefault((c,) + others, []).append((key,) + tuple(round(float(v), 5) for v in x[5:]))
            i += 1
            continue
        if w[:1] and w[0].lower() == "vdw":
            for r in L[i + 1:i + 1 + int(w[1])]:
                x = r.split()
                T["vdw"][tuple(sorted((x[0].lower(), x[1].lower())))] = (x[2].lower(),) + tuple(float(v) for v in x[3:])
            i += 1 + int(w[1])
            continue
        i += 1
    for q in list(T["dihed"]):
        if not T["dihed"][q] and T["s14"].get(q, (0, 0)) == (0, 0):   # no energy and no 1-4 pair: nothing
            del T["dihed"][q]
            T["s14"].pop(q, None)
        else:
            T["dihed"][q] = sorted(T["dihed"][q])
    return atoms, T


def same(a, b):
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    return abs(a - b) <= 1e-4 * max(1.0, abs(a))


def compare(ref, caps):
    (ra, rt), (ca, ct) = ref, caps
    out = []
    if len(ra) != len(ca):
        return [f"atoms {len(ra)} / {len(ca)}"]
    bad = [k for k in range(len(ra)) if ra[k][0] != ca[k][0] or not same(ra[k][2], ca[k][2])]
    if bad:
        k = bad[0]
        out.append(f"{len(bad)} atoms differ (first {k + 1}: {ra[k]} / {ca[k]})")
    mass = sorted({(ra[k][0], ra[k][1], ca[k][1]) for k in range(len(ra)) if k not in bad and not same(ra[k][1], ca[k][1])})
    if mass:   # reported, not counted: the reference's templates give some elements a mass that is not theirs
        out.append("mass only: " + ", ".join(f"{t} {a:g} (reference) / {b:g}" for t, a, b in mass))
    for kind in ("bond", "angle", "dihed", "inv", "vdw", "s14"):
        A, B = rt[kind], ct[kind]
        keys = set(A) | set(B)
        diff = []
        for key in sorted(keys, key=str):
            x, y = A.get(key), B.get(key)
            if x is None or y is None:
                diff.append((key, x, y))
                continue
            xs, ys = (x if isinstance(x, list) else [x]), (y if isinstance(y, list) else [y])
            xs, ys = [t if isinstance(t, tuple) else (t,) for t in xs], [t if isinstance(t, tuple) else (t,) for t in ys]
            if kind == "s14":
                xs, ys = [x], [y]
            if len(xs) != len(ys) or any(len(p) != len(q) or not all(same(u, v) for u, v in zip(p, q)) for p, q in zip(sorted(xs, key=str), sorted(ys, key=str))):
                diff.append((key, x, y))
        if diff:
            key, x, y = diff[0]
            out.append(f"{kind} {len(diff)} of {len(keys)} differ (first {key}: reference {x} · CAPS {y})")
    return out


def one(job):
    pid, name, smiles = job
    import caps
    d = os.path.join(R.WORK, pid)
    os.makedirs(d, exist_ok=True)
    try:
        doc = caps.polymer(smiles, dp=R.DP, chains=1, density=0.02, seed=1)
        xyz = os.path.join(d, pid + ".xyz")
        doc.save(xyz)
        lines = open(xyz).read().split("\n")
        lines[1] = f"{pid} {name}"
        open(xyz, "w").write("\n".join(lines))
    except Exception as e:
        return (pid, name, "build", str(e)[:160])
    try:
        R.reference(xyz, os.path.join(d, "ref"))
    except Exception as e:
        return (pid, name, "reference", str(e)[:200])
    ffp = R.CAPS_FF if os.path.exists(R.CAPS_FF) else os.path.join(R.ROOT, "data", "forcefields", R.CAPS_FF + ".json")
    out = os.path.join(d, "caps_dlpoly")
    r = subprocess.run([R.CAPS_CLI, "ff", "apply", xyz, "--ff", ffp, "--dlpoly", out], capture_output=True, text=True, timeout=600)
    if not os.path.exists(os.path.join(out, "FIELD")):
        return (pid, name, "caps", (r.stdout + r.stderr).strip()[-200:])
    try:
        diffs = compare(parse_field(os.path.join(d, "ref", "dlf_output1", "dl_poly.FIELD")), parse_field(os.path.join(out, "FIELD")))
    except Exception as e:
        return (pid, name, "parse", str(e)[:200])
    real = [x for x in diffs if not x.startswith("mass only")]
    return (pid, name, "ok" if not real else "differs", " · ".join(diffs) if diffs else "every atom and term the same")


def main():
    import json
    lib = json.load(open(os.path.join(R.ROOT, "data", "polymers", "library.json")))
    jobs = [(e["id"], e["name"], e["smiles"]) for e in lib["polymers"] if R.ONLY.lower() in e["name"].lower() or R.ONLY == e["id"]]
    with ProcessPoolExecutor(R.JOBS) as ex:
        rows = list(ex.map(one, jobs))
    for pid, name, st, msg in rows:
        if st != "reference":
            print(f"{pid} {name[:44]:44s} {st:9s} {msg}")
    ok = sum(r[2] == "ok" for r in rows)
    ref = sum(r[2] == "reference" for r in rows)
    print(f"\n{ok} of {len(rows) - ref} FIELD files the same as the reference's ({R.REF_FF} → CAPS {R.CAPS_FF}); {ref} the reference cannot do; work {R.WORK}")


if __name__ == "__main__":
    main()
