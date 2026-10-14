"""Which torsions does the chemistry CAPS is used for need, and which of them does a force field cover?

    python -m caps.census [--ff oplsaa2024] [--out DIR] [--only TEXT] [--no-react] [--no-copolymers] [--radius 2]

builds the chemistry CAPS is meant for: every polymer of the library as a DP-3 oligomer, every copolymer of the library
with each unit next to each other, the reactants of CAPS React's worked cases (curatives, crosslinkers, hardeners) and
their products cured to full conversion and relaxed (bench/ff/react_coverage.py's cure, relaxed=True). Each structure is typed with the
force field, and the census lists

  torsion classes    every proper torsion by the class names the force field looks torsions up with, and where its
                     parameters come from: "published" (the force field's own table), "CAPS fit" (a QM or AIMNet2 fit
                     CAPS added to the library, marked in its comment) or "missing"
  rotor environments every rotatable bond (not in a ring, not through a linear atom, a neighbour on both ends) keyed by
                     its chemical environment: the force-field types of the atoms within RADIUS bonds of the bond, ring
                     membership and cut bonds, hashed Weisfeiler-Lehman style so that the same environment in two
                     molecules is one key. One QM scan of the bond fits every torsion about it, so the number of
                     environments is the number of scans a library of this chemistry needs. The status of an
                     environment is the worst of its torsions'

and writes one capped fragment per environment (the atoms within RADIUS bonds, each cut bond capped by H, as
caps.torsionfit prepare cuts them): the set a QM campaign or an ML training set starts from.

Output in DIR (bench/ff/census/<ff> by default): census.json (classes, environments with an example each, structures),
summary.md (the counts and the campaign size) and fragments/<env>.mol2.
"""
import collections, hashlib, importlib.util, json, math, os, sys, tempfile

from . import core as caps
from .torsionfit import (CAP_BOND, arg, classes, element_of_mass, in_ring, min_image, neighbours, propers, read_box,
                         read_data, write_mol2)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
CAPS_FIT = ("CAPS fit", "CAPS auto-fit")
RANK = {"published": 0, "CAPS fit": 1, "missing": 2}
SCAN_POINTS = 24   # a 15 degree scan, as caps.torsionfit prepare writes it
VALENCE = {"H": 1, "C": 4, "N": 3, "O": 2, "S": 2, "Si": 4, "P": 3, "F": 1, "Cl": 1, "Br": 1, "I": 1}
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Br": 35, "I": 53}
SMALL_RING = 8     # rings up to this size are kept whole; larger cycles (a crosslinked network's) are cut
MAX_FRAGMENT = 80  # atoms; a fragment that would grow past this is cut where it stands and marked


# ------------------------------------------------------------------ the force field, with and without CAPS's fits
def ff_paths(ff, workdir):
    """the force field's file, and a copy without the terms CAPS fitted (their comment says so): the published table.
    The copy's relative typing path is made absolute, or CAPS cannot find the typing rules from the work directory"""
    full = caps._forcefield_path(ff)
    d = json.load(open(full))
    for sec in ("bonds", "angles", "dihedrals", "impropers"):
        d[sec] = [t for t in d.get(sec, []) if not any(m in t.get("comment", "") for m in CAPS_FIT)]
    if isinstance(d.get("typing"), str) and not os.path.isabs(d["typing"]):
        d["typing"] = os.path.normpath(os.path.join(os.path.dirname(full), d["typing"]))
    pub = os.path.join(workdir, "published-only.json")
    json.dump(d, open(pub, "w"))
    return full, pub


def missing_terms(report, cl):
    """the missing terms of an assignment: {kind: set of class tuples, each in its canonical direction}"""
    out = collections.defaultdict(set)
    for m in report.get("missing", []):
        w = m.split()
        key = tuple(x.split("(")[1].rstrip(")") if "(" in x else cl.get(x, x) for x in w[1:])
        out[w[0]].add(min(key, key[::-1]))
    return out


# ------------------------------------------------------------------ the chemistry
def library():
    lib = json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))
    return lib["polymers"], lib["copolymers"]


def react_cases():
    """the worked React cases of bench/ff/react_coverage.py (reactants and the cure), or none outside the repository"""
    path = os.path.join(ROOT, "bench", "ff", "react_coverage.py")
    if not os.path.exists(path):
        return [], None
    argv, sys.argv = sys.argv, [path]   # the script reads its own options from sys.argv on import
    try:
        spec = importlib.util.spec_from_file_location("react_coverage", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    finally:
        sys.argv = argv
    return mod.CASES, mod.cure


def structures(react=True, copolymers=True):
    """(kind, name, builder) for each structure of the census; a builder makes the Document when called"""
    polymers, copolys = library()
    by_id = {p["id"]: p["smiles"] for p in polymers}
    for p in polymers:
        yield "polymer", f"{p['id']} {p['name']}", lambda s=p["smiles"]: caps.polymer(s, dp=3, chains=1, density=0.02, seed=1)
    if copolymers:
        for c in copolys:
            units = [by_id.get(u, u) for u in c["units"]]
            pattern = "".join(chr(65 + k) for k in range(len(units)))
            yield "copolymer", f"{c['id']} {c['name']}", (lambda u=units, p=pattern:
                caps.polymer(u, dp=2 * len(u), chains=1, density=0.02, seed=1, pattern=p))
    if react:
        cases, cure = react_cases()
        seen = set()
        for _, _, mols, _ in cases:
            for smi, _ in mols:
                if smi not in seen:
                    seen.add(smi)
                    yield "reactant", smi, lambda s=smi: caps.build.smiles(s)
        for template, name, mols, box in cases:
            # cured as a real run cures (relaxed every cycle): the products' geometries are cut into QM fragments
            yield "product", f"{template}: {name}", lambda t=template, m=mols, b=box: cure(t, m, b, 1.0, relaxed=True)


# ------------------------------------------------------------------ environments
def ring_atoms(nb):
    return {x for x in range(len(nb)) for y in nb[x] if in_ring(nb, x, y)}


def linear(pos, nb, x, box=None):
    """an atom whose two bonds are collinear (an sp carbon): no torsion turns about it"""
    if len(nb[x]) != 2:
        return False
    a, c = nb[x]
    u = min_image([pos[a][k] - pos[x][k] for k in range(3)], box)
    v = min_image([pos[c][k] - pos[x][k] for k in range(3)], box)
    cosang = sum(u[k] * v[k] for k in range(3)) / (math.sqrt(sum(t * t for t in u)) * math.sqrt(sum(t * t for t in v)))
    return cosang < math.cos(math.radians(170))


def small_rings(nb, limit=SMALL_RING):
    """the smallest ring through each ring bond, up to `limit` atoms, as sets of atoms"""
    out = set()
    for u in range(len(nb)):
        for v in nb[u]:
            if u > v:
                continue
            prev, todo = {u: None}, [u]   # breadth-first from u to v without the bond u-v
            while todo and v not in prev:
                nxt = []
                for x in todo:
                    for y in nb[x]:
                        if y not in prev and not (x == u and y == v):
                            prev[y] = x
                            nxt.append(y)
                todo = nxt
                if len(prev) > 4000:
                    break
            if v in prev:
                path, x = [], v
                while x is not None:
                    path.append(x)
                    x = prev[x]
                if len(path) <= limit:
                    out.add(frozenset(path))
    return out


def cut_fragment(el, pos, nb, rings, b, c, radius, box=None):
    """the QM model of rotor b-c: the atoms within `radius` bonds, grown until no small ring is cut in two and no bond
    between two unsaturated atoms (fewer neighbours than their valence: a double, triple or aromatic bond) is cut, with
    every hydrogen of a kept atom; each remaining cut bond (a single bond) capped by H along the bond. A cut through a
    ring or a double bond puts two caps on one spot or leaves a radical, which no QM code agrees on. box (the cell
    edges): the fragment is made whole round b, each atom at the image nearest the atom it was reached from (a
    crosslinked network loops through the cell, so unwrapped positions leave some bonds a box length long).
    Returns keep, idx, elements, positions, bonds and a problem string ("" when the fragment is sound)."""
    unsat = lambda x: VALENCE.get(el[x], 0) > len(nb[x])
    dist = {b: 0, c: 0}
    todo = [b, c]
    while todo:
        x = todo.pop(0)
        for y in nb[x]:
            if y not in dist and dist[x] < radius:
                dist[y] = dist[x] + 1
                todo.append(y)
    keep, problem, grown = set(dist), "", True
    while grown and len(keep) <= MAX_FRAGMENT:
        grown = False
        for r in rings:
            if r & keep and not r <= keep:
                keep |= r
                grown = True
        for x in list(keep):
            for y in nb[x]:
                # a bond between two unsaturated atoms, or an atom two kept atoms bond to (two caps on one spot)
                if y not in keep and el[y] != "H" and ((unsat(x) and unsat(y)) or sum(z in keep for z in nb[y]) > 1):
                    keep.add(y)
                    grown = True
    if len(keep) > MAX_FRAGMENT:
        problem = f"over {MAX_FRAGMENT} atoms"
    keep |= {y for x in list(keep) for y in nb[x] if el[y] == "H"}
    whole, todo = {b: list(pos[b])}, [b]   # the fragment made whole round b
    while todo:
        x = todo.pop()
        for y in nb[x]:
            if y in keep and y not in whole:
                d = min_image([pos[y][k] - pos[x][k] for k in range(3)], box)
                whole[y] = [whole[x][k] + d[k] for k in range(3)]
                todo.append(y)
    keep = sorted(whole)
    idx = {a: k for k, a in enumerate(keep)}
    fel, fpos, fbonds = [el[a] for a in keep], [whole[a] for a in keep], []
    for a in keep:
        for y in nb[a]:
            if y in idx and idx[a] < idx[y]:
                fbonds.append((idx[a], idx[y]))
            elif y not in idx:
                v = min_image([pos[y][k] - pos[a][k] for k in range(3)], box)
                n = math.sqrt(sum(t * t for t in v))
                fel.append("H")
                fpos.append([whole[a][k] + v[k] / n * CAP_BOND.get(el[a], 1.09) for k in range(3)])
                fbonds.append((idx[a], len(fel) - 1))
    if not problem and sum(Z.get(x, 0) for x in fel) % 2:
        problem = "odd electron count"
    if not problem:
        close = min((math.dist(fpos[i], fpos[j]) for i in range(len(fpos)) for j in range(i + 1, len(fpos))), default=9)
        if close < 0.9:
            problem = f"two atoms {close:.2f} A apart"
    if not problem:
        crushed = smallest_sp3_angle(fel, fpos, fbonds)
        if crushed < 85:
            problem = f"an sp3 carbon angle of {crushed:.0f} deg"
    return keep, idx, fel, fpos, fbonds, problem


def smallest_sp3_angle(fel, fpos, fbonds):
    """the smallest bond angle at a four-coordinate carbon, three-membered rings left out (an epoxide's 60 degrees are
    real): under about 85 degrees the carbon is inverted or crushed, a structure no cure should give"""
    nb = [set() for _ in fel]
    for i, j in fbonds:
        nb[i].add(j)
        nb[j].add(i)
    worst = 180.0
    for c, e in enumerate(fel):
        if e != "C" or len(nb[c]) != 4:
            continue
        ring = sorted(nb[c])
        for x in range(4):
            for y in range(x + 1, 4):
                a, b = ring[x], ring[y]
                if b in nb[a]:
                    continue
                u = [fpos[a][k] - fpos[c][k] for k in range(3)]
                v = [fpos[b][k] - fpos[c][k] for k in range(3)]
                cosang = sum(u[k] * v[k] for k in range(3)) / (math.dist(fpos[a], fpos[c]) * math.dist(fpos[b], fpos[c]))
                worst = min(worst, math.degrees(math.acos(max(-1.0, min(1.0, cosang)))))
    return worst


def environment(types, nb, rings, b, c, radius):
    """the hash of rotor b-c's environment: the atoms within `radius` bonds of b or c, labelled by type, ring membership
    and the number of bonds the cut leaves (caps), refined over the subgraph (Weisfeiler-Lehman); b and c marked"""
    dist = {b: 0, c: 0}
    todo = [b, c]
    while todo:
        x = todo.pop(0)
        for y in nb[x]:
            if y not in dist and dist[x] < radius:
                dist[y] = dist[x] + 1
                todo.append(y)
    sub = {x: [y for y in nb[x] if y in dist] for x in dist}
    lab = {x: f"{'*' if x in (b, c) else ''}{types[x]}{'@r' if x in rings else ''}+{len(nb[x]) - len(sub[x])}" for x in dist}
    for _ in range(radius + 1):
        lab = {x: hashlib.sha1((lab[x] + "|" + ",".join(sorted(lab[y] for y in sub[x]))).encode()).hexdigest()[:16]
               for x in dist}
    key = "|".join(sorted(lab.values())) + "#" + "|".join(sorted((lab[b], lab[c])))
    return hashlib.sha1(key.encode()).hexdigest()[:12]


def census(ff="oplsaa2024", out=None, only="", react=True, copolymers=True, radius=2, log=print):
    out = out or os.path.join(ROOT, "bench", "ff", "census", ff)
    os.makedirs(os.path.join(out, "fragments"), exist_ok=True)
    work = tempfile.mkdtemp(prefix="caps-census-")
    full, pub = ff_paths(ff, work)
    cl = classes(full)
    tclass = {}        # torsion class -> status, count, structures
    envs = {}          # environment hash -> record
    other = collections.defaultdict(lambda: {"missing": set(), "CAPS fit": set()})   # bond and angle classes
    rows = []
    for kind, name, build in structures(react, copolymers):
        if only and only not in name:
            continue
        try:
            doc = build()
            r_pub = doc.field.assign(pub)
            r_full = doc.field.assign(full)
            data = os.path.join(work, "s.data")
            doc.save(data)
            types_file, masses, pos, bonds, _ = read_data(data)
            box = read_box(data)
        except Exception as e:   # the build, the cure or the assignment failed: listed, not counted
            rows.append({"kind": kind, "name": name, "error": str(e)[:300]})
            log(f"FAIL {kind}: {name}: {str(e)[:120]}")
            continue
        atoms = r_full.get("atoms", [])
        types = [a.get("type") or "" for a in atoms] if len(atoms) == len(masses) else types_file
        # the element of each typed atom (a data file's mass can stand for a type that is not the atom's element)
        el = ([a.get("el") or element_of_mass(m) for a, m in zip(atoms, masses)] if len(atoms) == len(masses)
              else [element_of_mass(m) for m in masses])
        nb = neighbours(len(el), bonds)
        rings = ring_atoms(nb)
        small = small_rings(nb)
        miss_full, miss_pub = missing_terms(r_full, cl), missing_terms(r_pub, cl)
        for k in ("bond", "angle"):
            other[k]["missing"] |= miss_full[k]
            other[k]["CAPS fit"] |= miss_pub[k] - miss_full[k]

        def status(key):
            return "missing" if key in miss_full["dihedral"] else "CAPS fit" if key in miss_pub["dihedral"] else "published"

        by_bond = collections.defaultdict(set)
        n_tors = 0
        for q in propers(nb):
            if any(not types[x] for x in q):
                continue
            key = tuple(cl.get(types[x], types[x]) for x in q)
            key = min(key, key[::-1])
            t = tclass.setdefault(key, {"status": status(key), "count": 0, "structures": set()})
            t["count"] += 1
            t["structures"].add(name)
            by_bond[tuple(sorted(q[1:3]))].add(key)
            n_tors += 1
        n_rot = 0
        for (b, c), keys in by_bond.items():
            if in_ring(nb, b, c) or len(nb[b]) < 2 or len(nb[c]) < 2 or linear(pos, nb, b, box) or linear(pos, nb, c, box):
                continue
            n_rot += 1
            h = environment(types, nb, rings, b, c, radius)
            st = max((status(k) for k in keys), key=RANK.get)
            e = envs.get(h)
            if e is None:
                keep, idx, fel, fpos, fbonds, problem = cut_fragment(el, pos, nb, small, b, c, radius, box)
                central = [idx[b] + 1, idx[c] + 1]   # the rotatable bond in the fragment's numbering
                write_mol2(os.path.join(out, "fragments", f"{h}.mol2"), fel, fpos, fbonds,
                           name=f"env {h} bond {central[0]}-{central[1]} from {name}")
                try:
                    smi = doc.fragment_smiles([x for x in keep if el[x] != "H"])   # 0-based; it adds the H itself
                except Exception:
                    smi = ""
                heavy = (any(el[x] != "H" for x in nb[b] if x != c) and any(el[x] != "H" for x in nb[c] if x != b))
                e = envs[h] = {"env": h, "status": st, "central": [types[b], types[c]], "central_atoms": central,
                               "classes": sorted("-".join(k) for k in keys), "count": 0, "structures": set(),
                               "example": name, "fragment_smiles": smi, "fragment_atoms": len(fel),
                               "fragment_heavy": sum(1 for x in fel if x != "H"), "heavy_rotor": heavy,
                               "elements": sorted(set(fel)), "problem": problem}
            e["count"] += 1
            e["structures"].add(name)
            if RANK[st] > RANK[e["status"]]:
                e["status"] = st
        rows.append({"kind": kind, "name": name, "atoms": len(el), "untyped": int(r_full.get("untyped", 0)),
                     "torsions": n_tors, "rotors": n_rot, "complete": bool(r_full.get("complete"))})
        log(f"ok   {kind}: {name}: {len(el)} atoms, {n_tors} torsions, {n_rot} rotors")
    return write(out, ff, full, radius, rows, tclass, envs, other)


# ------------------------------------------------------------------ report
def write(out, ff, full, radius, rows, tclass, envs, other):
    ok = [r for r in rows if "error" not in r]
    cs = collections.Counter(t["status"] for t in tclass.values())
    es = collections.Counter(e["status"] for e in envs.values())
    eh = collections.Counter(e["status"] for e in envs.values() if e["heavy_rotor"])
    by_kind = collections.Counter(r["kind"] for r in ok)
    heavy = sorted(e["fragment_heavy"] for e in envs.values())
    med = heavy[len(heavy) // 2] if heavy else 0
    need = sum(1 for e in envs.values() if e["status"] != "published")
    res = {
        "forcefield": ff, "file": os.path.relpath(full, ROOT), "radius": radius, "structures": rows,
        "torsion_classes": [{"class": "-".join(k), "status": v["status"], "count": v["count"],
                             "structures": sorted(v["structures"])[:5]} for k, v in sorted(tclass.items())],
        "environments": sorted(({**e, "structures": sorted(e["structures"])[:5]} for e in envs.values()),
                               key=lambda e: (-RANK[e["status"]], -e["count"], e["env"])),
        "bonds_angles": {k: {s: sorted("-".join(x) for x in v[s]) for s in v} for k, v in other.items()},
    }
    json.dump(res, open(os.path.join(out, "census.json"), "w"), indent=1)
    L = [f"# Torsion census: {ff}", "",
         f"Structures typed: {len(ok)} ({', '.join(f'{n} {k}' for k, n in sorted(by_kind.items()))}); "
         f"failed: {len(rows) - len(ok)}. Environment radius: {radius} bonds.", "",
         "| | published | CAPS fit | missing | total |", "|---|---|---|---|---|",
         f"| torsion classes | {cs['published']} | {cs['CAPS fit']} | {cs['missing']} | {len(tclass)} |",
         f"| rotor environments | {es['published']} | {es['CAPS fit']} | {es['missing']} | {len(envs)} |",
         f"| of them heavy-atom rotors | {eh['published']} | {eh['CAPS fit']} | {eh['missing']} | {sum(eh.values())} |",
         *(f"| {k} classes | - | {len(other[k]['CAPS fit'])} | {len(other[k]['missing'])} | - |" for k in ("bond", "angle")),
         "", "## Campaign size", "",
         f"- One relaxed scan per rotor environment, {SCAN_POINTS} points each (15 degrees).",
         f"- Library (environments the published table does not cover): {need} scans, {need * SCAN_POINTS} QM points.",
         f"- Training set (every environment, published ones included): {len(envs)} scans, "
         f"{len(envs) * SCAN_POINTS} QM points before distorted conformers are added.",
         f"- Fragment size: median {med} heavy atoms, largest {heavy[-1] if heavy else 0}.",
         f"- Fragments with a problem (left out of a campaign): {sum(1 for e in envs.values() if e['problem'])} "
         f"({', '.join(f'{n} {p}' for p, n in collections.Counter(e['problem'].split(' (')[0] for e in envs.values() if e['problem']).most_common())}).",
         f"- Environments seen once: {sum(1 for e in envs.values() if e['count'] == 1)} (many are chain ends or one "
         f"crosslink site; the ones seen often are the backbones).",
         "", "## Missing torsion classes", ""]
    L += [f"- {'-'.join(k)} ({v['count']}x; {', '.join(sorted(v['structures'])[:3])})"
          for k, v in sorted(tclass.items()) if v["status"] == "missing"] or ["- none"]
    L += ["", "## Failed structures", ""]
    L += [f"- {r['kind']}: {r['name']}: {r['error']}" for r in rows if "error" in r] or ["- none"]
    open(os.path.join(out, "summary.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:18]))
    return res


def main():
    census(ff=arg("--ff", "oplsaa2024"), out=arg("--out"), only=arg("--only", ""), react="--no-react" not in sys.argv,
           copolymers="--no-copolymers" not in sys.argv, radius=int(arg("--radius", "2")))


if __name__ == "__main__":
    main()
