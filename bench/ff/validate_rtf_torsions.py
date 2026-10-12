#!/usr/bin/env python3
"""Check which torsions CAPS gives a force field that lists one torsion per bond (CHARMM19's toph19: "a single
dihedral per bond rather than multiple torsions"; "torsion_selection": "one_per_bond") against the DIHE records of the
force field's own topology.

Each residue is built and typed as validate_rtf_types.py does (its mol2 files, RTF_KEEP); CAPS then assigns the force
field (caps ff apply) and the torsions of its GROMACS topology are read back. For every central bond inside the residue
the topology's torsions and CAPS's are compared: the same bonds (none where the topology lists none), the same number
about each bond, and the same four atoms (an atom of a neighbouring residue counts as its link name, -C or +N).

usage: validate_rtf_torsions.py RTF FF.json RULES.json [validate_rtf_types.py options: --united --cap-links ...] [--all]
"""
import collections, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CAPS = os.path.join(os.path.dirname(os.path.dirname(HERE)), "build", "cli", "caps")
rtf, ffjson, rules = sys.argv[1:4]
opts = [a for a in sys.argv[4:] if a != "--all"]
work = tempfile.mkdtemp()
subprocess.run([sys.executable, os.path.join(HERE, "validate_rtf_types.py"), rtf, ffjson, rules] + opts,
               env={**os.environ, "RTF_KEEP": work}, capture_output=True, text=True)

# DIHE records and atom names per residue (RESI only; patches aside)
res, cur = {}, None
for raw in open(rtf, errors="replace"):
    w = raw.split("!")[0].split()
    if not w:
        continue
    k = w[0].upper()[:4]
    if k in ("RESI", "PRES"):
        cur = {"atoms": [], "dihe": []} if k == "RESI" else None
        if cur is not None:
            res[w[1]] = cur
    elif cur is not None and k == "ATOM":
        cur["atoms"].append(w[1].upper())
    elif cur is not None and k in ("DIHE", "DIHED"):
        for q in range(1, len(w) - 3, 4):
            cur["dihe"].append(tuple(x.upper() for x in w[q:q + 4]))

import fnmatch, json
DIH = json.load(open(ffjson))["dihedrals"]


def zero_term(tq):
    """the force field's term for these four types (None: an atom of the neighbouring residue) is zero: CAPS writes no
    torsion for it, the energy is the same. The most specific matching entry, as CHARMM chooses (fewest wildcards)."""
    best = None
    for d in DIH:
        m = d["match"]
        for mm in (m, m[::-1]):
            if all(t is None and x == "*" or t is not None and fnmatch.fnmatchcase(t, x) for t, x in zip(tq, mm)) or \
               all(t is None or fnmatch.fnmatchcase(t, x) for t, x in zip(tq, mm)) and tq[0] is None and tq[3] is None:
                w = sum(x == "*" for x in m)
                if best is None or w < best[0]:
                    best = (w, d)
    if best is None:
        return False
    p = best[1]["params"]
    return all(p[1 + 3 * q] == 0 for q in range(int(p[0])))


zero = []
bonds_right = bonds_total = n_right = n_total = quad_right = quad_total = 0
extra, short, long_, other_ends, missing_bond = collections.Counter(), [], [], [], []
for name, r in res.items():
    mol2 = os.path.join(work, re.sub(r"[^A-Za-z0-9_]", "_", name) + ".mol2")
    if not os.path.exists(mol2) or not os.path.exists(mol2 + ".types"):
        continue
    names = r["atoms"]
    stem = mol2 + ".gmx"
    p = subprocess.run([CAPS, "ff", "apply", mol2, "--ff", ffjson, "--typing", rules, "--charges", "gasteiger", "--allow-missing",
                        "--gromacs", stem], capture_output=True, text=True)
    if not os.path.exists(stem + ".itp"):
        continue
    sec, caps = None, set()
    for l in open(stem + ".itp"):   # proper torsions: [ dihedrals ] functions 1, 3, 9 (impropers 2, 4)
        s = l.split(";")[0].strip()
        m = re.match(r"\[\s*(\S+)\s*\]", s)
        if m:
            sec = m.group(1)
            continue
        w = s.split()
        if sec == "dihedrals" and len(w) >= 5 and w[4] in ("1", "3", "9"):
            caps.add(tuple(int(x) - 1 for x in w[:4]))
    label = lambda i: names[i] if i < len(names) else "*ext"
    own = {a for a in names}
    ref = collections.defaultdict(set)
    for q in r["dihe"]:
        if q[1] in own and q[2] in own:
            ref[frozenset(q[1:3])].add(tuple("*ext" if x[:1] in "+-" else x for x in q))
    got = collections.defaultdict(set)
    for q in caps:
        if q[1] < len(names) and q[2] < len(names):
            got[frozenset((names[q[1]], names[q[2]]))].add(tuple(label(i) for i in q))
    for b in set(ref) | set(got):
        bonds_total += 1
        a, c = ref.get(b, set()), got.get(b, set())
        if not a:
            extra[name] += 1
            missing_bond.append((name, "CAPS only", sorted(c)))
            continue
        if not c:
            ty = [x.strip() for x in open(mol2 + ".types")]
            if all(zero_term(tuple(ty[names.index(x)] if x in names else None for x in q)) for q in a):
                zero.append((name, sorted(a)))
                bonds_right += 1
                continue
            missing_bond.append((name, "topology only", sorted(a)))
            continue
        bonds_right += 1
        n_total += 1
        canon = lambda s: {q if q <= q[::-1] else q[::-1] for q in s}
        if len(a) == len(c):
            n_right += 1
        else:
            (short if len(c) < len(a) else long_).append((name, sorted(a), sorted(c)))
        quad_total += len(a)
        same = canon(a) & canon(c)
        quad_right += len(same)
        if len(same) < len(a):
            other_ends.append((name, sorted(canon(a) - same), sorted(canon(c) - same)))
print(f"{os.path.basename(rtf)}: central bonds {bonds_right}/{bonds_total} the same (a torsion in both, or a zero term), "
      f"the number of torsions about them {n_right}/{n_total}, the same four atoms {quad_right}/{quad_total}")
if zero:
    print(f"  {len(zero)} bonds whose listed torsion has a zero term (CAPS writes none; the same energy): "
          + ", ".join(sorted({'-'.join(q[0]) for _, q in zero})))
for name, how, q in missing_bond:
    print(f"  {name}: {how}: " + ", ".join("-".join(x) for x in q))
for name, a, c in short + long_:
    print(f"  {name}: topology {len(a)} ({', '.join('-'.join(x) for x in a)}) / CAPS {len(c)}")
if "--all" in sys.argv:
    for name, a, c in other_ends:
        print(f"  {name}: other ends: topology {', '.join('-'.join(x) for x in a)} / CAPS {', '.join('-'.join(x) for x in c)}")
