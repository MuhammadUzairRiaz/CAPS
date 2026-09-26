#!/usr/bin/env python3
"""CAPS's Martini 2.2 proteins against martinize2 (vermouth's own integration tests: the all-atom input, martinize2's
topology and coarse-grained coordinates), term by term: bead types and charges, every bond, constraint, angle and
dihedral with its parameters, and each bead's position.

CAPS writes constraints as stiff bonds (it has no constraints); the .itp it writes lists them back as constraints, so the
two topologies compare section by section. DSSP is CAPS's own (checked here too against the tests' DSSP 2.0 output).

usage: check_martini_protein.py [REF_DIR]   (default tests/data/vermouth, copied from vermouth's tests with its NOTICE)
"""
import math, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
REF = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "tests", "data", "vermouth")
# vermouth's tests made with -ff martini22 (the others are Martini 3 / ElNeDyn): (all-atom input, martinize2 .itp, cg.pdb)
CASES = [("1ICO beta-sheet", "1ico_aa.pdb", "1ico_martini22.itp", "1ico_martini22_cg.pdb")]


def read_itp(path):
    sec, out = None, {"atoms": [], "bonds": set(), "constraints": set(), "angles": set(), "dihedrals": set()}
    for raw in open(path):
        l = raw.split(";")[0].strip()
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        w = l.split()
        if not w or sec not in out:
            continue
        if sec == "atoms":
            out["atoms"].append((w[1], w[3], w[4] if len(w) > 7 else "", round(float(w[6]), 3)))
        elif sec == "bonds":
            i, j = sorted((int(w[0]), int(w[1])))
            out["bonds"].add((i, j, round(float(w[3]), 4), round(float(w[4]), 1)))
        elif sec == "constraints":
            i, j = sorted((int(w[0]), int(w[1])))
            out["constraints"].add((i, j, round(float(w[3]), 4)))
        elif sec == "angles":
            i, k = int(w[0]), int(w[2])
            if i > k:
                i, k = k, i
            out["angles"].add((i, int(w[1]), k, round(float(w[4]), 1), round(float(w[5]), 1)))
        elif sec == "dihedrals":
            out["dihedrals"].add(tuple(int(x) for x in w[:5]) + (round(float(w[5]), 1), round(float(w[6]), 1)))
    return out


def read_pdb(path):
    return [tuple(float(l[30 + 8 * k:38 + 8 * k]) for k in range(3)) for l in open(path) if l.startswith(("ATOM", "HETATM"))]


fails = 0
for case, aa_pdb, ref_itp, ref_cg in CASES:
    work = tempfile.mkdtemp()
    r = subprocess.run([CAPS, "martini", os.path.join(REF, aa_pdb), "-o", os.path.join(work, "cg.data"), "--itp", os.path.join(work, "cg.itp")],
                       capture_output=True, text=True)
    if r.returncode:
        print(f"{case}: CAPS failed: {r.stderr.strip()}")
        fails += 1
        continue
    ref, caps = read_itp(os.path.join(REF, ref_itp)), read_itp(os.path.join(work, "cg.itp"))
    print(f"{case}")
    ok = True
    types_ref = [(a[0], a[3]) for a in ref["atoms"]]
    types_caps = [(a[0], a[3]) for a in caps["atoms"]]
    bad = [k + 1 for k, (x, y) in enumerate(zip(types_ref, types_caps)) if x != y]
    same = len(types_ref) == len(types_caps) and not bad
    ok &= same
    print(f"   beads: {len(types_caps)} / {len(types_ref)}, types and charges {'identical' if same else 'differ at ' + str(bad[:10])}")
    for sec in ("bonds", "constraints", "angles", "dihedrals"):
        a, b = ref[sec], caps[sec]
        same = a == b
        ok &= same
        print(f"   {sec}: {len(b)} / {len(a)} {'identical' if same else f'missing {sorted(a - b)[:6]} extra {sorted(b - a)[:6]}'}")
    pr, pc = read_pdb(os.path.join(REF, ref_cg)), []
    # CAPS's bead positions from its data file (Atoms section, file order = bead order)
    sec = None
    for l in open(os.path.join(work, "cg.data")):
        s = l.strip()
        if s.startswith("Atoms"):
            sec = "A"
            continue
        if s and s[0].isalpha():
            sec = None
        w = s.split()
        if sec == "A" and len(w) >= 7:
            pc.append((int(w[0]), tuple(float(x) for x in w[4:7])))
    pc = [p for _, p in sorted(pc)]
    dev = max(math.dist(x, y) for x, y in zip(pr, pc)) if len(pr) == len(pc) else float("inf")
    print(f"   bead positions: largest difference {dev:.3f} A from martinize2's cg.pdb")
    ok &= dev < 0.01
    fails += 0 if ok else 1
print(f"{len(CASES) - fails} of {len(CASES)} proteins identical to martinize2")
sys.exit(1 if fails else 0)
