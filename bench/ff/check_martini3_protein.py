#!/usr/bin/env python3
"""CAPS's Martini 3 proteins against martinize2 (vermouth's own Martini 3 integration tests: the all-atom input,
martinize2's command, topology and coarse-grained coordinates), term by term: bead types and charges, bonds (the
rubber band included), constraints, angles, dihedrals, impropers, virtual sites, exclusions, and each bead's position.

The reference .itp is read on its default path (#ifndef FLEXIBLE: constraints, not the FLEXIBLE bonds). CAPS writes its
constraints as stiff bonds and lists them back as constraints. DSSP is compared on its own (the references were made
with the DSSP their vermouth version called); the topology is built from the reference's secondary structure, so the
comparison is of the model, not of DSSP.

usage: check_martini3_protein.py [REF_DIR]   (default ~/vermouth-ref/tests-m3, the downloaded tests)
"""
import math, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
REF = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/vermouth-ref/tests-m3")
CASES = ["tier-0/dipro-termini", "tier-1/hst5", "tier-1/1UBQ", "tier-1/lysozyme", "tier-1/EN_chain"]


def caps_args(command):
    """martinize2's options as CAPS's"""
    w = command.split()
    out, ss, i = [], None, 0
    while i < len(w):
        a = w[i]
        nxt = w[i + 1] if i + 1 < len(w) else None
        if a == "-dssp":
            ss = "dssp"
        elif a == "-ss":
            ss = nxt
            i += 1
        elif a == "-nt":
            out += ["--nt"]
        elif a == "-noscfix":
            out += ["--noscfix"]
        elif a == "-elastic":
            out += ["--elastic"]
        elif a in ("-ef", "-el", "-eu", "-ea", "-ep", "-em", "-ermd", "-eunit"):
            out += ["--" + a[1:], nxt]
            i += 1
        elif a == "-id-regions":
            out += ["--idr", nxt]
            i += 1
        elif a == "-cys":
            out += ["--cys", nxt]
            i += 1
        elif a in ("-f", "-x", "-o", "-ff", "-ignore", "-p", "-maxwarn"):
            i += 1
        i += 1
    return out, ss


def read_itp(path):
    sec, skip, out = None, [], {k: [] for k in ("atoms", "bonds", "constraints", "angles", "dihedrals", "impropers", "virtual_sitesn", "exclusions")}
    header_ss = None
    lines = open(path).read().splitlines()
    for k, raw in enumerate(lines):
        if raw.startswith("; was used for the full system") or raw.startswith("; The following sequence of secondary"):
            for nxt in lines[k + 1:k + 3]:
                m = re.match(r"^;\s*([A-Z~ ]+)\s*$", nxt)
                if m and len(m.group(1).strip()) > 1:
                    header_ss = m.group(1).strip()
                    break
        l = raw.split(";")[0].strip()
        if l.startswith("#ifdef"):
            skip.append("FLEXIBLE" in l)
            continue
        if l.startswith("#ifndef"):
            skip.append(False)
            continue
        if l.startswith("#else"):
            skip[-1] = not skip[-1]
            continue
        if l.startswith("#endif"):
            skip.pop()
            continue
        if any(skip) or not l:
            continue
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        w = l.split()
        if sec == "atoms":
            out["atoms"].append((w[1], w[3], w[4], round(float(w[6]) if len(w) > 6 else 0.0, 3)))
        elif sec == "bonds":
            i, j = sorted((int(w[0]), int(w[1])))
            out["bonds"].append((i, j, round(float(w[3]), 3), round(float(w[4]), 1)))
        elif sec == "constraints":
            i, j = sorted((int(w[0]), int(w[1])))
            out["constraints"].append((i, j, round(float(w[3]), 3)))
        elif sec == "angles":
            out["angles"].append((int(w[0]), int(w[1]), int(w[2]), int(w[3]), round(float(w[4]), 1), round(float(w[5]), 1)))
        elif sec == "dihedrals":
            key = tuple(int(x) for x in w[:4])
            if w[4] == "2":
                out["impropers"].append(key + (round(float(w[5]), 1), round(float(w[6]), 1)))
            else:
                out["dihedrals"].append(key + (int(w[4]), round(float(w[5]), 1), round(float(w[6]), 1), int(w[7]) if len(w) > 7 else 1))
        elif sec == "virtual_sitesn":
            if "--" in w:   # martinize writes "site func atoms"
                w = [x for x in w if x != "--"]
            out["virtual_sitesn"].append((int(w[0]), int(w[1]), tuple(sorted(int(x) for x in w[2:]))))
        elif sec == "exclusions":
            a = int(w[0])
            for b in w[1:]:
                out["exclusions"].append(tuple(sorted((a, int(b)))))
    return out, header_ss


def read_pdb(path):
    return [tuple(float(l[30 + 8 * k:38 + 8 * k]) for k in range(3)) for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def multiset_diff(a, b):
    from collections import Counter
    ca, cb = Counter(a), Counter(b)
    return sorted((ca - cb).elements()), sorted((cb - ca).elements())


fails = 0
for case in CASES:
    d = os.path.join(REF, case)
    cmd = open(os.path.join(d, "martinize2", "command")).read()
    args, ss = caps_args(cmd)
    ref, ref_ss = read_itp(os.path.join(d, "martinize2", "molecule_0.itp"))
    work = tempfile.mkdtemp()
    print(f"{case}: martinize2 {' '.join(cmd.split()[1:])}")
    ok = True
    if ss == "dssp":
        r = subprocess.run([CAPS, "dssp", os.path.join(d, "aa.pdb")], capture_output=True, text=True)
        mine = r.stdout.split("\n")[0].strip()
        if ref_ss:
            same = mine == ref_ss
            diff = [k + 1 for k, (x, y) in enumerate(zip(mine, ref_ss)) if x != y]
            print(f"   DSSP: {'identical to the reference' if same else f'differs ({len(mine)} / {len(ref_ss)} letters) at ' + str(diff[:12])}")
            if not same:
                print(f"      CAPS      {mine}\n      reference {ref_ss}")
        ss_arg = ["--ss", ref_ss or mine]
    elif ss is None:
        ss_arg = ["--ss", "none"]
    else:
        ss_arg = ["--ss", ss]
    r = subprocess.run([CAPS, "martini", os.path.join(d, "aa.pdb"), "--martini", "3", "-o", os.path.join(work, "cg.data"),
                        "--itp", os.path.join(work, "cg.itp")] + ss_arg + args, capture_output=True, text=True)
    if r.returncode:
        print(f"   CAPS failed: {r.stderr.strip()}")
        fails += 1
        continue
    mine, _ = read_itp(os.path.join(work, "cg.itp"))
    ta = [(a[0], a[3]) for a in ref["atoms"]]
    tb = [(a[0], a[3]) for a in mine["atoms"]]
    bad = [k + 1 for k, (x, y) in enumerate(zip(ta, tb)) if x != y]
    same = len(ta) == len(tb) and not bad
    ok &= same
    print(f"   beads: {len(tb)} / {len(ta)}, types and charges {'identical' if same else 'differ at ' + str(bad[:8]) + ' ' + str([(ta[k-1], tb[k-1]) for k in bad[:4]])}")
    for sec in ("bonds", "constraints", "angles", "dihedrals", "impropers", "virtual_sitesn", "exclusions"):
        missing, extra = multiset_diff(ref[sec], mine[sec])
        same = not missing and not extra
        ok &= same
        if ref[sec] or mine[sec]:
            print(f"   {sec}: {len(mine[sec])} / {len(ref[sec])} {'identical' if same else f'missing {missing[:5]} extra {extra[:5]}'}")
    pr = read_pdb(os.path.join(d, "martinize2", "cg.pdb"))
    pc, sec, box = [], None, [0.0, 0.0, 0.0]
    for l in open(os.path.join(work, "cg.data")):
        s = l.strip()
        for k, ax in enumerate("xyz"):
            if s.endswith(f"{ax}lo {ax}hi"):
                lo, hi = map(float, s.split()[:2])
                box[k] = hi - lo
        if s.startswith("Atoms"):
            sec = "A"
            continue
        if s and s[0].isalpha():
            sec = None
        w = s.split()
        if sec == "A" and len(w) >= 7:
            img = [int(x) for x in w[7:10]] if len(w) >= 10 else [0, 0, 0]   # unwrap: the data file keeps image flags
            pc.append((int(w[0]), tuple(float(w[4 + k]) + img[k] * box[k] for k in range(3))))
    pc = [p for _, p in sorted(pc)]
    dev = max(math.dist(x, y) for x, y in zip(pr, pc)) if len(pr) == len(pc) else float("inf")
    print(f"   bead positions: largest difference {dev:.3f} A from martinize2's cg.pdb")
    ok &= dev < 0.01
    fails += 0 if ok else 1
print(f"{len(CASES) - fails} of {len(CASES)} Martini 3 proteins identical to martinize2")
sys.exit(1 if fails else 0)
