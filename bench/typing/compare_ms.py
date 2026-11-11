#!/usr/bin/env python3
"""PCFF, COMPASS and CVFF atom types (and force-field charges) from CAPS against Materials Studio's own typing.

The molecules are those of bench/typing/gaff_reference_set.smi, built by CAPS (compare_antechamber.py prepare writes
NAME.caps.mol2 and NAME.mol2). Materials Studio 23.1 types each one with Forcite (bench/typing/ms_types.pl, run by
RunMatScript) and writes NAME.ms_FF.types, "type charge" per atom in file order.

  compare_ms.py MOLDIR MSDIR [--show 60]

MOLDIR holds NAME.caps.mol2; MSDIR the NAME.ms_FF.types files. CAPS types each molecule with pcff-frc, compass-frc and
cvff and the two are compared atom by atom (types; charges within 0.001 e). COMPASS is compared with both COMPASS III and
COMPASS26: the COMPASS CAPS ships is the published 1998 subset, whose type names the later versions keep.
"""
import argparse, os, subprocess

PAIRS = [("pcff-frc", "pcff"), ("cvff", "cvff"), ("compass-frc", "COMPASSIII"), ("compass-frc", "COMPASS26")]


def ms_types(p):
    t, q = [], []
    for l in open(p):
        w = l.split()
        if w:
            t.append(w[0]); q.append(float(w[1]) if len(w) > 1 else float("nan"))
    return t, q


def caps_types(caps, mol2, ff, out):
    r = subprocess.run([caps, "ff", "type", mol2, "--ff", ff, "-o", out], capture_output=True, text=True)
    if not os.path.exists(out):
        return None
    return [w.split()[-1] for w in open(out) if w.strip() and not w.startswith("#")]


def caps_charges(caps, mol2, ff, types, out):
    """CAPS's force-field charges (bond increments) for the same molecule"""
    r = subprocess.run([caps, "ff", "apply", mol2, "--ff", ff, "--types", types, "--charges", "auto", "-o", out],
                       capture_output=True, text=True)
    if r.returncode or not os.path.exists(out):
        return None
    q, sec = [], False
    for l in open(out):
        if l.split("#")[0].strip() == "Atoms":
            sec = True; continue
        if sec:
            w = l.split()
            if not w:
                if q: break
                continue
            q.append((int(w[0]), float(w[3])))
    return [c for _, c in sorted(q)]


ap = argparse.ArgumentParser()
ap.add_argument("moldir"); ap.add_argument("msdir")
ap.add_argument("--caps", default="build/cli/caps"); ap.add_argument("--show", type=int, default=60)
a = ap.parse_args()
work = os.path.join(a.msdir, "caps")
os.makedirs(work, exist_ok=True)
names = sorted(f[:-10] for f in os.listdir(a.moldir) if f.endswith(".caps.mol2"))
for capsff, msff in PAIRS:
    ok = tot = untyped = qok = qtot = 0
    full, diffs, ms_missing = 0, [], 0
    for n in names:
        msp = os.path.join(a.msdir, f"{n}.ms_{msff}.types")
        if not os.path.exists(msp):
            ms_missing += 1
            continue
        mt, mq = ms_types(msp)
        mol2 = os.path.join(a.moldir, n + ".caps.mol2")
        tf = os.path.join(work, f"{n}.{capsff}.types")
        ct = caps_types(a.caps, mol2, capsff, tf)
        if ct is None or len(ct) != len(mt):
            diffs.append((n, "-", "CAPS failed" if ct is None else f"{len(ct)} vs {len(mt)} atoms", ""))
            continue
        cq = caps_charges(a.caps, mol2, capsff, tf, os.path.join(work, f"{n}.{capsff}.data")) if "?" not in ct else None
        allok = True
        for i, (x, y) in enumerate(zip(mt, ct)):
            tot += 1
            if y in ("?", "(untyped)"):
                untyped += 1; allok = False
                diffs.append((n, i + 1, x, "untyped"))
            elif x.lower() == y.lower():
                ok += 1
            else:
                allok = False
                diffs.append((n, i + 1, x, y))
            if cq and i < len(cq) and mq[i] == mq[i]:
                qtot += 1
                qok += abs(cq[i] - mq[i]) < 1e-3
        full += allok
    print(f"{capsff} vs Materials Studio {msff}: {ok}/{tot} atoms ({100 * ok / max(tot, 1):.1f} %), {untyped} untyped by CAPS, "
          f"{full}/{len(names) - ms_missing} molecules fully right; charges {qok}/{qtot} within 0.001 e")
    for d in diffs[:a.show]:
        print(f"   {d[0]:26s} atom {str(d[1]):4s} MS {d[2]!s:6s} CAPS {d[3]}")
