#!/usr/bin/env python3
"""Compare CAPS's automatic typing with Materials Studio's, as recorded in msi2lmp data files (the Masses comments
give each numeric type's force-field type). An independent check: these structures were typed by Accelrys'
own rules, not by DL_FIELD's templates.

usage: compare_msi_types.py FF.json RULES.json DATA [DATA ...]
"""
import collections, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
ff, rules, files = sys.argv[1], sys.argv[2], sys.argv[3:]
grand = [0, 0]
for path in files:
    text = open(path).read()
    label = {int(m.group(1)): m.group(2) for m in re.finditer(r"^\s*(\d+)\s+[\d.]+\s+#\s*(\S+)", text.split("Masses", 1)[1].split("Coeffs")[0].split("Atoms")[0], re.M)}
    atoms = []
    sec = text.split("Atoms", 1)[1].split("\n\n", 2)[1]
    for line in sec.splitlines():
        w = line.split()
        if len(w) >= 7:
            atoms.append((int(w[0]), int(w[2])))
    atoms.sort()
    expect = [label[t] for _, t in atoms]
    out = tempfile.mktemp(suffix=".types")
    r = subprocess.run([CAPS, "ff", "type", path, "--ff", ff, "--typing", rules, "-o", out], capture_output=True, text=True)
    if not os.path.exists(out):
        print(f"{os.path.basename(path)}: not typed: {(r.stderr or r.stdout).strip()[-200:]}")
        continue
    got = [x.strip() for x in open(out)]
    conf = collections.Counter((e, g) for e, g in zip(expect, got) if e != g)
    ok = sum(e == g for e, g in zip(expect, got))
    grand[0] += ok
    grand[1] += len(expect)
    print(f"{os.path.basename(path)}: {ok}/{len(expect)} atoms as Materials Studio" +
          ("" if not conf else "  differences: " + ", ".join(f"{e}→{g}×{k}" for (e, g), k in conf.most_common())))
print(f"total {grand[0]}/{grand[1]} ({100.0 * grand[0] / max(grand[1], 1):.1f} %)")
