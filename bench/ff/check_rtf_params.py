#!/usr/bin/env python3
"""Check that a CAPS CHARMM-format library parameterises every residue of its own topology file completely.

Each RESI of the .rtf (no lone pairs, no bonds to other residues) is written as a mol2 file with the topology's charges
and handed to caps ff apply, with the topology's own atom types (default) or with CAPS's automatic typing (--auto). A
residue passes when:
  * its charges sum to an integer, the one its RESI line gives;
  * every bond, angle and dihedral has a parameter of its own (none missing, none taken from an analogous type);
  * every improper the residue's IMPR records list is there (the same four atoms; one whose parameter is zero adds
    nothing and is not expected).
Impropers CAPS adds beyond those are counted apart: CAPS applies a parameter-file improper wherever its types match,
CHARMM only where the residue's IMPR records ask for one.

usage: check_rtf_params.py TOP.rtf FF.json [--auto RULES.json] [--only RESI,...] [--all]
"""
import collections, json, os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
rtf, ffjson = sys.argv[1], sys.argv[2]
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
only = [x for x in (arg("--only") or "").split(",") if x]
auto = arg("--auto")
FF = json.load(open(ffjson))
types_of = {t["name"]: t for t in FF["atom_types"]}
import fnmatch


def zero_improper(types4):
    """the library's improper for these four types (in the IMPR order or reversed) has a zero force constant: no term"""
    for e in reversed(FF.get("impropers", [])):
        for ty in (types4, types4[::-1]):
            if all(fnmatch.fnmatchcase(t, m) for t, m in zip(ty, e["match"])):
                k = e["params"][0] if e.get("style") == "harmonic" else e["params"][1]
                return k == 0
    return False


def residues(path):
    res, cur = [], None
    for line in open(path, errors="replace"):
        w = line.split("!")[0].split()
        if not w:
            continue
        k = w[0].upper()
        if k in ("RESI", "PRES"):
            cur = {"name": w[1], "charge": float(w[2]) if len(w) > 2 else 0.0, "atoms": [], "bonds": [], "impr": []} if k == "RESI" else None
            if cur is not None:
                res.append(cur)
        elif cur is None:
            continue
        elif k == "ATOM" and len(w) >= 4:
            cur["atoms"].append((w[1], w[2], float(w[3])))
        elif k in ("BOND", "DOUB", "DOUBLE", "TRIP", "TRIPLE", "AROM"):
            cur["bonds"] += list(zip(w[1::2], w[2::2]))
        elif k[:4] in ("IMPR", "IMPH"):
            cur["impr"] += [tuple(w[q:q + 4]) for q in range(1, len(w) - 3, 4)]
        elif k == "END":
            break
    return res


work = tempfile.mkdtemp()
counts, failures, extra = collections.Counter(), [], {}
for r in residues(rtf):
    if only and r["name"] not in only:
        continue
    names = [a.upper() for a, _, _ in r["atoms"]]
    if any(t.upper().startswith(("LP", "DUM")) or t not in types_of for _, t, _ in r["atoms"]):
        counts["skipped: lone pairs, dummies or types the library lacks"] += 1
        continue
    if any(x[:1] in "+-" for b in r["bonds"] for x in b) or len(names) < 2:
        counts["skipped: bonds to other residues or a single atom"] += 1
        continue
    idx = {a: i for i, a in enumerate(names)}
    bonds = sorted({tuple(sorted((idx[a.upper()], idx[b.upper()]))) for a, b in r["bonds"]})
    stem = os.path.join(work, re.sub(r"[^A-Za-z0-9_]", "_", r["name"]))
    with open(stem + ".mol2", "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{r['name']}\n{len(names)} {len(bonds)} 1\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, (a, t, q) in enumerate(r["atoms"]):
            f.write(f"{i + 1} {a} {i * 1.5:.3f} {(i % 3) * 0.7:.3f} 0 {types_of[t]['element']} 1 {r['name'][:4]} {q:.4f}\n")
        f.write("@<TRIPOS>BOND\n" + "".join(f"{k + 1} {i + 1} {j + 1} un\n" for k, (i, j) in enumerate(bonds)))
    with open(stem + ".types", "w") as f:
        f.write("".join(t + "\n" for _, t, _ in r["atoms"]))
    cmd = [CAPS, "ff", "apply", stem + ".mol2", "--ff", ffjson, "--charges", "keep", "--allow-missing", "--list"]
    cmd += ["--typing", auto] if auto else ["--types", stem + ".types"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = p.stdout
    why = []
    q = sum(x for _, _, x in r["atoms"])
    if abs(q - round(q)) > 1e-3 or abs(q - r["charge"]) > 1e-3:
        why.append(f"charges sum to {q:+.4f}, RESI says {r['charge']:+g}")
    if p.returncode:
        why.append("caps: " + ((p.stderr or out).strip().splitlines() or ["?"])[-1][:100])
    if "missing parameters" in out:
        why.append("missing: " + "; ".join(re.findall(r"^  (.+)$", out.split("missing parameters", 1)[1], re.M)[:4]))
    if m := re.search(r"(\d+) interactions \(\d+ kinds\) have no entry", out):
        why.append(f"{m.group(1)} terms by analogy")
    got = {frozenset(map(int, l.split()[1:5])) for l in out.splitlines() if l.startswith("improper ")}
    tmap = {a.upper(): t for a, t, _ in r["atoms"]}
    want = {frozenset(idx[x.upper()] + 1 for x in quad) for quad in r["impr"]
            if all(x.upper() in idx for x in quad) and not zero_improper([tmap[x.upper()] for x in quad])}
    if not p.returncode and want - got:
        why.append(f"impropers: {len(want - got)} of the topology's missing")
    if not p.returncode and got - want:
        extra[r["name"]] = len(got - want)
    counts["checked"] += 1
    if why:
        failures.append((r["name"], why))
    else:
        counts["complete"] += 1

print(f"{os.path.basename(rtf)} with {os.path.basename(ffjson)} ({'CAPS typing ' + os.path.basename(auto) if auto else 'the topology types'}): "
      f"{counts['complete']} of {counts['checked']} residues complete (integer charge, every bonded term, every improper the topology lists)"
      + "".join(f" · {k} {v}" for k, v in counts.items() if k.startswith("skipped")))
if extra:
    print(f"  extra impropers (a parameter-file improper the residue's IMPR records do not list): {sum(extra.values())} in {len(extra)} residues"
          + (": " + ", ".join(f"{k} {v}" for k, v in extra.items()) if "--all" in sys.argv else ""))
for name, why in failures[: (None if "--all" in sys.argv else 30)]:
    print(f"  {name}: " + " | ".join(why))
sys.exit(1 if failures else 0)
