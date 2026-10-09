#!/usr/bin/env python3
"""Check CAPS's automatic atom typing against DL_FIELD's molecule templates.

Every template in a DL_FIELD library (.sf MOLECULE blocks) lists its atoms with DL_FIELD atom types and its bonds;
the library's ATOM_TYPE table maps each DL_FIELD type to the force-field key the parameters use. This script writes
each template as a mol2 file WITHOUT bond orders (CAPS must perceive them from the connectivity, as for a PDB or
xyz input), types it with `caps ff type`, and compares the result with the template's keys.

When DL_FIELD itself was run on the template (validate_family.py work directory, --dl DIR), its re-typing
("Atom N - X changed to Y" in dl_field.output) is applied to the expected keys, since that is what DL_FIELD assigns.

usage: validate_typing.py LIB FF.json RULES.json [--dl DIR] [--only a,b] [--show N] [--all] [--apply] [--charges types|gasteiger] [--exclude REGEX]
"""
import collections, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import template_pdb

LIB = os.environ.get("DLFIELD", os.path.expanduser("~/project/dl_f_4.13")) + "/lib"
CAPS = os.path.join(ROOT, "build", "cli", "caps")
lib, ffjson, rules = sys.argv[1], sys.argv[2], sys.argv[3]
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
dl_dir = arg("--dl", "")
show = int(arg("--show", "40"))
only = [x for x in arg("--only", "").split(",") if x]

sf = os.path.join(LIB, lib + ".sf")
types, keys, mols = template_pdb.parse_sf(sf)
key_of = dict(template_pdb.ATOM_KEY)
work = tempfile.mkdtemp()


def usable(name, m):
    """Whole molecules only: every CONNECT partner is in the template (no links to neighbouring residues)."""
    labels = {a for a, _ in m["atoms"]}
    if len(m["atoms"]) < 1:
        return False
    if any(b not in labels for js in m["connect"].values() for b in js):
        return False
    text = (name + " " + m["remark"]).lower()
    # chain residues and terminal residues are fragments with open valences, not molecules
    return "link" not in text and not re.search(r"_residue|_(nh3|coo|s)$", name.lower())


def element(t):
    e = types.get(t, "")
    e = re.sub(r"[^A-Za-z]", "", e)
    return e[:1].upper() + e[1:2].lower() if len(e) > 1 and e[:2].capitalize() in template_pdb.COV else e[:1].upper()


def write_mol2(path, name, m):
    labels = [a for a, _ in m["atoms"]]
    idx = {a: i for i, a in enumerate(labels)}
    bonds = set()
    for a, js in m["connect"].items():
        for b in js:
            if a in idx and b in idx and a != b:
                bonds.add((min(idx[a], idx[b]), max(idx[a], idx[b])))
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{name}\n{len(labels)} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, (lab, t) in enumerate(m["atoms"]):
            f.write(f"{i + 1} {lab} {i * 1.5:.3f} 0.000 0.000 {types.get(t, 'Du')} 1 MOL 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for k, (i, j) in enumerate(sorted(bonds)):
            f.write(f"{k + 1} {i + 1} {j + 1} un\n")
    return len(labels)


def retyped(name, dl_types):
    """DL_FIELD's own re-typing of this template, when it was run on it."""
    out = os.path.join(dl_dir, name, "dl_field.output") if dl_dir else ""
    t = list(dl_types)
    if out and os.path.exists(out):
        for line in open(out, errors="replace"):
            mm = re.search(r"Atom\s+(\d+)\s*-\s*(\S+)\s+changed to\s+(\S+)", line)
            if mm and int(mm.group(1)) <= len(t):
                t[int(mm.group(1)) - 1] = mm.group(3)
    return t


names = only or [n for n, m in mols.items() if n in keys and usable(n, m)]
if "--exclude" in sys.argv:   # e.g. water-model templates (TIP4P / TIP5P virtual sites cannot be typed from bonds)
    names = [n for n in names if not re.search(arg("--exclude", ""), n)]
total = right = 0
wrong_mols, confusion, failed, apply_problems = [], collections.Counter(), [], []
for name in names:
    m = mols[name]
    if any(t not in types for _, t in m["atoms"]):
        failed.append((name, "atom type missing from ATOM_TYPE"))
        continue
    path = os.path.join(work, re.sub(r"[^A-Za-z0-9_.-]", "_", name) + ".mol2")
    write_mol2(path, name, m)
    expect = [key_of.get(t, "?") for t in retyped(name, [t for _, t in m["atoms"]])]
    out = path + ".types"
    r = subprocess.run([CAPS, "ff", "type", path, "--ff", ffjson, "--typing", rules, "-o", out, "--explain"],
                       capture_output=True, text=True)
    if not os.path.exists(out):
        failed.append((name, (r.stderr or r.stdout).strip().splitlines()[-1][:120] if (r.stderr or r.stdout).strip() else "?"))
        continue
    got = [x.strip() for x in open(out)]
    if "--apply" in sys.argv:
        # end to end: the automatically typed molecule must get every parameter
        ra = subprocess.run([CAPS, "ff", "apply", path, "--ff", ffjson, "--typing", rules, "--charges", arg("--charges", "types"), "--allow-missing"],
                            capture_output=True, text=True)
        miss = re.findall(r"^  (.+)$", ra.stdout.split("missing parameters", 1)[1], re.M) if "missing parameters" in ra.stdout else []
        if ra.returncode or miss:
            apply_problems.append((name, ra.returncode, miss[:6], (ra.stderr.strip().splitlines() or [""])[-1][:100]))
    bad = [(i, a[0], a[1], expect[i], got[i]) for i, a in enumerate(m["atoms"]) if expect[i] != got[i]]
    total += len(expect)
    right += len(expect) - len(bad)
    for b in bad:
        confusion[(b[3], b[4])] += 1
    if bad:
        wrong_mols.append((name, bad, r.stdout))

print(f"{lib}: {len(names)} templates, {right}/{total} atoms typed as DL_FIELD's templates "
      f"({100.0 * right / max(total, 1):.1f} %), {len(names) - len(wrong_mols) - len(failed)} molecules fully right, "
      f"{len(wrong_mols)} with differences, {len(failed)} not run")
if confusion:
    print("expected → CAPS (atoms):")
    for (e, g), k in confusion.most_common(show):
        print(f"  {e:>6} → {g:<6} {k}")
for name, why in failed[:min(show, 5)]:
    print(f"  not run: {name}: {why}")
if "--apply" in sys.argv:
    print(f"end to end (ff apply with automatic types): {len(names) - len(failed) - len(apply_problems)} molecules fully parameterised, "
          f"{len(apply_problems)} with missing terms or errors")
    for name, rc, miss, err in apply_problems[:show]:
        print(f"  {name}: " + ("missing " + "; ".join(miss) if miss else err))
if "--all" in sys.argv:
    for name, bad, log in wrong_mols:
        print(f"--- {name}")
        for i, lab, dlt, e, g in bad:
            print(f"   {i + 1:3d} {lab:6s} {dlt:24s} expected {e:6s} got {g}")
