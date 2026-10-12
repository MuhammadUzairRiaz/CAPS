#!/usr/bin/env python3
"""Check CAPS's automatic atom typing against DL_FIELD's molecule templates.

Every template in a DL_FIELD library (.sf MOLECULE blocks) lists its atoms with DL_FIELD atom types and its bonds;
the library's ATOM_TYPE table maps each DL_FIELD type to the force-field key the parameters use. This script writes
each template as a mol2 file WITHOUT bond orders (CAPS must perceive them from the connectivity, as for a PDB or
xyz input), types it with `caps ff type`, and compares the result with the template's keys.

When DL_FIELD itself was run on the template (validate_family.py work directory, --dl DIR), its re-typing
("Atom N - X changed to Y" in dl_field.output) is applied to the expected keys, since that is what DL_FIELD assigns.

--cap-open: residue templates whose only fault is an open valence (an sp3 C with three partners, an ether or hydroxyl O
with one, by the DL type's name: a sugar's anomeric C or glycosidic O) are typed with H on each open valence; the atoms
the cap retypes (the capped O, the H already on a capped C) are left out of the comparison.

usage: validate_typing.py LIB FF.json RULES.json [--dl DIR] [--only a,b] [--show N] [--all] [--apply] [--charges types|gasteiger] [--exclude REGEX] [--cap-open]
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
work = arg("--keep", "") or tempfile.mkdtemp()   # --keep DIR: the mol2 and type files stay there
os.makedirs(work, exist_ok=True)


UNITED = False   # united-atom rules: a CH2 site has two partners, nothing is open


def open_atoms(m):
    """{atom label: missing partners} for sp3 C and ether / hydroxyl O the template leaves open"""
    if UNITED:
        return {}
    nbs = collections.defaultdict(set)
    for a, js in m["connect"].items():
        for b in js:
            if a != b:
                nbs[a].add(b); nbs[b].add(a)
    # under-saturated atoms that a double bond can satisfy (furan's C=C, a carbonyl's C=O) are not open: a maximum
    # matching on the bonds between them; whatever stays unmatched is open
    full = {"C": 4, "N": 3, "O": 2, "S": 2}
    el = {a: element(t) for a, t in m["atoms"]}
    low = {a for a in el if el[a] in full and len(nbs[a]) < full[el[a]]}
    mate = {}

    def augment(a, seen):
        for b in nbs[a]:
            if b in low and b not in seen:
                seen.add(b)
                if b not in mate or augment(mate[b], seen):
                    mate[a], mate[b] = b, a
                    return True
        return False

    for a in sorted(low, key=lambda a: ({"O": 0, "S": 1, "N": 2}.get(el[a], 3), a)):   # a lone O / S / N takes its double bond first
        if a not in mate:
            augment(a, {a})
    out = {}
    for a, t in m["atoms"]:
        if a in mate:
            continue
        tl = t.lower()
        if element(t) == "C" and re.search(r"aliphatic|alkane|sp3|sialic|sugar|furan|pyran", tl) and len(nbs[a]) < 4:
            out[a] = 4 - len(nbs[a])
        elif element(t) == "O" and re.search(r"ether|hydroxyl|alcohol|glycosid|sugar", tl) and len(nbs[a]) < 2:
            out[a] = 2 - len(nbs[a])
        elif element(t) == "O" and re.match(r"ol_", tl) and len(nbs[a]) < 2:
            out[a] = 2 - len(nbs[a])   # a nucleotide's O3' linking to the next phosphate
        elif element(t) == "P" and re.search(r"phosph|po4", tl) and len(nbs[a]) < 4:
            out[a] = 4 - len(nbs[a])   # a nucleotide's phosphate linking to the previous sugar: capped as a methyl ester
        elif element(t) == "S" and re.search(r"disulph|disulf", tl) and len(nbs[a]) < 2:
            out[a] = 2 - len(nbs[a])   # half a disulfide bridge
        elif element(t) == "S" and len(nbs[a]) == 3 and all(element(dict(m["atoms"]).get(b, "")) == "O" for b in nbs[a]):
            out[a] = 1   # an SO3 group that a sugar O carries (sulfate ester)
    return out


def has_ring_o(m):
    """does the template hold a ring with an O (a sugar ring)?"""
    nbs = collections.defaultdict(set)
    for a, js in m["connect"].items():
        for b in js:
            if a != b:
                nbs[a].add(b); nbs[b].add(a)
    el = dict(m["atoms"])
    for o in (a for a, t in m["atoms"] if element(t) == "O" and len(nbs[a]) == 2):
        x, y = list(nbs[o])
        seen, front = {x}, [x]   # is y reachable from x without passing o?
        while front:
            front = [c for f in front for c in nbs[f] if c != o and c not in seen and not seen.add(c)]
        if y in seen:
            return True
    return False


def usable(name, m, allow_open=False):
    """Whole molecules only: every CONNECT partner is in the template (no links to neighbouring residues)."""
    labels = {a for a, _ in m["atoms"]}
    if len(m["atoms"]) < 1:
        return False
    if any(b not in labels for js in m["connect"].values() for b in js):
        return False
    if m.get("links"):   # an improper through a neighbouring residue's atom: a chain residue with an open valence
        return False
    # an open valence the template leaves for a neighbouring residue (a sugar's glycosidic C and O): an sp3 carbon with
    # fewer than four partners, an ether or hydroxyl O with fewer than two, by the DL type's own name
    if open_atoms(m) and not allow_open:
        return False
    text = (name + " " + m["remark"]).lower()
    # chain residues and terminal residues are fragments with open valences, not molecules
    return "link" not in text and not re.search(r"_residue|_(nh3|coo|s)$", name.lower())


def element(t):
    e = types.get(t, "")
    e = re.sub(r"[^A-Za-z]", "", e)
    return e[:1].upper() + e[1:2].lower() if len(e) > 1 and e[:2].capitalize() in template_pdb.COV else e[:1].upper()


# a sugar ring to cap the open valence of a group that only exists on a sugar (an aglycone, an O-methyl, a sulfate):
# oxan-2-yl, its C1' carrying the open atom (an O) or an O carrying it (a C or S)
OXANYL = [("C1", "C", ["O5", "C2", "H1"]), ("O5", "O", ["C5"]), ("C5", "C", ["C4", "H51", "H52"]), ("C4", "C", ["C3", "H41", "H42"]),
          ("C3", "C", ["C2", "H31", "H32"]), ("C2", "C", ["H21", "O2"]), ("O2", "O", ["HO2"])] + \
         [(h, "H", []) for h in ("H1", "H51", "H52", "H41", "H42", "H31", "H32", "H21", "HO2")]


def write_mol2(path, name, m, caps=None, sugar=False):
    """caps: {label: n} adds n H on that atom (open valences); sugar: one open valence per atom gets an oxan-2-yl ring
    (through an O when the open atom is not itself an O) instead"""
    m = {"atoms": list(m["atoms"]), "connect": {a: list(js) for a, js in m["connect"].items()}}
    for a, n in (caps or {}).items():
        if element(dict(m["atoms"])[a]) == "P":   # O-CH3
            pre = f"CAP{a}_"
            m["atoms"] += [(pre + "O", "__cap_O"), (pre + "C", "__cap_C")] + [(pre + f"H{k}", "__cap_H") for k in range(3)]
            m["connect"][pre + "O"] = [pre + "C"]
            m["connect"][pre + "C"] = [pre + f"H{k}" for k in range(3)]
            m["connect"].setdefault(a, []).append(pre + "O")
            n -= 1
        elif sugar:
            pre = f"CAP{a}_"
            for lab, el, js in OXANYL:
                m["atoms"].append((pre + lab, "__cap_" + el))
                m["connect"][pre + lab] = [pre + j for j in js]
            if element(dict(m["atoms"])[a]) == "O":
                m["connect"].setdefault(a, []).append(pre + "C1")
            else:
                m["atoms"].append((pre + "OL", "__cap_O"))
                m["connect"][pre + "OL"] = [pre + "C1"]
                m["connect"].setdefault(a, []).append(pre + "OL")
            n -= 1
        for k in range(n):
            h = f"HCAP{a}_{k}"
            m["atoms"].append((h, "__cap_H"))
            m["connect"].setdefault(a, []).append(h)
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
            el = t[6:] if t.startswith('__cap_') else types.get(t, 'Du')
            if UNITED and not t.startswith("__cap_") and element(t) == "C" and t in template_pdb.ATOM_MASS:
                # a united-atom site (the library's element column may hold the site name, CE1): carbon, named by its
                # hydrogens (CH2) as CAPS reads united-atom sites
                el = "C"
                nh = round((template_pdb.ATOM_MASS[t] - 12.011) / 1.008)
                if 1 <= nh <= 4:
                    el = f"CH{nh}"   # the mol2 type column, which CAPS reads as the site name
            f.write(f"{i + 1} {lab} {i * 1.5:.3f} 0.000 0.000 {el} 1 MOL 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for k, (i, j) in enumerate(sorted(bonds)):
            f.write(f"{k + 1} {i + 1} {j + 1} un\n")
    return len(labels)


def symmetry_classes(m):
    """atoms the molecule's graph cannot tell apart (pyridine's C3 and C5): colour refinement on element and neighbours"""
    labels = [a for a, _ in m["atoms"]]
    nbs = collections.defaultdict(set)
    for a, js in m["connect"].items():
        for b in js:
            if a != b:
                nbs[a].add(b); nbs[b].add(a)
    col = {a: element(t) for a, t in m["atoms"]}
    for _ in range(len(labels)):
        new = {a: (col[a], tuple(sorted(col[b] for b in nbs[a] if b in col))) for a in labels}
        ids = {v: k for k, v in enumerate(sorted(set(new.values()), key=repr))}
        new = {a: ids[new[a]] for a in labels}
        if len(set(new.values())) == len(set(col.values())):
            col = new
            break
        col = new
    return [col[a] for a in labels]


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


UNITED = bool(__import__("json").load(open(rules)).get("united_atom"))
CAP = "--cap-open" in sys.argv
# --deviations FILE: atoms where the library's template contradicts the force field's own typer (antechamber for GAFF),
# each with its evidence: reported as explained, not counted right
DEV = {}
if "--deviations" in sys.argv:
    for d in __import__("json").load(open(arg("--deviations", "")))["deviations"]:
        if d["lib"] == lib and d.get("ff", os.path.basename(ffjson)[:-5]) == os.path.basename(ffjson)[:-5]:
            for a in d["atoms"]:
                DEV[(d["template"], a)] = (d["dl"], d["caps"], d["reason"])
# a CAPS type whose every parameter class is another type's (equivalence vdw / bond / angle / dihedral all X, as
# misc's CTf3 / CTf1 for CT, which only carry a molecule's own charges) is that type under another name
ALIAS = {}
for t in __import__("json").load(open(ffjson))["atom_types"]:
    eq = t.get("equivalence") or {}
    if eq and len(set(eq.values())) == 1 and set(eq) >= {"vdw", "bond"}:
        ALIAS[t["name"]] = next(iter(eq.values()))
# --overlay 'REGEX=RULES.json;...': templates whose name matches take that typing overlay too (a water model: the same
# structure, a different model, which the user picks)
OVERLAYS = [tuple(x.split("=", 1)) for x in arg("--overlay", "").split(";") if x]
names = only or [n for n, m in mols.items() if n in keys and usable(n, m, CAP)]
if "--exclude" in sys.argv:   # e.g. water-model templates (TIP4P / TIP5P virtual sites cannot be typed from bonds)
    names = [n for n in names if not re.search(arg("--exclude", ""), n)]
total = right = sym_total = alias_total = 0
explained = []
wrong_mols, confusion, failed, apply_problems = [], collections.Counter(), [], []
for name in names:
    m = mols[name]
    if any(t not in types for _, t in m["atoms"]):
        failed.append((name, "atom type missing from ATOM_TYPE"))
        continue
    path = os.path.join(work, re.sub(r"[^A-Za-z0-9_.-]", "_", name) + ".mol2")
    opened = open_atoms(m) if CAP else {}
    sugar_cap = bool(opened) and not has_ring_o(m)   # a group that only exists on a sugar: capped with one
    write_mol2(path, name, m, opened, sugar_cap)
    over = [f for pat, f in OVERLAYS if re.fullmatch(pat, name)]
    expect = [key_of.get(t, "?") for t in retyped(name, [t for _, t in m["atoms"]])]
    # atoms a cap retypes: a capped O itself, and the H already on a capped C
    skip = {a for a in opened if element(dict(m["atoms"])[a]) == "O" and not sugar_cap}
    for a in ([] if sugar_cap else opened):
        for b in m["connect"].get(a, []):
            if element(dict(m["atoms"]).get(b, "")) == "H":
                skip.add(b)
        skip |= {x for x, js in m["connect"].items() if a in js and element(dict(m["atoms"]).get(x, "")) == "H"}
    out = path + ".types"
    r = subprocess.run([CAPS, "ff", "type", path, "--ff", ffjson, "--typing", ",".join([rules] + over), "-o", out, "--explain"],
                       capture_output=True, text=True)
    if not os.path.exists(out):
        failed.append((name, (r.stderr or r.stdout).strip().splitlines()[-1][:120] if (r.stderr or r.stdout).strip() else "?"))
        continue
    got = [x.strip() for x in open(out)][:len(m["atoms"])]
    if "--apply" in sys.argv:
        # end to end: the automatically typed molecule must get every parameter
        ra = subprocess.run([CAPS, "ff", "apply", path, "--ff", ffjson, "--typing", rules, "--charges", arg("--charges", "types"), "--allow-missing"],
                            capture_output=True, text=True)
        miss = re.findall(r"^  (.+)$", ra.stdout.split("missing parameters", 1)[1], re.M) if "missing parameters" in ra.stdout else []
        if ra.returncode or miss:
            apply_problems.append((name, ra.returncode, miss[:6], (ra.stderr.strip().splitlines() or [""])[-1][:100]))
    # a key that the library gives a symmetry-equivalent atom of the same molecule (pyridine's C5_2 for C3) is the same
    # type under another name: no typer can tell those atoms apart
    sym = symmetry_classes(m)
    same_class = collections.defaultdict(set)
    for i, c in enumerate(sym):
        same_class[c].add(expect[i])
    by_alias = [i for i, a in enumerate(m["atoms"]) if expect[i] != got[i] and a[0] not in skip and ALIAS.get(got[i]) == expect[i]]
    alias_total += len(by_alias)
    by_symmetry = [i for i, a in enumerate(m["atoms"]) if expect[i] != got[i] and a[0] not in skip and i not in by_alias
                   and got[i] in same_class[sym[i]]]
    by_symmetry += by_alias
    sym_total += len(by_symmetry)
    bad = [(i, a[0], a[1], expect[i], got[i]) for i, a in enumerate(m["atoms"])
           if expect[i] != got[i] and a[0] not in skip and i not in by_symmetry]
    dev = [b for b in bad if DEV.get((name, b[1]), (None, None))[:2] == (b[3], b[4])]
    for b in dev:
        explained.append((name, b[1], b[3], b[4], DEV[(name, b[1])][2]))
    bad = [b for b in bad if b not in dev]
    total += len(expect) - len(skip) - len(dev)
    right += len(expect) - len(skip) - len(bad) - len(dev)
    for b in bad:
        confusion[(b[3], b[4])] += 1
    if bad:
        wrong_mols.append((name, bad, r.stdout))

print(f"{lib}: {len(names)} templates, {right}/{total} atoms typed as DL_FIELD's templates "
      f"({100.0 * right / max(total, 1):.1f} %), {len(names) - len(wrong_mols) - len(failed)} molecules fully right, "
      f"{len(wrong_mols)} with differences, {len(failed)} not run" +
      (f"; {sym_total - alias_total} atoms carry the key of a symmetry-equivalent atom (the library numbers equivalent ring positions apart)" if sym_total > alias_total else "") +
      (f"; {alias_total} atoms carry a CAPS name for the key (every parameter class equivalent to it)" if alias_total else ""))
if explained:
    print(f"{len(explained)} atoms where the library's template contradicts the force field's own typer (not counted above):")
    for e in sorted({(x[0], x[2], x[3], x[4]) for x in explained}):
        print(f"  {e[0]}: library {e[1]}, CAPS {e[2]}: {e[3]}")
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
