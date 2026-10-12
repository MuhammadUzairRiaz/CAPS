#!/usr/bin/env python3
"""PCFF, COMPASS and CVFF (from their .frc files) on the polymer library: one table per force field, one row per polymer.

For every oligomer (a directory of .xyz files, one per polymer of data/polymers/library.json) and each of pcff-frc,
compass-frc and cvff-frc:

  · typing: atoms CAPS types and the atoms no rule matches; for COMPASS, why (the chemistry compass_published.frc has
    no type for, named by the atom's PCFF type);
  · missing terms: with the .frc file's explicit tables only and with its automatic parameters (caps ff apply
    --allow-missing lists them);
  · against DL_FIELD (compare_polymers.py --keep work directories): types, charges, every bond, angle, torsion and
    out-of-plane term compared by its atoms, and the LAMMPS energies per term. DL_FIELD writes no class II cross terms
    and no out-of-plane terms for PCFF / COMPASS: CAPS's file is run again with those terms set to zero, so the terms
    both programs write are compared on equal footing;
  · against msi2lmp (check_msi2lmp.py --polymers --json): the LAMMPS energy per term on msi2lmp's data file and on
    CAPS's for the same typed structure.

usage: class2_report.py --olig DIR --dlf FF=DIR[,FF=DIR] --msi MSI.json --out DIR [--only P0]
Writes OUT/{pcff,compass,cvff}_polymers.csv, OUT/class2_polymers.json and OUT/class2_summary.csv; prints the totals.
Every difference from DL_FIELD's parameters is matched against EXPLAINED (its cause, checked against the .frc file);
"dl_term_groups_unexplained" counts those that are not.
"""
import collections, csv, json, math, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import compare_polymers   # noqa: E402  (its readers of DL_FIELD's and CAPS's LAMMPS files)
CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FFS = {"pcff-frc": "pcff", "compass-frc": "compass_published", "cvff-frc": "cvff"}   # CAPS library → msi2lmp's .frc
CROSS = ("BondBond", "BondAngle", "MiddleBondTorsion", "EndBondTorsion", "AngleTorsion", "AngleAngleTorsion", "BondBond13", "AngleAngle")
TOL = 2e-4   # relative, for coefficients DL_FIELD writes with six decimals
DL_TOL = 1e-4   # relative, for DL_FIELD's LAMMPS energies (its coefficients and mixed pair terms carry six digits)

# why COMPASS (compass_published.frc) leaves an atom untyped, by the atom's PCFF type (pcff.frc's description)
COMPASS_GAP = [
    (r"^c=|^c=1|^c=2", "non-aromatic C=C carbon (no alkene type)"),
    (r"^(cl|f|br|i)$", "halogen (no halogen types)"),
    (r"^c_0$", "aldehyde or ketone C (no aldehyde or ketone type: the aldehyde is the H-capped end of an ester or amide oligomer)"),
    (r"^n$|^hn$|^nn$|^n2$|^hn2$|^n3n$|^nh$|^n_1|^n_2|^na$", "amine or N-H amide N, or its H (only n3m, an amide N without H)"),
    (r"^ct$|^nt$|^nz$", "nitrile C or N (no nitrile types)"),
    (r"^sf$|^s'$|^s$|^sp$|^s1$|^sc$", "sulfur (no sulfide or sulfone types)"),
    (r"^np$|^nb$|^n=|^npc$|^nr$", "aromatic or sp2 N (no pyridine, pyrrole or aniline N)"),
    (r"^cp$|^c5$|^cs$", "aromatic C bonded to a heteroatom or in a five-membered ring"),
]


# every difference from DL_FIELD's parameters found on the polymer library, with its cause (checked against the .frc file;
# for explicit-table terms msi2lmp assigns what CAPS assigns). Matched against "kind types: ..." of compare_terms.
EXPLAINED = {
    "pcff-frc": [
        (r"^torsion cp-cp-cp-cp:", "DL_FIELD drops the V1 term (8.3667) of pcff.frc's cp cp cp cp torsion"),
        (r"^torsion (\S*o_1-c_1-o_2-\S*|\S*-o_2-c_1-o_1|\S*-o_2-c_1-cp|cp-c_1-o_2-\S*|\S*-c_1-o_2-ho2):",
         "ester torsions: DL_FIELD's own values; CAPS pcff.frc's torsion_3 entries (as msi2lmp)"),
        (r"^torsion \S*(-n_2|-n|n_2-|n-)\S*: DL_FIELD \[\(0.158, 3, 0\)\]",
         "C-C-C-N torsion: DL_FIELD the generic cff91_auto * c_ c_ * (0.158) where pcff.frc lists c c c n / n_2 in torsion_3"),
        (r"^torsion \S*(cl|f)\S*: DL_FIELD \[\(0.158, 3, 0\)\]",
         "C-C-C-halogen torsion: DL_FIELD the generic cff91_auto * c_ c_ * (0.158) where pcff.frc lists c c c cl / f (V3 0.15)"),
        (r"^(bond|angle) \S*(c3m|c3h|o3e)\S*:",
         "three-membered ring: DL_FIELD cff91_auto values (ring angles 60 deg); pcff.frc's equivalences give c3m / c3h the class II "
         "terms of c and o3e those of o (as msi2lmp)"),
        (r"^bond sio-hsi:", "DL_FIELD the cff91_auto quadratic bond; pcff.frc's quartic sio-hsi (by equivalence, as msi2lmp)"),
        (r"^angle \S*ct\S*:", "nitrile C next to sp3 C: pcff.frc has no quartic c c ct angle (msi2lmp: none); DL_FIELD takes c c c's; "
                               "CAPS the cff91_auto entries (*7 c_ c_ 46.0, *2 c_ h_ 44.0)"),
        (r"^torsion cl-c=-c=-\S*:", "polychloroprene: follows from DL_FIELD's c=2 typing of the C=C carbon (listed typing deviation)"),
    ],
    "cvff-frc": [
        (r"^out-of-plane at \S+: DL_FIELD \[.*\] CAPS \[\]",
         "out-of-plane term at a centre cvff.frc's #out_of_plane does not list (ester / ketone c', substituted cp, trisubstituted c=, "
         "amide or aniline n): DL_FIELD adds one, msi2lmp and CAPS none"),
        (r"^torsion \S*-c'-n-\S*:", "amide torsions: DL_FIELD the cvff_auto * c'_ n_ * barrier spread over the bond (0.8) where cvff.frc "
                                    "lists c c' n c / hn and o' c' n c / hn (as msi2lmp)"),
        (r"^torsion \S+: DL_FIELD \[\((0.1778|0.195), 3, 0\)\]",
         "wildcard torsion (* c c *, * c o *) in a three-membered ring: DL_FIELD divides the barrier by the torsions present (8, 2), "
         "msi2lmp and CAPS by (n_j - 1)(n_k - 1) (9, 3)"),
        (r"^torsion (cp-cp-o-\S*|\S*-o-cp-cp):", "aryl ether: DL_FIELD the cvff_auto * cp_ o_ * (1.8) where cvff.frc lists * cp o * (1.5)"),
        (r"^torsion (hn-n2-cp-cp|cp-cp-n2-hn):", "aniline N: DL_FIELD the cvff_auto * cp_ n_ * (2.25) where cvff.frc lists * cp n2 * (10.0 over 4)"),
        (r"^torsion cp-cp-c-cp:", "cvff.frc's cp cp c cp (ver 1.9, 0.675, n 4), missing from DL_FIELD's CVFF"),
        (r"^Lennard-Jones c: DL_FIELD \(0.039, 3.875\) CAPS \(0.16, 3.475\)",
         "quaternary sp3 C (type c): DL_FIELD gives it cg's Lennard-Jones; CAPS cvff.frc's own c entry (listed parameter deviation)"),
    ],
    "compass-frc": [
        (r"^torsion \S*(o2e|c4o)\S*:", "C-C-O ether and alcohol torsions: DL_FIELD's COMPASS has entries compass_published.frc does not; CAPS the published "
                                  "file's c4 c4 entry (as msi2lmp)"),
        (r"^bond o2z-h1o:", "silanol end of PDMS: no o2z-h1o bond in compass_published.frc (missing term)"),
    ],
}


def explain(ff, key):
    return next((why for pat, why in EXPLAINED.get(ff, []) if re.search(pat, key)), None)


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def read_data(path):
    """a LAMMPS data file: labels of the atom types, atoms (type, charge), topology and coefficients per section"""
    d = {"label": {}, "atoms": {}, "topo": collections.defaultdict(list), "coef": collections.defaultdict(dict), "style": {}}
    sec = None
    for raw in open(path):
        body = raw.split("#")[0].strip()
        if not body:
            continue
        if re.match(r"^[A-Z][A-Za-z0-9 ]+$", body):
            sec = body.replace(" Coeffs", "_C").strip()
            if "#" in raw and sec.endswith("_C"):
                d["style"][sec[:-2]] = raw.split("#")[1].split()[0]
            continue
        w = body.split()
        if sec == "Masses":
            d["label"][int(w[0])] = raw.split("#")[1].split()[0] if "#" in raw else w[0]
        elif sec == "Atoms" and len(w) >= 7:
            d["atoms"][int(w[0])] = (int(w[2]), float(w[3]))
        elif sec in ("Bonds", "Angles", "Dihedrals", "Impropers"):
            d["topo"][sec].append((int(w[1]), tuple(int(x) for x in w[2:])))
        elif sec and sec.endswith("_C"):
            style = None
            vals = w[1:]
            if vals and re.match(r"^[a-z]", vals[0]):
                style, vals = vals[0], vals[1:]
            d["coef"][sec[:-2]][int(w[0])] = (style, [float(x) for x in vals if not re.match(r"^[a-z]", x)])
    return d


def lab(d, a):
    return d["label"][d["atoms"][a][0]]


def terms(d, cls):
    """every bonded term keyed by its atoms (orientation-free), with its parameters in one normal form"""
    out = {"bond": {}, "angle": {}, "torsion": collections.defaultdict(list), "oop": collections.defaultdict(list)}
    nb = collections.defaultdict(set)
    for t, (i, j) in d["topo"]["Bonds"]:
        nb[i].add(j)
        nb[j].add(i)
    for t, (i, j) in d["topo"]["Bonds"]:
        st, v = d["coef"]["Bond"][t]
        st = st or d["style"].get("Bond", "")
        if st == "zero" or not v:
            continue
        out["bond"][frozenset((i, j))] = (v[1], v[0]) if cls == "class1" else tuple(v[:4])   # (r0, K) or (r0, K2, K3, K4)
    for t, (i, j, k) in d["topo"]["Angles"]:
        st, v = d["coef"]["Angle"][t]
        out["angle"][(min(i, k), j, max(i, k))] = (v[1], v[0]) if cls == "class1" else tuple(v[:4])   # (θ0, K) or (θ0, K2, K3, K4)
    for t, a in d["topo"]["Dihedrals"]:
        st, v = d["coef"]["Dihedral"][t]
        st = st or d["style"].get("Dihedral", "")
        key = a if a[1] < a[2] else a[::-1]
        if st == "harmonic":                 # K [1 + d cos nφ]
            parts = [(v[0], int(v[2]), 0.0 if v[1] > 0 else 180.0)]
        elif st == "charmm":                 # K [1 + cos(nφ - d)]
            parts = [(v[0], int(v[1]), v[2] % 360)]
        elif st == "fourier":                # m (K n d)…
            parts = [(v[1 + 3 * q], int(v[2 + 3 * q]), v[3 + 3 * q] % 360) for q in range(int(v[0]))]
        elif st == "class2":                 # Σ K_n [1 - cos(nφ - φ_n)] = K_n [1 + cos(nφ - φ_n - 180)]
            parts = [(v[2 * q], q + 1, (v[2 * q + 1] + 180) % 360) for q in range(3)]
        else:
            raise ValueError(f"dihedral style {st}")
        out["torsion"][key] += [p for p in parts if abs(p[0]) > 1e-9]
    for t, a in d["topo"]["Impropers"]:
        st, v = d["coef"]["Improper"][t]
        centre = next((x for x in a if all(y in nb[x] for y in a if y != x)), a[1])
        if abs(v[0]) > 1e-9:   # a zero term (DL_FIELD writes some) is no term
            out["oop"][centre].append(tuple(round(x, 4) for x in v[:3]) if cls == "class1" else (round(v[0], 4), round(v[1], 3)))
    return out


def close(x, y):
    return len(x) == len(y) and all(abs(a - b) <= TOL * max(1.0, abs(a)) for a, b in zip(x, y))


def compare_terms(dl, cp, cls, oop=True):
    """differences between DL_FIELD's and CAPS's terms, grouped by the atom types of the term"""
    A, B = terms(dl, cls), terms(cp, cls)
    diff = collections.Counter()
    names = lambda d, key: "-".join(lab(d, a) for a in key)
    for kind in ("bond", "angle"):
        for key in set(A[kind]) | set(B[kind]):
            x, y = A[kind].get(key), B[kind].get(key)
            if x is None or y is None or not close(x, y):
                k = sorted(key) if kind == "bond" else key
                diff[f"{kind} {names(cp, k)}: DL_FIELD {fmt(x)} CAPS {fmt(y)}"] += 1
    for key in set(A["torsion"]) | set(B["torsion"]):
        x, y = sorted(A["torsion"].get(key, [])), sorted(B["torsion"].get(key, []))
        if len(x) != len(y) or not all(close(p, q) for p, q in zip(x, y)):
            diff[f"torsion {names(cp, key)}: DL_FIELD {fmt(x)} CAPS {fmt(y)}"] += 1
    if oop:
        for c in set(A["oop"]) | set(B["oop"]):
            x, y = sorted(A["oop"].get(c, [])), sorted(B["oop"].get(c, []))
            if x != y:
                diff[f"out-of-plane at {lab(cp, c)}: DL_FIELD {fmt(x)} CAPS {fmt(y)}"] += 1
    return diff


def fmt(x):
    if x is None:
        return "none"
    if isinstance(x, list):
        return "[" + ", ".join(fmt(p) for p in x) + "]"
    return "(" + ", ".join(f"{v:.4g}" for v in x) + ")"


def lammps_energy(dirn, script, data_from=None, data_to=None):
    txt = open(os.path.join(dirn, script)).read()
    if data_from:
        txt = txt.replace(data_from, data_to)
    txt = re.sub(r"(?m)^(thermo_style|thermo_modify|run) .*$", "", txt)
    txt += ("\nthermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong\nrun 0\n"
            "print \"$(ebond:%.10f) $(eangle:%.10f) $(edihed:%.10f) $(eimp:%.10f) $(evdwl:%.10f) $(ecoul:%.10f) $(elong:%.10f)\" file e.txt\n")
    open(os.path.join(dirn, "report.in"), "w").write(txt)
    if os.path.exists(os.path.join(dirn, "e.txt")):
        os.remove(os.path.join(dirn, "e.txt"))
    run([LMP, "-in", "report.in", "-log", "none", "-screen", "none"], cwd=dirn, timeout=600)
    if not os.path.exists(os.path.join(dirn, "e.txt")):
        return None
    v = [float(x) for x in open(os.path.join(dirn, "e.txt")).read().split()]
    return dict(zip(("bond", "angle", "dihedral", "improper", "vdw", "coul", "long"), v))


def without_cross_terms(src, dst):
    """CAPS's class II data file with every cross term and the out-of-plane (Wilson) term set to zero"""
    out, sec = [], ""
    for line in open(src).read().splitlines():
        body = line.split("#")[0].strip()
        if body and re.match(r"^[A-Z][A-Za-z0-9 ]+$", body):
            sec = body
        elif body and (sec.replace(" Coeffs", "") in CROSS or sec == "Improper Coeffs"):
            w = body.split()
            line = " ".join([w[0]] + ["0.0" if sec != "Improper Coeffs" or k == 0 else x for k, x in enumerate(w[1:])])
        out.append(line)
    open(dst, "w").write("\n".join(out) + "\n")


def caps_typing(xyz, ff, work):
    """CAPS's type of every atom ('' where no rule matches)"""
    out = os.path.join(work, "types.txt")
    if os.path.exists(out):
        os.remove(out)
    run([CAPS, "ff", "type", xyz, "--ff", os.path.join(ROOT, "data", "forcefields", ff + ".json"), "-o", out])
    n = sum(1 for l in open(xyz).read().splitlines()[2:] if l.split())
    types = [("" if l.strip() == "?" else l.strip()) for l in open(out).read().split("\n")][:n] if os.path.exists(out) else [""] * n
    return types + [""] * (n - len(types))


def missing_terms(xyz, ff, explicit, work):
    """the interactions caps ff apply finds no parameters for, with the .frc's automatic parameters (explicit=False) or
    with its explicit tables only"""
    lib = os.path.join(ROOT, "data", "forcefields", ff + ".json")
    if explicit:
        j = json.load(open(lib))
        for k in [k for k in j if k.startswith("auto_")]:
            j.pop(k)
        if "typing" in j and not os.path.isabs(j["typing"]):
            j["typing"] = os.path.normpath(os.path.join(ROOT, "data", "forcefields", j["typing"]))
        lib = os.path.join(work, ff + ".explicit.json")
        json.dump(j, open(lib, "w"))
    r = run([CAPS, "ff", "apply", xyz, "--ff", lib, "--charges", "auto", "-o", os.path.join(work, "m.data")])
    text = r.stdout + r.stderr
    m = re.search(r"no parameters for (\d+) interaction\(s\):\n((?:  .*\n?)*)", text)
    if m:
        return [l.strip() for l in m.group(2).splitlines() if l.strip()], ""
    if r.returncode:
        return None, (text.strip().splitlines() or ["?"])[-1]
    return [], ""


def main():
    olig = arg("--olig", "")
    outdir = arg("--out", "")
    dlf = dict(x.split("=") for x in arg("--dlf", "").split(",") if x)
    msi = {}
    for path in arg("--msi", "").split(","):
        if path:
            for r in json.load(open(path)):
                msi[tuple(r["case"].split())] = r
    lib = json.load(open(os.path.join(ROOT, "data", "polymers", "library.json")))["polymers"]
    only = arg("--only", "")
    work = os.path.join(outdir, ".work")
    os.makedirs(work, exist_ok=True)
    rows = []
    for p in lib:
        xyz = os.path.join(olig, p["id"] + ".xyz")
        if not os.path.exists(xyz) or (only and only not in p["id"]):
            continue
        # the PCFF types name the chemistry COMPASS lacks
        pcff_types = caps_typing(xyz, "pcff-frc", work)
        for ff, frc in FFS.items():
            cls = "class1" if ff == "cvff-frc" else "class2"
            row = {"id": p["id"], "polymer": p["name"], "forcefield": ff}
            types = caps_typing(xyz, ff, work)
            row["atoms"] = len(types)
            row["typed"] = sum(1 for t in types if t)
            row["untyped"] = row["atoms"] - row["typed"]
            if row["untyped"]:
                why = collections.Counter()
                for t, pt in zip(types, pcff_types):
                    if not t:
                        why[next((r for pat, r in COMPASS_GAP if re.search(pat, pt)), f"PCFF type {pt}")] += 1
                row["untyped_reason"] = "; ".join(f"{k} ×{v}" for k, v in why.most_common())
            else:
                w = os.path.join(work, p["id"] + "-" + ff)
                os.makedirs(w, exist_ok=True)
                m_explicit, _ = missing_terms(xyz, ff, True, w)
                m_auto, err = missing_terms(xyz, ff, False, w)
                row["missing_explicit"] = len(m_explicit) if m_explicit is not None else None
                row["missing"] = len(m_auto) if m_auto is not None else None
                row["missing_terms"] = m_auto or []
                row["auto_filled"] = sorted(set(m_explicit or []) - set(m_auto or []))
            # DL_FIELD
            dd = os.path.join(dlf.get(ff, ""), f"{p['id']}__{ff}")
            dl_data, cp_data = os.path.join(dd, "dlfield", "lammps1.data"), os.path.join(dd, "caps", "system.data")
            if os.path.exists(dl_data) and not os.path.exists(cp_data) and not row["untyped"]:
                # compare_polymers.py stops at a missing term: CAPS's file with that term zero (reported in missing_terms)
                os.makedirs(os.path.dirname(cp_data), exist_ok=True)
                run([CAPS, "ff", "apply", xyz, "--ff", os.path.join(ROOT, "data", "forcefields", ff + ".json"), "--charges", "auto",
                     "--allow-missing", "-o", cp_data, "--lammps-input", os.path.join(dd, "caps", "system.in")])
            if not os.path.exists(dl_data) or not read_data(dl_data)["atoms"]:
                so = os.path.join(dd, "dlfield", "dl_field.stdout")
                err = [l.strip() for l in open(so).read().splitlines() if l.strip().startswith("Error") or "not available" in l] if os.path.exists(so) else []
                row["dlfield"] = "failed: " + (err[0] if err else "no output")
            elif not os.path.exists(cp_data):
                row["dlfield"] = "CAPS wrote no data file (a missing term; see missing_terms)"
            else:
                D, C = read_data(dl_data), read_data(cp_data)
                ids = sorted(C["atoms"])
                row["dl_types_same"] = sum(1 for i in ids if lab(D, i) == lab(C, i))
                row["dl_type_diffs"] = sorted({f"{lab(D, i)}→{lab(C, i)}" for i in ids if lab(D, i) != lab(C, i)})
                row["dl_charge_max_diff"] = max(abs(D["atoms"][i][1] - C["atoms"][i][1]) for i in ids)
                row["dl_term_diffs"] = dict(compare_terms(D, C, cls, oop=(cls == "class1")))
                # each atom's Lennard-Jones parameters (pair_coeff i i of DL_FIELD's input, CAPS's Pair Coeffs)
                lj_d = compare_polymers.lj_per_atom(dl_data, os.path.join(dd, "dlfield", "lammps.in"))
                lj_c = compare_polymers.lj_per_atom(cp_data, os.path.join(dd, "caps", "system.in"))
                for i, (x, y) in enumerate(zip(lj_d, lj_c)):
                    if x and y and not close(x, y) and not (x[0] == 0 and y[0] == 0):   # no well: its sigma means nothing
                        k = f"Lennard-Jones {lab(C, ids[i])}: DL_FIELD {fmt(x)} CAPS {fmt(y)}"
                        row["dl_term_diffs"][k] = row["dl_term_diffs"].get(k, 0) + 1
                ed = lammps_energy(os.path.join(dd, "dlfield"), "lammps.in")
                if cls == "class2":
                    without_cross_terms(cp_data, os.path.join(dd, "caps", "nocross.data"))
                    ec = lammps_energy(os.path.join(dd, "caps"), "system.in", "system.data", "nocross.data")
                else:
                    ec = lammps_energy(os.path.join(dd, "caps"), "system.in")
                row["dl_energy"], row["dl_energy_caps"] = ed, ec
                if ed and ec:
                    keys = ("bond", "angle", "dihedral") + (("improper",) if cls == "class1" else ()) + ("vdw",)
                    row["dl_rel"] = {k: abs(ed[k] - ec[k]) / max(1.0, abs(ed[k])) for k in keys}
            # msi2lmp
            m = msi.get((p["id"], frc))
            if m is None:
                row["msi2lmp"] = "not run"
            elif "error" in m:
                row["msi2lmp"] = m["error"]
            else:
                row["msi2lmp"] = "ok" if m["ok"] else "DIFFERS"
                row["msi_rel"] = {k: abs(m["caps"][k] - m["msi2lmp"][k]) / max(1.0, abs(m["msi2lmp"][k])) for k in m["msi2lmp"]}
                row["msi_energy"], row["msi_energy_caps"] = m["msi2lmp"], m["caps"]
                row["msi_note"] = m.get("note", "")
            rows.append(row)
            print(f"{p['id']} {ff:12s} typed {row['typed']}/{row['atoms']} missing {row.get('missing')} "
                  f"DL {row.get('dlfield') or (str(row['dl_types_same']) + ' types, ' + str(sum(row['dl_term_diffs'].values())) + ' terms differ')} "
                  f"msi2lmp {row['msi2lmp']}", flush=True)
    for r in rows:
        d = r.get("dl_term_diffs") or {}
        r["dl_unexplained"] = sorted(k for k in d if not explain(r["forcefield"], k))
        r["dl_explanations"] = sorted({explain(r["forcefield"], k) for k in d} - {None})
    json.dump(rows, open(os.path.join(outdir, "class2_polymers.json"), "w"), indent=1)
    write_tables(rows, outdir)
    shutil.rmtree(work, ignore_errors=True)
    return rows


E_DL = ("bond", "angle", "dihedral", "improper", "vdw")
E_MSI = ("ebond", "eangle", "edihed", "eimp", "evdwl", "ecoul")


def write_tables(rows, outdir):
    head = ["id", "polymer", "atoms", "typed", "untyped", "untyped_reason", "missing_terms_explicit", "missing_terms",
            "auto_parameters", "dlfield", "dl_types_same", "dl_type_diffs", "dl_charge_max_diff", "dl_term_groups_differ",
            "dl_term_groups_unexplained"] + [f"dl_rel_{k}" for k in E_DL] + ["msi2lmp"] + [f"msi_rel_{k}" for k in E_MSI] + ["msi_note"]
    summary = []
    for ff in FFS:
        mine = [r for r in rows if r["forcefield"] == ff]
        with open(os.path.join(outdir, f"{ff.split('-')[0]}_polymers.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(head)
            for r in mine:
                dl = r.get("dl_rel") or {}
                ms = r.get("msi_rel") or {}
                w.writerow([r["id"], r["polymer"], r["atoms"], r["typed"], r["untyped"], r.get("untyped_reason", ""),
                            r.get("missing_explicit", ""), r.get("missing", ""), len(r.get("auto_filled", [])),
                            r.get("dlfield", "compared"), r.get("dl_types_same", ""), " ".join(r.get("dl_type_diffs", [])),
                            f"{r['dl_charge_max_diff']:.1e}" if "dl_charge_max_diff" in r else "",
                            len(r.get("dl_term_diffs") or {}), len(r.get("dl_unexplained") or [])] +
                           [f"{dl[k]:.1e}" if k in dl else "" for k in E_DL] + [r["msi2lmp"] if r["msi2lmp"] in ("ok", "DIFFERS") else "not run: " + r["msi2lmp"][:120]] +
                           [f"{ms[k]:.1e}" if k in ms else "" for k in E_MSI] + [r.get("msi_note", "")])
        typed = [r for r in mine if not r["untyped"]]
        dlc = [r for r in mine if "dl_types_same" in r]
        msi = [r for r in mine if r["msi2lmp"] in ("ok", "DIFFERS")]
        # DL_FIELD writes its coefficients with six significant digits: terms equal to 1e-4 (relative) are the same parameters
        dl_eq = [r for r in dlc if r.get("dl_rel") and max(r["dl_rel"].values()) < DL_TOL]
        summary.append({
            "forcefield": ff, "polymers": len(mine), "fully_typed": len(typed), "untyped_polymers": len(mine) - len(typed),
            "atoms_typed": sum(r["typed"] for r in mine), "atoms": sum(r["atoms"] for r in mine),
            "complete_no_missing_terms": sum(1 for r in typed if r.get("missing") == 0),
            "with_auto_parameters": sum(1 for r in typed if r.get("auto_filled")),
            "dlfield_compared": len(dlc), "dlfield_failed": sum(1 for r in mine if str(r.get("dlfield", "")).startswith("failed")),
            "dl_atoms_same_type": sum(r["dl_types_same"] for r in dlc), "dl_atoms": sum(r["atoms"] for r in dlc),
            "dl_polymers_same_types": sum(1 for r in dlc if r["dl_types_same"] == r["atoms"]),
            "dl_polymers_same_terms": sum(1 for r in dlc if not r["dl_term_diffs"]),
            "dl_polymers_energy_equal": len(dl_eq),
            "dl_polymers_differ_explained": sum(1 for r in dlc if r["dl_term_diffs"] and not r["dl_unexplained"]),
            "dl_polymers_unexplained": sum(1 for r in dlc if r["dl_unexplained"]),
            "dl_max_charge_diff": max((r["dl_charge_max_diff"] for r in dlc), default=0),
            "msi2lmp_compared": len(msi), "msi2lmp_agree": sum(1 for r in msi if r["msi2lmp"] == "ok"),
            "msi2lmp_exact": sum(1 for r in msi if max(r["msi_rel"].values()) < 1e-6),
            "msi2lmp_max_rel_exact": max((max(r["msi_rel"].values()) for r in msi if max(r["msi_rel"].values()) < 1e-6), default=0),
        })
    with open(os.path.join(outdir, "class2_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0]))
        w.writeheader()
        w.writerows(summary)
    for s in summary:
        print(json.dumps(s))


if __name__ == "__main__":
    main()
