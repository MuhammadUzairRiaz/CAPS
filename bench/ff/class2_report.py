#!/usr/bin/env python3
"""PCFF, COMPASS and CVFF (from their .frc files) on the polymer library: one table per force field, one row per polymer.

For every oligomer (a directory of .xyz files, one per polymer of data/polymers/library.json) and each of pcff-frc,
compass-frc and cvff-frc:

  · typing: atoms CAPS types and the atoms no rule matches; for COMPASS, why (the chemistry compass_published.frc has
    no type for, named by the atom's PCFF type);
  · missing terms: with the .frc file's explicit tables only and with its automatic parameters (caps ff apply
    --allow-missing lists them);
  · against msi2lmp (check_msi2lmp.py --polymers --json): the LAMMPS energy per term on msi2lmp's data file and on
    CAPS's for the same typed structure.

usage: class2_report.py --olig DIR --msi MSI.json --out DIR [--only P0]
Writes OUT/{pcff,compass,cvff}_polymers.csv, OUT/class2_polymers.json and OUT/class2_summary.csv; prints the totals.
"""
import collections, csv, json, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
FFS = {"pcff-frc": "pcff", "compass-frc": "compass_published", "cvff-frc": "cvff"}   # CAPS library → msi2lmp's .frc
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


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


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
                  f"msi2lmp {row['msi2lmp']}", flush=True)
    json.dump(rows, open(os.path.join(outdir, "class2_polymers.json"), "w"), indent=1)
    write_tables(rows, outdir)
    shutil.rmtree(work, ignore_errors=True)
    return rows


E_MSI = ("ebond", "eangle", "edihed", "eimp", "evdwl", "ecoul")


def write_tables(rows, outdir):
    head = ["id", "polymer", "atoms", "typed", "untyped", "untyped_reason", "missing_terms_explicit", "missing_terms",
            "auto_parameters", "msi2lmp"] + [f"msi_rel_{k}" for k in E_MSI] + ["msi_note"]
    summary = []
    for ff in FFS:
        mine = [r for r in rows if r["forcefield"] == ff]
        with open(os.path.join(outdir, f"{ff.split('-')[0]}_polymers.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(head)
            for r in mine:
                ms = r.get("msi_rel") or {}
                w.writerow([r["id"], r["polymer"], r["atoms"], r["typed"], r["untyped"], r.get("untyped_reason", ""),
                            r.get("missing_explicit", ""), r.get("missing", ""), len(r.get("auto_filled", []))] +
                           [r["msi2lmp"] if r["msi2lmp"] in ("ok", "DIFFERS") else "not run: " + r["msi2lmp"][:120]] +
                           [f"{ms[k]:.1e}" if k in ms else "" for k in E_MSI] + [r.get("msi_note", "")])
        typed = [r for r in mine if not r["untyped"]]
        msi = [r for r in mine if r["msi2lmp"] in ("ok", "DIFFERS")]
        summary.append({
            "forcefield": ff, "polymers": len(mine), "fully_typed": len(typed), "untyped_polymers": len(mine) - len(typed),
            "atoms_typed": sum(r["typed"] for r in mine), "atoms": sum(r["atoms"] for r in mine),
            "complete_no_missing_terms": sum(1 for r in typed if r.get("missing") == 0),
            "with_auto_parameters": sum(1 for r in typed if r.get("auto_filled")),
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
