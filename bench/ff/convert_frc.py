#!/usr/bin/env python3
"""CVFF, PCFF and COMPASS for CAPS from their original BIOVIA (Accelrys) force-field files, the .frc files LAMMPS's
msi2lmp tool reads (lammps/tools/msi2lmp/frc_files: cvff.frc, pcff.frc, compass_published.frc).

The explicit parameters are read as msi2lmp reads them, so that CAPS assigns what msi2lmp assigns:
  · the first section of each kind (#quartic_bond, #torsion_3 …; for CVFF #quadratic_bond, #torsion_1 …), its
    entries until the next '#', types cut to 4 characters;
  · an entry for the same ordered types with a higher version replaces the earlier one in place;
  · a missing right-hand half (bond-angle, end-bond-torsion, angle-torsion) repeats the left-hand one;
  · per term: the atoms' own types first, then their equivalences (CAPS "equivalence": "fallback"); exact entries
    before wildcard ones, the first in the file winning (written here in the order CAPS's last-match rule needs);
  · class II cross terms each looked up on their own ("cross_rules"), reference lengths and angles from the assigned
    bonds and angles; Wilson out-of-plane and angle-angle terms as msi2lmp builds them ("oop_scheme": "msi2lmp");
  · CVFF: a torsion with a wildcard end spreads its barrier over the torsions about the central bond.

Materials Studio's automatic parameters (the *_auto sections, by #auto_equivalence) are kept as the fallback CAPS
reports as automatic when no explicit entry exists (msi2lmp stops there instead). Charges: #bond_increments.

usage: convert_frc.py [FRC_DIR]   (default ~/lammps/tools/msi2lmp/frc_files) → data/forcefields/{cvff,pcff,compass}-frc.json
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/lammps/tools/msi2lmp/frc_files")


def sections(path):
    """Every section of the file: (keyword, label, [data lines]) in file order."""
    out, cur = [], None
    for raw in open(path, encoding="latin-1"):
        line = raw.rstrip("\n")
        if line.startswith("#"):
            w = line.replace("(", " (").split()
            cur = (w[0], " ".join(w[1:]), [])
            out.append(cur)
            continue
        if cur is None:
            continue
        s = line.strip()
        if not s or s[0] in "!>@":
            continue
        cur[2].append(s.split())
    return out


def first(secs, keyword, label=None):
    """msi2lmp reads the first section with this keyword (a label narrows it to the *_auto sections)."""
    for k, lab, lines in secs:
        if k == keyword and (label is None or lab.split()[0].endswith(label)):
            return lines
    return []


def entries(lines, nmem, npar):
    """(version, reference, types, params) as msi2lmp stores them: higher versions replace, missing halves repeat."""
    out = []
    for w in lines:
        try:
            ver, ref = float(w[0]), w[1]
        except ValueError:
            continue
        types = [t[:4] for t in w[2:2 + nmem]]
        vals = []
        for x in w[2 + nmem:2 + nmem + npar]:
            try:
                vals.append(float(x))
            except ValueError:
                break
        if not vals and npar:
            continue
        given = len(vals)
        while len(vals) < npar:
            vals.append(vals[len(vals) - given])
        for e in out:
            if e[2] == types:
                if ver > e[0]:
                    e[0], e[1], e[3] = ver, ref, vals
                break
        else:
            out.append([ver, ref, types, vals])
    return out


def wild(t):
    return t.startswith("*")


def ordered(es, make):
    """Rules in the order CAPS's last-match lookup needs for msi2lmp's choice: wildcard entries first, exact entries
    last, each group reversed (the first entry in the file is the last rule, so it wins)."""
    w = [e for e in es if any(wild(t) for t in e[2])]
    x = [e for e in es if not any(wild(t) for t in e[2])]
    return [make(e) for e in w[::-1] + x[::-1]]


def rule(e, params, style=None, name=None):
    # a literal * or ? in a type name (CVFF's o*, h* water) is escaped: only a leading * is msi2lmp's wildcard
    r = {"name": name or "-".join(e[2]), "match": ["*" if wild(t) else t.replace("\\", "\\\\").replace("*", "\\*").replace("?", "\\?")
                                                   for t in e[2]]}
    if style:
        r["style"] = style
    r["params"] = params
    r["comment"] = f"ver {e[0]:g}, ref {e[1]}"
    return r


def convert(frc, cls, name, version, references, typing, out_name, src=None):
    secs = sections(os.path.join(src or SRC, frc))
    auto = "_auto"
    atoms = entries(first(secs, "#atom_types"), 1, 0)
    # atom types: mass, element (columns after the type)
    info = {}
    for w in first(secs, "#atom_types"):
        try:
            float(w[0])
        except ValueError:
            continue
        info.setdefault(w[2][:4], (float(w[3]), w[4], " ".join(w[6:])))
    eq = {e[2][0]: e for e in entries(first(secs, "#equivalence"), 6, 0)}
    eqw = {w[2][:4]: [t[:4] for t in w[2:8]] for w in first(secs, "#equivalence") if len(w) >= 8}
    aeq = {w[2][:4]: [t[:4] for t in w[2:12]] for w in first(secs, "#auto_equivalence") if len(w) >= 12}
    types = []
    for t in dict.fromkeys(e[2][0] for e in atoms):
        mass, el, comment = info[t]
        q = {}
        if t in eqw:
            e = eqw[t]
            q.update(vdw=e[1], bond=e[2], angle=e[3], dihedral=e[4], improper=e[5])
        if t in aeq:
            a = aeq[t]
            q.update(increment=a[2], auto_bond=a[3], auto_angle_end=a[4], auto_angle_apex=a[5], auto_torsion_end=a[6],
                     auto_torsion_center=a[7], auto_oop_end=a[8], auto_oop_center=a[9])
        types.append({"name": t, "element": el if el not in ("*", "") else "", "mass": mass, "description": comment,
                      "equivalence": {k: v for k, v in q.items() if v != t}})
    ff = {"format": "caps-forcefield", "format_version": 1, "name": name, "version": version,
          "source": f"{frc} (BIOVIA / Accelrys, as distributed with LAMMPS's msi2lmp tool)" if src is None else
                    f"{frc} (INTERFACE force field 1.5, Heinz group, bionanostructures.com)", "references": references,
          "units": "real", "mixing": "sixthpower" if cls == 2 else "geometric",
          "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 12,
          "equivalence": "fallback"}
    if cls == 2:
        ff["styles"] = {"pair": "lj/class2/coul/long", "bond": "class2", "angle": "class2", "dihedral": "class2", "improper": "class2"}
        ff["oop_scheme"] = "msi2lmp"
        ff["improper_order"] = "center2_sorted"
    else:
        ff["styles"] = {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "harmonic", "improper": "cvff"}
        ff["wildcard_torsion_scaling"] = "msi2lmp"
        ff["oop_scheme"] = "msi2lmp"
        ff["improper_order"] = "center2_sorted"
    ff["atom_types"] = types
    # van der Waals: 9-6 (r, eps) or 12-6 (A, B)
    if cls == 2:
        ff["pairs"] = ordered(entries(first(secs, "#nonbond"), 1, 2), lambda e: rule(e, [e[3][1], e[3][0]]))
    else:
        def ab(e):
            A, B = e[3]
            return rule(e, [B * B / (4 * A), (A / B) ** (1 / 6)] if A and B else [0.0, 0.0])
        ff["pairs"] = ordered(entries(first(secs, "#nonbond"), 1, 2), ab)
    if cls == 2:
        ff["bonds"] = ordered(entries(first(secs, "#quartic_bond"), 2, 4), lambda e: rule(e, e[3], "class2"))
        ff["angles"] = ordered(entries(first(secs, "#quartic_angle"), 3, 4), lambda e: rule(e, e[3], "class2"))
        ff["dihedrals"] = ordered(entries(first(secs, "#torsion_3"), 4, 6), lambda e: rule(e, e[3], "class2"))
        ff["impropers"] = ordered(entries(first(secs, "#wilson_out_of_plane"), 4, 2), lambda e: rule(e, e[3], "class2"))
        cr = {}
        for g, key, nm, npar in (("bb", "#bond-bond", 3, 1), ("ba", "#bond-angle", 3, 2), ("mbt", "#middle_bond-torsion_3", 4, 3),
                                 ("ebt", "#end_bond-torsion_3", 4, 6), ("at", "#angle-torsion_3", 4, 6),
                                 ("aat", "#angle-angle-torsion_1", 4, 1), ("bb13", "#bond-bond_1_3", 4, 1), ("aa", "#angle-angle", 4, 1)):
            cr[g] = ordered(entries(first(secs, key), nm, npar), lambda e: rule(e, e[3]))
        ff["cross_rules"] = cr
    else:
        ff["bonds"] = ordered(entries(first(secs, "#quadratic_bond"), 2, 2), lambda e: rule(e, [e[3][1], e[3][0]], "harmonic"))
        ff["angles"] = ordered(entries(first(secs, "#quadratic_angle"), 3, 2), lambda e: rule(e, [e[3][1], e[3][0]], "harmonic"))

        def t1(e, style="harmonic"):
            K, n, phi0 = e[3]
            if phi0 not in (0.0, 180.0):
                raise SystemExit(f"{frc}: torsion {e[2]} has phi0 {phi0}: not a CVFF form")
            return rule(e, [K, 1 if phi0 == 0 else -1, int(n)], style)
        ff["dihedrals"] = ordered(entries(first(secs, "#torsion_1"), 4, 3), t1)
        ff["impropers"] = ordered(entries(first(secs, "#out_of_plane"), 4, 3), lambda e: t1(e, "cvff"))
    # Materials Studio's automatic parameters (by #auto_equivalence), used only when no explicit entry exists
    def auto_rule(e, params, style):
        r = rule(e, params, style, "auto " + " ".join(e[2]))
        r["comment"] = "auto, " + r["comment"]
        return r
    ab_ = entries(first(secs, "#quadratic_bond", auto), 2, 2)
    aa_ = entries(first(secs, "#quadratic_angle", auto), 3, 2)
    at_ = entries(first(secs, "#torsion_1", auto), 4, 3)
    if cls == 2:
        ff["auto_bonds"] = ordered(ab_, lambda e: auto_rule(e, [e[3][0], e[3][1], 0, 0], "class2"))
        ff["auto_angles"] = ordered(aa_, lambda e: auto_rule(e, [e[3][0], e[3][1], 0, 0], "class2"))

        def t1c2(e):
            # Kphi [1 + cos(nφ − φ0)] = Kphi [1 − cos(nφ − (φ0 + 180))]: LAMMPS's class2 form, K and φ at position n
            K, n, phi0 = e[3]
            p = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            if int(n) not in (1, 2, 3):
                return auto_rule(e, [1, K, int(n), phi0], "fourier")
            p[2 * (int(n) - 1)] = K
            p[2 * (int(n) - 1) + 1] = (phi0 + 180.0) % 360.0
            return auto_rule(e, p, "class2")
        ff["auto_dihedrals"] = ordered(at_, t1c2)
    else:
        ff["auto_bonds"] = ordered(ab_, lambda e: auto_rule(e, [e[3][1], e[3][0]], "harmonic"))
        ff["auto_angles"] = ordered(aa_, lambda e: auto_rule(e, [e[3][1], e[3][0]], "harmonic"))
        ff["auto_dihedrals"] = ordered(at_, lambda e: auto_rule(e, [e[3][0], 1 if e[3][2] == 0 else -1, int(e[3][1])], "harmonic"))
    ff["bond_increments"] = ordered(entries(first(secs, "#bond_increments"), 2, 2), lambda e: rule(e, e[3]))
    ff["notes"] = [f"converted by CAPS from {frc}: explicit parameters as LAMMPS's msi2lmp assigns them; Materials Studio's "
                   "automatic parameters where no explicit entry exists (reported as automatic)"]
    if cls == 1:
        ff["notes"].append("CVFF as msi2lmp writes it for LAMMPS: harmonic bonds and angles, no cross terms, no Morse bonds")
    # the typing rules of the CAPS file for this force field, less those for types the .frc file does not have
    tj = json.load(open(os.path.normpath(os.path.join(ROOT, "data", "forcefields", typing))))
    have = {t["name"] for t in types}
    # moltemplate spells the characters its names cannot hold: c3prime is c3', o1=star is o1=*
    renamed = 0
    for r in tj["rules"]:
        t = r["type"].replace("prime", "'").replace("star", "*")
        if t != r["type"] and t in have:
            r["type"] = t
            renamed += 1
    dropped = sorted({r["type"] for r in tj["rules"] if r["type"] not in have})
    if dropped or renamed:
        tj["rules"] = [r for r in tj["rules"] if r["type"] in have]
        tj["forcefield"] = name
        tj["description"] = tj.get("description", "") + f" (for {frc}: rules for {', '.join(dropped)}, which it does not define, left out)"
        typing = "../typing/" + out_name.replace(".json", ".typing.json")
        json.dump(tj, open(os.path.normpath(os.path.join(ROOT, "data", "forcefields", typing)), "w"), indent=1, ensure_ascii=False)
    ff["typing"] = typing
    ff["status"] = "validated"
    dst = os.path.join(ROOT, "data", "forcefields", out_name)
    json.dump(ff, open(dst, "w"), indent=1, ensure_ascii=False)
    n = {k: len(v) for k, v in ff.items() if isinstance(v, list) and k not in ("references", "special_lj", "special_coul", "notes")}
    if "cross_rules" in ff:
        n["cross"] = {k: len(v) for k, v in ff["cross_rules"].items()}
    print(f"{frc} -> {os.path.relpath(dst, ROOT)}: {n}")


convert("pcff.frc", 2, "PCFF (pcff.frc, full class II)", "cff91 / pcff.frc 4.0",
        ["H. Sun, S. J. Mumby, J. R. Maple, A. T. Hagler, J. Am. Chem. Soc. 116, 2978 (1994)",
         "H. Sun, Macromolecules 28, 701 (1995)"], "../typing/pcff.typing.json", "pcff-frc.json")
convert("cvff.frc", 1, "CVFF (cvff.frc)", "cvff.frc 2.4",
        ["P. Dauber-Osguthorpe et al., Proteins 4, 31 (1988)"], "../typing/cvff.typing.json", "cvff-frc.json")
IFF = os.path.expanduser("~/iff-ref/INTERFACE_FF_1_5/FORCE_FIELDS")
if os.path.isdir(IFF):   # the INTERFACE force field: PCFF and CVFF with the inorganic phases (clays, silica, metals, cement …)
    convert("pcff_interface_v1_5.frc", 2, "INTERFACE (IFF 1.5, PCFF)", "IFF 1.5 on pcff.frc",
            ["H. Heinz, T.-J. Lin, R. K. Mishra, F. S. Emami, Langmuir 29, 1754 (2013)"], "../typing/pcff.typing.json", "iff-pcff.json", IFF)
    convert("cvff_interface_v1_5.frc", 1, "INTERFACE (IFF 1.5, CVFF)", "IFF 1.5 on cvff.frc",
            ["H. Heinz, T.-J. Lin, R. K. Mishra, F. S. Emami, Langmuir 29, 1754 (2013)"], "../typing/cvff.typing.json", "iff-cvff.json", IFF)
convert("compass_published.frc", 2, "COMPASS (compass_published.frc, full class II)", "compass_published.frc 1.1",
        ["H. Sun, J. Phys. Chem. B 102, 7338 (1998)"], "../typing/compass-published-moltemplate.typing.json", "compass-frc.json")


# the library catalogue: the .frc conversions listed before the other CVFF / PCFF / COMPASS files
CHECKED = {
    "iff-pcff": "the parameters msi2lmp assigns from pcff_interface_v1_5.frc for models of IFF's own database (pyrophyllite, "
                "kaolinite, mica, montmorillonite, cristobalite and hydroxylated silica, gold and aluminium surfaces, "
                "hydroxyapatite, gypsum, a hydrated C3A surface, PEO), energies term by term in LAMMPS (bench/ff/check_msi2lmp.py)",
    "iff-cvff": "the parameters msi2lmp assigns from cvff_interface_v1_5.frc for models of IFF's own database (clays, silica, "
                "metals, hydroxyapatite, PEO), energies term by term in LAMMPS (bench/ff/check_msi2lmp.py)",
    "pcff-frc": "the parameters msi2lmp assigns from pcff.frc for LAMMPS's msi2lmp test structures (water, ethane, benzene, "
                "naphthalene, a carbon nanotube, hydroxyapatite), energies term by term in LAMMPS; CAPS's own energies and "
                "forces equal LAMMPS's (bench/ff/check_msi2lmp.py, check_data_lammps.py)",
    "cvff-frc": "the parameters msi2lmp assigns from cvff.frc for LAMMPS's msi2lmp test structures (crambin, nylon, a carbon "
                "nanotube, water, aromatics): equal, apart from msi2lmp counting a wildcard torsion's neighbours by type; CAPS's "
                "own energies and forces equal LAMMPS's (bench/ff/check_msi2lmp.py, check_data_lammps.py)",
    "compass-frc": "the parameters msi2lmp assigns from compass_published.frc for LAMMPS's msi2lmp test structures (ethane, "
                   "benzene, naphthalene, a carbon nanotube, hydrogen): equal (msi2lmp leaves out bond-bond-1-3 terms without a cp "
                   "type); CAPS's own energies and forces equal LAMMPS's (bench/ff/check_msi2lmp.py, check_data_lammps.py)",
}
cat_p = os.path.join(ROOT, "data", "forcefields", "catalogue.json")
cat = json.load(open(cat_p))
fs = [e for e in cat["forcefields"] if e["id"] not in CHECKED]
at = next((k for k, e in enumerate(fs) if e["id"] in ("pcff", "compass", "cvff")), len(fs))
new = []
for fid in ("pcff-frc", "cvff-frc", "compass-frc", "iff-pcff", "iff-cvff"):
    if not os.path.exists(os.path.join(ROOT, "data", "forcefields", fid + ".json")):
        continue
    ff = json.load(open(os.path.join(ROOT, "data", "forcefields", fid + ".json")))
    new.append({"id": fid, "name": ff["name"], "version": ff["version"], "references": ff["references"],
                "origin": "INTERFACE .frc" if fid.startswith("iff") else "BIOVIA .frc",
                "source_file": ff["source"].split(" ")[0], "status": "validated", "notes": "validated against " + CHECKED[fid],
                "file": fid + ".json",
                "typing": {"rules": ff["typing"].replace("../", ""),
                           "evidence": ("the typing rules of CAPS's PCFF / CVFF library for organic atoms; inorganic types from the "
                                        "model's .car file (IFF's model database)") if fid.startswith("iff") else
                                       "the typing rules of CAPS's " + ff["name"].split(" ")[0] + " library, for the types the .frc file defines"},
                "counts": {k: len(ff[k]) for k in ("atom_types", "pairs", "bonds", "angles", "dihedrals", "impropers")}})
cat["forcefields"] = fs[:at] + new + fs[at:]
json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)
print(f"catalogue: {', '.join(e['id'] for e in new)} listed before the other CVFF / PCFF / COMPASS files")
