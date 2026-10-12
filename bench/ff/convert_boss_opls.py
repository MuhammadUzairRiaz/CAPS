#!/usr/bin/env python3
"""A CAPS force-field file from Jorgensen's own OPLS-AA parameter files in BOSS format: the .par (atom types with their
charges and Lennard-Jones, the Fourier torsions) and the .sb (bond stretching, angle bending), as published with
W. L. Jorgensen et al., J. Phys. Chem. B 128, 250 (2024), supporting information jp3c06602_si_002 / _003.

  convert_boss_opls.py FILE.par FILE.sb OUT.json [--name NAME] [--typing ../typing/ID.typing.json]

Read:
  .par  "TYP AN AT CHARGE SIGMA EPSILON comment": type number, atomic number, the bond class AT that names the type in
        the .sb and in the torsions, its charge, sigma (A) and epsilon (kcal/mol); after "The following are the Fourier
        coefficients": "type V1 V2 V3 V4" and the four classes in columns 48-60 (X, Y, Z wildcards); a line marked
        "improper" is an improper torsion with its centre second
  .sb   bonds "A -B  K r0" (columns 1-5), a blank line, angles "A -B -C  K theta0" (columns 1-8)
CAPS's form (real units, as LAMMPS's opls / cvff / harmonic styles take them):
  types      opls_NNN with its own charge; the bond class names it in bonds, angles, torsions and impropers
  pairs      epsilon, sigma; geometric mixing; 1-4 x 0.5 for both
  dihedrals  opls V1 V2 V3 V4 as given
  impropers  cvff K = V2 / 2, d = -1, n = 2 (V2 / 2 [1 - cos 2 phi]); the centre (second in BOSS) written first, the
             other three in index order, as LigParGen writes them
Several lines may name the same classes: the first is the general one (Jorgensen's convention, as moltemplate's
oplsaa2lt.py records it), the later ones belong to particular molecules (butane, benzoic acids ...) and are listed in
the notes, not applied. "allenes improper" lines name four nearly collinear atoms and are left out (as moltemplate does).
"""
import json, re, sys

PAR, SB, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
ELEMENT = ["X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
           "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
           "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba"]
MASS = {"H": 1.008, "He": 4.0026, "Li": 6.941, "Be": 9.012, "B": 10.811, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
        "Ne": 20.18, "Na": 22.99, "Mg": 24.305, "Al": 26.982, "Si": 28.086, "P": 30.974, "S": 32.065, "Cl": 35.453, "Ar": 39.948,
        "K": 39.098, "Ca": 40.078, "Fe": 55.845, "Cu": 63.546, "Zn": 65.38, "Br": 79.904, "Kr": 83.798, "Rb": 85.468, "Sr": 87.62,
        "I": 126.904, "Xe": 131.293, "Cs": 132.905, "Ba": 137.327, "X": 0.0}


def cls(t):
    """a class as a CAPS pattern: X / Y / Z match any class; a ? is BOSS's one-character wildcard in its two-character
    class field ("C?" is C, CT, CA, C* ...: classes have at most two characters, so CAPS's glob C*); a literal * in a
    class name (C*, N*) is escaped"""
    t = t.strip()
    if t in ("X", "Y", "Z"):
        return "*"
    return t.replace("*", "\\*").replace("?", "*")


types, pairs, notes = [], [], []
lines = open(PAR, errors="replace").read().splitlines()
k = 0
while k < len(lines) and "The following are the Fourier coefficients" not in lines[k]:
    l = lines[k]
    k += 1
    if l.startswith("#") or len(l) < 8 or not l[5:7].strip().isdigit():
        continue
    w = l.split()
    try:
        num, an, at, q, sig, eps = int(w[0]), int(w[1]), w[2], float(w[3]), float(w[4]), float(w[5])
    except (ValueError, IndexError):
        continue
    el = ELEMENT[an] if 0 <= an < len(ELEMENT) else "X"
    name = "opls_%03d" % num
    types.append({"name": name, "element": el, "mass": MASS.get(el, 0.0), "charge": q, "description": " ".join(w[6:]),
                  "equivalence": {"bond": at, "angle": at, "dihedral": at, "improper": at}})
    pairs.append({"name": name, "match": [name], "params": [eps, sig]})

# Where the first line of a group belongs to one molecule, the general line further down is applied instead (the file
# itself marks a line "repeated here so taken first": first wins otherwise)
GENERAL = {("F", "CT", "CT", "HC"): "alkyl fluoride",                 # first: trifluoroethanol
           ("CK", "N*", "CT", "OS"): "imidazoles, indoles, purines",  # first: a nucleoside's
           ("CK", "NA", "CT", "CT"): "imidazoles, indoles, purines"}
dihedrals, impropers, seen_d, seen_i, alternatives = [], [], {}, {}, []
for l in lines[k + 1:]:
    if l.startswith("#") or l.startswith("Type") or len(l) < 48:
        continue
    w = l[:47].split()
    if len(w) != 5:
        continue
    defn = l[47:60].strip()
    if defn in ("Dummy", "", "Harmonic Rest"):
        continue
    if defn.endswith("P in"):   # a line whose comment runs into the class columns
        defn = defn[:-3]
    names = [x.strip() for x in defn.split("-")]
    if len(names) != 4:
        continue
    V = [float(x) for x in w[1:5]]
    comment = l[60:].strip()
    if "improper" in comment:
        if comment == "allenes improper":
            continue
        key = tuple(names)
        if key in seen_i:
            continue
        seen_i[key] = True
        # the centre (second in BOSS) first, the other three in index order: as LigParGen, Jorgensen's own tool, writes
        # them (moltemplate's oplsaa2lt.py found its cenIsortJKL order reproduces LigParGen's LAMMPS files)
        impropers.append({"name": "-".join(names), "match": [cls(x) for x in (names[1], names[0], names[2], names[3])], "style": "cvff", "params": [V[1] / 2, -1, 2],
                          "comment": f"improper, BOSS type {w[0]}" + (f": {comment}" if comment else "")})
        continue
    key = tuple(names) if tuple(names) <= tuple(names[::-1]) else tuple(names[::-1])
    general = next((c for k, c in GENERAL.items() if key in (k, k[::-1])), None)
    if key in seen_d:
        if general and comment.startswith(general) and not seen_d[key]["comment"].split(": ", 1)[-1].startswith(general):
            old = seen_d[key]
            alternatives.append(f"{old['name']} {old['params']} ({old['comment'].split(': ', 1)[-1]})")
            old.update({"name": "-".join(names), "match": [cls(x) for x in names], "params": V, "comment": f"BOSS type {w[0]}: {comment}"})
        elif V != seen_d[key]["params"]:
            alternatives.append(f"{'-'.join(names)} {V} ({comment})")
        continue
    d = {"name": "-".join(names), "match": [cls(x) for x in names], "style": "opls", "params": V,
         "comment": f"BOSS type {w[0]}" + (f": {comment}" if comment else "")}
    seen_d[key] = d
    dihedrals.append(d)
# CAPS takes the last rule that matches: X-wildcard rules first, then ?-wildcard ones, then fully named rules
dihedrals.sort(key=lambda d: 0 if "*" in d["match"] else 1 if any("*" in m.replace("\\*", "") for m in d["match"]) else 2)

bonds, angles, seen_b, seen_a = [], [], {}, {}
section = "bonds"
for l in open(SB, errors="replace").read().splitlines():
    if l.startswith("*") or l.startswith("#") or l.startswith('"""'):
        continue
    if not l.strip():
        if section == "bonds" and bonds:
            section = "angles"
        continue
    if section == "bonds":
        names, rest = [x.strip() for x in l[:5].split("-")], l[5:].split()
        if len(names) != 2 or len(rest) < 2:
            continue
        key = tuple(sorted(names))
        if key in seen_b:
            continue
        seen_b[key] = True
        bonds.append({"name": "-".join(names), "match": [cls(x) for x in names], "style": "harmonic",
                      "params": [float(rest[0]), float(rest[1])], **({"comment": " ".join(rest[2:])} if rest[2:] else {})})
    else:
        names, rest = [x.strip() for x in l[:8].split("-")], l[8:].split()
        if len(names) != 3 or len(rest) < 2:
            continue
        key = tuple(names) if names[0] <= names[2] else tuple(names[::-1])
        if key in seen_a:
            continue
        seen_a[key] = True
        angles.append({"name": "-".join(names), "match": [cls(x) for x in names], "style": "harmonic",
                       "params": [float(rest[0]), float(rest[1])], **({"comment": " ".join(rest[2:])} if rest[2:] else {})})

if arg("--supplement"):   # terms the files lack, each from the OPLS-AA paper that defines it (Jorgensen's files win)
    sup = json.load(open(arg("--supplement")))
    added = []
    for kind, lst, seen in (("bonds", bonds, seen_b), ("angles", angles, seen_a)):
        for t in sup.get(kind, []):
            m = t["match"]
            key = tuple(sorted(m)) if kind == "bonds" else (tuple(m) if m[0] <= m[2] else tuple(m[::-1]))
            if key in seen:
                continue
            seen[key] = True
            lst.append({"name": "-".join(m), "match": [cls(x) for x in m], "style": "harmonic", "params": t["params"],
                        "comment": f"{t['source']}: {t['reason']}"})
            added.append("-".join(m))
    if added:
        notes.append(f"{len(added)} terms the 2024 files lack, each from the OPLS-AA paper that defines it (bench/ff/open/oplsaa2024_supplement.json): "
                     + ", ".join(added))
if alternatives:
    notes.append(f"{len(alternatives)} torsions the .par gives again for particular molecules (the first, general line is "
                 "applied): " + "; ".join(alternatives))
notes.append("'allenes improper' lines (four nearly collinear atoms) are left out")
doc = {"format": "caps-forcefield", "format_version": 1, "name": arg("--name", "OPLS-AA 2024 (Jorgensen)"),
       "source": "BOSS parameter files of W. L. Jorgensen et al., J. Phys. Chem. B 128, 250 (2024), SI jp3c06602_si_002 (.par) "
                 "and jp3c06602_si_003 (.sb)", "units": "real",
       "references": ["W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)",
                      "W. L. Jorgensen et al., J. Phys. Chem. B 128, 250 (2024)"],
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"},
       "mixing": "geometric", "special_lj": [0, 0, 0.5], "special_coul": [0, 0, 0.5], "cutoff": 12,
       "improper_order": "center1_sorted", "equivalence": "replace",
       "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": impropers,
       "notes": ["converted by CAPS (bench/ff/convert_boss_opls.py) from Jorgensen's published BOSS files"] + notes,
       "version": "2024 (parameter files of October 2023)"}
if arg("--typing"):
    doc["typing"] = arg("--typing")
json.dump(doc, open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"{OUT}: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} dihedrals, {len(impropers)} impropers, "
      f"{len(alternatives)} molecule-specific torsions listed")
