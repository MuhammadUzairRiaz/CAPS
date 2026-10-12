#!/usr/bin/env python3
"""A CAPS force-field file from AMBER's own parameter files (AmberTools dat/leap/parm: gaff.dat, gaff2.dat, parm99.dat,
frcmod.*), which AmberTools' LICENSE states "have been placed in the public domain by their authors".

  convert_amber_dat.py OUT.json --name NAME --dat gaff2.dat [--frcmod frcmod.ff03 ...] [--reference "..."] [--version "..."]
                     [--typing ../typing/ID.typing.json] [--scale14 LJ COUL]

The .dat sections, in order and separated by blank lines: atom types (name, mass), the hydrophilic list, bonds, angles,
dihedrals, impropers, the 10-12 hydrogen-bond terms, nonbonded equivalences, and the nonbonded block (MOD4 RE: R*, epsilon);
a frcmod names its sections (MASS, BOND, ANGL, DIHE, IMPR, NONB) and overrides the .dat entry for the same types. CAPS's
form:

  pairs      epsilon, sigma = 2 R* / 2^(1/6)
  bonds      harmonic K (AMBER's K, energy K (r - r0)^2) and r0
  angles     harmonic K theta0
  dihedrals  fourier: the number of terms, then (PK / IDIVF, |PN|, PHASE) for each; a negative PN means more terms follow
  impropers  cvff K d n with d = cos(PHASE) (AMBER's improper is K [1 + cos(n phi - phase)], phase 0 or 180)
  X          wildcard, written * (the name keeps X)
  1-4        Lennard-Jones x 0.5, Coulomb x 1/1.2 (AMBER's SCNB 2, SCEE 1.2)
"""
import json, math, re, sys

arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
many = lambda k: [sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == k]
OUT = sys.argv[1]
SIG = 2 ** (1 / 6)
ELEMENTS = [("H", 1.008), ("C", 12.011), ("N", 14.007), ("O", 15.999), ("F", 18.998), ("Na", 22.99), ("Mg", 24.305), ("P", 30.974),
            ("S", 32.06), ("Cl", 35.45), ("K", 39.098), ("Ca", 40.078), ("Fe", 55.845), ("Cu", 63.546), ("Zn", 65.38),
            ("Br", 79.904), ("I", 126.904), ("Li", 6.94), ("Rb", 85.468), ("Cs", 132.905), ("Si", 28.085), ("B", 10.81),
            ("Se", 78.971), ("Ag", 107.868), ("Au", 196.967), ("Pt", 195.08), ("He", 4.0026), ("Ne", 20.18), ("Ar", 39.948)]

types, pairs, bonds, angles, dihedrals, impropers, notes = {}, {}, {}, {}, {}, {}, []
equiv = {}


def element_of(mass):
    return "X" if mass < 0.5 else min(ELEMENTS, key=lambda e: abs(e[1] - mass))[0]


def wild(t):
    return "*" if t.upper() == "X" else t


def keyfields(line, n):
    """the n type names of a bonded line ("c -c -c " fixed width 2 per name, '-' between) and the rest of the line"""
    width = 3 * n - 1
    names = [t.strip() for t in line[:width].split("-")]
    return names, line[width:].split()


def add_dihedral(names, rest, store, carry):
    """one dihedral line: IDIVF PK PHASE PN; a negative PN continues the same key on the next line"""
    idivf, pk, phase, pn = float(rest[0]), float(rest[1]), float(rest[2]), float(rest[3])
    k = tuple(names) if tuple(names) <= tuple(names[::-1]) else tuple(names[::-1])
    if k not in carry:   # the first line of this key in this file replaces any earlier file's terms
        store[k] = {"name": "-".join(names), "match": [wild(x) for x in names], "style": "fourier", "terms": []}
        carry.add(k)
    terms = store[k]["terms"]
    n = int(abs(round(pn)))
    terms[:] = [t for t in terms if t[1] != n] + [(pk / idivf, n, phase)]
    return pn < 0


def add_improper(names, rest):
    pk, phase, pn = float(rest[0]), float(rest[1]), float(rest[2])
    d = round(math.cos(math.radians(phase)))
    impropers[tuple(names)] = {"name": "-".join(names), "match": [wild(x) for x in names], "style": "cvff",
                               "params": [pk, d, int(abs(round(pn)))]}


def read_dat(path):
    L = open(path, errors="replace").read().splitlines()
    i = 1   # title
    while i < len(L) and L[i].strip():   # atom types
        w = L[i].split()
        if len(w) >= 2 and re.match(r"^[\d.]+$", w[1]):
            types[w[0]] = {"name": w[0], "element": element_of(float(w[1])), "mass": float(w[1])}
        i += 1
    i += 2   # the hydrophilic list
    dih_seen = set()
    for section in ("bonds", "angles", "dihedrals", "impropers"):
        n = {"bonds": 2, "angles": 3, "dihedrals": 4, "impropers": 4}[section]
        more, last = False, None
        while i < len(L) and L[i].strip():
            names, rest = keyfields(L[i], n)
            if section == "dihedrals" and more and not any(names):   # a further term of the last key (GLYCAM leaves the names blank)
                names = last
            last = names
            if section == "bonds":
                k = tuple(sorted(names))
                bonds[k] = {"name": "-".join(names), "match": [wild(x) for x in names], "style": "harmonic", "params": [float(rest[0]), float(rest[1])]}
            elif section == "angles":
                k = tuple(names) if names[0] <= names[2] else tuple(names[::-1])
                angles[k] = {"name": "-".join(names), "match": [wild(x) for x in names], "style": "harmonic", "params": [float(rest[0]), float(rest[1])]}
            elif section == "dihedrals":
                more = add_dihedral(names, rest, dihedrals, dih_seen)
            else:
                add_improper(names, rest)
            i += 1
        i += 1
    # 10-12 hydrogen bonds (unused in AMBER's current force fields), then the nonbonded equivalences
    while i < len(L) and L[i].strip():
        i += 1
    i += 1
    while i < len(L) and L[i].strip():
        w = L[i].split()
        for t in w[1:]:
            equiv[t] = w[0]
        i += 1
    i += 1
    while i < len(L) and not L[i].strip().startswith("END"):   # MOD4 RE, then type R* eps
        w = L[i].split()
        if len(w) >= 3 and re.match(r"^-?[\d.]+$", w[1]) and re.match(r"^-?[\d.]+$", w[2]):
            pairs[w[0]] = {"name": w[0], "match": [w[0]], "params": [float(w[2]), 2 * float(w[1]) / SIG]}
        i += 1


def read_frcmod(path):
    sec, dih_seen = None, set()
    for line in open(path, errors="replace").read().splitlines()[1:]:
        head = line.strip().upper()[:4]
        if head in ("MASS", "BOND", "ANGL", "DIHE", "IMPR", "NONB", "HBON", "CMAP"):
            sec = head
            continue
        if not line.strip():
            continue
        if sec == "MASS":
            w = line.split()
            types[w[0]] = {"name": w[0], "element": element_of(float(w[1])), "mass": float(w[1])}
        elif sec == "BOND":
            names, rest = keyfields(line, 2)
            bonds[tuple(sorted(names))] = {"name": "-".join(names), "match": [wild(x) for x in names], "style": "harmonic", "params": [float(rest[0]), float(rest[1])]}
        elif sec == "ANGL":
            names, rest = keyfields(line, 3)
            k = tuple(names) if names[0] <= names[2] else tuple(names[::-1])
            angles[k] = {"name": "-".join(names), "match": [wild(x) for x in names], "style": "harmonic", "params": [float(rest[0]), float(rest[1])]}
        elif sec == "DIHE":
            names, rest = keyfields(line, 4)
            add_dihedral(names, rest, dihedrals, dih_seen)
        elif sec == "IMPR":
            names, rest = keyfields(line, 4)
            add_improper(names, rest)
        elif sec == "NONB":
            w = line.split()
            pairs[w[0]] = {"name": w[0], "match": [w[0]], "params": [float(w[2]), 2 * float(w[1]) / SIG]}


for i, a in enumerate(sys.argv):   # in the order given, each overriding the last (leap: parm99, frcmod.ff03, then GLYCAM)
    if a == "--dat":
        read_dat(sys.argv[i + 1])
    elif a == "--frcmod":
        read_frcmod(sys.argv[i + 1])
for t, base in equiv.items():   # nonbonded equivalences: the type takes its base type's Lennard-Jones
    if t not in pairs and base in pairs:
        pairs[t] = {"name": t, "match": [t], "params": list(pairs[base]["params"]), "comment": f"Lennard-Jones of {base} (equivalence)"}
# 1-4 scaling: AMBER's 0.5 / 1/1.2 unless the source sets its own (GLYCAM: SCEE = SCNB = 1, --scale14 1 1)
LJ14, COUL14 = (1 / float(sys.argv[sys.argv.index("--scale14") + 1]), 1 / float(sys.argv[sys.argv.index("--scale14") + 2])) \
    if "--scale14" in sys.argv else (0.5, 1 / 1.2)
for d in dihedrals.values():
    t = sorted(d.pop("terms"), key=lambda x: x[1])
    d["params"] = [len(t)] + [v for term in t for v in term]
missing = sorted(t for t in types if t not in pairs)
if missing:
    notes.append(f"types without Lennard-Jones in the source: {', '.join(missing)}")
doc = {"format": "caps-forcefield", "format_version": 1, "name": arg("--name"), "source": ", ".join(p.split("/")[-1] for p in many("--dat") + many("--frcmod")),
       "references": many("--reference"), "units": "real",
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "cvff", "special": "amber"},
       "mixing": "arithmetic", "special_lj": [0, 0, LJ14], "special_coul": [0, 0, COUL14], "cutoff": 10,
       "improper_order": "center3_sorted", "equivalence": "replace", "improper_matched_order": True, "improper_reversible": True,
       "improper_max_neighbours": 3,
       "atom_types": [types[t] for t in sorted(types)], "pairs": list(pairs.values()), "bonds": list(bonds.values()),
       "angles": list(angles.values()), "dihedrals": list(dihedrals.values()), "impropers": list(impropers.values()),
       "notes": [f"converted by CAPS (bench/ff/convert_amber_dat.py) from {', '.join(p.split('/')[-1] for p in many('--dat') + many('--frcmod'))} "
                 "(AmberTools dat/leap/parm, public domain)"] + notes,
       "version": arg("--version", "")}
if arg("--typing"):   # CAPS's typing rules for this library (relative to the library file)
    doc["typing"] = arg("--typing")
json.dump(doc, open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"{OUT}: {len(types)} types, {len(pairs)} LJ, {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} dihedrals, {len(impropers)} impropers")
