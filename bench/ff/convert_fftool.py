#!/usr/bin/env python3
"""A CAPS force-field file from a force field in fftool's format (Padua's .ff files: CL&P's il.ff, its des.ff).

  convert_fftool.py FILE.ff [FILE.ff ...] OUT.json --name NAME [--reference "..."] [--version "..."] [--typing ../typing/ID.typing.json]

The format (units kJ/mol, A, degrees): ATOMS "name class mass charge lj sigma epsilon" (class names the type in the
bonded terms), BONDS / ANGLES "i j [k] harm|cons x0 k" with energy k/2 (x - x0)^2 ('cons': a constrained bond or angle,
kept flexible here with its constant), DIHEDRALS "i j k l opls V1 V2 V3 V4", IMPROPER the same with the central atom
third. CAPS's form, real units: pairs epsilon / 4.184, sigma; bonds and angles harmonic K = k / 2 / 4.184; dihedrals opls
V / 4.184; impropers cvff K = V2 / 2 / 4.184, d = -1, n = 2 (V2 / 2 [1 - cos 2 phi]), centre third; geometric mixing,
1-4 x 0.5 (OPLS-AA, which CL&P follows). Later files add to and override earlier ones.
"""
import json, sys

files = [a for a in sys.argv[1:] if a.endswith(".ff")]
OUT = next(a for a in sys.argv[1:] if a.endswith(".json") and sys.argv[sys.argv.index(a) - 1] not in ("--typing",))
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
many = lambda k: [sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == k]
KJ = 4.184
ELEMENTS = [("H", 1.008), ("Li", 6.94), ("B", 10.81), ("C", 12.011), ("N", 14.007), ("O", 15.999), ("F", 18.998), ("Na", 22.99),
            ("Mg", 24.305), ("Al", 26.982), ("Si", 28.085), ("P", 30.974), ("S", 32.06), ("Cl", 35.45), ("K", 39.098), ("Ca", 40.078),
            ("Zn", 65.38), ("Br", 79.904), ("I", 126.904), ("Cs", 132.905)]


def element_of(m):
    return min(ELEMENTS, key=lambda e: abs(e[1] - m))[0] if m > 0.5 else "X"


types, pairs, bonds, angles, dihedrals, impropers, notes = {}, {}, {}, {}, {}, {}, []
for path in files:
    sec = None
    for raw in open(path, errors="replace"):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        w = line.split()
        if w[0] in ("ATOMS", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER", "IMPROPERS"):
            sec = w[0]
            continue
        if sec == "ATOMS" and len(w) >= 7 and w[4] == "lj":
            name, cls, m, q, sig, eps = w[0], w[1], float(w[2]), float(w[3]), float(w[5]), float(w[6])
            types[name] = {"name": name, "element": element_of(m), "mass": m, "charge": q,
                           "equivalence": {"bond": cls, "angle": cls, "dihedral": cls, "improper": cls}}
            pairs[name] = {"name": name, "match": [name], "params": [eps / KJ, sig]}
        elif sec == "BONDS" and len(w) >= 5:
            k = tuple(sorted(w[:2]))
            bonds[k] = {"name": "-".join(w[:2]), "match": w[:2], "style": "harmonic", "params": [float(w[4]) / 2 / KJ, float(w[3])],
                        **({"comment": "constrained in the source"} if w[2] == "cons" else {})}
        elif sec == "ANGLES" and len(w) >= 6:
            k = tuple(w[:3]) if w[0] <= w[2] else tuple(w[2::-1])
            angles[k] = {"name": "-".join(w[:3]), "match": w[:3], "style": "harmonic", "params": [float(w[5]) / 2 / KJ, float(w[4])],
                         **({"comment": "constrained in the source"} if w[3] == "cons" else {})}
        elif sec == "DIHEDRALS" and len(w) >= 9 and w[4] == "opls":
            k = tuple(w[:4]) if tuple(w[:4]) <= tuple(w[3::-1]) else tuple(w[3::-1])
            dihedrals[k] = {"name": "-".join(w[:4]), "match": w[:4], "style": "opls", "params": [float(x) / KJ for x in w[5:9]]}
        elif sec in ("IMPROPER", "IMPROPERS") and len(w) >= 9 and w[4] == "opls":
            impropers[tuple(w[:4])] = {"name": "-".join(w[:4]), "match": w[:4], "style": "cvff", "params": [float(w[6]) / 2 / KJ, -1, 2]}
doc = {"format": "caps-forcefield", "format_version": 1, "name": arg("--name"), "source": ", ".join(f.split("/")[-1] for f in files),
       "references": many("--reference"), "units": "real",
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"},
       "mixing": "geometric", "special_lj": [0, 0, 0.5], "special_coul": [0, 0, 0.5], "cutoff": 12,
       "improper_order": "center3_sorted", "equivalence": "replace",
       "atom_types": list(types.values()), "pairs": list(pairs.values()), "bonds": list(bonds.values()), "angles": list(angles.values()),
       "dihedrals": list(dihedrals.values()), "impropers": list(impropers.values()),
       "notes": [f"converted by CAPS (bench/ff/convert_fftool.py) from {', '.join(f.split('/')[-1] for f in files)}"] + notes,
       "version": arg("--version", "")}
if arg("--typing"):
    doc["typing"] = arg("--typing")
json.dump(doc, open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"{OUT}: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} dihedrals, {len(impropers)} impropers")
