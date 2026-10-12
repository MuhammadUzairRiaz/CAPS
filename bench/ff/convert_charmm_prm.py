#!/usr/bin/env python3
"""A CAPS force-field file from CHARMM's own parameter files (the MacKerell lab's toppar distribution).

  convert_charmm_prm.py OUT.json --name NAME --prm par_all36_prot.prm [--prm toppar_water_ions.str ...]
                        [--rtf top_all36_prot.rtf ...] [--reference "..."] [--source "..."] [--version "..."]

Reads the CHARMM36-format sections of each file in order (ATOMS / MASS, BONDS, ANGLES, DIHEDRALS, IMPROPER, NONBONDED,
NBFIX; a stream file's read para block too) and writes CAPS's format:

  pairs      epsilon = |eps|, sigma = 2 Rmin/2 / 2^(1/6); the 1-4 values (CHARMM's optional second pair) as lj14
  NBFIX      pair-specific entries, sigma = Rmin / 2^(1/6)
  bonds      harmonic K (CHARMM's K, energy K (r - r0)^2) and r0
  angles     harmonic, or charmm with the Urey-Bradley K and S0
  dihedrals  fourier: the number of terms, then K n delta for each multiplicity CHARMM lists for the same four types
  impropers  harmonic K psi0 (multiplicity 0), else fourier 1 K n delta (CHARMM's cosine improper)
  X          wildcard, written * (CAPS matches by glob; the name keeps X)

Masses from MASS records (the rtf, or the prm's ATOMS section); the element is the one whose standard mass is nearest.
CMAP is not represented (noted in the file). Later files override earlier entries with the same types, as CHARMM reads
them; a duplicate in one file (CHARMM keeps the last) is noted.
"""
import json, math, re, sys

arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
many = lambda k: [sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == k]
OUT = sys.argv[1]
SIG = 2 ** (1 / 6)
ELEMENTS = [("H", 1.008), ("He", 4.0026), ("Li", 6.94), ("B", 10.81), ("C", 12.011), ("N", 14.007), ("O", 15.999), ("F", 18.998),
            ("Ne", 20.18), ("Na", 22.99), ("Mg", 24.305), ("Al", 26.982), ("Si", 28.085), ("P", 30.974), ("S", 32.06),
            ("Cl", 35.45), ("K", 39.098), ("Ca", 40.078), ("Fe", 55.845), ("Cu", 63.546), ("Zn", 65.38), ("Br", 79.904),
            ("Rb", 85.468), ("Cs", 132.905), ("I", 126.904), ("Ba", 137.327), ("Cd", 112.414), ("Li", 6.941), ("Ar", 39.948),
            ("Se", 78.971), ("Mn", 54.938), ("Co", 58.933), ("Ni", 58.693), ("Hg", 200.59), ("Sr", 87.62), ("Pt", 195.08),
            ("Au", 196.967), ("Ag", 107.868), ("Pb", 207.2), ("Ra", 226.0), ("Fr", 223.0), ("Tl", 204.38), ("Be", 9.012)]


def element_of(name, mass, given=None):
    if given:
        return given.capitalize()
    if mass < 0.5:
        return "X"   # lone pair / dummy site
    return min(ELEMENTS, key=lambda e: abs(e[1] - mass))[0]


def lines(path):
    """logical lines without comments, '-' continuations joined, with the section they belong to"""
    buf = ""
    for raw in open(path, errors="replace"):
        s = raw.split("!")[0].rstrip()
        if s.endswith("-") and not re.match(r"^\s*-?\d", s.split()[-1] if s.split() else ""):
            buf += s[:-1] + " "
            continue
        s = buf + s
        buf = ""
        if s.strip():
            yield s.strip()


# CHARMM reads a keyword by its first four letters (BOND, ANGL / THET, DIHE / PHI, IMPR / IMPH, NONB / NBON ...)
SECTIONS = {"ATOM": "atoms", "BOND": "bonds", "ANGL": "angles", "THET": "angles", "DIHE": "dihedrals", "PHI": "dihedrals",
            "IMPR": "impropers", "IMPH": "impropers", "CMAP": "cmap", "NONB": "nonbonded", "NBON": "nonbonded", "NBFI": "nbfix",
            "HBON": "hbond", "END": None}


def section(w):
    """the section a header line opens, or None for a data line"""
    h = w[0].upper()
    k = h[:4] if h[:4] in SECTIONS else (h if h in SECTIONS else None)
    if k is None:
        return False, None
    # a header stands alone, except NONBONDED / NBFIX / CMAP / HBOND, which carry options
    if len(w) == 1 or k in ("NONB", "NBON", "CMAP", "HBON", "NBFI", "END"):
        if k in ("NONB", "NBON") and len(w) > 1 and re.match(r"^-?[\d.]+$", w[1]):
            return False, None   # a type named NONB... with numbers: data
        return True, SECTIONS[k]
    return False, None


types, pairs, nbfix, bonds, angles, dihedrals, impropers, notes = {}, {}, {}, {}, {}, {}, {}, []
cmap_seen = False


def wild(t):
    return "*" if t.upper() == "X" else t.replace("*", "\\*").replace("?", "\\?")   # a literal * or ? escaped (CAPS patterns are globs)


def read_masses(path):
    for s in lines(path):
        w = s.split()
        if w[0].upper() == "MASS" and len(w) >= 4:
            el = w[4] if len(w) >= 5 and re.fullmatch(r"[A-Za-z]{1,2}", w[4]) else None
            m = float(w[3])
            types[w[2]] = {"name": w[2], "element": element_of(w[2], m, el), "mass": m}


def read_params(path):
    global cmap_seen
    sec, in_stream_para = None, not path.endswith(".str")
    seen = {k: set() for k in ("bonds", "angles", "dihedral_keys", "impropers")}
    dih_file = {}
    for s in lines(path):
        w = s.split()
        head = w[0].upper()
        if path.endswith(".str") and head == "READ":
            in_stream_para = len(w) > 1 and w[1].upper().startswith("PARA")
            continue
        if not in_stream_para:
            if head == "MASS":
                el = w[4] if len(w) >= 5 and re.fullmatch(r"[A-Za-z]{1,2}", w[4]) else None
                types.setdefault(w[2], {"name": w[2], "element": element_of(w[2], float(w[3]), el), "mass": float(w[3])})
            continue
        is_header, new = section(w)
        if is_header:
            sec = new
            if sec == "cmap":
                cmap_seen = True
            continue
        if head == "END" or head == "RETURN":
            sec = None
            continue
        if sec == "atoms" and head == "MASS":
            el = w[4] if len(w) >= 5 and re.fullmatch(r"[A-Za-z]{1,2}", w[4]) else None
            types.setdefault(w[2], {"name": w[2], "element": element_of(w[2], float(w[3]), el), "mass": float(w[3])})
        elif sec == "bonds" and len(w) >= 4:
            key = tuple(sorted(w[:2]))
            if key in seen["bonds"]:
                notes.append(f"{path.split('/')[-1]}: bond {'-'.join(w[:2])} listed twice; the last is used, as CHARMM does")
            seen["bonds"].add(key)
            bonds[key] = {"name": "-".join(w[:2]), "match": [wild(x) for x in w[:2]], "style": "harmonic", "params": [float(w[2]), float(w[3])]}
        elif sec == "angles" and len(w) >= 5:
            key = tuple(w[:3]) if w[0] <= w[2] else tuple(w[:3][::-1])
            p = [float(w[3]), float(w[4])]
            e = {"name": "-".join(w[:3]), "match": [wild(x) for x in w[:3]], "style": "harmonic", "params": p}
            if len(w) >= 7 and float(w[5]) != 0:
                e["style"], e["params"] = "charmm", p + [float(w[5]), float(w[6])]
            angles[key] = e
        elif sec == "dihedrals" and len(w) >= 7:
            k4 = tuple(w[:4]) if tuple(w[:4]) <= tuple(w[:4][::-1]) else tuple(w[:4][::-1])
            K, n, d = float(w[4]), int(float(w[5])), float(w[6])
            if k4 not in dih_file:   # the first listing in this file replaces any earlier file's terms
                dih_file[k4] = []
                dihedrals[k4] = {"name": "-".join(w[:4]), "match": [wild(x) for x in w[:4]], "style": "fourier", "terms": []}
            terms = dihedrals[k4]["terms"]
            terms[:] = [t for t in terms if t[1] != n] + [(K, n, d)]   # the same multiplicity again: the last wins
        elif sec == "impropers" and len(w) >= 7:
            key = tuple(w[:4])
            K, n, d = float(w[4]), int(float(w[5])), float(w[6])
            # multiplicity 0: harmonic K (psi - psi0)^2; otherwise CHARMM's cosine improper K [1 + cos(n psi - d)]
            impropers[key] = {"name": "-".join(w[:4]), "match": [wild(x) for x in w[:4]],
                              **({"style": "harmonic", "params": [K, d]} if n == 0 else {"style": "fourier", "params": [1, K, n, d]})}
        elif sec == "nonbonded" and len(w) >= 4 and not head.startswith(("CUTNB", "NBXMOD", "!")) and re.match(r"^-?[\d.]+$", w[1]):
            eps, rmin2 = abs(float(w[2])), float(w[3])
            e = {"name": w[0], "match": [w[0]], "params": [eps, 2 * rmin2 / SIG]}
            if len(w) >= 7:
                e["lj14"] = [abs(float(w[5])), 2 * float(w[6]) / SIG]
            pairs[w[0]] = e
        elif sec == "nbfix" and len(w) >= 4:
            key = tuple(sorted(w[:2]))
            e = {"name": "-".join(w[:2]), "match": list(w[:2]), "params": [abs(float(w[2])), float(w[3]) / SIG]}
            if len(w) >= 6:
                e["lj14"] = [abs(float(w[4])), float(w[5]) / SIG]
            nbfix[key] = e


for r in many("--rtf"):
    read_masses(r)
for p in many("--prm"):
    read_params(p)

for d in dihedrals.values():
    t = sorted(d.pop("terms"), key=lambda x: x[1])
    d["params"] = [len(t)] + [v for term in t for v in term]
# entries naming a type no file read defines (a stream's cross terms with another CHARMM file: na_combined, pyridines)
# can never apply in this library: left out, counted
foreign = set()
for coll in (bonds, angles, dihedrals, impropers, nbfix):
    for k in list(coll):
        out = [x for x in coll[k]["match"] if x != "*" and x not in types]
        if out:
            foreign |= set(out)
            del coll[k]
for t in list(pairs):
    if t not in types:
        foreign.add(t)
        del pairs[t]
if foreign:
    notes.append(f"entries naming types of other CHARMM files ({', '.join(sorted(foreign))}) left out: they apply only when those files are read too")
doc = {"format": "caps-forcefield", "format_version": 1, "name": arg("--name"), "source": arg("--source", ", ".join(many("--prm"))),
       "references": many("--reference"), "units": "real",
       "styles": {"pair": "lj/charmmfsw/coul/long", "bond": "harmonic", "angle": "charmm", "dihedral": "charmmfsw", "improper": "harmonic", "special": "charmm"},
       "mixing": "arithmetic", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 12,
       "pair_settings": {"lj_inner": 10, "model_cutoff": True, "lj_modifier": "charmm-force-switch"},
       "improper_order": "center1_sorted", "equivalence": "replace", "improper_matched_order": True, "improper_reversible": True,
       "improper_all_explicit": True, "improper_max_neighbours": 3,
       "atom_types": [types[t] for t in sorted(types)], "pairs": list(pairs.values()) + list(nbfix.values()),
       "bonds": list(bonds.values()), "angles": list(angles.values()), "dihedrals": list(dihedrals.values()), "impropers": list(impropers.values()),
       "notes": [f"converted by CAPS (bench/ff/convert_charmm_prm.py) from {', '.join(p.split('/')[-1] for p in many('--prm'))}"
                 + (", masses from " + ", ".join(r.split('/')[-1] for r in many("--rtf")) if many("--rtf") else "")]
                + (["CMAP corrections in the source are not represented"] if cmap_seen else []) + notes,
       "version": arg("--version", "")}
json.dump(doc, open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"{OUT}: {len(types)} types, {len(pairs)} LJ + {len(nbfix)} NBFIX, {len(bonds)} bonds, {len(angles)} angles, "
      f"{len(dihedrals)} dihedrals, {len(impropers)} impropers; {len(notes)} notes")
