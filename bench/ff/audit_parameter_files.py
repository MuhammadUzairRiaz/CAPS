#!/usr/bin/env python3
"""Consistency audit of every parameter file CAPS ships: data/forcefields/*.json, data/typing/**/*.typing.json,
data/martini/*.json and data/potentials/catalogue.json.

For each force field:
  - every type name a pair, bond, angle, dihedral, improper, automatic or bond-increment entry names is defined: it
    matches (as the loader's glob, core/src/ffdef.cpp glob_match) an atom type or an equivalence name of the file, or
    of the file it extends
  - no two entries of a term with the same style and the same type sequence (either direction for bonds, angles and
    dihedrals) carry different parameters
  - every atom type's mass is its element's standard atomic weight (IUPAC 2021, abridged), a united-atom group's
    (element + n H), or a declared shell / massless site
  - parameters in the physical range of the file's units (real: kcal/mol, A, degrees) for their style
  - the parameter count of each entry is one the loader reads for that style
  - every typing rule names a type its force field defines; every catalogue entry's file exists
  - Martini 3 molecules: every bead type defined, each molecule at an integer charge
  - potentials: each catalogued file exists and names the elements the catalogue lists

usage: audit_parameter_files.py [--ff ID] [--verbose]
Prints one line per finding: error (wrong data), dup (one entry given twice with different values), range (outside
the physical range for its style; each read against its source), unref (a term naming a type the file does not define:
the source file's own unused entries, kept as converted), warn. Exit status 1 when an error is found.
"""
import glob, gzip, json, math, os, re, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FFDIR = os.path.join(ROOT, "data", "forcefields")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
ONLY = arg("--ff", "")
VERBOSE = "--verbose" in sys.argv

# IUPAC 2021 standard atomic weights (abridged / conventional values), and the radioactive elements' usual mass numbers
MASS = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "Ar": 39.95, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996,
    "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630,
    "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41,
    "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33,
    "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145.0, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25,
    "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49,
    "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59,
    "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209.0, "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0,
    "Ac": 227.0, "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237.0, "Pu": 244.0, "Am": 243.0, "D": 2.014,
}

findings = []   # (level, ff, message)
# masses that differ from the standard atomic weight in the source file itself, kept as the source gives them
SOURCE_MASS = {("gaff2-2.11-ambertools", "Cu"): "gaff211.dat (AmberTools) gives Cu 63.010"}


def note(level, ff, msg):
    findings.append((level, ff, msg))


def glob_match(pattern, text):
    """core/src/ffdef.cpp glob_match: '*' any run, '?' one character, '\\' makes the next character literal"""
    rx = ""
    i = 0
    while i < len(pattern):
        c = pattern[i]
        if c == "\\" and i + 1 < len(pattern):
            rx += re.escape(pattern[i + 1]); i += 2; continue
        rx += ".*" if c == "*" else "." if c == "?" else re.escape(c)
        i += 1
    return re.fullmatch(rx, text, re.S) is not None


def load(path):
    """the file and the file it builds on: its "extends", or the library an overlay names in its notes
    ("optional overlay on oplsaa2024.json", applied with caps ff apply --overlay)"""
    d = json.load(open(path))
    base = None
    if "extends" in d:
        base = load(os.path.join(os.path.dirname(path), d["extends"]))
    else:
        m = re.search(r"overlay on ([\w.-]+\.json)", json.dumps(d.get("notes", "")))
        if m and os.path.exists(os.path.join(os.path.dirname(path), m.group(1))):
            base = load(os.path.join(os.path.dirname(path), m.group(1)))
    return d, base


STYLE_NPARAMS = {   # parameter counts the loader reads for each (term, style); None: any count of that form
    ("bonds", "harmonic"): {2}, ("bonds", "class2"): {4}, ("bonds", "morse"): {3}, ("bonds", "fene"): {2, 4},
    ("bonds", "harmonic/order"): {2}, ("bonds", "gromos"): {2},
    ("angles", "harmonic"): {2}, ("angles", "class2"): {4}, ("angles", "charmm"): {4}, ("angles", "cosine"): {1},
    ("angles", "cosine/squared"): {2}, ("angles", "sdk"): {2}, ("angles", "screened"): {4},
    ("dihedrals", "opls"): {4}, ("dihedrals", "class2"): {6}, ("dihedrals", "harmonic"): {3},
    ("dihedrals", "charmm"): {4}, ("dihedrals", "charmmfsw"): {4}, ("dihedrals", "cistrans"): {2},
    ("impropers", "cvff"): {3}, ("impropers", "class2"): {2}, ("impropers", "harmonic"): {2},
    ("impropers", "fourier"): {4, 5}, ("impropers", "planar"): {2}, ("impropers", "umbrella"): {2},
}


def term_style(d, base, term, e):
    st = (d.get("styles") or (base[0].get("styles") if base else None) or {})
    key = {"pairs": "pair", "bonds": "bond", "angles": "angle", "dihedrals": "dihedral", "impropers": "improper",
           "auto_bonds": "bond", "auto_angles": "angle", "auto_dihedrals": "dihedral"}.get(term)
    return e.get("style") or st.get(key, "")


def check_ff(path):
    fid = os.path.basename(path)[:-5]
    d, base = load(path)
    units = d.get("units") or (base[0].get("units") if base else None)
    types = list(d.get("atom_types", []))
    if base:
        types += base[0].get("atom_types", [])
        if base[1]:
            types += base[1][0].get("atom_types", [])
    names = set()
    for t in types:
        names.add(t["name"])
        for v in (t.get("equivalence") or {}).values():
            if isinstance(v, str):
                names.add(v)
    # names a term may use besides atom types: bead names, virtual sites, shells
    for v in d.get("virtual_sites", []) or []:
        names.add(v.get("type", ""))

    known = {}

    def defined(p):
        if p not in known:
            if not any(c in p for c in "*?\\"):
                known[p] = p in names
            else:
                known[p] = any(glob_match(p, n) for n in names)
        return known[p]

    # 1. masses
    seen_type = {}
    for t in d.get("atom_types", []):
        n, el, m = t["name"], t.get("element", ""), t.get("mass")
        if n in seen_type:
            same = {k: v for k, v in t.items() if k != "description"} == {k: v for k, v in seen_type[n].items() if k != "description"}
            note("error" if not same else "warn", fid, f"atom type {n} defined twice" + ("" if same else " with different values"))
        seen_type[n] = t
        if m is None or not el:
            continue   # coarse-grained beads and unnamed sites
        desc = (t.get("description", "") + " " + n).lower()
        if el == "X" and m == 0:
            continue   # a massless site (GROMACS's MW / LP / dummy types)
        if "relative free energy" in desc:
            continue   # the .frc files' null atom (nu)
        if (fid, n) in SOURCE_MASS:
            note("warn", fid, f"type {n} mass {m}: {SOURCE_MASS[fid, n]}, kept")
            continue
        if el not in MASS:
            if el in ("X", "L", "Lp", "D") and any(w in desc for w in ("dummy", "lone pair", "electron pair", "null atom", "site", " m ", "relative free energy")):
                continue   # a declared pseudo-atom (BOSS / GROMACS dummies, .frc lone pairs and null atoms)
            note("error", fid, f"type {n}: element '{el}' is not an element (mass {m}, '{t.get('description', '')[:50]}')")
            continue
        if any(w in desc for w in ("lone pair", "with lone pairs", "electron pair")):
            continue   # an O whose mass is shared with its lone-pair sites (IFF's TIP5P-like O, 14.9994 + 2 x 0.5)
        ref = MASS[el]
        lp = [o.get("mass", 0) for o in d.get("atom_types", []) if o.get("element") in ("Lp", "L") and 0 < o.get("mass", 0) < 1]
        if el == "O" and any(abs(m + 2 * x - ref) <= 0.002 for x in lp):
            continue   # IFF's TIP5P-like O: 14.9994 + two lone-pair sites of 0.5
        if m == 0 or m < 0.5:
            if not any(w in desc for w in ("shell", "virtual", "massless", "dummy", "site", "lone pair", "drude", "lp")):
                note("warn", fid, f"type {n} ({el}) mass {m}: zero or tiny without a stated shell / site")
            continue
        if abs(m - ref) <= max(0.02, 0.0015 * ref):
            continue
        # united-atom CHn / NHn / OH (TraPPE, OPLS-UA, ...)
        nh = round((m - ref) / 1.008)
        if 1 <= nh <= 4 and abs(m - ref - nh * 1.008) <= 0.03:
            continue
        # a core-shell split: a core and a shell of the same element whose masses add up to the atom's
        if m < ref and any(o is not t and o.get("element") == el and abs(o.get("mass", 0) + m - ref) <= max(0.02, 0.0015 * ref)
                           for o in d.get("atom_types", [])):
            continue
        desc = (t.get("description", "") + " " + n).lower()
        if "shell" in desc and m < 0.15 * ref:
            continue   # a light shell (the paper's own shell mass) whose core keeps the element's mass
        lvl = "error" if abs(m - ref) > 0.5 else "warn"
        note(lvl, fid, f"type {n} ({el}) mass {m}, standard atomic weight {ref}")

    # 2. undefined types, duplicates, parameter counts and ranges
    for term in ("pairs", "bonds", "angles", "dihedrals", "impropers", "auto_bonds", "auto_angles", "auto_dihedrals",
                 "bond_increments"):
        groups = defaultdict(list)
        for e in d.get(term, []) or []:
            match = e.get("match", [])
            for p in match:
                if p in ("*", "X", "x", ""):
                    continue
                if not defined(p):
                    note("unref", fid, f"{term} {e.get('name')}: type '{p}' is not defined")
            st = term_style(d, base, term, e)
            pr = e.get("params", [])
            if not all(isinstance(x, (int, float)) and math.isfinite(x) for x in pr):
                note("error", fid, f"{term} {e.get('name')}: a parameter is not a finite number: {pr}")
                continue
            want = STYLE_NPARAMS.get((term.replace("auto_", ""), st))
            if want and len(pr) not in want and not (st == "fourier"):
                note("error", fid, f"{term} {e.get('name')}: style {st} with {len(pr)} parameters (the loader reads {sorted(want)})")
            if term.replace("auto_", "") == "dihedrals" and st == "fourier" and (len(pr) - 1) % 3 != 0 and len(pr) != 4:
                note("error", fid, f"{term} {e.get('name')}: fourier with {len(pr)} parameters (m, then K n d per term)")
            if term in ("bonds", "angles", "dihedrals"):
                key = (st, min(tuple(match), tuple(reversed(match))))
            else:
                key = (st, tuple(match))
            groups[key].append(e)
            if units == "real":
                check_range(fid, term.replace("auto_", ""), st, e, d)
        for key, es in groups.items():
            ps = {json.dumps([round(x, 6) if isinstance(x, float) else x for x in e["params"]]) for e in es}
            if len(ps) > 1 and term != "bond_increments":
                contexts = {json.dumps({k: v for k, v in e.items() if k not in ("params", "comment", "name", "fit")}, sort_keys=True) for e in es}
                if len(contexts) == len(es):
                    continue   # the entries differ by context (contact, molecule, ...)
                note("dup", fid, f"{term} {'-'.join(key[1])} ({key[0]}): {len(es)} entries with different parameters: " +
                     "; ".join(f"{e['params']}" + (f" [{e.get('comment','')[:60]}]" if e.get('comment') else "") for e in es[:4]))
            elif len(es) > 1 and VERBOSE:
                note("info", fid, f"{term} {'-'.join(key[1])}: {len(es)} identical entries")
    return d


CG = re.compile(r"martini|sdk|cooke|mw-molinero")
METALS = set("Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Cs Ba La Ce Hf Ta W "
             "Re Os Ir Pt Au Hg Tl Pb Bi U".split())


def check_range(fid, term, st, e, d):
    p = e["params"]
    nm = f"{term} {e.get('name')}"
    bad = lambda msg: note("range", fid, f"{nm} {st}: {msg} {p}")
    el = {t["name"]: t.get("element", "") for t in d.get("atom_types", [])}
    if CG.search(fid) and term in ("bonds", "angles"):
        return   # coarse-grained beads: bond lengths of 3-6 A and soft constants are the model
    if term == "bonds" and st == "class2" and len(p) >= 2 and p[0] == 0.0 and p[1] > 0:
        return   # a core-shell spring (r0 0, k2 the paper's spring constant)
    if term == "bonds" and len(p) >= 2 and p[1] < 0.5 and any(x in ("LP", "DM", "EP", "LPH") or x.startswith("LP") for x in e["match"]):
        return   # a massless site held near its atom (TIP4F's M, BOSS dummies)
    if term == "pairs" and st.startswith("lj") and len(p) >= 2 and p[0] < 1e-20:
        return   # a repulsion-only pair written through the epsilon-sigma form (A/r^9 with epsilon -> 0)
    if term == "pairs" and len(p) >= 2 and any(el.get(x, "") in METALS for x in e["match"]) and 0 <= p[0] <= 40:
        return   # metal Lennard-Jones well depths (Heinz's metals, the .frc files' W, Mo, Pt) reach 20-30 kcal/mol
    if term == "pairs":
        if st.startswith("lj") or st in ("",) or st.startswith("lj/"):
            if len(p) >= 2:
                eps, sig = p[0], p[1]
                if eps < 0 or eps > 5:
                    bad("epsilon outside 0-5 kcal/mol")
                if eps > 0 and not (0.3 <= sig <= 8.0):
                    bad("sigma outside 0.3-8 A")
        elif st.startswith("buck") and len(p) >= 3:
            A, rho, C = p[0], p[1], p[2]
            if A < 0 or not (0.02 <= rho <= 1.5) or C < 0:
                bad("A < 0, rho outside 0.02-1.5 A or C < 0")
    elif term == "bonds":
        if st in ("harmonic", "harmonic/order") and len(p) >= 2:
            if not (0 <= p[0] <= 5000) or not (0.5 <= p[1] <= 4.0):
                bad("K outside 0-5000 kcal/mol/A^2 or r0 outside 0.5-4 A")
        elif st == "class2" and len(p) >= 4:
            if not (0.5 <= p[0] <= 4.0) or p[1] <= 0:
                bad("r0 outside 0.5-4 A or K2 <= 0")
        elif st == "morse" and len(p) >= 3:
            if p[0] <= 0 or not (0.3 <= p[1] <= 10) or not (0.5 <= p[2] <= 4):
                bad("D <= 0, alpha outside 0.3-10 or r0 outside 0.5-4")
    elif term == "angles":
        if st in ("harmonic", "charmm", "sdk", "cosine/squared", "screened") and len(p) >= 2:
            if not (0 <= p[0] <= 1500) or not (0 < p[1] <= 180):
                bad("K outside 0-1500 or theta0 outside 0-180")
        elif st == "class2" and len(p) >= 4:
            if not (0 < p[0] <= 180) or p[1] < 0:
                bad("theta0 outside 0-180 or K2 < 0")
    elif term == "dihedrals":
        if st == "fourier":
            m = int(p[0])
            if m != (len(p) - 1) // 3:
                bad("term count m disagrees with the parameters")
            for k in range(m):
                K, n, dd = p[1 + 3 * k: 4 + 3 * k]
                if abs(K) > 60 or n != int(n) or not (0 <= n <= 6) or not (-360 <= dd <= 360):
                    bad(f"term {k}: |K| > 60, n not 0-6 or delta outside -360..360")
        elif st == "opls":
            if any(abs(x) > 60 for x in p):
                bad("|V| > 60 kcal/mol")
        elif st in ("charmm", "charmmfsw") and len(p) >= 3:
            if abs(p[0]) > 60 or p[1] != int(p[1]) or not (0 <= p[1] <= 6):
                bad("|K| > 60 or n not 0-6")
        elif st == "harmonic" and len(p) >= 3:
            if abs(p[0]) > 60 or p[1] not in (-1, 1) or not (0 <= p[2] <= 6):
                bad("|K| > 60, d not +-1 or n not 0-6")
        elif st == "class2" and len(p) >= 6:
            if any(abs(p[i]) > 60 for i in (0, 2, 4)) or any(not (-360 <= p[i] <= 360) for i in (1, 3, 5)):
                bad("|V| > 60 or phase outside -360..360")
    elif term == "impropers":
        if st == "cvff" and len(p) >= 3:
            if p[1] not in (-1, 1) or not (0 <= p[2] <= 6) or abs(p[0]) > 100:
                bad("d not +-1, n not 0-6 or |K| > 100")
        elif st == "harmonic" and len(p) >= 2:
            if not (0 <= p[0] <= 500) or not (-180 <= p[1] <= 180):
                bad("K outside 0-500 or chi0 outside -180..180")


def check_typing():
    cat = json.load(open(os.path.join(FFDIR, "catalogue.json")))["forcefields"]
    for e in cat:
        f = os.path.join(FFDIR, e["file"])
        if not os.path.exists(f):
            local = "never distributed" in e.get("license", "") or "Not distributed" in e.get("notes", "")
            note("warn" if local else "error", e["id"], f"catalogue file {e['file']} missing" + (" (built locally, never distributed)" if local else ""))
            continue
        if ONLY and e["id"] != ONLY:
            continue
        rules = (e.get("typing") or {}).get("rules")
        if not rules:
            continue
        rp = os.path.join(os.path.dirname(FFDIR), rules)
        if not os.path.exists(rp):
            note("error", e["id"], f"typing rules {rules} missing")
            continue
        d, base = load(f)
        defined = {t["name"] for t in d.get("atom_types", [])}
        if base:
            defined |= {t["name"] for t in base[0].get("atom_types", [])}
            if base[1]:
                defined |= {t["name"] for t in base[1][0].get("atom_types", [])}
        tr = json.load(open(rp))
        unknown = set(tr.get("unknown_types", []) or [])
        missing = defaultdict(int)
        for r in tr.get("rules", []):
            t = r["type"]
            if t not in defined and t not in unknown:
                missing[t] += 1
        for t, n in sorted(missing.items()):
            note("error", e["id"], f"typing rules ({rules}) assign '{t}' ({n} rules), not an atom type of {e['file']}")
        c = e.get("counts") or {}
        for k in ("atom_types", "pairs", "bonds", "angles", "dihedrals", "impropers"):
            if k == "pairs" and d.get("pair_table"):
                continue   # pairs in a separate table (Martini 3's nonbonded TSV)
            if k in c and c[k] != len(d.get(k, []) or []):
                note("warn", e["id"], f"catalogue counts {k} {c[k]}, the file has {len(d.get(k, []) or [])}")


def check_martini():
    d = json.load(open(os.path.join(ROOT, "data", "martini", "martini3-molecules.json")))
    beads = set()
    with gzip.open(os.path.join(ROOT, "data", "martini", "martini3-nonbonded.tsv.gz"), "rt") as f:
        for line in f:
            w = line.split("\t")
            if len(w) >= 2 and not line.startswith("#"):
                beads.add(w[0]); beads.add(w[1])
    ff = json.load(open(os.path.join(FFDIR, "martini3.json")))
    beads |= {t["name"] for t in ff["atom_types"]}
    for name, m in d["molecules"].items():
        q = sum(a.get("charge", 0.0) or 0.0 for a in m["atoms"])
        if abs(q - round(q)) > 1e-3:
            note("error", "martini3-molecules", f"{name}: charges sum to {q:.4f}, not an integer")
        for a in m["atoms"]:
            if a["type"] not in beads:
                note("error", "martini3-molecules", f"{name}: bead type {a['type']} not defined")
                break
        n = len(m["atoms"])
        for k in ("bonds", "constraints", "angles", "dihedrals"):
            for b in m.get(k, []):
                idx = b.get("atoms", b) if isinstance(b, dict) else b
                ii = [x for x in idx if isinstance(x, int)][: {"bonds": 2, "constraints": 2, "angles": 3, "dihedrals": 4}[k]]
                if any(x < 0 or x >= n for x in ii):
                    note("error", "martini3-molecules", f"{name}: {k} entry {b} indexes past its {n} atoms")
                    break


def check_potentials():
    pd = os.path.join(ROOT, "data", "potentials")
    cat = json.load(open(os.path.join(pd, "catalogue.json")))["potentials"]
    for p in cat:
        f = os.path.join(pd, p["file"])
        if not os.path.exists(f):
            note("error", "potentials", f"{p['id']}: {p['file']} missing")
            continue
        txt = open(f, errors="replace").read()
        body = "\n".join(l for l in txt.splitlines() if not l.lstrip().startswith("#"))
        for el in p.get("elements", []):
            if not re.search(r"(?<![A-Za-z])" + re.escape(el) + r"(?![a-z])", body):
                note("warn", "potentials", f"{p['id']}: element {el} not found in {p['file']}")
        m = re.search(r"UNITS:\s*(\w+)", txt)
        if m and p.get("units") and m.group(1) != p["units"]:
            note("error", "potentials", f"{p['id']}: file says units {m.group(1)}, the catalogue {p['units']}")


def main():
    for f in sorted(glob.glob(os.path.join(FFDIR, "*.json"))):
        if f.endswith("catalogue.json"):
            continue
        if ONLY and os.path.basename(f)[:-5] != ONLY:
            continue
        check_ff(f)
    check_typing()
    if not ONLY:
        check_martini()
        check_potentials()
    order = {"error": 0, "dup": 1, "range": 2, "unref": 3, "warn": 4, "info": 5}
    for lvl, ff, msg in sorted(findings, key=lambda x: (order[x[0]], x[1])):
        print(f"{lvl:5s} {ff}: {msg}")
    counts = defaultdict(int)
    for lvl, _, _ in findings:
        counts[lvl] += 1
    print("summary:", dict(counts))
    return 1 if counts["error"] else 0


if __name__ == "__main__":
    sys.exit(main())
