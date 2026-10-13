#!/usr/bin/env python3
"""A CAPS force-field file from CHARMM's own parameter files (the MacKerell lab's toppar distribution).

  convert_charmm_prm.py OUT.json --name NAME --prm par_all36_prot.prm [--prm toppar_water_ions.str ...]
                        [--rtf top_all36_prot.rtf ...] [--reference "..."] [--source "..."] [--version "..."]
                        [--one-torsion-per-bond] [--typing ../typing/ID.typing.json] [--note "..." ...]
                        [--impropers-from-rtf]

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
CMAP is not represented (noted in the file).

--impropers-from-rtf: CHARMM puts an improper only where a residue's IMPR record lists one, while CAPS (which types by
SMARTS, without residues) would put one wherever a parameter entry matches. The topology's IMPR records become the
library's "improper_placement": for every three-connected centre of every residue (no links to other residues), its
environment (its type and its three neighbours' types) and the impropers the residue lists there (types in IMPR
order, centre position). Where residues disagree for one environment the most common listing is kept (counted in the
notes). Parameter entries that no IMPR record uses are listed as unlisted and never applied; at a centre whose
environment no residue shows, the other entries apply as their types match. Later files override earlier entries with the same types, as CHARMM reads
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
    # a united atom (CHARMM19's CH2E, SH1E): the element its name starts with, plus whole hydrogens
    for sym in sorted({e for e, _ in ELEMENTS}, key=len, reverse=True):
        if name.upper().startswith(sym.upper()):
            m_el = dict(ELEMENTS)[sym]
            if any(abs(mass - m_el - n * 1.008) < 0.02 for n in range(0, 5)):
                return sym
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


nb_settings = set()   # the NONBONDED header's keywords (GEOM: geometric combination; E14FAC x: 1-4 Coulomb factor)
nb_values = {}


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
        if (is_header and new == "nonbonded") or (sec == "nonbonded" and head in ("CUTNB", "NBXMOD")):
            up = [x.upper() for x in w]
            nb_settings.update(up)
            for k, x in enumerate(up[:-1]):
                if re.match(r"^-?[\d.]+$", up[k + 1]):
                    nb_values[x] = float(up[k + 1])
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
# CHARMM's wildcard Lennard-Jones entries (param19: C% = C and one more character, N* = every N type): each gives its
# parameters to the types it matches that have no entry of their own
import fnmatch
for w in [t for t in pairs if any(c in t for c in "*%")]:
    pat = w.replace("%", "?")
    for t in types:
        if t not in pairs and fnmatch.fnmatchcase(t, pat):
            pairs[t] = {**pairs[w], "name": t, "match": [t], "comment": f"Lennard-Jones of the wildcard entry {w}"}
    del pairs[w]
for t in list(pairs):
    if t not in types:
        foreign.add(t)
        del pairs[t]
if foreign:
    notes.append(f"entries naming types of other CHARMM files ({', '.join(sorted(foreign))}) left out: they apply only when those files are read too")
# OPLS in CHARMM's format (the Jorgensen lab's OPLS-AA/M files): the file's NONBONDED header says GEOM (geometric
# combination of sigma and epsilon) and E14FAC (1-4 Coulomb factor); 1-4 Lennard-Jones from the second (1-4) columns; no
# force switch (OPLS cuts plainly); impropers as the IMPR records write them, the out-of-plane atom first, the centre
# second, V2 [1 - cos 2 phi]
OPLS = "GEOM" in nb_settings
if OPLS:
    notes.append("OPLS settings from the file's NONBONDED header: geometric combination, 1-4 Coulomb x "
                 f"{nb_values.get('E14FAC', 1.0):g}, 1-4 Lennard-Jones from its 1-4 columns")
doc = {"format": "caps-forcefield", "format_version": 1, "name": arg("--name"), "source": arg("--source", ", ".join(p.split("/")[-1] for p in many("--prm"))),
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
                + (["CMAP corrections in the source are not represented"] if cmap_seen else []) + notes + many("--note"),
       "version": arg("--version", "")}
if arg("--typing"):
    doc["typing"] = arg("--typing")
if "--impropers-from-rtf" in sys.argv:
    import collections, itertools
    imp_rules = doc["impropers"]

    def charmm_lookup(ty):
        """the entry CHARMM uses for these four types: without wildcards (forward / backward) first, the last listed"""
        for wild in (False, True):
            for e in reversed(imp_rules):
                if ("*" in e["match"]) != wild:
                    continue
                for t in (ty, ty[::-1]):
                    if all(fnmatch.fnmatchcase(a, m) for a, m in zip(t, e["match"])):
                        return e["name"]
        return None

    instances, used, n_res, odd, reordered = [], set(), 0, [], []
    DEPTH = 16

    def tree(A, nbr, c, d, parent=None):
        """a centre's environment d bonds out, as CAPS reads it (ffdef.cpp improper_environment): the type, then the sorted
        environments of its neighbours other than the one it came from"""
        if d == 0:
            return A[c]
        return A[c] + "(" + ",".join(sorted(tree(A, nbr, x, d - 1, c) for x in nbr[c] if x != parent)) + ")"
    for rtf in many("--rtf"):
        res, cur = [], None
        for s_ in lines(rtf):
            w_ = s_.split()
            k_ = w_[0].upper()
            if k_ in ("RESI", "PRES"):
                cur = {"name": w_[1], "atoms": {}, "bonds": [], "impr": []} if k_ == "RESI" else None
                if cur is not None:
                    res.append(cur)
            elif cur is None:
                continue
            elif k_ == "ATOM" and len(w_) >= 3:
                cur["atoms"][w_[1].upper()] = w_[2]
            elif k_ in ("BOND", "DOUB", "DOUBLE", "TRIP", "TRIPLE", "AROM"):
                cur["bonds"] += [(a.upper(), b.upper()) for a, b in zip(w_[1::2], w_[2::2])]
            elif k_[:4] in ("IMPR", "IMPH"):
                cur["impr"] += [tuple(x.upper() for x in w_[q:q + 4]) for q in range(1, len(w_) - 3, 4)]
        for r in res:
            A = r["atoms"]
            if any(x[:1] in "+-" for b in r["bonds"] for x in b) or any(t not in types for t in A.values()):
                continue
            n_res += 1
            nbr = collections.defaultdict(set)
            for a, b in r["bonds"]:
                nbr[a].add(b), nbr[b].add(a)
            at = collections.defaultdict(list)
            for q in r["impr"]:
                if not all(x in A for x in q):
                    continue
                cs = [i for i in range(4) if all(q[j] in nbr[q[i]] for j in range(4) if j != i)]
                if len(cs) != 1 or len(nbr[q[cs[0]]]) != 3:
                    odd.append(f"{r['name']} {' '.join(q)}")
                    continue
                ty = [A[x] for x in q]
                at[q[cs[0]]].append((tuple(ty), cs[0]))
                hit = charmm_lookup(ty)
                if hit is None:   # an order no entry has (jul17's SAH): the entry CAPS then takes (ffdef.cpp, the same search)
                    others = sorted(ty[:cs[0]] + ty[cs[0] + 1:])
                    hit = next((h for cp in (0, 3) for p in itertools.permutations(others)
                                for h in [charmm_lookup(list(p[:cp]) + [ty[cs[0]]] + list(p[cp:]))] if h), None)
                    if hit:
                        reordered.append(f"{r['name']} {' '.join(q)}")
                used.add(hit)
            for c in A:
                if len(nbr[c]) == 3:
                    instances.append(([tree(A, nbr, c, d) for d in range(1, DEPTH + 1)], tuple(sorted(at.get(c, [])))))
    # the most common listing for each depth-1 environment (a centre and its three neighbours); then, deeper (the
    # neighbours' neighbours and so on), an exception wherever the residues showing a deeper environment list otherwise
    chosen, exceptions, parent_of, site_of, depth_used, left = {}, collections.defaultdict(list), {}, {}, 1, 0
    for d in range(1, DEPTH + 1):
        votes = collections.defaultdict(collections.Counter)
        for keys, listing in instances:
            votes[keys[d - 1]][listing] += 1
            site_of[keys[d - 1]] = keys[0]
            if d > 1:
                parent_of[keys[d - 1]] = keys[d - 2]
        for k, cnt in votes.items():
            top = max(cnt.values())
            best = [x for x, v in cnt.items() if v == top]
            inherit = chosen.get(parent_of.get(k)) if d > 1 else None
            pick = inherit if inherit in best else max(best, key=len)
            if d > 1 and pick != inherit:
                exceptions[site_of[k]].append((d, k, pick))
                depth_used = max(depth_used, d)
            chosen[k] = pick
        left = sum(1 for cnt in votes.values() if len(cnt) > 1)
        if not left:
            break
    as_json = lambda listing: [{"types": list(t), "centre": c} for t, c in listing]
    sites = []
    for k1 in sorted(k for k in chosen if "(" in k and k.count("(") == 1):
        best = chosen[k1]
        centre, nbs = k1[:k1.index("(")], k1[k1.index("(") + 1:-1].split(",")
        if not best and not exceptions.get(k1):   # a centre without impropers: kept only where some entry could put one
            if not any(charmm_lookup([centre] + list(q)) or charmm_lookup(list(q) + [centre]) for q in itertools.permutations(nbs)):
                continue
        site = {"centre": centre, "neighbours": nbs, "impropers": as_json(best)}
        if exceptions.get(k1):
            site["exceptions"] = [{"depth": d, "environment": k, "impropers": as_json(x)} for d, k, x in sorted(exceptions[k1])]
        sites.append(site)
    unl = [e["name"] for e in imp_rules if e["name"] not in used]
    doc["improper_placement"] = {"sites": sites, "unlisted_entries": unl}
    n_exc = sum(len(v) for v in exceptions.values())
    doc["notes"].append(f"impropers where the topology's IMPR records put them: {len(sites)} centre environments (a centre's type and "
                        f"its neighbours') from {n_res} residues of {', '.join(r.split('/')[-1] for r in many('--rtf'))}; where "
                        f"residues list one environment differently, the most common listing, with {n_exc} exception{'' if n_exc == 1 else 's'} for deeper "
                        f"environments (up to {depth_used} bonds out) that list otherwise"
                        + (f" ({left} environments still listed differently at {DEPTH} bonds)" if left else "")
                        + f"; the {len(unl)} parameter entries no IMPR record uses are never applied")
    if reordered:
        doc["notes"].append("IMPR records in an atom order no parameter entry has (CHARMM would stop on them); CAPS takes the "
                            "atoms in the first order an entry has (the centre first or last): " + "; ".join(reordered))
    if odd:
        doc["notes"].append("IMPR records whose centre is not a three-connected atom bonded to the other three, left out: " + "; ".join(odd))
if "--one-torsion-per-bond" in sys.argv:   # CHARMM19: "a single dihedral per bond rather than multiple torsions"
    doc["torsion_selection"] = "one_per_bond"
    doc["torsion_single_h"] = ["[#7X3;$(*-[#6]=[#8X1])]"]   # toph19: one torsion on an amide NH2, every H elsewhere
    doc["notes"].append("one torsion per central bond, as the topology lists them (aromatic ring bonds none)")
if not OPLS and abs(nb_values.get("E14FAC", 1.0) - 1.0) > 1e-9:   # CHARMM19: 1-4 Coulomb x 0.4
    doc["special_coul"] = [0, 0, nb_values["E14FAC"]]
    doc["notes"].append(f"1-4 Coulomb x {nb_values['E14FAC']:g} (the file's E14FAC)")
if OPLS and many("--rtf"):
    # OPLS gives each type its charge: the charge the type carries in the topology's residues (the most common one where
    # residues differ, listed)
    import collections as _c
    seen_q = _c.defaultdict(_c.Counter)
    for rtf in many("--rtf"):
        for s_ in lines(rtf):
            w_ = s_.split()
            if w_[0].upper() == "ATOM" and len(w_) >= 4 and re.match(r"^-?[\d.]+$", w_[3]):
                seen_q[w_[2]][round(float(w_[3]), 4)] += 1
    varied = []
    for t, cnt in seen_q.items():
        if t in types:
            q, _ = cnt.most_common(1)[0]
            types[t]["charge"] = q
            if len(cnt) > 1:
                varied.append(f"{t} {dict(cnt)}")
    doc["atom_types"] = [types[t] for t in sorted(types)]
    if varied:
        doc["notes"].append("types whose charge differs between residues (the most common taken): " + "; ".join(varied))
if OPLS:
    doc.update({"styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
                "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, nb_values.get("E14FAC", 1.0)],
                "improper_order": "center2_sorted"})
    # 1-4 columns that are every epsilon times one factor with the same sigma are a 1-4 scale (OPLS's 0.5)
    # types found only in residues of fewer than four atoms (TIP3P water, ions) are never in a 1-4 pair
    size_of, cur_ = {}, None
    for rtf in many("--rtf"):
        for s_ in lines(rtf):
            w_ = s_.split()
            if w_[0].upper() in ("RESI", "PRES"):
                cur_ = [] if w_[0].upper() == "RESI" else None
            elif w_[0].upper() == "ATOM" and cur_ is not None and len(w_) >= 3:
                cur_.append(w_[2])
                for t_ in cur_:
                    size_of[t_] = max(size_of.get(t_, 0), len(cur_))
    in14 = lambda t_: size_of.get(t_, 4) >= 4
    ratios = {round(p_["lj14"][0] / p_["params"][0], 6) for p_ in doc["pairs"] if "lj14" in p_ and p_["params"][0] > 0
              and abs(p_["lj14"][1] - p_["params"][1]) < 1e-6 and in14(p_["name"])}
    if len(ratios) == 1 and all("lj14" in p_ or p_["params"][0] == 0 for p_ in doc["pairs"] if in14(p_["name"])):
        f = ratios.pop()
        for p_ in doc["pairs"]:
            p_.pop("lj14", None)
        doc["special_lj"] = [0, 0, f]
        doc["notes"].append(f"the 1-4 Lennard-Jones columns are every epsilon x {f:g} with the same sigma: written as that 1-4 scale")
    doc.pop("pair_settings", None)
json.dump(doc, open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"{OUT}: {len(types)} types, {len(pairs)} LJ + {len(nbfix)} NBFIX, {len(bonds)} bonds, {len(angles)} angles, "
      f"{len(dihedrals)} dihedrals, {len(impropers)} impropers; {len(notes)} notes")
