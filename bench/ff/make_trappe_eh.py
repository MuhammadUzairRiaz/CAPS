#!/usr/bin/env python3
"""TraPPE-EH (explicit hydrogen) for CAPS from the TraPPE group's own parameter database (trappe.oit.umn.edu), model by
model, as bench/ff/make_trappe.py does for TraPPE-UA.

  make_trappe_eh.py SITE_DIR [--check] [--out DATA_DIR]

TraPPE-EH charges are not transferable (each model's own fit), so every model has its own types: one per set of
symmetry-equivalent sites of the model (its graph's automorphism classes), named <element><TraPPE type>_<model>[a-z].
Each type takes the Lennard-Jones terms of its TraPPE site type (vdw equivalence: C1_eh, H2_eh ...) and the model's own
charge. CAPS's typing rule for a type is the whole molecule written as SMARTS from that site (element and number of
partners for every atom, bonds of any order), so a rule matches only its own molecule.

  sites      the H of the alkanes, amines and amides is TraPPE-EH's interaction site on the C-H bond, 0.55 A from the C
             (the model's "M"): CAPS takes the H atom itself as that site, the C-M bond length 0.55 A
  pairs      epsilon / k_B (K) x k_B, sigma; Lorentz-Berthelot; 1-2, 1-3 and 1-4 pairs excluded (sites four or more bonds
             apart interact); 14 A cut-off, tail corrections, Ewald (the database's own notes)
  bonds      fixed lengths: stiff harmonic (600 kcal/mol/A^2, to be constrained); rigid models K 0 (CAPS's mark of a
             rigid molecule)
  bends      k_theta / 2 (theta - theta0)^2; a bend the model holds rigid in a flexible model is stiff (100 kcal/mol/rad^2),
             in a rigid model K 0
  torsions   every form of the database as Fourier terms of phi (as make_trappe.py); a torsion a model does not list: zero
  charges    where the database gives symmetry-equivalent sites different charges (phenol's ring, isoindole's H, 1,4-
             dihydroxybenzene, 1,3,5-trichlorobenzene: an in-plane OH or the fit's rounding) the class takes their mean,
             each model's total unchanged; the largest change is listed

--check: every EH model typed by CAPS's rules as the model's own sites, and every bonded term CAPS assigns (GROMACS
topology) the model's own (same atoms, same parameters).
"""
import collections, itertools, json, math, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
SITE = sys.argv[1]
OUT = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(ROOT, "data")
FID = "trappe-eh-umn"
KB = 0.0019872043   # kcal/mol/K
KBOND, KSTIFF = 600.0, 100.0
ELEM = {"M": "H"}
MASS = {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999, "S": 32.06, "F": 18.998, "Cl": 35.453, "Br": 79.904}
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "S": 16, "Cl": 17, "Br": 35}

pubs = {p["articleid"]: p for p in json.load(open(os.path.join(SITE, "publications.json")))["reg_pubs"]}


def cite(ref):
    p = pubs.get(ref)
    if not p:
        return f"TraPPE reference {ref}"
    a = p.get("authors") or []
    names = ", ".join(f"{x['initials']} {x['name']}" for x in a[:3]) + (" et al." if len(a) > 3 else "")
    return f"{names}, {p['journal']} {p['vol']}, {p['pages'].split('-')[0]} ({p['year']})"


# ---------------------------------------------------------------- the models
mols, CORRECTED, DROPPED = [], [], []
for e in json.load(open(os.path.join(SITE, "select.json")))["search"]:
    if e["family"] != "EH":
        continue
    d = json.load(open(os.path.join(SITE, "molecules", e["molecule_ID"] + ".json")))
    m = d["molname"][0]
    d["id"] = m["molecule_ID"]
    d["name"] = re.sub(r"<[^>]+>", "", m["name"])
    d["rigid"] = m["lrigid"] == "1"
    d["stretch"] = [x for x in d.get("stretch") or [] if x and x.get("A") is not None]
    d["bend"] = [x for x in d.get("bend") or [] if x and x.get("A") is not None]
    d["torsion_type"] = [T for T in d.get("torsion_type") or [] if T and T.get("torsion")]
    if d["name"] == "2-chlorothiophene":   # the database gives the ring S the furan O's site type (6); its bonds, SMILES and charges are thiophene's
        for a in d["nonbonded"]:
            if a["natom"] == "1" and a["nonbonded_ID"] == "6":
                a.update(nonbonded_ID="7", atomtype="C(aro)-[S](aro)-X(aro)", prtype="S", epsilon="158.00", sigma="3.55")
                CORRECTED.append("2-chlorothiophene: its site 1 (the ring S) carries the furan O's type 6 in the database; read as the thiophene S, type 7")
    d["atoms"] = {int(a["natom"]): a for a in d["nonbonded"]}
    d["el"] = {k: ELEM.get(a["prtype"], a["prtype"]) for k, a in d["atoms"].items()}
    adj = collections.defaultdict(set)
    for s in d["stretch"]:
        i, j = int(s["A"]), int(s["B"])
        adj[i].add(j)
        adj[j].add(i)
    d["adj"] = adj
    # the database's own slips: a bend or torsion listed with its atoms out of order (ethane's 1-2-3 for the M on C1):
    # read as the one bonded path through the same atoms
    for b in d["bend"]:
        q = [int(b[x]) for x in "ABC"]
        if not (q[0] in adj[q[1]] and q[2] in adj[q[1]]):
            fix = [p for p in itertools.permutations(q) if p[0] in adj[p[1]] and p[2] in adj[p[1]] and p[0] < p[2]]
            if len(fix) == 1:
                b["A"], b["B"], b["C"] = map(str, fix[0])
                CORRECTED.append(f"{d['name']}: bend {'-'.join(map(str, q))} read as {'-'.join(map(str, fix[0]))}")
            else:
                DROPPED.append(f"{d['name']}: bend {'-'.join(map(str, q))} (no bonded path through these atoms)")
                b["drop"] = True
    d["bend"] = [b for b in d["bend"] if not b.get("drop")]
    for T in d["torsion_type"]:
        for t in T["torsion"]:
            q = [int(t[x]) for x in "ABCD"]
            path = lambda p: p[1] in adj[p[0]] and p[2] in adj[p[1]] and p[3] in adj[p[2]] and len(set(p)) == 4
            if not path(q):
                fix = [p for p in itertools.permutations(q) if path(p) and p[0] < p[3]]
                if len(fix) == 1:
                    t["A"], t["B"], t["C"], t["D"] = map(str, fix[0])
                    CORRECTED.append(f"{d['name']}: torsion {'-'.join(map(str, q))} read as {'-'.join(map(str, fix[0]))}")
                else:
                    DROPPED.append(f"{d['name']}: torsion {'-'.join(map(str, q))} (no bonded path through these atoms)")
                    t["drop"] = True
        T["torsion"] = [t for t in T["torsion"] if not t.get("drop")]
    mols.append(d)


def classes(d):
    """symmetry classes of the model's sites: colour refinement on element and bonds (exact for these small molecules;
    checked below against every automorphism-related charge)"""
    col = {k: d["el"][k] for k in d["atoms"]}
    while True:
        sig = {k: (col[k], tuple(sorted(col[j] for j in d["adj"][k]))) for k in col}
        keys = {v: i for i, v in enumerate(sorted(set(sig.values()), key=str))}
        new = {k: keys[sig[k]] for k in col}
        if len(set(new.values())) == len(set(col.values())):
            return new
        col = new


def smarts_from(d, root):
    """the whole molecule as SMARTS, written from root: [#Z;Xn] atoms, ~ bonds, ring closures"""
    seen, order, parent = set(), [], {}
    closures = {}   # (i, j) back edge -> digit
    free = list(range(1, 100))

    def atom(k):
        return f"[#{Z[d['el'][k]]};X{len(d['adj'][k])}]"

    # depth-first tree, children in a fixed order
    tree = collections.defaultdict(list)
    back = []

    def dfs(k):
        seen.add(k)
        order.append(k)
        for j in sorted(d["adj"][k]):
            if j == parent.get(k):
                continue
            if j in seen:
                if (j, k) not in back:
                    back.append((k, j))
                continue
            parent[j] = k
            tree[k].append(j)
            dfs(j)
    dfs(root)
    # ring-closure digits: opened at the earlier atom, closed at the later one
    pos = {k: i for i, k in enumerate(order)}
    opens, closes = collections.defaultdict(list), collections.defaultdict(list)
    for a, b in back:
        first, last = (a, b) if pos[a] < pos[b] else (b, a)
        opens[first].append(last)
    out = []
    used = {}

    def emit(k):
        s = atom(k)
        for x in closes[k]:
            s += "~" + (str(x) if x < 10 else f"%{x}")
            free.append(x)
            free.sort()
        for last in opens[k]:
            x = free.pop(0)
            s += "~" + (str(x) if x < 10 else f"%{x}")
            closes[last].append(x)
        kids = tree[k]
        for c in kids[:-1]:
            s += "(~" + emit(c) + ")"
        if kids:
            s += "~" + emit(kids[-1])
        return s
    return emit(root)


# ---------------------------------------------------------------- types, charges, pairs, rules
lj = {}          # TraPPE site type id -> (prtype, eps K, sigma, atomtype, ref)
types, rules, notes_dev = [], [], []
tname = {}       # (model id, site) -> CAPS type
for d in mols:
    cl = classes(d)
    d["class"] = cl
    by = collections.defaultdict(list)
    for k in sorted(d["atoms"]):
        by[cl[k]].append(k)
    letters = iter("abcdefghijklmnopqrstuvwxyz")
    for c, ks in sorted(by.items(), key=lambda kv: min(kv[1])):
        a0 = d["atoms"][ks[0]]
        ids = {d["atoms"][k]["nonbonded_ID"] for k in ks}
        if len(ids) > 1:   # the same site in TraPPE's eyes under two type ids (equal parameters): the first
            pars = {(d["atoms"][k]["epsilon"], d["atoms"][k]["sigma"]) for k in ks}
            assert len(pars) == 1, (d["name"], ks)
        nb = a0["nonbonded_ID"]
        lj.setdefault(nb, (a0["prtype"], float(a0["epsilon"]), float(a0["sigma"]), a0["atomtype"], a0["ref_ID"]))
        qs = [float(d["atoms"][k]["charge"]) for k in ks]
        q = sum(qs) / len(qs)
        if max(qs) - min(qs) > 1e-9:
            notes_dev.append((d["name"], ks, qs, q))
        el = d["el"][ks[0]]
        name = f"{el}{nb}_{d['id']}{next(letters)}"
        for k in ks:
            tname[(d["id"], k)] = name
        types.append({"name": name, "element": el, "mass": MASS[el], "charge": round(q, 6),
                      "description": f"{d['name']}: {a0['atomtype']} (TraPPE-EH type {nb}); {cite(a0['q_ref_ID'] or a0['ref_ID'])}",
                      "equivalence": {"vdw": f"{el}{nb}_eh"}})
        rules.append({"type": name, "priority": 0, "smarts": smarts_from(d, ks[0]), "description": f"{d['name']}, sites {', '.join(map(str, ks))}"})
vdw_types = [{"name": f"{ELEM.get(p, p)}{nb}_eh", "element": ELEM.get(p, p), "mass": MASS[ELEM.get(p, p)], "charge": 0.0,
              "description": f"TraPPE-EH site type {nb}: {at} (Lennard-Jones class); {cite(ref)}"} for nb, (p, e, s, at, ref) in sorted(lj.items(), key=lambda x: int(x[0]))]
pairs = [{"name": f"{ELEM.get(p, p)}{nb}_eh", "match": [f"{ELEM.get(p, p)}{nb}_eh"], "params": [e * KB, s], "comment": f"eps/kB {e} K, sigma {s} A"}
         for nb, (p, e, s, at, ref) in sorted(lj.items(), key=lambda x: int(x[0]))]


# ---------------------------------------------------------------- torsion forms (as make_trappe.py)
def cos_power(n):
    out = collections.defaultdict(float)
    for k in range(n + 1):
        out[abs(n - 2 * k)] += math.comb(n, k) / 2 ** n
    return out


def fourier(form, t):
    c = lambda k: float(t[k]) * KB
    terms, const = [], 0.0
    if form == "opls":
        const = c("c0")
        terms = [(c("c1"), 1, 0.0), (c("c2"), 2, 180.0), (c("c3"), 3, 0.0)]
    elif form in ("C3", "C4", "C6"):
        const = c("c0")
        for n in range(1, int(form[1]) + 1):
            terms.append((c(f"c{n}"), n, 0.0))
            const -= c(f"c{n}")
    elif form == "aro":
        terms = [(c("e0"), 2, 180.0 - float(t["e1"]))]
    elif form == "EH":
        terms = [(c("c1"), 3, 180.0)]
    else:
        raise ValueError(form)
    terms = [(k, n, dd % 360.0) for k, n, dd in terms if abs(k) > 1e-12]
    if abs(const) > 1e-12:
        terms.append((const, 0, 90.0))
    return terms


def site_u(form, t, phi):
    g = lambda k: float(t[k])
    r = math.radians
    if form == "opls":
        return g("c0") + g("c1") * (1 + math.cos(r(phi))) + g("c2") * (1 - math.cos(r(2 * phi))) + g("c3") * (1 + math.cos(r(3 * phi)))
    if form in ("C3", "C4", "C6"):
        return g("c0") + sum(g(f"c{n}") * math.cos(r(n * phi)) for n in range(1, int(form[1]) + 1))
    if form == "aro":
        return g("e0") * (1 - math.cos(r(2 * phi) + r(g("e1"))))
    if form == "EH":
        return g("c1") * (1 - math.cos(r(3 * phi)))


worst = 0.0
for d in mols:
    for T in d["torsion_type"]:
        for t in T["torsion"]:
            fs = fourier(T["torsion_type"], t)
            for phi in range(0, 360, 5):
                u = sum(k * (1 + math.cos(math.radians(n * phi - dd))) for k, n, dd in fs)
                worst = max(worst, abs(u - site_u(T["torsion_type"], t, phi) * KB))
print(f"{len(mols)} TraPPE-EH models, {len(lj)} site types, {len(types)} CAPS types; torsion forms as Fourier terms: "
      f"largest difference {worst:.1e} kcal/mol")


# ---------------------------------------------------------------- bonded terms from each model's own lists
def key_sym(k):
    return min(tuple(k), tuple(k)[::-1])


BONDS, ANGLES, DIHEDRALS = {}, {}, {}
conflicts = []


def put(store, key, val, what):
    if key in store and store[key][0] != val[0]:
        conflicts.append((what, key, store[key][0], val[0]))
    store.setdefault(key, val)


for d in mols:
    T_ = lambda k: tname[(d["id"], k)]
    ref = lambda x: cite(x.get("ref_ID"))
    for s in d["stretch"]:
        k = key_sym((T_(int(s["A"])), T_(int(s["B"]))))
        put(BONDS, k, ((0.0 if d["rigid"] else KBOND, float(s["distance"])), ("rigid model: run it as a rigid body" if d["rigid"] else "fixed length in TraPPE") + f" ({ref(s)})"), d["name"])
    listed_b = set()
    for b in d["bend"]:
        q = (int(b["A"]), int(b["B"]), int(b["C"]))
        listed_b.add(key_sym(q))
        kt = float(b["ktheta"] or 0)
        K = kt / 2 * KB if kt > 0 else (0.0 if d["rigid"] else KSTIFF)
        put(ANGLES, key_sym(tuple(T_(x) for x in q)), ((round(K, 12), float(b["theta"])),
            (f"k_theta {kt:g} K/rad^2" if kt > 0 else "rigid model: run it as a rigid body" if d["rigid"] else "rigid in TraPPE (stiff here)") + f"; {ref(b)}"), d["name"])
    listed = set()
    for T in d["torsion_type"]:
        for t in T["torsion"]:
            q = (int(t["A"]), int(t["B"]), int(t["C"]), int(t["D"]))
            listed.add(key_sym(q))
            fs = tuple((round(a, 10), n, round(dd, 8)) for a, n, dd in fourier(T["torsion_type"], t))
            put(DIHEDRALS, key_sym(tuple(T_(x) for x in q)), (fs, f"TraPPE {T['torsion_type']} form; {ref(t)}"), d["name"])
    for j in d["adj"]:   # torsions the model does not list: none in TraPPE
        for k_ in d["adj"][j]:
            for i in d["adj"][j] - {k_}:
                for l in d["adj"][k_] - {j, i}:
                    if key_sym((i, j, k_, l)) not in listed:   # the database lists one of a set of equivalent torsions: theirs
                        DIHEDRALS.setdefault(key_sym(tuple(T_(x) for x in (i, j, k_, l))), ((), "no torsion in the TraPPE model"))
    # bends the model does not list (a rigid ring's): none in TraPPE
    for j in d["adj"]:
        for i in d["adj"][j]:
            for k_ in d["adj"][j]:
                if i < k_ and key_sym((i, j, k_)) not in listed_b:   # an equivalent listed bend's, else none
                    ANGLES.setdefault(key_sym((T_(i), T_(j), T_(k_))), ((0.0, 120.0), "not in the TraPPE model (rigid): none"))
if conflicts:
    print(f"{len(conflicts)} terms with two parameter sets in one model (first kept):")
    for c in conflicts[:20]:
        print("  ", c)
bonds = [{"name": "-".join(k), "match": list(k), "style": "harmonic", "params": list(v[0]), "comment": v[1]} for k, v in sorted(BONDS.items())]
angles = [{"name": "-".join(k), "match": list(k), "style": "harmonic", "params": list(v[0]), "comment": v[1]} for k, v in sorted(ANGLES.items())]
dihedrals = [{"name": "-".join(k), "match": list(k), "style": "fourier",
              "params": [len(v[0])] + [x for f in v[0] for x in f] if v[0] else [1, 0.0, 1, 0.0], "comment": v[1]} for k, v in sorted(DIHEDRALS.items())]
big = max(notes_dev, key=lambda x: max(abs(q - x[3]) for q in x[2])) if notes_dev else None
NOTES = ["built by CAPS (bench/ff/make_trappe_eh.py) from the TraPPE database's own EH models (trappe.oit.umn.edu), each term with its paper",
         "charges are each model's own (not transferable): one CAPS type per set of symmetry-equivalent sites of a model",
         "alkane, amine and amide H: the H atom is TraPPE-EH's C-H bond-centre site (C-M 0.55 A)",
         "1-2, 1-3 and 1-4 pairs excluded; 14 A cut-off with tail corrections and Ewald sums, as the database asks",
         "fixed bond lengths: stiff harmonic bonds (600 kcal/mol/A^2), to be constrained; rigid models carry K 0 (hold them rigid)"]
if conflicts:
    NOTES.append("terms the database gives two values in one model, for symmetry-equivalent atoms (CAPS keeps the first; ring bends of a "
                 "flexible ring within 1 deg, rigid models' geometry only, and the second OH of 1,3- and 1,4-dihydroxybenzene listed "
                 "rigid at 109.5 / 109.6 deg where every other phenolic C-O-H has 113 deg, 70420 K): " + "; ".join(
                     f"{m} {'-'.join(k)} {v1[1]:g} / {v2[1]:g} deg" for m, k, v1, v2 in conflicts))
NOTES += ["corrected: " + c for c in CORRECTED] + ["left out: " + c for c in DROPPED]
print(f"{len(CORRECTED)} corrections of the database's lists, {len(DROPPED)} entries left out")
for c in CORRECTED + DROPPED:
    print("  ", c)
if notes_dev:
    NOTES.append("symmetry-equivalent sites with different charges in the database take their mean: " + "; ".join(
        f"{n} sites {','.join(map(str, ks))} ({', '.join(f'{q:+.3f}' for q in qs)} -> {m:+.4f})" for n, ks, qs, m in notes_dev))
refs = sorted({cite(a["ref_ID"]) for d in mols for a in d["nonbonded"]} | {cite(a["q_ref_ID"]) for d in mols for a in d["nonbonded"] if a.get("q_ref_ID")})
doc = {"format": "caps-forcefield", "format_version": 1, "name": "TraPPE-EH", "units": "real",
       "source": "the TraPPE database (trappe.oit.umn.edu): the Siepmann group's tables of their published parameters, model by model (explicit-hydrogen family)",
       "references": refs, "version": "TraPPE database, retrieved 2026-10",
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "harmonic"},
       "mixing": "arithmetic", "special_lj": [0, 0, 0], "special_coul": [0, 0, 0], "cutoff": 14,
       "atom_types": vdw_types + types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": [],
       "notes": NOTES, "typing": f"../typing/{FID}.typing.json"}
ff_path = os.path.join(OUT, "forcefields", FID + ".json")
ty_path = os.path.join(OUT, "typing", FID + ".typing.json")
json.dump(doc, open(ff_path, "w"), indent=1, ensure_ascii=False)
json.dump({"format": "caps-typing", "version": 1, "forcefield": "TraPPE-EH", "unknown_types": "untyped", "bonds": "defined",
           "description": "CAPS rules for the TraPPE-EH models, one whole-molecule SMARTS per symmetry class (bench/ff/make_trappe_eh.py)",
           "rules": rules}, open(ty_path, "w"), indent=1, ensure_ascii=False)
print(f"{ff_path}: {len(doc['atom_types'])} types, {len(bonds)} bonds, {len(angles)} bends, {len(dihedrals)} torsions")
if big:
    print(f"largest charge averaging: {big[0]} sites {big[1]} {big[2]} -> {big[3]:.4f}")


# ---------------------------------------------------------------- reference structures and the checks
def mol2(d, path):
    """the model's atoms and bonds; positions from its molfile where its atoms map onto the model's"""
    n = max(d["atoms"])
    ks = sorted(d["atoms"])
    idx = {k: i + 1 for i, k in enumerate(ks)}
    bonds = [(idx[int(s["A"])], idx[int(s["B"])]) for s in d["stretch"]]
    xyz = {i + 1: (i * 1.3, (i % 3) * 0.9, (i % 2) * 0.7) for i in range(len(ks))}
    mf = os.path.join(SITE, "molfiles", d["molname"][0]["molfile"])
    if os.path.exists(mf):
        L = open(mf).read().splitlines()
        if len(L) > 4 and L[3][:3].strip().isdigit() and int(L[3][:3]) == len(ks):
            na, nbd = len(ks), int(L[3][3:6])
            pos = [tuple(float(x) * 1.5 for x in L[4 + k].split()[:3]) for k in range(na)]
            el = [L[4 + k].split()[3] for k in range(na)]
            mb = {tuple(sorted((int(L[4 + na + b][:3]) - 1, int(L[4 + na + b][3:6]) - 1))) for b in range(nbd)}
            sb = {tuple(sorted((i - 1, j - 1))) for i, j in bonds}
            want = [d["el"][k] for k in ks]
            mapping = [-1] * na

            def place(k, used):
                if k == na:
                    return True
                for x in range(na):
                    if x in used or el[x] != want[k]:
                        continue
                    if all(((tuple(sorted((k, j))) in sb) == (tuple(sorted((x, mapping[j]))) in mb)) for j in range(k)):
                        mapping[k] = x
                        if place(k + 1, used | {x}):
                            return True
                mapping[k] = -1
                return False
            if len(mb) == len(sb) and place(0, frozenset()):
                xyz = {k + 1: pos[mapping[k]] for k in range(na)}
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{d['name']}\n{len(ks)} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, k in enumerate(ks):
            e = d["el"][k]
            f.write(f"{i + 1} {e}{i + 1} {xyz[i + 1][0]:.4f} {xyz[i + 1][1]:.4f} {xyz[i + 1][2]:.4f} {e} 1 MOL 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for b, (i, j) in enumerate(bonds):
            f.write(f"{b + 1} {i} {j} 1\n")
    return ks


if "--check" in sys.argv:
    work = tempfile.mkdtemp()
    right = total = 0
    bad = []
    for d in mols:
        p = os.path.join(work, f"m{d['id']}.mol2")
        ks = mol2(d, p)
        out = p + ".types"
        r = subprocess.run([CAPS, "ff", "type", p, "--ff", ff_path, "-o", out], capture_output=True, text=True)
        if not os.path.exists(out):
            bad.append((d["name"], "failed: " + (r.stderr or r.stdout).strip()[-160:]))
            continue
        got = [x.strip() for x in open(out)][:len(ks)]
        want = [tname[(d["id"], k)] for k in ks]
        total += len(want)
        diff = [(k, w, g) for k, w, g in zip(ks, want, got) if w != g]
        right += len(want) - len(diff)
        if diff:
            bad.append((d["name"], ", ".join(f"{k} {w}->{g}" for k, w, g in diff[:6])))
    print(f"typing: {right}/{total} sites as the TraPPE database types them; {len(bad)} models differ")
    for n_, b in bad:
        print("  ", n_, ":", b)

    def itp_terms(path):
        sec, out = None, {"bonds": {}, "angles": {}, "dihedrals": collections.defaultdict(list), "pairs": set()}
        for l in open(path):
            l = l.split(";")[0].strip()
            m = re.match(r"\[\s*(\S+)\s*\]", l)
            if m:
                sec = m.group(1)
                continue
            w = l.split()
            if not w:
                continue
            if sec == "bonds":
                out["bonds"][key_sym((int(w[0]), int(w[1])))] = (float(w[3]) * 10, float(w[4]) / 4.184 / 100 / 2 if len(w) > 4 else 0.0)
            elif sec == "angles":
                out["angles"][key_sym((int(w[0]), int(w[1]), int(w[2])))] = (float(w[4]), float(w[5]) / 4.184 / 2)
            elif sec == "dihedrals" and w[4] in ("1", "9"):
                out["dihedrals"][key_sym((int(w[0]), int(w[1]), int(w[2]), int(w[3])))].append((round(float(w[6]) / 4.184, 6), int(w[7]), round(float(w[5]) % 360, 4)))
            elif sec == "pairs":
                out["pairs"].add(key_sym((int(w[0]), int(w[1]))))
        return out
    ok = 0
    bad = []
    for d in mols:
        p = os.path.join(work, f"m{d['id']}.mol2")
        ks = sorted(d["atoms"])
        I = {k: i + 1 for i, k in enumerate(ks)}
        stem = p + ".gmx"
        r = subprocess.run([CAPS, "ff", "apply", p, "--ff", ff_path, "--charges", "types", "--gromacs", stem], capture_output=True, text=True)
        if not os.path.exists(stem + ".itp"):
            bad.append((d["name"], (r.stderr or r.stdout).strip()[-200:]))
            continue
        got = itp_terms(stem + ".itp")
        errs = []
        for s in d["stretch"]:
            g = got["bonds"].get(key_sym((I[int(s["A"])], I[int(s["B"])])))
            if not g or abs(g[0] - float(s["distance"])) > 1e-6:
                errs.append(f"bond {s['A']}-{s['B']} {g}")
        for b in d["bend"]:
            kt = float(b["ktheta"] or 0)
            K = kt / 2 * KB if kt > 0 else (0.0 if d["rigid"] else KSTIFF)
            g = got["angles"].get(key_sym((I[int(b["A"])], I[int(b["B"])], I[int(b["C"])])))
            if not g or abs(g[0] - float(b["theta"])) > 1e-6 or abs(g[1] - K) > 1e-6 * max(1, K):
                errs.append(f"bend {b['A']}-{b['B']}-{b['C']} {g} vs {(float(b['theta']), K)}")
        for T in d["torsion_type"]:
            for t in T["torsion"]:
                q = key_sym((I[int(t["A"])], I[int(t["B"])], I[int(t["C"])], I[int(t["D"])]))
                want = sorted((round(a, 6), n, round(dd % 360, 4)) for a, n, dd in fourier(T["torsion_type"], t))
                g = sorted(got["dihedrals"].get(q, []))
                g = [x for x in g if abs(x[0]) > 1e-9]
                if len(g) != len(want) or any(abs(a[0] - b[0]) > 1e-5 or a[1] != b[1] or abs(a[2] - b[2]) > 1e-3 for a, b in zip(g, want)):
                    errs.append(f"torsion {t['A']}-{t['B']}-{t['C']}-{t['D']} {g} vs {want}")
        if got["pairs"]:   # written with epsilon 0 and fudgeQQ 0: no 1-4 interaction
            itp = open(stem + ".itp").read()
            if re.search(r"^\s*\d+\s+\d+\s+1\s+\S+\s+(?!0\s*$)\S+\s*$", itp.split("[ pairs ]")[1].split("[")[0], re.M) or \
                    not re.search(r"^\s*1\s+2\s+no\s+\S+\s+0\s*$", open(stem + ".top").read(), re.M):
                errs.append("1-4 pairs interact (TraPPE excludes them)")
        if errs:
            bad.append((d["name"], "; ".join(errs[:3]) + (" ..." if len(errs) > 3 else "")))
        else:
            ok += 1
    print(f"bonded terms: {ok}/{len(mols)} models get exactly their own bonds, bends and torsions")
    for n_, b in bad:
        print("  ", n_, ":", b)


# --dl-map: the reference keys of DL_FIELD's TraPPE-EH templates for CAPS's types, by matching each template's graph onto
# the database model of the same molecule (elements and bonds only; CAPS's typing plays no part), for
# bench/ff/validate_typing.py --map trappe-eh-umn
if "--dl-map" in sys.argv:
    sys.path.insert(0, HERE)
    import template_pdb
    lib = os.path.join(os.environ.get("DLFIELD", os.path.expanduser("~/project/dl_f_4.13")), "lib", "TRAPPE_EH.sf")
    tel_of, keys, tmpl = template_pdb.parse_sf(lib)   # DL type -> element
    by_name = {d["name"]: d for d in mols}
    cmap = {}
    for tn, t in tmpl.items():
        d = by_name.get(tn)
        if not d:
            print("no database model for template", tn)
            continue
        labels = [a for a, _ in t["atoms"]]
        tadj = collections.defaultdict(set)
        for a, js in t["connect"].items():
            for b in js:
                if a != b:
                    tadj[a].add(b); tadj[b].add(a)
        tel = {a: tel_of[ty] for a, ty in t["atoms"]}
        ks = sorted(d["atoms"])
        mapping = {}

        def place(i, used):
            if i == len(labels):
                return True
            a = labels[i]
            for k in ks:
                if k in used or d["el"][k] != tel[a]:
                    continue
                if all((b in tadj[a]) == (mapping[b] in d["adj"][k]) for b in labels[:i]):
                    mapping[a] = k
                    if place(i + 1, used | {k}):
                        return True
            mapping.pop(a, None)
            return False
        if not place(0, frozenset()):
            print("template", tn, "does not match the database model")
            continue
        for a, ty in t["atoms"]:
            cmap.setdefault(tname[(d["id"], mapping[a])], template_pdb.ATOM_KEY[ty])
    mp = os.path.join(ROOT, "bench", "typing", "reference_type_maps.json")
    allm = json.load(open(mp))
    allm[FID] = {"lib": "TRAPPE_EH", "map": cmap}
    json.dump(allm, open(mp, "w"), indent=1)
    print(f"{mp}: {len(cmap)} CAPS types mapped to TRAPPE_EH keys")
