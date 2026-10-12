#!/usr/bin/env python3
"""TraPPE-UA for CAPS from the TraPPE group's own parameter database (trappe.oit.umn.edu, the Siepmann group's
tables of their published parameters, each with its paper), transcribed molecule by molecule.

  make_trappe.py SITE_DIR [--check] [--out DATA_DIR]

SITE_DIR holds what the site serves: select.json (the list of models), molecules/ID.json (each model's sites with
their TraPPE type, epsilon / k_B, sigma, charge, and its bonds, bends and torsions by type and number, each with its
reference), publications.json. CAPS keeps its own typing rules (one per TraPPE site type, SMARTS on the united-atom
structure) and builds the force field from the models:

  pairs      epsilon / k_B (K) x k_B, sigma; Lorentz-Berthelot; charges per site type; 1-4 excluded (the models that
             count 1-4 Coulomb x 0.5 are listed: CAPS's 1-4 rule is per force field)
  bonds      TraPPE's fixed lengths, as stiff harmonic bonds (K 600 kcal/mol/A^2; constrain them in MD)
  bends      k_theta / 2 (theta - theta0)^2: harmonic K = k_theta / 2 x k_B
  torsions   every functional form of the site written exactly as Fourier terms of phi (IUPAC: trans 180 deg); the
             site's psi = phi - 180 deg where it uses psi
  terms      each bend and torsion keyed by the CAPS types of its atoms as the models list them; a centre (bend) or
             central bond (torsion) with one parameter set in every model it appears in gets a wildcard rule

--check: every UA model's sites typed by CAPS's rules must be the site's types, and every bonded term CAPS assigns
the site's (same atoms, same parameters).
"""
import collections, glob, json, math, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
SITE = sys.argv[1]
OUT = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(ROOT, "data")
KB = 0.0019872043   # kcal/mol/K

# ---------------------------------------------------------------- CAPS's typing rules, one per TraPPE site type
# (id, priority, SMARTS on the united-atom structure, implicit H counts as perceived from the CHn sites)
RULES = [
    (1, 0, "[CX4H4]"),
    (2, 0, "[CX4H3]"), (3, 0, "[CX4H2]"), (4, 0, "[CX4H1]"), (5, 0, "[CX4H0]"),
    # alcohols
    (6, 1, "[OX2H1][CX4]"), (7, 1, "[H][OX2H1][CX4]"),
    (8, 2, "[CX4H3][OX2H1]"), (9, 2, "[CX4H2][OX2H1]"), (10, 2, "[CX4H1][OX2H1]"), (11, 2, "[CX4H0][OX2H1]"),
    # aldehydes
    (12, 2, "[CX4H3][CX3H1]=[OX1]"), (13, 2, "[CX4H2][CX3H1]=[OX1]"), (21, 2, "[CX3H1](=[OX1])[#6]"), (22, 2, "[OX1]=[CX3H1][#6]"),
    # ethers and the alkyl side of esters
    (16, 2, "[CX4H3][OX2][#6]"), (17, 2, "[CX4H2][OX2][#6]"), (18, 2, "[CX4H1][OX2][#6]"), (19, 2, "[CX4H0][OX2][#6]"),
    (20, 1, "[OX2]([CX4])[CX4]"),
    # ketones
    (23, 2, "[CX3](=[OX1])([#6])[#6]"), (24, 2, "[OX1]=[CX3]([#6])[#6]"),
    # alkenes: CH2=, CH= and C= next to an sp3 carbon, or (conjugated) to an sp2 one
    (29, 1, "[CX3H2]=[#6]"), (30, 1, "[CX3H1]=[#6]"), (31, 1, "[CX3H0]=[#6]"),
    (138, 2, "[CX3H1](=[#6])[#6X3]"), (139, 2, "[CX3H0](=[#6])[#6X3]"),
    # nitriles
    (47, 2, "[CX2](#[NX1])[CX4]"), (48, 2, "[NX1]#[CX2][CX4]"), (49, 3, "[CX4H3][CX2]#[NX1]"), (50, 3, "[CX4H2][CX2]#[NX1]"),
    # nitro
    (53, 2, "[OX1]~[NX3](~[OX1])[#6]"), (54, 2, "[NX3](~[OX1])(~[OX1])[#6]"), (59, 3, "[c][NX3](~[OX1])~[OX1]"),
    # aromatics
    (88, 1, "[cH1]"), (89, 2, "[cH0]-[CX4]"), (90, 2, "[cH0;$(*(:c)(:c)-c),$(*(:c)(:c):c)]"),
    (60, 2, "[nX2;r6]"), (61, 3, "[cH1;r6]:[nX2]"), (62, 4, "[cH1;r6](:[nX2]):[nX2]"),
    (79, 2, "[sX2;r5]"), (80, 3, "[cH1;r5]:[sX2]"),
    # thiols, sulfides, disulfides
    (63, 2, "[SX2H1][CX4]"), (64, 2, "[H][SX2][CX4]"),
    (65, 3, "[CX4H3][SX2H1]"), (66, 3, "[CX4H2][SX2H1]"), (67, 3, "[CX4H1][SX2H1]"), (68, 3, "[CX4H0][SX2H1]"),
    (69, 2, "[SX2]([CX4])[CX4]"), (70, 3, "[CX4H3][SX2][CX4]"), (71, 3, "[CX4H2][SX2][CX4]"),
    (74, 3, "[SX2]([CX4])[SX2]"), (75, 4, "[CX4H3][SX2][SX2]"), (76, 4, "[CX4H2][SX2][SX2]"),
    # esters
    (99, 3, "[CX3](=[OX1])[OX2][#6]"), (100, 3, "[OX1]=[CX3][OX2][#6]"), (101, 3, "[OX2]([#6])[CX3]=[OX1]"),
    # cycloalkanes: CH2 in a ring of 5, 6, 7 or more
    (141, 1, "[CX4H2;r5]"), (140, 1, "[CX4H2;r6]"), (110, 1, "[CX4H2;R;!r3;!r4;!r5;!r6]"),
    # cyclic ethers: tetrahydrofuran, 1,3-dioxolane, tetrahydropyran, 1,4- and 1,3-dioxane, 1,3,5-trioxane
    (144, 3, "[OX2;r5]([CX4H2])[CX4H2]"), (142, 3, "[CX4H2;r5][OX2;r5]"), (143, 3, "[CX4H2;r5;$(*[CX4H2][OX2])]"),
    (145, 4, "[OX2;r5;$(*[CX4H2][OX2])]"), (146, 4, "[CX4H2;r5]([OX2])[OX2]"), (147, 4, "[CX4H2;r5;$(*[OX2][CX4H2][OX2])]"),
    (148, 3, "[OX2;r6]([CX4H2])[CX4H2]"), (149, 3, "[CX4H2;r6][OX2;r6]"), (150, 3, "[CX4H2;r6;$(*[CX4H2][OX2])]"),
    (151, 4, "[OX2;r6;$(*[CX4H2][CX4H2][OX2])]"), (152, 4, "[CX4H2;r6;$(*[OX2])$(*[CX4H2][OX2])]"),
    (153, 4, "[OX2;r6;$(*[CX4H2][OX2])]"), (154, 4, "[CX4H2;r6]([OX2])[OX2]"), (155, 4, "[CX4H2;r6;$(*[OX2][CX4H2][OX2]);!$(*([OX2])[OX2])]"),
    (156, 4, "[CX4H2;r6;$(*([CX4H2][OX2])[CX4H2][OX2])]"),
    (157, 5, "[OX2;r6;$(*[CX4H2][OX2][CX4H2][OX2])]"), (158, 5, "[CX4H2;r6;$(*([OX2])[OX2][CX4H2][OX2])]"),
    # organophosphorus: DMMP, sarin, soman
    (115, 3, "[PX4](=[OX1])([CX4H3])([OX2])[OX2]"), (114, 3, "[OX1]=[PX4]([CX4H3])([OX2])[OX2]"), (113, 3, "[OX2]([CX4H3])[PX4]"),
    (111, 3, "[CX4H3][OX2][PX4]"), (112, 3, "[CX4H3][PX4]([OX2])[OX2]"),
    (117, 3, "[PX4](=[OX1])F"), (116, 3, "[F][PX4]"), (118, 3, "[OX1]=[PX4]F"), (119, 3, "[OX2]([CX4])[PX4]F"), (120, 4, "[CX4H3][PX4]F"),
    (121, 4, "[CX4H3][CX4H1][OX2][PX4]"), (122, 4, "[CX4H1]([CX4H3])([CX4H3])[OX2][PX4]"), (123, 4, "[CX4H1]([CX4H3])([CX4H0])[OX2][PX4]"),
    (124, 4, "[CX4H3][CX4H0][CX4H1][OX2][PX4]"), (125, 4, "[CX4H0]([CX4H3])([CX4H3])([CX4H3])[CX4H1][OX2][PX4]"),
    # TraPPE-small (rigid): carbon dioxide, nitrogen, oxygen, ammonia, hydrogen sulfide, ethylene oxide
    (86, 3, "[CX2](=[OX1])=[OX1]"), (87, 3, "[OX1]=[CX2]=[OX1]"), (84, 3, "[NX1]#[NX1]"), (108, 3, "[OX1]=[OX1]"),
    (126, 3, "[NX3H3]"), (127, 3, "[H][NX3H3]"), (159, 3, "[SX2H2]"), (160, 3, "[H][SX2H2]"),
    (106, 4, "[OX2;r3]([CX4H2])[CX4H2]"), (107, 4, "[CX4H2;r3][OX2;r3]"),
]


# the massless charge sites of the small models, as weighted sums of their atoms (CAPS's virtual sites): N2's and O2's
# midway between the atoms; NH3's on the axis 0.08 A from N towards the H, H2S's on the bisector 0.3 A from S
def _axis_w(d_site, bond, angle_deg, nh):
    """weights (centre, each H) of a site d_site from the centre towards the centroid of nh H atoms"""
    a = math.radians(angle_deg)
    cos_ax = math.sqrt((1 + 2 * math.cos(a)) / 3) if nh == 3 else math.cos(a / 2)
    f = d_site / (bond * cos_ax)
    return [1 - f] + [f / nh] * nh


VSITES = [
    (85, "[NX1]#[NX1]", [0, 1], [0.5, 0.5], "N2's centre charge, midway between the atoms"),
    (109, "[OX1]=[OX1]", [0, 1], [0.5, 0.5], "O2's centre charge, midway between the atoms"),
    (128, "[NX3]([H])([H])[H]", [0, 1, 2, 3], _axis_w(0.08, 1.0124, 106.68, 3), "NH3's charge on the axis, 0.08 A from N towards the H"),
    (161, "[SX2]([H])[H]", [0, 1, 2], _axis_w(0.3, 1.34, 92.0, 2), "H2S's charge on the bisector, 0.3 A from S"),
]
MASS = {"M": 0.0, "CH4": 16.043, "CH3": 15.035, "CH2": 14.027, "CH": 13.019, "C": 12.011, "O": 15.999, "N": 14.007, "S": 32.06, "P": 30.974,
        "F": 18.998, "H": 1.008}
ELEM = {"CH4": "C", "CH3": "C", "CH2": "C", "CH": "C", "C": "C"}
NAMEH = {"CH4": "CH4", "CH3": "CH3", "CH2": "CH2", "CH": "CH", "C": "C"}

sel = json.load(open(os.path.join(SITE, "select.json")))["search"]
mols = []
for e in sel:
    d = json.load(open(os.path.join(SITE, "molecules", e["molecule_ID"] + ".json")))
    m = d["molname"][0]
    if m["family"] not in ("UA", "small"):
        continue
    # the site's own gaps: a row with no site type (some esters) makes the model incomplete
    d["incomplete"] = any(a.get("nonbonded_ID") is None for a in d["nonbonded"])
    d["nonbonded"] = [a for a in d["nonbonded"] if a.get("nonbonded_ID") is not None]
    for k in ("stretch", "bend"):
        d[k] = [x for x in (d.get(k) or []) if x and x.get("A") is not None]
    d["torsion_type"] = [T for T in (d.get("torsion_type") or []) if T and T.get("torsion")]
    for T in d["torsion_type"]:
        T["torsion"] = [t for t in T["torsion"] if t and t.get("A") is not None]
    d["name"] = re.sub(r"<[^>]+>", "", m["name"])
    keep = sorted(int(a["natom"]) for a in d["nonbonded"] if a["prtype"] != "M")
    d["idx"] = {k: i + 1 for i, k in enumerate(keep)}   # the model's atom number -> the reference structure's
    d["nsites"] = sum(a["prtype"] == "M" for a in d["nonbonded"])
    if d["name"] == "carbon dioxide":   # the database's bond list gives 1-2 and 2-3 (O-O); carbon 1 holds both oxygens
        d["stretch"] = [x if (x["A"], x["B"]) != ("2", "3") else {**x, "A": "1", "B": "3"} for x in d["stretch"]]
    d["id"] = m["molecule_ID"]
    mols.append(d)
pubs = {p["articleid"]: p for p in json.load(open(os.path.join(SITE, "publications.json")))["reg_pubs"]}


def cite(ref):
    p = pubs.get(ref)
    if not p:
        return f"TraPPE reference {ref}"
    a = p.get("authors") or []
    names = ", ".join(f"{x['initials']} {x['name']}" for x in a[:3]) + (" et al." if len(a) > 3 else "")
    return f"{names}, {p['journal']} {p['vol']}, {p['pages'].split('-')[0]} ({p['year']})"


# which models CAPS can hold: no massless charge sites (M), no united CFn sites
def usable(d):
    """CAPS holds the model: no united CFn sites, and no quadrupole charge sites out of the ring plane (TraPPE's
    benzene, toluene, pyridine, pyrimidine and nitrobenzene with explicit quadrupoles)"""
    quad = any(m and m.get("message_ID") == "4" for m in d.get("messages") or [])
    return not quad and all(a["prtype"] in MASS or a["prtype"] in ELEM for a in d["nonbonded"])


types = {}
for d in mols:
    for a in d["nonbonded"]:
        i = int(a["nonbonded_ID"])
        types.setdefault(i, a)
rule_ids = {r[0] for r in RULES} | {v[0] for v in VSITES}
print(f"{len(mols)} TraPPE-UA models; {len(types)} site types; rules for {len(rule_ids & set(types))}; without a rule: "
      + ", ".join(f"{i} {types[i]['atomtype']}" for i in sorted(set(types) - rule_ids)))

# ---------------------------------------------------------------- library: types and pairs
def tname(i):
    a = types[i]
    return f"{a['prtype']}_{i}"


lib_types, lib_pairs = [], []
for i in sorted(types):
    a = types[i]
    if a["prtype"] not in MASS:
        continue
    lib_types.append({"name": tname(i), "element": "" if a["prtype"] == "M" else ELEM.get(a["prtype"], a["prtype"]), "mass": MASS[a["prtype"]], "charge": float(a["charge"]),
                      "description": f"{a['atomtype']} ({a['functionalgroup']}); {cite(a['ref_ID'])}"})
    lib_pairs.append({"name": tname(i), "match": [tname(i)], "params": [float(a["epsilon"]) * KB, float(a["sigma"])]})
rules = [{"type": tname(i), "priority": p, "smarts": sm, "description": types[i]["atomtype"] + f" ({types[i]['functionalgroup']})"}
         for i, p, sm in RULES if i in types]
# the charge sites CAPS adds (element 0, named by their type)
rules += [{"type": tname(i), "priority": 5, "smarts": "[#0]", "atom_name": tname(i), "description": c} for i, sm, fr, w, c in VSITES if i in types]


# ---------------------------------------------------------------- reference structures (the site's own models)
def mol2(d, path):
    """the model's sites and bonds (united atoms named CH3, CH2 ... so CAPS reads their hydrogens)"""
    I = d["idx"]
    at = {I[int(a["natom"])]: a for a in d["nonbonded"] if a["prtype"] != "M"}
    n = max(at)
    bonds = [(I[int(s["A"])], I[int(s["B"])]) for s in d["stretch"] if int(s["A"]) in I and int(s["B"]) in I]
    # positions from the model's own drawing (its molfile, atoms in the model's order) where it has one: a cis isomer
    # drawn cis; 2D bond lengths of 1 scaled to 1.5 A
    xyz = {k: (k * 1.3, (k % 3) * 0.9, 0.0) for k in range(1, n + 1)}
    mf = os.path.join(SITE, "molfiles", d["molname"][0]["molfile"])
    if os.path.exists(mf):
        L = open(mf).read().splitlines()
        if len(L) > 4 and L[3][:3].strip().isdigit() and int(L[3][:3]) == n:
            nbonds = int(L[3][3:6])
            pos = [tuple(float(x) * 1.5 for x in L[4 + k].split()[:3]) for k in range(n)]
            el = [L[4 + k].split()[3] for k in range(n)]
            mb = {tuple(sorted((int(L[4 + n + b][:3]) - 1, int(L[4 + n + b][3:6]) - 1))) for b in range(nbonds)}
            sb = {tuple(sorted((i - 1, j - 1))) for i, j in bonds}
            want = [ELEM.get(at[k + 1]["prtype"], at[k + 1]["prtype"]) for k in range(n)]
            # the drawing's atoms in another order: map them onto the model's by element and bonds
            mapping = [-1] * n

            def place(k, used):
                if k == n:
                    return True
                for x in range(n):
                    if x in used or el[x] != want[k]:
                        continue
                    if all(((tuple(sorted((k, j))) in sb) == (tuple(sorted((x, mapping[j]))) in mb)) for j in range(k)):
                        mapping[k] = x
                        if place(k + 1, used | {x}):
                            return True
                mapping[k] = -1
                return False
            if len(mb) == len(sb) and place(0, frozenset()):
                xyz = {k + 1: pos[mapping[k]] for k in range(n)}
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{d['name']}\n{n} {len(bonds)} 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n")
        for k in range(1, n + 1):
            p = at[k]["prtype"]
            # mol2's type column is the name CAPS reads (CH3, CH2 ...: united atoms with their hydrogens)
            f.write(f"{k} {ELEM.get(p, p)}{k} {xyz[k][0]:.4f} {xyz[k][1]:.4f} {xyz[k][2]:.4f} {NAMEH.get(p, p)} 1 MOL 0.0\n")
        f.write("@<TRIPOS>BOND\n")
        for k, (i, j) in enumerate(bonds):
            f.write(f"{k + 1} {i} {j} un\n")


def write_library(path_ff, path_ty):
    doc = {"format": "caps-forcefield", "format_version": 1, "name": "TraPPE-UA", "units": "real",
           "source": "the TraPPE database (trappe.oit.umn.edu): the Siepmann group's tables of their published parameters, model by model",
           "references": sorted({cite(a["ref_ID"]) for d in mols for a in d["nonbonded"]}),
           "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "harmonic"},
           "mixing": "arithmetic", "special_lj": [0, 0, 0], "special_coul": [0, 0, 0], "cutoff": 14, "equivalence": "replace",
           "coul14_rules": [{"in_molecule": ["[OX2][CX4;!R][CX4;!R][OX2]", "[OX2H1][CX4][CX4][CX4][OX2H1]", "[SX2][SX2]",
                                             "[CX3H2]=[CX3][CX3](=[OX1])[OX2]"], "coul14": 0.5,
                             "comment": "TraPPE's glycol ethers, diols, disulfides and (meth)acrylates count 1-4 Coulomb x 0.5"}],
           "virtual_sites": [{"smarts": sm, "type": tname(i), "from": fr, "w": [round(x, 12) for x in w], "comment": c}
                             for i, sm, fr, w, c in VSITES if i in types],
           "atom_types": lib_types, "pairs": lib_pairs, "bonds": BONDS, "angles": WILD_ANGLES + ANGLES, "dihedrals": DIHEDRALS, "impropers": IMPROPERS,
           "notes": NOTES, "version": "TraPPE database, retrieved 2026-10", "typing": "../typing/trappe-ua-umn.typing.json"}
    json.dump(doc, open(path_ff, "w"), indent=1, ensure_ascii=False)
    json.dump({"format": "caps-typing", "version": 1, "forcefield": "TraPPE-UA", "unknown_types": "untyped", "united_atom": True,
               "description": "CAPS rules for the TraPPE-UA site types (bench/ff/make_trappe.py), checked on the TraPPE database's own models",
               "rules": rules}, open(path_ty, "w"), indent=1, ensure_ascii=False)


# ---------------------------------------------------------------- torsions: every site form as Fourier terms of phi
def cos_power(n):
    """cos^n(phi) as sum_m b_m cos(m phi)"""
    out = collections.defaultdict(float)
    for k in range(n + 1):
        out[abs(n - 2 * k)] += math.comb(n, k) / 2 ** n
    return out


def fourier(form, t):
    """[(K, n, d)] with sum K (1 + cos(n phi - d)) = the site's u(phi) (kcal/mol); a constant as n = 0, d = 90 (K)"""
    c = lambda k: float(t[k]) * KB
    terms, const = [], 0.0
    if form == "opls":
        const = c("c0")
        terms = [(c("c1"), 1, 0.0), (c("c2"), 2, 180.0), (c("c3"), 3, 0.0)]
    elif form in ("C3", "C4", "C6"):
        nmax = int(form[1])
        const = c("c0")
        for n in range(1, nmax + 1):
            terms.append((c(f"c{n}"), n, 0.0))
            const -= c(f"c{n}")
    elif form == "aro":
        terms = [(c("e0"), 2, 180.0 - float(t["e1"]))]
    elif form == "phosphorus":
        f = float(t["f"])
        const = c("c0")
        terms = [(c("c1"), 1, -f), (c("c2"), 2, 180.0 - 2 * f), (c("c3"), 3, -3 * f)]
    elif form == "EH":
        terms = [(c("c1"), 3, 180.0)]
    elif form == "pfc":   # sum c_n cos^n(psi), psi = phi - 180: cos psi = -cos phi
        b = collections.defaultdict(float)
        for n in range(8):
            for m, w in cos_power(n).items():
                b[m] += c(f"c{n}") * (-1) ** n * w
        const = b.pop(0, 0.0)
        for m, v in sorted(b.items()):
            terms.append((v, m, 0.0))
            const -= v
    else:
        raise ValueError(form)
    terms = [(k, n, d % 360.0) for k, n, d in terms if abs(k) > 1e-12]
    if abs(const) > 1e-12:
        terms.append((const, 0, 90.0))
    return terms


def site_u(form, t, phi):
    """the site's own formula (K), for checking the conversion"""
    g = lambda k: float(t[k])
    r = math.radians
    if form == "opls":
        return g("c0") + g("c1") * (1 + math.cos(r(phi))) + g("c2") * (1 - math.cos(r(2 * phi))) + g("c3") * (1 + math.cos(r(3 * phi)))
    if form in ("C3", "C4", "C6"):
        return g("c0") + sum(g(f"c{n}") * math.cos(r(n * phi)) for n in range(1, int(form[1]) + 1))
    if form == "aro":
        return g("e0") * (1 - math.cos(r(2 * phi) + r(g("e1"))))
    if form == "phosphorus":
        f = g("f")
        return g("c0") + g("c1") * (1 + math.cos(r(phi + f))) + g("c2") * (1 - math.cos(r(2 * (phi + f)))) + g("c3") * (1 + math.cos(r(3 * (phi + f))))
    if form == "EH":
        return g("c1") * (1 - math.cos(r(3 * phi)))
    if form == "pfc":
        x = math.cos(r(phi - 180))
        return sum(g(f"c{n}") * x ** n for n in range(8))


worst = 0.0
for d in mols:
    for T in d["torsion_type"]:
        if T["torsion_type"] == "cistrans":
            continue
        for t in T["torsion"]:
            fs = fourier(T["torsion_type"], t)
            for phi in range(0, 360, 5):
                u = sum(k * (1 + math.cos(math.radians(n * phi - dd))) for k, n, dd in fs)
                worst = max(worst, abs(u - site_u(T["torsion_type"], t, phi) * KB))
print(f"torsion forms as Fourier terms: largest difference from the site's formulas {worst:.1e} kcal/mol")

# ---------------------------------------------------------------- bonded terms from the models' own lists
def key_sym(k):
    return min(tuple(k), tuple(k)[::-1])


obs = {"bond": collections.defaultdict(set), "angle": collections.defaultdict(set), "dihedral": collections.defaultdict(set)}
rigid_models, coul14_models = [], []
seen_in = {"bond": collections.defaultdict(set), "angle": collections.defaultdict(set)}   # key -> {rigid?} of the models it is in
cistrans, cistrans_ref = collections.defaultdict(lambda: [None, None]), {}
for d in mols:
    if d["incomplete"] or not usable(d):
        continue
    ty = {int(a["natom"]): tname(int(a["nonbonded_ID"])) for a in d["nonbonded"]}
    if d["molname"][0]["lrigid"] == "1":
        rigid_models.append(d["name"])
    if any(m and m.get("message_ID") == "3" for m in d.get("messages") or []):
        coul14_models.append(d["name"])
    real = lambda *xs: all(int(x) in d["idx"] for x in xs)
    for s_ in d["stretch"]:
        if real(s_["A"], s_["B"]):
            obs["bond"][key_sym((ty[int(s_["A"])], ty[int(s_["B"])]))].add((float(s_["distance"]), s_["ref_ID"]))
            seen_in["bond"][key_sym((ty[int(s_["A"])], ty[int(s_["B"])]))].add(d["molname"][0]["lrigid"] == "1")
    for b in d["bend"]:
        if not real(b["A"], b["B"], b["C"]):
            continue
        k = (ty[int(b["A"])], ty[int(b["B"])], ty[int(b["C"])])
        obs["angle"][key_sym(k)].add((float(b["theta"]), float(b["ktheta"] or 0), b["ref_ID"]))
        seen_in["angle"][key_sym(k)].add(d["molname"][0]["lrigid"] == "1")
    listed = set()
    for T in d["torsion_type"]:
        for t in T["torsion"]:
            q = (int(t["A"]), int(t["B"]), int(t["C"]), int(t["D"]))
            listed.add(key_sym(q))
            k = tuple(ty[x] for x in q)
            if T["torsion_type"] == "cistrans":   # d0/2 (psi - psi0)^2; cis psi0 180, trans psi0 0 (phi 0 / 180)
                cis = abs(float(t["phi0"]) - 180.0) < 1e-6
                cistrans[k][0 if cis else 1] = float(t["d0"]) / 2 * KB
                cistrans_ref[k] = t["ref_ID"]
                continue
            fs = tuple((round(a, 10), n, round(dd, 8)) for a, n, dd in fourier(T["torsion_type"], t))
            obs["dihedral"][key_sym(k)].add((fs, T["torsion_type"], t["ref_ID"], d["name"]))
    # a torsion the model does not list has none in TraPPE (a linear nitrile's, a diene's second end ...): a zero term
    nbr = collections.defaultdict(set)
    for s_ in d["stretch"]:
        if real(s_["A"], s_["B"]):
            nbr[int(s_["A"])].add(int(s_["B"]))
            nbr[int(s_["B"])].add(int(s_["A"]))
    for j in nbr:
        for k_ in nbr[j]:
            for i in nbr[j] - {k_}:
                for l in nbr[k_] - {j, i}:
                    if key_sym((i, j, k_, l)) not in listed:
                        obs["dihedral"][key_sym((ty[i], ty[j], ty[k_], ty[l]))].add(((), "none", None, d["name"]))
# where models disagree on a torsion the listed one stands, and a model that omits it is named (the database's gap)
lost = collections.defaultdict(list)
for k, v in obs["dihedral"].items():
    sets = {x[0] for x in v}
    if len(sets) > 1:
        keep = [x for x in v if x[0]]
        if len({x[0] for x in keep}) > 1:
            print(f"  torsion {'-'.join(k)}: models disagree: {sorted({(x[0], x[3]) for x in keep})}")
        for x in v:
            if not x[0]:
                lost[x[3]].append("-".join(k))
        obs["dihedral"][k] = {x for x in v if x[0]}
for kind in ("bond", "angle"):
    for k, v in obs[kind].items():
        if len({x[:-1] for x in v}) > 1:
            print(f"  {kind} {'-'.join(k)} has {len(v)} parameter sets in different models: {sorted(v)}")
KBOND = 600.0   # stiff harmonic for TraPPE's fixed lengths
# a term found only in TraPPE's rigid models (CO2, N2, O2 ...) takes K 0: CAPS's mark of a rigid molecule (apply_rigid_geometry,
# sorption / adsorption, LAMMPS fix rigid), as the converted library had it
rigid_only = lambda kind, k: seen_in[kind].get(k) == {True}
BONDS = [{"name": "-".join(k), "match": list(k), "style": "harmonic", "params": [0.0 if rigid_only("bond", k) else KBOND, sorted(v)[0][0]],
          "comment": (f"rigid model: run it as a rigid body" if rigid_only("bond", k) else "fixed length in TraPPE") + f" ({cite(sorted(v)[0][1])})"}
         for k, v in sorted(obs["bond"].items())]
ANGLES = []
for k, v in sorted(obs["angle"].items()):
    th, kt, ref = sorted(v)[0]
    ANGLES.append({"name": "-".join(k), "match": list(k), "style": "harmonic",
                   "params": [kt / 2 * KB if kt > 0 else (0.0 if rigid_only("angle", k) else 100.0), th],
                   "comment": (f"k_theta {kt:g} K/rad^2" if kt > 0 else "rigid model: run it as a rigid body" if rigid_only("angle", k)
                               else "rigid in TraPPE (stiff here; hold the molecule rigid)") + f"; {cite(ref)}"})
DIHEDRALS = []
for k, v in sorted(obs["dihedral"].items()):
    fs, form, ref, _ = sorted(v)[0]
    if not fs:
        DIHEDRALS.append({"name": "-".join(k), "match": list(k), "style": "fourier", "params": [1, 0.0, 1, 0.0],
                          "comment": "no torsion in the TraPPE models that have it"})
        continue
    DIHEDRALS.append({"name": "-".join(k), "match": list(k), "style": "fourier", "params": [len(fs)] + [x for f in fs for x in f],
                      "comment": f"TraPPE {form} form; {cite(ref)}"})
# wildcard rules: a bend centre, or a torsion's central bond, with one parameter set in every model it appears in
# written before the exact rules (the last matching rule wins): an exact rule (a rigid model's K 0) is never overridden
WILD_ANGLES = []
by_centre, centre_rigid = collections.defaultdict(set), collections.defaultdict(set)
for k, v in obs["angle"].items():
    by_centre[k[1]] |= {x[:2] for x in v}
    centre_rigid[k[1]] |= seen_in["angle"].get(k, set())
for c_, v in sorted(by_centre.items()):
    if len(v) == 1:
        th, kt = next(iter(v))
        WILD_ANGLES.append({"name": f"X-{c_}-X", "match": ["*", c_, "*"], "style": "harmonic",
                       "params": [kt / 2 * KB if kt > 0 else (0.0 if centre_rigid[c_] == {True} else 100.0), th],
                       "comment": "every bend at this centre in the TraPPE models (CAPS's wildcard for other neighbours)"})
by_bond = collections.defaultdict(set)
for k, v in obs["dihedral"].items():
    by_bond[key_sym(k[1:3])] |= {x[0] for x in v}   # a zero (no torsion) set counts as its own
for b, v in sorted(by_bond.items()):
    if len(v) == 1:
        fs = next(iter(v)) or ((0.0, 1, 0.0),)
        DIHEDRALS.append({"name": f"X-{b[0]}-{b[1]}-X", "match": ["*", b[0], b[1], "*"], "style": "fourier", "params": [len(fs)] + [x for f in fs for x in f],
                          "comment": "every torsion about this bond in the TraPPE models (CAPS's wildcard for other ends)"})
DEVIATIONS = {"2-hydroxyethyl acrylate": "its CH2-CH2-OH bend next to the ether CH2 is the ether bend (112 deg, 50300 K) where 2-methoxyethan-1-ol, "
                                        "the same chemistry, has the alcohol bend (109.47 deg, 50400 K); CAPS gives the latter",
              "cyclododecane": "its torsion list is numbered for a ten-membered ring (8-9-10-1, 9-10-1-2, 10-1-2-3) and leaves two torsions out",
              "n-butyl methacrylate": "its bend list names 2-3-5 (atom 5 is bonded to 2, not 3) for 2-3-6",
              "ethylbenzene": "its bend list gives 2-3-4 twice and leaves 3-4-5 out (the rigid ring's bends)",
              "naphthalene": "its bend list leaves two ring bends out (a rigid model)",
              "sarin": "its bend list leaves the O-P-F bend out", "soman": "its bend list leaves the O-P-F bend out",
              "dimethyl methylphosphonate": "its bend list leaves the second P-O-CH3 bend (1-4-7) out, and it lists no torsions"}
CORRECTED = {"carbon dioxide": "its bond list gives 2-3 (O-O) for 1-3 (the second C=O); read as 1-3"}
for m_, ks in lost.items():
    DEVIATIONS.setdefault(m_, "its torsion list leaves out " + ", ".join(sorted(set(ks))) + ", which the other TraPPE models with these types list")
IMPROPERS = []
NOTES = ["built by CAPS (bench/ff/make_trappe.py) from the TraPPE database's own models (trappe.oit.umn.edu), each term with its paper",
         "rigid models (" + ", ".join(sorted(set(rigid_models))) + ") have no intramolecular terms in TraPPE: their angles are stiff here; hold them rigid",
         "1-4 Coulomb x 0.5 in TraPPE for " + ", ".join(sorted(set(coul14_models))) + " (coul14_rules); every other 1-4 pair excluded",
         "TraPPE fixes bond lengths: bonds are stiff harmonic (600 kcal/mol/A^2), to be constrained in MD"]
for k, (kc, kt) in sorted(cistrans.items()):
    DIHEDRALS.insert(0, {"name": "-".join(k) + " cis/trans", "match": list(k), "style": "cistrans", "params": [kc, kt],
                         "comment": f"d0 / 2 (psi - psi0)^2 about the isomer the structure has (cis or trans); {cite(cistrans_ref[k])}"})
    obs["dihedral"].pop(key_sym(k), None)
    DIHEDRALS[:] = [x for x in DIHEDRALS if not (x["match"] == list(k) and x["style"] != "cistrans")]
work = tempfile.mkdtemp()
ff_path = os.path.join(OUT, "forcefields", "trappe-ua-umn.json")
ty_path = os.path.join(OUT, "typing", "trappe-ua-umn.typing.json")
write_library(ff_path, ty_path)

if "--check" in sys.argv:
    right = total = 0
    bad = []
    for d in mols:
        if d["incomplete"] or not usable(d):
            continue
        p = os.path.join(work, f"m{d['id']}.mol2")
        mol2(d, p)
        out = p + ".types"
        r = subprocess.run([CAPS, "ff", "type", p, "--ff", ff_path, "-o", out], capture_output=True, text=True)
        if not os.path.exists(out):
            bad.append((d["name"], "failed: " + (r.stderr or r.stdout).strip()[-120:]))
            continue
        got = [x.strip() for x in open(out)]
        want = [tname(int(a["nonbonded_ID"])) for a in sorted(d["nonbonded"], key=lambda a: int(a["natom"])) if a["prtype"] != "M"]
        sites = sorted(tname(int(a["nonbonded_ID"])) for a in d["nonbonded"] if a["prtype"] == "M")
        if sorted(got[len(want):]) != sites:
            bad.append((d["name"], f"charge sites {got[len(want):]} where the model has {sites}"))
        got = got[:len(want)]
        total += len(want)
        diff = [(k + 1, w, g) for k, (w, g) in enumerate(zip(want, got)) if w != g]
        right += len(want) - len(diff)
        if diff:
            bad.append((d["name"], ", ".join(f"{k} {w}->{g}" for k, w, g in diff)))
    print(f"typing: {right}/{total} sites as the TraPPE database types them; {len(bad)} models differ")
    for n_, b in bad:
        print("  ", n_, ":", b)

    # bonded terms: CAPS's GROMACS topology of each model against the model's own list (same atoms, same parameters)
    def itp_terms(path):
        sec, out = None, {"bonds": {}, "angles": {}, "dihedrals": collections.defaultdict(list)}
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
                out["bonds"][key_sym((int(w[0]), int(w[1])))] = float(w[3]) * 10
            elif sec == "angles":
                out["angles"][key_sym((int(w[0]), int(w[1]), int(w[2])))] = (float(w[4]), float(w[5]) / 4.184 / 2)
            elif sec == "dihedrals" and w[4] == "2":
                q = (int(w[0]), int(w[1]), int(w[2]), int(w[3]))
                out["dihedrals"][key_sym(q)].append(("quadratic", round(float(w[6]) / 4.184 / 2, 6), round(float(w[5]) % 360, 4)))
            elif sec == "dihedrals" and w[4] in ("1", "9"):
                q = (int(w[0]), int(w[1]), int(w[2]), int(w[3]))
                out["dihedrals"][key_sym(q)].append((round(float(w[6]) / 4.184, 6), int(w[7]), round(float(w[5]) % 360, 4)))
        return out
    nb_ok = nb_tot = 0
    bad = []
    for d in mols:
        if d["incomplete"] or not usable(d) or d["name"] in DEVIATIONS:
            continue
        p = os.path.join(work, f"m{d['id']}.mol2")
        stem = p + ".gmx"
        r = subprocess.run([CAPS, "ff", "apply", p, "--ff", ff_path, "--charges", "types", "--gromacs", stem], capture_output=True, text=True)
        if not os.path.exists(stem + ".itp"):
            bad.append((d["name"], (r.stderr or r.stdout).strip()[-200:]))
            continue
        got = itp_terms(stem + ".itp")
        half = any(m and m.get("message_ID") == "3" for m in d.get("messages") or [])
        itp = open(stem + ".itp").read()
        top = open(stem + ".top").read()
        fudge = re.search(r"\[\s*defaults\s*\]\s*\n(?:;.*\n)*\s*\S+\s+\S+\s+\S+\s+\S+\s+(\S+)", top)
        # 1-4 Coulomb x 0.5: [ pairs ] function 2 with qq scale 0.5, or every pair by fudgeQQ 0.5
        if half != (bool(re.search(r"^\s*\d+\s+\d+\s+2\s+0\.5\s", itp, re.M)) or (fudge is not None and abs(float(fudge.group(1)) - 0.5) < 1e-9)):
            bad.append((d["name"], f"1-4 Coulomb x 0.5 in TraPPE: {half}, in CAPS: {not half}"))
            continue
        I = d["idx"]
        tyd = {int(a["natom"]): tname(int(a["nonbonded_ID"])) for a in d["nonbonded"]}
        stiff = lambda b: 0.0 if rigid_only("angle", key_sym(tuple(tyd[int(b[x])] for x in "ABC"))) else 100.0   # as the library writes it
        want_a = {key_sym((I[int(b["A"])], I[int(b["B"])], I[int(b["C"])])): (float(b["theta"]), float(b["ktheta"] or 0) / 2 * KB if float(b["ktheta"] or 0) > 0 else stiff(b))
                  for b in d["bend"] if all(int(b[x]) in I for x in "ABC")}
        want_d = {}
        for T in d["torsion_type"]:
            for t in T["torsion"]:
                q = key_sym((I[int(t["A"])], I[int(t["B"])], I[int(t["C"])], I[int(t["D"])]))
                if T["torsion_type"] == "cistrans":
                    want_d[q] = [("quadratic", round(float(t["d0"]) / 2 * KB, 6), round((float(t["phi0"]) + 180) % 360, 4))]
                else:
                    want_d[q] = sorted((round(a, 6), n, round(dd % 360, 4)) for a, n, dd in fourier(T["torsion_type"], t))
        errs = []
        for k, (th, kk) in want_a.items():
            g = got["angles"].get(k)
            if not g or abs(g[0] - th) > 1e-6 or abs(g[1] - kk) > 1e-6 * max(1, kk):
                errs.append(f"bend {k} {g} vs {(th, kk)}")
        for k in got["angles"]:
            if k not in want_a:
                errs.append(f"extra bend {k}")
        for k, v in want_d.items():
            g = sorted(got["dihedrals"].get(k, []))
            if v and v[0][0] == "quadratic":
                if g != v and not (len(g) == 1 and g[0][0] == "quadratic" and abs(g[0][1] - v[0][1]) < 1e-5 and abs(g[0][2] - v[0][2]) < 1e-3):
                    errs.append(f"torsion {k} {g} vs {v}")
            elif len(g) != len(v) or any(abs(a[0] - b[0]) > 1e-5 or a[1] != b[1] or abs(a[2] - b[2]) > 1e-3 for a, b in zip(g, v)):
                errs.append(f"torsion {k} {g} vs {v}")
        for k in got["dihedrals"]:
            if k not in want_d:
                errs.append(f"extra torsion {k}")
        for s_ in d["stretch"]:
            if not (int(s_["A"]) in I and int(s_["B"]) in I):
                continue
            k = key_sym((I[int(s_["A"])], I[int(s_["B"])]))
            if abs(got["bonds"].get(k, -1) - float(s_["distance"])) > 1e-6:
                errs.append(f"bond {k}")
        nb_tot += 1
        if errs:
            bad.append((d["name"], "; ".join(errs[:4]) + (" ..." if len(errs) > 4 else "")))
        else:
            nb_ok += 1
    print(f"bonded terms: {nb_ok}/{nb_tot} models get exactly their own bonds, bends and torsions; the database's own errors listed apart:")
    for m_, why in DEVIATIONS.items():
        print(f"   deviation {m_}: {why}")
    for m_, why in CORRECTED.items():
        print(f"   corrected {m_}: {why}")
    for n_, b in bad:
        print("  ", n_, ":", b)
