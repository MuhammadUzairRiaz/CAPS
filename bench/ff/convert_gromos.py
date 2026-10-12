#!/usr/bin/env python3
"""GROMOS 54A7 for CAPS from the files GROMACS distributes (share/top/gromos54a7.ff, LGPL 2.1).

  convert_gromos.py FF_DIR OUT.json [--typing ../typing/ID.typing.json]

GROMOS is a building-block force field: each residue's topology (aminoacids.rtp) names its atoms' types and charges
and every bonded term by code (gb_N, ga_N, gd_N, gi_N, #defined in ffbonded.itp); the termini (aminoacids.n.tdb /
c.tdb) patch the first and last residue. CAPS keeps the same building blocks and applies them as pdb2gmx does:

  types      ffnonbonded.itp [ atomtypes ] C6 / C12 (the diagonal), atomtypes.atp masses; every other pair from
             [ nonbond_params ] (GROMOS's C12 depends on the pair: a pair-specific entry each); 1-4 from [ pairtypes ],
             which are the geometric means of each type's own 1-4 C6 / C12 (checked here), so per-type lj14 values
  codes      gb_N: V = 1/4 kb (b^2 - b0^2)^2; ga_N: 1/2 k (cos theta - cos theta0)^2; gd_N: k (1 + cos(n phi - phi_s));
             gi_N: 1/2 k (xi - xi0)^2 (kept in GROMACS units: nm, kJ/mol, degrees; the core converts)
  residues   every [ residue ] of the .rtp; each amino acid (a -C and a +N link) also with each N terminus (NH3+, NH2,
             GLY-NH3+, PRO-NH2+, PRO-NH: GLY and PRO take their own) and each C terminus (COO-, COOH), and both: the patch's
             replace / add / delete applied, the residue's terms on deleted or absent atoms dropped (-C at the N end, +N at
             the C end), the patch's own terms in place of the residue's on the same atoms
  1-4        every pair three bonds apart is a 1-4 pair (nrexcl 3, gen-pairs from the bonds), less the residue's [ exclusions ]
"""
import json, math, os, re, sys

D = sys.argv[1]
OUT = sys.argv[2]
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
KJ = 4.184
SIG = lambda c6, c12: (c12 / c6) ** (1 / 6) if c6 > 0 and c12 > 0 else 0.0
EPS = lambda c6, c12: c6 * c6 / (4 * c12) if c6 > 0 and c12 > 0 else 0.0


def sections(path):
    sec = None
    for raw in open(path, errors="replace"):
        l = raw.split(";")[0].strip()
        if not l:
            continue
        m = re.match(r"\[\s*([^\]]+?)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        yield sec, l.split(), l


# ---- non-bonded
mass = {}
for l in open(os.path.join(D, "atomtypes.atp")):
    w = l.split(";")[0].split()
    if len(w) >= 2:
        mass[w[0]] = float(w[1])
atomtypes, nbp, p14 = {}, {}, {}
for sec, w, _ in sections(os.path.join(D, "ffnonbonded.itp")):
    if sec == "atomtypes" and len(w) >= 7:
        atomtypes[w[0]] = (int(w[1]), float(w[5]), float(w[6]), w[4])
    elif sec == "nonbond_params" and len(w) >= 5:
        nbp[tuple(sorted(w[:2]))] = (float(w[3]), float(w[4]))
    elif sec == "pairtypes" and len(w) >= 5:
        p14[tuple(sorted(w[:2]))] = (float(w[3]), float(w[4]))
# per-type 1-4 C6 / C12 from the diagonal; every listed pair must be their geometric mean
c14 = {a: p14[(a, a)] for a in atomtypes if (a, a) in p14}
worst = 0.0
for (a, b), (c6, c12) in p14.items():
    if a in c14 and b in c14:
        g6, g12 = math.sqrt(c14[a][0] * c14[b][0]), math.sqrt(c14[a][1] * c14[b][1])
        worst = max(worst, abs(g6 - c6) / max(c6, 1e-30), abs(g12 - c12) / max(c12, 1e-30))
# chloroform's own pairs (CChl, CLChl, HChl) are not geometric means, but a five-atom molecule has no 1-4 pair
worst = max([0.0] + [v for v in [0.0]])
for (a, b), (c6, c12) in p14.items():
    if a in c14 and b in c14 and not ({a, b} <= {"CChl", "CLChl", "HChl"}):
        g6, g12 = math.sqrt(c14[a][0] * c14[b][0]), math.sqrt(c14[a][1] * c14[b][1])
        worst = max(worst, abs(g6 - c6) / max(c6, 1e-30), abs(g12 - c12) / max(c12, 1e-30))
print(f"pairtypes: {len(p14)}, largest departure from the geometric mean of the diagonal {worst:.1e} (chloroform's own pairs aside)")
if worst > 1e-3:
    raise SystemExit("1-4 pairtypes are not geometric per type")
EL = {0: "X", 1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na", 12: "Mg", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 20: "Ca",
      26: "Fe", 29: "Cu", 30: "Zn", 35: "Br"}
types, pairs = [], []
for name, (z, c6, c12, ptype) in atomtypes.items():
    if ptype != "A":
        continue
    types.append({"name": name, "element": EL.get(z, "X"), "mass": mass.get(name, 0.0)})
    e = {"name": name, "match": [name], "params": [EPS(c6, c12) / KJ, SIG(c6, c12) * 10]}
    if name in c14:
        e["lj14"] = [EPS(*c14[name]) / KJ, SIG(*c14[name]) * 10]
    pairs.append(e)
for (a, b), (c6, c12) in sorted(nbp.items()):
    if a in atomtypes and b in atomtypes and a != b:
        pairs.append({"name": f"{a}-{b}", "match": [a, b], "params": [EPS(c6, c12) / KJ, SIG(c6, c12) * 10]})

# ---- bonded codes
codes = {}
for raw in open(os.path.join(D, "ffbonded.itp")):
    m = re.match(r"#define\s+(g[badi]_\d+)\s+(.*)", raw)
    if m:
        v = [float(x) for x in m.group(2).split(";")[0].split()]
        codes[m.group(1)] = v
# [ bondtypes ] / [ angletypes ] / [ dihedraltypes ] for special bonds (disulfide, heme): kept by type
spec = {"bonds": [], "angles": [], "dihedrals": []}
for sec, w, _ in sections(os.path.join(D, "ffbonded.itp")):
    if sec in ("bondtypes", "angletypes", "dihedraltypes") and w and w[-1] in codes:
        spec[{"bondtypes": "bonds", "angletypes": "angles", "dihedraltypes": "dihedrals"}[sec]].append(w[:-2] + [w[-1]])


# ---- residues
def sections_with_headers(path):
    """as sections(), with each header yielded too (w None)"""
    sec = None
    for raw in open(path, errors="replace"):
        l = raw.split(";")[0].strip()
        if not l:
            continue
        m = re.match(r"\[\s*([^\]]+?)\s*\]", l)
        if m:
            sec = m.group(1)
            yield sec, None, l
            continue
        yield sec, l.split(), l


def read_rtp(path):
    res, cur = {}, None
    for sec, w, l in sections_with_headers(path):
        if w is None:
            if sec not in ("atoms", "bonds", "angles", "dihedrals", "impropers", "exclusions", "bondedtypes"):
                cur = {"name": sec, "atoms": [], "bonds": [], "angles": [], "dihedrals": [], "impropers": [], "exclusions": []}
                res[sec] = cur
            continue
        if sec in ("atoms", "bonds", "angles", "dihedrals", "impropers", "exclusions"):
            if cur is None:
                continue
            if sec == "atoms":
                cur["atoms"].append([w[0], w[1], float(w[2])])
            elif sec == "exclusions":
                cur["exclusions"].append(w[:2])
            else:
                n = {"bonds": 2, "angles": 3, "dihedrals": 4, "impropers": 4}[sec]
                if len(w) > n and w[n] in codes:
                    cur[sec].append(w[:n] + [w[n]])
                elif len(w) > n:
                    raise SystemExit(f"{cur['name']}: {sec} line with an unknown code: {l}")
                else:   # no code: GROMACS's [ ...types ] by atom type (none in 54A7's residues)
                    raise SystemExit(f"{cur['name']}: {sec} line without a code: {l}")
    return res


def read_tdb(path):
    out, cur, sec = {}, None, None
    for raw in open(path, errors="replace"):
        l = raw.split(";")[0].strip()
        if not l:
            continue
        m = re.match(r"\[\s*([^\]]+?)\s*\]", l)
        if m:
            k = m.group(1)
            if k in ("replace", "add", "delete", "bonds", "angles", "dihedrals", "impropers"):
                sec = k
                pending = None
            else:
                cur = {"replace": [], "add": [], "delete": [], "bonds": [], "angles": [], "dihedrals": [], "impropers": []}
                out[k] = cur
                sec = None
            continue
        w = l.split()
        if sec == "replace":
            if not re.match(r"^-?[\d.]+$", w[2]):   # old new type mass charge
                cur["replace"].append({"old": w[0], "new": w[1], "type": w[2], "charge": float(w[4])})
            else:                                   # old type mass charge
                cur["replace"].append({"old": w[0], "new": w[0], "type": w[1], "charge": float(w[3])})
        elif sec == "add":
            if w[0].isdigit():
                cur["add"].append({"n": int(w[0]), "geom": int(w[1]), "name": w[2], "control": w[3:]})
            else:
                cur["add"][-1].update({"type": w[0], "charge": float(w[2])})
        elif sec == "delete":
            cur["delete"].append(w[0])
        elif sec in ("bonds", "angles", "dihedrals", "impropers"):
            n = {"bonds": 2, "angles": 3, "dihedrals": 4, "impropers": 4}[sec]
            cur[sec].append(w[:n] + [w[n]])
    return out


rtp = read_rtp(os.path.join(D, "aminoacids.rtp"))
ntdb = read_tdb(os.path.join(D, "aminoacids.n.tdb"))
ctdb = read_tdb(os.path.join(D, "aminoacids.c.tdb"))
KINDS = ("bonds", "angles", "dihedrals", "impropers")


def patch(r, p, end):
    """residue r with terminus patch p at end 'N' or 'C' (pdb2gmx's tdb semantics)"""
    atoms = [list(a) for a in r["atoms"]]
    rename = {}
    gone = set(p["delete"])
    atoms = [a for a in atoms if a[0] not in gone]
    for a in p["add"]:
        names = [a["name"]] if a["n"] == 1 else [f"{a['name']}{k}" for k in range(1, a["n"] + 1)]
        for nm in names:
            if nm not in [x[0] for x in atoms]:
                atoms.append([nm, a["type"], a["charge"]])
    for rp in p["replace"]:
        hit = [a for a in atoms if a[0] == rp["old"]]
        if rp["old"] != rp["new"] and any(a[0] == rp["new"] for a in atoms) and hit:
            atoms = [a for a in atoms if a[0] != rp["new"]]   # the added atom of the same name: the renamed one stands
        for a in atoms:
            if a[0] == rp["old"]:
                a[0], a[1], a[2] = rp["new"], rp["type"], rp["charge"]
                rename[rp["old"]] = rp["new"]
    names = {a[0] for a in atoms}
    link = "-" if end == "N" else "+"
    out = {"name": r["name"], "atoms": atoms, "exclusions": [[rename.get(x, x) for x in e] for e in r["exclusions"]]}
    for k in KINDS:
        terms = []
        for t in r[k]:
            at = [rename.get(x, x) for x in t[:-1]]
            if any(x[:1] == link for x in at) or any(x[:1] not in "+-" and x not in names for x in at):
                continue
            terms.append(at + [t[-1]])
        own = [list(t) for t in p[k]]
        keys = {frozenset(t[:-1]) for t in own}
        terms = [t for t in terms if frozenset(t[:-1]) not in keys] + own
        out[k] = terms
    out["exclusions"] = [e for e in out["exclusions"] if all(x in names for x in e)]
    return out


def links(r):
    ls = {x for k in KINDS for t in r[k] for x in t[:-1] if x[:1] in "+-"}
    return any(x[:1] == "-" for x in ls), any(x[:1] == "+" for x in ls)


residues = {}
for name, r in rtp.items():
    residues[name] = {k: r[k] for k in ("atoms", "exclusions") + KINDS}
    prev, nxt = links(r)
    if not (prev and nxt and any(a[0] == "CA" for a in r["atoms"]) and any(a[0] == "N" for a in r["atoms"])):
        continue   # not an amino acid in a chain
    nlist = ["PRO-NH2+", "PRO-NH"] if name == "PRO" else (["GLY-NH3+", "NH2"] if name == "GLY" else ["NH3+", "NH2"])
    for nt in nlist:
        residues[f"{name}:{nt}"] = patch(r, ntdb[nt], "N")
    for ct in ("COO-", "COOH"):
        residues[f"{name}::{ct}"] = patch(r, ctdb[ct], "C")
        for nt in nlist:
            residues[f"{name}:{nt}:{ct}"] = patch(patch(r, ntdb[nt], "N"), ctdb[ct], "C")
canon = {t["name"].upper(): t["name"] for t in types}   # GROMACS matches type names whatever their case (CCHL, CChl)
for v in residues.values():
    v.pop("name", None)
    for a in v["atoms"]:
        a[1] = canon.get(a[1].upper(), a[1])

doc = {"format": "caps-forcefield", "format_version": 1, "name": "GROMOS 54A7 (GROMACS)",
       "source": "gromos54a7.ff (GROMACS share/top: ffnonbonded.itp, ffbonded.itp, atomtypes.atp, aminoacids.rtp / .n.tdb / .c.tdb)",
       "references": ["N. Schmid et al., Eur. Biophys. J. 40, 843 (2011)", "C. Oostenbrink et al., J. Comput. Chem. 25, 1656 (2004)"],
       "units": "real", "version": "54A7 (GROMACS gromos54a7.ff)",
       "styles": {"pair": "lj/cut/coul/long", "bond": "gromos", "angle": "cosine/squared", "dihedral": "charmm", "improper": "harmonic"},
       "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 14,
       "atom_types": types, "pairs": pairs, "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
       "residue_templates": {"format": "caps-residue-templates", "units": "GROMACS (nm, kJ/mol, degrees)",
                             "codes": codes, "special": spec, "residues": residues,
                             "nrexcl": 3, "termini": {"n": sorted(ntdb), "c": sorted(ctdb)}},
       "notes": ["converted by CAPS (bench/ff/convert_gromos.py) from GROMACS's gromos54a7.ff (LGPL 2.1)",
                 "a building-block force field: residues CAPS recognises take their topology from aminoacids.rtp and the termini patches, "
                 "as pdb2gmx applies them; the Lennard-Jones of every type pair from nonbond_params, 1-4 from pairtypes (geometric per type)"]}
if arg("--typing"):
    doc["typing"] = arg("--typing")
json.dump(doc, open(OUT, "w"), separators=(",", ":"), ensure_ascii=False)   # compact: 420 residue variants
print(f"{OUT}: {len(types)} types, {len(pairs)} pair entries, {len(codes)} codes, {len(rtp)} residues ({len(residues)} with termini)")
