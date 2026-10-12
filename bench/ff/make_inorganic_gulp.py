#!/usr/bin/env python3
"""Ionic potential sets for CAPS from their papers, as the UCL Database of Published Interatomic Potential Parameters
(mcc.chem.ucl.ac.uk/Potentials) gives each paper's set in GULP's library format; one CAPS library per paper.

  make_inorganic_gulp.py LIB_DIR [ID ...] [--out DATA_DIR]

GULP's units and forms (eV, A) to CAPS's (real units, kcal/mol):
  species     each label (Mg2, O1 ...) core and, with a shell, shell: types LABELc / LABELs, the label alone without one;
              the charges as the paper gives them; a shell carries 10 % of the ion's mass (adiabatic core-shell, as
              LAMMPS's core/shell examples), the core the rest
  buckingham  A exp(-r/rho) - C/r^6 between the two species named: CAPS buck pair (A, C x 23.060548), between their
              types only; lennard A/r^12 - B/r^6 as Lennard-Jones (eps = B^2/4A, sigma = (A/B)^(1/6)); morse as Morse
  spring      1/2 k2 r^2 + 1/24 k4 r^4 between a core and its shell: class II bond K2 = k2/2, K4 = k4/24 (r0 = 0)
  sets        the compounds a paper covers, from which labels its pair terms join: each label is typed in a structure
              holding the elements of its partners and none of the paper's other elements; where one composition has
              several sets (Lewis & Catlow's two MgO models), the first is the default and each other an overlay
              (data/typing/sets/ID/SET.typing.json); a paper whose blocks repeat a pair with other values (Catlow,
              Diller & Norgett's alkali halides: a C per compound) is split block by block into its compounds
Every short-range term acts on the species the paper names; Coulomb between all, a core and its own shell excluded.
"""
import collections, json, os, re, sys

LIB = sys.argv[1]
OUT = sys.argv[sys.argv.index("--out") + 1] if "--out" in sys.argv else os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data")
EV = 23.060548   # kcal/mol per eV
MASS = {"H": 1.008, "Li": 6.94, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Na": 22.99, "Mg": 24.305, "Al": 26.982,
        "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942,
        "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.63,
        "Br": 79.904, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Ag": 107.868, "Sn": 118.71, "I": 126.904,
        "Cs": 132.905, "Ba": 137.327, "La": 138.905, "Ce": 140.116, "Pr": 140.908, "Nd": 144.242, "Eu": 151.964, "Gd": 157.25,
        "Tb": 158.925, "Yb": 173.045, "Ta": 180.948, "Pb": 207.2, "Bi": 208.98, "Th": 232.038, "U": 238.029}

# the papers: id -> (file, name, citation, what it covers)
# blocks of a file that come from other works (diller.lib: AgCl from Catlow, Corish & Jacobs 1979, Na2S unpublished)
NOT_THIS_PAPER = {"catlow1977-halides": {"AgCl", "Na2S", "NaS"}}
PAPERS = {
    "catlow1977-halides": ("diller.lib", "Alkali halides (Catlow, Diller & Norgett 1977)",
                           "C. R. A. Catlow, K. M. Diller, M. J. Norgett, J. Phys. C: Solid State Phys. 10, 1395 (1977)"),
    "bush1994-oxides": ("bush.lib", "Binary oxides (Bush, Gale, Catlow & Battle 1994)",
                        "T. S. Bush, J. D. Gale, C. R. A. Catlow, P. D. Battle, J. Mater. Chem. 4, 831 (1994)"),
    "lewis-catlow1985-oxides": ("lewis_b.lib", "Binary oxides (Lewis & Catlow 1985)",
                                "G. V. Lewis, C. R. A. Catlow, J. Phys. C: Solid State Phys. 18, 1149 (1985)"),
    "matsui-akaogi1991-tio2": ("matsui-akaogi.lib", "TiO2 polymorphs (Matsui & Akaogi 1991)",
                               "M. Matsui, M. Akaogi, Mol. Simul. 6, 239 (1991)"),
    "cherry1995-perovskites": ("cherry.lib", "Perovskite oxides (Cherry, Islam & Catlow 1995)",
                               "M. Cherry, M. S. Islam, C. R. A. Catlow, J. Solid State Chem. 118, 125 (1995)"),
}


def element(label):
    m = re.match(r"([A-Z][a-z]?)", label)
    e = m.group(1)
    return e if e in MASS else label[0]


def parse(path):
    """species {(label, kind): charge}, terms [(section, keywords, [(label, kind) ...], params, block)]"""
    species, terms = {}, []
    sec, block = None, 0
    NEEDS = {"buckingham": 2, "lennard": 2, "morse": 2, "spring": 1, "three": 3}
    for raw in open(path, errors="replace"):
        line = raw.split("#")[0].rstrip()
        w = line.split()
        if not w:
            block += 1
            continue
        if w[0].lower() in ("species", "buckingham", "lennard", "morse", "spring", "three", "cutoffs", "polarisability", "torsion"):
            sec = w[0].lower()
            block += 1
            continue
        if sec == "species" and len(w) >= 3:
            species[(w[0], w[1][:4])] = float(w[2])
        elif sec in NEEDS:
            n = NEEDS[sec]
            ats, k = [], 0
            if sec == "spring":
                terms.append((sec, [], [(w[0], "core"), (w[0], "shel")], [float(x) for x in w[1:]], block))
                continue
            while len(ats) < n:
                lab = w[k]
                kind = w[k + 1][:4] if k + 1 < len(w) and w[k + 1][:4] in ("core", "shel", "bshe") else "core"
                ats.append((lab, kind))
                k += 2 if kind and k + 1 < len(w) and w[k + 1][:4] in ("core", "shel", "bshe") else 1
            nums = [float(x) for x in w[k:] if re.match(r"^-?[\d.]+(e-?\d+)?$", x, re.I)]
            terms.append((sec, [x for x in w[k:] if not re.match(r"^-?[\d.]", x)], ats, nums, block))
    return species, terms


def build(pid):
    fname, name, ref = PAPERS[pid]
    species, terms = parse(os.path.join(LIB, fname))
    # a paper that gives one pair different values in different blocks: one set per block (its compound)
    seen, per_block = {}, False
    for sec, kw, ats, p, blk in terms:
        if sec in ("buckingham", "lennard", "morse"):
            key = (sec, tuple(sorted(ats)))
            if key in seen and seen[key] != p:
                per_block = True
            seen.setdefault(key, p)
    if per_block:
        blocks = collections.defaultdict(list)
        for t in terms:
            if t[0] in ("buckingham", "lennard", "morse"):
                blocks[t[4]].append(t)
        sets = []
        for blk, ts in sorted(blocks.items()):
            labs = sorted({a[0] for t in ts for a in t[2]})
            q = lambda l: species.get((l, "core"), 0) + species.get((l, "shel"), 0)
            comp = "".join(sorted({element(l) for l in labs}, key=lambda e: (min(q(l) for l in labs if element(l) == e) < 0, e)))   # cations first: LiF, Na2S
            sets.append((comp + f"_{len([s for s in sets if s[0].startswith(comp)]) + 1}" if any(s[0].startswith(comp) for s in sets) else comp, labs, ts))
    else:
        sets = None
    return species, terms, sets, name, ref, fname


def write(pid):
    species, terms, sets, name, ref, fname = build(pid)
    types, pairs, bonds, rules, shells, overlays = {}, [], [], [], {}, {}
    all_elements = sorted({element(l) for l, _ in species})
    has_shell = {l for (l, k) in species if k == "shel"}

    def tn(label, kind, suffix=""):
        if suffix:   # a compound's own species: Li_LiF_c, Cl_NaCl_s
            return f"{label}{suffix}_{'c' if kind == 'core' else 's'}" if label in has_shell else f"{label}{suffix}"
        return f"{label}{'c' if kind == 'core' else 's'}" if label in has_shell else label

    def add_type(label, kind, suffix=""):
        t = tn(label, kind, suffix)
        if t in types:
            return t
        el = element(label)
        q = species.get((label, kind), 0.0)
        mass = MASS[el] * (0.1 if kind == "shel" else 0.9 if label in has_shell else 1.0)
        types[t] = {"name": t, "element": el, "mass": round(mass, 6), "charge": q, "description": f"{label} {'shell' if kind == 'shel' else 'core'}; {ref}"}
        if kind == "core" and label in has_shell:
            shells[t] = tn(label, "shel", suffix)
        return t

    def add_terms(ts, suffix=""):
        for sec, kw, ats, p, blk in ts:
            a, b = (add_type(l, k, suffix) for l, k in ats)
            if sec == "buckingham":
                A, rho, C = p[0], p[1], p[2] if len(p) > 2 else 0.0
                pairs.append({"name": f"{a}-{b}", "match": [a, b], "style": "buck", "params": [A * EV, rho, C * EV],
                              "comment": f"A {A:g} eV, rho {rho:g} A, C {C:g} eV A^6; cut-off {p[-1]:g} A in the paper's set"})
            elif sec == "lennard":
                A, B = p[0] * EV, p[1] * EV
                pairs.append({"name": f"{a}-{b}", "match": [a, b], "params": [B * B / (4 * A) if A else 0.0, (A / B) ** (1 / 6) if A and B else 0.0],
                              "comment": f"lennard A {p[0]:g} eV A^12, B {p[1]:g} eV A^6"})
            elif sec == "morse":
                pairs.append({"name": f"{a}-{b}", "match": [a, b], "style": "morse", "params": [p[0] * EV, p[1], p[2]],
                              "comment": f"Morse D {p[0]:g} eV, a {p[1]:g} 1/A, r0 {p[2]:g} A"})

    springs = [t for t in terms if t[0] == "spring"]
    if sets is None:
        add_terms([t for t in terms if t[0] in ("buckingham", "lennard", "morse")])
        for sec, kw, ats, p, blk in springs:
            c, s = tn(ats[0][0], "core"), tn(ats[0][0], "shel")
            add_type(ats[0][0], "core"), add_type(ats[0][0], "shel")
            bonds.append({"name": f"{c}-{s}", "match": [c, s], "style": "class2", "params": [0.0, p[0] / 2 * EV, 0.0, (p[1] / 24 * EV) if len(p) > 1 else 0.0],
                          "comment": f"core-shell spring k2 {p[0]:g} eV/A^2" + (f", k4 {p[1]:g} eV/A^4" if len(p) > 1 else "")})
        # compounds: each cation label with the anion labels it is joined to
        graph = collections.defaultdict(set)
        for sec, kw, ats, p, blk in terms:
            if sec in ("buckingham", "lennard", "morse") and element(ats[0][0]) != element(ats[1][0]):
                graph[ats[0][0]].add(ats[1][0])
                graph[ats[1][0]].add(ats[0][0])
        by_el = collections.defaultdict(set)
        for l, k in species:
            by_el[element(l)].add(l)
        if all(len(v) == 1 for v in by_el.values()):
            # one species per element (Bush et al.: one oxygen for every cation, made for mixed oxides as well): typed by
            # element alone
            for e, labs in sorted(by_el.items()):
                lab = next(iter(labs))
                rules.append({"type": tn(lab, "core"), "smarts": f"[#{ELZ[e]}]", "priority": 0, "description": f"{lab} ({e})"})
            graph = {}
        comps = collections.OrderedDict()   # composition -> list of label sets
        for lab in sorted(graph, key=lambda l: (species.get((l, "core"), 0) + species.get((l, "shel"), 0) < 0, l)):
            partners = graph[lab]
            labs = {lab} | partners
            comp = tuple(sorted({element(l) for l in labs}))
            comps.setdefault(comp, [])
            if labs not in comps[comp]:
                comps[comp].append(labs)
        for comp, alts in comps.items():
            for k, labs in enumerate(alts):
                rs = []
                for lab in sorted(labs):
                    rs.append({"type": tn(lab, "core"), "smarts": f"[#{ELZ[element(lab)]}]", "priority": len(comp) + (50 if k else 0),
                               "requires": [e for e in comp if e != element(lab)], "excludes": [e for e in all_elements if e not in comp],
                               "description": f"{lab} in {''.join(comp)} ({ref.split(',')[0].split('.')[-1].strip()} et al.)" if False else f"{lab} in {''.join(comp)}"})
                if k == 0:
                    rules += rs
                else:
                    overlays[f"{''.join(comp)}_{k + 1}"] = rs
    else:
        for comp, labs, ts in sets:
            if comp in NOT_THIS_PAPER.get(pid, ()):
                continue
            sfx = "_" + comp
            add_terms(ts, sfx)
            for lab in labs:
                for sec, kw, ats, p, blk in springs:
                    if ats[0][0] == lab:
                        c, s = tn(lab, "core", sfx), tn(lab, "shel", sfx)
                        add_type(lab, "core", sfx), add_type(lab, "shel", sfx)
                        if not any(b["name"] == f"{c}-{s}" for b in bonds):
                            bonds.append({"name": f"{c}-{s}", "match": [c, s], "style": "class2", "params": [0.0, p[0] / 2 * EV, 0.0, (p[1] / 24 * EV) if len(p) > 1 else 0.0],
                                          "comment": f"core-shell spring k2 {p[0]:g} eV/A^2"})
            els = sorted({element(l) for l in labs})
            for lab in labs:
                rules.append({"type": tn(lab, "core", sfx), "smarts": f"[#{ELZ[element(lab)]}]", "priority": len(els),
                              "requires": [e for e in els if e != element(lab)], "excludes": [e for e in all_elements if e not in els],
                              "description": f"{lab} in {comp}"})
    # every type a pair entry of its own (zero: the short-range terms are the paper's pairs only)
    own = [{"name": t, "match": [t], "params": [0.0, 0.0]} for t in types]
    doc = {"format": "caps-forcefield", "format_version": 1, "name": name, "source": f"{fname} (UCL Database of Published Interatomic Potential Parameters)",
           "references": [ref], "units": "real",
           "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
           "mixing": "geometric", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10, "equivalence": "replace",
           "atom_types": list(types.values()), "pairs": own + pairs, "bonds": bonds, "angles": [], "dihedrals": [], "impropers": [],
           "notes": [f"built by CAPS (bench/ff/make_inorganic_gulp.py) from {ref}, as the UCL potentials database gives the set ({fname})",
                     "shells carry 10 % of the ion's mass (adiabatic core-shell); the source's shells are massless"],
           "version": ref.split("(")[-1].rstrip(")"), "typing": f"../typing/{pid}.typing.json"}
    json.dump(doc, open(os.path.join(OUT, "forcefields", pid + ".json"), "w"), indent=1, ensure_ascii=False)
    ty = {"format": "caps-typing", "version": 1, "forcefield": name, "unknown_types": "untyped", "bonds": "defined",
          "description": f"CAPS rules for {name}: each species typed by the composition of the compound it belongs to (bench/ff/make_inorganic_gulp.py)",
          "rules": rules}
    if shells:
        ty["shells"] = shells
        # the shells CAPS adds carry their type as their name
        ty["rules"] += [{"type": sh, "smarts": f"[#{ELZ[types[sh]['element']]}]", "priority": 99, "atom_name": sh, "description": f"shell of {c}"}
                        for c, sh in shells.items()]
    json.dump(ty, open(os.path.join(OUT, "typing", pid + ".typing.json"), "w"), indent=1, ensure_ascii=False)
    for oname, rs in overlays.items():
        os.makedirs(os.path.join(OUT, "typing", "sets", pid), exist_ok=True)
        json.dump({"format": "caps-typing", "version": 1, "forcefield": f"{pid}: {oname}", "unknown_types": "untyped",
                   "description": f"Overlay on the {pid} rules: the {oname} set in place of the composition's first", "rules": rs},
                  open(os.path.join(OUT, "typing", "sets", pid, oname + ".typing.json"), "w"), indent=1)
    print(f"{pid}: {len(types)} types, {len(pairs)} pair terms, {len(bonds)} springs, {len(rules)} rules, {len(overlays)} overlay sets")


ELZ = {e: z for z, e in enumerate(["X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
                                   "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
                                   "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce",
                                   "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",
                                   "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"])}

if __name__ == "__main__":
    # checked against an independent transcription of the same paper (the others wait for their papers: Catlow, Diller
    # & Norgett's shell charges and springs, which the file itself says may not be the published ones; Lewis &
    # Catlow's file leaves the O-O C6 as a placeholder, Matsui & Akaogi's units, Cherry's four unconfirmed terms)
    for pid in [a for a in sys.argv[2:] if a in PAPERS] or ["bush1994-oxides"]:
        write(pid)
