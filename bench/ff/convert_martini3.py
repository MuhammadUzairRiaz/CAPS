#!/usr/bin/env python3
"""Martini 3 for CAPS, from the cgmartini.nl downloads (martini_v300.zip: martini_v3.0.0.itp and its solvent, ion, small
molecule, sugar, nucleobase and phospholipid files):

  data/martini/martini3-nonbonded.tsv.gz   every bead type with its mass ("T name mass") and every pair as the file gives it
                                           ("P a b sigma_nm epsilon_kJ", GROMACS sigma-epsilon form): 844 types, the whole
                                           matrix; CAPS reads the pairs of the types a structure uses
  data/martini/martini3-molecules.json     each molecule of the parameter files with its explicit topology (GROMACS units:
                                           nm, kJ/mol, degrees, functions as in the .itp); "#ifdef FLEXIBLE" blocks are
                                           taken (stiff bonds where the default has constraints: CAPS has no constraints)
  data/forcefields/martini3.json           the force field: bead types, the pair table, Martini 3's run settings
  data/typing/martini3.typing.json         beads typed by name; molecules built from the templates carry their topology

Run settings (Martini 3's recommended GROMACS parameters, Souza et al. 2021 and cgmartini.nl): Lennard-Jones cut at
1.1 nm and shifted to zero there (Potential-shift, no dispersion correction), reaction-field Coulomb with epsilon_r 15
and epsilon_rf infinite (0), cut at 1.1 nm; special_bonds 0 1 1 (nrexcl 1).

usage: convert_martini3.py [MARTINI_V300_DIR]   (default ~/vermouth-ref/cgmartini/m3/v300/martini_v300)
"""
import gzip, json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/vermouth-ref/cgmartini/m3/v300/martini_v300")
MAIN = os.path.join(SRC, "martini_v3.0.0.itp")
MOLFILES = ["martini_v3.0.0_solvents_v1.itp", "martini_v3.0.0_ions_v1.itp", "martini_v3.0.0_small_molecules_v1.itp",
            "martini_v3.0.0_sugars_v1.itp", "martini_v3.0.0_nucleobases_v1.itp", "martini_v3.0.0_phospholipids_v1.itp"]
DEFINES = {"FLEXIBLE"}


def lines_with_defines(path):
    """the file's lines with #ifdef / #ifndef / #else / #endif applied for DEFINES"""
    stack, out = [], []
    for raw in open(path, encoding="utf-8", errors="replace"):
        l = raw.strip()
        if l.startswith("#ifdef"):
            stack.append(l.split()[1] in DEFINES)
            continue
        if l.startswith("#ifndef"):
            stack.append(l.split()[1] not in DEFINES)
            continue
        if l.startswith("#else"):
            stack[-1] = not stack[-1]
            continue
        if l.startswith("#endif"):
            stack.pop()
            continue
        if l.startswith("#"):
            continue
        if all(stack):
            out.append(raw)
    return out


def nonbonded():
    sec, types, pairs = None, [], []
    for raw in lines_with_defines(MAIN):
        l = raw.split(";")[0].strip()
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        w = l.split()
        if sec == "atomtypes" and len(w) >= 2:
            types.append((w[0], float(w[1])))
        elif sec == "nonbond_params" and len(w) >= 5 and w[2] == "1":
            pairs.append((w[0], w[1], w[3], w[4]))
    return types, pairs


def molecules(path):
    """{name: molecule} from an .itp: atoms (name, type, charge, mass or None), and the terms by 0-based atom index"""
    out, cur, sec = {}, None, None
    insane = None
    for raw in lines_with_defines(path):
        # a lipid's building blocks, as insane writes them: head, linker, one letter per tail bead (C saturated, D
        # unsaturated); CAPS maps all-atom lipids with them
        m = re.match(r"^;@INSANE\s+(.*)", raw.strip())
        if m:
            f = dict(re.findall(r"(\w+)=([^,]*)", m.group(1)))
            insane = {"head": f.get("alhead", "").split(), "link": f.get("allink", "").split(), "tails": f.get("altail", "").split()}
            continue
        l = raw.split(";")[0].strip()
        if not l:
            continue
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1)
            continue
        w = l.split()
        if sec == "moleculetype":
            cur = {"atoms": [], "bonds": [], "constraints": [], "angles": [], "dihedrals": [], "exclusions": [], "vsites": [],
                   "nrexcl": int(w[1]) if len(w) > 1 else 1, "source": os.path.basename(path)}
            if insane:
                cur["insane"] = insane
                insane = None
            out[w[0]] = cur
            sec = None
            continue
        if cur is None:
            continue
        ix = lambda t: int(t) - 1
        if sec == "atoms":
            cur["atoms"].append({"type": w[1], "name": w[4], "charge": float(w[6]) if len(w) > 6 else 0.0,
                                 "mass": float(w[7]) if len(w) > 7 else None})
        elif sec == "bonds" and len(w) >= 5:
            cur["bonds"].append([ix(w[0]), ix(w[1]), int(w[2]), float(w[3]), float(w[4])])
        elif sec == "constraints" and len(w) >= 4:
            cur["constraints"].append([ix(w[0]), ix(w[1]), int(w[2]), float(w[3])])
        elif sec == "angles" and len(w) >= 6:
            cur["angles"].append([ix(w[0]), ix(w[1]), ix(w[2]), int(w[3]), float(w[4]), float(w[5])])
        elif sec == "dihedrals" and len(w) >= 7:
            cur["dihedrals"].append([ix(w[0]), ix(w[1]), ix(w[2]), ix(w[3]), int(w[4])] + [float(x) for x in w[5:]])
        elif sec == "exclusions" and len(w) >= 2:
            cur["exclusions"].append([ix(x) for x in w])
        elif sec == "virtual_sitesn" and len(w) >= 3:   # 1 centre of geometry, 2 centre of mass, 3 weighted (atom weight ...)
            fn = int(w[1])
            if fn == 3:
                at, wt = [ix(x) for x in w[2::2]], [float(x) for x in w[3::2]]
                cur["vsites"].append({"site": ix(w[0]), "from": at, "weights": [x / sum(wt) for x in wt]})
            else:
                at = [ix(x) for x in w[2:]]
                cur["vsites"].append({"site": ix(w[0]), "from": at, "weights": [1 / len(at)] * len(at) if fn == 1 else None})
        elif sec == "virtual_sites2" and len(w) >= 5:     # site i j 1 a: (1 − a) x_i + a x_j
            a = float(w[4])
            cur["vsites"].append({"site": ix(w[0]), "from": [ix(w[1]), ix(w[2])], "weights": [1 - a, a]})
        elif sec == "virtual_sites3" and len(w) >= 7:
            if int(w[4]) != 1:   # 3fd, 3fad, 3out: non-linear constructions CAPS does not have (nucleobases)
                cur["unsupported"] = f"virtual_sites3 function {w[4]}"
                continue
            a, b = float(w[5]), float(w[6])   # (1 − a − b) x_i + a x_j + b x_k
            cur["vsites"].append({"site": ix(w[0]), "from": [ix(w[1]), ix(w[2]), ix(w[3])], "weights": [1 - a - b, a, b]})
        elif sec in ("settles", "position_restraints"):
            continue
        elif sec not in (None, "bonds", "constraints", "angles", "dihedrals", "exclusions", "virtual_sitesn", "virtual_sites2", "virtual_sites3"):
            raise SystemExit(f"{path}: section [{sec}] is not handled")
    return out


if __name__ == "__main__":
    types, pairs = nonbonded()
    names = [t for t, _ in types]
    have = {tuple(sorted((a, b))) for a, b, _, _ in pairs}
    missing = sum(1 for i, a in enumerate(names) for b in names[i:] if tuple(sorted((a, b))) not in have)
    dd = os.path.join(ROOT, "data", "martini")
    os.makedirs(dd, exist_ok=True)
    table = os.path.join(dd, "martini3-nonbonded.tsv.gz")
    # mtime 0: the same table gives the same bytes
    with gzip.GzipFile(table, "wb", compresslevel=9, mtime=0) as raw, __import__("io").TextIOWrapper(raw, encoding="utf-8") as f:
        f.write("# Martini 3.0.0 non-bonded parameters (martini_v3.0.0.itp, cgmartini.nl): T name mass(g/mol); P a b sigma(nm) epsilon(kJ/mol)\n")
        for t, m in types:
            f.write(f"T {t} {m:g}\n")
        for a, b, sg, ep in pairs:
            f.write(f"P {a} {b} {sg} {ep}\n")
    mols = {}
    for fn in MOLFILES:
        for k, v in molecules(os.path.join(SRC, fn)).items():
            if k in mols:
                print(f"note: {k} in {fn} again (the later one kept)")
            mols[k] = v
    # nrexcl above 1 (GMY, RIB): the pairs within nrexcl bonds (bonds and constraints, as GROMACS counts them) written out as
    # exclusions, so every molecule reads with nrexcl 1 (special_bonds 0 1 1)
    for k, v in mols.items():
        if v["nrexcl"] <= 1:
            continue
        n = len(v["atoms"])
        nb = [set() for _ in range(n)]
        for b in v["bonds"] + v["constraints"]:
            nb[b[0]].add(b[1]); nb[b[1]].add(b[0])
        for i in range(n):
            seen, front = {i}, {i}
            for _ in range(v["nrexcl"]):
                front = {y for x in front for y in nb[x]} - seen
                seen |= front
            far = sorted(j for j in seen if j > i and j not in nb[i])
            if far:
                v["exclusions"].append([i] + far)
        v["nrexcl"] = 1
    for k in [k for k, v in mols.items() if v.get("unsupported")]:
        print(f"left out {k}: {mols[k]['unsupported']} (a non-linear virtual site)")
        del mols[k]
    # molecules whose beads are all in the particle file
    bad = [k for k, v in mols.items() if any(a["type"] not in set(names) for a in v["atoms"])]
    for k in bad:
        print(f"left out {k}: bead types not in martini_v3.0.0.itp")
        del mols[k]
    doc = {"format": "caps-martini-molecules", "version": 1, "model": "Martini 3.0.0",
           "source": "cgmartini.nl martini_v300.zip: " + ", ".join(MOLFILES) + " (FLEXIBLE variants)",
           "units": "GROMACS: nm, kJ/mol, degrees; functions as in the .itp", "molecules": mols}
    json.dump(doc, open(os.path.join(dd, "martini3-molecules.json"), "w"), ensure_ascii=False, separators=(",", ":"))
    ff = {"format": "caps-forcefield", "format_version": 1, "name": "Martini 3", "version": "3.0.0",
          "source": "martini_v3.0.0.itp and its molecule files (cgmartini.nl)",
          "references": ["P. C. T. Souza et al., Nat. Methods 18, 382 (2021)"],
          "units": "real", "styles": {"pair": "lj/cut/coul/reaction-field", "bond": "harmonic", "angle": "cosine/squared", "dihedral": "fourier",
                                      "improper": "harmonic"},
          "mixing": "arithmetic", "special_lj": [0, 1, 1], "special_coul": [0, 1, 1], "cutoff": 11.0,
          "pair_settings": {"coulomb": "reaction-field", "eps_rf": 0, "dielectric": 15.0, "lj_modifier": "potential-shift", "model_cutoff": True},
          "pair_table": "../martini/martini3-nonbonded.tsv.gz", "molecule_templates": "../martini/martini3-molecules.json",
          "torsion_terms": "if_defined", "angle_terms": "if_defined", "typing": "../typing/martini3.typing.json",
          "atom_types": [{"name": t, "element": "", "mass": m, "description": f"Martini 3 bead {t}"} for t, m in types],
          "pairs": [], "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
          "notes": [f"{len(types)} bead types and {len(pairs)} pairs from martini_v3.0.0.itp (the pair table)",
                    "run settings: LJ cut at 11 A and shifted (Potential-shift), reaction field epsilon_r 15, epsilon_rf infinite, special_bonds 0 1 1",
                    f"{len(mols)} molecules from the parameter files as templates with their topologies (FLEXIBLE variants: stiff bonds for constraints)"]}
    json.dump(ff, open(os.path.join(ROOT, "data", "forcefields", "martini3.json"), "w"), ensure_ascii=False, indent=1)
    rules = [{"type": t, "smarts": "*", "atom_name": t, "priority": 0} for t in names]
    # all-atom proteins become Martini 3 beads as martinize2 makes them (data/martini/martini3-protein.json, from
    # vermouth's martini3001 files by bench/ff/convert_vermouth_martini3.py)
    tdoc = {"format": "caps-typing", "version": 1, "forcefield": "Martini 3", "coarse_grained": True, "unknown_types": "untyped",
            "martini_protein": "../martini/martini3-protein.json",
            "martini_small_molecules": "../martini/martini3-small-molecules.json",
            "description": "Martini 3 beads typed by name (bench/ff/convert_martini3.py); all-atom proteins mapped as martinize2 maps them, small molecules by vermouth's mappings",
            "rules": rules}
    json.dump(tdoc, open(os.path.join(ROOT, "data", "typing", "martini3.typing.json"), "w"), ensure_ascii=False, separators=(",", ":"))
    cat_p = os.path.join(ROOT, "data", "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    entry = {"id": "martini3", "name": "Martini 3", "version": "3.0.0", "references": ff["references"], "origin": "cgmartini.nl",
             "source_file": "martini_v3.0.0.itp", "status": "converted",
             "notes": f"Martini 3: {len(types)} bead types, the full pair table, {len(mols)} molecule templates; reaction field (GROMACS export)",
             "file": "martini3.json", "counts": {"atom_types": len(types), "pairs": len(pairs), "bonds": 0, "angles": 0, "dihedrals": 0, "impropers": 0},
             "typing": {"rules": "typing/martini3.typing.json", "evidence": "beads typed by name; molecules from the Martini 3 parameter files"}}
    ids = [e["id"] for e in cat["forcefields"]]
    if "martini3" in ids:
        cat["forcefields"][ids.index("martini3")] = entry
    else:
        cat["forcefields"].insert(ids.index("martini22-proteins") + 1, entry)
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    print(f"{table}: {len(types)} types, {len(pairs)} pairs ({missing} type pairs without an entry), {os.path.getsize(table)} bytes")
    print(f"{len(mols)} molecules")
