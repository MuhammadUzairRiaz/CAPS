#!/usr/bin/env python3
"""Martini 3 small molecules for CAPS, from vermouth-martinize (github.com/marrink-lab/vermouth-martinize, Apache-2.0):
force_fields/martini3001/25-small_molecule_martini3.ff (the coarse-grained blocks), mappings/martini3001/*.charmm36.map
(which atoms make each bead, with vermouth's fractional weights) and force_fields/charmm/10-small_molecule_charmm.ff
(the all-atom CHARMM residues the maps are written for: their atoms' elements and bonds).

CAPS maps an all-atom molecule by graph: its atoms and bonds are matched to a CHARMM residue's (elements and bonds, the
names do not matter), and the matched atoms go into the beads the map puts them in. Written to
data/martini/martini3-small-molecules.json:
  molecules  per molecule: "aa" (the CHARMM residue: atom names, elements, bonds by index; lone pairs left out),
             "map" (per all-atom atom index, [[bead, weight], ...]; atoms the map does not list are in no bead),
             "block" (the Martini 3 block: beads with types, charges, masses; interactions as in the .ff)

The blocks are checked here against cgmartini's martini_v3.0.0_small_molecules_v1.itp (the library's templates,
data/martini/martini3-molecules.json): the same bead types, charges and bonded terms from two sources.

usage: convert_vermouth_martini3_small.py [VERMOUTH_REF_DIR]   (default ~/vermouth-ref, the downloaded files)
"""
import glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import convert_vermouth_martini3 as vm   # the .ff reader and .map reader

SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/vermouth-ref")


def element(ctype):
    t = ctype.upper()
    for pre, el in (("CL", "Cl"), ("BR", "Br"), ("LP", None), ("H", "H"), ("C", "C"), ("N", "N"), ("O", "O"), ("S", "S"),
                    ("P", "P"), ("F", "F"), ("I", "I")):
        if t.startswith(pre):
            return el
    raise SystemExit(f"no element for CHARMM type {ctype}")


NREF = {"bonds": 2, "constraints": 2, "angles": 3, "dihedrals": 4, "impropers": 4, "virtual_sites2": 3, "virtual_sites3": 4}


def zero_based_fixed(path, dst):
    """25-small_molecule_martini3.ff numbers the atoms of its blocks' terms from 0 ("0 1 1 0.330"), while vermouth's
    reader counts from 1 (and so takes "0" for the last atom); a block whose terms use atom 0 is renumbered from 1 here.
    Read this way, the 8 molecules also in cgmartini's martini_v3.0.0_small_molecules_v1.itp are identical to it."""
    lines = open(path).read().splitlines()
    blocks, cur = [], None
    for k, l in enumerate(lines):
        if l.strip().lower().startswith("[ moleculetype") or l.strip().lower() == "[moleculetype]":
            cur = [k]
            blocks.append(cur)
        elif cur is not None:
            cur.append(k)
    fixed, sec_of = [], {}
    for b in blocks:
        sec, refs = None, []
        for k in b:
            code = lines[k].split(";")[0].strip()
            if code.startswith("["):
                sec = code.strip("[] ").lower()
                continue
            w = code.split()
            if not w or w[0].startswith("#") or sec in (None, "moleculetype", "atoms", "meta"):
                continue
            n = NREF.get(sec, None)
            if sec == "virtual_sitesn":
                n = w.index("--") if "--" in w else len([x for x in w if not x.startswith("{")]) - 1
            elif sec == "exclusions":
                n = len([x for x in w if not x.startswith("{")])
            refs.append((k, n))
            sec_of[k] = sec
        if any(w == "0" for k, n in refs for w in lines[k].split(";")[0].split()[:n]):
            fixed.append(lines[b[1]].split()[0] if len(b) > 1 else "?")
            for k, n in refs:
                code, _, comment = lines[k].partition(";")
                w = code.split()
                w[:n] = [str(int(x) + 1) if x.isdigit() else x for x in w[:n]]
                if sec_of[k] == "virtual_sitesn" and "--" not in w:   # "site atoms… function": the function after "--"
                    w = w[:n] + ["--"] + w[n:]
                lines[k] = " ".join(w) + ((" ;" + comment) if comment else "")
    open(dst, "w").write("\n".join(lines) + "\n")
    return fixed


def main():
    cg = {"blocks": {}, "modifications": {}, "links": [], "variables": {}}
    import tempfile
    tmp = os.path.join(tempfile.mkdtemp(), "25-small_molecule_martini3.ff")
    renumbered = zero_based_fixed(os.path.join(SRC, "force_fields", "martini3001", "25-small_molecule_martini3.ff"), tmp)
    vm.read_ff(tmp, cg)
    aa = {"blocks": {}, "modifications": {}, "links": [], "variables": {}}
    vm.read_ff(os.path.join(SRC, "force_fields", "charmm", "10-small_molecule_charmm.ff"), aa)
    out = {"format": "caps-martini3-small-molecules", "version": 1,
           "source": "vermouth-martinize (github.com/marrink-lab/vermouth-martinize, Apache-2.0): force_fields/martini3001/"
                     "25-small_molecule_martini3.ff, mappings/martini3001/*.charmm36.map, force_fields/charmm/10-small_molecule_charmm.ff",
           "references": ["R. Alessandri et al., Adv. Theory Simul. 5, 2100391 (2022)", "P. C. T. Souza et al., Nat. Methods 18, 382 (2021)"],
           "units": "GROMACS: nm, kJ/mol, degrees; parameters as the .ff writes them", "molecules": {}}
    missing = []
    for path in sorted(glob.glob(os.path.join(SRC, "mappings", "martini3001", "*.charmm36.map"))):
        name, beads, atoms = vm.read_map(path)
        if name not in cg["blocks"]:
            continue   # amino acids (the protein model) and others
        if name not in aa["blocks"]:
            missing.append(name)
            continue
        ab = aa["blocks"][name]
        names = [a["atomname"] for a in ab["atoms"]]
        keep = [i for i, a in enumerate(ab["atoms"]) if element(a["atype"])]
        idx = {names[i]: k for k, i in enumerate(keep)}
        bonds = []
        for t in ("bonds", "constraints"):
            for it in ab["interactions"].get(t, []):
                a, b = it["atoms"]
                if a in idx and b in idx:
                    bonds.append(sorted((idx[a], idx[b])))
        mp = {}
        for an, ent in atoms:
            if an not in idx:
                raise SystemExit(f"{path}: atom {an} is not in the CHARMM residue {name}")
            mp[idx[an]] = ent
        blk = cg["blocks"][name]
        bead_names = [a["atomname"] for a in blk["atoms"]]
        for ent in mp.values():
            for b, _ in ent:
                if b not in bead_names:
                    raise SystemExit(f"{path}: bead {b} is not in the Martini 3 block {name}")
        empty = [b for b in bead_names if not any(b == e[0] and e[1] > 0 for ent in mp.values() for e in ent)]
        out["molecules"][name] = {
            "aa": {"atoms": [{"name": names[i], "element": element(ab["atoms"][i]["atype"])} for i in keep], "bonds": sorted(bonds)},
            "map": [[i, mp[i]] for i in sorted(mp)],
            "block": {"atoms": blk["atoms"], "interactions": blk["interactions"], "edges": blk["edges"]},
            "unmapped_beads": empty,   # beads no atom weighs into (virtual sites placed from the others)
        }
    # the blocks against cgmartini's own molecule files (the library's templates), up to the order of the beads: the
    # same bead types and charges, bonds and constraints (atom pair and length), angles, dihedrals and virtual sites
    from itertools import permutations
    tpl = json.load(open(os.path.join(ROOT, "data", "martini", "martini3-molecules.json")))["molecules"]
    same, differ = [], []
    for name, m in out["molecules"].items():
        t = tpl.get(name)
        if not t:
            continue
        order = [a["atomname"] for a in m["block"]["atoms"]]
        A = [(a["atype"], round(float(a.get("charge", 0)), 3)) for a in m["block"]["atoms"]]
        B = [(a["type"], round(float(a["charge"]), 3)) for a in t["atoms"]]
        def ff_terms(p):   # the .ff's terms with bead k renumbered p[k]
            out_ = []
            for t_ in ("bonds", "constraints"):
                for it in m["block"]["interactions"].get(t_, []):
                    if it["meta"].get("ifndef") == "FLEXIBLE":
                        continue
                    out_.append(("b", tuple(sorted(p[order.index(x)] for x in it["atoms"])), round(float(it["params"][1]), 3)))
            for t_ in ("angles", "dihedrals", "impropers"):
                for it in m["block"]["interactions"].get(t_, []):
                    ids = [p[order.index(x)] for x in it["atoms"]]
                    if ids[-1] < ids[0]:
                        ids = ids[::-1]
                    out_.append(("a", tuple(ids), tuple(round(float(x), 2) for x in it["params"][:3])))
            for it in m["block"]["interactions"].get("virtual_sitesn", []):
                ids = [p[order.index(x)] for x in it["atoms"]]
                out_.append(("v", ids[0], tuple(sorted(ids[1:]))))
            return sorted(out_)
        tpl_terms = []
        for x in t["bonds"] + t["constraints"]:
            tpl_terms.append(("b", tuple(sorted(x[:2])), round(float(x[3]), 3)))
        for x in t["angles"]:
            ids = list(x[:3])
            tpl_terms.append(("a", tuple(ids if ids[-1] >= ids[0] else ids[::-1]), tuple(round(float(v), 2) for v in x[3:6])))
        for x in t["dihedrals"]:
            ids = list(x[:4])
            tpl_terms.append(("a", tuple(ids if ids[-1] >= ids[0] else ids[::-1]), tuple(round(float(v), 2) for v in x[4:7])))
        for v in t["vsites"]:
            tpl_terms.append(("v", v["site"], tuple(sorted(v["from"]))))
        tpl_terms.sort()
        found = False
        n = len(A)
        if sorted(A) == sorted(B):
            for p in permutations(range(n)):
                if all(A[k] == B[p[k]] for k in range(n)) and ff_terms(p) == tpl_terms:
                    found = True
                    break
        (same if found else differ).append(name)
    out["renumbered_from_zero"] = renumbered
    out["checked_against_cgmartini"] = {"identical": same, "different": differ,
                                        "note": "molecules in both vermouth's .ff and cgmartini's martini_v3.0.0_small_molecules_v1.itp"}
    dst = os.path.join(ROOT, "data", "martini", "martini3-small-molecules.json")
    json.dump(out, open(dst, "w"), indent=1)
    print(f"{len(out['molecules'])} molecules -> {os.path.relpath(dst, ROOT)} ({os.path.getsize(dst) // 1024} kB); "
          f"no CHARMM residue: {missing}; in cgmartini's file too: identical {same}, different {differ}")


if __name__ == "__main__":
    main()
