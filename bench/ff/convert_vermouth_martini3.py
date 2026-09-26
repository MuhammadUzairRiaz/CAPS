#!/usr/bin/env python3
"""Martini 3 proteins for CAPS, from vermouth-martinize (github.com/marrink-lab/vermouth-martinize, Apache-2.0):
force_fields/martini3001/{05-aminoacids,10-modifications,15-links_M3,20-links_M3_IDR}.ff and
mappings/martini3001/*.charmm36.map + modifications.charmm36.mapping.

The .ff files are read as vermouth reads them (vermouth/ffinput.py): macros substituted in the text; blocks with their
atoms, interactions (dihedrals of function 2 become impropers) and the edges their bonds, constraints, angles and proper
dihedrals make (unless "edge": false); modifications with the attributes they replace; links with their nodes (the link's
own attributes applied to every node, "+" / "-" / ">" prefixes as the order), edges, non-edges, patterns, molecule meta,
interactions to add and interactions to remove, parameters measured from the structure (dihphase). CAPS applies the links
in file order with vermouth's matching rules (core/src/martini3_protein.cpp).

Attribute values: strings, numbers, booleans, null; a "A|B" choice becomes {"$choice": ["A", "B"]}.
A measured parameter becomes {"$effector": "dihphase", "keys": [...], "format": ".01f"}.

The mappings become, per residue, the heavy atoms of each bead with their weight (0 for "!" entries) and the weight of the
hydrogens on each heavy atom (the map's hydrogens are attached by CHARMM's naming: HB1 / HB2 on CB, HD1 on CD1 or ND1 ...);
the modification mappings add the terminal atoms (OXT, the extra hydrogens on N) to the backbone bead.

usage: convert_vermouth_martini3.py [VERMOUTH_REF_DIR]   (default ~/vermouth-ref, the downloaded files)
"""
import glob, json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/vermouth-ref")
FFDIR = os.path.join(SRC, "force_fields", "martini3001")
MAPDIR = os.path.join(SRC, "mappings", "martini3001")
FILES = ["05-aminoacids.ff", "10-modifications.ff", "15-links_M3.ff", "20-links_M3_IDR.ff"]
NATOMS = {"bonds": 2, "constraints": 2, "pairs": 2, "angles": 3, "dihedrals": 4, "impropers": 4, "cmap": 5,
          "virtual_sites2": 3, "virtual_sites3": 4, "virtual_sites4": 5}
EDGE_TYPES = ("bonds", "angles", "dihedrals", "cmap", "constraints")


def tokenize(line):
    """whitespace-separated tokens; a {...} JSON object is one token"""
    out, i = [], 0
    while i < len(line):
        if line[i].isspace():
            i += 1
            continue
        if line[i] == "{":
            depth, j, instr = 0, i, False
            while j < len(line):
                c = line[j]
                if c == '"' and line[j - 1] != "\\":
                    instr = not instr
                elif not instr and c == "{":
                    depth += 1
                elif not instr and c == "}":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            out.append(line[i:j + 1])
            i = j + 1
        else:
            j = i
            while j < len(line) and not line[j].isspace():
                j += 1
            out.append(line[i:j])
            i = j
    return out


def choice(v):
    return {"$choice": v.split("|")} if isinstance(v, str) and "|" in v else v


def attrs_of(token):
    a = json.loads(token)
    return {k: choice(v) for k, v in a.items()}


def split_key(key):
    n = 0
    while n < len(key) and key[n] in "+-><*":
        n += 1
    prefix, base = key[:n], key[n:]
    if not base or len(set(prefix)) > 1:
        raise SystemExit(f"bad atom key {key!r}")
    return prefix, base


def treat_prefix(ref, attrs):
    """vermouth's _treat_atom_prefix: the key's prefix is the order; the base is the atom name unless given"""
    prefix, base = split_key(ref)
    a = dict(attrs)
    if "order" not in a:
        a["order"] = (len(prefix) if prefix[:1] == "+" else -len(prefix) if prefix[:1] == "-" else prefix) if prefix else 0
        key = ref
    else:
        o = a["order"]
        key = ref if prefix else ((("+" if o > 0 else "-") * abs(o)) if isinstance(o, int) else o) + base
    a.setdefault("atomname", base)
    return key, a


def get_atoms(tokens, natoms):
    atoms = []
    while tokens:
        if tokens[0] == "--":
            tokens.pop(0)
            break
        if natoms is not None and len(atoms) >= natoms:
            break
        t = tokens.pop(0)
        if tokens and tokens[0].startswith("{"):
            atoms.append([t, attrs_of(tokens.pop(0))])
        else:
            atoms.append([t, {}])
    return atoms


def parse_params(tokens):
    out = []
    for t in tokens:
        if "(" in t and not t.startswith("(") and t.endswith(")"):
            name, arg = t.split("(", 1)
            arg = arg[:-1]
            fmt = None
            if "|" in arg:
                arg, fmt = arg.split("|")
            out.append({"$effector": name, "keys": [x.strip() for x in arg.split(",")], "format": fmt})
        else:
            out.append(t)
    return out


class Context:
    def __init__(self, kind, name=None):
        self.kind, self.name = kind, name
        self.nodes = {}            # key -> attrs (ordered)
        self.edges = []
        self.non_edges = []
        self.patterns = []
        self.interactions = {}
        self.removed = {}
        self.all_nodes = {}        # the link's own attributes (header lines)
        self.all_inter = {}        # section -> #meta attributes
        self.molmeta = {}
        self.features = []
        self.nrexcl = None

    def add_edge(self, a, b):
        for k in (a, b):
            self.nodes.setdefault(k, {})
        if [a, b] not in self.edges and [b, a] not in self.edges:
            self.edges.append([a, b])

    def link_atoms(self, atoms, section):
        refs = []
        for ref, at in atoms:
            a = dict(self.all_nodes)
            a.update(at)
            key, a = treat_prefix(ref, a)
            refs.append(key)
            if key in self.nodes:
                for k, v in a.items():
                    if k in self.nodes[key] and self.nodes[key][k] != v:
                        raise SystemExit(f"conflicting attribute {k} for {ref} in a {self.kind}")
                self.nodes[key].update(a)
            else:
                self.nodes[key] = a
        return refs

    def interaction(self, section, tokens, delete):
        natoms = NATOMS.get(section)
        atoms = get_atoms(tokens, natoms)
        if self.kind == "block":
            names = list(self.nodes)
            refs = []
            for ref, _ in atoms:
                if ref.isdigit():
                    ref = names[int(ref) - 1]
                if ref not in self.nodes:
                    raise SystemExit(f"no atom {ref} in block {self.name}")
                refs.append(ref)
        else:
            refs = self.link_atoms(atoms, section)
        meta = {}
        if tokens and tokens[-1].startswith("{"):
            meta = json.loads(tokens.pop())
        params = parse_params(tokens)
        full = dict(self.all_inter.get(section, {}))
        full.update(meta)
        entry = {"atoms": refs, "params": params, "meta": full}
        if delete:
            entry["atom_attrs"] = [a for _, a in atoms]
            self.removed.setdefault(section, []).append(entry)
        else:
            # proper and improper dihedrals share [ dihedrals ]: function 2 is an improper
            if section == "dihedrals" and params and params[0] == "2":
                section = "impropers"
            self.interactions.setdefault(section, []).append(entry)

    def finish(self):
        for t in EDGE_TYPES:
            for it in self.interactions.get(t, []):
                if it["meta"].get("edge", True):
                    for a, b in zip(it["atoms"][:-1], it["atoms"][1:]):
                        self.add_edge(a, b)


def read_ff(path, out):
    text = open(path).read()
    lines = [l.split(";")[0].rstrip() for l in text.splitlines()]
    macros = {}
    section = []
    cur = None

    def close():
        nonlocal cur
        if cur is None:
            return
        if cur.kind == "block":
            cur.finish()
            out["blocks"][cur.name] = {"atoms": list(cur.nodes.values()), "interactions": cur.interactions, "edges": cur.edges,
                                       "nrexcl": cur.nrexcl}
        elif cur.kind == "link":
            cur.finish()
            out["links"].append({"nodes": cur.nodes, "edges": cur.edges, "non_edges": cur.non_edges, "patterns": cur.patterns,
                                 "interactions": cur.interactions, "removed": cur.removed, "molmeta": cur.molmeta,
                                 "features": cur.features, "source": os.path.basename(path)})
        elif cur.kind == "modification":
            out["modifications"][cur.name] = {"atoms": cur.nodes, "edges": cur.edges}
        cur = None

    for raw in lines:
        l = raw.strip()
        if not l:
            continue
        m = re.match(r"^\[\s*(\S+)\s*\]$", l)
        if m:
            name = m.group(1).lower()
            if name in ("moleculetype", "link", "modification", "macros", "variables", "citations"):
                close()
                section = [name]
                if name == "link":
                    cur = Context("link")
                elif name == "modification":
                    cur = Context("modification")
                elif name == "moleculetype":
                    cur = Context("block")
            else:
                section = section[:1] + [name]
            continue
        # macros: $name replaced in the text (after the [ macros ] section itself)
        if section == ["macros"]:
            k, v = l.split(None, 1)
            macros[k] = v.strip()
            continue
        for k in sorted(macros, key=len, reverse=True):
            l = l.replace("$" + k, macros[k])
        top, sub = section[0], (section[1] if len(section) > 1 else None)
        tokens = tokenize(l)
        if top == "variables":
            k, v = tokens
            try:
                out["variables"][k] = json.loads(v)
            except ValueError:
                out["variables"][k] = v
            continue
        if top == "citations" or sub in ("citation", "citations", "meta", "info"):
            continue
        if top == "moleculetype" and sub is None:
            cur.name, cur.nrexcl = tokens[0], int(tokens[1])
            continue
        if top == "link" and sub in (None, "molmeta"):
            k, v = tokens
            val = json.loads(v)
            (cur.molmeta if sub == "molmeta" else cur.all_nodes)[k] = choice(val)
            continue
        if top == "modification" and sub is None:
            cur.name = l.strip()
            continue
        if sub == "atoms":
            if cur.kind == "block":
                attrs = attrs_of(tokens.pop()) if tokens[-1].startswith("{") else {}
                _, atype, resid, resname, name, cg = tokens[:6]
                a = {"atomname": name, "atype": atype, "resname": resname, "resid": int(resid), "charge_group": int(cg)}
                if len(tokens) > 6:
                    a["charge"] = float(tokens[6])
                if len(tokens) > 7:
                    a["mass"] = float(tokens[7])
                a.update(attrs)
                cur.nodes[name] = a
            else:
                ref, at = tokens[0], attrs_of(tokens[1])
                key, at = treat_prefix(ref, at)
                full = dict(cur.all_nodes)
                full.update(at)
                if cur.kind == "modification":
                    full.setdefault("PTM_atom", False)
                cur.nodes[key] = {**cur.nodes.get(key, {}), **full}
            continue
        if sub == "edges":
            atoms = get_atoms(tokens, 2)
            keys = []
            for ref, at in atoms:
                key, _ = treat_prefix(ref, at)
                keys.append(key)
            if cur.kind == "block":
                cur.add_edge(*keys)
            else:
                cur.add_edge(*keys)
            continue
        if sub == "non-edges":
            atoms = get_atoms(tokens, 2)
            k0, _ = treat_prefix(*atoms[0])
            _, a1 = treat_prefix(*atoms[1])
            full = dict(cur.all_nodes)
            full.update(a1)
            cur.non_edges.append([k0, full])
            continue
        if sub == "patterns":
            cur.patterns.append(get_atoms(tokens, None))
            continue
        if sub == "features":
            cur.features.extend(tokens)
            continue
        delete = sub.startswith("!")
        isec = sub[1:] if delete else sub
        if tokens[0] == "#meta":
            cur.all_inter.setdefault(isec, {}).update(json.loads(tokens[1]))
            continue
        cur.interaction(isec, tokens, delete)
    close()


# ------------------------------------------------------------------------------------------------ mappings
def heavy_of(h, heavies):
    """the heavy atom a CHARMM hydrogen name belongs to: HN on N, HA* on CA, else 'H' + the heavy atom's name after its
    element letter (HB1 on CB, HD21 on ND2, HH11 on NH1, HG1 on OG / OG1 / SG ...), the longest match"""
    if h in ("HN", "H", "HT1", "HT2", "HT3", "HN1", "HN2", "HN3"):
        return "N"
    best = None
    for x in heavies:
        suf = x[1:]
        if not suf and x == "N":
            continue
        if h[1:].startswith(suf) and (best is None or len(suf) > len(best[1:])):
            best = x
    if best is None and h.startswith("HA"):
        best = "CA"
    return best


def read_map(path):
    beads, atoms, name, sec = [], [], None, None
    for raw in open(path):
        l = raw.split(";")[0].strip()
        m = re.match(r"^\[\s*(\S+)\s*\]$", l)
        if m:
            sec = m.group(1).lower()
            continue
        if not l:
            continue
        w = l.split()
        if sec == "molecule":
            name = w[0]
        elif sec == "martini":
            beads = w
        elif sec == "atoms":
            # vermouth's weights (map_input._compute_weights): an atom listed n times for a bead among m non-null
            # entries weighs n / m in it; "!" entries weigh 0
            nonnull = [b for b in w[2:] if not b.startswith("!")]
            ent = [[b, nonnull.count(b) / len(nonnull)] for b in dict.fromkeys(nonnull)]
            ent += [[b[1:], 0.0] for b in dict.fromkeys(w[2:]) if b.startswith("!")]
            atoms.append([w[1], ent])
    return name, beads, atoms


def residue_mapping(path, block):
    name, _, atoms = read_map(path)
    heavy = [a for a, _ in atoms if not a.startswith("H")]
    # heavy atoms: [name, [[bead, weight], ...]]; h: the hydrogens on each heavy atom, [[bead, weight], ...] (every
    # hydrogen the map lists on that atom weighs the same; hydrogens the map weighs 0 contribute nothing)
    out = {"beads": [a["atomname"] for a in block["atoms"]], "heavy": [], "h": {}}
    hw = {}
    for a, ent in atoms:
        if a.startswith("H"):
            x = heavy_of(a, heavy)
            if x is None:
                raise SystemExit(f"{path}: no heavy atom for {a}")
            live = tuple(sorted((b, round(w, 12)) for b, w in ent if w > 0))
            hw.setdefault(x, set()).add(live)
        else:
            out["heavy"].append([a, ent])
    for x, s in hw.items():
        if len(s) > 1:
            raise SystemExit(f"{path}: the hydrogens on {x} weigh differently: {s}")
        out["h"][x] = [list(e) for e in next(iter(s))]
    return name, out


def mod_mappings(path):
    """terminal modifications: the atoms they add to beads (beyond the residue's own)"""
    out, cur, sec = {}, None, None
    for raw in open(path):
        l = raw.split(";")[0].strip()
        m = re.match(r"^\[\s*(.+?)\s*\]$", l)
        if m:
            sec = m.group(1).lower()
            continue
        if not l:
            continue
        w = l.split()
        if sec == "to blocks":
            cur = out.setdefault(w[0], [])
        elif sec == "atoms" or sec == "mapping":
            if len(w) >= 2 and not w[0].startswith("charmm") and cur is not None:
                cur.append([w[-2], w[-1]])
    return out


def main():
    out = {"format": "caps-vermouth-model", "version": 1, "model": "Martini 3 proteins (martini3001)",
           "source": "vermouth-martinize (github.com/marrink-lab/vermouth-martinize, Apache-2.0): force_fields/martini3001/"
                     "05-aminoacids.ff, 10-modifications.ff, 15-links_M3.ff, 20-links_M3_IDR.ff; mappings/martini3001/*.charmm36.map, "
                     "modifications.charmm36.mapping; dssp/dssp.py",
           "references": ["P. C. T. Souza et al., Nat. Methods 18, 382 (2021)",
                          "P. C. T. Souza et al., Martini 3 protein models (scfix, IDP: Nat. Commun. 2024)",
                          "P. C. Kroon et al., Martinize2 and Vermouth, eLife 12, RP90627 (2023)"],
           "units": "GROMACS: nm, kJ/mol, degrees; parameters as the .ff writes them",
           "variables": {"center_weight": "mass", "bb_atomname": "BB"}, "blocks": {}, "modifications": {}, "links": []}
    for f in FILES:
        read_ff(os.path.join(FFDIR, f), out)
    # the secondary-structure conversion is vermouth's, the same for every Martini version (dssp/dssp.py)
    m22 = json.load(open(os.path.join(ROOT, "data", "martini", "martini22-protein.json")))
    out["ss"] = m22["ss"]
    out["mapping"] = {}
    for path in sorted(glob.glob(os.path.join(MAPDIR, "*.charmm36.map"))):
        name, _, _ = read_map(path)
        if name not in out["blocks"] or out["blocks"][name]["atoms"][0].get("atomname") != "BB":
            continue   # small molecules (their own work) and residues without a Martini 3 block
        rname, mp = residue_mapping(path, out["blocks"][name])
        out["mapping"][rname] = mp
    mm = mod_mappings(os.path.join(MAPDIR, "modifications.charmm36.mapping"))
    out["mod_mapping"] = {k: v for k, v in mm.items() if k in out["modifications"]}
    # the hydrogens a modification mapping puts in a bead with weight 1, by the heavy atom they sit on (LYS-HZ3 and
    # LYS-LSN put CE's hydrogens in SC2, which the residue's own map weighs 0)
    out["mod_h"] = {}
    for k, v in out["mod_mapping"].items():
        # hydrogens attach to the residue's own heavy atoms (LYS-HZ3: LYS's), the termini's to the backbone's
        resn = k.split("-")[0]
        own = [a for a, _ in out["mapping"][resn]["heavy"]] if resn in out["mapping"] else ["N", "CA", "C", "O"]
        heavy = list(dict.fromkeys([a for a, _ in v if not a.startswith("H")] + own))
        hs = {}
        for a, b in v:
            if a.startswith("H"):
                x = heavy_of(a, heavy)
                if x is None:
                    raise SystemExit(f"modification {k}: no heavy atom for {a}")
                hs[x] = b
        out["mod_h"][k] = hs
    dst = os.path.join(ROOT, "data", "martini", "martini3-protein.json")
    json.dump(out, open(dst, "w"), indent=1)
    nlinks = len(out["links"])
    print(f"{len(out['blocks'])} blocks, {len(out['modifications'])} modifications, {nlinks} links, "
          f"{len(out['mapping'])} residue mappings -> {os.path.relpath(dst, ROOT)} ({os.path.getsize(dst) // 1024} kB)")


if __name__ == "__main__":
    main()
