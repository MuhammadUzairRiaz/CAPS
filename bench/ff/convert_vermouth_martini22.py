#!/usr/bin/env python3
"""Martini 2.2 proteins for CAPS, from vermouth-martinize (github.com/marrink-lab/vermouth-martinize, Apache-2.0):
vermouth/data/force_fields/martini22/05-aminoacids.ff (residue beads, types, charges, side-chain bonds, constraints,
angles and dihedrals; the backbone "links") and vermouth/data/mappings/*.charmm36.map (which atoms make each bead).

Written to data/martini/martini22-protein.json:
  residues   per residue: beads (name, type, charge), side-chain terms in GROMACS units (nm, kJ/mol, degrees), the atoms
             of each bead (heavy atoms by name; "h" counts the hydrogens the map lists on each heavy atom: those join the
             bead, the others are not part of its centre, as vermouth averages the listed atoms by mass)
  links      the backbone rules of the .ff as tables CAPS applies in order (later rules override earlier ones, as
             vermouth applies its links): bead types by secondary structure, backbone bonds / constraints, BBB angles,
             BBS and first SBB angles, helix dihedrals, the elastic bonds of extended regions, disulfide constraints,
             termini. Every number is checked against the line of the .ff it comes from.
  ss         vermouth's conversion of DSSP letters to Martini's (dssp.py: SS_CG and the helix patterns)

usage: convert_vermouth_martini22.py [VERMOUTH_REF_DIR]   (default ~/vermouth-ref, the downloaded files)
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/vermouth-ref")
FF = os.path.join(SRC, "force_fields", "martini22", "05-aminoacids.ff")
MODS = os.path.join(SRC, "force_fields", "martini22", "10-modifications.ff")
MAPS = os.path.join(SRC, "mappings")
ff_text = open(FF).read()
mods_text = open(MODS).read()


def source(line, text=ff_text):
    """the .ff must hold this line (whitespace-insensitive): every number below comes from one"""
    norm = lambda t: re.sub(r"\s+", " ", t.strip())
    if not any(norm(l).startswith(norm(line)) for l in text.splitlines()):
        raise SystemExit(f"not in the source: {line!r}")
    return True


def moleculetypes():
    out, cur, sec = {}, None, None
    for raw in ff_text.splitlines():
        l = raw.split(";")[0].strip()
        if not l or l.startswith("#"):
            continue
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1).lower()
            if sec == "link":
                cur = None
            continue
        w = l.split()
        if sec == "moleculetype":
            cur = {"beads": [], "bonds": [], "constraints": [], "angles": [], "dihedrals": []}
            out[w[0]] = cur
            sec = None
        elif cur is None:
            continue
        elif sec == "atoms":
            cur["beads"].append({"name": w[4], "type": w[1], "charge": float(w[6])})
        elif sec == "bonds":
            cur["bonds"].append([int(w[0]) - 1, int(w[1]) - 1, float(w[3]), float(w[4])])
        elif sec == "constraints":
            cur["constraints"].append([int(w[0]) - 1, int(w[1]) - 1, float(w[3])])
        elif sec == "angles":
            cur["angles"].append([int(w[0]) - 1, int(w[1]) - 1, int(w[2]) - 1, int(w[3]), float(w[4]), float(w[5])])
        elif sec == "dihedrals":
            cur["dihedrals"].append([int(w[0]) - 1, int(w[1]) - 1, int(w[2]) - 1, int(w[3]) - 1, int(w[4]), float(w[5]), float(w[6])] + ([int(w[7])] if len(w) > 7 else []))
    return out


def mapping(res):
    """{bead: [[heavy atom, listed hydrogens on it], ...]} from the .charmm36.map"""
    f = os.path.join(MAPS, res.lower() + ".charmm36.map")
    if not os.path.exists(f):
        return None
    beads, sec, last = {}, None, None
    for raw in open(f):
        l = raw.split(";")[0].strip()
        m = re.match(r"\[\s*(\S+)\s*\]", l)
        if m:
            sec = m.group(1).lower()
            continue
        w = l.split()
        if sec != "atoms" or len(w) < 3 or not w[0].isdigit():
            continue
        name, bead = w[1], w[2]
        if name.startswith("H"):
            if last is not None:
                last[1] += 1
            continue
        last = [name, 0]
        beads.setdefault(bead, []).append(last)
    return beads


SS_CG = {"1": "H", "2": "H", "3": "H", "H": "H", "G": "H", "I": "H", "B": "E", "E": "E", "T": "T", "S": "S", "C": "C", " ": "C", "-": "C"}
HELIX_PATTERNS = [[".H.", ".3."], [".HH.", ".33."], [".HHH.", ".333."], [".HHHH.", ".3333."], [".HHHHH.", ".13332."],
                  [".HHHHHH.", ".113322."], [".HHHHHHH.", ".1113222."], [".HHHH", ".1111"], ["HHHH.", "2222."]]

if __name__ == "__main__":
    mt = moleculetypes()
    residues = {}
    for res, t in mt.items():
        mp = mapping({"HSD": "hsd", "HSE": "hse", "HSP": "hsp"}.get(res, res))
        # neutral forms map as the charged residue (the same atoms; the protonation decides the name)
        if mp is None:
            mp = mapping({"ASP0": "asp", "GLU0": "glu", "LSN": "lys", "ARG0": "arg"}.get(res, res))
        if mp is None:
            print(f"{res}: no mapping in the source (left out)")
            continue
        t["atoms"] = mp
        residues[res] = t
    # residue modifications that map more atoms into a bead (mappings/modifications.mapping): LYS-HZ3 (a charged lysine)
    # and LYS-LSN (a neutral one) put CE's hydrogens into SC2 as well
    mtxt = open(os.path.join(MAPS, "modifications.mapping")).read()
    for block in re.split(r"\[\s*modification\s*\]", mtxt):
        m = re.search(r"\[\s*from blocks\s*\]\s*\n\s*(\S+)", block)
        mm = re.search(r"\[\s*mapping\s*\]\s*\n(.*?)(?:\n\s*\[|\Z)", block, re.S)
        if not m or not mm or not m.group(1).startswith("LYS-"):
            continue
        target = {"LYS-HZ3": "LYS", "LYS-LSN": "LSN"}.get(m.group(1))
        if target not in residues:
            continue
        pairs = [l.split() for l in mm.group(1).splitlines() if len(l.split()) == 2 and not l.startswith(";")]
        for atom, bead in pairs:
            if atom.startswith("H"):
                continue
            hs = sum(1 for a, b in pairs if a.startswith("H") and b == bead and a[1:2] == atom[1:2])   # HE1 HE2 on CE
            for e in residues[target]["atoms"].get(bead, []):
                if e[0] == atom and hs > e[1]:
                    e[1] = hs
                    print(f"{target} {bead}: {atom} with {hs} hydrogens ({m.group(1)} modification mapping)")
    # the backbone links (05-aminoacids.ff, ";;; Links"), each checked against its line
    source("BB +BB 1 0.350 1250"); source("-BB  BB  SC1  2 100 25"); source("SC1 BB +BB 2 100 25")
    source('BB {"cgsecstruct": "T|3|E", "replace": {"atype": "Nda"}, "modifications": null}')
    source('BB {"cgsecstruct": "2", "replace": {"atype": "Na"}, "modifications": null}')
    source('BB {"cgsecstruct": "1", "replace": {"atype": "Nd"}, "modifications": null}')
    source('BB {"cgsecstruct": "H|F", "replace": {"atype": "N0"}, "modifications": null}')
    source('BB {"cgsecstruct": "S", "replace": {"atype": "P4"}, "modifications": null}')
    source('BB {"cgsecstruct": "T|3|2|1|E", "replace": {"atype": "N0"}, "modifications": null}')
    source('BB {"cgsecstruct": "H|F", "replace": {"atype": "C5"}, "modifications": null}')
    source("BB +BB 1 0.310"); source("BB +BB 1 0.33"); source("-BB BB +BB 2 96 700"); source('-BB BB {"cgsecstruct": "E"} +BB 2 134 25')
    source('-BB BB {"cgsecstruct": "S"} +BB 2 130 20'); source("-BB BB +BB 2 98 100"); source("-BB BB +BB 2 134 25")
    source("-BB BB +BB 2 130 25"); source("-BB BB +BB 2 127 25"); source("-BB BB +BB 2 100 25"); source("-BB BB +BB 2 130 20")
    source("-BB BB +BB 2 127 20"); source("-BB BB +BB 2 100 20"); source('-BB BB {"cgsecstruct": "H|1|2|3", "resname": "PRO|HYP"} +BB 2 98 100')
    source("-BB BB +BB ++BB 1 -120 400 1"); source('BB ++BB 1 0.640 2500'); source('BB +++BB 1 0.970 2500'); source("SC1 >SC1 1 0.24")
    source('BB {"replace": {"atype": "Qa", "charge": -1}}', mods_text); source('BB {"replace": {"atype": "Qd", "charge": 1}}', mods_text)
    links = {
        "bb_type": [   # in order; a residue with a terminus modification keeps the terminus type ("modifications": null)
            {"res": "*", "ss": "T3E", "type": "Nda"}, {"res": "*", "ss": "2", "type": "Na"}, {"res": "*", "ss": "1", "type": "Nd"},
            {"res": "*", "ss": "HF", "type": "N0"},
            {"res": "ALA PRO HYP", "ss": "S", "type": "P4"}, {"res": "ALA PRO HYP", "ss": "T321E", "type": "N0"},
            {"res": "ALA PRO HYP", "ss": "HF", "type": "C5"}, {"res": "PRO", "ss": "2", "type": "Na"}],
        "bb_bond": {"coil": [0.350, 1250], "helix_helix_constraint": 0.310, "helix_other_constraint": 0.33, "helix": "H123",
                    "other": "SCTE"},
        # BBB angles in the .ff's order: [pattern, degrees, kJ/mol]; a rule applies when any of the three residues has
        # one of the codes (and, for "pro" / "nonpro", is or is not a proline); "C" also when no code is set; "mid_pro_helix":
        # the middle residue is a helical proline. The last rule that applies wins.
        "bbb_angle": [
            {"ss": "H123", "res": "*", "theta": 96, "k": 700}, {"ss": "E", "res": "*", "theta": 134, "k": 25},
            {"ss": "S", "res": "*", "theta": 130, "k": 20},
            {"ss": "H123", "res": "PRO HYP", "theta": 98, "k": 100}, {"ss": "E", "res": "PRO HYP", "theta": 134, "k": 25},
            {"ss": "S", "res": "PRO HYP", "theta": 130, "k": 25}, {"ss": "C", "res": "PRO HYP", "theta": 127, "k": 25},
            {"ss": "T", "res": "PRO HYP", "theta": 100, "k": 25},
            {"ss": "S", "res": "!PRO HYP", "theta": 130, "k": 20}, {"ss": "C", "res": "*", "theta": 127, "k": 20, "or_none": True},
            {"ss": "T", "res": "!PRO HYP", "theta": 100, "k": 20},
            {"ss": "H123", "res": "PRO HYP", "theta": 98, "k": 100, "middle_only": True}],
        "bbs_angle": [100, 25], "first_sbb_angle": [100, 25],
        "helix_dihedral": {"ss": "H123", "phi": -120, "k": 400, "n": 1},
        "extended_elastic": {"ss": "E", "short": [0.640, 2500], "long": [0.970, 2500]},
        "disulfide_constraint": 0.24,
        "termini": {"N": {"type": "Qd", "charge": 1.0}, "C": {"type": "Qa", "charge": -1.0}},
        # constraints become stiff bonds (CAPS has no constraints): 20000 kJ/mol/nm2, the value the MARTINI parameter files
        # of the moltemplate distribution give the same distances (aminoacids.prm: AC1-P51 20000 0.31, SC4-SC4 20000 0.27)
        "constraint_k": 20000,
    }
    doc = {"format": "caps-martini-protein", "version": 1, "model": "Martini 2.2 proteins",
           "source": "vermouth-martinize (github.com/marrink-lab/vermouth-martinize, Apache-2.0): force_fields/martini22/05-aminoacids.ff, "
                     "10-modifications.ff, mappings/*.charmm36.map, dssp/dssp.py",
           "references": ["S. J. Marrink et al., J. Phys. Chem. B 111, 7812 (2007)", "L. Monticelli et al., J. Chem. Theory Comput. 4, 819 (2008)",
                          "D. H. de Jong et al., J. Chem. Theory Comput. 9, 687 (2013)",
                          "P. C. Kroon et al., Martinize2 and Vermouth, eLife 12, RP90627 (2023)"],
           "units": "GROMACS: nm, kJ/mol, degrees; bonds and impropers 1/2 k, angles 1/2 k (cos - cos0)^2",
           "ss": {"dssp_to_cg": SS_CG, "helix_patterns": HELIX_PATTERNS},
           "residues": residues, "links": links}
    os.makedirs(os.path.join(ROOT, "data", "martini"), exist_ok=True)
    out = os.path.join(ROOT, "data", "martini", "martini22-protein.json")
    json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
    print(f"{out}: {len(residues)} residues ({', '.join(residues)})")
    # the force field: MARTINI 2.0's beads and pairs, typing that maps a protein and names the beads by type
    FFD = os.path.join(ROOT, "data", "forcefields")
    fid = "martini22-proteins"
    ff = {"format": "caps-forcefield", "format_version": 1, "name": "Martini 2.2 proteins", "version": "2.2",
          "source": "vermouth-martinize martini22 (05-aminoacids.ff, mappings) on the library's MARTINI 2.0", "references": doc["references"],
          "extends": "martini-moltemplate.json", "typing": f"../typing/{fid}.typing.json", "atom_types": [], "pairs": [], "bonds": [],
          "angles": [], "dihedrals": [], "impropers": [],
          "notes": ["an all-atom protein is mapped onto Martini 2.2 beads with martinize's explicit topology (data/martini/martini22-protein.json; "
                    "DSSP for the secondary structure); non-bonded terms are MARTINI 2.0's (the library's martini.lt)",
                    "constraints are stiff bonds (20000 kJ/mol/nm2): CAPS has no constraints",
                    "AC1 / AC2 (VAL, LEU, ILE side chains) have no non-bonded parameters in the library (martini_v2.2.itp gives them)",
                    "ring (S) beads weigh 54 in the library's MARTINI file (EMC); Martini 2.2's itp gives 45: masses change dynamics, not energies"]}
    json.dump(ff, open(os.path.join(FFD, fid + ".json"), "w"), ensure_ascii=False, indent=1)
    base_rules = json.load(open(os.path.join(ROOT, "data", "typing", "martini-moltemplate.typing.json")))["rules"]
    tdoc = {"format": "caps-typing", "version": 1, "forcefield": ff["name"], "coarse_grained": True, "unknown_types": "untyped",
            "martini_protein": "../martini/martini22-protein.json",
            "description": "Martini 2.2 proteins (bench/ff/convert_vermouth_martini22.py): an all-atom protein becomes beads named by their "
                           "Martini type, with martinize's topology; beads typed by name with MARTINI 2.0's rules",
            "rules": base_rules}
    json.dump(tdoc, open(os.path.join(ROOT, "data", "typing", fid + ".typing.json"), "w"), ensure_ascii=False, indent=1)
    cat_p = os.path.join(FFD, "catalogue.json")
    cat = json.load(open(cat_p))
    entry = {"id": fid, "name": ff["name"], "version": "2.2", "references": doc["references"], "origin": "vermouth-martinize",
             "source_file": "force_fields/martini22/05-aminoacids.ff", "status": "converted",
             "notes": "proteins mapped as martinize2 does (term by term identical on vermouth's martini22 test), DSSP included; VAL / LEU / ILE lack AC1 / AC2 pairs",
             "file": fid + ".json", "counts": {"atom_types": 0, "pairs": 0, "bonds": 0, "angles": 0, "dihedrals": 0, "impropers": 0},
             "typing": {"rules": f"typing/{fid}.typing.json", "evidence": "bench/ff/check_martini_protein.py: identical to martinize2 (types, charges, every term, positions)"}}
    ids = [e["id"] for e in cat["forcefields"]]
    if fid in ids:
        cat["forcefields"][ids.index(fid)] = entry
    else:
        cat["forcefields"].insert(ids.index("martini-aminoacids") + 1, entry)
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    print(f"{fid}: force field, typing and catalogue entry written")
