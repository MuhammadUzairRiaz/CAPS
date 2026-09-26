#!/usr/bin/env python3
"""MARTINI overlays from the EMC parameter files in the moltemplate distribution (force_fields/martini_original_format):
polymers (PEO, a PEGylated lipid), solvents (alkanes, alcohols, benzene, chloroform ...), surfactants (SDS, DPC, EO5 ...)
and sugars (glucose, sucrose, maltose ...). The library's MARTINI 2.0 file was converted from martini.prm, lipids.prm,
ions.prm and cholesterol.prm only; each other file defines its own numbered bead types (C11, SN01, P41 ...) whose
parameters differ between files, so each becomes its own force field that "extends" MARTINI 2.0 (loaded as one).

Units: nm and kJ/mol in the source. Bonds ½k(r − l0)² → K = k / 2 / 4.184 / 100, r0 = 10 l0; angles (MARTINI's
cosine-harmonic, as the library file) ½k(cos θ − cos θ0)² → K = k / 2 / 4.184; torsions Σ k [1 + cos(nφ − δ)] → fourier
with K = k / 4.184; impropers ½k(ξ − ξ0)² → harmonic K = k / 2 / 4.184. A new type's non-bonded parameters are its
"pair" equivalent's (a MARTINI 2.0 type). Bead templates: the file's TEMPLATES whose beads all have parameters (the
polymer file's PS gives STY / SCY bonds and angles but no masses or pairs: left out). Amino acids: their impropers are read as GROMACS type 2 with the centre second, as every row
names it (the bead bonded to the other three); ILE, LEU, PRO and VAL use AC1 / AC2 beads that no source file gives
masses or pairs for, and are left out.

usage: convert_emc_martini.py [MARTINI_ORIGINAL_FORMAT_DIR]
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FF = os.path.join(ROOT, "data", "forcefields")
TY = os.path.join(ROOT, "data", "typing")
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/moltemplate/moltemplate/force_fields/martini_original_format")
KJ = 4.184
SETS = {"polymers": "PEO and a PEGylated lipid", "solvents": "solvents (alkanes, alcohols, benzene, chloroform ...)",
        "surfactants": "surfactants (SDS, DPC, EO5 ...)", "sugars": "sugars (glucose, sucrose, maltose ...)",
        "aminoacids": "amino acids (single residues)"}


def sections(path):
    out, sec = {}, None
    for l in open(path):
        l = l.rstrip("\n")
        if l.startswith("ITEM\t"):
            sec = l.split("\t")[1]
            if sec == "END":
                sec = None
            else:
                out.setdefault(sec, [])
            continue
        if sec and l.strip() and not l.startswith("#"):
            out[sec].append(l.split())
    return out


def bead_types(smiles):
    return [re.sub(r"[+-]\d*$", "", b) for b in re.findall(r"\[([^\]]+)\]", smiles)]


def short(t):
    m = re.match(r"moltemplate @atom:(\S+)", t.get("description", ""))
    return m.group(1) if m else t["name"]


if __name__ == "__main__":
    base = json.load(open(os.path.join(FF, "martini-moltemplate.json")))
    base_full = {short(t): t["name"] for t in base["atom_types"]}
    base_rules = json.load(open(os.path.join(TY, "martini-moltemplate.typing.json")))["rules"]
    cat_p = os.path.join(FF, "catalogue.json")
    cat = json.load(open(cat_p))
    ids = {e["id"] for e in cat["forcefields"]}
    for key, what in SETS.items():
        sec = sections(os.path.join(SRC, key + ".prm"))
        new = {}
        for w in sec.get("MASS", []):
            if len(w) >= 5:
                new[w[0]] = {"mass": float(w[1]), "charge": float(w[4])}
        # rows: type pair bond angle torsion improper (the header also names an "incr" column the rows do not have:
        # polymers.prm's SN01 → pair SN0, bond SN0 is what lets the PEGylated lipid's SN01-SN02 bond take SN0-SN0)
        equiv = {w[0]: w[1:] for w in sec.get("EQUIVALENCE", []) if len(w) >= 6}

        def pat(t, kind):
            """a name in a bonded row: one of this file's names as is, a MARTINI 2.0 type by its glob"""
            if t in new:
                return t
            if t not in base_full:
                return None
            return {"bond": f"*_b{t}_a*_d*_i*", "angle": f"*_b*_a{t}_d*_i*", "dihedral": f"*_b*_a*_d{t}_i*", "improper": f"*_b*_a*_d*_i{t}"}[kind]

        types = []
        for t, v in sorted(new.items()):
            e = equiv.get(t, [t] * 6)
            pair = e[0]
            if pair not in base_full:
                continue
            eq = {"vdw": base_full[pair]}
            for k, kind in ((1, "bond"), (2, "angle"), (3, "dihedral"), (4, "improper")):
                n = e[k] if k < len(e) else t
                eq[kind] = base_full[n] if n in base_full and n not in new else n
            types.append({"name": t, "element": "", "mass": v["mass"], **({"charge": v["charge"]} if v["charge"] else {}),
                          "description": f"MARTINI {key}.prm bead {t} (non-bonded as {pair})", "equivalence": eq})
        known = set(base_full) | {t["name"] for t in types}
        bonds, angles, dihedrals, impropers, skipped = [], [], [], [], 0
        for w in sec.get("BOND", []):
            m = [pat(x, "bond") for x in w[:2]]
            if None in m or len(w) < 4:
                skipped += 1
                continue
            bonds.append({"name": f"{w[0]}-{w[1]}", "match": m, "style": "harmonic", "params": [round(float(w[2]) / 2 / KJ / 100, 8), round(float(w[3]) * 10, 6)]})
        for w in sec.get("ANGLE", []):
            m = [pat(x, "angle") for x in w[:3]]
            if None in m or len(w) < 5:
                skipped += 1
                continue
            angles.append({"name": "-".join(w[:3]), "match": m, "style": "cosine/squared", "params": [round(float(w[3]) / 2 / KJ, 8), float(w[4])]})
        for w in sec.get("TORSION", []):
            m = [pat(x, "dihedral") for x in w[:4]]
            terms = w[4:]
            if None in m or len(terms) < 3 or len(terms) % 3:
                skipped += 1
                continue
            p = [len(terms) // 3]
            for k in range(0, len(terms), 3):
                p += [round(float(terms[k]) / KJ, 8), int(float(terms[k + 1])), float(terms[k + 2])]
            dihedrals.append({"name": "-".join(w[:4]), "match": m, "style": "fourier", "params": p})
        for w in sec.get("IMPROPER", []):
            m = [pat(x, "improper") for x in w[:4]]
            if None in m or len(w) < 6:
                skipped += 1
                continue
            impropers.append({"name": "-".join(w[:4]), "match": m, "style": "harmonic", "params": [round(float(w[4]) / 2 / KJ, 8), float(w[5])]})
        templates = {}
        left_out = []
        for w in sec.get("TEMPLATES", []):
            if len(w) == 2:
                if set(bead_types(w[1])) <= known:
                    templates[w[0]] = w[1]
                else:
                    left_out.append(w[0])
        fid = f"martini-{key}"
        name = f"MARTINI 2.0 · {key}"
        doc = {"format": "caps-forcefield", "format_version": 1, "name": name, "version": "2.0",
               "source": f"{key}.prm (EMC MARTINI parameters, moltemplate distribution)", "references": [],
               "extends": "martini-moltemplate.json", "typing": f"../typing/{fid}.typing.json",
               "atom_types": types, "pairs": [], "bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": impropers,
               "bead_templates": templates,
               # impropers as GROMACS type 2 (the angle between planes i-j-k and j-k-l): every row names the centre second, its
               # three neighbours around it (HIS: P5 SC4 SP1 SP1, the ring bead bonded to the backbone bead first)
               **({"improper_order": "center2_sorted", "improper_matched_order": True, "improper_max_neighbours": 3} if impropers else {}),
               "notes": [f"converted by CAPS from {key}.prm (bench/ff/convert_emc_martini.py): nm, kJ/mol → Å, kcal/mol; bonds and impropers "
                         "½k → K; angles cosine-harmonic as MARTINI 2.0; new bead types take their pair equivalent's non-bonded terms"]
               + ([f"templates left out (beads without parameters in the file): {', '.join(left_out)}"] if left_out else [])
               + ([f"{skipped} bonded rows name types the file does not define and were skipped"] if skipped else [])}
        json.dump(doc, open(os.path.join(FF, fid + ".json"), "w"), ensure_ascii=False, indent=1)
        shadow = {t["name"] for t in types}   # this file's bead names replace MARTINI 2.0's of the same name
        rules = [r for r in base_rules if r.get("atom_name") not in shadow and r["type"] not in shadow] + [{"type": t["name"], "smarts": "*", "atom_name": t["name"], "priority": 0,
                                     "description": f"a MARTINI {key} bead named {t['name']}"} for t in types]
        tdoc = {"format": "caps-typing", "version": 1, "forcefield": name, "coarse_grained": True, "unknown_types": "untyped",
                "description": f"MARTINI 2.0's bead-name rules with the {key}.prm beads (bench/ff/convert_emc_martini.py)", "rules": rules}
        json.dump(tdoc, open(os.path.join(TY, fid + ".typing.json"), "w"), ensure_ascii=False, indent=1)
        entry = {"id": fid, "name": name, "version": "2.0", "references": ["S. J. Marrink et al., J. Phys. Chem. B 111, 7812 (2007)"],
                 "origin": "moltemplate", "source_file": f"martini_original_format/{key}.prm", "status": "converted",
                 "notes": f"MARTINI 2.0 with the {what} of {key}.prm ({len(templates)} bead templates); extends MARTINI 2.0",
                 "file": fid + ".json",
                 "counts": {"atom_types": len(types), "pairs": 0, "bonds": len(bonds), "angles": len(angles), "dihedrals": len(dihedrals), "impropers": len(impropers)},
                 "typing": {"rules": f"typing/{fid}.typing.json", "evidence": "MARTINI 2.0 bead-name rules plus this file's beads"}}
        if fid in ids:
            cat["forcefields"] = [entry if e["id"] == fid else e for e in cat["forcefields"]]
        else:
            k = next(i for i, e in enumerate(cat["forcefields"]) if e["id"] == "martini-moltemplate")
            cat["forcefields"].insert(k + 1, entry)
            ids.add(fid)
        print(f"{fid}: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} torsions, {len(impropers)} impropers; "
              f"{len(templates)} templates (left out: {', '.join(left_out) or 'none'}); {skipped} rows skipped")
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
