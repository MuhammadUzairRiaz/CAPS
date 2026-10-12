#!/usr/bin/env python3
"""Catalogue entries (data/forcefields/catalogue.json) for the libraries CAPS builds from the force fields' own published
files (bench/ff/open/build_*.sh): source, version, citation and licence for each, the Studio listing moved to it from the
entry it replaces, and the replaced entry removed.

  catalogue_open.py [FF_DIR]
"""
import json, os, sys

FF = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "..", "..", "data", "forcefields")
GAFF = "J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, D. A. Case, J. Comput. Chem. 25, 1157 (2004)"
AMBERTOOLS = {"source": "AmberTools 23.6, dat/leap/parm (ambermd.org)",
              "license": "public domain (AmberTools LICENSE: the parameter files in dat/leap/parm have been placed in the public domain by their authors)"}
CHARMM = {"source": "CHARMM force-field files toppar_c36_jul17.tgz (mackerell.umaryland.edu)",
          "license": "no licence stated with the files: built on the user's machine by bench/ff/open/build_charmm.sh, never distributed"}

# id: (name, version, references, source description, typing evidence, list (label, key, order) or None, replaces)
ENTRIES = {
    "gaff-1.81-ambertools": ("GAFF 1.81", "1.81 (gaff.dat, May 2017)", [GAFF], AMBERTOOLS,
                             "CAPS rules in the order of antechamber's ATOMTYPE_GFF.DEF (AmberTools 23.6), with its ring classes and "
                             "conjugated pairs (bench/typing/make_gaff_rules.py)", ("GAFF", "gaff", 8), ["gaff-amber16"]),
    "gaff2-2.11-ambertools": ("GAFF2 2.11", "2.11 (gaff211.dat, May 2016)", [GAFF], AMBERTOOLS,
                              "CAPS rules in the order of antechamber's ATOMTYPE_GFF2.DEF, without c5 / c6 (added in GAFF2 2.2)", None, []),
    "gaff2-2.2.20-ambertools": ("GAFF2", "2.2.20 (gaff2.dat, March 2021)", [GAFF], AMBERTOOLS,
                                "CAPS rules in the order of antechamber's ATOMTYPE_GFF2.DEF (AmberTools 23.6), with its ring classes and "
                                "conjugated pairs (bench/typing/make_gaff_rules.py)", ("GAFF2", "gaff2", 18), ["gaff-amber25"]),
    "amber-ff03-ambertools": ("AMBER ff03", "parm99.dat + frcmod.ff03",
                              ["Y. Duan et al., J. Comput. Chem. 24, 1999 (2003)", "J. Wang, P. Cieplak, P. A. Kollman, J. Comput. Chem. 21, 1049 (2000)"],
                              AMBERTOOLS, "CAPS rules in AMBER's own type names (bench/typing/make_amber_rules.py); AmberTools' residue "
                              "libraries all_amino03.lib, all_aminont94.lib, all_aminoct94.lib and all_nucleic94.lib typed as the "
                              "libraries type them: 2270 / 2270 atoms, 111 / 111 residues (bench/ff/validate_off_types.py)",
                              ("AMBER ff03", "amber", 17), ["amber"]),
}
cat_path = os.path.join(FF, "catalogue.json")
cat = json.load(open(cat_path))
fl = cat["forcefields"]
for fid, (name, version, refs, lic, evidence, listing, replaces) in ENTRIES.items():
    path = os.path.join(FF, fid + ".json")
    e = {"id": fid, "name": name, "version": version, "references": refs, "origin": "built from source", **lic,
         "status": "validated", "file": fid + ".json"}
    if os.path.exists(path):
        d = json.load(open(path))
        e["counts"] = {"atom_types": len(d["atom_types"]), **{k: len(d.get(k, [])) for k in ("pairs", "bonds", "angles", "dihedrals", "impropers")}}
    if evidence:
        e["typing"] = {"rules": f"typing/{fid}.typing.json", "evidence": evidence}
    if listing:
        e["list"] = {"label": listing[0], "key": listing[1], "order": listing[2]}
    fl[:] = [x for x in fl if x["id"] not in replaces and x["id"] != fid]
    fl.append(e)
json.dump(cat, open(cat_path, "w"), indent=1, ensure_ascii=False)
print(f"{cat_path}: {len(fl)} entries")
