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
    "amber-ff14sb-ambertools": ("AMBER ff14SB / OL15 / OL3", "ff14SB proteins, OL15 DNA, OL3 RNA (parm10.dat + frcmod.ff14SB + frcmod.DNA.OL15)",
                                ["J. A. Maier et al., J. Chem. Theory Comput. 11, 3696 (2015)", "M. Zgarbova et al., J. Chem. Theory Comput. 11, 5723 (2015)",
                                 "M. Zgarbova et al., J. Chem. Theory Comput. 7, 2886 (2011)"], AMBERTOOLS,
                                "CAPS rules in AMBER's own type names (bench/typing/make_amber_rules.py ff14sb); AmberTools' residue libraries "
                                "amino12, aminont12, aminoct12, DNA.OL15 and RNA typed as the libraries type them: 2288 / 2288 atoms, 111 / 111 "
                                "residues (bench/ff/validate_off_types.py)", ("AMBER ff14SB", "amber14sb", 16), []),
    "oplsaa-gromacs": ("OPLS-AA (GROMACS oplsaa.ff)", "OPLS-AA/L as GROMACS distributes it (oplsaa.ff)",
                       ["W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)",
                        "G. A. Kaminski, R. A. Friesner, J. Tirado-Rives, W. L. Jorgensen, J. Phys. Chem. B 105, 6474 (2001)"],
                       {"source": "GROMACS share/top/oplsaa.ff (ffnonbonded.itp, ffbonded.itp)", "license": "GNU LGPL 2.1 (part of GROMACS)"},
                       "CAPS's OPLS-AA rules in opls_NNN names plus GROMACS's protein types (bench/typing/make_opls_gromacs_rules.py): "
                       "aminoacids.rtp 504 / 504 atoms, 46 / 46 residues (PGLU's entry listed as a deviation), methanol.itp exact; "
                       "every parameter CAPS assigns equals GROMACS's own lookup (bench/ff/check_source_gromacs.py)", ("OPLS-AA/L (GROMACS)", "oplsaal", 6), []),
    "oplsaa2024": ("OPLS-AA 2024", "Jorgensen's BOSS parameter files of October 2023 (J. Phys. Chem. B 128, 250 (2024), SI)",
                   ["W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)",
                    "W. L. Jorgensen et al., J. Phys. Chem. B 128, 250 (2024)"],
                   {"source": "supporting information jp3c06602_si_002.txt (.par) and jp3c06602_si_003.txt (.sb), pubs.acs.org/doi/suppl/10.1021/acs.jpcb.3c06602",
                    "license": "published parameter tables of the paper's supporting information (moltemplate, MIT, redistributes them as public domain)"},
                   "CAPS's OPLS-AA rules by type number (bench/typing/make_oplsaa_rules.py) in opls_NNN names; on the polymer library "
                   "every term equals moltemplate's independent conversion of the same files except CT-CT-OH-HO, where that conversion "
                   "applies trifluoroethanol's torsion to every alcohol (bench/ff/convert_boss_opls.py applies the general line)",
                   ("OPLS-AA 2024", "oplsaa2024", 4), ["opls2020"]),
    "oplsaam-jorgensen": ("OPLS-AA/M", "OPLS-AA/M of June 2015 (the Jorgensen lab's CHARMM-format files)",
                          ["M. J. Robertson, J. Tirado-Rives, W. L. Jorgensen, J. Chem. Theory Comput. 11, 3499 (2015)",
                           "W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)"],
                          {"source": "top_opls_aam.inp and par_opls_aam.inp, jorgensenresearch.com/opls-aa-m",
                           "license": "CC BY 4.0 (the Jorgensen lab: free for anyone to use, including commercially, with attribution)"},
                          "CAPS's OPLS-AA rules in OPLS-AA/M's names plus the protein and terminal types (bench/typing/make_oplsaam_rules.py): "
                          "the file's own topology typed as it types it, 508 / 508 atoms, its 25 residues and 11 terminal patches "
                          "(bench/ff/validate_rtf_types.py --patches)", ("OPLS-AA/M", "opls_aam", 11), ["opls-aam"]),
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
    if listing:   # one entry per listing key: an older entry with this key keeps its file but leaves the list
        for x in fl:
            if x.get("list", {}).get("key") == listing[1]:
                x.pop("list")
    fl.append(e)
json.dump(cat, open(cat_path, "w"), indent=1, ensure_ascii=False)
print(f"{cat_path}: {len(fl)} entries")
