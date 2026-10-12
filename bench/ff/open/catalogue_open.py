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
    "clandp": ("CL&P (ionic liquids)", "il.ff of 2026/04/16 (github.com/paduagroup/clandp)",
               ["J. N. Canongia Lopes, A. A. H. Padua, Theor. Chem. Acc. 131, 1129 (2012)"],
               {"source": "il.ff, github.com/paduagroup/clandp", "license": "MIT (Copyright (c) 2019 Agilio Padua)"},
               "CAPS rules for CL&P's types by ion family (bench/typing/make_clandp_rules.py): CL&P's own ion files (.zmat, .xyz) "
               "typed as they type them, 1165 / 1165 atoms in 55 molecules (P111i4's two swapped labels listed; CT2 and CS, the same "
               "parameters, counted as one)", ("CL&P ionic liquids", "opls_cl_p", 12), ["opls-clp"]),
    "dreiding1990": ("DREIDING", "Mayo, Olafson, Goddard 1990 (Table I and the paper's rules)",
                     ["S. L. Mayo, B. D. Olafson, W. A. Goddard III, J. Phys. Chem. 94, 8897 (1990)"],
                     {"source": "the paper: J. Phys. Chem. 94, 8897 (1990), Table I (radii, angles, van der Waals) and its rules for "
                                "bonds, angles, torsions, inversions and hydrogen bonds (bench/ff/make_dreiding.py)",
                      "license": "published parameters of the paper"},
                     "CAPS rules for the paper's labels (element and hybridisation: _3, _R, _2, _1; H___A on N or O); LAMMPS "
                     "reproduces CAPS's energies and forces (bond constants by bond order, cosine-harmonic angles, the torsion rules, "
                     "planar inversions, hydrogen bonds)", ("DREIDING", "dreiding", 14), ["dreiding"]),
    "trappe-ua-foyer": ("TraPPE-UA (alkanes, alcohols)", "foyer trappe-ua.xml 0.0.2 (mosdef-hub/foyer, 2021-08-11)",
                        ["M. G. Martin, J. I. Siepmann, J. Phys. Chem. B 102, 2569 (1998)", "M. G. Martin, J. I. Siepmann, J. Phys. Chem. B 103, 4508 (1999)",
                         "B. Chen, J. J. Potoff, J. I. Siepmann, J. Phys. Chem. B 105, 3093 (2001)"],
                        {"source": "foyer/forcefields/xml/trappe-ua.xml, github.com/mosdef-hub/foyer (bench/ff/convert_foyer.py)",
                         "license": "MIT (Copyright (c) 2015 Vanderbilt University)"},
                        "CAPS rules for foyer's eleven united-atom types (alkane and alcohol CHn, O, H); its Ryckaert-Bellemans torsions "
                        "written exactly as Fourier terms; LAMMPS and GROMACS reproduce CAPS's energies. Torsions foyer's file lacks "
                        "(CHn-CH2-CH-OH and similar, C-C between branched centres) are reported missing",
                        ("TraPPE united-atom", "trappe_ua", 16), ["trappe-ua"]),
    "gromos54a7-gromacs": ("GROMOS 54A7", "54A7 as GROMACS distributes it (gromos54a7.ff)",
                           ["N. Schmid et al., Eur. Biophys. J. 40, 843 (2011)", "C. Oostenbrink et al., J. Comput. Chem. 25, 1656 (2004)"],
                           {"source": "GROMACS share/top/gromos54a7.ff (ffnonbonded.itp, ffbonded.itp, atomtypes.atp, aminoacids.rtp, .n.tdb, .c.tdb)",
                            "license": "GNU LGPL 2.1 (part of GROMACS)"},
                           "a building-block force field: CAPS recognises the residues of aminoacids.rtp (with the termini of the .tdb "
                           "files) in a structure and takes their types, charges and coded terms as pdb2gmx does; peptides of all 20 "
                           "amino acids, a proline N-terminus and a glycine zwitterion give pdb2gmx's topology term for term and its "
                           "energies (bench/ff/compare_pdb2gmx.py)", ("GROMOS 54A7", "gromos_54a7", 9), ["gromos-54a7"]),
    "trappe-ua-umn": ("TraPPE-UA", "the TraPPE database (trappe.oit.umn.edu), retrieved October 2026",
                      ["M. G. Martin, J. I. Siepmann, J. Phys. Chem. B 102, 2569 (1998)",
                       "B. L. Eggimann et al., J. Chem. Eng. Data 65, 1330 (2020)"],
                      {"source": "the TraPPE group's own parameter tables, model by model (trappe.oit.umn.edu), each term with its paper; "
                                 "transcribed by bench/ff/make_trappe.py",
                       "license": "published parameters of the TraPPE papers, cited per term (no file of the site is distributed)"},
                      "CAPS rules for 92 TraPPE-UA site types: the database's own models typed as it types them (617 / 617 sites); "
                      "every model's bonds, bends and torsions as it lists them (87 / 87 models; nine models whose own lists are "
                      "inconsistent named apart); 1-4 Coulomb x 0.5 for the glycol ethers, diols, disulfides and acrylates as their "
                      "models count it; LAMMPS and GROMACS reproduce CAPS's energies", ("TraPPE united-atom", "trappe_ua", 16), ["trappe-ua"]),
    "bush1994-oxides": ("Oxides (Bush et al. 1994)", "Bush, Gale, Catlow & Battle 1994 shell-model set",
                        ["T. S. Bush, J. D. Gale, C. R. A. Catlow, P. D. Battle, J. Mater. Chem. 4, 831 (1994)"],
                        {"source": "the paper's set as the UCL Database of Published Interatomic Potential Parameters gives it (bush.lib, "
                                   "mcc.chem.ucl.ac.uk/Potentials); transcribed by bench/ff/make_inorganic_gulp.py",
                         "license": "published parameters of the paper (no file of the database is distributed)"},
                        "one shell-model set for 20 cations with oxygen, typed by element (a shell on every core); every Buckingham "
                        "term equals an independent transcription of the same paper; LAMMPS reproduces CAPS's energies on MgO, "
                        "SrTiO3 and corundum (CORESHELL styles)", ("Oxides (Bush 1994)", "bush1994", 30), []),
    "clayff2004": ("ClayFF", "Cygan, Liang & Kalinichev 2004 (Tables 1 and 2)",
                   ["R. T. Cygan, J.-J. Liang, A. G. Kalinichev, J. Phys. Chem. B 108, 1255 (2004)"],
                   {"source": "the paper's Tables 1 (nonbond) and 2 (bonds, angles); bench/ff/make_clayff.py",
                    "license": "published parameters of the paper"},
                   "CAPS rules by coordination and neighbours type LAMMPS's ClayFF pyrophyllite as its file does (1280 / 1280 atoms); "
                   "LAMMPS reproduces CAPS's energies on it (metal-O-H bends through contacts)", ("ClayFF", "clayff", 26), ["inorganic-clay"]),
    "gan2020": ("GaN (Wonglakhon & Zahn 2020)", "2020 (Table 1)", ["T. Wonglakhon, D. Zahn, J. Phys.: Condens. Matter 32, 205401 (2020)"],
                {"source": "the paper's Table 1; bench/ff/make_gan.py", "license": "published parameters of the paper (the article is CC BY)"},
                "Ga and N ions with formal charges, Lennard-Jones Ga-Ga and N-N, Born-Mayer Ga-N; LAMMPS reproduces CAPS's energies on "
                "wurtzite", ("GaN", "gan", 31), ["inorganic-binary-misc"]),
    "catlow1977-halides": ("Alkali halides (Catlow, Diller & Norgett 1977)", "1977, Table 2 (set 1, rigid ions)",
                           ["C. R. A. Catlow, K. M. Diller, M. J. Norgett, J. Phys. C: Solid State Phys. 10, 1395 (1977)"],
                           {"source": "the paper's Table 2; bench/ff/make_catlow1977.py", "license": "published parameters of the paper"},
                           "16 halides of Li, Na, K, Rb with F, Cl, Br, I, each typed by the structure's composition; every pair term "
                           "equals an independent transcription of the same table; LAMMPS reproduces CAPS's energies on NaCl and LiF",
                           ("Alkali halides", "alkali_halides", 28), ["inorganic-binary-halides"]),
}
RTF = "the topology's own residues typed as it types them (bench/ff/validate_rtf_types.py --cap-links): "
for fid, name, version, refs, evidence, listing, replaces in (
        ("charmm36-prot-toppar", "CHARMM36 proteins", "toppar_c36_jul17 (par_all36_prot.prm, top_all36_prot.rtf, its streams)",
         ["R. B. Best et al., J. Chem. Theory Comput. 8, 3257 (2012)", "A. D. MacKerell Jr. et al., J. Phys. Chem. B 102, 3586 (1998)"],
         RTF + "top_all36_prot.rtf 391 / 391 atoms, 24 / 24 residues", ("CHARMM36", "charmm36", 7), ["charmm36-prot"]),
        ("charmm22-prot-toppar", "CHARMM22 proteins", "toppar_c36_jul17 (par_all22_prot.prm, top_all22_prot.rtf)",
         ["A. D. MacKerell Jr. et al., J. Phys. Chem. B 102, 3586 (1998)"],
         RTF + "top_all22_prot.rtf 381 / 381 atoms, 23 / 23 residues", ("CHARMM22 protein", "charmm22_prot", 21), ["charmm22-prot"]),
        ("charmm36-carb-toppar", "CHARMM36 carbohydrates", "toppar_c36_jul17 (par_all36_carb.prm, top_all36_carb.rtf, its streams)",
         ["O. Guvench et al., J. Chem. Theory Comput. 7, 3162 (2011)"],
         RTF + "top_all36_carb.rtf 1786 / 1786 atoms, 75 / 75 residues", ("CHARMM36 carbohydrate", "charmm36_carb", 22), ["charmm36-carb"]),
        ("cgenff-toppar", "CHARMM36 CGenFF", "toppar_c36_jul17 (par_all36_cgenff.prm, top_all36_cgenff.rtf: CGenFF 3.0.1 / 4.0)",
         ["K. Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010)"],
         RTF + "top_all36_cgenff.rtf 9872 / 9872 atoms, 519 / 519 residues (residues with lone pairs left out)",
         ("CHARMM36 CGenFF", "charmm36_cgenff", 23), ["cgenff"]),
        ("charmm36-lipid-toppar", "CHARMM36 lipids", "toppar_c36_jul17 (par_all36_lipid.prm, top_all36_lipid.rtf, its streams)",
         ["J. B. Klauda et al., J. Phys. Chem. B 114, 7830 (2010)"],
         RTF + "top_all36_lipid.rtf 4224 / 4224 atoms, 35 / 35 residues", ("CHARMM36 lipid", "charmm36_lipid", 24), ["charmm36-lipid"]),
        ("charmm36-nucl-toppar", "CHARMM36 nucleic acids", "toppar_c36_jul17 (par_all36_na.prm, top_all36_na.rtf, its streams)",
         ["K. Hart et al., J. Chem. Theory Comput. 8, 348 (2012)"],
         RTF + "top_all36_na.rtf 161 / 161 atoms, 5 / 5 residues", ("CHARMM36 nucleic acid", "charmm36_nucl", 25), ["charmm36-nucl"]),
        ("charmm19-toppar", "CHARMM19 (united atom)", "param19 / toph19 (toppar_c36_jul17)",
         ["E. Neria, S. Fischer, M. Karplus, J. Chem. Phys. 105, 1902 (1996)"],
         RTF + "toph19.inp 319 / 319 atoms, 36 / 36 residues (--united); one torsion per bond as toph19 lists them, 103 / 104 central "
         "bonds (bench/ff/validate_rtf_torsions.py; TRP's cross-ring terms left out)", ("CHARMM19", "charmm19", 20), ["charmm19"])):
    ENTRIES[fid] = (name, version, refs, CHARMM, evidence, listing, replaces)
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
