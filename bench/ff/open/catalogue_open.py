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
                    "license": "published parameter tables of the paper's supporting information"},
                   "CAPS's OPLS-AA rules by type number (bench/typing/make_oplsaa_rules.py) in opls_NNN names; where the .par gives a "
                   "torsion again for particular molecules, the first, general line is applied (bench/ff/convert_boss_opls.py)",
                   ("OPLS-AA 2024", "oplsaa2024", 4), ["opls2020"]),
    "clandp": ("CL&P (ionic liquids)", "il.ff of 2026/04/16 (github.com/paduagroup/clandp)",
               ["J. N. Canongia Lopes, A. A. H. Padua, Theor. Chem. Acc. 131, 1129 (2012)"],
               {"source": "il.ff, github.com/paduagroup/clandp", "license": "MIT (Copyright (c) 2019 Agilio Padua)"},
               "CAPS rules for CL&P's types by ion family (bench/typing/make_clandp_rules.py): CL&P's own ion files (.zmat, .xyz) "
               "typed as they type them, 1165 / 1165 atoms in 55 molecules (P111i4's two swapped labels listed; CT2 and CS, the same "
               "parameters, counted as one)", ("CL&P ionic liquids", "opls_cl_p", 12), ["opls-clp"]),
    "clandp-des": ("CL&P DES (deep-eutectic solvents)", "des.ff of 2025/12/19 (github.com/paduagroup/clandp)",
                   ["J. N. Canongia Lopes, A. A. H. Padua, Theor. Chem. Acc. 131, 1129 (2012)",
                    "W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)"],
                   {"source": "des.ff, github.com/paduagroup/clandp", "license": "MIT (Copyright (c) 2019 Agilio Padua)"},
                   "CAPS rules for des.ff's types by DES component (bench/typing/make_clandp_des_rules.py): CL&P's own molecule "
                   "files whose types des.ff defines (cholinium.xyz, the ammonium and phosphonium ions, Cl, Br, acetate) typed as "
                   "they type them, 362 / 362 atoms in 14 molecules, each at its integer charge; 29 DES components (cholinium, "
                   "ammonium, phosphonium, betaine, halides, carboxylates, urea, alcohols, diols, glycerol, carboxylic acids, "
                   "levulinic acid, ketones) typed whole and neutral or at their ion charge (bench/ff/check_clandp_des_components.py); "
                   "those des.ff's charges cannot make neutral (methanol, malonic and oxalic acid) left untyped",
                   ("CL&P deep-eutectic solvents", "opls_des", 13), ["opls-des"]),
    "oplsdes-doherty2018": ("OPLS-DES (Doherty & Acevedo 2018)", "2018 (SI Tables S1-S6)",
                            ["B. Doherty, O. Acevedo, J. Phys. Chem. B 122, 9982 (2018)",
                             "S. V. Sambasivarao, O. Acevedo, J. Chem. Theory Comput. 5, 1038 (2009)",
                             "W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)"],
                            {"source": "the paper and its Supporting Information (Tables S1-S6, Schemes S1-S2); "
                                       "bench/ff/make_opls_des_doherty2018.py",
                             "license": "published parameters of the paper"},
                            "CAPS rules for the paper's nine molecules as the SI's Schemes S1 and S2 type them, each rule only inside "
                            "its whole molecule (bench/typing/make_oplsdes_doherty2018_rules.py): choline, chloride, ethylene glycol, "
                            "glycerol, levulinic, malonic and oxalic acid, phenol, urea, 102 / 102 atoms by parameters (100 by name: "
                            "urea's cis H, HC, has HT's terms and is typed HT); choline +0.8, chloride -0.8, the HBDs neutral, no "
                            "missing term (bench/ff/check_oplsdes_doherty2018_components.py); LAMMPS and GROMACS reproduce CAPS's "
                            "energies on reline, glyceline and the other HBDs", ("OPLS-DES (Doherty & Acevedo 2018)", "oplsdes_doherty2018", 13), []),
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
    "trappe-ua-umn": ("TraPPE-UA", "the TraPPE database (trappe.oit.umn.edu), retrieved October 2026",
                      ['M. G. Martin, J. I. Siepmann, J. Phys. Chem. B 102, 2569 (1998)', 'M. G. Martin, J. I. Siepmann, J. Phys. Chem. B 103, 4508 (1999)', 'C. D. Wick, M. G. Martin, J. I. Siepmann, J. Phys. Chem. B 104, 8008 (2000)', 'B. Chen, J. J. Potoff, J. I. Siepmann, J. Phys. Chem. B 105, 3093 (2001)', 'J. M. Stubbs, J. J. Potoff, J. I. Siepmann, J. Phys. Chem. B 108, 17596 (2004)', 'S. J. Keasler, S. M. Charan, C. D. Wick, I. G. Economou, J. I. Siepmann, J. Phys. Chem. B 116, 11234 (2012)', 'B. L. Eggimann et al., J. Chem. Eng. Data 65, 1330 (2020)'],
                      {"source": "the TraPPE group's own parameter tables, model by model (trappe.oit.umn.edu), each term with its paper, and the transferable hydrocarbon bonds, bends and torsions of the TraPPE-UA papers' tables; transcribed by bench/ff/make_trappe.py",
                       "license": 'published parameters of the TraPPE papers, cited per term (no file of the site is distributed)'},
                      "CAPS rules for 106 TraPPE-UA site types: the database's own models typed as it types them (634 / 634 sites); every model's bonds, bends and torsions as it lists them (92 / 92 models; ten models whose own lists are inconsistent or incomplete named apart); the transferable terms of the hydrocarbon (1998, 1999, 2000), alcohol (2001), ether and glycol (2004) and cycloalkane-ring (2012) papers (109 bonds, 845 bends, 2776 torsions) equal every database model's term for the same types; 30 of the 114 library polymers parameterised whole as DP-4 chains (bench/ff/trappe_polymer_coverage.py); 1-4 Coulomb x 0.5 for the glycol ethers, diols, disulfides and acrylates as their models count it; LAMMPS and GROMACS reproduce CAPS's energies", ("TraPPE united-atom", "trappe_ua", 16), ["trappe-ua"]),
    "bush1994-oxides": ("Oxides (Bush et al. 1994)", "Bush, Gale, Catlow & Battle 1994 shell-model set",
                        ["T. S. Bush, J. D. Gale, C. R. A. Catlow, P. D. Battle, J. Mater. Chem. 4, 831 (1994)"],
                        {"source": "the paper's set as the UCL Database of Published Interatomic Potential Parameters gives it (bush.lib, "
                                   "mcc.chem.ucl.ac.uk/Potentials); transcribed by bench/ff/make_inorganic_gulp.py",
                         "license": "published parameters of the paper (no file of the database is distributed)"},
                        "one shell-model set for 20 cations with oxygen, typed by element (a shell on every core); every Buckingham "
                        "term equals an independent transcription of the same paper; LAMMPS reproduces CAPS's energies on MgO, "
                        "SrTiO3 and corundum (CORESHELL styles)", ("Oxides (Bush 1994)", "bush1994", 30), []),
    "clayff-cygan2004": ("ClayFF", "Cygan, Liang & Kalinichev 2004 (Tables 1 and 2, p 1259)",
                         ["R. T. Cygan, J.-J. Liang, A. G. Kalinichev, J. Phys. Chem. B 108, 1255 (2004)"],
                         {"source": "the paper's Table 1 (charges, Do, Ro) and Table 2 (O-H bonds, H-O-H and metal-O-H bends), "
                                    "p 1259, read from its page images; bench/ff/make_clayff_cygan2004.py",
                          "license": "published parameters of the paper"},
                         "CAPS rules by element, coordination and neighbours type the INTERFACE force field's clay models "
                         "(pyrophyllite, kaolinite, muscovite, five montmorillonites; 10652 / 10652 atoms) as an independent "
                         "reference typing of ClayFF does; a muscovite with Mg for every third octahedral Al 1241 / 1260, the other "
                         "19 the paper's obss (an O by a tetrahedral and an octahedral substitution, p 1259) where that reference has "
                         "obos; LAMMPS and GROMACS reproduce CAPS's energies (metal-O-H bends through contacts)", ("ClayFF", "clayff", 26), ["clayff2004", "inorganic-clay"]),
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
    "lewis-catlow1985-oxides": ("Oxides (Lewis & Catlow 1985)", "1985, Table 3 (shell model, method (c)(i)); O-O from Catlow 1977 (C 20.37, his potential 2's)",
                                ["G. V. Lewis, C. R. A. Catlow, J. Phys. C: Solid State Phys. 18, 1149 (1985)", "C. R. A. Catlow, Proc. R. Soc. Lond. A 353, 533 (1977)"],
                                {"source": "the papers' tables; bench/ff/make_lewis_catlow1985.py", "license": "published parameters of the papers"},
                                "20 oxides, each a shell-model set typed by composition (tetrahedral sites and second oxidation states as "
                                "overlays); every cation-O term equals an independent transcription (27 / 27); LAMMPS reproduces CAPS's "
                                "energies on MgO and CaO", ("Oxides (Lewis & Catlow 1985)", "lewis_catlow1985", 31), []),
    "catlow1977-uo2": ("UO2 (Catlow 1977)", "1977, Table 3 (potentials 1 and 2)", ["C. R. A. Catlow, Proc. R. Soc. Lond. A 353, 533 (1977)"],
                       {"source": "the paper's Table 3; bench/ff/make_catlow1977_uo2.py", "license": "published parameters of the paper"},
                       "shell model for U and O; potential 1 by default, potential 2 as an overlay; LAMMPS reproduces CAPS's energies "
                       "on fluorite UO2 (CORESHELL styles)", ("UO2 (Catlow 1977)", "catlow1977_uo2", 31), []),
    "matsui1991-tio2": ("TiO2 (Matsui & Akaogi 1991)", "1991, Table 1", ["M. Matsui, M. Akaogi, Mol. Simul. 6, 239 (1991)"],
                        {"source": "the paper's Table 1; bench/ff/make_matsui1991.py", "license": "published parameters of the paper"},
                        "rigid Ti and O ions with partial charges; Ti-Ti, Ti-O and O-O Buckingham terms from the per-atom parameters "
                        "(Ti-O 16957.7 eV, 0.194 A, 12.59 eV A^6, as the widely used transcription); LAMMPS reproduces CAPS's "
                        "energies on rutile", ("TiO2 (Matsui & Akaogi)", "matsui1991", 31), []),
    "sanders1984-sio2": ("SiO2 (Sanders, Leslie & Catlow 1984)", "1984, Table 1",
                         ["M. J. Sanders, M. Leslie, C. R. A. Catlow, J. Chem. Soc., Chem. Commun. 1271 (1984)"],
                         {"source": "the paper's Table 1; bench/ff/make_sanders1984.py", "license": "published parameters of the paper"},
                         "rigid Si4+, shell-model O, Buckingham Si-O and O-O, the O-Si-O bend on the O shells within 1.9 A of each Si "
                         "as the paper writes it (K_B (theta - theta0)^2; GULP-style files apply 1/2 K_B with the same number); LAMMPS "
                         "reproduces CAPS's energies on alpha-quartz", ("SiO2 (Sanders 1984)", "sanders1984", 31), []),
    "catlow1976-mgo": ("MgO (Catlow, Faux & Norgett 1976)", "1976, Table 1 (models 1 and 2)",
                       ["C. R. A. Catlow, I. D. Faux, M. J. Norgett, J. Phys. C: Solid State Phys. 9, 419 (1976)"],
                       {"source": "the paper's Table 1(b), 1(c); bench/ff/make_catlow1976_mgo.py", "license": "published parameters of the paper"},
                       "rigid Mg2+ and a shell-model O2-; model 1 by default, model 2 as an overlay; the perfect lattice at r0 2.106 A "
                       "gives -40.82 eV per MgO (the paper: 40.85); LAMMPS reproduces CAPS's energies on rock-salt MgO with either model",
                       ("MgO (Catlow 1976)", "catlow1976_mgo", 31), []),
    "catlow1985-rutile": ("Rutile (Catlow, Freeman & Royle 1985)", "1985, Table IVa",
                          ["C. R. A. Catlow, C. M. Freeman, R. L. Royle, Physica B+C 131, 1 (1985)"],
                          {"source": "the paper's Table IVa; bench/ff/make_catlow1985_rutile.py", "license": "published parameters of the paper"},
                          "Ti and O shells with harmonic and quartic core-shell springs (the quartic constant read as K w^4/24); the "
                          "paper's alpha-quartz set is sanders1984-sio2; LAMMPS reproduces CAPS's energies on rutile",
                          ("Rutile (Catlow 1985)", "catlow1985_rutile", 31), []),
    "khan1998-zro2": ("ZrO2 (Khan, Islam & Bates 1998)", "1998, Table 1", ["M. S. Khan, M. S. Islam, D. R. Bates, J. Mater. Chem. 8, 2299 (1998)"],
                      {"source": "the paper's Table 1; bench/ff/make_khan1998_zro2.py", "license": "published parameters of the paper"},
                      "Zr and O shells, formal charges, 10.254 A cut-off; cubic ZrO2 at a = 5.07 A gives -109.76 eV per ZrO2 (the paper's "
                      "lattice energy); LAMMPS reproduces CAPS's energies", ("ZrO2 (Khan 1998)", "khan1998_zro2", 31), []),
    "sayle2000-alkaline-earth-oxides": ("MgO, CaO, SrO, BaO rigid ions (Sayle & Watson 2000)", "2000, Table 1",
                                        ["D. C. Sayle, G. W. Watson, Phys. Chem. Chem. Phys. 2, 5491 (2000)"],
                                        {"source": "the paper's Table 1; bench/ff/make_sayle2000.py", "license": "published parameters of the paper"},
                                        "formal-charge rigid ions, one O for the four oxides and their interfaces, 10 A cut-off; LAMMPS "
                                        "reproduces CAPS's energies on MgO and BaO", ("Alkaline-earth oxides (Sayle 2000)", "sayle2000", 31), []),
    "sepliarsky2005-perovskites": ("Ferroelectric perovskites (Sepliarsky et al. 2005)", "2005, Table 1 (K(Nb,Ta)O3, (Ba,Sr)TiO3)",
                                   ["M. Sepliarsky, A. Asthagiri, S. R. Phillpot, M. G. Stachiotti, R. L. Migoni, Curr. Opin. Solid State "
                                    "Mater. Sci. 9, 107 (2005)"],
                                   {"source": "the paper's Table 1; bench/ff/make_sepliarsky2005.py", "license": "published parameters of the paper"},
                                   "every ion a core and a shell with anharmonic springs c2 w^2/2 + c4 w^4/24, Buckingham pairs, each set "
                                   "typed by composition; PbTiO3 left out (Table 1's Pb-Ti row makes it unbound); LAMMPS reproduces CAPS's "
                                   "energies on KNbO3, KTaO3, BaTiO3 and SrTiO3", ("Perovskites (Sepliarsky 2005)", "sepliarsky2005", 30), []),
    "bandura2003-tio2-water": ("TiO2-H2O (Bandura & Kubicki 2003)", "2003, Tables 2, 3, 5, 6, 9, 10",
                               ["A. V. Bandura, J. D. Kubicki, J. Phys. Chem. B 107, 11072 (2003)"],
                               {"source": "the paper's tables; bench/ff/make_bandura2003.py", "license": "published parameters of the paper"},
                               "Matsui & Akaogi's rigid ions with the paper's surface charges and pairs, sites typed by neighbours (TiV, "
                               "TiVI, TiVIoh, Ob, Os, Obh, Ooh), Morse O-H, Ti-O-H bends by distance, flexible SPC/E water; LAMMPS "
                               "reproduces CAPS's energies on a hydroxylated rutile (110) slab with water", ("TiO2-water (Bandura 2003)", "bandura2003", 31), []),
    "teter-glass": ("Oxide glasses, Teter partial charges", "Cormack 2002, Du 2009, Okhotnikov 2013 (each term with its table)",
                    ["A. N. Cormack, J. Du, T. R. Zeitler, Phys. Chem. Chem. Phys. 4, 3193 (2002)", "J. Du, J. Am. Ceram. Soc. 92, 87 (2009)",
                     "K. Okhotnikov, B. Stevensson, M. Eden, Phys. Chem. Chem. Phys. 15, 15041 (2013)",
                     "A. Tilocca, N. H. de Leeuw, A. N. Cormack, Phys. Rev. B 73, 104209 (2006)"],
                    {"source": "the papers' tables; bench/ff/make_glass.py", "license": "published parameters of the papers"},
                    "rigid ions with partial charges (O -1.2, Si 2.4, Al/Y/La/Lu/Sc 1.8, Na 0.6), Buckingham terms with O; Okhotnikov's "
                    "own Y-O and Tilocca's Teter set (Table III, with Ca) as overlays; LAMMPS reproduces CAPS's energies on quartz, "
                    "corundum and Y2O3, each overlay included", ("Glass (Teter partial charges)", "teter_glass", 30), ["inorganic-glass"]),
    "tilocca2006-sm": ("Silicate glasses, shell model (Tilocca 2006)", "2006, Table I (SM1; SM0 an overlay)",
                       ["A. Tilocca, N. H. de Leeuw, A. N. Cormack, Phys. Rev. B 73, 104209 (2006)"],
                       {"source": "the paper's Table I; bench/ff/make_glass.py", "license": "published parameters of the paper"},
                       "formal-charge Si, Na, Ca and shell-model O; SM1's screened O-Si-O bend (CAPS and DL_POLY; LAMMPS and GROMACS have "
                       "no such style), its energy on quartz equal to an independent evaluation of the paper's equation; SM0 (raw "
                       "harmonic bend) as an overlay, LAMMPS reproducing CAPS's energies on quartz", ("Glass (Tilocca shell model)", "tilocca2006_sm", 30), []),
    "vogel2003-lipo3": ("LiPO3 glass (Vogel 2003)", "2003, Table I", ["M. Vogel, arXiv:cond-mat/0307738v2 (2003)"],
                        {"source": "the paper's Table I; bench/ff/make_glass.py", "license": "published parameters of the paper"},
                        "Li, P, O with effective charges and Born-Mayer terms for every pair; LAMMPS reproduces CAPS's energies",
                        ("Glass (LiPO3, Vogel)", "vogel2003", 30), []),
    "xiao2011-caco3": ("CaCO3 (Xiao, Edwards & Grater 2011)", "2011, Table 1",
                       ["S. Xiao, S. A. Edwards, F. Grater, J. Phys. Chem. C 115, 20067 (2011)",
                        "W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, M. L. Klein, J. Chem. Phys. 79, 926 (1983)"],
                       {"source": "the paper's Table 1 and energy expression; bench/ff/make_xiao2011.py", "license": "published parameters of the paper"},
                       "Ca, carbonate C and O with the paper's charges and C12/C6 pairs (Ca-O's C6 0 written as 1e-12), the O-C-O "
                       "bend and the K (1 - cos 2 psi) plane term by distance (no C-O bond term, no exclusions), TIP3P water with "
                       "the paper's mineral-water pairs; LAMMPS reproduces CAPS's energies on calcite and aragonite",
                       ("CaCO3 (Xiao 2011)", "xiao2011", 30), []),
    "gunn2012-pyrochlore": ("Gd pyrochlores (Gunn et al. 2012)", "2012, Table 1",
                            ["D. S. D. Gunn, N. L. Allan, H. Foxhall, J. H. Harding, J. A. Purton, W. Smith, M. J. Stein, I. T. Todorov, "
                             "K. P. Travis, J. Mater. Chem. 22, 4675 (2012)"],
                            {"source": "the paper's Table 1; bench/ff/make_gunn2012.py", "license": "published parameters of the paper"},
                            "rigid Gd, Ti, Zr and O ions with formal charges, Buckingham terms with O, 12 A cut-off; LAMMPS reproduces "
                            "CAPS's energies on Gd2Ti2O7 and Gd2Zr2O7", ("Pyrochlores (Gunn 2012)", "gunn2012", 30), []),
    "la-perovskites": ("La perovskite oxides (Cherry 1995; Read 1999)", "Cherry 1995 and Read 1999, Table 1 each",
                       ["M. Cherry, M. S. Islam, C. R. A. Catlow, J. Solid State Chem. 118, 125 (1995)",
                        "M. S. D. Read, M. S. Islam, F. King, F. E. Hancock, J. Phys. Chem. B 103, 1558 (1999)"],
                       {"source": "the papers' Table 1; bench/ff/make_la_perovskites.py", "license": "published parameters of the papers"},
                       "shell models with formal charges: LaMO3 (M = Cr, Mn, Fe, Co) from Cherry, La2NiO4 and its 2+ dopants from "
                       "Read, each set typed by composition (Ni present or not); LAMMPS reproduces CAPS's energies on LaFeO3, LaCoO3 "
                       "and La2NiO4 (CORESHELL styles)", ("La perovskites", "la_perovskites", 30), []),
    "schroder1992-zeolite": ("Zeolites with bridging hydroxyls (Schroder et al. 1992)", "1992, Table 1",
                             ["K.-P. Schroder, J. Sauer, M. Leslie, C. R. A. Catlow, J. M. Thomas, Chem. Phys. Lett. 188, 320 (1992)"],
                             {"source": "the paper's Table 1 (read from the page image); bench/ff/make_schroder1992.py",
                              "license": "published parameters of the paper"},
                             "formal-charge Si and Al, shell-model framework O, hydroxyl O and H with the paper's fractional charges; the "
                             "O-H Morse as a pair term within 1.4 A with its Coulomb term kept; O-T-O bends k (theta - theta0)^2 by "
                             "distance; LAMMPS reproduces CAPS's energies on quartz and on Si(OH)4 / Al(OH)3",
                             ("Zeolites (Schroder 1992)", "schroder1992", 30), []),
    "ramsahye2005-nay": ("Zeolite NaY and chloroform (Ramsahye & Bell 2005)", "2005, Table 1",
                         ["N. A. Ramsahye, R. G. Bell, J. Phys. Chem. B 109, 4738 (2005)"],
                         {"source": "the paper's Table 1; bench/ff/make_ramsahye2005.py", "license": "published parameters of the paper"},
                         "Kramer-charge framework (Si, Al, O) with Na, Buckingham terms and O-T-O bends by distance; chloroform's "
                         "charges, bonds, bends and its Lennard-Jones pairs with itself, O and Na; LAMMPS reproduces CAPS's energies on "
                         "quartz and on chloroform with Na+, AlO4 and SiO4", ("Zeolite NaY (Ramsahye 2005)", "ramsahye2005", 30), []),
    "forester1997-silicalite": ("Silicalite and benzene (Forester & Smith 1997)", "1997, Tables 1 and 2",
                                ["T. R. Forester, W. Smith, J. Chem. Soc., Faraday Trans. 93, 3249 (1997)"],
                                {"source": "the paper's Tables 1 and 2; bench/ff/make_forester1997.py", "license": "published parameters of the paper"},
                                "Schroder and Sauer's silica shell model as the paper gives it (Si-O shell repulsion only, O-Si-O bends "
                                "on the shells) and the benzene-O shell pairs (rigid benzene); LAMMPS reproduces CAPS's energies on "
                                "quartz and on benzene with SiO4", ("Silicalite (Forester 1997)", "forester1997", 30), []),
    "hill-sauer1995": ("Aluminosilicates (Hill & Sauer 1995)", "1995, Table 8",
                       ["J.-R. Hill, J. Sauer, J. Phys. Chem. 99, 9536 (1995)", "J.-R. Hill, J. Sauer, J. Phys. Chem. 98, 1238 (1994)"],
                       {"source": "the paper's Table 8 (cross terms read from the page image); bench/ff/make_hill_sauer1995.py",
                        "license": "published parameters of the paper"},
                       "the ab initio CFF for silica and zeolite clusters: quartic bonds and angles, torsions, bond-bond, bond-angle, "
                       "angle-angle and angle-angle-torsion cross terms, the bridging O's out-of-plane term, A_i A_j / r^9 repulsion "
                       "and bond-increment charges; LAMMPS reproduces CAPS's energies on disilicic acid, a bridging-hydroxyl Si-O(H)-Al "
                       "cluster, (HO)3Si-O-Al(OH)3 and a cyclotetrasilicic acid", ("Aluminosilicates (Hill & Sauer)", "hill_sauer1995", 30),
                       ["inorganic-zeolite-hill-sauer"]),
    "trappe-eh-umn": ("TraPPE-EH", "the TraPPE database (trappe.oit.umn.edu), explicit-hydrogen family, retrieved October 2026",
                      ["N. Rai, J. I. Siepmann, J. Phys. Chem. B 111, 10790 (2007)", "B. L. Eggimann et al., J. Chem. Eng. Data 65, 1330 (2020)"],
                      {"source": "the TraPPE group's own parameter tables, model by model (trappe.oit.umn.edu), each term with its paper; "
                                 "transcribed by bench/ff/make_trappe_eh.py",
                       "license": "published parameters of the TraPPE papers, cited per term (no file of the site is distributed)"},
                      "the 65 TraPPE-EH models (alkanes, amines, amides, nitriles, nitro compounds, aromatics and heteroaromatics), each "
                      "with its own charges: one CAPS type per symmetry class of a model, typed by a whole-molecule SMARTS (848 / 848 "
                      "sites as the database types them); every model's bonds, bends and torsions as it lists them (59 / 65 models "
                      "exact, the rest the database's symmetric-position differences); LAMMPS and GROMACS reproduce CAPS's energies",
                      ("TraPPE explicit hydrogen", "trappe_eh", 17), ["trappe-eh"]),
    "omalley2013-alkanes": ("All-atom alkanes in silicalite (O'Malley & Catlow 2013)", "2013, Tables 2 and 3 (all-atom model)",
                            ["A. J. O'Malley, C. R. A. Catlow, Phys. Chem. Chem. Phys. 15, 19024 (2013)"],
                            {"source": "the paper's Tables 2 and 3 (read from the page images); bench/ff/make_zeolite_sorbates.py",
                             "license": "published parameters of the paper"},
                            "alkane C typed by its H count (charges +0.1 per H, each alkane neutral), H, and the formal-charge silica "
                            "framework (Buckingham Si-O and O-O within 10 A, O-Si-O bends by distance); the paper's torsion (Raj, Sastre "
                            "& Catlow 1999) not included; LAMMPS reproduces CAPS's energies on alpha-quartz and on alkanes with SiO4, "
                            "GROMACS on the alkanes", ("Alkanes in silicalite (O'Malley 2013)", "omalley2013", 30), ["misc"]),
    "plant2007-methanol": ("Methanol in zeolite NaY (Plant, Maurin & Bell 2007)", "2007, Tables 1-3",
                           ["D. F. Plant, G. Maurin, R. G. Bell, J. Phys. Chem. B 111, 2836 (2007)"],
                           {"source": "the paper's Tables 1-3; bench/ff/make_zeolite_sorbates.py", "license": "published parameters of the paper"},
                           "Kramer-charge NaY framework with Na (as Ramsahye & Bell), methanol's charges, bonds, bends, torsion and "
                           "pairs, and its pairs with the framework O and Na; LAMMPS reproduces CAPS's energies on methanol with Na+, "
                           "AlO4 and SiO4, GROMACS on methanol", ("Methanol in NaY (Plant 2007)", "plant2007", 30), []),
    "peguin2009-hfa134a": ("HFA134a (Peguin et al. 2009, FF no. 2)", "2009, Tables 1 and 2 (FF no. 2)",
                           ["R. P. S. Peguin, G. Kamath, J. J. Potoff, S. R. P. da Rocha, J. Phys. Chem. B 113, 178 (2009)"],
                           {"source": "the paper's Tables 1 (read from the page image) and 2; bench/ff/make_zeolite_sorbates.py",
                            "license": "published parameters of the paper"},
                           "CF3CH2F with the paper's FF no. 2 charges and Lennard-Jones terms, OPLS-AA bends and the fitted F-C-C-F and "
                           "F-C-C-H torsions, no intramolecular pairs; LAMMPS and GROMACS reproduce CAPS's energies",
                           ("HFA134a (Peguin 2009)", "peguin2009", 30), []),
}
RTF = "the topology's own residues typed as it types them (bench/ff/validate_rtf_types.py --cap-links): "
for fid, name, version, refs, evidence, listing, replaces in (
        ("cgenff-toppar", "CHARMM36 CGenFF", "toppar_c36_jul17 (par_all36_cgenff.prm, top_all36_cgenff.rtf: CGenFF 3.0.1 / 4.0)",
         ["K. Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010)"],
         RTF + "top_all36_cgenff.rtf 9872 / 9872 atoms, 519 / 519 residues (residues with lone pairs left out)",
         ("CHARMM36 CGenFF", "charmm36_cgenff", 23), ["cgenff"]),):
    ENTRIES[fid] = (name, version, refs, CHARMM, evidence, listing, replaces)
# CGenFF v. 2a5, the version published with the CGenFF paper, from its Supporting Information (bench/ff/open/build_cgenff_si.sh)
ENTRIES["cgenff-2a5-jcc2010"] = (
    "CGenFF 2a5 (J. Comput. Chem. 2010)", "2a5 (the paper's Supporting Information: par_all36_cgenff.prm, top_all36_cgenff.rtf)",
    ["K. Vanommeslaeghe, E. Hatcher, C. Acharya, S. Kundu, S. Zhong, J. Shim, E. Darian, O. Guvench, P. Lopes, I. Vorobyov, "
     "A. D. MacKerell Jr., J. Comput. Chem. 31, 671 (2010), doi:10.1002/jcc.21367"],
    {"source": "Supporting Information of the paper, jcc_21367_sm_suppinfo.tgz (doi:10.1002/jcc.21367): par_all36_cgenff.prm and "
               "top_all36_cgenff.rtf, CGenFF v. 2a5",
     "license": "parameter values as published in the Supporting Information of Vanommeslaeghe et al. 2010; CAPS's conversion"},
    "CAPS's CGenFF rules (bench/typing/make_cgenff_rules.py) without the types later releases added, with 2a5's phosphorus "
    "and sulfoxide types; the SI's top_all36_cgenff.rtf typed as it types itself for 7471 / 7505 atoms, 371 / 388 residues "
    "(bench/ff/validate_rtf_types.py; the other 34 atoms are 2a5 assignments that later releases revised, as CAPS's rules "
    "follow, or residues later releases dropped); every residue complete at its integer charge with no missing bonded term, "
    "388 / 388 with the topology's types and 386 / 388 with CAPS's (PRAL and GABD, where CAPS's types differ from 2a5's, each "
    "lack one improper, bench/ff/check_rtf_params.py); impropers only where the residues' IMPR records put them (the "
    "library's improper_placement, convert_charmm_prm.py --impropers-from-rtf): the 323 they list, 0 missing, 0 extra; "
    "LAMMPS and GROMACS reproduce CAPS's energies (Urey-Bradley, 1-4 Lennard-Jones)",
    ("CGenFF 2a5 (2010)", "cgenff_2a5", 23), [])
cat_path = os.path.join(FF, "catalogue.json")
cat = json.load(open(cat_path))
fl = cat["forcefields"]
# biomolecular force fields removed on 2026-10-14 (CAPS focuses on polymers, crystals and materials): their entries go
REMOVED = {"amber-ff03-ambertools", "amber-ff14sb-ambertools", "oplsaam-jorgensen", "gromos54a7-gromacs", "martini22-proteins",
           "charmm36-prot-toppar", "charmm22-prot-toppar", "charmm36-carb-toppar", "charmm36-lipid-toppar",
           "charmm36-nucl-toppar", "charmm19-toppar"}
fl[:] = [x for x in fl if x["id"] not in REMOVED]
for fid, (name, version, refs, lic, evidence, listing, replaces) in ENTRIES.items():
    path = os.path.join(FF, fid + ".json")
    e = {"id": fid, "name": name, "version": version, "references": refs, "origin": "built from source", **lic,
         "status": "validated", "file": fid + ".json"}
    if os.path.exists(path):
        d = json.load(open(path))
        e["counts"] = {"atom_types": len(d["atom_types"]), **{k: len(d.get(k, [])) for k in ("pairs", "bonds", "angles", "dihedrals", "impropers")}}
    else:   # a library built only on the user's machine (CHARMM's toppar) and absent here: the counts the catalogue had
        e.update({k: v for x in fl if x["id"] == fid for k, v in x.items() if k == "counts"})
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
