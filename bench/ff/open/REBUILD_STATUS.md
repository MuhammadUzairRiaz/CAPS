# Force-field libraries: provenance and rebuild status

Commit 860a126b stopped tracking 32 libraries (and their typing files) that had been converted from DL_FIELD 4.13's
`.par` / `.sf` files, whose licence forbids redistributing derived data. Each family is being rebuilt from its original
published source (paper, official parameter file or the authors' database). DL_FIELD stays on the developer's machine
as an accuracy reference for typing and energies only; no number in a rebuilt library comes from it.

State of 2026-10-13 (branch `ff-rebuild-rest`, commits 1053f7a2 and b2ee52e2; the inorganic sets of section 4 on branch
`inorg-papers`). "Local" means the source is in `~/project/ff-sources/`.

**2026-10-14, removed by decision: CAPS focuses on polymers, crystals and materials.** The biomolecular libraries
(those that build proteins, nucleic acids, lipids and carbohydrates) left the repository with their typing rules,
generators and checks: amber-ff03-ambertools, amber-ff14sb-ambertools, gromos54a7-gromacs, martini22-proteins,
oplsaam-jorgensen, and the typing of the locally built charmm36-prot / -lipid / -nucl / -carb-toppar,
charmm22-prot-toppar, charmm19-toppar and glycam06j-ambertools. GAFF / GAFF2, CGenFF and the CG models stay.

## 1. The 32 libraries removed in 860a126b

| removed library | was | rebuilt as | original source | local | licence of the source | status |
|---|---|---|---|---|---|---|
| amber | AMBER.par | amber-ff03-ambertools, amber-ff14sb-ambertools | AmberTools 23.6 `dat/leap/parm` (parm99, parm10, frcmod.ff03 / ff14SB / OL15) | yes (`ambertools-23.6`) | public domain | done (other session); removed 2026-10-14 by decision (biomolecular) |
| gaff-amber16 | AMBER16_gaff.par | gaff-1.81-ambertools | AmberTools `gaff.dat` | yes | public domain | done (other session) |
| gaff-amber25 | AMBER25_gaff.par | gaff2-2.2.20-ambertools (and gaff2-2.11) | AmberTools `gaff2.dat`, `gaff211.dat` | yes | public domain | done (other session) |
| charmm | CHARMM.par | charmm36-*-toppar | MacKerell lab `toppar_c36_jul17` | yes (`charmm-toppar-jul17`) | no licence stated: built locally, never committed | done (other session); removed 2026-10-14 by decision (biomolecular) |
| charmm19 | CHARMM19.par | charmm19-toppar | `toppar_c36_jul17` (param19 / toph19) | yes | as above | done (other session); removed 2026-10-14 by decision (biomolecular) |
| charmm22-prot | CHARMM22_prot.par | charmm22-prot-toppar | `toppar_c36_jul17` | yes | as above | done (other session); removed 2026-10-14 by decision (biomolecular) |
| charmm36-carb / -lipid / -nucl / -prot | CHARMM36_*.par | charmm36-carb/-lipid/-nucl/-prot-toppar | `toppar_c36_jul17` | yes | as above | done (other session); removed 2026-10-14 by decision (biomolecular) |
| cgenff | CHARMM36_cgenff.par | cgenff-2a5-jcc2010 (published); cgenff-toppar (local only) | v2a5: Supporting Information of Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010), doi:10.1002/jcc.21367 (`par_all36_cgenff.prm`, `top_all36_cgenff.rtf`); later versions: `toppar_c36_jul17` (CGenFF 4.x) | yes (`cgenff-jcc2010-si`, `charmm-toppar-jul17`) | SI: the paper's published parameter values; toppar: no licence stated, built locally, never committed | v2a5 built from the paper's SI and published (`bench/ff/open/build_cgenff_si.sh`; SI RTF 7471 / 7505 atoms, 371 / 388 residues as it types them; 388 / 388 residues complete; impropers only where the IMPR records put them, 323 / 323 with 0 extra); later versions (jul17 4.x) local only via toppar |
| pcff | PCFF.par | pcff-frc | `pcff.frc` | - | - | other session |
| compass | COMPASS.par | compass-frc | `compass_published.frc`; H. Sun, J. Phys. Chem. B 102, 7338 (1998) | - | - | other session |
| cvff | CVFF.par | cvff-frc | `cvff.frc` | - | - | other session |
| dreiding | DREIDING.par | dreiding1990 | Mayo, Olafson & Goddard, J. Phys. Chem. 94, 8897 (1990) | yes (`papers/j100389a010.pdf`) | published parameters | done (other session) |
| gromos-54a7 | GROMOS_G54A7.par | gromos54a7-gromacs | GROMACS `share/top/gromos54a7.ff` | yes (GROMACS) | LGPL 2.1 | done (other session); removed 2026-10-14 by decision (biomolecular) |
| opls2005 | OPLS2005.par | oplsaa2024 (and oplsaa-gromacs) | Jorgensen's BOSS files, J. Phys. Chem. B 128, 250 (2024) SI; GROMACS `oplsaa.ff` | yes | SI tables; LGPL 2.1 | other session (OPLS) |
| opls2020 | OPLS2020.par | oplsaa2024 | as above | yes | as above | other session (OPLS) |
| opls-aam | OPLS_AAM.par | oplsaam-jorgensen | Jorgensen lab `top_opls_aam.inp`, `par_opls_aam.inp` | yes (`opls-aam`) | CC BY 4.0 | done (other session); removed 2026-10-14 by decision (biomolecular) |
| opls-clp | OPLS_CL_P.par | clandp | Padua's `il.ff` (github.com/paduagroup/clandp) | yes (`clandp-dl` is the reference only) | MIT | other session (CL&P) |
| opls-des | OPLS_DES.par | oplsdes-doherty2018 (the set DL_FIELD's file held); clandp-des (CL&P's DES set) | oplsdes-doherty2018: Doherty & Acevedo, J. Phys. Chem. B 122, 9982 (2018) and its SI (Tables S1-S6, Schemes S1-S2), transcribed; Acevedo's own GROMACS files (github.com/orlandoacevedo/DES, GPL 3, local in `opls-des-dl/DES-master`) are an accuracy reference only, nothing taken from them. clandp-des: Padua's `des.ff` of 2025/12/19 (github.com/paduagroup/clandp) | yes (`papers/jp8b06647.pdf`, `papers/jp8b06647_si_001.pdf`, `papers/jp109914s.pdf` (its ref. 3); `clandp-dl/clandp-master/des.ff`) | published parameters; MIT | done (section 3) |
| trappe-ua | TRAPPE_UA.par | trappe-ua-umn, trappe-ua-foyer | TraPPE database (trappe.oit.umn.edu); foyer `trappe-ua.xml`; the transferable tables of Martin & Siepmann 1998, 1999, Wick, Martin & Siepmann 2000, Chen, Potoff & Siepmann 2001, Stubbs, Potoff & Siepmann 2004 and Keasler et al. 2012 | yes (`trappe-site`, `foyer-trappe`, `papers/jp972543+.pdf`, `jp984742e.pdf`, `jp001044x.pdf`, `jp003882x.pdf`, `jp049459w.pdf`, `jp302975c.pdf`) | published parameters; MIT (foyer) | done (other session); transferable terms **added here** |
| **trappe-eh** | TRAPPE_EH.par | **trappe-eh-umn** | TraPPE database, EH family (Rai & Siepmann, J. Phys. Chem. B 111, 10790 (2007) and the later EH papers) | yes (`trappe-site`) | published parameters, cited per term | **done here** |
| **misc** | MISC_FF.par | **omalley2013-alkanes, plant2007-methanol, ramsahye2005-nay (chloroform), peguin2009-hfa134a** | the four papers (see section 3) | yes (`papers/c3cp52653d.pdf`, `jp0674524.pdf`, `jp046958o.pdf`, `jp806213w.pdf`) | published parameters | **done here** (one torsion paper missing) |
| inorganic-binary-halides | INORGANIC_binary_halides.par | catlow1977-halides | Catlow, Diller & Norgett, J. Phys. C 10, 1395 (1977) | yes | published parameters | done |
| inorganic-binary-misc | INORGANIC_binary_misc.par | gan2020 | Wonglakhon & Zahn, J. Phys.: Condens. Matter 32, 205401 (2020) | yes | CC BY | done |
| inorganic-binary-oxides | INORGANIC_binary_oxides.par | bush1994-oxides, lewis-catlow1985-oxides, catlow1977-uo2, matsui1991-tio2, sanders1984-sio2, **catlow1976-mgo, catlow1985-rutile, khan1998-zro2, sayle2000-alkaline-earth-oxides, bandura2003-tio2-water** | the papers of DL_FIELD's refs 1-11 | yes | published parameters | done but Bi2O3: Walther & Zahn 2018 prints no parameters (section 4) |
| inorganic-clay | INORGANIC_clay.par | clayff2004 | Cygan, Liang & Kalinichev, J. Phys. Chem. B 108, 1255 (2004) | yes | published parameters | done (other session) |
| inorganic-glass | INORGANIC_glass.par | teter-glass, tilocca2006-sm, vogel2003-lipo3 | Cormack 2002, Du 2009, Okhotnikov 2013, Tilocca 2006, Vogel 2003 | yes | published parameters | done |
| inorganic-ternary-oxides | INORGANIC_ternary_oxides.par | xiao2011-caco3, gunn2012-pyrochlore, la-perovskites, **sepliarsky2005-perovskites** | Xiao 2011, Gunn 2012, Cherry 1995, Read 1999, Sepliarsky 2005 | yes | published parameters | done but PbTiO3: its Pb-Ti row is misprinted (section 4) |
| inorganic-zeolite | INORGANIC_zeolite.par | forester1997-silicalite, schroder1992-zeolite, ramsahye2005-nay, sanders1984-sio2 | the four papers | yes | published parameters | done |
| inorganic-zeolite-hill-sauer | INORGANIC_zeolite_Hill_Sauer.par | hill-sauer1995 | Hill & Sauer, J. Phys. Chem. 98, 1238 (1994); 99, 9536 (1995) | yes | published parameters | done |

The typing overlays and water-model files removed with them (`data/typing/*-tip3p*.typing.json`, `*-spc.typing.json`,
`data/typing/sets/*`) go with their families; the inorganic potential sets of `data/typing/sets/` are replaced by the
per-paper libraries above.

## 2. Libraries CAPS ships (tracked in `data/forcefields/`)

| library | provenance | licence of the source | status |
|---|---|---|---|
| pcff-frc, cvff-frc, compass-frc | BIOVIA `.frc` files (`pcff.frc`, `cvff.frc`, `compass_published.frc`) | handled by the class II session | other session |
| iff-pcff, iff-cvff | INTERFACE FF 1.5 `.frc` files (Heinz et al., Langmuir 29, 1754 (2013)) | as distributed by the Heinz group | original source |
| gaff-1.81-ambertools, gaff2-2.11-ambertools, gaff2-2.2.20-ambertools | AmberTools 23.6 `dat/leap/parm` | public domain | original source |
| oplsaa-gromacs | GROMACS `oplsaa.ff` | LGPL 2.1 | original source |
| oplsaa2024, oplsua2024 | Jorgensen's BOSS files (J. Phys. Chem. B 128, 250 (2024) SI) | published SI tables | original source (OPLS session) |
| **oplsua2024 (1984, 1986, 1990 terms)** | the hydrocarbon, alcohol and ether terms the 2024 file lacks, from Jorgensen, Madura & Swenson, J. Am. Chem. Soc. 106, 6638 (1984), Tables I and II; Jorgensen, J. Phys. Chem. 90, 1276 (1986), Tables II and III; Briggs, Matsui & Jorgensen, J. Comput. Chem. 11, 958 (1990), Table II and p. 960 (`bench/ff/open/oplsua2024_supplement.json`) | published parameters | original source |
| clandp | Padua's `il.ff` | MIT | original source |
| clandp-des | Padua's `des.ff` | MIT | original source |
| oplsdes-doherty2018 | Doherty & Acevedo 2018, the paper's SI tables | published parameters | original source (Acevedo's GPL 3 files are the reference only) |
| dreiding1990 | the 1990 paper | published parameters | original source |
| trappe-ua-umn, trappe-ua-foyer | TraPPE database; foyer `trappe-ua.xml` | published parameters; MIT | original source |
| **trappe-eh-umn** | TraPPE database, EH models | published parameters | original source (this branch) |
| **omalley2013-alkanes, plant2007-methanol, peguin2009-hfa134a** | the papers' tables | published parameters | original source (this branch) |
| bush1994-oxides | Bush et al. 1994 as the UCL database of published potentials transcribes it | published parameters | original paper's set, transcribed via the UCL database |
| catlow1977-halides, catlow1977-uo2, lewis-catlow1985-oxides, matsui1991-tio2, sanders1984-sio2, catlow1976-mgo, catlow1985-rutile, khan1998-zro2, sayle2000-alkaline-earth-oxides, sepliarsky2005-perovskites, bandura2003-tio2-water, gan2020, clayff2004, teter-glass, tilocca2006-sm, vogel2003-lipo3, xiao2011-caco3, gunn2012-pyrochlore, la-perovskites, schroder1992-zeolite, ramsahye2005-nay, forester1997-silicalite, hill-sauer1995 | each paper's tables | published parameters | original source |
| *-moltemplate (COMPASS published subset, GAFF, GAFF2, OPLS-AA 2008 / 2024, L-OPLS, DREIDING, TraPPE 1998, SDK, MARTINI 2.0 and its .prm sets, Dry MARTINI, Cooke-Deserno, graphene, SPC/E, TIP3P, mW and the water templates) | moltemplate's force-field files | MIT (moltemplate) | third-party conversion of the published sets |
| martini3 | `martini_v3.0.0.itp` (cgmartini.nl) | as distributed by the Martini group | original distribution |
| cgenff-toppar | CHARMM toppar (no licence) | - | built locally, never committed (`.gitignore`) |

## 3. Families rebuilt on this branch

| family | library | generator | typing against DL_FIELD | LAMMPS | GROMACS | ctest |
|---|---|---|---|---|---|---|
| MISC: alkanes in silicalite | omalley2013-alkanes | `bench/ff/make_zeolite_sorbates.py` | DL_FIELD's `alkanes1` is a hand-labelled template (Cp1 / Cs1 / Ct1 / H1 by the number of H) that DL_FIELD cannot type automatically; CAPS types CH3 / CH2 / CH / C / CH4 and H by the same rule | agrees (alpha-quartz; alkanes with SiO4) | agrees (alkanes) | `FFDef.ZeoliteSorbatesFromTheirPapers` |
| MISC: methanol in NaY | plant2007-methanol | same | `methanol2` 6 / 6 atoms | agrees (methanol, Na+, AlO4, SiO4) | agrees | same |
| MISC: chloroform in NaY | ramsahye2005-nay (existing) | `make_ramsahye2005.py` | `chloroform3` 5 / 5 atoms | agrees (existing case) | - | existing |
| MISC: HFA134a | peguin2009-hfa134a | `make_zeolite_sorbates.py` | `HFA_134a` 8 / 8 atoms | agrees | agrees | same |
| TraPPE-EH | trappe-eh-umn | `bench/ff/make_trappe_eh.py SITE_DIR --check` | 9 templates, 94 / 94 atoms (17 carry the key of a symmetry-equivalent position); the database's own 65 models 848 / 848 sites | agrees (flexible and rigid sets) | agrees | `FFDef.TrappeEhModelsTypedAsThemselves` |
| TraPPE-UA: transferable hydrocarbon terms | trappe-ua-umn | `bench/ff/make_trappe.py SITE_DIR --check`: 31 bonds, 132 bends, 213 torsions over CH3 / CH2 / CH / C (sp3), CH2= / CH= / C= and aromatic CH / R-C, each citing its table; every one equals the database models' term for the same types (21 bond, 32 bend, 46 torsion parameter sets; 6 rigid-model K 0 terms kept) | not compared (the database's own models: 634 / 634 sites, 92 / 92 models' bonded terms, ten models' own lists named apart); polymer library DP 4 (`bench/ff/trappe_polymer_coverage.py`): 18 / 114 complete and neutral, 4 before | agrees: polypropylene, cis- and trans-1,4-polybutadiene, polystyrene, polyisobutylene DP 4, energy terms 1.1e-7 to 3.7e-7, forces <= 7.5e-9 kcal/mol/A (exact styles); native styles 1.1e-7 to 3.7e-7 | agrees: the same five chains, one per periodic cell, 4.8e-6 to 3.1e-5 | `FFDef.TrappeUaTransferableHydrocarbonTerms` |
| TraPPE-UA: alcohol, ether, glycol and ring terms | trappe-ua-umn | `bench/ff/make_trappe.py SITE_DIR --check`: with the hydrocarbon terms 109 bonds, 845 bends, 2776 torsions; the new ones over the hydroxyl and ether O, H and the carbons bonded to them (Chen, Potoff & Siepmann 2001 Table 1; Stubbs, Potoff & Siepmann 2004 Table 1) and the five- and six-membered cycloalkane rings (Keasler et al. 2012 Tables 2, 3); every one equals the database models' term for the same types (45 bond, 68 bend, 87 torsion parameter sets; the one bend that differs is 2-hydroxyethyl acrylate's own slip, already listed) | not compared (the database's own models: 634 / 634 sites, 92 / 92 models' bonded terms); polymer library DP 4: 30 / 114 complete and neutral, 18 before | agrees: PVA (1-4 Coulomb x 0.5), PVME, poly(isobutyl vinyl ether), poly(vinyl cyclohexane), poly(vinyl cyclopentane) DP 4, energy terms 8.8e-8 to 4.7e-7, forces <= 1.5e-5 kcal/mol/A (exact styles); native styles bonded and vdW 8.8e-8 to 4.7e-7 | agrees: PVA and PVME with PME, PVCH, PVCP, 1.9e-5 to 4.3e-5 | `FFDef.TrappeUaOxygenAndRingTerms` |
| CL&P deep-eutectic solvents | clandp-des | `bench/ff/convert_fftool.py des.ff` (47 types, 16 bonds, 30 angles, 44 dihedrals, 6 impropers), `bench/typing/make_clandp_des_rules.py` (45 rules) | not comparable: DL_FIELD's OPLS_DES is Acevedo's set. Against Padua's own files (`validate_fftool_types.py`): the 14 whose types des.ff defines (cholinium.xyz, N1110, N1111, N2220, N2222, N4444, P1111, P66614, c1c1pyrr, c3c1pyrr, c4c1pyrr, Cl, Br, oac) 362 / 362 atoms, each at its integer charge; 29 DES components typed whole and neutral or at their ion charge, 7 that des.ff cannot make neutral left untyped (`check_clandp_des_components.py`) | agrees: reline and glyceline boxes (608 and 600 atoms, periodic) and a vacuum set of 10 components, energy terms 1.0e-7 to 3.8e-7, forces <= 1.9e-5 kcal/mol/A, virial <= 2.2e-6 (exact styles); native styles bonded and vdW 1.6e-8 to 3.8e-7 | agrees: reline and glyceline with PME 6.9e-6 and 3.3e-5, vacuum set 5.0e-6 | `FFDef.ClandpDesComponentsTypedWhole`, python_bindings |
| OPLS-UA hydrocarbons, alcohols and ethers (Jorgensen et al. 1984, 1986, 1990) | oplsua2024 | `bench/ff/make_oplsua2024.py` with `bench/ff/open/oplsua2024_supplement.json` (56 bends, 595 torsions added where the 2024 file has none, each citing its table; 33 / 420 for hydrocarbons, 23 / 175 for alcohols and ethers) | typing unchanged (the 2024 file's 29 united-atom types, 44 rules); polymer library, DP-4 chains (`polymer_coverage.py data/forcefields/oplsua2024.json`): 24 / 114 complete and neutral (was 1, polyethylene; 22 after the hydrocarbon terms), every hydrocarbon of the library, poly(vinyl alcohol) and polyoxytetramethylene; poly(ethylene oxide), poly(methylene oxide) and two epoxy resins lack O-C-C-O or acetal terms that no paper gives; 86 untyped (no united-atom type in the 2024 file) | agrees: polypropylene tetramer, 1,4-polybutadiene and cis-1,4-polyisoprene trimers, poly(vinyl alcohol) tetramer, polyoxytetramethylene trimer, ethanol / 1-propanol / 2-propanol / tert-butanol / diethyl ether; energy terms 1.7e-7 to 5.0e-7, forces <= 1.4e-5 kcal/mol/A, virial <= 4.9e-7, read back exact | agrees: the same six, 3.9e-6 to 4.2e-5 | `FFDef.OplsUa1984HydrocarbonTerms`, `FFDef.OplsUaAlcoholEtherTerms` |
| OPLS-DES (Doherty & Acevedo 2018) | oplsdes-doherty2018 | `bench/ff/make_opls_des_doherty2018.py` (54 types, 48 bonds, 68 angles, 66 dihedrals, 5 impropers, every term citing its SI table), `bench/typing/make_oplsdes_doherty2018_rules.py` (53 rules) | the SI's Schemes S1-S2: the paper's 9 molecules 102 / 102 atoms by parameters, 100 by name (urea's cis H, HC, typed HT: identical terms), choline +0.8, chloride -0.8, HBDs 0, no missing term; 12 other molecules left untyped (`check_oplsdes_doherty2018_components.py`). Against DL_FIELD's OPLS_DES templates and Acevedo's GROMACS files (references only, `compare_oplsdes_doherty2018.py`): see below | agrees: reline (16 ChCl + 32 urea) and glyceline (12 + 24) in 34 A boxes, a vacuum set of the other HBDs; energy terms 2.2e-8 to 6.4e-8, forces <= 9.9e-6 kcal/mol/A, virial <= 7.8e-7 (exact styles); native bonded and vdW the same | agrees: reline (48 + 96) and glyceline (24 + 48) in 38 A boxes with PME, 1.4e-6 and 4.0e-5; vacuum set 4.8e-5 | `FFDef.OplsDesDoherty2018FromThePaper` |
| binary oxides: MgO (Catlow, Faux & Norgett 1976) | catlow1976-mgo | `bench/ff/make_catlow1976_mgo.py` | `magnesium_oxide3` 2 / 2 atoms | agrees (rock-salt MgO, model 1 and the model 2 overlay) | - (Buckingham) | `FFDef.InorganicSetsFromTheirPapers` |
| binary oxides: rutile (Catlow, Freeman & Royle 1985) | catlow1985-rutile | `bench/ff/make_catlow1985_rutile.py` | `rutile1` 6 / 6 atoms | agrees (rutile, Ti and O shells) | - | same |
| binary oxides: ZrO2 (Khan, Islam & Bates 1998) | khan1998-zro2 | `bench/ff/make_khan1998_zro2.py` | `zirconium(iv)_oxide1` 6 / 6 atoms | agrees (cubic ZrO2) | - | same |
| binary oxides: MgO, CaO, SrO, BaO (Sayle & Watson 2000) | sayle2000-alkaline-earth-oxides | `bench/ff/make_sayle2000.py` | `magnesium_oxide4`, `calcium_oxide2`, `strontium_oxide2`, `barium_oxide2` 8 / 8 atoms | agrees (MgO, BaO) | - | same |
| binary oxides: TiO2 surfaces and water (Bandura & Kubicki 2003) | bandura2003-tio2-water | `bench/ff/make_bandura2003.py` (`--test-slab` writes `tests/data/rutile110_hydroxylated.pdb`) | DL_FIELD's `TiO2_hydroxylated1` template is empty: nothing to compare | agrees (hydroxylated rutile (110) slab with water) | - | `FFDef.BanduraKubickiSurfaceSites` |
| ternary oxides: K(Nb,Ta)O3, (Ba,Sr)TiO3 (Sepliarsky et al. 2005) | sepliarsky2005-perovskites | `bench/ff/make_sepliarsky2005.py` | `barium_titanate1`, `strontium_titanate1`, `potassium_tantalate1`, `potassium_niobate1` 40 / 40 atoms (and `lead_titanate1` 10 / 10 with the withheld PbTiO3 set) | agrees (cubic KNbO3, KTaO3, BaTiO3, SrTiO3) | - | `FFDef.InorganicSetsFromTheirPapers` |

Typing commands: `python3 bench/ff/validate_typing.py MISC_FF data/forcefields/ID.json data/typing/ID.typing.json --map ID --only TEMPLATE`
(INORGANIC_binary_oxides / INORGANIC_ternary_oxides for the oxides; a rebuilt shell model's shells are named through the map)
and `... TRAPPE_EH data/forcefields/trappe-eh-umn.json data/typing/trappe-eh-umn.typing.json --map trappe-eh-umn`; the
type maps are in `bench/typing/reference_type_maps.json` (the TraPPE-EH map is derived by matching each DL_FIELD
template's graph onto the database model, `make_trappe_eh.py --dl-map`, without CAPS's typing).

Notes on the rebuilt sets:

- O'Malley & Catlow 2013 print the C-H length as 1.054 A and the bond K in kJ mol-1 rad-2 (read as kJ mol-1 A-2 with
  DL_POLY's 1/2 K); the all-atom torsion is in their ref. 18, not in the paper: the library has a zero torsion until
  that paper is read. The united-atom model of their Table 1 is not included.
- Plant, Maurin & Bell 2007 share Ramsahye & Bell's framework (Table 1); every pair no table lists is Coulomb only.
- Peguin et al. 2009: Table 1's four torsion columns are read as C0..C3 of their eq. 3; eq. 1 has no intramolecular
  non-bonded energy, so all intramolecular pairs are excluded.
- TraPPE-EH: one CAPS type per symmetry class of each model (charges are not transferable). The database's slips are
  corrected and listed in the library's notes: 2-chlorothiophene's ring S carries furan's O type; ethane, propane,
  triethylamine and imidazole list bends or torsions with their atoms out of order; 17 entries that are no bonded path
  are left out (their symmetry partners carry the terms). Where symmetric sites have different charges (phenol, 1,4-
  dihydroxybenzene, 1,3,5-trichlorobenzene, isoindole) the class takes the mean (largest change 0.007 e).
- TraPPE-UA transferable terms (Martin & Siepmann, J. Phys. Chem. B 102, 2569 (1998), section 2; 103, 4508 (1999), Table 2;
  Wick, Martin & Siepmann, J. Phys. Chem. B 104, 8008 (2000), Tables 2-4): fixed bonds 1.54 / 1.33 / 1.40 A (stiff
  harmonic, as the library writes them), bends at CH2 / CH / C / CH= / C= and the rigid ring bends, torsions
  CHx-CH2-CH2-CHy, CHx-CH2-CH-CHy, CHx-CH-CH-CHy, CHx-CH2-C-CHy, CHx-CH2-CH=CHy, cis / trans CHx-CH=CH-CHy and the two
  aromatic-alkyl torsions. They come first in each list, so a database model's own term wins; none differs. Readings and
  disagreements: the 1999 Table 2 prints c1 335.03 K for CHx-CH2-CH2-CHy where the 1998 paper (eq 5), the 2000 Table 4 and
  the database give 355.03 K (used); the 2000 Table 4 prints d0 12400 / 13400 K for the cis / trans term d0/2 (phi - phi0)^2,
  the database's cis- and trans-2-butene 24800 / 26800 K in the same formula (the database's values used, so a polymer and
  2-butene get one term); the 2000 Table 4 names the second aromatic torsion CHx(aro)-C(aro)-CH(sp2)-CHx, which the
  database applies to isopropylbenzene's sp3 CH (so does CAPS). 2-ethylhexyl acrylate's model leaves out the torsion
  CH2-CH2-CH-CH2 (12-7-8-9) that the 1999 table defines: it was a zero term before and is the paper's now (the model is
  listed with the database's slips, hence 92 / 92). sp3 and sp2 carbons in three- and four-membered rings are no longer
  typed (no TraPPE-UA bends for them); no database model has one. Still incomplete on the polymer library for want of
  terms in these papers: polyisoprene (cis and trans) and poly(2-methylpentadiene) (torsions about CH=C(CH3) and CH2-C=),
  1,2-polybutadiene (CHx-CH-CH=CH2), the cycloalkyl polymers (ring CH2 types of the 2012 cyclic-alkane paper bonded to
  CH), and every polymer with chemistry outside TraPPE-UA (halogens, esters on a backbone CH, amides, ...). Poly(4-
  methylstyrene) is complete but its ring methyl takes the rigid toluene / xylene terms (K 0): run the ring with its
  methyl as a rigid body.
- TraPPE-UA alcohol, ether and ring terms (Chen, Potoff & Siepmann, J. Phys. Chem. B 105, 3093 (2001), Table 1; Stubbs,
  Potoff & Siepmann, J. Phys. Chem. B 108, 17596 (2004), Table 1; Keasler, Charan, Wick, Economou & Siepmann, J. Phys.
  Chem. B 116, 11234 (2012), Tables 2 and 3): bonds C-OH 1.43, O-H 0.945, C-O (ether) 1.41 A; bends C-C-C at the O-bonded
  carbons, C-C-OH 109.47 / 50400 K, C-C-O 112 / 50300 K, C-O-H 108.5 / 55400 K, C-O-C 112 / 60400 K; torsions CHx-CH2-CH2-O,
  O-CH2-CH2-O, CHx-CH2/CH/C-O-H, CHx-CHy-O-CHz; in five- and six-membered rings C-C-C 105.5 / 114 deg and the ring C-C-C-C
  torsions. No new Lennard-Jones or charge sites: the papers' sites are the database's types already. Readings: a 1999 C-C
  torsion with one O end takes the carbon term (the database's 2-butanol does); a ring CH or C carrying a substituent is
  typed as the linear CH / C (the 2012 paper fits CH2 rings), with the ring bend and torsions inside the ring and the linear
  terms for an end outside it; 1-4 Coulomb x 0.5 (the 2004 rule) now also for a chain with two O three carbons apart (PVA,
  the vinyl ethers). Disagreements: CHx-(CH)-(O)-H c0 215.96 K in the 2001 table, 215.89 K in the database's propan-2-ol
  and 2-butanol (a constant; the database's value used); 2-hydroxyethyl acrylate's CH2-CH2-OH bend (listed before). The
  2012 torsions are a series in psi = phi - 180 (their minima are the rings' own torsions; read in phi they put 15 kcal/mol
  into each cyclohexane torsion): the library had read them in phi for every ring model of the database (cyclopentane,
  cyclohexane, THF, oxane, the dioxanes, 1,3-dioxolane, 1,3,5-trioxane) and is corrected. Typing: ring CH2 types only for a
  CH2 in one ring; the linear CH2 never in a ring; a CH or C shared by two rings untyped (norbornane, decalin). Still
  incomplete: poly(methylene oxide) and poly(acetaldehyde) (O-C-O-C torsions, acetal bends), the glycidyl and alkylene
  oxide polyethers (O-CH2-CH-O), epoxidised isoprene (CHx-C-CH-CHy), ENB and vinylcyclohexene dioxide (sp2 or O-bearing ring
  carbons), the aryl ethers (no aromatic ether O in TraPPE-UA). The 2004 paper's intramolecular O...H 1-5 repulsion (1,2-diols
  and hydroxy ethers) is not written: none of the polymers here has such a pair.
- A core fix came with the O'Malley framework: a pair term with its own cut-off (Buckingham within 10 A under a 13 A
  global cut-off) no longer adds minus its shift beyond that cut-off (`core/src/field.cpp`).

Notes on the inorganic sets (terms against DL_FIELD's .par: `python3 bench/ff/compare_inorganic_par.py ID ...`):

- catlow1976-mgo: Table 1(b), 1(c). Mg-O and O-O equal DL_FIELD's `Mg3` / `O10` terms; the paper's model is a shell
  model with only the anion polarisable (Y- -2.62, K- 38.03 for model 1), DL_FIELD's `magnesium_oxide3` a rigid ion
  (O -2.0, no shell): the paper's model is written. At the perfect lattice (shells on their cores) both give the same
  energy, -40.82 eV per MgO at r0 = 2.106 A (the paper's cohesive energy 40.85 eV). Model 2 is an overlay; the breathing
  shell model has no CAPS form and is not written.
- catlow1985-rutile: Table IVa. Pairs, charges and harmonic springs equal DL_FIELD's `rutile1`; the footnote's quartic
  constant 1.1e5 (printed in eV/A^5) DL_FIELD keeps only as a comment (k4 0), CAPS writes it as K w^4/24 with K in
  eV/A^4: the paper gives no expression, the 1/24 is the anharmonic-shell convention of the same group's codes and of
  Sepliarsky's Table 1. Only shell displacements feel it (none at the perfect lattice). The paper's alpha-quartz set
  (Table Ia) is sanders1984-sio2, every number equal; DL_FIELD's `silicon_dioxide1` is that set.
- khan1998-zro2: Table 1; every term and charge equal DL_FIELD's (8 / 8). Cubic ZrO2 at a = 5.07 A: -109.76 eV per
  ZrO2, the paper's lattice energy (Table 2).
- sayle2000-alkaline-earth-oxides: Table 1; every term equals DL_FIELD's (10 / 10). Formal-charge rigid ions.
- sepliarsky2005-perovskites: Table 1 (signs from the rendered table; the PDF text layer drops the minus signs).
  K(Nb,Ta)O3 and (Ba,Sr)TiO3: 42 terms equal DL_FIELD's, charges and springs included (k4 as c4 of c4 w^4/24, as
  DL_FIELD's header writes it). Each cubic cell has its energy minimum near its lattice constant (BaTiO3 3.95-4.0 A,
  SrTiO3 3.905 A, KNbO3 3.97 A, KTaO3 3.99 A). PbTiO3 is withheld: section 2 states its Pb-Ti, Pb-O and Ti-O pairs
  are Rydberg terms (A + C r) exp(-r/rho). Read that way, Pb-O and Ti-O behave (without Pb-Ti, cubic PbTiO3 has its
  minimum near a = 3.85-3.9 A, the LDA lattice constant the set was fitted to), but the Pb-Ti row (0.09 + 12.57 r)
  exp(-r/2.4201) is +10.5 eV per pair at 3.4 A and makes the crystal unbound (+141 eV per formula unit at a = 3.97 A).
  DL_FIELD writes all three rows as Buckingham A exp(-r/rho) - C/r^6; then there is no minimum between 3.8 and 4.15 A,
  so that reading is wrong too. CAPS now evaluates Rydberg pairs (`rydberg` pair style, pair_func form 3; LAMMPS as
  `table spline` sections in caps_rydberg.table, which reproduced CAPS's energies and forces on tetragonal PbTiO3 to
  9e-6 / 2e-5 kcal/mol/A during the check), so the set can be written (`SEPLIARSKY_SETS=knt,bst,pt`) once the Pb-Ti
  row is read from the original papers (section 4).
- bandura2003-tio2-water: Tables 2, 3, 5, 6, 9, 10 (kcal/mol). DL_FIELD has no template and no terms for it. Surface
  sites by neighbours, CAPS's reading of Table 9's names: TiV a Ti with five oxide O (a water O on it not counted),
  TiVIoh a Ti holding a terminal OH, TiVI a Ti holding a bridging O, Os a 3-fold O between a TiV / TiVIoh and a TiVI
  (the in-plane O of the surface Ti layer); Ou keeps the bulk O's charge and terms. Table 10's Ti^c and O^d ("any
  type") are every Ti and every oxide O. Ti-O bonds carry no term and are dropped after typing (bond perception also
  links Ti to Ti 2.96 A apart, so the counts are of O neighbours, by recursive SMARTS); the Ti-O-H bends take the Ti
  within 2.6 A of the hydroxyl O; the Morse constant -A is dropped. Charges are Table 9's, so a hydroxylated slab is not
  exactly neutral (the paper's own scheme).

**OPLS-UA hydrocarbon terms (oplsua2024).** The 2024 OPLS parameter file keeps the united-atom types of the 1984
paper (66-77) but only some of their bonded terms: its UA torsions are the later "UA alkanes" line (C2-C2-C2-C2,
C2-CH-C2-C2 ...) and its jtr 5-14-91 hydrocarbon block, so polypropylene's CH-CH2 torsions, everything about a C=C in
polybutadiene and the polyisoprenes, and the bends at an sp3 carbon next to an sp2 one were missing. They come from
W. L. Jorgensen, J. D. Madura, C. J. Swenson, J. Am. Chem. Soc. 106, 6638 (1984), through a supplement the generator
expands over the classes each entry names, adding a term only where the 2024 file has none (the 2024 file wins; its own
terms are unchanged, polyethylene's energy is the same as before):
- Table II, isopentane (V1 1.526, V2 0.533, V3 -3.453, one function of the CH-CH2 rotation, phi = 0 with the CH2's
  carbon gauche to both carbons on the CH): written per torsion as V1 1.526, V2 -0.533, V3 1.7265 on each of the two
  C-CH-CH2-C torsions (they sit at phi +- 60, so V1' = V1, V2' = -V2, V3' = -V3/2 give the paper's function up to a
  constant), for every X-CH-C2-Y.
- Table II, 1-butene (V1 0.343, V2 -0.436, V3 -1.121) for C=C-C(sp3)-X; Table II, other n-alkanes (1.411, -0.271,
  3.145) for the remaining sp3-sp3 torsions. The paper's V0 are constants and left out.
- Table I gives bond angles only (the paper holds them fixed): C-C4-C 112.0 and C-C3-C 124.0 degrees, with the 2024
  file's force constants for UA bends at sp3 (63.0) and sp2 carbons (70.0). The 2024 block puts the C=C-C bends it
  lists at 118 degrees; those stay, so a polyisoprene has 118 at its C= and 124 at its CH=.
- Not in the paper, recorded as such in the supplement: torsions about C=C (the paper simulates the 2-butenes as rigid
  isomers) take the 2024 file's alkene line ??-CM-CM-?? (V2 14.0 per torsion); the torsion from a C= carbon's
  single-bonded substituent (isoprene's CH3-C(=)-CH2-CH2) is zero, as the 1-butene function carries the whole rotation
  on the C=C-C-C torsion.
- The paper fitted its torsions with fixed bonds and angles, no 1-4 Lennard-Jones and a weak C-C term beyond 1-4
  (sigma 4 A, epsilon 0.0074 kcal/mol); here they run with the 2024 file's flexible bonds and 1-4 x 0.5.
- `check_data_lammps.py` now takes an overlay's cut-off from the file it extends (the read-back check used 10 A for
  OPLS-UA's 12 A).

Alcohols and ethers, the same way, from W. L. Jorgensen, J. Phys. Chem. 90, 1276 (1986) and J. M. Briggs, T. Matsui,
W. L. Jorgensen, J. Comput. Chem. 11, 958 (1990), the papers of types 78-81, 106-110 (23 bends, 175 torsions):
- 1986 Table III: ethanol's C-O (0.834, -0.116, 0.747) for X-CH2-OH-H; 1-propanol's C1-C2 (0.702, -0.212, 3.060) for
  C-C-C-O, taken for ether O too (the 1990 paper turns only C-O bonds; the 2024 file gives OH and OS one C-C-C-O? line);
  2-propanol's C-O (V0 0.429, V1 0.784, V2 0.125, V3 -0.691, trans = hydroxyl H between the methyls) written per torsion
  as 0.784, -0.125, 0.3455 on each C-CH-O-H (the torsions sit at phi +- 60, as for isopentane); 2-methyl-2-propanol's
  V3 0.650 as 0.650 / 3 on each of its three C-C-O-H torsions.
- 1990 p. 960 (eq 2, EME and DEE): C-C-O-C (2.8827, -0.6507, 2.2183), one torsion per side so per torsion as printed.
- Bends: the papers' fixed angles (1986 Table II CCO 108.0, COH 108.5; 1990 Table II COC 112.0, OCC 112.0) with the
  2024 file's force constants (80, 55, 100, 80). The 2024 file's own C-C-O bends (109.5) stay.
- Not added: the 1-5 Lennard-Jones terms both papers fit with their torsions (1986: HO-CH3 sigma 2.6 A, epsilon 0.008;
  1990 eq 3: sigma 2.745 A, epsilon 0.008 kcal/mol). No paper has a torsion about a C-C bond between two oxygens
  (O-C-C-O of poly(ethylene oxide) and glycols) or acetal O-C-O terms; those stay missing. The 2024 file's all-atom
  OS-CT-CT-OS (V1 -0.55 only) works with H-C-C-O terms that a united-atom chain does not have, so it is not carried over.
- The 2024 file has no united-atom CH on an ether O; vinyl ethers stay untyped (the 1990 paper defines only CH3 and
  CH2 next to O).

**OPLS-DES (oplsdes-doherty2018).** Built from B. Doherty, O. Acevedo, J. Phys. Chem. B 122, 9982 (2018) and its SI
only. The paper: OPLS-AA form, geometric mixing, 1-4 pairs x 0.5, charges built for +-0.8 e ions (the SI prints them
scaled), choline chloride's charges and Lennard-Jones refitted to reline's ND/EPSR structure, choline's bonded terms from
OPLS-2009IL, the HBDs' bonded terms OPLS-AA's (ethylene glycol's torsions Szefczyk & Cordeiro's), the HBDs' non-bonded
terms tuned for DES mixtures only; GROMACS with a 16 A cut-off. The SI prints every term (Tables S1-S6), so no term
comes from another library. Readings: ethylene glycol's rows in OPLS-AA class names (CT, HC) and its O called OH are
Table S4's CG, HG, OG; levulinic acid's bend "HO-CB-CZ" is OH-CB-CZ; the X-c-X-X rows are impropers on c (centre second,
every triple of neighbours); the all-zero X-CM-X-X and X-CT-X-X are left out; ethylene glycol's block lists only
OG-CG-CG-OG and CG-CG-OG-HO, so its H-C-C-H, H-C-C-O and H-C-O-H torsions are zero (both references have none either,
and the SI's ref. 3, Szefczyk & Cordeiro, J. Phys. Chem. B 115, 3013 (2011), confirms it: its Table 2 holds only
OG-CG-CG-OG and CG-CG-OG-HO, and p. 3016 says the model has three torsions, one O-C-C-O and two C-C-O-H, fitted to
RHF/6-31G(d) scans; its values, 2.993 / -0.654 / 2.820 and 0.396 / -0.750 / 1.028 kcal/mol, are not the SI's rows,
3.887 / -1.192 / 3.206 and 0.413 / -0.754 / 1.028, and the SI's are kept as the set the paper ran);
urea's HC and HT carry identical terms and CAPS types both HT. The type-charge check takes the library's
`ion_charge_scale` 0.8 (new FFDef field) so the scaled ions count as whole.

Against the two references (on the reference's own geometry nudged by 0.08 A, `compare_oplsdes_doherty2018.py DES_DIR --dlfield`):
- DL_FIELD's OPLS_DES (its templates for the 8 molecules, 101 atoms): charges, Lennard-Jones, bonds, bends and torsions
  all equal; energies equal on ethylene glycol, glycerol and phenol. It differs only in impropers: none on choline's N
  (the SI has X-NA-X-X, V2 2.0: 6.08 kcal/mol on the test geometry), two carbonyl impropers on levulinic acid that the
  SI does not list (11.09 kcal/mol), and another atom order for the same impropers on malonic and oxalic acid and urea
  (energy differences 0.82, 0.03, 0.0005 kcal/mol).
- Acevedo's GROMACS files (GPL 3; nothing taken from them; 9 molecules, 102 atoms): charges and bends all equal, type
  names 100 / 102 (urea's HC), Lennard-Jones 94 / 102 (ethylene glycol, below); choline's every term equal.
  Where they differ from the SI, CAPS follows the SI: (1) no `[ pairs ]` section in any .itp, so the files have no 1-4
  Lennard-Jones or Coulomb although the paper scales them by 0.5; with CAPS's 1-4 list added, choline, levulinic acid
  and phenol match CAPS's total exactly (48.1794, 13.9829, 33.088 / 33.090 kcal/mol); (2) extra `[ exclusions ]` of
  1-5 O...H pairs in ethylene glycol, glycerol and oxalic acid (and H...O 1-4 pairs in malonic acid), which the paper
  does not describe; (3) ethylene glycol's epsilons are the SI's divided by 1.75 (CG 0.066 for 0.1155, HG 0.030 for
  0.0525, OG 0.170 for 0.2975; the files' own comments print the SI's values) and its C-O length is 1.43 A for the SI's
  1.410 (the comment says 1.410); (4) malonic acid's second CD-CT bond 1.552 A for 1.522 (the other copy and the
  comment 1.522); (5) oxalic acid's second HO-OH-CD-CD V2 2.871 kcal/mol (12.012 kJ/mol) for 5.5 (the other copy 23.012
  kJ/mol and the comment 5.50); (6) O-H force constants 463750 kJ/mol/nm^2 (554.2) beside 462750 (553) and phenol's C-C
  468.88 for 469 (unit rounding); (7) the same impropers with another atom order (urea 0.42 kcal/mol on the test geometry),
  and zero impropers on glycerol's CM and malonic acid's CT.
- Not in the paper and not added: impropers on phenol's ring and levulinic acid's carbonyls; Lennard-Jones or charges for
  any molecule outside the paper's set (the GPL repository's later additions, acetamide, ethylammonium and others, are
  not part of the 2018 paper).

## 4. Papers needed

| needed for | citation | DOI | state |
|---|---|---|---|
| O'Malley & Catlow's all-atom alkane torsion | N. Raj, G. Sastre, C. R. A. Catlow, J. Phys. Chem. B 103, 11007-11015 (1999) | 10.1021/jp991314f | read (`papers/jp991314f.pdf`): it has no torsion term (Table 3 lists bonds and bends only), so the torsion O'Malley & Catlow cite is not there. Still unresolved; omalley2013-alkanes keeps its zero torsion |
| methanol / faujasite as rigid ions (DL_FIELD's second reference for methanol2 and faujasite2) | A. J. O'Malley et al., Phys. Chem. Chem. Phys. 18, 17294 (2016) | to confirm | not local |
| PbTiO3's Pb-Ti term (Sepliarsky 2005 Table 1 row misprinted, section 3 notes) | M. Sepliarsky, R. E. Cohen, AIP Conf. Proc. 626, 36 (2002); M. Sepliarsky, Z. Wu, A. Asthagiri, R. E. Cohen, Ferroelectrics 301, 55 (2004) | 10.1063/1.1499550; 10.1080/00150190490454882 | both read (`papers/36_1_online.pdf`, `papers/Atomistic Model Potential for PbTiO3 and PMN ...pdf`): neither prints a parameter table. 2002 is an earlier set (Buckingham A-O, B-O, O-O, 20 parameters; 5-atom tetragonal cell -118.81 eV). 2004 is the review's set: Rydberg A-B, A-O, B-O and Buckingham O-O, 23 parameters, spring c2 w^2/2 + c4 w^4/24 (confirming the review's spring), but it gives no Rydberg expression. Its Fig. 1 cubic E(V) (minimum near 58.8 A^3, about +0.03 eV at 57.0 and +0.12 eV at 62.4 A^3) is matched by no reading of the Pb-Ti row: as printed the crystal is unbound; with Pb-Ti dropped, rho read as 1/2.4201 or the row as Buckingham, the minimum sits near 57.5 A^3 and E(62.4) - E(58.8) is +0.34 eV. PbTiO3 stays withheld; the parameter file the authors used (DL_POLY FIELD) or a later paper printing the set is needed. K(Nb,Ta)O3 and (Ba,Sr)TiO3 are done from the review (doi 10.1016/j.cossms.2006.05.002) |
| Bi2O3 (bismuth(iii)_oxide1) | M. Walther, D. Zahn, Eur. J. Inorg. Chem. 2015, 1178 (their ref. 11) | 10.1002/ejic.201403141 (to confirm) | Walther & Zahn, Chem. Phys. Lett. 691, 87 (2018) (`papers/1-s2.0-S0009261417310060-main.pdf`, doi 10.1016/j.cplett.2017.10.064) is read: it prints no parameters, only "the same set of force field parameters was employed [11]". Not built |
| the quartic core-shell expression of Catlow, Freeman & Royle's rutile (to confirm the K w^4/24 reading) | C. R. A. Catlow, R. James, Proc. R. Soc. Lond. A 384, 157 (1982); R. James, PhD thesis, University of London (1979) | 10.1098/rspa.1982.0153 | Catlow & James 1982 read (`papers/Catlow-DisorderTiO2x-1982.pdf`): its rutile model (Table 1: A12 656.74 eV, rho12 0.40431 A, C22 27.063, Y1 -35.863, Y2 -2.38856, K1 65974.0, K2 18.41) is the earlier, harmonic-only shell model that the 1985 paper replaced; it has no quartic term, so it neither confirms nor contradicts K w^4/24. catlow1985-rutile is unchanged; James's thesis is still the place to look |

The other sets of the former section 4 are built from their papers (now local in `~/project/ff-sources/papers/`):
Catlow, Faux & Norgett 1976 (doi 10.1088/0022-3719/9/3/008), Catlow, Freeman & Royle 1985 (10.1016/0378-4363(85)90134-2),
Bandura & Kubicki 2003 (10.1021/jp034093t), Khan, Islam & Bates 1998 (10.1039/a803917h), Sayle & Watson 2000
(10.1039/b006587k).

## 5. What is left

- PbTiO3 (Sepliarsky): the 2002 and 2004 papers print no parameters; the Pb-Ti row needs the authors' parameter file or a later paper that prints the set.
- The rutile quartic convention (K w^4/24): Catlow & James 1982 has no quartic term; R. James's 1979 thesis.
- Bi2O3 (Walther & Zahn): their 2015 paper.
- The O'Malley alkane torsion: Raj, Sastre & Catlow 1999 does not hold it; the source O'Malley & Catlow meant is unknown.
- Families handled by other sessions: OPLS (2005 / 2020 / 2024), PCFF, COMPASS, CVFF, ClayFF, CL&P.
- Acevedo's OPLS-DES is rebuilt from the paper as oplsdes-doherty2018 (section 3); its GPL 3 files stay a reference.
  Open readings, all recorded in the generator: the order of the outer atoms in the impropers (the SI writes only
  X-c-X-X) and the improper on choline's four-connected N (four triples, as the X-NA-X-X row reads). Ethylene glycol's
  hydrogen torsions are settled: zero, as in Szefczyk & Cordeiro, J. Phys. Chem. B 115, 3013 (2011), Table 2 and
  p. 3016 (read, `papers/jp109914s.pdf`), whose model has no torsion but O-C-C-O and C-C-O-H.
- `bench/typing/sweep_dlfield.sh` still sweeps the converted files where they exist on disk; the rebuilt libraries are
  checked with `validate_typing.py --map` as above.
