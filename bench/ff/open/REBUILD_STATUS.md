# Force-field libraries: provenance and rebuild status

Commit 860a126b stopped tracking 32 libraries (and their typing files) that had been converted from DL_FIELD 4.13's
`.par` / `.sf` files, whose licence forbids redistributing derived data. Each family is being rebuilt from its original
published source (paper, official parameter file or the authors' database). DL_FIELD stays on the developer's machine
as an accuracy reference for typing and energies only; no number in a rebuilt library comes from it.

State of 2026-10-13 (branch `ff-rebuild-rest`, commits 1053f7a2 and b2ee52e2; the inorganic sets of section 4 on branch
`inorg-papers`). "Local" means the source is in `~/project/ff-sources/`.

## 1. The 32 libraries removed in 860a126b

| removed library | was | rebuilt as | original source | local | licence of the source | status |
|---|---|---|---|---|---|---|
| amber | AMBER.par | amber-ff03-ambertools, amber-ff14sb-ambertools | AmberTools 23.6 `dat/leap/parm` (parm99, parm10, frcmod.ff03 / ff14SB / OL15) | yes (`ambertools-23.6`) | public domain | done (other session) |
| gaff-amber16 | AMBER16_gaff.par | gaff-1.81-ambertools | AmberTools `gaff.dat` | yes | public domain | done (other session) |
| gaff-amber25 | AMBER25_gaff.par | gaff2-2.2.20-ambertools (and gaff2-2.11) | AmberTools `gaff2.dat`, `gaff211.dat` | yes | public domain | done (other session) |
| charmm | CHARMM.par | charmm36-*-toppar | MacKerell lab `toppar_c36_jul17` | yes (`charmm-toppar-jul17`) | no licence stated: built locally, never committed | done (other session) |
| charmm19 | CHARMM19.par | charmm19-toppar | `toppar_c36_jul17` (param19 / toph19) | yes | as above | done (other session) |
| charmm22-prot | CHARMM22_prot.par | charmm22-prot-toppar | `toppar_c36_jul17` | yes | as above | done (other session) |
| charmm36-carb / -lipid / -nucl / -prot | CHARMM36_*.par | charmm36-carb/-lipid/-nucl/-prot-toppar | `toppar_c36_jul17` | yes | as above | done (other session) |
| cgenff | CHARMM36_cgenff.par | cgenff-toppar | `toppar_c36_jul17` (CGenFF 3.0.1 / 4.0) | yes | as above | done (other session) |
| pcff | PCFF.par | pcff-frc | `pcff.frc` | - | - | other session |
| compass | COMPASS.par | compass-frc | `compass_published.frc`; H. Sun, J. Phys. Chem. B 102, 7338 (1998) | - | - | other session |
| cvff | CVFF.par | cvff-frc | `cvff.frc` | - | - | other session |
| dreiding | DREIDING.par | dreiding1990 | Mayo, Olafson & Goddard, J. Phys. Chem. 94, 8897 (1990) | yes (`papers/j100389a010.pdf`) | published parameters | done (other session) |
| gromos-54a7 | GROMOS_G54A7.par | gromos54a7-gromacs | GROMACS `share/top/gromos54a7.ff` | yes (GROMACS) | LGPL 2.1 | done (other session) |
| opls2005 | OPLS2005.par | oplsaa2024 (and oplsaa-gromacs) | Jorgensen's BOSS files, J. Phys. Chem. B 128, 250 (2024) SI; GROMACS `oplsaa.ff` | yes | SI tables; LGPL 2.1 | other session (OPLS) |
| opls2020 | OPLS2020.par | oplsaa2024 | as above | yes | as above | other session (OPLS) |
| opls-aam | OPLS_AAM.par | oplsaam-jorgensen | Jorgensen lab `top_opls_aam.inp`, `par_opls_aam.inp` | yes (`opls-aam`) | CC BY 4.0 | done (other session) |
| opls-clp | OPLS_CL_P.par | clandp | Padua's `il.ff` (github.com/paduagroup/clandp) | yes (`clandp-dl` is the reference only) | MIT | other session (CL&P) |
| opls-des | OPLS_DES.par | oplsdes-doherty2018 (the set DL_FIELD's file held); clandp-des (CL&P's DES set) | oplsdes-doherty2018: Doherty & Acevedo, J. Phys. Chem. B 122, 9982 (2018) and its SI (Tables S1-S6, Schemes S1-S2), transcribed; Acevedo's own GROMACS files (github.com/orlandoacevedo/DES, GPL 3, local in `opls-des-dl/DES-master`) are an accuracy reference only, nothing taken from them. clandp-des: Padua's `des.ff` of 2025/12/19 (github.com/paduagroup/clandp) | yes (`papers/jp8b06647.pdf`, `papers/jp8b06647_si_001.pdf`; `clandp-dl/clandp-master/des.ff`) | published parameters; MIT | done (section 3) |
| trappe-ua | TRAPPE_UA.par | trappe-ua-umn, trappe-ua-foyer | TraPPE database (trappe.oit.umn.edu); foyer `trappe-ua.xml` | yes (`trappe-site`, `foyer-trappe`) | published parameters; MIT (foyer) | done (other session) |
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
| gaff-1.81-ambertools, gaff2-2.11-ambertools, gaff2-2.2.20-ambertools, amber-ff03-ambertools, amber-ff14sb-ambertools | AmberTools 23.6 `dat/leap/parm` | public domain | original source |
| oplsaa-gromacs, gromos54a7-gromacs | GROMACS `oplsaa.ff`, `gromos54a7.ff` | LGPL 2.1 | original source |
| oplsaa2024, oplsua2024 | Jorgensen's BOSS files (J. Phys. Chem. B 128, 250 (2024) SI) | published SI tables | original source (OPLS session) |
| oplsaam-jorgensen | Jorgensen lab CHARMM-format files | CC BY 4.0 | original source |
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
| martini22-proteins | vermouth-martinize `force_fields/martini22` | Apache 2.0 | original distribution |
| martini3 | `martini_v3.0.0.itp` (cgmartini.nl) | as distributed by the Martini group | original distribution |
| glycam*-ambertools, *-toppar | AmberTools GLYCAM (own copyright), CHARMM toppar (no licence) | - | built locally, never committed (`.gitignore`) |

## 3. Families rebuilt on this branch

| family | library | generator | typing against DL_FIELD | LAMMPS | GROMACS | ctest |
|---|---|---|---|---|---|---|
| MISC: alkanes in silicalite | omalley2013-alkanes | `bench/ff/make_zeolite_sorbates.py` | DL_FIELD's `alkanes1` is a hand-labelled template (Cp1 / Cs1 / Ct1 / H1 by the number of H) that DL_FIELD cannot type automatically; CAPS types CH3 / CH2 / CH / C / CH4 and H by the same rule | agrees (alpha-quartz; alkanes with SiO4) | agrees (alkanes) | `FFDef.ZeoliteSorbatesFromTheirPapers` |
| MISC: methanol in NaY | plant2007-methanol | same | `methanol2` 6 / 6 atoms | agrees (methanol, Na+, AlO4, SiO4) | agrees | same |
| MISC: chloroform in NaY | ramsahye2005-nay (existing) | `make_ramsahye2005.py` | `chloroform3` 5 / 5 atoms | agrees (existing case) | - | existing |
| MISC: HFA134a | peguin2009-hfa134a | `make_zeolite_sorbates.py` | `HFA_134a` 8 / 8 atoms | agrees | agrees | same |
| TraPPE-EH | trappe-eh-umn | `bench/ff/make_trappe_eh.py SITE_DIR --check` | 9 templates, 94 / 94 atoms (17 carry the key of a symmetry-equivalent position); the database's own 65 models 848 / 848 sites | agrees (flexible and rigid sets) | agrees | `FFDef.TrappeEhModelsTypedAsThemselves` |
| CL&P deep-eutectic solvents | clandp-des | `bench/ff/convert_fftool.py des.ff` (47 types, 16 bonds, 30 angles, 44 dihedrals, 6 impropers), `bench/typing/make_clandp_des_rules.py` (45 rules) | not comparable: DL_FIELD's OPLS_DES is Acevedo's set. Against Padua's own files (`validate_fftool_types.py`): the 14 whose types des.ff defines (cholinium.xyz, N1110, N1111, N2220, N2222, N4444, P1111, P66614, c1c1pyrr, c3c1pyrr, c4c1pyrr, Cl, Br, oac) 362 / 362 atoms, each at its integer charge; 29 DES components typed whole and neutral or at their ion charge, 7 that des.ff cannot make neutral left untyped (`check_clandp_des_components.py`) | agrees: reline and glyceline boxes (608 and 600 atoms, periodic) and a vacuum set of 10 components, energy terms 1.0e-7 to 3.8e-7, forces <= 1.9e-5 kcal/mol/A, virial <= 2.2e-6 (exact styles); native styles bonded and vdW 1.6e-8 to 3.8e-7 | agrees: reline and glyceline with PME 6.9e-6 and 3.3e-5, vacuum set 5.0e-6 | `FFDef.ClandpDesComponentsTypedWhole`, python_bindings |
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

**OPLS-DES (oplsdes-doherty2018).** Built from B. Doherty, O. Acevedo, J. Phys. Chem. B 122, 9982 (2018) and its SI
only. The paper: OPLS-AA form, geometric mixing, 1-4 pairs x 0.5, charges built for +-0.8 e ions (the SI prints them
scaled), choline chloride's charges and Lennard-Jones refitted to reline's ND/EPSR structure, choline's bonded terms from
OPLS-2009IL, the HBDs' bonded terms OPLS-AA's (ethylene glycol's torsions Szefczyk & Cordeiro's), the HBDs' non-bonded
terms tuned for DES mixtures only; GROMACS with a 16 A cut-off. The SI prints every term (Tables S1-S6), so no term
comes from another library. Readings: ethylene glycol's rows in OPLS-AA class names (CT, HC) and its O called OH are
Table S4's CG, HG, OG; levulinic acid's bend "HO-CB-CZ" is OH-CB-CZ; the X-c-X-X rows are impropers on c (centre second,
every triple of neighbours); the all-zero X-CM-X-X and X-CT-X-X are left out; ethylene glycol's block lists only
OG-CG-CG-OG and CG-CG-OG-HO, so its H-C-C-H, H-C-C-O and H-C-O-H torsions are zero (both references have none either);
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
  Open readings, all recorded in the generator: ethylene glycol's hydrogen torsions (zero: the SI lists none; OPLS-AA
  would give 0.300 / 0.468 / 0.352 kcal/mol), the order of the outer atoms in the impropers (the SI writes only
  X-c-X-X), and the improper on choline's four-connected N (four triples, as the X-NA-X-X row reads). Szefczyk &
  Cordeiro, J. Phys. Chem. B 115, 3013 (2011) (ethylene glycol's torsions) would settle the first; not local.
- `bench/typing/sweep_dlfield.sh` still sweeps the converted files where they exist on disk; the rebuilt libraries are
  checked with `validate_typing.py --map` as above.
