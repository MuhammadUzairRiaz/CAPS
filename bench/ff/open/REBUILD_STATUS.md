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
| opls-des | OPLS_DES.par | (clandp-des, in progress) | Doherty & Acevedo, J. Phys. Chem. B 122, 9982 (2018) | yes (`opls-des-dl` reference) | - | other session (DES) |
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
- Families handled by other sessions: OPLS (2005 / 2020 / 2024), PCFF, COMPASS, CVFF, ClayFF, CL&P / DES.
- `bench/typing/sweep_dlfield.sh` still sweeps the converted files where they exist on disk; the rebuilt libraries are
  checked with `validate_typing.py --map` as above.
