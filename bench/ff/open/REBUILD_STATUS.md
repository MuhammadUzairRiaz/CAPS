# Force-field libraries: provenance and rebuild status

Commit 860a126b stopped tracking 32 libraries (and their typing files) that had been converted from DL_FIELD 4.13's
`.par` / `.sf` files, whose licence forbids redistributing derived data. Each family is being rebuilt from its original
published source (paper, official parameter file or the authors' database). DL_FIELD stays on the developer's machine
as an accuracy reference for typing and energies only; no number in a rebuilt library comes from it.

State of 2026-10-13 (branch `ff-rebuild-rest`, commits 1053f7a2 and b2ee52e2). "Local" means the source is in `~/project/ff-sources/`.

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
| inorganic-binary-oxides | INORGANIC_binary_oxides.par | bush1994-oxides, lewis-catlow1985-oxides, catlow1977-uo2, matsui1991-tio2, sanders1984-sio2 | the papers of DL_FIELD's refs 1-6 | yes | published parameters | partly done: six sets wait on papers (section 4) |
| inorganic-clay | INORGANIC_clay.par | clayff2004 | Cygan, Liang & Kalinichev, J. Phys. Chem. B 108, 1255 (2004) | yes | published parameters | done (other session) |
| inorganic-glass | INORGANIC_glass.par | teter-glass, tilocca2006-sm, vogel2003-lipo3 | Cormack 2002, Du 2009, Okhotnikov 2013, Tilocca 2006, Vogel 2003 | yes | published parameters | done |
| inorganic-ternary-oxides | INORGANIC_ternary_oxides.par | xiao2011-caco3, gunn2012-pyrochlore, la-perovskites | Xiao 2011, Gunn 2012, Cherry 1995, Read 1999 | yes | published parameters | partly done: Sepliarsky's perovskites wait on the paper (section 4) |
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
| catlow1977-halides, catlow1977-uo2, lewis-catlow1985-oxides, matsui1991-tio2, sanders1984-sio2, gan2020, clayff2004, teter-glass, tilocca2006-sm, vogel2003-lipo3, xiao2011-caco3, gunn2012-pyrochlore, la-perovskites, schroder1992-zeolite, ramsahye2005-nay, forester1997-silicalite, hill-sauer1995 | each paper's tables | published parameters | original source |
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

Typing commands: `python3 bench/ff/validate_typing.py MISC_FF data/forcefields/ID.json data/typing/ID.typing.json --map ID --only TEMPLATE`
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

## 4. Papers needed (not local)

| needed for | citation | DOI |
|---|---|---|
| O'Malley & Catlow's all-atom alkane torsion | N. Raj, G. Sastre, C. R. A. Catlow, J. Phys. Chem. B 103, 11007-11015 (1999) | 10.1021/jp991434h (to confirm) |
| methanol / faujasite as rigid ions (DL_FIELD's second reference for methanol2 and faujasite2) | A. J. O'Malley et al., Phys. Chem. Chem. Phys. 18, 17294 (2016) | to confirm |
| BaTiO3, SrTiO3, PbTiO3, KTaO3, KNbO3 shell models | M. Sepliarsky, A. Asthagiri, S. R. Phillpot, M. G. Stachiotti, R. L. Migoni, Curr. Opin. Solid State Mater. Sci. 9, 107-113 (2005) | 10.1016/j.cossms.2006.01.001 (to confirm). `papers/1-s2.0-S1359028606000271-main.pdf` is another article of that issue (Kolasinski's editorial), not this one |
| MgO rigid-ion model (DL_FIELD's magnesium_oxide3) | C. R. A. Catlow, I. D. Faux, M. J. Norgett, J. Phys. C: Solid State Phys. 9, 419 (1976) | 10.1088/0022-3719/9/3/008 (to confirm) |
| rutile and alpha-quartz (rutile1, silicon_dioxide1) | C. R. A. Catlow, C. M. Freeman, R. L. Royle, Physica B+C 131, 1 (1985) | 10.1016/0378-4363(85)90134-3 (to confirm) |
| hydroxylated TiO2 surface (TiO2_hydroxylated1) | A. V. Bandura, J. D. Kubicki, J. Phys. Chem. B 107, 11072 (2003) | 10.1021/jp034093r (to confirm) |
| Bi2O3 (bismuth(iii)_oxide1) | M. Walther, D. Zahn, Chem. Phys. Lett. 691, 87 (2018); DL_FIELD's set came from a FIELD file the authors sent, so the paper may not hold every term | to confirm |
| ZrO2 (zirconium(iv)_oxide1) | M. S. Khan, M. S. Islam, D. R. Bates, J. Mater. Chem. 8, 2299 (1998) | to confirm |
| MgO, CaO, SrO, BaO rigid ions (magnesium_oxide4, calcium/strontium/barium_oxide2) | D. C. Sayle, G. W. Watson, Phys. Chem. Chem. Phys. 2, 5491 (2000) | to confirm |

## 5. What is left

- The binary-oxide and ternary-oxide sets of section 4, once their papers are here.
- The O'Malley alkane torsion (Raj, Sastre & Catlow 1999).
- Families handled by other sessions: OPLS (2005 / 2020 / 2024), PCFF, COMPASS, CVFF, ClayFF, CL&P / DES.
- `bench/typing/sweep_dlfield.sh` still sweeps the converted files where they exist on disk; the rebuilt libraries are
  checked with `validate_typing.py --map` as above.
