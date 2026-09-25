# CAPS force-field library

Built by `bench/ff/build_library.py` from the published parameter files. Status: **validated** — checked interaction by interaction and in energies / forces against their reference implementations and LAMMPS (evidence in each file's `validation` field); **converted** — parsed into CAPS rules, not yet compared; **pending** — the family's functional forms are not interpreted yet; **template** / **alias** — not a separate force field.

| id | force field | version / year | from | status | automatic typing |
|---|---|---|---|---|---|
| pcff | PCFF (diagonal terms) | cff91 family, 1994 | converted: PCFF.par | validated | yes |
| compass | COMPASS (diagonal terms) | 1998 | converted: COMPASS.par | validated | |
| cvff | CVFF | 1988 | converted: CVFF.par | validated | yes |
| opls2005 | OPLS-AA / OPLS 2005 | 2005 | converted: OPLS2005.par | validated | yes |
| opls2020 | OPLS 2020 bond / angle supplement | as distributed | converted: OPLS2020.par | converted | |
| opls-aam | OPLS-AA/M | 2015 | converted: OPLS_AAM.par | converted | |
| opls-clp | CL&P ionic liquids (OPLS-AA based) | 2004 onwards | converted: OPLS_CL_P.par | validated | |
| opls-des | OPLS-DES deep eutectic solvents | 2018 | converted: OPLS_DES.par | converted | |
| amber | AMBER | Cornell et al. 1995 family | converted: AMBER.par | validated | |
| gaff-amber16 | GAFF (AmberTools 16) | GAFF 1.8x, 2016 | converted: AMBER16_gaff.par | validated | yes |
| gaff-amber25 | GAFF (AmberTools 25) | 2025 | converted: AMBER25_gaff.par | validated | yes |
| charmm | CHARMM | | converted: CHARMM.par | validated | |
| charmm19 | CHARMM19 united atom | | converted: CHARMM19.par | validated | |
| charmm22-prot | CHARMM22 proteins | 1998 | converted: CHARMM22_prot.par | validated | |
| charmm36-carb | CHARMM36 carbohydrates | | converted: CHARMM36_carb.par | converted | |
| cgenff | CGenFF (CHARMM36) | | converted: CHARMM36_cgenff.par | validated | yes |
| charmm36-lipid | CHARMM36 lipids | 2010 | converted: CHARMM36_lipid.par | converted | |
| charmm36-nucl | CHARMM36 nucleic acids | | converted: CHARMM36_nucl.par | validated | |
| charmm36-prot | CHARMM36 proteins | 2012 | converted: CHARMM36_prot.par | validated | |
| dreiding | DREIDING | 1990 | converted: DREIDING.par | validated | |
| gromos-54a7 | GROMOS 54A7 | 2011 | converted: GROMOS_G54A7.par | validated | |
| trappe-ua | TraPPE-UA | 1998 | converted: TRAPPE_UA.par | validated | |
| trappe-eh | TraPPE-EH | 2007 | converted: TRAPPE_EH.par | validated | |
| misc | Miscellaneous | | converted: MISC_FF.par | converted | |
| inorganic-binary-halides | Inorganic: binary halides | | converted: INORGANIC_binary_halides.par | validated | |
| inorganic-binary-misc | Inorganic: binary misc | | converted: INORGANIC_binary_misc.par | validated | |
| inorganic-binary-oxides | Inorganic: binary oxides | | converted: INORGANIC_binary_oxides.par | validated | |
| inorganic-clay | Inorganic: clay | | converted: INORGANIC_clay.par | converted | |
| inorganic-glass | Inorganic: glass | | converted: INORGANIC_glass.par | converted | |
| inorganic-ternary-oxides | Inorganic: ternary oxides | | converted: INORGANIC_ternary_oxides.par | validated | |
| inorganic-zeolite | Inorganic: zeolite | | converted: INORGANIC_zeolite.par | validated | |
| inorganic-zeolite-hill-sauer | Inorganic: zeolite Hill Sauer | | converted: INORGANIC_zeolite_Hill_Sauer.par | converted | |
| compass-published-moltemplate | COMPASS (published subset, full class II) | 1998 | moltemplate: compass_published.lt | validated |  |
| gaff-moltemplate | GAFF (moltemplate) | GAFF 1.x | moltemplate: gaff.lt | converted |  |
| gaff2-moltemplate | GAFF2 (moltemplate) | GAFF 2 | moltemplate: gaff2.lt | validated |  |
| oplsaa2024-moltemplate | OPLS-AA (2024 parameter file) | 2024 | moltemplate: oplsaa2024.lt | validated |  |
| oplsaa2008-moltemplate | OPLS-AA (BOSS 4.8, 2008) | 2008 | moltemplate: oplsaa2008.lt | converted |  |
| loplsaa2024-moltemplate | L-OPLS overlay (on OPLS-AA 2024) | 2012 / 2015 | moltemplate: loplsaa2024.lt | converted |  |
| loplsaa2008-moltemplate | L-OPLS overlay (on OPLS-AA 2008) | 2012 / 2015 | moltemplate: loplsaa2008.lt | converted |  |
| dreiding-moltemplate | DREIDING (moltemplate) | 1990 | moltemplate: dreiding.lt | converted |  |
| trappe1998-moltemplate | TraPPE-UA alkanes | 1998 | moltemplate: trappe1998.lt | converted |  |
| sdk-moltemplate | SDK coarse-grained | 2007 / 2010 | moltemplate: sdk.lt | converted |  |
| martini-moltemplate | MARTINI 2.0 | 2007 | moltemplate: martini.lt | converted |  |
| drymartini-moltemplate | Dry MARTINI | 2015 | moltemplate: drymartini.lt | converted |  |
| cooke-deserno-moltemplate | Cooke–Deserno lipid model | 2005 | moltemplate: cooke_deserno_lipid.lt | converted |  |
| graphene-moltemplate | Graphene (LJ carbon) |  | moltemplate: graphene.lt | converted |  |
| spce-moltemplate | SPC/E water | 1987 | moltemplate: spce.lt | converted |  |
| tip3p-1983-moltemplate | TIP3P water (1983) | 1983 | moltemplate: tip3p_1983.lt | converted |  |
| tip3p-2004-moltemplate | TIP3P water (Ewald, 2004) | 2004 | moltemplate: tip3p_2004.lt | converted |  |
| mw-moltemplate | mW water | 2009 | moltemplate: watmw.lt | converted |  |
| graphite-moltemplate | graphite |  | moltemplate: graphite.lt | template |  |
| spc_oplsaa-moltemplate | spc_oplsaa |  | moltemplate: spc_oplsaa.lt | template |  |
| spc_oplsaa2008-moltemplate | spc_oplsaa2008 |  | moltemplate: spc_oplsaa2008.lt | template |  |
| spc_oplsaa2024-moltemplate | spc_oplsaa2024 |  | moltemplate: spc_oplsaa2024.lt | template |  |
| spce_oplsaa-moltemplate | spce_oplsaa |  | moltemplate: spce_oplsaa.lt | template |  |
| spce_oplsaa2024-moltemplate | spce_oplsaa2024 |  | moltemplate: spce_oplsaa2024.lt | template |  |
| spce_ice_rect8-moltemplate | spce_ice_rect8 |  | moltemplate: spce_ice_rect8.lt | template |  |
| spce_ice_rect16-moltemplate | spce_ice_rect16 |  | moltemplate: spce_ice_rect16.lt | template |  |
| spce_ice_rect32-moltemplate | spce_ice_rect32 |  | moltemplate: spce_ice_rect32.lt | template |  |
| tip3p_1983_oplsaa-moltemplate | tip3p_1983_oplsaa |  | moltemplate: tip3p_1983_oplsaa.lt | template |  |
| tip3p_1983_oplsaa2008-moltemplate | tip3p_1983_oplsaa2008 |  | moltemplate: tip3p_1983_oplsaa2008.lt | template |  |
| tip3p_1983_oplsaa2024-moltemplate | tip3p_1983_oplsaa2024 |  | moltemplate: tip3p_1983_oplsaa2024.lt | template |  |
| tip3p_2004_oplsaa-moltemplate | tip3p_2004_oplsaa |  | moltemplate: tip3p_2004_oplsaa.lt | template |  |
| tip3p_2004_oplsaa2024-moltemplate | tip3p_2004_oplsaa2024 |  | moltemplate: tip3p_2004_oplsaa2024.lt | template |  |
| tip5p_oplsaa-moltemplate | tip5p_oplsaa |  | moltemplate: tip5p_oplsaa.lt | template |  |
| tip5p_oplsaa2008-moltemplate | tip5p_oplsaa2008 |  | moltemplate: tip5p_oplsaa2008.lt | template |  |
| tip5p_oplsaa2024-moltemplate | tip5p_oplsaa2024 |  | moltemplate: tip5p_oplsaa2024.lt | template |  |
| oplsaa-moltemplate | oplsaa |  | moltemplate: oplsaa.lt | alias |  |
| loplsaa-moltemplate | loplsaa |  | moltemplate: loplsaa.lt | alias |  |
| graphene_rectangular-moltemplate | graphene_rectangular |  | moltemplate: graphene_rectangular.lt | alias |  |
| spce_amber-moltemplate | spce_amber |  | moltemplate: spce_amber.lt | alias |  |
| spce_dreiding-moltemplate | spce_dreiding |  | moltemplate: spce_dreiding.lt | alias |  |
| spce_hybrid-moltemplate | spce_hybrid |  | moltemplate: spce_hybrid.lt | alias |  |
| spce_more_comments-moltemplate | spce_more_comments |  | moltemplate: spce_more_comments.lt | alias |  |
| tip3p_1983_charmm-moltemplate | tip3p_1983_charmm |  | moltemplate: tip3p_1983_charmm.lt | alias |  |
| tip3p_1983_charmm_hybrid-moltemplate | tip3p_1983_charmm_hybrid |  | moltemplate: tip3p_1983_charmm_hybrid.lt | alias |  |
| tip3p_1983_hybrid-moltemplate | tip3p_1983_hybrid |  | moltemplate: tip3p_1983_hybrid.lt | alias |  |
| tip3p_2004_hybrid-moltemplate | tip3p_2004_hybrid |  | moltemplate: tip3p_2004_hybrid.lt | alias |  |
