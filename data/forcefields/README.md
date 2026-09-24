# CAPS force-field library

Built by `bench/ff/build_library.py` from moltemplate and DL_FIELD 4.13. Status: **validated** — checked interaction by interaction and in energies / forces against the source tool and LAMMPS (evidence in each file's `validation` field); **converted** — parsed into CAPS rules, not yet compared; **pending** — the family's functional forms are not interpreted yet; **template** / **alias** — not a separate force field.

| id | force field | version / year | from | status | automatic typing |
|---|---|---|---|---|---|
| pcff-dlfield | PCFF (DL_FIELD, diagonal terms) | cff91 family, 1994; DL_FIELD library v4.5 (2019) | dlfield: PCFF.par | validated | yes |
| compass-dlfield | COMPASS (DL_FIELD, diagonal terms) | 1998 | dlfield: COMPASS.par | validated |  |
| cvff-dlfield | CVFF (DL_FIELD) | 1988 | dlfield: CVFF.par | validated | yes |
| opls2005-dlfield | OPLS-AA / OPLS 2005 (DL_FIELD) | 2005 | dlfield: OPLS2005.par | validated | yes |
| opls2020-dlfield | OPLS 2020 bond / angle supplement (DL_FIELD) | as distributed with DL_FIELD 4.13 | dlfield: OPLS2020.par | converted |  |
| opls-aam-dlfield | OPLS-AA/M (DL_FIELD) | 2015 | dlfield: OPLS_AAM.par | converted |  |
| opls-clp-dlfield | CL&P ionic liquids (OPLS-AA based, DL_FIELD) | 2004 onwards | dlfield: OPLS_CL_P.par | validated |  |
| opls-des-dlfield | OPLS-DES deep eutectic solvents (DL_FIELD) | 2018 | dlfield: OPLS_DES.par | converted |  |
| amber-dlfield | AMBER (DL_FIELD) | Cornell et al. 1995 family | dlfield: AMBER.par | validated |  |
| gaff-amber16-dlfield | GAFF (AmberTools 16, DL_FIELD) | GAFF 1.8x, 2016 | dlfield: AMBER16_gaff.par | validated | yes |
| gaff-amber25-dlfield | GAFF (AmberTools 25, DL_FIELD) | 2025 | dlfield: AMBER25_gaff.par | validated | yes |
| charmm-dlfield | CHARMM (DL_FIELD) |  | dlfield: CHARMM.par | validated |  |
| charmm19-dlfield | CHARMM19 united atom (DL_FIELD) |  | dlfield: CHARMM19.par | validated |  |
| charmm22-prot-dlfield | CHARMM22 proteins (DL_FIELD) | 1998 | dlfield: CHARMM22_prot.par | validated |  |
| charmm36-carb-dlfield | CHARMM36 carbohydrates (DL_FIELD) |  | dlfield: CHARMM36_carb.par | converted |  |
| cgenff-dlfield | CGenFF (CHARMM36, DL_FIELD) |  | dlfield: CHARMM36_cgenff.par | validated | yes |
| charmm36-lipid-dlfield | CHARMM36 lipids (DL_FIELD) | 2010 | dlfield: CHARMM36_lipid.par | converted |  |
| charmm36-nucl-dlfield | CHARMM36 nucleic acids (DL_FIELD) |  | dlfield: CHARMM36_nucl.par | validated |  |
| charmm36-prot-dlfield | CHARMM36 proteins (DL_FIELD) | 2012 | dlfield: CHARMM36_prot.par | validated |  |
| dreiding-dlfield | DREIDING (DL_FIELD) | 1990 | dlfield: DREIDING.par | validated |  |
| gromos-54a7-dlfield | GROMOS 54A7 (DL_FIELD) | 2011 | dlfield: GROMOS_G54A7.par | validated |  |
| trappe-ua-dlfield | TraPPE-UA (DL_FIELD) | 1998 | dlfield: TRAPPE_UA.par | validated |  |
| trappe-eh-dlfield | TraPPE-EH (DL_FIELD) | 2007 | dlfield: TRAPPE_EH.par | validated |  |
| misc-dlfield | Miscellaneous (DL_FIELD) |  | dlfield: MISC_FF.par | converted |  |
| inorganic-binary-halides-dlfield | Inorganic: binary halides (DL_FIELD) |  | dlfield: INORGANIC_binary_halides.par | validated |  |
| inorganic-binary-misc-dlfield | Inorganic: binary misc (DL_FIELD) |  | dlfield: INORGANIC_binary_misc.par | validated |  |
| inorganic-binary-oxides-dlfield | Inorganic: binary oxides (DL_FIELD) |  | dlfield: INORGANIC_binary_oxides.par | validated |  |
| inorganic-clay-dlfield | Inorganic: clay (DL_FIELD) |  | dlfield: INORGANIC_clay.par | converted |  |
| inorganic-glass-dlfield | Inorganic: glass (DL_FIELD) |  | dlfield: INORGANIC_glass.par | converted |  |
| inorganic-ternary-oxides-dlfield | Inorganic: ternary oxides (DL_FIELD) |  | dlfield: INORGANIC_ternary_oxides.par | validated |  |
| inorganic-zeolite-dlfield | Inorganic: zeolite (DL_FIELD) |  | dlfield: INORGANIC_zeolite.par | validated |  |
| inorganic-zeolite-hill-sauer-dlfield | Inorganic: zeolite Hill Sauer (DL_FIELD) |  | dlfield: INORGANIC_zeolite_Hill_Sauer.par | converted |  |
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
