# CAPS force-field library

Built by `bench/ff/build_library.py` from moltemplate and the reference tool 4.13. Status: **validated** — checked interaction by interaction and in energies / forces against the source tool and LAMMPS (evidence in each file's `validation` field); **converted** — parsed into CAPS rules, not yet compared; **pending** — the family's functional forms are not interpreted yet; **template** / **alias** — not a separate force field.

| id | force field | version / year | from | status | automatic typing |
|---|---|---|---|---|---|
| pcff-reftool | PCFF (the reference tool, diagonal terms) | cff91 family, 1994; the reference tool library v4.5 (2019) | reftool: PCFF.par | validated | yes |
| compass-reftool | COMPASS (the reference tool, diagonal terms) | 1998 | reftool: COMPASS.par | validated |  |
| cvff-reftool | CVFF (the reference tool) | 1988 | reftool: CVFF.par | validated | yes |
| opls2005-reftool | OPLS-AA / OPLS 2005 (the reference tool) | 2005 | reftool: OPLS2005.par | validated | yes |
| opls2020-reftool | OPLS 2020 bond / angle supplement (the reference tool) | as distributed with the reference tool 4.13 | reftool: OPLS2020.par | converted |  |
| opls-aam-reftool | OPLS-AA/M (the reference tool) | 2015 | reftool: OPLS_AAM.par | converted |  |
| opls-clp-reftool | CL&P ionic liquids (OPLS-AA based, the reference tool) | 2004 onwards | reftool: OPLS_CL_P.par | validated |  |
| opls-des-reftool | OPLS-DES deep eutectic solvents (the reference tool) | 2018 | reftool: OPLS_DES.par | converted |  |
| amber-reftool | AMBER (the reference tool) | Cornell et al. 1995 family | reftool: AMBER.par | validated |  |
| gaff-amber16-reftool | GAFF (AmberTools 16, the reference tool) | GAFF 1.8x, 2016 | reftool: AMBER16_gaff.par | validated | yes |
| gaff-amber25-reftool | GAFF (AmberTools 25, the reference tool) | 2025 | reftool: AMBER25_gaff.par | validated | yes |
| charmm-reftool | CHARMM (the reference tool) |  | reftool: CHARMM.par | validated |  |
| charmm19-reftool | CHARMM19 united atom (the reference tool) |  | reftool: CHARMM19.par | validated |  |
| charmm22-prot-reftool | CHARMM22 proteins (the reference tool) | 1998 | reftool: CHARMM22_prot.par | validated |  |
| charmm36-carb-reftool | CHARMM36 carbohydrates (the reference tool) |  | reftool: CHARMM36_carb.par | converted |  |
| cgenff-reftool | CGenFF (CHARMM36, the reference tool) |  | reftool: CHARMM36_cgenff.par | validated | yes |
| charmm36-lipid-reftool | CHARMM36 lipids (the reference tool) | 2010 | reftool: CHARMM36_lipid.par | converted |  |
| charmm36-nucl-reftool | CHARMM36 nucleic acids (the reference tool) |  | reftool: CHARMM36_nucl.par | validated |  |
| charmm36-prot-reftool | CHARMM36 proteins (the reference tool) | 2012 | reftool: CHARMM36_prot.par | validated |  |
| dreiding-reftool | DREIDING (the reference tool) | 1990 | reftool: DREIDING.par | validated |  |
| gromos-54a7-reftool | GROMOS 54A7 (the reference tool) | 2011 | reftool: GROMOS_G54A7.par | validated |  |
| trappe-ua-reftool | TraPPE-UA (the reference tool) | 1998 | reftool: TRAPPE_UA.par | validated |  |
| trappe-eh-reftool | TraPPE-EH (the reference tool) | 2007 | reftool: TRAPPE_EH.par | validated |  |
| misc-reftool | Miscellaneous (the reference tool) |  | reftool: MISC_FF.par | converted |  |
| inorganic-binary-halides-reftool | Inorganic: binary halides (the reference tool) |  | reftool: INORGANIC_binary_halides.par | validated |  |
| inorganic-binary-misc-reftool | Inorganic: binary misc (the reference tool) |  | reftool: INORGANIC_binary_misc.par | validated |  |
| inorganic-binary-oxides-reftool | Inorganic: binary oxides (the reference tool) |  | reftool: INORGANIC_binary_oxides.par | validated |  |
| inorganic-clay-reftool | Inorganic: clay (the reference tool) |  | reftool: INORGANIC_clay.par | converted |  |
| inorganic-glass-reftool | Inorganic: glass (the reference tool) |  | reftool: INORGANIC_glass.par | converted |  |
| inorganic-ternary-oxides-reftool | Inorganic: ternary oxides (the reference tool) |  | reftool: INORGANIC_ternary_oxides.par | validated |  |
| inorganic-zeolite-reftool | Inorganic: zeolite (the reference tool) |  | reftool: INORGANIC_zeolite.par | validated |  |
| inorganic-zeolite-hill-sauer-reftool | Inorganic: zeolite Hill Sauer (the reference tool) |  | reftool: INORGANIC_zeolite_Hill_Sauer.par | converted |  |
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
