# CAPS force-field library

Built by `bench/ff/build_library.py` from the published parameter files. Status: **validated** — checked interaction by interaction and in energies / forces against their reference implementations and LAMMPS (evidence in each file's `validation` field); **converted** — parsed into CAPS rules, not yet compared; **pending** — the family's functional forms are not interpreted yet.

| id | force field | version / year | from | status | automatic typing |
|---|---|---|---|---|---|
| pcff | PCFF (diagonal terms) | cff91 family, 1994 | converted: PCFF.par | validated | yes |
| compass | COMPASS (diagonal terms) | 1998 | converted: COMPASS.par | validated | |
| cvff | CVFF | 1988 | converted: CVFF.par | validated | yes |
| opls2005 | OPLS-AA / OPLS 2005 | 2005 | converted: OPLS2005.par | validated | yes |
| opls2020 | OPLS 2020 bond / angle supplement | as distributed | converted: OPLS2020.par | converted | |
| opls-clp | CL&P ionic liquids (OPLS-AA based) | 2004 onwards | converted: OPLS_CL_P.par | validated | |
| opls-des | OPLS-DES deep eutectic solvents | 2018 | converted: OPLS_DES.par | converted | |
| gaff-amber16 | GAFF (AmberTools 16) | GAFF 1.8x, 2016 | converted: AMBER16_gaff.par | validated | yes |
| gaff-amber25 | GAFF (AmberTools 25) | 2025 | converted: AMBER25_gaff.par | validated | yes |
| charmm | CHARMM | | converted: CHARMM.par | validated | |
| cgenff | CGenFF (CHARMM36) | | converted: CHARMM36_cgenff.par | validated | yes |
| dreiding | DREIDING | 1990 | converted: DREIDING.par | validated | |
| trappe-ua | TraPPE-UA | 1998 | converted: TRAPPE_UA.par; CO2 added from Potoff & Siepmann 2001 (checked against gmso's TraPPE CO2 file, bench/ff/check_trappe_co2.py) | validated | |
| trappe-eh | TraPPE-EH | 2007 | converted: TRAPPE_EH.par | validated | |
| misc | Miscellaneous | | converted: MISC_FF.par | converted | |
| inorganic-binary-halides | Inorganic: binary halides | | converted: INORGANIC_binary_halides.par | validated | |
| inorganic-binary-misc | Inorganic: binary misc | | converted: INORGANIC_binary_misc.par | validated | |
| inorganic-binary-oxides | Inorganic: binary oxides | | converted: INORGANIC_binary_oxides.par | validated | |
| clayff-cygan2004 | ClayFF | Cygan, Liang & Kalinichev 2004 | the paper's Tables 1 and 2 (bench/ff/make_clayff_cygan2004.py) | validated | yes |
| inorganic-glass | Inorganic: glass | | converted: INORGANIC_glass.par | converted | |
| inorganic-ternary-oxides | Inorganic: ternary oxides | | converted: INORGANIC_ternary_oxides.par | validated | |
| inorganic-zeolite | Inorganic: zeolite | | converted: INORGANIC_zeolite.par | validated | |
| inorganic-zeolite-hill-sauer | Inorganic: zeolite Hill Sauer | | converted: INORGANIC_zeolite_Hill_Sauer.par | converted | |
