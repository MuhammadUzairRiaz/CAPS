# Torsion census: oplsaa2024

Structures typed: 211 (21 copolymer, 114 polymer, 38 product, 38 reactant); failed: 0. Environment radius: 2 bonds.

| | published | CAPS fit | missing | total |
|---|---|---|---|---|
| torsion classes | 250 | 77 | 17 | 344 |
| rotor environments | 1106 | 155 | 29 | 1290 |
| of them heavy-atom rotors | 858 | 138 | 28 | 1024 |
| bond classes | - | 0 | 0 | - |
| angle classes | - | 6 | 4 | - |

## Campaign size

- One relaxed scan per rotor environment, 24 points each (15 degrees).
- Library (environments the published table does not cover): 184 scans, 4416 QM points.
- Training set (every environment, published ones included): 1290 scans, 30960 QM points before distorted conformers are added.
- Fragment size: median 7 heavy atoms, largest 17.
- Environments seen once: 755 (many are chain ends or one crosslink site; the ones seen often are the backbones).

## Missing torsion classes

- C-CT-OH-HO (1x; CC(O)C(=O)O)
- C-OS-C-CB (2x; O=C1OC(=O)c2ccccc21)
- CA-CT-CT-CM (26x; C002 Styrene–butadiene rubber (SBR, 23.5 wt % styrene), C108 E-SBR (23.5 wt % styrene; butadiene 70 % trans, 15 % cis, 15 % vinyl), CC/C=C\CC(c1ccccc1)C/C=C\C)
- CA-NT-CT-CT (21x; epoxy_amine_primary: DGEBA + DDM (aromatic amine), epoxy_amine_primary: DGEBA + DDS (aromatic sulfone amine), epoxy_amine_primary: TGDDM + DDS)
- CA-NT-CT-CY (7x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- CA-OS-CT-CT (37x; epoxy_amine_primary: DGEBA + DDM (aromatic amine), epoxy_amine_primary: DGEBA + DDS (aromatic sulfone amine), epoxy_amine_primary: DGEBA + DETA)
- CM-CT-CT-CZ (3x; C003 Acrylonitrile–butadiene–styrene terpolymer (random, 1 : 1 : 2 units), C105 NBR (33 wt % acrylonitrile))
- CM-CT-CT-S (19x; polysulfide_allylic: BR + TESPD, polysulfide_allylic: NR + TESPT, sulfur_allylic: BR + H2S2)
- CT-CM-CT-S (15x; polysulfide_allylic: NR + TESPT, sulfur_allylic: NR + H2S2)
- CT-NT-CT-CY (7x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- CT-OS-CT-CY (2x; CCCCOCC1CO1, epoxy_amine_primary: butyl glycidyl ether + ethylenediamine)
- CY-CY-CT-NT (7x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- HC-CY-CT-NT (7x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- NT-CT-CT-OS (8x; CC(N)COCC(C)OCC(C)N, epoxy_amine_primary: DGEBA + Jeffamine D230 (n=2))
- NT-CT-CY-O$ (7x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- OH-C-CT-OH (1x; CC(O)C(=O)O)
- S-CT-CT-S (4x; sulfur_allylic: BR + H2S2)

## Failed structures

- none
