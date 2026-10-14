# Torsion census: oplsaa2024

Structures typed: 211 (21 copolymer, 114 polymer, 38 product, 38 reactant); failed: 0. Environment radius: 2 bonds.

| | published | CAPS fit | missing | total |
|---|---|---|---|---|
| torsion classes | 249 | 78 | 20 | 347 |
| rotor environments | 1122 | 173 | 36 | 1331 |
| of them heavy-atom rotors | 859 | 156 | 31 | 1046 |
| bond classes | - | 0 | 1 | - |
| angle classes | - | 6 | 5 | - |

## Campaign size

- One relaxed scan per rotor environment, 24 points each (15 degrees).
- Library (environments the published table does not cover): 209 scans, 5016 QM points.
- Training set (every environment, published ones included): 1331 scans, 31944 QM points before distorted conformers are added.
- Fragment size: median 8 heavy atoms, largest 85.
- Fragments with a problem (left out of a campaign): 7 (6 an sp3 carbon angle of 85 deg, 1 over 80 atoms).
- Environments seen once: 784 (many are chain ends or one crosslink site; the ones seen often are the backbones).

## Missing torsion classes

- C-CT-OH-HO (6x; CC(O)C(=O)O, ester_condensation: lactic acid (PLA))
- C-OS-C-CB (2x; O=C1OC(=O)c2ccccc21)
- CA-CT-CT-CM (26x; C002 Styrene–butadiene rubber (SBR, 23.5 wt % styrene), C108 E-SBR (23.5 wt % styrene; butadiene 70 % trans, 15 % cis, 15 % vinyl), CC/C=C\CC(c1ccccc1)C/C=C\C)
- CA-CT-CT-S (5x; sulfur_allylic: SBR + H2S2)
- CA-NT-CT-CT (7x; epoxy_amine_primary: DGEBA + DDM (aromatic amine), epoxy_amine_primary: DGEBA + DDS (aromatic sulfone amine), epoxy_amine_primary: TGDDM + DDS)
- CA-NT-CT-CY (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- CM-CT-CT-CZ (3x; C003 Acrylonitrile–butadiene–styrene terpolymer (random, 1 : 1 : 2 units), C105 NBR (33 wt % acrylonitrile))
- CM-CT-CT-S (19x; polysulfide_allylic: BR + TESPD, polysulfide_allylic: NR + TESPT, sulfur_allylic: BR + H2S2)
- CT-CM-CT-S (15x; polysulfide_allylic: NR + TESPT, sulfur_allylic: NR + H2S2)
- CT-NT-CT-CY (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- CT-OS-CT-CY (5x; CCCCOCC1CO1, epoxy_amine_primary: butyl glycidyl ether + ethylenediamine)
- CT-S-S-HS (2x; sulfur_allylic: BR + H2S2, sulfur_allylic: SBR + H2S2)
- CY-CY-CT-NT (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- HC-CY-CT-NT (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- HS-S-S-HS (5x; SS, sulfur_allylic: NR + H2S2, sulfur_allylic: SBR + H2S2)
- NT-CT-CT-OS (8x; CC(N)COCC(C)OCC(C)N, epoxy_amine_primary: DGEBA + Jeffamine D230 (n=2))
- NT-CT-CY-O$ (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- OH-C-CT-OH (5x; CC(O)C(=O)O, ester_condensation: lactic acid (PLA))
- OH-CT-C-OS (1x; ester_condensation: lactic acid (PLA))
- S-CT-CT-S (5x; sulfur_allylic: BR + H2S2, sulfur_allylic: NR + H2S2)

## Failed structures

- none
