# Torsion census: oplsaa2024

Structures typed: 211 (21 copolymer, 114 polymer, 38 product, 38 reactant); failed: 0. Environment radius: 2 bonds.

| | published | CAPS fit | missing | total |
|---|---|---|---|---|
| torsion classes | 256 | 75 | 14 | 345 |
| rotor environments | 1107 | 152 | 22 | 1281 |
| of them heavy-atom rotors | 861 | 135 | 20 | 1016 |
| bond classes | - | 0 | 0 | - |
| angle classes | - | 5 | 4 | - |

## Campaign size

- One relaxed scan per rotor environment, 24 points each (15 degrees).
- Library (environments the published table does not cover): 174 scans, 4176 QM points.
- Training set (every environment, published ones included): 1281 scans, 30744 QM points before distorted conformers are added.
- Fragment size: median 8 heavy atoms, largest 65.
- Fragments with a problem (left out of a campaign): 10 (7 two atoms 0.70 A apart, 3 two atoms 0.86 A apart).
- Environments seen once: 738 (many are chain ends or one crosslink site; the ones seen often are the backbones).

## Missing torsion classes

- C-CT-OH-HO (6x; CC(O)C(=O)O, ester_condensation: lactic acid (PLA))
- C-OS-C-CB (2x; O=C1OC(=O)c2ccccc21)
- CA-CT-CT-CM (26x; C002 Styrene–butadiene rubber (SBR, 23.5 wt % styrene), C108 E-SBR (23.5 wt % styrene; butadiene 70 % trans, 15 % cis, 15 % vinyl), CC/C=C\CC(c1ccccc1)C/C=C\C)
- CA-NT-CT-CT (7x; epoxy_amine_primary: DGEBA + DDM (aromatic amine), epoxy_amine_primary: DGEBA + DDS (aromatic sulfone amine), epoxy_amine_primary: TGDDM + DDS)
- CA-NT-CT-CY (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- CM-CT-CT-CZ (3x; C003 Acrylonitrile–butadiene–styrene terpolymer (random, 1 : 1 : 2 units), C105 NBR (33 wt % acrylonitrile))
- CT-NT-CT-CY (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- CT-OS-CT-CY (4x; CCCCOCC1CO1, epoxy_amine_primary: butyl glycidyl ether + ethylenediamine)
- CY-CY-CT-NT (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- HC-CY-CT-NT (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- NT-CT-CT-OS (8x; CC(N)COCC(C)OCC(C)N, epoxy_amine_primary: DGEBA + Jeffamine D230 (n=2))
- NT-CT-CY-O$ (10x; C(c1ccc(N(CC2CO2)CC2CO2)cc1)c1ccc(N(CC2CO2)CC2CO2)cc1, epoxy_amine_primary: TGDDM + DDS)
- OH-C-CT-OH (5x; CC(O)C(=O)O, ester_condensation: lactic acid (PLA))
- OH-CT-C-OS (1x; ester_condensation: lactic acid (PLA))

## Failed structures

- none
