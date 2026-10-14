# Pilot: GPU4PySCF (Alpha, A100) vs ORCA 6.0.1 (Barnard, 8 cores)

Level: wB97M-D3(BJ)/def2-TZVPPD, energy and gradient.

| fragment | atoms | dE GPU-ORCA (kcal/mol) | force RMS diff | force max diff (kcal/mol/A) | GPU s | ORCA s |
|---|---|---|---|---|---|---|
| 017cb631254b_0 | 14 | +0.097 | 0.019 | 0.051 | 39 | 80.77 |
| 04a86f7325d8_0 | 20 | +0.118 | 0.025 | 0.059 | 40 | 202.89 |
| 0c21215f9fbc_0 | 18 | +0.182 | 0.025 | 0.110 | 6 | 109.23 |
| 12f9bb3a1aa1_0 | 15 | +0.672 | 0.736 | 2.672 | 8 | 87.41 |
| 1432e8b58b5f_0 | 17 | +0.165 | 0.018 | 0.059 | 6 | 109.02 |
| 1e61496e0218_0 | 21 | +0.265 | 0.029 | 0.111 | 12 | 329.37 |
| 25f69b4b786e_0 | 19 | +0.146 | 0.025 | 0.113 | 7 | 157.95 |
| 26c2b9aa445d_0 | 27 | +0.247 | 0.211 | 0.879 | 14 | 539.84 |
| 28b4350afaa5_0 | 12 | +0.147 | 0.035 | 0.107 | 4 | 72.61 |
| 2d9204ae7408_0 | 21 | +0.260 | 0.287 | 1.186 | 13 | 431.41 |
| 5436f7288227_0 | 15 | +0.210 | 0.053 | 0.170 | 5 | 80.15 |
| 5c152cc1c0bc_0 | 21 | +0.162 | 0.022 | 0.076 | 6 | 155.3 |
| 6169ddde3018_0 | 21 | +0.173 | 0.019 | 0.047 | 7 | 160.23 |
| 80aadeb03aec_0 | 9 | +0.101 | 0.024 | 0.058 | 3 | 47.81 |
| 9acf7b60a1eb_0 | 8 | +0.096 | 0.024 | 0.067 | 9 | 35.95 |
| b5a4db729a66_0 | 16 | +0.134 | 0.023 | 0.051 | 6 | 112.67 |
| bafc1feb3458_0 | 20 | +0.242 | 0.211 | 0.829 | 9 | 298.74 |
| d6a9403e9f8e_0 | 20 | +0.415 | 0.052 | 0.186 | 12 | 184.98 |
| ecc5383e426f_0 | 20 | +0.163 | 0.020 | 0.081 | 6 | 120.43 |
| faae5059d9bf_0 | 23 | +0.187 | 0.034 | 0.147 | 9 | 225.23 |

Energy offset GPU - ORCA: mean +0.209, spread (SD) 0.129 kcal/mol over 20 points.
Forces: RMS difference median 0.025, largest component 2.672 kcal/mol/A.

## Relative conformer energies (GPU, kcal/mol)

- 1e61496e0218 (missing, sulfur_allylic: NR + H2S2): 0.00, 0.95, 0.15
- 2d9204ae7408 (CAPS fit, P071 Polyethersulfone): 0.00
- 12f9bb3a1aa1 (published, polysulfide_allylic: NR + TESPT): 1.90, 1.67, 0.00
- 5436f7288227 (published, P201 Polydimethylsiloxane (silicone rubb): 0.20, 0.00, 0.39
- b5a4db729a66 (CAPS fit, P097 Polychloroprene): 0.00, 3.98, 3.31
- 04a86f7325d8 (CAPS fit, P007 Polyvinylidene fluoride): 0.00, 2.81, 0.90
- 017cb631254b (published, P003 Polyvinyl fluoride): 0.00, 2.04, 0.92
- bafc1feb3458 (CAPS fit, P074 Polyaniline): 0.00
- 0c21215f9fbc (published, P017 Nylon 6,6): 0.45, 2.49, 0.00
- 28b4350afaa5 (missing, CC(O)C(=O)O): 4.16, 5.88, 0.00
- 25f69b4b786e (missing, C003 Acrylonitrile–butadiene–styrene ter): 0.00, 0.98, 3.12
- d6a9403e9f8e (missing, sulfur_allylic: BR + H2S2): 0.00, 51.23, 62.81
- 80aadeb03aec (CAPS fit, OC(=O)/C=C\C(=O)O): 0.00, 6.56, 7.03
- 6169ddde3018 (CAPS fit, C106 Butyl rubber (IIR, 2 mol % isoprene): 0.18, 0.38, 0.00
- 9acf7b60a1eb (published, CC(=O)O): 0.00, 0.19, 0.18
- 5c152cc1c0bc (published, P205 1,2-Polybutadiene (vinyl)): 0.00, 1.39, 2.98
- ecc5383e426f (published, P023 Poly(vinyl octyl ether)): 0.00, 0.63, 0.87
- 1432e8b58b5f (published, P034 Poly(butyl acrylate)): 0.00, 0.25, 8.72
- 26c2b9aa445d (CAPS fit, P073 Polysulfone): 0.00
- faae5059d9bf (published, P073 Polysulfone): 0.03, 0.07, 0.00

## Cost

- GPU: 10 s per conformer on average (SCF + gradient, one A100), 54 conformers.
- ORCA: 177 s per point on 8 cores.
- Full training set, 55,000 conformers: about 151 GPU-hours, 5 h of wall time on 32 A100s.
