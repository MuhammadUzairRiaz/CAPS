# CAPS's machine-learned potential

`caps-mace-v1.model` is CAPS's own MACE model. The C++ core evaluates it through libtorch (`caps ml eval FILE`, in a
build with `-DCAPS_WITH_TORCH=ON`), with no Python; `python -m caps.autofit` uses it to fit the bonded terms a force field
lacks, in place of a quantum-chemistry run (`--model caps-mace`, or `--model auto`, the default, for fragments of its
elements). It runs on a CPU: about 0.1 s per energy and force evaluation of a 20-atom fragment on an Apple M-series
Mac.

| | |
|---|---|
| File | `caps-mace-v1.model`: the compiled (TorchScript) MACE 0.3.16 model, float64, SHA-256 `ea149a99340660c91274f8974409c4319bdd46738ff252a02bac17870226d3f9` |
| Elements | H, C, N, O, F, Si, S, Cl; neutral, closed-shell molecules |
| Reference level | wB97M-D3(BJ)/def2-TZVPPD, density-fitted, GPU4PySCF 1.8.1, which reproduces ORCA 6.0.1 to 0.02-0.05 kcal/mol/A in force (`bench/ff/mlpilot`) |
| Training set | 48,420 conformers: relaxed torsion scans (AIMNet2, dihedral held) and distorted copies of 1,243 fragments from `caps.census`, covering the library polymers, copolymers, React curatives and cured networks. Assembled and audited by `bench/ff/mltrain/assemble.py` and `audit_trainset.py` |
| Model | hidden irreps 128x0e + 128x1o, r_max 5 A, 150 epochs (SWA from 120), one A100, 9.5 h |
| Test | 5,326 conformers of 131 fragments held out of training whole: energy RMSE 0.6 meV/atom, force RMSE 15.8 meV/A |
| Torsion profiles | 122 held-out fragments: RMSE median 0.15, 90th percentile 0.41, largest 0.92 kcal/mol; all within 1 kcal/mol; barrier error median 0.14 kcal/mol (AIMNet2 on the same conformers: 0.35 / 0.63 / 1.30, barrier 0.39) |

The test fragments come from the same chemistry as the training fragments (polymers, rubbers, curatives and crosslinks),
held out by fragment. Outside that chemistry, or for other elements, `--model auto` falls back to AIMNet2, and every
fitted term names the potential it came from.

Licence: the model is part of CAPS and under CAPS's licence (BSD 3-Clause, `LICENSE`). It was trained by CAPS on
CAPS's own DFT data with the MACE code (MIT, github.com/ACEsuit/mace); no third-party weights are in it. Running it
needs `mace-torch` (MIT), installed with pip, which is not part of CAPS.
