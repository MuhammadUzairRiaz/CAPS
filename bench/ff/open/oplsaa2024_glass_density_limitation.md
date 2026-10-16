# OPLS-AA 2024: glassy polymer densities (recorded 2026-10-16)

OPLS-AA 2024 under-predicts the room-temperature density of glassy polymers in CAPS. We record this as a limitation of
the force field and do not refit it: the evidence below rules out the cooling protocol and the backbone torsions as the
main cause, and a non-bonded refit to the polymers would undo the monomer liquid fits.

The run data lives in the paper's derivation folder, `docs/paper/runs/ffderive/` of the main checkout (not tracked).

## Densities at 298 K (g/cm³, 1 atm, about 12,000 atoms, DP 68-100)

| Polymer | Larsen 21 only | + slow anneal | Experiment | Gap after anneal | Melt vs Rodgers PVT |
|---|---|---|---|---|---|
| PMMA | 1.0830 | 1.0890 | 1.187 | -8.3 % | -4.2 to -5.1 % (525-475 K) |
| PS | 1.0031 | 1.0091 | 1.046 | -3.5 % | -3.0 to -3.2 % (525-475 K) |
| PVAc | 1.1376 | 1.1422 | 1.188 | -3.8 % | -1.5 to -1.7 % (450-400 K) |
| PET | 1.2501 | 1.2511 | 1.335 | -6.3 % | -1.9 % (600 K, inside the measured range) |

Experiment: `ffderive/polymers/polymer_reference.json`. Melt PVT: P. A. Rodgers, J. Appl. Polym. Sci. 48, 1061 (1993),
Table VIII (`ffderive/polymers/melt_pvt_reference.json`). Most melt holds lie 30-100 K above the measured range and use
the Tait V(0,t) extrapolated.

## What was ruled out

- **Cooling history.** `bench/ff/polymer_validate.py --anneal` (1 ns melt, then 25 K / 400 ps steps, 62.5 K/ns against
  Larsen's ~5,000 K/ns quench) gained only 0.1-0.6 %. Barnard jobs 39667881, 39667889, 39668041; results in
  `ffderive/polymers/anneal/`.
- **An unequilibrated melt.** The 1 ns melt holds at 525 K drift by 0.1 %/ns (PMMA, PS).
- **The atoms' Lennard-Jones terms as such.** The monomer liquids are right with the same atom types: methyl methacrylate
  +0.31 % (with `data/forcefields/oplsaa2024-lj-refit.json`), styrene +0.11 %, toluene +0.02 %, ethylbenzene -0.18 %
  (`ffderive/liquids/`).
- **The backbone torsions as the main cause.** MP2/cc-pVTZ//HF/6-31G** scans (`ffderive/torsions/pmma_ps_backbone_check/`):
  OPLS-AA 2024's own terms against QM, RMS 2.14 (PMMA CT-CT-CT-CT, 11 of 24 points), 1.62 (C-CT-CT-CT), 1.00 / 0.76
  (ester rotation), 1.00 (PS CT-CT-CT-CT) and 0.32 (PS CA-CT-CT-CT) kcal/mol. A separate PMMA dyad scan
  (`ffderive/pmma_dyad`, tf25) gave 1.86 / 1.78, and a joint refit of the shared classes only reached 1.87 / 1.48. These
  are real errors, but CT-CT-CT-CT is the generic alkane class (it fits polyethylene-like chains to 0.16-0.30), and
  torsion errors of this size change chain stiffness and Tg more than density. A PMMA- and PS-specific torsion refit was
  judged not worth it for density (expected gain about 1 %, not measured).

## What remains

The melts themselves are 3-5 % light for PMMA and PS, and about 2 % for PVAc and PET, while the monomers are right. The
gap is in how OPLS-AA packs these chains, not in any single term we could isolate. The rest of the glass gap beyond the
melt (PMMA about 3 %, PET about 4 %) follows from the higher simulated Tg and the glass line below it.

Report glassy-polymer densities from OPLS-AA 2024 with this caveat. Rubbery polymers (PE, PP, PCL) are within about 2 %.
