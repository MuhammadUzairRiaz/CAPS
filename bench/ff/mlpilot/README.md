# Pilot: GPU DFT for CAPS's own training set (2026-10-14)

20 census fragments with 54 conformers at wB97M-D3(BJ)/def2-TZVPPD, energies and forces:
GPU4PySCF 1.8.1 on Alpha (one A100 per task) against ORCA 6.0.1 on Barnard (8 cores, RIJCOSX, DefGrid3).
Files: make_pilot.py, run_gpu4pyscf.py + alpha_pilot.sbatch, orca_inputs.py, compare.py; the numbers are in
pilot/compare.md. Inputs, outputs, caches and logs ran in the horse workspace /data/horse/ws/USER-caps-ml.

## Result

- 16 of 19 fragments: forces agree to 0.018-0.053 kcal/mol/A RMS, energies to a near-constant offset of
  +0.10 to +0.27 kcal/mol (density fitting and grids; it cancels in relative energies).
- Cost: about 10 s per conformer on one A100 (SCF + gradient), against 36-431 s on 8 cores with ORCA. A
  training set of 55,000 conformers is about 150 GPU-hours.

## The three outliers were bad fragments, not the GPU code

The worst (12f9bb3a1aa1, a TESPT product) was rerun with ORCA without any integral approximation (NoRI), and on
the GPU with a denser grid (150, 974) and without density fitting:
GPU energy against NoRI +0.002 kcal/mol; GPU forces 0.72 kcal/mol/A RMS off whatever the grid or fitting,
concentrated on one C-H pair. The fragment file had two carbons written as hydrogens (the census took
elements from data-file masses). The other two (polyethersulfone, polyaniline) had two cap hydrogens 0.54 A
apart, where the radius cut went through a ring.

Over the whole census that was 4 fragments with a mislabelled element, 162 with colliding caps and 342 with an
odd electron count (a cut through a double or aromatic bond). caps.census now cuts its own fragments: small
rings kept whole, no bond between two unsaturated atoms cut, elements from the typed atoms, positions made
whole by minimum image in periodic products, and every fragment checked (even electron count, no atoms
closer than 0.9 A).

ORCA 6.0.1 refuses RIJK gradients for range-separated hybrids ("the range sep RIJK gradient has an unresolved
bug"), so the RIJK leg of the check could not run; NoRI is the reference.
