---
title: A DFT surface study
lede: Cut a Ti₃C₂ MXene sheet, terminate it with O, OH and F, validate it, place a nitrile on it, and write VASP case folders ready for a Slurm cluster.
time: 10 minutes
level: intermediate
uses: sheet · terminate · validate · adsorb-dft · vasp-set · vasp-jobs · vasp-analyze
---

## 1. A sheet

```bash
caps sheet --preset Ti3C2 -o Ti3C2.vasp
```

Presets come from `data/sheets/sheets.json` (each a list of layers: element, stacking position, height). To cut one complete layer from a bulk or MAX-phase structure instead, use `--from Ti3AlC2.cif --formula Ti3C2 --remove Al`.

## 2. Terminations

```bash
caps terminate Ti3C2.vasp --top O:0.5,OH:0.25,F:0.25 --supercell 3x3 --seed 1 -o mixed.vasp
```

Groups sit in fcc or hcp hollows, on top or bridge sites, at the height set by the bond length:

$$
h = \sqrt{d^2 - r^2}
$$

where d is the bond length to the nearest surface atoms and r the in-plane distance to them. Mixed faces are drawn by the fractions with the seed; the bottom face matches the top unless `--janus`.

## 3. Validate

```bash
caps validate mixed.vasp
```

A PASS/FAIL report: composition, one complete sheet and its layers, floating or isolated atoms, vacuum, close contacts, termination distances and coordination, O–H bonds, site classes and top/bottom asymmetry. It exits 1 on FAIL, so a batch script can stop on a broken slab.

## 4. A molecule on the slab

```bash
caps adsorb-dft mixed.vasp --smiles 'C#N' --out ads
```

CAPS finds anchors by SMARTS (here the nitrile N), places the molecule anchor-down, parallel and upright at several azimuths, lowers it to the contact limits and checks the distance to its periodic images. The `ads` folder holds `slab/`, `molecule/` and one `complex_*` folder per placement, all with the same settings so the adsorption energy is consistent:

$$
E_\mathrm{ads} = E_\mathrm{complex} - E_\mathrm{slab} - E_\mathrm{molecule}
$$

## 5. VASP case folders

```bash
caps vasp-set mixed.vasp --out vset --profile generic-slurm
```

The folder gets `POSCAR` (species grouped), `INCAR.cell`, `INCAR.relax` and `INCAR.static` with a reason beside every tag, `KPOINTS`, `job.slurm` (stages, resume, backup, a watchdog), `make_potcar.sh` and a validation report. POTCAR files are licensed by VASP, so CAPS never ships them: `--pp-dir` points at your PAW set, or run `make_potcar.sh` on the cluster.

## 6. On the cluster

```text
caps vasp-jobs submit|update|reset|cleanup|store    many cases at once
caps vasp-check · vasp-progress · vasp-health       what is running, converging or stuck
caps vasp-conv setup|collect CASE                   ENCUT and k-mesh convergence
caps vasp-analyze geom|wf|dos|cdd|bader|freq|md|summary
```

`vasp-analyze` reads the finished runs: geometry (heights, anchor distance, H-bonds, tilt), work function and its shift, grouped DOS, charge-density difference, Bader charges, the ν(C≡N) shift and AIMD stability; `summary` gathers them into `results_summary.csv`.

## Next

- [Theory: DFT surfaces](../theory/dft-surfaces.html): termination heights, one k-point density for every cell, the adsorption energy and the work function.
- Every DFT command takes `--help`, with each option's default and why.
