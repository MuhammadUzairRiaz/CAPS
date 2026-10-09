---
title: The caps command line
eyebrow: Command line · guide
lede: How the `caps` command is organised, the conventions every command shares, which command to reach for at each step of a project, and the scripts CAPS writes for you to run on a workstation or a cluster.
---

## First steps

Run `caps` with no arguments to see every command and its main options on one screen. The workflow commands (coarse-graining and DFT surfaces) also take `--help`, which prints every option with its default and the reason for that default:

```bash
caps                      # the usage: every command
caps cgfit --help         # one command's options, defaults and reasons
caps info ps.data          # any structure: atoms, molecules, bonds, cell, frames
```

`caps info` is the quickest check that a file reads as you expect. Sample structures and recipes ship in the program's `samples` folder.

> The force-field library, typing rules, presets and samples ship with the program and are found beside it. Set `CAPS_HOME` to a folder holding `data/` only to use another copy (a development checkout, say).

## Conventions every command shares

### Units

CAPS works in LAMMPS *real* units throughout: lengths in Å, time in fs (time steps) or ps (run lengths, where an option says so), temperature in K, pressure in atm, energy in kcal/mol, density in g/cm³, charge in e. Options name their unit in the usage when it is not one of these.

### Files

The format follows the extension, both ways:

| Kind | Read | Written |
|---|---|---|
| Structures | `.data` (LAMMPS), `.pdb`, `.xyz`, `.mol2`, `.sdf`/`.mol`, `.cif`, `.car`/`.mdf`, `.gro` + `.top`, `.vasp`/`POSCAR`/`CONTCAR`, AMBER `prmtop` with coordinates | the same, chosen by `-o` |
| Trajectories | LAMMPS dump (`.lammpstrj`), `.dcd`, `.xtc`, `.trr`, multi-frame `.xyz`, AMBER NetCDF / `mdcrd` | `caps frames` writes `.lammpstrj`, `.dcd`, `.xyz`, `.pdb`, `.gro`, `.trr` |
| Force fields | CAPS force-field JSON (`data/forcefields`), or a library name such as `opls2005`, `pcff`, `compass` | LAMMPS data + input, GROMACS `.top`/`.itp`/`.gro`/`.mdp`, DL_POLY |
| Recipes, pipelines | YAML or JSON | — |

A trajectory carries no bonds, so commands that read one take `--topology DATA` (the data file of the same atoms). `-o` names the output; a folder for the commands that write several files.

### Seeds and threads

Everything random takes `--seed` (default 1): the same command, input and seed give the same result on the same build. The commands that run in parallel (`md`, `run`, `pack`) take `--threads N` (default: all cores).

### Provenance

The recipe runner and the coarse-graining and DFT commands write `FILE.provenance.json` beside what they make: the program's version, the inputs' SHA-256, every step with its parameters, seed and the methods it follows. `caps provenance FILE --methods` turns it into a methods paragraph with numbered references, `--bibtex` into BibTeX, and `--compare OTHER` lists what differs between two runs.

### Output and exit codes

Commands print a short report and the files they wrote. `--json` (recipes, coarse-graining and DFT commands) prints the result as JSON for scripts. A command exits 0 on success and non-zero on an error. `caps run` distinguishes 2 (bad input), 3 (parameters missing from the force field) and 4 (a failed run); `caps validate` exits 1 when a slab fails a check, so it can stop a batch script.

## Which command for what

Each command links to its full entry in the [command reference](reference.html), which is generated from the program itself.

<!-- COMMAND TABLES -->

## Scripts CAPS writes for you to run

Long runs belong on a workstation or a cluster, not inside CAPS. For those, CAPS writes ready-to-run inputs and scripts, each with a comment block at the top saying how to run it and which variables it takes.

| File | Written by | Run it with |
|---|---|---|
| `STEM.in` + `STEM.data` | `caps ff apply … --lammps-input STEM.in`, `caps run` | `lmp -in STEM.in` |
| `STEM.top`, `.itp`, `.gro`, `.mdp` (+ `STEM_em.mdp`) | `caps ff apply … --gromacs STEM`, `caps run` | `gmx grompp -f STEM.mdp -c STEM.gro -p STEM.top -o run.tpr` then `gmx mdrun -deffnm run` |
| `react.in`, `*_pre.mol`, `*_post.mol`, `*_map.txt` | Python `Document.bond_react(...)` | `lmp -in react.in` (LAMMPS with REACTION) |
| `bonded/bonded.in` + tables | `caps cgfit bonded` | `include` it after `read_data`, with `-var BONDED bonded` |
| `ibi/run_ibi.sh`, `in.cg_run`, `in.cg_tg` | `caps cgfit ibi-start` | `bash ibi/run_ibi.sh 20` (20 iterations; `CAPS_LMP` picks the LAMMPS binary) |
| `melt/in.cg_equil` | `caps cgbuild` | `lmp -in in.cg_equil -var DATA melt.cg.data -var BONDED ../bonded -var PAIR ../lj -var OUT equil` |
| `in.cg_ppa` + `STEM.ends.in` | `caps ppa … --method lammps` | `lmp -in in.cg_ppa -var DATA equil.data -var OUT ppa`, then `caps ppa MAP equil.data ppa.lammpstrj` |
| `tension/run_tension.sh`, `in.cg_tensile_*` | `caps mech decks` | `DATA=equil.data bash run_tension.sh`, then `caps mech analyze` |
| `aa/in.backmap`, `pair_coeffs.in` | `caps backmap` | `lmp -in in.backmap -var DATA backmapped.data -var OUT relaxed`, then `caps backmap check` |
| `job.slurm`, `make_potcar.sh`, `INCAR.*`, `KPOINTS` | `caps vasp-set`, `caps adsorb-dft` | `bash make_potcar.sh` (your PAW set), `sbatch job.slurm`; `caps vasp-jobs` submits and tracks many |

> **Careful:** the LAMMPS inputs for coarse-grained models use `dihedral_style table/cut` (LAMMPS's EXTRA-MOLECULE package) and fix bond/react needs the REACTION package. `lmp -h` lists the packages a LAMMPS binary has.

## Putting commands together

### A shell loop over seeds

Independent replicas differ only in the seed; the provenance records it, so `caps provenance --compare` shows exactly that difference.

```bash
for s in 1 2 3; do
  caps grow -o ps_$s.data --chains 4 --dp 10 --density 0.5 --units '*CC(*)c1ccccc1' --seed $s
  caps relax ps_$s.data -o ps_${s}_min.data --quiet
  caps md ps_${s}_min.data -o ps_${s}_md.data --steps 2000 --barostat berendsen --seed $s --quiet
done
caps analyze ps_1_md.data --props density
```

### A recipe instead of a script

When the chain is fixed (build, type, grow, relax, equilibrate, run, analyse, export), write it once as a recipe and run it with `caps run`. The recipe's hash, every stage's parameters and the seed go into the provenance, so the whole study can be rebuilt; see [Tutorial 03](../tutorials/recipes.html).

### Pipelines over many trajectories

`caps pipeline FILE --steps STEPS.yaml` runs a chain of analysis steps (the Studio's pipeline, including Python steps) on one frame; `caps run PIPELINE.yaml --input 'runs/*/traj.lammpstrj' --csv results.csv` runs it over many inputs and writes one row per input.

## Command line, Python or Studio?

All three call the same core, so results agree. Use the **command line** for batch work, shell scripts and clusters; **Python** (`import caps`) when a study needs loops, conditions or your own analysis, or features with no command yet (CAPS React's fix bond/react export, polymer sweeps, χ by MD); the **Studio** to look at structures, build interactively, and read plots and the methods report. The Studio's Coarse-grain and DFT pages show the equivalent command line for what you set up, ready to copy.
