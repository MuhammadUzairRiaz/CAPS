---
title: Force fields and engine inputs
lede: Type a structure with OPLS-AA or PCFF, inspect the assignment, and write complete inputs for LAMMPS and GROMACS.
time: 5 minutes
level: first steps
uses: ff info · ff type · ff apply
---

## The library

CAPS ships its force fields as JSON files in `data/forcefields`: OPLS-AA (2005, 2020, 2024), PCFF and COMPASS (class II), CVFF, GAFF and GAFF2, AMBER, CHARMM, DREIDING, TraPPE, UFF, inorganic and clay sets, Martini and others. Commands accept a library name or a path.

```bash
caps ff info opls2005        # types, typing rules, styles and references
caps ff info pcff
```

## 1. Assign types

Typing rules are SMARTS patterns: each atom gets the type of the most specific rule it matches. `caps ff type` writes one type per atom and reports atoms left untyped or matched ambiguously.

```bash
caps ff type ps_md.data --ff opls2005 -o types.txt
```

```text
0 untyped, 0 ambiguous
wrote types.txt
```

`--explain` prints, for every atom, the rule that won and the ones it beat; use it when a type looks wrong.

## 2. Parameters and engine files

`caps ff apply` looks up every bond, angle, dihedral and improper and the pair coefficients, assigns charges (with `--charges auto`, the force field's own scheme: fixed charges per type or bond increments), and writes the files for each engine:

```bash
caps ff apply ps_md.data --ff opls2005 --types types.txt -o ps_opls.data \
  --lammps-input ps_opls.in --lammps-run npt --steps 5000 --gromacs ps_opls
```

```text
wrote ps_opls.in (the force field's own styles)
gromacs: CAPS's damped shifted force becomes PME in GROMACS (no DSF there)
gromacs: the cell is 23.8 Å across, less than 2 r_c (24 Å): GROMACS needs a larger cell (a supercell) or a shorter cut-off
wrote ps_opls.top, ps_opls.gro and ps_opls.mdp
```

The notes matter: this tutorial's cell is smaller than twice the cut-off, which GROMACS refuses. A production cell (thousands of atoms) passes. Run the results with:

```bash
lmp -in ps_opls.in
gmx grompp -f ps_opls.mdp -c ps_opls.gro -p ps_opls.top -o run.tpr && gmx mdrun -deffnm run
```

## Class II force fields

PCFF and COMPASS carry cross terms (bond–bond, bond–angle, angle–angle, middle-bond and end-bond torsion couplings) and 9-6 Lennard-Jones pairs. The same two commands handle them; the LAMMPS input uses `class2` styles:

```bash
caps ff type ps_md.data --ff pcff -o types_pcff.txt
caps ff apply ps_md.data --ff pcff --types types_pcff.txt -o ps_pcff.data --lammps-input ps_pcff.in
```

## Options worth knowing

| Option | Effect |
|---|---|
| `--charges auto\|keep\|types\|gasteiger` | where charges come from (`keep`: the file's) |
| `--kspace auto\|pppm\|ewald\|dsf\|cut` | long-range electrostatics in the LAMMPS input |
| `--lammps-style native\|exact` | the force field's own styles, or the styles that give exactly CAPS's energy (to check CAPS against LAMMPS) |
| `--overlay USER.json` | your own parameters on top of the library's (missing terms, refits) |
| `--dlpoly DIR` | DL_POLY CONTROL, FIELD and CONFIG |

## Next

- [Tutorial 04](crosslinking.html) assigns a force field from Python and keeps it through a cure.
- [Theory: Force fields](../theory/force-fields.html) for the functional forms and mixing rules.
