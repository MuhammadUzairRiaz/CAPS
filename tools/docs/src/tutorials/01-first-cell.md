---
title: Your first amorphous polymer cell
lede: Grow polystyrene chains in a periodic box, relax them, run molecular dynamics, and measure density, chain size and stiffness, using five commands.
time: 5 minutes
level: first steps
uses: grow · relax · md · analyze · check · render
---

## What you will make

A small atactic polystyrene melt: 4 chains of 10 repeat units (648 atoms), relaxed and run for 2 ps at 300 K and 1 atm. It is far too small and too short for production numbers; the point is the route, which is the same at any size.

## 1. Grow the chains

`caps grow` places repeat units one at a time from internal coordinates, choosing among trial torsions the one with the most room, so chains do not overlap even at a useful density. The repeat unit is a SMILES string with two `*` marking where units join.

```bash
caps grow -o ps.data --chains 4 --dp 10 --density 0.5 --units '*CC(*)c1ccccc1' --seed 1
```

Growing at 0.5 g/cm³ (about half polystyrene's density) leaves room; the barostat compresses the cell later. `--tacticity isotactic` or `syndiotactic` changes the stereochemistry, `--units` with several SMILES and `--sequence` makes copolymers.

## 2. Relax

Grown cells have strained bonds and close contacts. `caps relax` pushes overlapping atoms apart with capped forces, then minimises with L-BFGS until the largest force is below 0.5 kcal/mol/Å.

```bash
caps relax ps.data -o ps_min.data --iterations 500 --quiet
```

With no force field named, CAPS uses its default (GAFF for C and H, UFF otherwise) and writes it into the `.data` file, so LAMMPS can read the result directly.

## 3. Run molecular dynamics

```bash
caps md ps_min.data -o ps_md.data --steps 2000 --dt 1 --temp 300 --barostat berendsen \
  --dump ps.lammpstrj --every 200 --quiet
```

This is 2 ps of NPT dynamics with velocity Verlet, the Bussi thermostat (the default) and a Berendsen barostat, saving a frame every 200 steps. For a real study, compress and anneal first with `caps equilibrate --protocol larsen21` (the 21-step compression and decompression protocol), then run nanoseconds.

> Keep runs on a laptop short; for long runs, export the cell to LAMMPS or GROMACS ([Tutorial 02](force-fields.html)) and run it on a workstation or cluster.

## 4. Measure

```bash
caps analyze ps.lammpstrj --topology ps_md.data --props density,rg --json props.json
caps analyze ps.lammpstrj --topology ps_md.data --props rdf,ree,cn --pair C-C --inter --csv csv
```

A trajectory has no bonds, so `--topology` gives the data file of the same atoms. Each property prints its value, how it was computed and a standard error from block averages; `--csv` writes the curves (g(r), end-to-end distances, C\(_n\)) for plotting. The characteristic ratio comes from all internal distances along the backbones:

$$
C_n = \frac{\langle R^2(n) \rangle}{n\,\langle b^2 \rangle}, \qquad C_n = C_\infty\left(1 - \frac{a}{n}\right)
$$

With chains of 19 backbone bonds, CAPS notes that C\(_\infty\) is an extrapolation; longer chains give a value you can compare with experiment (about 10 for atactic polystyrene).

## 5. Check and look

```bash
caps check ps_md.data
caps render ps_md.data -o ps.png --style sticks --colour molecule --size 900x700
```

`caps check` reports counts, bonds, close contacts, charges and the cell; `caps render` draws a picture without opening the Studio. To look around interactively, open `ps_md.data` in CAPS Studio.

## Next

- [Tutorial 02](force-fields.html): the same cell with OPLS-AA, PCFF or COMPASS, written for LAMMPS and GROMACS.
- [Tutorial 03](recipes.html): the whole route as one recipe file, with a methods paragraph for your paper.
- [Theory: Growth](../theory/growth.html) and [Analysis](../theory/analysis.html) for the methods behind each step.
