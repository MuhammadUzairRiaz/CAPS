---
title: Analysing trajectories
lede: Density, structure, chain dimensions, dynamics, free volume and cohesive energy from any trajectory CAPS reads, with errors, curves for plotting, and figures.
time: 10 minutes
level: first steps
uses: analyze · rdf · chains · frames · render · pipeline
---

## One command, many properties

`caps analyze` computes any set of properties over the frames of a trajectory (from CAPS, LAMMPS, GROMACS or AMBER) and reports each with its method and a block-average error:

```bash
caps analyze ps.lammpstrj --topology ps_md.data --props density,rg,ree,cn --json props.json
caps analyze ps.lammpstrj --topology ps_md.data --props rdf,sq --pair C-C --inter --csv curves
```

| Property | What it gives |
|---|---|
| `density` | mass over cell volume |
| `rdf`, `sq`, `xray`, `neutron` | g(r) by element pair (`--pair C-C`, `--inter` for different molecules only), S(q), scattering |
| `rg`, `ree`, `cn`, `persistence` | radius of gyration, end-to-end distance, characteristic ratio C\(_n\) and C\(_\infty\), persistence length |
| `msd`, `diffusion`, `relaxation` | mean-square displacement, D from its slope, end-to-end relaxation |
| `ced`, `delta` | cohesive energy density and the Hildebrand solubility parameter |
| `ffv`, `psd` | fractional free volume and the pore-size distribution (`--probe`, `--grid`) |
| `crosslinks`, `entanglements` | network measures, primitive-path entanglements |
| `zprofile`, `adhesion`, `interaction`, `orientation` | interfaces (`--surface`, `--axis`) |

The diffusion coefficient comes from the Einstein relation in the mean-square displacement's linear range:

$$
D = \lim_{t\to\infty} \frac{\langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle}{6t}
$$

so give the time between frames: `--frame-ps 0.2` (or `--timestep-fs 1` when the dump records step numbers).

## Choosing frames

`--first`, `--last` and `--stride` pick frames; drop the equilibration part of a run with `--first`. `caps frames` writes a thinned or converted copy:

```bash
caps frames ps.lammpstrj thin.dcd --topology ps_md.data --stride 2
```

## Quick looks

```bash
caps rdf ps_md.data --pair C-C --inter       # g(r) of one frame
caps chains ps_md.data                       # backbones and mean-square internal distances
caps shape ps_md.data                        # per-molecule Rg and shape
caps render ps_md.data -o ps.png --style sticks --colour molecule
```

## Pipelines

For analysis you repeat, save the chain of steps (selections, modifiers, Python steps) as a pipeline, from the Studio or by hand, and run it on one frame or on many trajectories:

```bash
# steps.yaml: a pipeline saved from the Studio, or written by hand
caps pipeline ps_md.data --steps steps.yaml --out out
caps run steps.yaml --input 'runs/*/traj.lammpstrj' --csv results.csv
```

## Next

- [Theory: Analysis](../theory/analysis.html) gives every property's definition, defaults and tests.
- [Python](../python.html): `Document.analyze()` returns the same numbers as a dictionary.
