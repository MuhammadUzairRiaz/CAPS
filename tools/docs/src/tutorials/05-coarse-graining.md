---
title: Coarse-graining a polyester
lede: Map an all-atom poly(butylene succinate) run onto beads, invert bonded and non-bonded potentials, build a long-chain melt, measure entanglements, and put the atoms back.
time: 15 minutes
level: advanced
uses: cgmap · cgfit · cgbuild · ppa · cgdyn · mech · backmap
---

## The route

The coarse-graining commands follow the order of a study, and each writes the files the next one reads:

| Stage | Command | Reads | Writes |
|---|---|---|---|
| Map | `cgmap` | all-atom structure (+ dump) | `STEM.map.json`, `STEM.cg.data`, `STEM.cg.lammpstrj`, `types.json` |
| Bonded | `cgfit bonded` | maps + CG frames | `bonded/` tables, `bonded.in`, `bonded.json` |
| Non-bonded | `cgfit targets`, `ibi-start`/`ibi-step`/`ibi-run`, `fit`, `calibrate` | maps + frames, targets | pair tables, LAMMPS decks, `run_ibi.sh`, fits |
| Build | `cgbuild` | `bonded.json`, maps | `melt.cg.data`, `in.cg_equil` |
| Analyse | `ppa`, `mech`, `cgdyn` | CG frames | N\(_e\), stress–strain, g\(_1\)–g\(_3\) |
| Backmap | `backmap` | reference all-atom cell + CG frame | all-atom data, `in.backmap` |

## 0. An all-atom reference run

A small, short run is enough to follow the route; a real model needs an equilibrated melt with many frames.

```bash
caps grow -o pbs.data --chains 6 --dp 6 --density 0.6 --units '*OCCCCOC(=O)CCC(=O)*' --seed 1
caps relax pbs.data -o pbs_min.data --iterations 500 --quiet
caps md pbs_min.data -o pbs_md.data --steps 3000 --dt 1 --temp 450 --dump pbs.lammpstrj --every 100 --quiet
```

## 1. Map

The `ester-cut` preset cuts every ester C(=O)–O bond, so the fragments are the beads: B for the butanediol part, S for the succinate. Beads are named by rules; a fragment no rule names is reported, never merged silently.

```bash
caps cgmap pbs_md.data --preset ester-cut -o cg --dump pbs.lammpstrj
```

```text
    B S B S B S B S B S B S
    note: cut [CX3](=O)-[OX2;!H1]: 66 bonds
    trajectory: 31 frames → cg/pbs_md.cg.lammpstrj (0.0 s)
types: 2 bead, 1 bond, 2 angle, 1 dihedral
```

Several structures mapped together share one type list, so PBS, PBSA and PBAT models can share tables.

## 2. Bonded potentials

Boltzmann inversion of the bond, angle and dihedral distributions, with the Jacobians removed:

$$
U(r) = -k_\mathrm{B}T \ln\frac{P(r)}{r^2}, \qquad U(\theta) = -k_\mathrm{B}T \ln\frac{P(\theta)}{\sin\theta}, \qquad U(\varphi) = -k_\mathrm{B}T \ln P(\varphi)
$$

```bash
caps cgfit bonded cg/pbs_md.map.json cg/pbs_md.cg.lammpstrj -T 450 -o bonded
```

Each distribution is inverted only where it exceeds a share of its maximum, with slope-continuous walls outside. The two halves of the molecules are inverted separately as a check, and a NOT CONVERGED note says when they disagree; with this tutorial's 31 frames they do, which is the right warning.

## 3. Non-bonded potentials

Target radial distribution functions per bead pair, then iterative Boltzmann inversion:

$$
U_{i+1}(r) = U_i(r) + k_\mathrm{B}T \ln\frac{g_i(r)}{g_\mathrm{target}(r)}
$$

```bash
caps cgfit targets cg/pbs_md.map.json cg/pbs_md.cg.lammpstrj -T 450 -o targets
caps cgfit ibi-start cg/pbs_md.map.json --targets targets/targets.json --bonded bonded -T 450 -o ibi
```

`ibi-start` writes the starting tables, the LAMMPS decks and `run_ibi.sh`, which runs the loop on a cluster (`bash ibi/run_ibi.sh 20`). For small systems, `ibi-run` runs a few iterations with CAPS's own engine:

```bash
caps cgfit ibi-run cg/pbs_md.map.json cg/pbs_md.cg.data --targets targets/targets.json \
  --bonded bonded -T 450 --iterations 2 --steps 1000 -o ibirun
```

All pairs are updated jointly, with a pressure correction ramp. `cgfit fit` then fits Lennard-Jones 12-6, 9-6, Morse or Mie forms to the tables, and `cgfit calibrate` scales σ and ε towards a target density and T\(_g\).

## 4. Build a long-chain melt

```bash
caps cgbuild --units BS=B+S --dp 20 --chains 10 --density 1.2 --bonded bonded \
  --maps cg/pbs_md.map.json -T 450 -o melt cg/pbs_md.map.json cg/pbs_md.cg.lammpstrj
```

```text
  ⟨R²(n)⟩/n at n =  10: built 55.1 Å², all-atom 53.1 Å² (ratio 1.04)
wrote melt/melt.cg.data, .map.json, .internal.csv and melt/in.cg_equil
```

Chains are random walks drawn from the inverted bonded distributions, with real sequence statistics for copolymers (`--composition`, `--sequence markov`). Given the all-atom maps and frames, it compares the internal distances with the all-atom model's. `in.cg_equil` equilibrates the melt in LAMMPS: soft push-off, the model's pairs, a hot anneal, NPT.

## 5. Entanglements, mechanics, dynamics

```bash
caps ppa melt/melt.map.json melt/melt.cg.data -o ppa
caps cgdyn cg/pbs_md.map.json cg/pbs_md.cg.lammpstrj --dt 1 -o dyn
caps mech decks --rates 1e-6 --mode stress -o tension
```

`ppa` finds primitive paths (CAPS's own PPA, a LAMMPS deck with `--method lammps`, or Z1+) and reports the entanglement length by the estimators of Hoy, Foteinopoulou and Kröger (2009), for example the modified S-coil:

$$
N_e = (N - 1)\left(\frac{\langle L_\mathrm{pp}^2 \rangle}{\langle R^2 \rangle} - 1\right)^{-1}
$$

Run on an equilibrated melt; the unrelaxed build above only shows the output. `cgdyn` gives g\(_1\), g\(_2\), g\(_3\), the diffusion coefficient and τ\(_R\); `cgdyn timemap` maps CG time onto all-atom time. `mech decks` writes tension inputs; `mech analyze` gives the modulus, yield, softening and strain-hardening modulus.

## 6. Backmap

```bash
caps backmap cg/pbs_md.map.json pbs_md.data --cg melt/melt.map.json --frame melt/melt.cg.data -o aa
```

```text
400 beads → 4820 atoms · 4810 bonds, 8800 angles, 11580 dihedrals, 1200 impropers · net charge 0.0000 e
next: lmp -in aa/in.backmap -var DATA backmapped.data -var OUT relaxed (from aa), then caps backmap check aa/relaxed.data
```

Each bead gets a fragment of its class (head, inner or tail) from the reference cell, turned to its chain neighbours. The cut bonds and every term across them are restored with their own types and charges. `in.backmap` relaxes in stages, and `caps backmap check` compares bond lengths and angles with the force field's.

## Next

- The Studio's Coarse-grain page runs the same stages with a plot for each and the command line to copy.
- [Theory: Coarse-grained](../theory/coarse-grained.html) has every equation and its reference.
