---
title: Recipes, provenance and a methods paragraph
lede: Describe a whole study in one YAML file, run it with one command, and get the methods paragraph, BibTeX and a record that shows exactly how two runs differ.
time: 5 minutes
level: first steps
uses: run · provenance
---

## A recipe

A recipe lists the stages of a study: build, type, grow, relax, run, analyse, export. This one ships as `samples/ps_cell.yaml`:

```yaml
recipe: 1
name: ps_cell
build:
  polymer:
    smiles: "*CC(*)c1ccccc1"
    dp: 20
    chains: 10
    tacticity: atactic
type: { forcefield: gaff2 }
grow: { trials: 16, density: 0.50, seed: 20260923 }
relax: { method: lbfgs, fmax: 0.5 }
md: { ps: 1, temperature: 300, ensemble: nvt }
analyze: { properties: [density, rg] }
export: [lammps, gromacs]
```

## Run it

```bash
caps run samples/ps_cell.yaml --out rec
```

The folder gets the structure (`ps_cell.data`), the LAMMPS input, the GROMACS topology and coordinates, the properties (`ps_cell.properties.json`) and `ps_cell.data.provenance.json`. Exit codes tell a script what went wrong: 2 a bad recipe, 3 parameters missing from the force field, 4 a failed run.

## The methods paragraph

```bash
caps provenance rec/ps_cell.data --methods
caps provenance rec/ps_cell.data --bibtex > refs.bib
```

```text
The structure was made by the CAPS recipe "ps_cell" (SHA-256 51565067c3c3…), whose steps follow. Atom types were
assigned from GAFF2 (moltemplate) [1], with Gasteiger charges [2]. 10 chains of 20 repeat units were grown in a
periodic cell at an initial density of 0.5 g/cm³ with CAPS 0.1.0, each unit placed from internal coordinates [3] …
```

Every method CAPS used cites its source, and the references are checked against their DOIs. Edit the paragraph as you like; the facts in it come from the record, not from memory.

## What differs between two runs

```bash
caps run samples/ps_cell.yaml --out rec2 --seed 2
caps provenance rec/ps_cell.data --compare rec2/ps_cell.data
```

```text
step 3 grow.trials            rng                mt19937-64 · seed 20260923  vs  mt19937-64 · seed 2
step 5 dynamics.nvt           rng                mt19937-64 · seed 1  vs  mt19937-64 · seed 2
```

`--seed` overrides every seed in the recipe, so replicas need no edited copies.

## Bundles

`caps bundle FILE --steps S.json -o fig.caps-bundle.zip --include-input` packs a figure's data with its analysis pipeline, provenance and hashes; `caps reproduce fig.caps-bundle.zip` rebuilds the data from the input and compares the SHA-256 of every file. Put a bundle beside each figure of a paper.

## Next

- The [CLI guide](../cli/index.html#putting-commands-together) shows recipes beside shell loops and pipelines.
- `caps bench --quick` runs CAPS's built-in validation suite against reference values.
