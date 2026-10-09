---
title: Crosslinking a rubber network
lede: Cure natural rubber with a peroxide template inside CAPS, or export the same chemistry as a LAMMPS fix bond/react set with typed templates.
time: 10 minutes
level: intermediate
uses: Python · react · bond_react · export_engines
---

## Two ways to crosslink

**CAPS React** cures the structure itself, cycle by cycle: in each cycle it finds reactive pairs within a capture distance, makes the bonds the template describes, re-types the products and relaxes. It is quick and deterministic, and well suited to building a network at a target conversion.

**fix bond/react export** writes the same chemistry for LAMMPS's REACTION package: pre- and post-reaction templates and map files for every chemical environment found in your structure, typed with your force field, and an input that runs the cure under dynamics.

Both are driven from Python.

## 1. Build and type

```python
import caps

d = caps.polymer("*CC(C)=CC*", dp=10, chains=6, density=0.6, seed=1)   # cis-1,4-polyisoprene
d.relax(max_iterations=300)
d.field.assign("opls2005")
```

`caps.reaction_templates()` lists the built-in chemistry: C–C crosslinks, accelerated sulfur and peroxide cures of diene rubbers, silane polysulfides, epoxy–amine, epoxidised natural rubber with acids, anhydride–alcohol, esterification. `caps.reaction_template("peroxide_allylic")` prints one, with a comment on every atom.

## 2a. Export for LAMMPS fix bond/react

```python
r = d.bond_react("peroxide_allylic", "xl_lammps", forcefield="opls2005")
print(r["notes"])
```

```text
6 templates for 1 reactions cover 5572 of 26408 candidate pairs (within the survey distance) in this structure
reaction peroxide_allylic: 1802 rarer environments (20836 candidate pairs) have no template; raise the number of templates to cover them
stabilize_steps 200 on every reaction (the reacting atoms stay under nve/limit that long after a reaction)
```

The folder holds `react.data`, `react.in` and, per template, `_pre.mol`, `_post.mol` and `_map.txt`. Every template is checked against LAMMPS's rule that each atom must be reachable from an initiator without passing through an edge atom; hydrogens that change partner get a distance constraint in the map file. Options include `targets=[0.2, 0.4]` (stop and save at each crosslink density), `variants` (more templates for rarer environments), `type_groups` (atom types split by component, for blends) and `keep_chain_ids=True`.

```bash
cd xl_lammps && lmp -in react.in
```

## 2b. Or cure inside CAPS

```python
report = d.react("peroxide_allylic", cycles=5, per_cycle=4, target=0.3, seed=1)
print(report)
d.save("nr_xl.data")
```

```text
link 10 (cycle 5, peroxide_allylic): chain 3 unit 8 — direct bond — chain 6 unit 10
link 11 (cycle 5, peroxide_allylic): chain 2 unit 6 — direct bond — chain 6 unit 8
20 byproduct molecules removed
force field after the run: OPLS-AA / OPLS 2005 · complete
```

Stop at a conversion (`target`), or at a network measure: `crosslinks=("per_chain", 1.5)`, `("density", mol/m³)` or `("mc", g/mol)`. `between_chains=True` refuses links within a chain, and `keep_byproducts=True` keeps the H₂ or H₂O that leaves.

## 3. Engine inputs for the network

```python
d.export_engines("engines", stem="nr", run="npt", temperature=300, steps=5000)
```

writes `nr.data` and `nr.in` for LAMMPS and `nr.top`, `.itp`, `.gro`, `.mdp` and `nr_em.mdp` for GROMACS, with the force field re-assigned to the products.

## Next

- `d.react_summary()` returns the network numbers of the cure (links per chain, crosslink density, M\(_c\)); `caps analyze FILE --props crosslinks` counts sulfur bridges after a sulfur cure.
- `d.react_sites("peroxide_allylic")` counts each chain's reactive sites before you choose a target.
