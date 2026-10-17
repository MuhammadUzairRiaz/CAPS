# Force-field validation

Every converted force field in `data/forcefields` is checked against the program it came from, interaction by
interaction, and every CAPS energy form is checked against LAMMPS.

## Tools

| Script | What it checks |
|---|---|
| `check_forms.cpp` (+ `check_forms_cmp.py`) | each CAPS functional form with random coefficients against the LAMMPS style it mirrors (class II, Morse / GROMOS bonds, cosine angles, Urey–Bradley, umbrella inversion, Buckingham, Morse pair, separate 1-4 LJ) |
| `compare_moltemplate.py` | a moltemplate conversion against moltemplate's output and LAMMPS |
| `compare_class2.py` | class II (COMPASS) conversions: every cross term, bond increments, LAMMPS CLASS2 energies / forces |
| `compare_msi_types.py` | automatic typing against Materials Studio's types in msi2lmp data files |
| `check_data_lammps.py` | LAMMPS data files written by CAPS, run in LAMMPS: every energy term and force, for each force-field family, hybrid styles included |
| `validate_rtf_types.py` | automatic typing against a CHARMM topology's residues (CGenFF's own model compounds); `--apply` parameterises them |
| `convert_frc.py` | CVFF, PCFF, COMPASS and IFF 1.5 from their `.frc` files with msi2lmp's semantics; IFF's silica and metal typing rules and the silica charges IFF states |
| `check_msi2lmp.py` | the `.frc` conversions against msi2lmp's data files (LAMMPS's test structures, IFF's model database), term by term in LAMMPS; then IFF's rules on its own silica and metal models read as bare coordinates (types, charges, bonds) |
| `check_gromacs.py` | GROMACS topologies written by CAPS, run in GROMACS (grompp, a zero-step mdrun) |
| `convert_martini3.py` | Martini 3 from cgmartini.nl's martini_v300.zip: bead types and the full pair table (`data/martini/martini3-nonbonded.tsv.gz`), the molecules of its solvent, ion, small-molecule, sugar, nucleobase and phospholipid files as templates |
| `check_martini3.py` | Martini 3 boxes against GROMACS (reaction field, virtual sites, restricted bending) |
| `vermouth_ff.py` | vermouth-martinize's `.ff` and `.map` readers, for the Martini 3 small molecules |
| `convert_vermouth_martini3_small.py` | Martini 3 small molecules: vermouth's blocks and mappings, the CHARMM residues they are written for |
| `check_martini3_small.py` | every Martini 3 small molecule built all-atom, recognised by graph and mapped |

LAMMPS with CLASS2, MOFFF, CORESHELL, MANYBODY (Stillinger-Weber) and CG-SPICA (lj/sdk, angle sdk) is built at
`~/lammps/build-class2/lmp`.

### Coarse-grained and many-body forms (`check_data_lammps.py`, DSF and PME)

| Case | Result |
|---|---|
| MARTINI 2, its overlays, SDK, Cooke-Deserno lipids, mW water | cases kept in `check_data_lammps.py` under neutral aliases (martini, sdk, cooke-deserno, mw); they run once those models are rebuilt from their papers |

### Martini 3 (`check_martini3.py`, against GROMACS)

Martini 3 has no LAMMPS form (reaction-field Coulomb, virtual sites): CAPS refuses the LAMMPS export and GROMACS is the
reference. Boxes of Martini 3's own molecules are built from the templates, packed, written to a data file and
recognised again, exported with `caps ff apply --gromacs` and run with `gmx mdrun -rerun`:

| Case | Result |
|---|---|
| Water and ions (W, Na+, Cl-, Ca2+; Na+ and Cl- are both TQ5, told apart by charge) | energy 2e-6 (relative), forces 3e-4 kcal/mol/Å |
| Lipids in water (POPC, POPS, DOPE, Na+) | energy 3.7e-6, forces 3e-4 |
| Nucleobases (virtual sites, exclusions), TXE (restricted bending), BIM / PCRE (harmonic impropers) | energy 2.4e-6, forces 7.6e-3 |
| Every one of the 217 molecule templates once, in water | energy 2.9e-6, forces 2.5e-2 |

| All 43 Martini 3 small molecules, all-atom (built by CAPS) mapped by graph | energy 1.1e-6 |

Forces on beads held by constraints (stiff bonds in CAPS, 1e6 kJ/mol/nm²) carry the single-precision position error of
GROMACS (~0.012 kcal/mol/Å per stiff bond in a 9 nm box); a mapped structure that is not relaxed can have beads under
1e5 kcal/mol/Å, which single precision carries to ~4e-5 of the force. The check allows those and nothing else. The
.gro holds the virtual sites where CAPS places them (GROMACS's rerun takes coordinates as written).

### Martini 3 small molecules (`check_martini3_small.py`)

`data/martini/martini3-small-molecules.json` (`convert_vermouth_martini3_small.py`) holds vermouth's 43 Martini 3
small-molecule blocks, their mappings and the CHARMM residues the mappings name atoms of. CAPS matches an all-atom
molecule to a residue by graph (elements, bonds, hydrogen counts), so any atom names do. Two things in the sources:

* vermouth's `25-small_molecule_martini3.ff` numbers the atoms of its blocks' terms from 0 while vermouth's reader
  counts from 1 (so "0" is taken for the last atom and every term is shifted by one atom); CAPS reads those blocks from
  0 (and puts MIND's virtual-site function after `--`). Read so, the 8 molecules also in cgmartini's
  `martini_v3.0.0_small_molecules_v1.itp` are identical to it (bead types, charges, bonds, constraints, angles,
  dihedrals, virtual sites, up to bead order); read vermouth's way, TOLU, FURA and MIND are not.
* the beads go at the geometric centre of their mapped atoms, hydrogens included (vermouth weighs by mass;
  `caps martini --centre mass` does that). Every molecule built all-atom by CAPS from a SMILES written
  from its residue's graph, then mapped:

| Centre | Result |
|---|---|
| geometric (default) | 43 of 43 recognised and mapped; mapped bond lengths 0.32 Å rms from the model's on average, most aromatic molecules within 0.1 Å (INDO 0.04, NAFT 0.06, PYRM 0.05); saturated two-bead rings (THP, CHXE, DIOX, THF) about 1 Å (their model bond is longer than the mapped distance with either centre) |
| mass (vermouth's) | 43 of 43; 0.80 Å rms on average |

In GROMACS (`check_martini3.py`) the 43 mapped molecules together agree to 1.1e-6, anthracene's virtual sites built on
a virtual site included (exported as `virtual_sites2` / `virtual_sites3`).


## What the validation found in the source tools

- moltemplate COMPASS: bond-increment charges doubled with `Data Bond List`; up to 55 % of class II cross-term
  reference lengths / angles on the wrong bond or angle (canonical atom sort after matching).

## Automatic typing

| Rules | Check | Result |
|---|---|---|
| `data/typing/cvff.typing.json` | Materials Studio typing in the msi2lmp examples | 811 / 823 atoms |
| `data/typing/cgenff-toppar.typing.json` | the residues of `top_all36_cgenff.rtf` (`validate_rtf_types.py --cap-links`) | 9872 / 9872 atoms, 519 / 519 residues (residues with lone pairs left out) |

CAPS follows the cff91 definition for PCFF (`c5`: aromatic carbon in a five-membered ring). Materials Studio's own
examples are not consistent with one another (ethane and the cff97 phenylalanine use generic `c`, crambin `c1`/`c2`/`c3`; crambin and nylon type NH3+ `n4`/`hn`,
the cff97 example `n+`/`h+`); CAPS follows the majority. Where the libraries lack parameters for a correctly typed
group (PCFF nitroso N=O, charged-imidazole N–H; CVFF sulfate S–O−), CAPS reports the missing terms.

GAFF rules follow antechamber's definition files (`ATOMTYPE_GFF.DEF`, `ATOMTYPE_GFF2.DEF`) in order, and its
electron-withdrawing set for h1–h5 (N O S F Cl Br I, from antechamber's ring.c). CGenFF rules follow the type
definitions in the CGenFF topology. OPLS residue-specific keys (amino-acid side-chain CB atoms, sugars, nucleic bases) need a user overlay.
