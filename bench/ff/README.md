# Force-field validation

Every converted force field in `data/forcefields` is checked against the program it came from, interaction by
interaction, and every CAPS energy form is checked against LAMMPS.

## Tools

| Script | What it checks |
|---|---|
| `check_forms.cpp` (+ `check_forms_cmp.py`) | each CAPS functional form with random coefficients against the LAMMPS style it mirrors (class II, Morse / GROMOS bonds, cosine angles, Urey–Bradley, umbrella inversion, Buckingham, Morse pair, separate 1-4 LJ) |
| `compare_moltemplate.py` | a moltemplate conversion against moltemplate's output and LAMMPS |
| `compare_class2.py` | class II (COMPASS) conversions: every cross term, bond increments, LAMMPS CLASS2 energies / forces |
| `build_library.py` | builds the moltemplate entries of `data/forcefields` with each entry's version, references, status and evidence |
| `compare_msi_types.py` | automatic typing against Materials Studio's types in msi2lmp data files |
| `check_data_lammps.py` | LAMMPS data files written by CAPS, run in LAMMPS: every energy term and force, for each force-field family, hybrid styles included |
| `validate_rtf_types.py` | automatic typing against a CHARMM topology's residues (CGenFF's own model compounds); `--apply` parameterises them |
| `convert_frc.py` | CVFF, PCFF, COMPASS and IFF 1.5 from their `.frc` files with msi2lmp's semantics; IFF's silica and metal typing rules and the silica charges IFF states |
| `check_msi2lmp.py` | the `.frc` conversions against msi2lmp's data files (LAMMPS's test structures, IFF's model database), term by term in LAMMPS; then IFF's rules on its own silica and metal models read as bare coordinates (types, charges, bonds) |
| `check_gromacs.py` | GROMACS topologies written by CAPS, run in GROMACS (grompp, a zero-step mdrun) |
| `extend_inorganic_water.py` | charges and O-H / H-O-H terms for the moltemplate water models |
| `convert_martini3.py` | Martini 3 from cgmartini.nl's martini_v300.zip: bead types and the full pair table (`data/martini/martini3-nonbonded.tsv.gz`), the molecules of its solvent, ion, small-molecule, sugar, nucleobase and phospholipid files as templates |
| `check_martini3.py` | Martini 3 boxes against GROMACS (reaction field, virtual sites, restricted bending) |
| `vermouth_ff.py` | vermouth-martinize's `.ff` and `.map` readers, for the Martini 3 small molecules |
| `convert_vermouth_martini3_small.py` | Martini 3 small molecules: vermouth's blocks and mappings, the CHARMM residues they are written for |
| `check_martini3_small.py` | every Martini 3 small molecule built all-atom, recognised by graph and mapped |
| `convert_emc_martini.py` | MARTINI overlays (polymers, solvents, surfactants, sugars) from the EMC parameter files of the moltemplate distribution |

LAMMPS with CLASS2, MOFFF, CORESHELL, MANYBODY (Stillinger-Weber) and CG-SPICA (lj/sdk, angle sdk) is built at
`~/lammps/build-class2/lmp`.

### Coarse-grained and many-body forms (`check_data_lammps.py`, DSF and PME)

| Case | Result |
|---|---|
| MARTINI 2.0: DPPC, POPE, ions, water (lj/gromacs/coul/gromacs, dielectric 15, cosine/squared angles, 1-3 and 1-4 pairs kept) | energy 1e-7, forces 8e-9 kcal/mol/Å |
| MARTINI overlays: PEO with its torsions; sucrose, maltose, glucose | energy ≤ 4e-9, forces 8e-9 |
| MARTINI amino acids: HIS, PHE, TYR, TRP (harmonic impropers, GROMACS type-2 order), ARG, ASP, LYS, water | energy 1.3e-7, forces 7.8e-9 |

| SDK: DMPC / DMPE / water; C12E8 / water; SDS / Na+ / water at relative permittivity 80 (lj/sdk 9-6 / 12-4, angle sdk with its 1-3 repulsion) | energy ≤ 1.2e-7, forces ≤ 2e-6 |
| Cooke-Deserno lipids (cosine/squared, FENE, harmonic head-tail bond) | energy 8.5e-9, forces 7.9e-9 |
| mW water (Stillinger-Weber, `.sw` file) | energy 1.8e-10, forces 8e-9 |

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

The SDK beads for PEG and SDS come from LAMMPS's own SDK / SPICA examples (`examples/PACKAGES/cgspica`); their CM, CT,
CT2 and W terms are checked equal to the library's before they are added (`bench/typing/make_cg_rules.py`).

## What the validation found in the source tools

- moltemplate COMPASS: bond-increment charges doubled with `Data Bond List`; up to 55 % of class II cross-term
  reference lengths / angles on the wrong bond or angle (canonical atom sort after matching).

## Automatic typing

| Rules | Check | Result |
|---|---|---|
| `data/typing/cvff.typing.json` | Materials Studio typing in the msi2lmp examples | 811 / 823 atoms |
| `data/typing/cgenff-toppar.typing.json` | the residues of `top_all36_cgenff.rtf` (`validate_rtf_types.py --cap-links`) | 9872 / 9872 atoms, 519 / 519 residues (residues with lone pairs left out) |
| `data/typing/sdk-moltemplate.typing.json` (bead mapping) | all-atom DMPC, C12E8, SDS and n-alkanes C12–C30 | the model's own bead topologies (Shinoda's 13-bead DMPC; the LAMMPS examples' C12E8 and SDS); every chain length tiles |
| `data/typing/martini-*.typing.json` | the 94 bead templates of the MARTINI source files | all built, typed and parameterised (ILE, LEU, PRO, VAL not templated: their AC1 / AC2 beads have no parameters in the sources) |

CAPS follows the cff91 definition for PCFF (`c5`: aromatic carbon in a five-membered ring). Materials Studio's own
examples are not consistent with one another (ethane and the cff97 phenylalanine use generic `c`, crambin `c1`/`c2`/`c3`; crambin and nylon type NH3+ `n4`/`hn`,
the cff97 example `n+`/`h+`); CAPS follows the majority. Where the libraries lack parameters for a correctly typed
group (PCFF nitroso N=O, charged-imidazole N–H; CVFF sulfate S–O−), CAPS reports the missing terms.

GAFF rules follow antechamber's definition files (`ATOMTYPE_GFF.DEF`, `ATOMTYPE_GFF2.DEF`) in order, and its
electron-withdrawing set for h1–h5 (N O S F Cl Br I, from antechamber's ring.c). CGenFF rules follow the type
definitions in the CGenFF topology. OPLS residue-specific keys (amino-acid side-chain CB atoms, sugars, nucleic bases) need a user overlay.
