# Force-field validation

Every converted force field in `data/forcefields` is checked against the program it came from, interaction by
interaction, and every CAPS energy form is checked against LAMMPS.

## Tools

| Script | What it checks |
|---|---|
| `check_forms.cpp` (+ `check_forms_cmp.py`) | each CAPS functional form with random coefficients against the LAMMPS style it mirrors (class II, Morse / GROMOS bonds, cosine angles, Urey–Bradley, umbrella inversion, Buckingham, Morse pair, separate 1-4 LJ) |
| `compare_moltemplate.py` | a moltemplate conversion against moltemplate's output and LAMMPS |
| `compare_class2.py` | class II (COMPASS) conversions: every cross term, bond increments, LAMMPS CLASS2 energies / forces |
| `run_dlfield.py` | runs a private copy of DL_FIELD (`build/dlfield_work`) on a structure; `--dlpoly` keeps the DL_POLY FIELD, `--keys` shows library keys |
| `compare_dlfield.py` | a DL_FIELD conversion against DL_FIELD's LAMMPS output, with LAMMPS energies and forces |
| `compare_field.py` | a DL_FIELD conversion against DL_FIELD's DL_POLY FIELD file, term by term (`--from-config` for proteins, `--from-field` for ionic systems) |
| `validate_family.py` | a whole library on its own molecule templates (`template_pdb.py`) or ionic formula units (`ionic_pdb.py`) |
| `build_library.py` | builds `data/forcefields` with each entry's version, references, status and evidence |
| `validate_typing.py` | automatic typing of a library's own molecule templates (mol2 without bond orders) against the template types; `--apply` also parameterises each one |
| `compare_msi_types.py` | automatic typing against Materials Studio's types in msi2lmp data files |
| `check_data_lammps.py` | LAMMPS data files written by CAPS, run in LAMMPS: every energy term and force, for each force-field family, hybrid styles included |
| `validate_rtf_types.py` | automatic typing against a CHARMM topology's residues (CGenFF's own model compounds); `--apply` parameterises them |
| `convert_frc.py` | CVFF, PCFF, COMPASS and IFF 1.5 from their `.frc` files with msi2lmp's semantics; IFF's silica and metal typing rules and the silica charges IFF states |
| `check_msi2lmp.py` | the `.frc` conversions against msi2lmp's data files (LAMMPS's test structures, IFF's model database), term by term in LAMMPS; then IFF's rules on its own silica and metal models read as bare coordinates (types, charges, bonds) |
| `check_gromacs.py` | GROMACS topologies written by CAPS, run in GROMACS (grompp, a zero-step mdrun) |
| `extend_inorganic_water.py` | charges, shells and fixes for the inorganic and water force fields from their sources' templates |
| `convert_martini3.py` | Martini 3 from cgmartini.nl's martini_v300.zip: bead types and the full pair table (`data/martini/martini3-nonbonded.tsv.gz`), the molecules of its solvent, ion, small-molecule, sugar, nucleobase and phospholipid files as templates |
| `check_martini3.py` | Martini 3 boxes against GROMACS (reaction field, virtual sites, restricted bending) |
| `convert_vermouth_martini3.py` | Martini 3 proteins from vermouth-martinize (martini3001 force field and mappings) |
| `check_martini3_protein.py` | Martini 3 proteins against martinize2, term by term |
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
| Martini 2.2 proteins: AK helix (helix constraints, 96° / 700 angles, −120° / 400 dihedrals); an aromatic helix (W, Y, F, H rings, charged termini); 1ICO (VAL / ILE as AC1 / AC2, elastic bonds, disulfide) | energy ≤ 1.4e-7, forces ≤ 2.5e-8 |

| SDK: DMPC / DMPE / water; C12E8 / water; SDS / Na+ / water at relative permittivity 80 (lj/sdk 9-6 / 12-4, angle sdk with its 1-3 repulsion) | energy ≤ 1.2e-7, forces ≤ 2e-6 |
| Cooke-Deserno lipids (cosine/squared, FENE, harmonic head-tail bond) | energy 8.5e-9, forces 7.9e-9 |
| mW water (Stillinger-Weber, `.sw` file) | energy 1.8e-10, forces 8e-9 |

`check_martini_protein.py`: CAPS's Martini 2.2 protein of 1ICO (vermouth's martini22 integration test, `-dssp -cys auto`)
against martinize2 — 67 beads with identical types and charges, 67 bonds, 21 constraints, 65 angles, 4 dihedrals
identical, bead positions within 0.001 Å; CAPS's DSSP gives DSSP 2.0's letters on 1ICO, the Trp-cage (α and 3₁₀) and a
43-residue helical protein.

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

| Ubiquitin 1UBQ, all-atom mapped to Martini 3 (restricted-bending backbone, side-chain fix dihedrals) | energy 2.3e-6 |
| Lysozyme 3LZT, all-atom mapped to Martini 3 (tryptophan virtual sites, exclusions, impropers) | energy 6.1e-6 |
| Lysozyme as martinize2 wrote it: its cg.pdb and topol.top read by CAPS (`--topology`), elastic network 500 | energy 1.7e-5 |
| All 43 Martini 3 small molecules, all-atom (built by CAPS) mapped by graph | energy 1.1e-6 |

Forces on beads held by constraints (stiff bonds in CAPS, 1e6 kJ/mol/nm²) carry the single-precision position error of
GROMACS (~0.012 kcal/mol/Å per stiff bond in a 9 nm box); a mapped protein that is not relaxed has beads under
1e5 kcal/mol/Å, which single precision carries to ~4e-5 of the force. The check allows those and nothing else. The
.gro holds the virtual sites where CAPS places them (GROMACS's rerun takes coordinates as written).

`check_martini3_protein.py`: CAPS's Martini 3 proteins against martinize2 (-ff martini3001), vermouth's Martini 3
integration tests with their own commands — beads, types and charges, every bond (the rubber band included),
constraint, angle, dihedral, improper, virtual site and exclusion, bead positions:

| Test | Result |
|---|---|
| PRO-PRO, neutral termini, no secondary structure (vermouth 0.7.4) | identical, positions 0.001 Å |
| Histatin 5, a disordered region, hydrogens in the input (0.15.1) | identical, positions 0.001 Å |
| Ubiquitin 1UBQ, DSSP, side-chain fix, elastic network (0.0.1) | identical, DSSP identical, positions 0.001 Å |
| Lysozyme 3LZT, DSSP, elastic network 500, four disulfides, six tryptophans (0.15.1) | identical, DSSP identical, positions 0.001 Å |
| Insulin, two chains joined by disulfides, network per chain (0.7.3) | identical but for 2 elastic bonds (chain A 6–9 and 7–10) that vermouth 0.7.3 left out; DSSP differs in 5 letters (that DSSP ranks π-helices before α, DSSP 2.0 the other way; both are helix to Martini) |

The model is data (`data/martini/martini3-protein.json`, from vermouth's martini3001 files by
`convert_vermouth_martini3.py`): blocks, modifications, the 66 links with vermouth's matching rules, the mappings.

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
  `caps martini --martini 3 --centre mass` does that). Every molecule built all-atom by CAPS from a SMILES written
  from its residue's graph, then mapped:

| Centre | Result |
|---|---|
| geometric (default) | 43 of 43 recognised and mapped; mapped bond lengths 0.32 Å rms from the model's on average, most aromatic molecules within 0.1 Å (INDO 0.04, NAFT 0.06, PYRM 0.05); saturated two-bead rings (THP, CHXE, DIOX, THF) about 1 Å (their model bond is longer than the mapped distance with either centre) |
| mass (vermouth's) | 43 of 43; 0.80 Å rms on average |

### Martini 3 lipids

The 109 phospholipids of cgmartini's `martini_v3.0.0_phospholipids_v1.itp` (PC, PE, PS, PG, PA) carry insane's building
blocks (head, linker, one C / D letter per tail bead); the templates keep them. No Martini 3 atomistic lipid mapping is
published with the models, so CAPS maps by the building blocks (its own rule, `martini3_lipids`): the phosphate, the
glycerol with its two esters, the head group and each acyl chain are found by structure; a chain is split as evenly as
possible into the model's tail beads, a bead "D" when a double bond starts in it (this reproduces the models' letters
for oleoyl, linoleoyl, linolenoyl, arachidonoyl, EPA); of the templates whose letters match, the one nearest four
carbons per tail bead is the lipid. POPC, DPPC, DOPE, POPS, POPG, POPA and DAPC built from SMILES are recognised as
themselves; mapped bond lengths about 0.8 Å rms from the model's before relaxing.

In GROMACS (`check_martini3.py`) the 43 mapped molecules together agree to 1.1e-6, anthracene's virtual sites built on
a virtual site included (exported as `virtual_sites2` / `virtual_sites3`).

The SDK beads for PEG and SDS come from LAMMPS's own SDK / SPICA examples (`examples/PACKAGES/cgspica`); their CM, CT,
CT2 and W terms are checked equal to the library's before they are added (`bench/typing/make_cg_rules.py`).

## Results (DL_FIELD 4.13, FIELD-file comparison unless noted)

| Library | Result |
|---|---|
| PCFF, CVFF, COMPASS, OPLS 2005, GAFF 16 / 25 | also against DL_FIELD's LAMMPS output: energies and forces ≤ 6e-5 kcal/mol/Å |
| PCFF | 57 / 59 templates (open: one auto-torsion precedence case, two template-listed inversions) |
| CVFF | 20 / 20 |
| COMPASS | 7 / 7 |
| OPLS 2005 | 59 / 60 (aniline: DL_FIELD's template fixes an improper's atom order) |
| CL&P | 24 / 24 |
| AMBER | BPTI (6PTI), methyl glucoside |
| CHARMM | 48 / 48 |
| CHARMM22 proteins | 35 / 35 templates, SOD1 protein |
| CHARMM36 proteins | 4 / 4 templates, BPTI |
| CHARMM36 nucleic acids | DL_FIELD's nucleic-acid example |
| CGenFF | 60 / 60 |
| CHARMM19 | 3 / 3 templates, C18 alkane |
| DREIDING | 29 / 29 |
| GROMOS 54A7 | octane (the other templates are chain residues) |
| TraPPE-UA | 23 / 25 (the two are DL_FIELD's C-C-O-H torsion error) |
| TraPPE-EH | 9 / 9 |
| Inorganic halides, GaN, binary oxides, ternary oxides, zeolites | 14 / 14, 1 / 1, 17 / 17 (+ Al2O3, MgO shell models), 12 / 12, 4 / 4 |

Templates DL_FIELD itself cannot build (chain residues, some nucleotides) are skipped. Two kinds of difference are
reported but not counted as failures: *template policy* (DL_FIELD places impropers only where its molecule templates
list them; CAPS applies the rules at every planar centre they match) and *convention* (same rule and atoms, another
order of symmetric outer atoms; DL_FIELD's harmonic approximation of GROMOS bonds).

## What the validation found in the source tools

- moltemplate COMPASS: bond-increment charges doubled with `Data Bond List`; up to 55 % of class II cross-term
  reference lengths / angles on the wrong bond or angle (canonical atom sort after matching).
- DL_FIELD PCFF / COMPASS / CVFF: no class II cross terms (DL_POLY cannot evaluate them).
- DL_FIELD LAMMPS export: CHARMM 1-4 van der Waals dropped (`special_bonds lj 0 0 0` with dihedral weight 0),
  DREIDING inversions dropped, GROMOS not written; its DL_POLY FIELD file is complete.
- DL_FIELD GROMOS: quartic bonds written as a harmonic approximation.
- DL_FIELD TraPPE-UA: the C-C-O-H torsion takes terms from the O-C-C-O entry; the library's `CH2 CH2 OA HA` torsion is
  continued by rows keyed `CH3 CH2 OA HA`.
- DL_FIELD FIELD files: large fixed-width values run together (`7425.64200011530.500000`).
- DL_FIELD libraries: `buckinghsm` (read as Buckingham, as DL_FIELD reads four characters); wrong masses for
  `HC_benzyl` (12.0115, PCFF), `C_benzene` (1.00797, CVFF), `O_thioester` (15.0994), and in GAFF 2025
  `C_carboxylate`, `O_carboxylate`, `HCN_azetidine`; CAPS keeps the first, correct entry and notes the others.
- DL_FIELD GAFF 2025: one carbonyl improper written with the key's end atom moved; CAPS follows the key.
- DL_FIELD units: eV → kcal/mol 23.061, kJ → kcal 0.23901, K → kcal 0.0019872041 (rounded); CAPS uses the exact
  constants, so values agree to ~2e-5 relative.

## Automatic typing

| Rules | Check | Result |
|---|---|---|
| `data/typing/pcff.typing.json` | 92 DL_FIELD PCFF templates (chain-residue fragments left out) | 1387 / 1400 atoms, 87 molecules fully right |
| `data/typing/cvff.typing.json` | 23 DL_FIELD CVFF templates | 272 / 278 atoms (the rest: TIP3P / SPC water) |
| `data/typing/cvff.typing.json` | Materials Studio typing in the msi2lmp examples | 811 / 823 atoms |
| `data/typing/cgenff.typing.json` | the 494 model compounds of `top_all36_cgenff.rtf` (CGenFF 3.0.1) | 9024 / 9174 atoms, 417 residues fully right (97.7 % with the RTF's partial bond orders) |
| `data/typing/gaff-amber16.typing.json` | 32 DL_FIELD GAFF templates (water models aside) | 365 / 373 atoms |
| `data/typing/gaff-amber25.typing.json` | 27 DL_FIELD GAFF2 templates | 304 / 316 atoms |
| `data/typing/opls2005.typing.json` | 99 DL_FIELD OPLS 2005 templates | 1368 / 1398 atoms, 91 molecules fully right |
| `data/typing/sdk-moltemplate.typing.json` (bead mapping) | all-atom DMPC, C12E8, SDS and n-alkanes C12–C30 | the model's own bead topologies (Shinoda's 13-bead DMPC; the LAMMPS examples' C12E8 and SDS); every chain length tiles |
| `data/typing/martini-*.typing.json` | the 94 bead templates of the MARTINI source files | all built, typed and parameterised (ILE, LEU, PRO, VAL not templated: their AC1 / AC2 beads have no parameters in the sources) |

DL_FIELD's templates are not consistent with one another: its PCFF templates give the carbons of furan, oxazole and
indole `cp` but those of pyrrole, isoxazole and benzoxazole `c5`; neutral histidine carbons `ci` (the charged-ring type).
CAPS follows the cff91 definition (`c5`: aromatic carbon in a five-membered ring). Materials Studio's own examples differ
too (ethane and the cff97 phenylalanine use generic `c`, crambin `c1`/`c2`/`c3`; crambin and nylon type NH3+ `n4`/`hn`,
the cff97 example `n+`/`h+`); CAPS follows the majority. Where the libraries lack parameters for a correctly typed
group (PCFF nitroso N=O, charged-imidazole N–H; CVFF sulfate S–O−), CAPS reports the missing terms, as DL_FIELD does.

GAFF rules follow antechamber's definition files (`ATOMTYPE_GFF.DEF`, `ATOMTYPE_GFF2.DEF`) in order, and its
electron-withdrawing set for h1–h5 (N O S F Cl Br I, from antechamber's ring.c). DL_FIELD's GAFF templates depart from
antechamber for H on carbons next to O or S (hc where antechamber gives h1) and for azobenzene N (n2, where antechamber's
conjugation rule gives ne / nf); CAPS follows antechamber. CGenFF rules follow the type definitions in the CGenFF
topology; DL_FIELD's CGenFF templates often use CHARMM's generic HA for alkane H (CGenFF: HGA1–3), which CAPS does not
copy. OPLS residue-specific keys (amino-acid side-chain CB atoms, sugars, nucleic bases) need a user overlay.
