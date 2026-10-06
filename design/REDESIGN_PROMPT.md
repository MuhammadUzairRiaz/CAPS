# Design prompt v2: CAPS (Chain Assembly and Packing Suite)

> Paste everything below the line into your AI assistant or give it to your engineering team.

---

## Role

You are a principal architect and computational physicist. You have built
molecular dynamics engines, Monte Carlo polymer builders, and force-field
software used in published research. You know the internals of LAMMPS,
GROMACS, packmol, and Materials Studio (Amorphous Cell, Forcite). You write
production C++, modern Fortran, CUDA or SYCL, Python, and C#/.NET. You are
strict about scientific correctness: every equation you implement comes from
peer-reviewed literature, and you cite it.

## The problem

Building a simulation-ready polymer model today means chaining several
programs: a monomer and chain builder, packmol for packing, LAMMPS or GROMACS
for minimisation and relaxation, and DL_FIELD or moltemplate for atom typing
and parameters. Each step has its own formats and failure modes, growth often
uses a generic soft-sphere potential rather than the real force field, and
nothing records how a model was made.

## Vision

Design and build **CAPS (Chain Assembly and Packing
Suite)**. CAPS must be **self-contained**. It needs no external executable
for any step: building, typing, packing, minimisation, molecular dynamics,
equilibration, crosslinking, and analysis all run inside CAPS on its own
physics engines.

CAPS still **writes** input files for LAMMPS, GROMACS, and other engines, so
users can run long production simulations wherever they like. External
engines are used in development only to cross-validate CAPS, never at runtime.

Libraries linked at build time are allowed when their licence permits, for
example RDKit (BSD), FFTW or pocketFFT, Eigen, and oneTBB. Launching another
program as a subprocess is not.

The engines must be original implementations of published methods. They
should be good enough to publish in a peer-reviewed journal, with reproducible
benchmarks.

## Scientific integrity rules

These rules apply to every line of physics code.

1. **Cite the source.** Every equation, algorithm, and default parameter
   carries a literature reference in the code and the docs. It gives the
   paper, equation number, and any deviation from the paper.
2. **Invent nothing.** Force-field parameters come only from the published
   force field. Missing parameters are reported as errors, never guessed
   silently. Any estimated parameter is flagged in the output and the log.
3. **Use one unit system internally** and document it. Convert only at I/O
   boundaries, with unit tests for every conversion.
4. **Check energy and force parity before speed.** A kernel counts as correct
   only when its energies and forces match LAMMPS and GROMACS on the same
   system to a stated tolerance. Speed work starts only after that.
5. **Document approximations.** Cut-offs, tail corrections, and mixed
   precision are written into the provenance manifest of every output.
6. **Report results honestly.** Benchmarks report where CAPS is slower or
   less accurate, not only where it wins.

## The CAPS engines

Design each engine as an independent library with a clean C++ API and a
stable C ABI. For each one, specify its physics, algorithms, data structures,
parallel strategy, and validation tests.

### 1. CAPS Chem: chemistry and topology core
- SMILES and SMARTS parsing, ring perception, aromaticity, bond orders, and
  stereo. Use RDKit linked as a C++ library, or write a native implementation
  if you justify it.
- 3D embedding: distance geometry with ETKDG (Riniker and Landrum, J. Chem.
  Inf. Model. 2015).
- Polymer graph model: repeat units, head and tail atoms, tacticity,
  sequences (alternating, block, random, gradient), branches, stars, and
  networks.
- Charges come from the force field. Fallbacks are Gasteiger–Marsili
  (Tetrahedron 1980) and QEq (Rappé and Goddard, J. Phys. Chem. 1991), with
  neutrality enforced per molecule.

### 2. CAPS Field: force-field engine
- Class I force fields: OPLS-AA/UA (Jorgensen et al., J. Am. Chem. Soc.
  1996), GAFF and GAFF2, DREIDING, TraPPE, CHARMM, AMBER.
- Class II force fields: PCFF and COMPASS (Sun, J. Phys. Chem. B 1998),
  including the cross terms.
- Coarse-grained models: Kremer–Grest (J. Chem. Phys. 1990), MARTINI, SDK.
- Automatic atom typing with rule files of SMARTS patterns and chemical
  environments, and a report explaining why each atom got its type.
- A declarative, versioned parameter database (TOML or JSON schema) with a
  literature source for each parameter set. Import from moltemplate .lt, DL_FIELD
  .par, and GROMACS .itp files. Export to all of them.

### 3. CAPS Grow: amorphous chain growth
- Keep and extend the existing method: Theodorou–Suter growth (Macromolecules
  1985), RIS statistics (Flory 1969, Mattice and Suter 1994), Rosenbluth
  weighting (J. Chem. Phys. 1955), and Meirovitch scanning.
- Add configurational-bias growth (Siepmann and Frenkel, Mol. Phys. 1992)
  that uses the **real force-field energy** during growth, not only a soft
  sphere.
- Add all-atom growth that samples side-group torsions.
- Support ring backbones (PET, PC, PBAT) through fragment-based growth.
- Support polydispersity (Schulz–Zimm, log-normal, user histogram),
  branching, stars, solvents, dissolved gases, and oriented or nematic
  cells.
- Grow chains in parallel with spatial domain locking.

### 4. CAPS Pack: native packing, faster than packmol
- Baseline: packmol's method, which minimises an overlap penalty with GENCAN
  (Martínez et al., J. Comput. Chem. 2009). Re-derive the objective and
  document it.
- Rigid-body formulation: centre of mass plus a quaternion for each molecule,
  with analytic gradients.
- Optimisers: L-BFGS (Liu and Nocedal, Math. Program. 1989) and FIRE (Bitzek
  et al., Phys. Rev. Lett. 2006).
- Linked-cell and Verlet neighbour lists with Morton-ordered memory for
  cache efficiency. Multithreaded CPU kernels with SIMD, plus GPU kernels.
- Staged strategy: random sequential insertion at low density, overlap
  minimisation, then soft-core compression MD to the target density.
- Support region constraints for layers, slabs, spheres, cylinders,
  interfaces, and fixed molecules.
- Target: at least 10 times faster than packmol on 100,000 atoms or more,
  with an equal or better minimum-distance guarantee. Prove it with the
  benchmark tables below. Do not claim a speed-up before measuring it.

### 5. CAPS Relax: energy minimisation
- Steepest descent, Polak–Ribière conjugate gradient, L-BFGS, and FIRE.
- Soft-core and force-capped push-off for badly overlapping structures
  (Auhl et al., J. Chem. Phys. 2003).
- Minimisation of the box shape and volume under constant pressure.

### 6. CAPS Dynamics: molecular dynamics engine
- Integrators: velocity Verlet (Swope et al., J. Chem. Phys. 1982) and
  multiple time stepping with r-RESPA (Tuckerman, Berne and Martyna, J. Chem.
  Phys. 1992).
- Thermostats: Nosé–Hoover chains (Martyna, Klein and Tuckerman, J. Chem.
  Phys. 1992), the stochastic velocity-rescaling thermostat (Bussi, Donadio
  and Parrinello, J. Chem. Phys. 2007), and Langevin dynamics with the
  BAOAB splitting (Leimkuhler and Matthews, 2013).
- Barostats: Berendsen for early relaxation only (J. Chem. Phys. 1984), MTK
  (Martyna, Tobias and Klein, J. Chem. Phys. 1994), and stochastic cell
  rescaling (Bernetti and Bussi, J. Chem. Phys. 2020). Support isotropic,
  anisotropic, and triclinic boxes.
- Constraints: SHAKE (Ryckaert et al., J. Comput. Phys. 1977), RATTLE
  (Andersen, J. Comput. Phys. 1983), and LINCS (Hess et al., J. Comput. Chem.
  1997).
- Electrostatics: Ewald summation, smooth particle-mesh Ewald (Essmann et
  al., J. Chem. Phys. 1995), and damped shifted force or Wolf summation (Wolf
  et al., J. Chem. Phys. 1999) as a fast option.
- Van der Waals: cut-offs with long-range tail corrections, switching
  functions, and all standard mixing rules.
- Precision: double precision by default. Mixed-precision GPU kernels are
  allowed, with a deterministic mode for exact reproducibility.
- Target: correct NVE energy conservation, and speed competitive with LAMMPS
  for systems up to about one million atoms. Beating GROMACS at production
  MD is a later goal, not a first-release promise.

### 7. CAPS Equilibrate: published equilibration protocols
- The 21-step compression and decompression protocol (Larsen et al.,
  Macromolecules 2011).
- Push-off for long-chain melts (Auhl et al., J. Chem. Phys. 2003).
- Double-bridging and end-bridging Monte Carlo moves for long chains
  (Karayiannis, Mavrantzas and Theodorou, Phys. Rev. Lett. 2002).
- Simulated annealing cycles.
- Convergence criteria for density, energy, Rg, and internal-distance plots
  that decide automatically when a cell is equilibrated.

### 8. CAPS React: crosslinking and reactions
- Keep the existing atom-mapped reaction templates.
- Form bonds by distance and probability during MD, as in
  REACTER (Gissinger, Jensen and Wise, Polymer 2017), with relaxation
  after each step. Also implement the Polymatic cycle (Abbott, Hart and
  Colina, Theor. Chem. Acc. 2013).
- Retype atoms automatically after each reaction.
- Targets for degree of conversion, plus gel-point detection from cluster
  analysis and comparison with Flory–Stockmayer theory.

### 9. CAPS Analyze: property calculation
- Structure: density, radial distribution functions, structure factor,
  and X-ray and neutron scattering patterns.
- Chains: Rg, end-to-end distance, characteristic ratio C_n and C∞,
  persistence length, and entanglement analysis.
- Thermodynamics: cohesive energy density, Hildebrand solubility parameter,
  and Tg from density versus temperature cooling curves.
- Mechanics: elastic constants by the constant-strain method and by strain
  fluctuations, and uniaxial stress–strain curves.
- Dynamics: MSD, diffusion coefficients, and relaxation times.
- Free volume and pore size by probe insertion.

### 10. CAPS IO: file writers
- Write LAMMPS data and input files, GROMACS gro, top, itp and mdp files,
  and PDB, mol2, xyz, extended XYZ, CIF and DCD files. Keep the moltemplate
  and DL_FIELD formats.
- Input decks written for other engines must reproduce CAPS energies on the
  first step. This is tested automatically.

### 11. CAPS Studio: interactive 3D builder and visualiser

CAPS Studio is the main workspace of the desktop app. It must match or exceed
the builder and visualiser in Materials Studio and other modern tools such as
Avogadro 2, GaussView, and Maestro. Design an original interface, icons, and
naming. Do not copy another product's proprietary UI, artwork, or scripting
API.

**Architecture.** All editing logic lives in the C++ core as undoable
commands, not in the GUI. Every action a user takes with the mouse is also
available as a command in the CLI and Python API, so any structure built by
hand can be replayed, scripted, and tested. The renderer is a C++ library,
CAPS View, embedded in the Avalonia window.

**Building and sketching**
- Place atoms by clicking in 3D space. Pick elements from a full periodic
  table, with isotopes, formal charges, and valence states.
- Draw bonds by dragging between atoms. Cycle bond order through single,
  double, triple, aromatic, and partial bonds, and support coordinate bonds.
- Add hydrogens automatically with correct VSEPR geometry, and re-adjust them
  as the structure changes. Allow removing or hiding hydrogens.
- Change hybridisation and geometry: linear, trigonal, tetrahedral, square
  planar, trigonal bipyramidal, and octahedral.
- A fragment library: rings, functional groups, amino acids, nucleotides,
  sugars, ligands, solvents, ions, and common monomers. Users can save their
  own fragments.
- Attach fragments to open valences and fuse rings onto existing rings.
- Create structures from SMILES, InChI, or a chemical name, and from a 2D
  sketch panel that converts to 3D.
- Perceive bonds from coordinates using covalent radii (Cordero et al.,
  Dalton Trans. 2008) for files that carry no bonds.
- Keep several molecules in one document. Merge, split, copy, and paste
  structures between documents.

**Editing**
- Delete atoms, bonds, fragments, or selections. Replace an element and keep
  the geometry.
- Make and break bonds, and join two fragments at chosen atoms.
- Move, rotate, and translate atoms, fragments, or whole molecules, freely or
  along an axis.
- Set an exact bond length, angle, or dihedral, choosing which side of the
  molecule moves.
- Mirror, invert, align, and centre structures. Align one molecule onto
  another by RMSD superposition (Kabsch, Acta Cryst. A 1976).
- Unlimited undo and redo, with an editable history list.

**Measuring and monitoring**
- Measure distances, angles, dihedrals, and out-of-plane angles by clicking
  atoms.
- Pin measurements as live monitors that update while editing, minimising,
  or playing a trajectory.
- Detect and display hydrogen bonds, close contacts, and clashes.
- Report centre of mass, moments of inertia, molecular weight, formula,
  net charge, dipole moment from charges, and molecular volume and surface
  area.

**Stereochemistry**
- Assign R/S labels by Cahn–Ingold–Prelog rules, using the algorithm and
  validation set of Hanson et al. (J. Chem. Inf. Model. 2018). Assign E/Z
  labels to double bonds.
- Show labels in 3D, invert a stereocentre with one action, and set a
  wanted R/S or E/Z configuration directly.
- Detect axial and planar chirality, and warn about undefined stereocentres.
- For polymers, set tacticity (isotactic, syndiotactic, atactic, or a custom
  sequence) and check it along the chain.

**Selection**
- Select by click, box, or lasso, and by element, atom type, charge,
  residue, molecule, or fragment.
- Select by SMARTS pattern, by distance from other atoms, or by growing
  the selection along bonds.
- Invert selections and save named selection sets.

**Cleaning and optimisation**
- Clean up geometry instantly with a force field, as in Avogadro's
  auto-optimise mode, running in the background while the user edits.
- Minimise with constraints: fixed atoms, fixed distances, angles, and
  dihedrals, and restraints.
- Conformer search with ETKDG and minimisation, and torsion scans that
  plot energy against angle.
- All of this uses CAPS Field and CAPS Relax, so builder geometry matches
  the simulation engine exactly.

**Specialised builders**
- **Polymer builder:** define the repeat unit and its head and tail atoms.
  Build homopolymers and random, alternating, block, and gradient
  copolymers. Set tacticity, head-to-tail or head-to-head linkage, backbone
  torsions, and end caps. Build branched, star, comb, dendrimer, and network
  polymers. Send any chain straight to CAPS Grow.
- **Crystal builder:** all 230 space groups from the International Tables
  for Crystallography, Volume A. Enter lattice parameters and fractional
  coordinates, apply symmetry, find symmetry in an existing structure (for
  example with spglib, Togo et al.), and make supercells and primitive cells.
- **Surface and interface builder:** cleave surfaces by Miller indices, add
  vacuum slabs, and stack layers or interfaces with lattice matching.
- **Nanostructure builder:** graphene and other 2D sheets, nanotubes with a
  chosen chirality, nanoparticles and clusters, and pores.
- **Biomolecule builder:** peptides from a sequence with chosen secondary
  structure, and nucleic acids from a sequence.
- **Solvation and ions:** fill a box with solvent, and add ions to a chosen
  concentration or to neutralise charge. This uses CAPS Pack.

**Visualisation**
- Styles: ball and stick, stick, wireframe, space filling, polyhedra, and
  ribbons for biomolecules. Mix styles within one structure.
- Colour by element, atom type, charge, molecule, chain, residue, or any
  computed property.
- Labels for element, index, atom type, charge, R/S, and custom text.
- Show the unit cell, periodic images, and wrapped or unwrapped molecules.
- Surfaces: van der Waals, solvent-accessible, and Connolly surfaces
  (Connolly, J. Appl. Cryst. 1983), built with marching cubes (Lorensen and
  Cline, 1987). Map electrostatic potential or other properties onto them.
- Rendering quality: ray-cast sphere and cylinder impostors, ambient
  occlusion and edge outlines for depth (as in QuteMol, Tarini et al., IEEE
  TVCG 2006), depth cueing, orthographic and perspective views, and clip
  planes.
- Export high-resolution images, and movies of trajectories and rotations.
- Play MD trajectories from CAPS Dynamics with live plots linked to the
  3D view.

**Force-field view**
- Show atom types and charges in 3D. Explain why each atom got its type.
- Override atom types and charges by hand, with the change recorded in
  provenance.
- Highlight untyped atoms and missing parameters before a run starts.

**Validation while building**
- Warn in real time about wrong valences, overlapping atoms, unusual bond
  lengths or angles, unbalanced charge, and broken rings.
- Offer one-click fixes where the correct answer is clear.

**Productivity**
- Multiple documents in tabs, split views, and a project explorer.
- Drag and drop files to import them.
- Record any sequence of actions as a Python script and replay it. Optionally
  provide a Perl binding for users coming from Materials Studio's Perl
  scripting.
- Customisable keyboard shortcuts and mouse modes. Support trackpad
  gestures and 3D mice.
- Light and dark themes, scalable UI, colour-blind-safe palettes, and full
  keyboard access.

**File formats**
- Import and export: PDB, mmCIF, CIF, mol2, SDF and MOL, XYZ and extended
  XYZ, CML, SMILES, InChI, LAMMPS data, GROMACS gro and top, and trajectory
  formats (DCD, XTC, TRR, LAMMPS dump).

**Performance targets**
- Smooth 60 frames per second for 1,000,000 atoms on a mid-range GPU,
  using instancing and level of detail.
- Every editing action responds in under 16 ms for molecules up to 10,000
  atoms.
- Structures of 1,000,000 atoms load in under 5 seconds.

## Language roles

Use each language only where it clearly earns its place. Justify any
deviation.

| Layer | Technology | Responsibility |
|---|---|---|
| Engine core | C++20, CMake | All engines, data structures, neighbour lists, typing, growth, packing, MD loop, I/O |
| Numeric kernels | Modern Fortran 2018 with iso_c_binding | PME reciprocal space, Ewald sums, RIS transfer matrices, dense linear algebra. Called through a thin C ABI |
| Acceleration | OpenMP and SIMD on CPU. CUDA and HIP, or SYCL, on GPU | Nonbonded forces, PME charge spreading, packing penalty, growth trials |
| Scripting and API | Python 3.12+ via nanobind or pybind11, on PyPI and conda-forge | Jupyter API, YAML/JSON workflows, benchmark and figure scripts |
| Desktop GUI | C# / .NET 9 with Avalonia UI (MVVM) | Application shell, panels, dialogs, specialised builder forms, live energy and density plots, job control |
| 3D renderer | C++ with Vulkan, using MoltenVK on macOS, embedded in Avalonia | CAPS View: impostor rendering, surfaces, picking, trajectories, image and movie export |
| GUI-to-engine bridge | P/Invoke over a versioned C ABI, plus gRPC | Local runs, and remote runs on HPC clusters through a headless daemon |
| CLI | Single static C++ binary | All commands, scriptable for batch runs on clusters |
| Legacy text tooling | Python by default. Perl only to wrap existing Perl scripts | Converting force-field parameter files |
| Docs | Sphinx with Doxygen and Breathe, and LaTeX for the theory manual | API reference, theory with every equation and citation, tutorials |

## Architecture requirements

1. **Hexagonal layout.** A pure physics core with no I/O or UI. Adapters for
   formats, GUI, CLI, and Python sit around it.
2. **Plugin interfaces** for force fields, typers, growth energy models,
   integrators, thermostats, barostats, file formats, and reaction rules.
3. **Deterministic reproducibility.** The same seed, inputs, and version give
   bit-identical output on every OS in deterministic mode. Use a
   counter-based random number generator such as Philox.
4. **Provenance.** Every output carries a manifest with input hashes,
   version, parameters, approximations, citations, and hardware.
5. **Job model.** Every long operation is an async, cancellable job with
   progress events, structured logs, and resumable checkpoints.
6. **Performance by design.** Use structure-of-arrays layouts, spatial
   sorting, and neighbour lists with a Verlet buffer. Profile before
   optimising, and track performance in CI to catch regressions.

## Validation and benchmark suite

Build this suite as a first-class part of the product, not an afterthought.
It must run with one command and write every table below as Markdown, CSV,
and LaTeX, ready for a paper. Use the same hardware and pinned versions for
every tool. Repeat each run at least five times and report the mean and
standard deviation. Record the CPU, GPU, core count, compiler, and every tool
version. Publish the scripts and raw data with the paper.

Before publishing any comparison with Materials Studio, check its licence
agreement. Many commercial licences restrict publishing benchmarks.

**Table 1. Energy and force parity with reference engines.**

| System | Atoms | Force field | Term | CAPS (kcal/mol) | LAMMPS | GROMACS | Max relative error | Force RMSD |
|---|---|---|---|---|---|---|---|---|

**Table 2. Packing performance against packmol.**

| System | Molecules | Atoms | Target density | Tool | Threads/GPU | Wall time (s) | Minimum distance (Å) | Success rate | Speed-up |
|---|---|---|---|---|---|---|---|---|---|

**Table 3. Amorphous cell building against other builders.**

Compare against Materials Studio Amorphous Cell, Polymatic, PySIMM,
and RadonPy.

| Polymer | Chains × DP | Atoms | Tool | Build time | Density after build | Density after equilibration | C∞ | Rg (Å) | Close contacts |
|---|---|---|---|---|---|---|---|---|---|

**Table 4. Molecular dynamics throughput.**

| System | Atoms | Force field | Electrostatics | Hardware | CAPS (ns/day) | LAMMPS (ns/day) | GROMACS (ns/day) |
|---|---|---|---|---|---|---|---|

**Table 5. Energy conservation in NVE.**

| System | Time step (fs) | Constraints | Length (ns) | Drift (kT/ns/atom) | Reference engine drift |
|---|---|---|---|---|---|

**Table 6. Predicted properties against experiment.**

Cover at least PE, PP, PS, PMMA, PET, PC, PA6, and PDMS.

| Polymer | Force field | Density at 298 K (CAPS / exp.) | Tg in K (CAPS / exp.) | Solubility parameter (CAPS / exp.) | Bulk modulus (CAPS / exp.) | Experimental source |
|---|---|---|---|---|---|---|

**Table 7. Chain statistics against literature.**

| Polymer | Temperature (K) | C∞ (CAPS) | C∞ (literature) | Error (%) | Source |
|---|---|---|---|---|---|

**Table 8. Parallel scaling.**

| System | Atoms | Cores or GPUs | Wall time | Speed-up | Parallel efficiency (%) |
|---|---|---|---|---|---|

**Table 9. Cross-platform reproducibility.**

| Workflow | Linux | macOS | Windows | Bit-identical in deterministic mode |
|---|---|---|---|---|

**Table 10. 3D builder feature coverage.**

Compare against Materials Studio, Avogadro 2, GaussView, and Maestro. Base
every entry on the vendors' current documentation, and cite it.

| Feature group | Feature | CAPS | Materials Studio | Avogadro 2 | GaussView | Maestro |
|---|---|---|---|---|---|---|

**Table 11. 3D rendering and editing performance.**

| Structure | Atoms | Style | GPU | Frames per second | Load time (s) | Memory (MB) | Edit latency (ms) |
|---|---|---|---|---|---|---|---|

**Table 12. Builder accuracy.**

| Test set | Structures | Metric | CAPS result | Reference | Source |
|---|---|---|---|---|---|

Include at least these rows:
- R/S and E/Z assignment on the Hanson et al. (2018) CIP validation set.
- Bond-length and angle RMSD of cleaned geometries against experimental
  crystal structures.
- Space-group detection on a set of known crystal structures.

## Quality bar

- Keep golden outputs of every engine as a regression oracle.
- Use GoogleTest for C++, pFUnit for Fortran, pytest for Python, and xUnit
  for C#. Add property-based tests for graphs, typing, and periodic
  geometry.
- Add physics tests: analytic two-body and three-body energies, finite
  difference checks of every force, virial and pressure checks, NVE drift,
  and equipartition.
- CI on Linux, macOS (arm64 and x86), and Windows, with sanitizers,
  clang-tidy, coverage gates, and nightly benchmark runs.
- Packaging: conda-forge, PyPI wheels, a signed Windows MSI, a notarised
  macOS app, and Linux AppImage plus Flatpak.
- Choose the licence early. MIT, BSD-3, or LGPL-3 are typical. Because CAPS
  does not link LAMMPS or GROMACS, their GPL licences do not constrain it.

## Publication plan

Propose a publication strategy, for example:
1. A software paper describing CAPS and its architecture, for Journal of
   Open Source Software, SoftwareX, or Computer Physics Communications.
2. A methods paper on CAPS Pack and force-field-aware CAPS Grow with the
   benchmark tables, for J. Chem. Theory Comput. or J. Chem. Inf. Model.
3. An application paper validating polymer properties against experiment,
   for Macromolecules or Polymer.

## Migration strategy

Do not rewrite everything at once:

1. **Phase 0:** define the data model and C ABI, and the golden outputs the
   engines are tested against.
2. **Phase 1:** build CAPS Chem, CAPS Field, and CAPS IO. Pass Table 1 energy
   parity for single-point energies.
3. **Phase 2:** build CAPS Relax and CAPS Pack. Pass Table 2 against packmol.
4. **Phase 3:** build CAPS Grow and add force-field-aware growth. Pass Table 7.
5. **Phase 4:** build CAPS Dynamics and CAPS Equilibrate. Pass Tables 4 and 5.
6. **Phase 5:** build CAPS React and CAPS Analyze. Pass Tables 3 and 6.
7. **Phase 6:** build CAPS Studio, the 3D builder and visualiser, and the
   Avalonia GUI. Pass Tables 10, 11, and 12.
   Then release and publish.

Build the command layer behind CAPS Studio in Phase 1, together with CAPS
Chem, so every builder action is scriptable and tested long before the GUI
exists.

## What I want from you

Produce, in order:

1. An architecture document with component diagrams, the data model, the C
   ABI, and the threading and GPU model.
2. For each engine, a design specification with its equations, a
   literature reference for each, algorithms, complexity, parallel
   strategy, and validation tests.
3. A repository layout for a monorepo with CMake, dotnet, and a Python
   package, including build and CI configuration.
4. A module-by-module map of each CAPS engine to its language and
   dependencies.
5. The force-field database schema and the provenance manifest schema.
6. The design of the benchmark harness that produces Tables 1–12
   automatically.
7. A full CAPS Studio specification: the feature list grouped as above,
   the command architecture with undo and redo, the mouse and keyboard
   interaction map, the renderer design, and mockups of the main builder
   window and each specialised builder.
8. GUI information architecture and key screens for the rest of the
   Avalonia app.
9. A risk register covering numerical accuracy, performance targets,
   licensing, GPU portability, and interop between Fortran, C++, and C#.
10. A phased roadmap with measurable exit criteria, tied to the benchmark
   tables.

Ask clarifying questions only where the answer changes the architecture.
Otherwise state your assumptions and proceed.
