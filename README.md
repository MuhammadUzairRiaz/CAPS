# CAPS — Chain Assembly and Packing Suite

CAPS is a standalone polymer and soft-matter modelling suite: builders, packing, force fields, dynamics and
analysis, with a desktop Studio. It is independent software and lives in its own repository.

The product and engineering specification is `REDESIGN_PROMPT.md`. The screen designs are in `design/`
(open `design/index.html`); they are the reference for every Studio screen.

## What works today (v0.1)

| Part | Where | Status |
|---|---|---|
| Core data model, periodic (triclinic) cells, minimum image | `core/` (C++20) | working, tested |
| Readers: LAMMPS data (full / molecular / charge / atomic, image flags), LAMMPS text dump (x / xu / xs, ix iy iz, triclinic), GROMACS `.gro`, PDB (CRYST1, CONECT, models), XYZ / extended XYZ | `core/src/io_*.cpp` | working, tested |
| Writers: LAMMPS data (full), PDB, extended XYZ | `core/src/io_*.cpp` | working, tested |
| Bond perception (covalent radii, periodic cell list), whole molecules | `core/src/bonds.cpp` | working, tested |
| Analysis: molecule shape (Rg, gyration tensor, κ²), g(r) with cell list, distance / angle / dihedral | `core/src/analysis.cpp` | working, tested |
| Renderer: ball-and-stick, space filling, sticks, no-H, backbone; colour by element / molecule / type / distance; outlines, depth cue, picking | `core/src/render.cpp` | working |
| Figure export: PNG (dark, white or transparent with a real alpha channel) and SVG (no background shape when transparent) | `core/src/image.cpp`, `svg.cpp` | working, tested |
| Grow: all-atom polymer chains grown inside a periodic cubic cell. Any repeat unit written as SMILES with two attachment points (`*CC(*)c1ccccc1`); homopolymers and copolymers of up to eight units (alternating, block, random with a share per unit, gradient, pattern); tacticity by mirrored units; unit templates embedded and cleaned with GAFF2 (UFF for units GAFF2 cannot type, such as silicones); periodic contact checks against every placed atom, look-ahead, torsion refinement of near misses, back-tracking and restarts, and optional contact-scale stepping for crowded backbones. A curated polystyrene grower with GAFF names and Gasteiger charges | `core/src/polymer.cpp`, `grow.cpp` | working, tested |
| Crystals and surfaces: CIF reader (symmetry operations expanded, occupancies, bonds across the periodic cell), slabs cleaved along any (hkl) with every termination listed by the bonds it breaks per nm², rectangular surface cells (sheared ≤ 2 % when needed), supercells, vacuum, and passivation of dangling bonds (Si–OH, O–H, H). Library of fibre, cord and filler crystals: α-quartz, graphite, α-iron, β-brass, copper, zinc, zinc oxide, rutile, diamond, rock salt | `core/src/crystal.cpp`, `data/crystals` | working, tested |
| Space groups and the crystal builder: the 530 settings of the 230 space groups from their Hall symbols (International Tables Vol. B, Table A1.4.2.7) closed into their operations; a crystal from a space group, a lattice and an asymmetric unit (images within a tolerance merged, each site's multiplicity reported); sites snapped onto their special positions; the space group and asymmetric unit found for any crystal or CIF (the highest-order setting that maps it onto itself, origin shifts tried); primitive cells of A, B, C, I, F and R lattices; supercells | `core/src/spacegroup.cpp`, `core/src/spacegroups.inc` | working, tested |
| Peptides: all-atom peptides from a one-letter sequence (or FASTA) with a secondary structure per residue — backbone placed by NeRF from φ, ψ, ω (α-helix, β-strand, PPII, a seeded coil; proline at φ −65°), side chains and hydrogens superposed from each L-amino acid embedded from SMILES, side-chain charges at the pH, NH₃⁺ / NH₂ / acetyl and COO⁻ / COOH / N-methyl-amide termini, UFF clean-up; PDB atom names, residue names and numbers | `core/src/peptide.cpp` | working, tested |
| Solvation: a structure held at the centre of a periodic box (cubic, rectangular or its extent plus padding) with solvent and ions packed around it by CAPS Pack — water as TIP3P, SPC/E or TIP4P/2005 (the M-site charge on O, as LAMMPS' tip4p styles expect) or a swelling solvent built from SMILES (toluene, cyclohexane, n-hexane, chloroform, THF, acetone, methanol, ethanol, DMSO); the count from the density over the free volume (the box less the solute's van der Waals volume on a 0.5 Å grid); ions to neutralise, at a salt concentration (NaCl, KCl, LiCl, CaCl₂, MgCl₂, ZnCl₂) or as counts | `core/src/solvate.cpp` | working, tested |
| Appearance: styles per selection (ball and stick, stick, wireframe, space filling, coordination polyhedra — SiO₄, TiO₆ — and a ribbon along the backbones) chosen by pipeline expressions; colour by partial charge on a colour-blind-safe diverging ramp; atom labels (element, CIP R/S, type, charge); solvent-accessible, solvent-excluded and van der Waals surfaces by marching tetrahedra, coloured by the Coulomb potential of the charges or by the nearest atom, translucent | `core/src/appearance.cpp`, `core/src/render.cpp` | working, tested |
| Trajectory player data: LAMMPS logs read in both thermo styles (one-line tables and multi blocks, every run), per-frame series (time, density, volume, a chain's Rg and backbone end-to-end distance) with the log's columns matched to the frames by timestep, run starts as frames, positions smoothed over a window of frames | `core/src/trajectory.cpp` | working, tested |
| Torsion scans: a dihedral set by rotating the atoms on one side of its bond (ring bonds refused), rigid or with the four atoms held and the rest minimised at each step, with the document's force field or UFF; minima refined by a parabola and named trans, gauche±, anticlinal± or cis; the barrier; every point's geometry kept to show; a default backbone torsion | `core/src/torsion.cpp` | working, tested |
| Structure editing: atoms placed in the ideal direction for the parent's hybridisation, bonds drawn and broken, atoms deleted, elements and charges changed, missing hydrogens filled by valence (formal charge counted), a stereocentre inverted by swapping two branches, vinyl-polymer tacticity read as m/r dyads and mm/mr/rr triads and set iso- or syndiotactic, UFF clean-up of the picked atoms; selections by SMARTS, element, type, charge, distance, bonds or expression; undo and redo | `core/src/edit.cpp` | working, tested |
| Interactions and checks: hydrogen bonds by the Luzar–Chandler geometric criterion (D···A ≤ 3.5 Å, H–D···A ≤ 30°), close contacts below the van der Waals sum − 0.4 Å and clashes below 0.75 × the sum (1-2 and 1-3 pairs excluded), molecules cut by the cell edge, atoms over their valence or lacking hydrogens, the net charge — each with a fix (push apart with UFF, wrap, add H) | `core/src/interactions.cpp` | working, tested |
| Fragment library: 129 fragments (rings, functional groups, monomers — NR, BR, CR, NBR, IIR, silicone …, amino acids, bases, sugars, solvents, ions, rubber additives: sulfur, TMTD, MBTS, CBS, MBT, DPG, 6PPD, TMQ, stearic acid, dicumyl peroxide, Si69 silane, RFL resins) as SMILES with * attachment points; attached in place of a hydrogen of the picked atom, rolled clear of the structure, other points capped with H; whole molecules placed beside the structure | `core/src/edit.cpp`, `data/fragments` | working, tested |
| Import of files without a topology (XYZ, extended XYZ, PDB, CIF): bonds perceived from distances (covalent radii + a tolerance), read from the file or none; bond orders and aromaticity from valences; molecules by connectivity; molecules made whole across the cell; the cell kept or dropped; a preview on frame 0 with a ten-heavy-atom fragment | `core/src/import.cpp` | working, tested |
| Image and movie export: PNG at 8 or 16 bits per channel (16-bit keeps the supersampled average), dpi, sRGB, a provenance manifest (caps-image/1.0: source file and its sha256, frame, camera, render settings) in an iTXt chunk, labels and measurements laid over at full precision; movies as animated PNG (streamed) or PNG sequences of the frames or a turntable, MP4 through ffmpeg when installed | `core/src/image.cpp`, `capi/src/caps_c.cpp` | working, tested |
| Provenance: every document keeps the steps that produced it — engine, parameters, generator and seed, the papers behind each method (24 built-in BibTeX entries: UFF, QEq, GAFF, L-BFGS, FIRE, CSVR, stochastic cell rescaling, SPME, DSF, the 21-step protocol, PACKMOL, TIP4P/2005 …) and the approximations in force — with the inputs' sha256; saved beside the file as FILE.provenance.json (caps-manifest/1.0) and read back on open; runs compared step by step; `caps provenance FILE [--json | --bibtex | --methods | --compare OTHER]`; a methods paragraph for a paper written from the steps with numbered references | `core/src/provenance.cpp` | working, tested |
| Voids: a periodic distance field to the nearest Bondi surface, its maxima as the largest non-overlapping empty spheres, the share a probe reaches, spheres drawn translucent, PDB export (radius in the B-factor); neutron contrast by deuteration (every H, aliphatic H — a deuterated backbone —, aromatic H, H on O/N) | `core/src/voids.cpp`, `core/src/properties.cpp` | working, tested |
| Pores: a slit between graphene walls (1–3 AB sheets each, width H carbon centre to centre, periodic or with vacuum), a cylindrical channel carved from a crystal block (hexagonal lattices use the rectangular a × (a + 2b) cell; carved surface optionally passivated into the pore), or a framework supercell (zeolite, MOF); a fluid packed inside the pore only; walls are molecule 1; `caps pore` | `core/src/nano.cpp` | working, tested |
| Theory manual: 24 pages, one per implemented method (integrator, thermostats, barostats, electrostatics, minimisers, push-off, chain growth, packing, UFF, QEq, GAFF, Gasteiger, the 21-step protocol, structure factors, free volume, diffusion, elastic constants): summary, equation, symbols with CAPS's settings, when to use it, references, source and tests, and every deviation from the paper | `data/manual/manual.json` | written from the sources |
| Coarse-grained melts: Kremer–Grest bead-spring chains (M × N at ρσ³, optional bending k_θ) built as random walks, drawn, and written for LAMMPS in reduced units — soft push-off ramped with limited displacements, then FENE + WCA under Langevin at T = 1 (checked in LAMMPS: P ≈ 5.0 at ρσ³ = 0.85) | `core/src/kremer_grest.cpp` | working, tested |
| Large systems: level of detail in the renderer by distance from the focus (full atoms and bonds near, spheres without bonds, then points; ambient occlusion near only) — a million beads with bonds in 90 ms instead of 170 ms on the CPU —; per-tier counts; the memory a document holds; g(r) on 20 000 strided centres above that size | `core/src/render.cpp`, `core/src/analysis.cpp` | working, tested |
| Python bindings: the `caps` package (ctypes over libcaps) — `caps.open`, `caps.import_file`, `doc.provenance()`, `caps.provenance_file`, `caps.compare_provenance`, `caps.bibtex`, `caps.build.smiles / peptide / crystal / nano / solvate`, `doc.field.assign`, `doc.relax`, `doc.md`, `doc.edit` (every builder operation), undo, selections, tacticity, interactions, torsion scans, series, save and render; the library found through `CAPS_LIB` or beside the package | `data/python/caps/core.py` | working, tested (`tests/python`) |
| Interfaces: a polymer film grown onto a slab (any repeat units and sequences; film thickness, density, gap, vacuum or a periodic sandwich), the slab held in place while chains grow and in Relax — fibre–rubber and filler–rubber interfaces | `core/src/interface.cpp` | working, tested |
| Nanostructures: graphene sheets and H-capped flakes, (n, m) carbon nanotubes (periodic or capped), spheres, cubes, octahedra, cuboctahedra and periodic fibres cut from any CIF (passivated silica, metal clusters, glass fibres); a filler held in a periodic cell while polymer chains grow around it (filled rubber, CNT– and graphene–rubber composites) | `core/src/nano.cpp` | working, tested |
| Blends: two or more polymers (NR/BR, SBR/BR tyre compounds …) grown component after component into one cell, chain counts from weight fractions, mixed, two-slab or droplet start | `core/src/interface.cpp` (grow_blend) | working, tested |
| Polymer library (all 111 units grow: 82 at full contact limits, 29 with the contact scale stepped down): 111 repeat units as SMILES and 16 copolymer presets (natural rubber, ENR-25/50, high-cis BR, NBR, butyl, chloroprene, E-SBR, SBR, SAN, EVA …) | `data/polymers/library.json` | data |
| Molecule builder: SMILES parser and writer (chirality, E/Z, rings, brackets), 2D depiction, distance-bounds embedding (4D → 3D) with chirality, planarity and square-planar / trigonal-bipyramidal / octahedral centres, GAFF2 or UFF clean-up (UFF whenever the chosen force field cannot type every atom), conformers ranked by energy | `core/src/smiles.cpp`, `embed.cpp` | working, tested |
| Field (first slice): GAFF 1.81 typing for hydrocarbons (c3, ca, hc, ha) with the reason for each type; bonds, angles, Fourier torsions, impropers, LJ with arithmetic mixing, damped shifted force electrostatics or particle-mesh Ewald, AMBER 1-4 scaling; periodic neighbour list that includes images in cells narrower than twice the cut-off; multithreaded pair terms | `core/src/field.cpp` | working, tested against LAMMPS |
| UFF, the Universal Force Field (Rappé et al. 1992), for every element H–Lr: typed from elements and bonds (hybridisation from the steric number, conjugation, oxidation state, metal coordination geometry), UFF's bond, Fourier / periodic angle, torsion and inversion rules and full 1-4 van der Waals; clean-up in the molecule and polymer builders, and a force field for Field, Relax, Dynamics and LAMMPS export | `core/src/uff.cpp` | working, tested against LAMMPS |
| QEq charge equilibration (Rappé & Goddard 1991) for every element, Ohno–Klopman shielding, sparse conjugate gradients: `--charges qeq`, Field "QEq (every element)", and with UFF | `core/src/qeq.cpp` | working, tested |
| Particle-mesh Ewald: smooth PME (Essmann et al. 1995) — reciprocal part in Fortran 2018 (B-spline spreading, mixed-radix FFT, influence function, forces, virial tensor; the plain Ewald sum as a reference), run on the evaluator's threads; real space, self and exclusion terms in C++; agrees with LAMMPS's Ewald sum | `core/fortran/caps_kspace.f90`, `core/src/kspace.cpp` | working, tested against LAMMPS |
| Relax: steepest descent, Polak–Ribière CG, L-BFGS (m = 10) and FIRE; capped-force push-off; affine compression to a target density; isotropic box relaxation to a pressure; one frame recorded per stage | `core/src/relax.cpp` | working, tested |
| LAMMPS data export with the force field (coefficients, angles, dihedrals, impropers, velocities and the matching styles); LAMMPS dump export of whole trajectories | `core/src/relax.cpp`, `io_lammps.cpp` | working, tested |
| Dynamics: velocity Verlet; NVE, NVT (Bussi velocity rescaling, Langevin BAOAB) and isotropic NPT (stochastic cell rescaling, Berendsen); LJ tail corrections; velocities carried between runs; thermo log; multithreaded | `core/src/dynamics.cpp` | working, tested against LAMMPS |
| Equilibrate: the Larsen et al. 21-step compression / decompression protocol, simulated annealing cycles, MD push-off with capped LJ forces, custom plain-text protocols (NVT / NPT / NVE stages, temperature ramps, force caps); production blocks until density, energy and Rg converge | `core/src/equilibrate.cpp` | working, tested |
| Chain statistics: backbone detection (non-ring heavy-atom path), mean-square internal distances ⟨R²(n)⟩/(n⟨b²⟩), end-to-end distance | `core/src/analysis.cpp` | working, tested |
| Pack: rigid molecules packed into boxes, cubes, spheres, cylinders and half-spaces (inside / outside), fixed molecules, periodic or not; overlap penalty after Martínez et al. over centre + rotation, L-BFGS, multithreaded cell list, worst molecules moved between rounds; reads packmol input files; never returns a cell below tolerance | `core/src/pack.cpp` | working, tested, benchmarked against packmol |
| React: atom-mapped templates (form, break, move, delete); distance capture and probability (REACTER style); cycles of react → retype → minimise → optional dynamics (Polymatic cycle); conversion counted per reactive group; cluster analysis, gel point from the reduced weight-average mass, Flory–Stockmayer α_c. Rubber cures: sulfur (H–S–S–H donors inserted into the cell, allylic C–S bonds giving C–S–S–C crosslinks) and peroxide (allylic C–C) | `core/src/react.cpp` | working, tested |
| Analyze: density, g(r), S(q) (direct reciprocal-lattice sum + g(r) transform), X-ray and neutron scattering, Rg, end-to-end distance, C_n and C∞, persistence length, MSD, D (Einstein), end-to-end and segmental relaxation (KWW), cohesive energy density and δ, free volume by probe insertion, pore size distribution; interfaces: density profile along z (surface / film, first-layer peak), work of adhesion, chain orientation (S, Herman's f, P₂ against height), and steered-MD pull-out (interfacial shear strength, peak normal stress, work of separation); vulcanised networks: sulfur bridges by rank, pendant groups, crosslink density and Mc; block-average errors; JSON / CSV output | `core/src/properties.cpp` | working, tested against MDAnalysis and analytic cases |
| C ABI v16 for the GUI and other languages (`capi/include/caps_c.h` lists every call and the version that added it) | `capi/` | working |
| `caps` command line: info, render, shape, rdf, convert, build (SMILES → 3D), grow, pack, contacts, field, relax, md, equilibrate, chains, react, ff, analyze, elastic, tensile, tg, bench | `cli/` | working |
| Recipes: `caps run RECIPE.yaml` makes a structure end to end from one YAML (or JSON) file — build (polymer / molecule / file) → type (default, UFF, a library id such as `gaff2`, or a path; charges from the force field, Gasteiger or QEq) → grow → relax → md → equilibrate → analyze → export (LAMMPS data + input, GROMACS, PDB, xyz, mol2) with the provenance manifest beside the first file. One line per stage (`[k/n] stage  detail  done`, a bar while it runs) or `--json` lines; `--seed`, `--threads`, `--out DIR`. Exit codes 0 ok · 2 input · 3 missing parameters · 4 failed run. `samples/ps_cell.yaml` is the board's recipe | `core/src/recipe.cpp`, `core/src/yaml.cpp` | working, tested |
| Notebooks: `caps.polymer(smiles, dp, tacticity, seed)`, `caps.run(recipe)`, `doc.view(style=…)` — an interactive view drawn by the page itself (drag to rotate, wheel to zoom; no widget extension; the renderer's PNG where HTML is not shown), `doc.analyze([...])`, `doc.provenance.citations(fmt="bibtex")` / `.methods()`, `caps.sweep.run(…)` / `caps.sweep.result(condition)` (the Studio's Sweep folders, `results.json`), `caps.table(cells, ["density", "tg", "c_inf"])` (mean ± SD over seeds, HTML / text / CSV). `samples/ps_tacticity.ipynb` is the board's notebook | `data/python/caps/` | working, tested |
| Accessibility and system states: colour-vision check (Machado 2009 simulation, ΔE*ab pairs, the view previewed per deficiency; `core/src/colourvision.cpp`), Reduce motion (follows the OS; camera flights of 450 ms or cuts; F frames the selection), screen-reader names for every icon button, notices with severity (info · check · warning · blocks run), a manual update check against the GitHub releases, and `site/` — the static download page | `core/src/colourvision.cpp`, `studio/` | working, tested |
| Analyze focus pages: Glass transition (stepwise NPT cooling of the open structure, replicas by seed pooled with error bars, two-line fit, α glass / α melt, the cooling rate beside every Tg; `analyze: {properties: [tg]}` in recipes), Interface (density along the surface normal per surface element and for the film, gap, adsorbed layer, plateau, first-layer peak, work of adhesion — this frame and over the trajectory), Diffusion (log–log MSD with its fit window, local slope check, Einstein D in m²/s, Yeh–Hummer finite-size correction) | `core/src/properties.cpp`, `studio/` | working, tested |
| History & snapshots: every edit with its atom count and "you are here"; editing after an undo keeps the undone steps as a branch; named snapshots to restore or compare with the structure now (`caps_snapshot`) | `capi/src/caps_c.cpp`, `studio/` | working, tested |
| Bench: the built-in validation suite (forces, virial and PME against finite differences and Ewald, packing, MD throughput, NVE drift, thread scaling, reproducibility, rendering, the molecule builder) with Markdown, CSV and LaTeX export | `core/src/bench.cpp` | working |
| CAPS Studio (Avalonia, .NET 10), one page per design board: Start (quick start, builders, recent work), Studio (3D view, inspector, measurements, frames, figure export), Molecule builder (2D sketcher, SMILES, conformers), Polymer builder (library, units, sequences, tacticity), Fragment library (categories, search, 3D tiles, attach points ringed, attach to the picked atom or place, My fragments; the Studio sidebar's quick tiles), Analyze focus pages — Mechanics (E, ν, K, yield, the stress–strain curve with its fit window, the stiffness matrix, protocol: constant strain or fluctuations, rate, lateral pressure, axis), Scattering (X-ray and neutron I(q), peaks in q, d and 2θ for Cu/Mo/Co/Cr/Ag Kα, scattering lengths, isotope pattern, a measured curve overlaid from q or 2θ), Free volume (the largest voids drawn in the cell, Bondi FFV, accessible share, pore size distribution, voids as PDB) — with the calculations list beside them, Parameter sweep (tacticity × chain length × seed of any homopolymer, each run grown, relaxed, optionally run in NPT, analysed for density and Rg and saved with its provenance; a few at a time with the threads shared; the run grid coloured by state, results by condition as mean ± sd, cell sizes, pause and stop), Project home (the folder's structures with thumbnails, atoms, density, last step and status from provenance; a methods section written from the selected structure's provenance with replicas detected and numbered references, Copy, BibTeX; Share project as a zip with provenance, methods.txt and references.bib), Theory manual (method pages with the equation, symbols and CAPS settings, references, implementation and tests, deviations from the paper, the open structure's steps that used the method, Copy BibTeX), Compact layout below 1440 px (rail and toolbar view options as icons, the inspector and project as drawers over the view from side tabs, the curves dock folded to its tabs), First-run tour (five steps lighting the rail, toolbar, 3D view, inspector and curves with a callout and the real shortcuts; once with the first structure, again from ⌘K), Provenance viewer (the steps as a timeline with parameters, seeds and citations, compared with another run with the differing step marked, approximations, inputs with sha256, Export BibTeX, View JSON; from Jobs, the File tab and ⌘K), Export dialog (⌘E; Image: Screen, 4K, journal and poster presets, 8- or 16-bit PNG or SVG, background, supersampling, colour profile, labels and measurements, provenance in the metadata, copy to the clipboard; Movie: trajectory frames or a turntable, animated PNG, PNG sequence or MP4 with ffmpeg, frame rate and stride, progress and Stop), Import dialog (XYZ, PDB and CIF files: format and bonds detected, the first lines, atoms, cell and units, a bond-perception preview on a fragment, bonds perceive / file / none with the tolerance, bond orders, split, unwrap, cell, then a force field, charges and the checks), Macro recorder (every Studio action — open, build, edit, fragments, tacticity, field, relax, save — recorded as a line of Python using the caps package; script files in ~/.caps/macros; a selected literal promoted to a parameter of `macro(…)`; run with python3 and this Studio's library, output streamed, Stop), Split view (the open structure beside a second one, cameras synced, selection and frame linked by index, atoms / density / stereocentres / tacticity / Rg compared, the Paper theme a click away), Interactions & checks panel (H-bonds dashed and clashes red in the view, contacts, the validation list with Push apart / Wrap / Add H and Fix all, CSV), Selection & stereo panel (select by SMARTS, element, type, charge range, distance or along bonds with replace / add / subtract / invert, named sets, stereocentres and tacticity with the m/r dyad strip, make iso- or syndiotactic), Builder tools in the Studio (place atom, draw bond, delete, +H, invert, auto-clean, undo / redo with ⌘Z, the Element picker — a periodic table with charge and placement geometry, ⇧E), Torsion scan (four picked atoms or a backbone torsion, the energy curve with its scan points, conformers found — click one to show it —, the barrier, CSV), Trajectory player (speed, loop, smoothing, run starts on the timeline, a chain's end-to-end vector in the view, three plots linked to the frame with its value, a LAMMPS log beside the dump picked up), Appearance panel (style tiles for all atoms or an expression, layers, colour by charge with its legend, labels drawn over the view, surfaces, AO / outlines / depth cue), Solvation builder (box, solvent and water model, ions by concentration or to neutralise, counts from the free volume, a sparse sample in the preview, packing stages), Biomolecule builder (sequence grid coloured by structure, a structure applied to a selected range, φ ψ ω per structure, termini, protonation, a tube through CA in the preview), Crystal builder (space-group search, lattice fixed by the crystal system, asymmetric unit table, CIF import through Find symmetry, Apply symmetry, primitive cell, supercell drawn dashed around one unit cell), Surface & interface builder, Performance panel (level-of-detail tiers, measured frame time and the HUD, memory per atom, adaptive near tier to hold 16.7 ms, benchmark of the view; switched on for documents of 200 000 atoms or more; exports stay in full detail), Reaction template editor (the atom-mapped template drawn before and after with map numbers, reacting atoms ringed, formed and broken bonds coloured; bond changes; checks; a test on the open structure — reactive sites and matches within capture —; saved to ~/.caps/templates; Use in React), Coarse-grained builder (Kremer–Grest melts, LAMMPS deck), Nanostructure builder (sheet, nanotube, particle, pore: slit / cylinder / framework with a fluid), Blend builder, Grow, Pack, Relax, Dynamics, Equilibrate (convergence), React, Analyze (properties, mechanics, Tg), Field (typing and parameters), Jobs (runs with curves, log, provenance), Bench, Settings (theme, palettes, electrostatics, threads; Compute & remote: SSH hosts with scheduler, partition and work directory, a connection test that runs caps on the host through the SSH agent, the editable SLURM job template), File checks (every opened file checked, one action per problem), Analyze empty state (drop zone, recent files, samples), progressive open of large trajectories (frame 0 first, the other frames read in the background, Cancel keeps what was read), Analyze › Visualize (a pipeline of non-destructive steps run on every shown frame — expression selection, expand/invert/clear selection, delete, slice, colour coding by any property, assign colour, cluster analysis, coordination and g(r), compute property, wrap, replicate, histogram, spatial binning, create bonds, unwrap, molecule shape, topology distributions, displacements, smoothing, vectors and trajectory lines drawn as tubes and arrows, voids and pores with a probe sweep, Voronoi volumes (grid or radical), a mass-conserving density field with a slice — with a data inspector for particles, bonds, global attributes and data tables; the same JSON runs with `caps pipeline`), Analyze › Compare (two inputs side by side with linked cameras, the same pipeline on both, every attribute as A, B and B − A, a data table of both overlaid), Analyze › Batch (the Visualize pipeline over many inputs — a glob such as cells/seed*/PS_melt.data or picked files — on each input's first or last frame, a few at a time; failures kept as rows; mean ± s.d. and a per-input plot; batch/results.csv with the pipeline's sha256), Export › Figure bundle (a zip with the figure, the CSV behind every plot, the pipeline, the input and provenance.json with every file's sha256; `caps reproduce` rebuilds the data and compares the hashes), four viewports (top, front, left, perspective), a Colour-by gallery, Export › Data (LAMMPS data with force-field sections, LAMMPS dump, GROMACS .gro, PDB, extended XYZ, mol2 — the frame or the Visualize result — with a preview of the real file, its sections and size), Studio › Render (output-frame guide in the view, ambient occlusion, depth cue, antialias, overlays filled from live attributes — text label with [Title] [SourceFrame] [Particles] [Density] …, vertical colour legend, true scale bar, axis tripod — one image or an image sequence of the frames), Export › Figure (dark, white and transparent side by side; journal, poster and slide sizes in mm × dpi with the dpi written into the PNG; title, true scale bar and colour legend in ink that suits the background, also in SVG), keyboard walk of the 3D view with screen-reader announcements (↑ ↓ atoms, [ ] molecules, B bonds, Space select, M measure, ⌘⇧A announce), command palette (⌘K) | `studio/CapsStudio` | working |

Not built yet (see the roadmap in `REDESIGN_PROMPT.md`): branched, star and comb architectures; nucleic-acid strands; learning reaction templates from reactant / product pairs; per-atom constraints and ellipsoids in Pack; GPU kernels; double-bridging / end-bridging Monte Carlo; Nosé–Hoover chains and MTK; anisotropic cells; constraints (SHAKE / RATTLE / LINCS); r-RESPA; restraints and fixed atoms in Relax; the Vulkan viewport; coarse-grained dynamics inside CAPS and backmapping to all atoms; provenance for runs started from the command line (the Studio and Python record it).

## Build

Requirements: CMake ≥ 3.24, a C++20 compiler, a Fortran 2018 compiler (gfortran: `brew install gcc`, `apt install gfortran`),
zlib, and the .NET 10 SDK (`brew install dotnet`). The Fortran kernels need no Fortran runtime, so the built library
depends on nothing beyond the system. On Windows the native core builds with MinGW-w64 (MSYS2 UCRT64: gcc, g++,
gfortran, cmake, ninja), statically linked.

```bash
scripts/build.sh      # C++ core + tests, then the Studio + its self-test
scripts/studio.sh     # launch the Studio
scripts/studio.sh samples/ps_melt.lammpstrj samples/ps_melt.data
```

Command line:

```bash
build/cli/caps info samples/ps_melt.data
build/cli/caps render samples/ps_melt.data -o figure.png --bg transparent --size 1920x1080
build/cli/caps render samples/ps_melt.data -o figure.svg --bg white --colour molecule
build/cli/caps rdf samples/ps_melt.data --pair C-C --inter
```

Build a molecule and grow polymers (any repeat unit, `*` marks the head then the tail):

```bash
build/cli/caps build "N[C@@H](C)C(=O)O" -o alanine.mol2 --conformers 5 --ff data/forcefields/gaff-amber25-dlfield.json
build/cli/caps grow --units '[*]C/C=C(C)\C[*]' --chains 10 --dp 30 --density 0.5 -o natural_rubber.data
build/cli/caps grow --units '[*]C/C=C\C[*],*CC(*)c1ccccc1' --sequence random --weights 0.86,0.14 --chains 10 --dp 30 --density 0.5 -o sbr.data
build/cli/caps grow --units '*CC(*)(C)C(=O)OC' --chains 10 --dp 20 --density 0.5 --auto-scale -o pmma.data  # crowded backbones: contacts lowered as needed, then relax
build/cli/caps surface data/crystals/alpha-quartz.cif --hkl 0,0,1 --list                                     # terminations of quartz (001)
build/cli/caps surface data/crystals/alpha-quartz.cif --layers 3 --supercell 5,3 --passivate -o quartz_001.data
build/cli/caps interface data/crystals/alpha-quartz.cif --supercell 5,3 --layers 2 --passivate \
    --units '[*]C/C=C(C)\C[*]' --dp 12 --film 25 --film-density 0.9 -o nr_on_silica.data              # rubber film on a glass-fibre surface
build/cli/caps relax nr_on_silica.data --ff uff --fix-mol 1 -o nr_on_silica_relaxed.data                   # UFF, the surface held in place
build/cli/caps blend --components '[*]C/C=C(C)\C[*],[*]C/C=C\C[*]' --weights 0.7,0.3 --chains 8 --dp 30 -o nr_br.data  # NR/BR blend
build/cli/caps peptide AEAAAKEAAAKEAAAKA --helix --n-term ACE --c-term NME -o eaaak.pdb   # a capped α-helical peptide
build/cli/caps solvate eaaak.pdb --padding 8 --model TIP3P --salt NaCl --conc 0.15 -o eaaak_water.data   # peptide in 0.15 M NaCl
build/cli/caps solvate nr_cell.data --box 60,60,60 --solvent toluene --no-ions -o nr_toluene.data         # rubber swelling in toluene
build/cli/caps crystal --group Pnam --cell 7.40,4.93,2.534 --sites 'C1 C 0.038 0.065 0.25; H1 H 0.1848 0.0466 0.25; H2 H 0.0068 0.2811 0.25' --supercell 4,6,20 -o pe_crystal.data   # orthorhombic PE
build/cli/caps crystal data/crystals/rutile.cif --find-symmetry          # P 42/m n m (No. 136) and its two sites
build/cli/caps crystal --groups 14                                        # the settings of P2₁/c
build/cli/caps nano tube --n 10 --m 10 --length 25 --units '[*]C/C=C(C)\C[*]' --chains 12 --dp 20 --density 0.9 -o nr_cnt.data  # CNT in rubber
build/cli/caps nano particle data/crystals/alpha-quartz.cif --radius 10 --passivate --units '[*]C/C=C(C)\C[*]' --chains 30 -o nr_silica.data
build/cli/caps pull nr_on_silica_relaxed.data --distance 20 --rate 2 --csv pull.csv                           # pull-out: interfacial shear strength
build/cli/caps nano particle data/crystals/alpha-quartz.cif --shape fibre --radius 8 --passivate --units '[*]C/C=C(C)\C[*]' --chains 20 -o fibre.data
build/cli/caps pull fibre.data --axis z --distance 15 --rate 2                                                   # fibre pull-out (area 2πRL)
build/cli/caps analyze nr_on_silica_relaxed.data --props zprofile,adhesion                                  # film density along z, work of adhesion
build/cli/caps react nr_relaxed.data --template sulfur_allylic --insert SS --count 40 -o nr_vulcanised.data  # sulfur cure: C–S–S–C crosslinks
```

Particle-mesh Ewald instead of the damped shifted force (periodic cells):

```bash
build/cli/caps md cell.data -o out.data --steps 20000 --pme                 # β from erfc(β rc) = 1e-5, 1 Å grid, order 5
build/cli/caps md cell.data -o out.data --pme --ewald-rtol 1e-6 --pme-spacing 0.8 --pme-order 6
python3 bench/ff/check_data_lammps.py --pme                                   # against LAMMPS's Ewald sum
```

Grow a polystyrene cell with the curated grower:

```bash
build/cli/caps grow -o cell.data --chains 10 --dp 8 --density 0.4 --tacticity atactic --seed 1
build/cli/caps grow -o cell.data --chains 10 --dp 8 --density 0.7 --scale 0.75   # denser, needs minimisation later
```

What to expect from Grow today. Growth enforces contact limits (C–C 3.0, C–H 2.45, H–H 2.0 Å) with no
force field, so it cannot reach melt density on its own; Relax (below) compresses the grown cell.
Measured on this machine, five seeds each:

| Cell | Contact scale | Success | Time |
|---|---|---|---|
| 10 × DP 8, 0.4 g/cm³ | 1.0 | 5 / 5 | 0.7 s |
| 10 × DP 8, 0.5 g/cm³ | 1.0 | 4 / 5 | 1.9 s |
| 10 × DP 8, 0.7 g/cm³ | 0.75 | 1 / 1 tried | 4 s |
| 20 × DP 20, 0.4 g/cm³ | 1.0 | 3 / 5 | 10 s |
| 20 × DP 20, 0.5 g/cm³ | 0.85 | 2 / 5 | 18 s |

The Studio retries up to three consecutive seeds and says which one it used.

## Relax

Relax types the structure with GAFF, pushes overlaps apart with capped LJ forces, compresses stage by stage to a
target density (minimising after each stage) and minimises to a force tolerance. Grow a loose cell, then compress it:

```bash
build/cli/caps grow  -o grown.data --chains 20 --dp 20 --density 0.3 --seed 1
build/cli/caps relax grown.data -o melt.data --density 1.05             # L-BFGS, |F|max < 0.5 kcal/mol/Å
build/cli/caps relax grown.data -o cell.data --box-relax --pressure 1    # 0 K volume at 1 atm
build/cli/caps field melt.data                                           # types, terms, energy, pressure
```

The `.data` output carries the force field; it runs in LAMMPS with the styles written in its header
(`pair_style lj/cut/coul/dsf 0.2 10`, `bond_style harmonic`, `angle_style harmonic`, `dihedral_style fourier`,
`improper_style cvff`, `special_bonds amber`).

Measured on this machine (10-core Apple silicon), Grow followed by Relax to 1.05 g/cm³:

| Cell | Atoms | Relax | Final largest force |
|---|---|---|---|
| 10 × DP 8, grown at 0.4 g/cm³ | 1 300 | 1.8 s, 1 352 iterations | 0.47 kcal/mol/Å |
| 20 × DP 20, grown at 0.3 g/cm³ | 6 440 | 15 s, 2 648 iterations | 0.48 kcal/mol/Å |

Checked against LAMMPS (`scripts/check_lammps.sh`, needs `lmp`): bond, angle, dihedral, improper and vdW energies and
the pressure agree to printed precision, and per-atom forces to 2 × 10⁻⁵ kcal/mol/Å. The Coulomb energy differs by a
constant per pair because CAPS uses the Fennell–Gezelter form, which is zero at the cut-off; the forces are the same.

Limits of this version: a minimised cell at 1.05 g/cm³ is not an equilibrated melt. The 0 K pressure stays in the
thousands of atm, and chain conformations are those of the growth; run Dynamics (NPT) to equilibrate. Only C and H are
typed; other elements stop with a message naming the atom.

## Dynamics

```bash
build/cli/caps md melt.data -o hot.data  --steps 40000 --temp 500 --barostat crescale --pressure 1 --log hot.csv
build/cli/caps md hot.data  -o cold.data --steps 40000 --temp 300 --barostat crescale --dump cold.lammpstrj --every 1000
build/cli/caps md cold.data -o nve.data  --steps 10000 --thermostat none        # continues with the saved velocities
```

Units are LAMMPS `real` (fs, Å, kcal/mol, K, atm). The time step defaults to 1 fs with no constraints. Temperatures use
3N − 3 degrees of freedom. Pressure includes the kinetic term and the LJ tail correction. The `.data` output carries
velocities, so a run can continue in CAPS or in LAMMPS.

Checked against LAMMPS (`scripts/check_lammps.sh`): a 200-step NVE trajectory from the same positions and velocities
agrees to 10⁻⁵ Å (the precision of the files); with `pair_modify tail yes` the energies and pressure agree. An NPT run
at 500 K and 1 atm from the same start (CAPS stochastic cell rescaling, LAMMPS `fix npt` Nosé–Hoover, 40 ps each)
gives 0.87 and 0.86 g/cm³ over the second half, the same within the fluctuations of a 1 300-atom cell.

Speed on this machine (10-core Apple silicon, cut-off 10 Å, skin 2 Å, 1 fs):

| Cell | Atoms | 1 thread | 4 threads | 10 threads |
|---|---|---|---|---|
| 10 × DP 8 | 1 300 | 20 ns/day | 59 ns/day | 78 ns/day |
| 20 × DP 20 | 6 440 | | | 16 ns/day |

For comparison, this machine's serial LAMMPS build ran the 1 300-atom NPT cell at 16 ns/day.

Equilibration example (10 × DP 8, from Relax at 1.05 g/cm³): 40 ps NPT at 500 K gives 0.87 g/cm³; then 40 ps at 300 K
gives 0.976 g/cm³ at +99 atm mean pressure. These are short oligomers quenched quickly, not a converged glass.
Long-chain melts need ns-scale runs at high temperature; see Equilibrate below.

Studio self-checks (no window needed):

```bash
studio/CapsStudio/bin/Release/net10.0/CapsStudio --selftest samples out
studio/CapsStudio/bin/Release/net10.0/CapsStudio --screenshot out/studio.png samples/ps_melt.lammpstrj samples/ps_melt.data
```

## Equilibrate

```bash
build/cli/caps equilibrate melt.data --protocol larsen21 --print-protocol          # the 21 stages as text
build/cli/caps equilibrate melt.data -o eq.data --protocol larsen21 --tmax 600 --until-converged --log eq.csv
build/cli/caps equilibrate melt.data -o eq.data --protocol annealing --cycles 3 --tlow 300 --thigh 600
build/cli/caps equilibrate melt.data -o eq.data --protocol my_protocol.txt       # custom stages
build/cli/caps chains eq.data                                                   # internal distances
```

Protocols are lists of Dynamics stages. The same text format is used for printing, editing in the Studio and custom
files:

```
nvt 50 ps T 600                # 1 · heat
npt 50 ps T 300 P 986.92 atm   # 3 · compress 0.02 Pmax
npt 20 ps T 300 to 600 P 1 atm # a temperature ramp
nvt 5 ps T 300 cap 20          # LJ forces capped at 20 kcal/mol/Å
```

- **larsen21:** Larsen, Lin, Hart and Colina, Macromolecules 44, 6944 (2011). Seven heat / cool / compress cycles up
  to Pmax = 5 × 10⁴ bar, then decompression, then 800 ps of NPT at the final conditions: 1.56 ns in total. The high
  temperature is a parameter; choose it about 200–300 K above Tg. `--scale` shortens every stage, and the result is
  labelled as shortened, not as the published schedule.
- **annealing:** NPT ramps between two temperatures with holds.
- **pushoff:** NVT with LJ forces capped at 5 → 500 kcal/mol/Å, then uncapped NVT and NPT.

`--until-converged` continues with NPT blocks at the final conditions. It stops when two successive block-to-block
changes stay within 0.5 % in density, 0.005 kcal/mol per atom in potential energy, and 2 % in mean Rg. These checks
catch drift; they cannot prove that long chains have relaxed, which takes far longer than any MD run.

## React

```bash
build/cli/caps react x --list-templates                                     # built-in templates as text
build/cli/caps react melt.data -o network.data --template cc_crosslink --per-cycle 3 --cycles 20 --md-ps 5 --temp 500
build/cli/caps react mix.data -o epoxy.xyz --template epoxy_amine_primary --template epoxy_amine_secondary --no-relax --fa 2 --fb 4
build/cli/caps react mix.data -o out.data --template my_reaction.txt
```

A template is a small atom-mapped pattern with its edits:

```
reaction epoxy_amine_primary
atom 1 C ring3 H=2          # epoxide CH2
atom 2 O ring3 bonded 1
atom 3 C ring3 bonded 1 2
atom 4 N H=2                # primary amine
atom 5 H bonded 4
initiators 4 1              # the pair whose distance is tested
capture 4.5
form 4 1
break 1 2
move 5 2                    # proton from N to O
sites 1 2 3                 # conversion is counted per epoxide ring
```

Each cycle finds every match with its initiators inside the capture distance. It reacts the closest non-overlapping
matches, retypes and minimises the structure, and can run a short NVT stage. Conversion, clusters, the largest
cluster's mass fraction and the reduced weight-average mass (without the largest cluster) are recorded per cycle. The
gel point is reported where that reduced mass peaks; compare it with the Flory–Stockmayer α_c from `--fa --fb --ratio`.

Measured on the 1 300-atom polystyrene cell (`cc_crosslink`, 3 per cycle, 8 cycles, 24 crosslinks). After a full
minimisation of each result, the network is at +1 168 kcal/mol without dynamics between cycles and +498 kcal/mol with
5 ps at 500 K after each cycle. For reference, the uncrosslinked cell is at −69 kcal/mol, with 48 more hydrogens. In a
glass the chains cannot move to accommodate new bonds by minimisation alone; run dynamics between cycles, as the
Polymatic cycle does.

The C–C template uses the two leaving hydrogens as initiators, so only C–H bonds that point at each other react.

## Force fields

CAPS keeps force fields in its own JSON format (`caps-forcefield`, see `core/include/caps/ffdef.hpp`): atom types with
equivalences, parameter rules matched on type names (the last matching rule wins), bond increments, automatic
(auto-equivalence) parameters and every class II cross term. `data/forcefields/` holds the curated library built by
`bench/ff/build_library.py` from moltemplate and DL_FIELD 4.13; `catalogue.json` lists every source file with its
version / year, primary reference and status (validated, converted, pending, template, alias).

```bash
caps ff import-lt  ~/moltemplate/moltemplate/force_fields/gaff2.lt -o gaff2.json     # moltemplate
caps ff import-dlf ~/dl_f_4.13/lib/PCFF.par -o pcff.json                              # DL_FIELD (.par + .sf + .bci)
caps ff apply molecule.mol2 --ff data/forcefields/gaff-amber25-dlfield.json --list -o molecule.data
```

mol2 is the preferred input: its bonds, bond orders, atom types and charges are used as given (xyz needs bonds guessed
from distances). The evaluator has the class II forms (quartic bonds and angles, bond-bond, bond-angle, middle / end
bond-torsion, angle-torsion, angle-angle-torsion, bond-bond 1-3, Wilson out-of-plane with angle-angle), 9-6 LJ with
sixth-power mixing, and inversion/harmonic. Validation (scripts in `bench/ff/`):

| Check | Result |
|---|---|
| class II terms, random coefficients, vs LAMMPS CLASS2 | energies to 1e-5 kcal/mol, forces 6e-5 kcal/mol/Å |
| COMPASS (moltemplate) ester-ether, propylbenzene, siloxane vs LAMMPS | forces ≤ 7e-8 kcal/mol/Å; all cross-term references on their own bond / angle |
| GAFF2, OPLS-AA 2024 (moltemplate) vs moltemplate + LAMMPS | identical coefficients, forces ≤ 9e-8 kcal/mol/Å |
| PCFF, CVFF, COMPASS, OPLS 2005, GAFF (AmberTools 16 and 25) from DL_FIELD vs DL_FIELD's own LAMMPS output | every interaction identical; forces ≤ 6e-5 kcal/mol/Å, distorted geometries included |
| Morse / GROMOS bonds, cosine angles, Urey–Bradley, planar inversion, Buckingham, Morse pairs, 1-4 LJ vs LAMMPS (`bench/ff/check_forms.cpp`) | forces ≤ 5e-6 kcal/mol/Å |
| DL_FIELD libraries on their own molecule templates vs DL_FIELD's DL_POLY FIELD file (`bench/ff/validate_family.py`) | CGenFF 60/60, DREIDING 29/29, TraPPE-EH 9/9, CHARMM, CHARMM22, TraPPE-UA, GROMOS, CL&P, inorganic shell-model oxides — see `bench/ff/README.md` |

Findings recorded while validating: moltemplate's COMPASS doubles bond-increment charges when a molecule is written with
`Data Bond List`, and its canonical atom sort places up to 55 % of class II cross-term reference values on the wrong bond or
angle; DL_FIELD's PCFF / COMPASS / CVFF contain no class II cross terms (DL_POLY cannot evaluate them); a few DL_FIELD
library entries carry wrong masses (HC_benzyl 12.0115 in PCFF, C_benzene 1.00797 in CVFF, O_thioester 15.0994), which CAPS
reports and does not use. DL_FIELD's LAMMPS export leaves out CHARMM 1-4 van der Waals and DREIDING inversions and cannot write GROMOS, so those
families are checked against its DL_POLY FIELD file; its GROMOS bonds are a harmonic approximation (CAPS keeps the
quartic form); its TraPPE-UA C-C-O-H torsion takes terms from the O-C-C-O entry (CAPS uses the published alcohol
torsion). Pending: full PCFF / CVFF with class II cross terms from Accelrys `.frc` files.

### Automatic atom typing

Atoms get their force-field types from SMARTS rules (`core/include/caps/typing.hpp`, rule files in `data/typing/`).
CAPS first perceives the chemistry from the connectivity alone, so mol2, PDB and xyz all work: bond orders from
valences (a Kekulé structure found by search), formal charges (carboxylates, nitro groups, ammonium, iminium /
guanidinium / imidazolium cations), rings (SSSR) and Hückel aromaticity. Bond orders given in a mol2 file are kept.
Among the rules that match an atom, types another match `overrides` drop out, then the highest `priority` wins; rules
may refer to other types (`[%oh]`, as in foyer), and typing repeats until nothing changes.

```bash
caps ff type molecule.pdb --ff data/forcefields/pcff-dlfield.json --explain       # each atom's type and the rule behind it
caps ff apply molecule.xyz --ff data/forcefields/pcff-dlfield.json                  # types automatically, then parameterises
caps ff apply water.mol2 --ff data/forcefields/cvff-dlfield.json \
    --typing data/typing/cvff-dlfield.typing.json,data/typing/cvff-dlfield-tip3p.typing.json   # rules on top: TIP3P water
```

`--types FILE` still sets types by hand (all atoms, or `index type` lines for a few). A force field names its rules
file with `"typing"`: PCFF, CVFF, GAFF (AmberTools 16), GAFF2 (AmberTools 25), OPLS-AA 2005 and CGenFF from DL_FIELD
have rules (`data/typing/`; the GAFF, OPLS and CGenFF files are generated by `bench/typing/make_*_rules.py`, one
commented line per rule). Two additions beyond SMARTS serve GAFF and CGenFF: `{AR1}`..`{AR5}`, antechamber's ring
classes, and conjugated type pairs (`"pairs"`: GAFF's cc/cd, ce/cf, nc/nd ... alternate across double bonds; CGenFF's
CG2DC1/CG2DC2 with `"pair_mode": "double_same"`). A rules file marked `"ordered"` takes the first matching rule, as
antechamber does.

| Check | Result |
|---|---|
| PCFF: 92 DL_FIELD molecule templates, typed from connectivity alone (`bench/ff/validate_typing.py`) | 1387 / 1400 atoms (99.1 %) as DL_FIELD's templates; the rest are inconsistencies between DL_FIELD's own templates (furan, oxazole and indole ring carbons `cp`, pyrrole's `c5`) |
| CVFF: 23 DL_FIELD templates | 272 / 278 atoms (97.8 %); the rest are water, whose model (TIP3P / SPC) is a choice |
| CVFF: Materials Studio's typing of the msi2lmp examples — crambin, nylon, aromatics (`bench/ff/compare_msi_types.py`) | 811 / 823 atoms (98.5 %); crambin 642 / 642 |
| end to end: automatically typed templates through `ff apply` | every parameter found except where DL_FIELD's own library has none (nitroso N=O bonds, charged-imidazole N–H), the molecules DL_FIELD also refuses |
| CGenFF: the 494 model compounds of CGenFF 3.0.1's `top_all36_cgenff.rtf`, typed by its developers (`bench/ff/validate_rtf_types.py`) | 9024 / 9174 atoms (98.4 %) from bonds alone, 417 residues fully right; 97.7 % with the RTF's partial bond orders; 420 residues fully parameterised by DL_FIELD's CGenFF library |
| GAFF / GAFF2: rules after antechamber's ATOMTYPE_GFF(2).DEF, on DL_FIELD's templates | 97.9 % / 96.2 % of atoms; the rest are templates departing from antechamber (H on C–O / C–S as hc, azobenzene N) |
| OPLS-AA 2005 (DL_FIELD names), 99 templates | 1368 / 1398 atoms (97.9 %), 91 molecules fully right |
| 6440-atom polystyrene melt, GAFF2 | typed in 0.07 s, every atom as CAPS's GAFF builder typed it |

In the Studio, the Field page does the same interactively (`caps_field_*` in the C API, ABI 7): the chosen force field
is written into the structure, Relax, Dynamics, Equilibrate and LAMMPS data then use it, and while atoms are untyped or
parameters missing those runs are refused. Parameters can be imported (CAPS JSON, moltemplate `.lt`) or entered by hand;
hand-entered terms are counted as estimated in the report.

### LAMMPS data

`caps ff apply FILE --ff FF.json -o out.data --lammps-input out.in` (and Save in the Studio) writes every term CAPS
evaluates in the LAMMPS style with the same energy: harmonic, class II, Morse and GROMOS bonds; harmonic, class II,
CHARMM (Urey–Bradley), cosine and cosine/squared angles; Fourier and class II torsions; cvff, harmonic, class II,
inversion/harmonic and umbrella impropers; LJ 12-6 or 9-6, Buckingham and Morse pairs. A kind that mixes forms becomes a
hybrid style (DL_FIELD's PCFF: class II bonds and angles with Fourier torsions), and the class II cross-term sections get
`skip` lines for the other sub-styles' types. Pair coefficients are written for every i-j pair, so LAMMPS does no mixing.
The data file's header and the input script hold the matching LAMMPS commands. Separate 1-4 Lennard-Jones parameters
(CHARMM, GROMOS) have no exact LAMMPS form without switching; the writer refuses them rather than approximate.

| Check (`bench/ff/check_data_lammps.py`, LAMMPS `run 0` on the written files) | Result |
|---|---|
| PCFF and COMPASS from DL_FIELD (class II + Fourier, hybrid), CVFF, OPLS-AA, GAFF, GAFF2, DREIDING (umbrella), ionic crystals (Buckingham, periodic), a periodic polystyrene melt with GAFF2 and with PCFF, COMPASS with a class I overlay (hybrid in every kind, skip lines in every class II section) | 16 of 16: every energy term to ≤ 6 × 10⁻⁷ (relative), every force to ≤ 1.4 × 10⁻⁶ kcal/mol/Å |
| UFF: a mixed molecule set (P, S, Si, Pt, F, Cl; Fourier and cosine/periodic angles, umbrella and fourier impropers) and a periodic polystyrene melt | 2 of 2: energy terms ≤ 2 × 10⁻⁷ (relative), forces ≤ 2.8 × 10⁻⁶ kcal/mol/Å |
| CGenFF (separate 1-4 LJ) | refused with the reason |

Found on the way: CAPS's damped-shifted-force self energy left out the force-shift part of the shift constant that
LAMMPS `coul/dsf` includes (the r → 0 limit of the same pair potential). A constant, so forces were already equal;
Coulomb energies now equal LAMMPS's too (2.2 kcal/mol for favipiravir with PCFF's bond-increment charges).

UFF is built in rather than kept as JSON: `--ff uff` (CLI), "UFF · every element" (molecule builder) and "UFF (Rappé
1992) · every element" (Field) type any structure from its elements and bonds. The atomic parameters are the UFF table
as RDKit distributes it (BSD licence, `licenses/RDKit-UFF-parameters.txt`); CAPS writes its own typer and terms. UFF
runs without charges by default, as a clean-up force field; `--charges keep` uses the file's.

GAFF and DL_FIELD's OPLS carry no charges on their types: use `--charges gasteiger` (or charges from the file). Water
defaults to each force field's TIP3P (or CVFF's own); other water models are one-line overlay rules files.

## Analyze

`caps analyze TRAJ --topology DATA --props density,rdf,sq,rg,msd,diffusion,ced,ffv,... --frame-ps 0.5 [--json out.json] [--csv DIR]`
(and Analyze › Properties in the Studio) computes properties over the frames, each with its method, a block-average
standard error where it applies, notes when the input cannot support the number (short chains, sub-diffusive MSD, an
end-to-end vector that has not decorrelated) and the curves behind it.

| Group | Properties | Method |
|---|---|---|
| Structure | density, rdf, sq, xray, neutron | g(r) with a cell list; S(q) (Faber–Ziman, number, Cromer–Mann or neutron weights) as the exact sum over the cell's reciprocal lattice up to 4 Å⁻¹ (shells merged to ≥ 24 k-vectors), the g(r) transform with a Lorch window above |
| Chains | rg, ree, cn, persistence | mass-weighted Rg; backbone end-to-end; C_n = ⟨R²(n)⟩/(n⟨b²⟩) and C∞ extrapolated in 1/n; Flory projection and bond-correlation decay |
| Dynamics | msd, diffusion, relaxation | all time origins, system drift removed; D from the molecule-centre MSD over a chosen window, with the log-log slope; P2 bond and end-to-end autocorrelations with KWW fits |
| Thermo | ced, delta | (E isolated − E bulk)/V with the assigned force field (Field; GAFF of C and H otherwise), δ = √CED |
| Free volume | ffv, psd | probe insertion on a grid with Bondi radii (accessible fraction for several probes, Bondi FFV); Gelb–Gubbins pore sizes, largest cavity refined off the grid |

| Check (`bench/analyze/check_analyze.py` on a 100 ps polystyrene run, 201 frames) | Result |
|---|---|
| g(r) C–C against MDAnalysis `InterRDF` | largest difference 2.7 × 10⁻⁴ |
| √⟨Rg²⟩ against MDAnalysis `radius_of_gyration` | equal to 10⁻⁴ Å |
| atom and molecule-centre MSD, D against numpy (drift removed) and MDAnalysis `EinsteinMSD` | equal to 10⁻⁴ (relative) |
| S(q) against a numpy sum over all 5125 k-vectors of the cell | 6 × 10⁻⁷ |
| free fraction against Monte Carlo insertion (400 000 points) | 0.3976 vs 0.3964 ± 0.0023 |
| largest pore of a simple cubic lattice against 2(a√3/2 − r) | equal to 10⁻⁴ Å |

Unit tests add a Bragg peak of a simple cubic crystal, an ideal gas (S = 1), a single sphere's volume and D of a
random walk with drift. Not built yet: entanglements (primitive-path analysis).

### Visualize pipeline

`caps pipeline FILE [--topology DATA] --steps STEPS.json [--frame N] [--table NAME] [--particles EXPR]` runs the Studio's
Analyze › Visualize steps on one frame and prints each step's result, the global attributes and, with `--table`, a data
table as CSV. Steps are listed top first and run bottom to top; expressions use particle properties
(`Type == 2 && Position.Z > 13`, `Element == "O"`, `sqrt(Position.X^2 + Position.Y^2) < 10`).

```bash
caps pipeline samples/ps_melt.lammpstrj --topology samples/ps_melt.data \
  --steps '[{"type":"coordination","cutoff":1.25,"element_a":6,"element_b":1},{"type":"cluster"}]' --table clusters
```

### Python steps

A pipeline step can be a Python script: a function decorated with `@step` from the small `caps` package in
data/python reads `data.particles` (identifier, molecule, type, element, charge, selection, backbone flag, whole-molecule
positions and every property earlier steps made) and writes `data.attributes`, `data.tables`, new particle properties
and a selection, which the next steps and the data inspector see. CAPS runs it in a separate Python process
(`$CAPS_PYTHON`, else python3); numpy is used when the script imports it.

```python
import numpy as np
import caps
from caps.pipeline import step

@step(name="Backbone conformation")
def modify(frame, data):
    mol, bb = data.particles["Molecule Identifier"], data.particles["Backbone"]
    pos = data.particles.positions_unwrapped
    trans = []
    for m in np.unique(mol):
        idx = np.flatnonzero((mol == m) & (bb == 1))
        trans += list(np.abs(caps.geometry.dihedrals(pos[idx])) > 120)
    data.attributes["TransFraction"] = np.mean(trans)
```

### Figure bundles

`caps bundle FILE --steps STEPS.json --include-input -o NAME.caps-bundle.zip` writes a figure with the data behind it:
figure.png and figure.svg, data/*.csv for every table the pipeline makes, pipeline.json, the input structure,
provenance.json (CAPS version, the input's sha256, the pipeline's sha256 and the sha256 of every file) and a README.
`caps reproduce NAME.caps-bundle.zip` rebuilds the data files from the bundled input and pipeline and compares their
hashes with the record. The Studio writes the same bundles from Export › Figure bundle.

### Mechanics and the glass transition

| Command (and Analyze chip) | Method |
|---|---|
| `caps elastic FILE --method strain` (Cij strain) | Theodorou–Suter static constants: minimise at fixed cell, ±ε in each Voigt direction (pure strain), re-minimise, C_IJ = Δσ_I/Δε_J; averaged over `--configs` frames. The pair set is frozen at the minimum, so the truncated-LJ surface is smooth and the minimiser converges (L-BFGS, 10⁻⁴ kcal/mol/Å) |
| `caps elastic TRAJ --method fluct --temp T` (Cij fluct.) | stress fluctuations of an NVT trajectory (Lutsko 1989): ⟨Born⟩ − (V/kT) cov(σ) + NkT/V, block errors; the Born term by central differences of the second Piola–Kirchhoff virial stress under Lagrangian strain |
| `caps tensile DATA --axis x --rate 1e-3 --strain 0.2` (Stress–strain) | uniaxial deformation at constant engineering rate (as LAMMPS `fix deform erate`, remap x), lateral faces at 1 atm by per-axis Berendsen coupling or fixed; modulus and Poisson ratio from 0–2 %, 0.2 % offset yield, peak |
| `caps tg DATA --from 500 --to 200 --step 20 --ps 100` (Tg), `caps tg --fit TABLE.csv` | stepwise NPT cooling, density averaged over the second half of each hold, Tg = hinge of a continuous two-line fit of specific volume (bootstrap error), expansion coefficients above and below |

Every result reports Voigt, Reuss and Hill averages (E, K, G, ν, λ), the pre-stress of the cell, and notes (asymmetry,
unconverged minimisations, strain rates far above experiment, cooling rate). The evaluator now carries the full virial
tensor, and Dynamics reports the pressure tensor and supports per-axis pressure coupling and constant-rate deformation.

| Check | Result |
|---|---|
| virial tensor, every term type (bonds, angles, torsions, impropers, class II, Morse, Buckingham, UB, inversions, 1-4 LJ) against −∂E/∂ε | unit tests, all six components |
| virial tensor against LAMMPS `compute pressure NULL virial` (`bench/ff/check_data_lammps.py`, 16 force-field cases) | ≤ 1.6 × 10⁻⁷ relative |
| Born term, bonded (`bench/mechanics/check_born_lammps.py`) against LAMMPS `compute born/matrix numdiff` | 2 × 10⁻¹⁰ |
| Born term, bonded + LJ, against LAMMPS numdiff (step 10⁻⁷) | 2 × 10⁻⁹ |
| Born term of an fcc argon crystal against the lattice sum Σ(φ″ − φ′/r)X⁴/r² | 10⁻⁴ (unit test); LAMMPS numdiff also equals it |
| static constants of a Bravais lattice against Born + the Cauchy-stress terms | unit test |
| fluctuation constants of a crystal at 5 K against its static constants | within 6–10 % (unit test) |
| tensile modulus and Poisson ratio of a cold crystal against (C11 − C12)(C11 + 2C12)/(C11 + C12) and C12/(C11 + C12) | within 12 % / 0.1 |
| two-line fit of noisy synthetic data | Tg within 8 K |

Found on the way: LAMMPS's analytic `compute born/matrix` pair term (this build) disagrees with its own `numdiff`
mode and with the lattice sum: it sums the prefactor ½φ″ − φ′/r over a full neighbour list, which counts φ′/r twice
(argon: C11 3260 against 3400 kcal/mol). And DL_FIELD's library files carry element and mass typos in single atom-type
entries (AMBER25 `HC_alkyne → ha` as C 12.0115, CVFF `O_alcohol → oh` at 12.0115, CGenFF `NG2O1` at 12.007, TraPPE-EH
`H_3` at 12.0115, …); the importer now takes the element by majority over all entries of a type and replaces a mass
equal to another element's standard mass, and the library was rebuilt (masses only; energies and forces unchanged).

## Samples

`samples/ps_melt.*` is one built cell: 10 atactic polystyrene chains of 8 units with hydrogens in a 33 Å box
(1300 atoms, 1370 bonds, density 0.386 g/cm³ — loosely packed, not equilibrated). The `.lammpstrj` has three
frames in which frame *k* is frame 0 shifted by 0.5·*k* Å in x: a reader test pattern, not dynamics.

## Layout

```
core/       C++20 library: model, readers/writers, bonds, analysis, Grow, Pack, Field, Relax, Dynamics, Equilibrate, React, renderer
capi/       C ABI (libcaps) used by the Studio through P/Invoke
cli/        caps command-line tool
tests/      GoogleTest suite
studio/     Avalonia desktop app (C#, MVVM)
samples/    small input files used by tests and the Studio
design/     exported screen designs (boards, canvas, generator scripts)
scripts/    build, launch and LAMMPS cross-check helpers
```
