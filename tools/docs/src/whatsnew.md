---
title: What's new
eyebrow: Releases · changelog
lede: What changed in CAPS, newest first. Every force field CAPS ships is now built from its original paper or official distribution, CAPS can fill missing parameters with its own machine-learned model, and the Studio and the command line do the same things the same way.
---

## October 2026 update (version 0.1.0)

### Semicrystalline stacks and thermal conductivity

- `caps semicrystalline` builds the lamellar cells the papers on semicrystalline polyethylene simulate: a lamella of straight stems along z from the polymer crystal library (orthorhombic PE from Bunn 1939), and an amorphous layer between it and its periodic image holding loops, bridges and tails in the fractions you ask for, at the amorphous density. Every stem runs through the lamella, every chain ends in two tails, and none closes on itself. The cell carries the groups `crystal` and `amorphous`. A report gives the counts, the length of each segment, the crystallinity from the chain order parameter P2, the density and P2 profiles along the stack and the stem tilt; `--analyse` reads the same from any frame after MD. `--lammps-input` writes the cluster equilibration with the crystal tethered first. A united-atom force field (`--ff oplsua2024`) builds carbons only, for thicker layers at the same cost.
- Stems are found in any cell of the polymer crystal library, straight or helical (iPP, PEO, POM, PVDF, the nylons, PET, PLLA, natural rubber). The amorphous layers are grown for polyethylene so far.
- `caps thermal` computes the thermal conductivity by MD: Green–Kubo from the heat-flux autocorrelation, with the centroid form of the flux for angles, torsions and impropers, or reverse NEMD (Müller-Plathe). You say whether the material is 3D (bulk) or 2D (a sheet with vacuum along z; κ over an effective thickness and the sheet conductance); CAPS does not guess. On one frame the flux equals LAMMPS's `compute heat/flux` on `compute centroid/stress/atom`. On the Lennard-Jones liquid both methods agree with the published values. `--lammps-input` writes the long run for a cluster.
- In the Studio: **Build › Semicrystalline stack** and **Modules › Thermal conductivity**; in Python `build.semicrystalline()` and `Document.thermal()`; in recipes `build: { semicrystalline: … }` and a `thermal:` stage.

### Polymer crystals for DFT

- A library of polymer unit cells from the original structure papers: orthorhombic and monoclinic PE, iPP α, PEO, PVDF β and α, POM, nylon 6 and nylon 66, PET, PLLA, natural rubber, and PTFE in its room-temperature (15/7 helix) and low-temperature (13/6 helix) forms. Each cell carries its paper's lattice, density and page references. The PTFE chains are generated from the helix symmetry in the cells of Clark & Muus (1962), and the files say so.
- `caps vasp-relax` writes a four-step VASP relaxation chain (cutoff and cell steps, PBE with D3-BJ) with one batch script that stops when a step does not converge.
- `caps phono3py setup|estimate|check|analyze` prepares the lattice thermal conductivity from third-order force constants: supercells from the cell heights, full or random displacements, a dry run that sets the array time limits, a cost estimate before anything is submitted, and the analysis of κ along and across the chains. `--dimensionality 2d|3d` handles sheets with vacuum.

### Figures you can put in a paper

- Every plot in the Studio is drawn as a journal figure: a frame with inward ticks, round tick values, axis titles as quantity / unit, real log axes, legends that keep clear of the data, error bars and block averages where CAPS computes them, and fits drawn with their window. Tg shows the glass and melt fits with Tg marked.
- The DPD page shows the final domains in 3D; the adsorption, sorption and DPD pages report from real runs.
- The pictures in this manual and on the module pages are taken from real CAPS runs.

### Force fields from original sources only

Every force field in CAPS now comes from its original paper, its supporting information, or the distribution its authors publish. Each library file says where its numbers came from.

- Water from the original papers and their supporting information: SPC and SPC/E (Berendsen 1987), TIP3P and TIP4P (Jorgensen 1983), TIP3P for Ewald sums (Price & Brooks 2004), the CHARMM TIP3P (MacKerell 1998), TIP4P-Ew (Horn 2004), TIP4P/2005 and TIP4P/Ice (Abascal & Vega 2005), TIP5P (Mahoney & Jorgensen 2000), OPC and OPC3 (Izadi & Onufriev), the flexible SPC/Fw (Wu 2006) and mW (Molinero & Moore 2009). Each reproduces its paper's geometry and charges exactly, and its LAMMPS and GROMACS files give the same energy as CAPS.
- Gases: CO2 as EPM2 (Harris & Yung 1995) besides TraPPE, H2 (Buch 1994; Darkrim & Levesque 1998) and the noble gases (Rutkai et al. 2017). N2, O2, CH4, NH3 and H2S come from the TraPPE build.
- L-OPLS (Siu, Pluhackova & Böckmann 2012), now as an extension of OPLS-AA 2024.
- Coarse-grained: Martini 2 (Marrink 2007, the paper's interaction matrix, with the lipid angles of Marrink 2004), its polymer set for PEO/PEG (Lee et al. 2009) and polystyrene (Rossi et al. 2011), Dry Martini (Arnarez 2015), SDK (Shinoda 2007 and 2010) and the Cooke–Deserno lipid model, alongside Martini 3. All-atom PEO and polystyrene map onto the Martini 2 beads.
- Short names still work: `--ff tip4p2005`, `--ff spce`, `--ff opc`, `--ff martini` and so on resolve to these libraries.

### OPLS-AA 2024 and OPLS-UA 2024

- Sulfur-cured rubber: persulfide and polysulfide types, the S–S bonds, bends and torsions, and the allylic torsions next to the C=C, fitted to MP2/cc-pVTZ//HF/6-31G** scans. Natural rubber and polybutadiene with S1 to S8 bridges are complete in both force fields.
- OPLS-UA 2024 gains the sulfur torsions of Jorgensen's 1986 paper (Table III) and united-atom refits of the all-atom scans, and now builds on CAPS's own OPLS-AA 2024 file.
- Alkoxysilanes (silane coupling agents such as TESPT and TESPD, TEOS, methoxy and dialkoxy silanes) are complete: new types with charges from 1.27 × CM5 that keep any mix of alkoxy and alkyl groups on Si neutral, and Si–O torsions fitted to MP2 scans. TESPT and TESPD crosslink rubber in React with every term in place.
- Maleic anhydride, maleates and cured DGEBA are complete.
- Optional overlays: DLPNO-CCSD(T) torsions (`oplsaa2024-dlpno`), a Lennard-Jones refit from liquid data (`oplsaa2024-lj-refit`) and a rubber torsion refit. `caps ff overlays` lists them.
- All 135 polymers and copolymers of the library build, type, relax and run in LAMMPS and GROMACS with OPLS-AA 2024. The styrene–butadiene and acrylonitrile–butadiene junctions of SBR, NBR and ABS have new torsions fitted to MP2 scans.
- Known limit: glassy PMMA, PS, PVAc and PET come out 3–8 % below their measured densities with OPLS-AA 2024, even after a slow anneal. The note `bench/ff/open/oplsaa2024_glass_density_limitation.md` gives the evidence.

### CAPS ML: missing parameters filled by CAPS's own model

CAPS now carries its own machine-learned potential (CAPS MACE v1, trained on wB97M-D3(BJ)/def2-TZVPPD data, about 0.15 kcal/mol median error on torsion profiles of molecules it never saw). When a force field lacks a torsion, bend or bond for your structure, CAPS can fit it to relaxed scans with this model.

- `caps ff complete FILE --ff FF.json -o OVERLAY.json` fits what is missing and saves it as an overlay.
- `--complete ml` (on `caps ff apply`, `relax`, `md`, `equilibrate`, `react`, `tg` and `tensile`), `complete: ml` in a recipe, or the Studio's "Fill missing terms with CAPS ML" switch fills the gaps whenever a force field is assigned.
- Every parameter CAPS ML supplies is marked where you will see it: a `# CAPS ML ...` comment on its line in the LAMMPS data file, a `; CAPS ML ...` comment in the GROMACS `.itp`, a summary at the top of each file, and an entry in the provenance file.
- Charges and atom types are never guessed. A term the model cannot fit within its quality limits stays reported as missing, with the reason.

### Studio and command line in step

- React's options (crosslinks between chains, crosslink targets, byproducts, weights, automatic capture, sites per chain, removal of unreacted sulfur donors, the sulfur bridge census) work the same in the Studio, the command line and Python, through the same engine code.
- "Copy as command" on the Relax, Dynamics, Equilibrate, React, glass transition, tensile and force-field pages gives the `caps` command that does the same thing.
- The full-screen viewer draws on the GPU, as the main view does: about 2.5 ms per frame on a 14,000-atom blend, where it took 69 ms before.
- A clean install no longer offers force fields it does not ship.
- CAPS Studio tells you when a newer CAPS is out: once a day it reads the latest release's `release.json` from GitHub (nothing about you is sent) and shows a banner with What's new, Download, Later and Skip this update. After an update it shows what is new, with a button to each page. Settings › Files has the switch, Help has Check for updates now, `CAPS_NO_UPDATE_CHECK=1` turns it off, and `caps --check-update` does the same from the command line.

### Faster runs

- `caps tg` and `caps tensile` use the mixed-precision kernel of `caps md`: about twice as fast, with the same results.
- Chain growth is 10 to 60 times faster for every growth method (best-of-k trials, Rosenbluth with soft spheres or Lennard-Jones, look-ahead, stars, combs, slabs, oriented cells): 24 polyethylene chains of 100 units (14,448 atoms) grow in about 1 s instead of 64 s, and 15 polystyrene chains of 60 units in 3 s instead of 87 s. The same seed still gives the same cell.
- The grower copes with any mix of chain count and length: a chain that gets stuck relaxes its own contact limits a step at a time instead of restarting the cell, end caps are kept clear while the cell grows, and bulky ring units such as ethylidene norbornene grow without stalling. Across the polymer library at 0.3 and 1.1 g/cm³ (1 × 1 up to 200 × 2 and 50 × 5, up to 14,000 atoms), 2,955 of 3,240 cells now build and pass every structure check, up from 1,649. A cell denser than chain growth can fill now stops early and says to grow at a lower density and compress.

### Files CAPS writes

- LAMMPS coefficients are written with six decimals. GROMACS files use eight decimals for lengths in nm and energies in kJ/mol, because six decimals in nm is ten times coarser than in Å; charges keep six decimals and stay exactly neutral.
- GROMACS `.mdp` files use a fixed pair-list buffer, so a plain Lennard-Jones cut-off with tail corrections gives the same density as LAMMPS.
- `caps render --wrap` draws atoms folded into the cell.
