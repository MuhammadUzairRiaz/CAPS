---
title: What's new
eyebrow: Releases · changelog
lede: What changed in CAPS, newest first. Every force field CAPS ships is now built from its original paper or official distribution, CAPS can fill missing parameters with its own machine-learned model, and the Studio and the command line do the same things the same way.
---

## October 2026 update (version 0.1.0)

### Force fields from original sources only

Every force field in CAPS now comes from its original paper, its supporting information, or the distribution its authors publish. Each library file says where its numbers came from.

- Removed: every force field CAPS had converted from another program's bundled parameter files, including OPLS 2005, OPLS-AA 2008, L-OPLS 2008, the second copies of OPLS-AA 2024, GAFF, GAFF2, DREIDING, TraPPE and COMPASS, and the earlier Martini 2, Dry Martini, SDK, SPC/E, TIP3P, mW, graphene and Cooke–Deserno files. The ones in use are rebuilt below; `.lt` input and output still work.
- Water, rebuilt from the original papers and their supporting information: SPC and SPC/E (Berendsen 1987), TIP3P and TIP4P (Jorgensen 1983), TIP3P for Ewald sums (Price & Brooks 2004), the CHARMM TIP3P (MacKerell 1998), TIP4P-Ew (Horn 2004), TIP4P/2005 and TIP4P/Ice (Abascal & Vega 2005), TIP5P (Mahoney & Jorgensen 2000), OPC and OPC3 (Izadi & Onufriev), the flexible SPC/Fw (Wu 2006) and mW (Molinero & Moore 2009). Each reproduces its paper's geometry and charges exactly, and its LAMMPS and GROMACS files give the same energy as CAPS.
- Gases: CO2 as EPM2 (Harris & Yung 1995) besides TraPPE, H2 (Buch 1994; Darkrim & Levesque 1998) and the noble gases (Rutkai et al. 2017). N2, O2, CH4, NH3 and H2S come from the TraPPE build, as before.
- L-OPLS (Siu, Pluhackova & Böckmann 2012), now as an extension of OPLS-AA 2024.
- Coarse-grained: Martini 2 (Marrink 2007, the paper's interaction matrix), Dry Martini (Arnarez 2015), SDK (Shinoda 2007 and 2010) and the Cooke–Deserno lipid model, alongside Martini 3.
- Not yet rebuilt: the Martini 2 amino-acid, polymer and sugar sets, and graphene carbon (its 1973 source paper is still to be obtained).
- Short names still work: `--ff tip4p2005`, `--ff spce`, `--ff opc`, `--ff martini` and so on resolve to the rebuilt libraries.

### OPLS-AA 2024 and OPLS-UA 2024

- Sulfur-cured rubber: persulfide and polysulfide types, the S–S bonds, bends and torsions, and the allylic torsions next to the C=C, fitted to MP2/cc-pVTZ//HF/6-31G** scans. Natural rubber and polybutadiene with S1 to S8 bridges are complete in both force fields.
- OPLS-UA 2024 gains the sulfur torsions of Jorgensen's 1986 paper (Table III) and united-atom refits of the all-atom scans, and now builds on CAPS's own OPLS-AA 2024 file.
- Alkoxysilanes (silane coupling agents such as TESPT and TESPD, TEOS, methoxy and dialkoxy silanes) are typed, with charges from 1.27 × CM5 that keep any mix of alkoxy and alkyl groups on Si neutral.
- Maleic anhydride, maleates and cured DGEBA are complete.
- Optional overlays: DLPNO-CCSD(T) torsions (`oplsaa2024-dlpno`), a Lennard-Jones refit from liquid data (`oplsaa2024-lj-refit`) and a rubber torsion refit. `caps ff overlays` lists them.
- The 135 polymers and copolymers of the library build, type, relax and run in LAMMPS and GROMACS with OPLS-AA 2024 (131 without help; the other four are copolymer junctions that CAPS ML completes).
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

### Faster runs

- `caps tg` and `caps tensile` use the mixed-precision kernel of `caps md`: about twice as fast, with the same results.
- Chain growth is 10 to 60 times faster for every growth method (best-of-k trials, Rosenbluth with soft spheres or Lennard-Jones, look-ahead, stars, combs, slabs, oriented cells): 24 polyethylene chains of 100 units (14,448 atoms) grow in about 1 s instead of 64 s, and 15 polystyrene chains of 60 units in 3 s instead of 87 s. The same seed still gives the same cell.
- The grower copes with any mix of chain count and length: a chain that gets stuck relaxes its own contact limits a step at a time instead of restarting the cell, end caps are kept clear while the cell grows, and bulky ring units such as ethylidene norbornene grow without stalling. Across the polymer library at 0.3 and 1.1 g/cm³ (1 × 1 up to 200 × 2 and 50 × 5, up to 14,000 atoms), 2,955 of 3,240 cells now build and pass every structure check, up from 1,649. A cell denser than chain growth can fill now stops early and says to grow at a lower density and compress.

### Files CAPS writes

- LAMMPS coefficients are written with six decimals. GROMACS files use eight decimals for lengths in nm and energies in kJ/mol, because six decimals in nm is ten times coarser than in Å; charges keep six decimals and stay exactly neutral.
- GROMACS `.mdp` files use a fixed pair-list buffer, so a plain Lennard-Jones cut-off with tail corrections gives the same density as LAMMPS.
- `caps render --wrap` draws atoms folded into the cell.
