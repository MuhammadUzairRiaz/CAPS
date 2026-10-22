# Vitrimer and dynamic covalent network simulations: per-paper notes

These are the full extraction notes behind `literature_review.md`, one section per paper, for all 38 PDFs in
`~/project/ff-sources/papers/vitrimers/`. They were made on 2026-10-22 from text extracted from each PDF (PyMuPDF), so
equations that were images in the PDF are sometimes missing. Where that happened, the note says so. Every DOI was looked
up on Crossref, or on DataCite for arXiv DOIs, and title, authors, journal, volume, pages and year were compared. Numbers
are copied as printed. "[unverified]" marks anything that could not be checked, "[inference]" or "[derived]" marks the
note-taker's own arithmetic or reasoning, and "not stated" means the paper does not give the value. Where a local file is
a preprint, the published version was not compared line by line.

Groups A to F are the order in which the papers were read, not a ranking.

## Group A

---

### 1. Adhikari et al. 2026 (ChemRxiv), COLVARS-biased ReaxFF of DGEBA / oxalic acid

**Local file:** Adhikari2026_chemrxiv_COLVARS_ReaxFF.pdf

**Citation:** Y. Adhikari, U. Chapagai, S. Khanal, "COLVARS-Driven Reactive Molecular Dynamics of Vitrimer Polymerization and Thermomechanical Property Analysis", ChemRxiv preprint, posted 3 January 2026 (CC-BY 4.0, not peer reviewed). Affiliation: Reformix AI, Kathmandu.

**DOI:** 10.26434/chemrxiv-2026-k6fn6. Crossref: found, type posted-content, title and all three authors match, issued 2026-01-03. Published version: a Crossref bibliographic query on the title returned only the preprint; no journal version found as of 2026-10-22, so it is unpublished.

**Chemistry:** epoxy-acid ring opening (esterification) of DGEBA by oxalic acid, giving beta-hydroxy ester links that are the substrate for later transesterification. Oxalic acid is used as a model diacid "owing to its simple structure", not because it is used in real vitrimers. Catalyst: TBD is mentioned as the experimental catalyst but is NOT included; only the curing step is simulated. The exchange reaction itself (transesterification) is never simulated. Rubber relevance: none (glassy epoxy).

**Model level / force field:** all-atom ReaxFF, parameter set "CHNO-2017_weak" as modified by Vashisth et al. (ref 23, Vashisth 2018 JPCA), charges by fix qeq/reax every step. Reactant geometries pre-optimised by DFT B3LYP/def2-SVP in ORCA, gas phase; the dimer from stage 1 was re-optimised at the same level before replication.

**Exchange / reaction method:** biased MD with the Colvars module in LAMMPS (steered moving harmonic restraints, not metadynamics). Three distance CVs, each with a harmonic restraint V = 1/2 k [xi - xi0(t)]^2 whose centre is ramped linearly from initial to target over N_steps, plus harmonic walls active only outside bounds (Table 1, stage 1):

| CV | atom pair | initial centre (Å) | target (Å) | k (kcal/mol·Å²) | lower wall (Å) / k | upper wall (Å) / k |
|---|---|---|---|---|---|---|
| dist_CO_epoxide | O(epoxide)-C(epoxide) | 1.399 | 2.40 | 80.0 | 1.0 / 5.0 | 4.5 / 5.0 |
| dist_CO_acid | C(epoxide)-O(acid) | 5.459 | 1.43 | 800.0 | 1.0 / 5.0 | 6.5 / 5.0 |
| d_H_diff | H(acid)-O(epoxide) | 5.761 | 1.1 | 500.0 | -8.0 / 2.0 | 8.0 / 2.0 |

Ramp lengths (number of MD steps) are listed as a Table 1 column heading, but the values did not survive extraction [unverified]; the paper says the biasing duration was increased at later stages. Barrier / activation energy: not computed. Temperature 423 K, NVT Nose-Hoover, damping 100 fs, timestep 0.05 fs, starting density 0.5-0.7 g/cm³. Mechanism reported: SN2-like attack, then proton transfer. Transition-state snapshot: epoxide C-O 2.10 Å, forming C-O(acid) 1.77 Å; product C-O(ester) 1.31 Å, C-O(beta-hydroxy) 1.42 Å, O-H 1.01 Å (measured in Avogadro). Associative/dissociative: not applicable (curing only). Swap frequency, acceptance rule: not applicable.

**System size and run length:** hierarchical doubling 57 -> 114 -> 228 -> 456 atoms (final chain = 8 DGEBA + 8 oxalic acid, 1:1 stoichiometry). Bulk model: 8 such chains in a box (so 3648 atoms; my arithmetic, not stated), box edge labelled 36.22 Å in Fig. 11. Annealing: 300 -> 600 -> 300 K in 100 K steps, NVT ramp then NPT at each T, each 12.5 ps at 0.25 fs, Berendsen thermostat 100 fs and barostat 500 fs; two annealing cycles. Crosslink density: not stated.

**Software:** LAMMPS (version not stated) with Colvars and fix qeq/reax; ORCA (B3LYP/def2-SVP); pymbar timeseries; Avogadro.

**Properties and exact analysis:**
- Density: final 300 K NPT segment of anneal; pymbar detects equilibration onset and statistical inefficiency; SEM = sigma / sqrt(N_eff).
- Tg: NPT at 300 K, heat to 450 K, cool back to 300 K in 10 K steps; density vs T fitted with a cubic smoothing spline; Tg where the second derivative rho''(T) crosses zero (inflection). Details "in SI" (not in this file).
- Young's modulus: uniaxial tension along x, y, z separately, 20 % engineering strain in 100 ps (2 × 10^9 s^-1; the text also prints "2×10-9", a typo), 300 K, 0.25 fs, Nose-Hoover barostat on transverse directions; rolling average of stress; linear fit over 3-6 % strain.
- Poisson's ratio: linear fits of transverse vs axial strain, negative slopes of both transverse directions averaged.
- Not computed: stress relaxation, Arrhenius, Tv, creep, MSD, healing, topology.

**Key results:** density 1.2975 ± 0.0009 g·cm⁻³ at 300 K (compared with 1.215 ± 0.006 g·cm⁻³ for DGEBA/adipic acid from Kamble et al.); Tg at 349.15 (76 °C) vs ~42.8 °C reported for DGEBA/adipic acid; E = 4.29 (x), 2.67 (y), 6.13 (z), mean 4.36 ± 1.00 GPa; Poisson's ratio 0.37545 / 0.38780 / 0.34870, mean 0.371 ± 0.012. Caveats: Fig. 15 shows snapshots at 40 % strain though the protocol states 20 %; 3648 atoms with chains of only 8 repeat units, a single box, and 0.25 fs ReaxFF runs of 12.5 ps per anneal step are very short.

**CAPS relevance:** weak as a validation target (no exchange kinetics, tiny system, preprint). Useful pieces: (1) the Tg protocol (10 K cooling steps, spline, second-derivative zero) and modulus protocol (3-6 % strain window, three axes) are cheap to reproduce in CAPS's own engine on a DGEBA/diacid network built with OPLS-AA or PCFF; (2) the CV table gives a ready-made biasing recipe if CAPS ever exports Colvars input with a ReaxFF deck. For network building CAPS's fix bond/react path is the better route; the Kamble et al. DGEBA/adipic acid numbers (1.215 g/cm³, Tg ~42.8 °C) quoted here are a better benchmark than this paper's own values.

---

### 2. Ciarella, Sciortino, Ellenbroek 2018 (PRL), defects as a highway to stress relaxation

**Local file:** Ciarella2018_PRL_arXiv1804.01723.pdf (arXiv v1 text, dated 6 April 2018, with SI appendices A and B)

**Citation:** S. Ciarella, F. Sciortino, W. G. Ellenbroek, "Dynamics of Vitrimers: Defects as a Highway to Stress Relaxation", Physical Review Letters 121, 058003 (2018).

**DOI:** 10.1103/PhysRevLett.121.058003. Crossref: match (title, three authors, Phys. Rev. Lett. vol 121, issue 5, article 058003, published 2018-07-31). arXiv DOI 10.48550/arXiv.1804.01723 verified on DataCite, which records IsVersionOf -> 10.1103/physrevlett.121.058003. Published: yes. Note the local file is the arXiv preprint, so figures and wording may differ slightly from the PRL version.

**Chemistry:** generic. Two reactive end types (red/blue) stand for e.g. carboxyl and hydroxyl ends of ester vitrimers; the barrier stands in for catalyst concentration or temperature. Rubber relevance: none directly (star-polymer networks, loosely like tetra-PEG).

**Model level / force field:** coarse-grained bead-spring 8-arm star polymers. Units: length 1 nm (spring rest length, approx. one Kuhn segment of ~8 carbons), mass [m] = 100 u, energy [E] = kB·300 K, time [τ] = 6.33 ps. All beads: purely repulsive WCA, σ = 1 nm. Backbone: harmonic springs, rest length 1 nm (spring constant not stated). Reactive pair (red-blue only): generalised LJ 20-10, v = 4ε[(σ/r)^20 - (σ/r)^10], σ = 0.5 nm, ε = 100 kBT, rcut = 2.5σ.

**Exchange method:** Sciortino three-body swap potential (RevCross), v3 = λ ε v̂_ij v̂_ik with v̂ = 1 for r ≤ rmin and -v2b/ε for r > rmin, λ ≥ 1. Swap barrier to first approximation β ΔE_sw = βε(λ - 1) = 100(λ - 1). Continuous MD, no MC acceptance step, no swap interval (the rate comes from the potential). Associative, bond number conserved. Barriers scanned: β ΔE_sw from 0 to 90 (values in Figs. 2 and 4; plateau persists beyond reach for β ΔE_sw > 50). Temperature fixed (energy unit kB·300 K). Stress tensor includes a derived three-body virial term (SI eq. 46).

**System size and run length:** N ≈ 50000 beads, ≈ 1500 stars. Defect-free mixture (DFM): 900 stars with 8 A ends + 600 stars with 8 B ends (A-B only, primary loops impossible). Defect-allowing mixture (DAM): 950 stars with 7 A + 1 B, 550 stars with 7 B + 1 A. 4800 B ends all bonded, 2400 free A ends. Cubic box L = 40 nm, packing fraction φ ≈ 0.3 (2.2 × overlap concentration). 100 independent network topologies per mixture. Timestep and number of steps: not stated; simulated times reach tens of µs (τ_liq ≈ 20 µs).

**Software:** HOOMD-blue on GPUs with the authors' own module "RevCross" (later part of HOOMD-blue).

**Properties and exact analysis:**
- Stress relaxation: equilibrium Green-Kubo autocorrelation G(t) ≈ C(t) = (V/kBT) <σ_yz(t) σ_yz(0)> in nVT; no fit form applied (the authors avoid fitting a power law to the short-time Rouse regime). Solid vs liquid distinguished by the long-time limit C_∞ (G = C + G_eq - C_∞ for solids).
- Relaxation times: τ_net ≈ 5 ns (time to reach the elastic plateau); τ_liq defined by G(τ_liq)/G(0) ≡ 10^-4.
- Swap count normalised vs barrier (Fig. 4); topology: histogram of number of connections between connected stars; fraction of primary loops.
- Glass transition check (SI): compress DAM networks at β ΔE_sw = 0 from φ ≈ 0.3 to φ_max ≈ 0.61, decompress to 0.05, measure C(τ_liq) and bond decorrelation 1 - c_b(τ_liq); kinetic glass slowdown sets in at φ = 0.5. Tv analogue is the "topological glass transition" set by β ΔE_sw; no Arrhenius fit, no activation energy fitted.

**Key results:** DFM with no barrier liquid at τ_liq ≈ 20 µs; DAM τ_liq ≈ 0.1 τ_liq(DFM) ≈ 2 µs (ten times faster at a similar number of swaps); 1/5 to 1/3 of DAM bonds are primary loops; swap rate falls by two orders of magnitude as β ΔE_sw goes 0 -> 6; fluid-solid transition at a barrier of 10 kBT (Fig. 4 caption).

**CAPS relevance:** the standard CG vitrimer benchmark. First reproducible protocol: a bead-spring star or Kremer-Grest network with the 20-10 gLJ + three-body swap, measure Green-Kubo G(t) with the three-body virial included, and check (a) the solid-liquid crossover vs λ and (b) the 10x DAM/DFM gap. Requires CAPS's own engine to implement the three-body term and its virial (LAMMPS has no built-in RevCross; LAMMPS users instead use fix bond/swap or REACTER templates). A smaller first target is the Sciortino 2017 EPJE single-potential test (verified DOI 10.1140/epje/i2017-11496-5).

---

### 3. Holoman, Prajwal, Hocky, Truskett 2026 (COCIS review), simulating dynamic bonding

**Local file:** Holoman2026_COCIS_dynamic_bonding_review.pdf

**Citation:** T. R. Holoman, B. P. Prajwal, G. M. Hocky, T. M. Truskett, "Simulating dynamic bonding in soft materials", Current Opinion in Colloid & Interface Science 83, 102019 (2026). Open access CC BY-NC-ND.

**DOI:** 10.1016/j.cocis.2026.102019. Crossref: match (title, four authors, vol 83, article 102019, issued 2026-06).

**Scope:** short review of MD, MC, hybrid MD/MC and kinetic MC methods for associative and dissociative dynamic bonds; vitrimers, polymer-grafted nanoparticles, hydrogels, colloids and cytoskeletal networks. No rubber chemistry.

**Methods as categorised (Table 1 of the review):**
- RevCross three-body potential (MD, associative, HOOMD-blue). Eq. 1: V = λ Σ ε V3(r_ij) V3(r_ik); barrier ε(λ - 1); λ < 1 gives unphysical clustering. Limitation stated: monovalent only, weak at high density where collisions are incessant.
- Customised Tersoff bond-order bead-spring model (Liu and O'Connor, LAMMPS) for associative melts and extensional rheology.
- REACTER / fix bond/react (MD, associative and dissociative, LAMMPS): pre/post templates plus map file, bonding cutoff triggers the superimpose algorithm; since 2020 parallel, bond breaking, angle/dihedral constraints.
- Atomistic: QM/MM; ReaxFF (energy terms E_bond + E_over + E_under + E_lp + E_val + E_tors + E_vdW + E_Coulomb); accelerated ReaxFF with bond-boost restraint E_rest = F1{1 - exp[-F2 (r_ij - r12)^2]} (Vashisth), used for epoxy-amine, vitrimers and Diels-Alder CANs.
- MC: reactive canonical MC; explicit-bonding reaction ensemble MC (eb-RxMC, ESPResSo); multivalent associative bond-swap MC of Rao, Xia, Ni with P_acc = min[1, (p'/v'_u)/(p/v_u) exp(-β ΔG0)].
- Hybrid MD/MC: Hoy-Fredrickson sticky FENE bonds U_sb = U_FENE(r) - U_FENE(r0) - h, r0 ≈ 0.96σ; dybond (HOOMD-blue plugin, separate k_on and k_off); TILD in MATILDA.FT; associative swaps in MD/MC (Perego and Khabaz; Stukalin).
- Kinetic MC / agent-based for biopolymers: AFINES, Cytosim, aLENS, MEDYAN (Gillespie).
- Outlook: detailed vs global balance for multivalent moves, ML force fields for CANs, ML inverse design.

**Key primary simulation papers cited (DOIs verified on Crossref 2026-10-22 unless marked):**
1. Sciortino F., Three-body potential for simulating bond swaps in molecular dynamics, Eur. Phys. J. E 40, 3 (2017). 10.1140/epje/i2017-11496-5 [verified]
2. Ciarella S., Ellenbroek W. G., Associative bond swaps in molecular dynamics, SciPost Phys. 12, 128 (2022). 10.21468/SciPostPhys.12.4.128 [verified]
3. Gissinger J. R., Jensen B. D., Wise K. E., REACTER: a heuristic method for reactive molecular dynamics, Macromolecules 53, 9953-9961 (2020). 10.1021/acs.macromol.0c02012 [verified]; also Polymer 128, 211-217 (2017), 10.1016/j.polymer.2017.09.038 [verified] and Comput. Phys. Commun. 304, 109287 (2024), 10.1016/j.cpc.2024.109287 [verified]
4. Vashisth A. et al., Accelerated ReaxFF simulations for describing the reactive cross-linking of polymers, J. Phys. Chem. A 122, 6633-6642 (2018). 10.1021/acs.jpca.8b03826 [verified]
5. Zheng Y., Varshney V., Vashisth A., Accelerated ReaxFF simulations of vitrimers with dynamic covalent adaptive networks, Macromolecules 58, 4948-4958 (2025). 10.1021/acs.macromol.5c00501 [verified]
6. Perego A., Khabaz F., Volumetric and rheological properties of vitrimers: hybrid MD and MC, Macromolecules 53, 8406-8416 (2020). 10.1021/acs.macromol.0c01423 [verified]; follow-up Perego et al., Macromolecules 55, 7605-7613 (2022), 10.1021/acs.macromol.2c00588 [verified]
7. Rao P., Xia X., Ni R., A bond swap algorithm for simulating dynamically crosslinked polymers, J. Chem. Phys. 160 (2024). 10.1063/5.0186553 [verified; Crossref gives issue 6, article 061102, which the review does not print]
8. Singh G., Varshney V., Sundararaghavan V., Understanding creep in vitrimers (all-atom, REACTER), Polymer 313, 127667 (2024). 10.1016/j.polymer.2024.127667 [verified]; Singh and Sundararaghavan, self-healing of vitrimers, Chem. Phys. Lett. 760, 137966 (2020), 10.1016/j.cplett.2020.137966 [verified]
9. Zhao H. et al., MD of structural, mechanical and reprocessing properties of vitrimers, Macromolecules 55, 1091-1103 (2022). 10.1021/acs.macromol.1c02034 [verified]
10. Hoy R. S., Fredrickson G. H., Thermoreversible associating polymer networks I, J. Chem. Phys. 131, 224902 (2009). 10.1063/1.3268777 [verified]
Also cited and verified: Ciarella and Ellenbroek, Swap-driven self-adhesion and healing of vitrimers, Coatings 9, 114 (2019), 10.3390/coatings9020114; Sun Y. et al., ML force fields for CAN exchange, Macromolecules 56, 9003-9013 (2023), 10.1021/acs.macromol.3c01377 (DOI not printed in the review; found by Crossref query).

**CAPS relevance:** the review's Table 1 is a ready checklist for a CAPS "dynamic bond" module. For an all-atom builder with fix bond/react export the REACTER entries (Gissinger 2020/2024; Singh 2020/2024 vitrimer creep and healing) are the most direct templates to reproduce; for the own-engine CG path, RevCross (Sciortino 2017, Ciarella 2022) and the Hoy-Fredrickson/dybond MD/MC scheme are the two minimal algorithms to implement and cross-check.

---

### 4. Li et al. 2023/2024 (Advanced Science review), next-generation vitrimers by theory and simulation

**Local file:** Li2023_AdvSci_review.pdf

**Citation:** K. Li, N. V. Tran, Y. Pan, S. Wang, Z. Jin, G. Chen, S. Li, J. Zheng, X. J. Loh, Z. Li, "Next-Generation Vitrimers Design through Theoretical Understanding and Computational Simulations", Advanced Science 11 (5), 2302816 (2024; online 2023). Open access CC BY.

**DOI:** 10.1002/advs.202302816. Crossref: match (title, all ten authors, Adv. Sci. vol 11, issue 5, article 2302816, issued online 2023-12-07, print 2024-02). The file name says 2023; the volume year is 2024.

**Scope:** DFT, AIMD, classical MD (all-atom and CG), reactive FFs, hybrid MD-MC and ML potentials applied to vitrimers; properties (Tg, Tv, Arrhenius viscosity, stress relaxation, welding, creep, swelling); DFT for catalyst screening and functional-group design. Rubber relevance: disulfide vitrimers (Park 2020; Singh creep with 4-aminophenyl disulfide) only; no ENR/SBR/NR/PDMS simulations discussed.

**Simulation methods as categorised (review Table 1, system sizes as printed):**
- "Cut-off distance" atomistic MD: ≈10^0-10^1 nm, ≈10^1-10^2 ns (Yang 2015 Soft Matter: bonds made on proximity, topology accepted on new-bond energy; Singh 2020 adds explicit temperature-dependent reaction probabilities from experiment; Park 2020 probabilistic disulfide BER).
- Accelerated ReaxFF (bond boost, Miron-Fichthorn): ≈10^0-10^1 nm, ≈10^1 ps (Vashisth 2018; Kamble 2022 carbon-fibre vitrimer).
- Three-body potential CGMD: ≈10^1-10^2 nm, ≈10^1 µs (Sciortino 2017; Ciarella 2018).
- Hybrid CGMD-MC: ≈10^1-10^2 nm, ≈10^7 timesteps (Perego and Khabaz; Wu 2019; Stukalin 2013). Metropolis on ΔU_exchange = U_new(i,l) + U_new(j,k) - U_old(i,j) - U_old(k,l); origin in Pant-Theodorou end-bridging MC.
- Patchy particles (Smallenburg, Leibler, Sciortino 2013): Kern-Frenkel patches, event-driven MD, Wertheim theory.
- DFT: transition states by NEB / CI-NEB; functionals listed in review Table 2 (ωB97XD, M06, M06-2X, BP86, B3LYP).
- MLMD: Behler-Parrinello, q-SNAP, GAP; DeePMD-kit and others listed; challenges are data, transferability, cost.
- Network building: iterative crosslinking with relaxation, thermocycling anneal; 2D PBC with a free surface for welding.

**Property definitions summarised:** Tv defined experimentally at viscosity 10^12 Pa·s; in Perego and Khabaz 2020 Tg from volumetric (specific volume vs T̃) data and Tv from the thermal expansion coefficient vs T̃; bond lifetime Arrhenius-like at T̃ ≤ 1.7, temperature-insensitive at T̃ ≥ 1.8. Zhao 2022: stress relaxation after a step uniaxial strain of 300 %, self-healing efficiency vs ΔE_sw after triaxial cavitation. Welding (Yang 2016): modulus and yield strength of welded vs fresh samples vs time, T, degree of polymerisation, crosslink density.

**DFT numbers quoted:** transcarbamoylation Ea 35.5 kJ/mol at 90° OCNC torsion vs 70.7 kJ/mol relaxed (Fortman); Me-AE amine reactions ΔΔG0 < -95 kJ/mol, Me-AA -30 to 15 kJ/mol; catalyst order TBD ≳ DETO ≳ Zn(OAc)2 > 1-MI (Bhusal 2021); zinc acetate barrier 51.84 kcal/mol, zinc pivalate 51.36 kcal/mol (Shen 2022, LCP synthesis, not a vitrimer).

**Key primary simulation papers cited (no DOIs printed in the review; all found and verified by Crossref query 2026-10-22):**
1. Yang H. et al., A molecular dynamics study of bond exchange reactions in covalent adaptable networks, Soft Matter 11, 6305-6317 (2015). 10.1039/c5sm00942a [verified]
2. Smallenburg F., Leibler L., Sciortino F., Patchy particle model for vitrimers, Phys. Rev. Lett. 111, 188002 (2013). 10.1103/PhysRevLett.111.188002 [verified]
3. Perego A., Khabaz F., Macromolecules 53, 8406 (2020). 10.1021/acs.macromol.0c01423 [verified]; and Effect of bond exchange rate on dynamics and mechanics of vitrimers, J. Polym. Sci. 59, 2590-2602 (2021), 10.1002/pol.20210411 [verified]
4. Zhao H. et al., Macromolecules 55, 1091 (2022). 10.1021/acs.macromol.1c02034 [verified]
5. Ciarella S., Sciortino F., Ellenbroek W. G., PRL 121, 058003 (2018). 10.1103/PhysRevLett.121.058003 [verified]
6. Ciarella S., Biezemans R. A., Janssen L. M. C., Understanding, predicting, and tuning the fragility of vitrimeric polymers, PNAS 116, 25013-25022 (2019). 10.1073/pnas.1912571116 [verified]
7. Wu J.-B. et al., Dynamics and reaction kinetics of coarse-grained bulk vitrimers, Phys. Chem. Chem. Phys. 21, 13258-13267 (2019). 10.1039/c9cp01766f [verified]
8. Park C. et al., Enhanced self-healing performance of graphene oxide/vitrimer nanocomposites (disulfide BER, MD), Polymer 206, 122862 (2020). 10.1016/j.polymer.2020.122862 [verified]
9. Kamble M. et al., Reversing fatigue in carbon-fiber reinforced vitrimer composites (accelerated ReaxFF), Carbon 187, 108-114 (2022). 10.1016/j.carbon.2021.10.078 [verified]
10. Bhusal S. et al., Transesterification in vitrimer polymers using bifunctional catalysts (DFT), J. Phys. Chem. B 125, 2411-2424 (2021). 10.1021/acs.jpcb.0c10403 [verified]
Also cited: Singh and Sundararaghavan CPL 760, 137966 (2020) [verified above]; Perego and Khabaz, Macromol. Rapid Commun. 44, 2200313 (2023) [from review, unverified]; Yang et al. RSC Adv. 6, 22476 (2016) welding [from review, unverified].

**CAPS relevance:** most useful as a map of which atomistic protocols exist. The Yang 2015 cut-off BER scheme on DGEBA/tricarballylic acid and the Park 2020 disulfide BER are the closest all-atom analogues to what CAPS can build and export through fix bond/react; Perego and Khabaz 2020 gives the CG Tg/Tv-from-volumetric-data protocol that CAPS's own engine could reproduce as a first validation.

---

### 5. Rispo Constantinou et al. 2026 (Macromolecules), entropic phase separation in polymer-vitrimer melts

**Local file:** RispoConstantinou2026_Macromolecules_blends.pdf

**Citation:** A. A. Rispo Constantinou, B. Magyari, G. Ianniruberto, E. van Ruymbeke, D. J. Read, "Entropic Phase Separation in Polymer-Vitrimer Melts", Macromolecules 59 (19), 10615-10630 (2026). CC-BY 4.0. Received 27 March 2026, published 23 September 2026.

**DOI:** 10.1021/acs.macromol.6c00893. Crossref: match (title, five authors, vol 59, issue 19, pages 10615-10630, online 2026-09-23). A replication dataset is registered as 10.14428/dvn/aabcvq. An arXiv preprint is mentioned but its DOI is cut off in the text file.

**Chemistry:** generic associative vitrimer blended with a chemically identical homopolymer; "realistic" parameter case modelled on a polystyrene vitrimer with dioxaborolane cross-linker and stopper (Leibler-type dioxaborolane metathesis). No catalyst. Rubber relevance: none directly, but the result (vitrimer/thermoplastic blends phase separate entropically) matters for any blend of a vitrimer additive into a rubber or thermoplastic matrix.

**Model level / force field:** mean-field lattice theory (Semenov-Rubinstein / Dobrynin, Flory-Huggins base, nearest-neighbour same-chain correction with free parameter k) validated by Kremer-Grest bead-spring MD. WCA LJ truncated at rc = 2^(1/6)σ; FENE K = 30 ε/σ², R0 = 1.5σ, used for backbone bonds and for cross-links; fully flexible, no angle potential; χ = 0 by identical interactions.

**Exchange method:** hybrid MC/MD extending Cui et al. (ref 32). Two pathways: "bond swap" (metathesis of two bonded sticker pairs) and "bond shift" (free sticker + bonded pair). All participating stickers within rc = 1.12σ (WCA minimum). ΔU = Σ U_FENE(new) - Σ U_FENE(old); of the two possible new configurations only the lower-ΔU one is kept; Metropolis acceptance P_acc = exp(-ΔU/kBT*) for ΔU > 0. One topology change per sticker per MC cycle; satisfies global (not detailed) balance. Exchange attempts every 10^3 MD steps. Associative: cross-link number fixed. No barrier or activation energy (barrierless apart from FENE strain). T* = 1 ε/kB.

**System size and run length:** 2000 chains of N = 20 beads (40000 beads), ρ = 0.85 σ⁻³, box L ≫ 5 Rg. Vitrimer fraction ϕ = 0.1-0.9 in steps of 0.1; s = 5 stickers per vitrimer chain at beads 2, 6, 10, 14, 18; sticker conversion p varied in steps of 0.05. Protocol: soft-potential push-off 2 × 10^8 steps (δt 0.001 then 0.005 τ), Langevin τ_damp = 0.1τ; then WCA, 10^8 steps at δt = 0.001τ; cross-linking by LAMMPS fix bond/create every 10 steps, capture radius R_min = 1.4σ, P_bond = 0.0001, conversion monitored by an external Python script; then NVT 10^8 steps at δt = 0.01τ with permanent cross-links; production with exchanges: Langevin τ_damp = 0.5τ, δt = 0.01τ, 10^8 steps (t = 10^6 τ).

**Software:** LAMMPS (version not stated; fix bond/create) plus own Python MC exchange code, released as the PyPI package "MC-exchange"; data on the ReBond dataverse.

**Properties and exact analysis:**
- Phase stability, three ways, on 7³ = 343 sub-boxes (~5σ edge), grid shifted 10 times per axis (343 000 sub-domains): (1) histogram of vitrimer bead counts (bimodal = separated); (2) spatial composition correlation Σ[ϕ(i) - ϕ][ϕ(j) - ϕ] vs centre distance, phase separated if anticorrelation at large r reaches -0.1 (heuristic, box-size dependent); (3) density variance normalised to t = 0: < 1 stable, 1-1.5 ambiguous, > 1.5 unstable. Equilibrium reached by the half-way point (5 × 10^7 steps).
- Percolation: reduced molecular weight (all clusters except largest) maximum and inflection of largest-cluster size (Livraghi method), compared with Flory-Stockmayer p - q = 1/[(1 - q)(s - 1)].
- Theory: spinodal from zero of the Hessian determinant, plait point from eq. 13, binodal by common tangent with tie-line stepping, lower convex envelope check.
- Not computed: stress relaxation, Tg, Tv, modulus, MSD, healing.

**Key results:** fitted k = 0.48 ± 0.08 (Gaussian-chain estimate 0.39, or 0.88 with v0 = 1/ρ; from measured p,q 0.20 nearest-neighbour or 0.30 any same-chain; from P[l4] ≈ 0.3, k = 0.48); at ϕ = 0.5 gel point p = 0.4 (RMW max) to 0.475 (inflection); all clearly separated simulations lie inside the predicted unstable region. Asymptotic rule (p - q)s = 1 + sqrt(N/M): for N ≈ M, phase separation once two stickers per chain bond to other chains, vs one per chain for gelation. Realistic PS case: dp = 1500, s = 200, k ≈ 0.34 (175 kDa vitrimer, 156 kDa matrix); instability at much lower p; highly cross-linked blends separate into nearly pure phases.

**CAPS relevance:** a clean, fully specified CG protocol that CAPS's own engine can reproduce end to end: KG melt, fix bond/create-style cross-linking to a target p, then MC swap/shift moves with FENE ΔU and Metropolis, analysed by the sub-box density variance. The open MC-exchange package allows a direct cross-check. A good first validation target is the ϕ = 0.5 gel point p = 0.4-0.475 and the stable/unstable classification at ϕ = 0.5 vs p (Figs. 6-7).

---

### 6. Stukalin, Cai, Kumar, Leibler, Rubinstein 2013 (Macromolecules), self-healing of unentangled networks with reversible bonds

**Local file:** Stukalin2013_Macromolecules_selfheal.pdf (NIH author manuscript, PMC)

**Citation:** E. B. Stukalin, L.-H. Cai, N. A. Kumar, L. Leibler, M. Rubinstein, "Self-Healing of Unentangled Polymer Networks with Reversible Bonds", Macromolecules 46 (18), 7525-7541 (2013).

**DOI:** 10.1021/ma401111n. Crossref: match (title, five authors, vol 46, issue 18, pages 7525-7541, 2013).

**Chemistry:** generic pairwise reversible associations (motivated by Leibler's supramolecular hydrogen-bonded self-healing rubber); the stickers form and break (dissociative), with "hopping"/partner exchange as an emergent path. No catalyst. Rubber relevance: indirect; the model is for elastomers and gels above Tg, and the hybrid permanent-network + dangling-sticker picture fits a covalently cured rubber carrying reversible pendant groups.

**Model level / force field:** scaling theory plus hybrid MD/MC Kremer-Grest simulations. Dangling chains anchored at one end, sticker at the other. WCA LJ rc = 2^(1/6)σ_LJ; FENE R_FENE = 1.5σ_LJ (spring constant in an equation image, lost from the text file), b0 ≈ 0.96σ_LJ; T = 1.0 ε_LJ/kB; monomer density in an equation image, lost. Sticky bond potential U_SB = U_FENE minus a constant offset U_FENE(b0) + ε, where ε is the bond strength.

**Exchange method:** Metropolis MC formation/breaking of sticky bonds (dissociative, not a swap). Every τ_MC = τ_LJ each candidate pair (sticker pair closer than R_FENE) is selected once on average; bonded pair broken with Metropolis on U_SB, unbonded pair bonded with Metropolis; a sticker already bonded to another partner is left unchanged; previously selected pairs are not excluded within a cycle. MD: velocity Verlet, δt = 0.01τ_LJ, Langevin damping 10 τ_LJ. Bond strengths: ε-jumps εi = 0 or 4 kBT -> εf = 12 kBT (N = 10); εi = 0 -> εf = 10 kBT (N = 2). Theory regimes: kBT ln N < ε < 2 kBT ln N (anomalous Rouse sub-diffusion controls rebinding) and ε > 2 kBT ln N (hopping / partner exchange). Bare bond lifetime τ_b ∝ exp(ε/kBT) (eq. 10).

**System size and run length:** M = 2000 dangling chains of N = 10, and M = 1000 anchored dimers N = 2; 60 trajectories from 3 initial configurations; box L = 28.65σ_LJ; reaction/cut zone width R0 = 3.6σ_LJ for N = 10. Runs reach t/τ_LJ > 2 × 10^5 (hopping regime, N = 2). Healing protocol: equilibrate melt and anchor one end, switch on stickers, "ideal cut" in the xy mid-plane breaking all bridges with a confining wall potential, waiting period (no cross-plane bonds), shift halves by half the box in x, remove wall, heal.

**Software:** enhanced version of DL_POLY 2.12.

**Properties and exact analysis:**
- Open-sticker concentration c_open(t) after ε-jumps, fitted with the interpolation formula eq. 27 / A6 (two parameters: final c_open and constant C); hopping regime fitted with eq. A9 (rate constant and crossover concentration); waiting-period decay of open-sticker surface density σ_open = N_open/L² in the cut zone fitted with eq. A10 (C_w).
- Healing measure: surface density of bridges σ_b(t) across the interface (theory: direct bonding, then "stitching" by partner exchange, logarithmic growth at long times). No mechanical healing efficiency, no stress relaxation, Tg, modulus or Arrhenius fit in the simulations.

**Key results:** C = 3.2 ± 0.1 (same for both ε-jumps); exploration-volume prefactor B = 0.3, so A = B/C ≈ 0.1 (≈ 0.12 by direct counting); early bonding regime 0 < t/τ_LJ < 50 then hopping for N = 2; τ_b = 3200 τ_LJ at ε = 10 kBT; C_w = 1.6 ± 0.1. Theory: renormalised bond lifetime ≫ bare lifetime because of repeated returns to the old partner; full interfacial recovery takes τ_adh ≫ τ_eq; shorter waiting time gives faster bridge formation (healing vs self-adhesion).

**CAPS relevance:** defines the "bridges across a cut" healing metric and the cut / wait / heal protocol, which CAPS can reproduce with any bond-exchange engine (own MD/MC or LAMMPS fix bond/create + fix bond/break with probabilities). First reproducible target: N = 10 anchored chains, ε-jump 0 -> 12 kBT, recover C ≈ 3.2 from the c_open(t) fit, then run the cut/wait/heal sequence and count bridges vs waiting time. Several equations are images lost from the text file; re-read the PDF for FENE K, density and eqs. A6 to A10 before implementing.

---

### 7. Wilcoxson et al. 2026 (arXiv), pySNAP and DyBond-v2

**Local file:** Wilcoxson2026_arXiv2609.21204_pySNAP.pdf

**Citation:** T. A. Wilcoxson, T. Li, T. R. Holoman, G. Mitra, T. M. Truskett, G. M. Hocky, "Simulation and Network Assembly Pipelines for Dynamically Bonded Soft Materials", arXiv:2609.21204 [cond-mat.soft], v1 18 September 2026 (dated 21 September 2026).

**DOI:** 10.48550/arXiv.2609.21204 (DataCite: title and six authors match, year 2026, no related published version recorded). Crossref: arXiv DOIs are not in Crossref; a Crossref bibliographic query on the title returned no journal version (top hits were unrelated, plus the group's earlier Mitra et al. Soft Matter 2023 paper). Status: unpublished preprint as of 2026-10-22.

**Chemistry:** generic dynamic bonds (thia-Michael, hydrazone-type pairs, DNA binders); examples are patchy colloids, tetraPEG/nanocrystal gels and DNA-coated emulsion droplets. Not a vitrimer method: the paper states the scheme is purely dissociative and does not implement associative swaps. Rubber relevance: none.

**Model level / force field:** coarse-grained reduced-unit MD (σ, m, kBT, τ = σ sqrt(m/kBT)). Patchy colloid: patches held by FENE + WCA, octahedral arrangement by cosine-squared angle potential. Dynamic bonds: Harmonic, FENEWCA or Tether bonds (example harmonic k = 30.0, r0 = 1.0).

**Exchange / bonding method (DyBond-v2):** every n steps a sweep proposes bonds for each free reactive site (one partner within [r_min, r_max], nearest to r0), resolves conflicts so pairs are disjoint, accepts with P_on(r) = k_on n Δt exp(-U_b(r)/kBT) (T = instantaneous kinetic T), then breaks existing bonds with P_off = k_off n Δt; ε ≡ -ΔG/kBT = ln(k_on/k_off). Requires k n Δt ≤ 1. Multiple competing reactions share one occupancy pool, order randomised each sweep. Optional sigmoidal melting switch and Bell slip-bond P_off(r) = P0_off exp(a sqrt(U_b/kBT)). Intramolecular exclusion flag selfAvoidFlag. Dissociative only. No activation energy.

**System size and run length:** example A: 2000 colloids with 6 patches, η = 0.05, ε = 8 kBT, k_on = 2, 10, 50 τ⁻¹, trigger every 100 steps (every 10 for k_on = 50); example B: tPEG (four 6-bead arms) + NC mixtures from 400 NCs/2400 tPEGs to 2000 NCs/0 tPEGs, η = 0.05; Δt = 0.001τ, equilibration 5 × 10^5 steps (500τ), production 8 × 10^7 steps (8 × 10^4 τ). Example C: two droplets, D = 50σ to 200σ, 100 to 800 binding sites, bond length 1.97σ, Δt = 0.0005τ, 2 × 10^5 and 4 × 10^7 steps (100τ and 2 × 10^4 τ).

**Software:** HOOMD-blue with the DyBond-v2 C++/HIP plugin (CUDA and ROCm, Philox RNG); pySNAP = snap_simulate + snap_workflow (YAML parameter trees, DAG of stages, Slurm batching). Code on GitHub (hocky-research-group).

**Properties and exact analysis:** dynamic bond fraction vs time; bond correlation function (fraction of bonds present at t0 still intact after t); self-intermediate scattering F_s(k*, t) at the S(k) peak k* = 1.334σ⁻¹, fitted with stretched exponentials; static S(k) from 60 frames in [60000, 80000]τ; bond multiplicity per molecular pair; droplet patch outer radius (mean of five outermost bonded sites, final 50 frames). Validation: frozen four-reaction geometry vs exact partition function Z = 1 + K_AB + K_AC + K_DB + K_DC + K_AB K_DC + K_AC K_DB, populations within 0.02. No stress relaxation, Tg, Tv, healing or modulus.

**Key results:** stretched exponents 0.54, 0.62, 0.61 for increasing k_on; raising k_on from 2 to 10 more than doubles bonds formed in the first 50τ, 10 to 50 adds less than a fifth (reaction- to diffusion-limited crossover); patch outer radius r = c sqrt(D) with c = 1.53.

**CAPS relevance:** a well-specified, testable dissociative MD/MC bonding rule (P_on with Metropolis factor on the actual bond potential, P_off constant, disjoint matching) that CAPS's own engine could adopt for reversible (non-vitrimer) crosslinks such as hydrogen-bond or ionic rubbers. The frozen-geometry partition-function test (Appendix A) is the cheapest correctness check to copy into CAPS's test suite. For vitrimers it must be paired with a swap move, which the authors explicitly leave out.

## Group B

---

### 1. Arbe2023_ACSMacroLett_Tv_imine_PI

**Local file:** Arbe2023_ACSMacroLett_Tv_imine_PI.pdf

**Citation:** Arbe, A.; Alegría, A.; Colmenero, J.; Bhaumik, S.; Ntetsikas, K.; Hadjichristidis, N. Microscopic Evidence for the Topological Transition in Model Vitrimers. ACS Macro Lett. 2023, 12 (11), 1595-1601.

**DOI:** 10.1021/acsmacrolett.3c00586. Crossref: title, all six authors, journal, vol 12, issue 11, pages 1595-1601, published 2023-11-10 all match. Verified.

**Nature of paper:** Purely experimental. No simulation of any kind. Note: the structural probe is X-ray diffraction ("XR diffraction"), not neutron scattering. DSC and parallel-plate viscometry complete the set.

**Chemistry**
- Exchange chemistry: imine dynamic covalent bonds; two pathways named, transimination and imine metathesis. Literature activation energies quoted for these processes: 10 to 130 kJ/mol (cited review, not measured here).
- Polymer: α,ω-dialdehyde-functionalised polyisoprene (PI), microstructure about 65% 1,4-cis, 27% 1,4-trans, 8% 3,4.
- Crosslinker: tris(2-aminoethyl)amine (triamine). Small excess (≈5%) of primary amines left on purpose to help exchange. Unreacted material washed out with solvent.
- Catalyst: none (excess amine acts as promoter).
- Associative (transimination / metathesis, constant crosslink count).
- Rubber relevance: high. PI is the synthetic analogue of NR; the system is an NR-like vitrimer elastomer.

**Samples (Table 1)**

| Sample | Mn (kDa) | PDI | φt (triamine vol. fraction) | Tg (K) | D (Å) at 300 K | ⟨Re,o⟩ (Å) | βg (K⁻¹) | βv (K⁻¹) | βl (K⁻¹) |
|---|---|---|---|---|---|---|---|---|---|
| 2k-vit | 2.7 | 1.05 | 0.033 | 221.0 | 55 | 40 | 1.23 × 10⁻⁴ | 4.99 × 10⁻⁴ | 8.42 × 10⁻⁴ |
| 6k-vit | 6.3 | 1.04 | 0.014 | 216.5 | 72 | 61 | 1.51 × 10⁻⁴ | 4.33 × 10⁻⁴ | 7.18 × 10⁻⁴ |
| 11k-vit | 11.3 | 1.04 | 0.008 | 217.0 | 80 | 82 | 1.36 × 10⁻⁴ | 4.33 × 10⁻⁴ | 4.89 × 10⁻⁴ |

Triamine concentration 0.8% (11k-vit) to 3.3% (2k-vit). Reference Tg of standard PI with similar microstructure: 204-205 K.

**Model level / force field / exchange method / software:** not applicable (experiment).

**Properties and how they were analysed**
- Tg: DSC inflection point at about 1 K/min heating/cooling (Table 1).
- Viscosity Tv (Tv^η0): parallel-plate viscometry gives zero-shear viscosity η0(T); Arrhenius fit; extrapolated to η0 = 10¹² Pa·s. Ea(2k-vit) = 66 ± 20 kJ/mol, Ea(6k-vit) = 79 ± 7 kJ/mol. Tv^η0(2k-vit) = 272 ± 14 K, Tv^η0(6k-vit) = 274 ± 5 K. 11k-vit too viscous in the accessible range. Authors warn Arrhenius extrapolation to low T may not hold.
- Microscopic Tv ("microscopic dilatometry"): XRD at many T. Interchain peak Qc ≈ 1.3 Å⁻¹ (d = 2π/Qc ≈ 4.8 Å) and intercluster peak Qv ≈ 0.1 Å⁻¹ (D = 2π/Qv, 55 to 80 Å). Expansion coefficients from slopes of ln d(T) (βg below Tg, βv above Tg). ln D(T) is fitted with the same βg and βv slopes; Tv is where ln D(T) departs upward (extra mesoscale expansion), fitted above it with slope βl. Departure found at ≈260-270 K for all three samples; exact shaded Tv intervals are in Fig. 4 only [unverified numerically]. Intensity of the low-Q peak also drops above Tv (clearest signature for 11k-vit).
- DSC aging protocol for Tv: slow cooling at 0.25 K/min down to Tg, quench to 180 K, then heat at 20 K/min; compared with fast cooling (≈70 K/min). A broad endothermic excess appears, area proportional to crosslink density; standard DSC only shows a very weak feature above ≈280 K.
- Nanostructure: linkers aggregate into clusters of ≈20 triamines (≈50 chain ends) per cluster (SI estimate). Homogeneous linker spacing would be ≈19 Å (2k) to ≈35 Å (11k), much less than measured D.
- Stress relaxation, creep, modulus, MSD, healing: not stated.

**Key quantitative results**
- Tg 216.5 to 221.0 K (vitrimers) vs 204-205 K (plain PI).
- Tv ≈ 260-270 K by XRD; Tv^η0 = 272 ± 14 K and 274 ± 5 K by viscosity; agreement Tv^η0 ≈ Tv.
- Flow activation energies 66 ± 20 and 79 ± 7 kJ/mol.
- Chain expansion coefficient unchanged across Tv; only the mesoscale (cluster) spacing shows the transition.

**CAPS relevance (validation targets, no protocol to copy)**
- Best experimental anchor in this group for a rubber vitrimer: cis-PI (NR-like) + imine exchange, with Tg, Tv and Ea per Mn.
- Reproducible first in CAPS: an all-atom (OPLS-AA or GAFF) PI-dialdehyde + tris(2-aminoethyl)amine network built by fix bond/react (aldehyde + amine to imine + water), then a volume-temperature scan to reproduce the Tg trend (221 / 216.5 / 217 K vs 204 K for neat PI) and the interchain distance d ≈ 4.8 Å and its expansion coefficients βg, βv (compare with ln d(T) from the first S(Q) peak computed in CAPS).
- Harder second step: imine exchange through fix bond/react (transimination with a free amine) or a CG swap model, then check that the mesoscale cluster peak (D = 55 to 80 Å) appears and that its expansion steepens above Tv. Cluster formation (≈20 triamines per cluster) is a structural target a CG model could reproduce.

---

### 2. Ciarella2019_Coatings_adhesion

**Local file:** Ciarella2019_Coatings_adhesion.pdf

**Citation:** Ciarella, S.; Ellenbroek, W. G. Swap-Driven Self-Adhesion and Healing of Vitrimers. Coatings 2019, 9 (2), 114.

**DOI:** 10.3390/coatings9020114. Crossref: title, both authors, Coatings vol 9, issue 2, article 114, published 2019-02-12 match. Crossref also lists a preprint (10.20944/preprints201902.0065.v1). Verified.

**Chemistry**
- Generic associative exchange; binding groups are meant to represent e.g. carboxyl and hydroxyl ends, so bonds are esters and the swap is transesterification. Red and black end beads are distinct types.
- Polymer: 8-arm star polymers, arm = 5 beads. Mixture of NA = 950 stars with 7 red + 1 black end (7R1Bl) and NB = 550 stars with 1 red + 7 black (1R7Bl).
- Catalyst: not modelled explicitly (only through the swap barrier λ).
- Rubber relevance: generic elastomer network ("stretched and dried elastomer" visualisation); no specific rubber chemistry.

**Model level and force field**
- Coarse-grained bead-spring MD (star polymers), not Kremer-Grest FENE.
- Bonds within arms: harmonic springs, rest length L0 (unit of length). Spring constant: not stated.
- Non-bonded: WCA, σ = L0, cutoff 2^(1/6)σ.
- Reactive red-black bond: generalized short-ranged attractive Lennard-Jones (exponents not stated in this paper), binding energy ε_bond of order 100 kBT.
- Swap: Sciortino three-body potential (Eur. Phys. J. E 2017, 40, 3) with parameter λ; swap barrier ΔE_swap = ε_bond(λ − 1). λ = 1 is barrier-free.

**Exchange method**
- Algorithm: Sciortino three-body potential inside MD (no MC). Associative.
- Barrier: λ = 1 (ΔEsw = 0) for the healing runs; ΔEsw = 0.1ε, 0.5ε and 5ε also run (Figs 5d, 6, 8).
- Swap attempt frequency / acceptance rule: none; swaps are continuous dynamics. Time unit τsw = time for the average number of swaps per bond to reach 1 without barrier.
- Cutoff distances: not stated (in ref. 19, Ciarella et al. PRL 2018).
- Temperature: T = 1 (reduced).

**System size and run length**
- 1500 stars, 4800 bonds, 2400 unbound majority (red) ends; total N = 49,500 monomer beads per network (as printed). Packing fraction Φ ≈ 0.29 (below glass transition; 2.2× overlap concentration).
- M = 25 independent networks equilibrated with walls in z, periodic in x, y; 10⁷ timesteps of 4.2 × 10⁻⁵ τsw each.
- 300 adhesion simulations (pairs stacked, fully periodic), each giving two healing cuts. Healing followed up to tsh = 20 τsw (density profiles) and 4.2 τsw (bond profiles).

**Software:** HOOMD-blue with the authors' three-body swap implementation.

**Properties and how they were analysed**
- Mechanical healing: topology frozen at chosen tsh, then stress relaxation modulus Gij(t) by the stress autocorrelation method (Kriuchevskyi et al. 2017); long-time plateau = static shear modulus. Healing = Gxz and Gyz reaching Gxy (Gxy unaffected by the cut). Cumulative distribution Gc_ij = ∫ dG P(G) also shown.
- Interfacial bonds: B(z), number of star-star bonds crossing a plane at height z; cut at z ≈ 42σ.
- Bridges: Nbr = number of bonds between monomers from different halves; compared with Nbr(t) = (2/Lz) ∫ dz ρ1 ρ2 / ρ0².
- Interdiffusion: density ρ1(z,t) of stars originally in lower half, fitted with 1D diffusion equation from square-wave initial condition (Fourier series, slowest mode q = 2π/Lz ≈ 0.075, and Green function solution).
- Bond recovery near the cut: Nb,cut in 40.5σ ≤ z ≤ 43.5σ normalised by equilibrium bulk value; plotted as 1 − Nb,cut/Nb,equi.
- Stress relaxation of the bulk material: taken from ref. 19 (τrel ≈ 13.2 τsw).
- Swap location map: first 100 swaps after contact are uniform, no enhancement at the cut.

**Key quantitative results**
- Mechanical healing complete at tsh ≈ 0.15 τsw; components comparable by ≈0.1 τsw, indistinguishable before 1 τsw. Visual healing indistinguishable from bulk from tsh = 0.42 τsw.
- Star diffusion D ≈ 3.12 σ²/τsw (ΔEsw = 0); appreciable interdiffusion only after ≈5 τsw (30× slower than healing).
- τrel ≈ 13.2 τsw, so healing is about three orders of magnitude faster than flow.
- ΔEsw = 0.5ε already strongly suppresses bridge formation.
- Turning swaps off gives D = 0 (layered profile does not relax).

**CAPS relevance**
- Most direct protocol to reproduce first for a CG welding module: two equilibrated slabs (walls in z) brought into contact, then count bridges across the cut and B(z) per tsh, plus frozen-topology Green-Kubo G(t) components Gxz, Gyz vs Gxy. All three metrics are cheap and unit-free.
- Needs a three-body swap potential (Sciortino) in the CAPS engine or LAMMPS (pair_style sw-like custom). If CAPS only has fix bond/react or fix bond/swap, an MC-swap analogue with a tunable barrier could substitute, and the bridge count vs time and the 0.15 τsw / 5 τsw / 13.2 τsw time-scale separation are the validation targets.

---

### 3. Hou2026_CEJ_AI_hybrid_reactive

**Local file:** Hou2026_CEJ_AI_hybrid_reactive.pdf

**Citation:** Hou, J.-a.; Zhao, X.; Zhao, D.; Wang, X. Q.; Tam, L.-h.; Chow, C. L.; Lau, D. AI-accelerated hybrid reactive-nonreactive molecular dynamics simulation for polymer cross-linking. Chem. Eng. J. 2026, 549, 182583.

**DOI:** 10.1016/j.cej.2026.182583. Crossref: title, all seven authors, Chemical Engineering Journal vol 549, article 182583, issued 2026-12 match (paper says available online 30 September 2026). Verified.

**Chemistry**
- Not a vitrimer: irreversible epoxy-amine curing (no exchange).
- Polymer / monomers: DGEBA + triethylenetetramine (TETA); epoxide ring opening by amine N to give C-N bond and hydroxyl.
- Catalyst: none.
- Rubber relevance: none directly (thermoset epoxy for FRP). Relevant as a method for curing simulation.

**Model level and force field**
- Atomistic. ReaxFF (reactive) switched globally with CVFF (non-reactive). ReaxFF parameter set: not stated in the main text [unverified, may be in SI]. CVFF atom types reassigned at each switch from bonding environment; charges taken from converged ReaxFF QEq at the switching frame and held fixed during CVFF.
- Bond-order threshold BO_th and proximity cutoff r_cut: defined but values not stated in the main text [unverified].

**Reaction / switching method**
- ReaxFF does the chemistry. An XGBoost regressor predicts Preact (future normalised bond-formation activity) from four descriptors per 25 fs window (100 ReaxFF steps at 0.25 fs): Rbond (new-bond rate per candidate pair), ΔBO (bond-order fluctuation), Nclose (fraction of candidate pairs within r_cut), ΔE (potential energy change per atom).
- Switch ReaxFF to CVFF when Preact stays below Plow = 0.20; return when Preact > Phigh = 0.50 (hysteresis). In CVFF, every τmon = 7.5 ps (7500 CVFF steps at 1 fs) a copy is probed with ReaxFF for τprobe = 25 fs. Maximum cumulative CVFF time τmax = 500 ps. After CVFF to ReaxFF, minimisation and short NVT at 300 K.
- Standalone CVFF reference network made with LAMMPS fix bond/create (details in SI).
- Temperatures: dataset trajectories at 500 K; hybrid production runs at 300 K.

**System size and run length**
- Sizes: 301 atoms (5 DGEBA + 2 TETA), 3010, 15,050, 25,284, 50,267 atoms. Property analysis on the 50,267-atom system.
- Dataset protocol: minimisation, NPT 300 K, 1 atm, 250 ps; NVT heat 300 to 500 K in 100 ps; NVT 500 K for 500 ps; quench to 300 K in 100 ps; NPT 300 K 100 ps. Nosé-Hoover thermostat and barostat. 0.25 fs ReaxFF timestep.
- 500 trajectories for training/validation (350 of 301 atoms, 150 of 3010 atoms), 42,000 configurations each, 21,000,000 in total, ≈1 ns each; split 400 train / 100 validation by trajectory; plus 10 trajectories each for 15,050, 25,284 and 50,267 atoms in validation only (130 validation trajectories).
- Hybrid production run length: at least 11 ns (insets show 10 to 11 ns) [total not stated in text].

**Software:** LAMMPS named only for fix bond/create; ReaxFF and CVFF engine assumed LAMMPS [unverified]. Python conversion scripts (SI Code S1, S2). XGBoost, GridSearchCV with 10-fold trajectory-grouped CV. Hardware: Hygon C86 7285 32-core 2.5 GHz, 256 GB, 32 cores.

**Properties and how they were analysed**
- ML accuracy: R², RMSE, MAE on train and held-out trajectories; SHAP and leave-one-descriptor-out ablation.
- Crosslink degree vs time; gel point and largest connected cluster fraction (SI Fig. S4, Table S4).
- RDF for N-O, C-C, C-O, C-N.
- MSD (qualitative comparison, no D extracted).
- Density vs time.
- Tg: kink in density vs temperature from the hybrid run (cooling rate not stated in main text).
- Uniaxial tension stress-strain; Young's, bulk, shear modulus and Poisson ratio (method of extraction in SI [unverified]).
- Speed: elapsed time per MD step per core.

**Key quantitative results**
- XGBoost: train R² 0.9949, RMSE 0.0100, MAE 0.0076; validation R² 0.9616, RMSE 0.0395, MAE 0.0355 (Random Forest R² 0.9603, MLR 0.8890).
- With (0.20, 0.50): final crosslink degree 81.9%, elapsed-time reduction 83.1%; all threshold pairs gave 77.4 to 82.5% crosslinking and > 77.8% reduction.
- Per-step cost (s per step per core) CVFF / ReaxFF / hybrid: 301 atoms 0.00098 / 0.07247 / 0.01511 (79.1%); 50,267 atoms 0.09750 / 8.20001 / 1.23003 (85.0%).
- Gel point at 35.5% (hybrid) vs 39.3% (full ReaxFF) crosslinking; largest cluster 0.993 vs 0.973.
- Tg = 447.1 K (hybrid) vs 419.35 K experiment and 440.35 K earlier MD (ref. 62).
- Moduli (GPa), CVFF / ReaxFF / hybrid / experiment: Young's 4.43 ± 0.23 / 3.69 ± 0.19 / 3.82 ± 0.20 / 3.41 ± 0.27; bulk 4.35 ± 0.12 / 3.02 ± 0.08 / 3.11 ± 0.10 / 2.86 ± 0.23; shear 1.66 ± 0.08 / 1.34 ± 0.10 / 1.42 ± 0.06 / 1.15 ± 0.13; Poisson 0.36 / 0.35 / 0.35 / 0.36.

**CAPS relevance**
- Not a vitrimer validation source, but the DGEBA/TETA numbers (Tg 419.35 K experimental, moduli above, 70 to 85% crosslinking range quoted) are a ready benchmark for the CAPS fix bond/react epoxy-amine builder with PCFF or CVFF-class FFs.
- The switching idea (cheap classical FF between reactive bursts, periodic short reactive probes) maps onto a CAPS workflow of alternating fix bond/react cycles with plain MD relaxation; CAPS could reproduce the Tg-from-density-kink and modulus table first, with the hybrid ML layer as an optional later feature.

---

### 4. Lin2026_JCP_penetrant_chemrxiv

**Local file:** Lin2026_JCP_penetrant_chemrxiv.pdf (ChemRxiv preprint, posted 15 January 2026)

**Citation (preprint):** Lin, M.-H.; Lin, T.-W.; Sing, C. E. Molecular Simulation Study of Penetrant Diffusion in Vitrimer Networks. ChemRxiv 2026, preprint.

**DOI (preprint):** 10.26434/chemrxiv-2026-b6vsp. Crossref: title, three authors, posted-content 2026-01-15 match. Crossref relation "is-preprint-of" 10.1063/5.0322457.

**Published version:** Lin, M.-H.; Lin, T.-W.; Sing, C. E. Molecular simulation study of penetrant diffusion in vitrimer networks. J. Chem. Phys. 2026, 164 (14), 144908. DOI 10.1063/5.0322457. Crossref: title, authors, JCP vol 164, issue 14, article 144908, published 2026-04-09. Verified as published. Extraction below is from the preprint text; the published version may differ [unverified].

**Chemistry**
- Model of dynamic poly(n-butyl acrylate) (nBA) networks (mapped to experiments by Evans and coworkers). Linear chains with randomly placed reactive beads (A), 3-bead crosslinker molecules with two reactive ends (B), free excess reactive beads act as exchange initiators.
- Exchange: generic associative (Sciortino three-body). Catalyst not modelled.
- Rubber relevance: PnBA is an acrylic rubber; T near Tg and rubbery states probed. No ENR/SBR/NR/PDMS.

**Model level and force field**
- Coarse-grained semiflexible bead-spring (Kremer-Grest-like but harmonic bonds, not FENE). One bead = one nBA monomer. T* = 485 K; τ* = (mσ²/kBT*)^(1/2).
- Permanent bonds: harmonic, k = 2000, r0 = 1.
- Bending: kθ(1 − cos θ), kθ = 1.52 (matches PnBA Kuhn length).
- Dynamic A-B bond: generalized LJ with ε_AB = 30, cutoff r_cut = 1.2 (exponents unreadable in the text conversion; [unverified]).
- Three-body term: u_ijk = λ ε_AB ν3(r_ij) ν3(r_ik); ν3 = 1 for r ≤ r_min, −u_AB(r)/ε_AB beyond. Barrier βΔE_BE = ε_AB(λ − 1).
- Non-bonded: generalized LJ 4ε[((dα+dβ)/2r)^(2m) − ((dα+dβ)/2r)^m]. A-A and B-B: m = 100, ε = 30, cutoff 2^(1/100) (purely repulsive). Pairs with C: m = 6, ε = 1, cutoff 2.5. Penetrant P with others: m = 6, ε = 1, cutoff 2.5 × (dα + dP)/2.

**Exchange method**
- Sciortino three-body potential in MD (no MC, no acceptance rule). Associative, crosslink count conserved.
- λ = 1.1 (fast) and 1.25 (slow).
- Temperatures (reduced): T = 0.57 (near Tg), 0.63, 0.70.
- Number of excess reactive beads N_free = 40 (default) or 200 (about five per chain).

**System size and run length**
- nc = 40 chains of Nb = 30 beads; nx crosslinkers with nc·NA = 2nx. Crosslink density f_cross = nx/(nc·Nb − nx) = 0.053, 0.111, 0.200, 0.333.
- One penetrant of diameter dP = 2, 3, 4 (4.25 in SI).
- Protocol: NPT at P = 0 with A-B bonds treated as permanent harmonic and A-A, B-B as standard LJ; then short NVT at the equilibrated volume with dynamic bonds switched on; then NVT until the penetrant is Fickian, at least 3 × 10^? steps (exponent unreadable [unverified]). dt = 0.0005 τ*. Nosé-Hoover thermostat and barostat. Three independent trajectories per state.

**Software:** LAMMPS (fix names not stated; three-body exchange potential implementation not named).

**Properties and how they were analysed**
- Bond lifetime τ_bond: bond autocorrelation C_bond(t) = ⟨H(t):H(0)⟩/⟨H(0):H(0)⟩ of the dynamic-bond connectivity matrix, fitted with KWW stretched exponential exp[−(t/τ_bond)^βKWW].
- Segmental relaxation τ_α,K: second Legendre autocorrelation of Kuhn vectors (three consecutive bonds); τ at decay to 1/e.
- Penetrant hopping time τ_α,P: self-intermediate scattering function Fs(q*, t) at the first peak of the network S(q); τ at decay to 0.1 (text conversion reads "0.121", taken as 0.1 with reference 21 [unverified]).
- Diffusion D_P: Einstein relation, lim MSD/6t in the Fickian regime.
- Confinement parameter C = dP / mesh size; analysis uses X ∝ exp(−C²)/C (Cai et al.), with exp(−bC) (Dell and Schweizer) in SI. Plots of log(C·D_P) vs C², 1/τ_α,P vs C², and D_P/D_P,linear vs τ_α,K/τ_bond.
- Tg: T ≈ 0.57 taken as near Tg from earlier work (ref. 59), not measured here.
- Stress relaxation, modulus, healing: not computed.

**Key quantitative results**
- Confinement table C (dP = 2 / 3 / 4): f = 0.053: 0.46 / 0.69 / 0.92; 0.111: 0.58 / 0.86 / 1.15; 0.200: 0.66 / 0.98 / 1.31; 0.333: 0.73 / 1.10 / 1.47.
- For C² < 0.5 the exchange rate has little effect on D_P; for C² from about 0.48 to 1.21 and above, λ = 1.1 gives clearly higher D_P than λ = 1.25. Enhancement in 1/τ_α,P is much weaker than in D_P, so exchange mainly relaxes mesh confinement.
- N_free 40 to 200 lowers τ_bond by about one order of magnitude with τ_α,K nearly unchanged; with λ = 1.25 this gives more than half an order of magnitude higher D at C² ≈ 2.2.
- D_P/D_P,linear collapses vs τ_α,K/τ_bond; plateau onset at τ_α,K/τ_bond ≈ 1 (C = 0.92, 1.31, 1.47 shown).

**CAPS relevance**
- Clean CG vitrimer reference for transport: reproduce C_bond(t) KWW fits, Kuhn P2 relaxation and penetrant Fs(q*,t) for one state point (f = 0.111, λ = 1.25, T = 0.63) as a validation of a CAPS three-body exchange implementation.
- The analysis set (bond autocorrelation with KWW, P2 vector autocorrelation, Fs(q,t), MSD to D) is generic and worth adding to CAPS Analyze as vitrimer dynamics tools regardless of the model.

---

### 5. Sharma2026_RSCAdv_CF_disulfide

**Local file:** Sharma2026_RSCAdv_CF_disulfide.pdf

**Citation:** Sharma, H.; Kaushik, N.; Liu, S.; Anand, R.; Sahoo, N. G.; Yun, G. J.; Rana, S. Reprocessable carbon fiber composites via disulfide-exchange epoxy vitrimers: experimental and molecular simulation insights. RSC Adv. 2026, 16 (19), 17480-17496.

**DOI:** 10.1039/d6ra00851h. Crossref: title, all seven authors, RSC Advances vol 16, issue 19, pages 17480-17496, 2026 match. Verified.

**Chemistry**
- Exchange: aromatic disulfide metathesis, catalyst-free; plus hydroxy-ester linkages from the acid-epoxy cure (transesterification not discussed as active).
- Experimental matrix: DGEBA (340.41 g/mol) cured with 2,2'-dithiodibenzoic acid (DTBA, 306.36 g/mol) at 1, 2, 3 wt% (EPD-1/2/3) with 1 wt% imidazole initiator; cure 80 °C 2 h, 120 °C 2 h, 150 °C 5 h. Composite: 4 plies UD carbon fibre (3K, 200 g/m²) by VARIM, pre-cure 120 °C 2 h, post-cure 150 °C 5 h.
- Simulation matrix: DGEBA + AFD (aminophenyl disulfide), not DTBA. The paper does not comment on this mismatch; Fig. 1 labels "R2 reactive O atoms in AFD" (as printed).
- Rubber relevance: disulfide exchange is the chemistry of sulfur-vulcanised NR/SBR and polysulfide rubbers; the matrix itself is a glassy epoxy.

**Model level and force field**
- Atomistic, PCFF. Non-bonded cutoff 9.5 Å; Ewald accuracy 0.0001 kcal/mol; velocity Verlet; 1 fs.
- Carbon fibre represented by 4 parallel carbon sheets in a 47 × 49 × 129 Å box, fully periodic; epoxy placed in the voids.

**Crosslinking / exchange method**
- Crosslinking: script-based iterative method; 126 DGEBA and AFD molecules (as printed) in top and bottom voids; bonds made between reactive pairs within 8 Å; final crosslink density 90%.
- "Bond Exchange Reaction (BER) equilibration": disulfide bonds iteratively broken and reformed by a cyclic script during MD (details in ref. 50: Park et al., Polymer 2020, 206, 122862). Barrier, attempt frequency, acceptance rule, BER cutoff, number of cycles: not stated.
- Annealing: NPT 300 K 50 ps, heat to 600 K and 50 ps, cool to 300 K and 50 ps.
- Creep MD: PCFF cannot break bonds, and no BER is said to be active during creep runs; creep is therefore without exchange [as described; the paper still attributes T dependence to disulfide exchange].

**System size and run length**
- Atom count not stated. Equilibration 50 ps per stage. Tensile at 5 × 10⁸ s⁻¹ along z (perpendicular to sheets), x, y at 0 bar, strain range 0 to 0.08. Creep 80 MPa along z, NPT, 130 to 170 °C; creep run length not stated.

**Software:** LAMMPS (fix names not stated); PCFF.

**Properties and how they were analysed**
- Density: 1.27 g/cm³ (MD) vs 1.25 g/cm³ (experiment).
- Young's modulus: linear fit of virial stress vs engineering strain; 2.74 GPa (MD) vs 2.69 GPa (experiment, composite).
- Creep MD: strain vs time at 130, 140, 150, 160, 170 °C (qualitative).
- Experimental stress relaxation: 1% strain (TMA, three-point bending), 60 to 90 °C in 10 °C steps; τ = time to reach 1/e (36.7%) of initial stress; no stretched-exponential fit. τ = 480 s at 60 °C, 120 s at 90 °C.
- Arrhenius: ln τ = Ea/RT − ln A; linear fit y = 4763.9x − 11.95 reported, Ea = 48 kJ/mol. Check: a slope of 4763.9 K gives 4763.9 × 8.314 = 39.6 kJ/mol, and the two end points (480 s at 333.15 K, 120 s at 363.15 K) give about 46 kJ/mol, so the printed fit line and the printed Ea do not match exactly.
- Tg: tan δ peak (TMA 0.5 Hz, 0.1% strain, 3 °C/min); EPD-2 92 °C in text, 90.3 °C in Table 2 (TMA), 64 °C in Table 2 (DSC) and 64 °C in Table 3. Composite Tg 130 °C.
- Creep (matrix): 0.3 MPa for 5 min, 20 min recovery, 70 to 110 °C in 10 °C steps. Composite creep 130 to 170 °C: ≈0.5% strain after 5 min at 130 °C, ≈14% at 170 °C (stress not stated for composite).
- Healing efficiency: razor scratch, 80 °C for 20 and 40 min, optical microscope; η = (original width − healed width)/original width × 100%. EPD-1 84%, EPD-2 91%, EPD-3 86%.
- Others: Td5% (EPD-2 396 °C), gel content 98.3%, swelling 0.33%, storage modulus 55 GPa at 40 °C (TMA, as printed), flexural strength 18.5 MPa (Table 2), composite tensile strength 290 MPa, thermoforming at 140 °C, recovery at 100 °C, DMF degradation at 70 °C in 4 h.

**Key quantitative results**
- Ea = 48 kJ/mol (stress relaxation), τ 480 s to 120 s over 60 to 90 °C.
- Healing 91% (EPD-2, scratch width).
- MD E = 2.74 GPa vs 2.69 GPa; MD density 1.27 vs 1.25 g/cm³.

**CAPS relevance**
- Directly reproducible in CAPS with PCFF + fix bond/react: DGEBA network with an aromatic disulfide hardener, density and Young's modulus as targets (1.25 g/cm³, 2.69 GPa). Preferably model the actual DTBA acid-epoxy chemistry, which the paper did not do.
- Disulfide exchange could be added as a fix bond/react template (S-S + S-S to two new S-S) with a Metropolis or Arrhenius acceptance using Ea ≈ 46-48 kJ/mol; the experimental τ(T) is a weak validation target only, because MD times are far shorter.
- Weak paper for protocol detail (BER script, atom counts, run lengths missing); use as experimental data, not as a method template.

---

### 6. Teotonico2024_PolymChem_boronic

**Local file:** Teotonico2024_PolymChem_boronic.pdf

**Citation:** Teotonico, J.; Mantione, D.; Ballester-Bayarri, L.; Ximenis, M.; Sardon, H.; Ballard, N.; Ruipérez, F. A combined computational and experimental study of metathesis and nucleophile-mediated exchange mechanisms in boronic ester-containing vitrimers. Polym. Chem. 2024, 15 (3), 181-192.

**DOI:** 10.1039/d3py01047c. Crossref: title, all seven authors, Polymer Chemistry vol 15, issue 3, pages 181-192, 2024 match. Verified.

**Chemistry**
- Exchange: boronic ester (6-membered dioxaborinane) metathesis, either direct or nucleophile (MeOH, water, pendant OH) mediated. Associative.
- Model compounds: B1 2-(4-methylphenyl)-1,3,2-dioxaborinane, B2 5,5-dimethyl-2-phenyl-1,3,2-dioxaborinane, B4 product, B5/B6 crosslinker precursors.
- Network: n-butyl acrylate photopolymerised with 3 mol% B6 (dioxaborinane dimethacrylate/styrenic crosslinker), 1 wt% Darocure 1175; HEMA as internal nucleophilic catalyst at 0 (VT0), 0.6 (VT1), 1.5 (VT2) mol% relative to BA. Control: ethylene glycol dimethacrylate (non-dynamic).
- Rubber relevance: PBA network is an elastomer (Tg ≈ −41 °C). Boronic ester exchange is used in rubbers elsewhere, but no NR/SBR/ENR/PDMS here.

**Model level:** DFT only (no MD).
- ωB97XD / 6-31+G(d) geometries and frequencies (minima / TS check, ZPVE, thermal corrections to H and G at 298.15 K, harmonic approximation); single points ωB97XD / 6-311++G(2df,2p). Gaussian 16. Solvation model: not stated (gas phase assumed [unverified]).

**Exchange energetics (Gibbs free energies, kcal/mol)**
- MeOH-mediated path (B1 + B2 + MeOH, second MeOH assisting proton transfers): TS1 15.74, TS2 20.96, INT2 14.74, TS3 26.97, TS4 30.83, INT4 25.12, TS5 34.27, INT5 31.08, TS6 37.99, TS7 22.67, INT7 20.69; rate-determining step = final intramolecular metathesis of the cyclic dimer, barrier 17.68 kcal/mol (as stated; relative to INT7).
- Direct metathesis: [2+2] TS1 29.15 (rate-determining), INT1 17.58, TS2 32.58.

**Experimental kinetics / properties and how they were analysed**
- Model exchange: B1 + B2 equimolar, 0.5 M in CDCl3 at 60 °C, ¹H NMR of –OCH2– signals at 3.80 ppm (B2) and 3.82 ppm (B4). After 4 h without catalyst: B2/B4 = 9/1 in solution, 6.1/3.9 in bulk. MeOH 5 and 10 mol% little effect; 20 and 50 mol% much faster. Effective rate coefficient assuming equal forward and back rate coefficients (Fig. S31) extrapolates to a non-zero value at zero MeOH (direct metathesis contributes). 50% water (saturated) also faster.
- Hydrolytic stability of B4 at 60 °C, 4 h in DMSO (Table 1): MeOH 0.10 to 5 eq gives ≤1% to 3% alcoholysis; water 0.10 to 5 eq gives 2% to 28% hydrolysis.
- Stress relaxation: shear G(t), 8 mm plates, linear regime, several temperatures (representative 160 °C); fitted with a stretched exponential (KWW); τ vs 1/T Arrhenius (Fig. 6). Numerical activation energies: not stated in the text (Fig. 6 only) [unverified]. HEMA lowers both τ and Ea; control does not relax.
- Tg: DSC inflection, first heating at 20 °C/min, about −41 °C for all networks. DMTA single transition and rubbery plateau. Td5% > 335 °C, stable to 310 °C (TGA N2).
- Frequency sweeps: VT2 reaches a G′/G″ crossover at low frequency.
- Reprocessing: hot press 150 °C, 10 min, 8 MPa (methods also give 300 bar); properties retained.
- Healing efficiency, MD, Tv: not stated.

**Key quantitative results**
- DFT rate-determining barriers: 17.68 kcal/mol (MeOH-mediated) vs 29.15 kcal/mol (direct metathesis).
- Bulk exchange near equilibrium in 4 h at 60 °C without catalyst.
- Network Tg ≈ −41 °C; Td5% > 335 °C.

**CAPS relevance**
- Supplies barrier inputs for a boronic ester exchange rule: 17.68 vs 29.15 kcal/mol (73.97 vs 121.96 kJ/mol, check: × 4.184) can set the acceptance probability of a fix bond/react or MC-swap template with and without a pendant hydroxyl nearby (catalyst-on-demand rule).
- First reproducible target: an all-atom PBA network (OPLS-AA or GAFF; boron parameters would need care) built with fix bond/react at 3 mol% crosslinker; check Tg ≈ −41 °C (232 K) by volume-temperature scan. Exchange kinetics are too slow for direct MD comparison.

---

### 7. Yang2016_RSCAdv_welding

**Local file:** Yang2016_RSCAdv_welding.pdf

**Citation:** Yang, H.; Yu, K.; Mu, X.; Wei, Y.; Guo, Y.; Qi, H. J. Molecular dynamics studying on welding behavior in thermosetting polymers due to bond exchange reactions. RSC Adv. 2016, 6 (27), 22476-22487.

**DOI:** 10.1039/c5ra26128g. Crossref: title, all six authors, RSC Advances vol 6, issue 27, pages 22476-22487, 2016 match. Verified.

**Chemistry**
- Exchange: transesterification (hydroxyl + ester) in a Leibler-type epoxy-acid vitrimer. Catalyst not modelled explicitly (folded into the Arrhenius rate). Associative.
- Monomers: DGEBA + tricarballylic acid (triacid crosslinker).
- Rubber relevance: none (glassy epoxy, Tg 316 to 334 K).

**Model level and force field**
- Atomistic, PCFF. Nosé-Hoover thermostat and barostat. Timestep: not stated.

**Crosslinking and exchange method**
- Crosslinking: iterative distance criterion (Wu and Xu; Varshney et al. scheme as adapted in ref. 28, Yang et al. Soft Matter 2015, 6305-6317): active atoms closer than 4.5 Å are bonded, then minimisation and MD relaxation; 32 iterations to 97.22% conversion at 450 K, 1 atm.
- BER: in each iteration, if an active atom (hydroxyl O) is within a cutoff of an existing ester bond it is attached; the new and old bond energies are compared and the exchange is accepted if the new bond has lower energy. BER cutoff value: not stated. For temperature studies an extra MC test: probability from k = A exp(−Ea/RT) with Ea = 68.18 kJ/mol normalised to 1 at 450 K (rate at 420 K = 27% of 450 K); accept if a random number is smaller. Source of the 68.18 kJ/mol value: not stated [unverified].
- Each BER iteration corresponds to about 0.1 ns of MD (110 iterations = 11 ns, as printed).
- Temperatures: 450 K (main), 420 K.

**System size and run length**
- Each film: 90 DGEBA + 60 tricarballylic acid; two films total 10,860 atoms. Film height 36 Å held by repulsive walls; x, y = 42 Å × 46 Å after crosslinking; periodic in x, y only.
- Initial gap 1 Å between films set by rigid walls; middle wall removed after the first cross-interface BER. Constant volume, 300 BER iterations.
- 360 active units; 2850 BERs in 10 ns at 450 K, giving τ_BER = 1.26 ns (average time per BER per active unit).

**Software:** LAMMPS (fix names not stated; BER done by external iteration scripts).

**Properties and how they were analysed**
- Tg: change in thermal expansion coefficient from MD temperature scan; 330 K (welding films), 334 K (stress relaxation system), 316 K (80.96% DOP), 320 K (S2 lower crosslink density).
- Stress relaxation: stretch 30% in z at 1 × 10⁹ s⁻¹, 450 K, then BER iterations; σzz per iteration fitted with a single exponential, τ = 4.2 ns. Activation energy from τ at 420 and 450 K via shift factor aT = exp[(Ea/R)(1/(T+273) − 1/(Tr+273))] (as printed), Ea = 67.67 kJ/mol, close to the 68.18 kJ/mol input, so free-volume change has negligible effect.
- Welding strength: uniaxial tension of welded system at 450 K, 1 × 10⁹ s⁻¹ to 80% strain; initial Young's modulus and yield stress vs normalised welding time t/τ_BER; healed when stress-strain equals the fresh network.
- Interface bonds: number of bonds crossing the interface, normalised by saturation value, vs t/τ; compared with the Yu et al. 2015 lattice model n(t,0) = n0 + Σ Cn exp(−t/τn), N = 2 layers, K = 1/τ, τ = 4.2 ns.
- Penetration depth of active atoms vs time; compared with Stukalin hopping diffusion (DH ≈ R0²/τb, Dp = 2/3 DH, W ≈ (DH t)^(1/2)) and Leibler tube/Rouse model with a = 14.4 Å, S = 5, p = τb/(τb + τo).
- Sticker lifetimes per 10 iterations (Table 1): closed 745 to 971 ps, free 621 to 729 ps, p 0.54 to 0.58.
- Crosslink density ν = ρ/Mc: S1 1.04 g/cm³, Mc 340, 0.003; S2 1.02 g/cm³, Mc 680, 0.0015 (units printed as "cm3", presumably mol/cm³).

**Key quantitative results**
- Full welding at 8.73 τ_BER (110 iterations, 11 ns) at 450 K; modulus recovers to 600 MPa.
- 80.96% DOP needs 13.49 τ_BER (170 iterations) vs 8.73 τ_BER for 97.22%.
- S2 (lower crosslink density) welds in 6.35 τ_BER (80 iterations) vs S1 8.73 τ_BER; more BERs because active O : active C = 2:1 vs 1:1.
- 420 K vs 450 K: 1242 vs 2850 BERs after 100 iterations (7.94 τ_BER); penetration 20.97 Å vs 31.04 Å at 6.35 τ_BER; yield stress at 450 K nearly twice that at 420 K; full recovery needs 24.6 τ_BER (310 iterations) at 420 K vs 8.73 τ_BER at 450 K.
- Gap 0.5 Å: first BER after 0.2 τ_BER; 2 Å: 2.58 τ_BER (exponential increase).
- MD penetration depth matches Stukalin's theory; Leibler's model predicts less.

**CAPS relevance**
- The closest template to a CAPS atomistic workflow: PCFF + LAMMPS, iterative distance-based crosslinking, then an MC-gated transesterification exchange with an Arrhenius acceptance, all achievable with fix bond/react (ester + hydroxyl template, with a custom rate via the "prob" keyword scaled by exp(−Ea/RT)).
- First protocol to reproduce: the 10,860-atom two-film weld at 450 K; targets are τ_BER ≈ 1.26 ns, stress relaxation τ ≈ 4.2 ns with single-exponential fit, Ea recovered ≈ 67.67 kJ/mol from 420/450 K, and modulus recovery to 600 MPa after ≈11 ns. Bonds-across-interface count and penetration depth are easy analysis additions.
- Caveat: the energy-comparison acceptance rule has no detailed balance; a CAPS version should use a proper Metropolis criterion and report the difference.

---

### Cross-paper notes for CAPS

- Exchange implementations seen: Sciortino three-body potential in MD (Ciarella 2019, Lin 2026), iterative distance + bond-energy + Arrhenius MC script (Yang 2016), unspecified cyclic break/reform script (Sharma 2026), ReaxFF bond-order chemistry with ML-gated FF switching (Hou 2026, curing only), DFT barriers only (Teotonico 2024), experiment only (Arbe 2023).
- Relaxation fits: KWW for bond autocorrelation (Lin) and experimental G(t) (Teotonico); single exponential (Yang); 1/e time (Sharma); stress autocorrelation plateau of frozen topology (Ciarella).
- Tv: only Arbe 2023 measures it (XRD mesoscale expansion kink, and η0 = 10¹² Pa·s extrapolation). No simulation paper in this group computes Tv.
- Healing metrics: Gxz, Gyz vs Gxy plateau and bridge counts (Ciarella); modulus and yield recovery, bonds across interface, penetration depth (Yang); scratch width ratio (Sharma, experimental).

## Group C

---

### 1. Branham-Ferrari, Cheng, Sokolov, Simmons: reaction/diffusion competition in vitrimers

- **Local file:** BranhamFerrari2026_NatPhys_arXiv2509.19496.pdf (the text is the arXiv preprint with its SI)
- **Citation:** M. R. Branham-Ferrari, S. Cheng, A. P. Sokolov, D. S. Simmons, "Reaction-diffusion competition drives anomalous relaxation of vitrimers", Nature Physics (2026), published online 2026-10-15. Crossref lists no volume, issue or pages yet. Preprint: arXiv:2509.19496, "Reaction/Diffusion Competition Drives Anomalous Relaxation of Vitrimers", posted 2025-09-23.
- **DOI:** 10.1038/s41567-026-03468-0. **Crossref:** title, all four authors, journal and online date (2026-10-15) match; Crossref has no volume or pages yet. A Crossref bibliographic search found this published version. The arXiv DOI 10.48550/arXiv.2509.19496 is not in Crossref (DataCite DOI), but doi.org resolves it and the arXiv API returns the same title. **Status: published.** I used the preprint text; the published version may differ [unverified].
- **Chemistry:** The simulation is generic: telechelic strands whose end beads form exchangeable bonds, with no catalyst. The experiments use imine (transimination) vitrimers: amine-terminated PPG (poly(propylene glycol) bis(2-aminopropyl ether), Mn = 230, 400 and 4000 g/mol) crosslinked with benzene-1,3,5-tricarbaldehyde. Amine:aldehyde molar ratio is 1:0.75, and free NH2 drives the exchange. Free NH2 is given as "0.25*2 per molecule" in Methods and "0.5 per molecule" in the SI. The data of Soman et al. (boronic ester ethylene vitrimers) are re-analysed. Rubber relevance: none directly. PPG networks are rubbery soft networks, but there is no ENR, SBR, NR, disulfide or PDMS.
- **Model level:** Bead-spring Kremer-Grest. Force field: LJ 12-6 with σLJ = 1 and ε = 1. For non-bonded pairs between beads in exchangeable bonds, σLJ is reduced to 0.985, 0.96 or 0.92. Bonds are FENE with K = 30, R0 = 1.50, ε = 1.0 and σ = 1.0 for all bonds, exchangeable bonds included. Theory: a reaction-diffusion model and an activated-state model.
- **Exchange method:**
  - Algorithm: hybrid MC/MD using LAMMPS bond/swap (Metropolis MC bond exchange), restricted by group to the telechelic end beads. The exchange is associative (partner swap; the number of bonds is conserved).
  - Attempts: each core attempts one swap every Nevery steps, with Nevery/Ncores = 1 timestep/attempt. This gives a temperature-invariant attempt rate of about 15 τ/(attempt/bond) for 3000 exchangeable bonds.
  - Cutoff: attempts only on bonds whose beads are within 1.5 σ.
  - Acceptance: Metropolis.
  - Barrier: theory gives activation energy E* = 3.24 (LJ units) with l* = 0.995 for the standard parameters. It is obtained by minimising 2E_LJ(l) + 2E_FENE(l). It matches Arrhenius fits of τex, and smaller σLJ lowers it (Fig. 1c).
  - Activated-state geometry: bond length 0.99, non-bonded distance 0.99. Ground state: bond 0.96, non-bonded 1.22 in Fig. 1a, but the text says 1.12 σ. The largest displacement needed is 0.13 σ, against a Debye-Waller rattling displacement of 0.195 σ.
  - Catalyst: not modelled.
- **System and run:**
  - Size: 3,000 primary strands, each with 8 non-reactive beads plus one exchangeable bead at each end. That is 10 beads per strand and about 30,000 beads [inference].
  - Setup: one exchangeable bond per strand. The system starts as 3000 rings.
  - Timestep: 0.005 τLJ (τLJ is about 1 ps).
  - Thermostat and barostat: Nose-Hoover, damping 2 τLJ, P = 0.
  - Equilibration: 10^5 τLJ at T = 1.5 for the standard system, then 10^3 τLJ after reducing σLJ.
  - Cooling: PRESQ stepwise quench, annealing at each T for more than 10 τα, then long production runs. The temperature list is not stated numerically.
- **Software:** LAMMPS (version not stated; bond/swap). Analysis with AMDAT, the Simmons group toolkit. Experiments: AR2000ex rheometer, Novocontrol Alpha and Agilent 4291B dielectric analysers.
- **Properties and analysis:**
  - τα: self-intermediate scattering function at q = 7.07 (near the first S(q) peak). The long-time part is fitted with a stretched exponential, and τα is the time where the fit reaches 0.2.
  - Debye-Waller factor: MSD at about 1 τLJ.
  - τex: counts all exchange events, recursive ones included.
  - Bond decorrelation: fraction of exchangeable beads that still have their bonded partner from time 0. It is fitted with a stretched exponential at high T. At low T it is fitted with the two-step form f(t) = 1/2[1 + e^(-t/τ)] exp[-(t/τdiff)^β] (Eq. 3, SI Eq. 14).
  - τbond comes from these fits. Model: τbond ≈ τex + τdiff ≈ τex + cτα (Eq. 4).
  - τdiff: time at which segment MSD reaches 10^(1/2) σ (Fig. 3 caption). At low T, τbond converges to the time for segments to move about 2 σ.
  - Arrhenius fits of τex over the high-T range and the full range (Fig. 1c).
  - Experiments: τbond = 1/ωcross from the G'/G'' crossover (SAOS, 10^-3 to 10^2 rad/s, 8 mm plates, N2, dried at 423 K for 2 h). τα from BDS (10^-2 to 10^9 Hz) fitted with Havriliak-Negami functions (two HN terms for PPG4000, one for PPG230 and PPG400) plus a dc-conductivity term, with VFT fits for τα. Not computed: Tv, Tg, modulus, creep, healing.
- **Key results:** τex is Arrhenius and its barrier is below the apparent α barrier. Bond relaxation crosses over from Arrhenius at high T to super-Arrhenius at low T. A two-step decorrelation appears when τex < τα. Classification: displacement ratio x = r*ex/r*α; x ≥ 1 is type I (Eq. 1 holds), x < 1 is type II (reaction-diffusion). Fits of experimental PPG230 and PPG400 data to Eq. 4 give c ~ 10^7 to 10^10. For the simulation, Rouse theory predicts c ~ n^2.
- **CAPS relevance:** This is the simplest exchange protocol to reproduce first in the CAPS bead-spring engine. Use KG strands with FENE K = 30 and R0 = 1.5, an MC bond swap restricted to end beads within 1.5 σ, and Metropolis acceptance. Validation targets: E* = 3.24 ε for σLJ = 1, the trend of the barrier with σLJ, and the bond-autocorrelation analysis (fraction of original partners) with the two-step fit. It is also a useful test of an exported LAMMPS `fix bond/swap` deck against CAPS's own engine.

---

### 2. Golbabaei, Garner, Zhang: ML with cascade representation for vitrimer healing efficiency

- **Local file:** Golbabaei2026_DigitalDiscovery.pdf
- **Citation:** M. H. Golbabaei, H. Garner, N. Zhang, "A machine learning framework with cascade molecular representation for predicting self-healing efficiency in vitrimers", Digital Discovery 5 (10), 4413-4433 (2026).
- **DOI:** 10.1039/d6dd00171h. **Crossref:** title, three authors, journal, vol 5, issue 10, pages 4413-4433 and year 2026 all match.
- **Chemistry:** The dataset covers mostly epoxy-based vitrimers: transesterification, imine, disulfide, vinylogous urethane and so on, plus some PU. Features are epoxy, hardener, catalyst, additive, hardener ratio, catalyst wt%, additive wt%, Tg, healing temperature and healing time. Experimental validation used DGEBA cured with 2-aminophenyl disulfide (2-AFD) at a 1:1 molar ratio (epoxy:hardener), with no catalyst, cured at 160 °C for up to 6 h. Rubber relevance: none. This is aromatic disulfide epoxy, the same family as papers 5 and 6.
- **Model level:** Machine learning only, with no simulation. A feed-forward ANN works on RDKit descriptors.
- **Exchange method:** Not modelled.
- **Data and model:**
  - Literature screening: 75 studies screened. 9 were excluded for non-comparable healing methods and 18 for incompatible chemistry, leaving 48 papers and 181 vitrimer systems.
  - Imputation: KNN imputer, R2 0.99 on deliberately removed values. The other imputers: iterative RF 0.96, iterative GB 0.94, mean 0.17.
  - Preprocessing: outliers removed, under-represented ranges duplicated in the training data only, then StandardScaler and Yeo-Johnson PowerTransformer. Features with |PCC| > 0.8 were dropped, and the transforms C*cata and t*heal = -√theal were applied.
  - Descriptors: 208 RDKit descriptors (RDKit 2022.09.5). The "cascade" averages precursor descriptors recursively, up to depth 5.
  - PCA: 5 components per material, then 15 components overall.
  - Network: 15 inputs; hidden layers 256-128-64-32-16-8-4-2 with Leaky ReLU; linear output.
  - Training: RMSprop, learning rate 0.001, MAE loss, batch 24, at most 9000 epochs, early-stopping patience 300 (Table 4 lists 2027 epochs as "optimized"), seed 123, TensorFlow 2.16.1 / Keras. A web interface was built with Flask.
  - Test set: 41 samples, each with at least one constituent unseen in training.
- **Software:** RDKit, TensorFlow/Keras, scikit-learn-style baselines, Flask.
- **Properties and analysis:**
  - Target: healing efficiency taken from the literature as reported.
  - Their own experiment: η = (1 - w/w0) × 100, where w0 and w are the cut widths before and after healing. Images were taken with a KEYENCE VR-3100 at 80×, and widths measured in ImageJ.
  - Tg: DSC (TA Q20, N2, 10 to 200 °C at 5 °C/min), taken at the onset of the transition.
  - Metrics: MAE, signed mean error and R2.
- **Key results:**
  - Test set: MAE 2.76 %, R2 0.96, SME -1.33 %. Direct-SMILES subset (24 samples): MAE 1.91 %. Cascade-derived subset (17 samples): MAE 3.96 %, SME -2.71 %. Validation MAE was 6.70 %.
  - Baselines: KNN (k = 1) MAE 0.67 %, R2 0.987. RF MAE 2.29 %, R2 0.960. GBR MAE 2.360 %, R2 0.933. SVR MAE 5.66 %, R2 0.776. Dummy MAE 21.01 %, R2 -0.017. The authors prefer the ANN for its physically monotonic time response.
  - DGEBA/2-AFD experiment: Tg ≈ 72 °C. Healing at 90 °C for 1 h took the cut from 44.979 to 2.802 (units printed as "mm"; µm is more likely [unverified]). Measured 93.77 % against 96.88 % predicted.
  - Yang et al. (DGEBA/adipic acid/TBD/MWNT at 180 °C): 46.77 % after 3 min and 56.25 % after 1 h. Predicted: 45.74 % and 55.08 %.
  - Five external formulations: errors from -4.38 to +5.53 pp.
- **CAPS relevance:** There is no simulation protocol to reproduce. The paper gives experimental targets for the DGEBA/aromatic-disulfide system that papers 5 and 6 simulate: Tg ≈ 72 °C (DSC onset) for DGEBA:2-AFD at 1:1, and about 94 % healing at 90 °C after 1 h. A CAPS Tg run (volume-T kink) on a PCFF or GAFF DGEBA/2-AFD network could be checked against 72 °C, keeping in mind the cooling-rate offset.

---

### 3. Irigoyen, Matxain, Ruipérez: combined DFT and MD protocol for disulfide self-healing materials

- **Local file:** Irigoyen2022_FrontMater_disulfide.pdf
- **Citation:** M. Irigoyen, J. M. Matxain, F. Ruipérez, "Combined DFT and MD Simulation Protocol to Characterize Self-Healing Properties in Disulfide-Containing Materials: Polyurethanes and Polymethacrylates as Case Studies", Frontiers in Materials 9, 859482 (2022).
- **DOI:** 10.3389/fmats.2022.859482. **Crossref:** title, three authors, journal, vol 9, article 859482 and date 2022-03-21 all match.
- **Chemistry:**
  - Exchange: aromatic disulfide exchange by the (2+1) radical-mediated mechanism. The S-S bond cleaves into sulfenyl radicals, and a radical attacks a disulfide through a three-membered transition state. There is no catalyst; the reaction is catalyst-free.
  - Polyurethanes (Kim et al. 2017, "transparent elastomers"): aromatic-disulfide hard segments named after their isocyanate. IP is isophorone diisocyanate, HM is 4,4'-methylenebis(cyclohexyl isocyanate) and M is 4,4'-methylenebis(phenyl isocyanate).
  - Polymethacrylate networks (Takahashi et al. 2017): BTA uses the BiTEMPS bis(2,2,6,6-tetramethylpiperidin-1-yl) disulfide crosslinker, a sulfenamide-type S-N bond. ADSA uses a dialkyl disulfide diacrylate.
  - Rubber relevance: no NR, SBR, ENR or PDMS. The PUs are elastomers with aromatic disulfide crosslinks, the chemistry of disulfide self-healing rubbers (Rekondo-type poly(urea-urethane)). The aliphatic disulfide (ADSA) is a proxy for the dialkyl di- and polysulfide crosslinks of sulfur-vulcanised rubber. The BDE and ΔG values below are directly usable for disulfide rubber models.
- **Model level:**
  - DFT on one- and two-chain fragment models: ωB97XD/6-31+G(d,p) geometries and frequencies (ZPVE and 298 K thermal corrections to H and G), with single points at 6-311++G(2df,2p). Gas phase.
  - Ab initio BOMD: PBE with a DZP basis and RI, Nose-Hoover at 298 K. The other temperatures listed below were also run.
  - Classical all-atom MD: "GAFF ff14SB" AMBER force field (wording as printed) with RESP charges from HF/6-31G* ESP, in vacuum, NVT.
- **Exchange method:** Exchange is not simulated in MD. The paper instead computes three descriptors:
  - ρ, the probability of radical formation: the fraction of QMD steps with S-S longer than 2.3 Å (Eq. 1).
  - Rate constant: from the exchange barrier by Eyring/TST, k = (kBT/h) e^(-ΔG/RT) (Eq. 2).
  - ω, mobility: from the S-S radial distribution function, ω = I_i/(I_i + I_ii), with a reacting region R ≤ 4.5 Å and a neighbouring region 4.5 < R < 20 Å; R > 20 Å is ignored (Eq. 3). In the MD analysis the integration ranges are 3 to 4.5 Å and 4.5 to 20 Å.
  - ρMD: from the MD RDF, ρMD = I_II/(I_II + I_III), with region II 2.20 to 3.00 Å and region III 1.80 to 3.00 Å (Eq. 4).
  - Mechanism: associative radical-mediated exchange, with a dissociative first step (homolysis).
  - Barriers (DFT ΔG):
    - IP: 12.02 kcal/mol.
    - Earlier value for a para-urethane aromatic disulfide: 12.21 kcal/mol (Irigoyen 2019a).
    - Sulfenamides, used for BTA: 16-18 kcal/mol (Ruipérez 2018).
    - Barriers for HM and M are assumed similar to IP. ADSA's barrier was not computed.
  - BDE (kcal/mol, open/closed conformation):
    - IP: -/59.31
    - HM: -/47.34
    - M: 51.20/65.65
    - BTA: 33.77/-
    - ADSA: 58.19/63.13
    - Earlier reference values: about 48 for diphenyl disulfides, 33 for sulfenamides, 60 for aliphatic disulfides.
  - S-S bond length Re: 2.050 to 2.200 Å (Table 1).
  - Interchain interaction ΔH: M -49.48, BTA -31.71, ADSA -31.41 kcal/mol, "around 30 kcal/mol" in general.
  - Temperatures: IP, HM and M at 25, 40, 60 and 80 °C; BTA and ADSA at 25 and 120 °C.
- **System size and run:**
  - QMD: 40,000 a.u. (9.651 ps), timestep 40 a.u. The paper prints "(1.93 ps)"; it must be about 1 fs, since 40 a.u. ≈ 0.97 fs [inference; printed value inconsistent]. TURBOMOLE 6.6.
  - Classical MD: 16 chains placed at random, each holding two disulfide units with a short chain between them. Density reproduces experiment after equilibration.
  - MD protocol: minimisation of 100,000 steps (50,000 steepest descent, 50,000 conjugate gradient), then heating from 0 to 300 K over 200 ps, 2 fs timestep, SHAKE on X-H bonds. Production length and atom count are not stated.
- **Software:** Gaussian 16 (DFT), TURBOMOLE 6.6 (QMD), AMBER 14 with LEaP/AmberTools (MD).
- **Properties and analysis:**
  - BDE: energy of the disulfide minus the two sulfenyl fragments, so non-covalent interactions are included.
  - ρ from QMD (Table 3) and ρMD from the RDF (Table 4).
  - ω from the RDF (Table 5).
  - TST rate, from ΔG only.
  - Healing efficiency is not computed; it is compared qualitatively with experiment. No Tg, Tv, stress relaxation, modulus or MSD.
- **Key results:**
  - QMD ρ at 25 °C: IP 0.0117, HM 0.0000, M 0.0000, BTA 0.4750, ADSA 0.0000.
  - QMD ρ at higher T: BTA 0.2165 and ADSA 0.0000 at 120 °C. IP falls with T (0.0012 at 80 °C), which the authors put down to conformational bias.
  - ωave: HM 0.0806, IP 0.0613, M 0.0361. BTA 0.0313 at 25 °C and 0.0699 at 120 °C (the text quotes ωave 0.0506). ADSA 0.0004 at 25 °C and 0.0008 at 120 °C (text ωave 0.0006).
  - ρMD rises with T for IP: 0.0030, 0.0036, 0.0042, 0.0055 at 25, 40, 60, 80 °C.
  - The protocol qualitatively ranks experimental healing IP > HM > M, and BTA heals while ADSA does not.
  - Experimental numbers quoted (validation targets): IP heals a scratch fully in 2 h at room temperature, and in 30, 5 and 3 min at 40, 60 and 80 °C. A cut IP film lifted a 5 kg load after 4 h. BTA healed at 120 °C for 8 h, recovering 85 % of strength and 92 % of strain at break. The PU Tg is quoted as 120-130 °C.
  - The authors conclude the exchange barrier does not discriminate between materials; radical formation and mobility do.
- **CAPS relevance:**
  - Cheapest protocol to reproduce first: an all-atom GAFF (or OPLS-AA) amorphous cell of disulfide-containing chains in CAPS, then the S-S RDF integrals ω (3 to 4.5 Å over 4.5 to 20 Å) and ρMD (2.20 to 3.00 Å over 1.80 to 3.00 Å). No reactive chemistry is needed.
  - Parameters for a later CAPS `fix bond/react` disulfide-exchange template, or for a temperature-dependent acceptance in rubber vulcanisate models: ΔG values of 12.02 kcal/mol (aromatic) and 16-18 kcal/mol (sulfenamide), and BDEs of about 33 to 66 kcal/mol.
  - A DFT re-check in CAPS's own QM pipeline (ωB97XD at the same basis, ORCA or Gaussian) on the IP fragment would validate the 12.02 kcal/mol barrier.
  - The MD runs were in vacuum, so CAPS should use periodic bulk cells instead.

---

### 4. Oyarzún and Mognetti: efficient sampling of reversible cross-linking polymers

- **Local file:** Oyarzun2018_JCP_arXiv1803.00126.pdf (the text is arXiv v1, 28 Feb 2018)
- **Citation:** B. Oyarzún, B. M. Mognetti, "Efficient sampling of reversible cross-linking polymers: Self-assembly of single-chain polymeric nanoparticles", J. Chem. Phys. 148 (11), 114110 (2018).
- **DOI:** 10.1063/1.5020158. **Crossref:** title, both authors, journal, vol 148, issue 11, article 114110 and date 2018-03-20 all match. **Status:** arXiv:1803.00126 is published as this JCP paper. Singh 2020 cites it with the same DOI.
- **Chemistry:** Generic reversible dimerisation 2A ⇌ A2 between reactive monomers on one chain, giving single-chain polymeric nanoparticles (SCPNs) with supramolecular or reversible links such as H-bonds, metal-ligand and π stacking. There is no specific chemistry and no catalyst. Rubber relevance: none. It is a single chain in good or ideal solvent.
- **Model level:** Coarse-grained chains with fixed bond length σ. Monomers are either ideal (freely jointed, U = 0) or interact through a WCA potential (βU = 4[(σ/r)^12 - (σ/r)^6] + 1 for r ≤ 2^(1/6) σ). Bead-spring, but with rigid bonds and no FENE.
- **Exchange method:**
  - Algorithm: pure MC, with no dynamics. A bind/unbind move uses topological configurational-bias MC: the chain section inside a sphere Vm of radius Rm around the reactive site α is regrown with fixed-end Rosenbluth growth, guided by freely-jointed end-to-end distributions p(r, L).
  - Acceptance: Pacc(f→b) = min[1, Keq · nf · nf,Vm/(nν+1) · p(rSα, LSβ+1) p(rβE, LβE)/p(rSE, LSE) · WΓb/WΓf], and the reverse rule (Eqs. 10-11).
  - Chemistry parameter: Keq = qA2/(2 qA^2) (in units of σ^3), which can come from experiment or from atomistic/QM calculations (Eq. A5).
  - Other moves: single and double pivot, regrowth, and displacement of reacted complexes. Each move type is chosen with equal probability in each MC cycle.
  - Rm only affects efficiency (tested 3, 5, 7, 10, 15 and 100 σ); acceptance peaks at Rm = 5 to 10.
  - Bond type: reversible bonds that break and form (a dissociative, bind/unbind type, not a swap). The irreversible limit is emulated by rejecting unbinding.
  - No energy barrier or temperature beyond kBT units.
- **System and run:**
  - Chains: NC = 100 monomers, functionalisation f = NR/NC from 0 to 1, with regular (REG) or random (RAN) placement.
  - Keq values: 10 (Fig. 5), 100 (Table I, 20-monomer chain with 4 sites), 1000. Fig. 6 shows others [unverified].
  - Runs: 5 × 10^6 MC cycles of production after an equally long equilibration, about 5 h on a 2.3 GHz processor, averaged over 50 independent runs.
- **Software:** The authors' own code (not named).
- **Properties and analysis:**
  - ⟨Rg²⟩ from gyration-tensor eigenvalues. Reference values: ideal R0² = σ²NC/6 = 16.67 σ²; WCA R0² = 57.20 σ².
  - Fraction of bound sites ⟨nν⟩, loop-length distributions, and average connectivity maps ⟨νij⟩.
  - Validation against exact ideal-chain results: p0/pL = 1/(Keq nL), and two-loop partition functions (Appendix C).
  - No stress relaxation, Tv, Tg or healing.
- **Key results:**
  - Table I (theory vs simulation): p0,0 0.1350 vs 0.1347 ± 0.0004; p5,5 0.5343 vs 0.5335 ± 0.0008; p20,5 0.2001 vs 0.2000 ± 0.0005; p10,10 0.1306 vs 0.1318 ± 0.0004.
  - Ideal chains: Rg shrinks until a plateau near f ≈ 0.4.
  - WCA REG chains at high Keq: Rg is re-entrant in f (the f = 0.5 chain is larger than the f = 0.25 and f = 1 chains). The effect is absent for RAN chains and in the irreversible limit, where results depend on Rm, a sign of non-equilibrium.
- **CAPS relevance:** Low priority for a materials builder. The algorithm is useful if CAPS adds an MC topology sampler for reversible crosslinks, for example equilibrating supramolecular or ionic rubber networks before MD. The first check to reproduce is the exact Table I / Fig. 5 ideal-chain loop probabilities for any new bind/unbind move.

---

### 5. Singh and Sundararaghavan: MD self-healing of vitrimers with dynamic cross-linking

- **Local file:** Singh2020_CPL_selfheal_arXiv2005.02307.pdf (the text is arXiv v2, 24 Dec 2020)
- **Citation:** G. Singh, V. Sundararaghavan, "Modeling self-healing behavior of vitrimers using molecular dynamics with dynamic cross-linking capability", Chemical Physics Letters 760, 137966 (2020).
- **DOI:** 10.1016/j.cplett.2020.137966. **Crossref:** title, both authors, journal, vol 760, article 137966 and date 2020-12 all match. **Status:** arXiv:2005.02307 is published; the preprint states this itself.
- **Chemistry:**
  - Epoxy vitrimer: DGEBA cured with "4-aminophenyl disulfide (AFD)", aromatic disulfide exchange, no catalyst.
  - Isomer: the text also says the amine is at the ortho position, which would be 2-AFD, so the isomer is inconsistent.
  - Stoichiometry: epoxy:hardener 2:1.
  - Healing case: a (12,12) SWCNT is pulled out of the matrix.
  - Rubber relevance: aromatic disulfide exchange (same mechanism family as disulfide rubbers), but the matrix is a glassy epoxy.
- **Model level:** All-atom. CVFF with pair, bond, angle, dihedral and improper terms. LJ plus Coulomb with a cutoff (cutoff value and long-range method not stated).
- **Exchange method:**
  - Algorithm: LAMMPS fix bond/react (Gissinger) with pre- and post-reaction templates and a map, two-step reactions (templates in SI Fig. S6). Associative S-S metathesis that conserves the number of bonds.
  - Cutoff: 4.12 Å, twice the S-S bond length of 2.06 Å.
  - Probability: no barrier. The reaction probability is an empirical logistic function of T: p(T) = 1/[exp(-a(T - Tv) - ln(pTv/(1 - pTv))) + 1], with a = 2pTv/(w(1 - pTv)), w = 20 K, pTv = 0.9, and p(Tv - w) = 0.1.
  - Tv: assumed equal to Tg, 403 K.
  - Rates are deliberately accelerated relative to experiment.
  - Curing (also fix bond/react): C-N cutoff 3.5 Å for primary and 5.0 Å for secondary amine, probability 1, reaching up to 95 % crosslinking.
- **System and run:**
  - Size: 8×8×8 replicas of the unit cell, giving 2,048 DGEBA and 1,024 AFD, 68,608 atoms. CNT case: (12,12) SWCNT of 91 Å length.
  - Packing: compressed to 1.0 g/cm³ at 300 K and 1 atm (NPT).
  - Annealing: cycles between 1 K and 600 K, NPT at 1 bar, until density converges to 1.18 g/cc for the neat system and 1.19 g/cm³ with the CNT, at 1 K.
  - Heating scan: 150 ps at each T from 200 to 550 K in 5 K steps.
  - Healing: 503 K (Tv + 100 K) for 250 ps, then cool to 300 K.
  - Timestep, thermostat details and LAMMPS version are not stated.
- **Software:** LAMMPS (fix bond/react; version not stated).
- **Properties and analysis:**
  - Number of reactions vs T.
  - Volume-T curve normalised by V0 at 200 K, static vs dynamic S-S. Linear CTE = ΔV/(3V0).
  - Uniaxial tension at 2 × 10^8 s^-1 at several T.
  - Stiffness: linear regression within ε = 4 %, averaged over three directions.
  - Healing: stress-strain of damaged vs healed samples at 300 K, and modulus vs T for pristine, damaged and healed samples. No scalar healing efficiency is defined.
  - Tv is not detected; it is imposed through p(T).
  - No stress relaxation, MSD or Arrhenius analysis.
- **Key results:**
  - Static CTE below Tv: 41.1 ± 2.7 µK^-1 (literature range 45.0-65 µK^-1).
  - Dynamic S-S model: higher rubbery CTE above Tv.
  - Young's modulus at 300 K: 4.5 ± 0.21 GPa, against 2.6 GPa experimentally.
  - Exchange count: 53 % fewer exchanges in 150 ps for Tv = 200 K than for Tv = 403 K.
  - After healing, the modulus recovers to the pristine value below Tv, and the CNT hole fills.
- **CAPS relevance:**
  - This is the most direct template for CAPS's `fix bond/react` export. A DGEBA + AFD network can be built and cured in CAPS (PCFF, CVFF or GAFF), then exported with two templates: amine-epoxy curing at 3.5 and 5.0 Å, and S-S exchange at 4.12 Å with T-dependent `prob`.
  - Validation targets: density 1.18 g/cc, CTE 41.1 µK^-1 below Tv, E 4.5 GPa at 2 × 10^8 s^-1, and the shape of the reaction count vs T.
  - The logistic p(T) rule is ad hoc. A CAPS version could replace it with an Arrhenius probability built from a DFT barrier (paper 3: ΔG about 12 kcal/mol for aromatic disulfides).

---

### 6. Uyama, Kishimoto, Oya, Murashima, Koyanagi: disulfide bond healing in epoxy vitrimers (QC + MD)

- **Local file:** Uyama2026_Polymers_disulfide.pdf
- **Citation:** T. Uyama, N. Kishimoto, Y. Oya, T. Murashima, J. Koyanagi, "Molecular Mechanism of Disulfide Bond Healing and Network Repair in Epoxy Vitrimers Revealed by Quantum Chemical and Molecular Dynamics Simulations", Polymers 18 (7), 861 (2026).
- **DOI:** 10.3390/polym18070861. **Crossref:** title, five authors, journal, vol 18, issue 7, article 861 and date 2026-03-31 all match.
- **Chemistry:** DGEBA cured with 4-aminophenyl disulfide (APDS) at APDS:DGEBA = 1:2 (equimolar functional groups), with no catalyst. The focus is post-scission radical recombination (healing of broken S-S bonds), not exchange. Rubber relevance: aromatic disulfide, as in paper 3. The mechanochemical S-S scission and recombination idea carries over to disulfide-crosslinked rubbers, but the system is a glassy epoxy.
- **Model level:**
  - DFT: B3LYP/6-31+G(d) driven by GRRM20 (Global Reaction Route Mapping) through Gaussian 16.
  - Mulliken bond populations.
  - All-atom MD with PCFF.
- **Exchange or reaction method:**
  - DFT exchange: an APDS + APDS-fragment radical, run through GRRM. It finds one intermediate and two transition states involving three APDS units. Barriers are 150-180 kJ/mol, against a quoted S-S BDE of about 190 kJ/mol. Mechanism: associative radical exchange.
  - DFT recombination of two fragments, relative to the reformed S-S at 0 kJ/mol: S-C(ortho) adduct about 105 kJ/mol; S-C(para) adduct about 106 kJ/mol. The non-bonded state is about 71 kJ/mol below the para adduct and about 60 kJ/mol below the ortho adduct. The authors say these are energy differences, not barriers.
  - Mulliken populations: intrinsic S-C σ about 0.29; new S-C(ortho/para) 0.20-0.22; π-stacking 0.02-0.04.
  - MD curing: LAMMPS fix bond/react, distance rule 5 Å between epoxy and amine, with reaction-relaxation cycles under NPT. Final conversion about 80 %.
  - MD bond breaking: during tension, any bond except aromatic C-C and X-H bonds is deleted when it exceeds 0.64 × the GRRM maximum covalent length. The factor was calibrated so peak stress falls at a strain of about 0.1. Tested range "0.625 to 0.65 in increments of 0.05" (as printed).
  - MD healing: fix bond/react with five rule-based reactions, each triggered when the pair reaches the QC maximum covalent distance:
    1. direct S-S;
    2. S-C(ortho);
    3. S-C(para);
    4. ortho → S-S;
    5. para → ortho.
  - The healing reactions are directional (para → ortho → S-S), with no activation energies and no temperature dependence. The bond-breaking threshold and the reaction-probability values for bond/react are not stated.
- **System and run:**
  - Size: about 27,000 atoms. Three independent networks × tension along x, y and z gives 9 runs.
  - Mixing equilibration: NPT at 300 K and 1 atm.
  - Tension: strain rate 5.0 × 10^8 /s with lateral cross-section fixed, at 300 K. Damage is introduced by straining to 200 %.
  - Healing: damaged cell relaxed with NPT for 5 ns (300 K, 1 atm), then recombination; it saturates within about 5 ns.
  - Timestep and LAMMPS version are not stated.
- **Software:** LAMMPS (fix bond/react), GRRM20, Gaussian 16.
- **Properties and analysis:**
  - Stress-strain: Young's modulus and tensile strength.
  - Fraction of S-S bonds surviving vs strain.
  - Repair: the fraction of broken S-S bonds reformed vs time, split by pathway.
  - Stress-strain of undamaged, damaged (no healing) and healed cells. Healing efficiency is not expressed as a single ratio, only as an ordering.
  - No Tg, Tv, stress relaxation or MSD.
- **Key results:**
  - Mechanics: E ≈ 5.0 GPa, tensile strength ≈ 180 MPa (170-220 MPa across thresholds), strain at peak ≈ 0.1. Only S-S bonds break under tension.
  - Repair: 78 % of broken S-S bonds reform and 22 % stay unpaired. Of the reformed bonds, 64 % come from direct S-S and 36 % via the ortho pathway.
  - Modulus and strength order: undamaged > healed > damaged.
- **CAPS relevance:**
  - It shows a PCFF + fix bond/react workflow with bond scission, the same force-field family CAPS exports. Two reproducible elements:
    1. a selective bond-breaking criterion (scaled maximum covalent length per bond type);
    2. a multi-template bond/react healing step.
  - Validation targets: E ≈ 5.0 GPa and σmax ≈ 180 MPa at 5 × 10^8 /s with a fixed cross-section; S-S being the only broken bond type; about 78 % reformation within about 5 ns.
  - Check how CAPS handles PCFF class II cross terms when bonds are deleted. The paper does not state how they did it.

---

### 7. Zheng, Biswal, Guo, Thakolkaran, Kokane, Varshney, Kumar, Vashisth: MD-informed ML for vitrimer Tg

- **Local file:** Zheng2025_DigitalDiscovery_TgML.pdf
- **Citation:** Y. Zheng, A. K. Biswal, Y. Guo, P. Thakolkaran, Y. Kokane, V. Varshney, S. Kumar, A. Vashisth, "Toward sustainable polymer design: a molecular dynamics-informed machine learning approach for vitrimers", Digital Discovery 4 (9), 2559-2569 (2025).
- **DOI:** 10.1039/d5dd00239g. **Crossref:** title, eight authors, journal, vol 4, issue 9, pages 2559-2569 and year 2025 all match.
- **Chemistry:** Bifunctional transesterification vitrimers made from dicarboxylic acid + diepoxide at an equimolar ratio, with TBD catalyst in the experiments.
  - Hypothetical library: 1 million vitrimers sampled from 2.5 billion combinations of 50,000 acids × 50,000 epoxides (ZINC15).
  - Synthesizable library: 259 vitrimers from 37 acids × 7 epoxides sold by Sigma-Aldrich.
  - Synthesised: DGEBA + D,L-malic acid (MA) and DGEBA + 2,6-naphthalenedicarboxylic acid (NA), with 5 mol% TBD in 5 mL acetone, stirred at 140-170 °C. MA cured at 140 °C for 5 h; NA cured at 280 °C for 6 h under 1 MPa.
  - Rubber relevance: none.
- **Model level:** ML trained on all-atom MD Tg labels. The labels (8424 vitrimers) come from the earlier work (Zheng et al., Adv. Sci. 2025, 12, 2411385 [unverified]), calibrated against experiment with a Gaussian process. The MD force field, cell size, cooling rate and Tg fit are in the SI and the earlier paper. **Not stated in the main text.**
- **Exchange method:** None. Bond exchange is not simulated; only Tg is predicted.
- **Data and models:**
  - Split: 90/10 train/test.
  - Representations:
    - Morgan fingerprints: 2048 bit, radius 3, frequency-based, top 200 substructures reduced to 195.
    - 166 RDKit descriptors.
    - 620 Mordred descriptors.
    - Mol2vec: 300-dimensional.
    - SMILES of the repeating unit (epoxide opened by the acid).
    - Molecular graphs.
  - Models: LASSO, RF, SVR, XGBoost, FFNN, GNN (Hickey et al.) and a transformer (TransPolymer), 22 models in all. Five-fold cross-validation, plus SHAP analysis.
- **Software:** scikit-learn, XGBoost, PyTorch, RDKit, Mordred, Mol2vec. Code: github.com/vashisth-lab/VitrimerScreening, Zenodo 10.5281/zenodo.16391076 (not checked [unverified]).
- **Properties and analysis:** Tg only.
  - ML metrics: R2 and RMSE on the test set.
  - Screening: about 100 top and bottom candidates re-simulated with MD + GP calibration.
  - Experimental Tg: DSC at 5 °C/min.
  - Healing: only qualitative microscopy of cut and healed samples (temperature and time not stated).
- **Key results:**
  - Best single model: LASSO_Mordred, R2 0.76, RMSE 15.72 K.
  - Best ensemble (XGBoost_fp + XGBoost_Mordred + GNN + transformer): R2 0.78, RMSE 15.19 K.
  - MD-validated discovered vitrimers: low Tg 218-237 K, high Tg 510-549 K.
  - Experiments: MA vitrimer Tg 348 K (75 °C) against 364 K predicted; NA vitrimer Tg 395 K (122 °C) against 390 K predicted. The NA value is higher than any bifunctional transesterification vitrimer they found in the literature.
  - SHAP: rotatable bonds and FractionCSP3 lower Tg; aromatic, O/N and H-bonding substructures and VSA_EState2 raise it.
- **CAPS relevance:**
  - Targets for a CAPS Tg workflow (crosslinked network build, cooling scan, bilinear fit) on DGEBA + diacid transesterification networks: experimental Tg 348 K (MA) and 395 K (NA).
  - The open labelled dataset (8424 MD Tg values) could benchmark CAPS's own high-throughput Tg pipeline once the SI protocol (force field, cooling rate) is read [unverified].
  - The paper has no exchange chemistry.

---

### Cross-paper notes for CAPS

- **Same DGEBA/aromatic disulfide system in three papers.** Singh (CVFF, 68,608 atoms) and Uyama (PCFF, about 27,000 atoms) simulate it, and Golbabaei gives experimental Tg ≈ 72 °C and healing data for DGEBA/2-AFD at 1:1. This system is a natural first CAPS vitrimer benchmark: build and cure with `fix bond/react`, run a Tg scan and tension, then S-S exchange or recombination templates.
- **Barriers available for an Arrhenius bond/react probability.** Aromatic disulfide exchange ΔG 12.02 kcal/mol (ωB97XD, Irigoyen). Sulfenamide 16-18 kcal/mol. The GRRM B3LYP exchange path gives 150-180 kJ/mol (Uyama), much higher because its pathway and reference differ. The two disagree by a factor of about 3 (12 kcal/mol is about 50 kJ/mol); CAPS should not mix them without its own calculation.
- **Coarse-grained first step.** The Branham-Ferrari bond/swap KG model (E* = 3.24 ε analytically) is the cheapest exact validation for a CAPS bead-spring exchange engine.

## Group D

---

### 1. Chen et al. 2023, Polymer (NR with hydrogen-bond network)

- **Local file:** Chen2023_Polymer_NR_Hbond.pdf
- **Citation:** Qionghai Chen, Wanhui Huang, Liqun Zhang, Li Xi, Jun Liu. "Fully atomistic molecular dynamics simulation of chemically modified natural rubber with hydrogen-bonding network." Polymer 284 (2023) 126284.
- **DOI:** 10.1016/j.polymer.2023.126284. **Crossref:** verified. Title, all five authors, journal, volume 284, article 126284 and year 2023 all match.
- **Chemistry:** No covalent exchange chemistry. Sacrificial hydrogen bonds only. cis-1,4-polyisoprene (synthetic NR model) is modified with maleic anhydride (MAn) and 3-amino-1,2,4-triazole (ATA) (MAn-ATA-PIP, after Liu et al. Macromolecules 2016). No crosslinker, since chemical crosslinks are deliberately left out "to isolate the effect of hydrogen bonds". No catalyst. **Rubber relevance:** high (NR / polyisoprene).
- **Model level and force field:** fully atomistic, COMPASS (the paper cites the phosphazene and azide COMPASS papers, refs 44 and 45). PBC in 3D. The paper does not state COMPASS parameters, the cutoff or the electrostatics method.
- **Exchange method:** no reactive exchange. Hydrogen bonds form and break physically under the force field. A hydrogen bond is counted when the H-O or H-N distance is below 3.5 (the text gives no unit; Angstrom is implied). Barrier: not stated. Swap frequency: not applicable. Associative vs dissociative: not applicable (physical bonds).
- **System size and run length:**
  - 50 repeat units per chain (m + n = 50), with modified fraction α = 0, 1, 2, 4, 6, 8 and 10%.
  - The α = 1% case says "25 molecular chains are modified ... other 25 are not", which implies 50 chains.
  - Total atom count: not stated.
  - Equilibration: 20 ns NPT at 400 K and 0.1 MPa, cooled to 298 K over 10 ns, then 20 ns NPT at 298 K and 0.1 MPa. Timestep 1 fs. Nose-Hoover thermostat and barostat, velocity Verlet.
  - Healing runs: 1, 10 and 30 ns NVT at 298, 373 and 498 K.
- **Software:** LAMMPS (version not stated), OVITO for analysis.
- **Properties and how each was analysed:**
  - *H-bond content:* distance distribution f(r) (Eq. 4) for H-O and H-N pairs. Bonds are split into intra-chain and inter-chain by molecule ID.
  - *Density and non-bonded energy:* reported as a function of α.
  - *MSD:* atomic MSD (Eq. 5). No diffusion coefficient is fitted.
  - *Tg:* method given in the SI only, so not stated in the main text. The paper only says Tg rises with α.
  - *Uniaxial tension:* constant-volume deformation. The box is stretched in z at a constant engineering strain rate while x and y shrink. The paper writes the rate as "10−8 s−1" (almost certainly a typo for 10^8 s^-1; it is reported here exactly as printed). Stress is σT = (−3Pzz + ΣPii)/2.
  - *Mechanism analysis:* ΔE_non-bonded versus strain, and the backbone bond orientation <P2(cos θ)> (Eq. 6). Cracks are located from atom-number profiles along normalised z.
  - *Self-healing:*
    - Starting system: α = 6%.
    - Gap: a repulsive wall of 2500 H atoms normal to z leaves a gap of about 3 Å. The wall is deleted, then the system heals under NVT.
    - Gap closure: tracked by the recovery ratios ρm/ρm,average and ρn/ρn,average in the gap.
    - Test: triaxial stretch along z, with x and y held fixed and stress taken as Pzz.
    - **Healing efficiency:** η = S_self-healing / S_original, the ratio of the areas under the stress-strain curves (a toughness ratio).
  - *Stress relaxation:*
    - G(t) from equilibrium MD using the Green-Kubo relation, G(t) = V/(kBT) <σαβ(t)σαβ(0)>, averaged over σxy, σxz and σyz.
    - Uses the multi-tau correlator (Ramirez et al. 2010).
    - No KWW or Maxwell fit. G' and G'' come from piecewise-linear analytic Fourier integration of G(t) (Adeyemi, Zhu, Xi 2022). η* = G*/ω, and η0 is read from the low-ω plateau.
  - *Not studied:* Tv, Arrhenius analysis, creep and welding strength.
- **Key numbers:**
  - Inter-chain H bonds reach 83% of all H bonds at α = 10%.
  - The best mechanical properties occur at α = 6%.
  - ρm recovers 60% after 1 ns at 298 K, and fully after 1 ns at 498 K.
  - η = 44.7% after 1 ns at 498 K and 82.3% after 30 ns at 498 K.
  - At low ω, G' ∝ ω² and G'' ∝ ω.
- **CAPS relevance:** CAPS already has COMPASS, so this NR system is a direct candidate. Two protocols are cheap to copy:
  - **Healing protocol:** wall-gap healing followed by a triaxial test, with η as a toughness ratio.
  - **Viscoelastic analysis:** the Green-Kubo G(t) with a multi-tau correlator, then piecewise-linear G'/G''.

  The first validation target could be the trend of Tg and η0 against α, with η = 44.7% (1 ns) and 82.3% (30 ns) at 498 K as check values. A 50 × 50-unit PIP system is small enough to run.

---

### 2. Hamzehlou and Ruipérez 2022, J. Polym. Sci. (vinylogous urethane / urea transamination, DFT)

- **Local file:** Hamzehlou2022_JPolymSci_VU_DFT.pdf
- **Citation:** Shaghayegh Hamzehlou, Fernando Ruipérez. "Computational study of the transamination reaction in vinylogous acyls: Paving the way to design vitrimers with controlled exchange kinetics." Journal of Polymer Science 60(13) (2022) 1988-1999.
- **DOI:** 10.1002/pol.20220099. **Crossref:** verified. Title, both authors, journal, volume 60, pages 1988-1999 and year 2022 match.
- **Chemistry:**
  - **Reaction:** transamination of vinylogous acyls with a free amine (associative, vitrimer chemistry).
  - **Model pair:** n-butyl vinylogous urethane + benzylamine.
  - **Variants:** X = O, S, NH, CH2 (urethane, thiourethane, urea, amide); Y = O or S (carbonyl versus thiocarbonyl); R1 = n-butyl or tert-butyl; R2 = CH2Ph (benzylamine) or Ph (aniline).
  - **Conditions / catalysis:** protic (neutral or acidic, protonated amine), aprotic (basic), and aprotic with the Lewis acid B(CF3)3.
  - **Rubber relevance:** none directly. This is small-molecule chemistry; vinylogous urethane elastomers exist, but the paper does not model a polymer.
- **Model level and force field:** DFT. ωB97XD/6-31+G(d) geometry optimisations (Def2SVPP on Sn). Frequencies at the same level give ZPVE and H and G at 298.15 K (harmonic). Single points use 6-311++G(2df,2p) and Def2TZVPP. PCM solvent is benzene.
- **Exchange method:**
  - **Approach:** static reaction-path (PES) calculations, with intermediates and transition states TS1 to TS4 for the protic path.
  - **Barriers:**
    - Protic path: the rate-limiting step is TS1 (proton transfer to the vinylic α-carbon). ΔG barrier is 9.29 kcal/mol and the electronic barrier ΔE0 is 14.3 kcal/mol. The paper compares this with the experimental Ea of 14.1 ± 1.4 kcal/mol (59 ± 6 kJ/mol; Denissen et al. 2015).
    - Aprotic path: the TS1 barrier is ΔE0 = 33.99 kcal/mol (experimental range quoted as 24-41 kcal/mol).
    - With B(CF3)3: ΔE0 = 19.93 kcal/mol.
    - Experimental values the paper quotes for comparison: DBTDL-catalysed 10.8 ± 1.7 kcal/mol; protic 17.5 ± 2.9 kcal/mol.
    - Vinylogous urea (X = NH): the rate-limiting step is TS3, with a 12.56 kcal/mol barrier and ΔE0 = 12.97 kcal/mol. The experimental value is 12.92 kcal/mol (Denissen 2018).
  - **Not applicable** (no MD in this paper): swap frequency, cutoffs and acceptance rule.
  - **Catalysts:** treated explicitly as reacting species (protonated benzylamine, B(CF3)3).
  - **Mechanism type:** associative.
  - **Text-extraction caveat:** the minus signs are lost in the text file (for example, the reactant complex reads "8.11 kcal/mol" and is almost certainly −8.11). Signed values should be checked against the PDF [unverified signs].
- **System size and run length:** gas-phase / PCM molecules only. Run length is not applicable.
- **Software:** Gaussian 16 Rev. C.01.
- **Properties:**
  - **Table 2:** ΔG of TS1 to TS4 for 7 systems.
  - **Table 3:** step barriers for the same 7 systems.
  - **Table 1:** H-bond distances in the reactant complex.
  - **Not studied:** stress relaxation, Tv, Tg, modulus and healing.
- **Key numbers:**
  - Protic TS1 barrier: ΔG = 9.29 kcal/mol; ΔE0 = 14.3 kcal/mol (vs experimental 14.1 ± 1.4).
  - With tert-butyl as the leaving amine, the TS3 barrier is 1.78 kcal/mol (8.79 with n-butyl).
  - With aniline as the free amine, the TS1 barrier drops to 6.28 kcal/mol, but TS2 rises to 12.61 kcal/mol.
- **CAPS relevance:** Use it as a source of barrier inputs, not as an MD protocol. The 14.3 kcal/mol (protic VU) and 12.97 kcal/mol (vinylogous urea) values can feed an Arrhenius or ΔGa acceptance term in a CAPS bond/react or MC-swap exchange model. CAPS's ORCA path could reproduce the protic TS1 barrier on the n-butyl VU + benzylamine pair as a QM validation step.

---

### 3. Karatrantos et al. 2024, Polymers (review)

- **Local file:** Karatrantos2024_Polymers_review.pdf
- **Citation:** Argyrios V. Karatrantos, Olivier Couture, Channya Hesse, Daniel F. Schmidt. "Molecular Simulation of Covalent Adaptable Networks and Vitrimers: A Review." Polymers 16(10) (2024) 1373.
- **DOI:** 10.3390/polym16101373. **Crossref:** verified. Title, all four authors, Polymers volume 16, article 1373 and year 2024 (published 11 May 2024) match.

#### How the review classifies the simulation methods

1. **Atomistic MD (Section 2.1, Table 1).** Fixed-charge force fields with the reaction step bolted on:
   - **Distance cutoff ("cut-off" method) plus an energy test on the new bond.** PCFF, DGEBA/tricarballylic acid, about 5430 atoms (Yang 2015). Also used for welding (Yang 2016) and creep anisotropy (Hanzon 2017). The cutoff starts small and is increased for stability.
   - **Temperature-dependent Monte Carlo probability added to the cutoff step.** Yang 2016.
   - **Pre/post templates with fix bond/react.** CVFF, DGEBA/AFD disulfide, 68,608 atoms (Singh and Sundararaghavan 2020; Singh et al. creep). Uses a logistic temperature-dependent probability p(T) with p(Tv) = 0.9 and an S-S cutoff of 4.12 Å (twice the 2.06 Å S-S bond length).
   - **Accelerated ReaxFF bond boost.** DGEBA/adipic acid healing (Kamble 2022).
   - **Diels-Alder dissociative CANs with van 't Hoff equilibrium.** PCFF (Sun 2020; Yang 2017).
   - **Imine-exchange polyimine networks.** PCFF (Zhang 2022).
   - **Disulfide polyurea self-healing.** PCFF (Zheng 2021).
   - **Solvent decomposition by transesterification.** PCFF (Sun 2019).
   - **ML force field from DFT with metadynamics and umbrella sampling.** Imine exchange free-energy profiles (Sun 2023).
   - **MD plus Gaussian-process calibration of Tg inside an inverse-design workflow.** Zheng et al. arXiv 2023.
2. **Atomistic vitrimer nanocomposites (Section 2.2).** Graphene oxide in DGEBA/AFD, PCFF (Park 2020). Tg falls from 383 to 363 K with GO.
3. **Coarse-grained MD (Section 2.3):**
   - **Patchy Kern-Frenkel particles** with event-driven MD and Wertheim theory (Smallenburg 2013).
   - **Sciortino three-body swap potential:**
     - Pair term: LJ-like 20-10 potential between reactive sites, with σd = 5 Å, εd = 4 kBT and cutoff 2.5 σd.
     - Three-body term: λ Σ ε V3(rij)V3(rik). λ = 1 gives free swapping; λ ≫ 1 gives no swapping.
     - Used for star networks: stress relaxation is G(t) from the stress autocorrelation via Green-Kubo, plus Geq. A swap barrier β∆Esw ≤ 10 gives a second relaxation; ≥ 50 freezes the topology.
     - Used for the Zhao et al. KG vitrimers: BER lifetime is Arrhenius in ∆Esw and T, the best mechanics occur at ∆Esw = 100, and healing is best at ∆Esw = 10.
     - Used for nanorod and nanoparticle composites (N = 60 chains).
   - **Dynamic-bond viscoelasticity** (Wang 2024): the healing time is Arrhenius in T.
4. **Monte Carlo (Section 3):**
   - Dimerisation MC with an equilibrium constant.
   - The Rao/Xia/Ni detailed-balance bias for multivalent swaps, with acc = min[1, v_u1/(v_u0+1) exp(−β∆G0)]. The hybrid form multiplies by exp(−β∆Ga).
   - Linker-mediated (Lei 2020) and thermo-gelling (Xia 2022) vitrimers.
5. **Hybrid MD-MC on Kremer-Grest bead-spring chains (Section 4):**
   - **Lineage:** from Pant-Theodorou end-bridging.
   - **Wu 2019:** two relaxation regimes.
   - **Perego and Khabaz 2020:**
     - Tv from the change in the thermal expansion coefficient.
     - Bond autocorrelation C(t) fitted with a KWW stretched exponential, average θ = 0.84; τav = (τ/θ)Γ(1/θ).
     - τav is Arrhenius for T ≤ 1.7 (LJ units).
     - η0 from NEMD extrapolated to zero shear rate, with η0 = 5.09 exp(1.94/T).
   - **Perego 2021-2024:** bond exchange rate effects, creep, J(t) ~ t^1/2.
   - **Stukalin 2013:** sticky-bond MC healing across an "ideal cut".
   - **Xia, Kalow, Olvera de la Cruz 2023:** PDMS-like side-chain vitrimer, Tv = 0.6 from the bond relaxation time and Fs(q,t), with a sticky Rouse model comparison.
   - **Ricarte and Shanbhag 2021:** SSR versus inhomogeneous Rouse models.
   - **Karmakar et al. (arXiv 2024):** entropic cohesion.
   - **Cui 2023, 2024:** associative versus dissociative CANs with a modified FENE sticky potential, USB = UFENE − UFENE(b0) − H.
6. **Gaps the review names:**
   - No atomistic study includes the catalyst.
   - Entangled associating networks have not been simulated.
   - Crosslink density, semicrystallinity, polydispersity and vdW contrast are unexplored.

#### Most important primary simulation papers cited (review reference number in brackets)

Every DOI below was checked on Crossref on 2026-10-22, and authors, title, journal, volume, pages and year match unless noted. The review itself gives no DOIs; it prints only [CrossRef] links.

| Ref | Paper | Model / method | DOI (Crossref) |
|---|---|---|---|
| [58] | Yang, Yu, Mu, Shi, Wei, Guo, Qi. Soft Matter 11 (2015) 6305-6317 | Atomistic PCFF, cutoff BER, transesterification | 10.1039/c5sm00942a (verified) |
| [59] | Yang, Yu, Mu, Wei, Guo, Qi. RSC Adv. 6 (2016) 22476-22487 | Atomistic welding by BER | 10.1039/c5ra26128g (verified) |
| [60] | Hanzon, He, Yang, Shi, Yu. Soft Matter 13 (2017) 7061-7073 | Creep-induced anisotropy, PCFF | 10.1039/c7sm01174a (verified) |
| [64] | Singh, Sundararaghavan. Chem. Phys. Lett. 760 (2020) 137966 | CVFF, fix bond/react disulfide, p(T) model | 10.1016/j.cplett.2020.137966 (verified) |
| [65] | Zheng, Yang, Sun, Zhang, Guo. Polymer 212 (2021) 123111 | Disulfide self-healing polyurea, PCFF | 10.1016/j.polymer.2020.123111 (verified) |
| [66] | Kamble, Vashisth, Yang, ..., Koratkar. Carbon 187 (2022) 108-114 | Accelerated ReaxFF healing, DGEBA/adipic acid | 10.1016/j.carbon.2021.10.078 (verified) |
| [69] | Sun, Wan, Shen, He, Zhou, Wang, Yang, Shi. Macromolecules 56 (2023) 9003-9013 | ML force field plus metadynamics, imine exchange | 10.1021/acs.macromol.3c01377 (verified) |
| [71] | Park, Kim, Jung, Krishnakumar, Rana, Yun. Polymer 206 (2020) 122862 | GO/vitrimer nanocomposite healing, PCFF | 10.1016/j.polymer.2020.122862 (verified) |
| [76] | Stukalin, Cai, Kumar, Leibler, Rubinstein. Macromolecules 46 (2013) 7525-7541 | KG sticky-bond MC healing | 10.1021/ma401111n (verified) |
| [90] | Smallenburg, Leibler, Sciortino. Phys. Rev. Lett. 111 (2013) 188002 | Patchy particle vitrimer | 10.1103/physrevlett.111.188002 (verified) |
| [96] | Ciarella, Sciortino, Ellenbroek. Phys. Rev. Lett. 121 (2018) 058003 | Three-body swap, star network stress relaxation | 10.1103/physrevlett.121.058003 (verified) |
| [99] | Sciortino. Eur. Phys. J. E 40 (2017) 3 | Three-body swap potential | 10.1140/epje/i2017-11496-5 (verified) |
| [101] | Rao, Xia, Ni. J. Chem. Phys. 160 (2024) 061102 | Detailed-balance MC swap | 10.1063/5.0186553 (verified) |
| [110] | Zhao, Wei, Fang, Gao, Yue, Zhang, Ganesan, Meng, Liu. Macromolecules 55 (2022) 1091-1103 | KG + three-body, ∆Esw effects on mechanics and reprocessing | 10.1021/acs.macromol.1c02034 (verified) |
| [111] | Wu, Li, Liu, Qian, Lu. Phys. Chem. Chem. Phys. 21 (2019) 13258-13267 | CG hybrid MD-MC kinetics | 10.1039/c9cp01766f (verified) |
| [112] | Zhao, Duan, Li, Chen, Yue, Zhang, Ganesan, Liu. Macromolecules 56 (2023) 9336-9349 | Multiscale dynamics, TTSP/TESP | 10.1021/acs.macromol.3c01893 (verified) |
| [123] | Perego, Khabaz. Macromolecules 53 (2020) 8406-8416 | Hybrid MD-MC, Tv from expansion, KWW C(t), Arrhenius η0 | 10.1021/acs.macromol.0c01423 (verified) |
| [124] | Perego, Khabaz. J. Polym. Sci. 59 (2021) 2590-2602 | Bond exchange rate vs mechanics | 10.1002/pol.20210411 (verified) |
| [129] | Xia, Kalow, Olvera de la Cruz. Macromolecules 56 (2023) 8080-8093 | PDMS-like vitrimer, sticky Rouse | 10.1021/acs.macromol.3c01366 (verified) |
| [155] | Cui, Jiang, Shao, Zhang, Yang, Tang. Macromolecules 56 (2023) 772-784 | Associative vs dissociative CAN viscoelasticity | 10.1021/acs.macromol.2c02122 (verified) |

**Not verified on Crossref (taken from the review only):**
- [67] Singh et al. creep is cited as arXiv 2304.05518. It is now published; see section 5.
- [70] Zheng et al. inverse design is cited as arXiv 2312.03690 [from review, unverified].
- [105] Ciarella and Ellenbroek, Coatings 9 (2019) 114 [from review, unverified].
- [125] Perego et al., Macromolecules 55 (2022) 7605-7613 [from review, unverified].
- [140] Ricarte and Shanbhag, Macromolecules 54 (2021) 3304-3320 [from review, unverified].
- [143] Perego and Khabaz, Macromol. Rapid Commun. 44 (2023) 2200313 [from review, unverified].

**CAPS relevance:** The review is the map. For CAPS, the cheapest atomistic route is the fix bond/react template route ([64], [67]), because CAPS already exports fix bond/react. The cheapest validation of the analysis code is the Perego-Khabaz CG protocol [123]:
- KWW fit of C(t);
- Arrhenius fit of τav and η0;
- Tv from the kink in the thermal expansion coefficient.

---

### 4. Rao, Xia, Ni 2024, J. Chem. Phys. (bond swap MC algorithm; preprint arXiv:2311.02983)

- **Local file:** Rao2024_JCP_bondswap_arXiv2311.02983.pdf (the text is arXiv v2, dated 17 January 2024).
- **Citation (published version):** Peilin Rao, Xiuyang Xia, Ran Ni. "A bond swap algorithm for simulating dynamically crosslinked polymers." The Journal of Chemical Physics 160 (2024) 061102.
- **DOI:** 10.1063/5.0186553. **Crossref:** found with a bibliographic query. Title, three authors, J. Chem. Phys. volume 160, article 061102 and year 2024 (14 Feb 2024) match. **Published:** yes. The local file is the arXiv preprint, so any differences from the published text were not checked [unverified].
- **Chemistry:** generic pivot (P) / residue (R) model. The worked chemical example is dioxaborolane metathesis (Röttger et al. 2017). Rubber relevance: none directly, though the multivalent residue model maps onto polysulfide or multi-arm crosslinkers by analogy [my inference].
- **Model level:** coarse-grained ideal particles. A bond is an infinitely deep square well of length rb. No other interactions in the tests (ideal gas).
- **Exchange method:** Metropolis-Hastings MC bond swap with a bias term for detailed balance.
  - **Step 1:** pick a random pivot P0 (probability 1/NP).
  - **Step 2:** pick a random attacking residue R1 with free valence within rb (1/Natt).
  - **Step 3:** pick a random leaving residue R0 bonded to P0 (1/Nlv).
  - **Acceptance:** acc = min[1, v_u1/(v_u0 + 1) exp(−β∆G0)] (Eq. 3).
  - **For hybrid MD-MC:** acc = exp(−β∆Ga) min[1, (p'/v'u)/(p/vu) exp(−β∆G0)] (Eq. 18), where ∆Ga is a tunable activation barrier.
  - **Not stated:** the swap frequency in MD; the paper only says swaps run "at a random time interval".
  - **Mechanism:** associative, since valence is conserved.
  - **Catalyst:** not modelled.
- **System size and run length:**
  - Monovalent test: NP = 400, box 20 rb, ∆G/ε = −1, µB = µC = −0.916, 5×10^5 MC moves to equilibrate and 1×10^8 to sample.
  - Multivalent linker test: vrB = 2, vrC = 3, NP = 2400, box 30 rb, µ = −2.133, ∆G/ε = −1; same move counts.
  - Binary chain test: vp = vr = 2, NP = 800, T = 1, ∆G = 0, box 20 rb, 5×10^5 + 5×10^8 moves.
  - Ensemble for all tests: grand canonical MC.
- **Software:** own code (not named). LAMMPS is cited only as a parallel-scheme reference.
- **Properties and analysis:**
  - Equilibrium species fractions against chemical-equilibrium theory: Eq. 5 for monovalent and Eq. 7 with binomial coefficients for multivalent.
  - Chain and ring length distributions against linear-condensation theory: Eq. 9 for chains, and Eq. 17 for rings with configurational entropy from a χ3 / regularised-gamma approximation plus combinatorial ln(1/i).
  - Not studied: stress relaxation, Tv, Tg, modulus and healing.
- **Key results:** simulation and theory agree for f_PRC against kBT/ε and µC − µB, for multivalent star fractions, and for chain densities at βµR = −2. Ring densities differ slightly because of the ∆S_conf approximation.
- **CAPS relevance:** This is the acceptance rule for a CAPS MC-swap module (CG or united-atom) that has to handle multivalent sites, such as sulfur crosslinks of rank > 1 or tri- and tetra-functional crosslinkers. The first validation is to reproduce the monovalent diatomic GCMC test (NP = 400, ∆G = −1 ε, box 20 rb) against Eq. 5. The test is analytic and cheap, so it makes a good unit test before any polymer run.

---

### 5. Singh, Sundararaghavan, Varshney 2024, Polymer (creep in vitrimers; preprint arXiv:2304.05518)

- **Local file:** Singh2024_Polymer_creep_arXiv2304.05518.pdf (the text is arXiv v1, dated 11 April 2023, with SI).
- **Citation (published version):** Gurmeet Singh, Vikas Varshney, Veera Sundararaghavan. "Understanding creep in vitrimers: Insights from molecular dynamics simulations." Polymer 313 (2024) 127667. Crossref lists the authors in this order; the arXiv version has Singh, Sundararaghavan, Varshney.
- **DOI:** 10.1016/j.polymer.2024.127667. **Crossref:** found with a bibliographic query. Title, authors, Polymer volume 313, article 127667 and year 2024 (November) match. **Published:** yes. The local text is the preprint, so any changes in the published version were not checked [unverified].
- **Chemistry:** disulfide metathesis (aromatic S-S exchange). DGEBA cured with 4-aminophenyl disulfide (AFD) at a 2:1 ratio. Catalyst-free; the reaction probability p stands in for catalyst amount and temperature. **Rubber relevance:** disulfide exchange chemistry, applicable to disulfide rubbers, but the matrix is a glassy epoxy.
- **Model level and force field:** all-atom CVFF (pair, bond, angle, dihedral and improper terms). LJ and Coulomb with a 12 Å cutoff. Timestep 1 fs. Nose-Hoover thermostat and barostat.
- **Exchange method:** LAMMPS fix bond/react (REACTER), topology-based.
  - **Curing:** primary and secondary amine templates, with N-C cutoffs of 3.5 Å and 5.0 Å. DGEBA is built with the epoxide ring pre-opened, so curing is a one-step reaction. Cured to 95%.
  - **S-S exchange:** pre/post templates of 44 atoms each. Triggered when sulfur atoms of two disulfide pairs come within 4.12 Å. Fixed reaction probability p = 0.0, 0.01, 0.1, 0.2 or 1.0.
  - **Creep loop:** NPT for 5 ps under applied stress, then NVE for 4 ps. Reactions are allowed only once per NVE step (set by the fix bond/react frequency). The loop repeats.
  - **Barrier:** not used (probability only).
  - **Mechanism:** associative.
- **System size and run length:**
  - 8×8×8 repeat cells give 1024 DGEBA + 512 AFD and 68,608 atoms. Box 8.72 × 8.72 × 8.72 nm³.
  - Anneal: cycles between 1 K and 600 K, 50 ps at each, NPT at 1 bar.
  - Density: 1.184 g/cm³ at 1 K and 1.159 g/cm³ at 300 K.
  - Creep: 600 K under 150, 200, 300 or 500 MPa uniaxial true stress (σyy = σ0, σxx = σzz = 0, shear left unspecified so the box stays orthogonal). Total creep time is shown in ns on the figures, but no number is stated in the text.
- **Software:** LAMMPS (version not stated), fix bond/react, OVITO with Python. A commercial builder was used only to draw the monomers.
- **Properties and analysis:**
  - *Creep:* stretch ratio λ = l/l0 against time. Primary, secondary and tertiary regimes are identified by eye. Findley's power law is mentioned but not fitted.
  - *Free volume:* alpha-shape surface mesh in OVITO with a 3.5 Å probe radius; void fraction V_void/V_cell × 100. The largest void volume is normalised by the cross-section.
  - *Bond orientation:* mean projection of the S-S vector on the load axis (Eq. 5) and its probability distribution over time.
  - *Reaction counting:* number of reactions per step.
  - *Not studied:* stress relaxation, Tv, Tg and healing. Tv from volume-temperature data is in the earlier CPL 2020 paper, not this one.
- **Key numbers:**
  - At 500 MPa and 600 K with p = 1.0, all three creep stages appear.
  - The void fraction heals to about 10% before tertiary creep starts.
  - Secondary creep rate is about constant for 150-300 MPa at p = 1.
  - At p = 0.01 the vitrimer behaves like the epoxy; p = 0.1 already raises creep strain markedly.
  - The S-S projection clusters near 2.0 Å in the epoxy.
- **CAPS relevance:** This is the most directly reproducible protocol for CAPS, because it is LAMMPS fix bond/react with CVFF. CAPS can export both the curing and the S-S templates. Validation targets:
  - **First:** the cured, annealed density (1.159 g/cm³ at 300 K, 95% cure).
  - **Next:** the alternating NPT 5 ps / NVE 4 ps creep loop at 600 K, 500 MPa, comparing λ(t) with p = 0 versus p = 1.
  - **Size:** 68,608 atoms is above the CAPS 15k-atom limit. A 4×4×4 replica (8,576 atoms) is my suggestion for a reduced test [my inference].

---

### 6. Wick, Peters, Li 2025, J. Polym. Sci. (shape memory vitrimer framework)

- **Local file:** Wick2024_JPolymSci_SMvitrimer.pdf (accepted September 2024, published in the January 2025 issue).
- **Citation:** Collin D. Wick, Andrew J. Peters, Guoqiang Li. "A Molecular Dynamics Modeling Framework for Shape Memory Vitrimers." Journal of Polymer Science 63(1) (2025) 9-19.
- **DOI:** 10.1002/pol.20240616. **Crossref:** verified. Title, three authors, journal, volume 63, pages 9-19 and year 2025 match. The file name says 2024 because it was received, accepted and published online in 2024.
- **Chemistry:** UV-cured bisphenol A glycerolate dimethacrylate (BPAGMA). The network forms by free-radical chain polymerisation of the methacrylate end groups. The exchange chemistry is transesterification between alcohol and ester groups (each monomer has 2 OH and 2 ester groups). The catalyst is not modelled. Rubber relevance: none (glassy shape memory thermoset).
- **Model level and force field:**
  - All-atom OPLS(-AA), parameters from LigParGen using SMILES.
  - 9 Å truncation with PPPM ("particle mesh Ewald" in the text).
  - Nose-Hoover thermostat and barostat.
  - Timestep: 0.25 fs during compression, then 0.5 fs.
  - New bonds are first given a weak 5 kcal/mol/Å² force constant and minimised, then switched to the full OPLS constant and minimised again.
- **Network building:**
  - Initiator: 4.5% of alkene end groups are initiated (H added to make sp3).
  - Cycle: 10 ps NVT at 398 K (20,000 steps at 0.5 fs), then form a bond if the shortest radical-alkene or radical-radical distance is below r_bond, then do the two-stage minimisation.
  - r_bond values: 5 Å for intact networks (78-82% of monomers react); 8 Å for fragments (70-78%, versus 53-59% at 5 Å).
  - Fragments are made in a 100 Å box with a confining potential u_box = P(dx − s)² (s = 20 Å, P = 0.0001 kcal/mol/Å²). Eight fragments are replicated 2×2×2 and melted together at 398 K, 1 atm, 500 ps NPT.
- **Exchange method:** custom MC/MD cycle.
  - **Cycle steps:**
    1. 20 ps NPT at 429 K (Tg + 30 K), 1 atm, 0.5 fs. The 20 ps comes from about 15 ps of volume relaxation after a swap.
    2. Minimise; record E_old.
    3. Find all alcohol-ester pairs within 4.25 Å (the minimum of the R_long distribution). R_long is the longer of the O(alcohol)-C(sp2 ester) and H(alcohol)-O(sp3 ester) distances. Pick one pair with probability Pi ∝ exp(−βR²_i,long), where β = 0.2 Å⁻² (Pi < 0.03 near the cutoff).
    4. Two-stage minimisation; record E_new and compute the reverse-move probability.
    5. Accept with min[1, exp(−(E_new − E_old)/kBT)], T stated as 398 K.
  - **Run:** 8000 cycles for most systems; analysed systems go up to 10,000. Acceptance was above 90%.
  - **Barrier:** the paper explicitly says the method "does not consider the activation energy", so cycles are a proxy for time.
  - **Mechanism:** associative.
- **System size:**
  - Intact network: 343 monomers × 73 atoms = 25,039 atoms.
  - Fractured: 64 monomers (4672 atoms) per fragment; 8 fragments give 512 monomers and 37,376 atoms, with a box of about 75 Å after melting.
  - Five independent systems each.
- **Software:** LAMMPS (version not stated), LigParGen, Python 2 scripts for the reactions.
- **Properties and analysis:**
  - *Tg:* method as in their earlier paper (Wick, Peters, Li, Polymer 281 (2023) 126114); the SI is not in this text, so not stated here. Uncertainty is the SEM over 5 systems.
  - *Shape memory, programming:* compress 40% at 429 K over 4 ns. Anneal 1 ns at 600 K at fixed volume, then 2 ns at 429 K at fixed strain; this gives σ_relax.
  - *Shape memory, recovery:* cool to 298 K over 4 ns at fixed strain, then 4 ns at 1 atm. Reheat to 429 K over 4 ns, then hold 4 ns. σ_rec is the average over the last 0.5 ns. Values are the mean of 15 runs (3 directions × 5 systems).
  - *Tensile test:* 2 ns NpT at 298 K first. Then 50 steps of (0.01%/ps for 50 ps, then a 50 ps hold), reaching 28% strain, followed by a 2.5 ns hold at 28%.
  - **Recycling efficiency:** mean stress over 2.0-2.5 ns at 28% strain for the recycled system, divided by the same quantity for the intact system, × 100%.
  - *Topology:*
    - Distinct monomers bonded per monomer, counted along paths between ether O atoms bonded to aromatic C (EOAC).
    - Largest-network size.
    - <N_links>: minimum path length between monomer pairs. It is compared with a diamond network of the same size, <N_DN_links> = 3.1724 + 0.018629N − 3.2863×10⁻⁵N² + 2.4977×10⁻⁸N³.
    - Pearson correlation coefficients of these with σ_rec and with recycling efficiency.
  - *Not studied:* stress relaxation fits, Tv and MSD.
- **Key numbers:**
  - **Tg:** intact 399 ± 4 K; after 1000 cycles 402 ± 8 K; after 4000 cycles 403 ± 8 K. The paper quotes 348-431 K from experiment.
  - **Stresses:**
    - Intact: σ_relax = 62 ± 5 MPa, σ_rec = 51 ± 3 MPa (experiment: 28 and 13.4 MPa).
    - Recycled 10,000 cycles: σ_rec = 56 ± 4 MPa.
    - σ_rec against cycle count: correlation coefficient 0.999.
  - **Tensile:** intact plateau 170 MPa (experiment 37 MPa); about 120 MPa at the 28% hold.
  - **Recycling efficiency:** 47.5% after 500 cycles, 71% after 1000 (experiment best 69.5%), 100% at 10,000.
  - **Topology:** correlation of <N_sim>/<N_DN> with σ_rec is −0.95 when ≥ 70% of monomers are percolated.
- **CAPS relevance:** CAPS has OPLS-AA typing. Three parts can be copied:
  - Energy-Metropolis transesterification swaps with a distance-biased pair choice and two-stage minimisation. CAPS can do this with its own engine, or as fix bond/react with a probability.
  - The topology descriptors: largest network, and N_links versus the diamond reference.
  - The recycling-efficiency definition.

  First target: Tg of the intact network (399 ± 4 K), then recycling efficiency against cycles. The 25,039-atom system is above the CAPS 15k limit, so a reduced replica would be needed [my inference].

---

### 7. Zheng, Varshney, Vashisth 2025, Macromolecules (Accelerated ReaxFF vitrimers; preprint ChemRxiv)

- **Local file:** Zheng2025_Macromolecules_AccReaxFF_chemrxiv.pdf (ChemRxiv version, posted 4 February 2025).
- **Citations:**
  - Preprint: Yiwen Zheng, Vikas Varshney, Aniruddh Vashisth. "Accelerated ReaxFF Simulations of Vitrimers with Dynamic Covalent Adaptive Networks." ChemRxiv (2025), DOI 10.26434/chemrxiv-2025-5vptv.
  - Published: same authors and title, Macromolecules 58 (2025) 4948-4958, DOI 10.1021/acs.macromol.5c00501.
- **Crossref:** the preprint DOI is verified (posted-content, 2025-02-04). Its record has the relation "is-preprint-of 10.1021/acs.macromol.5c00501". The Macromolecules record (title, three authors, volume 58, pages 4948-4958, 2025) was found by a bibliographic query and matches. **Published:** yes. The local text is the preprint, so changes in the published version were not checked [unverified].
- **Chemistry:** transesterification between β-hydroxy esters (associative). The network is DGEBA cured with adipic acid (epoxy-carboxylic acid, 1:1). The catalyst is not modelled explicitly; the number of simultaneous bond boosts stands in for catalyst concentration. Rubber relevance: none (epoxy vitrimer).
- **Model level and force field:** ReaxFF, force field CHON2017_weak_bb (Vashisth et al.). Charges from EEM. Timestep 0.25 fs.
- **Exchange method:** Accelerated ReaxFF bond boost (Vashisth et al. J. Phys. Chem. A 2018, after Miron and Fichthorn).
  - **Trigger:** when labelled reactive atoms fall inside set distance windows, a one-sided restraint is added for N = 10,000 steps.
  - **Restraint:** E_rest = F1[1 − exp(−F2(R_ij − R12^t))]², with a target that moves linearly, R12^t = R0(1 − t/N) + R12 t/N. R12 is the equilibrium bond length for bond formation and 2.5 Å for dissociation.
  - **Transesterification parameters (Table 2):**

    | Pair | Window (Å) | R12 (Å) | F1 (kcal/mol) | F2 (Å⁻²) |
    |---|---|---|---|---|
    | Ce-Oe | 1.0-1.8 | 2.50 | 50 | 0.25 |
    | Ce-Oh | 5.0-8.0 | 1.40 | 205 | 1.00 |
    | Oe-Hh | 5.0-8.0 | 0.95 | 250 | 0.25 |
    | Oh-Hh | 0.8-1.1 | none | none | none |

    Oe and Oh are given the same type so the reaction can run in reverse.
  - **Curing parameters (Table 1):**

    | Pair | Window (Å) | R12 (Å) | F1 (kcal/mol) | F2 (Å⁻²) |
    |---|---|---|---|---|
    | Ce-Oe | 1.0-1.8 | 2.50 | 70 | 0.25 |
    | Ce-Oa | 3.0-10.0 | 1.40 | 210 | 0.95 |
    | Oe-Ha | 2.5-10.0 | 0.95 | 90 | 0.30 |
    | Oa-Ha | 0.8-1.1 | none | none | none |

  - **Parameter fitting:** F1 and F2 are chosen by Bayesian optimisation, minimising the peak restraint energy:
    - Gaussian-process surrogate with an expected-improvement acquisition function;
    - 20 random initial sets, then 100 iterations on each of 5 configurations;
    - monomer-pair runs, NVT with a Berendsen thermostat (100 fs), at 423 K for curing and 500 K for transesterification.
  - **Temperature dependence:** F1/F1(500 K) = F0 exp(−αEa/RT), with **Ea = 14.6 kcal/mol** and α = 0.2 chosen, giving a reaction-probability transition between 300 and 500 K.
    - The reaction probability is checked as successful BERs out of the possible 1000 (100 boosts × 10 boost periods, 25 ps).
    - Peak restraint energy for the BER: 15.96 kcal/mol, "comparable with" the 14.6 kcal/mol barrier.
  - **Not explicit:** the acceptance rule. A reaction counts as successful if two bonds break and two form.
- **Source of the 14.6 kcal/mol barrier (ref. 55):** Shusil Bhusal, Changjun Oh, Youngjong Kang, Vikas Varshney, Yixin Ren, Dhriti Nepal, Ajit Roy, Gary Kedziora. "Transesterification in Vitrimer Polymers Using Bifunctional Catalysts: Modeled with Solution-Phase Experimental Rates and Theoretical Analysis of Efficiency and Mechanisms." J. Phys. Chem. B 125 (2021) 2411-2424.
  - **DOI:** 10.1021/acs.jpcb.0c10403. **Crossref:** verified (all 8 authors, title, journal, volume 125, pages 2411-2424, 2021).
  - The Zheng text calls it "a previous DFT study". I did not read Bhusal itself, so whether 14.6 kcal/mol appears there exactly, and for which catalyst, is [unverified].
  - The same paper is reference 51 in Hamzehlou and Ruipérez 2022.
- **System size and run length:**
  - 150 adipic acid + 150 DGEBA in a 70 × 70 × 70 Å³ box. My count from the formulas (DGEBA C21H24O4, 49 atoms; adipic acid C6H10O4, 20 atoms) gives about 10,350 atoms [derived, not stated in the paper].
  - Curing: NVT at 423 K, 1 boost at a time, N = 10,000 steps of 0.25 fs, stopped after 15 ps with no new bond. Final cure ratio 92%.
  - Anneal: 423 K down to 100 K in 25 ps, then cycles from 100 K to 800 K in 100 K steps under NPT, 155 ps per cycle, until density converged (by the fifth cycle).
  - Creep: heat from 100 K to 800 K at a constant rate under constant x stress (NPT in y and z, PBC). Reference case: 50 MPa, 4 K/ps, 50 boosts. Other stress, rate and boost values are in the figures, not the text.
- **Software:** Amsterdam Modeling Suite AMS2024 (ReaxFF engine) and Multiwfn for free volume.
- **Properties and analysis:**
  - *Non-isothermal creep:* εx = lx(t)/lx(0) − 1 against T, and dε/dT.
  - *Reaction count:* number of successful BERs against T.
  - *Free volume:* fractional free volume from a voxel grid against vdW radii (Multiwfn), normalised to the 100 K value.
  - *Reactive-site distribution:* normalised mean nearest-neighbour distance of ester carbons, NND = (mean d_i)/(V/n)^(1/3).
  - *Bond integrity:* Ce-O RDF.
  - *Not studied:* stress relaxation fit, Arrhenius fit of relaxation times, MSD and healing efficiency.
  - *Tv:* not measured directly. The 300-500 K onset of BERs is said to agree with experimental Tv for transesterification vitrimers.
- **Key numbers:** cure ratio 92%; peak BER restraint 15.96 kcal/mol versus 14.6 kcal/mol from DFT; α = 0.2; BER onset at 300 K; large pores from about 600 K and rupture at 800 K (50 MPa, 4 K/ps, 50 boosts). Trends match Hubbard et al. 2021 experiments: more creep with higher force, slower heating and more catalyst.
- **CAPS relevance:** ReaxFF bond boost is outside the CAPS force fields (OPLS-AA, PCFF, COMPASS, GAFF). What carries over:
  - The Arrhenius-scaled reaction probability with Ea = 14.6 kcal/mol. It can become the p(T) in a fix bond/react transesterification export.
  - The non-isothermal creep protocol: 100 to 800 K at 4 K/ps under 50 MPa uniaxial, with strain and BER counting.
  - The NND and FFV analyses.

  The system (about 10k atoms) fits the CAPS 15k limit. First target: the DGEBA/adipic acid cure to about 92% with fix bond/react, then the creep onset near 300-500 K with a p(T) built from Ea = 14.6 kcal/mol.

## Group E

---

### 1. Yuan 2026, Materials: all-atom REACTER transesterification and nanoscale healing

**Local file:** Yuan2026_Materials_transesterification_healing.pdf

**Citation:** Yuan, H.; Knight, A.; Jiang, L.; Huang, L. Molecular Insights into Dynamic Relaxation in Vitrimers: Transesterification-Driven Self-Healing at the Atomistic Level. *Materials* 2026, 19 (17), 3702. Published 31 August 2026.

**DOI:** 10.3390/ma19173702

**Crossref:** verified. Title, all four authors, journal (Materials), volume 19, issue 17, article 3702, online 2026-08-31 all match the PDF.

**Chemistry**
- Exchange: transesterification, ester + hydroxyl, associative (explicitly: "associative exchange preserves the total bond count").
- Polymer: Bis-GMA (bisphenol A-glycidyl methacrylate) / 2-HEMA (2-hydroxyethyl methacrylate) methacrylate assembly, 1:2 molar ratio, made by free-radical-like vinyl coupling in REACTER.
- Crosslinker: Bis-GMA itself (two methacrylate groups).
- Catalyst: none modelled; exchange is gated by an imposed temperature-dependent acceptance probability.
- Rubber relevance: none directly (glassy dental/vinyl-ester type resin, Tg 394 K). Transferable as a workflow for ester-hydroxyl exchange in any polyester-type network.

**Model level and force field**
- All-atom, OPLS-AA (Jorgensen 1996 and OPLS/2020 for unsaturated hydrocarbons, alcohols, ethers). Full parameters in SI Table S1 [not seen, SI not in folder].
- LJ + Coulomb cutoff 14 Å, PPPM long-range electrostatics.
- All partial charges scaled by 0.8; one Bis-GMA hydroxyl atom type has a charge modification documented in SI [unverified].
- Monomers built in a commercial builder, packed with Packmol, visualised in OVITO 3.15.0.

**Exchange method and parameters**
- Algorithm: LAMMPS REACTER (fix bond/react) with pre/post templates and atom maps.
- Stage 1, polymerisation: vinyl C-C coupling templates (monomer-monomer and monomer-chain), 4.0 Å C-C cutoff, acceptance probability p = 1.0, NPT 300 K, 1 atm, 5 ns. p is stated to be a construction gate, not a rate constant.
- Stage 2, transesterification: separate template set, ester-hydroxyl, 2.82 Å reactive-pair cutoff. Acceptance probability is a logistic function of T (Eq. 1, reconstructed):
  p(T) = 1 / { exp[ -a (T - T_BER) - ln( p_TBER / (1 - p_TBER) ) ] + 1 }
  with T_BER = 400 K, p_TBER = 0.90, w = 20 K, a = 2 p_TBER / [w (1 - p_TBER)] = 0.90 K^-1.
- Barrier / activation energy: not computed. Authors state REACTER "does not calculate reaction free-energy barriers or absolute chemical rates".
- Attempt frequency/interval for fix bond/react: not stated (Nevery not given). Accepted-event counts reported per 100 ps (Fig. 2).
- Associative. Catalyst not treated. Polymerisation and exchange are run sequentially, not coupled.
- Sensitivity: T_BER shifted by -40 and +40 K only analytically (Fig. S2, Table S3), no extra trajectories.

**System size and run length**
- 300 Bis-GMA + 600 2-HEMA = 900 monomers. Final box about 6.9 nm (Fig. 3 caption). Atom count not stated.
- Polymerisation 5 ns; BER relaxation 5 ns at 400 K for topology figure; healing up to 5 ns; Tg cooling 600 K to 200 K in 20 K steps (time per step not stated).
- One network realisation only; no replicas ("controlled within-realization paired design").

**Software:** LAMMPS 2 Aug 2023 with REACTER; a commercial builder (monomers only); Packmol; OVITO 3.15.0.

**Properties and exactly how they were analysed**
- Density: uncured mixture 1.08 g cm^-3 at 300 K, 1 atm vs ideal-volume estimate 1.13 g cm^-3 (from pure densities 1.16 and 1.07); -4.3 %. After polymerisation 1.17 g cm^-3.
- Conversion/topology: fraction of monomers in at least one coupling (about 99 %); graph analysis of connected components: 247 clusters, lower bound 900 - 247 = 653 inter-component events; only 23 clusters > 300 atoms; motifs chain-like, cyclic, mixed. After 5 ns BER at 400 K: 247 -> 173 components.
- Tg: stepwise cooling 600 -> 200 K, 20 K increments, specific volume after equilibration at each T; alpha = (1/V)(dV/dT)_P; Tg from intersection of glassy and high-T linear fits of V(T). Tg ≈ 394 K.
- Tv: explicitly NOT determined. T_BER is the imposed schedule reference temperature, distinct from Tv.
- Modulus: uniaxial tension along x, engineering strain rate 1 x 10^-5 fs^-1 (1 x 10^10 s^-1), lateral faces at zero pressure; Young's modulus from initial linear region, called "apparent finite-rate tensile modulus". 300, 400, 450, 500 K with and without BER.
- Stress relaxation fit (KWW etc.): not done.
- Creep, MSD/diffusion: not done.
- Healing: fracture by stretching beyond failure, same fractured configuration for all runs, faces brought into contact under "gentle compression" (value not stated), relaxed up to 5 ns, then second tensile test. Healing efficiency eta = E_healed / E_original (apparent small-strain modulus). Also tracked vdW interaction energy, density, cumulative reaction count and rate vs time.
- Error bars: "descriptive variability markers", not standard deviations.

**Key quantitative results**
- Tg ≈ 394 K (within 383-433 K reported for related vinyl-ester thermosets).
- BER reduces apparent modulus by 0.51 GPa at 400 K and 0.11 GPa at 450 K relative to non-reactive control; converge at 500 K; identical at 300 K. Absolute moduli not given in text (figure only).
- Strain 0.5 reached in 50 ps at the applied rate.
- Healing at 300 K: eta ≈ 31.6 % (0.1 ns), 69 % (2 ns), 75.5 % (5 ns). At ≥ 400 K non-reactive: about 98 % within 5 ns. At 400 K, BER gives about 8 percentage points higher average recovery during 0-2 ns; converge later. At 500 K reactive and non-reactive nearly identical.
- Non-reactive vdW/density stabilisation: not reached in 5 ns at 300 K, about 3.8 ns at 400 K, within 0.2 ns at 500 K.

**Relevance for CAPS**
- Closest match in this group to the LAMMPS fix bond/react export path with OPLS-AA. First reproduction target: a two-stage REACTER workflow (vinyl polymerisation templates, then ester-hydroxyl exchange template with Arrhenius/logistic `prob` per temperature), 900-monomer box, Tg by stepwise cooling with bilinear fit, and constant-rate tension with BER on/off.
- Caveat: the logistic p(T) schedule is an imposed input, not chemistry; a CAPS implementation should expose `prob` per temperature and say clearly that it is not a rate.
- System is probably above the 15,000-atom cap: Bis-GMA (C29H36O8, 73 atoms) x 300 + 2-HEMA (C6H10O3, 19 atoms) x 600 ≈ 33,300 atoms [my estimate, atom count not stated in paper]. A CAPS reproduction would need to be scaled down (for example 100 + 200 monomers).

---

### 2. Sun 2026, JACS: polyacrylate-isocyanide vinylogous-urethane-like vitrimer

**Local file:** Sun2026_JACS_polyacrylate_vitrimer_MD.pdf

**Citation:** Sun, H.-L.; Drakopoulos, S. X.; Čamdžić, L.; Maguire, S. M.; Priestley, R. D.; Stache, E. E. Polyacrylate Vitrimer Network via In Situ Isocyanide Copolymerization: Synthesis and Molecular Dynamics. *J. Am. Chem. Soc.* 2026, 148 (10), 11047-11057. Published online 6 March 2026.

**DOI:** 10.1021/jacs.5c22502

**Crossref:** verified. Title, six authors, JACS 148 (10), 11047-11057, online 2026-03-06 (print 2026-03-18) match.

**Important:** this is an EXPERIMENTAL paper. "Molecular Dynamics" in the title (and the local file stem "_MD") means molecular relaxation dynamics measured by broadband dielectric spectroscopy, not MD simulation. No simulation in the main text. The acknowledgements thank someone for "assisting DFT calculations"; any DFT is in the SI [not seen].

**Chemistry**
- Exchange: transamination of enamine (vinylogous-urethane-like, sterically hindered) linkages, associative, activated by free amine.
- Polymer: poly(methyl acrylate) (PMA); also 2-methoxyethyl acrylate (MEA) and 2,2,2-trifluoroethyl acrylate (TFEA).
- Crosslinker: multifunctional isocyanides, mainly 1,4-diisocyanobutane (DICBu), radical copolymerised with acrylate; beta-imino ester tautomerises to enamine.
- Initiator: BAPO (phenylbis(2,4,6-trimethylbenzoyl)phosphine oxide) 0.1-0.6 mol %, white LED; PET-RAFT variant with Ir(ppy)3 0.1 mol %.
- Catalyst/activator: 1,4-butyl diamine (BD), 0.10-0.15 mol %. Trace water hydrolysis releases amine in BD-free samples.
- Rubber relevance: soft acrylate network (Tg about 19 °C by DMA tan delta), stretchy; elastomer-like at room temperature. Cites NR vitrimers for dielectric elastomers. No ENR/SBR/NR/PDMS directly.

**Model level / force field:** none (experimental).

**Exchange method (experimental):** associative amine-enamine exchange; Ea = 108.2 kJ/mol from Arrhenius fit of BDS bond-exchange relaxation times (annealed PMA-3%DICBu-0.15%BD), consistent with rheology stress relaxation times (SI Figs. S48, S49) [stress relaxation fit form in SI, not seen].

**System:** films from 4.0 mmol MA scale; compositions in Table 1.

**Software:** not stated.

**Properties and how analysed**
- Gel fraction in THF, swelling ratio in solvents; Mn and dispersity of hydrolysed networks by GPC vs polystyrene.
- Isocyanide incorporation from IR peak ratio C=O (ester) / C=C (enamine) using PMA-3.3%BI linear copolymer as standard.
- DMA: E' rubbery plateau reached at about 50 °C; tan delta peak 19 °C taken as Tg (agrees with DSC).
- Reprocessing: hot press 2 MPa, 60 °C (thin film); shred and compression mould 80 °C, 5.5 MPa, three cycles, tensile test each cycle (values in figure only).
- BDS: epsilon' vs T at 0.1 Hz; derivative permittivity epsilon''_der to remove dc conductivity; Havriliak-Negami fits to epsilon'' and epsilon''_der (shape parameters Table S5); relaxation maps; segmental alpha relaxation fitted with VFT (fragility D, T0); bond-exchange relaxation fitted with Arrhenius after annealing.
- Tv: not given as a number; footnote says a specific Tv value cannot be established from frequency-dependent BDS. The bond-exchange peak at about 0.2 Hz at 30 °C is described as near the onset/Tv region.
- Conductivity sigma0 (eq S8) vs BNN relation of bond-exchange Delta epsilon and tau (eq S9).
- Dielectric permittivity and loss tan delta at 10 kHz.

**Key quantitative results (validation targets)**
- Table 1 (x mol % DICBu, y mol % BAPO: Mn kDa, D, gel %): 2/0.3: 7.9, 2.3, 90; 3/0.3: 7.4, 2.2, 90; 5/0.3: 4.3, 2.0, 99; 10/0.3: 2.8, 2.3, 90; 15/0.3: 2.0, 1.7, 89; 5/0.1: 8.7, 1.8, 96; 5/0.6: 2.4, 1.4, 94.
- Incorporation 4.0 % for PMA-2%DICBu. Swelling < 180 % except DCM 215 %; gel fraction > 90 % in all solvents.
- IR: C=O 1650 cm^-1, C=C 1595 cm^-1, residual -NC 2148 cm^-1.
- Tg (DMA tan delta) 19 °C for PMA-3%DICBu-0.15%BD.
- epsilon' steps at about 10 °C (segmental) and 35 °C (bond exchange) at 0.1 Hz for PMA-2%DICBu-0.10%BD.
- Ea bond exchange = 108.2 kJ/mol (higher than classic vinylogous urethanes; attributed to steric hindrance).
- Annealing 50 °C overnight shifts segmental relaxation by about 15 K to higher T.
- epsilon' 3.5 to 7.5 (10 kHz, over T); tan delta_BDS < 0.1 at segmental peak, about 0.03 at higher T.
- beta relaxation -50 to -10 °C (carboxymethyl side-group rotation).

**Relevance for CAPS:** low for simulation. Useful only as an experimental Ea (108.2 kJ/mol) and Tg target for a transamination fix bond/react model of PMA networks. Nothing to reproduce first.

---

### 3. Mamat 2026, ACS Appl. Polym. Mater.: MD-MC bead-spring vitrimer, tunable exchange interval, creep and TTS failure

**Local file:** Mamat2026_ACSApplPolymMater_creep.pdf

**Citation:** Mamat, A.; Alrashdan, R.; Eakin, T.; Pandya, H.; Cuddihy, P.; Barz, K.; Cloitre, M.; Khabaz, F. Tuning Bond Exchange Rate to Control Long-Time Mobility and Creep in Vitrimer-like Networks. *ACS Appl. Polym. Mater.* 2026, 8 (16), 13492-13504. Published online 12 August 2026.

**DOI:** 10.1021/acsapm.6c00856

**Crossref:** verified. Title, eight authors, ACS Appl. Polym. Mater. 8 (16), 13492-13504, online 2026-08-12 (print 2026-08-28) match.

**Chemistry:** generic associative vitrimer (no specific chemistry). Tetra-functional crosslinkers + short linear oligomers; motivated by epoxy networks. Rubber relevance: generic coarse-grained elastomer-like network; no specific rubber.

**Model level and force field**
- Coarse-grained bead-spring (Kremer-Grest type), reduced LJ units (sigma, epsilon, m; tau = sqrt(m sigma^2 / epsilon)).
- Nonbonded: truncated and shifted LJ, rc = 2.5 sigma (attractive).
- Bonds: FENE, k = 30 epsilon/sigma^2, R0 = 1.5 sigma. (Text cites Kremer-Grest for FENE; whether the WCA part of standard KG FENE is included is not stated.)
- Network build: melt equilibrated at P = 0, T = 1; each crosslinker bonds to nearest available chain ends (two arms of one crosslinker cannot bond the two ends of the same chain) until all four arms and all chain ends are bonded; then 5 % of bonds removed (no chain loses both ends; free crosslinker sites act as "defects" that initiate swaps). FENE temporarily replaced by harmonic springs until mean bond length r = 0.96, then FENE restored; NPT then NVT at T = 1.

**Exchange method and parameters**
- Algorithm: hybrid MD-MC associative swap. Reaction A-B_i + B_j(free) <-> A-B_j + B_i(free): free crosslinker site j replaces bonded crosslinker site i at chain end k.
- Only free crosslinker sites (defects, 5 %) initiate moves.
- Attempt interval t_Exch = 100, 200, 500, 1000 (MD time units, reduced) and infinity (permanent network).
- Acceptance (Metropolis-Hastings, reconstructed from garbled text): P_acc = min[1, exp(-ΔU_Exch / T) N_f / N_r], N_f = eligible forward candidates in current topology, N_r = eligible reverse candidates in trial topology, same geometric cutoff and constraints.
- Geometric cutoff value for eligible partners: not stated.
- Barrier: no imposed barrier; ΔU from the force field. Apparent activation obtained from tau_bond(T).
- Associative; catalyst not modelled.

**System size and run length**
- Np = 2880 linear oligomers x 5 beads + Nc = 1440 tetrafunctional crosslinkers. Total beads 2880 x 5 + 1440 = 15,840 if a crosslinker is one bead [my arithmetic; the text mentions "central beads of the cross-linkers" and "cross-linker reacting bead", so crosslinkers may be multi-bead and the true total is not stated].
- Time step 0.005 tau. Five independent replicas per t_Exch.
- Quench from T = 1.0 in 16 stages of -0.05, 10^7 steps each, cooling rate 10^-6 (reduced).
- Lag times to at least 10^5 tau (van Hove at Δt = 10^4 and 10^5). Total production length not stated.

**Software:** not explicitly stated in text. LAMMPS papers (Plimpton 1995; Thompson 2022) are cited in the introduction only [LAMMPS use unverified]. TACC resources acknowledged.

**Properties and how analysed**
- Tg: V per bead vs T from the quench; Tg from change in slope of V(T), also from alpha = (1/V)(dV/dT)_P mapped back. Tg ≈ 0.43 (all t_Exch). Volumetric Tg, explicitly distinct from topology freezing.
- Relative volume V_R = (V - V_p)/V_p vs permanent network.
- g(r) of central crosslinker beads at T = 0.8 and 0.3.
- Exchange statistics: N_total accepted swaps and N_average accepted per MC application vs T.
- MSD of crosslinker central beads (NVT) at T = 0.40, 0.60, 0.80, 1.00; plateau, subdiffusive ~Δt^0.5, then ~Δt for fast exchange.
- Self van Hove G_s(r, Δt) at T = 0.80.
- Bond memory: h_b(t) occupancy, C(t) = <h_b(0) h_b(t)>_b / <h_b(0)>_b; tau_bond from C(tau_bond) = 1/e; Arrhenius fit ln tau_bond = ln tau_0 + E_app/T.
- Creep: constant shear stress sigma_0 = 0.1, J(t) = gamma(t)/sigma_0, at T = 0.40, 0.60, 0.80, 1.00.
- TTS: Tref = 0.6; horizontal a_T by manual shifting; creep also vertical b_T = (rho_ref T_ref)/(rho T) [reconstructed]; WLF fit to permanent-network a_T only.
- Stress relaxation G(t), KWW: not done. Tv: not determined here (earlier fast-exchange work did).

**Key quantitative results**
- Tg ≈ 0.43 (reduced).
- N_total spans more than 2 orders of magnitude from glassy to rubbery; ordered t_Exch 100 > 200 > 500 > 1000; N_average ~ linear in T and independent of t_Exch.
- E_app = 1.225 ± 0.090 (t_Exch 100), 1.134 ± 0.075 (200), 0.908 ± 0.031 (500), 0.995 ± 0.052 (1000), reduced units, R^2 = 0.979-0.996. tau_0 = 0.086 t_Exch.
- tau_bond = m(T) t_Exch + b(T): from T = 1.0 to 0.5, m 0.230 -> 0.599, b 0.164 -> 12.838.
- WLF (permanent network): MSD C1 = 2.001, C2 = 0.282; creep C1 = 1.560, C2 = 0.262.
- TTS works for permanent network, fails for all finite t_Exch (no single a_T).

**Relevance for CAPS**
- Best first target in this group for a Kremer-Grest module with an MC swap move in the CAPS MD engine: small parameter set (LJ 2.5 sigma, FENE 30/1.5, dt 0.005, 5 % defects, t_Exch scan, Metropolis-Hastings with N_f/N_r correction). Outputs to reproduce: tau_bond(T) Arrhenius with E_app about 0.9-1.2, Tg ≈ 0.43 from V(T), MSD crossover, J(t) at sigma_0 = 0.1.
- Watch the 15,000-atom cap: about 15,840 beads as published; halve the system (1440 chains, 720 crosslinkers) for local checks. The 10^7-step-per-stage quench is far beyond the 5,000-step local-check limit, so full runs belong on the cluster.
- Missing details to fix before reproducing: swap cutoff distance, thermostat/barostat type, crosslinker bead geometry.

---

### 4. Danikas 2026, Macromolecules: dioxaborolane PBMA vitrimers under temperature and pressure (experimental)

**Local file:** Danikas2026_Macromolecules_dioxaborolane_elastomer.pdf

**Citation:** Danikas, M.; Conti, M.; Mattsson, J.; Van Zee, N. J.; Nicolaÿ, R.; Ramos, L.; Floudas, G. Temperature and Pressure Effects on the Dynamics of Elastomeric Dioxaborolane Vitrimers. *Macromolecules* 2026, 59 (18), 10544-10556. Published online 5 September 2026.

**DOI:** 10.1021/acs.macromol.6c01360

**Crossref:** verified. Title, seven authors, Macromolecules 59 (18), 10544-10556, online 2026-09-05 (print 2026-09-22) match.

**Experimental only.** No simulation in the paper.

**Chemistry**
- Exchange: dioxaborolane metathesis (boronic ester exchange), catalyst-free, associative.
- Polymer: poly(n-butyl methacrylate) (PBMA) via RAFT; functional comonomer 5,6-dihydroxyhexyl methacrylate (DM, 9.2-12.0 mol %), condensed with phenylboronic acid to pendant dioxaborolanes (9.2-10.4 mol %).
- Crosslinker: 1,4-bis(4-methyl-1,3,2-dioxaborolan-2-yl)benzene (CL); mono-dioxaborolane by-product removed under vacuum at 120 °C.
- Samples: homopolymers H85, H290; precursors P85, P306, P348; vitrimers V85-3%, V306-1%, V348-5%, V348-9% (crosslink mol %).
- Rubber relevance: "elastomeric" PBMA vitrimers (Tg about 293-303 K, rubbery plateau above). Not ENR/SBR/NR/PDMS but a direct elastomer vitrimer benchmark.

**Methods and analysis**
- TM-DSC (Q2000): ±1.0 K amplitude, period/rate 20 s/10.0, 40 s/5.0, 60 s/3.3, 100 s/2.0, 150 s/1.3 K min^-1, 253-373 K. Tg = peak of d(c_p,reversing)/dT.
- TGA: N2, 303-773 K, 10 K min^-1.
- Rheology: AR-G2, 8 mm parallel plates, frequency sweeps 0.1-100 rad s^-1; master curves at Tref = 353 K by horizontal shift on tan delta then vertical shift; WLF below Tv, Arrhenius a_T = a0 exp(E*/RT) above.
- Tv definition: crossover from WLF to Arrhenius in a_T, located by chi^2 of WLF fit (chi^2 < 1.5 x 10^-3 good, ~3 x 10^-3 fails) and more sensitively by plotting [-d log a_T / d(1/T)]^(-1/2) vs 1/T (Arrhenius is a horizontal line). Explicitly different from Leibler's 10^12 Pa s viscosity definition.
- Dielectric spectroscopy: Novocontrol, 10^-2 to 3 x 10^6 Hz, 393 to 203 K in 5 K steps; pressure cell to at least 240 MPa. HN fits, tau_max from HN parameters; alpha and alpha-beta processes VFT (log f_inf fixed at 13.2), beta and gamma Arrhenius. Tg defined at tau = 100 s.
- Pressure: pressure-VFT with D_P fixed at 11.6; activation volume ΔV# = 2.303 RT (d log tau / dP)_T; Tg(P) = Tg(0)(1 + nu P / mu)^(1/nu); Tait EOS V(T,P) = V(T,0)[1 - 0.0894 ln(1 + P/B(T))]; Q_V/Q_P from density-modified VFT.
- van Gurp-Palmen plots for delta_min and G_N^0.

**Validation-target numbers**
- Precursor Mn (SEC): P85 22, P306 53, P348 63 kg/mol, Đ 1.2. Homopolymer Mn 16 and 48 kg/mol.
- TGA: V306-1% stable below 533 K; V348-5% 2.3 % loss below 533 K; V85-3% 3.8 % loss below 500 K.
- Tg increase for V348-9% about 20 K (TM-DSC). Rheology and DS Tg about 10 K below TM-DSC (40 s period) and within ±2 K of each other.
- Tg(0) at tau = 100 s: H290 295 K, P306 293 K, V348-5% 303 K. dTg/dP: 238 ± 7, 247 ± 5, 253 ± 3 K/GPa. mu: 1235 ± 85, 1179 ± 64, 1193 ± 69 MPa; nu: 7 ± 1, 5 ± 1, 7 ± 1.
- VFT (B K, T0 K), alpha-process: H290 2990 ± 60, 214 ± 2; P306 3070 ± 50, 209 ± 2; V306-1% 3100 ± 50, 212 ± 2; V348-5% 3050 ± 60, 220 ± 2; V348-9% 2763 ± 70, 235 ± 2; H85 3890 ± 140, 183 ± 5; P85 3290 ± 60, 200 ± 2; V85-3% 3100 ± 50, 212 ± 12. alpha-beta: H290 3790 ± 40, 177 ± 2; P306 3850 ± 30, 173 ± 2; V306-1% 3730 ± 30, 181 ± 1; V348-5% 3700 ± 30, 186 ± 1; V348-9% 4002 ± 60, 174 ± 3; H85 3840 ± 20, 172 ± 1; P85 3850 ± 15, 170 ± 1; V85-3% 3720 ± 40, 181 ± 2.
- gamma-process Ea 28 ± 2 kJ/mol (V306-1%).
- ΔV# alpha approaches repeat-unit molar volume about 140 cm^3 mol^-1 about 15 K above Tg; beta ΔV# about 25 cm^3/mol (about 18 %); vitrimer ΔV# about 20 % above precursor at 353 K.
- Q_V/Q_P about 0.75-0.85 (V348-5%).
- vGP (delta_min °, G_N^0 x 10^5 Pa): H290 32, 1.8; P306 34, 1.6; V306-1% 2.27, 1.1; V348-5% 0.12, 3.2; V348-9% 0.7, 8.3; H85 53, 3.9; P85 52, 3.3; V85-3% 0.46, 1.1. Crosslinks per chain CN_f: 3.1, 17.4, 31.3, 2.6; M_xx,theo 11.7, 3.0, 1.7, 3.6 kg/mol; N_xx 75, 19, 11, 23 (V306-1%, V348-5%, V348-9%, V85-3%). PBMA N_e about 125.
- Bond-exchange (flow) Ea from Arrhenius shift factors: 93 ± 5 (V306-1%), 94 ± 2 (V348-5%), 80 ± 2 (V348-9%) kJ/mol. Small-molecule dioxaborolane metathesis about 15.9 kJ/mol (cited). Not reported for V85-3% because of thermal stability.
- Tv/Tg about 1.4-1.5 (individual Tv values only in Fig. 14). Order Tg < T_merge(alpha-beta) < Tv.
- 8 h at 473 K, V348-5%: G'' at 10 rad/s rose by factor two, G' changed 4 %.

**Relevance for CAPS:** experimental validation set for an elastomeric methacrylate vitrimer: Tg (293-303 K), fragility, G_N^0 vs crosslink density (1.1 to 8.3 x 10^5 Pa), Ea ~ 80-94 kJ/mol, Tv/Tg ~ 1.4-1.5. A Kremer-Grest or OPLS-AA PBMA model could compare Tv/Tg ratio and plateau-modulus trend; absolute timescales are out of MD reach. Also a methods reference for locating Tv from a_T(T) in simulated TTS.

---

### 5. Rovigatti 2018, Macromolecules: three-body swap potential, binary patchy-particle vitrimer, q-independent mode

**Local file:** Rovigatti2018_Macromolecules_swap_particle_model.pdf

**Note:** the local file is a preprint / author manuscript (no journal header, title "Self and collective swap-driven dynamics in a particle model for vitrimers"). The published version's title differs slightly.

**Citation (published version):** Rovigatti, L.; Nava, G.; Bellini, T.; Sciortino, F. Self-Dynamics and Collective Swap-Driven Dynamics in a Particle Model for Vitrimers. *Macromolecules* 2018, 51 (3), 1232-1241. Published online 12 January 2018.

**DOI:** 10.1021/acs.macromol.7b02186

**Crossref:** verified for the published version (title as above, four authors match, Macromolecules 51 (3), 1232-1241, online 2018-01-12, print 2018-02-13). Title in local file differs from published title. Numbers below come from the preprint; equality with the published version is [unverified]. Preprint DOI (arXiv) not checked [unverified].

**Chemistry:** generic. Binary mixture: tetrafunctional A (f_A = 4, tetrahedral sites) and bifunctional B (f_B = 2, polar sites); only A-B bonds; near-stoichiometric with B excess. No real chemistry, no catalyst. Associative swap (bond number conserved). Rubber relevance: none directly.

**Model level and force field**
- Patchy particles (continuous version of Smallenburg-Leibler-Sciortino 2013 model).
- Core-core: WCA, 4 epsilon[(sigma/r)^12 - (sigma/r)^6] + epsilon for r ≤ 2^(1/6) sigma.
- Site-site A-B attraction, Stillinger-Weber-like (Eq. 2, reconstructed): V_ss(r) = 2 epsilon [ sigma_ss^4 / (2 r^4) - 1 ] exp[ sigma_ss / (r - r_c) + 2 ] for r ≤ r_c, else 0; sigma_ss = 0.4 sigma, r_c = 1.5 sigma_ss; well depth epsilon.
- Three-body swap term (Sciortino 2017): V_3b(r_ij, r_ik) = lambda epsilon V_3(r_ij) V_3(r_ik), V_3(r) = 1 for r ≤ sigma_ss, = -V_ss(r)/epsilon for sigma_ss < r < r_c. Swap activation energy = (lambda - 1) epsilon; lambda ≥ 1 required.
- Production uses lambda = 1 (zero swap barrier, fastest swapping).
- No FENE, no permanent bonds; bonds are site-site attractions that never break at k_B T = 0.03.

**Exchange method and parameters**
- Three-body potential in plain MD (no MC). Swap A_i B_k + B_j -> A_i B_j + B_k happens whenever a free B site approaches a bonded pair.
- Barrier (lambda - 1) epsilon = 0. No attempt frequency (continuous dynamics). Temperature k_B T = 0.03 epsilon (stated irrelevant as long as << epsilon). Since x < x_fb, all A sites always bonded; only free B sites swap.

**System size and run length**
- N = 800,000 particles; volume 128.16^3 sigma^3 (text "128.163"), rho sigma^3 = 0.38.
- x = N_A/N from 0.300 to 1/3 (Table 1: x, unreacted B sites, p_B, B monomers): 0.300, 160000, 0.857, 11450; 0.320, 64000, 0.941, 1893; 0.325, 40000, 0.963, 739; 0.330, 16000, 0.985, 121; 0.332, 6400, 0.994, 19; 0.3325, 4000, 0.996, 8.5; 0.333, 1600, 0.998, 2; 1/3, 0, 1.000, 0. Number of unreacted sites = 6N(x_fb - x).
- Velocity Verlet, dt = 0.003 t0 (t0 = sqrt(m sigma^2 / k_B T)); modified Andersen thermostat for equilibration, NVE production. Correlation functions to about 10^6 t0 (figure axes); total run length not stated.
- Percolation threshold x_perc = 1/7; all states far above.

**Software:** home-made GPU MD code (name not given).

**Properties and how analysed**
- Structure factors S_AA, S_AB, S_BB, S_NN.
- Bond autocorrelation C_b(t) = [n_b(t) - n_b(inf)] / [n_b(0) - n_b(inf)], n_b(inf) = f_A/N_B; tau_b from C_b(tau_b) = 1/e; tau_b ~ 1/(x_fb - x).
- MSD of A particles; plateau about 0.8 sigma; D from long-time MSD; D = D1(x_fb - x) + D2(x_fb - x)^2.
- F_self(q, t) fit A e^(-t/tau_f) + (1 - A) e^(-(t/tau_s)^beta), beta > 0.85 approaching 1 at low q; tau_self ~ q^-2; Gaussian approximation checked.
- F_coll(q, t), 0.04 ≤ q sigma ≤ 2, fit A e^(-t/tau_f) cos(Omega t) + (1 - A) e^(-(t/tau_s)^beta).
- No stress relaxation, modulus, Tg, Tv, healing.

**Key quantitative results**
- Slow collective time tau_s independent of q (q^0 mode); fast tau_f ~ q^-2 at small q.
- Sound speed 0.53 ± 0.03 sigma/t0.
- tau_s vs (x_fb - x) power law exponent about -1.2 (fig label -1.225); tau_b ~ (x_fb - x)^-1; tau_b at least twice tau_s.
- tau_s D_B collapses curves for all but lowest x.
- C_b(t) long-time apparent power law exponent about 1.5-2.

**Relevance for CAPS:** the reference implementation of the Sciortino three-body swap potential (pair_style sw/mod-like plus three-body term, available in LAMMPS as the swap potential used by later papers [LAMMPS availability unverified here]). For CAPS the reusable part is the potential definition (V_ss, V_3b, lambda) as an alternative to MC swap; the 800,000-particle system is far above the 15,000-particle cap, so only a small-N bond-lifetime check (tau_b vs x_fb - x, MSD plateau 0.8 sigma) is feasible. Not a first target compared with Mamat.

---

### Inconsistencies and caveats found

1. Sun 2026: title and file stem say "Molecular Dynamics" but the paper has no MD simulation; it is dielectric/rheological relaxation dynamics. DFT acknowledged but not in main text.
2. Rovigatti 2018: local file is a preprint whose title ("Self and collective swap-driven dynamics...") differs from the published Crossref title ("Self-Dynamics and Collective Swap-Driven Dynamics..."). Its Appendix cites "Ref. 35" for the swap method, but Ref. 35 is a supercooled-water paper; the method is Ref. 22 (Sciortino 2017). Box size printed as "128.163" means 128.16^3. Fig. 3 caption calls the function n_b(t) while the axis and text use C_b(t).
3. Danikas 2026: abstract says E ∼ 90 kJ/mol; body gives 93 ± 5, 94 ± 2, 80 ± 2 kJ/mol. Small-molecule metathesis Ea given as ∼15.9 kJ/mol in text and 16 kJ/mol in conclusions. Conclusions attribute the ∼90 kJ/mol to "zero-shear viscosity", whereas the body derives it from Arrhenius shift factors. Table 4 V85-3% M_xx is printed garbled ("3.63.6") in the text layer; 13/(2.6 + 1) = 3.6 kg/mol confirms 3.6.
4. Mamat 2026: software never named; total bead count not stated (15,840 only if each crosslinker is one bead); swap geometric cutoff not given; Crossref title contains a line break but otherwise matches.
5. Yuan 2026: atom count not stated (my estimate about 33,000, above the CAPS 15,000-atom cap); fix bond/react Nevery not stated; error bars are explicitly not standard deviations (single network realisation).

## Group F

Crossref checked 2026-10-22 with `https://api.crossref.org/works/<doi>`. arXiv DOIs are registered with DataCite, not Crossref, so the arXiv DOI was checked at `https://api.datacite.org/dois/<doi>`. Published versions of the preprints were searched with a Crossref `query.bibliographic` title query. "[derived]" marks a number I calculated from values in the paper; it is not printed in the paper.

---

### F1. Garcia Villatoro et al. 2026, Evolution of primary loops in vitrimers

**Local file:** `GarciaVillatoro2026_chemrxiv_primary_loops.pdf`

**Citation:** E. A. Garcia Villatoro, S. W. Fitz, J. Patro, A. Martell Monterroza, J. A. Johnson, M. Olvera de la Cruz, J. A. Kalow. "Evolution of primary loops in vitrimers." ChemRxiv, posted 9 September 2026 (preprint, not peer reviewed). CC-BY-NC 4.0.

**DOI:** 10.26434/chemrxiv.15008540/v1

**Crossref check:** found. Type posted-content, title matches, 7 authors in the same order, date 2026-09-09. Crossref lists the sixth author only as "Olvera de la Cruz" with no given name; the PDF has "Mónica Olvera de la Cruz". No journal, volume or pages (preprint). A Crossref title search found no published journal version as of 2026-10-22.

**Chemistry**
- Exchange: transalkylation of trialkylsulfonium salts (Du Prez poly(thioether) vitrimer chemistry). Strictly associative when the alkylating agent has a non-nucleophilic counterion.
- Network: trifunctional thiol A3 (A3H, plus a partially deuterated A3D as the mass label, mixed 1:1) and triethylene glycol divinyl ether (B2), crosslinked by photoinitiated thiol-ene. Photoinitiator DMPA (2,2-dimethoxy-2-phenylacetophenone), 365 nm, as low as 2 mW/cm². Vinyl excess 1.8 equiv (20 % excess per thiol).
- Activator: decyl brosylate or tripropylsulfonium brosylate, 0.15 equiv (5 mol % relative to thiol). It alkylates thioethers above 100 °C; exchange was activated at 150 °C and 170 °C.
- Solvents: diglyme (deswelling series), monoglyme (gel fraction), polyglyme Mn = 250 g/mol as a nonvolatile plasticizer (constant-concentration series). Concentrations 0.05 M to 1 M (1 M = dry).
- Disassembly for NDS: 1 M NaOH in 3:1 MeOH-H2O, 40 °C, 1 h, then LC-MS.
- Rubber relevance: no ENR, SBR, NR, PDMS or disulfide rubber. It is a thioether-sulfonium elastomeric gel/thermoset. The discussion cites transalkylation recycling of sulfur-cured rubber (Tang et al. 2018), the same C-S transalkylation family, and argues that "vitrimerized" elastomers will come out stiffer than the original because they have fewer loops.

**Model level and force field (simulation part)**
- Reactive coarse-grained MD based on Kremer-Grest. Beads interact through WCA with ε = σ = 1.0, cutoff 2^(1/6) σ. Bonds are FENE: U = -0.5 k R0² ln[1 - (r/R0)²], r ≤ R0. The main text does not give the values of k and R0 (not stated; possibly in the SI, which is not in the local file).
- No bending or angle potential ("no bending angles, relatively inextensible bonds, and hard-particle only interactions").
- Components: A3 stars (non-active), A3+ stars (pre-doped, one end with functionality three = sulfonium) and B2 linkers. A primary loop is a 10-bead ring (B2 bonded to two arms of the same A3). 6-bead rings are a B2 bonded twice to the same sulfonium end.
- Box: the text reads "V=200 3σ 3" (garbled in extraction, probably V = 200³ σ³) [unverified]. Number densities quoted: ρ = 0.0815 for the denser intermediate system (C3). Bead counts not stated.
- All-atom FF: none.

**Exchange method and parameters**
- Network formation: permanent diffusion-limited bonding between unbonded A3/A3+ and B2 end beads whenever they come within 2^(1/6) σ, run until >95 % conversion. The loop fraction of this network is the "before exchange" point.
- Exchange: LAMMPS `fix bond/react` with two reversible reaction templates. (1) A doped A3+ sulfonium end that carries a dangling end exchanges with any sulfide end of functionality ≤ 2. (2) Any sulfonium end without a dangling end exchanges with any sulfide end of functionality ≤ 2. Both include the reverse reaction. Associative.
- Initiation cutoff r_react ≤ 0.97 (σ), chosen to keep the FENE bonds stable while exchange runs.
- Reaction probability, barrier, Nevery and temperature of the simulations: not stated in the main text.
- The "+" on sulfonium beads carries no charge in the simulation.

**System size and run length:** three low-concentration systems (one in the sol regime, two near the sol-gel transition). Loops were followed to 3.5 × 10^8 timesteps without the 10-bead loop count saturating. Bead numbers and timestep: not stated.

**Software:** LAMMPS (`fix bond/react`); OVITO ring-finding was used to cross-check the authors' own loop-counting algorithm (SI Fig. S21).

**Properties computed and how**
- Primary loop fraction n1,3 (experiment): NDS. LC-MS ratio of H/H, H/D and D/D disassembly fragments; intermolecular links give all three, primary loops only H/H or D/D. Only thioether fragments were quantified; sulfonium fragments were too weak.
- Loop count vs time (simulation): number of 10-bead and 6-bead cycles from the authors' algorithm, checked with OVITO.
- Gel fraction (simulation): n_GCC / n_tot, beads in the largest connected component over all polymer beads (Eq. 3). The second-largest component was also tracked to confirm dissolution rather than splitting.
- Stress relaxation (experiment): dry samples at 130 to 180 °C fitted to a KWW stretched exponential, giving the average relaxation time ⟨τ⟩. β rose from 0.5 at 130 °C to 1.0 at 180 °C. Arrhenius fit of ⟨τ⟩ gives the flow activation energy. No Tv reported.
- Modulus: oscillatory shear time sweeps at 150 °C (storage modulus G') before and after activation. Two-tailed Student's t-test, p ≤ 0.05.
- Gel fraction and swelling ratio (experiment): soaking in monoglyme at 23 °C for 24 h.
- Tg, fragility, creep, MSD, fracture, healing efficiency: not computed.

**Key quantitative results**
- Ea (flow) = 123 kJ/mol from the Arrhenius fit of the KWW ⟨τ⟩ (Du Prez's uncatalysed transalkylation for comparison: 113 kJ/mol).
- Static networks (diglyme): primary loop fraction 0.02 at 1 M, 0.05 at 0.5 M, 0.32 at 0.05 M (still a sol).
- Deswell then activate at 150 °C for 3 h: all samples lose loops. Samples made at >0.5 M converge to 0.01. A 0.1 M sample went from 0.17 to 0.07 over 5 h.
- Constant concentration (polyglyme): above ~0.2 M loops decrease on activation; below 0.2 M loops increase and the gel turns into a sol. 0.5 M equilibrates in about 3 h at 150 °C. At 170 °C (50 min) both trends are stronger.
- G' at 0.5 M (150 °C): 41 → 71 kPa after equilibration, plateau within 30 min. Swelling ratio 1.51 → 1.61.
- Gel fraction: dry 0.90 → 0.96; 0.5 M deswelled 0.82 → 0.92.
- Recycling: 1 M and 0.5 M networks dissolve in 24 h at 150 °C; 0.3 M in 4 h. A gel re-forms after solvent removal and re-equilibration.
- Simulation: 10-bead loops first drop (converted to 6-bead loops, which saturate about an order of magnitude below), then rise steadily. At ρ = 0.0815 the largest cluster (~45 % of material after permanent crosslinking) shrinks once exchange starts (de-gelation). Conclusion: entropy alone is enough to explain loop growth at low concentration.

**Relevance for CAPS**
- Reproduce first: the simulation is a direct test case for the CAPS Kremer-Grest builder plus `fix bond/react` export. Build A3 stars and B2 linkers (WCA + FENE, no angles), crosslink by distance at 2^(1/6) σ to >95 %, then export the two reversible exchange templates with r_react = 0.97. Gaps that must be filled by choice (not stated): FENE k and R0 (standard KG 30 ε/σ² and 1.5 σ is the obvious guess), temperature, timestep, reaction probability, bead counts.
- Analysis CAPS should add: ring census by ring length (10-bead and 6-bead cycles) over a trajectory, largest and second-largest connected component fractions vs time. Both are pure topology tools and work for any crosslinked network, including vulcanized rubber models.
- Experimental benchmark for an exchange-activated KWW workflow: β 0.5 to 1.0 and Ea 123 kJ/mol.

---

### F2. Guo and Zhang 2026, Force-dependent fracture of dynamic polymer networks

**Local file:** `GuoZhang2026_chemrxiv_dynamic_fracture.pdf`

**Citation:** W. Guo, T. Zhang. "Force-Dependent Fracture of Dynamic Polymer Networks: Sacrificial Bond Rupture and Dispersion-Induced Physical Crosslinking." ChemRxiv, posted 3 September 2026 (preprint, not peer reviewed). CC-BY 4.0. School of Chemistry and Chemical Engineering, Shanghai Jiao Tong University.

**DOI:** 10.26434/chemrxiv.15008249/v1

**Crossref check:** found. Type posted-content, title matches exactly, authors Wenjing Guo and Tao Zhang, date 2026-09-03. No journal, volume or pages. A Crossref title search found no published version as of 2026-10-22.

**Chemistry**
- Generic. The weak linker represents S-S or Se-Se (disulfide or diselenide) bonds, which have lower dissociation energies than C-C. The physical-attraction part is motivated by Se···Se dispersion contacts and diselenide polyurethanes (PU-Se-0/50/100, the authors' own experimental series).
- Exchange in the model is dissociative-type: break and re-form of the weak linker bond only. There is no associative swap.
- No catalyst.
- Rubber relevance: generic elastomer network. Directly relevant to disulfide- or diselenide-crosslinked elastomers and polyurethanes; no ENR, SBR, NR or PDMS.

**Model level and force field (Table I, LJ units)**
- Kremer-Grest bead-spring, flexible, no angle potential. m = 1.0, ε = 1.0, σ = 1.0, τ = σ (m/ε)^(1/2).
- Nonbonded: LJ with r_cutoff = 2^(1/6) σ (purely repulsive) for all pairs except weak-linker/weak-linker pairs. Those use ε_linker = 1.0, 5.0, 10.0, 15.0, 20.0 and σ_linker = 1.0, 1.1, 1.2, 1.3, 1.4, with an extended cutoff of 4.0.
- Bonds: FENE, R0 = 1.5, K_backbone = 30, K_linker = 12, 18, 24, 30. The text writes the FENE term as -0.5 K R0² ln[1 - (r/R0)²]; no WCA-on-bonds term is mentioned.
- kBT/ε = 0.78; initial number density ρ = 0.513 m/σ³.
- No three-body swap potential, no all-atom force field.

**Exchange / rupture method and parameters**
- Bell model, applied to every bond (backbone and linker): k_r(F) = k_r(0) exp(γ0 F), with k_r(0) = ν0 exp(-Ea / kBT).
- Rupture probability per check: P_r = 1 - exp[-k_r(F) Δt_rxn] (≈ k_r Δt_rxn when small).
- γ0 = 1.025; ν0 = 4.0 τ⁻¹; Δt_rxn = 50 Δt = 0.25 τ.
- Ea,backbone = 104.55; Ea,linker = 73.05, 78.05, 88.05, 98.05, 104.55 (LJ energy units). SI Fig. S1 also uses 83.05 and Ea,backbone = ∞ (backbone rupture switched off).
- Regeneration: only the weak dynamic linker, with fixed probability P_f = 0.001 per reaction check. The distance criterion and partner rule for re-forming a bond are not stated.
- [derived] With these numbers a backbone bond needs γ0 F ≈ Ea/kBT - ln(100) ≈ 134.0 - 4.6, i.e. F ≈ 126 ε/σ, to reach P_r ≈ 0.01 per check. A FENE (K = 30, R0 = 1.5) bond carries that force at r ≈ 1.26 σ. The paper does not quote a threshold force.
- How the Bell rupture and regeneration were coded inside LAMMPS is not stated (no fix named). Network crosslinking used `fix bond/react` with a 1.5 σ distance criterion to ~99 % crosslinking.
- Temperature fixed at T = 0.78; no temperature series.

**System size and run length**
- 160 long chains of 10 beads, 160 short chains of 2 beads, 160 crosslinkers of 2 beads (one central linker bond each). N_backbone = 1920, N_linker = 320, total 2240 beads.
- Timestep 0.005 τ. Relaxation with `fix nve/limit` 2.5 × 10³ τ, NPT 2.4 × 10⁴ τ, NVT 2.4 × 10⁴ τ (Nosé-Hoover). Minimisation: steepest descent then conjugate gradient to relative tolerance 1 × 10⁻⁵.
- Uniaxial tension along x at constant volume, engineering strain rate 1.0 × 10⁻³ τ⁻¹ for 1.0 × 10⁴ τ to ε_eng = 10.0. Rate checks at 3 × 10⁻⁴ τ⁻¹ (≈3.33 × 10⁴ τ) and 3 × 10⁻³ τ⁻¹ (≈3.33 × 10³ τ).
- Five independent runs averaged per curve.

**Software:** LAMMPS; initial chains grown with an in-house Python script. Code and data "available from the corresponding author upon reasonable request".

**Properties computed and how**
- Stress-strain: σ_t = -P_xx (virial) vs engineering strain. Strain at break, tensile strength and small-strain modulus are read off the curves. No fracture energy integral is reported.
- Force-dependent rupture: cumulative ruptured linker bonds N_r^linker, ruptured backbone bonds N_r^backbone and regenerated linker bonds N_f^linker vs strain. F_max of any single linker and backbone bond vs strain (SI).
- Aggregation: weak-linker/weak-linker RDF g_ww(r) and cumulative coordination number N_ww(r) = 4π ρ_w ∫ g_ww(r') r'² dr'.
- Not computed: stress relaxation, Ea/Arrhenius, Tv, Tg, fragility, creep, MSD, healing efficiency, loops.

**Key quantitative results**
- Strong reference (Ea,linker = Ea,backbone = 104.55, K = 30) fails at strain ≈ 5.5. Ea,linker = 88.05: break strain ≈ 7.5. Ea,linker = 78.05 and 73.05: no complete failure by strain 10. The initial modulus barely changes.
- Lower K_linker (Ea,linker = 88.05) increases break strain and linker rupture and lowers backbone rupture.
- Regeneration (Ea,linker = 88.05, K_linker = 18, P_f = 0.001): more extensible than the reference but fails earlier than the non-regenerating network. Formation count stays below rupture count. The same ordering holds at all three strain rates.
- Dispersion (σ_linker = 1.2): first RDF peak 18.84 at ε_linker = 5.0 to 79.66 at 20.0. ε_linker = 5.0 raises tensile strength from ≈2.8 to 5.8 and lowers break strain from ≈9.2 to 7.7, with linker rupture down to about one tenth. ε_linker = 10.0: small-strain modulus about twice that at 5.0, break strain ≈4.0. ε_linker = 15 to 20: break strain ≈2.0.
- Dynamic crosslinker content: DCC-0 % fails near 5.5, DCC-100 % is strongest and stiffest but fails near 4, DCC-50 % has a peak strength close to DCC-100 % and extends beyond strain 8.

**Relevance for CAPS**
- Fully parameterised and small (2240 beads, ~5 × 10⁴ τ per run), so it fits the 15k-atom cap and short local checks. Reproduce first: the strong reference and one weak-linker case (Ea,linker = 88.05) under constant-volume tension, comparing break strain ≈5.5 vs ≈7.5.
- The Bell rupture has no named LAMMPS fix in the paper, so the cleanest route is to implement force-dependent rupture (k_r = ν0 exp(-Ea/kBT + γ0 F), checked every 50 steps) in the CAPS MD engine, and to export to LAMMPS only the parts that are standard (`fix bond/react` crosslinking, `fix deform` tension).
- Useful generic analyses to add: cumulative broken-bond counts per bond type vs strain, per-type maximum bond force, and selective pair RDF/CNF.
- Rubber link: a coarse model of disulfide-crosslinked elastomers (sacrificial S-S vs C-C backbone).

---

### F3. Ricarte 2026, Network Topology and Linear Viscoelasticity of Vitrimers (Perspective)

**Local file:** `Ricarte2026_Macromolecules_topology_LVE_chemrxiv.pdf` (the ChemRxiv preprint text)

**Citation (published):** R. G. Ricarte. "Network Topology and Linear Viscoelasticity of Vitrimers - A Polymer Physics Perspective." Macromolecules 2026, 59 (12), 6593-6619. Published online 2 June 2026 (Crossref issued date). FAMU-FSU College of Engineering. NSF DMR-2144007.

**Preprint:** R. Ricarte, ChemRxiv, posted 26 May 2026, CC-BY-NC-ND 4.0, DOI 10.26434/chemrxiv.15003901/v1.

**DOIs:** 10.1021/acs.macromol.6c00582 (journal); 10.26434/chemrxiv.15003901/v1 (preprint)

**Crossref check:** both found. Journal record: journal-article, Macromolecules, volume 59, issue 12, pages 6593-6619, single author Ralm G. Ricarte; relation has-preprint → 10.26434/chemrxiv.15003901/v1. Preprint record: posted-content dated 2026-05-26, author "Ralm Ricarte", relation is-preprint-of → 10.1021/acs.macromol.6c00582. The two are linked both ways. The local text is the preprint; page-level wording of the published version was not compared [unverified].

**Chemistry covered**
- Associative cross-link exchange (organometallic "associative" sense: bond formation before cleavage). Scheme 2B lists Ea^sm for esters, vinylogous urethanes (VU), imines, disulfides and boronic esters in kJ/mol and kBT at 393 K. The numbers are in a figure and are not in the extracted text (not stated in text). The text says esters and VUs are above the 10 kBT benchmark and imines, disulfides and boronic esters are below it.
- Catalysts: Lewis acids, Brønsted acids, strong organic bases, enzymes, neighbouring groups; VU and imine vitrimers need only excess amine; lithium salts shorten τ_xl by up to two orders of magnitude.
- Rubber relevance (strong): the history starts with Thiokol polysulfide rubbers (Stern and Tobolsky 1946, near-Maxwellian stress relaxation, Ea^rh ≈ 100 kJ/mol, matching S-S bond dissociation). Discusses polybutadiene (PB) vitrimers (boronic ester, the author's own work), PDMS vitrimers (boronic ester, imine, disulfide, boric acid), and vulcanization of commodity polyolefins, polydienes and polyesters as a route to vitrimers. Notes that disulfides can exchange dissociatively, so they do not strictly satisfy the vitrimer definition.

**Simulation approaches reviewed (Section 5.4.2 and Section 4.1)**
- Coarse-grained bead-spring MD with exchange added by either:
  - Hybrid MD/MC: after a set number of MD steps, find a free reactive site within a cutoff of a cross-link, compute ΔU_swap, accept with Metropolis P = min[1, exp(-ΔU_swap / kBT)]. Handles addition-elimination and metathesis. Caveats stated: it does not force a flat energy profile (ΔU_swap need not be zero), and MC steps can break detailed balance. Studies: Wu et al. 2019 (pioneer; small ΔU_swap → segmental relaxation dominates, large ΔU_swap → exchange dominates; KWW relaxation and Arrhenius τ_xl in the exchange regime), Perego and Khabaz (fast dynamics like permanent networks, WLF to Arrhenius crossover, creep), Cui/Tang (τ_xl⁰ temperature-independent for vitrimers; extra slow modes beyond sticky Rouse), Xia et al. (spatiotemporal heterogeneity, bimodal exchange, defect-rich regions stay active; gradient/block vitrimers; Rouse-like to reptation-like displacements with higher barrier; gelation threshold r ≥ 0.2 vs Flory-Stockmayer r = 0.17).
  - Three-body potential (Sciortino 2017): LJ-type pair bond, a repulsive three-body term proportional to λ (analogue of Ea^lcl) that cancels the energy of the second bond in a triplet, so the energy profile stays flat; addition-elimination only, needs excess reactive sites. Studies: Ciarella et al. 2018 (primary loops speed stress relaxation) and 2019 (density controls fragility), Zhao et al. (higher crosslink density and larger λ slow the slow regime), Liu and O'Connor (modified for associative plus dissociative, linear and nonlinear), Lin, Mei, Dutta, Schweizer, Sing 2025 (compared to SC-ECNLE; near-linear MSD; non-Arrhenius τ_xl near Tg).
- Kinetic Monte Carlo: Liu et al. 2024 for VU step-growth gelation parameterised by IR kinetics; deviations from Miller-Macosko at high conversion from intramolecular reactions.

**Theory approaches reviewed**
- Gelation: Flory-Stockmayer (p_c ≈ 1/N_w for vulcanization), Miller-Macosko (Cheng et al. added VU equilibrium constants), coherent state field theory (suggested).
- Rubber elasticity: affine (G_N = ν_xl RT = ρRT/M_xl), phantom (G_N = (1 - 2/f) ρRT/M_xl), Flory dangling-end correction.
- Swelling: Flory-Rehner, Bray-Merrill (crosslinked in solvent); χ assumed equal to the uncrosslinked polymer, which the author flags as possibly unjustified.
- Tg: Gibbs-DiMarzio; ECNLE and self-consistent SC-ECNLE (Tg ~ f_xl^0.5 for slow exchange; fast exchange plasticizes).
- Linear viscoelasticity: transient network theory (Green-Tobolsky, Maxwell G(t) = ν_xl RT exp(-t/τ_r)); KWW G(t) = ν_xl RT exp[-(t/τ_r)^β], ⟨τ_r⟩ = Γ(1/β) τ_r / β; generalized Maxwell model; sticky Rouse and sticky reptation; the Ricarte-Shanbhag inhomogeneous Rouse extension (needs full diagonalisation unless τ_xl N_xl² ≫ τ_seg N²); SC-ECNLE τ_xl = τ_xl^loc + t_s τ_seg^perm; renormalized bond lifetime; Evans and Sokolov τ_xl = A exp(ΔH_a/RT) + B τ_seg; Ge-Evans decoupling exponent ν (0.6 to 0.9); Carden-Sokolov steric factor p (10⁻¹¹ to 10⁻⁸); Ea^rh ≈ Ea^lcl + Ea^diff for nucleophile-driven exchange; Sciortino fluctuating-elasticity model for power-law relaxation; for dissociative networks also modified time-marching and phantom chain hopping.

**Recommended analyses of linear viscoelasticity (as argued in the Perspective)**
- Treat the spectrum as two regimes: fast segmental (WLF) and slow exchange-controlled (Arrhenius). Do TTS separately on each regime with its own horizontal shift factors (a_T^high WLF with prepolymer C1, C2; a_T^low Arrhenius) instead of one set; report vertical factors b_T = ρ_ref T_ref / (ρ T).
- Prefer a generalized Maxwell discrete spectrum {G_i, τ_i} and its moments over a single TNT or KWW fit. In the worked example (sticky Rouse, N = 100, N_xl = 10, Ea^lcl = 60 kJ/mol, T = 130 °C) KWW fitted to t > 0.001 τ_xl gives ⟨τ_r⟩ = 14 τ_xl, β = 0.65 and still relaxes too fast; a four-mode GMM overlaps the slow regime with first moment 17 τ_xl; τ_1 = N_xl² τ_xl = 100 τ_xl. GMM spectra are not unique.
- Estimate τ_xl from the low-frequency G'-G'' crossover (τ_xl ~ ω_cross⁻¹) or from BDS; note that for sticky Rouse the crossover lies below τ_xl⁻¹.
- Check terminal behaviour explicitly: G' ~ ω², G'' ~ ω, J ~ t in creep, exponential G(t). Report extended power laws when the terminal regime is not reached.
- Test Ea^rh against Ea^seg + Ea^lcl (Eq. 17) and plot ln(τ_xl/τ_seg) vs Tg/T (Eq. 18); experiments contradict both, so report Ea^rh vs crosslink density.
- Topology: compare G_N with ANT/PNT for the elastically effective fraction; gel fraction by Soxhlet; swelling with Flory-Rehner/Bray-Merrill with care about χ; Tg and its breadth (fast scanning calorimetry shows two events).
- G_N vs temperature distinguishes dissociative networks (G_N drops, thermal gel-sol at T_gel) from vitrimers (G_N constant or rising).

**Key quantitative statements**
- Polysulfide rubber Ea^rh ≈ 100 kJ/mol (Stern and Tobolsky). PB vitrimers: only about 50 % of crosslinks elastically effective at very high crosslink density; vulcanized vitrimers 20 to 100 % effective; boronic ester n-butyl methacrylate as low as 27 % of PNT.
- Boronic ester vitrimers Ea^rh 18 to 50 kJ/mol, imines 48 to 86 kJ/mol (Carden compilation); bis-boronic ester 30 to 60 kJ/mol; PDMS/PB blends ≈60 kJ/mol regardless of composition.
- Epoxy vitrimer fast regime Ea^eff ≈160 kJ/mol (Meng); 10 mol % excess amine cuts τ_xl by 6 to 8 times; disulfide PDMS G' ~ G'' ~ ω^0.2 over fifteen decades; Fang/Winter tan δ frequency-independent at T = 76 °C in an epoxy vitrimer.

**Relevance for CAPS**
- Not a primary data source; use it as the specification for a CAPS stress-relaxation / LVE analysis module: Maxwell, KWW (with ⟨τ⟩ = Γ(1/β)τ/β), discrete GMM spectrum with moments, two-regime TTS with separate shift factors, crossover-based τ_xl, Arrhenius/WLF fits with Ea and C1/C2, terminal-slope checks.
- Simulation methods to offer: hybrid MD/MC swap (Metropolis on ΔU) and the Sciortino three-body swap; both on Kremer-Grest networks. For rubber chemistry, disulfide (polysulfide rubber) and PB vulcanizate vitrimers are the closest targets.
- First thing to reproduce: the sticky Rouse G(t) of Fig. 7A (N = 100, N_xl = 10, Ea^lcl = 60 kJ/mol, 130 °C, τ_xl from Eq. 12) as a pure analysis test, then fit it with KWW and GMM to recover ⟨τ_r⟩ = 14 τ_xl, β = 0.65 and first moment 17 τ_xl.

---

### F4. Ciarella, Biezemans, Janssen 2019, Fragility of vitrimeric polymers

**Local file:** `Ciarella2019_PNAS_fragility_arXiv1910.00468.pdf` (arXiv:1910.00468v1 draft, compiled 2 October 2019)

**Citation (published):** S. Ciarella, R. A. Biezemans, L. M. C. Janssen. "Understanding, predicting, and tuning the fragility of vitrimeric polymers." Proc. Natl. Acad. Sci. U.S.A. 2019, 116 (50), 25013-25022. Eindhoven University of Technology.

**DOI:** 10.1073/pnas.1912571116

**Crossref check:** found. Journal-article, Proceedings of the National Academy of Sciences, volume 116, issue 50, pages 25013-25022, authors Simone Ciarella, Rutger A. Biezemans, Liesbeth M. C. Janssen, issued 2019-11-25. All match. The local file is the arXiv draft (PNAS placeholder DOI "pnas.XXXXXXXXXX"); numbers were not compared with the final PNAS text [unverified].

**Chemistry:** none specific. Generic vitrimer (bond-swapping end groups); the barrier ΔE_swap stands for catalyst control of exchange. No rubber chemistry; it is a star-polymer network model.

**Model level and force field**
- Coarse-grained star polymers, reduced units L0 = 1, ε = 1, m0 = 1, T in ε/kB.
- Each star: n = 25 segments (1 centre + 8 arms of 4 segments). Bonded neighbours: harmonic springs, equilibrium length L0 = 1; spring constant not stated.
- Non-neighbour segments: WCA with cutoff σ = 0.9 (L0).
- End segments of types A and B: A-B attraction via a generalized 20-10 Lennard-Jones potential, equilibrium length 2^(1/10) σ, depth ε ≫ kBT. Methods also gives r_min = 2^(1/10) σ ≈ 0.54, which implies a different σ for the A-B pair than the 0.9 WCA σ (see inconsistencies); the A-B σ value is not stated here (refers to Refs. 37, 42, 74).
- Three-body swap potential (Sciortino 2017): v3_ijk = λ ε v_ij^bond(r_ij) v_ik^bond(r_ik) over A-B-A and B-A-B triplets; ΔE_swap = ε(λ - 1). Runs at ΔE_swap = 0 (barrierless swaps, Tv → 0) and ΔE_swap = ∞ (no swaps, Tv → ∞). Equilibration always with ΔE_swap = 0.
- No all-atom force field.

**Exchange method and parameters**
- Algorithm: continuous three-body potential in MD (no MC). Associative, 1-to-1 bonding enforced, A-B-A / B-A-B swaps initiated by the unbonded majority ends.
- Barrier: ΔE_swap = ε(λ - 1), set to 0 or ∞; no physical kJ/mol value. No attempt frequency or acceptance rule (deterministic MD). Cutoffs of the 20-10 and three-body terms: not stated in this text.
- Temperatures: equilibration starts at T = 0.2, then stepwise cooling; Tg ranges from ≈10⁻¹ (high density) to ≈10⁻⁶ (lowest density).

**System size and run length**
- N1 = 600 stars with 7 A ends + 1 B end, N2 = 300 stars with 1 A + 7 B; 900 stars, N = 22500 segments; 2700 A-B bonds and 1800 free A ends when fully bonded.
- NVT, cubic periodic box; density set by box size (ρ from 0.91 to 2.89 shown).
- Production runs n_t = 10⁶ to 10⁸ steps, dt = 10⁻³; 2 to 5 independent runs per state point. Aging checked by comparing Fs(k,t) at the start and end of runs (Fig. S1).

**Software:** HOOMD-blue. MCT solved with Fuchs' algorithm (own code; not named).

**Properties computed and how**
- Self and collective intermediate scattering functions Fs(k,t), F(k,t) (monomer-averaged), static structure factor S(k), radial distribution function g(r) and its first-peak position r0.
- Structural relaxation time τ from Fs(k0, τ) = e⁻¹; k-resolved τ(k) and τc(k) from Fs and F.
- Operational Tg: the temperature where τ(Tg) = 10^4.5 τ0, with τ0 = 0.06 (m0 L0²/ε)^(1/2) at a reference state (T = 0.2; ρ = 1.58 in the Approach section, ρ = 1.5 in Methods).
- Fragility m = d log10 τ / d(Tg/T) at T = Tg, from polynomial fits to the Angell plots (simulation and MCT).
- Microscopic MCT (overdamped, site-averaged, form factor w(k) = 1) with simulated S(k) as the only input; 100 wavenumbers from k = 0.2 to 40; time step 10⁻⁶ doubling every 32 steps. Schematic MCT with S(k0;T) given by a ν-dependent growth law and C1 = C2 = 1, compared with the F2 model.
- Not computed: stress relaxation G(t), Arrhenius Ea of exchange, Tv, modulus, creep, fracture, healing, topology/loops.

**Key quantitative results**
- ρ ≥ 2.0: fragile (super-Arrhenius); 1.5 < ρ < 2.0: strong (Arrhenius); ρ < 1.5: superstrong (sub-Arrhenius), m ≈ 1 for ρ ≤ ρc ≈ 1.5; m rises monotonically above ρc. Fragility is k-independent (Fig. S2).
- Fragility along the glass line is practically the same for ΔE_swap = 0 and ∞ (Fig. S4); swaps shift the glass line to lower T and higher ρ.
- MCT predicts liquid behaviour at all simulated state points; Tg,MCT ≈ 0.5 Tg for the most fragile samples; MCT reproduces ρc ≈ 1.5.
- r0 ≈ L0 at low density (bond-dominated), ≈ σ = 0.9 at intermediate density, ≈ 0.8 L0 at the highest densities.
- Superstrong behaviour traces to a weak temperature growth of S(k0).

**Relevance for CAPS**
- Main use: a reference implementation of the Sciortino three-body swap and a clean glass-dynamics analysis chain. Reproduce first: a smaller star-polymer vitrimer (same 25-segment stars and 2:1 composition, fewer stars to stay under 15,000 beads) at one high and one low density, computing S(k), Fs(k0,t), τ and an Angell plot.
- Analyses worth adding to CAPS regardless of exchange method: ISF (self and collective), S(k), τ(k), operational Tg from a τ threshold, fragility index from an Angell plot.
- Whether LAMMPS ships a Sciortino swap potential was not checked here [unverified]; implementing it in the CAPS engine is the safe assumption.
- No rubber chemistry.

---

### F5. Shanbhag and Ricarte 2026, Fractional viscoelasticity in transient unentangled polymer networks

**Local file:** `ShanbhagRicarte2026_arXiv2608.04146_fractional.pdf` (arXiv:2608.04146v1 [cond-mat.soft], 4 Aug 2026; manuscript dated 6 August 2026)

**Citation:** S. Shanbhag, R. G. Ricarte. "Fractional Viscoelasticity in Transient Unentangled Polymer Networks." arXiv:2608.04146 (2026). Florida State University and FAMU-FSU College of Engineering.

**DOI:** 10.48550/arxiv.2608.04146

**Verification:** Crossref returns not found (arXiv DOIs are DataCite DOIs). DataCite: title matches, creators Shanbhag, Sachin and Ricarte, Ralm G., publication year 2026, v1 submitted 2026-08-04, updated 2026-08-06, no related identifiers. A Crossref title search found no published journal version as of 2026-10-22.

**Chemistry:** theory paper. The benchmark data are an imine-crosslinked polystyrene vitrimer, PS-v-8-1.2 (8 kDa precursor, diamine:aldehyde molar ratio 1.2, i.e. excess diamine; from the authors' earlier work, their Ref. 17). Other systems cited for power laws: dioxaborolane PB vitrimers (β = 0.4), imine PS vitrimers (β = 1.6), dioxazaborocane methacrylic vitrimers (β ~ 0.4 to 0.7), disulfide PDMS vitrimer (G* ~ ω^0.2 over more than 15 decades). Rubber relevance: PB and PDMS examples, Thiokol polysulfide history; no simulation of rubber chemistry.

**Model level:** analytical/numerical single-chain model, no particles, no force field. FIRM (fractional inhomogeneous Rouse model): overdamped generalized Langevin equation with a fractional Gaussian noise memory kernel K(t) = (t - t')^(-α) / Γ(1 - α), 0 < α ≤ 1, beads with different friction ζ_n (sticky beads ζ_x), linear springs k = 3kBT/b². Special cases: FRM (homogeneous), IHR (α = 1), SRM, and FSRM (two-population approximation, Eq. 23).

**Exchange method:** not simulated explicitly. Cross-links are sticky beads with friction ratio δ = ζ_x/ζ. With τ_x = τ exp(Ea/RT): δ = (N_x/N)(τ_x/τ)^α = ρ_x exp(α Ea/RT). The attempt time is set by the segmental time τ. Associative vs dissociative not distinguished by the model.

**System size / run length:** model chains N = 10 to 320 (N = 100, α = 0.8 for FRM; N = 50, ρ_x = 0.2, δ = 10³; N = 20, ρ_x = 0.1, δ = 10²; random cross-link placement averaged over 50 replicas). Experimental fit: N ≈ 11 Kuhn segments (720 g/mol per PS Kuhn segment), 2.1 effective cross-links per chain, ρ_x = 0.189.

**Software:** pyReSpect (relaxation spectrum and Prony fits), scipy.signal (partial fractions). Own code otherwise; data on request.

**Properties computed and how**
- Relaxation modulus: G(t)/G0 = (1/N) Σ_p E_α²(-(t/τ_p)^α) (one-parameter Mittag-Leffler), τ_p = λ_p^(-1/α) from the generalized eigenproblem C v_p = λ_p Γ v_p. Gives t^(-α/2) before and t^(-2α) after the terminal time τ_R = τ N^(2/α) (or τ_x N_x^(2/α) for crosslinked chains).
- Viscosity η = ∫G dt, approximated by η̂ = G0 τ I(α) ζ_R(2/α) N^(2/α - 1), I(α) ≈ (0.6488α² - 0.7638α + 0.3668)/(α - 0.5); η diverges for α ≤ 0.5.
- Stress relaxation fits to PS-v-8-1.2 at 130, 145, 160, 175, 190 °C: three fitting modes (α and Ea free; Ea only with α = 1, i.e. IHR; α only with Ea = 35 kJ/mol). Bare-chain parameters from unentangled PS13K at 104.7 °C, shifted with the Rouse molar mass dependence and WLF (C1 = 7.57, C2 = 115.7 K); G0 ≈ 1.54 MPa; τ = 7.2 ms (adjusting for Tg 87.6 °C vs 103 °C). Uncertainty from 100 replicas.
- Arrhenius analysis of τ_max and τ_min from the FIRM spectrum vs 1/T (linear regression), and E_WLF(T) = 2.303 R C1 C2 T² / (C2 + T - T_ref)² with T_ref = 418 K.
- End-to-end vector autocorrelation P(t), bead and centre-of-mass MSD (ρ_cm ∝ t^α), dielectric permittivity as a sum of Cole-Cole terms, G'(ω), G''(ω) via a fitted discrete spectrum.
- Comparison with the fractional Maxwell model (springpots); Prony-series approximation of the memory kernel (M = 10 modes for 10⁻² to 10², M = 12 for 10⁻⁶ to 10²).
- Not computed: Tv, Tg, fragility, creep, fracture, healing, topology analysis.

**Key quantitative results**
- Best-fit α = 0.85 ± 0.03 at all temperatures; α = 1 (IHR) gives clearly worse fits with an abrupt exponential terminal decay.
- Fitted Ea rises from ≈26 kJ/mol at 130 °C to 42 kJ/mol at 190 °C: Ea(T) = 155.5 - 5.1952 × 10⁴ / T (kJ/mol, K). With α = 1 the Ea values are 2 to 3 kJ/mol higher. Fixed Ea = 35 kJ/mol gives poor fits.
- Simplified system (N = 10, ρ_x = 0.2, α = 0.85): E_min ≈ 193 kJ/mol from τ_min; E_WLF ≈ 270 kJ/mol at 130 °C to 140 kJ/mol at 190 °C; E_max = 228 kJ/mol for constant Ea = 35; E_max = 106 kJ/mol with Ea(T), vs Eq. 31 estimate 108 kJ/mol at T = 433 K and the TTS value 117 ± 5 kJ/mol.
- FMM fit to FRM (N = 20, α = 0.4): c_a/G0 = 25.58 s^0.8, c_b/G0 = 0.14 s^0.2, a = 0.8, b = 0.2; shapes differ even with the same terminal slopes.

**Relevance for CAPS**
- No particle simulation to reproduce. Useful as an analysis option in the CAPS stress-relaxation module: fit G(t) from MD (e.g. Green-Kubo stress autocorrelation of a Kremer-Grest vitrimer) with Mittag-Leffler modes, compare α against the MSD exponent of the sticky beads measured in the same MD run (the paper proposes exactly this link).
- Reproduce first: the FRM G(t) for N = 100, α = 0.8 (power laws t^-0.4 and t^-1.6) and the viscosity scaling η ~ N^(2/α - 1); both need only an eigenproblem and Mittag-Leffler evaluation.

---

### Inconsistencies and gaps found

1. Ciarella 2019 (arXiv draft): the reference state for τ0 is ρ = 1.58, T = 0.2 in the Approach section but ρ = 1.5, T = 0.2 in Materials and Methods. Methods states the A-B bond length as 2^(1/10) σ with σ = 0.9 (≈0.96) but then r_min = 2^(1/10) σ ≈ 0.54, which needs σ_AB ≈ 0.5; the A-B σ is not given. Fig. 2(d) caption says the MCT axis is limited to Tg while the text says ~0.5 Tg,MCT. Harmonic spring constant not stated.
2. Garcia Villatoro 2026: Crossref drops the given name of M. Olvera de la Cruz. FENE k and R0, temperature, timestep, bead counts and reaction probability are not in the main text; the box volume is garbled in the text extraction. References [14] and [19] are the same paper. The text calls ref [6] (Ezzeddine, experimental PB networks) a prior simulation study. NDS is called "spectrometry" in the abstract and "spectroscopy" in Results A.
3. Guo and Zhang 2026: the Bell rupture and regeneration implementation in LAMMPS is not described (no fix named), and the re-formation partner rule/distance is not stated. The manuscript is dated 2 August 2026 but posted 3 September 2026. Figure 1 caption labels both rupture probabilities "Pr".
4. Ricarte 2026: preprint title uses an en dash, the journal record a box-drawing character; otherwise consistent (Macromolecules 59(12) 6593-6619, also cited that way by Garcia Villatoro). Scheme 2B activation energies exist only in the figure.
5. Shanbhag and Ricarte 2026: the arXiv DOI is not in Crossref (DataCite DOI; expected). Fig. 7 caption has a typo ("Ea and the only free parameter"). The FIRM Ea (26 to 42 kJ/mol, rising with T) differs from the 117 ± 5 kJ/mol TTS activation energy; the paper reconciles them through E_max ≈ E_min + a - 2b/T (my arithmetic check: 193 + 155.5 - 2 × 51952/433 ≈ 108.5 kJ/mol, matches the stated 108).
6. Numbers I re-checked [derived]: Ciarella composition gives 4500 A and 2700 B ends, so 2700 bonds and 1800 free A, as stated; Guo-Zhang bead totals 160 × 10 + 160 × 2 = 1920 backbone and 160 × 2 = 320 linker beads, as stated; Shanbhag Eq. 28 gives 26.6 kJ/mol at 130 °C and 43.3 kJ/mol at 190 °C (paper: ≈26 and 42); E_WLF at 130 and 190 °C computes to 268 and 139 kJ/mol (paper: ≈270 and 140).
