# Molecular simulation of vitrimers and self-healing dynamic covalent networks

Literature review for a CAPS bond-exchange module, written 2026-10-22. The module design built on it is in
`caps_module_proposal.md`. Full notes on each paper are in `paper_notes.md`: citation, Crossref check, chemistry, model,
force field, exchange method and parameters, system size, run length, software, how each property was analysed, and key
numbers. The PDFs are in `~/project/ff-sources/papers/vitrimers/`, and the README there maps file names to DOIs.

Conventions: numbers are copied as printed in the papers. "[unverified]" marks a value that could not be checked against
the paper, often because an equation was an image that did not survive text extraction, or because the value is in an SI
that was not available. When the local copy is a preprint, the published version's wording was not compared.

## 1. How the papers were found

- **Starting set.** An earlier pass had collected 28 PDFs.
- **Google Scholar.** The user's Scholar tab in Chrome was not reachable: it is outside the tab group this session can
  read. A search in a new tab of this session's own was answered with a "not a robot" check, so the tab was closed and
  Scholar was dropped.
- **Searches.** OpenAlex was searched by title and abstract for 2018 to 2026, newest first. Queries covered vitrimer MD,
  coarse-grained vitrimers, bond swap, dynamic covalent networks, disulfide and transesterification exchange, ReaxFF,
  covalent adaptable networks and ENR vitrimers. Crossref was used to check every DOI.
- **Downloads.** Ten more open-access PDFs were added: two straight from arXiv, eight through Chrome.
- **Paywalled.** 42 entries are listed in `~/project/ff-sources/papers/TO_DOWNLOAD.md` under "Vitrimers
  (2026-10-22)". This pass added twelve of them and removed Rovigatti 2018, which is
  now downloaded.

The final set is 38 papers: 31 primary simulation, theory, QM or ML papers, 4 reviews or perspectives, and 3 experimental
papers kept as validation targets (one of them, Sun 2026 in JACS, promises molecular dynamics in its title but has no
simulation).

### Papers by method

A paper is counted under its main method.

| Method | Count | Papers |
|---|---|---|
| Coarse-grained MD, three-body swap potential | 5 | Ciarella 2018 PRL; Ciarella and Ellenbroek 2019 Coatings; Ciarella, Biezemans, Janssen 2019 PNAS; Rovigatti 2018; Lin 2026 JCP |
| Coarse-grained MD with MC bond swaps (incl. LAMMPS bond/swap) | 3 | Branham-Ferrari 2026 Nat. Phys.; Rispo Constantinou 2026; Mamat 2026 |
| Coarse-grained MD, bind/unbind or rupture (dissociative) | 3 | Stukalin 2013; Wilcoxson 2026 (DyBond-v2); Guo and Zhang 2026 (Bell rupture) |
| Coarse-grained MD with fix bond/react exchange | 1 | Garcia Villatoro 2026 |
| Pure MC algorithms | 2 | Rao, Xia, Ni 2024; Oyarzun and Mognetti 2018 |
| All-atom fixed-charge MD with exchange | 7 | Yang 2016 (PCFF); Singh 2020, Singh 2024 (CVFF); Uyama 2026 (PCFF); Wick 2025 (OPLS-AA); Yuan 2026 (OPLS-AA); Sharma 2026 (PCFF, with experiment) |
| All-atom with physical (hydrogen) bonds only | 1 | Chen 2023 (COMPASS, NR) |
| Reactive force field (ReaxFF) | 3 | Zheng, Varshney, Vashisth 2025 (accelerated ReaxFF exchange); Adhikari 2026 and Hou 2026 (cure only) |
| QM barriers (DFT, ab initio MD) | 3 | Irigoyen 2022 (disulfide); Hamzehlou 2022 (vinylogous urethane); Teotonico 2024 (boronic ester) |
| ML property models | 2 | Golbabaei 2026 (healing efficiency); Zheng 2025 Digital Discovery (Tg) |
| Theory | 1 | Shanbhag and Ricarte 2026 (fractional Rouse) |
| Reviews and perspectives | 4 | Karatrantos 2024; Li 2023/2024; Holoman 2026; Ricarte 2026 |
| Experiment only (validation targets) | 3 | Arbe 2023 (imine polyisoprene); Danikas 2026 (dioxaborolane PBMA); Sun 2026 (polyacrylate) |

The table splits a few papers that cross categories. Rispo Constantinou and Stukalin are also theory papers. Branham-Ferrari,
Garcia Villatoro, Teotonico and Sharma also contain experiments. Irigoyen also runs classical MD.

## 2. The simulation methods

### 2.1 Coarse-grained, associative exchange

**Three-body swap potential.** The potential is Sciortino's (Eur. Phys. J. E 40, 3, 2017; paywalled, cited by every
CG paper here). It works on top of a short-ranged reactive pair potential, usually a generalised 20-10 Lennard-Jones term
between reactive end beads with a depth far above kT; Ciarella 2018 uses eps = 100 kT and sigma = 0.5 nm. The three-body
term

    v3 = lambda eps v_ij v_ik

cancels the energy of a second bond on the same site, so a swap costs eps(lambda - 1). Swaps happen in continuous MD
with no Monte Carlo step, and the barrier is set by lambda.

| Paper | System and size | Barrier and run | Main result |
|---|---|---|---|
| Ciarella 2018 PRL | About 1500 8-arm stars, about 50,000 beads, phi about 0.3 | beta dE_sw from 0 to 90; HOOMD-blue (RevCross) | Networks that allow primary loops relax ten times faster than defect-free ones, at a similar number of swaps |
| Ciarella and Ellenbroek 2019 Coatings | Two slabs of the same model, 300 adhesion runs | | Healing (0.15 tau_sw) is about three orders faster than flow (13.2 tau_sw) |
| Ciarella 2019 PNAS | 900 stars, 22,500 segments | | Density moves fragility from fragile to superstrong; mode-coupling theory reproduces this |
| Rovigatti 2018 | 800,000-particle binary patchy mixture | lambda = 1 | Collective relaxation time independent of wavevector; bond lifetime scales as 1/(x_fb - x) |
| Lin 2026 JCP | PnBA-like semiflexible chains | lambda = 1.1 and 1.25 at T* = 0.57 to 0.70; LAMMPS | Exchange speeds penetrant diffusion by relaxing mesh confinement, not by local hopping |

The method has two weaknesses (Holoman 2026; Ricarte 2026). It is monovalent: one partner per site. And it needs
excess free reactive ends.

**MC bond swap in hybrid MD/MC.** A move replaces bonds A1-A2 and B1-B2 by A1-B2 and B1-A2 when A1 and B1 are within a
cutoff, accepted by Metropolis on the change in bonded plus pair energy. This is LAMMPS `fix bond/swap` (Sides et al.
2004). Its manual page, read 2026-10-22, sets these requirements: molecule IDs, `special_bonds lj 0 1 1`, and no dihedral
or improper potential anywhere in the simulation.

- **Branham-Ferrari 2026** uses bond/swap on Kremer-Grest telechelic strands: FENE K = 30, R0 = 1.5, swaps between end
  beads within 1.5 sigma, about one attempt per timestep across the processors. From the geometry it obtains a swap
  barrier E* = 3.24 eps, and E* falls as the non-bonded sigma between exchanging beads is reduced (0.985, 0.96, 0.92).
  The main result: bond relaxation is Arrhenius at high T and turns super-Arrhenius once the exchange time falls below the
  alpha time. The bond decorrelation then fits a two-step form, and tau_bond is about tau_ex + c tau_alpha.
- **Rispo Constantinou 2026** uses its own Python MC (the PyPI package "MC-exchange") on 2000 Kremer-Grest chains of 20
  beads. Two move types, swap and shift, are tried every 1000 MD steps; the lower-energy outcome is kept, then accepted
  by Metropolis. The authors note this satisfies global, not detailed, balance.
- **Mamat 2026** scans the exchange interval t_exch from 100 to 1000 tau with a Metropolis-Hastings swap on a
  Kremer-Grest network that contains defects. Results:
  - volumetric Tg is 0.43, independent of t_exch;
  - tau_bond is Arrhenius with apparent activation 0.91 to 1.23 (reduced units);
  - tau_bond grows linearly with t_exch;
  - time-temperature superposition fails in creep and MSD for every finite t_exch.

**Detailed balance for multivalent sites.** Rao, Xia and Ni 2024 give an acceptance rule with a bias term:

    acc = min[1, v_u1/(v_u0 + 1) exp(-b dG0)]

For hybrid MD with an extra activation barrier it becomes

    acc = exp(-b dGa) min[1, (p'/v'u)/(p/vu) exp(-b dG0)]

The rule is checked against analytic grand-canonical equilibria. It is what a CAPS swap needs when a site can carry more
than one bond (multifunctional linkers, polysulfide sulfur).

### 2.2 Coarse-grained, dissociative bonds and rupture

**Stukalin 2013** (DL_POLY 2.12 hybrid MD/MC). Stickers bond and unbond by Metropolis on a sticky FENE bond:

    U_SB = U_FENE - U_FENE(b0) - eps

The protocol cuts the network, waits, and lets it heal; the measure is bridges across the cut. The paper also sets out
scaling theory for healing against self-adhesion.

**Wilcoxson 2026** (DyBond-v2, HOOMD-blue). Bonds form with probability

    P_on = k_on n dt exp(-U_b/kT)

and break with P_off = k_off n dt, with disjoint pair matching. It is checked against an exact four-reaction partition
function to within 0.02. The authors state that it has no associative swap move.

**Guo and Zhang 2026.** A 2240-bead Kremer-Grest network in which every bond can rupture with the Bell rate

    k_r(F) = nu0 exp(-Ea/kT) exp(gamma0 F)

Weak linkers that rupture and re-form delay backbone failure, and dispersion attraction between linkers trades
extensibility for strength. The model is fully specified, but the paper does not say how it was coded in LAMMPS.

**Oyarzun and Mognetti 2018.** A configurational-bias MC bind/unbind move for single-chain nanoparticles, checked against
exact loop statistics.

### 2.3 Coarse-grained, template exchange

**Garcia Villatoro 2026** runs Kremer-Grest A3 stars and B2 linkers with two reversible LAMMPS `fix bond/react`
templates for sulfonium transalkylation, initiated within r = 0.97 sigma. As networks exchange, primary loops grow at low
concentration and the gel dissolves; the paper's conclusion is that entropy alone drives this. The paper does not state
FENE constants, temperature, timestep or reaction probability.

### 2.4 All-atom, fixed-charge force fields

All six all-atom exchange papers use the same three steps:
1. build and cure a network by distance-triggered bonding;
2. find exchange candidates within a cutoff;
3. accept or reject the exchange, then relax locally.

They differ in the acceptance rule.

| Paper | Chemistry, FF | Exchange rule | Size | Notable output |
|---|---|---|---|---|
| Yang 2016 | DGEBA + tricarballylic acid, transesterification; PCFF | Distance cutoff; accept if the new bond energy is lower; extra MC test p ~ exp(-Ea/RT) with Ea = 68.18 kJ/mol normalised to 1 at 450 K (source of Ea not stated) | 10,860 atoms (two films), 450 and 420 K, one BER iteration about 0.1 ns | Welding complete at 8.73 tau_BER; stress relaxation single-exponential tau = 4.2 ns; recovered Ea 67.67 kJ/mol |
| Singh 2020 | DGEBA + aminophenyl disulfide; CVFF | fix bond/react S-S metathesis within 4.12 A; logistic p(T) centred on an assumed Tv = Tg = 403 K | 68,608 atoms | CTE 41.1 uK^-1 below Tv; E 4.5 GPa (2.6 GPa in experiment); healing after CNT pull-out |
| Singh 2024 | Same; CVFF, 12 A cutoff, 1 fs | fix bond/react with fixed p = 0 to 1; creep loop NPT 5 ps / NVE 4 ps | 68,608 atoms, density 1.159 g/cm3 at 300 K | All three creep stages at 600 K, 500 MPa, p = 1 |
| Uyama 2026 | DGEBA + 4-aminophenyl disulfide; PCFF, with DFT (B3LYP, GRRM) | Bond scission past 0.64 times the QM maximum covalent length; five healing templates (direct S-S, S-C ortho/para) | About 27,000 atoms | 78 % of broken S-S reform in about 5 ns; E about 5.0 GPa, strength about 180 MPa |
| Wick 2025 | BPAGMA shape-memory methacrylate, transesterification; OPLS-AA (LigParGen) | Pair within 4.25 A chosen with weight exp(-0.2 R^2); two-stage minimisation; Metropolis on energy with the reverse-move probability | 25,039 atoms intact, 37,376 recycled | Tg 399 +/- 4 K; recycling efficiency 47.5 % (500 cycles) to 100 % (10,000) |
| Yuan 2026 | Bis-GMA/HEMA, ester-hydroxyl exchange; OPLS-AA | REACTER polymerisation, then exchange with a logistic p(T) | About 33,000 atoms (agent estimate; not stated) | Tg 394 K; healing as E_healed / E_original |

Sharma 2026 uses PCFF for density and modulus of a disulfide epoxy (E 2.74 GPa against 2.69 GPa measured). It does not
describe its exchange script. Its simulated hardener (an aminophenyl disulfide) also differs from the experimental one
(DTBA).

Chen 2023 is the only atomistic rubber study with healing. It models cis-polyisoprene carrying MAn-ATA hydrogen-bonding
groups, with COMPASS, and no covalent exchange. Healing is the toughness ratio after a gap of about 3 A is closed: 44.7 %
after 1 ns and 82.3 % after 30 ns at 498 K. G(t) comes from Green-Kubo with a multi-tau correlator, and G' and G'' from
piecewise-linear Fourier integration.

### 2.5 Reactive and machine-learned force fields

**Zheng, Varshney and Vashisth 2025** (accelerated ReaxFF, CHON2017_weak_bb, AMS 2024) push transesterification along
with bond-boost restraints:

    E_rest = F1 {1 - exp[-F2 (r - R12)^2]}

F1 and F2 were fitted by Bayesian optimisation; their table is in the notes. Temperature enters as

    F1 / F1(500 K) = F0 exp(-alpha Ea / RT),  Ea = 14.6 kcal/mol,  alpha = 0.2

The 14.6 kcal/mol barrier is from Bhusal et al., J. Phys. Chem. B 125, 2411 (2021), doi 10.1021/acs.jpcb.0c10403
(paywalled; which catalyst it applies to is unverified). The protocol is creep under 50 MPa while heating from 100 to
800 K at 4 K/ps.

The ML force field for imine exchange (Sun et al., Macromolecules 56, 9003, 2023) and Kamble 2022 (accelerated ReaxFF
fatigue reversal) are paywalled or blocked. Adhikari 2026 and Hou 2026 apply ReaxFF to curing only.

CAPS has no ReaxFF engine. ML potentials are named as the open direction by Holoman, Li and Karatrantos, but no paper
here trains one for exchange except Sun 2023, which has not been read.

### 2.6 Barriers from quantum chemistry

| Chemistry | Value | Level | Paper |
|---|---|---|---|
| Aromatic disulfide exchange (radical, 2+1) | dG = 12.02 kcal/mol | wB97XD/6-311++G(2df,2p)//6-31+G(d,p), gas phase | Irigoyen 2022 |
| Disulfide BDE (aromatic, aliphatic, sulfenamide) | 33.8 to 65.7 kcal/mol | same | Irigoyen 2022 |
| Aromatic disulfide exchange path (GRRM) | 150 to 180 kJ/mol | B3LYP/6-31+G(d) | Uyama 2026 |
| Vinylogous urethane transamination, protic | dE0 14.3 kcal/mol (exp. 14.1 +/- 1.4); dG 9.29 | wB97XD, PCM benzene | Hamzehlou 2022 |
| Vinylogous urethane, aprotic / with B(CF3)3 | 33.99 / 19.93 kcal/mol | same | Hamzehlou 2022 |
| Vinylogous urea | 12.97 kcal/mol (exp. 12.92) | same | Hamzehlou 2022 |
| Boronic ester (dioxaborinane), MeOH-mediated | 17.68 kcal/mol | wB97XD/6-311++G(2df,2p)//6-31+G(d) | Teotonico 2024 |
| Boronic ester, direct [2+2] metathesis | 29.15 kcal/mol | same | Teotonico 2024 |
| Transesterification (beta-hydroxy ester) | 14.6 kcal/mol | DFT, as quoted by Zheng 2025 | Bhusal 2021 (not read) |

The two disulfide numbers disagree by about a factor of three. They follow different pathways from different
references, and Irigoyen concludes that the barrier alone does not rank healing ability; radical formation and chain
mobility do. Hamzehlou's text lost its minus signs in extraction, so the signed intermediate energies need checking in the
PDF.

## 3. Exchange parameters across papers

| Paper | Algorithm | Barrier or rate | Attempt frequency | Cutoff | T |
|---|---|---|---|---|---|
| Branham-Ferrari 2026 | LAMMPS bond/swap, Metropolis | E* = 3.24 eps (from geometry) | about 1 attempt per step over the processors; about 15 tau per attempt per bond | 1.5 sigma | stepwise quench, values in figures |
| Rispo Constantinou 2026 | Own MC (swap and shift), Metropolis on FENE dU | none beyond FENE strain | every 1000 MD steps, one move per sticker per cycle | 1.12 sigma | T* = 1 |
| Mamat 2026 | MC swap, Metropolis-Hastings with N_f/N_r | none imposed; E_app 0.91 to 1.23 from fits | t_exch = 100 to 1000 tau | not stated | reduced, Tg = 0.43 |
| Ciarella 2018 | Three-body potential | beta dE_sw = 0 to 90 | continuous | 2.5 sigma (20-10 term) | 300 K equivalent |
| Lin 2026 | Three-body potential | lambda = 1.1, 1.25 (eps_AB = 30) | continuous | 1.2 | 0.57 to 0.70 (reduced) |
| Stukalin 2013 | Metropolis bind/unbind on sticky FENE | eps = 10 to 12 kT | every tau_LJ | R_FENE | 1.0 |
| Wilcoxson 2026 | DyBond-v2 bind/unbind | eps = ln(k_on/k_off) | every n steps | [r_min, r_max] | instantaneous |
| Guo and Zhang 2026 | Bell rupture plus re-formation | Ea,linker 73 to 104.55 (LJ); gamma0 = 1.025 | every 50 steps (0.25 tau) | not stated | 0.78 |
| Garcia Villatoro 2026 | fix bond/react, two reversible templates | not stated | not stated | 0.97 sigma | not stated |
| Yang 2016 | Distance, energy test, Arrhenius MC | Ea 68.18 kJ/mol | one iteration about 0.1 ns | not stated | 450, 420 K |
| Singh 2020 | fix bond/react, logistic p(T) | Tv = 403 K, w = 20 K, p(Tv) = 0.9 | not stated | 4.12 A (S-S) | 200 to 550 K scan; healing at 503 K |
| Singh 2024 | fix bond/react, fixed p | p = 0, 0.01, 0.1, 0.2, 1 | once per NVE segment | 4.12 A | 600 K |
| Wick 2025 | Distance-biased choice, Metropolis on energy | none (cycles stand in for time) | one swap per 20 ps cycle | 4.25 A | 398 or 429 K |
| Zheng 2025 | Accelerated ReaxFF bond boost | Ea 14.6 kcal/mol scaled by alpha = 0.2 | boosts of 10,000 steps | windows per pair (notes) | 100 to 800 K |

## 4. How the properties were measured

### 4.1 Bond and exchange dynamics

- **Exchange time tau_ex** comes from counted events, recursive ones included (Branham-Ferrari 2026), or as the average
  time per exchange per active unit (Yang 2016: tau_BER = 1.26 ns at 450 K).
- **Bond-partner autocorrelation.** Three definitions are in use:
  - the fraction of sites still bonded to their time-zero partner (Branham-Ferrari);
  - the connectivity-matrix correlation <H(t):H(0)>/<H(0):H(0)> (Lin);
  - the occupancy correlation <h_b(0) h_b(t)>/<h_b(0)> (Mamat).

  The fits are KWW (Lin), the 1/e time (Mamat, Rovigatti), or the two-step form
  f = 1/2 [1 + exp(-t/tau)] exp[-(t/tau_diff)^beta] at low T (Branham-Ferrari).
- **Segmental and alpha relaxation.** Two measures appear: the self-intermediate scattering function at the first
  structure-factor peak (tau where Fs = 0.2 in Branham-Ferrari, 1/e in Ciarella 2019), and the P2 autocorrelation of Kuhn
  vectors (Lin).

### 4.2 Stress relaxation and viscoelasticity

- **Equilibrium Green-Kubo.** G(t) = (V/kT) <sigma_ab(t) sigma_ab(0)>, averaged over the three shear components; the
  three-body virial is needed when that potential is used (Ciarella 2018, SI). Ciarella 2018 does not fit G(t). It defines
  tau_liq by G(tau)/G(0) = 1e-4, and separates solid from liquid by the long-time plateau. Chen 2023 uses a multi-tau
  correlator and converts to G' and G''.
- **Frozen topology.** Ciarella 2019 freezes the topology at a chosen healing time and computes Gxz and Gyz (across the
  cut) against Gxy (unaffected) to judge mechanical healing.
- **Step strain.** Yang 2016 applies 30 % at 1e9 s^-1, holds, and fits sigma(t) with a single exponential. Zhao 2022 applies
  300 % uniaxial strain (paywalled, quoted by Li 2023).
- **Fit forms in experiments.** KWW is used by Teotonico 2024 and Garcia Villatoro 2026 (beta 0.5 to 1.0); the 1/e time by
  Sharma 2026.
- **Ricarte 2026 recommendations.** Fit a generalised Maxwell spectrum and report its moments instead of a single KWW.
  Shift the fast segmental regime and the slow exchange regime separately (WLF against Arrhenius). Check terminal slopes
  (G' ~ w^2, G'' ~ w). Estimate tau_xl from the low-frequency crossover. Shanbhag and Ricarte 2026 add Mittag-Leffler
  (fractional Rouse) fits, with the fractional exponent alpha to be compared with the MSD exponent of sticky beads from MD.

### 4.3 Arrhenius analysis and Tv

- **Arrhenius fits.** Applied to tau_ex, tau_bond, tau_rel or viscosity:
  - Branham-Ferrari fits the high-T range and the full range separately;
  - Mamat fits ln tau_bond = ln tau_0 + E_app/T;
  - Yang 2016 recovers Ea from shift factors at two temperatures.
- **Tv in simulation.** Four routes appear:
  - the kink in the thermal expansion coefficient, Perego and Khabaz 2020 (paywalled; described by Karatrantos and Li);
  - imposed rather than detected (Singh 2020 sets Tv = Tg in p(T));
  - inferred from the onset of exchange (Zheng 2025: 300 to 500 K);
  - in reduced units, from bond relaxation time and Fs(q,t) (Xia et al. 2023, via Karatrantos; paywalled).
- **Tv in experiment.** Arbe 2023 extrapolates viscosity to 1e12 Pa s and finds the same Tv from a mesoscale X-ray
  expansion kink. Danikas 2026 locates Tv from shift factors.
- **Tg in simulation.** Hinge of specific volume or density against T (most papers). Adhikari takes the zero of the second
  derivative of a spline. Ciarella 2019 uses an operational threshold, tau(Tg) = 10^4.5 tau0.

### 4.4 Healing, welding and recycling

| Definition | Paper |
|---|---|
| Bridges across the cut plane N_br(t); bond profile B(z) across the interface | Stukalin 2013; Ciarella 2019 |
| Bonds across the interface normalised by saturation; penetration depth of active atoms | Yang 2016 |
| Modulus and yield stress of the welded sample against the fresh one, against t/tau_BER | Yang 2016 |
| Toughness ratio: area under the healed stress-strain curve over the original | Chen 2023 |
| Modulus ratio E_healed / E_original | Yuan 2026 |
| Fraction of broken bonds reformed, by pathway | Uyama 2026 |
| Recycling efficiency: mean stress at 28 % strain (2.0 to 2.5 ns hold), recycled over intact | Wick 2025 |
| Density recovery in the gap rho_m / rho_m,average | Chen 2023 |
| Experimental scratch or cut width: (w0 - w) / w0 | Sharma 2026; Golbabaei 2026 |

The definitions do not agree with each other, so any CAPS healing efficiency has to be reported with its definition.

### 4.5 Creep

- Singh 2024 identifies primary, secondary and tertiary creep in the stretch ratio, and tracks void fraction (alpha-shape,
  3.5 A probe) and S-S orientation.
- Zheng 2025 records strain against temperature during heating under stress, with the exchange count and free volume.
- Mamat 2026 measures the compliance J(t) at sigma0 = 0.1.
- Perego and Khabaz 2023 (paywalled) report J(t) ~ t^1/2, as quoted by Karatrantos.

## 5. Rubber and elastomer relevance

No paper read here simulates exchange in ENR, SBR, NR or a sulfur-cured diene rubber. The closest are:

- **Chen 2023.** Atomistic NR with hydrogen-bonding grafts (COMPASS); physical bonds, healing measured. This is the only
  rubber simulation in the set.
- **Arbe 2023.** Imine-crosslinked polyisoprene, an NR analogue, as experiment. Values:
  - Tg 216.5 to 221 K (plain PI 204 to 205 K);
  - Tv 260 to 274 K;
  - flow Ea 66 +/- 20 and 79 +/- 7 kJ/mol.
- **Danikas 2026.** Elastomeric dioxaborolane PBMA vitrimers, experiment. Values:
  - Tg 293 to 303 K;
  - Tv/Tg about 1.4 to 1.5;
  - Ea 80 to 94 kJ/mol;
  - plateau modulus 1.1 to 8.3 x 10^5 Pa.
- **Irigoyen 2022.** Bond energies of aliphatic disulfides (58 to 63 kcal/mol) and an aromatic disulfide exchange barrier.
  These are the closest inputs for the disulfide and polysulfide crosslinks of vulcanised rubber. There is no barrier for
  aliphatic exchange.
- **Ricarte 2026.** Recalls the polysulfide rubber stress relaxation of Stern and Tobolsky (near-Maxwellian, Ea about
  100 kJ/mol, about the S-S dissociation energy). Also discusses vitrimers made by vulcanising polydienes, and polybutadiene
  and PDMS vitrimers.
- **Guo and Zhang 2026.** A coarse model of a sacrificial weak S-S or Se-Se linker in an elastomer network.
- **Lin 2026 and Teotonico 2024.** Acrylate rubbers (PnBA): CG exchange and boronic ester barriers.
- **Garcia Villatoro 2026.** Argues that vitrimerised elastomers will come out stiffer than the original rubber because
  they have fewer loops.

For the user's fibre-rubber composites there is no simulation precedent for exchange at a fibre-rubber interface. Two
pieces come nearest: Sharma 2026 (carbon fibre in a disulfide epoxy, mechanics only) and Singh 2020 (healing after CNT
pull-out in a disulfide epoxy). The experimental ENR target is Imbernon et al. 2015 (ENR with a disulfide diacid), which
is paywalled.

## 6. The most useful papers for CAPS

1. **Branham-Ferrari et al. 2026, Nature Physics.** The simplest complete protocol: Kremer-Grest, LAMMPS bond/swap, a
   swap barrier derived from the geometry (E* = 3.24 eps), and a bond-decorrelation analysis that separates reaction from
   diffusion control. It is the first validation target for a CAPS engine and its bond/swap export together.
2. **Rao, Xia and Ni 2024, J. Chem. Phys.** The acceptance rule that keeps detailed balance with multivalent sites, plus
   an analytic test. It defines how a CAPS swap must be accepted.
3. **Mamat et al. 2026, ACS Appl. Polym. Mater.** A clean Kremer-Grest MC-swap study with Tg, Arrhenius bond lifetimes,
   creep and an explicit failure of time-temperature superposition. Numbers to reproduce in reduced units.
4. **Ciarella, Sciortino and Ellenbroek 2018, PRL.** The standard three-body benchmark, and G(t) from Green-Kubo with the
   three-body virial. Shows topology defects (loops) controlling relaxation.
5. **Ciarella and Ellenbroek 2019, Coatings.** A complete healing and adhesion protocol with three cheap metrics (bridges,
   B(z), frozen-topology G components) and three separated time scales.
6. **Yang et al. 2016, RSC Advances.** The best-quantified atomistic study: PCFF transesterification welding, 10,860
   atoms (under the CAPS 15,000-atom limit), with tau_BER, a stress-relaxation tau, a recovered Ea and modulus recovery.
7. **Singh and Sundararaghavan 2020 and Singh et al. 2024.** Disulfide exchange in a DGEBA epoxy through LAMMPS fix
   bond/react with CVFF. Closest to CAPS's existing bond/react export; gives the creep protocol.
8. **Wick, Peters and Li 2025, J. Polym. Sci.** OPLS-AA, the force field CAPS supports best. Its Metropolis
   transesterification move includes the reverse-move probability. It defines recycling efficiency and topology
   descriptors.
9. **Chen et al. 2023, Polymer.** The only rubber (NR) healing simulation, with COMPASS, a gap-healing protocol, the
   toughness-ratio efficiency, and Green-Kubo G(t) through G' and G''.
10. **Ricarte 2026, Macromolecules (perspective).** The specification for the analysis side: two-regime superposition,
    generalised Maxwell spectra, terminal checks, crossover-based tau_xl, and when KWW misleads.

The barrier papers (Irigoyen 2022, Hamzehlou 2022, Teotonico 2024) are not protocols. They are the source of most numbers
in the proposed exchange library.

## 7. Gaps and open questions

- **No catalyst in any atomistic simulation** (also noted by Karatrantos 2024). A catalyst appears only as a rate factor
  or a number of boosts.
- **Time scale.** Every atomistic study speeds exchange up by orders of magnitude. Results are comparable with
  experiment only as trends or in units of the exchange time.
- **Detailed balance.** Energy-lowering acceptance (Yang 2016), logistic p(T) (Singh 2020, Yuan 2026) and keeping the
  lower-energy of two outcomes (Rispo Constantinou 2026) are not detailed-balance moves. Rao 2024 and Wick 2025 show how to
  do it properly.
- **Barriers.** Conflicting for disulfides; missing for aliphatic polysulfide, siloxane and retro Diels-Alder exchange in
  the papers read; one value (14.6 kcal/mol) traced to a paper not yet read.
- **Rubber.** No exchange simulation of ENR, SBR, sulfur-cured NR or PDMS in this set. PDMS and VU (Rukmani 2024) and
  PDMS-like side-chain vitrimers (Xia 2023) are in the paywalled list.
- **Interfaces.** No simulation of exchange across a fibre-rubber or filler-rubber interface.
- **Inconsistencies found in the papers.** Recorded in `paper_notes.md`:
  - the Ea of Sharma 2026 against its own fit;
  - non-bonded distances in Branham-Ferrari;
  - the isomer of the aminophenyl disulfide in Singh 2020;
  - the A-B sigma in Ciarella 2019 PNAS;
  - Danikas's three different Ea values.

## 8. Paywalled and blocked papers

`~/project/ff-sources/papers/TO_DOWNLOAD.md`, section "Vitrimers (2026-10-22)", lists them in three blocks:
- the CG and atomistic protocols that matter most (Perego and Khabaz 2020, 2022, 2023; Yang 2015; Sun 2023; Rukmani 2024;
  Kuhzadmohammadi and Zhang 2025; Zhao 2022; Wu 2019; Zheng 2021);
- a second block of methods and models (Sciortino 2017, Smallenburg 2013, Vermeersch 2024, Xia 2023 and others);
- a third block found in this pass (Bhusal 2021, Ma 2024, Li 2026 JPCB, Zhao 2024, Wang 2026, Wu and Wen 2026, Lee and
  Jeong 2026, Perego and Khabaz 2021, Pandya 2023, Cui 2023, Li 2024 Langmuir, Chen 2026).

The same file lists open-access papers whose scripted download was refused (Ciarella and Ellenbroek 2022 SciPost,
Mangialetto 2024, Kamble 2022, Xia 2024, Ricarte and Shanbhag 2021, Imbernon 2015).
