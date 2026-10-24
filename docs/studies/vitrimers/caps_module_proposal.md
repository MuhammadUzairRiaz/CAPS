# Dynamic bond exchange in CAPS: module proposal

Status: proposal only, written 2026-10-22. No CAPS code has been changed. The literature behind each choice is in
`literature_review.md` in this folder. The PDFs are in `~/project/ff-sources/papers/vitrimers/`, with the paper table in
that folder's README. Values marked [unverified] need checking against the paper or its SI before they go into a data
file.

Working name: **Exchange**. That means CLI `caps exchange`, a Studio page called Exchange, and a library in
`data/exchange/`. In a vitrimer, bonds swap partners while the network keeps most of its connectivity, and the name says
that.

## 1. What CAPS already has

The survey of the repository at a3f86bcd found these pieces to build on.

| Need | What exists | Where |
|---|---|---|
| Network building | `react()` with find, apply, retype and relax cycles; targets by conversion, crosslinks, Mc or density; byproducts; inter-chain links; sulfur census | `core/src/react.cpp`, `core/include/caps/react.hpp`, `caps react` |
| Bond breaking inside a template | `break a b`, `move h x`, `delete n` keywords in the template text format | `react.cpp:251-310` |
| Reacting during dynamics | `during_md`: MD in `md_ps` segments, reacting between segments, local capped minimisation of each new site | `react.cpp:1120-1162`, `1398-1420` |
| Exchange chemistries as schemes | `transesterification`, `carbamate_exchange` (the library entry calls it the exchange behind PU vitrimers), `schiff_base`, `disulfide` (formation only), `epoxy_acid`, `enr_pbs_direct`, `siloxane_condensation`, `diels_alder` (forward only) | `data/reactions/library.json` |
| LAMMPS reaction export | `write_bond_react` writes pre, post and map files per variant, and `stem.in` with `fix bond/react ... prob P seed`, stabilisation and crosslink-target chunking. It writes only distance constraints (no Arrhenius constraint), and there is no CLI entry | `core/src/bond_react.cpp`, C API `caps_bond_react_export`, Python `Document.bond_react` |
| Own MD engine | Velocity Verlet, r-RESPA, SHAKE/RATTLE, Bussi, Langevin, Nose-Hoover, C-rescale, MTK; deformation, SLLOD shear, pulling; per-step hooks `each_step` and `after_step` | `core/src/dynamics.cpp` |
| Coarse-grained models | Kremer-Grest builder and force field usable by CAPS's own MD, LAMMPS KG writer, DPD | `core/include/caps/kremer_grest.hpp`, `dpd.hpp` |
| Relaxation analysis | MSD, diffusion, end-to-end and P2 autocorrelation with KWW fits | `core/src/properties.cpp` (`kww()` at line 724) |
| Stress analysis | Green-Kubo viscosity from the stress autocorrelation, NEMD viscosity | `core/src/mechanics.cpp:86-113` |
| Mechanics | Static and fluctuation Cij; tensile with modulus, yield and peak; creep; pull-off work of adhesion | `mechanics.cpp`, `caps tensile`, `caps pull` |
| Tg | Cooling scans and hinge fits | `fit_bilinear`, `fit_cooling_scans`, `caps tg` |
| Interfaces | Two-phase interface builder with gap, film and density | `core/src/interface.cpp`, `caps interface` |
| Missing bonded terms | `caps ff complete` and `--complete ml` (CAPS MACE v1) for new bond, angle and torsion types made by a reaction | `core/include/caps/ffcomplete.hpp` |
| Recipes and cluster jobs | Fixed stage list (`build ... react, relax, cbmc, equilibrate, md, thermal, analyze, export`); any CLI command runs as a job | `core/src/recipe.cpp`, `core/src/cluster_job.cpp` |

## 2. What is missing

1. **Topology change inside one MD run.** `run_dynamics` copies the force field once and builds one `Evaluator`, so the
   bond, angle and torsion lists and the pair exclusions cannot change during a run. `during_md` works around this by
   restarting MD between reaction steps. That is fine for curing every few ps, but too coarse for an exchange attempt
   every few steps.
2. **Bond-swap and bind/unbind moves.** CAPS has configurational-bias and GCMC moves (`cbmc.cpp`, `sorption.cpp`), but no
   Monte Carlo move on bonds.
3. **Exchange schemes.** Several are missing: disulfide metathesis (S-S + S-S), transimination or imine metathesis,
   boronic ester metathesis, vinylogous urethane transamination, silanolate siloxane equilibration and retro Diels-Alder.
   None of the existing schemes carries a barrier or a source paper.
4. **Export for exchange.** There is no `fix bond/swap` deck for bead-spring vitrimers, and bond/react exports carry only
   a fixed `prob`. The LAMMPS bond/react manual summary lists an Arrhenius constraint for the map-file Constraints
   section (`arrhenius A n E_a seed`, k = A T^n exp(-E_a/k_B T)) [unverified: read the full manual page before writing
   it].
5. **Stress relaxation and Tv analysis.** CAPS lacks:
   - a step-strain G(t) protocol;
   - a multi-tau stress correlator, and G' and G'' from G(t);
   - a bond-partner autocorrelation;
   - an Arrhenius fit in core (one exists only in the Studio's `CurveFit.cs`);
   - Tv from a thermal-expansion kink.
6. **Healing.** There is no cut, crack, void or weld builder, no count of bridges across a plane, and no definition of
   healing efficiency.

## 3. Proposed design

### 3.1 Exchange chemistry library (`data/exchange/library.json`)

There is one entry per chemistry. Each entry points to a reaction scheme, either existing or new in
`data/reactions/library.json`, and adds the exchange data with a source for every number:

```
{ "id": "disulfide_aromatic",
  "mechanism": "associative",           // associative | dissociative | physical
  "scheme": "disulfide_metathesis",     // one template, used in both directions
  "barrier": { "value": 12.02, "unit": "kcal/mol", "kind": "DFT dG",
               "method": "wB97XD/6-311++G(2df,2p)//6-31+G(d,p)",
               "source": "10.3389/fmats.2022.859482" },
  "alternatives": [ { "value": "150-180", "unit": "kJ/mol", "kind": "DFT path (B3LYP, GRRM)",
                      "source": "10.3390/polym18070861" } ],
  "experiment": { "Ea_relax": { "value": 48, "unit": "kJ/mol", "source": "10.1039/d6ra00851h" } },
  "capture": { "value": 4.12, "unit": "A", "pair": "S..S", "source": "10.1016/j.cplett.2020.137966" },
  "notes": "Uyama and Irigoyen differ by about 3x; do not mix them without a CAPS calculation." }
```

These are the first entries and the source of each number. Details and caveats are in the review, section 5.

| Entry | Barrier input | Source | Status |
|---|---|---|---|
| Transesterification, beta-hydroxy ester (TBD or Zn) | 14.6 kcal/mol (DFT, as used by Zheng 2025) | Bhusal et al. 2021, 10.1021/acs.jpcb.0c10403 | value taken from Zheng 2025; read Bhusal to see which catalyst it applies to [unverified] |
| Transesterification, epoxy-acid (atomistic MC rate) | Ea = 68.18 kJ/mol | Yang et al. 2016, 10.1039/c5ra26128g | Yang 2016 does not say where the value comes from [unverified] |
| Disulfide, aromatic | dG = 12.02 kcal/mol | Irigoyen et al. 2022 | DFT, gas phase |
| Disulfide, aliphatic and polysulfide (rubber) | BDE 58 to 63 kcal/mol; no exchange barrier | Irigoyen et al. 2022 | barrier missing; needs a CAPS QM job |
| Sulfenamide (BiTEMPS type) | 16 to 18 kcal/mol | cited by Irigoyen 2022 (Ruiperez 2018) | primary paper not read [unverified] |
| Vinylogous urethane, protic | dE0 = 14.3 kcal/mol (exp. 14.1 +/- 1.4) | Hamzehlou and Ruiperez 2022, 10.1002/pol.20220099 | minus signs lost in the text extraction; check the PDF |
| Vinylogous urea | 12.97 kcal/mol (exp. 12.92) | same | |
| Boronic ester (dioxaborinane), MeOH-mediated | 17.68 kcal/mol | Teotonico et al. 2024, 10.1039/d3py01047c | |
| Boronic ester, direct metathesis | 29.15 kcal/mol | same | |
| Dioxaborolane (elastomer, experiment) | flow Ea 80 to 94 kJ/mol | Danikas et al. 2026, 10.1021/acs.macromol.6c01360 | experimental only |
| Imine (transimination) | no simulated barrier read yet; experimental flow Ea 66 and 79 kJ/mol (PI networks) | Arbe et al. 2023, 10.1021/acsmacrolett.3c00586 | the ML barrier in Sun et al. 2023 is paywalled |
| Siloxane, Diels-Alder | no value from any paper read here | | gap; Vermeersch 2024 (retro DA, ReaxFF) is paywalled |

The library holds only numbers taken from the original papers. No parameter is copied from another program's
parameter files.

### 3.2 Exchange engine in CAPS's own MD (core)

The evaluator gets a mutable-topology path, so MD no longer has to restart:

- **Topology updates.** `Topology::swap_partner(bond, old, new)`, `Topology::break_bond` and `Topology::make_bond` update
  the bond list, the angle and torsion lists generated from it, and the 1-2, 1-3 and 1-4 exclusions and scaled pairs.
  They then mark the neighbour list dirty. A swap leaves the number of bonds unchanged, so the arrays keep their size.
- **ExchangeMove.** An `ExchangeMove` interface is called from `DynamicsOptions::after_step` every `every` steps. Four
  move types:
  - **Pair swap (Metropolis).** For bead-spring networks. Bonds A1-A2 and B1-B2, with A1 and B1 within `capture`,
    become A1-B2 and B1-A2. dU is the change in bonded plus pair energy. The move is accepted with
    min(1, exp(-dU/kT)), times an optional Arrhenius factor exp(-Ea/kT). This is what LAMMPS `fix bond/swap` does
    (Sides et al. 2004), and what Branham-Ferrari et al. 2026, Rispo Constantinou et al. 2026 and Mamat et al. 2026 use.
    Mamat adds an N_f/N_r Metropolis-Hastings correction for the number of forward and reverse choices.
  - **Multivalent swap.** For crosslinkers with more than one reactive valence: tri- and tetra-functional linkers, and
    polysulfide bridges. It uses the Rao, Xia and Ni 2024 acceptance:
    acc = exp(-b dGa) min[1, (p'/v'u)/(p/vu) exp(-b dG0)].
  - **Bind/unbind (dissociative).** For hydrogen-bond, ionic or Diels-Alder links. It has two forms:
    - DyBond-v2 (Wilcoxson et al. 2026): P_on = k_on n dt exp(-U_b/kT) and P_off = k_off n dt, with disjoint matching.
    - Stukalin et al. 2013: a sticky FENE bond.
  - **Template exchange (all-atom).** A reaction template with `break` and `form` lines goes through the same match code
    as `react`. A short local relaxation of the site follows. The move is accepted on energy with Metropolis and an
    Arrhenius factor from the library barrier. CAPS should use the reversible form below and report the acceptance and
    reverse-move probabilities:
    - Yang et al. 2016 accepted a swap when the new bond energy was lower, which breaks detailed balance.
    - Wick et al. 2025 chose pairs with a distance bias and included the reverse-move probability, which keeps it.
- **Three-body swap potential.** This is Sciortino's (2017) potential, v3 = lambda eps v_ij v_ik. It goes into the CG
  evaluator as a pair style, with the three-body virial term that Ciarella et al. 2018 give. It is for star and
  telechelic models, where the barrier is eps(lambda - 1). It gives continuous dynamics without MC moves and is the
  standard coarse-grained benchmark.
- **Counters.** The thermo rows record attempts, accepted swaps, and broken and formed bonds per type.

**Charges.** An all-atom exchange restores the same atom types: ester stays ester and S-S stays S-S. So `charges keep` is
correct for transesterification and disulfide metathesis.

**Force-field terms.** Class II force fields (PCFF, COMPASS) need their cross terms regenerated for every touched angle
and torsion. CVFF, OPLS-AA and GAFF need only the regular terms. `--complete ml` fills terms missing from the force field.

### 3.3 LAMMPS export

- **Bead-spring.** The KG writer gains `fix bond/swap` decks. The manual page (read 2026-10-22) sets the requirements:
  molecule IDs, `special_bonds lj 0 1 1`, and no dihedral or improper potential in the simulation. The fix reports
  cumulative swaps performed and attempted. Bond/swap is therefore for bead-spring models only. The CAPS engine and the
  LAMMPS deck can then be checked against each other on the same system.
- **All-atom.** `write_bond_react` gains exchange templates:
  - pre and post templates with the same atom count and the same bond count;
  - the reverse direction generated from the same template;
  - `prob` from the library, or the Arrhenius constraint once it is verified;
  - a `rate_limit` option.
- **CLI.** A new `caps exchange export` entry, since bond/react export can be reached today only from the C API, Python
  and the Studio.

### 3.4 Protocols

Each protocol is a recipe stage or a `caps exchange` sub-command, so it runs as a cluster job with no extra code.

| Protocol | What it does | From |
|---|---|---|
| `exchange relax` | Equilibrium run with exchange; writes the bond-partner autocorrelation and the stress autocorrelation | Ciarella 2018; Lin 2026; Branham-Ferrari 2026; Mamat 2026 |
| `exchange step-strain` | Step strain (uniaxial or shear), then a hold with exchange on, recording sigma(t) | Yang 2016 (30 % at 1e9 s^-1, then hold); Zhao 2022 (300 % uniaxial) |
| `exchange tv` | Heating or cooling scan with exchange on; volume, CTE and swap rate at each T | Perego and Khabaz 2020 (via the Karatrantos review); Singh 2020; Mamat 2026 |
| `exchange creep` | Constant stress with exchange; isothermal, or heating at a set rate (non-isothermal) | Singh 2024 (NPT 5 ps / NVE 4 ps loop); Zheng 2025 (100 to 800 K at 4 K/ps, 50 MPa); Mamat 2026 (J(t)) |
| `exchange heal` | Cut, wall or gap; optional waiting time; wall removed; run with exchange; then a tensile or triaxial test | Stukalin 2013; Ciarella 2019; Chen 2023; Uyama 2026; Yuan 2026 |
| `exchange weld` | Two films with walls and an initial gap of 0.5 to 2 A; the wall is removed after the first cross-interface event | Yang 2016 |
| `exchange recycle` | Fragment the network, reassemble it by melting the fragments together, run exchange cycles, report recycling efficiency | Wick 2025 |

Four builders are needed:
- **Cut:** removes all bonds crossing a plane and inserts a repulsive wall.
- **Gap:** two slabs from `build_interface`, held apart by walls.
- **Crack or void:** removes atoms or beads in a region.
- **Fragments:** breaks the network into pieces and packs them.

### 3.5 Analysis

- **Bond-partner autocorrelation.** The fraction of exchangeable sites still bonded to their time-zero partner. Two fits:
  - KWW, or 1/e time (Mamat 2026);
  - the two-step form f(t) = 1/2 [1 + exp(-t/tau)] exp[-(t/tau_diff)^beta] (Branham-Ferrari 2026, Eq. 3).
- **Exchange time and Arrhenius analysis.** tau_ex comes from counted events. tau_ex, tau_bond and tau_alpha each get an
  Arrhenius fit with Ea and prefactor. A test checks for the crossover between Arrhenius and super-Arrhenius behaviour
  (Branham-Ferrari 2026). Mamat 2026 adds the linear law tau_bond = m(T) t_exch + b(T).
- **G(t).** Computed from the stress autocorrelation with a multi-tau correlator, over all three shear components, with
  the three-body virial where it applies. Then:
  - G' and G'' by piecewise-linear Fourier integration (as in Chen 2023);
  - tau_rel, from G(t) reaching a set fraction (Ciarella 2018 uses 1e-4 of G(0));
  - optional KWW, single-exponential and Maxwell fits.
- **Tv.** Three routes:
  - a hinge fit of the CTE against T (Perego and Khabaz 2020);
  - the temperature where tau_bond crosses tau_alpha;
  - shift factors a_T from time-temperature superposition (Danikas 2026; Mamat 2026 shows TTS fails once exchange is on,
    so this must report the failure).

  For comparison with experiment, the Arrhenius extrapolation of tau to the 1e12 Pa s viscosity criterion (Arbe 2023).
- **Healing.** Three kinds of output:
  - **Interface counts:** bridges across the cut plane N_br(t) and the bond profile B(z) (Ciarella 2019); bonds across
    the interface, normalised by the bulk value; penetration depth (Yang 2016).
  - **Repair by pathway:** the fraction of broken bonds reformed by each pathway (Uyama 2026).
  - **Efficiencies:** four definitions, each always reported with its definition:
    - toughness ratio, the area under the stress-strain curve (Chen 2023);
    - modulus ratio (Yuan 2026);
    - peak-stress ratio;
    - recovered stress at fixed strain (the recycling efficiency of Wick 2025).
- **Topology.** Four measures:
  - loop fraction (primary loops speed relaxation about tenfold in Ciarella 2018);
  - largest cluster;
  - gel point by reduced molecular weight (Rispo Constantinou 2026);
  - shortest-path links (Wick 2025).

All of this goes into core, so the CLI, Python and Studio share it. The Studio's Arrhenius code in `CurveFit.cs` should
then call the core fit.

## 4. Default protocol per chemistry

| Chemistry | Model level | Default protocol | Source |
|---|---|---|---|
| Generic associative | Kremer-Grest, FENE K = 30, R0 = 1.5; end beads exchange; swap within 1.5 sigma; Metropolis; about one attempt per step | `exchange relax` over T, then `exchange tv` | Branham-Ferrari et al. 2026 |
| Generic associative, network with defects | Kremer-Grest, LJ cutoff 2.5 sigma, FENE 30/1.5, dt 0.005; MC swaps every t_exch = 100 to 1000 tau | `exchange relax`, `exchange creep` (J(t) at sigma0 = 0.1) | Mamat et al. 2026 |
| Generic associative, star networks | 8-arm stars, three-body swap, 20-10 generalised LJ with eps = 100 kT and sigma = 0.5 nm | `exchange relax`, `exchange heal` | Ciarella et al. 2018, 2019 |
| Transesterification epoxy (DGEBA + acid) | All-atom PCFF; cured by react to about 92 to 97 %; template exchange with an Arrhenius factor | `exchange weld` at 450 K and 420 K, `exchange step-strain` | Yang 2016; barrier per 3.1 |
| Transesterification methacrylate (Bis-GMA/HEMA type) | All-atom OPLS-AA, REACTER polymerisation then ester-hydroxyl exchange | `exchange tv`, `exchange heal` | Yuan et al. 2026 |
| Disulfide aromatic (DGEBA + aminophenyl disulfide) | All-atom CVFF or PCFF; S...S capture 4.12 A | `exchange creep` (NPT 5 ps / NVE 4 ps loop), `exchange heal` | Singh 2020, 2024; Uyama 2026 |
| Transesterification shape-memory methacrylate | All-atom OPLS-AA; pairs within 4.25 A, chosen with weight exp(-0.2 R^2) | `exchange recycle` | Wick et al. 2025 |
| Rubber with hydrogen-bond side groups (NR) | All-atom COMPASS, physical bonds only | `exchange heal` with a wall gap of about 3 A, then a triaxial test | Chen et al. 2023 |
| Imine polyisoprene (NR-like vitrimer) | All-atom OPLS-AA or GAFF network; CG for the cluster structure | Tg and S(q) first, exchange later | Arbe et al. 2023 (experimental targets) |
| Boronic ester acrylate and methacrylate elastomers | All-atom, with boron parameters from `ff complete` | Tg check only at first | Teotonico et al. 2024; Danikas et al. 2026 |

## 5. Validation targets, in the order to reproduce them

1. **Rao, Xia and Ni 2024 monovalent test.** A unit test with no polymer: 400 pivots, box 20 r_b, dG = -1 eps. Species
   fractions are compared against their Eq. 5. This checks the acceptance rule and detailed balance.
2. **Wilcoxson et al. 2026 frozen four-reaction test.** A unit test: bond populations must fall within 0.02 of the exact
   partition function. This checks bind/unbind.
3. **Branham-Ferrari et al. 2026 barrier.** KG telechelic strands with 3000 exchangeable bonds. The Arrhenius fit of
   tau_ex should give E* = 3.24 eps for sigma_LJ = 1, and E* should fall when sigma_LJ is reduced to 0.985, 0.96 and 0.92.
   Run it in the CAPS engine and with the exported `fix bond/swap` deck.
4. **Mamat et al. 2026.** A KG network with MC swaps every t_exch = 100, 200, 500 and 1000. Targets:
   - Tg about 0.43 from V(T);
   - apparent activation E_app = 1.225, 1.134, 0.908 and 0.995 (reduced units) from Arrhenius fits of tau_bond;
   - TTS failing for every finite t_exch.

   The total bead count is not stated, and the agent estimated 15,840 [unverified]. A smaller cell may be needed.
5. **Rispo Constantinou et al. 2026 gel point.** 2000 chains of 20 beads, 5 stickers per chain, phi = 0.5. Target gel point
   p = 0.4 to 0.475. Their MC-exchange package (on PyPI) allows a direct comparison of swap statistics. 40,000 beads is
   above the 15,000 cap, so run a quarter-size cell and say so.
6. **Ciarella et al. 2018.** Defect-allowing networks should relax about ten times faster than defect-free ones, with the
   fluid-solid change near a 10 kT barrier. This needs the three-body potential. Their 50,000 beads must also be reduced.
7. **Ciarella and Ellenbroek 2019 healing times.** Mechanical healing at 0.15 tau_sw, interdiffusion after about 5 tau_sw,
   bulk relaxation at 13.2 tau_sw.
8. **Yang et al. 2016 weld (atomistic).** Two films, 10,860 atoms (under the limit), at 450 K. Targets:
   - tau_BER about 1.26 ns;
   - stress-relaxation tau of 4.2 ns from a single exponential;
   - Ea of 67.67 kJ/mol recovered from the 420 K and 450 K runs;
   - modulus back to 600 MPa after about 11 ns.
9. **Singh et al. 2020 and 2024 disulfide epoxy.** Targets:
   - density 1.159 g/cm3 at 300 K after 95 % cure;
   - CTE 41.1 uK^-1 below Tv;
   - E 4.5 GPa at 2e8 s^-1;
   - creep with p = 0 against p = 1 at 600 K and 500 MPa.

   Their 68,608-atom cell is over the 15,000-atom limit, so use a 4x4x4 replica (8,576 atoms) and say so.
10. **Chen et al. 2023 NR healing.** COMPASS. Healing efficiency (toughness ratio) of 44.7 % after 1 ns and 82.3 % after
    30 ns at 498 K. The atom count is not stated, so check that it fits the limit.
11. **Experimental anchors.** These cannot be reproduced in MD time; they serve for trends and Tg only:
    - imine PI: Tg 216.5 to 221 K, Tv about 260 to 274 K, flow Ea 66 and 79 kJ/mol (Arbe 2023);
    - dioxaborolane PBMA: Tg 293 to 303 K, Tv/Tg about 1.4 to 1.5, Ea 80 to 94 kJ/mol (Danikas 2026);
    - DGEBA/2-AFD: Tg about 72 C, and 94 % scratch healing at 90 C after 1 h (Golbabaei 2026);
    - DGEBA/DTBA: stress-relaxation Ea 48 kJ/mol as printed, 40 to 46 kJ/mol from its own data (Sharma 2026);
    - DGEBA with malic acid: Tg 348 K; with naphthalene diacid: Tg 395 K (Zheng 2025, Digital Discovery).

Local runs stay at 5,000 steps or fewer. Longer runs go to the cluster as jobs, at most 80 cores each.

## 6. Risks

- **Time scales.** Real exchange takes seconds to hours, but atomistic MD covers ns. Every atomistic paper read here
  speeds exchange up: a logistic p(T), boosted ReaxFF, or one swap per iteration. CAPS must report results in units of
  the exchange time (t/tau_BER, t/tau_sw) and state the speed-up. Absolute times cannot be compared with experiment.
- **Barrier disagreement.**
  - For aromatic disulfides, Irigoyen 2022 (dG about 50 kJ/mol) and Uyama 2026 (150 to 180 kJ/mol) differ threefold;
    they use different pathways and references.
  - Sharma 2026's printed Ea does not match its own fit.

  The library must keep every value with its method and let the user choose. CAPS's own QM path should recompute every
  barrier it uses.
- **Detailed balance.** Two published rules break it:
  - energy-lowering acceptance (Yang 2016);
  - the logistic p(T) of Singh 2020 and Yuan 2026.

  Multivalent swaps also need the Rao bias to keep it. CAPS defaults should be balanced moves.
- **Class II force fields.** Swapping bonds in PCFF or COMPASS changes the bond-bond and bond-angle cross terms. Uyama 2026
  does not say how they handled this. Regenerate the cross terms of all touched angles and torsions, and test energy
  conservation across a swap.
- **No catalyst in any atomistic study** (noted by Karatrantos 2024). Catalyst appears only as a rate factor.
- **System size.** Several reference systems exceed the 15,000-atom limit:
  - Singh: 68,608 atoms;
  - Wick: 25,039 and 37,376 atoms;
  - Hou: 50,267 atoms;
  - Yuan: about 33,000 atoms, the agent's estimate;
  - the Rispo Constantinou and Ciarella CG systems: 40,000 to 50,000 beads.

  Use reduced cells and report the difference.
- **LAMMPS limits.** `fix bond/swap` cannot be used with dihedrals or impropers, so it is for bead-spring models only. The
  bond/react Arrhenius constraint syntax still needs checking against the full manual page.
- **ReaxFF.** The accelerated ReaxFF route (Zheng 2025, Kamble 2022, Vermeersch 2024) needs a ReaxFF parameter set from
  the original papers. CAPS has no ReaxFF engine, so ReaxFF stays an export target at most.
- **Rubber chemistries are thin in the simulation literature.** No paper read here simulates exchange in ENR, SBR or a
  sulfur-cured diene rubber. Chen 2023 (NR with hydrogen bonds) is the only rubber simulation with healing. The rubber
  entries will start with experimental anchors only, and CAPS will have to compute their barriers. The anchors are
  Imbernon et al. 2015 (ENR with a disulfide diacid; paywalled), Arbe 2023 (PI imine) and Danikas 2026 (PBMA boronic).

## 7. Staged build plan

| Stage | Content | Done when |
|---|---|---|
| E0 Analysis | Bond-partner autocorrelation with KWW and two-step fits; core Arrhenius fit; multi-tau stress correlator, G(t), G'/G''; Tv by CTE hinge and by a_T; healing metrics (bridges, B(z), four efficiencies) | Unit tests on synthetic data pass, and the Studio's `CurveFit` uses the core fit |
| E1 Mutable topology | `Topology` swap, break and make, with exclusion and neighbour updates; an energy continuity test across a swap | NVE drift in a KG melt is unchanged by swaps |
| E2 CG exchange | Pair swap (Metropolis plus Arrhenius), multivalent swap, bind/unbind; `fix bond/swap` export | Targets 1 to 5 reproduced, and the CAPS engine and LAMMPS agree on tau_ex |
| E3 Three-body potential | Sciortino potential and its virial in the CG evaluator | Targets 6 and 7 |
| E4 Exchange library | `data/exchange/library.json` with the entries of 3.1; new schemes (disulfide metathesis, transimination, boronic ester, VU, siloxane equilibration, retro DA) | Every number has a DOI, and the schemes round-trip through `template_from_scheme` |
| E5 All-atom exchange | Template exchange move in the CAPS engine; bond/react exchange export with probability or Arrhenius constraint; `caps exchange export` | Target 8 in CAPS and in LAMMPS; target 9 on the reduced cell |
| E6 Healing and welding | Cut, gap, crack and fragment builders; `exchange heal`, `weld`, `recycle`, `creep`, `step-strain` and `tv` protocols, plus recipe stages | Targets 7, 8 and 10 |
| E7 Rubber vitrimers | ENR with a diacid or disulfide diacid, disulfide-cured NR and SBR, imine PI; barriers from CAPS QM jobs where papers have none | Tg and trend checks against experimental data scaled by Arrhenius |
| E8 Studio | Exchange page: chemistry picker showing barrier and source, protocol forms, live swap counter, G(t) and healing plots | The self-test covers one CG case and one all-atom case |

E0 should come first. It is useful on its own, since G(t), Arrhenius and Tv apply to any network CAPS builds today, and it
carries no risk. E1 is the one change to the engine's design that everything else depends on.

## 8. Rubber composites to study

These are the systems the module is for, in order. Each stays under 15,000 atoms all-atom; bead-spring models take the
long-time behaviour.

1. **ENR crosslinked with a diacid (transesterification).** The beta-hydroxy esters formed when an acid opens the ENR
   epoxide exchange with free OH. CAPS already makes this link (`enr_acid_ester`, `epoxy_acid`). With a disulfide diacid
   (the chemistry of Imbernon et al. 2015, paywalled) the network carries two exchange chemistries. First targets: Tg,
   G(t) and the Arrhenius energy of stress relaxation, against experimental trends.
2. **ENR vitrimer on a fibre surface.** A fibre surface carrying OH or ester groups (cellulose or natural fibre, PET cord)
   can bond into the ENR network by the same exchange. Interface welding (bridges, B(z)) and pull-out work before and after
   healing with `caps pull`. The user's composites use cellulose / natural fibres (decided 2026-10-24), so the
   surface is cellulose: its hydroxyls take part in transesterification with the ENR network's beta-hydroxy esters.
3. **Sulfur-cured NR and SBR (polysulfide exchange).** S-S bridges swap when hot. Stress relaxation and crack healing
   compared for mono-, di- and polysulfide crosslinks, using the existing sulfur cure and sulfur census. The exchange
   barrier of aliphatic polysulfides is not in any paper read here; it needs a CAPS QM calculation on model compounds.
4. **ENR/PBS through MAH.** The ester links made through MAH can swap; whether exchange makes the blend more compatible
   (phase stability as in Rispo Constantinou et al. 2026) and how strong the ENR/PBS interface becomes.

Imine polyisoprene (Arbe et al. 2023) and boronic-ester elastomers (Danikas et al. 2026) serve as experimental checks only.

## 9. Status

**E0 (analysis) is built** (2026-10-24), in `core/include/caps/exchange_analysis.hpp` and `core/src/exchange_analysis.cpp`:

- `fit_arrhenius`, `crossing_temperature`: Ea with standard errors (kJ/mol or reduced units), and the temperature where
  two relaxation times cross.
- `MultiTauCorrelator`, `stress_relaxation`, `dynamic_moduli`, `fit_kww`: G(t) from pressure-tensor samples, G'(ω) and
  G''(ω) by the exact piecewise-linear transform, the Green–Kubo integral to G's first zero, KWW fits.
- `bond_memory`, `fit_two_step`: partner autocorrelation, its 1/e time, KWW and the two-step form of Branham-Ferrari et
  al.; partner changes between frames.
- `interface_bonds`, `sides_from_positions`: bridges across a cut and B(z), periodic.
- `healing_efficiency`: toughness, peak-stress, modulus and stress-at-strain ratios, each named.

It is reachable as `caps exchange arrhenius | modulus | memory | bridges | heal`, as the C API entry
`caps_exchange_analyze` (ABI 74) and as Python `caps.exchange(op, ...)`. Unit tests: `Exchange.*` in
`tests/test_exchange.cpp` (10 tests on exact or analytic answers). The theory pages have four new entries
(`arrhenius-exchange`, `stress-relaxation-gk`, `bond-memory`, `healing-efficiency`).

Not in E0, by design: Tv from the thermal-expansion hinge uses the existing `fit_bilinear` (`caps tg --fit` on a volume
scan with exchange on). The Studio page waits for E8. The partner tables that `caps exchange memory` reads have to be
written from LAMMPS output by the user until CAPS's own engine exchanges bonds (E1, E2).

**E1 (mutable topology) is built** (2026-10-24), in `core/include/caps/topology_edit.hpp`, `core/src/topology_edit.cpp`,
`Evaluator::reload()` and `DynamicsOptions::topology`:

- `TopologyEditor` breaks and makes any bonds in one change (a swap of partners, or bonds broken or made alone) and
  rebuilds the force field around them from the force field's own terms by type: bonds (harmonic, FENE and the other
  extended forms, class II quartic), angles with Urey–Bradley terms, proper torsions, impropers (Fourier, harmonic,
  inversion), 1-4 pairs with their own scales, and exclusions by the rule the force field used (1-2 only for
  `special_bonds fene`, 1-2/1-3/1-4 for all-atom force fields). A change the tables cannot describe throws and leaves the
  structure and the force field as they were.
- `Evaluator::reload()` takes a force field changed in place.
- `DynamicsOptions::topology` calls a callback every `topology_every` steps with the run's own structure and force
  field. The callback may edit them, or put in a force field assigned again to the changed structure (the route for
  chemistry the tables cannot describe: new atom types, class II cross terms). `ctx.potential()` gives the energy for a
  Metropolis test. The energy the changes add is kept out of the conserved quantity.
- Tests (`Topology.*`, 6):
  - A swap of groups between two butanes and between two toluenes gives the same terms, exclusions, 1-4 pairs and
    energies, to 10⁻⁸, as GAFF assigned to the new structure from scratch. The same holds for Kremer–Grest chains with
    bending.
  - An impossible change leaves everything untouched.
  - Breaking a bond and making it again restores the energy.
  - NVE with 29 accepted Metropolis swaps (58 kcal/mol moved): the conserved quantity jitters 0.52 kcal/mol from row to
    row against 0.48 without swaps.
  - Rejected swaps restore the energy to 10⁻⁹.
  - The callback can reassign the force field.
  - Constraints are refused.

**E2 (exchange moves) is built** (2026-10-24), in `core/include/caps/exchange_moves.hpp`, `core/src/exchange_moves.cpp`
and `write_kg_bond_swap` in `core/src/kremer_grest.cpp`:

- `ExchangeEngine` makes Metropolis–Hastings moves of three kinds. Each acceptance carries the proposal ratio, so
  detailed balance holds:
  - swap: two bonds trade partners, as LAMMPS fix bond/swap, including its molecule-ID rule (a1/b1 and a2/b2 share
    labels; checked in the LAMMPS source);
  - shift: a pivot moves a bond between residues, with residue valences (the multivalent case);
  - bind/unbind: a binding constant K = exp(ε/k_B T).

  A barrier E_a slows every move without biasing it.
- ΔU comes from the atoms within three bonds of the change, before and after. It equals the whole system's ΔU, PME
  included (checked to 10⁻⁷). Force fields with many-body terms, hydrogen bonds or virtual sites fall back to the whole
  system's energy.
- `exchange_callback` plugs the engine into `DynamicsOptions::topology`. `caps exchange run` runs it from the command
  line.
- `caps exchange bond-swap` writes the LAMMPS deck for a Kremer–Grest system in either of LAMMPS's two uses:
  - chains: molecule ID = place along the chain;
  - stickers: one molecule ID, `newton on off`, no angles across sticker bonds.
- Tests (`ExchangeMoves.*`, 8):
  - frozen-atom sampling reproduces the Boltzmann weights of every topology (swap, shift at valence 1 and 2);
  - bind reproduces the four-site partition function Z = 1 + 4K + 2K² of Wilcoxson et al.;
  - local ΔU equals the whole system's (DSF, PME, bead-spring);
  - labels keep chain lengths, and without them lengths change;
  - swaps run during MD;
  - the LAMMPS deck has the right content.

**E3 (three-body swap potential) is built** (2026-10-24), as `SwapModel` (`ForceField::swap3`) and
`swap_model_terms` in `core/src/field.cpp`:

- Sciortino's potential as Ciarella et al. 2018 write it, with a shift option:
  - a generalised Lennard-Jones pair (20-10 by default) between complementary site classes;
  - the three-body term λ ε v̂ v̂ on every site with two partners.
- Sites are found by their own cell list. Bonded pairs are excluded. Each distance gets its virial term, so stress
  and G(t) include the triplets.
- Tests (`SwapModel.*`, 4):
  - the barrier ε(λ − 1) comes out exactly;
  - forces match finite differences, and the virial matches stretch and shear strain derivatives;
  - the cell-list sums equal a direct reference over 200 sites;
  - exclusions are honoured and a too-small cell is refused;
  - NVE with pairs forming and swapping conserves energy.

Not in E3: a command-line or Studio way to set up the model (sites and pairs are set on the force field from C++ for
now), and a builder for star or telechelic networks; both belong with E6 and E8.

Next: E4, the exchange chemistry library.

## 10. Chemistries the module covers

The user asked (2026-10-24) for every vitrimer and self-healing chemistry, not only the rubber systems of section 8.
Every family below maps onto the same three operations of the E1 editor, so the engine is common to all of them and the
chemistry lives in the exchange library (E4):

| Family | Mechanism | Change in the editor | Example systems | Barrier source (section 3.1 or to find) |
|---|---|---|---|---|
| Transesterification (beta-hydroxy ester, with Zn, TBD or self-catalysed) | associative | swap (alkoxy groups between esters, through a free OH) | epoxy-acid and epoxy-anhydride vitrimers, ENR + diacid, polyester vitrimers | Bhusal 2021 (14.6 kcal/mol); Yang 2016 (68.18 kJ/mol) |
| Transcarbamoylation (urethane exchange) | associative or dissociative | swap, or break then make | polyurethane vitrimers | to find (Li 2023 review quotes Fortman's 35.5 to 70.7 kJ/mol) |
| Transamination of vinylogous urethanes and ureas | associative | swap | VU elastomers, PDMS-VU | Hamzehlou 2022 |
| Transimination and imine metathesis | associative | swap | imine polyisoprene (NR-like), imine PS and PPG | experimental Ea (Arbe 2023); ML barrier Sun 2023 (paywalled) |
| Disulfide and polysulfide exchange (metathesis, radical 2+1, thiol-disulfide) | associative (radical step dissociative) | swap, or break then make | aromatic disulfide epoxies, sulfur-cured NR and SBR, disulfide PU | Irigoyen 2022; aliphatic barrier to compute |
| Diselenide exchange | associative (radical) | swap | diselenide PU | to find |
| Boronic ester and dioxaborolane metathesis | associative | swap | PBA, PBMA, PB, PDMS vitrimers | Teotonico 2024; Danikas 2026 (experimental) |
| Siloxane equilibration (silanolate) | associative | swap | PDMS vitrimers, silicone rubber | to find |
| Transalkylation (sulfonium, triazolium) | associative | swap | poly(thioether) vitrimers | Garcia Villatoro 2026 (experimental 123 kJ/mol) |
| Acetal and N,O-acetal exchange | associative | swap | bio-based CANs | to find |
| Diels–Alder (furan–maleimide) | dissociative | break, then make | DA self-healing networks | Vermeersch 2024 ReaxFF (paywalled) |
| Hindered urea and thiol-Michael | dissociative | break, then make | reversible PU and acrylates | to find |
| Hydrogen bonds, metal–ligand, ionic clusters | physical | none (force field), or bind/unbind stickers | Chen 2023 NR H-bond network, ionomers | force field only |

The rows marked "to find" need a paper read, or a CAPS QM calculation, before they go into the library. No value is
copied from another program's parameter files.
