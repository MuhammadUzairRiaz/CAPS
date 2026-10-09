#pragma once
// Analysis of finished VASP runs of 2D slabs and adsorption sets (the DFT surface & adsorption workbench): binding
// energies with a settings check, adsorption geometry, work function from LOCPOT, grouped DOS, charge-density
// difference, Bader transfer, nitrile stretch shift, AIMD stability, and a summary across surfaces — each a Json
// report (the CLI prints it, the Studio draws it) and files beside the run (csv / json).
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

// E_bind = E0(complex) − E0(slab) − E0(molecule) from 04_static (03_relax flagged "not final"), eV and kJ/mol, ranked,
// the spread of the best three, and every run's ENCUT/IVDW/ISMEAR/SIGMA/PREC/IDIPOL/LREAL/GGA compared with the slab's.
// Writes binding_energies.csv in the set.
Json binding_energies(const std::string& set_root);

// For each complex: lowest heavy atom and COM above the outer termination plane, the anchor (nitrile N, else the lowest
// heavy atom) to its nearest surface atom, O–H···N hydrogen bonds (H···N < 2.5 Å, angle > 120°), C–H···O/F contacts
// (< 2.7 Å), contacts (< 3.5 Å), tilt of the molecule's long axis. Writes geometry.json per complex, geometry_table.csv.
Json adsorption_geometry(const std::vector<std::string>& complex_dirs);
Json adsorption_geometry_of(const System& at, const System& slab_poscar);

// Work function from LOCPOT (planar average along z): vacuum plateaus > 4 Å from the atoms on both sides,
// φ = V_vac − E_F per face, plateau flatness; with `reference` (the slab alone) Δφ = φ_top − φ_top(reference).
Json work_function(const std::string& stage_dir, const std::string& reference_dir = "");

// DOS grouped: total, the surface metal, the outer termination layer, and (for a complex) the molecule and its anchor
// atom; DOS at E_F, HOMO-/LUMO-like molecular peaks. Writes dos.json and dos_groups.csv.
Json dos_analysis(const std::string& stage_dir, const std::string& slab_poscar = "");

// Δρ = ρ(complex) − ρ(slab part) − ρ(molecule part) at the same geometry (grids must agree), CHGCAR_diff, the
// plane-averaged profile and cumulative ΔQ(z), ΔQ at the mid-plane between surface and molecule (> 0: the slab side
// gained electrons). Reads <cx>/charge/04_static, <cx>/cdd_slab/04_static, <cx>/cdd_mol/04_static.
Json charge_density_difference(const std::string& complex_dir, const std::string& slab_poscar);

// Bader transfer: electrons gained by the molecule and by the slab (ZVAL from the POTCAR), the total (~0).
Json bader_transfer(const std::string& slab_poscar, const std::string& complex_poscar, const std::string& acf, const std::string& potcar = "");

// ν(C≡N) of the adsorbed and free molecule (2000–2500 cm⁻¹ window) and the shift, from <cx>/freq and <set>/molecule/freq.
Json nitrile_shift(const std::string& complex_dir);

// AIMD stability from seg_*: T(t) with its expected fluctuation T√(2/3N), energy drift (meV/atom/ps), COM and anchor
// height, closest heavy contact, contacts, broken bonds; production after `equil_ps`; verdict STAYS ADSORBED /
// DESORBS · UNSTABLE. Writes md_timeseries.csv and md_summary.txt.
Json md_analysis(const std::string& aimd_dir, double equil_ps = 1.0, double dt_fs = 1.0);

// Collects every surface's numbers (lattice constant, bonds, E_bind, best orientation, Bader, Δφ, DOS(E_F), CDD ΔQ,
// ν shift, AIMD verdict) into results_summary.csv / .txt: the hand-over to classical MD. `main` holds structures/<case>,
// adsorption/<case>_<molecule>, aimd/<case>_best.
Json summarize_results(const std::string& main);

}  // namespace caps
