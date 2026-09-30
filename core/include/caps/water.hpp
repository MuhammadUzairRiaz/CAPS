// Water models: the SPC and TIP families (and OPC), each with the geometry, charges and Lennard-Jones parameters of its
// publication, as LAMMPS's own tables give them (docs.lammps.org Howto_spc, Howto_tip3p, Howto_tip4p) where they do.
//
// Three-site models carry their charges on O and H. Four-site models (TIP4P, TIP4P-Ew, TIP4P/2005, TIP4P/Ice, OPC) put
// the negative charge on a massless site M on the H–O–H bisector, d_OM from O: CAPS adds M as an atom of its own after
// each water's hydrogens (element 0, name MW), placed at every evaluation as (1 − α) O + α/2 (H1 + H2) with
// α = d_OM / (r_OH cos(θ/2)), its force handed back to O and H in the same proportions (the linear virtual site of
// GROMACS's TIP4P, and where LAMMPS's tip4p pair styles put it). The LAMMPS files leave M out and give O its charge, with
// pair_style lj/cut/tip4p/long and the model's d_OM; the GROMACS files keep it as a virtual site.
//
// Rigid models are run with their bonds and angle constrained (SHAKE / RATTLE, as LAMMPS fix shake and GROMACS settles):
// their bond and angle constants here serve a flexible run only. SPC/Fw is flexible by design.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct WaterModel {
  std::string id, name, citation, bibkey;
  int sites = 3;                // 3 or 4
  double r_oh = 0, theta = 0;   // Å, degrees
  double q_h = 0;               // e
  double q_neg = 0;             // e: O's charge (3 sites) or M's (4 sites; O then carries none)
  double d_om = 0;              // Å, 4 sites
  double eps_o = 0, sigma_o = 0;   // O–O Lennard-Jones: kcal/mol, Å
  double eps_h = 0, sigma_h = 0;   // hydrogen Lennard-Jones (CHARMM's TIP3P), else 0
  bool rigid = true;
  double k_bond = 0, k_angle = 0;  // K (r − r0)², K (θ − θ0)² (LAMMPS harmonic): kcal/mol/Å², kcal/mol/rad²
  std::string k_source;         // where the bond and angle constants come from
  std::string note;             // what to know (a value that differs between sources …)
};

const std::vector<WaterModel>& water_models();
// By id ("tip4p2005") or name ("TIP4P/2005"), any case; throws std::invalid_argument naming the known ones.
const WaterModel& water_model(const std::string& id_or_name);
// M = (1 − α) O + α/2 (H1 + H2)
double water_m_alpha(const WaterModel& m);

// The waters of s: molecules of one O bonded to two H (and an M site, element 0), as {O, H1, H2, M or −1}.
std::vector<std::array<int64_t, 4>> find_waters(const System& s);

// Gives every water of s the model's geometry (O and the molecule's plane and orientation kept, the O–H lengths and H–O–H
// angle set) and charges; a 4-site model adds each water's M after its hydrogens, a 3-site one removes an M there was.
// Atoms after an insertion move up (bonds follow). Returns the number of waters; notes say what was done.
size_t apply_water_model(System& s, const WaterModel& m, std::vector<std::string>* notes = nullptr);

// The model's force field for the atoms `atoms` of s (its waters, after apply_water_model), in their order: types OW, HW
// (and MW), charges, masses (M none), the O–H bonds and H–O–H angle, M as a virtual site, every pair inside a water
// excluded. Throws FieldError when an atom is not part of a water.
ForceField water_forcefield(const System& s, const WaterModel& m, const std::vector<uint32_t>& atoms);

}  // namespace caps
