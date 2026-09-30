// CAPS solvation (design/boards/SolvationBuilder): a solute held at the centre of a periodic box, solvent molecules and
// ions packed around it with CAPS Pack.
//
//  Box        cubic (edge), rectangular (a, b, c) or the solute's extent plus a padding on every side.
//  Solvent    water (TIP3P, SPC/E or TIP4P/2005 geometry and charges; the TIP4P M-site charge sits on O, as LAMMPS'
//             tip4p pair styles expect) or an organic solvent built from SMILES (toluene, cyclohexane … for swelling
//             rubber). The count fills the free volume at the solvent's density: the box less the grid points inside
//             the solute's van der Waals spheres (Bondi radii).
//  Ions       none, neutralise the solute, a salt concentration (pairs = c · N_A · V_free, plus the neutralising ions),
//             or custom counts; each ion in water takes the place of one water molecule.
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "caps/pack.hpp"
#include "caps/system.hpp"

namespace caps {

struct SolventInfo {
  std::string id, name, smiles;   // water's smiles is "O"
  double density = 1.0;           // g/cm³ at 25 °C
  std::string use;                // what it is for
};
const std::vector<SolventInfo>& solvent_library();

struct SaltInfo {
  std::string id, cation, anion;  // "NaCl", "Na", "Cl"
  int zc = 1, za = -1;            // charges
};
const std::vector<SaltInfo>& salt_library();

struct SolvateOptions {
  int shape = 0;                  // 0 cubic, 1 rectangular, 2 solute plus padding
  double edge = 30;               // Å, cubic
  Vec3 edges{30, 30, 30};         // Å, rectangular
  double padding = 10;            // Å on each side, shape 2
  double tolerance = 2.0;         // Å between atoms of different molecules
  std::string solvent = "water";  // an id of solvent_library()
  std::string water_model = "TIP4P/2005";   // any of water_models() by name or id (SPC, SPC/E, SPC/Fw, TIP3P …, TIP4P/2005, OPC)
  double density = 0;             // g/cm³; 0: the solvent's
  int molecules = 0;              // solvent molecules; 0: from the density over the free volume
  int ion_mode = 2;               // 0 none, 1 neutralise, 2 concentration (and neutralise), 3 custom counts
  std::string salt = "NaCl";
  double concentration = 0.15;    // mol/L
  int cations = 0, anions = 0;    // custom counts
  uint64_t seed = 1;
  int threads = 0;
  std::function<bool(const PackProgress&)> progress;
};

struct SolvatePlan {
  Vec3 box{0, 0, 0};              // Å
  double box_volume = 0, solute_volume = 0, free_volume = 0;   // Å³
  int solute_atoms = 0;
  double solute_charge = 0;       // e, the sum of the solute's charges
  int solvent = 0, cations = 0, anions = 0;
  std::string solvent_name, cation, anion;
  double solvent_mass = 0;        // g/mol of one solvent molecule
  double concentration = 0;       // mol/L of the salt actually placed (pairs over the free volume)
  double density = 0;             // g/cm³ of the whole box when packed
  std::vector<std::string> notes;
};

// The box and the counts, without packing (milliseconds). solute may be null (pure solvent).
SolvatePlan solvate_plan(const System* solute, const SolvateOptions& o);

struct SolvateReport {
  SolvatePlan plan;
  PackReport pack;
  std::vector<std::string> notes;
};
// The solvated box: the solute (molecule 1, centred, unwrapped), then the solvent, then cations and anions. Throws
// PackError when the tolerance cannot be met.
System solvate(const System* solute, const SolvateOptions& o, SolvateReport* report = nullptr);

// One molecule of the solvent (water in the chosen model), with types, bonds and charges.
System solvent_molecule(const SolvateOptions& o, std::string* name = nullptr);

}  // namespace caps
