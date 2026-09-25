// CAPS recipes (design/boards/CommandLine): a structure made end to end from one YAML or JSON file — build, type, grow,
// relax, dynamics, equilibrate, analyse, export — each stage reported as it runs, the result written with its provenance
// manifest. The same stages as the Studio's builders and runs.
//
//   recipe: 1
//   name: ps_cell
//   build:  { polymer: { smiles: "*CC(*)c1ccccc1", dp: 40, chains: 20, tacticity: atactic } }
//           (or molecule: { smiles: CCO } · file: cell.data)
//   type:   { forcefield: gaff2 | uff | default | PATH.json, charges: types | gasteiger | qeq | keep }
//   grow:   { density: 0.5, contact_scale: auto, seed: 1 }
//   relax:  { method: lbfgs, fmax: 0.5, pushoff: true, target_density: 0 }
//   md:     { ps: 100, temperature: 300, ensemble: npt, pressure: 1, dt: 1 }
//   equilibrate: { protocol: larsen21, t_max: 600, t_final: 300, p_max: 50000 (bar), time_scale: 1 }
//   analyze: { properties: [density, rg] }
//   export: [lammps, gromacs, pdb, xyz, mol2]
//   electrostatics: dsf | pme          cutoff: 10          seed: 1
//
// Exit codes (the CLI's): 0 done · 2 the recipe or an input is wrong · 3 the force field lacks parameters · 4 a run failed.
#pragma once
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/json.hpp"
#include "caps/properties.hpp"
#include "caps/provenance.hpp"
#include "caps/system.hpp"

namespace caps {

struct RecipeError : std::runtime_error {
  int code;
  RecipeError(int c, const std::string& what) : std::runtime_error(what), code(c) {}
};

struct RecipeEvent {
  int stage = 0, stages = 0;   // 1-based
  std::string name;            // build, type, grow, relax, md, equilibrate, analyze, export
  std::string detail;
  std::string status;          // running, done, failed, waiting
  double fraction = 0;
};

struct RecipeOptions {
  std::string base_dir = ".";        // relative paths in the recipe
  std::string out_dir = ".";         // exported files
  std::string forcefield_dir;        // the library (catalogue.json) for force-field names
  long long seed = -1;               // ≥ 0 overrides every seed of the recipe
  int threads = 0;                   // 0: automatic
  std::function<void(const RecipeEvent&)> progress;
};

struct RecipeResult {
  std::string name;
  System system;
  std::vector<Property> properties;
  std::vector<std::string> files;
  Manifest manifest;
  std::shared_ptr<const ForceField> field;   // the typed force field (null when nothing typed it)
  std::string forcefield;                    // its name
};

// The stages the recipe will run, in order (checks the keys; throws RecipeError 2).
std::vector<std::string> recipe_stages(const Json& recipe);
RecipeResult run_recipe(const Json& recipe, const RecipeOptions& o);

}  // namespace caps
