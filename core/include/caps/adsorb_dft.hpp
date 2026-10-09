#pragma once
// A molecule on a 2D slab for DFT (the DFT surface & adsorption workbench, "adsorption set"): anchors found by SMARTS,
// orientations (anchor down, parallel, upright) × azimuths, the azimuths the outer layer's symmetry makes equivalent,
// the molecule lowered until the closest contact reaches dmin (any pair) or dheavy (heavy atoms), its distance to its
// own periodic images on the unwrapped molecule (with the smallest supercell that passes), and an independent audit of
// written complexes. The trio — slab, molecule in a box, complexes — is written with identical settings by vasp_set.
#include <map>
#include <string>
#include <vector>

#include "caps/slab2d.hpp"
#include "caps/system.hpp"

namespace caps {

struct AnchorDef { std::string name, smarts; };
// nitrile N, C=C, C=O oxygen, hydroxyl O, NH2 N, aromatic ring, thiol S, carboxylic O (data/sheets/anchors.json)
std::vector<AnchorDef> anchor_library(const std::string& data_dir);
// The anchors present in a molecule: name → the atoms of its first match (a connected group).
std::map<std::string, std::vector<int>> find_anchors(const System& mol, const std::vector<AnchorDef>& defs);

struct AdsorbSetOptions {
  std::vector<std::string> modes;          // "anchor:<name>", "parallel", "upright"; empty: every anchor found + parallel
  std::vector<double> azimuths = {0, 90, 180, 270};
  double dmin = 2.3;                       // Å, closest molecule–slab contact (any atoms)
  double dheavy = 3.0;                     // Å, closest non-H molecule atom to non-H slab atom
  double image_error = 2.5, image_warn = 4.0;   // Å, the unwrapped molecule to its own in-plane images
  double vacuum = 20.0;
  double box_padding = 8.0;                // Å around the free molecule
  std::vector<std::string> order = {"Ti", "C", "N", "O", "F", "H"};
  bool skip_equivalent = false;            // leave out azimuths equivalent by the surface's symmetry
};
struct AdsorbComplex {
  std::string name;                        // complex_<mode>_az<deg>
  std::string mode;
  double azimuth = 0;
  System system;                           // slab atoms first in each element block, then the molecule's
  std::vector<int> molecule;               // the molecule's atoms in `system`
  double any = 0, heavy = 0, image = 0;    // closest distances (Å)
  std::string equivalent_to;               // a complex this one equals by the surface's symmetry ("" none)
  std::vector<Finding> issues;
  std::string status = "ok";               // ok | FAIL
  int suggested_supercell = 0;             // smallest n × n (same shape) whose image distance passes, when this fails
};
struct AdsorbSet {
  System slab;                             // the slab as used (species grouped)
  System molecule;                         // in its box (Γ only)
  std::vector<AdsorbComplex> complexes;
  std::map<std::string, std::vector<int>> anchors;
  std::vector<double> symmetry_rotations;  // rotations (deg) about the placement point that map the outer layer onto itself
  std::vector<std::string> notes;
};
AdsorbSet build_adsorption_set(const System& slab, const System& molecule, const std::map<std::string, std::vector<int>>& anchors, const AdsorbSetOptions& o);

// Pre-screen (the force-field Adsorption Locator, caps/adsorption.hpp): UFF on slab + molecule, the molecule annealed as
// a rigid body above the top face (Monte Carlo simulated annealing), the lowest `keep` configurations turned into DFT
// complexes (complex_ff1 …, mode "prescreen") with the same checks. `energies` gets their UFF interaction (kcal/mol).
struct PrescreenOptions { int keep = 3; int cycles = 3; int steps = 20000; uint64_t seed = 1; double window = 10.0; };
std::vector<AdsorbComplex> prescreen_complexes(const System& slab, const System& molecule, const AdsorbSetOptions& o, const PrescreenOptions& p,
                                               std::vector<double>* energies = nullptr);

// The rotations (multiples of 30°, 0 < θ < 360) about the in-plane point c that map the outer layer of the top face
// (atoms within 0.5 Å of the highest heavy atom, by element) onto itself with the cell's periodicity.
std::vector<double> outer_layer_rotations(const System& slab, const Vec3& centre, double tol = 0.15);

// A relaxed slab (a CONTCAR) made one block along z (the largest z cluster must hold every atom), repeated na × nb,
// re-centred with the vacuum, species grouped.
System slab_from_relaxed(const System& contcar, int na, int nb, double vacuum, const std::vector<std::string>& order = {"Ti", "C", "N", "O", "F", "H"});

// Independent audit of written complexes (its own code path): the molecule = the atoms beyond the slab's composition
// per element, made whole along bonds (< 1.9 Å, minimum image), checked for atom count, contacts, image distance and
// intactness against the free molecule's bonds.
struct AuditRow { std::string name; int n_mol = 0; double any = 0, heavy = 0, image = 0; bool intact = false; bool ok = false; };
std::vector<AuditRow> audit_complexes(const System& slab, const System& molecule, const std::vector<std::pair<std::string, System>>& complexes,
                                      double dmin = 2.25, double dheavy = 2.95, double dimage = 4.0);

}  // namespace caps
