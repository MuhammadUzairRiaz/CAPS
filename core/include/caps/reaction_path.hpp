// CAPS reaction paths between two structures of the same atoms: gas-hop barriers in rubber, chain slip in crystals,
// vitrimer exchange barriers.
//
//  match_atoms     which atom of B is each atom of A when the files list them in another order: atoms of one element
//                  (and, when both structures have the same bonded environments, one environment, refined by the
//                  neighbours' environments) are assigned by the Hungarian method (Kuhn, Naval Res. Logist. Q. 2, 83
//                  (1955)) to minimise the summed squared distances, and, with the environments, then turned into the
//                  nearest map that keeps every bond (each molecule walked breadth-first, every atom given a free
//                  partner bonded to its placed neighbours' partners, the Hungarian one first, then the nearest;
//                  backtracking within a bounded number of tries, a note when none is found); a molecule is superposed first (Horn's quaternion
//                  fit, from the four principal-axis orientations, the assignment and the fit repeated until stable), a
//                  periodic cell uses minimum-image distances.
//  interpolate     the images between A and B (B in A's atom order, each atom at its image nearest A's; a molecule's B
//                  superposed on A first, so the path carries no rigid motion):
//                  linear     straight lines in Cartesian space;
//                  lst        each image minimises S = Σ_{i<j} (d_ij^κ − d_ij)² / d_ij⁴ + γ Σ_i |r_i − r_i^LI|², the
//                             distances d^κ interpolated linearly between the end points, γ = 10⁻⁶ bohr⁻⁴ (linear
//                             synchronous transit, Halgren & Lipscomb, Chem. Phys. Lett. 49, 225 (1977), as restated
//                             by Smidstrup et al. 2014);
//                  idpp       the same distance objective without the γ term, the images optimised together as a nudged
//                             elastic band on it (image dependent pair potential, Smidstrup, Pedersen, Stokbro &
//                             Jónsson, J. Chem. Phys. 140, 214106 (2014)).
//                  Pairs by the minimum image in periodic cells.
//  neb             climbing-image nudged elastic band (Henkelman, Uberuaga & Jónsson, J. Chem. Phys. 113, 9901 (2000))
//                  with the improved tangent (Henkelman & Jónsson, J. Chem. Phys. 113, 9978 (2000)): on each image the
//                  force perpendicular to the path plus a spring force k(|R_{i+1} − R_i| − |R_i − R_{i−1}|) τ along it;
//                  the highest image climbs (F − 2(F·τ)τ, no spring) once the band has settled; the images move by FIRE
//                  (Bitzek et al. 2006). Any energy: the force field (neb_forcefield) or a function of the coordinates.
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct AtomMatch {
  std::vector<int> map;          // map[i]: the atom of B that is atom i of A
  double rmsd = 0;               // Å, after the superposition (molecules) or by minimum image (cells)
  int reordered = 0;             // atoms whose index differs
  bool by_graph = false;         // bonded environments were used
  std::vector<std::string> notes;
};
AtomMatch match_atoms(const System& a, const System& b);
// B's atoms in A's order (positions and everything else from B), each placed at its periodic image nearest A's atom.
System reorder_atoms(const System& b, const std::vector<int>& map, const System& a);

struct PathOptions {
  std::string method = "idpp";   // linear | lst | idpp
  int images = 7;                // intermediate images
  double tolerance = 1e-3;       // objective force at which lst / idpp images stop (Å⁻³, as the objective's gradient)
  int max_steps = 3000;
};
struct PathReport {
  std::vector<double> min_distance;   // per image (ends included): the shortest interatomic distance, Å
  int steps = 0;
  std::vector<std::string> notes;
};
// The images from A to B (ends included); B must already be in A's atom order.
std::vector<std::vector<Vec3>> interpolate_path(const System& a, const System& b, const PathOptions& o = {}, PathReport* report = nullptr);

struct NebOptions {
  double spring = 5.0;           // energy per Å² (kcal/mol/Å² with a force field)
  bool climb = true;
  double fmax = 0.05;            // largest atomic NEB force at convergence (energy per Å)
  int max_steps = 2000;
  double climb_after = 0.5;      // the climbing image starts when the largest force falls below this (or at half the steps)
};
struct NebResult {
  std::vector<std::vector<Vec3>> images;   // ends included
  std::vector<double> energy, arc;          // per image; arc length along the path (Å)
  int climbing = -1;                        // the climbing image
  double forward = 0, reverse = 0;          // the highest intermediate image above A and above B
  double fmax = 0;
  int steps = 0;
  bool converged = false;
  std::vector<std::string> notes;
};
// image index, coordinates → energy, forces (−gradient)
using ImageEnergy = std::function<double(int image, const std::vector<Vec3>& x, std::vector<Vec3>& f)>;
NebResult neb(std::vector<std::vector<Vec3>> path, const ImageEnergy& energy, const NebOptions& o = {});
// The band with a force field (the ends' energies included); positions in Å, energies in kcal/mol.
NebResult neb_forcefield(const System& a, std::vector<std::vector<Vec3>> path, const ForceField& ff, const EnergyOptions& eo, const NebOptions& o = {});

}  // namespace caps
