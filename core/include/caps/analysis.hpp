#pragma once
#include <cstdint>
#include <utility>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct BondOptions {
  double tolerance = 0.45;     // Å added to the sum of covalent radii
  double min_distance = 0.4;   // closer pairs are reported, not bonded
  bool skip_hh = true;
};

// Distance-based bonds with periodic minimum image and a cell list. O(N).
std::vector<Bond> perceive_bonds(const System& s, const BondOptions& opt = {});

// Unwrap each molecule so it is whole (follows bonds with minimum-image steps).
void make_molecules_whole(System& s);
/// The positions with every molecule whole (each bonded atom at its partner's nearest image, walked from the molecule's
/// first atom); the structure is left as it is. Image flags written from these are consistent along every bond.
std::vector<Vec3> whole_positions(const System& s);

struct MoleculeShape {
  int molecule;
  int atoms;
  double mass;
  Vec3 com;
  double rg;
  double lambda[3];   // gyration tensor eigenvalues, ascending
  double kappa2;      // relative shape anisotropy
  Vec3 axis[3];       // the principal axes (unit vectors) for lambda[0..2]
};
std::vector<MoleculeShape> molecule_shapes(const System& s);

// Partial g(r) between element sets (0 = any). Pairs within one molecule are skipped when inter_only. Above
// kRdfMaxCentres centre atoms an evenly strided subset of them is used (all atoms still count as neighbours).
constexpr size_t kRdfMaxCentres = 20000;
std::vector<std::pair<double, double>> rdf(const System& s, int elem_a, int elem_b, double rmax, double dr, bool inter_only);

// Distance (2 atoms, Å), angle (3, degrees) or dihedral (4, degrees, IUPAC sign), minimum image.
double measure(const System& s, const std::vector<uint32_t>& idx);

// Backbone of each molecule with at least `min_atoms` heavy atoms: the longest path through heavy atoms that are not
// in rings (vinyl polymers), or through all heavy atoms when rings are part of the main chain. Atom indices in order.
std::vector<std::vector<uint32_t>> backbones(const System& s, int min_atoms = 4);

// Mean-square internal distances along the backbones (Theodorou and Suter, Macromolecules 1985): for each separation n,
// ⟨R²(n)⟩ / (n ⟨b²⟩), which rises towards C∞ for equilibrated chains. Positions must be whole molecules.
struct InternalDistances {
  std::vector<int> n;
  std::vector<double> ratio;    // ⟨R²(n)⟩ / (n ⟨b²⟩)
  double b2 = 0;                // mean squared backbone bond length, Å²
  double r2_end = 0;            // mean squared end-to-end distance of the backbones, Å²
  int chains = 0;
};
InternalDistances internal_distances(const System& s);

void symmetric_eigen3(const double A[3][3], double w[3], double V[3][3]);

}  // namespace caps
