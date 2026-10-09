// CAPS structure-based coarse-graining of a real polymer (design/boards/CoarseGrained "From a polymer"): an all-atom cell
// (a grown melt, or a trajectory of one) mapped to beads, with a bead model derived from it.
//
// Mapping schemes (centre of mass of each bead's atoms, hydrogens with their heavy atoms):
//   unit           one bead per repeat unit (Milano & Müller-Plathe, J. Phys. Chem. B 109, 18609 (2005) for PS; the
//                  common 1:1 mapping of vinyl polymers)
//   backbone_side  two beads per unit: its backbone atoms and its side group (Harmandaris, Adhikari, van der Vegt &
//                  Kremer, Macromolecules 39, 6708 (2006) for PS; moltemplate's "2-bead polymer" is the generic form);
//                  a unit without a side group is one bead
//   backbone_n     n backbone atoms per bead, side groups with their backbone atom (polyethylene 3:1 and the like)
//   rules          chemistry-aware beads (caps/cg_rules.hpp): bonds cut by SMARTS, fragment SMARTS or an atom → bead list,
//                  named by rules (polyesters: the ester-cut preset's B, S, A, T beads)
//
// The bead model (a starting point, as its report says):
//   bonds, angles   Boltzmann inversion of the mapped distributions over every bond and angle of the chains and every
//                   frame given, as harmonic terms: r0 = ⟨r⟩, k = k_B T / (2 var r) for E = k (r − r0)², the same for θ
//                   (Tschöp, Kremer, Batoulis, Bürger & Hahn, Acta Polym. 49, 61 (1998): direct inversion of bonded
//                   distributions). Types by bead types.
//   non-bonded      purely repulsive WCA with ε = k_B T and σ where the g(r) of the non-bonded bead pairs (other chains,
//                   or more than three bonds apart) first reaches 1/e — there the IBI starting potential −k_B T ln g(r)
//                   is k_B T, as the WCA is at r = σ — cut at 2^(1/6) σ; per type pair when there are several bead
//                   types (LAMMPS cosine/squared … wca), pooled otherwise. The
//                   attractions that give the right density and structure come from iterative Boltzmann inversion
//                   (Reith, Pütz & Müller-Plathe, J. Comput. Chem. 24, 1624 (2003)), which is not done here.
//   exclusions      1-2 and 1-3 bead pairs (LAMMPS special_bonds lj 0 0 1).
#pragma once
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "caps/cg_rules.hpp"
#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct CgMapOptions {
  std::string scheme = "unit";   // unit | backbone_side | backbone_n | rules (caps/cg_rules.hpp: bond cuts, fragments or a list)
  int per_bead = 3;              // backbone_n
  double temperature = 300;      // K: the inversion's k_B T
  CgRules rules;                 // the rules scheme
};

struct CgBondType { std::string a, b; double r0 = 0, k = 0, sd = 0; int count = 0; };
struct CgAngleType { std::string a, b, c; double theta0 = 0, k = 0, sd = 0; int count = 0; };   // θ0, sd in radians

struct CgMapResult {
  System beads;                            // the last frame mapped (beads named by type; element 0)
  std::vector<std::vector<Vec3>> frames;   // every frame mapped
  std::vector<Cell> cells;
  std::shared_ptr<ForceField> ff;          // the bead model for `beads`
  std::vector<CgBondType> bonds;
  std::vector<CgAngleType> angles;
  double sigma = 0, cut = 0, epsilon = 0;  // Å, Å, kcal/mol: the pooled repulsive size
  std::map<std::pair<std::string, std::string>, double> pair_cut;   // several bead types: each pair's own WCA cut (Å)
  double cut_max = 0;
  std::vector<std::pair<double, double>> gr;   // non-bonded bead–bead g(r), pooled
  std::vector<std::string> notes;
};

// aa: the all-atom structure (whole molecules; residues numbered by repeat unit for unit / backbone_side, as Grow
// numbers them); frames/cells: its trajectory (empty: aa's own positions).
CgMapResult cg_map(const System& aa, const CgMapOptions& o, const std::vector<std::vector<Vec3>>& frames = {},
                   const std::vector<Cell>& cells = {});

}  // namespace caps
