// CAPS powder indexing and direct-space structure solution.
//
//  Peaks        local maxima of the smoothed pattern above its background, positions by a parabola through the top three
//               points; a weaker peak where Kα2 of a stronger one falls is dropped.
//  Indexing     the successive dichotomy method (Boultif & Louër, J. Appl. Cryst. 24, 987 (1991)): boxes of cell lengths
//               (and the monoclinic angle) are kept while every observed Q = 1/d² lies, within its tolerance, in the range of
//               Q some hkl takes over the box, and halved until small; the survivors' cells are refined by least squares on Q
//               (linear in the reciprocal metric) and ranked by de Wolff's M20 = Q20 / (2 ε̄ N20) (J. Appl. Cryst. 1, 108
//               (1968)); Smith & Snyder's F_N (J. Appl. Cryst. 12, 60 (1979)) beside it. Cubic, tetragonal, hexagonal,
//               orthorhombic and monoclinic primitive cells (a centred lattice is found as its primitive or a larger cell).
//               As TREOR and DICVOL (Werner et al., J. Appl. Cryst. 18, 367 (1985)), up to a set number of lines may be left
//               unindexed (impurities).
//  Solution     simulated annealing (Kirkpatrick et al., Science 220, 671 (1983)) of a molecule's position, orientation and
//               torsions in a cell and space group against the measured pattern in direct space (David, Shankland &
//               Shankland, Chem. Commun. 931 (1998)): the profile Rwp with the scale and background solved at every trial
//               (refine.hpp), moves that bring atoms of different molecules within a contact distance refused.
#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "caps/refine.hpp"
#include "caps/system.hpp"

namespace caps {

struct PeakOptions {
  double min_height = 0.02;          // fraction of the strongest peak above the background
  int smooth = 2;                    // half-width of the moving average, points
  bool drop_kalpha2 = true;
  double lambda1 = 1.54056, lambda2 = 1.54439;   // Cu Kα1, Kα2 (Å), as powder.hpp
};
struct Peak { double two_theta = 0, height = 0; };
std::vector<Peak> find_peaks(const MeasuredPattern& p, const PeakOptions& o = {});

struct IndexOptions {
  double wavelength = 1.54056;       // Å (Cu Kα1, as powder.hpp)
  double tolerance = 0.03;           // degrees 2θ
  int lines = 20;                    // the first lines used
  int unindexed = 0;                 // lines allowed to stay unindexed
  double max_length = 20;            // Å
  double max_volume = 2500;          // Å³
  double min_beta = 90, max_beta = 125;   // monoclinic angle range, degrees
  std::vector<std::string> systems{"cubic", "tetragonal", "hexagonal", "orthorhombic", "monoclinic"};
  int solutions = 10;
  std::function<bool(const std::string&, double)> progress;
};

struct IndexSolution {
  std::string system;
  std::array<double, 6> cell{};      // a b c (Å) alpha beta gamma (degrees)
  double volume = 0;
  double M20 = 0, FN = 0;
  int indexed = 0;                   // of the lines used
  std::vector<std::array<int, 3>> hkl;   // per line used (0 0 0: unindexed)
  std::vector<double> two_theta_obs, two_theta_calc;
};
struct IndexResult {
  std::vector<double> lines;         // the 2θ positions used
  std::vector<IndexSolution> solutions;   // best M20 first
  std::vector<std::string> notes;
  double seconds = 0;
};
IndexResult index_pattern(const std::vector<double>& two_theta, const IndexOptions& o);

struct SolveOptions {
  std::string space_group = "P 1";   // symbol or number (spacegroup.hpp find_space_group)
  RefineOptions refine;              // radiation, wavelength, profile (U V W X Y, zero), background, range
  int steps = 6000;                  // annealing moves
  double t_start = 0.05, t_end = 1e-4;   // Rwp² "temperatures"
  double contact = 2.0;              // Å, the closest two atoms of different molecules may come
  bool torsions = true;
  int restarts = 2;
  uint64_t seed = 1;
  std::function<bool(const std::string&, double)> progress;
};
struct SolveResult {
  System crystal;                    // the cell with every symmetry copy of the molecule
  double rwp = 0;
  int accepted = 0, steps = 0, torsions = 0;
  std::vector<double> rwp_history;
  std::vector<std::string> notes;
};
// molecule: one molecule (bonds given); cell: a b c α β γ.
SolveResult solve_structure(const MeasuredPattern& data, const System& molecule, const std::array<double, 6>& cell, const SolveOptions& o);

}  // namespace caps
