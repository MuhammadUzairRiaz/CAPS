// CAPS torsion scan (design/boards/TorsionScan): the energy along a dihedral, rigid or with the rest relaxed, and the
// conformers it finds.
//
//  Scan        φ(a-b-c-d) set by rotating the atoms on d's side of the b–c bond about it (the bond must not be in a
//              ring), from `from` to `to` in `step`. Rigid: the energy of each rotated geometry. Relaxed: a, b, c and d
//              held and every other atom minimised (L-BFGS) at each step, starting from the previous step.
//  Conformers  local minima of the scanned curve (the full turn is periodic), refined by a parabola through the three
//              lowest points; named trans (|φ| ≥ 150°), anticlinal± (90–150°), gauche± (30–90°) or cis (< 30°).
#pragma once
#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct TorsionScanOptions {
  std::array<int, 4> atoms{-1, -1, -1, -1};
  double from = -180, to = 180, step = 15;   // degrees
  bool relax = false;
  double ftol = 0.1;                          // kcal/mol/Å for the relaxed scan
  EnergyOptions energy;
  std::function<bool(int, int)> progress;     // (done, total) → false stops
};

struct TorsionPoint {
  double phi = 0;            // degrees, as measured after setting
  double energy = 0;         // kcal/mol
  EnergyTerms terms;
  std::vector<Vec3> positions;
};

struct TorsionConformer {
  double phi = 0, energy = 0;   // degrees, kcal/mol relative to the lowest
  std::string state;             // trans, gauche+, gauche−, anticlinal+, anticlinal−, cis
};

struct TorsionScanResult {
  std::vector<TorsionPoint> points;
  std::vector<TorsionConformer> conformers;   // lowest first
  double minimum = 0, barrier = 0;            // kcal/mol: the lowest point, the highest above it
  double phi_start = 0;                       // the dihedral before the scan
  std::vector<uint32_t> moving;               // the rotated atoms
  std::vector<std::string> notes;
};

// Atoms on c's side of the b–c bond (c included); throws when b–c is in a ring or not a bond.
std::vector<uint32_t> moving_side(const System& s, int b, int c);
// Rotates the moving atoms about b→c so φ(a-b-c-d) becomes `phi` degrees.
void set_dihedral(System& s, const std::array<int, 4>& atoms, double phi, const std::vector<uint32_t>& moving);
// The IUPAC dihedral a-b-c-d in degrees (minimum image when the cell is valid).
double dihedral_angle(const System& s, const std::array<int, 4>& atoms);
std::string torsion_state(double phi);

TorsionScanResult torsion_scan(const System& s, const ForceField& ff, const TorsionScanOptions& o);

// A torsion to scan when none is chosen: heavy atoms around the middle bond of the longest backbone that is not in a
// ring ({-1, …} when there is none).
std::array<int, 4> default_torsion(const System& s);

}  // namespace caps
