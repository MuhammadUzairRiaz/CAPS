// CAPS force-field phonons: lattice dynamics of a crystal with its force field — dispersion along a path through the
// Brillouin zone, the density of states on a k-point mesh, the harmonic thermodynamics, and sound speeds.
//
//  Force constants  the cell repeated n1 × n2 × n3 times (each width at least twice the cut-off unless given), each atom
//                   of the cell displaced ±h along x, y, z, the forces on every atom of the supercell by central
//                   differences: Φ_aα,Jβ = −∂F_Jβ/∂u_aα (the pair list frozen at the reference, so the derivatives stay
//                   smooth). The acoustic sum rule Σ_J Φ_aα,Jβ = 0 is imposed on the self term.
//  Images           each pair takes the supercell image(s) nearest a: equally near ones share the weight (Wigner–Seitz).
//  Dynamical matrix D_aα,bβ(k) = Σ_J∈b Σ_images w Φ_aα,Jβ e^(i k·(r_J + L − r_a)) / √(m_a m_b), Hermitian; diagonalised as
//                   the real symmetric [[Re, −Im], [Im, Re]] (each eigenvalue twice). ν̃ = √λ / (2πc), 108.59 √λ cm⁻¹ for
//                   λ in kcal/mol/Å²/(g/mol); λ < 0 (not at a minimum) as negative wavenumbers.
//  DOS              Gaussian-broadened over a Γ-centred mesh, normalised to 3N modes per cell; per element from |e|².
//  Thermodynamics   quantum harmonic oscillators summed over the mesh, per cell: zero-point energy, U, F, S, C_v.
//  Sound speeds     the three acoustic branches' slopes dω/dk at a small k along a*, b* and c*.
// Born & Huang, Dynamical Theory of Crystal Lattices (1954); the supercell finite-displacement method as phonopy (Togo &
// Tanaka, Scr. Mater. 108, 1 (2015)).
#pragma once

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct PhononOptions {
  // the force field of any structure built from the cell (the supercell); required
  std::function<ForceField(const System&)> field_for;
  EnergyOptions energy;
  double step = 0.01;                          // Å, finite-difference displacement
  // relax the atoms first, in the supercell's periodic copies (the cell kept): a cell narrower than twice the cut-off
  // cannot be relaxed on its own, it sees one image of each neighbour; FIRE (Bitzek et al. 2006) on the cell's own
  // coordinates with each atom's force summed over its copies, until every one is below relax_ftol (kcal/mol/Å)
  bool relax = true;
  double relax_ftol = 1e-3;
  int relax_iterations = 20000;
  std::array<int, 3> supercell{0, 0, 0};       // repeats (0: each width ≥ 2 × cut-off)
  size_t max_atoms = 15000;                    // the automatic supercell stays below this
  // the path: points in fractional reciprocal coordinates with labels (empty: Γ X S Y Γ Z R of the cell's axes)
  std::vector<std::array<double, 3>> path;
  std::vector<std::string> labels;
  int path_points = 120;                       // over the whole path
  std::array<int, 3> mesh{0, 0, 0};            // DOS mesh (0: about 40 Å / cell length each)
  double dos_sigma = 5;                        // cm⁻¹, Gaussian
  double dos_step = 1;                         // cm⁻¹
  std::vector<double> temperatures{100, 200, 298.15, 300, 400, 500};   // K (0: the zero-point values)
  std::function<bool(const std::string&, double)> progress;   // false stops
};

struct PhononResult {
  std::array<int, 3> supercell{1, 1, 1};
  size_t supercell_atoms = 0;
  double static_energy = 0;                    // kcal/mol per cell at the reference (the copies divided out)
  double max_force = 0;                        // kcal/mol/Å at the reference: should be small (a minimum)
  double relax_shift = 0;                      // Å, the largest atom move of the relaxation
  int relax_steps = 0;
  std::vector<Vec3> positions;                 // the cell's atoms at the reference (relaxed)
  // dispersion: x along the path (Å⁻¹), the 3N wavenumbers at each point (ascending), ticks at the labelled points
  std::vector<double> path_x;
  std::vector<std::vector<double>> bands;
  std::vector<double> tick_x;
  std::vector<std::string> tick_labels;
  // DOS (states per cm⁻¹ per cell) and per element
  std::vector<double> dos_w, dos;
  std::vector<std::string> pdos_labels;
  std::vector<std::vector<double>> pdos;
  std::array<int, 3> mesh{1, 1, 1};
  // thermodynamics per cell (kcal/mol, cal/mol/K), at each temperature
  std::vector<double> T, zpe, U, F, S, Cv;
  // sound speeds (m/s) along a*, b*, c*: the three acoustic branches, slowest first
  std::array<std::array<double, 3>, 3> sound{};
  double min_wavenumber = 0;                   // the lowest over the mesh and path (negative: unstable)
  int imaginary = 0;                           // mesh modes below −10 cm⁻¹ (Γ's acoustic modes excluded)
  std::string field;
  std::vector<std::string> notes;
};

PhononResult phonons(const System& cell, const PhononOptions& o);

}  // namespace caps
