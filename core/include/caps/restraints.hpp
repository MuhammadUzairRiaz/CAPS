// CAPS restraints: harmonic and flat-bottomed penalties on distances, angles, torsions and positions, for dynamics and
// minimisation, written into LAMMPS and GROMACS inputs.
//
//  Energy     E = ½ k (q − q_lo)² below the lower bound, ½ k (q − q_hi)² above the upper one, 0 between (flat bottom);
//             lower = upper is the plain harmonic ½ k (q − q0)². q: a distance (Å, minimum image), an angle i-j-k
//             (degrees, k per rad²), a torsion i-j-k-l (IUPAC sign, degrees; the deviation from the centre of the bounds
//             taken on the circle) or an atom's displacement from its anchor over chosen axes (Å; the flat bottom a
//             sphere of radius `upper`). k ramps linearly to k_end over a run when k_end ≥ 0 (release in stages).
//  LAMMPS     fix restrain (bond K (r − r0)², angle K (θ − θ0)² with K = k/2; dihedral K [1 + cos(φ − d)], d = φ0 + 180°
//             in LAMMPS's own reckoning, which measures φ 180° from the IUPAC sign, K = k: equal to ½ k Δφ² near the
//             minimum), lbound for a lower bound alone, spring/self (½ K r², K = k) for position restraints on the
//             atoms where the data file puts them; Kstart → Kstop carries the ramp. Upper and two-sided flat bottoms
//             other than a distance's lower bound have no LAMMPS form and are left out, with a note.
//  GROMACS    in [ intermolecular_interactions ] by global atom number: bonds function 10 (low, up1, up2 = up1 + 100 nm,
//             k: the same piecewise ½ k form), angles function 1 (harmonic), dihedral_restraints (φ0 the centre, Δφ the
//             half-width: ½ k (|φ′| − Δφ)² outside); [ position_restraints ] function 1 (k_x k_y k_z) or 2 (flat-bottomed
//             sphere, g = 1) inside each molecule type, the anchors in a reference file for grompp -r. No ramps; a
//             flat-bottomed angle has no GROMACS form and is left out, with a note.
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct Restraint {
  enum class Kind { Distance, Angle, Torsion, Position };
  Kind kind = Kind::Distance;
  std::array<uint32_t, 4> atoms{0, 0, 0, 0};   // 0-based; 2, 3, 4 or 1 used
  double lower = 0, upper = 0;                 // Å or degrees; position: upper = flat-bottom radius (0: harmonic)
  double k = 10;                               // kcal/mol/Å² (distance, position) or kcal/mol/rad² (angle, torsion)
  double k_end = -1;                           // ≥ 0: k ramps to this over a run
  Vec3 anchor{0, 0, 0};                        // position: Å
  int axes = 7;                                // position: bit 0 x, 1 y, 2 z
};

int restraint_atom_count(Restraint::Kind k);
const char* to_string(Restraint::Kind k);

// The restrained coordinate now: Å, or degrees (torsion in (−180, 180]).
double restraint_value(const Restraint& r, const std::vector<double>& x, const Cell& cell);

struct RestraintTerms {
  double energy = 0;                            // kcal/mol
  double virial = 0;                            // Σ r·f of the internal restraints (positions excluded)
  double w[6] = {0, 0, 0, 0, 0, 0};             // xx yy zz xy xz yz
};
// Adds the restraint forces to f (3N, kcal/mol/Å) and returns the terms; ramp ∈ [0, 1] is the run's progress for k → k_end.
RestraintTerms add_restraint_forces(const std::vector<Restraint>& rs, const std::vector<double>& x, const Cell& cell, std::vector<double>& f, double ramp = 0);

// Text form, atoms numbered from 1, items separated by ';' or new lines:
//   distance I J R|LO:HI [K] [ramp KEND]            (K default 10 kcal/mol/Å²)
//   angle I J K THETA|LO:HI [K] [ramp KEND]          (degrees; K default 50 kcal/mol/rad²)
//   torsion I J K L PHI|LO:HI [K] [ramp KEND]        (degrees, IUPAC; K default 50 kcal/mol/rad²)
//   position SEL [K] [axes xyz] [radius R] [ramp KEND]   (SEL: 5, 1-20, 3,7-9; anchors where the atoms are in s)
// Throws std::invalid_argument naming the item it cannot read.
std::vector<Restraint> parse_restraints(const std::string& text, const System& s);
std::string restraints_text(const std::vector<Restraint>& rs);   // back to that form (positions with equal settings as ranges)
Json restraints_json(const std::vector<Restraint>& rs, const std::vector<double>* x = nullptr, const Cell* cell = nullptr);

// LAMMPS commands for the input (after read_data); energy_scale: 1 real units, 1/23.060549 metal. Notes for what is left out.
std::string lammps_restraint_lines(const std::vector<Restraint>& rs, const System& s, double energy_scale, std::vector<std::string>* notes);

// GROMACS: the [ intermolecular_interactions ] block (global numbering) to end the .top with; notes for what is left out.
std::string gromacs_intermolecular_restraints(const std::vector<Restraint>& rs, std::vector<std::string>* notes);
// [ position_restraints ] lines of the atoms first..last (local numbering from first), "" when none.
std::string gromacs_position_restraints(const std::vector<Restraint>& rs, uint32_t first, uint32_t last, std::vector<std::string>* notes);

}  // namespace caps
