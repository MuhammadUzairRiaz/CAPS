// Lattice tools for crystals and periodic structures (Materials Studio's Build › Crystals and Build › Surfaces, in
// CAPS's own form):
//   transform_cell      a new cell a' b' c' = (a b c)·M for any 3×3 matrix M (integer, or fractional such as ½ for a
//                       centred ↔ primitive change): the atoms of the old lattice that fall in the new cell; the count
//                       must come to N·|det M| (otherwise M does not map the lattice onto itself and it throws)
//   niggli_reduce       the Niggli-reduced cell of a lattice (Křivý & Gruber, Acta Cryst. A32, 297 (1976), with the
//                       relative ε of Grosse-Kunstleve, Sauter & Adams, Acta Cryst. A60, 1 (2004)); returns M
//   find_primitive_cell the primitive cell of the structure's own lattice (pure translations found from the atoms)
//   conventional_cell   the cell with the most symmetry: from the primitive, reduced cell, every cell of lattice vectors at
//                       90°/120° up to four times its volume, each tested by find_symmetry; the one with the most
//                       operations (then the smallest, then the most orthogonal) wins
//   vacuum_slab         a slab of a periodic structure along c: atoms made contiguous across the largest gap in c,
//                       then c lengthened (along its own direction) so the perpendicular height is thickness + vacuum
//   nanowire            a wire cut from a crystal along a lattice direction [uvw]: periodic along it, a cylinder or a
//                       prism of radius R across, vacuum around
#pragma once
#include <array>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

using Mat3 = std::array<std::array<double, 3>, 3>;   // M[r][c]; new vector j = Σ_i old_i · M[i][j] (columns)

System transform_cell(const System& s, const Mat3& m, double tolerance = 0.05);

// The smallest cell of the structure's own lattice: every translation that maps all atoms onto atoms of the same
// element (within `tolerance` Å) is found, and the three shortest that span the lattice at volume V/k become the cell.
// A supercell, or a centred cell, comes back primitive; k = 1 returns the structure unchanged.
System find_primitive_cell(const System& s, double tolerance = 0.1, int* k = nullptr);

struct NiggliResult {
  Mat3 m{};                 // old → reduced (columns as above)
  double a = 0, b = 0, c = 0, alpha = 0, beta = 0, gamma = 0;
  int iterations = 0;
};
NiggliResult niggli_reduce(const Cell& cell);
System niggli_cell(const System& s, NiggliResult* result = nullptr);

struct ConventionalResult {
  Mat3 m{};                 // reduced input → conventional
  std::string hm;           // space group found (Hermann–Mauguin) and its number
  int number = 0, operations = 0, candidates = 0;
};
System conventional_cell(const System& s, double tolerance = 0.1, ConventionalResult* result = nullptr);

struct SlabResult {
  double thickness = 0, vacuum = 0;
};
System vacuum_slab(const System& s, double vacuum, bool centre = true, SlabResult* result = nullptr);

struct WireOptions {
  std::array<int, 3> uvw{0, 0, 1};
  double radius = 10;          // Å (a prism: its inscribed radius)
  int repeats = 1;             // lattice periods along the wire
  std::string shape = "cylinder";   // cylinder | hexagonal | square
  double vacuum = 10;          // Å beyond the wire on each side
};
struct WireResult {
  double period = 0;           // Å, |u a + v b + w c|
  size_t atoms = 0;
};
System nanowire(const System& bulk, const WireOptions& o, WireResult* result = nullptr);

}  // namespace caps
