// CAPS space groups (design/boards/CrystalBuilder): the 530 settings of the 230 space groups from their Hall symbols
// (International Tables Vol. B, Table A1.4.2.7), the operations each generates, and a crystal built from a space
// group, a lattice and an asymmetric unit.
//
//  hall_operations   Hall (1981) symbols: lattice centring (with - for a centre of inversion at the origin), rotation
//                    matrices with their axes (default z, then x or a−b, then the body diagonal), translation symbols
//                    a b c n u v w d and screw subscripts, and an origin shift (vx vy vz)/12. The group is closed by
//                    multiplication, translations reduced modulo 1.
//  build_crystal     Sites expanded by the operations; images closer than the tolerance merged (the site's multiplicity
//                    is its number of distinct images); supercell repeats.
//  find_symmetry     The highest-order setting whose operations map the structure onto itself (same elements, within the
//                    tolerance), trying origin shifts; its asymmetric unit.
//  primitive_cell    The primitive cell of a centred lattice (A B C I F, R in hexagonal axes).
#pragma once
#include <array>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct SymOp {
  double R[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  double t[3] = {0, 0, 0};
};
SymOp parse_symop(std::string s);          // "x,y,z+1/2"
std::string symop_string(const SymOp& op);

struct SpaceGroupSetting {
  std::string key;     // "62", "5:b1", "227:2", "166:h"
  int number = 0;
  std::string hm;      // "P n m a", "F d -3 m:2"
  std::string hall;    // "-p 2ac 2n"
};
const std::vector<SpaceGroupSetting>& space_group_settings();
// By setting key ("227:2"), number ("62": the default setting, origin choice 2 where there are two, hexagonal axes for
// rhombohedral groups) or Hermann–Mauguin symbol with or without spaces ("Pnma", "P 21/c", "Fd-3m").
const SpaceGroupSetting* find_space_group(const std::string& query);
std::vector<SymOp> hall_operations(const std::string& hall);
std::string crystal_system(int number);   // triclinic … cubic

struct CrystalSite {
  std::string label;
  int element = 0;
  Vec3 frac{0, 0, 0};
  double occupancy = 1;
};

struct CrystalSpec {
  std::string space_group = "P 1";
  double a = 5, b = 5, c = 5, alpha = 90, beta = 90, gamma = 90;   // Å, degrees
  std::vector<CrystalSite> sites;
  std::array<int, 3> supercell{1, 1, 1};
  double tolerance = 0.01;   // Å: images closer than this are one atom
  std::string title;
};

struct CrystalReport {
  std::string key, hm, hall, system;
  int number = 0, operations = 0;
  size_t atoms_per_cell = 0;
  std::vector<int> multiplicity;   // per site
  double volume = 0;               // Å³ of the unit cell
  std::vector<std::string> notes;
};
System build_crystal(const CrystalSpec& spec, CrystalReport* report = nullptr);
// A cell from its parameters (Å, degrees): a along x, b in the xy plane.
Cell cell_parameters(double a, double b, double c, double alpha, double beta, double gamma);
// Sites moved onto their special positions: each site becomes the mean of its images that land within `snap` Å of it
// (its site-symmetry group), so x = 0.333 in P 63 m c becomes 1/3. `moved` counts the sites that changed.
CrystalSpec symmetrize_sites(const CrystalSpec& spec, double snap = 0.3, int* moved = nullptr);

struct SymmetryFound {
  std::string key, hm, system;
  int number = 0, operations = 0;
  Vec3 origin{0, 0, 0};              // fractional shift to the setting's origin
  std::vector<CrystalSite> sites;    // asymmetric unit in the setting's origin
};
SymmetryFound find_symmetry(const System& crystal, double tolerance = 0.1);

// The primitive cell of a centred lattice; centring 'A' 'B' 'C' 'I' 'F' or 'R' (hexagonal axes). 'P' returns a copy.
System primitive_cell(const System& crystal, char centring);
System supercell(const System& s, int nx, int ny, int nz);

}  // namespace caps
