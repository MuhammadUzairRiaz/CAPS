// CAPS core data model: one frame of atoms in a (possibly triclinic) periodic cell.
#pragma once
#include <array>
#include <cstdint>
#include <optional>
#include <memory>
#include <string>
#include <vector>

namespace caps {

using Vec3 = std::array<double, 3>;

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Vec3 operator*(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a);

// Periodic cell spanned by three edge vectors from an origin.
struct Cell {
  Vec3 origin{0, 0, 0};
  Vec3 a{0, 0, 0}, b{0, 0, 0}, c{0, 0, 0};
  std::array<bool, 3> periodic{true, true, true};
  bool valid() const;
  double volume() const;
  Vec3 to_fractional(const Vec3& r) const;
  Vec3 to_cartesian(const Vec3& f) const;
  Vec3 minimum_image(Vec3 d) const;   // shortest periodic image of a separation
  Vec3 wrap(const Vec3& r) const;     // position folded into the cell
};

struct Atom {
  int64_t id = 0;        // file identifier (1-based in most formats)
  int64_t mol = 0;       // molecule identifier, 0 = none
  int type = 0;          // numeric type, 0 = none
  std::string name;      // type or atom name from the file
  std::string resname;   // residue name (GROMACS, PDB)
  int64_t resid = 0;     // residue number (GROMACS, PDB, the peptide builder), 0 = none
  int element = 0;       // atomic number, 0 = unknown
  double charge = 0.0;
  Vec3 pos{0, 0, 0};     // Å, as read (may be unwrapped)
  std::array<int, 3> image{0, 0, 0};
};

struct Bond {
  uint32_t i, j;
  int order = 0;   // 0 unknown, 1 single, 2 double, 3 triple, 4 aromatic, 5 amide (Tripos mol2 codes)
};

struct TypeInfo {
  int type = 0;
  double mass = 0.0;
  std::string label;     // e.g. "c3" from "1 12.011 # c3"
};

// A molecule topology given term by term (a coarse-grained protein, as martinize writes one): parameterize uses these
// bonded terms, with their own parameters, in place of rule lookups; the non-bonded terms still come from the types.
// Units as LAMMPS real: bonds K (r − r0)² (kcal/mol/Å², Å); angles form 1 cosine/squared K (cos θ − cos θ0)²; dihedrals
// form 1 K [1 + cos(nφ − φ0)], form 2 harmonic K (ξ − ξ0)² on the i-j-k-l dihedral (radians). Valid only for the
// structure it was made with (natoms and the bond list must still match; parameterize checks).
struct ExplicitTopology {
  struct Bond { uint32_t i, j; double k, r0; std::string group; };
  struct Angle { uint32_t i, j, k; int form; double kt, theta0; std::string group; };
  struct Dihedral { uint32_t i, j, k, l; int form; double kd, phi0; int n; std::string group; };
  size_t natoms = 0;
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Dihedral> dihedrals;
  std::string source;   // "Martini 2.2 protein (martinize rules)"
};

struct System {
  std::string title;
  std::string source_format;
  int64_t timestep = -1;
  Cell cell;
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;
  std::vector<TypeInfo> types;
  std::vector<Vec3> velocities;   // Å/fs, one per atom, or empty
  bool bonds_from_file = false;
  bool has_charges = false;
  bool has_mol = false;
  bool unwrapped = false;   // positions are unwrapped (xu/yu/zu or image flags applied)
  std::vector<std::string> notes;  // what the reader inferred or skipped
  std::shared_ptr<const ExplicitTopology> topology;   // bonded terms given term by term, or null

  double mass_of(const Atom& a) const;
  double total_mass() const;          // g/mol
  double density() const;             // g/cm^3, 0 if no valid cell
  std::vector<std::vector<uint32_t>> neighbours() const;
  // Molecule index per atom (0..n-1): from mol ids when present, else bond-connected components.
  std::vector<int> molecules(int* count = nullptr) const;
};

// One trajectory: shared topology, positions per frame.
struct Trajectory {
  System topology;                    // frame 0 plus bonds/types
  std::vector<std::vector<Vec3>> positions;
  std::vector<Cell> cells;
  std::vector<int64_t> timesteps;
  size_t frames() const { return positions.size(); }
  System frame(size_t k) const;
};

}  // namespace caps
