// CAPS core data model: one frame of atoms in a (possibly triclinic) periodic cell.
#pragma once
#include <array>
#include <cstdint>
#include <optional>
#include <memory>
#include <map>
#include <string>
#include <vector>

namespace caps {

struct ForceField;

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
  uint32_t tags = 0;     // the System's tags this atom carries (bit k: tags[k]); see caps/tags.hpp
  int64_t chain = 0;     // the molecule it started in before reactions joined molecules (React sets it; 0: not recorded)
};

struct Bond {
  uint32_t i, j;
  int order = 0;   // 0 unknown, 1 single, 2 double, 3 triple, 4 aromatic, 5 amide (Tripos mol2 codes), 9 coordinate (dative)
};
// A coordinate (dative) bond, a ligand donor's lone pair to a metal (MDL V3000 bond type 9): connectivity, but no share in
// either atom's valence (hydrogens, bond orders and formal charges ignore it).
constexpr int kBondDative = 9;

struct TypeInfo {
  int type = 0;
  double mass = 0.0;
  std::string label;     // e.g. "c3" from "1 12.011 # c3"
};

// A molecule topology given term by term (a coarse-grained protein, as martinize writes one): parameterize uses these
// bonded terms, with their own parameters, in place of rule lookups; the non-bonded terms still come from the types.
// Units as LAMMPS real: bonds K (r − r0)² (kcal/mol/Å², Å); angles form 0 harmonic K (θ − θ0)², 1 cosine/squared
// K (cos θ − cos θ0)², 5 restricted bending K (cos θ − cos θ0)² / sin² θ; dihedrals form 1 (and 9, GROMACS's multiple form) K [1 + cos(nφ − φ0)], form 2
// harmonic K (ξ − ξ0)² on the i-j-k-l dihedral, form 4 periodic improper K [1 + cos(nφ − φ0)] (radians), form 11 combined
// bending–torsion sin³θ1 sin³θ2 Σ c_n cos^n φ; pairs: explicit LJ pairs; exclusions: pairs
// with no non-bonded interaction beyond the bonded ones; virtual sites: the centre of mass of their atoms (GROMACS
// virtual_sitesn 2), massless. Valid only for the structure it was made with (natoms and the bond list must still match;
// parameterize checks).
struct ExplicitTopology {
  struct Bond { uint32_t i, j; double k, r0; std::string group; };
  struct Angle { uint32_t i, j, k; int form; double kt, theta0; std::string group; };
  struct Dihedral { uint32_t i, j, k, l; int form; double kd, phi0; int n; std::string group; std::array<double, 5> c{}; };   // c: form 11's coefficients
  struct Pair { uint32_t i, j; double eps, sigma; };   // an explicit LJ pair (kcal/mol, Å)
  size_t natoms = 0;
  std::vector<Bond> bonds;
  std::vector<Angle> angles;
  std::vector<Dihedral> dihedrals;
  std::vector<std::pair<uint32_t, uint32_t>> exclusions;
  std::vector<Pair> pairs;
  struct VSite { uint32_t site; std::vector<uint32_t> from; std::vector<double> w; double c = 0; };   // w empty: the centre of mass; c: out of plane (1/Å)
  std::vector<VSite> vsites;
  std::vector<double> masses;   // per atom when the molecule sets them (NaN: the type's), else empty
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
  std::shared_ptr<const ForceField> forcefield;       // the force field the file carries for these atoms (AMBER prmtop), or null
  struct Tag { std::string name, colour; };           // named atom sets (design/boards/Tags): bit k of Atom::tags
  std::vector<Tag> tags;

  double mass_of(const Atom& a) const;
  double total_mass() const;          // g/mol
  double density() const;             // g/cm^3, 0 if no valid cell
  std::vector<std::vector<uint32_t>> neighbours() const;
  // Molecule index per atom (0..n-1): from mol ids when present, else bond-connected components.
  std::vector<int> molecules(int* count = nullptr) const;
};

// Which frames of a trajectory file to keep (file frame numbers from 0): first, first + stride, … up to last (inclusive).
// Frames not kept are skipped as the file is read, so a long trajectory is held at the size of what is kept.
struct FrameSelection {
  size_t first = 0, last = SIZE_MAX, stride = 1;
  bool all() const { return first == 0 && last == SIZE_MAX && stride <= 1; }
  bool wants(size_t k) const { return k >= first && k <= last && (k - first) % (stride ? stride : 1) == 0; }
  bool past(size_t k) const { return k >= last; }   // frame k was the last one wanted: stop reading
};

// One trajectory: shared topology, positions per frame.
struct Trajectory {
  System topology;                    // frame 0 plus bonds/types
  std::vector<std::vector<Vec3>> positions;
  std::vector<Cell> cells;
  std::vector<int64_t> timesteps;
  // per frame, when the file gives them: velocities (a LAMMPS dump's vx vy vz, as written: Å/fs in real units), and any
  // other per-atom numeric column by its name (fx, c_pe, v_stress …; |f| and |v| added from the components), atoms in
  // the topology's order
  std::vector<std::vector<Vec3>> velocities;
  std::map<std::string, std::vector<std::vector<float>>> columns;
  size_t frames() const { return positions.size(); }
  System frame(size_t k) const;
  // Reading with a selection: a reader sets `selection`, appends each frame it reads (positions, cell, timestep and
  // any velocities or columns) and calls admit(), which drops that frame when it is not wanted and returns false once
  // the selection's last frame has been read. frames_read counts the file's frames seen.
  FrameSelection selection;
  size_t frames_read = 0;
  bool admit();
  void drop_last();
  // Keeps the selected frames of those held (frames already in memory; the readers use admit instead).
  void select(const FrameSelection& sel);
};

// Held atoms (RelaxOptions::fixed, DynamicsOptions::fixed): 1 holds every coordinate; otherwise bits 2, 4, 8 hold x, y, z
// on their own (a substrate free to slide in its plane but not to leave it: 8).
inline bool holds_axis(char f, int k) { return f == 1 || ((f >> (k + 1)) & 1); }
inline bool holds_all(char f) { return f == 1 || (f & 14) == 14; }

}  // namespace caps
