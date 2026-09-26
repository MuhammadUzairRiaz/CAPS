// CAPS Field (first slice): GAFF typing and parameters for hydrocarbons, and the energy / force evaluator used by Relax.
#pragma once
#include <array>
#include <map>
#include <set>
#include <memory>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct FieldError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// AMBER functional forms (units real: kcal/mol, Å, degrees in the tables, radians inside).
struct BondTerm { uint32_t i, j; double k, r0; };                       // k (r − r0)²
struct AngleTerm { uint32_t i, j, k; double kt, theta0; };              // kt (θ − θ0)², j is the vertex
struct TorsionTerm { uint32_t i, j, k, l; double v; int n; double delta; };   // v [1 + cos(nφ − δ)]
struct PairType { double eps, sigma; };
struct HarmonicTorsion { uint32_t i, j, k, l; double k2, chi0; };          // k2 (χ − χ0)², χ the i-j-k-l dihedral (radians)

// Inversion (LAMMPS improper_style inversion/harmonic, DL_POLY / DL_FIELD "inversion"): c is the centre;
// E = (K/3) Σ (ω − ω0)² over the three angles ω between one bond from c and the plane of the other two (0 ≤ ω ≤ π/2).
// form 0: E = (K/3) Σ (ω − ω0)²; form 1 (DREIDING "plan", LAMMPS umbrella with ω0 = 0 over the three
// permutations at K/3): E = (K/3) Σ (1 − cos ω); form 2 (UFF P, As, Sb, Bi): E = (K/3) Σ [C0 + C1 cos ω + C2 cos 2ω] with
// C2 = 1, C1 = −4 cos ω0, C0 = −(C1 cos ω0 + C2 cos 2ω0) (LAMMPS improper fourier).
struct InversionTerm { uint32_t c, a, b, d; double kw, w0; int form = 0; };

// Other bonded forms. Bonds: form 1 Morse D [1 − e^(−α(r − r0))]² (a = D, b = α, c = r0); form 2 GROMOS quartic
// K (r² − r0²)² as LAMMPS bond gromos (K = ¼ of GROMOS's kb; a = K, b = r0). Angles: form 1 K (cos θ − cos θ0)² (LAMMPS cosine/squared; a = K, b = θ0);
// form 2 K (1 + cos θ) (LAMMPS cosine, linear centres); form 3 (UFF) K [C0 + C1 cos θ + C2 cos 2θ] with C2 = 1/(4 sin²θ0),
// C1 = −4 C2 cos θ0, C0 = C2 (2 cos²θ0 + 1) (LAMMPS fourier; a = K, b = θ0); forms 11–14 (UFF) K (1 − cos nθ)/n² with
// n = form − 10, n = 1 meaning K (1 + cos θ), plus UFF's wall e^(−20 (θ − θ0 + 0.25)) below 30° (LAMMPS cosine/periodic,
// without the wall; b = θ0); form 4 (SDK / SPICA, LAMMPS angle sdk) K (θ − θ0)² plus the repulsive part of the end
// atoms' SDK pair, cut at its minimum and shifted to zero there.
// Urey–Bradley: K (r13 − r0)² between the end atoms of an angle, counted as angle energy (CHARMM).
struct BondX { uint32_t i, j; int form; double a, b, c, d = 0; };   // form 3 FENE: −½ K R0² ln(1 − (r/R0)²) (a = K, b = R0) plus WCA ε, σ (c, d) below 2^(1/6) σ
struct AngleX { uint32_t i, j, k; int form; double a, b; };
struct UreyBradley { uint32_t i, k; double kub, r0; };

// Pair forms other than Lennard-Jones, per type pair: 1 Buckingham A e^(−r/ρ) − C/r⁶ (a = A, b = ρ, c = C);
// 2 Morse D0 [e^(−2α(r − r0)) − 2 e^(−α(r − r0))] (a = D0, b = α, c = r0);
// SDK / SPICA coarse-grained Lennard-Jones (LAMMPS lj/sdk; a = ε, b = σ), C ε [(σ/r)^m − (σ/r)^n]:
// 11 9-6 (C = 27/4), 12 12-4 (C = 3√3/2), 13 12-6 (C = 4), 14 12-5;
// 20 Lennard-Jones 12-6 with GROMACS's force switch from ForceField::lj_inner to the cut-off (LAMMPS lj/gromacs, MARTINI).
struct PairFunc { int form; double a, b, c; };
// 31 / 32 cosine-squared attraction (LAMMPS cosine/squared; Cooke–Deserno): −ε below σ, −ε cos²(π(r − σ) / 2(rc − σ)) up to
// the pair's cut-off rc (a = ε, b = σ, c = rc), 32 with WCA ε[(σ/r)¹² − 2(σ/r)⁶ + 1] below σ; rc = σ: WCA only.
constexpr int kPairSdk96 = 11, kPairSdk124 = 12, kPairSdk126 = 13, kPairSdk125 = 14, kPairGromacs = 20, kPairCos2 = 31, kPairCos2Wca = 32;

// Class II forms (COMPASS, PCFF), as LAMMPS bond / angle / dihedral / improper_style class2. Angles in radians.
struct Class2Bond { uint32_t i, j; double r0, k2, k3, k4; };               // K2 Δr² + K3 Δr³ + K4 Δr⁴
struct Class2Angle {                                                        // j is the vertex
  uint32_t i, j, k;
  double theta0, k2, k3, k4;                                               // K2 Δθ² + K3 Δθ³ + K4 Δθ⁴
  double bb_m, bb_r1, bb_r2;                                               // M (r_ij − r1)(r_jk − r2)
  double ba_n1, ba_n2, ba_r1, ba_r2;                                       // N1 (r_ij − r1) Δθ + N2 (r_jk − r2) Δθ
};
struct Class2Dihedral {
  uint32_t i, j, k, l;
  double k1, phi1, k2, phi2, k3, phi3;                                     // Σ Kn [1 − cos(nφ − φn)]
  double mbt[3], mbt_r2;                                                   // (r_jk − r2) Σ An cos nφ
  double ebt_b[3], ebt_c[3], ebt_r1, ebt_r3;                               // (r_ij − r1) Σ Bn cos nφ + (r_kl − r3) Σ Cn cos nφ
  double at_d[3], at_e[3], at_theta1, at_theta2;                           // (θ_ijk − θ1) Σ Dn cos nφ + (θ_jkl − θ2) Σ En cos nφ
  double aat_m, aat_theta1, aat_theta2;                                    // M (θ_ijk − θ1)(θ_jkl − θ2) cos φ
  double bb13_n, bb13_r1, bb13_r3;                                         // N (r_ij − r1)(r_kl − r3)
};
struct Class2Improper {                                                     // j is the centre
  uint32_t i, j, k, l;
  double kchi, chi0;                                                       // K (χ̄ − χ0)², χ̄ the mean of the three Wilson angles
  double m1, m2, m3, theta1, theta2, theta3;                               // M1 Δθ_ijk Δθ_kjl + M2 Δθ_ijk Δθ_ijl + M3 Δθ_ijl Δθ_kjl
};

// Stillinger–Weber (LAMMPS pair sw, one element: mW water, Molinero & Moore 2009): for r < aσ
//   φ2 = A ε [B (σ/r)^p − (σ/r)^q] exp(σ / (r − aσ))
//   φ3 = λ ε (cos θ_jik − cos θ0)² exp(γσ / (r_ij − aσ)) exp(γσ / (r_ik − aσ))   (each triplet j-i-k, i the vertex)
struct StillingerWeber {
  bool on = false;
  double eps = 0, sigma = 0, a = 0, lambda = 0, gamma = 0, cos0 = 0, A = 0, B = 0, p = 4, q = 0;
  std::vector<char> atom;   // per atom: takes part
};

struct ForceField {
  std::string name;                        // "GAFF 1.81 (hydrocarbon subset)"
  std::vector<std::string> atom_type;      // per atom: c3, ca, hc, ha
  std::vector<std::string> why;            // per atom: the environment that gave the type
  std::vector<int> type_index;             // per atom: index into type_names / lj
  std::vector<std::string> type_names;
  std::vector<PairType> lj;                // per type
  std::vector<double> charge;              // per atom, e
  std::vector<double> mass;                // per atom, g/mol
  std::vector<BondTerm> bonds;
  std::vector<AngleTerm> angles;
  std::vector<TorsionTerm> dihedrals;      // proper, one entry per Fourier term
  std::vector<TorsionTerm> impropers;      // AMBER order: central atom third
  std::vector<HarmonicTorsion> impropers_harmonic;
  std::vector<InversionTerm> inversions;
  std::vector<BondX> bonds_x;
  std::vector<AngleX> angles_x;
  std::vector<UreyBradley> urey_bradley;
  std::map<std::pair<int, int>, PairFunc> pair_func;       // non-LJ pair forms for type-index pairs (a ≤ b)
  // Separate 1-4 Lennard-Jones parameters per type (CHARMM ε14 / Rmin14, GROMOS C6/C12 1-4); empty: the normal
  // parameters scaled by lj14.
  std::vector<PairType> lj14_types;
  std::vector<Class2Bond> bonds2;
  std::vector<Class2Angle> angles2;
  std::vector<Class2Dihedral> dihedrals2;
  std::vector<Class2Improper> impropers2;
  std::string pair_form = "lj12-6";        // lj12-6: 4ε[(σ/r)¹² − (σ/r)⁶]; lj9-6: ε[2(σ/r)⁹ − 3(σ/r)⁶] (class II)
  std::string mixing = "arithmetic";       // arithmetic, geometric or sixthpower
  std::map<std::pair<int, int>, PairType> pair_override;   // explicit coefficients for type-index pairs (a ≤ b)
  std::vector<std::array<uint32_t, 2>> pairs14;
  std::vector<std::vector<uint32_t>> excluded;   // per atom, sorted: 1-2, 1-3 and 1-4 partners, left out of the pair list
  StillingerWeber sw;                         // many-body term (mW water); counted in the vdW energy
  double lj14 = 0.5, coul14 = 1.0 / 1.2;
  bool keep13 = false;                        // 1-3 pairs interact in full (MARTINI: special_bonds 0 1 1)
  std::set<std::pair<int, int>> excluded_type_pairs;   // type-index pairs (a ≤ b) that never interact (a held graphene sheet)
  // Coarse-grained electrostatics and cut-offs: dielectric (εr) divides every Coulomb term (MARTINI 15, SDK 80);
  // coul_gromacs: GROMACS's force switch from coul_inner to the cut-off instead of DSF / PME (MARTINI); lj_inner starts
  // the lj/gromacs switch; cutoff > 0 is the model's own cut-off (used whatever EnergyOptions says).
  bool coul_gromacs = false;
  double coul_inner = 0, lj_inner = 0, dielectric = 1, cutoff = 0;
  std::vector<std::string> notes;
};

// Types every atom of `s` (C and H only in this slice) and assigns GAFF parameters. Charges are kept from the
// system when it has them, otherwise Gasteiger–Marsili charges are computed. Throws FieldError for other elements
// or for a missing parameter.
ForceField assign_gaff(const System& s);

struct EnergyOptions {
  double cutoff = 10.0;          // Å, LJ and Coulomb
  double skin = 2.0;             // Å, neighbour-list margin (the LAMMPS default for units real)
  bool coulomb = true;           // electrostatics on
  // Long-range electrostatics: damped shifted force (Fennell & Gezelter, J. Chem. Phys. 2006), a cut-off method,
  // or smooth particle-mesh Ewald (Essmann et al. 1995; reciprocal part in Fortran), periodic cells only.
  enum class Electrostatics { DSF, PME } electrostatics = Electrostatics::DSF;
  double dsf_alpha = 0.2;        // Å⁻¹
  double ewald_rtol = 1e-5;      // PME: erfc(β rc) = ewald_rtol sets the Ewald coefficient β (as GROMACS)
  double pme_spacing = 1.0;      // PME: largest grid spacing, Å
  int pme_order = 5;             // PME: B-spline order (5, as LAMMPS pppm)
  double force_cap = 0.0;        // > 0: LJ becomes linear inside the radius where |F| reaches the cap (push-off)
  int threads = 0;               // worker threads for pair terms; 0 = one per hardware thread (at most 16)
  bool tail = true;              // LJ long-range tail corrections to energy and pressure (homogeneous fluid beyond rc)
  // Which terms (r-RESPA splits them): 1 bonded (bonds, angles, torsions, impropers, cross terms), 2 non-bonded (pairs,
  // 1-4 pairs, electrostatics of bonded partners, k-space, self and tail terms), 3 both
  int parts = 3;
};

struct EnergyTerms {
  double bond = 0, angle = 0, dihedral = 0, improper = 0, vdw = 0, coulomb = 0;   // vdw includes the tail term
  double virial = 0;             // Σ r·f, kcal/mol
  double w[6] = {0, 0, 0, 0, 0, 0};   // virial tensor Σ r_a f_b (symmetrised), xx yy zz xy xz yz; trace = virial
  double total() const { return bond + angle + dihedral + improper + vdw + coulomb; }
};

// Energy and forces for a force field over a periodic cell. The neighbour list is rebuilt when an atom has moved
// more than half the skin, or when the cell changes.
class ThreadPool;

class Evaluator {
 public:
  Evaluator(const ForceField& ff, const EnergyOptions& o);
  ~Evaluator();
  Evaluator(const Evaluator&) = delete;
  Evaluator& operator=(const Evaluator&) = delete;
  int threads() const;
  void set_options(const EnergyOptions& o);
  const EnergyOptions& options() const { return opt_; }
  // x: 3N positions (Å). f: 3N forces out (kcal/mol/Å). Returns the energy terms.
  EnergyTerms compute(const std::vector<double>& x, const Cell& cell, std::vector<double>& f);
  int list_builds() const { return builds_; }
  // Freezes the set of interacting pairs at configuration (x, cell): later calls use exactly the pairs inside the
  // cut-off there, whatever their new distance, and never rebuild the list (small deformations only). Makes the energy
  // a smooth function of strain for finite-difference second derivatives (the Born term without the impulsive
  // contribution of pairs crossing a truncated cut-off, as analytic Born terms are defined). unfreeze() ends it.
  void freeze_pairs(const std::vector<double>& x, const Cell& cell);
  void unfreeze() { frozen_ = false; }
  size_t pairs() const { return pi_.size(); }

 private:
  void build(const std::vector<double>& x, const Cell& cell);
  void cap_radii();
  const ForceField& ff_;
  EnergyOptions opt_;
  std::vector<double> x0_;
  Cell cell0_;
  std::vector<uint32_t> pi_, pj_;
  std::vector<double> shift_;        // 3 per pair: Cartesian image shift, d = x_j − x_i + shift
  std::vector<int32_t> ishift_;      // the same shifts as integer lattice coefficients
  Cell cells_;                       // the cell the Cartesian shifts were computed for
  std::vector<double> eps_, s6_, rcap2_, ecap_, fcap_;   // per type pair
  bool lj96_ = false;                // class II 9-6 LJ
  std::vector<uint8_t> form_;        // per type pair: 0 LJ, else the PairFunc form
  std::vector<double> gsw_;          // per type pair, 5 each: lj/gromacs switch coefficients
  std::vector<double> rmin2_, emin_; // per type pair: SDK minimum (r², energy), for the angle's 1-3 repulsion
  double excl_r2_ = 36.0;            // bonded partners closer than this (Å²) are the bonded image
  std::vector<char> skip_type_;      // per type pair: excluded by the force field (never in the pair list)
  std::vector<double> qeff_;         // charges / √εr: what every Coulomb term uses
  std::vector<double> pa_, pb_, pc_;
  std::vector<double> eps14_, s614_;  // separate 1-4 LJ, when the force field has them
  std::vector<double> type_count_;   // atoms per type, for the tail correction
  struct Excluded { uint32_t i, j; double factor; };
  std::vector<Excluded> excl_;       // 1-2, 1-3 (factor 0) and 1-4 (Coulomb factor) partners
  std::unique_ptr<ThreadPool> pool_;
  std::vector<std::vector<double>> tf_;   // per-worker force buffers
  int builds_ = 0;
  bool built_ = false;
  bool frozen_ = false;
  std::vector<uint8_t> inside_;      // per pair: inside the cut-off when frozen
};

// LJ parameters between two type indices, with the force field's mixing rule and explicit pairs.
PairType mixed_pair(const ForceField& ff, int a, int b);

// Pressure (atm) from a virial and volume at zero temperature.
double pressure_atm(double virial, double volume);

}  // namespace caps
