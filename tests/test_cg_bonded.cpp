// Bonded Boltzmann inversion (caps/cg_bonded.hpp): an ideal chain drawn from known bond, angle and dihedral potentials
// (no non-bonded terms, so its internal coordinates are independent and drawn exactly from r² e^(−U/kT), sin θ e^(−U/kT)
// and e^(−U/kT)) gives the potentials back within 0.1 k_B T where they are sampled; walls rise outside; the tables'
// files; weights; a bonded-IBI step that removes a deliberate error.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <random>

#include "caps/cg_bonded.hpp"

using namespace caps;

namespace {
constexpr double kB = 0.0019872067, T = 300, kT = kB * T, kPi = 3.14159265358979323846;

double Ubond(double r) { return 3.0 * (r - 5.0) * (r - 5.0) + 0.6 * (r - 5.0) * (r - 5.0) * (r - 5.0); }   // kcal/mol, anharmonic
double Uangle(double t) { return 0.0012 * (t - 130.0) * (t - 130.0); }                                       // degrees
double Udih(double p) { const double x = p * kPi / 180; return 0.8 * (1 + std::cos(x)) + 0.35 * (1 - std::cos(2 * x)) + 0.15 * std::sin(x); }   // asymmetric: catches a mirrored sign

// a uniform number in [0, 1) the same on every platform (std::uniform_real_distribution differs between libraries)
double unit(std::mt19937_64& g) { return double(g() >> 11) * (1.0 / 9007199254740992.0); }

// samples from a density on [lo, hi] by its tabulated CDF
struct Sampler {
  std::vector<double> x, cdf;
  Sampler(const std::function<double(double)>& w, double lo, double hi, int n = 200000) {
    x.resize(size_t(n)), cdf.resize(size_t(n));
    double c = 0;
    for (int i = 0; i < n; ++i) {
      x[size_t(i)] = lo + (hi - lo) * (i + 0.5) / n;
      c += w(x[size_t(i)]);
      cdf[size_t(i)] = c;
    }
    for (auto& v : cdf) v /= c;
  }
  double operator()(std::mt19937_64& g) const {
    const double u = unit(g);
    const size_t i = size_t(std::lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin());
    const double dx = x[1] - x[0];
    return x[std::min(i, x.size() - 1)] + (unit(g) - 0.5) * dx;
  }
};

Vec3 place(const Vec3& a, const Vec3& b, const Vec3& c, double r, double theta_deg, double phi_deg) {
  const double t = theta_deg * kPi / 180, p = phi_deg * kPi / 180;
  Vec3 bc = c - b;
  bc = bc * (1.0 / norm(bc));
  Vec3 n = cross(b - a, bc);
  n = n * (1.0 / norm(n));
  const Vec3 m = cross(n, bc);
  return c + bc * (-r * std::cos(t)) + m * (r * std::sin(t) * std::cos(p)) + n * (r * std::sin(t) * std::sin(p));
}

struct Ideal {
  CgTopology top;
  std::vector<std::vector<Vec3>> frames;
};

// molecules A-B-A-B: bonds A-B ×3, angles A-B-A and B-A-B, dihedral A-B-A-B; `shift` moves the bond minimum
Ideal ideal(int molecules, int frames, uint64_t seed, double shift = 0) {
  Ideal I;
  for (int m = 0; m < molecules; ++m) {
    const int o = 4 * m;
    for (int k = 0; k < 4; ++k) I.top.kind.push_back(k % 2 ? "B" : "A"), I.top.mol.push_back(m);
    I.top.bonds.push_back({o, o + 1}), I.top.bonds.push_back({o + 1, o + 2}), I.top.bonds.push_back({o + 2, o + 3});
    I.top.angles.push_back({o, o + 1, o + 2}), I.top.angles.push_back({o + 1, o + 2, o + 3});
    I.top.dihedrals.push_back({o, o + 1, o + 2, o + 3});
    I.top.chains.push_back({o, o + 1, o + 2, o + 3});
  }
  const Sampler rb([&](double r) { return r * r * std::exp(-Ubond(r - shift) / kT); }, 2.0, 8.0);
  const Sampler ta([&](double t) { return std::sin(t * kPi / 180) * std::exp(-Uangle(t) / kT); }, 0.0, 180.0);
  const Sampler pd([&](double p) { return std::exp(-Udih(p) / kT); }, -180.0, 180.0);
  std::mt19937_64 g(seed);
  for (int f = 0; f < frames; ++f) {
    std::vector<Vec3> pos;
    for (int m = 0; m < molecules; ++m) {
      const Vec3 a{0, 0, 0};
      const double r1 = rb(g), r2 = rb(g), r3 = rb(g), t1 = ta(g), t2 = ta(g), ph = pd(g);
      const Vec3 b{r1, 0, 0};
      const double t = t1 * kPi / 180;
      const Vec3 c = b + Vec3{-r2 * std::cos(t), r2 * std::sin(t), 0};
      const Vec3 d = place(a, b, c, r3, t2, ph);
      const Vec3 off{double(m % 20) * 30.0, double(m / 20) * 30.0, 0};
      pos.push_back(a + off), pos.push_back(b + off), pos.push_back(c + off), pos.push_back(d + off);
    }
    I.frames.push_back(pos);
  }
  return I;
}

CgTypes types() { return {{"A", "B"}, {"A-B"}, {"A-B-A", "B-A-B"}, {"A-B-A-B"}}; }

// largest |U_table − U_true| (k_B T) where U_true < cut k_B T above its minimum, both shifted to agree there on average
double worst(const CgBondedTable& t, const std::function<double(double)>& truth, double cut, double* where = nullptr) {
  double umin = 1e300;
  for (double x : t.x) umin = std::min(umin, truth(x));
  std::vector<std::pair<double, double>> v;
  std::vector<double> xs;
  for (size_t i = 0; i < t.x.size(); ++i)
    if (truth(t.x[i]) - umin < cut * kT && t.x[i] >= t.lo && t.x[i] <= t.hi) v.push_back({t.U[i], truth(t.x[i])}), xs.push_back(t.x[i]);
  double mean = 0;
  for (const auto& [a, b] : v) mean += a - b;
  mean /= double(v.size());
  double w = 0;
  for (size_t i = 0; i < v.size(); ++i) {
    const double d = std::fabs(v[i].first - v[i].second - mean) / kT;
    if (d > w) { w = d; if (where) *where = xs[i]; }
  }
  return w;
}
}  // namespace

TEST(CgBonded, IdealChainRoundTrip) {
  // the geometry first: the builder's dihedral is the one measured (an asymmetric potential would hide a mirror)
  const Vec3 a{0, 0, 0}, b{1.5, 0, 0}, c{2, 1.4, 0};
  for (double ph : {-150.0, -60.0, 0.0, 45.0, 170.0}) {
    const Vec3 d = place(a, b, c, 1.5, 110, ph);
    EXPECT_NEAR(cg_dihedral_deg(a, b, c, d), ph, 1e-9);
    EXPECT_NEAR(cg_angle_deg(b, c, d), 110, 1e-9);
  }
  const Ideal I = ideal(400, 1000, 11);   // 1.2 million bonds, 400 000 of each angle and the dihedral
  CgBondedOptions o;
  o.temperature = T;
  CgBondedAccumulator acc(types(), o);
  const int s = acc.add_system(I.top, 1, "ideal");
  acc.set_expected_frames(s, I.frames.size());
  for (const auto& f : I.frames) acc.add_frame(s, f, Cell{});
  const CgBondedResult r = invert_bonded(acc);
  ASSERT_EQ(r.bonds.size(), 1u);
  ASSERT_EQ(r.angles.size(), 2u);
  ASSERT_EQ(r.dihedrals.size(), 1u);
  double at = 0;
  EXPECT_LT(worst(r.bonds[0], Ubond, 2.5, &at), 0.1) << "bond, worst at " << at;
  for (const auto& t : r.angles) EXPECT_LT(worst(t, Uangle, 2.5, &at), 0.1) << t.key << " worst at " << at;
  std::printf("round trip (k_B T): bond %.3f · angles %.3f %.3f · dihedral %.3f · bond samples %ld · halves %.3f\n", worst(r.bonds[0], Ubond, 2.5),
              worst(r.angles[0], Uangle, 2.5), worst(r.angles[1], Uangle, 2.5), worst(r.dihedrals[0], Udih, 2.5), r.bonds[0].count, r.bonds[0].half_diff);
  EXPECT_LT(worst(r.dihedrals[0], Udih, 2.5, &at), 0.1) << "dihedral worst at " << at;
  // the well and the curvature of a parabola through it (E = k (r − r0)², k ≈ 3 kcal/mol/Å² near the minimum)
  EXPECT_NEAR(r.bonds[0].x0, 5.0, 0.05);
  EXPECT_NEAR(r.bonds[0].k_harmonic, 3.0, 0.6);
  EXPECT_NEAR(r.angles[0].x0, 130.0, 3.0);
  // walls: outside the sampled range the potential keeps rising, smoothly
  const auto& B = r.bonds[0];
  for (size_t i = 1; i < B.x.size(); ++i) {
    if (B.x[i] < B.lo) EXPECT_LT(B.U[i], B.U[i - 1] + 1e-9) << "rising to the left at " << B.x[i];
    if (B.x[i] > B.hi) EXPECT_GT(B.U[i], B.U[i - 1] - 1e-9) << "rising to the right at " << B.x[i];
    // value and slope continuous: each step matches the force at its ends to first order
    const double dx = B.x[i] - B.x[i - 1];
    EXPECT_NEAR(B.U[i] - B.U[i - 1], -0.5 * (B.F[i] + B.F[i - 1]) * dx, 0.02) << "continuous at " << B.x[i];
  }
  // well sampled: the halves agree
  EXPECT_LT(r.bonds[0].half_diff, 0.2);
  // a dihedral table covers less than a full turn (LAMMPS) and is periodic across ±180
  const auto& D = r.dihedrals[0];
  EXPECT_DOUBLE_EQ(D.x.front(), -180.0);
  EXPECT_LT(D.x.back(), 180.0);
  EXPECT_LT(std::fabs(D.U.front() - (D.U.back() + (D.U.back() - D.U[D.U.size() - 2]))), 0.05);
}

TEST(CgBonded, WeightsFilesAndRefinement) {
  // two systems: the second's bonds 0.3 Å longer; equal weights put the minimum halfway, weight 3 : 1 nearer the first
  const Ideal I1 = ideal(200, 100, 1), I2 = ideal(200, 400, 2, 0.3);   // the second has 4× the frames
  CgBondedOptions o;
  auto x0 = [&](double w1, double w2) {
    CgBondedAccumulator acc(types(), o);
    const int a = acc.add_system(I1.top, w1), b = acc.add_system(I2.top, w2);
    for (const auto& f : I1.frames) acc.add_frame(a, f, Cell{});
    for (const auto& f : I2.frames) acc.add_frame(b, f, Cell{});
    return invert_bonded(acc).bonds[0].x0;
  };
  EXPECT_NEAR(x0(1, 1), 5.15, 0.06) << "frame counts do not weigh: each system is normalised";
  EXPECT_LT(x0(3, 1), x0(1, 1) - 0.03);
  // files
  CgBondedAccumulator acc(types(), o);
  const int s = acc.add_system(I1.top);
  for (const auto& f : I1.frames) acc.add_frame(s, f, Cell{});
  const CgBondedResult r = invert_bonded(acc);
  const auto dir = std::filesystem::temp_directory_path() / "caps_cg_bonded_test";
  const auto files = write_bonded(r, types(), o, dir.string());
  EXPECT_GE(files.size(), 8u);
  std::ifstream in((dir / "bonded.in").string());
  std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_NE(text.find("bond_coeff     1 ${BONDED}/bonds.table A-B"), std::string::npos) << text;
  EXPECT_NE(text.find("angle_coeff    2 ${BONDED}/angles.table B-A-B"), std::string::npos);
  EXPECT_NE(text.find("dihedral_coeff 1 aat 1 170 180 ${BONDED}/dihedrals.table A-B-A-B"), std::string::npos);
  EXPECT_TRUE(std::filesystem::exists(dir / "gromacs" / "table_d0.xvg"));
  std::ifstream bj((dir / "bonded.json").string());
  std::string js((std::istreambuf_iterator<char>(bj)), std::istreambuf_iterator<char>());
  const CgBondedResult back = bonded_from_json(Json::parse(js));
  ASSERT_EQ(back.bonds.size(), 1u);
  EXPECT_EQ(back.bonds[0].U.size(), r.bonds[0].U.size());
  // refinement: tables from the shifted system, "CG run" = the true one; one full step (α = 1) moves the bond well back
  CgBondedAccumulator wrong(types(), o);
  const int w = wrong.add_system(I2.top);
  for (const auto& f : I2.frames) wrong.add_frame(w, f, Cell{});
  CgBondedResult target = invert_bonded(acc);       // what we want: I1's distributions
  CgBondedResult current = invert_bonded(wrong);   // tables that give I2's distributions (as if from a CG run)
  // a step U ← U + kT ln(P_run / P_target) with the run = I2 and the target = I1 applied to I2's own potential gives I1's
  CgBondedResult step = current;
  for (size_t k = 0; k < step.bonds.size(); ++k) step.bonds[k].P = target.bonds[k].P, step.bonds[k].Uinv = current.bonds[k].Uinv;
  std::map<std::string, double> res;
  const CgBondedResult fixed = refine_bonded(step, wrong, 1.0, &res);
  EXPECT_NEAR(fixed.bonds[0].x0, 5.0, 0.06);
  EXPECT_GT(res["bond A-B"], 0.1);
  std::filesystem::remove_all(dir);
}
