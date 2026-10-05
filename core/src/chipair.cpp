// χ from pair contacts (Fan, Olafson, Blanco & Hsu 1992); see chipair.hpp.
#include "caps/rng.hpp"
#include "caps/chipair.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <random>
#include <thread>

#include "caps/elements.hpp"
#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/typing.hpp"
#include "caps/molecule.hpp"
#include "caps/uff.hpp"

namespace caps {

namespace {

constexpr double kR = 0.0019872043;       // kcal/mol/K
constexpr double kCoul = 332.06371;       // kcal·Å/(mol·e²)
constexpr int kChunk = 256;               // samples per random stream: the same numbers for any thread count

uint64_t mix(uint64_t x) {   // splitmix64
  x += 0x9E3779B97F4A7C15ull;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
  return x ^ (x >> 31);
}

// A rigid molecule about its centroid, with per-atom Lennard-Jones type, charge and contact radius.
struct Mol {
  std::string name;
  std::vector<Vec3> x;
  std::vector<int> type;
  std::vector<double> q, rad;
  double reach = 0;   // largest centroid distance + radius
};

// Pair tables between two molecules.
struct Pair {
  size_t nx = 0, ny = 0;
  std::vector<double> eps4, s6, qq, sig2;   // 4ε, σ⁶, k q_i q_j, (contact σ)²
};

Pair make_pair_table(const Mol& X, const Mol& Y, const ForceField& ffc, int y_type_offset, double scale) {
  Pair p;
  p.nx = X.x.size(), p.ny = Y.x.size();
  const size_t n = p.nx * p.ny;
  p.eps4.resize(n), p.s6.resize(n), p.qq.resize(n), p.sig2.resize(n);
  for (size_t i = 0; i < p.nx; ++i)
    for (size_t j = 0; j < p.ny; ++j) {
      const size_t k = i * p.ny + j;
      const PairType pt = mixed_pair(ffc, X.type[i], Y.type[j] + y_type_offset);
      p.eps4[k] = 4 * pt.eps;
      p.s6[k] = std::pow(pt.sigma, 6);
      p.qq[k] = kCoul * X.q[i] * Y.q[j];
      const double s = scale * (X.rad[i] + Y.rad[j]);
      p.sig2[k] = s * s;
    }
  return p;
}

Vec3 rotate(const double q[4], const Vec3& v) {   // unit quaternion (w, x, y, z)
  const double w = q[0], x = q[1], y = q[2], z = q[3];
  return {(1 - 2 * (y * y + z * z)) * v[0] + 2 * (x * y - w * z) * v[1] + 2 * (x * z + w * y) * v[2],
          2 * (x * y + w * z) * v[0] + (1 - 2 * (x * x + z * z)) * v[1] + 2 * (y * z - w * x) * v[2],
          2 * (x * z - w * y) * v[0] + 2 * (y * z + w * x) * v[1] + (1 - 2 * (x * x + y * y)) * v[2]};
}

struct Rng {
  std::mt19937_64 g;
  caps::UniformReal<double> u{0.0, 1.0};
  explicit Rng(uint64_t s) : g(s) {}
  double operator()() { return u(g); }
  void quaternion(double q[4]) {   // uniform random rotation (Shoemake 1992)
    const double u1 = (*this)(), u2 = 2 * M_PI * (*this)(), u3 = 2 * M_PI * (*this)();
    const double a = std::sqrt(1 - u1), b = std::sqrt(u1);
    q[0] = a * std::sin(u2), q[1] = a * std::cos(u2), q[2] = b * std::sin(u3), q[3] = b * std::cos(u3);
  }
  Vec3 direction() {
    const double z = 2 * (*this)() - 1, phi = 2 * M_PI * (*this)(), r = std::sqrt(std::max(0.0, 1 - z * z));
    return {r * std::cos(phi), r * std::sin(phi), z};
  }
};

// Y rotated, then slid in along u from far away until the first atom pair touches (the largest root of
// |y_j + d u − x_i| = σ_ij over all pairs). Returns d and fills y (placed coordinates).
double contact(const Mol& X, const Mol& Y, const Pair& P, const double q[4], const Vec3& u, std::vector<Vec3>& y) {
  y.resize(Y.x.size());
  for (size_t j = 0; j < y.size(); ++j) y[j] = rotate(q, Y.x[j]);
  double d = -1e30;
  for (size_t i = 0; i < P.nx; ++i)
    for (size_t j = 0; j < P.ny; ++j) {
      const Vec3 D = y[j] - X.x[i];
      const double b = dot(D, u), disc = b * b - dot(D, D) + P.sig2[i * P.ny + j];
      if (disc >= 0) d = std::max(d, -b + std::sqrt(disc));
    }
  for (auto& v : y) v = v + u * d;
  return d;
}

double pair_energy(const Mol& X, const std::vector<Vec3>& y, const Pair& P) {
  double e = 0;
  for (size_t i = 0; i < P.nx; ++i)
    for (size_t j = 0; j < P.ny; ++j) {
      const size_t k = i * P.ny + j;
      const Vec3 D = y[j] - X.x[i];
      const double r2 = dot(D, D), sr = P.s6[k] / (r2 * r2 * r2);
      e += P.eps4[k] * (sr * sr - sr) + P.qq[k] / std::sqrt(r2);
    }
  return e;
}

template <class F>
void parallel(int n, int threads, F&& f) {
  const int nt = std::max(1, std::min(threads, n));
  std::atomic<int> next{0};
  std::vector<std::thread> pool;
  for (int t = 0; t < nt; ++t)
    pool.emplace_back([&] {
      for (int k; (k = next.fetch_add(1)) < n;) f(k);
    });
  for (auto& th : pool) th.join();
}

System as_system(const std::string& smiles) {
  BuildOptions b;
  b.forcefield = "uff";
  b.conformers = 4;
  return build_molecule(smiles, b).system;
}

}  // namespace

std::string capped_unit(const std::string& unit) {
  // a bare * or a bracket atom whose element is * ([*], [*:1], [1*]) becomes [H]
  std::string r;
  for (size_t i = 0; i < unit.size(); ++i) {
    if (unit[i] == '*') { r += "[H]"; continue; }
    if (unit[i] == '[') {
      const size_t e = unit.find(']', i);
      if (e != std::string::npos) {
        const std::string in = unit.substr(i + 1, e - i - 1);
        const size_t k = in.find_first_not_of("0123456789");
        if (k != std::string::npos && in[k] == '*') { r += "[H]"; i = e; continue; }
      }
    }
    r += unit[i];
  }
  return r;
}

ChiPairResult chi_by_contacts(const ChiPairOptions& o) {
  if (o.a_smiles.empty() || o.b_smiles.empty()) throw std::invalid_argument("χ by pair contacts needs two molecules (SMILES)");
  if (o.temperatures.empty()) throw std::invalid_argument("χ by pair contacts needs at least one temperature");
  auto tell = [&](const std::string& stage, double f) {
    if (o.progress && !o.progress(stage, f)) throw ChiPairCancelled();
  };
  const int threads = o.threads > 0 ? o.threads : std::max(1, int(std::thread::hardware_concurrency()));
  ChiPairResult r;
  tell("build", 0);
  const std::string sa = capped_unit(o.a_smiles), sb = capped_unit(o.b_smiles);
  const System A0 = as_system(sa), B0 = as_system(sb);

  // one force-field family for both: GAFF when it types both (C and H), else UFF with Gasteiger / QEq charges
  auto gaff_types = [](const System& s) {
    if (!std::all_of(s.atoms.begin(), s.atoms.end(), [](const Atom& a) { return a.element == 1 || a.element == 6; })) return false;
    try { assign_gaff(s); return true; } catch (const FieldError&) { return false; }
  };
  ForceField fa, fb;
  if (!o.forcefield.empty()) {
    const FFDef def = load_forcefield(o.forcefield);
    auto typed = [&](const System& s, const std::string& which) {
      const TypingResult tr = assign_types(s, def);
      if (tr.untyped) throw FieldError(which + ": " + std::to_string(tr.untyped) + " atoms match no typing rule of " + def.name);
      ParamReport rep;
      return parameterize(s, def, tr.types, "gasteiger", &rep, true);   // bonded gaps do not matter here
    };
    fa = typed(A0, sa), fb = typed(B0, sb);
    r.forcefield = def.name + ", Gasteiger charges";
  } else {
    const bool gaff = gaff_types(A0) && gaff_types(B0);
    std::string qa, qb;
    fa = gaff ? assign_gaff(A0) : uff_with_charges(A0, &qa);
    fb = gaff ? assign_gaff(B0) : uff_with_charges(B0, &qb);
    r.forcefield = gaff ? "GAFF (built-in C and H subset), Gasteiger charges" : "UFF, " + (qa == qb ? qa : qa + " / " + qb) + " charges";
  }
  ForceField ffc;   // both molecules' types in one table for the mixing rule
  ffc.mixing = fa.mixing;
  ffc.lj = fa.lj;
  ffc.lj.insert(ffc.lj.end(), fb.lj.begin(), fb.lj.end());
  const int off = int(fa.lj.size());

  auto fill = [](const System& s, const ForceField& ff, const std::string& smiles) {
    Mol m;
    Vec3 c{0, 0, 0};
    for (const auto& a : s.atoms) c = c + a.pos;
    c = c * (1.0 / double(s.atoms.size()));
    m.name = smiles;
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      m.x.push_back(s.atoms[i].pos - c);
      m.type.push_back(ff.type_index[i]);
      m.q.push_back(ff.charge[i]);
      m.rad.push_back(element(s.atoms[i].element).vdw);
      m.reach = std::max(m.reach, norm(m.x.back()) + m.rad.back());
    }
    return m;
  };
  const Mol A = fill(A0, fa, sa), B = fill(B0, fb, sb);
  r.a_name = sa, r.b_name = sb;
  r.a_atoms = int(A.x.size()), r.b_atoms = int(B.x.size());

  const size_t nT = o.temperatures.size();
  r.temperatures = o.temperatures;

  // energies of pairs in contact, Boltzmann averages per temperature with block errors
  auto sample = [&](const Mol& X, const Mol& Y, int xoff, int yoff, int kind, ChiPairKind& out) {
    Mol Xs = X;   // type indices of X and Y into ffc
    for (auto& t : Xs.type) t += xoff;
    const Mol& Ys = Y;
    const Pair P = make_pair_table(Xs, Ys, ffc, yoff, o.contact_scale);
    std::vector<double> E(size_t(std::max(10, o.samples)));
    const int chunks = int((E.size() + kChunk - 1) / kChunk);
    parallel(chunks, threads, [&](int c) {
      Rng g(mix(o.seed * 1000003ull + uint64_t(kind) * 7919ull + uint64_t(c)));
      std::vector<Vec3> y;
      double q[4];
      for (size_t s = size_t(c) * kChunk; s < std::min(E.size(), size_t(c + 1) * kChunk); ++s) {
        g.quaternion(q);
        const Vec3 u = g.direction();
        contact(Xs, Ys, P, q, u, y);
        E[s] = pair_energy(Xs, y, P);
      }
    });
    out.e_min = *std::min_element(E.begin(), E.end());
    double m = 0;
    for (double e : E) m += e;
    out.e_mean = m / double(E.size());
    const int nb = 10;
    const size_t per = E.size() / nb;
    for (double T : o.temperatures) {
      const double beta = 1 / (kR * T);
      auto avg = [&](size_t a, size_t b) {
        double w = 0, we = 0, lo = 1e300;
        for (size_t k = a; k < b; ++k) lo = std::min(lo, E[k]);
        for (size_t k = a; k < b; ++k) {
          const double x = std::exp(-beta * (E[k] - lo));
          w += x, we += x * E[k];
        }
        return we / w;
      };
      out.e_t.push_back(avg(0, E.size()));
      std::vector<double> blk;
      for (int k = 0; k < nb; ++k) blk.push_back(avg(size_t(k) * per, size_t(k + 1) * per));
      double bm = 0, s2 = 0;
      for (double x : blk) bm += x;
      bm /= nb;
      for (double x : blk) s2 += (x - bm) * (x - bm);
      out.e_t_error.push_back(std::sqrt(s2 / (nb - 1) / nb));
    }
    // histogram over the lowest 99 %
    std::vector<double> sorted = E;
    std::sort(sorted.begin(), sorted.end());
    const double lo = sorted.front(), hi = sorted[size_t(0.99 * double(sorted.size() - 1))];
    const int bins = 40;
    const double w = (hi - lo) / bins;
    if (w > 0) {
      out.hist_e.resize(bins), out.hist_p.assign(bins, 0);
      for (int k = 0; k < bins; ++k) out.hist_e[size_t(k)] = lo + (k + 0.5) * w;
      for (double e : E)
        if (e <= hi) out.hist_p[size_t(std::min(bins - 1, int((e - lo) / w)))] += 1.0 / double(E.size());
    }
  };

  // coordination: j molecules packed around one i at random contacts, none overlapping another
  auto coordination = [&](const Mol& X, const Mol& Y, int xoff, int yoff, int kind, ChiPairKind& out) {
    Mol Xs = X, Ys = Y;
    for (auto& t : Xs.type) t += xoff;
    const Pair P = make_pair_table(Xs, Ys, ffc, yoff, o.contact_scale);
    std::vector<double> sig2yy(Ys.x.size() * Ys.x.size());
    for (size_t i = 0; i < Ys.x.size(); ++i)
      for (size_t j = 0; j < Ys.x.size(); ++j) {
        const double s = o.contact_scale * (Ys.rad[i] + Ys.rad[j]);
        sig2yy[i * Ys.x.size() + j] = s * s;
      }
    std::vector<double> Z(size_t(std::max(2, o.pack_trials)));
    parallel(int(Z.size()), threads, [&](int t) {
      Rng g(mix(o.seed * 2000003ull + uint64_t(kind) * 104729ull + uint64_t(t)));
      std::vector<std::vector<Vec3>> placed;
      std::vector<Vec3> centres, y;
      double q[4];
      for (int miss = 0; miss < o.pack_misses;) {
        g.quaternion(q);
        const Vec3 u = g.direction();
        const double d = contact(Xs, Ys, P, q, u, y);
        const Vec3 cy = u * d;
        bool clash = false;
        for (size_t p = 0; p < placed.size() && !clash; ++p) {
          if (norm(centres[p] - cy) > 2 * Ys.reach) continue;
          for (size_t i = 0; i < y.size() && !clash; ++i)
            for (size_t j = 0; j < y.size(); ++j) {
              const Vec3 D = placed[p][j] - y[i];
              if (dot(D, D) < sig2yy[i * y.size() + j]) { clash = true; break; }
            }
        }
        if (clash) { ++miss; continue; }
        placed.push_back(y);
        centres.push_back(cy);
        miss = 0;
      }
      Z[size_t(t)] = double(placed.size());
    });
    double m = 0, s2 = 0;
    for (double z : Z) m += z;
    m /= double(Z.size());
    for (double z : Z) s2 += (z - m) * (z - m);
    out.z = m;
    out.z_error = std::sqrt(s2 / double(Z.size() - 1) / double(Z.size()));
  };

  r.aa.name = "A–A", r.ab.name = "A–B", r.ba.name = "B–A", r.bb.name = "B–B";
  tell("pair energies A–A", 0.05);
  sample(A, A, 0, 0, 1, r.aa);
  tell("pair energies A–B", 0.2);
  sample(A, B, 0, off, 2, r.ab);
  tell("pair energies B–B", 0.35);
  sample(B, B, off, off, 3, r.bb);
  r.ba.e_t = r.ab.e_t, r.ba.e_t_error = r.ab.e_t_error, r.ba.e_min = r.ab.e_min, r.ba.e_mean = r.ab.e_mean;   // the same pairs
  tell("coordination A–A", 0.5);
  coordination(A, A, 0, 0, 1, r.aa);
  tell("coordination A–B", 0.62);
  coordination(A, B, 0, off, 2, r.ab);
  tell("coordination B–A", 0.75);
  coordination(B, A, off, 0, 3, r.ba);
  tell("coordination B–B", 0.87);
  coordination(B, B, off, off, 4, r.bb);

  // χ(T) and its error (E_AB and E_BA are one set of samples: their errors add coherently)
  for (size_t k = 0; k < nT; ++k) {
    const double T = o.temperatures[k], RT = kR * T;
    const double eab = r.ab.e_t[k], eaa = r.aa.e_t[k], ebb = r.bb.e_t[k];
    const double de = 0.5 * (r.ab.z * eab + r.ba.z * eab - r.aa.z * eaa - r.bb.z * ebb);
    const double v = std::pow(0.5 * (r.ab.z + r.ba.z) * r.ab.e_t_error[k], 2) + std::pow(0.5 * r.aa.z * r.aa.e_t_error[k], 2) +
                     std::pow(0.5 * r.bb.z * r.bb.e_t_error[k], 2) + std::pow(0.5 * eab * r.ab.z_error, 2) + std::pow(0.5 * eab * r.ba.z_error, 2) +
                     std::pow(0.5 * eaa * r.aa.z_error, 2) + std::pow(0.5 * ebb * r.bb.z_error, 2);
    r.chi.push_back(de / RT);
    r.chi_error.push_back(std::sqrt(v) / RT);
  }
  // χ = A + B/T, least squares in 1/T
  if (nT >= 2) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t k = 0; k < nT; ++k) {
      const double x = 1 / o.temperatures[k];
      sx += x, sy += r.chi[k], sxx += x * x, sxy += x * r.chi[k];
    }
    const double n = double(nT), den = n * sxx - sx * sx;
    r.fit_b = den != 0 ? (n * sxy - sx * sy) / den : 0;
    r.fit_a = (sy - r.fit_b * sx) / n;
  } else {
    r.fit_a = 0, r.fit_b = r.chi[0] * o.temperatures[0];
  }
  // at the reporting temperature: from the fit (what the phase diagram uses), with the error of the nearest temperature
  {
    r.chi_report = r.fit_a + r.fit_b / o.report_temperature;
    size_t near = 0;
    for (size_t k = 1; k < nT; ++k)
      if (std::fabs(o.temperatures[k] - o.report_temperature) < std::fabs(o.temperatures[near] - o.report_temperature)) near = k;
    r.chi_report_error = r.chi_error[near];
  }
  char b[400];
  std::snprintf(b, sizeof b, "χ(%.0f K) = %.3f ± %.3f; χ(T) = %.4f + %.2f/T over %.0f–%.0f K", o.report_temperature, r.chi_report, r.chi_report_error,
                r.fit_a, r.fit_b, *std::min_element(o.temperatures.begin(), o.temperatures.end()), *std::max_element(o.temperatures.begin(), o.temperatures.end()));
  r.notes.push_back(b);
  std::snprintf(b, sizeof b, "Z: A–A %.2f, A–B %.2f, B–A %.2f, B–B %.2f (± %.2f); ⟨E⟩ at %.0f K: A–A %.3f, A–B %.3f, B–B %.3f kcal/mol", r.aa.z, r.ab.z, r.ba.z, r.bb.z,
                std::max({r.aa.z_error, r.ab.z_error, r.ba.z_error, r.bb.z_error}), o.temperatures[nT / 2], r.aa.e_t[nT / 2], r.ab.e_t[nT / 2], r.bb.e_t[nT / 2]);
  r.notes.push_back(b);
  std::snprintf(b, sizeof b, "%d contacts per pair of kinds, %d packed shells per coordination number; %s", o.samples, o.pack_trials, r.forcefield.c_str());
  r.notes.push_back(b);
  r.notes.push_back("rigid molecules (the lowest of 4 UFF conformers), repeat units capped with H, one molecule per lattice site: a screen, "
                    "not a measured χ (Fan et al., Macromolecules 25, 3667 (1992))");
  return r;
}

}  // namespace caps
