// CAPS adsorption locator (see adsorption.hpp).
#include "caps/rng.hpp"
#include "caps/adsorption.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/pairmodel.hpp"

namespace caps {

namespace {

constexpr double kB = 0.0019872043;   // kcal/(mol·K)
constexpr double kPi = 3.14159265358979323846;

Vec3 rotate(const Vec3& v, const Vec3& u, double t) {   // v about unit axis u by t (Rodrigues)
  const double c = std::cos(t), s = std::sin(t);
  return v * c + cross(u, v) * s + u * (dot(u, v) * (1 - c));
}

// The substrate (fixed atoms) and the mobile molecules, with the energies the locator counts.
struct Model {
  const System& s;
  const ForceField& ff;
  PairModel pm;
  const Cell& cell;
  bool per = false, ortho = false;
  double L[3] = {0, 0, 0};
  bool pdim[3] = {false, false, false};
  std::vector<std::vector<uint32_t>> mol;   // atoms of each mobile molecule
  std::vector<uint32_t> fixed;
  // bins of the fixed atoms, at least rc/2 wide, two either side (fewer than five: every bin once)
  int nb[3] = {1, 1, 1};
  std::vector<std::vector<uint32_t>> bins;
  std::vector<std::vector<int>> stencil;
  // A cut-off wider than half the cell's narrowest width reaches more than one image of an atom, so the minimum image
  // above is not enough: every periodic image within the cut-off counts instead. The substrate's images are ghost
  // copies binned in space (a mobile atom folded into the cell meets them at plain distances); adsorbate pairs, a
  // molecule with its own images included, are lattice sums.
  bool images = false;
  double w[3] = {0, 0, 0};           // the cell's widths (Å) across each pair of faces
  double wmin = 0;                   // the narrowest periodic one
  Vec3 rec[3];                       // fractional coordinate of a separation d along axis k: dot(rec[k], d)
  std::vector<Vec3> gpos;            // the substrate's images within the cut-off of the cell
  std::vector<uint32_t> gatom;       // ... and the atom each one copies
  Vec3 glo{0, 0, 0};
  double gw[3] = {1, 1, 1};
  int gn[3] = {1, 1, 1};
  std::vector<std::vector<uint32_t>> gbins;

  Model(const System& sys, const ForceField& f, const std::vector<int>& mobile, double rc, bool coul, double alpha)
      : s(sys), ff(f), pm(f, rc, coul, alpha), cell(sys.cell) {
    per = cell.valid();
    ortho = per && std::fabs(cell.a[1]) + std::fabs(cell.a[2]) + std::fabs(cell.b[0]) + std::fabs(cell.b[2]) + std::fabs(cell.c[0]) + std::fabs(cell.c[1]) < 1e-9;
    L[0] = cell.a[0], L[1] = cell.b[1], L[2] = cell.c[2];
    for (int d = 0; d < 3; ++d) pdim[d] = cell.periodic[size_t(d)];
    int nm = 0;
    const auto comp = sys.molecules(&nm);
    std::map<int, size_t> slot;
    for (size_t k = 0; k < mobile.size(); ++k) slot[mobile[k]] = k;
    mol.resize(mobile.size());
    for (uint32_t i = 0; i < sys.atoms.size(); ++i) {
      const auto it = slot.find(comp[i]);
      if (it == slot.end()) fixed.push_back(i);
      else mol[it->second].push_back(i);
    }
    if (per) {
      const double V = cell.volume(), det = dot(cell.a, cross(cell.b, cell.c));
      w[0] = V / norm(cross(cell.b, cell.c)), w[1] = V / norm(cross(cell.c, cell.a)), w[2] = V / norm(cross(cell.a, cell.b));
      rec[0] = cross(cell.b, cell.c) * (1.0 / det), rec[1] = cross(cell.c, cell.a) * (1.0 / det), rec[2] = cross(cell.a, cell.b) * (1.0 / det);
      wmin = 1e30;
      for (int d = 0; d < 3; ++d)
        if (pdim[d]) nb[d] = std::max(1, int(std::floor(2 * w[d] / rc))), wmin = std::min(wmin, w[d]);
      images = rc > 0.5 * wmin;
    }
    if (images) build_images(rc);
    bins.resize(size_t(nb[0] * nb[1] * nb[2]));
    for (uint32_t i : fixed) bins[size_t(bin_of(sys.atoms[i].pos))].push_back(i);
    stencil.resize(bins.size());
    for (int a = 0; a < nb[0]; ++a)
      for (int b = 0; b < nb[1]; ++b)
        for (int c = 0; c < nb[2]; ++c) {
          auto offs = [&](int d, int v) {
            std::vector<int> r;
            if (nb[d] < 5) for (int k = 0; k < nb[d]; ++k) r.push_back(k);
            else for (int dv = -2; dv <= 2; ++dv) r.push_back(((v + dv) % nb[d] + nb[d]) % nb[d]);
            return r;
          };
          for (int i : offs(0, a))
            for (int j : offs(1, b))
              for (int k : offs(2, c)) stencil[size_t((a * nb[1] + b) * nb[2] + c)].push_back((i * nb[1] + j) * nb[2] + k);
        }
  }
  void build_images(double rc) {
    // a bounded budget: a cut-off many times the cell's width would need millions of images
    double per_pair = 1;
    for (int d = 0; d < 3; ++d) if (pdim[d]) per_pair *= 1 + 2 * rc / w[d];
    const double est = per_pair * double(fixed.size());
    if (est > 4e6 || per_pair > 2e5) {
      char b[320];
      std::snprintf(b, sizeof b, "adsorption: a %.1f Å cut-off in a cell %.2f Å across needs about %.0f periodic images; use a shorter cut-off or a larger cell",
                    rc, wmin, std::max(est, per_pair));
      throw std::invalid_argument(b);
    }
    int n[3];
    double reach[3];
    for (int d = 0; d < 3; ++d) n[d] = pdim[d] ? int(std::ceil(rc / w[d])) : 0, reach[d] = pdim[d] ? rc / w[d] : 0;
    for (uint32_t b : fixed) {
      Vec3 f = cell.to_fractional(s.atoms[b].pos);
      for (int d = 0; d < 3; ++d) if (pdim[d]) f[size_t(d)] -= std::floor(f[size_t(d)]);
      for (int i = -n[0]; i <= n[0]; ++i)
        for (int j = -n[1]; j <= n[1]; ++j)
          for (int k = -n[2]; k <= n[2]; ++k) {
            const Vec3 g{f[0] + i, f[1] + j, f[2] + k};
            bool keep = true;   // within the cut-off of the cell (a point in it is at most rc / w away along each axis)
            for (int d = 0; d < 3 && keep; ++d)
              if (pdim[d]) keep = g[size_t(d)] > -reach[d] && g[size_t(d)] < 1 + reach[d];
            if (keep) gpos.push_back(cell.to_cartesian(g)), gatom.push_back(b);
          }
    }
    // bins at least rc wide over the images' box; a point outside the box searches the bins at its edge
    Vec3 hi{-1e30, -1e30, -1e30};
    glo = {1e30, 1e30, 1e30};
    for (const auto& p : gpos)
      for (size_t d = 0; d < 3; ++d) glo[d] = std::min(glo[d], p[d]), hi[d] = std::max(hi[d], p[d]);
    if (gpos.empty()) glo = {0, 0, 0};   // no substrate (a cluster search): one empty bin, and no 1e30 origin to overflow gbin
    for (int d = 0; d < 3; ++d) {
      const double ext = gpos.empty() ? 0 : hi[size_t(d)] - glo[size_t(d)];
      gn[d] = std::clamp(int(std::floor(ext / rc)), 1, 64);
      gw[d] = std::max(rc, ext / gn[d]);
    }
    gbins.assign(size_t(gn[0] * gn[1] * gn[2]), {});
    for (uint32_t g = 0; g < gpos.size(); ++g) {
      int id[3];
      gbin(gpos[g], id);
      gbins[size_t((id[0] * gn[1] + id[1]) * gn[2] + id[2])].push_back(g);
    }
  }
  void gbin(const Vec3& p, int id[3]) const {
    for (int d = 0; d < 3; ++d) id[d] = std::clamp(int(std::floor((p[size_t(d)] - glo[size_t(d)]) / gw[d])), 0, gn[d] - 1);
  }
  // pm.energy of atoms ia and ib over the lattice translations T with |d + T| inside the cut-off; skip0: T = 0 left out
  double lattice_sum(uint32_t ia, uint32_t ib, const Vec3& d, bool skip0) const {
    const double rc = pm.cutoff();
    int lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    for (int t = 0; t < 3; ++t)
      if (pdim[t]) {
        const double f = dot(rec[t], d);   // the separation's height across faces t, in widths: (f + i) w within ±rc
        lo[t] = int(std::ceil(-rc / w[t] - f)), hi[t] = int(std::floor(rc / w[t] - f));
      }
    double e = 0;
    for (int i = lo[0]; i <= hi[0]; ++i)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int k = lo[2]; k <= hi[2]; ++k) {
          if (skip0 && i == 0 && j == 0 && k == 0) continue;
          const Vec3 r = d + cell.a * double(i) + cell.b * double(j) + cell.c * double(k);
          e += pm.energy(ia, ib, dot(r, r));
        }
    return e;
  }
  Vec3 mi(Vec3 d) const {
    if (!per) return d;
    if (!ortho) return cell.minimum_image(d);
    for (int k = 0; k < 3; ++k)
      if (pdim[k]) d[size_t(k)] -= L[k] * std::round(d[size_t(k)] / L[k]);
    return d;
  }
  int bin_of(const Vec3& p) const {
    if (!per) return 0;
    const Vec3 f = cell.to_fractional(p);
    int id[3];
    for (int d = 0; d < 3; ++d) id[d] = std::min(nb[d] - 1, int((f[size_t(d)] - std::floor(f[size_t(d)])) * nb[d]));
    return (id[0] * nb[1] + id[1]) * nb[2] + id[2];
  }
  // molecule m at positions P (one per atom of m) with the substrate
  double with_fixed(size_t m, const std::vector<Vec3>& P) const {
    double e = 0;
    if (images) {
      for (size_t k = 0; k < mol[m].size(); ++k) {
        const Vec3 q = cell.wrap(P[k]);
        int id[3];
        gbin(q, id);
        for (int a = std::max(0, id[0] - 1); a <= std::min(gn[0] - 1, id[0] + 1); ++a)
          for (int b = std::max(0, id[1] - 1); b <= std::min(gn[1] - 1, id[1] + 1); ++b)
            for (int c = std::max(0, id[2] - 1); c <= std::min(gn[2] - 1, id[2] + 1); ++c)
              for (uint32_t g : gbins[size_t((a * gn[1] + b) * gn[2] + c)]) {
                const Vec3 d = gpos[g] - q;
                e += pm.energy(mol[m][k], gatom[g], dot(d, d));
              }
        if (!std::isfinite(e)) return e;
      }
      return e;
    }
    for (size_t k = 0; k < mol[m].size(); ++k) {
      const uint32_t a = mol[m][k];
      for (int bn : stencil[size_t(bin_of(P[k]))])
        for (uint32_t b : bins[size_t(bn)]) {
          const Vec3 d = mi(s.atoms[b].pos - P[k]);
          e += pm.energy(a, b, dot(d, d));
        }
      if (!std::isfinite(e)) return e;
    }
    return e;
  }
  double between(size_t m, const std::vector<Vec3>& P, size_t n, const std::vector<Vec3>& Q) const {
    double e = 0;
    for (size_t k = 0; k < mol[m].size(); ++k)
      for (size_t q = 0; q < mol[n].size(); ++q) {
        if (images) {
          e += lattice_sum(mol[m][k], mol[n][q], Q[q] - P[k], false);
          continue;
        }
        const Vec3 d = mi(Q[q] - P[k]);
        e += pm.energy(mol[m][k], mol[n][q], dot(d, d));
      }
    return e;
  }
  // molecule m at P with its own periodic images (each pair of images once); 0 while one image is all the cut-off reaches
  double self(size_t m, const std::vector<Vec3>& P) const {
    if (!images) return 0;
    double e = 0;
    for (size_t k = 0; k < mol[m].size(); ++k)
      for (size_t q = 0; q < mol[m].size(); ++q) e += lattice_sum(mol[m][k], mol[m][q], P[q] - P[k], true);
    return 0.5 * e;
  }
};

std::string formula(const System& s, const std::vector<uint32_t>& atoms) {   // Hill order
  std::map<std::string, int> n;
  const bool real = std::any_of(atoms.begin(), atoms.end(), [&](uint32_t a) { return s.atoms[a].element > 0; });
  for (uint32_t a : atoms)
    if (!real || s.atoms[a].element > 0) ++n[element(s.atoms[a].element).symbol];   // a model's massless sites are not in the formula
  std::string f;
  auto put = [&](const std::string& el) {
    const auto it = n.find(el);
    if (it == n.end()) return;
    f += el + (it->second > 1 ? std::to_string(it->second) : "");
    n.erase(it);
  };
  if (n.count("C")) { put("C"); put("H"); }
  for (const auto& [el, c] : std::map<std::string, int>(n)) put(el);
  return f;
}

std::vector<int> default_mobile(const System& s) {   // every molecule but the largest
  int nm = 0;
  const auto comp = s.molecules(&nm);
  std::vector<int> size(size_t(nm), 0);
  for (int c : comp) ++size[size_t(c)];
  const int big = int(std::max_element(size.begin(), size.end()) - size.begin());
  std::vector<int> out;
  for (int m = 0; m < nm; ++m) if (m != big) out.push_back(m);
  return out;
}

}  // namespace

std::vector<int> adsorbate_molecules(const System& s, int first_atom) {
  if (first_atom < 0) return default_mobile(s);
  int nm = 0;
  const auto comp = s.molecules(&nm);
  std::set<int> m;
  for (size_t i = size_t(first_atom); i < s.atoms.size(); ++i) m.insert(comp[i]);
  return {m.begin(), m.end()};
}

double adsorption_interaction(const System& s, const ForceField& ff, const std::vector<int>& mobile_in, double cutoff, bool coulomb, double dsf_alpha) {
  const auto mobile = mobile_in.empty() ? default_mobile(s) : mobile_in;
  Model M(s, ff, mobile, cutoff, coulomb, dsf_alpha);
  std::vector<std::vector<Vec3>> P(M.mol.size());
  for (size_t m = 0; m < M.mol.size(); ++m) for (uint32_t a : M.mol[m]) P[m].push_back(s.atoms[a].pos);
  double e = 0;
  for (size_t m = 0; m < M.mol.size(); ++m) {
    e += M.with_fixed(m, P[m]) + M.self(m, P[m]);
    for (size_t n = m + 1; n < M.mol.size(); ++n) e += M.between(m, P[m], n, P[n]);
  }
  return e;
}

void locate_adsorption(System& s, const ForceField& ff, const AdsorptionOptions& o, AdsorptionReport* report) {
  const auto t0 = std::chrono::steady_clock::now();
  AdsorptionReport rep;
  const size_t n = s.atoms.size();
  if (ff.type_index.size() != n || ff.charge.size() != n) throw std::invalid_argument("adsorption: the force field was assigned to another structure");
  if (!ff.pair_func.empty() || ff.hbond.on() || ff.sw.on)
    throw std::invalid_argument("adsorption: the energies need Lennard-Jones pairs (no other pair forms or hydrogen bonds)");
  // massless sites move with their rigid adsorbate (TraPPE N2's centre); on the fixed substrate they have no place
  for (const auto& v : ff.vsites)
    if (o.first_mobile_atom >= 0 && v.site < uint32_t(o.first_mobile_atom))
      throw std::invalid_argument("adsorption: virtual sites on the substrate are not handled");
  if (o.steps < 10 || o.cycles < 1 || o.t_high <= 0 || o.t_low <= 0 || o.t_low > o.t_high) throw std::invalid_argument("adsorption: at least 10 steps and one cycle, 0 < t_low ≤ t_high");
  const auto mobile = o.mobile.empty() ? adsorbate_molecules(s, o.first_mobile_atom) : o.mobile;
  if (mobile.empty()) throw std::invalid_argument("adsorption: no adsorbate molecule (the structure is one molecule)");
  if (!(o.cutoff > 0) || !std::isfinite(o.cutoff)) throw std::invalid_argument("adsorption: the cut-off must be a positive length (Å)");
  // the cut-off holds as asked: past half the cell's narrowest width, every periodic image within it counts
  Model M(s, ff, mobile, o.cutoff, o.coulomb, o.dsf_alpha);
  const size_t nm = M.mol.size();
  for (size_t m = 0; m < nm; ++m)
    if (M.mol[m].empty()) throw std::invalid_argument("adsorption: molecule " + std::to_string(mobile[m]) + " does not exist");
  if (M.images) {
    char b[240];
    std::snprintf(b, sizeof b, "the %.1f Å cut-off is wider than half the cell (%.2f Å across): every periodic image within it counts (%zu substrate images)",
                  o.cutoff, M.wmin, M.gpos.size());
    rep.notes.push_back(b);
    // a molecule within reach of its own images: the cell holds a dense layer, one molecule per cell
    double lmin = 1e30, span = 0;
    for (int i = -2; i <= 2; ++i)
      for (int j = -2; j <= 2; ++j)
        for (int k = -2; k <= 2; ++k) {
          if ((i == 0 && j == 0 && k == 0) || (i && !M.pdim[0]) || (j && !M.pdim[1]) || (k && !M.pdim[2])) continue;
          lmin = std::min(lmin, norm(s.cell.a * double(i) + s.cell.b * double(j) + s.cell.c * double(k)));
        }
    for (const auto& ms : M.mol)
      for (uint32_t x : ms)
        for (uint32_t y : ms) span = std::max(span, norm(s.atoms[y].pos - s.atoms[x].pos));
    if (lmin < o.cutoff + span) {
      int sc[3];
      for (int d = 0; d < 3; ++d) sc[d] = M.pdim[d] ? std::max(1, int(std::ceil(2 * o.cutoff / M.w[d]))) : 1;
      std::snprintf(b, sizeof b, "each adsorbate meets its own periodic images %.2f Å away (one molecule per cell, a dense layer); for one molecule alone make a %dx%dx%d supercell first",
                    lmin, sc[0], sc[1], sc[2]);
      rep.notes.push_back(b);
    }
  }
  if (M.fixed.empty()) rep.notes.push_back("no substrate: every molecule moves (a cluster search)");

  std::mt19937_64 rng(o.seed);
  caps::UniformReal<double> U(0.0, 1.0);
  std::vector<std::vector<Vec3>> P(nm);
  for (size_t m = 0; m < nm; ++m) for (uint32_t a : M.mol[m]) P[m].push_back(s.atoms[a].pos);
  const bool zband = o.z_hi > o.z_lo;
  auto inside = [&](const std::vector<Vec3>& Q) {
    if (!zband) return true;
    for (const auto& q : Q) if (q[2] < o.z_lo || q[2] > o.z_hi) return false;
    return true;
  };
  auto centre = [](const std::vector<Vec3>& Q) { Vec3 c{0, 0, 0}; for (const auto& q : Q) c = c + q; return c * (1.0 / double(Q.size())); };
  // with a z band only the in-plane directions fold: wrapping z would carry a molecule above the
  // cell top to the slab's other face, outside the band it was placed in
  Cell fold_cell = s.cell;
  if (zband) fold_cell.periodic[2] = false;
  auto fold = [&](std::vector<Vec3>& Q) {   // the centre back into the cell (whole molecule shifted)
    if (!M.per) return;
    const Vec3 c = centre(Q), d = fold_cell.wrap(c) - c;
    for (auto& q : Q) q = q + d;
  };
  auto random_axis = [&] {
    const double z = 2 * U(rng) - 1, t = 2 * kPi * U(rng), r = std::sqrt(std::max(0.0, 1 - z * z));
    return Vec3{r * std::cos(t), r * std::sin(t), z};
  };
  auto random_place = [&](const std::vector<Vec3>& Q) {   // a uniform place and a uniform orientation
    const Vec3 c = centre(Q);
    const Vec3 u = random_axis();
    const double ang = 2 * kPi * U(rng);
    Vec3 target;
    if (M.per) {
      Vec3 f{U(rng), U(rng), U(rng)};
      target = s.cell.to_cartesian(f);
      if (zband) target[2] = o.z_lo + (o.z_hi - o.z_lo) * U(rng);
    } else {
      Vec3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
      for (uint32_t b : M.fixed) for (int k = 0; k < 3; ++k) lo[size_t(k)] = std::min(lo[size_t(k)], s.atoms[b].pos[size_t(k)]), hi[size_t(k)] = std::max(hi[size_t(k)], s.atoms[b].pos[size_t(k)]);
      if (M.fixed.empty()) lo = c - Vec3{10, 10, 10}, hi = c + Vec3{10, 10, 10};
      for (int k = 0; k < 3; ++k) target[size_t(k)] = lo[size_t(k)] - 5 + (hi[size_t(k)] - lo[size_t(k)] + 10) * U(rng);
      if (zband) target[2] = o.z_lo + (o.z_hi - o.z_lo) * U(rng);
    }
    std::vector<Vec3> R(Q.size());
    for (size_t k = 0; k < Q.size(); ++k) R[k] = target + rotate(Q[k] - c, u, ang);
    return R;
  };
  auto mol_energy = [&](size_t m, const std::vector<Vec3>& Q) {
    double e = M.with_fixed(m, Q);
    if (!std::isfinite(e)) return e;
    e += M.self(m, Q);
    for (size_t k = 0; k < nm; ++k) if (k != m) e += M.between(m, Q, k, P[k]);
    return e;
  };
  auto total = [&] {
    double e = 0;
    for (size_t m = 0; m < nm; ++m) {
      e += M.with_fixed(m, P[m]) + M.self(m, P[m]);
      for (size_t k = m + 1; k < nm; ++k) e += M.between(m, P[m], k, P[k]);
    }
    return e;
  };
  // start: each molecule inserted in turn at the best of 50 random places
  if (o.randomise)
    for (size_t m = 0; m < nm; ++m) {
      std::vector<Vec3> best = P[m];
      double be = std::numeric_limits<double>::infinity();
      for (int t = 0; t < 50; ++t) {
        auto Q = random_place(P[m]);
        fold(Q);
        if (!inside(Q)) continue;
        const double e = mol_energy(m, Q);
        if (e < be) be = e, best = Q;
      }
      P[m] = best;
    }

  std::vector<AdsorptionConfig> finals;
  std::vector<std::vector<Vec3>> bestP = P;
  double best_e = total();
  std::vector<double> samples;
  long acc = 0, att = 0;
  const int every = std::max(1, o.steps / 100);
  const double ratio = std::log(o.t_high / o.t_low);
  bool stopped = false;
  for (int cyc = 0; cyc < o.cycles && !stopped; ++cyc) {
    if (cyc > 0) P = bestP;
    double cur = total();
    double cyc_best = cur;
    std::vector<std::vector<Vec3>> cyc_bestP = P;
    auto step = [&](double T, double dmax, double amax, double pjump, bool cold) {
      const size_t m = size_t(U(rng) * double(nm)) % nm;
      const double r = U(rng);
      std::vector<Vec3> Q = P[m];
      if (r < pjump) Q = random_place(P[m]);
      else if (r < pjump + 0.5 * (1 - pjump)) {   // translate
        const Vec3 d{(2 * U(rng) - 1) * dmax, (2 * U(rng) - 1) * dmax, (2 * U(rng) - 1) * dmax};
        for (auto& q : Q) q = q + d;
      } else {   // rotate about the centre
        const Vec3 c = centre(Q), u = random_axis();
        const double a = (2 * U(rng) - 1) * amax;
        for (auto& q : Q) q = c + rotate(q - c, u, a);
      }
      ++att;
      fold(Q);
      if (!inside(Q)) return;
      const double e_old = mol_energy(m, P[m]), e_new = mol_energy(m, Q);
      if (!std::isfinite(e_new)) return;
      const double dE = e_new - (std::isfinite(e_old) ? e_old : 1e300);
      if (dE <= 0 || (!cold && U(rng) < std::exp(-dE / (kB * T)))) {
        P[m] = std::move(Q);
        cur += std::isfinite(e_old) ? dE : 0;
        if (!std::isfinite(e_old)) cur = total();
        ++acc;
        if (cur < best_e) best_e = cur, bestP = P;
      }
    };
    for (int k = 0; k < o.steps; ++k) {
      const double f = 1.0 - double(k) / std::max(1, o.steps - 1);    // 1 hot … 0 cold
      const double T = o.t_low * std::exp(ratio * f);
      step(T, 0.2 + 2.8 * f, (5 + 55 * f) * kPi / 180, 0.3 * f, false);   // jumps fade out as it cools
      if (cur < cyc_best) cyc_best = cur, cyc_bestP = P;
      if (k >= o.steps * 4 / 5 && k % 10 == 0) samples.push_back(cur);
      if ((k + 1) % every == 0 && o.progress && !o.progress(cyc, k + 1, best_e)) { stopped = true; break; }
    }
    P = cyc_bestP, cur = cyc_best;   // the cycle's lowest configuration, polished
    {   // polish downhill with steps shrinking from 0.3 to 0.005 Å (and 10° to 0.2°)
      const int np = std::max(400, o.steps / 5);
      for (int k = 0; k < np; ++k) {
        const double f = std::pow(0.005 / 0.3, double(k) / (np - 1));
        step(0, 0.3 * f, 10 * f * kPi / 180, 0, true);
      }
    }
    cur = total();
    if (cur < best_e) best_e = cur, bestP = P;
    AdsorptionConfig c;
    c.energy = cur, c.cycle = cyc + 1;
    c.positions.resize(n);
    for (size_t i = 0; i < n; ++i) c.positions[i] = s.atoms[i].pos;
    for (size_t m = 0; m < nm; ++m) for (size_t k = 0; k < M.mol[m].size(); ++k) c.positions[M.mol[m][k]] = P[m][k];
    finals.push_back(std::move(c));
  }
  if (stopped) rep.notes.push_back("stopped before the last cycle ended");
  // the lowest configuration into the structure
  P = bestP;
  for (size_t m = 0; m < nm; ++m) for (size_t k = 0; k < M.mol[m].size(); ++k) s.atoms[M.mol[m][k]].pos = P[m][k];
  s.velocities.clear();
  {
    AdsorptionConfig b;
    b.energy = best_e, b.cycle = 0;
    b.positions.resize(n);
    for (size_t i = 0; i < n; ++i) b.positions[i] = s.atoms[i].pos;
    finals.push_back(std::move(b));
  }
  std::sort(finals.begin(), finals.end(), [](const auto& a, const auto& b) { return a.energy < b.energy; });
  for (const auto& c : finals) {
    if (!rep.configs.empty() && std::fabs(c.energy - rep.configs.back().energy) < 1e-6) continue;   // the same minimum twice
    rep.configs.push_back(c);
    if (int(rep.configs.size()) >= std::max(1, o.keep)) break;
  }
  // energies of the lowest configuration, by kind and per component
  double as = 0, aa = 0;
  std::vector<double> per(nm, 0.0);
  for (size_t m = 0; m < nm; ++m) {
    const double f = M.with_fixed(m, P[m]), own = M.self(m, P[m]);
    as += f, per[m] += f;
    aa += own, per[m] += own;   // with its own images, as the other adsorbates
    for (size_t k = m + 1; k < nm; ++k) {
      const double e = M.between(m, P[m], k, P[k]);
      aa += e, per[m] += e, per[k] += e;
    }
  }
  rep.adsorption_energy = as + aa, rep.adsorbate_substrate = as, rep.adsorbate_adsorbate = aa;
  std::map<std::string, std::pair<int, double>> comp;
  std::vector<std::string> order;
  for (size_t m = 0; m < nm; ++m) {
    const std::string f = formula(s, M.mol[m]);
    if (!comp.count(f)) order.push_back(f);
    auto& c = comp[f];
    ++c.first, c.second += per[m];
  }
  for (const auto& f : order) rep.components.push_back({f, comp[f].first, comp[f].second / comp[f].first});
  if (!samples.empty()) {
    const double lo = *std::min_element(samples.begin(), samples.end()), hi = *std::max_element(samples.begin(), samples.end());
    const int nb = 30;
    const double w = hi > lo ? (hi - lo) / nb : 1.0;
    for (int b = 0; b <= nb; ++b) rep.hist_edges.push_back(lo + b * w);
    rep.hist_counts.assign(size_t(nb), 0.0);
    for (double x : samples) rep.hist_counts[size_t(std::clamp(int((x - lo) / w), 0, nb - 1))] += 1;
  }
  rep.steps = att;
  rep.acceptance = att ? double(acc) / double(att) : 0;
  {
    char b[240];
    std::snprintf(b, sizeof b, "%zu adsorbate molecule%s on %zu fixed atoms · %d cycle%s of %d steps, %.0f → %.0f K · lowest %.3f kcal/mol",
                  nm, nm == 1 ? "" : "s", M.fixed.size(), o.cycles, o.cycles == 1 ? "" : "s", o.steps, o.t_high, o.t_low, best_e);
    rep.notes.insert(rep.notes.begin(), b);
    if (best_e > 0) {
      std::snprintf(b, sizeof b, "the lowest energy found is positive (%.1f kcal/mol): the adsorbates do not fit (a bulk cell has no room; a surface needs vacuum above it)", best_e);
      rep.notes.push_back(b);
    }
    rep.notes.push_back(std::string("energies: the force field's van der Waals") + (M.pm.coulomb() ? " and damped shifted force electrostatics" : "") +
                        " between molecules; rigid adsorbates (no deformation energy)");
  }
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (report) *report = std::move(rep);
}

}  // namespace caps
