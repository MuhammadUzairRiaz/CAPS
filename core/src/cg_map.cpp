// Structure-based coarse-graining (see cg_map.hpp).
#include "caps/cg_map.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/resolution.hpp"

namespace caps {

namespace {
constexpr double kB = 0.0019872067;   // kcal/(mol K)
constexpr double kPi = 3.14159265358979323846;
constexpr double kContact = 0.36787944117144233;   // 1/e: where −k_B T ln g(r) = k_B T

std::pair<double, double> mean_sd(const std::vector<double>& v) {
  if (v.empty()) return {0, 0};
  double m = 0;
  for (double x : v) m += x;
  m /= double(v.size());
  double ss = 0;
  for (double x : v) ss += (x - m) * (x - m);
  return {m, v.size() > 1 ? std::sqrt(ss / double(v.size() - 1)) : 0.0};
}
}  // namespace

CgMapResult cg_map(const System& aa, const CgMapOptions& o, const std::vector<std::vector<Vec3>>& frames_in, const std::vector<Cell>& cells_in) {
  const size_t n = aa.atoms.size();
  if (n == 0) throw std::invalid_argument("no atoms to coarse-grain");
  if (o.temperature <= 0) throw std::invalid_argument("give the temperature of the inversion");
  // ---- which bead each atom joins, and each bead's type
  std::vector<int> bead(n, -1);
  std::vector<std::string> btype;
  std::vector<int64_t> bmol;
  const auto nb = aa.neighbours();
  int nm = 0;
  const auto mol = aa.molecules(&nm);
  if (o.scheme == "unit" || o.scheme == "backbone_side") {
    bool numbered = false;
    for (const auto& a : aa.atoms) numbered = numbered || a.resid > 0;
    if (!numbered)
      throw std::invalid_argument("the atoms carry no repeat-unit numbers (a cell grown in CAPS numbers them): use n backbone atoms per bead instead");
    std::vector<char> backbone(n, 0);
    if (o.scheme == "backbone_side") {
      for (const auto& path : backbones(aa, 2))
        for (uint32_t i : path) backbone[i] = 1;
      for (uint32_t i = 0; i < n; ++i)   // hydrogens go with their heavy atom
        if (aa.atoms[i].element == 1)
          for (uint32_t w : nb[i]) if (backbone[w]) backbone[i] = 1;
    }
    std::map<std::tuple<int, int64_t, int>, int> key;
    for (uint32_t i = 0; i < n; ++i) {
      const auto k = std::make_tuple(mol[i], aa.atoms[i].resid, o.scheme == "backbone_side" ? (backbone[i] ? 0 : 1) : 0);
      auto [it, fresh] = key.emplace(k, int(btype.size()));
      if (fresh) {
        const std::string res = aa.atoms[i].resname.empty() ? std::string("U") : aa.atoms[i].resname;
        btype.push_back(o.scheme == "unit" ? res : res + (std::get<2>(k) ? "_S" : "_B"));
        bmol.push_back(mol[i]);
      }
      bead[i] = it->second;
    }
  } else if (o.scheme == "backbone_n") {
    if (o.per_bead < 1) throw std::invalid_argument("at least one backbone atom per bead");
    ResolutionReport rr;
    const System cg = coarse_grain(aa, o.per_bead, &rr);
    bead = rr.site_of;
    for (const auto& a : cg.atoms) btype.push_back("B" + std::to_string(o.per_bead)), bmol.push_back(a.mol - 1);
  } else if (o.scheme == "rules") {
    const CgMapping m = cg_mapping(aa, o.rules);
    bead = m.bead_of();
    btype = m.bead_kind;
    for (int mm : m.bead_mol) bmol.push_back(mm);
  } else {
    throw std::invalid_argument("mapping scheme: unit, backbone_side, backbone_n or rules");
  }
  const size_t nbead = btype.size();
  // masses: each bead's atoms'
  std::vector<double> bmass(nbead, 0.0);
  for (uint32_t i = 0; i < n; ++i) bmass[size_t(bead[i])] += element(aa.atoms[i].element).mass;
  // ---- bead topology from the atom bonds
  std::set<std::pair<uint32_t, uint32_t>> bbond;
  for (const auto& bd : aa.bonds) {
    uint32_t x = uint32_t(bead[bd.i]), y = uint32_t(bead[bd.j]);
    if (x == y) continue;
    if (x > y) std::swap(x, y);
    bbond.insert({x, y});
  }
  std::vector<std::vector<uint32_t>> bnb(nbead);
  for (const auto& [x, y] : bbond) bnb[x].push_back(y), bnb[y].push_back(x);
  std::vector<std::array<uint32_t, 3>> bangle;
  for (uint32_t j = 0; j < nbead; ++j)
    for (size_t a = 0; a < bnb[j].size(); ++a)
      for (size_t c = a + 1; c < bnb[j].size(); ++c) bangle.push_back({bnb[j][a], j, bnb[j][c]});
  // ---- every frame mapped (centres of mass, each bead unwrapped about its first atom)
  std::vector<std::vector<Vec3>> frames = frames_in;
  std::vector<Cell> cells = cells_in;
  if (frames.empty()) {
    std::vector<Vec3> p;
    for (const auto& a : aa.atoms) p.push_back(a.pos);
    frames.push_back(p);
    cells = {aa.cell};
  }
  if (cells.size() != frames.size()) cells.assign(frames.size(), aa.cell);
  std::vector<int> first(nbead, -1);
  for (uint32_t i = 0; i < n; ++i) if (first[size_t(bead[i])] < 0) first[size_t(bead[i])] = int(i);
  CgMapResult r;
  for (size_t f = 0; f < frames.size(); ++f) {
    if (frames[f].size() != n) throw std::invalid_argument("frame " + std::to_string(f + 1) + " has another number of atoms");
    const Cell& c = cells[f];
    std::vector<Vec3> com(nbead, Vec3{0, 0, 0});
    for (uint32_t i = 0; i < n; ++i) {
      const size_t b = size_t(bead[i]);
      const Vec3 ref = frames[f][size_t(first[b])];
      const Vec3 p = c.valid() ? ref + c.minimum_image(frames[f][i] - ref) : frames[f][i];
      com[b] = com[b] + p * element(aa.atoms[i].element).mass;
    }
    for (size_t b = 0; b < nbead; ++b) com[b] = com[b] * (1.0 / bmass[b]);
    r.frames.push_back(std::move(com));
    r.cells.push_back(c);
  }
  // ---- bonded distributions by type, inverted (harmonic)
  auto dist = [&](size_t f, uint32_t x, uint32_t y) {
    Vec3 d = r.frames[f][y] - r.frames[f][x];
    if (r.cells[f].valid()) d = r.cells[f].minimum_image(d);
    return d;
  };
  const double kT = kB * o.temperature;
  std::map<std::pair<std::string, std::string>, std::vector<double>> blen;
  for (const auto& [x, y] : bbond) {
    auto key = std::minmax(btype[x], btype[y]);
    for (size_t f = 0; f < r.frames.size(); ++f) blen[{key.first, key.second}].push_back(norm(dist(f, x, y)));
  }
  std::map<std::tuple<std::string, std::string, std::string>, std::vector<double>> bang;
  for (const auto& t : bangle) {
    std::string a = btype[t[0]], c = btype[t[2]];
    if (c < a) std::swap(a, c);
    for (size_t f = 0; f < r.frames.size(); ++f) {
      const Vec3 u = dist(f, t[1], t[0]), v = dist(f, t[1], t[2]);
      const double cs = std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0);
      bang[{a, btype[t[1]], c}].push_back(std::acos(cs));
    }
  }
  std::map<std::pair<std::string, std::string>, size_t> bond_ix;
  for (const auto& [k, v] : blen) {
    const auto [m, sd] = mean_sd(v);
    CgBondType bt{k.first, k.second, m, sd > 1e-6 ? kT / (2 * sd * sd) : 0.0, sd, int(v.size())};
    if (bt.k <= 0) throw std::invalid_argument("bond " + k.first + "–" + k.second + " has no spread to invert (a single rigid sample): map more chains or frames");
    bond_ix[k] = r.bonds.size();
    r.bonds.push_back(bt);
  }
  std::map<std::tuple<std::string, std::string, std::string>, size_t> angle_ix;
  for (const auto& [k, v] : bang) {
    const auto [m, sd] = mean_sd(v);
    if (sd <= 1e-6) continue;
    angle_ix[k] = r.angles.size();
    r.angles.push_back({std::get<0>(k), std::get<1>(k), std::get<2>(k), m, kT / (2 * sd * sd), sd, int(v.size())});
  }
  // ---- non-bonded bead–bead g(r), pooled over types and frames: the repulsive size
  {
    const double dr = 0.1;
    double rmax = 15.0;
    const Cell& c0 = r.cells.back();
    const bool periodic = c0.valid();
    if (periodic) rmax = std::min(rmax, 0.5 * std::min({norm(c0.a), norm(c0.b), norm(c0.c)}));
    const size_t nbin = size_t(rmax / dr);
    std::vector<double> h(nbin, 0.0);
    // the non-bonded pairs: other chains, or more than three bonds apart along one (intermolecular pairs alone carry
    // the chains' correlation hole, which keeps g(r) under 1 out to the chain size)
    std::vector<std::set<uint32_t>> near(nbead);
    for (uint32_t x = 0; x < nbead; ++x) {
      std::vector<uint32_t> front{x};
      for (int step = 0; step < 3; ++step) {
        std::vector<uint32_t> next;
        for (uint32_t u : front)
          for (uint32_t w : bnb[u])
            if (w != x && near[x].insert(w).second) next.push_back(w);
        front = next;
      }
    }
    double pairs = 0, vol = 0;
    std::map<std::pair<std::string, std::string>, std::pair<std::vector<double>, double>> hp;   // per type pair: histogram, pairs
    for (size_t f = 0; f < r.frames.size(); ++f) {
      const double V = r.cells[f].valid() ? r.cells[f].volume() : 0;
      for (size_t x = 0; x < nbead; ++x)
        for (size_t y = x + 1; y < nbead; ++y) {
          if (bmol[x] == bmol[y] && near[x].count(uint32_t(y))) continue;
          pairs += 1;
          auto key = std::minmax(btype[x], btype[y]);
          auto& e = hp[{key.first, key.second}];
          if (e.first.empty()) e.first.assign(nbin, 0.0);
          e.second += 1;
          const double d = norm(dist(f, uint32_t(x), uint32_t(y)));
          const size_t k = size_t(d / dr);   // d just under rmax can round to nbin
          if (d < rmax && k < nbin) h[k] += 1, e.first[k] += 1;
        }
      vol += V;
    }
    // each type pair's own cut where its g(r) first reaches 1, when it has the samples (several bead types)
    if (periodic && vol > 0 && hp.size() > 1) {
      const double V = vol / double(r.frames.size());
      for (const auto& [k, e] : hp) {
        const double per_frame = e.second / double(r.frames.size());
        double below = 0, cut = 0;
        std::vector<double> g(nbin);
        for (size_t b = 0; b < nbin; ++b) {
          const double r0 = b * dr, r1 = r0 + dr, shell = 4.0 / 3.0 * kPi * (r1 * r1 * r1 - r0 * r0 * r0);
          g[b] = e.first[b] / double(r.frames.size()) / (per_frame * shell / V);
        }
        for (size_t b = 1; b + 1 < nbin; ++b) {
          below += e.first[b - 1];
          if ((g[b - 1] + g[b] + g[b + 1]) / 3 >= kContact) { cut = (b + 0.5) * dr * std::pow(2.0, 1.0 / 6); break; }
        }
        if (cut > 0 && below >= 30) r.pair_cut[k] = cut;   // enough pairs inside it to trust the crossing
      }
    }
    if (periodic && pairs > 0 && vol > 0) {
      const double V = vol / double(r.frames.size()), per_frame = pairs / double(r.frames.size());
      for (size_t k = 0; k < nbin; ++k) {
        const double r0 = k * dr, r1 = r0 + dr, shell = 4.0 / 3.0 * kPi * (r1 * r1 * r1 - r0 * r0 * r0);
        const double g = h[k] / double(r.frames.size()) / (per_frame * shell / V);
        r.gr.push_back({r0 + dr / 2, g});
      }
      // σ where g, smoothed over three bins, first reaches 1/e: there −k_B T ln g(r) (the IBI starting potential) is
      // k_B T, which a WCA with ε = k_B T is at r = σ
      for (size_t k = 1; k + 1 < r.gr.size(); ++k) {
        const double g = (r.gr[k - 1].second + r.gr[k].second + r.gr[k + 1].second) / 3;
        if (g >= kContact) { r.cut = r.gr[k].first * std::pow(2.0, 1.0 / 6); break; }
      }
    }
    if (r.cut <= 0) throw std::invalid_argument(periodic ? "the non-bonded bead g(r) never reaches 1/e within half the box: map a denser or larger cell"
                                                         : "a periodic cell (a melt) is needed for the beads' repulsive size");
    r.sigma = r.cut / std::pow(2.0, 1.0 / 6);
    r.epsilon = kT;
  }
  // ---- the bead structure (the last frame) and its force field
  std::map<std::string, int> tix;
  std::vector<std::string> tnames;
  for (const auto& t : btype)
    if (tix.emplace(t, int(tnames.size())).second) tnames.push_back(t);
  std::vector<double> tmass(tnames.size(), 0.0), tcount(tnames.size(), 0.0);
  for (size_t b = 0; b < nbead; ++b) tmass[size_t(tix[btype[b]])] += bmass[b], tcount[size_t(tix[btype[b]])] += 1;
  for (size_t t = 0; t < tnames.size(); ++t) tmass[t] /= tcount[t];
  System& s = r.beads;
  s.title = (aa.title.empty() ? std::string("polymer") : aa.title) + " (coarse-grained, " + o.scheme + ")";
  s.cell = r.cells.back();
  s.has_mol = true;
  s.bonds_from_file = true;
  s.has_charges = true;
  s.unwrapped = true;
  for (size_t t = 0; t < tnames.size(); ++t) s.types.push_back({int(t + 1), tmass[t], tnames[t]});
  for (size_t b = 0; b < nbead; ++b) {
    Atom a;
    a.id = int64_t(b + 1);
    a.mol = bmol[b] + 1;
    a.type = tix[btype[b]] + 1;
    a.name = btype[b];
    a.element = 0;
    a.pos = r.frames.back()[b];
    s.atoms.push_back(a);
  }
  for (const auto& [x, y] : bbond) s.bonds.push_back({x, y, 1});
  auto F = std::make_shared<ForceField>();
  char b[240];
  std::snprintf(b, sizeof b, "Structure-based CG (%s) · Boltzmann-inverted bonds and angles, WCA σ = %.3g Å at %.4g K", o.scheme.c_str(), r.sigma, o.temperature);
  F->name = b;
  F->type_names = tnames;
  F->lj.assign(tnames.size(), PairType{r.epsilon, r.sigma});
  if (tnames.size() > 1) {   // a WCA per type pair, cut where that pair's g(r) reaches 1 (LAMMPS cosine/squared … wca)
    double maxcut = r.cut;
    for (size_t a = 0; a < tnames.size(); ++a)
      for (size_t c = a; c < tnames.size(); ++c) {
        auto key = std::minmax(tnames[a], tnames[c]);
        auto it = r.pair_cut.find({key.first, key.second});
        const double cut = it != r.pair_cut.end() ? it->second : r.cut;
        maxcut = std::max(maxcut, cut);
        F->pair_func[{int(a), int(c)}] = PairFunc{kPairCos2Wca, r.epsilon, cut, cut};
      }
    r.cut_max = maxcut;
  }
  F->pair_form = "lj12-6";
  F->mixing = "arithmetic";
  F->cutoff = std::max(r.cut, r.cut_max);
  F->lj_shift = true;
  F->lj14 = 1.0, F->coul14 = 1.0;   // special_bonds lj 0 0 1: 1-2 and 1-3 excluded below, 1-4 in full
  F->native_pair = "lj/cut";
  F->type_index.resize(nbead);
  F->atom_type.resize(nbead);
  F->why.assign(nbead, "mapped bead (" + o.scheme + ")");
  F->charge.assign(nbead, 0.0);
  F->mass.resize(nbead);
  for (size_t k = 0; k < nbead; ++k) {
    F->type_index[k] = tix[btype[k]];
    F->atom_type[k] = btype[k];
    F->mass[k] = tmass[size_t(tix[btype[k]])];
  }
  F->excluded.assign(nbead, {});
  double shortest_period = 1e30;
  for (const auto& [x, y] : bbond) {
    auto key = std::minmax(btype[x], btype[y]);
    const CgBondType& bt = r.bonds[bond_ix[{key.first, key.second}]];
    F->bonds.push_back({x, y, bt.k, bt.r0});
    F->excluded[x].push_back(y), F->excluded[y].push_back(x);
    const double mu = F->mass[x] * F->mass[y] / (F->mass[x] + F->mass[y]);
    shortest_period = std::min(shortest_period, 2 * kPi * std::sqrt(mu / (2 * bt.k)) * 48.88821291);   // fs
  }
  for (const auto& t : bangle) {
    std::string a = btype[t[0]], c = btype[t[2]];
    if (c < a) std::swap(a, c);
    auto it = angle_ix.find({a, btype[t[1]], c});
    if (it != angle_ix.end()) F->angles.push_back({t[0], t[1], t[2], r.angles[it->second].k, r.angles[it->second].theta0});
    F->excluded[t[0]].push_back(t[2]), F->excluded[t[2]].push_back(t[0]);
  }
  for (auto& e : F->excluded) {
    std::sort(e.begin(), e.end());
    e.erase(std::unique(e.begin(), e.end()), e.end());
  }
  // a step a twentieth of the fastest bond's period, in whole half femtoseconds
  F->native_timestep = std::max(0.5, std::floor(shortest_period / 20 / 0.5) * 0.5);
  std::snprintf(b, sizeof b, "%zu beads (%zu types) from %zu atoms · %zu bond and %zu angle types inverted from %zu frame(s) · σ = %.3f Å (g(r) of the non-bonded "
                "beads reaches 1/e at σ; cut %.3f Å = 2^(1/6) σ), ε = k_B T = %.4f kcal/mol · Δt ≤ %.1f fs",
                nbead, tnames.size(), n, r.bonds.size(), r.angles.size(), r.frames.size(), r.sigma, r.cut, r.epsilon, F->native_timestep);
  r.notes.push_back(b);
  if (!r.pair_cut.empty()) {
    std::string pc;
    for (const auto& [k, c] : r.pair_cut) { char q[80]; std::snprintf(q, sizeof q, "%s%s–%s %.2f Å", pc.empty() ? "" : ", ", k.first.c_str(), k.second.c_str(), c); pc += q; }
    r.notes.push_back("each type pair's WCA from its own g(r) (σ where it reaches 1/e), cut: " + pc + " (other pairs: the pooled " + std::to_string(r.cut).substr(0, 5) + " Å)");
  }
  r.notes.push_back("a structure-based starting model: bonded terms reproduce the mapped distributions; the one repulsive size does not give the melt's "
                    "density or g(r) — refine the non-bonded potential by iterative Boltzmann inversion (Reith, Pütz & Müller-Plathe 2003) for that");
  if (r.frames.size() == 1) r.notes.push_back("one frame: the distributions are over the chains of a single configuration — map an equilibrated trajectory for better statistics");
  F->notes = r.notes;
  r.ff = F;
  s.forcefield = F;
  return r;
}

}  // namespace caps
