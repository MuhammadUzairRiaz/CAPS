// CAPS grafting chains onto an existing structure (see graft_chains.hpp).
#include "caps/graft_chains.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

#include "caps/rng.hpp"

namespace caps {

namespace {
struct Env {
  int element = 0, hydrogens = 0;
  std::vector<int> neighbours;   // heavy neighbours' elements, sorted
  std::string resname;
  bool operator==(const Env& o) const { return element == o.element && hydrogens == o.hydrogens && neighbours == o.neighbours && resname == o.resname; }
};
Env env_of(const System& s, const std::vector<std::vector<uint32_t>>& nb, uint32_t i) {
  Env e;
  e.element = s.atoms[i].element;
  e.resname = s.atoms[i].resname;
  for (uint32_t j : nb[i]) {
    if (s.atoms[j].element == 1) ++e.hydrogens;
    else e.neighbours.push_back(s.atoms[j].element);
  }
  std::sort(e.neighbours.begin(), e.neighbours.end());
  return e;
}
}  // namespace

std::vector<uint32_t> equivalent_atoms(const System& s, uint32_t atom) {
  if (atom >= s.atoms.size()) throw std::invalid_argument("graft: atom " + std::to_string(atom + 1) + " is not in the structure");
  const auto nb = s.neighbours();
  const Env ref = env_of(s, nb, atom);
  std::vector<uint32_t> out;
  for (uint32_t i = 0; i < s.atoms.size(); ++i)
    if (env_of(s, nb, i) == ref) out.push_back(i);
  return out;
}

System graft_chains(const System& host, const ChainSpec& graft, const ChainGraftOptions& o, ChainGraftReport* report) {
  ChainGraftReport R;
  if (host.atoms.empty()) throw GrowError("graft: the host has no atoms");
  if (graft.units.empty()) throw GrowError("graft: the grafts' repeat units");
  if (o.sites.empty()) throw GrowError("graft: choose the site atoms");
  if (o.per_site < 1) throw GrowError("graft: at least one graft per site");
  const auto nb = host.neighbours();
  const size_t n = host.atoms.size();
  // candidates
  std::vector<uint32_t> cand;
  {
    std::set<uint32_t> seen;
    for (uint32_t a0 : o.sites) {
      if (a0 >= n) throw GrowError("graft: atom " + std::to_string(a0 + 1) + " is not in the structure");
      uint32_t a = a0;   // a hydrogen picked stands for the atom it sits on
      if (host.atoms[a].element == 1 && nb[a].size() == 1) a = nb[a][0];
      for (uint32_t e : o.equivalent ? equivalent_atoms(host, a) : std::vector<uint32_t>{a})
        if (seen.insert(e).second) cand.push_back(e);
    }
  }
  R.candidates = int(cand.size());
  // the cell the grafts grow in: the host's (orthorhombic), else a box with room around the host
  const bool periodic = host.cell.valid();
  std::array<double, 3> L{0, 0, 0};
  Vec3 shift{0, 0, 0};
  if (periodic) {
    const Cell& c = host.cell;
    if (std::fabs(c.a[1]) > 1e-6 || std::fabs(c.a[2]) > 1e-6 || std::fabs(c.b[0]) > 1e-6 || std::fabs(c.b[2]) > 1e-6 || std::fabs(c.c[0]) > 1e-6 || std::fabs(c.c[1]) > 1e-6)
      throw GrowError("graft: the host's cell is not orthorhombic");
    L = {c.a[0], c.b[1], c.c[2]};
    shift = Vec3{0, 0, 0} - c.origin;
  } else {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& a : host.atoms)
      for (int k = 0; k < 3; ++k) lo[size_t(k)] = std::min(lo[size_t(k)], a.pos[size_t(k)]), hi[size_t(k)] = std::max(hi[size_t(k)], a.pos[size_t(k)]);
    const double pad = o.padding > 0 ? o.padding : std::max(10.0, 1.3 * graft.dp + 6.0);
    for (int k = 0; k < 3; ++k) L[size_t(k)] = hi[size_t(k)] - lo[size_t(k)] + 2 * pad, shift[size_t(k)] = pad - lo[size_t(k)];
    R.notes.push_back("the host has no cell: grafted in a box with " + std::to_string(int(std::lround(pad))) + " Å of room around it");
  }
  Cell box;
  box.a = {L[0], 0, 0}, box.b = {0, L[1], 0}, box.c = {0, 0, L[2]};
  auto mi = [&](Vec3 d) { return box.minimum_image(d); };
  // the draw: kept with the probability, apart by min_spacing, up to max_grafts
  std::mt19937_64 rng(o.grow.seed * 0x9E3779B97F4A7C15ull + 31);
  caps::UniformReal<double> U(0, 1);
  std::vector<uint32_t> order = cand;
  caps::shuffle(order.begin(), order.end(), rng);
  std::vector<uint32_t> pick;
  for (uint32_t a : order) {
    if (o.max_grafts > 0 && int(pick.size()) * o.per_site >= o.max_grafts) break;
    if (U(rng) >= o.probability) continue;
    int hs = 0;
    for (uint32_t j : nb[a]) hs += host.atoms[j].element == 1;
    if (hs < o.per_site) {
      R.notes.push_back("atom " + std::to_string(a + 1) + " has " + std::to_string(hs) + " hydrogen(s), fewer than the grafts asked for: left out");
      continue;
    }
    bool ok = true;
    for (uint32_t p : pick)
      if (norm(mi(host.atoms[a].pos - host.atoms[p].pos)) < o.min_spacing) { ok = false; break; }
    if (ok) pick.push_back(a);
  }
  if (pick.empty()) throw GrowError("graft: no site was drawn (" + std::to_string(cand.size()) + " candidates; raise the probability or check the sites have hydrogens)");
  // each site's hydrogens with the most room, given up for the grafts
  std::vector<char> drop(n, 0);
  struct Site { uint32_t atom; Vec3 dir; };
  std::vector<Site> sites;
  for (uint32_t a : pick) {
    std::vector<std::pair<double, uint32_t>> hs;
    for (uint32_t h : nb[a]) {
      if (host.atoms[h].element != 1) continue;
      const Vec3 d = mi(host.atoms[h].pos - host.atoms[a].pos);
      const Vec3 probe = host.atoms[a].pos + d * (2.5 / std::max(1e-9, norm(d)));
      double room = 1e300;
      for (size_t j = 0; j < n; ++j)
        if (j != a && j != h) room = std::min(room, norm(mi(host.atoms[j].pos - probe)));
      hs.push_back({-room, h});
    }
    std::sort(hs.begin(), hs.end());
    for (int k = 0; k < o.per_site && k < int(hs.size()); ++k) {
      const uint32_t h = hs[size_t(k)].second;
      drop[h] = 1;
      const Vec3 d = mi(host.atoms[h].pos - host.atoms[a].pos);
      sites.push_back({a, d * (1.0 / norm(d))});
      R.sites.push_back(a);
    }
  }
  // the host without those hydrogens, in the box
  std::vector<int64_t> newi(n, -1);
  System sub = host;
  sub.atoms.clear(), sub.bonds.clear(), sub.velocities.clear();
  for (size_t i = 0; i < n; ++i)
    if (!drop[i]) {
      newi[i] = int64_t(sub.atoms.size());
      Atom at = host.atoms[i];
      at.pos = at.pos + shift;
      if (periodic)
        for (int k = 0; k < 3; ++k) at.pos[size_t(k)] -= L[size_t(k)] * std::floor(at.pos[size_t(k)] / L[size_t(k)]);
      sub.atoms.push_back(at);
    }
  for (const auto& b : host.bonds)
    if (newi[b.i] >= 0 && newi[b.j] >= 0) sub.bonds.push_back({uint32_t(newi[b.i]), uint32_t(newi[b.j]), b.order});
  sub.cell = box;
  GrowOptions g = o.grow;
  g.cell = L;
  g.box = 0, g.density = 0;
  g.substrate = &sub;
  g.chains = int(sites.size());
  g.anchors.clear();
  for (const auto& st : sites) g.anchors.push_back({uint32_t(newi[st.atom]), {st.dir[0], st.dir[1], st.dir[2]}});
  g.auto_scale = true;
  ChainSpec spec = graft;
  spec.chain_dp.clear();
  System s = grow_chains(spec, g, &R.grow);
  // positions back where the host was, molecules by connectivity
  for (auto& a : s.atoms) a.pos = a.pos - shift;
  if (periodic) s.cell = host.cell;
  else s.cell = Cell{};
  int nmol = 0;
  s.has_mol = false;   // by bonds, not the ids grow gave each chain
  const auto mol = s.molecules(&nmol);
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].mol = int64_t(mol[i]) + 1;
  s.has_mol = true;
  R.grafted = int(sites.size());
  s.title = "CAPS graft: " + std::to_string(R.grafted) + " grafts of " + std::to_string(graft.dp) + " units on " + (host.title.empty() ? std::string("the structure") : host.title);
  for (const auto& nt : R.grow.notes) R.notes.push_back(nt);
  if (report) *report = R;
  return s;
}

}  // namespace caps
