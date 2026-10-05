// CAPS functionalisation of fillers (see functionalize.hpp).
#include "caps/rng.hpp"
#include "caps/functionalize.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/edit.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;

const std::vector<std::pair<std::string, std::string>>& groups() {
  static const std::vector<std::pair<std::string, std::string>> g = {
      {"hydroxyl", "*O"},        {"carboxyl", "*C(=O)O"},  {"amine", "*N"},         {"methyl", "*C"},
      {"fluoro", "*F"},          {"phenyl", "*c1ccccc1"},  {"nitrophenyl", "*c1ccc(cc1)[N+](=O)[O-]"},
      {"amide", "*C(=O)N"},      {"ester", "*C(=O)OC"},    {"hydroxymethyl", "*CO"}, {"vinyl", "*C=C"},
      {"thiol", "*S"},           {"aminopropyl", "*CCCN"}, {"octadecylamide", "*C(=O)NCCCCCCCCCCCCCCCCCC"},
      {"peg3", "*OCCOCCOCCO"}};
  return g;
}

Vec3 unit(const Vec3& v) { const double n = norm(v); return n > 1e-12 ? v * (1 / n) : Vec3{0, 0, 0}; }

}  // namespace

const std::vector<std::string>& functional_group_names() {
  static const std::vector<std::string> n = [] { std::vector<std::string> v; for (const auto& [k, s] : groups()) v.push_back(k); return v; }();
  return n;
}

std::string functional_group_smiles(const std::string& name) {
  for (const auto& [k, s] : groups()) if (k == name) return s;
  if (name.find('*') == std::string::npos) throw std::invalid_argument("'" + name + "' is neither a preset group nor a SMILES with an attachment point *");
  return name;
}

FunctionalizeReport functionalize(System& s, const FunctionalizeOptions& o) {
  FunctionalizeReport rep;
  const std::string smi = functional_group_smiles(o.group);
  if (std::count(smi.begin(), smi.end(), '*') != 1) throw std::invalid_argument("the group needs exactly one attachment point *");
  const size_t n = s.atoms.size();
  const auto nb = s.neighbours();
  const bool per = s.cell.valid();
  auto mi = [&](Vec3 d) { return per ? s.cell.minimum_image(d) : d; };
  // which atoms can carry a group
  std::set<int> zs;
  {
    std::string t = o.elements;
    for (auto& c : t) if (c == ',' || c == ';') c = ' ';
    std::istringstream in(t);
    for (std::string e; in >> e;) {
      const int z = element_from_symbol(e);
      if (z <= 0) throw std::invalid_argument("unknown element '" + e + "'");
      zs.insert(z);
    }
  }
  std::vector<uint32_t> frame;
  for (uint32_t i = 0; i < n; ++i)
    if (s.atoms[i].element != 1 && (zs.empty() || zs.count(s.atoms[i].element))) frame.push_back(i);
  if (frame.empty()) {   // say which elements the structure has, so the filter can be fixed
    std::set<int> have;
    for (const auto& a : s.atoms) if (a.element > 1) have.insert(a.element);
    std::string list;
    for (int z : have) list += (list.empty() ? "" : ", ") + std::string(element(z).symbol);
    throw std::invalid_argument("no atom of the chosen elements: this structure has " + (list.empty() ? std::string("none but hydrogen") : list) +
                                " — clear 'On elements' (any) or name one of these");
  }
  auto hydrogen_of = [&](uint32_t i) { for (uint32_t q : nb[i]) if (s.atoms[q].element == 1) return int(q); return -1; };
  // the filler's centre, extents and the directions it spans across the cell (periodic: no "outward" there)
  Vec3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30}, c{0, 0, 0};
  for (uint32_t i : frame) {
    c = c + s.atoms[i].pos;
    for (int k = 0; k < 3; ++k) lo[size_t(k)] = std::min(lo[size_t(k)], s.atoms[i].pos[size_t(k)]), hi[size_t(k)] = std::max(hi[size_t(k)], s.atoms[i].pos[size_t(k)]);
  }
  c = c * (1.0 / double(frame.size()));
  bool spans[3] = {false, false, false};
  if (per) {
    const double L[3] = {norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)};
    for (int k = 0; k < 3; ++k) spans[k] = s.cell.periodic[size_t(k)] && hi[size_t(k)] - lo[size_t(k)] > 0.8 * L[k];
  }
  int axis = 2;   // band / helix: the longest direction (z on ties)
  for (int k = 0; k < 3; ++k) if (hi[size_t(k)] - lo[size_t(k)] > hi[size_t(axis)] - lo[size_t(axis)] + 0.5) axis = k;
  std::mt19937_64 rng(o.seed);
  caps::UniformReal<double> U(0.0, 1.0);
  // the outward normal of a three-connected site
  auto normal_of = [&](uint32_t i) {
    std::vector<Vec3> b;
    for (uint32_t q : nb[i]) b.push_back(unit(mi(s.atoms[q].pos - s.atoms[i].pos)));
    Vec3 nn{0, 0, 0};
    if (b.size() == 3) nn = cross(b[0], b[1]) + cross(b[1], b[2]) + cross(b[2], b[0]);
    nn = unit(nn);
    Vec3 r = s.atoms[i].pos - c;
    for (int k = 0; k < 3; ++k) if (spans[k]) r[size_t(k)] = 0;
    const Vec3 rh = unit(r);
    if (norm(nn) < 1e-6) nn = rh;
    double sg = o.side == "inner" ? -1.0 : 1.0;
    if (o.side == "both") sg = U(rng) < 0.5 ? -1.0 : 1.0;
    if (std::fabs(dot(nn, rh)) > 0.3) return nn * (dot(nn, rh) > 0 ? sg : -sg);   // away from (or toward) the centre
    return nn * ((nn[2] >= 0 ? 1.0 : -1.0) * sg);                                    // a flat layer: on top (or below)
  };
  // candidates by pattern
  const std::string pat = o.pattern;
  std::vector<uint32_t> cand;
  if (pat == "atoms") {
    for (uint32_t a : o.atoms) if (a < n && s.atoms[a].element != 1) cand.push_back(a);
  } else if (pat == "ends" || pat == "edges") {
    for (uint32_t i : frame) if (hydrogen_of(i) >= 0) cand.push_back(i);
    if (cand.empty()) throw std::invalid_argument(pat == "ends" ? "the tube has no hydrogen-capped ends (build it finite, not periodic)" : "the sheet has no hydrogen-capped edges (build a flake)");
  } else if (pat == "random" || pat == "all" || pat == "band" || pat == "helix") {
    for (uint32_t i : frame) {
      if (nb[i].size() != 3 || hydrogen_of(i) >= 0) continue;   // free sidewall sites: three bonds, no hydrogen
      if (pat == "band") {
        const double t = hi[size_t(axis)] > lo[size_t(axis)] ? (s.atoms[i].pos[size_t(axis)] - lo[size_t(axis)]) / (hi[size_t(axis)] - lo[size_t(axis)]) : 0;
        if (t < o.from - 1e-9 || t > o.to + 1e-9) continue;
      }
      if (pat == "helix") {
        if (o.pitch <= 0) throw std::invalid_argument("helix: the pitch must be positive");
        const double th = std::atan2(s.atoms[i].pos[1] - c[1], s.atoms[i].pos[0] - c[0]);
        const double want = 2 * kPi * (s.atoms[i].pos[2] - lo[2]) / o.pitch + o.phase * kPi / 180;
        double d = std::remainder(th - want, 2 * kPi);
        if (std::fabs(d) > 20 * kPi / 180) continue;
      }
      cand.push_back(i);
    }
  } else {
    throw std::invalid_argument("unknown pattern '" + pat + "' (random, all, band, helix, ends, edges, atoms)");
  }
  rep.eligible = cand.size();
  if (cand.empty()) throw std::invalid_argument("no site fits the pattern (sidewall sites have three bonds and no hydrogen)");
  if (pat == "random" || pat == "band") caps::shuffle(cand.begin(), cand.end(), rng);
  if (pat == "helix" || pat == "all") std::stable_sort(cand.begin(), cand.end(), [&](uint32_t a, uint32_t b) { return s.atoms[a].pos[size_t(axis)] < s.atoms[b].pos[size_t(axis)]; });
  size_t want = cand.size();
  if (o.count > 0) want = size_t(o.count);
  else if (pat == "random" || pat == "band") want = std::max<size_t>(1, size_t(std::lround(o.fraction * double(cand.size()))));
  // take them in turn, at least min_spacing apart
  std::vector<uint32_t> pick;
  for (uint32_t i : cand) {
    if (pick.size() >= want) break;
    bool ok = true;
    for (uint32_t j : pick) if (norm(mi(s.atoms[i].pos - s.atoms[j].pos)) < o.min_spacing) { ok = false; break; }
    if (ok) pick.push_back(i);
  }
  if (pick.size() < want && o.count > 0)
    rep.notes.push_back(std::to_string(pick.size()) + " of the " + std::to_string(want) + " groups asked for fit at " + std::to_string(o.min_spacing).substr(0, 4) + " Å apart");
  // directions first (the structure grows while grafting), hydrogens replaced at the end
  std::vector<Vec3> dirs;
  std::vector<char> del(n, 0);
  for (uint32_t i : pick) {
    const int h = hydrogen_of(i);
    if (h >= 0) { dirs.push_back(unit(mi(s.atoms[size_t(h)].pos - s.atoms[i].pos))); del[size_t(h)] = 1; }
    else dirs.push_back(normal_of(i));
  }
  const size_t before = s.atoms.size();
  for (size_t k = 0; k < pick.size(); ++k) attach_fragment(s, pick[k], smi, 0, false, &dirs[k]);
  del.resize(s.atoms.size(), 0);
  size_t removed = 0;
  for (char d : del) removed += d;
  if (removed) delete_atoms(s, del);
  rep.grafted = pick.size();
  rep.sites = pick;
  rep.added_atoms = s.atoms.size() + removed - before;
  rep.degree = 100.0 * double(pick.size()) / double(frame.size());
  char b[240];
  std::snprintf(b, sizeof b, "%zu %s group%s on %zu eligible site%s (%s) · %.2f per 100 framework atoms · %zu hydrogens replaced", pick.size(), o.group.c_str(),
                pick.size() == 1 ? "" : "s", cand.size(), cand.size() == 1 ? "" : "s", pat.c_str(), rep.degree, removed);
  rep.notes.insert(rep.notes.begin(), b);
  if (pat != "ends" && pat != "edges") rep.notes.push_back("the grafted atoms keep their sp² places: relax (Minimise) to let them pyramidalise before dynamics");
  return rep;
}

}  // namespace caps
