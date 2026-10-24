// CAPS Exchange, moves (see exchange_moves.hpp).
#include "caps/exchange_moves.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <unordered_map>

namespace caps {

namespace {

constexpr double kB = 0.0019872043;   // kcal/(mol·K)
using Pair = TopologyEditor::Pair;

Pair pair_of(uint32_t i, uint32_t j) { return i < j ? Pair{i, j} : Pair{j, i}; }

}  // namespace

ExchangeKind exchange_kind_from_string(const std::string& s) {
  if (s == "swap") return ExchangeKind::Swap;
  if (s == "shift") return ExchangeKind::Shift;
  if (s == "bind") return ExchangeKind::Bind;
  if (s == "template") return ExchangeKind::Template;
  throw std::invalid_argument("exchange kind must be swap, shift, bind or template (not '" + s + "')");
}

const char* to_string(ExchangeKind k) {
  return k == ExchangeKind::Swap ? "swap" : k == ExchangeKind::Shift ? "shift" : k == ExchangeKind::Bind ? "bind" : "template";
}

ExchangeEngine::ExchangeEngine(System& s, ForceField& ff, const ExchangeOptions& o, const EnergyOptions& eo)
    : s_(s), ff_(ff), o_(o), eo_(eo), ed_(s, ff), rng_(o.seed) {
  const size_t n = s.atoms.size();
  in_a_.assign(n, 0), in_b_.assign(n, 0);
  for (uint32_t a : o.group_a) {
    if (a >= n) throw std::invalid_argument("exchange: group A names an atom the structure does not have");
    in_a_[a] = 1;
  }
  for (uint32_t b : o.group_b) {
    if (b >= n) throw std::invalid_argument("exchange: group B names an atom the structure does not have");
    if (in_a_[b]) throw std::invalid_argument("exchange: an atom cannot be in both groups");
    in_b_[b] = 1;
  }
  if (o.kind == ExchangeKind::Template) {
    if (o.templates.empty()) throw std::invalid_argument("exchange: template moves need a template");
    for (const auto& t : o.templates)
      if (!t.move.empty() || !t.remove.empty())
        throw std::invalid_argument("exchange: template " + t.name + " moves a hydrogen or deletes atoms; run it in LAMMPS (caps exchange bond-react)");
  } else if (o.group_a.empty()) {
    throw std::invalid_argument("exchange: group A is empty");
  }
  if (o.kind == ExchangeKind::Shift && o.group_b.empty()) throw std::invalid_argument("exchange: a shift needs residues (group B)");
  if (!o.labels.empty() && o.labels.size() != n) throw std::invalid_argument("exchange: one label per atom");
  if (!(o.capture > 0)) throw std::invalid_argument("exchange: the capture distance must be positive");
  if (o.valence_a < 1 || o.valence_b < 1) throw std::invalid_argument("exchange: valences start at 1");
  eo_.threads = 1, eo_.mixed = false, eo_.alchemical = nullptr, eo_.parts = 3;
  st_.local = o.local;
  if (o.local && (ff.sw.on || ff.swap3.on() || ff.manybody.on() || ff.hbond.on() || !ff.vsites.empty())) {
    st_.local = false;
    st_.note = "the force field has many-body terms, hydrogen bonds or virtual sites: every move is weighed with the whole system's energy";
  }
  if (st_.local) {
    shell_ = std::make_unique<ForceField>(ff);
    auto& f = *shell_;
    f.atom_type.clear(), f.why.clear(), f.type_index.clear(), f.charge.clear(), f.mass.clear(), f.excluded.clear();
    f.bonds.clear(), f.angles.clear(), f.dihedrals.clear(), f.impropers.clear(), f.impropers_harmonic.clear(), f.inversions.clear();
    f.bonds_x.clear(), f.angles_x.clear(), f.urey_bradley.clear(), f.cbt.clear(), f.lj_pairs.clear();
    f.bonds2.clear(), f.angles2.clear(), f.dihedrals2.clear(), f.impropers2.clear();
    f.pairs14.clear(), f.pairs14_lj.clear(), f.pairs14_coul.clear(), f.impropers_dlpoly.clear(), f.notes.clear(), f.term_sources.clear();
  }
}

ExchangeEngine::~ExchangeEngine() = default;

double ExchangeEngine::distance(const std::vector<double>& x, const Cell& cell, uint32_t i, uint32_t j) const {
  Vec3 d{x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]};
  if (cell.valid()) d = cell.minimum_image(d);
  return norm(d);
}

double ExchangeEngine::cluster_energy(const std::vector<uint32_t>& atoms, const std::vector<double>& x, const Cell& cell) const {
  ForceField sub = *shell_;
  std::unordered_map<uint32_t, uint32_t> id;
  for (uint32_t k = 0; k < atoms.size(); ++k) id[atoms[k]] = k;
  auto in = [&](uint32_t a) { return id.count(a) > 0; };
  std::vector<double> xs;
  for (uint32_t a : atoms) {
    if (!ff_.atom_type.empty()) sub.atom_type.push_back(ff_.atom_type[a]);
    if (ff_.why.size() == ff_.type_index.size()) sub.why.push_back(ff_.why[a]);
    sub.type_index.push_back(ff_.type_index[a]);
    if (!ff_.charge.empty()) sub.charge.push_back(ff_.charge[a]);
    if (!ff_.mass.empty()) sub.mass.push_back(ff_.mass[a]);
    std::vector<uint32_t> ex;
    if (a < ff_.excluded.size())
      for (uint32_t e : ff_.excluded[a])
        if (in(e)) ex.push_back(id[e]);
    std::sort(ex.begin(), ex.end());
    sub.excluded.push_back(std::move(ex));
    xs.insert(xs.end(), {x[3 * a], x[3 * a + 1], x[3 * a + 2]});
  }
  // the terms with every atom in the cluster, renumbered
  auto keep2 = [&](const auto& from, auto& to) {
    for (auto t : from)
      if (in(t.i) && in(t.j)) t.i = id[t.i], t.j = id[t.j], to.push_back(t);
  };
  auto keep3 = [&](const auto& from, auto& to) {
    for (auto t : from)
      if (in(t.i) && in(t.j) && in(t.k)) t.i = id[t.i], t.j = id[t.j], t.k = id[t.k], to.push_back(t);
  };
  auto keep4 = [&](const auto& from, auto& to) {
    for (auto t : from)
      if (in(t.i) && in(t.j) && in(t.k) && in(t.l)) t.i = id[t.i], t.j = id[t.j], t.k = id[t.k], t.l = id[t.l], to.push_back(t);
  };
  keep2(ff_.bonds, sub.bonds), keep2(ff_.bonds_x, sub.bonds_x), keep2(ff_.bonds2, sub.bonds2), keep2(ff_.lj_pairs, sub.lj_pairs);
  keep3(ff_.angles, sub.angles), keep3(ff_.angles_x, sub.angles_x), keep3(ff_.angles2, sub.angles2);
  for (auto t : ff_.urey_bradley)
    if (in(t.i) && in(t.k)) t.i = id[t.i], t.k = id[t.k], sub.urey_bradley.push_back(t);
  keep4(ff_.dihedrals, sub.dihedrals), keep4(ff_.impropers, sub.impropers), keep4(ff_.impropers_harmonic, sub.impropers_harmonic);
  keep4(ff_.dihedrals2, sub.dihedrals2), keep4(ff_.impropers2, sub.impropers2), keep4(ff_.cbt, sub.cbt);
  for (auto v : ff_.inversions)
    if (in(v.c) && in(v.a) && in(v.b) && in(v.d)) v.c = id[v.c], v.a = id[v.a], v.b = id[v.b], v.d = id[v.d], sub.inversions.push_back(v);
  const bool per = ff_.per_pair14();
  for (size_t k = 0; k < ff_.pairs14.size(); ++k) {
    const auto& p = ff_.pairs14[k];
    if (!in(p[0]) || !in(p[1])) continue;
    sub.pairs14.push_back({id[p[0]], id[p[1]]});
    if (per) sub.pairs14_lj.push_back(ff_.pairs14_lj[k]), sub.pairs14_coul.push_back(k < ff_.pairs14_coul.size() ? ff_.pairs14_coul[k] : ff_.coul14);
  }
  Evaluator ev(sub, eo_);
  std::vector<double> f;
  return ev.compute(xs, cell, f).total();
}

double ExchangeEngine::apply(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::vector<double>& x, const Cell& cell,
                             const std::function<double()>& potential) {
  if (!st_.local) {
    if (!potential) throw std::invalid_argument("exchange: the whole system's energy is needed (no local energies for this force field)");
    const double e0 = potential();
    ed_.change(broken, made);
    return potential() - e0;
  }
  // the cluster: atoms within three bonds of the change, in the graph before it together with the bonds it makes
  std::set<uint32_t> seen, front;
  for (const auto& p : broken) front.insert(p[0]), front.insert(p[1]);
  for (const auto& p : made) front.insert(p[0]), front.insert(p[1]);
  seen = front;
  const auto& nb = ed_.neighbours();
  for (int depth = 0; depth < 3; ++depth) {
    std::set<uint32_t> next;
    for (uint32_t a : front) {
      for (uint32_t b : nb[a])
        if (!seen.count(b)) next.insert(b);
      for (const auto& p : made)
        for (int k = 0; k < 2; ++k)
          if (p[k] == a && !seen.count(p[1 - k])) next.insert(p[1 - k]);
    }
    for (uint32_t b : next) seen.insert(b);
    front = std::move(next);
  }
  const std::vector<uint32_t> atoms(seen.begin(), seen.end());
  const double e0 = cluster_energy(atoms, x, cell);
  ed_.change(broken, made);
  return cluster_energy(atoms, x, cell) - e0;
}

bool ExchangeEngine::decide(const std::vector<Pair>& broken, const std::vector<Pair>& made, double log_ratio, double beta,
                            const std::vector<double>& x, const Cell& cell, const std::function<double()>& potential) {
  double du;
  try {
    du = apply(broken, made, x, cell, potential);
  } catch (const TopologyError&) {
    ++st_.unparameterised;
    return false;
  }
  const double l = log_ratio - beta * du;
  const double acc = std::exp(-beta * o_.barrier) * (l >= 0 ? 1.0 : std::exp(l));
  if (std::uniform_real_distribution<double>(0, 1)(rng_) < acc) {
    ++st_.accepted;
    st_.energy += du;
    return true;
  }
  ed_.change(made, broken);   // undone: the editor restores the broken bonds' own parameters
  return false;
}

double ExchangeEngine::delta(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::vector<double>& x, const Cell& cell,
                             const std::function<double()>& potential) {
  const double du = apply(broken, made, x, cell, potential);
  ed_.change(made, broken);
  return du;
}

bool ExchangeEngine::try_swap(const std::vector<double>& x, const Cell& cell, double beta, const std::function<double()>& potential) {
  ++st_.attempts;
  auto u = [&] { return std::uniform_real_distribution<double>(0, 1)(rng_); };
  auto dynamic_bonds = [&] {   // bonds with both atoms in group A
    std::vector<Pair> d;
    const auto& nb = ed_.neighbours();
    for (uint32_t i = 0; i < nb.size(); ++i)
      if (in_a_[i])
        for (uint32_t j : nb[i])
          if (i < j && in_a_[j]) d.push_back({i, j});
    return d;
  };
  auto dyn = dynamic_bonds();
  if (dyn.size() < 2) { ++st_.no_candidate; return false; }
  const Pair pick = dyn[std::min(dyn.size() - 1, size_t(u() * double(dyn.size())))];
  const bool flip = u() < 0.5;
  const uint32_t a1 = flip ? pick[1] : pick[0], a2 = flip ? pick[0] : pick[1];
  const double cap = o_.capture;
  // the swap x1–x2 + y1–y2 → x1–y2 + y1–x2 is allowed in the current topology
  auto valid = [&](uint32_t x1, uint32_t x2, uint32_t y1, uint32_t y2) {
    if (y1 == x1 || y1 == x2 || y2 == x1 || y2 == x2) return false;
    if (!ed_.bonded(x1, x2) || !ed_.bonded(y1, y2)) return false;
    if (!o_.labels.empty() && (o_.labels[x1] != o_.labels[y1] || o_.labels[x2] != o_.labels[y2])) return false;   // as LAMMPS: both ends
    if (ed_.bonded(x1, y2) || ed_.bonded(y1, x2)) return false;
    return distance(x, cell, x1, x2) < cap && distance(x, cell, y1, y2) < cap && distance(x, cell, x1, y2) < cap && distance(x, cell, y1, x2) < cap;
  };
  auto candidates = [&](uint32_t x1, uint32_t x2) {
    std::vector<Pair> c;
    for (const auto& b : dyn)
      for (int o = 0; o < 2; ++o) {
        const uint32_t y1 = o ? b[1] : b[0], y2 = o ? b[0] : b[1];
        if (valid(x1, x2, y1, y2)) c.push_back({y1, y2});
      }
    return c;
  };
  const auto cand = candidates(a1, a2);
  if (cand.empty()) { ++st_.no_candidate; return false; }
  const Pair c = cand[std::min(cand.size() - 1, size_t(u() * double(cand.size())))];
  const uint32_t b1 = c[0], b2 = c[1];
  // proposal weights: the four ways of picking the same swap (each bond first, each way round), before it and, for the
  // way back, after it; the bond count and the 1/2 for the direction are the same both ways
  auto weight = [&](uint32_t x1, uint32_t x2, uint32_t y1, uint32_t y2) { return valid(x1, x2, y1, y2) ? 1.0 / double(candidates(x1, x2).size()) : 0.0; };
  const double qf = weight(a1, a2, b1, b2) + weight(b1, b2, a1, a2) + weight(a2, a1, b2, b1) + weight(b2, b1, a2, a1);
  ed_.swap(a1, a2, b1, b2);
  dyn = dynamic_bonds();
  const double qr = weight(a1, b2, b1, a2) + weight(b1, a2, a1, b2) + weight(b2, a1, a2, b1) + weight(a2, b1, b2, a1);
  ed_.swap(a1, b2, b1, a2);
  if (!(qf > 0) || !(qr > 0)) { ++st_.no_candidate; return false; }
  return decide({pair_of(a1, a2), pair_of(b1, b2)}, {pair_of(a1, b2), pair_of(b1, a2)}, std::log(qr / qf), beta, x, cell, potential);
}

bool ExchangeEngine::try_shift(const std::vector<double>& x, const Cell& cell, double beta, const std::function<double()>& potential) {
  ++st_.attempts;
  auto u = [&] { return std::uniform_real_distribution<double>(0, 1)(rng_); };
  auto residues_of = [&](uint32_t p) {
    std::vector<uint32_t> r;
    for (uint32_t b : ed_.neighbours()[p])
      if (in_b_[b]) r.push_back(b);
    return r;
  };
  auto load = [&](uint32_t r) {   // bonds a residue holds to pivots
    int k = 0;
    for (uint32_t b : ed_.neighbours()[r]) k += in_a_[b];
    return k;
  };
  auto attackers = [&](uint32_t p) {
    std::vector<uint32_t> a;
    for (uint32_t r : o_.group_b)
      if (!ed_.bonded(p, r) && load(r) < o_.valence_b && distance(x, cell, p, r) < o_.capture) a.push_back(r);
    return a;
  };
  std::vector<uint32_t> pivots;   // pivots holding a residue: a shift keeps their number
  for (uint32_t p : o_.group_a)
    if (!residues_of(p).empty()) pivots.push_back(p);
  if (pivots.empty()) { ++st_.no_candidate; return false; }
  const uint32_t p = pivots[std::min(pivots.size() - 1, size_t(u() * double(pivots.size())))];
  const auto lv = residues_of(p);
  const uint32_t r0 = lv[std::min(lv.size() - 1, size_t(u() * double(lv.size())))];
  const auto att = attackers(p);
  if (att.empty()) { ++st_.no_candidate; return false; }
  const uint32_t r1 = att[std::min(att.size() - 1, size_t(u() * double(att.size())))];
  // the way back: the same pivot and leaving count; its attackers after the move
  try {
    ed_.change({pair_of(p, r0)}, {pair_of(p, r1)});
  } catch (const TopologyError&) {
    ++st_.unparameterised;
    return false;
  }
  const double back = double(attackers(p).size());
  ed_.change({pair_of(p, r1)}, {pair_of(p, r0)});
  return decide({pair_of(p, r0)}, {pair_of(p, r1)}, std::log(double(att.size()) / back), beta, x, cell, potential);
}

bool ExchangeEngine::try_template(const std::vector<double>& x, const Cell& cell, double beta, const std::function<double()>& potential) {
  ++st_.attempts;
  auto u = [&] { return std::uniform_real_distribution<double>(0, 1)(rng_); };
  const auto& t = o_.templates[std::min(o_.templates.size() - 1, size_t(u() * double(o_.templates.size())))];
  // the matches at the run's positions (the structure's positions are the run's: the dynamics callback copies them in)
  for (size_t i = 0; i < s_.atoms.size(); ++i) s_.atoms[i].pos = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
  const auto matches = find_matches(s_, t);
  if (matches.empty()) { ++st_.no_candidate; return false; }
  const auto& m = matches[std::min(matches.size() - 1, size_t(u() * double(matches.size())))];
  auto at = [&](int map) {
    for (size_t q = 0; q < t.atoms.size(); ++q)
      if (t.atoms[q].map == map) return m.atoms[q];
    throw std::invalid_argument("exchange: template " + t.name + " names map atom " + std::to_string(map) + " it does not have");
  };
  std::vector<Pair> broken, made;
  for (const auto& b : t.brk) broken.push_back(pair_of(at(b.first), at(b.second)));
  for (const auto& f : t.form) made.push_back(pair_of(at(f.first), at(f.second)));
  // the way back is the same template on the products: the ratio of the matches before and after
  try {
    ed_.change(broken, made);
  } catch (const TopologyError&) {
    ++st_.unparameterised;
    return false;
  }
  const double back = double(find_matches(s_, t).size());
  ed_.change(made, broken);
  if (!(back > 0)) { ++st_.no_candidate; return false; }
  return decide(broken, made, std::log(double(matches.size()) / back), beta, x, cell, potential);
}

bool ExchangeEngine::try_bind(const std::vector<double>& x, const Cell& cell, double beta, const std::function<double()>& potential) {
  ++st_.attempts;
  auto u = [&] { return std::uniform_real_distribution<double>(0, 1)(rng_); };
  const bool self = o_.group_b.empty();
  std::vector<uint32_t> sites = o_.group_a;
  sites.insert(sites.end(), o_.group_b.begin(), o_.group_b.end());
  const double ns = double(sites.size());
  auto complementary = [&](uint32_t i, uint32_t j) { return self ? (in_a_[i] && in_a_[j]) : ((in_a_[i] && in_b_[j]) || (in_b_[i] && in_a_[j])); };
  auto valence = [&](uint32_t i) { return in_a_[i] ? o_.valence_a : o_.valence_b; };
  auto dynamic = [&](uint32_t i) {
    std::vector<uint32_t> d;
    for (uint32_t j : ed_.neighbours()[i])
      if (complementary(i, j)) d.push_back(j);
    return d;
  };
  auto free_partners = [&](uint32_t i) {
    std::vector<uint32_t> c;
    for (uint32_t j : sites)
      if (j != i && complementary(i, j) && !ed_.bonded(i, j) && int(dynamic(j).size()) < valence(j) && distance(x, cell, i, j) < o_.capture) c.push_back(j);
    return c;
  };
  // the probability of proposing to bind i and j, and to break them: either site picked first, then bind or break (½ each)
  auto q_bind = [&](uint32_t i, uint32_t j) {
    double q = 0;
    for (const auto& [a, b] : {std::pair<uint32_t, uint32_t>{i, j}, {j, i}}) {
      if (int(dynamic(a).size()) >= valence(a)) continue;
      const auto c = free_partners(a);
      if (std::find(c.begin(), c.end(), b) != c.end()) q += 1.0 / double(c.size());
    }
    return q / (2 * ns);
  };
  auto q_break = [&](uint32_t i, uint32_t j) { return (1.0 / double(dynamic(i).size()) + 1.0 / double(dynamic(j).size())) / (2 * ns); };
  const uint32_t i = sites[std::min(sites.size() - 1, size_t(u() * ns))];
  const double lk = beta * o_.binding;   // ln K
  if (u() < 0.5) {   // bind
    if (int(dynamic(i).size()) >= valence(i)) { ++st_.no_candidate; return false; }
    const auto c = free_partners(i);
    if (c.empty()) { ++st_.no_candidate; return false; }
    const uint32_t j = c[std::min(c.size() - 1, size_t(u() * double(c.size())))];
    const double qf = q_bind(i, j);
    try {
      ed_.change({}, {pair_of(i, j)});
    } catch (const TopologyError&) {
      ++st_.unparameterised;
      return false;
    }
    const double qr = q_break(i, j);
    ed_.change({pair_of(i, j)}, {});
    return decide({}, {pair_of(i, j)}, std::log(qr / qf) + lk, beta, x, cell, potential);
  }
  const auto d = dynamic(i);   // break
  if (d.empty()) { ++st_.no_candidate; return false; }
  const uint32_t j = d[std::min(d.size() - 1, size_t(u() * double(d.size())))];
  const double qf = q_break(i, j);
  ed_.change({pair_of(i, j)}, {});
  const double qr = q_bind(i, j);
  ed_.change({}, {pair_of(i, j)});
  if (!(qr > 0)) { ++st_.no_candidate; return false; }
  return decide({pair_of(i, j)}, {}, std::log(qr / qf) - lk, beta, x, cell, potential);
}

bool ExchangeEngine::run(const std::vector<double>& x, const Cell& cell, double temperature, const std::function<double()>& potential) {
  if (!(temperature > 0)) throw std::invalid_argument("exchange: the temperature must be positive");
  const double beta = 1.0 / (kB * temperature);
  bool any = false;
  for (int t = 0; t < o_.attempts; ++t) {
    const bool kept = o_.kind == ExchangeKind::Swap       ? try_swap(x, cell, beta, potential)
                      : o_.kind == ExchangeKind::Shift    ? try_shift(x, cell, beta, potential)
                      : o_.kind == ExchangeKind::Template ? try_template(x, cell, beta, potential)
                                                          : try_bind(x, cell, beta, potential);
    any = any || kept;
  }
  return any;
}

std::function<bool(TopologyContext&)> exchange_callback(const ExchangeOptions& o, const EnergyOptions& eo, std::shared_ptr<ExchangeStats> stats) {
  auto engine = std::make_shared<std::unique_ptr<ExchangeEngine>>();
  return [o, eo, stats, engine](TopologyContext& c) {
    if (!*engine) *engine = std::make_unique<ExchangeEngine>(c.system, c.field, o, eo);
    const bool changed = (*engine)->run(c.x, c.cell, c.temperature, c.potential);
    if (stats) *stats = (*engine)->stats();
    return changed;
  };
}

std::vector<int64_t> chain_position_labels(const System& s, bool from_both_ends) {
  const size_t n = s.atoms.size();
  std::vector<std::vector<uint32_t>> nb(n);
  for (const auto& b : s.bonds)
    if (b.i < n && b.j < n && b.i != b.j) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  std::vector<int64_t> label(n, 0);
  std::vector<char> done(n, 0);
  for (uint32_t start = 0; start < n; ++start) {
    if (done[start]) continue;
    // the component, and whether it is a linear chain
    std::vector<uint32_t> comp{start}, stack{start};
    done[start] = 1;
    while (!stack.empty()) {
      const uint32_t a = stack.back();
      stack.pop_back();
      for (uint32_t b : nb[a])
        if (!done[b]) done[b] = 1, comp.push_back(b), stack.push_back(b);
    }
    std::vector<uint32_t> ends;
    bool linear = true;
    for (uint32_t a : comp) {
      if (nb[a].size() > 2) linear = false;
      if (nb[a].size() <= 1) ends.push_back(a);
    }
    if (!linear || (comp.size() > 1 && ends.size() != 2)) continue;
    uint32_t cur = *std::min_element(ends.begin(), ends.end()), prev = cur;
    const int64_t len = int64_t(comp.size());
    for (int64_t k = 1; k <= len; ++k) {
      label[cur] = from_both_ends ? std::min(k, len + 1 - k) : k;
      uint32_t next = cur;
      for (uint32_t b : nb[cur])
        if (b != prev) { next = b; break; }
      if (next == cur) break;
      prev = cur, cur = next;
    }
  }
  return label;
}

}  // namespace caps
