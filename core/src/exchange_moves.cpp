// CAPS Exchange, moves (see exchange_moves.hpp).
#include "caps/exchange_moves.hpp"

#include "caps/elements.hpp"

#include <algorithm>
#include <cmath>
#include <map>
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
    for (const auto& t : o.templates) {
      if (!t.remove.empty()) throw std::invalid_argument("exchange: template " + t.name + " deletes atoms; run it in LAMMPS (caps exchange bond-react)");
      if (!t.move.empty() && o.relax_steps <= 0)
        throw std::invalid_argument("exchange: template " + t.name + " moves a hydrogen, which needs relaxed moves (relax_steps > 0; --relax N)");
    }
    for (auto& t : o_.templates)
      if (t.min_path <= 0) t.min_path = 4;   // any molecule, the network's own included
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
  if (o.relax_steps > 0 && !st_.local)
    throw std::invalid_argument("exchange: relaxed moves need local energies (" + (st_.note.empty() ? std::string("local is off") : st_.note) + ")");
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

void ExchangeEngine::edit(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::vector<TopologyEditor::Retype>& retype) {
  if (!o_.parameterize) {
    ed_.change(broken, made, retype);
    return;
  }
  // the neighbourhood after the change, six bonds deep, typed again by the run's force field; atoms within four bonds of
  // the change have every neighbour in it, so their types and terms are the structure's own. The atoms within three bonds
  // take those types (what a fresh assignment gives them: a secondary carbon that becomes an ester's stays a CH), and
  // the terms of combinations the structure never had are learned from it
  const size_t n = s_.atoms.size();
  auto nb = ed_.neighbours();
  for (const auto& p : broken) {
    nb[p[0]].erase(std::remove(nb[p[0]].begin(), nb[p[0]].end(), p[1]), nb[p[0]].end());
    nb[p[1]].erase(std::remove(nb[p[1]].begin(), nb[p[1]].end(), p[0]), nb[p[1]].end());
  }
  for (const auto& p : made) nb[p[0]].push_back(p[1]), nb[p[1]].push_back(p[0]);
  std::vector<int> depth(n, -1);
  std::vector<uint32_t> order;
  for (const auto& v : {broken, made})
    for (const auto& p : v)
      for (uint32_t a : p)
        if (depth[a] < 0) depth[a] = 0, order.push_back(a);
  for (const auto& r : retype)
    if (depth[r.atom] < 0) depth[r.atom] = 0, order.push_back(r.atom);
  std::map<Pair, int> bond_order;
  for (const auto& b : s_.bonds) bond_order[pair_of(b.i, b.j)] = b.order;
  // past six bonds only through double and aromatic bonds: a cut between two unsaturated atoms would not type
  for (size_t h = 0; h < order.size(); ++h) {
    const uint32_t a = order[h];
    for (uint32_t b : nb[a]) {
      if (depth[b] >= 0) continue;
      const auto it = bond_order.find(pair_of(a, b));
      if (depth[a] < 6 || (it != bond_order.end() && it->second != 1)) depth[b] = depth[a] + 1, order.push_back(b);
    }
  }
  std::unordered_map<uint32_t, uint32_t> id;
  System frag;
  frag.cell = s_.cell;
  std::vector<char> complete;
  for (uint32_t a : order) {
    id[a] = uint32_t(frag.atoms.size());
    frag.atoms.push_back(s_.atoms[a]);
    complete.push_back(depth[a] <= 4);
  }
  for (uint32_t a : order)
    for (uint32_t b : nb[a]) {
      if (a < b && id.count(b)) {
        const auto it = bond_order.find(pair_of(a, b));
        frag.bonds.push_back({id[a], id[b], it != bond_order.end() ? it->second : 1});
      } else if (!id.count(b)) {
        // a cut bond: one hydrogen cap at 1.09 Å towards the atom left out, so the edge still types (as ff complete cuts)
        Vec3 d = s_.atoms[b].pos - s_.atoms[a].pos;
        if (s_.cell.valid()) d = s_.cell.minimum_image(d);
        const double len = norm(d);
        Atom h;
        h.element = 1;
        h.mol = s_.atoms[a].mol;
        h.pos = s_.atoms[a].pos + (len > 1e-6 ? d * (1.09 / len) : Vec3{1.09, 0, 0});
        frag.bonds.push_back({id[a], uint32_t(frag.atoms.size()), 1});
        frag.atoms.push_back(h);
        complete.push_back(0);
      }
    }
  ForceField typed;
  try {
    typed = o_.parameterize(frag);
  } catch (const std::exception& e) {
    // the assignment could not type the neighbourhood (or fill what it lacks): the change stays one it cannot describe
    if (st_.note.empty()) st_.note = std::string("terms for a new combination of types could not be learned: ") + e.what();
    throw TopologyError(e.what());
  }
  // a type the run's force field has for no atom yet (an ester's tertiary carbon that becomes a tertiary alcohol's)
  // joins it with its non-bonded parameters from the same assignment
  std::map<std::string, int> ours;
  for (size_t k = 0; k < ff_.type_names.size(); ++k) ours.emplace(ff_.type_names[k], int(k));
  for (size_t k = 0; k < order.size(); ++k) {
    if (depth[order[k]] > 3) continue;
    const int fi = typed.type_index[k];
    if (fi < 0 || size_t(fi) >= typed.type_names.size() || size_t(fi) >= typed.lj.size()) throw TopologyError("the run's assignment left an atom untyped");
    const std::string& name = typed.type_names[fi];
    if (ours.count(name)) continue;
    if (!ff_.pair_override.empty() || !ff_.pair_func.empty() || (!ff_.lj14_types.empty() && size_t(fi) >= typed.lj14_types.size()))
      throw TopologyError("the change makes type " + name + ", new to the run, and the force field gives pairs by type pair");
    for (ForceField* f : {&ff_, shell_.get()}) {
      if (!f) continue;
      f->type_names.push_back(name);
      f->lj.push_back(typed.lj[fi]);
      if (!f->lj14_types.empty()) f->lj14_types.push_back(typed.lj14_types[fi]);
      if (!f->type_part.empty()) f->type_part.push_back(f->type_part[size_t(ff_.type_index[order[k]])]);
    }
    ours.emplace(name, int(ff_.type_names.size()) - 1);
    ++st_.new_types;
  }
  st_.learned += int64_t(ed_.learn(frag, typed, complete));
  std::map<uint32_t, double> role_charge;   // charges follow roles, so the total stays
  for (const auto& r : retype) role_charge[r.atom] = r.charge;
  std::vector<TopologyEditor::Retype> rt;
  unretype_.clear();
  for (size_t k = 0; k < order.size(); ++k) {
    const uint32_t a = order[k];
    if (depth[a] > 3) continue;
    const int fi = typed.type_index[id[a]];
    const std::string name = fi >= 0 && size_t(fi) < typed.type_names.size() ? typed.type_names[fi] : std::string();
    const auto it = ours.find(name);
    const auto rc = role_charge.find(a);
    const double q = rc != role_charge.end() ? rc->second : a < ff_.charge.size() ? ff_.charge[a] : 0.0;
    if (it->second == ff_.type_index[a] && (ff_.charge.empty() || q == ff_.charge[a])) continue;
    rt.push_back({a, it->second, name, q});
    unretype_.push_back({a, ff_.type_index[a], a < ff_.atom_type.size() ? ff_.atom_type[a] : std::string(), a < ff_.charge.size() ? ff_.charge[a] : 0.0});
  }
  ed_.change(broken, made, rt);
}

double ExchangeEngine::distance(const std::vector<double>& x, const Cell& cell, uint32_t i, uint32_t j) const {
  Vec3 d{x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]};
  if (cell.valid()) d = cell.minimum_image(d);
  return norm(d);
}

double ExchangeEngine::cluster_energy(const std::vector<uint32_t>& atoms, const std::vector<double>& x, const Cell& cell) const {
  const ForceField sub = cluster_field(atoms);
  std::vector<double> xs;
  for (uint32_t a : atoms) xs.insert(xs.end(), {x[3 * a], x[3 * a + 1], x[3 * a + 2]});
  Evaluator ev(sub, eo_);
  std::vector<double> f;
  return ev.compute(xs, cell, f).total();
}

ForceField ExchangeEngine::cluster_field(const std::vector<uint32_t>& atoms) const {
  ForceField sub = *shell_;
  std::unordered_map<uint32_t, uint32_t> id;
  for (uint32_t k = 0; k < atoms.size(); ++k) id[atoms[k]] = k;
  auto in = [&](uint32_t a) { return id.count(a) > 0; };
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
  return sub;
}

double ExchangeEngine::apply_relaxed(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::vector<double>& x, const Cell& cell) {
  const auto& nb = ed_.neighbours();
  const size_t n = s_.atoms.size();
  // movable: the atoms of the change, the atoms bonded to them (before and after) and the hydrogens on those
  std::set<uint32_t> mov;
  for (const auto& p : broken) mov.insert(p[0]), mov.insert(p[1]);
  for (const auto& p : made) mov.insert(p[0]), mov.insert(p[1]);
  for (int d = 0; d < o_.relax_depth; ++d)
    for (uint32_t a : std::vector<uint32_t>(mov.begin(), mov.end()))
      for (uint32_t b : nb[a]) mov.insert(b);
  for (uint32_t a : std::vector<uint32_t>(mov.begin(), mov.end()))
    for (uint32_t b : nb[a])
      if (s_.atoms[b].element == 1) mov.insert(b);
  // held: every atom within three bonds of a movable one (each bonded term holding a movable atom), and every atom within
  // the cut-off of one (its non-bonded partners)
  std::set<uint32_t> seen = mov, front = mov;
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
  const double rc = (ff_.cutoff > 0 ? ff_.cutoff : eo_.cutoff) + 0.5;
  for (uint32_t i = 0; i < n; ++i) {
    if (seen.count(i)) continue;
    for (uint32_t a : mov)
      if (distance(x, cell, a, i) < rc) { seen.insert(i); break; }
  }
  const std::vector<uint32_t> atoms(seen.begin(), seen.end());
  std::unordered_map<uint32_t, uint32_t> id;
  std::vector<double> xs;
  for (uint32_t k = 0; k < atoms.size(); ++k) id[atoms[k]] = k, xs.insert(xs.end(), {x[3 * atoms[k]], x[3 * atoms[k] + 1], x[3 * atoms[k] + 2]});
  std::vector<char> movable(atoms.size(), 0);
  for (uint32_t a : mov) movable[id.at(a)] = 1;
  // steepest descent on the movable atoms, the step grown on success and halved on failure, at most 0.2 Å per step
  auto relax = [&](const ForceField& sub, std::vector<double>& q) {
    Evaluator ev(sub, eo_);
    std::vector<double> f, ft, trial;
    double e = ev.compute(q, cell, f).total(), step = 0.02;
    for (int it = 0; it < o_.relax_steps; ++it) {
      double fmax = 0;
      for (size_t k = 0; k < atoms.size(); ++k)
        if (movable[k]) fmax = std::max(fmax, std::sqrt(f[3 * k] * f[3 * k] + f[3 * k + 1] * f[3 * k + 1] + f[3 * k + 2] * f[3 * k + 2]));
      if (fmax < 0.1) break;
      trial = q;
      const double s = std::min(step, 0.2) / fmax;
      for (size_t k = 0; k < atoms.size(); ++k)
        if (movable[k])
          for (int c = 0; c < 3; ++c) trial[3 * k + c] += s * f[3 * k + c];
      const double et = ev.compute(trial, cell, ft).total();
      if (et < e) q.swap(trial), f.swap(ft), e = et, step *= 1.2;
      else step *= 0.5;
    }
    return e;
  };
  // both sides relaxed alike, so ΔU compares two minima and not a minimum with a thermal state
  std::vector<double> before = xs;
  const double e0 = relax(cluster_field(atoms), before);
  edit(broken, made, retype_);
  // a moved hydrogen at a covalent bond length from its new partner, on the side it came from (React's placement)
  for (const auto& [h, p] : moves_) {
    const uint32_t hi = id.at(h), pi = id.at(p);
    Vec3 d{xs[3 * hi] - xs[3 * pi], xs[3 * hi + 1] - xs[3 * pi + 1], xs[3 * hi + 2] - xs[3 * pi + 2]};
    if (cell.valid()) d = cell.minimum_image(d);
    double len = norm(d);
    if (len < 1e-6) d = {1, 0, 0}, len = 1;
    const double r0 = element(s_.atoms[h].element).covalent + element(s_.atoms[p].element).covalent;
    for (int c = 0; c < 3; ++c) xs[3 * hi + c] = xs[3 * pi + c] + d[c] * (r0 / len);
  }
  const double e = relax(cluster_field(atoms), xs);
  pending_.clear();
  for (uint32_t a : mov) {
    const uint32_t k = id.at(a);
    pending_.push_back({a, Vec3{xs[3 * k], xs[3 * k + 1], xs[3 * k + 2]}});
  }
  ++st_.relaxed;
  return e - e0;
}

double ExchangeEngine::apply(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::vector<double>& x, const Cell& cell,
                             const std::function<double()>& potential) {
  pending_.clear();
  if (o_.relax_steps > 0) return apply_relaxed(broken, made, x, cell);
  if (!st_.local) {
    if (!potential) throw std::invalid_argument("exchange: the whole system's energy is needed (no local energies for this force field)");
    const double e0 = potential();
    edit(broken, made, retype_);
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
  edit(broken, made, retype_);
  return cluster_energy(atoms, x, cell) - e0;
}

bool ExchangeEngine::decide(const std::vector<Pair>& broken, const std::vector<Pair>& made, double log_ratio, double beta,
                            const std::vector<double>& x, const Cell& cell, const std::function<double()>& potential) {
  double du;
  try {
    du = apply(broken, made, x, cell, potential);
  } catch (const TopologyError& e) {
    if (st_.note.empty()) st_.note = std::string("a change the force field cannot describe: ") + e.what();
    ++st_.unparameterised;
    return false;
  }
  const double l = log_ratio - beta * du;
  const double acc = std::exp(-beta * o_.barrier) * (l >= 0 ? 1.0 : std::exp(l));
  if (std::uniform_real_distribution<double>(0, 1)(rng_) < acc) {
    ++st_.accepted;
    st_.energy += du;
    if (!pending_.empty()) {
      if (!xw_) throw std::logic_error("exchange: a relaxed move needs the run's positions to write to (run with a mutable x)");
      for (const auto& [a, r] : pending_)
        for (int c = 0; c < 3; ++c) (*xw_)[3 * a + c] = r[c];
      pending_.clear();
    }
    return true;
  }
  pending_.clear();
  ed_.change(made, broken, unretype_);   // undone: the editor restores the broken bonds' own parameters
  return false;
}

double ExchangeEngine::delta(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::vector<double>& x, const Cell& cell,
                             const std::function<double()>& potential) {
  const double du = apply(broken, made, x, cell, potential);
  ed_.change(made, broken, unretype_);
  pending_.clear();
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
  // a moved hydrogen leaves every partner it has for its new one
  moves_.clear();
  for (const auto& [h, p] : t.move) {
    const uint32_t hi = at(h), pi = at(p);
    for (uint32_t q : ed_.neighbours()[hi])
      if (q != pi) broken.push_back(pair_of(hi, q));
    made.push_back(pair_of(hi, pi));
    moves_.push_back({hi, pi});
  }
  for (auto* v : {&broken, &made}) std::sort(v->begin(), v->end()), v->erase(std::unique(v->begin(), v->end()), v->end());
  // the way back is the same template on the products: the ratio of the matches before and after
  System prod = s_;
  auto has = [](const std::vector<Pair>& v, uint32_t i, uint32_t j) { return std::find(v.begin(), v.end(), pair_of(i, j)) != v.end(); };
  prod.bonds.erase(std::remove_if(prod.bonds.begin(), prod.bonds.end(), [&](const Bond& b) { return has(broken, b.i, b.j); }), prod.bonds.end());
  for (const auto& p : made) prod.bonds.push_back({p[0], p[1], 1});
  const auto backs = find_matches(prod, t);
  if (backs.empty()) { ++st_.no_candidate; return false; }
  // types follow roles: the back match that undoes this one names each role's atom after the change
  std::set<Pair> heavy_broken;
  for (const auto& b : t.brk) heavy_broken.insert(pair_of(at(b.first), at(b.second)));
  std::vector<std::vector<uint32_t>> after(prod.atoms.size());
  for (const auto& b : prod.bonds) after[b.i].push_back(b.j), after[b.j].push_back(b.i);
  retype_.clear(), unretype_.clear();
  for (const auto& r : backs) {
    std::set<Pair> undo;
    for (const auto& f : t.form) undo.insert(pair_of(r.atoms[size_t(std::find_if(t.atoms.begin(), t.atoms.end(), [&](const TemplateAtom& a) { return a.map == f.first; }) - t.atoms.begin())],
                                                      r.atoms[size_t(std::find_if(t.atoms.begin(), t.atoms.end(), [&](const TemplateAtom& a) { return a.map == f.second; }) - t.atoms.begin())]));
    if (undo != heavy_broken) continue;
    std::map<uint32_t, uint32_t> from;   // atom after the change → the atom whose type and charge it takes
    std::map<uint32_t, size_t> role_before, role_after;
    for (size_t q = 0; q < m.atoms.size(); ++q) role_before[m.atoms[q]] = q, role_after[r.atoms[q]] = q;
    for (size_t q = 0; q < m.atoms.size(); ++q) from[r.atoms[q]] = m.atoms[q];
    // atoms that leave the pattern take the types of atoms that join it, paired through their neighbours' roles
    std::vector<uint32_t> gained;
    for (uint32_t b : r.atoms)
      if (!role_before.count(b)) gained.push_back(b);
    std::vector<char> used(gained.size(), 0);
    for (uint32_t a : m.atoms) {
      if (role_after.count(a)) continue;
      for (size_t g = 0; g < gained.size() && !from.count(a); ++g) {
        if (used[g] || s_.atoms[gained[g]].element != s_.atoms[a].element) continue;
        bool linked = false;
        for (uint32_t c : after[a])
          if (const auto ra = role_after.find(c); ra != role_after.end())
            for (uint32_t d : ed_.neighbours()[gained[g]]) linked = linked || d == m.atoms[ra->second];
        if (linked) from[a] = gained[g], used[g] = 1;
      }
    }
    for (const auto& [a, b] : from) {
      if (a == b || (ff_.type_index[a] == ff_.type_index[b] && (ff_.charge.empty() || ff_.charge[a] == ff_.charge[b]))) continue;
      auto snap = [&](uint32_t x, uint32_t as) {
        return TopologyEditor::Retype{as, ff_.type_index[x], x < ff_.atom_type.size() ? ff_.atom_type[x] : std::string(), x < ff_.charge.size() ? ff_.charge[x] : 0.0};
      };
      retype_.push_back(snap(b, a)), unretype_.push_back(snap(a, a));
    }
    break;
  }
  const double back = double(backs.size());
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

bool ExchangeEngine::run(std::vector<double>& x, const Cell& cell, double temperature, const std::function<double()>& potential) {
  xw_ = &x;
  try {
    const bool any = run(static_cast<const std::vector<double>&>(x), cell, temperature, potential);
    xw_ = nullptr;
    return any;
  } catch (...) {
    xw_ = nullptr;
    throw;
  }
}

bool ExchangeEngine::run(const std::vector<double>& x, const Cell& cell, double temperature, const std::function<double()>& potential) {
  if (o_.relax_steps > 0 && xw_ != &x) throw std::invalid_argument("exchange: relaxed moves write positions: run them with a mutable x");
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
