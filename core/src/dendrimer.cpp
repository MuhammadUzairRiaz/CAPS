// Dendrimers and hyperbranched molecules from AB_n units (dendrimer.hpp).
#include "caps/dendrimer.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>

#include "caps/edit.hpp"
#include "caps/elements.hpp"
#include "caps/molecule.hpp"
#include "caps/relax.hpp"
#include "caps/ring_edit.hpp"
#include "caps/rng.hpp"
#include "caps/uff.hpp"

namespace caps {

const std::vector<DendrimerPreset>& dendrimer_presets() {
  static const std::vector<DendrimerPreset> p = {
      {"pamam", "PAMAM (ethylenediamine core, amidoamine units, NH₂ ends)", "[*:2]N([*:2])CCN([*:2])[*:2]", "[*:1]CCC(=O)NCCN([*:2])[*:2]", "", 1},
      {"ppi", "PPI (1,4-diaminobutane core, propylene imine units, NH₂ ends)", "[*:2]N([*:2])CCCCN([*:2])[*:2]", "[*:1]CCCN([*:2])[*:2]", "", 0},
      {"bismpa", "bis-MPA polyester (trimethylolpropane core, OH ends)", "CCC(CO[*:2])(CO[*:2])CO[*:2]", "[*:1]C(=O)C(C)(CO[*:2])CO[*:2]", "[*:2]O", 0},
  };
  return p;
}

namespace {

uint64_t draw(std::mt19937_64& g, uint64_t n) { return random_u64(g) % n; }

// A unit, core or end as a 3D piece: its sites made hydrogens (they leave when joined); the A site and the B sites
struct Piece {
  System s;
  int a = -1;
  std::vector<int> b;
};
Piece make_piece(const std::string& smiles, const char* what) {
  MolGraph g = parse_smiles(smiles);
  Piece p;
  std::vector<int> a_sites, b_sites;
  for (size_t k = 0; k < g.atoms.size(); ++k) {
    auto& at = g.atoms[k];
    if (at.element != 0) continue;
    if (at.map == 1) a_sites.push_back(int(k));
    else if (at.map == 2) b_sites.push_back(int(k));
    else throw std::invalid_argument(std::string(what) + " '" + smiles + "': mark sites [*:1] (A) or [*:2] (B)");
    at.element = 1, at.hcount = 0, at.bracket = true;
  }
  for (int k : a_sites) {
    int deg = 0;
    for (const auto& bd : g.bonds) deg += bd.a == k || bd.b == k;
    if (deg != 1) throw std::invalid_argument(std::string(what) + ": a site has exactly one bond");
  }
  add_hydrogens(g);
  p.s = molecule_system(g, embed(g));
  for (auto& x : p.s.atoms) x.mol = 1;
  p.s.has_mol = true;
  p.a = a_sites.empty() ? -1 : a_sites[0];
  p.b = b_sites;
  if (a_sites.size() > 1) throw std::invalid_argument(std::string(what) + " '" + smiles + "' has more than one A site [*:1]");
  return p;
}

std::string hill(const System& s, double* mass) {
  std::map<std::string, int> n;
  *mass = 0;
  for (const auto& a : s.atoms) { ++n[element(a.element).symbol]; *mass += element(a.element).mass; }
  std::string f;
  for (const char* first : {"C", "H"})
    if (n.count(first)) { f += first + (n[first] > 1 ? std::to_string(n[first]) : ""); n.erase(first); }
  for (const auto& [el, c] : n) f += el + (c > 1 ? std::to_string(c) : "");
  return f;
}

}  // namespace

DendrimerTree perfect_tree(int core_sites, int n, int shells) {
  DendrimerTree t;
  std::vector<std::pair<int, int>> open;   // (unit, slot); unit −1 the core
  for (int k = 0; k < core_sites; ++k) open.push_back({-1, k});
  for (int sh = 0; sh < shells; ++sh) {
    std::vector<std::pair<int, int>> next;
    for (auto [u, sl] : open) {
      const int id = int(t.parent.size());
      t.parent.push_back(u), t.slot.push_back(sl);
      for (int k = 0; k < n; ++k) next.push_back({id, k});
    }
    open = std::move(next);
  }
  return t;
}

DendrimerTree sma_tree(int units, int n, int core_sites, uint64_t seed) {
  if (units < 1) throw std::invalid_argument("a hyperbranched molecule needs at least one unit");
  std::mt19937_64 g(seed);
  DendrimerTree t;
  std::vector<std::pair<int, int>> open;
  if (core_sites > 0) {
    for (int k = 0; k < core_sites; ++k) open.push_back({-1, k});
  } else {   // the focal unit
    t.parent.push_back(-1), t.slot.push_back(0);
    for (int k = 0; k < n; ++k) open.push_back({0, k});
  }
  while (int(t.parent.size()) < units) {
    const size_t pick = size_t(draw(g, open.size()));
    const auto [u, sl] = open[pick];
    open[pick] = open.back();
    open.pop_back();
    const int id = int(t.parent.size());
    t.parent.push_back(u), t.slot.push_back(sl);
    for (int k = 0; k < n; ++k) open.push_back({id, k});
  }
  return t;
}

DendrimerTree flory_tree(int units, int n, uint64_t seed) {
  if (units < 1) throw std::invalid_argument("a hyperbranched molecule needs at least one unit");
  std::mt19937_64 g(seed);
  const int N = units;
  // each B site reacted with probability p, the sequence redrawn until N − 1 sites carry a unit (the size condition);
  // p is free (the conditioned tree does not depend on it), this one makes the condition likeliest
  const double p = N > 1 ? double(N - 1) / (double(n) * N) : 0.0;
  std::vector<char> bits(size_t(N) * size_t(n));
  for (int tries = 0;; ++tries) {
    if (tries > 1000000) throw std::runtime_error("flory_tree: no sequence of the right size");
    long used = 0;
    for (auto& b : bits) { b = uniform01(g) < p ? 1 : 0; used += b; }
    if (used == N - 1) break;
  }
  // cycle lemma: the rotation that starts after the first minimum of the partial sums of (children − 1)
  long sum = 0, best = 1;
  int at = 0;
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < n; ++k) sum += bits[size_t(i) * n + k];
    sum -= 1;
    if (sum < best) best = sum, at = i;
  }
  std::vector<int> order(static_cast<size_t>(N));
  for (int i = 0; i < N; ++i) order[size_t(i)] = (at + 1 + i) % N;
  // preorder: each unit takes the next open slot of the last unit that still has one
  DendrimerTree t;
  std::vector<std::pair<int, int>> pending;
  for (int j = 0; j < N; ++j) {
    if (j == 0) t.parent.push_back(-1), t.slot.push_back(0);
    else {
      if (pending.empty()) throw std::logic_error("flory_tree: the rotation is not a tree");
      const auto [u, sl] = pending.back();
      pending.pop_back();
      t.parent.push_back(u), t.slot.push_back(sl);
    }
    const size_t src = size_t(order[size_t(j)]);
    for (int k = n - 1; k >= 0; --k)
      if (bits[src * n + k]) pending.push_back({j, k});
  }
  if (!pending.empty()) throw std::logic_error("flory_tree: open slots left");
  return t;
}

BranchCounts branch_counts(const DendrimerTree& t, int n) {
  BranchCounts c;
  const size_t N = t.parent.size();
  std::vector<int> used(N, 0), shell(N, 0);
  for (size_t k = 0; k < N; ++k) {
    if (t.parent[k] >= 0) {
      ++used[size_t(t.parent[k])];
      shell[k] = shell[size_t(t.parent[k])] + 1;
    } else {
      shell[k] = 1;
    }
  }
  c.by_used.assign(size_t(n) + 1, 0);
  for (size_t k = 0; k < N; ++k) {
    ++c.by_used[size_t(used[k])];
    if (size_t(shell[k]) > c.per_shell.size()) c.per_shell.resize(size_t(shell[k]), 0);
    ++c.per_shell[size_t(shell[k] - 1)];
  }
  if (n == 2) {
    c.terminal = c.by_used[0], c.linear = c.by_used[1], c.dendritic = c.by_used[2];
    const double D = c.dendritic, L = c.linear, T = c.terminal;
    // Frey's DB counts linear units as the defects; a tree of terminal units only (one shell on a core, PAMAM G0 or
    // PPI G1) has none, so it is perfectly branched like any other tree with L = 0
    c.db_frey = L == 0 ? (D + T > 0 ? 1.0 : NAN) : 2 * D / (2 * D + L);
    c.db_frechet = D + T + L > 0 ? (D + T) / (D + T + L) : NAN;
  }
  return c;
}

System build_dendrimer(const DendrimerSpec& spec0, DendrimerReport* report) {
  DendrimerSpec spec = spec0;
  DendrimerReport rep;
  if (!spec.preset.empty()) {
    const auto& ps = dendrimer_presets();
    const auto it = std::find_if(ps.begin(), ps.end(), [&](const DendrimerPreset& p) { return p.name == spec.preset; });
    if (it == ps.end()) throw std::invalid_argument("unknown dendrimer preset '" + spec.preset + "' (pamam, ppi, bismpa)");
    if (spec.core.empty()) spec.core = it->core;
    if (spec.unit.empty()) spec.unit = it->unit;
    if (spec.focal.empty()) spec.focal = it->focal;
    if (spec.generation >= 0) spec.shells = spec.generation + it->shell_offset;
    rep.notes.push_back(it->title);
  }
  if (spec.unit.empty()) throw std::invalid_argument("name the branching unit (a SMILES with one [*:1] and its [*:2] sites)");
  const Piece unit = make_piece(spec.unit, "the unit");
  if (unit.a < 0) throw std::invalid_argument("the unit needs its A site [*:1]");
  if (unit.b.empty()) throw std::invalid_argument("the unit needs at least one B site [*:2]");
  const int n = int(unit.b.size());
  Piece core;
  const bool has_core = !spec.core.empty() && spec.growth != "flory";
  if (!spec.core.empty() && spec.growth == "flory") rep.notes.push_back("Flory growth is the AB_n polycondensation without a core: the core is left out");
  if (has_core) {
    core = make_piece(spec.core, "the core");
    if (core.a >= 0) throw std::invalid_argument("the core has B sites [*:2] only");
    if (core.b.empty()) throw std::invalid_argument("the core needs B sites [*:2]");
  }
  Piece end;
  if (!spec.end.empty()) {
    end = make_piece(spec.end, "the end group");
    if (end.a < 0 || !end.b.empty()) throw std::invalid_argument("the end group has one A site [*:1] and no B sites");
  }
  DendrimerTree tree;
  if (spec.growth == "perfect") {
    if (!has_core) throw std::invalid_argument("a perfect dendrimer grows on a core");
    if (spec.shells < 1 || spec.shells > 10) throw std::invalid_argument("1 to 10 shells");
    tree = perfect_tree(int(core.b.size()), n, spec.shells);
  } else if (spec.growth == "sma") {
    tree = sma_tree(spec.units, n, has_core ? int(core.b.size()) : 0, spec.seed);
  } else if (spec.growth == "flory") {
    tree = flory_tree(spec.units, n, spec.seed);
  } else {
    throw std::invalid_argument("growth perfect, sma or flory");
  }
  const size_t N = tree.parent.size();
  // the size before building (each join loses two site atoms)
  size_t open_left = (has_core ? core.b.size() : 0) + N * size_t(n) - (has_core ? N : N - 1);
  const size_t estimate = (has_core ? core.s.atoms.size() : 0) + N * unit.s.atoms.size() - 2 * (has_core ? N : N - 1) +
                          (spec.end.empty() ? 0 : open_left * (end.s.atoms.size() - 2));
  if (estimate > 15000) throw std::invalid_argument("about " + std::to_string(estimate) + " atoms: CAPS builds molecules of up to 15000 atoms");

  // assembly: the tracked site atoms follow every deletion (each join removes one atom of the structure)
  System s;
  std::vector<std::vector<int>> sites(N);   // each unit's B-site atoms in s (−1 once used)
  std::vector<int> core_sites;
  auto shift = [&](int gone) {
    for (auto& v : sites) for (auto& x : v) if (x > gone) --x;
    for (auto& x : core_sites) if (x > gone) --x;
  };
  auto map_piece = [&](const RingEditResult& r, const Piece& p, int q) { return int(r.added[size_t(q < p.a ? q : q - 1)]); };
  size_t first = 0;
  if (has_core) {
    s = core.s;
    core_sites = core.b;
  } else {   // the focal unit as built, its A site kept as a hydrogen
    s = unit.s;
    sites[0] = unit.b;
    first = 1;
  }
  if (!has_core && !spec.focal.empty()) {   // the focal unit's A site capped
    const Piece f = make_piece(spec.focal, "the focal group");
    if (f.a >= 0 || f.b.size() != 1) throw std::invalid_argument("the focal group has one site [*:2]");
    const int la = unit.a;
    join_pieces(s, f.s, uint32_t(la), uint32_t(f.b[0]));
    shift(la);
    rep.notes.push_back("focal group " + spec.focal + " on the first unit");
  }
  for (size_t k = first; k < N; ++k) {
    int& site = tree.parent[k] < 0 ? core_sites[size_t(tree.slot[k])] : sites[size_t(tree.parent[k])][size_t(tree.slot[k])];
    const int la = site;
    site = -1;
    const auto r = join_pieces(s, unit.s, uint32_t(la), uint32_t(unit.a));
    shift(la);
    for (int q : unit.b) sites[k].push_back(map_piece(r, unit, q));
  }
  // end groups on every B site left
  int ends = 0;
  if (!spec.end.empty()) {
    std::vector<int*> open;
    for (auto& x : core_sites) if (x >= 0) open.push_back(&x);
    for (auto& v : sites) for (auto& x : v) if (x >= 0) open.push_back(&x);
    for (int* x : open) {
      const int la = *x;
      *x = -1;
      join_pieces(s, end.s, uint32_t(la), uint32_t(end.a));
      shift(la);
      ++ends;
    }
  }
  for (auto& a : s.atoms) a.mol = 1;
  s.has_mol = true;
  if (spec.relax && spec.relax_steps > 0) {
    RelaxOptions ro;
    ro.field = std::make_shared<ForceField>(assign_uff(s));
    ro.ftol = 0.5, ro.max_iterations = spec.relax_steps, ro.pushoff = true;
    relax(s, ro);
    rep.notes.push_back("relaxed with UFF (" + std::to_string(spec.relax_steps) + " steps at most)");
  }
  rep.units = int(N);
  rep.n = n;
  rep.core_sites = has_core ? int(core.b.size()) : 0;
  rep.ends = ends;
  rep.counts = branch_counts(tree, n);
  rep.shells = int(rep.counts.per_shell.size());
  if (n == 2) rep.db_limit = spec.growth == "sma" ? 2.0 / 3.0 : spec.growth == "flory" ? 0.5 : 1.0;
  rep.formula = hill(s, &rep.mass);
  if (report) *report = rep;
  return s;
}

}  // namespace caps
