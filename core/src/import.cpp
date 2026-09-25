#include "caps/import.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {

std::string cell_text(const Cell& c) {
  if (!c.valid()) return "";
  const auto len = [](const Vec3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); };
  const double a = len(c.a), b = len(c.b), cc = len(c.c);
  const bool ortho = std::fabs(c.a[1]) + std::fabs(c.a[2]) + std::fabs(c.b[0]) + std::fabs(c.b[2]) + std::fabs(c.c[0]) + std::fabs(c.c[1]) < 1e-6;
  char buf[160];
  std::snprintf(buf, sizeof buf, "%.4g × %.4g × %.4g Å, %s", a, b, cc, ortho ? "orthogonal" : "triclinic");
  return buf;
}

}  // namespace

// Bond orders of a structure without hydrogens (a PDB of heavy atoms), from its geometry: hybridisation from bond
// angles (or the bond length for terminal atoms), planar six-membered rings of sp² atoms aromatic, then double and
// triple bonds paired between the remaining sp² and sp atoms, shortest bonds first.
void orders_from_geometry(System& s) {
  const size_t n = s.atoms.size();
  std::vector<std::vector<std::pair<uint32_t, size_t>>> nb(n);   // (neighbour, bond index)
  for (size_t k = 0; k < s.bonds.size(); ++k) nb[s.bonds[k].i].push_back({s.bonds[k].j, k}), nb[s.bonds[k].j].push_back({s.bonds[k].i, k});
  auto vec = [&](uint32_t a, uint32_t b) { Vec3 d = s.atoms[b].pos - s.atoms[a].pos; return s.cell.valid() ? s.cell.minimum_image(d) : d; };
  auto len = [&](uint32_t a, uint32_t b) { return norm(vec(a, b)); };
  std::vector<int> hyb(n, 3);   // 1 sp, 2 sp², 3 sp³
  for (uint32_t i = 0; i < n; ++i) {
    const auto& N = nb[i];
    const int z = s.atoms[i].element;
    if (N.size() >= 4 || z == 1) continue;
    if (N.size() >= 2) {
      double sum = 0;
      int k = 0;
      for (size_t a = 0; a < N.size(); ++a)
        for (size_t b = a + 1; b < N.size(); ++b) {
          const Vec3 u = vec(i, N[a].first), v = vec(i, N[b].first);
          sum += std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / 3.14159265358979323846;
          ++k;
        }
      const double mean = sum / k;
      hyb[i] = mean > 155 ? 1 : mean > 115 ? 2 : 3;
      if (N.size() == 3 && (z == 7)) hyb[i] = mean > 117 ? 2 : 3;   // planar N (amides, aromatic N)
    } else if (N.size() == 1) {   // terminal: the bond length against the single-bond sum of covalent radii
      const uint32_t j = N[0].first;
      const double ratio = len(i, j) / (element(z).covalent + element(s.atoms[j].element).covalent);
      hyb[i] = ratio < 0.83 ? 1 : ratio < 0.93 ? 2 : 3;
    }
  }
  for (auto& b : s.bonds) b.order = 1;
  // flat rings of five or six atoms through atoms with at most three bonds (shortest path back to the bond's start)
  auto flat_rings = [&](size_t size, auto&& member_ok) {
    std::vector<std::vector<uint32_t>> rings;
    for (size_t k = 0; k < s.bonds.size(); ++k) {
      const uint32_t u = s.bonds[k].i, v = s.bonds[k].j;
      if (!member_ok(u) || !member_ok(v)) continue;
      std::vector<int> prev(n, -1), dist(n, -1);
      std::deque<uint32_t> q{v};
      dist[v] = 0;
      while (!q.empty()) {
        const uint32_t a = q.front();
        q.pop_front();
        if (a == u || dist[a] >= int(size) - 1) continue;
        for (auto [w, bk] : nb[a]) {
          if (bk == k || !member_ok(w) || dist[w] >= 0) continue;
          dist[w] = dist[a] + 1, prev[w] = int(a);
          q.push_back(w);
        }
      }
      if (dist[u] != int(size) - 1) continue;
      std::vector<uint32_t> ring{u};
      for (int a = prev[u]; a >= 0; a = prev[size_t(a)]) ring.push_back(uint32_t(a));
      std::vector<Vec3> p;
      for (uint32_t a : ring) p.push_back(s.atoms[ring[0]].pos + vec(ring[0], a));
      Vec3 c{0, 0, 0};
      for (const auto& x : p) c = c + x * (1.0 / double(p.size()));
      Vec3 nrm = cross(p[1] - p[0], p[size - 2] - p[0]);
      if (norm(nrm) < 1e-9) continue;
      nrm = nrm * (1.0 / norm(nrm));
      bool flat = true;
      for (const auto& x : p) flat = flat && std::fabs(dot(x - c, nrm)) < 0.15;
      if (flat) rings.push_back(ring);
    }
    return rings;
  };
  // five-rings read ~108° angles as sp³: a flat one is conjugated (imidazole, pyrrole, furan, thiophene)
  for (const auto& ring : flat_rings(5, [&](uint32_t a) { return nb[a].size() <= 3 && s.atoms[a].element != 1; }))
    for (uint32_t a : ring) hyb[a] = 2;
  // aromatic six-rings: every atom sp², flat, and none carrying an exocyclic double bond to a terminal atom (C=O)
  auto exo_double = [&](uint32_t a) {
    for (auto [w, bk] : nb[a])
      if (nb[w].size() == 1 && hyb[w] == 2) return true;
    return false;
  };
  std::vector<char> in_arom(n, 0);
  for (const auto& ring : flat_rings(6, [&](uint32_t a) { return hyb[a] == 2; })) {
    if (std::any_of(ring.begin(), ring.end(), exo_double)) continue;
    for (size_t a = 0; a < ring.size(); ++a) {
      const uint32_t x = ring[a], y = ring[(a + 1) % ring.size()];
      for (auto [w, bk] : nb[x]) if (w == y) s.bonds[bk].order = 4;
      in_arom[x] = 1;
    }
  }
  // double and triple bonds: shortest candidate bonds first, each atom taking one
  std::vector<size_t> cand;
  for (size_t k = 0; k < s.bonds.size(); ++k) {
    const auto& b = s.bonds[k];
    if (b.order == 4 || in_arom[b.i] || in_arom[b.j]) continue;
    if ((hyb[b.i] <= 2 && hyb[b.j] <= 2)) cand.push_back(k);
  }
  std::sort(cand.begin(), cand.end(), [&](size_t a, size_t b) {
    return len(s.bonds[a].i, s.bonds[a].j) / (element(s.atoms[s.bonds[a].i].element).covalent + element(s.atoms[s.bonds[a].j].element).covalent) <
           len(s.bonds[b].i, s.bonds[b].j) / (element(s.atoms[s.bonds[b].i].element).covalent + element(s.atoms[s.bonds[b].j].element).covalent);
  });
  std::vector<int> used(n, 0);   // π bonds each atom still wants: sp 2, sp² 1
  for (uint32_t i = 0; i < n; ++i) used[i] = hyb[i] == 1 ? 2 : hyb[i] == 2 ? 1 : 0;
  for (uint32_t i = 0; i < n; ++i)   // a planar N with three bonds keeps its lone pair (amides, pyrrole-type N)
    if (s.atoms[i].element == 7 && nb[i].size() == 3) used[i] = 0;
  for (size_t k : cand) {
    auto& b = s.bonds[k];
    if (used[b.i] <= 0 || used[b.j] <= 0) continue;
    const int pi = hyb[b.i] == 1 && hyb[b.j] == 1 ? 2 : 1;
    b.order = 1 + pi;
    used[b.i] -= pi, used[b.j] -= pi;
  }
}

namespace {

// Bond orders from valences (all hydrogens explicit): 4 for aromatic bonds, else the Kekulé order. Without any hydrogen
// (heavy atoms only) the orders come from the geometry instead.
void write_orders(System& s) {
  if (s.bonds.empty()) return;
  if (std::none_of(s.atoms.begin(), s.atoms.end(), [](const Atom& a) { return a.element == 1; }) && s.atoms.size() > 2) {
    orders_from_geometry(s);
    s.notes.push_back("bond orders from the geometry (no hydrogens in the file): angles, flat sp² rings aromatic");
    return;
  }
  const Perception p = perceive(s);
  for (auto& b : s.bonds) {
    const auto& nb = p.nb[b.i];
    for (size_t k = 0; k < nb.size(); ++k)
      if (nb[k] == b.j) { b.order = p.arom_bond[b.i][k] ? 4 : p.order[b.i][k]; break; }
  }
}

}  // namespace

void apply_import(Trajectory& t, const ImportOptions& o, size_t bonds_in_file) {
  System& top = t.topology;
  if (!o.use_cell) {
    top.cell = Cell{};
    for (auto& c : t.cells) c = Cell{};
  }
  if (o.bonds == ImportOptions::None) {
    top.bonds.clear();
    top.bonds_from_file = false;
  } else if (o.bonds == ImportOptions::FromFile) {
    if (bonds_in_file == 0) {
      top.bonds.clear();
      top.notes.push_back("the file gives no bonds; none were perceived (Bonds: read from file)");
    }
  } else {
    BondOptions bo;
    bo.tolerance = std::clamp(o.tolerance, 0.0, 1.5);
    System f0 = t.frame(0);
    if (!o.use_cell) f0.cell = Cell{};
    top.bonds = perceive_bonds(f0, bo);
    top.bonds_from_file = false;
    char buf[120];
    std::snprintf(buf, sizeof buf, "%zu bonds perceived from distances (covalent radii + %.2f Å)", top.bonds.size(), bo.tolerance);
    top.notes.erase(std::remove_if(top.notes.begin(), top.notes.end(), [](const std::string& n) { return n.find("bonds perceived from distances") != std::string::npos; }), top.notes.end());
    top.notes.push_back(buf);
  }
  if (o.bond_orders && !top.bonds.empty()) {
    const bool have = std::any_of(top.bonds.begin(), top.bonds.end(), [](const Bond& b) { return b.order > 0; });
    if (!have) {
      try { write_orders(top); } catch (const std::exception& e) { top.notes.push_back(std::string("bond orders not assigned: ") + e.what()); }
    }
  }
  if (o.split) {
    System probe = top;
    probe.has_mol = false;
    for (auto& a : probe.atoms) a.mol = 0;
    int n = 0;
    const auto mol = probe.molecules(&n);
    for (size_t i = 0; i < top.atoms.size(); ++i) top.atoms[i].mol = mol[i] + 1;
    top.has_mol = true;
  }
  if (o.unwrap && !top.bonds.empty()) {
    for (size_t k = 0; k < t.frames(); ++k) {
      System f = t.frame(k);
      if (!f.cell.valid()) continue;
      make_molecules_whole(f);
      for (size_t i = 0; i < f.atoms.size(); ++i) t.positions[k][i] = f.atoms[i].pos;
    }
    for (size_t i = 0; i < top.atoms.size(); ++i) top.atoms[i].pos = t.positions[0][i];
    if (top.cell.valid()) top.unwrapped = true;
  }
}

Trajectory import_file(const std::string& path, const std::string& topology_path, const ImportOptions& o) {
  Trajectory t = open_file(path, topology_path);
  const size_t in_file = t.topology.bonds_from_file ? t.topology.bonds.size() : 0;
  apply_import(t, o, in_file);
  return t;
}

ImportPreview import_preview(const std::string& path, const ImportOptions& o, int heavy) {
  ImportPreview r;
  r.file = inspect_file(path, "", 8);
  OpenProgress one;
  one.max_frames = 1;
  Trajectory t = open_file(path, "", one);
  r.bonds_in_file = t.topology.bonds_from_file ? t.topology.bonds.size() : 0;
  apply_import(t, o, r.bonds_in_file);
  const System& s = t.topology;
  r.atoms = s.atoms.size();
  r.bonds = s.bonds.size();
  int nmol = 0;
  s.molecules(&nmol);
  r.molecules = size_t(nmol);
  for (const auto& b : s.bonds) {
    if (b.order == 1) ++r.single;
    else if (b.order == 2) ++r.dbl;
    else if (b.order == 3) ++r.triple;
    else if (b.order == 4) ++r.aromatic;
  }
  r.cell = cell_text(s.cell);
  r.notes = s.notes;

  // fragment: breadth-first over heavy atoms from the first carbon (else the first heavy atom), then their hydrogens
  const auto nb = s.neighbours();
  int start = -1;
  for (size_t i = 0; i < s.atoms.size() && start < 0; ++i) if (s.atoms[i].element == 6) start = int(i);
  for (size_t i = 0; i < s.atoms.size() && start < 0; ++i) if (s.atoms[i].element > 1) start = int(i);
  if (start < 0 && !s.atoms.empty()) start = 0;
  std::vector<int> map(s.atoms.size(), -1);
  std::vector<uint32_t> pick;
  if (start >= 0) {
    std::deque<uint32_t> q{uint32_t(start)};
    map[start] = 0;
    pick.push_back(start);
    while (!q.empty() && int(pick.size()) < heavy) {
      const uint32_t i = q.front();
      q.pop_front();
      for (uint32_t j : nb[i]) {
        if (map[j] >= 0 || s.atoms[j].element == 1 || int(pick.size()) >= heavy) continue;
        map[j] = int(pick.size());
        pick.push_back(j);
        q.push_back(j);
      }
    }
    r.fragment_heavy = int(pick.size());
    const size_t nheavy = pick.size();
    for (size_t k = 0; k < nheavy; ++k)
      for (uint32_t j : nb[pick[k]])
        if (map[j] < 0 && s.atoms[j].element == 1) { map[j] = int(pick.size()); pick.push_back(j); }
  }
  System& f = r.fragment;
  f.title = "import preview";
  f.types = s.types;
  // positions made whole about the first atom (minimum image), then centred
  Vec3 c{0, 0, 0};
  for (size_t k = 0; k < pick.size(); ++k) {
    Atom a = s.atoms[pick[k]];
    if (k > 0 && s.cell.valid()) {
      const Vec3 d = s.cell.minimum_image(a.pos - s.atoms[pick[0]].pos);
      a.pos = s.atoms[pick[0]].pos + d;
    }
    a.id = int64_t(k + 1);
    a.mol = 1;
    f.atoms.push_back(a);
    c = c + a.pos;
  }
  if (!f.atoms.empty()) {
    c = c * (1.0 / double(f.atoms.size()));
    for (auto& a : f.atoms) a.pos = a.pos - c;
  }
  for (const auto& b : s.bonds)
    if (map[b.i] >= 0 && map[b.j] >= 0) f.bonds.push_back({uint32_t(map[b.i]), uint32_t(map[b.j]), b.order});
  f.bonds_from_file = true;
  f.has_mol = true;
  return r;
}

}  // namespace caps
