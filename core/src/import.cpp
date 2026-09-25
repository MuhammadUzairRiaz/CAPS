#include "caps/import.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>

#include "caps/analysis.hpp"
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

// Bond orders from valences (all hydrogens explicit): 4 for aromatic bonds, else the Kekulé order.
void write_orders(System& s) {
  if (s.bonds.empty()) return;
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
