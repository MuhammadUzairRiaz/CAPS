#include "caps/molinfo.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <set>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/molecule.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {

double mono_mass(int z) {
  switch (z) {
    case 1: return 1.00782503207; case 3: return 7.0160045; case 5: return 11.0093054; case 6: return 12.0;
    case 7: return 14.0030740048; case 8: return 15.99491461956; case 9: return 18.99840322; case 11: return 22.9897692809;
    case 12: return 23.985041700; case 13: return 26.98153863; case 14: return 27.9769265325; case 15: return 30.97376163;
    case 16: return 31.97207100; case 17: return 34.96885268; case 19: return 38.96370668; case 20: return 39.96259098;
    case 22: return 47.9479463; case 26: return 55.9349375; case 29: return 62.9295975; case 30: return 63.9291422;
    case 35: return 78.9183371; case 53: return 126.904473;
    default: return std::nan("");
  }
}

// eigenvalues of a symmetric 3 × 3 matrix, ascending (Jacobi)
std::array<double, 3> eig3(double a[3][3]) {
  double m[3][3];
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) m[i][j] = a[i][j];
  for (int sweep = 0; sweep < 50; ++sweep) {
    double off = std::fabs(m[0][1]) + std::fabs(m[0][2]) + std::fabs(m[1][2]);
    if (off < 1e-14) break;
    for (int p = 0; p < 2; ++p)
      for (int q = p + 1; q < 3; ++q) {
        if (std::fabs(m[p][q]) < 1e-300) continue;
        const double th = 0.5 * std::atan2(2 * m[p][q], m[q][q] - m[p][p]);
        const double c = std::cos(th), s = std::sin(th);
        for (int k = 0; k < 3; ++k) {
          const double mkp = m[k][p], mkq = m[k][q];
          m[k][p] = c * mkp - s * mkq;
          m[k][q] = s * mkp + c * mkq;
        }
        for (int k = 0; k < 3; ++k) {
          const double mpk = m[p][k], mqk = m[q][k];
          m[p][k] = c * mpk - s * mqk;
          m[q][k] = s * mpk + c * mqk;
        }
      }
  }
  std::array<double, 3> e{m[0][0], m[1][1], m[2][2]};
  std::sort(e.begin(), e.end());
  return e;
}

}  // namespace

MoleculeInfo molecule_info(const System& s0, uint32_t atom) {
  if (atom >= s0.atoms.size()) throw std::out_of_range("atom " + std::to_string(atom + 1) + " is not in the structure");
  System s = s0;
  if (s.cell.valid() && !s.unwrapped) make_molecules_whole(s);
  const auto mol = s.molecules();
  MoleculeInfo r;
  r.molecule = mol[atom];
  std::vector<uint32_t> idx;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) if (mol[i] == r.molecule) idx.push_back(i);
  std::vector<int> local(s.atoms.size(), -1);
  for (size_t k = 0; k < idx.size(); ++k) local[idx[k]] = int(k);
  // the molecule as its own system (for perception)
  System m;
  m.cell = s.cell;
  for (uint32_t i : idx) m.atoms.push_back(s.atoms[i]);
  for (const auto& b : s.bonds)
    if (local[b.i] >= 0 && local[b.j] >= 0) m.bonds.push_back({uint32_t(local[b.i]), uint32_t(local[b.j]), b.order});
  m.cell = Cell{};
  r.atoms = int(m.atoms.size());
  r.bonds = int(m.bonds.size());
  // formula (Hill), masses, DBE
  std::map<int, int> count;
  bool mono_ok = true;
  for (const auto& a : m.atoms) {
    count[a.element]++;
    r.mass += element(a.element).mass;
    const double mm = mono_mass(a.element);
    if (std::isnan(mm)) mono_ok = false; else r.monoisotopic += mm;
  }
  if (!mono_ok) r.monoisotopic = std::nan("");
  auto add = [&](int z) { if (count.count(z)) { r.formula += element(z).symbol; if (count[z] > 1) r.formula += std::to_string(count[z]); count.erase(z); } };
  if (count.count(6)) { add(6); add(1); }
  std::vector<std::pair<std::string, int>> rest;
  for (const auto& [z, n] : count) rest.push_back({element(z).symbol, n});
  std::sort(rest.begin(), rest.end());
  for (const auto& [sym, n] : rest) { r.formula += sym; if (n > 1) r.formula += std::to_string(n); }
  int c = 0, h = 0, n = 0, x = 0;
  for (const auto& a : m.atoms) {
    const int z = a.element;
    if (z == 6 || z == 14) ++c; else if (z == 1) ++h; else if (z == 7 || z == 15) ++n; else if (z == 9 || z == 17 || z == 35 || z == 53) ++x;
  }
  r.dbe = c - (h + x) / 2.0 + n / 2.0 + 1;
  // shape: inertia about the centre of mass, Rg
  double M = 0;
  Vec3 com{0, 0, 0};
  for (const auto& a : m.atoms) { const double w = element(a.element).mass; com = com + a.pos * w; M += w; }
  com = com * (1.0 / M);
  double I[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, rg2 = 0;
  for (const auto& a : m.atoms) {
    const double w = element(a.element).mass;
    const Vec3 d = a.pos - com;
    const double r2 = dot(d, d);
    rg2 += w * r2;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) I[i][j] += w * ((i == j ? r2 : 0) - d[i] * d[j]);
  }
  const auto e = eig3(I);
  for (int k = 0; k < 3; ++k) r.inertia[k] = e[size_t(k)];
  r.inertia_defect = e[2] - e[0] - e[1];
  r.rg = std::sqrt(rg2 / M);
  // charges and dipole (debye; about the centre of mass, which is exact for a neutral molecule)
  Vec3 mu{0, 0, 0};
  for (const auto& a : m.atoms) { r.net_charge += a.charge; mu = mu + (a.pos - com) * a.charge; r.has_charges = r.has_charges || std::fabs(a.charge) > 1e-9; }
  r.dipole = r.has_charges ? norm(mu) * 4.80320471 : std::nan("");
  // SMILES from the perceived bonds, rings, rotatable bonds
  const Perception p = perceive(m);
  r.rings = int(p.rings.size());
  MolGraph g;
  std::vector<int> gi(m.atoms.size(), -1);
  for (uint32_t i = 0; i < m.atoms.size(); ++i) {
    if (m.atoms[i].element == 1 && p.nb[i].size() == 1 && m.atoms[p.nb[i][0]].element != 1) continue;   // written as implicit H
    MolAtom a;
    a.element = m.atoms[i].element;
    a.aromatic = p.aromatic[i];
    a.charge = p.charge[i];
    const bool organic = a.element == 5 || a.element == 6 || a.element == 7 || a.element == 8 || a.element == 9 || a.element == 15 ||
                         a.element == 16 || a.element == 17 || a.element == 35 || a.element == 53;
    a.bracket = !organic || a.charge != 0;
    a.hcount = a.bracket ? p.hcount[i] : -1;
    gi[i] = int(g.atoms.size());
    g.atoms.push_back(a);
  }
  for (uint32_t i = 0; i < m.atoms.size(); ++i)
    for (size_t k = 0; k < p.nb[i].size(); ++k) {
      const uint32_t j = p.nb[i][k];
      if (j <= i || gi[i] < 0 || gi[j] < 0) continue;
      MolBond b;
      b.a = gi[i], b.b = gi[j];
      b.order = p.arom_bond[i][k] ? 4 : p.order[i][k];
      g.bonds.push_back(b);
      // rotatable: acyclic single bond, both ends with another heavy neighbour
      auto heavy_nb = [&](uint32_t a) { int n2 = 0; for (uint32_t w : p.nb[a]) n2 += m.atoms[w].element != 1; return n2; };
      if (!p.ring_bond[i][k] && p.order[i][k] == 1 && !p.arom_bond[i][k] && heavy_nb(i) > 1 && heavy_nb(j) > 1) {
        bool triple = false;
        for (uint32_t a : {i, j}) for (int o : p.order[a]) triple = triple || o == 3;
        if (!triple) ++r.rotatable;
      }
    }
  g.heavy = int(g.atoms.size());
  try { r.smiles = write_smiles(g); } catch (const std::exception&) { r.smiles = ""; }
  return r;
}

SasaResult sasa(const System& s, double probe, int points) {
  SasaResult r;
  r.probe = probe;
  r.points = std::max(8, points);
  const size_t n = s.atoms.size();
  r.area.assign(n, 0.0);
  if (!n) return r;
  // golden-spiral points on the unit sphere
  std::vector<Vec3> sp;
  const double golden = M_PI * (3 - std::sqrt(5.0));
  for (int k = 0; k < r.points; ++k) {
    const double y = 1 - 2 * (k + 0.5) / r.points, rad = std::sqrt(1 - y * y), th = golden * k;
    sp.push_back({std::cos(th) * rad, y, std::sin(th) * rad});
  }
  std::vector<double> R(n);
  double rmax = 0;
  for (size_t i = 0; i < n; ++i) { R[i] = element(s.atoms[i].element).vdw + probe; rmax = std::max(rmax, R[i]); }
  const bool pbc = s.cell.valid();
  // cell list on the atom positions (bins of 2 rmax)
  const double bin = 2 * rmax;
  Vec3 lo{1e300, 1e300, 1e300};
  for (const auto& a : s.atoms) for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], a.pos[k]);
  std::map<std::array<long, 3>, std::vector<uint32_t>> cells;
  auto key = [&](const Vec3& p) { return std::array<long, 3>{long(std::floor((p[0] - lo[0]) / bin)), long(std::floor((p[1] - lo[1]) / bin)), long(std::floor((p[2] - lo[2]) / bin))}; };
  std::vector<std::vector<uint32_t>> near(n);
  if (!pbc) {
    for (uint32_t i = 0; i < n; ++i) cells[key(s.atoms[i].pos)].push_back(i);
    for (uint32_t i = 0; i < n; ++i) {
      const auto k0 = key(s.atoms[i].pos);
      for (long a = -1; a <= 1; ++a) for (long b = -1; b <= 1; ++b) for (long c = -1; c <= 1; ++c) {
        auto it = cells.find({k0[0] + a, k0[1] + b, k0[2] + c});
        if (it == cells.end()) continue;
        for (uint32_t j : it->second) {
          if (j == i) continue;
          const Vec3 d = s.atoms[j].pos - s.atoms[i].pos;
          if (dot(d, d) < (R[i] + R[j]) * (R[i] + R[j])) near[i].push_back(j);
        }
      }
    }
  } else {
    // periodic: all pairs by minimum image within the cutoff (cells in fractional space would be faster; SASA is
    // interactive on the structures the page is meant for)
    for (uint32_t i = 0; i < n; ++i)
      for (uint32_t j = i + 1; j < n; ++j) {
        const Vec3 d = s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos);
        if (dot(d, d) < (R[i] + R[j]) * (R[i] + R[j])) { near[i].push_back(j); near[j].push_back(i); }
      }
  }
  for (uint32_t i = 0; i < n; ++i) {
    int free = 0;
    for (const auto& u : sp) {
      const Vec3 p = s.atoms[i].pos + u * R[i];
      bool buried = false;
      for (uint32_t j : near[i]) {
        const Vec3 d = pbc ? s.cell.minimum_image(p - s.atoms[j].pos) : p - s.atoms[j].pos;
        if (dot(d, d) < R[j] * R[j]) { buried = true; break; }
      }
      free += !buried;
    }
    r.area[i] = 4 * M_PI * R[i] * R[i] * free / r.points;
    r.total += r.area[i];
  }
  return r;
}

std::string repeat_unit_smiles(const System& s0, uint32_t head, uint32_t tail) {
  if (head >= s0.atoms.size() || tail >= s0.atoms.size()) throw std::out_of_range("pick two atoms of the structure");
  if (head == tail) throw std::invalid_argument("the head and the tail must be two different atoms");
  System s = s0;
  if (s.cell.valid() && !s.unwrapped) make_molecules_whole(s);
  const auto mol = s.molecules();
  if (mol[head] != mol[tail]) throw std::invalid_argument("the head and the tail are in different molecules");
  std::vector<uint32_t> idx;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) if (mol[i] == mol[head]) idx.push_back(i);
  std::vector<int> local(s.atoms.size(), -1);
  for (size_t k = 0; k < idx.size(); ++k) local[idx[k]] = int(k);
  System m;
  for (uint32_t i : idx) m.atoms.push_back(s.atoms[i]);
  for (const auto& b : s.bonds)
    if (local[b.i] >= 0 && local[b.j] >= 0) m.bonds.push_back({uint32_t(local[b.i]), uint32_t(local[b.j]), b.order});
  const Perception p = perceive(m);
  // the attachment points: a picked hydrogen itself, else a hydrogen of the picked atom (not the other's)
  auto point = [&](uint32_t at, int avoid) -> uint32_t {
    const uint32_t a = uint32_t(local[at]);
    if (m.atoms[a].element == 1) return a;
    for (uint32_t w : p.nb[a]) if (m.atoms[w].element == 1 && int(w) != avoid) return w;
    throw std::invalid_argument("atom " + std::to_string(at + 1) + " (" + element(m.atoms[a].element).symbol + ") has no hydrogen to become the attachment point");
  };
  const uint32_t hp = point(head, -1), tp = point(tail, int(hp));
  if (hp == tp) throw std::invalid_argument("the head and the tail give the same attachment point");
  // the graph: the head's * first, then the other atoms; the two points become *, other hydrogens implicit
  MolGraph g;
  std::vector<int> gi(m.atoms.size(), -1);
  std::vector<uint32_t> order{hp};
  for (uint32_t i = 0; i < m.atoms.size(); ++i) if (i != hp) order.push_back(i);
  for (uint32_t i : order) {
    const bool dummy = i == hp || i == tp;
    if (!dummy && m.atoms[i].element == 1 && p.nb[i].size() == 1 && m.atoms[p.nb[i][0]].element != 1) continue;   // implicit H
    MolAtom a;
    a.element = dummy ? 0 : m.atoms[i].element;
    a.aromatic = !dummy && p.aromatic[i];
    a.charge = dummy ? 0 : p.charge[i];
    const int z = a.element;
    const bool organic = z == 0 || z == 5 || z == 6 || z == 7 || z == 8 || z == 9 || z == 15 || z == 16 || z == 17 || z == 35 || z == 53;
    a.bracket = !organic || a.charge != 0;
    // a bracket atom states its hydrogens: the perceived count less the one that became *
    int h = dummy ? 0 : p.hcount[i];
    if (!dummy) for (uint32_t w : p.nb[i]) if (w == hp || w == tp) --h;
    a.hcount = a.bracket ? std::max(0, h) : -1;
    gi[i] = int(g.atoms.size());
    g.atoms.push_back(a);
  }
  for (uint32_t i = 0; i < m.atoms.size(); ++i)
    for (size_t k = 0; k < p.nb[i].size(); ++k) {
      const uint32_t j = p.nb[i][k];
      if (j <= i || gi[i] < 0 || gi[j] < 0) continue;
      MolBond b;
      b.a = gi[i], b.b = gi[j];
      b.order = p.arom_bond[i][k] ? 4 : p.order[i][k];
      g.bonds.push_back(b);
    }
  g.heavy = int(g.atoms.size());
  return write_smiles(g);
}

}  // namespace caps
