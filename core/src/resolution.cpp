// CAPS model resolution: see caps/resolution.hpp.
#include "caps/resolution.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {
std::vector<std::vector<uint32_t>> adjacency(const System& s) {
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  return nb;
}
// A type per label with its mass (the site's mass: a united CH2 weighs C + 2 H).
int type_of(System& out, std::map<std::string, int>& types, const std::string& label, double mass) {
  auto it = types.find(label);
  if (it != types.end()) return it->second;
  TypeInfo t;
  t.type = int(out.types.size()) + 1;
  t.label = label;
  t.mass = mass;
  out.types.push_back(t);
  return types[label] = t.type;
}
}  // namespace

ResolutionReport all_atom_summary(const System& s) {
  ResolutionReport r;
  r.sites = int(s.atoms.size());
  for (const auto& a : s.atoms) {
    r.mass += element(a.element).mass;
    r.hydrogens += a.element == 1 ? 1 : 0;
  }
  return r;
}

System united_atom(const System& s, ResolutionReport* rep) {
  const size_t n = s.atoms.size();
  const auto nb = adjacency(s);
  std::vector<int> host(n, -1);   // the carbon a hydrogen folds into
  for (size_t i = 0; i < n; ++i)
    if (s.atoms[i].element == 1 && nb[i].size() == 1 && s.atoms[nb[i][0]].element == 6) host[i] = int(nb[i][0]);
  System out;
  out.title = s.title.empty() ? "united-atom" : s.title + " (united-atom)";
  out.cell = s.cell;
  out.has_mol = s.has_mol;
  out.unwrapped = s.unwrapped;
  out.bonds_from_file = true;
  std::map<std::string, int> types;
  ResolutionReport r;
  r.site_of.assign(n, -1);
  std::vector<int> hcount(n, 0);
  std::vector<double> hcharge(n, 0.0);
  for (size_t i = 0; i < n; ++i)
    if (host[i] >= 0) ++hcount[size_t(host[i])], hcharge[size_t(host[i])] += s.atoms[i].charge;
  for (size_t i = 0; i < n; ++i) {
    if (host[i] >= 0) continue;
    Atom a = s.atoms[i];
    a.id = int64_t(out.atoms.size() + 1);
    const int h = hcount[i];
    const double mass = element(a.element).mass + h * element(1).mass;
    std::string label = element(a.element).symbol;
    if (a.element == 6 && h > 0) label = h == 1 ? "CH" : "CH" + std::to_string(h);
    a.name = label;
    a.charge += hcharge[i];
    a.type = type_of(out, types, label, mass);
    r.site_of[i] = int(out.atoms.size());
    out.atoms.push_back(a);
    r.mass += mass;
    r.hydrogens += a.element == 1 ? 1 : 0;
  }
  for (size_t i = 0; i < n; ++i)
    if (host[i] >= 0) r.site_of[i] = r.site_of[size_t(host[i])];
  for (const auto& b : s.bonds) {
    if (host[b.i] >= 0 || host[b.j] >= 0) continue;
    out.bonds.push_back({uint32_t(r.site_of[b.i]), uint32_t(r.site_of[b.j]), b.order});
  }
  r.sites = int(out.atoms.size());
  int folded = 0;
  for (size_t i = 0; i < n; ++i) folded += host[i] >= 0 ? 1 : 0;
  r.notes.push_back(std::to_string(folded) + " hydrogens on carbon folded into their carbons; " + std::to_string(r.hydrogens) + " polar hydrogens kept");
  if (rep) *rep = std::move(r);
  return out;
}

System coarse_grain(const System& s, int per_bead, ResolutionReport* rep) {
  if (per_bead < 1) throw std::invalid_argument("coarse_grain: at least one backbone atom per bead");
  const size_t n = s.atoms.size();
  const auto nb = adjacency(s);
  int nm = 0;
  const auto mol = s.molecules(&nm);
  // bead of each backbone atom, then every other atom takes the bead of the nearest backbone atom along bonds
  std::vector<int> bead(n, -1);
  int beads = 0;
  std::vector<int> mol_of_bead;
  std::vector<char> has(size_t(std::max(nm, 1)), 0);
  for (const auto& path : backbones(s, 2)) {
    for (size_t k = 0; k < path.size(); ++k) {
      if (k % size_t(per_bead) == 0) {
        // a short tail joins the bead before it (fewer than half a bead)
        if (k > 0 && path.size() - k < size_t(per_bead + 1) / 2) { bead[path[k]] = beads - 1; continue; }
        mol_of_bead.push_back(mol[path[k]]);
        ++beads;
      }
      bead[path[k]] = beads - 1;
    }
    has[size_t(mol[path[0]])] = 1;
  }
  std::vector<uint32_t> queue;
  for (uint32_t i = 0; i < n; ++i) if (bead[i] >= 0) queue.push_back(i);
  for (size_t h = 0; h < queue.size(); ++h)
    for (uint32_t w : nb[queue[h]])
      if (bead[w] < 0) bead[w] = bead[queue[h]], queue.push_back(w);
  std::vector<int> mol_bead(size_t(std::max(nm, 1)), -1);   // molecules without a backbone: one bead
  for (uint32_t i = 0; i < n; ++i) {
    if (bead[i] >= 0) continue;
    int& b = mol_bead[size_t(mol[i])];
    if (b < 0) b = beads++, mol_of_bead.push_back(mol[i]);
    bead[i] = b;
  }
  // centre of mass per bead, unwrapped about the bead's first atom (minimum image in a periodic cell)
  std::vector<double> mass(size_t(beads), 0.0), charge(size_t(beads), 0.0);
  std::vector<Vec3> com(size_t(beads), Vec3{0, 0, 0});
  std::vector<int> first(size_t(beads), -1);
  const bool periodic = s.cell.valid();
  for (uint32_t i = 0; i < n; ++i) {
    const auto b = size_t(bead[i]);
    if (first[b] < 0) first[b] = int(i);
    Vec3 p = s.atoms[i].pos;
    if (periodic) p = s.atoms[size_t(first[b])].pos + s.cell.minimum_image(p - s.atoms[size_t(first[b])].pos);
    const double m = element(s.atoms[i].element).mass;
    com[b] = com[b] + p * m;
    mass[b] += m;
    charge[b] += s.atoms[i].charge;
  }
  System out;
  out.title = s.title.empty() ? "coarse-grained" : s.title + " (coarse-grained)";
  out.cell = s.cell;
  out.has_mol = true;
  out.bonds_from_file = true;
  std::map<std::string, int> types;
  ResolutionReport r;
  for (int b = 0; b < beads; ++b) {
    Atom a;
    a.id = b + 1;
    a.element = 6;   // drawn as a carbon-sized sphere; the type carries the bead's mass
    a.mol = int64_t(mol_of_bead[size_t(b)] + 1);
    a.pos = com[size_t(b)] * (1.0 / std::max(1e-12, mass[size_t(b)]));
    a.charge = charge[size_t(b)];
    char label[32];
    std::snprintf(label, sizeof label, "B%.0f", std::round(mass[size_t(b)]));
    a.name = label;
    a.type = type_of(out, types, label, mass[size_t(b)]);
    out.atoms.push_back(a);
    r.mass += mass[size_t(b)];
  }
  std::map<std::pair<int, int>, bool> seen;
  for (const auto& bd : s.bonds) {
    int x = bead[bd.i], y = bead[bd.j];
    if (x == y) continue;
    if (x > y) std::swap(x, y);
    if (seen.emplace(std::pair{x, y}, true).second) out.bonds.push_back({uint32_t(x), uint32_t(y), 1});
  }
  r.sites = beads;
  r.site_of = bead;
  r.notes.push_back(std::to_string(per_bead) + " backbone atoms per bead, centre of mass; " + std::to_string(beads) + " beads");
  if (rep) *rep = std::move(r);
  return out;
}

namespace {

// Rotation (3×3, row major) that best maps points p onto q (both relative to their own centres): Horn's quaternion.
std::array<double, 9> best_rotation(const std::vector<Vec3>& p, const std::vector<Vec3>& q) {
  double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (size_t k = 0; k < p.size(); ++k)
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) S[a][b] += p[k][a] * q[k][b];
  double N[4][4] = {
      {S[0][0] + S[1][1] + S[2][2], S[1][2] - S[2][1], S[2][0] - S[0][2], S[0][1] - S[1][0]},
      {S[1][2] - S[2][1], S[0][0] - S[1][1] - S[2][2], S[0][1] + S[1][0], S[2][0] + S[0][2]},
      {S[2][0] - S[0][2], S[0][1] + S[1][0], -S[0][0] + S[1][1] - S[2][2], S[1][2] + S[2][1]},
      {S[0][1] - S[1][0], S[2][0] + S[0][2], S[1][2] + S[2][1], -S[0][0] - S[1][1] + S[2][2]}};
  // Jacobi: the eigenvector of the largest eigenvalue is the quaternion
  double V[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  for (int sweep = 0; sweep < 50; ++sweep) {
    double off = 0;
    for (int i = 0; i < 4; ++i) for (int j = i + 1; j < 4; ++j) off += N[i][j] * N[i][j];
    if (off < 1e-20) break;
    for (int i = 0; i < 4; ++i)
      for (int j = i + 1; j < 4; ++j) {
        if (std::fabs(N[i][j]) < 1e-15) continue;
        const double th = 0.5 * std::atan2(2 * N[i][j], N[j][j] - N[i][i]);
        const double c = std::cos(th), sn = std::sin(th);
        for (int k = 0; k < 4; ++k) {   // columns i, j
          const double a = N[k][i], b = N[k][j];
          N[k][i] = c * a - sn * b, N[k][j] = sn * a + c * b;
        }
        for (int k = 0; k < 4; ++k) {   // rows i, j
          const double a = N[i][k], b = N[j][k];
          N[i][k] = c * a - sn * b, N[j][k] = sn * a + c * b;
        }
        for (int k = 0; k < 4; ++k) {
          const double a = V[k][i], b = V[k][j];
          V[k][i] = c * a - sn * b, V[k][j] = sn * a + c * b;
        }
      }
  }
  int m = 0;
  for (int i = 1; i < 4; ++i) if (N[i][i] > N[m][m]) m = i;
  const double w = V[0][m], x = V[1][m], y = V[2][m], z = V[3][m];
  return {w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y),
          2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x),
          2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z};
}

Vec3 turn(const std::array<double, 9>& R, const Vec3& v) {
  return {R[0] * v[0] + R[1] * v[1] + R[2] * v[2], R[3] * v[0] + R[4] * v[1] + R[5] * v[2], R[6] * v[0] + R[7] * v[1] + R[8] * v[2]};
}

// the smallest rotation taking unit vector a onto unit vector b
std::array<double, 9> align(Vec3 a, Vec3 b) {
  a = a * (1 / norm(a)), b = b * (1 / norm(b));
  const Vec3 v = cross(a, b);
  const double c = dot(a, b), s = norm(v);
  if (s < 1e-12) return c > 0 ? std::array<double, 9>{1, 0, 0, 0, 1, 0, 0, 0, 1} : std::array<double, 9>{-1, 0, 0, 0, -1, 0, 0, 0, 1};
  const double k = (1 - c) / (s * s);
  return {c + k * v[0] * v[0], k * v[0] * v[1] - v[2], k * v[0] * v[2] + v[1],
          k * v[0] * v[1] + v[2], c + k * v[1] * v[1], k * v[1] * v[2] - v[0],
          k * v[0] * v[2] - v[1], k * v[1] * v[2] + v[0], c + k * v[2] * v[2]};
}

}  // namespace

System backmap(const System& all_atom, const System& beads_in, int per_bead, BackmapReport* rep) {
  System aa = all_atom;
  if (aa.cell.valid() && !aa.unwrapped) make_molecules_whole(aa);
  ResolutionReport rr;
  const System ref = coarse_grain(aa, per_bead, &rr);
  if (beads_in.atoms.size() != ref.atoms.size())
    throw std::invalid_argument("the bead file has " + std::to_string(beads_in.atoms.size()) + " beads; this structure coarse-grains to " +
                                std::to_string(ref.atoms.size()) + " at " + std::to_string(per_bead) + " backbone atoms per bead");
  // the moved beads as chains: the reference's bonds and molecules, made whole across the cell
  System beads = beads_in;
  beads.bonds = ref.bonds;
  for (size_t k = 0; k < beads.atoms.size(); ++k) beads.atoms[k].mol = ref.atoms[k].mol;
  beads.has_mol = true;
  if (beads.cell.valid() && !beads.unwrapped) make_molecules_whole(beads);
  const auto nb = ref.neighbours();
  const size_t nbead = ref.atoms.size();
  std::vector<std::array<double, 9>> R(nbead);
  double turn2 = 0;
  for (size_t k = 0; k < nbead; ++k) {
    const Vec3 o = ref.atoms[k].pos, n = beads.atoms[k].pos;
    if (nb[k].empty()) { R[k] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; continue; }
    if (nb[k].size() == 1) {
      R[k] = align(ref.atoms[nb[k][0]].pos - o, beads.atoms[nb[k][0]].pos - n);
    } else {
      std::vector<Vec3> p, q;
      for (uint32_t w : nb[k]) p.push_back(ref.atoms[w].pos - o), q.push_back(beads.atoms[w].pos - n);
      R[k] = best_rotation(p, q);
    }
    const double tr = std::clamp((R[k][0] + R[k][4] + R[k][8] - 1) / 2, -1.0, 1.0);
    turn2 += std::pow(std::acos(tr) * 180 / 3.14159265358979323846, 2);
  }
  System out = aa;
  for (size_t i = 0; i < out.atoms.size(); ++i) {
    const int k = rr.site_of[i];
    out.atoms[i].pos = turn(R[size_t(k)], aa.atoms[i].pos - ref.atoms[size_t(k)].pos) + beads.atoms[size_t(k)].pos;
  }
  out.cell = beads.cell.valid() ? beads.cell : aa.cell;
  out.unwrapped = true;
  out.velocities.clear();
  double worst = 0;
  for (const auto& b : out.bonds)
    if (rr.site_of[b.i] != rr.site_of[b.j]) worst = std::max(worst, norm(out.atoms[b.j].pos - out.atoms[b.i].pos));
  if (rep) {
    rep->beads = int(nbead);
    rep->atoms = int(out.atoms.size());
    rep->rms_turn = std::sqrt(turn2 / double(std::max<size_t>(1, nbead)));
    rep->worst_bond = worst;
    char b[240];
    std::snprintf(b, sizeof b, "%d atoms carried by %d beads (%d backbone atoms each) · rotations %.1f° rms · longest bond between beads %.2f Å before relaxing",
                  rep->atoms, rep->beads, per_bead, rep->rms_turn, worst);
    rep->notes.push_back(b);
  }
  return out;
}

}  // namespace caps
