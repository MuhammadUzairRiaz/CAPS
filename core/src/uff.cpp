// CAPS UFF: typing and parameters of the Universal Force Field for every element (see caps/uff.hpp).
#include "caps/uff.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <set>

#include "caps/elements.hpp"
#include "caps/qeq.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {

struct UffParam {
  const char* label;
  double r1, theta0, x1, D1, zeta, Z1, V1, U1, Xi, hard, radius;
};

const UffParam kUff[] = {
#include "uff_params.inc"
};
constexpr int kUffCount = int(sizeof(kUff) / sizeof(kUff[0]));

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kG = 332.06;         // UFF bond / angle force-constant prefactor, kcal Å / mol
constexpr double kLambda = 0.1332;    // bond-order correction

const UffParam* by_label(const std::string& l) {
  for (const auto& p : kUff)
    if (l == p.label) return &p;
  return nullptr;
}

std::string label_symbol(const char* l) {
  std::string s(l, std::min<size_t>(2, std::strlen(l)));
  if (s.size() == 2 && s[1] == '_') s.resize(1);
  return s;
}

std::string label_element(const char* l) {   // leading letters: "Be3+2" → Be, "C_3" → C
  std::string s(1, l[0]);
  if (l[1] >= 'a' && l[1] <= 'z') s += l[1];
  return s;
}
char label_hyb(const char* l) { return std::strlen(l) > 2 ? l[2] : 0; }
int label_ox(const char* l) {
  const char* p = std::strchr(l, '+');
  return p && std::isdigit(static_cast<unsigned char>(p[1])) ? p[1] - '0' : -1;
}

// Group (1–18) of a main-group element; 0 for d- and f-block.
int main_group(int z) {
  if (z == 1) return 1;
  if (z == 2) return 18;
  if (z >= 5 && z <= 10) return z + 8;
  if (z >= 13 && z <= 18) return z;
  if (z >= 31 && z <= 36) return z - 18;
  if (z >= 49 && z <= 54) return z - 36;
  if (z >= 81 && z <= 86) return z - 68;
  return 0;
}
bool alkali_like(int z) {   // groups 1 and 2 below H
  return z == 3 || z == 4 || z == 11 || z == 12 || z == 19 || z == 20 || z == 37 || z == 38 || z == 55 || z == 56 || z == 87 || z == 88;
}
bool group16(int z) { return z == 8 || z == 16 || z == 34 || z == 52 || z == 84; }

// Hybridisation codes: 0 none (s), 1 sp, 2 sp2, 3 sp3, 4 square planar, 5 trigonal bipyramidal, 6 octahedral.
struct AtomState {
  int z = 0, deg = 0, hyb = 0, valence = 0, charge = 0, lone_pairs = 0;
  bool pi = false;          // has a double, triple or aromatic bond
  bool promoted = false;    // lone-pair N / O conjugated with a neighbouring π system (amide N, ester O, aniline N)
  bool resonant = false;    // UFF "R": aromatic or conjugated sp2 C, N, O, S
  std::string label, why;
};

double bond_order_value(const Perception& p, uint32_t a, size_t k) { return p.arom_bond[a][k] ? 1.5 : double(p.order[a][k]); }

// Is the angle between two neighbours of c (from positions) near 180°? Used for square planar and axial pairs.
double cos_angle(const System& s, uint32_t c, uint32_t a, uint32_t b) {
  const Vec3 u = s.atoms[a].pos - s.atoms[c].pos, v = s.atoms[b].pos - s.atoms[c].pos;
  const double nu = norm(u), nv = norm(v);
  if (nu < 1e-6 || nv < 1e-6) return 2.0;   // no geometry
  return dot(u, v) / (nu * nv);
}

std::vector<AtomState> type_atoms(const System& s, const Perception& p) {
  const size_t n = s.atoms.size();
  std::vector<AtomState> A(n);
  for (size_t i = 0; i < n; ++i) {
    auto& a = A[i];
    a.z = s.atoms[i].element;
    a.deg = int(p.nb[i].size());
    a.charge = p.charge[i];
    for (size_t k = 0; k < p.nb[i].size(); ++k) {
      a.valence += p.order[i][k];
      if (p.order[i][k] >= 2 || p.arom_bond[i][k]) a.pi = true;
    }
    if (p.aromatic[i]) a.pi = true;
  }
  // hybridisation of the main-group atoms from the steric number; metals from their coordination
  for (size_t i = 0; i < n; ++i) {
    auto& a = A[i];
    const int g = main_group(a.z);
    if (a.z == 1) { a.hyb = 0; continue; }
    if (g >= 13) {
      a.lone_pairs = std::max(0, (g - 10) - a.charge - a.valence) / 2;
      const int sn = a.deg + a.lone_pairs;
      if (p.aromatic[i]) a.hyb = 2;
      else if (a.deg == 0) a.hyb = 0;
      else if (sn <= 2) a.hyb = 1;
      else if (sn == 3) a.hyb = 2;
      else if (sn == 4) a.hyb = 3;
      else if (sn == 5) a.hyb = 5;
      else a.hyb = 6;
      continue;
    }
    // d- and f-block: the coordination geometry
    if (a.deg <= 0) a.hyb = 0;
    else if (a.deg <= 2) a.hyb = 1;
    else if (a.deg == 3) a.hyb = 2;
    else if (a.deg == 4) {
      bool trans = false, geometry = false;
      for (size_t x = 0; x < p.nb[i].size(); ++x)
        for (size_t y = x + 1; y < p.nb[i].size(); ++y) {
          const double c = cos_angle(s, uint32_t(i), p.nb[i][x], p.nb[i][y]);
          if (c <= 1.0) geometry = true;
          if (c < -0.85) trans = true;
        }
      if (geometry) a.hyb = trans ? 4 : 3;
      else {   // no coordinates: the element's own table geometry (Pd, Pt, Ni, Au are square planar)
        a.hyb = 3;
        for (const auto& q : kUff)
          if (element_from_symbol(label_symbol(q.label)) == a.z && label_hyb(q.label) == '4') a.hyb = 4;
      }
    } else if (a.deg == 5) a.hyb = 5;
    else a.hyb = 6;
  }
  // lone-pair N and O next to a π system are sp2 (conjugated), as in amides, esters, anilines and phenols
  for (size_t i = 0; i < n; ++i) {
    auto& a = A[i];
    if ((a.z != 7 && a.z != 8) || a.hyb != 3 || a.lone_pairs < 1 || a.pi) continue;
    for (uint32_t j : p.nb[i])
      if (A[j].pi && (A[j].hyb == 2 || A[j].hyb == 1)) { a.hyb = 2; a.promoted = true; break; }
  }
  auto conj = [&](uint32_t x) { return A[x].pi || A[x].promoted; };
  for (size_t i = 0; i < n; ++i) {
    auto& a = A[i];
    if (a.hyb != 2 || !(a.z == 6 || a.z == 7 || a.z == 8 || a.z == 16)) continue;
    if (p.aromatic[i] || a.promoted) { a.resonant = true; continue; }
    if (!a.pi) continue;
    for (size_t k = 0; k < p.nb[i].size() && !a.resonant; ++k) {
      const uint32_t b = p.nb[i][k];
      if (p.order[i][k] == 1 && !p.arom_bond[i][k]) {
        if (conj(b)) a.resonant = true;
      } else {   // the π partner: conjugated when it has another single bond to a π or lone-pair system
        for (size_t m = 0; m < p.nb[b].size(); ++m) {
          const uint32_t c = p.nb[b][m];
          if (c != i && p.order[b][m] == 1 && !p.arom_bond[b][m] && conj(c)) { a.resonant = true; break; }
        }
      }
    }
  }
  // labels: the exact UFF key when the table has it, otherwise the closest label of the element
  static const char* hyb_name[] = {"", "sp", "sp2", "sp3", "square planar", "trigonal bipyramidal", "octahedral"};
  for (size_t i = 0; i < n; ++i) {
    auto& a = A[i];
    if (a.z <= 0 || a.z > 103) continue;
    if (a.z == 1) {
      bool bridge = a.deg == 2;
      for (uint32_t j : p.nb[i]) bridge = bridge && A[j].z == 5;
      a.label = bridge ? "H_b" : "H_";
      a.why = bridge ? "H bridging two B" : "H";
      continue;
    }
    const char want = a.resonant ? 'R' : a.hyb == 0 ? 0 : char('0' + a.hyb);
    const bool metal = main_group(a.z) == 0 || alkali_like(a.z);
    const int ox = metal && a.charge > 0 ? a.charge : a.valence + std::max(0, a.charge);
    auto hyb_num = [](char h) { return h == 'R' ? 2.0 : h == 0 ? 0.0 : double(h - '0'); };
    double best = 1e9;
    const UffParam* pick = nullptr;
    for (const auto& q : kUff) {
      if (element_from_symbol(label_symbol(q.label)) != a.z) continue;
      if (!std::strcmp(q.label, "O_3_z") || !std::strcmp(q.label, "P_3+q")) continue;   // zeolite O, metal-bound P
      const char h = label_hyb(q.label);
      double cost = 0;
      if (want != h) {
        if (h == 0) cost = 0;   // the element's generic label (halogens, alkali metals)
        else if (want == 0) cost = 5;
        else if ((want == 'R' && h == '2') || (want == '2' && h == 'R')) cost = 5;
        else cost = 10 * std::fabs(hyb_num(want) - hyb_num(h));
      }
      const int lo = label_ox(q.label);
      if (lo >= 0) cost += (metal ? 0.1 : 1.0) * std::abs(lo - ox);   // a metal's bond count is not its oxidation state
      if (cost < best) { best = cost; pick = &q; }
    }
    if (!pick) continue;
    a.label = pick->label;
    const std::string sym = element(a.z).symbol;
    a.why = sym + (a.hyb ? std::string(" ") + hyb_name[a.hyb] : std::string()) + (a.resonant ? " (aromatic or conjugated)" : "") +
            (a.promoted ? " (lone pair conjugated)" : "") + ", " + std::to_string(a.deg) + " neighbours";
    if (best >= 1) a.why += "; closest UFF label";
  }
  return A;
}

double rest_length(double bo, const UffParam& a, const UffParam& b) {
  const double ri = a.r1, rj = b.r1;
  const double rbo = -kLambda * (ri + rj) * std::log(bo);
  const double d = std::sqrt(a.Xi) - std::sqrt(b.Xi);
  const double ren = ri * rj * d * d / (a.Xi * ri + b.Xi * rj);
  return ri + rj + rbo - ren;
}

double angle_k(double theta0, double bo12, double bo23, const UffParam& a, const UffParam& b, const UffParam& c) {
  const double ct = std::cos(theta0);
  const double r12 = rest_length(bo12, a, b), r23 = rest_length(bo23, b, c);
  const double r13 = std::sqrt(r12 * r12 + r23 * r23 - 2 * r12 * r23 * ct);
  const double inner = 3 * r12 * r23 * (1 - ct * ct) - r13 * r13 * ct;
  return 2 * kG * a.Z1 * c.Z1 / std::pow(r13, 5) * inner;
}

double equation17(double bo, const UffParam& a, const UffParam& b) { return 5 * std::sqrt(a.U1 * b.U1) * (1 + 4.18 * std::log(bo)); }

}  // namespace

ForceField default_forcefield(const System& s) {
  const bool ch = std::all_of(s.atoms.begin(), s.atoms.end(), [](const Atom& a) { return a.element == 1 || a.element == 6; });
  if (ch) {
    try {
      return assign_gaff(s);
    } catch (const FieldError&) {   // the built-in GAFF knows sp3 and aromatic carbon only (not C=C of rubbers)
    }
  }
  return assign_uff(s);
}

bool is_uff(const std::string& name) {
  std::string f = name.substr(name.find_last_of("/\\") == std::string::npos ? 0 : name.find_last_of("/\\") + 1);
  for (auto& ch : f) ch = char(std::tolower(static_cast<unsigned char>(ch)));
  return f == "uff" || f == "uff.json";
}

FFDef uff_definition() {
  FFDef d;
  d.name = "UFF";
  d.version = "1992";
  d.source = "Rappé, Casewit, Colwell, Goddard & Skiff, J. Am. Chem. Soc. 114, 10024 (1992); parameter table from RDKit (BSD)";
  d.references = {"A. K. Rappé, C. J. Casewit, K. S. Colwell, W. A. Goddard III, W. M. Skiff, J. Am. Chem. Soc. 114, 10024 (1992)"};
  for (const auto& q : kUff) {
    if (!std::strcmp(q.label, "O_3_z") || !std::strcmp(q.label, "P_3+q")) continue;
    FFType t;
    t.name = q.label;
    t.element = element_from_symbol(label_symbol(q.label));
    t.mass = element(t.element).mass;
    char b[160];
    std::snprintf(b, sizeof b, "r1 %.3f Å · θ0 %.2f° · x1 %.3f Å · D1 %.3f kcal/mol · Z1 %.3f · χ %.3f", q.r1, q.theta0, q.x1, q.D1, q.Z1, q.Xi);
    t.description = b;
    t.source = "UFF";
    d.types.push_back(t);
  }
  return d;
}

bool qeq_parameters(int z, double& chi, double& J, double& radius) {
  for (const auto& q : kUff) {
    if (element_from_symbol(label_symbol(q.label)) != z) continue;
    chi = q.Xi;
    J = 2 * q.hard;   // the table holds the hardness η = J / 2
    radius = q.radius;
    return true;
  }
  return false;
}

int uff_label_count() { return kUffCount; }
bool uff_has_label(const std::string& l) { return by_label(l) != nullptr; }

std::vector<std::string> uff_types(const System& s, std::vector<std::string>* why) {
  const Perception p = perceive(s);
  const auto A = type_atoms(s, p);
  std::vector<std::string> out;
  if (why) why->clear();
  for (const auto& a : A) {
    out.push_back(a.label);
    if (why) why->push_back(a.why);
  }
  return out;
}

ForceField assign_uff(const System& s, const UffOptions& o) {
  const size_t n = s.atoms.size();
  const Perception p = perceive(s);
  const auto A0 = type_atoms(s, p);
  ForceField ff;
  ff.name = "UFF (Rappé et al. 1992)";
  ff.mixing = "geometric";
  ff.lj14 = 1.0;
  ff.coul14 = 1.0;
  std::vector<const UffParam*> P(n, nullptr);
  std::map<std::string, int> tindex;
  int fallback = 0;
  auto A = A0;
  int by_hand = 0;
  for (size_t i = 0; i < n && i < o.labels.size(); ++i) {
    if (o.labels[i].empty() || o.labels[i] == A[i].label) continue;
    const UffParam* q = by_label(o.labels[i]);
    if (!q) throw FieldError("atom " + std::to_string(i + 1) + ": " + o.labels[i] + " is not a UFF label");
    A[i].label = q->label;
    A[i].why = "set by hand";
    ++by_hand;
  }
  for (size_t i = 0; i < n; ++i) {
    if (A[i].label.empty())
      throw FieldError("atom " + std::to_string(i + 1) + " has no element UFF covers (Z = " + std::to_string(A[i].z) + ")");
    P[i] = by_label(A[i].label);
    ff.atom_type.push_back(A[i].label);
    ff.why.push_back(A[i].why);
    if (A[i].why.find("closest") != std::string::npos) ++fallback;
    auto it = tindex.find(A[i].label);
    if (it == tindex.end()) {
      it = tindex.emplace(A[i].label, int(ff.type_names.size())).first;
      ff.type_names.push_back(A[i].label);
      ff.lj.push_back({P[i]->D1, P[i]->x1 / std::pow(2.0, 1.0 / 6.0)});
    }
    ff.type_index.push_back(it->second);
    ff.mass.push_back(element(A[i].z).mass);
    ff.charge.push_back(o.keep_charges && s.has_charges ? s.atoms[i].charge : 0.0);
  }
  auto bo = [&](uint32_t a, uint32_t b) {
    for (size_t k = 0; k < p.nb[a].size(); ++k)
      if (p.nb[a][k] == b) return bond_order_value(p, a, k);
    return 1.0;
  };

  // bonds
  for (uint32_t i = 0; i < n; ++i)
    for (uint32_t j : p.nb[i]) {
      if (j <= i) continue;
      const double r0 = rest_length(bo(i, j), *P[i], *P[j]);
      ff.bonds.push_back({i, j, kG * P[i]->Z1 * P[j]->Z1 / (r0 * r0 * r0), r0});
    }

  // rings of three and four for the sp2 special cases
  std::vector<char> in3(n, 0), in4(n, 0);
  for (const auto& r : p.rings) {
    if (r.size() == 3) for (uint32_t a : r) in3[a] = 1;
    if (r.size() == 4) for (uint32_t a : r) in4[a] = 1;
  }
  int skipped_centres = 0;
  // angles
  for (uint32_t j = 0; j < n; ++j) {
    const auto& nb = p.nb[j];
    if (nb.size() < 2) continue;
    if (nb.size() > 6) { ++skipped_centres; continue; }
    auto add = [&](uint32_t i, uint32_t k, int order, double theta0) {
      const double K = angle_k(theta0, bo(i, j), bo(j, k), *P[i], *P[j], *P[k]);
      ff.angles_x.push_back({i, j, k, order == 0 ? 3 : 10 + order, K, theta0});
    };
    const double t0 = P[j]->theta0 * kDeg;
    if (A[j].hyb == 5 && nb.size() == 5) {
      // trigonal bipyramid: the most nearly opposite pair is axial
      double most = 3;
      size_t ax1 = 0, ax2 = 1;
      for (size_t x = 0; x < 5; ++x)
        for (size_t y = x + 1; y < 5; ++y) {
          const double c = cos_angle(s, j, nb[x], nb[y]);
          if (c < most) { most = c; ax1 = x; ax2 = y; }
        }
      for (size_t x = 0; x < 5; ++x)
        for (size_t y = x + 1; y < 5; ++y) {
          const bool ax_x = x == ax1 || x == ax2, ax_y = y == ax1 || y == ax2;
          if (ax_x && ax_y) add(nb[x], nb[y], 2, t0);
          else if (!ax_x && !ax_y) add(nb[x], nb[y], 3, t0);
          else add(nb[x], nb[y], 0, t0);
        }
      continue;
    }
    for (size_t x = 0; x < nb.size(); ++x)
      for (size_t y = x + 1; y < nb.size(); ++y) {
        const uint32_t i = nb[x], k = nb[y];
        int order = 0;
        double th = t0;
        switch (A[j].hyb) {
          case 1: order = 1; break;
          case 2:
            order = 3;
            if (in3[j]) {
              if (in3[i] != in3[k]) { order = 0; th = 150 * kDeg; }
              else if (in3[i] && in3[k]) { order = 0; th = 60 * kDeg; }
            } else if (in4[j]) {
              if (in4[i] != in4[k]) { order = 0; th = 135 * kDeg; }
              else if (in4[i] && in4[k]) { order = 0; th = 90 * kDeg; }
            }
            break;
          case 4: case 6: order = 4; break;
          default: order = 0; break;
        }
        add(i, k, order, th);
      }
  }

  // torsions about bonds between sp2 / sp3 atoms that are not terminal and not in a triple bond
  auto in_triple = [&](uint32_t a) {
    for (size_t k = 0; k < p.nb[a].size(); ++k)
      if (p.order[a][k] == 3 && !p.arom_bond[a][k]) return true;
    return false;
  };
  auto sp23 = [&](uint32_t a) { return A[a].hyb == 2 || A[a].hyb == 3; };
  for (uint32_t b = 0; b < n; ++b)
    for (uint32_t c : p.nb[b]) {
      if (c <= b) continue;
      if (p.nb[b].size() < 2 || p.nb[c].size() < 2 || in_triple(b) || in_triple(c) || !sp23(b) || !sp23(c)) continue;
      const double order = bo(b, c);
      const bool single = std::fabs(order - 1.0) < 1e-9;
      std::vector<std::array<uint32_t, 4>> quads;
      for (uint32_t a : p.nb[b])
        if (a != c)
          for (uint32_t d : p.nb[c])
            if (d != b && d != a) quads.push_back({a, b, c, d});
      if (quads.empty()) continue;
      for (const auto& q : quads) {
        const bool end_sp2 = A[q[0]].hyb == 2 || A[q[3]].hyb == 2;
        double V;
        int nn;
        double cos_term;
        const int zb = A[b].z, zc = A[c].z;
        if (A[b].hyb == 3 && A[c].hyb == 3) {
          V = std::sqrt(P[b]->V1 * P[c]->V1);
          nn = 3;
          cos_term = -1;
          if (single && group16(zb) && group16(zc)) {
            const double v2 = zb == 8 ? 2.0 : 6.8, v3 = zc == 8 ? 2.0 : 6.8;
            V = std::sqrt(v2 * v3);
            nn = 2;
          }
        } else if (A[b].hyb == 2 && A[c].hyb == 2) {
          V = equation17(order, *P[b], *P[c]);
          nn = 2;
          cos_term = 1;
        } else {
          V = 1.0;
          nn = 6;
          cos_term = 1;
          if (single) {
            if ((A[b].hyb == 3 && group16(zb) && !group16(zc)) || (A[c].hyb == 3 && group16(zc) && !group16(zb))) {
              V = equation17(order, *P[b], *P[c]);
              nn = 2;
              cos_term = -1;
            } else if (end_sp2) {
              V = 2.0;
              nn = 3;
              cos_term = -1;
            }
          }
        }
        V /= double(quads.size());
        if (std::fabs(V) < 1e-12) continue;
        // ½ V [1 − cos(nφ0) cos(nφ)] = v [1 + cos(nφ − δ)] with δ = 0 when cos(nφ0) = −1, π when +1
        ff.dihedrals.push_back({q[0], q[1], q[2], q[3], V / 2, nn, cos_term < 0 ? 0.0 : kPi});
      }
    }

  // inversions
  for (uint32_t c = 0; c < n; ++c) {
    if (p.nb[c].size() != 3) continue;
    const int z = A[c].z;
    const auto& nb = p.nb[c];
    if ((z == 6 || z == 7 || z == 8) && A[c].hyb == 2) {
      bool to_o2 = false;
      if (z == 6)
        for (uint32_t x : nb) to_o2 = to_o2 || (A[x].z == 8 && A[x].hyb == 2);
      ff.inversions.push_back({c, nb[0], nb[1], nb[2], to_o2 ? 50.0 : 6.0, 0.0, 1});
    } else if (z == 15 || z == 33 || z == 51 || z == 83) {
      const double w0 = (z == 15 ? 84.4339 : z == 33 ? 86.9735 : z == 51 ? 87.7047 : 90.0) * kDeg;
      const double c2 = 1, c1 = -4 * std::cos(w0), c0 = -(c1 * std::cos(w0) + c2 * std::cos(2 * w0));
      ff.inversions.push_back({c, nb[0], nb[1], nb[2], 22.0 / (c0 + c1 + c2), w0, 2});
    }
  }

  // exclusions: 1-2 and 1-3 pairs; 1-4 pairs interact in full
  ff.excluded.assign(n, {});
  for (uint32_t i = 0; i < n; ++i) {
    for (uint32_t j : p.nb[i]) {
      ff.excluded[i].push_back(j);
      for (uint32_t k : p.nb[j])
        if (k != i) ff.excluded[i].push_back(k);
    }
    std::sort(ff.excluded[i].begin(), ff.excluded[i].end());
    ff.excluded[i].erase(std::unique(ff.excluded[i].begin(), ff.excluded[i].end()), ff.excluded[i].end());
  }

  std::map<std::string, int> counts;
  for (const auto& t : ff.atom_type) counts[t]++;
  std::string c;
  for (const auto& [t, k] : counts) c += (c.empty() ? "" : ", ") + std::to_string(k) + " " + t;
  ff.notes.push_back("typed " + std::to_string(n) + " atoms: " + c);
  if (by_hand) ff.notes.push_back(std::to_string(by_hand) + " atoms have UFF labels set by hand");
  if (fallback) ff.notes.push_back(std::to_string(fallback) + " atoms took the closest UFF label of their element");
  if (skipped_centres) ff.notes.push_back(std::to_string(skipped_centres) + " centres with more than six neighbours have no angle terms");
  ff.notes.push_back(std::to_string(ff.bonds.size()) + " bonds, " + std::to_string(ff.angles_x.size()) + " angles, " +
                     std::to_string(ff.dihedrals.size()) + " torsion terms, " + std::to_string(ff.inversions.size()) + " inversions");
  if (o.qeq) {
    QEqReport qr;
    ff.charge = qeq_charges(s, QEqOptions{}, &qr);
    ff.notes.push_back(qr.notes.front());
  } else {
    ff.notes.push_back(o.keep_charges && s.has_charges ? "charges from the structure" : "no charges (UFF clean-up)");
  }
  for (const auto& note : p.notes) ff.notes.push_back(note);
  return ff;
}

bool uff_vdw(int z, double& x, double& d) {
  if (z <= 0) return false;
  const std::string sym = element(z).symbol;
  for (const auto& p : kUff)
    if (label_element(p.label) == sym) { x = p.x1, d = p.D1; return true; }
  return false;
}

}  // namespace caps
