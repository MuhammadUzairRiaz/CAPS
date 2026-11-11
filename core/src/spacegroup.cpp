// CAPS space groups (see caps/spacegroup.hpp).
#include "caps/spacegroup.hpp"

#include <array>
#include <cstdlib>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iterator>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/crystal.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {

struct Row { const char* key; const char* hm; const char* hall; };
const Row kRows[] = {
#include "spacegroups.inc"
};

constexpr int kT = 24;   // translations in 1/24 (covers 1/2, 1/3, 1/4, 1/6, 1/8 and the 1/12 origin shifts)

struct IOp {
  int R[3][3];
  int t[3];
  bool operator<(const IOp& o) const {
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) if (R[i][j] != o.R[i][j]) return R[i][j] < o.R[i][j];
    for (int i = 0; i < 3; ++i) if (t[i] != o.t[i]) return t[i] < o.t[i];
    return false;
  }
};

int mod(int a) { return ((a % kT) + kT) % kT; }

IOp mul(const IOp& a, const IOp& b) {
  IOp c{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      c.R[i][j] = 0;
      for (int k = 0; k < 3; ++k) c.R[i][j] += a.R[i][k] * b.R[k][j];
    }
    int v = a.t[i];
    for (int k = 0; k < 3; ++k) v += a.R[i][k] * b.t[k];
    c.t[i] = mod(v);
  }
  return c;
}

IOp identity() {
  IOp o{};
  for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) o.R[i][j] = i == j; o.t[i] = 0; }
  return o;
}

using M3 = std::array<std::array<int, 3>, 3>;

// Rotation matrices of Hall's table (proper rotations of order N about an axis)
M3 rotation(int n, char axis, char principal) {
  const M3 I{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  if (n == 1) return I;
  if (axis == '*') return M3{{{0, 0, 1}, {1, 0, 0}, {0, 1, 0}}};   // 3 along the body diagonal
  if (axis == '\'' || axis == '"') {   // 2-fold along a−b (') or a+b (") relative to the principal axis
    const bool dbl = axis == '"';
    if (principal == 'x') return dbl ? M3{{{-1, 0, 0}, {0, 0, 1}, {0, 1, 0}}} : M3{{{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}}};
    if (principal == 'y') return dbl ? M3{{{0, 0, 1}, {0, -1, 0}, {1, 0, 0}}} : M3{{{0, 0, -1}, {0, -1, 0}, {-1, 0, 0}}};
    return dbl ? M3{{{0, 1, 0}, {1, 0, 0}, {0, 0, -1}}} : M3{{{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}}};
  }
  if (axis == 'x') {
    switch (n) {
      case 2: return M3{{{1, 0, 0}, {0, -1, 0}, {0, 0, -1}}};
      case 3: return M3{{{1, 0, 0}, {0, 0, -1}, {0, 1, -1}}};
      case 4: return M3{{{1, 0, 0}, {0, 0, -1}, {0, 1, 0}}};
      case 6: return M3{{{1, 0, 0}, {0, 1, -1}, {0, 1, 0}}};
    }
  } else if (axis == 'y') {
    switch (n) {
      case 2: return M3{{{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}}};
      case 3: return M3{{{-1, 0, 1}, {0, 1, 0}, {-1, 0, 0}}};
      case 4: return M3{{{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}}};
      case 6: return M3{{{0, 0, 1}, {0, 1, 0}, {-1, 0, 1}}};
    }
  } else {
    switch (n) {
      case 2: return M3{{{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}}};
      case 3: return M3{{{0, -1, 0}, {1, -1, 0}, {0, 0, 1}}};
      case 4: return M3{{{0, -1, 0}, {1, 0, 0}, {0, 0, 1}}};
      case 6: return M3{{{1, -1, 0}, {1, 0, 0}, {0, 0, 1}}};
    }
  }
  throw std::invalid_argument("Hall symbol: no " + std::to_string(n) + "-fold rotation about " + std::string(1, axis));
}

std::vector<IOp> hall_int_ops(const std::string& hall_in) {
  std::string hall = hall_in;
  for (char& c : hall) c = char(std::tolower(static_cast<unsigned char>(c)));
  // origin shift (vx vy vz), in twelfths
  int shift[3] = {0, 0, 0};
  if (const auto p = hall.find('('); p != std::string::npos) {
    std::istringstream v(hall.substr(p + 1, hall.find(')', p) - p - 1));
    for (int& x : shift) v >> x;
    hall = hall.substr(0, p);
  }
  std::istringstream in(hall);
  std::string tok;
  in >> tok;
  if (tok.empty()) throw std::invalid_argument("empty Hall symbol");
  const bool centric = tok[0] == '-';
  const char lat = centric ? tok[1] : tok[0];
  std::vector<std::array<int, 3>> centring = {{0, 0, 0}};
  const int h = kT / 2, t3 = kT / 3;
  switch (lat) {
    case 'p': break;
    case 'a': centring.push_back({0, h, h}); break;
    case 'b': centring.push_back({h, 0, h}); break;
    case 'c': centring.push_back({h, h, 0}); break;
    case 'i': centring.push_back({h, h, h}); break;
    case 'r': centring.push_back({2 * t3, t3, t3}); centring.push_back({t3, 2 * t3, 2 * t3}); break;
    case 's': centring.push_back({t3, t3, 2 * t3}); centring.push_back({2 * t3, 2 * t3, t3}); break;
    case 't': centring.push_back({t3, 2 * t3, t3}); centring.push_back({2 * t3, t3, 2 * t3}); break;
    case 'f': centring.push_back({0, h, h}); centring.push_back({h, 0, h}); centring.push_back({h, h, 0}); break;
    default: throw std::invalid_argument("Hall symbol: unknown lattice " + std::string(1, lat));
  }
  std::vector<IOp> gens;
  for (const auto& c : centring) {
    IOp o = identity();
    for (int i = 0; i < 3; ++i) o.t[i] = c[i];
    gens.push_back(o);
  }
  if (centric) {
    IOp o = identity();
    for (int i = 0; i < 3; ++i) o.R[i][i] = -1;
    gens.push_back(o);
  }
  int position = 0, previous = 0;
  char principal = 'z';
  while (in >> tok) {
    size_t k = 0;
    const bool improper = tok[k] == '-';
    if (improper) ++k;
    if (k >= tok.size() || !std::isdigit(static_cast<unsigned char>(tok[k]))) throw std::invalid_argument("Hall symbol: rotation expected in '" + tok + "'");
    const int n = tok[k++] - '0';
    char axis = 0;
    int screw = 0;
    int t[3] = {0, 0, 0};
    for (; k < tok.size(); ++k) {
      const char c = tok[k];
      if (std::isdigit(static_cast<unsigned char>(c))) screw = c - '0';
      else if (c == 'x' || c == 'y' || c == 'z' || c == '\'' || c == '"' || c == '*') axis = c;
      else if (c == 'a') t[0] += h;
      else if (c == 'b') t[1] += h;
      else if (c == 'c') t[2] += h;
      else if (c == 'n') { t[0] += h; t[1] += h; t[2] += h; }
      else if (c == 'u') t[0] += kT / 4;
      else if (c == 'v') t[1] += kT / 4;
      else if (c == 'w') t[2] += kT / 4;
      else if (c == 'd') { t[0] += kT / 4; t[1] += kT / 4; t[2] += kT / 4; }
      else throw std::invalid_argument("Hall symbol: unknown '" + std::string(1, c) + "' in '" + tok + "'");
    }
    if (!axis) {   // Hall's default axes
      if (position == 0) axis = 'z';
      else if (position == 1 && n == 2) axis = (previous == 2 || previous == 4) ? 'x' : '\'';
      else if (position == 2 && n == 3) axis = '*';
      else axis = 'z';
    }
    if (position == 0 && (axis == 'x' || axis == 'y' || axis == 'z')) principal = axis;
    const M3 R = rotation(n, axis, principal);
    IOp o{};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) o.R[i][j] = improper ? -R[i][j] : R[i][j];
    if (screw && n > 1) {   // screw along the rotation axis
      const int s = kT * screw / n;
      if (axis == 'x') t[0] += s; else if (axis == 'y') t[1] += s; else t[2] += s;
    }
    for (int i = 0; i < 3; ++i) o.t[i] = mod(t[i]);
    gens.push_back(o);
    previous = n;
    ++position;
  }
  // closure
  std::set<IOp> group = {identity()};
  std::vector<IOp> frontier = {identity()};
  while (!frontier.empty()) {
    std::vector<IOp> next;
    for (const auto& a : frontier)
      for (const auto& g : gens) {
        const IOp c = mul(g, a);
        if (group.insert(c).second) next.push_back(c);
      }
    frontier = std::move(next);
    if (group.size() > 400) throw std::invalid_argument("Hall symbol " + hall_in + " does not close");
  }
  std::vector<IOp> ops(group.begin(), group.end());
  // the origin shift V: S' = V S V⁻¹, t' = t + (I − R) v
  if (shift[0] || shift[1] || shift[2])
    for (auto& o : ops)
      for (int i = 0; i < 3; ++i) {
        int v = 0;
        for (int k = 0; k < 3; ++k) v += ((i == k) - o.R[i][k]) * shift[k] * (kT / 12);
        o.t[i] = mod(o.t[i] + v);
      }
  // identity first, then the rest as generated
  std::stable_sort(ops.begin(), ops.end(), [](const IOp& a, const IOp& b) {
    auto isid = [](const IOp& o) { for (int i = 0; i < 3; ++i) { if (o.t[i]) return false; for (int j = 0; j < 3; ++j) if (o.R[i][j] != (i == j)) return false; } return true; };
    return isid(a) && !isid(b);
  });
  return ops;
}

std::string squash(std::string s) {
  std::string out;
  for (char c : s) if (!std::isspace(static_cast<unsigned char>(c)) && c != '_') out += char(std::tolower(static_cast<unsigned char>(c)));
  return out;
}

Vec3 wrap01(Vec3 f) {
  for (auto& x : f) { x -= std::floor(x); if (x >= 1 - 1e-9) x = 0; }
  return f;
}

Cell cell_from(double a, double b, double c, double al, double be, double ga) {
  constexpr double d = M_PI / 180;
  al *= d; be *= d; ga *= d;
  Cell cell;
  cell.a = {a, 0, 0};
  cell.b = {b * std::cos(ga), b * std::sin(ga), 0};
  const double cx = c * std::cos(be), cy = c * (std::cos(al) - std::cos(be) * std::cos(ga)) / std::sin(ga);
  cell.c = {cx, cy, std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy))};
  return cell;
}

double frac_distance(const Cell& cell, Vec3 df) {
  for (auto& x : df) x -= std::round(x);
  return norm(cell.a * df[0] + cell.b * df[1] + cell.c * df[2]);
}

Vec3 apply_op(const SymOp& op, const Vec3& f) {
  Vec3 r;
  for (int i = 0; i < 3; ++i) r[i] = op.R[i][0] * f[0] + op.R[i][1] * f[1] + op.R[i][2] * f[2] + op.t[i];
  return r;
}

void finish_types(System& s) {
  std::map<int, int> type_of;
  s.types.clear();
  for (auto& a : s.atoms) {
    auto it = type_of.find(a.element);
    if (it == type_of.end()) {
      it = type_of.emplace(a.element, int(type_of.size()) + 1).first;
      TypeInfo ti;
      ti.type = it->second;
      ti.mass = element(a.element).mass;
      ti.label = element(a.element).symbol;
      s.types.push_back(ti);
    }
    a.type = it->second;
  }
}

}  // namespace

Cell cell_parameters(double a, double b, double c, double alpha, double beta, double gamma) { return cell_from(a, b, c, alpha, beta, gamma); }

// ---------------------------------------------------------------- operations

SymOp parse_symop(std::string s) {
  SymOp op;
  for (auto& r : op.R) for (double& x : r) x = 0;
  s.erase(std::remove_if(s.begin(), s.end(), [](char c) { return std::isspace(static_cast<unsigned char>(c)) || c == '\'' || c == '"'; }), s.end());
  for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  std::vector<std::string> comp;
  std::stringstream ss(s);
  for (std::string x; std::getline(ss, x, ',');) comp.push_back(x);
  if (comp.size() != 3) throw CrystalError("symmetry operation '" + s + "' does not have three components");
  auto fraction = [](const std::string& f) {
    const auto p = f.find('/');
    return p == std::string::npos ? std::stod(f) : std::stod(f.substr(0, p)) / std::stod(f.substr(p + 1));
  };
  for (int row = 0; row < 3; ++row) {
    const std::string& c = comp[size_t(row)];
    size_t i = 0;
    while (i < c.size()) {
      double sign = 1;
      if (c[i] == '+' || c[i] == '-') { sign = c[i] == '-' ? -1 : 1; ++i; }
      size_t j = i;
      while (j < c.size() && c[j] != '+' && c[j] != '-') ++j;
      const std::string term = c.substr(i, j - i);
      i = j;
      if (term.empty()) continue;
      const auto xyz = term.find_first_of("xyz");
      if (xyz != std::string::npos) {
        std::string coef = term.substr(0, xyz);
        if (!coef.empty() && coef.back() == '*') coef.pop_back();
        op.R[row][term[xyz] - 'x'] += sign * (coef.empty() ? 1.0 : fraction(coef));
      } else {
        op.t[row] += sign * fraction(term);
      }
    }
  }
  return op;
}

std::string symop_string(const SymOp& op) {
  std::string out;
  for (int i = 0; i < 3; ++i) {
    std::string c;
    for (int k = 0; k < 3; ++k) {
      const double v = op.R[i][k];
      if (std::fabs(v) < 1e-9) continue;
      c += (v < 0 ? "-" : (c.empty() ? "" : "+"));
      c += char('x' + k);
    }
    const double t = op.t[i] - std::floor(op.t[i] + 1e-9);
    if (t > 1e-6) {
      int best = 0;
      for (int d : {2, 3, 4, 6, 8, 12, 24}) if (std::fabs(t * d - std::round(t * d)) < 1e-6) { best = d; break; }
      if (best) {
        c += "+" + std::to_string(int(std::lround(t * best))) + "/" + std::to_string(best);
      } else {   // not a rational of the usual denominators (an origin shift): written as a decimal
        char v[32];
        std::snprintf(v, sizeof v, "+%.6f", t);
        c += v;
      }
    }
    out += (i ? "," : "") + (c.empty() ? std::string("0") : c);
  }
  return out;
}

const std::vector<SpaceGroupSetting>& space_group_settings() {
  static const std::vector<SpaceGroupSetting> all = [] {
    std::vector<SpaceGroupSetting> v;
    for (const auto& r : kRows) v.push_back({r.key, std::atoi(r.key), r.hm, r.hall});
    return v;
  }();
  return all;
}

const SpaceGroupSetting* find_space_group(const std::string& query) {
  const auto& all = space_group_settings();
  const std::string q = squash(query);
  if (q.empty()) return nullptr;
  for (const auto& s : all) if (squash(s.key) == q) return &s;
  if (std::all_of(q.begin(), q.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) {
    const int n = std::atoi(q.c_str());
    const SpaceGroupSetting* first = nullptr;
    for (const auto& s : all)
      if (s.number == n) {
        if (!first) first = &s;
        if (s.key.size() > 2 && s.key.substr(s.key.size() - 2) == ":2") return &s;   // origin choice 2 when there are two
      }
    return first;
  }
  const SpaceGroupSetting* plain = nullptr;
  for (const auto& s : all) {
    std::string hm = squash(s.hm);
    const auto colon = hm.find(':');
    const std::string bare = colon == std::string::npos ? hm : hm.substr(0, colon);
    if (hm == q) return &s;
    if (bare == q && !plain) {
      plain = &s;
      if (colon != std::string::npos && hm.substr(colon) == ":1") {   // prefer origin choice 2 and hexagonal axes
        for (const auto& t : all) if (t.number == s.number && squash(t.hm) == bare + ":2") return &t;
      }
    }
  }
  if (plain) return plain;
  // short monoclinic symbols ("P 21/c" for "P 1 21/c 1"): the first setting, unique axis b
  for (const auto& s : all) {
    if (s.number < 3 || s.number > 15) continue;
    std::istringstream in(s.hm.substr(0, s.hm.find(':')));
    std::vector<std::string> t{std::istream_iterator<std::string>(in), std::istream_iterator<std::string>()};
    if (t.size() != 4) continue;
    std::string shortened = t[0];
    for (size_t k = 1; k < 4; ++k) if (t[k] != "1") shortened += t[k];
    if (squash(shortened) == q) return &s;
  }
  return nullptr;
}

std::vector<SymOp> hall_operations(const std::string& hall) {
  static std::mutex m;
  static std::map<std::string, std::vector<SymOp>> cache;
  std::lock_guard<std::mutex> lock(m);
  if (auto it = cache.find(hall); it != cache.end()) return it->second;
  std::vector<SymOp> out;
  for (const auto& o : hall_int_ops(hall)) {
    SymOp s;
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) s.R[i][j] = o.R[i][j]; s.t[i] = double(o.t[i]) / kT; }
    out.push_back(s);
  }
  return cache[hall] = out;
}

std::string crystal_system(int n) {
  if (n <= 2) return "triclinic";
  if (n <= 15) return "monoclinic";
  if (n <= 74) return "orthorhombic";
  if (n <= 142) return "tetragonal";
  if (n <= 167) return "trigonal";
  if (n <= 194) return "hexagonal";
  return "cubic";
}

// ---------------------------------------------------------------- crystals

System build_crystal(const CrystalSpec& spec, CrystalReport* report) {
  const SpaceGroupSetting* sg = find_space_group(spec.space_group);
  if (!sg) throw CrystalError("unknown space group '" + spec.space_group + "'");
  if (spec.a <= 0 || spec.b <= 0 || spec.c <= 0) throw CrystalError("the cell lengths must be positive");
  const auto ops = hall_operations(sg->hall);
  System s;
  s.title = spec.title.empty() ? sg->hm : spec.title;
  s.source_format = "crystal";
  s.cell = cell_from(spec.a, spec.b, spec.c, spec.alpha, spec.beta, spec.gamma);
  if (!(s.cell.volume() > 1e-6)) throw CrystalError("the cell angles do not make a cell");
  CrystalReport R;
  R.key = sg->key; R.hm = sg->hm; R.hall = sg->hall; R.number = sg->number; R.system = crystal_system(sg->number);
  R.operations = int(ops.size());
  std::vector<Vec3> frac;
  std::vector<int> zs;
  std::vector<std::string> labels;
  int overlaps = 0;
  for (size_t k = 0; k < spec.sites.size(); ++k) {
    const auto& site = spec.sites[k];
    if (site.element <= 0) throw CrystalError("site " + site.label + ": no element");
    int mult = 0;
    const size_t start = frac.size();
    for (const auto& op : ops) {
      const Vec3 f = wrap01(apply_op(op, site.frac));
      bool dup = false;
      for (size_t q = 0; q < frac.size() && !dup; ++q)
        if (frac_distance(s.cell, f - frac[q]) < std::max(spec.tolerance, 1e-4)) {
          dup = true;
          if (q < start) ++overlaps;   // lands on another site's atom
        }
      if (dup) continue;
      frac.push_back(f);
      zs.push_back(site.element);
      labels.push_back(site.label.empty() ? element(site.element).symbol : site.label);
      ++mult;
    }
    R.multiplicity.push_back(mult);
  }
  if (overlaps) R.notes.push_back(std::to_string(overlaps) + " image(s) fall on atoms of another site: check the coordinates");
  for (size_t i = 0; i < frac.size(); ++i) {
    Atom a;
    a.id = int64_t(i + 1);
    a.mol = 1;
    a.element = zs[i];
    a.name = labels[i];
    a.pos = s.cell.to_cartesian(frac[i]);
    s.atoms.push_back(a);
  }
  finish_types(s);
  s.has_mol = false;   // molecules follow the perceived bonds, as for a CIF; mol = 1 above is only a placeholder
  s.unwrapped = true;
  R.atoms_per_cell = s.atoms.size();
  R.volume = s.cell.volume();
  s.symmetry = {sg->hm, sg->number, ops};
  if (spec.supercell[0] * spec.supercell[1] * spec.supercell[2] > 1) s = supercell(s, spec.supercell[0], spec.supercell[1], spec.supercell[2]);
  s.bonds = crystal_bonds(s);
  char b[200];
  std::snprintf(b, sizeof b, "%s (No. %d, %s) · %d operations · %zu atoms per cell · V %.2f Å³", sg->hm.c_str(), sg->number, R.system.c_str(), R.operations,
                R.atoms_per_cell, R.volume);
  s.notes.insert(s.notes.begin(), b);
  if (report) *report = R;
  return s;
}

CrystalSpec symmetrize_sites(const CrystalSpec& spec, double snap, int* moved) {
  const SpaceGroupSetting* sg = find_space_group(spec.space_group);
  if (!sg) throw CrystalError("unknown space group '" + spec.space_group + "'");
  const auto ops = hall_operations(sg->hall);
  const Cell cell = cell_from(spec.a, spec.b, spec.c, spec.alpha, spec.beta, spec.gamma);
  CrystalSpec out = spec;
  int count = 0;
  for (auto& site : out.sites) {
    const Vec3 start = site.frac;
    for (int pass = 0; pass < 3; ++pass) {
      // the images that land near the site itself are its site-symmetry group; their mean is the special position
      Vec3 sum{0, 0, 0};
      int k = 0;
      for (const auto& op : ops) {
        Vec3 d = apply_op(op, site.frac) - site.frac;
        for (int i = 0; i < 3; ++i) d[i] -= std::round(d[i]);
        if (frac_distance(cell, d) < snap) sum = sum + d, ++k;
      }
      if (k > 0) site.frac = site.frac + sum * (1.0 / k);
    }
    for (int i = 0; i < 3; ++i) {
      site.frac[i] -= std::floor(site.frac[i]);
      if (std::fabs(site.frac[i]) < 1e-9 || std::fabs(site.frac[i] - 1) < 1e-9) site.frac[i] = 0;
    }
    if (frac_distance(cell, site.frac - start - Vec3{std::round(site.frac[0] - start[0]), std::round(site.frac[1] - start[1]), std::round(site.frac[2] - start[2])}) > 1e-6) ++count;
  }
  if (moved) *moved = count;
  return out;
}

System supercell(const System& s, int nx, int ny, int nz) {
  if (!s.cell.valid()) throw CrystalError("a supercell needs a cell");
  nx = std::max(1, nx); ny = std::max(1, ny); nz = std::max(1, nz);
  System out = s;
  out.atoms.clear();
  out.bonds.clear();
  out.velocities.clear();
  const size_t na = s.atoms.size();
  int64_t max_mol = 0;
  for (const auto& a : s.atoms) max_mol = std::max(max_mol, a.mol);
  out.atoms.reserve(na * size_t(nx) * ny * nz);
  int64_t id = 0, image = 0;
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k, ++image)
        for (const auto& a : s.atoms) {
          Atom b = a;
          b.pos = a.pos + s.cell.a * i + s.cell.b * j + s.cell.c * k;
          b.id = ++id;
          if (a.mol > 0) b.mol = a.mol + image * max_mol;   // each copy its own molecules
          out.atoms.push_back(b);
        }
  out.cell.a = s.cell.a * nx;
  out.cell.b = s.cell.b * ny;
  out.cell.c = s.cell.c * nz;
  // bonds perceived from distances (a crystal, where a small cell can bond an atom to several images of one neighbour) are
  // perceived again in the new cell; a topology (a file's, a builder's) is copied into every image
  if (s.bonds.empty() || !s.bonds_from_file) {
    out.bonds = crystal_bonds(out);
    return out;
  }
  // the structure's own bonds, in every copy: a bond that crosses the cell's wall (its minimum image) joins the
  // neighbouring copy (across the supercell's own wall, periodically), so chains and networks stay as they were
  const int n[3] = {nx, ny, nz};
  out.bonds.reserve(s.bonds.size() * size_t(nx) * ny * nz);
  std::vector<std::array<int, 3>> shift(s.bonds.size());
  for (size_t q = 0; q < s.bonds.size(); ++q) {
    const Vec3 fi = s.cell.to_fractional(s.atoms[s.bonds[q].i].pos), fj = s.cell.to_fractional(s.atoms[s.bonds[q].j].pos);
    for (int c = 0; c < 3; ++c) shift[q][size_t(c)] = int(std::lround(fj[size_t(c)] - fi[size_t(c)]));
  }
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k) {
        const int here[3] = {i, j, k};
        const size_t base = size_t((i * ny + j) * nz + k) * na;
        for (size_t q = 0; q < s.bonds.size(); ++q) {
          int t[3];
          for (int c = 0; c < 3; ++c) t[c] = ((here[c] - shift[q][size_t(c)]) % n[c] + n[c]) % n[c];   // the partner lies at j − shift·cell
          const size_t other = size_t((t[0] * ny + t[1]) * nz + t[2]) * na;
          out.bonds.push_back({uint32_t(base + s.bonds[q].i), uint32_t(other + s.bonds[q].j), s.bonds[q].order});
        }
      }
  return out;
}

System primitive_cell(const System& s, char centring) {
  const Vec3 a = s.cell.a, b = s.cell.b, c = s.cell.c;
  Vec3 p[3];
  switch (std::toupper(static_cast<unsigned char>(centring))) {
    case 'A': p[0] = a; p[1] = (b + c) * 0.5; p[2] = (c - b) * 0.5; break;
    case 'B': p[0] = (a + c) * 0.5; p[1] = b; p[2] = (c - a) * 0.5; break;
    case 'C': p[0] = (a + b) * 0.5; p[1] = (b - a) * 0.5; p[2] = c; break;
    case 'I': p[0] = (b + c - a) * 0.5; p[1] = (a + c - b) * 0.5; p[2] = (a + b - c) * 0.5; break;
    case 'F': p[0] = (b + c) * 0.5; p[1] = (a + c) * 0.5; p[2] = (a + b) * 0.5; break;
    case 'R': p[0] = (a * 2 + b + c) * (1.0 / 3); p[1] = (b + c - a) * (1.0 / 3); p[2] = (c - a - b * 2) * (1.0 / 3); break;
    default: return s;
  }
  System out = s;
  out.cell.a = p[0]; out.cell.b = p[1]; out.cell.c = p[2];
  out.atoms.clear();
  out.bonds.clear();
  out.velocities.clear();
  std::vector<Vec3> frac;
  for (const auto& at : s.atoms) {
    const Vec3 f = wrap01(out.cell.to_fractional(at.pos));
    bool dup = false;
    for (size_t q = 0; q < frac.size() && !dup; ++q)
      dup = out.atoms[q].element == at.element && frac_distance(out.cell, f - frac[q]) < 0.05;
    if (dup) continue;
    frac.push_back(f);
    Atom b = at;
    b.pos = out.cell.to_cartesian(f);
    b.id = int64_t(out.atoms.size() + 1);
    out.atoms.push_back(b);
  }
  out.bonds = crystal_bonds(out);
  out.symmetry = {};
  out.notes.insert(out.notes.begin(), "primitive cell of the " + std::string(1, char(std::toupper(static_cast<unsigned char>(centring)))) + "-centred lattice · " +
                                          std::to_string(out.atoms.size()) + " atoms (from " + std::to_string(s.atoms.size()) + ")");
  return out;
}

// ---------------------------------------------------------------- finding the space group

namespace {

// Atoms in fractional coordinates, binned so that the nearest atom of an element to a point is looked for only in the
// point's bin and its neighbours: every bin is at least `reach` Å across between its faces, so an atom closer than
// `reach` lies at most one bin away along each axis.
struct Frame {
  Cell cell;
  std::vector<Vec3> f;
  std::vector<int> el;
  int nb[3] = {1, 1, 1};
  std::vector<uint32_t> start, list;

  Frame(const Cell& c, std::vector<Vec3> frac, std::vector<int> elements, double reach) : cell(c), f(std::move(frac)), el(std::move(elements)) {
    const Vec3 e[3] = {c.a, c.b, c.c};
    const double V = std::fabs(dot(e[0], cross(e[1], e[2])));
    for (int d = 0; d < 3; ++d) {
      const double recip = norm(cross(e[(d + 1) % 3], e[(d + 2) % 3])) / V;   // 1 / spacing of the planes
      nb[d] = int(std::clamp(std::floor(1.0 / std::max(reach * recip, 1e-9)), 1.0, 64.0));
    }
    const size_t cap = std::max<size_t>(8, 8 * f.size());
    while (size_t(nb[0]) * size_t(nb[1]) * size_t(nb[2]) > cap) {
      int* m = std::max_element(nb, nb + 3);
      *m = std::max(1, *m / 2);
    }
    for (auto& x : f) x = wrap01(x);
    const size_t bins = size_t(nb[0]) * size_t(nb[1]) * size_t(nb[2]);
    start.assign(bins + 1, 0);
    std::vector<size_t> home(f.size());
    for (size_t i = 0; i < f.size(); ++i) home[i] = bin(f[i]), ++start[home[i] + 1];
    for (size_t b = 0; b < bins; ++b) start[b + 1] += start[b];
    list.assign(f.size(), 0);
    std::vector<uint32_t> fill(start.begin(), start.end() - 1);
    for (size_t i = 0; i < f.size(); ++i) list[fill[home[i]]++] = uint32_t(i);
  }
  int axis_bin(double x, int d) const { return std::min(nb[d] - 1, std::max(0, int(std::floor(x * nb[d])))); }
  size_t bin(const Vec3& w) const { return (size_t(axis_bin(w[0], 0)) * size_t(nb[1]) + size_t(axis_bin(w[1], 1))) * size_t(nb[2]) + size_t(axis_bin(w[2], 2)); }
  // the atom of element e nearest the point among the neighbouring bins (−1 if there is none there), and its distance;
  // the first one closer than `enough` when that is given
  long nearest(const Vec3& g, int e, double* dist, double enough = -1) const {
    const Vec3 w = wrap01(g);
    int lo[3], cnt[3];
    for (int d = 0; d < 3; ++d) {
      if (nb[d] <= 3) lo[d] = 0, cnt[d] = nb[d];
      else lo[d] = axis_bin(w[d], d) - 1, cnt[d] = 3;
    }
    double best = 1e300;
    long hit = -1;
    for (int x = 0; x < cnt[0]; ++x)
      for (int y = 0; y < cnt[1]; ++y)
        for (int z = 0; z < cnt[2]; ++z) {
          const size_t b = (size_t((lo[0] + x + nb[0]) % nb[0]) * size_t(nb[1]) + size_t((lo[1] + y + nb[1]) % nb[1])) * size_t(nb[2]) +
                           size_t((lo[2] + z + nb[2]) % nb[2]);
          for (uint32_t k = start[b]; k < start[b + 1]; ++k) {
            const uint32_t j = list[k];
            if (el[j] != e) continue;
            const double r = frac_distance(cell, w - f[j]);
            if (r < best) best = r, hit = long(j);
            if (best < enough) {
              if (dist) *dist = best;
              return hit;
            }
          }
        }
    if (dist) *dist = best;
    return hit;
  }
  bool has(const Vec3& g, int e, double tol) const {
    double r = 1e300;
    nearest(g, e, &r, tol);
    return r < tol;
  }
};

using IM = std::array<int, 9>;   // an integer rotation, row by row
const IM kIdentity{1, 0, 0, 0, 1, 0, 0, 0, 1};

IM int_rotation(const SymOp& o) {
  IM m;
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) m[size_t(3 * i + j)] = int(std::lround(o.R[i][j]));
  return m;
}

IM im_mul(const IM& a, const IM& b) {
  IM c{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) c[size_t(3 * i + j)] += a[size_t(3 * i + k)] * b[size_t(3 * k + j)];
  return c;
}

Vec3 rotate(const IM& R, const Vec3& v) {
  Vec3 r;
  for (int i = 0; i < 3; ++i) r[size_t(i)] = R[size_t(3 * i)] * v[0] + R[size_t(3 * i + 1)] * v[1] + R[size_t(3 * i + 2)] * v[2];
  return r;
}

// a setting as the search uses it: its operations, their rotations, and operations whose rotations generate its point group
struct Setting {
  const SpaceGroupSetting* sg = nullptr;
  std::vector<SymOp> ops;
  std::vector<IM> rot;
  std::vector<size_t> gens;
};

// every distinct setting, the highest order first (table order among equals)
const std::vector<Setting>& search_order() {
  static const std::vector<Setting> all = [] {
    std::vector<Setting> v;
    std::set<std::string> seen;
    for (const auto& sg : space_group_settings()) {
      if (!seen.insert(sg.hall).second) continue;
      Setting S;
      S.sg = &sg;
      S.ops = hall_operations(sg.hall);
      std::set<IM> group = {kIdentity};
      for (size_t k = 0; k < S.ops.size(); ++k) {
        S.rot.push_back(int_rotation(S.ops[k]));
        if (group.count(S.rot[k])) continue;
        S.gens.push_back(k);
        group.insert(S.rot[k]);
        for (bool grew = true; grew;) {   // closed under products
          grew = false;
          const std::vector<IM> cur(group.begin(), group.end());
          for (const auto& a : cur)
            for (const auto& b : cur) grew = group.insert(im_mul(a, b)).second || grew;
        }
      }
      v.push_back(std::move(S));
    }
    std::stable_sort(v.begin(), v.end(), [](const Setting& x, const Setting& y) { return x.ops.size() > y.ops.size(); });
    return v;
  }();
  return all;
}

// (R, t) maps every atom onto an atom of its element within tol
bool maps_onto(const Frame& F, const IM& R, const Vec3& t, double tol) {
  for (size_t i = 0; i < F.f.size(); ++i)
    if (!F.has(rotate(R, F.f[i]) + t, F.el[i], tol)) return false;
  return true;
}

// The origin shifts s with (I − R_k) s ≡ d_k (mod 1) for every generator k. The stacked integer matrix A is brought to
// diagonal form by unimodular row and column operations, D = U A V, so that D y ≡ U d (mod 1) with s = V y: a nonzero
// D_t gives the |D_t| solutions y_t = (c_t + k) / D_t, a zero one leaves y_t free (a polar direction, taken as 0).
std::vector<Vec3> origin_solutions(const std::vector<IM>& R, const std::vector<Vec3>& d) {
  const size_t m = 3 * R.size();
  std::vector<std::array<long long, 3>> A(m);
  std::vector<double> c(m);
  for (size_t k = 0; k < R.size(); ++k)
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) A[3 * k + size_t(i)][size_t(j)] = (i == j) - R[k][size_t(3 * i + j)];
      c[3 * k + size_t(i)] = d[k][size_t(i)];
    }
  long long V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  int rank = 0;
  for (int t = 0; t < 3 && size_t(t) < m; ++t) {
    bool pivot = false;
    for (;;) {
      size_t pi = m;
      int pj = -1;
      long long pv = 0;
      for (size_t i = size_t(t); i < m; ++i)
        for (int j = t; j < 3; ++j)
          if (A[i][size_t(j)] && (!pv || std::llabs(A[i][size_t(j)]) < pv)) pv = std::llabs(A[i][size_t(j)]), pi = i, pj = j;
      if (pi == m) break;
      pivot = true;
      std::swap(A[size_t(t)], A[pi]);
      std::swap(c[size_t(t)], c[pi]);
      if (pj != t) {
        for (auto& row : A) std::swap(row[size_t(t)], row[size_t(pj)]);
        for (auto& row : V) std::swap(row[t], row[pj]);
      }
      const long long p = A[size_t(t)][size_t(t)];
      bool clean = true;
      for (size_t i = size_t(t) + 1; i < m; ++i) {
        const long long q = A[i][size_t(t)] / p;
        if (q) {
          for (int j = 0; j < 3; ++j) A[i][size_t(j)] -= q * A[size_t(t)][size_t(j)];
          c[i] -= double(q) * c[size_t(t)];
        }
        clean = clean && !A[i][size_t(t)];
      }
      for (int j = t + 1; j < 3; ++j) {
        const long long q = A[size_t(t)][size_t(j)] / p;
        if (q) {
          for (auto& row : A) row[size_t(j)] -= q * row[size_t(t)];
          for (auto& row : V) row[j] -= q * row[t];
        }
        clean = clean && !A[size_t(t)][size_t(j)];
      }
      if (clean) break;
    }
    if (!pivot) break;
    rank = t + 1;
  }
  long long D[3] = {1, 1, 1};
  size_t total = 1;
  for (int t = 0; t < rank; ++t) D[t] = A[size_t(t)][size_t(t)], total *= size_t(std::llabs(D[t]));
  std::vector<Vec3> out;
  if (total > 512) return out;   // no crystallographic operations give this many
  for (size_t idx = 0; idx < total; ++idx) {
    double y[3] = {0, 0, 0};
    size_t rest = idx;
    for (int t = 0; t < rank; ++t) {
      const size_t a = size_t(std::llabs(D[t]));
      y[t] = (c[size_t(t)] + double(rest % a)) / double(D[t]);
      rest /= a;
    }
    Vec3 s;
    for (int i = 0; i < 3; ++i) s[size_t(i)] = double(V[i][0]) * y[0] + double(V[i][1]) * y[1] + double(V[i][2]) * y[2];
    // the w come from atom positions, so a zero component can come out as 0.9999997: taken as 0, so that the origin
    // and the shifted operations read 0, not 1 or x+0.999999
    for (auto& x : s) if (std::fabs(x - std::round(x)) < 1e-6) x = 0;
    out.push_back(wrap01(s));
  }
  return out;
}

struct Hit {
  const Setting* set = nullptr;
  Vec3 origin{0, 0, 0};
  double shift = 1e300;   // Å from the structure's own origin
};

constexpr size_t kChoiceBudget = 512;     // translation choices for a setting's generators
constexpr size_t kOriginBudget = 20000;   // origin shifts tried per setting
constexpr size_t kMaxOps = 6144;          // operations a supercell's group may list in its own axes

// The highest-order setting whose operations, moved to some origin s, all map the structure onto itself. Every
// operation (R, w) of the structure takes the anchor atom (the rarest element) onto an atom of its element, so
// W(R) = { f_j − R f_anchor } that fit is every translation that goes with R. A setting fits at origin s when each
// operation (R, t) has t + (I − R) s in W(R); s is solved from the setting's generators for each choice of their
// translations. Of the settings of the highest order that fit, the smallest shift wins (a structure already in a
// standard setting keeps it: F d -3 m:2 at s = 0 rather than :1 at 1/8 1/8 1/8). `rep` limits the rotations to those
// that keep a supercell of rep[0] × rep[1] × rep[2] of this cell (D⁻¹ R D integral).
Hit search(const Frame& F, double tol, const std::array<int, 3>& rep) {
  const Cell& cell = F.cell;
  const size_t n = F.f.size();
  double G[3][3];
  const Vec3 e[3] = {cell.a, cell.b, cell.c};
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) G[i][j] = dot(e[i], e[j]);
  const double gscale = std::max({G[0][0], G[1][1], G[2][2]});
  std::map<IM, bool> allowed;   // RᵀGR = G, and the supercell's lattice kept
  auto ok = [&](const IM& R) {
    if (auto it = allowed.find(R); it != allowed.end()) return it->second;
    bool good = true;
    for (int i = 0; i < 3 && good; ++i)
      for (int j = 0; j < 3 && good; ++j) {
        if ((R[size_t(3 * i + j)] * rep[size_t(j)]) % rep[size_t(i)]) good = false;
        double v = 0;
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l) v += R[size_t(3 * k + i)] * G[k][l] * R[size_t(3 * l + j)];
        if (std::fabs(v - G[i][j]) > 0.01 * gscale) good = false;
      }
    return allowed[R] = good;
  };
  std::map<int, std::vector<size_t>> by;
  for (size_t i = 0; i < n; ++i) by[F.el[i]].push_back(i);
  const auto& rare = std::min_element(by.begin(), by.end(), [](const auto& x, const auto& y) { return x.second.size() < y.second.size(); })->second;
  const size_t anchor = rare[0];
  std::map<IM, std::vector<Vec3>> W;
  std::function<const std::vector<Vec3>&(const IM&)> translations = [&](const IM& R) -> const std::vector<Vec3>& {
    if (auto it = W.find(R); it != W.end()) return it->second;
    std::vector<Vec3> w;
    const Vec3 r0 = rotate(R, F.f[anchor]);
    if (R == kIdentity) {   // the pure translations, every one tested
      for (size_t j : rare) {
        const Vec3 t = wrap01(F.f[j] - r0);
        if (std::any_of(w.begin(), w.end(), [&](const Vec3& u) { return frac_distance(cell, t - u) < tol; })) continue;
        if (maps_onto(F, R, t, tol)) w.push_back(t);
      }
    } else {   // a coset of the pure translations: the first that fits, plus each of them
      const std::vector<Vec3> T = translations(kIdentity);
      for (size_t j : rare) {
        const Vec3 t = wrap01(F.f[j] - r0);
        if (!maps_onto(F, R, t, tol)) continue;
        for (const auto& p : T) w.push_back(wrap01(t + p));
        break;
      }
    }
    std::stable_sort(w.begin(), w.end(), [&](const Vec3& x, const Vec3& y) { return frac_distance(cell, x) < frac_distance(cell, y); });
    return W[R] = std::move(w);
  };
  auto member = [&](const IM& R, const Vec3& t) {
    for (const auto& w : translations(R))
      if (frac_distance(cell, t - w) < tol) return true;
    return maps_onto(F, R, t, tol);
  };
  Hit best;
  size_t order = 0;
  for (const auto& S : search_order()) {
    if (order && S.ops.size() < order) break;
    if (!std::all_of(S.rot.begin(), S.rot.end(), ok)) continue;
    bool possible = true;
    for (size_t k = 0; k < S.ops.size() && possible; ++k) possible = !translations(S.rot[k]).empty();
    // the centring translations do not depend on the origin
    for (size_t k = 0; k < S.ops.size() && possible; ++k)
      if (S.rot[k] == kIdentity) possible = member(kIdentity, Vec3{S.ops[k].t[0], S.ops[k].t[1], S.ops[k].t[2]});
    if (!possible) continue;
    std::vector<IM> gr;
    std::vector<std::vector<Vec3>> choice;
    size_t combos = 1;
    for (size_t g : S.gens) {
      gr.push_back(S.rot[g]);
      choice.push_back(translations(S.rot[g]));
      combos *= choice.back().size();
    }
    if (combos > kChoiceBudget) {
      // many pure translations (a cell repeated other than along its axes): the shortest translation of each coset
      const auto& T = translations(kIdentity);
      combos = 1;
      for (auto& list : choice) {
        std::vector<Vec3> reps;
        for (const auto& w : list) {
          const bool same = std::any_of(reps.begin(), reps.end(), [&](const Vec3& u) {
            return std::any_of(T.begin(), T.end(), [&](const Vec3& p) { return frac_distance(cell, w - u - p) < tol; });
          });
          if (!same) reps.push_back(w);
        }
        list = std::move(reps);
        combos *= list.size();
      }
      if (combos > kChoiceBudget) continue;
    }
    std::set<std::array<long, 3>> tried;
    size_t budget = kOriginBudget;
    auto attempt = [&](const Vec3& s) {   // true when this setting can do no better
      std::array<long, 3> key;
      for (int i = 0; i < 3; ++i) key[size_t(i)] = std::lround(s[size_t(i)] * 4096) % 4096;
      if (!tried.insert(key).second || !budget) return false;
      --budget;
      for (size_t k = 0; k < S.ops.size(); ++k) {
        if (S.rot[k] == kIdentity) continue;
        const Vec3 rs = rotate(S.rot[k], s);
        if (!member(S.rot[k], Vec3{S.ops[k].t[0] + s[0] - rs[0], S.ops[k].t[1] + s[1] - rs[1], S.ops[k].t[2] + s[2] - rs[2]})) return false;
      }
      const double sh = frac_distance(cell, s);
      if (sh < best.shift - 1e-6) best = {&S, s, sh};
      return sh < 1e-6;
    };
    if (!attempt({0, 0, 0})) {
      bool done = false;
      std::vector<Vec3> d(choice.size());
      for (size_t c = 0; c < combos && !done && budget; ++c) {
        size_t rest = c;
        for (size_t k = 0; k < choice.size(); ++k) {
          const Vec3& w = choice[k][rest % choice[k].size()];
          rest /= choice[k].size();
          const SymOp& op = S.ops[S.gens[k]];
          d[k] = {w[0] - op.t[0], w[1] - op.t[1], w[2] - op.t[2]};
        }
        for (const auto& s : origin_solutions(gr, d))
          if ((done = attempt(s))) break;
      }
    }
    if (best.set && !order) order = S.ops.size();
  }
  return best;
}

}  // namespace

SymmetryFound find_symmetry(const System& s, double tol) {
  if (!s.cell.valid()) throw CrystalError("finding symmetry needs a cell");
  const size_t n = s.atoms.size();
  if (n == 0) throw CrystalError("no atoms");
  if (n > 2000) throw CrystalError("more than 2000 atoms: find symmetry on the unit cell");
  std::vector<Vec3> f(n);
  std::vector<int> el(n);
  for (size_t i = 0; i < n; ++i) f[i] = wrap01(s.cell.to_fractional(s.atoms[i].pos)), el[i] = s.atoms[i].element;
  const Frame whole(s.cell, f, el, tol);
  // pure translations along the axes: the input repeats a smaller cell, and the group is searched in that cell
  std::array<int, 3> rep{1, 1, 1};
  {
    std::map<int, std::vector<size_t>> by;
    for (size_t i = 0; i < n; ++i) by[el[i]].push_back(i);
    const auto& rare = std::min_element(by.begin(), by.end(), [](const auto& x, const auto& y) { return x.second.size() < y.second.size(); })->second;
    std::vector<Vec3> pure;
    for (size_t j : rare) {
      const Vec3 t = wrap01(f[j] - f[rare[0]]);
      if (frac_distance(s.cell, t) >= tol && maps_onto(whole, kIdentity, t, tol)) pure.push_back(t);
    }
    const Vec3 axes[3] = {s.cell.a, s.cell.b, s.cell.c};
    for (int d = 0; d < 3 && !pure.empty(); ++d) {
      double step = 1;
      for (const auto& t : pure) {
        Vec3 off = t;
        off[size_t(d)] = 0;
        if (frac_distance(s.cell, off) >= tol) continue;   // not along this axis
        const double x = std::min(t[size_t(d)], 1 - t[size_t(d)]);
        if (x * norm(axes[d]) >= tol) step = std::min(step, x);
      }
      const int k = int(std::lround(1 / step));
      Vec3 t{0, 0, 0};
      t[size_t(d)] = 1.0 / std::max(k, 1);
      if (k >= 2 && maps_onto(whole, kIdentity, t, tol)) rep[size_t(d)] = k;
    }
  }
  const size_t copies = size_t(rep[0]) * size_t(rep[1]) * size_t(rep[2]);
  // the group of frame F (rep copies of it make the input), as the result in the input's coordinates
  bool too_many = false;
  auto result = [&](const Frame& F, const std::array<int, 3>& k, const std::vector<size_t>& atom_of) {
    const Hit h = search(F, tol, k);
    SymmetryFound r;
    r.cell = F.cell;
    r.supercell = k;
    std::vector<SymOp> ops;   // in F's coordinates: t′ = t + (I − R) s
    if (h.set) {
      r.key = h.set->sg->key, r.hm = h.set->sg->hm, r.number = h.set->sg->number, r.system = crystal_system(r.number);
      r.operations = int(h.set->ops.size());
      r.origin = h.origin;
      for (size_t q = 0; q < h.set->ops.size(); ++q) {
        SymOp o = h.set->ops[q];
        const Vec3 rs = rotate(h.set->rot[q], h.origin);
        for (int i = 0; i < 3; ++i) o.t[i] += h.origin[size_t(i)] - rs[size_t(i)];
        ops.push_back(o);
      }
    } else {
      r.key = "1", r.hm = "P 1", r.number = 1, r.system = "triclinic", r.operations = 1;
      ops.push_back(SymOp{});
    }
    // asymmetric unit: one atom per orbit, in the setting's origin
    std::vector<char> done(F.f.size(), 0);
    for (size_t i = 0; i < F.f.size(); ++i) {
      if (done[i]) continue;
      for (const auto& o : ops) {
        double d = 1e300;
        const long j = F.nearest(apply_op(o, F.f[i]), F.el[i], &d);
        if (j >= 0 && d < tol) done[size_t(j)] = 1;
      }
      done[i] = 1;
      const Atom& a = s.atoms[atom_of[i]];
      r.sites.push_back({a.name.empty() ? element(a.element).symbol : a.name, a.element, wrap01(F.f[i] - r.origin), 1});
    }
    // every operation in the input's coordinates, f = D⁻¹ f_F: R′ = D⁻¹ R D, t′ = D⁻¹ t, with each translation of F's cell
    const size_t N = size_t(k[0]) * size_t(k[1]) * size_t(k[2]);
    too_many = ops.size() * N > kMaxOps;
    if (too_many) return r;
    for (const auto& o : ops) {
      SymOp b;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) b.R[i][j] = o.R[i][j] * k[size_t(j)] / k[size_t(i)];
        b.t[i] = o.t[i] / k[size_t(i)];
      }
      for (int x = 0; x < k[0]; ++x)
        for (int y = 0; y < k[1]; ++y)
          for (int z = 0; z < k[2]; ++z) {
            SymOp c = b;
            c.t[0] += double(x) / k[0], c.t[1] += double(y) / k[1], c.t[2] += double(z) / k[2];
            for (double& t : c.t) t -= std::floor(t + 1e-9);
            r.ops.push_back(c);
          }
    }
    return r;
  };
  if (copies > 1) {
    Cell c = s.cell;
    c.a = s.cell.a * (1.0 / rep[0]), c.b = s.cell.b * (1.0 / rep[1]), c.c = s.cell.c * (1.0 / rep[2]);
    std::vector<Vec3> g;
    std::vector<int> ge;
    std::vector<size_t> atom_of;
    for (size_t i = 0; i < n; ++i) {
      const Vec3 h = wrap01(Vec3{f[i][0] * rep[0], f[i][1] * rep[1], f[i][2] * rep[2]});
      bool dup = false;
      for (size_t q = 0; q < g.size() && !dup; ++q) dup = ge[q] == el[i] && frac_distance(c, h - g[q]) < tol;
      if (dup) continue;
      g.push_back(h), ge.push_back(el[i]), atom_of.push_back(i);
    }
    if (g.size() * copies == n) {
      const Frame small(c, g, ge, tol);
      SymmetryFound r = result(small, rep, atom_of);
      if (!too_many) return r;
      // too many operations to list in the supercell's axes: its own group there instead
    }
  }
  std::vector<size_t> self(n);
  for (size_t i = 0; i < n; ++i) self[i] = i;
  return result(whole, {1, 1, 1}, self);
}

// ---------------------------------------------------------------- the space group kept with the structure

Symmetry symmetry_of(const System& s, double tol) {
  const SymmetryFound f = find_symmetry(s, tol);
  Symmetry g;
  g.name = f.hm, g.number = f.number;
  g.ops = f.ops;
  return g;
}

bool symmetry_holds(const System& s, const Symmetry& g, double tol, double* worst) {
  if (worst) *worst = 0;
  if (!s.cell.valid() || g.ops.empty()) return false;
  const size_t n = s.atoms.size();
  std::vector<Vec3> f(n);
  std::vector<int> el(n);
  for (size_t i = 0; i < n; ++i) f[i] = s.cell.to_fractional(s.atoms[i].pos), el[i] = s.atoms[i].element;
  std::map<int, std::vector<size_t>> by;
  for (size_t i = 0; i < n; ++i) by[el[i]].push_back(i);
  const Frame F(s.cell, f, el, std::max(tol, 1e-6));
  double w = 0;
  for (const auto& op : g.ops)
    for (size_t i = 0; i < n; ++i) {
      const Vec3 h = apply_op(op, f[i]);
      double best = 1e300;
      F.nearest(h, el[i], &best, worst ? -1 : tol);   // the exact nearest only when the worst image is wanted
      if (best >= tol) {   // nothing within tol in the neighbouring bins: the exact distance from every atom of the element
        if (!worst) return false;
        for (size_t j : by[el[i]]) best = std::min(best, frac_distance(s.cell, h - f[j]));
      }
      w = std::max(w, best);
      if (w > tol && !worst) return false;
    }
  if (worst) *worst = w;
  return w <= tol;
}

std::vector<size_t> asymmetric_unit(const System& s, const Symmetry& g, double tol) {
  const size_t n = s.atoms.size();
  std::vector<Vec3> f(n);
  for (size_t i = 0; i < n; ++i) f[i] = s.cell.to_fractional(s.atoms[i].pos);
  std::vector<char> done(n, 0);
  std::vector<size_t> keep;
  for (size_t i = 0; i < n; ++i) {
    if (done[i]) continue;
    keep.push_back(i);
    for (const auto& op : g.ops) {
      const Vec3 h = apply_op(op, f[i]);
      for (size_t j = 0; j < n; ++j)
        if (!done[j] && s.atoms[j].element == s.atoms[i].element && frac_distance(s.cell, h - f[j]) < tol) done[j] = 1;
    }
    done[i] = 1;
  }
  return keep;
}

Cell symmetric_cell(const Cell& c, const Symmetry& g) {
  const Vec3 e[3] = {c.a, c.b, c.c};
  double G[3][3], A[3][3] = {};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) G[i][j] = dot(e[i], e[j]);
  // the distinct rotation parts (the point group); Rᵀ G R averaged over them
  std::vector<const SymOp*> rots;
  for (const auto& op : g.ops) {
    bool seen = false;
    for (const auto* r : rots) {
      bool same = true;
      for (int i = 0; i < 3 && same; ++i)
        for (int j = 0; j < 3 && same; ++j) same = std::fabs(r->R[i][j] - op.R[i][j]) < 1e-9;
      seen = seen || same;
    }
    if (!seen) rots.push_back(&op);
  }
  if (rots.empty()) return c;
  for (const auto* r : rots)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double v = 0;
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l) v += r->R[k][i] * G[k][l] * r->R[l][j];
        A[i][j] += v / double(rots.size());
      }
  const double a = std::sqrt(A[0][0]), b = std::sqrt(A[1][1]), cc = std::sqrt(A[2][2]);
  auto ang = [](double x) { return std::acos(std::clamp(x, -1.0, 1.0)) * 180 / M_PI; };
  Cell out = cell_from(a, b, cc, ang(A[1][2] / (b * cc)), ang(A[0][2] / (a * cc)), ang(A[0][1] / (a * b)));
  out.origin = c.origin;
  out.periodic = c.periodic;
  return out;
}

void set_cell_symmetric(System& s, double a, double b, double c, double alpha, double beta, double gamma) {
  Cell want = cell_from(a, b, c, alpha, beta, gamma);
  if (!(want.volume() > 1e-9)) throw CrystalError("these cell parameters give no volume");
  if (!s.symmetry.ops.empty()) want = symmetric_cell(want, s.symmetry);
  want.origin = s.cell.origin;
  want.periodic = s.cell.periodic;
  const Cell old = s.cell;
  for (auto& at : s.atoms) at.pos = want.to_cartesian(old.to_fractional(at.pos));
  s.cell = want;
}

size_t move_atom_symmetric(System& s, size_t i, const Vec3& pos, double tol) {
  if (i >= s.atoms.size()) throw std::out_of_range("move: no atom " + std::to_string(i + 1));
  const auto& ops = s.symmetry.ops;
  if (ops.empty() || !s.cell.valid() || !symmetry_holds(s, s.symmetry, tol)) {
    s.atoms[i].pos = pos;
    return 1;
  }
  const size_t n = s.atoms.size();
  std::vector<Vec3> f(n);
  for (size_t k = 0; k < n; ++k) f[k] = s.cell.to_fractional(s.atoms[k].pos);
  const Vec3 fi = f[i], d = s.cell.to_fractional(pos) - fi;
  // the site-symmetry group: the move projected onto what it leaves invariant
  Vec3 dp{0, 0, 0};
  int ns = 0;
  for (const auto& op : ops)
    if (frac_distance(s.cell, apply_op(op, fi) - fi) < tol) {
      SymOp r = op;
      for (double& t : r.t) t = 0;
      dp = dp + apply_op(r, d);
      ++ns;
    }
  dp = dp * (1.0 / std::max(1, ns));
  std::vector<char> done(n, 0);
  size_t moved = 0;
  for (const auto& op : ops) {
    const Vec3 h = apply_op(op, fi);
    for (size_t j = 0; j < n; ++j) {
      if (done[j] || s.atoms[j].element != s.atoms[i].element || frac_distance(s.cell, h - f[j]) >= tol) continue;
      Vec3 shift = f[j] - h;
      for (auto& x : shift) x = std::round(x);
      s.atoms[j].pos = s.cell.to_cartesian(apply_op(op, fi + dp) + shift);
      done[j] = 1, ++moved;
      break;
    }
  }
  return moved;
}

}  // namespace caps
