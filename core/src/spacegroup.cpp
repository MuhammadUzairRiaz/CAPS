// CAPS space groups (see caps/spacegroup.hpp).
#include "caps/spacegroup.hpp"

#include <cstdlib>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
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
      int best = 1;
      for (int d : {2, 3, 4, 6, 8, 12, 24}) if (std::fabs(t * d - std::round(t * d)) < 1e-6) { best = d; break; }
      c += "+" + std::to_string(int(std::lround(t * best))) + "/" + std::to_string(best);
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
  s.has_mol = true;
  s.unwrapped = true;
  R.atoms_per_cell = s.atoms.size();
  R.volume = s.cell.volume();
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
  out.notes.insert(out.notes.begin(), "primitive cell of the " + std::string(1, char(std::toupper(static_cast<unsigned char>(centring)))) + "-centred lattice · " +
                                          std::to_string(out.atoms.size()) + " atoms (from " + std::to_string(s.atoms.size()) + ")");
  return out;
}

SymmetryFound find_symmetry(const System& s, double tol) {
  if (!s.cell.valid()) throw CrystalError("finding symmetry needs a cell");
  const size_t n = s.atoms.size();
  if (n == 0) throw CrystalError("no atoms");
  if (n > 2000) throw CrystalError("more than 2000 atoms: find symmetry on the unit cell");
  std::vector<Vec3> f(n);
  for (size_t i = 0; i < n; ++i) f[i] = wrap01(s.cell.to_fractional(s.atoms[i].pos));
  // metric tensor: an operation belongs only if RᵀGR = G
  double G[3][3];
  const Vec3 e[3] = {s.cell.a, s.cell.b, s.cell.c};
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) G[i][j] = dot(e[i], e[j]);
  const double gscale = std::max({G[0][0], G[1][1], G[2][2]});
  auto metric_ok = [&](const SymOp& op) {
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double v = 0;
        for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) v += op.R[k][i] * G[k][l] * op.R[l][j];
        if (std::fabs(v - G[i][j]) > 0.01 * gscale) return false;
      }
    return true;
  };
  // the rarest element anchors the origin candidates
  std::map<int, std::vector<size_t>> by;
  for (size_t i = 0; i < n; ++i) by[s.atoms[i].element].push_back(i);
  const auto& rare = std::min_element(by.begin(), by.end(), [](const auto& x, const auto& y) { return x.second.size() < y.second.size(); })->second;
  std::vector<Vec3> shifts = {{0, 0, 0}};
  for (size_t j : rare)
    for (int h = 0; h < 8; ++h) {
      Vec3 m = (f[rare[0]] + f[j]) * 0.5;
      m[0] += (h & 1) * 0.5; m[1] += (h >> 1 & 1) * 0.5; m[2] += (h >> 2 & 1) * 0.5;
      shifts.push_back(wrap01(m));
    }
  for (size_t j : rare) shifts.push_back(f[j]);
  auto fits = [&](const SymOp& op, const Vec3& sh) {
    // the operation in the structure's coordinates: t' = t + (I − R) s
    SymOp o = op;
    for (int i = 0; i < 3; ++i) { double v = 0; for (int k = 0; k < 3; ++k) v += ((i == k) - op.R[i][k]) * sh[k]; o.t[i] += v; }
    for (size_t i = 0; i < n; ++i) {
      const Vec3 g = apply_op(o, f[i]);
      bool hit = false;
      for (size_t j : by.at(s.atoms[i].element))
        if (frac_distance(s.cell, g - f[j]) < tol) { hit = true; break; }
      if (!hit) return false;
    }
    return true;
  };
  // settings from the highest order down; the first that fits is the symmetry
  std::vector<const SpaceGroupSetting*> order;
  std::set<std::string> seen;
  for (const auto& sg : space_group_settings()) if (seen.insert(sg.hall).second) order.push_back(&sg);
  std::stable_sort(order.begin(), order.end(), [](const SpaceGroupSetting* x, const SpaceGroupSetting* y) {
    return hall_operations(x->hall).size() > hall_operations(y->hall).size();
  });
  for (const auto* sg : order) {
    const auto ops = hall_operations(sg->hall);
    if (!std::all_of(ops.begin(), ops.end(), metric_ok)) continue;
    for (const auto& sh : shifts) {
      bool ok = true;
      for (size_t k = 1; k < ops.size() && ok; ++k) ok = fits(ops[k], sh);
      if (!ok) continue;
      SymmetryFound r;
      r.key = sg->key; r.hm = sg->hm; r.number = sg->number; r.system = crystal_system(sg->number);
      r.operations = int(ops.size());
      r.origin = sh;
      // asymmetric unit: one atom per orbit, in the setting's origin
      std::vector<char> done(n, 0);
      for (size_t i = 0; i < n; ++i) {
        if (done[i]) continue;
        const Vec3 fi = wrap01(f[i] - sh);
        for (const auto& op : ops) {
          const Vec3 g = wrap01(apply_op(op, fi) + sh);
          for (size_t j : by.at(s.atoms[i].element)) if (frac_distance(s.cell, g - f[j]) < tol) done[j] = 1;
        }
        r.sites.push_back({s.atoms[i].name.empty() ? element(s.atoms[i].element).symbol : s.atoms[i].name, s.atoms[i].element, fi, 1});
      }
      return r;
    }
  }
  SymmetryFound r;
  r.key = "1"; r.hm = "P 1"; r.number = 1; r.system = "triclinic"; r.operations = 1;
  for (size_t i = 0; i < n; ++i) r.sites.push_back({s.atoms[i].name, s.atoms[i].element, f[i], 1});
  return r;
}

}  // namespace caps
