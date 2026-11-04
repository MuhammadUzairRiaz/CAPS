// CAPS restraints (see restraints.hpp).
#include "caps/restraints.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>
#include <stdexcept>

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180;
constexpr double kKJ = 4.184;   // kJ per kcal

Vec3 at(const std::vector<double>& x, uint32_t i) { return {x[3 * i], x[3 * i + 1], x[3 * i + 2]}; }
Vec3 mi(const Cell& c, Vec3 v) { return c.valid() ? c.minimum_image(v) : v; }

// signed distance outside [lo, hi] (0 inside)
double outside(double q, double lo, double hi) { return q < lo ? q - lo : q > hi ? q - hi : 0.0; }

double ramped(const Restraint& r, double ramp) { return r.k_end >= 0 ? r.k + (r.k_end - r.k) * std::clamp(ramp, 0.0, 1.0) : r.k; }

std::string g(double v) {
  char b[40];
  std::snprintf(b, sizeof b, "%.6f", v);
  std::string s = b;
  while (s.size() > 1 && s.back() == '0') s.pop_back();
  if (!s.empty() && s.back() == '.') s.pop_back();
  return s == "-0" ? "0" : s;
}
}  // namespace

int restraint_atom_count(Restraint::Kind k) {
  switch (k) {
    case Restraint::Kind::Distance: return 2;
    case Restraint::Kind::Angle: return 3;
    case Restraint::Kind::Torsion: return 4;
    case Restraint::Kind::Position: return 1;
  }
  return 0;
}

const char* to_string(Restraint::Kind k) {
  switch (k) {
    case Restraint::Kind::Distance: return "distance";
    case Restraint::Kind::Angle: return "angle";
    case Restraint::Kind::Torsion: return "torsion";
    case Restraint::Kind::Position: return "position";
  }
  return "?";
}

namespace {
struct Geo {
  double q = 0;                       // value (Å or rad)
  std::array<Vec3, 4> grad{};         // dq/dx_a
  std::array<Vec3, 4> rel{};          // positions relative to the first atom (minimum-image chain)
};

bool geometry(const Restraint& r, const std::vector<double>& x, const Cell& c, Geo& out) {
  const auto& a = r.atoms;
  switch (r.kind) {
    case Restraint::Kind::Distance: {
      const Vec3 d = mi(c, at(x, a[1]) - at(x, a[0]));
      const double len = norm(d);
      if (len < 1e-12) return false;
      out.q = len;
      out.grad[1] = d * (1 / len), out.grad[0] = d * (-1 / len);
      out.rel[1] = d;
      return true;
    }
    case Restraint::Kind::Angle: {
      const Vec3 u = mi(c, at(x, a[0]) - at(x, a[1])), v = mi(c, at(x, a[2]) - at(x, a[1]));
      const double lu = norm(u), lv = norm(v);
      if (lu < 1e-12 || lv < 1e-12) return false;
      const double cs = std::clamp(dot(u, v) / (lu * lv), -1.0, 1.0);
      const double sn = std::sqrt(std::max(1e-16, 1 - cs * cs));
      out.q = std::acos(cs);
      const Vec3 uh = u * (1 / lu), vh = v * (1 / lv);
      out.grad[0] = (vh - uh * cs) * (-1 / (sn * lu));
      out.grad[2] = (uh - vh * cs) * (-1 / (sn * lv));
      out.grad[1] = (out.grad[0] + out.grad[2]) * -1.0;
      out.rel[1] = Vec3{0, 0, 0} - u, out.rel[2] = v - u;   // relative to atom i
      return true;
    }
    case Restraint::Kind::Torsion: {   // IUPAC sign; gradient of Blondel & Karplus, J. Comput. Chem. 17, 1132 (1996)
      const Vec3 b1 = mi(c, at(x, a[1]) - at(x, a[0])), b2 = mi(c, at(x, a[2]) - at(x, a[1])), b3 = mi(c, at(x, a[3]) - at(x, a[2]));
      const Vec3 m = cross(b1, b2), nn = cross(b2, b3);
      const double m2 = dot(m, m), n2 = dot(nn, nn), lb2 = norm(b2);
      if (m2 < 1e-12 || n2 < 1e-12 || lb2 < 1e-9) return false;
      out.q = std::atan2(lb2 * dot(b1, nn), dot(m, nn));
      const Vec3 gi = m * (-lb2 / m2), gl = nn * (lb2 / n2);
      const double s1 = dot(b1, b2) / (lb2 * lb2), s3 = dot(b3, b2) / (lb2 * lb2);
      out.grad[0] = gi, out.grad[3] = gl;
      out.grad[1] = gi * (s1 - 1) - gl * s3;
      out.grad[2] = gl * (s3 - 1) - gi * s1;
      out.rel[1] = b1, out.rel[2] = b1 + b2, out.rel[3] = b1 + b2 + b3;
      return true;
    }
    case Restraint::Kind::Position: {
      Vec3 d = mi(c, at(x, a[0]) - r.anchor);
      for (int k = 0; k < 3; ++k)
        if (!((r.axes >> k) & 1)) d[size_t(k)] = 0;
      const double len = norm(d);
      out.q = len;
      out.grad[0] = len > 1e-12 ? d * (1 / len) : Vec3{0, 0, 0};
      return true;
    }
  }
  return false;
}
}  // namespace

double restraint_value(const Restraint& r, const std::vector<double>& x, const Cell& cell) {
  Geo gm;
  if (!geometry(r, x, cell, gm)) return std::nan("");
  return r.kind == Restraint::Kind::Angle || r.kind == Restraint::Kind::Torsion ? gm.q / kDeg : gm.q;
}

RestraintTerms add_restraint_forces(const std::vector<Restraint>& rs, const std::vector<double>& x, const Cell& cell, std::vector<double>& f, double ramp) {
  RestraintTerms t;
  const size_t n = x.size() / 3;
  for (const auto& r : rs) {
    const int na = restraint_atom_count(r.kind);
    bool ok = true;
    for (int q = 0; q < na; ++q) ok = ok && r.atoms[size_t(q)] < n;
    if (!ok) continue;
    const double k = ramped(r, ramp);
    if (r.kind == Restraint::Kind::Position) {
      // harmonic: ½ k |d|² over the axes; flat-bottomed sphere: ½ k (|d| − R)² beyond R
      Vec3 d = mi(cell, at(x, r.atoms[0]) - r.anchor);
      for (int a = 0; a < 3; ++a)
        if (!((r.axes >> a) & 1)) d[size_t(a)] = 0;
      const double len = norm(d);
      const double R = std::max(0.0, r.upper);
      if (len <= R || len < 1e-300) continue;
      const double ex = len - R;
      t.energy += 0.5 * k * ex * ex;
      for (int a = 0; a < 3; ++a) f[3 * r.atoms[0] + size_t(a)] -= k * ex * d[size_t(a)] / len;
      continue;
    }
    Geo gm;
    if (!geometry(r, x, cell, gm)) continue;
    double dev;
    if (r.kind == Restraint::Kind::Distance) {
      dev = outside(gm.q, r.lower, r.upper);
    } else if (r.kind == Restraint::Kind::Angle) {
      dev = outside(gm.q, r.lower * kDeg, r.upper * kDeg);
    } else {   // on the circle: the deviation from the centre, less the half-width
      const double c = 0.5 * (r.lower + r.upper) * kDeg, hw = 0.5 * std::fabs(r.upper - r.lower) * kDeg;
      double d = gm.q - c;
      d -= 2 * kPi * std::round(d / (2 * kPi));
      dev = std::fabs(d) > hw ? (d > 0 ? d - hw : d + hw) : 0.0;
    }
    if (dev == 0) continue;
    t.energy += 0.5 * k * dev * dev;
    const double dE = k * dev;   // dE/dq
    for (int a = 0; a < na; ++a) {
      const Vec3 fa = gm.grad[size_t(a)] * (-dE);
      for (int c = 0; c < 3; ++c) f[3 * r.atoms[size_t(a)] + size_t(c)] += fa[size_t(c)];
      const Vec3& p = gm.rel[size_t(a)];
      t.virial += dot(p, fa);
      t.w[0] += p[0] * fa[0], t.w[1] += p[1] * fa[1], t.w[2] += p[2] * fa[2];
      t.w[3] += 0.5 * (p[0] * fa[1] + p[1] * fa[0]), t.w[4] += 0.5 * (p[0] * fa[2] + p[2] * fa[0]), t.w[5] += 0.5 * (p[1] * fa[2] + p[2] * fa[1]);
    }
  }
  return t;
}

namespace {
std::vector<uint32_t> parse_selection(const std::string& s, size_t n) {
  std::vector<uint32_t> out;
  std::stringstream ss(s);
  std::string part;
  while (std::getline(ss, part, ',')) {
    if (part.empty()) continue;
    const auto dash = part.find('-', 1);
    long a = std::stol(part.substr(0, dash)), b = dash == std::string::npos ? a : std::stol(part.substr(dash + 1));
    if (a < 1 || b < a || size_t(b) > n) throw std::invalid_argument("atoms '" + part + "' outside 1.." + std::to_string(n));
    for (long i = a; i <= b; ++i) out.push_back(uint32_t(i - 1));
  }
  return out;
}

// "R" or "LO:HI"
void parse_bounds(const std::string& v, double& lo, double& hi) {
  const auto c = v.find(':');
  if (c == std::string::npos) lo = hi = std::stod(v);
  else lo = std::stod(v.substr(0, c)), hi = std::stod(v.substr(c + 1));
  if (hi < lo) std::swap(lo, hi);
}
}  // namespace

std::vector<Restraint> parse_restraints(const std::string& text, const System& s) {
  std::vector<Restraint> out;
  std::string norm_text = text;
  for (auto& ch : norm_text) if (ch == '\n') ch = ';';
  std::stringstream items(norm_text);
  std::string item;
  const size_t n = s.atoms.size();
  while (std::getline(items, item, ';')) {
    std::stringstream ws(item);
    std::vector<std::string> w;
    for (std::string t; ws >> t;) w.push_back(t);
    if (w.empty() || w[0][0] == '#') continue;
    try {
      Restraint r;
      const std::string kind = w[0];
      size_t p = 1;
      auto atom = [&](const std::string& t) {
        const long i = std::stol(t);
        if (i < 1 || size_t(i) > n) throw std::invalid_argument("atom " + t + " outside 1.." + std::to_string(n));
        return uint32_t(i - 1);
      };
      // the options after the numbers: [K] [ramp KEND] [axes xyz] [radius R]
      auto tail = [&](double kdef) {
        r.k = kdef;
        bool kset = false;
        while (p < w.size()) {
          if (w[p] == "ramp" && p + 1 < w.size()) r.k_end = std::stod(w[p + 1]), p += 2;
          else if (w[p] == "axes" && p + 1 < w.size()) {
            r.axes = 0;
            for (char ch : w[p + 1]) r.axes |= ch == 'x' ? 1 : ch == 'y' ? 2 : ch == 'z' ? 4 : throw std::invalid_argument("axes: x, y, z");
            p += 2;
          } else if (w[p] == "radius" && p + 1 < w.size()) r.upper = std::stod(w[p + 1]), p += 2;
          else if (!kset) r.k = std::stod(w[p]), kset = true, ++p;
          else throw std::invalid_argument("'" + w[p] + "' not understood");
        }
        if (r.k < 0 || (r.k_end < 0 && r.k_end != -1)) throw std::invalid_argument("force constants are not negative");
      };
      if (kind == "distance" || kind == "angle" || kind == "torsion") {
        r.kind = kind == "distance" ? Restraint::Kind::Distance : kind == "angle" ? Restraint::Kind::Angle : Restraint::Kind::Torsion;
        const int na = restraint_atom_count(r.kind);
        if (w.size() < size_t(na) + 2) throw std::invalid_argument("needs " + std::to_string(na) + " atoms and a value");
        for (int q = 0; q < na; ++q) r.atoms[size_t(q)] = atom(w[p++]);
        parse_bounds(w[p++], r.lower, r.upper);
        tail(r.kind == Restraint::Kind::Distance ? 10.0 : 50.0);
        if (r.kind == Restraint::Kind::Distance && r.lower < 0) throw std::invalid_argument("a distance is not negative");
        if (r.kind == Restraint::Kind::Angle && (r.lower < 0 || r.upper > 180)) throw std::invalid_argument("an angle lies in 0..180°");
        if (r.kind == Restraint::Kind::Torsion && r.upper - r.lower >= 360) throw std::invalid_argument("a torsion's bounds span less than 360°");
        out.push_back(r);
      } else if (kind == "position") {
        if (w.size() < 2) throw std::invalid_argument("needs atoms");
        const auto sel = parse_selection(w[p++], n);
        r.kind = Restraint::Kind::Position;
        tail(10.0);
        if (r.axes == 0) throw std::invalid_argument("axes: at least one");
        for (uint32_t i : sel) {
          Restraint q = r;
          q.atoms[0] = i;
          q.anchor = s.atoms[i].pos;
          out.push_back(q);
        }
      } else {
        throw std::invalid_argument("'" + kind + "': distance, angle, torsion or position");
      }
    } catch (const std::invalid_argument& e) {
      throw std::invalid_argument("restraint '" + item + "': " + e.what());
    } catch (const std::out_of_range&) {
      throw std::invalid_argument("restraint '" + item + "': a number out of range");
    }
  }
  return out;
}

std::string restraints_text(const std::vector<Restraint>& rs) {
  std::ostringstream o;
  for (size_t idx = 0; idx < rs.size(); ++idx) {
    const auto& r = rs[idx];
    o << to_string(r.kind);
    if (r.kind == Restraint::Kind::Position) {   // runs of atoms with the same settings as ranges: 3-10,14
      auto same = [&](const Restraint& q) { return q.kind == r.kind && q.k == r.k && q.k_end == r.k_end && q.axes == r.axes && q.upper == r.upper; };
      std::string sel;
      size_t j = idx;
      while (j < rs.size() && same(rs[j])) {
        size_t e = j;
        while (e + 1 < rs.size() && same(rs[e + 1]) && rs[e + 1].atoms[0] == rs[e].atoms[0] + 1) ++e;
        sel += (sel.empty() ? "" : ",") + std::to_string(rs[j].atoms[0] + 1) + (e > j ? "-" + std::to_string(rs[e].atoms[0] + 1) : "");
        j = e + 1;
      }
      o << " " << sel;
      idx = j - 1;
    } else {
      for (int q = 0; q < restraint_atom_count(r.kind); ++q) o << " " << r.atoms[size_t(q)] + 1;
    }
    if (r.kind != Restraint::Kind::Position) o << " " << (r.lower == r.upper ? g(r.lower) : g(r.lower) + ":" + g(r.upper));
    o << " " << g(r.k);
    if (r.kind == Restraint::Kind::Position) {
      if (r.axes != 7) o << " axes " << ((r.axes & 1) ? "x" : "") << ((r.axes & 2) ? "y" : "") << ((r.axes & 4) ? "z" : "");
      if (r.upper > 0) o << " radius " << g(r.upper);
    }
    if (r.k_end >= 0) o << " ramp " << g(r.k_end);
    o << "\n";
  }
  return o.str();
}

Json restraints_json(const std::vector<Restraint>& rs, const std::vector<double>* x, const Cell* cell) {
  Json a = Json::array();
  for (const auto& r : rs) {
    Json j = Json::object(), at = Json::array();
    j["kind"] = std::string(to_string(r.kind));
    for (int q = 0; q < restraint_atom_count(r.kind); ++q) at.push_back(double(r.atoms[size_t(q)]));
    j["atoms"] = at, j["lower"] = r.lower, j["upper"] = r.upper, j["k"] = r.k, j["k_end"] = r.k_end, j["axes"] = double(r.axes);
    if (r.kind == Restraint::Kind::Position) {
      Json an = Json::array();
      for (double v : r.anchor) an.push_back(v);
      j["anchor"] = an;
    }
    if (x && cell) {
      const double v = restraint_value(r, *x, *cell);
      if (std::isfinite(v)) j["value"] = v;
      std::vector<double> f(x->size(), 0.0);
      j["energy"] = add_restraint_forces({r}, *x, *cell, f).energy;
    }
    a.push_back(j);
  }
  return a;
}

std::string lammps_restraint_lines(const std::vector<Restraint>& rs, const System& s, double energy_scale, std::vector<std::string>* notes) {
  std::ostringstream o;
  std::vector<std::string> parts;
  int skipped_flat = 0, skipped_anchor = 0, skipped_fbpos = 0;
  for (const auto& r : rs) {
    const double ks = r.k * energy_scale, ke = (r.k_end >= 0 ? r.k_end : r.k) * energy_scale;
    const auto& a = r.atoms;
    switch (r.kind) {
      case Restraint::Kind::Distance:
        if (r.lower == r.upper)
          parts.push_back("bond " + std::to_string(a[0] + 1) + " " + std::to_string(a[1] + 1) + " " + g(ks / 2) + " " + g(ke / 2) + " " + g(r.lower) + " " + g(r.lower));
        else if (r.upper >= 1e6)   // a lower bound alone
          parts.push_back("lbound " + std::to_string(a[0] + 1) + " " + std::to_string(a[1] + 1) + " " + g(ks / 2) + " " + g(ke / 2) + " " + g(r.lower) + " " + g(r.lower));
        else ++skipped_flat;
        break;
      case Restraint::Kind::Angle:
        if (r.lower == r.upper)
          parts.push_back("angle " + std::to_string(a[0] + 1) + " " + std::to_string(a[1] + 1) + " " + std::to_string(a[2] + 1) + " " + g(ks / 2) + " " + g(ke / 2) + " " + g(r.lower));
        else ++skipped_flat;
        break;
      case Restraint::Kind::Torsion:
        if (r.lower == r.upper) {
          double p0 = r.lower + 180;   // LAMMPS's φ is 180° from the IUPAC sign
          p0 -= 360 * std::round(p0 / 360);
          parts.push_back("dihedral " + std::to_string(a[0] + 1) + " " + std::to_string(a[1] + 1) + " " + std::to_string(a[2] + 1) + " " + std::to_string(a[3] + 1) + " " +
                          g(ks) + " " + g(ke) + " " + g(p0));
        } else ++skipped_flat;
        break;
      case Restraint::Kind::Position: break;
    }
  }
  if (!parts.empty()) {
    o << "\n# restraints (CAPS ½ k Δ²: bond and angle K = k/2; dihedral K [1 − cos(φ − φ0)] ≈ ½ K Δφ², φ0 written + 180° for LAMMPS's reckoning)\n";
    o << "fix             restraints all restrain";
    for (const auto& p : parts) o << " &\n                " << p;
    o << "\nfix_modify      restraints energy yes\n";
  }
  // positions: one group and spring/self per (k, k_end, axes); the anchors must be where the data file puts the atoms
  std::map<std::tuple<double, double, int>, std::vector<uint32_t>> springs;
  for (const auto& r : rs) {
    if (r.kind != Restraint::Kind::Position) continue;
    if (r.upper > 0) { ++skipped_fbpos; continue; }
    const uint32_t i = r.atoms[0];
    if (i >= s.atoms.size() || norm(mi(s.cell, s.atoms[i].pos - r.anchor)) > 1e-6) { ++skipped_anchor; continue; }
    springs[{r.k, r.k_end, r.axes}].push_back(i);
  }
  int id = 0;
  bool ramped = false;
  for (const auto& [key, atoms] : springs) {
    const auto [k, kend, axes] = key;
    ++id;
    const std::string grp = "posres" + std::to_string(id);
    o << (id == 1 ? "\n# position restraints: ½ K r² from where each atom starts (fix spring/self)\n" : "");
    o << "group           " << grp << " id";
    for (size_t q = 0; q < atoms.size(); ++q) o << (q % 16 == 0 && q ? " &\n               " : "") << " " << atoms[q] + 1;
    o << "\n";
    std::string dir = std::string((axes & 1) ? "x" : "") + ((axes & 2) ? "y" : "") + ((axes & 4) ? "z" : "");
    o << "fix             " << grp << " " << grp << " spring/self " << g(k * energy_scale) << (dir == "xyz" ? std::string() : " " + dir) << "\n"
      << "fix_modify      " << grp << " energy yes\n";
    if (kend >= 0) ramped = true;
  }
  if (notes) {
    if (skipped_flat) notes->push_back(std::to_string(skipped_flat) + " flat-bottomed restraint(s) left out: LAMMPS's fix restrain has only harmonic bonds, angles and dihedrals and a distance's lower bound");
    if (skipped_fbpos) notes->push_back(std::to_string(skipped_fbpos) + " flat-bottomed position restraint(s) left out: fix spring/self is harmonic");
    if (skipped_anchor) notes->push_back(std::to_string(skipped_anchor) + " position restraint(s) left out: their anchors are not where the atoms are written, and fix spring/self tethers each atom to its start");
    if (ramped) notes->push_back("position restraints keep their first force constant in LAMMPS (fix spring/self does not ramp)");
    bool dihedral = false;
    for (const auto& r : rs) dihedral = dihedral || (r.kind == Restraint::Kind::Torsion && r.lower == r.upper);
    if (dihedral) notes->push_back("torsion restraints use LAMMPS's cosine form K [1 − cos(φ − φ0)], equal to CAPS's ½ k Δφ² near φ0 and softer far from it");
  }
  return o.str();
}

std::string gromacs_intermolecular_restraints(const std::vector<Restraint>& rs, std::vector<std::string>* notes) {
  std::ostringstream bonds, angles, dihres;
  int skipped = 0;
  bool ramped = false;
  char b[256];
  for (const auto& r : rs) {
    if (r.k_end >= 0) ramped = true;
    const auto& a = r.atoms;
    switch (r.kind) {
      case Restraint::Kind::Distance: {
        const double up1 = std::min(r.upper, 1e5) / 10;
        std::snprintf(b, sizeof b, "%7u %7u  10  %.8f  %.8f  %.8f  %.6f\n", a[0] + 1, a[1] + 1, r.lower / 10, up1, up1 + 100.0, r.k * kKJ * 100);
        bonds << b;
        break;
      }
      case Restraint::Kind::Angle:
        if (r.lower != r.upper) { ++skipped; break; }
        std::snprintf(b, sizeof b, "%7u %7u %7u  1  %.6f  %.6f\n", a[0] + 1, a[1] + 1, a[2] + 1, r.lower, r.k * kKJ);
        angles << b;
        break;
      case Restraint::Kind::Torsion: {
        double c = 0.5 * (r.lower + r.upper);
        c -= 360 * std::round(c / 360);
        std::snprintf(b, sizeof b, "%7u %7u %7u %7u  1  %.6f  %.6f  %.6f\n", a[0] + 1, a[1] + 1, a[2] + 1, a[3] + 1, c, 0.5 * (r.upper - r.lower), r.k * kKJ);
        dihres << b;
        break;
      }
      case Restraint::Kind::Position: break;
    }
  }
  std::string out;
  if (!bonds.str().empty() || !angles.str().empty() || !dihres.str().empty()) {
    out = "\n[ intermolecular_interactions ]\n";
    if (!bonds.str().empty()) out += "; CAPS distance restraints: ½ k (r − low)² below low, ½ k (r − up1)² above up1 (linear past up2)\n[ bonds ]\n;    ai      aj funct  low (nm)  up1  up2  k (kJ/mol/nm²)\n" + bonds.str();
    if (!angles.str().empty()) out += "; CAPS angle restraints: ½ k (θ − θ0)²\n[ angles ]\n;    ai      aj      ak funct  theta0 (deg)  k (kJ/mol/rad²)\n" + angles.str();
    if (!dihres.str().empty()) out += "; CAPS torsion restraints: ½ k (|φ − φ0| − Δφ)² outside ±Δφ\n[ dihedral_restraints ]\n;    ai      aj      ak      al funct  phi0 (deg)  dphi (deg)  k (kJ/mol/rad²)\n" + dihres.str();
  }
  if (notes) {
    if (skipped) notes->push_back(std::to_string(skipped) + " flat-bottomed angle restraint(s) left out: GROMACS's angle restraints are a cosine between two vectors");
    if (ramped) notes->push_back("restraints keep their first force constant in GROMACS (no ramp)");
  }
  return out;
}

std::string gromacs_position_restraints(const std::vector<Restraint>& rs, uint32_t first, uint32_t last, std::vector<std::string>* notes) {
  std::ostringstream o;
  char b[200];
  int partial_fb = 0;
  for (const auto& r : rs) {
    if (r.kind != Restraint::Kind::Position || r.atoms[0] < first || r.atoms[0] > last) continue;
    const double k = r.k * kKJ * 100;
    if (r.upper > 0) {
      if (r.axes != 7) { ++partial_fb; continue; }
      std::snprintf(b, sizeof b, "%7u  2  1  %.8f  %.6f\n", r.atoms[0] - first + 1, r.upper / 10, k);
    } else {
      std::snprintf(b, sizeof b, "%7u  1  %.6f  %.6f  %.6f\n", r.atoms[0] - first + 1, (r.axes & 1) ? k : 0.0, (r.axes & 2) ? k : 0.0, (r.axes & 4) ? k : 0.0);
    }
    o << b;
  }
  if (notes && partial_fb) notes->push_back(std::to_string(partial_fb) + " flat-bottomed position restraint(s) on some axes only left out (GROMACS's sphere covers all three)");
  const std::string body = o.str();
  return body.empty() ? std::string() : "\n[ position_restraints ]\n; CAPS: ½ k |r − r_ref|² (funct 1: kx ky kz, kJ/mol/nm²) or flat-bottomed sphere (funct 2: g r k)\n" + body;
}

}  // namespace caps
