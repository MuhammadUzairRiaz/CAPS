// CAPS backmapping of chemistry-specific coarse-grained melts (see caps/cg_backmap.hpp).
#include "caps/cg_backmap.hpp"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/superpose.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;

bool starts_alpha(const std::string& s) {
  const auto p = s.find_first_not_of(" \t\r");
  return p != std::string::npos && std::isalpha(static_cast<unsigned char>(s[p]));
}
std::string trim(std::string s) {
  const auto c = s.find('#');
  std::string body = c == std::string::npos ? s : s.substr(0, c);
  body.erase(0, body.find_first_not_of(" \t\r"));
  body.erase(body.find_last_not_of(" \t\r") + 1);
  return body;
}
}  // namespace

LammpsFull read_lammps_full(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::invalid_argument("cannot open " + path);
  LammpsFull d;
  std::vector<std::string> L;
  for (std::string l; std::getline(f, l);) L.push_back(l);
  if (L.empty()) throw std::invalid_argument(path + " is empty");
  d.title = L[0];
  double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0}, tilt[3] = {0, 0, 0};
  size_t i = 1;
  // header
  for (; i < L.size(); ++i) {
    const std::string s = trim(L[i]);
    if (s.empty()) continue;
    if (starts_alpha(s)) break;
    std::istringstream ss(s);
    double a, b, c;
    std::string w1, w2, w3;
    if (s.find("xlo") != std::string::npos) { ss >> lo[0] >> hi[0]; continue; }
    if (s.find("ylo") != std::string::npos) { ss >> lo[1] >> hi[1]; continue; }
    if (s.find("zlo") != std::string::npos) { ss >> lo[2] >> hi[2]; continue; }
    if (s.find("xy") != std::string::npos) { ss >> tilt[0] >> tilt[1] >> tilt[2]; continue; }
    ss >> a >> w1 >> w2;
    (void)b, (void)c, (void)w3;
    const int n = int(a);
    if (w1 == "atom" && w2 == "types") d.atom_types = n;
    else if (w1 == "bond" && w2 == "types") d.bond_types = n;
    else if (w1 == "angle" && w2 == "types") d.angle_types = n;
    else if (w1 == "dihedral" && w2 == "types") d.dihedral_types = n;
    else if (w1 == "improper" && w2 == "types") d.improper_types = n;
  }
  d.cell.origin = {lo[0], lo[1], lo[2]};
  d.cell.a = {hi[0] - lo[0], 0, 0};
  d.cell.b = {tilt[0], hi[1] - lo[1], 0};
  d.cell.c = {tilt[1], tilt[2], hi[2] - lo[2]};
  // sections
  while (i < L.size()) {
    const std::string head = L[i];
    const std::string name = trim(head);
    ++i;
    std::vector<std::string> lines;
    for (; i < L.size(); ++i) {
      if (trim(L[i]).empty()) { if (!lines.empty() && i + 1 < L.size() && starts_alpha(L[i + 1])) { ++i; break; } continue; }
      if (starts_alpha(L[i])) break;
      lines.push_back(L[i]);
    }
    if (name == "Masses") d.masses = lines;
    else if (name.size() > 7 && name.compare(name.size() - 6, 6, "Coeffs") == 0) {
      std::string h = head;
      h.erase(h.find_last_not_of(" \t\r") + 1);
      d.coeffs.push_back({h, lines});
    } else if (name == "Atoms") {
      for (const auto& l : lines) {
        std::istringstream ss(trim(l));
        LammpsFull::Atom a;
        double x, y, z;
        int ix = 0, iy = 0, iz = 0;
        ss >> a.id >> a.mol >> a.type >> a.q >> x >> y >> z;
        if (!(ss >> ix >> iy >> iz)) ix = iy = iz = 0;
        a.x = Vec3{x, y, z} + d.cell.a * double(ix) + d.cell.b * double(iy) + d.cell.c * double(iz);
        d.atoms.push_back(a);
      }
      std::sort(d.atoms.begin(), d.atoms.end(), [](const auto& p, const auto& q) { return p.id < q.id; });
    } else if (name == "Bonds" || name == "Angles" || name == "Dihedrals" || name == "Impropers") {
      auto& v = name == "Bonds" ? d.bonds : name == "Angles" ? d.angles : name == "Dihedrals" ? d.dihedrals : d.impropers;
      const size_t k = name == "Bonds" ? 2 : name == "Angles" ? 3 : 4;
      for (const auto& l : lines) {
        std::istringstream ss(trim(l));
        int64_t id;
        LammpsFull::Term t;
        ss >> id >> t.type;
        for (size_t q = 0; q < k; ++q) { int64_t x; ss >> x; t.ids.push_back(x); }
        v.push_back(t);
      }
    }
  }
  if (d.atoms.empty()) throw std::invalid_argument(path + ": no Atoms section (atom_style full is read)");
  return d;
}

void write_lammps_full(const LammpsFull& d, const std::string& path, bool pair_coeffs) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot write " + path);
  char b[256];
  f << (d.title.empty() ? "CAPS" : d.title) << "\n\n";
  f << d.atoms.size() << " atoms\n" << d.bonds.size() << " bonds\n" << d.angles.size() << " angles\n" << d.dihedrals.size() << " dihedrals\n" << d.impropers.size() << " impropers\n\n";
  f << d.atom_types << " atom types\n" << d.bond_types << " bond types\n" << d.angle_types << " angle types\n" << d.dihedral_types << " dihedral types\n" << d.improper_types << " improper types\n\n";
  const Cell& c = d.cell;
  std::snprintf(b, sizeof b, "%.6f %.6f xlo xhi\n%.6f %.6f ylo yhi\n%.6f %.6f zlo zhi\n", c.origin[0], c.origin[0] + c.a[0], c.origin[1], c.origin[1] + c.b[1], c.origin[2], c.origin[2] + c.c[2]);
  f << b;
  if (std::fabs(c.b[0]) + std::fabs(c.c[0]) + std::fabs(c.c[1]) > 1e-12) { std::snprintf(b, sizeof b, "%.6f %.6f %.6f xy xz yz\n", c.b[0], c.c[0], c.c[1]); f << b; }
  f << "\nMasses\n\n";
  for (const auto& m : d.masses) f << m << "\n";
  for (const auto& [h, lines] : d.coeffs) {
    if (!pair_coeffs && (h.rfind("Pair Coeffs", 0) == 0 || h.rfind("PairIJ Coeffs", 0) == 0)) continue;
    f << "\n" << h << "\n\n";
    for (const auto& l : lines) f << l << "\n";
  }
  f << "\nAtoms  # full\n\n";
  for (const auto& a : d.atoms) { std::snprintf(b, sizeof b, "%lld %lld %d %.6f %.6f %.6f %.6f\n", (long long)a.id, (long long)a.mol, a.type, a.q, a.x[0], a.x[1], a.x[2]); f << b; }
  auto terms = [&](const char* name, const std::vector<LammpsFull::Term>& v) {
    if (v.empty()) return;
    f << "\n" << name << "\n\n";
    for (size_t k = 0; k < v.size(); ++k) {
      f << k + 1 << " " << v[k].type;
      for (int64_t x : v[k].ids) f << " " << x;
      f << "\n";
    }
  };
  terms("Bonds", d.bonds);
  terms("Angles", d.angles);
  terms("Dihedrals", d.dihedrals);
  terms("Impropers", d.impropers);
}

// ------------------------------------------------------------------------------------------------ library
namespace {
std::string role_of(size_t p, size_t L) { return L == 1 ? "single" : p == 0 ? "head" : p + 1 == L ? "tail" : "inner"; }

// the two chain neighbours that orient a bead: i−1 and i+1, else the two on its one side
std::pair<int, int> orienting(size_t p, size_t L) {
  if (L < 2) return {-1, -1};
  if (p > 0 && p + 1 < L) return {int(p) - 1, int(p) + 1};
  if (p == 0) return {1, L > 2 ? 2 : -1};
  return {int(p) - 1, p >= 2 ? int(p) - 2 : -1};
}
}  // namespace

BackmapLibrary backmap_library(const LammpsFull& aa, const Json& mj, int conformers) {
  BackmapLibrary lib;
  // a System for the mapping (ids, masses by type, bonds)
  System s;
  std::map<int, double> tmass;
  for (const auto& m : aa.masses) { std::istringstream ss(trim(m)); int t; double w; if (ss >> t >> w) tmass[t] = w; }
  std::map<int64_t, uint32_t> idx;
  for (const auto& a : aa.atoms) {
    Atom x;
    x.id = a.id, x.mol = a.mol, x.type = a.type, x.charge = a.q, x.pos = a.x;
    idx[a.id] = uint32_t(s.atoms.size());
    s.atoms.push_back(x);
  }
  for (const auto& [t, w] : tmass) s.types.push_back({t, w, ""});
  for (const auto& b : aa.bonds) s.bonds.push_back({idx.at(b.ids[0]), idx.at(b.ids[1]), 1});
  s.cell = aa.cell;
  const CgMapping m = cg_mapping_from_json(mj, s);
  auto mass = [&](uint32_t a) { const auto it = tmass.find(s.atoms[a].type); return it != tmass.end() ? it->second : 0.0; };
  // bead centres (whole: positions are unwrapped already, bead atoms made whole about the first)
  std::vector<Vec3> com(m.beads());
  for (size_t b = 0; b < m.beads(); ++b) {
    const Vec3 ref = s.atoms[m.bead_atoms[b].front()].pos;
    Vec3 c{0, 0, 0};
    double w = 0;
    for (uint32_t a : m.bead_atoms[b]) {
      const Vec3 p = ref + (aa.cell.valid() ? aa.cell.minimum_image(s.atoms[a].pos - ref) : s.atoms[a].pos - ref);
      c = c + p * mass(a), w += mass(a);
    }
    com[b] = c * (1.0 / w);
  }
  // bead → (chain, position); class per bead = (class key, role)
  std::vector<int> chain_of(m.beads(), -1), pos_of(m.beads(), -1), cls(m.beads(), -1);
  std::map<std::pair<std::string, std::string>, int> class_ix;
  for (size_t c = 0; c < m.chains.size(); ++c)
    for (size_t p = 0; p < m.chains[c].size(); ++p) {
      const int b = m.chains[c][p];
      chain_of[size_t(b)] = int(c), pos_of[size_t(b)] = int(p);
      const std::string role = role_of(p, m.chains[c].size());
      const auto key = std::make_pair(m.bead_class.empty() || m.bead_class[size_t(b)].empty() ? m.bead_kind[size_t(b)] : m.bead_class[size_t(b)], role);
      auto it = class_ix.find(key);
      if (it == class_ix.end()) {
        BackmapLibrary::ClassInfo ci;
        ci.key = key.first, ci.kind = m.bead_kind[size_t(b)], ci.role = role, ci.atoms = int(m.bead_atoms[size_t(b)].size());
        for (uint32_t a : m.bead_atoms[size_t(b)]) ci.mass += mass(a);
        it = class_ix.emplace(key, int(lib.classes.size())).first;
        lib.classes.push_back(ci);
      }
      cls[size_t(b)] = it->second;
    }
  // fragments: atoms in id order; an instance joins its class when its types match the first one's
  std::map<int, int> skipped;
  for (size_t c = 0; c < m.chains.size(); ++c)
    for (size_t p = 0; p < m.chains[c].size(); ++p) {
      const int b = m.chains[c][p];
      auto& ci = lib.classes[size_t(cls[size_t(b)])];
      if (int(ci.instances.size()) >= conformers) continue;
      auto at = m.bead_atoms[size_t(b)];
      std::sort(at.begin(), at.end(), [&](uint32_t x, uint32_t y) { return s.atoms[x].id < s.atoms[y].id; });
      BackmapLibrary::Fragment fr;
      const Vec3 ref = s.atoms[at.front()].pos;
      for (uint32_t a : at) {
        const Vec3 pa = ref + (aa.cell.valid() ? aa.cell.minimum_image(s.atoms[a].pos - ref) : s.atoms[a].pos - ref);
        fr.type.push_back(s.atoms[a].type), fr.q.push_back(s.atoms[a].charge), fr.rel.push_back(pa - com[size_t(b)]);
      }
      if (!ci.instances.empty() && fr.type != ci.instances.front().type) { skipped[cls[size_t(b)]]++; continue; }
      const auto [n1, n2] = orienting(p, m.chains[c].size());
      auto dir = [&](int q) {
        if (q < 0) return Vec3{0, 0, 0};
        const Vec3 d = com[size_t(m.chains[c][size_t(q)])] - com[size_t(b)];
        return aa.cell.valid() ? aa.cell.minimum_image(d) : d;
      };
      fr.dirs = {dir(n1), dir(n2)};
      ci.instances.push_back(std::move(fr));
    }
  for (const auto& [k, n] : skipped) lib.notes.push_back(std::to_string(n) + " beads of class " + lib.classes[size_t(k)].kind + " (" + lib.classes[size_t(k)].role + ") list their atoms in another order: not used as conformers");
  // (kind, role) → the most common class
  std::map<std::pair<std::string, std::string>, std::map<int, int>> votes;
  for (size_t b = 0; b < m.beads(); ++b) votes[{m.bead_kind[b], lib.classes[size_t(cls[b])].role}][cls[b]]++;
  for (const auto& [kr, v] : votes) {
    lib.by_kind_role[kr] = std::max_element(v.begin(), v.end(), [](auto& x, auto& y) { return x.second < y.second; })->first;
    if (v.size() > 1) lib.notes.push_back("several fragment classes for " + kr.first + " (" + kr.second + "): the most common is used");
  }
  // templates: every term of the chains by the classes of the consecutive beads it spans and each atom's role
  std::vector<int> bead_of(s.atoms.size(), -1), role(s.atoms.size(), -1);
  for (size_t b = 0; b < m.beads(); ++b) {
    auto at = m.bead_atoms[b];
    std::sort(at.begin(), at.end(), [&](uint32_t x, uint32_t y) { return s.atoms[x].id < s.atoms[y].id; });
    for (size_t r = 0; r < at.size(); ++r) bead_of[at[r]] = int(b), role[at[r]] = int(r);
  }
  std::set<std::string> seen;
  int across_chains = 0;
  auto add = [&](int kind, const LammpsFull::Term& t) {
    std::vector<int> beads;
    for (int64_t id : t.ids) beads.push_back(bead_of[idx.at(id)]);
    const int ch = chain_of[size_t(beads[0])];
    int pmin = 1 << 30, pmax = -1;
    for (int b : beads) {
      if (chain_of[size_t(b)] != ch) { ++across_chains; return; }
      pmin = std::min(pmin, pos_of[size_t(b)]), pmax = std::max(pmax, pos_of[size_t(b)]);
    }
    BackmapLibrary::Template tp;
    tp.kind = kind, tp.type = t.type;
    for (int p = pmin; p <= pmax; ++p) tp.classes.push_back(cls[size_t(m.chains[size_t(ch)][size_t(p)])]);
    for (size_t k = 0; k < t.ids.size(); ++k) {
      const uint32_t a = idx.at(t.ids[k]);
      tp.atoms.push_back({pos_of[size_t(bead_of[a])] - pmin, role[a]});
    }
    std::string key = std::to_string(kind) + ":" + std::to_string(tp.type) + ":";
    for (int c : tp.classes) key += std::to_string(c) + ",";
    key += ":";
    for (const auto& [o, r] : tp.atoms) key += std::to_string(o) + "." + std::to_string(r) + ",";
    if (seen.insert(key).second) lib.templates.push_back(std::move(tp));
  };
  for (const auto& t : aa.bonds) add(0, t);
  for (const auto& t : aa.angles) add(1, t);
  for (const auto& t : aa.dihedrals) add(2, t);
  for (const auto& t : aa.impropers) add(3, t);
  if (across_chains) lib.notes.push_back(std::to_string(across_chains) + " terms join two chains (crosslinks): not templated");
  lib.reference = aa;
  lib.reference.atoms.clear(), lib.reference.bonds.clear(), lib.reference.angles.clear(), lib.reference.dihedrals.clear(), lib.reference.impropers.clear();
  return lib;
}

// ------------------------------------------------------------------------------------------------ placement
BackmapResult backmap_fragments(const BackmapLibrary& lib, const CgTopology& t, const std::vector<Vec3>& beads, const Cell& cell, uint64_t seed) {
  if (beads.size() != t.beads()) throw std::invalid_argument("the bead positions do not match the CG topology");
  BackmapResult R;
  R.data = lib.reference;
  R.data.cell = cell;
  R.data.title = "CAPS backmapped from " + std::to_string(t.beads()) + " coarse-grained beads";
  std::mt19937_64 g(seed);
  auto unit01 = [&]() { return double(g() >> 11) * (1.0 / 9007199254740992.0); };
  int64_t next_id = 1;
  std::vector<double> cut;
  // per chain: the classes (reversed when the ends only match that way)
  for (size_t c = 0; c < t.chains.size(); ++c) {
    std::vector<int> chain = t.chains[c];
    const size_t L = chain.size();
    auto classes_of = [&](const std::vector<int>& ch, std::vector<int>& out) {
      out.clear();
      for (size_t p = 0; p < ch.size(); ++p) {
        const auto it = lib.by_kind_role.find({t.kind[size_t(ch[p])], role_of(p, ch.size())});
        if (it == lib.by_kind_role.end()) return false;
        out.push_back(it->second);
      }
      return true;
    };
    std::vector<int> cls;
    if (!classes_of(chain, cls)) {
      std::reverse(chain.begin(), chain.end());
      if (!classes_of(chain, cls)) {
        std::string k = t.kind[size_t(chain.back())] + " … " + t.kind[size_t(chain.front())];
        throw std::invalid_argument("chain " + std::to_string(c + 1) + " (" + k + "): the reference has no fragment for one of its beads in that place (a kind at a chain end the reference chains never end with)");
      }
      ++R.reversed_chains;
    }
    // atoms of each bead
    std::vector<std::vector<int64_t>> ids(L);
    for (size_t p = 0; p < L; ++p) {
      const auto& ci = lib.classes[size_t(cls[p])];
      const auto& fr = ci.instances[size_t(unit01() * double(ci.instances.size())) % ci.instances.size()];
      // the target bead's orienting neighbours (whole along the chain by minimum images)
      const auto [n1, n2] = orienting(p, L);
      auto dir = [&](int q) {
        if (q < 0) return Vec3{0, 0, 0};
        const Vec3 d = beads[size_t(chain[size_t(q)])] - beads[size_t(chain[p])];
        return cell.valid() ? cell.minimum_image(d) : d;
      };
      const std::array<Vec3, 2> td = {dir(n1), dir(n2)};
      std::array<std::array<double, 3>, 3> R3{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
      if (n1 >= 0 && n2 >= 0) {
        const Superposition fit = superpose({Vec3{0, 0, 0}, td[0], td[1]}, {Vec3{0, 0, 0}, fr.dirs[0], fr.dirs[1]});
        R3 = fit.rot;
      } else {
        // one neighbour (a dimer) or none: a random turn, the neighbour direction matched when there is one
        const double u1 = unit01(), u2 = unit01(), u3 = unit01();
        const double q0 = std::sqrt(1 - u1) * std::sin(2 * kPi * u2), q1 = std::sqrt(1 - u1) * std::cos(2 * kPi * u2), q2 = std::sqrt(u1) * std::sin(2 * kPi * u3), q3 = std::sqrt(u1) * std::cos(2 * kPi * u3);
        R3 = {{{1 - 2 * (q2 * q2 + q3 * q3), 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)},
               {2 * (q1 * q2 + q0 * q3), 1 - 2 * (q1 * q1 + q3 * q3), 2 * (q2 * q3 - q0 * q1)},
               {2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), 1 - 2 * (q1 * q1 + q2 * q2)}}};
      }
      for (size_t k = 0; k < fr.type.size(); ++k) {
        LammpsFull::Atom a;
        a.id = next_id++;
        a.mol = int64_t(c + 1);
        a.type = fr.type[k];
        a.q = fr.q[k];
        const Vec3& v = fr.rel[k];
        const Vec3 rv{R3[0][0] * v[0] + R3[0][1] * v[1] + R3[0][2] * v[2], R3[1][0] * v[0] + R3[1][1] * v[1] + R3[1][2] * v[2], R3[2][0] * v[0] + R3[2][1] * v[1] + R3[2][2] * v[2]};
        a.x = beads[size_t(chain[p])] + rv;
        R.charge += a.q;
        ids[p].push_back(a.id);
        R.data.atoms.push_back(a);
      }
    }
    // chains made whole (the beads may be wrapped): each bead's atoms follow the minimum image from the one before
    for (size_t p = 1; p < L; ++p) {
      if (!cell.valid()) break;
      const Vec3 d = beads[size_t(chain[p])] - beads[size_t(chain[p - 1])];
      const Vec3 shift = cell.minimum_image(d) - d;
      if (norm(shift) < 1e-9) continue;
      // shift every later bead of the chain
      for (size_t q = p; q < L; ++q)
        for (int64_t id : ids[q]) { auto& a = R.data.atoms[size_t(id - 1)]; a.x = a.x + shift; }
    }
    // terms by template
    std::vector<char> linked(L > 0 ? L - 1 : 0, 0);
    for (const auto& tp : lib.templates) {
      const size_t w = tp.classes.size();
      if (w > L) continue;
      for (size_t p0 = 0; p0 + w <= L; ++p0) {
        bool match = true;
        for (size_t k = 0; k < w && match; ++k) match = tp.classes[k] == cls[p0 + k];
        if (!match) continue;
        LammpsFull::Term term;
        term.type = tp.type;
        for (const auto& [o, r] : tp.atoms) term.ids.push_back(ids[p0 + size_t(o)][size_t(r)]);
        if (tp.kind == 0) {
          R.data.bonds.push_back(term);
          if (tp.atoms[0].first != tp.atoms[1].first) {
            const size_t lo = p0 + size_t(std::min(tp.atoms[0].first, tp.atoms[1].first));
            linked[lo] = 1;
            cut.push_back(norm(R.data.atoms[size_t(term.ids[1] - 1)].x - R.data.atoms[size_t(term.ids[0] - 1)].x));
          }
        } else if (tp.kind == 1) R.data.angles.push_back(term);
        else if (tp.kind == 2) R.data.dihedrals.push_back(term);
        else R.data.impropers.push_back(term);
      }
    }
    for (char l : linked) R.terms_unmatched += !l;
  }
  R.beads = int(t.beads());
  if (!cut.empty()) {
    R.cut_min = *std::min_element(cut.begin(), cut.end()), R.cut_max = *std::max_element(cut.begin(), cut.end());
    double s = 0;
    for (double x : cut) s += x;
    R.cut_mean = s / double(cut.size());
  }
  if (R.terms_unmatched) R.notes.push_back(std::to_string(R.terms_unmatched) + " bonded bead pairs got no bond: their class sequence does not occur in the reference");
  if (R.reversed_chains) R.notes.push_back(std::to_string(R.reversed_chains) + " chains placed in reverse order (their ends matched the reference that way)");
  for (const auto& n : lib.notes) R.notes.push_back(n);
  return R;
}

std::vector<std::string> style_lines_of(const std::string& input) {
  std::ifstream f(input);
  if (!f) throw std::invalid_argument("cannot open " + input);
  std::vector<std::string> out;
  static const char* keys[] = {"pair_style", "pair_modify", "kspace_style", "kspace_modify", "special_bonds", "bond_style", "angle_style", "dihedral_style", "improper_style", "dielectric"};
  for (std::string l; std::getline(f, l);) {
    const std::string s = trim(l);
    for (const char* k : keys)
      if (s.rfind(k, 0) == 0 && (s.size() == std::strlen(k) || s[std::strlen(k)] == ' ')) { out.push_back(s); break; }
  }
  return out;
}

std::string backmap_deck(const std::vector<std::string>& style, const BackmapDeckOptions& o) {
  auto pick = [&](const char* k) {
    std::string out;
    for (const auto& s : style) if (s.rfind(k, 0) == 0) out += s + "\n";
    return out;
  };
  char b[1200];
  std::string d = "# CAPS backmapping relaxation (LAMMPS): 1 bonded terms only, 2 a soft-core push-off, 3 the full force field, 4 a short NPT\n"
                  "#   lmp -in in.backmap -var DATA backmapped.data -var OUT relaxed\n";
  std::snprintf(b, sizeof b, "variable T index %g\nvariable P index %g\nvariable SOFT index %lld\nvariable NPT index %lld\nvariable SEED index 4928459\n", o.T, o.P,
                (long long)o.soft_steps, (long long)o.npt_steps);
  d += b;
  d += "units real\natom_style full\nboundary p p p\n";
  d += pick("bond_style") + pick("angle_style") + pick("dihedral_style") + pick("improper_style") + pick("special_bonds");
  d += "read_data ${DATA}\nneighbor 2.0 bin\nneigh_modify delay 0 every 1 check yes one 10000\nthermo_style custom step temp press density pe ebond eangle edihed evdwl ecoul\nthermo 1000\n";
  d += "# --- 1 bonded terms only: the fragments joined at the restored bonds\npair_style zero 8.0\npair_coeff * *\nmin_style cg\nmin_modify dmax 0.05\nminimize 1.0e-4 1.0e-6 5000 50000\n";
  std::snprintf(b, sizeof b, "# --- 2 soft-core push-off (overlaps between fragments)\npair_style soft 2.5\npair_coeff * * 0.0\nvariable A equal ramp(0.0,%g)\n"
                "fix push all adapt 1 pair soft a * * v_A\ntimestep %g\nvelocity all create ${T} ${SEED} mom yes rot yes dist gaussian\n"
                "fix lim all nve/limit 0.05\nfix lang all langevin ${T} ${T} 100.0 ${SEED}\nrun ${SOFT}\nunfix push\nunfix lim\nunfix lang\n",
                o.soft_max, o.dt_soft);
  d += b;
  // stage 3 minimises with Ewald when the reference uses PPPM (a PPPM grid fails on the large first moves of a fresh
  // structure), stage 4 runs with the reference's own k-space
  std::string ks = pick("kspace_style"), ks_min = ks;
  if (ks.find("pppm") != std::string::npos) {
    std::istringstream ss(ks);
    std::string w, style, acc;
    ss >> w >> style >> acc;
    ks_min = "kspace_style ewald " + (acc.empty() ? std::string("1.0e-4") : acc) + "\n";
  }
  d += "# --- 3 the full force field (the reference's pair style and coefficients)\n" + pick("pair_style") + "include pair_coeffs.in\n" + pick("pair_modify") + ks_min + pick("dielectric");
  d += "minimize 1.0e-4 1.0e-6 5000 50000\nwrite_data ${OUT}.minimized.data\nreset_timestep 0\n";
  if (ks_min != ks) d += ks + pick("kspace_modify");
  else d += pick("kspace_modify");
  std::snprintf(b, sizeof b, "# --- 4 a short NPT\ntimestep %g\nfix md all npt temp ${T} ${T} 100.0 iso ${P} ${P} 1000.0\nrun ${NPT}\nwrite_data ${OUT}.data\n", o.dt);
  d += b;
  return d;
}

std::vector<BackmapCheckRow> backmap_check(const LammpsFull& d) {
  std::vector<BackmapCheckRow> out;
  std::map<int64_t, size_t> at;
  for (size_t k = 0; k < d.atoms.size(); ++k) at[d.atoms[k].id] = k;
  auto vec = [&](int64_t a, int64_t b) { const Vec3 v = d.atoms[at.at(b)].x - d.atoms[at.at(a)].x; return d.cell.valid() ? d.cell.minimum_image(v) : v; };
  auto coeffs = [&](const char* prefix, std::string* style) {
    std::map<int, double> ref;
    for (const auto& [h, lines] : d.coeffs) {
      if (h.rfind(prefix, 0) != 0) continue;
      const auto c = h.find('#');
      *style = c == std::string::npos ? "" : trim(h.substr(c + 1));
      for (const auto& l : lines) { std::istringstream ss(trim(l)); int t; double k, x0; if (ss >> t >> k >> x0) ref[t] = x0; }
    }
    return ref;
  };
  std::string bs, as;
  const auto br = coeffs("Bond Coeffs", &bs), ar = coeffs("Angle Coeffs", &as);
  std::map<int, std::vector<double>> bd, ad;
  for (const auto& t : d.bonds) bd[t.type].push_back(norm(vec(t.ids[0], t.ids[1])));
  for (const auto& t : d.angles) {
    const Vec3 u = vec(t.ids[1], t.ids[0]), v = vec(t.ids[1], t.ids[2]);
    ad[t.type].push_back(std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / kPi);
  }
  auto rows = [&](const char* kind, const std::map<int, std::vector<double>>& m, const std::map<int, double>& ref, bool harmonic) {
    for (const auto& [t, v] : m) {
      BackmapCheckRow r;
      r.kind = kind, r.type = t, r.count = int(v.size());
      const auto it = ref.find(t);
      if (!harmonic || it == ref.end()) { r.ref = std::nan(""); out.push_back(r); continue; }
      r.ref = it->second;
      double s = 0;
      for (double x : v) s += std::fabs(x - r.ref), r.max_dev = std::max(r.max_dev, std::fabs(x - r.ref));
      r.mean_dev = s / double(v.size());
      out.push_back(r);
    }
  };
  rows("bond", bd, br, bs.empty() || bs == "harmonic");
  rows("angle", ad, ar, as.empty() || as == "harmonic");
  return out;
}

}  // namespace caps
