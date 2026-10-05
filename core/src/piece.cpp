#include "caps/piece.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/json.hpp"

namespace caps {

std::string piece_json(const System& s, const std::vector<size_t>& atoms, const std::string& name) {
  Vec3 c{0, 0, 0};
  std::vector<size_t> list;
  for (size_t i : atoms) if (i < s.atoms.size()) list.push_back(i);
  if (list.empty()) throw std::runtime_error("nothing to copy");
  // the piece whole: each atom at the image nearest the first one (a molecule across a wall stays in one piece)
  std::vector<Vec3> pos(list.size());
  for (size_t k = 0; k < list.size(); ++k) {
    const Vec3 d = s.atoms[list[k]].pos - s.atoms[list[0]].pos;
    pos[k] = s.atoms[list[0]].pos + (s.cell.valid() ? s.cell.minimum_image(d) : d);
    c = c + pos[k];
  }
  c = c * (1.0 / double(list.size()));
  std::map<size_t, size_t> at;
  std::map<int64_t, int> mols;
  const auto molidx = s.molecules();
  std::map<int, std::string> label;
  std::map<int, double> mass;
  for (const auto& t : s.types) label[t.type] = t.label, mass[t.type] = t.mass;
  std::map<std::string, int> el;
  Json a = Json::array();
  for (size_t k = 0; k < list.size(); ++k) {
    const Atom& x = s.atoms[list[k]];
    at[list[k]] = k;
    const int64_t m = s.has_mol ? x.mol : int64_t(molidx[list[k]]) + 1;
    if (!mols.count(m)) { const int next = int(mols.size()) + 1; mols[m] = next; }
    Json o = Json::object();
    o["el"] = element(x.element).symbol;
    o["x"] = pos[k][0] - c[0], o["y"] = pos[k][1] - c[1], o["z"] = pos[k][2] - c[2];
    o["q"] = x.charge;
    o["type"] = x.type;
    o["label"] = label.count(x.type) ? label[x.type] : "";
    o["mass"] = mass.count(x.type) ? mass[x.type] : s.mass_of(x);
    o["name"] = x.name;
    o["res"] = x.resname;
    o["mol"] = mols[m];
    a.push_back(o);
    ++el[element(x.element).symbol];
  }
  Json b = Json::array();
  for (const auto& bd : s.bonds) {
    auto i = at.find(bd.i), j = at.find(bd.j);
    if (i == at.end() || j == at.end()) continue;
    Json e = Json::array();
    e.push_back(double(i->second)), e.push_back(double(j->second)), e.push_back(bd.order);
    b.push_back(e);
  }
  std::string formula;   // Hill order
  auto part = [&](const std::string& e) { if (el.count(e)) { formula += e + (el[e] > 1 ? std::to_string(el[e]) : ""); el.erase(e); } };
  if (el.count("C")) { part("C"); part("H"); }
  for (const auto& [e, n] : el) formula += e + (n > 1 ? std::to_string(n) : "");
  Json j = Json::object();
  j["format"] = "caps-piece";
  j["name"] = name.empty() ? formula : name;
  j["formula"] = formula;
  j["molecules"] = double(mols.size());
  j["atoms"] = a;
  j["bonds"] = b;
  return j.dump_exact();
}

StampResult stamp_piece(System& s, const std::string& piece, const Vec3& at, const Vec3& axis, double degrees, double clear) {
  const Json j = Json::parse(piece);
  if (!j.is_object() || j.text("format") != "caps-piece" || !j.has("atoms")) throw std::runtime_error("not a CAPS piece");
  const auto& A = j["atoms"].items();
  if (A.empty()) throw std::runtime_error("the piece has no atoms");
  // turned about its centre (Rodrigues), then moved to the point
  const double n = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  const Vec3 u = n > 1e-12 ? axis * (1.0 / n) : Vec3{0, 0, 1};
  const double th = degrees * M_PI / 180.0, co = std::cos(th), si = std::sin(th);
  std::vector<Vec3> p;
  for (const auto& x : A) {
    const Vec3 v{x.num("x", 0), x.num("y", 0), x.num("z", 0)};
    const double ud = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
    const Vec3 cr{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
    p.push_back(v * co + cr * si + u * (ud * (1 - co)) + at);
  }
  // the structure's atoms near enough to matter
  double rad = 0;
  for (const auto& q : p) { const Vec3 d = q - at; rad = std::max(rad, std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])); }
  const double reach = rad + 12.0 + clear;
  auto dvec = [&](const Vec3& a, const Vec3& b) { const Vec3 d = a - b; return s.cell.valid() ? s.cell.minimum_image(d) : d; };
  std::vector<Vec3> near;
  for (const auto& x : s.atoms) { const Vec3 d = dvec(x.pos, at); if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= reach * reach) near.push_back(at + d); }
  auto closest = [&](const Vec3& shift) {
    double m2 = 1e300;
    for (const auto& q : p)
      for (const auto& x : near) {
        const Vec3 d = q + shift - x;
        m2 = std::min(m2, d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      }
    return std::sqrt(m2);
  };
  StampResult r;
  Vec3 best{0, 0, 0};
  double got = closest(best);
  if (got < clear && !near.empty()) {
    // the shortest move that clears it: shells of 0.25 Å, 96 directions each (a Fibonacci sphere); no room anywhere
    // within 12 Å (a large piece in a dense cell): the roomiest place found, its close contacts said (Minimise clears them)
    bool found = false;
    Vec3 roomiest = best;
    double room = got;
    for (double dist = 0.25; dist <= 12.0 + 1e-9 && !found; dist += 0.25)
      for (int k = 0; k < 96; ++k) {
        const double zz = 1.0 - 2.0 * (k + 0.5) / 96.0, rr = std::sqrt(1 - zz * zz), ph = k * M_PI * (3.0 - std::sqrt(5.0));
        const Vec3 sh{rr * std::cos(ph) * dist, rr * std::sin(ph) * dist, zz * dist};
        const double c = closest(sh);
        if (c >= clear) { best = sh, got = c, found = true; break; }
        if (c > room) room = c, roomiest = sh;
      }
    if (!found) best = roomiest, got = room;
  }
  r.shift = std::sqrt(best[0] * best[0] + best[1] * best[1] + best[2] * best[2]);
  r.closest = near.empty() ? 0.0 : got;
  // types by label, molecules after the structure's
  int64_t mol0 = 0;
  const auto molidx = s.molecules();
  for (size_t i = 0; i < s.atoms.size(); ++i) mol0 = std::max(mol0, s.has_mol ? s.atoms[i].mol : int64_t(molidx[i]) + 1);
  if (!s.has_mol) for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].mol = int64_t(molidx[i]) + 1;
  int64_t id = 0;
  for (const auto& x : s.atoms) id = std::max(id, x.id);
  const size_t base = s.atoms.size();
  for (size_t k = 0; k < A.size(); ++k) {
    const Json& x = A[k];
    Atom a;
    a.element = element_from_symbol(x.text("el", "C"));
    a.pos = p[k] + best;
    a.charge = x.num("q", 0);
    a.name = x.text("name");
    a.resname = x.text("res");
    a.mol = mol0 + int64_t(x.num("mol", 1));
    a.id = ++id;
    const std::string lab = x.text("label");
    int type = 0;
    if (!lab.empty()) for (const auto& t : s.types) if (t.label == lab) type = t.type;
    if (type == 0 && (!s.types.empty() || x.num("type", 0) > 0)) {
      for (const auto& t : s.types) type = std::max(type, t.type);
      ++type;
      TypeInfo ti;
      ti.type = type, ti.mass = x.num("mass", element(a.element).mass), ti.label = lab.empty() ? element(a.element).symbol : lab;
      s.types.push_back(ti);
    }
    a.type = type;
    r.added.push_back(s.atoms.size());
    s.atoms.push_back(a);
  }
  if (j.has("bonds"))
    for (const auto& b : j["bonds"].items())
      if (b.size() >= 2) {
        const size_t bi = size_t(b[0].number()), bj = size_t(b[1].number());
        if (bi < A.size() && bj < A.size()) s.bonds.push_back({uint32_t(base + bi), uint32_t(base + bj), b.size() > 2 ? int(b[2].number()) : 0});
      }
  s.has_mol = true;
  if (std::any_of(A.begin(), A.end(), [](const Json& x) { return std::abs(x.num("q", 0)) > 0; })) s.has_charges = true;
  return r;
}

}  // namespace caps
