// caps ff complete (see caps/ffcomplete.hpp): detect what a force field lacks for a structure and fit it with CAPS's
// own machine-learned potential, in the core, with no Python.
#include "caps/ffcomplete.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "caps/elements.hpp"
#include "caps/field.hpp"
#include "caps/mlpot.hpp"
#include "caps/typing.hpp"

namespace caps {
namespace {

constexpr double kEvKcal = 23.0605478;
constexpr double kPi = 3.14159265358979323846;

using Nb = std::vector<std::vector<int>>;

Nb neighbours(size_t n, const std::vector<Bond>& bonds) {
  Nb nb(n);
  for (const auto& b : bonds) {
    nb[b.i].push_back(static_cast<int>(b.j));
    nb[b.j].push_back(static_cast<int>(b.i));
  }
  return nb;
}

// the atoms on c's side of the bond b-c (c reached from itself without going through b)
std::set<int> side(const Nb& nb, int c, int b) {
  std::set<int> seen{c};
  std::vector<int> todo{c};
  while (!todo.empty()) {
    const int x = todo.back();
    todo.pop_back();
    for (int y : nb[x])
      if (y != b && !seen.count(y)) { seen.insert(y); todo.push_back(y); }
  }
  return seen;
}

bool in_ring(const Nb& nb, int b, int c) { return side(nb, c, b).count(b) > 0; }

// the smallest ring through each ring bond, up to `limit` atoms
std::vector<std::set<int>> small_rings(const Nb& nb, size_t limit = 8) {
  std::set<std::set<int>> out;
  for (size_t u0 = 0; u0 < nb.size(); ++u0) {
    const int u = static_cast<int>(u0);
    for (int v : nb[u]) {
      if (u > v) continue;
      std::map<int, int> prev{{u, -1}};
      std::vector<int> layer{u};
      while (!layer.empty() && !prev.count(v) && prev.size() < 4000) {
        std::vector<int> next;
        for (int x : layer)
          for (int y : nb[x])
            if (!prev.count(y) && !(x == u && y == v)) { prev[y] = x; next.push_back(y); }
        layer = next;
      }
      if (!prev.count(v)) continue;
      std::set<int> path;
      for (int x = v; x != -1; x = prev[x]) path.insert(x);
      if (path.size() <= limit) out.insert(path);
    }
  }
  return {out.begin(), out.end()};
}

int valence(int z) {
  switch (z) {
    case 1: case 9: case 17: case 35: case 53: return 1;
    case 8: case 16: return 2;
    case 7: case 15: return 3;
    case 6: case 14: return 4;
    default: return 0;
  }
}

double cap_length(int z) {   // X-H, as caps.torsionfit caps a cut bond
  switch (z) {
    case 7: return 1.01; case 8: return 0.96; case 16: return 1.34; case 14: return 1.48; case 15: return 1.42;
    default: return 1.09;
  }
}

struct Fragment {
  System sys;                   // non-periodic, bonds of unknown order (as a mol2 with "un" bonds)
  std::map<int, int> idx;       // structure atom -> fragment atom
  std::vector<int> z;
  std::string problem;
};

// the fragment round `seed` (structure atoms), `radius` bonds out, grown as caps.census grows it: whole rings of up to
// 8 atoms, no bond between two unsaturated atoms cut, no atom left that two kept atoms bond to; every H of a kept atom;
// each remaining cut bond capped by H; whole across the periodic cell; checked for an even electron count and contacts
Fragment cut(const System& s, const Nb& nb, const std::vector<std::set<int>>& rings, const std::vector<int>& seed, int radius) {
  auto z = [&](int a) { return s.atoms[a].element; };
  auto unsat = [&](int a) { return valence(z(a)) > static_cast<int>(nb[a].size()); };
  std::map<int, int> dist;
  std::vector<int> todo;
  for (int a : seed) { dist[a] = 0; todo.push_back(a); }
  for (size_t k = 0; k < todo.size(); ++k) {
    const int x = todo[k];
    for (int y : nb[x])
      if (!dist.count(y) && dist[x] < radius) { dist[y] = dist[x] + 1; todo.push_back(y); }
  }
  std::set<int> keep;
  for (const auto& [a, d] : dist) keep.insert(a);
  Fragment f;
  for (bool grown = true; grown && keep.size() <= 80;) {
    grown = false;
    for (const auto& r : rings) {
      bool any = false, all = true;
      for (int a : r) { any |= keep.count(a) > 0; all &= keep.count(a) > 0; }
      if (any && !all) { keep.insert(r.begin(), r.end()); grown = true; }
    }
    for (int x : std::vector<int>(keep.begin(), keep.end()))
      for (int y : nb[x]) {
        if (keep.count(y) || z(y) == 1) continue;
        int kept = 0;
        for (int v : nb[y]) kept += keep.count(v) > 0;
        if ((unsat(x) && unsat(y)) || kept > 1) { keep.insert(y); grown = true; }
      }
  }
  if (keep.size() > 80) f.problem = "the fragment would grow past 80 atoms";
  for (int x : std::vector<int>(keep.begin(), keep.end()))
    for (int y : nb[x])
      if (z(y) == 1) keep.insert(y);
  std::map<int, Vec3> whole{{seed[0], s.atoms[seed[0]].pos}};
  std::vector<int> st{seed[0]};
  while (!st.empty()) {
    const int x = st.back();
    st.pop_back();
    for (int y : nb[x])
      if (keep.count(y) && !whole.count(y)) {
        Vec3 d = s.atoms[y].pos - s.atoms[x].pos;
        if (s.cell.valid()) d = s.cell.minimum_image(d);
        whole[y] = whole[x] + d;
        st.push_back(y);
      }
  }
  for (const auto& [a, p] : whole) {
    Atom at;
    at.element = z(a);
    at.pos = p;
    at.id = static_cast<int64_t>(f.sys.atoms.size() + 1);
    at.name = element(z(a)).symbol;
    f.idx[a] = static_cast<int>(f.sys.atoms.size());
    f.sys.atoms.push_back(at);
    f.z.push_back(z(a));
  }
  for (const auto& [a, p] : whole)
    for (int y : nb[a]) {
      if (f.idx.count(y)) {
        if (f.idx[a] < f.idx[y]) f.sys.bonds.push_back({static_cast<uint32_t>(f.idx[a]), static_cast<uint32_t>(f.idx[y]), 0});
      } else {
        Vec3 d = s.atoms[y].pos - s.atoms[a].pos;
        if (s.cell.valid()) d = s.cell.minimum_image(d);
        Atom h;
        h.element = 1;
        h.pos = p + d * (cap_length(z(a)) / norm(d));
        h.id = static_cast<int64_t>(f.sys.atoms.size() + 1);
        h.name = "H";
        f.sys.bonds.push_back({static_cast<uint32_t>(f.idx[a]), static_cast<uint32_t>(f.sys.atoms.size()), 0});
        f.sys.atoms.push_back(h);
        f.z.push_back(1);
      }
    }
  if (f.problem.empty()) {
    int electrons = 0;
    for (int e : f.z) electrons += e;
    if (electrons % 2) f.problem = "odd electron count";
  }
  for (size_t i = 0; i < f.sys.atoms.size() && f.problem.empty(); ++i)
    for (size_t j = i + 1; j < f.sys.atoms.size(); ++j)
      if (norm(f.sys.atoms[i].pos - f.sys.atoms[j].pos) < 0.9) { f.problem = "two atoms closer than 0.9 A"; break; }
  return f;
}

// a missing interaction as the parameter report gives it: "dihedral opls_465(C) opls_467(OS) ..."
struct Missing {
  std::string kind;
  std::vector<std::string> types, classes;
};

Missing parse_missing(const std::string& line) {
  std::istringstream in(line);
  Missing m;
  in >> m.kind;
  std::string w;
  while (in >> w) {
    const auto p = w.find('(');
    if (p == std::string::npos) { m.types.push_back(w); m.classes.push_back(w); continue; }
    m.types.push_back(w.substr(0, p));
    m.classes.push_back(w.substr(p + 1, w.size() - p - 2));
  }
  return m;
}

std::vector<std::string> canonical(const std::vector<std::string>& v) {
  std::vector<std::string> r(v.rbegin(), v.rend());
  return std::min(v, r);
}

std::string join(const std::vector<std::string>& v) {
  std::string s;
  for (size_t k = 0; k < v.size(); ++k) s += (k ? "-" : "") + v[k];
  return s;
}

double six(double x) { return std::round(x * 1e6) / 1e6; }   // coefficients are written to six decimals

std::string fmt(double x, int p = 2) {
  char b[64];
  std::snprintf(b, sizeof b, "%.*f", p, x);
  return b;
}

// a small dense weighted least-squares solve, (A^T W A + diag(ridge)) x = A^T W y
std::vector<double> lsq(const std::vector<std::vector<double>>& A, const std::vector<double>& y, const std::vector<double>& w,
                        const std::vector<double>& ridge) {
  const size_t m = A.empty() ? 0 : A[0].size();
  std::vector<std::vector<double>> M(m, std::vector<double>(m, 0.0));
  std::vector<double> b(m, 0.0);
  for (size_t r = 0; r < A.size(); ++r)
    for (size_t i = 0; i < m; ++i) {
      b[i] += w[r] * A[r][i] * y[r];
      for (size_t j = 0; j < m; ++j) M[i][j] += w[r] * A[r][i] * A[r][j];
    }
  for (size_t i = 0; i < m; ++i) M[i][i] += ridge[i];
  for (size_t c = 0; c < m; ++c) {
    size_t p = c;
    for (size_t r = c + 1; r < m; ++r) if (std::fabs(M[r][c]) > std::fabs(M[p][c])) p = r;
    std::swap(M[c], M[p]);
    std::swap(b[c], b[p]);
    if (std::fabs(M[c][c]) < 1e-14) continue;
    for (size_t r = c + 1; r < m; ++r) {
      const double f = M[r][c] / M[c][c];
      for (size_t q = c; q < m; ++q) M[r][q] -= f * M[c][q];
      b[r] -= f * b[c];
    }
  }
  std::vector<double> x(m, 0.0);
  for (size_t c = m; c-- > 0;) {
    double s = b[c];
    for (size_t q = c + 1; q < m; ++q) s -= M[c][q] * x[q];
    x[c] = std::fabs(M[c][c]) < 1e-14 ? 0.0 : s / M[c][c];
  }
  return x;
}

double dihedral(const std::vector<Vec3>& x, int a, int b, int c, int d) {
  const Vec3 b1 = x[b] - x[a], b2 = x[c] - x[b], b3 = x[d] - x[c];
  const Vec3 m = cross(b1, b2), n = cross(b2, b3);
  return std::atan2(norm(b2) * dot(b1, n), dot(m, n));
}

double angle(const std::vector<Vec3>& x, int i, int j, int k) {
  const Vec3 u = x[i] - x[j], v = x[k] - x[j];
  return std::acos(std::max(-1.0, std::min(1.0, dot(u, v) / (norm(u) * norm(v)))));
}

Vec3 rotate(const Vec3& p, const Vec3& origin, const Vec3& axis, double t) {   // Rodrigues; axis of unit length
  const Vec3 v = p - origin;
  const double c = std::cos(t), s = std::sin(t);
  return origin + v * c + cross(axis, v) * s + axis * (dot(axis, v) * (1 - c));
}

}  // namespace

CompleteReport complete_forcefield(const System& s, const FFDef& ff, const CompleteOptions& opt,
                                   const std::function<void(const std::string&)>& log) {
  auto say = [&](const std::string& m) { if (log) log(m); };
  CompleteReport rep;
  const TypingResult tr = assign_types(s, ff);
  if (tr.untyped) throw std::runtime_error(std::to_string(tr.untyped) + " atom(s) match no typing rule of " + ff.name);
  const auto& types = tr.types;
  ParamReport pr;
  (void)parameterize(s, ff, types, "types", &pr, true);
  rep.missing_before = pr.missing;
  if (pr.missing.empty()) {
    rep.complete = true;
    return rep;
  }

  // a type's class for one kind of term (the name the force field looks that term up with)
  std::map<std::string, const FFType*> by_name;
  for (const auto& t : ff.types) by_name[t.name] = &t;
  auto cls = [&](const std::string& type, const std::string& kind) {
    auto it = by_name.find(type);
    if (it == by_name.end()) return type;
    auto e = it->second->equiv.find(kind);
    return e != it->second->equiv.end() && !e->second.empty() ? e->second : type;
  };

  const Nb nb = neighbours(s.atoms.size(), s.bonds);
  const auto rings = small_rings(nb);
  ml::Potential pot(opt.model);
  std::unique_ptr<ml::Potential> pot2;   // a second copy, for the backward torsion scan run alongside the forward one
  const std::string label = "CAPS MACE v1, " + opt.model.substr(opt.model.find_last_of("/\\") + 1);

  // one term per class, in fitting order (bonds, angles, dihedrals); other kinds are reported
  std::map<std::string, std::vector<Missing>> by_kind;
  std::set<std::pair<std::string, std::vector<std::string>>> seen;
  for (const auto& line : pr.missing) {
    Missing m = parse_missing(line);
    if (m.kind != "bond" && m.kind != "angle" && m.kind != "dihedral") {
      rep.terms.push_back({m.kind, m.classes, "not_fitted", {}, 0, 0, 0, m.kind + " terms are not fitted by caps ff complete"});
      continue;
    }
    if (seen.insert({m.kind, canonical(m.classes)}).second) by_kind[m.kind].push_back(m);
  }
  FFDef& ov = rep.overlay;
  ov.name = "CAPS ff complete";
  ov.units = ff.units;
  FFDef current = ff;   // the force field with every term fitted so far, for the MM energy of the next scans

  EnergyOptions eo;
  eo.cutoff = ff.cutoff > 0 ? ff.cutoff : 10.0;
  auto mm_energy = [&](const Fragment& f, const std::vector<std::string>& ftypes, const std::vector<Vec3>& x) {
    System fs = f.sys;
    for (size_t k = 0; k < x.size(); ++k) fs.atoms[k].pos = x[k];
    ParamReport r;
    const ForceField fld = parameterize(fs, current, ftypes, "types", &r, true);
    Evaluator ev(fld, eo);
    std::vector<double> flat(3 * x.size()), g;
    for (size_t k = 0; k < x.size(); ++k)
      for (int d = 0; d < 3; ++d) flat[3 * k + d] = x[k][d];
    return ev.compute(flat, fs.cell, g).total();
  };

  // an occurrence (structure atoms) of a missing term, by the report's types in either direction
  auto locate = [&](const Missing& m) -> std::vector<int> {
    auto match = [&](const std::vector<int>& q) {
      std::vector<std::string> t;
      for (int a : q) t.push_back(types[a]);
      std::vector<std::string> r(t.rbegin(), t.rend());
      return t == m.types || r == m.types;
    };
    for (size_t b0 = 0; b0 < nb.size(); ++b0) {
      const int b = static_cast<int>(b0);
      for (int c : nb[b]) {
        if (m.kind == "bond" && b < c && match({b, c})) return {b, c};
        if (m.kind == "angle")
          for (int k : nb[b])
            if (k != c && match({c, b, k})) return {c, b, k};
        if (m.kind == "dihedral")
          for (int a : nb[b])
            for (int d : nb[c])
              if (a != c && d != b && a != d && match({a, b, c, d})) return {a, b, c, d};
      }
    }
    return {};
  };

  // the fragment that keeps the types of the term's atoms and of every neighbour of the seed atoms ("" when found)
  auto fragment_for = [&](const std::vector<int>& q, const std::vector<int>& seed, Fragment& f, std::vector<std::string>& ftypes) {
    std::set<int> check(q.begin(), q.end());
    for (int a : seed) check.insert(nb[a].begin(), nb[a].end());
    std::string why = "no fragment keeps the structure's types";
    for (int r = opt.radius; r <= opt.radius + 4; ++r) {
      f = cut(s, nb, rings, seed, r);
      if (!f.problem.empty()) { why = f.problem; continue; }
      if (!pot.covers(f.z)) return std::string("an element outside CAPS's model (it covers H C N O F Si S Cl): needs QM");
      const TypingResult ft = assign_types(f.sys, ff);
      bool same = !ft.untyped;
      for (int a : check) same = same && f.idx.count(a) && ft.types[f.idx.at(a)] == types[a];
      if (same) { ftypes = ft.types; return std::string(); }
    }
    return why;
  };

  auto relaxed = [&](const Fragment& f, const std::vector<Vec3>& x, std::vector<ml::Hold> holds) {
    for (auto& h : holds) h.target = ml::hold_value(h, x);
    return ml::relax(pot, f.z, x, holds, opt.fmax, opt.max_steps);
  };
  auto positions = [](const Fragment& f) {
    std::vector<Vec3> x;
    for (const auto& a : f.sys.atoms) x.push_back(a.pos);
    return x;
  };

  // ------------------------------------------------------------------ bonds and angles
  for (const std::string kind : {"bond", "angle"}) {
    for (const Missing& m : by_kind[kind]) {
      CompleteTerm term{kind, m.classes, "not_fitted", {}, 0, 0, 0, ""};
      const auto q = locate(m);
      Fragment f;
      std::vector<std::string> ft;
      const std::string why = q.empty() ? std::string("no occurrence found") : fragment_for(q, q, f, ft);
      if (!why.empty()) {
        term.why = why;
        rep.terms.push_back(term);
        say(kind + " " + join(m.classes) + ": not fitted (" + why + ")");
        continue;
      }
      term.fragment_atoms = static_cast<int>(f.z.size());
      std::vector<int> fq;
      for (int a : q) fq.push_back(f.idx[a]);
      const std::vector<Vec3> x0 = relaxed(f, positions(f), {}).pos;
      // every occurrence of the class in the fragment (a class met twice is fitted as the sum of its terms)
      const auto key = canonical(m.classes);
      const Nb fnb = neighbours(f.z.size(), f.sys.bonds);
      std::vector<std::vector<int>> occ;
      for (size_t j0 = 0; j0 < fnb.size(); ++j0) {
        const int j = static_cast<int>(j0);
        for (int i : fnb[j]) {
          if (kind == "bond") {
            if (j < i && canonical({cls(ft[j], "bond"), cls(ft[i], "bond")}) == key) occ.push_back({j, i});
          } else {
            for (int k : fnb[j])
              if (i < k && canonical({cls(ft[i], "angle"), cls(ft[j], "angle"), cls(ft[k], "angle")}) == key) occ.push_back({i, j, k});
          }
        }
      }
      // bonds: r_e +- 0.06 A in 0.02 A steps; angles: theta_e +- 12 deg in 3 deg steps
      const int lim = kind == "bond" ? 3 : 4;
      std::vector<std::vector<double>> A;
      std::vector<double> y, w;
      std::vector<CompleteTerm::Point> pts;
      for (int st = -lim; st <= lim; ++st) {
        std::vector<Vec3> x = x0;
        if (kind == "bond") {   // both atoms moved half the change along the bond
          const double dlt = 0.02 * st;
          const Vec3 u = (x[fq[1]] - x[fq[0]]) * (1.0 / norm(x[fq[1]] - x[fq[0]]));
          x[fq[0]] = x[fq[0]] - u * (dlt / 2);
          x[fq[1]] = x[fq[1]] + u * (dlt / 2);
        } else {                 // the third atom turned about the vertex, in the angle's plane
          Vec3 nrm = cross(x[fq[0]] - x[fq[1]], x[fq[2]] - x[fq[1]]);
          nrm = nrm * (1.0 / norm(nrm));
          x[fq[2]] = rotate(x[fq[2]], x[fq[1]], nrm, 3.0 * st * kPi / 180.0);
        }
        ml::Hold h;
        h.n = static_cast<int>(fq.size());
        for (size_t k = 0; k < fq.size(); ++k) h.idx[k] = fq[k];
        const auto r = relaxed(f, x, {h});
        const double emm = mm_energy(f, ft, r.pos);
        double s1 = 0, s2 = 0;
        for (const auto& o : occ) {
          const double v = kind == "bond" ? norm(r.pos[o[0]] - r.pos[o[1]]) : angle(r.pos, o[0], o[1], o[2]);
          s1 += v;
          s2 += v * v;
        }
        A.push_back({s2, s1, 1.0});
        y.push_back(r.energy * kEvKcal - emm);
        w.push_back(1.0);
        const double v = kind == "bond" ? norm(r.pos[fq[0]] - r.pos[fq[1]]) : angle(r.pos, fq[0], fq[1], fq[2]) * 180.0 / kPi;
        pts.push_back({v, r.energy * kEvKcal, emm, r.converged});
      }
      {
        double a = 1e300, b = 1e300;
        for (const auto& p : pts) { a = std::min(a, p.e_ml); b = std::min(b, p.e_mm); }
        for (auto& p : pts) { p.e_ml -= a; p.e_mm -= b; }
        term.scan = pts;
      }
      // E = K sum (x - x0)^2 + const = K sum x^2 - 2 K x0 sum x + const
      const auto sol = lsq(A, y, w, {0, 0, 0});
      const double K = sol[0], x0f = K != 0 ? -sol[1] / (2 * K) : std::nan("");
      double rss = 0;
      for (size_t r = 0; r < A.size(); ++r) {
        const double e = A[r][0] * sol[0] + A[r][1] * sol[1] + sol[2] - y[r];
        rss += e * e;
      }
      term.rms = std::sqrt(rss / static_cast<double>(A.size()));
      term.points = static_cast<int>(A.size());
      const double xv = kind == "bond" ? x0f : x0f * 180.0 / kPi;
      term.params = {K, xv};
      const bool sane = kind == "angle" ? (K >= 20 && K <= 250 && xv >= 90 && xv <= 150) : (K >= 100 && K <= 1500 && xv >= 0.8 && xv <= 2.6);
      if (!sane || term.rms > opt.stretch_rms_max) {
        term.status = "rejected";
        term.why = !sane ? "outside the physical range" : "a poor fit (RMS " + fmt(term.rms) + " kcal/mol)";
        rep.terms.push_back(term);
        say(kind + " " + join(m.classes) + ": K " + fmt(K, 1) + ", x0 " + fmt(xv, 3) + " rejected: " + term.why);
        continue;
      }
      FFRule rule;
      rule.name = join(m.classes);
      rule.match = m.classes;
      rule.style = "harmonic";
      rule.params = {six(K), six(xv)};
      rule.comment = "CAPS ff complete (" + label + "): " + std::to_string(term.points) + "-point relaxed scan of the " + kind +
                     " in a " + std::to_string(term.fragment_atoms) + "-atom fragment, RMS " + fmt(term.rms) + " kcal/mol";
      (kind == "bond" ? ov.bonds : ov.angles).push_back(rule);
      current = ff;
      merge_forcefield(current, ov);
      term.status = "fitted";
      rep.terms.push_back(term);
      say(kind + " " + join(m.classes) + ": K " + fmt(K, 3) + ", x0 " + fmt(xv, 4) + ", RMS " + fmt(term.rms));
    }
  }

  // ------------------------------------------------------------------ torsions: every missing class about the bond together
  std::set<std::vector<std::string>> missing_dih, done;
  for (const Missing& m : by_kind["dihedral"]) missing_dih.insert(canonical(m.classes));
  for (const Missing& m : by_kind["dihedral"]) {
    if (done.count(canonical(m.classes))) continue;
    const auto q = locate(m);
    Fragment f;
    std::vector<std::string> ft;
    std::string why = q.empty() ? "no occurrence found" : "";
    if (why.empty() && in_ring(nb, q[1], q[2])) why = "the central bond is in a ring (not scanned)";
    if (why.empty()) why = fragment_for(q, {q[1], q[2]}, f, ft);
    if (!why.empty()) {
      rep.terms.push_back({"dihedral", m.classes, "not_fitted", {}, 0, 0, 0, why});
      done.insert(canonical(m.classes));
      say("dihedral " + join(m.classes) + ": not fitted (" + why + ")");
      continue;
    }
    const Nb fnb = neighbours(f.z.size(), f.sys.bonds);
    const int b = f.idx[q[1]], c = f.idx[q[2]];
    std::vector<std::vector<std::string>> group;   // the missing classes about b-c, fitted together
    std::map<std::vector<std::string>, std::vector<std::array<int, 4>>> occ;
    for (int a : fnb[b])
      for (int d : fnb[c]) {
        if (a == c || d == b || a == d) continue;
        const auto k = canonical({cls(ft[a], "dihedral"), cls(ft[b], "dihedral"), cls(ft[c], "dihedral"), cls(ft[d], "dihedral")});
        if (!missing_dih.count(k) || done.count(k)) continue;
        if (!occ.count(k)) group.push_back(k);
        occ[k].push_back({a, b, c, d});
      }
    if (group.empty()) {
      rep.terms.push_back({"dihedral", m.classes, "not_fitted", {}, 0, 0, static_cast<int>(f.z.size()), "its class is not found about the bond in the fragment"});
      done.insert(canonical(m.classes));
      continue;
    }
    const std::array<int, 4> quad = {f.idx[q[0]], b, c, f.idx[q[3]]};
    // the smaller side turns: c's side by +delta changes the dihedral by +delta, b's side by -delta
    std::set<int> moving = side(fnb, c, b);
    const std::set<int> other = side(fnb, b, c);
    int sign = 1;
    if (moving.size() > other.size()) { moving = other; sign = -1; }
    const std::vector<Vec3> x0 = relaxed(f, positions(f), {}).pos;
    const double phi0 = dihedral(x0, quad[0], quad[1], quad[2], quad[3]);
    std::vector<std::set<int>> near(f.z.size());   // 1-2 and 1-3 partners, for the clash test
    for (size_t i = 0; i < f.z.size(); ++i)
      for (int j : fnb[i]) {
        near[i].insert(j);
        for (int k : fnb[j]) near[i].insert(k);
      }
    // sequential scans round the full circle from the relaxed start, forward and backward: each point starts from the
    // last relaxed one turned to the next angle (small turns, few clashes), and at each angle the lower of the two
    // relaxed energies is kept, so a scan that followed one hydrogen-bond pattern past where another is lower does not
    // leave a jump in the profile
    std::vector<std::vector<double>> A;
    std::vector<double> y, e_ml;
    std::vector<CompleteTerm::Point> pts;
    std::vector<std::vector<double>> rows;
    std::vector<double> ys;
    const int npts = (360 + opt.step - 1) / opt.step;
    // the two directions run at once, each with its own copy of the model
    std::vector<std::vector<int>> have(2, std::vector<int>(npts, 0));
    std::vector<std::vector<ml::RelaxResult>> got(2, std::vector<ml::RelaxResult>(npts));
    if (!pot2) pot2 = std::make_unique<ml::Potential>(opt.model);
    auto chain = [&](int dir, const ml::Potential& P, std::vector<int>& hv, std::vector<ml::RelaxResult>& gv) {
      std::vector<Vec3> last = x0;
      for (int k = 0; k < npts; ++k) {
        const int slot = ((dir * k) % npts + npts) % npts;
        const double target = phi0 + dir * k * opt.step * kPi / 180.0;
        std::vector<Vec3> x = last;
        Vec3 axis = x[c] - x[b];
        axis = axis * (1.0 / norm(axis));
        const double delta = target - dihedral(x, quad[0], quad[1], quad[2], quad[3]);
        for (int a : moving) x[a] = rotate(x[a], x[c], axis, sign * delta);
        bool clash = false;
        for (size_t i = 0; i < x.size() && !clash; ++i)
          for (size_t j = i + 1; j < x.size(); ++j)
            if (!near[i].count(static_cast<int>(j)) && norm(x[i] - x[j]) < 1.0) { clash = true; break; }
        if (clash) continue;
        ml::Hold h;
        h.n = 4;
        for (int q4 = 0; q4 < 4; ++q4) h.idx[q4] = quad[q4];
        h.target = ml::hold_value(h, x);
        const auto r = ml::relax(P, f.z, x, {h}, opt.fmax, opt.max_steps);
        last = r.pos;
        gv[slot] = r;
        hv[slot] = 1;
      }
    };
    std::thread back([&] { chain(-1, *pot2, have[1], got[1]); });
    chain(1, pot, have[0], got[0]);
    back.join();
    for (int slot = 0; slot < npts; ++slot) {
      if (!have[0][slot] && !have[1][slot]) continue;
      const int pick = !have[1][slot] ? 0 : !have[0][slot] ? 1 : (got[1][slot].energy < got[0][slot].energy ? 1 : 0);
      const auto& r = got[pick][slot];
      const double emm = mm_energy(f, ft, r.pos);
      std::vector<double> row;
      for (const auto& kk : group) {   // CAPS's OPLS form: 1/2 V1 (1 + cos) + 1/2 V2 (1 - cos 2) + 1/2 V3 (1 + cos 3) + 1/2 V4 (1 - cos 4)
        double s1 = 0, s2 = 0, s3 = 0, s4 = 0;
        for (const auto& o : occ[kk]) {
          const double p = dihedral(r.pos, o[0], o[1], o[2], o[3]);
          s1 += 0.5 * (1 + std::cos(p));
          s2 += 0.5 * (1 - std::cos(2 * p));
          s3 += 0.5 * (1 + std::cos(3 * p));
          s4 += 0.5 * (1 - std::cos(4 * p));
        }
        row.insert(row.end(), {s1, s2, s3, s4});
      }
      row.push_back(1.0);
      rows.push_back(row);
      ys.push_back(r.energy * kEvKcal - emm);
      pts.push_back({dihedral(r.pos, quad[0], quad[1], quad[2], quad[3]) * 180.0 / kPi, r.energy * kEvKcal, emm, r.converged});
    }
    {
      double a0 = 1e300, b0 = 1e300;
      for (const auto& p : pts) { a0 = std::min(a0, p.e_ml); b0 = std::min(b0, p.e_mm); }
      for (size_t i = 0; i < pts.size(); ++i) {
        pts[i].e_ml -= a0;
        pts[i].e_mm -= b0;
        pts[i].used = pts[i].e_ml <= opt.window && pts[i].e_mm <= opt.window;
        if (pts[i].used) {
          A.push_back(rows[i]);
          y.push_back(ys[i]);
          e_ml.push_back(pts[i].e_ml);
        }
      }
    }
    for (const auto& k : group) done.insert(k);
    if (A.size() < 6) {
      for (const auto& k : group)
        rep.terms.push_back({"dihedral", k, "not_fitted", {}, 0, static_cast<int>(A.size()), static_cast<int>(f.z.size()), "too few scan points without a clash"});
      continue;
    }
    // points over 10 kcal/mol above the scan's lowest weigh 0.2 (as caps.torsionfit); a light ridge on the V terms
    const double emin = *std::min_element(e_ml.begin(), e_ml.end());
    std::vector<double> w;
    for (double e : e_ml) w.push_back(e - emin > 10 ? 0.2 : 1.0);
    // classes about one bond share terms the scan cannot tell apart (their cos 3 phi terms are the same function when
    // their dihedrals differ by 120 degrees): a ridge of 0.05 then splits them evenly, the smallest V that fit the
    // profile; one class alone keeps a ridge of 0.001
    std::vector<double> ridge(A[0].size(), group.size() > 1 ? 0.05 : 0.001);
    ridge.back() = 0;
    const auto sol = lsq(A, y, w, ridge);
    double rss = 0;
    for (size_t r = 0; r < A.size(); ++r) {
      double e = -y[r];
      for (size_t k = 0; k < sol.size(); ++k) e += A[r][k] * sol[k];
      rss += e * e;
    }
    const double rms = std::sqrt(rss / static_cast<double>(A.size()));
    for (size_t g = 0; g < group.size(); ++g) {
      CompleteTerm term{"dihedral", group[g], "fitted", {sol[4 * g], sol[4 * g + 1], sol[4 * g + 2], sol[4 * g + 3]}, rms,
                        static_cast<int>(A.size()), static_cast<int>(f.z.size()), "", pts};
      if (rms > opt.torsion_rms_max) {
        term.status = "rejected";
        term.why = "a poor fit (RMS " + fmt(rms) + " kcal/mol)";
        rep.terms.push_back(term);
        say("dihedral " + join(group[g]) + ": rejected (" + term.why + ")");
        continue;
      }
      FFRule rule;
      rule.name = join(group[g]);
      rule.match = group[g];
      rule.style = "opls";
      rule.params = {six(term.params[0]), six(term.params[1]), six(term.params[2]), six(term.params[3])};
      rule.comment = "CAPS ff complete (" + label + "): " + std::to_string(term.points) + "-point relaxed scan about the " +
                     group[g][1] + "-" + group[g][2] + " bond in a " + std::to_string(term.fragment_atoms) + "-atom fragment, " +
                     std::to_string(group.size()) + " class(es) fitted together, RMS " + fmt(rms) + " kcal/mol";
      ov.dihedrals.push_back(rule);
      rep.terms.push_back(term);
      say("dihedral " + join(group[g]) + ": V = [" + fmt(term.params[0], 4) + ", " + fmt(term.params[1], 4) + ", " +
          fmt(term.params[2], 4) + ", " + fmt(term.params[3], 4) + "], RMS " + fmt(rms));
    }
    current = ff;
    merge_forcefield(current, ov);
  }

  // ------------------------------------------------------------------ the check
  FFDef merged = ff;
  merge_forcefield(merged, ov);
  ParamReport after;
  (void)parameterize(s, merged, types, "types", &after, true);
  rep.missing_after = after.missing;
  rep.complete = after.missing.empty();
  return rep;
}

}  // namespace caps
