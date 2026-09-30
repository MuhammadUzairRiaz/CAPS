// CAPS DL_POLY export (see dlpoly.hpp).
#include "caps/dlpoly.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>

#include "caps/elements.hpp"

namespace caps {

namespace {

constexpr double kDeg = 180.0 / 3.14159265358979323846;

std::string fmt(const char* f, double a) { char b[64]; std::snprintf(b, sizeof b, f, a); return b; }

// one line of a term: key, atom indices (1-based within the molecule), parameters
std::string line(const std::string& key, std::initializer_list<int> idx, std::initializer_list<double> p) {
  char b[64];
  std::string s = key;
  s.resize(std::max<size_t>(s.size(), 5), ' ');
  for (int i : idx) { std::snprintf(b, sizeof b, " %6d", i); s += b; }
  for (double v : p) { std::snprintf(b, sizeof b, " %16.8f", v); s += b; }
  return s;
}

bool near(double a, double b) { return std::fabs(a - b) < 1e-9; }

}  // namespace

std::vector<std::string> write_dlpoly(const System& s, const ForceField& ff, const std::string& dir, const DlpolyOptions& o) {
  const size_t n = s.atoms.size();
  if (ff.type_index.size() != n) throw std::invalid_argument("DL_POLY: the force field was assigned to another structure");
  if (!ff.vsites.empty()) throw std::invalid_argument("DL_POLY: virtual sites are not written");
  if (ff.sw.on || ff.hbond.on()) throw std::invalid_argument("DL_POLY: Stillinger–Weber and DREIDING hydrogen-bond terms are not written");
  if (ff.manybody.on())
    throw std::invalid_argument("DL_POLY: the " + ff.manybody.style + " potential is a LAMMPS file; DL_POLY takes its own TERSOFF / METAL tables, which CAPS does not convert. Export to LAMMPS");
  if (!ff.lj_pairs.empty()) throw std::invalid_argument("DL_POLY: explicit [ pairs ] Lennard-Jones terms are not written");
  if (ff.keep13) throw std::invalid_argument("DL_POLY: force fields where 1-3 pairs interact (special_bonds 0 1 1) are not written");
  std::vector<std::string> notes;
  std::map<std::string, size_t> dropped;   // what DL_POLY cannot express, with counts

  // molecules (connected components) in the order of their first atom; atoms keep their order within a molecule
  int nm = 0;
  const auto comp = s.molecules(&nm);
  std::vector<std::vector<uint32_t>> mols(static_cast<size_t>(nm));
  {
    std::vector<int> rank(static_cast<size_t>(nm), -1);
    int next = 0;
    for (uint32_t i = 0; i < n; ++i)
      if (rank[size_t(comp[i])] < 0) rank[size_t(comp[i])] = next++;
    for (uint32_t i = 0; i < n; ++i) mols[size_t(rank[size_t(comp[i])])].push_back(i);
  }
  std::vector<int> molof(n), local(n);
  bool reordered = false;
  {
    uint32_t k = 0;
    for (size_t m = 0; m < mols.size(); ++m)
      for (size_t a = 0; a < mols[m].size(); ++a) {
        molof[mols[m][a]] = int(m), local[mols[m][a]] = int(a) + 1;
        reordered |= mols[m][a] != k++;
      }
  }
  if (reordered) notes.push_back("CONFIG lists the atoms molecule by molecule (DL_POLY needs each molecule's atoms together)");

  // ---- term lines per molecule
  struct Terms { std::vector<std::string> bonds, angles, dihedrals, inversions; };
  std::vector<Terms> T(mols.size());
  auto M = [&](uint32_t a) -> Terms& { return T[size_t(molof[a])]; };
  auto L = [&](uint32_t a) { return local[a]; };

  for (const auto& b : ff.bonds) M(b.i).bonds.push_back(line("harm", {L(b.i), L(b.j)}, {2 * b.k, b.r0}));
  for (const auto& b : ff.bonds2) M(b.i).bonds.push_back(line("quar", {L(b.i), L(b.j)}, {2 * b.k2, b.r0, 3 * b.k3, 4 * b.k4}));
  for (const auto& b : ff.bonds_x) {
    if (b.form == 1) M(b.i).bonds.push_back(line("mors", {L(b.i), L(b.j)}, {b.a, b.c, b.b}));   // E0 [(1 − e^(−k(r−r0)))² − 1]: the same up to a constant
    else if (b.form == 3) {
      M(b.i).bonds.push_back(line("fene", {L(b.i), L(b.j)}, {b.a, b.b, 0.0}));
      if (b.c > 0) ++dropped["the WCA part of FENE bonds"];
    } else throw std::invalid_argument("DL_POLY: bond form " + std::to_string(b.form) + " has no DL_POLY equivalent");
  }
  for (const auto& u : ff.urey_bradley) { M(u.i).bonds.push_back(line("harm", {L(u.i), L(u.k)}, {2 * u.kub, u.r0})); ++dropped["Urey–Bradley terms written as harmonic bonds between the angle's end atoms"]; }
  for (const auto& a : ff.angles) M(a.i).angles.push_back(line("harm", {L(a.i), L(a.j), L(a.k)}, {2 * a.kt, a.theta0 * kDeg}));
  for (const auto& a : ff.angles2) {
    M(a.i).angles.push_back(line("quar", {L(a.i), L(a.j), L(a.k)}, {2 * a.k2, a.theta0 * kDeg, 3 * a.k3, 4 * a.k4}));
    if (a.bb_m != 0) ++dropped["class II bond–bond terms"];
    if (a.ba_n1 != 0 || a.ba_n2 != 0) ++dropped["class II bond–angle terms"];
  }
  for (const auto& a : ff.angles_x) {
    const std::initializer_list<int> ix = {L(a.i), L(a.j), L(a.k)};
    if (a.form == 1) M(a.i).angles.push_back(line("hcos", ix, {2 * a.a, a.b * kDeg}));
    else if (a.form == 2) M(a.i).angles.push_back(line("cos", ix, {a.a, 0.0, 1.0}));
    else if (a.form >= 11 && a.form <= 14) {
      const int m = a.form - 10;
      M(a.i).angles.push_back(m == 1 ? line("cos", ix, {a.a, 0.0, 1.0}) : line("cos", ix, {a.a / (m * m), 180.0, double(m)}));
    } else throw std::invalid_argument("DL_POLY: angle form " + std::to_string(a.form) + " has no DL_POLY equivalent");
  }

  // dihedrals: the 1-4 scale factors go on one term per 1-4 pair
  // each 1-4 pair's scales (its own part's when force fields with different 1-4 scalings were merged by group)
  std::map<std::pair<uint32_t, uint32_t>, std::pair<double, double>> p14;
  for (size_t k = 0; k < ff.pairs14.size(); ++k) {
    const auto& p = ff.pairs14[k];
    p14[{std::min(p[0], p[1]), std::max(p[0], p[1])}] = {k < ff.pairs14_coul.size() ? ff.pairs14_coul[k] : ff.coul14, k < ff.pairs14_lj.size() ? ff.pairs14_lj[k] : ff.lj14};
  }
  std::set<std::pair<uint32_t, uint32_t>> scaled;
  auto scale14 = [&](uint32_t i, uint32_t l) -> std::pair<double, double> {
    const auto key = std::make_pair(std::min(i, l), std::max(i, l));
    const auto it = p14.find(key);
    if (it == p14.end() || scaled.count(key)) return {0.0, 0.0};
    scaled.insert(key);
    return it->second;
  };
  // Fourier terms of one quadruple (same atom order): "cos3" when they are OPLS's three, else one "cos" line each
  auto write_torsions = [&](uint32_t i, uint32_t j, uint32_t k, uint32_t l, const std::vector<const TorsionTerm*>& ts, bool proper) {
    const auto sc = proper ? scale14(i, l) : std::make_pair(0.0, 0.0);
    double A[4] = {0, 0, 0, 0};
    bool cos3 = !ts.empty();
    for (const auto* t : ts) {
      const bool ok = (t->n == 1 && near(t->delta, 0)) || (t->n == 2 && near(std::fabs(t->delta), 3.14159265358979323846)) || (t->n == 3 && near(t->delta, 0));
      if (!ok || A[t->n] != 0) { cos3 = false; break; }
      A[t->n] = 2 * t->v;
    }
    Terms& D = M(i);
    if (cos3) {
      D.dihedrals.push_back(line("cos3", {L(i), L(j), L(k), L(l)}, {A[1], A[2], A[3], sc.first, sc.second}));
      return;
    }
    bool first = true;
    for (const auto* t : ts) {
      D.dihedrals.push_back(line("cos", {L(i), L(j), L(k), L(l)}, {t->v, t->delta * kDeg, double(t->n), first ? sc.first : 0.0, first ? sc.second : 0.0}));
      first = false;
    }
    if (ts.empty()) D.dihedrals.push_back(line("cos", {L(i), L(j), L(k), L(l)}, {0.0, 0.0, 1.0, sc.first, sc.second}));
  };
  auto group_terms = [&](const std::vector<TorsionTerm>& list, bool proper) {
    std::map<std::array<uint32_t, 4>, std::vector<const TorsionTerm*>> g;
    std::vector<std::array<uint32_t, 4>> order;
    for (const auto& t : list) {
      const std::array<uint32_t, 4> q = {t.i, t.j, t.k, t.l};
      auto it = g.find(q);
      if (it == g.end()) { order.push_back(q); g[q] = {&t}; } else it->second.push_back(&t);
    }
    for (const auto& q : order) write_torsions(q[0], q[1], q[2], q[3], g[q], proper);
  };
  group_terms(ff.dihedrals, true);
  for (const auto& d : ff.dihedrals2) {
    const auto sc = scale14(d.i, d.l);
    const double K[3] = {d.k1, d.k2, d.k3}, P[3] = {d.phi1, d.phi2, d.phi3};
    bool first = true;
    for (int m = 0; m < 3; ++m) {   // K [1 − cos(mφ − φm)] = K [1 + cos(mφ − φm − 180°)]
      if (K[m] == 0 && !(first && m == 2)) continue;
      M(d.i).dihedrals.push_back(line("cos", {L(d.i), L(d.j), L(d.k), L(d.l)}, {K[m], std::fmod(P[m] * kDeg + 180.0, 360.0), double(m + 1), first ? sc.first : 0.0, first ? sc.second : 0.0}));
      first = false;
    }
    bool cross = d.mbt_r2 != 0 && (d.mbt[0] != 0 || d.mbt[1] != 0 || d.mbt[2] != 0);
    for (int m = 0; m < 3; ++m) cross |= d.ebt_b[m] != 0 || d.ebt_c[m] != 0 || d.at_d[m] != 0 || d.at_e[m] != 0;
    cross |= d.aat_m != 0 || d.bb13_n != 0;
    if (cross) ++dropped["class II torsion cross terms (middle/end bond–torsion, angle–torsion, angle–angle–torsion, bond–bond 1-3)"];
  }
  if (ff.impropers_dlpoly_set && ff.impropers_dlpoly_ok) {   // DL-derived files: as their DL_POLY files write them
    group_terms(ff.impropers_dlpoly, false);
    if (ff.improper_written == "center2" && !ff.impropers.empty()) notes.push_back("impropers as this force field's DL_POLY files write them (centre third), which differs from its LAMMPS order (centre second)");
  } else if (ff.improper_written == "center2") {   // held i C j l (their LAMMPS order); their DL_POLY files differ
    {
      std::vector<TorsionTerm> imp = ff.impropers;
      for (auto& t : imp) t = {t.k, t.i, t.j, t.l, t.v, t.n, t.delta};
      group_terms(imp, false);
    }
    if (!ff.impropers.empty()) notes.push_back("impropers as this force field's DL_POLY files write them (centre third), which differs from its LAMMPS order (centre second)");
  } else {
    group_terms(ff.impropers, false);
  }
  for (const auto& h : ff.impropers_harmonic)
    M(h.i).dihedrals.push_back(line("harm", {L(h.i), L(h.j), L(h.k), L(h.l)}, {2 * h.k2, h.chi0 * kDeg, 0.0, 0.0, 0.0}));
  for (const auto& v : ff.inversions) {
    if (v.form != 0) throw std::invalid_argument("DL_POLY: planar (DREIDING) and UFF inversions are not written");
    M(v.c).inversions.push_back(line("harm", {L(v.c), L(v.a), L(v.b), L(v.d)}, {2 * v.kw, v.w0 * kDeg}));
  }
  for (const auto& v : ff.impropers2) {   // class II out-of-plane as DL_POLY's harmonic inversion, the centre first
    M(v.j).inversions.push_back(line("harm", {L(v.j), L(v.i), L(v.k), L(v.l)}, {2 * v.kchi, v.chi0 * kDeg}));
    if (v.m1 != 0 || v.m2 != 0 || v.m3 != 0) ++dropped["class II angle–angle terms"];
  }
  if (!ff.impropers2.empty()) notes.push_back("class II out-of-plane terms are written as DL_POLY harmonic inversions (DL_POLY's inversion angle, not the mean of the three Wilson angles)");
  if (!ff.cbt.empty()) dropped["combined bending–torsion terms"] += ff.cbt.size();
  // 1-4 pairs with no torsion term: a zero-amplitude dihedral carries their scale factors
  if (!p14.empty()) {
    const auto nb = s.neighbours();
    for (const auto& [pr, sc] : p14) {
      if (scaled.count(pr)) continue;
      bool done = false;
      for (uint32_t j : nb[pr.first]) {
        for (uint32_t k : nb[j])
          if (k != pr.first && std::find(nb[k].begin(), nb[k].end(), pr.second) != nb[k].end()) {
            write_torsions(pr.first, j, k, pr.second, {}, true);
            done = true;
            break;
          }
        if (done) break;
      }
    }
  }
  if (!ff.lj14_types.empty()) notes.push_back("DL_POLY scales the ordinary Lennard-Jones for 1-4 pairs: the force field's separate 1-4 parameters are not written");
  for (const auto& [what, count] : dropped) notes.push_back(std::to_string(count) + " " + what + (what.find("written as") == std::string::npos ? " left out (DL_POLY has no such term)" : ""));

  // ---- molecular types: consecutive identical molecules share one
  auto type_name = [&](uint32_t a) {
    std::string t = ff.type_index[a] >= 0 && size_t(ff.type_index[a]) < ff.type_names.size() ? ff.type_names[size_t(ff.type_index[a])] : element(s.atoms[a].element).symbol;
    for (char& c : t) if (c == ' ') c = '_';
    return t.substr(0, 8);
  };
  auto mass_of = [&](uint32_t a) { return a < ff.mass.size() ? ff.mass[a] : s.mass_of(s.atoms[a]); };
  auto charge_of = [&](uint32_t a) { return a < ff.charge.size() ? ff.charge[a] : s.atoms[a].charge; };
  auto body = [&](size_t m) {
    std::string t;
    char b[160];
    std::snprintf(b, sizeof b, "atoms %zu\n", mols[m].size());
    t += b;
    for (uint32_t a : mols[m]) {
      std::string nm8 = type_name(a);
      nm8.resize(std::max<size_t>(nm8.size(), 8), ' ');
      std::snprintf(b, sizeof b, "%s  %12.6f  %10.6f   1   0\n", nm8.c_str(), mass_of(a), charge_of(a));
      t += b;
    }
    auto sec = [&](const char* name, const std::vector<std::string>& v) {
      if (v.empty()) return;
      t += std::string(name) + " " + std::to_string(v.size()) + "\n";
      for (const auto& x : v) t += x + "\n";
    };
    sec("bonds", T[m].bonds);
    sec("angles", T[m].angles);
    sec("dihedrals", T[m].dihedrals);
    sec("inversions", T[m].inversions);
    return t;
  };
  struct MolType { std::string body; int count; size_t first; };
  std::vector<MolType> types;
  for (size_t m = 0; m < mols.size(); ++m) {
    std::string b = body(m);
    if (o.group_molecules && !types.empty() && types.back().body == b) ++types.back().count;
    else types.push_back({std::move(b), 1, m});
  }

  // ---- vdw: every pair of the types present
  std::map<std::string, int> tname_to_index;
  for (uint32_t a = 0; a < n; ++a) tname_to_index.emplace(type_name(a), ff.type_index[a]);
  std::vector<std::pair<std::string, int>> tlist(tname_to_index.begin(), tname_to_index.end());
  std::vector<std::string> vdw;
  const bool lj96 = ff.pair_form == "lj9-6";
  for (size_t x = 0; x < tlist.size(); ++x)
    for (size_t y = x; y < tlist.size(); ++y) {
      const int a = tlist[x].second, b = tlist[y].second;
      std::string l = tlist[x].first;
      l.resize(std::max<size_t>(l.size(), 8), ' ');
      l += "  " + tlist[y].first;
      l.resize(std::max<size_t>(l.size(), 18), ' ');
      char buf[160];
      const auto pf = ff.pair_func.find({std::min(a, b), std::max(a, b)});
      if (pf != ff.pair_func.end()) {
        if (pf->second.form == 1) std::snprintf(buf, sizeof buf, "  buck  %14.6f %12.6f %14.6f", pf->second.a, pf->second.b, pf->second.c);
        else if (pf->second.form == 2) std::snprintf(buf, sizeof buf, "  mors  %14.6f %12.6f %12.6f", pf->second.a, pf->second.c, pf->second.b);
        else throw std::invalid_argument("DL_POLY: pair form " + std::to_string(pf->second.form) + " has no DL_POLY equivalent");
      } else {
        const bool off = ff.excluded_type_pairs.count({std::min(a, b), std::max(a, b)}) > 0;
        const PairType p = off ? PairType{0, 0} : mixed_pair(ff, a, b);
        if (lj96) std::snprintf(buf, sizeof buf, "  nm    %14.6e %6d %6d %12.6f", p.eps, 9, 6, p.sigma);
        else std::snprintf(buf, sizeof buf, "  lj    %14.6e %12.6f", p.eps, p.sigma);
      }
      vdw.push_back(l + buf);
    }

  namespace fs = std::filesystem;
  fs::create_directories(dir);
  {
    std::ofstream f(fs::path(dir) / "FIELD");
    if (!f) throw std::runtime_error("cannot write " + (fs::path(dir) / "FIELD").string());
    f << o.title << " · " << ff.name << "\nunits kcal\nmolecular types " << types.size() << "\n";
    for (size_t t = 0; t < types.size(); ++t) {
      f << "molecule " << (t + 1) << "\nnummols " << types[t].count << "\n" << types[t].body << "finish\n";
    }
    f << "vdw " << vdw.size() << "\n";
    for (const auto& v : vdw) f << v << "\n";
    f << "close\n";
  }
  // CONFIG: positions about the cell centre, folded into the cell
  const bool per = s.cell.valid();
  int imcon = 0;
  if (per) {
    const auto& c = s.cell;
    const bool ortho = std::fabs(c.a[1]) + std::fabs(c.a[2]) + std::fabs(c.b[0]) + std::fabs(c.b[2]) + std::fabs(c.c[0]) + std::fabs(c.c[1]) < 1e-9;
    imcon = !ortho ? 3 : near(c.a[0], c.b[1]) && near(c.a[0], c.c[2]) ? 1 : 2;
  }
  {
    std::ofstream f(fs::path(dir) / "CONFIG");
    char b[200];
    f << o.title << "\n";
    std::snprintf(b, sizeof b, "%10d%10d%10zu\n", 0, imcon, n);
    f << b;
    if (per)
      for (const Vec3& v : {s.cell.a, s.cell.b, s.cell.c}) {
        std::snprintf(b, sizeof b, "%20.10f%20.10f%20.10f\n", v[0], v[1], v[2]);
        f << b;
      }
    const Vec3 centre = per ? s.cell.origin + (s.cell.a + s.cell.b + s.cell.c) * 0.5 : Vec3{0, 0, 0};
    size_t idx = 0;
    for (const auto& mol : mols)
      for (uint32_t a : mol) {
        Vec3 p = s.atoms[a].pos;
        if (per) p = s.cell.minimum_image(p - centre);
        std::snprintf(b, sizeof b, "%-8s%10zu\n%20.10f%20.10f%20.10f\n", type_name(a).c_str(), ++idx, p[0], p[1], p[2]);
        f << b;
      }
  }
  // CONTROL: a generic NVT run (edit to taste)
  {
    const bool charged = std::any_of(ff.charge.begin(), ff.charge.end(), [](double q) { return q != 0; });
    std::ofstream f(fs::path(dir) / "CONTROL");
    f << o.title << " · generic NVT run written by CAPS (edit before use)\n\n"
      << "temperature      " << fmt("%.2f", o.temperature) << "\n"
      << "ensemble nvt hoover 0.1\n"
      << "steps            " << o.steps << "\n"
      << "equilibration    0\n"
      << "timestep         " << fmt("%.6f", o.timestep_fs / 1000.0) << "\n"
      << "cutoff           " << fmt("%.2f", o.cutoff) << "\n"
      << "rvdw             " << fmt("%.2f", o.cutoff) << "\n"
      << (charged ? (per ? "ewald precision  1.0e-6\n" : "coul\n") : "no elec\n")
      << "print            100\nstats            100\n"
      << "trajectory       0 1000 0\n"
      << "job time         1000000\nclose time       100\n\nfinish\n";
  }
  if (per) {
    const double V = s.cell.volume();
    const double w = std::min({V / norm(cross(s.cell.b, s.cell.c)), V / norm(cross(s.cell.c, s.cell.a)), V / norm(cross(s.cell.a, s.cell.b))});
    if (o.cutoff > 0.5 * w) notes.push_back("the cut-off " + fmt("%.1f", o.cutoff) + " Å exceeds half the cell's narrowest width (" + fmt("%.1f", 0.5 * w) + " Å): shorten it in CONTROL");
  }
  notes.insert(notes.begin(), std::to_string(types.size()) + " molecular type" + (types.size() == 1 ? "" : "s") + " for " + std::to_string(mols.size()) + " molecule" +
                                  (mols.size() == 1 ? "" : "s") + " · " + std::to_string(vdw.size()) + " van der Waals pairs");
  return notes;
}

}  // namespace caps
