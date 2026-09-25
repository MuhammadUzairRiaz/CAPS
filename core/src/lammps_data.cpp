// LAMMPS data files carrying a CAPS force field: every term the evaluator knows, written with the LAMMPS style that
// computes the same energy. Kinds that mix forms (class II bonds with class I torsions, as in DL_FIELD's PCFF; CVFF
// with harmonic impropers; Morse and harmonic bonds) become hybrid styles, and the class II cross-term sections get
// "skip" lines for the types of the other sub-styles, as LAMMPS reads them.
//
// Checked term by term, energies and forces, against LAMMPS for each force-field family (bench/ff/check_data_lammps.py).
#include <cstdio>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "caps/relax.hpp"

namespace caps {

namespace {

constexpr double R2D = 57.29577951308232;

std::string num(std::initializer_list<double> v) {
  std::string r;
  char b[40];
  for (double x : v) {
    std::snprintf(b, sizeof b, " %.10g", x);
    r += b;
  }
  return r;
}

// One kind of bonded term (bonds, angles, dihedrals or impropers) with its types, sub-styles and class II sections.
struct Kind {
  Kind(std::string n, std::vector<std::string> c) : name(std::move(n)), cross(std::move(c)) {}
  std::string name;                        // Bond, Angle, Dihedral, Improper
  std::vector<std::string> cross;          // class II cross-term sections, in LAMMPS's names
  std::vector<std::string> styles;         // sub-styles in first-seen order
  struct Type {
    std::string style, coef;
    std::vector<std::string> cross;        // per cross section (class2 only)
    std::string label;
  };
  std::vector<Type> types;
  std::map<std::string, int> ids;
  std::vector<int> term_type;
  std::vector<std::vector<uint32_t>> term_atoms;

  void add(const std::string& style, const std::string& coef, std::vector<std::string> cross_coef, std::vector<uint32_t> atoms,
           const std::string& label) {
    std::string key = style + "|" + coef;
    for (const auto& c : cross_coef) key += "|" + c;
    auto [it, fresh] = ids.emplace(key, int(types.size()) + 1);
    if (fresh) {
      types.push_back({style, coef, std::move(cross_coef), label});
      if (std::find(styles.begin(), styles.end(), style) == styles.end()) styles.push_back(style);
    }
    term_type.push_back(it->second);
    term_atoms.push_back(std::move(atoms));
  }
  bool hybrid() const { return styles.size() > 1; }
  bool has(const std::string& s) const { return std::find(styles.begin(), styles.end(), s) != styles.end(); }
  std::string style_line() const {
    std::string r = hybrid() ? "hybrid" : "";
    for (const auto& s : styles) r += (r.empty() ? "" : " ") + s;
    return r;
  }
};

struct Layout {
  Kind bonds{"Bond", {}}, angles{"Angle", {"BondBond", "BondAngle"}},
      dihedrals{"Dihedral", {"MiddleBondTorsion", "EndBondTorsion", "AngleTorsion", "AngleAngleTorsion", "BondBond13"}},
      impropers{"Improper", {"AngleAngle"}};
  std::string pair_base;                   // lj/cut or lj/class2
  bool pair_hybrid = false;                // hybrid/overlay (9-6 LJ, or Buckingham / Morse pairs)
  std::set<std::string> pair_styles;       // sub-styles used in PairIJ lines
  bool periodic = true;                    // a cell: tail corrections apply (CAPS adds none without a volume)
};

Layout build(const System& s, const ForceField& ff) {
  if (!ff.lj14_types.empty())
    throw FieldError(ff.name + ": separate 1-4 Lennard-Jones parameters (CHARMM, GROMOS) have no exact LAMMPS form without switching "
                     "(lj/charmm/coul/*); LAMMPS data for them is not written");
  Layout L;
  const auto& T = ff.atom_type;
  auto lab = [&](std::initializer_list<uint32_t> a) {
    std::string r;
    for (uint32_t x : a) r += (r.empty() ? "" : "-") + T[x];
    return r;
  };
  // bonds; structure bonds without a term still go in (style zero) so LAMMPS excludes the same 1-2 / 1-3 / 1-4 pairs
  std::set<std::pair<uint32_t, uint32_t>> have;
  auto mark = [&](uint32_t i, uint32_t j) { have.insert({std::min(i, j), std::max(i, j)}); };
  for (const auto& b : ff.bonds) { L.bonds.add("harmonic", num({b.k, b.r0}), {}, {b.i, b.j}, lab({b.i, b.j})); mark(b.i, b.j); }
  for (const auto& b : ff.bonds2) { L.bonds.add("class2", num({b.r0, b.k2, b.k3, b.k4}), {}, {b.i, b.j}, lab({b.i, b.j})); mark(b.i, b.j); }
  for (const auto& b : ff.bonds_x) {
    if (b.form == 1) L.bonds.add("morse", num({b.a, b.b, b.c}), {}, {b.i, b.j}, lab({b.i, b.j}));
    else if (b.form == 2) L.bonds.add("gromos", num({b.a, b.b}), {}, {b.i, b.j}, lab({b.i, b.j}));
    else throw FieldError("bond form " + std::to_string(b.form) + " has no LAMMPS style");
    mark(b.i, b.j);
  }
  for (const auto& b : s.bonds)
    if (b.i != b.j && !have.count({std::min(b.i, b.j), std::max(b.i, b.j)})) {
      L.bonds.add("zero", "", {}, {b.i, b.j}, lab({b.i, b.j}) + " (no term)");
      mark(b.i, b.j);
    }
  // angles: Urey–Bradley terms join their angle (angle_style charmm)
  std::map<std::pair<uint32_t, uint32_t>, const UreyBradley*> ub;
  for (const auto& u : ff.urey_bradley) ub[{std::min(u.i, u.k), std::max(u.i, u.k)}] = &u;
  std::set<std::pair<uint32_t, uint32_t>> ub_used;
  for (const auto& a : ff.angles) {
    const auto key = std::make_pair(std::min(a.i, a.k), std::max(a.i, a.k));
    auto it = ub.find(key);
    if (!ub.empty()) {
      const double kub = it == ub.end() ? 0 : it->second->kub, rub = it == ub.end() ? 0 : it->second->r0;
      if (it != ub.end()) ub_used.insert(key);
      L.angles.add("charmm", num({a.kt, a.theta0 * R2D, kub, rub}), {}, {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
    } else {
      L.angles.add("harmonic", num({a.kt, a.theta0 * R2D}), {}, {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
    }
  }
  if (ub_used.size() != ub.size()) throw FieldError("a Urey–Bradley term has no angle to join (LAMMPS angle_style charmm needs one)");
  for (const auto& a : ff.angles2)
    L.angles.add("class2", num({a.theta0 * R2D, a.k2, a.k3, a.k4}), {num({a.bb_m, a.bb_r1, a.bb_r2}), num({a.ba_n1, a.ba_n2, a.ba_r1, a.ba_r2})},
                 {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
  for (const auto& a : ff.angles_x) {
    if (a.form == 1) L.angles.add("cosine/squared", num({a.a, a.b * R2D}), {}, {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
    else if (a.form == 2) L.angles.add("cosine", num({a.a}), {}, {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
    else if (a.form == 3) {
      const double s0 = std::sin(a.b), c0 = std::cos(a.b);
      const double C2 = 1 / (4 * std::max(s0 * s0, 1e-8)), C1 = -4 * C2 * c0, C0 = C2 * (2 * c0 * c0 + 1);
      L.angles.add("fourier", num({a.a, C0, C1, C2}), {}, {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
    } else if (a.form > 10 && a.form <= 16) {
      // cosine/periodic: E = (2/n²) C [1 − B (−1)ⁿ cos nθ]; C = K/2, B = (−1)ⁿ, or B = 1 for n = 1 (1 + cos θ)
      const int n = a.form - 10;
      const int B = n == 1 ? 1 : (n % 2 ? -1 : 1);
      L.angles.add("cosine/periodic", num({a.a / 2}) + " " + std::to_string(B) + " " + std::to_string(n), {}, {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
    }
    else throw FieldError("angle form " + std::to_string(a.form) + " has no LAMMPS style");
  }
  // dihedrals: the Fourier terms of one atom quadruple make one dihedral_style fourier entry
  std::map<std::array<uint32_t, 4>, std::vector<const TorsionTerm*>> quad;
  std::vector<std::array<uint32_t, 4>> order;
  for (const auto& t : ff.dihedrals) {
    const std::array<uint32_t, 4> k{t.i, t.j, t.k, t.l};
    auto& v = quad[k];
    if (v.empty()) order.push_back(k);
    v.push_back(&t);
  }
  for (const auto& k : order) {
    std::string c = " " + std::to_string(quad[k].size());
    for (const auto* t : quad[k]) {
      if (t->n < 0) throw FieldError("a torsion with multiplicity < 0 has no LAMMPS fourier form");
      c += num({t->v}) + " " + std::to_string(t->n) + num({t->delta * R2D});
    }
    L.dihedrals.add("fourier", c, {}, {k[0], k[1], k[2], k[3]}, lab({k[0], k[1], k[2], k[3]}));
  }
  for (const auto& d : ff.dihedrals2)
    L.dihedrals.add("class2", num({d.k1, d.phi1 * R2D, d.k2, d.phi2 * R2D, d.k3, d.phi3 * R2D}),
                    {num({d.mbt[0], d.mbt[1], d.mbt[2], d.mbt_r2}),
                     num({d.ebt_b[0], d.ebt_b[1], d.ebt_b[2], d.ebt_c[0], d.ebt_c[1], d.ebt_c[2], d.ebt_r1, d.ebt_r3}),
                     num({d.at_d[0], d.at_d[1], d.at_d[2], d.at_e[0], d.at_e[1], d.at_e[2], d.at_theta1 * R2D, d.at_theta2 * R2D}),
                     num({d.aat_m, d.aat_theta1 * R2D, d.aat_theta2 * R2D}), num({d.bb13_n, d.bb13_r1, d.bb13_r3})},
                    {d.i, d.j, d.k, d.l}, lab({d.i, d.j, d.k, d.l}));
  // impropers
  for (const auto& t : ff.impropers) {
    const double c = std::cos(t.delta);
    if (std::fabs(std::fabs(c) - 1) > 1e-6)
      throw FieldError("an improper torsion with phase " + std::to_string(t.delta * R2D) + "° has no LAMMPS cvff form (phase 0 or 180° only)");
    L.impropers.add("cvff", num({t.v}) + (c > 0 ? " 1 " : " -1 ") + std::to_string(t.n), {}, {t.i, t.j, t.k, t.l}, lab({t.i, t.j, t.k, t.l}));
  }
  for (const auto& t : ff.impropers_harmonic)
    L.impropers.add("harmonic", num({t.k2, t.chi0 * R2D}), {}, {t.i, t.j, t.k, t.l}, lab({t.i, t.j, t.k, t.l}));
  for (const auto& t : ff.impropers2)
    L.impropers.add("class2", num({t.kchi, t.chi0 * R2D}), {num({t.m1, t.m2, t.m3, t.theta1 * R2D, t.theta2 * R2D, t.theta3 * R2D})},
                    {t.i, t.j, t.k, t.l}, lab({t.i, t.j, t.k, t.l}));
  for (const auto& v : ff.inversions) {
    if (v.form == 0) {
      L.impropers.add("inversion/harmonic", num({v.kw, v.w0 * R2D}), {}, {v.c, v.a, v.b, v.d}, lab({v.c, v.a, v.b, v.d}));
    } else if (v.form == 1) {
      if (std::fabs(v.w0) > 1e-12) throw FieldError("a planar inversion with omega0 ≠ 0 has no exact LAMMPS umbrella form");
      // LAMMPS umbrella: the centre first, the angle between the last bond and the plane of the other two;
      // the three permutations at K/3 give CAPS's average over the three bonds
      const uint32_t o[3] = {v.a, v.b, v.d};
      for (int p = 0; p < 3; ++p)
        L.impropers.add("umbrella", num({v.kw / 3, 0.0}), {}, {v.c, o[(p + 1) % 3], o[(p + 2) % 3], o[p]}, lab({v.c, v.a, v.b, v.d}));
    } else if (v.form == 2) {
      // improper fourier: the centre first, ω between the I-L axis and the I-J-K plane; "all" sums the three
      // permutations, each at K/3
      const double C2 = 1, C1 = -4 * std::cos(v.w0), C0 = -(C1 * std::cos(v.w0) + C2 * std::cos(2 * v.w0));
      L.impropers.add("fourier", num({v.kw / 3, C0, C1, C2}) + " 1", {}, {v.c, v.a, v.b, v.d}, lab({v.c, v.a, v.b, v.d}));
    } else {
      throw FieldError("inversion form " + std::to_string(v.form) + " has no LAMMPS style");
    }
  }
  // pairs: the sub-styles the i-j pairs use (LAMMPS rejects a hybrid sub-style that no pair uses)
  L.pair_base = ff.pair_form == "lj9-6" ? "lj/class2" : "lj/cut";
  const int nt = int(ff.type_names.size());
  for (int a2 = 0; a2 < nt; ++a2)
    for (int b2 = a2; b2 < nt; ++b2) {
      auto it = ff.pair_func.find({a2, b2});
      L.pair_styles.insert(it == ff.pair_func.end() ? L.pair_base : it->second.form == 1 ? "buck" : it->second.form == 2 ? "morse" : "?");
    }
  if (L.pair_styles.count("?")) throw FieldError("a pair form has no LAMMPS style");
  L.pair_hybrid = ff.pair_form == "lj9-6" || !ff.pair_func.empty();
  L.periodic = s.cell.valid();
  return L;
}

// PME is used (and written) only for periodic cells, as the evaluator does.
bool pme(const EnergyOptions& e, const Layout& L) { return e.electrostatics == EnergyOptions::Electrostatics::PME && L.periodic; }

// The LAMMPS commands (after units / atom_style) that reproduce CAPS's energy with the data file.
std::vector<std::string> style_lines(const Layout& L, const ForceField& ff, const EnergyOptions& e) {
  std::vector<std::string> r;
  char b[400];
  if (!L.pair_hybrid) {
    if (e.coulomb && pme(e, L)) std::snprintf(b, sizeof b, "pair_style lj/cut/coul/long %.6g", e.cutoff);
    else if (e.coulomb) std::snprintf(b, sizeof b, "pair_style lj/cut/coul/dsf %.6g %.6g", e.dsf_alpha, e.cutoff);
    else std::snprintf(b, sizeof b, "pair_style lj/cut %.6g", e.cutoff);
    r.push_back(b);
  } else {
    std::string p = "pair_style hybrid/overlay";
    for (const auto& st : L.pair_styles) {
      std::snprintf(b, sizeof b, " %s %.6g", st.c_str(), e.cutoff);
      p += b;
    }
    if (e.coulomb && pme(e, L)) {
      std::snprintf(b, sizeof b, " coul/long %.6g", e.cutoff);
      p += b;
    } else if (e.coulomb) {
      std::snprintf(b, sizeof b, " coul/dsf %.6g %.6g", e.dsf_alpha, e.cutoff);
      p += b;
    }
    r.push_back(p);
  }
  // CAPS: with tail corrections the potentials are truncated at the cut-off (plus the tail when there is a cell);
  // without, they are shifted to zero there
  if (!e.tail) r.push_back("pair_modify shift yes");
  else if (L.periodic) r.push_back("pair_modify tail yes");
  for (const Kind* k : {&L.bonds, &L.angles, &L.dihedrals, &L.impropers}) {
    if (k->types.empty()) continue;
    std::string nm = k->name;
    for (auto& c : nm) c = char(std::tolower(static_cast<unsigned char>(c)));
    r.push_back(nm + "_style " + k->style_line());
  }
  std::snprintf(b, sizeof b, "special_bonds lj 0 0 %.10g coul 0 0 %.10g", ff.lj14, ff.coul14);
  r.push_back(b);
  if (e.coulomb && pme(e, L)) {
    // CAPS's PME with its own β; LAMMPS's Ewald sum to the same accuracy reaches the same total electrostatics
    std::snprintf(b, sizeof b, "kspace_style ewald %.3g", std::max(1e-12, e.ewald_rtol * 0.01));
    r.push_back(b);
  }
  return r;
}

// Commands that must follow read_data (hybrid pair coefficients the data file cannot hold).
std::vector<std::string> after_read(const Layout& L, const EnergyOptions& e) {
  std::vector<std::string> r;
  if (L.pair_hybrid && e.coulomb) r.push_back(pme(e, L) ? "pair_coeff * * coul/long" : "pair_coeff * * coul/dsf");
  return r;
}

}  // namespace

void write_lammps_data_ff(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& path) {
  const Layout L = build(s, ff);
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  char buf[512];
  const std::vector<const Kind*> kinds = {&L.bonds, &L.angles, &L.dihedrals, &L.impropers};

  out << "CAPS 0.1 · " << (s.title.empty() ? "structure" : s.title) << " · " << ff.name << " · units real · atom_style full\n\n";
  out << s.atoms.size() << " atoms\n";
  const char* plural[] = {"bonds", "angles", "dihedrals", "impropers"};
  for (int k = 0; k < 4; ++k) out << kinds[size_t(k)]->term_type.size() << " " << plural[k] << "\n";
  out << "\n" << ff.type_names.size() << " atom types\n";
  const char* tnames[] = {"bond", "angle", "dihedral", "improper"};
  for (int k = 0; k < 4; ++k) out << kinds[size_t(k)]->types.size() << " " << tnames[k] << " types\n";
  out << "\n";
  const Cell& c = s.cell;
  Vec3 lo = c.origin, a = c.a, b = c.b, cc = c.c;
  if (!c.valid()) { lo = {-50, -50, -50}; a = {100, 0, 0}; b = {0, 100, 0}; cc = {0, 0, 100}; }
  std::snprintf(buf, sizeof buf, "%.8f %.8f xlo xhi\n%.8f %.8f ylo yhi\n%.8f %.8f zlo zhi\n", lo[0], lo[0] + a[0], lo[1], lo[1] + b[1], lo[2],
                lo[2] + cc[2]);
  out << buf;
  if (std::fabs(b[0]) + std::fabs(cc[0]) + std::fabs(cc[1]) > 0) {
    std::snprintf(buf, sizeof buf, "%.8f %.8f %.8f xy xz yz\n", b[0], cc[0], cc[1]);
    out << buf;
  }
  out << "\n# LAMMPS commands for these coefficients (units real, atom_style full; CAPS evaluates Coulomb as damped shifted force):\n";
  for (const auto& l : style_lines(L, ff, e)) out << "#   " << l << "\n";
  out << "#   read_data <this file>\n";
  for (const auto& l : after_read(L, e)) out << "#   " << l << "\n";
  out << "# Pair coefficients: every i-j pair, " << ff.mixing << " mixing applied (nothing is left to LAMMPS's mixing).\n";

  out << "\nMasses\n\n";
  for (size_t t = 0; t < ff.type_names.size(); ++t) {
    double m = 0;
    for (size_t i = 0; i < s.atoms.size(); ++i)
      if (ff.type_index[i] == int(t)) { m = ff.mass[i]; break; }
    std::snprintf(buf, sizeof buf, "%zu %.6f  # %s\n", t + 1, m, ff.type_names[t].c_str());
    out << buf;
  }
  // every i-j pair, mixed by the force field's rule (and its explicit pairs): nothing is left to LAMMPS's mixing
  out << "\nPairIJ Coeffs  # " << (L.pair_hybrid ? std::string("hybrid/overlay") : e.coulomb ? std::string(pme(e, L) ? "lj/cut/coul/long" : "lj/cut/coul/dsf") : L.pair_base) << "\n\n";
  for (size_t a2 = 0; a2 < ff.type_names.size(); ++a2)
    for (size_t b2 = a2; b2 < ff.type_names.size(); ++b2) {
      auto it = ff.pair_func.find({int(a2), int(b2)});
      std::string coef, style = L.pair_base;
      if (it != ff.pair_func.end()) {
        style = it->second.form == 1 ? "buck" : "morse";
        coef = num({it->second.a, it->second.b, it->second.c});
      } else {
        const PairType pt = mixed_pair(ff, int(a2), int(b2));
        coef = num({pt.eps, pt.sigma});
      }
      out << a2 + 1 << " " << b2 + 1 << (L.pair_hybrid ? " " + style : "") << coef << "  # " << ff.type_names[a2] << " " << ff.type_names[b2] << "\n";
    }
  for (const Kind* k : kinds) {
    if (k->types.empty()) continue;
    out << "\n" << k->name << " Coeffs  # " << (k->hybrid() ? std::string("hybrid") : k->style_line()) << "\n\n";
    for (size_t t = 0; t < k->types.size(); ++t)
      out << t + 1 << (k->hybrid() ? " " + k->types[t].style : "") << k->types[t].coef << "  # " << k->types[t].label << "\n";
    if (!k->has("class2")) continue;
    for (size_t sec = 0; sec < k->cross.size(); ++sec) {
      out << "\n" << k->cross[sec] << " Coeffs\n\n";
      for (size_t t = 0; t < k->types.size(); ++t) {
        const auto& ty = k->types[t];
        if (ty.style == "class2") out << t + 1 << (k->hybrid() ? " class2" : "") << ty.cross[sec] << "  # " << ty.label << "\n";
        else out << t + 1 << " skip  # " << ty.label << " (" << ty.style << ")\n";
      }
    }
  }
  out << "\nAtoms  # full\n\n";
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const auto& at = s.atoms[i];
    const Vec3 fr = c.valid() ? c.to_fractional(at.pos) : Vec3{0, 0, 0};
    int im[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) im[k] = c.valid() && c.periodic[k] ? int(std::floor(fr[k])) : 0;
    const Vec3 w = at.pos - (c.a * im[0] + c.b * im[1] + c.c * im[2]);
    std::snprintf(buf, sizeof buf, "%zu %lld %d %.8f %.10f %.10f %.10f %d %d %d\n", i + 1, static_cast<long long>(at.mol), ff.type_index[i] + 1,
                  ff.charge[i], w[0], w[1], w[2], im[0], im[1], im[2]);
    out << buf;
  }
  if (s.velocities.size() == s.atoms.size()) {
    out << "\nVelocities\n\n";
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      std::snprintf(buf, sizeof buf, "%zu %.8f %.8f %.8f\n", i + 1, s.velocities[i][0], s.velocities[i][1], s.velocities[i][2]);
      out << buf;
    }
  }
  const char* sections[] = {"Bonds", "Angles", "Dihedrals", "Impropers"};
  for (int k = 0; k < 4; ++k) {
    const Kind* kd = kinds[size_t(k)];
    if (kd->term_type.empty()) continue;
    out << "\n" << sections[k] << "\n\n";
    for (size_t q = 0; q < kd->term_type.size(); ++q) {
      out << q + 1 << " " << kd->term_type[q];
      for (uint32_t x : kd->term_atoms[q]) out << " " << x + 1;
      out << "\n";
    }
  }
}

void write_lammps_input(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& data_path, const std::string& path,
                        int64_t held_mol) {
  const Layout L = build(s, ff);
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << "# LAMMPS input for " << data_path << " (" << ff.name << "), written by CAPS: reproduces CAPS's energy and forces\n";
  // a structure without a cell sits in a 100 Å box (as in the data file): periodic, but too large for images to interact
  out << "units real\natom_style full\nboundary p p p\n";
  for (const auto& l : style_lines(L, ff, e)) out << l << "\n";
  out << "read_data " << data_path << "\n";
  for (const auto& l : after_read(L, e)) out << l << "\n";
  char b[160];
  std::snprintf(b, sizeof b, "neighbor %.3g bin\ncomm_modify cutoff %.3g\n", e.skin, e.cutoff + e.skin + 2.0);
  out << b;
  if (held_mol > 0)
    out << "# molecule " << held_mol << " (the surface or filler) held in place, as in CAPS: no velocity, no force\n"
        << "group held molecule " << held_mol << "\n"
        << "velocity held set 0.0 0.0 0.0\n"
        << "fix held_in_place held setforce 0.0 0.0 0.0\n";
  out << "thermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong press\n"
         "thermo_modify format float %.10f\n"
         "run 0\n";
}

}  // namespace caps
