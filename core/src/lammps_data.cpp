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

#include "caps/io.hpp"
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
  bool sdk = false;                        // SDK / SPICA pairs (lj/sdk): no tail correction in LAMMPS
  bool gromacs = false;                    // MARTINI: lj/gromacs(/coul/gromacs), one style for every pair
  bool cos2 = false;                       // cosine/squared pairs (Cooke–Deserno): no shift, no tail
  std::vector<std::string> sw_types;       // per atom type: its Stillinger–Weber element name, or NULL (pair_style sw)
};

Layout build(const System& s, const ForceField& ff) {
  if (!ff.vsites.empty())
    throw FieldError(ff.name + ": virtual sites (Martini 3's tryptophan, ...) have no LAMMPS form; export to GROMACS instead");
  {   // exclusions LAMMPS can make: bonded 1-2 pairs, and 1-3 / 1-4 pairs when their scaling is 0
    const auto nb = s.neighbours();
    std::set<std::pair<uint32_t, uint32_t>> can;
    for (const auto& b : s.bonds) can.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
    if (!ff.keep13)
      for (uint32_t j = 0; j < nb.size(); ++j)
        for (uint32_t i : nb[j])
          for (uint32_t k : nb[j])
            if (i < k) can.insert({i, k});
    for (uint32_t i = 0; i < ff.excluded.size(); ++i)
      for (uint32_t j : ff.excluded[i])
        if (i < j && !can.count({i, j})) {
          bool is14 = false;
          for (const auto& p : ff.pairs14) is14 = is14 || (std::min(p[0], p[1]) == i && std::max(p[0], p[1]) == j);
          if (!is14)
            throw FieldError(ff.name + ": explicit exclusions between atoms that are not bonded (Martini 3's aromatic side chains, ...) have no "
                             "LAMMPS form; export to GROMACS instead");
        }
  }
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
    else if (b.form == 3) L.bonds.add("fene", num({b.a, b.b, b.c, b.d}), {}, {b.i, b.j}, lab({b.i, b.j}));
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
    else if (a.form == 4) L.angles.add("sdk", num({a.a, a.b * R2D}), {}, {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
    else if (a.form == 5) L.angles.add("cosine/squared/restricted", num({a.a, a.b * R2D}), {}, {a.i, a.j, a.k}, lab({a.i, a.j, a.k}));
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
      const int f = it == ff.pair_func.end() ? 0 : it->second.form;
      L.pair_styles.insert(f == 0 ? L.pair_base : f == 1 ? "buck" : f == 2 ? "morse" : f >= kPairSdk96 && f <= kPairSdk125 ? "lj/sdk" : f == kPairGromacs ? "lj/gromacs" : f == kPairCos2 || f == kPairCos2Wca ? "cosine/squared" : "?");
    }
  if (L.pair_styles.count("?")) throw FieldError("a pair form has no LAMMPS style");
  L.pair_hybrid = ff.pair_form == "lj9-6" || !ff.pair_func.empty();
  L.sdk = L.pair_styles.count("lj/sdk") > 0;
  L.cos2 = L.pair_styles.count("cosine/squared") > 0;
  if (L.pair_styles.count("lj/gromacs")) {
    if (L.pair_styles.size() > 1) throw FieldError("lj/gromacs with other pair forms has no LAMMPS style");
    L.gromacs = true;
    L.pair_hybrid = false;
    L.pair_base = "lj/gromacs";
  }
  if (ff.sw.on) {   // Stillinger–Weber overlays the pair terms (its types' Lennard-Jones is zero)
    L.sw_types.assign(size_t(nt), "NULL");
    for (size_t i = 0; i < ff.type_index.size(); ++i)
      if (i < ff.sw.atom.size() && ff.sw.atom[i]) L.sw_types[size_t(ff.type_index[i])] = ff.type_names[size_t(ff.type_index[i])];
    L.pair_hybrid = true;
  }
  L.periodic = s.cell.valid();
  return L;
}

// PME is used (and written) only for periodic cells, as the evaluator does.
bool pme(const EnergyOptions& e, const Layout& L) { return e.electrostatics == EnergyOptions::Electrostatics::PME && L.periodic && !L.gromacs; }

// The LAMMPS commands (after units / atom_style) that reproduce CAPS's energy with the data file.
std::vector<std::string> style_lines(const Layout& L, const ForceField& ff, const EnergyOptions& e) {
  std::vector<std::string> r;
  char b[400];
  if (L.gromacs) {   // MARTINI: the GROMACS switch for LJ (and Coulomb), inner and outer radii
    if (e.coulomb && ff.coul_gromacs)
      std::snprintf(b, sizeof b, "pair_style lj/gromacs/coul/gromacs %.6g %.6g %.6g %.6g", ff.lj_inner, e.cutoff, std::max(ff.coul_inner, 1e-6), e.cutoff);
    else if (e.coulomb) throw FieldError("lj/gromacs needs the GROMACS Coulomb form");
    else std::snprintf(b, sizeof b, "pair_style lj/gromacs %.6g %.6g", ff.lj_inner, e.cutoff);
    r.push_back(b);
  } else if (!L.pair_hybrid) {
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
    if (!L.sw_types.empty()) p += " sw";
    r.push_back(p);
  }
  // CAPS: with tail corrections the potentials are truncated at the cut-off (plus the tail when there is a cell);
  // without, they are shifted to zero there
  // lj/gromacs is zero at the cut-off by itself; lj/sdk has no tail correction (CAPS adds none for it either)
  if (L.gromacs || (L.cos2 && L.pair_styles.size() == 1)) {
  } else if (!e.tail) r.push_back("pair_modify shift yes");
  else if (L.periodic && L.sdk) {
    if (L.pair_styles.size() > 1) throw FieldError("SDK pairs with other Lennard-Jones pairs and tail corrections have no LAMMPS form (lj/sdk has no tail)");
  } else if (L.periodic) r.push_back("pair_modify tail yes");
  if (ff.dielectric != 1) {
    std::snprintf(b, sizeof b, "dielectric %.10g", ff.dielectric);
    r.push_back(b);
  }
  for (const Kind* k : {&L.bonds, &L.angles, &L.dihedrals, &L.impropers}) {
    if (k->types.empty()) continue;
    std::string nm = k->name;
    for (auto& c : nm) c = char(std::tolower(static_cast<unsigned char>(c)));
    r.push_back(nm + "_style " + k->style_line());
  }
  std::snprintf(b, sizeof b, "special_bonds lj 0 %d %.10g coul 0 %d %.10g", ff.keep13 ? 1 : 0, ff.lj14, ff.keep13 ? 1 : 0, ff.coul14);
  r.push_back(b);
  if (e.coulomb && pme(e, L)) {
    // CAPS's PME with its own β; LAMMPS's Ewald sum to the same accuracy reaches the same total electrostatics
    std::snprintf(b, sizeof b, "kspace_style ewald %.3g", std::max(1e-12, e.ewald_rtol * 0.01));
    r.push_back(b);
  }
  return r;
}

std::string clean_title(const std::string& t, const std::string& ffname) { return export_title(t, ffname); }

// Every i-j pair, mixed by the force field's rule (and its explicit pairs): "i j [style] coefficients  # A B".
std::vector<std::string> pair_lines(const Layout& L, const ForceField& ff) {
  std::vector<std::string> r;
  for (size_t a2 = 0; a2 < ff.type_names.size(); ++a2)
    for (size_t b2 = a2; b2 < ff.type_names.size(); ++b2) {
      auto it = ff.pair_func.find({int(a2), int(b2)});
      std::string coef, style = L.pair_base;
      const int f = it != ff.pair_func.end() ? it->second.form : 0;
      if (f >= kPairSdk96 && f <= kPairSdk125) {
        static const char* nm[] = {"lj9_6", "lj12_4", "lj12_6", "lj12_5"};
        style = "lj/sdk";
        coef = std::string(" ") + nm[f - kPairSdk96] + num({it->second.a, it->second.b});
      } else if (f == kPairGromacs) {
        coef = num({it->second.a, it->second.b});
      } else if (f == kPairCos2 || f == kPairCos2Wca) {
        style = "cosine/squared";
        coef = num({it->second.a, it->second.b, it->second.c}) + (f == kPairCos2Wca ? " wca" : "");
      } else if (it != ff.pair_func.end()) {
        style = it->second.form == 1 ? "buck" : "morse";
        coef = num({it->second.a, it->second.b, it->second.c});
      } else {
        const PairType pt = mixed_pair(ff, int(a2), int(b2));
        coef = num({pt.eps, pt.sigma});
      }
      r.push_back(std::to_string(a2 + 1) + " " + std::to_string(b2 + 1) + (L.pair_hybrid ? " " + style : "") + coef + "  # " + ff.type_names[a2] + " " +
                  ff.type_names[b2]);
    }
  return r;
}

// The Stillinger–Weber parameter file next to a data file: STEM.sw.
std::string sw_path(const std::string& data_path) {
  const size_t slash = data_path.find_last_of('/'), dot = data_path.find_last_of('.');
  return (dot != std::string::npos && (slash == std::string::npos || dot > slash) ? data_path.substr(0, dot) : data_path) + ".sw";
}

// Commands that must follow read_data (hybrid pair coefficients the data file cannot hold).
std::vector<std::string> after_read(const Layout& L, const EnergyOptions& e, const std::string& data_path,
                                    const std::set<std::pair<int, int>>& ff_excl = {}) {
  std::vector<std::string> r;
  if (L.pair_hybrid && e.coulomb) r.push_back(pme(e, L) ? "pair_coeff * * coul/long" : "pair_coeff * * coul/dsf");
  for (const auto& [a, b] : ff_excl) r.push_back("neigh_modify exclude type " + std::to_string(a + 1) + " " + std::to_string(b + 1));
  if (!L.sw_types.empty()) {
    std::string l = "pair_coeff * * sw " + sw_path(data_path);
    for (const auto& t : L.sw_types) l += " " + t;
    r.push_back(l);
  }
  return r;
}

// LAMMPS's Stillinger–Weber file: one entry per element triplet, CAPS's single parameter set in each.
void write_sw_file(const Layout& L, const ForceField& ff, const std::string& path) {
  std::set<std::string> el(L.sw_types.begin(), L.sw_types.end());
  el.erase("NULL");
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << "# Stillinger-Weber parameters written by CAPS for " << ff.name << " (units real: kcal/mol, A)\n"
      << "# el1 el2 el3  epsilon sigma a lambda gamma costheta0 A B p q tol\n";
  const auto& S = ff.sw;
  char b[400];
  for (const auto& i : el)
    for (const auto& j : el)
      for (const auto& k : el) {
        std::snprintf(b, sizeof b, "%s %s %s  %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g 0.0\n", i.c_str(), j.c_str(), k.c_str(), S.eps,
                      S.sigma, S.a, S.lambda, S.gamma, S.cos0, S.A, S.B, S.p, S.q);
        out << b;
      }
}

}  // namespace

// A title read back from a CAPS file carries the old header line: keep the description only.
std::string export_title(std::string t, const std::string& ffname) {
  auto erase_all = [&](const std::string& x) {
    if (x.empty()) return;
    for (size_t k; (k = t.find(x)) != std::string::npos;) t.erase(k, x.size());
  };
  for (const char* junk : {"CAPS 0.1 · ", "CAPS · ", " · atom_style full", " · units real"}) erase_all(junk);
  erase_all(" · " + ffname);
  while (!t.empty() && (t.back() == ' ' || t.back() == '\n' || t.back() == '\r')) t.pop_back();
  while (!t.empty() && t.front() == ' ') t.erase(t.begin());
  return t.empty() ? "structure" : t;
}

std::string write_lammps_data_or_structure(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& path) {
  try {
    write_lammps_data_ff(s, ff, e, path);
    return "";
  } catch (const FieldError& x) {
    write_lammps_data(s, path);
    return std::string("structure only, no coefficients (") + x.what() + ")";
  }
}

void write_lammps_data_ff(const System& s, const ForceField& ff, const EnergyOptions& e0, const std::string& path, bool pair_coeffs) {
  EnergyOptions e = e0;
  if (ff.cutoff > 0) e.cutoff = ff.cutoff;   // the model's own cut-off (MARTINI)
  // no charges, no Coulomb term (LAMMPS refuses an Ewald sum on an uncharged system; the energy is the same)
  if (std::all_of(ff.charge.begin(), ff.charge.end(), [](double q) { return q == 0; })) e.coulomb = false;
  if (ff.lj_shift) e.tail = false;   // Martini 3: shifted at the cut-off
  if (ff.coul_rf && e.coulomb)
    throw FieldError(ff.name + ": reaction-field Coulomb (Martini 3) has no LAMMPS pair style; export to GROMACS instead");
  const Layout L = build(s, ff);
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  if (!L.sw_types.empty()) write_sw_file(L, ff, sw_path(path));
  char buf[512];
  const std::vector<const Kind*> kinds = {&L.bonds, &L.angles, &L.dihedrals, &L.impropers};

  // the LAMMPS data format: one title line, the counts, the box, then the sections (atom_style full, units real)
  out << "CAPS · " << clean_title(s.title, ff.name) << " · " << ff.name << "\n\n";
  out << s.atoms.size() << " atoms\n";
  const char* plural[] = {"bonds", "angles", "dihedrals", "impropers"};
  for (int k = 0; k < 4; ++k)
    if (!kinds[size_t(k)]->term_type.empty()) out << kinds[size_t(k)]->term_type.size() << " " << plural[k] << "\n";
  out << "\n" << ff.type_names.size() << " atom types\n";
  const char* tnames[] = {"bond", "angle", "dihedral", "improper"};
  for (int k = 0; k < 4; ++k)
    if (!kinds[size_t(k)]->types.empty()) out << kinds[size_t(k)]->types.size() << " " << tnames[k] << " types\n";
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
  out << "\nMasses\n\n";
  for (size_t t = 0; t < ff.type_names.size(); ++t) {
    double m = 0;
    for (size_t i = 0; i < s.atoms.size(); ++i)
      if (ff.type_index[i] == int(t)) { m = ff.mass[i]; break; }
    std::snprintf(buf, sizeof buf, "%zu %.6f  # %s\n", t + 1, m, ff.type_names[t].c_str());
    out << buf;
  }
  // every i-j pair, mixed by the force field's rule (and its explicit pairs): nothing is left to LAMMPS's mixing
  if (pair_coeffs) {
    out << "\nPairIJ Coeffs  # " << (L.pair_hybrid ? std::string("hybrid/overlay") : L.gromacs ? std::string(e.coulomb ? "lj/gromacs/coul/gromacs" : "lj/gromacs")
                                      : e.coulomb ? std::string(pme(e, L) ? "lj/cut/coul/long" : "lj/cut/coul/dsf") : L.pair_base) << "\n\n";
    for (const auto& l : pair_lines(L, ff)) out << l << "\n";
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

void write_lammps_input(const System& s, const ForceField& ff, const EnergyOptions& e0, const std::string& data_path, const std::string& path,
                        int64_t held_mol, bool pair_coeffs, const LammpsRun& run) {
  EnergyOptions e = e0;
  if (ff.cutoff > 0) e.cutoff = ff.cutoff;   // the model's own cut-off (MARTINI)
  // no charges, no Coulomb term (LAMMPS refuses an Ewald sum on an uncharged system; the energy is the same)
  if (std::all_of(ff.charge.begin(), ff.charge.end(), [](double q) { return q == 0; })) e.coulomb = false;
  if (ff.lj_shift) e.tail = false;   // Martini 3: shifted at the cut-off
  if (ff.coul_rf && e.coulomb)
    throw FieldError(ff.name + ": reaction-field Coulomb (Martini 3) has no LAMMPS pair style; export to GROMACS instead");
  const Layout L = build(s, ff);
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  char b[400];
  out << "# LAMMPS input written by CAPS: " << clean_title(s.title, ff.name) << " · " << ff.name << "\n";
  out << "# the same force field and cut-offs CAPS uses (energies and forces checked against LAMMPS: bench/ff/check_data_lammps.py)\n\n";
  // a structure without a cell sits in a 100 Å box (as in the data file): periodic, but too large for images to interact
  out << "units           real\natom_style      full\nboundary        p p p\n\n";
  auto aligned = [](const std::string& l) {   // "keyword       arguments", as the rest of the script
    const size_t sp = l.find(' ');
    if (sp == std::string::npos || sp >= 16) return l;
    return l.substr(0, sp) + std::string(16 - sp, ' ') + l.substr(sp + 1);
  };
  for (const auto& l : style_lines(L, ff, e)) out << aligned(l) << "\n";
  out << "\nread_data       " << data_path << "\n";
  if (pair_coeffs) {
    out << "\n# pair coefficients: every type pair, " << ff.mixing << " mixing applied by CAPS (nothing left to LAMMPS's mixing)\n";
    for (const auto& l : pair_lines(L, ff)) out << "pair_coeff      " << l << "\n";
  }
  for (const auto& l : after_read(L, e, data_path, ff.excluded_type_pairs)) out << aligned(l) << "\n";
  std::snprintf(b, sizeof b, "\nneighbor        %.3g bin\nneigh_modify    delay 0 every 1 check yes\ncomm_modify     cutoff %.3g\n", e.skin, e.cutoff + e.skin + 2.0);
  out << b;
  if (held_mol > 0)
    out << "\n# molecule " << held_mol << " (the surface or filler) held in place, as in CAPS: no velocity, no force\n"
        << "group           held molecule " << held_mol << "\n"
        << "velocity        held set 0.0 0.0 0.0\n"
        << "fix             held_in_place held setforce 0.0 0.0 0.0\n";
  using K = LammpsRun::Kind;
  switch (run.kind) {
    case K::Check:
      out << "\nthermo_style custom step pe ebond eangle edihed eimp evdwl ecoul elong press\n"
             "thermo_modify format float %.10f\n"
             "run 0\n";
      return;
    case K::None:
      return;
    default:
      break;
  }
  const std::string mobile = held_mol > 0 ? "mobile" : "all";
  if (held_mol > 0) out << "group           mobile subtract all held\n";
  std::snprintf(b, sizeof b, "\nthermo          %d\nthermo_style    custom step temp press pe ke etotal density vol\n", std::max(1, run.thermo_every));
  out << b;
  if (run.minimize_first || run.kind == K::Minimize)
    out << "\n# 1. energy minimisation\nmin_style       cg\nminimize        1.0e-4 1.0e-6 5000 50000\nreset_timestep  0\n";
  if (run.kind == K::Minimize) {
    out << "\nwrite_data      minimized.data\n";
    return;
  }
  const bool npt = run.kind == K::NPT;
  out << "\n# 2. " << (npt ? "NPT" : "NVT") << " molecular dynamics (Nosé–Hoover)\n";
  std::snprintf(b, sizeof b, "velocity        %s create %.6g %llu mom yes rot yes dist gaussian\ntimestep        %.6g\n", mobile.c_str(), run.temperature,
                static_cast<unsigned long long>(run.seed), run.dt);
  out << b;
  if (npt)
    std::snprintf(b, sizeof b, "fix             integrate %s npt temp %.6g %.6g %.6g iso %.6g %.6g %.6g\n", mobile.c_str(), run.temperature, run.temperature,
                  run.tdamp, run.pressure, run.pressure, run.pdamp);
  else
    std::snprintf(b, sizeof b, "fix             integrate %s nvt temp %.6g %.6g %.6g\n", mobile.c_str(), run.temperature, run.temperature, run.tdamp);
  out << b;
  std::snprintf(b, sizeof b, "dump            traj all custom %d traj.lammpstrj id mol type q xu yu zu\ndump_modify     traj sort id\nrun             %lld\n",
                std::max(1, run.dump_every), static_cast<long long>(run.steps));
  out << b;
  out << "\nwrite_data      final.data\nwrite_restart   final.restart\n";
}

}  // namespace caps
