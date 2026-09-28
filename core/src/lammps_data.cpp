// LAMMPS data files carrying a CAPS force field: every term the evaluator knows, written with the LAMMPS style that
// computes the same energy. Kinds that mix forms (class II bonds with class I torsions, as in DL_FIELD's PCFF; CVFF
// with harmonic impropers; Morse and harmonic bonds) become hybrid styles, and the class II cross-term sections get
// "skip" lines for the types of the other sub-styles, as LAMMPS reads them.
//
// Checked term by term, energies and forces, against LAMMPS for each force-field family (bench/ff/check_data_lammps.py).
#include <filesystem>
#include <cstdio>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "caps/manybody.hpp"
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
  bool force_hybrid = false;               // the user asked for hybrid styles
  bool hybrid() const { return force_hybrid || styles.size() > 1; }
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
  bool charmm = false;                     // CHARMM (native): lj/charmmfsw, 1-4 pairs through dihedral charmmfsw weights
  bool hbond = false;                      // DREIDING hydrogen bonds: hbond/dreiding/lj overlaid on the pair style
  bool coreshell = false;                  // core-shell pairs (a shell on its core): the long-range Coulomb as .../cs
  bool cs_buck = false;                    // …with only Buckingham (and empty) pairs: one CORESHELL style for all of them
  std::vector<std::string> sw_types;       // per atom type: its Stillinger–Weber element name, or NULL (pair_style sw)
  std::vector<std::string> mb_types;       // per atom type: its element in the many-body potential file, or NULL
  // how the files are written (exact CAPS styles, or the force field's native ones)
  bool native = false, force_hybrid = false;
  std::string coul = "none";               // resolved Coulomb: none, dsf, long (with kspace), cut
  std::string kspace;                      // pppm or ewald when coul == long
  double kspace_accuracy = 1e-4;
  std::string pair_combined;               // native, one pair style for every pair: lj/cut/coul/long …
  std::vector<std::string> notes;
};

Layout build(const System& s, const ForceField& ff, const LammpsStyle& st = {}) {
  if (!ff.vsites.empty())
    throw FieldError(ff.name + ": virtual sites (Martini 3's tryptophan, ...) have no LAMMPS form; export to GROMACS instead");
  if (!ff.cbt.empty())
    throw FieldError(ff.name + ": combined bending–torsion dihedrals (Martini 3 polymers) have no LAMMPS form; export to GROMACS instead");
  if (!ff.lj_pairs.empty())
    throw FieldError(ff.name + ": explicit Lennard-Jones pairs (Martini 3 polymers) have no LAMMPS form; export to GROMACS instead");
  {   // exclusions LAMMPS can make: bonded 1-2 pairs, and 1-3 / 1-4 pairs when their scaling is 0
    const auto nb = s.neighbours();
    std::set<std::pair<uint32_t, uint32_t>> can;
    for (const auto& b : s.bonds) can.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
    if (!ff.keep13)
      for (uint32_t j = 0; j < nb.size(); ++j)
        for (uint32_t i : nb[j])
          for (uint32_t k : nb[j])
            if (i < k) can.insert({i, k});
    std::set<std::pair<uint32_t, uint32_t>> p14;   // looked up, not scanned: a 180 000-atom melt has half a million 1-4 pairs
    for (const auto& p : ff.pairs14) p14.insert({std::min(p[0], p[1]), std::max(p[0], p[1])});
    for (uint32_t i = 0; i < ff.excluded.size(); ++i)
      for (uint32_t j : ff.excluded[i])
        if (i < j && !can.count({i, j})) {
          const bool is14 = p14.count({i, j}) > 0;
          if (!is14)
            throw FieldError(ff.name + ": explicit exclusions between atoms that are not bonded (Martini 3's aromatic side chains, ...) have no "
                             "LAMMPS form; export to GROMACS instead");
        }
  }
  // CHARMM in its own styles: lj/charmmfsw with the 1-4 pairs (their own ε14, σ14) computed by dihedral charmmfsw
  const bool charmm = ff.lj_fsw && st.native;
  if (ff.lj_fsw && !st.native)
    throw FieldError(ff.name + ": CHARMM's force-switched Lennard-Jones and 1-4 terms are written in CHARMM's own LAMMPS styles only "
                     "(lj/charmmfsw, dihedral charmmfsw): choose the force field's own styles");
  if (charmm && (ff.lj14 != 1 || ff.coul14 != 1 || ff.keep13))
    throw FieldError(ff.name + ": dihedral charmmfsw computes 1-4 pairs in full; other 1-4 scaling has no LAMMPS form with lj/charmmfsw");
  if (!ff.lj14_types.empty() && !charmm)
    throw FieldError(ff.name + ": separate 1-4 Lennard-Jones parameters (GROMOS) have no exact LAMMPS form; export to GROMACS instead");
  Layout L;
  L.native = st.native;
  L.charmm = charmm;
  L.force_hybrid = st.hybrid;
  for (Kind* k : {&L.bonds, &L.angles, &L.dihedrals, &L.impropers}) k->force_hybrid = st.hybrid;
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
  // (not between atoms of a many-body potential: special_bonds would take the pair out of the neighbour list it reads)
  auto many_body = [&](uint32_t i) {
    if (!ff.manybody.on() || i >= ff.type_index.size()) return false;
    const size_t t = size_t(ff.type_index[i]);
    return t < ff.manybody.element.size() && !ff.manybody.element[t].empty();
  };
  for (const auto& b : s.bonds)
    if (b.i != b.j && !have.count({std::min(b.i, b.j), std::max(b.i, b.j)}) && !(many_body(b.i) && many_body(b.j))) {
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
  // native styles: CVFF's harmonic (one term K [1 + d cos(nφ)]), OPLS (½K1(1+cos φ) + ½K2(1−cos 2φ) + ½K3(1+cos 3φ) + ½K4(1−cos 4φ): a term v(1+cos(nφ−δ)) is
  // K_n = 2v with δ 0 for odd n, 180° for even n) and CHARMM (one dihedral line per term, K n d with d a whole degree,
  // weight 0: the 1-4 pairs come from special_bonds); a quadruple that does not fit stays a Fourier sum
  const std::string nd = st.native ? ff.native_dihedral : "";
  // a constant term (n = 0: TraPPE's c0, v(1 + cos δ)) has no opls coefficient; it moves no atom, so it is left out and
  // its total over the structure's torsions reported (LAMMPS's dihedral energy is lower by exactly that)
  double opls_constant = 0;
  auto opls_of = [&](const std::vector<const TorsionTerm*>& v, std::string& coef) {
    double K[4] = {0, 0, 0, 0}, c0 = 0;
    for (const auto* t : v) {
      if (t->n == 0) { c0 += t->v * (1 + std::cos(t->delta)); continue; }
      if (t->n < 1 || t->n > 4) return false;
      const double want = (t->n % 2) ? 0.0 : 180.0, d = std::fmod(std::fabs(t->delta * R2D), 360.0);
      if (std::fabs(d - want) > 1e-6 && std::fabs(d - want - 360) > 1e-6) return false;
      K[t->n - 1] += 2 * t->v;
    }
    coef = num({K[0], K[1], K[2], K[3]});
    opls_constant += c0;
    return true;
  };
  size_t fourier_left = 0;
  std::vector<std::pair<std::array<uint32_t, 4>, std::string>> native_lines;   // (quadruple, style|coef)
  if (L.charmm) {
    // CHARMM: every 1-4 pair through its torsions' weights (special_bonds charmm leaves them out of the pair list): the
    // first line of each quadruple carries 1 / (quadruples sharing its end atoms), so a pair two torsions reach (a
    // six-membered ring's) counts once; a 1-4 pair with no torsion term gets a zero one to carry it
    std::set<std::pair<uint32_t, uint32_t>> p14;
    for (const auto& p : ff.pairs14) p14.insert({std::min(p[0], p[1]), std::max(p[0], p[1])});
    std::map<std::pair<uint32_t, uint32_t>, int> reach;
    for (const auto& k : order) ++reach[{std::min(k[0], k[3]), std::max(k[0], k[3])}];
    std::vector<std::vector<uint32_t>> nbl(s.atoms.size());
    for (const auto& bd : s.bonds) { nbl[bd.i].push_back(bd.j); nbl[bd.j].push_back(bd.i); }
    static const TorsionTerm kZero{0, 0, 0, 0, 0.0, 1, 0.0};
    for (const auto& pr : p14) {
      if (reach.count(pr)) continue;
      bool placed = false;   // a path i-j-k-l for the pair
      for (uint32_t j : nbl[pr.first]) {
        for (uint32_t k : nbl[j])
          if (k != pr.first && std::find(nbl[k].begin(), nbl[k].end(), pr.second) != nbl[k].end() && k != pr.second && j != pr.second) {
            const std::array<uint32_t, 4> q{pr.first, j, k, pr.second};
            order.push_back(q);
            quad[q].push_back(&kZero);
            ++reach[pr];
            placed = true;
            break;
          }
        if (placed) break;
      }
      if (!placed) throw FieldError(ff.name + ": a 1-4 pair with no bonded path for its torsion (charmmfsw needs one)");
    }
    for (const auto& k : order) {
      const auto pr = std::make_pair(std::min(k[0], k[3]), std::max(k[0], k[3]));
      const double w = p14.count(pr) ? 1.0 / reach[pr] : 0.0;
      bool first = true;
      for (const auto* t : quad[k]) {
        const double dd = t->delta * R2D;
        if (t->n < 0 || std::fabs(dd - std::round(dd)) > 1e-9)
          throw FieldError(ff.name + ": a torsion with a non-integer phase or negative multiplicity has no dihedral charmmfsw form");
        L.dihedrals.add("charmmfsw", num({t->v}) + " " + std::to_string(t->n) + " " + std::to_string(((std::lround(dd) % 360) + 360) % 360) + num({first ? w : 0.0}),
                        {}, {k[0], k[1], k[2], k[3]}, lab({k[0], k[1], k[2], k[3]}));
        first = false;
      }
    }
    order.clear();
  }
  for (const auto& k : order) {
    const auto& v = quad[k];
    std::string c;
    if (nd == "opls" && opls_of(v, c)) {
      L.dihedrals.add("opls", c, {}, {k[0], k[1], k[2], k[3]}, lab({k[0], k[1], k[2], k[3]}));
      continue;
    }
    if (nd == "class2") {
      // a class II force field's plain torsions (the diagonal PCFF / COMPASS files) as its own dihedral class2 with no
      // cross terms: v [1 + cos(nφ − δ)] = K_n [1 − cos(nφ − φ_n)], K_n = v, φ_n = δ + 180°, n 1 to 3
      double K[3] = {0, 0, 0}, P[3] = {0, 0, 0};
      bool fits = true;
      for (const auto* t : v) {
        if (t->n < 1 || t->n > 3 || K[t->n - 1] != 0) { fits = false; break; }
        K[t->n - 1] = t->v;
        P[t->n - 1] = std::fmod(t->delta * R2D + 180.0 + 720.0, 360.0);
      }
      if (fits) {
        L.dihedrals.add("class2", num({K[0], P[0], K[1], P[1], K[2], P[2]}),
                        {num({0, 0, 0, 0}), num({0, 0, 0, 0, 0, 0, 0, 0}), num({0, 0, 0, 0, 0, 0, 0, 0}), num({0, 0, 0}), num({0, 0, 0})},
                        {k[0], k[1], k[2], k[3]}, lab({k[0], k[1], k[2], k[3]}));
        continue;
      }
    }
    if (nd == "harmonic" && v.size() == 1 && v[0]->n >= 0 && std::fabs(std::fabs(std::cos(v[0]->delta)) - 1) < 1e-9) {
      // CVFF as msi2lmp writes it: K [1 + d cos(nφ)], d = +1 (phase 0) or −1 (phase 180°)
      L.dihedrals.add("harmonic", num({v[0]->v}) + (std::cos(v[0]->delta) > 0 ? " 1 " : " -1 ") + std::to_string(v[0]->n), {},
                      {k[0], k[1], k[2], k[3]}, lab({k[0], k[1], k[2], k[3]}));
      continue;
    }
    if (nd == "charmm") {
      bool ok = true;
      for (const auto* t : v) ok = ok && t->n >= 0 && std::fabs(t->delta * R2D - std::round(t->delta * R2D)) < 1e-9;
      if (ok) {
        for (const auto* t : v)
          L.dihedrals.add("charmm", num({t->v}) + " " + std::to_string(t->n) + " " + std::to_string(((std::lround(t->delta * R2D) % 360) + 360) % 360) + " 0.0", {},
                          {k[0], k[1], k[2], k[3]}, lab({k[0], k[1], k[2], k[3]}));
        continue;
      }
    }
    std::string f = " " + std::to_string(v.size());
    for (const auto* t : v) {
      if (t->n < 0) throw FieldError("a torsion with multiplicity < 0 has no LAMMPS fourier form");
      f += num({t->v}) + " " + std::to_string(t->n) + num({t->delta * R2D});
    }
    L.dihedrals.add("fourier", f, {}, {k[0], k[1], k[2], k[3]}, lab({k[0], k[1], k[2], k[3]}));
    ++fourier_left;
  }
  if (opls_constant != 0) {
    char nb[200];
    std::snprintf(nb, sizeof nb, "torsion constants (n = 0 terms, e.g. TraPPE's c0) have no opls form and are left out: forces unchanged, "
                  "the dihedral energy lower by %.10g kcal/mol", opls_constant);
    L.notes.push_back(nb);
  }
  if (!nd.empty() && nd != "fourier" && fourier_left)
    L.notes.push_back(std::to_string(fourier_left) + " dihedral types have terms the " + nd + " style cannot hold: written as fourier (a hybrid dihedral style)");
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
      if (st.native && ff.native_improper == "fourier")   // UFF's own form: K [1 − cos ω] as improper fourier over all three
        L.impropers.add("fourier", num({v.kw / 3, 1.0, -1.0, 0.0}) + " 1", {}, {v.c, v.a, v.b, v.d}, lab({v.c, v.a, v.b, v.d}));
      else
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
  if (ff.manybody.on()) {   // a literature many-body potential overlays the pair terms (its types' Lennard-Jones among themselves is zero)
    if (L.gromacs) throw FieldError(ff.name + ": lj/gromacs beside a many-body potential has no LAMMPS style");
    L.mb_types.assign(size_t(nt), "NULL");
    for (size_t t = 0; t < size_t(nt) && t < ff.manybody.element.size(); ++t)
      if (!ff.manybody.element[t].empty()) L.mb_types[t] = ff.manybody.element[t];
    L.pair_hybrid = true;
  }
  L.periodic = s.cell.valid();
  if (st.hybrid && !L.gromacs) L.pair_hybrid = true;
  // shell models: a core and its shell bonded at (nearly) the same place; coul/long would evaluate their excluded pair's
  // Ewald correction at r = 0 (nan), CORESHELL's coul/long/cs does it stably
  for (const auto& b : s.bonds) {
    const Vec3 d = s.cell.valid() ? s.cell.minimum_image(s.atoms[b.j].pos - s.atoms[b.i].pos) : s.atoms[b.j].pos - s.atoms[b.i].pos;
    if (dot(d, d) < 0.25 * 0.25) { L.coreshell = true; break; }
  }
  // a shell model's pairs all Buckingham (the others empty): CORESHELL's buck/coul/long/cs, or born/coul/dsf/cs with DSF
  // (Buckingham is Born-Mayer-Huggins with σ = 0, D = 0), so no sub-style meets a core and its shell at r = 0
  if (L.coreshell && !L.gromacs && L.sw_types.empty()) {
    bool all = true, any_buck = false;
    for (int a2 = 0; a2 < nt && all; ++a2)
      for (int b2 = a2; b2 < nt && all; ++b2) {
        auto it = ff.pair_func.find({a2, b2});
        if (it != ff.pair_func.end()) { all = it->second.form == 1; any_buck = true; }
        else all = mixed_pair(ff, a2, b2).eps == 0;
      }
    if (all && any_buck) { L.cs_buck = true; L.pair_hybrid = false; }
  }
  L.hbond = ff.hbond.on();
  if (L.hbond) {
    if (L.gromacs || L.charmm) throw FieldError(ff.name + ": DREIDING hydrogen bonds with this pair style have no LAMMPS form");
    L.pair_hybrid = true;   // hbond/dreiding/lj overlays the Lennard-Jones and Coulomb pairs
  }
  return L;
}

// The native layout's Coulomb and pair style: the force field's long-range sum (PPPM unless asked otherwise) when it is
// long-range, charged and periodic; otherwise what the user chose. One pair style for every pair when it can be.
void resolve_native(Layout& L, const ForceField& ff, const LammpsStyle& st, bool charged) {
  const bool ff_long = ff.native_pair.find("long") != std::string::npos || ff.native_pair.find("coul") == std::string::npos;
  std::string c = st.coulomb;
  if (!charged) c = "none";
  else if (c == "auto") c = ff_long && L.periodic ? "pppm" : L.periodic ? "dsf" : "cut";
  if ((c == "pppm" || c == "ewald") && !L.periodic) {
    L.notes.push_back("no periodic cell: Coulomb cut off, not " + c);
    c = "cut";
  }
  L.coul = c == "pppm" || c == "ewald" ? "long" : c;
  if (L.coul == "long") L.kspace = c;
  L.kspace_accuracy = st.kspace_accuracy;
  if (L.cs_buck) {   // a shell model: one CORESHELL style for every pair
    L.pair_combined = L.coul == "long" ? "buck/coul/long/cs" : "born/coul/dsf/cs";
    if (L.coul == "cut" || L.coul == "none") L.notes.push_back("shell model without long-range Coulomb: born/coul/dsf/cs (DSF)");
    if (L.coul != "long") L.coul = "dsf";
    L.pair_styles = {L.pair_combined};
    L.pair_hybrid = st.hybrid || !L.mb_types.empty();
    return;
  }
  if (L.charmm) {   // CHARMM: the force-switched LJ with its long-range sum, else CHARMM's force-shifted Coulomb
    L.pair_combined = std::string("lj/charmmfsw/coul/") + (L.coul == "long" ? "long" : "charmmfsh");
    if (L.coul == "dsf" || L.coul == "cut") L.notes.push_back("Coulomb: CHARMM's force shift (coul/charmmfsh), not " + c);
    L.pair_styles = {L.pair_combined};
    L.pair_hybrid = st.hybrid || !L.mb_types.empty();
    return;
  }
  // SDK / SPICA: every pair lj/sdk, with its long-range Coulomb as the one style lj/sdk/coul/long (as its own inputs write
  // it); without charges lj/sdk alone. Other Coulomb forms stay an overlay.
  if (L.sdk && L.pair_styles.size() == 1 && L.sw_types.empty() && (L.coul == "long" || L.coul == "none")) {
    L.pair_combined = L.coul == "long" ? "lj/sdk/coul/long" : "lj/sdk";
    L.pair_styles = {L.pair_combined};
    L.pair_hybrid = st.hybrid || L.hbond || !L.mb_types.empty();
    return;
  }
  // one style for all pairs: lj/cut or lj/class2 with its Coulomb (lj/class2 has no DSF form)
  if (ff.pair_func.empty() && L.sw_types.empty() && !L.gromacs && !(L.pair_base == "lj/class2" && L.coul == "dsf")) {
    L.pair_combined = L.pair_base + (L.coul == "none" ? "" : "/coul/" + L.coul);
    L.pair_styles = {L.pair_combined};
    L.pair_hybrid = st.hybrid || L.hbond || !L.mb_types.empty();   // one style (lj/class2/coul/long too): hybrid only when asked or overlaid
  }
}

// Pair-style arguments as force-field input files write them: "12.000000" (fixed, six decimals, as the dsf α too)
std::string fmt_args(std::initializer_list<double> v) {
  std::string r;
  char b[40];
  for (double x : v) {
    std::snprintf(b, sizeof b, " %.6f", x);
    r += b;
  }
  return r;
}

// A k-space accuracy as "1.0e-4": one decimal, the exponent without its leading zero
std::string fmt_accuracy(double x) {
  char b[40];
  std::snprintf(b, sizeof b, "%.1e", x);
  std::string s = b;
  const auto e = s.find('e');
  if (e != std::string::npos) {
    std::size_t k = e + 1;
    if (k < s.size() && (s[k] == '-' || s[k] == '+')) ++k;
    while (k + 1 < s.size() && s[k] == '0') s.erase(k, 1);
    if (s[e + 1] == '+') s.erase(e + 1, 1);
  }
  return s;
}

// DREIDING hydrogen bonds: "hbond/dreiding/lj power r_in r_out angle", and one pair_coeff per (donor type, acceptor type):
// the lower type index first, the flag i or j naming the donor, the hydrogen type, ε σ n.
std::string hbond_style(const ForceField& ff) {
  char b[160];
  std::snprintf(b, sizeof b, "hbond/dreiding/lj %d %.6g %.6g %.6g", ff.hbond.power, ff.hbond.inner, ff.hbond.outer, ff.hbond.angle_deg);
  return b;
}
std::vector<std::string> hbond_lines(const ForceField& ff) {
  std::vector<std::string> r;
  char b[240];
  for (const auto& [key, p] : ff.hbond.param) {
    const int d = key.first, a = key.second, h = ff.hbond.htype.at(key);
    std::snprintf(b, sizeof b, "%d %d hbond/dreiding/lj %d %s %.10g %.10g %d  # donor %s, acceptor %s, hydrogen %s", std::min(d, a) + 1,
                  std::max(d, a) + 1, h + 1, d <= a ? "i" : "j", p[0], p[1], int(p[2]), ff.type_names[size_t(d)].c_str(),
                  ff.type_names[size_t(a)].c_str(), ff.type_names[size_t(h)].c_str());
    r.push_back(b);
  }
  return r;
}

// PME is used (and written) only for periodic cells, as the evaluator does.
bool pme(const EnergyOptions& e, const Layout& L) { return e.electrostatics == EnergyOptions::Electrostatics::PME && L.periodic && !L.gromacs; }

// The LAMMPS commands (after units / atom_style) that reproduce CAPS's energy with the data file.
// The many-body style as pair_style names it, with its arguments (airebo 3.0 1 1).
std::string mb_style_word(const ForceField& ff) { return ff.manybody.style + (ff.manybody.args.empty() ? "" : " " + ff.manybody.args); }

std::vector<std::string> style_lines(const Layout& L, const ForceField& ff, const EnergyOptions& e) {
  std::vector<std::string> r;
  char b[400];
  if (L.native && !L.pair_combined.empty()) {
    // the force field's own form: one pair style (hybrid when asked), its long-range sum by PPPM or Ewald
    const std::string args = L.charmm ? fmt_args({ff.lj_inner, e.cutoff}) : L.coul == "dsf" ? fmt_args({e.dsf_alpha, e.cutoff}) : fmt_args({e.cutoff});
    const std::string mb = L.mb_types.empty() ? "" : " " + mb_style_word(ff);   // a many-body potential overlaid
    if (L.hbond || !mb.empty()) {   // DREIDING: the hydrogen bond overlaid on the Lennard-Jones and Coulomb pairs
      r.push_back("pair_style hybrid/overlay " + (L.hbond ? hbond_style(ff) + " " : std::string()) + L.pair_combined + args + mb);
      if (e.tail && L.periodic) r.push_back("pair_modify pair " + L.pair_combined + " tail yes");
      else if (!e.tail) r.push_back("pair_modify pair " + L.pair_combined + " shift yes");
    } else {
      r.push_back("pair_style " + std::string(L.pair_hybrid ? "hybrid " : "") + L.pair_combined + args);
      if (L.charmm) {}   // switched to zero at the cut-off: no tail, no shift
      else if (L.sdk) { if (!e.tail) r.push_back("pair_modify shift yes"); }   // lj/sdk has no tail correction (nor has CAPS for it)
      else if (e.tail && L.periodic) r.push_back("pair_modify tail yes");
      else if (!e.tail) r.push_back("pair_modify shift yes");
    }
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
    // LAMMPS's own keyword where the force field has one: amber (lj 0 0 0.5, coul 0 0 5/6 exactly), dreiding (0 0 1)
    if (L.charmm)   // 1-2, 1-3, 1-4 all out of the pair list: dihedral charmmfsw adds the 1-4 pairs
      r.push_back("special_bonds charmm");
    else if (ff.native_special == "amber" && !ff.keep13 && ff.lj14 == 0.5 && std::fabs(ff.coul14 - 5.0 / 6.0) < 1e-9)
      r.push_back("special_bonds amber");
    else if (ff.native_special == "dreiding" && !ff.keep13 && ff.lj14 == 1 && ff.coul14 == 1)
      r.push_back("special_bonds dreiding");
    else {
      std::snprintf(b, sizeof b, "special_bonds lj 0.0 %s %.6f coul 0.0 %s %.6f", ff.keep13 ? "1.0" : "0.0", ff.lj14, ff.keep13 ? "1.0" : "0.0", ff.coul14);
      r.push_back(b);
    }
    if (L.coul == "long") r.push_back("kspace_style " + L.kspace + " " + fmt_accuracy(L.kspace_accuracy));
    return r;
  }
  if (L.gromacs) {   // MARTINI: the GROMACS switch for LJ (and Coulomb), inner and outer radii
    if (e.coulomb && ff.coul_gromacs)
      std::snprintf(b, sizeof b, "pair_style lj/gromacs/coul/gromacs %.6g %.6g %.6g %.6g", ff.lj_inner, e.cutoff, std::max(ff.coul_inner, 1e-6), e.cutoff);
    else if (e.coulomb) throw FieldError("lj/gromacs needs the GROMACS Coulomb form");
    else std::snprintf(b, sizeof b, "pair_style lj/gromacs %.6g %.6g", ff.lj_inner, e.cutoff);
    r.push_back(b);
  } else if (L.cs_buck) {   // a shell model (CORESHELL): Buckingham with its Coulomb, stable for a core on its shell
    if (e.coulomb && pme(e, L)) std::snprintf(b, sizeof b, "pair_style buck/coul/long/cs %.6g", e.cutoff);
    else std::snprintf(b, sizeof b, "pair_style born/coul/dsf/cs %.6g %.6g", e.dsf_alpha, e.cutoff);
    r.push_back(b);
  } else if (!L.pair_hybrid) {
    if (e.coulomb && pme(e, L)) std::snprintf(b, sizeof b, "pair_style lj/cut/coul/long %.6g", e.cutoff);
    else if (e.coulomb) std::snprintf(b, sizeof b, "pair_style lj/cut/coul/dsf %.6g %.6g", e.dsf_alpha, e.cutoff);
    else std::snprintf(b, sizeof b, "pair_style lj/cut %.6g", e.cutoff);
    r.push_back(b);
  } else {
    std::string p = "pair_style hybrid/overlay";
    if (L.hbond) p += " " + hbond_style(ff);
    for (const auto& st : L.pair_styles) {
      std::snprintf(b, sizeof b, " %s %.6g", st.c_str(), e.cutoff);
      p += b;
    }
    if (e.coulomb && pme(e, L)) {
      std::snprintf(b, sizeof b, L.coreshell ? " coul/long/cs %.6g" : " coul/long %.6g", e.cutoff);
      p += b;
    } else if (e.coulomb) {
      std::snprintf(b, sizeof b, " coul/dsf %.6g %.6g", e.dsf_alpha, e.cutoff);
      p += b;
    }
    if (!L.sw_types.empty()) p += " sw";
    if (!L.mb_types.empty()) p += " " + mb_style_word(ff);
    r.push_back(p);
  }
  // CAPS: with tail corrections the potentials are truncated at the cut-off (plus the tail when there is a cell);
  // without, they are shifted to zero there
  // lj/gromacs is zero at the cut-off by itself; lj/sdk has no tail correction (CAPS adds none for it either)
  if (L.gromacs || (L.cos2 && L.pair_styles.size() == 1)) {
  } else if (L.cs_buck && L.pair_combined.find("dsf") != std::string::npos) {   // born/coul/dsf/cs: no tail correction in LAMMPS
    if (!e.tail) r.push_back("pair_modify shift yes");
  } else if (L.hbond) {   // the tail or shift for the Lennard-Jones sub-style only
    if (!e.tail) r.push_back("pair_modify pair " + L.pair_base + " shift yes");
    else if (L.periodic) r.push_back("pair_modify pair " + L.pair_base + " tail yes");
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
  if (L.native)   // the force field's own form of the line (MARTINI: lj 0.0 1.0 1.0)
    std::snprintf(b, sizeof b, "special_bonds lj 0.0 %s %.6f coul 0.0 %s %.6f", ff.keep13 ? "1.0" : "0.0", ff.lj14, ff.keep13 ? "1.0" : "0.0", ff.coul14);
  else
    std::snprintf(b, sizeof b, "special_bonds lj 0 %d %.10g coul 0 %d %.10g", ff.keep13 ? 1 : 0, ff.lj14, ff.keep13 ? 1 : 0, ff.coul14);
  r.push_back(b);
  if (e.coulomb && pme(e, L)) {
    // CAPS's PME with its own β; LAMMPS's Ewald sum to the same accuracy reaches the same total electrostatics
    if (L.native) std::snprintf(b, sizeof b, "kspace_style %s %s", L.kspace.c_str(), fmt_accuracy(L.kspace_accuracy).c_str());
    else std::snprintf(b, sizeof b, "kspace_style ewald %.3g", std::max(1e-12, e.ewald_rtol * 0.01));
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
      std::string coef, style = L.pair_combined.empty() ? L.pair_base : L.pair_combined;
      const int f = it != ff.pair_func.end() ? it->second.form : 0;
      if (f >= kPairSdk96 && f <= kPairSdk125) {
        static const char* nm[] = {"lj9_6", "lj12_4", "lj12_6", "lj12_5"};
        style = L.pair_combined.empty() ? "lj/sdk" : L.pair_combined;   // lj/sdk/coul/long when that is the one style
        coef = std::string(" ") + nm[f - kPairSdk96] + num({it->second.a, it->second.b});
      } else if (f == kPairGromacs) {
        coef = num({it->second.a, it->second.b});
      } else if (f == kPairCos2 || f == kPairCos2Wca) {
        style = "cosine/squared";
        coef = num({it->second.a, it->second.b, it->second.c}) + (f == kPairCos2Wca ? " wca" : "");
      } else if (L.cs_buck) {   // Buckingham (A ρ C), or as Born–Mayer–Huggins (A ρ σ=0 C D=0); empty pairs zero
        const bool born = L.pair_combined.find("born") != std::string::npos;
        const double A = it != ff.pair_func.end() ? it->second.a : 0.0, rho = it != ff.pair_func.end() ? it->second.b : 1.0,
                     C = it != ff.pair_func.end() ? it->second.c : 0.0;
        style = L.pair_combined;
        coef = born ? num({A, rho, 0.0, C, 0.0}) : num({A, rho, C});
      } else if (it != ff.pair_func.end()) {
        style = it->second.form == 1 ? "buck" : "morse";
        coef = num({it->second.a, it->second.b, it->second.c});
      } else if (L.charmm) {   // ε σ ε14 σ14, both pairs mixed by the force field's rule
        const PairType pt = mixed_pair(ff, int(a2), int(b2));
        PairType p14 = pt;
        if (!ff.lj14_types.empty()) {
          ForceField f14;
          f14.mixing = ff.mixing;
          f14.lj = ff.lj14_types;
          p14 = mixed_pair(f14, int(a2), int(b2));
        }
        coef = num({pt.eps, pt.sigma, p14.eps, p14.sigma});
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
// The many-body potential file next to a data file.
std::string mb_path(const ForceField& ff, const std::string& data_path) {
  const size_t slash = data_path.find_last_of('/');
  return (slash == std::string::npos ? std::string() : data_path.substr(0, slash + 1)) + manybody_file_name(ff.manybody);
}

std::vector<std::string> after_read(const Layout& L, const EnergyOptions& e, const std::string& data_path,
                                    const std::set<std::pair<int, int>>& ff_excl = {}, const ForceField* ff = nullptr) {
  std::vector<std::string> r;
  if (L.pair_hybrid && e.coulomb && L.pair_combined.empty())
    r.push_back(pme(e, L) ? (L.coreshell ? "pair_coeff * * coul/long/cs" : "pair_coeff * * coul/long") : "pair_coeff * * coul/dsf");
  for (const auto& [a, b] : ff_excl) r.push_back("neigh_modify exclude type " + std::to_string(a + 1) + " " + std::to_string(b + 1));
  if (!L.sw_types.empty()) {
    std::string l = "pair_coeff * * sw " + sw_path(data_path);
    for (const auto& t : L.sw_types) l += " " + t;
    r.push_back(l);
  }
  if (!L.mb_types.empty() && ff) {   // every type mapped: its element in the file, or NULL
    std::string l = "pair_coeff * * " + ff->manybody.style + " " + mb_path(*ff, data_path);
    for (const auto& t : L.mb_types) l += " " + t;
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

// A title read back from a CAPS file carries the old header line: keep the description only.
}  // namespace

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

namespace {

// The energy options and layout a pair of LAMMPS files is written with (the same for the data file and the input).
Layout prepare(const System& s, const ForceField& ff, EnergyOptions& e, const LammpsStyle& st) {
  if (ff.cutoff > 0) e.cutoff = ff.cutoff;   // the model's own cut-off (MARTINI)
  const bool charged = !std::all_of(ff.charge.begin(), ff.charge.end(), [](double q) { return q == 0; });
  // no charges, no Coulomb term (LAMMPS refuses an Ewald sum on an uncharged system; the energy is the same)
  if (!charged) e.coulomb = false;
  if (ff.lj_shift || ff.lj_fsw) e.tail = false;   // Martini 3: shifted at the cut-off; CHARMM: switched
  if (ff.coul_rf && e.coulomb)
    throw FieldError(ff.name + ": reaction-field Coulomb (GROMOS, Martini 3) has no LAMMPS pair style; export to GROMACS instead");
  Layout L = build(s, ff, st);
  if (st.native) {
    // a coarse-grained model's own cut-off belongs to the model (MARTINI 12 Å, SDK 15 Å); otherwise the user's, else the
    // force field file's
    if (ff.cutoff > 0) {
      if (st.cutoff > 0 && std::fabs(st.cutoff - ff.cutoff) > 1e-9)
        L.notes.push_back(ff.name + "'s cut-off is part of the model: " + std::to_string(ff.cutoff) + " Å kept");
    } else if (st.cutoff > 0) e.cutoff = st.cutoff;
    else if (ff.native_cutoff > 0) e.cutoff = ff.native_cutoff;
    if (st.tail >= 0) e.tail = st.tail == 1;
    resolve_native(L, ff, st, charged && e.coulomb);
    e.coulomb = L.coul != "none";
    if (L.coul == "long") e.electrostatics = EnergyOptions::Electrostatics::PME;
    else if (L.coul == "dsf") e.electrostatics = EnergyOptions::Electrostatics::DSF;
  }
  // a shell model in CAPS's own styles: the CORESHELL style that matches the Coulomb in use
  if (L.cs_buck && !L.native) L.pair_combined = e.coulomb && pme(e, L) ? "buck/coul/long/cs" : "born/coul/dsf/cs";
  if (L.cs_buck && L.pair_combined == "born/coul/dsf/cs" && e.tail && L.periodic)
    L.notes.push_back("born/coul/dsf/cs has no tail correction: LAMMPS's van der Waals leaves out CAPS's tail term (use --no-tail, or Ewald / PPPM)");
  return L;
}

}  // namespace

bool lammps_metal_units(const ForceField& ff, const LammpsStyle& st) {
  if (st.units != "auto" && st.units != "real" && st.units != "metal") throw FieldError("LAMMPS units " + st.units + ": real, metal or auto");
  const bool need = ff.manybody.on() && ff.manybody.metal_only;
  if (need && st.units == "real")
    throw FieldError(ff.manybody.style + " is read by LAMMPS in metal units only: write this system in metal units (eV, ps, bar)");
  return st.units == "metal" || (st.units == "auto" && need);
}

ForceField forcefield_in_metal_units(const ForceField& ff0) {
  constexpr double f = 1.0 / 23.060549;   // kcal/mol → eV, as LAMMPS converts potential files
  ForceField ff = ff0;
  for (auto& p : ff.lj) p.eps *= f;
  for (auto& p : ff.lj14_types) p.eps *= f;
  for (auto& [k, p] : ff.pair_override) p.eps *= f;
  for (auto& [k, p] : ff.pair_func) {
    if (p.form == 1) p.a *= f, p.c *= f;   // Buckingham A, C
    else p.a *= f;                         // Morse D0; SDK, lj/gromacs and cosine/squared ε
  }
  for (auto& t : ff.lj_pairs) t.eps *= f;
  for (auto& t : ff.bonds) t.k *= f;
  for (auto& t : ff.angles) t.kt *= f;
  for (auto& t : ff.dihedrals) t.v *= f;
  for (auto& t : ff.impropers) t.v *= f;
  for (auto& t : ff.impropers_dlpoly) t.v *= f;
  for (auto& t : ff.impropers_harmonic) t.k2 *= f;
  for (auto& t : ff.inversions) t.kw *= f;
  for (auto& t : ff.bonds_x) {
    t.a *= f;                        // Morse D, GROMOS K, FENE K
    if (t.form == 3) t.c *= f;       // FENE's WCA ε
  }
  for (auto& t : ff.angles_x) t.a *= f;
  for (auto& t : ff.urey_bradley) t.kub *= f;
  for (auto& t : ff.cbt) for (double& a : t.a) a *= f;
  for (auto& t : ff.bonds2) t.k2 *= f, t.k3 *= f, t.k4 *= f;
  for (auto& t : ff.angles2) t.k2 *= f, t.k3 *= f, t.k4 *= f, t.bb_m *= f, t.ba_n1 *= f, t.ba_n2 *= f;
  for (auto& t : ff.dihedrals2) {
    t.k1 *= f, t.k2 *= f, t.k3 *= f, t.aat_m *= f, t.bb13_n *= f;
    for (int n = 0; n < 3; ++n) t.mbt[n] *= f, t.ebt_b[n] *= f, t.ebt_c[n] *= f, t.at_d[n] *= f, t.at_e[n] *= f;
  }
  for (auto& t : ff.impropers2) t.kchi *= f, t.m1 *= f, t.m2 *= f, t.m3 *= f;
  for (auto& [k, p] : ff.hbond.param) p[0] *= f;   // ε
  ff.sw.eps *= f;
  return ff;
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

void write_lammps_data_ff(const System& s, const ForceField& ff0, const EnergyOptions& e0, const std::string& path, bool pair_coeffs,
                          const LammpsStyle& st) {
  const bool metal = lammps_metal_units(ff0, st);
  const ForceField mff = metal ? forcefield_in_metal_units(ff0) : ForceField{};
  const ForceField& ff = metal ? mff : ff0;
  EnergyOptions e = e0;
  const Layout L = prepare(s, ff, e, st);
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  if (!L.sw_types.empty()) write_sw_file(L, ff, sw_path(path));
  if (!L.mb_types.empty()) write_manybody_file(ff.manybody, std::filesystem::path(path).parent_path().string());
  char buf[512];
  const std::vector<const Kind*> kinds = {&L.bonds, &L.angles, &L.dihedrals, &L.impropers};

  // the LAMMPS data format: one title line, the counts, the box, then the sections (atom_style full, units real)
  out << "CAPS · " << clean_title(s.title, ff.name) << " · " << ff.name << (metal ? " · units metal (eV)" : "") << "\n\n";
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
  if (!c.valid()) {   // no cell: a box around the atoms, 25 Å of vacuum on every side (more than any cut-off), at least 100 Å
    Vec3 mn{1e30, 1e30, 1e30}, mx{-1e30, -1e30, -1e30};
    for (const auto& at : s.atoms)
      for (int k = 0; k < 3; ++k) mn[k] = std::min(mn[k], at.pos[k]), mx[k] = std::max(mx[k], at.pos[k]);
    if (s.atoms.empty()) mn = {0, 0, 0}, mx = {0, 0, 0};
    Vec3 w;
    for (int k = 0; k < 3; ++k) {
      w[k] = std::max(100.0, mx[k] - mn[k] + 50.0);
      lo[k] = 0.5 * (mn[k] + mx[k]) - 0.5 * w[k];
    }
    a = {w[0], 0, 0}, b = {0, w[1], 0}, cc = {0, 0, w[2]};
  }
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
    out << "\nPairIJ Coeffs  # " << (!L.pair_combined.empty() ? (L.hbond ? "hybrid/overlay" : L.pair_hybrid ? "hybrid" : L.pair_combined) : L.pair_hybrid ? std::string("hybrid/overlay") : L.gromacs ? std::string(e.coulomb ? "lj/gromacs/coul/gromacs" : "lj/gromacs")
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
  // molecule ids: the file's own, else the bonded fragments (a PDB's chain ids are not molecules)
  std::vector<int> frag;
  if (!s.has_mol) frag = s.molecules();
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const auto& at = s.atoms[i];
    const Vec3 fr = c.valid() ? c.to_fractional(at.pos) : Vec3{0, 0, 0};
    int im[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) im[k] = c.valid() && c.periodic[k] ? int(std::floor(fr[k])) : 0;
    const Vec3 w = at.pos - (c.a * im[0] + c.b * im[1] + c.c * im[2]);
    std::snprintf(buf, sizeof buf, "%zu %lld %d %.8f %.10f %.10f %.10f %d %d %d\n", i + 1, static_cast<long long>(s.has_mol ? at.mol : int64_t(frag[i]) + 1), ff.type_index[i] + 1,
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

// fix shake on the bonds CAPS constrains: to hydrogen by mass (and water's H–O–H angle by its type), or every bond type
std::string shake_fix(const System& s, const ForceField& ff, const Layout& L, ConstraintMode mode, const std::string& group) {
  if (mode == ConstraintMode::None || L.bonds.types.empty()) return {};
  char b[160];
  std::string what;
  if (mode == ConstraintMode::AllBonds) {
    what = "b";
    for (size_t t = 1; t <= L.bonds.types.size(); ++t) what += " " + std::to_string(t);
  } else {
    double mh = 0;
    for (size_t i = 0; i < s.atoms.size(); ++i)
      if (s.atoms[i].element == 1) mh = std::max(mh, ff.mass[i]);
    if (mh <= 0) return {};
    std::snprintf(b, sizeof b, "m %.4g", mh + 0.05);
    what = b;
    const auto nb = s.neighbours();
    std::set<int> water;
    for (size_t k = 0; k < L.angles.term_atoms.size(); ++k) {
      const auto& t = L.angles.term_atoms[k];
      if (t.size() != 3) continue;
      const uint32_t o = t[1];
      if (s.atoms[o].element == 8 && nb[o].size() == 2 && s.atoms[t[0]].element == 1 && s.atoms[t[2]].element == 1) water.insert(L.angles.term_type[k]);
    }
    if (!water.empty()) {
      what += " a";
      for (int t : water) what += " " + std::to_string(t);
    }
  }
  std::snprintf(b, sizeof b, "fix             hold_bonds %s shake 1.0e-6 100 0 ", group.c_str());
  return b + what + "\n";
}

std::string lammps_shake_fix(const System& s, const ForceField& ff0, const EnergyOptions& e0, ConstraintMode mode, const std::string& group,
                             const LammpsStyle& st) {
  const bool metal = lammps_metal_units(ff0, st);
  const ForceField mff = metal ? forcefield_in_metal_units(ff0) : ForceField{};
  const ForceField& ff = metal ? mff : ff0;
  EnergyOptions e = e0;
  const Layout L = prepare(s, ff, e, st);
  return shake_fix(s, ff, L, mode, group);
}

void write_lammps_input(const System& s, const ForceField& ff0, const EnergyOptions& e0, const std::string& data_path, const std::string& path,
                        int64_t held_mol, bool pair_coeffs, const LammpsRun& run, const LammpsStyle& st, std::vector<std::string>* notes) {
  const bool metal = lammps_metal_units(ff0, st);
  const ForceField mff = metal ? forcefield_in_metal_units(ff0) : ForceField{};
  const ForceField& ff = metal ? mff : ff0;
  const double tu = metal ? 1e-3 : 1.0;          // fs → ps
  const double pu = metal ? 1.01325 : 1.0;       // atm → bar
  EnergyOptions e = e0;
  const Layout L = prepare(s, ff, e, st);
  if (notes) notes->insert(notes->end(), L.notes.begin(), L.notes.end());
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  char b[400];
  out << "# LAMMPS input written by CAPS: " << clean_title(s.title, ff.name) << " · " << ff.name << "\n";
  if (st.native)
    out << "# " << ff.name << " in its own LAMMPS styles" << (st.hybrid ? " (hybrid form)" : "") << "\n\n";
  else
    out << "# the same force field and cut-offs CAPS uses (energies and forces checked against LAMMPS: bench/ff/check_data_lammps.py)\n\n";
  // a structure without a cell sits in a 100 Å box (as in the data file): periodic, but too large for images to interact
  if (metal)
    out << "# units metal: energies in eV (every parameter of the force field divided by 23.060549, LAMMPS's own factor), time in ps,\n"
           "# pressure in bar" << (ff.manybody.on() ? ", as LAMMPS reads " + ff.manybody.style + "'s file" : std::string()) << "\n";
  out << (metal ? "units           metal\n" : "units           real\n") << "atom_style      full\nboundary        p p p\n";
  {
    char t[64];   // in fs (ps in metal units), with the other settings at the top, as force-field input files give it
    std::snprintf(t, sizeof t, "timestep        %.6g\n\n", lammps_timestep(run, ff) * tu);
    out << t;
  }
  auto aligned = [](const std::string& l) {   // "keyword       arguments", as the rest of the script
    const size_t sp = l.find(' ');
    if (sp == std::string::npos || sp >= 16) return l;
    return l.substr(0, sp) + std::string(16 - sp, ' ') + l.substr(sp + 1);
  };
  // the long-range solver after read_data: LAMMPS sets PPPM up for the box it is defined with (a triclinic cell needs it
  // defined after the cell is read)
  std::vector<std::string> kspace;
  for (const auto& l : style_lines(L, ff, e)) {
    if (l.rfind("kspace_style", 0) == 0) kspace.push_back(l);
    else out << aligned(l) << "\n";
  }
  out << "\nread_data       " << data_path << "\n";
  for (const auto& l : kspace) out << aligned(l) << "\n";
  if (pair_coeffs) {
    out << "\n# pair coefficients: every type pair, " << ff.mixing << " mixing applied by CAPS (nothing left to LAMMPS's mixing)\n";
    for (const auto& l : pair_lines(L, ff)) out << "pair_coeff      " << l << "\n";
  }
  if (L.hbond) {
    out << "\n# DREIDING hydrogen bonds: donor-acceptor type pairs, the hydrogen type\n";
    for (const auto& l : hbond_lines(ff)) out << "pair_coeff      " << l << "\n";
  }
  if (!L.mb_types.empty())
    out << "\n# " << ff.manybody.style << " (" << ff.manybody.file.substr(ff.manybody.file.find_last_of('/') + 1) << (ff.manybody.citation.empty() ? "" : "; " + ff.manybody.citation)
        << ") for its elements, Lennard-Jones between them and the rest; "
        << (metal ? (ff.manybody.units == "metal" ? "the file in metal units, as this input" : "LAMMPS converts the file from real to metal units")
                  : (ff.manybody.units == "metal" ? "LAMMPS converts the file from metal to real units" : "the file is in real units")) << "\n";
  if (!L.mb_types.empty() && !L.bonds.term_type.empty())
    out << "# (LAMMPS warns of a many-body potential beside bonds: its atoms have none, the special_bonds exclusions act on the other groups only)\n";
  for (const auto& l : after_read(L, e, data_path, ff.excluded_type_pairs, &ff)) out << aligned(l) << "\n";
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
  {
    std::snprintf(b, sizeof b, "\n# 1. energy minimisation\nmin_style       cg\nminimize        1.0e-4 %s 5000 50000\nreset_timestep  0\n", metal ? "4.34e-8" : "1.0e-6");
    out << b;   // the force tolerance 1e-6 kcal/mol/Å, in eV/Å in metal units
  }
  if (run.kind == K::Minimize) {
    out << "\nwrite_data      minimized.data\n";
    return;
  }
  const bool npt = run.kind == K::NPT;
  out << "\n# 2. " << (npt ? "NPT" : "NVT") << " molecular dynamics (Nosé–Hoover)\n";
  if (run.constraints != ConstraintMode::None) {
    const std::string line = shake_fix(s, ff, L, run.constraints, mobile);
    if (!line.empty())
      out << "# " << (run.constraints == ConstraintMode::AllBonds ? "every bond" : "bonds to hydrogen (and water's angle)") << " held at its length, as in CAPS (before the velocities: the temperature counts the constraints)\n" << line;
  }
  if (L.coreshell) {
    // a shell model (LAMMPS CORESHELL): the thermostat sees the ions' centre-of-mass motion, not the core-shell vibration
    std::set<int> cores, shells;
    for (const auto& bd : s.bonds) {
      const Vec3 d = s.cell.valid() ? s.cell.minimum_image(s.atoms[bd.j].pos - s.atoms[bd.i].pos) : s.atoms[bd.j].pos - s.atoms[bd.i].pos;
      if (dot(d, d) >= 0.25 * 0.25) continue;
      const bool j_shell = ff.mass[bd.j] < ff.mass[bd.i];
      shells.insert(ff.type_index[j_shell ? bd.j : bd.i] + 1);
      cores.insert(ff.type_index[j_shell ? bd.i : bd.j] + 1);
    }
    std::string gc, gs;
    for (int t : cores) gc += " " + std::to_string(t);
    for (int t : shells) gs += " " + std::to_string(t);
    out << "group           cores type" << gc << "\ngroup           shells type" << gs << "\n"
        << "comm_modify     vel yes\ncompute         CSequ all temp/cs cores shells\nthermo_modify   temp CSequ\n";
    std::snprintf(b, sizeof b, "velocity        %s create %.6g %llu dist gaussian mom yes rot no bias yes temp CSequ\n", mobile.c_str(), run.temperature,
                  static_cast<unsigned long long>(run.seed));
  } else {
    std::snprintf(b, sizeof b, "velocity        %s create %.6g %llu mom yes rot yes dist gaussian\n", mobile.c_str(), run.temperature,
                  static_cast<unsigned long long>(run.seed));
  }
  out << b;
  if (npt)
    std::snprintf(b, sizeof b, "fix             integrate %s npt temp %.6g %.6g %.6g iso %.6g %.6g %.6g\n", mobile.c_str(), run.temperature, run.temperature,
                  run.tdamp * tu, run.pressure * pu, run.pressure * pu, run.pdamp * tu);
  else
    std::snprintf(b, sizeof b, "fix             integrate %s nvt temp %.6g %.6g %.6g\n", mobile.c_str(), run.temperature, run.temperature, run.tdamp * tu);
  out << b;
  if (L.coreshell) out << "fix_modify      integrate temp CSequ\n";

  std::snprintf(b, sizeof b, "dump            traj all custom %d traj.lammpstrj id mol type q xu yu zu\ndump_modify     traj sort id\nrun             %lld\n",
                std::max(1, run.dump_every), static_cast<long long>(run.steps));
  out << b;
  out << "\nwrite_data      final.data\nwrite_restart   final.restart\n";
}

double lammps_timestep(const LammpsRun& run, const ForceField& ff) {
  if (run.dt > 0) return run.dt;
  if (ff.native_timestep > 0) return ff.native_timestep;
  return run.constraints != ConstraintMode::None ? 2.0 : 0.5;
}

}  // namespace caps
