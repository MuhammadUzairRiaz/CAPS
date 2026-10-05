#include "caps/water.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace caps {

namespace {
constexpr double kBoltzmannKcal = 0.0019872043;   // kcal/mol/K (CODATA 2018 k_B N_A)
constexpr double kPi = 3.14159265358979323846;

std::string lower(std::string s) {
  for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
}  // namespace

const std::vector<WaterModel>& water_models() {
  // constants for a flexible run of a rigid model: LAMMPS Howto_tip3p's (TIP3P and TIP4P families), moltemplate's
  // spce.lt (SPC, SPC/E); only SPC/Fw's are the model's own (Wu, Tepper & Voth: k/2 (r − r0)², so K = k/2)
  const std::string tip_k = "LAMMPS Howto_tip3p (K 450 kcal/mol/Å², 55 kcal/mol/rad²); the model is rigid";
  const std::string spc_k = "moltemplate spce.lt (K 600 kcal/mol/Å², 75 kcal/mol/rad²); the model is rigid";
  static const std::vector<WaterModel> m = {
      {"spc", "SPC", "H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, J. Hermans, in Intermolecular Forces (Reidel, 1981), p. 331",
       "berendsen1981", 3, 1.0, 109.47, 0.41, -0.82, 0, 0.1553, 3.166, 0, 0, true, 600, 75, spc_k, ""},
      {"spce", "SPC/E", "H. J. C. Berendsen, J. R. Grigera, T. P. Straatsma, J. Phys. Chem. 91, 6269 (1987)", "berendsen1987", 3, 1.0, 109.47,
       0.4238, -0.8476, 0, 0.1553, 3.166, 0, 0, true, 600, 75, spc_k, ""},
      {"spcfw", "SPC/Fw", "Y. Wu, H. L. Tepper, G. A. Voth, J. Chem. Phys. 124, 024503 (2006)", "wu2006", 3, 1.012, 113.24, 0.41, -0.82, 0,
       0.1554253, 3.165492, 0, 0, false, 1059.162 / 2, 75.90 / 2,
       "the model's own: k_b 1059.162 kcal/mol/Å², k_θ 75.90 kcal/mol/rad² for k/2 (x − x0)²", "flexible: run without constraints"},
      {"tip3p", "TIP3P", "W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, M. L. Klein, J. Chem. Phys. 79, 926 (1983)", "jorgensen1983",
       3, 0.9572, 104.52, 0.417, -0.834, 0, 0.1521, 3.1507, 0, 0, true, 450, 55, tip_k, ""},
      {"tip3p-charmm", "TIP3P (CHARMM)", "A. D. MacKerell Jr. et al., J. Phys. Chem. B 102, 3586 (1998)", "mackerell1998", 3, 0.9572, 104.52,
       0.417, -0.834, 0, 0.1521, 3.1507, 0.0460, 0.4, true, 450, 55, tip_k, "hydrogens with Lennard-Jones (ε 0.046 kcal/mol, σ 0.4 Å), as CHARMM uses it"},
      {"tip3p-ew", "TIP3P-Ew", "D. J. Price, C. L. Brooks III, J. Chem. Phys. 121, 10096 (2004)", "price2004", 3, 0.9572, 104.52, 0.415,
       -0.830, 0, 0.102, 3.188, 0, 0, true, 450, 55, tip_k,
       "charges −0.830 / 0.415 as moltemplate's tip3p_2004.lt gives them (LAMMPS's Howto_tip3p lists the TIP3P charges −0.834 / 0.417 for it)"},
      {"tip4p", "TIP4P", "W. L. Jorgensen, J. Chandrasekhar, J. D. Madura, R. W. Impey, M. L. Klein, J. Chem. Phys. 79, 926 (1983)", "jorgensen1983",
       4, 0.9572, 104.52, 0.520, -1.040, 0.15, 0.1550, 3.1536, 0, 0, true, 450, 55, tip_k, ""},
      {"tip4p-ew", "TIP4P-Ew", "H. W. Horn, W. C. Swope, J. W. Pitera, J. D. Madura, T. J. Dick, G. L. Hura, T. Head-Gordon, J. Chem. Phys. 120, 9665 (2004)",
       "horn2004", 4, 0.9572, 104.52, 0.52422, -1.04844, 0.125, 0.16275, 3.16435, 0, 0, true, 450, 55, tip_k, ""},
      {"tip4p2005", "TIP4P/2005", "J. L. F. Abascal, C. Vega, J. Chem. Phys. 123, 234505 (2005)", "abascal2005", 4, 0.9572, 104.52, 0.5564,
       -1.1128, 0.1546, 93.2 * kBoltzmannKcal, 3.1589, 0, 0, true, 450, 55, tip_k, "ε/k_B = 93.2 K"},
      {"tip4p-ice", "TIP4P/Ice", "J. L. F. Abascal, E. Sanz, R. García Fernández, C. Vega, J. Chem. Phys. 122, 234511 (2005)", "abascal2005ice", 4,
       0.9572, 104.52, 0.5897, -1.1794, 0.1577, 106.1 * kBoltzmannKcal, 3.1668, 0, 0, true, 450, 55, tip_k, "ε/k_B = 106.1 K; made for ice and the melting point"},
      {"tip5p", "TIP5P", "M. W. Mahoney, W. L. Jorgensen, J. Chem. Phys. 112, 8910 (2000)", "mahoney2000", 5, 0.9572, 104.52, 0.241, -0.241, 0,
       0.16, 3.12, 0, 0, true, 450, 55, tip_k, "two lone pairs 0.70 Å from O at 109.47°; GROMACS only (LAMMPS has no five-site water)", 0.70, 109.47},
      {"opc", "OPC", "S. Izadi, R. Anandakrishnan, A. V. Onufriev, J. Phys. Chem. Lett. 5, 3863 (2014)", "izadi2014", 4, 0.8724, 103.6, 0.6791,
       -1.3582, 0.1594, 0.2128, 3.1660, 0, 0, true, 450, 55, tip_k, ""},
  };
  return m;
}

const WaterModel& water_model(const std::string& id_or_name) {
  const std::string k = lower(id_or_name);
  for (const auto& m : water_models())
    if (lower(m.id) == k || lower(m.name) == k) return m;
  std::string known;
  for (const auto& m : water_models()) known += (known.empty() ? "" : ", ") + m.name;
  throw std::invalid_argument("unknown water model '" + id_or_name + "' (" + known + ")");
}

double water_m_alpha(const WaterModel& m) {
  if (m.sites != 4) return 0;
  return m.d_om / (m.r_oh * std::cos(m.theta * kPi / 360.0));
}

std::array<double, 2> water_lp_coefficients(const WaterModel& m) {
  if (m.sites != 5) return {0, 0};
  const double ht = m.theta * kPi / 360.0, hl = m.theta_l * kPi / 360.0;
  return {-m.d_ol * std::cos(hl) / (2 * m.r_oh * std::cos(ht)), m.d_ol * std::sin(hl) / (m.r_oh * m.r_oh * std::sin(2 * ht))};
}

std::vector<std::array<int64_t, 5>> find_waters(const System& s) {
  const size_t n = s.atoms.size();
  std::vector<std::vector<uint32_t>> nb(n);
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  int nmol = 0;
  const auto mol = s.molecules(&nmol);
  std::vector<std::vector<uint32_t>> members(size_t(std::max(nmol, 0)));
  for (size_t i = 0; i < n; ++i) if (mol[i] >= 0 && mol[i] < nmol) members[size_t(mol[i])].push_back(uint32_t(i));
  std::vector<std::array<int64_t, 5>> out;
  for (const auto& mem : members) {
    if (mem.size() < 3 || mem.size() > 5) continue;
    int64_t o = -1;
    std::vector<int64_t> h, x;
    for (auto i : mem) {
      const int el = s.atoms[i].element;
      if (el == 8 && o < 0) o = i;
      else if (el == 1) h.push_back(i);
      else if (el == 0) x.push_back(i);
      else { o = -2; break; }
    }
    if (o < 0 || h.size() != 2 || mem.size() != 3 + x.size()) continue;
    auto bonded = [&](int64_t a, int64_t b) { return std::find(nb[size_t(a)].begin(), nb[size_t(a)].end(), uint32_t(b)) != nb[size_t(a)].end(); };
    if (!bonded(o, h[0]) || !bonded(o, h[1])) continue;
    out.push_back({o, h[0], h[1], x.size() > 0 ? x[0] : -1, x.size() > 1 ? x[1] : -1});
  }
  std::sort(out.begin(), out.end());
  return out;
}

size_t apply_water_model(System& s, const WaterModel& m, std::vector<std::string>* notes) {
  if (s.topology) throw std::invalid_argument("a structure with a term-by-term topology keeps its own water: set the model in its topology file");
  const auto waters = find_waters(s);
  if (waters.empty()) return 0;
  if (!s.has_mol) {   // molecule ids from the bonds, so an M site (bonded to nothing) stays in its water's molecule
    int nm = 0;
    const auto mol = s.molecules(&nm);
    for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].mol = int64_t(mol[i]) + 1;
    s.has_mol = true;
  }
  const double half = m.theta * kPi / 360.0;
  // geometry: O, the plane and the bisector kept
  for (const auto& w : waters) {
    const Vec3 o = s.atoms[size_t(w[0])].pos;
    Vec3 d1 = s.atoms[size_t(w[1])].pos - o, d2 = s.atoms[size_t(w[2])].pos - o;
    if (s.cell.valid()) d1 = s.cell.minimum_image(d1), d2 = s.cell.minimum_image(d2);
    Vec3 u = d1 * (1.0 / std::max(1e-12, norm(d1))) + d2 * (1.0 / std::max(1e-12, norm(d2)));
    if (norm(u) < 1e-8) u = Vec3{0, 0, 1};
    u = u * (1.0 / norm(u));
    Vec3 v = d1 - u * dot(d1, u);
    if (norm(v) < 1e-8) v = std::fabs(u[0]) < 0.9 ? cross(u, Vec3{1, 0, 0}) : cross(u, Vec3{0, 1, 0});
    v = v * (1.0 / norm(v));
    s.atoms[size_t(w[1])].pos = o + (u * std::cos(half) + v * std::sin(half)) * m.r_oh;
    s.atoms[size_t(w[2])].pos = o + (u * std::cos(half) - v * std::sin(half)) * m.r_oh;
    s.atoms[size_t(w[0])].charge = m.sites >= 4 ? 0.0 : m.q_neg;
    s.atoms[size_t(w[1])].charge = s.atoms[size_t(w[2])].charge = m.q_h;
    for (int k : {3, 4}) if (w[size_t(k)] >= 0) s.atoms[size_t(w[size_t(k)])].charge = m.q_neg;
    s.atoms[size_t(w[0])].name = "OW", s.atoms[size_t(w[1])].name = "HW1", s.atoms[size_t(w[2])].name = "HW2";
  }
  s.has_charges = true;
  // extra sites: M (4 sites) or two lone pairs (5 sites) added after each water's last hydrogen, the old ones removed
  const int want_x = m.sites - 3;
  const bool want_m = want_x > 0;
  const bool have_m = std::any_of(waters.begin(), waters.end(), [](const auto& w) { return w[3] >= 0; });
  size_t added = 0, removed = 0;
  if (want_m || have_m) {
    std::vector<char> is_m(s.atoms.size(), 0);
    std::vector<int64_t> after(s.atoms.size(), -1);   // water index whose M follows this atom
    for (size_t k = 0; k < waters.size(); ++k) {
      for (int j : {3, 4}) if (waters[k][size_t(j)] >= 0) is_m[size_t(waters[k][size_t(j)])] = 1;
      after[size_t(std::max({waters[k][0], waters[k][1], waters[k][2]}))] = int64_t(k);
    }
    int mtype = 0;
    if (want_m && !s.types.empty()) {
      const std::string xl = want_x == 2 ? "LP" : "MW";
      for (const auto& t : s.types) if (t.label == xl) mtype = t.type;
      if (!mtype) {
        for (const auto& t : s.types) mtype = std::max(mtype, t.type);
        s.types.push_back({++mtype, 0.0, xl});
      }
    }
    const double alpha = water_m_alpha(m);
    const auto lp = water_lp_coefficients(m);
    std::vector<Atom> atoms;
    std::vector<int64_t> map(s.atoms.size(), -1);
    std::vector<Vec3> vel;
    const bool has_v = s.velocities.size() == s.atoms.size();
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      if (is_m[i]) { ++removed; continue; }
      map[i] = int64_t(atoms.size());
      atoms.push_back(s.atoms[i]);
      if (has_v) vel.push_back(s.velocities[i]);
      if (want_m && after[i] >= 0) {
        const auto& w = waters[size_t(after[i])];
        const Vec3 o = s.atoms[size_t(w[0])].pos, h1 = s.atoms[size_t(w[1])].pos, h2 = s.atoms[size_t(w[2])].pos;
        Atom a = s.atoms[size_t(w[0])];
        a.element = 0, a.charge = m.q_neg, a.type = mtype;
        a.id = 0;
        if (want_x == 1) {
          a.name = "MW";
          a.pos = o * (1 - alpha) + (h1 + h2) * (alpha / 2);
          atoms.push_back(a);
        } else {
          Vec3 r1 = h1 - o, r2 = h2 - o;
          if (s.cell.valid()) r1 = s.cell.minimum_image(r1), r2 = s.cell.minimum_image(r2);
          for (int sgn : {1, -1}) {
            a.name = sgn > 0 ? "LP1" : "LP2";
            a.pos = o + (r1 + r2) * lp[0] + cross(r1, r2) * (sgn * lp[1]);
            atoms.push_back(a);
          }
        }
        for (int j = 0; j < want_x; ++j)
          if (has_v) vel.push_back(s.velocities[size_t(w[0])]);
        added += size_t(want_x);
      }
    }
    std::vector<Bond> bonds;
    for (const auto& b : s.bonds)
      if (map[b.i] >= 0 && map[b.j] >= 0) bonds.push_back({uint32_t(map[b.i]), uint32_t(map[b.j]), b.order});
    s.atoms = std::move(atoms);
    s.bonds = std::move(bonds);
    if (has_v) s.velocities = std::move(vel);
    for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].id = int64_t(i) + 1;
  }
  if (notes) {
    notes->push_back(std::to_string(waters.size()) + " waters as " + m.name + " (" + m.citation + ")" +
                     (added ? (want_x == 2 ? ", two lone pairs added to each" : ", an M site added to each") : removed ? ", their extra sites removed" : ""));
    if (!m.note.empty()) notes->push_back(m.name + ": " + m.note);
  }
  return waters.size();
}

ForceField water_forcefield(const System& s, const WaterModel& m, const std::vector<uint32_t>& atoms) {
  const auto waters = find_waters(s);
  std::vector<int64_t> local(s.atoms.size(), -1);
  for (size_t k = 0; k < atoms.size(); ++k) local[atoms[k]] = int64_t(k);
  ForceField F;
  F.name = m.name + " water";
  F.type_names = {"OW", "HW"};
  F.lj = {{m.eps_o, m.sigma_o}, {m.eps_h, m.sigma_h}};
  if (m.sites == 4) F.type_names.push_back("MW"), F.lj.push_back({0, 0});
  if (m.sites == 5) F.type_names.push_back("LP"), F.lj.push_back({0, 0});
  const auto lp = water_lp_coefficients(m);
  const size_t n = atoms.size();
  F.atom_type.assign(n, ""), F.why.assign(n, ""), F.type_index.assign(n, -1), F.charge.assign(n, 0), F.mass.assign(n, 0);
  F.excluded.assign(n, {});
  F.pair_form = "lj12-6", F.mixing = "arithmetic";
  F.lj14 = 0.5, F.coul14 = 1.0 / 1.2;   // no 1-4 pairs in a water: these are never used
  const double alpha = water_m_alpha(m);
  for (const auto& w : waters) {
    int64_t o = local[size_t(w[0])], h1 = local[size_t(w[1])], h2 = local[size_t(w[2])], mm = w[3] >= 0 ? local[size_t(w[3])] : -1,
            l2 = w[4] >= 0 ? local[size_t(w[4])] : -1;
    if (o < 0 && h1 < 0 && h2 < 0 && mm < 0 && l2 < 0) continue;
    const int have = (mm >= 0) + (l2 >= 0);
    if (o < 0 || h1 < 0 || h2 < 0 || have != m.sites - 3)
      throw FieldError(m.name + ": a water is split between groups, or " +
                       (have < m.sites - 3 ? "lacks its extra sites (apply the model first)" : "has extra sites this model does not"));
    auto set = [&](int64_t i, int t, double q, double mass, const std::string& why) {
      F.type_index[size_t(i)] = t, F.atom_type[size_t(i)] = F.type_names[size_t(t)], F.charge[size_t(i)] = q, F.mass[size_t(i)] = mass, F.why[size_t(i)] = why;
    };
    set(o, 0, m.sites >= 4 ? 0.0 : m.q_neg, 15.9994, m.name + " oxygen");
    set(h1, 1, m.q_h, 1.008, m.name + " hydrogen");
    set(h2, 1, m.q_h, 1.008, m.name + " hydrogen");
    F.bonds.push_back({uint32_t(o), uint32_t(h1), m.k_bond, m.r_oh});
    F.bonds.push_back({uint32_t(o), uint32_t(h2), m.k_bond, m.r_oh});
    F.angles.push_back({uint32_t(h1), uint32_t(o), uint32_t(h2), m.k_angle, m.theta * kPi / 180.0});   // θ0 in radians
    std::vector<uint32_t> all{uint32_t(o), uint32_t(h1), uint32_t(h2)};
    if (m.sites == 4) {
      set(mm, 2, m.q_neg, 0.0, m.name + " M site (" + std::to_string(m.d_om) + " Å from O)");
      F.vsites.push_back({uint32_t(mm), {uint32_t(o), uint32_t(h1), uint32_t(h2)}, {1 - alpha, alpha / 2, alpha / 2}});
      all.push_back(uint32_t(mm));
    }
    if (m.sites == 5) {   // the two lone pairs, on either side of the plane
      for (auto [l, sgn] : {std::pair<int64_t, double>{mm, 1.0}, {l2, -1.0}}) {
        set(l, 2, m.q_neg, 0.0, m.name + " lone pair (" + std::to_string(m.d_ol) + " Å from O)");
        F.vsites.push_back({uint32_t(l), {uint32_t(o), uint32_t(h1), uint32_t(h2)}, {1 - 2 * lp[0], lp[0], lp[0]}, sgn * lp[1]});
        all.push_back(uint32_t(l));
      }
    }
    std::sort(all.begin(), all.end());
    for (auto i : all)
      for (auto j : all) if (i != j) F.excluded[i].push_back(j);
  }
  for (size_t k = 0; k < n; ++k)
    if (F.type_index[k] < 0) throw FieldError(m.name + ": atom " + std::to_string(atoms[k] + 1) + " is not part of a water");
  F.notes.push_back(m.name + " (" + m.citation + ")");
  if (m.rigid) F.notes.push_back(m.name + " is rigid: run it with its bonds and angle constrained (SHAKE / RATTLE); the bond and angle constants (" + m.k_source + ") serve a flexible run only");
  if (!m.note.empty()) F.notes.push_back(m.name + ": " + m.note);
  return F;
}

}  // namespace caps
