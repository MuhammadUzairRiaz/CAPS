// GROMACS files carrying a CAPS force field: a topology (.top) with every term in the GROMACS function that has the
// same energy, the coordinates (.gro) and the run parameters (.mdp) for the same cut-offs. Units: kcal/mol → kJ/mol,
// Å → nm; CAPS's harmonic k (r − r0)² is GROMACS's ½ kb (r − b0)², so kb = 2k.
//
// What makes the energies agree term by term:
//   · every i-j Lennard-Jones pair is written in [nonbond_params] with CAPS's mixing (and overrides) applied, so no
//     GROMACS combination rule is used;
//   · every 1-4 pair is written in [pairs] with its scaled σ, ε (so fudgeLJ is 1 and CHARMM's separate 1-4 parameters
//     are kept); fudgeQQ is CAPS's 1-4 Coulomb scale;
//   · nrexcl is 0 and [exclusions] lists CAPS's excluded partners exactly (1-2, 1-3 and 1-4 pairs);
//   · tail corrections: DispCorr AllEnerPres (dispersion and repulsion), else the LJ potential is shifted at the
//     cut-off (vdw-modifier Potential-shift), as CAPS.
// Electrostatics: GROMACS has no damped shifted force; periodic cells get PME with CAPS's ewald-rtol, which is what CAPS's
// PME computes and what DSF approximates. Without a cell the molecule sits in a box 2 r_c larger than it is, with
// plain cut-off Coulomb (reaction field, ε_rf = 1).
//
// Checked against CAPS with gmx grompp, mdrun -rerun and gmx energy (bench/ff/check_gromacs.py).
#include <cctype>
#include <cstdio>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <map>
#include <queue>
#include <set>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/relax.hpp"

namespace caps {

namespace {

constexpr double R2D = 57.29577951308232;
constexpr double KJ = 4.184;

std::string fmt(const char* f, double a) {
  char b[64];
  std::snprintf(b, sizeof b, f, a);
  return b;
}

// GROMACS names: no spaces, at most 16 characters
std::string clean(std::string s) {
  for (auto& c : s)
    if (c == ' ' || c == ';' || c == '[' || c == ']') c = '_';
  if (s.empty()) s = "X";
  return s.substr(0, 16);
}

// A type name as GROMACS writes it: moltemplate's OPLS-AA names carry their bonded classes (135_bCT_aCT_dCT_iCT), which
// the itp's bonded lines never need, so the number is kept, spelt as GROMACS's own oplsaa.ff does (opls_135); a name
// starting with a digit gets a letter in front, as grompp reads such a name as a number.
std::string type_label(const std::string& full) {
  std::string s = full;
  const size_t k = s.find("_b");
  if (k != std::string::npos && k > 0 && s.find("_a", k) != std::string::npos && s.find("_d", k) != std::string::npos) s = s.substr(0, k);
  s = clean(s);
  if (!s.empty() && std::isdigit(static_cast<unsigned char>(s[0]))) s = (s.size() <= 11 ? "opls_" : "t") + s;
  return clean(s);
}

// An atom name for the itp and gro (5 characters): the structure's own name when it is one (a PDB name such as CA or
// HB2; not a force-field type a typed file carries as the name), otherwise the element with its count in the molecule
// (C1, H14; C* past 9999)
std::string atom_label(const Atom& a, const std::string& type, std::map<int, int>& count) {
  const std::string sym = element(a.element).symbol;
  const std::string& nm = a.name;
  if (!nm.empty() && nm != type && nm.size() <= 5 && nm.find('_') == std::string::npos && std::isalpha(static_cast<unsigned char>(nm[0])))
    return clean(nm);
  const int k = ++count[a.element];
  const std::string num = std::to_string(k);
  return sym.size() + num.size() <= 5 ? sym + num : sym + "*";
}

bool periodic_bonds(const System& s, std::vector<Vec3>& pos) {
  // makes molecules whole across the cell; true when a bond still spans it afterwards (an infinite network)
  const size_t n = s.atoms.size();
  pos.resize(n);
  for (size_t i = 0; i < n; ++i) pos[i] = s.atoms[i].pos - s.cell.origin;
  if (!s.cell.valid()) return false;
  const auto nb = s.neighbours();
  std::vector<char> seen(n, 0);
  for (size_t r = 0; r < n; ++r) {
    if (seen[r]) continue;
    seen[r] = 1;
    std::queue<uint32_t> q;
    q.push(uint32_t(r));
    while (!q.empty()) {
      const uint32_t a = q.front();
      q.pop();
      for (uint32_t b : nb[a])
        if (!seen[b]) {
          seen[b] = 1;
          pos[b] = pos[a] + s.cell.minimum_image(pos[b] - pos[a]);
          q.push(b);
        }
    }
  }
  for (const auto& b : s.bonds) {
    const Vec3 d = pos[b.j] - pos[b.i], m = s.cell.minimum_image(d);
    if (norm(d - m) > 1e-6) return true;
  }
  return false;
}

// MARTINI 2's pairs as LAMMPS lj/gromacs: Lennard-Jones with GROMACS's force switch from lj_inner to the cut-off
bool gromacs_switch(const ForceField& ff) {
  if (ff.pair_func.empty()) return false;
  for (const auto& [k, p] : ff.pair_func)
    if (p.form != kPairGromacs) return false;
  return ff.lj_inner > 0;
}

}  // namespace

std::vector<std::string> gromacs_notes(const System& s, const ForceField& ff, const EnergyOptions& e) {
  std::vector<std::string> notes;
  const size_t n = s.atoms.size();
  if (ff.atom_type.size() != n) throw FieldError("the force field does not cover every atom");
  if (ff.pair_form != "lj12-6") throw FieldError(ff.name + ": the 9-6 Lennard-Jones form (class II) has no GROMACS function");
  // lj/gromacs (MARTINI 2's shifted Lennard-Jones) is GROMACS's own force switch; the other pair forms have no Verlet form
  for (const auto& [k, p] : ff.pair_func)
    if (p.form != kPairGromacs)
      throw FieldError(ff.name + ": " + (p.form == 1 ? std::string("Buckingham") : p.form == 2 ? std::string("Morse") : p.form >= kPairSdk96 && p.form <= kPairSdk125 ? std::string("SDK / SPICA") : std::string("cosine-squared")) +
                       " pairs have no GROMACS form in the Verlet scheme; export to LAMMPS instead");
  if (gromacs_switch(ff)) {
    bool charged = false;
    for (double q : ff.charge) charged = charged || q != 0;
    if (charged && ff.coul_gromacs && e.coulomb)
      throw FieldError(ff.name + ": MARTINI 2's shifted Coulomb (ε_r " + fmt("%g", ff.dielectric) +
                       ") has no form in GROMACS's Verlet scheme (the group scheme that had it is gone): export this charged system to LAMMPS, which keeps it exactly");
    notes.push_back("Lennard-Jones with GROMACS's force switch from " + fmt("%g", ff.lj_inner / 10) + " nm to the cut-off (vdw-modifier Force-switch), MARTINI 2's lj/gromacs");
  }
  if (ff.sw.on) throw FieldError(ff.name + ": the Stillinger–Weber three-body term (mW water) has no GROMACS form");
  if (ff.manybody.on())
    throw FieldError(ff.name + ": GROMACS has no " + ff.manybody.style + " (or any other many-body potential read from a file); export this system to LAMMPS, or give the crystal a force field of pair terms (IFF, INTERFACE, UFF)");
  if (!ff.bonds2.empty() || !ff.angles2.empty() || !ff.dihedrals2.empty() || !ff.impropers2.empty())
    throw FieldError(ff.name + ": class II terms (COMPASS, PCFF) have no GROMACS functions; export to LAMMPS instead");
  if (!ff.inversions.empty()) throw FieldError(ff.name + ": inversion (umbrella) terms (DREIDING, UFF) have no GROMACS function; export to LAMMPS instead");
  if (ff.hbond.on()) throw FieldError(ff.name + ": DREIDING's hydrogen-bond term (hbond/dreiding/lj) has no GROMACS function; export to LAMMPS instead");
  if (ff.lj_fsw)
    notes.push_back("Lennard-Jones with GROMACS's force switch from " + std::to_string(ff.lj_inner / 10).substr(0, 4) + " nm (CHARMM36's GROMACS setting); "
                    "its polynomial differs slightly from CHARMM's own switch (CAPS, LAMMPS lj/charmmfsw), and GROMACS's 1-4 pairs carry "
                    "no switch offset: about 0.1 % of the Lennard-Jones energy");
  for (const auto& b : ff.bonds_x)
    if (b.form != 1 && b.form != 2) throw FieldError("bond form " + std::to_string(b.form) + " has no GROMACS function");
  for (const auto& a : ff.angles_x)
    if (a.form != 1 && a.form != 5) throw FieldError("angle form " + std::to_string(a.form) + " has no GROMACS function");
  for (const auto& t : ff.dihedrals)
    if (t.n < 0) throw FieldError("a torsion with multiplicity < 0 has no GROMACS form");
  for (const auto& t : ff.impropers)
    if (t.n < 0) throw FieldError("an improper with multiplicity < 0 has no GROMACS form");
  std::vector<Vec3> pos;
  if (periodic_bonds(s, pos)) notes.push_back("bonds cross the cell (an infinite network): periodic-molecules = yes");
  if (!s.cell.valid()) notes.push_back("no periodic cell: the molecule is centred in a box 2 r_c larger than it is");
  const bool any_charge = std::any_of(ff.charge.begin(), ff.charge.end(), [](double q) { return q != 0; });
  if (!any_charge) {
  } else if (ff.coul_rf && e.coulomb) notes.push_back("reaction-field Coulomb (ε_r " + fmt("%g", ff.dielectric) + ", ε_rf " + (ff.eps_rf == 0 ? std::string("∞") : fmt("%g", ff.eps_rf)) + "), as in CAPS");
  else if (s.cell.valid() && e.coulomb && e.electrostatics != EnergyOptions::Electrostatics::PME)
    notes.push_back("CAPS's damped shifted force becomes PME in GROMACS (no DSF there)");
  if (s.cell.valid()) {
    const Cell& c = s.cell;
    const double w = std::min({c.a[0], c.b[1], c.c[2]});
    if (w < 2 * e.cutoff)
      notes.push_back("the cell is " + fmt("%.1f", w) + " Å across, less than 2 r_c (" + fmt("%.0f", 2 * e.cutoff) +
                      " Å): GROMACS needs a larger cell (a supercell) or a shorter cut-off; CAPS sums the further images");
  }
  if (!s.cell.valid() && e.coulomb && !ff.coul_rf) notes.push_back("no cell: plain cut-off Coulomb in GROMACS, damped shifted force in CAPS");
  if (!e.tail) notes.push_back("1-4 Lennard-Jones is unshifted in GROMACS: a constant offset from CAPS's shifted 1-4 terms, no force difference");
  return notes;
}

std::string gromacs_mdp(const System& s, const ForceField& ff, const EnergyOptions& e0) {
  EnergyOptions e = e0;
  if (std::all_of(ff.charge.begin(), ff.charge.end(), [](double q) { return q == 0; })) e.coulomb = false;   // no charges: no Coulomb term
  if (ff.cutoff > 0) e.cutoff = ff.cutoff;
  const bool cell = s.cell.valid();
  const bool pme = cell && e.coulomb;
  std::ostringstream m;
  m << "; non-bonded settings matching CAPS ("
    << (e.coulomb ? (ff.coul_rf ? "reaction field" : e.electrostatics == EnergyOptions::Electrostatics::PME ? "PME" : "damped shifted force") : "no Coulomb")
    << ", r_c " << e.cutoff << " Å)\n";
  m << "cutoff-scheme            = Verlet\n";
  m << "pbc                      = xyz\n";
  m << "rvdw                     = " << e.cutoff / 10 << "\n";
  m << "rcoulomb                 = " << e.cutoff / 10 << "\n";
  m << "vdwtype                  = Cut-off\n";
  const bool tail = e.tail && !ff.lj_shift && !ff.lj_fsw && !gromacs_switch(ff);
  if (ff.lj_fsw || gromacs_switch(ff)) {   // CHARMM (as CHARMM-GUI writes CHARMM36) or MARTINI 2's lj/gromacs: the force switch from lj_inner
    m << "vdw-modifier             = Force-switch\n";
    m << "rvdw-switch              = " << ff.lj_inner / 10 << "\n";
  } else {
    m << "vdw-modifier             = " << (tail ? "None" : "Potential-shift") << "\n";
  }
  m << "DispCorr                 = " << (tail && cell ? "AllEnerPres" : "no") << "\n";
  if (ff.dielectric != 1 && e.coulomb) m << "epsilon-r                = " << ff.dielectric << "\n";
  if (e.coulomb && ff.coul_rf) {   // Martini 3: reaction field, ε_rf 0 meaning infinite
    m << "coulombtype              = Reaction-Field\n";
    m << "epsilon-rf               = " << ff.eps_rf << "\n";
  } else if (!e.coulomb) {
    m << "coulombtype              = Cut-off       ; no charges: no Coulomb energy\n";
  } else if (pme) {
    if (e.electrostatics != EnergyOptions::Electrostatics::PME)
      m << "; CAPS runs this with damped shifted force (α " << e.dsf_alpha << " Å⁻¹); GROMACS has no DSF, so PME: the Ewald sum DSF approximates\n";
    m << "coulombtype              = PME\n";
    m << "ewald-rtol               = " << fmt("%.3g", e.ewald_rtol) << "\n";
    m << "fourier-spacing          = " << fmt("%.4g", std::min(e.pme_spacing, 1.2) / 10) << "\n";
    m << "pme-order                = " << std::clamp(e.pme_order, 4, 12) << "\n";
  } else {
    m << "; no periodic cell in CAPS: the molecule is in a box 2 r_c larger, plain cut-off Coulomb (CAPS: damped shifted force)\n";
    m << "coulombtype              = Reaction-Field\n";
    m << "epsilon-rf               = 1\n";
  }
  m << "constraints              = none        ; CAPS keeps every bond flexible\n";
  return m.str();
}

std::vector<std::string> write_gromacs(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& stem) {
  const std::vector<std::string> notes = gromacs_notes(s, ff, e);
  const size_t n = s.atoms.size();

  const int nt = int(ff.type_names.size());
  std::ofstream top(stem + ".top");
  if (!top) throw std::runtime_error("cannot write " + stem + ".top");
  char b[512];
  const std::string title = export_title(s.title, ff.name);
  top << "; CAPS 0.1 · " << title << " · " << ff.name << "\n";
  top << "; kJ/mol and nm; every Lennard-Jones pair written with " << ff.mixing << " mixing applied, 1-4 pairs with their scaled σ, ε\n\n";
  // GROMOS writes its van der Waals as C6 = 4εσ⁶ and C12 = 4εσ¹² (comb-rule 1); the others as σ, ε (comb-rule 2).
  // Either way every pair is listed, so the numbers GROMACS uses are the same.
  const bool c6c12 = ff.native_gromacs_lj == "c6c12";
  auto lj_pair = [&](double sigma, double eps) {   // the two numbers for σ (Å), ε (kcal/mol)
    char t[80];
    if (c6c12) {
      const double s6 = std::pow(sigma / 10, 6);
      std::snprintf(t, sizeof t, "%.10g %.10g", 4 * eps * KJ * s6, 4 * eps * KJ * s6 * s6);
    } else {
      std::snprintf(t, sizeof t, "%.10g %.10g", sigma / 10, eps * KJ);
    }
    return std::string(t);
  };
  top << "[ defaults ]\n; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ\n";
  std::snprintf(b, sizeof b, "  1        %d          no         1.0      %.10g\n\n", c6c12 ? 1 : 2, ff.coul14);
  top << b;

  // atom types: the mass of the first atom of each type; charges per atom
  std::vector<int> first(size_t(nt), -1);
  for (size_t i = 0; i < n; ++i)
    if (first[size_t(ff.type_index[i])] < 0) first[size_t(ff.type_index[i])] = int(i);
  std::vector<std::string> tname(static_cast<size_t>(nt));
  std::set<std::string> used;
  for (int t = 0; t < nt; ++t) {
    std::string nm = type_label(ff.type_names[size_t(t)]);
    if (used.count(nm)) nm = clean(ff.type_names[size_t(t)]);
    for (int k = 2; used.count(nm); ++k) nm = clean(ff.type_names[size_t(t)]).substr(0, 12) + "_" + std::to_string(k);
    used.insert(nm);
    tname[size_t(t)] = nm;
  }
  top << (c6c12 ? "[ atomtypes ]\n; name  at.num  mass  charge  ptype  C6 (kJ/mol nm⁶)  C12 (kJ/mol nm¹²)\n"
                 : "[ atomtypes ]\n; name  at.num  mass  charge  ptype  sigma (nm)  epsilon (kJ/mol)\n");
  for (int t = 0; t < nt; ++t) {
    const int a = first[size_t(t)];
    const int z = a >= 0 ? s.atoms[size_t(a)].element : 0;
    const double mass = a >= 0 ? ff.mass[size_t(a)] : 0;
    std::snprintf(b, sizeof b, "%-16s %3d %12.6f  0.0  A %s\n", tname[size_t(t)].c_str(), z, mass, lj_pair(ff.lj[size_t(t)].sigma, ff.lj[size_t(t)].eps).c_str());
    top << b;
  }
  top << (c6c12 ? "\n[ nonbond_params ]\n; i  j  func  C6  C12\n" : "\n[ nonbond_params ]\n; i  j  func  sigma (nm)  epsilon (kJ/mol)\n");
  for (int a = 0; a < nt; ++a)
    for (int c = a; c < nt; ++c) {
      PairType p = mixed_pair(ff, a, c);
      if (auto it = ff.pair_func.find({a, c}); it != ff.pair_func.end() && it->second.form == kPairGromacs) p = {it->second.a, it->second.b};   // lj/gromacs: ε, σ
      std::snprintf(b, sizeof b, "%-16s %-16s 1 %s\n", tname[size_t(a)].c_str(), tname[size_t(c)].c_str(), lj_pair(p.sigma, p.eps).c_str());
      top << b;
    }

  // The molecule sections, as lines tagged with their atoms: each written with molecule-local numbering, so identical
  // molecules (a melt's chains, a solvent) become one [ moleculetype ] in STEM.itp, listed with their count in [ molecules ].
  std::vector<Vec3> pos;
  const bool periodic_mol = periodic_bonds(s, pos);
  // virtual sites where their atoms put them (GROMACS's rerun takes the coordinates as written; CAPS places the sites
  // at every evaluation)
  for (const auto& v : ff.vsites) {
    Vec3 c{0, 0, 0};
    double wt = 0;
    const Vec3 ref = pos[v.from[0]];
    for (size_t k = 0; k < v.from.size(); ++k) {
      const Vec3 d = s.cell.valid() ? s.cell.minimum_image(pos[v.from[k]] - ref) : pos[v.from[k]] - ref;
      c = c + (ref + d) * v.w[k];
      wt += v.w[k];
    }
    if (wt > 0) pos[v.site] = c * (1 / wt);
  }
  int nmol = 0;
  const auto mol = s.molecules(&nmol);
  enum Sec { ATOMS, BONDS, PAIRS, ANGLES, DIHEDRALS, VSITES, VSITES2, VSITES3, EXCLUSIONS, NSEC };
  static const char* sec_head[NSEC] = {"[ atoms ]\n; nr  type  resnr  residue  atom  cgnr  charge  mass\n", "[ bonds ]\n; i  j  func  parameters\n",
                                       "[ pairs ]\n; i  j  func  sigma (nm)  epsilon (kJ/mol), scaled\n", "[ angles ]\n; i  j  k  func  parameters\n",
                                       "[ dihedrals ]\n; i  j  k  l  func  parameters\n",
                                       "[ virtual_sitesn ]\n; site  func  atom weight ... (3: weighted centre)\n",
                                       "[ virtual_sites2 ]\n; site  i  j  func  a   (x = (1 - a) x_i + a x_j)\n",
                                       "[ virtual_sites3 ]\n; site  i  j  k  func  a  b   (x = (1 - a - b) x_i + a x_j + b x_k)\n", "[ exclusions ]\n"};
  struct Line { Sec sec; std::vector<uint32_t> atoms; std::string tail; };
  std::vector<Line> lines;
  auto add = [&](Sec sec, std::vector<uint32_t> at, const std::string& tail) { lines.push_back({sec, std::move(at), tail}); };

  // bonds (structure bonds without a term: function 5, a connection, so tools see the molecule)
  std::set<std::pair<uint32_t, uint32_t>> have;
  auto mark = [&](uint32_t i, uint32_t j) { have.insert({std::min(i, j), std::max(i, j)}); };
  for (const auto& t : ff.bonds) {
    std::snprintf(b, sizeof b, " 1 %.10g %.10g", t.r0 / 10, 2 * t.k * KJ * 100);
    add(BONDS, {t.i, t.j}, b);
    mark(t.i, t.j);
  }
  for (const auto& t : ff.bonds_x) {
    if (t.form == 1) std::snprintf(b, sizeof b, " 3 %.10g %.10g %.10g", t.c / 10, t.a * KJ, t.b * 10);   // Morse: b0, D, β
    else if (t.form == 2) std::snprintf(b, sizeof b, " 2 %.10g %.10g", t.b / 10, t.a * KJ * 1e4);   // GROMOS quartic
    else throw FieldError("bond form " + std::to_string(t.form) + " has no GROMACS function");
    add(BONDS, {t.i, t.j}, b);
    mark(t.i, t.j);
  }
  for (const auto& t : s.bonds)
    if (t.i != t.j && !have.count({std::min(t.i, t.j), std::max(t.i, t.j)})) {
      add(BONDS, {t.i, t.j}, " 5");
      mark(t.i, t.j);
    }
  // 1-4 pairs, each with its scaled coefficients
  {
    ForceField f14;
    f14.mixing = ff.mixing;
    f14.lj = ff.lj14_types;
    for (const auto& p : ff.pairs14) {
      const int ta = ff.type_index[p[0]], tb = ff.type_index[p[1]];
      const PairType q = ff.lj14_types.empty() ? mixed_pair(ff, ta, tb) : mixed_pair(f14, ta, tb);
      add(PAIRS, {p[0], p[1]}, " 1 " + lj_pair(q.sigma, ff.lj14 * q.eps));
    }
  }
  // explicit LJ pairs (Martini 3 polymers): their own σ, ε (GROMACS adds fudgeQQ × the pair's Coulomb, as CAPS does)
  for (const auto& p : ff.lj_pairs) {
    add(PAIRS, {p.i, p.j}, " 1 " + lj_pair(p.sigma, p.eps));
  }
  // angles (Urey–Bradley terms join their angle: GROMACS function 5)
  std::map<std::pair<uint32_t, uint32_t>, const UreyBradley*> ub;
  for (const auto& u : ff.urey_bradley) ub[{std::min(u.i, u.k), std::max(u.i, u.k)}] = &u;
  size_t ub_used = 0;
  for (const auto& a : ff.angles) {
    auto it = ub.find({std::min(a.i, a.k), std::max(a.i, a.k)});
    if (it != ub.end()) {
      ++ub_used;
      std::snprintf(b, sizeof b, " 5 %.10g %.10g %.10g %.10g", a.theta0 * R2D, 2 * a.kt * KJ, it->second->r0 / 10, 2 * it->second->kub * KJ * 100);
    } else {
      std::snprintf(b, sizeof b, " 1 %.10g %.10g", a.theta0 * R2D, 2 * a.kt * KJ);
    }
    add(ANGLES, {a.i, a.j, a.k}, b);
  }
  if (ub_used != ub.size()) throw FieldError("a Urey–Bradley term has no angle to join (GROMACS angle function 5 needs one)");
  for (const auto& a : ff.angles_x) {
    if (a.form == 1) std::snprintf(b, sizeof b, " 2 %.10g %.10g", a.b * R2D, 2 * a.a * KJ);   // K (cos θ − cos θ0)² = ½ kθ (…)², GROMOS-96 angle
    else if (a.form == 5) std::snprintf(b, sizeof b, " 10 %.10g %.10g", a.b * R2D, 2 * a.a * KJ);   // restricted bending
    else throw FieldError("angle form " + std::to_string(a.form) + " has no GROMACS function");
    add(ANGLES, {a.i, a.j, a.k}, b);
  }
  // virtual sites: function 3, each atom with its weight (the centre of mass, or a fixed combination, exactly)
  for (const auto& v : ff.vsites) {
    std::vector<uint32_t> at{v.site};
    at.insert(at.end(), v.from.begin(), v.from.end());
    double sum = 0, lo = 0;
    for (double x : v.w) sum += x, lo = std::min(lo, x);
    // a linear combination of two or three atoms (weights summing to 1, perhaps negative: a site outside its atoms)
    // as GROMACS's own linear sites; virtual_sitesn takes positive weights only
    // (two atoms always so: a site built on it may be a virtual_sites3, which GROMACS allows only on lower functions)
    if ((v.from.size() == 2 || (v.from.size() == 3 && lo < 0)) && std::fabs(sum - 1) < 1e-9) {
      add(v.from.size() == 2 ? VSITES2 : VSITES3, at,
          v.from.size() == 2 ? " 1 " + fmt("%.12g", v.w[1]) : " 1 " + fmt("%.12g", v.w[1]) + " " + fmt("%.12g", v.w[2]));
      continue;
    }
    if (lo < 0) throw FieldError("a virtual site with negative weights on " + std::to_string(v.from.size()) + " atoms has no GROMACS form");
    std::string w;
    for (double x : v.w) w += " " + fmt("%.12g", x);
    add(VSITES, at, w);
  }
  // proper torsions: every Fourier term of a quadruple on consecutive lines (function 9)
  std::map<std::array<uint32_t, 4>, std::vector<const TorsionTerm*>> quad;
  std::vector<std::array<uint32_t, 4>> order;
  for (const auto& t : ff.dihedrals) {
    const std::array<uint32_t, 4> k{t.i, t.j, t.k, t.l};
    auto& v = quad[k];
    if (v.empty()) order.push_back(k);
    v.push_back(&t);
  }
  for (const auto& k : order)
    for (const auto* t : quad[k]) {
      if (t->n < 0) throw FieldError("a torsion with multiplicity < 0 has no GROMACS form");
      // GROMOS's own proper dihedral is function 1, one term per quadruple; several terms need function 9
      std::snprintf(b, sizeof b, " %d %.10g %.10g %d", c6c12 && quad[k].size() == 1 ? 1 : 9, t->delta * R2D, t->v * KJ, t->n);
      add(DIHEDRALS, {k[0], k[1], k[2], k[3]}, b);
    }
  // impropers: periodic (function 4, AMBER order with the centre third) and harmonic (function 2)
  std::map<std::array<uint32_t, 4>, int> imp_count;
  for (const auto& t : ff.impropers) imp_count[{t.i, t.j, t.k, t.l}]++;
  for (const auto& t : ff.impropers) {
    if (t.n < 0) throw FieldError("an improper with multiplicity < 0 has no GROMACS form");
    const int fn = imp_count[{t.i, t.j, t.k, t.l}] > 1 ? 9 : 4;   // several terms on one quadruple: function 9
    std::snprintf(b, sizeof b, " %d %.10g %.10g %d", fn, t.delta * R2D, t.v * KJ, t.n);
    add(DIHEDRALS, {t.i, t.j, t.k, t.l}, b);
  }
  for (const auto& t : ff.impropers_harmonic) {
    std::snprintf(b, sizeof b, " 2 %.10g %.10g", t.chi0 * R2D, 2 * t.k2 * KJ);
    add(DIHEDRALS, {t.i, t.j, t.k, t.l}, b);
  }
  for (const auto& t : ff.cbt) {   // combined bending–torsion: k 1, the coefficients in kJ/mol
    std::snprintf(b, sizeof b, " 11 1 %.10g %.10g %.10g %.10g %.10g", t.a[0] * KJ, t.a[1] * KJ, t.a[2] * KJ, t.a[3] * KJ, t.a[4] * KJ);
    add(DIHEDRALS, {t.i, t.j, t.k, t.l}, b);
  }
  // exclusions: CAPS's own list (nrexcl 0)
  for (uint32_t i = 0; i < ff.excluded.size(); ++i) {
    std::vector<uint32_t> at{i};
    for (uint32_t j : ff.excluded[i])
      if (j > i) at.push_back(j);
    if (at.size() > 1) add(EXCLUSIONS, at, "");
  }

  // molecules: contiguous atom ranges in order, else one molecule type for the whole system
  std::vector<uint32_t> first_atom(static_cast<size_t>(nmol), UINT32_MAX), last_atom(static_cast<size_t>(nmol), 0), count(static_cast<size_t>(nmol), 0);
  for (uint32_t i = 0; i < n; ++i) {
    const auto m = size_t(mol[i]);
    first_atom[m] = std::min(first_atom[m], i);
    last_atom[m] = std::max(last_atom[m], i);
    ++count[m];
  }
  bool contiguous = true;
  for (int m = 0; m < nmol && contiguous; ++m) contiguous = last_atom[size_t(m)] - first_atom[size_t(m)] + 1 == count[size_t(m)];
  for (const auto& l : lines)
    for (uint32_t a : l.atoms) contiguous = contiguous && mol[a] == mol[l.atoms[0]];
  std::vector<int> mol_of(n);
  std::vector<uint32_t> start;
  if (contiguous) {
    std::vector<int> byfirst(static_cast<size_t>(nmol));
    for (int m = 0; m < nmol; ++m) byfirst[size_t(m)] = m;
    std::sort(byfirst.begin(), byfirst.end(), [&](int a, int c) { return first_atom[size_t(a)] < first_atom[size_t(c)]; });
    std::vector<int> rank(static_cast<size_t>(nmol));
    for (int r = 0; r < nmol; ++r) rank[size_t(byfirst[size_t(r)])] = r, start.push_back(first_atom[size_t(byfirst[size_t(r)])]);
    for (uint32_t i = 0; i < n; ++i) mol_of[i] = rank[size_t(mol[i])];
  } else {
    nmol = 1;
    start = {0};
    std::fill(mol_of.begin(), mol_of.end(), 0);
  }
  std::vector<std::vector<const Line*>> by_mol(static_cast<size_t>(nmol));
  for (const auto& l : lines) by_mol[size_t(mol_of[l.atoms[0]])].push_back(&l);
  // atom names numbered within each molecule type, the same in the itp and the gro
  std::vector<std::string> aname(n);
  {
    std::map<int, int> per_element;
    for (uint32_t i = 0; i < n; ++i) {
      if (i == 0 || mol_of[i] != mol_of[i - 1]) per_element.clear();
      aname[i] = atom_label(s.atoms[i], ff.type_names[size_t(ff.type_index[i])], per_element);
    }
  }
  // each molecule's sections with local numbering; identical text is one molecule type
  auto body_of = [&](int m) {
    const uint32_t s0 = start[size_t(m)];
    const uint32_t s1 = m + 1 < nmol ? start[size_t(m) + 1] : uint32_t(n);
    std::ostringstream o;
    long res0 = -1;
    o << sec_head[ATOMS];
    for (uint32_t i = s0; i < s1; ++i) {
      const Atom& a = s.atoms[i];
      long res = a.resid > 0 ? long(a.resid) : 1;
      if (res0 < 0) res0 = res;
      const std::string rn = clean(a.resname.empty() ? "MOL" : a.resname).substr(0, 5);
      const std::string& an = aname[i];
      std::snprintf(b, sizeof b, "%7u %-16s %6ld %-5s %-5s %7u %14.10f %12.6f\n", i - s0 + 1, tname[size_t(ff.type_index[i])].c_str(), res - res0 + 1, rn.c_str(),
                    an.c_str(), i - s0 + 1, ff.charge[i], ff.mass[i]);
      o << b;
    }
    for (int sec = BONDS; sec < NSEC; ++sec) {
      bool head = false;
      for (const Line* l : by_mol[size_t(m)]) {
        if (l->sec != sec) continue;
        if (!head) o << "\n" << sec_head[sec], head = true;
        if (sec == VSITES) {   // virtual_sitesn: the site, function 3, then each atom and its weight
          std::istringstream ws(l->tail);
          o << std::setw(7) << (l->atoms[0] - s0 + 1) << "  3";
          for (size_t k = 1; k < l->atoms.size(); ++k) {
            double wk = 0;
            ws >> wk;
            o << std::setw(7) << (l->atoms[k] - s0 + 1) << " " << fmt("%.12g", wk);
          }
          o << "\n";
          continue;
        }
        for (uint32_t a : l->atoms) o << std::setw(7) << (a - s0 + 1);
        o << l->tail << "\n";
      }
    }
    return o.str();
  };
  std::vector<std::string> kinds;              // molecule-type bodies
  std::vector<std::pair<int, int>> runs;       // (kind, count) in order
  for (int m = 0; m < nmol; ++m) {
    const std::string body = body_of(m);
    int k = int(std::find(kinds.begin(), kinds.end(), body) - kinds.begin());
    if (k == int(kinds.size())) kinds.push_back(body);
    if (!runs.empty() && runs.back().first == k) ++runs.back().second;
    else runs.push_back({k, 1});
  }
  // STEM.itp: the molecule types
  const std::string itp_name = std::filesystem::path(stem + ".itp").filename().string();
  std::ofstream itp(stem + ".itp");
  if (!itp) throw std::runtime_error("cannot write " + stem + ".itp");
  itp << "; CAPS 0.1 · molecule types of " << title << " · " << ff.name << " · kJ/mol and nm\n";
  std::vector<std::string> kname(kinds.size());
  for (size_t k = 0; k < kinds.size(); ++k) {
    kname[k] = contiguous ? (kinds.size() == 1 ? std::string("MOL") : "MOL" + std::to_string(k + 1)) : std::string("SYSTEM");
    int copies = 0;
    for (const auto& r : runs) if (size_t(r.first) == k) copies += r.second;
    itp << "\n[ moleculetype ]\n; name  nrexcl   (" << copies << (copies == 1 ? " molecule" : " molecules") << "; exclusions listed below)\n"
        << kname[k] << "  0\n\n" << kinds[k];
  }
  itp.close();
  top << "\n#include \"" << itp_name << "\"\n";
  top << "\n[ system ]\n" << title << "\n\n[ molecules ]\n; name  count\n";
  for (const auto& r : runs) top << kname[size_t(r.first)] << "  " << r.second << "\n";
  top.close();

  // coordinates: nm, 8 decimals (GROMACS reads the precision from the first line)
  std::ofstream gro(stem + ".gro");
  if (!gro) throw std::runtime_error("cannot write " + stem + ".gro");
  gro << title << "\n" << n << "\n";
  Vec3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
  for (const auto& p : pos)
    for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], p[k]); hi[k] = std::max(hi[k], p[k]); }
  Vec3 shift{0, 0, 0};
  double box[9] = {0};
  if (s.cell.valid()) {
    const Cell& c = s.cell;
    const double v[9] = {c.a[0], c.b[1], c.c[2], c.a[1], c.a[2], c.b[0], c.b[2], c.c[0], c.c[1]};
    for (int k = 0; k < 9; ++k) box[k] = v[k] / 10;
  } else {
    // a box 2 r_c larger than the molecule, the molecule centred
    for (int k = 0; k < 3; ++k) {
      box[k] = (hi[k] - lo[k] + 2 * e.cutoff + 2) / 10;
      shift[k] = box[k] * 5 - 0.5 * (lo[k] + hi[k]);
    }
  }
  for (size_t i = 0; i < n; ++i) {
    const Atom& a = s.atoms[i];
    const long res = (a.resid > 0 ? long(a.resid) : s.has_mol && a.mol > 0 ? long(a.mol) : long(mol[i] + 1)) % 100000;
    const std::string rn = clean(a.resname.empty() ? "MOL" : a.resname).substr(0, 5);
    const std::string& an = aname[i];
    const Vec3 p = pos[i] + shift;
    std::snprintf(b, sizeof b, "%5ld%-5s%5s%5ld%13.8f%13.8f%13.8f\n", res, rn.c_str(), an.c_str(), long((i + 1) % 100000), p[0] / 10, p[1] / 10, p[2] / 10);
    gro << b;
  }
  std::snprintf(b, sizeof b, "%12.8f %12.8f %12.8f", box[0], box[1], box[2]);
  gro << b;
  if (std::fabs(box[3]) + std::fabs(box[4]) + std::fabs(box[5]) + std::fabs(box[6]) + std::fabs(box[7]) + std::fabs(box[8]) > 1e-12) {
    std::snprintf(b, sizeof b, " %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f", box[3], box[4], box[5], box[6], box[7], box[8]);
    gro << b;
  }
  gro << "\n";
  gro.close();

  // run parameters: a single-point evaluation; the Studio replaces the run section with the ensemble
  std::ofstream mdp(stem + ".mdp");
  if (!mdp) throw std::runtime_error("cannot write " + stem + ".mdp");
  mdp << "; GROMACS run parameters written by CAPS: " << ff.name << "\n";
  mdp << "integrator               = md\nnsteps                   = 0\ndt                       = " << (ff.native_timestep > 0 ? ff.native_timestep : 0.5) / 1000 << "\n";
  mdp << "nstcalcenergy            = 1\nnstenergy                = 1\n";
  mdp << gromacs_mdp(s, ff, e);
  if (periodic_mol) mdp << "periodic-molecules       = yes         ; bonds cross the cell: an infinite network\n";
  mdp.close();

  return notes;
}

}  // namespace caps
