// Random-coefficient checks of CAPS's extra functional forms against LAMMPS: writes one LAMMPS data + input per
// case and the CAPS energy / forces for it.
#include <cmath>
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include "caps/field.hpp"
#include "caps/grow.hpp"
using namespace caps;
constexpr double D = 3.14159265358979323846 / 180;

int main(int argc, char** argv) {
  const std::string which = argv[1];
  GrowOptions g; g.chains = 3; g.dp = 4; g.density = 0.30; g.seed = 5;
  System s = grow(g);
  ForceField base = assign_gaff(s);
  std::mt19937_64 rng(argc > 2 ? atoi(argv[2]) : 1);
  std::uniform_real_distribution<double> u(0, 1);
  auto r = [&](double lo, double hi) { return lo + (hi - lo) * u(rng); };
  std::normal_distribution<double> nd(0, 0.08);
  std::vector<double> x;
  for (auto& a : s.atoms) { for (int k = 0; k < 3; ++k) a.pos[k] += nd(rng); x.insert(x.end(), a.pos.begin(), a.pos.end()); }
  ForceField ff = base;
  ff.bonds.clear(); ff.angles.clear(); ff.dihedrals.clear(); ff.impropers.clear();
  for (auto& q : ff.charge) q = 0;
  const size_t nt = ff.lj.size();
  for (auto& p : ff.lj) p.eps = 0;   // pair energy off unless the case sets it
  FILE* d = fopen("case.data", "w");
  FILE* in = fopen("in.case", "w");
  const auto nb = s.neighbours();
  std::string styles, coeffs, terms;
  char b[512];
  int nterm = 0;
  std::string section;
  fprintf(in, "units real\natom_style full\n");
  if (which == "morse" || which == "gromos") {
    for (const auto& bd : base.bonds) {
      const double r0 = bd.r0 + r(-0.1, 0.1);
      if (which == "morse") { ff.bonds_x.push_back({bd.i, bd.j, 1, r(50, 120), r(1.5, 2.5), r0}); snprintf(b, sizeof b, "%d %.10g %.10g %.10g\n", nterm + 1, ff.bonds_x.back().a, ff.bonds_x.back().b, r0); }
      else { ff.bonds_x.push_back({bd.i, bd.j, 2, r(50, 200), r0, 0}); snprintf(b, sizeof b, "%d %.10g %.10g\n", nterm + 1, ff.bonds_x.back().a, r0); }
      coeffs += b; snprintf(b, sizeof b, "%d %d %u %u\n", nterm + 1, nterm + 1, bd.i + 1, bd.j + 1); terms += b; ++nterm;
    }
    fprintf(in, "bond_style %s\n", which.c_str()); section = "Bond";
  } else if (which == "cossq" || which == "charmm") {
    for (const auto& an : base.angles) {
      const double t0 = an.theta0 / D + r(-5, 5), K = r(30, 80);
      if (which == "cossq") { ff.angles_x.push_back({an.i, an.j, an.k, 1, K, t0 * D}); snprintf(b, sizeof b, "%d %.10g %.10g\n", nterm + 1, K, t0); }
      else {
        const double kub = r(5, 40), rub = r(2.0, 2.6);
        ff.angles.push_back({an.i, an.j, an.k, K, t0 * D}); ff.urey_bradley.push_back({an.i, an.k, kub, rub});
        snprintf(b, sizeof b, "%d %.10g %.10g %.10g %.10g\n", nterm + 1, K, t0, kub, rub);
      }
      coeffs += b; snprintf(b, sizeof b, "%d %d %u %u %u\n", nterm + 1, nterm + 1, an.i + 1, an.j + 1, an.k + 1); terms += b; ++nterm;
    }
    fprintf(in, "angle_style %s\n", which == "cossq" ? "cosine/squared" : "charmm"); section = "Angle";
  } else if (which == "umbrella") {
    for (uint32_t c = 0; c < s.atoms.size(); ++c) {
      if (nb[c].size() != 3) continue;
      const double K = r(10, 60);
      ff.inversions.push_back({c, nb[c][0], nb[c][1], nb[c][2], K, 0, 1});
      // LAMMPS umbrella: I is the centre, the angle between IL and the IJK plane; three permutations at K/3
      const uint32_t o[3] = {nb[c][0], nb[c][1], nb[c][2]};
      for (int p = 0; p < 3; ++p) {
        snprintf(b, sizeof b, "%d %.10g 0.0\n", nterm + 1, K / 3); coeffs += b;
        snprintf(b, sizeof b, "%d %d %u %u %u %u\n", nterm + 1, nterm + 1, c + 1, o[(p + 1) % 3] + 1, o[(p + 2) % 3] + 1, o[p] + 1); terms += b; ++nterm;
      }
    }
    fprintf(in, "improper_style umbrella\n"); section = "Improper";
  } else if (which == "buck" || which == "pmorse" || which == "lj14") {
    // pair forms: every type pair
    for (size_t a = 0; a < nt; ++a)
      for (size_t c = a; c < nt; ++c) {
        if (which == "buck") ff.pair_func[{int(a), int(c)}] = {1, r(2000, 20000), r(0.25, 0.32), r(5, 50)};
        else if (which == "pmorse") ff.pair_func[{int(a), int(c)}] = {2, r(0.05, 0.2), r(1.2, 2.0), r(3.0, 4.0)};
      }
    if (which == "lj14") {
      ff.lj14_types.clear();
      for (size_t a = 0; a < nt; ++a) ff.lj14_types.push_back({r(0.02, 0.1), r(2.5, 3.5)});
      ff.lj14 = 1.0;
      ff.pairs14 = base.pairs14;
    } else {
      ff.lj14 = 1.0;   // special_bonds lj 0 0 1 in LAMMPS
    }
    fprintf(in, "%s\n", which == "buck" ? "pair_style buck 10.0\npair_modify tail yes" : which == "pmorse" ? "pair_style morse 10.0\npair_modify shift yes" : "pair_style lj/charmm/coul/charmm 9.99 10.0");
  }
  // CAPS
  EnergyOptions eo; eo.coulomb = false; eo.tail = which == "buck" || which == "lj14"; eo.cutoff = 10.0;
  if (which == "pmorse") eo.tail = false;
  Evaluator ev(ff, eo);
  std::vector<double> f;
  EnergyTerms e = ev.compute(x, s.cell, f);
  // LAMMPS data
  fprintf(d, "case %s\n\n%zu atoms\n%zu atom types\n", which.c_str(), s.atoms.size(), nt);
  const bool dih = which == "lj14";
  const bool pairs = which == "buck" || which == "pmorse";   // bonds (zero energy) so LAMMPS excludes 1-2 / 1-3 pairs
  if (pairs) fprintf(d, "%zu bonds\n1 bond types\n", base.bonds.size());
  if (section == "Bond") fprintf(d, "%d bonds\n%d bond types\n", nterm, nterm);
  if (section == "Angle") fprintf(d, "%d angles\n%d angle types\n", nterm, nterm);
  if (section == "Improper") fprintf(d, "%d impropers\n%d improper types\n", nterm, nterm);
  if (dih) fprintf(d, "%zu dihedrals\n1 dihedral types\n", base.pairs14.size());
  const Cell& c = s.cell;
  fprintf(d, "%.8f %.8f xlo xhi\n%.8f %.8f ylo yhi\n%.8f %.8f zlo zhi\n\nMasses\n\n", c.origin[0], c.origin[0] + c.a[0], c.origin[1], c.origin[1] + c.b[1], c.origin[2], c.origin[2] + c.c[2]);
  for (size_t t = 0; t < nt; ++t) fprintf(d, "%zu 12.0\n", t + 1);
  if (!section.empty()) fprintf(d, "\n%s Coeffs\n\n%s", section.c_str(), coeffs.c_str());
  if (which == "buck" || which == "pmorse" || which == "lj14") {
    fprintf(d, "\nPairIJ Coeffs\n\n");
    ForceField f14; f14.mixing = ff.mixing; f14.lj = ff.lj14_types;
    for (size_t a = 0; a < nt; ++a)
      for (size_t cc = a; cc < nt; ++cc) {
        if (which == "lj14") { PairType p = mixed_pair(f14, int(a), int(cc)); fprintf(d, "%zu %zu 0.0 1.0 %.10g %.10g\n", a + 1, cc + 1, p.eps, p.sigma); }
        else { auto pf = ff.pair_func.at({int(a), int(cc)}); fprintf(d, "%zu %zu %.10g %.10g %.10g\n", a + 1, cc + 1, pf.a, pf.b, pf.c); }
      }
  }
  if (dih) fprintf(d, "\nDihedral Coeffs\n\n1 0.0 1 0 1.0\n");
  fprintf(d, "\nAtoms\n\n");
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    Vec3 fr = c.to_fractional(s.atoms[i].pos); int im[3];
    for (int k = 0; k < 3; ++k) im[k] = int(std::floor(fr[k]));
    Vec3 w = s.atoms[i].pos - (c.a * im[0] + c.b * im[1] + c.c * im[2]);
    fprintf(d, "%zu 1 %d 0.0 %.10f %.10f %.10f %d %d %d\n", i + 1, ff.type_index[i] + 1, w[0], w[1], w[2], im[0], im[1], im[2]);
  }
  if (!section.empty()) fprintf(d, "\n%ss\n\n%s", section.c_str(), terms.c_str());
  if (pairs) {
    fprintf(d, "\nBonds\n\n");
    for (size_t k = 0; k < base.bonds.size(); ++k) fprintf(d, "%zu 1 %u %u\n", k + 1, base.bonds[k].i + 1, base.bonds[k].j + 1);
  }
  if (dih) {
    // one dihedral per 1-4 pair (a path i-j-k-l), K = 0: only its 1-4 LJ counts
    fprintf(d, "\nDihedrals\n\n");
    int k = 0;
    std::set<std::pair<uint32_t, uint32_t>> want(base.pairs14.begin() == base.pairs14.end() ? std::set<std::pair<uint32_t, uint32_t>>() : std::set<std::pair<uint32_t, uint32_t>>());
    for (const auto& p : base.pairs14) want.insert({p[0], p[1]});
    std::set<std::pair<uint32_t, uint32_t>> done;
    for (uint32_t i = 0; i < s.atoms.size(); ++i)
      for (uint32_t j : nb[i]) for (uint32_t kk : nb[j]) { if (kk == i) continue; for (uint32_t l : nb[kk]) {
        if (l == j || l == i) continue;
        auto key = std::make_pair(std::min(i, l), std::max(i, l));
        if (!want.count(key) || done.count(key)) continue;
        done.insert(key);
        fprintf(d, "%d 1 %u %u %u %u\n", ++k, i + 1, j + 1, kk + 1, l + 1);
      } }
    fprintf(in, "dihedral_style charmm\nspecial_bonds lj 0 0 0 coul 0 0 0\n");
  } else {
    fprintf(in, "special_bonds lj 0 0 %g coul 0 0 0\n", which == "buck" || which == "pmorse" ? 1.0 : 0.0);
  }
  if (pairs) { fprintf(in, "bond_style zero\n"); }
  fprintf(in, "read_data case.data\ncomm_modify cutoff 12.0\n");
  if (pairs) fprintf(in, "bond_coeff *\n");
  fprintf(in, "thermo_style custom ebond eangle edihed eimp evdwl\nthermo_modify format float %%.10f\ndump dd all custom 1 lmp.dump id fx fy fz\ndump_modify dd sort id format float %%.10f\nrun 0\n");
  fclose(d); fclose(in);
  FILE* o = fopen("caps.txt", "w");
  fprintf(o, "%.10f %.10f %.10f %.10f %.10f\n", e.bond, e.angle, e.dihedral, e.improper, e.vdw);
  for (size_t i = 0; i < f.size() / 3; ++i) fprintf(o, "%zu %.10f %.10f %.10f\n", i + 1, f[3 * i], f[3 * i + 1], f[3 * i + 2]);
  fclose(o);
}
