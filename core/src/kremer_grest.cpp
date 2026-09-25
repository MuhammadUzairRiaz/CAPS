#include "caps/kremer_grest.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <stdexcept>

namespace caps {

System kremer_grest(const KgOptions& o, KgReport* rep) {
  if (o.chains < 1 || o.beads < 2) throw std::invalid_argument("at least one chain of two beads");
  if (o.density <= 0.05 || o.density > 1.2) throw std::invalid_argument("reduced density ρσ³ between 0.05 and 1.2");
  const double L = std::cbrt(double(o.chains) * o.beads / o.density);
  System s;
  s.title = "Kremer–Grest melt";
  s.cell.a = {L, 0, 0};
  s.cell.b = {0, L, 0};
  s.cell.c = {0, 0, L};
  s.types.push_back({1, 1.0, "B1"});
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0, 1);
  std::normal_distribution<double> G(0, 1);
  auto unit = [&]() {
    Vec3 v{G(rng), G(rng), G(rng)};
    const double n = std::sqrt(dot(v, v));
    return v * (1 / n);
  };
  double r2sum = 0;
  for (int c = 0; c < o.chains; ++c) {
    Vec3 p{U(rng) * L, U(rng) * L, U(rng) * L};
    Vec3 prev{0, 0, 0};
    const uint32_t first = uint32_t(s.atoms.size());
    for (int k = 0; k < o.beads; ++k) {
      if (k > 0) {
        // the next direction: bond angle above 60° (no back-folding); with k_θ, Boltzmann weight exp(−k_θ (1 − cos θ)) on the
        // angle between consecutive bonds (θ = 0: straight), by rejection
        Vec3 d;
        for (int tries = 0;; ++tries) {
          d = unit();
          if (k == 1) break;
          const double cosb = dot(d, prev);       // 1: straight on
          if (cosb < -0.5) continue;               // folds back past 120° between bonds
          if (o.k_theta <= 0 || U(rng) < std::exp(-o.k_theta * (1 - cosb)) || tries > 1000) break;
        }
        p = p + d * o.bond;
        prev = d;
      }
      Atom a;
      a.element = 6;
      a.type = 1;
      a.name = "B1";
      a.mol = c + 1;
      a.pos = p;   // unwrapped along the chain
      s.atoms.push_back(a);
      if (k > 0) s.bonds.push_back({uint32_t(s.atoms.size() - 2), uint32_t(s.atoms.size() - 1), 1});
    }
    const Vec3 ee = s.atoms.back().pos - s.atoms[first].pos;
    r2sum += dot(ee, ee) / (o.beads - 1);
  }
  s.has_mol = true;
  s.bonds_from_file = true;
  s.unwrapped = true;
  if (rep) {
    rep->box = L;
    rep->mean_r2 = r2sum / o.chains;
    // closest non-bonded pair on a sample of beads (cell list would be exact; a sample is enough for the report)
    double best = 1e9;
    const size_t n = s.atoms.size(), step = std::max<size_t>(1, n / 400);
    for (size_t i = 0; i < n; i += step)
      for (size_t j = 0; j < n; ++j) {
        if (j == i || (s.atoms[i].mol == s.atoms[j].mol && (j + 1 == i || i + 1 == j))) continue;
        const Vec3 d = s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos);
        best = std::min(best, std::sqrt(dot(d, d)));
      }
    rep->closest = best;
  }
  return s;
}

void write_kg_lammps(const System& s, const KgOptions& o, const std::string& stem, double pushoff_steps, double run_steps) {
  const bool angles = o.k_theta > 0;
  const double L = s.cell.a[0];
  {
    std::ofstream f(stem + ".data");
    if (!f) throw std::runtime_error("cannot write " + stem + ".data");
    size_t nang = 0;
    for (size_t i = 1; i + 1 < s.atoms.size(); ++i) nang += angles && s.atoms[i - 1].mol == s.atoms[i].mol && s.atoms[i + 1].mol == s.atoms[i].mol;
    char b[256];
    f << "Kremer-Grest melt: " << o.chains << " chains x " << o.beads << " beads, rho sigma^3 = " << o.density << " (CAPS, reduced units)\n\n";
    f << s.atoms.size() << " atoms\n" << s.bonds.size() << " bonds\n";
    if (angles) f << nang << " angles\n";
    f << "1 atom types\n1 bond types\n";
    if (angles) f << "1 angle types\n";
    std::snprintf(b, sizeof b, "\n0 %.6f xlo xhi\n0 %.6f ylo yhi\n0 %.6f zlo zhi\n\nMasses\n\n1 1.0\n\nAtoms # angle\n\n", L, L, L);
    f << b;
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      const Vec3 w = s.cell.wrap(s.atoms[i].pos);
      const Vec3 f3 = s.cell.to_fractional(s.atoms[i].pos);
      const int ix = int(std::floor(f3[0])), iy = int(std::floor(f3[1])), iz = int(std::floor(f3[2]));
      std::snprintf(b, sizeof b, "%zu %lld 1 %.6f %.6f %.6f %d %d %d\n", i + 1, (long long)s.atoms[i].mol, w[0], w[1], w[2], ix, iy, iz);
      f << b;
    }
    f << "\nBonds\n\n";
    for (size_t k = 0; k < s.bonds.size(); ++k) f << k + 1 << " 1 " << s.bonds[k].i + 1 << " " << s.bonds[k].j + 1 << "\n";
    if (angles) {
      f << "\nAngles\n\n";
      size_t k = 0;
      for (size_t i = 1; i + 1 < s.atoms.size(); ++i)
        if (s.atoms[i - 1].mol == s.atoms[i].mol && s.atoms[i + 1].mol == s.atoms[i].mol) f << ++k << " 1 " << i << " " << i + 1 << " " << i + 2 << "\n";
    }
  }
  std::ofstream f(stem + ".in");
  if (!f) throw std::runtime_error("cannot write " + stem + ".in");
  const std::string data = stem.substr(stem.find_last_of("/\\") == std::string::npos ? 0 : stem.find_last_of("/\\") + 1) + ".data";
  const std::string angle = angles ? "angle_style cosine\nangle_coeff 1 " + std::to_string(o.k_theta) + "\n" : std::string();
  const std::string base = data.substr(0, data.size() - 5);
  char b[4096];
  std::snprintf(b, sizeof b,
                "# Kremer-Grest melt written by CAPS: %d x %d beads, rho sigma^3 = %.3f, k_theta = %.3f\n"
                "# Kremer & Grest, J. Chem. Phys. 92, 5057 (1990); push-off after Auhl et al., J. Chem. Phys. 119, 12718 (2003)\n"
                "units lj\natom_style angle\nboundary p p p\nread_data %s\n"
                "bond_style fene\nbond_coeff 1 30.0 1.5 1.0 1.0\nspecial_bonds fene\n%s"
                "comm_modify cutoff 2.0\nneighbor 0.4 bin\nneigh_modify every 1 delay 0 check yes\ntimestep 0.01\n"
                "velocity all create 1.0 %llu dist gaussian\nfix lang all langevin 1.0 1.0 1.0 %llu\nthermo_style custom step temp pe ebond press\nthermo 1000\n\n"
                "# 1. push-off: a soft repulsion ramped up, displacements limited, so random-walk overlaps open without breaking bonds\n"
                "pair_style soft 1.122462\npair_coeff * * 0.0\nvariable prefactor equal ramp(0,60)\nfix push all adapt 1 pair soft a * * v_prefactor\n"
                "fix move all nve/limit 0.02\nrun %.0f\nunfix push\n\n"
                "# 2. Kremer-Grest: WCA repulsion (LJ cut at 2^(1/6) sigma, shifted); still limited for a while, then free\n"
                "pair_style lj/cut 1.122462\npair_modify shift yes\npair_coeff * * 1.0 1.0 1.122462\n"
                "run 2000\nunfix move\nfix move all nve\n"
                "dump d all custom 10000 %s.lammpstrj id mol type xu yu zu\nrun %.0f\nwrite_data %s_equilibrated.data\n",
                o.chains, o.beads, o.density, o.k_theta, data.c_str(), angle.c_str(), (unsigned long long)(o.seed + 1), (unsigned long long)(o.seed + 2),
                pushoff_steps, base.c_str(), run_steps, base.c_str());
  f << b;
}

}  // namespace caps
