#include "caps/rng.hpp"
#include "caps/kremer_grest.hpp"
#include "caps/exchange_moves.hpp"

#include <algorithm>
#include <array>
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
  caps::UniformReal<double> U(0, 1);
  caps::Normal<double> G(0, 1);
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

ForceField kremer_grest_forcefield(const System& s, const KgOptions& o) {
  const KgUnits u = kg_units(o);
  const bool mapped = u.eps > 0;
  const double sg = mapped ? o.sigma : 1.0, ep = mapped ? u.eps : 1.0, ms = mapped ? o.bead_mass : 1.0;
  const size_t n = s.atoms.size();
  ForceField F;
  char b[200];
  if (mapped) std::snprintf(b, sizeof b, "Kremer–Grest (FENE + WCA) · σ = %.4g Å, ε = k_B·%.4g K, m = %.4g g/mol", o.sigma, o.temperature, o.bead_mass);
  else std::snprintf(b, sizeof b, "Kremer–Grest (FENE + WCA) · reduced units");
  F.name = b;
  F.type_names = {"KG"};
  F.lj = {{ep, sg}};
  F.pair_form = "lj12-6";
  F.mixing = "geometric";
  F.cutoff = std::pow(2.0, 1.0 / 6) * sg;   // WCA: the model's own cut-off, shifted to zero there
  F.lj_shift = true;
  F.keep13 = true, F.lj14 = 1.0, F.coul14 = 1.0;   // special_bonds fene: 1-3 and 1-4 pairs in full
  F.type_index.assign(n, 0);
  F.atom_type.assign(n, "KG");
  F.why.assign(n, "Kremer–Grest bead");
  F.charge.assign(n, 0.0);
  F.mass.assign(n, ms);
  F.excluded.assign(n, {});
  for (const auto& bd : s.bonds) {
    F.bonds_x.push_back({bd.i, bd.j, 3, 30.0 * ep / (sg * sg), 1.5 * sg, ep, sg});
    F.excluded[bd.i].push_back(bd.j), F.excluded[bd.j].push_back(bd.i);
  }
  for (auto& e : F.excluded) std::sort(e.begin(), e.end());
  if (o.k_theta > 0) {
    const auto nb = s.neighbours();
    for (uint32_t j = 0; j < nb.size(); ++j)
      for (size_t a = 0; a < nb[j].size(); ++a)
        for (size_t c = a + 1; c < nb[j].size(); ++c) F.angles_x.push_back({nb[j][a], j, nb[j][c], 2, o.k_theta * ep, 0.0});   // K (1 + cos θ)
  }
  F.native_timestep = 0.01 * (mapped ? u.tau_fs : 48.88821291);   // Δt = 0.01 τ
  F.notes.push_back(mapped ? "Kremer–Grest mapped as given: σ = " + std::to_string(o.sigma) + " Å, ε = k_B T at " + std::to_string(o.temperature) + " K, m = " +
                                 std::to_string(o.bead_mass) + " g/mol; τ = " + std::to_string(u.tau_fs / 1000) + " ps"
                           : "Kremer–Grest in reduced units, read by CAPS as σ = 1 Å, ε = 1 kcal/mol, m = 1 g/mol: T* = 1 is 503.2 K and τ = 48.9 fs; the LAMMPS export is in units lj");
  return F;
}

KgUnits kg_units(const KgOptions& o) {
  KgUnits u;
  if (!(o.sigma > 0 && o.temperature > 0 && o.bead_mass > 0)) return u;
  constexpr double kB = 0.0019872067;          // kcal/(mol K), LAMMPS units real
  constexpr double tau_real = 48.88821291;     // fs: √(g/mol · Å² / (kcal/mol)), LAMMPS's time unit in units real
  u.eps = kB * o.temperature;
  u.tau_fs = o.sigma * std::sqrt(o.bead_mass / u.eps) * tau_real;
  u.density = o.density * o.bead_mass / 6.02214076e23 / std::pow(o.sigma * 1e-8, 3);
  return u;
}

void write_kg_lammps(const System& s, const KgOptions& o, const std::string& stem, double pushoff_steps, double run_steps) {
  const bool angles = o.k_theta > 0;
  const KgUnits u = kg_units(o);
  const bool real = u.eps > 0;
  // reduced units: σ = ε = m = τ = 1; real units: Å, kcal/mol, g/mol, fs
  const double sg = real ? o.sigma : 1.0, ep = real ? u.eps : 1.0, ms = real ? o.bead_mass : 1.0, tu = real ? u.tau_fs : 1.0, T = real ? o.temperature : 1.0;
  const double L = s.cell.a[0] * sg;
  {
    std::ofstream f(stem + ".data");
    if (!f) throw std::runtime_error("cannot write " + stem + ".data");
    size_t nang = 0;
    for (size_t i = 1; i + 1 < s.atoms.size(); ++i) nang += angles && s.atoms[i - 1].mol == s.atoms[i].mol && s.atoms[i + 1].mol == s.atoms[i].mol;
    char b[256];
    f << "Kremer-Grest melt: " << o.chains << " chains x " << o.beads << " beads, rho sigma^3 = " << o.density << (real ? " (CAPS, units real: sigma = " + std::to_string(o.sigma) + " A)" : std::string(" (CAPS, reduced units)")) << "\n\n";
    f << s.atoms.size() << " atoms\n" << s.bonds.size() << " bonds\n";
    if (angles) f << nang << " angles\n";
    f << "1 atom types\n1 bond types\n";
    if (angles) f << "1 angle types\n";
    std::snprintf(b, sizeof b, "\n0 %.6f xlo xhi\n0 %.6f ylo yhi\n0 %.6f zlo zhi\n\nMasses\n\n1 %.6g\n\nAtoms # angle\n\n", L, L, L, ms);
    f << b;
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      const Vec3 w = s.cell.wrap(s.atoms[i].pos) * sg;
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
  char num[64];
  auto g = [&](double x) {   // "30.0", "1.5", "0.71539441": whole numbers keep their ".0"
    std::snprintf(num, sizeof num, "%.8g", x);
    std::string r = num;
    if (r.find_first_of(".e") == std::string::npos) r += ".0";
    return r;
  };
  const std::string angle = angles ? "angle_style cosine\nangle_coeff 1 " + g(o.k_theta * ep) + "\n" : std::string();
  const std::string base = data.substr(0, data.size() - 5);
  const double rc = 1.122462048309373 * sg;   // 2^(1/6) σ
  std::string head = "# Kremer-Grest melt written by CAPS: " + std::to_string(o.chains) + " x " + std::to_string(o.beads) + " beads, rho sigma^3 = " + g(o.density) +
                     ", k_theta = " + g(o.k_theta) + "\n# Kremer & Grest, J. Chem. Phys. 92, 5057 (1990); push-off after Auhl et al., J. Chem. Phys. 119, 12718 (2003)\n";
  if (real)
    head += "# units real, mapped as given: sigma = " + g(o.sigma) + " A, epsilon = k_B T = " + g(u.eps) + " kcal/mol (T = " + g(o.temperature) + " K), m = " + g(o.bead_mass) +
            " g/mol, tau = sigma sqrt(m/epsilon) = " + g(u.tau_fs) + " fs, density " + g(u.density) + " g/cm3\n";
  char b[4096];
  std::snprintf(b, sizeof b,
                "%s"
                "units %s\natom_style angle\nboundary p p p\nread_data %s\n"
                "bond_style fene\nbond_coeff 1 %s %s %s %s\nspecial_bonds fene\n%s"
                "comm_modify cutoff %s\nneighbor %s bin\nneigh_modify every 1 delay 0 check yes\ntimestep %s\n"
                "velocity all create %s %llu dist gaussian\nfix lang all langevin %s %s %s %llu\nthermo_style custom step temp pe ebond press\nthermo 1000\n\n"
                "# 1. push-off: a soft repulsion ramped up, displacements limited, so random-walk overlaps open without breaking bonds\n"
                "pair_style soft %s\npair_coeff * * 0.0\nvariable prefactor equal ramp(0,%s)\nfix push all adapt 1 pair soft a * * v_prefactor\n"
                "fix move all nve/limit %s\nrun %.0f\nunfix push\n\n"
                "# 2. Kremer-Grest: WCA repulsion (LJ cut at 2^(1/6) sigma, shifted); still limited for a while, then free\n"
                "pair_style lj/cut %s\npair_modify shift yes\npair_coeff * * %s %s %s\n"
                "run 2000\nunfix move\nfix move all nve\n"
                "dump d all custom 10000 %s.lammpstrj id mol type xu yu zu\nrun %.0f\nwrite_data %s_equilibrated.data pair ij\n",
                head.c_str(), real ? "real" : "lj", data.c_str(), g(30.0 * ep / (sg * sg)).c_str(), g(1.5 * sg).c_str(), g(ep).c_str(), g(sg).c_str(), angle.c_str(),
                g(2.0 * sg).c_str(), g(0.4 * sg).c_str(), g(0.01 * tu).c_str(), g(T).c_str(), (unsigned long long)(o.seed + 1), g(T).c_str(), g(T).c_str(), g(1.0 * tu).c_str(),
                (unsigned long long)(o.seed + 2), g(rc).c_str(), g(60.0 * ep).c_str(), g(0.02 * sg).c_str(), pushoff_steps, g(rc).c_str(), g(ep).c_str(), g(sg).c_str(), g(rc).c_str(),
                base.c_str(), run_steps, base.c_str());
  f << b;
}


void write_kg_bond_swap(const System& s, const KgOptions& o, const BondSwapExport& e, const std::string& stem) {
  const Cell& c = s.cell;
  if (!c.valid() || c.a[1] != 0 || c.a[2] != 0 || c.b[0] != 0 || c.b[2] != 0 || c.c[0] != 0 || c.c[1] != 0)
    throw std::invalid_argument("bond/swap export: an orthogonal periodic cell is needed");
  if (e.every < 1 || !(e.fraction > 0 && e.fraction <= 1) || !(e.cutoff > 0)) throw std::invalid_argument("bond/swap export: every ≥ 1, 0 < fraction ≤ 1, cutoff > 0");
  const size_t n = s.atoms.size();
  std::vector<char> in(n, e.group.empty() ? 1 : 0);
  for (uint32_t a : e.group) {
    if (a >= n) throw std::invalid_argument("bond/swap export: the group names an atom the structure does not have");
    in[a] = 1;
  }
  const KgUnits u = kg_units(o);
  const bool real = u.eps > 0;
  const double sg = real ? o.sigma : 1.0, ep = real ? u.eps : 1.0, ms = real ? o.bead_mass : 1.0, tu = real ? u.tau_fs : 1.0, T = real ? o.temperature : 1.0;
  std::vector<int64_t> mol(n, 2);
  if (e.stickers) {
    for (size_t i = 0; i < n; ++i)
      if (in[i]) mol[i] = 1;
  } else {
    mol = chain_position_labels(s, e.from_both_ends);
    for (size_t i = 0; i < n; ++i)
      if (mol[i] == 0) throw std::invalid_argument("bond/swap export: chain mode needs linear chains (atom " + std::to_string(i + 1) + " is not on one)");
  }
  // angles from the bond graph (none across a sticker–sticker bond in sticker mode)
  std::vector<std::vector<uint32_t>> nb(n);
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  auto sticker_bond = [&](uint32_t i, uint32_t j) { return e.stickers && in[i] && in[j]; };
  std::vector<std::array<uint32_t, 3>> angles;
  if (o.k_theta > 0)
    for (uint32_t j = 0; j < n; ++j)
      for (size_t a = 0; a < nb[j].size(); ++a)
        for (size_t b = a + 1; b < nb[j].size(); ++b)
          if (!sticker_bond(nb[j][a], j) && !sticker_bond(j, nb[j][b])) angles.push_back({nb[j][a], j, nb[j][b]});
  char b[512];
  {
    std::ofstream f(stem + ".data");
    if (!f) throw std::runtime_error("cannot write " + stem + ".data");
    f << "Kremer-Grest system for fix bond/swap (CAPS" << (real ? ", units real" : ", reduced units") << "; molecule IDs: "
      << (e.stickers ? "1 for every sticker" : "place along the chain") << ")\n\n";
    f << n << " atoms\n" << s.bonds.size() << " bonds\n" << angles.size() << " angles\n2 atom types\n1 bond types\n1 angle types\n";
    std::snprintf(b, sizeof b, "\n0 %.6f xlo xhi\n0 %.6f ylo yhi\n0 %.6f zlo zhi\n\nMasses\n\n1 %.6g\n2 %.6g\n\nAtoms # angle\n\n", c.a[0] * sg, c.b[1] * sg, c.c[2] * sg, ms, ms);
    f << b;
    for (size_t i = 0; i < n; ++i) {
      const Vec3 w = c.wrap(s.atoms[i].pos) * sg;
      const Vec3 f3 = c.to_fractional(s.atoms[i].pos);
      std::snprintf(b, sizeof b, "%zu %lld %d %.6f %.6f %.6f %d %d %d\n", i + 1, (long long)mol[i], in[i] ? 2 : 1, w[0], w[1], w[2], int(std::floor(f3[0])),
                    int(std::floor(f3[1])), int(std::floor(f3[2])));
      f << b;
    }
    f << "\nBonds\n\n";
    for (size_t k = 0; k < s.bonds.size(); ++k) f << k + 1 << " 1 " << s.bonds[k].i + 1 << " " << s.bonds[k].j + 1 << "\n";
    if (!angles.empty()) {
      f << "\nAngles\n\n";
      for (size_t k = 0; k < angles.size(); ++k) f << k + 1 << " 1 " << angles[k][0] + 1 << " " << angles[k][1] + 1 << " " << angles[k][2] + 1 << "\n";
    }
  }
  std::ofstream f(stem + ".in");
  if (!f) throw std::runtime_error("cannot write " + stem + ".in");
  const std::string data = stem.substr(stem.find_last_of("/\\") == std::string::npos ? 0 : stem.find_last_of("/\\") + 1) + ".data";
  const std::string base = data.substr(0, data.size() - 5);
  auto g = [&](double x) {
    char num[64];
    std::snprintf(num, sizeof num, "%.8g", x);
    std::string r = num;
    if (r.find_first_of(".e") == std::string::npos) r += ".0";
    return r;
  };
  f << "# Bond exchange in a Kremer-Grest system, written by CAPS: LAMMPS fix bond/swap (MC package)\n"
       "# Sides, Grest, Stevens & Plimpton, J. Polym. Sci. B 42, 199 (2004); Kremer & Grest, J. Chem. Phys. 92, 5057 (1990)\n"
       "# a swap A1-A2 + B1-B2 -> A1-B2 + B1-A2 needs all four atoms in group exch, A1 and B1 with the same molecule ID, A2 and\n"
       "# B2 likewise, and the four distances A1-A2, B1-B2, A1-B2, B1-A2 below the cutoff; Metropolis on pair, bond and angle energy\n";
  if (e.stickers) f << "# sticker mode: every sticker has molecule ID 1; bonds stored on both atoms\nnewton on off\n";
  else f << "# chain mode: molecule ID = place along the chain" << (e.from_both_ends ? ", counted from the nearer end" : "") << ", so swaps keep chain lengths\n";
  f << "units " << (real ? "real" : "lj") << "\natom_style angle\nboundary p p p\nread_data " << data << "\n";
  f << "bond_style fene\nbond_coeff 1 " << g(30.0 * ep / (sg * sg)) << " " << g(1.5 * sg) << " " << g(ep) << " " << g(sg) << "\nspecial_bonds fene\n";
  f << "angle_style cosine\nangle_coeff 1 " << g(o.k_theta * ep) << "\n";
  f << "pair_style lj/cut " << g(1.122462048309373 * sg) << "\npair_modify shift yes\npair_coeff * * " << g(ep) << " " << g(sg) << "\n";
  f << "neighbor " << g(0.4 * sg) << " bin\nneigh_modify every 1 delay 0 check yes\ncomm_modify cutoff " << g(std::max(2.0, e.cutoff + 0.5) * sg) << "\n";
  f << "group exch type 2\ntimestep " << g(0.01 * tu) << "\n";
  f << "velocity all create " << g(T) << " " << (unsigned long long)(e.seed + 1) << " dist gaussian\n";
  f << "fix move all nve\nfix lang all langevin " << g(T) << " " << g(T) << " " << g(1.0 * tu) << " " << (unsigned long long)(e.seed + 2) << "\n";
  f << "fix swap exch bond/swap " << e.every << " " << g(e.fraction) << " " << g(e.cutoff * sg) << " " << (unsigned long long)(e.seed + 3) << "\n";
  f << "thermo_style custom step temp pe ebond f_swap[1] f_swap[2]\nthermo 1000\n";
  f << "dump d all custom 10000 " << base << "_swap.lammpstrj id mol type xu yu zu\n";
  std::snprintf(b, sizeof b, "run %.0f\nwrite_data %s_swapped.data\n", e.run_steps, base.c_str());
  f << b;
}

}  // namespace caps
