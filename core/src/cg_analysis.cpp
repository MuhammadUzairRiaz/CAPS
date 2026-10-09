// CAPS entanglement, mechanics and dynamics analyses of coarse-grained melts (see caps/cg_analysis.hpp).
#include "caps/cg_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "caps/entangle.hpp"

namespace caps {

namespace {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// each chain's beads made whole along the chain (minimum image from its first bead)
std::vector<std::vector<Vec3>> whole_chains(const CgTopology& t, const std::vector<Vec3>& pos, const Cell& cell) {
  std::vector<std::vector<Vec3>> out;
  const bool pbc = cell.valid();
  for (const auto& ch : t.chains) {
    std::vector<Vec3> p(ch.size());
    p[0] = pos[size_t(ch[0])];
    for (size_t i = 1; i < ch.size(); ++i) {
      const Vec3 d = pos[size_t(ch[i])] - pos[size_t(ch[i - 1])];
      p[i] = p[i - 1] + (pbc ? cell.minimum_image(d) : d);
    }
    out.push_back(std::move(p));
  }
  return out;
}

double lin_slope(const std::vector<double>& x, const std::vector<double>& y, double* intercept = nullptr, double* rms = nullptr) {
  const size_t n = x.size();
  if (n < 2) return kNaN;
  double mx = 0, my = 0;
  for (size_t i = 0; i < n; ++i) mx += x[i], my += y[i];
  mx /= double(n), my /= double(n);
  double sxy = 0, sxx = 0;
  for (size_t i = 0; i < n; ++i) sxy += (x[i] - mx) * (y[i] - my), sxx += (x[i] - mx) * (x[i] - mx);
  const double b = sxx > 0 ? sxy / sxx : kNaN, a = my - b * mx;
  if (intercept) *intercept = a;
  if (rms) {
    double e = 0;
    for (size_t i = 0; i < n; ++i) e += (y[i] - a - b * x[i]) * (y[i] - a - b * x[i]);
    *rms = std::sqrt(e / double(n));
  }
  return b;
}
}  // namespace

// ------------------------------------------------------------------------------------------------ Z1
void write_z1_config(const CgTopology& t, const std::vector<Vec3>& pos, const Cell& cell, const std::string& path) {
  if (!cell.valid()) throw std::invalid_argument("Z1 needs a periodic box");
  if (std::fabs(cell.b[0]) > 1e-9 || std::fabs(cell.c[0]) > 1e-9 || std::fabs(cell.c[1]) > 1e-9) throw std::invalid_argument("Z1: an orthogonal box is needed");
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot write " + path);
  char b[160];
  f << t.chains.size() << "\n";
  std::snprintf(b, sizeof b, "%.6f %.6f %.6f\n", cell.a[0], cell.b[1], cell.c[2]);
  f << b;
  for (size_t k = 0; k < t.chains.size(); ++k) f << (k ? " " : "") << t.chains[k].size();
  f << "\n";
  for (const auto& ch : whole_chains(t, pos, cell))
    for (const auto& p : ch) {
      std::snprintf(b, sizeof b, "%.6f %.6f %.6f\n", p[0] - cell.origin[0], p[1] - cell.origin[1], p[2] - cell.origin[2]);
      f << b;
    }
}

CgChainPaths read_z1_paths(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::invalid_argument("cannot open " + path);
  std::vector<std::string> lines;
  for (std::string l; std::getline(f, l);) {
    if (l.find_first_not_of(" \t\r") != std::string::npos) lines.push_back(l);
  }
  if (lines.size() < 3) throw std::invalid_argument(path + ": too short for Z1 paths");
  const int nch = std::stoi(lines[0]);
  auto ints = [](const std::string& s) { std::istringstream ss(s); std::vector<double> v; double x; while (ss >> x) v.push_back(x); return v; };
  CgChainPaths P;
  size_t li = 2;
  std::vector<int> counts;
  const auto third = ints(lines[2]);
  const bool all_on_third = int(third.size()) == nch && nch > 1;
  if (all_on_third) { for (double x : third) counts.push_back(int(x)); li = 3; }
  for (int c = 0; c < nch; ++c) {
    int n;
    if (all_on_third) n = counts[size_t(c)];
    else {
      if (li >= lines.size()) throw std::invalid_argument(path + ": ends before chain " + std::to_string(c + 1));
      n = int(ints(lines[li++]).at(0));
    }
    std::vector<Vec3> nodes;
    for (int k = 0; k < n; ++k) {
      if (li >= lines.size()) throw std::invalid_argument(path + ": ends inside chain " + std::to_string(c + 1));
      const auto v = ints(lines[li++]);
      if (v.size() < 3) throw std::invalid_argument(path + ": a node line without x y z");
      nodes.push_back({v[0], v[1], v[2]});
    }
    double L = 0;
    for (size_t k = 1; k < nodes.size(); ++k) L += norm(nodes[k] - nodes[k - 1]);
    P.lpp.push_back(L);
    P.r2.push_back(nodes.size() > 1 ? dot(nodes.back() - nodes.front(), nodes.back() - nodes.front()) : 0.0);
    P.z.push_back(std::max(0, n - 2));
    P.beads.push_back(0);
  }
  return P;
}

CgChainPaths paths_from_positions(const CgTopology& t, const std::vector<Vec3>& start, const std::vector<Vec3>& after, const Cell& cell) {
  const auto a = whole_chains(t, start, cell), b = whole_chains(t, after, cell);
  CgChainPaths P;
  for (size_t c = 0; c < a.size(); ++c) {
    double L = 0;
    for (size_t k = 1; k < b[c].size(); ++k) L += norm(b[c][k] - b[c][k - 1]);
    const Vec3 R = a[c].back() - a[c].front();
    P.lpp.push_back(L);
    P.r2.push_back(dot(R, R));
    P.beads.push_back(int(a[c].size()));
  }
  return P;
}

CgChainPaths caps_ppa(const CgTopology& t, const std::vector<Vec3>& pos, const Cell& cell, double sigma, int max_steps) {
  System s;
  s.cell = cell;
  const auto whole = whole_chains(t, pos, cell);
  std::vector<std::vector<uint32_t>> bb;
  for (size_t c = 0; c < t.chains.size(); ++c) {
    std::vector<uint32_t> ids;
    for (size_t k = 0; k < t.chains[c].size(); ++k) {
      Atom a;
      a.id = int64_t(s.atoms.size() + 1);
      a.mol = int64_t(c + 1);
      a.element = 6;   // a bead drawn as carbon; PPA uses only positions and bonds
      a.pos = whole[c][k];
      ids.push_back(uint32_t(s.atoms.size()));
      s.atoms.push_back(a);
      if (k > 0) s.bonds.push_back({ids[k - 1], ids[k], 1});
    }
    bb.push_back(ids);
  }
  s.has_mol = true;
  PrimitivePathOptions o;
  o.sigma = sigma;
  o.max_steps = max_steps;
  const PrimitivePaths pp = primitive_paths(s, bb, o);
  CgChainPaths P;
  P.lpp = pp.lpp, P.r2 = pp.r2;
  P.steps = pp.steps, P.converged = pp.converged, P.sigma = pp.sigma;
  for (const auto& c : t.chains) P.beads.push_back(int(c.size()));
  return P;
}

CgEntanglement entanglement_of(const CgChainPaths& p) {
  CgEntanglement e;
  const size_t n = p.lpp.size();
  if (!n) return e;
  e.chains = int(n);
  double N = 0, r2 = 0, l = 0, l2 = 0, z = 0;
  for (size_t c = 0; c < n; ++c) N += c < p.beads.size() ? p.beads[c] : 0, r2 += p.r2[c], l += p.lpp[c], l2 += p.lpp[c] * p.lpp[c];
  e.N = N / double(n), e.r2 = r2 / double(n), e.lpp = l / double(n), e.lpp2 = l2 / double(n);
  e.a_pp = e.lpp > 0 ? e.r2 / e.lpp : 0;
  if (e.N > 1 && e.lpp > 0) e.ne_s_coil = (e.N - 1) * e.r2 / (e.lpp * e.lpp);
  if (e.N > 1 && e.r2 > 0 && e.lpp2 / e.r2 > 1) e.ne_mod_s_coil = (e.N - 1) / (e.lpp2 / e.r2 - 1);
  if (!p.z.empty() && p.z.size() == n) {
    for (int k : p.z) z += k;
    e.z = z / double(n);
    if (e.N > 1) e.ne_s_kink = e.N * (e.N - 1) / (e.z * (e.N - 1) + e.N);
    if (e.z > 0) e.ne_mod_s_kink = e.N / e.z;
  }
  return e;
}

CgMultiEstimate multi_estimators(const std::vector<CgEntanglement>& sets, const std::vector<double>& internal, double l0sq) {
  CgMultiEstimate m;
  std::vector<double> N, Z, Y;
  for (const auto& s : sets) {
    N.push_back(s.N);
    if (s.z >= 0) Z.push_back(s.z);
    if (l0sq > 0 && s.N > 1) Y.push_back(s.lpp * s.lpp / ((s.N - 1) * l0sq));
  }
  if (sets.size() < 2) { m.note = "the M-estimators need two chain lengths or more"; return m; }
  if (Z.size() == sets.size()) {
    const double g = lin_slope(N, Z);
    if (g > 0) m.ne_m_kink = 1.0 / g;
  } else m.note = "M-kink needs ⟨Z⟩ (Z1+) for every chain length";
  if (Y.size() == sets.size() && !internal.empty() && l0sq > 0) {
    const double slope = lin_slope(N, Y);
    // C(x)/x with C(x) = ⟨R²(x − 1 bonds)⟩ / ((x − 1) l₀²) = internal[x − 2] / l₀²; the first x where it falls to the slope
    if (slope > 0)
      for (size_t n = 1; n < internal.size(); ++n) {
        const double x = double(n + 1), cx = internal[n - 1] / l0sq;
        if (cx / x <= slope) { m.ne_m_coil = x; break; }
      }
    if (m.ne_m_coil == 0) m.note += (m.note.empty() ? "" : "; ") + std::string("M-coil: C(x)/x does not fall to d(⟨L_pp⟩²/R²_RW)/dN within the chains measured");
  }
  return m;
}

std::string ppa_ends(const CgTopology& t) {
  std::string out;
  std::vector<int> ids;
  for (const auto& c : t.chains) { ids.push_back(c.front() + 1); if (c.size() > 1) ids.push_back(c.back() + 1); }
  for (size_t k = 0; k < ids.size(); k += 200) {
    out += "group ends id";
    for (size_t j = k; j < std::min(ids.size(), k + 200); ++j) out += " " + std::to_string(ids[j]);
    out += "\n";
  }
  return out;
}

std::string ppa_deck(const CgTopology& t, double sigma, const std::string& ends_file) {
  (void)t;
  char b[2400];
  std::snprintf(b, sizeof b,
                "# CAPS primitive-path analysis (LAMMPS; Everaers et al., Science 303, 823 (2004)): chain ends held, pairs within a\n"
                "# chain off, FENE bonds pulling with zero rest length, chains kept from crossing by a WCA of diameter %.3f Å (ε = 1\n"
                "# kcal/mol, the bonds' K = 30 ε/σ² — much stiffer bonds would pull chains through each other); quenched\n"
                "#   lmp -in in.cg_ppa -var DATA equil.data -var OUT ppa;  caps ppa MAP equil.data ppa.lammpstrj …\n"
                "variable SIGMA index %.4f\nvariable QUENCH index 50000\nunits real\natom_style full\nboundary p p p\nspecial_bonds lj 0 1 1\nread_data ${DATA}\n"
                "angle_style none\ndihedral_style none\nimproper_style none\n"
                "# Everaers et al.'s bonds, as CAPS's own PPA: FENE with K = 30 ε/σ², R0 = 1.5 σ and no rest length (its own WCA off)\n"
                "bond_style fene\nbond_coeff * $(30.0/v_SIGMA^2) $(1.5*v_SIGMA) 0.0 ${SIGMA}\n"
                "pair_style lj/cut $(v_SIGMA*2.0^(1.0/6.0))\npair_coeff * * 1.0 ${SIGMA}\npair_modify shift yes\n"
                "neigh_modify exclude molecule/intra all one 10000\ncomm_modify cutoff 30.0\ninclude %s\n"
                "thermo_style custom step pe ebond epair\nthermo 1000\ntimestep 1.0\n"
                "# an overdamped quench with limited steps (the bonds pull hard at first), the ends held after the thermostat's forces\n"
                "fix lim all nve/limit 0.05\nfix cold all langevin 0.001 0.001 100.0 4928459\nfix hold ends setforce 0.0 0.0 0.0\nvelocity all set 0.0 0.0 0.0\n"
                "run ${QUENCH}\nunfix lim\nunfix cold\n"
                "min_style cg\nmin_modify dmax 0.1\nminimize 1.0e-8 1.0e-10 100000 1000000\n"
                "write_dump all custom ${OUT}.lammpstrj id mol type xu yu zu modify sort id\n",
                sigma, sigma, ends_file.c_str());
  return b;
}

// ------------------------------------------------------------------------------------------------ tension
std::string tension_deck(const CgTensionDeck& d) {
  if (d.mode != "stress" && d.mode != "volume") throw std::invalid_argument("tension mode: stress or volume");
  const char* special = d.exclude <= 2 ? "0 1 1" : d.exclude == 3 ? "0 0 1" : "0 0 0";
  const double steps = d.max_strain / (d.rate * d.dt);
  const long long nsteps = (long long)std::llround(steps);
  const long long dump = std::max(1LL, (long long)std::llround(steps / std::max(1, d.frames)));
  const long long pe = d.print_every > 0 ? d.print_every : std::max(1LL, nsteps / 500);
  char b[5200];
  std::snprintf(b, sizeof b,
                "# CAPS uniaxial tension (%s mode), z at an engineering rate of %g /fs to strain %g (LAMMPS)\n"
                "#   lmp -in in.cg_tensile -var DATA equil.data -var BONDED ../bonded -var PAIR ../lj -var OUT %s_%g\n"
                "# σ = −(Pzz − (Pxx + Pyy)/2) in MPa and its parts (stress/atom per term, summed): %s\n"
                "variable T index %g\nvariable P index %g\nvariable DT index %g\nvariable RATE index %g\nvariable STEPS index %lld\nvariable DUMP index %lld\nvariable PRINT index %lld\n"
                "units real\natom_style full\nboundary p p p\nspecial_bonds lj %s\nread_data ${DATA}\ninclude ${BONDED}/bonded.in\ninclude ${PAIR}/pair.in\n"
                "neighbor 3.0 bin\nneigh_modify delay 0 every 1 check yes\ntimestep ${DT}\n"
                "velocity all create ${T} 4928459 mom yes rot yes dist gaussian\n"
                "# the start relaxed at P first (no stress before the pull)\n"
                "variable EQ index 50000\nfix eq all npt temp ${T} ${T} $(100.0*dt) iso ${P} ${P} $(1000.0*dt)\nrun ${EQ}\nunfix eq\nreset_timestep 0\n"
                "variable lz0 equal $(lz)\nvariable strain equal (lz-v_lz0)/v_lz0\n"
                "variable mpa equal 0.101325\n"
                "variable sig equal -(pzz-0.5*(pxx+pyy))*v_mpa\n"
                "compute sb all stress/atom NULL bond\ncompute sa all stress/atom NULL angle\ncompute sd all stress/atom NULL dihedral\n"
                "compute sp all stress/atom NULL pair\ncompute sk all stress/atom NULL ke\n"
                "compute rb all reduce sum c_sb[1] c_sb[2] c_sb[3]\ncompute ra all reduce sum c_sa[1] c_sa[2] c_sa[3]\ncompute rd all reduce sum c_sd[1] c_sd[2] c_sd[3]\n"
                "compute rp all reduce sum c_sp[1] c_sp[2] c_sp[3]\ncompute rk all reduce sum c_sk[1] c_sk[2] c_sk[3]\n"
                "# per-atom stress is −pressure × volume: the axial stress of a part is (s_zz − (s_xx + s_yy)/2) / V\n"
                "variable sbond equal (c_rb[3]-0.5*(c_rb[1]+c_rb[2]))/vol*v_mpa\nvariable sangle equal (c_ra[3]-0.5*(c_ra[1]+c_ra[2]))/vol*v_mpa\n"
                "variable sdih equal (c_rd[3]-0.5*(c_rd[1]+c_rd[2]))/vol*v_mpa\nvariable spair equal (c_rp[3]-0.5*(c_rp[1]+c_rp[2]))/vol*v_mpa\n"
                "variable ske equal (c_rk[3]-0.5*(c_rk[1]+c_rk[2]))/vol*v_mpa\n"
                "%s"
                "thermo_style custom step temp press density pxx pyy pzz lx ly lz\nthermo ${PRINT}\n"
                "fix out all print ${PRINT} \"$(v_strain:%%.6f) $(v_sig:%%.4f) $(v_sbond:%%.4f) $(v_sangle:%%.4f) $(v_sdih:%%.4f) $(v_spair:%%.4f) $(v_ske:%%.4f) $(density:%%.5f)\" "
                "file ${OUT}.stress_strain.dat title \"# strain sigma_MPa bond angle dihedral pair kinetic density\" screen no\n"
                "dump frames all custom ${DUMP} ${OUT}.lammpstrj id mol type xu yu zu\ndump_modify frames sort id\n"
                "log ${OUT}.log\nrun ${STEPS}\nwrite_data ${OUT}.data\n",
                d.mode.c_str(), d.rate, d.max_strain, d.mode.c_str(), d.rate,
                d.mode == "stress" ? "lateral axes at P (npt, x and y uncoupled)" : "volume held (x and y shrink with z)", d.T, d.P, d.dt, d.rate, nsteps, dump, pe, special,
                d.mode == "stress" ? "fix md all npt temp ${T} ${T} $(100.0*dt) x ${P} ${P} $(1000.0*dt) y ${P} ${P} $(1000.0*dt) couple none\n"
                                     "fix pull all deform 1 z erate ${RATE} remap x\n"
                                   : "fix md all nvt temp ${T} ${T} $(100.0*dt)\nfix pull all deform 1 z erate ${RATE} x volume y volume remap x\n");
  return b;
}

CgStressStrain read_stress_strain(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::invalid_argument("cannot open " + path);
  CgStressStrain c;
  for (std::string l; std::getline(f, l);) {
    if (l.empty() || l[0] == '#') continue;
    std::istringstream ss(l);
    std::vector<double> v;
    double x;
    while (ss >> x) v.push_back(x);
    if (v.size() < 2) continue;
    c.strain.push_back(v[0]), c.stress.push_back(v[1]);
    auto at = [&](size_t k) { return k < v.size() ? v[k] : kNaN; };
    c.bond.push_back(at(2)), c.angle.push_back(at(3)), c.dihedral.push_back(at(4)), c.pair.push_back(at(5)), c.kinetic.push_back(at(6)), c.density.push_back(at(7));
  }
  if (c.strain.size() < 5) throw std::invalid_argument(path + ": fewer than five stress lines");
  return c;
}

CgTensionAnalysis analyse_tension(const CgStressStrain& c, double fit_strain, double hardening_from) {
  CgTensionAnalysis a;
  const size_t n = c.strain.size();
  // modulus: linear fit from 0 to fit_strain
  std::vector<double> x, y;
  for (size_t i = 0; i < n; ++i) if (c.strain[i] <= fit_strain) x.push_back(c.strain[i]), y.push_back(c.stress[i]);
  if (x.size() >= 3) a.modulus = lin_slope(x, y);
  // smoothed stress: ±0.5 % strain (as CAPS's own tensile analysis), at least one neighbour each side
  std::vector<double> sm(n);
  for (size_t i = 0; i < n; ++i) {
    double s = 0;
    int k = 0;
    for (size_t j = 0; j < n; ++j)
      if (std::fabs(c.strain[j] - c.strain[i]) <= 0.005 || (j + 1 >= i && j <= i + 1)) s += c.stress[j], ++k;
    sm[i] = s / k;
  }
  // yield: the running maximum when the smoothed stress first falls 2 % below it (past the linear part)
  size_t iy = 0, imax = 0;
  double runmax = -1e300;
  for (size_t i = 0; i < n; ++i) {
    if (c.strain[i] < fit_strain) continue;
    if (sm[i] > runmax) runmax = sm[i], imax = i;
    else if (sm[i] < runmax - 0.02 * std::fabs(runmax)) { iy = imax; break; }
  }
  if (iy) {
    a.yield_stress = sm[iy], a.yield_strain = c.strain[iy];
    size_t imin = iy;
    for (size_t i = iy; i < n; ++i) {
      if (sm[i] < sm[imin]) imin = i;
      if (i > imin && sm[i] > sm[imin] + 0.25 * (a.yield_stress - sm[imin])) break;   // past the minimum, rising again
    }
    a.softening = a.yield_stress - sm[imin];
    a.min_strain = c.strain[imin];
  } else a.note = "no stress maximum (no yield) in the range";
  // hardening: σ against λ² − 1/λ past the softening minimum (or from the strain given)
  const double from = hardening_from >= 0 ? hardening_from : (iy ? a.min_strain : 0.3);
  std::vector<double> gx, gy;
  for (size_t i = 0; i < n; ++i)
    if (c.strain[i] >= from) {
      const double lam = 1 + c.strain[i];
      gx.push_back(lam * lam - 1 / lam), gy.push_back(c.stress[i]);
    }
  if (gx.size() >= 5) {
    double rms;
    a.hardening_modulus = lin_slope(gx, gy, nullptr, &rms);
    a.hardening_rms = rms;
    a.hardening_from = from;
  } else a.note += (a.note.empty() ? "" : "; ") + std::string("too little strain past the softening for a hardening fit");
  return a;
}

CgOrientation orientation_of(const CgTopology& t, const std::vector<Vec3>& pos, const Cell& cell, double probe) {
  CgOrientation o;
  double p2 = 0;
  long nb = 0;
  const bool pbc = cell.valid();
  for (const auto& [i, j] : t.bonds) {
    Vec3 d = pos[size_t(j)] - pos[size_t(i)];
    if (pbc) d = cell.minimum_image(d);
    const double c = d[2] / norm(d);
    p2 += 1.5 * c * c - 0.5, ++nb;
  }
  o.p2 = nb ? p2 / double(nb) : 0;
  double rx = 0, ry = 0, rz = 0;
  const auto whole = whole_chains(t, pos, cell);
  for (const auto& ch : whole) {
    const Vec3 R = ch.back() - ch.front();
    rx += R[0] * R[0], ry += R[1] * R[1], rz += R[2] * R[2];
  }
  o.ree_anisotropy = (rx + ry) > 0 ? rz / (0.5 * (rx + ry)) : 0;
  if (probe > 0 && pbc) {
    // a 1 Å grid: the share of points farther than `probe` from every bead (cell list)
    const double L[3] = {cell.a[0], cell.b[1], cell.c[2]};
    int ng[3], nc[3];
    for (int k = 0; k < 3; ++k) ng[k] = std::max(1, int(L[k])), nc[k] = std::max(1, int(L[k] / std::max(probe, 1.0)));
    std::vector<std::vector<uint32_t>> bins(size_t(nc[0] * nc[1] * nc[2]));
    std::vector<Vec3> w(pos.size());
    auto cidx = [&](int x, int y, int z) { return size_t(((x % nc[0] + nc[0]) % nc[0]) + nc[0] * (((y % nc[1] + nc[1]) % nc[1]) + nc[1] * ((z % nc[2] + nc[2]) % nc[2]))); };
    for (size_t i = 0; i < pos.size(); ++i) {
      w[i] = cell.wrap(pos[i]);
      int q[3];
      for (int k = 0; k < 3; ++k) q[k] = std::min(nc[k] - 1, int((w[i][size_t(k)] - cell.origin[size_t(k)]) / L[k] * nc[k]));
      bins[cidx(q[0], q[1], q[2])].push_back(uint32_t(i));
    }
    long empty = 0, total = 0;
    const int reach = int(std::ceil(probe / (std::min({L[0] / nc[0], L[1] / nc[1], L[2] / nc[2]})))) ;
    for (int gx = 0; gx < ng[0]; ++gx)
      for (int gy = 0; gy < ng[1]; ++gy)
        for (int gz = 0; gz < ng[2]; ++gz) {
          const Vec3 p{cell.origin[0] + (gx + 0.5) * L[0] / ng[0], cell.origin[1] + (gy + 0.5) * L[1] / ng[1], cell.origin[2] + (gz + 0.5) * L[2] / ng[2]};
          const int c0 = int((p[0] - cell.origin[0]) / L[0] * nc[0]), c1 = int((p[1] - cell.origin[1]) / L[1] * nc[1]), c2 = int((p[2] - cell.origin[2]) / L[2] * nc[2]);
          bool hit = false;
          for (int dx = -reach; dx <= reach && !hit; ++dx)
            for (int dy = -reach; dy <= reach && !hit; ++dy)
              for (int dz = -reach; dz <= reach && !hit; ++dz)
                for (uint32_t i : bins[cidx(c0 + dx, c1 + dy, c2 + dz)]) {
                  const Vec3 d = cell.minimum_image(w[i] - p);
                  if (dot(d, d) < probe * probe) { hit = true; break; }
                }
          empty += !hit, ++total;
        }
    o.void_fraction = double(empty) / double(total);
  }
  return o;
}

// ------------------------------------------------------------------------------------------------ dynamics
CgDynamics cg_dynamics(const CgTopology& t, const std::vector<std::vector<Vec3>>& frames, const std::vector<Cell>& cells, const std::vector<double>& times) {
  CgDynamics D;
  const size_t F = frames.size();
  if (F < 3) throw std::invalid_argument("dynamics need three frames or more");
  // chains whole in each frame, then unwrapped in time (each chain's first bead follows the minimum image between frames)
  std::vector<std::vector<std::vector<Vec3>>> W(F);
  for (size_t f = 0; f < F; ++f) W[f] = whole_chains(t, frames[f], f < cells.size() ? cells[f] : Cell{});
  for (size_t f = 1; f < F; ++f) {
    const Cell& c = f < cells.size() ? cells[f] : Cell{};
    if (!c.valid()) continue;
    for (size_t ch = 0; ch < W[f].size(); ++ch) {
      const Vec3 jump = W[f][ch][0] - W[f - 1][ch][0];
      const Vec3 shift = c.minimum_image(jump) - jump;
      for (auto& p : W[f][ch]) p = p + shift;
    }
  }
  const size_t nch = t.chains.size();
  std::vector<std::vector<Vec3>> com(F, std::vector<Vec3>(nch));
  for (size_t f = 0; f < F; ++f)
    for (size_t ch = 0; ch < nch; ++ch) {
      Vec3 s{0, 0, 0};
      for (const auto& p : W[f][ch]) s = s + p;
      com[f][ch] = s * (1.0 / double(W[f][ch].size()));
    }
  // lags: every frame step up to 10, then logarithmic
  std::vector<size_t> lags;
  for (size_t k = 1; k < F; k = k < 10 ? k + 1 : size_t(std::ceil(double(k) * 1.2))) lags.push_back(k);
  double ree0 = 0;
  long nr = 0;
  for (size_t f = 0; f < F; ++f)
    for (size_t ch = 0; ch < nch; ++ch) { const Vec3 R = W[f][ch].back() - W[f][ch].front(); ree0 += dot(R, R), ++nr; }
  ree0 /= double(nr);
  for (size_t lag : lags) {
    double g1 = 0, g2 = 0, g3 = 0, p1 = 0, rr = 0;
    long n1 = 0, n3 = 0, nb = 0;
    for (size_t f = 0; f + lag < F; ++f) {
      const size_t h = f + lag;
      for (size_t ch = 0; ch < nch; ++ch) {
        const auto& a = W[f][ch];
        const auto& b = W[h][ch];
        const size_t L = a.size(), lo = L / 4, hi = L - L / 4;
        const Vec3 dc = com[h][ch] - com[f][ch];
        for (size_t k = lo; k < std::max(hi, lo + 1) && k < L; ++k) {
          const Vec3 d = b[k] - a[k];
          g1 += dot(d, d);
          const Vec3 e = d - dc;
          g2 += dot(e, e);
          ++n1;
        }
        g3 += dot(dc, dc), ++n3;
        for (size_t k = 0; k + 1 < L; ++k) {
          const Vec3 u = a[k + 1] - a[k], v = b[k + 1] - b[k];
          p1 += dot(u, v) / (norm(u) * norm(v)), ++nb;
        }
        rr += dot(a.back() - a.front(), b.back() - b.front());
      }
    }
    D.t.push_back(times[lag] - times[0]);
    D.g1.push_back(g1 / double(n1)), D.g2.push_back(g2 / double(n1)), D.g3.push_back(g3 / double(n3));
    D.p1.push_back(p1 / double(nb)), D.ree.push_back(rr / double(n3) / ree0);
  }
  // D from g₃'s last decade (or the last half of the lags)
  {
    std::vector<double> x, y;
    const double tmax = D.t.back();
    for (size_t k = 0; k < D.t.size(); ++k) if (D.t[k] >= tmax / 10) x.push_back(D.t[k]), y.push_back(D.g3[k]);
    if (x.size() >= 2) D.D = lin_slope(x, y) / 6;
  }
  for (size_t k = 1; k < D.t.size(); ++k)
    if (D.ree[k] <= std::exp(-1.0) && D.ree[k - 1] > std::exp(-1.0)) {
      const double w = (D.ree[k - 1] - std::exp(-1.0)) / (D.ree[k - 1] - D.ree[k]);
      D.tau_R = D.t[k - 1] + w * (D.t[k] - D.t[k - 1]);
      break;
    }
  // τ_e: g₁'s local exponent first below 3/8 after it has been near ½
  bool rouse = false;
  for (size_t k = 1; k + 1 < D.t.size(); ++k) {
    const double a = std::log(D.g1[k + 1] / D.g1[k - 1]) / std::log(D.t[k + 1] / D.t[k - 1]);
    if (a > 0.42 && a < 0.65) rouse = true;
    if (rouse && a < 0.375) { D.tau_e = D.t[k]; break; }
  }
  std::string n;
  if (D.tau_R == 0) n += "the end-to-end correlation does not fall to 1/e within the run (τ_R longer than it)";
  if (D.tau_e == 0) n += std::string(n.empty() ? "" : "; ") + "no reptation crossover in g₁ (unentangled, or the run is short)";
  D.note = n;
  return D;
}

double time_mapping(const std::vector<double>& ta, const std::vector<double>& ga, const std::vector<double>& tc, const std::vector<double>& gc, double* spread, int* points) {
  // t at a given g by log–log interpolation on a monotone part of each curve
  auto t_at = [](const std::vector<double>& t, const std::vector<double>& g, double v) {
    for (size_t k = 1; k < g.size(); ++k)
      if ((g[k - 1] <= v && g[k] >= v) && g[k] > g[k - 1] && t[k - 1] > 0) {
        const double w = std::log(v / g[k - 1]) / std::log(g[k] / g[k - 1]);
        return std::exp(std::log(t[k - 1]) + w * (std::log(t[k]) - std::log(t[k - 1])));
      }
    return kNaN;
  };
  double lo = 0, hi = 1e300;
  for (const auto* g : {&ga, &gc}) {
    double gmin = 1e300, gmax = 0;
    for (double v : *g) if (v > 0) gmin = std::min(gmin, v), gmax = std::max(gmax, v);
    lo = std::max(lo, gmin), hi = std::min(hi, gmax);
  }
  if (!(hi > lo)) throw std::invalid_argument("the two g₁ curves do not overlap in g₁");
  std::vector<double> r;
  for (int k = 0; k <= 40; ++k) {
    const double v = std::exp(std::log(lo) + (std::log(hi) - std::log(lo)) * k / 40.0);
    const double a = t_at(ta, ga, v), c = t_at(tc, gc, v);
    if (std::isfinite(a) && std::isfinite(c) && c > 0) r.push_back(std::log(a / c));
  }
  if (r.empty()) throw std::invalid_argument("no common g₁ values to map the times");
  double m = 0;
  for (double x : r) m += x;
  m /= double(r.size());
  double v = 0;
  for (double x : r) v += (x - m) * (x - m);
  if (spread) *spread = std::sqrt(v / double(r.size()));
  if (points) *points = int(r.size());
  return std::exp(m);
}

}  // namespace caps
