// CAPS coarse-grained polymer builder (see caps/cg_build.hpp).
#include "caps/cg_build.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <set>
#include <stdexcept>

#include "caps/polymer.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kAvogadro = 6.02214076e23;

double unit01(std::mt19937_64& g) { return double(g() >> 11) * (1.0 / 9007199254740992.0); }

// draws from a bonded table's distribution with its Jacobian put back (r², sin θ; dihedrals as they are)
struct Drawer {
  std::vector<double> x, cdf;
  double bin = 0;
  explicit Drawer(const CgBondedTable& T) {
    const size_t n = T.xh.size();
    if (n < 2) throw std::invalid_argument("the " + T.key + " table has no distribution");
    bin = T.xh[1] - T.xh[0];
    double c = 0;
    for (size_t i = 0; i < n; ++i) {
      const double jac = T.kind == 0 ? T.xh[i] * T.xh[i] : T.kind == 1 ? std::sin(T.xh[i] * kPi / 180.0) : 1.0;
      c += std::max(0.0, T.P[i]) * jac;
      x.push_back(T.xh[i]);
      cdf.push_back(c);
    }
    if (c <= 0) throw std::invalid_argument("the " + T.key + " table has an empty distribution");
    for (auto& v : cdf) v /= c;
  }
  double operator()(std::mt19937_64& g) const {
    const double u = unit01(g);
    const size_t i = std::min(x.size() - 1, size_t(std::lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin()));
    return x[i] + (unit01(g) - 0.5) * bin;
  }
};

Vec3 unit_vector(std::mt19937_64& g) {
  const double z = 2 * unit01(g) - 1, phi = 2 * kPi * unit01(g), s = std::sqrt(std::max(0.0, 1 - z * z));
  return {s * std::cos(phi), s * std::sin(phi), z};
}

// the next bead from the last three (or fewer), bond r, angle θ (degrees) at the last bead, dihedral φ (IUPAC)
Vec3 place(const Vec3* a, const Vec3& b, const Vec3& c, double r, double theta, double phi, std::mt19937_64& g) {
  const double t = theta * kPi / 180.0;
  Vec3 bc = c - b;
  bc = bc * (1.0 / norm(bc));
  Vec3 n;
  if (a) n = cross(b - *a, bc);
  else {   // no bead before b: any direction perpendicular to bc
    Vec3 q = unit_vector(g);
    n = cross(q, bc);
  }
  const double nn = norm(n);
  if (nn < 1e-9) { Vec3 q = std::fabs(bc[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}; n = cross(q, bc); }
  n = n * (1.0 / norm(n));
  const Vec3 m = cross(n, bc);
  const double p = phi * kPi / 180.0;
  return c + bc * (-r * std::cos(t)) + m * (r * std::sin(t) * std::cos(p)) + n * (r * std::sin(t) * std::sin(p));
}

const CgBondedTable* table(const std::vector<CgBondedTable>& v, const std::string& key) {
  for (const auto& t : v) if (t.key == key && t.sampled && !t.P.empty()) return &t;
  return nullptr;
}
}  // namespace

CgBuildResult build_cg_polymer(const CgBuildOptions& o, const CgBondedResult& bonded, const std::map<std::string, double>& masses) {
  if (o.units.empty()) throw std::invalid_argument("give the repeat units as bead sequences (BS = B S …)");
  if (o.dp < 1 || o.chains < 1) throw std::invalid_argument("dp and chains must be at least 1");
  if (o.density <= 0) throw std::invalid_argument("give the melt's density (g/cm³)");
  for (const auto& u : o.units)
    for (const auto& b : u.beads)
      if (!masses.count(b)) throw std::invalid_argument("no mass for bead " + b + " (give the maps, or the masses)");
  CgBuildResult R;
  std::mt19937_64 g(o.seed);
  // composition (units in order; missing ones 0)
  std::vector<double> share;
  double tot = 0;
  for (const auto& u : o.units) share.push_back(o.composition.count(u.name) ? o.composition.at(u.name) : (o.composition.empty() ? 1.0 : 0.0)), tot += share.back();
  if (tot <= 0) throw std::invalid_argument("the composition gives no unit a share");
  for (auto& s : share) s /= tot;
  auto pick = [&](const std::vector<double>& w) {
    double u = unit01(g), c = 0, t = 0;
    for (double x : w) t += x;
    for (size_t k = 0; k < w.size(); ++k) { c += w[k] / t; if (u < c) return int(k); }
    return int(w.size()) - 1;
  };
  auto index_of = [&](const std::string& n) {
    for (size_t k = 0; k < o.units.size(); ++k) if (o.units[k].name == n) return int(k);
    throw std::invalid_argument("no unit " + n);
  };
  // chain lengths
  std::vector<int> dps = o.lengths == "monodisperse" || o.lengths.empty() ? std::vector<int>(size_t(o.chains), o.dp)
                                                                          : draw_chain_lengths(o.lengths, double(o.dp), o.pdi, o.chains, o.seed);
  R.chain_dp = dps;
  // sequences
  std::vector<std::vector<int>> seq;
  for (int c = 0; c < o.chains; ++c) {
    const int n = dps[size_t(c)];
    std::vector<int> s;
    if (o.sequence == "bernoulli") for (int i = 0; i < n; ++i) s.push_back(pick(share));
    else if (o.sequence == "markov") {
      s.push_back(pick(share));
      for (int i = 1; i < n; ++i) {
        const auto row = o.markov.find(o.units[size_t(s.back())].name);
        if (row == o.markov.end()) { s.push_back(pick(share)); continue; }
        std::vector<double> w(o.units.size(), 0.0);
        for (const auto& [to, p] : row->second) w[size_t(index_of(to))] = p;
        s.push_back(pick(w));
      }
    } else if (o.sequence == "block") {
      if (o.blocks.empty()) throw std::invalid_argument("block sequences need the block lengths");
      for (size_t b = 0; int(s.size()) < n; ++b)
        for (int k = 0; k < o.blocks[b % o.blocks.size()] && int(s.size()) < n; ++k) s.push_back(int(b % o.units.size()));
    } else if (o.sequence == "gradient") {
      if (o.units.size() < 2) throw std::invalid_argument("a gradient needs two units or more");
      for (int i = 0; i < n; ++i) {
        const double p0 = n > 1 ? 1.0 - double(i) / double(n - 1) : 0.5;
        std::vector<double> w(o.units.size(), 0.0);
        double rest = 0;
        for (size_t k = 1; k < o.units.size(); ++k) rest += share[k];
        for (size_t k = 1; k < o.units.size(); ++k) w[k] = (1 - p0) * (rest > 0 ? share[k] / rest : 1.0 / double(o.units.size() - 1));
        w[0] = p0;
        s.push_back(pick(w));
      }
    } else if (o.sequence == "alternating") for (int i = 0; i < n; ++i) s.push_back(i % int(o.units.size()));
    else if (o.sequence == "pattern") {
      if (o.pattern.empty()) throw std::invalid_argument("give the pattern (letters A, B, … for the units in order)");
      for (int i = 0; i < n; ++i) {
        const int k = o.pattern[size_t(i) % o.pattern.size()] - 'A';
        if (k < 0 || k >= int(o.units.size())) throw std::invalid_argument("the pattern names a unit the list does not have");
        s.push_back(k);
      }
    } else throw std::invalid_argument("sequence: bernoulli, markov, block, gradient, alternating or pattern");
    for (int k : s) R.units_drawn[o.units[size_t(k)].name]++;
    seq.push_back(std::move(s));
  }
  // beads, the box
  double mass = 0;
  std::vector<std::vector<std::string>> kinds;
  for (const auto& s : seq) {
    std::vector<std::string> k;
    for (int u : s) for (const auto& b : o.units[size_t(u)].beads) k.push_back(b), mass += masses.at(b);
    kinds.push_back(std::move(k));
  }
  const double volume = mass / kAvogadro / o.density * 1e24;   // Å³
  R.box = std::cbrt(volume);
  // samplers per bonded type
  std::map<std::string, Drawer> bond_d, angle_d, dih_d;
  std::set<std::string> uniform_dih;
  auto drawer = [&](std::map<std::string, Drawer>& cache, const std::vector<CgBondedTable>& v, const std::string& key, bool required, const char* what) -> const Drawer* {
    auto it = cache.find(key);
    if (it != cache.end()) return &it->second;
    const CgBondedTable* T = table(v, key);
    if (!T) {
      if (required) throw std::invalid_argument(std::string("the bonded set has no ") + what + " " + key);
      return nullptr;
    }
    return &cache.emplace(key, Drawer(*T)).first->second;
  };
  // walks
  System& S = R.beads;
  S.title = "coarse-grained melt";
  S.cell.a = {R.box, 0, 0}, S.cell.b = {0, R.box, 0}, S.cell.c = {0, 0, R.box};
  S.has_mol = true;
  CgMapping& M = R.topology;
  M.position = "built";
  double contour = 0, ree2 = 0;
  for (size_t c = 0; c < kinds.size(); ++c) {
    const auto& k = kinds[c];
    std::vector<Vec3> p;
    p.push_back({unit01(g) * R.box, unit01(g) * R.box, unit01(g) * R.box});
    double L = 0;
    for (size_t i = 1; i < k.size(); ++i) {
      const double r = (*drawer(bond_d, bonded.bonds, cg_key({k[i - 1], k[i]}), true, "bond"))(g);
      L += r;
      if (i == 1) { p.push_back(p[0] + unit_vector(g) * r); continue; }
      const double th = (*drawer(angle_d, bonded.angles, cg_key({k[i - 2], k[i - 1], k[i]}), true, "angle"))(g);
      double ph = 360 * unit01(g) - 180;
      if (i >= 3) {
        const std::string key = cg_key({k[i - 3], k[i - 2], k[i - 1], k[i]});
        if (const Drawer* d = drawer(dih_d, bonded.dihedrals, key, false, "dihedral")) {
          ph = (*d)(g);
          // the table's sign is for the canonical key's direction: a reversed key reads the same dihedral (φ(abcd) = φ(dcba))
        } else uniform_dih.insert(key);
      }
      p.push_back(place(i >= 3 ? &p[i - 3] : nullptr, p[i - 2], p[i - 1], r, th, ph, g));
    }
    contour += L;
    const double re = norm(p.back() - p.front());
    R.ree.push_back(re);
    ree2 += re * re;
    std::vector<int> chain;
    for (size_t i = 0; i < k.size(); ++i) {
      Atom a;
      a.id = int64_t(S.atoms.size() + 1);
      a.mol = int64_t(c + 1);
      a.name = k[i];
      a.element = 0;
      a.pos = p[i];
      chain.push_back(int(S.atoms.size()));
      M.bead_kind.push_back(k[i]);
      M.bead_mol.push_back(int(c));
      M.bead_mass.push_back(masses.at(k[i]));
      S.atoms.push_back(a);
      if (i > 0) M.bonds.push_back({chain[i - 1], chain[i]}), S.bonds.push_back({uint32_t(chain[i - 1]), uint32_t(chain[i]), 1});
      if (i > 1) M.angles.push_back({chain[i - 2], chain[i - 1], chain[i]});
      if (i > 2) M.dihedrals.push_back({chain[i - 3], chain[i - 2], chain[i - 1], chain[i]});
    }
    M.chains.push_back(chain);
  }
  R.ree_rms = std::sqrt(ree2 / double(kinds.size()));
  R.contour = contour / double(kinds.size());
  // types (sorted kinds), masses
  std::set<std::string> ks(M.bead_kind.begin(), M.bead_kind.end());
  int ti = 0;
  for (const auto& k : ks) S.types.push_back({++ti, masses.at(k), k});
  for (auto& a : S.atoms) a.type = int(std::distance(ks.begin(), ks.find(a.name))) + 1;
  // internal distances of the walks
  CgTopology top;
  top.kind = M.bead_kind, top.mol = M.bead_mol, top.bonds = M.bonds, top.chains = M.chains;
  std::vector<Vec3> pos;
  for (const auto& a : S.atoms) pos.push_back(a.pos);
  R.internal = cg_internal_distances(top, {pos}, {Cell{}});
  for (size_t c = 0; c < seq.size() && c < 3; ++c) {
    std::string q;
    for (int u : seq[c]) q += (q.empty() ? "" : " ") + o.units[size_t(u)].name;
    R.sequences.push_back(q);
  }
  char b[320];
  std::snprintf(b, sizeof b, "%d chains, %zu beads, box %.1f Å at %.3f g/cm³ · ⟨R_ee²⟩^½ %.1f Å, contour %.0f Å · L / R_ee %.2f", o.chains, S.atoms.size(), R.box, o.density,
                R.ree_rms, R.contour, R.box / R.ree_rms);
  R.notes.push_back(b);
  if (R.box < R.ree_rms)
    R.notes.push_back("WARNING: the box is smaller than the chains (L < ⟨R_ee²⟩^½): chains meet their own periodic images — use more chains of the same length (a larger box) so that L ≥ R_ee, better 1.5 R_ee");
  if (!uniform_dih.empty()) {
    std::string u;
    for (const auto& k : uniform_dih) u += (u.empty() ? "" : ", ") + k;
    R.notes.push_back("no dihedral table for " + u + ": those dihedrals drawn uniformly");
  }
  R.notes.push_back("internal coordinates are drawn independently (no correlation between neighbouring angles or dihedrals) and the walks overlap: push them apart and anneal (in.cg_equil), then check ⟨R²(n)⟩/n against the all-atom model");
  return R;
}

std::vector<double> cg_internal_distances(const CgTopology& t, const std::vector<std::vector<Vec3>>& frames, const std::vector<Cell>& cells, std::vector<double>* counts) {
  size_t nmax = 0;
  for (const auto& c : t.chains) nmax = std::max(nmax, c.size());
  std::vector<double> sum(nmax, 0.0), cnt(nmax, 0.0);
  for (size_t f = 0; f < frames.size(); ++f) {
    const Cell& cell = f < cells.size() ? cells[f] : Cell{};
    const bool pbc = cell.valid();
    for (const auto& ch : t.chains) {
      // whole along the chain
      std::vector<Vec3> p(ch.size());
      p[0] = frames[f][size_t(ch[0])];
      for (size_t i = 1; i < ch.size(); ++i) {
        const Vec3 d = frames[f][size_t(ch[i])] - frames[f][size_t(ch[i - 1])];
        p[i] = p[i - 1] + (pbc ? cell.minimum_image(d) : d);
      }
      for (size_t i = 0; i < p.size(); ++i)
        for (size_t j = i + 1; j < p.size(); ++j) {
          const Vec3 d = p[j] - p[i];
          sum[j - i] += dot(d, d), cnt[j - i] += 1;
        }
    }
  }
  std::vector<double> out;
  for (size_t n = 1; n < nmax; ++n) out.push_back(cnt[n] > 0 ? sum[n] / cnt[n] / double(n) : 0.0);
  if (counts) counts->assign(cnt.begin() + (nmax > 0 ? 1 : 0), cnt.end());
  return out;
}

std::vector<std::string> write_cg_build(const CgBuildResult& r, const CgTypes& types, const std::string& stem) {
  std::vector<std::string> files;
  write_cg_lammps_data(r.topology, r.beads, types, stem + ".cg.data");
  files.push_back(stem + ".cg.data");
  Json j = Json::object();
  j["format"] = "caps-cg-map";
  j["version"] = 1;
  j["position"] = "built";
  j["atoms"] = 0;
  Json beads = Json::array();
  for (size_t b = 0; b < r.topology.bead_kind.size(); ++b) {
    Json e = Json::object();
    e["kind"] = r.topology.bead_kind[b];
    e["mol"] = r.topology.bead_mol[b];
    e["mass"] = r.topology.bead_mass[b];
    beads.push_back(e);
  }
  j["beads"] = beads;
  Json chains = Json::array();
  for (const auto& c : r.topology.chains) {
    Json x = Json::array();
    for (int b : c) x.push_back(b);
    chains.push_back(x);
  }
  j["chains"] = chains;
  Json bonds = Json::array();
  for (const auto& [x, y] : r.topology.bonds) {
    Json e = Json::array();
    e.push_back(x);
    e.push_back(y);
    bonds.push_back(e);
  }
  j["bonds"] = bonds;
  std::ofstream f(stem + ".map.json");
  f << j.dump(0) << "\n";
  files.push_back(stem + ".map.json");
  return files;
}

std::string cg_equil_deck(const CgEquilOptions& o) {
  const char* special = o.exclude <= 2 ? "0 1 1" : o.exclude == 3 ? "0 0 1" : "0 0 0";
  char b[4096];
  std::snprintf(b, sizeof b,
                "# CAPS equilibration of a built coarse-grained melt (LAMMPS; dihedral_style table/cut needs EXTRA-MOLECULE)\n"
                "#   lmp -in in.cg_equil -var DATA melt.cg.data -var BONDED ../bonded -var PAIR ../lj -var OUT equil\n"
                "# 1 push-off: the walks overlap; a soft repulsion ramped up pushes them apart while the bonded terms keep each chain's\n"
                "#   local structure (Auhl, Everaers, Grest, Kremer & Plimpton, J. Chem. Phys. 119, 12718 (2003)); steps limited in length\n"
                "# 2 the model's pairs, a minimisation, a hot anneal (chains relax faster), cooling, NPT\n"
                "variable T index %g\nvariable P index %g\nvariable TA index %g\nvariable DT index %g\nvariable SOFT index %g\nvariable SEED index 4928459\n"
                "variable PUSH index %lld\nvariable RELAX index %lld\nvariable ANNEAL index %lld\nvariable COOL index %lld\nvariable NPT index %lld\nvariable DUMP index %lld\n"
                "units real\natom_style full\nboundary p p p\nspecial_bonds lj %s\nread_data ${DATA}\ninclude ${BONDED}/bonded.in\n"
                "neighbor 3.0 bin\nneigh_modify delay 0 every 1 check yes\ntimestep ${DT}\n"
                "thermo_style custom step temp press density pe ebond eangle edihed epair\nthermo 5000\n"
                "log ${OUT}.log\n"
                "# --- 1 push-off\n"
                "pair_style soft 5.0\npair_coeff * * 0.0\nvariable A equal ramp(0.0,${SOFT})\nfix push all adapt 1 pair soft a * * v_A\n"
                "velocity all create ${TA} ${SEED} mom yes rot yes dist gaussian\n"
                "fix lim all nve/limit 0.1\nfix lang all langevin ${TA} ${TA} $(100.0*dt) ${SEED}\n"
                "dump push all custom ${DUMP} ${OUT}.pushoff.lammpstrj id mol type xu yu zu\ndump_modify push sort id\n"
                "run ${PUSH}\nunfix push\nunfix lim\nunfix lang\nundump push\n"
                "# --- 2 the model's pairs: first with steps still limited (overlaps the soft pairs left), then free\ninclude ${PAIR}/pair.in\n"
                "fix lim all nve/limit 0.05\nfix lang all langevin ${TA} ${TA} $(100.0*dt) ${SEED}\nrun ${RELAX}\nunfix lim\nunfix lang\nreset_timestep 0\n"
                "fix md all nvt temp ${TA} ${TA} $(100.0*dt)\n"
                "dump eq all custom ${DUMP} ${OUT}.lammpstrj id mol type xu yu zu\ndump_modify eq sort id\n"
                "run ${ANNEAL}\nunfix md\n"
                "fix md all nvt temp ${TA} ${T} $(100.0*dt)\nrun ${COOL}\nunfix md\n"
                "fix md all npt temp ${T} ${T} $(100.0*dt) iso ${P} ${P} $(1000.0*dt)\nrun ${NPT}\n"
                "write_data ${OUT}.data\n"
                "# Optional, homopolymers with one bond type: double bridging (fix bond/swap) relaxes long chains faster; LAMMPS's\n"
                "# restrictions for fix bond/swap apply (see its documentation) — not used here.\n",
                o.T, o.P, o.anneal_T, o.dt, o.soft_max, static_cast<long long>(o.pushoff_steps), static_cast<long long>(o.pushoff_steps / 3), static_cast<long long>(o.anneal_steps),
                static_cast<long long>(o.cool_steps), static_cast<long long>(o.npt_steps), static_cast<long long>(o.dump_every), special);
  return b;
}

}  // namespace caps
