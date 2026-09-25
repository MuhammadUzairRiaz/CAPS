// χ from the energy of mixing by MD (see chimd.hpp).
#include "caps/chimd.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "caps/dynamics.hpp"
#include "caps/properties.hpp"
#include "caps/uff.hpp"
#include "caps/molecule.hpp"
#include "caps/pack.hpp"
#include "caps/relax.hpp"

namespace caps {
namespace {

constexpr double kR = 0.0019872043;   // kcal/mol/K
constexpr double kNAc = 0.602214076;  // g/cm³ ↔ (g/mol)/Å³

double total_mass(const System& s) {
  double m = 0;
  for (const auto& a : s.atoms) m += s.mass_of(a);
  return m;
}

// mean and standard error from block averages (5 blocks)
std::pair<double, double> block_mean(const std::vector<double>& v) {
  if (v.empty()) return {0, 0};
  double m = 0;
  for (double x : v) m += x;
  m /= double(v.size());
  const int nb = 5;
  if (int(v.size()) < 2 * nb) return {m, 0};
  const size_t per = v.size() / nb;
  std::vector<double> b(nb, 0.0);
  for (int k = 0; k < nb; ++k) {
    for (size_t i = 0; i < per; ++i) b[size_t(k)] += v[size_t(k) * per + i];
    b[size_t(k)] /= double(per);
  }
  double s2 = 0;
  for (double x : b) s2 += (x - m) * (x - m);
  return {m, std::sqrt(s2 / (nb - 1) / nb)};
}

System solvent_molecule_from(const std::string& smiles) {
  BuildOptions b;
  b.forcefield = "uff";
  return build_molecule(smiles, b).system;
}

}  // namespace

ChiMdResult chi_by_md(const ChiMdOptions& o) {
  ChiMdResult r;
  if (o.polymer.units.empty()) throw std::invalid_argument("χ by MD needs component A (a polymer)");
  if (!o.b_polymer && o.solvent_smiles.empty()) throw std::invalid_argument("χ by MD needs component B: a solvent SMILES or a second polymer");
  auto tell = [&](const std::string& stage, double f) {
    if (o.progress && !o.progress(stage, f)) throw ChiMdCancelled();
  };
  const System guest = o.b_polymer ? System{} : solvent_molecule_from(o.solvent_smiles);
  const double m_guest = o.b_polymer ? 0 : total_mass(guest);

  // relax with push-off, compress towards a liquid density, NPT; e = E/V averaged over the production rows
  auto measure = [&](System& s, ChiMdCell& cell, const std::string& stage) {
    tell(stage + ": relax", 0);
    RelaxOptions rl;
    rl.target_density = 0.85;
    rl.ftol = 2.0;
    rl.max_iterations = 2000;
    rl.energy.cutoff = o.cutoff;
    rl.energy.threads = o.threads;
    relax(s, rl);
    DynamicsOptions d;
    d.dt = o.dt;
    d.temperature = o.temperature;
    d.thermostat = Thermostat::Bussi;
    d.barostat = Barostat::CRescale;
    d.pressure = o.pressure;
    d.new_velocities = true;
    d.seed = o.seed;
    d.frame_every = 0;
    d.energy.cutoff = o.cutoff;
    d.energy.threads = o.threads;
    const int64_t eq = std::max<int64_t>(1, std::llround(o.eq_ps * 1000 / o.dt)), prod = std::max<int64_t>(1, std::llround(o.prod_ps * 1000 / o.dt));
    d.steps = eq + prod;
    d.thermo_every = int(std::max<int64_t>(1, prod / 200));
    d.frame_every = int(std::max<int64_t>(1, prod / 10));   // ten frames of the production run for the cohesive energy
    std::vector<double> E, V;
    Trajectory t;
    d.progress = [&](const ThermoRow& row) {
      if (row.step > eq && row.volume > 0) E.push_back(row.potential), V.push_back(row.volume);
      tell(stage + ": NPT", double(row.step) / double(d.steps));
      return true;
    };
    d.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
      if (step <= eq) return;
      std::vector<Vec3> p(x.size() / 3);
      for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
      t.positions.push_back(std::move(p));
      t.cells.push_back(c);
      t.timesteps.push_back(step);
    };
    run_dynamics(s, d);
    cell.atoms = int(s.atoms.size());
    int nm = 0;
    s.molecules(&nm);
    cell.molecules = nm;
    cell.energy = block_mean(E).first;
    cell.volume = block_mean(V).first;
    cell.density = cell.volume > 0 ? total_mass(s) / (kNAc * cell.volume) : 0;
    // cohesive energy density: (Σ E of each molecule alone − E of the cell) / V over the frames, the same force field
    tell(stage + ": cohesive energy", 0);
    t.topology = s;
    t.topology.unwrapped = false;   // the analysis makes every molecule whole (inserted ones may straddle the boundary)
    const ForceField ff = default_forcefield(s);
    AnalyzeOptions ao;
    ao.ff = &ff;
    ao.energy.cutoff = o.cutoff;
    ao.energy.threads = o.threads;
    ao.threads = o.threads;
    const auto props = analyze(t, {"ced"}, ao);
    const Property& ced = props.front();
    if (std::isnan(ced.value)) throw std::runtime_error("no cohesive energy for " + cell.name + (ced.notes.empty() ? "" : ": " + ced.notes.front()));
    constexpr double kJcm3 = 1.4393e-4;   // J/cm³ → kcal/mol/Å³
    cell.e_density = ced.value * kJcm3;
    cell.e_error = std::isnan(ced.error) ? 0 : ced.error * kJcm3;
  };

  auto grow_cell = [&](const ChainSpec& spec, int chains, double density, uint64_t seed) {
    GrowOptions g;
    g.chains = chains;
    g.density = density;
    g.seed = seed;
    g.auto_scale = true;
    return grow_chains(spec, g);
  };
  auto insert_into = [&](const System& host, int count, uint64_t seed) {
    PackOptions po;
    po.seed = seed;
    po.threads = o.threads;
    return insert_molecules(host, guest, count, po);
  };

  // A alone
  const int na = std::max(2, o.chains);
  r.a.name = "A";
  System A = grow_cell(o.polymer, na, 0.4, o.seed);
  const double mass_a = total_mass(A);
  measure(A, r.a, "A alone");

  // B alone
  int nb = 0;
  System B;
  r.b.name = "B";
  if (o.b_polymer) {
    nb = std::max(2, o.chains_b);
    B = grow_cell(o.polymer_b, nb, 0.4, o.seed + 1);
  } else {
    nb = o.solvent_molecules > 0 ? o.solvent_molecules : std::max(10, int(std::lround(mass_a / m_guest)));
    // the solvent packed into a periodic cube at half a liquid's density (NPT brings it to its own)
    const double L = std::cbrt(nb * m_guest / (0.5 * kNAc));
    PackOptions po;
    po.cell.a = {L, 0, 0}, po.cell.b = {0, L, 0}, po.cell.c = {0, 0, L};
    po.periodic = true;
    po.seed = o.seed + 1;
    po.threads = o.threads;
    PackItem it;
    it.name = "solvent";
    it.molecule = guest;
    it.count = nb;
    tell("B alone: pack", 0);
    B = pack({it}, po);
  }
  measure(B, r.b, "B alone");

  // the mixture: half of each
  const int ma = std::max(1, na / 2), mb = std::max(1, nb / 2);
  r.mix.name = "mixture";
  System M;
  if (o.b_polymer) {
    BlendOptions bo;
    bo.density = 0.4;
    bo.grow.seed = o.seed + 2;
    BlendComponent ca, cb;
    ca.spec = o.polymer, ca.chains = ma;
    cb.spec = o.polymer_b, cb.chains = mb;
    M = grow_blend({ca, cb}, bo);
  } else {
    // A's chains at a low density in a cell with room for B, then B into the free space
    const double m_mix = mass_a * double(ma) / na + m_guest * mb;
    GrowOptions g;
    g.chains = ma;
    g.box = std::cbrt(m_mix / (0.35 * kNAc));
    g.seed = o.seed + 2;
    g.auto_scale = true;
    System host = grow_chains(o.polymer, g);
    tell("mixture: pack", 0);
    M = insert_into(host, mb, o.seed + 3);
  }
  measure(M, r.mix, "mixture");

  // volume fractions from the pure cells' volumes per molecule; V_ref
  const double va = r.a.volume / na, vb = r.b.volume / nb;
  r.phi_a = ma * va / (ma * va + mb * vb);
  const double pb = 1 - r.phi_a;
  if (o.b_polymer) {
    const double ua = r.a.volume / (na * std::max(1, o.polymer.dp)), ub = r.b.volume / (nb * std::max(1, o.polymer_b.dp));
    r.v_ref = std::sqrt(ua * ub);
  } else {
    r.v_ref = vb;
  }
  const double RT = kR * o.temperature;
  // e_density is the cohesive energy density (positive): mixing costs energy when the mixture holds together less
  r.de_mix = r.phi_a * r.a.e_density + pb * r.b.e_density - r.mix.e_density;
  r.chi = r.v_ref * r.de_mix / (RT * r.phi_a * pb);
  const double err = std::sqrt(r.mix.e_error * r.mix.e_error + r.phi_a * r.phi_a * r.a.e_error * r.a.e_error + pb * pb * r.b.e_error * r.b.e_error);
  r.chi_error = r.v_ref * err / (RT * r.phi_a * pb);
  char b[320];
  std::snprintf(b, sizeof b, "χ = %.3f ± %.3f at %.0f K (φ_A %.3f, V_ref %.1f Å³, Δe_mix %.3g kcal/mol/Å³): enthalpic, from %.0f ps of NPT per cell",
                r.chi, r.chi_error, o.temperature, r.phi_a, r.v_ref, r.de_mix, o.prod_ps);
  r.notes.push_back(b);
  for (const ChiMdCell* c : {&r.a, &r.b, &r.mix}) {
    std::snprintf(b, sizeof b, "%s: %d atoms, %d molecules, ρ %.3f g/cm³, CED %.1f ± %.1f J/cm³", c->name.c_str(), c->atoms, c->molecules, c->density,
                  c->e_density / 1.4393e-4, c->e_error / 1.4393e-4);
    r.notes.push_back(b);
  }
  r.notes.push_back("force field: the built-in typing (GAFF for C and H, UFF otherwise), the same for all three cells; no entropic part");
  if (r.chi_error > 0.5 * std::fabs(r.chi)) r.notes.push_back("the error is large next to χ: run longer (prod_ps) or with more molecules");
  return r;
}

}  // namespace caps
