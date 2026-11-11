// CAPS solvation: solvent and ions packed around a solute (see caps/solvate.hpp).
#include "caps/solvate.hpp"
#include "caps/water.hpp"

#include <cctype>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <mutex>
#include <stdexcept>
#include <utility>

#include "caps/elements.hpp"
#include "caps/molecule.hpp"

namespace caps {

namespace {

constexpr double kAvogadro = 6.02214076e23;

void finish_types_by_element(System& s) {
  std::map<int, int> type_of;
  s.types.clear();
  for (auto& a : s.atoms) {
    auto it = type_of.find(a.element);
    if (it == type_of.end()) {
      it = type_of.emplace(a.element, int(type_of.size()) + 1).first;
      TypeInfo ti;
      ti.type = it->second;
      ti.mass = element(a.element).mass;
      ti.label = element(a.element).symbol;
      s.types.push_back(ti);
    }
    a.type = it->second;
  }
}

double molar_mass(const System& s) {
  double m = 0;
  for (const auto& a : s.atoms) m += element(a.element).mass;
  return m;
}

const SolventInfo* find_solvent(const std::string& id) {
  for (const auto& s : solvent_library()) if (s.id == id) return &s;
  return nullptr;
}

const SaltInfo* find_salt(const std::string& id) {
  for (const auto& s : salt_library()) if (s.id == id) return &s;
  return nullptr;
}

System ion(int z, int charge, const std::string& resname) {
  System s;
  Atom a;
  a.id = 1, a.mol = 1, a.element = z, a.charge = charge, a.name = element(z).symbol, a.resname = resname;
  s.atoms.push_back(a);
  s.has_mol = true;
  finish_types_by_element(s);
  return s;
}

// Where the solute's atoms go. A solute with a cell can be a network bonded through that cell (a sheet from caps
// surface, a crystal's chain): no unwrapping makes it whole, its image flags are arbitrary, and its bonds close only
// under the cell's period. Along such an axis the box keeps the solute's period and the atoms are folded into it; along
// the others the solute is made whole (a molecule written wrapped is joined again) and the solvent fills the space.
struct SoluteLayout {
  std::vector<Vec3> pos;                              // whole, folded into the period along a periodic axis
  std::array<bool, 3> periodic{false, false, false};  // x, y, z: the bonded network repeats along it
  Vec3 period{0, 0, 0};                               // Å, the repeat along a periodic axis
  Vec3 lo{0, 0, 0}, hi{0, 0, 0};                      // extent; along a periodic axis the period itself
};

SoluteLayout solute_layout(const System& s) {
  SoluteLayout L;
  const size_t n = s.atoms.size();
  L.pos.resize(n);
  for (size_t i = 0; i < n; ++i) L.pos[i] = s.atoms[i].pos;
  const Cell& c = s.cell;
  bool whole = true;   // every bond already at its nearest image: the positions are used as they are
  if (c.valid())
    for (const auto& b : s.bonds) {
      const Vec3 d = s.atoms[b.j].pos - s.atoms[b.i].pos;
      if (norm(d - c.minimum_image(d)) > 1e-6) { whole = false; break; }
    }
  if (!whole) {
    // walk the bonds at the nearest image (a finite molecule comes out whole)
    const auto nb = s.neighbours();
    std::vector<char> done(n, 0);
    std::vector<uint32_t> stack;
    for (uint32_t root = 0; root < n; ++root) {
      if (done[root]) continue;
      done[root] = 1;
      stack.push_back(root);
      while (!stack.empty()) {
        const uint32_t i = stack.back();
        stack.pop_back();
        for (uint32_t j : nb[i]) {
          if (done[j]) continue;
          L.pos[j] = L.pos[i] + c.minimum_image(L.pos[j] - L.pos[i]);
          done[j] = 1;
          stack.push_back(j);
        }
      }
    }
    // a bond the walk reached at another image than its nearest closes a ring around the cell: the network repeats
    // along the cell vectors of that lattice translation
    std::array<bool, 3> along{false, false, false};
    for (const auto& b : s.bonds) {
      const Vec3 walked = L.pos[b.j] - L.pos[b.i];
      const Vec3 f = c.to_fractional(c.origin + (walked - c.minimum_image(walked)));
      for (int k = 0; k < 3; ++k)
        if (std::fabs(f[k]) > 0.5) along[k] = true;
    }
    const Vec3 vec[3] = {c.a, c.b, c.c};
    std::string names;
    int count = 0;
    for (int k = 0; k < 3; ++k)
      if (along[k]) names += std::string(count++ ? " and " : "") + "abc"[k];
    if (count == 3)
      throw std::invalid_argument("the solute is bonded through its cell along a, b and c (a crystal or a network that fills its cell): "
                                  "there is no outside to solvate; pack the solvent into its pores with caps pack (the solute as a fixed structure)");
    for (int k = 0; k < 3; ++k) {
      if (!along[k]) continue;
      // the box is orthogonal, so the period must lie along its own axis (a along x, b along y, c along z)
      const double len = norm(vec[k]);
      for (int m = 0; m < 3; ++m)
        if (m != k && std::fabs(vec[k][m]) > 1e-4 * std::max(1.0, len))
          throw std::invalid_argument("the solute is bonded through its cell along " + names + " (a sheet or chain that repeats), and solvate keeps "
                                      "that period only for a cell vector along its own axis (a along x, b along y, c along z); this cell's " +
                                      std::string(1, "abc"[k]) + " is not along " + "xyz"[k] + ": solvate a finite piece or an orthogonal cell");
      L.periodic[k] = true;
      L.period[k] = len;
      for (auto& p : L.pos) p[k] -= len * std::floor((p[k] - c.origin[k]) / len);
    }
  }
  L.lo = {1e300, 1e300, 1e300}, L.hi = {-1e300, -1e300, -1e300};
  for (const auto& p : L.pos)
    for (int k = 0; k < 3; ++k) L.lo[k] = std::min(L.lo[k], p[k]), L.hi[k] = std::max(L.hi[k], p[k]);
  for (int k = 0; k < 3; ++k)
    if (L.periodic[k]) L.lo[k] = c.origin[k], L.hi[k] = c.origin[k] + L.period[k];
  return L;
}

}  // namespace

const std::vector<SolventInfo>& solvent_library() {
  static const std::vector<SolventInfo> lib = {
      {"water", "Water", "O", 0.997, "hydrated fillers, latex, polar polymers"},
      {"toluene", "Toluene", "Cc1ccccc1", 0.862, "swelling of NR, BR and SBR (Flory–Rehner crosslink density)"},
      {"cyclohexane", "Cyclohexane", "C1CCCCC1", 0.774, "swelling of non-polar rubbers"},
      {"hexane", "n-Hexane", "CCCCCC", 0.655, "extraction, swelling of EPDM"},
      {"chloroform", "Chloroform", "ClC(Cl)Cl", 1.479, "swelling, polymer solutions"},
      {"thf", "Tetrahydrofuran", "C1CCOC1", 0.883, "GPC solvent, polymer solutions"},
      {"acetone", "Acetone", "CC(C)=O", 0.784, "extraction of additives"},
      {"methanol", "Methanol", "CO", 0.787, "polar solvent"},
      {"ethanol", "Ethanol", "CCO", 0.785, "polar solvent"},
      {"dmso", "Dimethyl sulfoxide", "CS(C)=O", 1.095, "polar aprotic solvent"},
  };
  return lib;
}

const std::vector<SaltInfo>& salt_library() {
  static const std::vector<SaltInfo> lib = {
      {"NaCl", "Na", "Cl", 1, -1}, {"KCl", "K", "Cl", 1, -1}, {"LiCl", "Li", "Cl", 1, -1},
      {"CaCl2", "Ca", "Cl", 2, -1}, {"MgCl2", "Mg", "Cl", 2, -1}, {"ZnCl2", "Zn", "Cl", 2, -1},
  };
  return lib;
}

System solvent_molecule(const SolvateOptions& o, std::string* name) {
  const SolventInfo* sv = find_solvent(o.solvent);
  if (!sv) throw std::invalid_argument("unknown solvent '" + o.solvent + "'");
  if (sv->id == "water") {
    // the model's geometry and charges (caps/water.hpp); a four- or five-site model's negative charge sits on O until
    // the model is applied with its extra sites (the Field's water model, or the edit water_model)
    const WaterModel& wm = water_model(o.water_model);
    const double r = wm.r_oh, theta = wm.theta, qo = wm.sites == 5 ? 2 * wm.q_neg : wm.q_neg;   // TIP5P: both lone pairs' charge
    const double h = theta * 3.14159265358979323846 / 360.0;
    System w;
    const Vec3 p[3] = {{0, 0, 0}, {r * std::sin(h), r * std::cos(h), 0}, {-r * std::sin(h), r * std::cos(h), 0}};
    const char* names[3] = {"OW", "HW1", "HW2"};
    for (int k = 0; k < 3; ++k) {
      Atom a;
      a.id = k + 1, a.mol = 1, a.element = k ? 1 : 8, a.charge = k ? -qo / 2 : qo, a.name = names[k], a.resname = "SOL", a.pos = p[k];
      w.atoms.push_back(a);
    }
    w.bonds = {{0, 1, 1}, {0, 2, 1}};
    w.has_mol = true;
    finish_types_by_element(w);
    if (name) *name = "water (" + wm.name + ")";
    return w;
  }
  // shared by every caller (the Studio plans on the UI thread while a preview solvates on a worker), so the map is
  // only touched under the lock; the molecule is built outside it and the result is copied out while still held
  static std::map<std::string, System> cache;
  static std::mutex cache_mutex;
  {
    std::lock_guard<std::mutex> lock(cache_mutex);
    auto it = cache.find(sv->smiles);
    if (it != cache.end()) {
      if (name) *name = sv->name;
      return it->second;
    }
  }
  System built;
  {
    BuildOptions bo;
    bo.forcefield = "uff";
    bo.seed = 3;
    System m = build_molecule(sv->smiles, bo).system;
    for (auto& a : m.atoms) a.resname = sv->id.substr(0, 3), a.charge = 0;
    for (auto& c : m.atoms[0].resname) c = char(std::toupper(static_cast<unsigned char>(c)));
    for (auto& a : m.atoms) a.resname = m.atoms[0].resname;
    m.has_mol = true;
    finish_types_by_element(m);
    built = std::move(m);
  }
  std::lock_guard<std::mutex> lock(cache_mutex);
  auto it = cache.emplace(sv->smiles, std::move(built)).first;   // keeps the first entry if another thread got there first
  if (name) *name = sv->name;
  return it->second;
}

SolvatePlan solvate_plan(const System* solute, const SolvateOptions& o) {
  SolvatePlan P;
  const SolventInfo* sv = find_solvent(o.solvent);
  if (!sv) throw std::invalid_argument("unknown solvent '" + o.solvent + "'");
  const bool has_solute = solute && !solute->atoms.empty();
  SoluteLayout lay;
  if (has_solute) lay = solute_layout(*solute);
  const Vec3 lo = lay.lo, hi = lay.hi;
  const Vec3 span = has_solute ? hi - lo : Vec3{0, 0, 0};
  switch (o.shape) {
    case 0: P.box = {o.edge, o.edge, o.edge}; break;
    case 1: P.box = o.edges; break;
    default:
      if (!has_solute) throw std::invalid_argument("a box around the solute needs a solute");
      for (int k = 0; k < 3; ++k) P.box[k] = span[k] + 2 * o.padding;
  }
  // a solute bonded through its cell keeps its period: the edge and the padding apply to the other axes only, and a
  // box given in full must repeat it
  std::string per_axes;
  for (int k = 0; k < 3; ++k) {
    if (!lay.periodic[k]) continue;
    if (o.shape == 1 && std::fabs(o.edges[k] - lay.period[k]) > 0.01) {
      char b[200];
      std::snprintf(b, sizeof b, "the solute repeats every %.4f Å along %c (bonded through its cell): the box must keep that length along %c", lay.period[k],
                    "xyz"[k], "xyz"[k]);
      throw std::invalid_argument(b);
    }
    P.box[k] = lay.period[k];
    per_axes += std::string(per_axes.empty() ? "" : " and ") + "xyz"[k];
  }
  if (!per_axes.empty()) {
    std::string lens;
    for (int k = 0; k < 3; ++k)
      if (lay.periodic[k]) {
        char b[32];
        std::snprintf(b, sizeof b, "%s%.3f", lens.empty() ? "" : " × ", lay.period[k]);
        lens += b;
      }
    P.notes.push_back("the solute is bonded through its cell along " + per_axes + ": the box keeps its period there (" + lens +
                      " Å) and the solvent fills the other axes");
  }
  for (int k = 0; k < 3; ++k)
    if (!(P.box[k] > 2 * o.tolerance)) throw std::invalid_argument("the box is too small for the tolerance");
  if (has_solute)
    for (int k = 0; k < 3; ++k)
      if (!lay.periodic[k] && span[k] + o.tolerance > P.box[k])
        throw std::invalid_argument(std::string("the solute is ") + std::to_string(int(std::ceil(span[k]))) + " Å along " + "xyz"[k] +
                                    ": the box must be larger (or use the solute plus padding)");
  P.box_volume = P.box[0] * P.box[1] * P.box[2];

  // the solute's van der Waals volume on a 0.5 Å grid (periodic, the solute centred)
  if (has_solute) {
    P.solute_atoms = int(solute->atoms.size());
    // about 0.5 Å, each axis divided evenly so the grid repeats with the box
    const int nx = std::max(1, int(std::ceil(P.box[0] / 0.5))), ny = std::max(1, int(std::ceil(P.box[1] / 0.5))), nz = std::max(1, int(std::ceil(P.box[2] / 0.5)));
    const double hx = P.box[0] / nx, hy = P.box[1] / ny, hz = P.box[2] / nz, h = std::min({hx, hy, hz});
    std::vector<char> in(size_t(nx) * ny * nz, 0);
    const Vec3 shift = P.box * 0.5 - (lo + hi) * 0.5;
    for (size_t ai = 0; ai < solute->atoms.size(); ++ai) {
      const auto& a = solute->atoms[ai];
      P.solute_charge += a.charge;
      const double r = element(a.element).vdw;
      const Vec3 c = lay.pos[ai] + shift;
      const int rx = int(std::ceil(r / h)) + 1;
      const int cx = int(std::floor(c[0] / hx)), cy = int(std::floor(c[1] / hy)), cz = int(std::floor(c[2] / hz));
      for (int i = cx - rx; i <= cx + rx; ++i)
        for (int j = cy - rx; j <= cy + rx; ++j)
          for (int k = cz - rx; k <= cz + rx; ++k) {
            const Vec3 g{(i + 0.5) * hx, (j + 0.5) * hy, (k + 0.5) * hz};
            if (norm(g - c) > r) continue;
            const int ii = ((i % nx) + nx) % nx, jj = ((j % ny) + ny) % ny, kk = ((k % nz) + nz) % nz;
            in[(size_t(ii) * ny + jj) * nz + kk] = 1;
          }
    }
    P.solute_volume = double(std::count(in.begin(), in.end(), 1)) * P.box_volume / double(in.size());
  }
  P.free_volume = std::max(0.0, P.box_volume - P.solute_volume);

  const System one = solvent_molecule(o, &P.solvent_name);
  P.solvent_mass = molar_mass(one);
  const double rho = o.density > 0 ? o.density : sv->density;
  P.solvent = o.molecules > 0 ? o.molecules : int(std::lround(rho * P.free_volume * 1e-24 * kAvogadro / P.solvent_mass));

  // ions
  const SaltInfo* salt = find_salt(o.salt);
  if (o.ion_mode != 0 && !salt) throw std::invalid_argument("unknown salt '" + o.salt + "'");
  if (salt) P.cation = salt->cation, P.anion = salt->anion;
  if (o.ion_mode == 3) {
    P.cations = std::max(0, o.cations), P.anions = std::max(0, o.anions);
  } else if (o.ion_mode == 1 || o.ion_mode == 2) {
    int pairs = 0;
    if (o.ion_mode == 2) pairs = int(std::lround(o.concentration * kAvogadro * P.free_volume * 1e-27));
    P.cations = pairs;
    P.anions = pairs * (salt->zc / -salt->za);
    const long q = std::lround(P.solute_charge);
    if (std::fabs(P.solute_charge - double(q)) > 0.05)
      P.notes.push_back("the solute's charge is " + std::to_string(P.solute_charge).substr(0, 6) + " e, not a whole number: neutralised to the nearest integer");
    if (q > 0) P.anions += int((q + (-salt->za) - 1) / (-salt->za));
    if (q < 0) P.cations += int((-q + salt->zc - 1) / salt->zc);
    const long left = q + long(P.cations) * salt->zc + long(P.anions) * salt->za;
    if (left != 0) P.notes.push_back("the box keeps a net charge of " + std::to_string(left) + " e (the ions' charges do not divide the solute's)");
  }
  if (P.cations + P.anions > 0 && sv->id == "water") P.solvent = std::max(0, P.solvent - P.cations - P.anions);   // each ion displaces a water
  if (P.solvent == 0 && P.cations + P.anions == 0 && !has_solute) throw std::invalid_argument("nothing to pack");
  const int units = salt ? std::min(P.cations, P.anions / std::max(1, salt->zc / -salt->za)) : 0;   // formula units of the salt
  P.concentration = P.free_volume > 0 ? double(units) / (kAvogadro * P.free_volume * 1e-27) : 0;
  double mass = P.solvent * P.solvent_mass;
  if (has_solute) mass += molar_mass(*solute);
  if (salt) mass += P.cations * element(element_from_symbol(salt->cation)).mass + P.anions * element(element_from_symbol(salt->anion)).mass;
  P.density = mass / kAvogadro / (P.box_volume * 1e-24);
  return P;
}

System solvate(const System* solute, const SolvateOptions& o, SolvateReport* report) {
  SolvateReport R;
  R.plan = solvate_plan(solute, o);
  const SolvatePlan& P = R.plan;
  const bool has_solute = solute && !solute->atoms.empty();
  std::vector<PackItem> items;
  const Region inside{Region::InsideBox, {0, 0, 0}, P.box, 0, 0, {}};
  if (has_solute) {
    // placed as the plan measured it: whole (folded into the period along a periodic axis), its extent centred in the
    // box; bonds through the old cell stay at their nearest image in the box, which keeps that period
    const SoluteLayout lay = solute_layout(*solute);
    const Vec3 shift = P.box * 0.5 - (lay.lo + lay.hi) * 0.5;
    PackItem s;
    s.name = "solute";
    s.molecule = *solute;
    s.molecule.cell = Cell{};
    for (size_t i = 0; i < s.molecule.atoms.size(); ++i) {
      s.molecule.atoms[i].pos = lay.pos[i] + shift;
      s.molecule.atoms[i].mol = 1;
    }
    s.fixed = true;
    items.push_back(std::move(s));
  }
  if (P.solvent > 0) {
    PackItem w;
    w.name = P.solvent_name;
    w.molecule = solvent_molecule(o);
    w.count = P.solvent;
    w.regions = {inside};
    items.push_back(std::move(w));
  }
  if (P.cations > 0 || P.anions > 0) {
    const SaltInfo* salt = find_salt(o.salt);
    if (P.cations > 0) {
      PackItem c;
      c.name = salt->cation;
      c.molecule = ion(element_from_symbol(salt->cation), salt->zc, salt->cation);
      c.count = P.cations;
      c.regions = {inside};
      items.push_back(std::move(c));
    }
    if (P.anions > 0) {
      PackItem a;
      a.name = salt->anion;
      a.molecule = ion(element_from_symbol(salt->anion), salt->za, salt->anion);
      a.count = P.anions;
      a.regions = {inside};
      items.push_back(std::move(a));
    }
  }
  PackOptions po;
  po.tolerance = o.tolerance;
  po.cell.a = {P.box[0], 0, 0}, po.cell.b = {0, P.box[1], 0}, po.cell.c = {0, 0, P.box[2]};
  po.periodic = true;
  po.seed = o.seed;
  po.threads = o.threads;
  po.progress = o.progress;
  System s = pack(items, po, &R.pack);
  s.title = (has_solute ? (solute->title.empty() ? std::string("solute") : solute->title) + " in " : std::string()) + P.solvent_name;
  char b[240];
  std::snprintf(b, sizeof b, "%.1f × %.1f × %.1f Å · %d %s · %d %s · %d %s · %.3f g/cm³", P.box[0], P.box[1], P.box[2], P.solvent, P.solvent_name.c_str(), P.cations,
                P.cation.empty() ? "cations" : (P.cation + "+").c_str(), P.anions, P.anion.empty() ? "anions" : (P.anion + "−").c_str(), P.density);
  s.notes.insert(s.notes.begin(), b);
  for (const auto& n : P.notes) s.notes.push_back(n);
  R.notes = P.notes;
  if (report) *report = R;
  return s;
}

}  // namespace caps
