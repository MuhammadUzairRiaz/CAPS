// CAPS solvation: solvent and ions packed around a solute (see caps/solvate.hpp).
#include "caps/solvate.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>

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

// The solute's extent (unwrapped positions).
void extent(const System& s, Vec3& lo, Vec3& hi) {
  lo = {1e300, 1e300, 1e300}, hi = {-1e300, -1e300, -1e300};
  for (const auto& a : s.atoms)
    for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], a.pos[k]), hi[k] = std::max(hi[k], a.pos[k]);
}

}  // namespace

const std::vector<SolventInfo>& solvent_library() {
  static const std::vector<SolventInfo> lib = {
      {"water", "Water", "O", 0.997, "hydrated fillers, latex, biomolecules"},
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
    double r = 0.9572, theta = 104.52, qo = -1.1128;
    if (o.water_model == "TIP3P") qo = -0.834;
    else if (o.water_model == "SPC/E") r = 1.0, theta = 109.47, qo = -0.8476;
    else if (o.water_model != "TIP4P/2005") throw std::invalid_argument("water model TIP3P, SPC/E or TIP4P/2005, not '" + o.water_model + "'");
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
    if (name) *name = "water (" + o.water_model + ")";
    return w;
  }
  static std::map<std::string, System> cache;
  auto it = cache.find(sv->smiles);
  if (it == cache.end()) {
    BuildOptions bo;
    bo.forcefield = "uff";
    bo.seed = 3;
    System m = build_molecule(sv->smiles, bo).system;
    for (auto& a : m.atoms) a.resname = sv->id.substr(0, 3), a.charge = 0;
    for (auto& c : m.atoms[0].resname) c = char(std::toupper(static_cast<unsigned char>(c)));
    for (auto& a : m.atoms) a.resname = m.atoms[0].resname;
    m.has_mol = true;
    finish_types_by_element(m);
    it = cache.emplace(sv->smiles, m).first;
  }
  if (name) *name = sv->name;
  return it->second;
}

SolvatePlan solvate_plan(const System* solute, const SolvateOptions& o) {
  SolvatePlan P;
  const SolventInfo* sv = find_solvent(o.solvent);
  if (!sv) throw std::invalid_argument("unknown solvent '" + o.solvent + "'");
  const bool has_solute = solute && !solute->atoms.empty();
  Vec3 lo{0, 0, 0}, hi{0, 0, 0};
  if (has_solute) extent(*solute, lo, hi);
  const Vec3 span = has_solute ? hi - lo : Vec3{0, 0, 0};
  switch (o.shape) {
    case 0: P.box = {o.edge, o.edge, o.edge}; break;
    case 1: P.box = o.edges; break;
    default:
      if (!has_solute) throw std::invalid_argument("a box around the solute needs a solute");
      for (int k = 0; k < 3; ++k) P.box[k] = span[k] + 2 * o.padding;
  }
  for (int k = 0; k < 3; ++k)
    if (!(P.box[k] > 2 * o.tolerance)) throw std::invalid_argument("the box is too small for the tolerance");
  if (has_solute)
    for (int k = 0; k < 3; ++k)
      if (span[k] + o.tolerance > P.box[k])
        throw std::invalid_argument(std::string("the solute is ") + std::to_string(int(std::ceil(span[k]))) + " Å along " + "xyz"[k] +
                                    ": the box must be larger (or use the solute plus padding)");
  P.box_volume = P.box[0] * P.box[1] * P.box[2];

  // the solute's van der Waals volume on a 0.5 Å grid (periodic, the solute centred)
  if (has_solute) {
    P.solute_atoms = int(solute->atoms.size());
    const double h = 0.5;
    const int nx = std::max(1, int(std::ceil(P.box[0] / h))), ny = std::max(1, int(std::ceil(P.box[1] / h))), nz = std::max(1, int(std::ceil(P.box[2] / h)));
    std::vector<char> in(size_t(nx) * ny * nz, 0);
    const Vec3 shift = P.box * 0.5 - (lo + hi) * 0.5;
    for (const auto& a : solute->atoms) {
      P.solute_charge += a.charge;
      const double r = element(a.element).vdw;
      const Vec3 c = a.pos + shift;
      const int rx = int(std::ceil(r / h)) + 1;
      const int cx = int(std::floor(c[0] / h)), cy = int(std::floor(c[1] / h)), cz = int(std::floor(c[2] / h));
      for (int i = cx - rx; i <= cx + rx; ++i)
        for (int j = cy - rx; j <= cy + rx; ++j)
          for (int k = cz - rx; k <= cz + rx; ++k) {
            const Vec3 g{(i + 0.5) * h, (j + 0.5) * h, (k + 0.5) * h};
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
  const Region inside{Region::InsideBox, {0, 0, 0}, P.box};
  if (has_solute) {
    PackItem s;
    s.name = "solute";
    s.molecule = *solute;
    s.molecule.cell = Cell{};
    for (auto& a : s.molecule.atoms) a.mol = 1;
    s.fixed = true;
    s.center = true;
    s.position = P.box * 0.5;
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
