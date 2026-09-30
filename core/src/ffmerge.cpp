// CAPS force fields by group (see ffmerge.hpp).
#include "caps/ffmerge.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>

namespace caps {

namespace {

bool same(double a, double b) { return std::fabs(a - b) <= 1e-9 * std::max(1.0, std::fabs(a)); }

std::string fmt(double x) {
  char b[40];
  std::snprintf(b, sizeof b, "%.4g", x);
  return b;
}

}  // namespace

ForceField merge_forcefields(size_t natoms, const std::vector<FFPart>& parts, const MergeOptions& o, std::vector<std::string>* notes) {
  if (parts.empty()) throw FieldError("no parts to merge");
  if (o.eps_rule != "geometric" && o.eps_rule != "arithmetic") throw FieldError("cross ε rule: geometric or arithmetic");
  if (o.sigma_rule != "arithmetic" && o.sigma_rule != "geometric" && o.sigma_rule != "sixthpower") throw FieldError("cross σ rule: arithmetic, geometric or sixthpower");
  std::vector<std::string> said;
  auto note = [&](const std::string& s) { said.push_back(s); };
  // every atom in exactly one part
  std::vector<int> owner(natoms, -1);
  for (size_t p = 0; p < parts.size(); ++p) {
    const auto& P = parts[p];
    if (!P.ff) throw FieldError("part " + P.tag + " has no force field");
    if (P.ff->type_index.size() != P.atoms.size()) throw FieldError("part " + P.tag + ": its force field is for " + std::to_string(P.ff->type_index.size()) + " atoms, not " + std::to_string(P.atoms.size()));
    for (uint32_t a : P.atoms) {
      if (a >= natoms) throw FieldError("part " + P.tag + ": atom " + std::to_string(a) + " is not in the structure");
      if (owner[a] >= 0) throw FieldError("atom " + std::to_string(a + 1) + " is in two parts (" + parts[size_t(owner[a])].tag + ", " + P.tag + ")");
      owner[a] = int(p);
    }
  }
  for (size_t i = 0; i < natoms; ++i)
    if (owner[i] < 0) throw FieldError("atom " + std::to_string(i + 1) + " is in no part: give every atom a force field");
  // a part without 1-4 pairs (a water model: three bonds never separate two of its atoms) sets no 1-4 scaling
  auto has14 = [](const ForceField& F) { return !F.dihedrals.empty() || !F.dihedrals2.empty() || !F.pairs14.empty() || !F.cbt.empty(); };
  size_t first = 0;   // the settings everybody shares come from the first part a force field types (not a many-body group)
  while (first + 1 < parts.size() && (parts[first].ff->manybody.on() || !has14(*parts[first].ff))) ++first;
  if (parts[first].ff->manybody.on() || !has14(*parts[first].ff)) {
    first = 0;
    while (first + 1 < parts.size() && parts[first].ff->manybody.on()) ++first;
  }
  const ForceField& F0 = *parts[first].ff;
  // what one simulation holds for everybody
  bool mixed_forms = false, own14 = false;
  for (const auto& P : parts) {
    const ForceField& F = *P.ff;
    if (F.manybody.on()) {   // a many-body group: no 1-4 pairs, and the non-bonded settings of the rest
      if (F.pair_form != F0.pair_form) mixed_forms = true;   // its cross Lennard-Jones is 12-6 (UFF)
      continue;
    }
    if (has14(F) && F.keep13 != F0.keep13)
      throw FieldError(F.name + " and " + F0.name + " treat 1-3 pairs differently (one keeps them): they cannot share one simulation");
    if (has14(F) && (!same(F.lj14, F0.lj14) || !same(F.coul14, F0.coul14))) {
      const std::string why = F.name + " scales 1-4 pairs by LJ " + fmt(F.lj14) + ", Coulomb " + fmt(F.coul14) + ", " + F0.name + " by LJ " + fmt(F0.lj14) +
                              ", Coulomb " + fmt(F0.coul14);
      if (o.scaling14 == "first") note("1-4 scaling: " + why + "; the first part's (" + fmt(F0.lj14) + ", " + fmt(F0.coul14) + ") is used for every part, as asked");
      else if (o.scaling14 == "refuse") throw FieldError(why + ": choose force fields of one family, keep each part's own 1-4 scaling (scaling14: own), or take the first part's for all (scaling14: first)");
      else own14 = true;
    }
    if (F.pair_form != F0.pair_form) mixed_forms = true;
    if (F.dielectric != F0.dielectric || F.coul_rf != F0.coul_rf || F.lj_shift != F0.lj_shift || F.lj_fsw != F0.lj_fsw || F.coul_gromacs != F0.coul_gromacs)
      throw FieldError(F.name + " and " + F0.name + " treat the non-bonded terms differently (dielectric, reaction field or switching): they cannot share one simulation");
    if (!F.lj14_types.empty() || F.hbond.on()) throw FieldError(F.name + ": separate 1-4 Lennard-Jones types (CHARMM) and DREIDING hydrogen bonds are not merged with other force fields");
  }
  if (mixed_forms && o.cross96 != "rmin")
    throw FieldError("one part uses Lennard-Jones 9-6 (class II) and another 12-6: the cross pairs need one form. Use force fields of one class (IFF-PCFF covers minerals with PCFF polymers), or give the 9-6 sites a 12-6 form with the same well depth and minimum (cross96: rmin)");
  int sw_parts = 0, mb_parts = 0;
  for (const auto& P : parts) sw_parts += P.ff->sw.on, mb_parts += P.ff->manybody.on();
  if (sw_parts > 1) throw FieldError("Stillinger–Weber in more than one part: one parameter set only");
  if (mb_parts > 1) throw FieldError("a many-body potential file in more than one group: one per system (a file that covers every element, e.g. SiC.tersoff for Si and C, serves several)");
  if (mb_parts && sw_parts) throw FieldError("a many-body potential file together with the mW Stillinger–Weber water: not written together");

  ForceField M;
  M.pair_form = mixed_forms ? "lj9-6" : F0.pair_form;
  M.mixing = F0.mixing;
  M.lj14 = F0.lj14, M.coul14 = F0.coul14, M.keep13 = F0.keep13;
  M.dielectric = F0.dielectric, M.coul_rf = F0.coul_rf, M.eps_rf = F0.eps_rf, M.lj_shift = F0.lj_shift, M.lj_fsw = F0.lj_fsw, M.coul_gromacs = F0.coul_gromacs;
  M.coul_inner = F0.coul_inner, M.lj_inner = F0.lj_inner, M.cutoff = F0.cutoff;
  M.atom_type.assign(natoms, ""), M.why.assign(natoms, ""), M.type_index.assign(natoms, 0), M.charge.assign(natoms, 0), M.mass.assign(natoms, 0);
  M.excluded.assign(natoms, {});
  // the engine styles: kept only when every part declares the same
  bool same_native = true;
  for (const auto& P : parts)
    if (!P.ff->manybody.on()) same_native = same_native && P.ff->native_pair == F0.native_pair && P.ff->native_dihedral == F0.native_dihedral && P.ff->native_improper == F0.native_improper &&
                  P.ff->native_special == F0.native_special;
  if (same_native) M.native_pair = F0.native_pair, M.native_dihedral = F0.native_dihedral, M.native_improper = F0.native_improper, M.native_special = F0.native_special;
  M.native_timestep = F0.native_timestep;
  for (const auto& P : parts) if (P.ff->native_timestep > 0) M.native_timestep = M.native_timestep > 0 ? std::min(M.native_timestep, P.ff->native_timestep) : P.ff->native_timestep;

  // types: each part's in turn; a name that another part already uses gets the part's tag
  std::vector<std::vector<int>> tmap(parts.size());
  std::vector<int> type_part;                 // merged type → part
  std::vector<int> pair14_part;               // merged 1-4 pair → part
  std::vector<bool> type_is96;                // merged type's own form is 9-6
  std::set<std::string> taken;
  for (size_t p = 0; p < parts.size(); ++p) {
    const ForceField& F = *parts[p].ff;
    for (size_t t = 0; t < F.type_names.size(); ++t) {
      std::string name = F.type_names[t];
      if (taken.count(name)) name += "_" + (parts[p].tag.empty() ? std::to_string(p + 1) : parts[p].tag);
      while (taken.count(name)) name += "'";
      taken.insert(name);
      tmap[p].push_back(int(M.type_names.size()));
      M.type_names.push_back(name);
      M.lj.push_back(F.lj[t]);
      type_part.push_back(int(p));
      type_is96.push_back(F.pair_form == "lj9-6");
    }
  }
  // per atom
  for (size_t p = 0; p < parts.size(); ++p) {
    const auto& P = parts[p];
    const ForceField& F = *P.ff;
    for (size_t k = 0; k < P.atoms.size(); ++k) {
      const uint32_t g = P.atoms[k];
      M.type_index[g] = tmap[p][size_t(F.type_index[k])];
      M.atom_type[g] = M.type_names[size_t(M.type_index[g])];
      M.why[g] = k < F.why.size() ? F.why[k] : "";
      M.charge[g] = F.charge[k];
      M.mass[g] = F.mass[k];
      for (uint32_t x : F.excluded[k]) M.excluded[g].push_back(P.atoms[x]);
      std::sort(M.excluded[g].begin(), M.excluded[g].end());
    }
    auto A = [&](uint32_t local) { return P.atoms[local]; };
    for (auto t : F.bonds) { t.i = A(t.i), t.j = A(t.j); M.bonds.push_back(t); }
    for (auto t : F.angles) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k); M.angles.push_back(t); }
    for (auto t : F.dihedrals) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k), t.l = A(t.l); M.dihedrals.push_back(t); }
    for (auto t : F.impropers) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k), t.l = A(t.l); M.impropers.push_back(t); }
    for (auto t : F.impropers_harmonic) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k), t.l = A(t.l); M.impropers_harmonic.push_back(t); }
    for (auto t : F.inversions) { t.c = A(t.c), t.a = A(t.a), t.b = A(t.b), t.d = A(t.d); M.inversions.push_back(t); }
    for (auto t : F.bonds_x) { t.i = A(t.i), t.j = A(t.j); M.bonds_x.push_back(t); }
    for (auto t : F.angles_x) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k); M.angles_x.push_back(t); }
    for (auto t : F.urey_bradley) { t.i = A(t.i), t.k = A(t.k); M.urey_bradley.push_back(t); }
    for (auto t : F.cbt) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k), t.l = A(t.l); M.cbt.push_back(t); }
    for (auto t : F.lj_pairs) { t.i = A(t.i), t.j = A(t.j); M.lj_pairs.push_back(t); }
    for (auto t : F.bonds2) { t.i = A(t.i), t.j = A(t.j); M.bonds2.push_back(t); }
    for (auto t : F.angles2) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k); M.angles2.push_back(t); }
    for (auto t : F.dihedrals2) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k), t.l = A(t.l); M.dihedrals2.push_back(t); }
    for (auto t : F.impropers2) { t.i = A(t.i), t.j = A(t.j), t.k = A(t.k), t.l = A(t.l); M.impropers2.push_back(t); }
    for (auto t : F.pairs14) {
      t[0] = A(t[0]), t[1] = A(t[1]);
      M.pairs14.push_back(t);
      pair14_part.push_back(int(p));
    }
    for (auto v : F.vsites) {
      v.site = A(v.site);
      for (auto& f : v.from) f = A(f);
      M.vsites.push_back(v);
    }
    if (F.manybody.on()) {   // the element of each merged type ("" for the other groups' types: NULL)
      M.manybody = F.manybody;
      M.manybody.element.clear();
      M.manybody.entry.clear();
      for (size_t t = 0; t < F.type_names.size(); ++t) {
        const size_t mt = size_t(tmap[p][t]);
        if (M.manybody.element.size() < mt + 1) M.manybody.element.resize(mt + 1), M.manybody.entry.resize(mt + 1);
        M.manybody.element[mt] = F.manybody.element[t];
        if (t < F.manybody.entry.size()) M.manybody.entry[mt] = F.manybody.entry[t];
      }
    }
    if (F.sw.on) {
      M.sw = F.sw;
      M.sw.atom.assign(natoms, 0);
      for (size_t k = 0; k < F.sw.atom.size() && k < P.atoms.size(); ++k) M.sw.atom[P.atoms[k]] = F.sw.atom[k];
    }
    // the part's own pair terms: explicit pairs and other forms, remapped; its excluded type pairs
    for (const auto& [ab, pt] : F.pair_override) M.pair_override[{tmap[p][size_t(ab.first)], tmap[p][size_t(ab.second)]}] = pt;
    for (const auto& [ab, pf] : F.pair_func) M.pair_func[{tmap[p][size_t(ab.first)], tmap[p][size_t(ab.second)]}] = pf;
    for (const auto& [a, b] : F.excluded_type_pairs) {
      const int x = tmap[p][size_t(a)], y = tmap[p][size_t(b)];
      M.excluded_type_pairs.insert({std::min(x, y), std::max(x, y)});
    }
    for (const auto& nt : F.notes) M.notes.push_back((parts.size() > 1 ? P.tag + ": " : "") + nt);
  }
  // unlike pairs inside a part whose own mixing (or form) differs from the merged one: written out with the part's rule
  const int nt = int(M.type_names.size());
  for (size_t p = 0; p < parts.size(); ++p) {
    const ForceField& F = *parts[p].ff;
    const bool own_form_differs = F.pair_form != M.pair_form;
    if (F.mixing == M.mixing && !own_form_differs) continue;
    for (size_t a = 0; a < F.type_names.size(); ++a)
      for (size_t b = a; b < F.type_names.size(); ++b) {
        const int x = tmap[p][a], y = tmap[p][b];
        const std::pair<int, int> key{std::min(x, y), std::max(x, y)};
        if (M.pair_func.count(key)) continue;
        const PairType pt = mixed_pair(F, int(a), int(b));
        if (own_form_differs) M.pair_func[key] = PairFunc{kPairSdk126, pt.eps, pt.sigma, 0};   // 4ε[(σ/r)¹² − (σ/r)⁶] in a 9-6 system
        else if (a != b) M.pair_override[key] = pt;
      }
  }
  // cross pairs between parts: the chosen rule, or the pair given; 9-6 sites in a 12-6 cross pair keep ε and r_min
  std::map<std::pair<std::string, std::string>, PairType> given;
  for (const auto& c : o.explicit_pairs) given[{std::min(c.a, c.b), std::max(c.a, c.b)}] = {c.eps, c.sigma};
  int cross = 0, given_used = 0;
  for (int a = 0; a < nt; ++a)
    for (int b = a + 1; b < nt; ++b) {
      if (type_part[size_t(a)] == type_part[size_t(b)]) continue;
      ++cross;
      const std::pair<int, int> key{a, b};
      auto site = [&](int t) {   // the site as 12-6 when the cross pair is 12-6
        PairType s = M.lj[size_t(t)];
        if (mixed_forms && type_is96[size_t(t)]) s.sigma /= std::pow(2.0, 1.0 / 6);   // 9-6 σ is r_min
        return s;
      };
      PairType pt;
      if (auto it = given.find({std::min(M.type_names[size_t(a)], M.type_names[size_t(b)]), std::max(M.type_names[size_t(a)], M.type_names[size_t(b)])}); it != given.end()) {
        pt = it->second;
        ++given_used;
      } else {
        const PairType A = site(a), B = site(b);
        if (o.sigma_rule == "sixthpower") {
          const double s6a = std::pow(A.sigma, 6), s6b = std::pow(B.sigma, 6);
          pt.sigma = std::pow(0.5 * (s6a + s6b), 1.0 / 6);
          pt.eps = s6a + s6b > 0 ? 2 * std::sqrt(A.eps * B.eps) * std::pow(A.sigma, 3) * std::pow(B.sigma, 3) / (s6a + s6b) : 0;
        } else {
          pt.sigma = o.sigma_rule == "geometric" ? std::sqrt(A.sigma * B.sigma) : 0.5 * (A.sigma + B.sigma);
          pt.eps = o.eps_rule == "geometric" ? std::sqrt(A.eps * B.eps) : 0.5 * (A.eps + B.eps);
        }
      }
      if (mixed_forms) M.pair_func[key] = PairFunc{kPairSdk126, pt.eps, pt.sigma, 0};
      else M.pair_override[key] = pt;
    }
  if (M.manybody.on()) {
    M.manybody.element.resize(M.type_names.size());
    if (!M.manybody.entry.empty()) M.manybody.entry.resize(M.type_names.size());
  }
  // names and notes
  std::string name;
  for (size_t p = 0; p < parts.size(); ++p) name += (p ? " + " : "") + parts[p].ff->name + (parts.size() > 1 ? " (" + parts[p].tag + ")" : "");
  M.name = name;
  if (parts.size() > 1) {
    note(std::to_string(parts.size()) + " force fields by group; " + std::to_string(cross) + " cross type pairs by ε " +
         (o.sigma_rule == "sixthpower" ? std::string("and σ sixth-power") : o.eps_rule + ", σ " + o.sigma_rule) + " mixing" +
         (given_used ? " (" + std::to_string(given_used) + " given explicitly)" : "") + ", written out explicitly");
    if (mixed_forms)
      note("9-6 and 12-6 parts: every pair with a 12-6 site is 4ε[(σ/r)¹² − (σ/r)⁶] (LAMMPS lj/sdk lj12_6 beside lj/class2); a 9-6 site enters a cross pair with its own ε and minimum r_min (σ₁₂ = r_min / 2^(1/6))");
  }
  for (auto& s : said) M.notes.push_back(s);
  if (notes) *notes = said;
  if (own14) {   // each part's 1-4 pairs by its own scaling
    M.pairs14_lj.resize(M.pairs14.size()), M.pairs14_coul.resize(M.pairs14.size());
    for (size_t k = 0; k < M.pairs14.size(); ++k) {
      const ForceField& F = *parts[size_t(pair14_part[k])].ff;
      M.pairs14_lj[k] = F.lj14, M.pairs14_coul[k] = F.coul14;
    }
    M.type_part = type_part;
    for (const auto& P : parts) M.part14.push_back({P.ff->lj14, P.ff->coul14});
    std::string list;
    for (const auto& P : parts)
      if (has14(*P.ff)) list += (list.empty() ? "" : ", ") + P.ff->name + " LJ " + fmt(P.ff->lj14) + " / Coulomb " + fmt(P.ff->coul14);
    note("1-4 pairs scaled by each part's own force field (" + list + "): LAMMPS writes a pair sub-style per part with its own special weights, GROMACS each 1-4 pair with its own scaling");
  }
  return M;
}

}  // namespace caps
