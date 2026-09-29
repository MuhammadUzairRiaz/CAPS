// CAPS reactions as LAMMPS fix bond/react templates: see caps/bond_react.hpp.
#include "caps/bond_react.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {

// A reaction site cut from the structure: the atoms within `radius` bonds of the atoms whose bonds change.
struct Site {
  int reaction = 0;
  Match match;
  std::vector<uint32_t> atoms;      // structure indices, in template order
  std::string signature;            // the chemical environment (types, bonds, the pattern's roles)
};

std::vector<int> changed_maps(const ReactionTemplate& t) {
  std::set<int> m{t.init_a, t.init_b};
  for (auto [a, b] : t.form) m.insert(a), m.insert(b);
  for (auto [a, b] : t.brk) m.insert(a), m.insert(b);
  for (auto [a, b] : t.move) m.insert(a), m.insert(b);
  for (int a : t.remove) m.insert(a);
  for (int a : t.byproduct) m.insert(a);
  return {m.begin(), m.end()};
}

uint32_t atom_of(const ReactionTemplate& t, const Match& m, int map) {
  for (size_t q = 0; q < t.atoms.size(); ++q)
    if (t.atoms[q].map == map) return m.atoms[q];
  throw ReactError("internal: map " + std::to_string(map) + " not in the match");
}

// The fragment around a match and a label of its environment: Weisfeiler–Lehman refinement of (type, role) labels over
// the fragment's bonds, the sorted final labels as the signature.
Site cut(const System&, const std::vector<std::vector<uint32_t>>& nb, const std::vector<std::string>& type, const ReactionTemplate& t,
         const Match& m, int reaction, int radius) {
  Site site;
  site.reaction = reaction;
  site.match = m;
  std::map<uint32_t, int> dist;
  std::vector<uint32_t> q;
  for (int map : changed_maps(t)) {
    const uint32_t a = atom_of(t, m, map);
    if (dist.emplace(a, 0).second) q.push_back(a);
  }
  for (size_t h = 0; h < q.size(); ++h) {
    const uint32_t u = q[h];
    if (dist[u] >= radius) continue;
    for (uint32_t v : nb[u])
      if (dist.emplace(v, dist[u] + 1).second) q.push_back(v);
  }
  // template order: the pattern atoms first (by map number), then by distance and index
  std::map<uint32_t, int> role;
  for (size_t k = 0; k < t.atoms.size(); ++k) role[m.atoms[k]] = t.atoms[k].map;
  std::vector<uint32_t> order(q.begin(), q.end());
  std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
    const int ra = role.count(a) ? role[a] : 1 << 20, rb = role.count(b) ? role[b] : 1 << 20;
    if (ra != rb) return ra < rb;
    if (dist[a] != dist[b]) return dist[a] < dist[b];
    return a < b;
  });
  site.atoms = order;
  std::map<uint32_t, std::string> label;
  for (uint32_t a : order)
    label[a] = type[a] + "|" + (role.count(a) ? "p" + std::to_string(role[a]) : "d" + std::to_string(dist[a]));
  const std::set<uint32_t> in(order.begin(), order.end());
  for (int it = 0; it < 4; ++it) {
    std::map<uint32_t, std::string> next;
    for (uint32_t a : order) {
      std::vector<std::string> nl;
      for (uint32_t v : nb[a]) nl.push_back(in.count(v) ? label[v] : "edge");
      std::sort(nl.begin(), nl.end());
      std::string x = label[a] + "(";
      for (const auto& l : nl) x += l + ",";
      next[a] = std::to_string(std::hash<std::string>{}(x + ")"));
    }
    label = std::move(next);
  }
  std::vector<std::string> all;
  for (uint32_t a : order) all.push_back(label[a]);
  std::sort(all.begin(), all.end());
  std::string sig = std::to_string(order.size()) + ":";
  for (const auto& l : all) sig += l + ";";
  site.signature = sig;
  return site;
}

std::string clean(std::string n) {
  for (auto& c : n)
    if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
  return n;
}

}  // namespace

BondReactReport write_bond_react(const System& s0, const std::vector<ReactionTemplate>& templates,
                                 const std::function<std::shared_ptr<const ForceField>(const System&)>& field, const std::string& dir,
                                 const std::string& stem, const BondReactOptions& o) {
  if (templates.empty()) throw ReactError("no reaction templates");
  if (!field) throw ReactError("fix bond/react templates need the structure's force field: assign one first");
  if (!o.weights.empty() && o.weights.size() != templates.size()) throw ReactError("one weight per template");
  BondReactReport rep;
  std::filesystem::create_directories(dir);
  System s = s0;
  if (s.cell.valid()) make_molecules_whole(s);   // bonded atoms on the same side: consistent image flags, whole templates
  const size_t n = s.atoms.size();
  const auto ff0 = field(s);
  if (!ff0) throw ReactError("the force field cannot describe the structure");
  // stages: the structure, then — for a reaction whose groups only a first reaction makes (the acid a maleic anhydride
  // leaves when an OH opens it) — copies of it with another reaction applied at one site; each reaction's templates are
  // cut from the first stage where it has pairs, typed there
  struct Stage { System sys; std::shared_ptr<const ForceField> ff; std::vector<std::vector<uint32_t>> nb; size_t offset = 0; std::string after; };
  std::vector<Stage> stages;
  stages.push_back({s, ff0, s.neighbours(), 0, ""});
  struct Group { std::vector<Site> sites; };
  std::vector<std::map<std::string, Group>> groups(templates.size());
  std::vector<int> stage_of(templates.size(), -1);
  auto survey = [&](size_t k, int st) {
    const Stage& S = stages[size_t(st)];
    const double cap = o.survey_capture > 0 ? o.survey_capture : std::max(templates[k].capture, 8.0);
    int found = 0;
    for (const auto& m : find_matches(S.sys, templates[k], int(k), {}, cap)) {
      Site site = cut(S.sys, S.nb, S.ff->atom_type, templates[k], m, int(k), o.radius);
      groups[k][site.signature].sites.push_back(std::move(site));
      ++found;
    }
    if (found) stage_of[k] = st;
    rep.candidates += found;
    return found;
  };
  for (size_t k = 0; k < templates.size(); ++k) survey(k, 0);
  for (size_t k = 0; k < templates.size(); ++k) {
    if (stage_of[k] >= 0) continue;
    for (size_t k2 = 0; k2 < templates.size() && stage_of[k] < 0; ++k2) {
      if (k2 == k || stage_of[k2] != 0 || groups[k2].empty()) continue;
      // the most frequent environment of the first reaction, applied once
      const Group* best = nullptr;
      for (const auto& [sig, g] : groups[k2]) if (!best || g.sites.size() > best->sites.size()) best = &g;
      Stage S;
      S.sys = s;
      std::vector<ReactionTemplate> one{templates[k2]};
      Match m = best->sites.front().match;
      m.reaction = 0;
      apply_matches(S.sys, one, {m}, o.keep_byproducts);
      if (S.sys.cell.valid()) make_molecules_whole(S.sys);
      try {
        S.ff = field(S.sys);
      } catch (const std::exception& e) {
        throw ReactError(std::string("the force field cannot describe the product of ") + templates[k2].name + ": " + e.what());
      }
      S.nb = S.sys.neighbours();
      S.after = templates[k2].name;
      stages.push_back(std::move(S));
      if (survey(k, int(stages.size()) - 1))
        rep.notes.push_back("reaction " + templates[k].name + ": its groups appear only after " + templates[k2].name + " — templates cut from a copy where that reaction has happened");
      else
        stages.pop_back();
    }
  }

  // the variants: the most frequent environments of each reaction, each with its reacted copy
  struct Variant { Site site; int stage = 0; System post; std::vector<int64_t> post_of; std::string name; };   // post_of: pre atom → post index (−1 deleted)
  std::vector<Variant> variants;
  for (size_t k = 0; k < templates.size(); ++k) {
    std::vector<std::pair<size_t, const Group*>> by;
    for (const auto& [sig, g] : groups[k]) by.push_back({g.sites.size(), &g});
    std::stable_sort(by.begin(), by.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    if (by.empty()) {
      rep.notes.push_back("reaction " + templates[k].name + ": no reactive pair in the structure (nor after another reaction) within the survey distance — no template written");
      continue;
    }
    const Stage& S = stages[size_t(stage_of[k])];
    const size_t sn = S.sys.atoms.size();
    for (size_t v = 0; v < by.size() && int(v) < o.max_variants; ++v) {
      Variant var;
      var.site = by[v].second->sites.front();
      var.stage = stage_of[k];
      var.name = clean(templates[k].name) + "_" + std::to_string(v + 1);
      // the reaction applied at that site; deleted atoms (and byproducts, unless kept) leave the copy — they stay in the
      // post template with their pre-reaction types as DeleteIDs — the rest keep their order
      var.post = S.sys;
      std::vector<ReactionTemplate> one{templates[k]};
      one[0].keep_charges = false;
      Match m = var.site.match;
      m.reaction = 0;
      apply_matches(var.post, one, {m}, o.keep_byproducts);
      std::vector<char> dead(sn, 0);
      for (int map : templates[k].remove) dead[atom_of(templates[k], var.site.match, map)] = 1;
      if (!o.keep_byproducts)
        for (int map : templates[k].byproduct) dead[atom_of(templates[k], var.site.match, map)] = 1;
      var.post_of.assign(sn, -1);
      int64_t next = 0;
      for (size_t i = 0; i < sn; ++i)
        if (!dead[i]) var.post_of[i] = next++;
      if (size_t(next) != var.post.atoms.size()) throw ReactError("internal: the reacted copy has " + std::to_string(var.post.atoms.size()) + " atoms, expected " + std::to_string(next));
      // LAMMPS updates types and charges only inside the template, and only on atoms whose 1-2 and 1-3 neighbours are all
      // in it ("landlocked"): the template reaches three bonds past the farthest atom whose type or charge the reaction changes
      {
        std::shared_ptr<const ForceField> ffp;
        try {
          ffp = field(var.post);
        } catch (const std::exception& e) {
          throw ReactError(std::string("the force field cannot describe a reacted site of ") + templates[k].name + ": " + e.what() +
                           " — add the missing parameters (Force field › Fill gaps)");
        }
        std::map<uint32_t, int> dist;
        std::vector<uint32_t> q;
        for (int map : changed_maps(templates[k])) {
          const uint32_t a = atom_of(templates[k], var.site.match, map);
          if (dist.emplace(a, 0).second) q.push_back(a);
        }
        for (size_t h = 0; h < q.size(); ++h)
          if (dist[q[h]] < 12)
            for (uint32_t w : S.nb[q[h]])
              if (dist.emplace(w, dist[q[h]] + 1).second) q.push_back(w);
        int far = 0;
        for (size_t i = 0; i < sn; ++i) {
          if (var.post_of[i] < 0) continue;
          const size_t j = size_t(var.post_of[i]);
          if (ffp->atom_type[j] == S.ff->atom_type[i] && std::fabs(ffp->charge[j] - S.ff->charge[i]) < 1e-6) continue;
          auto it = dist.find(uint32_t(i));
          if (it == dist.end()) throw ReactError("the reaction " + templates[k].name + " changes atom " + std::to_string(i + 1) + " more than 12 bonds away (type or charge)");
          far = std::max(far, it->second);
        }
        const int need = std::max(o.radius, far + 3);
        if (need > o.radius) var.site = cut(S.sys, S.nb, S.ff->atom_type, templates[k], var.site.match, int(k), need);
      }
      rep.covered += int(by[v].first);
      rep.variants.push_back({templates[k].name, var.name, int(by[v].first), int(var.site.atoms.size()), 0, 0});
      variants.push_back(std::move(var));
    }
    if (int(by.size()) > o.max_variants) {
      size_t left = 0;
      for (size_t v = size_t(o.max_variants); v < by.size(); ++v) left += by[v].first;
      rep.notes.push_back("reaction " + templates[k].name + ": " + std::to_string(by.size() - size_t(o.max_variants)) + " rarer environments (" +
                          std::to_string(left) + " candidate pairs) have no template; raise the number of templates to cover them");
    }
  }
  if (variants.empty()) throw ReactError("no reactive pair of any template in the structure: nothing to write");

  // the union: the structure, then each reacted copy — one force field over all of it numbers every type the reactions
  // create, and the data file keeps only the structure's atoms
  System u = s;
  auto append = [&](const System& x) {
    const uint32_t off = uint32_t(u.atoms.size());
    int64_t mol_off = 0;
    for (const auto& a : u.atoms) mol_off = std::max(mol_off, a.mol);
    for (auto a : x.atoms) { a.mol += mol_off; u.atoms.push_back(a); }
    for (auto b : x.bonds) { b.i += off; b.j += off; u.bonds.push_back(b); }
    return size_t(off);
  };
  for (size_t k = 1; k < stages.size(); ++k) stages[k].offset = append(stages[k].sys);
  std::vector<size_t> offset;
  for (const auto& v : variants) offset.push_back(append(v.post));
  u.has_charges = false;
  u.velocities.clear();
  std::shared_ptr<const ForceField> ffu;
  try {
    ffu = field(u);
  } catch (const std::exception& e) {
    throw ReactError(std::string("the force field cannot describe a reacted site: ") + e.what() + " — add the missing parameters (Force field › Fill gaps)");
  }
  if (!ffu) throw ReactError("the force field cannot describe a reacted site");
  for (size_t i = 0; i < n; ++i)
    if (ffu->atom_type[i] != ff0->atom_type[i])
      throw ReactError("internal: the structure's types changed when the reacted copies were added (atom " + std::to_string(i + 1) + ")");
  const LammpsTerms terms = lammps_terms(u, *ffu, o.energy, o.style);

  const std::string base = (std::filesystem::path(dir) / stem).string();
  LammpsStyle st = o.style;
  st.write_atoms = n;
  write_lammps_data_ff(u, *ffu, o.energy, base + ".data", true, st);
  rep.files.push_back(base + ".data");

  // the molecule templates and maps
  auto unwrap = [&](const System& x, uint32_t i, const Vec3& ref) {
    const Vec3 d = x.cell.valid() ? x.cell.minimum_image(x.atoms[i].pos - ref) : x.atoms[i].pos - ref;
    return ref + d;
  };
  std::vector<std::string> molecule_lines, react_lines;
  double wmax = 0;
  for (double w : o.weights) wmax = std::max(wmax, w);
  for (size_t v = 0; v < variants.size(); ++v) {
    const Variant& var = variants[v];
    const Stage& S = stages[size_t(var.stage)];
    const ReactionTemplate& t = templates[size_t(var.site.reaction)];
    const auto& pre = var.site.atoms;
    std::map<uint32_t, int> id;   // structure index → template id (1-based)
    for (size_t k = 0; k < pre.size(); ++k) id[pre[k]] = int(k + 1);
    // post: the same atoms in the same order (a deleted atom keeps its pre-reaction type and place: LAMMPS removes it)
    std::vector<int64_t> post_u(pre.size(), -1);   // index in u of each template atom after the reaction
    for (size_t k = 0; k < pre.size(); ++k)
      if (var.post_of[pre[k]] >= 0) post_u[k] = int64_t(offset[v]) + var.post_of[pre[k]];
    std::set<int> deleted;
    for (int map : t.remove) deleted.insert(id.at(atom_of(t, var.site.match, map)));
    if (!o.keep_byproducts)
      for (int map : t.byproduct) deleted.insert(id.at(atom_of(t, var.site.match, map)));
    // edge atoms: a bond to an atom outside the template
    std::vector<int> edge;
    for (size_t k = 0; k < pre.size(); ++k)
      for (uint32_t w : S.nb[pre[k]])
        if (!id.count(w)) { edge.push_back(int(k + 1)); break; }
    rep.variants[v].edge = int(edge.size());
    rep.variants[v].deleted = int(deleted.size());
    const Vec3 ref = S.sys.atoms[atom_of(t, var.site.match, t.init_a)].pos;

    auto write_mol = [&](const std::string& path, bool after) {
      // the template atoms as indices of u, and their terms among themselves
      std::vector<int64_t> at(pre.size());
      for (size_t k = 0; k < pre.size(); ++k) at[k] = after && post_u[k] >= 0 ? post_u[k] : int64_t(S.offset + pre[k]);
      std::map<int64_t, int> tid;
      for (size_t k = 0; k < at.size(); ++k) tid[at[k]] = int(k + 1);
      auto collect = [&](const std::vector<std::pair<int, std::vector<uint32_t>>>& all) {
        std::vector<std::pair<int, std::vector<int>>> out;
        for (const auto& [ty, xs] : all) {
          std::vector<int> ids;
          for (uint32_t x : xs) {
            auto it = tid.find(int64_t(x));
            if (it == tid.end()) break;
            ids.push_back(it->second);
          }
          if (ids.size() == xs.size()) out.push_back({ty, ids});
        }
        return out;
      };
      const auto B = collect(terms.bonds), A = collect(terms.angles), D = collect(terms.dihedrals), I = collect(terms.impropers);
      std::ofstream f(path);
      if (!f) throw std::runtime_error("cannot write " + path);
      char b[256];
      f << "# CAPS " << (after ? "post" : "pre") << "-reaction template: " << t.name << " (" << var.name << "), " << ffu->name << "\n\n";
      f << pre.size() << " atoms\n";
      if (!B.empty()) f << B.size() << " bonds\n";
      if (!A.empty()) f << A.size() << " angles\n";
      if (!D.empty()) f << D.size() << " dihedrals\n";
      if (!I.empty()) f << I.size() << " impropers\n";
      f << "\nCoords\n\n";
      for (size_t k = 0; k < pre.size(); ++k) {
        const Vec3 p = after && post_u[k] >= 0 ? unwrap(var.post, uint32_t(var.post_of[pre[k]]), ref) : unwrap(S.sys, pre[k], ref);
        std::snprintf(b, sizeof b, "%zu %.6f %.6f %.6f\n", k + 1, p[0], p[1], p[2]);
        f << b;
      }
      f << "\nTypes\n\n";
      for (size_t k = 0; k < pre.size(); ++k) f << k + 1 << " " << terms.atom_type[size_t(at[k])] << "  # " << ffu->atom_type[size_t(at[k])] << "\n";
      f << "\nCharges\n\n";
      for (size_t k = 0; k < pre.size(); ++k) {
        std::snprintf(b, sizeof b, "%zu %.6f\n", k + 1, ffu->charge[size_t(at[k])]);
        f << b;
      }
      auto section = [&](const char* name, const std::vector<std::pair<int, std::vector<int>>>& xs) {
        if (xs.empty()) return;
        f << "\n" << name << "\n\n";
        for (size_t q = 0; q < xs.size(); ++q) {
          f << q + 1 << " " << xs[q].first;
          for (int x : xs[q].second) f << " " << x;
          f << "\n";
        }
      };
      section("Bonds", B);
      section("Angles", A);
      section("Dihedrals", D);
      section("Impropers", I);
    };
    const std::string pre_path = base + "_" + var.name + "_pre.mol", post_path = base + "_" + var.name + "_post.mol",
                      map_path = base + "_" + var.name + "_map.txt";
    write_mol(pre_path, false);
    write_mol(post_path, true);
    {
      std::ofstream f(map_path);
      if (!f) throw std::runtime_error("cannot write " + map_path);
      f << "# CAPS map: " << t.name << " (" << var.name << ")\n\n";
      f << edge.size() << " edgeIDs\n";
      f << pre.size() << " equivalences\n";
      if (!deleted.empty()) f << deleted.size() << " deleteIDs\n";
      f << "\nInitiatorIDs\n\n" << id.at(atom_of(t, var.site.match, t.init_a)) << "\n" << id.at(atom_of(t, var.site.match, t.init_b)) << "\n";
      if (!edge.empty()) {
        f << "\nEdgeIDs\n\n";
        for (int e : edge) f << e << "\n";
      }
      if (!deleted.empty()) {
        f << "\nDeleteIDs\n\n";
        for (int d : deleted) f << d << "\n";
      }
      f << "\nEquivalences\n\n";
      for (size_t k = 0; k < pre.size(); ++k) f << k + 1 << " " << k + 1 << "\n";
    }
    rep.files.push_back(pre_path);
    rep.files.push_back(post_path);
    rep.files.push_back(map_path);
    const std::string fn = std::filesystem::path(pre_path).filename().string(), fp = std::filesystem::path(post_path).filename().string(),
                      fm = std::filesystem::path(map_path).filename().string();
    molecule_lines.push_back("molecule        " + var.name + "_pre " + fn);
    molecule_lines.push_back("molecule        " + var.name + "_post " + fp);
    const size_t k = size_t(var.site.reaction);
    double prob = t.probability * (wmax > 0 ? o.weights[k] / wmax : 1.0);
    char b[512];
    std::snprintf(b, sizeof b, "  react %s all %d 0.0 %.4g %s_pre %s_post %s prob %.4g %llu%s", var.name.c_str(), o.nevery, o.rmax > 0 ? std::min(o.rmax, t.capture) : t.capture, var.name.c_str(),
                  var.name.c_str(), fm.c_str(), prob, static_cast<unsigned long long>(o.seed + v), o.between_chains ? " molecule inter" : "");
    react_lines.push_back(b);
  }

  // the input: CAPS's own for this force field (styles, pair coefficients in the data file), then the fix
  {
    LammpsRun run;
    run.kind = LammpsRun::Kind::None;
    run.temperature = o.temperature;
    run.dt = o.timestep;
    std::vector<std::string> notes;
    write_lammps_input(u, *ffu, o.energy, std::filesystem::path(base + ".data").filename().string(), base + ".in", 0, false, run, o.style, &notes);
    std::ifstream in(base + ".in");
    std::stringstream ss;
    ss << in.rdbuf();
    std::string text = ss.str(), out;
    std::istringstream lines(text);
    const double dt = lammps_timestep(run, *ffu);
    for (std::string line; std::getline(lines, line);) {
      if (line.rfind("read_data", 0) == 0) {
        // room for the bonds, angles, dihedrals and special neighbours the reactions add
        line += " extra/bond/per/atom 4 extra/angle/per/atom 12 extra/dihedral/per/atom 36 extra/improper/per/atom 12 extra/special/per/atom 48";
      }
      out += line + "\n";
    }
    char b[512];
    out += "\n# ---- reactions (fix bond/react, REACTER): templates cut from this structure by CAPS, typed with " + ffu->name + "\n";
    for (const auto& l : molecule_lines) out += l + "\n";
    out += "\nfix             rxns all bond/react stabilization yes statted_grp " + std::to_string(o.stabilization_xmax).substr(0, 5) + " &\n";
    for (size_t k = 0; k < react_lines.size(); ++k) out += react_lines[k] + (k + 1 < react_lines.size() ? " &\n" : "\n");
    std::snprintf(b, sizeof b,
                  "\nvelocity        all create %.2f %llu dist gaussian\n"
                  "fix             nvt statted_grp_REACT nvt temp %.2f %.2f %.1f\n"
                  "fix             rescale bond_react_MASTER_group temp/rescale 1 %.2f %.2f 10 1\n"
                  "timestep        %.4g\n"
                  "thermo_style    custom step temp pe press density",
                  o.temperature, static_cast<unsigned long long>(o.seed), o.temperature, o.temperature, 100 * dt, o.temperature, o.temperature, dt);
    out += b;
    for (size_t k = 0; k < react_lines.size(); ++k) out += " f_rxns[" + std::to_string(k + 1) + "]";
    out += "\nthermo          " + std::to_string(std::max(1, o.nevery)) + "\n";
    out += "dump            traj all custom " + std::to_string(std::max<int64_t>(1, o.steps / 50)) + " " + stem + "_react.lammpstrj id mol type q x y z\n";
    out += "run             " + std::to_string(o.steps) + "\n";
    out += "write_data      " + stem + "_reacted.data\n";
    std::ofstream f(base + ".in");
    f << out;
    rep.files.push_back(base + ".in");
    for (const auto& note : notes) rep.notes.push_back(note);
  }
  char b[256];
  std::snprintf(b, sizeof b, "%zu templates for %zu reactions cover %d of %d candidate pairs (within the survey distance) in this structure", variants.size(),
                templates.size(), rep.covered, rep.candidates);
  rep.notes.insert(rep.notes.begin(), b);
  if (o.between_chains)
    rep.notes.push_back("molecule inter: fix bond/react resets molecule ids to the bonded pieces after each reaction (reset_mol_ids, its default), so "
                        "initiators must sit in pieces not yet joined — a crosslinker cannot close back on its own chain, but two chains already "
                        "joined take no second link (a tree-like network; CAPS's own runs count chains and allow it)");
  return rep;
}

// ---------------------------------------------------------------------------------------------------------------
// fix bond/react templates read back

namespace {

struct MolFile {
  std::vector<int> type;
  std::vector<std::string> label;              // type labels, when the file names types instead of numbering them
  std::map<int, double> mass;                  // Masses section (type → mass)
  std::vector<double> atom_mass;               // per atom, when the file gives masses per atom
  std::vector<std::pair<int, int>> bonds;      // 0-based
  int natoms = 0;
};

bool is_section(const std::string& w) {
  static const std::set<std::string> k = {"Coords", "Types", "Charges", "Bonds", "Angles", "Dihedrals", "Impropers", "Masses", "Molecules",
                                          "Special", "Diameters", "Fragments", "Shake", "Body"};
  return k.count(w) > 0;
}

std::vector<std::string> words(const std::string& line) {
  std::string l = line.substr(0, line.find('#'));
  std::istringstream is(l);
  std::vector<std::string> w;
  for (std::string x; is >> x;) w.push_back(x);
  return w;
}

MolFile read_mol(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReactError("cannot read " + path);
  MolFile m;
  std::string line, section;
  std::getline(in, line);   // title
  while (std::getline(in, line)) {
    const auto w = words(line);
    if (w.empty()) continue;
    if (is_section(w[0])) { section = w[0]; continue; }
    if (section.empty()) {
      if (w.size() >= 2 && w[1] == "atoms") m.natoms = std::stoi(w[0]);
      continue;
    }
    if (section == "Types" && w.size() >= 2) {
      const int i = std::stoi(w[0]) - 1;
      if (i < 0) throw ReactError(path + ": atom ids start at 1");
      if (m.type.size() <= size_t(i)) m.type.resize(size_t(i) + 1, 0), m.label.resize(size_t(i) + 1);
      char* end = nullptr;
      const long t = std::strtol(w[1].c_str(), &end, 10);
      if (end && *end == 0) m.type[size_t(i)] = int(t);
      else m.label[size_t(i)] = w[1];
    } else if (section == "Masses" && w.size() >= 2) {
      const int i = std::stoi(w[0]) - 1;
      if (m.atom_mass.size() <= size_t(i)) m.atom_mass.resize(size_t(i) + 1, 0);
      m.atom_mass[size_t(i)] = std::stod(w[1]);
    } else if (section == "Bonds" && w.size() >= 4) {
      m.bonds.push_back({std::stoi(w[2]) - 1, std::stoi(w[3]) - 1});
    }
  }
  if (m.natoms <= 0) m.natoms = int(m.type.size());
  if (int(m.type.size()) < m.natoms) m.type.resize(size_t(m.natoms), 0), m.label.resize(size_t(m.natoms));
  return m;
}

std::map<int, double> data_masses(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReactError("cannot read " + path);
  std::map<int, double> m;
  std::string line;
  bool on = false;
  while (std::getline(in, line)) {
    const auto w = words(line);
    if (w.empty()) continue;
    if (w[0] == "Masses") { on = true; continue; }
    if (on) {
      if (!std::isdigit(static_cast<unsigned char>(w[0][0]))) break;
      m[std::stoi(w[0])] = std::stod(w[1]);
    }
  }
  if (m.empty()) throw ReactError(path + " has no Masses section");
  return m;
}

int element_of_mass(double m) {
  int best = 0;
  double d = 1e9;
  for (int z = 1; z <= 100; ++z) {
    const double dz = std::fabs(element(z).mass - m);
    if (dz < d) d = dz, best = z;
  }
  return d < 0.6 ? best : 0;
}

}  // namespace

ReactionTemplate read_bond_react(const std::string& pre_path, const std::string& post_path, const std::string& map_path, const std::string& masses_from,
                                 const std::string& name, double capture, std::vector<std::string>* notes) {
  const MolFile pre = read_mol(pre_path), post = read_mol(post_path);
  std::map<int, double> masses;
  if (!masses_from.empty()) masses = data_masses(masses_from);
  // the map
  std::vector<int> init, edge, del;
  std::map<int, int> eq;   // pre id → post id (0-based)
  {
    std::ifstream in(map_path);
    if (!in) throw ReactError("cannot read " + map_path);
    std::string line, section;
    std::getline(in, line);
    while (std::getline(in, line)) {
      const auto w = words(line);
      if (w.empty()) continue;
      if (w[0] == "InitiatorIDs" || w[0] == "EdgeIDs" || w[0] == "DeleteIDs" || w[0] == "Equivalences" || w[0] == "CreateIDs" || w[0] == "ChiralIDs" ||
          w[0] == "Constraints") {
        section = w[0];
        if ((section == "CreateIDs") ) throw ReactError(map_path + ": atoms created by the reaction (CreateIDs) are not supported by CAPS templates");
        if (notes && (section == "ChiralIDs" || section == "Constraints")) notes->push_back(section + " in the map are not carried into the CAPS template");
        continue;
      }
      if (section.empty() || !std::isdigit(static_cast<unsigned char>(w[0][0]))) continue;
      const int a = std::stoi(w[0]) - 1;
      if (section == "InitiatorIDs") init.push_back(a);
      else if (section == "EdgeIDs") edge.push_back(a);
      else if (section == "DeleteIDs") del.push_back(a);
      else if (section == "Equivalences" && w.size() >= 2) eq[a] = std::stoi(w[1]) - 1;
    }
  }
  if (init.size() != 2) throw ReactError(map_path + ": needs two InitiatorIDs");
  const int n = pre.natoms;
  if (int(eq.size()) != n) throw ReactError(map_path + ": every pre-reaction atom needs an equivalence (" + std::to_string(eq.size()) + " of " + std::to_string(n) + ")");
  // elements
  std::vector<int> el(static_cast<size_t>(n), 0);
  for (int i = 0; i < n; ++i) {
    double m = 0;
    if (size_t(i) < pre.atom_mass.size() && pre.atom_mass[size_t(i)] > 0) m = pre.atom_mass[size_t(i)];
    else if (pre.type[size_t(i)] > 0 && masses.count(pre.type[size_t(i)])) m = masses.at(pre.type[size_t(i)]);
    el[size_t(i)] = m > 0 ? element_of_mass(m) : 0;
    if (!el[size_t(i)] && !pre.label[size_t(i)].empty()) {   // a type label: its leading letters, when they name an element
      std::string sym;
      for (char c : pre.label[size_t(i)]) { if (!std::isalpha(static_cast<unsigned char>(c))) break; sym += c; }
      for (size_t k = std::min<size_t>(2, sym.size()); k > 0 && !el[size_t(i)]; --k) {
        std::string x = sym.substr(0, k);
        x[0] = char(std::toupper(static_cast<unsigned char>(x[0])));
        if (k == 2) x[1] = char(std::tolower(static_cast<unsigned char>(x[1])));
        el[size_t(i)] = element_from_symbol(x);
      }
    }
    if (!el[size_t(i)])
      throw ReactError("atom " + std::to_string(i + 1) + " of " + pre_path + ": no element (give the data file whose Masses number these types)");
  }
  // bonds before and after, in pre ids
  std::map<int, int> back;
  for (auto [a, b] : eq) back[b] = a;
  auto key = [](int a, int b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
  std::set<std::pair<int, int>> before, after;
  for (auto [a, b] : pre.bonds) before.insert(key(a, b));
  const std::set<int> dead(del.begin(), del.end());
  for (auto [a, b] : post.bonds) {
    if (!back.count(a) || !back.count(b)) throw ReactError(post_path + ": a bond to an atom with no equivalence");
    after.insert(key(back[a], back[b]));
  }
  std::vector<std::vector<int>> nb(static_cast<size_t>(n));
  for (auto [a, b] : before) nb[size_t(a)].push_back(b), nb[size_t(b)].push_back(a);
  std::vector<std::pair<int, int>> form, brk, move;
  for (const auto& e : after) if (!before.count(e) && !dead.count(e.first) && !dead.count(e.second)) form.push_back(e);
  for (const auto& e : before) if (!after.count(e) && !dead.count(e.first) && !dead.count(e.second)) brk.push_back(e);
  // a hydrogen that leaves one partner for another: a move
  for (int h = 0; h < n; ++h) {
    if (el[size_t(h)] != 1 || dead.count(h)) continue;
    std::vector<size_t> f, b;
    for (size_t k = 0; k < form.size(); ++k) if (form[k].first == h || form[k].second == h) f.push_back(k);
    for (size_t k = 0; k < brk.size(); ++k) if (brk[k].first == h || brk[k].second == h) b.push_back(k);
    if (f.size() != 1 || b.size() != 1) continue;
    const int to = form[f[0]].first == h ? form[f[0]].second : form[f[0]].first;
    move.push_back({h, to});
    form.erase(form.begin() + long(f[0]));
    brk.erase(brk.begin() + long(b[0]));
  }
  // the pattern: the atoms whose bonds change, and the heavy atoms within two bonds of them
  std::set<int> changed(init.begin(), init.end());
  for (auto [a, b] : form) changed.insert(a), changed.insert(b);
  for (auto [a, b] : brk) changed.insert(a), changed.insert(b);
  for (auto [a, b] : move) {
    changed.insert(a), changed.insert(b);
    for (int x : nb[size_t(a)]) changed.insert(x);   // the moved atom's old partner
  }
  for (int d : del) changed.insert(d);
  std::map<int, int> dist;
  std::vector<int> q;
  for (int c : changed) { dist[c] = 0; q.push_back(c); }
  for (size_t h = 0; h < q.size(); ++h) {
    const int u = q[h];
    if (dist[u] >= 2) continue;
    for (int v : nb[size_t(u)])
      if (el[size_t(v)] != 1 && !dist.count(v)) { dist[v] = dist[u] + 1; q.push_back(v); }
  }
  // post components with no changed-chemistry partner in the network: byproducts (water, H2 …) — small pieces cut off
  std::vector<int> byproduct;
  {
    std::vector<std::vector<int>> pn(static_cast<size_t>(n));
    for (const auto& e : after) pn[size_t(e.first)].push_back(e.second), pn[size_t(e.second)].push_back(e.first);
    std::vector<int> comp(static_cast<size_t>(n), -1);
    int nc = 0;
    std::vector<int> size;
    for (int i = 0; i < n; ++i) {
      if (comp[size_t(i)] >= 0 || dead.count(i)) continue;
      std::vector<int> st{i};
      comp[size_t(i)] = nc;
      int sz = 0;
      while (!st.empty()) {
        const int u = st.back();
        st.pop_back();
        ++sz;
        for (int v : pn[size_t(u)]) if (comp[size_t(v)] < 0) { comp[size_t(v)] = nc; st.push_back(v); }
      }
      size.push_back(sz);
      ++nc;
    }
    const std::set<int> edges(edge.begin(), edge.end());
    for (int i = 0; i < n; ++i) {
      if (dead.count(i)) continue;
      const int c = comp[size_t(i)];
      bool open = false;   // a component reaching the edge is part of the network
      for (int j = 0; j < n; ++j) if (comp[size_t(j)] == c && edges.count(j)) open = true;
      bool cut_off = false;   // it lost a bond in the reaction
      for (auto [a, b] : brk) if (comp[size_t(a)] == c || comp[size_t(b)] == c) cut_off = true;
      for (auto [a, b] : move) if (comp[size_t(a)] == c) cut_off = true;
      if (!open && cut_off && size[size_t(c)] <= 8) byproduct.push_back(i);
    }
    for (int b : byproduct) changed.insert(b), dist.emplace(b, 0);
  }
  std::set<int> in;
  for (auto [a, d] : dist) in.insert(a);
  for (int c : changed) in.insert(c);
  // order: from the initiators through pattern bonds, so each atom is bonded to one defined before it
  std::vector<int> order;
  std::set<int> placed;
  for (int s0 : init) {
    if (placed.count(s0)) continue;
    std::vector<int> bfs{s0};
    placed.insert(s0);
    for (size_t h = 0; h < bfs.size(); ++h)
      for (int v : nb[size_t(bfs[h])])
        if (in.count(v) && !placed.count(v)) { placed.insert(v); bfs.push_back(v); }
    order.insert(order.end(), bfs.begin(), bfs.end());
  }
  for (int a : in)
    if (!placed.count(a)) throw ReactError("pre-reaction atom " + std::to_string(a + 1) + " is not bonded to an initiator's piece of the template");
  std::map<int, int> map_of;
  for (size_t k = 0; k < order.size(); ++k) map_of[order[k]] = int(k + 1);
  const std::set<int> edges(edge.begin(), edge.end());
  ReactionTemplate t;
  t.name = name.empty() ? std::filesystem::path(map_path).stem().string() : name;
  for (int a : order) {
    TemplateAtom x;
    x.map = map_of[a];
    x.element = el[size_t(a)];
    int h = 0;
    for (int v : nb[size_t(a)]) h += el[size_t(v)] == 1;
    if (!edges.count(a)) {
      x.degree = int(nb[size_t(a)].size());
      x.h_min = x.h_max = h;
    } else if (h > 0) {
      x.h_min = h;
    }
    for (int v : nb[size_t(a)]) {   // three-membered rings (epoxides)
      for (int w : nb[size_t(v)])
        if (w != a && std::find(nb[size_t(a)].begin(), nb[size_t(a)].end(), w) != nb[size_t(a)].end()) x.ring3 = true;
    }
    for (int v : nb[size_t(a)])
      if (map_of.count(v) && map_of[v] < x.map) x.bonded.push_back(map_of[v]);
    t.atoms.push_back(x);
  }
  for (auto& a : t.atoms)
    for (int m : std::vector<int>(a.bonded))
      for (auto& b : t.atoms) if (b.map == m) b.bonded.push_back(a.map);
  t.init_a = map_of[init[0]];
  t.init_b = map_of[init[1]];
  t.capture = capture > 0 ? capture : 5.0;
  t.min_path = 0;
  for (auto [a, b] : form) t.form.push_back({map_of[a], map_of[b]});
  for (auto [a, b] : brk) t.brk.push_back({map_of[a], map_of[b]});
  for (auto [a, b] : move) t.move.push_back({map_of[a], map_of[b]});
  for (int d : del) t.remove.push_back(map_of[d]);
  for (int b : byproduct) t.byproduct.push_back(map_of[b]);
  t.sites = {t.init_a};
  // through text, so the parser checks it as it checks any template
  auto parsed = parse_templates(template_text(t));
  if (notes) {
    notes->push_back("pattern: " + std::to_string(t.atoms.size()) + " of the " + std::to_string(n) + " pre-reaction atoms (the changing atoms and the heavy atoms within two bonds), by element, hydrogens and bonds — LAMMPS matches by atom type");
    if (!byproduct.empty()) notes->push_back(std::to_string(byproduct.size()) + " atoms leave as a byproduct molecule (kept or removed by the run)");
  }
  return parsed.front();
}

std::string template_text(const ReactionTemplate& t) {
  std::ostringstream o;
  o << "reaction " << clean(t.name) << "\n";
  for (const auto& a : t.atoms) {
    o << "atom " << a.map << " " << element(a.element).symbol;
    if (a.ring3) o << " ring3";
    if (a.not_aromatic) o << " not_aromatic";
    if (a.aromatic) o << " aromatic";
    if (a.h_min >= 0 && a.h_min == a.h_max) o << " H=" << a.h_min;
    else {
      if (a.h_min >= 0) o << " H>=" << a.h_min;
      if (a.h_max >= 0) o << " H<=" << a.h_max;
    }
    if (a.degree >= 0) o << " degree=" << a.degree;
    std::vector<int> earlier;
    for (int m : a.bonded) if (m < a.map) earlier.push_back(m);
    if (!earlier.empty()) {
      o << " bonded";
      for (int m : earlier) o << " " << m;
    }
    o << "\n";
  }
  o << "initiators " << t.init_a << " " << t.init_b << "\n";
  o << "capture " << t.capture << "\nprobability " << t.probability << "\nmin_path " << t.min_path << "\n";
  for (auto [a, b] : t.form) o << "form " << a << " " << b << "\n";
  for (auto [a, b] : t.brk) o << "break " << a << " " << b << "\n";
  for (auto [a, b] : t.move) o << "move " << a << " " << b << "\n";
  if (!t.remove.empty()) { o << "delete"; for (int d : t.remove) o << " " << d; o << "\n"; }
  if (!t.byproduct.empty()) { o << "byproduct"; for (int d : t.byproduct) o << " " << d; o << "\n"; }
  if (!t.sites.empty()) { o << "sites"; for (int d : t.sites) o << " " << d; o << "\n"; }
  if (t.weight != 1.0) o << "weight " << t.weight << "\n";
  if (t.keep_charges) o << "charges keep\n";
  return o.str();
}

}  // namespace caps
