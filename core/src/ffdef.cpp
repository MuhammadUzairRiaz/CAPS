// CAPS force-field definitions: JSON format, moltemplate import, and parameter assignment.
#include "caps/ffdef.hpp"
#include "caps/charges.hpp"
#include "caps/resolution.hpp"
#include "caps/typing.hpp"
#include "caps/qeq.hpp"

#include <filesystem>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <map>
#include <set>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/grow.hpp"
#include "caps/json.hpp"

namespace caps {

namespace {

Json typing_json(const std::vector<TypingRule>& rules) {
  Json a = Json::array();
  for (const auto& r : rules) {
    Json o = Json::object();
    o["type"] = r.type;
    o["smarts"] = r.smarts;
    if (r.priority) o["priority"] = double(r.priority);
    if (!r.overrides.empty()) {
      Json ov = Json::array();
      for (const auto& x : r.overrides) ov.push_back(x);
      o["overrides"] = ov;
    }
    if (!r.description.empty()) o["description"] = r.description;
    a.push_back(o);
  }
  return a;
}

std::vector<TypingRule> typing_from(const Json& a, const std::string& where) {
  std::vector<TypingRule> out;
  for (const auto& o : a.items()) {
    TypingRule r;
    r.type = o.text("type");
    r.smarts = o.text("smarts");
    if (r.type.empty() || r.smarts.empty()) throw FFError(where + ": a typing rule needs \"type\" and \"smarts\"");
    r.description = o.text("description");
    if (o.has("priority")) r.priority = int(o["priority"].number());
    if (o.has("overrides"))
      for (const auto& x : o["overrides"].items()) r.overrides.push_back(x.str());
    for (auto [key, vec] : {std::make_pair("requires", &r.needs_elements), std::make_pair("excludes", &r.no_elements)})
      if (o.has(key))
        for (const auto& x : o[key].items()) {
          const int z = element_from_symbol(x.str());
          if (z <= 0) throw FFError(where + ": " + key + ": unknown element " + x.str());
          vec->push_back(z);
        }
    r.atom_name = o.text("atom_name");
    out.push_back(std::move(r));
  }
  return out;
}

}  // namespace

void load_typing(FFDef& ff, const std::string& path) {
  {
    std::error_code ec;
    const std::string key = std::filesystem::weakly_canonical(path, ec).string();
    if (std::find(ff.typing_files.begin(), ff.typing_files.end(), key) != ff.typing_files.end()) return;
    ff.typing_files.push_back(key);
  }
  std::ifstream in(path);
  if (!in) throw FFError("cannot open typing rules " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  Json j;
  try { j = Json::parse(ss.str()); } catch (const JsonError& e) { throw FFError(path + ": " + e.what()); }
  if (j.text("format") != "caps-typing") throw FFError(path + ": not a CAPS typing-rules file");
  auto rules = typing_from(j["rules"], path);
  std::set<std::string> known;
  for (const auto& t : ff.types) known.insert(t.name);
  // A rule may name a type by its short label: an alias, or the leading label of moltemplate's encoded names
  // ("135" for 135_bCT_aCT_dCT_iCT in OPLS-AA, "c4" for c4~pc4~bc4~… in COMPASS). Used only when unambiguous.
  {
    std::map<std::string, std::string> shortname;
    std::set<std::string> ambiguous;
    auto add = [&](const std::string& k, const std::string& full) {
      if (k.empty() || k == full || known.count(k)) return;
      auto [it, fresh] = shortname.emplace(k, full);
      if (!fresh && it->second != full) ambiguous.insert(k);
    };
    for (const auto& t : ff.types) {
      for (const auto& a : t.aliases) add(a, t.name);
      const size_t b = t.name.find("_b"), w = t.name.find('~');
      if (b != std::string::npos && b > 0) add(t.name.substr(0, b), t.name);
      if (w != std::string::npos && w > 0) add(t.name.substr(0, w), t.name);
    }
    auto resolve = [&](std::string& n) {
      if (known.count(n)) return;
      auto it = shortname.find(n);
      if (it != shortname.end() && !ambiguous.count(n)) n = it->second;
    };
    for (auto& r : rules) {
      resolve(r.type);
      for (auto& o : r.overrides) resolve(o);
    }
  }
  std::string unknown;
  for (const auto& r : rules)
    if (!known.count(r.type) && unknown.find(" " + r.type + ",") == std::string::npos) unknown += " " + r.type + ",";
  // "unknown_types": "untyped" — rules shared with a larger version of the force field (antechamber's ordered GAFF table
  // on moltemplate's GAFF): a rule for a type this file lacks stays in its place, so it still stops later, more general
  // rules from taking those atoms, and the atoms it matches are reported untyped
  if (!unknown.empty() && j.text("unknown_types") != "untyped")
    throw FFError(path + ": rules for types not in " + ff.name + ":" + unknown.substr(0, unknown.size() - 1));
  ff.typing.insert(ff.typing.end(), rules.begin(), rules.end());
  ff.typing_ordered = ff.typing_ordered || (j.has("ordered") && j["ordered"].boolean());
  ff.typing_unknown_untyped = ff.typing_unknown_untyped || j.text("unknown_types") == "untyped";
  ff.united_atom = ff.united_atom || (j.has("united_atom") && j["united_atom"].boolean());
  ff.keep_defined_bonds = ff.keep_defined_bonds || j.text("bonds") == "defined";
  if (j.has("shells") && j["shells"].is_object())
    for (const auto& [core, shell] : j["shells"].members()) ff.shells[core] = shell.str();
  if (j.has("variants") && j["variants"].is_object()) {
    std::set<std::string> have;
    for (const auto& ty : ff.types) have.insert(ty.name);
    for (const auto& [base, list] : j["variants"].members()) {
      std::vector<std::string> v;
      for (const auto& x : list.items())
        if (have.count(x.str())) v.push_back(x.str());
      if (!v.empty() && have.count(base)) ff.type_variants[base] = v;
    }
  }
  if (j.has("analogies") && j["analogies"].is_object()) {
    for (const auto& [ty, list] : j["analogies"].members())
      for (const auto& x : list.items()) ff.analogies[ty].push_back(x.str());
    ff.analogy_source = j.text("analogy_source", "parameters of analogous types");
  }
  ff.bond_k_per_order = j.num("bond_k_per_order", ff.bond_k_per_order);
  ff.bond_conjugated_single = j.num("conjugated_single_order", ff.bond_conjugated_single);
  if (j.text("pair_mode") == "double_same") ff.typing_pairs_double_same = true;
  else if (!j.text("pair_mode").empty() && j.text("pair_mode") != "double_differs")
    throw FFError(path + ": pair_mode is double_differs (GAFF) or double_same (CGenFF)");
  if (j.has("pairs"))
    for (const auto& pr : j["pairs"].items()) {
      if (pr.size() != 2) throw FFError(path + ": each entry of \"pairs\" names two types");
      for (int k = 0; k < 2; ++k)
        if (!known.count(pr[k].str()) && j.text("unknown_types") != "untyped") throw FFError(path + ": pair type " + pr[k].str() + " not in " + ff.name);
      ff.typing_pairs.push_back({pr[0].str(), pr[1].str()});
    }
  ff.typing_source = path;
}

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
}  // namespace

const FFType* FFDef::type(const std::string& n) const {
  for (const auto& t : types)
    if (t.name == n) return &t;
  return nullptr;
}

// '*' any run, '?' one character, '\' makes the next character literal (type names such as PCFF's "h*").
bool glob_match(const std::string& pattern, const std::string& text) {
  size_t p = 0, t = 0, star = std::string::npos, mark = 0;
  while (t < text.size()) {
    if (p + 1 < pattern.size() && pattern[p] == '\\' && pattern[p + 1] == text[t]) { p += 2; ++t; }
    else if (p < pattern.size() && pattern[p] != '\\' && pattern[p] != '*' && (pattern[p] == '?' || pattern[p] == text[t])) { ++p; ++t; }
    else if (p < pattern.size() && pattern[p] == '*') { star = p++; mark = t; }
    else if (star != std::string::npos) { p = star + 1; t = ++mark; }
    else return false;
  }
  while (p < pattern.size() && pattern[p] == '*') ++p;
  return p == pattern.size();
}

std::string glob_escape(const std::string& name) {
  std::string r;
  for (char c : name) {
    if (c == '*' || c == '?' || c == '\\') r += '\\';
    r += c;
  }
  return r;
}

// ---------------------------------------------------------------------------------------------------------------
// JSON

namespace {

Json rules_json(const std::vector<FFRule>& rules) {
  Json a = Json::array();
  for (const auto& r : rules) {
    Json o = Json::object();
    if (!r.name.empty()) o["name"] = r.name;
    Json m = Json::array();
    for (const auto& x : r.match) m.push_back(x);
    o["match"] = m;
    if (!r.style.empty()) o["style"] = r.style;
    Json p = Json::array();
    for (double v : r.params) p.push_back(v);
    o["params"] = p;
    for (const auto& [g, v] : r.cross) {
      Json q = Json::array();
      for (double x : v) q.push_back(x);
      o[g] = q;
    }
    if (!r.comment.empty()) o["comment"] = r.comment;
    a.push_back(o);
  }
  return a;
}

std::vector<FFRule> rules_from(const Json& a) {
  std::vector<FFRule> out;
  for (const auto& o : a.items()) {
    FFRule r;
    r.name = o.text("name");
    for (const auto& m : o["match"].items()) r.match.push_back(m.str());
    r.style = o.text("style");
    for (const auto& v : o["params"].items()) r.params.push_back(v.number());
    // named parameter groups: class II cross terms, separate 1-4 van der Waals, alternative forms (DREIDING)
    static const std::set<std::string> plain = {"name", "match", "style", "params", "comment"};
    for (const auto& [g, v] : o.members())
      if (!plain.count(g) && v.is_array())
        for (const auto& x : v.items()) r.cross[g].push_back(x.number());
    r.comment = o.text("comment");
    out.push_back(std::move(r));
  }
  return out;
}

}  // namespace

void save_forcefield(const FFDef& ff, const std::string& path) {
  Json j = Json::object();
  j["format"] = "caps-forcefield";
  j["format_version"] = 1;
  j["name"] = ff.name;
  if (!ff.version.empty()) j["version"] = ff.version;
  j["source"] = ff.source;
  Json refs = Json::array();
  for (const auto& r : ff.references) refs.push_back(r);
  j["references"] = refs;
  j["units"] = ff.units;
  Json st = Json::object();
  st["pair"] = ff.pair_style;
  st["bond"] = ff.bond_style;
  st["angle"] = ff.angle_style;
  st["dihedral"] = ff.dihedral_style;
  st["improper"] = ff.improper_style;
  j["styles"] = st;
  j["mixing"] = ff.mixing;
  Json sl = Json::array(), sc = Json::array();
  for (int k = 0; k < 3; ++k) { sl.push_back(ff.special_lj[k]); sc.push_back(ff.special_coul[k]); }
  j["special_lj"] = sl;
  j["special_coul"] = sc;
  j["cutoff"] = ff.cutoff;
  j["improper_order"] = ff.improper_order;
  j["equivalence"] = ff.equivalence;
  if (ff.improper_matched_order) j["improper_matched_order"] = true;
  if (ff.improper_reversible) j["improper_reversible"] = true;
  if (ff.improper_all_explicit) j["improper_all_explicit"] = true;
  if (ff.improper_max_neighbours) j["improper_max_neighbours"] = ff.improper_max_neighbours;
  if (ff.wildcard_torsion_scaling != "none") j["wildcard_torsion_scaling"] = ff.wildcard_torsion_scaling;
  Json types = Json::array();
  for (const auto& t : ff.types) {
    Json o = Json::object();
    o["name"] = t.name;
    o["element"] = t.element > 0 ? element(t.element).symbol : "";
    o["mass"] = t.mass;
    if (!std::isnan(t.charge)) o["charge"] = t.charge;
    if (!t.description.empty()) o["description"] = t.description;
    if (!t.smarts.empty()) o["smarts"] = t.smarts;
    if (t.priority) o["priority"] = double(t.priority);
    if (!t.overrides.empty()) {
      Json ov = Json::array();
      for (const auto& x : t.overrides) ov.push_back(x);
      o["overrides"] = ov;
    }
    if (!t.source.empty()) o["source"] = t.source;
    if (!t.equiv.empty()) {
      Json e = Json::object();
      for (const auto& [k, v] : t.equiv) e[k] = v;
      o["equivalence"] = e;
    }
    if (!t.aliases.empty()) {
      Json a = Json::array();
      for (const auto& x : t.aliases) a.push_back(x);
      o["aliases"] = a;
    }
    types.push_back(o);
  }
  j["atom_types"] = types;
  j["pairs"] = rules_json(ff.pairs);
  j["bonds"] = rules_json(ff.bonds);
  j["angles"] = rules_json(ff.angles);
  j["dihedrals"] = rules_json(ff.dihedrals);
  j["impropers"] = rules_json(ff.impropers);
  if (!ff.bond_increments.empty()) j["bond_increments"] = rules_json(ff.bond_increments);
  if (!ff.auto_bonds.empty()) j["auto_bonds"] = rules_json(ff.auto_bonds);
  if (!ff.auto_angles.empty()) j["auto_angles"] = rules_json(ff.auto_angles);
  if (!ff.auto_dihedrals.empty()) j["auto_dihedrals"] = rules_json(ff.auto_dihedrals);
  if (!ff.auto_impropers.empty()) j["auto_impropers"] = rules_json(ff.auto_impropers);
  if (!ff.cross_rules.empty()) {
    Json c = Json::object();
    for (const auto& [g, rules] : ff.cross_rules) c[g] = rules_json(rules);
    j["cross_rules"] = c;
  }
  if (!ff.oop_scheme.empty()) j["oop_scheme"] = ff.oop_scheme;
  if (!ff.typing.empty()) j["typing"] = typing_json(ff.typing);
  if (ff.typing_ordered) j["typing_ordered"] = true;
  if (ff.typing_pairs_double_same) j["typing_pair_mode"] = "double_same";
  if (!ff.typing_pairs.empty()) {
    Json a = Json::array();
    for (const auto& [x, y] : ff.typing_pairs) {
      Json pr = Json::array();
      pr.push_back(x);
      pr.push_back(y);
      a.push_back(pr);
    }
    j["typing_pairs"] = a;
  }
  Json notes = Json::array();
  for (const auto& n : ff.notes) notes.push_back(n);
  j["notes"] = notes;
  std::ofstream out(path);
  if (!out) throw FFError("cannot write " + path);
  out << j.dump(1) << "\n";
}

FFDef load_forcefield(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw FFError("cannot open " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  Json j;
  try { j = Json::parse(ss.str()); } catch (const JsonError& e) { throw FFError(path + ": " + e.what()); }
  if (j.text("format") != "caps-forcefield") throw FFError(path + ": not a CAPS force-field file");
  FFDef ff;
  ff.name = j.text("name");
  ff.version = j.text("version");
  ff.source = j.text("source");
  if (j.has("references"))
    for (const auto& r : j["references"].items()) ff.references.push_back(r.str());
  ff.units = j.text("units", "real");
  if (j.has("styles")) {
    const Json& st = j["styles"];
    ff.pair_style = st.text("pair", ff.pair_style);
    ff.bond_style = st.text("bond", ff.bond_style);
    ff.angle_style = st.text("angle", ff.angle_style);
    ff.dihedral_style = st.text("dihedral", ff.dihedral_style);
    ff.improper_style = st.text("improper", ff.improper_style);
  }
  ff.mixing = j.text("mixing", ff.mixing);
  if (j.has("special_lj"))
    for (int k = 0; k < 3; ++k) ff.special_lj[k] = j["special_lj"][k].number();
  if (j.has("special_coul"))
    for (int k = 0; k < 3; ++k) ff.special_coul[k] = j["special_coul"][k].number();
  ff.cutoff = j.num("cutoff", ff.cutoff);
  ff.improper_order = j.text("improper_order", ff.improper_order);
  ff.equivalence = j.text("equivalence", ff.equivalence);
  ff.improper_matched_order = j.has("improper_matched_order") && j["improper_matched_order"].boolean();
  ff.improper_reversible = j.has("improper_reversible") && j["improper_reversible"].boolean();
  ff.improper_all_explicit = j.has("improper_all_explicit") && j["improper_all_explicit"].boolean();
  ff.improper_max_neighbours = int(j.num("improper_max_neighbours", 0));
  ff.wildcard_torsion_scaling = j.text("wildcard_torsion_scaling", ff.wildcard_torsion_scaling);
  for (const auto& o : j["atom_types"].items()) {
    FFType t;
    t.name = o["name"].str();
    t.element = element_from_symbol(o.text("element"));
    t.mass = o.num("mass", 0);
    if (o.has("charge")) t.charge = o["charge"].number();
    t.description = o.text("description");
    t.smarts = o.text("smarts");
    if (o.has("priority")) t.priority = int(o["priority"].number());
    if (o.has("overrides"))
      for (const auto& x : o["overrides"].items()) t.overrides.push_back(x.str());
    t.source = o.text("source");
    if (o.has("equivalence"))
      for (const auto& [k, v] : o["equivalence"].members()) t.equiv[k] = v.str();
    if (o.has("aliases"))
      for (const auto& x : o["aliases"].items()) t.aliases.push_back(x.str());
    ff.types.push_back(std::move(t));
  }
  for (auto [key, vec] : {std::make_pair("pairs", &ff.pairs), std::make_pair("bonds", &ff.bonds), std::make_pair("angles", &ff.angles),
                          std::make_pair("dihedrals", &ff.dihedrals), std::make_pair("impropers", &ff.impropers),
                          std::make_pair("bond_increments", &ff.bond_increments), std::make_pair("auto_bonds", &ff.auto_bonds),
                          std::make_pair("auto_angles", &ff.auto_angles), std::make_pair("auto_dihedrals", &ff.auto_dihedrals),
                          std::make_pair("auto_impropers", &ff.auto_impropers)})
    if (j.has(key)) *vec = rules_from(j[key]);
  if (j.has("cross_rules"))
    for (const auto& [g, rules] : j["cross_rules"].members()) ff.cross_rules[g] = rules_from(rules);
  ff.oop_scheme = j.text("oop_scheme");
  if (j.has("notes"))
    for (const auto& n : j["notes"].items()) ff.notes.push_back(n.str());
  // an overlay naming the force field it extends (L-OPLS on OPLS-AA): the base first, this file's types and terms on
  // top (merge_forcefield; later terms win), one complete force field under this file's name
  if (j.has("extends")) {
    const std::filesystem::path bp(j["extends"].str());
    FFDef base = load_forcefield((bp.is_absolute() ? bp : std::filesystem::path(path).parent_path() / bp).lexically_normal().string());
    if (j.has("typing")) {   // this file's own rules replace the base's
      base.typing.clear();
      base.typing_pairs.clear();
      base.analogies.clear();
      base.typing_files.clear();
    }
    merge_forcefield(base, ff);
    base.name = ff.name;
    base.version = ff.version;
    base.source = ff.source + " on " + base.source;
    base.references.insert(base.references.end(), ff.references.begin(), ff.references.end());
    ff = std::move(base);
  }
  if (j.has("typing")) {
    if (j["typing"].is_string()) {
      // relative to the force-field file (either separator: Windows paths use backslashes)
      const std::filesystem::path f(j["typing"].str());
      load_typing(ff, (f.empty() || f.is_absolute() ? f : std::filesystem::path(path).parent_path() / f).lexically_normal().string());
    } else {
      ff.typing = typing_from(j["typing"], path);
      ff.typing_source = path;
      ff.typing_ordered = j.has("typing_ordered") && j["typing_ordered"].boolean();
      ff.typing_pairs_double_same = j.text("typing_pair_mode") == "double_same";
      if (j.has("typing_pairs"))
        for (const auto& pr : j["typing_pairs"].items()) ff.typing_pairs.push_back({pr[0].str(), pr[1].str()});
    }
  }
  for (const auto& t : ff.types)
    if (!t.smarts.empty()) ff.typing.push_back({t.name, t.smarts, t.description, t.overrides, t.priority});
  return ff;
}

void merge_forcefield(FFDef& base, const FFDef& ov) {
  for (const auto& t : ov.types) {
    auto it = std::find_if(base.types.begin(), base.types.end(), [&](const FFType& x) { return x.name == t.name; });
    if (it != base.types.end()) *it = t;
    else base.types.push_back(t);
  }
  auto append = [](std::vector<FFRule>& a, const std::vector<FFRule>& b) { a.insert(a.end(), b.begin(), b.end()); };
  append(base.pairs, ov.pairs);
  append(base.bonds, ov.bonds);
  append(base.angles, ov.angles);
  append(base.dihedrals, ov.dihedrals);
  append(base.impropers, ov.impropers);
  append(base.bond_increments, ov.bond_increments);
  append(base.auto_bonds, ov.auto_bonds);
  append(base.auto_angles, ov.auto_angles);
  append(base.auto_dihedrals, ov.auto_dihedrals);
  append(base.auto_impropers, ov.auto_impropers);
  for (const auto& [g, rules] : ov.cross_rules) append(base.cross_rules[g], rules);
  base.typing.insert(base.typing.begin(), ov.typing.begin(), ov.typing.end());   // the overlay's rules come first
  base.typing_pairs.insert(base.typing_pairs.end(), ov.typing_pairs.begin(), ov.typing_pairs.end());
  base.notes.push_back("merged overlay: " + (ov.name.empty() ? std::string("(unnamed)") : ov.name));
}

// ---------------------------------------------------------------------------------------------------------------
// moltemplate import

namespace {

std::string strip_prefix(const std::string& w, const char* pre) {
  const std::string p = pre;
  return w.rfind(p, 0) == 0 ? w.substr(p.size()) : w;
}

bool numeric(const std::string& w) {
  if (w.empty()) return false;
  char* end = nullptr;
  std::strtod(w.c_str(), &end);
  return end && *end == 0;
}

int lt_element_from_mass(double m) {
  if (m < 0.5) return 0;   // massless sites (virtual sites, TIP4P M)
  return element_from_mass(m, 0.1);
}

}  // namespace

FFDef import_moltemplate(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw FFError("cannot open " + path);
  FFDef ff;
  ff.source = "moltemplate: " + path;
  std::string line, section, top;
  int depth = 0;
  std::map<std::string, std::string> alias;              // short → full
  std::map<std::string, double> masses, charges;
  std::map<std::string, std::string> desc, charge_desc;
  std::map<std::string, FFRule> coeff[5];                 // bond, angle, dihedral, improper; [4] unused
  std::vector<std::pair<std::string, std::vector<std::string>>> bytype[4];
  std::vector<std::pair<std::string, std::vector<std::string>>> pair_lines;   // "a b" → params (with style)
  std::map<std::string, std::string> pair_comment;
  std::string improper_plugin;
  std::vector<FFRule> increments;
  int lineno = 0;
  auto bad = [&](const std::string& why) { return FFError(path + ":" + std::to_string(lineno) + ": " + why); };

  while (std::getline(in, line)) {
    ++lineno;
    std::string comment;
    if (auto h = line.find('#'); h != std::string::npos) {
      comment = line.substr(h + 1);
      line = line.substr(0, h);
    }
    // type descriptions in header comments, e.g. "#  @atom:c3a  C  "aromatic carbon" (ver=1.0, ref=1)"
    if (line.find_first_not_of(" \t") == std::string::npos && comment.find("@atom:") != std::string::npos && depth <= 1) {
      std::istringstream cs(comment);
      std::string w;
      cs >> w;
      if (w.rfind("@atom:", 0) == 0) {
        std::string rest;
        std::getline(cs, rest);
        const auto a = rest.find_first_not_of(" \t");
        if (a != std::string::npos && !desc.count(strip_prefix(w, "@atom:"))) desc[strip_prefix(w, "@atom:")] = rest.substr(a);
      }
    }
    std::istringstream ls(line);
    std::vector<std::string> t;
    for (std::string w; ls >> w;) t.push_back(w);
    if (t.empty()) continue;
    auto trim = [](std::string s) {
      const auto a = s.find_first_not_of(" \t\r");
      const auto b = s.find_last_not_of(" \t\r");
      return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
    };
    const std::string joined = trim(line);
    if (joined.rfind("replace{", 0) == 0 || joined.rfind("replace {", 0) == 0) {
      std::string body = joined.substr(joined.find('{') + 1);
      body = body.substr(0, body.find('}'));
      std::istringstream bs(body);
      std::string a, b;
      bs >> a >> b;
      alias[strip_prefix(a, "@atom:")] = strip_prefix(b, "@atom:");
      continue;
    }
    if (joined.find("write_once(") == 0 || joined.find("write(") == 0) {
      const auto q1 = joined.find('"'), q2 = joined.find('"', q1 + 1);
      section = joined.substr(q1 + 1, q2 - q1 - 1);
      if (section.find("Impropers By Type") != std::string::npos) {
        const auto p1 = section.find('('), p2 = section.find(')');
        improper_plugin = p1 != std::string::npos ? section.substr(p1 + 1, p2 - p1 - 1) : "";
      }
      depth = 2;
      continue;
    }
    if (depth == 0 && t.back() == "{") {
      top = t[0];
      ff.name = top;
      depth = 1;
      continue;
    }
    if (t[0] == "}") {
      if (depth == 2) { section.clear(); depth = 1; }
      else depth = 0;
      continue;
    }
    if (section.empty()) continue;
    const std::string cmt = trim(comment);
    if (section == "In Init") {
      if (t[0] == "pair_style") {
        size_t k = 1;
        if (t.size() > 1 && t[1] == "hybrid") ++k;
        if (k < t.size()) ff.pair_style = t[k];
        for (size_t i = k + 1; i < t.size(); ++i)
          if (numeric(t[i])) ff.cutoff = std::stod(t[i]);   // last numeric is the (outer) cut-off
      } else if (t[0] == "bond_style" || t[0] == "angle_style" || t[0] == "dihedral_style" || t[0] == "improper_style") {
        const std::string st = t.size() > 2 && t[1] == "hybrid" ? t[2] : t.size() > 1 ? t[1] : "";
        (t[0] == "bond_style" ? ff.bond_style : t[0] == "angle_style" ? ff.angle_style : t[0] == "dihedral_style" ? ff.dihedral_style
                                                                                                   : ff.improper_style) = st;
      } else if (t[0] == "special_bonds") {
        for (size_t i = 1; i < t.size(); ++i) {
          if (t[i] == "amber") { ff.special_lj[2] = 0.5; ff.special_coul[2] = 0.8333333333; }
          else if (t[i] == "charmm") { ff.special_lj[2] = 0; ff.special_coul[2] = 0; }
          else if (t[i] == "dreiding") { ff.special_lj[2] = 1; ff.special_coul[2] = 1; }
          else if ((t[i] == "lj/coul" || t[i] == "lj" || t[i] == "coul") && i + 3 < t.size()) {
            for (int k = 0; k < 3; ++k) {
              const double v = std::stod(t[i + 1 + k]);
              if (t[i] != "coul") ff.special_lj[k] = v;
              if (t[i] != "lj") ff.special_coul[k] = v;
            }
            i += 3;
          }
        }
      } else if (t[0] == "pair_modify") {
        for (size_t i = 1; i + 1 < t.size(); ++i)
          if (t[i] == "mix") ff.mixing = t[i + 1];
      }
    } else if (section == "Data Masses") {
      if (t.size() >= 2 && t[0].rfind("@atom:", 0) == 0) {
        masses[strip_prefix(t[0], "@atom:")] = std::stod(t[1]);
        if (!cmt.empty() && !desc.count(strip_prefix(t[0], "@atom:"))) desc[strip_prefix(t[0], "@atom:")] = cmt;
      }
    } else if (section == "In Charges") {
      // set type @atom:81 charge 0.265 # C - C2 | CH2 IN ETHANOL
      if (t.size() >= 5 && t[0] == "set" && t[1] == "type") {
        charges[strip_prefix(t[2], "@atom:")] = std::stod(t[4]);
        if (!cmt.empty()) charge_desc[strip_prefix(t[2], "@atom:")] = cmt;
      }
    } else if (section == "In Settings") {
      const std::string& cmd = t[0];
      const int kind = cmd == "bond_coeff" ? 0 : cmd == "angle_coeff" ? 1 : cmd == "dihedral_coeff" ? 2 : cmd == "improper_coeff" ? 3 : -1;
      if (cmd == "pair_coeff" && t.size() >= 4) {
        std::vector<std::string> v(t.begin() + 3, t.end());
        pair_lines.push_back({strip_prefix(t[1], "@atom:") + " " + strip_prefix(t[2], "@atom:"), v});
        pair_comment[pair_lines.back().first] = cmt;
      } else if (kind >= 0 && t.size() >= 3) {
        FFRule r;
        const char* pre = kind == 0 ? "@bond:" : kind == 1 ? "@angle:" : kind == 2 ? "@dihedral:" : "@improper:";
        r.name = strip_prefix(t[1], pre);
        size_t k = 2;
        if (!numeric(t[2])) { r.style = t[2]; k = 3; }
        for (; k < t.size(); ++k) {
          if (!numeric(t[k])) throw bad("expected a number in '" + joined + "'");
          r.params.push_back(std::stod(t[k]));
        }
        r.comment = cmt;
        // class II cross terms arrive as extra coeff lines for the same type: "angle_coeff @angle:X bb M r1 r2"
        static const std::set<std::string> cross = {"bb", "ba", "mbt", "ebt", "at", "aat", "bb13", "aa"};
        if (cross.count(r.style)) {
          FFRule& base = coeff[kind][r.name];
          base.name = r.name;
          base.cross[r.style] = r.params;
          if (base.comment.empty()) base.comment = cmt;
        } else {
          FFRule& base = coeff[kind][r.name];
          r.cross = base.cross;
          base = r;
        }
      }
    } else if (section == "Data Charge By Bond") {
      // @atom:a @atom:b δa δb
      if (t.size() >= 4 && numeric(t[2]) && numeric(t[3])) {
        FFRule r;
        r.match = {strip_prefix(t[0], "@atom:"), strip_prefix(t[1], "@atom:")};
        r.params = {std::stod(t[2]), std::stod(t[3])};
        r.comment = cmt;
        increments.push_back(r);
      }
    } else if (section.find("By Type") != std::string::npos) {
      const int kind = section.find("Bonds") == 5 ? 0 : section.find("Angles") == 5 ? 1 : section.find("Dihedrals") == 5 ? 2
                       : section.find("Impropers") == 5 ? 3 : -1;
      if (kind < 0 || t.size() < 3) continue;
      const char* pre = kind == 0 ? "@bond:" : kind == 1 ? "@angle:" : kind == 2 ? "@dihedral:" : "@improper:";
      std::vector<std::string> pats;
      for (size_t i = 1; i < t.size(); ++i) {
        if (t[i].rfind("@bond:", 0) == 0) continue;   // bond-type constraints are not used in these files
        pats.push_back(strip_prefix(t[i], "@atom:"));
      }
      bytype[kind].push_back({strip_prefix(t[0], pre), pats});
    }
  }

  auto full = [&](const std::string& n) { auto it = alias.find(n); return it == alias.end() ? n : it->second; };
  // Types: every name with a mass, under its full name.
  std::set<std::string> seen;
  for (const auto& [n, m] : masses) {
    FFType t;
    t.name = full(n);
    if (!seen.insert(t.name).second) continue;
    t.mass = m;
    t.element = lt_element_from_mass(m);
    if (auto it = charges.find(n); it != charges.end()) t.charge = it->second;
    else if (auto it2 = charges.find(t.name); it2 != charges.end()) t.charge = it2->second;
    std::string d = desc.count(n) ? desc[n] : desc.count(t.name) ? desc[t.name] : "";
    if (charge_desc.count(n)) d = charge_desc[n] + (d.empty() ? "" : " · " + d);
    if (t.name != n) d = "moltemplate @atom:" + n + (d.empty() ? "" : " · " + d);
    t.description = d;
    ff.types.push_back(t);
  }
  // Pair coefficients: self and explicit pairs (styles kept when hybrid).
  for (const auto& [ab, v] : pair_lines) {
    std::istringstream is(ab);
    std::string a, b;
    is >> a >> b;
    FFRule r;
    a = full(a);
    b = full(b);
    r.match = a == b ? std::vector<std::string>{a} : std::vector<std::string>{a, b};
    size_t k = 0;
    if (!v.empty() && !numeric(v[0])) { r.style = v[0]; k = 1; }
    for (; k < v.size(); ++k)
      if (numeric(v[k])) r.params.push_back(std::stod(v[k]));
    r.comment = pair_comment[ab];
    ff.pairs.push_back(r);
  }
  // Bonded rules: the By Type patterns in file order, each with its coefficients.
  const char* kinds[] = {"bond", "angle", "dihedral", "improper"};
  std::vector<FFRule>* out[] = {&ff.bonds, &ff.angles, &ff.dihedrals, &ff.impropers};
  for (int k = 0; k < 4; ++k) {
    int missing = 0;
    for (const auto& [name, pats] : bytype[k]) {
      auto it = coeff[k].find(name);
      if (it == coeff[k].end()) { ++missing; continue; }
      FFRule r = it->second;
      r.match.clear();
      for (const auto& p : pats) r.match.push_back(full(p));
      out[k]->push_back(r);
    }
    if (missing) ff.notes.push_back(std::to_string(missing) + " " + kinds[k] + " rules without coefficients were skipped");
  }
  for (auto& r : increments)
    for (auto& m : r.match) m = full(m);
  ff.bond_increments = increments;
  if (improper_plugin == "cenIsortJKL.py" || improper_plugin == "cenIsortJKL") ff.improper_order = "center1_sorted";
  else if (improper_plugin == "cenJsortIKL.py" || improper_plugin == "cenJsortIKL") ff.improper_order = "center2_sorted";
  else if (improper_plugin == "gaff_imp.py" || improper_plugin == "gaff_imp" || improper_plugin.empty()) ff.improper_order = "center3_sorted";
  else ff.notes.push_back("improper plugin '" + improper_plugin + "' is treated as centre-third, sorted");
  ff.notes.push_back("converted from moltemplate by CAPS; rules keep moltemplate's order (the last matching rule wins)");
  if (ff.types.empty()) throw FFError(path + ": no atom types (Data Masses) found");
  return ff;
}

// ---------------------------------------------------------------------------------------------------------------
// Parameter assignment

namespace {

bool match_fw(const std::vector<std::string>& pat, const std::vector<const std::string*>& ty) {
  if (pat.size() != ty.size()) return false;
  for (size_t i = 0; i < pat.size(); ++i)
    if (!glob_match(pat[i], *ty[i])) return false;
  return true;
}

// Last rule whose patterns match the types forwards or backwards (*reversed tells which; forwards is preferred).
// wild: 0 any rule, 1 only rules without wildcards, 2 only rules with a wildcard
const FFRule* last_match(const std::vector<FFRule>& rules, const std::vector<const std::string*>& ty, bool reversible,
                         bool* reversed = nullptr, int wild = 0) {
  std::vector<const std::string*> rev(ty.rbegin(), ty.rend());
  for (auto it = rules.rbegin(); it != rules.rend(); ++it) {
    if (wild) {
      const bool w = std::find(it->match.begin(), it->match.end(), std::string("*")) != it->match.end();
      if ((wild == 1) == w) continue;
    }
    if (match_fw(it->match, ty)) { if (reversed) *reversed = false; return &*it; }
    if (reversible && match_fw(it->match, rev)) { if (reversed) *reversed = true; return &*it; }
  }
  return nullptr;
}

// A class II cross-term group, zeros when the rule does not give it.
std::vector<double> group(const FFRule& r, const char* g, size_t n) {
  auto it = r.cross.find(g);
  std::vector<double> v = it == r.cross.end() ? std::vector<double>() : it->second;
  if (it != r.cross.end() && v.size() != n)
    throw FFError("class2 " + std::string(g) + " terms of " + r.name + " need " + std::to_string(n) + " coefficients");
  v.resize(n, 0.0);
  return v;
}


}  // namespace

// Bond-order variants (FFDef::type_variants): each conjugated system of atoms whose type has variants takes the
// combination that makes every bond's harmonic constant bond_k_per_order x its perceived order (aromatic bonds 1.5),
// found by depth-first search with the bonds to already-placed atoms checked at each step.
void refine_bond_order_variants(const System& s, const Perception& p, const FFDef& ff, std::vector<std::string>& types, std::vector<std::string>& why) {
  if (ff.type_variants.empty() || ff.bond_k_per_order <= 0) return;
  std::map<std::string, const FFType*> byname;
  for (const auto& ty : ff.types) byname[ty.name] = &ty;
  auto bond_name = [&](const std::string& t) {
    auto it = byname.find(t);
    if (it == byname.end()) return t;
    auto e = it->second->equiv.find("bond");
    return e == it->second->equiv.end() ? t : e->second;
  };
  std::map<std::pair<std::string, std::string>, double> kcache;
  auto k_of = [&](const std::string& a, const std::string& b) {
    auto key = std::make_pair(a, b);
    auto it = kcache.find(key);
    if (it != kcache.end()) return it->second;
    const std::string A = bond_name(a), B = bond_name(b);
    const FFRule* r = last_match(ff.bonds, {&A, &B}, true);
    const double k = r && !r->params.empty() ? r->params[0] : -1;
    kcache[key] = k;
    kcache[{b, a}] = k;
    return k;
  };
  const size_t n = s.atoms.size();
  std::vector<std::vector<std::string>> cand(n);
  std::vector<char> var(n, 0);
  for (size_t i = 0; i < n; ++i) {
    cand[i] = {types[i]};
    auto it = ff.type_variants.find(types[i]);
    if (it != ff.type_variants.end()) {
      cand[i].insert(cand[i].end(), it->second.begin(), it->second.end());
      var[i] = 1;
    }
  }
  // a single bond between two atoms that each carry a double or aromatic bond is conjugated: some files give it 1.5
  std::vector<char> unsat(n, 0);
  for (size_t i = 0; i < n; ++i)
    for (size_t k = 0; k < p.nb[i].size(); ++k)
      if (p.arom_bond[i][k] || p.order[i][k] >= 2) unsat[i] = 1;
  auto expected = [&](uint32_t i, size_t k) {
    if (p.arom_bond[i][k]) return 1.5;
    const int o = p.order[i][k];
    if (o == 1 && unsat[i] && unsat[p.nb[i][k]]) return ff.bond_conjugated_single;
    return double(o);
  };
  auto fits = [&](const std::string& a, const std::string& b, double order) {
    const double k = k_of(a, b);
    return k >= 0 && std::fabs(k - ff.bond_k_per_order * order) <= 1e-6 * std::max(1.0, k);
  };
  std::vector<char> seen(n, 0);
  for (uint32_t start = 0; start < n; ++start) {
    if (!var[start] || seen[start]) continue;
    std::vector<uint32_t> comp{start};   // the conjugated system: variable atoms joined by bonds
    seen[start] = 1;
    for (size_t q = 0; q < comp.size(); ++q)
      for (uint32_t j : p.nb[comp[q]])
        if (var[j] && !seen[j]) { seen[j] = 1; comp.push_back(j); }
    std::vector<char> placed(n, 0);
    std::vector<std::string> cur = types;
    size_t steps = 0;
    std::function<bool(size_t)> dfs = [&](size_t q) -> bool {
      if (q == comp.size()) return true;
      if (++steps > 200000) return false;
      const uint32_t a = comp[q];
      for (const auto& c : cand[a]) {
        cur[a] = c;
        bool ok = true;
        for (size_t k = 0; k < p.nb[a].size() && ok; ++k) {
          const uint32_t b = p.nb[a][k];
          if (var[b] && !placed[b]) continue;   // checked when b is placed
          ok = fits(cur[a], cur[b], expected(a, k));
        }
        if (!ok) continue;
        placed[a] = 1;
        if (dfs(q + 1)) return true;
        placed[a] = 0;
      }
      cur[a] = types[a];
      return false;
    };
    const bool solved = dfs(0);
    for (uint32_t a : comp) {
      if (solved && cur[a] != types[a]) {
        types[a] = cur[a];
        if (a < why.size()) why[a] += " · variant " + cur[a] + " for its bond orders";
      } else if (!solved && a < why.size()) {
        why[a] += " · no variant combination gives every bond its order's constant";
      }
    }
  }
}

ForceField parameterize(const System& s, const FFDef& def, const std::vector<std::string>& types_in, const std::string& charges,
                        ParamReport* rep_out, bool allow_missing) {
  ParamReport rep;
  const size_t n = s.atoms.size();
  if (types_in.size() != n) throw FFError("one force-field type per atom is needed");
  ForceField ff;
  ff.name = def.name;
  ff.mixing = def.mixing;
  ff.lj14 = def.special_lj[2];
  ff.coul14 = def.special_coul[2];
  if (def.special_lj[0] != 0 || def.special_lj[1] != 0 || def.special_coul[0] != 0 || def.special_coul[1] != 0)
    throw FFError(def.name + ": 1-2 or 1-3 non-bonded scaling other than 0 is not supported yet");
  // resolve names: full names, or a moltemplate short name recorded in the description
  std::map<std::string, const FFType*> byname;
  for (const auto& t : def.types) byname[t.name] = &t;
  for (const auto& t : def.types)
    for (const auto& a : t.aliases) byname.emplace(a, &t);
  for (const auto& t : def.types)
    if (t.description.rfind("moltemplate @atom:", 0) == 0) {
      std::istringstream is(t.description.substr(18));
      std::string shortn;
      is >> shortn;
      byname.emplace(shortn, &t);
    }
  std::vector<std::string> T(n);
  std::vector<const FFType*> FT(n);
  std::set<std::string> unknown;
  for (size_t i = 0; i < n; ++i) {
    auto it = byname.find(types_in[i]);
    if (it == byname.end()) { unknown.insert(types_in[i]); continue; }
    FT[i] = it->second;
    T[i] = it->second->name;
  }
  if (!unknown.empty()) {
    std::string l;
    for (const auto& u : unknown) l += (l.empty() ? "" : ", ") + u;
    throw FFError("types not in " + def.name + ": " + l);
  }
  ff.atom_type = T;
  ff.why.assign(n, "type given");
  // type indices in first-use order
  std::map<std::string, int> tix;
  for (size_t i = 0; i < n; ++i) {
    auto [it, fresh] = tix.emplace(T[i], int(ff.type_names.size()));
    if (fresh) ff.type_names.push_back(T[i]);
    ff.type_index.push_back(it->second);
    // a type without a mass: the element's, plus the hydrogens a united-atom site carries (a CH4 site weighs 16.04)
    const int uh = std::max(0, united_atom_hydrogens(s.atoms[i].name, s.atoms[i].element));
    ff.mass.push_back(FT[i]->mass > 0 ? FT[i]->mass : element(s.atoms[i].element).mass + uh * element(1).mass);
  }
  // Names used for each kind of lookup: the type's own, or its equivalent for that kind.
  auto names_for = [&](const char* kind) {
    std::vector<std::string> v(n);
    for (size_t i = 0; i < n; ++i) {
      auto it = FT[i]->equiv.find(kind);
      v[i] = it == FT[i]->equiv.end() ? T[i] : it->second;
    }
    return v;
  };
  const auto Nb = names_for("bond"), Na = names_for("angle"), Nd = names_for("dihedral"), Ni = names_for("improper"),
             Nq = names_for("increment");
  // second-choice equivalents (DL_FIELD "bond2_" etc.): the first names where a type has none
  auto names2 = [&](const char* kind, const std::vector<std::string>& first) {
    std::vector<std::string> v = first;
    bool any = false;
    for (size_t i = 0; i < n; ++i)
      if (auto it = FT[i]->equiv.find(kind); it != FT[i]->equiv.end()) { v[i] = it->second; any = true; }
    return any ? v : std::vector<std::string>();
  };
  const auto Nb2 = names2("bond2", Nb), Na2 = names2("angle2", Na), Nd2 = names2("dihedral2", Nd), Ni2 = names2("improper2", Ni);
  const bool fallback = def.equivalence == "fallback";
  const auto Ab = names_for("auto_bond"), Aae = names_for("auto_angle_end"), Aaa = names_for("auto_angle_apex"),
             Ate = names_for("auto_torsion_end"), Atc = names_for("auto_torsion_center");
  int n_auto = 0;
  // automatic parameters, by position-dependent names (only when the explicit rules gave nothing)
  auto auto_lookup = [&](const std::vector<FFRule>& rules, std::vector<const std::string*> names, bool* rev) -> const FFRule* {
    if (rules.empty()) return nullptr;
    const FFRule* r = last_match(rules, names, true, rev);
    if (r) ++n_auto;
    return r;
  };
  // last matching rule over the equivalent names (msi2lmp "fallback": the own names first)
  auto lookup = [&](const std::vector<FFRule>& rules, const std::vector<std::string>& N, std::initializer_list<uint32_t> atoms,
                    bool* rev = nullptr) -> const FFRule* {
    std::vector<const std::string*> raw, eq;
    for (uint32_t a : atoms) { raw.push_back(&T[a]); eq.push_back(&N[a]); }
    if (fallback)
      if (const FFRule* r = last_match(rules, raw, true, rev)) return r;
    return last_match(rules, eq, true, rev);
  };
  // with second-choice names (DL_FIELD): exact rules by the first names, exact by the second, then wildcard rules
  auto lookup12 = [&](const std::vector<FFRule>& rules, const std::vector<std::string>& N1, const std::vector<std::string>& N2,
                      std::initializer_list<uint32_t> atoms, bool* rev = nullptr) -> const FFRule* {
    if (N2.empty()) return lookup(rules, N1, atoms, rev);
    std::vector<const std::string*> a1, a2;
    for (uint32_t a : atoms) { a1.push_back(&N1[a]); a2.push_back(&N2[a]); }
    for (int wild : {1, 2})
      for (const auto* names : {&a1, &a2})
        if (const FFRule* r = last_match(rules, *names, true, rev, wild)) return r;
    return nullptr;
  };
  auto shown = [&](const std::vector<std::string>& N, std::initializer_list<uint32_t> atoms) {
    std::string r;
    for (uint32_t a : atoms) r += (r.empty() ? "" : " ") + T[a] + (N[a] != T[a] ? "(" + N[a] + ")" : "");
    return r;
  };
  // By analogy (FFDef::analogies, parmchk2's approach): one atom's type replaced by an analogue, then two, most similar
  // first; the term is kept and listed as estimated with the types that stood in.
  auto analog = [&](const std::vector<FFRule>& rules, const std::vector<std::string>& N, std::initializer_list<uint32_t> atoms, bool* rev,
                    std::string* used) -> const FFRule* {
    if (def.analogies.empty() || rules.empty()) return nullptr;
    const std::vector<uint32_t> at(atoms);
    std::vector<std::string> cur;
    for (uint32_t a : at) cur.push_back(N[a]);
    std::vector<const std::vector<std::string>*> alt(at.size(), nullptr);
    bool any = false;
    for (size_t k = 0; k < at.size(); ++k)
      if (auto it = def.analogies.find(T[at[k]]); it != def.analogies.end()) alt[k] = &it->second, any = true;
    if (!any) return nullptr;
    auto attempt = [&]() -> const FFRule* {
      std::vector<const std::string*> p;
      for (const auto& x : cur) p.push_back(&x);
      const FFRule* r = last_match(rules, p, true, rev);
      if (r && used) {
        *used = "";
        for (const auto& x : cur) *used += (used->empty() ? "" : " ") + x;
      }
      return r;
    };
    for (size_t k = 0; k < at.size(); ++k)
      if (alt[k])
        for (const auto& x : *alt[k]) {
          const std::string keep = cur[k];
          cur[k] = x;
          if (const FFRule* r = attempt()) return r;
          cur[k] = keep;
        }
    for (size_t k = 0; k < at.size(); ++k)
      for (size_t m = k + 1; m < at.size(); ++m)
        if (alt[k] && alt[m])
          for (const auto& x : *alt[k])
            for (const auto& y : *alt[m]) {
              const std::string kk = cur[k], mm = cur[m];
              cur[k] = x, cur[m] = y;
              if (const FFRule* r = attempt()) return r;
              cur[k] = kk, cur[m] = mm;
            }
    return nullptr;
  };
  auto estimated = [&](const std::string& what) {
    ++rep.estimated_terms;
    if (std::find(rep.estimated.begin(), rep.estimated.end(), what) == rep.estimated.end()) rep.estimated.push_back(what);
  };
  // LJ: last self rule matching each type (by its vdW equivalent); explicit pairs as overrides
  std::vector<std::string> vdw_name;
  for (const auto& tn : ff.type_names) {
    const FFType* ft = byname.at(tn);
    auto it = ft->equiv.find("vdw");
    vdw_name.push_back(it == ft->equiv.end() ? tn : it->second);
  }
  std::vector<PairType> lj14;
  bool any14 = false;
  for (size_t ti = 0; ti < ff.type_names.size(); ++ti) {
    const std::string& tn = ff.type_names[ti];
    const FFRule* r = nullptr;
    for (auto it = def.pairs.rbegin(); it != def.pairs.rend() && !r; ++it)
      if (it->match.size() == 1 && (glob_match(it->match[0], vdw_name[ti]) || (fallback && glob_match(it->match[0], tn)))) r = &*it;
    if (!r || r->params.size() < 2) {
      // a type that only has explicit pair entries (Buckingham ion pairs) has no self term of its own
      bool explicit_pairs = false;
      for (const auto& q : def.pairs)
        explicit_pairs = explicit_pairs || (q.match.size() == 2 && (glob_match(q.match[0], vdw_name[ti]) || glob_match(q.match[1], vdw_name[ti])));
      if (!explicit_pairs) rep.missing.push_back("pair " + tn);
      ff.lj.push_back({0, 0});
      lj14.push_back({0, 0});
      continue;
    }
    // any Lennard-Jones variant (lj/cut, lj/charmm/coul/long, ... as moltemplate names hybrid sub-styles) mixes
    if (!r->style.empty() && r->style.rfind("lj", 0) != 0)
      throw FFError("pair style '" + r->style + "' on a single type (" + r->name + ") has no mixing rule; give it per pair");
    ff.lj.push_back({r->params[0], r->params[1]});
    if (auto c14 = r->cross.find("lj14"); c14 != r->cross.end() && c14->second.size() >= 2) {
      lj14.push_back({c14->second[0], c14->second[1]});
      any14 = true;
    } else {
      lj14.push_back({r->params[0], r->params[1]});
    }
    rep.used["pair " + tn]++;
  }
  if (any14) ff.lj14_types = lj14;
  for (const auto& r : def.pairs) {
    if (r.match.size() != 2 || r.params.size() < 2) continue;
    for (size_t a = 0; a < ff.type_names.size(); ++a)
      for (size_t b = a; b < ff.type_names.size(); ++b)
        if ((glob_match(r.match[0], vdw_name[a]) && glob_match(r.match[1], vdw_name[b])) ||
            (glob_match(r.match[1], vdw_name[a]) && glob_match(r.match[0], vdw_name[b]))) {
          if (r.style == "buck" || r.style == "born" || r.style == "morse") {
            if (r.params.size() < 3) throw FFError("pair " + r.name + ": " + r.style + " needs 3 parameters");
            ff.pair_func[{int(a), int(b)}] = {r.style == "morse" ? 2 : 1, r.params[0], r.params[1], r.params[2]};
          } else if (r.style.empty() || r.style.rfind("lj", 0) == 0) {
            ff.pair_override[{int(a), int(b)}] = {r.params[0], r.params[1]};
          } else {
            throw FFError("pair style '" + r.style + "' (" + r.name + ") is not supported yet");
          }
        }
  }
  if (def.pair_style.find("class2") != std::string::npos) {
    ff.pair_form = "lj9-6";
    ff.mixing = "sixthpower";   // LAMMPS lj/class2 always mixes by the sixth-power rule
  }

  // charges
  if (charges == "types" && !def.bond_increments.empty()) {
    // bond increments: each bond moves charge between its two atoms
    ff.charge.assign(n, 0.0);
    for (size_t i = 0; i < n; ++i)
      if (!std::isnan(FT[i]->charge)) ff.charge[i] = FT[i]->charge;
    int unmatched = 0;
    for (const auto& b : s.bonds) {
      bool rev = false;
      const FFRule* r = lookup(def.bond_increments, Nq, {b.i, b.j}, &rev);
      // then the bond equivalents, as COMPASS's tools match increments (c43, c44, c4o are c4 for bonds: c4-c43 is c4-c4)
      if (!r) r = lookup(def.bond_increments, Nb, {b.i, b.j}, &rev);
      if (!r || r->params.size() < 2) {
        if (!unmatched++) rep.missing.push_back("bond increment " + shown(Nq, {b.i, b.j}));
        continue;
      }
      ff.charge[b.i] += rev ? r->params[1] : r->params[0];
      ff.charge[b.j] += rev ? r->params[0] : r->params[1];
    }
    rep.notes.push_back("charges from bond increments" + (unmatched ? " (" + std::to_string(unmatched) + " bonds without an increment)" : std::string()));
  } else if (charges == "types") {
    for (size_t i = 0; i < n; ++i) {
      if (std::isnan(FT[i]->charge)) throw FFError(def.name + " has no charge for type " + T[i] + "; use charges \"keep\" or \"gasteiger\"");
      ff.charge.push_back(FT[i]->charge);
    }
  } else if (charges == "keep") {
    for (const auto& a : s.atoms) ff.charge.push_back(a.charge);
  } else if (charges == "gasteiger") {
    std::vector<char> arom(n, 0);
    for (size_t i = 0; i < n; ++i) arom[i] = FT[i]->description.find("aromatic") != std::string::npos;
    ff.charge = gasteiger_ch(s, arom);
    rep.notes.push_back("Gasteiger–Marsili charges (H, C, N, O, halogens, sp³ S)");
  } else if (charges == "qeq") {
    QEqReport qr;
    ff.charge = qeq_charges(s, QEqOptions{}, &qr);
    rep.notes.push_back(qr.notes.front());
  } else {
    throw FFError("charges must be \"types\", \"keep\", \"gasteiger\" or \"qeq\"");
  }

  const auto nb = s.neighbours();
  // reference values for class II cross terms given as separate rules: the assigned bond lengths and angles
  std::map<std::pair<uint32_t, uint32_t>, double> ref_r0;
  std::map<std::tuple<uint32_t, uint32_t, uint32_t>, double> ref_t0;   // (min end, vertex, max end) → θ0 (rad)
  auto r0_of = [&](uint32_t a, uint32_t b) {
    auto it = ref_r0.find({std::min(a, b), std::max(a, b)});
    return it == ref_r0.end() ? 0.0 : it->second;
  };
  auto t0_of = [&](uint32_t a, uint32_t v, uint32_t b) {
    auto it = ref_t0.find({std::min(a, b), v, std::max(a, b)});
    return it == ref_t0.end() ? 0.0 : it->second;
  };
  auto cross_lookup = [&](const char* g, const std::vector<std::string>& N, std::initializer_list<uint32_t> atoms, bool* rev) -> const FFRule* {
    auto it = def.cross_rules.find(g);
    if (it == def.cross_rules.end()) return nullptr;
    return lookup(it->second, N, atoms, rev);
  };
  auto cross_params = [&](const FFRule* r, size_t nparam) {
    std::vector<double> v = r ? r->params : std::vector<double>();
    v.resize(nparam, 0.0);
    return v;
  };
  auto missing = [&](const std::string& what) {
    if (std::find(rep.missing.begin(), rep.missing.end(), what) == rep.missing.end()) rep.missing.push_back(what);
  };
  // bonds
  for (const auto& b : s.bonds) {
    const FFRule* r = lookup12(def.bonds, Nb, Nb2, {b.i, b.j});
    if (!r) r = auto_lookup(def.auto_bonds, {&Ab[b.i], &Ab[b.j]}, nullptr);
    if (std::string u; !r && (r = analog(def.bonds, Nb, {b.i, b.j}, nullptr, &u))) estimated("bond " + shown(Nb, {b.i, b.j}) + " as " + u);
    if (!r) { missing("bond " + shown(Nb, {b.i, b.j})); continue; }
    const std::string st = r->style.empty() ? def.bond_style : r->style;
    if (st == "class2" && r->params.size() >= 4) {
      ff.bonds2.push_back({b.i, b.j, r->params[0], r->params[1], r->params[2], r->params[3]});
      ref_r0[{std::min(b.i, b.j), std::max(b.i, b.j)}] = r->params[0];
    }
    else if (st == "morse" && r->params.size() >= 3) ff.bonds_x.push_back({b.i, b.j, 1, r->params[0], r->params[1], r->params[2]});
    else if (st == "gromos" && r->params.size() >= 2) ff.bonds_x.push_back({b.i, b.j, 2, r->params[0], r->params[1], 0});
    else if (st != "harmonic" || r->params.size() < 2) throw FFError("bond style '" + st + "' (" + r->name + ") is not supported yet");
    else {
      ff.bonds.push_back({b.i, b.j, r->params[0], r->params[1]});
      ref_r0[{std::min(b.i, b.j), std::max(b.i, b.j)}] = r->params[1];
    }
    rep.used["bond " + r->name]++;
  }
  // angles
  int trans_skipped = 0;
  std::set<std::pair<uint32_t, uint32_t>> ex12, ex13;
  for (const auto& b : s.bonds) ex12.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
  for (uint32_t j = 0; j < n; ++j)
    for (size_t x = 0; x < nb[j].size(); ++x)
      for (size_t y = x + 1; y < nb[j].size(); ++y) {
        uint32_t i = nb[j][x], k = nb[j][y];
        if (i > k) std::swap(i, k);
        ex13.insert({i, k});
        bool rev = false;
        const FFRule* r = lookup12(def.angles, Na, Na2, {i, j, k}, &rev);
        // octahedral / square-planar centres: a rule with θ0 = 90° also matches the trans pairs (θ ≈ 180°), which
        // must carry no angle term (as DL_FIELD writes them); decided from the structure's geometry
        if (r && (r->style.empty() ? def.angle_style : r->style) == "harmonic" && r->params.size() >= 2 && std::fabs(r->params[1] - 90.0) < 5.0) {
          const Vec3 u = s.atoms[i].pos - s.atoms[j].pos, v = s.atoms[k].pos - s.atoms[j].pos;
          const double c = dot(u, v) / (norm(u) * norm(v));
          if (c < std::cos(150.0 * kDeg)) { ++trans_skipped; continue; }
        }
        if (!r) r = auto_lookup(def.auto_angles, {&Aae[i], &Aaa[j], &Aae[k]}, &rev);
        if (std::string u; !r && (r = analog(def.angles, Na, {i, j, k}, &rev, &u))) estimated("angle " + shown(Na, {i, j, k}) + " as " + u);
        if (!r) { missing("angle " + shown(Na, {i, j, k})); continue; }
        const std::string st = r->style.empty() ? def.angle_style : r->style;
        if (st == "cosine") {   // K (1 + cos θ), minimum at 180°
          ff.angles_x.push_back({i, j, k, 2, r->params.at(0), kPi});
        } else if (st == "cosine/squared") {   // K (cos θ − cos θ0)²
          if (r->params.size() < 2) throw FFError("cosine/squared angle " + r->name + " needs K theta0");
          ff.angles_x.push_back({i, j, k, 1, r->params[0], r->params[1] * kDeg});
        } else if (st == "charmm") {   // harmonic + Urey–Bradley: K (θ − θ0)² + Kub (r13 − r_ub)²
          const auto& p = r->params;
          if (p.size() < 4) throw FFError("charmm angle " + r->name + " needs K theta0 K_ub r_ub");
          ff.angles.push_back({i, j, k, p[0], p[1] * kDeg});
          if (p[2] != 0) ff.urey_bradley.push_back({i, k, p[2], p[3]});
        } else if (st == "quartic") {   // LAMMPS angle_style quartic: the class II angle without cross terms
          const auto& p = r->params;
          if (p.size() < 4) throw FFError("quartic angle " + r->name + " needs theta0 K2 K3 K4");
          ff.angles2.push_back({i, j, k, p[0] * kDeg, p[1], p[2], p[3], 0, 0, 0, 0, 0, 0, 0});
          ref_t0[{std::min(i, k), j, std::max(i, k)}] = p[0] * kDeg;
        } else if (st == "class2") {
          if (r->params.size() < 4) throw FFError("class2 angle " + r->name + " needs theta0 K2 K3 K4");
          const auto& p = r->params;
          auto bb = group(*r, "bb", 3), ba = group(*r, "ba", 4);
          const uint32_t a0 = rev ? k : i, a2 = rev ? i : k;   // oriented as the rule
          if (!def.cross_rules.empty()) {
            // separate bond-bond and bond-angle rules (.frc): their own match, references from the assigned bonds
            bool rb = false, ra = false;
            const auto m = cross_params(cross_lookup("bb", Na, {a0, j, a2}, &rb), 1);
            const auto na = cross_params(cross_lookup("ba", Na, {a0, j, a2}, &ra), 2);
            const double r1 = r0_of(a0, j), r2 = r0_of(j, a2);
            bb = {m[0], r1, r2};
            ba = {ra ? na[1] : na[0], ra ? na[0] : na[1], r1, r2};
          }
          ff.angles2.push_back({a0, j, a2, p[0] * kDeg, p[1], p[2], p[3], bb[0], bb[1], bb[2], ba[0], ba[1], ba[2], ba[3]});
          ref_t0[{std::min(a0, a2), j, std::max(a0, a2)}] = p[0] * kDeg;
        } else if (st != "harmonic" || r->params.size() < 2) {
          throw FFError("angle style '" + st + "' (" + r->name + ") is not supported yet");
        } else {
          ff.angles.push_back({i, j, k, r->params[0], r->params[1] * kDeg});
          ref_t0[{std::min(i, k), j, std::max(i, k)}] = r->params[1] * kDeg;
        }
        rep.used["angle " + r->name]++;
      }
  // dihedrals (one per i-j-k-l with i < l, as moltemplate's canonical order)
  std::set<std::pair<uint32_t, uint32_t>> p14;
  for (const auto& b : s.bonds)
    for (int dir = 0; dir < 2; ++dir) {
      const uint32_t j = dir ? b.j : b.i, k = dir ? b.i : b.j;
      for (uint32_t i : nb[j]) {
        if (i == k) continue;
        for (uint32_t l : nb[k]) {
          if (l == j || l == i || i > l) continue;   // each dihedral once, in the direction with i < l
          const auto key = std::make_pair(std::min(i, l), std::max(i, l));
          if (!ex12.count(key) && !ex13.count(key)) p14.insert(key);
          bool rev = false;
          const FFRule* r = lookup12(def.dihedrals, Nd, Nd2, {i, j, k, l}, &rev);
          if (!r) r = auto_lookup(def.auto_dihedrals, {&Ate[i], &Atc[j], &Atc[k], &Ate[l]}, &rev);
          if (std::string u; !r && (r = analog(def.dihedrals, Nd, {i, j, k, l}, &rev, &u))) estimated("dihedral " + shown(Nd, {i, j, k, l}) + " as " + u);
          if (!r) { missing("dihedral " + shown(Nd, {i, j, k, l})); continue; }
          const std::string st = r->style.empty() ? def.dihedral_style : r->style;
          const auto& p = r->params;
          // msi2lmp / DL_FIELD (CVFF): a wildcard end spreads the barrier over the torsions about the central bond
          double scale = 1;
          // automatic (cvff_auto) torsions already hold per-torsion values: not scaled
          if (def.wildcard_torsion_scaling == "msi2lmp" && r->match.size() == 4 && r->name.rfind("auto", 0) != 0) {
            const uint32_t end0 = rev ? k : j, end3 = rev ? j : k;   // central atom next to the rule's first / last position
            if (r->match[0] == "*") scale /= std::max<size_t>(1, nb[end0].size() - 1);
            if (r->match[3] == "*") scale /= std::max<size_t>(1, nb[end3].size() - 1);
          }
          auto term = [&](double v, int nn, double d) { if (v != 0) ff.dihedrals.push_back({i, j, k, l, v * scale, nn, d * kDeg}); };
          if (st == "opls") {
            // ½K1(1+cosφ) + ½K2(1−cos2φ) + ½K3(1+cos3φ) + ½K4(1−cos4φ)
            if (p.size() < 4) throw FFError("opls dihedral " + r->name + " needs 4 coefficients");
            term(0.5 * p[0], 1, 0); term(0.5 * p[1], 2, 180); term(0.5 * p[2], 3, 0); term(0.5 * p[3], 4, 180);
          } else if (st == "fourier") {
            const int m = int(p.at(0));
            for (int q = 0; q < m; ++q) term(p.at(1 + 3 * q), int(p.at(2 + 3 * q)), p.at(3 + 3 * q));
          } else if (st == "charmm" || st == "charmmfsw") {
            term(p.at(0), int(p.at(1)), p.at(2));
            if (p.size() > 3 && p[3] != 1.0 && p[3] != 0.0) rep.notes.push_back("charmm dihedral weights other than 0 / 1 are ignored");
          } else if (st == "harmonic") {
            term(p.at(0), int(p.at(2)), p.at(1) >= 0 ? 0 : 180);
          } else if (st == "class2") {
            if (p.size() < 6) throw FFError("class2 dihedral " + r->name + " needs K1 phi1 K2 phi2 K3 phi3");
            const auto mbt = group(*r, "mbt", 4), ebt = group(*r, "ebt", 8), at = group(*r, "at", 8), aat = group(*r, "aat", 3),
                       bb13 = group(*r, "bb13", 3);
            Class2Dihedral d{};
            if (rev) { d.i = l; d.j = k; d.k = j; d.l = i; }
            else { d.i = i; d.j = j; d.k = k; d.l = l; }
            d.k1 = p[0]; d.phi1 = p[1] * kDeg; d.k2 = p[2]; d.phi2 = p[3] * kDeg; d.k3 = p[4]; d.phi3 = p[5] * kDeg;
            for (int q = 0; q < 3; ++q) {
              d.mbt[q] = mbt[q]; d.ebt_b[q] = ebt[q]; d.ebt_c[q] = ebt[3 + q]; d.at_d[q] = at[q]; d.at_e[q] = at[3 + q];
            }
            d.mbt_r2 = mbt[3]; d.ebt_r1 = ebt[6]; d.ebt_r3 = ebt[7];
            d.at_theta1 = at[6] * kDeg; d.at_theta2 = at[7] * kDeg;
            d.aat_m = aat[0]; d.aat_theta1 = aat[1] * kDeg; d.aat_theta2 = aat[2] * kDeg;
            d.bb13_n = bb13[0]; d.bb13_r1 = bb13[1]; d.bb13_r3 = bb13[2];
            if (!def.cross_rules.empty()) {
              // separate cross-term rules (.frc), each matched on its own; end-specific halves swap when it matches
              // the other way round; references from the assigned bonds and angles
              bool r_ = false;
              const auto A = cross_params(cross_lookup("mbt", Nd, {d.i, d.j, d.k, d.l}, &r_), 3);
              bool re = false;
              const auto B = cross_params(cross_lookup("ebt", Nd, {d.i, d.j, d.k, d.l}, &re), 6);
              bool rt = false;
              const auto D = cross_params(cross_lookup("at", Nd, {d.i, d.j, d.k, d.l}, &rt), 6);
              const auto M = cross_params(cross_lookup("aat", Nd, {d.i, d.j, d.k, d.l}, &r_), 1);
              const auto N13 = cross_params(cross_lookup("bb13", Nd, {d.i, d.j, d.k, d.l}, &r_), 1);
              for (int q = 0; q < 3; ++q) {
                d.mbt[q] = A[q];
                d.ebt_b[q] = re ? B[3 + q] : B[q];
                d.ebt_c[q] = re ? B[q] : B[3 + q];
                d.at_d[q] = rt ? D[3 + q] : D[q];
                d.at_e[q] = rt ? D[q] : D[3 + q];
              }
              d.mbt_r2 = r0_of(d.j, d.k);
              d.ebt_r1 = d.bb13_r1 = r0_of(d.i, d.j);
              d.ebt_r3 = d.bb13_r3 = r0_of(d.k, d.l);
              d.at_theta1 = d.aat_theta1 = t0_of(d.i, d.j, d.k);
              d.at_theta2 = d.aat_theta2 = t0_of(d.j, d.k, d.l);
              d.aat_m = M[0];
              d.bb13_n = N13[0];
            }
            ff.dihedrals2.push_back(d);
          } else {
            throw FFError("dihedral style '" + st + "' (" + r->name + ") is not supported yet");
          }
          rep.used["dihedral " + r->name]++;
        }
      }
    }
  for (const auto& q : p14) ff.pairs14.push_back({q.first, q.second});
  // impropers: every centre with three or more neighbours, triples sorted by index; last matching rule wins over all
  // orderings of the three (as moltemplate's canonical ordering does)
  const int cpos = def.improper_order == "center1_sorted" ? 0 : def.improper_order == "center2_sorted" ? 1 : 2;
  if (def.oop_scheme == "msi2lmp") {
    // msi2lmp / Discover class II out-of-plane terms: a three-connected centre B with neighbours (A, C, D) in bond order
    // gets the Wilson term (the rule may match any order of A, C, D) and the three angle-angle couplings; a centre with
    // more neighbours gets angle-angle terms for every triple.
    auto aa_rules = def.cross_rules.find("aa");
    auto aa_one = [&](const std::vector<std::string>& N, uint32_t a, uint32_t b, uint32_t c2, uint32_t d2, bool& found) {
      // K for the pair of angles a-b-c2 and c2-b-d2: the entry (a b c2 d2) or (d2 b c2 a)
      if (aa_rules == def.cross_rules.end()) return 0.0;
      for (auto t4 : {std::array<uint32_t, 4>{a, b, c2, d2}, std::array<uint32_t, 4>{d2, b, c2, a}}) {
        std::vector<const std::string*> ty = {&N[t4[0]], &N[t4[1]], &N[t4[2]], &N[t4[3]]};
        if (const FFRule* r = last_match(aa_rules->second, ty, false)) { found = true; return r->params.empty() ? 0.0 : r->params[0]; }
      }
      return 0.0;
    };
    auto aa_terms = [&](uint32_t A, uint32_t B, uint32_t C, uint32_t D, double m[3]) {
      // M1: ABC & CBD (shared C), M2: ABC & ABD (shared A), M3: ABD & CBD (shared D); own names first, then equivalents
      for (int pass = 0; pass < 2; ++pass) {
        const std::vector<std::string>& N = pass == 0 ? T : Ni;
        bool found = false;
        m[0] = aa_one(N, A, B, C, D, found);
        m[1] = aa_one(N, D, B, A, C, found);
        m[2] = aa_one(N, A, B, D, C, found);
        if (found) return;
      }
    };
    for (uint32_t c = 0; c < n; ++c) {
      const auto& nbc = nb[c];
      if (nbc.size() == 3) {
        const uint32_t A = nbc[0], C = nbc[1], D = nbc[2];
        // msi2lmp's order: permutations a b c d, a b d c, d b a c, d b c a, c b a d, c b d a (outer atoms a c d), each
        // exact before wildcard; own names before the out-of-plane equivalents; never read backwards
        double K = 0, chi0 = 0;
        bool hit = false;
        for (int pass = 0; pass < 2 && !hit; ++pass) {
          const std::vector<std::string>& N = pass == 0 ? T : Ni;
          const uint32_t perm[6][3] = {{A, C, D}, {A, D, C}, {D, A, C}, {D, C, A}, {C, A, D}, {C, D, A}};
          for (const auto& pm : perm) {
            std::vector<const std::string*> ty = {&N[pm[0]], &N[c], &N[pm[1]], &N[pm[2]]};
            if (const FFRule* r = last_match(def.impropers, ty, false)) {
              hit = true;
              K = r->params.size() > 0 ? r->params[0] : 0;
              chi0 = r->params.size() > 1 ? r->params[1] : 0;
              rep.used["improper " + r->name]++;
              break;
            }
          }
        }
        double m[3] = {0, 0, 0};
        aa_terms(A, c, C, D, m);
        ff.impropers2.push_back({A, c, C, D, K, chi0 * kDeg, m[0], m[1], m[2], t0_of(A, c, C), t0_of(A, c, D), t0_of(C, c, D)});
      } else if (nbc.size() > 3) {
        for (size_t x = 0; x + 2 < nbc.size(); ++x)
          for (size_t y = x + 1; y + 1 < nbc.size(); ++y)
            for (size_t z = y + 1; z < nbc.size(); ++z) {
              const uint32_t A = nbc[x], C = nbc[y], D = nbc[z];
              double m[3] = {0, 0, 0};
              aa_terms(A, c, C, D, m);
              ff.impropers2.push_back({A, c, C, D, 0.0, 0.0, m[0], m[1], m[2], t0_of(A, c, C), t0_of(A, c, D), t0_of(C, c, D)});
            }
      }
    }
  }
  for (uint32_t c = 0; def.oop_scheme != "msi2lmp" && c < n; ++c) {
    const auto& nbc = nb[c];
    if (nbc.size() < 3) continue;
    if (def.improper_max_neighbours > 0 && int(nbc.size()) > def.improper_max_neighbours) continue;
    for (size_t x = 0; x < nbc.size(); ++x)
      for (size_t y = x + 1; y < nbc.size(); ++y)
        for (size_t z = y + 1; z < nbc.size(); ++z) {
          std::array<uint32_t, 3> o = {nbc[x], nbc[y], nbc[z]};
          // DL_FIELD / msi2lmp try the outer atoms in the order of the structure's bond list; others sort them
          const bool by_type = def.improper_order == "center3_type_sorted";
          if (by_type)   // AMBER (tleap) and DL_FIELD: outer atoms by type name, ties in bond-list order
            std::stable_sort(o.begin(), o.end(), [&](uint32_t a, uint32_t b) { return T[a] < T[b]; });
          else if (!def.improper_matched_order)
            std::sort(o.begin(), o.end());
          const FFRule* best = nullptr;
          std::array<uint32_t, 3> matched = o;   // the outer atoms in the order the rule matched them
          // find the last rule (in file order) that matches any ordering of the outer atoms
          for (int pass = fallback ? 0 : 1; pass < 3 && !best; ++pass) {
            if (pass == 2 && Ni2.empty()) break;
            const std::vector<std::string>& N = pass == 0 ? T : pass == 1 ? Ni : Ni2;
            for (size_t ri = def.impropers.size(); ri-- > 0 && !best;) {
              const auto& r = def.impropers[ri];
              if (r.match.size() != 4) continue;
              static constexpr int kPerm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
              for (const auto& pi : kPerm) {
                const std::array<uint32_t, 3> pm = {o[pi[0]], o[pi[1]], o[pi[2]]};
                std::vector<const std::string*> ty(4);
                int q = 0;
                for (int pos = 0; pos < 4; ++pos) ty[pos] = pos == cpos ? &N[c] : &N[pm[q++]];
                bool hit = match_fw(r.match, ty);
                if (!hit && def.improper_reversible) {
                  std::vector<const std::string*> rv(ty.rbegin(), ty.rend());
                  hit = match_fw(r.match, rv);
                }
                if (hit) { best = &r; matched = pm; break; }
              }
            }
          }
          if (!best) continue;   // impropers are only where a rule asks for one
          // CHARMM: every explicit rule that matches applies, each in its own atom order (its topologies list e.g.
          // both N1-C5-C2-H1 and N1-C2-C5-H1 for histidine); other force fields use the one best rule
          std::vector<std::pair<const FFRule*, std::array<uint32_t, 3>>> picks = {{best, matched}};
          auto explicit_rule = [](const FFRule& r) { return std::find(r.match.begin(), r.match.end(), std::string("*")) == r.match.end(); };
          if (def.improper_all_explicit && explicit_rule(*best)) {
            static constexpr int kPerm2[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
            for (size_t ri = def.impropers.size(); ri-- > 0;) {
              const auto& r = def.impropers[ri];
              if (&r == best || r.match.size() != 4 || !explicit_rule(r)) continue;
              for (const auto& pi : kPerm2) {
                const std::array<uint32_t, 3> pm = {o[pi[0]], o[pi[1]], o[pi[2]]};
                std::vector<const std::string*> ty(4);
                int q = 0;
                for (int pos = 0; pos < 4; ++pos) ty[pos] = pos == cpos ? &Ni[c] : &Ni[pm[q++]];
                if (match_fw(r.match, ty)) { picks.push_back({&r, pm}); break; }
              }
            }
          }
          for (const auto& pk : picks) {
            const FFRule* best = pk.first;
            const std::array<uint32_t, 3> matched = pk.second;
            {
            // a zero term adds nothing (class II impropers with K = 0 still carry angle-angle terms)
            bool zero = !best->params.empty() && best->params[0] == 0 && best->style != "fourier";
            for (const auto& [g, v] : best->cross)
              for (double x : v) zero = zero && (g != "aa" || x == 0 || &x - v.data() >= 3);
            if (zero) continue;
          }
            uint32_t q4[4];
            int q = 0;
            const std::string st = best->style.empty() ? def.improper_style : best->style;
            // class II angle-angle terms depend on the order of the outer atoms: keep the matched order
            const auto& outer = st == "class2" || (def.improper_matched_order && !by_type) ? matched : o;
            for (int pos = 0; pos < 4; ++pos) q4[pos] = pos == cpos ? c : outer[q++];
            const auto& p = best->params;
            if (st == "fourier") {   // m (K n d)×m, each K [1 + cos(nφ − d)] on the i-j-k-l torsion
              const int m = p.empty() ? 0 : int(p[0]);
              for (int q = 0; q < m; ++q)
                if (p.at(1 + 3 * q) != 0)
                  ff.impropers.push_back({q4[0], q4[1], q4[2], q4[3], p.at(1 + 3 * q), int(p.at(2 + 3 * q)), p.at(3 + 3 * q) * kDeg});
            } else if (st == "planar" || st == "umbrella") {
              // DREIDING planar inversion K (1 − cos ω) averaged over the three bonds (LAMMPS umbrella, ω0 = 0)
              if (p.size() < 1) throw FFError("planar inversion " + best->name + " needs K");
              if (p.size() > 1 && p[1] != 0) rep.notes.push_back("planar inversion " + best->name + ": omega0 ignored (DREIDING planar form)");
              ff.inversions.push_back({c, outer[0], outer[1], outer[2], p[0], 0.0, 1});
            } else if (st == "inversion") {
              if (p.size() < 2) throw FFError("inversion " + best->name + " needs K omega0");
              ff.inversions.push_back({c, outer[0], outer[1], outer[2], p[0], p[1] * kDeg});
            } else if (st == "class2") {
              if (cpos != 1) throw FFError("class2 impropers need the centre in position 2 (improper_order center2_sorted)");
              if (p.size() < 2) throw FFError("class2 improper " + best->name + " needs K chi0");
              const auto aa = group(*best, "aa", 6);
              ff.impropers2.push_back({q4[0], q4[1], q4[2], q4[3], p[0], p[1] * kDeg, aa[0], aa[1], aa[2], aa[3] * kDeg, aa[4] * kDeg,
                                       aa[5] * kDeg});
            } else if (st == "cvff") {
              ff.impropers.push_back({q4[0], q4[1], q4[2], q4[3], p.at(0), int(p.at(2)), p.at(1) >= 0 ? 0.0 : kPi});
            } else if (st == "harmonic") {
              ff.impropers_harmonic.push_back({q4[0], q4[1], q4[2], q4[3], p.at(0), p.at(1) * kDeg});
            } else {
              throw FFError("improper style '" + st + "' (" + best->name + ") is not supported yet");
            }
            rep.used["improper " + best->name]++;
          }
        }
  }
  // exclusions: 1-2, 1-3, 1-4 partners
  ff.excluded.assign(n, {});
  auto add_ex = [&](uint32_t a, uint32_t b) { ff.excluded[a].push_back(b); ff.excluded[b].push_back(a); };
  for (const auto& p : ex12) add_ex(p.first, p.second);
  for (const auto& p : ex13) if (!ex12.count(p)) add_ex(p.first, p.second);
  for (const auto& p : p14) add_ex(p.first, p.second);
  for (auto& e : ff.excluded) {
    std::sort(e.begin(), e.end());
    e.erase(std::unique(e.begin(), e.end()), e.end());
  }
  ff.notes.push_back(def.name + ": " + std::to_string(ff.type_names.size()) + " types, " + std::to_string(ff.bonds.size()) + " bonds, " +
                     std::to_string(ff.angles.size()) + " angles, " + std::to_string(ff.dihedrals.size()) + " torsion terms, " +
                     std::to_string(ff.impropers.size() + ff.impropers_harmonic.size()) + " impropers" +
                     (ff.bonds2.size() + ff.angles2.size() + ff.dihedrals2.size() + ff.impropers2.size()
                          ? "; class II: " + std::to_string(ff.bonds2.size()) + " bonds, " + std::to_string(ff.angles2.size()) + " angles, " +
                                std::to_string(ff.dihedrals2.size()) + " dihedrals, " + std::to_string(ff.impropers2.size()) + " impropers"
                          : std::string()));
  if (trans_skipped) rep.notes.push_back(std::to_string(trans_skipped) + " trans angles (≈180°) at centres with 90° rules carry no angle term");
  if (n_auto) rep.notes.push_back(std::to_string(n_auto) + " interactions use automatic (auto-equivalence) parameters");
  if (rep.estimated_terms)
    rep.notes.push_back(std::to_string(rep.estimated_terms) + " interactions (" + std::to_string(rep.estimated.size()) + " kinds) have no entry in " + def.name +
                        " and use the parameters of analogous types (" + def.analogy_source + "): estimated, listed in the report");
  if (def.pair_style.find("charmm") != std::string::npos)
    rep.notes.push_back("the force field uses " + def.pair_style + " (switched LJ); CAPS evaluates plain LJ truncated at the cut-off");
  if (def.pair_style.find("long") != std::string::npos)
    rep.notes.push_back("the force field uses long-range (Ewald / PPPM) electrostatics; CAPS uses the damped shifted force sum");
  if (!rep.missing.empty() && !allow_missing) {
    std::string l;
    for (size_t k = 0; k < rep.missing.size() && k < 30; ++k) l += "\n  " + rep.missing[k];
    if (rep.missing.size() > 30) l += "\n  … " + std::to_string(rep.missing.size() - 30) + " more";
    if (rep_out) *rep_out = rep;
    throw FFError(def.name + ": no parameters for " + std::to_string(rep.missing.size()) + " interaction(s):" + l);
  }
  if (rep_out) *rep_out = std::move(rep);
  return ff;
}

std::string prepare_for_forcefield(System& s, const FFDef& ff, std::string& charges) {
  std::string note;
  if (!ff.shells.empty()) {
    std::set<std::string> shell_names;
    for (const auto& [c, sh] : ff.shells) shell_names.insert(sh);
    bool have = false;
    for (const auto& a : s.atoms) have = have || shell_names.count(a.name);
    if (!have) {
      const TypingResult tr = assign_types(s, ff);
      const size_t n = s.atoms.size();
      int added = 0;
      for (size_t i = 0; i < n && i < tr.types.size(); ++i) {
        auto it = ff.shells.find(tr.types[i]);
        if (it == ff.shells.end()) continue;
        Atom sh = s.atoms[i];
        sh.id = int64_t(s.atoms.size() + 1);
        sh.name = it->second;
        s.atoms.push_back(sh);
        s.bonds.push_back({uint32_t(i), uint32_t(s.atoms.size() - 1), 1});
        ++added;
      }
      if (added) {
        s.bonds_from_file = true;
        note = ff.name + " is a shell model: " + std::to_string(added) + " shells added, one on each core, bonded to it by the core-shell spring";
      }
    }
  }
  if (ff.keep_defined_bonds) {
    const TypingResult tr = assign_types(s, ff);
    std::vector<Bond> kept;
    for (const auto& b : s.bonds) {
      const std::string& ti = tr.types[b.i];
      const std::string& tj = tr.types[b.j];
      auto bname = [&](const std::string& t) {
        const FFType* ft = ff.type(t);
        if (!ft) return t;
        auto it = ft->equiv.find("bond");
        return it == ft->equiv.end() ? t : it->second;
      };
      const std::string a = bname(ti), c = bname(tj);
      const bool defined = !ti.empty() && !tj.empty() && last_match(ff.bonds, {&a, &c}, true) != nullptr;
      if (defined) kept.push_back(b);
    }
    const size_t dropped = s.bonds.size() - kept.size();
    if (dropped) {
      s.bonds = std::move(kept);
      s.bonds_from_file = true;
      note += (note.empty() ? "" : "; ") + ff.name + " has no bonds between ions: " + std::to_string(dropped) +
              " neighbour bonds dropped, " + std::to_string(s.bonds.size()) + " kept (the terms the force field defines)";
    }
  }
  if (!ff.united_atom) return note;
  bool ch = false;
  for (const auto& b : s.bonds)
    ch = ch || (s.atoms[b.i].element == 6 && s.atoms[b.j].element == 1) || (s.atoms[b.j].element == 6 && s.atoms[b.i].element == 1);
  if (!ch) return "";
  std::string how;
  if (charges == "auto" || charges == "gasteiger") {   // on the all-atom structure, where the method is defined
    ChargeReport q;
    try {
      q = compute_charges(s, "gasteiger");
      how = "Gasteiger–Marsili";
    } catch (const std::exception&) {
      if (charges == "gasteiger") throw;
      q = compute_charges(s, "qeq");
      how = "QEq";
    }
    for (size_t i = 0; i < s.atoms.size() && i < q.q.size(); ++i) s.atoms[i].charge = q.q[i];
    charges = "keep";
  }
  ResolutionReport rep;
  const size_t before = s.atoms.size();
  s = united_atom(s, &rep);
  return ff.name + " is united-atom: " + std::to_string(before - s.atoms.size()) + " hydrogens on carbon folded into their carbons (" +
         std::to_string(s.atoms.size()) + " sites)" + (how.empty() ? "" : "; " + how + " charges computed on the all-atom structure and summed into each site");
}

}  // namespace caps
