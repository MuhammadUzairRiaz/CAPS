#include "caps/tags.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "caps/json.hpp"

namespace caps {

int tag_index(const System& s, const std::string& name) {
  for (size_t k = 0; k < s.tags.size(); ++k)
    if (s.tags[k].name == name) return int(k);
  return -1;
}

std::vector<size_t> tag_atoms(const System& s, int k) {
  std::vector<size_t> out;
  if (k < 0 || size_t(k) >= kMaxTags) return out;
  const uint32_t bit = 1u << k;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].tags & bit) out.push_back(i);
  return out;
}

int tag_edit(System& s, const std::string& name, const std::string& colour, const std::vector<size_t>& atoms, const std::string& op) {
  if (name.empty()) throw std::runtime_error("a tag needs a name");
  int k = tag_index(s, name);
  if (k < 0) {
    if (op == "remove") return -1;
    if (s.tags.size() >= kMaxTags) throw std::runtime_error("a structure holds at most 32 tags");
    s.tags.push_back({name, colour.empty() ? "#F0A83C" : colour});
    k = int(s.tags.size() - 1);
  } else if (!colour.empty()) {
    s.tags[size_t(k)].colour = colour;
  }
  const uint32_t bit = 1u << k;
  if (op == "set")
    for (auto& a : s.atoms) a.tags &= ~bit;
  for (size_t i : atoms) {
    if (i >= s.atoms.size()) continue;
    if (op == "remove") s.atoms[i].tags &= ~bit;
    else s.atoms[i].tags |= bit;
  }
  return k;
}

void tag_delete(System& s, int k) {
  if (k < 0 || size_t(k) >= s.tags.size()) return;
  const uint32_t low = (1u << k) - 1u;   // the tags below k stay; those above move down one bit
  for (auto& a : s.atoms) a.tags = (a.tags & low) | ((a.tags >> 1) & ~low);
  s.tags.erase(s.tags.begin() + k);
}

std::string tag_group_name(const std::string& name) {
  std::string g;
  for (unsigned char c : name) g += std::isalnum(c) || c == '_' ? char(c) : '_';
  if (g.empty() || std::isdigit(static_cast<unsigned char>(g[0]))) g = "t_" + g;
  return g;
}

std::string id_ranges(const std::vector<size_t>& atoms) {
  std::ostringstream o;
  for (size_t k = 0; k < atoms.size();) {
    size_t e = k;
    while (e + 1 < atoms.size() && atoms[e + 1] == atoms[e] + 1) ++e;
    if (k) o << ' ';
    if (e > k) o << atoms[k] + 1 << ':' << atoms[e] + 1;
    else o << atoms[k] + 1;
    k = e + 1;
  }
  return o.str();
}

namespace {
std::vector<size_t> parse_ranges(const std::string& text, size_t n) {
  std::vector<size_t> out;
  std::istringstream in(text);
  std::string w;
  while (in >> w) {
    const auto c = w.find(':');
    try {
      const long a = std::stol(w.substr(0, c)), b = c == std::string::npos ? a : std::stol(w.substr(c + 1));
      for (long i = a; i <= b; ++i)
        if (i >= 1 && size_t(i) <= n) out.push_back(size_t(i - 1));
    } catch (...) {}
  }
  return out;
}
}  // namespace

std::string tags_json(const System& s) {
  Json j = Json::object();
  j["format"] = "caps-tags";
  j["atoms"] = double(s.atoms.size());
  Json a = Json::array();
  for (size_t k = 0; k < s.tags.size(); ++k) {
    Json t = Json::object();
    t["name"] = s.tags[k].name;
    t["colour"] = s.tags[k].colour;
    t["atoms"] = id_ranges(tag_atoms(s, int(k)));
    a.push_back(t);
  }
  j["tags"] = a;
  return j.dump(1);
}

bool tags_from_json(System& s, const std::string& json) {
  Json j;
  try { j = Json::parse(json); } catch (...) { return false; }
  if (!j.is_object() || j.text("format") != "caps-tags" || size_t(j.num("atoms", -1)) != s.atoms.size() || !j.has("tags")) return false;
  s.tags.clear();
  for (auto& a : s.atoms) a.tags = 0;
  for (size_t k = 0; k < j["tags"].size() && k < kMaxTags; ++k) {
    const Json& t = j["tags"][k];
    if (t.text("name").empty()) continue;
    tag_edit(s, t.text("name"), t.text("colour"), parse_ranges(t.text("atoms"), s.atoms.size()), "set");
  }
  return true;
}

std::string tags_sidecar(const std::string& structure_path) { return structure_path + ".tags.json"; }

void write_tags(const System& s, const std::string& structure_path) {
  const auto p = tags_sidecar(structure_path);
  if (s.tags.empty()) {
    std::error_code ec;
    std::filesystem::remove(p, ec);
    return;
  }
  std::ofstream(p) << tags_json(s);
}

bool read_tags(System& s, const std::string& structure_path) {
  std::ifstream in(tags_sidecar(structure_path));
  if (!in) return false;
  std::stringstream ss;
  ss << in.rdbuf();
  return tags_from_json(s, ss.str());
}

std::string lammps_tag_groups(const System& s) {
  if (s.tags.empty()) return "";
  std::ostringstream o;
  o << "\n# tags from CAPS (each a group: compute, dump or fix on it)\n";
  for (size_t k = 0; k < s.tags.size(); ++k) {
    const auto atoms = tag_atoms(s, int(k));
    if (atoms.empty()) continue;
    const std::string head = "group           " + tag_group_name(s.tags[k].name) + " id";
    std::string line = head;
    std::istringstream r(id_ranges(atoms));
    std::string w;
    while (r >> w) {
      if (line.size() + w.size() > 200) { o << line << " &\n"; line = "               "; }
      line += ' ' + w;
    }
    o << line << "\n";
  }
  return o.str();
}

}  // namespace caps
