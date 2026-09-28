// CAPS pipelines as YAML (design/boards/SavePipeline): plain text to diff, review and version. A small subset —
//
//   caps_pipeline: 1
//   name: PS melt · structure report
//   source:
//     file: PS_melt.lammpstrj
//     topology: PS_melt.data
//   steps:                                   # in the order they run
//     - unwrap: {}
//     - cluster: {mode: cutoff, cutoff: 3.3, unit: molecules}
//     - select_expression: {expression: "Type == 2 && Position.Z > 13", enabled: false}
//
// Values are numbers, true/false, quoted or bare strings and flow lists [a, b, c].
#include <cctype>
#include <cstdio>
#include <sstream>
#include <stdexcept>

#include "caps/pipeline.hpp"

namespace caps {

namespace {

std::string yaml_scalar(const Json& v) {
  char b[40];
  switch (v.kind()) {
    case Json::Number: std::snprintf(b, sizeof b, "%.10g", v.number()); return b;
    case Json::Bool: return v.boolean() ? "true" : "false";
    case Json::String: {
      const std::string& s = v.str();
      const bool bare = !s.empty() && s.find_first_of(":,{}[]#\"'&*!|>%@`") == std::string::npos && s.front() != ' ' && s.back() != ' ' &&
                        s != "true" && s != "false" && !(std::isdigit(static_cast<unsigned char>(s[0])) || s[0] == '-' || s[0] == '.');
      if (bare) return s;
      std::string q = "\"";
      for (char c : s) { if (c == '"' || c == '\\') q += '\\'; q += c; }
      return q + "\"";
    }
    case Json::Array: {
      std::string out = "[";
      for (size_t k = 0; k < v.size(); ++k) out += (k ? ", " : "") + yaml_scalar(v[k]);
      return out + "]";
    }
    default: return "null";
  }
}

struct Flow {
  const std::string& s;
  size_t p = 0;
  void ws() { while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p]))) ++p; }
  [[noreturn]] void fail(const std::string& what) { throw std::invalid_argument("pipeline YAML: " + what + " near \"" + s.substr(p, 20) + "\""); }
  Json value() {
    ws();
    if (p >= s.size()) fail("value expected");
    if (s[p] == '[') {
      ++p;
      Json a = Json::array();
      ws();
      if (p < s.size() && s[p] == ']') { ++p; return a; }
      for (;;) {
        a.push_back(value());
        ws();
        if (p < s.size() && s[p] == ',') { ++p; continue; }
        if (p < s.size() && s[p] == ']') { ++p; return a; }
        fail("',' or ']' expected");
      }
    }
    if (s[p] == '{') return mapping();
    if (s[p] == '"' || s[p] == '\'') {
      const char q = s[p++];
      std::string out;
      while (p < s.size() && s[p] != q) {
        if (s[p] == '\\' && p + 1 < s.size()) ++p;
        out += s[p++];
      }
      if (p >= s.size()) fail("unterminated string");
      ++p;
      return Json(out);
    }
    const size_t b = p;
    while (p < s.size() && s[p] != ',' && s[p] != '}' && s[p] != ']') ++p;
    std::string t = s.substr(b, p - b);
    while (!t.empty() && std::isspace(static_cast<unsigned char>(t.back()))) t.pop_back();
    if (t == "true") return Json(true);
    if (t == "false") return Json(false);
    if (t == "null" || t == "~") return Json();
    char* end = nullptr;
    const double d = std::strtod(t.c_str(), &end);
    if (!t.empty() && end && *end == '\0') return Json(d);
    return Json(t);
  }
  Json mapping() {
    ws();
    if (p >= s.size() || s[p] != '{') fail("'{' expected");
    ++p;
    Json o = Json::object();
    ws();
    if (p < s.size() && s[p] == '}') { ++p; return o; }
    for (;;) {
      ws();
      const size_t b = p;
      while (p < s.size() && s[p] != ':') ++p;
      if (p >= s.size()) fail("':' expected");
      std::string key = s.substr(b, p - b);
      while (!key.empty() && std::isspace(static_cast<unsigned char>(key.back()))) key.pop_back();
      ++p;
      o[key] = value();
      ws();
      if (p < s.size() && s[p] == ',') { ++p; continue; }
      if (p < s.size() && s[p] == '}') { ++p; return o; }
      fail("',' or '}' expected");
    }
  }
};

std::string trim_line(std::string s) {
  const auto h = s.find(" #");
  if (h != std::string::npos && s.find('"') == std::string::npos) s = s.substr(0, h);
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
  return s;
}

}  // namespace

std::string pipeline_to_yaml(const Pipeline& p, const std::string& name, const std::string& file, const std::string& topology) {
  std::string out = "caps_pipeline: 1\n";
  if (!name.empty()) out += "name: " + yaml_scalar(Json(name)) + "\n";
  if (!file.empty()) {
    out += "source:\n  file: " + yaml_scalar(Json(file)) + "\n";
    if (!topology.empty()) out += "  topology: " + yaml_scalar(Json(topology)) + "\n";
  }
  if (!p.branch.empty()) out += "branch: " + yaml_scalar(Json(p.branch)) + "        # the branch a run shows (steps of others are skipped)\n";
  out += "steps:                 # in the order they run\n";
  for (size_t k = p.steps.size(); k-- > 0;) {
    const auto& s = p.steps[k];
    out += "  - " + s.type + ": {";
    bool first = true;
    for (const auto& [key, v] : s.params.members()) {
      out += (first ? "" : ", ") + key + ": " + yaml_scalar(v);
      first = false;
    }
    if (!s.enabled) out += std::string(first ? "" : ", ") + "enabled: false";
    out += "}\n";
  }
  if (!p.outputs.empty()) {
    out += "outputs:\n";
    for (const auto& o : p.outputs) {
      if (o.kind == "render") {
        char b[80];
        std::snprintf(b, sizeof b, ", size: [%d, %d]}\n", o.width, o.height);
        out += "  - render: {file: " + yaml_scalar(Json(o.path)) + b;
      } else {
        out += "  - " + o.kind + ": " + (o.what.empty() ? "" : o.what + " ") + "-> " + o.path + "\n";
      }
    }
  }
  return out;
}

Pipeline pipeline_from_yaml(const std::string& text, std::string* name, std::string* file, std::string* topology) {
  std::istringstream in(text);
  std::string line;
  Pipeline p;
  std::vector<PipelineStep> in_order;
  bool steps = false, source = false, seen = false, outputs = false;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const std::string t = trim_line(line);
    if (t.find_first_not_of(' ') == std::string::npos || t[t.find_first_not_of(' ')] == '#') continue;
    const bool indented = t[0] == ' ';
    const std::string body = t.substr(t.find_first_not_of(' '));
    if (!indented) {
      steps = source = outputs = false;
      const auto c = body.find(':');
      const std::string key = body.substr(0, c), val = c == std::string::npos ? "" : trim_line(body.substr(c + 1));
      if (key == "caps_pipeline") seen = true;
      else if (key == "name" && name) {
        if (val.empty()) *name = "";
        else { Flow f{val}; const Json v = f.value(); *name = v.kind() == Json::String ? v.str() : val; }
      }
      else if (key == "branch" && !val.empty()) { Flow f{val}; const Json v = f.value(); p.branch = v.kind() == Json::String ? v.str() : val; }
      else if (key == "source") source = true;
      else if (key == "steps") steps = true;
      else if (key == "outputs") outputs = true;
      continue;
    }
    if (outputs) {   // - table: Ree -> ree.csv · - render: {file: view.png, size: [w, h]}
      if (body.rfind("- ", 0) != 0) throw std::invalid_argument("pipeline YAML: an output is \"- kind: what -> file\", not \"" + body + "\"");
      const std::string item = body.substr(2);
      const auto c = item.find(':');
      if (c == std::string::npos) throw std::invalid_argument("pipeline YAML: an output needs a kind: \"" + body + "\"");
      PipelineOutput o;
      o.kind = trim_line(item.substr(0, c));
      std::string rest = trim_line(item.substr(c + 1));
      while (!rest.empty() && rest.front() == ' ') rest.erase(rest.begin());
      if (!rest.empty() && rest.front() == '{') {
        Flow f{rest};
        const Json m = f.mapping();
        o.path = m.text("file", "");
        o.what = m.text("what", "");
        if (m.has("size") && m["size"].is_array() && m["size"].size() == 2) o.width = int(m["size"][0].number()), o.height = int(m["size"][1].number());
      } else {
        const auto arrow = rest.find("->");
        if (arrow == std::string::npos) throw std::invalid_argument("pipeline YAML: an output is \"kind: what -> file\": \"" + body + "\"");
        o.what = trim_line(rest.substr(0, arrow));
        std::string f = rest.substr(arrow + 2);
        while (!f.empty() && f.front() == ' ') f.erase(f.begin());
        o.path = trim_line(f);
      }
      if (o.path.empty()) throw std::invalid_argument("pipeline YAML: the output has no file: \"" + body + "\"");
      p.outputs.push_back(std::move(o));
      continue;
    }
    if (source) {
      const auto c = body.find(':');
      if (c == std::string::npos) continue;
      const std::string key = body.substr(0, c), val = trim_line(body.substr(c + 1));
      Flow f{val};
      const Json v = val.empty() ? Json("") : f.value();
      if (key == "file" && file) *file = v.kind() == Json::String ? v.str() : val;
      if (key == "topology" && topology) *topology = v.kind() == Json::String ? v.str() : val;
      continue;
    }
    if (steps) {
      if (body.rfind("- ", 0) != 0) throw std::invalid_argument("pipeline YAML: a step is \"- type: {…}\", not \"" + body + "\"");
      const std::string item = body.substr(2);
      const auto c = item.find(':');
      PipelineStep s;
      s.type = trim_line(item.substr(0, c == std::string::npos ? item.size() : c));
      std::string rest = c == std::string::npos ? "{}" : trim_line(item.substr(c + 1));
      if (rest.empty()) rest = "{}";
      Flow f{rest};
      const Json params = f.mapping();
      s.params = Json::object();
      for (const auto& [k, v] : params.members()) {
        if (k == "enabled" && v.kind() == Json::Bool) s.enabled = v.boolean();
        else s.params[k] = v;
      }
      in_order.push_back(std::move(s));
    }
  }
  if (!seen) throw std::invalid_argument("not a CAPS pipeline (no caps_pipeline: line)");
  p.steps.assign(in_order.rbegin(), in_order.rend());   // listed top first, run bottom to top
  return p;
}

}  // namespace caps
