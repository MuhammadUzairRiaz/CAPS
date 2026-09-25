#include "caps/yaml.hpp"

#include <cctype>
#include <cstdlib>
#include <stdexcept>
#include <vector>

namespace caps {

namespace {

struct Line {
  int indent;
  std::string text;   // without indentation and comment
  int number;         // 1-based
};

std::string strip_comment(const std::string& s) {
  char quote = 0;
  for (size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (quote) { if (c == quote) quote = 0; continue; }
    if (c == '"' || c == '\'') quote = c;
    else if (c == '#' && (i == 0 || std::isspace(static_cast<unsigned char>(s[i - 1])))) return s.substr(0, i);
  }
  return s;
}

std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}

[[noreturn]] void fail(int line, const std::string& why) { throw std::invalid_argument("line " + std::to_string(line) + ": " + why); }

Json scalar(const std::string& raw, int line) {
  const std::string s = trim(raw);
  if (s.empty() || s == "~" || s == "null") return Json();
  if ((s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\'')) {
    if (s.size() < 2) fail(line, "unclosed quote");
    std::string out;
    for (size_t i = 1; i + 1 < s.size(); ++i) {
      if (s.front() == '"' && s[i] == '\\' && i + 2 < s.size()) {
        const char n = s[++i];
        out += n == 'n' ? '\n' : n == 't' ? '\t' : n;
      } else if (s.front() == '\'' && s[i] == '\'' && i + 2 < s.size() && s[i + 1] == '\'') {
        out += '\'';
        ++i;
      } else {
        out += s[i];
      }
    }
    return Json(out);
  }
  if (s == "true" || s == "yes") return Json(true);
  if (s == "false" || s == "no") return Json(false);
  char* end = nullptr;
  const double v = std::strtod(s.c_str(), &end);
  if (end && *end == 0 && !s.empty() && (std::isdigit(static_cast<unsigned char>(s[0])) || s[0] == '-' || s[0] == '+' || s[0] == '.')) return Json(v);
  return Json(s);
}

// Flow collections: { k: v, … } and [a, …], nested.
Json flow(const std::string& s, size_t& i, int line);
void skip(const std::string& s, size_t& i) { while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i; }
std::string flow_token(const std::string& s, size_t& i) {
  skip(s, i);
  std::string out;
  if (i < s.size() && (s[i] == '"' || s[i] == '\'')) {
    const char q = s[i];
    out += s[i++];
    while (i < s.size() && s[i] != q) out += s[i++];
    if (i < s.size()) out += s[i++];
    return out;
  }
  while (i < s.size() && s[i] != ',' && s[i] != '}' && s[i] != ']' && s[i] != ':') out += s[i++];
  // a colon not followed by a space belongs to the scalar (e.g. a SMILES [*:1])
  while (i < s.size() && s[i] == ':' && i + 1 < s.size() && !std::isspace(static_cast<unsigned char>(s[i + 1])) && s[i + 1] != '}' && s[i + 1] != ',') {
    out += s[i++];
    while (i < s.size() && s[i] != ',' && s[i] != '}' && s[i] != ']' && s[i] != ':') out += s[i++];
  }
  return out;
}
Json flow_value(const std::string& s, size_t& i, int line) {
  skip(s, i);
  if (i < s.size() && (s[i] == '{' || s[i] == '[')) return flow(s, i, line);
  return scalar(flow_token(s, i), line);
}
Json flow(const std::string& s, size_t& i, int line) {
  skip(s, i);
  if (s[i] == '[') {
    ++i;
    Json a = Json::array();
    skip(s, i);
    if (i < s.size() && s[i] == ']') { ++i; return a; }
    for (;;) {
      a.push_back(flow_value(s, i, line));
      skip(s, i);
      if (i >= s.size()) fail(line, "unclosed [");
      if (s[i] == ']') { ++i; return a; }
      if (s[i] != ',') fail(line, "expected , or ] in a list");
      ++i;
    }
  }
  ++i;   // '{'
  Json m = Json::object();
  skip(s, i);
  if (i < s.size() && s[i] == '}') { ++i; return m; }
  for (;;) {
    const std::string key = trim(flow_token(s, i));
    skip(s, i);
    if (i >= s.size() || s[i] != ':') fail(line, "expected key: value in { }");
    ++i;
    m[scalar(key, line).is_string() ? scalar(key, line).str() : key] = flow_value(s, i, line);
    skip(s, i);
    if (i >= s.size()) fail(line, "unclosed {");
    if (s[i] == '}') { ++i; return m; }
    if (s[i] != ',') fail(line, "expected , or } in a mapping");
    ++i;
  }
}

Json value_of(const std::string& raw, int line) {
  const std::string s = trim(raw);
  if (!s.empty() && (s[0] == '{' || s[0] == '[')) {
    size_t i = 0;
    Json v = flow(s, i, line);
    skip(s, i);
    if (i != s.size()) fail(line, "text after a flow collection");
    return v;
  }
  return scalar(s, line);
}

// The first ": " (or ":" at the end) outside quotes and brackets separates key and value.
size_t key_colon(const std::string& s) {
  char quote = 0;
  int depth = 0;
  for (size_t i = 0; i < s.size(); ++i) {
    const char c = s[i];
    if (quote) { if (c == quote) quote = 0; continue; }
    if (c == '"' || c == '\'') quote = c;
    else if (c == '[' || c == '{') ++depth;
    else if (c == ']' || c == '}') --depth;
    else if (c == ':' && depth == 0 && (i + 1 == s.size() || s[i + 1] == ' ')) return i;
  }
  return std::string::npos;
}

Json block(const std::vector<Line>& L, size_t& k, int indent);

// A mapping entry "key: value" (value on the line, or a nested block below).
void entry(const std::vector<Line>& L, size_t& k, int indent, const std::string& text, int number, Json& map) {
  const size_t c = key_colon(text);
  if (c == std::string::npos) fail(number, "expected key: value");
  const std::string key = trim(text.substr(0, c));
  const std::string rest = trim(text.substr(c + 1));
  const Json kj = scalar(key, number);
  const std::string name = kj.is_string() ? kj.str() : key;
  if (map.has(name)) fail(number, "the key '" + name + "' appears twice");
  if (!rest.empty()) { map[name] = value_of(rest, number); return; }
  if (k < L.size() && (L[k].indent > indent || (L[k].indent == indent && L[k].text.rfind("- ", 0) == 0))) map[name] = block(L, k, L[k].indent);
  else map[name] = Json();
}

Json block(const std::vector<Line>& L, size_t& k, int indent) {
  if (k >= L.size()) return Json();
  if (L[k].text == "-" || L[k].text.rfind("- ", 0) == 0) {
    Json list = Json::array();
    while (k < L.size() && L[k].indent == indent && (L[k].text == "-" || L[k].text.rfind("- ", 0) == 0)) {
      const Line& ln = L[k++];
      const std::string item = trim(ln.text.substr(1));
      if (item.empty()) { list.push_back(k < L.size() && L[k].indent > indent ? block(L, k, L[k].indent) : Json()); continue; }
      if (key_colon(item) != std::string::npos && item[0] != '{' && item[0] != '[' && item[0] != '"' && item[0] != '\'') {
        // "- key: value" opens a mapping whose other keys sit at the item's indentation
        Json m = Json::object();
        const int inner = indent + 2;
        entry(L, k, inner, item, ln.number, m);
        while (k < L.size() && L[k].indent == inner && L[k].text.rfind("- ", 0) != 0) {
          const Line& e = L[k++];
          entry(L, k, inner, e.text, e.number, m);
        }
        list.push_back(std::move(m));
      } else {
        list.push_back(value_of(item, ln.number));
      }
    }
    return list;
  }
  Json map = Json::object();
  while (k < L.size() && L[k].indent == indent) {
    const Line& ln = L[k++];
    if (ln.text.rfind("- ", 0) == 0) fail(ln.number, "a list item where a key was expected");
    entry(L, k, indent, ln.text, ln.number, map);
  }
  if (k < L.size() && L[k].indent > indent) fail(L[k].number, "unexpected indentation");
  return map;
}

}  // namespace

Json yaml_parse(const std::string& text) {
  std::vector<Line> L;
  int number = 0;
  size_t start = 0;
  while (start <= text.size()) {
    size_t end = text.find('\n', start);
    if (end == std::string::npos) end = text.size();
    std::string raw = text.substr(start, end - start);
    ++number;
    start = end + 1;
    if (!raw.empty() && raw.back() == '\r') raw.pop_back();
    if (raw.find('\t') != std::string::npos && raw.find_first_not_of(" \t") != std::string::npos && raw.find('\t') < raw.find_first_not_of(" \t"))
      fail(number, "tabs are not allowed for indentation");
    const std::string body = strip_comment(raw);
    if (trim(body).empty() || trim(body) == "---") { if (end == text.size()) break; continue; }
    int indent = 0;
    while (indent < int(body.size()) && body[size_t(indent)] == ' ') ++indent;
    L.push_back({indent, trim(body), number});
    if (end == text.size()) break;
  }
  if (L.empty()) return Json::object();
  size_t k = 0;
  Json out = block(L, k, L[0].indent);
  if (k < L.size()) fail(L[k].number, "unexpected indentation");
  return out;
}

}  // namespace caps
