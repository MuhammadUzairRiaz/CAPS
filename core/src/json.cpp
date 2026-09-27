#include "caps/json.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace caps {

double Json::number() const {
  if (kind_ != Number) throw JsonError("expected a number");
  return d_;
}
bool Json::boolean() const {
  if (kind_ != Bool) throw JsonError("expected true or false");
  return b_;
}
const std::string& Json::str() const {
  if (kind_ != String) throw JsonError("expected a string");
  return s_;
}
const std::vector<Json>& Json::items() const {
  if (kind_ != Array) throw JsonError("expected an array");
  return a_;
}
std::vector<Json>& Json::items() {
  if (kind_ != Array) throw JsonError("expected an array");
  return a_;
}
const std::vector<std::pair<std::string, Json>>& Json::members() const {
  if (kind_ != Object) throw JsonError("expected an object");
  return o_;
}
bool Json::has(const std::string& k) const {
  if (kind_ != Object) return false;
  for (const auto& m : o_)
    if (m.first == k) return true;
  return false;
}
const Json& Json::operator[](const std::string& k) const {
  for (const auto& m : members())
    if (m.first == k) return m.second;
  throw JsonError("missing key '" + k + "'");
}
Json& Json::operator[](const std::string& k) {
  if (kind_ == Null) kind_ = Object;
  if (kind_ != Object) throw JsonError("not an object");
  for (auto& m : o_)
    if (m.first == k) return m.second;
  o_.push_back({k, Json()});
  return o_.back().second;
}
size_t Json::size() const { return kind_ == Array ? a_.size() : kind_ == Object ? o_.size() : 0; }
void Json::push_back(Json v) {
  if (kind_ == Null) kind_ = Array;
  if (kind_ != Array) throw JsonError("not an array");
  a_.push_back(std::move(v));
}
// a number, or true / false as 1 / 0 (options written from Python or by hand say passivate: true)
double Json::num(const std::string& k, double def) const {
  if (!has(k)) return def;
  const Json& v = (*this)[k];
  if (v.is_number()) return v.number();
  if (v.kind() == Bool) return v.boolean() ? 1.0 : 0.0;
  return def;
}
std::string Json::text(const std::string& k, const std::string& def) const {
  return has(k) && (*this)[k].is_string() ? (*this)[k].str() : def;
}

namespace {

struct Parser {
  const std::string& t;
  size_t p = 0;
  int line = 1;
  [[noreturn]] void fail(const std::string& why) { throw JsonError("JSON line " + std::to_string(line) + ": " + why); }
  void ws() {
    while (p < t.size()) {
      const char c = t[p];
      if (c == '\n') { ++line; ++p; }
      else if (c == ' ' || c == '\t' || c == '\r') ++p;
      else break;
    }
  }
  Json value() {
    ws();
    if (p >= t.size()) fail("unexpected end");
    const char c = t[p];
    if (c == '{') return object();
    if (c == '[') return array();
    if (c == '"') return Json(string());
    if (t.compare(p, 4, "true") == 0) { p += 4; return Json(true); }
    if (t.compare(p, 5, "false") == 0) { p += 5; return Json(false); }
    if (t.compare(p, 4, "null") == 0) { p += 4; return Json(); }
    if (c == '-' || (c >= '0' && c <= '9')) {
      char* end = nullptr;
      const double d = std::strtod(t.c_str() + p, &end);
      if (end == t.c_str() + p) fail("bad number");
      p = size_t(end - t.c_str());
      return Json(d);
    }
    fail(std::string("unexpected '") + c + "'");
  }
  std::string string() {
    ++p;
    std::string s;
    while (p < t.size() && t[p] != '"') {
      char c = t[p++];
      if (c == '\n') ++line;
      if (c != '\\') { s += c; continue; }
      if (p >= t.size()) fail("bad escape");
      c = t[p++];
      switch (c) {
        case 'n': s += '\n'; break;
        case 't': s += '\t'; break;
        case 'r': s += '\r'; break;
        case 'b': s += '\b'; break;
        case 'f': s += '\f'; break;
        case 'u': {
          if (p + 4 > t.size()) fail("bad \\u escape");
          unsigned cp = unsigned(std::strtoul(t.substr(p, 4).c_str(), nullptr, 16));
          p += 4;
          if (cp >= 0xD800 && cp < 0xDC00 && p + 6 <= t.size() && t[p] == '\\' && t[p + 1] == 'u') {
            const unsigned lo = unsigned(std::strtoul(t.substr(p + 2, 4).c_str(), nullptr, 16));
            p += 6;
            cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
          }
          if (cp < 0x80) s += char(cp);
          else if (cp < 0x800) { s += char(0xC0 | (cp >> 6)); s += char(0x80 | (cp & 0x3F)); }
          else if (cp < 0x10000) { s += char(0xE0 | (cp >> 12)); s += char(0x80 | ((cp >> 6) & 0x3F)); s += char(0x80 | (cp & 0x3F)); }
          else { s += char(0xF0 | (cp >> 18)); s += char(0x80 | ((cp >> 12) & 0x3F)); s += char(0x80 | ((cp >> 6) & 0x3F)); s += char(0x80 | (cp & 0x3F)); }
          break;
        }
        default: s += c;
      }
    }
    if (p >= t.size()) fail("unterminated string");
    ++p;
    return s;
  }
  Json array() {
    ++p;
    Json a = Json::array();
    ws();
    if (p < t.size() && t[p] == ']') { ++p; return a; }
    for (;;) {
      a.push_back(value());
      ws();
      if (p < t.size() && t[p] == ',') { ++p; continue; }
      if (p < t.size() && t[p] == ']') { ++p; return a; }
      fail("expected ',' or ']'");
    }
  }
  Json object() {
    ++p;
    Json o = Json::object();
    ws();
    if (p < t.size() && t[p] == '}') { ++p; return o; }
    for (;;) {
      ws();
      if (p >= t.size() || t[p] != '"') fail("expected a key");
      const std::string k = string();
      ws();
      if (p >= t.size() || t[p] != ':') fail("expected ':'");
      ++p;
      o[k] = value();
      ws();
      if (p < t.size() && t[p] == ',') { ++p; continue; }
      if (p < t.size() && t[p] == '}') { ++p; return o; }
      fail("expected ',' or '}'");
    }
  }
};

void quote(std::string& out, const std::string& s) {
  out += '"';
  for (unsigned char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\t': out += "\\t"; break;
      case '\r': out += "\\r"; break;
      default:
        if (c < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); out += b; }
        else out += char(c);
    }
  }
  out += '"';
}

}  // namespace

Json Json::parse(const std::string& text) {
  Parser ps{text};
  Json v = ps.value();
  ps.ws();
  if (ps.p != text.size()) ps.fail("trailing characters");
  return v;
}

void Json::dump_to(std::string& out, int indent, int depth) const {
  auto nl = [&](int d) {
    if (indent <= 0) return;
    out += '\n';
    out.append(size_t(indent * d), ' ');
  };
  switch (kind_) {
    case Null: out += "null"; break;
    case Bool: out += b_ ? "true" : "false"; break;
    case Number: {
      char b[32];
      if (std::isfinite(d_) && d_ == std::floor(d_) && std::fabs(d_) < 1e15) std::snprintf(b, sizeof b, "%.0f", d_);
      else std::snprintf(b, sizeof b, "%.10g", d_);
      out += b;
      break;
    }
    case String: quote(out, s_); break;
    case Array: {
      // short arrays of scalars stay on one line
      bool flat = a_.size() <= 8;
      for (const auto& v : a_) flat = flat && v.kind_ != Array && v.kind_ != Object;
      out += '[';
      for (size_t i = 0; i < a_.size(); ++i) {
        if (i) out += flat ? ", " : ",";
        if (!flat) nl(depth + 1);
        a_[i].dump_to(out, indent, depth + 1);
      }
      if (!flat && !a_.empty()) nl(depth);
      out += ']';
      break;
    }
    case Object: {
      out += '{';
      for (size_t i = 0; i < o_.size(); ++i) {
        if (i) out += ',';
        nl(depth + 1);
        quote(out, o_[i].first);
        out += indent > 0 ? ": " : ":";
        o_[i].second.dump_to(out, indent, depth + 1);
      }
      if (!o_.empty()) nl(depth);
      out += '}';
      break;
    }
  }
}

std::string Json::dump(int indent) const {
  std::string s;
  dump_to(s, indent, 0);
  return s;
}

}  // namespace caps
