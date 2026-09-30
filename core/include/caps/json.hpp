// Minimal JSON value, reader and writer for CAPS data files (force fields, typing overrides).
#pragma once
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace caps {

struct JsonError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

class Json {
 public:
  enum Kind { Null, Bool, Number, String, Array, Object };
  Json() = default;
  Json(std::nullptr_t) {}
  Json(bool b) : kind_(Bool), b_(b) {}
  Json(double d) : kind_(Number), d_(d) {}
  Json(int i) : kind_(Number), d_(i) {}
  Json(const char* s) : kind_(String), s_(s) {}
  Json(std::string s) : kind_(String), s_(std::move(s)) {}
  static Json array() { Json j; j.kind_ = Array; return j; }
  static Json object() { Json j; j.kind_ = Object; return j; }

  Kind kind() const { return kind_; }
  bool is_null() const { return kind_ == Null; }
  bool is_number() const { return kind_ == Number; }
  bool is_string() const { return kind_ == String; }
  bool is_array() const { return kind_ == Array; }
  bool is_object() const { return kind_ == Object; }

  double number() const;
  bool boolean() const;
  const std::string& str() const;
  const std::vector<Json>& items() const;
  std::vector<Json>& items();
  // Object members keep insertion order.
  const std::vector<std::pair<std::string, Json>>& members() const;
  bool has(const std::string& k) const;
  const Json& operator[](const std::string& k) const;   // throws when missing
  Json& operator[](const std::string& k);               // inserts on an object
  const Json& operator[](size_t i) const { return items().at(i); }
  size_t size() const;
  void push_back(Json v);

  // Convenience with defaults.
  double num(const std::string& k, double def) const;
  std::string text(const std::string& k, const std::string& def = "") const;

  static Json parse(const std::string& text);
  std::string dump(int indent = 2) const;
  // Numbers with 17 significant digits: read back exactly (a force field file); dump() keeps 10 for reports
  std::string dump_exact(int indent = 0) const;

 private:
  Kind kind_ = Null;
  bool b_ = false;
  double d_ = 0;
  std::string s_;
  std::vector<Json> a_;
  std::vector<std::pair<std::string, Json>> o_;
  void dump_to(std::string& out, int indent, int depth) const;
};

}  // namespace caps
