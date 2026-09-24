#pragma once
#include <string>
#include <vector>

namespace caps {

inline std::string trim(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  const auto e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

inline std::string strip_comment(const std::string& s, std::string* comment = nullptr) {
  const auto p = s.find('#');
  if (p == std::string::npos) return s;
  if (comment) *comment = s.substr(p + 1);
  return s.substr(0, p);
}

inline std::vector<std::string> split(const std::string& s) {
  std::vector<std::string> out;
  size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n')) ++i;
    const size_t j = i;
    while (i < s.size() && !(s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n')) ++i;
    if (i > j) out.push_back(s.substr(j, i - j));
  }
  return out;
}

inline std::string lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(c >= 'A' && c <= 'Z' ? c + 32 : c);
  return s;
}

}  // namespace caps
