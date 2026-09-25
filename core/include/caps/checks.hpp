// CAPS file checks (design board VisProblems): every opened structure is checked, and each finding says what was
// found, what CAPS did about it and what the user can change. Rules: data is never silently dropped; reading goes on
// where it can; one action per finding.
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct FileCheck {
  std::string level;    // "pass", "note", "warn", "error"
  std::string title;    // "Atom counts agree"
  std::string detail;   // what was found and what was done
  std::string action;   // one suggested action id: "wrap", "field", "relax", "" (none)
};

// Checks the trajectory as read: atom counts and ids, bonds (count, longest, stretched), close contacts, charges,
// elements, the cell and atoms outside it, frame spacing, and every note the reader left.
std::vector<FileCheck> file_checks(const Trajectory& t);
std::string file_checks_json(const std::vector<FileCheck>& checks);
std::string file_checks_text(const std::vector<FileCheck>& checks, const std::string& title);

}  // namespace caps
