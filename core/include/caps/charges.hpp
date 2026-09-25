// Partial charges (design/boards/Charges): computed without touching the structure, described by group, and applied on
// request. Methods: gasteiger (Gasteiger–Marsili, 6 iterations; H C N O F Cl Br I and sp³ S), qeq (Rappé–Goddard, every
// element), file (a .chg file: one charge per atom per line, or "index charge"; # comments), keep (the charges the
// structure carries). AM1-BCC and RESP need external programs and are not computed here.
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct ChargeGroup {
  std::string name;   // "C aromatic", "H on sp³ C", …
  int n = 0;
  double mean = 0, lo = 0, hi = 0;
};
struct ChargeReport {
  std::string method;
  std::vector<double> q;
  double net = 0, max_abs = 0;
  int formal = 0;                       // the perceived formal charges' sum: what the net charge should be
  std::vector<ChargeGroup> groups;      // largest first
  std::vector<double> edges, counts;    // histogram, 21 bins over ±max|q|
  std::vector<std::string> notes;
};

ChargeReport compute_charges(const System& s, const std::string& method, const std::string& path = "");
ChargeReport describe_charges(const System& s, const std::vector<double>& q, const std::string& method);
std::vector<double> read_charge_file(const std::string& path, size_t atoms);

}  // namespace caps
