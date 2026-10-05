#include "caps/charges.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/grow.hpp"
#include "caps/qeq.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {

// the atom's environment: element with hybridisation or aromaticity; hydrogens by the atom they sit on
std::string group_of(const System& s, const Perception& p, uint32_t i) {
  const int z = s.atoms[i].element;
  const std::string el = element(z).symbol;
  auto hyb = [&](uint32_t a) -> std::string {
    if (p.aromatic[a]) return "aromatic";
    const size_t k = p.nb[a].size();
    int maxo = 1;
    for (int o : p.order[a]) maxo = std::max(maxo, o);
    if (s.atoms[a].element == 6) return maxo == 3 || (k == 2 && maxo == 2) ? "sp" : maxo == 2 ? "sp²" : "sp³";
    if (s.atoms[a].element == 8) return k == 1 && maxo == 2 ? "carbonyl" : k == 1 ? "terminal" : "sp³";
    if (s.atoms[a].element == 7) return maxo == 3 ? "sp" : maxo == 2 ? "sp²" : "sp³";
    return "";
  };
  if (z == 1) {
    if (p.nb[i].empty()) return "H";
    const uint32_t a = p.nb[i][0];
    const std::string h = hyb(a);
    return "H on " + (h.empty() ? "" : h + " ") + element(s.atoms[a].element).symbol;
  }
  const std::string h = hyb(i);
  return h.empty() ? el : el + " " + h;
}

}  // namespace

std::vector<double> read_charge_file(const std::string& path, size_t atoms) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot read " + path);
  std::vector<double> q(atoms, std::nan(""));
  std::string line;
  size_t next = 0;
  int lineno = 0;
  while (std::getline(f, line)) {
    ++lineno;
    if (const auto h = line.find('#'); h != std::string::npos) line.resize(h);
    std::istringstream is(line);
    std::vector<double> v;
    double x;
    while (is >> x) v.push_back(x);
    if (v.empty()) continue;
    size_t i = next;
    double c = v.back();
    if (v.size() >= 2) i = size_t(std::llround(v[0])) - 1;   // "index charge", 1-based
    if (i >= atoms) throw std::runtime_error(path + ":" + std::to_string(lineno) + ": atom " + std::to_string(i + 1) + " is not in the structure");
    q[i] = c;
    next = i + 1;
  }
  const auto missing = std::count_if(q.begin(), q.end(), [](double c) { return std::isnan(c); });
  if (missing) throw std::runtime_error(path + " gives no charge for " + std::to_string(missing) + " of " + std::to_string(atoms) + " atoms");
  return q;
}

ChargeReport describe_charges(const System& s, const std::vector<double>& q, const std::string& method) {
  ChargeReport r;
  r.method = method;
  r.q = q;
  const Perception p = perceive(s);
  for (int c : p.charge) r.formal += c;
  std::map<std::string, std::vector<double>> by;
  for (uint32_t i = 0; i < s.atoms.size() && i < q.size(); ++i) {
    r.net += q[i];
    r.max_abs = std::max(r.max_abs, std::fabs(q[i]));
    by[group_of(s, p, i)].push_back(q[i]);
  }
  for (const auto& [name, v] : by) {
    ChargeGroup g;
    g.name = name;
    g.n = int(v.size());
    double sum = 0;
    g.lo = 1e300, g.hi = -1e300;
    for (double x : v) sum += x, g.lo = std::min(g.lo, x), g.hi = std::max(g.hi, x);
    g.mean = sum / double(v.size());
    r.groups.push_back(g);
  }
  std::sort(r.groups.begin(), r.groups.end(), [](const ChargeGroup& a, const ChargeGroup& b) { return a.n > b.n; });
  const double m = r.max_abs > 0 ? r.max_abs : 1;
  const int bins = 21;
  for (int b = 0; b <= bins; ++b) r.edges.push_back(-m + 2 * m * b / bins);
  r.counts.assign(size_t(bins), 0);
  for (double x : q) r.counts[size_t(std::clamp(int((x + m) / (2 * m) * bins), 0, bins - 1))] += 1;
  return r;
}

ChargeReport compute_charges(const System& s, const std::string& method, const std::string& path) {
  std::vector<double> q;
  std::vector<std::string> notes;
  if (method == "gasteiger") {
    const Perception p = perceive(s);
    std::vector<char> arom(s.atoms.size(), 0);
    for (size_t i = 0; i < arom.size(); ++i) arom[i] = p.aromatic[i] ? 1 : 0;
    q = gasteiger_ch(s, arom);
    notes.push_back("Gasteiger–Marsili, 6 iterations; the net charge follows from the formal charges");
  } else if (method == "qeq") {
    QEqReport rep;
    q = qeq_charges(s, QEqOptions{}, &rep);
    notes.push_back("QEq (Rappé & Goddard 1991), Ohno–Klopman shielding");
  } else if (method == "file") {
    q = read_charge_file(path, s.atoms.size());
    notes.push_back("from " + path);
  } else if (method == "keep") {
    for (const auto& a : s.atoms) q.push_back(a.charge);
  } else {
    throw std::invalid_argument("charges: method gasteiger, qeq, file or keep (AM1-BCC and RESP need AmberTools or a QM program; import their .chg)");
  }
  auto r = describe_charges(s, q, method);
  r.notes = notes;
  return r;
}


std::vector<int> equivalent_atoms(const System& s) {
  const size_t n = s.atoms.size();
  const auto nb = s.neighbours();
  std::vector<int64_t> lab(n);
  for (size_t i = 0; i < n; ++i) lab[i] = s.atoms[i].element;
  size_t classes = 0;
  for (int it = 0; it < 64; ++it) {   // refine: a label is the old label and the sorted labels of the neighbours
    std::map<std::pair<int64_t, std::vector<int64_t>>, int64_t> ids;
    std::vector<int64_t> nl(n);
    for (size_t i = 0; i < n; ++i) {
      std::vector<int64_t> ns;
      for (uint32_t k : nb[i]) ns.push_back(lab[k]);
      std::sort(ns.begin(), ns.end());
      auto key = std::make_pair(lab[i], ns);
      auto f = ids.find(key);
      nl[i] = f != ids.end() ? f->second : (ids[key] = int64_t(ids.size()));
    }
    lab.swap(nl);
    if (ids.size() == classes) break;
    classes = ids.size();
  }
  return std::vector<int>(lab.begin(), lab.end());
}

void adjust_charges(const System& s, std::vector<double>& q, const ChargeAdjust& a, int formal, std::vector<std::string>* notes) {
  const size_t n = std::min(q.size(), s.atoms.size());
  if (a.average) {
    const auto cls = equivalent_atoms(s);
    std::map<int, std::pair<double, int>> sum;
    for (size_t i = 0; i < n; ++i) sum[cls[i]].first += q[i], sum[cls[i]].second += 1;
    size_t shared = 0;
    for (size_t i = 0; i < n; ++i) {
      const auto& p = sum[cls[i]];
      if (p.second > 1) ++shared;
      q[i] = p.first / p.second;
    }
    if (notes) notes->push_back("averaged over " + std::to_string(sum.size()) + " classes of equivalent atoms (" + std::to_string(shared) + " atoms share a class)");
  }
  if (a.scale != 1.0) {
    for (size_t i = 0; i < n; ++i) q[i] *= a.scale;
    if (notes) notes->push_back("scaled × " + std::to_string(a.scale).substr(0, 6));
  }
  if (a.neutralise == "even" || a.neutralise == "proportional") {
    const double target = std::isnan(a.target) ? formal * a.scale : a.target;
    double tot = 0, wsum = 0;
    for (size_t i = 0; i < n; ++i) tot += q[i], wsum += std::fabs(q[i]);
    const double dq = target - tot;
    if (a.neutralise == "proportional" && wsum > 0)
      for (size_t i = 0; i < n; ++i) q[i] += dq * std::fabs(q[i]) / wsum;
    else
      for (size_t i = 0; i < n; ++i) q[i] += dq / double(n);
    if (notes) notes->push_back("total brought to " + std::to_string(target).substr(0, 7) + " e (" + a.neutralise + ", a shift of " +
                                std::to_string(dq).substr(0, 9) + " e in all)");
  }
}

}  // namespace caps
