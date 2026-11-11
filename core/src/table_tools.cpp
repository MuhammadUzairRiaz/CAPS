// Study-table tools: histograms, multi-column order, an SD file of the rows (table_tools.hpp).
#include "caps/table_tools.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "caps/io.hpp"

namespace caps {

namespace {
// the q-quantile by linear interpolation between order statistics (sorted input)
double quantile(const std::vector<double>& v, double q) {
  const double h = (double(v.size()) - 1) * q;
  const size_t lo = size_t(std::floor(h)), hi = std::min(v.size() - 1, lo + 1);
  return v[lo] + (h - double(lo)) * (v[hi] - v[lo]);
}
}  // namespace

Histogram histogram(const std::vector<double>& values, const std::string& rule, int bins) {
  Histogram H;
  std::vector<double> v;
  for (double x : values) {
    if (std::isfinite(x)) v.push_back(x);
    else ++H.left_out;
  }
  H.n = int(v.size());
  if (v.empty()) throw std::invalid_argument("no finite values to bin");
  std::sort(v.begin(), v.end());
  const double lo = v.front(), hi = v.back(), range = hi - lo;
  const double n = double(v.size());
  int k = 0;
  if (rule == "bins") {
    if (bins < 1) throw std::invalid_argument("bins: give a number of bins of 1 or more");
    k = bins, H.rule = std::to_string(bins) + " bins";
  } else if (rule == "fd" || rule == "scott") {
    double h = 0;
    if (rule == "fd") {
      h = 2 * (quantile(v, 0.75) - quantile(v, 0.25)) / std::cbrt(n);
      H.rule = "Freedman–Diaconis";
    } else {
      const double mean = std::accumulate(v.begin(), v.end(), 0.0) / n;
      double ss = 0;
      for (double x : v) ss += (x - mean) * (x - mean);
      h = v.size() > 1 ? 3.49 * std::sqrt(ss / (n - 1)) / std::cbrt(n) : 0;
      H.rule = "Scott";
    }
    // a tiny spread next to one far value asks for billions of bins (and range/h past INT_MAX): cap the
    // count in double before any int cast and fall back to Sturges above it
    constexpr double kMaxBins = 10000;
    const double want = (h > 0 && range > 0) ? std::ceil(range / h - 1e-12) : 0;
    if (h > 0 && range > 0 && want <= kMaxBins) k = std::max(1, int(want));
    else if (h > 0 && range > 0) H.rule += " gives over 10000 bins: Sturges", k = 0;
    else H.rule += " gives no width: Sturges", k = 0;
  } else if (rule != "sturges") {
    throw std::invalid_argument("bin rule fd, scott, sturges or bins");
  }
  if (k == 0) {
    k = int(std::ceil(std::log2(n))) + 1;
    if (rule == "sturges") H.rule = "Sturges";
  }
  if (range <= 0) k = 1;
  const double w = range > 0 ? range / k : 1.0;
  const double base = range > 0 ? lo : lo - 0.5;
  H.width = w;
  for (int b = 0; b <= k; ++b) H.edges.push_back(base + w * b);
  H.counts.assign(size_t(k), 0);
  for (double x : v) {
    int b = int(std::floor((x - base) / w));
    b = std::clamp(b, 0, k - 1);   // the last bin closed on the right
    ++H.counts[size_t(b)];
  }
  return H;
}

std::vector<size_t> sort_rows(const std::vector<std::vector<double>>& keys, const std::vector<bool>& descending) {
  const size_t n = keys.empty() ? 0 : keys[0].size();
  for (const auto& k : keys)
    if (k.size() != n) throw std::invalid_argument("every key column needs a value for each row");
  std::vector<size_t> order(n);
  std::iota(order.begin(), order.end(), size_t(0));
  std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    for (size_t k = 0; k < keys.size(); ++k) {
      const double x = keys[k][a], y = keys[k][b];
      const bool fx = std::isfinite(x), fy = std::isfinite(y);
      if (fx != fy) return fx;   // not finite last
      if (!fx || x == y) continue;
      const bool desc = k < descending.size() && descending[k];
      return desc ? x > y : x < y;
    }
    return false;
  });
  return order;
}

int write_sdf_table(const std::vector<SdfRow>& rows, const std::string& path, std::vector<std::string>* notes) {
  // every record is built first: a row whose structure cannot be read leaves the output file as it was (no empty file)
  std::ostringstream out;
  int written = 0, skipped = 0;
  for (const auto& r : rows) {
    if (r.structure.empty()) { ++skipped; continue; }
    if (!std::filesystem::exists(r.structure))
      throw std::invalid_argument("cannot open " + r.structure + " (the structure column holds structure file paths, not SMILES)");
    Trajectory t = open_file(r.structure);
    System s = t.frames() ? t.frame(t.frames() - 1) : t.topology;
    if (!r.title.empty()) s.title = r.title;
    for (const auto& [name, value] : r.data)
      if (name.find('>') != std::string::npos || name.find('<') != std::string::npos)
        throw std::invalid_argument("a data item name cannot hold < or >: '" + name + "'");
    write_sdf_record(s, out, r.data);
    ++written;
  }
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot write " + path);
  f << out.str();
  if (!f) throw std::runtime_error("cannot write " + path);
  if (notes && skipped) notes->push_back(std::to_string(skipped) + " row(s) without a structure left out");
  return written;
}

}  // namespace caps
