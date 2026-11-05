#pragma once
// Study-table tools: histograms with a bin rule, rows ordered by several columns, and an SD file of the rows that have a
// structure, each with its columns as data items.
//
//  histogram  bins of equal width over [min, max]; the width from Freedman & Diaconis (h = 2 IQR n^(−1/3), Z.
//             Wahrscheinlichkeitstheorie verw. Gebiete 57, 453 (1981)), Scott (h = 3.49 s n^(−1/3), Biometrika 66, 605
//             (1979)), Sturges (k = ⌈log₂ n⌉ + 1 bins, J. Am. Stat. Assoc. 21, 65 (1926)) or a number of bins; the last
//             bin closed on the right. Values that are not finite are left out (and counted).
//  sort_rows  a stable order by the first key, ties by the next …, each ascending or descending; a value that is not
//             finite goes last whichever the direction.
//  write_sdf_table  one record per row with a structure file (its last frame), the row's values as data items
//             (>  <name>, Dalby et al., J. Chem. Inf. Comput. Sci. 32, 244 (1992)); rows without one are skipped.
#include <string>
#include <utility>
#include <vector>

namespace caps {

struct Histogram {
  std::vector<double> edges;   // bins + 1
  std::vector<int> counts;
  double width = 0;
  std::string rule;            // what set the bins
  int n = 0, left_out = 0;     // values binned, values not finite
};
// rule: "fd" (Freedman–Diaconis), "scott", "sturges" or "bins" (with bins > 0). A rule that gives no width (all
// values equal, IQR 0) falls back to Sturges, said in `rule`.
Histogram histogram(const std::vector<double>& values, const std::string& rule = "fd", int bins = 0);

// keys[k][r]: the value of key column k in row r; descending[k] for each key. Returns the rows in order.
std::vector<size_t> sort_rows(const std::vector<std::vector<double>>& keys, const std::vector<bool>& descending);

struct SdfRow {
  std::string title;                                    // the record's first line
  std::string structure;                                // a structure file CAPS reads; empty: skipped
  std::vector<std::pair<std::string, std::string>> data;
};
// Returns the records written; notes the rows skipped.
int write_sdf_table(const std::vector<SdfRow>& rows, const std::string& path, std::vector<std::string>* notes = nullptr);

}  // namespace caps
