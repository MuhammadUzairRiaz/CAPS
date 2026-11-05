// CAPS geometry at every match of a SMARTS pattern (see pattern_geometry.hpp).
#include "caps/pattern_geometry.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

#include "caps/typing.hpp"

namespace caps {

std::vector<int> smarts_map_numbers(const std::string& t, std::string* stripped) {
  std::vector<int> maps;
  std::string out;
  size_t i = 0;
  while (i < t.size()) {
    const char c = t[i];
    if (c == '[') {   // one bracket atom; its map is ":n" just before the closing bracket (nested brackets in $() kept)
      int depth = 0;
      size_t j = i;
      for (; j < t.size(); ++j) {
        if (t[j] == '[') ++depth;
        else if (t[j] == ']' && --depth == 0) break;
      }
      if (j >= t.size()) throw std::invalid_argument("SMARTS: an unclosed [");
      std::string body = t.substr(i + 1, j - i - 1);
      int map = 0;
      size_t k = body.size();
      while (k > 0 && std::isdigit(static_cast<unsigned char>(body[k - 1]))) --k;
      if (k > 0 && k < body.size() && body[k - 1] == ':') {
        map = std::stoi(body.substr(k));
        body = body.substr(0, k - 1);
      }
      maps.push_back(map);
      out += "[" + body + "]";
      i = j + 1;
      continue;
    }
    // the organic subset and wildcards outside brackets
    if (c == 'C' && i + 1 < t.size() && t[i + 1] == 'l') { maps.push_back(0), out += "Cl", i += 2; continue; }
    if (c == 'B' && i + 1 < t.size() && t[i + 1] == 'r') { maps.push_back(0), out += "Br", i += 2; continue; }
    if (std::string("BCNOPSFIbcnops*Aa").find(c) != std::string::npos) { maps.push_back(0), out += c, ++i; continue; }
    out += c, ++i;   // bonds, branches, ring closures
  }
  if (stripped) *stripped = out;
  return maps;
}

PatternGeometry pattern_geometry(const Trajectory& t, const std::vector<size_t>& frames, const std::string& smarts, double bin) {
  PatternGeometry R;
  if (frames.empty()) throw std::invalid_argument("pattern geometry: no frames");
  std::string plain;
  const auto maps = smarts_map_numbers(smarts, &plain);
  const Smarts sm(plain);
  if (sm.size() != maps.size()) throw std::invalid_argument("pattern geometry: could not tell the pattern's atoms apart (" + std::to_string(maps.size()) + " read, " +
                                                            std::to_string(sm.size()) + " in the pattern); write each marked atom in brackets");
  // the points: map numbers in increasing order, each the pattern atoms carrying it
  std::map<int, std::vector<size_t>> groups;
  for (size_t k = 0; k < maps.size(); ++k)
    if (maps[k] > 0) groups[maps[k]].push_back(k);
  const int np = int(groups.size());
  if (np < 2 || np > 4) throw std::invalid_argument("pattern geometry: mark 2, 3 or 4 points with atom-map numbers ([S:1][S:2]); this pattern marks " + std::to_string(np));
  R.kind = np == 2 ? "distance" : np == 3 ? "angle" : "torsion";
  R.unit = np == 2 ? "Å" : "°";
  std::vector<std::vector<size_t>> point_atoms;
  for (const auto& [m, ks] : groups) point_atoms.push_back(ks);
  // the matches, each once
  const System& s0 = t.topology;
  const Perception p = perceive(s0);
  std::set<std::vector<std::vector<uint32_t>>> seen;
  for (uint32_t a = 0; a < s0.atoms.size(); ++a)
    for (const auto& e : sm.embeddings(s0, p, a)) {
      std::vector<std::vector<uint32_t>> pts;
      for (const auto& ks : point_atoms) {
        std::vector<uint32_t> at;
        for (size_t k : ks) at.push_back(e[k]);
        pts.push_back(at);
      }
      auto key = pts;
      for (auto& g : key) std::sort(g.begin(), g.end());
      auto rev = key;
      std::reverse(rev.begin(), rev.end());
      if (seen.count(key) || seen.count(rev)) continue;
      seen.insert(key);
      R.matches.push_back(pts);
    }
  if (R.matches.empty()) { R.notes.push_back("the pattern matches nothing in the structure (an unwritten bond means single: write ~ for any bond, : for aromatic)"); return R; }
  // the values
  for (size_t f : frames) {
    const auto& P = t.positions[f];
    const Cell& cell = f < t.cells.size() ? t.cells[f] : s0.cell;
    auto mi = [&](Vec3 d) { return cell.valid() ? cell.minimum_image(d) : d; };
    double fsum = 0;
    for (const auto& pts : R.matches) {
      // each point's centroid, all made whole from the first point's first atom
      const Vec3 ref = P[pts[0][0]];
      std::vector<Vec3> x;
      for (const auto& g : pts) {
        Vec3 c{0, 0, 0};
        for (uint32_t i : g) c = c + mi(P[i] - ref);
        x.push_back(c * (1.0 / double(g.size())));
      }
      double v;
      if (np == 2) v = norm(x[1] - x[0]);
      else if (np == 3) {
        const Vec3 u = x[0] - x[1], w = x[2] - x[1];
        v = std::acos(std::clamp(dot(u, w) / (norm(u) * norm(w)), -1.0, 1.0)) * 180 / M_PI;
      } else {
        const Vec3 b1 = x[1] - x[0], b2 = x[2] - x[1], b3 = x[3] - x[2];
        v = std::atan2(norm(b2) * dot(b1, cross(b2, b3)), dot(cross(b1, b2), cross(b2, b3))) * 180 / M_PI;
      }
      R.values.push_back(v);
      fsum += v;
    }
    R.frame_mean.push_back(fsum / double(R.matches.size()));
  }
  double s = 0, s2 = 0;
  R.min = 1e300, R.max = -1e300;
  for (double v : R.values) s += v, s2 += v * v, R.min = std::min(R.min, v), R.max = std::max(R.max, v);
  const double n = double(R.values.size());
  R.mean = s / n;
  R.sd = std::sqrt(std::max(0.0, s2 / n - R.mean * R.mean));
  const double w = bin > 0 ? bin : np == 2 ? 0.01 : 1.0;
  double lo = std::floor(R.min / w) * w, hi = std::floor(R.max / w) * w + w;
  if (np == 4) lo = -180, hi = 180;
  const int nb = std::max(1, int(std::lround((hi - lo) / w)));
  R.hist_y.assign(size_t(nb), 0.0);
  for (int k = 0; k < nb; ++k) R.hist_x.push_back(lo + (k + 0.5) * w);
  for (double v : R.values) R.hist_y[size_t(std::clamp(int((v - lo) / w), 0, nb - 1))] += 1.0 / (n * w);
  R.notes.push_back(std::to_string(R.matches.size()) + " matches over " + std::to_string(frames.size()) + " frames");
  return R;
}

}  // namespace caps
