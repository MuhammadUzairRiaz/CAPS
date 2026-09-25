// CAPS visualize pipeline (see caps/pipeline.hpp).
#include "caps/pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <numeric>
#include <set>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "cell_list.hpp"

namespace caps {

namespace {

// CAPS chain palette (as the renderer's molecule colours) and the continuous maps.
const unsigned kCat[10] = {0xF0A83C, 0x6CC4D8, 0xDE775D, 0x9B7AD5, 0x7DC884, 0xD6AC5C, 0xE9ECEF, 0x2271DB, 0xC77DBA, 0x8FB8A8};
const unsigned kViridis[9] = {0x440154, 0x472D7B, 0x3B528B, 0x2C728E, 0x21918C, 0x28AE80, 0x5EC962, 0xADDC30, 0xFDE725};
const unsigned kDiverge[5] = {0x2166AC, 0x92C5DE, 0xF7F7F7, 0xF4A582, 0xB2182B};
constexpr unsigned kSelected = 0xE5484D;

unsigned lerp_rgb(unsigned a, unsigned b, double t) {
  auto ch = [&](int s) { return unsigned(std::lround(((a >> s) & 255) + (double(((b >> s) & 255)) - ((a >> s) & 255)) * t)) & 255; };
  return (ch(16) << 16) | (ch(8) << 8) | ch(0);
}
unsigned ramp(const unsigned* stops, int n, double t) {
  if (!std::isfinite(t)) return 0x808080;
  t = std::clamp(t, 0.0, 1.0) * (n - 1);
  const int k = std::min(n - 2, int(t));
  return lerp_rgb(stops[k], stops[k + 1], t - k);
}
unsigned parse_rgb(const std::string& s, unsigned def) {
  std::string t = s;
  if (!t.empty() && t[0] == '#') t = t.substr(1);
  if (t.size() != 6) return def;
  try { return unsigned(std::stoul(t, nullptr, 16)); } catch (...) { return def; }
}
std::string fmt(const char* f, double v) {
  char b[64];
  std::snprintf(b, sizeof b, f, v);
  return b;
}
bool flag(const Json& p, const std::string& k, bool def) { return p.has(k) && p[k].kind() == Json::Bool ? p[k].boolean() : def; }

// Keep only the particles with keep[i]: atoms, bonds, properties, selection, colours and origins shrink together.
void compact(PipelineState& st, const std::vector<char>& keep) {
  const size_t n = st.system.atoms.size();
  std::vector<int64_t> map(n, -1);
  size_t m = 0;
  for (size_t i = 0; i < n; ++i) if (keep[i]) map[i] = int64_t(m++);
  std::vector<Atom> atoms;
  atoms.reserve(m);
  for (size_t i = 0; i < n; ++i) if (keep[i]) atoms.push_back(st.system.atoms[i]);
  std::vector<Bond> bonds;
  for (const auto& b : st.system.bonds)
    if (keep[b.i] && keep[b.j]) bonds.push_back({uint32_t(map[b.i]), uint32_t(map[b.j]), b.order});
  auto shrink = [&](auto& v) {
    auto w = v;
    w.clear();
    for (size_t i = 0; i < n; ++i) if (keep[i]) w.push_back(v[i]);
    v = std::move(w);
  };
  shrink(st.origin);
  shrink(st.selected);
  shrink(st.colour);
  for (auto& [k, v] : st.props) shrink(v);
  if (st.system.velocities.size() == n) shrink(st.system.velocities);
  st.system.atoms = std::move(atoms);
  st.system.bonds = std::move(bonds);
}

std::vector<std::vector<uint32_t>> adjacency(const System& s) {
  std::vector<std::vector<uint32_t>> adj(s.atoms.size());
  for (const auto& b : s.bonds) { adj[b.i].push_back(b.j); adj[b.j].push_back(b.i); }
  return adj;
}

// Pairs closer than cutoff (minimum image), each once.
template <class F>
void for_pairs(const System& s, double cutoff, F&& f) {
  if (s.atoms.empty() || cutoff <= 0) return;
  const Grid g(s, cutoff);
  const double c2 = cutoff * cutoff;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    const Vec3& fi = g.frac[i];
    g.for_neighbour_bins(g.bin(fi, 0), g.bin(fi, 1), g.bin(fi, 2), [&](const std::vector<uint32_t>& bin) {
      for (uint32_t j : bin) {
        if (j <= i) continue;
        const Vec3 d = g.sep(i, j);
        const double r2 = dot(d, d);
        if (r2 < c2) f(i, j, std::sqrt(r2), d);
      }
    });
  }
}

double auto_distance(const PipelineState& st, const Vec3& n) {
  const Cell& c = st.system.cell;
  if (c.valid()) return dot(n, c.origin + (c.a + c.b + c.c) * 0.5);
  Vec3 com{0, 0, 0};
  for (const auto& a : st.system.atoms) com = com + a.pos;
  return st.system.atoms.empty() ? 0 : dot(n, com * (1.0 / st.system.atoms.size()));
}

// ---------------------------------------------------------------- steps

void step_select_expression(PipelineState& st, const Json& p, StepStatus& out) {
  const std::string e = p.text("expression", "");
  const auto v = evaluate_expression(st, e);
  size_t k = 0;
  for (size_t i = 0; i < v.size(); ++i) { st.selected[i] = v[i] != 0 && std::isfinite(v[i]); k += st.selected[i]; }
  out.summary = std::to_string(k) + " selected";
  if (k == 0) out.level = "warning";
  st.set_attribute("ExpressionSelection.count", double(k));
}

void step_expand_selection(PipelineState& st, const Json& p, StepStatus& out) {
  const std::string mode = p.text("mode", "bonds");
  const int iters = std::max(1, int(p.num("iterations", 1)));
  const size_t before = st.selected_count();
  if (mode == "cutoff") {
    const double rc = p.num("cutoff", 3.0);
    for (int it = 0; it < iters; ++it) {
      auto next = st.selected;
      for_pairs(st.system, rc, [&](uint32_t i, uint32_t j, double, const Vec3&) {
        if (st.selected[i]) next[j] = 1;
        if (st.selected[j]) next[i] = 1;
      });
      st.selected = std::move(next);
    }
  } else {
    const auto adj = adjacency(st.system);
    for (int it = 0; it < iters; ++it) {
      auto next = st.selected;
      for (size_t i = 0; i < adj.size(); ++i)
        if (st.selected[i]) for (uint32_t j : adj[i]) next[j] = 1;
      st.selected = std::move(next);
    }
  }
  out.summary = std::to_string(before) + " → " + std::to_string(st.selected_count()) + " selected";
}

void step_delete_selected(PipelineState& st, const Json&, StepStatus& out) {
  std::vector<char> keep(st.selected.size());
  size_t k = 0;
  for (size_t i = 0; i < keep.size(); ++i) { keep[i] = !st.selected[i]; k += st.selected[i]; }
  compact(st, keep);
  std::fill(st.selected.begin(), st.selected.end(), 0);
  out.summary = std::to_string(k) + " deleted";
  if (k == 0) out.level = "warning";
}

void step_slice(PipelineState& st, const Json& p, StepStatus& out) {
  Vec3 n{0, 0, 1};
  if (p.has("normal") && p["normal"].is_array() && p["normal"].size() == 3) n = {p["normal"][0].number(), p["normal"][1].number(), p["normal"][2].number()};
  const double len = norm(n);
  if (len < 1e-12) throw std::invalid_argument("the normal is zero");
  n = n * (1.0 / len);
  const double d = p.has("distance") && p["distance"].is_number() ? p["distance"].number() : auto_distance(st, n);
  const double w = std::max(0.0, p.num("width", 12.0));
  const bool invert = flag(p, "invert", false), select_only = flag(p, "select_only", false);
  std::vector<char> in(st.system.atoms.size());
  size_t k = 0;
  for (size_t i = 0; i < in.size(); ++i) {
    const double s = dot(n, st.system.atoms[i].pos) - d;
    const bool inside = w > 0 ? std::fabs(s) <= w / 2 : s <= 0;
    in[i] = inside != invert;
    k += in[i];
  }
  if (select_only) {
    st.selected = in;
    out.summary = std::to_string(k) + " selected in the slab";
  } else {
    compact(st, in);
    out.summary = std::to_string(k) + " kept · " + fmt("%.1f Å", w) + " slab at " + fmt("%.1f Å", d);
  }
}

void step_colour_coding(PipelineState& st, const Json& p, StepStatus& out) {
  const std::string prop = p.text("property", "Molecule");
  std::vector<double> v;
  if (!property_values(st, prop, v)) throw std::invalid_argument("no property " + prop);
  const bool only_sel = flag(p, "only_selected", false), lighten = flag(p, "lighten_h", true);
  const size_t n = v.size();
  std::string mode = p.text("mode", "auto");
  std::set<double> distinct;
  bool integral = true;
  for (size_t i = 0; i < n; ++i) {
    if (only_sel && !st.selected[i]) continue;
    integral &= std::fabs(v[i] - std::round(v[i])) < 1e-9;
    if (distinct.size() <= 64) distinct.insert(v[i]);
  }
  if (mode == "auto") mode = integral && distinct.size() <= 64 ? "categorical" : "continuous";
  PipelineLegend L;
  L.property = prop;
  if (mode == "categorical") {
    std::map<double, unsigned> col;
    int k = 0;
    for (double x : distinct) col[x] = kCat[k++ % 10];
    for (size_t i = 0; i < n; ++i) {
      if (only_sel && !st.selected[i]) continue;
      unsigned c = col[v[i]];
      if (lighten && st.system.atoms[i].element == 1) c = lerp_rgb(c, 0xFFFFFF, 0.55);
      st.colour[i] = c;
    }
    for (auto& [x, c] : col)
      if (L.entries.size() < 24) L.entries.push_back({prop == "Element" ? std::string(element(int(x)).symbol) : fmt("%g", x), c});
    out.summary = std::to_string(distinct.size()) + " values · categorical";
  } else {
    double lo = 1e300, hi = -1e300;
    for (size_t i = 0; i < n; ++i) {
      if ((only_sel && !st.selected[i]) || !std::isfinite(v[i])) continue;
      lo = std::min(lo, v[i]);
      hi = std::max(hi, v[i]);
    }
    if (p.has("start") && p["start"].is_number()) lo = p["start"].number();
    if (p.has("end") && p["end"].is_number()) hi = p["end"].number();
    if (lo > hi) lo = hi = 0;
    const bool div = p.text("map", "viridis") == "diverging";
    if (div && lo < 0 && hi > 0) { const double m = std::max(-lo, hi); lo = -m; hi = m; }
    const double span = hi - lo > 1e-12 ? hi - lo : 1.0;
    for (size_t i = 0; i < n; ++i) {
      if (only_sel && !st.selected[i]) continue;
      const double t = (v[i] - lo) / span;
      st.colour[i] = div ? ramp(kDiverge, 5, t) : ramp(kViridis, 9, t);
    }
    L.continuous = true;
    L.lo = lo;
    L.hi = hi;
    L.map = div ? "diverging" : "viridis";
    out.summary = fmt("%.4g", lo) + " → " + fmt("%.4g", hi) + (div ? " · diverging" : " · viridis");
  }
  st.legend = L;
  st.has_legend = true;
}

void step_assign_colour(PipelineState& st, const Json& p, StepStatus& out) {
  const unsigned c = parse_rgb(p.text("colour", "#E5484D"), kSelected);
  size_t k = 0;
  for (size_t i = 0; i < st.colour.size(); ++i)
    if (st.selected[i]) { st.colour[i] = c; ++k; }
  if (!flag(p, "keep_selection", false)) std::fill(st.selected.begin(), st.selected.end(), 0);
  out.summary = std::to_string(k) + " coloured";
  if (k == 0) out.level = "warning";
}

void step_cluster(PipelineState& st, const Json& p, StepStatus& out) {
  const System& s = st.system;
  const size_t n = s.atoms.size();
  const bool only_sel = flag(p, "only_selected", false);
  const std::string mode = p.text("mode", "bonds");
  const double rc = p.num("cutoff", 3.2);
  // union–find over the chosen links
  std::vector<uint32_t> parent(n);
  std::iota(parent.begin(), parent.end(), 0u);
  std::function<uint32_t(uint32_t)> find = [&](uint32_t x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; };
  auto unite = [&](uint32_t a, uint32_t b) {
    if (only_sel && (!st.selected[a] || !st.selected[b])) return;
    a = find(a); b = find(b);
    if (a != b) parent[std::max(a, b)] = std::min(a, b);
  };
  const bool heavy = flag(p, "heavy_only", false);
  if (mode == "cutoff")
    for_pairs(s, rc, [&](uint32_t i, uint32_t j, double, const Vec3&) {
      if (!heavy || (s.atoms[i].element != 1 && s.atoms[j].element != 1)) unite(i, j);
    });
  else for (const auto& b : s.bonds) unite(b.i, b.j);
  if (p.text("unit", "atoms") == "molecules") {   // whole molecules: every atom joins its molecule's cluster
    const auto mol = s.molecules();
    std::map<int, uint32_t> first;
    for (uint32_t i = 0; i < n; ++i) {
      if (only_sel && !st.selected[i]) continue;
      auto [it, fresh] = first.emplace(mol[i], i);
      if (!fresh) unite(it->second, i);
    }
  }
  std::map<uint32_t, std::vector<uint32_t>> groups;
  for (uint32_t i = 0; i < n; ++i)
    if (!only_sel || st.selected[i]) groups[find(i)].push_back(i);
  std::vector<std::vector<uint32_t>> cl;
  for (auto& [r, g] : groups) cl.push_back(std::move(g));
  if (flag(p, "sort_by_size", true))
    std::stable_sort(cl.begin(), cl.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
  auto& id = st.props["Cluster"];
  id.assign(n, 0);
  DataTable t;
  t.name = "clusters";
  t.title = "Cluster list";
  t.columns = {"Cluster", "Size", "Mass (g/mol)", "Rg (Å)", "COM.X", "COM.Y", "COM.Z"};
  const auto adj = mode == "cutoff" ? std::vector<std::vector<uint32_t>>{} : adjacency(s);
  for (size_t c = 0; c < cl.size(); ++c) {
    const auto& g = cl[c];
    for (uint32_t i : g) id[i] = double(c + 1);
    // positions made whole from the first atom by minimum-image steps
    std::vector<Vec3> pos;
    pos.reserve(g.size());
    const Vec3 r0 = s.atoms[g[0]].pos;
    double m = 0;
    Vec3 com{0, 0, 0};
    for (uint32_t i : g) {
      const Vec3 r = s.cell.valid() ? r0 + s.cell.minimum_image(s.atoms[i].pos - r0) : s.atoms[i].pos;
      const double mi = s.mass_of(s.atoms[i]);
      pos.push_back(r);
      com = com + r * mi;
      m += mi;
    }
    if (m > 0) com = com * (1.0 / m);
    double rg2 = 0;
    for (size_t k = 0; k < g.size(); ++k) { const Vec3 d = pos[k] - com; rg2 += s.mass_of(s.atoms[g[k]]) * dot(d, d); }
    const double rg = m > 0 ? std::sqrt(rg2 / m) : 0;
    t.rows.push_back({double(c + 1), double(g.size()), m, rg, com[0], com[1], com[2]});
  }
  st.tables.push_back(std::move(t));
  st.set_attribute("ClusterAnalysis.cluster_count", double(cl.size()));
  st.set_attribute("ClusterAnalysis.largest_size", cl.empty() ? 0.0 : double(cl.front().size()));
  out.summary = std::to_string(cl.size()) + " clusters" + (cl.empty() ? "" : " · largest " + std::to_string(cl.front().size()));
  if (flag(p, "colour", false))
    for (size_t i = 0; i < n; ++i)
      if (id[i] > 0) st.colour[i] = kCat[(int(id[i]) - 1) % 10];
}

void step_coordination(PipelineState& st, const Json& p, StepStatus& out) {
  const System& s = st.system;
  const size_t n = s.atoms.size();
  const double rc = std::max(0.5, p.num("cutoff", 3.2));
  const int bins = std::clamp(int(p.num("bins", 200)), 10, 5000);
  const int ea = int(p.num("element_a", 0)), eb = int(p.num("element_b", 0));
  const bool only_sel = flag(p, "only_selected", false);
  auto is_a = [&](uint32_t i) { return (!only_sel || st.selected[i]) && (ea == 0 || s.atoms[i].element == ea); };
  auto is_b = [&](uint32_t i) { return (!only_sel || st.selected[i]) && (eb == 0 || s.atoms[i].element == eb); };
  auto& cn = st.props["Coordination"];
  cn.assign(n, 0);
  std::vector<double> hist(size_t(bins), 0);
  const double dr = rc / bins;
  const bool inter = flag(p, "inter_only", false);
  const auto mol = inter ? s.molecules() : std::vector<int>{};
  for_pairs(s, rc, [&](uint32_t i, uint32_t j, double r, const Vec3&) {
    if (inter && mol[i] == mol[j]) return;
    const bool ab = is_a(i) && is_b(j), ba = is_a(j) && is_b(i);
    if (ab) cn[i] += 1;
    if (ba) cn[j] += 1;
    const int k = std::min(bins - 1, int(r / dr));
    hist[size_t(k)] += (ab ? 1 : 0) + (ba ? 1 : 0);
  });
  size_t na = 0, nb = 0;
  for (uint32_t i = 0; i < n; ++i) { na += is_a(i); nb += is_b(i); }
  DataTable t;
  t.name = "rdf";
  t.title = "Radial distribution g(r)";
  t.columns = {"r (Å)", "g(r)"};
  const double vol = s.cell.valid() ? s.cell.volume() : 0;
  for (int k = 0; k < bins; ++k) {
    const double r0 = k * dr, r1 = r0 + dr;
    const double shell = 4.0 / 3.0 * M_PI * (r1 * r1 * r1 - r0 * r0 * r0);
    const double ideal = vol > 0 && na && nb ? double(na) * double(nb) / vol * shell : 0;
    t.rows.push_back({r0 + dr / 2, ideal > 0 ? hist[size_t(k)] / ideal : 0});
  }
  st.tables.push_back(std::move(t));
  double mean = 0;
  for (uint32_t i = 0; i < n; ++i) if (is_a(i)) mean += cn[i];
  mean = na ? mean / na : 0;
  st.set_attribute("CoordinationAnalysis.mean", mean);
  out.summary = fmt("mean %.2f", mean) + " within " + fmt("%.2f Å", rc) + " · g(r) " + std::to_string(bins) + " bins";
  if (vol <= 0) { out.level = "warning"; out.summary += " · no cell: g(r) not normalised"; }
}

void step_compute_property(PipelineState& st, const Json& p, StepStatus& out) {
  std::string name = p.text("name", "Custom");
  if (name.empty()) throw std::invalid_argument("the property needs a name");
  const auto v = evaluate_expression(st, p.text("expression", "0"));
  const bool only_sel = flag(p, "only_selected", false);
  const size_t n = v.size();
  auto& atoms = st.system.atoms;
  size_t k = 0;
  auto write = [&](auto&& set) { for (size_t i = 0; i < n; ++i) if (!only_sel || st.selected[i]) { set(i, v[i]); ++k; } };
  if (name == "Position.X") write([&](size_t i, double x) { atoms[i].pos[0] = x; });
  else if (name == "Position.Y") write([&](size_t i, double x) { atoms[i].pos[1] = x; });
  else if (name == "Position.Z") write([&](size_t i, double x) { atoms[i].pos[2] = x; });
  else if (name == "Charge") { write([&](size_t i, double x) { atoms[i].charge = x; }); st.system.has_charges = true; }
  else if (name == "Selection") write([&](size_t i, double x) { st.selected[i] = x != 0; });
  else {
    auto& dst = st.props[name];
    if (dst.size() != n) dst.assign(n, 0);
    write([&](size_t i, double x) { dst[i] = x; });
  }
  double lo = 1e300, hi = -1e300;
  for (double x : v) if (std::isfinite(x)) { lo = std::min(lo, x); hi = std::max(hi, x); }
  out.summary = name + " · " + (lo <= hi ? fmt("%.4g", lo) + " … " + fmt("%.4g", hi) : std::string("no finite values"));
}

void step_wrap(PipelineState& st, const Json&, StepStatus& out) {
  if (!st.system.cell.valid()) { out.level = "warning"; out.summary = "no cell: nothing to wrap"; return; }
  size_t k = 0;
  for (auto& a : st.system.atoms) {
    const Vec3 w = st.system.cell.wrap(a.pos);
    if (norm(w - a.pos) > 1e-9) ++k;
    a.pos = w;
  }
  st.system.unwrapped = false;
  out.summary = std::to_string(k) + " moved into the cell";
}

void step_replicate(PipelineState& st, const Json& p, StepStatus& out) {
  const int nx = std::clamp(int(p.num("nx", 2)), 1, 20), ny = std::clamp(int(p.num("ny", 2)), 1, 20), nz = std::clamp(int(p.num("nz", 1)), 1, 20);
  if (!st.system.cell.valid()) throw std::invalid_argument("replicate needs a cell");
  const size_t n = st.system.atoms.size(), images = size_t(nx) * ny * nz;
  if (n * images > 5'000'000) throw std::invalid_argument("more than 5 million particles");
  System& s = st.system;
  int64_t maxid = 0, maxmol = 0;
  for (const auto& a : s.atoms) { maxid = std::max(maxid, a.id); maxmol = std::max(maxmol, a.mol); }
  const auto atoms0 = s.atoms;
  const auto bonds0 = s.bonds;
  const auto origin0 = st.origin;
  const auto sel0 = st.selected;
  const auto col0 = st.colour;
  const auto props0 = st.props;
  s.atoms.clear(); s.bonds.clear(); st.origin.clear(); st.selected.clear(); st.colour.clear();
  for (auto& [k, v] : st.props) v.clear();
  size_t img = 0;
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
      for (int k = 0; k < nz; ++k, ++img) {
        const Vec3 shift = s.cell.a * i + s.cell.b * j + s.cell.c * k;
        const uint32_t base = uint32_t(img * n);
        for (size_t q = 0; q < n; ++q) {
          Atom a = atoms0[q];
          a.pos = a.pos + shift;
          a.id += int64_t(img) * maxid;
          if (a.mol) a.mol += int64_t(img) * maxmol;
          s.atoms.push_back(a);
          st.origin.push_back(origin0[q]);
          st.selected.push_back(sel0[q]);
          st.colour.push_back(col0[q]);
          for (auto& [name, v] : st.props) v.push_back(props0.at(name)[q]);
        }
        for (const auto& b : bonds0) s.bonds.push_back({base + b.i, base + b.j, b.order});
      }
  if (flag(p, "adjust_cell", true)) { s.cell.a = s.cell.a * nx; s.cell.b = s.cell.b * ny; s.cell.c = s.cell.c * nz; }
  s.velocities.clear();
  out.summary = std::to_string(nx) + " × " + std::to_string(ny) + " × " + std::to_string(nz) + " · " + std::to_string(s.atoms.size()) + " particles";
}

void step_histogram(PipelineState& st, const Json& p, StepStatus& out) {
  const std::string prop = p.text("property", "Position.Z");
  std::vector<double> v;
  if (!property_values(st, prop, v)) throw std::invalid_argument("no property " + prop);
  const bool only_sel = flag(p, "only_selected", false);
  const int bins = std::clamp(int(p.num("bins", 40)), 1, 10000);
  double lo = 1e300, hi = -1e300;
  size_t used = 0;
  for (size_t i = 0; i < v.size(); ++i) {
    if ((only_sel && !st.selected[i]) || !std::isfinite(v[i])) continue;
    lo = std::min(lo, v[i]); hi = std::max(hi, v[i]); ++used;
  }
  if (p.has("start") && p["start"].is_number()) lo = p["start"].number();
  if (p.has("end") && p["end"].is_number()) hi = p["end"].number();
  if (!used) { out.level = "warning"; out.summary = "no values"; return; }
  if (hi - lo < 1e-12) hi = lo + 1;
  std::vector<double> h(size_t(bins), 0);
  const double w = (hi - lo) / bins;
  for (size_t i = 0; i < v.size(); ++i) {
    if ((only_sel && !st.selected[i]) || !std::isfinite(v[i]) || v[i] < lo || v[i] > hi) continue;
    h[size_t(std::min(bins - 1, int((v[i] - lo) / w)))] += 1;
  }
  DataTable t;
  t.name = "histogram";
  t.title = "Histogram · " + prop;
  t.columns = {prop, "Count"};
  for (int k = 0; k < bins; ++k) t.rows.push_back({lo + (k + 0.5) * w, h[size_t(k)]});
  st.tables.push_back(std::move(t));
  out.summary = std::to_string(used) + " values in " + std::to_string(bins) + " bins";
}

void step_binning(PipelineState& st, const Json& p, StepStatus& out) {
  const System& s = st.system;
  const int axis = std::clamp(int(p.num("axis", 2)), 0, 2);
  const int bins = std::clamp(int(p.num("bins", 50)), 1, 10000);
  const std::string red = p.text("reduction", "mean");
  const std::string prop = p.text("property", red == "density" ? "Mass" : "Charge");
  std::vector<double> v;
  if (!property_values(st, prop, v)) throw std::invalid_argument("no property " + prop);
  double lo, hi, area = 0;
  if (s.cell.valid()) {
    const Vec3 edge = axis == 0 ? s.cell.a : axis == 1 ? s.cell.b : s.cell.c;
    lo = s.cell.origin[axis];
    hi = lo + edge[axis];
    area = s.cell.volume() / std::max(1e-12, edge[axis]);
  } else {
    lo = 1e300; hi = -1e300;
    for (const auto& a : s.atoms) { lo = std::min(lo, a.pos[axis]); hi = std::max(hi, a.pos[axis]); }
    if (s.atoms.empty() || hi - lo < 1e-9) { lo = 0; hi = 1; }
  }
  const double w = (hi - lo) / bins;
  std::vector<double> sum(size_t(bins), 0), cnt(size_t(bins), 0);
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    double x = s.atoms[i].pos[axis];
    if (s.cell.valid() && s.cell.periodic[size_t(axis)]) x = lo + std::fmod(std::fmod(x - lo, hi - lo) + (hi - lo), hi - lo);
    if (x < lo || x > hi) continue;
    const size_t k = size_t(std::min(bins - 1, int((x - lo) / w)));
    sum[k] += v[i];
    cnt[k] += 1;
  }
  DataTable t;
  t.name = "binning";
  const char* ax = axis == 0 ? "x" : axis == 1 ? "y" : "z";
  t.title = std::string("Profile along ") + ax + " · " + (red == "density" ? "density" : prop);
  t.columns = {std::string(ax) + " (Å)", red == "density" ? "Density (g/cm³)" : red == "sum" ? prop + " (sum)" : prop + " (mean)", "Count"};
  for (int k = 0; k < bins; ++k) {
    double y = sum[size_t(k)];
    if (red == "mean") y = cnt[size_t(k)] > 0 ? y / cnt[size_t(k)] : 0;
    else if (red == "density") y = area > 0 ? y / (area * w) * 1.66053906660 : 0;
    t.rows.push_back({lo + (k + 0.5) * w, y, cnt[size_t(k)]});
  }
  st.tables.push_back(std::move(t));
  out.summary = std::to_string(bins) + " bins along " + ax + " · " + (red == "density" ? std::string("density") : red + " of " + prop);
  if (red == "density" && area <= 0) { out.level = "warning"; out.summary += " · no cell: density not normalised"; }
}

void step_create_bonds(PipelineState& st, const Json& p, StepStatus& out) {
  System& s = st.system;
  const size_t before = s.bonds.size();
  const bool only_sel = flag(p, "only_selected", false);
  std::vector<Bond> made;
  if (p.text("mode", "perceive") == "cutoff") {
    const double rc = p.num("cutoff", 1.6);
    for_pairs(s, rc, [&](uint32_t i, uint32_t j, double r, const Vec3&) {
      if (r < 0.4 || (only_sel && (!st.selected[i] || !st.selected[j]))) return;
      made.push_back({i, j, 0});
    });
  } else {
    BondOptions o;
    o.tolerance = p.num("tolerance", 0.45);
    for (const auto& b : perceive_bonds(s, o))
      if (!only_sel || (st.selected[b.i] && st.selected[b.j])) made.push_back(b);
  }
  if (flag(p, "replace", false)) s.bonds.clear();
  std::set<std::pair<uint32_t, uint32_t>> have;
  for (const auto& b : s.bonds) have.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
  for (const auto& b : made)
    if (have.insert({std::min(b.i, b.j), std::max(b.i, b.j)}).second) s.bonds.push_back(b);
  out.summary = std::to_string(s.bonds.size()) + " bonds (" + (s.bonds.size() >= before ? "+" : "") + std::to_string(long(s.bonds.size()) - long(before)) + ")";
}

void step_unwrap(PipelineState& st, const Json&, StepStatus& out) {
  if (!st.system.cell.valid()) { out.level = "warning"; out.summary = "no cell: nothing to unwrap"; return; }
  const auto before = st.system.atoms;
  make_molecules_whole(st.system);
  st.system.unwrapped = true;
  size_t k = 0;
  for (size_t i = 0; i < before.size(); ++i) k += norm(before[i].pos - st.system.atoms[i].pos) > 1e-6;
  out.summary = std::to_string(k) + " moved · molecules whole";
}

void step_molecule_shape(PipelineState& st, const Json&, StepStatus& out) {
  System whole = st.system;
  if (!whole.unwrapped && whole.cell.valid()) make_molecules_whole(whole);
  const auto shapes = molecule_shapes(whole);
  const auto mol = whole.molecules();
  const size_t n = whole.atoms.size();
  auto& rg = st.props["MoleculeRg"];
  auto& k2 = st.props["MoleculeKappa2"];
  auto& as = st.props["MoleculeAsphericity"];
  rg.assign(n, 0); k2.assign(n, 0); as.assign(n, 0);
  DataTable t;
  t.name = "molecules";
  t.title = "Molecule shape";
  t.columns = {"Molecule", "Atoms", "Mass (g/mol)", "Rg (Å)", "κ²", "Asphericity (Å²)", "COM.X", "COM.Y", "COM.Z"};
  double mrg = 0, mk2 = 0;
  for (size_t m = 0; m < shapes.size(); ++m) {
    const auto& sh = shapes[m];
    const double b = sh.lambda[2] - 0.5 * (sh.lambda[0] + sh.lambda[1]);
    t.rows.push_back({double(m + 1), double(sh.atoms), sh.mass, sh.rg, sh.kappa2, b, sh.com[0], sh.com[1], sh.com[2]});
    mrg += sh.rg;
    mk2 += sh.kappa2;
  }
  for (size_t i = 0; i < n; ++i) {
    const auto& sh = shapes[size_t(mol[i])];
    rg[i] = sh.rg;
    k2[i] = sh.kappa2;
    as[i] = sh.lambda[2] - 0.5 * (sh.lambda[0] + sh.lambda[1]);
  }
  st.tables.push_back(std::move(t));
  const double nm = std::max<size_t>(1, shapes.size());
  st.set_attribute("MoleculeShape.mean_rg", mrg / nm);
  st.set_attribute("MoleculeShape.mean_kappa2", mk2 / nm);
  out.summary = std::to_string(shapes.size()) + " molecules · mean Rg " + fmt("%.2f Å", mrg / nm) + " · κ² " + fmt("%.3f", mk2 / nm);
}

void step_topology(PipelineState& st, const Json& p, StepStatus& out) {
  const System& s = st.system;
  const int bins = std::clamp(int(p.num("bins", 60)), 5, 2000);
  const auto adj = adjacency(s);
  auto sep = [&](uint32_t a, uint32_t b) { const Vec3 d = s.atoms[b].pos - s.atoms[a].pos; return s.cell.valid() ? s.cell.minimum_image(d) : d; };
  std::vector<double> len, ang, dih;
  for (const auto& b : s.bonds) len.push_back(norm(sep(b.i, b.j)));
  for (uint32_t j = 0; j < adj.size(); ++j)
    for (size_t x = 0; x < adj[j].size(); ++x)
      for (size_t y = x + 1; y < adj[j].size(); ++y) {
        const Vec3 u = sep(j, adj[j][x]), v = sep(j, adj[j][y]);
        ang.push_back(std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / M_PI);
      }
  for (const auto& b : s.bonds)
    for (uint32_t a : adj[b.i])
      for (uint32_t d : adj[b.j]) {
        if (a == b.j || d == b.i || a == d) continue;
        const Vec3 b1 = sep(a, b.i), b2 = sep(b.i, b.j), b3 = sep(b.j, d);
        const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
        const Vec3 m1 = cross(n1, b2 * (1.0 / std::max(1e-12, norm(b2))));
        dih.push_back(std::atan2(dot(m1, n2), dot(n1, n2)) * 180 / M_PI);
      }
  auto hist = [&](const std::vector<double>& v, const char* name, const char* title, const char* col, double lo, double hi) {
    DataTable t;
    t.name = name;
    t.title = title;
    t.columns = {col, "Count"};
    if (v.empty()) { st.tables.push_back(std::move(t)); return; }
    if (lo >= hi) {
      lo = *std::min_element(v.begin(), v.end());
      hi = *std::max_element(v.begin(), v.end());
      if (hi - lo < 1e-9) { lo -= 0.05; hi += 0.05; }
    }
    std::vector<double> h(size_t(bins), 0);
    const double w = (hi - lo) / bins;
    for (double x : v) h[size_t(std::clamp(int((x - lo) / w), 0, bins - 1))] += 1;
    for (int k = 0; k < bins; ++k) t.rows.push_back({lo + (k + 0.5) * w, h[size_t(k)]});
    st.tables.push_back(std::move(t));
  };
  hist(len, "bonds", "Bond lengths", "Length (Å)", 0, 0);
  hist(ang, "angles", "Bond angles", "Angle (°)", 0, 180);
  hist(dih, "dihedrals", "Dihedral angles", "Dihedral (°)", -180, 180);
  double ml = 0;
  for (double x : len) ml += x;
  st.set_attribute("Topology.mean_bond", len.empty() ? 0 : ml / len.size());
  out.summary = std::to_string(len.size()) + " bonds · " + std::to_string(ang.size()) + " angles · " + std::to_string(dih.size()) + " dihedrals";
}

// Positions of trajectory frame k as the pipeline sees particles: by origin index, molecules whole.
std::vector<Vec3> frame_positions(const PipelineState& st, size_t k) {
  System f = st.traj->frame(k);
  if (!f.unwrapped && f.cell.valid()) make_molecules_whole(f);
  std::vector<Vec3> out(st.origin.size());
  for (size_t i = 0; i < out.size(); ++i) {
    const int o = st.origin[i];
    out[i] = o >= 0 && size_t(o) < f.atoms.size() ? f.atoms[size_t(o)].pos : st.system.atoms[i].pos;
  }
  return out;
}

void step_displacements(PipelineState& st, const Json& p, StepStatus& out) {
  if (!st.traj || st.traj->frames() < 2) { out.level = "warning"; out.summary = "one frame: no displacements"; return; }
  const std::string ref = p.text("reference", "first");
  size_t rf = 0;
  if (ref == "previous") rf = st.frame > 0 ? size_t(st.frame - 1) : 0;
  else if (ref == "frame") rf = size_t(std::clamp(int(p.num("frame", 0)), 0, int(st.traj->frames()) - 1));
  const auto r0 = frame_positions(st, rf);
  const auto r1 = frame_positions(st, size_t(st.frame));
  const size_t n = st.system.atoms.size();
  auto& dm = st.props["Displacement"];
  auto& dx = st.props["Displacement.X"];
  auto& dy = st.props["Displacement.Y"];
  auto& dz = st.props["Displacement.Z"];
  dm.assign(n, 0); dx.assign(n, 0); dy.assign(n, 0); dz.assign(n, 0);
  // unwrapped displacements: whole molecules drift through the boundary, so follow each atom by minimum image
  // between the two frames only when they are consecutive; else compare the whole-molecule positions as they are
  double msd = 0, mx = 0;
  for (size_t i = 0; i < n; ++i) {
    Vec3 d = r1[i] - r0[i];
    if (!st.system.unwrapped && st.system.cell.valid() && ref == "previous") d = st.system.cell.minimum_image(d);
    dx[i] = d[0]; dy[i] = d[1]; dz[i] = d[2];
    dm[i] = norm(d);
    msd += dot(d, d);
    mx = std::max(mx, dm[i]);
  }
  msd = n ? msd / n : 0;
  st.set_attribute("Displacements.msd", msd);
  st.set_attribute("Displacements.max", mx);
  out.summary = "vs frame " + std::to_string(rf) + " · MSD " + fmt("%.3g Å²", msd) + " · max " + fmt("%.2f Å", mx);
  if (!st.system.unwrapped && ref != "previous") { out.level = "warning"; out.summary += " · wrapped input: long runs need unwrapped coordinates"; }
}

void step_smooth(PipelineState& st, const Json& p, StepStatus& out) {
  if (!st.traj || st.traj->frames() < 2) { out.level = "warning"; out.summary = "one frame: nothing to average"; return; }
  const int w = std::clamp(int(p.num("window", 5)), 1, 1001);
  const int f0 = std::max(0, st.frame - w / 2), f1 = std::min(int(st.traj->frames()) - 1, st.frame + w / 2);
  const size_t n = st.system.atoms.size();
  std::vector<Vec3> sum(n, Vec3{0, 0, 0});
  const auto here = frame_positions(st, size_t(st.frame));
  for (int f = f0; f <= f1; ++f) {
    const auto r = frame_positions(st, size_t(f));
    for (size_t i = 0; i < n; ++i) {
      Vec3 d = r[i] - here[i];
      if (st.system.cell.valid() && !st.system.unwrapped) d = st.system.cell.minimum_image(d);
      sum[i] = sum[i] + d;
    }
  }
  const double k = 1.0 / (f1 - f0 + 1);
  for (size_t i = 0; i < n; ++i) st.system.atoms[i].pos = here[i] + sum[i] * k;
  out.summary = "frames " + std::to_string(f0) + "–" + std::to_string(f1) + " averaged";
}

struct StepDef {
  const char* type;
  const char* title;
  const char* about;
  void (*run)(PipelineState&, const Json&, StepStatus&);
};

const StepDef kSteps[] = {
    {"select_expression", "Expression selection", "select where an expression is true", step_select_expression},
    {"invert_selection", "Invert selection", "selected ↔ not selected", [](PipelineState& st, const Json&, StepStatus& o) {
       for (auto& x : st.selected) x = !x;
       o.summary = std::to_string(st.selected_count()) + " selected";
     }},
    {"clear_selection", "Clear selection", "nothing selected", [](PipelineState& st, const Json&, StepStatus& o) {
       std::fill(st.selected.begin(), st.selected.end(), 0);
       o.summary = "selection cleared";
     }},
    {"expand_selection", "Expand selection", "by bonds or distance", step_expand_selection},
    {"delete_selected", "Delete selected", "remove the selected particles", step_delete_selected},
    {"slice", "Slice", "slab by normal and width", step_slice},
    {"colour_coding", "Colour coding", "any property, categorical or continuous", step_colour_coding},
    {"assign_colour", "Assign colour", "to the selection", step_assign_colour},
    {"cluster", "Cluster analysis", "by bonds or cutoff · sizes, Rg", step_cluster},
    {"coordination", "Coordination & RDF", "neighbours within a cutoff, g(r)", step_coordination},
    {"compute_property", "Compute property", "expression per particle", step_compute_property},
    {"wrap", "Wrap into cell", "fold positions into the cell", step_wrap},
    {"replicate", "Replicate", "periodic images", step_replicate},
    {"histogram", "Histogram", "distribution of a property", step_histogram},
    {"binning", "Spatial binning", "1-D profile along an axis", step_binning},
    {"create_bonds", "Create bonds", "perceived from distances, or by cutoff", step_create_bonds},
    {"unwrap", "Unwrap", "molecules made whole across the boundary", step_unwrap},
    {"molecule_shape", "Molecule shape", "Rg, κ², asphericity per molecule", step_molecule_shape},
    {"topology", "Topology distributions", "bond lengths, angles, dihedrals", step_topology},
    {"displacements", "Displacements", "vs a reference frame, MSD", step_displacements},
    {"smooth", "Smooth trajectory", "positions averaged over a window of frames", step_smooth},
};

}  // namespace

// ---------------------------------------------------------------- state

double PipelineState::attribute(const std::string& name, double def) const {
  for (const auto& [k, v] : attributes) if (k == name) return v;
  return def;
}
void PipelineState::set_attribute(const std::string& name, double v) {
  for (auto& [k, x] : attributes) if (k == name) { x = v; return; }
  attributes.emplace_back(name, v);
}
size_t PipelineState::selected_count() const { return size_t(std::count(selected.begin(), selected.end(), char(1))); }

std::vector<std::array<std::string, 3>> pipeline_step_catalogue() {
  std::vector<std::array<std::string, 3>> out;
  for (const auto& s : kSteps) out.push_back({s.type, s.title, s.about});
  return out;
}
std::string step_title(const std::string& type) {
  for (const auto& s : kSteps) if (type == s.type) return s.title;
  return type;
}

Pipeline pipeline_from_json(const Json& j) {
  Pipeline p;
  const Json* arr = &j;
  if (j.is_object() && j.has("steps")) arr = &j["steps"];
  if (!arr->is_array()) throw std::invalid_argument("a pipeline is a list of steps");
  for (const auto& s : arr->items()) {
    if (!s.is_object() || !s.has("type")) throw std::invalid_argument("every step needs a type");
    PipelineStep st;
    st.type = s["type"].str();
    st.enabled = !s.has("enabled") || s["enabled"].kind() != Json::Bool || s["enabled"].boolean();
    st.params = Json::object();
    for (const auto& [k, v] : s.members())
      if (k != "type" && k != "enabled") st.params[k] = v;
    p.steps.push_back(std::move(st));
  }
  return p;
}

Json pipeline_to_json(const Pipeline& p) {
  Json arr = Json::array();
  for (const auto& s : p.steps) {
    Json o = Json::object();
    o["type"] = s.type;
    o["enabled"] = s.enabled;
    for (const auto& [k, v] : s.params.members()) o[k] = v;
    arr.push_back(std::move(o));
  }
  Json j = Json::object();
  j["steps"] = std::move(arr);
  return j;
}

PipelineState run_pipeline(const System& frame, const Pipeline& p, int frame_index, int64_t timestep, const Trajectory* traj) {
  PipelineState st;
  st.traj = traj;
  st.system = frame;
  const size_t n = frame.atoms.size();
  st.origin.resize(n);
  std::iota(st.origin.begin(), st.origin.end(), 0);
  st.selected.assign(n, 0);
  st.colour.assign(n, kNoColour);
  st.frame = frame_index;
  st.timestep = timestep;
  st.steps.resize(p.steps.size());
  for (size_t k = p.steps.size(); k-- > 0;) {   // bottom to top
    const auto& step = p.steps[k];
    StepStatus& out = st.steps[k];
    out.type = step.type;
    out.title = step_title(step.type);
    if (!step.enabled) { out.level = "off"; out.summary = "off"; continue; }
    const StepDef* def = nullptr;
    for (const auto& d : kSteps) if (step.type == d.type) def = &d;
    if (!def) { out.level = "error"; out.summary = "unknown step " + step.type; continue; }
    try {
      def->run(st, step.params, out);
    } catch (const std::exception& e) {
      out.level = "error";
      out.summary = e.what();
    }
  }
  // global attributes of the result, ahead of the ones the steps added
  const System& s = st.system;
  int nmol = 0;
  s.molecules(&nmol);
  double q = 0;
  for (const auto& a : s.atoms) q += a.charge;
  std::vector<std::pair<std::string, double>> base = {
      {"SourceFrame", double(frame_index)}, {"Timestep", double(timestep)}, {"Particles", double(s.atoms.size())},
      {"Bonds", double(s.bonds.size())}, {"Molecules", double(nmol)}, {"CellVolume", s.cell.valid() ? s.cell.volume() : 0.0},
      {"Mass", s.total_mass()}, {"Density", s.density()}, {"TotalCharge", q}, {"Selected", double(st.selected_count())}};
  base.insert(base.end(), st.attributes.begin(), st.attributes.end());
  st.attributes = std::move(base);
  return st;
}

DataTable pipeline_series(const Trajectory& traj, const Pipeline& p, int stride, bool wrap, const std::function<bool(int, int)>& progress) {
  DataTable t;
  t.name = "series";
  t.title = "Time series";
  stride = std::max(1, stride);
  const int n = int(traj.frames());
  int done = 0, total = (n + stride - 1) / stride;
  for (int f = 0; f < n; f += stride) {
    System s = traj.frame(size_t(f));
    if (!s.unwrapped && s.cell.valid()) make_molecules_whole(s);
    if (wrap && s.cell.valid()) for (auto& a : s.atoms) a.pos = s.cell.wrap(a.pos);
    const auto st = run_pipeline(s, p, f, f < int(traj.timesteps.size()) ? traj.timesteps[size_t(f)] : 0, &traj);
    if (t.columns.empty()) {
      t.columns = {"Frame"};
      for (const auto& [k, v] : st.attributes) if (k != "SourceFrame") t.columns.push_back(k);
    }
    std::vector<double> row = {double(f)};
    for (size_t c = 1; c < t.columns.size(); ++c) row.push_back(st.attribute(t.columns[c], std::nan("")));
    t.rows.push_back(std::move(row));
    if (progress && !progress(++done, total)) break;
  }
  return t;
}

// ---------------------------------------------------------------- properties

std::vector<std::string> property_names(const PipelineState& st) {
  std::vector<std::string> names = {"Identifier", "Index", "Molecule", "Type", "Element", "Mass", "Charge", "Position.X", "Position.Y", "Position.Z",
                                    "Selection", "DistanceToCOM"};
  for (const auto& [k, v] : st.props) names.push_back(k);
  return names;
}

bool property_values(const PipelineState& st, const std::string& name, std::vector<double>& out) {
  const System& s = st.system;
  const size_t n = s.atoms.size();
  out.resize(n);
  auto each = [&](auto&& f) { for (size_t i = 0; i < n; ++i) out[i] = f(i); return true; };
  if (auto it = st.props.find(name); it != st.props.end()) { out = it->second; out.resize(n, 0); return true; }
  if (name == "Identifier" || name == "Particle Identifier" || name == "ParticleIdentifier") return each([&](size_t i) { return double(s.atoms[i].id); });
  if (name == "Index") return each([&](size_t i) { return double(i); });
  if (name == "Molecule" || name == "Molecule Identifier" || name == "MoleculeIdentifier") {
    if (s.has_mol) return each([&](size_t i) { return double(s.atoms[i].mol); });
    const auto mol = s.molecules();
    return each([&](size_t i) { return double(mol[i] + 1); });
  }
  if (name == "Type" || name == "Particle Type" || name == "ParticleType") return each([&](size_t i) { return double(s.atoms[i].type); });
  if (name == "Element") return each([&](size_t i) { return double(s.atoms[i].element); });
  if (name == "Mass") return each([&](size_t i) { return s.mass_of(s.atoms[i]); });
  if (name == "Charge") return each([&](size_t i) { return s.atoms[i].charge; });
  if (name == "Position.X") return each([&](size_t i) { return s.atoms[i].pos[0]; });
  if (name == "Position.Y") return each([&](size_t i) { return s.atoms[i].pos[1]; });
  if (name == "Position.Z") return each([&](size_t i) { return s.atoms[i].pos[2]; });
  if (name == "Selection") return each([&](size_t i) { return st.selected.size() == n ? double(st.selected[i]) : 0.0; });
  if (name == "DistanceToCOM") {
    System whole = s;
    if (!whole.unwrapped) make_molecules_whole(whole);
    const auto shapes = molecule_shapes(whole);
    const auto mol = whole.molecules();
    return each([&](size_t i) { return norm(whole.atoms[i].pos - shapes[size_t(mol[i])].com); });
  }
  out.clear();
  return false;
}

// ---------------------------------------------------------------- results

Json pipeline_result_json(const PipelineState& st) {
  Json j = Json::object();
  Json attrs = Json::array();
  for (const auto& [k, v] : st.attributes) {
    Json a = Json::object();
    a["name"] = k;
    a["value"] = v;
    attrs.push_back(std::move(a));
  }
  j["attributes"] = std::move(attrs);
  Json steps = Json::array();
  for (const auto& s : st.steps) {
    Json o = Json::object();
    o["type"] = s.type;
    o["title"] = s.title;
    o["summary"] = s.summary;
    o["level"] = s.level;
    steps.push_back(std::move(o));
  }
  j["steps"] = std::move(steps);
  Json tables = Json::array();
  for (const auto& t : st.tables) {
    Json o = Json::object();
    o["name"] = t.name;
    o["title"] = t.title;
    Json cols = Json::array();
    for (const auto& c : t.columns) cols.push_back(c);
    o["columns"] = std::move(cols);
    Json rows = Json::array();
    for (const auto& r : t.rows) {
      Json row = Json::array();
      for (double x : r) row.push_back(std::isfinite(x) ? x : 0.0);
      rows.push_back(std::move(row));
    }
    o["rows"] = std::move(rows);
    tables.push_back(std::move(o));
  }
  j["tables"] = std::move(tables);
  if (st.has_legend) {
    Json l = Json::object();
    l["property"] = st.legend.property;
    l["continuous"] = st.legend.continuous;
    l["lo"] = st.legend.lo;
    l["hi"] = st.legend.hi;
    l["map"] = st.legend.map;
    Json e = Json::array();
    for (const auto& [label, rgb] : st.legend.entries) {
      Json o = Json::object();
      o["label"] = label;
      char b[8];
      std::snprintf(b, sizeof b, "#%06X", rgb & 0xFFFFFF);
      o["colour"] = std::string(b);
      e.push_back(std::move(o));
    }
    l["entries"] = std::move(e);
    j["legend"] = std::move(l);
  }
  Json props = Json::array();
  for (const auto& p : property_names(st)) props.push_back(p);
  j["properties"] = std::move(props);
  j["particles"] = double(st.system.atoms.size());
  j["bonds"] = double(st.system.bonds.size());
  j["selected"] = double(st.selected_count());
  return j;
}

Json particles_json(const PipelineState& st, const std::string& filter, size_t offset, size_t count) {
  const System& s = st.system;
  std::vector<size_t> rows;
  const auto keep = evaluate_expression(st, filter);
  for (size_t i = 0; i < s.atoms.size(); ++i) if (keep[i] != 0 && std::isfinite(keep[i])) rows.push_back(i);
  std::vector<std::string> cols = {"Identifier", "Molecule", "Type", "Name", "Element", "Charge", "Position.X", "Position.Y", "Position.Z"};
  for (const auto& [k, v] : st.props) cols.push_back(k);
  Json j = Json::object();
  Json c = Json::array();
  for (const auto& x : cols) c.push_back(x);
  j["columns"] = std::move(c);
  j["total"] = double(rows.size());
  j["in_frame"] = double(s.atoms.size());
  j["offset"] = double(offset);
  std::vector<double> mol;
  property_values(st, "Molecule", mol);
  Json out = Json::array();
  char b[48];
  for (size_t r = offset; r < rows.size() && r < offset + count; ++r) {
    const size_t i = rows[r];
    const Atom& a = s.atoms[i];
    Json row = Json::array();
    row.push_back(std::to_string(a.id));
    std::snprintf(b, sizeof b, "%.0f", mol[i]); row.push_back(std::string(b));
    row.push_back(std::to_string(a.type));
    row.push_back(a.name);
    row.push_back(std::string(element(a.element).symbol));
    std::snprintf(b, sizeof b, "%+.4f", a.charge); row.push_back(std::string(b));
    for (int k = 0; k < 3; ++k) { std::snprintf(b, sizeof b, "%.3f", a.pos[k]); row.push_back(std::string(b)); }
    for (const auto& [k, v] : st.props) { std::snprintf(b, sizeof b, "%.4g", v[i]); row.push_back(std::string(b)); }
    Json o = Json::object();
    o["index"] = double(i);
    o["origin"] = double(st.origin[i]);
    o["selected"] = st.selected[i] != 0;
    o["cells"] = std::move(row);
    out.push_back(std::move(o));
  }
  j["rows"] = std::move(out);
  return j;
}

Json bonds_json(const PipelineState& st, size_t offset, size_t count) {
  const System& s = st.system;
  Json j = Json::object();
  Json c = Json::array();
  for (const char* x : {"Bond", "A", "B", "Names", "Length (Å)", "Order"}) c.push_back(x);
  j["columns"] = std::move(c);
  j["total"] = double(s.bonds.size());
  j["offset"] = double(offset);
  Json out = Json::array();
  char b[48];
  for (size_t k = offset; k < s.bonds.size() && k < offset + count; ++k) {
    const auto& bd = s.bonds[k];
    const Atom &a = s.atoms[bd.i], &bb = s.atoms[bd.j];
    const Vec3 d = s.cell.valid() ? s.cell.minimum_image(bb.pos - a.pos) : bb.pos - a.pos;
    Json row = Json::array();
    row.push_back(std::to_string(k + 1));
    row.push_back(std::to_string(a.id));
    row.push_back(std::to_string(bb.id));
    row.push_back((a.name.empty() ? element(a.element).symbol : a.name) + std::string("–") + (bb.name.empty() ? element(bb.element).symbol : bb.name));
    std::snprintf(b, sizeof b, "%.4f", norm(d)); row.push_back(std::string(b));
    row.push_back(std::to_string(bd.order));
    Json o = Json::object();
    o["cells"] = std::move(row);
    out.push_back(std::move(o));
  }
  j["rows"] = std::move(out);
  return j;
}

}  // namespace caps
