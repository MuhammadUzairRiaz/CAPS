// CAPS visualize pipeline (see caps/pipeline.hpp).
#include "caps/pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <filesystem>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/bundle.hpp"
#include "caps/elements.hpp"
#include "caps/entangle.hpp"
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
  out.summary = std::to_string(k) + " of " + std::to_string(v.size()) + " selected";
  // what the numeric types mean here, when the expression uses them (Type 2 → ca)
  if (e.find("Type") != std::string::npos) {
    std::map<int, std::string> names;
    for (const auto& a : st.system.atoms) if (a.type > 0 && !a.name.empty() && !names.count(a.type)) names[a.type] = a.name;
    for (const auto& t : st.system.types) if (!t.label.empty()) names[t.type] = t.label;
    std::string legend;
    for (const auto& [t, n] : names) legend += (legend.empty() ? "" : ", ") + std::to_string(t) + " " + n;
    if (!legend.empty() && names.size() <= 12) out.summary += " · types " + legend;
  }
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
    for (double x : distinct) col[x] = prop == "Element" ? element(int(x)).rgb : kCat[k++ % 10];   // elements keep their own colours
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
  // which molecules each cluster holds (molecule ids from 1)
  const auto molid = s.molecules();
  t.label_column = "Molecules";
  for (const auto& g : cl) {
    std::set<int> ms;
    for (uint32_t i : g) ms.insert(molid[i] + 1);
    std::string txt;
    int shown = 0;
    for (int m : ms) {
      if (shown++ == 8) { txt += ", …"; break; }
      txt += (txt.empty() ? "" : ", ") + std::to_string(m);
    }
    t.labels.push_back(txt);
  }
  // cutoff sweep (cutoff mode): how the cluster count and the largest cluster move with the cutoff
  if (mode == "cutoff" && flag(p, "sweep", true)) {
    DataTable sw;
    sw.name = "cluster_sweep";
    sw.title = "Cutoff sweep";
    sw.columns = {"Cutoff (Å)", "Clusters", "Largest (molecules)", "Largest (atoms)"};
    const bool whole = p.text("unit", "atoms") == "molecules";
    for (double f : {rc - 0.2, rc, rc + 0.3, rc + 0.7}) {
      if (f <= 0) continue;
      std::vector<uint32_t> par(n);
      std::iota(par.begin(), par.end(), 0u);
      std::function<uint32_t(uint32_t)> fd = [&](uint32_t x) { while (par[x] != x) x = par[x] = par[par[x]]; return x; };
      auto un = [&](uint32_t a, uint32_t b) {
        if (only_sel && (!st.selected[a] || !st.selected[b])) return;
        a = fd(a), b = fd(b);
        if (a != b) par[std::max(a, b)] = std::min(a, b);
      };
      for_pairs(s, f, [&](uint32_t i, uint32_t j, double, const Vec3&) { if (!heavy || (s.atoms[i].element != 1 && s.atoms[j].element != 1)) un(i, j); });
      if (whole) {
        std::map<int, uint32_t> first;
        for (uint32_t i = 0; i < n; ++i) { auto [it, fresh] = first.emplace(molid[i], i); if (!fresh) un(it->second, i); }
      }
      std::map<uint32_t, std::set<int>> mols;
      std::map<uint32_t, int> atoms;
      for (uint32_t i = 0; i < n; ++i) if (!only_sel || st.selected[i]) mols[fd(i)].insert(molid[i]), ++atoms[fd(i)];
      size_t big = 0;
      int big_atoms = 0;
      for (const auto& [r, m] : mols) if (m.size() > big) big = m.size(), big_atoms = atoms[r];
      sw.rows.push_back({f, double(mols.size()), double(big), double(big_atoms)});
    }
    st.tables.push_back(std::move(sw));
  }
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
  const double rmax = std::max(rc, p.num("rmax", rc));   // g(r) out to rmax; coordination counted within the cutoff
  const int bins = std::clamp(int(p.num("bins", 200)), 10, 5000);
  const int ea = int(p.num("element_a", 0)), eb = int(p.num("element_b", 0));
  const bool only_sel = flag(p, "only_selected", false);
  auto is_a = [&](uint32_t i) { return (!only_sel || st.selected[i]) && (ea == 0 || s.atoms[i].element == ea); };
  auto is_b = [&](uint32_t i) { return (!only_sel || st.selected[i]) && (eb == 0 || s.atoms[i].element == eb); };
  auto& cn = st.props["Coordination"];
  cn.assign(n, 0);
  std::vector<double> hist(size_t(bins), 0);
  const double dr = rmax / bins;
  const bool inter = flag(p, "inter_only", false);
  const auto mol = inter ? s.molecules() : std::vector<int>{};
  for_pairs(s, rmax, [&](uint32_t i, uint32_t j, double r, const Vec3&) {
    if (inter && mol[i] == mol[j]) return;
    const bool ab = is_a(i) && is_b(j), ba = is_a(j) && is_b(i);
    if (r <= rc) {
      if (ab) cn[i] += 1;
      if (ba) cn[j] += 1;
    }
    const int k = std::min(bins - 1, int(r / dr));
    hist[size_t(k)] += (ab ? 1 : 0) + (ba ? 1 : 0);
  });
  // averaged over the trajectory (every `every`-th frame; this frame's selection filter does not apply to the others)
  int frames = 1;
  if (flag(p, "average_frames", false) && st.traj && st.traj->frames() > 1) {
    const int every = std::max(1, int(p.num("every", 1)));
    for (size_t f = 0; f < st.traj->frames(); f += size_t(every)) {
      if (int(f) == st.frame) continue;
      System o = st.traj->frame(f);
      if (o.atoms.size() != n) continue;
      auto oa = [&](uint32_t i) { return ea == 0 || o.atoms[i].element == ea; };
      auto ob = [&](uint32_t i) { return eb == 0 || o.atoms[i].element == eb; };
      for_pairs(o, rmax, [&](uint32_t i, uint32_t j, double r, const Vec3&) {
        if (inter && mol[i] == mol[j]) return;
        const int k = std::min(bins - 1, int(r / dr));
        hist[size_t(k)] += (oa(i) && ob(j) ? 1 : 0) + (oa(j) && ob(i) ? 1 : 0);
      });
      ++frames;
    }
    for (auto& h : hist) h /= frames;
  }
  st.set_attribute("CoordinationAnalysis.frames", double(frames));
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
  double mx = 0;
  for (uint32_t i = 0; i < n; ++i) if (is_a(i)) mx = std::max(mx, cn[i]);
  st.set_attribute("CoordinationAnalysis.max", mx);
  out.summary = fmt("mean %.2f", mean) + fmt(" · max %.0f", mx) + " within " + fmt("%.2f Å", rc) + " · g(r) to " + fmt("%.3g Å", rmax) +
                (frames > 1 ? " over " + std::to_string(frames) + " frames" : "");
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

// Positions folded back by whole lattice vectors: the shift of each atom as integer images (for the tables)
std::array<double, 3> lattice_shift(const Cell& c, const Vec3& from, const Vec3& to) {
  const Vec3 f = c.to_fractional(to) - c.to_fractional(from);
  return {std::round(f[0]), std::round(f[1]), std::round(f[2])};
}

// bonds still longer than half the narrowest cell width (a molecule left split across a face)
size_t split_bonds(const System& s) {
  const double half = 0.5 * std::min({norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)});
  size_t n = 0;
  for (const auto& b : s.bonds) n += norm(s.atoms[b.i].pos - s.atoms[b.j].pos) > half ? 1 : 0;
  return n;
}

// mode molecules: each molecule made whole and moved by lattice vectors so its centre of mass lies in the cell (no
// bond is cut, as for pictures and for engines that want whole molecules); atoms: every atom folded in
void step_wrap_molecules(PipelineState& st, StepStatus& out) {
  System& s = st.system;
  if (!s.unwrapped) make_molecules_whole(s);
  int nm = 0;
  const auto mol = s.molecules(&nm);
  std::vector<Vec3> com(size_t(nm), Vec3{0, 0, 0});
  std::vector<double> mass(size_t(nm), 0);
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const double m = s.mass_of(s.atoms[i]);
    com[size_t(mol[i])] = com[size_t(mol[i])] + s.atoms[i].pos * m;
    mass[size_t(mol[i])] += m;
  }
  std::vector<Vec3> shift(size_t(nm), Vec3{0, 0, 0});
  int moved = 0;
  for (int m = 0; m < nm; ++m) {
    const Vec3 c = mass[size_t(m)] > 0 ? com[size_t(m)] * (1.0 / mass[size_t(m)]) : Vec3{0, 0, 0};
    shift[size_t(m)] = s.cell.wrap(c) - c;
    moved += norm(shift[size_t(m)]) > 1e-9;
  }
  DataTable t;
  t.name = "outside";
  t.title = "Image shifts · molecules folded in by their centre of mass";
  t.columns = {"Molecule", "Atoms", "na", "nb", "nc"};
  std::vector<int> count(size_t(nm), 0);
  for (int m : mol) ++count[size_t(m)];
  for (int m = 0; m < nm && t.rows.size() < 200; ++m)
    if (norm(shift[size_t(m)]) > 1e-9) {
      const auto n = lattice_shift(s.cell, {0, 0, 0}, shift[size_t(m)]);
      t.rows.push_back({double(m + 1), double(count[size_t(m)]), n[0], n[1], n[2]});
    }
  size_t outside = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    s.atoms[i].pos = s.atoms[i].pos + shift[size_t(mol[i])];
    outside += norm(s.cell.wrap(s.atoms[i].pos) - s.atoms[i].pos) > 1e-9;
  }
  s.unwrapped = false;   // folded by molecule: follow atoms between frames by minimum image
  st.tables.push_back(std::move(t));
  st.set_attribute("Wrap.molecules_moved", double(moved));
  st.set_attribute("Wrap.atoms_outside", double(outside));
  st.set_attribute("Wrap.bonds_across_faces", double(split_bonds(s)));
  out.summary = std::to_string(moved) + " of " + std::to_string(nm) + " molecules folded in whole · " + std::to_string(outside) + " atoms stick out of the cell";
}

void step_wrap(PipelineState& st, const Json& p, StepStatus& out) {
  if (!st.system.cell.valid()) { out.level = "warning"; out.summary = "no cell: nothing to wrap"; return; }
  const std::string mode = p.text("mode", "atoms");
  if (mode == "molecules") return step_wrap_molecules(st, out);
  if (mode != "atoms") throw std::invalid_argument("wrap mode is atoms or molecules");
  System& s = st.system;
  const auto mol = s.molecules();
  std::set<int> crossing;
  DataTable t;   // the first atoms outside the cell: where they are, their image, where they fold to
  t.name = "outside";
  t.title = "Image flags · atoms outside the cell";
  t.columns = {"Identifier", "Molecule", "x (Å)", "y (Å)", "z (Å)", "ix", "iy", "iz", "wrapped x", "wrapped y", "wrapped z"};
  size_t k = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    auto& a = s.atoms[i];
    const Vec3 w = s.cell.wrap(a.pos);
    if (norm(w - a.pos) > 1e-9) {
      ++k;
      crossing.insert(mol[i]);
      if (t.rows.size() < 200) {
        const Vec3 f = s.cell.to_fractional(a.pos);
        t.rows.push_back({double(a.id), double(a.mol ? a.mol : mol[i] + 1), a.pos[0], a.pos[1], a.pos[2], std::floor(f[0]), std::floor(f[1]), std::floor(f[2]), w[0], w[1], w[2]});
      }
    }
    a.pos = w;
  }
  // bonds that now reach across a face (drawn hidden rather than as long sticks)
  const double half = 0.5 * std::min({norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)});
  size_t across = 0;
  for (const auto& b : s.bonds) across += norm(s.atoms[b.i].pos - s.atoms[b.j].pos) > half ? 1 : 0;
  s.unwrapped = false;
  st.tables.push_back(std::move(t));
  st.set_attribute("Wrap.atoms_outside", double(k));
  st.set_attribute("Wrap.molecules_crossing", double(crossing.size()));
  st.set_attribute("Wrap.bonds_across_faces", double(across));
  int nm = 0;
  s.molecules(&nm);
  out.summary = std::to_string(k) + " atoms outside folded in · " + std::to_string(crossing.size()) + " of " + std::to_string(nm) + " molecules cross faces · " +
                std::to_string(across) + " bonds across faces";
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
  // stacked by a category (Type, Molecule, Element …): one count column per value, labelled by the type's name
  const std::string by = p.text("stack_by", "");
  std::vector<double> cat;
  std::map<double, std::vector<double>> per;
  if (!by.empty() && by != "none" && property_values(st, by, cat)) {
    for (size_t i = 0; i < v.size(); ++i) {
      if ((only_sel && !st.selected[i]) || !std::isfinite(v[i]) || v[i] < lo || v[i] > hi) continue;
      auto& col = per[cat[i]];
      col.resize(size_t(bins), 0);
      col[size_t(std::min(bins - 1, int((v[i] - lo) / w)))] += 1;
    }
    for (const auto& [c, col] : per) {
      std::string name = by + " " + fmt("%g", c);
      if (by == "Type")
        for (const auto& ti : st.system.types) if (ti.type == int(c) && !ti.label.empty()) name = ti.label;
      if (by == "Element") name = element(int(c)).symbol;
      t.columns.push_back(name);
    }
  }
  for (int k = 0; k < bins; ++k) {
    std::vector<double> row = {lo + (k + 0.5) * w, h[size_t(k)]};
    for (const auto& [c, col] : per) row.push_back(col[size_t(k)]);
    t.rows.push_back(std::move(row));
  }
  st.tables.push_back(std::move(t));
  out.summary = std::to_string(used) + " values in " + std::to_string(bins) + " bins" + (per.empty() ? "" : " · stacked by " + by + " (" + std::to_string(per.size()) + ")");
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
  const bool only_sel = flag(p, "only_selected", false), inter = flag(p, "inter_only", false);
  const auto mol = inter ? s.molecules() : std::vector<int>{};
  std::vector<Bond> made;
  const std::string mode = p.text("mode", "perceive");
  if (mode == "cutoff" || mode == "pairs") {
    // one cutoff, or a cutoff per element pair ("C-C": 1.70, "C-H": 1.25; 0 or absent: never bonded)
    std::map<std::pair<int, int>, double> pair_rc;
    double rmax = p.num("cutoff", 1.6);
    auto add_pair = [&](const std::string& k, double v) {
      const auto dash = k.find_first_of("-–");
      if (dash == std::string::npos || v <= 0) return;
      const size_t skip = k.compare(dash, 3, "–") == 0 ? 3 : 1;
      const int a = element_from_symbol(k.substr(0, dash)), b = element_from_symbol(k.substr(dash + skip));
      if (!a || !b) return;
      pair_rc[{std::min(a, b), std::max(a, b)}] = v;
    };
    if (mode == "pairs" && p.has("pairs")) {
      if (p["pairs"].is_object()) {
        for (const auto& [k, v] : p["pairs"].members()) if (v.is_number()) add_pair(k, v.number());
      } else if (p["pairs"].is_string()) {   // "C-C 1.70, C-H 1.25"
        std::string t = p["pairs"].str();
        for (char& c : t) if (c == ',' || c == ';' || c == ':' || c == '=') c = ' ';
        std::istringstream in(t);
        std::string k;
        double v;
        while (in >> k >> v) add_pair(k, v);
      }
      rmax = 0;
      for (const auto& [k, v] : pair_rc) rmax = std::max(rmax, v);
    }
    if (rmax > 0)
      for_pairs(s, rmax, [&](uint32_t i, uint32_t j, double r, const Vec3&) {
        if (r < 0.4 || (only_sel && (!st.selected[i] || !st.selected[j])) || (inter && mol[i] == mol[j])) return;
        if (mode == "pairs") {
          const int a = s.atoms[i].element, b = s.atoms[j].element;
          auto it = pair_rc.find({std::min(a, b), std::max(a, b)});
          if (it == pair_rc.end() || r > it->second) return;
        }
        made.push_back({i, j, 0});
      });
  } else {
    BondOptions o;
    o.tolerance = p.num("tolerance", 0.45);
    for (const auto& b : perceive_bonds(s, o))
      if ((!only_sel || (st.selected[b.i] && st.selected[b.j])) && (!inter || mol[b.i] != mol[b.j])) made.push_back(b);
  }
  auto key = [](const Bond& b) { return std::pair{std::min(b.i, b.j), std::max(b.i, b.j)}; };
  // a file with bonds keeps them: the new bonds are a check against them, listed, never silently merged
  const bool keep = flag(p, "keep_file", true) && before > 0 && !flag(p, "replace", false);
  if (keep) {
    std::set<std::pair<uint32_t, uint32_t>> file, cut;
    for (const auto& b : s.bonds) file.insert(key(b));
    for (const auto& b : made) cut.insert(key(b));
    size_t both = 0, cut_only = 0, file_only = 0;
    for (const auto& k : cut) (file.count(k) ? both : cut_only)++;
    for (const auto& k : file) file_only += cut.count(k) ? 0 : 1;
    double lo = 1e300, hi = 0;
    for (const auto& b : s.bonds) {
      Vec3 d = s.atoms[b.j].pos - s.atoms[b.i].pos;
      if (s.cell.valid()) d = s.cell.minimum_image(d);
      lo = std::min(lo, norm(d)), hi = std::max(hi, norm(d));
    }
    for (const auto& k : cut)   // spurious cutoff bonds in red
      if (!file.count(k)) {
        Vec3 d = s.atoms[k.second].pos - s.atoms[k.first].pos;
        if (s.cell.valid()) d = s.cell.minimum_image(d);
        st.segments.push_back({s.atoms[k.first].pos, s.atoms[k.first].pos + d, 0xE5484D, 0.12, false});
      }
    st.set_attribute("CreateBonds.in_both", double(both));
    st.set_attribute("CreateBonds.cutoff_only", double(cut_only));
    st.set_attribute("CreateBonds.topology_only", double(file_only));
    st.set_attribute("CreateBonds.topology_bonds", double(file.size()));
    if (!s.bonds.empty()) st.set_attribute("CreateBonds.length_min", lo), st.set_attribute("CreateBonds.length_max", hi);
    out.summary = std::to_string(cut.size()) + " bonds by " + mode + " · " + std::to_string(file.size()) + " in the file · " + std::to_string(cut_only) + " extra · " +
                  std::to_string(file_only) + " missing";
    if (cut_only || file_only) out.level = "warning";
    return;
  }
  if (flag(p, "replace", false)) s.bonds.clear();
  std::set<std::pair<uint32_t, uint32_t>> have;
  for (const auto& b : s.bonds) have.insert(key(b));
  for (const auto& b : made)
    if (have.insert(key(b)).second) s.bonds.push_back(b);
  out.summary = std::to_string(s.bonds.size()) + " bonds (" + (s.bonds.size() >= before ? "+" : "") + std::to_string(long(s.bonds.size()) - long(before)) + ")";
}

std::vector<Vec3> frame_positions(const PipelineState& st, size_t k);

// method bonds: molecules made whole by following their bonds (minimum image along each bond); images: the file's
// image flags (LAMMPS ix iy iz), as read; nojump: each atom followed from frame 0 (made whole) through the frames by
// its shortest step, so nothing jumps by a box length (MSD and diffusion)
void step_unwrap(PipelineState& st, const Json& p, StepStatus& out) {
  if (!st.system.cell.valid()) { out.level = "warning"; out.summary = "no cell: nothing to unwrap"; return; }
  System& s = st.system;
  const std::string method = p.text("method", "bonds");
  const auto before = s.atoms;
  std::string how;
  if (method == "bonds") {
    make_molecules_whole(s);
    how = "molecules whole along their bonds";
  } else if (method == "images") {
    size_t flagged = 0;
    for (const auto& a : s.atoms) flagged += a.image[0] || a.image[1] || a.image[2];
    if (!s.unwrapped)
      for (auto& a : s.atoms) a.pos = a.pos + s.cell.a * a.image[0] + s.cell.b * a.image[1] + s.cell.c * a.image[2];
    how = flagged ? std::to_string(flagged) + " atoms with image flags" + (s.unwrapped ? " (applied when the file was read)" : "")
                  : "no image flags in the file";
    if (!flagged && !s.unwrapped) {   // the fallback: follow the bonds
      make_molecules_whole(s);
      how += ": molecules whole along their bonds instead";
    }
  } else if (method == "nojump") {
    if (!st.traj || st.traj->frames() < 2) {
      make_molecules_whole(s);
      how = "one frame: molecules whole along their bonds";
    } else {
      auto raw = [&](size_t k) {
        const System f = st.traj->frame(k);
        std::vector<Vec3> r(st.origin.size());
        for (size_t i = 0; i < r.size(); ++i) {
          const int o = st.origin[i];
          r[i] = o >= 0 && size_t(o) < f.atoms.size() ? f.atoms[size_t(o)].pos : s.atoms[i].pos;
        }
        return r;
      };
      std::vector<Vec3> x = frame_positions(st, 0), prev = raw(0);
      for (size_t k = 1; k <= size_t(std::max(0, st.frame)); ++k) {
        auto cur = raw(k);
        for (size_t i = 0; i < x.size(); ++i) x[i] = x[i] + s.cell.minimum_image(cur[i] - prev[i]);
        prev = std::move(cur);
      }
      for (size_t i = 0; i < x.size(); ++i) s.atoms[i].pos = x[i];
      how = "atoms followed through " + std::to_string(st.frame + 1) + " frames without jumps";
    }
  } else throw std::invalid_argument("unwrap method is bonds, images or nojump");
  s.unwrapped = true;
  DataTable t;   // what moved, by how many cells
  t.name = "images";
  t.title = "Image shifts · atoms moved by unwrapping";
  t.columns = {"Identifier", "Molecule", "na", "nb", "nc"};
  size_t k = 0;
  for (size_t i = 0; i < before.size(); ++i)
    if (norm(before[i].pos - s.atoms[i].pos) > 1e-6) {
      ++k;
      if (t.rows.size() < 200) {
        const auto n = lattice_shift(s.cell, before[i].pos, s.atoms[i].pos);
        t.rows.push_back({double(s.atoms[i].id), double(s.atoms[i].mol), n[0], n[1], n[2]});
      }
    }
  const size_t split = split_bonds(s);
  st.tables.push_back(std::move(t));
  st.set_attribute("Unwrap.atoms_moved", double(k));
  st.set_attribute("Unwrap.bonds_split", double(split));
  out.summary = std::to_string(k) + " moved · " + how;
  if (split) {
    out.level = "warning";
    out.summary += " · " + std::to_string(split) + " bonds still span half the cell" + (method == "images" ? " (inconsistent image flags: unwrap along bonds)" : "");
  }
}

Vec3 unit_or_zero(const Vec3& v) { const double n = norm(v); return n > 0 ? v * (1 / n) : v; }

// Largest eigenvalue and its vector of a symmetric 3 × 3 matrix (Jacobi)
std::pair<double, Vec3> largest_eigen(double A[3][3]) {
  double V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (int sweep = 0; sweep < 50; ++sweep) {
    int p = 0, q = 1;
    for (int i = 0; i < 3; ++i)
      for (int j = i + 1; j < 3; ++j)
        if (std::fabs(A[i][j]) > std::fabs(A[p][q])) p = i, q = j;
    if (std::fabs(A[p][q]) < 1e-14) break;
    const double th = 0.5 * std::atan2(2 * A[p][q], A[q][q] - A[p][p]), c = std::cos(th), sn = std::sin(th);
    for (int k = 0; k < 3; ++k) {
      const double akp = A[k][p], akq = A[k][q];
      A[k][p] = c * akp - sn * akq, A[k][q] = sn * akp + c * akq;
    }
    for (int k = 0; k < 3; ++k) {
      const double apk = A[p][k], aqk = A[q][k];
      A[p][k] = c * apk - sn * aqk, A[q][k] = sn * apk + c * aqk;
    }
    for (int k = 0; k < 3; ++k) {
      const double vkp = V[k][p], vkq = V[k][q];
      V[k][p] = c * vkp - sn * vkq, V[k][q] = sn * vkp + c * vkq;
    }
  }
  int im = 0;
  for (int i = 1; i < 3; ++i) if (A[i][i] > A[im][im]) im = i;
  return {A[im][im], Vec3{V[0][im], V[1][im], V[2][im]}};
}

// Freeze property: the property as it was at a reference frame (the steps below this one run on that frame), matched
// to this frame's particles by identifier — colour by starting height to see flow, or keep frame-0 clusters
void step_freeze(PipelineState& st, const Json& p, StepStatus& out) {
  const std::string prop = p.text("property", "Position.Z");
  const std::string name = p.text("output", prop + " frozen");
  const int nframes = st.traj ? int(st.traj->frames()) : 1;
  const int ref = std::clamp(int(p.num("frame", 0)), 0, std::max(0, nframes - 1));
  std::vector<double> refv;
  std::vector<int64_t> refid;
  if (ref == st.frame || !st.traj) {   // this frame is the reference
    if (!property_values(st, prop, refv)) throw std::invalid_argument("no property " + prop);
    for (const auto& a : st.system.atoms) refid.push_back(a.id);
  } else {
    Pipeline up;   // the steps below this one, as they run on the reference frame
    if (st.pipeline)
      for (size_t k = st.step_index + 1; k < st.pipeline->steps.size(); ++k) up.steps.push_back(st.pipeline->steps[k]);
    const PipelineState rs = run_pipeline(st.traj->frame(size_t(ref)), up, ref, 0, st.traj);
    if (!property_values(rs, prop, refv)) throw std::invalid_argument("no property " + prop + " at frame " + std::to_string(ref));
    for (const auto& a : rs.system.atoms) refid.push_back(a.id);
  }
  std::unordered_map<int64_t, double> by_id;
  for (size_t i = 0; i < refid.size(); ++i) by_id[refid[i]] = refv[i];
  auto& dst = st.props[name];
  dst.assign(st.system.atoms.size(), std::nan(""));
  size_t found = 0;
  for (size_t i = 0; i < dst.size(); ++i)
    if (auto it = by_id.find(st.system.atoms[i].id); it != by_id.end()) dst[i] = it->second, ++found;
  out.summary = prop + " at frame " + std::to_string(ref) + " → " + name + (found < dst.size() ? " · " + std::to_string(dst.size() - found) + " particles not there" : "");
  if (found < dst.size()) out.level = "warning";
}

// Chain orientation per atom (design/boards/PipelineSteps "Chain orientation", "Crystallinity (polymer)"): the backbone
// chord through each backbone atom (its neighbours i−1 → i+1), P₂ of its angle to the axis (the director, or x, y, z),
// and crystalline where at least `neighbours` other chords with midpoints within `radius` are aligned within `angle`
// degrees; hydrogens and side groups take their backbone atom's values
void step_orientation(PipelineState& st, const Json& p, StepStatus& out) {
  System whole = st.system;
  if (!whole.unwrapped && whole.cell.valid()) make_molecules_whole(whole);
  const auto bb = backbones(whole, 3);
  const size_t n = whole.atoms.size();
  std::vector<Vec3> u, mid;
  std::vector<uint32_t> at;   // the atom each chord is centred on
  for (const auto& b : bb)
    for (size_t i = 1; i + 1 < b.size(); ++i) {
      const Vec3 d = whole.atoms[b[i + 1]].pos - whole.atoms[b[i - 1]].pos;
      const double l = norm(d);
      if (l < 1e-9) continue;
      u.push_back(d * (1 / l));
      mid.push_back(whole.atoms[b[i]].pos);
      at.push_back(b[i]);
    }
  if (u.empty()) { out.level = "warning"; out.summary = "no chain backbones of three or more heavy atoms"; return; }
  double Q[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (const Vec3& v : u)
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) Q[a][b] += (1.5 * v[a] * v[b] - (a == b ? 0.5 : 0)) / double(u.size());
  auto [S, director] = largest_eigen(Q);
  const std::string axis = p.text("axis", "director");
  Vec3 ax = director;
  if (axis == "x") ax = {1, 0, 0};
  else if (axis == "y") ax = {0, 1, 0};
  else if (axis == "z") ax = {0, 0, 1};
  else if (axis != "director") throw std::invalid_argument("axis is director, x, y or z");
  const double radius = std::max(1.0, p.num("radius", 5.0)), cosa = std::cos(std::clamp(p.num("angle", 10.0), 0.0, 90.0) * M_PI / 180);
  const int need = std::max(1, int(p.num("neighbours", 8)));
  // chords binned by midpoint (fractional in a cell, else a box around them) for the neighbour search
  const Cell& c = whole.cell;
  const bool per = c.valid();
  Vec3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
  for (const Vec3& m : mid)
    for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], m[k]), hi[k] = std::max(hi[k], m[k]);
  std::array<int, 3> nb{1, 1, 1};
  const double width[3] = {per ? 1.0 / std::fabs(dot(c.a, unit_or_zero(cross(c.b, c.c)))) : 0, per ? 1.0 / std::fabs(dot(c.b, unit_or_zero(cross(c.c, c.a)))) : 0,
                           per ? 1.0 / std::fabs(dot(c.c, unit_or_zero(cross(c.a, c.b)))) : 0};
  for (int k = 0; k < 3; ++k) nb[k] = std::max(1, int((per ? 1.0 / width[k] : hi[k] - lo[k] + 1e-6) / radius));
  auto bin = [&](const Vec3& m, int k) {
    double f = per ? c.to_fractional(m)[k] : (m[k] - lo[k]) / (hi[k] - lo[k] + 1e-6);
    if (per) f -= std::floor(f);
    return std::clamp(int(f * nb[k]), 0, nb[k] - 1);
  };
  std::vector<std::vector<uint32_t>> grid(size_t(nb[0]) * nb[1] * nb[2]);
  std::vector<std::array<int, 3>> where(u.size());
  for (size_t i = 0; i < u.size(); ++i) {
    where[i] = {bin(mid[i], 0), bin(mid[i], 1), bin(mid[i], 2)};
    grid[(size_t(where[i][0]) * nb[1] + where[i][1]) * nb[2] + where[i][2]].push_back(uint32_t(i));
  }
  std::vector<double> p2(n, std::nan("")), cr(n, 0), chordp2(u.size()), chordcr(u.size());
  size_t crystalline = 0;
  double f_axis = 0;
  for (size_t i = 0; i < u.size(); ++i) {
    const double ct = dot(u[i], ax);
    chordp2[i] = 1.5 * ct * ct - 0.5;
    f_axis += chordp2[i];
    int aligned = 0;
    std::set<size_t> seen;
    for (int dx = -1; dx <= 1 && aligned < need; ++dx)
      for (int dy = -1; dy <= 1 && aligned < need; ++dy)
        for (int dz = -1; dz <= 1 && aligned < need; ++dz) {
          int bx = where[i][0] + dx, by = where[i][1] + dy, bz = where[i][2] + dz;
          if (per) bx = (bx + nb[0]) % nb[0], by = (by + nb[1]) % nb[1], bz = (bz + nb[2]) % nb[2];
          else if (bx < 0 || by < 0 || bz < 0 || bx >= nb[0] || by >= nb[1] || bz >= nb[2]) continue;
          const size_t cellid = (size_t(bx) * nb[1] + by) * nb[2] + bz;
          if (!seen.insert(cellid).second) continue;
          for (uint32_t j : grid[cellid]) {
            if (j == i) continue;
            const Vec3 d = per ? c.minimum_image(mid[j] - mid[i]) : mid[j] - mid[i];
            if (dot(d, d) <= radius * radius && std::fabs(dot(u[i], u[j])) >= cosa && ++aligned >= need) break;
          }
        }
    chordcr[i] = aligned >= need;
    crystalline += aligned >= need;
    p2[at[i]] = chordp2[i];
    cr[at[i]] = chordcr[i];
  }
  // the rest of each molecule: its nearest backbone atom's values (along bonds)
  const auto nbrs = whole.neighbours();
  std::vector<uint32_t> q;
  std::vector<char> done(n, 0);
  for (size_t i = 0; i < n; ++i) if (!std::isnan(p2[i])) done[i] = 1, q.push_back(uint32_t(i));
  for (size_t h = 0; h < q.size(); ++h)
    for (uint32_t w : nbrs[q[h]])
      if (!done[w]) done[w] = 1, p2[w] = p2[q[h]], cr[w] = cr[q[h]], q.push_back(w);
  for (size_t i = 0; i < n; ++i) if (!done[i]) p2[i] = 0;
  st.props["Orientation"] = p2;
  st.props["Crystalline"] = cr;
  const double frac = double(crystalline) / double(u.size());
  st.set_attribute("Orientation.S", S);
  st.set_attribute("Orientation.P2_axis", f_axis / double(u.size()));
  st.set_attribute("Orientation.director_x", director[0]);
  st.set_attribute("Orientation.director_y", director[1]);
  st.set_attribute("Orientation.director_z", director[2]);
  st.set_attribute("Crystallinity.fraction", frac);
  char buf[200];
  std::snprintf(buf, sizeof buf, "S %.3f · ⟨P₂⟩ along %s %.3f · %.1f %% of %zu chords crystalline", S, axis.c_str(), f_axis / double(u.size()), 100 * frac, u.size());
  out.summary = buf;
}

// Affine transformation (shear, strain, rotation): positions and/or the cell mapped by the 3 × 3 matrix (row-major,
// x' = M x) and shifted by the translation; `strain` (exx eyy ezz) is the diagonal shortcut, about the cell centre
void step_affine(PipelineState& st, const Json& p, StepStatus& out) {
  double M[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  if (p.has("matrix") && p["matrix"].is_array()) {
    if (p["matrix"].size() != 9) throw std::invalid_argument("the matrix needs nine numbers (row by row)");
    for (int k = 0; k < 9; ++k) M[k / 3][k % 3] = p["matrix"][size_t(k)].number();
  }
  if (p.has("strain") && p["strain"].is_array() && p["strain"].size() == 3)
    for (int k = 0; k < 3; ++k) M[k][k] *= 1 + p["strain"][size_t(k)].number();
  Vec3 t{0, 0, 0};
  if (p.has("translation") && p["translation"].is_array() && p["translation"].size() == 3) t = {p["translation"][0].number(), p["translation"][1].number(), p["translation"][2].number()};
  const std::string target = p.text("target", "all");
  if (target != "all" && target != "particles" && target != "cell") throw std::invalid_argument("target is all, particles or cell");
  System& s = st.system;
  // about the cell's centre (or the particles' centre), so a strain keeps the sample in place
  Vec3 o{0, 0, 0};
  if (s.cell.valid()) o = s.cell.origin + (s.cell.a + s.cell.b + s.cell.c) * 0.5;
  else if (!s.atoms.empty()) { for (const auto& a : s.atoms) o = o + a.pos; o = o * (1.0 / double(s.atoms.size())); }
  auto apply = [&](const Vec3& v) { return Vec3{M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2], M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2], M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2]}; };
  const bool only_sel = flag(p, "only_selected", false);
  size_t moved = 0;
  if (target != "cell")
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      if (only_sel && !st.selected[i]) continue;
      s.atoms[i].pos = o + apply(s.atoms[i].pos - o) + t;
      ++moved;
    }
  const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  if (target != "particles" && s.cell.valid() && !only_sel) {
    const Vec3 corner = o + apply(s.cell.origin - o) + t;
    s.cell.a = apply(s.cell.a), s.cell.b = apply(s.cell.b), s.cell.c = apply(s.cell.c);
    s.cell.origin = corner;
  }
  st.set_attribute("AffineTransformation.volume_ratio", det);
  char buf[160];
  std::snprintf(buf, sizeof buf, "%zu particles%s · volume × %.4f", moved, target != "particles" && s.cell.valid() && !only_sel ? " and the cell" : "", det);
  out.summary = buf;
}

// The k nearest neighbours of every particle (vectors from it, shortest first), periodic in a cell: particles binned
// by a radius from the number density, widened for a particle that finds fewer than k
std::vector<std::vector<Vec3>> nearest_neighbours(const System& s, int k) {
  const size_t n = s.atoms.size();
  std::vector<std::vector<Vec3>> out(n);
  if (n < 2) return out;
  const Cell& c = s.cell;
  const bool per = c.valid();
  Vec3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
  for (const auto& a : s.atoms)
    for (int d = 0; d < 3; ++d) lo[d] = std::min(lo[d], a.pos[d]), hi[d] = std::max(hi[d], a.pos[d]);
  double vol = per ? c.volume() : std::max(1.0, (hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1));
  const double rk = std::cbrt(3.0 * k * vol / (4 * M_PI * double(n)));
  double width[3];
  for (int d = 0; d < 3; ++d) width[d] = per ? 1.0 / norm(cross(d == 0 ? c.b : d == 1 ? c.c : c.a, d == 0 ? c.c : d == 1 ? c.a : c.b)) * c.volume() : hi[d] - lo[d] + 1e-6;
  const double bin = std::max(1.0, 0.8 * rk);
  int nb[3];
  for (int d = 0; d < 3; ++d) nb[d] = std::max(1, int(width[d] / bin));
  auto frac = [&](const Vec3& x, int d) {
    double f = per ? c.to_fractional(x)[d] : (x[d] - lo[d]) / width[d];
    if (per) f -= std::floor(f);
    return std::clamp(int(f * nb[d]), 0, nb[d] - 1);
  };
  std::vector<std::vector<uint32_t>> grid(size_t(nb[0]) * nb[1] * nb[2]);
  std::vector<std::array<int, 3>> where(n);
  for (size_t i = 0; i < n; ++i) {
    where[i] = {frac(s.atoms[i].pos, 0), frac(s.atoms[i].pos, 1), frac(s.atoms[i].pos, 2)};
    grid[(size_t(where[i][0]) * nb[1] + where[i][1]) * nb[2] + where[i][2]].push_back(uint32_t(i));
  }
  for (size_t i = 0; i < n; ++i) {
    std::vector<std::pair<double, Vec3>> cand;
    for (int reach = 1;; ++reach) {
      cand.clear();
      std::set<size_t> seen;
      const bool all = 2 * reach + 1 >= nb[0] && 2 * reach + 1 >= nb[1] && 2 * reach + 1 >= nb[2];
      for (int dx = -reach; dx <= reach; ++dx)
        for (int dy = -reach; dy <= reach; ++dy)
          for (int dz = -reach; dz <= reach; ++dz) {
            int b[3] = {where[i][0] + dx, where[i][1] + dy, where[i][2] + dz};
            bool skip = false;
            for (int d = 0; d < 3; ++d) {
              if (per) b[d] = ((b[d] % nb[d]) + nb[d]) % nb[d];
              else if (b[d] < 0 || b[d] >= nb[d]) skip = true;
            }
            if (skip) continue;
            const size_t id = (size_t(b[0]) * nb[1] + b[1]) * nb[2] + b[2];
            if (!seen.insert(id).second) continue;
            for (uint32_t j : grid[id]) {
              if (j == i) continue;
              const Vec3 v = per ? c.minimum_image(s.atoms[j].pos - s.atoms[i].pos) : s.atoms[j].pos - s.atoms[i].pos;
              cand.push_back({dot(v, v), v});
            }
          }
      // enough, and the k-th is closer than the searched shell's inner reach (else a nearer one may lie beyond)
      if (int(cand.size()) >= k) {
        std::partial_sort(cand.begin(), cand.begin() + k, cand.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        double shell = 1e30;
        for (int d = 0; d < 3; ++d) shell = std::min(shell, reach * width[d] / nb[d]);
        if (all || std::sqrt(cand[size_t(k - 1)].first) <= shell) break;
      } else if (all) break;
    }
    const size_t m = std::min<size_t>(size_t(k), cand.size());
    std::partial_sort(cand.begin(), cand.begin() + long(m), cand.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    for (size_t q = 0; q < m; ++q) out[i].push_back(cand[q].second);
  }
  return out;
}

// One CNA signature per neighbour: common neighbours (within rc of both), bonds among them, the largest connected set
// of those bonds (Honeycutt & Andersen 1987; a-CNA's local cutoff, Stukowski 2012)
std::array<int, 3> cna_signature(const std::vector<Vec3>& nb, size_t j, double rc2) {
  std::vector<size_t> common;
  for (size_t k = 0; k < nb.size(); ++k)
    if (k != j) { const Vec3 d = nb[k] - nb[j]; if (dot(d, d) < rc2) common.push_back(k); }
  std::vector<std::pair<size_t, size_t>> bonds;
  for (size_t a = 0; a < common.size(); ++a)
    for (size_t b = a + 1; b < common.size(); ++b) {
      const Vec3 d = nb[common[a]] - nb[common[b]];
      if (dot(d, d) < rc2) bonds.push_back({a, b});
    }
  // largest set of bonds connected through shared atoms
  std::vector<int> comp(bonds.size(), -1);
  int longest = 0;
  for (size_t b0 = 0; b0 < bonds.size(); ++b0) {
    if (comp[b0] >= 0) continue;
    comp[b0] = int(b0);
    std::vector<size_t> q{b0};
    int size = 0;
    while (!q.empty()) {
      const size_t b = q.back();
      q.pop_back();
      ++size;
      for (size_t o = 0; o < bonds.size(); ++o)
        if (comp[o] < 0 && (bonds[o].first == bonds[b].first || bonds[o].first == bonds[b].second || bonds[o].second == bonds[b].first || bonds[o].second == bonds[b].second))
          comp[o] = int(b0), q.push_back(o);
    }
    longest = std::max(longest, size);
  }
  return {int(common.size()), int(bonds.size()), longest};
}

// Adaptive common neighbour analysis (Stukowski, Modelling Simul. Mater. Sci. Eng. 20, 045021 (2012)): FCC, HCP, BCC,
// icosahedral or other, per particle → Structure Type (0 other, 1 FCC, 2 HCP, 3 BCC, 4 ICO)
void step_cna(PipelineState& st, const Json& p, StepStatus& out) {
  const System& s = st.system;
  const size_t n = s.atoms.size();
  const auto nn = nearest_neighbours(s, 14);
  auto& type = st.props["Structure Type"];
  type.assign(n, 0);
  const bool only_sel = flag(p, "only_selected", false);
  std::array<size_t, 5> count{0, 0, 0, 0, 0};
  const double f = (1 + std::sqrt(2.0)) / 2;
  for (size_t i = 0; i < n; ++i) {
    if (only_sel && !st.selected[i]) continue;
    const auto& all = nn[i];
    int t = 0;
    if (all.size() >= 12) {
      std::vector<Vec3> n12(all.begin(), all.begin() + 12);
      double mean = 0;
      for (const Vec3& v : n12) mean += norm(v);
      const double rc = f * mean / 12;
      int s421 = 0, s422 = 0, s555 = 0;
      for (size_t j = 0; j < 12; ++j) {
        const auto sig = cna_signature(n12, j, rc * rc);
        if (sig == std::array<int, 3>{4, 2, 1}) ++s421;
        else if (sig == std::array<int, 3>{4, 2, 2}) ++s422;
        else if (sig == std::array<int, 3>{5, 5, 5}) ++s555;
      }
      if (s421 == 12) t = 1;
      else if (s421 == 6 && s422 == 6) t = 2;
      else if (s555 == 12) t = 4;
    }
    if (t == 0 && all.size() >= 14) {
      double m8 = 0, m6 = 0;
      for (size_t j = 0; j < 8; ++j) m8 += norm(all[j]);
      for (size_t j = 8; j < 14; ++j) m6 += norm(all[j]);
      const double rc = f * 0.5 * (2 / std::sqrt(3.0) * m8 / 8 + m6 / 6);
      int s444 = 0, s666 = 0;
      for (size_t j = 0; j < 14; ++j) {
        const auto sig = cna_signature(all, j, rc * rc);
        if (sig == std::array<int, 3>{4, 4, 4}) ++s444;
        else if (sig == std::array<int, 3>{6, 6, 6}) ++s666;
      }
      if (s444 == 6 && s666 == 8) t = 3;
    }
    type[i] = t;
    ++count[size_t(t)];
  }
  static const char* names[5] = {"Other", "FCC", "HCP", "BCC", "ICO"};
  DataTable tab;
  tab.name = "structures";
  tab.title = "Common neighbour analysis · particles per structure";
  tab.columns = {"Type", "Count", "Fraction"};
  const double tot = double(std::max<size_t>(1, only_sel ? st.selected_count() : n));
  std::string sum;
  for (int t = 0; t < 5; ++t) {
    tab.rows.push_back({double(t), double(count[size_t(t)]), count[size_t(t)] / tot});
    st.set_attribute(std::string("CommonNeighborAnalysis.counts.") + names[t], double(count[size_t(t)]));
    if (count[size_t(t)]) sum += (sum.empty() ? "" : " · ") + std::string(names[t]) + " " + std::to_string(count[size_t(t)]);
  }
  st.tables.push_back(std::move(tab));
  out.summary = sum.empty() ? "no particles" : sum;
}

// Centrosymmetry parameter (Kelchner, Plimpton & Hamilton, Phys. Rev. B 58, 11085 (1998)), as LAMMPS computes it: the
// sum of the N/2 smallest |r_i + r_j|² over pairs of the N nearest neighbours (0 in a perfect centrosymmetric lattice)
void step_centrosymmetry(PipelineState& st, const Json& p, StepStatus& out) {
  const int N = int(p.num("neighbours", 12));
  if (N < 2 || N > 32 || N % 2) throw std::invalid_argument("the neighbour count is even, 2 to 32 (12 for FCC, 8 for BCC)");
  const size_t n = st.system.atoms.size();
  const auto nn = nearest_neighbours(st.system, N);
  auto& csp = st.props["Centrosymmetry"];
  csp.assign(n, 0);
  double mx = 0, mean = 0;
  for (size_t i = 0; i < n; ++i) {
    const auto& v = nn[i];
    std::vector<double> sums;
    for (size_t a = 0; a < v.size(); ++a)
      for (size_t b = a + 1; b < v.size(); ++b) { const Vec3 q = v[a] + v[b]; sums.push_back(dot(q, q)); }
    const size_t half = std::min(sums.size(), v.size() / 2);
    std::partial_sort(sums.begin(), sums.begin() + long(half), sums.end());
    csp[i] = std::accumulate(sums.begin(), sums.begin() + long(half), 0.0);
    mx = std::max(mx, csp[i]);
    mean += csp[i];
  }
  mean /= double(std::max<size_t>(1, n));
  st.set_attribute("Centrosymmetry.mean", mean);
  st.set_attribute("Centrosymmetry.max", mx);
  char buf[120];
  std::snprintf(buf, sizeof buf, "%d neighbours · mean %.3g Å² · max %.3g Å²", N, mean, mx);
  out.summary = buf;
}

void step_molecule_shape(PipelineState& st, const Json& p, StepStatus& out) {
  System whole = st.system;
  if (!whole.unwrapped && whole.cell.valid()) make_molecules_whole(whole);
  const auto shapes = molecule_shapes(whole);
  const auto mol = whole.molecules();
  const size_t n = whole.atoms.size();
  auto& rg = st.props["MoleculeRg"];
  auto& k2 = st.props["MoleculeKappa2"];
  auto& as = st.props["MoleculeAsphericity"];
  rg.assign(n, 0); k2.assign(n, 0); as.assign(n, 0);
  // end-to-end distance of each molecule's backbone (0 for molecules without one)
  std::vector<double> ree(shapes.size(), 0.0);
  for (const auto& path : backbones(whole, 4)) ree[size_t(mol[path.front()])] = norm(whole.atoms[path.back()].pos - whole.atoms[path.front()].pos);
  DataTable t;
  t.name = "molecules";
  t.title = "Molecule shape";
  t.columns = {"Molecule", "Atoms", "Mass (g/mol)", "Rg (Å)", "Ree (Å)", "Ree²/Rg²", "b (Å²)", "c (Å²)", "κ²", "λ₁", "λ₂", "λ₃", "COM.X", "COM.Y", "COM.Z"};
  double mrg = 0, mk2 = 0;
  const bool glyphs = flag(p, "glyphs", true);
  for (size_t m = 0; m < shapes.size(); ++m) {
    const auto& sh = shapes[m];
    const double b = sh.lambda[2] - 0.5 * (sh.lambda[0] + sh.lambda[1]);   // asphericity
    const double c = sh.lambda[1] - sh.lambda[0];                            // acylindricity
    t.rows.push_back({double(m + 1), double(sh.atoms), sh.mass, sh.rg, ree[m], sh.rg > 0 ? ree[m] * ree[m] / (sh.rg * sh.rg) : 0.0, b, c, sh.kappa2,
                      sh.lambda[0], sh.lambda[1], sh.lambda[2], sh.com[0], sh.com[1], sh.com[2]});
    mrg += sh.rg;
    mk2 += sh.kappa2;
    if (glyphs && sh.atoms > 2)   // principal axes ±√(3λ); the longest drawn thick
      for (int k = 0; k < 3; ++k) {
        const double half = std::sqrt(3 * std::max(0.0, sh.lambda[k]));
        if (half < 0.05) continue;
        st.segments.push_back({sh.com - sh.axis[k] * half, sh.com + sh.axis[k] * half, k == 2 ? 0xF0A83Cu : 0x6CC4D8u, k == 2 ? 0.28 : 0.14, false});
      }
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
        // IUPAC sign: φ = atan2(b̂₂·(n₁×n₂), n₁·n₂)
        dih.push_back(std::atan2(dot(cross(n1, n2), b2 * (1.0 / std::max(1e-12, norm(b2)))), dot(n1, n2)) * 180 / M_PI);
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
  // ranges by kind of bond (C–C, C–H …) and of angle (C–C–C …)
  std::map<std::string, std::pair<double, double>> brange, arange;
  auto sym = [&](uint32_t i) { return std::string(element(s.atoms[i].element).symbol); };
  auto widen = [](std::pair<double, double>& r, double v) { r.first = std::min(r.first, v), r.second = std::max(r.second, v); };
  for (const auto& b : s.bonds) {
    std::string a = sym(b.i), c = sym(b.j);
    if (c < a || (a == "H" && c != "H")) std::swap(a, c);
    if (c == "H" && a != "H") {}   // heavy first: C–H
    auto [it, fresh] = brange.emplace(a + "–" + c, std::pair{1e300, -1e300});
    widen(it->second, norm(sep(b.i, b.j)));
  }
  for (uint32_t j = 0; j < adj.size(); ++j)
    for (size_t x = 0; x < adj[j].size(); ++x)
      for (size_t y = x + 1; y < adj[j].size(); ++y) {
        std::string a = sym(adj[j][x]), c = sym(adj[j][y]);
        if (c < a) std::swap(a, c);
        const Vec3 u = sep(j, adj[j][x]), v = sep(j, adj[j][y]);
        auto [it, fresh] = arange.emplace(a + "–" + sym(j) + "–" + c, std::pair{1e300, -1e300});
        widen(it->second, std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / M_PI);
      }
  DataTable rt;
  rt.name = "ranges";
  rt.title = "Ranges by kind";
  rt.columns = {"#", "min", "max"};
  rt.label_column = "Kind";
  int kk = 0;
  for (const auto& [k, r] : brange) rt.rows.push_back({double(++kk), r.first, r.second}), rt.labels.push_back(k + " (Å)");
  for (const auto& [k, r] : arange) rt.rows.push_back({double(++kk), r.first, r.second}), rt.labels.push_back(k + " (°)");
  st.tables.push_back(std::move(rt));
  // backbone dihedrals: trans |φ| > 120°, gauche± otherwise; the backbone coloured by the state of each dihedral
  int nt = 0, ngp = 0, ngm = 0;
  const bool colour = flag(p, "colour_states", true);
  for (const auto& path : backbones(s, 4))
    for (size_t k = 0; k + 3 < path.size(); ++k) {
      const Vec3 b1 = sep(path[k], path[k + 1]), b2 = sep(path[k + 1], path[k + 2]), b3 = sep(path[k + 2], path[k + 3]);
      const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
      const double phi = std::atan2(dot(cross(n1, n2), b2 * (1.0 / std::max(1e-12, norm(b2)))), dot(n1, n2)) * 180 / M_PI;
      const int state = std::fabs(phi) > 120 ? 0 : phi > 0 ? 1 : 2;
      (state == 0 ? nt : state == 1 ? ngp : ngm)++;
      if (colour) {
        const unsigned c = state == 0 ? 0xF0A83Cu : state == 1 ? 0x6CC4D8u : 0x9B7AD5u;
        st.colour[path[k + 1]] = c, st.colour[path[k + 2]] = c;
      }
    }
  st.set_attribute("Topology.backbone_trans", double(nt));
  st.set_attribute("Topology.backbone_gauche_plus", double(ngp));
  st.set_attribute("Topology.backbone_gauche_minus", double(ngm));
  if (nt + ngp + ngm > 0) st.set_attribute("Topology.trans_fraction", double(nt) / double(nt + ngp + ngm));
  double ml = 0;
  for (double x : len) ml += x;
  st.set_attribute("Topology.mean_bond", len.empty() ? 0 : ml / len.size());
  out.summary = std::to_string(len.size()) + " bonds · " + std::to_string(ang.size()) + " angles · " + std::to_string(dih.size()) + " dihedrals" +
                (nt + ngp + ngm > 0 ? " · backbone t " + std::to_string(nt) + " · g+ " + std::to_string(ngp) + " · g− " + std::to_string(ngm) : "");
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
  std::vector<Vec3> dv(n);
  Vec3 drift{0, 0, 0};
  double mtot = 0;
  for (size_t i = 0; i < n; ++i) {
    Vec3 d = r1[i] - r0[i];
    if (!st.system.unwrapped && st.system.cell.valid() && ref == "previous") d = st.system.cell.minimum_image(d);
    dv[i] = d;
    const double m = st.system.mass_of(st.system.atoms[i]);
    drift = drift + d * m;
    mtot += m;
  }
  drift = mtot > 0 ? drift * (1.0 / mtot) : drift;
  const bool no_drift = flag(p, "subtract_drift", false);   // the system's centre-of-mass motion removed
  double msd = 0, mx = 0;
  for (size_t i = 0; i < n; ++i) {
    const Vec3 d = no_drift ? dv[i] - drift : dv[i];
    dx[i] = d[0]; dy[i] = d[1]; dz[i] = d[2];
    dm[i] = norm(d);
    msd += dot(d, d);
    mx = std::max(mx, dm[i]);
  }
  msd = n ? msd / n : 0;
  st.set_attribute("Displacements.msd", msd);
  st.set_attribute("Displacements.max", mx);
  st.set_attribute("Displacements.drift", norm(drift));
  {   // per molecule: mean and largest atom displacement (heavy atoms), and the centre-of-mass shift
    const auto mol = st.system.molecules();
    std::map<int, std::array<double, 7>> acc;   // sum |d| heavy, n heavy, max, Σm d (3), Σm
    for (size_t i = 0; i < n; ++i) {
      auto& a = acc[mol[i]];
      const Vec3 d = no_drift ? dv[i] - drift : dv[i];
      const double m = st.system.mass_of(st.system.atoms[i]);
      if (st.system.atoms[i].element != 1) a[0] += norm(d), a[1] += 1, a[2] = std::max(a[2], norm(d));
      a[3] += m * d[0], a[4] += m * d[1], a[5] += m * d[2], a[6] += m;
    }
    DataTable t;
    t.name = "displacements";
    t.title = "Displacements per molecule";
    t.columns = {"Molecule", "mean |d| (Å)", "max |d| (Å)", "COM shift (Å)"};
    for (const auto& [m, a] : acc) {
      const Vec3 c = a[6] > 0 ? Vec3{a[3], a[4], a[5]} * (1.0 / a[6]) : Vec3{0, 0, 0};
      t.rows.push_back({double(m + 1), a[1] > 0 ? a[0] / a[1] : 0.0, a[2], norm(c)});
    }
    st.tables.push_back(std::move(t));
  }
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
  // what averaging does to geometry: C–C bond lengths stored and averaged, and how far the atoms moved
  auto cc = [&](const std::function<Vec3(size_t)>& pos) {
    double sm = 0, s2 = 0;
    int c = 0;
    for (const auto& b : st.system.bonds) {
      if (st.system.atoms[b.i].element != 6 || st.system.atoms[b.j].element != 6) continue;
      Vec3 d = pos(b.j) - pos(b.i);
      if (st.system.cell.valid()) d = st.system.cell.minimum_image(d);
      const double l = norm(d);
      sm += l, s2 += l * l, ++c;
    }
    const double m = c ? sm / c : 0;
    return std::pair{m, c ? std::sqrt(std::max(0.0, s2 / c - m * m)) : 0.0};
  };
  const auto before = cc([&](size_t i) { return here[i]; });
  double rms = 0;
  for (size_t i = 0; i < n; ++i) {
    const Vec3 shift = sum[i] * k;
    rms += dot(shift, shift);
    st.system.atoms[i].pos = here[i] + shift;
  }
  rms = n ? std::sqrt(rms / n) : 0;
  const auto after = cc([&](size_t i) { return st.system.atoms[i].pos; });
  st.set_attribute("Smooth.window_frames", double(f1 - f0 + 1));
  st.set_attribute("Smooth.cc_stored_mean", before.first), st.set_attribute("Smooth.cc_stored_sd", before.second);
  st.set_attribute("Smooth.cc_averaged_mean", after.first), st.set_attribute("Smooth.cc_averaged_sd", after.second);
  st.set_attribute("Smooth.rms_shift", rms);
  st.set_attribute("Smoothed", 1);   // smoothed frames are labelled: bond and angle analysis should not read them unasked
  out.summary = "frames " + std::to_string(f0) + "–" + std::to_string(f1) + " averaged · C–C " + fmt("%.3f", before.first) + " ± " + fmt("%.3f", before.second) + " → " +
                fmt("%.3f", after.first) + " ± " + fmt("%.3f Å", after.second) + " · RMS shift " + fmt("%.3f Å", rms);
}

void step_vectors(PipelineState& st, const Json& p, StepStatus& out) {
  const std::string what = p.text("property", "end_to_end");
  const double scale = p.num("scale", 1.0), radius = std::max(0.02, p.num("radius", 0.3));
  System whole = st.system;
  if (!whole.unwrapped && whole.cell.valid()) make_molecules_whole(whole);
  const auto mol = whole.molecules();
  DataTable t;
  t.name = "vectors";
  size_t count = 0;
  double sum = 0, sum2 = 0;
  auto add = [&](const Vec3& a, const Vec3& b, unsigned rgb) {
    st.segments.push_back({a, b, rgb, radius, true});
    const double l = norm(b - a);
    sum += l;
    sum2 += l * l;
    ++count;
  };
  if (what == "end_to_end") {
    t.title = "End-to-end vectors";
    t.columns = {"Molecule", "|R| (Å)", "R.X", "R.Y", "R.Z", "angle to z (°)"};
    const bool flip = flag(p, "flip", false);
    double Q[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, p2z = 0;
    for (const auto& bb : backbones(whole)) {
      if (bb.size() < 2) continue;
      Vec3 a = whole.atoms[bb.front()].pos, b = whole.atoms[bb.back()].pos;
      if (flip) std::swap(a, b);
      const Vec3 r = b - a;
      const int m = mol[bb.front()];
      add(a, b, kCat[m % 10]);
      const double len = norm(r);
      const double cz = len > 0 ? r[2] / len : 0;
      t.rows.push_back({double(m + 1), len, r[0], r[1], r[2], std::acos(std::fabs(cz)) * 180 / M_PI});
      if (len > 0)   // the order tensor of the unit vectors (head–tail symmetric)
        for (int x = 0; x < 3; ++x)
          for (int y = 0; y < 3; ++y) Q[x][y] += 1.5 * r[size_t(x)] * r[size_t(y)] / (len * len) - (x == y ? 0.5 : 0);
      p2z += 1.5 * cz * cz - 0.5;
    }
    st.set_attribute("Vectors.mean_ree", count ? sum / count : 0);
    st.set_attribute("Vectors.mean_ree2", count ? sum2 / count : 0);
    std::string order;
    if (count > 0) {
      for (auto& row : Q) for (double& q : row) q /= double(count);
      double w[3], V[3][3];
      symmetric_eigen3(Q, w, V);
      st.set_attribute("Vectors.order_S", w[2]);
      st.set_attribute("Vectors.director_x", V[0][2]), st.set_attribute("Vectors.director_y", V[1][2]), st.set_attribute("Vectors.director_z", V[2][2]);
      st.set_attribute("Vectors.P2_z", p2z / double(count));
      order = " · S " + fmt("%.3f", w[2]) + " · P₂(z) " + fmt("%.3f", p2z / double(count));
    }
    out.summary = std::to_string(count) + " chains · mean |R| " + fmt("%.2f Å", count ? sum / count : 0) + " · ⟨R²⟩ " + fmt("%.1f Å²", count ? sum2 / count : 0) + order;
  } else if (what == "dipole") {
    t.title = "Molecular dipoles";
    t.columns = {"Molecule", "|μ| (D)", "μ.X (e·Å)", "μ.Y", "μ.Z"};
    const auto shapes = molecule_shapes(whole);
    std::vector<Vec3> mu(shapes.size(), Vec3{0, 0, 0});
    for (size_t i = 0; i < whole.atoms.size(); ++i) mu[size_t(mol[i])] = mu[size_t(mol[i])] + (whole.atoms[i].pos - shapes[size_t(mol[i])].com) * whole.atoms[i].charge;
    for (size_t m = 0; m < shapes.size(); ++m) {
      add(shapes[m].com, shapes[m].com + mu[m] * scale, kCat[m % 10]);
      t.rows.push_back({double(m + 1), norm(mu[m]) * 4.80320, mu[m][0], mu[m][1], mu[m][2]});
    }
    out.summary = std::to_string(shapes.size()) + " molecules · arrows × " + fmt("%g", scale) + " Å per e·Å";
    if (!whole.has_charges) { out.level = "warning"; out.summary += " · no charges in the file"; }
  } else if (what == "displacement" || what == "velocity") {
    const bool disp = what == "displacement";
    std::vector<double> dx, dy, dz;
    if (disp && (!property_values(st, "Displacement.X", dx) || !property_values(st, "Displacement.Y", dy) || !property_values(st, "Displacement.Z", dz)))
      throw std::invalid_argument("add a Displacements step below this one");
    if (!disp && st.system.velocities.size() != st.system.atoms.size()) throw std::invalid_argument("the file has no velocities");
    const bool only_sel = st.selected_count() > 0;
    for (size_t i = 0; i < st.system.atoms.size() && count < 50000; ++i) {
      if (only_sel && !st.selected[i]) continue;
      const Vec3 d = disp ? Vec3{dx[i], dy[i], dz[i]} : st.system.velocities[i] * 1000.0;   // velocities in Å/ps
      if (norm(d) < 1e-9) continue;
      const Vec3 r = st.system.atoms[i].pos;
      add(disp ? r - d * scale : r, disp ? r : r + d * scale, st.colour[i] != kNoColour ? st.colour[i] : 0xF5A524);
    }
    out.summary = std::to_string(count) + " arrows" + (only_sel ? " (selected)" : "") + " · mean " + fmt("%.3g", count ? sum / count / std::max(1e-12, scale) : 0) + (disp ? " Å" : " Å/ps");
  } else {
    throw std::invalid_argument("unknown vector " + what + " (end_to_end, dipole, displacement, velocity)");
  }
  if (!t.rows.empty()) st.tables.push_back(std::move(t));
}

void step_trajectory_lines(PipelineState& st, const Json& p, StepStatus& out) {
  if (!st.traj || st.traj->frames() < 2) { out.level = "warning"; out.summary = "one frame: no paths"; return; }
  const int nf = int(st.traj->frames());
  const int from = std::clamp(int(p.num("from", 0)), 0, nf - 1);
  const int to = std::clamp(int(p.num("to", nf - 1)), from, nf - 1);
  const int stride = std::max(1, int(p.num("stride", std::max(1, (to - from) / 200))));
  const double radius = std::max(0.02, p.num("radius", 0.12));
  const bool centres = p.text("particles", "centres") == "centres";
  // what to trace: molecule centres of mass, or the selected particles (by their frame index)
  const System& s0 = st.system;
  const auto mol = s0.molecules();
  std::vector<int> picks;
  if (!centres)
    for (size_t i = 0; i < s0.atoms.size() && picks.size() < 2000; ++i)
      if (st.selected[i]) picks.push_back(st.origin[i]);
  if (!centres && picks.empty()) throw std::invalid_argument("select the particles to trace (or trace molecule centres)");
  int nm = 0;
  st.traj->topology.molecules(&nm);
  const size_t tracks = centres ? size_t(nm) : picks.size();
  std::vector<Vec3> prev(tracks), path(tracks);
  bool first = true;
  size_t segs = 0;
  std::vector<Vec3> start(tracks);
  std::vector<double> travelled(tracks, 0.0);
  for (int f = from; f <= to; f += stride) {
    System fr = st.traj->frame(size_t(f));
    if (!fr.unwrapped && fr.cell.valid()) make_molecules_whole(fr);
    std::vector<Vec3> now(tracks);
    if (centres) {
      const auto shapes = molecule_shapes(fr);
      for (size_t m = 0; m < tracks && m < shapes.size(); ++m) now[m] = shapes[m].com;
    } else {
      for (size_t k = 0; k < tracks; ++k) now[k] = fr.atoms[size_t(picks[k])].pos;
    }
    for (size_t k = 0; k < tracks; ++k) {
      if (first) { path[k] = now[k]; start[k] = now[k]; }
      else {
        Vec3 d = now[k] - prev[k];
        if (fr.cell.valid() && !fr.unwrapped) d = fr.cell.minimum_image(d);   // follow across the boundary
        const Vec3 next = path[k] + d;
        travelled[k] += norm(d);
        st.segments.push_back({path[k], next, kCat[(centres ? k : size_t(mol[size_t(std::max(0, picks[k]))])) % 10], radius, false});
        path[k] = next;
        ++segs;
      }
      prev[k] = now[k];
    }
    first = false;
  }
  {   // how far each path went, and how far it got: net / path near 0 means the walk turned back on itself
    DataTable t;
    t.name = "paths";
    t.title = centres ? "Centre-of-mass paths" : "Particle paths";
    t.columns = {centres ? "Molecule" : "Particle", "path length (Å)", "net shift (Å)", "net / path"};
    for (size_t k = 0; k < tracks; ++k) {
      const double net = norm(path[k] - start[k]);
      t.rows.push_back({double(centres ? k + 1 : size_t(picks[k] + 1)), travelled[k], net, travelled[k] > 0 ? net / travelled[k] : 0.0});
    }
    st.tables.push_back(std::move(t));
  }
  out.summary = std::to_string(tracks) + (centres ? " molecule centres" : " particles") + " · frames " + std::to_string(from) + "–" + std::to_string(to) +
                (stride > 1 ? " every " + std::to_string(stride) : "") + " · " + std::to_string(segs) + " segments";
}

// Primitive paths (Everaers et al. 2004, see entangle.hpp): each chain pulled tight between its fixed ends without
// crossing another, drawn as a line through the chain; N_e and the tube step as attributes, one row per chain.
void step_primitive_paths(PipelineState& st, const Json& p, StepStatus& out) {
  const double radius = std::max(0.02, p.num("radius", 0.3));
  System whole = st.system;
  if (!whole.unwrapped && whole.cell.valid()) make_molecules_whole(whole);
  std::vector<std::vector<uint32_t>> bb;
  for (auto& b : backbones(whole))
    if (b.size() >= 3) bb.push_back(std::move(b));
  if (bb.empty()) { out.level = "warning"; out.summary = "no chains (backbones of at least 3 heavy atoms)"; return; }
  PrimitivePathOptions po;
  po.max_steps = std::max(1000, int(p.num("max_steps", 200000)));
  const PrimitivePaths pp = primitive_paths(whole, bb, po);
  const EntanglementEstimate e = entanglement_estimate(pp);
  const auto mol = st.system.molecules();
  DataTable t;
  t.name = "primitive_paths";
  t.title = "Primitive paths";
  t.columns = {"Molecule", "backbone bonds", "R (Å)", "L_pp (Å)", "L_pp / R"};
  for (size_t c = 0; c < pp.paths.size(); ++c) {
    // drawn where the chain is drawn: the path moved with the chain's first backbone atom
    const uint32_t a0 = bb[c].front();
    const Vec3 shift = st.system.atoms[a0].pos - whole.atoms[a0].pos;
    const unsigned rgb = kCat[size_t(mol[a0]) % 10];
    const auto& q = pp.paths[c];
    for (size_t k = 1; k < q.size(); ++k)
      if (norm(q[k] - q[k - 1]) > 1e-6) st.segments.push_back({q[k - 1] + shift, q[k] + shift, rgb, radius, false});
    const double R = std::sqrt(pp.r2[c]);
    t.rows.push_back({double(mol[a0] + 1), pp.bonds[c], R, pp.lpp[c], R > 0 ? pp.lpp[c] / R : 1.0});
  }
  st.tables.push_back(std::move(t));
  st.set_attribute("PrimitivePath.chains", e.chains);
  st.set_attribute("PrimitivePath.Lpp", e.lpp);
  st.set_attribute("PrimitivePath.a_pp", e.a_pp);
  if (e.ne_mscoil > 0) {
    st.set_attribute("PrimitivePath.Ne", e.ne_mscoil);
    st.set_attribute("PrimitivePath.Ne_coil", e.ne_coil);
    st.set_attribute("PrimitivePath.Z", e.z);
  }
  char b[200];
  if (e.ne_mscoil > 0)
    std::snprintf(b, sizeof b, "%d chains · ⟨L_pp⟩ %.1f Å · N_e %.0f bonds (Z %.1f per chain) · a_pp %.1f Å%s", e.chains, e.lpp, e.ne_mscoil, e.z, e.a_pp,
                  pp.converged ? "" : " · not converged");
  else
    std::snprintf(b, sizeof b, "%d chains · paths straight: not entangled with each other", e.chains);
  out.summary = b;
  if (!pp.converged) out.level = "warning";
  // the paths run inside the chains: unless asked, the particles are left out so the paths show
  if (!flag(p, "show_chains", false)) compact(st, std::vector<char>(st.system.atoms.size(), 0));
}

// A periodic grid over the cell: n[k] points along each edge, spacing about h.
struct CellGrid {
  Cell cell;
  int n[3] = {1, 1, 1};
  double voxel = 0;
  CellGrid(const Cell& c, double h) : cell(c) {
    for (int k = 0; k < 3; ++k) n[k] = std::clamp(int(std::ceil(norm(k == 0 ? c.a : k == 1 ? c.b : c.c) / h)), 4, 400);
    voxel = c.volume() / (double(n[0]) * n[1] * n[2]);
  }
  size_t size() const { return size_t(n[0]) * n[1] * n[2]; }
  size_t index(int i, int j, int k) const { return (size_t(i) * n[1] + j) * n[2] + k; }
  Vec3 point(int i, int j, int k) const { return cell.origin + cell.a * ((i + 0.5) / n[0]) + cell.b * ((j + 0.5) / n[1]) + cell.c * ((k + 0.5) / n[2]); }
  // every grid point within r of p (minimum image), with its separation vector
  template <class F>
  void around(const Vec3& p, double r, F&& f) const {
    const Vec3 fr = cell.to_fractional(p);
    const double v = cell.volume();
    const double w[3] = {v / norm(cross(cell.b, cell.c)), v / norm(cross(cell.c, cell.a)), v / norm(cross(cell.a, cell.b))};
    int lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
      const double span = r / w[k] * n[k];
      lo[k] = int(std::floor(fr[k] * n[k] - 0.5 - span));
      hi[k] = int(std::ceil(fr[k] * n[k] - 0.5 + span));
      if (hi[k] - lo[k] >= n[k]) { lo[k] = 0; hi[k] = n[k] - 1; }
    }
    const double r2 = r * r;
    for (int i = lo[0]; i <= hi[0]; ++i)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int k = lo[2]; k <= hi[2]; ++k) {
          const int ii = ((i % n[0]) + n[0]) % n[0], jj = ((j % n[1]) + n[1]) % n[1], kk = ((k % n[2]) + n[2]) % n[2];
          const Vec3 d = cell.minimum_image(point(ii, jj, kk) - p);
          const double d2 = dot(d, d);
          if (d2 <= r2) f(index(ii, jj, kk), d2);
        }
  }
};

double vdw_radius(const Atom& a) { return element(a.element).vdw; }

// Accessible grid points for a probe, and their connected voids (face neighbours, periodic), largest first.
std::vector<char> accessible(const System& s, const CellGrid& g, double probe) {
  std::vector<char> ok(g.size(), 1);
  for (const auto& a : s.atoms) g.around(a.pos, vdw_radius(a) + probe, [&](size_t k, double) { ok[k] = 0; });
  return ok;
}
std::vector<std::vector<size_t>> components(const CellGrid& g, const std::vector<char>& ok) {
  std::vector<int> label(g.size(), -1);
  std::vector<std::vector<size_t>> comps;
  std::vector<size_t> stack;
  for (size_t start = 0; start < g.size(); ++start) {
    if (!ok[start] || label[start] >= 0) continue;
    comps.emplace_back();
    auto& c = comps.back();
    label[start] = int(comps.size() - 1);
    stack.assign(1, start);
    while (!stack.empty()) {
      const size_t q = stack.back();
      stack.pop_back();
      c.push_back(q);
      const int i = int(q / (size_t(g.n[1]) * g.n[2])), j = int((q / g.n[2]) % g.n[1]), k = int(q % g.n[2]);
      const int nb[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
      for (const auto& d : nb) {
        const size_t r = g.index((i + d[0] + g.n[0]) % g.n[0], (j + d[1] + g.n[1]) % g.n[1], (k + d[2] + g.n[2]) % g.n[2]);
        if (ok[r] && label[r] < 0) { label[r] = label[q]; stack.push_back(r); }
      }
    }
  }
  std::sort(comps.begin(), comps.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
  return comps;
}

void step_voids(PipelineState& st, const Json& p, StepStatus& out) {
  const System& s = st.system;
  if (!s.cell.valid()) throw std::invalid_argument("voids need a periodic cell");
  const double probe = std::max(0.0, p.num("probe", 1.4));
  const CellGrid g(s.cell, std::clamp(p.num("grid", 0.5), 0.2, 2.0));
  const auto ok = accessible(s, g, probe);
  const auto comps = components(g, ok);
  size_t acc = 0;
  for (char c : ok) acc += c;
  const double share = double(acc) / g.size();
  DataTable t;
  t.name = "voids";
  t.title = "Voids · probe " + fmt("%.2g Å", probe);
  t.columns = {"Void", "Volume (Å³)", "Points"};
  for (size_t c = 0; c < comps.size() && c < 500; ++c) t.rows.push_back({double(c + 1), comps[c].size() * g.voxel, double(comps[c].size())});
  st.tables.push_back(std::move(t));
  DataTable sw;   // how accessibility falls as the probe grows
  sw.name = "probe_sweep";
  sw.title = "Accessible volume against probe radius";
  sw.columns = {"Probe (Å)", "Accessible (%)", "Voids", "Largest share (%)"};
  for (double r : {0.0, 0.5, 1.0, 1.4, 2.0, 2.5}) {
    const auto o2 = accessible(s, g, r);
    const auto c2 = components(g, o2);
    size_t a2 = 0;
    for (char c : o2) a2 += c;
    sw.rows.push_back({r, 100.0 * a2 / g.size(), double(c2.size()), a2 ? 100.0 * (c2.empty() ? 0 : c2.front().size()) / a2 : 0});
  }
  st.tables.push_back(std::move(sw));
  st.set_attribute("Voids.accessible_fraction", share);
  st.set_attribute("Voids.count", double(comps.size()));
  st.set_attribute("Voids.largest_volume", comps.empty() ? 0 : comps.front().size() * g.voxel);
  if (flag(p, "show", true)) {   // void points coloured by void, at most about 6 000
    size_t total = 0;
    for (const auto& c : comps) total += c.size();
    const size_t every = std::max<size_t>(1, total / 6000);
    for (size_t c = 0; c < comps.size(); ++c)
      for (size_t q = 0; q < comps[c].size(); q += every) {
        const size_t idx = comps[c][q];
        const int i = int(idx / (size_t(g.n[1]) * g.n[2])), j = int((idx / g.n[2]) % g.n[1]), k = int(idx % g.n[2]);
        const Vec3 r = g.point(i, j, k);
        st.segments.push_back({r, r, kCat[c % 10], 0.14, false});
      }
  }
  out.summary = fmt("%.1f %% accessible", 100 * share) + " · " + std::to_string(comps.size()) + " voids · probe " + fmt("%.2g Å", probe) + " · grid " +
                std::to_string(g.n[0]) + "×" + std::to_string(g.n[1]) + "×" + std::to_string(g.n[2]);
}

void step_voronoi(PipelineState& st, const Json& p, StepStatus& out) {
  const System& s = st.system;
  if (!s.cell.valid()) throw std::invalid_argument("Voronoi volumes need a periodic cell");
  const bool radical = p.text("method", "grid") == "radical";
  const CellGrid g(s.cell, std::clamp(p.num("grid", 0.5), 0.15, 2.0));
  const size_t n = s.atoms.size();
  std::vector<double> best(g.size(), 1e300);
  std::vector<int> owner(g.size(), -1);
  // every atom claims the points within 6 Å; the few points farther from all atoms are settled one by one
  for (size_t a = 0; a < n; ++a) {
    const double w = radical ? vdw_radius(s.atoms[a]) : 0;
    g.around(s.atoms[a].pos, 6.0, [&](size_t k, double d2) {
      const double score = d2 - w * w;
      if (score < best[k]) { best[k] = score; owner[k] = int(a); }
    });
  }
  for (int i = 0; i < g.n[0]; ++i)
    for (int j = 0; j < g.n[1]; ++j)
      for (int k = 0; k < g.n[2]; ++k) {
        const size_t q = g.index(i, j, k);
        if (owner[q] >= 0) continue;
        const Vec3 pt = g.point(i, j, k);
        for (size_t a = 0; a < n; ++a) {
          const Vec3 d = s.cell.minimum_image(s.atoms[a].pos - pt);
          const double w = radical ? vdw_radius(s.atoms[a]) : 0;
          const double score = dot(d, d) - w * w;
          if (score < best[q]) { best[q] = score; owner[q] = int(a); }
        }
      }
  auto& vol = st.props["AtomicVolume"];
  vol.assign(n, 0);
  for (size_t k = 0; k < g.size(); ++k) if (owner[k] >= 0) vol[size_t(owner[k])] += g.voxel;
  std::map<std::string, std::vector<double>> by;
  for (size_t a = 0; a < n; ++a) by[s.atoms[a].name.empty() ? element(s.atoms[a].element).symbol : s.atoms[a].name].push_back(vol[a]);
  DataTable t;
  t.name = "voronoi";
  t.title = "Voronoi volume by type";
  t.columns = {"Type #", "Atoms", "mean (Å³)", "median (Å³)", "max (Å³)"};
  int k = 0;
  std::string names;
  for (auto& [name, v] : by) {
    std::sort(v.begin(), v.end());
    double m = 0;
    for (double x : v) m += x;
    t.rows.push_back({double(++k), double(v.size()), m / v.size(), v[v.size() / 2], v.back()});
    names += (names.empty() ? "" : ", ") + std::to_string(k) + " " + name;
  }
  st.tables.push_back(std::move(t));
  double sum = 0;
  for (double x : vol) sum += x;
  st.set_attribute("Voronoi.sum", sum);
  st.set_attribute("Voronoi.mean", n ? sum / n : 0);
  out.summary = std::string(radical ? "radical" : "nearest atom") + " · sum " + fmt("%.0f Å³", sum) + " of " + fmt("%.0f Å³", s.cell.volume()) + " · types " + names;
}

void step_density_field(PipelineState& st, const Json& p, StepStatus& out) {
  const System& s = st.system;
  if (!s.cell.valid()) throw std::invalid_argument("a density field needs a periodic cell");
  const double sigma = std::clamp(p.num("sigma", 1.5), 0.3, 6.0);
  const CellGrid g(s.cell, std::clamp(p.num("grid", 0.8), 0.2, 3.0));
  std::vector<double> rho(g.size(), 0);
  const double norm3 = 1.0 / std::pow(2 * M_PI * sigma * sigma, 1.5);
  // each atom's Gaussian is normalised on the grid, so the field holds exactly the cell's mass
  std::vector<std::pair<size_t, double>> kernel;
  for (const auto& a : s.atoms) {
    const double m = s.mass_of(a) * 1.66053906660;   // g/mol per Å³ → g/cm³
    kernel.clear();
    double sum = 0;
    g.around(a.pos, 3 * sigma, [&](size_t k, double d2) { const double w = std::exp(-d2 / (2 * sigma * sigma)); kernel.push_back({k, w}); sum += w; });
    if (sum <= 0) continue;
    for (const auto& [k, w] : kernel) rho[k] += m * w / (sum * g.voxel);
  }
  (void)norm3;
  double mean = 0, empty = 0;
  for (double x : rho) { mean += x; empty += x < 0.05; }
  mean /= g.size();
  st.set_attribute("DensityField.mean", mean);
  st.set_attribute("DensityField.empty_fraction", empty / g.size());
  const int axis = std::clamp(int(p.num("axis", 2)), 0, 2);
  DataTable prof;
  prof.name = "density_profile";
  const char* ax = axis == 0 ? "x" : axis == 1 ? "y" : "z";
  prof.title = std::string("Density profile along ") + ax;
  prof.columns = {std::string(ax) + " (fraction of the cell)", "Density (g/cm³)"};
  const int na = g.n[axis];
  std::vector<double> layer(size_t(na), 0);
  for (int i = 0; i < g.n[0]; ++i)
    for (int j = 0; j < g.n[1]; ++j)
      for (int k = 0; k < g.n[2]; ++k) layer[size_t(axis == 0 ? i : axis == 1 ? j : k)] += rho[g.index(i, j, k)];
  const double per = double(g.size()) / na;
  for (int q = 0; q < na; ++q) prof.rows.push_back({(q + 0.5) / na, layer[size_t(q)] / per});
  st.tables.push_back(std::move(prof));
  // the slice at the chosen position, as points coloured by density
  const double pos = std::clamp(p.num("position", 0.5), 0.0, 1.0);
  const int q = std::clamp(int(pos * na), 0, na - 1);
  double hi = 0;
  for (double x : rho) hi = std::max(hi, x);
  const double r = 0.45 * std::cbrt(g.voxel);
  for (int i = 0; i < g.n[0]; ++i)
    for (int j = 0; j < g.n[1]; ++j)
      for (int k = 0; k < g.n[2]; ++k) {
        if ((axis == 0 ? i : axis == 1 ? j : k) != q) continue;
        const double x = rho[g.index(i, j, k)];
        const Vec3 pt = g.point(i, j, k);
        st.segments.push_back({pt, pt, ramp(kViridis, 9, hi > 0 ? x / hi : 0), r, false});
      }
  st.legend = PipelineLegend{};
  st.legend.property = "Density (g/cm³)";
  st.legend.continuous = true;
  st.legend.lo = 0;
  st.legend.hi = hi;
  st.has_legend = true;
  out.summary = "mean " + fmt("%.3f g/cm³", mean) + " (cell " + fmt("%.3f", s.density()) + ") · " + fmt("%.0f %%", 100 * empty / g.size()) + " below 0.05 · σ " +
                fmt("%.2g Å", sigma);
}

void step_msd(PipelineState& st, const Json& p, StepStatus& out) {
  if (!st.traj || st.traj->frames() < 3) { out.level = "warning"; out.summary = "fewer than three frames: no MSD"; return; }
  const Trajectory& tr = *st.traj;
  const int nf = int(tr.frames());
  const int max_lag = std::clamp(int(p.num("max_lag", nf / 2)), 1, nf - 1);
  const bool heavy = flag(p, "heavy_only", true);
  // the lags the diffusion fit uses (frames); default: the last three quarters
  const int fit_lo = std::clamp(int(p.num("fit_from", std::max(1, max_lag / 4))), 1, max_lag), fit_hi = std::clamp(int(p.num("fit_to", max_lag)), fit_lo, max_lag);
  const int every = std::max(1, int(p.num("every", 1)));
  std::vector<uint32_t> atoms;
  for (uint32_t i = 0; i < tr.topology.atoms.size(); ++i)
    if ((!heavy || tr.topology.atoms[i].element != 1) && i % every == 0) atoms.push_back(i);
  int nm = 0;
  tr.topology.molecules(&nm);
  // continuous paths: frames made whole, then followed from frame to frame by minimum image
  std::vector<std::vector<Vec3>> ra(static_cast<size_t>(nf)), rc(static_cast<size_t>(nf));
  for (int f = 0; f < nf; ++f) {
    System s = tr.frame(size_t(f));
    if (!s.unwrapped && s.cell.valid()) make_molecules_whole(s);
    ra[size_t(f)].resize(atoms.size());
    for (size_t k = 0; k < atoms.size(); ++k) ra[size_t(f)][k] = s.atoms[atoms[k]].pos;
    const auto shapes = molecule_shapes(s);
    rc[size_t(f)].resize(shapes.size());
    for (size_t m = 0; m < shapes.size(); ++m) rc[size_t(f)][m] = shapes[m].com;
    if (f > 0 && s.cell.valid() && !s.unwrapped) {
      for (size_t k = 0; k < atoms.size(); ++k) ra[size_t(f)][k] = ra[size_t(f - 1)][k] + s.cell.minimum_image(ra[size_t(f)][k] - ra[size_t(f - 1)][k]);
      for (size_t m = 0; m < shapes.size() && m < rc[size_t(f - 1)].size(); ++m)
        rc[size_t(f)][m] = rc[size_t(f - 1)][m] + s.cell.minimum_image(rc[size_t(f)][m] - rc[size_t(f - 1)][m]);
    }
  }
  auto msd = [&](const std::vector<std::vector<Vec3>>& r, int lag) {
    double sum = 0;
    size_t cnt = 0;
    const int stride = std::max(1, (nf - lag) / 100);   // up to about 100 time origins
    for (int t0 = 0; t0 + lag < nf; t0 += stride)
      for (size_t k = 0; k < r[size_t(t0)].size() && k < r[size_t(t0 + lag)].size(); ++k) {
        const Vec3 d = r[size_t(t0 + lag)][k] - r[size_t(t0)][k];
        sum += dot(d, d);
        ++cnt;
      }
    return cnt ? sum / cnt : 0.0;
  };
  const double dts = tr.timesteps.size() >= 2 ? double(tr.timesteps.back() - tr.timesteps.front()) / (nf - 1) : 1.0;
  DataTable t;
  t.name = "msd";
  t.title = "Mean-square displacement";
  t.columns = {"Lag (timesteps)", "MSD atoms (Å²)", "MSD centres (Å²)", "Lag (frames)"};
  std::vector<double> x, y;
  for (int lag = 1; lag <= max_lag; ++lag) {
    const double a = msd(ra, lag), c = msd(rc, lag);
    t.rows.push_back({lag * dts, a, c, double(lag)});
    if (lag >= fit_lo && lag <= fit_hi) { x.push_back(lag * dts); y.push_back(c); }
  }
  st.tables.push_back(std::move(t));
  // D from the centres of mass: MSD = 6 D τ over the fit lags
  double D = 0;
  if (x.size() >= 2) {
    double mx = 0, my = 0;
    for (size_t k = 0; k < x.size(); ++k) { mx += x[k]; my += y[k]; }
    mx /= x.size(); my /= x.size();
    double sxy = 0, sxx = 0;
    for (size_t k = 0; k < x.size(); ++k) { sxy += (x[k] - mx) * (y[k] - my); sxx += (x[k] - mx) * (x[k] - mx); }
    D = sxx > 0 ? sxy / sxx / 6 : 0;   // Å² per timestep unit
  }
  const double fs = p.num("timestep_fs", 0);
  st.set_attribute("MSD.D_centres", D);
  st.set_attribute("MSD.fit_from", double(fit_lo)), st.set_attribute("MSD.fit_to", double(fit_hi));
  if (fs > 0) st.set_attribute("MSD.D_centres_cm2s", D / fs * 0.1);   // Å²/fs = 0.1 cm²/s
  out.summary = std::to_string(atoms.size()) + " atoms, " + std::to_string(nm) + " centres · lags to " + std::to_string(max_lag) + " frames · D " +
                (fs > 0 ? fmt("%.3g cm²/s", D / fs * 0.1) : fmt("%.3g Å²/timestep", D));
}

void step_scatter(PipelineState& st, const Json& p, StepStatus& out) {
  const std::string xn = p.text("x", "DistanceToCOM"), yn = p.text("y", "Charge");
  std::vector<double> x, y;
  if (!property_values(st, xn, x)) throw std::invalid_argument("no property " + xn);
  if (!property_values(st, yn, y)) throw std::invalid_argument("no property " + yn);
  const bool only_sel = flag(p, "only_selected", false);
  size_t n = 0;
  for (size_t i = 0; i < x.size(); ++i) n += (!only_sel || st.selected[i]) && std::isfinite(x[i]) && std::isfinite(y[i]);
  const size_t every = std::max<size_t>(1, n / 20000);
  DataTable t;
  t.name = "scatter";
  t.title = "Scatter · " + yn + " against " + xn;
  t.columns = {xn, yn};
  t.points = true;
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  size_t used = 0, seen = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    if ((only_sel && !st.selected[i]) || !std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
    sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; syy += y[i] * y[i]; sxy += x[i] * y[i]; ++used;
    if (seen++ % every == 0) t.rows.push_back({x[i], y[i]});
  }
  st.tables.push_back(std::move(t));
  double r = 0;
  if (used > 1) {
    const double cxy = sxy - sx * sy / used, cxx = sxx - sx * sx / used, cyy = syy - sy * sy / used;
    r = cxx > 0 && cyy > 0 ? cxy / std::sqrt(cxx * cyy) : 0;
  }
  st.set_attribute("Scatter.pearson_r", r);
  out.summary = std::to_string(used) + " points · Pearson r " + fmt("%.3f", r);
}

// ---------------------------------------------------------------- Python steps

std::string python_package_dir(const Json& p) {
  if (p.has("path") && p["path"].is_string() && !p["path"].str().empty()) return p["path"].str();
  if (const char* e = std::getenv("CAPS_PYTHON_PATH"); e && *e) return e;
  for (const char* c : {"data/python", "../data/python", "../../data/python", "../share/caps/data/python"})
    if (std::filesystem::exists(std::filesystem::path(c) / "caps" / "runner.py")) return std::filesystem::absolute(c).string();
  return "";
}

void step_python(PipelineState& st, const Json& p, StepStatus& out) {
  // a file, or the code typed in the step (written to a script named by its hash, so an edit is a new file)
  std::string script = p.text("file", "");
  const std::string code = p.text("code", "");
  if (script.empty() && !code.empty()) {
    const auto path = std::filesystem::temp_directory_path() / ("caps_step_" + sha256_hex(code).substr(0, 16) + ".py");
    if (!std::filesystem::exists(path)) {
      std::ofstream f(path);
      f << code << (code.back() == '\n' ? "" : "\n");
    }
    script = path.string();
  }
  if (script.empty()) throw std::invalid_argument("choose a Python file with an @step function, or type the step");
  if (!std::filesystem::exists(script)) throw std::invalid_argument("no file " + script);
  const std::string pkg = python_package_dir(p);
  if (pkg.empty()) throw std::runtime_error("the caps Python package was not found (set CAPS_PYTHON_PATH to data/python)");
  const char* py_env = std::getenv("CAPS_PYTHON");
#ifdef _WIN32
  const std::string python = py_env && *py_env ? py_env : "python";
#else
  const std::string python = py_env && *py_env ? py_env : "python3";
#endif
  // the frame, molecules whole, as the script sees it
  System whole = st.system;
  if (!whole.unwrapped && whole.cell.valid()) make_molecules_whole(whole);
  const size_t n = whole.atoms.size();
  std::vector<double> backbone(n, 0);
  for (const auto& bb : backbones(whole)) for (uint32_t i : bb) backbone[i] = 1;
  std::vector<double> mol;
  property_values(st, "Molecule", mol);
  Json parts = Json::object();
  auto column = [&](const std::string& name, auto&& f) {
    Json a = Json::array();
    for (size_t i = 0; i < n; ++i) a.push_back(f(i));
    parts[name] = std::move(a);
  };
  column("Particle Identifier", [&](size_t i) { return Json(double(whole.atoms[i].id)); });
  column("Molecule Identifier", [&](size_t i) { return Json(mol[i]); });
  column("Particle Type", [&](size_t i) { return Json(double(whole.atoms[i].type)); });
  column("Element", [&](size_t i) { return Json(std::string(element(whole.atoms[i].element).symbol)); });
  column("Charge", [&](size_t i) { return Json(whole.atoms[i].charge); });
  column("Selection", [&](size_t i) { return Json(double(st.selected[i])); });
  column("Backbone", [&](size_t i) { return Json(backbone[i]); });
  column("Position", [&](size_t i) {
    Json r = Json::array();
    for (int k = 0; k < 3; ++k) r.push_back(whole.atoms[i].pos[k]);
    return r;
  });
  for (const auto& [name, v] : st.props) column(name, [&](size_t i) { return Json(v[i]); });
  Json in = Json::object();
  in["count"] = double(n);
  in["frame"] = double(st.frame);
  in["particles"] = std::move(parts);
  Json bonds = Json::array();
  for (const auto& b : whole.bonds) { Json pr = Json::array(); pr.push_back(double(b.i)); pr.push_back(double(b.j)); bonds.push_back(std::move(pr)); }
  in["bonds"] = std::move(bonds);
  if (whole.cell.valid()) {
    Json c = Json::array();
    for (const Vec3& e : {whole.cell.a, whole.cell.b, whole.cell.c}) { Json r = Json::array(); for (int k = 0; k < 3; ++k) r.push_back(e[k]); c.push_back(std::move(r)); }
    in["cell"] = std::move(c);
  }
  Json attrs = Json::object();
  for (const auto& [k, v] : st.attributes) attrs[k] = v;
  in["attributes"] = std::move(attrs);
  const auto dir = std::filesystem::temp_directory_path();
  const auto tag = std::to_string(reinterpret_cast<uintptr_t>(&st)) + "_" + std::to_string(st.frame);
  const auto fin = dir / ("caps_step_in_" + tag + ".json"), fout = dir / ("caps_step_out_" + tag + ".json");
  { std::ofstream f(fin); f << in.dump(0); }
  std::filesystem::remove(fout);
  // run it with the package on PYTHONPATH (the child inherits the environment)
  std::string old = std::getenv("PYTHONPATH") ? std::getenv("PYTHONPATH") : "";
#ifdef _WIN32
  _putenv_s("PYTHONPATH", (pkg + (old.empty() ? "" : ";" + old)).c_str());
  const std::string cmd = "\"\"" + python + "\" -m caps.runner \"" + script + "\" \"" + fin.string() + "\" \"" + fout.string() + "\" 2>&1\"";
  FILE* pipe = _popen(cmd.c_str(), "r");
#else
  setenv("PYTHONPATH", (pkg + (old.empty() ? "" : ":" + old)).c_str(), 1);
  const std::string cmd = "'" + python + "' -m caps.runner '" + script + "' '" + fin.string() + "' '" + fout.string() + "' 2>&1";
  FILE* pipe = popen(cmd.c_str(), "r");
#endif
  std::string log;
  if (pipe) {
    char buf[512];
    while (std::fgets(buf, sizeof buf, pipe)) log += buf;
#ifdef _WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif
  }
#ifdef _WIN32
  _putenv_s("PYTHONPATH", old.c_str());
#else
  if (old.empty()) unsetenv("PYTHONPATH"); else setenv("PYTHONPATH", old.c_str(), 1);
#endif
  std::filesystem::remove(fin);
  out.output = log.size() > 8000 ? "…" + log.substr(log.size() - 8000) : log;   // the console, kept on errors too
  if (!std::filesystem::exists(fout)) {   // Python itself failed: say which error (the traceback's last line)
    std::string last = log;
    while (!last.empty() && (last.back() == '\n' || last.back() == '\r')) last.pop_back();
    if (const auto nl = last.rfind('\n'); nl != std::string::npos) last = last.substr(nl + 1);
    throw std::runtime_error("Python did not run" + (log.empty() ? std::string(" (is ") + python + " installed?)" : ": " + last.substr(0, 300)));
  }
  Json res;
  { std::ifstream f(fout); res = Json::parse(std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>())); }
  std::filesystem::remove(fout);
  if (!res.has("ok") || !res["ok"].boolean()) {
    std::string where;
    if (res.has("trace")) {   // the script's own line
      const std::string& tr = res["trace"].str();
      out.output += tr;
      const auto at = tr.rfind(std::filesystem::path(script).filename().string());
      if (at != std::string::npos) where = " · " + tr.substr(at, tr.find('\n', at) - at);
    }
    throw std::runtime_error(res.text("error", "the step failed") + where);
  }
  size_t na = 0, np = 0, nt = 0;
  if (res.has("attributes"))
    for (const auto& [k, v] : res["attributes"].members())
      if (v.is_number()) {   // new or changed ones (values that went through JSON unchanged are not counted)
        const double before = st.attribute(k, std::nan(""));
        if (!(std::fabs(before - v.number()) <= 1e-9 * std::max(1.0, std::fabs(before)))) { st.set_attribute(k, v.number()); ++na; }
      }
  if (res.has("properties"))
    for (const auto& [k, v] : res["properties"].members()) {
      if (!v.is_array() || v.size() != n) continue;
      auto& dst = st.props[k];
      dst.resize(n);
      for (size_t i = 0; i < n; ++i) dst[i] = v[i].number();
      ++np;
    }
  if (res.has("tables"))
    for (const auto& t : res["tables"].items()) {
      DataTable d;
      d.name = t.text("name", "python");
      d.title = res.text("name", "Python") + " · " + d.name;
      for (const auto& c : t["columns"].items()) d.columns.push_back(c.str());
      for (const auto& r : t["rows"].items()) {
        std::vector<double> row;
        for (const auto& x : r.items()) row.push_back(x.number());
        d.rows.push_back(std::move(row));
      }
      st.tables.push_back(std::move(d));
      ++nt;
    }
  if (res.has("selection") && res["selection"].size() == n)
    for (size_t i = 0; i < n; ++i) st.selected[i] = res["selection"][i].number() != 0;
  out.title = res.text("name", "Python step");
  out.summary = std::to_string(na) + " attribute" + (na == 1 ? "" : "s") + ", " + std::to_string(np) + " propert" + (np == 1 ? "y" : "ies") + ", " +
                std::to_string(nt) + " table" + (nt == 1 ? "" : "s") + (log.empty() ? "" : [&] { const auto nl = std::count(log.begin(), log.end(), '\n'); return " · printed " + std::to_string(nl) + (nl == 1 ? " line" : " lines"); }());
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
    {"freeze_property", "Freeze property", "a property's values at a reference frame, on every frame", step_freeze},
    {"orientation", "Chain orientation", "backbone chords: P₂ per atom, S, director, local crystallinity", step_orientation},
    {"affine_transform", "Affine transformation", "strain, shear or rotate particles and cell", step_affine},
    {"cna", "Common neighbour analysis", "adaptive CNA: FCC, HCP, BCC, icosahedral, other", step_cna},
    {"centrosymmetry", "Centrosymmetry", "Kelchner's parameter from the N nearest neighbours", step_centrosymmetry},
    {"molecule_shape", "Molecule shape", "Rg, κ², asphericity per molecule", step_molecule_shape},
    {"topology", "Topology distributions", "bond lengths, angles, dihedrals", step_topology},
    {"displacements", "Displacements", "vs a reference frame, MSD", step_displacements},
    {"smooth", "Smooth trajectory", "positions averaged over a window of frames", step_smooth},
    {"vectors", "Vectors", "end-to-end, dipoles, displacements, velocities as arrows", step_vectors},
    {"trajectory_lines", "Trajectory lines", "paths of molecule centres or particles", step_trajectory_lines},
    {"primitive_paths", "Primitive paths", "chains pulled tight without crossing: entanglements, N_e", step_primitive_paths},
    {"voids", "Voids & pores", "accessible volume for a probe, voids by size", step_voids},
    {"voronoi", "Voronoi volumes", "volume per atom on a grid", step_voronoi},
    {"density_field", "Density field", "Gaussian-smoothed mass density, profile, slice", step_density_field},
    {"msd", "Mean-square displacement", "MSD(τ) of atoms and chain centres, diffusion", step_msd},
    {"scatter", "Scatter plot", "one property against another", step_scatter},
    {"python", "Python step", "a script with an @step function, run in Python", step_python},
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
  st.pipeline = &p;
  for (size_t k = p.steps.size(); k-- > 0;) {   // bottom to top
    const auto& step = p.steps[k];
    StepStatus& out = st.steps[k];
    st.step_index = k;
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
                                    "Selection", "DistanceToCOM", "Monomer", "MoleculeCOM.X", "MoleculeCOM.Y", "MoleculeCOM.Z"};
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
  if (name == "Monomer" || name == "Residue") return each([&](size_t i) { return double(s.atoms[i].resid); });   // Grow numbers units from 1
  if (name == "MoleculeCOM.X" || name == "MoleculeCOM.Y" || name == "MoleculeCOM.Z") {
    // the molecule's centre of mass placed by the particle's own image: Position − MoleculeCOM is its minimum-image offset
    const int c = name.back() - 'X';
    System whole = s;
    if (!whole.unwrapped) make_molecules_whole(whole);
    const auto shapes = molecule_shapes(whole);
    const auto mol = whole.molecules();
    return each([&](size_t i) { return s.atoms[i].pos[size_t(c)] + shapes[size_t(mol[i])].com[size_t(c)] - whole.atoms[i].pos[size_t(c)]; });
  }
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
    if (!s.output.empty()) o["output"] = s.output;
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
    o["points"] = t.points;
    if (!t.labels.empty()) {
      Json labels = Json::array();
      for (const auto& l : t.labels) labels.push_back(l);
      o["labels"] = std::move(labels);
      o["label_column"] = t.label_column;
    }
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
