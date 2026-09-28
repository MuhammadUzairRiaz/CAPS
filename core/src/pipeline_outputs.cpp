// CAPS pipeline outputs (see pipeline.hpp): the files a saved pipeline writes — tables as CSV and SVG plots, the global
// attributes, a render of the view and the last grid.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include "caps/appearance.hpp"
#include "caps/bundle.hpp"
#include "caps/pipeline.hpp"
#include "caps/render.hpp"

namespace caps {

namespace {

std::string esc(const std::string& s) {
  std::string o;
  for (char c : s) {
    if (c == '<') o += "&lt;";
    else if (c == '>') o += "&gt;";
    else if (c == '&') o += "&amp;";
    else if (c == '"') o += "&quot;";
    else o += c;
  }
  return o;
}

// about five round ticks between lo and hi
std::vector<double> ticks(double lo, double hi) {
  std::vector<double> t;
  if (!(hi > lo)) return {lo};
  const double raw = (hi - lo) / 5, mag = std::pow(10.0, std::floor(std::log10(raw))), r = raw / mag;
  const double step = (r < 1.5 ? 1 : r < 3 ? 2 : r < 7 ? 5 : 10) * mag;
  for (double x = std::ceil(lo / step) * step; x <= hi + 1e-9 * step; x += step) t.push_back(std::fabs(x) < 1e-12 * step ? 0.0 : x);
  return t;
}

void write_text(const std::filesystem::path& p, const std::string& s) {
  std::ofstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + p.string());
  f << s;
}

}  // namespace

std::string table_svg(const DataTable& t, int width, int height) {
  const double ml = 64, mr = 16, mt = 30, mb = 46;
  const double W = width - ml - mr, H = height - mt - mb;
  double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
  for (const auto& r : t.rows) {
    if (r.empty() || !std::isfinite(r[0])) continue;
    x0 = std::min(x0, r[0]), x1 = std::max(x1, r[0]);
    for (size_t c = 1; c < r.size(); ++c)
      if (std::isfinite(r[c])) y0 = std::min(y0, r[c]), y1 = std::max(y1, r[c]);
  }
  if (!(x1 >= x0)) x0 = 0, x1 = 1;
  if (!(y1 >= y0)) y0 = 0, y1 = 1;
  if (x1 == x0) x0 -= 0.5, x1 += 0.5;
  if (y1 == y0) y0 -= 0.5, y1 += 0.5;
  const double pad = 0.05 * (y1 - y0);
  y0 -= pad, y1 += pad;
  auto X = [&](double x) { return ml + (x - x0) / (x1 - x0) * W; };
  auto Y = [&](double y) { return mt + (1 - (y - y0) / (y1 - y0)) * H; };
  static const char* colours[] = {"#D98A1E", "#3B82C4", "#4E9A52", "#C24D6B", "#7A5FC0", "#2A9D8F"};
  char b[256];
  std::string s;
  std::snprintf(b, sizeof b, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"11\">\n",
                width, height, width, height);
  s += b;
  s += "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
  s += "<text x=\"" + std::to_string(int(ml)) + "\" y=\"18\" font-size=\"13\" fill=\"#222\">" + esc(t.title.empty() ? t.name : t.title) + "</text>\n";
  // axes and ticks
  std::snprintf(b, sizeof b, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#888\"/>\n", ml, mt, W, H);
  s += b;
  for (double x : ticks(x0, x1)) {
    std::snprintf(b, sizeof b, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#888\"/><text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" fill=\"#444\">%g</text>\n",
                  X(x), mt + H, X(x), mt + H + 4, X(x), mt + H + 16, x);
    s += b;
  }
  for (double y : ticks(y0, y1)) {
    std::snprintf(b, sizeof b, "<line x1=\"%.1f\" y1=\"%.1f\" x2=\"%.1f\" y2=\"%.1f\" stroke=\"#888\"/><text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" fill=\"#444\">%g</text>\n",
                  ml - 4, Y(y), ml, Y(y), ml - 6, Y(y) + 4, y);
    s += b;
  }
  if (!t.columns.empty()) {
    std::snprintf(b, sizeof b, "<text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\" fill=\"#222\">", ml + W / 2, height - 10);
    s += b + esc(t.columns[0]) + "</text>\n";
  }
  if (t.columns.size() == 2) {
    std::snprintf(b, sizeof b, "<text transform=\"translate(14,%.1f) rotate(-90)\" text-anchor=\"middle\" fill=\"#222\">", mt + H / 2);
    s += b + esc(t.columns[1]) + "</text>\n";
  }
  // the series
  const size_t ncol = t.columns.empty() ? (t.rows.empty() ? 0 : t.rows[0].size()) : t.columns.size();
  for (size_t c = 1; c < ncol; ++c) {
    const char* col = colours[(c - 1) % 6];
    if (t.points) {
      for (const auto& r : t.rows)
        if (c < r.size() && std::isfinite(r[0]) && std::isfinite(r[c])) {
          std::snprintf(b, sizeof b, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"2\" fill=\"%s\"/>\n", X(r[0]), Y(r[c]), col);
          s += b;
        }
    } else {
      std::string d;
      bool pen = false;
      for (const auto& r : t.rows) {
        if (c >= r.size() || !std::isfinite(r[0]) || !std::isfinite(r[c])) { pen = false; continue; }
        std::snprintf(b, sizeof b, "%s%.1f %.1f ", pen ? "L" : "M", X(r[0]), Y(r[c]));
        d += b;
        pen = true;
      }
      s += "<path d=\"" + d + "\" fill=\"none\" stroke=\"" + col + "\" stroke-width=\"1.6\"/>\n";
    }
    if (ncol > 2 && c < t.columns.size()) {   // a legend line per series
      std::snprintf(b, sizeof b, "<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"end\" fill=\"%s\">", ml + W - 6, mt + 14 + 13.0 * double(c - 1), col);
      s += b + esc(t.columns[c]) + "</text>\n";
    }
  }
  s += "</svg>\n";
  return s;
}

RenderOptions pipeline_render_options(const PipelineState& st, const RenderOptions& base) {
  RenderOptions r = base;
  r.colours = st.colour;
  for (size_t i = 0; i < r.colours.size() && i < st.selected.size(); ++i) if (st.selected[i]) r.colours[i] = 0xE5484D;
  r.segments = st.segments;
  r.highlight.clear();
  r.focus = -1;
  for (const auto& pm : st.meshes) if (pm.mesh) r.meshes.push_back({pm.mesh.get(), pm.rgb, pm.opacity});
  const size_t n = st.system.atoms.size();
  if (auto it = st.props.find("Radius"); it != st.props.end() && it->second.size() == n) r.radius.assign(it->second.begin(), it->second.end());
  if (auto it = st.props.find("Transparency"); it != st.props.end() && it->second.size() == n) r.transparency.assign(it->second.begin(), it->second.end());
  if (r.colour_by == ColourBy::Property) property_values(st, "DistanceToCOM", r.property);
  return r;
}

std::vector<std::string> write_pipeline_outputs(const PipelineState& st, const Pipeline& p, const std::string& dir, const Camera& cam) {
  std::vector<std::string> log;
  if (p.outputs.empty()) return log;
  const std::filesystem::path root = dir.empty() ? std::filesystem::path(".") : std::filesystem::path(dir);
  std::filesystem::create_directories(root);
  auto table = [&](const std::string& name) -> const DataTable* {
    for (const auto& t : st.tables) if (t.name == name) return &t;
    return nullptr;
  };
  for (const auto& o : p.outputs) {
    const std::filesystem::path path = root / o.path;
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    try {
      if (o.kind == "table" || o.kind == "plot") {
        const DataTable* t = table(o.what);
        if (!t) {
          std::string names;
          for (const auto& x : st.tables) names += (names.empty() ? "" : ", ") + x.name;
          log.push_back("not written: " + o.path + " — no table \"" + o.what + "\" (tables: " + (names.empty() ? "none" : names) + ")");
          continue;
        }
        write_text(path, o.kind == "table" ? table_csv(*t) : table_svg(*t));
      } else if (o.kind == "attributes") {
        std::string csv = "attribute,value\n";
        char b[64];
        for (const auto& [k, v] : st.attributes) {
          std::snprintf(b, sizeof b, "%.10g", v);
          csv += (k.find(',') != std::string::npos ? "\"" + k + "\"" : k) + "," + b + "\n";
        }
        write_text(path, csv);
      } else if (o.kind == "render") {
        RenderOptions base;
        base.width = std::clamp(o.width, 64, 8192);
        base.height = std::clamp(o.height, 64, 8192);
        const RenderOptions r = pipeline_render_options(st, base);
        const std::string ext = path.extension().string();
        if (ext == ".svg") write_text(path, render_svg(st.system, cam, r));
        else {
          Renderer renderer;
          write_png(renderer.render(st.system, cam, r), path.string());
        }
      } else if (o.kind == "grid") {
        if (!st.grid) { log.push_back("not written: " + o.path + " — no step made a grid (Density field)"); continue; }
        write_grid(*st.grid, st.system, path.string());
      } else {
        log.push_back("not written: " + o.path + " — unknown output kind \"" + o.kind + "\" (table, plot, attributes, render, grid)");
        continue;
      }
      log.push_back("wrote " + path.string());
    } catch (const std::exception& e) {
      log.push_back("not written: " + o.path + " — " + e.what());
    }
  }
  return log;
}

}  // namespace caps
