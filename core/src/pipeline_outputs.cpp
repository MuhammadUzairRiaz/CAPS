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

void write_text(const std::filesystem::path& p, const std::string& s) {
  std::ofstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + p.string());
  f << s;
}

// a round major step (1, 2, 2.5 or 5 times a power of ten) for about `target` intervals, and its minor step
double nice_step(double span, int target) {
  const double raw = span / std::max(1, target), mag = std::pow(10.0, std::floor(std::log10(std::max(1e-300, raw)))), f = raw / mag;
  return (f < 1.5 ? 1 : f < 2.25 ? 2 : f < 3.5 ? 2.5 : f < 7.5 ? 5 : 10) * mag;
}
double minor_of(double step) {
  const double f = step / std::pow(10.0, std::floor(std::log10(step)));
  return step / (std::fabs(f - 2) < 1e-9 ? 4 : 5);
}
std::vector<double> ticks(double lo, double hi, double step) {
  std::vector<double> t;
  for (double k = std::ceil(lo / step - 1e-9); k * step <= hi + 1e-9 * step; k += 1) t.push_back(std::fabs(k * step) < 1e-9 * step ? 0.0 : k * step);
  return t;
}
// a tick value with the decimals its step needs, a true minus sign
std::string tick_text(double v, double step) {
  char b[64];
  const double a = std::fabs(v);
  if (a != 0 && (a >= 1e5 || a < 1e-3)) std::snprintf(b, sizeof b, "%.1e", v);
  else {
    int dec = 0;
    while (dec < 6 && std::fabs(step * std::pow(10.0, dec) - std::round(step * std::pow(10.0, dec))) > 1e-6 * std::pow(10.0, dec)) ++dec;
    std::snprintf(b, sizeof b, "%.*f", dec, v);
  }
  std::string s;
  for (const char* p = b; *p; ++p) s += *p == '-' ? std::string("−") : std::string(1, *p);
  return s;
}
// "Time (ps)" → "Time / ps" (an axis title as quantity / unit)
std::string axis_title(const std::string& c) {
  const auto o = c.rfind(" ("), e = c.rfind(')');
  if (o == std::string::npos || e != c.size() - 1 || e < o) return c;
  return c.substr(0, o) + " / " + c.substr(o + 2, e - o - 2);
}
// width of a string in em of the font (Helvetica-like proportions; digits and most letters about 0.56 em)
double text_em(const std::string& s) {
  double w = 0;
  for (size_t i = 0; i < s.size(); ++i) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if ((c & 0xC0) == 0x80) continue;   // a UTF-8 continuation byte
    w += c == ' ' || c == '.' || c == ',' || c == 'i' || c == 'l' || c == 'j' || c == 't' || c == 'f' || c == 'r' || c == '(' || c == ')' || c == '/' ? 0.3
         : c == 'm' || c == 'w' || c == 'M' || c == 'W' ? 0.85 : c >= 'A' && c <= 'Z' ? 0.68 : 0.56;
  }
  return w;
}

}  // namespace

std::string table_svg(const DataTable& t, int width, int height) {
  // a journal figure: sizes in points (1/72 in), 8 pt text, a closed frame with inward major and minor ticks on all
  // four sides at round values, axis titles as quantity / unit, Okabe and Ito's colours, open markers on series of up
  // to 150 points, no legend frame
  const double fs = 8.0, fsT = 7.5, tick = 3.5, tminor = 2.0, lw = 0.6;
  static const char* colours[] = {"#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"};
  static const char* dashes[] = {"", "4 2", "1.2 1.6", "5 2 1.2 2"};
  const size_t ncol = t.columns.empty() ? (t.rows.empty() ? 0 : t.rows[0].size()) : t.columns.size();
  double x0 = 1e300, x1 = -1e300, y0 = 1e300, y1 = -1e300;
  std::vector<size_t> count(ncol, 0);
  for (const auto& r : t.rows) {
    if (r.empty() || !std::isfinite(r[0])) continue;
    x0 = std::min(x0, r[0]), x1 = std::max(x1, r[0]);
    for (size_t c = 1; c < r.size() && c < ncol; ++c)
      if (std::isfinite(r[c])) y0 = std::min(y0, r[c]), y1 = std::max(y1, r[c]), ++count[c];
  }
  if (!(x1 >= x0)) x0 = 0, x1 = 1;
  if (!(y1 >= y0)) y0 = 0, y1 = 1;
  if (x1 == x0) x0 -= 0.5, x1 += 0.5;
  if (y1 == y0) { const double d = y0 == 0 ? 1 : std::fabs(y0) * 0.05; y0 -= d, y1 += d; }
  const double ypad = 0.04 * (y1 - y0);
  if (!(y0 >= 0 && y0 < 0.3 * (y1 - y0))) y0 -= ypad;   // a quantity that comes near zero starts at zero
  else y0 = 0;
  y1 += ypad;
  const double xs = nice_step(x1 - x0, 5), ys = nice_step(y1 - y0, 5);
  x0 = std::floor(x0 / minor_of(xs) + 1e-9) * minor_of(xs), x1 = std::ceil(x1 / minor_of(xs) - 1e-9) * minor_of(xs);
  y0 = std::floor(y0 / minor_of(ys) + 1e-9) * minor_of(ys), y1 = std::ceil(y1 / minor_of(ys) - 1e-9) * minor_of(ys);
  double ytw = 0;
  for (double y : ticks(y0, y1, ys)) ytw = std::max(ytw, text_em(tick_text(y, ys)) * fsT);
  const bool ytitle = ncol == 2 && t.columns.size() == 2;
  const double ml = 2 + (ytitle ? fs * 1.3 : 0) + ytw + 3.5, mr = 6, mt = 6, mb = 3.5 + fsT * 1.2 + (t.columns.empty() ? 0 : fs * 1.4) + 2;
  const double W = width - ml - mr, H = height - mt - mb;
  auto X = [&](double x) { return ml + (x - x0) / (x1 - x0) * W; };
  auto Y = [&](double y) { return mt + (1 - (y - y0) / (y1 - y0)) * H; };
  char b[512];
  std::string s;
  std::snprintf(b, sizeof b, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.3fin\" height=\"%.3fin\" viewBox=\"0 0 %d %d\" font-family=\"Arial, Helvetica, sans-serif\" font-size=\"%g\">\n",
                width / 72.0, height / 72.0, width, height, fs);
  s += b;
  s += "<title>" + esc(t.title.empty() ? t.name : t.title) + "</title>\n";
  s += "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
  // the series, clipped to the frame
  std::snprintf(b, sizeof b, "<clipPath id=\"frame\"><rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\"/></clipPath>\n<g clip-path=\"url(#frame)\">\n", ml, mt, W, H);
  s += b;
  for (size_t c = 1; c < ncol; ++c) {
    const char* col = colours[(c - 1) % 7];
    const bool marks = t.points || count[c] <= 150;
    if (!t.points) {
      std::string d;
      bool pen = false;
      for (const auto& r : t.rows) {
        if (c >= r.size() || !std::isfinite(r[0]) || !std::isfinite(r[c])) { pen = false; continue; }
        std::snprintf(b, sizeof b, "%s%.2f %.2f ", pen ? "L" : "M", X(r[0]), Y(r[c]));
        d += b;
        pen = true;
      }
      const char* da = dashes[(c - 1) % 4];
      s += "<path d=\"" + d + "\" fill=\"none\" stroke=\"" + col + "\" stroke-width=\"1\" stroke-linejoin=\"round\"" + (*da ? std::string(" stroke-dasharray=\"") + da + "\"" : "") + "/>\n";
    }
    if (marks)
      for (const auto& r : t.rows)
        if (c < r.size() && std::isfinite(r[0]) && std::isfinite(r[c])) {
          std::snprintf(b, sizeof b, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%s\" fill=\"%s\" stroke=\"%s\" stroke-width=\"0.6\"/>\n", X(r[0]), Y(r[c]),
                        t.points ? "1.2" : "1.5", t.points ? col : "white", col);
          s += b;
        }
  }
  s += "</g>\n";
  // the frame and the ticks, inward on all four sides
  std::string tk;
  auto line = [&](double ax, double ay, double bx, double by) { std::snprintf(b, sizeof b, "M%.2f %.2fL%.2f %.2f", ax, ay, bx, by); tk += b; };
  for (double x : ticks(x0, x1, minor_of(xs))) line(X(x), mt + H, X(x), mt + H - tminor), line(X(x), mt, X(x), mt + tminor);
  for (double y : ticks(y0, y1, minor_of(ys))) line(ml, Y(y), ml + tminor, Y(y)), line(ml + W, Y(y), ml + W - tminor, Y(y));
  for (double x : ticks(x0, x1, xs)) line(X(x), mt + H, X(x), mt + H - tick), line(X(x), mt, X(x), mt + tick);
  for (double y : ticks(y0, y1, ys)) line(ml, Y(y), ml + tick, Y(y)), line(ml + W, Y(y), ml + W - tick, Y(y));
  std::snprintf(b, sizeof b, "\" fill=\"none\" stroke=\"black\" stroke-width=\"%.2f\"/>\n", lw * 0.8);
  s += "<path d=\"" + tk + b;
  std::snprintf(b, sizeof b, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"none\" stroke=\"black\" stroke-width=\"%.2f\"/>\n", ml, mt, W, H, lw);
  s += b;
  for (double x : ticks(x0, x1, xs)) {
    std::snprintf(b, sizeof b, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"middle\" font-size=\"%g\">", X(x), mt + H + 3.5 + fsT * 0.8, fsT);
    s += b + tick_text(x, xs) + "</text>\n";
  }
  for (double y : ticks(y0, y1, ys)) {
    std::snprintf(b, sizeof b, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"end\" font-size=\"%g\">", ml - 3.5, Y(y) + fsT * 0.35, fsT);
    s += b + tick_text(y, ys) + "</text>\n";
  }
  if (!t.columns.empty()) {
    std::snprintf(b, sizeof b, "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"middle\">", ml + W / 2, double(height) - 3);
    s += b + esc(axis_title(t.columns[0])) + "</text>\n";
  }
  if (ytitle) {
    std::snprintf(b, sizeof b, "<text transform=\"translate(%.2f,%.2f) rotate(-90)\" text-anchor=\"middle\">", 2 + fs * 0.8, mt + H / 2);
    s += b + esc(axis_title(t.columns[1])) + "</text>\n";
  }
  // a legend without a frame, top right, for more than one series
  if (ncol > 2)
    for (size_t c = 1; c < ncol && c < t.columns.size(); ++c) {
      const char* col = colours[(c - 1) % 7];
      const double yy = mt + tick + 3 + fsT * 1.25 * double(c - 1) + fsT * 0.5;
      const double tx = ml + W - tick - 3 - text_em(axis_title(t.columns[c])) * fsT;
      const char* da = dashes[(c - 1) % 4];
      if (!t.points) {
        std::snprintf(b, sizeof b, "<path d=\"M%.2f %.2fh14\" stroke=\"%s\" stroke-width=\"1\"%s%s%s/>\n", tx - 18, yy, col, *da ? " stroke-dasharray=\"" : "", da, *da ? "\"" : "");
        s += b;
      }
      if (t.points || count[c] <= 150) {
        std::snprintf(b, sizeof b, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"1.5\" fill=\"%s\" stroke=\"%s\" stroke-width=\"0.6\"/>\n", tx - 11, yy, t.points ? col : "white", col);
        s += b;
      }
      std::snprintf(b, sizeof b, "<text x=\"%.2f\" y=\"%.2f\" font-size=\"%g\">", tx, yy + fsT * 0.35, fsT);
      s += b + esc(axis_title(t.columns[c])) + "</text>\n";
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
