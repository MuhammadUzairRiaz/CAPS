// CAPS readers for measured data and grids (see caps/measured.hpp).
#include "caps/measured.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"

namespace caps {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kBohr = 0.529177210903;   // Å

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
std::string upper(std::string s) {
  for (char& c : s) c = char(std::toupper(static_cast<unsigned char>(c)));
  return s;
}
std::vector<std::string> lines_of(const std::string& text) {
  std::vector<std::string> out;
  std::istringstream in(text);
  for (std::string l; std::getline(in, l);) {
    if (!l.empty() && l.back() == '\r') l.pop_back();
    out.push_back(l);
  }
  return out;
}
std::string extension(const std::string& name) {
  const auto dot = name.find_last_of('.');
  const auto slash = name.find_last_of("/\\");
  if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return "";
  return upper(name.substr(dot + 1));
}
// every number in s (blanks, commas and semicolons separate; a sign starts a new number, as PAC allows)
std::vector<double> numbers(const std::string& s) {
  std::vector<double> v;
  const char* p = s.c_str();
  while (*p) {
    while (*p == ' ' || *p == '\t' || *p == ',' || *p == ';') ++p;
    if (!*p) break;
    if (*p == '?') { v.push_back(kNaN); ++p; continue; }
    char* e = nullptr;
    const double x = std::strtod(p, &e);
    if (e == p) { ++p; continue; }
    v.push_back(x);
    p = e;
  }
  return v;
}

// ---------------------------------------------------------------- JCAMP-DX

// The ASDF pseudo-digits (McDonald & Wilks 1988): SQZ @A–I (+0…+9) a–i (−1…−9); DIF % J–R (+0…+9) j–r (−1…−9);
// DUP S–Z s (1…9).
bool sqz(char c, int& d) {
  if (c == '@') { d = 0; return true; }
  if (c >= 'A' && c <= 'I') { d = c - 'A' + 1; return true; }
  if (c >= 'a' && c <= 'i') { d = -(c - 'a' + 1); return true; }
  return false;
}
bool dif(char c, int& d) {
  if (c == '%') { d = 0; return true; }
  if (c >= 'J' && c <= 'R') { d = c - 'J' + 1; return true; }
  if (c >= 'j' && c <= 'r') { d = -(c - 'j' + 1); return true; }
  return false;
}
bool dup(char c, int& d) {
  if (c >= 'S' && c <= 'Z') { d = c - 'S' + 1; return true; }
  if (c == 's') { d = 9; return true; }
  return false;
}
// A pseudo-digit's value with the plain digits after it: the sign from the pseudo-digit, then |d| and the digits.
double with_digits(int d, const std::string& s, size_t& i) {
  const bool neg = d < 0;
  std::string m = std::to_string(std::abs(d));
  while (i < s.size() && (std::isdigit(static_cast<unsigned char>(s[i])) || s[i] == '.')) m += s[i++];
  const double v = std::strtod(m.c_str(), nullptr);
  return neg ? -v : v;
}
bool asdf_line(const std::string& s) {
  for (char c : s) {
    int d;
    if (c == 'E' || c == 'e') continue;   // also an exponent: decides nothing alone
    if (sqz(c, d) || dif(c, d) || dup(c, d)) return true;
  }
  return false;
}

// One XYDATA line: the abscissa and the ordinates; `ends_dif` whether the last ordinate is in DIF form.
void decode_xydata_line(const std::string& raw, double& x, std::vector<double>& ys, bool& ends_dif) {
  std::string s = raw;
  if (const auto c = s.find("$$"); c != std::string::npos) s = s.substr(0, c);
  ys.clear();
  ends_dif = false;
  if (!asdf_line(s)) {   // AFFN / PAC
    const auto v = numbers(s);
    if (v.empty()) throw std::invalid_argument("JCAMP-DX: an empty data line");
    x = v[0];
    ys.assign(v.begin() + 1, v.end());
    return;
  }
  size_t i = 0;
  while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
  {   // the abscissa, AFFN without exponent
    size_t j = i;
    if (j < s.size() && (s[j] == '+' || s[j] == '-')) ++j;
    while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j] == '.')) ++j;
    if (j == i) throw std::invalid_argument("JCAMP-DX: a data line without its abscissa: " + raw);
    x = std::strtod(s.substr(i, j - i).c_str(), nullptr);
    i = j;
  }
  double last = kNaN, diff = 0;
  bool in_dif = false;
  while (i < s.size()) {
    const char c = s[i];
    int d;
    if (c == ' ' || c == '\t' || c == ',') { ++i; continue; }
    if (c == '?') {
      ys.push_back(kNaN), last = kNaN, in_dif = false, ++i;
      continue;
    }
    if (sqz(c, d)) {
      ++i;
      last = with_digits(d, s, i), in_dif = false;
      ys.push_back(last);
    } else if (dif(c, d)) {
      ++i;
      diff = with_digits(d, s, i);
      last = last + diff, in_dif = true;
      ys.push_back(last);
    } else if (dup(c, d)) {
      ++i;
      const int count = int(with_digits(d, s, i));
      for (int k = 1; k < count; ++k) {
        if (in_dif) last += diff;
        ys.push_back(last);
      }
    } else if (c == '+' || c == '-' || c == '.' || std::isdigit(static_cast<unsigned char>(c))) {   // PAC / AFFN within the line
      size_t j = i + 1;
      while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j] == '.')) ++j;
      last = std::strtod(s.substr(i, j - i).c_str(), nullptr), in_dif = false;
      ys.push_back(last);
      i = j;
    } else {
      throw std::invalid_argument(std::string("JCAMP-DX: an unexpected character '") + c + "' in a data line");
    }
  }
  ends_dif = in_dif;
}

MeasuredCurve parse_jcamp(const std::string& text) {
  MeasuredCurve r;
  r.format = "JCAMP-DX";
  struct Ldr { std::string label, value; std::vector<std::string> data; };
  std::vector<Ldr> ldrs;
  for (const auto& l : lines_of(text)) {
    const std::string t = trim(l);
    if (t.rfind("##", 0) == 0) {
      const auto eq = t.find('=');
      std::string label = eq == std::string::npos ? t.substr(2) : t.substr(2, eq - 2), key;
      for (char c : upper(label))
        if (c != ' ' && c != '-' && c != '/' && c != '_') key += c;   // labels compare without blanks, hyphens, slashes, underscores
      std::string value = eq == std::string::npos ? "" : trim(t.substr(eq + 1));
      if (const auto cm = value.find("$$"); cm != std::string::npos) value = trim(value.substr(0, cm));
      ldrs.push_back({key, value, {}});
    } else if (!ldrs.empty() && !t.empty()) {
      ldrs.back().data.push_back(t);
    }
  }
  // the first block that holds data, with its own parameters
  size_t start = 0, data = ldrs.size();
  for (size_t k = 0; k < ldrs.size(); ++k) {
    if (ldrs[k].label == "TITLE") start = k;
    if (ldrs[k].label == "XYDATA" || ldrs[k].label == "XYPOINTS") { data = k; break; }
  }
  if (data == ldrs.size()) throw std::invalid_argument("JCAMP-DX: no ##XYDATA or ##XYPOINTS (peak tables and NMR blocks are not read)");
  double xf = 1, yf = 1, firstx = kNaN, lastx = kNaN, deltax = kNaN;
  long npoints = -1;
  auto num = [](const std::string& v) { return std::strtod(v.c_str(), nullptr); };
  for (size_t k = start; k < ldrs.size(); ++k) {
    const auto& L = ldrs[k];
    if (L.label == "TITLE") r.title = L.value;
    else if (L.label == "XUNITS") r.xunits = L.value;
    else if (L.label == "YUNITS") r.yunits = L.value;
    else if (L.label == "XFACTOR") xf = num(L.value);
    else if (L.label == "YFACTOR") yf = num(L.value);
    else if (L.label == "FIRSTX") firstx = num(L.value);
    else if (L.label == "LASTX") lastx = num(L.value);
    else if (L.label == "DELTAX") deltax = num(L.value);
    else if (L.label == "NPOINTS") npoints = std::lround(num(L.value));
    else if (L.label == "END") break;
  }
  int blocks = 0;
  for (const auto& L : ldrs) blocks += L.label == "TITLE";
  if (blocks > 1) r.notes.push_back(std::to_string(blocks) + " blocks in the file: the first with data read");
  const auto& D = ldrs[data];
  const std::string form = upper(D.value);
  if (D.label == "XYPOINTS") {
    std::string all;
    for (const auto& l : D.data) all += " " + (l.find("$$") == std::string::npos ? l : l.substr(0, l.find("$$")));
    const auto v = numbers(all);
    if (v.size() % 2) throw std::invalid_argument("JCAMP-DX: ##XYPOINTS needs x, y pairs");
    for (size_t k = 0; k + 1 < v.size(); k += 2) r.x.push_back(v[k] * xf), r.y.push_back(v[k + 1] * yf);
    return r;
  }
  if (form.find("X++(Y..Y)") == std::string::npos) throw std::invalid_argument("JCAMP-DX: ##XYDATA=" + D.value + " is not read (only (X++(Y..Y)))");
  std::vector<double> ys, line;
  double x = 0, first_line_x = kNaN;
  bool prev_dif = false;
  int checks = 0, bad_checks = 0;
  for (const auto& l : D.data) {
    bool ends_dif = false;
    decode_xydata_line(l, x, line, ends_dif);
    if (std::isnan(first_line_x)) first_line_x = x;
    size_t from = 0;
    if (prev_dif && !line.empty() && !ys.empty()) {   // the Y check value repeats the previous line's last ordinate
      ++checks;
      if (!(std::fabs(line[0] - ys.back()) <= 1e-9 * std::max(1.0, std::fabs(ys.back())))) ++bad_checks;
      from = 1;
    }
    ys.insert(ys.end(), line.begin() + long(from), line.end());
    prev_dif = ends_dif;
  }
  if (ys.empty()) throw std::invalid_argument("JCAMP-DX: no ordinates");
  if (bad_checks) r.notes.push_back(std::to_string(bad_checks) + " of " + std::to_string(checks) + " Y check values disagree: the data may be damaged");
  const size_t n = ys.size();
  if (npoints > 0 && size_t(npoints) != n) r.notes.push_back("##NPOINTS says " + std::to_string(npoints) + ", the table holds " + std::to_string(n));
  if (std::isnan(firstx)) firstx = first_line_x * xf;
  double step;
  if (!std::isnan(lastx) && n > 1) step = (lastx - firstx) / double(n - 1);
  else if (!std::isnan(deltax)) step = deltax;
  else throw std::invalid_argument("JCAMP-DX: neither ##LASTX nor ##DELTAX gives the abscissa spacing");
  for (size_t k = 0; k < n; ++k) r.x.push_back(firstx + double(k) * step), r.y.push_back(ys[k] * yf);
  return r;
}

// ---------------------------------------------------------------- XRDML

// the position of an element's start tag "<tag" or "<prefix:tag" (followed by blank, '>' or '/') from `from`
size_t find_tag(const std::string& t, const std::string& tag, size_t from, size_t to = std::string::npos) {
  for (size_t p = t.find(tag, from); p != std::string::npos && p < to; p = t.find(tag, p + 1)) {
    const size_t q = p + tag.size();
    if (q >= t.size() || !(t[q] == ' ' || t[q] == '>' || t[q] == '/' || t[q] == '\n' || t[q] == '\t' || t[q] == '\r')) continue;
    if (p > 0 && t[p - 1] == '<') return p - 1;
    if (p > 1 && t[p - 1] == ':') {
      size_t c = p - 1;
      while (c > 0 && (std::isalnum(static_cast<unsigned char>(t[c - 1])) || t[c - 1] == '_')) --c;
      if (c > 0 && t[c - 1] == '<') return c - 1;
    }
  }
  return std::string::npos;
}
// the text inside the element that starts at `at`, and where the element ends
std::string element_text(const std::string& t, size_t at, size_t* end = nullptr) {
  const size_t gt = t.find('>', at);
  if (gt == std::string::npos) throw std::invalid_argument("XRDML: an unclosed tag");
  if (t[gt - 1] == '/') { if (end) *end = gt + 1; return ""; }
  const size_t close = t.find("</", gt);
  if (close == std::string::npos) throw std::invalid_argument("XRDML: an unclosed element");
  if (end) *end = t.find('>', close) + 1;
  return t.substr(gt + 1, close - gt - 1);
}
std::string attribute(const std::string& t, size_t at, const std::string& name) {
  const size_t gt = t.find('>', at);
  const std::string head = t.substr(at, gt - at);
  const auto p = head.find(name + "=\"");
  if (p == std::string::npos) return "";
  const auto q = head.find('"', p + name.size() + 2);
  return head.substr(p + name.size() + 2, q - p - name.size() - 2);
}

MeasuredCurve parse_xrdml(const std::string& t) {
  MeasuredCurve r;
  r.format = "XRDML";
  r.xunits = "2θ (°)";
  size_t scans = 0;
  for (size_t p = find_tag(t, "scan", 0); p != std::string::npos; p = find_tag(t, "scan", p + 5)) ++scans;
  const size_t dp = find_tag(t, "dataPoints", 0);
  if (dp == std::string::npos) throw std::invalid_argument("XRDML: no <dataPoints>");
  // the end of the first dataPoints element: its closing tag
  size_t dp_end = t.find("dataPoints>", t.find('>', dp));
  if (dp_end == std::string::npos) throw std::invalid_argument("XRDML: <dataPoints> is not closed");
  std::vector<double> x;
  double start = kNaN, end = kNaN;
  for (size_t p = find_tag(t, "positions", dp, dp_end); p != std::string::npos; p = find_tag(t, "positions", p + 10, dp_end)) {
    if (attribute(t, p, "axis") != "2Theta") continue;
    const size_t pe = t.find("positions>", t.find('>', p));
    if (const size_t lp = find_tag(t, "listPositions", p, pe); lp != std::string::npos) {
      x = numbers(element_text(t, lp));
    } else {
      const size_t s = find_tag(t, "startPosition", p, pe), e = find_tag(t, "endPosition", p, pe);
      if (s != std::string::npos && e != std::string::npos)
        start = std::strtod(element_text(t, s).c_str(), nullptr), end = std::strtod(element_text(t, e).c_str(), nullptr);
      else
        throw std::invalid_argument("XRDML: 2Theta is held fixed (a scan of another axis)");
    }
    break;
  }
  if (x.empty() && std::isnan(start)) throw std::invalid_argument("XRDML: no 2Theta positions in the first scan");
  size_t yt = find_tag(t, "counts", dp, dp_end);
  const bool counts = yt != std::string::npos;
  if (!counts) yt = find_tag(t, "intensities", dp, dp_end);
  if (yt == std::string::npos) throw std::invalid_argument("XRDML: neither <counts> nor <intensities>");
  r.y = numbers(element_text(t, yt));
  const size_t n = r.y.size();
  if (n < 2) throw std::invalid_argument("XRDML: fewer than two data points");
  if (x.empty())   // the schema: step = (end − start)/(points − 1)
    for (size_t k = 0; k < n; ++k) x.push_back(start + (end - start) * double(k) / double(n - 1));
  if (x.size() != n) throw std::invalid_argument("XRDML: " + std::to_string(x.size()) + " positions for " + std::to_string(n) + " values");
  r.x = x;
  r.yunits = counts ? "counts (not corrected for attenuators)" : "intensities (counts × attenuation factor)";
  if (const size_t ct = find_tag(t, "commonCountingTime", dp, dp_end); ct != std::string::npos)
    r.notes.push_back("counting time " + trim(element_text(t, ct)) + " s per step");
  if (scans > 1) r.notes.push_back(std::to_string(scans) + " scans in the file: the first read");
  return r;
}

// ---------------------------------------------------------------- GSAS

MeasuredCurve parse_gsas(const std::string& text) {
  MeasuredCurve r;
  r.format = "GSAS";
  r.xunits = "2θ (°)";
  const auto L = lines_of(text);
  size_t b = 0;
  while (b < L.size() && L[b].rfind("BANK", 0) != 0) ++b;
  if (b == L.size()) throw std::invalid_argument("GSAS: no BANK record");
  if (!L.empty() && b > 0) r.title = trim(L[0]);
  std::istringstream hs(L[b]);
  std::string word, bintyp, type;
  long ibank = 0, nchan = 0, nrec = 0;
  double c1 = 0, c2 = 0, c3 = 0, c4 = 0;
  hs >> word >> ibank >> nchan >> nrec >> bintyp >> c1 >> c2 >> c3 >> c4;
  if (!(hs >> type)) type = "STD";
  type = upper(type), bintyp = upper(bintyp);
  if (nchan <= 0) throw std::invalid_argument("GSAS: the BANK record gives no points");
  if (bintyp.rfind("CONS", 0) != 0) throw std::invalid_argument("GSAS: BINTYP " + bintyp + " is not read (constant-wavelength CONS only)");
  int banks = 0;
  for (const auto& l : L) banks += l.rfind("BANK", 0) == 0;
  if (banks > 1) r.notes.push_back(std::to_string(banks) + " banks in the file: bank " + std::to_string(ibank) + " read");
  auto field = [](const std::string& line, size_t at, size_t w) { return at < line.size() ? trim(line.substr(at, w)) : std::string(); };
  size_t k = b + 1;
  if (type == "STD" || type == "ESD") {
    for (; k < L.size() && long(r.y.size()) < nchan; ++k) {
      if (L[k].rfind("BANK", 0) == 0) break;
      std::string line = L[k];
      if (!line.empty() && line.back() == '\x1a') line.pop_back();
      if (type == "STD") {   // 10 × (I2, F6.0): counters, intensity
        for (int f = 0; f < 10 && long(r.y.size()) < nchan; ++f) {
          const std::string nc = field(line, size_t(8 * f), 2), yo = field(line, size_t(8 * f + 2), 6);
          if (yo.empty()) continue;
          const double y = std::strtod(yo.c_str(), nullptr);
          const int ncv = nc.empty() ? 0 : std::atoi(nc.c_str());
          r.y.push_back(y), r.sigma.push_back(std::sqrt(std::max(0.0, y / (ncv > 0 ? ncv : 1))));
        }
      } else {   // 5 × (F8, F8): intensity, its σ
        for (int f = 0; f < 5 && long(r.y.size()) < nchan; ++f) {
          const std::string yo = field(line, size_t(16 * f), 8), e = field(line, size_t(16 * f + 8), 8);
          if (yo.empty()) continue;
          r.y.push_back(std::strtod(yo.c_str(), nullptr)), r.sigma.push_back(std::strtod(e.c_str(), nullptr));
        }
      }
    }
    for (long i = 0; i < long(r.y.size()); ++i) r.x.push_back((c1 + double(i) * c2) / 100.0);   // centidegrees
  } else if (type == "FXY" || type == "FXYE") {
    for (; k < L.size() && long(r.y.size()) < nchan; ++k) {
      if (L[k].rfind("BANK", 0) == 0) break;
      const auto v = numbers(L[k]);
      if (v.size() < 2) continue;
      r.x.push_back(v[0] / 100.0), r.y.push_back(v[1]);
      if (type == "FXYE" && v.size() > 2) r.sigma.push_back(v[2]);
    }
    if (r.sigma.size() != r.y.size()) r.sigma.clear();
  } else {
    throw std::invalid_argument("GSAS: TYPE " + type + " is not read (STD, ESD, FXY and FXYE are)");
  }
  if (long(r.y.size()) != nchan) r.notes.push_back("the BANK record says " + std::to_string(nchan) + " points, " + std::to_string(r.y.size()) + " read");
  r.yunits = "counts";
  return r;
}

// ---------------------------------------------------------------- UXD

MeasuredCurve parse_uxd(const std::string& text) {
  MeasuredCurve r;
  r.format = "UXD";
  r.xunits = "2θ (°)";
  double start = kNaN, step = kNaN;
  int mode = 0;   // 0 header, 1 single-column data, 2 two-column data, 3 a later range
  int ranges = 0;
  for (const auto& raw : lines_of(text)) {
    const std::string l = trim(raw);
    if (l.empty() || l[0] == ';') continue;
    if (l[0] == '_') {
      const auto eq = l.find('=');
      const std::string key = upper(trim(l.substr(0, eq)));
      if (key == "_DRIVE") ++ranges;
      if ((mode == 1 || mode == 2) && key == "_DRIVE") mode = 3;   // only the first range is read
      if (mode == 3) continue;
      if (eq != std::string::npos) {
        const std::string v = trim(l.substr(eq + 1));
        if (key == "_START") start = std::strtod(v.c_str(), nullptr);
        else if (key == "_STEPSIZE") step = std::strtod(v.c_str(), nullptr);
        else if (key == "_SAMPLE" && r.title.empty()) r.title = v;
        continue;
      }
      if (key.rfind("_2THETACOUNTS", 0) == 0 || key.rfind("_2THETACPS", 0) == 0) mode = 2, r.yunits = key.find("CPS") != std::string::npos ? "counts per second" : "counts";
      else if (key.rfind("_COUNT", 0) == 0 || key.rfind("_CPS", 0) == 0) mode = 1, r.yunits = key.find("CPS") != std::string::npos ? "counts per second" : "counts";
      continue;
    }
    if (mode != 1 && mode != 2) continue;
    if (!(std::isdigit(static_cast<unsigned char>(l[0])) || l[0] == '-' || l[0] == '+' || l[0] == '.')) continue;
    const auto v = numbers(l);
    if (mode == 1)
      for (double y : v) r.y.push_back(y);
    else if (v.size() >= 2)
      r.x.push_back(v[0]), r.y.push_back(v[1]);
  }
  if (r.y.empty()) throw std::invalid_argument("UXD: no _COUNTS, _CPS or _2THETACOUNTS data");
  if (r.x.empty()) {
    if (std::isnan(start) || std::isnan(step)) throw std::invalid_argument("UXD: single-column data need _START and _STEPSIZE");
    for (size_t k = 0; k < r.y.size(); ++k) r.x.push_back(start + double(k) * step);
  }
  if (ranges > 1) r.notes.push_back(std::to_string(ranges) + " ranges in the file: the first read");
  return r;
}

// ---------------------------------------------------------------- Rigaku *TYPE text

MeasuredCurve parse_rigaku(const std::string& text) {
  MeasuredCurve r;
  r.format = "Rigaku";
  r.xunits = "2θ (°)";
  double start = kNaN, stop = kNaN, step = kNaN;
  long count = -1;
  int groups = 0;
  bool in_group = false, reading = false, done = false;
  for (const auto& raw : lines_of(text)) {
    const std::string l = trim(raw);
    if (l.empty()) continue;
    if (l[0] == '*') {
      const auto eq = l.find('=');
      const std::string key = upper(trim(l.substr(0, eq)));
      const std::string v = eq == std::string::npos ? "" : trim(l.substr(eq + 1));
      if (key == "*BEGIN") { ++groups; in_group = !done; continue; }
      if (key == "*END") { if (reading) done = true; reading = false, in_group = false; continue; }
      if (!in_group || done) continue;
      if (key == "*START") start = std::strtod(v.c_str(), nullptr);
      else if (key == "*STOP") stop = std::strtod(v.c_str(), nullptr);
      else if (key == "*STEP") step = std::strtod(v.c_str(), nullptr);
      else if (key == "*COUNT") count = std::atol(v.c_str()), reading = true;
      continue;
    }
    if (reading)
      for (double y : numbers(l)) r.y.push_back(y);
  }
  if (r.y.empty()) throw std::invalid_argument("Rigaku: no *COUNT data in a *BEGIN … *END group");
  if (std::isnan(start) || std::isnan(step)) throw std::invalid_argument("Rigaku: the group gives no *START and *STEP");
  if (count >= 0 && long(r.y.size()) != count) r.notes.push_back("*COUNT says " + std::to_string(count) + ", " + std::to_string(r.y.size()) + " values read");
  for (size_t k = 0; k < r.y.size(); ++k) r.x.push_back(start + double(k) * step);
  if (!std::isnan(stop) && std::fabs(r.x.back() - stop) > 0.5 * std::fabs(step)) r.notes.push_back("the last point is not at *STOP");
  if (groups > 1) r.notes.push_back(std::to_string(groups) + " groups in the file: the first read");
  r.yunits = "as written (counts or counts per second; the file does not say)";
  return r;
}

// ---------------------------------------------------------------- columns

MeasuredCurve parse_columns(const std::string& text) {
  MeasuredCurve r;
  r.format = "columns";
  bool sig = true;
  for (const auto& raw : lines_of(text)) {
    if (raw.empty() || raw[0] == '#' || raw[0] == '!') continue;
    std::string line = raw;
    for (char& c : line) if (c == ',' || c == ';' || c == '\t') c = ' ';
    std::istringstream is(line);
    double a, b, s;
    if (!(is >> a >> b)) continue;
    r.x.push_back(a), r.y.push_back(b);
    if (is >> s) r.sigma.push_back(s);
    else sig = false;
  }
  if (!sig || r.sigma.size() != r.y.size()) r.sigma.clear();
  if (r.x.empty()) throw std::invalid_argument("no x y columns found");
  return r;
}

}  // namespace

MeasuredCurve parse_measured(const std::string& text, const std::string& name) {
  const std::string ext = extension(name);
  const std::string head = text.substr(0, 1024);
  size_t a = 0;
  while (a < text.size() && (std::isspace(static_cast<unsigned char>(text[a])) || text[a] == '\xef' || text[a] == '\xbb' || text[a] == '\xbf')) ++a;
  const std::string lead = text.substr(a, 16);
  bool uxd = false;
  for (const auto& l : lines_of(text.substr(0, 4096))) {
    const std::string t = trim(l);
    if (t.empty() || t[0] == ';') continue;
    uxd = t.rfind("_FILEVERSION", 0) == 0;
    break;
  }
  bool bank = false;
  for (const auto& l : lines_of(text.substr(0, 8192))) bank = bank || l.rfind("BANK ", 0) == 0;
  MeasuredCurve r;
  if (lead.rfind("##", 0) == 0 || ext == "JDX" || ext == "DX" || ext == "JCAMP") r = parse_jcamp(text);
  else if (head.find("www.xrdml.com") != std::string::npos || ext == "XRDML") r = parse_xrdml(text);
  else if (uxd || ext == "UXD") r = parse_uxd(text);
  else if (lead.rfind("*TYPE", 0) == 0) r = parse_rigaku(text);
  else if (bank || ext == "GSA" || ext == "GSS" || ext == "FXYE" || ext == "GSAS") r = parse_gsas(text);
  else r = parse_columns(text);
  if (r.x.size() != r.y.size()) throw std::invalid_argument(r.format + ": " + std::to_string(r.x.size()) + " abscissae for " + std::to_string(r.y.size()) + " ordinates");
  return r;
}

MeasuredCurve read_measured(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::invalid_argument("cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  const std::string ext = extension(path);
  if (ext == "RAW" && b.str().rfind("RAW", 0) == 0)
    throw std::invalid_argument(path + ": a Bruker binary .raw is not read (no public specification): export it as .uxd, .xrdml or columns");
  if (ext == "XYD") throw std::invalid_argument(path + ": STOE .xyd is not read (no public description of the format): export it as columns");
  return parse_measured(b.str(), path);
}

// ---------------------------------------------------------------- Gaussian cube

GridField parse_cube(const std::string& text, System* atoms) {
  std::istringstream in(text);
  std::string t1, t2, line;
  std::getline(in, t1), std::getline(in, t2);
  GridField g;
  g.name = trim(t1).empty() ? "cube" : trim(t1);
  g.unit = "as in the file";
  auto need = [&](std::istream& s, auto& v) { if (!(s >> v)) throw std::invalid_argument("cube: the header ends early"); };
  long natoms = 0;
  double o[3];
  std::getline(in, line);
  long nval = 1;
  {
    std::istringstream h(line);
    need(h, natoms), need(h, o[0]), need(h, o[1]), need(h, o[2]);
    long nv;
    if (h >> nv && nv > 0) nval = nv;
  }
  Vec3 ax[3];
  bool angstrom = false;
  for (int k = 0; k < 3; ++k) {
    std::getline(in, line);
    std::istringstream h(line);
    long n;
    need(h, n), need(h, ax[k][0]), need(h, ax[k][1]), need(h, ax[k][2]);
    if (n == 0) throw std::invalid_argument("cube: a zero voxel count");
    if (n < 0) angstrom = true;
    g.n[k] = int(std::labs(n));
  }
  const double unit = angstrom ? 1.0 : kBohr;
  for (auto& v : ax) v = v * unit;
  const Vec3 first{o[0] * unit, o[1] * unit, o[2] * unit};
  System s;
  for (long a = 0; a < std::labs(natoms); ++a) {
    std::getline(in, line);
    std::istringstream h(line);
    double z, q, x, y, w;
    need(h, z), need(h, q), need(h, x), need(h, y), need(h, w);
    Atom at;
    at.element = int(std::lround(z));
    at.pos = {x * unit, y * unit, w * unit};
    s.atoms.push_back(at);
  }
  if (natoms < 0) {   // NMO and the orbital numbers (10I5 records): NMO + 1 integers
    long nmo;
    need(in, nmo);
    for (long k = 0; k < nmo; ++k) {
      long m;
      need(in, m);
    }
    nval = std::max(1L, nmo);
  }
  const size_t N = size_t(g.n[0]) * size_t(g.n[1]) * size_t(g.n[2]);
  g.values.resize(N);
  const std::string rest((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  const char* p = rest.c_str();
  for (size_t q = 0; q < N; ++q)
    for (long v = 0; v < nval; ++v) {
      char* e = nullptr;
      const double x = std::strtod(p, &e);
      if (e == p) throw std::invalid_argument("cube: " + std::to_string(q) + " of " + std::to_string(N) + " values read");
      if (v == 0) g.values[q] = x;
      p = e;
    }
  g.cell.a = ax[0] * double(g.n[0]), g.cell.b = ax[1] * double(g.n[1]), g.cell.c = ax[2] * double(g.n[2]);
  g.cell.origin = first - (ax[0] + ax[1] + ax[2]) * 0.5;
  if (atoms) {
    s.cell = g.cell;
    s.title = g.name;
    *atoms = s;
  }
  return g;
}

GridField read_cube(const std::string& path, System* atoms) {
  std::ifstream in(path);
  if (!in) throw std::invalid_argument("cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  return parse_cube(b.str(), atoms);
}

Json measured_json(const MeasuredCurve& c) {
  Json j = Json::object();
  auto arr = [](const std::vector<double>& v) { Json a = Json::array(); for (double x : v) a.push_back(std::isnan(x) ? Json() : Json(x)); return a; };
  j["format"] = c.format, j["title"] = c.title, j["xunits"] = c.xunits, j["yunits"] = c.yunits;
  j["x"] = arr(c.x), j["y"] = arr(c.y), j["sigma"] = arr(c.sigma);
  Json n = Json::array();
  for (const auto& x : c.notes) n.push_back(x);
  j["notes"] = n;
  return j;
}

}  // namespace caps
