// Release manifests and the update decision (caps/update.hpp).
#include "caps/update.hpp"

#include "caps/config.hpp"
#include "caps/json.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace caps {

const char* const kReleaseManifestUrl = "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest/download/release.json";
namespace {
const char* const kReleasesPage = "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest";
const char* const kDownloadBase = "https://github.com/MuhammadUzairRiaz/CAPS/releases/latest/download/";

std::string lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
std::string trim(const std::string& s) {
  size_t a = 0, b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return s.substr(a, b - a);
}
// days since 1970-01-01 of a proleptic Gregorian date (Howard Hinnant's days_from_civil)
long long days_from_civil(long long y, unsigned m, unsigned d) {
  y -= m <= 2;
  const long long era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return era * 146097 + static_cast<long long>(doe) - 719468;
}
}  // namespace

ReleaseManifest parse_release_manifest(const std::string& json) {
  Json j;
  try {
    j = Json::parse(json);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("release manifest: not JSON (") + e.what() + ")");
  }
  if (!j.is_object()) throw std::runtime_error("release manifest: not a JSON object");
  const std::string format = j.text("format");
  if (format.rfind("caps-release/1", 0) != 0) throw std::runtime_error("release manifest: format is '" + format + "', not caps-release/1");
  ReleaseManifest m;
  m.version = trim(j.text("version"));
  m.commit = lower(trim(j.text("commit")));
  m.date = trim(j.text("date"));
  m.title = j.text("title");
  m.notes_url = j.text("notes_url");
  if (m.version.empty() || m.commit.empty() || m.date.empty()) throw std::runtime_error("release manifest: version, commit and date are required");
  if (j.has("highlights") && j["highlights"].is_array())
    for (const auto& h : j["highlights"].items()) {
      if (h.is_string()) { m.highlights.push_back({h.str(), "", ""}); continue; }
      if (!h.is_object()) continue;
      ReleaseHighlight r{h.text("title"), h.text("detail"), h.text("page")};
      if (!r.title.empty()) m.highlights.push_back(r);
    }
  if (j.has("assets") && j["assets"].is_object())
    for (const auto& [key, a] : j["assets"].members()) {
      if (a.is_string()) m.assets[key] = {a.str(), ""};
      else if (a.is_object() && !a.text("file").empty()) m.assets[key] = {a.text("file"), lower(a.text("sha256"))};
    }
  return m;
}

std::string release_manifest_json(const ReleaseManifest& m) {
  Json j = Json::object();
  j["format"] = "caps-release/1";
  j["version"] = m.version;
  j["commit"] = m.commit;
  j["date"] = m.date;
  j["title"] = m.title;
  Json hs = Json::array();
  for (const auto& h : m.highlights) {
    Json o = Json::object();
    o["title"] = h.title;
    if (!h.detail.empty()) o["detail"] = h.detail;
    if (!h.page.empty()) o["page"] = h.page;
    hs.push_back(o);
  }
  j["highlights"] = hs;
  j["notes_url"] = m.notes_url;
  Json as = Json::object();
  for (const auto& [k, a] : m.assets) {
    Json o = Json::object();
    o["file"] = a.file;
    if (!a.sha256.empty()) o["sha256"] = a.sha256;
    as[k] = o;
  }
  j["assets"] = as;
  return j.dump(1);
}

int compare_versions(const std::string& a, const std::string& b) {
  auto split = [](const std::string& v, std::vector<long long>& nums, std::string& pre) {
    std::string s = trim(v);
    if (!s.empty() && (s[0] == 'v' || s[0] == 'V')) s.erase(0, 1);
    const auto dash = s.find_first_of("-+");
    if (dash != std::string::npos) {
      if (s[dash] == '-') pre = s.substr(dash + 1);
      const auto plus = pre.find('+');   // build metadata never orders
      if (plus != std::string::npos) pre.erase(plus);
      s.erase(dash);
    }
    std::stringstream ss(s);
    std::string part;
    while (std::getline(ss, part, '.')) {
      long long n = 0;
      for (char c : part) {
        if (!std::isdigit(static_cast<unsigned char>(c))) break;
        n = n * 10 + (c - '0');
      }
      nums.push_back(n);
    }
  };
  std::vector<long long> na, nb;
  std::string pa, pb;
  split(a, na, pa);
  split(b, nb, pb);
  const size_t n = std::max(na.size(), nb.size());
  for (size_t k = 0; k < n; ++k) {
    const long long x = k < na.size() ? na[k] : 0, y = k < nb.size() ? nb[k] : 0;
    if (x != y) return x < y ? -1 : 1;
  }
  if (pa == pb) return 0;
  if (pa.empty()) return 1;   // 0.2.0 > 0.2.0-rc1
  if (pb.empty()) return -1;
  return pa < pb ? -1 : 1;
}

long long parse_iso_time(const std::string& raw) {
  const std::string s = trim(raw);
  int y = 0, mo = 0, d = 0, h = 0, mi = 0, se = 0;
  if (s.size() < 10 || std::sscanf(s.c_str(), "%4d-%2d-%2d", &y, &mo, &d) != 3) return -1;
  if (mo < 1 || mo > 12 || d < 1 || d > 31) return -1;
  size_t k = 10;
  long long offset = 0;
  if (k < s.size() && (s[k] == 'T' || s[k] == ' ')) {
    if (std::sscanf(s.c_str() + k + 1, "%2d:%2d", &h, &mi) != 2) return -1;
    k += 6;
    if (k < s.size() && s[k] == ':') {
      if (std::sscanf(s.c_str() + k + 1, "%2d", &se) != 1) return -1;
      k += 3;
    }
    if (k < s.size() && s[k] == '.')
      for (++k; k < s.size() && std::isdigit(static_cast<unsigned char>(s[k])); ++k) {}
    if (k < s.size() && (s[k] == '+' || s[k] == '-')) {
      int oh = 0, om = 0;
      if (std::sscanf(s.c_str() + k + 1, "%2d:%2d", &oh, &om) < 1 && std::sscanf(s.c_str() + k + 1, "%2d%2d", &oh, &om) < 1) return -1;
      offset = (s[k] == '+' ? 1 : -1) * (oh * 3600LL + om * 60LL);
    }
  }
  return days_from_civil(y, static_cast<unsigned>(mo), static_cast<unsigned>(d)) * 86400LL + h * 3600LL + mi * 60LL + se - offset;
}

bool same_commit(const std::string& a, const std::string& b) {
  const std::string x = lower(trim(a)), y = lower(trim(b));
  if (x.size() < 7 || y.size() < 7 || x == "unknown" || y == "unknown") return false;
  const size_t n = std::min(x.size(), y.size());
  return x.compare(0, n, y, 0, n) == 0;
}

BuildIdentity this_build() { return {version_string(), commit_string(), commit_date_string()}; }

UpdateVerdict decide_update(const BuildIdentity& local, const ReleaseManifest& remote) {
  const int v = compare_versions(remote.version, local.version);
  if (v > 0) return UpdateVerdict::Update;
  if (v < 0) return UpdateVerdict::LocalNewer;
  if (same_commit(local.commit, remote.commit)) return UpdateVerdict::Current;
  const long long tr = parse_iso_time(remote.date), tl = parse_iso_time(local.date);
  if (tr < 0 || tl < 0) return UpdateVerdict::Unknown;
  return tr > tl ? UpdateVerdict::Update : UpdateVerdict::LocalNewer;
}

const char* verdict_name(UpdateVerdict v) {
  switch (v) {
    case UpdateVerdict::Update: return "update";
    case UpdateVerdict::Current: return "current";
    case UpdateVerdict::LocalNewer: return "newer";
    default: return "unknown";
  }
}

bool update_checks_disabled() {
  const char* e = std::getenv("CAPS_NO_UPDATE_CHECK");
  if (!e) return false;
  const std::string v = lower(trim(e));
  return !(v.empty() || v == "0" || v == "false" || v == "no" || v == "off");
}

std::string update_feed() {
  const char* f = std::getenv("CAPS_UPDATE_FEED");
  return f && *f ? std::string(f) : std::string(kReleaseManifestUrl);
}

std::string fetch_release_manifest(const std::string& feed, int timeout_s) {
  const bool url = feed.rfind("https://", 0) == 0 || feed.rfind("http://", 0) == 0;
  if (!url) {
    std::ifstream f(feed, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read " + feed);
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
  }
  // the URL is ours or the tester's CAPS_UPDATE_FEED; it goes to curl as one quoted word
  for (char c : feed)
    if (c == '"' || c == '\'' || c == '`' || c == '$' || c == '\\' || std::isspace(static_cast<unsigned char>(c)) || c == '%' || c == '!')
      throw std::runtime_error("the update feed URL has characters CAPS does not pass on");
  const std::string proto = feed.rfind("https://", 0) == 0 ? "=https" : "=http,https";
  std::string cmd = "curl -q -sSfL --proto " + proto + " --proto-redir =https --max-redirs 5 --connect-timeout " +
                    std::to_string(std::max(1, timeout_s - 1)) + " --max-time " + std::to_string(std::max(1, timeout_s)) +
                    " -A CAPS -H \"Accept: application/json, application/octet-stream\" \"" + feed + "\"";
#ifdef _WIN32
  cmd += " 2>NUL";
  FILE* p = _popen(cmd.c_str(), "rb");
#else
  cmd += " 2>/dev/null";
  FILE* p = popen(cmd.c_str(), "r");
#endif
  if (!p) throw std::runtime_error("cannot start curl");
  std::string out;
  char buf[4096];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof buf, p)) > 0 && out.size() < (1u << 20)) out.append(buf, n);
#ifdef _WIN32
  const int rc = _pclose(p);
#else
  const int rc = pclose(p);
#endif
  if (rc != 0 || out.empty()) throw std::runtime_error("offline, or the release manifest could not be fetched");
  return out;
}

std::string platform_key() {
#if defined(__APPLE__)
#if defined(__aarch64__) || defined(__arm64__)
  return "macos-arm64";
#else
  return "macos-x86_64";
#endif
#elif defined(_WIN32)
  return "windows-x64";
#else
#if defined(__aarch64__)
  return "linux-aarch64";
#else
  return "linux-x86_64";
#endif
#endif
}

std::string download_url(const ReleaseManifest& m, const std::string& key) {
  auto it = m.assets.find(key);
  if (it == m.assets.end()) {
    // "linux-x86_64-deb" without a deb: the platform's main file
    const auto dash = key.find_last_of('-');
    if (dash != std::string::npos && key.find('-') != dash) it = m.assets.find(key.substr(0, dash));
  }
  if (it != m.assets.end() && !it->second.file.empty()) return kDownloadBase + it->second.file;
  return m.notes_url.empty() ? std::string(kReleasesPage) : m.notes_url;
}

}  // namespace caps
