// CAPS live output (caps/live.hpp): files that grow while a run goes, each line flushed.
#include "caps/live.hpp"

#include "caps/config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <ctime>
#include <stdexcept>

#ifndef CAPS_VERSION_STRING
#define CAPS_VERSION_STRING "0.1.0"
#endif
#ifndef CAPS_COMMIT_STRING
#define CAPS_COMMIT_STRING "unknown"
#endif

namespace caps {

const char* version_string() { return CAPS_VERSION_STRING; }
const char* commit_string() { return CAPS_COMMIT_STRING; }

namespace {
int positive_int(const char* s) {
  if (!s || !*s) return 0;
  char* end = nullptr;
  const long v = std::strtol(s, &end, 10);
  return end && *end == '\0' && v >= 1 && v <= 100000 ? int(v) : 0;
}

std::unique_ptr<std::ofstream> open_file(const std::string& path) {
  const std::filesystem::path parent = std::filesystem::path(path).parent_path();
  std::error_code ec;
  if (!parent.empty()) std::filesystem::create_directories(parent, ec);
  auto f = std::make_unique<std::ofstream>(path, std::ios::out | std::ios::trunc);
  if (!*f) throw std::runtime_error("cannot write " + path);
  return f;
}

std::string iso_now() {
  const std::time_t t = std::time(nullptr);
  char b[32];
  std::strftime(b, sizeof b, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  return b;
}

// a finite number for JSON (NaN and ±∞ would make the line unreadable)
double fin(double v) { return std::isfinite(v) ? v : 0.0; }
}  // namespace

int env_threads() {
  if (const int n = positive_int(std::getenv("CAPS_THREADS"))) return n;
  return positive_int(std::getenv("SLURM_CPUS_PER_TASK"));
}

std::string default_progress_file() {
  const char* job = std::getenv("SLURM_JOB_ID");
  const char* caps_job = std::getenv("CAPS_JOB");
  return (job && *job) || (caps_job && *caps_job) ? "progress.jsonl" : "";
}

std::string thermo_header() {
  return "step,time_ps,stage,target_K,temperature_K,potential,kinetic,total,conserved,pressure_atm,"
         "pxx_atm,pyy_atm,pzz_atm,pxy_atm,pxz_atm,pyz_atm,lx_A,ly_A,lz_A,volume_A3,density_g_cm3";
}

std::string thermo_csv_row(const ThermoRow& r, const std::string& stage) {
  std::string st = stage;
  std::replace(st.begin(), st.end(), ',', ';');
  std::replace(st.begin(), st.end(), '\n', ' ');
  char b[512];
  std::snprintf(b, sizeof b, "%lld,%.4f,%s,%.2f,%.3f,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f,%.2f,%.5f",
                static_cast<long long>(r.step), r.time_ps, st.c_str(), r.target_temperature, r.temperature, r.potential, r.kinetic, r.total,
                r.conserved, r.pressure, r.p[0], r.p[1], r.p[2], r.p[3], r.p[4], r.p[5], r.lx, r.ly, r.lz, r.volume, r.density);
  return b;
}

void LiveOutput::open_thermo(const std::string& path) {
  thermo_ = open_file(path);
  *thermo_ << thermo_header() << "\n";
  thermo_->flush();
}

void LiveOutput::open_frames(const std::string& path, const System& topology) {
  frames_ = open_file(path);
  topology_ = topology;
  nframes_ = 0;
}

void LiveOutput::open_progress(const std::string& path) {
  progress_ = open_file(path);
  Json j = Json::object();
  j["command"] = command_;
  j["state"] = std::string("running");
  j["start"] = iso_now();
  *progress_ << j.dump(0) << "\n";
  progress_->flush();
}

void LiveOutput::thermo(const ThermoRow& r, const std::string& stage) {
  if (!thermo_) return;
  *thermo_ << thermo_csv_row(r, stage) << "\n";
  thermo_->flush();
}

void LiveOutput::frame(const std::vector<double>& x, const Cell& c, int64_t step) {
  if (!frames_) return;
  std::vector<Vec3> p(x.size() / 3);
  for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
  write_lammps_dump_frame(*frames_, topology_, p, c, step);
  frames_->flush();
  ++nframes_;
}

void LiveOutput::progress(const std::string& stage, int stage_index, int stages, double fraction, const ThermoRow* r, const Json& extra) {
  if (!progress_) return;
  Json j = Json::object();
  j["command"] = command_;
  j["stage"] = stage;
  j["stage_index"] = double(stage_index);
  j["stages"] = double(stages);
  if (fraction >= 0) j["fraction"] = fin(std::min(1.0, fraction));
  if (r) {
    j["step"] = double(r->step);
    j["time_ps"] = fin(r->time_ps);
    j["T"] = fin(r->temperature);
    if (r->target_temperature > 0) j["T_target"] = fin(r->target_temperature);
    j["P"] = fin(r->pressure);
    j["rho"] = fin(r->density);
    j["pe"] = fin(r->potential);
    j["ke"] = fin(r->kinetic);
    j["etotal"] = fin(r->total);
    j["volume"] = fin(r->volume);
    j["lx"] = fin(r->lx), j["ly"] = fin(r->ly), j["lz"] = fin(r->lz);
  }
  if (extra.is_object())
    for (const auto& [k, v] : extra.members()) j[k] = v;
  *progress_ << j.dump(0) << "\n";
  progress_->flush();
}

void LiveOutput::done(bool ok, const std::string& note) {
  if (!progress_) return;
  Json j = Json::object();
  j["command"] = command_;
  j["state"] = std::string(ok ? "finished" : "failed");
  j["end"] = iso_now();
  if (!note.empty()) j["note"] = note;
  *progress_ << j.dump(0) << "\n";
  progress_->flush();
}

std::vector<Json> read_progress_lines(const std::string& text, size_t* consumed) {
  std::vector<Json> out;
  size_t pos = 0, used = 0;
  while (pos < text.size()) {
    const size_t nl = text.find('\n', pos);
    if (nl == std::string::npos) break;   // a line still being written
    const std::string line = text.substr(pos, nl - pos);
    pos = used = nl + 1;
    if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
    try {
      out.push_back(Json::parse(line));
    } catch (const std::exception&) {
      // a damaged line (a file truncated and rewritten): skipped
    }
  }
  if (consumed) *consumed = used;
  return out;
}

}  // namespace caps
