// CAPS jobs on a cluster or workstation (caps/cluster_job.hpp): profiles, the job script, folders, the state file,
// SLURM output and the poll a viewer reads.
#include "caps/cluster_job.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/config.hpp"
#include "caps/live.hpp"

namespace caps {

namespace fs = std::filesystem;

namespace {
std::string slurp(const std::string& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}
void write_file(const std::string& p, const std::string& text) {
  std::ofstream f(p, std::ios::binary | std::ios::trunc);
  if (!f) throw std::runtime_error("cannot write " + p);
  f << text;
}
std::string iso_now() {
  const std::time_t t = std::time(nullptr);
  char b[32];
  std::strftime(b, sizeof b, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  return b;
}
std::string trim(std::string s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const auto b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}
std::string str_or(const Json& j, const char* k, const std::string& d) { return j.has(k) && j[k].is_string() ? j[k].str() : d; }
double num_or(const Json& j, const char* k, double d) { return j.has(k) && j[k].is_number() ? j[k].number() : d; }
bool bool_or(const Json& j, const char* k, bool d) { return j.has(k) && j[k].kind() == Json::Bool ? j[k].boolean() : d; }

// "2G" per core × cores → "32gb" (PBS asks for the node's total)
std::string total_mem(const std::string& per_cpu, int cpus) {
  if (per_cpu.empty()) return "";
  double v = std::atof(per_cpu.c_str());
  const char u = char(std::toupper(static_cast<unsigned char>(per_cpu.back())));
  if (u == 'M') v /= 1024;
  if (u == 'T') v *= 1024;
  if (v <= 0) return "";
  char b[32];
  std::snprintf(b, sizeof b, "%.0fgb", std::ceil(v * std::max(1, cpus)));
  return b;
}

// the next free <kind>-<n> in a folder
std::string next_job_name(const fs::path& parent, const std::string& kind) {
  int n = 0;
  std::error_code ec;
  if (fs::exists(parent, ec))
    for (const auto& e : fs::directory_iterator(parent, ec)) {
      const std::string name = e.path().filename().string();
      if (name.rfind(kind + "-", 0) != 0) continue;
      const std::string tail = name.substr(kind.size() + 1);
      if (!tail.empty() && std::all_of(tail.begin(), tail.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); }))
        n = std::max(n, std::atoi(tail.c_str()));
    }
  return kind + "-" + std::to_string(n + 1);
}

// "caps md …" → exec "$CAPS" md … (cmd.sh, run in the scratch folder with CAPS set by job.sh)
std::string command_script(const std::string& command) {
  // every caps that starts a command (after &&, || or ;) is the installed program; one command runs in place (exec)
  const std::string c = trim(command);
  std::string out;
  bool at_start = true, chained = false;
  for (size_t i = 0; i < c.size();) {
    if (at_start) {
      size_t j = i;
      while (j < c.size() && c[j] == ' ') ++j;
      if (c.compare(j, 5, "caps ") == 0 || (j + 4 == c.size() && c.compare(j, 4, "caps") == 0)) { out += c.substr(i, j - i) + "\"$CAPS\""; i = j + 4; at_start = false; continue; }
      at_start = false;
    }
    if (c.compare(i, 2, "&&") == 0 || c.compare(i, 2, "||") == 0) { out += c.substr(i, 2); i += 2; at_start = chained = true; continue; }
    if (c[i] == ';') { out += ';'; ++i; at_start = chained = true; continue; }
    out += c[i++];
  }
  return "#!/bin/bash\n# the command of this CAPS job (cmd.txt as written; caps is the installed program, $CAPS)\n" + std::string(chained ? "" : "exec ") + out + "\n";
}

void write_task(const fs::path& dir, const std::string& title, const std::string& job, const std::string& command,
                const std::vector<std::string>& inputs, const std::string& array) {
  fs::create_directories(dir);
  for (const auto& in : inputs) {
    if (in.empty()) continue;
    if (!fs::exists(in)) throw std::invalid_argument("input not found: " + in);
    const fs::path to = dir / fs::path(in).lexically_normal().filename();
    std::error_code ec;
    if (fs::is_directory(in)) {   // a folder the command reads (a coarse-graining stage's tables): copied whole
      fs::copy(in, to, fs::copy_options::recursive | fs::copy_options::overwrite_existing);
      continue;
    }
    // a hard link where the file system allows it (a resumed job's trajectory can be large), else a copy
    fs::remove(to, ec);
    fs::create_hard_link(in, to, ec);
    if (ec) fs::copy_file(in, to, fs::copy_options::overwrite_existing);
  }
  write_file((dir / "cmd.txt").string(), trim(command) + "\n");
  write_file((dir / "cmd.sh").string(), command_script(command));
  fs::permissions(dir / "cmd.sh", fs::perms::owner_exec | fs::perms::group_exec, fs::perm_options::add);
  Json st = Json::object();
  st["title"] = title;
  st["job"] = job;
  st["state"] = std::string("created");
  st["created"] = iso_now();
  if (!array.empty()) st["array"] = array;
  write_job_state(dir.string(), st);
}
}  // namespace

// ------------------------------------------------------------------------------------------------- profiles
Json host_profile_json(const HostProfile& h) {
  Json j = Json::object();
  j["format"] = std::string("caps-host");
  j["version"] = 1.0;
  j["name"] = h.name;
  j["ssh"] = h.ssh;
  j["port"] = double(h.port);
  j["scheduler"] = h.scheduler;
  j["account"] = h.account;
  j["partition"] = h.partition;
  j["modules"] = h.modules;
  j["build_modules"] = h.build_modules;
  j["root"] = h.root;
  j["scratch"] = h.scratch;
  j["ws_filesystem"] = h.ws_filesystem;
  j["ws_days"] = double(h.ws_days);
  j["scratch_var"] = h.scratch_var;
  j["cpus"] = double(h.cpus);
  j["mem_per_cpu"] = h.mem_per_cpu;
  j["time"] = h.time;
  j["mail_type"] = h.mail_type;
  j["mail_user"] = h.mail_user;
  j["keep_workspace"] = h.keep_workspace;
  j["node_cores"] = double(h.node_cores);
  j["node_mem"] = h.node_mem;
  if (!h.job_template.empty()) j["job_template"] = h.job_template;
  return j;
}

HostProfile host_profile_from_json(const Json& j) {
  HostProfile h;
  if (!j.is_object()) return h;
  h.name = str_or(j, "name", h.name);
  h.ssh = str_or(j, "ssh", h.ssh);
  h.port = int(num_or(j, "port", h.port));
  h.scheduler = str_or(j, "scheduler", h.scheduler);
  h.account = str_or(j, "account", h.account);
  h.partition = str_or(j, "partition", h.partition);
  h.modules = str_or(j, "modules", h.modules);
  h.build_modules = str_or(j, "build_modules", h.build_modules);
  h.root = str_or(j, "root", h.root);
  h.scratch = str_or(j, "scratch", h.scratch);
  h.ws_filesystem = str_or(j, "ws_filesystem", h.ws_filesystem);
  h.ws_days = int(num_or(j, "ws_days", h.ws_days));
  h.scratch_var = str_or(j, "scratch_var", h.scratch_var);
  h.cpus = int(num_or(j, "cpus", h.cpus));
  h.mem_per_cpu = str_or(j, "mem_per_cpu", h.mem_per_cpu);
  h.time = str_or(j, "time", h.time);
  h.mail_type = str_or(j, "mail_type", h.mail_type);
  h.mail_user = str_or(j, "mail_user", h.mail_user);
  h.keep_workspace = bool_or(j, "keep_workspace", h.keep_workspace);
  h.node_cores = int(num_or(j, "node_cores", h.node_cores));
  h.node_mem = str_or(j, "node_mem", h.node_mem);
  h.job_template = str_or(j, "job_template", h.job_template);
  if (h.scheduler != "slurm" && h.scheduler != "pbs" && h.scheduler != "direct") h.scheduler = "slurm";
  if (h.scratch != "workspace" && h.scratch != "env" && h.scratch != "none") h.scratch = "none";
  if (h.cpus < 1) h.cpus = 1;
  return h;
}

std::vector<std::pair<std::string, HostProfile>> host_presets() {
  std::vector<std::pair<std::string, HostProfile>> out;
  HostProfile slurm;
  slurm.name = "SLURM cluster";
  slurm.modules = "module purge\nmodule load GCC";
  slurm.build_modules = "module purge\nmodule load GCC CMake";
  out.push_back({"slurm", slurm});
  HostProfile ws = slurm;   // SLURM sites with workspaces (ws_allocate / ws_release)
  ws.name = "SLURM with workspaces";
  ws.scratch = "workspace";
  out.push_back({"slurm-workspace", ws});
  HostProfile pbs;
  pbs.name = "PBS cluster";
  pbs.scheduler = "pbs";
  pbs.modules = "module purge\nmodule load gcc";
  pbs.build_modules = pbs.modules + " cmake";
  pbs.scratch = "env";
  pbs.scratch_var = "$TMPDIR";
  out.push_back({"pbs", pbs});
  HostProfile ws2;
  ws2.name = "Workstation";
  ws2.scheduler = "direct";
  ws2.cpus = 8;
  ws2.time.clear();
  ws2.mem_per_cpu.clear();
  out.push_back({"workstation", ws2});
  return out;
}

std::string expand_home(const std::string& path) {
  if (path.empty() || path[0] != '~') return path;
  const char* home = std::getenv("HOME");
#ifdef _WIN32
  if (!home) home = std::getenv("USERPROFILE");
#endif
  return (home ? std::string(home) : std::string(".")) + path.substr(1);
}

std::string default_host_profile_path() { return expand_home("~/CAPS/host.json"); }

HostProfile load_host_profile(const std::string& path) {
  if (path.empty() || !fs::exists(path)) return HostProfile{};
  return host_profile_from_json(Json::parse(slurp(path)));
}

void save_host_profile(const HostProfile& h, const std::string& path) {
  const fs::path p(path);
  if (!p.parent_path().empty()) fs::create_directories(p.parent_path());
  write_file(path, host_profile_json(h).dump(2) + "\n");
}

std::string sanitize_title(const std::string& title) {
  std::string s;
  for (char c : trim(title)) s += std::isalnum(static_cast<unsigned char>(c)) || c == '.' || c == '_' || c == '-' ? c : '_';
  while (!s.empty() && s[0] == '.') s.erase(0, 1);   // no hidden folders
  if (s.size() > 80) s.resize(80);
  return s.empty() ? std::string("job") : s;
}

// ------------------------------------------------------------------------------------------------- the job script
namespace {
const char* kSlurmHead =
    "#!/bin/bash\n"
    "#SBATCH --job-name=caps-{job}\n"
    "#SBATCH --account={account}\n"
    "#SBATCH --partition={partition}\n"
    "#SBATCH --nodes=1\n"
    "#SBATCH --ntasks=1\n"
    "#SBATCH --cpus-per-task={cpus}\n"
    "#SBATCH --mem-per-cpu={mem}\n"
    "#SBATCH --time={time}\n"
    "#SBATCH --chdir={jobdir}\n"
    "#SBATCH --output=slurm-%j.out\n"
    "#SBATCH --error=slurm-%j.err\n"
    "#SBATCH --signal=B:USR1@600\n"
    "{array}\n"
    "{mail_lines}\n";
const char* kPbsHead =
    "#!/bin/bash\n"
    "#PBS -N caps-{job}\n"
    "#PBS -A {account}\n"
    "#PBS -q {partition}\n"
    "#PBS -l select=1:ncpus={cpus}:mem={mem_total}\n"
    "#PBS -l walltime={time}\n"
    "#PBS -o {jobdir}/pbs.out\n"
    "#PBS -e {jobdir}/pbs.err\n"
    "{array}\n"
    "{mail_lines}\n";
const char* kDirectHead = "#!/bin/bash\n# a workstation: run in the background by caps job submit (nohup), no queue\n";
const char* kBody =
    "# CAPS job {title} / {job} (written by caps job new, CAPS {version}); edit freely, the Studio's Reset restores it.\n"
    "# CAPS runs threads, not MPI ranks: one node, one task, {cpus} cores.\n"
    "set -u\n"
    "JOBDIR={jobdir}\n"
    "TASK=${SLURM_ARRAY_TASK_ID:-${PBS_ARRAY_INDEX:-${CAPS_TASK:-}}}\n"
    "if [ -n \"$TASK\" ] && [ -f \"$JOBDIR/tasks.txt\" ]; then JOBDIR=$(sed -n \"$((TASK + 1))p\" \"$JOBDIR/tasks.txt\"); fi\n"
    "TITLE=$(basename \"$(dirname \"$JOBDIR\")\"); JOB=$(basename \"$JOBDIR\")\n"
    "JOBKEY=${SLURM_JOB_ID:-${PBS_JOBID:-local-$$}}; JOBKEY=${JOBKEY%%.*}\n"
    "export OMP_NUM_THREADS=1\n"
    "export CAPS_THREADS=${SLURM_CPUS_PER_TASK:-${NCPUS:-{cpus}}}\n"
    "export CAPS_JOB=$JOB\n"
    "export CAPS={caps}\n"
    "KEEP={keep_workspace}\n"
    "{modules}\n"
    "cd \"$JOBDIR\" || exit 1\n"
    "START=$(date -u +%Y-%m-%dT%H:%M:%SZ); END=\"\"; SCR=\"\"\n"
    "state() {   # state STATE [EXIT]: caps-job.json, written whole\n"
    "  printf '{\"title\":\"%s\",\"job\":\"%s\",\"state\":\"%s\",\"slurm_job\":\"%s\",\"array_task\":\"%s\",\"node\":\"%s\",\"scratch\":\"%s\",'"
    "'\"start\":\"%s\",\"end\":\"%s\",\"exit\":%s,\"threads\":%s,\"scheduler\":\"{scheduler}\",\"log\":\"run.log\",\"thermo\":\"thermo.csv\",\"progress\":\"progress.jsonl\",'"
    "'\"frames\":\"traj.lammpstrj\"}\\n' \"$TITLE\" \"$JOB\" \"$1\" \"$JOBKEY\" \"$TASK\" \"$(hostname)\" \"$SCR\" \"$START\" \"$END\" \"${2:-null}\" \"$CAPS_THREADS\" "
    "> \"$JOBDIR/caps-job.json.tmp\" && mv \"$JOBDIR/caps-job.json.tmp\" \"$JOBDIR/caps-job.json\"\n"
    "}\n"
    "{scratch_setup}\n"
    "state running\n"
    "for f in \"$JOBDIR\"/* ; do case \"$(basename \"$f\")\" in out|slurm-*|pbs.*) ;; *) cp -a \"$f\" \"$SCR/\" ;; esac; done\n"
    "cd \"$SCR\" || { state failed 1; exit 1; }\n"
    "copy_back() { [ \"$SCR\" = \"$JOBDIR/out\" ] && return 0; mkdir -p \"$JOBDIR/out\" && cp -a \"$SCR/.\" \"$JOBDIR/out/\"; }\n"
    "finish() {   # finish STATE EXIT: copy back, release the scratch only after a good copy\n"
    "  END=$(date -u +%Y-%m-%dT%H:%M:%SZ)\n"
    "  if copy_back; then\n"
    "    if [ \"$KEEP\" != 1 ] && [ -n \"$RELEASE\" ]; then eval \"$RELEASE\" && SCR=\"$JOBDIR/out\"; fi\n"
    "  else\n"
    "    echo \"caps-job: copying $SCR back failed: the scratch is kept\"\n"
    "  fi\n"
    "  state \"$1\" \"$2\"\n"
    "}\n"
    "on_limit() {   # USR1 (the time limit near) or TERM (the limit, or caps job cancel): a checkpoint, the results back\n"
    "  if [ -e \"$JOBDIR/cancel-requested\" ]; then echo \"caps-job: cancelled\"; else echo \"caps-job: time limit near: asking caps for a checkpoint\"; fi\n"
    "  kill -USR1 \"$PID\" 2>/dev/null; wait \"$PID\"\n"
    "  if [ -e \"$JOBDIR/cancel-requested\" ]; then finish cancelled 130; exit 130; fi\n"
    "  finish timeout 124; exit 124\n"
    "}\n"
    "trap on_limit USR1 TERM\n"
    "{command} > run.log 2>&1 &\n"
    "PID=$!\n"
    "tail -n +1 -f run.log --pid=\"$PID\" 2>/dev/null &\n"
    "wait \"$PID\"; RC=$?\n"
    "trap - USR1 TERM\n"
    "if [ \"$RC\" -eq 0 ]; then finish finished 0; else finish failed \"$RC\"; fi\n"
    "exit \"$RC\"\n";
}  // namespace

std::string default_job_template() { return std::string(kSlurmHead) + kBody; }

std::string job_template_for(const std::string& scheduler) {
  if (scheduler == "pbs") return std::string(kPbsHead) + kBody;
  if (scheduler == "direct") return std::string(kDirectHead) + kBody;
  return default_job_template();
}

std::string scratch_setup(const HostProfile& h) {
  const std::string name = "\"caps-$JOBKEY${TASK:+-$TASK}\"";
  if (h.scratch == "workspace") {
    const std::string fsopt = h.ws_filesystem.empty() ? "" : " -F " + h.ws_filesystem;
    return "# scratch: a workspace (ws_allocate prints its path)\n"
           "SCR=$(ws_allocate" + fsopt + " " + name + " " + std::to_string(std::max(1, h.ws_days)) + " 2>/dev/null | tail -n 1)\n"
           "if [ -z \"$SCR\" ] || [ ! -d \"$SCR\" ]; then echo \"caps-job: ws_allocate gave no workspace\"; state failed 1; exit 1; fi\n"
           "RELEASE=\"ws_release" + fsopt + " caps-$JOBKEY${TASK:+-$TASK}\"";
  }
  if (h.scratch == "env") {
    const std::string var = h.scratch_var.empty() ? "$TMPDIR" : h.scratch_var;
    return "# scratch: a folder under " + var + "\n"
           "SCR=\"" + var + "/caps-$JOBKEY${TASK:+-$TASK}\"\n"
           "mkdir -p \"$SCR\" || { echo \"caps-job: cannot make $SCR\"; state failed 1; exit 1; }\n"
           "RELEASE=\"rm -rf '$SCR'\"";
  }
  return "# no scratch: the run happens in out/ of the job folder\n"
         "SCR=\"$JOBDIR/out\"; mkdir -p \"$SCR\"\n"
         "RELEASE=\"\"";
}

std::string fill_job_template(const std::string& tmpl, const std::map<std::string, std::string>& values) {
  std::stringstream in(tmpl);
  std::string line, out;
  while (std::getline(in, line)) {
    bool empty_directive = false;
    std::string filled;
    for (size_t i = 0; i < line.size();) {
      if (line[i] == '{') {
        const size_t e = line.find('}', i);
        if (e != std::string::npos) {
          const std::string key = line.substr(i + 1, e - i - 1);
          if (const auto it = values.find(key); it != values.end()) {
            filled += it->second;
            if (it->second.empty()) empty_directive = true;
            i = e + 1;
            continue;
          }
        }
      }
      filled += line[i++];
    }
    const bool directive = line.rfind("#SBATCH", 0) == 0 || line.rfind("#PBS", 0) == 0;
    if (directive && empty_directive) continue;            // a directive without its value
    if (trim(filled).empty() && !trim(line).empty()) continue;   // a placeholder line that stayed empty
    out += filled + "\n";
  }
  return out;
}

namespace {
std::map<std::string, std::string> template_values(const HostProfile& h, const JobRequest& r, const std::string& title, const std::string& job,
                                                   const std::string& jobdir, int array_count) {
  const int cpus = r.cpus > 0 ? r.cpus : h.cpus;
  const std::string mem = !r.mem.empty() ? r.mem : h.mem_per_cpu;
  std::map<std::string, std::string> v;
  v["title"] = title;
  v["job"] = job;
  v["jobdir"] = jobdir;
  v["account"] = h.account;
  v["partition"] = !r.partition.empty() ? r.partition : h.partition;
  v["cpus"] = std::to_string(cpus);
  v["mem"] = mem;
  v["mem_total"] = total_mem(mem, cpus);
  v["time"] = !r.time.empty() ? r.time : h.time;
  std::string mail;
  if (!h.mail_user.empty()) {
    if (h.scheduler == "pbs") mail = "#PBS -m " + (h.mail_type.empty() ? std::string("ae") : h.mail_type) + "\n#PBS -M " + h.mail_user;
    else if (h.scheduler == "slurm") mail = "#SBATCH --mail-type=" + (h.mail_type.empty() ? std::string("END,FAIL") : h.mail_type) + "\n#SBATCH --mail-user=" + h.mail_user;
  }
  v["mail_lines"] = mail;
  v["modules"] = h.modules;
  v["scratch_setup"] = scratch_setup(h);
  v["command"] = "bash ./cmd.sh";
  v["keep_workspace"] = (r.keep_workspace >= 0 ? r.keep_workspace == 1 : h.keep_workspace) ? "1" : "0";
  v["array"] = array_count <= 0 ? "" : h.scheduler == "pbs" ? "#PBS -J 0-" + std::to_string(array_count - 1) : h.scheduler == "slurm" ? "#SBATCH --array=0-" + std::to_string(array_count - 1) : "";
  v["caps"] = "\"" + (!r.caps_bin.empty() ? expand_home(r.caps_bin) : expand_home(h.root) + "/bin/caps") + "\"";
  v["version"] = std::string(version_string()) + " " + commit_string();
  v["scheduler"] = h.scheduler;
  return v;
}

std::string script_template(const HostProfile& h) { return !h.job_template.empty() ? h.job_template : job_template_for(h.scheduler); }
}  // namespace

JobFolder make_job(const HostProfile& h, const JobRequest& r) {
  if (trim(r.command).empty()) throw std::invalid_argument("the job needs a command (caps …)");
  const std::string title = sanitize_title(r.title.empty() ? "job" : r.title);
  const std::string kind = sanitize_title(r.kind.empty() ? "run" : r.kind);
  const fs::path tdir = fs::path(expand_home(h.root)) / title;
  fs::create_directories(tdir);
  JobFolder f;
  f.job = next_job_name(tdir, kind);
  const fs::path dir = tdir / f.job;
  f.dir = fs::absolute(dir).lexically_normal().string();
  write_task(dir, title, f.job, r.command, r.inputs, "");
  f.script = (dir / "job.sh").string();
  write_file(f.script, fill_job_template(script_template(h), template_values(h, r, title, f.job, f.dir, 0)));
  fs::permissions(f.script, fs::perms::owner_exec, fs::perm_options::add);
  return f;
}

JobFolder make_array_job(const HostProfile& h, const JobRequest& r, const std::vector<ArrayTask>& tasks) {
  if (tasks.empty()) throw std::invalid_argument("the array list has no tasks");
  const std::string kind = sanitize_title(r.kind.empty() ? "run" : r.kind);
  const fs::path root = expand_home(h.root);
  const fs::path adir_parent = root / "_arrays";
  fs::create_directories(adir_parent);
  JobFolder f;
  f.job = next_job_name(adir_parent, kind);
  const fs::path adir = adir_parent / f.job;
  fs::create_directories(adir);
  f.dir = fs::absolute(adir).lexically_normal().string();
  std::string list;
  for (size_t k = 0; k < tasks.size(); ++k) {
    const std::string title = sanitize_title(tasks[k].title);
    const fs::path tdir = root / title;
    fs::create_directories(tdir);
    const std::string job = next_job_name(tdir, kind);
    const fs::path dir = fs::absolute(tdir / job).lexically_normal();
    write_task(dir, title, job, tasks[k].command, r.inputs, "_arrays/" + f.job);
    f.task_dirs.push_back(dir.string());
    list += dir.string() + "\n";
  }
  write_file((adir / "tasks.txt").string(), list);
  f.script = (adir / "job.sh").string();
  write_file(f.script, fill_job_template(script_template(h), template_values(h, r, "array", f.job, f.dir, int(tasks.size()))));
  fs::permissions(f.script, fs::perms::owner_exec, fs::perm_options::add);
  Json st = Json::object();
  st["job"] = f.job;
  st["state"] = std::string("created");
  st["tasks"] = double(tasks.size());
  write_job_state(f.dir, st);
  return f;
}

JobFolder make_resume_job(const HostProfile& h, const std::string& stopped_dir, const JobRequest& overrides) {
  const Json st = read_job_state(stopped_dir);
  const fs::path out = fs::path(stopped_dir) / "out";
  fs::path from = out;
  const std::string scr = str_or(st, "scratch", "");
  std::error_code ec;
  if (!fs::exists(out / "resume.txt") && !scr.empty() && fs::exists(fs::path(scr) / "resume.txt", ec)) from = scr;   // the workspace was kept
  if (!fs::exists(from / "resume.txt"))
    throw std::invalid_argument("the job left no resume.txt: it finished, or it stopped before its first checkpoint (submit it again)");
  std::ifstream rf(from / "resume.txt");
  std::string line;
  std::getline(rf, line);
  JobRequest r = overrides;
  r.command = trim(line);
  if (r.title.empty()) r.title = str_or(st, "title", fs::path(stopped_dir).parent_path().filename().string());
  if (r.kind.empty() || r.kind == "run") {
    const std::string job = str_or(st, "job", fs::path(stopped_dir).filename().string());
    r.kind = job.substr(0, job.rfind('-'));
  }
  static const std::set<std::string> skip = {"caps-job.json", "caps-job.json.tmp", "cmd.txt", "cmd.sh", "job.sh", "run.log", "resume.txt"};
  for (const auto& e : fs::directory_iterator(from, ec))
    if (e.is_regular_file(ec) && !skip.count(e.path().filename().string()) && e.path().filename().string().rfind("slurm-", 0) != 0)
      r.inputs.push_back(e.path().string());
  JobFolder f = make_job(h, r);
  Json c = Json::object();
  c["resumed_by"] = f.dir;
  update_job_state(stopped_dir, c);
  Json n = Json::object();
  n["resumes"] = std::string(fs::absolute(stopped_dir).lexically_normal().string());
  update_job_state(f.dir, n);
  return f;
}

std::vector<ArrayTask> read_array_list(const std::string& text) {
  std::vector<ArrayTask> out;
  std::stringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') continue;
    const size_t sp = line.find_first_of(" \t");
    if (sp == std::string::npos) throw std::invalid_argument("array list: \"" + line + "\" needs a title and a command");
    out.push_back({line.substr(0, sp), trim(line.substr(sp))});
  }
  return out;
}

// ------------------------------------------------------------------------------------------------- state file
Json read_job_state(const std::string& dir) {
  const fs::path p = fs::path(dir) / "caps-job.json";
  if (!fs::exists(p)) return Json::object();
  try {
    return Json::parse(slurp(p.string()));
  } catch (const std::exception&) {
    return Json::object();   // being rewritten: the next look reads it
  }
}

void write_job_state(const std::string& dir, const Json& j) {
  const fs::path p = fs::path(dir) / "caps-job.json";
  const fs::path tmp = fs::path(dir) / "caps-job.json.tmp";
  write_file(tmp.string(), j.dump(1) + "\n");
  fs::rename(tmp, p);
}

void update_job_state(const std::string& dir, const Json& changes) {
  Json j = read_job_state(dir);
  if (!j.is_object()) j = Json::object();
  if (changes.is_object())
    for (const auto& [k, v] : changes.members()) j[k] = v;
  write_job_state(dir, j);
}

// ------------------------------------------------------------------------------------------------- SLURM output
QueueInfo parse_squeue_line(const std::string& line) {
  QueueInfo q;
  const std::string t = trim(line);
  if (t.empty()) return q;
  std::vector<std::string> f;
  std::stringstream ss(t);
  std::string x;
  while (std::getline(ss, x, '|')) f.push_back(trim(x));
  if (f.empty() || f[0].empty()) return q;
  q.found = true;
  q.state = f[0];
  if (f.size() > 1) q.reason = f[1] == "None" ? "" : f[1];
  if (f.size() > 2) q.elapsed = f[2];
  if (f.size() > 3) q.limit = f[3];
  if (f.size() > 4) q.node = f[4];
  return q;
}

AcctInfo parse_sacct(const std::string& text) {
  AcctInfo a;
  std::stringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    std::vector<std::string> f;
    std::stringstream ss(line);
    std::string x;
    while (std::getline(ss, x, '|')) f.push_back(trim(x));
    if (f.size() < 3 || f[0].empty()) continue;
    if (!a.found) {   // the job's own line comes first; its steps (batch, extern) follow and carry the memory
      a.found = true;
      a.state = f[0].substr(0, f[0].find(' '));   // "CANCELLED by 123" → CANCELLED
      a.exit_code = f[1];
      a.elapsed = f[2];
    }
    if (f.size() > 3 && !f[3].empty()) a.max_rss = f[3];
  }
  return a;
}

std::string job_state_from_slurm(const std::string& s) {
  if (s == "PENDING" || s == "CONFIGURING" || s == "REQUEUED") return "pending";
  if (s == "RUNNING" || s == "COMPLETING" || s == "SUSPENDED") return "running";
  if (s == "COMPLETED") return "finished";
  if (s == "TIMEOUT" || s == "DEADLINE") return "timeout";
  if (s == "CANCELLED") return "cancelled";
  if (s == "OUT_OF_MEMORY") return "failed";
  if (s.empty()) return "";
  return "failed";   // FAILED, NODE_FAIL, BOOT_FAIL, PREEMPTED
}

int64_t last_frame_offset(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return -1;
  f.seekg(0, std::ios::end);
  const int64_t size = f.tellg();
  const std::string key = "ITEM: TIMESTEP";
  const int64_t chunk = 1 << 20;
  int64_t end = size;
  while (end > 0) {
    const int64_t start = std::max<int64_t>(0, end - chunk);
    const int64_t len = std::min<int64_t>(size, end + int64_t(key.size())) - start;   // overlap: a key across the cut
    std::string buf(size_t(len), '\0');
    f.clear();
    f.seekg(start);
    f.read(&buf[0], len);
    if (const size_t at = buf.rfind(key); at != std::string::npos) return start + int64_t(at);
    end = start;
  }
  return -1;
}

double eta_seconds(double fraction, double elapsed_s) {
  if (!(fraction > 0) || fraction >= 1 || !(elapsed_s > 0)) return -1;
  return elapsed_s * (1 - fraction) / fraction;
}

double slurm_duration_seconds(const std::string& s) {
  const std::string t = trim(s);
  if (t.empty() || t == "UNLIMITED" || t == "INVALID" || t == "NOT_SET") return -1;
  double days = 0;
  std::string rest = t;
  if (const size_t d = t.find('-'); d != std::string::npos) days = std::atof(t.substr(0, d).c_str()), rest = t.substr(d + 1);
  std::vector<double> parts;
  std::stringstream ss(rest);
  std::string x;
  while (std::getline(ss, x, ':')) parts.push_back(std::atof(x.c_str()));
  double secs = 0;
  if (parts.size() == 3) secs = parts[0] * 3600 + parts[1] * 60 + parts[2];
  else if (parts.size() == 2) secs = parts[0] * 60 + parts[1];
  else if (parts.size() == 1) secs = parts[0] * 60;   // minutes (sbatch's form)
  else return -1;
  return days * 86400 + secs;
}

// ------------------------------------------------------------------------------------------------- scaling
std::vector<int> scaling_threads(int cores) {
  std::set<int> t;
  const int top = cores > 0 ? cores : 32;
  for (int p = 1; p < top; p *= 2) t.insert(p);
  if (cores > 0) t.insert(std::max(1, cores / 2)), t.insert(cores);
  else t.insert(32);
  return {t.begin(), t.end()};
}

Json scaling_table(const std::vector<Json>& bench) {
  std::vector<std::pair<double, double>> rows;   // threads, ns/day
  std::string node, commit;
  for (const auto& b : bench) {
    if (!b.is_object() || !b.has("threads") || !b.has("ns_per_day")) continue;
    rows.push_back({b["threads"].number(), b["ns_per_day"].number()});
    if (node.empty()) node = str_or(b, "node", "");
    if (commit.empty()) commit = str_or(b, "commit", "");
  }
  std::sort(rows.begin(), rows.end());
  Json t = Json::object(), arr = Json::array();
  t["node"] = node, t["commit"] = commit;
  if (!rows.empty()) {
    const double t0 = rows[0].first, r0 = rows[0].second;
    double best = 0;
    int best_threads = 0;
    for (const auto& [th, nd] : rows) {
      Json r = Json::object();
      r["threads"] = th, r["ns_per_day"] = nd;
      const double speedup = r0 > 0 ? nd / r0 : 0, eff = th > 0 && t0 > 0 ? speedup / (th / t0) : 0;
      r["speedup"] = speedup, r["efficiency"] = eff;
      arr.push_back(r);
      if (eff >= 0.7 && nd > best) best = nd, best_threads = int(th);   // the fastest that still uses its cores well
    }
    if (best_threads > 0) t["suggested_threads"] = double(best_threads);
  }
  t["rows"] = arr;
  return t;
}

std::string scaling_text(const Json& t) {
  std::string s = "threads   ns/day   speedup  efficiency\n";
  char b[128];
  for (const auto& r : t["rows"].items()) {
    std::snprintf(b, sizeof b, "%7.0f %8.3f %9.2f %10.0f %%\n", r["threads"].number(), r["ns_per_day"].number(), r["speedup"].number(), 100 * r["efficiency"].number());
    s += b;
  }
  if (t.has("suggested_threads"))
    s += "suggested: " + std::to_string(int(t["suggested_threads"].number())) + " threads (the fastest at 70 % parallel efficiency or more)\n";
  return s;
}

// ------------------------------------------------------------------------------------------------- poll
std::string job_live_dir(const std::string& dir, const Json& st) {
  const std::string state = str_or(st, "state", "");
  const std::string scr = str_or(st, "scratch", "");
  std::error_code ec;
  if (!scr.empty() && fs::is_directory(scr, ec) && (state == "running" || state == "timeout" || state == "failed")) return scr;
  if (fs::is_directory(fs::path(dir) / "out", ec)) return (fs::path(dir) / "out").string();
  return dir;
}

namespace {
// bytes of a file from an offset (a file that shrank — rewritten — is read from the start); next: the new offset
std::string read_from(const std::string& path, int64_t offset, size_t max_bytes, int64_t* next, int64_t* size_out) {
  std::ifstream f(path, std::ios::binary);
  if (!f) { *next = offset; *size_out = -1; return ""; }
  f.seekg(0, std::ios::end);
  const int64_t size = f.tellg();
  *size_out = size;
  if (offset > size || offset < 0) offset = 0;
  int64_t start = offset;
  if (size - start > int64_t(max_bytes)) start = size - int64_t(max_bytes);   // far behind: the newest part
  std::string buf(size_t(size - start), '\0');
  f.seekg(start);
  if (!buf.empty()) f.read(&buf[0], int64_t(buf.size()));
  *next = size;
  return buf;
}
}  // namespace

Json job_poll(const PollRequest& p, const std::string& squeue_line, const std::string& sacct_text) {
  Json out = Json::object();
  const Json st = read_job_state(p.dir);
  out["dir"] = p.dir;
  out["job_state"] = st;
  const QueueInfo q = parse_squeue_line(squeue_line);
  Json qj = Json::object();
  qj["found"] = q.found;
  if (q.found) {
    qj["state"] = q.state, qj["reason"] = q.reason, qj["elapsed"] = q.elapsed, qj["limit"] = q.limit, qj["node"] = q.node;
    qj["elapsed_s"] = slurm_duration_seconds(q.elapsed), qj["limit_s"] = slurm_duration_seconds(q.limit);
  }
  out["queue"] = qj;
  const AcctInfo a = parse_sacct(sacct_text);
  if (a.found) {
    Json aj = Json::object();
    aj["state"] = a.state, aj["exit_code"] = a.exit_code, aj["elapsed"] = a.elapsed, aj["max_rss"] = a.max_rss;
    out["acct"] = aj;
  }
  // the state: the queue's while it is there, else the script's, else the accounting's
  std::string state = q.found ? job_state_from_slurm(q.state) : str_or(st, "state", "");
  if (!q.found && a.found && (state.empty() || state == "submitted" || state == "running" || state == "pending")) state = job_state_from_slurm(a.state);
  out["state"] = state;
  const std::string live = job_live_dir(p.dir, st);
  out["live_dir"] = live;
  int64_t next = 0, size = 0;
  out["log"] = read_from((fs::path(live) / str_or(st, "log", "run.log")).string(), p.log_offset, p.max_bytes, &next, &size);
  out["log_offset"] = double(next);
  std::string prog = read_from((fs::path(live) / str_or(st, "progress", "progress.jsonl")).string(), p.progress_offset, p.max_bytes, &next, &size);
  size_t used = 0;
  // a chunk cut at the front (far behind) starts inside a line: from the first newline
  if (p.progress_offset > 0 && size >= 0 && next - int64_t(prog.size()) > p.progress_offset) {
    const size_t nl = prog.find('\n');
    prog = nl == std::string::npos ? "" : prog.substr(nl + 1);
  }
  const auto lines = read_progress_lines(prog, &used);
  Json pl = Json::array();
  for (const auto& l : lines) pl.push_back(l);
  out["progress"] = pl;
  out["progress_offset"] = double(next - int64_t(prog.size() - used));
  // the scheduler's error file (slurm-<id>.err, pbs.err)
  std::string err_path;
  const std::string sj = str_or(st, "slurm_job", "");
  if (!sj.empty() && fs::exists(fs::path(p.dir) / ("slurm-" + sj + ".err"))) err_path = (fs::path(p.dir) / ("slurm-" + sj + ".err")).string();
  else if (fs::exists(fs::path(p.dir) / "pbs.err")) err_path = (fs::path(p.dir) / "pbs.err").string();
  if (!err_path.empty()) {
    out["err"] = read_from(err_path, p.err_offset, p.max_bytes, &next, &size);
    out["err_offset"] = double(next);
  }
  // the live folder's files
  Json files = Json::array();
  std::error_code ec;
  if (fs::is_directory(live, ec))
    for (const auto& e : fs::directory_iterator(live, ec)) {
      if (!e.is_regular_file(ec)) continue;
      Json fj = Json::object();
      fj["name"] = e.path().filename().string();
      fj["bytes"] = double(e.file_size(ec));
      files.push_back(fj);
    }
  out["files"] = files;
  const std::string frames = (fs::path(live) / str_or(st, "frames", "traj.lammpstrj")).string();
  out["last_frame_offset"] = double(last_frame_offset(frames));
  // the latest progress line's fraction and the ETA
  for (auto it = pl.items().rbegin(); it != pl.items().rend(); ++it)
    if (it->has("fraction") && (*it)["fraction"].is_number()) {
      out["fraction"] = (*it)["fraction"].number();
      if (q.found) out["eta_s"] = eta_seconds((*it)["fraction"].number(), slurm_duration_seconds(q.elapsed));
      break;
    }
  return out;
}

}  // namespace caps
