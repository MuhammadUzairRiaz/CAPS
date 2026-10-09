// CAPS jobs on a cluster or workstation (caps/cluster_job.hpp): the host profile round trip (old files keep defaults),
// the job script for each scheduler and scratch kind (valid bash, empty directives left out), folders and arrays, the
// state file, SLURM's squeue / sacct lines, durations and the ETA, the last dump frame's offset, and the poll a viewer
// reads (offsets, a progress line cut off mid-write, the live folder while it runs).
#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/cluster_job.hpp"

using namespace caps;
namespace fs = std::filesystem;

namespace {
std::string slurp(const fs::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}
void spit(const fs::path& p, const std::string& s) {
  std::ofstream f(p, std::ios::binary);
  f << s;
}
fs::path tmpdir(const char* name) {
  auto d = fs::temp_directory_path() / name;
  fs::remove_all(d);
  fs::create_directories(d);
  return d;
}
bool bash_ok(const fs::path& script) {
  if (std::system("command -v bash > /dev/null 2>&1") != 0) return true;
  return std::system(("bash -n '" + script.string() + "'").c_str()) == 0;
}
}  // namespace

TEST(ClusterJob, ProfileRoundTripAndOldFiles) {
  HostProfile h;
  h.name = "my cluster", h.ssh = "me@login.example.org", h.account = "proj1", h.scheduler = "pbs", h.scratch = "env";
  h.scratch_var = "$TMPDIR", h.cpus = 48, h.mem_per_cpu = "3G", h.keep_workspace = true, h.node_cores = 96;
  const HostProfile b = host_profile_from_json(Json::parse(host_profile_json(h).dump(0)));
  EXPECT_EQ(b.ssh, h.ssh);
  EXPECT_EQ(b.scheduler, "pbs");
  EXPECT_EQ(b.cpus, 48);
  EXPECT_TRUE(b.keep_workspace);
  EXPECT_EQ(b.node_cores, 96);
  // a file from before a key existed, and nonsense values: the defaults
  const HostProfile o = host_profile_from_json(Json::parse(R"({"account":"x","scheduler":"lsf","scratch":"tape","cpus":0})"));
  EXPECT_EQ(o.account, "x");
  EXPECT_EQ(o.scheduler, "slurm");
  EXPECT_EQ(o.scratch, "none");
  EXPECT_EQ(o.cpus, 1);
  EXPECT_EQ(o.root, "~/CAPS");
  // presets carry no account, user or address
  for (const auto& [key, p] : host_presets()) {
    EXPECT_TRUE(p.account.empty()) << key;
    EXPECT_TRUE(p.ssh.empty()) << key;
    EXPECT_TRUE(p.mail_user.empty()) << key;
  }
}

TEST(ClusterJob, TemplateDropsEmptyDirectives) {
  const std::string t = fill_job_template("#!/bin/bash\n#SBATCH --account={account}\n#SBATCH --time={time}\n{mail_lines}\necho {title}\n{unknown}\n",
                                          {{"account", ""}, {"time", "01:00:00"}, {"mail_lines", ""}, {"title", "PS"}});
  EXPECT_EQ(t.find("--account"), std::string::npos);
  EXPECT_NE(t.find("#SBATCH --time=01:00:00"), std::string::npos);
  EXPECT_NE(t.find("echo PS"), std::string::npos);
  EXPECT_NE(t.find("{unknown}"), std::string::npos);   // not ours: left as written
  EXPECT_EQ(t.find("\n\n"), std::string::npos);        // the empty mail line is gone
}

TEST(ClusterJob, ScriptsForEachSchedulerAndScratch) {
  const auto root = tmpdir("caps_cj_scripts");
  for (const std::string sched : {"slurm", "pbs", "direct"})
    for (const std::string scr : {"workspace", "env", "none"}) {
      HostProfile h;
      h.root = (root / (sched + "_" + scr)).string();
      h.scheduler = sched, h.scratch = scr, h.account = "acct", h.ws_filesystem = scr == "workspace" ? "fsA" : "";
      h.mail_user = "someone@example.org";
      JobRequest r;
      r.title = "PS melt (atactic)", r.kind = "md", r.command = "caps md structure.data -o out.data --steps 1000";
      const JobFolder f = make_job(h, r);
      const std::string s = slurp(f.script);
      EXPECT_TRUE(bash_ok(f.script)) << sched << " " << scr << "\n" << s;
      EXPECT_NE(f.dir.find("PS_melt__atactic_"), std::string::npos) << f.dir;
      EXPECT_EQ(fs::path(f.dir).filename().string(), "md-1");
      EXPECT_NE(slurp(fs::path(f.dir) / "cmd.sh").find("exec \"$CAPS\" md structure.data"), std::string::npos);
      EXPECT_EQ(read_job_state(f.dir)["state"].str(), "created");
      if (sched == "slurm") {
        EXPECT_NE(s.find("#SBATCH --cpus-per-task=16"), std::string::npos);
        EXPECT_NE(s.find("#SBATCH --account=acct"), std::string::npos);
        EXPECT_NE(s.find("#SBATCH --signal=B:USR1@600"), std::string::npos);
        EXPECT_NE(s.find("#SBATCH --mail-user=someone@example.org"), std::string::npos);
      }
      if (sched == "pbs") {
        EXPECT_NE(s.find("#PBS -l select=1:ncpus=16:mem=32gb"), std::string::npos) << s;
        EXPECT_EQ(s.find("#SBATCH"), std::string::npos);
      }
      if (sched == "direct") EXPECT_EQ(s.find("#SBATCH"), std::string::npos);
      if (scr == "workspace") {
        EXPECT_NE(s.find("ws_allocate -F fsA"), std::string::npos);
        EXPECT_NE(s.find("ws_release -F fsA"), std::string::npos);
      }
      if (scr == "none") EXPECT_EQ(s.find("ws_allocate"), std::string::npos);
      // the scratch is released only after the copy back succeeded
      const size_t copy = s.find("if copy_back; then"), rel = s.find("eval \"$RELEASE\"");
      EXPECT_LT(copy, rel);
      EXPECT_NE(s.find("export CAPS_THREADS="), std::string::npos);
      EXPECT_NE(s.find("export OMP_NUM_THREADS=1"), std::string::npos);
      // a second job of the same kind: the next number
      EXPECT_EQ(fs::path(make_job(h, r).dir).filename().string(), "md-2");
    }
}

TEST(ClusterJob, ArrayJobs) {
  const auto root = tmpdir("caps_cj_array");
  spit(root / "in.data", "x");
  HostProfile h;
  h.root = (root / "CAPS").string();
  JobRequest r;
  r.kind = "md";
  r.inputs = {(root / "in.data").string()};
  const auto tasks = read_array_list("# seeds\nPS_s1 caps md in.data -o o.data --seed 1\n\nPS_s2\tcaps md in.data -o o.data --seed 2\n");
  ASSERT_EQ(tasks.size(), 2u);
  EXPECT_EQ(tasks[1].title, "PS_s2");
  const JobFolder f = make_array_job(h, r, tasks);
  ASSERT_EQ(f.task_dirs.size(), 2u);
  const std::string s = slurp(f.script);
  EXPECT_NE(s.find("#SBATCH --array=0-1"), std::string::npos);
  EXPECT_TRUE(bash_ok(f.script));
  EXPECT_EQ(slurp(fs::path(f.dir) / "tasks.txt"), f.task_dirs[0] + "\n" + f.task_dirs[1] + "\n");
  EXPECT_TRUE(fs::exists(fs::path(f.task_dirs[1]) / "in.data"));
  EXPECT_NE(slurp(fs::path(f.task_dirs[1]) / "cmd.txt").find("--seed 2"), std::string::npos);
  EXPECT_EQ(read_job_state(f.task_dirs[0])["array"].str(), "_arrays/md-1");
  EXPECT_THROW(read_array_list("lonely\n"), std::invalid_argument);
}

TEST(ClusterJob, StateFile) {
  const auto d = tmpdir("caps_cj_state");
  EXPECT_TRUE(read_job_state(d.string()).is_object());
  Json j = Json::object();
  j["state"] = std::string("running");
  j["scratch"] = std::string("/scratch/x");
  write_job_state(d.string(), j);
  Json c = Json::object();
  c["state"] = std::string("finished");
  c["exit"] = 0.0;
  update_job_state(d.string(), c);
  const Json r = read_job_state(d.string());
  EXPECT_EQ(r["state"].str(), "finished");
  EXPECT_EQ(r["scratch"].str(), "/scratch/x");
  spit(d / "caps-job.json", "{\"state\":\"runn");   // being rewritten: an empty object, never a throw
  EXPECT_TRUE(read_job_state(d.string()).is_object());
}

TEST(ClusterJob, SlurmLines) {
  const QueueInfo q = parse_squeue_line("RUNNING|None|1-02:03:04|2-00:00:00|n1234\n");
  ASSERT_TRUE(q.found);
  EXPECT_EQ(q.state, "RUNNING");
  EXPECT_EQ(q.reason, "");
  EXPECT_EQ(q.node, "n1234");
  EXPECT_FALSE(parse_squeue_line("").found);
  EXPECT_EQ(parse_squeue_line("PENDING|Priority|0:00|1:00:00|").reason, "Priority");
  const AcctInfo a = parse_sacct("TIMEOUT|0:0|01:00:12|\nCANCELLED|0:15|01:00:13|1234M\nCOMPLETED|0:0|01:00:12|88M\n");
  ASSERT_TRUE(a.found);
  EXPECT_EQ(a.state, "TIMEOUT");
  EXPECT_EQ(a.max_rss, "88M");
  EXPECT_EQ(parse_sacct("CANCELLED by 4242|0:0|00:01:00|\n").state, "CANCELLED");
  EXPECT_EQ(job_state_from_slurm("PENDING"), "pending");
  EXPECT_EQ(job_state_from_slurm("COMPLETED"), "finished");
  EXPECT_EQ(job_state_from_slurm("TIMEOUT"), "timeout");
  EXPECT_EQ(job_state_from_slurm("OUT_OF_MEMORY"), "failed");
  EXPECT_DOUBLE_EQ(slurm_duration_seconds("1-02:03:04"), 93784);
  EXPECT_DOUBLE_EQ(slurm_duration_seconds("03:04"), 184);
  EXPECT_DOUBLE_EQ(slurm_duration_seconds("UNLIMITED"), -1);
  EXPECT_DOUBLE_EQ(eta_seconds(0.25, 100), 300);
  EXPECT_DOUBLE_EQ(eta_seconds(0, 100), -1);
}

TEST(ClusterJob, LastFrameOffset) {
  const auto d = tmpdir("caps_cj_frames");
  std::string dump;
  for (int k = 0; k < 3; ++k) dump += "ITEM: TIMESTEP\n" + std::to_string(100 * k) + "\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\nITEM: ATOMS id mol type xu yu zu\n1 1 1 0 0 0\n";
  spit(d / "t.lammpstrj", dump);
  const int64_t off = last_frame_offset((d / "t.lammpstrj").string());
  EXPECT_EQ(dump.substr(size_t(off), 18), "ITEM: TIMESTEP\n200");
  // a frame across the 1 MiB read window
  std::string big(1 << 20, 'x');
  spit(d / "big.lammpstrj", big + dump);
  EXPECT_EQ(last_frame_offset((d / "big.lammpstrj").string()), int64_t(big.size()) + off);
  EXPECT_EQ(last_frame_offset((d / "none.lammpstrj").string()), -1);
}

TEST(ClusterJob, PollReadsWhatIsNew) {
  const auto d = tmpdir("caps_cj_poll");
  const auto scratch = d / "scratch";
  fs::create_directories(scratch);
  Json st = Json::object();
  st["title"] = std::string("PS"), st["job"] = std::string("md-1"), st["state"] = std::string("running");
  st["scratch"] = scratch.string(), st["slurm_job"] = std::string("77");
  write_job_state(d.string(), st);
  spit(scratch / "run.log", "step 100\nstep 200\n");
  spit(scratch / "progress.jsonl", "{\"state\":\"running\"}\n{\"step\":100,\"fraction\":0.25}\n{\"step\":2");   // the last line half written
  spit(d / "slurm-77.err", "warning: x\n");
  PollRequest p;
  p.dir = d.string();
  Json r = job_poll(p, "RUNNING|None|00:10:00|01:00:00|n1", "");
  EXPECT_EQ(r["state"].str(), "running");
  EXPECT_EQ(r["live_dir"].str(), scratch.string());
  EXPECT_EQ(r["log"].str(), "step 100\nstep 200\n");
  ASSERT_EQ(r["progress"].items().size(), 2u);
  EXPECT_DOUBLE_EQ(r["fraction"].number(), 0.25);
  EXPECT_DOUBLE_EQ(r["eta_s"].number(), 1800);
  EXPECT_EQ(r["err"].str(), "warning: x\n");
  // the next look: only the new log bytes, and the progress line finished since
  std::ofstream(scratch / "run.log", std::ios::app) << "step 300\n";
  std::ofstream(scratch / "progress.jsonl", std::ios::app) << "00,\"fraction\":0.5}\n";
  p.log_offset = int64_t(r["log_offset"].number());
  p.progress_offset = int64_t(r["progress_offset"].number());
  r = job_poll(p, "RUNNING|None|00:20:00|01:00:00|n1", "");
  EXPECT_EQ(r["log"].str(), "step 300\n");
  ASSERT_EQ(r["progress"].items().size(), 1u);
  EXPECT_DOUBLE_EQ(r["progress"].items()[0]["step"].number(), 200);
  // the job ended: off the queue, the accounting says why; out/ is read
  fs::create_directories(d / "out");
  spit(d / "out" / "run.log", "done\n");
  st["state"] = std::string("finished");
  write_job_state(d.string(), st);
  p.log_offset = 0;
  r = job_poll(p, "", "COMPLETED|0:0|00:30:00|\n");
  EXPECT_EQ(r["state"].str(), "finished");
  EXPECT_EQ(r["live_dir"].str(), (d / "out").string());
  EXPECT_EQ(r["log"].str(), "done\n");
  // the scheduler's word when the script never ran (killed in the queue)
  st["state"] = std::string("submitted");
  write_job_state(d.string(), st);
  EXPECT_EQ(job_poll(p, "", "CANCELLED by 1|0:0|00:00:00|\n")["state"].str(), "cancelled");
}

TEST(ClusterJob, ChainedCommandsAndFolderInputs) {
  const auto root = tmpdir("caps_cj_chain");
  fs::create_directories(root / "bonded");
  spit(root / "bonded" / "bonds.table", "x");
  HostProfile h;
  h.root = (root / "CAPS").string();
  JobRequest r;
  r.title = "T", r.kind = "mech";
  r.command = "caps tensile s.data -o a.data && caps elastic s.data --json e.json; caps info a.data";
  r.inputs = {(root / "bonded").string()};
  const JobFolder f = make_job(h, r);
  const std::string sh = slurp(fs::path(f.dir) / "cmd.sh");
  EXPECT_EQ(sh.find("exec "), std::string::npos);   // several commands: none replaces the shell
  EXPECT_NE(sh.find("\"$CAPS\" tensile s.data -o a.data && \"$CAPS\" elastic"), std::string::npos) << sh;
  EXPECT_NE(sh.find("; \"$CAPS\" info a.data"), std::string::npos) << sh;
  EXPECT_TRUE(fs::exists(fs::path(f.dir) / "bonded" / "bonds.table"));   // a folder input, copied whole
  EXPECT_TRUE(bash_ok(fs::path(f.dir) / "cmd.sh"));
}

TEST(ClusterJob, TitlesAreSafeFolderNames) {
  EXPECT_EQ(sanitize_title("PBS DP-25/40"), "PBS_DP-25_40");
  EXPECT_EQ(sanitize_title("../../etc"), "_.._etc");
  EXPECT_EQ(sanitize_title("  "), "job");
}
