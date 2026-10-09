// CAPS jobs on a SLURM cluster: CAPS's own engines run as batch jobs, from the command line on the cluster (caps job …)
// or from the Studio (which runs the same commands over ssh), with one folder layout, one job script and one state file.
//
//   host profile   ~/CAPS/host.json: account, partition, the modules for jobs and for building, the root folder, the
//                  scratch kind (a workspace from ws_allocate, a variable such as $SCRATCH, or none), default cores,
//                  memory, time limit, optional mail
//   folder         <root>/<title>/<kind>-<n>/: the inputs, cmd.txt (the caps command), job.sh, caps-job.json, the SLURM
//                  output files, out/ (the scratch folder copied back); an array job keeps its job.sh and tasks.txt in
//                  <root>/_arrays/<kind>-<n>/, one task folder per line
//   job script     from a template with {placeholders}: one node, one task, --cpus-per-task cores (CAPS uses threads,
//                  not MPI ranks), CAPS_THREADS from the job, the scratch made and recorded in caps-job.json, the inputs
//                  copied there, the command run in the background with its output in run.log, USR1 ten minutes before
//                  the time limit asking caps for a checkpoint, the scratch copied back to out/ and released only when
//                  the copy succeeded
//   caps-job.json  {title, job, state (created | submitted | running | finished | failed | timeout | cancelled),
//                  slurm_job, array_task, node, scratch, start, end, exit, log, thermo, progress, frames}
//   poll           one JSON object with everything a viewer needs since its last look (job_poll): the queue line, the
//                  state file, new bytes of run.log and progress.jsonl from offsets, the scratch listing, the byte
//                  offset of the last dump frame
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "caps/json.hpp"

namespace caps {

// Any user's own cluster or workstation: nothing here is tied to one site. Presets (host_presets) fill the fields for
// common set-ups; the user edits them.
struct HostProfile {
  std::string name;                 // the Studio's label for the host
  std::string ssh;                  // user@login-node (the Studio's connection; not used on the cluster itself)
  int port = 22;
  std::string scheduler = "slurm";  // slurm | pbs (PBS Pro / Torque) | direct (a workstation: no queue, run in the background)
  std::string account, partition;
  std::string modules;              // the job's module lines ("module purge; module load GCC")
  std::string build_modules;        // for building caps
  std::string root = "~/CAPS";      // job folders and bin/caps
  std::string scratch = "none";     // workspace (ws_allocate) | env (a variable such as $SCRATCH or $TMPDIR) | none
  std::string ws_filesystem;        // ws_allocate -F (empty: the site's default)
  int ws_days = 10;
  std::string scratch_var = "$SCRATCH";   // scratch = env: the variable (or a path)
  int cpus = 16;
  std::string mem_per_cpu = "2G";
  std::string time = "24:00:00";
  std::string mail_type, mail_user;  // optional (no default address)
  bool keep_workspace = false;
  int node_cores = 0;               // detected by the Studio's host test (sinfo); 0 unknown
  std::string node_mem;
  std::string job_template;          // empty: the built-in template
};
Json host_profile_json(const HostProfile& h);
// Starting points the user picks and edits: "slurm" (a generic SLURM cluster), "slurm-workspace" (SLURM with
// ws_allocate workspaces), "pbs", "workstation"
std::vector<std::pair<std::string, HostProfile>> host_presets();
HostProfile host_profile_from_json(const Json& j);   // missing keys keep their defaults (old settings load)
HostProfile load_host_profile(const std::string& path);   // a missing file: the defaults
void save_host_profile(const HostProfile& h, const std::string& path);
std::string default_host_profile_path();             // $HOME/CAPS/host.json
std::string expand_home(const std::string& path);    // a leading ~ → $HOME

// The built-in job script; placeholders {title} {job} {jobdir} {account} {partition} {cpus} {mem} {time} {mail_lines}
// {modules} {scratch_setup} {command} {keep_workspace} {array} {caps} {version}. A #SBATCH line whose placeholder is
// empty is left out.
std::string default_job_template();               // SLURM
std::string job_template_for(const std::string& scheduler);   // slurm | pbs | direct
std::string fill_job_template(const std::string& tmpl, const std::map<std::string, std::string>& values);
// The scratch lines for a profile (SCR, RELEASE)
std::string scratch_setup(const HostProfile& h);

std::string sanitize_title(const std::string& title);   // [A-Za-z0-9._-], others '_'; never empty

struct JobRequest {
  std::string title, kind = "run";
  std::string command;              // "caps md structure.data -o out.data …" (a leading caps becomes the full path)
  std::vector<std::string> inputs;  // files copied into the job folder
  int cpus = 0;                     // 0: the profile's
  std::string mem, time, partition; // empty: the profile's
  int keep_workspace = -1;          // -1: the profile's
  std::string caps_bin;             // empty: <root>/bin/caps
};
struct ArrayTask { std::string title, command; };

struct JobFolder {
  std::string dir;                  // the job folder (an array: <root>/_arrays/<kind>-<n>)
  std::string job;                  // <kind>-<n>
  std::string script;               // dir/job.sh
  std::vector<std::string> task_dirs;   // an array's task folders, in task order
};
// Makes the folder(s), copies the inputs, writes cmd.txt, job.sh and caps-job.json (state created)
JobFolder make_job(const HostProfile& h, const JobRequest& r);
JobFolder make_array_job(const HostProfile& h, const JobRequest& r, const std::vector<ArrayTask>& tasks);
// A job that stopped near its time limit continued as a new job of the same title and kind: its results (out/) as the
// inputs, its resume.txt as the command; the stopped job's state names the new one (resumed_by)
JobFolder make_resume_job(const HostProfile& h, const std::string& stopped_dir, const JobRequest& overrides);
std::vector<ArrayTask> read_array_list(const std::string& text);   // "title<TAB or spaces>caps …" per line, # comments

// caps-job.json
Json read_job_state(const std::string& dir);   // {} when missing
void write_job_state(const std::string& dir, const Json& j);   // whole (a temporary file renamed)
void update_job_state(const std::string& dir, const Json& changes);

// SLURM output
struct QueueInfo { bool found = false; std::string state, reason, elapsed, limit, node; };
QueueInfo parse_squeue_line(const std::string& line);   // "%T|%r|%M|%l|%N"
struct AcctInfo { bool found = false; std::string state, exit_code, elapsed, max_rss; };
AcctInfo parse_sacct(const std::string& text);          // "State|ExitCode|Elapsed|MaxRSS" lines: the job's own line
std::string job_state_from_slurm(const std::string& slurm_state);   // RUNNING → running, COMPLETED → finished …

// The byte offset of the last "ITEM: TIMESTEP" in a dump (-1: none): the last frame is read from there
int64_t last_frame_offset(const std::string& path);
// Seconds left from the fraction done and the seconds elapsed (-1: unknown)
double eta_seconds(double fraction, double elapsed_s);
double slurm_duration_seconds(const std::string& s);   // "1-02:03:04", "02:03:04", "03:04" → seconds; -1 unknown

// A thread-scaling check (caps job scaling): the same short MD at several thread counts on one node. The counts to try
// on a node of `cores`: powers of two below it, half the node and the whole node (0 cores: 1 … 32). The table from the
// measurements (each task's bench.json): threads, ns/day, the speedup over the fewest threads and the parallel
// efficiency (speedup / thread ratio), sorted by threads.
std::vector<int> scaling_threads(int cores);
Json scaling_table(const std::vector<Json>& bench);
std::string scaling_text(const Json& table);

// Everything a viewer needs since its last look, as JSON (caps job poll; the Studio parses this one object)
struct PollRequest { std::string dir; int64_t log_offset = 0, progress_offset = 0, err_offset = 0; size_t max_bytes = 262144; };
Json job_poll(const PollRequest& p, const std::string& squeue_line, const std::string& sacct_text);
// Where a job's live files are: the scratch folder while it runs, else out/ (else the folder itself)
std::string job_live_dir(const std::string& dir, const Json& state);

}  // namespace caps
