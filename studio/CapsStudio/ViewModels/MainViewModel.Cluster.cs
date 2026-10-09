using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Cluster jobs (any user's own SLURM or PBS cluster, or a workstation): the host's profile written as its
/// ~/CAPS/host.json, caps installed there by a short build job, runs sent with caps job new --submit and followed with
/// caps job poll over one shared ssh connection — the same commands a user types on the cluster, so a job started by
/// hand can be followed here too (Find my jobs). Keys come from the SSH agent only; nothing personal is stored beyond
/// what the user typed into the host's fields.</summary>
public sealed partial class MainViewModel
{
    // ---------------------------------------------------------------- the host's profile
    internal static string SchedKey(RemoteHost h) => h.Scheduler switch { "PBS" => "pbs", "none" => "direct", _ => "slurm" };

    /// <summary>The host's fields as caps job reads them (caps/cluster_job.hpp host_profile_json).</summary>
    internal static JsonObject HostProfileJson(RemoteHost h) => new()
    {
        ["format"] = "caps-host", ["version"] = 1, ["name"] = h.Name, ["ssh"] = Target(h), ["port"] = h.Port, ["scheduler"] = SchedKey(h),
        ["account"] = h.Account, ["partition"] = h.Partition, ["modules"] = h.Modules, ["build_modules"] = h.BuildModules, ["root"] = h.Root,
        ["scratch"] = h.Scratch, ["ws_filesystem"] = h.WsFilesystem, ["ws_days"] = h.WsDays, ["scratch_var"] = h.ScratchVar, ["cpus"] = h.Cpus,
        ["mem_per_cpu"] = h.MemPerCpu, ["time"] = h.TimeLimit, ["mail_type"] = h.MailType, ["mail_user"] = h.MailUser, ["keep_workspace"] = h.KeepWorkspace,
        ["node_cores"] = h.NodeCores, ["node_mem"] = h.NodeMem, ["job_template"] = h.JobScript,
    };

    /// <summary>A path under the host's root as a shell word: "~/CAPS/bin/caps" → "$HOME/CAPS/bin/caps" (expanded by the
    /// host's shell), anything else single-quoted.</summary>
    internal static string RootPath(RemoteHost h, string sub = "")
    {
        var root = (h.Root.Length > 0 ? h.Root : "~/CAPS").TrimEnd('/');
        var p = sub.Length > 0 ? root + "/" + sub : root;
        if (p == "~") return "\"$HOME\"";
        return p.StartsWith("~/", StringComparison.Ordinal) ? "\"$HOME/" + p[2..].Replace("\"", "\\\"").Replace("$", "\\$") + "\"" : Q(p);
    }
    internal static string CapsOnHost(RemoteHost h) => RootPath(h, "bin/caps");

    /// <summary>One ssh connection shared by the polls (ControlMaster; not on Windows, whose OpenSSH has none).</summary>
    private static string[] SharedConnection()
    {
        if (OperatingSystem.IsWindows()) return [];
        var dir = Path.Combine(Path.GetTempPath(), "caps-ssh");
        try { Directory.CreateDirectory(dir); } catch { return []; }
        return ["-o", "ControlMaster=auto", "-o", "ControlPersist=10m", "-o", "ControlPath=" + Path.Combine(dir, "%C")];
    }

    /// <summary>The profile written on the host (root/host.json): the CLI there and the Studio read the same file.</summary>
    private async Task WriteHostProfile(RemoteHost h)
    {
        var local = Path.Combine(RemoteFolder, "host-" + string.Concat(h.Name.Where(char.IsLetterOrDigit)) + ".json");
        Directory.CreateDirectory(RemoteFolder);
        File.WriteAllText(local, HostProfileJson(h).ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
        var mk = await Tool("ssh", SshArgs(h, $"mkdir -p {RootPath(h)} && cd {RootPath(h)} && pwd"), 30000);
        if (mk.Code != 0) throw new InvalidOperationException("ssh: " + FirstLine(mk.Err, mk.Code));
        var up = await Tool("scp", ScpArgs(h, [local], $"{Target(h)}:{mk.Out.Split('\n').Last().Trim()}/host.json"), 60000);
        if (up.Code != 0) throw new InvalidOperationException("scp: " + FirstLine(up.Err, up.Code));
    }
    private static string FirstLine(string err, int code) => err.Length > 0 ? err.Split('\n')[0] : $"exit {code}";

    // ---------------------------------------------------------------- presets and the host's fields
    public List<string> HostPresetNames { get; } = ["SLURM cluster", "SLURM with workspaces", "PBS cluster", "Workstation"];
    private static readonly string[] PresetKeys = ["slurm", "slurm-workspace", "pbs", "workstation"];

    /// <summary>Fills the host's job fields from a generic preset (the core's host_presets); the user then adds their
    /// own account, modules and limits. Name, hostname, user and port are kept.</summary>
    public void ApplyHostPreset(int index)
    {
        if (_host?.Host is not { } h || index < 0 || index >= PresetKeys.Length) return;
        try
        {
            var all = JsonNode.Parse(CapsDocument.HostPresets())!.AsArray();
            var p = all.OfType<JsonObject>().FirstOrDefault(x => (string?)x["key"] == PresetKeys[index])?["profile"] as JsonObject;
            if (p == null) return;
            h.Scheduler = (string?)p["scheduler"] switch { "pbs" => "PBS", "direct" => "none", _ => "SLURM" };
            h.Modules = (string?)p["modules"] ?? "";
            h.BuildModules = (string?)p["build_modules"] ?? "";
            h.Scratch = (string?)p["scratch"] ?? "none";
            h.WsFilesystem = (string?)p["ws_filesystem"] ?? "";
            h.WsDays = (int?)p["ws_days"] ?? 10;
            h.ScratchVar = (string?)p["scratch_var"] ?? "$SCRATCH";
            h.Cpus = (int?)p["cpus"] ?? 16;
            h.MemPerCpu = (string?)p["mem_per_cpu"] ?? "";
            h.TimeLimit = (string?)p["time"] ?? "";
            HostChanged();
            RaiseClusterFields();
            Status = $"{HostPresetNames[index]}: fill in your account, modules and limits";
        }
        catch (Exception e) { Status = "Preset: " + e.Message; }
    }

    private void RaiseClusterFields()
    {
        foreach (var n in new[] { nameof(HostAccount), nameof(HostModules), nameof(HostBuildModules), nameof(HostRoot), nameof(HostScratchIndex), nameof(HostWsFilesystem),
                                  nameof(HostWsDays), nameof(HostScratchVar), nameof(HostCpus), nameof(HostMem), nameof(HostTime), nameof(HostMailType), nameof(HostMailUser),
                                  nameof(HostKeepWorkspace), nameof(HostNodeText), nameof(HostUsesWorkspace), nameof(HostUsesScratchVar), nameof(HostJobScript),
                                  nameof(HostScheduler), nameof(HostIsQueue) })
            Raise(n);
    }
    private void SetHost(Action<RemoteHost> f) { if (_host?.Host is { } h) { f(h); HostChanged(); RaiseClusterFields(); } }
    public string HostAccount { get => H.Account; set => SetHost(h => h.Account = (value ?? "").Trim()); }
    public string HostModules { get => H.Modules; set => SetHost(h => h.Modules = value ?? ""); }
    public string HostBuildModules { get => H.BuildModules; set => SetHost(h => h.BuildModules = value ?? ""); }
    public string HostRoot { get => H.Root; set => SetHost(h => h.Root = (value ?? "").Trim() is { Length: > 0 } v ? v : "~/CAPS"); }
    public static readonly string[] ScratchKinds = ["Workspace (ws_allocate)", "A folder under a variable ($SCRATCH, $TMPDIR)", "None: run in the job folder"];
    public int HostScratchIndex
    {
        get => H.Scratch switch { "workspace" => 0, "env" => 1, _ => 2 };
        set => SetHost(h => h.Scratch = value switch { 0 => "workspace", 1 => "env", _ => "none" });
    }
    public bool HostUsesWorkspace => H.Scratch == "workspace";
    public bool HostUsesScratchVar => H.Scratch == "env";
    public bool HostIsQueue => H.Scheduler != "none";
    public string HostWsFilesystem { get => H.WsFilesystem; set => SetHost(h => h.WsFilesystem = (value ?? "").Trim()); }
    public decimal HostWsDays { get => H.WsDays; set => SetHost(h => h.WsDays = (int)Math.Clamp(value, 1, 365)); }
    public string HostScratchVar { get => H.ScratchVar; set => SetHost(h => h.ScratchVar = (value ?? "").Trim()); }
    public decimal HostCpus { get => H.Cpus; set => SetHost(h => h.Cpus = (int)Math.Clamp(value, 1, 4096)); }
    public string HostMem { get => H.MemPerCpu; set => SetHost(h => h.MemPerCpu = (value ?? "").Trim()); }
    public string HostTime { get => H.TimeLimit; set => SetHost(h => h.TimeLimit = (value ?? "").Trim()); }
    public string HostMailType { get => H.MailType; set => SetHost(h => h.MailType = (value ?? "").Trim()); }
    public string HostMailUser { get => H.MailUser; set => SetHost(h => h.MailUser = (value ?? "").Trim()); }
    public bool HostKeepWorkspace { get => H.KeepWorkspace; set => SetHost(h => h.KeepWorkspace = value); }
    public string HostNodeText => H.NodeCores > 0 ? $"{H.NodeCores} cores per node{(H.NodeMem.Length > 0 ? " · " + H.NodeMem + " MB" : "")} (from the host test)" : "node size: run Test connection";
    /// <summary>The host's job script: its own when edited, else the built-in one for its scheduler (Reset clears it).</summary>
    public string HostJobScript
    {
        get => H.JobScript.Length > 0 ? H.JobScript : SafeTemplate(SchedKey(H));
        set { if (_host?.Host is { } h && value != null && value != HostJobScript) { h.JobScript = value.Replace("\r\n", "\n"); Changed("Job script"); Raise(); } }
    }
    public void ResetHostJobScript() { if (_host?.Host is { } h) { h.JobScript = ""; Changed("Job script"); Raise(nameof(HostJobScript)); } }
    private static string SafeTemplate(string sched) { try { return CapsDocument.JobTemplate(sched); } catch { return ""; } }

    /// <summary>The Studio's own build: "0.1.0 (commit a1b2c3d4e5)".</summary>
    internal static (string Version, string Commit) StudioBuild()
    {
        try
        {
            var j = JsonNode.Parse(CapsDocument.BuildInfo())!;
            return ((string?)j["version"] ?? "", (string?)j["commit"] ?? "unknown");
        }
        catch { return ("", "unknown"); }
    }

    /// <summary>"caps 0.1.0 (commit X)" → X ("" when the line names none: an older caps).</summary>
    internal static string CommitOf(string versionLine)
    {
        var at = versionLine.IndexOf("(commit ", StringComparison.Ordinal);
        return at < 0 ? "" : versionLine[(at + 8)..].TrimEnd(')', ' ', '\n');
    }

    /// <summary>The host test's script: caps's version (the installed one, else PATH), the node size, the scheduler and
    /// whether workspaces exist. One line each, keyed.</summary>
    internal static string HostProbeCommand(RemoteHost h) =>
        $"echo \"version:$({CapsOnHost(h)} --version 2>/dev/null || caps --version 2>/dev/null || echo none)\"; " +
        "echo \"node:$(sinfo -h -o '%c %m' 2>/dev/null | sort -n | tail -n 1)\"; " +
        "command -v sbatch >/dev/null 2>&1 && echo sched:slurm; command -v qsub >/dev/null 2>&1 && echo sched:pbs; " +
        "command -v ws_allocate >/dev/null 2>&1 && echo ws:yes; command -v git >/dev/null 2>&1 && echo git:yes; command -v cmake >/dev/null 2>&1 && echo cmake:yes";

    internal static Dictionary<string, string> ParseProbe(string text)
    {
        var d = new Dictionary<string, string>();
        foreach (var line in text.Split('\n'))
        {
            var c = line.IndexOf(':');
            if (c > 0) d[line[..c].Trim()] = line[(c + 1)..].Trim();
        }
        return d;
    }

    // ---------------------------------------------------------------- installing caps on the host
    public const string SourceRepository = "https://github.com/MuhammadUzairRiaz/CAPS.git";

    /// <summary>The build script run as a short job (never on a login node when there is a queue): the source at the
    /// Studio's commit (a clone of the public repository), cmake with the host's build modules, bin/caps linked; its
    /// state in root/build.state, its output in root/build.log.</summary>
    internal static string BuildScript(RemoteHost h, string commit)
    {
        var root = RootPath(h);
        var checkout = commit.Length > 0 && commit != "unknown" ? $"git checkout --quiet {Q(commit)} 2>/dev/null || git checkout --quiet main" : "git checkout --quiet main && git pull --quiet --ff-only";
        return "#!/bin/bash\n# CAPS: build the caps command line on this host (written by CAPS Studio's Install / update caps)\n" +
               "set -u\n" +
               $"ROOT={root}\n" +
               "echo building > \"$ROOT/build.state\"\n" +
               "exec > \"$ROOT/build.log\" 2>&1\n" +
               "fail() { echo \"failed: $1\"; echo \"failed $1\" > \"$ROOT/build.state\"; exit 1; }\n" +
               (h.BuildModules.Trim().Length > 0 ? h.BuildModules.Replace("\r\n", "\n").Trim() + "\n" : "") +
               "mkdir -p \"$ROOT/src\" \"$ROOT/bin\" || fail mkdir\n" +
               "cd \"$ROOT/src\" || fail cd\n" +
               $"[ -d CAPS/.git ] || git clone --quiet {SourceRepository} CAPS || fail clone\n" +
               "cd CAPS && git fetch --quiet origin || fail fetch\n" +
               $"{checkout} || fail checkout\n" +
               "cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCAPS_BUILD_TESTS=OFF || fail cmake\n" +
               "cmake --build build -j \"${SLURM_CPUS_PER_TASK:-${NCPUS:-4}}\" --target caps_cli || fail build\n" +
               "ln -sf \"$ROOT/src/CAPS/build/cli/caps\" \"$ROOT/bin/caps\" || fail link\n" +
               "V=$(\"$ROOT/bin/caps\" --version) || fail run\n" +
               "echo \"$V\"; echo \"ok $V\" > \"$ROOT/build.state\"\n";
    }

    /// <summary>How the build script is started: a short job on the queue, or in the background on a workstation.</summary>
    internal static string BuildSubmit(RemoteHost h)
    {
        var root = RootPath(h);
        var acct = h.Account.Length > 0 ? (h.Scheduler == "PBS" ? $" -A {Q(h.Account)}" : $" --account={Q(h.Account)}") : "";
        var part = h.Partition.Length > 0 ? (h.Scheduler == "PBS" ? $" -q {Q(h.Partition)}" : $" --partition={Q(h.Partition)}") : "";
        return h.Scheduler switch
        {
            "SLURM" => $"cd {root} && sbatch --parsable --job-name=caps-build --nodes=1 --ntasks=1 --cpus-per-task=8 --time=00:45:00{acct}{part} --output=build-%j.out build-caps.sh",
            "PBS" => $"cd {root} && qsub -N caps-build -l select=1:ncpus=8 -l walltime=00:45:00{acct}{part} build-caps.sh",
            _ => $"cd {root} && (nohup bash build-caps.sh > /dev/null 2>&1 & echo pid:$!)",
        };
    }

    private bool _installing;
    public bool HostInstallIdle { get => !_installing; private set => Set(ref _installing, !value); }
    private string _installText = "";
    public string HostInstallText { get => _installText; private set => Set(ref _installText, value); }

    /// <summary>Install / update caps: the build script goes up and runs as a job; Jobs follows it (build.state).</summary>
    public async Task InstallCapsOnHost()
    {
        if (_host?.Host is not { } h || _installing) return;
        if (h.Hostname.Length == 0) { HostInstallText = "Enter the hostname first"; return; }
        HostInstallIdle = false;
        HostInstallText = "Sending the build…";
        var commit = StudioBuild().Commit;
        var k = _jobCounters["Install"] = _jobCounters.GetValueOrDefault("Install") + 1;
        var job = new Job { Id = $"install-{k}", Kind = "Install", Module = 10, Title = $"Install caps · on {h.Name}", Document = "caps " + commit, Atoms = 0 };
        job.Remote = new RemoteRun { Host = h.Name, Scheduler = h.Scheduler, Mode = "install", Local = Path.Combine(RemoteFolder, $"install-{k}-{DateTime.Now:yyyyMMdd-HHmmss}") };
        job.Status = "queued";
        Jobs.Insert(0, job);
        Raise(nameof(HasJobs)); Raise(nameof(JobsSummary));
        try
        {
            await WriteHostProfile(h);
            Directory.CreateDirectory(job.Remote.Local);
            var script = Path.Combine(job.Remote.Local, "build-caps.sh");
            File.WriteAllText(script, BuildScript(h, commit));
            var root = (await Tool("ssh", SshArgs(h, $"cd {RootPath(h)} && pwd"), 30000)).Out.Split('\n').Last().Trim();
            job.Remote.Dir = root;
            var up = await Tool("scp", ScpArgs(h, [script], $"{Target(h)}:{root}/build-caps.sh"), 60000);
            if (up.Code != 0) throw new InvalidOperationException("scp: " + FirstLine(up.Err, up.Code));
            var sub = await Tool("ssh", SshArgs(h, BuildSubmit(h)), 60000);
            if (sub.Code != 0 || sub.Out.Length == 0) throw new InvalidOperationException("submit: " + FirstLine(sub.Err, sub.Code));
            job.Remote.JobId = sub.Out.Split('\n').Last().Split(';')[0].Trim();
            job.Status = "running";
            job.Add($"Building caps {commit} on {h.Name} ({(h.Scheduler == "none" ? "in the background" : h.Scheduler + " job " + job.Remote.JobId)}): the source from {SourceRepository}");
            HostInstallText = "Building on the host: Jobs follows it (a few minutes)";
            StartRemotePoll();
        }
        catch (Exception e)
        {
            job.Status = "failed";
            job.Error = e.Message;
            job.Add("Could not start the build: " + e.Message);
            job.Ended = DateTime.Now;
            HostInstallText = "Could not start the build: " + e.Message;
        }
        finally { HostInstallIdle = true; SaveJobs(); Raise(nameof(JobsSummary)); }
    }

    /// <summary>The build's state (build.state) and the end of its log.</summary>
    private async Task PollInstall(Job j, RemoteHost h, RemoteRun r)
    {
        var st = await Tool("ssh", [.. SharedConnection(), .. SshArgs(h, $"cat {RootPath(h, "build.state")} 2>/dev/null; echo ---; tail -n 15 {RootPath(h, "build.log")} 2>/dev/null")], 30000);
        var parts = st.Out.Split("---", 2);
        var state = parts[0].Trim();
        if (parts.Length > 1 && parts[1].Trim().Length > 0) { j.AppendLive("", ""); j.SetQueue(parts[1].Trim()); }
        if (state.StartsWith("ok", StringComparison.Ordinal))
        {
            j.Status = "done";
            j.Progress = 1;
            j.Ended = DateTime.Now;
            h.CapsVersion = state[2..].Trim();
            j.Add($"Installed: {h.CapsVersion} in {r.Dir}/bin/caps");
            Changed("caps installed");
            if (_host?.Host == h) { HostInstallText = "Installed: " + h.CapsVersion; _ = TestHost(); }
        }
        else if (state.StartsWith("failed", StringComparison.Ordinal))
        {
            j.Status = "failed";
            j.Ended = DateTime.Now;
            j.Error = "The build " + state + (parts.Length > 1 ? ":\n" + parts[1].Trim() : "");
            j.Add(j.Error);
            j.Suggestion = "Check the build modules (a C++20 compiler, gfortran and CMake ≥ 3.24) in Settings › Compute & remote, then Install again.";
            j.SuggestModule = 10;
        }
        SaveJobs();
    }

    // ---------------------------------------------------------------- the run's settings for a cluster
    private string _rTitle = "";
    /// <summary>The job's title (its folder on the host): the structure's name, editable.</summary>
    public string RemoteTitle { get => _rTitle.Length > 0 ? _rTitle : DefaultRemoteTitle; set { _rTitle = value ?? ""; Raise(); Raise(nameof(RemoteCommand)); } }
    private string DefaultRemoteTitle => string.Concat(Title.Replace(" (unsaved)", "").Select(c => char.IsLetterOrDigit(c) || c is '-' or '_' or '.' ? c : '_')).Trim('_') is { Length: > 0 } t ? t : "structure";
    private int _rCpus;
    public decimal RemoteCpus
    {
        get => _rCpus > 0 ? _rCpus : RemoteHostNow?.Cpus ?? 16;
        set { var cap = RemoteHostNow?.NodeCores is > 0 and var n ? n : 4096; _rCpus = (int)Math.Clamp(value, 1, cap); Raise(); Raise(nameof(RemoteCommand)); }
    }
    private string _rMem = "", _rTime = "";
    public string RemoteMem { get => _rMem.Length > 0 ? _rMem : RemoteHostNow?.MemPerCpu ?? ""; set { _rMem = (value ?? "").Trim(); Raise(); Raise(nameof(RemoteCommand)); } }
    public string RemoteTime { get => _rTime.Length > 0 ? _rTime : RemoteHostNow?.TimeLimit ?? ""; set { _rTime = (value ?? "").Trim(); Raise(); Raise(nameof(RemoteCommand)); } }
    private double _rFramePs = 10, _rThermoPs = 1;
    /// <summary>A frame of the run every this many ps on the host (0: none), for the live Structure view.</summary>
    public decimal RemoteFramePs { get => (decimal)_rFramePs; set { _rFramePs = (double)Math.Max(0, value); Raise(); } }
    public decimal RemoteThermoPs { get => (decimal)_rThermoPs; set { _rThermoPs = (double)Math.Max(0.001m, value); Raise(); } }
    private bool? _rKeep;
    public bool RemoteKeepWorkspace { get => _rKeep ?? RemoteHostNow?.KeepWorkspace ?? false; set { _rKeep = value; Raise(); Raise(nameof(RemoteCommand)); } }
    private RemoteHost? RemoteHostNow => RunWhereIndex > 0 && RunWhereIndex <= _settings.Hosts.Count ? _settings.Hosts[RunWhereIndex - 1] : null;
    public string RemoteHostLabel => RemoteHostNow is { } h ? $"{(h.Scheduler == "none" ? "no queue" : h.Scheduler)} · {(h.NodeCores > 0 ? h.NodeCores + " cores per node" : "node size unknown")} · CAPS runs threads, not MPI ranks: one node, one task" : "";

    /// <summary>What will run on the host: the exact command a user could type there.</summary>
    public string RemoteCommand
    {
        get
        {
            if (RemoteHostNow is not { } h) return "";
            var inv = CultureInfo.InvariantCulture;
            var opts = $"--cpus {RemoteCpus.ToString("0", inv)}" + (RemoteMem.Length > 0 ? $" --mem {RemoteMem}" : "") + (RemoteTime.Length > 0 ? $" --time {RemoteTime}" : "") +
                       (RemoteKeepWorkspace ? " --keep-workspace" : "");
            return $"{h.Root}/bin/caps job new --title {RemoteTitle} --kind <run> --input structure.caps.data,recipe.json {opts} --submit -- caps run recipe.json --out .";
        }
    }
    private void RaiseRemoteRun()
    {
        foreach (var n in new[] { nameof(RemoteCpus), nameof(RemoteMem), nameof(RemoteTime), nameof(RemoteKeepWorkspace), nameof(RemoteCommand), nameof(RemoteHostLabel), nameof(RemoteTitle) }) Raise(n);
    }

    /// <summary>The recipe's stages write frames and thermo rows on the host (the live views read them).</summary>
    private void AddLiveOutput(JsonObject recipe)
    {
        if (recipe["md"] is JsonObject md)
        {
            if (_rFramePs > 0) md["frame_ps"] = _rFramePs;
            md["thermo_every"] = Math.Max(1, (int)Math.Round(_rThermoPs * 1000 / Math.Max(0.1, (double?)md["dt"] ?? 1.0)));
        }
        if (recipe["equilibrate"] is JsonObject eq)
        {
            if (_rFramePs > 0) eq["frame_ps"] = _rFramePs;
            eq["thermo_ps"] = _rThermoPs;
        }
    }

    // ---------------------------------------------------------------- sending a run with caps job
    /// <summary>Uploads the inputs to a staging folder under the root, writes the host's profile, and makes and submits the
    /// job there with caps job new; the job folder and the scheduler's id come back.</summary>
    private async Task SendJob(Job job, RemoteHost h, string kind, string[] files, string command)
    {
        var r = job.Remote!;
        r.Mode = "job";
        await WriteHostProfile(h);
        var stage = $".uploads/{job.Id}-{DateTime.Now:yyyyMMddHHmmss}";
        var mk = await Tool("ssh", SshArgs(h, $"mkdir -p {RootPath(h, stage)} && cd {RootPath(h, stage)} && pwd"), 30000);
        if (mk.Code != 0) throw new InvalidOperationException("ssh: " + FirstLine(mk.Err, mk.Code));
        var dir = mk.Out.Split('\n').Last().Trim();
        var folders = files.Any(f => Directory.Exists(Path.Combine(r.Local, f)));
        var up = await Tool("scp", ScpArgs(h, files.Select(f => Path.Combine(r.Local, f)), $"{Target(h)}:{dir}/", folders), 300000);
        if (up.Code != 0) throw new InvalidOperationException("scp: " + FirstLine(up.Err, up.Code));
        job.Add($"Uploaded {string.Join(", ", files)} to {h.Name}");
        var inv = CultureInfo.InvariantCulture;
        var title = string.Concat(RemoteTitle.Select(c => char.IsLetterOrDigit(c) || c is '-' or '_' or '.' ? c : '_'));
        var opts = $"--title {Q(title)} --kind {Q(kind.ToLowerInvariant())} --input {Q(string.Join(",", files.Select(f => dir + "/" + f)))} --cpus {RemoteCpus.ToString("0", inv)}" +
                   (RemoteMem.Length > 0 ? $" --mem {Q(RemoteMem)}" : "") + (RemoteTime.Length > 0 ? $" --time {Q(RemoteTime)}" : "") +
                   (RemoteKeepWorkspace ? " --keep-workspace" : "") + (h.Partition.Length > 0 ? $" --partition {Q(h.Partition)}" : "");
        var cmd = $"{CapsOnHost(h)} job new --host-profile {RootPath(h, "host.json")} {opts} --json --submit -- {command}; RC=$?; rm -rf {Q(dir)}; exit $RC";
        var sub = await Tool("ssh", SshArgs(h, cmd), 120000);
        var line = sub.Out.Split('\n').LastOrDefault(l => l.TrimStart().StartsWith('{')) ?? "";
        JsonObject? res = null;
        try { res = JsonNode.Parse(line) as JsonObject; } catch { }
        if (res == null)
            throw new InvalidOperationException(sub.Out.Contains("No such file", StringComparison.Ordinal) || sub.Err.Contains("No such file", StringComparison.Ordinal) || sub.Code == 127
                ? $"caps is not installed in {h.Root}/bin on {h.Name}: Settings › Compute & remote › Install / update caps"
                : "caps job new: " + FirstLine(sub.Err.Length > 0 ? sub.Err : sub.Out, sub.Code));
        r.Dir = (string?)res["dir"] ?? "";
        if (res["submitted"]?.GetValue<bool>() != true) throw new InvalidOperationException("the scheduler refused the job: " + ((string?)res["error"] ?? "").Trim());
        r.JobId = (string?)res["id"] ?? "";
        job.Status = "queued";
        job.Add($"Submitted · {(h.Scheduler == "none" ? "process" : h.Scheduler + " job")} {r.JobId} · {r.Dir}");
        job.Add($"On the host: {h.Root}/bin/caps job status {r.Dir}");
    }

    // ---------------------------------------------------------------- following a job
    /// <summary>One look at a cluster job: caps job poll on the host (the queue, the state file, what is new in the log
    /// and the progress lines, the files, the last frame), through the shared connection. An unreachable host never
    /// marks a job failed.</summary>
    private async Task PollJob(Job j, RemoteHost h, RemoteRun r)
    {
        if (r.Checking) return;
        r.Checking = true;
        r.LastPoll = DateTime.Now;
        try
        {
            var cmd = $"{CapsOnHost(h)} job poll {Q(r.Dir)} --host-profile {RootPath(h, "host.json")} --log-offset {r.LogOffset} --progress-offset {r.ProgressOffset} --err-offset {r.ErrOffset}";
            var res = await Tool("ssh", [.. SharedConnection(), .. SshArgs(h, cmd)], 60000);
            if (res.Code != 0 || JsonNode.Parse(res.Out.Split('\n').Last(l => l.Length > 0)) is not JsonObject p)
            {
                j.Add($"{r.Host} not reached: {FirstLine(res.Err, res.Code)} (the job keeps running there)");
                return;
            }
            ApplyPoll(j, r, p);
            var state = (string?)p["state"] ?? "";
            if (state is "finished" or "failed" or "timeout" or "cancelled")
            {
                j.Add(state switch
                {
                    "finished" => $"Finished on {r.Host}",
                    "timeout" => $"Stopped at its time limit on {r.Host}: a checkpoint is in the results; Resume continues it",
                    "cancelled" => $"Cancelled on {r.Host}",
                    _ => $"Failed on {r.Host}" + ((string?)p["acct"]?["exit_code"] is { Length: > 0 } x ? $" (exit {x})" : ""),
                });
                if (state == "cancelled") { j.Status = "cancelled"; j.Ended = DateTime.Now; }
                else await BringBack(j, h, r, state == "finished");
                if (state == "timeout" && j.Status == "failed") { j.Status = "stopped"; j.Error = ""; }
            }
            SaveJobs();
            Raise(nameof(JobsSummary)); Raise(nameof(ComputeText));
        }
        catch (Exception e) { j.Add($"{r.Host} not reached: {e.Message} (the job keeps running there)"); }
        finally { r.Checking = false; }
    }

    /// <summary>What a poll says, onto the job: the queue panel, the log, the curves and the progress.</summary>
    internal static void ApplyPoll(Job j, RemoteRun r, JsonObject p)
    {
        var q = p["queue"] as JsonObject;
        var found = q?["found"]?.GetValue<bool>() == true;
        r.QueueState = found ? (string?)q!["state"] ?? "" : "";
        r.Reason = found ? (string?)q!["reason"] ?? "" : "";
        r.Node = found ? (string?)q!["node"] ?? "" : r.Node;
        r.Elapsed = found ? (string?)q!["elapsed"] ?? "" : r.Elapsed;
        r.Limit = found ? (string?)q!["limit"] ?? "" : r.Limit;
        if (p["acct"] is JsonObject a) r.Acct = $"{(string?)a["state"]} · exit {(string?)a["exit_code"]} · {(string?)a["elapsed"]}{((string?)a["max_rss"] is { Length: > 0 } m ? " · MaxRSS " + m : "")}";
        r.LiveDir = (string?)p["live_dir"] ?? r.LiveDir;
        r.LastFrameOffset = (long?)p["last_frame_offset"] ?? -1;
        r.LogOffset = (long?)p["log_offset"] ?? r.LogOffset;
        r.ProgressOffset = (long?)p["progress_offset"] ?? r.ProgressOffset;
        if (p["err_offset"] != null) r.ErrOffset = (long?)p["err_offset"] ?? r.ErrOffset;
        r.Eta = (double?)p["eta_s"] ?? -1;
        var state = (string?)p["state"] ?? "";
        if (state == "pending" || state == "submitted") j.Status = "queued";
        else if (state == "running") j.Status = "running";
        j.AppendLive((string?)p["log"] ?? "", (string?)p["err"] ?? "");
        // thermo from the progress lines: T and ρ against time (the two curve panels)
        var any = false;
        if (p["progress"] is JsonArray lines)
            foreach (var l in lines.OfType<JsonObject>())
            {
                if ((double?)l["fraction"] is { } f) j.Progress = f;
                if ((int?)l["stages"] is > 0 and var n) { j.Stages = n; j.Stage = (int?)l["stage_index"] ?? 0; }
                if (l["time_ps"] is null) continue;
                var t = (double?)l["time_ps"] ?? 0;
                if ((double?)l["T"] is { } temp) j.A.Add((t, temp));
                if ((double?)l["rho"] is { } rho) j.B.Add((t, rho));
                any = true;
            }
        if (any)
        {
            if (j.CurveA.Length == 0) { j.CurveA = "Temperature"; j.AxisA = "K"; j.CurveB = "Density"; j.AxisB = "g/cm³"; j.AxisX = "time / ps"; }
            if (j.A.Count > 20000) j.A.RemoveRange(0, j.A.Count - 20000);
            if (j.B.Count > 20000) j.B.RemoveRange(0, j.B.Count - 20000);
            j.CurvesChanged();
        }
        var files = p["files"] is JsonArray fa ? fa.OfType<JsonObject>().Select(x => new RemoteFile((string?)x["name"] ?? "", (long?)x["bytes"] ?? 0)).Where(x => x.Name.Length > 0).ToList() : new();
        var lines2 = new List<string>
        {
            $"state      {state}{(r.QueueState.Length > 0 ? $"  ({r.QueueState}{(r.Reason.Length > 0 ? ", " + r.Reason : "")})" : "")}",
        };
        if (r.Node.Length > 0) lines2.Add($"node       {r.Node}");
        if (r.Elapsed.Length > 0) lines2.Add($"elapsed    {r.Elapsed}{(r.Limit.Length > 0 ? " of " + r.Limit : "")}");
        if (r.Eta >= 0) lines2.Add($"ETA        {TimeSpan.FromSeconds(r.Eta):d\\.hh\\:mm\\:ss}");
        if (r.Acct.Length > 0 && !found) lines2.Add($"accounting {r.Acct}");
        lines2.Add($"files      {r.LiveDir}");
        foreach (var f in files.OrderByDescending(f => f.Bytes).Take(6)) lines2.Add($"           {f.Name}  {(f.Bytes < 1 << 20 ? (f.Bytes / 1024.0).ToString("F1", CultureInfo.InvariantCulture) + " kB" : f.Size)}");
        j.SetQueue(string.Join("\n", lines2));
    }

    /// <summary>Asks the host now (Refresh); the timer asks every 30 s for the job shown in Jobs, every 5 min otherwise.</summary>
    public Task RefreshJob(Job? j) => j?.Remote is { } r && _settings.Hosts.FirstOrDefault(x => x.Name == r.Host) is { } h
        ? r.Mode switch { "job" => PollJob(j, h, r), "install" => PollInstall(j, h, r), _ => CheckRemote(j) }
        : Task.CompletedTask;

    private bool PollDue(Job j)
    {
        if (j.Remote is not { } r) return false;
        var shown = IsJobs && SelectedJob == j;
        return (DateTime.Now - r.LastPoll).TotalSeconds >= (shown ? 25 : 290);
    }

    /// <summary>Cancel on the host (after the page's confirmation): caps job cancel.</summary>
    public async Task CancelRemote(Job? j)
    {
        if (j?.Remote is not { Mode: "job" } r || _settings.Hosts.FirstOrDefault(x => x.Name == r.Host) is not { } h) return;
        var res = await Tool("ssh", [.. SharedConnection(), .. SshArgs(h, $"{CapsOnHost(h)} job cancel {Q(r.Dir)} --host-profile {RootPath(h, "host.json")}")], 60000);
        j.Add(res.Code == 0 ? $"Cancel sent to {r.Host}" : $"Cancel: {FirstLine(res.Err.Length > 0 ? res.Err : res.Out, res.Code)}");
        await PollJob(j, h, r);
    }

    /// <summary>A job stopped near its time limit, again from its checkpoint (caps job resume on the host): a new job of
    /// the same title and kind, followed here like the first.</summary>
    public async Task ResumeRemote(Job? j)
    {
        if (j?.Remote is not { Mode: "job" } r || _settings.Hosts.FirstOrDefault(x => x.Name == r.Host) is not { } h) return;
        try
        {
            var res = await Tool("ssh", [.. SharedConnection(), .. SshArgs(h, $"{CapsOnHost(h)} job resume {Q(r.Dir)} --host-profile {RootPath(h, "host.json")} --submit --json")], 120000);
            var line = res.Out.Split('\n').LastOrDefault(l => l.TrimStart().StartsWith('{')) ?? "";
            if (JsonNode.Parse(line) is not JsonObject o || o["submitted"]?.GetValue<bool>() != true)
            {
                j.Add("Resume: " + FirstLine(res.Err.Length > 0 ? res.Err : res.Out, res.Code));
                return;
            }
            var k = _jobCounters[j.Kind] = _jobCounters.GetValueOrDefault(j.Kind) + 1;
            var job = new Job { Id = $"{j.Kind.ToLowerInvariant()}-{k}", Kind = j.Kind, Module = j.Module, Title = j.Title.Replace(" · resumed", "") + " · resumed", Document = j.Document, Atoms = j.Atoms };
            job.Remote = new RemoteRun
            {
                Host = r.Host, Scheduler = r.Scheduler, Mode = "job", Stem = r.Stem, Dir = (string?)o["dir"] ?? "", JobId = (string?)o["id"] ?? "",
                Local = Path.Combine(RemoteFolder, $"{job.Id}-{DateTime.Now:yyyyMMdd-HHmmss}"),
            };
            Directory.CreateDirectory(job.Remote.Local);
            foreach (var f in new[] { "structure.caps.data" })   // the structure that went up, for the latest-frame view
                if (File.Exists(Path.Combine(r.Local, f))) File.Copy(Path.Combine(r.Local, f), Path.Combine(job.Remote.Local, f), true);
            job.Status = "queued";
            job.Add($"Resumes {j.Id} from its checkpoint · {h.Scheduler} {job.Remote.JobId} · {job.Remote.Dir}");
            j.Add($"Resumed as {job.Id} ({job.Remote.Dir})");
            Jobs.Insert(0, job);
            SelectedJob = job;
            Raise(nameof(HasJobs)); Raise(nameof(JobsSummary));
            SaveJobs();
            StartRemotePoll();
        }
        catch (Exception e) { j.Add($"{r.Host} not reached: {e.Message}"); }
    }

    /// <summary>The terminal commands that show the same job by hand: ssh, then the live log.</summary>
    public string TerminalCommands(Job? j)
    {
        if (j?.Remote is not { } r || _settings.Hosts.FirstOrDefault(x => x.Name == r.Host) is not { } h) return "";
        var port = h.Port != 22 ? $" -p {h.Port}" : "";
        return $"ssh{port} {Target(h)}\n{h.Root}/bin/caps job status {r.Dir}\n{h.Root}/bin/caps job tail {r.Dir} -f\n";
    }

    /// <summary>The run's latest frame from the host (from the last frame's byte offset only), opened over the structure
    /// that went up.</summary>
    public async Task LatestFrame(Job? j)
    {
        if (j?.Remote is not { Mode: "job" } r || _settings.Hosts.FirstOrDefault(x => x.Name == r.Host) is not { } h) return;
        if (r.LastFrameOffset < 0) { Status = "No frame yet: frames are written every " + RemoteFramePs.ToString("0.###", CultureInfo.InvariantCulture) + " ps"; return; }
        var local = Path.Combine(r.Local, "latest.lammpstrj");
        var cmd = $"tail -c +{r.LastFrameOffset + 1} {Q(r.LiveDir + "/traj.lammpstrj")}";
        var res = await Tool("ssh", [.. SharedConnection(), .. SshArgs(h, cmd)], 120000);
        if (res.Code != 0 || !res.Out.StartsWith("ITEM: TIMESTEP", StringComparison.Ordinal)) { Status = "Latest frame: " + FirstLine(res.Err, res.Code); return; }
        File.WriteAllText(local, res.Out + "\n");
        var top = Path.Combine(r.Local, "structure.caps.data");
        var step = res.Out.Split('\n').Skip(1).FirstOrDefault()?.Trim() ?? "";
        j.Add($"Latest frame (step {step}) from {r.Host}");
        Open(local, File.Exists(top) ? top : null);
    }

    /// <summary>Jobs started by hand on the host (caps job new on the cluster): listed with caps job list and followed here.</summary>
    public async Task FindHostJobs()
    {
        if (_host?.Host is not { } h) return;
        HostTestText = "Looking for jobs…";
        try
        {
            var res = await Tool("ssh", SshArgs(h, $"{CapsOnHost(h)} job list --json --host-profile {RootPath(h, "host.json")}"), 60000);
            if (res.Code != 0 || JsonNode.Parse(res.Out.Split('\n').Last(l => l.Length > 0)) is not JsonArray all) { HostTestText = "caps job list: " + FirstLine(res.Err, res.Code); return; }
            var known = Jobs.Where(x => x.Remote?.Mode == "job").Select(x => x.Remote!.Dir).ToHashSet();
            var added = 0;
            foreach (var s in all.OfType<JsonObject>())
            {
                var dir = (string?)s["dir"] ?? "";
                if (dir.Length == 0 || known.Contains(dir)) continue;
                var name = (string?)s["job"] ?? "job";
                var kind = name.Split('-')[0];
                var state = (string?)s["state"] ?? "";
                var job = new Job
                {
                    Id = $"{name}@{(string?)s["title"]}", Kind = char.ToUpperInvariant(kind[0]) + kind[1..], Module = kind switch { "md" or "dynamics" => 3, "equilibrate" => 4, "relax" => 2, "react" => 6, _ => 0 },
                    Title = $"{(string?)s["title"]} / {name} · on {h.Name}", Document = (string?)s["title"] ?? "", Atoms = 0,
                };
                job.Remote = new RemoteRun { Host = h.Name, Scheduler = h.Scheduler, Mode = "job", Dir = dir, JobId = (string?)s["slurm_job"] ?? "", Local = Path.Combine(RemoteFolder, $"{name}-{DateTime.Now:yyyyMMdd-HHmmss}-{added}"), Stem = "md" };
                Directory.CreateDirectory(job.Remote.Local);
                job.Status = state is "finished" ? "done" : state is "failed" ? "failed" : state is "cancelled" ? "cancelled" : state is "running" ? "running" : "queued";
                job.Add($"Found on {h.Name}: {dir} ({state})");
                Jobs.Insert(0, job);
                added++;
            }
            HostTestText = added == 0 ? "No new jobs on the host" : $"{added} job{(added == 1 ? "" : "s")} from the host added to Jobs";
            if (added > 0) { Raise(nameof(HasJobs)); Raise(nameof(JobsSummary)); SaveJobs(); StartRemotePoll(); }
        }
        catch (Exception e) { HostTestText = e.Message; }
    }
}
