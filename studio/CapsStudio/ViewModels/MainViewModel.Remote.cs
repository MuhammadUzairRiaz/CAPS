using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>A row of the host list: this machine or a remote host, with its last known state.</summary>
public sealed class HostRow : ObservableObject
{
    public HostRow(RemoteHost? host, string name, string detail) { Host = host; _name = name; _detail = detail; }
    public RemoteHost? Host { get; }
    private string _name, _detail, _state = "not tested";
    private int _level;   // 0 unknown, 1 ok, 2 error
    public string Name { get => _name; set => Set(ref _name, value); }
    public string Detail { get => _detail; set => Set(ref _detail, value); }
    public string State { get => _state; set => Set(ref _state, value); }
    public int Level { get => _level; set { if (Set(ref _level, value)) { Raise(nameof(IsOk)); Raise(nameof(IsError)); } } }
    public bool IsOk => _level == 1;
    public bool IsError => _level == 2;
    public bool IsLocal => Host == null;
}

/// <summary>Compute &amp; remote (design/boards/RemoteCompute): this machine and SSH hosts (credentials stay in the SSH
/// agent), a connection test that runs caps on the host, the batch-job template and what moves with a job.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] Schedulers = ["SLURM", "PBS", "none"];
    public ObservableCollection<HostRow> Hosts { get; } = new();
    private HostRow? _host;

    private void LoadHosts()
    {
        Hosts.Clear();
        Hosts.Add(new HostRow(null, "This machine", $"local · {Environment.ProcessorCount} threads · CPU") { State = "ready", Level = 1 });
        foreach (var h in _settings.Hosts) Hosts.Add(new HostRow(h, h.Name, HostDetail(h)));
        SelectedHost = Hosts.Count > 1 ? Hosts[1] : Hosts[0];
        Raise(nameof(HostsStatus));
    }

    private static string HostDetail(RemoteHost h) =>
        $"{(h.User.Length > 0 ? h.User + "@" : "")}{(h.Hostname.Length > 0 ? h.Hostname : "no hostname")}:{h.Port} · SSH · {(h.Scheduler == "none" ? "no scheduler" : h.Scheduler)}";

    public HostRow? SelectedHost
    {
        get => _host;
        set
        {
            if (!Set(ref _host, value)) return;
            foreach (var n in new[] { nameof(HostEditable), nameof(HostName), nameof(HostHostname), nameof(HostUser), nameof(HostPort), nameof(HostScheduler),
                                      nameof(HostPartition), nameof(HostWorkDir), nameof(HostTitle), nameof(HostTestText) }) Raise(n);
            RaiseClusterFields();
        }
    }
    public bool HostEditable => _host?.Host != null;
    public string HostTitle => _host?.Name ?? "";
    private RemoteHost H => _host?.Host ?? new RemoteHost();
    private void HostChanged()
    {
        if (_host?.Host is not { } h) return;
        _host.Name = h.Name;
        _host.Detail = HostDetail(h);
        Raise(nameof(RunWhereChoices));
        _host.State = "not tested";
        _host.Level = 0;
        foreach (var n in new[] { nameof(HostTitle), nameof(HostName), nameof(HostHostname), nameof(HostUser), nameof(HostPort), nameof(HostScheduler), nameof(HostPartition), nameof(HostWorkDir) })
            Raise(n);
        Changed("Host");
    }
    public string HostName { get => H.Name; set { if (_host?.Host is { } h && value != null) { h.Name = value; HostChanged(); } } }
    public string HostHostname { get => H.Hostname; set { if (_host?.Host is { } h && value != null) { h.Hostname = value.Trim(); HostChanged(); } } }
    public string HostUser { get => H.User; set { if (_host?.Host is { } h && value != null) { h.User = value.Trim(); HostChanged(); } } }
    public decimal HostPort { get => H.Port; set { if (_host?.Host is { } h) { h.Port = (int)Math.Clamp(value, 1, 65535); HostChanged(); } } }
    public string HostScheduler { get => H.Scheduler; set { if (_host?.Host is { } h && value != null) { h.Scheduler = value; HostChanged(); } } }
    public string HostPartition { get => H.Partition; set { if (_host?.Host is { } h && value != null) { h.Partition = value.Trim(); HostChanged(); } } }
    public string HostWorkDir { get => H.WorkDir; set { if (_host?.Host is { } h && value != null) { h.WorkDir = value.Trim(); HostChanged(); } } }
    public string HostsStatus => $"{Hosts.Count(h => h.IsOk)} of {Hosts.Count} hosts reachable";

    public void AddHost()
    {
        var h = new RemoteHost { Name = $"host-{_settings.Hosts.Count + 1}" };
        _settings.Hosts.Add(h);
        var row = new HostRow(h, h.Name, HostDetail(h));
        Hosts.Add(row);
        SelectedHost = row;
        Changed("Host added");
        Raise(nameof(RunWhereChoices)); Raise(nameof(HasHosts)); Raise(nameof(RunWhereIndex));
    }

    public void RemoveHost()
    {
        if (_host?.Host is not { } h) return;
        _settings.Hosts.Remove(h);
        Hosts.Remove(_host);
        SelectedHost = Hosts.LastOrDefault();
        Changed("Host removed");
        Raise(nameof(RunWhereChoices)); Raise(nameof(HasHosts)); Raise(nameof(RunWhereIndex));
    }

    private string _hostTest = "";
    public string HostTestText { get => _hostTest; private set => Set(ref _hostTest, value); }
    private bool _hostTesting;
    public bool HostIdle { get => !_hostTesting; private set => Set(ref _hostTesting, !value); }

    /// <summary>ssh in batch mode (the agent's keys, no password prompt) and run `caps --version` on the host.</summary>
    public async Task TestHost()
    {
        if (_host?.Host is not { } h || _hostTesting) return;
        if (h.Hostname.Length == 0) { HostTestText = "Enter the hostname first"; return; }
        HostIdle = false;
        HostTestText = "Connecting…";
        var row = _host;
        try
        {
            var psi = new ProcessStartInfo("ssh") { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false };
            foreach (var a in new[] { "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", "-p", h.Port.ToString(CultureInfo.InvariantCulture),
                                      h.User.Length > 0 ? $"{h.User}@{h.Hostname}" : h.Hostname, HostProbeCommand(h) })
                psi.ArgumentList.Add(a);
            using var p = Process.Start(psi) ?? throw new InvalidOperationException("cannot start ssh");
            var outTask = p.StandardOutput.ReadToEndAsync();
            var errTask = p.StandardError.ReadToEndAsync();
            var done = await Task.Run(() => p.WaitForExit(15000));
            if (!done) { try { p.Kill(); } catch { } throw new TimeoutException("no answer in 15 s"); }
            var output = (await outTask).Trim();
            var error = (await errTask).Trim();
            var probe = ParseProbe(output);
            if (probe.TryGetValue("node", out var node) && node.Split(' ', StringSplitOptions.RemoveEmptyEntries) is { Length: >= 1 } nf && int.TryParse(nf[0].TrimEnd('+'), out var cores))
            {
                h.NodeCores = cores;
                h.NodeMem = nf.Length > 1 ? nf[1].TrimEnd('+') : "";
                RaiseClusterFields();
            }
            var version = probe.GetValueOrDefault("version", "none");
            var extras = new List<string>();
            if (h.NodeCores > 0) extras.Add($"{h.NodeCores} cores per node");
            if (probe.ContainsKey("sched:slurm")) extras.Add("SLURM");
            else if (probe.ContainsKey("sched:pbs")) extras.Add("PBS");
            else extras.Add("no queue found");
            if (probe.ContainsKey("ws")) extras.Add("workspaces (ws_allocate)");
            if (p.ExitCode == 0 && version.StartsWith("caps "))
            {
                h.CapsVersion = version;
                var mine = StudioBuild().Commit;
                var theirs = CommitOf(version);
                var same = theirs.Length == 0 || mine == "unknown" || theirs == "unknown" || mine.StartsWith(theirs, StringComparison.Ordinal) || theirs.StartsWith(mine, StringComparison.Ordinal);
                row.State = same ? "connected" : "caps differs";
                row.Level = same ? 1 : 2;
                HostTestText = $"Reachable · {version} on the host · {string.Join(" · ", extras)}" +
                               (same ? "" : $". The Studio is {mine}: Install / update caps builds the same one there");
            }
            else if (p.ExitCode == 0)
            {
                row.State = "no caps";
                row.Level = 2;
                HostTestText = $"Reachable · {string.Join(" · ", extras)} · caps is not installed there yet: Install / update caps builds it in {h.Root}/bin" +
                               (probe.ContainsKey("git") && probe.ContainsKey("cmake") ? "" : " (it needs git and cmake: load them with the build modules)");
            }
            else
            {
                row.State = "unreachable";
                row.Level = 2;
                HostTestText = "ssh failed: " + (error.Length > 0 ? error.Split('\n')[0] : $"exit {p.ExitCode}") + " (keys come from your SSH agent; no passwords are asked or stored)";
            }
        }
        catch (Exception e)
        {
            row.State = "unreachable";
            row.Level = 2;
            HostTestText = e.Message;
        }
        finally
        {
            HostIdle = true;
            Raise(nameof(HostsStatus));
            if (row.State == "unreachable")
                Notify(new Notice
                {
                    Key = "host." + h.Name, Severity = "warning", Icon = "server", Title = $"{(h.Name.Length > 0 ? h.Name : h.Hostname)} unreachable",
                    Body = $"{HostTestText}. Local runs are unaffected; jobs already on the host keep running there.",
                    Primary = "Retry now", OnPrimary = () => _ = TestHost(),
                    Secondary = "Details", OnSecondary = () => { SettingsTab = 3; SetModule(10); },
                });
            else Dismiss("host." + h.Name);
        }
    }

    // ---------------------------------------------------------------- remote jobs (design/boards/RemoteCompute, Jobs)
    // A Dynamics or Equilibrate run sent to a host: the structure and a recipe of the page's settings go up with scp, the
    // job template (filled) is submitted to the host's scheduler over ssh, its state is polled, and the result comes back
    // with its provenance. Keys stay in the SSH agent; nothing is stored but the host's name and the job's folder.
    public List<string> RunWhereChoices => new[] { "This machine" }.Concat(_settings.Hosts.Select(h => h.Name.Length > 0 ? h.Name : h.Hostname)).ToList();
    private int _runWhere;
    public int RunWhereIndex { get => Math.Min(_runWhere, _settings.Hosts.Count); set { if (Set(ref _runWhere, Math.Clamp(value, 0, _settings.Hosts.Count))) { Raise(nameof(RunsRemote)); _rCpus = 0; _rMem = _rTime = ""; _rKeep = null; RaiseRemoteRun(); } } }
    public bool RunsRemote => RunWhereIndex > 0;
    public bool HasHosts => _settings.Hosts.Count > 0;
    public static string RemoteFolder => AppSettings.Override != null ? Path.Combine(Path.GetDirectoryName(AppSettings.Override)!, "caps-remote") : Path.Combine(AppSettings.Folder, "remote");
    private Avalonia.Threading.DispatcherTimer? _remotePoll;

    /// <summary>ssh / scp with the agent's keys (batch mode: no prompts); exit code, stdout and stderr.</summary>
    private static async Task<(int Code, string Out, string Err)> Tool(string exe, IEnumerable<string> args, int timeoutMs = 60000)
    {
        var psi = new ProcessStartInfo(exe) { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false };
        foreach (var a in args) psi.ArgumentList.Add(a);
        using var p = Process.Start(psi) ?? throw new InvalidOperationException("cannot start " + exe);
        var o = p.StandardOutput.ReadToEndAsync();
        var e = p.StandardError.ReadToEndAsync();
        if (!await Task.Run(() => p.WaitForExit(timeoutMs))) { try { p.Kill(true); } catch { } throw new TimeoutException($"{exe}: no answer in {timeoutMs / 1000} s"); }
        return (p.ExitCode, (await o).Trim(), (await e).Trim());
    }
    private static string Target(RemoteHost h) => h.User.Length > 0 ? $"{h.User}@{h.Hostname}" : h.Hostname;
    private static string[] SshArgs(RemoteHost h, string command) =>
        ["-o", "BatchMode=yes", "-o", "ConnectTimeout=10", "-p", h.Port.ToString(CultureInfo.InvariantCulture), Target(h), command];
    private static IEnumerable<string> ScpArgs(RemoteHost h, IEnumerable<string> from, string to, bool recursive = false) =>
        new[] { "-B", "-q", "-P", h.Port.ToString(CultureInfo.InvariantCulture) }.Concat(recursive ? ["-r"] : Array.Empty<string>()).Concat(from).Append(to);
    private static string Q(string s) => "'" + s.Replace("'", "'\\''") + "'";   // a single-quoted shell word

    /// <summary>The page's run as a recipe for caps run on the host: the structure from the file, typed with the Field's
    /// force field, the run, and the result exported.</summary>
    private JsonObject RemoteRecipe(string kind, string stem)
    {
        var inv = CultureInfo.InvariantCulture;
        var r = new JsonObject { ["recipe"] = 1, ["name"] = stem, ["build"] = new JsonObject { ["file"] = "structure.caps.data" }, ["cutoff"] = _relaxCutoff, ["seed"] = _mdSeed };
        // the force field as assigned here, read whole on the host (groups, water models, every parameter kept); React
        // changes the structure, so it types again from the library
        var ff = Field.Assigned && Field.FfIndex >= 0 && Field.FfIndex < Field.Library.Count ? Field.Library[Field.FfIndex].Id : "default";
        r["type"] = Field.Assigned && kind != "React" ? new JsonObject { ["file"] = "structure.ff.json" } : new JsonObject { ["forcefield"] = ff };
        string[] th = ["bussi", "langevin", "nose-hoover"], ba = ["crescale", "berendsen", "mtk"], cons = ["none", "h-bonds", "all-bonds"], solver = ["shake", "lincs"];
        var common = new JsonObject
        {
            ["thermostat"] = th[Math.Clamp(_mdThermostat, 0, 2)], ["barostat"] = ba[Math.Clamp(_mdBarostat, 0, 2)], ["tau_t"] = _mdTauT, ["tau_p"] = _mdTauP,
            ["constraints"] = cons[Math.Clamp(_mdConstraints, 0, 2)], ["constraint_solver"] = solver[Math.Clamp(_mdConstraintSolver, 0, 1)], ["seed"] = _mdSeed,
        };
        if (kind == "Relax")
        {
            var rx = new JsonObject
            {
                ["method"] = _relaxMethod switch { 0 => "sd", 1 => "cg", 3 => "fire", _ => "lbfgs" }, ["fmax"] = _relaxFtol, ["max_iterations"] = _relaxIterations,
                ["pushoff"] = _relaxPushoff,
            };
            if (_relaxCompress) rx["target_density"] = _relaxDensity;
            if (_relaxPushoff && _relaxPushoffMd) { rx["pushoff_md_ps"] = (double)_relaxRampPs; rx["pushoff_cap"] = (double)_relaxCap; rx["pushoff_temperature"] = (double)_relaxPushoffT; }
            r["relax"] = rx;
        }
        else if (kind == "React")
        {
            r["react"] = new JsonObject
            {
                ["templates"] = _rxText, ["cycles"] = _rxCycles, ["per_cycle"] = _rxPerCycle, ["target"] = _rxTarget, ["relax"] = _rxRelax,
                ["md_ps"] = _rxMdPs, ["temperature"] = _rxTemp, ["seed"] = _rxSeed,
            };
            if (_rxCapture > 0) ((JsonObject)r["react"]!)["capture"] = _rxCapture;
        }
        else if (kind == "Dynamics")
        {
            var md = (JsonObject)common.DeepClone();
            md["ensemble"] = new[] { "nve", "nvt", "npt", "nph" }[Math.Clamp(_mdEnsemble, 0, 3)];
            md["ps"] = _mdSteps * _mdDt / 1000.0;
            md["dt"] = _mdDt;
            md["temperature"] = _mdTemp;
            md["pressure"] = _mdPressure;
            if (RespaSteps > 1) md["respa"] = RespaSteps;
            r["md"] = md;
        }
        else
        {
            var eq = (JsonObject)common.DeepClone();
            eq["protocol_text"] = _eqText;
            eq["until_converged"] = _eqUntil;
            eq["block_ps"] = _eqBlock;
            eq["max_blocks"] = _eqMaxBlocks;
            r["equilibrate"] = eq;
        }
        r["export"] = new JsonArray("lammps", "pdb");
        _ = inv;
        return r;
    }

    /// <summary>Sends the open structure and the page's run (Dynamics or Equilibrate) to the chosen host.</summary>
    public async Task SubmitRemote(string kind)
    {
        var grow = kind == "Grow";   // a new cell from the Grow recipe: no structure goes up
        if ((_doc == null && !grow) || RunWhereIndex == 0) return;
        var h = _settings.Hosts[RunWhereIndex - 1];
        if (h.Hostname.Length == 0) { Status = "The host has no hostname: set it in Settings › Compute & remote"; return; }
        var k = _jobCounters[kind] = _jobCounters.GetValueOrDefault(kind) + 1;
        var id = $"{kind.ToLowerInvariant()}-{k}";
        var growRecipe = grow ? GrowRecipe() : "";
        var stem = grow ? growRecipe.Split('\n').FirstOrDefault(l => l.StartsWith("name:"))?[5..].Trim().Trim('"') ?? "cell"
                 : string.Concat(Title.Replace(" (unsaved)", "").Where(char.IsLetterOrDigit).Take(24)) is { Length: > 0 } t ? t : "structure";
        var job = new Job
        {
            Id = id, Kind = kind, Module = kind switch { "Dynamics" => 3, "Equilibrate" => 4, "Relax" => 2, "React" => 6, _ => 0 }, Title = $"{kind} · on {h.Name}",
            Document = grow ? "new cell" : Title, Atoms = grow ? 0 : _doc!.Summary().Atoms, Provenance = Manifest(kind),
        };
        job.Status = "queued";
        job.Remote = new RemoteRun { Host = h.Name, Scheduler = h.Scheduler, Local = Path.Combine(RemoteFolder, $"{id}-{DateTime.Now:yyyyMMdd-HHmmss}"), Stem = stem };
        Jobs.Insert(0, job);
        SelectedJob = job;
        Raise(nameof(HasJobs)); Raise(nameof(JobsSummary));
        try
        {
            var local = job.Remote.Local;
            Directory.CreateDirectory(local);
            var recipeFile = grow ? "recipe.yaml" : "recipe.json";
            if (grow) File.WriteAllText(Path.Combine(local, recipeFile), growRecipe);
            else
            {
                _doc!.Save(Path.Combine(local, "structure.caps.data"));
                var recipe = RemoteRecipe(kind, stem);
                AddLiveOutput(recipe);
                File.WriteAllText(Path.Combine(local, recipeFile), recipe.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
            }
            if (!_relaxCoulomb && !grow) job.Add("Note: recipes always include Coulomb terms; the host's run has them although the page has them off");
            var files = grow ? new List<string> { recipeFile } : new List<string> { "structure.caps.data", recipeFile };
            if (!grow && Field.Assigned && kind != "React") files.AddRange(SaveForceFieldForHost(local));
            await SendJob(job, h, kind, files.ToArray(), $"caps run {recipeFile} --out .");
            Status = $"{kind} sent to {h.Name} · {(h.Scheduler == "none" ? "process" : h.Scheduler)} {job.Remote.JobId} · Jobs follows it";
            StartRemotePoll();
        }
        catch (Exception e)
        {
            job.Status = "failed";
            job.Error = e.Message;
            job.Add("Could not send the job: " + e.Message);
            job.Suggestion = "Test the host in Settings › Compute & remote (keys come from your SSH agent; Install / update caps puts caps on the host).";
            job.SuggestModule = 10;
            job.Ended = DateTime.Now;
            Status = "Could not send the job: " + e.Message;
        }
        SaveJobs();
        Raise(nameof(JobsSummary)); Raise(nameof(ComputeText));
    }

    /// <summary>The assigned force field written whole for a host (structure.ff.json), with a many-body potential's file
    /// beside it (its path made local); the files to upload.</summary>
    private List<string> SaveForceFieldForHost(string local)
    {
        var files = new List<string> { "structure.ff.json" };
        var path = Path.Combine(local, "structure.ff.json");
        _doc!.FieldSave(path);
        var j = JsonNode.Parse(File.ReadAllText(path))!.AsObject();
        if (j["manybody"] is JsonObject mb)
            foreach (var key in new[] { "file", "file2" })
                if ((string?)mb[key] is { Length: > 0 } f && File.Exists(f))
                {
                    var name = Path.GetFileName(f);
                    File.Copy(f, Path.Combine(local, name), true);
                    mb[key] = name;
                    files.Add(name);
                }
        File.WriteAllText(path, j.ToJsonString());
        return files;
    }

    /// <summary>The job folder prepared locally (the recipe, a structure) goes up to the host with job.sh and is submitted.</summary>
    private async Task SendPrepared(Job job, RemoteHost h, string id, string[] files, string recipeFile)
    {
        var local = job.Remote!.Local;
        var script = _settings.JobTemplate.Replace("{job}", id).Replace("{partition}", h.Partition).Replace("{recipe}", recipeFile + " --out out");
        job.Add($"Prepared {local}: {string.Join(", ", files.Where(f => f != "job.sh"))}, job.sh");
        // the host's folder (the template's {workdir} is expanded there: $USER, ~)
        var mk = await Tool("ssh", SshArgs(h, $"mkdir -p \"{h.WorkDir}/{id}\" && cd \"{h.WorkDir}/{id}\" && pwd"), 30000);
        if (mk.Code != 0) throw new InvalidOperationException("ssh: " + (mk.Err.Length > 0 ? mk.Err.Split('\n')[0] : $"exit {mk.Code}"));
        var dir = mk.Out.Split('\n').Last().Trim();
        job.Remote.Dir = dir;
        script = script.Replace("{workdir}/" + id, dir).Replace("{workdir}", Path.GetDirectoryName(dir.Replace('\\', '/'))?.Replace('\\', '/') ?? dir);
        if (h.Scheduler == "none" && !script.Contains("caps run")) script += $"\ncaps run {recipeFile} --out out\n";
        File.WriteAllText(Path.Combine(local, "job.sh"), script.Replace("\r\n", "\n"));
        var up = await Tool("scp", ScpArgs(h, files.Select(f => Path.Combine(local, f)), $"{Target(h)}:{dir}/"), 120000);
        if (up.Code != 0) throw new InvalidOperationException("scp: " + (up.Err.Length > 0 ? up.Err.Split('\n')[0] : $"exit {up.Code}"));
        job.Add($"Uploaded to {h.Name}:{dir}");
        var submit = h.Scheduler switch
        {
            "SLURM" => $"cd {Q(dir)} && sbatch --parsable job.sh",
            "PBS" => $"cd {Q(dir)} && qsub job.sh",
            _ => $"cd {Q(dir)} && (nohup bash job.sh > caps.log 2>&1 & echo $!)",
        };
        var sub = await Tool("ssh", SshArgs(h, submit), 60000);
        if (sub.Code != 0 || sub.Out.Length == 0) throw new InvalidOperationException("submit: " + (sub.Err.Length > 0 ? sub.Err.Split('\n')[0] : $"exit {sub.Code}"));
        job.Remote.JobId = sub.Out.Split('\n').Last().Split(';')[0].Trim();
        job.Status = "running";
        job.Add($"Submitted · {(h.Scheduler == "none" ? "process" : h.Scheduler + " job")} {job.Remote.JobId}");
    }

    /// <summary>A recipe (caps run) sent to a host as its own job: the sweep's cells on a cluster.</summary>
    internal async Task<Job> SendRecipe(RemoteHost h, string kind, string title, string stem, string recipeJson)
    {
        var k = _jobCounters[kind] = _jobCounters.GetValueOrDefault(kind) + 1;
        var id = $"{kind.ToLowerInvariant()}-{k}";
        var job = new Job { Id = id, Kind = kind, Module = 0, Title = $"{title} · on {h.Name}", Document = "new cell", Atoms = 0, Provenance = Manifest(kind) };
        job.Status = "queued";
        job.Remote = new RemoteRun { Host = h.Name, Scheduler = h.Scheduler, Local = Path.Combine(RemoteFolder, $"{id}-{DateTime.Now:yyyyMMdd-HHmmss}"), Stem = stem };
        Jobs.Insert(0, job);
        Raise(nameof(HasJobs)); Raise(nameof(JobsSummary));
        try
        {
            Directory.CreateDirectory(job.Remote.Local);
            File.WriteAllText(Path.Combine(job.Remote.Local, "recipe.json"), recipeJson);
            var keep = _rTitle;
            _rTitle = stem;
            try { await SendJob(job, h, kind, ["recipe.json"], "caps run recipe.json --out ."); }
            finally { _rTitle = keep; }
            StartRemotePoll();
        }
        catch (Exception e)
        {
            job.Status = "failed";
            job.Error = e.Message;
            job.Add("Could not send the job: " + e.Message);
            job.Suggestion = "Test the host in Settings › Compute & remote (keys come from your SSH agent; Install / update caps puts caps on the host).";
            job.SuggestModule = 10;
            job.Ended = DateTime.Now;
        }
        SaveJobs();
        Raise(nameof(JobsSummary)); Raise(nameof(ComputeText));
        return job;
    }

    private void StartRemotePoll()
    {
        if (_remotePoll != null) return;
        // every 30 s for the job shown in Jobs, every 5 min for the others (PollDue); the earlier plain scripts every minute
        _remotePoll = new Avalonia.Threading.DispatcherTimer { Interval = TimeSpan.FromSeconds(30) };
        _remotePoll.Tick += (_, _) =>
        {
            foreach (var j in Jobs.Where(j => j.IsRemote && (j.IsRunning || j.IsQueued)).ToList())
                if (j.Remote!.Mode.Length == 0 ? (DateTime.Now - j.Remote.LastPoll).TotalSeconds >= 55 : PollDue(j))
                {
                    if (j.Remote.Mode.Length == 0) j.Remote.LastPoll = DateTime.Now;
                    _ = RefreshJob(j);
                }
        };
        _remotePoll.Start();
    }

    /// <summary>Asks the host where the job stands; when it has ended, brings the result and its log back.</summary>
    public async Task CheckRemote(Job? j)
    {
        if (j?.Remote is not { } r || !(j.IsRunning || j.IsQueued) || r.Checking) return;
        var h = _settings.Hosts.FirstOrDefault(x => x.Name == r.Host);
        if (h == null) { j.Add($"The host {r.Host} is no longer in Settings"); return; }
        if (r.Mode == "job") { await PollJob(j, h, r); return; }
        if (r.Mode == "install") { await PollInstall(j, h, r); return; }
        r.Checking = true;
        try
        {
            var probe = r.Scheduler switch
            {
                "SLURM" => $"squeue -h -j {Q(r.JobId)} -o %T 2>/dev/null",
                "PBS" => $"qstat {Q(r.JobId)} >/dev/null 2>&1 && echo RUNNING",
                _ => $"kill -0 {Q(r.JobId)} 2>/dev/null && echo RUNNING",
            };
            var st = await Tool("ssh", SshArgs(h, probe), 30000);
            var state = st.Out.Trim();
            if (state.Length > 0)
            {
                if (state != r.LastState) j.Add($"{r.Host}: {state.ToLowerInvariant()}");
                r.LastState = state;
                return;
            }
            await BringBack(j, h, r);
        }
        catch (Exception e) { j.Add($"{r.Host} not reached: {e.Message} (the job keeps running there; checked again in a minute)"); }
        finally { r.Checking = false; }
    }

    /// <summary>A job that has ended: the result, the log and the small outputs come back; a large trajectory stays on
    /// the host until asked for.</summary>
    private async Task BringBack(Job j, RemoteHost h, RemoteRun r, bool finished = false)
    {
        {
            var outDir = Path.Combine(r.Local, "out");
            Directory.CreateDirectory(outDir);
            var ls = await Tool("ssh", SshArgs(h, ListOutCommand(r.Dir)), 60000);
            var (now, later) = CopyPlan(ParseListing(ls.Out), r.Stem, LargeRemoteBytes);
            var back = (Code: 0, Out: "", Err: "");
            if (now.Count > 0)
                back = await Tool("scp", ScpArgs(h, now.Select(f => $"{Target(h)}:{r.Dir}/out/{f.Name}"), outDir + "/"), CopyTimeoutMs(now.Sum(f => f.Bytes)));
            await Tool("scp", ScpArgs(h, [$"{Target(h)}:{r.Dir}/{(r.Mode == "job" ? "slurm-*" : "*.log")}"], r.Local), 60000);
            var result = Path.Combine(outDir, r.Stem + ".data");
            j.Ended = DateTime.Now;
            // a structure, or (a Glass replica: the scan runs on a copy) its properties
            if (back.Code == 0 && (File.Exists(result) || File.Exists(Path.Combine(outDir, r.Stem + ".properties.json")) || finished))
            {
                j.Status = "done";
                j.Progress = 1;
                r.OnHost = later;
                j.Add($"Finished on {r.Host}; {now.Count} file{(now.Count == 1 ? "" : "s")} brought back to {outDir}");
                foreach (var f in later) j.Add($"{f.Name} ({f.Size}) left on {r.Host}: copy it back from the outputs, whole or every 10th / 100th frame");
                RemoteOutputs(j);
                if (r.Batch.Length > 0 && Jobs.Where(x => x.Remote?.Batch == r.Batch).All(x => !x.IsRunning)) CollectGlassReplicas(r.Batch);
            }
            else
            {
                var log = File.Exists(Path.Combine(outDir, "run.log")) ? File.ReadAllText(Path.Combine(outDir, "run.log"))
                        : Directory.Exists(r.Local) ? Directory.GetFiles(r.Local, "*.log").Select(File.ReadAllText).FirstOrDefault() ?? "" : "";
                j.Status = "failed";
                j.Error = "No result came back" + (log.Length > 0 ? ":\n" + string.Join("\n", log.Split('\n').TakeLast(8)) : $" ({(back.Err.Length > 0 ? back.Err.Split('\n')[0] : ls.Code != 0 ? "no out folder on the host" : "no " + r.Stem + ".data in out")})");
                j.Add(j.Error);
            }
            SaveJobs();
            Raise(nameof(JobsSummary)); Raise(nameof(ComputeText));
        }
    }

    // ---------------------------------------------------------------- copying results back
    /// <summary>Outputs above this size stay on the host when a job ends (copied on request, whole or thinned).</summary>
    public const long LargeRemoteBytes = 200L << 20;
    private static readonly string[] TrajectoryExts = [".lammpstrj", ".dump", ".dcd", ".xtc", ".trr", ".nc", ".mdcrd"];

    /// <summary>The out folder's files with their sizes, one "bytes TAB name" per line (POSIX sh on the host).</summary>
    internal static string ListOutCommand(string dir) =>
        $"cd {Q(dir + "/out")} && for f in *; do [ -f \"$f\" ] && printf '%s\\t%s\\n' \"$(wc -c < \"$f\" | tr -d ' ')\" \"$f\"; done";

    internal static List<RemoteFile> ParseListing(string text) =>
        text.Split('\n').Select(l => l.TrimEnd('\r').Split('\t', 2)).Where(p => p.Length == 2 && long.TryParse(p[0].Trim(), out _))
            .Select(p => new RemoteFile(p[1], long.Parse(p[0].Trim(), CultureInfo.InvariantCulture))).ToList();

    /// <summary>Which files come back when the job ends (the result always; others up to large) and which wait.</summary>
    internal static (List<RemoteFile> Now, List<RemoteFile> Later) CopyPlan(IEnumerable<RemoteFile> files, string stem, long large)
    {
        var now = new List<RemoteFile>();
        var later = new List<RemoteFile>();
        foreach (var f in files) (f.Name == stem + ".data" || f.Bytes <= large ? now : later).Add(f);
        return (now, later);
    }

    /// <summary>scp / rsync time allowed: two minutes, or the size at 1 MB/s when longer.</summary>
    internal static int CopyTimeoutMs(long bytes) => (int)Math.Min(int.MaxValue, Math.Max(120000L, bytes / 1048576L * 1000L));

    internal static bool IsTrajectoryFile(string name) => TrajectoryExts.Contains(Path.GetExtension(name).ToLowerInvariant());

    /// <summary>The thinned copy's name: every Nth frame, in a format caps frames writes (.xtc as .trr, AMBER as .dcd).</summary>
    internal static string ThinnedName(string name, int every)
    {
        var ext = Path.GetExtension(name).ToLowerInvariant();
        var outExt = ext switch { ".xtc" => ".trr", ".nc" or ".mdcrd" => ".dcd", _ => ext };
        return Path.GetFileNameWithoutExtension(name) + $".every{every}" + outExt;
    }

    /// <summary>Thins a trajectory on the host with caps frames (the frames not kept are passed over as it is read).</summary>
    internal static string ThinCommand(string dir, string name, int every, string? topology) =>
        $"cd {Q(dir + "/out")} && caps frames {Q(name)} {Q(ThinnedName(name, every))} --stride {every}" + (topology != null ? $" --topology {Q(topology)}" : "");

    /// <summary>The job's outputs: the result brought back, trajectories copied back, and actions for files left on the host.</summary>
    private void RemoteOutputs(Job j)
    {
        if (j.Remote is not { } r) return;
        j.Outputs.Clear();
        var outDir = Path.Combine(r.Local, "out");
        var result = Path.Combine(outDir, r.Stem + ".data");
        if (File.Exists(result)) j.Outputs.Add(new JobOutput($"Result · {r.Stem}.data", "cube", () => Open(result)));
        if (r.Batch.Length > 0 && File.Exists(Path.Combine(outDir, r.Stem + ".properties.json")))
            j.Outputs.Add(new JobOutput("Glass transition · pool the replicas back", "chart", () => CollectGlassReplicas(r.Batch)));
        if (Directory.Exists(outDir))
            foreach (var t in Directory.GetFiles(outDir).Where(f => IsTrajectoryFile(f)).OrderBy(f => f))
            {
                var path = t;
                j.Outputs.Add(new JobOutput($"Trajectory · {Path.GetFileName(t)}", "play", () => Open(path, File.Exists(result) ? result : null)));
            }
        foreach (var f in r.OnHost)
        {
            var file = f;
            j.Outputs.Add(new JobOutput($"Copy {f.Name} back · {f.Size}", "download", () => _ = CopyBack(j, file, 0)));
            if (IsTrajectoryFile(f.Name))
            {
                j.Outputs.Add(new JobOutput($"Copy every 10th frame of {f.Name}", "download", () => _ = CopyBack(j, file, 10)));
                j.Outputs.Add(new JobOutput($"Copy every 100th frame of {f.Name}", "download", () => _ = CopyBack(j, file, 100)));
            }
        }
    }

    /// <summary>Brings a file left on the host back: whole (rsync --partial, so a broken copy resumes when asked again;
    /// scp when rsync is missing) or, for a trajectory, every Nth frame (thinned on the host by caps frames, then copied).</summary>
    public async Task CopyBack(Job j, RemoteFile f, int every)
    {
        if (j.Remote is not { } r || r.Checking) return;
        var h = _settings.Hosts.FirstOrDefault(x => x.Name == r.Host);
        if (h == null) { j.Add($"The host {r.Host} is no longer in Settings"); return; }
        var outDir = Path.Combine(r.Local, "out");
        Directory.CreateDirectory(outDir);
        r.Checking = true;
        try
        {
            var name = f.Name;
            long bytes = f.Bytes;
            if (every > 1)
            {
                j.Add($"Thinning {f.Name} on {r.Host}: every {every}th frame…");
                var top = File.Exists(Path.Combine(outDir, r.Stem + ".data")) ? r.Stem + ".data" : null;
                var thin = await Tool("ssh", SshArgs(h, ThinCommand(r.Dir, f.Name, every, top)), Math.Max(600000, CopyTimeoutMs(f.Bytes) / 4));
                if (thin.Code != 0)
                {
                    var why = (thin.Err + "\n" + thin.Out).Contains("usage", StringComparison.OrdinalIgnoreCase) ? "the host's caps has no frames command: update caps there" : thin.Err.Split('\n')[0];
                    j.Add($"Could not thin {f.Name}: {why}");
                    return;
                }
                name = ThinnedName(f.Name, every);
                bytes = f.Bytes / every;
                j.Add($"{r.Host}: {thin.Out.Split('\n').LastOrDefault()}");
            }
            j.Add($"Copying {name} back (~{new RemoteFile(name, bytes).Size})…");
            var src = $"{Target(h)}:{r.Dir}/out/{name}";
            var port = h.Port.ToString(CultureInfo.InvariantCulture);
            (int Code, string Out, string Err) copy;
            try { copy = await Tool("rsync", ["-t", "--partial", "-e", $"ssh -p {port} -o BatchMode=yes", src, outDir + "/"], CopyTimeoutMs(bytes)); }
            catch (Exception e) when (e is not TimeoutException) { copy = (127, "", e.Message); }
            if (copy.Code != 0) copy = await Tool("scp", ScpArgs(h, [src], outDir + "/"), CopyTimeoutMs(bytes));
            if (copy.Code != 0 || !File.Exists(Path.Combine(outDir, name)))
            {
                j.Add($"Could not copy {name}: {(copy.Err.Length > 0 ? copy.Err.Split('\n')[0] : $"exit {copy.Code}")} (ask again to resume)");
                return;
            }
            if (every <= 1) r.OnHost.RemoveAll(x => x.Name == f.Name);
            j.Add($"{name} is in {outDir}");
            RemoteOutputs(j);
            SaveJobs();
        }
        catch (Exception e) { j.Add($"{r.Host} not reached: {e.Message}"); }
        finally { r.Checking = false; }
    }

    /// <summary>Opens the structure a remote job brought back.</summary>
    public void OpenRemoteResult(Job? j)
    {
        if (j?.Remote is not { } r) return;
        var result = Path.Combine(r.Local, "out", r.Stem + ".data");
        if (!File.Exists(result)) { Status = "No result yet"; return; }
        Open(result);
    }

    /// <summary>The batch-job template ({job}, {partition}, {workdir}, {recipe} are filled per job).</summary>
    public string JobTemplate { get => _settings.JobTemplate; set { if (value != null && value != _settings.JobTemplate) { _settings.JobTemplate = value; Raise(); Changed("Job template"); } } }
    public void ResetJobTemplate() { JobTemplate = RemoteHost.DefaultTemplate; }
}
