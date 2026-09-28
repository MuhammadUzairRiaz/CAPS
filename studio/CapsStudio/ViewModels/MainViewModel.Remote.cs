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
                                      h.User.Length > 0 ? $"{h.User}@{h.Hostname}" : h.Hostname, "caps --version || echo caps-not-found" })
                psi.ArgumentList.Add(a);
            using var p = Process.Start(psi) ?? throw new InvalidOperationException("cannot start ssh");
            var outTask = p.StandardOutput.ReadToEndAsync();
            var errTask = p.StandardError.ReadToEndAsync();
            var done = await Task.Run(() => p.WaitForExit(15000));
            if (!done) { try { p.Kill(); } catch { } throw new TimeoutException("no answer in 15 s"); }
            var output = (await outTask).Trim();
            var error = (await errTask).Trim();
            if (p.ExitCode == 0 && output.StartsWith("caps "))
            {
                row.State = "connected";
                row.Level = 1;
                HostTestText = $"Reachable · {output} on the host";
            }
            else if (p.ExitCode == 0)
            {
                row.State = "no caps";
                row.Level = 2;
                HostTestText = "Reachable, but caps is not on the host's PATH: install the CAPS command-line tools there";
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
    public int RunWhereIndex { get => Math.Min(_runWhere, _settings.Hosts.Count); set { if (Set(ref _runWhere, Math.Clamp(value, 0, _settings.Hosts.Count))) Raise(nameof(RunsRemote)); } }
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
        var r = new JsonObject { ["recipe"] = 1, ["name"] = stem, ["build"] = new JsonObject { ["file"] = "structure.data" }, ["cutoff"] = _relaxCutoff, ["seed"] = _mdSeed };
        var ff = Field.Assigned && Field.FfIndex >= 0 && Field.FfIndex < Field.Library.Count ? Field.Library[Field.FfIndex].Id : "default";
        r["type"] = new JsonObject { ["forcefield"] = ff };
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
                _doc!.Save(Path.Combine(local, "structure.data"));
                File.WriteAllText(Path.Combine(local, recipeFile), RemoteRecipe(kind, stem).ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
            }
            if (!_relaxCoulomb && !grow) job.Add("Note: recipes always include Coulomb terms; the host's run has them although the page has them off");
            var script = _settings.JobTemplate.Replace("{job}", id).Replace("{partition}", h.Partition).Replace("{recipe}", recipeFile + " --out out");
            job.Add($"Prepared {local}: {(grow ? "" : "structure.data, ")}{recipeFile}, job.sh");
            // the host's folder (the template's {workdir} is expanded there: $USER, ~)
            var mk = await Tool("ssh", SshArgs(h, $"mkdir -p \"{h.WorkDir}/{id}\" && cd \"{h.WorkDir}/{id}\" && pwd"), 30000);
            if (mk.Code != 0) throw new InvalidOperationException("ssh: " + (mk.Err.Length > 0 ? mk.Err.Split('\n')[0] : $"exit {mk.Code}"));
            var dir = mk.Out.Split('\n').Last().Trim();
            job.Remote.Dir = dir;
            script = script.Replace("{workdir}/" + id, dir).Replace("{workdir}", Path.GetDirectoryName(dir.Replace('\\', '/'))?.Replace('\\', '/') ?? dir);
            if (h.Scheduler == "none" && !script.Contains("caps run")) script += $"\ncaps run {recipeFile} --out out\n";
            File.WriteAllText(Path.Combine(local, "job.sh"), script.Replace("\r\n", "\n"));
            var files = grow ? new[] { recipeFile, "job.sh" } : new[] { "structure.data", recipeFile, "job.sh" };
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
            Status = $"{kind} sent to {h.Name} · {h.Scheduler} {job.Remote.JobId} · Jobs follows it";
            StartRemotePoll();
        }
        catch (Exception e)
        {
            job.Status = "failed";
            job.Error = e.Message;
            job.Add("Could not send the job: " + e.Message);
            job.Suggestion = "Test the host in Settings › Compute & remote (keys come from your SSH agent; caps must be on the host's PATH).";
            job.SuggestModule = 10;
            job.Ended = DateTime.Now;
            Status = "Could not send the job: " + e.Message;
        }
        SaveJobs();
        Raise(nameof(JobsSummary)); Raise(nameof(ComputeText));
    }

    private void StartRemotePoll()
    {
        if (_remotePoll != null) return;
        _remotePoll = new Avalonia.Threading.DispatcherTimer { Interval = TimeSpan.FromSeconds(60) };
        _remotePoll.Tick += (_, _) => { foreach (var j in Jobs.Where(j => j.IsRemote && j.IsRunning).ToList()) _ = CheckRemote(j); };
        _remotePoll.Start();
    }

    /// <summary>Asks the host where the job stands; when it has ended, brings the result and its log back.</summary>
    public async Task CheckRemote(Job? j)
    {
        if (j?.Remote is not { } r || !j.IsRunning || r.Checking) return;
        var h = _settings.Hosts.FirstOrDefault(x => x.Name == r.Host);
        if (h == null) { j.Add($"The host {r.Host} is no longer in Settings"); return; }
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
            // ended: the results, the log and the recipe's provenance come back
            Directory.CreateDirectory(r.Local);
            var back = await Tool("scp", ScpArgs(h, [$"{Target(h)}:{r.Dir}/out"], r.Local, true), 600000);
            await Tool("scp", ScpArgs(h, [$"{Target(h)}:{r.Dir}/*.log"], r.Local), 60000);
            var result = Path.Combine(r.Local, "out", r.Stem + ".data");
            j.Ended = DateTime.Now;
            if (back.Code == 0 && File.Exists(result))
            {
                j.Status = "done";
                j.Progress = 1;
                j.Add($"Finished on {r.Host}; the result is in {Path.Combine(r.Local, "out")}");
                j.Outputs.Add(new JobOutput($"Result · {r.Stem}.data", "cube", () => Open(result)));
            }
            else
            {
                var log = Directory.Exists(r.Local) ? Directory.GetFiles(r.Local, "*.log").Select(File.ReadAllText).FirstOrDefault() ?? "" : "";
                j.Status = "failed";
                j.Error = "No result came back" + (log.Length > 0 ? ":\n" + string.Join("\n", log.Split('\n').TakeLast(8)) : $" ({(back.Err.Length > 0 ? back.Err.Split('\n')[0] : "no out folder on the host")})");
                j.Add(j.Error);
            }
            SaveJobs();
            Raise(nameof(JobsSummary)); Raise(nameof(ComputeText));
        }
        catch (Exception e) { j.Add($"{r.Host} not reached: {e.Message} (the job keeps running there; checked again in a minute)"); }
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
