using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Where a remote job runs: the host, its scheduler and job id, the folder there and the one here.</summary>
public sealed class RemoteRun
{
    public string Host { get; set; } = "";
    public string Scheduler { get; set; } = "SLURM";
    public string JobId { get; set; } = "";
    public string Dir { get; set; } = "";
    public string Local { get; set; } = "";
    public string Stem { get; set; } = "structure";
    public string LastState { get; set; } = "";
    /// <summary>Jobs sent together (a Glass scan's replicas): pooled when they are back.</summary>
    public string Batch { get; set; } = "";
    public bool Checking { get; set; }
    /// <summary>Outputs left on the host when the job came back (large trajectories): copied on request.</summary>
    public List<RemoteFile> OnHost { get; set; } = new();
    public JsonObject Json() => new()
    {
        ["host"] = Host, ["scheduler"] = Scheduler, ["job_id"] = JobId, ["dir"] = Dir, ["local"] = Local, ["stem"] = Stem, ["batch"] = Batch,
        ["on_host"] = new JsonArray(OnHost.Select(f => (JsonNode)new JsonObject { ["name"] = f.Name, ["bytes"] = f.Bytes }).ToArray()),
    };
    public static RemoteRun From(JsonObject o) => new()
    {
        Host = (string?)o["host"] ?? "", Scheduler = (string?)o["scheduler"] ?? "SLURM", JobId = (string?)o["job_id"] ?? "", Dir = (string?)o["dir"] ?? "",
        Local = (string?)o["local"] ?? "", Stem = (string?)o["stem"] ?? "structure", Batch = (string?)o["batch"] ?? "",
        OnHost = o["on_host"] is JsonArray a ? a.OfType<JsonObject>().Select(f => new RemoteFile((string?)f["name"] ?? "", (long?)f["bytes"] ?? 0)).Where(f => f.Name.Length > 0).ToList() : new(),
    };
}

/// <summary>A file in a remote job's out folder and its size.</summary>
public sealed record RemoteFile(string Name, long Bytes)
{
    public string Size => Bytes >= 1L << 30 ? (Bytes / 1073741824.0).ToString("F1", CultureInfo.InvariantCulture) + " GB"
                        : (Bytes / 1048576.0).ToString("F1", CultureInfo.InvariantCulture) + " MB";
}

/// <summary>One line of a job's log.</summary>
public sealed record JobLine(string Time, string Text);

/// <summary>A key and value of a job's provenance manifest.</summary>
public sealed record JobFact(string Key, string Value);

/// <summary>A run (Grow, Relax, Dynamics, Equilibrate, Pack, React, Analyze): its progress, log, curves and the
/// provenance needed to repeat it. Kept for the session and in ~/.caps/jobs.json (without the curves).</summary>
public sealed class Job : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));

    public string Id { get; init; } = "";
    public string Kind { get; init; } = "";
    /// <summary>The structure whose job folder holds it in the project tree (the session only).</summary>
    public ProjectItem? Item { get; set; }
    /// <summary>What the job left: the result, the trajectory, charts, the report …</summary>
    public ObservableCollection<JobOutput> Outputs { get; } = new();
    private bool _expanded;
    public bool Expanded { get => _expanded; set { _expanded = value; Raise(nameof(Expanded)); } }
    public string FolderLabel => $"{Title} · {Started:HH:mm}";
    public int Module { get; init; }
    public string Title { get; init; } = "";
    public string Document { get; init; } = "";
    public long Atoms { get; init; }
    public DateTime Started { get; set; } = DateTime.Now;
    private DateTime? _ended;
    public DateTime? Ended { get => _ended; set { _ended = value; Raise(nameof(Ended)); Raise(nameof(Duration)); } }

    private string _status = "running";
    /// <summary>queued | running | done | failed | cancelled | stopped (finished without meeting its target)</summary>
    public string Status
    {
        get => _status;
        set
        {
            _status = value;
            foreach (var n in new[] { nameof(Status), nameof(IsRunning), nameof(IsQueued), nameof(IsFailed), nameof(IsDone), nameof(IsQuiet), nameof(StatusText), nameof(ShowPause), nameof(CanCancel),
                                      nameof(CanCheckRemote), nameof(CanOpenRemote), nameof(Where), nameof(HasFix) }) Raise(n);
        }
    }
    public bool IsRunning => _status == "running";
    public bool IsQueued => _status == "queued";
    private bool _paused;
    /// <summary>A running job held at its next report (Pause); nothing is lost.</summary>
    public bool Paused { get => _paused; set { _paused = value; Raise(nameof(Paused)); Raise(nameof(StatusText)); Raise(nameof(ShowPause)); } }
    public bool ShowPause => IsRunning && !_paused && Kind is "Dynamics" or "Equilibrate" or "Relax";
    public bool CanCancel => (IsRunning || IsQueued) && !IsRemote;
    public bool IsFailed => _status == "failed";
    public bool IsDone => _status == "done";
    public bool IsQuiet => _status is "cancelled" or "stopped";
    public string StatusText => _paused && _status == "running" ? "paused" : _status;

    private double _progress;
    public double Progress { get => _progress; set { _progress = Math.Clamp(value, 0, 1); Raise(nameof(Progress)); Raise(nameof(ProgressWidth)); } }
    /// <summary>For the list's thin progress line (0 … 1 of its width, drawn by the view).</summary>
    public double ProgressWidth => _progress;
    private int _stage, _stages;
    public int Stage { get => _stage; set { _stage = value; Raise(nameof(Stage)); Raise(nameof(StageText)); } }
    public int Stages { get => _stages; set { _stages = value; Raise(nameof(Stages)); Raise(nameof(StageText)); } }
    public string StageText => _stages > 0 ? $"stage {_stage} of {_stages}" : "";

    public ObservableCollection<JobLine> Log { get; } = new();
    public List<JobFact> Provenance { get; init; } = new();
    public string Error { get; set; } = "";
    /// <summary>What to do about a failure the Studio recognises, and the page that does it (-1: none).</summary>
    public string Suggestion { get; set; } = "";
    public int SuggestModule { get; set; } = -1;
    public bool HasSuggestion => Suggestion.Length > 0;
    public string SuggestButton => SuggestModule switch { 7 => "Open Field", 2 => "Open Relax", 3 => "Open Dynamics", 5 => "Open Pack", _ => "Open" };
    /// <summary>Restart with fix (design/boards/FailedJob): a change the Studio makes and runs again — md-dt (half the time
    /// step, continue from the checkpoint), rx-gentle (fewer reactions a cycle, longer MD after each, capture ≤ 4 Å),
    /// pack-box (every edge 10 % longer); "" when there is none. FixText says what it changes.</summary>
    public string Fix { get; set; } = "";
    public string FixText { get; set; } = "";
    public bool HasFix => Fix.Length > 0 && Status == "failed";
    public string Subtitle => $"{Document} · {Atoms.ToString("N0", CultureInfo.InvariantCulture)} atoms · started {Started:HH:mm}";
    public string Where => Remote is { } r ? $"{Id} · {r.Host} · {(r.Scheduler == "none" ? "process" : r.Scheduler)} {r.JobId}" : $"{Id} · Local · {Environment.ProcessorCount} threads";
    /// <summary>A job sent to a host (Settings › Compute &amp; remote); null for this machine.</summary>
    public RemoteRun? Remote { get; set; }
    public bool IsRemote => Remote != null;
    public bool CanCheckRemote => IsRemote && IsRunning;
    public bool CanOpenRemote => IsRemote && IsDone;
    public string Duration => Ended is { } e ? Seconds(e - Started) : Seconds(DateTime.Now - Started) + " so far";

    private static string Seconds(TimeSpan t) => t.TotalSeconds < 90 ? $"{t.TotalSeconds:F1} s" : t.TotalMinutes < 90 ? $"{t.TotalMinutes:F1} min" : $"{t.TotalHours:F1} h";

    // curves (two panels): the view reads these when CurvesChanged fires
    public string CurveA { get; set; } = "";
    public string CurveB { get; set; } = "";
    public string AxisA { get; set; } = "";
    public string AxisB { get; set; } = "";
    public string AxisX { get; set; } = "";
    public List<(double X, double Y)> A { get; } = new();
    public List<(double X, double Y)> B { get; } = new();
    public bool HasCurves => A.Count > 0 || B.Count > 0;
    public void CurvesChanged() => Raise(nameof(HasCurves));
    public void Tick() => Raise(nameof(Duration));

    public void Add(string text)
    {
        foreach (var line in text.Split('\n', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
        {
            if (Log.Count > 0 && Log[^1].Text == line) continue;
            Log.Add(new JobLine(DateTime.Now.ToString("HH:mm:ss", CultureInfo.InvariantCulture), line));
            if (Log.Count > 400) Log.RemoveAt(0);
        }
    }

    public JsonObject ToJson(bool full)
    {
        var o = new JsonObject
        {
            ["id"] = Id, ["kind"] = Kind, ["module"] = Module, ["title"] = Title, ["document"] = Document, ["atoms"] = Atoms,
            ["started"] = Started.ToString("o"), ["ended"] = Ended?.ToString("o"), ["status"] = _status, ["error"] = Error,
            ["provenance"] = new JsonObject(Provenance.Select(f => KeyValuePair.Create(f.Key, (JsonNode?)f.Value))),
        };
        o["log"] = new JsonArray(Log.TakeLast(full ? 400 : 60).Select(l => (JsonNode)new JsonObject { ["t"] = l.Time, ["text"] = l.Text }).ToArray());
        if (Remote != null) o["remote"] = Remote.Json();
        return o;
    }

    public static Job FromJson(JsonNode n)
    {
        var j = new Job
        {
            Id = (string?)n["id"] ?? "", Kind = (string?)n["kind"] ?? "", Module = (int?)n["module"] ?? 8, Title = (string?)n["title"] ?? "",
            Document = (string?)n["document"] ?? "", Atoms = (long?)n["atoms"] ?? 0,
            Started = DateTime.TryParse((string?)n["started"], CultureInfo.InvariantCulture, DateTimeStyles.RoundtripKind, out var s) ? s : DateTime.Now,
            Error = (string?)n["error"] ?? "",
        };
        if (DateTime.TryParse((string?)n["ended"], CultureInfo.InvariantCulture, DateTimeStyles.RoundtripKind, out var e)) j.Ended = e;
        var st = (string?)n["status"] ?? "done";
        if (n["remote"] is JsonObject ro) j.Remote = RemoteRun.From(ro);
        j.Status = st == "running" && j.Remote == null ? "stopped" : st;   // a local job running when the Studio closed did not finish; a remote one goes on
        if (n["provenance"] is JsonObject p) foreach (var kv in p) j.Provenance.Add(new JobFact(kv.Key, (string?)kv.Value ?? ""));
        if (n["log"] is JsonArray log) foreach (var l in log) if (l != null) j.Log.Add(new JobLine((string?)l["t"] ?? "", (string?)l["text"] ?? ""));
        j.Progress = st == "done" ? 1 : 0;
        return j;
    }
}

/// <summary>Jobs (design/boards/Jobs, FailedJob): every run the Studio starts, with its log, curves and provenance.</summary>
public sealed partial class MainViewModel
{
    public bool IsJobs => _module == 11;
    public ObservableCollection<Job> Jobs { get; } = new();
    private Job? _job;
    public Job? SelectedJob { get => _job; set { if (Set(ref _job, value)) { Raise(nameof(HasSelectedJob)); JobSelected?.Invoke(); } } }
    public bool HasSelectedJob => _job != null;
    public bool HasJobs => Jobs.Count > 0;
    public event Action? JobSelected;
    public event Action? JobCurvesChanged;
    public string JobsSummary
    {
        get
        {
            var r = Jobs.Count(j => j.IsRunning);
            var f = Jobs.Count(j => j.IsFailed);
            var q = Jobs.Count(j => j.IsQueued);
            return Jobs.Count == 0 ? "no jobs" : $"{r} running{(q > 0 ? $" · {q} queued" : "")} · {Jobs.Count(j => j.IsDone)} done{(f > 0 ? $" · {f} failed" : "")}";
        }
    }

    private readonly Dictionary<string, Job> _live = new();
    private readonly Dictionary<string, int> _jobCounters = new();
    private bool _jobsHooked;

    public static string JobsFile => AppSettings.Override != null ? Path.Combine(Path.GetDirectoryName(AppSettings.Override)!, "caps-jobs.json") : Path.Combine(AppSettings.Folder, "jobs.json");

    /// <summary>Starts watching the runs (once) and reads the saved history.</summary>
    public void HookJobs()
    {
        if (_jobsHooked) return;
        _jobsHooked = true;
        try
        {
            if (File.Exists(JobsFile) && JsonNode.Parse(File.ReadAllText(JobsFile)) is JsonArray arr)
                foreach (var n in arr) if (n != null) Jobs.Add(Job.FromJson(n));
            foreach (var j in Jobs)
                if (int.TryParse(j.Id.Split('-').LastOrDefault(), out var k)) _jobCounters[j.Kind] = Math.Max(_jobCounters.GetValueOrDefault(j.Kind), k);
        }
        catch { /* a broken history is left behind */ }
        RefreshJobsShown();
        foreach (var j in Jobs.Where(j => j.IsRemote && j.Status == "done")) RemoteOutputs(j);   // results brought back, files left on the host
        if (Jobs.Any(j => j.IsRemote && j.IsRunning)) Avalonia.Threading.Dispatcher.UIThread.Post(StartRemotePoll);   // remote jobs sent before
        PropertyChanged += OnRunProperty;
        Analyze.PropertyChanged += (_, e) =>
        {
            if (e.PropertyName == nameof(AnalyzeViewModel.Working)) Track("Analyze", Analyze.Working, 1, "Analyze · properties");
            if (e.PropertyName == nameof(AnalyzeViewModel.Log) && _live.TryGetValue("Analyze", out var j)) j.Add(Analyze.Log);
            if (e.PropertyName is nameof(AnalyzeViewModel.SurfaceMolecules) or nameof(AnalyzeViewModel.AxisIndex)) Raise(nameof(IfSetupText));
        };
        ThermoChanged += OnThermo;
        RelaxCurvesChanged += OnRelaxCurves;
        Raise(nameof(HasJobs));
        Raise(nameof(JobsSummary));
    }

    private void OnRunProperty(object? s, PropertyChangedEventArgs e)
    {
        if (e.PropertyName is nameof(MdRunning) or nameof(EqRunning) or nameof(Relaxing)) Raise(nameof(CanPause));
        if (e.PropertyName == nameof(RecipeRunning) && !RecipeRunning && _live.TryGetValue("Recipe", out var rj)) rj.Add(RecipeLog.Trim() + (Status.Length > 0 ? "\n" + Status : ""));
        switch (e.PropertyName)
        {
            case nameof(Growing): Track("Grow", Growing, 0, "Grow · amorphous cell"); break;
            case nameof(RecipeRunning): Track("Recipe", RecipeRunning, 0, "Recipe"); break;
            case nameof(Relaxing): Track("Relax", Relaxing, 2, "Relax · " + Minimisers[Math.Clamp(RelaxMethod, 0, Minimisers.Length - 1)]); break;
            case nameof(MdRunning): Track("Dynamics", MdRunning, 3, "Dynamics · " + Ensembles[Math.Clamp(MdEnsemble, 0, Ensembles.Length - 1)]); break;
            case nameof(EqRunning): Track("Equilibrate", EqRunning, 4, "Equilibrate · " + Protocols[Math.Clamp(EqProtocol, 0, Protocols.Length - 1)]); break;
            case nameof(Packing): Track("Pack", Packing, 5, "Pack · molecules into a box"); break;
            case nameof(Reacting): Track("React", Reacting, 6, "React · crosslinking"); break;
            case nameof(BenchRunning): Track("Bench", BenchRunning, 12, "Bench · validation suite"); break;
            case nameof(BenchProgress): Line("Bench", BenchProgress); break;
            case nameof(GrowLog): Line("Grow", GrowLog); if (_live.TryGetValue("Grow", out var g)) { g.Progress = GrowProgress; AnnounceProgress(g); } break;
            case nameof(RelaxLog): Line("Relax", RelaxLog); break;
            case nameof(MdLog): Line("Dynamics", MdLog); break;
            case nameof(PackLog): Line("Pack", PackLog); break;
            case nameof(RxLog): Line("React", RxLog); break;
            case nameof(EqLog):
                Line("Equilibrate", EqLog);
                if (_live.TryGetValue("Equilibrate", out var q))
                {
                    var m = System.Text.RegularExpressions.Regex.Match(EqLog, @"stage (\d+) of (\d+)");
                    if (m.Success) { q.Stage = int.Parse(m.Groups[1].Value); q.Stages = int.Parse(m.Groups[2].Value); q.Progress = (q.Stage - 1.0) / Math.Max(1, q.Stages); AnnounceProgress(q); }
                }
                break;
        }
    }

    private void Line(string kind, string text) { if (_live.TryGetValue(kind, out var j)) j.Add(text); }

    private void Track(string kind, bool running, int module, string title)
    {
        if (running && !_live.ContainsKey(kind))
        {
            Job job;
            if (_starting is { } q && q.Kind == kind)   // a queued job starting: the same entry runs
            {
                job = q;
                _starting = null;
                job.Status = "running";
                job.Started = DateTime.Now;
                Jobs.Remove(job);
            }
            else
            {
                var k = _jobCounters[kind] = _jobCounters.GetValueOrDefault(kind) + 1;
                var s = _doc?.Summary();
                job = new Job
                {
                    Id = $"{kind.ToLowerInvariant()}-{k}", Kind = kind, Module = module, Title = title,
                    Document = kind == "Grow" || kind == "Pack" ? "new cell" : Title, Atoms = kind == "Grow" || kind == "Pack" ? 0 : s?.Atoms ?? 0,
                    Provenance = Manifest(kind),
                };
            }
            if (kind is "Dynamics" or "Equilibrate") { job.CurveA = "Density"; job.AxisA = "density (g/cm³)"; job.CurveB = "Temperature"; job.AxisB = "temperature (K)"; job.AxisX = "time (ps)"; }
            else if (kind == "Relax") { job.CurveA = "Energy"; job.AxisA = "E (kcal/mol)"; job.CurveB = "Largest force"; job.AxisB = "log₁₀ |F|max"; job.AxisX = "iteration"; }
            else if (kind == "React") { job.CurveA = "Energy after each cycle"; job.AxisA = "E (kcal/mol)"; job.CurveB = "Force before the next cycle"; job.AxisB = "log₁₀ |F|max"; job.AxisX = "reaction cycle"; }
            _live[kind] = job;
            Jobs.Insert(0, job);
            if (kind is not ("Grow" or "Pack")) AttachJob(job, _activeItem);   // the structure it runs on
            job.Add($"{title} started on {job.Document}");
            SelectedJob = job;
        }
        else if (!running && _live.Remove(kind, out var job))
        {
            job.Ended = DateTime.Now;
            var last = job.Log.Count > 0 ? job.Log[^1].Text : "";
            var all = string.Join("\n", job.Log.TakeLast(6).Select(l => l.Text));
            if (all.Contains("ancelled")) job.Status = "cancelled";
            else if (all.Contains("Could not") || all.Contains("failed") || (kind == "Pack" && all.Contains("could not pack")))
            {
                job.Status = "failed";
                job.Error = string.Join("\n", job.Log.SkipWhile(l => !l.Text.StartsWith("Could not")).Select(l => l.Text));
                if (job.Error.Length == 0) job.Error = job.Log.LastOrDefault(l => l.Text.Contains("failed"))?.Text ?? last;
                (job.Suggestion, job.SuggestModule) = Suggest(kind, job.Error);
                (job.Fix, job.FixText) = FixFor(kind, job.Error);
                if (kind == "Dynamics" && MdCanContinue)   // FailedJob: nothing is lost — the checkpoint is intact
                    (job.Suggestion, job.SuggestModule) = ($"{MdCheckpointText}. {MdContinueLabel} on the Dynamics page" +
                        (job.Suggestion.Length > 0 ? " — " + char.ToLowerInvariant(job.Suggestion[0]) + job.Suggestion[1..] : "."), 3);
            }
            else if (all.Contains("not reached") || all.Contains("did not pass")) job.Status = "stopped";
            else job.Status = "done";
            if (job.Status == "done") job.Progress = 1;
            // the reader hears it once, politely (AccessibilityMap: "Growth finished, 20 chains, no close contacts")
            Announcement = $"{job.Title} {(job.Status == "done" ? "finished" : job.Status)}" + (last.Length > 0 ? ". " + last : ".");
            if (_doc != null && job.Atoms == 0) job.Provenance.Add(new JobFact("result", $"{_doc.Summary().Atoms:N0} atoms"));
            AttachJob(job, _activeItem);   // Grow and Pack: the cell they made
            BuildJobOutputs(job);
            SaveJobs();
            if (_settings.NotifyRuns && job.Ended is { } end && end - job.Started >= TimeSpan.FromMinutes(1))
                SystemNotify("CAPS · " + job.Title + " " + (job.Status == "done" ? "finished" : job.Status), last.Length > 0 ? last : job.Duration);
            if (_live.Count == 0) { ResumeRun(); Avalonia.Threading.Dispatcher.UIThread.Post(() => _ = StartNextQueued()); }
        }
        Raise(nameof(HasJobs));
        Raise(nameof(JobsSummary));
    }

    /// <summary>Recovery for the failures CAPS explains in its messages.</summary>
    internal (string Fix, string Text) FixFor(string kind, string error)
    {
        var e = error.ToLowerInvariant();
        var inv = CultureInfo.InvariantCulture;
        if (kind == "Dynamics" && (e.Contains("nan") || e.Contains("blew up") || e.Contains("moved") || e.Contains("unstable")) && _mdDt > 0.25)
            return ("md-dt", string.Format(inv, "Δt {0:0.##} → {1:0.##} fs, {2}", _mdDt, Math.Max(0.25, _mdDt / 2), MdCanContinue ? "from the last checkpoint" : "from the start"));
        if (kind == "React" && (e.Contains("failed at cycle") || e.Contains("force") || e.Contains("moved")))
        {
            var md = Math.Max(_rxMdPs * 4, 5);
            return ("rx-gentle", string.Format(inv, "{0} → {1} reactions a cycle, MD after each {2:0.#} → {3:0.#} ps{4}, from the structure kept after the last completed cycle",
                _rxPerCycle, Math.Max(1, _rxPerCycle / 2), _rxMdPs, md, _rxCapture > 4 ? string.Format(inv, ", capture {0:0.#} → 4 Å", _rxCapture) : ""));
        }
        if (kind == "Pack" && (e.Contains("could not pack") || e.Contains("tolerance")) && IsCapsPack(_packText))
            return ("pack-box", string.Format(inv, "box {0:0.#} × {1:0.#} × {2:0.#} → {3:0.#} × {4:0.#} × {5:0.#} Å (regions keep their numbers)", _packX, _packY, _packZ, _packX * 1.1, _packY * 1.1, _packZ * 1.1));
        return ("", "");
    }

    /// <summary>Restart with fix: makes the job's change and runs it again (the new run is a new job).</summary>
    public async Task RestartWithFix(Job job)
    {
        if (!job.HasFix || !Idle) { Status = job.HasFix ? "Wait for the running job to finish" : "This failure has no automatic fix"; return; }
        var inv = CultureInfo.InvariantCulture;
        switch (job.Fix)
        {
            case "md-dt":
                MdDtD = (decimal)Math.Max(0.25, _mdDt / 2);
                Status = $"Restarting Dynamics with Δt {_mdDt.ToString("0.##", inv)} fs";
                if (MdCanContinue) await ContinueMd(); else await RunMd();
                break;
            case "rx-gentle":
                RxPerCycleD = Math.Max(1, _rxPerCycle / 2);
                RxMdPsD = (decimal)Math.Max(_rxMdPs * 4, 5);
                RxRelax = true;
                if (_rxCapture > 4) RxCaptureD = 4;
                Status = $"Restarting React: {_rxPerCycle} reactions a cycle, {_rxMdPs.ToString("0.#", inv)} ps after each";
                await RunReact();
                break;
            case "pack-box":
                var (x0, y0, z0) = (_packX, _packY, _packZ);
                PackXD = (decimal)(x0 * 1.1); PackYD = (decimal)(y0 * 1.1); PackZD = (decimal)(z0 * 1.1);
                string F(double v) => v.ToString("0.###", inv);
                var text = _packText;
                // the cell line and every whole-cell box take the new edges
                text = System.Text.RegularExpressions.Regex.Replace(text, @"^(\s*cell\s+)\S+\s+\S+\s+\S+", m => m.Groups[1].Value + $"{F(_packX)} {F(_packY)} {F(_packZ)}", System.Text.RegularExpressions.RegexOptions.Multiline);
                text = text.Replace($"to {F(x0)} {F(y0)} {F(z0)}", $"to {F(_packX)} {F(_packY)} {F(_packZ)}");
                PackText = text;
                Status = $"Restarting Pack in a {F(_packX)} × {F(_packY)} × {F(_packZ)} Å box";
                await RunPack();
                break;
        }
    }

    private static (string, int) Suggest(string kind, string error)
    {
        var e = error.ToLowerInvariant();
        if (e.Contains("no parameters") || e.Contains("types c and h only") || e.Contains("untyped") || e.Contains("missing"))
            return ("Assign a force field in Field (GAFF2 types O, N, S and more); Relax, Dynamics and Equilibrate then use it.", 7);
        if (e.Contains("nan") || e.Contains("blew up") || e.Contains("moved") || e.Contains("unstable"))
            return (kind == "Dynamics" || kind == "Equilibrate"
                ? "Relax the structure first (forces from overlaps make the first steps unstable), or halve the timestep." : "", 2);
        if (kind == "React" && e.Contains("failed at cycle"))
            return ("Nothing is lost: the structure after the last completed cycle is kept — React again to continue from it (fewer reactions per cycle, or more relaxation, if the same cycle fails).", 6);
        if (e.Contains("could not pack") || e.Contains("tolerance"))
            return ("Use a larger box or fewer molecules; Pack never returns a cell with contacts closer than the tolerance.", 5);
        if (e.Contains("periodic cell"))
            return ("This step needs a periodic cell: open a structure with a box, or build one with Grow or Pack.", -1);
        return ("", -1);
    }

    private void OnThermo()
    {
        foreach (var kind in new[] { "Dynamics", "Equilibrate" })
        {
            if (!_live.TryGetValue(kind, out var j)) continue;
            j.A.Clear();
            j.B.Clear();
            foreach (var r in _thermo) { if (r.Density > 0) j.A.Add((r.TimePs, r.Density)); j.B.Add((r.TimePs, r.Temperature)); }
            if (kind == "Dynamics" && _thermo.Count > 0 && _mdSteps > 0) { j.Progress = (double)_thermo[^1].Step / _mdSteps; AnnounceProgress(j); }
            j.CurvesChanged();
            if (j == _job) JobCurvesChanged?.Invoke();
        }
    }

    /// <summary>React jobs: the energy and the largest force after each cycle's relaxation (FailedJob: the force before the failure).</summary>
    private void OnReactCurves()
    {
        if (!_live.TryGetValue("React", out var j)) return;
        j.A.Clear();
        j.B.Clear();
        foreach (var r in _rxRows)
        {
            if (r.Energy != 0) j.A.Add((r.Cycle, r.Energy));
            if (r.MaxForce > 0) j.B.Add((r.Cycle, Math.Log10(r.MaxForce)));
        }
        if (_rxCycles > 0) { j.Progress = Math.Min(1, (double)(_rxRows.Count > 0 ? _rxRows[^1].Cycle : 0) / _rxCycles); AnnounceProgress(j); }
        j.CurvesChanged();
        if (j == _job) JobCurvesChanged?.Invoke();
    }

    private void OnRelaxCurves()
    {
        if (!_live.TryGetValue("Relax", out var j)) return;
        j.A.Clear(); j.A.AddRange(_relaxEnergy);
        j.B.Clear(); j.B.AddRange(_relaxForce);
        j.CurvesChanged();
        if (j == _job) JobCurvesChanged?.Invoke();
    }

    /// <summary>What a run needs to be repeated: the program, the input, the force field, the settings.</summary>
    private List<JobFact> Manifest(string kind)
    {
        var inv = CultureInfo.InvariantCulture;
        var f = new List<JobFact>
        {
            new("caps", $"{typeof(MainViewModel).Assembly.GetName().Version?.ToString(3)} · ABI {Native.AbiVersion()}"),
        };
        if (_doc != null && kind is not ("Grow" or "Pack"))
        {
            f.Add(new("input", File.Exists(_doc.Path) ? RecentFiles.Tilde(Path.GetFullPath(_doc.Path)) : _doc.Path));
            if (File.Exists(_doc.Path))
                try
                {
                    using var fs = File.OpenRead(_doc.Path);
                    var h = Convert.ToHexString(SHA256.HashData(fs)).ToLowerInvariant();
                    f.Add(new("sha256", h[..4] + "…" + h[^4..]));
                }
                catch { }
            f.Add(new("frame", _frame.ToString(inv)));
        }
        f.Add(new("ff", Field.Assigned ? Field.ForceFieldName : "built-in GAFF (C and H) or UFF (other elements)"));
        switch (kind)
        {
            case "Relax":
                f.Add(new("minimiser", Minimisers[Math.Clamp(RelaxMethod, 0, Minimisers.Length - 1)]));
                f.Add(new("electrostatics", RelaxCoulomb ? $"{ElectrostaticsText} · rc {_relaxCutoff:F0} Å" : "off"));
                f.Add(new("vdW", $"cut {_relaxCutoff:F0} Å"));
                break;
            case "Dynamics":
                f.Add(new("electrostatics", ElectrostaticsText));
                f.Add(new("ensemble", Ensembles[Math.Clamp(MdEnsemble, 0, Ensembles.Length - 1)]));
                f.Add(new("timestep", $"{_mdDt:F2} fs · {MdStepsD:N0} steps"));
                f.Add(new("temperature", $"{_mdTemp:F0} K"));
                f.Add(new("seed", _mdSeed.ToString(inv)));
                break;
            case "Equilibrate":
                f.Add(new("protocol", Protocols[Math.Clamp(EqProtocol, 0, Protocols.Length - 1)]));
                f.Add(new("protocol sha256", Convert.ToHexString(SHA256.HashData(System.Text.Encoding.UTF8.GetBytes(EqText))).ToLowerInvariant()[..8]));
                break;
            case "Grow":
                f.Add(new("chains", $"{GrowChains} × DP {GrowDp} · {Tacticities[Math.Clamp(GrowTacticity, 0, Tacticities.Length - 1)].ToLowerInvariant()}"));
                f.Add(new("seed", GrowSeed.ToString(inv)));
                break;
        }
        f.Add(new("precision", "double"));
        f.Add(new("threads", _settings.Threads == 0 ? $"auto ({Math.Min(16, Environment.ProcessorCount)})" : _settings.Threads.ToString(inv)));
        return f;
    }

    private void SaveJobs()
    {
        try
        {
            Directory.CreateDirectory(Path.GetDirectoryName(JobsFile)!);
            var arr = new JsonArray(Jobs.Where(j => (!j.IsRunning || j.IsRemote) && !j.IsQueued).Take(40).Select(j => (JsonNode)j.ToJson(false)).ToArray());
            File.WriteAllText(JobsFile, arr.ToJsonString(JsonOut));
        }
        catch { }
    }

    // ---------------------------------------------------------------- pause (design/boards/Jobs "Pause")
    // The run's progress callback waits at this gate while paused: the core holds its state, nothing is lost, and
    // Resume carries on from the same step. Cancel opens the gate so the run can stop.
    private readonly ManualResetEventSlim _runGate = new(true);
    private bool _runPaused;
    public bool RunPaused { get => _runPaused; private set { if (Set(ref _runPaused, value)) { Raise(nameof(PauseLabel)); foreach (var j in _live.Values) j.Paused = value; } } }
    public string PauseLabel => _runPaused ? "Resume" : "Pause";
    public bool CanPause => MdRunning || EqRunning || Relaxing;

    /// <summary>On the run's thread: holds while paused (a cancel still gets through).</summary>
    private void WaitIfPaused(CancellationToken token)
    {
        while (!_runGate.Wait(200)) if (token.IsCancellationRequested) return;
    }

    public void PauseRun()
    {
        if (!CanPause || _runPaused) return;
        _runGate.Reset();
        RunPaused = true;
        Status = "Paused at the run's next report · Resume carries on from the same step";
    }

    public void ResumeRun()
    {
        _runGate.Set();
        if (_runPaused) { RunPaused = false; Status = CanPause ? "Resumed" : Status; }
    }

    public void TogglePause() { if (_runPaused) ResumeRun(); else PauseRun(); }

    // ---------------------------------------------------------------- queue (design/boards/Jobs "queued")
    // Dynamics and Equilibrate runs queued while another run goes, with their settings as they were when queued; each
    // starts when the one before it ends, on the structure it was queued on.
    private sealed record QueuedRun(Job Job, CapsDocument? Doc, Func<Task> Start);   // Doc null: the run makes a new structure
    private readonly List<QueuedRun> _queue = new();
    private Job? _starting;
    public int QueuedCount => _queue.Count;
    public string QueueText => _queue.Count == 0 ? "" : $"{_queue.Count} queued";

    private void EnqueueNew(string kind, int module, string title, Func<Task> start) => Enqueue(kind, module, title, start, newStructure: true);

    private void Enqueue(string kind, int module, string title, Func<Task> start, bool newStructure = false, string? note = null)
    {
        if (_doc == null && !newStructure) return;
        var k = _jobCounters[kind] = _jobCounters.GetValueOrDefault(kind) + 1;
        var job = new Job
        {
            Id = $"{kind.ToLowerInvariant()}-{k}", Kind = kind, Module = module, Title = title, Document = newStructure ? "new cell" : Title,
            Atoms = newStructure ? 0 : _doc!.Summary().Atoms, Provenance = Manifest(kind),
        };
        job.Status = "queued";
        job.Add(note ?? (newStructure ? $"{title} queued; it builds a new structure when the current run ends, with the settings as they are now"
                                      : $"{title} queued on {job.Document}; it starts when the current run ends, with the settings as they are now"));
        _queue.Add(new QueuedRun(job, newStructure ? null : _doc, start));
        Jobs.Insert(0, job);
        Raise(nameof(HasJobs)); Raise(nameof(JobsSummary)); Raise(nameof(QueuedCount)); Raise(nameof(QueueText));
        Status = $"{title} queued · {_queue.Count} in the queue";
        if (Idle) _ = StartNextQueued();
    }

    public void QueueMd()
    {
        var opts = MdOptions(null);
        var ens = MdEnsemble;
        Enqueue("Dynamics", 3, "Dynamics · " + Ensembles[Math.Clamp(ens, 0, Ensembles.Length - 1)], () => RunMd(null, opts, ens));
    }

    public void QueueEquilibrate()
    {
        var text = EqText;
        var opts = EqOptions(null);
        var target = EqTargetCurve();
        var until = opts.UntilConverged != 0;
        Enqueue("Equilibrate", 4, "Equilibrate · " + Protocols[Math.Clamp(EqProtocol, 0, Protocols.Length - 1)], () => RunEquilibrate(text, until, opts, target));
    }

    public void QueueRelax()
    {
        var opts = RelaxOptions();
        Enqueue("Relax", 2, "Relax · " + Minimisers[Math.Clamp(RelaxMethod, 0, Minimisers.Length - 1)], () => Relax(opts));
    }

    /// <summary>Pack queued with its input as it is now; the cell it packs becomes a new structure.</summary>
    public void QueuePack()
    {
        string text;
        try { text = PackTextToRun(); } catch (Exception e) { Status = "Could not queue the packing: " + e.Message; return; }
        var baseDir = PackBaseDir;
        EnqueueNew("Pack", 5, "Pack · molecules into a box", () => RunPack(text, baseDir));
    }

    /// <summary>The glass transition queued on this structure: its replicas run when the current run ends, with the Glass
    /// page's settings as they are then (its cooling scan reads them when it starts).</summary>
    public void QueueGlass()
    {
        var reps = GtReplicas;
        Enqueue("Glass", 47, $"Glass transition · {reps} replica{(reps == 1 ? "" : "s")}", RunGlass,
                note: $"Glass transition queued on {Title}; it starts when the current run ends and uses the Glass page's settings at that moment");
    }

    public void QueueReact()
    {
        var opts = ReactOptions();
        var text = _rxText;
        Enqueue("React", 6, "React · crosslinking", () => RunReact(opts, text));
    }

    /// <summary>Grow queued as its recipe (the settings as they are now): it builds a new cell when its turn comes, with
    /// the files the recipe exports in the queue folder.</summary>
    public void QueueGrow()
    {
        var text = GrowRecipe();
        var dir = Path.Combine(QueueFolder, $"grow-{DateTime.Now:yyyyMMdd-HHmmss}");
        var named = text.Split('\n').FirstOrDefault(l => l.StartsWith("name:"))?[5..].Trim().Trim('"') ?? "";
        var label = named.Length > 0 ? named + " cell" : "grown cell";
        EnqueueNew("Recipe", 0, "Grow · " + label + " (recipe)", () => { Directory.CreateDirectory(dir); return RunRecipeText(text, label, dir); });
    }
    public static string QueueFolder => AppSettings.Override != null ? Path.Combine(Path.GetDirectoryName(AppSettings.Override)!, "caps-queue") : Path.Combine(AppSettings.Folder, "queue");

    private async Task StartNextQueued()
    {
        while (_queue.Count > 0 && Idle)
        {
            var q = _queue[0];
            _queue.RemoveAt(0);
            Raise(nameof(QueuedCount)); Raise(nameof(QueueText));
            if (q.Doc != null && !ReferenceEquals(q.Doc, _doc))
            {
                q.Job.Status = "cancelled";
                q.Job.Add("Not started: the structure it was queued on is no longer the open one");
                q.Job.Ended = DateTime.Now;
                continue;
            }
            _starting = q.Job;
            await q.Start();
            if (q.Job.IsQueued)   // the page refused (a field missing …): it never started
            {
                _starting = null;
                q.Job.Status = "failed";
                q.Job.Error = "Could not start: " + Status;
                q.Job.Add(q.Job.Error);
                q.Job.Ended = DateTime.Now;
                continue;
            }
            return;   // the next one starts when this run ends
        }
        Raise(nameof(JobsSummary));
    }

    public void CancelJob(Job? j)
    {
        if (j is { IsQueued: true })
        {
            _queue.RemoveAll(q => q.Job == j);
            j.Status = "cancelled";
            j.Add("Removed from the queue");
            j.Ended = DateTime.Now;
            Raise(nameof(QueuedCount)); Raise(nameof(QueueText)); Raise(nameof(JobsSummary));
            return;
        }
        switch (j?.Kind)
        {
            case "Relax": CancelRelax(); break;
            case "Dynamics": CancelMd(); break;
            case "Equilibrate": CancelEquilibrate(); break;
            case "Pack": CancelPack(); break;
            case "React": CancelReact(); break;
            case "Grow": CancelGrow(); break;
            case "Analyze": Analyze.Cancel(); break;
            case "Bench": CancelBench(); break;
        }
    }

    public void ClearJobs()
    {
        foreach (var j in Jobs.Where(j => !j.IsRunning && !j.IsQueued).ToList()) Jobs.Remove(j);
        SelectedJob = Jobs.FirstOrDefault();
        SaveJobs();
        Raise(nameof(HasJobs));
        Raise(nameof(JobsSummary));
    }

    public string JobJson(Job j) => j.ToJson(true).ToJsonString(JsonOut);
    private static readonly JsonSerializerOptions JsonOut = new() { WriteIndented = true, Encoder = System.Text.Encodings.Web.JavaScriptEncoder.UnsafeRelaxedJsonEscaping };

    /// <summary>A system notification (D13): macOS Notification Centre through osascript, notify-send on Linux; nothing on
    /// Windows (the status bar says it). Never blocks and never fails the run.</summary>
    internal static void SystemNotify(string title, string text)
    {
        try
        {
            static string Esc(string s) => s.Replace("\\", "\\\\").Replace("\"", "\\\"");
            System.Diagnostics.ProcessStartInfo? psi = null;
            if (OperatingSystem.IsMacOS())
            {
                psi = new System.Diagnostics.ProcessStartInfo("osascript") { UseShellExecute = false, CreateNoWindow = true };
                psi.ArgumentList.Add("-e");
                psi.ArgumentList.Add($"display notification \"{Esc(text)}\" with title \"{Esc(title)}\"");
            }
            else if (OperatingSystem.IsLinux())
            {
                psi = new System.Diagnostics.ProcessStartInfo("notify-send") { UseShellExecute = false, CreateNoWindow = true };
                psi.ArgumentList.Add(title);
                psi.ArgumentList.Add(text);
            }
            if (psi != null) System.Diagnostics.Process.Start(psi)?.Dispose();
        }
        catch { /* no notifier on this system */ }
    }
}
