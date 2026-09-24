using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Security.Cryptography;
using System.Text.Json;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

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
    public int Module { get; init; }
    public string Title { get; init; } = "";
    public string Document { get; init; } = "";
    public long Atoms { get; init; }
    public DateTime Started { get; init; } = DateTime.Now;
    private DateTime? _ended;
    public DateTime? Ended { get => _ended; set { _ended = value; Raise(nameof(Ended)); Raise(nameof(Duration)); } }

    private string _status = "running";
    /// <summary>running | done | failed | cancelled | stopped (finished without meeting its target)</summary>
    public string Status
    {
        get => _status;
        set
        {
            _status = value;
            foreach (var n in new[] { nameof(Status), nameof(IsRunning), nameof(IsFailed), nameof(IsDone), nameof(IsQuiet), nameof(StatusText) }) Raise(n);
        }
    }
    public bool IsRunning => _status == "running";
    public bool IsFailed => _status == "failed";
    public bool IsDone => _status == "done";
    public bool IsQuiet => _status is "cancelled" or "stopped";
    public string StatusText => _status == "stopped" ? "stopped" : _status;

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
    public string Subtitle => $"{Document} · {Atoms.ToString("N0", CultureInfo.InvariantCulture)} atoms · started {Started:HH:mm}";
    public string Where => $"{Id} · Local · {Environment.ProcessorCount} threads";
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
        j.Status = st == "running" ? "stopped" : st;   // a job running when the Studio closed did not finish
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
            return Jobs.Count == 0 ? "no jobs" : $"{r} running · {Jobs.Count(j => j.IsDone)} done{(f > 0 ? $" · {f} failed" : "")}";
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
        PropertyChanged += OnRunProperty;
        Analyze.PropertyChanged += (_, e) =>
        {
            if (e.PropertyName == nameof(AnalyzeViewModel.Working)) Track("Analyze", Analyze.Working, 1, "Analyze · properties");
            if (e.PropertyName == nameof(AnalyzeViewModel.Log) && _live.TryGetValue("Analyze", out var j)) j.Add(Analyze.Log);
        };
        ThermoChanged += OnThermo;
        RelaxCurvesChanged += OnRelaxCurves;
        Raise(nameof(HasJobs));
        Raise(nameof(JobsSummary));
    }

    private void OnRunProperty(object? s, PropertyChangedEventArgs e)
    {
        switch (e.PropertyName)
        {
            case nameof(Growing): Track("Grow", Growing, 0, "Grow · amorphous cell"); break;
            case nameof(Relaxing): Track("Relax", Relaxing, 2, "Relax · " + Minimisers[Math.Clamp(RelaxMethod, 0, Minimisers.Length - 1)]); break;
            case nameof(MdRunning): Track("Dynamics", MdRunning, 3, "Dynamics · " + Ensembles[Math.Clamp(MdEnsemble, 0, Ensembles.Length - 1)]); break;
            case nameof(EqRunning): Track("Equilibrate", EqRunning, 4, "Equilibrate · " + Protocols[Math.Clamp(EqProtocol, 0, Protocols.Length - 1)]); break;
            case nameof(Packing): Track("Pack", Packing, 5, "Pack · molecules into a box"); break;
            case nameof(Reacting): Track("React", Reacting, 6, "React · crosslinking"); break;
            case nameof(BenchRunning): Track("Bench", BenchRunning, 12, "Bench · validation suite"); break;
            case nameof(BenchProgress): Line("Bench", BenchProgress); break;
            case nameof(GrowLog): Line("Grow", GrowLog); if (_live.TryGetValue("Grow", out var g)) g.Progress = GrowProgress; break;
            case nameof(RelaxLog): Line("Relax", RelaxLog); break;
            case nameof(MdLog): Line("Dynamics", MdLog); break;
            case nameof(PackLog): Line("Pack", PackLog); break;
            case nameof(RxLog): Line("React", RxLog); break;
            case nameof(EqLog):
                Line("Equilibrate", EqLog);
                if (_live.TryGetValue("Equilibrate", out var q))
                {
                    var m = System.Text.RegularExpressions.Regex.Match(EqLog, @"stage (\d+) of (\d+)");
                    if (m.Success) { q.Stage = int.Parse(m.Groups[1].Value); q.Stages = int.Parse(m.Groups[2].Value); q.Progress = (q.Stage - 1.0) / Math.Max(1, q.Stages); }
                }
                break;
        }
    }

    private void Line(string kind, string text) { if (_live.TryGetValue(kind, out var j)) j.Add(text); }

    private void Track(string kind, bool running, int module, string title)
    {
        if (running && !_live.ContainsKey(kind))
        {
            var k = _jobCounters[kind] = _jobCounters.GetValueOrDefault(kind) + 1;
            var s = _doc?.Summary();
            var job = new Job
            {
                Id = $"{kind.ToLowerInvariant()}-{k}", Kind = kind, Module = module, Title = title,
                Document = kind == "Grow" || kind == "Pack" ? "new cell" : Title, Atoms = kind == "Grow" || kind == "Pack" ? 0 : s?.Atoms ?? 0,
                Provenance = Manifest(kind),
            };
            if (kind is "Dynamics" or "Equilibrate") { job.CurveA = "Density"; job.AxisA = "density (g/cm³)"; job.CurveB = "Temperature"; job.AxisB = "temperature (K)"; job.AxisX = "time (ps)"; }
            else if (kind == "Relax") { job.CurveA = "Energy"; job.AxisA = "E (kcal/mol)"; job.CurveB = "Largest force"; job.AxisB = "log₁₀ |F|max"; job.AxisX = "iteration"; }
            _live[kind] = job;
            Jobs.Insert(0, job);
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
                if (job.Error.Length == 0) job.Error = last;
                (job.Suggestion, job.SuggestModule) = Suggest(kind, job.Error);
            }
            else if (all.Contains("not reached") || all.Contains("did not pass")) job.Status = "stopped";
            else job.Status = "done";
            if (job.Status == "done") job.Progress = 1;
            if (_doc != null && job.Atoms == 0) job.Provenance.Add(new JobFact("result", $"{_doc.Summary().Atoms:N0} atoms"));
            SaveJobs();
        }
        Raise(nameof(HasJobs));
        Raise(nameof(JobsSummary));
    }

    /// <summary>Recovery for the failures CAPS explains in its messages.</summary>
    private static (string, int) Suggest(string kind, string error)
    {
        var e = error.ToLowerInvariant();
        if (e.Contains("no parameters") || e.Contains("types c and h only") || e.Contains("untyped") || e.Contains("missing"))
            return ("Assign a force field in Field (GAFF2 types O, N, S and more); Relax, Dynamics and Equilibrate then use it.", 7);
        if (e.Contains("nan") || e.Contains("blew up") || e.Contains("moved") || e.Contains("unstable"))
            return (kind == "Dynamics" || kind == "Equilibrate"
                ? "Relax the structure first (forces from overlaps make the first steps unstable), or halve the timestep." : "", 2);
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
            if (kind == "Dynamics" && _thermo.Count > 0 && _mdSteps > 0) j.Progress = (double)_thermo[^1].Step / _mdSteps;
            j.CurvesChanged();
            if (j == _job) JobCurvesChanged?.Invoke();
        }
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
        f.Add(new("ff", Field.Assigned ? Field.ForceFieldName : "built-in GAFF (C and H)"));
        switch (kind)
        {
            case "Relax":
                f.Add(new("minimiser", Minimisers[Math.Clamp(RelaxMethod, 0, Minimisers.Length - 1)]));
                f.Add(new("electrostatics", RelaxCoulomb ? $"DSF · rc {_relaxCutoff:F0} Å" : "off"));
                f.Add(new("vdW", $"cut {_relaxCutoff:F0} Å"));
                break;
            case "Dynamics":
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
            var arr = new JsonArray(Jobs.Where(j => !j.IsRunning).Take(40).Select(j => (JsonNode)j.ToJson(false)).ToArray());
            File.WriteAllText(JobsFile, arr.ToJsonString(JsonOut));
        }
        catch { }
    }

    public void CancelJob(Job? j)
    {
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
        foreach (var j in Jobs.Where(j => !j.IsRunning).ToList()) Jobs.Remove(j);
        SelectedJob = Jobs.FirstOrDefault();
        SaveJobs();
        Raise(nameof(HasJobs));
        Raise(nameof(JobsSummary));
    }

    public string JobJson(Job j) => j.ToJson(true).ToJsonString(JsonOut);
    private static readonly JsonSerializerOptions JsonOut = new() { WriteIndented = true, Encoder = System.Text.Encodings.Web.JavaScriptEncoder.UnsafeRelaxedJsonEscaping };
}
