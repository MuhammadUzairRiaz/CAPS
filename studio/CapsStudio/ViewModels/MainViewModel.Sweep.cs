using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One run of a sweep: its condition, seed, state and results.</summary>
public sealed class SweepRun : ObservableObject
{
    public string Tacticity { get; init; } = "";
    public int TacticityCode { get; init; }
    public int Dp { get; init; }
    public ulong Seed { get; init; }
    public string Path { get; init; } = "";
    public string Condition => $"{Tacticity} · DP {Dp}";
    private string _status = "queued", _detail = "";
    public string Status { get => _status; set { if (Set(ref _status, value)) { Raise(nameof(Fill)); Raise(nameof(Dim)); Raise(nameof(Tip)); } } }
    public string Detail { get => _detail; set { if (Set(ref _detail, value)) Raise(nameof(Tip)); } }
    public double Density { get; set; } = double.NaN;
    /// <summary>Set by the worker when the run ends; applied on the UI side after the sweep (no dispatcher round-trip).</summary>
    public string? FinalStatus { get; set; }
    public string FinalDetail { get; set; } = "";
    public double Rg { get; set; } = double.NaN;
    public double Tg { get; set; } = double.NaN;
    public string SeedText => Seed.ToString(CultureInfo.InvariantCulture);
    public IBrush Fill => Tokens.Brush(_status switch { "done" => "OkB", "running" => "AccB", "failed" => "ErrB", _ => "DimB" });
    public double Dim => _status is "queued" or "stopped" ? 0.45 : 1.0;
    public string Tip => $"{Condition} · seed {Seed}: {_status}" + (_detail.Length > 0 ? " · " + _detail : "");
}

/// <summary>A grid cell: one condition (tacticity × DP) and its seeds.</summary>
public sealed record SweepCell(int Dp, ObservableCollection<SweepRun> Runs);
/// <summary>A grid row: one tacticity across the chain lengths.</summary>
public sealed record SweepRow(string Tacticity, ObservableCollection<SweepCell> Cells);
/// <summary>Results of one condition over its seeds.</summary>
public sealed record SweepResult(string Condition, string Density, string Rg, string Tg, string Seeds);

/// <summary>Parameter sweep (design/boards/ParameterSweep): a grid of runs — tacticity × chain length × seed — each grown,
/// relaxed, optionally run in NPT, analysed (density, Rg) and saved with its provenance, a few at a time on this machine;
/// results by condition as mean ± sd over the seeds.</summary>
public sealed partial class MainViewModel
{
    public bool IsSweep => _module == 43;
    public static readonly string[] SweepTacticityNames = ["atactic", "isotactic", "syndiotactic"];
    public ObservableCollection<SweepRow> SweepRows { get; } = new();
    public ObservableCollection<SweepResult> SweepResults { get; } = new();
    public ObservableCollection<ProvRow> SweepSizes { get; } = new();
    private readonly List<SweepRun> _sweepRuns = new();
    private bool _swIso = true, _swSyn = true, _swAta = true, _swNpt, _swRunning, _swPaused;
    private string _swDps = "10, 20, 40", _swSeeds = "1, 2, 3", _swFolder = "", _swSummary = "", _swError = "";
    private decimal _swChains = 10, _swDensity = 0.5m, _swNptPs = 20, _swTemp = 300, _swConcurrency = 2;
    private FilmPolymer? _swPolymer;
    private CancellationTokenSource? _swCancel;

    public bool SweepIso { get => _swIso; set { if (Set(ref _swIso, value)) PlanSweep(); } }
    public bool SweepSyn { get => _swSyn; set { if (Set(ref _swSyn, value)) PlanSweep(); } }
    public bool SweepAta { get => _swAta; set { if (Set(ref _swAta, value)) PlanSweep(); } }
    public string SweepDps { get => _swDps; set { if (Set(ref _swDps, value ?? "")) PlanSweep(); } }
    public string SweepSeeds { get => _swSeeds; set { if (Set(ref _swSeeds, value ?? "")) PlanSweep(); } }
    public decimal SweepChains { get => _swChains; set { if (Set(ref _swChains, Math.Clamp(Math.Round(value), 1, 500))) PlanSweep(); } }
    public decimal SweepDensity { get => _swDensity; set => Set(ref _swDensity, Math.Clamp(value, 0.05m, 1.5m)); }
    public bool SweepNpt { get => _swNpt; set { if (Set(ref _swNpt, value)) RaiseSweep(); } }
    // each run equilibrated by the Equilibrate page's protocol, and its Tg by stepwise NPT cooling (Analyze's Tg)
    private bool _swEquil, _swTg;
    public bool SweepEquilibrate { get => _swEquil; set { if (Set(ref _swEquil, value)) RaiseSweep(); } }
    public bool SweepTg { get => _swTg; set { if (Set(ref _swTg, value)) RaiseSweep(); } }
    public decimal SweepNptPs { get => _swNptPs; set => Set(ref _swNptPs, Math.Clamp(value, 1, 100000)); }
    public decimal SweepTemperature { get => _swTemp; set => Set(ref _swTemp, Math.Clamp(value, 1, 2000)); }
    public decimal SweepConcurrency { get => _swConcurrency; set => Set(ref _swConcurrency, Math.Clamp(Math.Round(value), 1, 8)); }
    public FilmPolymer? SweepPolymer { get => _swPolymer ?? SweepPolymers.FirstOrDefault(); set { if (Set(ref _swPolymer, value)) PlanSweep(); } }
    /// <summary>Every homopolymer of the library (rubbers first).</summary>
    public ObservableCollection<FilmPolymer> SweepPolymers { get; } = new();
    public string SweepFolder { get => _swFolder; set => Set(ref _swFolder, value ?? ""); }
    public string SweepSummary { get => _swSummary; private set => Set(ref _swSummary, value); }
    public string SweepError { get => _swError; private set => Set(ref _swError, value); }
    public bool SweepRunning { get => _swRunning; private set { if (Set(ref _swRunning, value)) { Raise(nameof(SweepIdle)); Raise(nameof(SweepPauseText)); } } }
    public bool SweepIdle => !_swRunning;
    public string SweepPauseText => _swPaused ? "Resume sweep" : "Pause sweep";
    public string SweepTitle => $"Sweep · tacticity × chain length · {SweepPolymer?.Name.Split(" (")[0] ?? "polymer"}";
    public string SweepSubtitle => (_swCombine == 0
        ? $"{_sweepRuns.Count} runs = {SweepTacticities().Count} tacticit{(SweepTacticities().Count == 1 ? "y" : "ies")} × {ParseInts(_swDps).Count} DP × {ParseInts(_swSeeds).Count} seeds"
        : $"{_sweepRuns.Count} runs · {SweepCombines[_swCombine].ToLowerInvariant()} × {ParseInts(_swSeeds).Count} seeds") + $" · Grow → Relax{(_swEquil ? " → Equilibrate" : "")}{(_swNpt ? " → NPT" : "")} → Analyze{(_swTg ? " → Tg" : "")}";
    public string SweepDpHeader(int k) => k < ParseInts(_swDps).Count ? "DP " + ParseInts(_swDps)[k] : "";
    public ObservableCollection<string> SweepDpHeaders { get; } = new();

    private List<(string Name, int Code)> SweepTacticities()
    {
        var l = new List<(string, int)>();
        if (_swIso) l.Add(("isotactic", 1));
        if (_swSyn) l.Add(("syndiotactic", 2));
        if (_swAta) l.Add(("atactic", 0));
        return l;
    }
    private static List<int> ParseInts(string s) => s.Split([',', ' ', ';'], StringSplitOptions.RemoveEmptyEntries)
        .Select(x => int.TryParse(x, NumberStyles.Integer, CultureInfo.InvariantCulture, out var v) ? v : -1).Where(v => v > 0).Distinct().Take(12).ToList();

    public void OpenSweep()
    {
        LoadSurface();   // the polymer library
        if (_swFolder.Length == 0)
            SweepFolder = System.IO.Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "CAPS", "sweeps",
                                                 "sweep-" + DateTime.Now.ToString("yyyyMMdd-HHmm", CultureInfo.InvariantCulture));
        // tacticity needs a stereocentre in the backbone: polystyrene first, else the library's first polymer
        if (SweepPolymers.Count == 0)
            foreach (var p in PolymerLibrary.Where(p => !p.Copolymer).OrderByDescending(p => p.Rubber))
                SweepPolymers.Add(new FilmPolymer(p.Name, new JsonObject
                {
                    ["units"] = new JsonArray(new JsonObject { ["name"] = p.Name, ["smiles"] = p.Smiles }),
                    ["sequence"] = "homopolymer",
                }.ToJsonString()));
        _swPolymer ??= SweepPolymers.FirstOrDefault(p => p.Name == "Polystyrene") ?? SweepPolymers.FirstOrDefault();
        Raise(nameof(SweepPolymer));
        SetModule(43);
        if (!_swRunning) PlanSweep();
    }

    // Combine as: every tacticity × chain length (full grid); one factor at a time from the first of each (the first
    // tacticity at every length, every tacticity at the first length); or paired (the i-th tacticity with the i-th length)
    public static readonly string[] SweepCombines = ["Full grid", "One factor at a time", "Paired (i-th with i-th)"];
    private int _swCombine;
    public int SweepCombine { get => _swCombine; set { if (Set(ref _swCombine, Math.Clamp(value, 0, 2))) PlanSweep(); } }

    /// <summary>The runs for the parameters as combined (nothing runs yet).</summary>
    private void PlanSweep()
    {
        if (_swRunning) return;
        _sweepRuns.Clear();
        SweepRows.Clear();
        SweepDpHeaders.Clear();
        var dps = ParseInts(_swDps);
        var seeds = ParseInts(_swSeeds);
        foreach (var dp in dps) SweepDpHeaders.Add("DP " + dp);
        var ti = -1;
        foreach (var (name, code) in SweepTacticities())
        {
            ++ti;
            var cells = new ObservableCollection<SweepCell>();
            var di = -1;
            foreach (var dp in dps)
            {
                ++di;
                var runs = new ObservableCollection<SweepRun>();
                var planned = _swCombine switch { 1 => ti == 0 || di == 0, 2 => ti == di, _ => true };
                foreach (var seed in planned ? seeds : [])
                {
                    var run = new SweepRun { Tacticity = name, TacticityCode = code, Dp = dp, Seed = (ulong)seed,
                                             Path = System.IO.Path.Combine(_swFolder, $"{name}_dp{dp}_seed{seed}.data") };
                    runs.Add(run);
                    _sweepRuns.Add(run);
                }
                cells.Add(new SweepCell(dp, runs));
            }
            SweepRows.Add(new SweepRow(name, cells));
        }
        SweepSizes.Clear();
        if (SweepPolymer is { } poly)
            foreach (var dp in dps)
            {
                var atoms = "—";
                try
                {
                    var spec = JsonNode.Parse(poly.Spec)!.AsObject();
                    spec["dp"] = dp;
                    var p = JsonNode.Parse(CapsDocument.ChainPreview(spec.ToJsonString(), 1))!;
                    if (p["atoms"]?.GetValue<double>() is double a) atoms = ((int)(a * (double)_swChains)).ToString("N0", CultureInfo.InvariantCulture);
                }
                catch { }
                SweepSizes.Add(new ProvRow(dp.ToString(CultureInfo.InvariantCulture), atoms, (seeds.Count * SweepTacticities().Count).ToString(CultureInfo.InvariantCulture)));
            }
        RaiseSweep();
    }

    private void RaiseSweep()
    {
        int n(string s) => _sweepRuns.Count(r => r.Status == s);
        SweepSummary = $"{n("done")} done · {n("failed")} failed · {n("running")} running · {n("queued")} queued";
        foreach (var p in new[] { nameof(SweepTitle), nameof(SweepSubtitle) }) Raise(p);
        SweepResults.Clear();
        var inv = CultureInfo.InvariantCulture;
        static string Ms(IList<double> v, string f) => v.Count == 0 ? "—" : v.Count == 1 ? v[0].ToString(f, CultureInfo.InvariantCulture)
            : $"{v.Average().ToString(f, CultureInfo.InvariantCulture)} ± {Math.Sqrt(v.Sum(x => (x - v.Average()) * (x - v.Average())) / (v.Count - 1)).ToString(f, CultureInfo.InvariantCulture)}";
        foreach (var g in _sweepRuns.GroupBy(r => r.Condition))
        {
            var done = g.Where(r => r.Status == "done").ToList();
            var failed = g.Count(r => r.Status == "failed");
            SweepResults.Add(new SweepResult(g.Key, Ms(done.Where(r => double.IsFinite(r.Density)).Select(r => r.Density).ToList(), "0.000"),
                                             Ms(done.Where(r => double.IsFinite(r.Rg)).Select(r => r.Rg).ToList(), "0.00"),
                                             Ms(done.Where(r => double.IsFinite(r.Tg)).Select(r => r.Tg).ToList(), "0"),
                                             $"{done.Count} / {g.Count()}" + (failed > 0 ? $" · {failed} failed" : "")));
        }
    }

    // where the runs go: this machine (a few at a time) or a host from Settings › Compute & remote, one recipe job per cell
    private int _swHost;
    public List<string> SweepHosts => RunWhereChoices;
    public int SweepHostIndex { get => Math.Min(_swHost, _settings.Hosts.Count); set { if (Set(ref _swHost, Math.Clamp(value, 0, _settings.Hosts.Count))) RaiseSweep(); } }

    /// <summary>One cell as a caps run recipe (JSON): grow, relax, the optional equilibration and NPT, analysis, export.</summary>
    private string SweepRecipe(FilmPolymer poly, SweepRun run)
    {
        var spec = JsonNode.Parse(poly.Spec)!.AsObject();
        var units = new JsonArray();
        foreach (var u in (spec["units"] as JsonArray) ?? new JsonArray()) if (u?["smiles"]?.GetValue<string>() is string smi) units.Add(smi);
        var build = new JsonObject { ["units"] = units, ["sequence"] = (string?)spec["sequence"] ?? "homopolymer", ["dp"] = run.Dp, ["chains"] = (int)_swChains,
                                     ["tacticity"] = run.Tacticity.ToLowerInvariant() };
        var r = new JsonObject
        {
            ["recipe"] = 1, ["name"] = System.IO.Path.GetFileNameWithoutExtension(run.Path), ["seed"] = run.Seed,
            ["build"] = new JsonObject { ["polymer"] = build }, ["type"] = new JsonObject { ["forcefield"] = "default" },
            ["grow"] = new JsonObject { ["density"] = (double)_swDensity, ["seed"] = run.Seed, ["contact_scale"] = "auto", ["curve"] = true },
            ["relax"] = new JsonObject { ["method"] = "lbfgs", ["fmax"] = 1.0 },
        };
        if (_swEquil) r["equilibrate"] = new JsonObject { ["protocol_text"] = _eqText, ["seed"] = run.Seed };
        if (_swNpt) r["md"] = new JsonObject { ["ensemble"] = "npt", ["ps"] = (double)_swNptPs, ["temperature"] = (double)_swTemp, ["seed"] = run.Seed };
        var props = new JsonArray("density", "rg");
        if (_swTg) props.Add("tg");
        r["analyze"] = new JsonObject { ["properties"] = props };
        r["export"] = new JsonArray("lammps");
        return r.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true });
    }

    /// <summary>Runs every queued run, a few at a time; each saves its cell with provenance.</summary>
    public async Task RunSweep()
    {
        if (_swRunning || SweepPolymer is not { } poly) return;
        if (_sweepRuns.Count == 0) { SweepError = "No runs: choose at least one tacticity, chain length and seed"; return; }
        if (SweepHostIndex > 0)
        {   // a host: each queued cell goes as its own recipe job; Jobs follows them and brings the results back
            if (_swEquil && _eqText.Trim().Length == 0) { SweepError = "Equilibrate needs a protocol: choose one on the Equilibrate page"; return; }
            var h = _settings.Hosts[SweepHostIndex - 1];
            SweepError = "";
            SweepRunning = true;
            // one array job: a task per cell (caps job new --array), each followed in Jobs and fetched when it ends
            var todo = _sweepRuns.Where(r => r.Status is "queued" or "stopped" or "failed").ToList();
            try
            {
                var jobs = await SendRecipeArray(h, "Sweep", todo.Select(run => ($"Sweep · {run.Condition} · seed {run.Seed}", System.IO.Path.GetFileNameWithoutExtension(run.Path), SweepRecipe(poly, run))).ToList());
                for (var i = 0; i < todo.Count && i < jobs.Count; ++i) { todo[i].Status = "running"; todo[i].Detail = $"on {h.Name} · {jobs[i].Remote?.JobId}"; }
            }
            catch (Exception e)
            {
                SweepError = "Could not send the sweep: " + e.Message;   // the host is unreachable: the runs stay queued
            }
            RaiseSweep();
            SweepRunning = false;
            RaiseSweep();
            Status = $"Sweep sent to {h.Name} · Jobs follows the runs and fetches each cell when it finishes";
            return;
        }
        Directory.CreateDirectory(_swFolder);
        SweepError = "";
        SweepRunning = true;
        _swPaused = false;
        Raise(nameof(SweepPauseText));
        _swCancel = new CancellationTokenSource();
        var token = _swCancel.Token;
        var conc = (int)_swConcurrency;
        var threads = Math.Max(1, Environment.ProcessorCount / conc);
        var gate = new SemaphoreSlim(conc);
        var tasks = new List<Task>();
        var (npt, ps, temp, chains, density) = (_swNpt, (double)_swNptPs, (double)_swTemp, (int)_swChains, (double)_swDensity);
        var (equil, tg, eqText) = (_swEquil, _swTg, _eqText);
        if (equil && eqText.Trim().Length == 0) { SweepError = "Equilibrate needs a protocol: choose one on the Equilibrate page"; SweepRunning = false; return; }
        var eqOpts = EqOptions(false);
        foreach (var run in _sweepRuns.Where(r => r.Status is "queued" or "stopped" or "failed").ToList())
        {
            await gate.WaitAsync();
            while (_swPaused && !token.IsCancellationRequested) await Task.Delay(200);
            if (token.IsCancellationRequested) { gate.Release(); break; }
            run.Status = "running";
            run.Detail = "growing";
            RaiseSweep();
            tasks.Add(Task.Run(() =>
            {
                try
                {
                    var spec = JsonNode.Parse(poly.Spec)!.AsObject();
                    spec["dp"] = run.Dp;
                    var g = new CapsGrowOpts { Chains = chains, Dp = run.Dp, Tacticity = run.TacticityCode, Seed = run.Seed, Density = density, ContactScale = -1.0, Curve = 1 };
                    var (doc, _) = CapsDocument.GrowChains(spec.ToJsonString(), g, (_, _, _) => !token.IsCancellationRequested, run.Condition);
                    using (doc)
                    {
                        Post(() => run.Detail = "relaxing");
                        doc.Relax(new CapsRelaxOpts { Method = 2, Ftol = 1.0, MaxIterations = 3000, Pushoff = 1, Cutoff = 10, Coulomb = 1, Threads = threads },
                                  (_, _, _, _, _, _) => !token.IsCancellationRequested);
                        if (equil)
                        {
                            Post(() => run.Detail = "equilibrating");
                            doc.Equilibrate(eqText, eqOpts with { Threads = threads, Seed = run.Seed }, (_, _, _, _) => !token.IsCancellationRequested);
                        }
                        if (npt)
                        {
                            Post(() => run.Detail = $"NPT {ps:0} ps");
                            doc.Md(new CapsMdOpts { Dt = 1, Steps = (long)(ps * 1000), Temperature = temp, Thermostat = 1, TauT = 100, Barostat = 1, Pressure = 1, TauP = 1000,
                                                    NewVelocities = 1, Seed = run.Seed, ThermoEvery = 1000, FrameEvery = 0, Cutoff = 10, Coulomb = 1, Tail = 1, Threads = threads },
                                   (_, _) => !token.IsCancellationRequested);
                        }
                        Post(() => run.Detail = "analysing");
                        var ao = new CapsAnalyzeOpts { First = 0, Last = -1, Stride = 1, Blocks = 5, Threads = threads };
                        var json = doc.Analyze("density,rg", ao, null);
                        double Value(string id)
                        {
                            foreach (var p in (JsonArray)JsonNode.Parse(json)!["properties"]!)
                                if (p!["id"]?.GetValue<string>() == id && p["value"]?.GetValue<double>() is double v) return v;
                            return double.NaN;
                        }
                        run.Density = Value("density");
                        run.Rg = Value("rg");
                        if (tg)
                        {   // a copy is cooled step by step (the document keeps its state); Tg from the two-line fit of v(T)
                            Post(() => run.Detail = "Tg scan");
                            json = doc.Analyze("tg", ao, (_, _) => !token.IsCancellationRequested);
                            run.Tg = Value("tg");
                        }
                        doc.Save(run.Path);
                    }
                    run.FinalDetail = $"ρ {run.Density:0.000} g/cm³ · Rg {run.Rg:0.00} Å" + (double.IsFinite(run.Tg) ? $" · Tg {run.Tg:0} K" : "");
                    run.FinalStatus = "done";
                    Post(() => { run.Status = "done"; run.Detail = run.FinalDetail; RaiseSweep(); });
                }
                catch (Exception e)
                {
                    var stopped = token.IsCancellationRequested;
                    run.FinalDetail = stopped ? "" : e.Message;
                    run.FinalStatus = stopped ? "stopped" : "failed";
                    Post(() => { run.Status = run.FinalStatus ?? run.Status; run.Detail = run.FinalDetail; RaiseSweep(); });
                }
                finally { gate.Release(); }
            }));
        }
        await Task.WhenAll(tasks);
        foreach (var r in _sweepRuns.Where(r => r.FinalStatus != null)) { r.Status = r.FinalStatus!; r.Detail = r.FinalDetail; r.FinalStatus = null; }
        WriteSweepResults(poly);
        SweepRunning = false;
        RaiseSweep();
        Status = $"Sweep finished · {SweepSummary} · {RecentFiles.Tilde(_swFolder)}";
    }

    /// <summary>results.json beside the cells (caps-sweep/1, as caps.sweep in Python writes and reads it).</summary>
    private void WriteSweepResults(FilmPolymer poly)
    {
        try
        {
            var units = new JsonArray();
            foreach (var u in (JsonNode.Parse(poly.Spec)?["units"] as JsonArray) ?? new JsonArray())
                if (u?["smiles"]?.GetValue<string>() is string smi) units.Add(smi);
            var runs = new JsonArray();
            foreach (var r in _sweepRuns)
            {
                var props = new JsonObject();
                if (double.IsFinite(r.Density)) props["density"] = new JsonObject { ["value"] = r.Density, ["unit"] = "g/cm³", ["name"] = "Density" };
                if (double.IsFinite(r.Rg)) props["rg"] = new JsonObject { ["value"] = r.Rg, ["unit"] = "Å", ["name"] = "Radius of gyration" };
                if (double.IsFinite(r.Tg)) props["tg"] = new JsonObject { ["value"] = r.Tg, ["unit"] = "K", ["name"] = "Glass transition temperature" };
                runs.Add(new JsonObject
                {
                    ["tacticity"] = r.Tacticity, ["dp"] = r.Dp, ["seed"] = r.Seed, ["file"] = System.IO.Path.GetFileName(r.Path),
                    ["status"] = r.Status, ["properties"] = props,
                });
            }
            File.WriteAllText(System.IO.Path.Combine(_swFolder, "results.json"),
                              new JsonObject { ["format"] = "caps-sweep/1", ["units"] = units, ["runs"] = runs }.ToJsonString(new System.Text.Json.JsonSerializerOptions { WriteIndented = true }));
        }
        catch (Exception e) { SweepError = "results.json: " + e.Message; }
    }

    private static void Post(Action a) => Avalonia.Threading.Dispatcher.UIThread.Post(a);

    public void PauseSweep() { _swPaused = !_swPaused; Raise(nameof(SweepPauseText)); }
    public void StopSweep() => _swCancel?.Cancel();

    /// <summary>More seeds for every condition (the next free numbers).</summary>
    public void AddSweepSeeds()
    {
        if (_swRunning) return;
        var seeds = ParseInts(_swSeeds);
        var next = (seeds.Count == 0 ? 0 : seeds.Max()) + 1;
        var kept = _sweepRuns.ToList();
        SweepSeeds = string.Join(", ", seeds.Append(next));
        foreach (var r in _sweepRuns)   // keep what already ran
            if (kept.FirstOrDefault(k => k.Condition == r.Condition && k.Seed == r.Seed) is { } old)
            {
                r.Status = old.Status;
                r.Detail = old.Detail;
                r.Density = old.Density;
                r.Rg = old.Rg;
            }
        RaiseSweep();
    }

    public void OpenSweepRun(SweepRun r)
    {
        if (r.Status == "done" && File.Exists(r.Path)) { SetModule(8); Open(r.Path, null); }
        else Status = r.Tip;
    }
}
