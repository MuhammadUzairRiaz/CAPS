using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One barrier the literature gives for an exchange chemistry, as the picker lists it.</summary>
public sealed record ExchangeBarrierRow(string Value, string What, string Source, bool Read)
{
    public string Status => Read ? "read" : "quoted";
}

/// <summary>An exchange chemistry of data/exchange/library.json.</summary>
public sealed record ExchangeChemistryRow(string Id, string Name, string Family, string Mechanism, string Move, string Template, double? Barrier,
                                          string BarrierWhy, double? Capture, string CaptureNote, ExchangeBarrierRow[] Barriers, string Notes)
{
    public override string ToString() => $"{Name} · {Family}";
}

/// <summary>Modules › Bond exchange (core exchange_moves / exchange_protocols; caps exchange): vitrimers and self-healing
/// networks. A chemistry from the library (its move, barrier and capture, every barrier with its source), moves during
/// CAPS's own MD with the accepted moves counted as they come, the protocols (step strain with G(t), Tv scan, heal run),
/// the cuts a healing test starts from, and the healing efficiency of a healed against an intact pull.</summary>
public sealed partial class MainViewModel
{
    public bool IsExchange => _module == 83;
    public void OpenExchange() { SetModule(83); if (ExChemistries.Count == 0) LoadExchangeLibrary(); }

    private bool _exRunning;
    public bool ExRunning { get => _exRunning; private set { if (Set(ref _exRunning, value)) RaiseBusy(); } }
    private CancellationTokenSource? _exCts;
    public void CancelExchange() => _exCts?.Cancel();

    // ------------------------------------------------------------------ the chemistry
    public ObservableCollection<ExchangeChemistryRow> ExChemistries { get; } = new();
    private ExchangeChemistryRow? _exChem;
    public ExchangeChemistryRow? ExChemistry
    {
        get => _exChem;
        set
        {
            if (!Set(ref _exChem, value) || value == null) return;
            // the chemistry's move, barrier and capture: the user may change each afterwards
            ExKind = value.Template.Length > 0 ? 3 : value.Move switch { "shift" => 1, "bind" => 2, _ => 0 };
            ExBarrier = (decimal)(value.Barrier ?? 0);
            if (value.Capture is { } c) ExCapture = (decimal)c;
            Raise(nameof(ExChemHasBarrier)); Raise(nameof(ExChemSummary));
        }
    }
    public bool ExChemHasBarrier => _exChem?.Barrier != null;
    public string ExChemSummary => _exChem == null ? "" : _exChem.Mechanism + " · " + (_exChem.Template.Length > 0 ? "template " + _exChem.Template + " (all-atom, relaxed)" : _exChem.Move + " moves")
        + (_exChem.Barrier is { } b ? string.Format(CultureInfo.InvariantCulture, " · default barrier {0:0.##} kcal/mol", b) : " · no barrier read yet");

    public void LoadExchangeLibrary()
    {
        ExChemistries.Clear();
        try
        {
            var path = Paths.Data is { } d ? Path.Combine(d, "exchange", "library.json") : "";
            var r = JsonNode.Parse(CapsDocument.ExchangeLibrary(path))!;
            if (r["ok"]?.GetValue<bool>() != true) { ExError = (string?)r["error"] ?? "the exchange library could not be read"; return; }
            foreach (var c in (JsonArray)r["chemistries"]!)
            {
                double? N(JsonNode? x) => x == null ? null : x.GetValue<double>();
                var bars = ((JsonArray)c!["barriers"]!).Select(b => new ExchangeBarrierRow(
                    $"{(string?)b!["value"]} {(string?)b["unit"]}", string.Join(" · ", new[] { (string?)b["kind"], (string?)b["method"], (string?)b["conditions"] }.Where(s => !string.IsNullOrEmpty(s))),
                    "doi:" + (string?)b["source"], (string?)b["status"] == "read")).ToArray();
                ExChemistries.Add(new ExchangeChemistryRow((string)c["id"]!, (string?)c["name"] ?? "", (string?)c["family"] ?? "", (string?)c["mechanism"] ?? "",
                    (string?)c["move"] ?? "", (string?)c["template"] ?? "", N(c["default_kcal"]), (string?)c["default_why"] ?? "", N(c["capture"]),
                    c["capture"] == null ? "" : $"{(string?)c["capture_pair"]} · doi:{(string?)c["capture_source"]}", bars, (string?)c["notes"] ?? ""));
            }
            ExChemistry = ExChemistries.FirstOrDefault(c => c.Id == "transesterification") ?? ExChemistries.FirstOrDefault();
        }
        catch (Exception e) { ExError = e.Message; }
    }

    // ------------------------------------------------------------------ the moves
    public static readonly string[] ExKinds = ["Swap partners", "Shift (valences)", "Bind / break", "Template (all-atom)"];
    private static readonly string[] ExKindIds = ["swap", "shift", "bind", "template"];
    private int _exKind;
    private decimal _exCapture = 1.3m, _exBarrier, _exBinding = 10, _exValence = 2, _exAttempts = 10, _exEvery = 10, _exRelax = 1000;
    private bool _exUseSelection, _exNoMoves;
    public int ExKind { get => _exKind; set { if (Set(ref _exKind, Math.Clamp(value, 0, 3))) { Raise(nameof(ExIsBind)); Raise(nameof(ExIsTemplate)); } } }
    public bool ExIsBind => _exKind is 1 or 2;
    public bool ExIsTemplate => _exKind == 3;
    public decimal ExCapture { get => _exCapture; set => Set(ref _exCapture, Math.Clamp(value, 0.5m, 10m)); }
    public decimal ExBarrier { get => _exBarrier; set => Set(ref _exBarrier, Math.Clamp(value, 0m, 100m)); }
    public decimal ExBinding { get => _exBinding; set => Set(ref _exBinding, Math.Clamp(value, -100m, 100m)); }
    public decimal ExValence { get => _exValence; set => Set(ref _exValence, Math.Clamp(Math.Round(value), 1, 8)); }
    public decimal ExAttempts { get => _exAttempts; set => Set(ref _exAttempts, Math.Clamp(Math.Round(value), 1, 10000)); }
    public decimal ExEvery { get => _exEvery; set => Set(ref _exEvery, Math.Clamp(Math.Round(value), 1, 100000)); }
    public decimal ExRelax { get => _exRelax; set => Set(ref _exRelax, Math.Clamp(Math.Round(value), 0, 100000)); }
    public bool ExUseSelection { get => _exUseSelection; set => Set(ref _exUseSelection, value); }
    /// <summary>The network as it is (no moves): the reference a protocol's exchange is compared with.</summary>
    public bool ExNoMoves { get => _exNoMoves; set => Set(ref _exNoMoves, value); }

    // ------------------------------------------------------------------ the protocol
    public static readonly string[] ExProtocols = ["Exchange during MD", "Step strain · G(t)", "Tv scan", "Heal a cut"];
    private static readonly string[] ExProtocolIds = ["run", "step_strain", "tv", "heal_run"];
    public static readonly string[] ExAxes = ["x", "y", "z"];
    private int _exProtocol, _exAxis = 2;
    private decimal _exT = 450, _exDt = 1, _exSteps = 2000, _exStrain = 0.05m, _exPs = 2, _exTFrom = 300, _exTTo = 500, _exTStep = 50, _exPsPerT = 1, _exChunks = 10, _exPlane = 0.5m;
    private bool _exNpt;
    public int ExProtocol { get => _exProtocol; set { if (Set(ref _exProtocol, Math.Clamp(value, 0, 3))) { Raise(nameof(ExIsRun)); Raise(nameof(ExIsStrain)); Raise(nameof(ExIsTv)); Raise(nameof(ExIsHeal)); Raise(nameof(ExStepsHere)); } } }
    public bool ExIsRun => _exProtocol == 0;
    public bool ExIsStrain => _exProtocol == 1;
    public bool ExIsTv => _exProtocol == 2;
    public bool ExIsHeal => _exProtocol == 3;
    public int ExAxis { get => _exAxis; set => Set(ref _exAxis, Math.Clamp(value, 0, 2)); }
    public decimal ExT { get => _exT; set => Set(ref _exT, Math.Clamp(value, 1, 5000)); }
    public decimal ExDt { get => _exDt; set { if (Set(ref _exDt, Math.Clamp(value, 0.01m, 10))) Raise(nameof(ExStepsHere)); } }
    public decimal ExSteps { get => _exSteps; set { if (Set(ref _exSteps, Math.Clamp(Math.Round(value), 10, 100_000_000))) Raise(nameof(ExStepsHere)); } }
    public decimal ExStrain { get => _exStrain; set => Set(ref _exStrain, Math.Clamp(value, -0.5m, 2m)); }
    public decimal ExPs { get => _exPs; set { if (Set(ref _exPs, Math.Clamp(value, 0.01m, 100000))) Raise(nameof(ExStepsHere)); } }
    public decimal ExTFrom { get => _exTFrom; set { if (Set(ref _exTFrom, Math.Clamp(value, 1, 5000))) Raise(nameof(ExStepsHere)); } }
    public decimal ExTTo { get => _exTTo; set { if (Set(ref _exTTo, Math.Clamp(value, 1, 5000))) Raise(nameof(ExStepsHere)); } }
    public decimal ExTStep { get => _exTStep; set { if (Set(ref _exTStep, Math.Clamp(value, 1, 2000))) Raise(nameof(ExStepsHere)); } }
    public decimal ExPsPerT { get => _exPsPerT; set { if (Set(ref _exPsPerT, Math.Clamp(value, 0.01m, 100000))) Raise(nameof(ExStepsHere)); } }
    public decimal ExChunks { get => _exChunks; set => Set(ref _exChunks, Math.Clamp(Math.Round(value), 1, 1000)); }
    public decimal ExPlane { get => _exPlane; set => Set(ref _exPlane, Math.Clamp(value, 0.01m, 0.99m)); }
    public bool ExNpt { get => _exNpt; set => Set(ref _exNpt, value); }
    /// <summary>The MD steps the protocol takes, for the local-run guide.</summary>
    public string ExStepsHere
    {
        get
        {
            var ps = _exProtocol switch
            {
                0 => _exSteps * _exDt / 1000,
                2 => (Math.Floor(Math.Abs(_exTTo - _exTFrom) / _exTStep) + 1) * _exPsPerT,
                _ => _exPs,
            };
            return string.Format(CultureInfo.InvariantCulture, "{0:N0} MD steps ({1:0.###} ps)", Math.Round(ps * 1000 / _exDt), ps);
        }
    }

    // ------------------------------------------------------------------ the live counter and the results
    private long _exAccepted, _exTried;
    private string _exError = "", _exProgress = "", _exResult = "", _exFieldAfter = "", _exLearned = "—";
    public long ExAccepted { get => _exAccepted; private set { if (Set(ref _exAccepted, value)) Raise(nameof(ExRatio)); } }
    public long ExTried { get => _exTried; private set { if (Set(ref _exTried, value)) Raise(nameof(ExRatio)); } }
    public string ExRatio => _exTried > 0 ? (100.0 * _exAccepted / _exTried).ToString("0.0", CultureInfo.InvariantCulture) + " %" : "—";
    public string ExLearned { get => _exLearned; private set => Set(ref _exLearned, value); }
    public string ExError { get => _exError; private set { if (Set(ref _exError, value)) Raise(nameof(ExHasError)); } }
    public bool ExHasError => _exError.Length > 0;
    public string ExProgress { get => _exProgress; private set => Set(ref _exProgress, value); }
    public string ExResult { get => _exResult; private set => Set(ref _exResult, value); }
    public string ExFieldAfter { get => _exFieldAfter; private set => Set(ref _exFieldAfter, value); }
    /// <summary>Accepted moves against time, as they come (t in ps).</summary>
    public List<(double X, double Y)> ExLive { get; } = new();
    /// <summary>The protocol's own curve: G(t), v(T) with the rate, or bridges against time.</summary>
    public (double X, double Y)[] ExCurve { get; private set; } = [];
    public (double X, double Y)[] ExCurve2 { get; private set; } = [];
    public string ExCurveKind { get; private set; } = "";
    public event Action? ExChanged;

    private string ExchangeJson()
    {
        var o = new JsonObject
        {
            ["protocol"] = ExProtocolIds[_exProtocol], ["exchange"] = !_exNoMoves, ["kind"] = ExKindIds[_exKind],
            ["capture"] = (double)_exCapture, ["barrier"] = (double)_exBarrier, ["binding"] = (double)_exBinding,
            ["valence_a"] = (int)_exValence, ["valence_b"] = (int)_exValence, ["attempts"] = (int)_exAttempts, ["every"] = (int)_exEvery,
            ["relax"] = (int)_exRelax, ["seed"] = ExSeed.Take(), ["temperature"] = (double)_exT, ["dt"] = (double)_exDt,
            ["steps"] = (long)_exSteps, ["axis"] = ExAxes[_exAxis], ["strain"] = (double)_exStrain, ["ps"] = (double)_exPs,
            ["t_from"] = (double)_exTFrom, ["t_to"] = (double)_exTTo, ["t_step"] = (double)_exTStep, ["ps_per_step"] = (double)_exPsPerT,
            ["chunks"] = (int)_exChunks, ["plane"] = (double)_exPlane, ["npt"] = _exNpt,
        };
        if (_exChem != null)
        {
            o["chemistry"] = _exChem.Id;
            o["library"] = Paths.Data is { } d ? Path.Combine(d, "exchange", "library.json") : "";
            o["reactions"] = Paths.Reactions ?? "";
        }
        if (_exUseSelection && _selection.Count > 0) o["sites"] = new JsonArray(_selection.Select(i => (JsonNode)i).ToArray());
        return o.ToJsonString();
    }

    private SeedChoice? _exSeedC;
    public SeedChoice ExSeed => _exSeedC ??= new SeedChoice();

    public async Task RunExchange()
    {
        if (!Idle || _doc == null) { if (_doc == null) Status = "Open a structure first"; return; }
        var doc = _doc;
        var json = ExchangeJson();
        var inv = CultureInfo.InvariantCulture;
        ExRunning = true;
        ExError = "";
        ExAccepted = ExTried = 0;
        ExLive.Clear();
        _exCts = new CancellationTokenSource();
        var ct = _exCts.Token;
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var total = (double)(_exProtocol == 0 ? _exSteps * _exDt / 1000 : _exPs);
        Status = "Bond exchange: running…";
        try
        {
            var res = await Task.Run(() => doc.ExchangeRun(json, (stage, f) =>
            {
                // "A accepted of T tried": the counter and its curve follow every report
                var parts = stage.Split(' ');
                long.TryParse(parts.ElementAtOrDefault(0), NumberStyles.Integer, inv, out var acc);
                long.TryParse(parts.ElementAtOrDefault(3), NumberStyles.Integer, inv, out var tried);
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    if (!ExRunning) return;
                    ExAccepted = acc; ExTried = tried;
                    ExLive.Add((f * total, acc));
                    ExProgress = $"{f * 100:0} % · {acc} accepted of {tried} tried · {sw.Elapsed.TotalSeconds:0} s";
                    ExChanged?.Invoke();
                });
                return !ct.IsCancellationRequested;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { ExError = (string?)r["error"] ?? "the run failed"; Status = "Bond exchange: " + ExError; return; }
            var st = r["stats"]!;
            double S(string k) => st[k]?.GetValue<double>() ?? 0;
            ExAccepted = (long)S("accepted"); ExTried = (long)S("attempts");
            ExLearned = $"{S("learned"):0} learned · {S("new_types"):0} new types";
            double[] A(JsonNode? a) => a is JsonArray v ? v.Select(x => x?.GetValue<double>() ?? double.NaN).ToArray() : [];
            var proto = (string?)r["protocol"] ?? "run";
            ExCurveKind = proto;
            ExCurve2 = [];
            var lines = new List<string> { (string?)r["summary"] ?? "" };
            switch (proto)
            {
                case "run":
                    ExCurve = ((JsonArray)r["swaps"]!).Select(p => (p![0]!.GetValue<double>(), p[1]!.GetValue<double>())).ToArray();
                    break;
                case "step_strain":
                    ExCurve = A(r["t"]).Zip(A(r["g"])).Where(p => p.First > 0).Select(p => (p.First, p.Second)).ToArray();
                    lines.Add(r["kww_mean"] is null ? "the stress did not decay enough for a KWW fit: hold longer, or raise the temperature"
                        : string.Format(inv, "KWW: τ {0:0.###} ps, β {1:0.###}, mean relaxation time {2:0.###} ps", r["kww_tau"]!.GetValue<double>(), r["kww_beta"]!.GetValue<double>(), r["kww_mean"]!.GetValue<double>()));
                    break;
                case "tv":
                    ExCurve = A(r["T"]).Zip(A(r["v"])).Select(p => (p.First, p.Second)).ToArray();
                    ExCurve2 = A(r["T"]).Zip(A(r["rate"])).Select(p => (p.First, p.Second)).ToArray();
                    lines.Add(r["tv"] is null ? "no hinge in the specific volume: " + ((string?)r["tv_note"] ?? "")
                        : string.Format(inv, "hinge of the specific volume at {0:0} ± {1:0} K (Tv with exchange, Tg without)", r["tv"]!.GetValue<double>(), r["tv_err"]!.GetValue<double>()));
                    if (r["ea"] is not null) lines.Add(string.Format(inv, "exchange rate: Arrhenius E_a {0:0.0} ± {1:0.0} kJ/mol", r["ea"]!.GetValue<double>(), r["ea_err"]!.GetValue<double>()));
                    break;
                default:
                    ExCurve = A(r["t"]).Zip(A(r["bridges"])).Select(p => (p.First, p.Second)).ToArray();
                    ExCurve2 = A(r["t"]).Zip(A(r["crossing"])).Select(p => (p.First, p.Second)).ToArray();
                    lines.Add(ExCurve.Length > 0 ? string.Format(inv, "{0:0} bonds bridge the cut at the end", ExCurve[^1].Y) : "");
                    break;
            }
            if (!string.IsNullOrEmpty((string?)st["note"])) lines.Add((string)st["note"]!);
            ExResult = string.Join("\n", lines.Where(l => l.Length > 0));
            ExFieldAfter = "Force field after the run: " + ((string?)r["field_after"] ?? "");
            AfterRun(doc, " · exchanged");
            Status = $"Bond exchange: {ExAccepted} accepted of {ExTried} tried";
            ExChanged?.Invoke();
        }
        catch (Exception e) { ExError = e.Message; Status = "Bond exchange: " + e.Message; }
        finally { ExRunning = false; ExProgress = ""; }
    }

    // ------------------------------------------------------------------ a healing test's start
    public static readonly string[] ExCuts = ["Cut", "Gap (weld)", "Crack (void)", "Fragments (recycle)"];
    private static readonly string[] ExCutIds = ["cut", "gap", "crack", "fragment"];
    private int _exCut;
    private decimal _exGap = 3, _exRadius = 4;
    public int ExCut { get => _exCut; set { if (Set(ref _exCut, Math.Clamp(value, 0, 3))) { Raise(nameof(ExCutIsGap)); Raise(nameof(ExCutIsCrack)); Raise(nameof(ExCutHasPlane)); } } }
    public bool ExCutIsGap => _exCut == 1;
    public bool ExCutIsCrack => _exCut == 2;
    public bool ExCutHasPlane => _exCut is 0 or 1;
    public decimal ExGap { get => _exGap; set => Set(ref _exGap, Math.Clamp(value, 0.1m, 100)); }
    public decimal ExRadius { get => _exRadius; set => Set(ref _exRadius, Math.Clamp(value, 0.5m, 100)); }

    public void BuildExchangeStart()
    {
        if (!Idle || _doc == null) { if (_doc == null) Status = "Open a structure first"; return; }
        ExError = "";
        try
        {
            var json = new JsonObject
            {
                ["op"] = ExCutIds[_exCut], ["axis"] = ExAxes[_exAxis], ["plane"] = (double)_exPlane, ["gap"] = (double)_exGap, ["radius"] = (double)_exRadius, ["axes"] = "xyz",
            }.ToJsonString();
            var r = JsonNode.Parse(_doc.ExchangeBuild(json))!;
            if (r["ok"]?.GetValue<bool>() != true) { ExError = (string?)r["error"] ?? "could not cut"; return; }
            ExFieldAfter = "Force field after the cut: " + ((string?)r["field_after"] ?? "");
            AfterRun(_doc, " · cut");
            ExProtocol = _exCut == 3 ? 0 : 3;   // a cut heals; fragments rejoin by exchange
            Status = "Healing test: " + ((string?)r["summary"] ?? "");
        }
        catch (Exception e) { ExError = e.Message; }
    }

    // ------------------------------------------------------------------ healing efficiency
    private string _exHealed = "", _exIntact = "", _exHealText = "";
    public string ExHealed { get => _exHealed; set { if (Set(ref _exHealed, value ?? "")) Raise(nameof(ExHealedName)); } }
    public string ExIntact { get => _exIntact; set { if (Set(ref _exIntact, value ?? "")) Raise(nameof(ExIntactName)); } }
    public string ExHealedName => _exHealed.Length > 0 ? Path.GetFileName(_exHealed) : "healed pull (strain, stress)";
    public string ExIntactName => _exIntact.Length > 0 ? Path.GetFileName(_exIntact) : "intact pull (strain, stress)";
    public string ExHealText { get => _exHealText; private set => Set(ref _exHealText, value); }
    public (double X, double Y)[] ExHealCurve { get; private set; } = [];
    public (double X, double Y)[] ExIntactCurve { get; private set; } = [];

    /// <summary>A stress–strain table (two columns of numbers; lines of text skipped), as caps tensile writes it.</summary>
    private static (double X, double Y)[] ReadCurve(string path)
    {
        var inv = CultureInfo.InvariantCulture;
        var pts = new List<(double, double)>();
        foreach (var line in File.ReadLines(path))
        {
            var f = line.Split([',', ' ', '\t', ';'], StringSplitOptions.RemoveEmptyEntries);
            if (f.Length >= 2 && double.TryParse(f[0], NumberStyles.Float, inv, out var x) && double.TryParse(f[1], NumberStyles.Float, inv, out var y)) pts.Add((x, y));
        }
        return pts.ToArray();
    }

    public void CompareHealing()
    {
        ExError = "";
        if (_exHealed.Length == 0 || _exIntact.Length == 0) { ExError = "Choose the healed and the intact stress–strain tables"; return; }
        try
        {
            ExHealCurve = ReadCurve(_exHealed);
            ExIntactCurve = ReadCurve(_exIntact);
            JsonObject C((double X, double Y)[] c) => new() { ["strain"] = new JsonArray(c.Select(p => (JsonNode)p.X).ToArray()), ["stress"] = new JsonArray(c.Select(p => (JsonNode)p.Y).ToArray()) };
            var r = JsonNode.Parse(CapsDocument.ExchangeAnalyze(new JsonObject { ["op"] = "heal", ["healed"] = C(ExHealCurve), ["intact"] = C(ExIntactCurve) }.ToJsonString()))!;
            if (r["ok"]?.GetValue<bool>() != true) { ExError = (string?)r["error"] ?? "could not compare"; return; }
            var p = r["properties"]![0]!;
            var lines = new List<string> { string.Format(CultureInfo.InvariantCulture, "toughness ratio {0:0.#} %", p["value"]!.GetValue<double>()) };
            if (p["extra"] is JsonObject ex)
                foreach (var (k, v) in ex) lines.Add(string.Format(CultureInfo.InvariantCulture, "{0} {1:0.###}", k, v!.GetValue<double>()));
            ExHealText = string.Join("\n", lines);
            ExChanged?.Invoke();
        }
        catch (Exception e) { ExError = e.Message; }
    }

    /// <summary>The page's run as the caps command line runs it (caps exchange …, the same engine and options).</summary>
    public string ExchangeCommand()
    {
        var inv = CultureInfo.InvariantCulture;
        var head = _exProtocol switch
        {
            0 => $"caps exchange run IN.data -o OUT.data --steps {_exSteps:0}",
            1 => string.Format(inv, "caps exchange step-strain IN.data -o OUT.data --axis {0} --strain {1} --ps {2}", ExAxes[_exAxis], _exStrain, _exPs),
            2 => string.Format(inv, "caps exchange tv IN.data -o OUT.data --from {0} --to {1} --step {2} --ps {3}", _exTFrom, _exTTo, _exTStep, _exPsPerT),
            _ => string.Format(inv, "caps exchange heal-run IN.data -o OUT.data --axis {0} --plane {1} --ps {2} --chunks {3}{4}", ExAxes[_exAxis], _exPlane, _exPs, _exChunks, _exNpt ? " --npt" : ""),
        };
        var moves = _exNoMoves ? " --no-exchange" : string.Format(inv, "{0} --kind {1} --capture {2} --barrier {3} --attempts {4} --every {5}{6}{7}",
            _exChem != null ? " --chemistry " + _exChem.Id : "", ExKindIds[_exKind], _exCapture, _exBarrier, _exAttempts, _exEvery,
            ExIsBind ? string.Format(inv, " --binding {0} --valence {1}", _exBinding, _exValence) : "", ExIsTemplate ? $" --relax {_exRelax:0}" : "");
        return head + moves + string.Format(inv, " --temp {0} --dt {1} --ff FF.json", _exT, _exDt);
    }
}
