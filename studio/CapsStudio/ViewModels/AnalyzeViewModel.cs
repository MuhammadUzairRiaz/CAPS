using System.Collections.ObjectModel;
using System.Globalization;
using System.Text;
using System.Text.Json;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One calculation that can be switched on (a property id of caps analyze), or one that is not built yet.</summary>
public sealed class CalcChip : ObservableObject
{
    public string Id { get; init; } = "";
    public string Label { get; init; } = "";
    public string Tip { get; init; } = "";
    public bool Available { get; init; } = true;
    private bool _on;
    public bool IsOn { get => _on; set => Set(ref _on, value && Available); }
    private bool _active;
    /// <summary>Shown on the focused page it belongs to (Mechanics, Scattering, Free volume).</summary>
    public bool Active { get => _active; set => Set(ref _active, value); }
}

public sealed record CalcGroup(string Name, CalcChip[] Chips);

/// <summary>An experimental or literature range for one property of one material (data/reference/polymers.json).</summary>
public sealed record RefValue(double Lo, double Hi, string Unit, string Condition, string Source);

public sealed record RefMaterial(string Id, string Name, Dictionary<string, RefValue> Values)
{
    public override string ToString() => Name;
}

/// <summary>A property result as a card: value ± error, method, notes, further numbers and the comparison with experiment.</summary>
public sealed class ResultCard
{
    public string Id { get; init; } = "";
    public string Name { get; init; } = "";
    public double Value { get; init; } = double.NaN;
    public double Error { get; init; } = double.NaN;
    public string Unit { get; init; } = "";
    public string Method { get; init; } = "";
    public string[] Notes { get; init; } = [];
    public (string Key, double Value)[] Extra { get; init; } = [];
    public RefValue? Ref { get; set; }

    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    public bool HasValue => !double.IsNaN(Value);
    private bool HasError => HasValue && !double.IsNaN(Error) && Error > 0;
    /// <summary>Decimals that show the error to two significant figures (the value is rounded to match).</summary>
    private int ErrorDecimals => Math.Clamp(1 - (int)Math.Floor(Math.Log10(Error)), 0, 8);
    // shown in the display units chosen in Settings › Units (the stored value keeps CAPS's units)
    private double F => DisplayUnits.For(Unit).Factor;
    private int ShownDecimals => Math.Clamp(1 - (int)Math.Floor(Math.Log10(Error * F)), 0, 8);
    public string ValueText => !HasValue ? "—" : HasError ? (Value * F).ToString("F" + ShownDecimals, Inv) : Num(Value * F);
    public string ErrorText => HasError ? "± " + (Error * F).ToString("F" + ShownDecimals, Inv) : "";
    public string UnitText => DisplayUnits.For(Unit).Unit.Replace("^½", "½");
    public string NotesText => string.Join("\n", Notes);
    public bool HasNotes => Notes.Length > 0;
    public string ExtraText => string.Join("\n", Extra.Select(e => $"{e.Key}: {Num(e.Value)}"));
    public bool HasExtra => Extra.Length > 0;
    public bool HasRef => Ref != null;
    public string RefText => Ref == null ? "" :
        $"{(Ref.Source.StartsWith("Flory", StringComparison.Ordinal) ? "lit." : "exp.")} {(Ref.Lo == Ref.Hi ? Num(Ref.Lo) : Num(Ref.Lo) + "–" + Num(Ref.Hi))} {Ref.Unit.Replace("^½", "½")}".TrimEnd() +
        (Ref.Condition.Length > 0 ? " · " + Ref.Condition : "");
    public string RefSource => Ref?.Source ?? "";
    /// <summary>Within the range (allowing the error), outside it, or not compared.</summary>
    public string Verdict
    {
        get
        {
            if (Ref == null || !HasValue) return "";
            var e = double.IsNaN(Error) ? 0 : Error;
            var uncertain = e > 0.3 * Math.Abs(Value);
            if (Value + e >= Ref.Lo && Value - e <= Ref.Hi)
                return uncertain ? string.Format(Inv, "overlaps the range, but the error is {0:F0} % of the value", 100 * e / Math.Abs(Value)) : "within the range";
            var d = Value < Ref.Lo ? (Ref.Lo - Value) / Math.Max(Math.Abs(Ref.Lo), 1e-12) : (Value - Ref.Hi) / Math.Max(Math.Abs(Ref.Hi), 1e-12);
            return string.Format(Inv, "{0:F1} % {1} the range", 100 * d, Value < Ref.Lo ? "below" : "above");
        }
    }
    public IBrush VerdictBrush => Verdict.StartsWith("within", StringComparison.Ordinal) ? Good : Warn;   // "overlaps … error" stays amber
    internal static IBrush Good => Tokens.Brush("OkB");
    internal static IBrush Warn => Tokens.Brush("WarnB");

    public static string Num(double v)
    {
        if (double.IsNaN(v)) return "—";
        var a = Math.Abs(v);
        if (a == 0) return "0";
        if (a >= 1e5 || a < 1e-3) return v.ToString("0.###e0", Inv);
        var digits = Math.Max(0, 3 - (int)Math.Floor(Math.Log10(a)));   // four significant figures
        return Math.Round(v, Math.Min(digits, 6)).ToString("0.######", Inv);
    }
}

/// <summary>A curve behind a result.</summary>
public sealed record SeriesItem(string Property, string Label, string XLabel, string YLabel, double[] X, double[] Y, bool LogLog, double? RefY)
{
    /// <summary>Draw the data as points (with Overlay as the line through them: a fit or a smoothed curve).</summary>
    public bool Markers { get; init; }
    public double[]? OverlayX { get; init; }
    public double[]? OverlayY { get; init; }
    public override string ToString() => $"{Property} · {Label}";
}

/// <summary>Analyze › Properties: choose calculations, frames and groups, run them on the open trajectory, and read the
/// results as cards compared with experiment, with the curves behind them. Export CSV and a LaTeX table.</summary>
public sealed class AnalyzeViewModel : ObservableObject
{
    private readonly Func<CapsDocument?> _doc;
    private readonly Action<string> _status;
    private readonly Action<bool> _running;
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public AnalyzeViewModel(Func<CapsDocument?> doc, Action<string> status, Action<bool> running)
    {
        _doc = doc;
        _status = status;
        _running = running;
        Groups =
        [
            new("Structure", [Chip("density", "Density", on: true), Chip("rdf", "RDF", on: true), Chip("sq", "S(q)"), Chip("xray", "X-ray"), Chip("neutron", "Neutron")]),
            new("Chains", [Chip("rg", "Rg", on: true), Chip("ree", "Ree"), Chip("cn", "Cn, C∞"), Chip("persistence", "Persistence"), Chip("orientation", "Orientation"),
                Soon("Entanglements", "Primitive-path analysis (Z1-type) is not built yet")]),
            new("Thermo", [Chip("ced", "CED"), Chip("delta", "δ"), TgChip]),
            new("Mechanics", [StrainChip, FluctChip, TensileChip]),
            new("Dynamics", [Chip("msd", "MSD"), Chip("diffusion", "D"), Chip("relaxation", "Relaxation")]),
            new("Free volume", [Chip("ffv", "Probe insertion"), Chip("psd", "Pore size")]),
            new("Interface", [Chip("zprofile", "z profile"), Chip("adhesion", "Adhesion"), PullShearChip, PullNormalChip]),
            new("Rubber network", [Chip("crosslinks", "Crosslink density")]),
        ];
        LoadReferences();
        PullShearChip.PropertyChanged += (_, _) => Raise(nameof(PullOn));
        PullNormalChip.PropertyChanged += (_, _) => Raise(nameof(PullOn));
    }

    // protocols (their settings show when switched on)
    public CalcChip TgChip { get; } = Chip("tg", "Tg", tip: "Glass transition from a stepwise NPT cooling run of the current frame (a copy: the document is not changed)");
    public CalcChip StrainChip { get; } = Chip("cij_strain", "Cij strain", tip: "Static elastic constants: minimise, strain ±ε in each direction, re-minimise (Theodorou & Suter)");
    public CalcChip FluctChip { get; } = Chip("cij_run", "Cij fluct.", tip: "Elastic constants from stress fluctuations: an NVT run of the current frame, the stress sampled at every step (Lutsko; Clavier et al.)");
    public CalcChip PullShearChip { get; } = Chip("pull_shear", "Pull · shear", tip: "Steered MD: the film dragged along x over the held surface (molecule 1); interfacial shear strength and work");
    public CalcChip PullNormalChip { get; } = Chip("pull_normal", "Pull · normal", tip: "Steered MD: the film pulled off the held surface along +z (needs vacuum above the film); peak normal stress and work of separation");
    public CalcChip TensileChip { get; } = Chip("tensile", "Stress–strain", tip: "Uniaxial deformation MD of the current frame: modulus, Poisson ratio, yield");

    private static CalcChip Chip(string id, string label, bool on = false, string tip = "") => new() { Id = id, Label = label, IsOn = on, Tip = tip.Length > 0 ? tip : Tips.GetValueOrDefault(id, "") };
    private static CalcChip Soon(string label, string tip) => new() { Id = "", Label = label, Tip = tip, Available = false };

    private static readonly Dictionary<string, string> Tips = new()
    {
        ["density"] = "Mass over volume, averaged over the frames",
        ["rdf"] = "g(r) of the chosen pair, averaged over the frames; first peak and coordination number",
        ["sq"] = "Total structure factor: direct reciprocal-lattice sum at low q, g(r) transform above",
        ["xray"] = "X-ray I(q) with Cromer–Mann form factors (Faber–Ziman)",
        ["neutron"] = "Neutron S(q) with coherent scattering lengths (Faber–Ziman)",
        ["rg"] = "Radius of gyration of the chains, √⟨Rg²⟩",
        ["ree"] = "Backbone end-to-end distance, √⟨R²⟩",
        ["cn"] = "Characteristic ratio C_n and C∞ extrapolated in 1/n",
        ["persistence"] = "Persistence length: Flory projection and bond-correlation decay",
        ["ced"] = "Cohesive energy density (E isolated − E bulk)/V with the force field (Field, else GAFF for C and H, UFF otherwise)",
        ["delta"] = "Hildebrand solubility parameter δ = √CED",
        ["cij_fluct"] = "Elastic constants from stress fluctuations of the saved frames of an NVT run",
        ["msd"] = "Mean-square displacement of atoms and molecule centres, all time origins, drift removed",
        ["diffusion"] = "Diffusion coefficient from the linear part of the molecule-centre MSD (Einstein)",
        ["relaxation"] = "End-to-end and segmental (P2) autocorrelations with KWW fits",
        ["ffv"] = "Free volume by probe insertion on a grid: accessible fraction and Bondi FFV",
        ["psd"] = "Pore size distribution: largest atom-free sphere containing each free point",
        ["orientation"] = "Nematic order S of backbone chords, director, Herman's f along z, local crystallinity, and P₂ against height (orientation near a surface)",
        ["crosslinks"] = "Sulfur bridges (mono-, di-, polysulfidic), pendant groups, crosslink density ν and strand mass Mc = ρ/2ν",
        ["zprofile"] = "Mass density along z for the surface (molecule 1) and the film: first-layer peak and the film's own density",
        ["adhesion"] = "Work of adhesion −(E all − E surface − E film)/area between the surface (molecule 1) and the film, with the force field",
    };

    public CalcGroup[] Groups { get; }

    // ---------------------------------------------------------------- source
    private int _first, _last = -1, _stride = 1, _pair;
    private double _framePs, _timestepFs = 1, _probe, _grid = 0.4, _fitFrom = 0.2, _fitTo = 0.5;
    private bool _inter;
    public decimal FirstD { get => _first; set => Set(ref _first, (int)Math.Max(0, value), nameof(FirstD)); }
    /// <summary>Last frame, −1 for the last one.</summary>
    public decimal LastD { get => _last; set => Set(ref _last, (int)Math.Max(-1, value), nameof(LastD)); }
    public decimal StrideD { get => _stride; set => Set(ref _stride, (int)Math.Max(1, value), nameof(StrideD)); }
    /// <summary>Time between frames in ps; 0 takes it from the timesteps in the file × the MD time step.</summary>
    public decimal FramePsD { get => (decimal)_framePs; set => Set(ref _framePs, (double)Math.Max(0, value), nameof(FramePsD)); }
    public decimal TimestepFsD { get => (decimal)_timestepFs; set => Set(ref _timestepFs, (double)Math.Max(0.01m, value), nameof(TimestepFsD)); }
    public int PairIndex { get => _pair; set => Set(ref _pair, value); }
    public bool InterOnly { get => _inter; set => Set(ref _inter, value); }
    private double _qmax = 25, _dq = 0.02, _qDirect = 4;
    private int _deuterate;
    public decimal QmaxD { get => (decimal)_qmax; set => Set(ref _qmax, (double)Math.Clamp(value, 1m, 40m), nameof(QmaxD)); }
    public decimal DqD { get => (decimal)_dq; set => Set(ref _dq, (double)Math.Clamp(value, 0.002m, 0.2m), nameof(DqD)); }
    public decimal QDirectD { get => (decimal)_qDirect; set => Set(ref _qDirect, (double)Math.Clamp(value, 0m, 10m), nameof(QDirectD)); }
    /// <summary>Neutron contrast: 0 none, 1 every H → D, 2 aliphatic H (d-backbone), 3 aromatic H (d-ring), 4 H on O/N.</summary>
    public int Deuterate { get => _deuterate; set => Set(ref _deuterate, Math.Clamp(value, 0, 4)); }
    public decimal ProbeD { get => (decimal)_probe; set => Set(ref _probe, (double)Math.Max(0, value), nameof(ProbeD)); }
    public decimal GridD { get => (decimal)_grid; set => Set(ref _grid, (double)Math.Clamp(value, 0.1m, 2m), nameof(GridD)); }
    public decimal FitFromD { get => (decimal)_fitFrom; set => Set(ref _fitFrom, (double)Math.Clamp(value, 0m, 0.95m), nameof(FitFromD)); }
    public decimal FitToD { get => (decimal)_fitTo; set => Set(ref _fitTo, (double)Math.Clamp(value, 0.05m, 1m), nameof(FitToD)); }
    // protocol settings
    private double _fluctPs = 100, _eqPs = 20;
    /// <summary>Unsampled NPT run before the tensile pull and the first cooling hold (ps).</summary>
    public decimal EqPsD { get => (decimal)_eqPs; set => Set(ref _eqPs, (double)Math.Max(0, value), nameof(EqPsD)); }
    public decimal FluctPsD { get => (decimal)_fluctPs; set => Set(ref _fluctPs, (double)Math.Max(1, value), nameof(FluctPsD)); }
    private double _tgFrom = 500, _tgTo = 200, _tgStep = 20, _tgPs = 100, _cijStrain = 1e-4, _fluctT = 300, _tensRate = 1e-3, _tensMax = 0.1, _tensT = 300;
    private int _cijConfigs = 1, _tensAxis;
    private bool _tensFixed;
    public decimal TgFromD { get => (decimal)_tgFrom; set => Set(ref _tgFrom, (double)value, nameof(TgFromD)); }
    public decimal TgToD { get => (decimal)_tgTo; set => Set(ref _tgTo, (double)value, nameof(TgToD)); }
    public decimal TgStepD { get => (decimal)_tgStep; set => Set(ref _tgStep, (double)Math.Max(1, value), nameof(TgStepD)); }
    public decimal TgPsD { get => (decimal)_tgPs; set => Set(ref _tgPs, (double)Math.Max(1, value), nameof(TgPsD)); }
    public decimal CijStrainD { get => (decimal)_cijStrain; set => Set(ref _cijStrain, (double)Math.Clamp(value, 1e-6m, 0.01m), nameof(CijStrainD)); }
    public decimal CijConfigsD { get => _cijConfigs; set => Set(ref _cijConfigs, (int)Math.Max(1, value), nameof(CijConfigsD)); }
    public decimal FluctTD { get => (decimal)_fluctT; set => Set(ref _fluctT, (double)Math.Max(1, value), nameof(FluctTD)); }
    public decimal TensRateD { get => (decimal)_tensRate; set => Set(ref _tensRate, (double)Math.Max(1e-6m, value), nameof(TensRateD)); }
    public decimal TensMaxD { get => (decimal)_tensMax; set => Set(ref _tensMax, (double)Math.Clamp(value, 0.005m, 2m), nameof(TensMaxD)); }
    public decimal TensTD { get => (decimal)_tensT; set => Set(ref _tensT, (double)Math.Max(1, value), nameof(TensTD)); }
    public int TensAxis { get => _tensAxis; set => Set(ref _tensAxis, value); }
    public bool TensFixed { get => _tensFixed; set => Set(ref _tensFixed, value); }
    public static readonly string[] Axes = ["x", "y", "z"];
    public string TensRateText => string.Format(Inv, "{0:0.##e0} s⁻¹", _tensRate * 1e12);
    protected override void OnChanged(string? name) { if (name == nameof(TensRateD)) Raise(nameof(TensRateText)); }

    // pull test (interfaces): distance and rate travel in the tensile fields when the tensile run is off
    private double _pullDist = 10, _pullRate = 2, _pullT = 300, _pullEq = 5;
    public decimal PullDistD { get => (decimal)_pullDist; set => Set(ref _pullDist, (double)Math.Clamp(value, 0.5m, 200m), nameof(PullDistD)); }
    public decimal PullRateD { get => (decimal)_pullRate; set => Set(ref _pullRate, (double)Math.Clamp(value, 0.01m, 100m), nameof(PullRateD)); }
    public decimal PullTD { get => (decimal)_pullT; set => Set(ref _pullT, (double)Math.Max(1, value), nameof(PullTD)); }
    public decimal PullEqD { get => (decimal)_pullEq; set => Set(ref _pullEq, (double)Math.Max(0, value), nameof(PullEqD)); }
    public bool PullOn => PullShearChip.IsOn || PullNormalChip.IsOn;
    private int _pullAxis;
    /// <summary>Shear axis: x, y, or z for pull-out along a fibre.</summary>
    public int PullAxis { get => _pullAxis; set => Set(ref _pullAxis, value); }
    public static readonly string[] PullAxes = ["x (slab)", "y (slab)", "z (fibre pull-out)"];

    /// <summary>Seed of the runs (cooling scans, pulls): replicas differ only in it.</summary>
    public ulong MechSeed { get; set; } = 1;

    public CapsMechOpts MechOptions()
    {
        var pull = PullOn && !TensileChip.IsOn;
        return new CapsMechOpts
        {
            Configurations = _cijConfigs, Strain = _cijStrain,
            Temperature = pull ? _pullT : FluctChip.IsOn && !TensileChip.IsOn ? _fluctT : TensileChip.IsOn ? _tensT : _fluctT,
            Axis = pull ? _pullAxis : _tensAxis, Rate = pull ? _pullRate : _tensRate, MaxStrain = pull ? _pullDist : _tensMax, LateralFixed = _tensFixed ? 1 : 0,
            TStart = _tgFrom, TEnd = _tgTo, TStep = _tgStep, PsPerStep = _tgPs, RunPs = _fluctPs,
            EquilibratePs = pull ? (_pullEq > 0 ? _pullEq : -1) : _eqPs > 0 ? _eqPs : -1,
            Seed = MechSeed,
        };
    }

    public static readonly string[] Pairs = ["all – all", "C – C", "C – H", "H – H", "C – O", "C – N", "O – H"];
    private static readonly (int A, int B)[] PairElements = [(0, 0), (6, 6), (6, 1), (1, 1), (6, 8), (6, 7), (8, 1)];

    private string _sourceText = "No trajectory open";
    public string SourceText { get => _sourceText; private set => Set(ref _sourceText, value); }
    private string _framesText = "";
    public string FramesText { get => _framesText; private set => Set(ref _framesText, value); }

    /// <summary>Called when the document or frame count changes.</summary>
    public void OnDocument(string title, CapsSummary? s)
    {
        if (s is not CapsSummary x) { SourceText = "No trajectory open"; FramesText = ""; return; }
        SourceText = title;
        FramesText = x.Frames > 1 ? string.Format(Inv, "{0:N0} frames · {1:N0} atoms", x.Frames, x.Atoms) : string.Format(Inv, "one structure · {0:N0} atoms · dynamics need frames", x.Atoms);
    }

    // ---------------------------------------------------------------- references
    public ObservableCollection<RefMaterial> References { get; } = new();
    private int _refIndex;
    public int RefIndex { get => _refIndex; set { if (Set(ref _refIndex, value)) ApplyReferences(); } }
    private RefMaterial? SelectedRef => _refIndex > 0 && _refIndex < References.Count ? References[_refIndex] : null;

    private void LoadReferences()
    {
        References.Add(new RefMaterial("", "No comparison", new()));
        var f = Paths.References;
        if (f != null)
        {
            try
            {
                using var js = JsonDocument.Parse(File.ReadAllText(f));
                foreach (var m in js.RootElement.GetProperty("materials").EnumerateArray())
                {
                    var vals = new Dictionary<string, RefValue>();
                    foreach (var v in m.GetProperty("values").EnumerateObject())
                        vals[v.Name] = new RefValue(v.Value.GetProperty("lo").GetDouble(), v.Value.GetProperty("hi").GetDouble(), Str(v.Value, "unit"),
                            Str(v.Value, "condition"), Str(v.Value, "source"));
                    References.Add(new RefMaterial(Str(m, "id"), Str(m, "name"), vals));
                }
            }
            catch (Exception e) { _status("Cannot read " + f + ": " + e.Message); }
        }
    }

    private static string Str(JsonElement e, string k) => e.TryGetProperty(k, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString()! : "";

    private void ApplyReferences()
    {
        var m = SelectedRef;
        var cards = Results.ToList();
        Results.Clear();
        foreach (var c in cards)
        {
            c.Ref = m != null && m.Values.TryGetValue(c.Id, out var r) ? r : null;
            Results.Add(c);
        }
    }

    // ---------------------------------------------------------------- run
    private bool _working;
    public bool Working { get => _working; private set { if (Set(ref _working, value)) { Raise(nameof(NotWorking)); _running(value); } } }
    public bool NotWorking => !_working;
    private string _log = "Choose calculations and run them on the open trajectory.";
    public string Log { get => _log; private set => Set(ref _log, value); }
    private double _progress;
    public double Progress { get => _progress; private set => Set(ref _progress, value); }
    private CancellationTokenSource? _cancel;

    public ObservableCollection<ResultCard> Results { get; } = new();
    public ObservableCollection<SeriesItem> Curves { get; } = new();
    private int _curveIndex = -1;
    public int CurveIndex { get => _curveIndex; set { if (Set(ref _curveIndex, value)) Raise(nameof(Curve)); } }
    public SeriesItem? Curve => _curveIndex >= 0 && _curveIndex < Curves.Count ? Curves[_curveIndex] : null;
    public bool HasResults => Results.Count > 0;
    private string _runInfo = "";
    public string RunInfo { get => _runInfo; private set => Set(ref _runInfo, value); }
    private string _json = "";

    public string[] SelectedIds => Groups.SelectMany(g => g.Chips).Where(c => c.IsOn && c.Available).Select(c => c.Id).ToArray();

    public CapsAnalyzeOpts Options()
    {
        var (a, b) = PairElements[Math.Clamp(_pair, 0, PairElements.Length - 1)];
        return new CapsAnalyzeOpts
        {
            First = _first, Last = _last < 0 ? -1 : _last, Stride = _stride, FramePs = _framePs, TimestepFs = _timestepFs, Blocks = 5,
            ElemA = a, ElemB = b, InterOnly = _inter ? 1 : 0, FitFrom = _fitFrom, FitTo = _fitTo, Probe = _probe, Grid = _grid,
            Qmax = _qmax, Dq = _dq, QDirect = _qDirect, Deuterate = _deuterate,
        };
    }

    /// <summary>The result cards again, so they show the display units chosen in Settings.</summary>
    public void RefreshDisplayUnits()
    {
        var cards = Results.ToList();
        Results.Clear();
        foreach (var c in cards) Results.Add(c);
    }

    public async Task Run()
    {
        var doc = _doc();
        if (doc == null || Working) return;
        var ids = SelectedIds;
        if (ids.Length == 0) { Log = "Switch on at least one calculation."; return; }
        Working = true;
        Progress = 0;
        _cancel = new CancellationTokenSource();
        var token = _cancel.Token;
        var o = Options();
        var sw = System.Diagnostics.Stopwatch.StartNew();
        Log = "Starting…";
        _status("Analyzing " + string.Join(", ", ids) + "…");
        try
        {
            if (FluctChip.IsOn && TensileChip.IsOn && Math.Abs(_fluctT - _tensT) > 1e-9)
                Log = "Note: the fluctuation constants and the tensile run share one temperature; the tensile temperature is used";
            var mo = MechOptions();
            var json = await Task.Run(() => doc.Analyze(string.Join(",", ids), o, mo, (what, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    if (!_working) return;
                    Progress = f;
                    Log = string.Format(Inv, "{0} · {1:F0} % · {2:F0} s", what, 100 * f, sw.Elapsed.TotalSeconds);
                });
                return !token.IsCancellationRequested;
            }));
            Load(json);
            Log = string.Format(Inv, "{0} properties in {1:F1} s", Results.Count, sw.Elapsed.TotalSeconds);
            _status("Analysis finished · " + RunInfo);
        }
        catch (Exception e)
        {
            var cancelled = e.Message.Contains("cancelled");
            Log = cancelled ? "Cancelled." : "Could not analyse: " + e.Message;
            _status(cancelled ? "Analysis cancelled" : "Analysis failed — see the Analyze page");
        }
        finally { Working = false; Progress = 0; }
    }

    public void Cancel() => _cancel?.Cancel();

    /// <summary>Reads the JSON report of caps_analyze into cards and curves.</summary>
    public void Load(string json)
    {
        _json = json;
        Results.Clear();
        Curves.Clear();
        if (json.Length == 0) { Raise(nameof(HasResults)); return; }
        using var js = JsonDocument.Parse(json);
        var root = js.RootElement;
        RunInfo = string.Format(Inv, "{0} of {1} frames · {2:N0} atoms", root.GetProperty("frames").GetInt32(), root.GetProperty("of").GetInt32(), root.GetProperty("atoms").GetInt32());
        var m = SelectedRef;
        foreach (var p in root.GetProperty("properties").EnumerateArray())
        {
            var id = Str(p, "id");
            var extra = new List<(string, double)>();
            if (p.TryGetProperty("extra", out var ex))
                foreach (var e in ex.EnumerateObject())
                    extra.Add((e.Name, e.Value.ValueKind == JsonValueKind.Number ? e.Value.GetDouble() : double.NaN));
            var notes = p.TryGetProperty("notes", out var ns) ? ns.EnumerateArray().Select(n => n.GetString() ?? "").ToArray() : [];
            var card = new ResultCard
            {
                Id = id, Name = Str(p, "name"), Unit = Str(p, "unit"), Method = Str(p, "method"), Notes = notes, Extra = extra.ToArray(),
                Value = Dbl(p, "value"), Error = Dbl(p, "error"),
                Ref = m != null && m.Values.TryGetValue(id, out var r) ? r : null,
            };
            Results.Add(card);
            if (p.TryGetProperty("series", out var ss))
            {
                var list = new List<SeriesItem>();
                foreach (var s in ss.EnumerateArray())
                {
                    var x = s.GetProperty("x").EnumerateArray().Select(v => v.ValueKind == JsonValueKind.Number ? v.GetDouble() : double.NaN).ToArray();
                    var y = s.GetProperty("y").EnumerateArray().Select(v => v.ValueKind == JsonValueKind.Number ? v.GetDouble() : double.NaN).ToArray();
                    double? refY = id is "rdf" or "sq" or "xray" or "neutron" ? 1.0 : id is "relaxation" ? 0.0 : null;
                    list.Add(new SeriesItem(card.Name, Str(s, "label"), Str(s, "x_label"), Str(s, "y_label"), x, y, id == "msd", refY));
                }
                if (id == "tg" && list.Count >= 2)          // specific volume as points, the two-line fit through them
                    list = [list[0] with { Markers = true, OverlayX = list[1].X, OverlayY = list[1].Y, Label = "specific volume and two-line fit" }, .. list.Skip(2)];
                else if (id == "tensile_modulus" && list.Count >= 2)   // raw stress as points, the smoothed curve through them
                    list = [list[1] with { Markers = true, OverlayX = list[0].X, OverlayY = list[0].Y, Label = "stress–strain (points) and smoothed" }, .. list.Skip(2)];
                foreach (var it in list) Curves.Add(it);
            }
        }
        Raise(nameof(HasResults));
        CurveIndex = -1;
        CurveIndex = Curves.Count == 0 ? -1 : Math.Max(0, Curves.ToList().FindIndex(c => c.X.Length >= 2));
    }

    private static double Dbl(JsonElement e, string k) => e.TryGetProperty(k, out var v) && v.ValueKind == JsonValueKind.Number ? v.GetDouble() : double.NaN;

    // ---------------------------------------------------------------- export
    /// <summary>Writes results.csv (one row per property, with the comparison), one CSV per curve, results.tex (a booktabs
    /// table) and results.json into a folder. Returns the number of files.</summary>
    public int Export(string dir)
    {
        Directory.CreateDirectory(dir);
        var n = 0;
        var sb = new StringBuilder("id,property,value,error,unit,reference_lo,reference_hi,reference_source,method\n");
        foreach (var c in Results)
            sb.Append(Csv(c.Id)).Append(',').Append(Csv(c.Name)).Append(',').Append(c.HasValue ? c.Value.ToString("R", Inv) : "").Append(',')
              .Append(double.IsNaN(c.Error) ? "" : c.Error.ToString("R", Inv)).Append(',').Append(Csv(c.Unit)).Append(',')
              .Append(c.Ref?.Lo.ToString("R", Inv) ?? "").Append(',').Append(c.Ref?.Hi.ToString("R", Inv) ?? "").Append(',')
              .Append(Csv(c.Ref?.Source ?? "")).Append(',').Append(Csv(c.Method)).Append('\n');
        File.WriteAllText(Path.Combine(dir, "results.csv"), sb.ToString());
        n++;
        var used = new HashSet<string>();
        foreach (var s in Curves)
        {
            var name = Slug(s.Property + "_" + s.Label);
            while (!used.Add(name)) name += "_";
            var w = new StringBuilder();
            w.Append(Csv(s.XLabel)).Append(',').Append(Csv(s.YLabel)).Append('\n');
            for (var k = 0; k < s.X.Length; k++) w.Append(s.X[k].ToString("R", Inv)).Append(',').Append(s.Y[k].ToString("R", Inv)).Append('\n');
            File.WriteAllText(Path.Combine(dir, name + ".csv"), w.ToString());
            n++;
        }
        File.WriteAllText(Path.Combine(dir, "results.tex"), Latex());
        n++;
        if (_json.Length > 0) { File.WriteAllText(Path.Combine(dir, "results.json"), _json); n++; }
        return n;
    }

    public string Latex()
    {
        var sb = new StringBuilder();
        sb.Append("% CAPS Analyze · ").Append(SourceText).Append(" · ").Append(RunInfo).Append("\n");
        sb.Append("\\begin{tabular}{llll}\n\\toprule\nProperty & Simulation & Experiment & Unit \\\\\n\\midrule\n");
        foreach (var c in Results.Where(c => c.HasValue))
        {
            var sim = c.ValueText + (c.ErrorText.Length == 0 ? "" : " $\\pm$ " + c.ErrorText[2..]);
            var exp = c.Ref == null ? "--" : (c.Ref.Lo == c.Ref.Hi ? ResultCard.Num(c.Ref.Lo) : ResultCard.Num(c.Ref.Lo) + "--" + ResultCard.Num(c.Ref.Hi));
            sb.Append(Tex(c.Name)).Append(" & ").Append(sim).Append(" & ").Append(exp).Append(" & ").Append(Tex(c.Unit)).Append(" \\\\\n");
        }
        sb.Append("\\bottomrule\n\\end{tabular}\n");
        var sources = Results.Where(c => c.Ref != null).Select(c => c.Ref!.Source).Distinct().ToArray();
        if (sources.Length > 0) sb.Append("% experiment: ").Append(string.Join("; ", sources)).Append('\n');
        return sb.ToString();
    }

    private static string Csv(string s) => s.IndexOfAny([',', '"', '\n']) >= 0 ? "\"" + s.Replace("\"", "\"\"") + "\"" : s;
    private static string Slug(string s) => new(s.Select(ch => char.IsLetterOrDigit(ch) ? char.ToLowerInvariant(ch) : '_').ToArray());
    private static string Tex(string s) => s.Replace("\\", "\\textbackslash{}").Replace("&", "\\&").Replace("%", "\\%").Replace("_", "\\_")
        .Replace("#", "\\#").Replace("^", "\\^{}").Replace("³", "$^3$").Replace("²", "$^2$").Replace("⁻¹", "$^{-1}$").Replace("⁻⁵", "$^{-5}$")
        .Replace("½", "$^{1/2}$").Replace("√", "$\\sqrt{}$").Replace("⟨", "$\\langle$").Replace("⟩", "$\\rangle$").Replace("Å", "\\AA{}")
        .Replace("δ", "$\\delta$").Replace("∞", "$_\\infty$").Replace("τ", "$\\tau$");
}
