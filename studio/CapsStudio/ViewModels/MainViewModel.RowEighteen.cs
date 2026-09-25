using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;
using System.Threading.Tasks;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record ModelRow(string Label, string A, string B, bool Strong = false);
public sealed record MonomerChoice(string Name, string Short, string Smiles);
public sealed record SolventRow(string Name, string V, string Delta, string Chi, string Predicted, string Known, bool Mismatch, bool Agrees, double Bar, bool Good, bool Bad)
{
    public bool Unknown => !Mismatch && !Agrees;
}
public sealed record MeshRow(string Axis, string L, string Grid, string Spacing);
public sealed record EwaldRow(string Tolerance, string Beta, string BetaRc, bool Current);
public sealed record ResultRow(string Quantity, string Value);

/// <summary>Row 18 of the design (Polydispersity, Copolymer, SolventScreen, Tacticity, BlendPhase, Electrostatics): the
/// statistics behind what Grow builds, and the small analytic models beside them, each saying what it assumes.</summary>
public partial class MainViewModel
{
    private static string Thin(double v, string f = "#,0") => v.ToString(f, Inv).Replace(",", " ");
    private string? DefaultCleanForceField() =>
        CleanChoices.FirstOrDefault(c => c.File != null && Path.GetFileNameWithoutExtension(c.File) == _settings.ForceField)?.File
        ?? CleanChoices.FirstOrDefault(c => c.File != null)?.File;

    /// <summary>Grow's polymer as a spec (the built-in polystyrene as a styrene repeat unit when none was chosen).</summary>
    private JsonObject GrowSpecObject()
    {
        if (_growSpec != null) return JsonNode.Parse(_growSpec)!.AsObject();
        var o = new JsonObject
        {
            ["units"] = new JsonArray(new JsonObject { ["name"] = "styrene", ["smiles"] = "[*]CC([*])C1=CC=CC=C1" }),
            ["sequence"] = "homopolymer",
            ["dp"] = _growDp,
        };
        if (DefaultCleanForceField() is { } ff) o["forcefield"] = ff;
        return o;
    }
    private string GrowShortName => _growSpec == null ? "PS" : _growSpecName;

    // ---------------------------------------------------------------- Grow › Polydispersity (design/boards/Polydispersity)
    public bool IsPolydispersity => _module == 59;
    public static readonly string[] PdDistIds = ["schulz-zimm", "flory", "poisson", "monodisperse"];
    public string[] PdDistNames { get; } = ["Schulz–Zimm", "Most probable (Flory)", "Poisson", "Monodisperse"];
    public string[] PdMatchNames { get; } = ["draw, then report the sample", "best of 100 draws (closest Đ)"];
    private int _pdDist, _pdMatch, _pdSeed = 2026;
    private decimal _pdNn = 40, _pdPdi = 1.10m;
    public int PdDist { get => _pdDist; set { if (Set(ref _pdDist, Math.Clamp(value, 0, 3))) { Raise(nameof(PdShowPdi)); PdRecompute(); } } }
    public int PdMatch { get => _pdMatch; set { if (Set(ref _pdMatch, Math.Clamp(value, 0, 1))) PdRecompute(); } }
    public decimal PdNn { get => _pdNn; set { if (Set(ref _pdNn, Math.Clamp(value, 2, 100000))) PdRecompute(); } }
    public decimal PdPdi { get => _pdPdi; set { if (Set(ref _pdPdi, Math.Clamp(value, 1.001m, 10))) PdRecompute(); } }
    public bool PdShowPdi => _pdDist == 0;
    public ObservableCollection<ModelRow> PdRows { get; } = new();
    public string[] PdLengths { get; private set; } = [];
    private int[] _pdLengths = [];
    private string _pdTitle = "", _pdK = "", _pdStatus = "", _pdM0 = "", _pdNote = "", _pdSeedText = "", _pdCellTitle = "";
    public string PdTitle { get => _pdTitle; private set => Set(ref _pdTitle, value); }
    public string PdK { get => _pdK; private set => Set(ref _pdK, value); }
    public string PdStatus { get => _pdStatus; private set => Set(ref _pdStatus, value); }
    public string PdM0 { get => _pdM0; private set => Set(ref _pdM0, value); }
    public string PdNote { get => _pdNote; private set => Set(ref _pdNote, value); }
    public string PdSeedText { get => _pdSeedText; private set => Set(ref _pdSeedText, value); }
    public string PdCellTitle { get => _pdCellTitle; private set => Set(ref _pdCellTitle, value); }
    public (double X, double Y)[] PdNumber { get; private set; } = [];
    public (double X, double Y)[] PdWeight { get; private set; } = [];
    public event Action? PdChanged;
    private double _pdSampleNn, _pdSamplePdi;
    /// <summary>The document grown with these lengths (its chains are coloured short → long on the page).</summary>
    public CapsDocument? PdGrown { get; private set; }
    public int[] PdGrownLengths { get; private set; } = [];

    public void OpenPolydispersity() { SetModule(59); PdRecompute(); }
    public void PdRedraw() { _pdSeed++; PdRecompute(); }

    private double RepeatUnitMass()
    {
        var (_, mass) = GrowChainSize();
        return _growDp > 0 ? (mass - 2 * 1.008) / _growDp : 104.15;
    }

    public void PdRecompute()
    {
        try
        {
            var m0 = RepeatUnitMass();
            var j = new JsonObject
            {
                ["distribution"] = PdDistIds[_pdDist], ["nn"] = (double)_pdNn, ["pdi"] = (double)_pdPdi, ["count"] = _growChains,
                ["seed"] = _pdSeed, ["m0"] = m0, ["best_of"] = _pdMatch == 1 ? 100 : 1,
            };
            var r = JsonNode.Parse(CapsDocument.ChainLengths(j.ToJsonString()))!;
            if ((bool?)r["ok"] != true) { PdNote = (string?)r["error"] ?? "cannot draw"; return; }
            _pdLengths = (r["lengths"] as JsonArray ?? []).Select(x => (int)(double)x!).ToArray();
            PdLengths = _pdLengths.OrderBy(x => x).Select(x => x.ToString(Inv)).ToArray();
            Raise(nameof(PdLengths));
            var s = r["sample"]!;
            var t = r["target"]!;
            double D(JsonNode n, string k) => (double?)n[k] ?? double.NaN;
            _pdSampleNn = D(s, "nn");
            _pdSamplePdi = D(s, "pdi");
            PdRows.Clear();
            PdRows.Add(new("Nₙ", D(t, "nn").ToString("0.#", Inv), _pdSampleNn.ToString("0.0", Inv)));
            PdRows.Add(new("Mₙ (g/mol)", Thin(D(t, "mn")), Thin(D(s, "mn"))));
            PdRows.Add(new("Mw (g/mol)", Thin(D(t, "mw")), Thin(D(s, "mw"))));
            PdRows.Add(new("Đ = Mw/Mₙ", D(t, "pdi").ToString("0.00", Inv), _pdSamplePdi.ToString("0.000", Inv), true));
            PdRows.Add(new("Shortest · longest", "—", $"{D(s, "min"):0} · {D(s, "max"):0}"));
            PdTitle = PdDistNames[_pdDist] + " distribution";
            var k = D(r, "k");
            PdK = _pdDist switch
            {
                0 => $"k = 1/(Đ − 1) = {k.ToString(k >= 100 ? "0" : "0.##", Inv)}",
                1 => $"Đ → 2 − 1/Nₙ = {D(t, "pdi"):0.000}",
                2 => $"Đ ≈ 1 + 1/Nₙ = {D(t, "pdi"):0.000}",
                _ => "every chain Nₙ units",
            };
            PdSeedText = $"seed {_pdSeed}";
            PdStatus = $"{GrowShortName} · {_growChains} chains · Σ N = {D(s, "sum"):0}";
            PdM0 = $"M₀ = {m0.ToString("0.00", Inv)} g/mol";
            PdNote = _pdMatch == 1
                ? $"Kept draw {(int)D(r, "kept_draw") + 1} of 100, the one with Đ closest to the target; the provenance records the sample values and that a best-of draw was used."
                : $"With only {_growChains} chains the sample Đ differs from the target; CAPS reports the sample values in the provenance, never the target as if achieved.";
            var c = r["curve"]!;
            var n = (c["n"] as JsonArray ?? []).Select(x => (double)x!).ToArray();
            PdNumber = n.Zip((c["number"] as JsonArray ?? []).Select(x => (double)x!)).ToArray();
            PdWeight = n.Zip((c["weight"] as JsonArray ?? []).Select(x => (double)x!)).ToArray();
            PdCellTitle = PdGrown != null && PdGrown == _doc ? $"Cell · {PdGrownLengths.Length} chains, short → long" : "Cell · grow with these lengths to see it";
            PdChanged?.Invoke();
        }
        catch (Exception e) { PdNote = e.Message; }
    }

    // Grow uses per-chain lengths once they are sent here
    private int[]? _growChainDp;
    private string _growChainDpText = "";
    public bool GrowPolydisperse => _growChainDp != null;
    public string GrowDispersityChip => _growChainDp == null ? "Monodisperse" : "Polydisperse · " + PdDistNames[_pdDist];
    public string GrowDispersityText => _growChainDp == null ? "Đ 1.00" : $"Đ {_pdSamplePdi.ToString("0.000", Inv)}";
    public string GrowLengthsText => _growChainDpText;

    public void UsePdLengths()
    {
        if (_pdLengths.Length == 0) return;
        _growChainDp = _pdLengths.ToArray();
        _growChainDpText = string.Format(Inv, "{0} lengths drawn from {1} (Nₙ {2:0.#}{3}, seed {4}{5}) · sample Nₙ {6:0.0}, Đ {7:0.000}, {8}–{9} units",
            _growChainDp.Length, PdDistNames[_pdDist], _pdNn, _pdDist == 0 ? $", Đ {_pdPdi:0.00}" : "", _pdSeed, _pdMatch == 1 ? ", best of 100 draws" : "",
            _pdSampleNn, _pdSamplePdi, _growChainDp.Min(), _growChainDp.Max());
        RaiseGrowDispersity();
        SetModule(0);
        Status = "Grow uses these chain lengths · " + _growChainDpText;
    }
    public void ClearPdLengths()
    {
        _growChainDp = null;
        _growChainDpText = "";
        RaiseGrowDispersity();
    }
    private void RaiseGrowDispersity()
    {
        Raise(nameof(GrowPolydisperse)); Raise(nameof(GrowDispersityChip)); Raise(nameof(GrowDispersityText)); Raise(nameof(GrowLengthsText));
    }
    /// <summary>Called by Grow: the spec with the per-chain lengths (and P_m) when set; null keeps Grow's own path.</summary>
    private string? GrowSpecWithStatistics(string? spec)
    {
        var poly = _growChainDp != null && _growChainDp.Length == _growChains;
        var stereo = _growTact == 0 && _growStereo != null;
        if (!poly && !stereo) return spec;
        var sj = spec != null ? JsonNode.Parse(spec)!.AsObject() : GrowSpecObject();
        sj["dp"] = _growDp;
        if (poly)
        {
            sj["chain_dp"] = new JsonArray(_growChainDp!.Select(x => (JsonNode)x).ToArray());
            sj["chain_lengths"] = _growChainDpText;
        }
        if (stereo)
        {
            var (kind, pm, a, b) = _growStereo!.Value;
            if (kind == 1) { sj["p_mr"] = a; sj["p_rm"] = b; sj["stereo"] = $"first-order Markov, P(r after m) {a:0.###}, P(m after r) {b:0.###}"; }
            else { sj["pm"] = pm; sj["stereo"] = $"Bernoulli, P_m {pm:0.###}"; }
        }
        return sj.ToJsonString();
    }
    private void AfterGrowStatistics(CapsDocument doc)
    {
        if (_growChainDp != null && _growChainDp.Length == _growChains) { PdGrown = doc; PdGrownLengths = _growChainDp.ToArray(); }
    }

    // ---------------------------------------------------------------- Studio › Copolymer builder (design/boards/Copolymer)
    public bool IsCopolymer => _module == 60;
    public MonomerChoice[] CoMonomers { get; } =
    [
        new("styrene", "S", "[*]CC([*])C1=CC=CC=C1"),
        new("methyl methacrylate", "MMA", "[*]CC(C)(C(=O)OC)[*]"),
        new("1,3-butadiene", "Bd", "[*]CC=CC[*]"),
        new("acrylonitrile", "AN", "[*]CC(C#N)[*]"),
        new("butyl acrylate", "BA", "[*]CC(C(=O)OCCCC)[*]"),
        new("vinyl acetate", "VAc", "[*]CC(OC(C)=O)[*]"),
    ];
    public string[] CoMonomerNames => CoMonomers.Select(m => m.Name).ToArray();
    // commonly quoted free-radical reactivity ratios (Odian, Principles of Polymerization, table 6-2): inputs, not data CAPS vouches for
    private static readonly (string A, string B, double R1, double R2, string Note)[] CoPresets =
    [
        ("S", "MMA", 0.52, 0.46, "S/MMA free-radical copolymerisation"),
        ("Bd", "S", 1.35, 0.58, "butadiene/styrene (SBR) at 50 °C"),
        ("Bd", "AN", 0.30, 0.02, "butadiene/acrylonitrile (NBR) at 40 °C"),
        ("S", "AN", 0.40, 0.04, "styrene/acrylonitrile (SAN) at 60 °C"),
    ];
    private int _coM1, _coM2 = 1, _coDp = 80, _coSeed = 12;
    private decimal _coR1 = 0.52m, _coR2 = 0.46m, _coF1 = 0.50m;
    public int CoM1 { get => _coM1; set { if (Set(ref _coM1, Math.Clamp(value, 0, CoMonomers.Length - 1))) { CoPreset(); CoRecompute(true); } } }
    public int CoM2 { get => _coM2; set { if (Set(ref _coM2, Math.Clamp(value, 0, CoMonomers.Length - 1))) { CoPreset(); CoRecompute(true); } } }
    public decimal CoR1 { get => _coR1; set { if (Set(ref _coR1, Math.Clamp(value, 0, 100))) CoRecompute(true); } }
    public decimal CoR2 { get => _coR2; set { if (Set(ref _coR2, Math.Clamp(value, 0, 100))) CoRecompute(true); } }
    public decimal CoF1 { get => _coF1; set { if (Set(ref _coF1, Math.Clamp(value, 0, 1))) CoRecompute(true); } }
    public decimal CoDp { get => _coDp; set { if (Set(ref _coDp, (int)Math.Clamp(value, 4, 2000))) CoRecompute(true); } }
    private string _coPresetNote = "", _coTitle = "", _coNote = "", _coStatus = "", _coFooter = "", _coError = "";
    public string CoPresetNote { get => _coPresetNote; private set => Set(ref _coPresetNote, value); }
    public string CoTitle { get => _coTitle; private set => Set(ref _coTitle, value); }
    public string CoNote { get => _coNote; private set => Set(ref _coNote, value); }
    public string CoStatus { get => _coStatus; private set => Set(ref _coStatus, value); }
    public string CoFooter { get => _coFooter; private set => Set(ref _coFooter, value); }
    public string CoError { get => _coError; private set => Set(ref _coError, value); }
    public string CoM1Name => CoMonomers[_coM1].Name;
    public string CoM2Name => CoMonomers[_coM2].Name;
    public string CoBuildText => $"Build {_growChains} chains";
    public ObservableCollection<ModelRow> CoRows { get; } = new();
    public int[] CoSequence { get; private set; } = [];
    public (double X, double Y)[] CoCurve { get; private set; } = [];
    public double CoF1Model { get; private set; }
    public double? CoAzeotrope { get; private set; }
    public event Action? CoChanged;
    private CapsDocument? _coDoc;
    public CapsDocument? CoDoc { get => _coDoc; private set => Set(ref _coDoc, value); }
    private bool _coBuilding;
    public bool CoBuilding { get => _coBuilding; private set => Set(ref _coBuilding, value); }

    public void OpenCopolymer() { SetModule(60); CoPreset(); CoRecompute(true); }
    public void CoRedraw() { _coSeed++; CoRecompute(true); }

    private void CoPreset()
    {
        var a = CoMonomers[_coM1].Short;
        var b = CoMonomers[_coM2].Short;
        foreach (var p in CoPresets)
        {
            if (p.A == a && p.B == b) { _coR1 = (decimal)p.R1; _coR2 = (decimal)p.R2; }
            else if (p.A == b && p.B == a) { _coR1 = (decimal)p.R2; _coR2 = (decimal)p.R1; }
            else continue;
            Raise(nameof(CoR1)); Raise(nameof(CoR2));
            CoPresetNote = $"Reactivity ratios are inputs; the values shown are commonly quoted for {p.Note}. Check your source.";
            return;
        }
        CoPresetNote = a == b ? "Both monomers are the same: choose two different ones." : $"No commonly quoted pair here for {a}/{CoMonomers[_coM2].Short}: enter reactivity ratios from your source.";
    }

    private JsonObject CoSpec()
    {
        var o = new JsonObject
        {
            ["units"] = new JsonArray(
                new JsonObject { ["name"] = CoMonomers[_coM1].Short, ["smiles"] = CoMonomers[_coM1].Smiles },
                new JsonObject { ["name"] = CoMonomers[_coM2].Short, ["smiles"] = CoMonomers[_coM2].Smiles }),
            ["sequence"] = "terminal",
            ["dp"] = _coDp,
            ["weights"] = new JsonArray((double)_coF1, 1 - (double)_coF1),
            ["r1"] = (double)_coR1,
            ["r2"] = (double)_coR2,
        };
        if (DefaultCleanForceField() is { } ff) o["forcefield"] = ff;
        return o;
    }
    private string CoName => $"P({CoMonomers[_coM1].Short}-co-{CoMonomers[_coM2].Short})";

    public void CoRecompute(bool rebuild)
    {
        try
        {
            CoError = "";
            var r = JsonNode.Parse(CapsDocument.Copolymer(new JsonObject { ["r1"] = (double)_coR1, ["r2"] = (double)_coR2, ["f1"] = (double)_coF1, ["dp"] = _coDp, ["seed"] = _coSeed }.ToJsonString()))!;
            if ((bool?)r["ok"] != true) { CoError = (string?)r["error"] ?? "cannot compute"; return; }
            double D(JsonNode n, string k) => (double?)n[k] ?? double.NaN;
            CoSequence = (r["sequence"] as JsonArray ?? []).Select(x => (int)(double)x!).ToArray();
            var ch = r["chain"]!;
            var (a, b) = (CoMonomers[_coM1].Short, CoMonomers[_coM2].Short);
            CoF1Model = D(r, "F1");
            CoAzeotrope = r["azeotrope"] is JsonValue v ? (double)v : null;
            CoRows.Clear();
            CoRows.Add(new($"F₁ ({CoMonomers[_coM1].Name})", CoF1Model.ToString("0.000", Inv), D(ch, "F1").ToString("0.000", Inv)));
            CoRows.Add(new($"Mean run, {a}", D(r, "run_a").ToString("0.00", Inv), D(ch, "run_a").ToString("0.00", Inv)));
            CoRows.Add(new($"Mean run, {b}", D(r, "run_b").ToString("0.00", Inv), D(ch, "run_b").ToString("0.00", Inv)));
            CoRows.Add(new("Longest run", "—", D(ch, "longest").ToString("0", Inv)));
            var az = CoAzeotrope is double z
                ? $"Azeotrope at f₁ = (1 − r₂)/(2 − r₁ − r₂) = {z:0.000}."
                : "No azeotrope: F₁ ≠ f₁ at every feed (one ratio above 1, the other below).";
            CoNote = $"Mayo–Lewis: F₁ = (r₁f₁² + f₁f₂) / (r₁f₁² + 2f₁f₂ + r₂f₂²) = {CoF1Model:0.000}. {az} Low-conversion model: composition drift is not included. Mayo & Lewis, J. Am. Chem. Soc. 66, 1594 (1944).";
            CoTitle = $"{CoName} · {_coDp} units";
            CoStatus = $"{_coDp} units · {D(ch, "a"):0} {a} · {D(ch, "b"):0} {b}";
            CoFooter = $"sequence from the terminal model, seed {_coSeed}";
            var c = r["curve"]!;
            CoCurve = (c["f1"] as JsonArray ?? []).Select(x => (double)x!).Zip((c["F1"] as JsonArray ?? []).Select(x => (double)x!)).ToArray();
            Raise(nameof(CoM1Name)); Raise(nameof(CoM2Name)); Raise(nameof(CoBuildText));
            CoChanged?.Invoke();
            if (rebuild) _ = BuildCoChain();
        }
        catch (Exception e) { CoError = e.Message; }
    }

    private int _coBuildTicket;
    /// <summary>One chain of the current sequence (Grow's chain 0 at this seed), for the view.</summary>
    public async Task BuildCoChain()
    {
        if (_coM1 == _coM2) return;
        var ticket = ++_coBuildTicket;
        var spec = CoSpec().ToJsonString();
        var seed = (ulong)_coSeed;
        var seq = CoSequence;
        CoBuilding = true;
        try
        {
            var (doc, _) = await Task.Run(() => CapsDocument.GrowChains(spec, new CapsGrowOpts { Chains = 1, Dp = 0, Tacticity = 0, Seed = seed, Density = 0.02, ContactScale = 0.8, Curve = 1 }, null, "copolymer"));
            if (ticket != _coBuildTicket) { doc.Dispose(); return; }
            // atoms coloured by their unit: amber for M₁, blue for M₂ (the design's category pair)
            var res = doc.AtomResidues();
            doc.SetAtomValues(res.Select(k => k >= 1 && k <= seq.Length ? (seq[k - 1] == 0 ? 1.0 : 0.0) : 0.5).ToArray(), 1);
            var old = _coDoc;
            CoDoc = doc;
            old?.Dispose();
        }
        catch (Exception e) { if (ticket == _coBuildTicket) CoError = "Chain: " + e.Message; }
        finally { if (ticket == _coBuildTicket) CoBuilding = false; }
    }

    /// <summary>Grow builds this copolymer (the terminal-model sequence per chain) and runs.</summary>
    public async Task BuildCopolymerChains()
    {
        if (_coM1 == _coM2 || _growing) return;
        _growSpec = CoSpec().ToJsonString();
        _growSpecName = CoName;
        _growDp = _coDp;
        Raise(nameof(GrowDpD)); Raise(nameof(GrowHasSpec)); Raise(nameof(GrowComponentName)); Raise(nameof(GrowAtomsText)); Raise(nameof(GrowEstimate)); Raise(nameof(GrowCommand));
        SetModule(0);
        await Grow();
    }

    // ---------------------------------------------------------------- Analyze › Solvent screen (design/boards/SolventScreen)
    public bool IsSolventScreen => _module == 61;
    private JsonObject? _solventData;
    private string[] _ssNames = [];
    /// <summary>Loaded on first use so the polymer list is there before the page binds its selection.</summary>
    public string[] SsPolymerNames { get { LoadSolvents(); return _ssNames; } }
    private int _ssPolymer;
    private decimal _ssDelta = 18.6m, _ssT = 298.15m;
    public int SsPolymer
    {
        get => _ssPolymer;
        set
        {
            if (!Set(ref _ssPolymer, Math.Clamp(value, 0, Math.Max(0, SsPolymerNames.Length - 1)))) return;
            if (SsPolymerJson() is { } p) { _ssDelta = (decimal)((double?)p["delta"] ?? 18.6); Raise(nameof(SsDelta)); }
            SsRecompute();
        }
    }
    public decimal SsDelta { get => _ssDelta; set { if (Set(ref _ssDelta, Math.Clamp(value, 5, 60))) SsRecompute(); } }
    public decimal SsT { get => _ssT; set { if (Set(ref _ssT, Math.Clamp(value, 100, 700))) SsRecompute(); } }
    public ObservableCollection<SolventRow> SsRows { get; } = new();
    private string _ssChartTitle = "", _ssChartNote = "", _ssExample = "", _ssFailText = "", _ssFailChip = "", _ssKnownHeader = "", _ssInputs = "", _ssStatus = "";
    public string SsChartTitle { get => _ssChartTitle; private set => Set(ref _ssChartTitle, value); }
    public string SsChartNote { get => _ssChartNote; private set => Set(ref _ssChartNote, value); }
    public string SsExample { get => _ssExample; private set => Set(ref _ssExample, value); }
    public string SsFailText { get => _ssFailText; private set => Set(ref _ssFailText, value); }
    public string SsFailChip { get => _ssFailChip; private set => Set(ref _ssFailChip, value); }
    public string SsKnownHeader { get => _ssKnownHeader; private set => Set(ref _ssKnownHeader, value); }
    public string SsInputs { get => _ssInputs; private set => Set(ref _ssInputs, value); }
    public string SsStatus { get => _ssStatus; private set => Set(ref _ssStatus, value); }
    public bool SsHasMismatch => SsRows.Any(r => r.Mismatch);
    /// <summary>Where the χ = 0.5 guide sits along the bars (0–1).</summary>
    public double SsHalfAt { get; private set; }

    private JsonObject? SsPolymerJson() => (_solventData?["polymers"] as JsonArray)?.OfType<JsonObject>().ElementAtOrDefault(_ssPolymer);

    private void LoadSolvents()
    {
        if (_solventData != null) return;
        try
        {
            if (Paths.Solvents is { } p) _solventData = JsonNode.Parse(File.ReadAllText(p))!.AsObject();
            _ssNames = (_solventData?["polymers"] as JsonArray ?? []).OfType<JsonObject>().Select(x => (string?)x["name"] ?? "").Append("other polymer (enter δ)").ToArray();
            if (SsPolymerJson() is { } p0) _ssDelta = (decimal)((double?)p0["delta"] ?? 18.6);
        }
        catch (Exception e) { SsFailText = "solvents.json: " + e.Message; }
    }

    public void OpenSolventScreen()
    {
        LoadSolvents();
        SetModule(61);
        SsRecompute();
    }

    public void SsRecompute()
    {
        if (_solventData == null) return;
        var solvents = _solventData["solvents"] as JsonArray ?? [];
        var poly = SsPolymerJson();
        var polyName = (string?)poly?["name"] ?? "this polymer";
        var polyShort = (string?)poly?["short"] ?? "polymer";
        var known = poly?["known"] as JsonObject;
        var t = (double)_ssT;
        var dp = (double)_ssDelta;
        var r = JsonNode.Parse(CapsDocument.SolventChi(new JsonObject { ["delta_polymer"] = dp, ["t"] = t, ["solvents"] = solvents.DeepClone() }.ToJsonString()))!;
        var chis = (r["solvents"] as JsonArray ?? []).OfType<JsonObject>().ToArray();
        var rows = solvents.OfType<JsonObject>().Zip(chis).Select(p => (S: p.First, C: p.Second)).OrderBy(p => (double)p.C["chi"]!).ToArray();
        var max = rows.Length > 0 ? Math.Max(1.0, rows.Max(p => (double)p.C["chi"]!)) : 1;
        SsRows.Clear();
        foreach (var (s, c) in rows)
        {
            var name = (string?)s["name"] ?? "";
            var chi = (double)c["chi"]!;
            var pred = (string?)c["predicted"] ?? "";
            var k = (string?)known?[name] ?? "—";
            var knownGood = k.Contains("good", StringComparison.Ordinal);
            var knownBad = k.Contains("non-solvent", StringComparison.Ordinal);
            var mismatch = (pred == "solvent" && knownBad) || (pred == "non-solvent" && knownGood);
            var agrees = (pred == "solvent" && knownGood) || (pred == "non-solvent" && knownBad) || (pred == "borderline" && k.Contains('Θ'));
            SsRows.Add(new(name, ((double)s["v"]!).ToString("0.0", Inv), ((double)s["delta"]!).ToString("0.0", Inv), chi.ToString("0.00", Inv), pred, k,
                mismatch, agrees, Math.Clamp(chi / max, 0.02, 1), pred == "solvent", pred == "non-solvent"));
        }
        SsHalfAt = 0.5 / max;
        var rt = (double?)r["rt"] ?? 8.314462618 * t;
        SsChartTitle = $"χ against {polyName}";
        SsChartNote = $"δ({polyShort}) = {dp:0.0} MPa½ · {t:0.##} K";
        var tol = solvents.OfType<JsonObject>().FirstOrDefault(s => (string?)s["name"] == "Toluene") ?? solvents.OfType<JsonObject>().FirstOrDefault();
        if (tol != null)
        {
            var v = (double)tol["v"]!;
            var dd = (double)tol["delta"]! - dp;
            SsExample = string.Format(Inv, "{0}: {1:0.0} × {2:0.00} / {3:0} + 0.34 = {4:0.000}", ((string?)tol["name"] ?? "").ToLowerInvariant(), v, dd * dd, rt, v * dd * dd / rt + 0.34);
        }
        var bad = SsRows.Where(x => x.Mismatch).ToArray();
        SsFailChip = bad.Length == 0 ? "no mismatch" : bad.Length == 1 ? "1 mismatch" : $"{bad.Length} mismatches";
        SsKnownHeader = known == null ? "KNOWN" : $"KNOWN FOR {polyShort.ToUpperInvariant()}";
        SsFailText = known == null
            ? "No known behaviour is on file for this polymer, so the estimate cannot be checked here. Hildebrand parameters ignore polarity and hydrogen bonding: treat the ranking as a first guess."
            : bad.Length == 0
                ? $"Every prediction here agrees with the known behaviour of {polyName}, but Hildebrand parameters ignore polarity and hydrogen bonding; polar solvents can still be misjudged."
                : $"Hildebrand parameters ignore polarity and hydrogen bonding. Here the estimate calls {string.Join(", ", bad.Select(b => b.Name + (b.Predicted == "solvent" ? " a solvent" : " a non-solvent")))}, against what is known for {polyName}. CAPS shows this check beside every screen; Hansen parameters or a direct χ from MD are the next step when polarity matters.";
        SsInputs = $"Inputs: {solvents.Count} values from {(string?)_solventData["source"] ?? "the solvents file"}. δ({polyShort}) is an input ({(string?)poly?["range"] ?? "enter the value you cite"}).";
        SsStatus = $"{solvents.Count} solvents · {t:0.##} K";
        Raise(nameof(SsHasMismatch)); Raise(nameof(SsHalfAt));
    }

    // ---------------------------------------------------------------- Studio › Tacticity statistics (design/boards/Tacticity)
    public bool IsTacticityStats => _module == 62;
    public string[] TsModelNames { get; } = ["Bernoulli", "Markov (first order)"];
    private int _tsModel, _tsDp = 200, _tsSeed = 99;
    private decimal _tsPm = 0.50m, _tsPmr = 0.50m, _tsPrm = 0.50m;
    public int TsModel { get => _tsModel; set { if (Set(ref _tsModel, Math.Clamp(value, 0, 1))) { Raise(nameof(TsMarkov)); Raise(nameof(TsBernoulli)); _ = TsBuild(); } } }
    public bool TsMarkov => _tsModel == 1;
    public bool TsBernoulli => _tsModel == 0;
    public decimal TsPm { get => _tsPm; set { if (Set(ref _tsPm, Math.Clamp(value, 0, 1))) _ = TsBuild(); } }
    public decimal TsPmr { get => _tsPmr; set { if (Set(ref _tsPmr, Math.Clamp(value, 0, 1))) _ = TsBuild(); } }
    public decimal TsPrm { get => _tsPrm; set { if (Set(ref _tsPrm, Math.Clamp(value, 0, 1))) _ = TsBuild(); } }
    public decimal TsDp { get => _tsDp; set { if (Set(ref _tsDp, (int)Math.Clamp(value, 5, 2000))) _ = TsBuild(); } }
    public ObservableCollection<ModelRow> TsTriads { get; } = new();
    private string _tsLabel = "", _tsDyadTitle = "", _tsSeedText = "", _tsMeso = "", _tsRacemo = "", _tsStatus = "", _tsChainHeader = "", _tsError = "", _tsMeasured = "", _tsFit = "";
    public string TsLabel { get => _tsLabel; private set => Set(ref _tsLabel, value); }
    public string TsDyadTitle { get => _tsDyadTitle; private set => Set(ref _tsDyadTitle, value); }
    public string TsSeedText { get => _tsSeedText; private set => Set(ref _tsSeedText, value); }
    public string TsMeso { get => _tsMeso; private set => Set(ref _tsMeso, value); }
    public string TsRacemo { get => _tsRacemo; private set => Set(ref _tsRacemo, value); }
    public string TsStatus { get => _tsStatus; private set => Set(ref _tsStatus, value); }
    public string TsChainHeader { get => _tsChainHeader; private set => Set(ref _tsChainHeader, value); }
    public string TsError { get => _tsError; private set => Set(ref _tsError, value); }
    public string TsMeasured { get => _tsMeasured; set => Set(ref _tsMeasured, value ?? ""); }
    public string TsFit { get => _tsFit; private set => Set(ref _tsFit, value); }
    public string[] TsPentadNames { get; private set; } = [];
    public double[] TsModelPentads { get; private set; } = [];
    public double[] TsChainPentads { get; private set; } = [];
    public int[] TsDyads { get; private set; } = [];
    public int[] TsFirstUnits { get; private set; } = [];
    private CapsDocument? _tsDoc;
    public CapsDocument? TsDoc { get => _tsDoc; private set => Set(ref _tsDoc, value); }
    private bool _tsBuilding;
    public bool TsBuilding { get => _tsBuilding; private set => Set(ref _tsBuilding, value); }
    public event Action? TsChanged;
    /// <summary>Grow's stereo model once applied: kind 0 Bernoulli (pm), 1 Markov (a = P(r after m), b = P(m after r)).</summary>
    private (int Kind, double Pm, double A, double B)? _growStereo;
    public string GrowStereoText => _growStereo is { } g ? g.Kind == 1 ? $"Markov · P(r|m) {g.A:0.###} · P(m|r) {g.B:0.###}" : $"Bernoulli · P_m {g.Pm:0.###}" : "";
    public bool GrowHasStereo => _growStereo != null;

    public void OpenTacticityStats() { SetModule(62); _ = TsBuild(); }
    public void TsRedraw() { _tsSeed++; _ = TsBuild(); }

    private JsonObject TsModelJson() => _tsModel == 1
        ? new JsonObject { ["model"] = "markov", ["p_mr"] = (double)_tsPmr, ["p_rm"] = (double)_tsPrm }
        : new JsonObject { ["model"] = "bernoulli", ["pm"] = (double)_tsPm };

    private int _tsTicket;
    public async Task TsBuild()
    {
        var ticket = ++_tsTicket;
        var spec = GrowSpecObject();
        spec["dp"] = _tsDp;
        spec.Remove("chain_dp");
        if (_tsModel == 1) { spec["p_mr"] = (double)_tsPmr; spec["p_rm"] = (double)_tsPrm; }
        else spec["pm"] = (double)_tsPm;
        var specText = spec.ToJsonString();
        var seed = (ulong)_tsSeed;
        TsBuilding = true;
        TsError = "";
        try
        {
            var (doc, _) = await Task.Run(() => CapsDocument.GrowChains(specText, new CapsGrowOpts { Chains = 1, Dp = 0, Tacticity = 0, Seed = seed, Density = 0.02, ContactScale = 0.8, Curve = 1 }, null, "tacticity"));
            if (ticket != _tsTicket) { doc.Dispose(); return; }
            var tac = JsonNode.Parse(doc.Tacticity())!;
            var dyads = (tac["chains"] as JsonArray)?.OfType<JsonObject>().Select(c => (string?)c["dyads"] ?? "").FirstOrDefault() ?? "";
            var old = _tsDoc;
            TsDoc = doc;
            old?.Dispose();
            var q = TsModelJson();
            if (dyads.Length > 0) q["dyads"] = dyads;
            var r = JsonNode.Parse(CapsDocument.Stereo(q.ToJsonString()))!;
            if ((bool?)r["ok"] != true) { TsError = (string?)r["error"] ?? "cannot compute"; return; }
            var m = r["model"]!;
            double D(JsonNode? n, string k) => (double?)n?[k] ?? double.NaN;
            TsPentadNames = (r["names"] as JsonArray ?? []).Select(x => (string?)x ?? "").ToArray();
            TsModelPentads = (m["pentads"] as JsonArray ?? []).Select(x => (double)x!).ToArray();
            var ch = r["chain"];
            TsChainPentads = (ch?["pentads"] as JsonArray ?? []).Select(x => (double)x!).ToArray();
            TsDyads = dyads.Select(c => c == 'm' ? 0 : 1).ToArray();
            TsTriads.Clear();
            string C(string k) => ch == null ? "—" : D(ch, k).ToString("0.000", Inv);
            TsTriads.Add(new("mm", D(m, "mm").ToString("0.000", Inv), C("mm")));
            TsTriads.Add(new("mr", D(m, "mr").ToString("0.000", Inv), C("mr")));
            TsTriads.Add(new("rr", D(m, "rr").ToString("0.000", Inv), C("rr")));
            var nTriads = ch == null ? 0 : (int)D(ch, "triads");
            TsChainHeader = $"CHAIN · {nTriads}";
            TsDyadTitle = dyads.Length == 0 ? "Dyad sequence · this polymer has no stereocentres" : $"Dyad sequence · {dyads.Length} dyads";
            TsSeedText = $"seed {_tsSeed}";
            TsMeso = $"m (meso) · {(ch == null ? 0 : D(ch, "m")):0}";
            TsRacemo = $"r (racemo) · {(ch == null ? 0 : D(ch, "r")):0}";
            var pm = D(m, "pm");
            TsLabel = pm > 0.9 ? "isotactic" : pm < 0.1 ? "syndiotactic" : "atactic";
            TsStatus = _tsModel == 1
                ? $"{GrowShortName} · DP {_tsDp} · P(r|m) = {_tsPmr:0.##} · P(m|r) = {_tsPrm:0.##}"
                : $"{GrowShortName} · DP {_tsDp} · P_m = {_tsPm:0.##}";
            // the first eight units framed in the view
            var res = doc.AtomResidues();
            TsFirstUnits = Enumerable.Range(0, res.Length).Where(i => res[i] >= 1 && res[i] <= 8).ToArray();
            TsChanged?.Invoke();
        }
        catch (Exception e) { if (ticket == _tsTicket) TsError = e.Message; }
        finally { if (ticket == _tsTicket) TsBuilding = false; }
    }

    /// <summary>Fits Bernoulli and Markov models to measured pentads (ten numbers in the order shown).</summary>
    public void TsFitMeasured()
    {
        var parts = _tsMeasured.Split([',', ' ', ';', '\t'], StringSplitOptions.RemoveEmptyEntries);
        var vals = parts.Select(p => double.TryParse(p, NumberStyles.Float, Inv, out var v) ? v : double.NaN).ToArray();
        if (vals.Length != 10 || vals.Any(v => !double.IsFinite(v) || v < 0)) { TsFit = "Enter ten pentad fractions in the order mmmm, mmmr, rmmr, mmrr, mmrm, rmrm, rmrr, rrrr, rrrm, mrrm."; return; }
        var r = JsonNode.Parse(CapsDocument.Stereo(new JsonObject { ["measured"] = new JsonArray(vals.Select(v => (JsonNode)v).ToArray()) }.ToJsonString()))!;
        var f = r["fit"]!;
        double D(JsonNode? n, string k) => (double?)n?[k] ?? double.NaN;
        var b = f["bernoulli"];
        var mk = f["markov"];
        var ok = Math.Abs(D(f, "mm_rr") - D(f, "mr2_4")) < 0.1 * Math.Max(1e-6, D(f, "mr2_4"));
        TsFit = string.Format(Inv, "Bernoulli: P_m = {0:0.000} (rms {1:0.0000}). Markov: P(r|m) = {2:0.000}, P(m|r) = {3:0.000} (rms {4:0.0000}). Test mm·rr = {5:0.0000} vs (mr/2)² = {6:0.0000}: {7}",
            D(b, "pm"), D(f, "bernoulli_rms"), D(mk, "p_mr"), D(mk, "p_rm"), D(f, "markov_rms"), D(f, "mm_rr"), D(f, "mr2_4"),
            ok ? "Bernoulli is consistent." : "Bernoulli fails; use the Markov values.");
        if (ok) { _tsModel = 0; _tsPm = (decimal)Math.Round(D(b, "pm"), 3); }
        else { _tsModel = 1; _tsPmr = (decimal)Math.Round(D(mk, "p_mr"), 3); _tsPrm = (decimal)Math.Round(D(mk, "p_rm"), 3); }
        Raise(nameof(TsModel)); Raise(nameof(TsMarkov)); Raise(nameof(TsBernoulli)); Raise(nameof(TsPm)); Raise(nameof(TsPmr)); Raise(nameof(TsPrm));
        _ = TsBuild();
    }

    /// <summary>Grow draws its atactic chains from this model.</summary>
    public void TsApply()
    {
        _growStereo = _tsModel == 1 ? (1, double.NaN, (double)_tsPmr, (double)_tsPrm) : (0, (double)_tsPm, double.NaN, double.NaN);
        _growTact = 0;
        Raise(nameof(GrowTacticity)); Raise(nameof(GrowStereoText)); Raise(nameof(GrowHasStereo)); Raise(nameof(GrowTacticityName));
        SetModule(0);
        Status = "Grow draws atactic chains from " + GrowStereoText;
    }
    public void ClearGrowStereo() { _growStereo = null; Raise(nameof(GrowStereoText)); Raise(nameof(GrowHasStereo)); }

    // ---------------------------------------------------------------- Analyze › Blend phase diagram (design/boards/BlendPhase)
    public bool IsBlendPhase => _module == 63;
    private decimal _bpNa = 100, _bpNb = 200, _bpA = -0.02m, _bpB = 15, _bpT = 300;
    public decimal BpNa { get => _bpNa; set { if (Set(ref _bpNa, Math.Clamp(value, 1, 1e6m))) BpRecompute(); } }
    public decimal BpNb { get => _bpNb; set { if (Set(ref _bpNb, Math.Clamp(value, 1, 1e6m))) BpRecompute(); } }
    public decimal BpA { get => _bpA; set { if (Set(ref _bpA, Math.Clamp(value, -10, 10))) BpRecompute(); } }
    public decimal BpB { get => _bpB; set { if (Set(ref _bpB, Math.Clamp(value, -1e5m, 1e5m))) BpRecompute(); } }
    public decimal BpT { get => _bpT; set { if (Set(ref _bpT, Math.Clamp(value, 1, 5000))) BpRecompute(); } }
    public ObservableCollection<ResultRow> BpResults { get; } = new();
    private string _bpNote = "", _bpStatus = "";
    public string BpNote { get => _bpNote; private set => Set(ref _bpNote, value); }
    public string BpStatus { get => _bpStatus; private set => Set(ref _bpStatus, value); }
    public (double X, double Y)[] BpBinodal { get; private set; } = [];
    public (double X, double Y)[] BpSpinodal { get; private set; } = [];
    public (double X, double Y)? BpCritical { get; private set; }
    public double BpTNow => (double)_bpT;
    public event Action? BpChanged;

    public void OpenBlendPhase() { SetModule(63); BpRecompute(); }

    public void BpRecompute()
    {
        try
        {
            var r = JsonNode.Parse(CapsDocument.BlendPhase(new JsonObject { ["na"] = (double)_bpNa, ["nb"] = (double)_bpNb, ["a"] = (double)_bpA, ["b"] = (double)_bpB, ["t"] = (double)_bpT }.ToJsonString()))!;
            double D(JsonNode? n, string k) => (double?)n?[k] ?? double.NaN;
            (double, double)[] Curve(JsonNode? c) => c == null ? [] : (c["phi"] as JsonArray ?? []).Select(x => (double)x!).Zip((c["t"] as JsonArray ?? []).Select(x => (double)x!)).ToArray();
            BpBinodal = Curve(r["binodal"]);
            BpSpinodal = Curve(r["spinodal_curve"]);
            var tc = r["tc"] is JsonValue tv ? (double)tv : double.NaN;
            BpCritical = double.IsFinite(tc) ? (D(r, "phi_c"), tc) : null;
            string Pair(JsonNode? a) => a is JsonArray p && p.Count == 2 ? $"{(double)p[0]!:0.0000} · {(double)p[1]!:0.0000}" : "one phase";
            BpResults.Clear();
            BpResults.Add(new("χ_c = ½ (1/√N_A + 1/√N_B)²", D(r, "chi_c").ToString("0.00000", Inv)));
            BpResults.Add(new("φ_c = √N_B / (√N_A + √N_B)", D(r, "phi_c").ToString("0.0000", Inv)));
            BpResults.Add(new("T_c = B / (χ_c − A)", double.IsFinite(tc) ? tc.ToString("0.0", Inv) + " K" : "none (χ never reaches χ_c)"));
            BpResults.Add(new($"χ at {_bpT:0.#} K", D(r, "chi_t").ToString("0.0000", Inv)));
            BpResults.Add(new($"Coexisting at {_bpT:0.#} K", Pair(r["coexist"])));
            BpResults.Add(new($"Spinodal at {_bpT:0.#} K", Pair(r["spinodal"])));
            var kind = (string?)r["kind"];
            BpNote = kind == "ucst"
                ? "Two phases below the binodal. χ(T) = A + B/T gives an upper critical solution temperature because B > 0."
                : kind == "lcst" ? "Two phases above the binodal. χ(T) = A + B/T with B < 0 gives a lower critical solution temperature."
                : "B = 0: χ does not depend on temperature, so there is no critical temperature; the blend is one phase when χ < χ_c.";
            BpStatus = $"N_A = {_bpNa:0} · N_B = {_bpNb:0}";
            BpChanged?.Invoke();
        }
        catch (Exception e) { BpNote = e.Message; }
    }

    public string BlendCsv()
    {
        var sb = new System.Text.StringBuilder("# CAPS blend phase diagram (Flory–Huggins), chi = A + B/T\n");
        sb.Append(string.Format(Inv, "# N_A {0}, N_B {1}, A {2}, B {3} K\ncurve,phi_A,T_K\n", _bpNa, _bpNb, _bpA, _bpB));
        foreach (var (x, y) in BpBinodal) sb.Append(string.Format(Inv, "binodal,{0:0.######},{1:0.###}\n", x, y));
        foreach (var (x, y) in BpSpinodal) sb.Append(string.Format(Inv, "spinodal,{0:0.######},{1:0.###}\n", x, y));
        return sb.ToString();
    }

    // ---------------------------------------------------------------- Dynamics › Electrostatics (design/boards/Electrostatics)
    public bool IsElectrostatics => _module == 64;
    public string[] EsMethods { get; } = ["PME", "Damped shifted force"];
    public string[] EsTolerances { get; } = ["1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"];
    private int _esMethod, _esTol = 2;
    private decimal _esOrder = 4, _esCutoff = 12, _esSpacing = 1.2m, _esBox = 45.3m;
    public int EsMethod { get => _esMethod; set { if (Set(ref _esMethod, Math.Clamp(value, 0, 1))) { Raise(nameof(EsPme)); EsRecompute(); } } }
    public bool EsPme => _esMethod == 0;
    public int EsTolerance { get => _esTol; set { if (Set(ref _esTol, Math.Clamp(value, 0, 5))) EsRecompute(); } }
    public decimal EsOrder { get => _esOrder; set { if (Set(ref _esOrder, Math.Clamp(value, 3, 10))) EsRecompute(); } }
    public decimal EsCutoff { get => _esCutoff; set { if (Set(ref _esCutoff, Math.Clamp(value, 4, 30))) EsRecompute(); } }
    public decimal EsSpacing { get => _esSpacing; set { if (Set(ref _esSpacing, Math.Clamp(value, 0.3m, 3))) EsRecompute(); } }
    public decimal EsBox { get => _esBox; set { if (Set(ref _esBox, Math.Clamp(value, 5, 10000))) EsRecompute(); } }
    public bool EsNoCell => _doc == null || _doc.Summary().CellValid == 0;
    public ObservableCollection<EwaldRow> EsTable { get; } = new();
    public ObservableCollection<MeshRow> EsMesh { get; } = new();
    private string _esBeta = "", _esMeshText = "", _esFft = "", _esStatus = "", _esStatusRight = "", _esWarn = "";
    public string EsBeta { get => _esBeta; private set => Set(ref _esBeta, value); }
    public string EsMeshText { get => _esMeshText; private set => Set(ref _esMeshText, value); }
    public string EsFft { get => _esFft; private set => Set(ref _esFft, value); }
    public string EsStatus { get => _esStatus; private set => Set(ref _esStatus, value); }
    public string EsStatusRight { get => _esStatusRight; private set => Set(ref _esStatusRight, value); }
    public string EsWarn { get => _esWarn; private set => Set(ref _esWarn, value); }
    public (double X, double Y)[] EsCurve { get; private set; } = [];
    public double EsTolValue => Math.Pow(10, -3 - _esTol);
    public event Action? EsChanged;

    public void OpenElectrostatics()
    {
        _esMethod = _settings.Electrostatics == 1 ? 0 : 1;
        _esTol = Math.Clamp(-(int)Math.Round(Math.Log10(_settings.EwaldRtol)) - 3, 0, 5);
        _esOrder = _settings.PmeOrder;
        _esSpacing = (decimal)_settings.PmeSpacing;
        _esCutoff = (decimal)_relaxCutoff;
        foreach (var n in new[] { nameof(EsMethod), nameof(EsPme), nameof(EsTolerance), nameof(EsOrder), nameof(EsSpacing), nameof(EsCutoff), nameof(EsNoCell) }) Raise(n);
        SetModule(64);
        EsRecompute();
    }

    public void EsRecompute()
    {
        try
        {
            var j = new JsonObject { ["cutoff"] = (double)_esCutoff, ["tolerance"] = EsTolValue, ["spacing"] = (double)_esSpacing, ["order"] = (int)_esOrder };
            if (EsNoCell) j["edges"] = new JsonArray((double)_esBox, (double)_esBox, (double)_esBox);
            var r = JsonNode.Parse(CapsDocument.EwaldParams(EsNoCell ? null : _doc, j.ToJsonString()))!;
            double D(JsonNode? n, string k) => (double?)n?[k] ?? double.NaN;
            var beta = D(r, "beta");
            EsTable.Clear();
            foreach (var t in (r["table"] as JsonArray ?? []).OfType<JsonObject>())
            {
                var tol = D(t, "tolerance");
                EsTable.Add(new("10" + Sup((int)Math.Round(Math.Log10(tol))), D(t, "beta").ToString("0.0000", Inv), D(t, "beta_rc").ToString("0.000", Inv), Math.Abs(tol - EsTolValue) < 1e-15));
            }
            EsMesh.Clear();
            var edges = (r["edges"] as JsonArray ?? []).Select(x => (double)x!).ToArray();
            var mesh = (r["mesh"] as JsonArray ?? []).Select(x => (int)(double)x!).ToArray();
            var sp = (r["spacing"] as JsonArray ?? []).Select(x => (double)x!).ToArray();
            for (var k = 0; k < Math.Min(3, mesh.Length); ++k) EsMesh.Add(new("xyz"[k].ToString(), edges[k].ToString("0.0", Inv), mesh[k].ToString(Inv), sp[k].ToString("0.0000", Inv)));
            EsBeta = beta.ToString("0.0000", Inv) + " Å⁻¹";
            EsMeshText = mesh.Length == 3 ? $"{mesh[0]} × {mesh[1]} × {mesh[2]}" : "—";
            EsFft = "FFT sizes: factors 2·3·5·7";
            var cube = edges.Length == 3 && Math.Abs(edges[0] - edges[1]) < 1e-6 && Math.Abs(edges[1] - edges[2]) < 1e-6;
            EsStatus = edges.Length == 3 ? (cube ? $"box {edges[0]:0.0} Å cubic" : $"box {edges[0]:0.0} × {edges[1]:0.0} × {edges[2]:0.0} Å") + (EsNoCell ? " (entered: no periodic cell open)" : "") : "";
            EsStatusRight = _esMethod == 0 ? $"β = {beta:0.0000} Å⁻¹ · mesh {(mesh.Length == 3 && cube ? mesh[0] + "³" : EsMeshText)}" : "damped shifted force · no mesh";
            EsWarn = (bool?)r["fits"] == false ? $"The cut-off is more than half the shortest box edge: use a cut-off ≤ {edges.Min() / 2:0.0} Å or a larger box." : "";
            var c = r["curve"]!;
            EsCurve = (c["r"] as JsonArray ?? []).Select(x => (double)x!).Zip((c["erfc"] as JsonArray ?? []).Select(x => (double)x!)).ToArray();
            EsChanged?.Invoke();
        }
        catch (Exception e) { EsWarn = e.Message; }
    }
    private static string Sup(int n)
    {
        const string d = "⁰¹²³⁴⁵⁶⁷⁸⁹";
        var s = Math.Abs(n).ToString(Inv).Select(ch => d[ch - '0']);
        return (n < 0 ? "⁻" : "") + string.Concat(s);
    }

    public void EsApply()
    {
        _settings.Electrostatics = _esMethod == 0 ? 1 : 0;
        _settings.EwaldRtol = EsTolValue;
        _settings.PmeSpacing = (double)_esSpacing;
        _settings.PmeOrder = (int)_esOrder;
        _relaxCutoff = (double)_esCutoff;
        ApplyElectrostatics();
        Raise(nameof(SetElectrostatics)); Raise(nameof(PmeOn)); Raise(nameof(SetEwaldExponent)); Raise(nameof(SetPmeSpacing)); Raise(nameof(SetPmeOrder));
        Raise(nameof(RelaxCutoffD)); Raise(nameof(ElectrostaticsText));
        Changed("Electrostatics");
        Status = $"Electrostatics: {(_esMethod == 0 ? "PME" : "damped shifted force")} · cut-off {_esCutoff:0.#} Å" + (_esMethod == 0 ? $" · β {EsBeta} · mesh {EsMeshText}" : "");
    }
    public void EsReset()
    {
        // the Studio's defaults (Settings › Force fields): damped shifted force; for PME 1e-5, 1.0 Å, order 5; 10 Å cut-off
        _esMethod = 1; _esTol = 2; _esOrder = 5; _esSpacing = 1.0m; _esCutoff = 10;
        foreach (var n in new[] { nameof(EsMethod), nameof(EsPme), nameof(EsTolerance), nameof(EsOrder), nameof(EsSpacing), nameof(EsCutoff) }) Raise(n);
        EsRecompute();
    }
}
