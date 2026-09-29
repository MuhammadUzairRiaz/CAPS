using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A solvent of the library: id, name, density and what it is for.</summary>
public sealed record SolventItem(string Id, string Name, double Density, string Use)
{
    public override string ToString() => Name;
}

/// <summary>A packing stage of the solvation progress list.</summary>
public sealed class SolvStage : ObservableObject
{
    private string _title;
    public SolvStage(string title) { _title = title; }
    public string Title { get => _title; set => Set(ref _title, value); }
    private int _state;   // 0 pending, 1 running, 2 done, 3 failed
    public int State { get => _state; set { if (Set(ref _state, value)) { Raise(nameof(StateText)); Raise(nameof(IsRunning)); Raise(nameof(IsDone)); Raise(nameof(IsFailed)); } } }
    public string StateText => _state switch { 1 => "running", 2 => "done", 3 => "failed", _ => "pending" };
    public bool IsRunning => _state == 1;
    public bool IsDone => _state == 2;
    public bool IsFailed => _state == 3;
}

/// <summary>Solvation builder (design/boards/SolvationBuilder): the open structure held at the centre of a periodic
/// box, solvent and ions packed around it with CAPS Pack; the counts from the free volume, the density and the salt
/// concentration; a sparse sample of the solvent in the preview.</summary>
public sealed partial class MainViewModel
{
    public bool IsSolvation => _module == 31;
    public List<SolventItem> Solvents { get; } = new();
    public List<string> Salts { get; } = new();
    public static readonly string[] SolvShapes = ["Cubic", "Rectangular", "Padded"];   // padded: the solute's extent plus the padding
    public static readonly string[] WaterModels = ["TIP4P/2005", "TIP3P", "SPC/E"];
    public SolvStage[] SolvStages { get; } = [new("Random sequential insertion"), new("Overlap minimisation · L-BFGS"), new("Verify min. distance ≥ 2.0 Å")];

    public void OpenSolvation()
    {
        if (Solvents.Count == 0)
        {
            try
            {
                var lib = JsonNode.Parse(CapsDocument.SolventLibrary())!;
                foreach (var s in lib["solvents"]!.AsArray())
                    Solvents.Add(new SolventItem(s!["id"]!.GetValue<string>(), s["name"]!.GetValue<string>(), s["density"]!.GetValue<double>(), s["use"]!.GetValue<string>()));
                foreach (var s in lib["salts"]!.AsArray()) Salts.Add(s!["id"]!.GetValue<string>());
            }
            catch (Exception e) { SolvError = "Cannot read the solvent library: " + e.Message; }
            _solvent = Solvents.FirstOrDefault();
            _salt = Salts.FirstOrDefault() ?? "NaCl";
            Raise(nameof(Solvents)); Raise(nameof(Salts)); Raise(nameof(SolvSolvent)); Raise(nameof(SolvSalt));
        }
        _solvUseSolute = _doc != null;
        if (_doc != null && _solvShape == 0) _solvShape = 2;   // around a structure: its extent plus padding
        Raise(nameof(SolvUseSolute)); Raise(nameof(SolvShape));
        RaiseSolv();
        SetModule(31);
        SolvRefresh();
    }

    // ---------------------------------------------------------------- solute and box

    private bool _solvUseSolute;
    public bool SolvHasDocument => _doc != null;
    public bool SolvUseSolute { get => _solvUseSolute && _doc != null; set { if (Set(ref _solvUseSolute, value)) { if (!value && _solvShape == 2) SolvShape = 0; RaiseSolv(); SolvRefresh(); } } }
    public string SolvSoluteText => _doc == null ? "No structure open: a box of pure solvent" : $"{Title}";

    private int _solvShape;
    public int SolvShape { get => _solvShape; set { if (Set(ref _solvShape, Math.Clamp(value, 0, 2))) { RaiseSolv(); SolvRefresh(); } } }
    public bool SolvIsCubic => _solvShape == 0;
    public bool SolvIsRect => _solvShape == 1;
    public bool SolvIsPadding => _solvShape == 2;
    private decimal _solvEdge = 30, _solvA = 30, _solvB = 30, _solvC = 30, _solvPad = 10, _solvTol = 2.0m;
    public decimal SolvEdge { get => _solvEdge; set { if (Set(ref _solvEdge, Math.Clamp(value, 8, 300))) SolvRefresh(); } }
    public decimal SolvA { get => _solvA; set { if (Set(ref _solvA, Math.Clamp(value, 8, 300))) SolvRefresh(); } }
    public decimal SolvB { get => _solvB; set { if (Set(ref _solvB, Math.Clamp(value, 8, 300))) SolvRefresh(); } }
    public decimal SolvC { get => _solvC; set { if (Set(ref _solvC, Math.Clamp(value, 8, 300))) SolvRefresh(); } }
    public decimal SolvPadding { get => _solvPad; set { if (Set(ref _solvPad, Math.Clamp(value, 2, 100))) SolvRefresh(); } }
    public decimal SolvTolerance { get => _solvTol; set { if (Set(ref _solvTol, Math.Clamp(value, 1, 4))) { SolvStages[2].State = 0; SolvStages[2].Title = SolvVerifyTitle; Raise(nameof(SolvVerifyTitle)); SolvRefresh(); } } }
    public string SolvVerifyTitle => $"Verify min. distance ≥ {_solvTol:0.0} Å";

    // ---------------------------------------------------------------- solvent

    private SolventItem? _solvent;
    public SolventItem? SolvSolvent { get => _solvent; set { if (value != null && Set(ref _solvent, value)) { _solvDensity = 0; Raise(nameof(SolvDensity)); RaiseSolv(); SolvRefresh(); } } }
    public bool SolvIsWater => _solvent?.Id == "water";
    private string _waterModel = "TIP4P/2005";
    public string SolvWaterModel { get => _waterModel; set { if (value != null && Set(ref _waterModel, value)) SolvRefresh(); } }
    private decimal _solvDensity;
    /// <summary>g/cm³; shows the solvent's own density until edited.</summary>
    public decimal SolvDensity { get => _solvDensity > 0 ? _solvDensity : (decimal)(_solvent?.Density ?? 1.0); set { if (Set(ref _solvDensity, Math.Clamp(value, 0.1m, 3m))) SolvRefresh(); } }
    private string _solvMolecules = "auto";
    public string SolvMolecules { get => _solvMolecules; set { if (Set(ref _solvMolecules, value ?? "auto")) SolvRefresh(); } }
    public string SolvSolventUse => _solvent?.Use ?? "";

    // ---------------------------------------------------------------- ions

    private int _ionMode = 2;
    public int SolvIonMode { get => _ionMode; set { if (Set(ref _ionMode, Math.Clamp(value, 0, 3))) { RaiseSolv(); SolvRefresh(); } } }
    public bool SolvIonNone { get => _ionMode == 0; set { if (value) SolvIonMode = 0; } }
    public bool SolvIonNeutral { get => _ionMode == 1; set { if (value) SolvIonMode = 1; } }
    public bool SolvIonConc { get => _ionMode == 2; set { if (value) SolvIonMode = 2; } }
    public bool SolvIonCustom { get => _ionMode == 3; set { if (value) SolvIonMode = 3; } }
    public bool SolvShowSalt => _ionMode != 0;
    private string _salt = "NaCl";
    public string SolvSalt { get => _salt; set { if (value != null && Set(ref _salt, value)) SolvRefresh(); } }
    private decimal _conc = 0.15m, _cations, _anions;
    public decimal SolvConcentration { get => _conc; set { if (Set(ref _conc, Math.Clamp(value, 0, 5))) SolvRefresh(); } }
    public decimal SolvCations { get => _cations; set { if (Set(ref _cations, Math.Clamp(Math.Round(value), 0, 5000))) SolvRefresh(); } }
    public decimal SolvAnions { get => _anions; set { if (Set(ref _anions, Math.Clamp(Math.Round(value), 0, 5000))) SolvRefresh(); } }

    private void RaiseSolv()
    {
        foreach (var n in new[] { nameof(SolvIsCubic), nameof(SolvIsRect), nameof(SolvIsPadding), nameof(SolvIsWater), nameof(SolvIonNone), nameof(SolvIonNeutral), nameof(SolvIonConc),
                                  nameof(SolvIonCustom), nameof(SolvShowSalt), nameof(SolvSolventUse), nameof(SolvHasDocument), nameof(SolvSoluteText), nameof(SolvUseSolute), nameof(SolvVerifyTitle) })
            Raise(n);
    }

    private string SolvOptions(int? sample = null) => new JsonObject
    {
        ["shape"] = _solvShape, ["edge"] = (double)_solvEdge, ["edges"] = new JsonArray((double)_solvA, (double)_solvB, (double)_solvC), ["padding"] = (double)_solvPad,
        ["tolerance"] = (double)_solvTol, ["solvent"] = _solvent?.Id ?? "water", ["water_model"] = _waterModel,
        ["density"] = _solvDensity > 0 ? (double)_solvDensity : 0,
        ["molecules"] = sample ?? (int.TryParse(_solvMolecules, out var m) && m > 0 ? m : 0),
        ["ion_mode"] = _ionMode, ["salt"] = _salt, ["concentration"] = (double)_conc, ["cations"] = (int)_cations, ["anions"] = (int)_anions, ["seed"] = SolvSeed.Take(),
    }.ToJsonString();

    // ---------------------------------------------------------------- plan, preview, run

    private string _solvCounts = "", _solvError = "", _solvLog = "", _solvBox = "", _solvHud = "";
    public string SolvCounts { get => _solvCounts; private set => Set(ref _solvCounts, value); }
    public string SolvBoxText { get => _solvBox; private set => Set(ref _solvBox, value); }
    public string SolvHud { get => _solvHud; private set => Set(ref _solvHud, value); }
    public string SolvError { get => _solvError; private set { if (Set(ref _solvError, value)) Raise(nameof(SolvHasError)); } }
    public bool SolvHasError => _solvError.Length > 0;
    public string SolvLog { get => _solvLog; private set => Set(ref _solvLog, value); }
    private string _cationLabel = "Na⁺", _anionLabel = "Cl⁻";
    public string SolvCationLabel { get => _cationLabel; private set => Set(ref _cationLabel, value); }
    public string SolvAnionLabel { get => _anionLabel; private set => Set(ref _anionLabel, value); }
    private int _solvPlanned;
    private CapsDocument? _solvDoc;
    public CapsDocument? SolvDoc { get => _solvDoc; private set { if (Set(ref _solvDoc, value)) SolvViewChanged?.Invoke(); } }
    public event Action? SolvViewChanged;
    private int _solvTicket;

    private static string Sup(string el, int z) => el + (Math.Abs(z) == 2 ? "²" : "") + (z > 0 ? "⁺" : "⁻");

    /// <summary>The counts (at once) and a preview with a sparse sample of the solvent (in the background).</summary>
    public void SolvRefresh()
    {
        if (_module != 31 || SolvRunning) return;
        var solute = SolvUseSolute ? _doc : null;
        try
        {
            var p = JsonNode.Parse(CapsDocument.SolvatePlan(solute, SolvOptions()))!;
            if (p["ok"]?.GetValue<bool>() != true) { SolvError = p["error"]?.GetValue<string>() ?? "cannot plan"; SolvCounts = "—"; return; }
            SolvError = "";
            int nsol = (int)p["solvent"]!.GetValue<double>(), nc = (int)p["cations"]!.GetValue<double>(), na = (int)p["anions"]!.GetValue<double>();
            _solvPlanned = nsol;
            var zc = _salt is "CaCl2" or "MgCl2" or "ZnCl2" ? 2 : 1;
            SolvCationLabel = Sup(p["cation"]!.GetValue<string>(), zc);
            SolvAnionLabel = Sup(p["anion"]!.GetValue<string>(), -1);
            var name = p["solvent_name"]!.GetValue<string>();
            SolvCounts = string.Format(CultureInfo.InvariantCulture, "{0:N0} {1} · {2} {3} · {4} {5}", nsol, name.Split(' ')[0].ToLowerInvariant(), nc, SolvCationLabel, na, SolvAnionLabel);
            var box = p["box"]!.AsArray().Select(x => x!.GetValue<double>()).ToArray();
            SolvBoxText = string.Format(CultureInfo.InvariantCulture, "{0:0.0} × {1:0.0} × {2:0.0} Å · free {3:N0} of {4:N0} Å³ · {5:0.000} g/cm³ · {6:0.000} mol/L",
                box[0], box[1], box[2], p["free_volume"]!.GetValue<double>(), p["box_volume"]!.GetValue<double>(), p["density"]!.GetValue<double>(), p["concentration"]!.GetValue<double>());
            var notes = p["notes"]!.AsArray().Select(x => x!.GetValue<string>()).ToList();
            SolvLog = notes.Count > 0 ? string.Join("\n", notes) : (solute != null ? $"Solute charge {p["solute_charge"]!.GetValue<double>():+0.##;−0.##;0} e" : "");
        }
        catch (Exception e) { SolvError = e.Message; return; }
        // preview: the same box with at most 150 solvent molecules
        var ticket = ++_solvTicket;
        var sample = Math.Min(150, Math.Max(1, _solvPlanned));
        var opts = SolvOptions(sample);
        SolvHud = _solvPlanned > sample ? $"Preview · {sample} of {_solvPlanned:N0} solvent molecules shown" : "Preview";
        Task.Run(() =>
        {
            try { var (d, _) = CapsDocument.Solvate(solute, opts, null, "preview"); return (Doc: (CapsDocument?)d, Error: (string?)null); }
            catch (Exception e) { return (Doc: (CapsDocument?)null, Error: (string?)e.Message); }
        }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            var (doc, err) = t.Result;
            if (ticket != _solvTicket) { doc?.Dispose(); return; }
            if (err != null || doc == null) { SolvError = err ?? "cannot preview"; return; }
            var old = _solvDoc;
            SolvDoc = doc;
            old?.Dispose();
        }));
    }

    private bool _solvRunning;
    public bool SolvRunning { get => _solvRunning; private set { if (Set(ref _solvRunning, value)) Raise(nameof(SolvIdle)); } }
    public bool SolvIdle => !_solvRunning;
    private double _solvProgress;
    public double SolvProgress { get => _solvProgress; private set => Set(ref _solvProgress, value); }
    private string _solvStageChip = "3 stages";
    public string SolvStageChip { get => _solvStageChip; private set => Set(ref _solvStageChip, value); }

    /// <summary>Packs the full box and opens it as the Studio document.</summary>
    public async Task Solvate()
    {
        if (SolvRunning) return;
        SolvRunning = true;
        foreach (var s in SolvStages) s.State = 0;
        SolvStages[0].State = 1;
        SolvProgress = 0.05;
        SolvStageChip = "stage 1 of 3";
        var solute = SolvUseSolute ? _doc : null;
        var opts = SolvOptions();
        var title = (solute != null ? Title + " in " : "") + (_solvent?.Name.ToLowerInvariant() ?? "solvent");
        try
        {
            var (doc, rep) = await Task.Run(() => CapsDocument.Solvate(solute, opts, (stage, loop, loops, dmin, bad) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    for (var k = 0; k < SolvStages.Length; k++) SolvStages[k].State = k < stage ? 2 : k == stage ? 1 : 0;
                    SolvStageChip = $"stage {stage + 1} of 3";
                    SolvProgress = stage == 0 ? 0.1 : stage == 1 ? 0.2 + 0.7 * Math.Min(1.0, loops > 0 ? (double)loop / loops : 0) : 0.95;
                    Status = string.Format(CultureInfo.InvariantCulture, "Solvating · {0} · min. distance {1:0.00} Å · {2} molecules too close", SolvStages[stage].Title, dmin, bad);
                });
                return true;
            }, title));
            foreach (var s in SolvStages) s.State = 2;
            SolvProgress = 1;
            SolvStageChip = "done";
            SolvLog = rep;
            Show(doc, title);
            Record($"doc = caps.build.solvate({(solute != null ? "doc" : "None")}, **{opts})");
            GrownUnsaved = true;
            Status = "Solvated · " + (rep.Split('\n').FirstOrDefault() ?? "");
            SetModule(8);
        }
        catch (Exception e)
        {
            foreach (var s in SolvStages) if (s.State == 1) s.State = 3;
            SolvError = e.Message;
            Status = "Could not solvate: " + e.Message;
        }
        finally { SolvRunning = false; }
    }
}
