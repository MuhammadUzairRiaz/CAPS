using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Build › Semicrystalline stack (core semicrystal.hpp; caps semicrystalline): a crystal lamella of straight stems
/// and an amorphous layer of loops, bridges and tails, grown at the amorphous density; the report (counts, lengths,
/// crystallinity, tilt, density and P2 profiles) and the view coloured by the crystal and amorphous groups.
/// Modules › Thermal conductivity (core thermal.hpp; caps thermal): Green–Kubo or reverse NEMD on the open structure,
/// 2D or 3D said by the user, and the LAMMPS inputs for the long run on a cluster.</summary>
public sealed partial class MainViewModel
{
    public bool IsSemicrystalline => _module == 81;
    public bool IsThermal => _module == 82;
    public void OpenSemicrystalline() => SetModule(81);
    public void OpenThermal() => SetModule(82);

    // ---- the builder
    private decimal _scNx = 4, _scNy = 6, _scCrystal = 5, _scAmorphous = 4, _scLoops = 0.5m, _scBridges = 0.2m, _scTails = 0.3m, _scDensity = 0.855m;
    private int _scReentry;
    private bool _scUnited, _scRelax = true, _scRunning;
    private string _scError = "", _scProgress = "", _scSummary = "", _scNotes = "";
    public static readonly string[] ScReentries = ["Adjacent re-entry (mostly neighbouring stems)", "Random re-entry (switchboard)"];
    public decimal ScNx { get => _scNx; set { if (Set(ref _scNx, Math.Clamp(Math.Round(value), 1, 20))) Raise(nameof(ScEstimate)); } }
    public decimal ScNy { get => _scNy; set { if (Set(ref _scNy, Math.Clamp(Math.Round(value), 1, 20))) Raise(nameof(ScEstimate)); } }
    public decimal ScCrystal { get => _scCrystal; set { if (Set(ref _scCrystal, Math.Clamp(value, 0.5m, 50))) Raise(nameof(ScEstimate)); } }
    public decimal ScAmorphous { get => _scAmorphous; set { if (Set(ref _scAmorphous, Math.Clamp(value, 0.5m, 50))) Raise(nameof(ScEstimate)); } }
    public decimal ScLoops { get => _scLoops; set => Set(ref _scLoops, Math.Clamp(value, 0, 1)); }
    public decimal ScBridges { get => _scBridges; set => Set(ref _scBridges, Math.Clamp(value, 0, 1)); }
    public decimal ScTails { get => _scTails; set => Set(ref _scTails, Math.Clamp(value, 0, 1)); }
    public decimal ScDensity { get => _scDensity; set { if (Set(ref _scDensity, Math.Clamp(value, 0.3m, 1.2m))) Raise(nameof(ScEstimate)); } }
    public int ScReentry { get => _scReentry; set => Set(ref _scReentry, Math.Clamp(value, 0, 1)); }
    public bool ScUnited { get => _scUnited; set { if (Set(ref _scUnited, value)) Raise(nameof(ScEstimate)); } }
    public bool ScRelax { get => _scRelax; set => Set(ref _scRelax, value); }
    public bool ScRunning { get => _scRunning; private set { if (Set(ref _scRunning, value)) RaiseBusy(); } }
    public string ScError { get => _scError; private set { if (Set(ref _scError, value)) Raise(nameof(ScHasError)); } }
    public bool ScHasError => _scError.Length > 0;
    public string ScProgress { get => _scProgress; private set => Set(ref _scProgress, value); }
    public string ScSummary { get => _scSummary; private set => Set(ref _scSummary, value); }
    public string ScNotes { get => _scNotes; private set => Set(ref _scNotes, value); }
    public SeedChoice ScSeed => _scSeedC ??= new SeedChoice();
    private SeedChoice? _scSeedC;
    public (double X, double Y)[] ScRho { get; private set; } = [];
    public (double X, double Y)[] ScP2 { get; private set; } = [];
    public (double X, double Y)[] ScBridgeHist { get; private set; } = [];
    public event Action? ScChanged;

    /// <summary>The atoms the cell will hold, before it is built (PE: 3 atoms per CH2 all-atom, 1 united atom).</summary>
    public string ScEstimate
    {
        get
        {
            const double a = 7.40, b = 4.93, c = 2.534;
            var nc = Math.Max(2, Math.Round((double)_scCrystal * 10 / c));
            var stems = 2 * (double)_scNx * (double)_scNy;
            var crys = stems * 2 * nc;
            var am = (double)_scDensity * a * (double)_scNx * b * (double)_scNy * (double)_scAmorphous * 10 / 1.66053906660 / 14.027;
            var atoms = (_scUnited ? 1 : 3) * (crys + am);
            return string.Format(CultureInfo.InvariantCulture, "{0:0} stems · {1:0} cells along c · about {2:N0} atoms{3}", stems, nc, atoms,
                atoms > 15000 ? " — above the 15,000-atom limit: fewer repeats, thinner layers or united atom" : "");
        }
    }

    public async Task BuildSemicrystalline()
    {
        if (!Idle) return;
        var inv = CultureInfo.InvariantCulture;
        var json = new JsonObject
        {
            ["polymer"] = "pe", ["nx"] = (int)_scNx, ["ny"] = (int)_scNy, ["crystal_nm"] = (double)_scCrystal, ["amorphous_nm"] = (double)_scAmorphous,
            ["loops"] = (double)_scLoops, ["bridges"] = (double)_scBridges, ["tails"] = (double)_scTails, ["reentry"] = _scReentry == 0 ? "adjacent" : "random",
            ["density"] = (double)_scDensity, ["seed"] = ScSeed.Take(), ["united_atom"] = _scUnited, ["relax"] = _scRelax,
        }.ToJsonString();
        ScRunning = true;
        ScError = "";
        Status = "Building the lamellar stack…";
        try
        {
            var title = $"PE lamellar stack {(int)_scNx}×{(int)_scNy}";
            var (doc, rep) = await Task.Run(() => CapsDocument.Semicrystalline(json, f =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => ScProgress = $"growing the amorphous layer · {f * 100:0} %");
                return true;
            }, title));
            ScRunning = false;
            Show(doc, title);
            GrownUnsaved = true;
            SetModule(81);
            ShowScReport(rep);
            ColourByPhase();
            Status = "Semicrystalline stack built: " + ScSummary.Split('\n')[0];
        }
        catch (Exception e) { ScError = e.Message; Status = "Semicrystalline: " + e.Message; }
        finally { ScRunning = false; ScProgress = ""; }
    }

    /// <summary>The report of the open structure (any frame of a lamellar cell, after MD too).</summary>
    public void AnalyseSemicrystalline()
    {
        if (_doc == null) { Status = "Open a structure first"; return; }
        try { ShowScReport(_doc.SemicrystallineReport()); }
        catch (Exception e) { ScError = e.Message; }
    }

    private void ShowScReport(string rep)
    {
        var inv = CultureInfo.InvariantCulture;
        var r = JsonNode.Parse(rep)!;
        double D(string k) => r[k]?.GetValue<double>() ?? 0;
        double[] A(string k) => r[k] is JsonArray a ? a.Select(x => x!.GetValue<double>()).ToArray() : [];
        var segs = Math.Max(1, D("loops") + D("bridges") + D("tails"));
        double Mean(double[] v) => v.Length > 0 ? v.Average() : 0;
        ScSummary = string.Format(inv,
            "{0:N0} atoms · {1:0} stems · {2:0} chains\n{3:0} loops ({4:0} %) · {5:0} bridges ({6:0} %) · {7:0} tails ({8:0} %)\n" +
            "mean CH2: loops {9:0.0} · bridges {10:0.0} · tails {11:0.0}\ncrystallinity {12:0.000} (P2 > {13:0.00}) · {14:0.000} by mass · stem tilt {15:0.0} ± {16:0.0}°\n" +
            "density {17:0.000} g/cm³ · amorphous layer {18:0.000} g/cm³",
            D("atoms"), D("stems"), D("chains"), D("loops"), 100 * D("loops") / segs, D("bridges"), 100 * D("bridges") / segs, D("tails"), 100 * D("tails") / segs,
            Mean(A("loop_units")), Mean(A("bridge_units")), Mean(A("tail_units")), D("crystallinity_p2"), D("p2_cut"), D("crystallinity_mass"), D("tilt_deg"), D("tilt_sd"),
            D("density"), D("amorphous_density"));
        ScNotes = r["notes"] is JsonArray n ? string.Join(" · ", n.Select(x => (string?)x ?? "")) : "";
        var z = A("z");
        ScRho = z.Zip(A("rho")).Select(t => (t.First, t.Second)).ToArray();
        ScP2 = z.Zip(A("p2")).Select(t => (t.First, t.Second)).ToArray();
        var b = A("bridge_units");
        ScBridgeHist = b.Length == 0 ? [] : b.GroupBy(x => Math.Floor(x / 10) * 10).OrderBy(g => g.Key).Select(g => (g.Key + 5, (double)g.Count())).ToArray();
        ScChanged?.Invoke();
    }

    /// <summary>The view coloured by the structure's tags (crystal, amorphous): the Visualize pipeline's colour coding.</summary>
    public void ColourByPhase()
    {
        var row = PipelineRows.FirstOrDefault(r => r.Type == "colour_coding");
        if (row == null)
        {
            PipeSelected = PipelineRows.FirstOrDefault();
            AddStep("colour_coding");
            row = PipelineRows.First(r => r.Type == "colour_coding");
        }
        row.Params["property"] = "Tags";
        row.Params["mode"] = "categorical";
        row.Params["map"] = "";
        row.Params["lighten_h"] = false;
        PipeSelected = null;
        PipeSelected = row;
        Status = "Colour by tag: crystal and amorphous";
    }

    // ---- thermal conductivity
    private int _thDim = -1, _thMethod, _thAxis = 2;
    private decimal _thT = 300, _thEq = 2, _thPs = 3, _thDt = 1, _thEvery = 5, _thCorr = 1, _thSlabs = 20, _thSwap = 100, _thThick, _thLmpPs = 2000;
    private bool _thRunning;
    private string _thError = "", _thProgress = "", _thResult = "", _thNotes = "";
    public static readonly string[] ThDims = ["3D bulk (κ along x, y, z)", "2D sheet (vacuum along z: in-plane κ, sheet conductance)"];
    public static readonly string[] ThMethods = ["Green–Kubo (heat-flux autocorrelation, NVE)", "Reverse NEMD (Müller-Plathe exchanges)"];
    public static readonly string[] ThAxes = ["x", "y", "z"];
    /// <summary>−1 until the user says 2D or 3D: CAPS does not guess it.</summary>
    public int ThDim { get => _thDim; set { if (Set(ref _thDim, Math.Clamp(value, -1, 1))) { Raise(nameof(ThIs2D)); Raise(nameof(ThCanRun)); } } }
    public bool ThIs2D => _thDim == 1;
    public bool ThCanRun => _thDim >= 0;
    public int ThMethod { get => _thMethod; set { if (Set(ref _thMethod, Math.Clamp(value, 0, 1))) Raise(nameof(ThNemd)); } }
    public bool ThNemd => _thMethod == 1;
    public int ThAxis { get => _thAxis; set => Set(ref _thAxis, Math.Clamp(value, 0, 2)); }
    public decimal ThT { get => _thT; set => Set(ref _thT, Math.Clamp(value, 1, 5000)); }
    public decimal ThEq { get => _thEq; set => Set(ref _thEq, Math.Clamp(value, 0, 100000)); }
    public decimal ThPs { get => _thPs; set => Set(ref _thPs, Math.Clamp(value, 0.1m, 100000)); }
    public decimal ThDt { get => _thDt; set => Set(ref _thDt, Math.Clamp(value, 0.1m, 5)); }
    public decimal ThEvery { get => _thEvery; set => Set(ref _thEvery, Math.Clamp(Math.Round(value), 1, 1000)); }
    public decimal ThCorr { get => _thCorr; set => Set(ref _thCorr, Math.Clamp(value, 0.01m, 1000)); }
    public decimal ThSlabs { get => _thSlabs; set => Set(ref _thSlabs, Math.Clamp(Math.Round(value / 2) * 2, 6, 200)); }
    public decimal ThSwap { get => _thSwap; set => Set(ref _thSwap, Math.Clamp(Math.Round(value), 1, 100000)); }
    public decimal ThThick { get => _thThick; set => Set(ref _thThick, Math.Clamp(value, 0, 100)); }
    public decimal ThLmpPs { get => _thLmpPs; set => Set(ref _thLmpPs, Math.Clamp(value, 1, 1000000)); }
    public bool ThRunning { get => _thRunning; private set { if (Set(ref _thRunning, value)) RaiseBusy(); } }
    public string ThError { get => _thError; private set { if (Set(ref _thError, value)) Raise(nameof(ThHasError)); } }
    public bool ThHasError => _thError.Length > 0;
    public string ThProgress { get => _thProgress; private set => Set(ref _thProgress, value); }
    public string ThResult { get => _thResult; private set => Set(ref _thResult, value); }
    public string ThNotes { get => _thNotes; private set => Set(ref _thNotes, value); }
    public string ThSteps => string.Format(CultureInfo.InvariantCulture, "{0:N0} steps here", Math.Round((_thEq + _thPs) * 1000 / _thDt));
    public (double X, double Y)[][] ThCurves { get; private set; } = [];
    public event Action? ThChanged;

    private string ThermalJson() => new JsonObject
    {
        ["dimensionality"] = _thDim == 1 ? 2 : 3, ["thickness"] = (double)_thThick, ["method"] = _thMethod == 0 ? "gk" : "nemd", ["temperature"] = (double)_thT,
        ["equilibrate_ps"] = (double)_thEq, ["ps"] = (double)_thPs, ["dt"] = (double)_thDt, ["sample_every"] = (int)_thEvery, ["corr_ps"] = (double)_thCorr,
        ["axis"] = _thAxis, ["slabs"] = (int)_thSlabs, ["swap_every"] = (int)_thSwap, ["lammps_ps"] = (double)_thLmpPs, ["lammps_equilibrate_ps"] = 100.0,
    }.ToJsonString();

    public async Task RunThermal()
    {
        if (!Idle || _doc == null) { if (_doc == null) Status = "Open a structure first"; return; }
        if (_thDim < 0) { ThError = "Say whether the material is 3D (bulk) or 2D (a sheet with vacuum along z)"; return; }
        var doc = _doc;
        var json = ThermalJson();
        ThRunning = true;
        ThError = "";
        Status = "Thermal conductivity: running…";
        try
        {
            var res = await Task.Run(() => doc.Thermal(json, f =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => ThProgress = $"{f * 100:0} %");
                return true;
            }));
            var r = JsonNode.Parse(res)!;
            ThResult = ((string?)r["summary"] ?? "").Trim();
            ThNotes = r["notes"] is JsonArray n ? string.Join(" · ", n.Select(x => (string?)x ?? "")) : "";
            double[] A(JsonNode? a) => a is JsonArray v ? v.Select(x => x!.GetValue<double>()).ToArray() : [];
            if ((string?)r["method"] == "gk")
            {
                var t = A(r["t_ps"]);
                ThCurves = ((JsonArray)r["running"]!).Select(run => t.Zip(A(run)).Select(p => (p.First, p.Second)).ToArray()).Where(c => c.Length > 0).ToArray();
            }
            else ThCurves = [A(r["slab_x"]).Zip(A(r["slab_t"])).Select(p => (p.First, p.Second)).ToArray()];
            Status = "Thermal conductivity: " + ThResult;
            ThChanged?.Invoke();
        }
        catch (Exception e) { ThError = e.Message; Status = "Thermal conductivity: " + e.Message; }
        finally { ThRunning = false; ThProgress = ""; }
    }

    public string WriteThermalLammps(string path)
    {
        if (_doc == null) return "Open a structure first";
        if (_thDim < 0) { ThError = "Say whether the material is 3D (bulk) or 2D (a sheet with vacuum along z)"; return ThError; }
        try { return _doc.ThermalLammps(ThermalJson(), path); }
        catch (Exception e) { ThError = e.Message; return "Could not write: " + e.Message; }
    }
}
