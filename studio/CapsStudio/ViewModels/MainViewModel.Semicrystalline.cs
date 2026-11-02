using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Build › Semicrystalline stack (core semicrystal.hpp; caps semicrystalline): a crystal lamella of the polymer
/// crystal library's cells (the ones data/crystals/polymers/semicrystalline.json lists) and an amorphous layer of loops,
/// bridges and tails grown from the crystal's repeat unit at the amorphous density, relaxed with the polymer's force
/// field; the report (counts, lengths, crystallinity, tilt, density and P2 profiles) and the view coloured by the crystal
/// and amorphous groups.
/// Modules › Thermal conductivity (core thermal.hpp; caps thermal): Green–Kubo or reverse NEMD on the open structure,
/// 2D or 3D said by the user, and the LAMMPS inputs for the long run on a cluster.</summary>
public sealed partial class MainViewModel
{
    public bool IsSemicrystalline => _module == 81;
    public bool IsThermal => _module == 82;
    public void OpenSemicrystalline() => SetModule(81);
    public void OpenThermal() => SetModule(82);

    // ---- the builder
    /// <summary>A crystal the builder takes, with its defaults (semicrystalline.json) and its cell (polymers.json).</summary>
    public sealed record ScPolymer(string Id, string Name, string Source, string Forcefield, int Nx, int Ny, double Crystal, double Amorphous,
                                   double Loops, double Bridges, double Tails, double Density, string DensitySource,
                                   double Area, double Height, double AtomsPerCell, double MassPerAtom)
    {
        public override string ToString() => Name;
    }
    private static ScPolymer[]? _scPolymers;
    /// <summary>The crystals that build (data/crystals/polymers/semicrystalline.json), PE first.</summary>
    public static ScPolymer[] ScPolymers => _scPolymers ??= LoadScPolymers();
    private static ScPolymer[] LoadScPolymers()
    {
        try
        {
            if (Paths.Data is not { } d) return [];
            var dir = Path.Combine(d, "crystals", "polymers");
            var lib = JsonNode.Parse(File.ReadAllText(Path.Combine(dir, "polymers.json")))!["crystals"]!.AsArray();
            var def = JsonNode.Parse(File.ReadAllText(Path.Combine(dir, "semicrystalline.json")))!["polymers"]!.AsArray();
            var list = new List<ScPolymer>();
            foreach (var e in def)
            {
                var id = (string)e!["id"]!;
                var c = lib.FirstOrDefault(x => (string?)x!["id"] == id);
                if (c == null) continue;
                var cell = c["cell"]!;
                double G(string k) => cell[k]!.GetValue<double>();
                double a = G("a"), b = G("b"), cc = G("c");
                double ca = Math.Cos(G("alpha") * Math.PI / 180), cb = Math.Cos(G("beta") * Math.PI / 180), cg = Math.Cos(G("gamma") * Math.PI / 180);
                var vol = a * b * cc * Math.Sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg);
                // the face across the chains: the two edges other than the chain axis
                var axis = (string?)c["chain_axis"] ?? "c";
                var area = axis switch
                {
                    "a" => b * cc * Math.Sqrt(1 - ca * ca),
                    "b" => a * cc * Math.Sqrt(1 - cb * cb),
                    _ => a * b * Math.Sqrt(1 - cg * cg),
                };
                var atoms = c["atoms"]!.GetValue<double>();
                var mass = c["density"]!.GetValue<double>() * vol / 1.66053906660;
                var rep = e["repeats"]!.AsArray();
                list.Add(new ScPolymer(id, (string)c["name"]!, (string?)c["source"] ?? "", (string)e["forcefield"]!, rep[0]!.GetValue<int>(), rep[1]!.GetValue<int>(),
                    e["crystal_nm"]!.GetValue<double>(), e["amorphous_nm"]!.GetValue<double>(), e["loops"]!.GetValue<double>(), e["bridges"]!.GetValue<double>(),
                    e["tails"]!.GetValue<double>(), e["density"]!.GetValue<double>(), (string?)e["density_source"] ?? "", area, vol / area, atoms, mass / atoms));
            }
            return list.ToArray();
        }
        catch { return []; }
    }
    public static readonly string[] ScForcefields = ["The polymer's own (automatic)", "OPLS-AA 2024", "PCFF", "COMPASS", "None (built-in GAFF / UFF for the relaxation)"];
    private int _scPolymer, _scFf;
    private decimal _scNx = 4, _scNy = 6, _scCrystal = 7.6m, _scAmorphous = 7.6m, _scLoops = 0.35m, _scBridges = 0.10m, _scTails = 0.55m, _scDensity = 0.855m;
    private int _scReentry;
    private bool _scUnited, _scRelax = true, _scRunning;
    private string _scError = "", _scProgress = "", _scSummary = "", _scNotes = "";
    public static readonly string[] ScReentries = ["Adjacent re-entry (mostly neighbouring stems)", "Random re-entry (switchboard)"];
    /// <summary>The crystal: choosing one fills its defaults (lateral repeats, thicknesses, fractions, density).</summary>
    public int ScPolymerIndex
    {
        get => _scPolymer;
        set
        {
            if (!Set(ref _scPolymer, Math.Clamp(value, 0, Math.Max(0, ScPolymers.Length - 1)))) return;
            if (ScPolymers.Length > 0)
            {
                var p = ScPolymers[_scPolymer];
                ScNx = p.Nx; ScNy = p.Ny; ScCrystal = (decimal)p.Crystal; ScAmorphous = (decimal)p.Amorphous;
                ScLoops = (decimal)p.Loops; ScBridges = (decimal)p.Bridges; ScTails = (decimal)p.Tails; ScDensity = (decimal)p.Density;
            }
            Raise(nameof(ScCrystalText)); Raise(nameof(ScEstimate)); Raise(nameof(ScFfText));
        }
    }
    private ScPolymer? ScCurrent => ScPolymers.Length > 0 ? ScPolymers[Math.Clamp(_scPolymer, 0, ScPolymers.Length - 1)] : null;
    public string ScCrystalText => ScCurrent is { } p ? $"{p.Source}\nAmorphous density: {p.DensitySource}" : "The polymer crystal library (data/crystals/polymers) was not found";
    public int ScFf { get => _scFf; set { if (Set(ref _scFf, Math.Clamp(value, 0, ScForcefields.Length - 1))) Raise(nameof(ScFfText)); } }
    /// <summary>The library id of the chosen force field, or "" for none.</summary>
    private string ScFfId => _scFf switch { 0 => ScCurrent?.Forcefield ?? "", 1 => "oplsaa2024", 2 => "pcff", 3 => "compass", _ => "" };
    public string ScFfText => ScFfId switch
    {
        "" => "The relaxation uses the built-in GAFF / UFF; assign a force field on the Force field page afterwards",
        var f => $"Relaxed and assigned with {FfName(f)}" + (_scFf == 0 ? " (chosen for this polymer: it types every atom of the cell)" : ""),
    };
    private static string FfName(string id) => id switch { "oplsaa2024" => "OPLS-AA 2024", "pcff" => "PCFF", "compass" => "COMPASS", _ => id };
    private static string? FfFile(string id) => id.Length == 0 || Paths.Data is not { } d ? null
        : Path.Combine(d, "forcefields", (id switch { "pcff" => "pcff-frc", "compass" => "compass-frc", _ => id }) + ".json");
    public decimal ScNx { get => _scNx; set { if (Set(ref _scNx, Math.Clamp(Math.Round(value), 1, 20))) Raise(nameof(ScEstimate)); } }
    public decimal ScNy { get => _scNy; set { if (Set(ref _scNy, Math.Clamp(Math.Round(value), 1, 20))) Raise(nameof(ScEstimate)); } }
    public decimal ScCrystal { get => _scCrystal; set { if (Set(ref _scCrystal, Math.Clamp(value, 0.5m, 50))) Raise(nameof(ScEstimate)); } }
    public decimal ScAmorphous { get => _scAmorphous; set { if (Set(ref _scAmorphous, Math.Clamp(value, 0.5m, 50))) Raise(nameof(ScEstimate)); } }
    public decimal ScLoops { get => _scLoops; set => Set(ref _scLoops, Math.Clamp(value, 0, 1)); }
    public decimal ScBridges { get => _scBridges; set => Set(ref _scBridges, Math.Clamp(value, 0, 1)); }
    public decimal ScTails { get => _scTails; set => Set(ref _scTails, Math.Clamp(value, 0, 1)); }
    public decimal ScDensity { get => _scDensity; set { if (Set(ref _scDensity, Math.Clamp(value, 0.3m, 3m))) Raise(nameof(ScEstimate)); } }
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

    /// <summary>The atoms the cell will hold, before it is built (from the crystal's cell and the amorphous density).</summary>
    public string ScEstimate
    {
        get
        {
            if (ScCurrent is not { } p) return "";
            var nc = Math.Max(2, Math.Round((double)_scCrystal * 10 / p.Height));
            var cells = (double)_scNx * (double)_scNy;
            var crys = cells * nc * p.AtomsPerCell;
            var am = (double)_scDensity * cells * p.Area * (double)_scAmorphous * 10 / 1.66053906660 / p.MassPerAtom;
            var atoms = crys + am;
            if (_scUnited) atoms *= 0.4;   // heavy atoms only: about two in five for these polymers
            return string.Format(CultureInfo.InvariantCulture, "{0:0} cells along the chains · about {1:N0} atoms", nc, atoms);
        }
    }

    public async Task BuildSemicrystalline()
    {
        if (!Idle) return;
        var inv = CultureInfo.InvariantCulture;
        var poly = ScCurrent;
        var ffId = ScFfId;
        var ffFile = FfFile(ffId);
        var json = new JsonObject
        {
            ["polymer"] = poly?.Id ?? "pe", ["data_dir"] = Paths.Data ?? "", ["nx"] = (int)_scNx, ["ny"] = (int)_scNy, ["crystal_nm"] = (double)_scCrystal, ["amorphous_nm"] = (double)_scAmorphous,
            ["loops"] = (double)_scLoops, ["bridges"] = (double)_scBridges, ["tails"] = (double)_scTails, ["reentry"] = _scReentry == 0 ? "adjacent" : "random",
            ["density"] = (double)_scDensity, ["seed"] = ScSeed.Take(), ["united_atom"] = _scUnited, ["relax"] = _scRelax,
        };
        if (ffFile != null && !_scUnited) json["forcefield"] = ffFile;
        var text = json.ToJsonString();
        ScRunning = true;
        ScError = "";
        Status = "Building the lamellar stack…";
        try
        {
            var title = $"{poly?.Name ?? "PE"} lamellar stack {(int)_scNx}×{(int)_scNy}";
            var (doc, rep) = await Task.Run(() => CapsDocument.Semicrystalline(text, f =>
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
            // the cell carries the force field it was relaxed with (the Force field page shows it; exports use it)
            if (ffFile != null && !_scUnited)
                try { doc.FieldAssign(ffFile, null, 0); }
                catch (Exception fe) { ScNotes += " · " + FfName(ffId) + " not assigned: " + fe.Message; }
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
            "mean backbone atoms: loops {9:0.0} · bridges {10:0.0} · tails {11:0.0}\ncrystallinity {12:0.000} (P2 > {13:0.00}) · {14:0.000} by mass · stem tilt {15:0.0} ± {16:0.0}°\n" +
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
    /// <summary>The curves are Green–Kubo running integrals (else the NEMD slab temperatures).</summary>
    public bool ThGreenKubo { get; private set; }
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
            ThGreenKubo = (string?)r["method"] == "gk";
            if (ThGreenKubo)
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
