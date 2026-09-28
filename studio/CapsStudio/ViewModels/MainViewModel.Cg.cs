using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Coarse-grained melts (design/boards/CoarseGrained): Kremer–Grest bead-spring chains built as random walks at a
/// reduced density, drawn, and written for LAMMPS with the push-off and the FENE + WCA run.</summary>
public sealed partial class MainViewModel
{
    public bool IsCg => _module == 44;
    private decimal _cgChains = 50, _cgBeads = 100, _cgDensity = 0.85m, _cgKTheta = 0, _cgSeed = 1, _cgPush = 20000, _cgRun = 100000;
    public decimal CgChains { get => _cgChains; set { if (Set(ref _cgChains, Math.Clamp(Math.Round(value), 1, 20000))) CgPreview(); } }
    public decimal CgBeads { get => _cgBeads; set { if (Set(ref _cgBeads, Math.Clamp(Math.Round(value), 2, 100000))) CgPreview(); } }
    public decimal CgDensity { get => _cgDensity; set { if (Set(ref _cgDensity, Math.Clamp(value, 0.06m, 1.2m))) { CgPreview(); Raise(nameof(CgMapText)); } } }
    public decimal CgKTheta { get => _cgKTheta; set { if (Set(ref _cgKTheta, Math.Clamp(value, 0, 50))) CgPreview(); } }
    public decimal CgSeed { get => _cgSeed; set { if (Set(ref _cgSeed, Math.Max(1, Math.Round(value)))) CgPreview(); } }
    public decimal CgPushSteps { get => _cgPush; set => Set(ref _cgPush, Math.Clamp(Math.Round(value), 100, 1e8m)); }
    public decimal CgRunSteps { get => _cgRun; set => Set(ref _cgRun, Math.Clamp(Math.Round(value), 100, 1e10m)); }
    public string CgBoxText => ((double)(_cgChains * _cgBeads / _cgDensity) is var v && v > 0 ? Math.Cbrt(v) : 0).ToString("0.00", CultureInfo.InvariantCulture);
    private string _cgHud = "", _cgInfo = "", _cgError = "";
    public string CgHud { get => _cgHud; private set => Set(ref _cgHud, value); }
    public string CgInfo { get => _cgInfo; private set => Set(ref _cgInfo, value); }
    public string CgError { get => _cgError; private set => Set(ref _cgError, value); }
    public string CgFooter => $"{_cgChains * _cgBeads:N0} beads · L = {CgBoxText} σ (from ρσ³ = {_cgDensity.ToString("0.00", CultureInfo.InvariantCulture)})";
    private CapsDocument? _cgDoc;
    public CapsDocument? CgDoc { get => _cgDoc; private set { var old = _cgDoc; if (Set(ref _cgDoc, value)) { CgViewChanged?.Invoke(); old?.Dispose(); } } }
    public event Action? CgViewChanged;
    private int _cgTicket;

    public string CgOptions()
    {
        var o = new JsonObject { ["chains"] = (int)_cgChains, ["beads"] = (int)_cgBeads, ["density"] = (double)_cgDensity, ["k_theta"] = (double)_cgKTheta, ["seed"] = (int)_cgSeed };
        if (_cgUnits == 1) { o["sigma"] = (double)_cgSigma; o["temperature"] = (double)_cgTemp; o["bead_mass"] = (double)_cgMass; }
        return o.ToJsonString();
    }

    /// <summary>Units of the LAMMPS files: reduced (ε, σ, m), or real with the user's mapping — σ in Å, the temperature that
    /// sets ε = k_B T, the bead mass — every length, energy, mass and time scaled exactly (τ = σ √(m/ε)).</summary>
    public static readonly string[] CgUnitModes = ["Reduced (ε, σ, m)", "Real · mapped by σ, T and bead mass"];
    private int _cgUnits;
    private decimal _cgSigma = 5.0m, _cgTemp = 300m, _cgMass = 50m;   // round placeholders: the mapping is the user's
    public int CgUnits { get => _cgUnits; set { if (Set(ref _cgUnits, Math.Clamp(value, 0, 1))) { Raise(nameof(CgMapped)); Raise(nameof(CgMapText)); } } }
    public bool CgMapped => _cgUnits == 1;
    public decimal CgSigma { get => _cgSigma; set { if (Set(ref _cgSigma, Math.Clamp(value, 0.5m, 50m))) Raise(nameof(CgMapText)); } }
    public decimal CgTemp { get => _cgTemp; set { if (Set(ref _cgTemp, Math.Clamp(value, 1m, 5000m))) Raise(nameof(CgMapText)); } }
    public decimal CgMass { get => _cgMass; set { if (Set(ref _cgMass, Math.Clamp(value, 1m, 100000m))) Raise(nameof(CgMapText)); } }
    /// <summary>ε, τ and the mass density the mapping gives (LAMMPS units real: k_B = 0.0019872067 kcal/mol/K, 48.888 fs per time unit).</summary>
    public string CgMapText
    {
        get
        {
            var eps = 0.0019872067 * (double)_cgTemp;
            var tau = (double)_cgSigma * Math.Sqrt((double)_cgMass / eps) * 48.88821291;
            var rho = (double)_cgDensity * (double)_cgMass / 6.02214076e23 / Math.Pow((double)_cgSigma * 1e-8, 3);
            return string.Format(CultureInfo.InvariantCulture, "ε = k_B T = {0:0.####} kcal/mol · τ = σ√(m/ε) = {1:0.##} ps · Δt = 0.01 τ = {2:0.#} fs · ρ = {3:0.###} g/cm³ at ρσ³ = {4:0.00}. The mapping is yours to choose (e.g. from Everaers et al. 2020's Kuhn-scale mapping of a polymer); CAPS rescales exactly.",
                                 eps, tau / 1000, tau * 0.01, rho, _cgDensity);
        }
    }

    // ---------------------------------------------------------------- the model: Kremer–Grest or MARTINI
    /// <summary>0 Kremer–Grest (generic bead-spring, reduced or mapped units), 1 MARTINI (chemistry-specific beads, 4 heavy
    /// atoms per bead, its own force field file).</summary>
    private int _cgModel;
    public int CgModel { get => _cgModel; set { if (Set(ref _cgModel, Math.Clamp(value, 0, 1))) { Raise(nameof(CgIsKg)); Raise(nameof(CgIsMartini)); Raise(nameof(CgBuildText)); Raise(nameof(CgModelChip)); } } }
    public bool CgIsKg => _cgModel == 0;
    public bool CgIsMartini => _cgModel == 1;
    public string CgBuildText => _cgModel == 1 ? "Build MARTINI melt" : "Build CG melt";
    public string CgModelChip => _cgModel == 1 ? "MARTINI" : "Kremer–Grest";
    /// <summary>The MARTINI force fields of the catalogue with bead typing (MARTINI 2.0 and its polymer, lipid … extensions,
    /// Martini 3), read from data/forcefields/catalogue.json: the Field page lists one entry per force field, the melt
    /// builder needs the MARTINI 2 files too.</summary>
    private List<FfEntry>? _mtFfs;
    public List<FfEntry> MtForceFields => _mtFfs ??= LoadMartiniForceFields();
    private static List<FfEntry> LoadMartiniForceFields()
    {
        var r = new List<FfEntry>();
        if (Paths.ForceFields is not { } dir) return r;
        try
        {
            using var js = System.Text.Json.JsonDocument.Parse(System.IO.File.ReadAllText(System.IO.Path.Combine(dir, "catalogue.json")));
            foreach (var e in js.RootElement.GetProperty("forcefields").EnumerateArray())
            {
                string S(string k) => e.TryGetProperty(k, out var v) && v.ValueKind == System.Text.Json.JsonValueKind.String ? v.GetString()! : "";
                var id = S("id");
                if (!id.StartsWith("martini", StringComparison.Ordinal) || id is "martini-aminoacids" or "martini22-proteins") continue;   // proteins: the Biomolecule builder
                if (!e.TryGetProperty("typing", out var t) || t.ValueKind != System.Text.Json.JsonValueKind.Object) continue;
                r.Add(new FfEntry(id, S("name"), S("version"), S("status"), System.IO.Path.Combine(dir, S("file")), true));
            }
        }
        catch { }
        // MARTINI 2.0 first: the polymer mappings in the examples are MARTINI 2
        return r.OrderBy(x => x.Id == "martini-moltemplate" ? 0 : x.Id == "martini-polymers" ? 1 : 2).ToList();
    }
    private int _mtFf;
    public int MtFfIndex { get => _mtFf; set => Set(ref _mtFf, value); }
    public static readonly string[] MtExamples = ["[C1] · 4 CH₂ per bead (polyethylene-like; MARTINI 2)", "[SN0] · poly(ethylene oxide) (Lee et al. 2009; MARTINI 2 polymers)",
                                                   "[C1][C3] · two apolar beads per unit (write your own mapping)"];
    private string _mtRepeat = "[C1]";
    public string MtRepeat { get => _mtRepeat; set => Set(ref _mtRepeat, value ?? ""); }
    public int MtExample { get => -1; set { if (value >= 0 && value < MtExamples.Length) { MtRepeat = MtExamples[value].Split(' ')[0]; if (value == 1) SelectMartini("martini-polymers"); } } }
    private void SelectMartini(string id) { var k = MtForceFields.FindIndex(e => e.Id == id); if (k >= 0) MtFfIndex = k; }
    private decimal _mtRepeats = 20, _mtChains = 40, _mtDensity = 0.9m, _mtSeed = 1;
    public decimal MtRepeats { get => _mtRepeats; set => Set(ref _mtRepeats, Math.Clamp(Math.Round(value), 1, 5000)); }
    public decimal MtChains { get => _mtChains; set => Set(ref _mtChains, Math.Clamp(Math.Round(value), 1, 20000)); }
    public decimal MtDensity { get => _mtDensity; set => Set(ref _mtDensity, Math.Clamp(value, 0.1m, 3m)); }
    public decimal MtSeed { get => _mtSeed; set => Set(ref _mtSeed, Math.Max(1, Math.Round(value))); }
    private bool _mtBuilding;
    public bool MtBuilding { get => _mtBuilding; private set { if (Set(ref _mtBuilding, value)) Raise(nameof(CgIdle)); } }
    public bool CgIdle => !_mtBuilding;

    /// <summary>Builds the MARTINI melt (packed, assigned, compressed) and opens it in the Studio, ready to export.</summary>
    public async Task BuildMartini()
    {
        var ffs = MtForceFields;
        if (_mtFf < 0 || _mtFf >= ffs.Count) { CgError = "Choose a MARTINI force field"; return; }
        var ff = ffs[_mtFf];
        var opts = new JsonObject { ["forcefield"] = ff.File, ["repeat"] = _mtRepeat.Trim(), ["repeats"] = (int)_mtRepeats, ["chains"] = (int)_mtChains,
                                    ["density"] = (double)_mtDensity, ["seed"] = (int)_mtSeed }.ToJsonString();
        CgError = "";
        MtBuilding = true;
        Status = $"Building the MARTINI melt · {_mtChains:0} chains of {_mtRepeats:0} × {_mtRepeat.Trim()}";
        try
        {
            var (doc, rep) = await Task.Run(() => CapsDocument.MartiniMelt(opts, "MARTINI melt"));
            Show(doc, $"MARTINI_{_mtRepeat.Trim().Trim('[', ']')}_{_mtChains:0}x{_mtRepeats:0}");
            GrownUnsaved = true;
            var j = JsonNode.Parse(rep)!;
            Status = string.Format(CultureInfo.InvariantCulture, "MARTINI melt built · {0} beads at {1:0.000} g/cm³ with {2} · export it from Export (GROMACS or LAMMPS)",
                                   (double?)j["beads"], (double?)j["density"], (string?)j["forcefield"]);
            SetModule(8);
        }
        catch (Exception e) { CgError = e.Message; Status = "Could not build the MARTINI melt: " + e.Message; }
        finally { MtBuilding = false; }
    }

    public void OpenCg()
    {
        SetModule(44);
        CgPreview();
    }

    /// <summary>The melt as built (instant: random walks), for the preview.</summary>
    public void CgPreview()
    {
        Raise(nameof(CgBoxText));
        Raise(nameof(CgFooter));
        if (!IsCg) return;
        var ticket = ++_cgTicket;
        var opts = CgOptions();
        Task.Run(() =>
        {
            try { var (d, r) = CapsDocument.KgBuild(opts, "Kremer–Grest melt"); return (Doc: (CapsDocument?)d, Report: r, Error: (string?)null); }
            catch (Exception e) { return (Doc: (CapsDocument?)null, Report: "", Error: (string?)e.Message); }
        }).ContinueWith(t => Avalonia.Threading.Dispatcher.UIThread.Post(() =>
        {
            var (doc, rep, err) = t.Result;
            if (ticket != _cgTicket) { doc?.Dispose(); return; }
            CgError = err ?? "";
            if (doc == null) return;
            CgDoc = doc;
            var inv = CultureInfo.InvariantCulture;
            var j = JsonNode.Parse(rep)!;
            CgHud = $"{_cgChains:0} × {_cgBeads:0} beads · bond 0.97 σ";
            CgInfo = $"⟨R²⟩/(N−1) = {j["r2_per_bond"]!.GetValue<double>().ToString("0.00", inv)} σ² per bond · closest pair {j["closest"]!.GetValue<double>().ToString("0.00", inv)} σ — random walks overlap until the push-off opens them";
        }));
    }

    /// <summary>The melt as the Studio document.</summary>
    public void BuildCg()
    {
        try
        {
            var (doc, _) = CapsDocument.KgBuild(CgOptions(), "Kremer–Grest melt");
            Show(doc, $"KG_melt_{_cgChains:0}x{_cgBeads:0}");
            GrownUnsaved = true;
            Status = "Kremer–Grest melt built · export it for LAMMPS from Builders › Coarse-grained (σ = 1 Å in CAPS files)";
            SetModule(8);
        }
        catch (Exception e) { CgError = e.Message; }
    }

    public string ExportCg(string stem)
    {
        if (_cgDoc == null) return "Nothing built";
        _cgDoc.KgLammps(CgOptions(), stem, (double)_cgPush, (double)_cgRun);
        return $"Wrote {System.IO.Path.GetFileName(stem)}.data and {System.IO.Path.GetFileName(stem)}.in · run: lmp -in {System.IO.Path.GetFileName(stem)}.in";
    }
}
