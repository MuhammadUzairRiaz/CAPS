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
    public decimal CgSeed { get => _cgSeed; set => CgSeedChoice.Value = value; }
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
        var o = new JsonObject { ["chains"] = (int)_cgChains, ["beads"] = (int)_cgBeads, ["density"] = (double)_cgDensity, ["k_theta"] = (double)_cgKTheta, ["seed"] = CgSeedChoice.Take() };
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
    public int CgModel { get => _cgModel; set { if (Set(ref _cgModel, Math.Clamp(value, 0, 2))) { Raise(nameof(CgIsKg)); Raise(nameof(CgIsMartini)); Raise(nameof(CgIsMapped)); Raise(nameof(CgBuildText)); Raise(nameof(CgModelChip)); } } }
    public bool CgIsKg => _cgModel == 0;
    public bool CgIsMartini => _cgModel == 1;
    public bool CgIsMapped => _cgModel == 2;
    public string CgBuildText => _cgModel switch { 1 => "Build MARTINI melt", 2 => "Build mapped CG melt", _ => "Build CG melt" };
    public string CgModelChip => _cgModel switch { 1 => "MARTINI", 2 => "From a polymer", _ => "Kremer–Grest" };

    // ---- From a polymer (structure-based, cg_map.hpp): an all-atom reference melt of the chosen polymer mapped to beads
    public static readonly string[] MpSchemes =
    [
        "1 bead per repeat unit (centre of mass)",
        "2 beads per unit: backbone + side group",
        "n backbone atoms per bead",
    ];
    private static readonly string[] MpSchemeIds = ["unit", "backbone_side", "backbone_n"];
    public static readonly string[] MpSchemeHelp =
    [
        "Each repeat unit becomes one bead at its centre of mass — the simplest chemistry-specific model (e.g. polystyrene 1:1, Milano & Müller-Plathe 2005). Fast; the side group's shape is lost.",
        "Each unit becomes a backbone bead (its backbone atoms) and a side-group bead (the rest) — the '2-bead polymer' of moltemplate's examples, and the polystyrene model of Harmandaris et al. 2006. Keeps where the side group sits; a unit without a side group stays one bead.",
        "Every n backbone atoms make a bead, side groups with them — for chains without big side groups (polyethylene 3:1: three CH₂ per bead).",
    ];
    private int _mpPolymer, _mpScheme;
    private decimal _mpPerBead = 3, _mpChains = 10, _mpDp = 20, _mpDensity = 1.0m, _mpTemp = 450;
    private string _mpLog = "";
    public int MpPolymer { get => _mpPolymer; set => Set(ref _mpPolymer, Math.Max(0, value)); }
    public int MpScheme { get => _mpScheme; set { if (Set(ref _mpScheme, Math.Clamp(value, 0, 2))) { Raise(nameof(MpHelp)); Raise(nameof(MpIsN)); } } }
    public string MpHelp => MpSchemeHelp[_mpScheme];
    public bool MpIsN => _mpScheme == 2;
    public decimal MpPerBead { get => _mpPerBead; set => Set(ref _mpPerBead, Math.Clamp(Math.Round(value), 1, 20)); }
    public decimal MpChains { get => _mpChains; set => Set(ref _mpChains, Math.Clamp(Math.Round(value), 2, 500)); }
    public decimal MpDp { get => _mpDp; set => Set(ref _mpDp, Math.Clamp(Math.Round(value), 3, 500)); }
    public decimal MpDensity { get => _mpDensity; set => Set(ref _mpDensity, Math.Clamp(value, 0.3m, 3m)); }
    public decimal MpTemp { get => _mpTemp; set => Set(ref _mpTemp, Math.Clamp(value, 50m, 2000m)); }
    public string MpLog { get => _mpLog; private set => Set(ref _mpLog, value); }
    private string MpOptions() => new JsonObject
    {
        ["scheme"] = MpSchemeIds[_mpScheme], ["per_bead"] = (int)_mpPerBead, ["temperature"] = (double)_mpTemp,
        ["chains"] = (int)_mpChains, ["dp"] = (int)_mpDp, ["density"] = (double)_mpDensity, ["seed"] = MtSeedChoice.Take(),
    }.ToJsonString();

    /// <summary>The report as text: the notes, then each bond and angle type with its inverted parameters.</summary>
    private static string MpReportText(string json)
    {
        try
        {
            var r = JsonNode.Parse(json)!;
            var inv = CultureInfo.InvariantCulture;
            var sb = new System.Text.StringBuilder();
            foreach (var n in (JsonArray)r["notes"]!) sb.Append(n!.GetValue<string>()).Append('\n');
            foreach (var b in (JsonArray)r["bonds"]!)
                sb.Append(string.Format(inv, "bond {0}: r₀ {1:0.###} Å, σ(r) {2:0.###} Å → k {3:0.###} kcal/mol/Å² ({4} samples)\n", b!["types"], b["r0"]!.GetValue<double>(),
                                        b["sd"]!.GetValue<double>(), b["k"]!.GetValue<double>(), b["count"]!.GetValue<double>()));
            foreach (var a in (JsonArray)r["angles"]!)
                sb.Append(string.Format(inv, "angle {0}: θ₀ {1:0.#}°, σ(θ) {2:0.#}° → k {3:0.###} kcal/mol/rad² ({4} samples)\n", a!["types"], a["theta0"]!.GetValue<double>(),
                                        a["sd"]!.GetValue<double>(), a["k"]!.GetValue<double>(), a["count"]!.GetValue<double>()));
            return sb.ToString().TrimEnd();
        }
        catch { return json; }
    }

    /// <summary>Builds the chosen polymer's all-atom reference melt, maps it, and opens the beads with their model.</summary>
    public async Task BuildMappedCg()
    {
        if (!Idle) return;
        var list = CgBackmapPolymers;
        if (_mpPolymer >= list.Count) { CgError = "Choose the polymer"; return; }
        var poly = list[_mpPolymer];
        var name = $"{poly.Name.Split(" (")[0]}_CG_{MpSchemeIds[_mpScheme]}";
        CgError = "";
        Status = $"Growing and compressing an all-atom {poly.Name} melt, then mapping it to beads…";
        try
        {
            var (spec, opts) = (poly.Spec, MpOptions());
            var (d, rep) = await Task.Run(() => CapsDocument.CgFromPolymer(spec, opts, name));
            Show(d, name);
            GrownUnsaved = true;
            MpLog = MpReportText(rep);
            Status = "Mapped CG melt built · its bead model is assigned (Field shows it); relax, then run or export it";
        }
        catch (Exception e) { CgError = e.Message; Status = "Could not build the mapped model: " + e.Message; }
    }

    /// <summary>The open all-atom structure (every frame of its trajectory) mapped to beads.</summary>
    public async Task MapOpenStructure()
    {
        if (_doc == null || !Idle) { CgError = "Open an all-atom cell (a grown melt, or its trajectory) first"; return; }
        var doc = _doc;
        var name = Title.Replace(" (unsaved)", "") + $" (CG, {MpSchemeIds[_mpScheme]})";
        CgError = "";
        try
        {
            var opts = MpOptions();
            var (d, rep) = await Task.Run(() => doc.CgMap(opts, name));
            Show(d, name);
            GrownUnsaved = true;
            MpLog = MpReportText(rep);
            Status = "Mapped to beads with their model · the all-atom structure stays open in the project";
        }
        catch (Exception e) { CgError = e.Message; Status = "Could not map: " + e.Message; }
    }
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
    public decimal MtSeed { get => _mtSeed; set => MtSeedChoice.Value = value; }
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
                                    ["density"] = (double)_mtDensity, ["seed"] = MtSeedChoice.Take() }.ToJsonString();
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

    // Backmap to all atoms (design/boards/CoarseGrained): one repeat unit per bead, then GAFF2 (else GAFF) on the result
    private int _cgBackUnit;
    private List<FilmPolymer>? _cgBackPolymers;
    /// <summary>Every homopolymer of the library (the unit each bead stands for), polystyrene first.</summary>
    public List<FilmPolymer> CgBackmapPolymers
    {
        get
        {
            if (_cgBackPolymers != null) return _cgBackPolymers;
            LoadPolymerLibrary();
            _cgBackPolymers = PolymerLibrary.Where(p => !p.Copolymer && p.Smiles.Count(c => c == '*') == 2)
                .OrderBy(p => p.Name == "Polystyrene" ? 0 : 1).ThenBy(p => p.Name)
                .Select(p => new FilmPolymer(p.Name, new JsonObject
                {
                    ["units"] = new JsonArray(new JsonObject { ["name"] = p.Name, ["smiles"] = p.Smiles }),
                    ["sequence"] = "homopolymer",
                }.ToJsonString())).ToList();
            return _cgBackPolymers;
        }
    }
    public int CgBackmapUnit { get => _cgBackUnit; set => Set(ref _cgBackUnit, Math.Max(0, value)); }
    private string _cgBackLog = "";
    public string CgBackmapLog { get => _cgBackLog; private set => Set(ref _cgBackLog, value); }
    public async Task BackmapCg()
    {
        if (_doc == null || !Idle) return;
        if (_cgBackUnit >= CgBackmapPolymers.Count) { CgError = "Choose the repeat unit each bead stands for"; return; }
        var poly = CgBackmapPolymers[_cgBackUnit];
        var doc = _doc;
        var name = Title.Replace(" (unsaved)", "") + $" → {poly.Name.Split(" (")[0]}, all atoms";
        CgError = "";
        Status = "Backmapping the beads to " + poly.Name + " and relaxing…";
        try
        {
            var (d, report) = await Task.Run(() => doc.KgBackmap(poly.Spec, "{}", name));
            Show(d, name);
            GrownUnsaved = true;
            CgBackmapLog = report;
            var ff = Field.Library.ToList().FindIndex(x => x.Id == "gaff2-moltemplate");
            if (ff < 0) ff = Field.Library.ToList().FindIndex(x => x.Id == "gaff-amber25");
            if (ff >= 0) { Field.FfIndex = ff; await Field.Assign(); }
            Status = "Backmapped: " + report.Split('\n').FirstOrDefault() + (ff >= 0 ? " · " + Field.ForceFieldName + " assigned" : "");
        }
        catch (Exception e) { CgError = "Backmap: " + e.Message; Status = CgError; }
    }

    public string ExportCg(string stem)
    {
        if (_cgDoc == null) return "Nothing built";
        _cgDoc.KgLammps(CgOptions(), stem, (double)_cgPush, (double)_cgRun);
        return $"Wrote {System.IO.Path.GetFileName(stem)}.data and {System.IO.Path.GetFileName(stem)}.in · run: lmp -in {System.IO.Path.GetFileName(stem)}.in";
    }
}
