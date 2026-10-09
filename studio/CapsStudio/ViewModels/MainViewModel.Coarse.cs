using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A stage of the coarse-grain workflow strip (Mapping → … → Analyse).</summary>
public sealed class CgwStageItem : ObservableObject
{
    public int Index { get; init; }
    public string Title { get; init; } = "";
    public string Sub { get; init; } = "";
    public string Number => (Index + 1).ToString(CultureInfo.InvariantCulture);
    private bool _done, _active;
    public bool Done { get => _done; set { if (Set(ref _done, value)) Raise(nameof(DotBrush)); } }
    public bool Active { get => _active; set { if (Set(ref _active, value)) { Raise(nameof(DotBrush)); Raise(nameof(TitleBrush)); Raise(nameof(BackBrush)); } } }
    public IBrush DotBrush => Tokens.Brush(_done ? "OkB" : _active ? "AccB" : "DimB");
    public IBrush TitleBrush => Tokens.Brush(_active ? "TextB" : "MutedB");
    public IBrush BackBrush => Tokens.Brush(_active ? "Bg3B" : "Bg1B");
}
/// <summary>A system given to a coarse-graining command: its map.json and the frames (or dump and log) that follow it.</summary>
public sealed record CgwSysRow(string Map, string[] Frames)
{
    public string MapName => Path.GetFileName(Map);
    public string FramesText => string.Join(" · ", Frames.Select(Path.GetFileName));
}
/// <summary>A bonded table of `caps cgfit bonded`.</summary>
public sealed record CgwBondedRow(string Kind, string Key, string Min, string Range, string Halves, string Samples, bool Flag)
{
    public IBrush FlagBrush => Tokens.Brush(Flag ? "ErrB" : "TextB");
}

/// <summary>The coarse-grain workflow (module 80): chemistry-aware mapping, bonded Boltzmann inversion, non-bonded IBI with its
/// fit and calibration, the CG melt, its equilibration deck, and the analyses (entanglements, tension, dynamics, backmapping) —
/// every action a `caps` command through caps_cg_run (caps/cg_commands), so the Studio, the CLI and Python write the same
/// files; the command line of the last action is shown to copy. Everything shown comes from the replies or the files they wrote.</summary>
public partial class MainViewModel
{
    public bool IsCgw => _module == 80;
    private static readonly CultureInfo CgInv = CultureInfo.InvariantCulture;

    // ---------------------------------------------------------------- equations (typeset by MathView)
    public const string CgwMapTex = @"\mathbf{R}_{I}=\frac{\sum_{i\in I} m_{i}\,\mathbf{r}_{i}}{\sum_{i\in I} m_{i}}";
    public const string CgwBondedTex = @"U(r)=-k_{B}T\,\ln\frac{P(r)}{r^{2}}\qquad U(\theta)=-k_{B}T\,\ln\frac{P(\theta)}{\sin\theta}\qquad U(\varphi)=-k_{B}T\,\ln P(\varphi)";
    public const string CgwIbiTex = @"U_{i+1}(r)=U_{i}(r)+k_{B}T\,\ln\frac{g_{i}(r)}{g_{\mathrm{target}}(r)}";
    public const string CgwBuildTex = @"\frac{\langle R^{2}(n)\rangle}{n}\qquad \frac{L}{\langle R_{ee}^{2}\rangle^{1/2}}\ge 1";
    public const string CgwPpaTex = @"M_{e}=N_{e}\,\bar{m}_{\mathrm{bead}}\qquad Z=\frac{N}{N_{e}}\qquad a_{pp}=\frac{\langle R_{ee}^{2}\rangle}{\langle L_{pp}\rangle}";
    public const string CgwTensionTex = @"\sigma=-\left(P_{zz}-\frac{P_{xx}+P_{yy}}{2}\right)\qquad \sigma=G_{R}\left(\lambda^{2}-\frac{1}{\lambda}\right)+\sigma_{0}";
    public const string CgwDynTex = @"g_{1}(t)=\langle|\mathbf{r}_{i}(t)-\mathbf{r}_{i}(0)|^{2}\rangle\qquad D=\lim_{t\to\infty}\frac{g_{3}(t)}{6t}";

    // ---------------------------------------------------------------- common
    public ObservableCollection<CgwStageItem> CgwStages { get; } = new()
    {
        new CgwStageItem { Index = 0, Title = "Mapping", Sub = "atoms → beads", Active = true },
        new CgwStageItem { Index = 1, Title = "Bonded", Sub = "Boltzmann inversion" },
        new CgwStageItem { Index = 2, Title = "Non-bonded", Sub = "IBI · fit · calibrate" },
        new CgwStageItem { Index = 3, Title = "Build CG melt", Sub = "chains as random walks" },
        new CgwStageItem { Index = 4, Title = "Equilibrate", Sub = "LAMMPS deck" },
        new CgwStageItem { Index = 5, Title = "Analyse", Sub = "PPA · tension · dynamics" },
    };
    private int _cgwStage;
    public int CgwStage
    {
        get => _cgwStage;
        set
        {
            if (!Set(ref _cgwStage, Math.Clamp(value, 0, 5))) return;
            foreach (var s in CgwStages) s.Active = s.Index == _cgwStage;
            foreach (var n in new[] { nameof(CgwOnMap), nameof(CgwOnBonded), nameof(CgwOnNonBonded), nameof(CgwOnBuild), nameof(CgwOnEquil), nameof(CgwOnAnalyse), nameof(CgwPrimaryText), nameof(CgwPrimaryIcon) }) Raise(n);
            CgwChanged?.Invoke();
        }
    }
    /// <summary>The header's primary action: the current stage's main command.</summary>
    public string CgwPrimaryText => _cgwStage switch { 0 => "Map", 1 => "Invert bonded", 2 => "g(r) targets", 3 => "Build melt", 4 => "Copy LAMMPS line", _ => "Entanglements" };
    public string CgwPrimaryIcon => _cgwStage switch { 0 => "atom", 1 => "chart", 2 => "chart", 3 => "grow", 4 => "copy", _ => "chart" };
    public async Task CgwPrimary()
    {
        switch (_cgwStage)
        {
            case 0: await CgwRunMap(); break;
            case 1: await CgwRunBonded(); break;
            case 2: await CgwRunTargets(); break;
            case 3: await CgwBuild(); break;
            case 5: await CgwRunPpa(); break;
        }
    }
    public bool CgwOnMap => _cgwStage == 0;
    public bool CgwOnBonded => _cgwStage == 1;
    public bool CgwOnNonBonded => _cgwStage == 2;
    public bool CgwOnBuild => _cgwStage == 3;
    public bool CgwOnEquil => _cgwStage == 4;
    public bool CgwOnAnalyse => _cgwStage == 5;

    public string CgwRoot
    {
        get => _settings.CgwRoot.Length > 0 ? _settings.CgwRoot : Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "CAPS", "cg");
        set
        {
            if (_settings.CgwRoot == value) return;
            _settings.CgwRoot = value ?? "";
            Raise();
            Changed("Coarse-grain folder");
            CgwDefaults(true);
        }
    }
    /// <summary>The folder the last file was picked from (pickers start there).</summary>
    public string CgwPickDir
    {
        get => _settings.CgwPickDir;
        set { if (value == null || _settings.CgwPickDir == value) return; _settings.CgwPickDir = value; try { _settings.Save(); } catch { } }
    }
    private string _cgwCli = "", _cgwStatus = "", _cgwReport = "";
    public string CgwCli { get => _cgwCli; private set { if (Set(ref _cgwCli, value)) Raise(nameof(HasCgwCli)); } }
    public bool HasCgwCli => _cgwCli.Length > 0;
    public string CgwStatus { get => _cgwStatus; private set => Set(ref _cgwStatus, value); }
    /// <summary>The readable report of the last action (its "text").</summary>
    public string CgwReport { get => _cgwReport; private set { if (Set(ref _cgwReport, value)) Raise(nameof(HasCgwReport)); } }
    public bool HasCgwReport => _cgwReport.Length > 0;
    private bool _cgwBusy;
    public bool CgwBusy { get => _cgwBusy; private set { if (Set(ref _cgwBusy, value)) Raise(nameof(CgwIdle)); } }
    public bool CgwIdle => !_cgwBusy;
    public event Action? CgwChanged;

    /// <summary>Runs a coarse-graining command off the UI thread; records its command line and report; null on failure (CgwStatus says why).</summary>
    public async Task<JsonObject?> Cgw(string command, JsonObject args)
    {
        if (DftDataDir.Length == 0) { CgwStatus = "CAPS data folder not found (data/cg)"; return null; }
        CgwBusy = true;
        CgwStatus = $"Running {command} …";
        try
        {
            var json = args.ToJsonString();
            var text = await Task.Run(() => CapsDocument.CgRun(command, json, DftDataDir));
            var r = JsonNode.Parse(text) as JsonObject;
            CgwCli = (string?)r?["command"] ?? "";
            CgwReport = ((string?)r?["text"] ?? "").TrimEnd();
            CgwStatus = "";
            return r;
        }
        catch (Exception e) { CgwStatus = $"{command}: {e.Message}"; return null; }
        finally { CgwBusy = false; }
    }
    private string CgwSub(string name) => Path.Combine(CgwRoot, name);
    private static double CNum(JsonNode? n, double def = double.NaN)
    {
        if (n is not JsonValue v) return def;
        if (v.TryGetValue<double>(out var d)) return d;
        return v.TryGetValue<string>(out var s) && double.TryParse(s, NumberStyles.Float, CgInv, out d) ? d : def;
    }
    private static string CTxt(JsonNode? n) => n is JsonValue v ? (v.TryGetValue<string>(out var s) ? s : v.TryGetValue<double>(out var d) ? d.ToString("0.###", CgInv) : v.ToJsonString()) : n?.ToJsonString() ?? "";
    private static string CF(double v, string f = "0.000") => double.IsFinite(v) ? v.ToString(f, CgInv) : "–";
    private static JsonArray Inputs(IEnumerable<string> xs) => new(xs.Select(x => (JsonNode?)JsonValue.Create(x)).ToArray());
    private static IEnumerable<string> SysInputs(IEnumerable<CgwSysRow> rows) => rows.SelectMany(r => new[] { r.Map }.Concat(r.Frames));
    private static string Num(decimal d) => d.ToString("0.######", CgInv);

    /// <summary>Reads a CSV (or whitespace table, '#' comments) into its header and numeric rows.</summary>
    internal static (string[] Head, List<double[]> Rows) ReadTable(string path)
    {
        var head = Array.Empty<string>();
        var rows = new List<double[]>();
        if (!File.Exists(path)) return (head, rows);
        foreach (var raw in File.ReadLines(path))
        {
            var l = raw.Trim();
            if (l.Length == 0) continue;
            if (l[0] == '#') { if (head.Length == 0) head = l.TrimStart('#').Trim().Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries); continue; }
            var parts = l.Contains(',') ? l.Split(',') : l.Split((char[]?)null, StringSplitOptions.RemoveEmptyEntries);
            var v = new double[parts.Length];
            var numeric = true;
            for (var i = 0; i < parts.Length; ++i)
                if (!double.TryParse(parts[i], NumberStyles.Float, CgInv, out v[i])) { numeric = false; break; }
            if (!numeric) { if (rows.Count == 0) head = parts.Select(p => p.Trim()).ToArray(); continue; }
            rows.Add(v);
        }
        return (head, rows);
    }

    public void OpenCgw()
    {
        CgwDefaults(false);
        SetModule(80);
        Avalonia.Threading.Dispatcher.UIThread.Post(() => { Raise(nameof(CgwRoot)); Raise(nameof(CgwPreset)); CgwChanged?.Invoke(); });
    }
    /// <summary>The rule sets of data/cg/mapping_rules.json (read when first shown, so the combo has its items before its selection).</summary>
    private void CgwLoadPresets()
    {
        try
        {
            var rules = Path.Combine(DftDataDir, "cg", "mapping_rules.json");
            if (DftDataDir.Length > 0 && File.Exists(rules))
                foreach (var p in (JsonNode.Parse(File.ReadAllText(rules))?["presets"] as JsonArray ?? []).OfType<JsonObject>())
                    _cgwPresets.Add((string?)p["id"] ?? "");
        }
        catch (Exception e) { CgwStatus = "mapping rules: " + e.Message; }
        _cgwPresets.Add(CgwCustom);
        if (!_cgwPresets.Contains(_cgwPreset)) _cgwPreset = _cgwPresets[0];
    }
    /// <summary>The output folders under the working folder (kept when the user chose others, unless the root moves).</summary>
    private void CgwDefaults(bool force)
    {
        if (force || _cgwMapOut.Length == 0) CgwMapOut = CgwSub("cg");
        if (force || _cgwBondedOut.Length == 0) CgwBondedOut = CgwSub("bonded");
        if (force || _cgwIbiOut.Length == 0) CgwIbiOut = CgwSub("ibi");
        if (force || _cgwBuildOut.Length == 0) CgwBuildOut = CgwSub("melt");
        if (force || _cgwHistory.Length == 0) CgwHistory = Path.Combine(CgwSub("calib"), "calibration.json");
    }

    // ---------------------------------------------------------------- 1 Mapping
    public const string CgwCustom = "custom (cut / names)";
    private readonly ObservableCollection<string> _cgwPresets = new();
    public ObservableCollection<string> CgwPresets { get { if (_cgwPresets.Count == 0) CgwLoadPresets(); return _cgwPresets; } }
    public ObservableCollection<string> CgwInputs { get; } = new();
    private string _cgwPreset = "ester-cut", _cgwCut = "", _cgwNames = "", _cgwMapOut = "", _cgwMapDump = "", _cgwSequence = "", _cgwMapNotes = "", _cgwTypes = "";
    public string CgwPreset { get => _cgwPreset; set { if (value != null && Set(ref _cgwPreset, value)) Raise(nameof(CgwIsCustom)); } }
    public bool CgwIsCustom => _cgwPreset == CgwCustom;
    /// <summary>Custom rules: bond SMARTS cut (comma-separated) and naming rules NAME=SMARTS.</summary>
    public string CgwCut { get => _cgwCut; set => Set(ref _cgwCut, value ?? ""); }
    public string CgwNames { get => _cgwNames; set => Set(ref _cgwNames, value ?? ""); }
    public string CgwMapOut { get => _cgwMapOut; set => Set(ref _cgwMapOut, value ?? ""); }
    /// <summary>An all-atom LAMMPS dump of the (single) structure, mapped frame by frame (optional).</summary>
    public string CgwMapDump { get => _cgwMapDump; set { if (Set(ref _cgwMapDump, value ?? "")) Raise(nameof(HasCgwMapDump)); } }
    public bool HasCgwMapDump => _cgwMapDump.Length > 0;
    public ObservableCollection<DftTableRow> CgwMapRows { get; } = new();
    public string CgwSequence { get => _cgwSequence; private set => Set(ref _cgwSequence, value); }
    public string CgwMapNotes { get => _cgwMapNotes; private set => Set(ref _cgwMapNotes, value); }
    /// <summary>The shared type list (types.json) every later stage uses.</summary>
    public string CgwTypes { get => _cgwTypes; set => Set(ref _cgwTypes, value ?? ""); }
    public string[] CgwKindCats { get; private set; } = [];
    public List<(string Name, (double X, double Y)[] Pts)> CgwKindBars { get; } = new();

    public void CgwAddInputs(IEnumerable<string> files) { foreach (var f in files) if (!CgwInputs.Contains(f)) CgwInputs.Add(f); }

    public async Task CgwRunMap()
    {
        if (CgwInputs.Count == 0) { CgwStatus = "Add the all-atom structures to map"; return; }
        var a = new JsonObject { ["inputs"] = Inputs(CgwInputs), ["o"] = _cgwMapOut };
        if (CgwIsCustom)
        {
            if (_cgwCut.Length == 0) { CgwStatus = "Custom rules: give the bonds to cut (SMARTS)"; return; }
            a["cut"] = _cgwCut;
            if (_cgwNames.Length > 0) a["names"] = _cgwNames;
        }
        else a["preset"] = _cgwPreset;
        if (_cgwMapDump.Length > 0) a["dump"] = _cgwMapDump;
        var r = await Cgw("cgmap", a);
        if (r == null) return;
        CgwMapRows.Clear();
        CgwKindBars.Clear();
        var sys = (r["systems"] as JsonArray ?? []).OfType<JsonObject>().ToList();
        var kinds = sys.SelectMany(s => (s["kinds"] as JsonObject ?? []).Select(k => k.Key)).Distinct().OrderBy(k => k, StringComparer.Ordinal).ToArray();
        CgwKindCats = kinds;
        var notes = new List<string>();
        CgwSystems.Clear();
        foreach (var s in sys)
        {
            var fr = string.Join(" · ", (s["fractions"] as JsonObject ?? []).Select(k => $"{k.Key} {CF(100 * CNum(k.Value), "0.0")} %"));
            CgwMapRows.Add(new DftTableRow([Path.GetFileName((string?)s["input"] ?? ""), CF(CNum(s["atoms"]), "0"), CF(CNum(s["beads"]), "0"), CF(CNum(s["chains"]), "0"),
                CTxt(s["beads_per_chain"]), s["beads_per_unit"] is JsonNode bu ? CF(CNum(bu), "0.###") : "–", fr, CNum(s["mass_error"]).ToString("0.0e+0", CgInv)]));
            var kd = s["kinds"] as JsonObject ?? [];
            CgwKindBars.Add((Path.GetFileNameWithoutExtension((string?)s["input"] ?? ""), kinds.Select((k, i) => ((double)i, kd[k] is JsonNode n ? CNum(n, 0) : 0)).ToArray()));
            foreach (var n in (s["notes"] as JsonArray ?? [])) notes.Add($"{Path.GetFileName((string?)s["input"] ?? "")}: {(string?)n}");
            var frames = s["trajectory"]?["out"] is JsonNode tr ? new[] { (string?)tr ?? "" } : new[] { (string?)s["data"] ?? "" };
            CgwSystems.Add(new CgwSysRow((string?)s["map"] ?? "", frames));
        }
        CgwSequence = sys.Count > 0 && (sys[0]["sequences"] as JsonArray)?.FirstOrDefault() is JsonNode q ? (string?)q ?? "" : "";
        CgwMapNotes = string.Join("\n", notes);
        var typesFile = Path.Combine(_cgwMapOut, "types.json");
        if (File.Exists(typesFile)) CgwTypes = typesFile;
        CgwBuildMaps = string.Join(",", CgwSystems.Select(x => x.Map));
        if (CgwBackRefMap.Length == 0 && CgwSystems.Count > 0) CgwBackRefMap = CgwSystems[0].Map;
        if (CgwBackRefData.Length == 0 && CgwInputs.Count > 0 && CgwInputs[0].EndsWith(".data", StringComparison.OrdinalIgnoreCase)) CgwBackRefData = CgwInputs[0];
        CgwStages[0].Done = true;
        CgwStatus = $"Mapped {sys.Count} system(s): {string.Join(", ", CgwMapRows.Select(x => $"{x.Cells[0]} {x.Cells[2]} beads"))}";
        CgwChanged?.Invoke();
    }

    // ---------------------------------------------------------------- 2 Bonded
    /// <summary>The all-atom reference systems (map + mapped frames): bonded inversion, g(r) targets, the melt's reference.</summary>
    public ObservableCollection<CgwSysRow> CgwSystems { get; } = new();
    private decimal _cgwT = 300;
    private string _cgwBondedOut = "", _cgwTimestep = "", _cgwNotConverged = "", _cgwBondedNotes = "";
    public decimal CgwT { get => _cgwT; set => Set(ref _cgwT, Math.Clamp(value, 1, 5000)); }
    public string CgwBondedOut { get => _cgwBondedOut; set => Set(ref _cgwBondedOut, value ?? ""); }
    public ObservableCollection<CgwBondedRow> CgwBondedRows { get; } = new();
    private CgwBondedRow? _cgwBondedSel;
    public CgwBondedRow? CgwBondedSel { get => _cgwBondedSel; set { if (Set(ref _cgwBondedSel, value)) { CgwLoadBondedPlot(); CgwChanged?.Invoke(); } } }
    public string CgwTimestep { get => _cgwTimestep; private set => Set(ref _cgwTimestep, value); }
    /// <summary>The NOT CONVERGED note of the inversion (the two halves of the samples disagree by more than 1 kT).</summary>
    public string CgwNotConverged { get => _cgwNotConverged; private set { if (Set(ref _cgwNotConverged, value)) Raise(nameof(HasCgwNotConverged)); } }
    public bool HasCgwNotConverged => _cgwNotConverged.Length > 0;
    public string CgwBondedNotes { get => _cgwBondedNotes; private set => Set(ref _cgwBondedNotes, value); }
    private JsonObject? _cgwBondedJson;
    public (double X, double Y)[] CgwPlotP { get; private set; } = [];
    public (double X, double Y)[] CgwPlotU { get; private set; } = [];
    public string CgwPlotUnit { get; private set; } = "";
    /// <summary>The sampled range (P above 0.1 % of its peak, padded) and the potential up to 12 kT above its minimum: the walls reach far higher.</summary>
    public (double Lo, double Hi, double UMax)? CgwPlotRange { get; private set; }

    public void CgwAddSystem(string map, IEnumerable<string> frames) => CgwSystems.Add(new CgwSysRow(map, frames.ToArray()));

    public async Task CgwRunBonded()
    {
        if (CgwSystems.Count == 0) { CgwStatus = "Add a system: its map.json and its mapped frames"; return; }
        var a = new JsonObject { ["inputs"] = Inputs(new[] { "bonded" }.Concat(SysInputs(CgwSystems))), ["T"] = (double)_cgwT, ["o"] = _cgwBondedOut };
        if (_cgwTypes.Length > 0) a["types"] = _cgwTypes;
        var r = await Cgw("cgfit", a);
        if (r == null) return;
        CgwBondedRows.Clear();
        var unsampled = 0;
        foreach (var t in (r["tables"] as JsonArray ?? []).OfType<JsonObject>())
        {
            if ((bool?)t["sampled"] != true) { unsampled++; continue; }
            var kind = (string?)t["kind"] ?? "";
            var u = kind == "bond" ? " Å" : " °";
            var hd = CNum(t["half_diff_kT"]);
            CgwBondedRows.Add(new CgwBondedRow(kind, (string?)t["key"] ?? "", CF(CNum(t["x0"]), "0.###") + u, $"{CF(CNum(t["lo"]), "0.##")}–{CF(CNum(t["hi"]), "0.##")}{u}",
                CF(hd, "0.00"), CF(CNum(t["count"]), "0"), hd > 1));
        }
        CgwTimestep = r["timestep_fs"] is JsonNode ts ? $"Suggested time step: {CF(CNum(ts), "0.0")} fs" : "";
        _cgwBondedJson = null;
        try { _cgwBondedJson = JsonNode.Parse(File.ReadAllText(Path.Combine(_cgwBondedOut, "bonded.json"))) as JsonObject; } catch (Exception e) { CgwStatus = "bonded.json: " + e.Message; }
        var notes = (_cgwBondedJson?["notes"] as JsonArray ?? []).Select(n => (string?)n ?? "").ToList();
        CgwNotConverged = string.Join("\n", notes.Where(n => n.Contains("NOT CONVERGED")));
        CgwBondedNotes = (unsampled > 0 ? $"{unsampled} terms of the type list had no samples (no table)" : "")
                         + string.Concat(notes.Where(n => !n.Contains("NOT CONVERGED") && !n.StartsWith("no samples")).Select(n => "\n" + n));
        CgwBondedSel = CgwBondedRows.FirstOrDefault(x => x.Flag) ?? CgwBondedRows.FirstOrDefault();
        CgwLoadBondedPlot();
        CgwStages[1].Done = true;
        CgwStatus = $"Bonded: {CgwBondedRows.Count} tables · {CgwBondedRows.Count(x => x.Flag)} with halves > 1 kT apart";
        CgwChanged?.Invoke();
    }

    private void CgwLoadBondedPlot()
    {
        CgwPlotP = []; CgwPlotU = [];
        if (_cgwBondedSel == null || _cgwBondedJson == null) return;
        var group = _cgwBondedSel.Kind == "bond" ? "bonds" : _cgwBondedSel.Kind == "angle" ? "angles" : "dihedrals";
        var t = (_cgwBondedJson[group] as JsonArray ?? []).OfType<JsonObject>().FirstOrDefault(x => (string?)x["key"] == _cgwBondedSel.Key);
        if (t == null) return;
        static double[] Arr(JsonNode? n) => (n as JsonArray ?? []).Select(v => CNum(v)).ToArray();
        var xh = Arr(t["xh"]); var p = Arr(t["P"]); var x = Arr(t["x"]); var u = Arr(t["U"]);
        CgwPlotP = xh.Zip(p).Select(z => (z.First, z.Second)).ToArray();
        CgwPlotU = x.Zip(u).Where(z => double.IsFinite(z.Second)).Select(z => (z.First, z.Second)).ToArray();
        CgwPlotUnit = _cgwBondedSel.Kind == "bond" ? "r (Å)" : _cgwBondedSel.Kind == "angle" ? "θ (°)" : "φ (°)";
        CgwPlotRange = null;
        if (CgwPlotP.Length > 0 && CgwPlotU.Length > 0)
        {
            var pmax = CgwPlotP.Max(q => q.Y);
            var on = CgwPlotP.Where(q => q.Y > 1e-3 * pmax).ToArray();
            var pad = _cgwBondedSel.Kind == "bond" ? 1.0 : 15.0;
            var lo = on.Min(q => q.X) - pad; var hi = on.Max(q => q.X) + pad;
            if (_cgwBondedSel.Kind != "bond") { lo = Math.Max(lo, CgwPlotU.Min(q => q.X)); hi = Math.Min(hi, CgwPlotU.Max(q => q.X)); }
            var umin = CgwPlotU.Where(q => q.X >= lo && q.X <= hi).Select(q => q.Y).DefaultIfEmpty(0).Min();
            CgwPlotRange = (lo, hi, umin + 12 * 0.0019872 * (double)_cgwT);
        }
    }

    // ---------------------------------------------------------------- 3 Non-bonded
    private string _cgwIbiOut = "", _cgwPressure = "1", _cgwTargets = "", _cgwPairsPath = "", _cgwIbiScript = "", _cgwIbiLoop = "", _cgwStepOut = "", _cgwStepSummary = "";
    private decimal _cgwRmax = 15, _cgwDr = 0.05m, _cgwRc = 15;
    public string CgwIbiOut { get => _cgwIbiOut; set => Set(ref _cgwIbiOut, value ?? ""); }
    public decimal CgwRmax { get => _cgwRmax; set => Set(ref _cgwRmax, Math.Clamp(value, 1, 100)); }
    public decimal CgwDr { get => _cgwDr; set => Set(ref _cgwDr, Math.Clamp(value, 0.001m, 1)); }
    public decimal CgwRc { get => _cgwRc; set => Set(ref _cgwRc, Math.Clamp(value, 1, 100)); }
    /// <summary>Each system's all-atom pressure (atm, comma-separated): the IBI pressure target.</summary>
    public string CgwPressure { get => _cgwPressure; set => Set(ref _cgwPressure, value ?? ""); }
    public string CgwTargets { get => _cgwTargets; set => Set(ref _cgwTargets, value ?? ""); }
    /// <summary>The current tables (pairs.json): ibi-step and fit start from them.</summary>
    public string CgwPairsPath { get => _cgwPairsPath; set { if (Set(ref _cgwPairsPath, value ?? "")) { CgwLoadPairs(); CgwChanged?.Invoke(); } } }
    public string CgwIbiScript { get => _cgwIbiScript; private set { if (Set(ref _cgwIbiScript, value)) Raise(nameof(HasCgwIbiScript)); } }
    public bool HasCgwIbiScript => _cgwIbiScript.Length > 0;
    /// <summary>The loop's command line (from the ibi-start report): it runs on a cluster.</summary>
    public string CgwIbiLoop { get => _cgwIbiLoop; private set => Set(ref _cgwIbiLoop, value); }
    public ObservableCollection<CgwSysRow> CgwStepRows { get; } = new();
    public string CgwStepOut { get => _cgwStepOut; set => Set(ref _cgwStepOut, value ?? ""); }
    public ObservableCollection<DftTableRow> CgwStepPairs { get; } = new();
    public string CgwStepSummary { get => _cgwStepSummary; private set => Set(ref _cgwStepSummary, value); }
    public ObservableCollection<string> CgwPairKeys { get; } = new();
    private string _cgwPairKey = "";
    public string CgwPairKey { get => _cgwPairKey; set { if (Set(ref _cgwPairKey, value ?? "")) CgwChanged?.Invoke(); } }
    /// <summary>g(r) per system (gr.csv): system → pair key → points.</summary>
    public Dictionary<string, Dictionary<string, (double X, double Y)[]>> CgwGr { get; } = new();
    private Dictionary<string, (double X, double Y)[]> _cgwPairU = new(), _cgwFitU = new();
    public (double X, double Y)[] CgwPairU => _cgwPairU.TryGetValue(_cgwPairKey, out var p) ? p : [];
    public (double X, double Y)[] CgwFitU => _cgwFitU.TryGetValue(_cgwPairKey, out var p) ? p : [];

    public async Task CgwRunTargets()
    {
        if (CgwSystems.Count == 0) { CgwStatus = "Add the systems (Bonded stage): maps and their frames"; return; }
        var a = new JsonObject { ["inputs"] = Inputs(new[] { "targets" }.Concat(SysInputs(CgwSystems))), ["T"] = (double)_cgwT, ["rmax"] = (double)_cgwRmax, ["dr"] = (double)_cgwDr, ["o"] = _cgwIbiOut };
        if (_cgwPressure.Length > 0) a["pressure"] = _cgwPressure;
        if (_cgwTypes.Length > 0) a["types"] = _cgwTypes;
        var r = await Cgw("cgfit", a);
        if (r == null) return;
        var files = (r["files"] as JsonArray ?? []).Select(f => (string?)f ?? "").ToList();
        CgwTargets = files.FirstOrDefault(f => f.EndsWith("targets.json")) ?? Path.Combine(_cgwIbiOut, "targets.json");
        CgwGr.Clear();
        foreach (var f in files.Where(f => f.EndsWith(".gr.csv")))
        {
            var (head, rows) = ReadTable(f);
            var d = new Dictionary<string, (double X, double Y)[]>();
            for (var c = 1; c < head.Length; ++c)
                if (rows.Any(x => c < x.Length && x[c] != 0)) d[head[c]] = rows.Where(x => c < x.Length).Select(x => (x[0], x[c])).ToArray();
            CgwGr[Path.GetFileName(f)[..^".gr.csv".Length]] = d;
        }
        CgwRefreshPairKeys();
        CgwStatus = $"Targets: {CgwGr.Count} system(s) · {CgwPairKeys.Count} pairs with g(r)";
        CgwChanged?.Invoke();
    }
    private void CgwRefreshPairKeys()
    {
        var keys = CgwGr.Values.SelectMany(d => d.Keys).Concat(_cgwPairU.Keys).Distinct().OrderBy(k => k, StringComparer.Ordinal).ToList();
        CgwPairKeys.Clear();
        foreach (var k in keys) CgwPairKeys.Add(k);
        if (!keys.Contains(_cgwPairKey)) CgwPairKey = keys.FirstOrDefault() ?? "";
        else Raise(nameof(CgwPairKey));
    }
    private static Dictionary<string, (double X, double Y)[]> ReadPairs(string path)
    {
        var d = new Dictionary<string, (double X, double Y)[]>();
        if (!File.Exists(path)) return d;
        try
        {
            var j = JsonNode.Parse(File.ReadAllText(path)) as JsonObject;
            var r0 = CNum(j?["r0"], 0); var dr = CNum(j?["dr"], 0.05);
            foreach (var p in (j?["pairs"] as JsonArray ?? []).OfType<JsonObject>())
                d[(string?)p["key"] ?? ""] = (p["U"] as JsonArray ?? []).Select((u, i) => (r0 + i * dr, CNum(u))).Where(x => double.IsFinite(x.Item2)).ToArray();
        }
        catch { }
        return d;
    }
    private void CgwLoadPairs()
    {
        _cgwPairU = ReadPairs(_cgwPairsPath);
        CgwRefreshPairKeys();
    }

    public async Task CgwIbiStart()
    {
        if (_cgwTargets.Length == 0) { CgwStatus = "Make the g(r) targets first"; return; }
        var a = new JsonObject { ["inputs"] = Inputs(["ibi-start"]), ["targets"] = _cgwTargets, ["bonded"] = _cgwBondedOut, ["rc"] = (double)_cgwRc, ["o"] = _cgwIbiOut };
        if (_cgwTypes.Length > 0) a["types"] = _cgwTypes;
        var r = await Cgw("cgfit", a);
        if (r == null) return;
        var files = (r["files"] as JsonArray ?? []).Select(f => (string?)f ?? "").ToList();
        CgwIbiScript = files.FirstOrDefault(f => f.EndsWith("run_ibi.sh")) ?? "";
        var loopLine = ((string?)r["text"] ?? "").Split('\n').FirstOrDefault(l => l.Contains("loop: "));
        CgwIbiLoop = loopLine != null ? loopLine[(loopLine.IndexOf("loop: ", StringComparison.Ordinal) + 6)..].Trim() : "";
        var p0 = files.FirstOrDefault(f => f.EndsWith("pairs.json")) ?? Path.Combine(_cgwIbiOut, "it000", "pairs.json");
        CgwPairsPath = p0;
        CgwStepOut = CgwNextIteration(p0);
        CgwStatus = $"IBI start: it000 tables for {_cgwPairU.Count} pairs · the loop runs on a cluster";
        CgwChanged?.Invoke();
    }
    /// <summary>itNNN beside the current tables' folder (it000 → it001).</summary>
    private string CgwNextIteration(string pairs)
    {
        var d = Path.GetDirectoryName(pairs) ?? _cgwIbiOut;
        var name = Path.GetFileName(d);
        if (name.StartsWith("it") && int.TryParse(name[2..], out var n)) return Path.Combine(Path.GetDirectoryName(d) ?? _cgwIbiOut, $"it{n + 1:000}");
        return Path.Combine(_cgwIbiOut, "it001");
    }

    public async Task CgwIbiStep()
    {
        if (CgwStepRows.Count == 0) { CgwStatus = "Add each system's map, CG dump and log from the last iteration"; return; }
        if (_cgwPairsPath.Length == 0 || _cgwTargets.Length == 0) { CgwStatus = "The step needs the targets and the current pairs.json"; return; }
        var a = new JsonObject { ["inputs"] = Inputs(new[] { "ibi-step" }.Concat(SysInputs(CgwStepRows))), ["targets"] = _cgwTargets, ["pairs"] = _cgwPairsPath, ["o"] = _cgwStepOut };
        if (_cgwTypes.Length > 0) a["types"] = _cgwTypes;
        var r = await Cgw("cgfit", a);
        if (r == null) return;
        CgwStepPairs.Clear();
        foreach (var p in (r["pairs"] as JsonArray ?? []).OfType<JsonObject>())
            CgwStepPairs.Add(new DftTableRow([(string?)p["key"] ?? "", CF(CNum(p["residual"]), "0.0000"), CF(CNum(p["max_dev"]), "0.000")]));
        CgwStepSummary = $"residual {CF(CNum(r["residual"]), "0.0000")}" + (r["pressure"] is JsonNode pr ? " · pressure " + pr.ToJsonString() : "");
        var next = Path.Combine(_cgwStepOut, "pairs.json");
        if (File.Exists(next)) { CgwPairsPath = next; CgwStepOut = CgwNextIteration(next); }
        CgwChanged?.Invoke();
    }

    public static readonly string[] CgwForms = ["lj126", "lj96", "morse", "mie"];
    private string _cgwForm = "lj126", _cgwFits = "", _cgwTgFile = "", _cgwCalNext = "";
    public string CgwForm { get => _cgwForm; set => Set(ref _cgwForm, value ?? "lj126"); }
    public string CgwFits { get => _cgwFits; set => Set(ref _cgwFits, value ?? ""); }
    public async Task CgwFit()
    {
        if (_cgwPairsPath.Length == 0) { CgwStatus = "Fit needs the tables (pairs.json of the last IBI iteration)"; return; }
        var o = CgwSub("fit_" + _cgwForm);
        var a = new JsonObject { ["inputs"] = Inputs(["fit"]), ["pairs"] = _cgwPairsPath, ["form"] = _cgwForm, ["o"] = o };
        if (_cgwTypes.Length > 0) a["types"] = _cgwTypes;
        var r = await Cgw("cgfit", a);
        if (r == null) return;
        CgwFits = Path.Combine(o, "fits.json");
        _cgwFitU = ReadPairs(Path.Combine(o, "tables", "pairs.json"));
        CgwStages[2].Done = true;
        CgwStatus = $"Fit ({_cgwForm}): {_cgwFitU.Count} pairs · {CgwFits}";
        CgwChanged?.Invoke();
    }

    private string _cgwHistory = "";
    private decimal _cgwSSigma = 1, _cgwSEps = 1, _cgwDensity, _cgwTg, _cgwTargetDensity, _cgwTargetTg;
    public string CgwHistory { get => _cgwHistory; set => Set(ref _cgwHistory, value ?? ""); }
    public decimal CgwSSigma { get => _cgwSSigma; set => Set(ref _cgwSSigma, Math.Clamp(value, 0.1m, 10)); }
    public decimal CgwSEps { get => _cgwSEps; set => Set(ref _cgwSEps, Math.Clamp(value, 0.01m, 100)); }
    public decimal CgwDensity { get => _cgwDensity; set => Set(ref _cgwDensity, Math.Clamp(value, 0, 10)); }
    public decimal CgwTg { get => _cgwTg; set => Set(ref _cgwTg, Math.Clamp(value, 0, 3000)); }
    public decimal CgwTargetDensity { get => _cgwTargetDensity; set => Set(ref _cgwTargetDensity, Math.Clamp(value, 0, 10)); }
    public decimal CgwTargetTg { get => _cgwTargetTg; set => Set(ref _cgwTargetTg, Math.Clamp(value, 0, 3000)); }
    public string CgwTgFile { get => _cgwTgFile; set => Set(ref _cgwTgFile, value ?? ""); }
    public string CgwCalNext { get => _cgwCalNext; private set => Set(ref _cgwCalNext, value); }

    /// <summary>T_g from a cooling scan's "T density" lines; the value goes into the calibration's reported T_g.</summary>
    public async Task CgwRunTg()
    {
        if (!File.Exists(_cgwTgFile)) { CgwStatus = "Pick the cooling scan (T density per line)"; return; }
        var r = await Cgw("cgfit", new JsonObject { ["inputs"] = Inputs(["tg", _cgwTgFile]) });
        if (r == null) return;
        var tg = CNum(r["tg"]);
        if (double.IsFinite(tg)) CgwTg = (decimal)Math.Round(tg, 1);
        CgwStatus = $"T_g {CF(tg, "0.0")} ± {CF(CNum(r["tg_err"]), "0.0")} K";
    }
    public async Task CgwCalibrate()
    {
        if (_cgwFits.Length == 0) { CgwStatus = "Calibrate scales a fit: run Fit first (or pick fits.json)"; return; }
        var step = 1;
        var calib = Path.GetDirectoryName(_cgwHistory) ?? CgwSub("calib");
        while (Directory.Exists(Path.Combine(calib, $"step{step}"))) step++;
        var a = new JsonObject
        {
            ["inputs"] = Inputs(["calibrate"]), ["fits"] = _cgwFits, ["history"] = _cgwHistory, ["s_sigma"] = (double)_cgwSSigma, ["s_eps"] = (double)_cgwSEps,
            ["o"] = Path.Combine(calib, $"step{step}"),
        };
        if (_cgwDensity > 0) a["density"] = (double)_cgwDensity;
        if (_cgwTg > 0) a["tg"] = (double)_cgwTg;
        if (_cgwTargetDensity > 0) a["target_density"] = (double)_cgwTargetDensity;
        if (_cgwTargetTg > 0) a["target_tg"] = (double)_cgwTargetTg;
        if (_cgwTypes.Length > 0) a["types"] = _cgwTypes;
        Directory.CreateDirectory(calib);
        var r = await Cgw("cgfit", a);
        if (r == null) return;
        CgwCalNext = $"next: s_σ {CF(CNum(r["s_sigma"]), "0.0000")} · s_ε {CF(CNum(r["s_eps"]), "0.0000")}  (step{step}: pair.in and decks for them)";
        CgwStatus = CgwCalNext;
    }

    // ---------------------------------------------------------------- 4 Build CG melt
    public static readonly string[] CgwSequences = ["bernoulli", "markov", "block", "gradient", "alternating", "pattern"];
    public static readonly string[] CgwLengthKinds = ["monodisperse", "schulz-zimm", "flory", "poisson", "log-normal"];
    private string _cgwUnits = "BS=B+S", _cgwComposition = "", _cgwSeq = "bernoulli", _cgwSeqParam = "", _cgwLengths = "monodisperse", _cgwBuildMaps = "", _cgwBuildOut = "", _cgwBuildWarning = "";
    private int _cgwDp = 100, _cgwChains = 10;
    private decimal _cgwPdi = 1, _cgwBuildDensity = 1.2m;
    private bool _cgwBuildRef;
    public string CgwUnits { get => _cgwUnits; set => Set(ref _cgwUnits, value ?? ""); }
    public string CgwComposition { get => _cgwComposition; set => Set(ref _cgwComposition, value ?? ""); }
    public string CgwSeq { get => _cgwSeq; set { if (Set(ref _cgwSeq, value ?? "bernoulli")) { Raise(nameof(CgwSeqHasParam)); Raise(nameof(CgwSeqParamLabel)); } } }
    public bool CgwSeqHasParam => _cgwSeq is "markov" or "block" or "pattern";
    public string CgwSeqParamLabel => _cgwSeq switch { "markov" => "Transitions FROM>TO:p", "block" => "Block lengths (cycling the units)", "pattern" => "Pattern (A, B … for the units)", _ => "" };
    public string CgwSeqParam { get => _cgwSeqParam; set => Set(ref _cgwSeqParam, value ?? ""); }
    public int CgwDp { get => _cgwDp; set => Set(ref _cgwDp, Math.Clamp(value, 1, 100000)); }
    public int CgwChains { get => _cgwChains; set => Set(ref _cgwChains, Math.Clamp(value, 1, 1000000)); }
    public string CgwLengths { get => _cgwLengths; set { if (Set(ref _cgwLengths, value ?? "monodisperse")) Raise(nameof(CgwDrawn)); } }
    public bool CgwDrawn => _cgwLengths != "monodisperse";
    public decimal CgwPdi { get => _cgwPdi; set => Set(ref _cgwPdi, Math.Clamp(value, 1, 20)); }
    public decimal CgwBuildDensity { get => _cgwBuildDensity; set => Set(ref _cgwBuildDensity, Math.Clamp(value, 0.01m, 10)); }
    public string CgwBuildMaps { get => _cgwBuildMaps; set => Set(ref _cgwBuildMaps, value ?? ""); }
    public string CgwBuildOut { get => _cgwBuildOut; set => Set(ref _cgwBuildOut, value ?? ""); }
    /// <summary>Compare ⟨R²(n)⟩/n with the all-atom reference systems (the Bonded stage's list).</summary>
    public bool CgwBuildRef { get => _cgwBuildRef; set => Set(ref _cgwBuildRef, value); }
    public ObservableCollection<Row> CgwBuildFacts { get; } = new();
    public string CgwBuildWarning { get => _cgwBuildWarning; private set { if (Set(ref _cgwBuildWarning, value)) Raise(nameof(HasCgwBuildWarning)); } }
    public bool HasCgwBuildWarning => _cgwBuildWarning.Length > 0;
    public (double X, double Y)[] CgwR2Built { get; private set; } = [];
    public (double X, double Y)[] CgwR2Aa { get; private set; } = [];

    public async Task CgwBuild()
    {
        if (_cgwUnits.Length == 0) { CgwStatus = "Give the repeat units as NAME=BEAD+BEAD"; return; }
        var bonded = _cgwBondedOut;
        if (!File.Exists(Path.Combine(bonded, "bonded.json"))) { CgwStatus = "The melt draws its chains from bonded.json: run the Bonded stage (or set its folder)"; return; }
        var a = new JsonObject
        {
            ["units"] = _cgwUnits, ["sequence"] = _cgwSeq, ["dp"] = _cgwDp, ["chains"] = _cgwChains, ["lengths"] = _cgwLengths, ["density"] = (double)_cgwBuildDensity,
            ["bonded"] = bonded, ["o"] = _cgwBuildOut,
        };
        if (_cgwComposition.Length > 0) a["composition"] = _cgwComposition;
        if (CgwSeqHasParam && _cgwSeqParam.Length > 0) a[_cgwSeq == "markov" ? "markov" : _cgwSeq == "block" ? "blocks" : "pattern"] = _cgwSeqParam;
        if (CgwDrawn) a["pdi"] = (double)_cgwPdi;
        if (_cgwBuildMaps.Length > 0) a["maps"] = _cgwBuildMaps;
        if (_cgwTypes.Length > 0) a["types"] = _cgwTypes;
        if (_cgwBuildRef && CgwSystems.Count > 0) a["inputs"] = Inputs(SysInputs(CgwSystems));
        var r = await Cgw("cgbuild", a);
        if (r == null) return;
        CgwBuildFacts.Clear();
        var box = CNum(r["box"]); var ree = CNum(r["ree_rms"]);
        CgwBuildFacts.Add(new Row("Beads", CF(CNum(r["beads"]), "0")));
        CgwBuildFacts.Add(new Row("Box edge L", CF(box, "0.0") + " Å"));
        CgwBuildFacts.Add(new Row("R_ee (rms)", CF(ree, "0.0") + " Å"));
        CgwBuildFacts.Add(new Row("L / R_ee", CF(box / ree, "0.00")));
        CgwBuildFacts.Add(new Row("Mean contour length", CF(CNum(r["contour"]), "0.0") + " Å"));
        CgwBuildWarning = string.Join("\n", ((string?)r["text"] ?? "").Split('\n').Where(l => l.Contains("WARNING")).Select(l => l.Trim().Replace("note: ", "")));
        var (head, rows) = ReadTable(Path.Combine(_cgwBuildOut, "melt.internal.csv"));
        CgwR2Built = rows.Where(x => x.Length > 1).Select(x => (x[0], x[1])).ToArray();
        CgwR2Aa = head.Length > 2 ? rows.Where(x => x.Length > 2 && x[2] > 0).Select(x => (x[0], x[2])).ToArray() : [];
        var map = Path.Combine(_cgwBuildOut, "melt.map.json");
        var data = Path.Combine(_cgwBuildOut, "melt.cg.data");
        CgwPpaRows.Clear();
        if (File.Exists(map)) CgwPpaRows.Add(new CgwSysRow(map, [data]));
        if (File.Exists(map)) { CgwBackCg = map; CgwDynMap = map; CgwMechMap = map; }
        CgwLoadEquil();
        CgwStages[3].Done = true;
        CgwStatus = $"Melt: {CF(CNum(r["beads"]), "0")} beads in a {CF(box, "0.0")} Å box";
        CgwChanged?.Invoke();
    }

    // ---------------------------------------------------------------- 5 Equilibrate (the deck, run elsewhere)
    private string _cgwEquilDeck = "", _cgwEquilCommand = "", _cgwEquilPath = "";
    public string CgwEquilDeck { get => _cgwEquilDeck; private set { if (Set(ref _cgwEquilDeck, value)) Raise(nameof(HasCgwEquil)); } }
    public bool HasCgwEquil => _cgwEquilDeck.Length > 0;
    public string CgwEquilCommand { get => _cgwEquilCommand; private set => Set(ref _cgwEquilCommand, value); }
    public string CgwEquilPath { get => _cgwEquilPath; private set => Set(ref _cgwEquilPath, value); }
    public void CgwLoadEquil()
    {
        var p = Path.Combine(_cgwBuildOut, "in.cg_equil");
        CgwEquilPath = File.Exists(p) ? p : "";
        CgwEquilDeck = File.Exists(p) ? File.ReadAllText(p) : "";
        // the deck says how it is run in its own header (a comment line "lmp -in in.cg_equil …")
        var line = _cgwEquilDeck.Split('\n').Take(10).FirstOrDefault(l => l.TrimStart('#', ' ').StartsWith("lmp ", StringComparison.Ordinal));
        CgwEquilCommand = line == null ? "" : $"cd {_cgwBuildOut} && {line.TrimStart('#', ' ').Trim()}";
        if (_cgwEquilDeck.Length > 0) CgwStages[4].Done = true;
    }

    // ---------------------------------------------------------------- 6 Analyse
    public ObservableCollection<CgwSysRow> CgwPpaRows { get; } = new();
    public static readonly string[] CgwPpaMethods = ["caps", "z1", "both", "lammps"];
    private string _cgwPpaMethod = "caps", _cgwPpaMulti = "";
    public string CgwPpaMethod { get => _cgwPpaMethod; set => Set(ref _cgwPpaMethod, value ?? "caps"); }
    public ObservableCollection<DftTableRow> CgwPpaTable { get; } = new();
    public string CgwPpaMulti { get => _cgwPpaMulti; private set => Set(ref _cgwPpaMulti, value); }
    public async Task CgwRunPpa()
    {
        if (CgwPpaRows.Count == 0) { CgwStatus = "Add the CG systems: map.json and frames (an equilibrated melt)"; return; }
        var r = await Cgw("ppa", new JsonObject { ["inputs"] = Inputs(SysInputs(CgwPpaRows)), ["method"] = _cgwPpaMethod, ["o"] = CgwSub("ppa") });
        if (r == null) return;
        CgwPpaTable.Clear();
        foreach (var s in (r["systems"] as JsonArray ?? []).OfType<JsonObject>())
        {
            var n = CNum(s["N"]); var ne = CNum(s["ne_mod_s_coil"]); var m = CNum(s["bead_mass"]); var z = CNum(s["z"], -1);
            CgwPpaTable.Add(new DftTableRow([(string?)s["system"] ?? "", CF(n, "0"), CF(CNum(s["chains"]), "0"), CF(CNum(s["lpp"]), "0.0"), CF(CNum(s["a_pp"]), "0.0"),
                CF(CNum(s["ne_s_coil"]), "0.0"), CF(ne, "0.0"), z >= 0 ? CF(CNum(s["ne_s_kink"]), "0.0") : "–", z >= 0 ? CF(CNum(s["ne_mod_s_kink"]), "0.0") : "–",
                ne > 0 ? CF(ne * m, "0") : "–", ne > 0 ? CF(n / ne, "0.00") : "–"]));
        }
        CgwPpaMulti = r["ne_m_kink"] != null || r["ne_m_coil"] != null
            ? $"over the chain lengths: M-kink N_e {CF(CNum(r["ne_m_kink"]), "0.0")} · M-coil N_e {CF(CNum(r["ne_m_coil"]), "0.0")}" : "M-kink and M-coil need two chain lengths or more";
        CgwStages[5].Done = true;
        CgwStatus = $"Entanglements: {CgwPpaTable.Count} system(s)";
        CgwChanged?.Invoke();
    }

    public static readonly string[] CgwMechModes = ["both", "stress", "volume"];
    private string _cgwRates = "1e-6,1e-7", _cgwMechMode = "both", _cgwMechFiles = "", _cgwStressFile = "", _cgwMechMap = "", _cgwMechDump = "";
    private decimal _cgwMaxStrain = 3, _cgwMechT = 300, _cgwMechP = 1, _cgwMechDt = 10, _cgwProbe;
    private int _cgwMechFrames = 60;
    public string CgwRates { get => _cgwRates; set => Set(ref _cgwRates, value ?? ""); }
    public string CgwMechMode { get => _cgwMechMode; set => Set(ref _cgwMechMode, value ?? "both"); }
    public decimal CgwMaxStrain { get => _cgwMaxStrain; set => Set(ref _cgwMaxStrain, Math.Clamp(value, 0.01m, 20)); }
    public decimal CgwMechT { get => _cgwMechT; set => Set(ref _cgwMechT, Math.Clamp(value, 1, 3000)); }
    public decimal CgwMechP { get => _cgwMechP; set => Set(ref _cgwMechP, Math.Clamp(value, 0, 100000)); }
    public decimal CgwMechDt { get => _cgwMechDt; set => Set(ref _cgwMechDt, Math.Clamp(value, 0.1m, 100)); }
    public int CgwMechFrames { get => _cgwMechFrames; set => Set(ref _cgwMechFrames, Math.Clamp(value, 1, 10000)); }
    public string CgwMechFiles { get => _cgwMechFiles; private set => Set(ref _cgwMechFiles, value); }
    public string CgwStressFile { get => _cgwStressFile; set => Set(ref _cgwStressFile, value ?? ""); }
    public string CgwMechMap { get => _cgwMechMap; set => Set(ref _cgwMechMap, value ?? ""); }
    public string CgwMechDump { get => _cgwMechDump; set => Set(ref _cgwMechDump, value ?? ""); }
    public decimal CgwProbe { get => _cgwProbe; set => Set(ref _cgwProbe, Math.Clamp(value, 0, 10)); }
    public ObservableCollection<Row> CgwMechFacts { get; } = new();
    public List<(string Name, (double X, double Y)[] Pts)> CgwStress { get; } = new();

    public async Task CgwMechDecks()
    {
        var o = CgwSub("tension");
        var r = await Cgw("mech", new JsonObject
        {
            ["inputs"] = Inputs(["decks"]), ["rates"] = _cgwRates, ["mode"] = _cgwMechMode, ["max_strain"] = (double)_cgwMaxStrain, ["T"] = (double)_cgwMechT,
            ["P"] = (double)_cgwMechP, ["dt"] = (double)_cgwMechDt, ["frames"] = _cgwMechFrames, ["o"] = o,
        });
        if (r == null) return;
        var files = (r["files"] as JsonArray ?? []).Select(f => Path.GetFileName((string?)f ?? "")).ToList();
        CgwMechFiles = $"{o}: {string.Join(" · ", files)}";
        CgwStatus = $"Tension decks: {files.Count(f => f.StartsWith("in."))} runs (run_tension.sh on a cluster)";
    }
    public async Task CgwMechAnalyze()
    {
        if (!File.Exists(_cgwStressFile)) { CgwStatus = "Pick a STEM.stress_strain.dat written by a tension run"; return; }
        var ins = new List<string> { "analyze", _cgwStressFile };
        if (File.Exists(_cgwMechMap) && File.Exists(_cgwMechDump)) { ins.Add(_cgwMechMap); ins.Add(_cgwMechDump); }
        var a = new JsonObject { ["inputs"] = Inputs(ins), ["o"] = Path.Combine(CgwSub("tension"), "analysis") };
        if (_cgwProbe > 0) a["probe"] = (double)_cgwProbe;
        var r = await Cgw("mech", a);
        if (r == null) return;
        CgwMechFacts.Clear();
        CgwMechFacts.Add(new Row("Modulus", CF(CNum(r["modulus"]), "0.0") + " MPa"));
        CgwMechFacts.Add(new Row("Yield", $"{CF(CNum(r["yield_stress"]), "0.0")} MPa at strain {CF(CNum(r["yield_strain"]), "0.000")}"));
        CgwMechFacts.Add(new Row("Softening", CF(CNum(r["softening"]), "0.0") + " MPa"));
        CgwMechFacts.Add(new Row("Hardening G_R", CF(CNum(r["hardening_modulus"]), "0.00") + " MPa"));
        CgwStress.Clear();
        var (_, rows) = ReadTable(_cgwStressFile);
        string[] names = ["σ", "bond", "angle", "dihedral", "pair", "kinetic"];
        for (var c = 1; c <= 6; ++c)
        {
            var pts = rows.Where(x => x.Length > c && double.IsFinite(x[c])).Select(x => (x[0], x[c])).ToArray();
            if (pts.Length > 0) CgwStress.Add((names[c - 1], pts));
        }
        CgwStatus = "Tension analysed: " + Path.GetFileName(_cgwStressFile);
        CgwChanged?.Invoke();
    }

    private string _cgwDynMap = "", _cgwDynDump = "", _cgwAaCsv = "", _cgwCgCsv = "", _cgwTimeMap = "";
    private decimal _cgwDynDt = 10;
    public string CgwDynMap { get => _cgwDynMap; set => Set(ref _cgwDynMap, value ?? ""); }
    public string CgwDynDump { get => _cgwDynDump; set => Set(ref _cgwDynDump, value ?? ""); }
    /// <summary>fs per timestep of the dump.</summary>
    public decimal CgwDynDt { get => _cgwDynDt; set => Set(ref _cgwDynDt, Math.Clamp(value, 0.01m, 1000)); }
    public ObservableCollection<Row> CgwDynFacts { get; } = new();
    public List<(string Name, (double X, double Y)[] Pts)> CgwDyn { get; } = new();
    public string CgwAaCsv { get => _cgwAaCsv; set => Set(ref _cgwAaCsv, value ?? ""); }
    public string CgwCgCsv { get => _cgwCgCsv; set => Set(ref _cgwCgCsv, value ?? ""); }
    public string CgwTimeMap { get => _cgwTimeMap; private set => Set(ref _cgwTimeMap, value); }

    public async Task CgwRunDyn()
    {
        if (!File.Exists(_cgwDynMap) || !File.Exists(_cgwDynDump)) { CgwStatus = "Dynamics: pick the map.json and the CG dump"; return; }
        var o = CgwSub("dynamics");
        var r = await Cgw("cgdyn", new JsonObject { ["inputs"] = Inputs([_cgwDynMap, _cgwDynDump]), ["dt"] = (double)_cgwDynDt, ["o"] = o });
        if (r == null) return;
        CgwDynFacts.Clear();
        CgwDynFacts.Add(new Row("D", $"{CNum(r["D"]).ToString("0.###e+0", CgInv)} Å²/fs"));
        var tr = CNum(r["tau_R"]); var te = CNum(r["tau_e"]);
        CgwDynFacts.Add(new Row("τ_R", tr > 0 ? CF(tr / 1000, "0.###") + " ps" : "not reached"));
        CgwDynFacts.Add(new Row("τ_e", te > 0 ? CF(te / 1000, "0.###") + " ps" : "not reached"));
        CgwDyn.Clear();
        var csv = Path.Combine(o, "dynamics.csv");
        var (_, rows) = ReadTable(csv);
        foreach (var (name, c) in new[] { ("g₁", 1), ("g₂", 2), ("g₃", 3) })
            CgwDyn.Add((name, rows.Where(x => x.Length > c && x[0] > 0 && x[c] > 0).Select(x => (Math.Log10(x[0]), x[c])).ToArray()));
        if (_cgwCgCsv.Length == 0) CgwCgCsv = csv;
        CgwStatus = "Dynamics: " + csv;
        CgwChanged?.Invoke();
    }
    public async Task CgwRunTimeMap()
    {
        if (!File.Exists(_cgwAaCsv) || !File.Exists(_cgwCgCsv)) { CgwStatus = "Time mapping: the AA and the CG dynamics.csv"; return; }
        var r = await Cgw("cgdyn", new JsonObject { ["inputs"] = Inputs(["timemap", _cgwAaCsv, _cgwCgCsv]) });
        if (r == null) return;
        CgwTimeMap = $"t_AA = {CNum(r["factor"]).ToString("0.###", CgInv)} × t_CG · ln s spread ±{CF(CNum(r["spread"]), "0.00")}";
    }

    private string _cgwBackRefMap = "", _cgwBackRefData = "", _cgwBackCg = "", _cgwBackFrame = "", _cgwBackInput = "", _cgwBackCheckData = "";
    public string CgwBackRefMap { get => _cgwBackRefMap; set => Set(ref _cgwBackRefMap, value ?? ""); }
    public string CgwBackRefData { get => _cgwBackRefData; set => Set(ref _cgwBackRefData, value ?? ""); }
    public string CgwBackCg { get => _cgwBackCg; set => Set(ref _cgwBackCg, value ?? ""); }
    public string CgwBackFrame { get => _cgwBackFrame; set => Set(ref _cgwBackFrame, value ?? ""); }
    public string CgwBackInput { get => _cgwBackInput; set => Set(ref _cgwBackInput, value ?? ""); }
    public string CgwBackCheckData { get => _cgwBackCheckData; set => Set(ref _cgwBackCheckData, value ?? ""); }
    public ObservableCollection<Row> CgwBackFacts { get; } = new();
    public ObservableCollection<DftTableRow> CgwBackCheck { get; } = new();
    public async Task CgwBackmap()
    {
        if (!File.Exists(_cgwBackRefMap) || !File.Exists(_cgwBackRefData) || !File.Exists(_cgwBackCg) || !File.Exists(_cgwBackFrame))
        { CgwStatus = "Backmap: the reference map and all-atom data, the CG map and its frame"; return; }
        var o = CgwSub("backmap");
        var a = new JsonObject { ["inputs"] = Inputs([_cgwBackRefMap, _cgwBackRefData]), ["cg"] = _cgwBackCg, ["frame"] = _cgwBackFrame, ["o"] = o };
        if (File.Exists(_cgwBackInput)) a["input"] = _cgwBackInput;
        var r = await Cgw("backmap", a);
        if (r == null) return;
        CgwBackFacts.Clear();
        CgwBackFacts.Add(new Row("Atoms", CF(CNum(r["atoms"]), "0")));
        CgwBackFacts.Add(new Row("Net charge", CF(CNum(r["charge"]), "0.0000") + " e"));
        CgwBackFacts.Add(new Row("Restored cut bonds", "mean " + CF(CNum(r["cut_mean"]), "0.00") + " Å (before relaxation)"));
        CgwBackFacts.Add(new Row("Unmatched terms", CTxt(r["unmatched"])));
        var data = Path.Combine(o, "backmapped.data");
        if (File.Exists(data)) CgwBackCheckData = data;
        CgwStatus = "Backmapped: " + data;
    }
    public async Task CgwBackmapCheck()
    {
        if (!File.Exists(_cgwBackCheckData)) { CgwStatus = "Check: pick an all-atom data file"; return; }
        var r = await Cgw("backmap", new JsonObject { ["inputs"] = Inputs(["check", _cgwBackCheckData]) });
        if (r == null) return;
        CgwBackCheck.Clear();
        foreach (var x in (r["rows"] as JsonArray ?? []).OfType<JsonObject>())
        {
            var bond = (string?)x["kind"] == "bond";
            var u = bond ? " Å" : " °";
            CgwBackCheck.Add(new DftTableRow([(string?)x["kind"] ?? "", CTxt(x["type"]), CTxt(x["count"]), x["ref"] is JsonNode rf ? CF(CNum(rf), "0.000") + u : "not harmonic",
                CF(CNum(x["mean_dev"]), "0.000"), CF(CNum(x["max_dev"]), "0.000")]));
        }
        CgwStatus = $"Check: {CgwBackCheck.Count} term types";
    }
}
