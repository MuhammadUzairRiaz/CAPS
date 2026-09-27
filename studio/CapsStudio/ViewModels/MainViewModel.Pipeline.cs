using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>One step of the pipeline strip (design/boards/PipelineGrow). State: done | warn | optional | next | todo.</summary>
public sealed record PipelineStep(int Number, string Name, string Detail, string State, int Module, bool Current)
{
    public string Mark => Number.ToString(CultureInfo.InvariantCulture);
    public bool IsDone => State == "done" && !Current;
    public bool IsWarn => State == "warn" && !Current;
    public bool IsOptional => State == "optional" && !Current;
    public bool IsNext => State == "next" && !Current;
    public bool IsTodo => State == "todo" && !Current;
    public bool ShowNumber => Number > 0 && !IsDone && !IsWarn;
    public bool HasDetail => Detail.Length > 0;
}

public sealed record EngineFile(string Name, string What, string Size);

/// <summary>The pipeline in order (design/boards/PipelineGrow, ForceFieldStep, ExportCenter): a strip over every page
/// that says where the project is (Build, Grow, Force field, Pack, Relax, Equilibrate, Dynamics, Analyze), the force
/// field chosen in Grow and assigned when growing (or packing) finishes, and the Export center that writes LAMMPS and
/// GROMACS files from any step.</summary>
public sealed partial class MainViewModel
{
    // ---------------------------------------------------------------- where the project is
    private readonly HashSet<string> _pipeDone = new();
    private string _pipeBuild = "";
    public ObservableCollection<PipelineStep> PipelineSteps { get; } = new();
    public bool ShowPipelineStrip => _doc != null;

    /// <summary>A new document: what made it (Grow, Pack or a file) is marked by the caller after this.</summary>
    private void PipelineNewDocument(string title)
    {
        _pipeDone.Clear();
        _pipeBuild = title.Replace(" (unsaved)", "");
        RefreshSteps();
    }

    private void MarkPipeline(string step, string? buildDetail = null)
    {
        _pipeDone.Add(step);
        if (buildDetail != null) _pipeBuild = buildDetail;
        if (_activeItem != null && step is "Grow" or "Pack") _activeItem.Origin = step == "Grow" ? "Polymer cell" : "Packing";
        RefreshSteps();
    }

    public void RefreshSteps()
    {
        Raise(nameof(ShowPipelineStrip));
        PipelineSteps.Clear();
        if (_doc == null) return;
        var inv = CultureInfo.InvariantCulture;
        long atoms;
        try { atoms = _doc.Summary().Atoms; } catch (ObjectDisposedException) { return; }   // closed under us: no strip
        var ff = Field.Assigned ? Field.ForceFieldName : "";
        var ffState = !Field.Assigned ? "todo" : Field.Complete ? "done" : "warn";
        var ffDetail = !Field.Assigned ? "not assigned" : Field.Complete ? ShortFf(ff) + (_pipeAutoFf ? " · auto" : "") : "incomplete";
        // the active structure's state, no order imposed: each module works on whichever structure is active
        string Did(string step, string done) => _pipeDone.Contains(step) ? done : "not yet";
        var made = _activeItem?.Origin is { Length: > 0 } o ? o : "";
        var rows = new List<(string Name, string Detail, string State, int Module, bool Current)>
        {
            ("Structure", Shorten(Title.Replace(" (unsaved)", ""), 26), "done", -1, false),
            ("Made by", made.Length > 0 ? made : Shorten(_pipeBuild, 22), "done", made == "Polymer cell" ? 0 : made == "Packing" ? 5 : 13, false),
            ("Force field", ffDetail, ffState, 7, _module == 7),
            // optional in CAPS: a built cell can go straight to LAMMPS or GROMACS
            ("Minimise", Did("Relax", "minimised").Replace("not yet", "optional"), _pipeDone.Contains("Relax") ? "done" : "todo", 2, _module == 2),
            ("Equilibrate", Did("Equilibrate", "equilibrated").Replace("not yet", "optional"), _pipeDone.Contains("Equilibrate") ? "done" : "todo", 4, _module == 4),
            ("Dynamics", Did("Dynamics", "run").Replace("not yet", "optional"), _pipeDone.Contains("Dynamics") ? "done" : "todo", 3, _module == 3),
            ("Export", Did("Export", "written"), _pipeDone.Contains("Export") ? "done" : "todo", 68, _module == 68),
        };
        foreach (var r in rows) PipelineSteps.Add(new PipelineStep(0, r.Name, r.Detail, r.State, r.Module, r.Current));
        if (_activeItem != null) UpdateItemInfo(_activeItem);
        Raise(nameof(PipelineNextLabel));
        Raise(nameof(HasPipelineNext));
    }

    private static string ShortFf(string name) => Shorten(name, 22);
    private static string Shorten(string s, int n) => s.Length <= n ? s : s[..(n - 1)] + "…";

    // no "next step": modules are independent (kept for the command palette and older callers)
    public string PipelineNextLabel => "";
    public bool HasPipelineNext => PipelineNextLabel.Length > 0;
    public void GoPipelineNext()
    {
        if (PipelineSteps.FirstOrDefault(s => s.State is "next" or "warn") is { Current: false } s) GoPipelineStep(s);
    }
    public void GoPipelineStep(PipelineStep s)
    {
        if (s.Module < 0) return;
        if (s.Module == 13) { SetModule(13); return; }
        SetModule(s.Module);
    }

    /// <summary>The rail's Build and Export entries light up for their pages.</summary>
    public bool IsBuildRail => _module is 13 or 14 or 15 or 16 or 29 or 30 or 31 or 44 or 60;
    public bool IsExportRail => _module is 68 or 21 or 18 or 26 or 54;

    // ---------------------------------------------------------------- the force field in Grow and Pack
    private bool _growAssignFf = true, _pipeAutoFf;
    /// <summary>Assign Grow's own force field choice as soon as growing finishes (Pack has its own).</summary>
    public bool GrowAssignField { get => _growAssignFf; set => Set(ref _growAssignFf, value); }
    public string GrowFieldLine => !Field.Assigned ? (_growAssignFf ? "Assigned when growing finishes" : "Not assigned: choose one in the Force field step")
        : Field.Complete ? $"{Field.ForceFieldName}: every atom typed, every term found"
        : $"{Field.ForceFieldName}: {Field.UntypedText}, {Field.MissingText}";
    public bool GrowFieldOk => Field.Assigned && Field.Complete;
    public bool GrowFieldWarn => Field.Assigned && !Field.Complete;

    private void RaiseGrowField() { Raise(nameof(GrowFieldLine)); Raise(nameof(GrowFieldOk)); Raise(nameof(GrowFieldWarn)); }

    // ---------------------------------------------------------------- Export center (module 68)
    public bool IsExportCenter => _module == 68;
    public static readonly string[] EngineRuns = ["Check (single point)", "Minimise", "NVT", "NPT"];
    private static readonly string[] EngineRunIds = ["check", "minimize", "nvt", "npt"];
    private bool _engLammps = true, _engGromacs = true, _engMinFirst = true, _engBusy;
    private int _engRun = 3;
    private double _engTemp = 300, _engPress = 1, _engDt = 0;
    private long _engSteps = 100000;
    private string _engFolder = "", _engStem = "system", _engError = "", _engPreviewName = "";
    public ObservableCollection<EngineFile> EngineLammpsFiles { get; } = new();
    public ObservableCollection<EngineFile> EngineGromacsFiles { get; } = new();
    public ObservableCollection<string> EngineFileNames { get; } = new();
    public ObservableCollection<PreviewLine> EnginePreviewLines { get; } = new();
    public ObservableCollection<Row> EngineChecks { get; } = new();
    public ObservableCollection<string> EngineNotes { get; } = new();
    private readonly Dictionary<string, List<string>> _engHeads = new();

    public bool EngineLammps { get => _engLammps; set { if (Set(ref _engLammps, value)) { Raise(nameof(EngineLammpsRefused)); RefreshEngines(); } } }
    public bool EngineGromacs { get => _engGromacs; set { if (Set(ref _engGromacs, value)) { Raise(nameof(EngineGromacsRefused)); RefreshEngines(); } } }
    public int EngineRun { get => _engRun; set { if (Set(ref _engRun, Math.Clamp(value, 0, 3))) { Raise(nameof(EngineIsMd)); RefreshEngines(); } } }
    public bool EngineIsMd => _engRun >= 2;
    public bool EngineMinimiseFirst { get => _engMinFirst; set { if (Set(ref _engMinFirst, value)) RefreshEngines(); } }
    private int _engConstraints;
    public int EngineConstraints { get => _engConstraints; set { if (Set(ref _engConstraints, Math.Clamp(value, 0, 2))) RefreshEngines(); } }
    public double EngineTemperature { get => _engTemp; set { if (Set(ref _engTemp, value)) RefreshEngines(); } }
    public double EnginePressure { get => _engPress; set { if (Set(ref _engPress, value)) RefreshEngines(); } }
    public double EngineDt { get => _engDt; set { if (Set(ref _engDt, value)) RefreshEngines(); } }
    public long EngineSteps { get => _engSteps; set { if (Set(ref _engSteps, value)) RefreshEngines(); } }
    public decimal? EngineTemperatureD { get => (decimal)_engTemp; set => EngineTemperature = Math.Clamp((double)(value ?? 300m), 1, 5000); }
    public decimal? EnginePressureD { get => (decimal)_engPress; set => EnginePressure = Math.Clamp((double)(value ?? 1m), 0, 1e5); }
    public decimal? EngineDtD { get => (decimal)_engDt; set => EngineDt = Math.Clamp((double)(value ?? 0m), 0, 50); }
    public decimal? EngineStepsD { get => _engSteps; set => EngineSteps = (long)Math.Clamp(value ?? 100000m, 0m, 1_000_000_000m); }
    public string EngineFolder { get => _engFolder; set { if (Set(ref _engFolder, value)) Raise(nameof(EngineCanWrite)); } }
    public string EngineStem { get => _engStem; set { if (Set(ref _engStem, value)) RefreshEngines(); } }
    private string _engLammpsError = "";
    /// <summary>Why LAMMPS cannot take this force field (GROMOS's reaction field …): shown in the LAMMPS card; GROMACS still exports.</summary>
    public string EngineLammpsError { get => _engLammpsError; private set { if (Set(ref _engLammpsError, value)) Raise(nameof(EngineLammpsRefused)); } }
    public bool EngineLammpsRefused => _engLammpsError.Length > 0 && _engLammps;
    private string _engGromacsError = "";
    /// <summary>Why GROMACS cannot take this force field (class II …): shown in the GROMACS card; LAMMPS still exports.</summary>
    public string EngineGromacsError { get => _engGromacsError; private set { if (Set(ref _engGromacsError, value)) Raise(nameof(EngineGromacsRefused)); } }
    public bool EngineGromacsRefused => _engGromacsError.Length > 0 && _engGromacs;
    public string EngineError { get => _engError; private set { if (Set(ref _engError, value)) { Raise(nameof(EngineHasError)); Raise(nameof(EngineCanWrite)); } } }
    public bool EngineHasError => _engError.Length > 0;
    public bool EngineBusy { get => _engBusy; private set { if (Set(ref _engBusy, value)) Raise(nameof(EngineCanWrite)); } }
    public bool EngineCanWrite => !_engBusy && !EngineHasError && _engFolder.Length > 0 && (_engLammps || _engGromacs);
    public string EngineWriteLabel => $"Write {EngineLammpsFiles.Count + EngineGromacsFiles.Count} files";
    public string EngineSummary { get; private set; } = "";
    public string EnginePreviewName
    {
        get => _engPreviewName;
        set { if (Set(ref _engPreviewName, value)) ShowEngineHead(); }
    }

    public void OpenExportCenter()
    {
        if (_doc == null) { Status = "Open or build a structure first"; return; }
        if (_engFolder.Length == 0)
        {
            // next to the file when there is one, else Documents/CAPS/<title>_export (a grown cell is not saved yet)
            var p = _doc.Path;
            var name = string.Concat(Title.Replace(" (unsaved)", "").Select(ch => char.IsLetterOrDigit(ch) || ch is '_' or '-' or '.' ? ch : '_'));
            EngineFolder = p is { Length: > 0 } && System.IO.File.Exists(p)
                ? System.IO.Path.Combine(System.IO.Path.GetDirectoryName(p) ?? "", System.IO.Path.GetFileNameWithoutExtension(p) + "_export")
                : System.IO.Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.MyDocuments), "CAPS", (name.Length > 0 ? name : "structure") + "_export");
        }
        SetModule(68);
        RefreshEngines();
    }

    // LAMMPS styles: the force field's own (its dihedral style, long-range Coulomb by PPPM …) or CAPS-exact (the energy
    // CAPS computes, for checking); hybrid; the long-range sum; the cut-off (0: the force field's) and k-space accuracy
    public static readonly string[] EngineStyleModes = ["Force field's own styles", "CAPS-exact (verification)"];
    public static readonly string[] EngineCoulombModes = ["Automatic", "PPPM", "Ewald", "Damped shifted force", "Plain cut-off"];
    private static readonly string[] EngineCoulombIds = ["auto", "pppm", "ewald", "dsf", "cut"];
    private int _engStyle, _engCoulomb;
    private bool _engHybrid;
    private double _engCutoff, _engKspace = 1e-4;
    public int EngineStyle { get => _engStyle; set { if (Set(ref _engStyle, Math.Clamp(value, 0, 1))) { Raise(nameof(EngineNative)); RefreshEngines(); } } }
    public bool EngineNative => _engStyle == 0;
    public bool EngineHybrid { get => _engHybrid; set { if (Set(ref _engHybrid, value)) RefreshEngines(); } }
    public int EngineCoulomb { get => _engCoulomb; set { if (Set(ref _engCoulomb, Math.Clamp(value, 0, 4))) { Raise(nameof(EngineKspaceVisible)); RefreshEngines(); } } }
    public bool EngineKspaceVisible => _engStyle == 0 && _engCoulomb <= 2;
    public decimal? EngineCutoffD { get => (decimal)_engCutoff; set { var v = Math.Clamp((double)(value ?? 0m), 0, 50); if (Math.Abs(v - _engCutoff) > 1e-12) { _engCutoff = v; Raise(); RefreshEngines(); } } }
    public string EngineKspaceText { get => _engKspace.ToString("0.##E+0", CultureInfo.InvariantCulture); set { if (double.TryParse(value, NumberStyles.Float, CultureInfo.InvariantCulture, out var v) && v > 0 && v < 1) { _engKspace = v; Raise(); RefreshEngines(); } } }

    private string EngineOptions(bool preview) => new JsonObject
    {
        ["lammps"] = _engLammps, ["gromacs"] = _engGromacs, ["stem"] = _engStem, ["run"] = EngineRunIds[_engRun],
        ["minimize_first"] = _engMinFirst, ["temperature"] = _engTemp, ["pressure"] = _engPress, ["dt"] = _engDt, ["steps"] = _engSteps, ["constraints"] = _engConstraints switch { 1 => "h-bonds", 2 => "all-bonds", _ => "none" },
        ["lammps_styles"] = _engStyle == 0 ? "native" : "exact", ["hybrid"] = _engHybrid, ["coulomb"] = EngineCoulombIds[_engCoulomb],
        ["cutoff"] = _engCutoff, ["kspace_accuracy"] = _engKspace,
        ["preview"] = preview, ["head_lines"] = preview ? 60 : 0,
    }.ToJsonString();

    private int _engGen;
    /// <summary>Writes the files to a scratch folder off the UI thread: the list, each file's head and the checks.</summary>
    public void RefreshEngines()
    {
        if (_doc == null || !IsExportCenter) return;
        var gen = ++_engGen;
        var doc = _doc;
        var opts = EngineOptions(true);
        EngineBusy = true;
        Task.Run(() =>
        {
            string json;
            try { json = doc.ExportEngines("", opts); } catch (Exception e) { json = new JsonObject { ["ok"] = false, ["error"] = e.Message }.ToJsonString(); }
            Avalonia.Threading.Dispatcher.UIThread.Post(() => { if (gen == _engGen) LoadEngines(json); });
        });
    }

    /// <summary>The same on the calling thread (scripts and the self-test).</summary>
    public void RefreshEnginesNow()
    {
        if (_doc == null) return;
        ++_engGen;
        string json;
        try { json = _doc.ExportEngines("", EngineOptions(true)); } catch (Exception e) { json = new JsonObject { ["ok"] = false, ["error"] = e.Message }.ToJsonString(); }
        LoadEngines(json);
    }

    private void LoadEngines(string json)
    {
        EngineBusy = false;
        EngineLammpsFiles.Clear(); EngineGromacsFiles.Clear(); EngineChecks.Clear(); EngineNotes.Clear(); _engHeads.Clear();
        var inv = CultureInfo.InvariantCulture;
        var j = JsonNode.Parse(json)!;
        if (j["ok"]?.GetValue<bool>() != true)
        {
            EngineError = (string?)j["error"] ?? "export failed";
            EngineFileNames.Clear(); EnginePreviewLines.Clear();
            EngineSummary = ""; Raise(nameof(EngineSummary)); Raise(nameof(EngineWriteLabel));
            return;
        }
        EngineError = "";
        EngineGromacsError = (string?)j["gromacs_error"] ?? "";
        EngineLammpsError = (string?)j["lammps_error"] ?? "";
        long total = 0;
        var names = new List<string>();
        foreach (var f in (JsonArray)j["files"]!)
        {
            var name = (string?)f!["name"] ?? "";
            var bytes = (long)((double?)f["bytes"] ?? 0);
            total += bytes;
            var ef = new EngineFile(name, (string?)f["what"] ?? "", FileSize(bytes));
            if (name.EndsWith(".data") || name.EndsWith(".in")) EngineLammpsFiles.Add(ef); else EngineGromacsFiles.Add(ef);
            _engHeads[name] = ((JsonArray?)f["head"] ?? new JsonArray()).Select(l => (string?)l ?? "").ToList();
            names.Add(name);
        }
        // the file tabs: kept in order; the one shown stays when it is still written
        if (!names.SequenceEqual(EngineFileNames)) { EngineFileNames.Clear(); foreach (var n in names) EngineFileNames.Add(n); }
        if (!names.Contains(_engPreviewName)) _engPreviewName = names.FirstOrDefault(n => n.EndsWith(".in")) ?? names.FirstOrDefault() ?? "";
        Raise(nameof(EnginePreviewName));
        ShowEngineHead();
        var c = j["checks"]!;
        double N(string k) => (double?)c[k] ?? 0;
        EngineChecks.Add(new Row("Atoms typed", $"{N("typed").ToString("N0", inv)} of {N("atoms").ToString("N0", inv)} · {N("types"):0} types"));
        EngineChecks.Add(new Row("Every term found", $"{N("bonds").ToString("N0", inv)} bonds · {N("angles").ToString("N0", inv)} angles · {N("dihedrals").ToString("N0", inv)} torsion terms · {N("impropers").ToString("N0", inv)} impropers"));
        EngineChecks.Add(new Row("Net charge", $"{N("net_charge").ToString("0.000", inv)} e · {ChargeWord((string?)c["charges"] ?? "")}"));
        EngineChecks.Add(new Row("Pair coefficients", $"all {N("type_pairs"):0} type pairs, mixing applied by CAPS"));
        if (N("density") > 0) EngineChecks.Add(new Row("Density", N("density").ToString("0.000", inv) + " g/cm³" + (_pipeDone.Contains("Relax") || _pipeDone.Contains("Equilibrate") || _pipeDone.Contains("Dynamics") ? "" : " · not relaxed yet")));
        foreach (var n in (JsonArray)j["notes"]!) EngineNotes.Add((string?)n ?? "");
        EngineSummary = $"{(string?)c["forcefield"]} · {names.Count} files · {FileSize(total)}";
        Raise(nameof(EngineSummary));
        Raise(nameof(EngineWriteLabel));
    }

    private static string ChargeWord(string c) => c switch { "automatic" => "the force field's, else Gasteiger", "types" => "the force field's", "gasteiger" => "Gasteiger–Marsili", "keep" => "the file's", "qeq" => "QEq", _ => c };
    private static string FileSize(long b) => b >= 1 << 20 ? (b / 1048576.0).ToString("0.0", CultureInfo.InvariantCulture) + " MB"
        : b >= 10240 ? (b / 1024.0).ToString("0", CultureInfo.InvariantCulture) + " KB" : (b / 1024.0).ToString("0.0", CultureInfo.InvariantCulture) + " KB";

    private void ShowEngineHead()
    {
        EnginePreviewLines.Clear();
        if (!_engHeads.TryGetValue(_engPreviewName, out var lines)) return;
        var k = 0;
        foreach (var l in lines) EnginePreviewLines.Add(new PreviewLine((++k).ToString(CultureInfo.InvariantCulture), l));
    }

    public async Task WriteEngines()
    {
        if (_doc == null || !EngineCanWrite) return;
        var doc = _doc;
        var dir = _engFolder;
        var opts = EngineOptions(false);
        EngineBusy = true;
        try
        {
            var json = await Task.Run(() => doc.ExportEngines(dir, opts));
            var j = JsonNode.Parse(json)!;
            if (j["ok"]?.GetValue<bool>() != true) { EngineError = (string?)j["error"] ?? "export failed"; Status = "Export failed: " + EngineError; return; }
            var n = ((JsonArray)j["files"]!).Count;
            Status = $"Wrote {n} files to {dir}" + (_engLammps ? $" · lmp -in {_engStem}.in" : "") + (_engGromacs ? $" · gmx grompp -f {_engStem}.mdp -c {_engStem}.gro -p {_engStem}.top" : "");
            MarkPipeline("Export");
            Record($"doc.export_engines({PyStr(dir)}, stem={PyStr(_engStem)}, lammps={(_engLammps ? "True" : "False")}, gromacs={(_engGromacs ? "True" : "False")}, run={PyStr(EngineRunIds[_engRun])})");
        }
        catch (Exception e) { EngineError = e.Message; Status = "Export failed: " + e.Message; }
        finally { EngineBusy = false; }
    }
}
