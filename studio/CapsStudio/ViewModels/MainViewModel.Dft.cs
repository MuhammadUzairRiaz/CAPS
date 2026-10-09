using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A species and its fraction on one face of a slab (2D sheets page).</summary>
public sealed class TermRow : ObservableObject
{
    private string _species = "O";
    private decimal _fraction = 1;
    public string Species { get => _species; set => Set(ref _species, value ?? "O"); }
    public decimal Fraction { get => _fraction; set => Set(ref _fraction, Math.Clamp(value, 0, 1000)); }
}
/// <summary>A validator finding: level, text, the atoms it is about (0-based).</summary>
public sealed record DftFinding(string Level, string Text, int[] Atoms)
{
    public IBrush Brush => Tokens.Brush(Level == "ERROR" ? "ErrB" : "WarnB");
}
/// <summary>A surface site for the site map: in-plane position, kind, the atom beneath, r, the termination's height.</summary>
public sealed record DftSite(double X, double Y, string Kind, string Face, string Beneath, double R, double H);
/// <summary>A complex of an adsorption set.</summary>
public sealed record DftComplexRow(string Name, string Mode, double Azimuth, double Any, double Heavy, double Image, string Status, string EquivalentTo, string Issues, int Suggested)
{
    public string AnyText => Any.ToString("0.00", CultureInfo.InvariantCulture);
    public string HeavyText => Heavy.ToString("0.00", CultureInfo.InvariantCulture);
    public string ImageText => Image.ToString("0.00", CultureInfo.InvariantCulture);
    public bool Equivalent => EquivalentTo.Length > 0;
    public string EquivalentText => Equivalent ? "≡ " + EquivalentTo.Replace("complex_", "") : "";
    public double RowOpacity => Equivalent ? 0.55 : 1.0;
    public IBrush StatusBrush => Tokens.Brush(Status == "ok" ? "OkB" : "ErrB");
}
/// <summary>An anchor a molecule can bind through (checkbox: used as an orientation).</summary>
public sealed class DftAnchorRow : ObservableObject
{
    public string Name { get; init; } = "";
    public string Atoms { get; init; } = "";
    private bool _on = true;
    public bool On { get => _on; set => Set(ref _on, value); }
}
/// <summary>An INCAR line with the reason for it.</summary>
public sealed record IncarRow(string File, string Key, string Value, string Why);
/// <summary>A VASP case folder in the run monitor.</summary>
public sealed record DftRunRow(string Case, string Path, string State, string Done, string Steps, string Energy, string Force, string Eta, string Flags, string Actions, double[] Fmax)
{
    public bool Flagged => Flags.Length > 0;
    public IBrush FlagBrush => Tokens.Brush(Flagged ? "ErrB" : State.StartsWith("running") ? "AccB" : "OkB");
}
/// <summary>A row of a results table: name and values as text.</summary>
public sealed record DftTableRow(string[] Cells);

/// <summary>The DFT surface &amp; adsorption workbench (design: CAPS_DFT_surface_workbench): 2D sheets and terminations,
/// adsorption sets, VASP calculation sets, the run monitor and the results — every action a `caps` command
/// (caps/dft_commands.hpp) through caps_dft_run, so the Studio, the CLI and Python give the same files; the command line
/// of the last action is shown to copy.</summary>
public partial class MainViewModel
{
    public bool IsSheets => _module == 75;
    public bool IsAdsorbDft => _module == 76;
    public bool IsDftJob => _module == 77;
    public bool IsDftRuns => _module == 78;
    public bool IsDftResults => _module == 79;
    public static readonly CultureInfo DftInv = CultureInfo.InvariantCulture;

    // ---------------------------------------------------------------- common
    public string DftDataDir => Paths.Data ?? "";
    public string DftRoot
    {
        get => _settings.DftRoot.Length > 0 ? _settings.DftRoot : Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "CAPS", "dft");
        set { if (_settings.DftRoot == value) return; _settings.DftRoot = value ?? ""; Raise(); Changed("DFT folder"); }
    }
    private string _dftCli = "", _dftStatus = "";
    /// <summary>The command line of the last DFT action (copyable: the same result from a terminal).</summary>
    public string DftCli { get => _dftCli; private set { if (Set(ref _dftCli, value)) Raise(nameof(HasDftCli)); } }
    public bool HasDftCli => _dftCli.Length > 0;
    public string DftStatus { get => _dftStatus; private set => Set(ref _dftStatus, value); }
    private bool _dftBusy;
    public bool DftBusy { get => _dftBusy; private set { if (Set(ref _dftBusy, value)) Raise(nameof(DftIdle)); } }
    public bool DftIdle => !_dftBusy;
    public event Action? DftChanged;

    /// <summary>Runs a workbench command (off the UI thread), records its command line; null on failure (DftStatus says why).</summary>
    public async Task<JsonObject?> Dft(string command, JsonObject args)
    {
        if (DftDataDir.Length == 0) { DftStatus = "CAPS data folder not found (data/sheets)"; return null; }
        DftBusy = true;
        try
        {
            var text = await Task.Run(() => CapsDocument.DftRun(command, args.ToJsonString(), DftDataDir));
            var r = JsonNode.Parse(text) as JsonObject;
            DftCli = (string?)r?["command"] ?? "";
            return r;
        }
        catch (Exception e) { DftStatus = $"{command}: {e.Message}"; return null; }
        finally { DftBusy = false; }
    }
    internal JsonObject? DftNow(string command, JsonObject args)
    {
        try
        {
            var r = JsonNode.Parse(CapsDocument.DftRun(command, args.ToJsonString(), DftDataDir)) as JsonObject;
            DftCli = (string?)r?["command"] ?? "";
            return r;
        }
        catch (Exception e) { DftStatus = $"{command}: {e.Message}"; return null; }
    }
    private string DftWork(params string[] parts) { var d = Path.Combine(new[] { DftRoot }.Concat(parts).ToArray()); Directory.CreateDirectory(d); return d; }
    private static double DNum(JsonNode? n, double def = double.NaN) => n is JsonValue v && v.TryGetValue<double>(out var d) ? d : def;
    private static string F(double v, string f = "0.000") => double.IsFinite(v) ? v.ToString(f, DftInv) : "–";

    // ---------------------------------------------------------------- 2D sheets & terminations (module 75)
    public ObservableCollection<string> SheetPresets { get; } = new();
    public ObservableCollection<string> TermSpecies { get; } = new();
    public static readonly string[] SiteKinds = ["fcc", "hcp", "top", "bridge"];
    public ObservableCollection<TermRow> TopTerms { get; } = new() { new TermRow { Species = "O", Fraction = 1 } };
    public ObservableCollection<TermRow> BottomTerms { get; } = new() { new TermRow { Species = "O", Fraction = 1 } };
    private string _sheetPreset = "Ti3C2", _sheetFrom = "", _sheetFormula = "", _sheetRemove = "", _siteTop = "fcc", _siteBottom = "fcc", _sheetOrder = "Ti,C,N,O,F,H";
    private decimal _sheetA, _sheetVacuum = 20, _sheetSeed = 1;
    private int _sheetNa = 1, _sheetNb = 1;
    private bool _sheetJanus;
    public string SheetPreset { get => _sheetPreset; set { if (value != null) Set(ref _sheetPreset, value); } }
    /// <summary>A structure to cut one layer out of ("" = the preset).</summary>
    public string SheetFrom { get => _sheetFrom; set { if (Set(ref _sheetFrom, value ?? "")) Raise(nameof(SheetFromPreset)); } }
    public bool SheetFromPreset => _sheetFrom.Length == 0;
    public string SheetFormula { get => _sheetFormula; set => Set(ref _sheetFormula, value ?? ""); }
    public string SheetRemove { get => _sheetRemove; set => Set(ref _sheetRemove, value ?? ""); }
    public decimal SheetA { get => _sheetA; set => Set(ref _sheetA, Math.Clamp(value, 0, 100)); }
    public string SiteTop { get => _siteTop; set => Set(ref _siteTop, value ?? "fcc"); }
    public string SiteBottom { get => _siteBottom; set => Set(ref _siteBottom, value ?? "fcc"); }
    public bool SheetJanus { get => _sheetJanus; set => Set(ref _sheetJanus, value); }
    public int SheetNa { get => _sheetNa; set => Set(ref _sheetNa, Math.Clamp(value, 1, 20)); }
    public int SheetNb { get => _sheetNb; set => Set(ref _sheetNb, Math.Clamp(value, 1, 20)); }
    public decimal SheetVacuum { get => _sheetVacuum; set => Set(ref _sheetVacuum, Math.Clamp(value, 0, 200)); }
    public decimal SheetSeed { get => _sheetSeed; set => Set(ref _sheetSeed, Math.Clamp(Math.Round(value), 0, int.MaxValue)); }
    public string SheetOrder { get => _sheetOrder; set => Set(ref _sheetOrder, value ?? ""); }
    public ObservableCollection<DftFinding> SheetFindings { get; } = new();
    public ObservableCollection<DftSite> SheetSites { get; } = new();
    public ObservableCollection<Row> SheetInfo { get; } = new();
    private string _sheetStatus = "", _sheetReport = "", _slabPath = "", _sheetPath = "";
    public string SheetValidation { get => _sheetStatus; private set { if (Set(ref _sheetStatus, value)) { Raise(nameof(SheetPass)); Raise(nameof(SheetBadge)); } } }
    public bool SheetPass => _sheetStatus == "PASS";
    public IBrush SheetBadge => Tokens.Brush(_sheetStatus == "PASS" ? "OkB" : _sheetStatus == "FAIL" ? "ErrB" : "DimB");
    public string SheetReport { get => _sheetReport; private set => Set(ref _sheetReport, value); }
    /// <summary>The terminated slab last built (a POSCAR in the DFT folder).</summary>
    public string SlabPath { get => _slabPath; private set { if (Set(ref _slabPath, value)) Raise(nameof(HasSlab)); } }
    public bool HasSlab => _slabPath.Length > 0 && File.Exists(_slabPath);
    public CapsDocument? SlabDoc { get; private set; }
    public string SitesNote { get; private set; } = "";
    /// <summary>The cell of the bare sheet (for the site map): a and b in-plane vectors.</summary>
    public (double Ax, double Ay, double Bx, double By) SheetCell { get; private set; }

    public void OpenSheets()
    {
        if (SheetPresets.Count == 0 && DftDataDir.Length > 0)
        {
            try
            {
                foreach (var s in JsonNode.Parse(File.ReadAllText(Path.Combine(DftDataDir, "sheets", "sheets.json")))!["sheets"]!.AsArray()) SheetPresets.Add((string?)s!["name"] ?? "");
                foreach (var t in JsonNode.Parse(File.ReadAllText(Path.Combine(DftDataDir, "sheets", "terminations.json")))!["terminations"]!.AsArray()) TermSpecies.Add((string?)t!["name"] ?? "");
            }
            catch (Exception e) { DftStatus = "DFT data: " + e.Message; }
            if (!SheetPresets.Contains(_sheetPreset) && SheetPresets.Count > 0) _sheetPreset = SheetPresets[0];
        }
        SetModule(75);
        Avalonia.Threading.Dispatcher.UIThread.Post(() => Raise(nameof(SheetPreset)));
    }
    public void AddTerm(bool top) => (top ? TopTerms : BottomTerms).Add(new TermRow { Species = "OH", Fraction = 0.25m });
    public void RemoveTerm(TermRow r) { TopTerms.Remove(r); BottomTerms.Remove(r); }
    private static string Fractions(IEnumerable<TermRow> rows) => string.Join(",", rows.Where(r => r.Fraction > 0).Select(r => $"{r.Species}:{r.Fraction.ToString("0.####", DftInv)}"));

    /// <summary>Builds the sheet (preset or cut from a file), its sites, the terminated slab, and validates it.</summary>
    public async Task BuildSlab()
    {
        var work = DftWork("work");
        _sheetPath = Path.Combine(work, "sheet.vasp");
        var sa = new JsonObject { ["o"] = _sheetPath, ["order"] = _sheetOrder };
        if (_sheetFrom.Length > 0)
        {
            sa["from"] = _sheetFrom;
            if (_sheetFormula.Length > 0) sa["formula"] = _sheetFormula;
            if (_sheetRemove.Length > 0) sa["remove"] = _sheetRemove;
        }
        else
        {
            sa["preset"] = _sheetPreset;
            if (_sheetA > 0) sa["a"] = (double)_sheetA;
        }
        var sheet = await Dft("sheet", sa);
        if (sheet == null) return;
        var label = $"{(string?)sheet["formula"]}_{string.Join("", TopTerms.Select(t => t.Species))}";
        SlabPath = "";
        var ta = TerminateArgs(Path.Combine(work, label + ".vasp"));
        var r = await Dft("terminate", ta);
        if (r == null) return;
        SlabPath = (string?)r["wrote"] ?? "";
        LoadValidation(r["validation"] as JsonObject, (string?)r["report"] ?? "");
        SheetInfo.Clear();
        SheetInfo.Add(new Row("Formula", (string?)r["formula"] ?? ""));
        SheetInfo.Add(new Row("Atoms", ((double?)r["atoms"] ?? 0).ToString("0", DftInv)));
        SheetInfo.Add(new Row("a (per sheet cell)", F(DNum(sheet["a"]), "0.0000") + " Å"));
        SheetInfo.Add(new Row("Site in-plane r", F(DNum(r["r_site"]), "0.0000") + " Å"));
        foreach (var n in (sheet["notes"] as JsonArray ?? [])) SheetInfo.Add(new Row("Note", (string?)n ?? ""));
        await LoadSites();
        try { SlabDoc?.Dispose(); SlabDoc = CapsDocument.Open(SlabPath); } catch (Exception e) { DftStatus = "preview: " + e.Message; }
        DftStatus = $"{(string?)r["formula"]} · {SheetValidation}";
        DftCli = (string?)r["command"] ?? DftCli;
        DftChanged?.Invoke();
    }

    private JsonObject TerminateArgs(string output)
    {
        var a = new JsonObject
        {
            ["inputs"] = new JsonArray(_sheetPath), ["top"] = Fractions(TopTerms), ["site"] = _siteTop, ["supercell"] = $"{_sheetNa}x{_sheetNb}",
            ["vacuum"] = (double)_sheetVacuum, ["seed"] = (double)_sheetSeed, ["order"] = _sheetOrder, ["o"] = output,
        };
        if (_sheetJanus) { a["janus"] = true; a["bottom"] = Fractions(BottomTerms); a["site_bottom"] = _siteBottom; }
        return a;
    }

    private void LoadValidation(JsonObject? v, string report)
    {
        SheetFindings.Clear();
        SheetValidation = (string?)v?["status"] ?? "";
        SheetReport = report;
        foreach (var f in (v?["findings"] as JsonArray ?? []).OfType<JsonObject>())
            SheetFindings.Add(new DftFinding((string?)f["level"] ?? "", (string?)f["text"] ?? "", (f["atoms"] as JsonArray ?? []).Select(x => (int)DNum(x, 0)).ToArray()));
    }

    private async Task LoadSites()
    {
        SheetSites.Clear();
        var r = await Dft("sites", new JsonObject { ["inputs"] = new JsonArray(_sheetPath), ["species"] = TopTerms.FirstOrDefault()?.Species ?? "O" });
        if (r == null) return;
        try
        {
            var s = CapsDocument.Open(_sheetPath);
            var sum = s.Summary();
            SheetCell = (sum.CellA, 0, -sum.CellB / 2, sum.CellB * Math.Sqrt(3) / 2);
            s.Dispose();
        }
        catch { }
        foreach (var x in (r["sites"] as JsonArray ?? []).OfType<JsonObject>())
        {
            var p = x["pos"] as JsonArray;
            SheetSites.Add(new DftSite(DNum(p?[0]), DNum(p?[1]), (string?)x["kind"] ?? "", (string?)x["face"] ?? "", (string?)x["beneath"] ?? "", DNum(x["r"]), DNum(x["h"])));
        }
        SitesNote = r["bond"] is JsonNode b ? $"{TopTerms.FirstOrDefault()?.Species}: bond {F(DNum(b), "0.00")} Å ({(string?)r["bond_source"]}); h = √(d² − r²) above the {(string?)r["surface_element"]} plane" : "";
        Raise(nameof(SitesNote));
    }

    /// <summary>Writes the slab as POSCAR / CIF (by extension), or a VASP set into a folder (vasp_set).</summary>
    public async Task<string?> ExportSlab(string path, bool vaspSet)
    {
        if (!HasSlab) { DftStatus = "Build a slab first"; return null; }
        var a = TerminateArgs(vaspSet ? Path.Combine(path, "POSCAR") : path);
        if (vaspSet)
        {
            a["vasp_set"] = path;
            a["kdens"] = (double)_kdens;
            a["profile"] = DftProfile;
            if (_settings.PotcarDir.Length > 0) a["pp_dir"] = _settings.PotcarDir;
        }
        var r = await Dft("terminate", a);
        if (r == null) return null;
        DftStatus = vaspSet ? $"VASP set written to {path} ({F(DNum(r["vasp_set"]?["encut"]), "0")} eV, k {string.Join("×", (r["vasp_set"]?["kpoints"] as JsonArray ?? []).Select(x => x?.ToString()))})" : $"Wrote {path}";
        return path;
    }

    // ---------------------------------------------------------------- adsorption set (module 76)
    private string _adsSlab = "", _adsSmiles = "C/C=C/CCC(C#N)C/C=C/C", _adsMolFile = "", _adsAzimuths = "0,90,180,270", _adsOut = "", _adsSymmetry = "";
    private bool _adsFromRelaxed, _adsParallel = true, _adsUpright, _adsSkipEq;
    private int _adsNa = 5, _adsNb = 5;
    private decimal _adsDmin = 2.3m, _adsDheavy = 3.0m;
    public string AdsSlab { get => _adsSlab; set => Set(ref _adsSlab, value ?? ""); }
    public bool AdsFromRelaxed { get => _adsFromRelaxed; set => Set(ref _adsFromRelaxed, value); }
    public int AdsNa { get => _adsNa; set => Set(ref _adsNa, Math.Clamp(value, 1, 20)); }
    public int AdsNb { get => _adsNb; set => Set(ref _adsNb, Math.Clamp(value, 1, 20)); }
    public string AdsSmiles { get => _adsSmiles; set => Set(ref _adsSmiles, value ?? ""); }
    public string AdsMolFile { get => _adsMolFile; set => Set(ref _adsMolFile, value ?? ""); }
    public string AdsAzimuths { get => _adsAzimuths; set => Set(ref _adsAzimuths, value ?? ""); }
    public bool AdsParallel { get => _adsParallel; set => Set(ref _adsParallel, value); }
    public bool AdsUpright { get => _adsUpright; set => Set(ref _adsUpright, value); }
    public bool AdsSkipEquivalent { get => _adsSkipEq; set => Set(ref _adsSkipEq, value); }
    private int _adsPrescreen;
    /// <summary>Also the N lowest configurations of the force-field Adsorption Locator (UFF annealing) as complexes.</summary>
    public int AdsPrescreen { get => _adsPrescreen; set => Set(ref _adsPrescreen, Math.Clamp(value, 0, 20)); }
    public decimal AdsDmin { get => _adsDmin; set => Set(ref _adsDmin, Math.Clamp(value, 0.5m, 10)); }
    public decimal AdsDheavy { get => _adsDheavy; set => Set(ref _adsDheavy, Math.Clamp(value, 0.5m, 10)); }
    public string AdsOut { get => _adsOut; set => Set(ref _adsOut, value ?? ""); }
    public string AdsSymmetry { get => _adsSymmetry; private set => Set(ref _adsSymmetry, value); }
    public ObservableCollection<DftAnchorRow> AdsAnchors { get; } = new();
    public ObservableCollection<DftComplexRow> AdsComplexes { get; } = new();
    public ObservableCollection<Row> AdsAudit { get; } = new();
    public bool HasAdsAudit => AdsAudit.Count > 0;
    public CapsDocument? AdsPreviewDoc { get; private set; }
    private string _adsSummary = "";
    public string AdsSummary { get => _adsSummary; private set => Set(ref _adsSummary, value); }

    public void OpenAdsorbDft()
    {
        if (_adsSlab.Length == 0 && HasSlab) AdsSlab = SlabPath;
        if (_adsOut.Length == 0) AdsOut = Path.Combine(DftRoot, "adsorption", "set");
        SetModule(76);
    }

    /// <summary>Builds the adsorption set (anchors found on the first run; then the checked ones) and writes the trio.</summary>
    public async Task BuildAdsorptionSet()
    {
        if (!File.Exists(_adsSlab)) { DftStatus = "Choose the slab (a POSCAR or a relaxed CONTCAR)"; return; }
        var a = new JsonObject
        {
            ["inputs"] = new JsonArray(_adsSlab), ["supercell"] = $"{_adsNa}x{_adsNb}", ["azimuths"] = _adsAzimuths, ["dmin"] = (double)_adsDmin, ["dheavy"] = (double)_adsDheavy,
            ["out"] = _adsOut, ["kdens"] = (double)_kdens, ["profile"] = DftProfile,
        };
        if (_adsFromRelaxed) a["from_relaxed"] = true;
        if (_adsSkipEq) a["skip_equivalent"] = true;
        if (_adsPrescreen > 0) a["prescreen"] = _adsPrescreen;
        if (_adsMolFile.Length > 0) a["molecule"] = _adsMolFile; else a["smiles"] = _adsSmiles;
        if (_settings.PotcarDir.Length > 0) a["pp_dir"] = _settings.PotcarDir;
        if (AdsAnchors.Count > 0)
        {
            var modes = AdsAnchors.Where(x => x.On).Select(x => "anchor:" + x.Name).ToList();
            if (_adsParallel) modes.Add("parallel");
            if (_adsUpright) modes.Add("upright");
            if (modes.Count == 0) { DftStatus = "Choose at least one orientation"; return; }
            a["modes"] = string.Join(",", modes);
        }
        if (Directory.Exists(_adsOut)) Directory.Delete(_adsOut, true);   // a set is rebuilt whole: no stale complexes
        var r = await Dft("adsorb-dft", a);
        if (r == null) return;
        var on = AdsAnchors.ToDictionary(x => x.Name, x => x.On);
        AdsAnchors.Clear();
        foreach (var (k, v) in (r["anchors"] as JsonObject ?? []))
            AdsAnchors.Add(new DftAnchorRow { Name = k, Atoms = string.Join(" ", (v as JsonArray ?? []).Select(x => "#" + ((int)DNum(x, 0) + 1))), On = !on.TryGetValue(k, out var o) || o });
        AdsComplexes.Clear();
        foreach (var c in (r["complexes"] as JsonArray ?? []).OfType<JsonObject>())
            AdsComplexes.Add(new DftComplexRow((string?)c["name"] ?? "", (string?)c["mode"] ?? "", DNum(c["azimuth"]), DNum(c["any"]), DNum(c["heavy"]), DNum(c["image"]), (string?)c["status"] ?? "",
                (string?)c["equivalent_to"] ?? "", string.Join("\n", (c["issues"] as JsonArray ?? []).Select(i => $"[{(string?)i?["level"]}] {(string?)i?["text"]}")), (int)DNum(c["suggested_supercell"], 0)));
        AdsSymmetry = string.Join(" ", (r["notes"] as JsonArray ?? []).Select(n => (string?)n ?? "")) + (r["prescreen_note"] is JsonNode pn ? " Pre-screen: " + (string?)pn + "." : "");
        AdsAudit.Clear();
        foreach (var x in (r["audit"] as JsonArray ?? []).OfType<JsonObject>())
            AdsAudit.Add(new Row(((string?)x["name"] ?? "").Replace("complex_", ""), $"{(int)DNum(x["n_mol"], 0)} atoms · any {F(DNum(x["any"]), "0.00")} · heavy {F(DNum(x["heavy"]), "0.00")} · image {F(DNum(x["image"]), "0.00")} · {((bool?)x["intact"] == true ? "intact" : "BROKEN")} · {((bool?)x["ok"] == true ? "ok" : "FAIL")}"));
        Raise(nameof(HasAdsAudit));
        var ok = AdsComplexes.Count(c => c.Status == "ok");
        AdsSummary = $"{AdsComplexes.Count} complexes of a {(int)DNum(r["molecule_atoms"], 0)}-atom molecule · {ok} pass · {AdsComplexes.Count(c => c.Equivalent)} equivalent by symmetry · written to {_adsOut}";
        DftStatus = AdsSummary;
        if (AdsComplexes.FirstOrDefault(c => c.Status == "ok") is { } first) PreviewComplex(first);
        DftChanged?.Invoke();
    }

    public void PreviewComplex(DftComplexRow c)
    {
        var p = Path.Combine(_adsOut, c.Name, "POSCAR");
        if (!File.Exists(p)) { DftStatus = $"{c.Name} was not written ({c.Issues.Split('\n').FirstOrDefault()})"; return; }
        try { AdsPreviewDoc?.Dispose(); AdsPreviewDoc = CapsDocument.Open(p); DftChanged?.Invoke(); } catch (Exception e) { DftStatus = "preview: " + e.Message; }
    }

    // ---------------------------------------------------------------- VASP set designer (module 77)
    public ObservableCollection<string> DftProfiles { get; } = new();
    public static readonly string[] DftStageSets = ["slab: 01_cell 02_cell2 03_relax 04_static", "fixed cell: 03_relax 04_static"];
    public static readonly string[] DftStaticOuts = ["none: energy + DOS", "pot: + LOCPOT (work function)", "all: + CHGCAR, AECCAR (Bader, CDD)"];
    public static readonly string[] DftDipoles = ["auto (faces differ)", "on", "off"];
    public static readonly string[] DftIvdws = ["12 · D3-BJ", "11 · D3 zero damping", "13 · D4", "0 · none"];
    private string _jobStructure = "", _jobLabel = "", _jobOut = "", _jobProfileText = "", _jobMesh = "", _jobPreviewFile = "", _jobPreviewText = "", _jobStore = "";
    private int _jobStages, _jobStaticOut = 1, _jobDipole, _jobIvdw, _jobFile;
    private decimal _kdens = 45, _jobEncut;
    private bool _jobGamma;
    public string JobStructure { get => _jobStructure; set { if (Set(ref _jobStructure, value ?? "")) UpdateMesh(); } }
    public string JobLabel { get => _jobLabel; set => Set(ref _jobLabel, value ?? ""); }
    public string JobOut { get => _jobOut; set => Set(ref _jobOut, value ?? ""); }
    public int JobStages { get => _jobStages; set => Set(ref _jobStages, value); }
    public int JobStaticOut { get => _jobStaticOut; set => Set(ref _jobStaticOut, value); }
    public int JobDipole { get => _jobDipole; set => Set(ref _jobDipole, value); }
    public int JobIvdw { get => _jobIvdw; set => Set(ref _jobIvdw, value); }
    public bool JobGamma { get => _jobGamma; set { if (Set(ref _jobGamma, value)) UpdateMesh(); } }
    public decimal JobEncut { get => _jobEncut; set => Set(ref _jobEncut, Math.Clamp(value, 0, 3000)); }
    /// <summary>k-point density (Å): N_i = ceil(kdens / |a_i|), 1 along the vacuum; the same density in every cell.</summary>
    public decimal Kdens { get => _kdens; set { if (Set(ref _kdens, Math.Clamp(Math.Round(value), 5, 200))) UpdateMesh(); } }
    public string JobMesh { get => _jobMesh; private set => Set(ref _jobMesh, value); }
    public string DftProfile
    {
        get => _settings.DftProfile.Length > 0 ? _settings.DftProfile : "slurm-workspace";
        set { if (value == null || _settings.DftProfile == value) return; _settings.DftProfile = value; Raise(); Changed("Cluster profile"); LoadProfileText(); }
    }
    public string JobProfileText { get => _jobProfileText; private set => Set(ref _jobProfileText, value); }
    /// <summary>Your licensed PAW directory (potpaw_PBE): POTCARs are assembled from it, never shipped with CAPS.</summary>
    public string PotcarDir { get => _settings.PotcarDir; set { if (_settings.PotcarDir == (value ?? "")) return; _settings.PotcarDir = value ?? ""; Raise(); Changed("PAW directory"); } }
    public ObservableCollection<IncarRow> JobIncar { get; } = new();
    public ObservableCollection<string> JobFiles { get; } = new();
    public int JobFile { get => _jobFile; set { if (Set(ref _jobFile, value)) ShowJobFile(); } }
    public string JobPreviewText { get => _jobPreviewText; private set => Set(ref _jobPreviewText, value); }
    public ObservableCollection<Row> JobFacts { get; } = new();
    public ObservableCollection<string> JobStageChips { get; } = new();
    /// <summary>Long-term storage path for the project: written to <project>/.store; jobs create stage folders there.</summary>
    public string JobStore { get => _jobStore; set => Set(ref _jobStore, value ?? ""); }
    private string _jobPreviewDir = "";

    public void OpenDftJob()
    {
        if (DftProfiles.Count == 0 && DftDataDir.Length > 0)
            try { foreach (var p in JsonNode.Parse(File.ReadAllText(Path.Combine(DftDataDir, "dft", "clusters.json")))!["profiles"]!.AsArray()) DftProfiles.Add((string?)p!["name"] ?? ""); } catch { }
        if (_jobStructure.Length == 0 && HasSlab) JobStructure = SlabPath;
        if (_jobOut.Length == 0) JobOut = Path.Combine(DftRoot, "structures", "case");
        LoadProfileText();
        UpdateMesh();
        SetModule(77);
    }

    private void LoadProfileText()
    {
        try
        {
            var p = JsonNode.Parse(File.ReadAllText(Path.Combine(DftDataDir, "dft", "clusters.json")))!["profiles"]!.AsArray().OfType<JsonObject>().FirstOrDefault(x => (string?)x["name"] == DftProfile);
            if (p == null) { JobProfileText = ""; return; }
            var rules = string.Join("\n", (p["parallel"] as JsonArray ?? []).OfType<JsonObject>().Select(r =>
                $"  {((int)DNum(r["max_atoms"], 0) == 0 ? "larger" : "≤ " + (int)DNum(r["max_atoms"], 0) + " atoms")}: {(int)DNum(r["ranks"], 0)} ranks · KPAR {(int)DNum(r["kpar"], 0)} · NCORE {(int)DNum(r["ncore"], 0)} · {(string?)r["time"]}"));
            JobProfileText = $"{(string?)p["about"]}\naccount {(string?)p["account"]} · modules {(string?)p["modules"]}\n{rules}\nno e-mail lines are written; edit data/dft/clusters.json for your cluster";
        }
        catch (Exception e) { JobProfileText = e.Message; }
    }

    private void UpdateMesh()
    {
        if (!File.Exists(_jobStructure)) { JobMesh = ""; return; }
        try
        {
            using var d = CapsDocument.Open(_jobStructure);
            var s = d.Summary();
            int K(double a) => _jobGamma ? 1 : Math.Max(1, (int)Math.Ceiling((double)_kdens / a));
            JobMesh = $"{K(s.CellA)} × {K(s.CellB)} × 1   (|a| {s.CellA.ToString("0.000", DftInv)} Å, |b| {s.CellB.ToString("0.000", DftInv)} Å, Gamma-centred)";
        }
        catch { JobMesh = ""; }
    }

    private JsonObject JobArgs(string outDir) => new()
    {
        ["inputs"] = new JsonArray(_jobStructure), ["out"] = outDir, ["label"] = _jobLabel.Length > 0 ? _jobLabel : Path.GetFileName(_jobOut.TrimEnd('/', '\\')),
        ["stages"] = _jobStages == 0 ? "slab" : "fixed", ["static_out"] = _jobStaticOut switch { 0 => "none", 2 => "all", _ => "pot" },
        ["dipole"] = _jobDipole switch { 1 => "on", 2 => "off", _ => "auto" }, ["ivdw"] = _jobIvdw switch { 1 => 11, 2 => 13, 3 => 0, _ => 12 },
        ["kdens"] = (double)_kdens, ["profile"] = DftProfile, ["gamma"] = _jobGamma,
    };

    /// <summary>The set written into a scratch folder: every file shown, the INCAR lines with their reasons.</summary>
    public async Task PreviewJob()
    {
        if (!File.Exists(_jobStructure)) { DftStatus = "Choose the structure"; return; }
        _jobPreviewDir = Path.Combine(Path.GetTempPath(), "caps-dft-preview", Path.GetFileName(_jobOut.TrimEnd('/', '\\')));
        if (Directory.Exists(_jobPreviewDir)) Directory.Delete(_jobPreviewDir, true);
        var a = JobArgs(_jobPreviewDir);
        if (_jobEncut > 0) a["encut"] = (double)_jobEncut;
        if (PotcarDir.Length > 0) a["pp_dir"] = PotcarDir;
        var r = await Dft("vasp-set", a);
        if (r == null) return;
        var set = r["set"] as JsonObject;
        JobFiles.Clear();
        foreach (var f in (set?["files"] as JsonArray ?? [])) JobFiles.Add((string?)f ?? "");
        JobIncar.Clear();
        foreach (var file in new[] { "INCAR.cell", "INCAR.relax", "INCAR.static" })
        {
            var p = Path.Combine(_jobPreviewDir, file);
            if (!File.Exists(p)) continue;
            foreach (var line in File.ReadAllLines(p).Skip(1))
            {
                var eq = line.IndexOf('=');
                if (eq < 0) continue;
                var hash = line.IndexOf("   # ", StringComparison.Ordinal);
                JobIncar.Add(new IncarRow(file.Replace("INCAR.", ""), line[..eq].Trim(), (hash > eq ? line[(eq + 1)..hash] : line[(eq + 1)..]).Trim(), hash > 0 ? line[(hash + 5)..] : ""));
            }
        }
        JobFacts.Clear();
        JobFacts.Add(new Row("ENCUT", $"{F(DNum(set?["encut"]), "0")} eV"));
        JobFacts.Add(new Row("k-mesh", string.Join(" × ", (set?["kpoints"] as JsonArray ?? []).Select(x => x?.ToString()))));
        JobFacts.Add(new Row("Parallel", $"{(int)DNum(set?["ranks"], 0)} ranks · KPAR {(int)DNum(set?["kpar"], 0)} · NCORE {(int)DNum(set?["ncore"], 0)} · {(string?)set?["time"]}"));
        JobFacts.Add(new Row("Validation", (string?)r["validation"] ?? ""));
        foreach (var n in (set?["notes"] as JsonArray ?? [])) JobFacts.Add(new Row("Note", (string?)n ?? ""));
        JobStageChips.Clear();
        foreach (var s in _jobStages == 0 ? new[] { "01_cell", "02_cell2", "03_relax", "04_static" } : new[] { "03_relax", "04_static" }) JobStageChips.Add(s);
        _jobFile = 0;
        ShowJobFile();
        Avalonia.Threading.Dispatcher.UIThread.Post(() => Raise(nameof(JobFile)));   // after the list is in the box
        DftStatus = $"Preview of {JobFiles.Count} files (not written yet): Write set puts them in {_jobOut}";
    }

    private void ShowJobFile()
    {
        if (_jobFile < 0 || _jobFile >= JobFiles.Count) { JobPreviewText = ""; return; }
        var p = Path.Combine(_jobPreviewDir, JobFiles[_jobFile]);
        JobPreviewText = JobFiles[_jobFile] == "POTCAR" ? "(POTCAR: your licensed PAW data — not shown)" : File.Exists(p) ? File.ReadAllText(p) : "";
    }

    public async Task WriteJob()
    {
        var a = JobArgs(_jobOut);
        if (_jobEncut > 0) a["encut"] = (double)_jobEncut;
        if (PotcarDir.Length > 0) a["pp_dir"] = PotcarDir;
        var r = await Dft("vasp-set", a);
        if (r == null) return;
        DftStatus = $"VASP set written to {_jobOut} · submit: cd {_jobOut} && sbatch job.slurm";
    }

    /// <summary>The project's long-term storage: <root>/.store holds the path; jobs then create their stage folders there.</summary>
    public void SetStore(string projectRoot)
    {
        if (_jobStore.Length == 0) { DftStatus = "Give the storage path first"; return; }
        File.WriteAllText(Path.Combine(projectRoot, ".store"), _jobStore + "\n");
        DftStatus = $"{Path.Combine(projectRoot, ".store")} → {_jobStore}: new jobs below it store their stages there (move existing ones with Runs › Store)";
    }

    // ---------------------------------------------------------------- run monitor (module 78)
    private string _runsRoot = "";
    public string RunsRoot { get => _runsRoot; set => Set(ref _runsRoot, value ?? ""); }
    public ObservableCollection<DftRunRow> DftRuns { get; } = new();
    private DftRunRow? _dftRun;
    public DftRunRow? SelectedDftRun { get => _dftRun; set { if (Set(ref _dftRun, value)) { Raise(nameof(HasDftRun)); LoadRunStages(); DftChanged?.Invoke(); } } }
    public bool HasDftRun => _dftRun != null;
    public ObservableCollection<string> DftRunStages { get; } = new();
    private string _dftRunStage = "", _dftPlan = "";
    public string DftRunStage { get => _dftRunStage; set => Set(ref _dftRunStage, value ?? ""); }
    public string DftPlan { get => _dftPlan; private set => Set(ref _dftPlan, value); }

    public void OpenDftRuns()
    {
        if (_runsRoot.Length == 0) RunsRoot = DftRoot;
        SetModule(78);
    }

    /// <summary>Every case folder (one with a job.slurm) below the root: progress, and for adsorption sets the health flags.</summary>
    public async Task RefreshRuns()
    {
        if (!Directory.Exists(_runsRoot)) { DftStatus = "Choose the project folder"; return; }
        var cases = new List<string>();
        void Walk(string d, int depth)
        {
            if (File.Exists(Path.Combine(d, "job.slurm"))) { cases.Add(d); return; }
            if (depth > 4) return;
            try { foreach (var s in Directory.GetDirectories(d).OrderBy(x => x, StringComparer.Ordinal)) if (!Path.GetFileName(s).StartsWith('.')) Walk(s, depth + 1); } catch { }
        }
        Walk(_runsRoot, 0);
        var prog = await Dft("vasp-progress", new JsonObject { ["inputs"] = new JsonArray(cases.Select(c => (JsonNode)c).ToArray()) });
        // health of every adsorption set (a folder with slab/, molecule/, complex_*)
        var flags = new Dictionary<string, (string Text, string Action)>();
        foreach (var set in cases.Select(c => Path.GetDirectoryName(c)!).Distinct().Where(s => Directory.Exists(Path.Combine(s, "slab")) && Directory.Exists(Path.Combine(s, "molecule"))))
        {
            var h = await Dft("vasp-health", new JsonObject { ["inputs"] = new JsonArray(set) });
            foreach (var row in (h?["health"]?["complexes"] as JsonArray ?? []).OfType<JsonObject>())
            {
                var fl = (row["flags"] as JsonArray ?? []).OfType<JsonObject>().ToList();
                if (fl.Count > 0)
                    flags[Path.Combine(set, (string?)row["complex"] ?? "")] = (string.Join(" · ", fl.Select(f => $"{(string?)f["kind"]} {(string?)f["text"]}")), string.Join("\n", fl.Select(f => $"{(string?)f["kind"]}: {(string?)f["action"]}")));
            }
        }
        DftRuns.Clear();
        foreach (var c in (prog?["cases"] as JsonArray ?? []).OfType<JsonObject>())
        {
            var path = (string?)c["case"] ?? "";
            var stages = (c["stages"] as JsonArray ?? []).Count;
            var done = (c["done"] as JsonArray ?? []).Count;
            var fm = (c["fmax"] as JsonArray ?? []).Select(x => DNum(x)).Where(double.IsFinite).ToArray();
            var last = (c["last_E0"] as JsonArray ?? []).Select(x => DNum(x)).LastOrDefault(double.NaN);
            var eta = DNum(c["eta_seconds"]);
            var state = (string?)c["state"] ?? "";
            if (c["running"] is JsonNode rn) state = $"running {(string?)rn}";
            var (ft, fa) = flags.TryGetValue(path, out var x) ? x : ("", "");
            if (c["note"] is JsonNode note && ft.Length == 0) (ft, fa) = ("NOT CONVERGED", (string?)note ?? "");
            DftRuns.Add(new DftRunRow(Path.GetRelativePath(_runsRoot, path), path, state, $"{done}/{stages}", ((int)DNum(c["ionic_steps"], 0)).ToString(DftInv), F(last, "0.0000"),
                F(DNum(c["max_force"]), "0.000"), double.IsFinite(eta) ? $"~{eta / 60:0} min" : "", ft, fa, fm));
        }
        DftStatus = $"{DftRuns.Count} case folders · {DftRuns.Count(r => r.State.StartsWith("running"))} running · {DftRuns.Count(r => r.Flagged)} flagged";
        DftChanged?.Invoke();
    }

    private void LoadRunStages()
    {
        DftRunStages.Clear();
        if (_dftRun == null) return;
        try { foreach (var d in Directory.GetDirectories(_dftRun.Path).Select(Path.GetFileName).Where(n => n != null && System.Text.RegularExpressions.Regex.IsMatch(n, @"^(0[1-4]_\w+|seg_\d+)$")).OrderBy(n => n)) DftRunStages.Add(d!); } catch { }
        DftRunStage = DftRunStages.LastOrDefault() ?? "";
    }

    public async Task ResubmitRun()
    {
        if (_dftRun == null) return;
        var r = await Dft("vasp-jobs", new JsonObject { ["inputs"] = new JsonArray("submit", _dftRun.Path), ["run"] = true });
        var sub = r?["submit"] as JsonArray ?? [];
        DftStatus = sub.Count == 0 ? $"{_dftRun.Case}: finished or already queued — not submitted" : $"{_dftRun.Case}: sbatch exit {(int)DNum(sub[0]?["exit"], -1)} (on a cluster login node; here run the command shown)";
    }

    public async Task ResetRunStage()
    {
        if (_dftRun == null || _dftRunStage.Length == 0) return;
        var r = await Dft("vasp-jobs", new JsonObject { ["inputs"] = new JsonArray("reset", _dftRun.Path, _dftRunStage) });
        if (r != null) { DftStatus = $"{_dftRunStage} put aside as {(string?)r["moved_to"]}: submit again to redo it"; LoadRunStages(); }
    }

    /// <summary>Cleanup: first a dry run (the plan), the second press removes.</summary>
    public async Task CleanupRuns(bool keepBest)
    {
        var apply = DftPlan.StartsWith("Cleanup plan", StringComparison.Ordinal);
        var a = new JsonObject { ["inputs"] = new JsonArray("cleanup", _runsRoot), ["keep_best"] = keepBest };
        if (apply) a["yes"] = true;
        var r = await Dft("vasp-jobs", a);
        if (r == null) return;
        var items = (r["plan"]?["items"] as JsonArray ?? []).Count;
        var gb = DNum(r["plan"]?["bytes"], 0) / 1e9;
        DftPlan = apply ? $"Removed {items} files ({gb:0.00} GB)." : $"Cleanup plan: {items} files ({gb:0.00} GB) of finished stages — press Clean up again to remove them";
        if (apply) DftPlan = "";
        DftStatus = apply ? $"Removed {items} files, {gb:0.00} GB" : DftPlan;
    }

    public void OpenRunStructure()
    {
        if (_dftRun == null) return;
        foreach (var st in new[] { _dftRunStage }.Concat(DftRunStages.Reverse()))
            foreach (var f in new[] { "CONTCAR", "POSCAR.in" })
            {
                var p = Path.Combine(_dftRun.Path, st, f);
                if (File.Exists(p) && new FileInfo(p).Length > 0) { Open(p); return; }
            }
        var poscar = Path.Combine(_dftRun.Path, "POSCAR");
        if (File.Exists(poscar)) Open(poscar);
    }

    // ---------------------------------------------------------------- results (module 79)
    private string _resultsSet = "", _resultsAimd = "", _resultsMain = "";
    public string ResultsSet { get => _resultsSet; set => Set(ref _resultsSet, value ?? ""); }
    public string ResultsAimd { get => _resultsAimd; set => Set(ref _resultsAimd, value ?? ""); }
    public string ResultsMain { get => _resultsMain; set => Set(ref _resultsMain, value ?? ""); }
    public ObservableCollection<DftTableRow> BindRows { get; } = new();
    public ObservableCollection<DftTableRow> GeomRows { get; } = new();
    public ObservableCollection<DftTableRow> DftSummaryRows { get; } = new();
    public ObservableCollection<Row> ResultsFacts { get; } = new();
    public (double X, double Y)[] BindBars { get; private set; } = [];
    public string[] BindNames { get; private set; } = [];
    public Dictionary<string, (double X, double Y)[]> DosGroups { get; } = new();
    public ObservableCollection<string> DosGroupNames { get; } = new();
    private string _dosGroup = "total";
    public string DosGroup { get => _dosGroup; set { if (Set(ref _dosGroup, value ?? "total")) DftChanged?.Invoke(); } }
    public (double X, double Y)[] WfProfile { get; private set; } = [];
    public (double X, double Y)[] CddProfile { get; private set; } = [];
    public (double X, double Y)[] CddCumulative { get; private set; } = [];
    public (double X, double Y)[] MdHeight { get; private set; } = [];
    public (double X, double Y)[] MdTemperature { get; private set; } = [];
    private string _resultsNote = "";
    public string ResultsNote { get => _resultsNote; private set => Set(ref _resultsNote, value); }

    public void OpenDftResults()
    {
        if (_resultsSet.Length == 0) ResultsSet = AdsOut.Length > 0 ? AdsOut : DftRoot;
        if (_resultsMain.Length == 0) ResultsMain = DftRoot;
        SetModule(79);
    }

    /// <summary>Every analysis of a set that its finished runs allow; what is missing is said, not guessed.</summary>
    public async Task RefreshResults()
    {
        var notes = new List<string>();
        ResultsFacts.Clear();
        BindRows.Clear(); GeomRows.Clear(); DosGroups.Clear(); DosGroupNames.Clear();
        BindBars = []; BindNames = []; WfProfile = []; CddProfile = []; CddCumulative = []; MdHeight = []; MdTemperature = [];
        string best = "";
        if (Directory.Exists(_resultsSet))
        {
            var b = await Dft("vasp-bind", new JsonObject { ["inputs"] = new JsonArray(_resultsSet) });
            var bind = b?["binding"] as JsonObject;
            if (bind?["error"] is JsonNode e) notes.Add("Binding: " + (string?)e);
            else if (bind != null)
            {
                var rows = (bind["complexes"] as JsonArray ?? []).OfType<JsonObject>().ToList();
                foreach (var r in rows)
                    BindRows.Add(new DftTableRow([((string?)r["complex"] ?? "").Replace("complex_", ""), F(DNum(r["E_bind"]), "0.000"), F(DNum(r["E_bind_kJmol"]), "0.0"), (bool?)r["final"] == true ? "final" : "not final"]));
                BindBars = rows.Select((r, i) => ((double)i, DNum(r["E_bind"]))).ToArray();
                BindNames = rows.Select(r => ((string?)r["complex"] ?? "").Replace("complex_", "")).ToArray();
                best = (string?)bind["best"] ?? "";
                if (best.Length > 0) ResultsFacts.Add(new Row("Best", $"{best.Replace("complex_", "")} · {F(DNum(bind["best_E_bind"]), "0.000")} eV · spread of the best three {F(DNum(bind["spread_best3"]), "0.000")} eV"));
                foreach (var w in (bind["settings_warnings"] as JsonArray ?? [])) notes.Add("Settings: " + (string?)w);
            }
            var cx = Directory.Exists(_resultsSet) ? Directory.GetDirectories(_resultsSet, "complex_*").OrderBy(x => x).ToArray() : [];
            if (cx.Length > 0)
            {
                var g = await Dft("vasp-analyze", new JsonObject { ["inputs"] = new JsonArray(new JsonNode[] { "geom" }.Concat(cx.Select(c => (JsonNode)c)).ToArray()) });
                foreach (var r in (g?["result"]?["complexes"] as JsonArray ?? []).OfType<JsonObject>())
                    GeomRows.Add(new DftTableRow([((string?)r["complex"] ?? "").Replace("complex_", ""), F(DNum(r["h_min"]), "0.00"), F(DNum(r["h_com"]), "0.00"), $"{F(DNum(r["d_anchor"]), "0.00")} {(string?)r["d_anchor_partner"]}",
                        ((int)DNum(r["Hbond_OH_anchor"], 0)).ToString(DftInv), ((int)DNum(r["CH_X_contacts"], 0)).ToString(DftInv), ((int)DNum(r["n_contacts"], 0)).ToString(DftInv), F(DNum(r["tilt_deg"]), "0")]));
            }
            if (best.Length > 0)
            {
                var bd = Path.Combine(_resultsSet, best);
                var st = Path.Combine(bd, "04_static");
                if (File.Exists(Path.Combine(st, "DOSCAR")))
                {
                    var d = await Dft("vasp-analyze", new JsonObject { ["inputs"] = new JsonArray("dos", st) });
                    var res = d?["result"] as JsonObject;
                    var en = (res?["energy"] as JsonArray ?? []).Select(x => DNum(x)).ToArray();
                    foreach (var (k, v) in (res?["groups"] as JsonObject ?? []))
                    {
                        var y = (v as JsonArray ?? []).Select(x => DNum(x)).ToArray();
                        DosGroups[k] = en.Zip(y).Where(p => p.First > -10 && p.First < 5).Select(p => (p.First, p.Second)).ToArray();
                        DosGroupNames.Add(k);
                    }
                    if (res != null) ResultsFacts.Add(new Row("DOS at E_F", $"{F(DNum(res["DOS_at_EF_total"]), "0.00")} states/eV" + (res["mol_HOMO_like"] is JsonNode h ? $" · molecule HOMO-like {F(DNum(h), "0.00")} eV" : "") + (res["mol_LUMO_like"] is JsonNode l ? $" · LUMO-like {F(DNum(l), "0.00")} eV" : "")));
                }
                else notes.Add("DOS: no DOSCAR in the best complex's 04_static yet");
                var chg = Path.Combine(bd, "charge", "04_static");
                var slabSt = Path.Combine(_resultsSet, "slab", "04_static");
                if (File.Exists(Path.Combine(chg, "LOCPOT")))
                {
                    var w = await Dft("vasp-analyze", new JsonObject { ["inputs"] = new JsonArray("wf", chg), ["reference"] = slabSt });
                    var res = w?["result"] as JsonObject;
                    WfProfile = (res?["profile"] as JsonArray ?? []).OfType<JsonArray>().Select(p => (DNum(p[0]), DNum(p[1]))).ToArray();
                    if (res != null) ResultsFacts.Add(new Row("Work function", $"φ top {F(DNum(res["phi_top"]), "0.00")} eV · bottom {F(DNum(res["phi_bottom"]), "0.00")} eV" + (res["delta_phi"] is JsonNode dp ? $" · Δφ {F(DNum(dp), "+0.00;−0.00")} eV" : "") + $" · plateau flatness {F(DNum(res["plateau_flatness_top"]), "0.000")} eV"));
                }
                else notes.Add("Work function: no LOCPOT in <best>/charge/04_static (Runs › Derived set › charge)");
                if (File.Exists(Path.Combine(bd, "cdd_slab", "04_static", "CHGCAR")))
                {
                    var c = await Dft("vasp-analyze", new JsonObject { ["inputs"] = new JsonArray("cdd", bd) });
                    var res = c?["result"] as JsonObject;
                    var prof = (res?["profile"] as JsonArray ?? []).OfType<JsonArray>().ToList();
                    CddProfile = prof.Select(p => (DNum(p[0]), DNum(p[1]))).ToArray();
                    CddCumulative = prof.Select(p => (DNum(p[0]), DNum(p[2]))).ToArray();
                    if (res != null) ResultsFacts.Add(new Row("Charge transfer (CDD)", $"ΔQ at the mid-plane {F(DNum(res["dQ_at_mid"]), "+0.000;−0.000")} e (> 0: the slab side gained) · CHGCAR_diff written"));
                }
                else notes.Add("CDD: run the cdd set of the best complex first");
                var acf = Path.Combine(chg, "ACF.dat");
                if (File.Exists(acf))
                {
                    var q = await Dft("vasp-analyze", new JsonObject { ["inputs"] = new JsonArray("bader", Path.Combine(_resultsSet, "slab", "POSCAR"), Path.Combine(chg, "POSCAR.in"), acf, Path.Combine(bd, "POTCAR")) });
                    if (q?["result"] is JsonObject res) ResultsFacts.Add(new Row("Bader", $"molecule {F(DNum(res["molecule_gained"]), "+0.000;−0.000")} e · slab {F(DNum(res["slab_gained"]), "+0.000;−0.000")} e · total {F(DNum(res["total"]), "+0.000;−0.000")} e"));
                }
                else notes.Add("Bader: no ACF.dat in <best>/charge/04_static (run bader on its CHGCAR there)");
                if (Directory.Exists(Path.Combine(bd, "freq")))
                {
                    var f = await Dft("vasp-analyze", new JsonObject { ["inputs"] = new JsonArray("freq", bd) });
                    if (f?["result"] is JsonObject res && res["shift"] is JsonNode sh) ResultsFacts.Add(new Row("ν(C≡N)", $"{F(DNum(res["nu_CN_adsorbed"]), "0")} cm⁻¹ adsorbed · {F(DNum(res["nu_CN_free"]), "0")} free · shift {F(DNum(sh), "+0;−0")} cm⁻¹"));
                    else if (f?["result"]?["error"] is JsonNode fe) notes.Add("Frequencies: " + (string?)fe);
                }
            }
        }
        else notes.Add("Choose the adsorption set's folder");
        if (Directory.Exists(_resultsAimd) && Directory.GetDirectories(_resultsAimd, "seg_*").Length > 0)
        {
            var m = await Dft("vasp-analyze", new JsonObject { ["inputs"] = new JsonArray("md", _resultsAimd) });
            if (m?["result"] is JsonObject res)
            {
                var ts = (res["timeseries"] as JsonArray ?? []).OfType<JsonArray>().ToList();
                MdHeight = ts.Select(r => (DNum(r[0]), DNum(r[4]))).ToArray();
                MdTemperature = ts.Select(r => (DNum(r[0]), DNum(r[1]))).ToArray();
                ResultsFacts.Add(new Row("AIMD", $"{(string?)res["verdict"]} · T {F(DNum(res["T_mean"]), "0")} ± {F(DNum(res["T_std"]), "0")} K (expected ± {F(DNum(res["T_std_expected"]), "0")}) · drift {F(DNum(res["drift_meV_atom_ps"]), "+0.000;−0.000")} meV/atom/ps · height {F(DNum(res["height_mean"]), "0.00")} Å"));
            }
        }
        DftSummaryRows.Clear();
        if (Directory.Exists(Path.Combine(_resultsMain, "structures")))
        {
            var s = await Dft("vasp-analyze", new JsonObject { ["inputs"] = new JsonArray("summary", _resultsMain) });
            foreach (var r in (s?["result"]?["rows"] as JsonArray ?? []).OfType<JsonObject>())
                DftSummaryRows.Add(new DftTableRow([(string?)r["case"] ?? "", F(DNum(r["a"]), "0.0000"), F(DNum(r["E_bind"]), "0.000"), ((string?)r["best_complex"] ?? "–").Replace("complex_", ""),
                    F(DNum(r["bader_molecule"]), "+0.000;−0.000"), F(DNum(r["delta_phi"]), "+0.00;−0.00"), F(DNum(r["nu_shift"]), "+0;−0"), (string?)r["aimd_verdict"] ?? "–"]));
        }
        ResultsNote = string.Join("\n", notes);
        DftStatus = $"Results of {Path.GetFileName(_resultsSet.TrimEnd('/', '\\'))}: {BindRows.Count} complexes with energies";
        if (DosGroupNames.Count > 0 && !DosGroupNames.Contains(_dosGroup)) DosGroup = DosGroupNames[0];
        DftChanged?.Invoke();
    }

    /// <summary>A derived set made from the selected run (charge, cdd, freq) or an AIMD set from its CONTCAR.</summary>
    public async Task MakeDerived(string kind)
    {
        if (_dftRun == null) { DftStatus = "Choose a run"; return; }
        var a = new JsonObject { ["inputs"] = new JsonArray(kind, _dftRun.Path), ["profile"] = DftProfile };
        if (PotcarDir.Length > 0) a["pp_dir"] = PotcarDir;
        if (kind == "aimd")
        {
            var cont = Path.Combine(_dftRun.Path, "03_relax", "CONTCAR");
            if (!File.Exists(cont)) { DftStatus = "AIMD starts from 03_relax/CONTCAR: not there yet"; return; }
            a["inputs"] = new JsonArray("aimd", cont);
            a["out"] = Path.Combine(Path.GetDirectoryName(_dftRun.Path)!, "aimd_" + Path.GetFileName(_dftRun.Path));
        }
        var r = await Dft("vasp-derived", a);
        if (r != null) DftStatus = $"{kind} set: {r["made"]?.ToJsonString()}";
    }
}
