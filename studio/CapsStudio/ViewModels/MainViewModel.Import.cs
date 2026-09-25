using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>Import (design/boards/ImportDialog): a structure file without a topology (XYZ, extended XYZ, PDB, CIF) opened
/// with the choices it leaves open — where bonds come from and the tolerance, bond orders and aromaticity, molecules by
/// connectivity, unwrapping, the cell — previewed on the first frame and a ten-carbon fragment; then a force field,
/// charges and the validation checks.</summary>
public sealed partial class MainViewModel
{
    /// <summary>Files whose bonds, molecules and cell the import dialog settles (the others open through Open file).</summary>
    public static bool NeedsImport(string path) =>
        Path.GetExtension(path).ToLowerInvariant() is ".xyz" or ".extxyz" or ".pdb" or ".ent" or ".cif";

    private bool _importOpen;
    public bool ImportOpen { get => _importOpen; set => Set(ref _importOpen, value); }
    private string _importPath = "";
    public string ImportTitle => "Import " + Path.GetFileName(_importPath);
    public ObservableCollection<PreviewLine> ImportHead { get; } = new();
    public ObservableCollection<string> ImportNotes { get; } = new();
    public event Action? ImportFragmentChanged;
    public CapsDocument? ImportFragmentDoc { get; private set; }

    private string _impFormat = "", _impBondsChip = "", _impAtoms = "", _impCell = "", _impUnits = "", _impFragText = "", _impSummary = "", _impError = "";
    private bool _impBondsWarn, _impBusy;
    public string ImportFormatChip { get => _impFormat; private set => Set(ref _impFormat, value); }
    public string ImportBondsChip { get => _impBondsChip; private set => Set(ref _impBondsChip, value); }
    public bool ImportBondsWarn { get => _impBondsWarn; private set => Set(ref _impBondsWarn, value); }
    public string ImportAtoms { get => _impAtoms; private set => Set(ref _impAtoms, value); }
    public string ImportCell { get => _impCell; private set => Set(ref _impCell, value); }
    public string ImportUnits { get => _impUnits; private set => Set(ref _impUnits, value); }
    public string ImportFragmentText { get => _impFragText; private set => Set(ref _impFragText, value); }
    public string ImportSummary { get => _impSummary; private set => Set(ref _impSummary, value); }
    public string ImportError { get => _impError; private set { if (Set(ref _impError, value)) Raise(nameof(ImportHasError)); } }
    public bool ImportHasError => _impError.Length > 0;
    public bool ImportBusy { get => _impBusy; private set { if (Set(ref _impBusy, value)) Raise(nameof(ImportIdle)); } }
    public bool ImportIdle => !_impBusy;

    // choices
    private int _impBonds;
    private decimal _impTol = 0.45m;
    private bool _impOrders = true, _impSplit = true, _impUnwrap = true, _impCellOn = true, _impChecks = true;
    public int ImportBondMode { get => _impBonds; set { if (Set(ref _impBonds, value)) { Raise(nameof(ImportPerceives)); RefreshImport(); } } }
    public bool ImportPerceives => _impBonds == 0;
    public decimal ImportTolerance { get => _impTol; set { if (Set(ref _impTol, Math.Clamp(value, 0m, 1.5m))) RefreshImport(); } }
    public bool ImportOrders { get => _impOrders; set { if (Set(ref _impOrders, value)) RefreshImport(); } }
    public bool ImportSplit { get => _impSplit; set { if (Set(ref _impSplit, value)) RefreshImport(); } }
    public bool ImportUnwrap { get => _impUnwrap; set { if (Set(ref _impUnwrap, value)) RefreshImport(); } }
    public bool ImportUseCell { get => _impCellOn; set { if (Set(ref _impCellOn, value)) RefreshImport(); } }
    public bool ImportChecks { get => _impChecks; set => Set(ref _impChecks, value); }
    public static readonly string[] ImportRadii = ["Covalent · Cordero 2008"];

    /// <summary>"Keep as read", then the force-field library (typing with it after import).</summary>
    public IEnumerable<string> ImportForceFields => new[] { "None · keep as read" }.Concat(Field.Library.Select(e => "Type with " + e.Name));
    private int _impFf, _impCharges;
    public int ImportForceField { get => _impFf; set { if (Set(ref _impFf, value)) Raise(nameof(ImportChargesEnabled)); } }
    public int ImportCharges { get => _impCharges; set => Set(ref _impCharges, value); }
    public bool ImportChargesEnabled => _impFf > 0;
    public static string[] ImportChargeModes => FieldViewModel.ChargeModes;

    private string ImportOptionsJson() => JsonSerializer.Serialize(new
    {
        bonds = _impBonds switch { 1 => "file", 2 => "none", _ => "perceive" },
        tolerance = (double)_impTol,
        bond_orders = _impOrders,
        split = _impSplit,
        unwrap = _impUnwrap,
        use_cell = _impCellOn,
    });

    public void ShowImport(string path)
    {
        _importPath = path;
        Raise(nameof(ImportTitle));
        Raise(nameof(ImportForceFields));
        _impBonds = 0;
        Raise(nameof(ImportBondMode));
        Raise(nameof(ImportPerceives));
        _firstImportPreview = true;
        ImportError = "";
        ImportOpen = true;
        RefreshImport();
    }

    private int _impGen;
    private bool _firstImportPreview;
    private Task? _importTask;
    /// <summary>The choices applied to frame 0 (off the UI thread; the newest request wins).</summary>
    private void RefreshImport()
    {
        if (!_importOpen || _importPath.Length == 0) return;
        var gen = ++_impGen;
        var (path, opts) = (_importPath, ImportOptionsJson());
        ImportBusy = true;
        _importTask = Task.Run(() =>
        {
            string? json = null, error = null;
            CapsDocument? frag = null;
            try
            {
                json = CapsDocument.ImportPreview(path, opts);
                if (JsonNode.Parse(json)?["ok"]?.GetValue<bool>() == true) frag = CapsDocument.ImportFragment(path, opts);
            }
            catch (Exception e) { error = e.Message; }
            Avalonia.Threading.Dispatcher.UIThread.Post(() =>
            {
                if (gen != _impGen) { frag?.Dispose(); return; }
                ImportBusy = false;
                FillImport(json, error, frag);
            });
        });
    }

    /// <summary>Waits for the preview in flight (self-test and screenshots).</summary>
    public void WaitImport()
    {
        for (var i = 0; i < 400 && (_impBusy || _importTask is { IsCompleted: false }); i++) { Avalonia.Threading.Dispatcher.UIThread.RunJobs(); Thread.Sleep(10); }
        Avalonia.Threading.Dispatcher.UIThread.RunJobs();
    }

    private void FillImport(string? json, string? error, CapsDocument? frag)
    {
        var inv = CultureInfo.InvariantCulture;
        var j = json != null ? JsonNode.Parse(json) : null;
        if (j == null || j["ok"]?.GetValue<bool>() != true)
        {
            ImportError = error ?? j?["error"]?.GetValue<string>() ?? "cannot read the file";
            frag?.Dispose();
            return;
        }
        ImportError = "";
        var inFile = (int)(j["bonds_in_file"]?.GetValue<double>() ?? 0);
        if (_firstImportPreview)
        {
            _firstImportPreview = false;
            if (inFile > 0) { _impBonds = 1; Raise(nameof(ImportBondMode)); Raise(nameof(ImportPerceives)); frag?.Dispose(); RefreshImport(); return; }
        }
        var fmt = j["format"]?.GetValue<string>() ?? "";
        var head = (JsonArray)j["head"]!;
        var extended = fmt == "xyz" && head.Count > 1 && (head[1]?.GetValue<string>() ?? "").Contains("Lattice=");
        ImportFormatChip = (extended ? "Extended XYZ" : j["format_name"]?.GetValue<string>() ?? fmt) + " detected";
        ImportBondsWarn = inFile == 0;
        ImportBondsChip = inFile == 0 ? "no bonds in file" : $"{inFile:N0} bonds in file";
        ImportHead.Clear();
        var k = 0;
        foreach (var l in head.Take(8)) ImportHead.Add(new PreviewLine((++k).ToString(inv), l?.GetValue<string>() ?? ""));
        if (head.Count > 8) ImportHead.Add(new PreviewLine("", "…"));
        var atoms = j["atoms"]!.GetValue<double>();
        ImportAtoms = atoms.ToString("N0", inv);
        var cell = j["cell"]?.GetValue<string>() ?? "";
        ImportCell = cell.Length > 0 ? cell : _impCellOn ? "none in the file" : "dropped";
        var units = j["units"]?.GetValue<string>() ?? "";
        ImportUnits = units.Length > 0 ? units : fmt == "xyz" ? "Å (XYZ convention)" : "Å";
        var bonds = j["bonds"]!.GetValue<double>();
        var mols = j["molecules"]!.GetValue<double>();
        var parts = new List<string> { $"{bonds:N0} bonds", $"{mols:N0} molecule{(mols == 1 ? "" : "s")}" };
        foreach (var (key, name) in new[] { ("double", "double"), ("triple", "triple"), ("aromatic", "aromatic") })
            if (j[key]?.GetValue<double>() is > 0 and var n) parts.Add($"{n:N0} {name}");
        ImportSummary = string.Join(" · ", parts);
        ImportFragmentText = $"{j["fragment_bonds"]!.GetValue<double>():N0} bonds {(_impBonds == 0 ? "perceived" : _impBonds == 1 ? "from the file" : "")} in the preview fragment · {j["fragment_heavy"]!.GetValue<double>():N0} heavy atoms";
        ImportNotes.Clear();
        foreach (var n in (JsonArray)j["notes"]!) ImportNotes.Add(n?.GetValue<string>() ?? "");
        var old = ImportFragmentDoc;
        ImportFragmentDoc = frag;
        ImportFragmentChanged?.Invoke();
        old?.Dispose();
    }

    public void CloseImport()
    {
        ImportOpen = false;
        ++_impGen;
        var old = ImportFragmentDoc;
        ImportFragmentDoc = null;
        ImportFragmentChanged?.Invoke();
        old?.Dispose();
    }

    /// <summary>Opens the whole file with the choices as a new document, then the force field and the checks.</summary>
    public async Task ConfirmImport()
    {
        if (_importPath.Length == 0 || _impBusy) return;
        var (path, opts) = (_importPath, ImportOptionsJson());
        ImportBusy = true;
        CapsDocument doc;
        try { doc = await Task.Run(() => CapsDocument.Import(path, null, opts)); }
        catch (Exception e) { ImportBusy = false; ImportError = e.Message; return; }
        ImportBusy = false;
        CloseImport();
        if (!IsStudio) SetModule(8);
        Show(doc, Path.GetFileName(path));
        if (_doc != doc) return;
        Remember(path, null);
        var args = new List<string> { PyStr(Path.GetFullPath(path)) };
        if (_impBonds != 0) args.Add($"bonds=\"{(_impBonds == 1 ? "file" : "none")}\"");
        if (_impTol != 0.45m) args.Add("tolerance=" + _impTol.ToString(CultureInfo.InvariantCulture));
        if (!_impOrders) args.Add("bond_orders=False");
        if (!_impSplit) args.Add("split=False");
        if (!_impUnwrap) args.Add("unwrap=False");
        if (!_impCellOn) args.Add("use_cell=False");
        Record($"doc = caps.import_file({string.Join(", ", args)})");
        if (_impFf > 0 && _impFf - 1 < Field.Library.Count)
        {
            Field.FfIndex = _impFf - 1;
            Field.ChargeMode = _impCharges;
            await Field.Assign();
        }
        if (_impChecks) OpenChecks();
        Status = $"Imported {Path.GetFileName(path)} · {ImportSummary}";
        var sum = doc.Summary();
        if (_impBonds == 0 && sum.BondsFromFile == 0 && sum.Bonds > 0)
            Notify(new Notice
            {
                Key = "import.bonds", Severity = "check", Icon = "link", Title = "Bonds were inferred",
                Body = $"{Path.GetFileName(path)} had no bond records. {sum.Bonds:N0} bonds were perceived from covalent radii (tolerance {_impTol.ToString(CultureInfo.InvariantCulture)} Å); review them before typing.",
                Primary = "Review bonds", OnPrimary = OpenChecks,
                Secondary = "Undo import", OnSecondary = CloseDocument,
            });
    }
}
