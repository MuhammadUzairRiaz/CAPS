using System.Linq;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A force field of the CAPS library (data/forcefields/catalogue.json).</summary>
public sealed record FfEntry(string Id, string Name, string Version, string Status, string File, bool AutoTyping, string Key = "")
{
    public string Label => (Key.Length > 0 ? $"{Name}  [{Key}]" : Name) + (AutoTyping ? "" : " · types from the file");
    public override string ToString() => Label;
}

/// <summary>Type names as shown: moltemplate's OPLS-AA names carry their bonded classes (135_bCT_aCT_dCT_iCT), shown as
/// the number and its class (135 · CT); the full name stays in the files and tooltips.</summary>
public static class FieldNames
{
    private static readonly System.Text.RegularExpressions.Regex Opls = new(@"^(\d+)_b([^_]+)_a[^_]+_d[^_]+_i[^_]+$");
    public static string Short(string name) => Opls.Match(name) is { Success: true } m ? $"{m.Groups[1].Value} · {m.Groups[2].Value}" : name;
}

/// <summary>One atom of the typing report.</summary>
public sealed class FieldAtomRow
{
    public int Index { get; init; }            // 0-based
    public string Atom { get; init; } = "";    // element + number, as in the design (C1, H14)
    public string Element { get; init; } = "";
    public string Type { get; init; } = "";
    public string Charge { get; init; } = "";
    public string Rule { get; init; } = "";
    public string Source { get; init; } = "";
    public string Description { get; init; } = "";
    public int Priority { get; init; }
    public string Why { get; init; } = "";
    public string[] Candidates { get; init; } = [];
    public bool Overridden { get; init; }
    public bool Untyped => Type.Length == 0;
    public string TypeShown => Untyped ? "?" : FieldNames.Short(Type);
    public IBrush TypeBrush => Untyped ? Bad : Overridden ? Sel : Text;
    public IBrush RuleBrush => Untyped ? Bad : Dim;
    internal static IBrush Bad => Tokens.Brush("ErrB");
    internal static IBrush Sel => Tokens.Brush("SelB");
    internal static IBrush Text => Tokens.Brush("TextB");
    internal static IBrush Dim => Tokens.Brush("MutedB");
}

/// <summary>A type present in the structure, with the viewer's colour.</summary>
public sealed record FieldTypeSwatch(string Name, int Count, IBrush Colour, string Description)
{
    public string Label => $"{FieldNames.Short(Name)} · {Count}";
}

/// <summary>A term or atom that keeps the structure from running.</summary>
public sealed class FieldMissingRow
{
    public string Kind { get; init; } = "";      // bond, angle, dihedral, improper, pair, "bond increment", "atom type"
    public string Title { get; init; } = "";
    public string Detail { get; init; } = "";
    public string Types { get; init; } = "";     // lookup names, for a rule entered by hand
    public int Atom { get; init; } = -1;         // untyped atom, 0-based
    public bool CanEnter => Kind is "bond" or "angle" or "dihedral" or "improper" or "pair";
}

/// <summary>CAPS Field: choose a force field, see why every atom has its type, override types, and complete missing
/// parameters (imported, or entered by hand and flagged estimated).</summary>
public sealed partial class FieldViewModel : ObservableObject
{
    private readonly Func<CapsDocument?> _doc;
    private readonly Action<string> _status;
    private readonly Action _changed;   // types or charges in the document changed: re-render, refresh the inspector
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public FieldViewModel(Func<CapsDocument?> doc, Action<string> status, Action changed)
    {
        _doc = doc;
        _status = status;
        _changed = changed;
        LoadCatalogue();
    }

    // ---------------------------------------------------------------- library
    public ObservableCollection<FfEntry> Library { get; } = new();
    public string LibraryNote { get; private set; } = "";
    private int _ffIndex = -1;
    public int FfIndex { get => _ffIndex; set { if (Set(ref _ffIndex, value)) Raise(nameof(FfNote)); } }
    public FfEntry? Selected => _ffIndex >= 0 && _ffIndex < Library.Count ? Library[_ffIndex] : null;
    public string FfNote => Selected is { } e ? $"{e.Version} · {e.Status}" + (e.AutoTyping ? " · automatic typing" : " · types must be the atom names in the file") : "";
    /// <summary>Index 0 is automatic (core mode 4); the others are core modes 0–3 in order.</summary>
    public static readonly string[] ChargeModes = ["Automatic (force field, else Gasteiger)", "From the force field", "Gasteiger–Marsili", "Keep the file's charges", "QEq (every element)"];
    private static int CoreCharges(int ui) => ui == 0 ? 4 : ui - 1;
    private int _chargeMode;
    public int ChargeMode { get => _chargeMode; set => Set(ref _chargeMode, value); }

    /// <summary>Kept for callers that rank force fields: the list's own order (catalogue "list.order") now decides.</summary>
    internal static readonly string[] PolymerFirst = ["pcff-frc", "compass-frc", "opls2005", "oplsaa2024-moltemplate"];

    private void LoadCatalogue()
    {
        var dir = FindLibrary();
        if (dir == null) { LibraryNote = "The force-field library (data/forcefields) was not found next to the Studio."; return; }
        try
        {
            using var js = JsonDocument.Parse(File.ReadAllText(Path.Combine(dir, "catalogue.json")));
            var list = new List<(int Order, FfEntry E)>();
            foreach (var e in js.RootElement.GetProperty("forcefields").EnumerateArray())
            {
                // one entry per force field: the catalogue marks the ones listed ("list": label, key, order); the others
                // (other copies of the same force field) stay loadable by id for recipes and scripts
                if (!e.TryGetProperty("list", out var l) || l.ValueKind != JsonValueKind.Object) continue;
                if (!e.TryGetProperty("file", out var f) || f.ValueKind != JsonValueKind.String) continue;
                var typing = e.TryGetProperty("typing", out var t) && t.ValueKind == JsonValueKind.Object;
                list.Add((l.TryGetProperty("order", out var o) ? o.GetInt32() : 999,
                          new FfEntry(Str(e, "id"), Str(l, "label"), Str(e, "version"), Str(e, "status"), Path.Combine(dir, f.GetString()!), typing, Str(l, "key"))));
            }
            foreach (var (_, x) in list.OrderBy(x => x.Order)) Library.Add(x);
            // UFF (built into the core): every element, typed from bonds, hybridisation and oxidation state — before the inorganic ones
            var inorganic = Library.ToList().FindIndex(x => x.Id.StartsWith("inorganic", StringComparison.Ordinal));
            Library.Insert(inorganic < 0 ? Library.Count : inorganic, new FfEntry("uff", "UFF · every element", "1992", "validated", "uff", true, "uff"));
            LibraryNote = $"{Library.Count} force fields · {Library.Count(x => x.AutoTyping)} with automatic typing";
            FfIndex = Library.Count > 0 ? 0 : -1;
        }
        catch (Exception ex) { LibraryNote = "Cannot read the force-field catalogue: " + ex.Message; }
    }

    // installed layouts and the development tree (see Paths)
    private static string? FindLibrary() => Paths.ForceFields;

    private static string Str(JsonElement e, string k) => e.TryGetProperty(k, out var v) && v.ValueKind == JsonValueKind.String ? v.GetString()! : "";

    // ---------------------------------------------------------------- state
    private bool _assigned, _complete, _working;
    public bool Assigned { get => _assigned; private set { if (Set(ref _assigned, value)) { Raise(nameof(NotAssigned)); Raise(nameof(ShowReport)); Raise(nameof(Incomplete)); } } }
    public bool NotAssigned => !_assigned;
    public bool ShowReport => _assigned;
    public bool Complete { get => _complete; private set { if (Set(ref _complete, value)) Raise(nameof(Incomplete)); } }
    public bool Incomplete => _assigned && !_complete;
    public bool Working { get => _working; private set { if (Set(ref _working, value)) Raise(nameof(NotWorking)); } }
    public bool NotWorking => !_working;

    private string _ffName = "", _typedText = "", _untypedText = "", _missingText = "", _estimatedText = "", _chargeText = "", _footerTyper = "",
                   _footerDb = "", _energyText = "", _log = "Choose a force field and assign it. CAPS types every atom from its rules, sets the charges and looks up every parameter; nothing is guessed.";
    private bool _chargeOdd;
    public string ForceFieldName { get => _ffName; private set => Set(ref _ffName, value); }
    public string TypedText { get => _typedText; private set => Set(ref _typedText, value); }
    public string UntypedText { get => _untypedText; private set => Set(ref _untypedText, value); }
    public string MissingText { get => _missingText; private set => Set(ref _missingText, value); }
    public string EstimatedText { get => _estimatedText; private set => Set(ref _estimatedText, value); }
    public string ChargeText { get => _chargeText; private set => Set(ref _chargeText, value); }
    public bool ChargeOdd { get => _chargeOdd; private set => Set(ref _chargeOdd, value); }
    public string FooterTyper { get => _footerTyper; private set => Set(ref _footerTyper, value); }
    public string FooterDb { get => _footerDb; private set => Set(ref _footerDb, value); }
    public string EnergyText { get => _energyText; private set => Set(ref _energyText, value); }
    public string Log { get => _log; private set => Set(ref _log, value); }
    private bool _hasUntyped, _hasMissing, _hasEstimated;
    public bool HasUntyped { get => _hasUntyped; private set => Set(ref _hasUntyped, value); }
    public bool HasMissing { get => _hasMissing; private set => Set(ref _hasMissing, value); }
    public bool HasEstimated { get => _hasEstimated; private set => Set(ref _hasEstimated, value); }

    /// <summary>One line for the Relax / Dynamics / Equilibrate inspectors.</summary>
    public string RunLine => !_assigned ? "Force field: built-in GAFF for C and H, UFF for other elements · assign another in Field"
        : _complete ? $"Force field: {_ffName} (from Field) · complete"
        : $"Force field: {_ffName} (from Field) · incomplete: runs are blocked until Field is complete";

    // ---------------------------------------------------------------- atoms table
    private readonly List<FieldAtomRow> _all = new();
    public ObservableCollection<FieldAtomRow> Rows { get; } = new();
    private int _filterMode;   // 0 all, 1 untyped, 2 overridden
    private string _filterText = "";
    public int FilterMode { get => _filterMode; set { if (Set(ref _filterMode, value)) { RaiseFilter(); ApplyFilter(); } } }
    public bool FilterAll { get => _filterMode == 0; set { if (value) FilterMode = 0; } }
    public bool FilterUntyped { get => _filterMode == 1; set { if (value) FilterMode = 1; } }
    public bool FilterOverridden { get => _filterMode == 2; set { if (value) FilterMode = 2; } }
    private void RaiseFilter() { Raise(nameof(FilterAll)); Raise(nameof(FilterUntyped)); Raise(nameof(FilterOverridden)); }
    public string FilterText { get => _filterText; set { if (Set(ref _filterText, value)) ApplyFilter(); } }
    private string _rowsNote = "";
    public string RowsNote { get => _rowsNote; private set => Set(ref _rowsNote, value); }

    private void ApplyFilter()
    {
        var words = _filterText.Split(' ', StringSplitOptions.RemoveEmptyEntries);
        IEnumerable<FieldAtomRow> q = _all;
        if (_filterMode == 1) q = q.Where(r => r.Untyped);
        if (_filterMode == 2) q = q.Where(r => r.Overridden);
        foreach (var w in words)
            q = q.Where(r => r.Atom.Equals(w, StringComparison.OrdinalIgnoreCase) || r.Type.Equals(w, StringComparison.Ordinal) ||
                             r.Element.Equals(w, StringComparison.OrdinalIgnoreCase) || r.Rule.Contains(w, StringComparison.Ordinal) ||
                             r.Charge.StartsWith(w, StringComparison.Ordinal));
        var sel = _selected?.Index;
        Rows.Clear();
        foreach (var r in q) Rows.Add(r);
        RowsNote = Rows.Count == _all.Count ? $"{_all.Count:N0} atoms" : $"{Rows.Count:N0} of {_all.Count:N0} atoms";
        if (sel is int s) Selected0 = Rows.FirstOrDefault(r => r.Index == s);
    }

    // ---------------------------------------------------------------- selected atom: why, override
    private FieldAtomRow? _selected;
    public FieldAtomRow? SelectedRow { get => _selected; set { if (Set(ref _selected, value)) OnSelected(); } }
    private FieldAtomRow? Selected0 { set { _selected = value; Raise(nameof(SelectedRow)); OnSelected(); } }
    public event Action<int>? AtomSelected;   // the viewer highlights it
    private string _whyTitle = "Pick an atom", _whyText = "Select a row (or click an atom in the 3D view) to see the rule that typed it, the rules it also matched, and to set its type by hand.";
    public string WhyTitle { get => _whyTitle; private set => Set(ref _whyTitle, value); }
    public string WhyText { get => _whyText; private set => Set(ref _whyText, value); }
    public ObservableCollection<string> OverrideTypes { get; } = new();
    private string? _overrideType;
    public string? OverrideType { get => _overrideType; set => Set(ref _overrideType, value); }
    public bool HasSelection => _selected != null;
    public bool SelectedOverridden => _selected?.Overridden == true;
    private readonly Dictionary<string, (string El, string Desc)> _fftypes = new();

    private void OnSelected()
    {
        Raise(nameof(HasSelection));
        Raise(nameof(SelectedOverridden));
        OverrideTypes.Clear();
        if (_selected is not { } r) return;
        AtomSelected?.Invoke(r.Index);
        WhyTitle = r.Untyped ? $"Why {r.Atom} has no type" : $"Why {r.Atom} is {r.Type}";
        WhyText = Explain(r);
        foreach (var (name, (el, _)) in _fftypes.OrderBy(x => x.Key, StringComparer.Ordinal))
            if (el == r.Element || el == "" || el == "?") OverrideTypes.Add(name);
        OverrideType = r.Untyped ? OverrideTypes.FirstOrDefault() : r.Type;
    }

    /// <summary>Why an atom has its type, without changing the page's selection (the Studio inspector asks too).</summary>
    public string Explain(int index) => _all.FirstOrDefault(r => r.Index == index) is { } row ? Explain(row) : "";

    private string Explain(FieldAtomRow r)
    {
        var t = new System.Text.StringBuilder();
        if (r.Overridden) t.Append($"Set to {r.Type} by hand (the rules gave {(r.Candidates.Length > 0 ? string.Join(", ", r.Candidates) : "no type")}).");
        else if (r.Untyped) t.Append(r.Candidates.Length > 0 ? $"Rules matched {string.Join(", ", r.Candidates)}, but none of these is a type of the force field." : "No typing rule matches this atom's element and bonding. Set a type by hand, or check the structure (missing hydrogens, wrong bonds).");
        else if (r.Source == "file") t.Append($"The force field has no typing rules: {r.Type} is the atom name in the file.");
        else
        {
            t.Append($"Rule {r.Type}");
            if (r.Priority != 0) t.Append($" (priority {r.Priority})");
            t.Append($" matched {r.Rule}");
            if (r.Description.Length > 0) t.Append($": {r.Description.TrimEnd('.')}");
            t.Append('.');
            var others = r.Candidates.Where(c => c != r.Type).ToArray();
            if (others.Length > 0) t.Append($" Rules for {string.Join(", ", others)} also matched and lost on priority, an override or file order.");
        }
        if (_fftypes.TryGetValue(r.Type, out var d) && d.Desc.Length > 0 && !r.Untyped) t.Append($" Type {r.Type}: {d.Desc.TrimEnd('.')}.");
        if (r.Charge.Length > 0) t.Append($" Charge {r.Charge} e ({ChargeModes[_chargeMode].ToLowerInvariant()}).");
        return t.ToString();
    }

    /// <summary>Select the row of an atom picked in the viewer.</summary>
    public void SelectAtom(int index)
    {
        if (index < 0 || index >= _all.Count) return;
        if (Rows.All(r => r.Index != index)) { _filterText = ""; Raise(nameof(FilterText)); FilterMode = 0; ApplyFilter(); }
        SelectedRow = Rows.FirstOrDefault(r => r.Index == index);
    }

    // ---------------------------------------------------------------- types in 3D, missing terms
    public ObservableCollection<FieldTypeSwatch> Swatches { get; } = new();
    public ObservableCollection<FieldMissingRow> Missing { get; } = new();
    public ObservableCollection<string> Entered { get; } = new();
    public bool HasEntered => Entered.Count > 0;
    public ObservableCollection<string> Notes { get; } = new();
    private Dictionary<string, string> _styles = new();

    private FieldMissingRow? _missingSel;
    public FieldMissingRow? SelectedMissing
    {
        get => _missingSel;
        set
        {
            if (!Set(ref _missingSel, value)) return;
            Raise(nameof(CanEnterRule));
            if (value == null) return;
            if (value.Atom >= 0) { SelectAtom(value.Atom); return; }
            if (!value.CanEnter) return;
            EntryKind = value.Kind;
            EntryTypes = value.Types;
            EntryStyle = _styles.TryGetValue(value.Kind, out var st) ? st : "";
            EntryParams = "";
        }
    }
    public bool CanEnterRule => _missingSel?.CanEnter == true;
    private string _entryKind = "", _entryTypes = "", _entryStyle = "", _entryParams = "";
    public string EntryKind { get => _entryKind; set { if (Set(ref _entryKind, value)) Raise(nameof(EntryHint)); } }
    public string EntryTypes { get => _entryTypes; set => Set(ref _entryTypes, value); }
    public string EntryStyle { get => _entryStyle; set { if (Set(ref _entryStyle, value)) Raise(nameof(EntryHint)); } }
    public string EntryParams { get => _entryParams; set => Set(ref _entryParams, value); }
    public string EntryHint => ParamHint(_entryKind, _entryStyle);

    /// <summary>The parameter order of each style, in CAPS / LAMMPS units (kcal/mol, Å, degrees).</summary>
    private static string ParamHint(string kind, string style) => (kind, style) switch
    {
        ("pair", _) => "ε (kcal/mol)  σ (Å)",
        ("bond", "class2") => "r0 (Å)  K2  K3  K4",
        ("bond", "morse") => "D (kcal/mol)  α (1/Å)  r0 (Å)",
        ("bond", _) => "K (kcal/mol/Å²)  r0 (Å)   ·   E = K (r − r0)²",
        ("angle", "class2" or "quartic") => "θ0 (°)  K2  K3  K4",
        ("angle", "charmm") => "K  θ0 (°)  K_UB  r_UB (Å)",
        ("angle", _) => "K (kcal/mol/rad²)  θ0 (°)   ·   E = K (θ − θ0)²",
        ("dihedral", "opls") => "V1  V2  V3  V4 (kcal/mol)",
        ("dihedral", "class2") => "K1 φ1  K2 φ2  K3 φ3",
        ("dihedral", "charmm") => "K  n  d (°)",
        ("dihedral", _) => "m, then K n d for each of the m terms   ·   E = Σ K [1 + cos(nφ − d)]",
        ("improper", "cvff") => "K  d (±1)  n",
        ("improper", "class2" or "inversion") => "K  χ0 (°)",
        ("improper", "harmonic") => "K  χ0 (°)",
        ("improper", _) => "m, then K n d for each term",
        _ => "",
    };

    // ---------------------------------------------------------------- actions
    private async Task Do(string what, Func<CapsDocument, bool> action)
    {
        var doc = _doc();
        if (doc == null) { Log = "Open or build a structure first."; return; }
        if (doc.LongRunning) { _status("A run is using this structure: assign or change the force field when it finishes"); return; }
        Working = true;
        try
        {
            var complete = await Task.Run(() => action(doc));
            if (!ReferenceEquals(_doc(), doc)) return;   // another structure is open now: this assignment went with the old one
            LoadReport(doc);
            _changed();
            _ = CheckCoverage(auto: true);   // when this force field cannot describe the structure: why, and which can
            _status(complete ? $"{what} · {_ffName}: complete" : $"{what} · {_ffName}: {UntypedText}, {MissingText}");
        }
        catch (ObjectDisposedException) { }   // the structure was closed while it was being typed
        catch (Exception ex) { Log = what + " failed: " + ex.Message; _status(Log); }
        finally { Working = false; }
    }

    /// <summary>The macro recorder's hook: each assignment as a line of Python.</summary>
    public Action<string>? Recorder { get; set; }

    /// <summary>UFF for every element (the fallback a blocked run offers).</summary>
    public async void UseUff()
    {
        var k = Library.ToList().FindIndex(e => e.Id == "uff");
        if (k < 0) return;
        FfIndex = k;
        await Assign();
    }

    public Task Assign()
    {
        if (Selected is not { } e) { Log = "Choose a force field."; return Task.CompletedTask; }
        var mode = _chargeMode;
        Recorder?.Invoke($"doc.field.assign(\"{e.File.Replace("\\", "/")}\", charges=\"{(mode switch { 1 => "forcefield", 2 => "gasteiger", 3 => "keep", 4 => "qeq", _ => "auto" })}\")");
        return Do("Assigned", d => d.FieldAssign(e.File, null, CoreCharges(mode)));
    }

    public Task ApplyOverride()
    {
        if (_selected is not { } r || string.IsNullOrEmpty(_overrideType)) return Task.CompletedTask;
        var (i, t) = (r.Index, _overrideType);
        return Do($"{r.Atom} set to {t}", d => d.FieldOverride(i, t));
    }

    public Task ResetOverride()
    {
        if (_selected is not { } r) return Task.CompletedTask;
        var i = r.Index;
        return Do($"{r.Atom} back to its rule", d => d.FieldOverride(i, null));
    }

    public Task AddRule()
    {
        var (k, t, s, p) = (_entryKind, _entryTypes, _entryStyle, _entryParams);
        if (k.Length == 0 || p.Trim().Length == 0) { Log = "Choose a missing term and give its parameters."; return Task.CompletedTask; }
        return Do($"Entered {k} {t} (estimated)", d => d.FieldAddRule(k, t, s == _styles.GetValueOrDefault(k) ? "" : s, p));
    }

    public Task Import(string path) => Do($"Imported {Path.GetFileName(path)}", d => d.FieldImport(path));
    /// <summary>Fill gaps: another force field's terms only where this one has none (OPLS-AA 2024 from OPLS 2005 …),
    /// each borrowed term listed in the report.</summary>
    public Task FillGaps(string path) => Do($"Missing terms filled from {Path.GetFileName(path)} (only where the force field has none)", d => d.FieldFillGaps(path));
    public Task RemoveRules() => Do("Removed imported and entered parameters", d => d.FieldRemoveRules());

    public async Task Clear()
    {
        var doc = _doc();
        if (doc == null) return;
        await Task.Run(doc.FieldClear);
        Reset();
        _changed();
        _status("Force-field assignment cleared; the file's own types and charges are back");
    }

    public void SaveTypes(string path)
    {
        _doc()?.FieldTypesFile(path);
        _status($"Saved the types to {path} (caps ff apply --types)");
    }

    /// <summary>A new structure: the assignment belongs to the old one.</summary>
    public void Reset()
    {
        Assigned = false;
        Complete = false;
        _all.Clear();
        Rows.Clear();
        Swatches.Clear();
        Missing.Clear();
        Entered.Clear();
        Raise(nameof(HasEntered));
        Notes.Clear();
        _fftypes.Clear();
        SelectedRow = null;
        ForceFieldName = "";
        EnergyText = "";
        Log = "Choose a force field and assign it. CAPS types every atom from its rules, sets the charges and looks up every parameter; nothing is guessed.";
        Raise(nameof(RunLine));
    }

    /// <summary>Reads the core's report after any change (also after Relax / Dynamics, which keep the assignment).</summary>
    public void LoadReport(CapsDocument doc)
    {
        var json = doc.FieldReport();
        if (json.Length == 0) { Reset(); return; }
        using var js = JsonDocument.Parse(json);
        var r = js.RootElement;
        ForceFieldName = Str(r, "forcefield");
        var typed = r.GetProperty("typed").GetDouble();
        var untyped = r.GetProperty("untyped").GetDouble();
        var missing = r.GetProperty("missing").EnumerateArray().Select(x => x.GetString()!).ToList();
        var estimated = r.GetProperty("estimated").GetDouble();
        var imported = r.GetProperty("imported").GetDouble();
        TypedText = $"{typed.ToString("N0", Inv)} typed";
        UntypedText = $"{untyped:N0} untyped";
        MissingText = missing.Count == 1 ? "1 missing term" : $"{missing.Count:N0} missing terms";
        var filled = r.TryGetProperty("filled", out var fl) ? fl.GetDouble() : 0;
        EstimatedText = $"{estimated:N0} estimated" + (imported > 0 ? $" · {imported:N0} imported" : "") + (filled > 0 ? $" · {filled:N0} filled from another force field" : "");
        HasUntyped = untyped > 0;
        HasMissing = missing.Count > 0;
        HasEstimated = estimated > 0;
        if (r.GetProperty("has_charges").GetBoolean())
        {
            var q = r.GetProperty("net_charge").GetDouble();
            ChargeOdd = Math.Abs(q - Math.Round(q)) > 0.01;
            ChargeText = $"net charge {q.ToString("+0.000;−0.000;0.000", Inv)} e" + (ChargeOdd ? " · not a whole charge" : "");
        }
        else { ChargeText = "charges after typing"; ChargeOdd = false; }
        _fftypes.Clear();
        foreach (var t in r.GetProperty("fftypes").EnumerateArray()) _fftypes[Str(t, "name")] = (Str(t, "el"), Str(t, "desc"));
        _styles = r.GetProperty("styles").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString() ?? "");

        _all.Clear();
        var counter = new Dictionary<string, int>();
        foreach (var a in r.GetProperty("atoms").EnumerateArray())
        {
            var el = Str(a, "el");
            counter[el] = counter.GetValueOrDefault(el) + 1;
            _all.Add(new FieldAtomRow
            {
                Index = a.GetProperty("i").GetInt32() - 1,
                Atom = el + counter[el].ToString(Inv),
                Element = el,
                Type = Str(a, "type"),
                Charge = a.TryGetProperty("q", out var q) ? q.GetDouble().ToString("+0.0000;−0.0000;0.0000", Inv) : "",
                Rule = Str(a, "rule"),
                Source = Str(a, "src"),
                Description = Str(a, "desc"),
                Priority = a.TryGetProperty("prio", out var p) ? p.GetInt32() : 0,
                Why = Str(a, "why"),
                Candidates = a.GetProperty("cands").EnumerateArray().Select(x => x.GetString()!).ToArray(),
                Overridden = a.GetProperty("ov").GetBoolean(),
            });
        }
        ApplyFilter();

        Swatches.Clear();
        foreach (var u in r.GetProperty("used").EnumerateArray())
            Swatches.Add(new FieldTypeSwatch(Str(u, "name"), u.GetProperty("count").GetInt32(), new SolidColorBrush(Color.Parse(Str(u, "colour"))), Str(u, "desc")));

        Missing.Clear();
        foreach (var row in _all.Where(x => x.Untyped).Take(200))
            Missing.Add(new FieldMissingRow { Kind = "atom type", Title = $"Atom type · {row.Atom}", Detail = "No typing rule matched: set its type by hand", Atom = row.Index });
        foreach (var m in missing.Take(500))
        {
            if (m.StartsWith("charges:", StringComparison.Ordinal))   // the physics check, not a parameter to enter
            {
                Missing.Add(new FieldMissingRow { Kind = "charges", Title = "Charges do not balance", Detail = m["charges:".Length..].Trim() });
                continue;
            }
            // "bond c3(c) ca": the kind, then each type with its lookup name in parentheses when different
            var kind = m.StartsWith("bond increment", StringComparison.Ordinal) ? "bond increment" : m.Split(' ')[0];
            var rest = m[kind.Length..].Trim();
            var names = rest.Split(' ', StringSplitOptions.RemoveEmptyEntries)
                            .Select(w => w.Contains('(') && w.EndsWith(')') ? w[(w.IndexOf('(') + 1)..^1] : w).ToArray();
            Missing.Add(new FieldMissingRow
            {
                Kind = kind,
                Title = $"{char.ToUpperInvariant(kind[0])}{kind[1..]} · {string.Join(" – ", rest.Split(' ', StringSplitOptions.RemoveEmptyEntries).Select(FieldNames.Short))}",
                Detail = kind == "bond increment" ? $"Not in {ForceFieldName}: charges need it (or use Gasteiger charges)" : $"Not in {ForceFieldName}",
                Types = string.Join(' ', names),
            });
        }
        Entered.Clear();
        foreach (var x in r.GetProperty("entered").EnumerateArray()) Entered.Add(x.GetString()! + " · estimated");
        foreach (var x in r.GetProperty("imported_files").EnumerateArray()) Entered.Add("imported: " + Path.GetFileName(x.GetString()!));
        if (r.TryGetProperty("filled_terms", out var ft))
            foreach (var x in ft.EnumerateArray()) Entered.Add(x.GetString()!.Replace("filled: ", "") + " · filled (the force field has none)");
        if (r.TryGetProperty("by_analogy", out var an))
            foreach (var x in an.EnumerateArray())
                Entered.Add(string.Join(' ', x.GetString()!.Split(' ').Select(FieldNames.Short)) + " · by analogy (estimated)");
        Raise(nameof(HasEntered));
        Notes.Clear();
        foreach (var x in r.GetProperty("notes").EnumerateArray()) Notes.Add(x.GetString()!);

        var src = Path.GetFileName(Str(r, "typing"));
        FooterTyper = $"Typer: {(src.Length > 0 ? src : "atom names")} · {r.GetProperty("rules").GetInt32():N0} SMARTS rules · charges {ChargeModes[_chargeMode].ToLowerInvariant()}";
        var refs = r.GetProperty("references").EnumerateArray().Select(x => x.GetString()!).FirstOrDefault() ?? "";
        FooterDb = $"Parameters: {ForceFieldName} {Str(r, "version")}".TrimEnd() + (refs.Length > 0 ? " · " + (refs.Length > 90 ? refs[..90] + "…" : refs) : "");
        if (r.TryGetProperty("energy", out var en))
            EnergyText = string.Format(Inv, "Energy of the shown frame (kcal/mol): bond {0:F1} · angle {1:F1} · dihedral {2:F1} · improper {3:F1} · vdW {4:F1} · Coulomb {5:F1} · total {6:F1}",
                en.GetProperty("bond").GetDouble(), en.GetProperty("angle").GetDouble(), en.GetProperty("dihedral").GetDouble(), en.GetProperty("improper").GetDouble(),
                en.GetProperty("vdw").GetDouble(), en.GetProperty("coulomb").GetDouble(), en.GetProperty("total").GetDouble());
        else EnergyText = "";
        Complete = r.GetProperty("complete").GetBoolean();
        Assigned = true;
        Log = Complete ? $"{ForceFieldName}: every atom typed and every parameter found. Relax, Dynamics, Equilibrate and LAMMPS data use it now."
                       : $"{ForceFieldName}: {UntypedText}, {MissingText}. Relax and Dynamics are blocked until these are complete; CAPS does not guess parameters.";
        Raise(nameof(RunLine));
    }
}
