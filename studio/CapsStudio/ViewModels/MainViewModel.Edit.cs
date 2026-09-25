using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>An element of the picker's periodic table.</summary>
public sealed record PeriodicCell(int Z, string Symbol, int Row, int Column, bool Common, string Colour)
{
    public string Tip => $"{Symbol}, atomic number {Z}";
    public double X => (Column - 1) * 41;
    public double Y => (Row - 1) * 41 + (Row >= 8 ? 10 : 0);
    public string ZText => Z.ToString(System.Globalization.CultureInfo.InvariantCulture);
}

/// <summary>Studio editing (the Main board's builder tools, design/boards/ElementPicker): place atoms, draw bonds,
/// delete, change elements, fill hydrogens, invert a centre, clean up with UFF, undo and redo.</summary>
public sealed partial class MainViewModel
{
    // 0 select, 1 place atom, 2 draw bond, 3 delete
    private int _editTool;
    public int EditTool { get => _editTool; set { if (Set(ref _editTool, value)) { _bondFirst = -1; Raise(nameof(IsPlaceTool)); Raise(nameof(IsBondTool)); Raise(nameof(IsDeleteTool)); Raise(nameof(EditHint)); } } }
    public bool IsPlaceTool => _editTool == 1;
    public bool IsBondTool => _editTool == 2;
    public bool IsDeleteTool => _editTool == 3;
    private int _bondFirst = -1;
    public string EditHint => _editTool switch
    {
        1 => $"Place {_buildElement}: click an atom to bond one to it",
        2 => _bondFirst < 0 ? "Draw bond: click the first atom" : $"Draw bond: click the atom to bond to {_bondFirst + 1}",
        3 => "Delete: click an atom",
        _ => "",
    };

    private string _buildElement = "C";
    private int _buildGeometry, _buildCharge = 1, _buildIsotope;
    /// <summary>The element new atoms get (the toolbar's element button, the Element picker).</summary>
    public string BuildElement { get => _buildElement; set { if (value != null && Set(ref _buildElement, value)) { Raise(nameof(BuildElementColour)); Raise(nameof(BuildElementText)); Raise(nameof(EditHint)); RaisePicker(); } } }
    public string BuildElementColour => ElementColour(_buildElement);
    /// <summary>0 from the neighbours, 1 sp³, 2 sp², 3 sp.</summary>
    public int BuildGeometry { get => _buildGeometry; set { if (Set(ref _buildGeometry, Math.Clamp(value, 0, 3))) Raise(nameof(BuildElementText)); } }
    /// <summary>Formal charge of placed atoms: 0 → −1, 1 → 0, 2 → +1, 3 → +2.</summary>
    public int BuildCharge { get => _buildCharge; set => Set(ref _buildCharge, Math.Clamp(value, 0, 3)); }
    public int BuildIsotope { get => _buildIsotope; set => Set(ref _buildIsotope, value); }
    public string BuildElementText => _buildGeometry switch { 1 => "sp³", 2 => "sp²", 3 => "sp", _ => "auto" };
    private static string ElementColour(string sym) => $"#{Native.ElementColour(Native.ElementNumber(sym)):X6}";

    private string _editError = "";
    public string EditError { get => _editError; private set => Set(ref _editError, value); }
    public ObservableCollection<string> EditHistory { get; } = new();
    public bool CanUndo => EditHistory.Count > 0;
    private bool _canRedo;
    public bool CanRedo { get => _canRedo; private set => Set(ref _canRedo, value); }

    private JsonNode? RunEdit(object spec)
    {
        if (_doc == null) return null;
        var json = System.Text.Json.JsonSerializer.Serialize(spec);
        var r = JsonNode.Parse(_doc.Edit(json))!;
        if (r["ok"]?.GetValue<bool>() != true)
        {
            EditError = r["error"]?.GetValue<string>() ?? "the edit failed";
            Status = EditError;
            return null;
        }
        EditError = "";
        AfterEdit(r["what"]?.GetValue<string>() ?? "Edit");
        return r;
    }

    private void AfterEdit(string what)
    {
        _selection.Clear();
        _labelTexts = null;
        GrownUnsaved = true;
        RefreshSummary();
        RefreshSelection();
        RefreshHistory();
        Status = what;
        RenderRequested?.Invoke();
    }

    private void RefreshHistory()
    {
        EditHistory.Clear();
        if (_doc == null) { Raise(nameof(CanUndo)); CanRedo = false; return; }
        var h = JsonNode.Parse(_doc.History())!;
        foreach (var x in h["undo"]!.AsArray()) EditHistory.Add(x!.GetValue<string>());
        CanRedo = h["redo"]!.AsArray().Count > 0;
        Raise(nameof(CanUndo));
    }

    private int ChargeValue => _buildCharge - 1;

    /// <summary>A click in the view with a builder tool active.</summary>
    public void ToolClick(int atom)
    {
        switch (_editTool)
        {
            case 1:
                if (atom < 0) RunEdit(new { op = "add_atom", to = -1, element = _buildElement, charge = ChargeValue });
                else RunEdit(new { op = "add_atom", to = atom, element = _buildElement, geometry = _buildGeometry switch { 1 => 3, 2 => 2, 3 => 1, _ => 0 }, charge = ChargeValue });
                break;
            case 2:
                if (atom < 0) return;
                if (_bondFirst < 0) { _bondFirst = atom; Pick(atom); Raise(nameof(EditHint)); return; }
                if (atom != _bondFirst) RunEdit(new { op = "bond", i = _bondFirst, j = atom, order = 1 });
                _bondFirst = -1;
                Raise(nameof(EditHint));
                break;
            case 3:
                if (atom >= 0) RunEdit(new { op = "delete", atoms = new[] { atom } });
                break;
        }
    }

    /// <summary>The picked atoms become the chosen element (Element picker: Replace selected atom).</summary>
    public void ReplacePickedElement()
    {
        if (_selection.Count == 0) { Status = "Pick the atoms to change first"; return; }
        RunEdit(new { op = "element", atoms = _selection.ToArray(), element = _buildElement });
    }

    public void AddHydrogensAll() => RunEdit(_selection.Count > 0 ? new { op = "add_h", atoms = (object)_selection.ToArray() } : new { op = "add_h", atoms = (object)"" });
    public void DeletePicked() { if (_selection.Count > 0) RunEdit(new { op = "delete", atoms = _selection.ToArray() }); }
    public void InvertPicked()
    {
        if (_selection.Count != 1) { Status = "Pick one tetrahedral centre to invert"; return; }
        RunEdit(new { op = "invert", centre = _selection[0] });
    }

    public async Task AutoClean()
    {
        if (_doc == null || Busy) return;
        var doc = _doc;
        var atoms = _selection.Count > 0 ? _selection.ToArray() : null;
        Status = "Cleaning up with UFF…";
        var text = await Task.Run(() => doc.Edit(atoms == null ? "{\"op\":\"clean\"}" : System.Text.Json.JsonSerializer.Serialize(new { op = "clean", atoms })));
        var r = JsonNode.Parse(text)!;
        if (r["ok"]?.GetValue<bool>() != true) { EditError = r["error"]?.GetValue<string>() ?? "clean-up failed"; Status = EditError; return; }
        AfterEdit(r["what"]!.GetValue<string>());
    }

    public void UndoEdit(bool redo)
    {
        if (_doc == null) return;
        if (!_doc.Undo(redo)) { Status = redo ? "Nothing to redo" : "Nothing to undo"; return; }
        AfterEdit(redo ? "Redone" : "Undone");
    }

    // ---------------------------------------------------------------- the Element picker

    private bool _pickerOpen;
    public bool ElementPickerOpen { get => _pickerOpen; set { if (Set(ref _pickerOpen, value)) RaisePicker(); } }
    private static readonly string[] Syms = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr " +
        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu " +
        "Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr " +
        "Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").Split(' ');
    private static readonly HashSet<string> Common = ["H", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "Na", "K", "Ca", "Li", "I", "Zn", "Ti"];
    private List<PeriodicCell>? _table;
    public List<PeriodicCell> PeriodicTable => _table ??= Syms.Select((s, i) => { var (r, c) = TablePos(i + 1); return new PeriodicCell(i + 1, s, r, c, Common.Contains(s), ElementColour(s)); }).ToList();

    private static (int Row, int Col) TablePos(int z)
    {
        if (z == 1) return (1, 1);
        if (z == 2) return (1, 18);
        foreach (var (start, period) in new[] { (3, 2), (11, 3) })
            if (z >= start && z < start + 8) { var i = z - start; return (period, i < 2 ? i + 1 : i + 11); }
        foreach (var (start, period) in new[] { (19, 4), (37, 5) })
            if (z >= start && z < start + 18) return (period, z - start + 1);
        foreach (var (start, period) in new[] { (55, 6), (87, 7) })
        {
            var i = z - start;
            if (i < 0 || i >= 32) continue;
            if (i < 2) return (period, i + 1);
            if (i < 17) return (period + 3, i + 1);   // the f-block rows below the table
            return (period, i - 13);
        }
        return (1, 1);
    }

    private void RaisePicker()
    {
        foreach (var n in new[] { nameof(PickerSymbol), nameof(PickerName), nameof(PickerLine), nameof(PickerRadius), nameof(PickerIsotopes), nameof(PickerPlaceText) }) Raise(n);
    }

    public string PickerSymbol => _buildElement;
    public string PickerName => ElementTitles.TryGetValue(_buildElement, out var n) ? n : _buildElement;
    public string PickerLine
    {
        get
        {
            var z = Native.ElementNumber(_buildElement);
            return $"Z = {z} · {Native.ElementMass(z).ToString("0.###", CultureInfo.InvariantCulture)} u";
        }
    }
    public string PickerRadius
    {
        get
        {
            var z = Native.ElementNumber(_buildElement);
            return $"covalent {Native.ElementCovalent(z).ToString("0.00", CultureInfo.InvariantCulture)} Å · van der Waals {Native.ElementVdw(z).ToString("0.00", CultureInfo.InvariantCulture)} Å";
        }
    }
    public string PickerIsotopes => Isotopes.TryGetValue(_buildElement, out var t) ? t : "natural abundance";
    public string PickerPlaceText => $"Place {_buildElement}";

    private static readonly Dictionary<string, string> Isotopes = new()
    {
        ["H"] = "¹H 99.985 % · ²H 0.015 %", ["C"] = "¹²C 98.93 % · ¹³C 1.07 %", ["N"] = "¹⁴N 99.636 % · ¹⁵N 0.364 %",
        ["O"] = "¹⁶O 99.757 % · ¹⁷O 0.038 % · ¹⁸O 0.205 %", ["S"] = "³²S 94.99 % · ³⁴S 4.25 %", ["Cl"] = "³⁵Cl 75.76 % · ³⁷Cl 24.24 %",
        ["Br"] = "⁷⁹Br 50.69 % · ⁸¹Br 49.31 %", ["Si"] = "²⁸Si 92.23 % · ²⁹Si 4.68 % · ³⁰Si 3.09 %",
    };
    private static readonly Dictionary<string, string> ElementTitles = new()
    {
        ["H"] = "Hydrogen", ["He"] = "Helium", ["Li"] = "Lithium", ["Be"] = "Beryllium", ["B"] = "Boron", ["C"] = "Carbon", ["N"] = "Nitrogen", ["O"] = "Oxygen",
        ["F"] = "Fluorine", ["Ne"] = "Neon", ["Na"] = "Sodium", ["Mg"] = "Magnesium", ["Al"] = "Aluminium", ["Si"] = "Silicon", ["P"] = "Phosphorus", ["S"] = "Sulfur",
        ["Cl"] = "Chlorine", ["Ar"] = "Argon", ["K"] = "Potassium", ["Ca"] = "Calcium", ["Ti"] = "Titanium", ["Fe"] = "Iron", ["Cu"] = "Copper", ["Zn"] = "Zinc",
        ["Br"] = "Bromine", ["Ag"] = "Silver", ["Sn"] = "Tin", ["I"] = "Iodine", ["Au"] = "Gold", ["Pt"] = "Platinum", ["Co"] = "Cobalt", ["Ni"] = "Nickel",
        ["Mn"] = "Manganese", ["Cr"] = "Chromium", ["Pb"] = "Lead", ["Hg"] = "Mercury", ["Se"] = "Selenium", ["Ge"] = "Germanium", ["Ga"] = "Gallium", ["As"] = "Arsenic",
    };

    /// <summary>Jump in the picker by typing a symbol ("Cl").</summary>
    public bool PickerJump(string text)
    {
        var sym = PeriodicTable.FirstOrDefault(c => string.Equals(c.Symbol, text, StringComparison.OrdinalIgnoreCase))?.Symbol;
        if (sym == null) return false;
        BuildElement = sym;
        return true;
    }
}
