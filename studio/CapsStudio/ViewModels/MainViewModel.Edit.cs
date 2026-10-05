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
    // 0 select, 1 place atom, 2 draw bond, 3 delete, 4 lasso select, 5 translate
    private int _editTool;
    public int EditTool
    {
        get => _editTool;
        set
        {
            if (!Set(ref _editTool, value)) return;
            _bondFirst = -1;
            foreach (var n in new[] { nameof(IsPlaceTool), nameof(IsBondTool), nameof(IsDeleteTool), nameof(IsLassoTool), nameof(IsMoveTool), nameof(IsRotateTool), nameof(EditHint) }) Raise(n);
        }
    }
    public bool IsPlaceTool => _editTool == 1;
    public bool IsBondTool => _editTool == 2;
    public bool IsDeleteTool => _editTool == 3;
    private int _bondFirst = -1;
    public string EditHint => _editTool switch
    {
        1 => $"Place {_buildElement}: click an atom to bond one to it",
        2 => _bondFirst < 0 ? "Draw bond: click the first atom" : $"Draw bond: click the atom to bond to {_bondFirst + 1}",
        3 => "Delete: click an atom",
        4 => "Lasso: drag around atoms to select them (⇧ adds to the selection)",
        5 => "Move: drag the selection (or the molecule under the cursor) in the view plane; release to place it",
        6 => "Rotate: drag the selection (or the molecule under the cursor) about the view axis through its centre; hold X, Y or Z for that axis, ⇧ for 15° steps",
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

    /// <summary>A run (dynamics, minimisation, a protocol …) holds the structure: a change now would be lost, so it waits.</summary>
    private bool EditLocked()
    {
        if (_doc?.LongRunning != true) return false;
        Status = "A run is using this structure: change it when the run finishes (or Stop it)";
        return true;
    }

    private JsonNode? RunEdit(object spec)
    {
        if (_doc == null) return null;
        if (EditLocked()) return null;
        var json = System.Text.Json.JsonSerializer.Serialize(spec);
        var r = JsonNode.Parse(_doc.Edit(json))!;
        if (r["ok"]?.GetValue<bool>() != true)
        {
            EditError = r["error"]?.GetValue<string>() ?? "the edit failed";
            Status = EditError;
            return null;
        }
        EditError = "";
        RecordEdit(json);
        var what = r["what"]?.GetValue<string>() ?? "Edit";
        // H autopilot: the edited atoms' hydrogens put right (missing added, surplus removed), its own undo step
        if (_settings.HAutopilot && HFixAfter(json, r) is { Length: > 0 } hreg)
        {
            var h = JsonNode.Parse(_doc.Edit(System.Text.Json.JsonSerializer.Serialize(new { op = "fix_h", atoms = hreg })))!;
            if (h["ok"]?.GetValue<bool>() == true)
            {
                int ha = (int)((double?)h["h_added"] ?? 0), hr = (int)((double?)h["h_removed"] ?? 0);
                HAutoText = $"+{ha} H, −{hr} H after the last edit";
                what += $" · H autopilot +{ha} H, −{hr} H";
            }
            else HAutoText = "H autopilot · the hydrogens were right";
        }
        if (_settings.AutoClean && CleanAfter(json, r) is { Length: > 0 } region)
        {
            // auto-clean (A): the edited atoms and their nearest neighbours relaxed with UFF, as its own undo step
            var c = JsonNode.Parse(_doc.Edit(System.Text.Json.JsonSerializer.Serialize(new { op = "clean", atoms = region })))!;
            if (c["ok"]?.GetValue<bool>() == true) { what += " · auto-cleaned (UFF)"; Record($"doc.edit(op=\"clean\", atoms=[{string.Join(", ", region)}])"); }
        }
        AfterEdit(what);
        return r;
    }

    // the sketch edits whose hydrogens the autopilot puts right: their atoms and the atoms bonded to them
    private static readonly HashSet<string> HFixedOps = ["add_atom", "bond", "unbond", "element", "charge", "attach", "fuse_ring", "delete"];
    private int[] HFixAfter(string json, JsonNode reply)
    {
        if (_doc == null || JsonNode.Parse(json) is not JsonObject spec || !HFixedOps.Contains((string?)spec["op"] ?? "")) return [];
        var core = new HashSet<int>();
        foreach (var x in reply["added"] as JsonArray ?? []) core.Add((int)(double)x!);
        foreach (var k in new[] { "to", "i", "j", "target" })
            if (spec[k] is JsonValue v && v.TryGetValue<int>(out var a) && a >= 0) core.Add(a);
        if ((string?)spec["op"] != "delete" && spec["atoms"] is JsonArray aa) foreach (var x in aa) if (x is JsonValue v && v.TryGetValue<int>(out var a)) core.Add(a);
        var n = (int)_doc.Summary().Atoms;
        core.RemoveWhere(a => a >= n);
        if (core.Count == 0 && (string?)spec["op"] == "add_atom") core.Add(n - 1);
        var region = new HashSet<int>(core);
        foreach (var a in core) foreach (var (i, d) in _doc.Neighbours(a, 4)) if (d < 1.9) region.Add(i);   // its bonded neighbours
        return region.Where(a => a < n).OrderBy(a => a).ToArray();
    }
    private string _hAutoText = "";
    public string HAutoText { get => _hAutoText; private set => Set(ref _hAutoText, value); }
    public bool HAutopilotOn => _settings.HAutopilot;
    public void ToggleHAutopilot()
    {
        _settings.HAutopilot = !_settings.HAutopilot;
        _settings.Save();
        Raise(nameof(HAutopilotOn));
        HAutoText = _settings.HAutopilot ? "H autopilot on" : "";
        Status = _settings.HAutopilot ? "H autopilot on: after each sketch edit, missing hydrogens are added and surplus ones removed" : "H autopilot off";
    }

    // the builder edits that change bonding: their atoms (and new ones) with the six nearest atoms of each
    private static readonly HashSet<string> CleanedOps = ["add_atom", "bond", "unbond", "element", "charge", "attach", "fuse_ring", "add_h"];
    private int[] CleanAfter(string json, JsonNode reply)
    {
        if (_doc == null || JsonNode.Parse(json) is not JsonObject spec || !CleanedOps.Contains((string?)spec["op"] ?? "")) return [];
        var core = new HashSet<int>();
        foreach (var x in reply["added"] as JsonArray ?? []) core.Add((int)(double)x!);
        foreach (var k in new[] { "to", "i", "j", "target" })
            if (spec[k] is JsonValue v && v.TryGetValue<int>(out var a) && a >= 0) core.Add(a);
        if (spec["atoms"] is JsonArray aa) foreach (var x in aa) if (x is JsonValue v && v.TryGetValue<int>(out var a)) core.Add(a);
        var n = (int)_doc.Summary().Atoms;
        core.RemoveWhere(a => a >= n);
        if (core.Count == 0 && (string?)spec["op"] == "add_atom") core.Add(n - 1);
        var region = new HashSet<int>(core);
        foreach (var a in core)
            foreach (var (i, _) in _doc.Neighbours(a, 6)) region.Add(i);
        return region.Where(a => a < n).OrderBy(a => a).ToArray();
    }

    /// <summary>The window opens the full-screen view (true: follow the running job's live snapshots).</summary>
    public event Action<bool>? FullViewRequested;
    public bool AutoCleanOn => _settings.AutoClean;
    public string AutoCleanText => _settings.AutoClean ? "Auto-clean on · UFF" : "";
    public void ToggleAutoClean()
    {
        _settings.AutoClean = !_settings.AutoClean;
        _settings.Save();
        Raise(nameof(AutoCleanOn)); Raise(nameof(AutoCleanText));
        Status = _settings.AutoClean ? "Auto-clean on: every builder edit is followed by a UFF clean-up of the atoms it touched (A turns it off)"
                                     : "Auto-clean off";
    }

    private void AfterEdit(string what)
    {
        _selection.Clear();
        _labelTexts = null;
        GrownUnsaved = true;
        RefreshSummary();
        RefreshSelection();
        RefreshLayers();
        RefreshHistory();
        if (_selOpen) RefreshStereo();
        if (_ixOpen && _ixLive) RunInteractions();
        // the core keeps the selection when the atoms stay the same (tacticity, clean-up), clears it otherwise
        try { SelectedCount = (int)(JsonNode.Parse(_doc!.SelectionJson())!["count"]?.GetValue<double>() ?? 0); } catch { SelectedCount = 0; }
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
        if (_histOpen) RefreshHistoryPanel();
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
                if (atom != _bondFirst) RunEdit(new { op = "bond", i = _bondFirst, j = atom, order = BondOrderCode });
                _bondFirst = -1;
                Raise(nameof(EditHint));
                break;
            case 3:
                if (atom >= 0) RunEdit(new { op = "delete", atoms = new[] { atom } });
                break;
            case 4:   // a click, not a lasso: that atom alone
            case 5:
                if (atom >= 0) Pick(atom);
                break;
        }
    }

    /// <summary>The picked atoms become the chosen element (Element picker: Replace selected atom).</summary>
    public void ReplacePickedElement()
    {
        var sel = SelectionAtoms();
        if (sel.Length == 0) { Status = "Pick the atoms to change first"; return; }
        RunEdit(new { op = "element", atoms = sel, element = _buildElement });
    }

    public void AddHydrogensAll() { var sel = SelectionAtoms(); RunEdit(sel.Length > 0 ? new { op = "add_h", atoms = (object)sel } : new { op = "add_h", atoms = (object)"" }); }
    // the Modify toolbar (Materials Studio's Modify Element / Bond Type / Hybridization, on the picked atoms)
    public static readonly string[] QuickElements = ["H", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br"];
    /// <summary>The picked atoms made this element (their type follows).</summary>
    public void ModifyElementPicked(string symbol)
    {
        var sel = SelectionAtoms();
        if (sel.Length == 0) { Status = "Select atoms (click, ⇧ click adds, double-click a molecule), then an element"; return; }
        RunEdit(new { op = "element", atoms = sel, element = symbol });
    }
    /// <summary>The bond between the two picked atoms made single, double or triple (made if they were not bonded).</summary>
    public void BondOrderPicked(int order)
    {
        if (_selection.Count != 2) { Status = "Pick the two atoms of a bond (⇧ click), then its order"; return; }
        RunEdit(new { op = "bond", i = _selection[0], j = _selection[1], order });
    }
    public void BreakBondPicked()
    {
        if (_selection.Count != 2) { Status = "Pick the two atoms of a bond (⇧ click), then Break"; return; }
        RunEdit(new { op = "unbond", i = _selection[0], j = _selection[1] });
    }
    public static readonly (string Id, string Name)[] Geometries =
        [("linear", "Linear (sp)"), ("trigonal", "Trigonal planar (sp²)"), ("tetrahedral", "Tetrahedral (sp³)"), ("square_planar", "Square planar"),
         ("trigonal_bipyramidal", "Trigonal bipyramidal"), ("square_pyramidal", "Square pyramidal"), ("octahedral", "Octahedral")];
    /// <summary>The picked atom's neighbours placed at the ideal directions of a geometry (each keeps its bond length).</summary>
    public void GeometryPicked(string geometry)
    {
        if (_selection.Count != 1) { Status = "Pick one atom, then its geometry"; return; }
        RunEdit(new { op = "set_coordination", atom = _selection[0], geometry });
    }

    /// <summary>Delete (key, toolbar, bar, ring): the whole selection.</summary>
    public void DeletePicked() => DeleteSelection();
    /// <summary>A benzene ring fused onto the bond between the two picked atoms (each needs a hydrogen on that side).</summary>
    public void FuseRingPicked()
    {
        if (_selection.Count != 2) { Status = "Pick the two atoms of a bond (⇧ click), then Fuse ring"; return; }
        RunEdit(new { op = "fuse_ring", i = _selection[0], j = _selection[1] });
    }

    // exact geometry (design/boards/InteractionMap "Set exact value"): the picked bond, angle or dihedral made the value
    // typed in the live monitor; the side of the last-picked atom moves
    private string _measureTarget = "";
    public string MeasureTarget { get => _measureTarget; set => Set(ref _measureTarget, value); }
    public bool CanSetMeasure => _selection.Count is >= 2 and <= 4 && _doc != null;
    public void SetMeasured()
    {
        if (!CanSetMeasure) { Status = "Pick 2, 3 or 4 atoms (⇧ click): a bond, an angle or a dihedral"; return; }
        if (!double.TryParse(_measureTarget.Replace("Å", "").Replace("°", "").Trim(), System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out var v))
        { Status = "Type the value to set (Å for a bond, degrees for an angle or dihedral)"; return; }
        var picks = _selection.ToArray();
        if (RunEdit(new { op = "set_geometry", atoms = picks, value = v }) == null) return;
        for (var k = 0; k < picks.Length; ++k) Pick(picks[k], k > 0);   // the monitor shows the new value
    }
    /// <summary>The selection (or every atom) turned about a Cartesian axis through its centre.</summary>
    public void RotateSelection(int axis, double degrees)
    {
        if (_doc == null) return;
        var ax = new double[3];
        ax[axis] = 1;
        var sel = SelectionAtoms();
        RunEdit(sel.Length > 0 ? new { op = "rotate", atoms = (object)sel, axis = ax, degrees }
                                     : new { op = "rotate", atoms = (object)Enumerable.Range(0, (int)_doc.Summary().Atoms).ToArray(), axis = ax, degrees });
    }
    /// <summary>The selection (or every atom) reflected through the plane normal to a Cartesian axis: its mirror image.</summary>
    public void MirrorSelection(int axis)
    {
        if (_doc == null) return;
        var n = new double[3];
        n[axis] = 1;
        var sel = SelectionAtoms();
        RunEdit(sel.Length > 0 ? new { op = "mirror", atoms = (object)sel, normal = n }
                                     : new { op = "mirror", atoms = (object)Enumerable.Range(0, (int)_doc.Summary().Atoms).ToArray(), normal = n });
    }
    /// <summary>The picked stereocentre made R or S (inverted only when it is the other).</summary>
    public void MakePicked(string rs)
    {
        if (_selection.Count != 1) { Status = "Pick one stereocentre, then make it " + rs; return; }
        RunEdit(new { op = "set_rs", centre = _selection[0], to = rs });
    }

    public void InvertPicked()
    {
        if (_selection.Count != 1) { Status = "Pick one tetrahedral centre to invert"; return; }
        RunEdit(new { op = "invert", centre = _selection[0] });
    }

    public async Task AutoClean()
    {
        if (_doc == null || Busy) return;
        var doc = _doc;
        var sel = SelectionAtoms();
        var atoms = sel.Length > 0 ? sel : null;
        if (atoms == null && doc.AtomStates() is var st && st.Any(v => (v & 4) != 0))   // locked layers are held
            atoms = Enumerable.Range(0, st.Length).Where(i => (st[i] & 4) == 0).ToArray();
        Status = "Cleaning up with UFF…";
        var text = await Task.Run(() => doc.Edit(atoms == null ? "{\"op\":\"clean\"}" : System.Text.Json.JsonSerializer.Serialize(new { op = "clean", atoms })));
        var r = JsonNode.Parse(text)!;
        if (r["ok"]?.GetValue<bool>() != true) { EditError = r["error"]?.GetValue<string>() ?? "clean-up failed"; Status = EditError; return; }
        Record(atoms == null ? "doc.edit(op=\"clean\")" : $"doc.edit(op=\"clean\", atoms=[{string.Join(", ", atoms)}])");
        AfterEdit(r["what"]!.GetValue<string>());
    }

    public void UndoEdit(bool redo)
    {
        if (_doc == null) return;
        if (EditLocked()) return;
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
