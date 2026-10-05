using System.Globalization;
using System.Text;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>One action of the ring menu: its icon, its word and what it does.</summary>
public sealed record RingItem(string Icon, string Word, string Id);

/// <summary>The selection bar and ring menu (design/boards/SelectionBar): one selection, whatever made it — clicked atoms,
/// a lasso, a box, a double-clicked molecule, an expression, a named set — with a bar floating over it (count, what it is,
/// how it was made, the actions) and a ring of the same actions under the pointer on a held right button. Hide and ghost
/// are view states of atoms (the structure, its exports and calculations keep every atom).</summary>
public sealed partial class MainViewModel
{
    /// <summary>The whole selection: the core's (lasso, box, expression, set) and the clicked atoms, sorted, without ghosted
    /// or hidden atoms.</summary>
    public int[] SelectionAtoms()
    {
        if (_doc == null) return [];
        var set = new SortedSet<int>(_selection);
        if (_selCount > 0)
            try
            {
                if (JsonNode.Parse(_doc.SelectionJson())?["indices"] is JsonArray a)
                    foreach (var x in a) set.Add((int)((double?)x ?? -1));
            }
            catch (Exception) { }
        set.Remove(-1);
        if (set.Count == 0) return [];
        var st = _doc.AtomStates();
        return set.Where(i => i < st.Length && st[i] == 0).ToArray();
    }

    private int _barCount;
    private string _barWhat = "", _barHow = "";
    private int[] _barAtoms = [];
    public int SelBarCount { get => _barCount; private set { if (Set(ref _barCount, value)) { Raise(nameof(HasSelBar)); Raise(nameof(SelBarCountText)); } } }
    public string SelBarCountText => _barCount.ToString("N0", CultureInfo.InvariantCulture);
    /// <summary>What the selection is and how it was made: "molecule 4 · lasso".</summary>
    public string SelBarWhat { get => _barWhat; private set => Set(ref _barWhat, value); }
    public string SelBarHow { get => _barHow; private set => Set(ref _barHow, value); }
    public bool HasSelBar => _barCount > 0 && _doc != null && IsStudio && !Busy;
    /// <summary>A spread of the selected atoms (at most 4000) the view projects to place the bar.</summary>
    public int[] SelBarSample
    {
        get
        {
            if (_barAtoms.Length <= 4000) return _barAtoms;
            var step = _barAtoms.Length / 4000.0;
            return Enumerable.Range(0, 4000).Select(k => _barAtoms[(int)(k * step)]).ToArray();
        }
    }
    public event Action? SelectionBarChanged;
    private string _pickHow = "click";

    private bool _barQueued;
    /// <summary>Many changes in one turn (a count, then how it was made) make one refresh of the bar.</summary>
    private void QueueSelBar()
    {
        if (_barQueued) return;
        _barQueued = true;
        Avalonia.Threading.Dispatcher.UIThread.Post(() => { _barQueued = false; RefreshSelBar(); RefreshViewStates(); }, Avalonia.Threading.DispatcherPriority.Background);
    }

    /// <summary>The bar's count, description and origin, after any change of the selection.</summary>
    public void RefreshSelBar()
    {
        if (HasTags) RefreshTags();   // each chip's selected count
        if (_brushOpen) RefillBrush();   // the histograms' selected share
        _barAtoms = SelectionAtoms();
        SelBarCount = _barAtoms.Length;
        if (_barAtoms.Length == 0) { SelBarWhat = SelBarHow = ""; LightLayers(_barAtoms); SelectionBarChanged?.Invoke(); return; }
        SelBarWhat = DescribeAtoms(_barAtoms);
        var how = new List<string>();
        if (_selCount > 0 && _selHud.Length > 0) how.Add(_selHud);
        else if (_selCount > 0) how.Add("selection");
        if (_selection.Count > 0) how.Add(_pickHow);
        SelBarHow = string.Join(" + ", how.Distinct());
        LightLayers(_barAtoms);
        SelectionBarChanged?.Invoke();
    }

    /// <summary>"molecule 4", "3 molecules", "12 C atoms", "C H · 5 molecules" — what a set of atoms is.</summary>
    private string DescribeAtoms(int[] atoms)
    {
        if (_doc == null || atoms.Length == 0) return "";
        try
        {
            var mol = _doc.MoleculeIds();
            var mols = new HashSet<long>();
            foreach (var a in atoms) if (a < mol.Length) mols.Add(mol[a]);
            // a whole molecule, or part of one
            if (mols.Count == 1)
            {
                var m = mols.First();
                var whole = mol.Count(x => x == m) == atoms.Length;
                return whole ? $"molecule {m}" : $"{Elements(atoms)} in molecule {m}";
            }
            return atoms.Length <= 20000 ? $"{Elements(atoms)} · {mols.Count:N0} molecules" : $"{mols.Count:N0} molecules";
        }
        catch (Exception) { return $"{atoms.Length:N0} atoms"; }
    }

    private string Elements(int[] atoms)
    {
        var el = new SortedDictionary<string, int>(StringComparer.Ordinal);
        foreach (var a in atoms.Take(20000)) { var s = _doc!.Atom(a).ElementSymbol; el[s] = el.GetValueOrDefault(s) + 1; }
        return el.Count == 1 ? $"{el.First().Value:N0} {el.First().Key}" : string.Join(" ", el.OrderByDescending(p => p.Value).Take(4).Select(p => p.Key));
    }

    /// <summary>Clicking a new atom (not adding) starts a new selection: the core's is cleared too.</summary>
    public void StartSelection(int hit, string how)
    {
        _pickHow = how;
        if (_selCount > 0 && _doc != null) { _doc.Select("{\"mode\":\"none\"}"); SelectedCount = 0; SelectHud = ""; }
        Pick(hit, false);
    }

    /// <summary>Double-click: the atom's molecule; with ⌥ every atom of its type (force-field type, else the file's type).</summary>
    public void SelectLike(int atom, bool sameType, bool add)
    {
        if (_doc == null || atom < 0) return;
        JsonObject o;
        string how;
        if (sameType)
        {
            var a = _doc.Atom(atom);
            o = new JsonObject { ["mode"] = "type", ["pattern"] = a.Type.ToString(CultureInfo.InvariantCulture), ["op"] = add ? "add" : "replace" };
            how = $"type {a.Type} ({a.ElementSymbol})";
        }
        else
        {
            o = new JsonObject { ["mode"] = "molecule", ["atoms"] = new JsonArray(atom), ["op"] = add ? "add" : "replace" };
            how = "double-click";
        }
        var r = JsonNode.Parse(_doc.Select(o.ToJsonString()))!;
        if (r["ok"]?.GetValue<bool>() != true) { Status = r["error"]?.GetValue<string>() ?? "cannot select"; return; }
        if (!add) { _selection.Clear(); RefreshSelection(); }
        SelectedCount = (int)r["count"]!.GetValue<double>();
        SelectHud = how;
        RefreshSelBar();
        Status = $"{SelectedCount:N0} atoms selected · {how}";
        RenderRequested?.Invoke();
    }

    /// <summary>Esc: nothing selected.</summary>
    public void ClearAllSelection()
    {
        _selection.Clear();
        RefreshSelection();
        if (_selCount > 0) ClearDocSelection();
        RefreshSelBar();
    }

    // ---------------------------------------------------------------- view states

    private int _hiddenCount, _ghostCount;
    public int HiddenCount { get => _hiddenCount; private set { if (Set(ref _hiddenCount, value)) { Raise(nameof(HasHiddenAtoms)); Raise(nameof(HiddenChip)); } } }
    public int GhostCount { get => _ghostCount; private set { if (Set(ref _ghostCount, value)) { Raise(nameof(HasHiddenAtoms)); Raise(nameof(HiddenChip)); } } }
    public bool HasHiddenAtoms => _hiddenCount + _ghostCount > 0;
    public string HiddenChip => string.Join(" · ", new[] { _hiddenCount > 0 ? $"{_hiddenCount:N0} hidden" : "", _ghostCount > 0 ? $"{_ghostCount:N0} ghosted" : "" }.Where(x => x.Length > 0));

    private void RefreshViewStates()
    {
        if (_doc == null) { HiddenCount = GhostCount = 0; return; }
        var st = _doc.AtomStates();
        HiddenCount = st.Count(v => (v & 3) == 2);
        GhostCount = st.Count(v => (v & 3) == 1);
        Raise(nameof(LayerChip));
    }

    private void SetStates(int[]? atoms, int state, string what)
    {
        if (_doc == null) return;
        _doc.SetAtomState(atoms, state);
        Record(atoms == null ? $"doc.set_atom_state(None, {state})" : $"doc.set_atom_state({atoms.Length} atoms, {state})");
        _selection.Clear();
        RefreshSelection();
        if (_selCount > 0) ClearDocSelection();
        RefreshViewStates();
        RefreshLayerStates();
        Raise(nameof(LayerChip));
        RefreshSelBar();
        Status = what;
        RenderRequested?.Invoke();
    }

    public void HideSelection()
    {
        var a = SelectionAtoms();
        if (a.Length == 0) { Status = "Select atoms to hide (click, double-click a molecule, ⌥ drag a box)"; return; }
        SetStates(a, 2, $"{a.Length:N0} atoms hidden · Show all brings them back (they stay in the structure and its exports)");
    }

    public void GhostSelection()
    {
        var a = SelectionAtoms();
        if (a.Length == 0) { Status = "Select atoms to ghost"; return; }
        SetStates(a, 1, $"{a.Length:N0} atoms ghosted: faint, and clicks go through them");
    }

    public void ShowOnlySelection()
    {
        var a = SelectionAtoms();
        if (a.Length == 0 || _doc == null) { Status = "Select the atoms to keep in view"; return; }
        _doc.SetAtomState(null, 2);
        SetStates(a, 0, $"Only {a.Length:N0} atoms shown · Show all brings the rest back");
        RefreshViewStates();
    }

    public void ShowAllAtoms()
    {
        if (!HasHiddenAtoms) { Status = "Every atom is shown"; return; }
        SetStates(null, 0, "Every atom shown");
    }

    /// <summary>The selected atoms as XYZ text (element and Å), for the clipboard.</summary>
    public string SelectionXyz()
    {
        var a = SelectionAtoms();
        if (a.Length == 0 || _doc == null) return "";
        var inv = CultureInfo.InvariantCulture;
        var sb = new StringBuilder();
        sb.Append(a.Length.ToString(inv)).Append('\n').Append(Title).Append(" · ").Append(SelBarWhat).Append('\n');
        foreach (var i in a)
        {
            var x = _doc.Atom(i);
            sb.Append(x.ElementSymbol).Append(' ').Append(x.X.ToString("0.00000", inv)).Append(' ').Append(x.Y.ToString("0.00000", inv)).Append(' ').Append(x.Z.ToString("0.00000", inv)).Append('\n');
        }
        return sb.ToString();
    }

    public void DeleteSelection()
    {
        var a = SelectionAtoms();
        if (a.Length == 0) { Status = "Select atoms to delete"; return; }
        RunEdit(new { op = "delete", atoms = a });
    }

    /// <summary>The ring's actions: on the selection, or (nothing selected) on the view.</summary>
    public IReadOnlyList<RingItem> RingItems() => _barCount > 0
        ? [new("eye", "Hide", "hide"), new("layers", "Ghost", "ghost"), new("tag", "Tag", "tag"), new("copy", "Copy", "copy"),
           new("wand", "Clean up", "clean"), new("frame", "Frame", "frame"), new("filter", "Show only", "only"), new("xcircle", "Delete", "delete")]
        : [new("eye", "Show all", "showall"), new("cursor", "Select all", "all"), new("frame", "Fit view", "fit"), new("mirror", "Invert", "invert")];

    /// <summary>Runs a ring or bar action; "copy" returns the text to put on the clipboard.</summary>
    public async Task<string?> RunSelectionAction(string id)
    {
        switch (id)
        {
            case "hide": HideSelection(); break;
            case "ghost": GhostSelection(); break;
            case "tag": TagSelection(); break;
            case "copy":
                var t = SelectionXyz();
                CopyToTray();   // the piece in the clipboard tray, to place as a stamp
                if (t.Length > 0) Status = $"{_barCount:N0} atoms copied: in the clipboard tray (pick it to place it) and as XYZ text";
                return t;
            case "clean": await AutoClean(); break;
            case "frame": FrameSelection(); break;
            case "only": ShowOnlySelection(); break;
            case "delete": DeleteSelection(); break;
            case "showall": ShowAllAtoms(); break;
            case "all":
                if (_doc == null) break;
                var r = JsonNode.Parse(_doc.Select("{\"mode\":\"all\",\"op\":\"replace\"}"))!;
                SelectedCount = (int)(r["count"]?.GetValue<double>() ?? 0);
                SelectHud = "all";
                RefreshSelBar();
                RenderRequested?.Invoke();
                break;
            case "fit": FrameSelection(); break;
            case "invert": RunSelect("invert"); RefreshSelBar(); break;
        }
        return null;
    }
}
