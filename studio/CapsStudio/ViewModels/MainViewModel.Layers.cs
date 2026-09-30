using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>One row of the layers list: a kind of molecule (its molecules under it) or one molecule.</summary>
public sealed class LayerRow : ObservableObject
{
    public LayerRow(string name, string sub, int atoms, string colour, int level, bool isKind, long[] mols, double[] z)
    {
        Name = name; Sub = sub; Atoms = atoms; Colour = colour; Level = level; IsKind = isKind; Molecules = mols; Z = z;
    }
    public string Name { get; }
    public string Sub { get; }
    public int Atoms { get; }
    public string AtomsText => Atoms.ToString("N0", CultureInfo.InvariantCulture);
    public string Colour { get; }
    public Avalonia.Media.IBrush Brush => Avalonia.Media.Brush.Parse(Colour);
    public int Level { get; }
    public Avalonia.Thickness Indent => new(8 + 14 * Level, 0, 4, 0);
    public bool IsKind { get; }
    /// <summary>The molecule ids the row stands for.</summary>
    public long[] Molecules { get; }
    /// <summary>Its atoms along z (8 bins, peak 1).</summary>
    public double[] Z { get; }
    public bool HasChildren { get; init; }
    public bool IsMore { get; init; }
    private bool _expanded;
    public bool Expanded { get => _expanded; set { if (Set(ref _expanded, value)) Raise(nameof(Chevron)); } }
    public string Chevron => _expanded ? "chev" : "chevr";
    private string _state = "shown";
    /// <summary>shown, ghost, hidden or mixed.</summary>
    public string State { get => _state; set { if (Set(ref _state, value)) { Raise(nameof(IsGhost)); Raise(nameof(IsHidden)); Raise(nameof(GhostChip)); Raise(nameof(HiddenChip)); Raise(nameof(RowOpacity)); } } }
    public bool IsGhost => _state == "ghost";
    public bool IsHidden => _state == "hidden";
    // the words on kind rows; a molecule row shows its state by its eye and its dimming
    public bool GhostChip => IsKind && IsGhost;
    public bool HiddenChip => IsKind && IsHidden;
    public double RowOpacity => _state == "hidden" ? 0.5 : _state == "ghost" ? 0.75 : 1;
    private bool _locked;
    public bool Locked { get => _locked; set => Set(ref _locked, value); }
    private int _selected;
    /// <summary>How many of its atoms are selected (lit when any, strongly when all).</summary>
    public int Selected { get => _selected; set { if (Set(ref _selected, value)) { Raise(nameof(IsLit)); Raise(nameof(IsFullyLit)); } } }
    public bool IsLit => _selected > 0;
    public bool IsFullyLit => _selected > 0 && _selected >= Atoms;
    public bool Held { get; init; }
}

/// <summary>Layers (design/boards/Layers): the structure as what it is made of — kinds of molecule (the slab, the film and
/// its chains) with each molecule under its kind — each with a profile along z, its atom count and its state. The eye
/// cycles shown → ghost (faint, not pickable) → hidden; the pin locks a layer against picking and edits. Nothing leaves the
/// structure, its exports or calculations. Selecting in the view lights the rows; a row selects its atoms.</summary>
public sealed partial class MainViewModel
{
    private static readonly string[] LayerColours = ["#6CC4D8", "#E0A060", "#9B7AD5", "#7CC784", "#E07A5F", "#D4B85A", "#8FA8C8", "#C98BB9"];
    public ObservableCollection<LayerRow> LayerRows { get; } = new();
    private readonly List<LayerRow> _layerKinds = new();
    private readonly Dictionary<LayerRow, List<LayerRow>> _layerChildren = new();
    private readonly HashSet<string> _expandedKinds = new();
    private long[] _layerMolIds = [];
    public bool HasLayers => LayerRows.Count > 0;
    private string _layerAxis = "";
    public string LayerCaption => _layerAxis.Length == 0 ? "" : $"z along the {_layerAxis} · atoms · state";

    /// <summary>The layers of the document, rebuilt (states and selections are refreshed in place by RefreshLayerStates).</summary>
    public void RefreshLayers()
    {
        LayerRows.Clear();
        _layerKinds.Clear();
        _layerChildren.Clear();
        if (_doc == null) { Raise(nameof(HasLayers)); return; }
        try
        {
            var j = JsonNode.Parse(_doc.LayersJson())!;
            _layerAxis = (string?)j["z_axis"] ?? "";
            _layerMolIds = _doc.MoleculeIds();
            var k = 0;
            foreach (var kind in j["kinds"]!.AsArray())
            {
                var colour = LayerColours[k++ % LayerColours.Length];
                var mols = kind!["molecules"]!.AsArray();
                var name = (string?)kind["name"] ?? "";
                var elements = (string?)kind["elements"] ?? "";
                double[] Z(JsonNode m) => m["z"]!.AsArray().Select(v => (double)v!).ToArray();
                if (mols.Count == 1)
                {
                    var m = mols[0]!;
                    var id = (long)(double)m["id"]!;
                    var held = (bool?)m["held"] == true;
                    var row = new LayerRow($"{name}", $"molecule {id} · {elements}{(held ? " · held" : "")}", (int)(double)m["atoms"]!, colour, 0, true, [id], Z(m)) { Held = held };
                    Apply(row, m);
                    _layerKinds.Add(row);
                    continue;
                }
                // a kind: its profile sums its molecules'
                var zs = new double[8];
                foreach (var m in mols) { var z = Z(m!); for (var b = 0; b < 8; ++b) zs[b] += z[b]; }
                var peak = zs.Max();
                if (peak > 0) for (var b = 0; b < 8; ++b) zs[b] /= peak;
                var kindRow = new LayerRow($"{mols.Count:N0} × {name}", elements, (int)(double)kind["atoms"]!, colour, 0, true,
                                           mols.Select(m => (long)(double)m!["id"]!).ToArray(), zs) { HasChildren = true };
                var children = new List<LayerRow>();
                foreach (var m in mols.Take(400))
                {
                    var id = (long)(double)m!["id"]!;
                    var child = new LayerRow($"molecule {id}", elements, (int)(double)m["atoms"]!, colour, 1, false, [id], Z(m)) { Held = (bool?)m["held"] == true };
                    Apply(child, m);
                    children.Add(child);
                }
                if (mols.Count > 400) children.Add(new LayerRow($"… {mols.Count - 400:N0} more (select the kind, or use a query)", "", 0, colour, 1, false, [], new double[8]) { IsMore = true });
                _layerChildren[kindRow] = children;
                AggregateKind(kindRow, children);
                kindRow.Expanded = _expandedKinds.Contains(name) || (mols.Count <= 12 && !_expandedKinds.Contains("-" + name));
                _layerKinds.Add(kindRow);
            }
        }
        catch (Exception e) { Status = "Layers: " + e.Message; }
        Flatten();
        Raise(nameof(LayerCaption));
    }

    private static void Apply(LayerRow r, JsonNode m)
    {
        r.State = (string?)m["state"] ?? "shown";
        r.Locked = (bool?)m["locked"] == true;
        r.Selected = (int)((double?)m["selected"] ?? 0);
    }

    private static void AggregateKind(LayerRow kind, List<LayerRow> children)
    {
        var real = children.Where(c => !c.IsMore).ToList();
        var states = real.Select(c => c.State).Distinct().ToList();
        kind.State = states.Count == 1 ? states[0] : "mixed";
        kind.Locked = real.Count > 0 && real.All(c => c.Locked);
        kind.Selected = real.Sum(c => c.Selected);
    }

    private void Flatten()
    {
        LayerRows.Clear();
        foreach (var k in _layerKinds)
        {
            LayerRows.Add(k);
            if (k.Expanded && _layerChildren.TryGetValue(k, out var ch)) foreach (var c in ch) LayerRows.Add(c);
        }
        Raise(nameof(HasLayers));
    }

    /// <summary>States, locks and selected counts from the core, without rebuilding the list (after a hide, a lock, a pick).</summary>
    public void RefreshLayerStates()
    {
        if (_doc == null || _layerKinds.Count == 0) return;
        try
        {
            var j = JsonNode.Parse(_doc.LayersJson())!;
            var byId = new Dictionary<long, JsonNode>();
            foreach (var kind in j["kinds"]!.AsArray()) foreach (var m in kind!["molecules"]!.AsArray()) byId[(long)(double)m!["id"]!] = m;
            if (byId.Count != _layerKinds.Sum(k => k.Molecules.Length)) { RefreshLayers(); return; }   // the molecules changed
            foreach (var k in _layerKinds)
            {
                if (_layerChildren.TryGetValue(k, out var ch))
                {
                    foreach (var c in ch) if (!c.IsMore && byId.TryGetValue(c.Molecules[0], out var m)) Apply(c, m);
                    if (ch.Count > 0 && ch[^1].IsMore)
                    {
                        // beyond the listed children: the kind from all its molecules
                        var all = k.Molecules.Select(id => byId[id]).ToList();
                        var states = all.Select(m => (string?)m["state"]).Distinct().ToList();
                        k.State = states.Count == 1 ? states[0] ?? "shown" : "mixed";
                        k.Locked = all.All(m => (bool?)m["locked"] == true);
                        k.Selected = all.Sum(m => (int)((double?)m["selected"] ?? 0));
                    }
                    else AggregateKind(k, ch);
                }
                else if (byId.TryGetValue(k.Molecules[0], out var m)) Apply(k, m);
            }
        }
        catch (Exception) { }
    }

    /// <summary>The rows lit by the selection (from the selected atoms and the molecule ids: no call into the core).</summary>
    private void LightLayers(int[] selected)
    {
        if (_layerKinds.Count == 0) return;
        var per = new Dictionary<long, int>();
        var ids = _layerMolIds;
        foreach (var a in selected) if (a < ids.Length) per[ids[a]] = per.GetValueOrDefault(ids[a]) + 1;
        foreach (var k in _layerKinds)
        {
            if (_layerChildren.TryGetValue(k, out var ch)) foreach (var c in ch) if (!c.IsMore) c.Selected = per.GetValueOrDefault(c.Molecules[0]);
            k.Selected = k.Molecules.Sum(m => per.GetValueOrDefault(m));
        }
    }

    public void ToggleLayer(LayerRow r)
    {
        if (!r.HasChildren) return;
        r.Expanded = !r.Expanded;
        var key = r.Name[(r.Name.IndexOf('×') + 2)..];
        _expandedKinds.Remove(key); _expandedKinds.Remove("-" + key);
        _expandedKinds.Add(r.Expanded ? key : "-" + key);
        Flatten();
    }

    /// <summary>The atoms of a layer row (its molecules).</summary>
    private int[] LayerAtoms(LayerRow r)
    {
        if (r.Molecules.Length == 0) return [];
        var want = r.Molecules.ToHashSet();
        var ids = _layerMolIds;
        var list = new List<int>();
        for (var i = 0; i < ids.Length; ++i) if (want.Contains(ids[i])) list.Add(i);
        return list.ToArray();
    }

    /// <summary>The eye: shown → ghost → hidden → shown (mixed goes to shown).</summary>
    public void CycleLayer(LayerRow r)
    {
        if (_doc == null || r.IsMore) return;
        var next = r.State switch { "shown" => 1, "ghost" => 2, _ => 0 };
        var atoms = LayerAtoms(r);
        _doc.SetAtomState(atoms, next);
        Record($"doc.set_atom_state({atoms.Length} atoms of {r.Name}, {next})");
        AfterLayerChange($"{r.Name}: {(next == 0 ? "shown" : next == 1 ? "ghosted (faint, not pickable)" : "hidden")}");
    }

    /// <summary>The pin: locked against picking and edits (the held slab, say), or free.</summary>
    public void LockLayer(LayerRow r)
    {
        if (_doc == null || r.IsMore) return;
        var atoms = LayerAtoms(r);
        _doc.SetAtomLock(atoms, !r.Locked);
        Record($"doc.lock({atoms.Length} atoms of {r.Name}, {(!r.Locked ? "True" : "False")})");
        AfterLayerChange(r.Locked ? $"{r.Name} free again" : $"{r.Name} locked: not picked, not edited (clean-up holds it)");
    }

    private void AfterLayerChange(string what)
    {
        // what was selected in a layer that left the view (or was locked) leaves the selection
        var st = _doc!.AtomStates();
        if (_selection.RemoveAll(i => i < st.Length && st[i] != 0) > 0) RefreshSelection();
        RefreshViewStates();
        RefreshLayerStates();
        Raise(nameof(LayerChip));
        QueueSelBar();
        Status = what;
        RenderRequested?.Invoke();
    }

    /// <summary>A row selects its atoms in the view (⇧ adds).</summary>
    public void SelectLayer(LayerRow r, bool add)
    {
        if (_doc == null || r.IsMore) return;
        var st = _doc.AtomStates();
        var atoms = LayerAtoms(r).Where(i => i >= st.Length || st[i] == 0).ToArray();
        if (atoms.Length == 0) { Status = $"{r.Name} is hidden, ghosted or locked: show or free it to select it"; return; }
        var res = JsonNode.Parse(_doc.Select(new JsonObject
        {
            ["mode"] = "indices", ["atoms"] = new JsonArray(atoms.Select(a => (JsonNode)a).ToArray()), ["op"] = add ? "add" : "replace",
        }.ToJsonString()))!;
        if (!add) { _selection.Clear(); RefreshSelection(); }
        SelectedCount = (int)((double?)res["count"] ?? 0);
        SelectHud = "layer";
        RenderRequested?.Invoke();
    }

    /// <summary>The HUD's reminder of what is out of view: "h-BN slab hidden · molecule 5 ghosted".</summary>
    public string LayerChip
    {
        get
        {
            var parts = new List<string>();
            foreach (var k in _layerKinds)
            {
                if (k.State is "hidden" or "ghost") { parts.Add($"{k.Name} {(k.State == "hidden" ? "hidden" : "ghosted")}"); continue; }
                if (_layerChildren.TryGetValue(k, out var ch))
                    parts.AddRange(ch.Where(c => c.State is "hidden" or "ghost").Select(c => $"{c.Name} {(c.State == "hidden" ? "hidden" : "ghosted")}"));
            }
            if (parts.Count == 0) return HiddenChip;
            return parts.Count <= 2 ? string.Join(" · ", parts) : $"{parts[0]} · {parts[1]} · {parts.Count - 2} more";
        }
    }
}
