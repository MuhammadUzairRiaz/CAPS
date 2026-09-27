using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One atom of the typing example: where it sits (head, body or tail unit, which monomer), what the rules
/// gave it and what was assigned by hand.</summary>
public sealed class UtAtom : ObservableObject
{
    public int Index { get; init; }
    public string Label { get; init; } = "";        // C12
    public string Element { get; init; } = "";
    public string Role { get; init; } = "";         // Head · A, Body · B, Tail · A
    public int RoleKind { get; init; }              // 0 head, 1 body, 2 tail
    public string Neighbours { get; init; } = "";
    private string _auto = "", _assigned = "";
    public string Auto { get => _auto; set { if (Set(ref _auto, value)) Raise(nameof(Shown)); } }
    /// <summary>Set by hand ("" = the rules' type stands).</summary>
    public string Assigned { get => _assigned; set { if (Set(ref _assigned, value)) { Raise(nameof(Shown)); Raise(nameof(IsHand)); } } }
    public string Shown => FieldNames.Short(_assigned.Length > 0 ? _assigned : _auto.Length > 0 ? _auto : "—");
    public string Type => _assigned.Length > 0 ? _assigned : _auto;
    public bool IsHand => _assigned.Length > 0;
    public bool IsUntyped => Type.Length == 0;
}

/// <summary>A type of the force field, for the manual assignment list.</summary>
public sealed record UtType(string Name, string Element, string Charge, string Description)
{
    public string Short => FieldNames.Short(Name);
    /// <summary>The description as a chemist reads it: "3-phenyl pyrrole C3", without the file's bookkeeping.</summary>
    public string Text
    {
        get
        {
            var d = Description;
            var bar = d.IndexOf('|');
            if (bar >= 0) d = d[(bar + 1)..];
            else if (d.StartsWith("moltemplate @atom:", StringComparison.Ordinal) && d.IndexOf('·') is var dot and >= 0) d = d[(dot + 1)..];
            return d.Replace("\"", "").Trim();
        }
    }
    internal (int, string) SortKey => (int.TryParse(new string(Name.TakeWhile(char.IsDigit).ToArray()), out var k) ? k : int.MaxValue, Name);
}

/// <summary>Force field › Type by hand: the repeat unit as head, body and tail of a short chain (for a copolymer, every
/// unit and every junction between two units), shown in 3D. Types come from the force field's own rules where they
/// apply; any atom can be given a type by hand — with every atom of the same chemical environment — and the result is
/// carried to the whole structure by each atom's environment (core/example_typing.hpp).</summary>
public sealed partial class MainViewModel
{
    public bool IsUnitTyping => _module == 69;
    private string _utSequence = "";
    private List<string> _utUnitNames = [];
    public ObservableCollection<UtAtom> UtAtoms { get; } = new();
    public ObservableCollection<UtType> UtTypes { get; } = new();
    private readonly List<UtType> _utAllTypes = new();
    private CapsDocument? _utDoc;
    public CapsDocument? UtDoc { get => _utDoc; private set { var old = _utDoc; if (Set(ref _utDoc, value)) old?.Dispose(); } }
    private string _utStatus = "", _utFilter = "", _utSource = "", _utFfName = "";
    public string UtStatus { get => _utStatus; private set => Set(ref _utStatus, value); }
    public string UtSource { get => _utSource; private set => Set(ref _utSource, value); }
    public string UtFfName { get => _utFfName; private set => Set(ref _utFfName, value); }
    public string UtFilter { get => _utFilter; set { if (Set(ref _utFilter, value)) FilterUtTypes(); } }
    private bool _utEquivalents = true;
    /// <summary>Assign to every atom with the same surroundings (to the learning radius), not only the picked one.</summary>
    public bool UtEquivalents { get => _utEquivalents; set => Set(ref _utEquivalents, value); }
    private UtAtom? _utSel;
    public UtAtom? UtSelected { get => _utSel; set { if (Set(ref _utSel, value)) { FilterUtTypes(); Raise(nameof(UtHighlights)); Raise(nameof(UtSelectedText)); } } }
    public string UtSelectedText => _utSel == null ? "Pick an atom in the view or the table" : $"{_utSel.Label} · {_utSel.Role} · bonded to {_utSel.Neighbours}";
    private UtType? _utTypeSel;
    public UtType? UtTypeSelected { get => _utTypeSel; set => Set(ref _utTypeSel, value); }
    public int[] UtHighlights => _utSel == null || _utDoc == null ? [] : (_utEquivalents ? _utDoc.EquivalentAtoms(_utSel.Index, 3) : [_utSel.Index]).Take(4).ToArray();
    public int UtUntyped => UtAtoms.Count(a => a.IsUntyped);
    public string UtCounts => $"{UtAtoms.Count} atoms · {UtAtoms.Count(a => a.IsHand)} by hand · {UtUntyped} untyped";
    public event Action? UtChanged;

    /// <summary>Opens the page for the force field chosen in Field and the polymer the structure was grown from (or the
    /// polymer builder's units).</summary>
    public void OpenUnitTyping()
    {
        SetModule(69);
        BuildUnitExample();
    }

    /// <summary>The example chain: head, body, tail of one unit; for several units a sequence holding every ordered pair
    /// of units once (a de Bruijn sequence), so every junction is in it.</summary>
    public void BuildUnitExample()
    {
        UtAtoms.Clear();
        _utAllTypes.Clear();
        UtTypes.Clear();
        var ff = Field.FfIndex >= 0 && Field.FfIndex < Field.Library.Count ? Field.Library[Field.FfIndex] : null;
        if (ff == null || ff.Id == "uff") { UtStatus = "Choose a force field with named types in Field first (UFF types itself)"; return; }
        UtFfName = ff.Label;
        JsonObject spec;
        try { spec = (JsonObject)JsonNode.Parse(_growSpec ?? PolySpecJson())!; }
        catch (Exception e) { UtStatus = "No polymer to build the example from: " + e.Message; return; }
        var units = spec["units"] as JsonArray;
        var n = units?.Count ?? 0;
        if (n == 0) { UtStatus = "The polymer has no repeat units"; return; }
        spec.Remove("architecture");
        spec.Remove("chain_dp");
        if (n == 1) { spec["sequence"] = "homopolymer"; spec["dp"] = 3; }
        else
        {
            var letters = string.Concat(Enumerable.Range(0, n).Select(k => (char)('A' + k)));
            var seq = DeBruijnPairs(letters);
            spec["sequence"] = "pattern";
            spec["pattern"] = seq;
            spec["dp"] = seq.Length;
        }
        var names = units!.Select(u => (string?)u?["name"] ?? "").ToList();
        // residue k (1-based, in chain order) is unit _utSequence[k − 1]
        _utSequence = n == 1 ? new string('A', 3) : (string?)spec["pattern"] ?? "";
        _utUnitNames = names;
        UtSource = n == 1 ? $"{names[0]} · head, body and tail" : $"{string.Join(" / ", names)} · {spec["dp"]} units covering every junction ({(string?)spec["pattern"]})";
        try
        {
            var o = new CapsGrowOpts { Chains = 1, Dp = 0, Tacticity = _growTact, Seed = 1, Density = 0.02, ContactScale = -0.8, Curve = 1 };   // negative: steps down where crowded
            var (doc, _) = CapsDocument.GrowChains(spec.ToJsonString(), o, null, "example");
            UtDoc = doc;
        }
        catch (Exception e) { UtStatus = "Could not build the example: " + e.Message; return; }
        AutoTypeUnitExample();
    }

    /// <summary>Every ordered pair of letters once, as consecutive units (de Bruijn B(k, 2), closed back to its start).</summary>
    internal static string DeBruijnPairs(string letters)
    {
        var k = letters.Length;
        var a = new int[2 * k];
        var seq = new List<int>();
        void Db(int t, int p)
        {
            if (t > 2) { if (2 % p == 0) for (var j = 1; j <= p; j++) seq.Add(a[j]); }
            else
            {
                a[t] = a[t - p];
                Db(t + 1, p);
                for (var j = a[t - p] + 1; j < k; j++) { a[t] = j; Db(t + 1, t); }
            }
        }
        Db(1, 1);
        seq.Add(seq[0]);
        return string.Concat(seq.Select(x => letters[x]));
    }

    /// <summary>The force field's own rules on the example (what they cannot type stays for the hand).</summary>
    public void AutoTypeUnitExample()
    {
        if (_utDoc == null) return;
        var ff = Field.Library[Field.FfIndex];
        try { _utDoc.FieldAssign(ff.File, null, 0); }
        catch (Exception e) { UtStatus = "The rules could not type the example (" + e.Message + "): give the types by hand"; }
        LoadUnitRows();
    }

    private void LoadUnitRows()
    {
        if (_utDoc == null) return;
        var hand = UtAtoms.Where(a => a.IsHand).ToDictionary(a => a.Index, a => a.Assigned);
        UtAtoms.Clear();
        var rep = _utDoc.FieldReport();
        var types = new List<string>();
        if (rep.Length > 0 && JsonNode.Parse(rep) is JsonObject r)
        {
            if (r["atoms"] is JsonArray aa) foreach (var a in aa) types.Add((string?)a?["type"] ?? "");
            if (_utAllTypes.Count == 0 && r["fftypes"] is JsonArray ft)
                foreach (var t in ft.OfType<JsonObject>())
                    _utAllTypes.Add(new UtType((string?)t["name"] ?? "", (string?)t["el"] ?? "",
                                               t["q"] is JsonValue q ? ((double)q).ToString("+0.000;−0.000", CultureInfo.InvariantCulture) : "", (string?)t["desc"] ?? ""));
        }
        var s = _utDoc.Summary();
        var lastRes = 1L;
        var props = new List<JsonObject>();
        for (var i = 0; i < s.Atoms; i++)
        {
            var p = JsonNode.Parse(_utDoc.AtomProperties(i)) as JsonObject ?? new JsonObject();
            props.Add(p);
            lastRes = Math.Max(lastRes, (long?)p["residue"] ?? 1);
        }
        for (var i = 0; i < props.Count; i++)
        {
            var p = props[i];
            var res = (long?)p["residue"] ?? 1;
            var kind = res <= 1 ? 0 : res >= lastRes ? 2 : 1;
            var unit = (string?)p["resname"] ?? "";
            var r0 = (int)res - 1;
            if (r0 >= 0 && r0 < _utSequence.Length && _utSequence[r0] - 'A' is var u && u >= 0 && u < _utUnitNames.Count)
                unit = _utUnitNames.Count > 1 ? $"{_utSequence[r0]} {_utUnitNames[u]}" : _utUnitNames[u];
            var nb = p["neighbours"] is JsonArray na ? string.Join(" ", na.Select(x => (string?)x?["element"] + ((int?)x?["index"] + 1))) : "";
            var el = (string?)p["element"] ?? "";
            UtAtoms.Add(new UtAtom
            {
                Index = i, Label = el + (i + 1), Element = el, Role = new[] { "Head", "Body", "Tail" }[kind] + (unit.Length > 0 ? " · " + unit : ""), RoleKind = kind,
                Neighbours = nb, Auto = i < types.Count ? types[i] : "", Assigned = hand.GetValueOrDefault(i, ""),
            });
        }
        FilterUtTypes();
        RaiseUt();
        UtStatus = UtUntyped == 0 ? "Every atom of the example has a type: change any by hand, then apply to the structure"
                                  : $"{UtUntyped} atoms of the example have no type from the rules: pick each and give it one";
    }

    private void RaiseUt() { Raise(nameof(UtCounts)); Raise(nameof(UtUntyped)); Raise(nameof(UtHighlights)); UtChanged?.Invoke(); }

    private void FilterUtTypes()
    {
        UtTypes.Clear();
        var el = _utSel?.Element ?? "";
        var q = _utFilter.Trim();
        var used = UtAtoms.Select(a => a.Type).Where(t => t.Length > 0).ToHashSet();
        // the types the example already carries first (the unit's chemistry), then the rest of the force field by number
        foreach (var t in _utAllTypes.OrderBy(t => used.Contains(t.Name) ? 0 : 1).ThenBy(t => t.SortKey))
        {
            if (el.Length > 0 && !string.Equals(t.Element, el, StringComparison.OrdinalIgnoreCase)) continue;
            if (q.Length > 0 && !(t.Name.Contains(q, StringComparison.OrdinalIgnoreCase) || t.Text.Contains(q, StringComparison.OrdinalIgnoreCase))) continue;
            UtTypes.Add(t);
            if (UtTypes.Count >= 400) break;
        }
    }

    /// <summary>The atom picked in the example's 3D view.</summary>
    public void PickUnitAtom(int index) => UtSelected = UtAtoms.FirstOrDefault(a => a.Index == index);

    /// <summary>The chosen type on the picked atom (and every atom of the same environment).</summary>
    public void AssignUnitType()
    {
        if (_utSel == null || _utTypeSel == null || _utDoc == null) { UtStatus = "Pick an atom and a type"; return; }
        var targets = _utEquivalents ? _utDoc.EquivalentAtoms(_utSel.Index, 3) : [_utSel.Index];
        foreach (var i in targets)
            if (UtAtoms.FirstOrDefault(a => a.Index == i) is { } a && a.Element == _utSel.Element) a.Assigned = _utTypeSel.Name;
        RaiseUt();
        UtStatus = $"{FieldNames.Short(_utTypeSel.Name)} on {targets.Length} atom{(targets.Length == 1 ? "" : "s")} of the example ({_utSel.Role})";
    }

    public void ClearUnitType()
    {
        if (_utSel == null || _utDoc == null) return;
        var targets = _utEquivalents ? _utDoc.EquivalentAtoms(_utSel.Index, 3) : [_utSel.Index];
        foreach (var i in targets) if (UtAtoms.FirstOrDefault(a => a.Index == i) is { } a) a.Assigned = "";
        RaiseUt();
    }

    /// <summary>The example's types (rules + hand) carried to every atom of the structure with the same environment.</summary>
    public void ApplyUnitTyping()
    {
        if (_doc == null || _utDoc == null) { UtStatus = "Open or grow the structure to type first"; return; }
        if (UtUntyped > 0) { UtStatus = $"{UtUntyped} atoms of the example have no type yet"; return; }
        var ff = Field.Library[Field.FfIndex];
        try
        {
            // the structure in this force field (its rules first; what they leave is filled from the example)
            try { _doc.FieldAssign(ff.File, null, 0); } catch { /* untyped atoms are fine here: the example types them */ }
            var types = new JsonArray(UtAtoms.OrderBy(a => a.Index).Select(a => (JsonNode)a.Type).ToArray());
            var (complete, report) = _doc.FieldTypeByExample(_utDoc, types.ToJsonString());
            var r = JsonNode.Parse(report) as JsonObject;
            var set = (int?)(double?)r?["set"] ?? 0;
            var unmatched = (int?)(double?)r?["unmatched"] ?? 0;
            var conflicts = r?["conflicts"] as JsonArray;
            Field.LoadReport(_doc);
            RenderRequested?.Invoke();
            UtStatus = $"{set:N0} atoms of the structure typed from the example (environments to {(int?)(double?)r?["radius"] ?? 3} bonds)" +
                       (unmatched > 0 ? $" · {unmatched:N0} atoms have surroundings the example lacks (branch points, crosslinks): type them in Field" : "") +
                       (conflicts is { Count: > 0 } ? " · " + (string?)conflicts[0] : "") +
                       (complete ? " · every parameter found" : " · some parameters are missing: see Field");
            Status = UtStatus;
        }
        catch (Exception e) { UtStatus = "Could not apply: " + e.Message; }
    }
}
