using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>One atom property as a histogram: every atom's bin, the selected atoms' bins, and the brushed range.</summary>
public sealed class BrushHist : ObservableObject
{
    public string Key { get; init; } = "";        // the column: z, charge, bonds, hybrid, distance:tag:NAME …
    public string Title { get; init; } = "";
    public string Unit { get; init; } = "";
    public bool Categorical { get; init; }
    public string[] Categories { get; init; } = [];   // categorical: the bars' names (hybrid: sp, sp2, sp3, other)
    public double[] Values = [];
    public int[] Bins = [];
    public int[] Selected = [];
    public double Lo, Hi;
    public event Action? Changed;
    /// <summary>Asked again: the selection follows the brushes.</summary>
    public Action? Brushed;

    public double? BrushFrom { get; private set; }
    public double? BrushTo { get; private set; }
    public bool HasBrush => BrushFrom != null;
    public string RangeText => !HasBrush ? Categorical ? string.Join(" · ", Categories) : $"{F(Lo)} – {F(Hi)}" + (Unit.Length > 0 ? " " + Unit : "")
        : Categorical ? string.Join(" · ", Categories.Where((_, k) => k >= BrushFrom && k <= BrushTo))
        : $"{F(BrushFrom!.Value)} – {F(BrushTo!.Value)}" + (Unit.Length > 0 ? " " + Unit : "");
    public string CountText { get; private set; } = "";
    private static string F(double v) => v.ToString(Math.Abs(v) >= 100 ? "0" : Math.Abs(v) >= 10 ? "0.0" : "0.0##", CultureInfo.InvariantCulture);

    public double FractionOf(double v, bool upper) => Categorical ? (v + (upper ? 1 : 0)) / Math.Max(1, Bins.Length) : Hi > Lo ? Math.Clamp((v - Lo) / (Hi - Lo), 0, 1) : 0;
    public double ValueAt(double f) => Lo + f * (Hi - Lo);

    public void SetBrush(double? from, double? to)
    {
        BrushFrom = from; BrushTo = to;
        Raise(nameof(RangeText)); Raise(nameof(HasBrush));
        Changed?.Invoke();
        Brushed?.Invoke();
    }
    public void BrushCategories(double f0, double f1)
    {
        var n = Math.Max(1, Bins.Length);
        SetBrush(Math.Clamp(Math.Floor(f0 * n), 0, n - 1), Math.Clamp(Math.Floor(Math.Min(f1, 0.9999) * n), 0, n - 1));
    }

    /// <summary>Bins from the values (numeric: 40 across the range; categorical: one per category).</summary>
    public void Fill(double[] values, bool[]? selected)
    {
        Values = values;
        var finite = values.Where(double.IsFinite).Where(v => v < 1e299).ToArray();
        if (Categorical) { Lo = 0; Hi = Categories.Length; }
        else { Lo = finite.Length > 0 ? finite.Min() : 0; Hi = finite.Length > 0 ? finite.Max() : 1; if (Hi - Lo < 1e-9) { Lo -= 0.5; Hi += 0.5; } }
        var n = Categorical ? Categories.Length : 40;
        Bins = new int[n];
        Selected = new int[n];
        for (var i = 0; i < values.Length; ++i)
        {
            var v = values[i];
            if (!double.IsFinite(v) || v >= 1e299) continue;
            var k = Categorical ? CategoryOf(v) : (int)Math.Clamp(Math.Floor((v - Lo) / (Hi - Lo) * n), 0, n - 1);
            if (k < 0) continue;
            ++Bins[k];
            if (selected != null && i < selected.Length && selected[i]) ++Selected[k];
        }
        CountText = selected == null ? "" : $"{selected.Count(x => x):N0} selected";
        Raise(nameof(RangeText)); Raise(nameof(CountText));
        Changed?.Invoke();
    }
    // hybrid: 1 sp, 2 sp2, 3 sp3 → bars 0, 1, 2; 0 (none) → bar 3
    private int CategoryOf(double v) => Key == "hybrid" ? (v >= 1 && v <= 3 ? (int)v - 1 : 3) : (int)v;

    /// <summary>The brush as query text: "z 7..11", "(hybrid sp2 or hybrid sp3)", "within 6 of tag filler and not within 2 of tag filler".</summary>
    public string? QueryText()
    {
        if (BrushFrom is not { } a || BrushTo is not { } b) return null;
        var inv = CultureInfo.InvariantCulture;
        string R(double v) => Math.Round(v, 3).ToString("0.###", inv);
        if (Key == "hybrid")
        {
            var parts = Enumerable.Range((int)a, (int)b - (int)a + 1).Where(k => k <= 2).Select(k => "hybrid " + new[] { "sp", "sp2", "sp3" }[k]).ToList();
            return parts.Count == 0 ? "none" : parts.Count == 1 ? parts[0] : "(" + string.Join(" or ", parts) + ")";
        }
        if (Key.StartsWith("distance:", StringComparison.Ordinal))
        {
            var of = Key[9..].StartsWith("tag:", StringComparison.Ordinal) ? $"tag \"{Key[13..]}\"" : $"index {int.Parse(Key[9..], inv) + 1}";
            return a <= Lo + 1e-9 ? $"within {R(b)} of {of}" : $"(within {R(b)} of {of} and not within {R(a)} of {of})";
        }
        if (Key == "bonds") return $"bonds {Math.Ceiling(a - 1e-9).ToString(inv)}..{Math.Floor(b + 1e-9).ToString(inv)}";
        return $"{Key} {R(a)}..{R(b)}";
    }
}

/// <summary>Brush to select (design/boards/BrushSelect): each property of the atoms as a histogram; dragging across one
/// selects that range, across several intersects them; the view rings the match live and the line under the histograms
/// is the same selection as query text — editable, and reusable in the query bar, macros and recipes.</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<BrushHist> BrushHists { get; } = new();
    private bool _brushOpen;
    public bool BrushOpen
    {
        get => _brushOpen;
        set
        {
            if (!Set(ref _brushOpen, value)) return;
            if (value) { AppearanceOpen = false; SelectionOpen = false; InteractionsOpen = false; HistoryOpen = false; StatesOpen = false; LensOpen = false; LodOpen = false; OpenBrushDefaults(); }
            Raise(nameof(ShowBrushPanel));
            Raise(nameof(ShowStudioTabs));
        }
    }
    /// <summary>In the inspector's place, while no other panel holds it.</summary>
    public bool ShowBrushPanel => IsStudio && _brushOpen && _doc != null && !_appOpen && !_selOpen && !_ixOpen && !_lodOpen && !_histOpen && !_lensOpen && !_stOpen;
    private string _brushQuery = "", _brushCount = "", _brushError = "";
    public string BrushQuery { get => _brushQuery; set => Set(ref _brushQuery, value ?? ""); }
    public string BrushCount { get => _brushCount; private set => Set(ref _brushCount, value); }
    public string BrushError { get => _brushError; private set { if (Set(ref _brushError, value)) Raise(nameof(HasBrushError)); } }
    public bool HasBrushError => _brushError.Length > 0;

    /// <summary>The properties Add more offers (Distance to… adds the tags and the picked atom).</summary>
    public static readonly (string Key, string Title, string Unit)[] BrushChoices =
    [
        ("z", "Height z", "Å"), ("x", "Position x", "Å"), ("y", "Position y", "Å"), ("charge", "Partial charge", "e"),
        ("bonds", "Bonds per atom", ""), ("hybrid", "Hybridisation", ""), ("mass", "Mass", "g/mol"), ("molecule", "Molecule", ""),
    ];

    private void OpenBrushDefaults()
    {
        if (BrushHists.Count == 0) foreach (var k in new[] { "z", "charge", "bonds", "hybrid" }) AddBrushHist(k, refill: false);
        RefillBrush();
    }

    public void OpenBrush() { SetModule(8); BrushOpen = true; }

    public void AddBrushHist(string key, bool refill = true)
    {
        if (BrushHists.Any(h => h.Key == key)) return;
        var c = BrushChoices.FirstOrDefault(x => x.Key == key);
        string title = c.Title ?? key, unit = c.Unit ?? "";
        if (key.StartsWith("distance:tag:", StringComparison.Ordinal)) { title = "Distance to " + key[13..]; unit = "Å"; }
        else if (key.StartsWith("distance:", StringComparison.Ordinal)) { title = "Distance to atom " + (int.Parse(key[9..], CultureInfo.InvariantCulture) + 1); unit = "Å"; }
        var h = new BrushHist
        {
            Key = key, Title = title, Unit = unit, Categorical = key == "hybrid", Categories = key == "hybrid" ? ["sp", "sp2", "sp3", "other"] : [],
        };
        h.Brushed = ApplyBrushes;
        BrushHists.Add(h);
        if (refill) RefillBrush();
    }

    public void RemoveBrushHist(BrushHist h)
    {
        BrushHists.Remove(h);
        ApplyBrushes();
    }

    /// <summary>Distance to…: the picked atom (the first selected), or a tag.</summary>
    public void AddBrushDistance(string? tag)
    {
        if (tag != null) { AddBrushHist("distance:tag:" + tag); return; }
        var a = SelectionAtoms();
        if (a.Length == 0) { BrushError = "Pick an atom first (or choose a tag): distances are measured from it"; return; }
        AddBrushHist("distance:" + a[0].ToString(CultureInfo.InvariantCulture));
    }

    /// <summary>Every histogram filled again from the structure, with the selection's share.</summary>
    public void RefillBrush()
    {
        if (!_brushOpen) return;
        BrushError = "";
        if (_doc == null) { foreach (var h in BrushHists) h.Fill([], null); BrushCount = ""; return; }
        bool[]? sel = null;
        var picked = SelectionAtoms();
        var n = 0;
        try
        {
            foreach (var h in BrushHists)
            {
                var v = _doc.AtomColumn(h.Key);
                n = v.Length;
                if (sel == null) { sel = new bool[n]; foreach (var i in picked) if (i >= 0 && i < n) sel[i] = true; }
                h.Fill(v, sel);
            }
        }
        catch (Exception e) { BrushError = e.Message; }
        BrushCount = $"{picked.Length:N0} of {(n > 0 ? n : _doc.Summary().Atoms):N0}";
    }

    /// <summary>The brushes intersected: the query text, and the selection it makes.</summary>
    private void ApplyBrushes()
    {
        var parts = BrushHists.Select(h => h.QueryText()).Where(q => q != null).ToList();
        BrushQuery = string.Join(" and ", parts);
        if (parts.Count == 0) { ClearAllSelection(); RefillBrush(); RenderRequested?.Invoke(); return; }
        RunBrushQuery();
    }

    /// <summary>The query line run (as typed or as brushed): the selection replaced.</summary>
    public void RunBrushQuery()
    {
        if (_doc == null || _brushQuery.Trim().Length == 0) return;
        try
        {
            var r = JsonNode.Parse(_doc.Select(new JsonObject { ["mode"] = "query", ["pattern"] = _brushQuery, ["op"] = "replace" }.ToJsonString()))!;
            if ((bool?)r["ok"] == false) { BrushError = (string?)r["error"] ?? "the query did not run"; return; }
            _selection.Clear();
            RefreshSelection();
            SelectedCount = (int)((double?)r["count"] ?? 0);
            SelectHud = "brush";
            BrushError = "";
        }
        catch (Exception e) { BrushError = e.Message; return; }
        RefreshSelBar();
        RefillBrush();
        RenderRequested?.Invoke();
    }

    /// <summary>The brushes cleared (the selection with them).</summary>
    public void ClearBrushes()
    {
        foreach (var h in BrushHists) if (h.HasBrush) { h.Brushed = null; h.SetBrush(null, null); h.Brushed = ApplyBrushes; }
        BrushQuery = "";
        ClearAllSelection();
        RefillBrush();
        RenderRequested?.Invoke();
    }
}
