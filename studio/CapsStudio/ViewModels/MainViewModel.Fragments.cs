using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A fragment of the library: a SMILES with * attachment points (or a whole molecule without them).</summary>
public sealed class FragmentItem : ObservableObject
{
    public FragmentItem(string category, string name, string smiles, string meta) { Category = category; Name = name; Smiles = smiles; Meta = meta; }
    public string Category { get; }
    public string Name { get; }
    public string Smiles { get; }
    public string Meta { get; }
    public int Points => Smiles.Count(c => c == '*');
    public bool Attaches => Points > 0;
    public string Tip => Attaches ? $"{Name} · {Smiles} · click to attach to the picked atom" : $"{Name} · {Smiles} · click to place beside the structure";
    public string Short => Name.Length > 11 ? Name[..10] + "…" : Name;
    public string Icon => Category switch
    {
        "Rings" => Name.Contains("Cyclo") ? "hex" : "ring",
        "Monomers" => "grow", "Solvents" => "flask", "Ions" => "atom", "Amino acids" => "link", "Rubber additives" => "flask",
        _ => Name is "Methyl" or "Ethyl" ? "atom" : "link",
    };
    /// <summary>The SMILES with each * written as a hydrogen (a stand-in the builder can embed).</summary>
    public string SmilesH
    {
        get
        {
            var sb = new System.Text.StringBuilder();
            for (var k = 0; k < Smiles.Length; k++)
                if (Smiles[k] != '*') sb.Append(Smiles[k]);
                else sb.Append(k > 0 && Smiles[k - 1] == '[' ? "H" : "[H]");
            return sb.ToString();
        }
    }
    /// <summary>Indices (among the written atoms) of the * points: the atoms to ring in the preview.</summary>
    public int[] StarIndices
    {
        get
        {
            var list = new List<int>();
            var atom = 0;
            for (var k = 0; k < Smiles.Length; k++)
            {
                var c = Smiles[k];
                if (c == '[') { var close = Smiles.IndexOf(']', k); if (Smiles.Substring(k, close - k).Contains('*')) list.Add(atom); atom++; k = close; continue; }
                if (c == '*') { list.Add(atom++); continue; }
                if (c == '%') { k += 2; continue; }
                if (c is 'C' && k + 1 < Smiles.Length && Smiles[k + 1] == 'l') { atom++; k++; continue; }
                if (c is 'B' && k + 1 < Smiles.Length && Smiles[k + 1] == 'r') { atom++; k++; continue; }
                if ("BCNOPSFIbcnops".Contains(c)) atom++;
            }
            return list.ToArray();
        }
    }
}

/// <summary>A category of the library with its count.</summary>
public sealed record FragmentCategory(string Name, int Count)
{
    public string CountText => Count.ToString(CultureInfo.InvariantCulture);
}

/// <summary>Fragment library (design/boards/FragmentLibrary): rings, functional groups, monomers (NR, BR, CR, silicone …),
/// amino acids, bases, sugars, solvents, ions and rubber additives (cure agents, accelerators, antioxidants, silane
/// coupling agents), plus the user's own; attach to the picked atom or place beside the structure.</summary>
public sealed partial class MainViewModel
{
    public bool IsFragments => _module == 35;
    private List<FragmentItem>? _fragments;
    public ObservableCollection<FragmentCategory> FragmentCategories { get; } = new();
    public ObservableCollection<FragmentItem> FragmentTiles { get; } = new();
    public ObservableCollection<FragmentItem> QuickFragments { get; } = new();
    public event Action? FragmentChanged;

    private List<FragmentItem> Fragments
    {
        get
        {
            if (_fragments != null) return _fragments;
            _fragments = new();
            try
            {
                var path = Paths.Fragments;
                if (path != null && File.Exists(path))
                    foreach (var f in JsonNode.Parse(File.ReadAllText(path))!["fragments"]!.AsArray())
                        _fragments.Add(new FragmentItem(f!["category"]!.GetValue<string>(), f["name"]!.GetValue<string>(), f["smiles"]!.GetValue<string>(), f["meta"]!.GetValue<string>()));
            }
            catch (Exception e) { FragmentError = "Cannot read the fragment library: " + e.Message; }
            foreach (var m in _settings.MyFragments) _fragments.Add(new FragmentItem("My fragments", m.Name, m.Smiles, "yours"));
            return _fragments;
        }
    }

    /// <summary>The sidebar's tiles: common builder pieces.</summary>
    public void LoadQuickFragments()
    {
        QuickFragments.Clear();
        foreach (var n in new[] { "Benzene", "Cyclohexane", "Ester", "Amide", "Methyl", "Vinyl", "Water", "Sodium", "Styrene" })
            if (Fragments.FirstOrDefault(f => f.Name == n) is { } item) QuickFragments.Add(item);
    }

    private string _fragCategory = "Functional groups", _fragQuery = "", _fragError = "";
    public string FragmentCategory { get => _fragCategory; set { if (value != null && Set(ref _fragCategory, value)) FilterFragments(); } }
    public string FragmentQuery { get => _fragQuery; set { if (Set(ref _fragQuery, value ?? "")) FilterFragments(); } }
    public string FragmentError { get => _fragError; private set => Set(ref _fragError, value); }
    public string FragmentTitle => _fragQuery.Trim().Length > 0 ? $"“{_fragQuery.Trim()}”" : _fragCategory;
    public string FragmentCountText => FragmentTiles.Count.ToString(CultureInfo.InvariantCulture);
    public string FragmentStatus => $"{Fragments.Count:N0} fragments · {_settings.MyFragments.Count} yours";

    public void OpenFragments(string? query = null)
    {
        FragmentCategories.Clear();
        FragmentCategories.Add(new FragmentCategory("All", Fragments.Count));
        foreach (var g in Fragments.GroupBy(f => f.Category)) FragmentCategories.Add(new FragmentCategory(g.Key, g.Count()));
        if (!FragmentCategories.Any(c => c.Name == "My fragments")) FragmentCategories.Add(new FragmentCategory("My fragments", 0));
        if (query != null) _fragQuery = query;
        Raise(nameof(FragmentQuery));
        SetModule(35);
        FilterFragments();
    }

    private void FilterFragments()
    {
        FragmentTiles.Clear();
        var q = _fragQuery.Trim();
        IEnumerable<FragmentItem> hits = q.Length > 0
            ? Fragments.Where(f => f.Name.Contains(q, StringComparison.OrdinalIgnoreCase) || f.Smiles.Contains(q, StringComparison.Ordinal) || f.Meta.Contains(q, StringComparison.OrdinalIgnoreCase))
            : _fragCategory == "All" ? Fragments : Fragments.Where(f => f.Category == _fragCategory);
        foreach (var f in hits) FragmentTiles.Add(f);
        Raise(nameof(FragmentTitle));
        Raise(nameof(FragmentCountText));
        Raise(nameof(FragmentStatus));
        if (_fragSel == null || !FragmentTiles.Contains(_fragSel)) SelectedFragment = FragmentTiles.FirstOrDefault();
        FragmentChanged?.Invoke();
    }

    private FragmentItem? _fragSel;
    public FragmentItem? SelectedFragment
    {
        get => _fragSel;
        set
        {
            if (!Set(ref _fragSel, value)) return;
            foreach (var n in new[] { nameof(FragmentPointsChip), nameof(FragmentAttachText), nameof(FragmentHasPoints), nameof(FragmentDefinition) }) Raise(n);
            FragmentChanged?.Invoke();
        }
    }
    public bool FragmentHasPoints => _fragSel?.Attaches == true;
    public string FragmentPointsChip => _fragSel == null ? "" : _fragSel.Points == 0 ? "a whole molecule" : $"{_fragSel.Points} attach point{(_fragSel.Points == 1 ? "" : "s")}";
    public string FragmentAttachText => _fragSel?.Attaches == true ? "Attach to selection" : "Place beside structure";
    public string FragmentDefinition => _fragSel == null ? "" : $"{_fragSel.Smiles}\n{_fragSel.Meta}";

    private bool _fragReplaceH = true, _fragClean = true;
    public bool FragmentReplaceH { get => _fragReplaceH; set => Set(ref _fragReplaceH, value); }
    public bool FragmentClean { get => _fragClean; set => Set(ref _fragClean, value); }

    /// <summary>Attaches the fragment to the picked atom (or places a whole molecule) and returns to the Studio.</summary>
    public async Task UseFragment(FragmentItem? f = null, bool back = true)
    {
        f ??= _fragSel;
        if (f == null || _doc == null) return;
        var doc = _doc;
        string json;
        if (f.Attaches)
        {
            if (_selection.Count == 0) { FragmentError = "Pick the atom to attach to in the Studio first (click it)"; Status = FragmentError; return; }
            json = JsonSerializer.Serialize(new { op = "attach", target = _selection[0], smiles = f.Smiles, name = f.Name, replace_h = _fragReplaceH ? 1 : 0, clean = _fragClean ? 1 : 0 });
        }
        else json = JsonSerializer.Serialize(new { op = "place", smiles = f.Smiles, name = f.Name, resname = f.Name.Length >= 3 ? f.Name[..3].ToUpperInvariant() : "MOL" });
        Status = f.Attaches ? $"Attaching {f.Name}…" : $"Placing {f.Name}…";
        var text = await Task.Run(() => doc.Edit(json));
        var r = JsonNode.Parse(text)!;
        if (r["ok"]?.GetValue<bool>() != true) { FragmentError = r["error"]?.GetValue<string>() ?? "cannot attach"; Status = FragmentError; return; }
        FragmentError = "";
        AfterEdit(r["what"]!.GetValue<string>());
        if (back && IsFragments) SetModule(8);
    }

    public void SaveFragmentCopy()
    {
        if (_fragSel == null) return;
        var name = _fragSel.Category == "My fragments" ? _fragSel.Name : _fragSel.Name + " (copy)";
        _settings.MyFragments.Add(new MyFragment { Name = name, Smiles = _fragSel.Smiles });
        Changed("Fragment saved");
        _fragments = null;
        OpenFragments();
        FragmentCategory = "My fragments";
    }

    /// <summary>A fragment typed as SMILES (with * points) into My fragments.</summary>
    public void AddMyFragment(string name, string smiles)
    {
        if (smiles.Trim().Length == 0) return;
        _settings.MyFragments.Add(new MyFragment { Name = name.Trim().Length > 0 ? name.Trim() : smiles.Trim(), Smiles = smiles.Trim() });
        Changed("Fragment added");
        _fragments = null;
        OpenFragments();
        FragmentCategory = "My fragments";
    }
}
