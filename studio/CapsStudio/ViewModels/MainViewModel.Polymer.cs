using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A repeat unit or copolymer preset from data/polymers/library.json.</summary>
public sealed record LibraryEntry(string Id, string Name, string Smiles, bool Rubber, bool Copolymer, JsonObject? Preset)
{
    public string Kind => Copolymer ? "copolymer" : Rubber ? "rubber" : "";
    public bool HasKind => Kind.Length > 0;
}

/// <summary>One repeat unit of the chain being built (A, B, C …).</summary>
public sealed class PolyUnit : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public Action? Changed;

    private string _letter = "A";
    public string Letter { get => _letter; set { _letter = value; Raise(nameof(Letter)); } }
    private string _name = "", _smiles = "", _info = "", _error = "";
    private decimal _weight = 1, _block = 10;
    public string Name { get => _name; set { _name = value; Raise(nameof(Name)); Changed?.Invoke(); } }
    public string Smiles
    {
        get => _smiles;
        set { _smiles = value; Raise(nameof(Smiles)); Describe(); Changed?.Invoke(); }
    }
    public string Info { get => _info; private set { _info = value; Raise(nameof(Info)); } }
    public string Error { get => _error; private set { _error = value; Raise(nameof(Error)); Raise(nameof(HasError)); Raise(nameof(Ok)); } }
    public bool HasError => _error.Length > 0;
    public bool Ok => _error.Length == 0 && _smiles.Length > 0;
    public int Stereocentres { get; private set; }
    public decimal Weight { get => _weight; set { _weight = Math.Max(0, value); Raise(nameof(Weight)); Changed?.Invoke(); } }
    public decimal Block { get => _block; set { _block = Math.Max(1, value); Raise(nameof(Block)); Changed?.Invoke(); } }

    private void Describe()
    {
        try
        {
            var r = JsonNode.Parse(CapsDocument.UnitInfo(_smiles))!;
            if ((bool?)r["ok"] == true)
            {
                Stereocentres = (int?)r["stereocentres"] ?? 0;
                Info = string.Format(CultureInfo.InvariantCulture, "{0} · {1:F2} g/mol · head {2} → tail {3}{4}", Sub((string?)r["formula"] ?? ""),
                    (double?)r["mass"] ?? 0, r["head_element"], r["tail_element"], Stereocentres > 0 ? " · stereocentre" : "");
                Error = "";
            }
            else { Info = ""; Error = (string?)r["error"] ?? "not a repeat unit"; }
        }
        catch (Exception e) { Info = ""; Error = e.Message; }
    }

    internal static string Sub(string f) => string.Concat(f.Select(c => c is >= '0' and <= '9' ? (char)('₀' + (c - '0')) : c));
}

/// <summary>The polymer builder (design/boards/PolymerBuilder): repeat units from the library or SMILES, the sequence
/// (homopolymer, alternating, block, random, gradient, pattern) over any number of units, tacticity, DP; one chain
/// previewed; sent to Grow for a cell.</summary>
public sealed partial class MainViewModel
{
    public bool IsPolymer => _module == 13;
    public static readonly string[] SequenceKinds = ["Homopolymer", "Alternating", "Block", "Random", "Gradient", "Pattern"];
    private static readonly string[] SequenceIds = ["homopolymer", "alternating", "block", "random", "gradient", "pattern"];

    public List<LibraryEntry> PolymerLibrary { get; } = new();
    public ObservableCollection<LibraryEntry> LibraryShown { get; } = new();
    private string _libQuery = "";
    private bool _libRubber;
    public string LibraryQuery { get => _libQuery; set { if (Set(ref _libQuery, value)) FilterLibrary(); } }
    public bool LibraryRubberOnly { get => _libRubber; set { if (Set(ref _libRubber, value)) FilterLibrary(); } }
    public string LibraryCount => $"{PolymerLibrary.Count(e => !e.Copolymer)} repeat units · {PolymerLibrary.Count(e => e.Copolymer)} copolymers";

    public ObservableCollection<PolyUnit> PolyUnits { get; } = new();
    private int _polySeq;
    public int PolySequence { get => _polySeq; set { if (Set(ref _polySeq, value)) { RaisePolyKinds(); PolyChanged(); } } }
    public bool PolyShowWeights => _polySeq == 3;
    public bool PolyShowBlocks => _polySeq == 2;
    public bool PolyShowPattern => _polySeq == 5;
    private string _polyPattern = "AB";
    public string PolyPattern { get => _polyPattern; set { if (Set(ref _polyPattern, value)) PolyChanged(); } }
    private string _polyName = "Polystyrene";
    public string PolyName { get => _polyName; set { if (Set(ref _polyName, value)) PolyChanged(); } }
    private string _polyPreview = "", _polyStrip = "", _polyError = "";
    public string PolyPreview { get => _polyPreview; private set => Set(ref _polyPreview, value); }
    public string PolyError { get => _polyError; private set { if (Set(ref _polyError, value)) Raise(nameof(PolyHasError)); } }
    public bool PolyHasError => _polyError.Length > 0;
    public int[] PolyStripUnits { get; private set; } = [];
    public event Action? PolyStripChanged;
    private CapsDocument? _polyDoc;
    public CapsDocument? PolyDoc { get => _polyDoc; private set => Set(ref _polyDoc, value); }
    private bool _polyBuilding;
    public bool PolyBuilding { get => _polyBuilding; private set => Set(ref _polyBuilding, value); }
    public bool PolyCanAdd => PolyUnits.Count < 8;

    private void RaisePolyKinds() { Raise(nameof(PolyShowWeights)); Raise(nameof(PolyShowBlocks)); Raise(nameof(PolyShowPattern)); }

    public void LoadPolymerLibrary()
    {
        if (PolymerLibrary.Count > 0) return;
        try
        {
            var path = Paths.Polymers;
            if (path == null) return;
            var j = JsonNode.Parse(File.ReadAllText(path))!;
            var byId = new Dictionary<string, JsonObject>();
            foreach (var p in (j["polymers"] as JsonArray ?? []).OfType<JsonObject>())
            {
                var id = (string?)p["id"] ?? "";
                byId[id] = p;
                var tags = (p["tags"] as JsonArray)?.Select(x => (string?)x).ToList() ?? [];
                PolymerLibrary.Add(new LibraryEntry(id, (string?)p["name"] ?? id, (string?)p["smiles"] ?? "", tags.Contains("rubber"), false, null));
            }
            foreach (var c in (j["copolymers"] as JsonArray ?? []).OfType<JsonObject>())
            {
                var tags = (c["tags"] as JsonArray)?.Select(x => (string?)x).ToList() ?? [];
                PolymerLibrary.Add(new LibraryEntry((string?)c["id"] ?? "", (string?)c["name"] ?? "", "", tags.Contains("rubber"), true, c));
            }
            _libById = byId;
        }
        catch (Exception e) { Status = "Polymer library: " + e.Message; }
        FilterLibrary();
        Raise(nameof(LibraryCount));
        if (PolyUnits.Count == 0) AddPolyUnit("Styrene", "*CC(*)c1ccccc1");
    }
    private Dictionary<string, JsonObject> _libById = new();

    private void FilterLibrary()
    {
        LibraryShown.Clear();
        var q = _libQuery.Trim().ToLowerInvariant();
        foreach (var e in PolymerLibrary)
        {
            if (_libRubber && !e.Rubber) continue;
            if (q.Length > 0 && !e.Name.ToLowerInvariant().Contains(q) && !e.Smiles.ToLowerInvariant().Contains(q) && !e.Id.ToLowerInvariant().Contains(q)) continue;
            LibraryShown.Add(e);
        }
    }

    public void AddPolyUnit(string name = "", string smiles = "*CC*")
    {
        if (PolyUnits.Count >= 8) return;
        var u = new PolyUnit { Letter = ((char)('A' + PolyUnits.Count)).ToString(), Changed = PolyChanged };
        u.Name = name;
        u.Smiles = smiles;
        PolyUnits.Add(u);
        Raise(nameof(PolyCanAdd));
        PolyChanged();
    }

    public void RemovePolyUnit(PolyUnit u)
    {
        if (PolyUnits.Count <= 1) return;
        PolyUnits.Remove(u);
        for (var k = 0; k < PolyUnits.Count; ++k) PolyUnits[k].Letter = ((char)('A' + k)).ToString();
        Raise(nameof(PolyCanAdd));
        PolyChanged();
    }

    /// <summary>A library entry: a repeat unit goes into the selected slot (or A), a copolymer preset replaces all.</summary>
    public void UseLibrary(LibraryEntry e, PolyUnit? slot)
    {
        if (!e.Copolymer)
        {
            var target = slot ?? PolyUnits.FirstOrDefault();
            if (target == null) AddPolyUnit(e.Name, e.Smiles);
            else { target.Name = e.Name; target.Smiles = e.Smiles; }
            if (PolyUnits.Count == 1) _polyName = e.Name;
            Raise(nameof(PolyName));
            PolyChanged();
            return;
        }
        var p = e.Preset!;
        PolyUnits.Clear();
        var units = (p["units"] as JsonArray ?? []).Select(x => (string?)x ?? "").ToList();
        var names = (p["unit_names"] as JsonArray)?.Select(x => (string?)x ?? "").ToList();
        for (var k = 0; k < units.Count; ++k)
        {
            var u = units[k];
            var (nm, smi) = _libById.TryGetValue(u, out var lp) ? ((string?)lp["name"] ?? u, (string?)lp["smiles"] ?? "") : (names != null && k < names.Count ? names[k] : u, u);
            AddPolyUnit(nm, smi);
        }
        var seq = Array.IndexOf(SequenceIds, (string?)p["sequence"] ?? "homopolymer");
        _polySeq = Math.Max(0, seq);
        RaisePolyKinds();
        Raise(nameof(PolySequence));
        if (p["weights"] is JsonArray w)
            for (var k = 0; k < Math.Min(w.Count, PolyUnits.Count); ++k) PolyUnits[k].Weight = (decimal)((double?)w[k] ?? 1);
        if (p["blocks"] is JsonArray b)
            for (var k = 0; k < Math.Min(b.Count, PolyUnits.Count); ++k) PolyUnits[k].Block = (decimal)((double?)b[k] ?? 10);
        _polyName = e.Name;
        Raise(nameof(PolyName));
        PolyChanged();
    }

    /// <summary>The chain as the core reads it (caps_chain_preview / caps_grow_chains).</summary>
    public string PolySpecJson(int? dp = null)
    {
        var o = new JsonObject
        {
            ["units"] = new JsonArray(PolyUnits.Select(u => (JsonNode)new JsonObject { ["name"] = u.Name, ["smiles"] = u.Smiles }).ToArray()),
            ["sequence"] = SequenceIds[Math.Clamp(_polySeq, 0, 5)],
            ["dp"] = dp ?? _growDp,
            ["weights"] = new JsonArray(PolyUnits.Select(u => (JsonNode)(double)u.Weight).ToArray()),
            ["blocks"] = new JsonArray(PolyUnits.Select(u => (JsonNode)(int)u.Block).ToArray()),
            ["pattern"] = _polyPattern,
        };
        // unit templates cleaned with the default force field when it types them (GAFF2 unless Settings says otherwise)
        if (CleanChoices.FirstOrDefault(c => c.File != null && Path.GetFileNameWithoutExtension(c.File) == _settings.ForceField)?.File is { } clean) o["forcefield"] = clean;
        else if (CleanChoices.FirstOrDefault(c => c.File != null)?.File is { } first) o["forcefield"] = first;
        return o.ToJsonString();
    }

    private void PolyChanged()
    {
        if (PolyUnits.Count == 0) return;
        if (PolyUnits.FirstOrDefault(u => !u.Ok) is { } bad) { PolyError = $"{bad.Letter}: {(bad.HasError ? bad.Error : "no SMILES")}"; return; }
        try
        {
            var r = JsonNode.Parse(CapsDocument.ChainPreview(PolySpecJson(), 1))!;
            if ((bool?)r["ok"] != true) { PolyError = (string?)r["error"] ?? "cannot build this chain"; return; }
            PolyError = "";
            PolyStripUnits = (r["sequence"] as JsonArray ?? []).Select(x => (int?)x ?? 0).ToArray();
            var counts = PolyUnits.Select((u, k) => PolyStripUnits.Count(x => x == k)).ToArray();
            PolyPreview = string.Format(CultureInfo.InvariantCulture, "{0} · {1:N0} g/mol · {2:N0} atoms per chain · {3}",
                PolyUnit.Sub((string?)r["formula"] ?? ""), (double?)r["mass"] ?? 0, (int?)r["atoms"] ?? 0,
                string.Join(" · ", PolyUnits.Select((u, k) => $"{u.Letter} {counts[k]}")));
            _polyAtoms = (int?)r["atoms"] ?? 0;
            _polyMass = (double?)r["mass"] ?? 0;
            PolyStripChanged?.Invoke();
        }
        catch (Exception e) { PolyError = e.Message; }
    }
    private int _polyAtoms;
    private double _polyMass;

    /// <summary>One chain (the preview) grown on its own.</summary>
    public async Task BuildPolyPreview()
    {
        if (PolyHasError || _polyBuilding) return;
        PolyBuilding = true;
        var spec = PolySpecJson(Math.Min(_growDp, 60));
        var tact = _growTact;
        try
        {
            var (doc, _) = await Task.Run(() => CapsDocument.GrowChains(spec, new CapsGrowOpts { Chains = 1, Dp = 0, Tacticity = tact, Seed = 1, Density = 0.02, ContactScale = 0.8, Curve = 1 }, null, "chain"));
            var old = _polyDoc;
            PolyDoc = doc;
            old?.Dispose();
        }
        catch (Exception e) { PolyError = "Preview: " + e.Message; }
        finally { PolyBuilding = false; }
    }

    /// <summary>One chain of the spec opened as the Studio document.</summary>
    public async Task BuildPolymerInStudio()
    {
        if (PolyHasError || !Idle) return;
        var spec = PolySpecJson();
        var tact = _growTact;
        PolyBuilding = true;
        try
        {
            var (doc, _) = await Task.Run(() => CapsDocument.GrowChains(spec, new CapsGrowOpts { Chains = 1, Dp = 0, Tacticity = tact, Seed = 1, Density = 0.02, ContactScale = 0.8, Curve = 1 }, null, "chain"));
            Show(doc, (_polyName.Length > 0 ? _polyName : "polymer") + $" · 1 chain × {_growDp}");
            GrownUnsaved = true;
            SetModule(8);
        }
        catch (Exception e) { PolyError = e.Message; }
        finally { PolyBuilding = false; }
    }

    // ---- Grow uses the chain spec once one is sent
    private string? _growSpec;
    private string _growSpecName = "PS";
    private int _growSpecAtoms;
    private double _growSpecMass;
    public bool GrowHasSpec => _growSpec != null;
    public string GrowComponentName => _growSpec == null ? "PS " + GrowTacticityName : _growSpecName;

    public void SendPolymerToGrow()
    {
        if (PolyHasError) return;
        _growSpec = PolySpecJson();
        _growSpecName = PolyUnits.Count == 1 ? (PolyUnits[0].Name.Length > 0 ? PolyUnits[0].Name : "polymer") : (_polyName.Length > 0 ? _polyName : "copolymer");
        _growSpecAtoms = _polyAtoms;
        _growSpecMass = _polyMass;
        Raise(nameof(GrowHasSpec)); Raise(nameof(GrowComponentName)); Raise(nameof(GrowAtomsText)); Raise(nameof(GrowEstimate)); Raise(nameof(GrowCommand));
        SetModule(0);
        Status = $"Grow builds {_growSpecName} now · set chains, density and contact scale, then Grow";
    }

    public void UsePolystyreneInGrow()
    {
        _growSpec = null;
        Raise(nameof(GrowHasSpec)); Raise(nameof(GrowComponentName)); Raise(nameof(GrowAtomsText)); Raise(nameof(GrowEstimate)); Raise(nameof(GrowCommand));
    }

    /// <summary>Atoms and mass per chain for the chosen polymer at the current DP (polystyrene when none).</summary>
    private (long Atoms, double Mass) GrowChainSize()
    {
        if (_growSpec == null) return (16L * _growDp + 2, _growDp * (8 * 12.011 + 8 * 1.008) + 2 * 1.008);
        try
        {
            var spec = JsonNode.Parse(_growSpec)!.AsObject();
            spec["dp"] = _growDp;
            var r = JsonNode.Parse(CapsDocument.ChainPreview(spec.ToJsonString(), (ulong)_growSeed))!;
            return ((int?)r["atoms"] ?? 0, (double?)r["mass"] ?? 0);
        }
        catch { return (_growSpecAtoms, _growSpecMass); }
    }
}
