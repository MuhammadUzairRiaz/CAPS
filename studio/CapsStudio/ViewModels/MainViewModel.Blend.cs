using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Text.Json.Nodes;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>One component of a blend: a library polymer, its weight share and chain length.</summary>
public sealed class BlendRow : INotifyPropertyChanged
{
    public event PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(n));
    public Action? Changed;
    public string Dot { get; init; } = "#F5A524";
    private LibraryEntry? _polymer;
    private decimal _weight = 50, _dp = 20, _density = 1.0m, _count = 8;
    private string _chains = "";
    /// <summary>The component's density (g/cm³), for volume shares.</summary>
    public decimal Density { get => _density; set { _density = Math.Clamp(value, 0.3m, 5m); Raise(nameof(Density)); Changed?.Invoke(); } }
    /// <summary>Chains of this component (chain-count mode).</summary>
    public decimal Count { get => _count; set { _count = Math.Clamp(Math.Round(value), 1, 10000); Raise(nameof(Count)); Changed?.Invoke(); } }
    public LibraryEntry? Polymer { get => _polymer; set { _polymer = value; Raise(nameof(Polymer)); Changed?.Invoke(); } }
    public decimal Weight { get => _weight; set { _weight = Math.Clamp(value, 0, 100); Raise(nameof(Weight)); Changed?.Invoke(); } }
    public decimal Dp { get => _dp; set { _dp = Math.Clamp(Math.Round(value), 2, 2000); Raise(nameof(Dp)); Changed?.Invoke(); } }
    public string ChainsText { get => _chains; set { _chains = value; Raise(nameof(ChainsText)); } }
    private FfEntry? _ff;
    /// <summary>The component's own force field (Field · by group after the build); the first entry: choose later.</summary>
    public FfEntry? ForceField { get => _ff; set { _ff = value; Raise(nameof(ForceField)); } }
    internal double Mass;   // g/mol of one chain
    /// <summary>A chain made in the polymer builder (its units, sequence and shares), in place of a library entry's.</summary>
    public JsonObject? Custom { get; set; }
}

/// <summary>Blend builder (design/boards/BlendBuilder): polymer blends such as NR/BR and SBR/BR tyre compounds, grown
/// component after component into one cell, mixed or as two slabs.</summary>
public sealed partial class MainViewModel
{
    public bool IsBlend => _module == 16;
    public ObservableCollection<BlendRow> BlendRows { get; } = new();
    public List<LibraryEntry> BlendLibrary { get; private set; } = new();
    private static readonly string[] BlendDots = ["#F5A524", "#5B8DEF", "#E07A5F", "#9B7BD6", "#7DC884"];

    private int _blendReturn = 13;
    /// <summary>Cancel: back where the blend was opened from (the polymer builder), not out of the builders.</summary>
    public void CancelBlend() => SetModule(_blendReturn is 16 or < 0 ? 13 : _blendReturn);

    public void OpenBlend()
    {
        if (_module != 16) _blendReturn = _module == 8 ? 13 : _module;
        LoadPolymerLibrary();
        // homopolymers, then copolymers (library and yours), then chains taken from the builder
        var built = BlendLibrary.Where(e => e.Id.StartsWith("builder:", StringComparison.Ordinal)).ToList();
        BlendLibrary = [.. PolymerLibrary.Where(p => !p.Copolymer), .. PolymerLibrary.Where(p => p.Copolymer && p.Preset != null), .. built];
        Raise(nameof(BlendLibrary));
        if (BlendRows.Count == 0)
        {
            AddBlendRow(BlendLibrary.FirstOrDefault(p => p.Name.Contains("natural rubber", StringComparison.OrdinalIgnoreCase)), 70);
            AddBlendRow(BlendLibrary.FirstOrDefault(p => p.Name.StartsWith("Cis-1,4-Polybutadiene", StringComparison.OrdinalIgnoreCase)), 30);
        }
        SetModule(16);
        BlendRecount();
    }

    private List<FfEntry>? _blendFfs;
    /// <summary>A force field per component: "choose later" first, then the library.</summary>
    public List<FfEntry> BlendForceFields => _blendFfs ??= [new FfEntry("", "Force field: choose later in Field", "", "", "", true), .. Field.Library];

    public void AddBlendRow(LibraryEntry? p = null, decimal weight = 50)
    {
        if (BlendRows.Count >= 5) return;
        var r = new BlendRow { Polymer = p ?? BlendLibrary.FirstOrDefault(), Weight = weight, Dot = BlendDots[BlendRows.Count % BlendDots.Length], ForceField = BlendForceFields[0] };
        r.Changed = BlendRecount;
        BlendRows.Add(r);
        BlendRecount();
    }
    public void RemoveBlendRow(BlendRow r)
    {
        if (BlendRows.Count <= 1) return;
        BlendRows.Remove(r);
        BlendRecount();
    }

    private decimal _blendChains = 8, _blendDensity = 0.5m;
    private int _blendMorph, _blendMode;
    /// <summary>How the composition is given: 0 weight %, 1 volume % (each component's density), 2 chain counts.</summary>
    public static readonly string[] BlendModes = ["Weight %", "Volume %", "Chain count"];
    public int BlendMode
    {
        get => _blendMode;
        set { if (Set(ref _blendMode, Math.Clamp(value, 0, 2))) { Raise(nameof(BlendByShare)); Raise(nameof(BlendByVolume)); Raise(nameof(BlendByCount)); Raise(nameof(BlendShareLabel)); BlendRecount(); } }
    }
    public bool BlendByShare => _blendMode != 2;
    public bool BlendByVolume => _blendMode == 1;
    public bool BlendByCount => _blendMode == 2;
    public string BlendShareLabel => _blendMode == 1 ? "vol %" : "wt %";
    /// <summary>A component's weight share as given (volume shares times densities).</summary>
    private static double WeightShare(BlendRow r, int mode) => mode == 1 ? (double)r.Weight * (double)r.Density : (double)r.Weight;
    public decimal BlendChains { get => _blendChains; set { if (Set(ref _blendChains, Math.Clamp(Math.Round(value), 1, 1000))) BlendRecount(); } }
    public decimal BlendDensity { get => _blendDensity; set => Set(ref _blendDensity, Math.Clamp(value, 0.1m, 1.5m)); }
    public int BlendMorph { get => _blendMorph; set { if (Set(ref _blendMorph, value)) { Raise(nameof(BlendMorphText)); } } }
    // growth as Grow's methods: best of k, Rosenbluth (soft spheres), Rosenbluth with UFF Lennard-Jones
    public static readonly string[] BlendGrowMethods = ["Best of k trials", "Rosenbluth · soft spheres", "Rosenbluth · UFF LJ (CBMC)"];
    private int _blendGrowMethod;
    private decimal _blendGrowTemp = 450;
    public int BlendGrowMethod { get => _blendGrowMethod; set { if (Set(ref _blendGrowMethod, Math.Clamp(value, 0, 2))) Raise(nameof(BlendGrowBiased)); } }
    public bool BlendGrowBiased => _blendGrowMethod > 0;
    public decimal BlendGrowTemp { get => _blendGrowTemp; set => Set(ref _blendGrowTemp, Math.Clamp(value, 100, 2000)); }
    public bool BlendMixed { get => _blendMorph == 0; set { if (value) BlendMorph = 0; } }
    public bool BlendSlabs { get => _blendMorph == 1; set { if (value) BlendMorph = 1; } }
    public bool BlendDroplet { get => _blendMorph == 2; set { if (value) BlendMorph = 2; } }
    public string BlendMorphText => _blendMorph == 1 ? "two-slab start" : _blendMorph == 2 ? "droplet start" : "mixed start";
    private string _blendSummary = "", _blendLog = "", _blendError = "";
    public string BlendSummary { get => _blendSummary; private set => Set(ref _blendSummary, value); }
    public string BlendLog { get => _blendLog; private set => Set(ref _blendLog, value); }
    public string BlendError { get => _blendError; private set { if (Set(ref _blendError, value)) Raise(nameof(BlendHasError)); } }
    public bool BlendHasError => _blendError.Length > 0;
    private bool _blendBuilding;
    public bool BlendBuilding { get => _blendBuilding; private set { if (Set(ref _blendBuilding, value)) Raise(nameof(BlendIdle)); } }
    public bool BlendIdle => !_blendBuilding;

    private readonly Dictionary<string, JsonObject> _builtChains = new();

    /// <summary>A component's chain for the core: a repeat unit as a homopolymer; a copolymer from its preset (units,
    /// sequence, unit shares, blocks, pattern), as the polymer builder reads it; a chain made in the builder as made.</summary>
    private string SpecOf(LibraryEntry p, int dp)
    {
        if (_builtChains.TryGetValue(p.Id, out var made))
        {
            var o = (JsonObject)made.DeepClone();
            o["dp"] = dp;
            return o.ToJsonString();
        }
        if (!p.Copolymer || p.Preset == null)
            return new JsonObject
            {
                ["units"] = new JsonArray(new JsonObject { ["name"] = p.Name, ["smiles"] = p.Smiles }),
                ["sequence"] = "homopolymer",
                ["dp"] = dp,
            }.ToJsonString();
        var pr = p.Preset;
        var ids = (pr["units"] as JsonArray ?? []).Select(x => (string?)x ?? "").ToList();
        var names = (pr["unit_names"] as JsonArray)?.Select(x => (string?)x ?? "").ToList();
        var units = new JsonArray();
        for (var k = 0; k < ids.Count; ++k)
        {
            var (nm, smi) = _libById.TryGetValue(ids[k], out var lp) ? ((string?)lp["name"] ?? ids[k], (string?)lp["smiles"] ?? "") : (names != null && k < names.Count ? names[k] : ids[k], ids[k]);
            units.Add(new JsonObject { ["name"] = nm, ["smiles"] = smi });
        }
        var spec = new JsonObject { ["units"] = units, ["sequence"] = (string?)pr["sequence"] ?? "random", ["dp"] = dp };
        foreach (var key in new[] { "weights", "blocks", "pattern" }) if (pr[key] is { } v) spec[key] = v.DeepClone();
        return spec.ToJsonString();
    }

    /// <summary>The polymer builder's chain as it stands (your own units, sequence and shares) becomes a blend component.</summary>
    public void AddBuilderChainToBlend()
    {
        if (PolyUnits.Count == 0 || PolyUnits.Any(u => u.Smiles.Trim().Length == 0)) { BlendError = "Make the chain in the polymer builder first (its repeat units and sequence)"; return; }
        var spec = JsonNode.Parse(PolySpecJson(20))!.AsObject();
        spec.Remove("dp");
        var name = (PolyName is { Length: > 0 } n ? n : string.Join("-", PolyUnits.Select(u => u.Name))) + " (builder)";
        var id = $"builder:{_builtChains.Count + 1}";
        _builtChains[id] = spec;
        var entry = new LibraryEntry(id, name, PolyUnits[0].Smiles, false, PolyUnits.Count > 1, null, true);
        BlendLibrary = [.. BlendLibrary, entry];
        Raise(nameof(BlendLibrary));
        if (BlendRows.Count >= 5) { BlendError = "A blend takes up to five components: remove one to add the builder's chain"; return; }
        AddBlendRow(entry, 30);
    }

    /// <summary>Chain counts from the weight shares and each chain's mass, the first component setting the scale.</summary>
    private void BlendRecount()
    {
        BlendError = "";
        foreach (var r in BlendRows)
        {
            if (r.Polymer == null) continue;
            try
            {
                var j = JsonNode.Parse(CapsDocument.ChainPreview(SpecOf(r.Polymer, (int)r.Dp), 1))!;
                r.Mass = (double?)j["mass"] ?? 0;
            }
            catch { r.Mass = 0; }
        }
        if (BlendRows.Count == 0 || BlendRows[0].Mass <= 0) return;
        var first = BlendRows[0];
        double total = 0;
        var counts = new List<int>();
        foreach (var r in BlendRows)
        {
            var n = _blendMode == 2 ? (int)r.Count
                : r == first ? (int)_blendChains
                : Math.Max(1, (int)Math.Round((double)_blendChains * first.Mass / Math.Max(1e-9, r.Mass) * WeightShare(r, _blendMode) / Math.Max(1e-9, WeightShare(first, _blendMode))));
            counts.Add(n);
            total += n * r.Mass;
        }
        double vtotal = 0;
        for (int k = 0; k < BlendRows.Count; ++k) vtotal += counts[k] * BlendRows[k].Mass / (double)BlendRows[k].Density;
        for (int k = 0; k < BlendRows.Count; ++k)
        {
            var r = BlendRows[k];
            r.ChainsText = _blendMode == 1
                ? string.Format(CultureInfo.InvariantCulture, "{0} · {1:F1} vol % · {2:F1} wt %", counts[k], 100 * counts[k] * r.Mass / (double)r.Density / vtotal, 100 * counts[k] * r.Mass / total)
                : string.Format(CultureInfo.InvariantCulture, "{0} · {1:F1} wt %", counts[k], 100 * counts[k] * r.Mass / total);
        }
        BlendSummary = string.Format(CultureInfo.InvariantCulture, "{0} components · {1} chains · {2:N0} g/mol in the cell", BlendRows.Count, counts.Sum(), total);
    }

    /// <summary>The last blend built: its document and each component's molecules.</summary>
    private (CapsDocument? Doc, List<(string Name, string Molecules, bool Crystal, string[] Elements)> Groups) _blendGroups = (null, new());

    /// <summary>Field · by group: the parts of the open structure CAPS knows (the held filler and the rest, a blend's
    /// components).</summary>
    private List<(string Name, string Molecules, bool Crystal, string[] Elements)> FieldGroupSuggestions()
    {
        var r = new List<(string, string, bool, string[])>();
        if (_doc == null) return r;
        long held = 0;
        try { held = _doc.HeldMolecule(); } catch { }
        if (held > 0)
        {
            // the filler's elements (a literature potential for them, when the library has one)
            var els = new HashSet<string>();
            try
            {
                var ids = _doc.MoleculeIds();
                for (var i = 0; i < ids.Length; ++i) if (ids[i] == held) els.Add(_doc.Atom(i).ElementSymbol);
            }
            catch { }
            r.Add(("filler", held.ToString(CultureInfo.InvariantCulture), true, els.ToArray()));
            r.Add(("matrix", "rest", false, []));
            return r;
        }
        if (ReferenceEquals(_blendGroups.Doc, _doc) && _blendGroups.Groups.Count > 0) return new(_blendGroups.Groups);
        return r;
    }

    public async Task BuildBlend()
    {
        if (BlendBuilding || BlendRows.Any(r => r.Polymer == null)) return;
        if (_blendMorph == 1 && BlendRows.Count != 2) { BlendError = "The two-slab start takes exactly two components"; return; }
        BlendBuilding = true;
        var comps = new JsonArray(BlendRows.Select(r => (JsonNode)new JsonObject
        {
            ["spec"] = JsonNode.Parse(SpecOf(r.Polymer!, (int)r.Dp)),
            ["weight"] = WeightShare(r, _blendMode),
            ["chains"] = _blendMode == 2 ? (int)r.Count : 0,
        }).ToArray());
        var opts = new JsonObject { ["components"] = comps, ["chains"] = (int)_blendChains, ["density"] = (double)_blendDensity, ["morphology"] = _blendMorph == 1 ? "slabs" : _blendMorph == 2 ? "droplet" : "mixed",
                                    ["method"] = GrowMethodIds[_blendGrowMethod], ["temperature"] = (double)_blendGrowTemp }.ToJsonString();
        var name = string.Join(" / ", BlendRows.Select(r => r.Polymer!.Name.Split(" (")[0]));
        try
        {
            var g = new CapsGrowOpts { Seed = (ulong)BlendSeed.Take(), ContactScale = 1.0, Curve = 1 };
            var (doc, rep) = await Task.Run(() => CapsDocument.GrowBlend(opts, g, (d, t, r) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => Status = $"Growing the blend · {d} of {t} chains · {r} restarts");
                return true;
            }, "blend"));
            // each component's molecules, for a force field per component (Field · by group)
            var parts = new List<(string, string, bool, string[])>();
            foreach (var line in rep.Split('\n'))
            {
                var m = System.Text.RegularExpressions.Regex.Match(line, @"^component (\d+) · molecules (\d+)-(\d+)$");
                if (m.Success && int.Parse(m.Groups[1].Value) - 1 is var k && k < BlendRows.Count)
                    parts.Add((BlendRows[k].Polymer!.Name.Split(" (")[0], $"{m.Groups[2].Value}-{m.Groups[3].Value}", false, []));
            }
            _blendGroups = (doc, parts);
            Show(doc, name + " blend");
            GrownUnsaved = true;
            BlendLog = rep;
            Status = "Blend built · compress it in Relax (target density), then equilibrate";
            // a force field per component: the groups set up in Field, assigned when every component has one
            var chosen = BlendRows.Select(r => r.ForceField is { Id.Length: > 0 } f ? f : null).ToList();
            if (chosen.Any(f => f != null) && parts.Count == BlendRows.Count)
            {
                Field.Groups.Clear();
                for (int k = 0; k < parts.Count; ++k)
                {
                    var grp = Field.AddGroup(parts[k].Item1, parts[k].Item2);
                    if (chosen[k] is { } f) grp.FfIndex = Field.Library.IndexOf(f);
                }
                Field.GroupMode = true;
                if (chosen.All(f => f != null))
                {
                    await Field.AssignGroups();
                    SetModule(7);   // Field: the typing report of every group
                    return;
                }
                Status = "Blend built · choose the remaining components' force fields in Field · by group";
            }
            SetModule(8);
        }
        catch (Exception e) { BlendError = e.Message; Status = "Could not build the blend: " + e.Message; }
        finally { BlendBuilding = false; }
    }
}
