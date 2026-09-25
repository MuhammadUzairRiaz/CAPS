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
    private decimal _weight = 50, _dp = 20;
    private string _chains = "";
    public LibraryEntry? Polymer { get => _polymer; set { _polymer = value; Raise(nameof(Polymer)); Changed?.Invoke(); } }
    public decimal Weight { get => _weight; set { _weight = Math.Clamp(value, 0, 100); Raise(nameof(Weight)); Changed?.Invoke(); } }
    public decimal Dp { get => _dp; set { _dp = Math.Clamp(Math.Round(value), 2, 2000); Raise(nameof(Dp)); Changed?.Invoke(); } }
    public string ChainsText { get => _chains; set { _chains = value; Raise(nameof(ChainsText)); } }
    internal double Mass;   // g/mol of one chain
}

/// <summary>Blend builder (design/boards/BlendBuilder): polymer blends such as NR/BR and SBR/BR tyre compounds, grown
/// component after component into one cell, mixed or as two slabs.</summary>
public sealed partial class MainViewModel
{
    public bool IsBlend => _module == 16;
    public ObservableCollection<BlendRow> BlendRows { get; } = new();
    public List<LibraryEntry> BlendLibrary { get; private set; } = new();
    private static readonly string[] BlendDots = ["#F5A524", "#5B8DEF", "#E07A5F", "#9B7BD6", "#7DC884"];

    public void OpenBlend()
    {
        LoadPolymerLibrary();
        if (BlendLibrary.Count == 0)
        {
            BlendLibrary = PolymerLibrary.Where(p => !p.Copolymer).ToList();
            Raise(nameof(BlendLibrary));
        }
        if (BlendRows.Count == 0)
        {
            AddBlendRow(BlendLibrary.FirstOrDefault(p => p.Name.Contains("natural rubber", StringComparison.OrdinalIgnoreCase)), 70);
            AddBlendRow(BlendLibrary.FirstOrDefault(p => p.Name.StartsWith("Cis-1,4-Polybutadiene", StringComparison.OrdinalIgnoreCase)), 30);
        }
        SetModule(16);
        BlendRecount();
    }

    public void AddBlendRow(LibraryEntry? p = null, decimal weight = 50)
    {
        if (BlendRows.Count >= 5) return;
        var r = new BlendRow { Polymer = p ?? BlendLibrary.FirstOrDefault(), Weight = weight, Dot = BlendDots[BlendRows.Count % BlendDots.Length] };
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
    private int _blendMorph;
    public decimal BlendChains { get => _blendChains; set { if (Set(ref _blendChains, Math.Clamp(Math.Round(value), 1, 1000))) BlendRecount(); } }
    public decimal BlendDensity { get => _blendDensity; set => Set(ref _blendDensity, Math.Clamp(value, 0.1m, 1.5m)); }
    public int BlendMorph { get => _blendMorph; set { if (Set(ref _blendMorph, value)) { Raise(nameof(BlendMorphText)); } } }
    public bool BlendMixed { get => _blendMorph == 0; set { if (value) BlendMorph = 0; } }
    public bool BlendSlabs { get => _blendMorph == 1; set { if (value) BlendMorph = 1; } }
    public string BlendMorphText => _blendMorph == 1 ? "two-slab start" : "mixed start";
    private string _blendSummary = "", _blendLog = "", _blendError = "";
    public string BlendSummary { get => _blendSummary; private set => Set(ref _blendSummary, value); }
    public string BlendLog { get => _blendLog; private set => Set(ref _blendLog, value); }
    public string BlendError { get => _blendError; private set { if (Set(ref _blendError, value)) Raise(nameof(BlendHasError)); } }
    public bool BlendHasError => _blendError.Length > 0;
    private bool _blendBuilding;
    public bool BlendBuilding { get => _blendBuilding; private set { if (Set(ref _blendBuilding, value)) Raise(nameof(BlendIdle)); } }
    public bool BlendIdle => !_blendBuilding;

    private static string SpecOf(LibraryEntry p, int dp) => new JsonObject
    {
        ["units"] = new JsonArray(new JsonObject { ["name"] = p.Name, ["smiles"] = p.Smiles }),
        ["sequence"] = "homopolymer",
        ["dp"] = dp,
    }.ToJsonString();

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
            var n = r == first ? (int)_blendChains
                : Math.Max(1, (int)Math.Round((double)_blendChains * first.Mass / Math.Max(1e-9, r.Mass) * (double)r.Weight / Math.Max(1e-9, (double)first.Weight)));
            counts.Add(n);
            total += n * r.Mass;
        }
        for (int k = 0; k < BlendRows.Count; ++k)
        {
            var r = BlendRows[k];
            r.ChainsText = string.Format(CultureInfo.InvariantCulture, "{0} · {1:F1} wt %", counts[k], 100 * counts[k] * r.Mass / total);
        }
        BlendSummary = string.Format(CultureInfo.InvariantCulture, "{0} components · {1} chains · {2:N0} g/mol in the cell", BlendRows.Count, counts.Sum(), total);
    }

    public async Task BuildBlend()
    {
        if (BlendBuilding || BlendRows.Any(r => r.Polymer == null)) return;
        if (_blendMorph == 1 && BlendRows.Count != 2) { BlendError = "The two-slab start takes exactly two components"; return; }
        BlendBuilding = true;
        var comps = new JsonArray(BlendRows.Select(r => (JsonNode)new JsonObject
        {
            ["spec"] = JsonNode.Parse(SpecOf(r.Polymer!, (int)r.Dp)),
            ["weight"] = (double)r.Weight,
        }).ToArray());
        var opts = new JsonObject { ["components"] = comps, ["chains"] = (int)_blendChains, ["density"] = (double)_blendDensity, ["morphology"] = _blendMorph == 1 ? "slabs" : "mixed" }.ToJsonString();
        var name = string.Join(" / ", BlendRows.Select(r => r.Polymer!.Name.Split(" (")[0]));
        try
        {
            var g = new CapsGrowOpts { Seed = 1, ContactScale = 1.0, Curve = 1 };
            var (doc, rep) = await Task.Run(() => CapsDocument.GrowBlend(opts, g, (d, t, r) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => Status = $"Growing the blend · {d} of {t} chains · {r} restarts");
                return true;
            }, "blend"));
            Show(doc, name + " blend");
            GrownUnsaved = true;
            BlendLog = rep;
            Status = "Blend built · compress it in Relax (target density), then equilibrate";
            SetModule(8);
        }
        catch (Exception e) { BlendError = e.Message; Status = "Could not build the blend: " + e.Message; }
        finally { BlendBuilding = false; }
    }
}
