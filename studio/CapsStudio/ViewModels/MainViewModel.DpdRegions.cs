using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>A region of a DPD start: a droplet, rod, slab, shell or tube placed in the box; species name it to start
/// there, and a region no species names starts empty.</summary>
public sealed class DpdRegionRow : ObservableObject
{
    public static readonly string[] Shapes = ["droplet", "rod", "slab", "shell", "tube"];
    public static readonly string[] Axes = ["x", "y", "z"];
    private string _name = "core", _centre = "0.5 0.5 0.5";
    private int _shape, _axis = 2;
    private decimal _radius = 3, _inner = 0, _thickness = 3, _priority;
    public string Name { get => _name; set => Set(ref _name, (value ?? "").Trim()); }
    public int Shape { get => _shape; set { if (Set(ref _shape, Math.Clamp(value, 0, Shapes.Length - 1))) { Raise(nameof(HasInner)); Raise(nameof(IsSlab)); Raise(nameof(HasAxis)); Raise(nameof(IsTube)); } } }
    public string Centre { get => _centre; set => Set(ref _centre, value ?? ""); }
    public decimal Radius { get => _radius; set => Set(ref _radius, Math.Clamp(value, 0.1m, 500)); }
    public decimal Inner { get => _inner; set => Set(ref _inner, Math.Clamp(value, 0, 500)); }
    public decimal Thickness { get => _thickness; set => Set(ref _thickness, Math.Clamp(value, 0.1m, 1000)); }
    public int Axis { get => _axis; set => Set(ref _axis, Math.Clamp(value, 0, 2)); }
    public decimal Priority { get => _priority; set => Set(ref _priority, Math.Clamp(Math.Round(value), -99, 99)); }
    public bool HasInner => _shape is 3 or 4;
    public bool IsSlab => _shape == 2;
    public bool HasAxis => _shape is 1 or 2 or 4;
    public bool IsTube => _shape == 4;

    public JsonObject ToJson()
    {
        var inv = CultureInfo.InvariantCulture;
        var c = _centre.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries).Select(x => double.Parse(x, inv)).ToArray();
        if (c.Length != 3) throw new FormatException($"region {_name}: the centre takes three box fractions, e.g. 0.5 0.5 0.5");
        if (HasInner && _inner >= _radius) throw new FormatException($"region {_name}: the inner radius must be below the outer one");
        return new JsonObject
        {
            ["name"] = _name, ["shape"] = Shapes[_shape], ["centre"] = new JsonArray(c[0], c[1], c[2]), ["radius"] = (double)_radius,
            ["inner"] = (double)_inner, ["thickness"] = (double)_thickness, ["axis"] = _axis, ["priority"] = (int)_priority,
        };
    }
}

/// <summary>Mesoscale (DPD) starts composed from regions (core dpd.hpp DpdRegion): a core–shell droplet, an interphase
/// round a filler rod, a compatibiliser along an interface, a void; amphiphiles laid head on the region's surface.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] DpdOrients = ["random walk", "head on surface, tail in", "head on surface, tail out"];
    private decimal _dpdBox;
    public ObservableCollection<DpdRegionRow> DpdRegions { get; } = new();
    /// <summary>Box edge in r_c; 0 sizes it from the beads and the density.</summary>
    public decimal DpdBox { get => _dpdBox; set => Set(ref _dpdBox, Math.Clamp(value, 0, 1000)); }
    public bool DpdHasRegions => DpdRegions.Count > 0;
    public void AddDpdRegion()
    {
        DpdRegions.Add(DpdRegions.Count == 0 ? new DpdRegionRow { Name = "core", Radius = 3, Priority = 1 }
                                             : new DpdRegionRow { Name = "shell", Shape = 3, Inner = 3, Radius = 5 });
        Raise(nameof(DpdHasRegions));
    }
    public void RemoveDpdRegion(DpdRegionRow r) { DpdRegions.Remove(r); Raise(nameof(DpdHasRegions)); }
    /// <summary>The region names a species may start in: "" (outside every region) first.</summary>
    public string[] DpdRegionNames => new[] { "" }.Concat(DpdRegions.Select(r => r.Name)).ToArray();

    /// <summary>The regions, the box and each species' region and orientation into the run's JSON.</summary>
    private void AddDpdRegions(JsonObject json, JsonArray species)
    {
        if (_dpdBox > 0) json["box"] = (double)_dpdBox;
        if (DpdRegions.Count == 0) return;
        var names = DpdRegions.Select(r => r.Name).ToHashSet();
        json["regions"] = new JsonArray(DpdRegions.Select(r => (JsonNode)r.ToJson()).ToArray());
        for (var i = 0; i < DpdSpecies.Count && i < species.Count; ++i)
        {
            var row = DpdSpecies[i];
            if (row.Region.Length == 0) continue;
            if (!names.Contains(row.Region)) throw new FormatException($"{row.Name}: no region named {row.Region}");
            species[i]!["region"] = row.Region;
            if (row.Orient > 0) species[i]!["orient"] = row.Orient == 1 ? "in" : "out";
        }
    }
}
