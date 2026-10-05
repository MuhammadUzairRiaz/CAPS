using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia.Media;

namespace CapsStudio.ViewModels;

/// <summary>A probe: a point, plane, axis or ellipsoid made from atoms, drawn in the view.</summary>
public sealed class ProbeItem : ObservableObject
{
    public string Name { get; init; } = "";
    public string Kind { get; init; } = "point";
    public int[] Atoms { get; init; } = [];
    public string Colour { get; init; } = "#6CC4D8";
    public IBrush Brush => SolidColorBrush.Parse(Colour);
    public string KindText => Kind switch { "point" => "centre", "plane" => "plane", "axis" => "axis", _ => "ellipsoid" };
    public JsonObject Json => new() { ["kind"] = Kind, ["atoms"] = new JsonArray(Atoms.Select(a => (JsonNode)a).ToArray()) };
}

/// <summary>A measurement pinned as a card: its value now, and a sparkline over the trajectory's frames.</summary>
public sealed class ProbeCard : ObservableObject
{
    public string Title { get; init; } = "";
    public string Sub { get; init; } = "";
    public string Unit { get; init; } = "Å";
    public double[] Values { get; init; } = [];
    public string ValueText { get; set; } = "";
    public double[] Spark => Values.Length < 2 ? [] : Normalised(Values);
    public bool HasSpark => Values.Length > 1;
    private static double[] Normalised(double[] v)
    {
        double lo = v.Min(), hi = v.Max();
        return hi - lo < 1e-12 ? v.Select(_ => 0.5).ToArray() : v.Select(x => (x - lo) / (hi - lo)).ToArray();
    }
}

/// <summary>Probes (design/boards/Probes): from the selection, a point (the mass-weighted centre), a best-fit plane, the
/// long axis or the enclosing ellipsoid; drawn in the view; measured against each other — a centre's height above a plane,
/// an axis's angle to a plane, a plane's flatness — each measurement a card with its value and a sparkline over the frames,
/// sent to Analyze as a curve with one click.</summary>
public sealed partial class MainViewModel
{
    public ObservableCollection<ProbeItem> Probes { get; } = new();
    public ObservableCollection<ProbeCard> ProbeCards { get; } = new();
    public bool HasProbes => Probes.Count > 0;
    public bool HasProbeCards => ProbeCards.Count > 0;
    private static readonly string[] ProbeColours = ["#6CC4D8", "#F0A83C", "#9B7BD6", "#7CC784", "#E07A5F", "#E6C35C"];

    /// <summary>The selection as a probe: "point", "plane", "axis" or "ellipsoid".</summary>
    public void MakeProbe(string kind)
    {
        var a = SelectionAtoms();
        if (_doc == null || a.Length == 0) { Status = "Select the atoms the probe is made from"; return; }
        if (kind is "plane" or "axis" or "ellipsoid" && a.Length < 3) { Status = $"A {kind} needs at least three atoms"; return; }
        var what = DescribeAtoms(a);
        var p = new ProbeItem { Kind = kind, Atoms = a, Colour = ProbeColours[Probes.Count % ProbeColours.Length], Name = $"{new ProbeItem { Kind = kind }.KindText} · {what}" };
        Probes.Add(p);
        DrawProbes();
        Status = $"Probe {p.Name} · measure it against another probe from its chip";
    }

    public void RemoveProbe(ProbeItem p)
    {
        Probes.Remove(p);
        DrawProbes();
    }

    private void DrawProbes()
    {
        Raise(nameof(HasProbes));
        if (_doc == null) return;
        try
        {
            _doc.SetProbes(new JsonArray(Probes.Select(p => (JsonNode)new JsonObject
            {
                ["kind"] = p.Kind, ["atoms"] = new JsonArray(p.Atoms.Select(i => (JsonNode)i).ToArray()),
                ["rgb"] = Convert.ToInt32(p.Colour.TrimStart('#'), 16),
            }).ToArray()).ToJsonString());
        }
        catch (Exception e) { Status = "Probes: " + e.Message; }
        RenderRequested?.Invoke();
    }

    /// <summary>The measurements a probe offers against another (or alone): (label, measure, the other probe).</summary>
    public List<(string Label, string Measure, ProbeItem? Other)> ProbeMeasures(ProbeItem a)
    {
        var list = new List<(string, string, ProbeItem?)>();
        foreach (var b in Probes.Where(b => b != a))
        {
            if (a.Kind is "point" || b.Kind is "point" || (a.Kind == "plane") != (b.Kind == "plane") || a.Kind == b.Kind)
            {
                var d = b.Kind == "plane" ? $"Height above {b.Name}" : a.Kind == "plane" ? $"Height of {b.Name} above it" : $"Distance to {b.Name}";
                list.Add((d, "distance", b));
            }
            if (a.Kind != "point" && b.Kind != "point") list.Add(($"Angle to {b.Name}", "angle", b));
        }
        if (a.Kind == "plane") list.Add(("Flatness (rms from the plane)", "rms", null));
        if (a.Kind == "ellipsoid") list.Add(("Size (the longest semi-axis)", "size", null));
        return list;
    }

    /// <summary>A measurement pinned: its value now (the frame shown) and every frame's for the sparkline.</summary>
    public void PinProbeMeasure(ProbeItem a, string measure, ProbeItem? b)
    {
        if (_doc == null) return;
        try
        {
            var q = new JsonObject { ["a"] = a.Json, ["measure"] = measure };
            if (b != null) q["b"] = b.Json;
            var r = JsonNode.Parse(_doc.ProbeSeries(q.ToJsonString()))!;
            var v = (r["values"] as JsonArray ?? []).Select(x => (double?)x ?? double.NaN).ToArray();
            var unit = (string?)r["unit"] ?? "Å";
            var now = v.Length == 0 ? double.NaN : v[Math.Clamp(Frame, 0, v.Length - 1)];
            var inv = CultureInfo.InvariantCulture;
            var title = measure switch
            {
                "distance" when b?.Kind == "plane" => $"{a.Name} ↑ {b.Name}",
                "distance" => $"{a.Name} ↔ {b?.Name}",
                "angle" => $"{a.Name} ∠ {b?.Name}",
                "rms" => $"{a.Name} flatness",
                _ => $"{a.Name} size",
            };
            var sub = measure switch
            {
                "distance" when b?.Kind == "plane" => "centre → plane, along its normal",
                "distance" => "between centres, minimum image",
                "angle" => a.Kind == "plane" && b?.Kind == "plane" ? "between the planes" : (a.Kind == "plane") != (b?.Kind == "plane") ? "axis to the plane" : "between the axes",
                "rms" => "rms distance of its atoms from the plane",
                _ => "the enclosing ellipsoid's longest semi-axis",
            } + (v.Length > 1 ? $" · {v.Length:N0} frames" : "");
            ProbeCards.Add(new ProbeCard { Title = title, Sub = sub, Unit = unit, Values = v, ValueText = now.ToString(unit == "°" ? "0.0" : "0.00", inv) + (unit == "°" ? "°" : " Å") });
            Raise(nameof(HasProbeCards));
        }
        catch (Exception e) { Status = "Could not measure: " + e.Message; }
    }

    public void RemoveProbeCard(ProbeCard c) { ProbeCards.Remove(c); Raise(nameof(HasProbeCards)); }

    /// <summary>A card to Analyze as a curve over the frames.</summary>
    public void ProbeCardToAnalyze(ProbeCard c)
    {
        if (c.Values.Length == 0) return;
        var x = Enumerable.Range(0, c.Values.Length).Select(i => (double)i).ToArray();
        Analyze.Curves.Add(new SeriesItem("Probe", c.Title, "frame", $"{c.Title} ({c.Unit})", x, c.Values, false, null) { Markers = c.Values.Length < 30 });
        Analyze.CurveIndex = Analyze.Curves.Count - 1;
        SetModule(1);
        Status = $"{c.Title} in Analyze as a curve over {c.Values.Length:N0} frame{(c.Values.Length == 1 ? "" : "s")}";
    }
}
