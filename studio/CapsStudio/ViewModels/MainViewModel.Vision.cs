using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;
using Avalonia.Media;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

public sealed record VisionSwatch(string Label, IBrush Fill, IBrush Ink);
public sealed record VisionRow(string Vision, VisionSwatch[] Swatches);
public sealed record VisionPalette(string Name, string Count, VisionRow[] Rows);
public sealed record VisionPairRow(string Palette, string Vision, string Pair, string De);

/// <summary>Settings › Colour vision (design/boards/ColourVision): the element, chain and status colours as seen with
/// protanopia, deuteranopia and tritanopia (Machado et al. 2009), the pairs that fall below ΔE*ab 12, and a preview of
/// the 3D view with each deficiency.</summary>
public partial class MainViewModel
{
    private static readonly (string Label, int Z)[] VisionElements = [("C", 6), ("H", 1), ("O", 8), ("N", 7), ("Si", 14), ("S", 16), ("Na", 11), ("Cl", 17)];
    public ObservableCollection<VisionPalette> VisionPalettes { get; } = new();
    public ObservableCollection<VisionPairRow> VisionPairs { get; } = new();
    private string _visionSummary = "";
    public string VisionSummary { get => _visionSummary; private set => Set(ref _visionSummary, value); }
    public bool IsColourVision => _module == 46;
    public bool HasVisionPairs => VisionPairs.Count > 0;
    private int _visionPreview;
    /// <summary>The view as seen with a deficiency (0 normal, 1 protanopia, 2 deuteranopia, 3 tritanopia); session only.</summary>
    public int VisionPreview
    {
        get => _visionPreview;
        set
        {
            if (!Set(ref _visionPreview, Math.Clamp(value, 0, 3))) return;
            try { _doc?.SetVision(_visionPreview); } catch { /* closed */ }
            RenderRequested?.Invoke();
            Status = _visionPreview == 0 ? "View: normal colour vision" : $"View simulated as {VisionNames[_visionPreview].ToLowerInvariant()} (Machado et al. 2009) · exports are not affected";
        }
    }
    private static readonly string[] VisionNames = ["Normal", "Protanopia", "Deuteranopia", "Tritanopia"];

    public void OpenColourVision()
    {
        CheckVision();
        SetModule(46);
    }

    /// <summary>Recomputes the swatches and pairs for the current palette and theme.</summary>
    public void CheckVision()
    {
        static string Hex(uint c) => "#" + (c & 0xFFFFFF).ToString("X6", CultureInfo.InvariantCulture);
        string Tok(string key) => Tokens.Brush(key) is ISolidColorBrush b ? $"#{b.Color.R:X2}{b.Color.G:X2}{b.Color.B:X2}" : "#808080";
        var input = new JsonObject
        {
            ["Elements"] = new JsonObject { ["labels"] = new JsonArray(VisionElements.Select(e => (JsonNode)e.Label).ToArray()),
                                            ["colours"] = new JsonArray(VisionElements.Select(e => (JsonNode)Hex(Native.ElementColour(e.Z))).ToArray()) },
            ["Chains"] = new JsonObject { ["labels"] = new JsonArray(Enumerable.Range(1, 8).Select(k => (JsonNode)k.ToString(CultureInfo.InvariantCulture)).ToArray()),
                                          ["colours"] = new JsonArray(Enumerable.Range(0, 8).Select(k => (JsonNode)Hex(Native.CategoryColour(k))).ToArray()) },
            ["Status"] = new JsonObject { ["labels"] = new JsonArray("ok", "running", "error", "select"),
                                          ["colours"] = new JsonArray(Tok("OkB"), Tok("AccB"), Tok("ErrB"), Tok("SelB")) },
        };
        JsonNode r;
        try { r = JsonNode.Parse(CapsDocument.VisionCheck(input.ToJsonString(), 12.0))!; }
        catch (Exception e) { VisionSummary = "Could not check the palettes: " + e.Message; return; }
        VisionPalettes.Clear();
        foreach (var p in (JsonArray)r["palettes"]!)
        {
            var labels = ((JsonArray)p!["labels"]!).Select(l => (string)l!).ToArray();
            var rows = new[] { ("Normal", "normal"), ("Protanopia", "protanopia"), ("Deuteranopia", "deuteranopia"), ("Tritanopia", "tritanopia") }
                .Select(v => new VisionRow(v.Item1, ((JsonArray)p[v.Item2]!).Select((c, i) =>
                {
                    var col = Color.Parse((string)c!);
                    var light = 0.2126 * col.R + 0.7152 * col.G + 0.0722 * col.B > 140;
                    return new VisionSwatch(labels[i], new SolidColorBrush(col), new SolidColorBrush(light ? Color.Parse("#16191C") : Colors.White));
                }).ToArray())).ToArray();
            var name = (string)p["name"]! == "Chains" ? "Chains (Grow)" : (string)p["name"]!;
            VisionPalettes.Add(new VisionPalette(name, $"{labels.Length} colours", rows));
        }
        VisionPairs.Clear();
        var pairs = (JsonArray)r["pairs"]!;
        foreach (var q in pairs.Take(12))
        {
            var pal = (string)q!["palette"]! == "Chains" ? "Chains (Grow)" : (string)q["palette"]!;
            var vis = (string)q["vision"]!;
            VisionPairs.Add(new VisionPairRow(pal, char.ToUpperInvariant(vis[0]) + vis[1..], $"{q["a"]} ↔ {q["b"]}",
                                              q["de"]!.GetValue<double>().ToString("0.0", CultureInfo.InvariantCulture)));
        }
        Raise(nameof(HasVisionPairs));
        VisionSummary = pairs.Count == 0 ? "No pair falls below ΔE*ab 12" : $"{pairs.Count} found" + (pairs.Count > 12 ? " · lowest 12 shown" : "");
    }
}
