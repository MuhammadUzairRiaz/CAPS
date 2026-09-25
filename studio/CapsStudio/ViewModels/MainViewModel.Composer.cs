using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json.Nodes;
using Avalonia.Media;
using Avalonia.Media.Imaging;
using CapsStudio.Interop;
using CapsStudio.Views;

namespace CapsStudio.ViewModels;

public sealed record ComposerCheck(string Text, bool Ok);

/// <summary>Export › Figure composer (design/boards/FigureComposer): panels a, b, c … laid out on a journal page — the 3D
/// view or any curve already in the project — with line, colour and text settings checked at print size, exported as
/// SVG (vector plots), PNG and TIFF at the chosen dpi, or PDF (the page as an image at that dpi).</summary>
public partial class MainViewModel
{
    public bool IsComposer => _module == 54;
    public static readonly string[] ComposerPresets = ["Single column · 3.5 in", "Double column · 7.0 in", "Full page · 7 × 9 in", "Slide · 13.33 × 7.5 in"];
    private static readonly (double W, double H)[] PresetSizes = [(3.5, 2.6), (7.0, 4.2), (7.0, 9.0), (13.33, 7.5)];
    public static readonly string[] ComposerLayouts = ["1 × 1", "1 × 2", "2 × 1", "2 × 2", "1 × 3", "3 × 1"];
    private static readonly (int R, int C)[] LayoutShapes = [(1, 1), (1, 2), (2, 1), (2, 2), (1, 3), (3, 1)];
    public static readonly string[] ComposerFonts = ["Arial", "IBM Plex Sans", "Times New Roman"];
    public static readonly string[] ComposerColourNames = ["amber · v2", "blue · v2", "red · v2", "violet · v2", "green · v2", "ink"];
    private static readonly string[] ComposerColourHex = ["#E8893A", "#3B7DD8", "#DE775D", "#9B7AD5", "#5FAF6E", "#1E2226"];

    public FigureComposition Composition { get; } = new();
    public ObservableCollection<string> ComposerSources { get; } = new();
    public ObservableCollection<ComposerCheck> ComposerChecks { get; } = new();
    private int _cmpPreset = 1, _cmpLayout = 3, _cmpPanel, _cmpFont, _cmpDpi = 600;
    private readonly int[] _cmpSource = new int[9], _cmpColour = new int[9];
    private readonly double[] _cmpLine = Enumerable.Repeat(1.0, 9).ToArray();
    public event Action? ComposerChanged;

    public int ComposerPreset { get => _cmpPreset; set { if (Set(ref _cmpPreset, Math.Clamp(value, 0, 3))) { (Composition.WidthIn, Composition.HeightIn) = PresetSizes[_cmpPreset]; RaiseComposer(); } } }
    public decimal ComposerWidth { get => (decimal)Composition.WidthIn; set { Composition.WidthIn = (double)Math.Clamp(value, 1m, 40m); RaiseComposer(); } }
    public decimal ComposerHeight { get => (decimal)Composition.HeightIn; set { Composition.HeightIn = (double)Math.Clamp(value, 1m, 40m); RaiseComposer(); } }
    public decimal ComposerDpi { get => _cmpDpi; set { _cmpDpi = (int)Math.Clamp(value, 72, 1200); RaiseComposer(); } }
    public string ComposerOutput => $"{(int)Math.Round(Composition.WidthIn * _cmpDpi)} × {(int)Math.Round(Composition.HeightIn * _cmpDpi)} px";
    public int ComposerLayout { get => _cmpLayout; set { if (Set(ref _cmpLayout, Math.Clamp(value, 0, 5))) { (Composition.Rows, Composition.Cols) = LayoutShapes[_cmpLayout]; _cmpPanel = Math.Min(_cmpPanel, PanelCount - 1); RaiseComposer(); } } }
    private int PanelCount => Composition.Rows * Composition.Cols;
    public int ComposerPanel { get => _cmpPanel; set { if (Set(ref _cmpPanel, Math.Clamp(value, 0, PanelCount - 1))) RaisePanel(); } }
    public string ComposerPanelTitle => $"PANEL {(char)('A' + _cmpPanel)}";
    public int ComposerSource { get => _cmpSource[_cmpPanel]; set { _cmpSource[_cmpPanel] = Math.Max(0, value); RaiseComposer(); } }
    public int ComposerColour { get => _cmpColour[_cmpPanel]; set { _cmpColour[_cmpPanel] = Math.Clamp(value, 0, ComposerColourHex.Length - 1); RaiseComposer(); } }
    public decimal ComposerLine { get => (decimal)_cmpLine[_cmpPanel]; set { _cmpLine[_cmpPanel] = (double)Math.Clamp(value, 0.1m, 6m); RaiseComposer(); } }
    public int ComposerFont { get => _cmpFont; set { if (Set(ref _cmpFont, Math.Clamp(value, 0, 2))) { Composition.Font = ComposerFonts[_cmpFont]; RaiseComposer(); } } }
    public decimal ComposerFontPt { get => (decimal)Composition.FontPt; set { Composition.FontPt = (double)Math.Clamp(value, 4m, 24m); RaiseComposer(); } }
    public bool ComposerLabels { get => Composition.Labels; set { Composition.Labels = value; RaiseComposer(); } }
    public string ComposerSummary => $"{PanelCount} panel{(PanelCount == 1 ? "" : "s")} · {ComposerOutput} at {_cmpDpi} dpi";

    public void OpenComposer()
    {
        RefreshComposerSources();
        // a first layout from what the project has: the view, then the first curves
        if (_cmpSource.All(s => s == 0))
        {
            _cmpSource[0] = _doc != null ? 1 : 0;
            // the curves with data, first come first
            var withData = Enumerable.Range(2, Math.Max(0, ComposerSources.Count - 2))
                .Where(i => ComposerSources[i].StartsWith("Dynamics log") || ComposerSources[i].StartsWith("Bars") ||
                            Analyze.Curves.FirstOrDefault(c => $"{c.Property} · {c.Label}" == ComposerSources[i]) is { } c && c.X.Length > 1)
                .ToList();
            for (var k = 1; k < 4 && k - 1 < withData.Count; ++k) _cmpSource[k] = withData[k - 1];
            for (var k = 0; k < 9; ++k) _cmpColour[k] = k % 2 == 0 ? 1 : 0;
        }
        (Composition.WidthIn, Composition.HeightIn) = PresetSizes[_cmpPreset];
        (Composition.Rows, Composition.Cols) = LayoutShapes[_cmpLayout];
        SetModule(54);
        RaiseComposer();
    }

    /// <summary>What a panel can show: nothing, the 3D view, every curve of the last analysis, the dynamics log.</summary>
    public void RefreshComposerSources()
    {
        ComposerSources.Clear();
        ComposerSources.Add("Empty");
        ComposerSources.Add("3D view · " + (Title.Length > 0 ? Title : "no structure"));
        foreach (var c in Analyze.Curves) ComposerSources.Add($"{c.Property} · {c.Label}");
        if (_thermo.Count > 1)
            foreach (var t in new[] { "temperature", "density", "potential energy" }) ComposerSources.Add("Dynamics log · " + t);
        var unit = Analyze.Results.Where(r => double.IsFinite(r.Value)).GroupBy(r => r.Unit).OrderByDescending(g => g.Count()).FirstOrDefault();
        if (unit != null && unit.Count() > 1) ComposerSources.Add($"Bars · results in {(unit.Key.Length > 0 ? unit.Key : "one unit")}");
    }

    private void RaisePanel()
    {
        foreach (var n in new[] { nameof(ComposerPanelTitle), nameof(ComposerSource), nameof(ComposerColour), nameof(ComposerLine) }) Raise(n);
        ComposerChanged?.Invoke();
    }

    public void RaiseComposer()
    {
        BuildPanels(null);
        foreach (var n in new[] { nameof(ComposerWidth), nameof(ComposerHeight), nameof(ComposerDpi), nameof(ComposerOutput), nameof(ComposerSummary), nameof(ComposerPanelTitle),
                                  nameof(ComposerSource), nameof(ComposerColour), nameof(ComposerLine), nameof(ComposerFontPt), nameof(ComposerLabels) }) Raise(n);
        CheckComposer();
        ComposerChanged?.Invoke();
    }

    /// <summary>Fills the composition's panels from the sources; the 3D view is rendered at pixelsPerInch (null: preview).</summary>
    public void BuildPanels(double? pixelsPerInch)
    {
        var inv = CultureInfo.InvariantCulture;
        Composition.Panels.Clear();
        for (var k = 0; k < PanelCount; ++k)
        {
            var src = _cmpSource[k] < ComposerSources.Count ? ComposerSources[_cmpSource[k]] : "Empty";
            var p = new FigurePanel { Source = src, LinePt = _cmpLine[k], Colour = Color.Parse(ComposerColourHex[_cmpColour[k]]) };
            if (src.StartsWith("3D view") && _doc != null)
            {
                p.Kind = "view";
                var r = Composition.PanelRect(k, 1);
                var ppi = pixelsPerInch ?? 110;
                var (w, h) = ((int)Math.Clamp(r.Width / 72 * ppi, 32, 8000), (int)Math.Clamp((r.Height - Composition.FontPt - 3) / 72 * ppi, 32, 8000));
                try
                {
                    var o = new CapsRenderOpts { Width = w, Height = h, Supersample = 2, Background = 1, Style = _style, ColourBy = _colour, Outlines = 1, DepthCue = 1, ShowCell = 1,
                                                 Highlight0 = -1, Highlight1 = -1, Highlight2 = -1, Highlight3 = -1 };
                    var rgba = new byte[w * h * 4];
                    _doc.Render(Camera, o, rgba);
                    p.View = ToBitmap(rgba, w, h);
                }
                catch { p.View = null; }
            }
            else if (src.StartsWith("Dynamics log · "))
            {
                p.Kind = "curve";
                var what = src["Dynamics log · ".Length..];
                p.Points = _thermo.Select(t => (t.TimePs, what == "temperature" ? t.Temperature : what == "density" ? t.Density : t.Potential)).ToArray();
                (p.XLabel, p.YLabel) = ("time (ps)", what == "temperature" ? "T (K)" : what == "density" ? "ρ (g/cm³)" : "E (kcal/mol)");
            }
            else if (src.StartsWith("Bars · "))
            {
                p.Kind = "bars";
                var unit = Analyze.Results.Where(r => double.IsFinite(r.Value)).GroupBy(r => r.Unit).OrderByDescending(g => g.Count()).First();
                p.Bars = unit.Select(r => (r.Name.Length > 12 ? r.Name[..12] : r.Name, r.Value)).ToArray();
                p.YLabel = unit.Key;
            }
            else if (Analyze.Curves.FirstOrDefault(c => $"{c.Property} · {c.Label}" == src) is { } c)
            {
                p.Kind = "curve";
                p.Points = c.X.Zip(c.Y).Where(q => double.IsFinite(q.First) && double.IsFinite(q.Second)).ToArray();
                (p.XLabel, p.YLabel, p.RefY) = (c.XLabel, c.YLabel, c.RefY);
            }
            Composition.Panels.Add(p);
        }
    }

    private void CheckComposer()
    {
        ComposerChecks.Clear();
        var thin = Composition.Panels.Where(p => p.Kind is "curve" && p.LinePt < 0.5).Select((p, i) => (char)('a' + Composition.Panels.IndexOf(p))).ToList();
        ComposerChecks.Add(new ComposerCheck(thin.Count == 0 ? "Lines ≥ 0.5 pt at print size" : $"Lines under 0.5 pt in panel {string.Join(", ", thin)}", thin.Count == 0));
        ComposerChecks.Add(new ComposerCheck(Composition.FontPt >= 5 ? $"Text {Composition.FontPt.ToString("0.#", CultureInfo.InvariantCulture)} pt (≥ 5 pt, axis labels 1 pt smaller)" : "Text under 5 pt: too small in print", Composition.FontPt >= 5));
        // the panels' colours, checked for colour-vision deficiencies
        var cols = Composition.Panels.Where(p => p.HasData && p.Kind != "view").Select(p => p.Colour).Distinct().ToList();
        var cvdOk = true;
        if (cols.Count > 1)
            try
            {
                var input = new JsonObject { ["Panels"] = new JsonObject { ["colours"] = new JsonArray(cols.Select(c => (JsonNode)$"#{c.R:X2}{c.G:X2}{c.B:X2}").ToArray()) } };
                cvdOk = ((JsonArray)JsonNode.Parse(CapsDocument.VisionCheck(input.ToJsonString(), 12))!["pairs"]!).Count == 0;
            }
            catch { }
        ComposerChecks.Add(new ComposerCheck(cvdOk ? "Palette v2 · colour-vision safe" : "Two panel colours are hard to tell apart for colour-blind readers", cvdOk));
        var empty = Composition.Panels.Select((p, i) => (p, i)).Where(x => !x.p.HasData).Select(x => (char)('a' + x.i)).ToList();
        ComposerChecks.Add(new ComposerCheck(empty.Count == 0 ? "Every panel has data" : $"Panel {string.Join(", ", empty)} has no data yet", empty.Count == 0));
    }

    public string ExportComposerSvg(string path)
    {
        BuildPanels(_cmpDpi);
        var svg = Composition.ToSvg(p =>
        {
            if (p.View == null) return null;
            using var ms = new MemoryStream();
            p.View.Save(ms);
            return Convert.ToBase64String(ms.ToArray());
        });
        File.WriteAllText(path, svg);
        BuildPanels(null);
        return $"Wrote {Path.GetFileName(path)} · plots as vectors, the 3D panel as an image";
    }

    /// <summary>PNG, TIFF or PDF of the page at the chosen dpi.</summary>
    public string ExportComposerRaster(string path)
    {
        BuildPanels(_cmpDpi);
        var (w, h) = ((int)Math.Round(Composition.WidthIn * _cmpDpi), (int)Math.Round(Composition.HeightIn * _cmpDpi));
        using var rtb = new RenderTargetBitmap(new Avalonia.PixelSize(w, h), new Avalonia.Vector(96, 96));
        using (var ctx = rtb.CreateDrawingContext()) Composition.Draw(ctx, _cmpDpi / 72.0);
        var ext = Path.GetExtension(path).ToLowerInvariant();
        string note;
        if (ext == ".png") { rtb.Save(path); FigureFiles.SetPngDpi(path, _cmpDpi); note = "PNG"; }
        else
        {
            var rgba = new byte[w * h * 4];
            unsafe { fixed (byte* p = rgba) rtb.CopyPixels(new Avalonia.PixelRect(0, 0, w, h), (IntPtr)p, rgba.Length, w * 4); }
            for (var i = 0; i < rgba.Length; i += 4) (rgba[i], rgba[i + 2]) = (rgba[i + 2], rgba[i]);   // BGRA → RGBA
            if (ext == ".pdf") { FigureFiles.WritePdf(path, rgba, w, h, Composition.WidthPt, Composition.HeightPt); note = "PDF (the page as an image)"; }
            else { FigureFiles.WriteTiff(path, rgba, w, h, _cmpDpi); note = "TIFF"; }
        }
        BuildPanels(null);
        return $"Wrote {Path.GetFileName(path)} · {note} · {w} × {h} px at {_cmpDpi} dpi";
    }
}
