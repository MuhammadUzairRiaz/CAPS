using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Avalonia;
using Avalonia.Media;
using Avalonia.Media.Imaging;

namespace CapsStudio.Views;

/// <summary>One panel of a composed figure: the 3D view (a bitmap and its pixels) or a plot (a SciFigure: a curve of the
/// project, the dynamics log, or bars over results).</summary>
public sealed class FigurePanel
{
    public string Source = "";                 // shown in the inspector
    public string Kind = "empty";              // empty | view | curve | bars
    public SciFigure? Figure;
    public Bitmap? View;                        // the 3D panel, rendered at the export size
    public byte[]? ViewRgba;                    // its pixels (RGBA, top row first) for the PDF
    public int ViewW, ViewH;
    public double LinePt = 1.0;
    public Color Colour = Color.Parse(SciStyle.PaperPalette[0]);
    public bool HasData => Kind == "view" ? View != null : Figure != null && Figure.Series.Any(s => s.X.Length > 0);
}

/// <summary>A multi-panel figure (design/boards/FigureComposer) drawn at any scale: on screen, into a bitmap for PNG
/// and TIFF, as vector SVG or as a vector PDF. Sizes are in points (1/72 in); panel labels a, b, c … bold. The plot
/// panels are drawn the way every figure CAPS writes for print is drawn (SciPlot in the Paper style: a closed frame,
/// inward ticks on four sides, round values, quantity / unit titles), in the panel's colour and line weight.</summary>
public sealed class FigureComposition
{
    public double WidthIn = 7.0, HeightIn = 4.2;
    public int Rows = 2, Cols = 2;
    public string Font = "Arial";
    public double FontPt = 7;
    public bool Labels = true;
    public List<FigurePanel> Panels = new();

    public double WidthPt => WidthIn * 72;
    public double HeightPt => HeightIn * 72;

    public Rect PanelRect(int k, double scale)
    {
        const double pad = 6;   // pt between panels and at the page edge
        var w = (WidthPt - pad * (Cols + 1)) / Cols;
        var h = (HeightPt - pad * (Rows + 1)) / Rows;
        var r = k / Cols;
        var c = k % Cols;
        return new Rect((pad + c * (w + pad)) * scale, (pad + r * (h + pad)) * scale, w * scale, h * scale);
    }

    /// <summary>The panel's area below its letter, in points.</summary>
    private Rect Inner(int k) => PanelRect(k, 1).Deflate(new Thickness(0, Labels ? FontPt + 3 : 0, 0, 0));

    /// <summary>The print style for a panel: its colour first in the palette, the composer's text size.</summary>
    private SciStyle StyleOf(FigurePanel p)
    {
        var hex = $"#{p.Colour.R:X2}{p.Colour.G:X2}{p.Colour.B:X2}";
        var paper = SciStyle.Paper;
        return new SciStyle
        {
            Print = true, Light = true, Ink = paper.Ink, Soft = paper.Soft, Bg = paper.Bg, Rule = paper.Rule, Font = new FontFamily(Font),
            Palette = [hex, .. SciStyle.PaperPalette.Where(c => !c.Equals(hex, StringComparison.OrdinalIgnoreCase))], BasePt = FontPt + 1,
        };
    }

    /// <summary>The panel's figure with its line weight (1 = the print default of 1 pt for a curve).</summary>
    private static SciFigure Weighted(FigurePanel p) => p.Figure! with { Series = p.Figure!.Series.Select(s => s with { Width = s.Width * p.LinePt }).ToArray() };

    /// <summary>Everything but the 3D panels drawn into a sink (points); view(k, rect) draws a 3D panel.</summary>
    private void Compose(ISciSink g, Action<FigurePanel, Rect> view)
    {
        var ink = Colors.Black;
        var dim = Color.Parse("#6B7178");
        for (var k = 0; k < Panels.Count && k < Rows * Cols; ++k)
        {
            var p = Panels[k];
            var r = PanelRect(k, 1);
            if (Labels) g.Text(((char)('a' + k)).ToString(), new Point(r.X, r.Y), FontPt + 1, false, true, ink);
            var inner = Inner(k);
            if (p.Kind == "view" && p.View != null)
            {
                var s = Math.Min(inner.Width / p.View.Size.Width, inner.Height / p.View.Size.Height);
                var (w, h) = (p.View.Size.Width * s, p.View.Size.Height * s);
                view(p, new Rect(inner.X + (inner.Width - w) / 2, inner.Y + (inner.Height - h) / 2, w, h));
                continue;
            }
            if (p.HasData)
            {
                g.PushTransform(Matrix.CreateTranslation(inner.X, inner.Y));
                SciPlot.DrawTo(g, Weighted(p), inner.Size, StyleOf(p));
                g.Pop();
                continue;
            }
            var box = inner.Deflate(2);
            g.Rect(box, null, 0, Color.Parse("#D9D5CF"), 0.75, [4, 3]);
            var text = p.Kind == "empty" ? "empty panel" : "no data yet";
            g.Text(text, new Point(box.Center.X - g.Measure(text, FontPt, false, false) / 2, box.Center.Y - g.Metrics(FontPt).Height / 2), FontPt, false, false, dim);
        }
    }

    /// <summary>Draws the page; scale = device pixels per point.</summary>
    public void Draw(DrawingContext ctx, double scale)
    {
        ctx.FillRectangle(Brushes.White, new Rect(0, 0, WidthPt * scale, HeightPt * scale));
        using var t = ctx.PushTransform(Matrix.CreateScale(scale, scale));
        Compose(new ScreenSink(ctx, new FontFamily(Font)), (p, r) => ctx.DrawImage(p.View!, new Rect(p.View!.Size), r));
    }

    /// <summary>SVG with the plots as vectors (the 3D panel as an embedded PNG); sizes in points.</summary>
    public string ToSvg(Func<FigurePanel, string?> viewPngBase64)
    {
        var s = new SvgSink(WidthPt, HeightPt);
        Compose(s, (p, r) =>
        {
            if (viewPngBase64(p) is { } b64)
                s.Raw($"<image x=\"{SciExport.F(r.X)}\" y=\"{SciExport.F(r.Y)}\" width=\"{SciExport.F(r.Width)}\" height=\"{SciExport.F(r.Height)}\" href=\"data:image/png;base64,{b64}\"/>\n");
        });
        return s.Document();
    }

    /// <summary>A one-page vector PDF: the plots as paths and text (fonts embedded), the 3D panel as an image.</summary>
    public byte[] ToPdf()
    {
        var pdf = new PdfSink(WidthPt, HeightPt);
        Compose(pdf, (p, r) => { if (p.ViewRgba != null) pdf.Image(p.ViewRgba, p.ViewW, p.ViewH, r); });
        return pdf.Document();
    }
}
