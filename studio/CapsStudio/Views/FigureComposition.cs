using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Avalonia;
using Avalonia.Media;
using Avalonia.Media.Imaging;

namespace CapsStudio.Views;

/// <summary>One panel of a composed figure: the 3D view (a bitmap), a curve (x, y), or bars (label, value).</summary>
public sealed class FigurePanel
{
    public string Source = "";                 // shown in the inspector
    public string Kind = "empty";              // empty | view | curve | bars
    public (double X, double Y)[] Points = [];
    public (string Label, double Value)[] Bars = [];
    public string XLabel = "", YLabel = "";
    public double? RefY;
    public Bitmap? View;                        // the 3D panel, rendered at the export size
    public double LinePt = 1.0;
    public Color Colour = Color.Parse("#E8893A");
    public bool HasData => Kind == "view" ? View != null : Kind == "curve" ? Points.Length > 1 : Kind == "bars" && Bars.Length > 0;
}

/// <summary>A multi-panel figure (design/boards/FigureComposer) drawn at any scale — on screen, into a bitmap for PNG,
/// TIFF and PDF, or as SVG with the plots as vectors. Sizes are in points (1/72 in); panel labels a, b, c … bold.</summary>
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

    private static string Fmt(double v) => Math.Abs(v) >= 1000 ? v.ToString("0", CultureInfo.InvariantCulture)
        : Math.Abs(v) >= 10 ? v.ToString("0.#", CultureInfo.InvariantCulture) : v.ToString("0.##", CultureInfo.InvariantCulture);

    private static (double Lo, double Hi) RangeX(IEnumerable<double> v)
    {
        var a = v.Where(double.IsFinite).ToArray();
        if (a.Length == 0) return (0, 1);
        var (lo, hi) = (a.Min(), a.Max());
        return hi - lo < 1e-12 ? (lo, lo + 1) : (lo, hi);
    }

    private static (double Lo, double Hi) Range(IEnumerable<double> v, bool fromZero)
    {
        var a = v.Where(double.IsFinite).ToArray();
        if (a.Length == 0) return (0, 1);
        double lo = fromZero ? Math.Min(0, a.Min()) : a.Min(), hi = a.Max();
        if (hi - lo < 1e-12) hi = lo + 1;
        var pad = (hi - lo) * 0.06;
        return (fromZero && lo >= 0 ? 0 : lo - pad, hi + pad);
    }

    /// <summary>Draws the page; scale = device pixels per point.</summary>
    public void Draw(DrawingContext ctx, double scale)
    {
        ctx.FillRectangle(Brushes.White, new Rect(0, 0, WidthPt * scale, HeightPt * scale));
        var tf = new Typeface(Font);
        var tfb = new Typeface(Font, FontStyle.Normal, FontWeight.Bold);
        var ink = new SolidColorBrush(Color.Parse("#1E2226"));
        var dim = new SolidColorBrush(Color.Parse("#6B7178"));
        var grid = new Pen(new SolidColorBrush(Color.Parse("#E7E4DF")), Math.Max(0.5, 0.5 * scale));
        var axis = new Pen(new SolidColorBrush(Color.Parse("#555B62")), Math.Max(0.5, 0.6 * scale));
        FormattedText Text(string s, bool bold = false, double pt = 0, IBrush? b = null) =>
            new(s, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, bold ? tfb : tf, (pt > 0 ? pt : FontPt) * scale, b ?? ink);
        for (var k = 0; k < Panels.Count && k < Rows * Cols; ++k)
        {
            var p = Panels[k];
            var r = PanelRect(k, scale);
            if (Labels) ctx.DrawText(Text(((char)('a' + k)).ToString(), true, FontPt + 1), new Point(r.X, r.Y));
            var inner = r.Deflate(new Thickness(0, (FontPt + 3) * scale, 0, 0));
            switch (p.Kind)
            {
                case "view" when p.View != null:
                {
                    var src = new Rect(p.View.Size);
                    var s = Math.Min(inner.Width / src.Width, inner.Height / src.Height);
                    var dst = new Rect(inner.X + (inner.Width - src.Width * s) / 2, inner.Y + (inner.Height - src.Height * s) / 2, src.Width * s, src.Height * s);
                    ctx.DrawImage(p.View, src, dst);
                    break;
                }
                case "curve" when p.Points.Length > 1:
                case "bars" when p.Bars.Length > 0:
                {
                    var L = inner.X + 26 * scale; var T = inner.Y + 4 * scale;
                    var W = inner.Right - L - 4 * scale; var H = inner.Bottom - T - 16 * scale;
                    if (W < 10 || H < 10) break;
                    var bars = p.Kind == "bars";
                    var (xlo, xhi) = bars ? (0.0, p.Bars.Length) : RangeX(p.Points.Select(q => q.X));
                    var (ylo, yhi) = bars ? Range(p.Bars.Select(q => q.Value), true) : Range(p.Points.Select(q => q.Y).Concat(p.RefY is double ry ? [ry] : []), false);
                    double X(double v) => L + (v - xlo) / (xhi - xlo) * W;
                    double Y(double v) => T + H - (v - ylo) / (yhi - ylo) * H;
                    foreach (var t in new[] { ylo, (ylo + yhi) / 2, yhi })
                    {
                        ctx.DrawLine(grid, new Point(L, Y(t)), new Point(L + W, Y(t)));
                        var ft = Text(Fmt(t), pt: FontPt - 1, b: dim);
                        ctx.DrawText(ft, new Point(L - ft.Width - 3 * scale, Y(t) - ft.Height / 2));
                    }
                    ctx.DrawLine(axis, new Point(L, T + H), new Point(L + W, T + H));
                    ctx.DrawLine(axis, new Point(L, T), new Point(L, T + H));
                    var pen = new Pen(new SolidColorBrush(p.Colour), p.LinePt * scale, lineJoin: PenLineJoin.Round);
                    if (bars)
                    {
                        var fill = new SolidColorBrush(p.Colour, 0.8);
                        for (var i = 0; i < p.Bars.Length; ++i)
                        {
                            var x0 = X(i + 0.18); var x1 = X(i + 0.82);
                            ctx.FillRectangle(fill, new Rect(x0, Y(p.Bars[i].Value), x1 - x0, Y(ylo) - Y(p.Bars[i].Value)));
                            var lt = Text(p.Bars[i].Label, pt: FontPt - 1, b: dim);
                            ctx.DrawText(lt, new Point(X(i + 0.5) - lt.Width / 2, T + H + 3 * scale));
                        }
                    }
                    else
                    {
                        foreach (var t in new[] { xlo, (xlo + xhi) / 2, xhi })
                        {
                            var ft = Text(Fmt(t), pt: FontPt - 1, b: dim);
                            ctx.DrawText(ft, new Point(Math.Clamp(X(t) - ft.Width / 2, L, L + W - ft.Width), T + H + 3 * scale));
                        }
                        if (p.RefY is double rv && rv > ylo && rv < yhi)
                            ctx.DrawLine(new Pen(dim, 0.5 * scale, new DashStyle([4, 3], 0)), new Point(L, Y(rv)), new Point(L + W, Y(rv)));
                        var geo = new StreamGeometry();
                        using (var g = geo.Open())
                        {
                            g.BeginFigure(new Point(X(p.Points[0].X), Y(p.Points[0].Y)), false);
                            foreach (var q in p.Points.Skip(1)) g.LineTo(new Point(X(q.X), Y(q.Y)));
                            g.EndFigure(false);
                        }
                        ctx.DrawGeometry(null, pen, geo);
                    }
                    var yl = Text(p.YLabel, pt: FontPt - 1, b: dim);
                    ctx.DrawText(yl, new Point(L + 2 * scale, T - yl.Height + 2 * scale));
                    if (!bars)
                    {
                        var xl = Text(p.XLabel, pt: FontPt - 1, b: dim);
                        ctx.DrawText(xl, new Point(L + W - xl.Width, T + H - xl.Height - 1 * scale));
                    }
                    break;
                }
                default:
                {
                    var box = new Pen(new SolidColorBrush(Color.Parse("#D9D5CF")), 0.75 * scale, new DashStyle([4, 3], 0));
                    ctx.DrawRectangle(null, box, inner.Deflate(2 * scale), 4 * scale, 4 * scale);
                    var ft = Text(p.Kind == "empty" ? "empty panel" : "no data yet", pt: FontPt, b: dim);
                    ctx.DrawText(ft, new Point(inner.Center.X - ft.Width / 2, inner.Center.Y - ft.Height / 2));
                    break;
                }
            }
        }
    }

    /// <summary>SVG with the plots as vectors (the 3D panel as an embedded PNG); sizes in points.</summary>
    public string ToSvg(Func<FigurePanel, string?> viewPngBase64)
    {
        var inv = CultureInfo.InvariantCulture;
        string F(double v) => v.ToString("0.##", inv);
        string Esc(string s) => s.Replace("&", "&amp;").Replace("<", "&lt;").Replace(">", "&gt;");
        string Hex(Color c) => $"#{c.R:X2}{c.G:X2}{c.B:X2}";
        var sb = new StringBuilder();
        sb.Append($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{F(WidthIn)}in\" height=\"{F(HeightIn)}in\" viewBox=\"0 0 {F(WidthPt)} {F(HeightPt)}\" font-family=\"{Esc(Font)}, sans-serif\" font-size=\"{F(FontPt)}\">\n");
        sb.Append($"<rect width=\"{F(WidthPt)}\" height=\"{F(HeightPt)}\" fill=\"#FFFFFF\"/>\n");
        for (var k = 0; k < Panels.Count && k < Rows * Cols; ++k)
        {
            var p = Panels[k];
            var r = PanelRect(k, 1);
            if (Labels) sb.Append($"<text x=\"{F(r.X)}\" y=\"{F(r.Y + FontPt + 1)}\" font-weight=\"bold\" font-size=\"{F(FontPt + 1)}\">{(char)('a' + k)}</text>\n");
            var inner = r.Deflate(new Thickness(0, FontPt + 3, 0, 0));
            if (p.Kind == "view" && viewPngBase64(p) is { } b64 && p.View != null)
            {
                var s = Math.Min(inner.Width / p.View.Size.Width, inner.Height / p.View.Size.Height);
                var (w, h) = (p.View.Size.Width * s, p.View.Size.Height * s);
                sb.Append($"<image x=\"{F(inner.X + (inner.Width - w) / 2)}\" y=\"{F(inner.Y + (inner.Height - h) / 2)}\" width=\"{F(w)}\" height=\"{F(h)}\" href=\"data:image/png;base64,{b64}\"/>\n");
                continue;
            }
            if (!p.HasData) { sb.Append($"<text x=\"{F(inner.Center.X)}\" y=\"{F(inner.Center.Y)}\" text-anchor=\"middle\" fill=\"#6B7178\">no data</text>\n"); continue; }
            double L = inner.X + 26, T = inner.Y + 4, W = inner.Right - L - 4, H = inner.Bottom - T - 16;
            var bars = p.Kind == "bars";
            var (xlo, xhi) = bars ? (0.0, p.Bars.Length) : RangeX(p.Points.Select(q => q.X));
            var (ylo, yhi) = bars ? Range(p.Bars.Select(q => q.Value), true) : Range(p.Points.Select(q => q.Y), false);
            double X(double v) => L + (v - xlo) / (xhi - xlo) * W;
            double Y(double v) => T + H - (v - ylo) / (yhi - ylo) * H;
            foreach (var t in new[] { ylo, (ylo + yhi) / 2, yhi })
                sb.Append($"<line x1=\"{F(L)}\" y1=\"{F(Y(t))}\" x2=\"{F(L + W)}\" y2=\"{F(Y(t))}\" stroke=\"#E7E4DF\" stroke-width=\"0.5\"/><text x=\"{F(L - 3)}\" y=\"{F(Y(t) + 2)}\" text-anchor=\"end\" font-size=\"{F(FontPt - 1)}\" fill=\"#6B7178\">{Fmt(t)}</text>\n");
            sb.Append($"<path d=\"M{F(L)},{F(T)} V{F(T + H)} H{F(L + W)}\" fill=\"none\" stroke=\"#555B62\" stroke-width=\"0.6\"/>\n");
            if (bars)
                for (var i = 0; i < p.Bars.Length; ++i)
                    sb.Append($"<rect x=\"{F(X(i + 0.18))}\" y=\"{F(Y(p.Bars[i].Value))}\" width=\"{F(X(i + 0.82) - X(i + 0.18))}\" height=\"{F(Y(ylo) - Y(p.Bars[i].Value))}\" fill=\"{Hex(p.Colour)}\" fill-opacity=\"0.8\"/>" +
                              $"<text x=\"{F(X(i + 0.5))}\" y=\"{F(T + H + FontPt)}\" text-anchor=\"middle\" font-size=\"{F(FontPt - 1)}\" fill=\"#6B7178\">{Esc(p.Bars[i].Label)}</text>\n");
            else
            {
                foreach (var t in new[] { xlo, (xlo + xhi) / 2, xhi })
                    sb.Append($"<text x=\"{F(X(t))}\" y=\"{F(T + H + FontPt)}\" text-anchor=\"middle\" font-size=\"{F(FontPt - 1)}\" fill=\"#6B7178\">{Fmt(t)}</text>\n");
                sb.Append($"<polyline fill=\"none\" stroke=\"{Hex(p.Colour)}\" stroke-width=\"{F(p.LinePt)}\" stroke-linejoin=\"round\" points=\"");
                sb.Append(string.Join(" ", p.Points.Select(q => $"{F(X(q.X))},{F(Y(q.Y))}")));
                sb.Append("\"/>\n");
                sb.Append($"<text x=\"{F(L + W)}\" y=\"{F(T + H - 2)}\" text-anchor=\"end\" font-size=\"{F(FontPt - 1)}\" fill=\"#6B7178\">{Esc(p.XLabel)}</text>\n");
            }
            sb.Append($"<text x=\"{F(L + 2)}\" y=\"{F(T - 1)}\" font-size=\"{F(FontPt - 1)}\" fill=\"#6B7178\">{Esc(p.YLabel)}</text>\n");
        }
        sb.Append("</svg>\n");
        return sb.ToString();
    }
}
