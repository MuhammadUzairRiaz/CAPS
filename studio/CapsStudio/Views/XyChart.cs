using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>One series of an XyChart. Kind: "line", "dash", "dots" (points on a baseline strip when Strip is set), "bars"
/// (grouped by category, x = category index). Colour: a theme token name ("AccB") or a #RRGGBB.</summary>
public sealed record ChartSeries(string Name, (double X, double Y)[] Points, string Kind = "line", string Colour = "AccB", double Width = 2, bool Legend = true);
/// <summary>A labelled point (a critical point, the current composition).</summary>
public sealed record ChartMarker(double X, double Y, string Colour = "TextB", string Label = "", double Radius = 4);
/// <summary>A dashed horizontal or vertical guide line.</summary>
public sealed record ChartGuide(double At, bool Vertical, string Colour = "DimB", string Label = "");

/// <summary>A small chart for the polymer-statistics pages (design/boards row 18): lines, dashed lines, dots, grouped bars,
/// markers, guide lines and a legend, in the theme's tokens.</summary>
public sealed class XyChart : Control
{
    private ChartSeries[] _series = [];
    private ChartMarker[] _markers = [];
    private ChartGuide[] _guides = [];
    public string XLabel { get; set; } = "";
    public string YLabel { get; set; } = "";
    public double? XMin { get; set; }
    public double? XMax { get; set; }
    public double? YMin { get; set; }
    public double? YMax { get; set; }
    /// <summary>Log-scaled y axis (decades).</summary>
    public bool LogY { get; set; }
    /// <summary>Category names for bar charts (x = 0, 1, …).</summary>
    public string[] Categories { get; set; } = [];
    public bool ShowLegend { get; set; } = true;
    public string Empty { get; set; } = "no data";

    public void Set(ChartSeries[] series, ChartMarker[]? markers = null, ChartGuide[]? guides = null)
    {
        _series = series;
        _markers = markers ?? [];
        _guides = guides ?? [];
        InvalidateVisual();
    }

    private static IBrush B(string c) => c.StartsWith('#') ? new SolidColorBrush(Color.Parse(c)) : Tokens.Brush(c);

    private static double Nice(double span, int ticks)
    {
        var raw = span / Math.Max(1, ticks);
        var mag = Math.Pow(10, Math.Floor(Math.Log10(Math.Max(1e-300, raw))));
        var f = raw / mag;
        return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 2.5 ? 2.5 : f <= 5 ? 5 : 10) * mag;
    }
    private static string Fmt(double v, double step)
    {
        if (v == 0) return "0";
        var dec = 0;   // enough decimals to write the step exactly (0.25 → 2)
        while (dec < 6 && Math.Abs(step * Math.Pow(10, dec) - Math.Round(step * Math.Pow(10, dec))) > 1e-6 * Math.Pow(10, dec)) dec++;
        return v.ToString(dec == 0 ? "0" : "0." + new string('0', dec), CultureInfo.InvariantCulture);
    }

    public override void Render(DrawingContext ctx)
    {
        var b = Bounds;
        var legendH = ShowLegend && _series.Any(s => s.Legend && s.Name.Length > 0) ? 22.0 : 0;
        const double L = 50, R = 14, T = 16;
        var B0 = 34 + legendH;
        var w = b.Width - L - R;
        var h = b.Height - T - B0;
        if (w < 30 || h < 30) return;
        var tf = new Typeface(Tokens.Mono);
        var ui = new Typeface(Tokens.Sans);
        void Text(string s, double x, double y, IBrush brush, double size = 10, bool right = false, bool centre = false, Typeface? face = null)
        {
            var ft = new FormattedText(s, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, face ?? tf, size, brush);
            ctx.DrawText(ft, new Point(right ? x - ft.Width : centre ? x - ft.Width / 2 : x, y - ft.Height / 2));
        }
        var dim = Tokens.Brush("DimB");
        var pts = _series.Where(s => s.Kind != "dots" || s.Points.All(p => double.IsFinite(p.Y))).SelectMany(s => s.Points)
            .Concat(_markers.Select(m => (m.X, m.Y))).Where(p => double.IsFinite(p.X) && double.IsFinite(p.Y) && (!LogY || p.Y > 0)).ToArray();
        var bars = _series.Where(s => s.Kind == "bars").ToArray();
        if (pts.Length == 0 && Categories.Length == 0)
        {
            Text(Empty, L + w / 2, T + h / 2, dim, 11, centre: true, face: ui);
            return;
        }
        double xmin, xmax, ymin, ymax;
        if (bars.Length > 0)
        {
            xmin = -0.5; xmax = Math.Max(Categories.Length, bars.Max(s => s.Points.Length)) - 0.5;
        }
        else
        {
            xmin = XMin ?? pts.Min(p => p.X); xmax = XMax ?? pts.Max(p => p.X);
        }
        if (xmax - xmin < 1e-12) xmax = xmin + 1;
        var ys = pts.Select(p => LogY ? Math.Log10(p.Y) : p.Y).Concat(_guides.Where(g => !g.Vertical).Select(g => LogY ? Math.Log10(Math.Max(1e-300, g.At)) : g.At)).ToArray();
        ymin = YMin ?? (bars.Length > 0 ? 0 : ys.Length > 0 ? ys.Min() : 0);
        ymax = YMax ?? (ys.Length > 0 ? ys.Max() : 1);
        if (!YMax.HasValue) ymax += (ymax - ymin) * 0.08;
        if (ymax - ymin < 1e-12) ymax = ymin + 1;
        double X(double v) => L + (v - xmin) / (xmax - xmin) * w;
        double Y(double v) => T + h - ((LogY ? Math.Log10(Math.Max(1e-300, v)) : v) - ymin) / (ymax - ymin) * h;

        // grid and ticks
        var grid = new Pen(Tokens.Brush("Bg3B"), 1);
        if (LogY)
        {
            for (var d = Math.Ceiling(ymin); d <= ymax + 1e-9; d += Math.Max(1, Math.Round((ymax - ymin) / 4)))
            {
                var yy = T + h - (d - ymin) / (ymax - ymin) * h;
                ctx.DrawLine(grid, new Point(L, yy), new Point(L + w, yy));
                Text("1e" + d.ToString("0", CultureInfo.InvariantCulture), L - 6, yy, dim, right: true);
            }
        }
        else
        {
            var ys0 = Nice(ymax - ymin, 4);
            for (var t = Math.Ceiling(ymin / ys0 - 1e-9) * ys0; t <= ymax + 1e-9 * ys0; t += ys0)
            {
                ctx.DrawLine(grid, new Point(L, Y(t)), new Point(L + w, Y(t)));
                Text(Fmt(t, ys0), L - 6, Y(t), dim, right: true);
            }
        }
        if (bars.Length > 0)
        {
            for (var i = 0; i < Categories.Length; ++i) Text(Categories[i], X(i), T + h + 12, dim, 9.5, centre: true);
        }
        else
        {
            var xs = Nice(xmax - xmin, Math.Max(2, (int)(w / 90)));
            for (var t = Math.Ceiling(xmin / xs - 1e-9) * xs; t <= xmax + 1e-9 * xs; t += xs) Text(Fmt(t, xs), X(t), T + h + 12, dim, centre: true);
        }
        ctx.DrawLine(new Pen(Tokens.Brush("LineB"), 1), new Point(L, T + h), new Point(L + w, T + h));
        if (XLabel.Length > 0) Text(XLabel, L + w, T + h + 26, dim, 10, right: true, face: ui);
        if (YLabel.Length > 0) Text(YLabel, L + 4, T - 6, dim, 10, face: ui);

        using (ctx.PushClip(new Rect(L - 6, T - 6, w + 12, h + 12)))
        {
            foreach (var g in _guides)
            {
                var pen = new Pen(B(g.Colour), 1, new DashStyle([4, 4], 0));
                if (g.Vertical) { ctx.DrawLine(pen, new Point(X(g.At), T), new Point(X(g.At), T + h)); if (g.Label.Length > 0) Text(g.Label, X(g.At) + 4, T + 6, B(g.Colour), 9.5); }
                else { ctx.DrawLine(pen, new Point(L, Y(g.At)), new Point(L + w, Y(g.At))); if (g.Label.Length > 0) Text(g.Label, L + w - 4, Y(g.At) - 8, B(g.Colour), 9.5, right: true); }
            }
            for (var k = 0; k < bars.Length; ++k)
            {
                var s = bars[k];
                var gw = 0.72 / bars.Length;
                foreach (var p in s.Points)
                {
                    if (!double.IsFinite(p.Y) || p.Y <= 0) continue;
                    var x0 = X(p.X - 0.36 + k * gw) + 1;
                    var x1 = X(p.X - 0.36 + (k + 1) * gw) - 1;
                    ctx.FillRectangle(B(s.Colour), new Rect(x0, Y(p.Y), Math.Max(1, x1 - x0), Math.Max(0, Y(ymin) - Y(p.Y))), 2);
                }
            }
            foreach (var s in _series.Where(s => s.Kind is "line" or "dash"))
            {
                var ok = s.Points.Where(p => double.IsFinite(p.X) && double.IsFinite(p.Y) && (!LogY || p.Y > 0)).ToArray();
                if (ok.Length < 2) continue;
                var geo = new StreamGeometry();
                using (var g = geo.Open())
                {
                    g.BeginFigure(new Point(X(ok[0].X), Y(ok[0].Y)), false);
                    foreach (var p in ok.Skip(1)) g.LineTo(new Point(X(p.X), Y(p.Y)));
                    g.EndFigure(false);
                }
                ctx.DrawGeometry(null, new Pen(B(s.Colour), s.Width, s.Kind == "dash" ? new DashStyle([5, 3], 0) : null, lineJoin: PenLineJoin.Round), geo);
            }
            foreach (var s in _series.Where(s => s.Kind == "dots"))
                foreach (var p in s.Points)
                    ctx.DrawEllipse(B(s.Colour), null, new Point(X(p.X), double.IsFinite(p.Y) ? Y(p.Y) : T + h - 8), s.Width, s.Width);
            foreach (var m in _markers)
            {
                ctx.DrawEllipse(B(m.Colour), new Pen(Tokens.Brush("Bg1B"), 1.5), new Point(X(m.X), Y(m.Y)), m.Radius, m.Radius);
                if (m.Label.Length > 0) Text(m.Label, X(m.X) + m.Radius + 4, Y(m.Y) - 8, B(m.Colour), 9.5);
            }
        }
        if (legendH > 0)
        {
            var x = L;
            var y = b.Height - 10;
            foreach (var s in _series.Where(s => s.Legend && s.Name.Length > 0))
            {
                if (s.Kind == "dash") ctx.DrawLine(new Pen(B(s.Colour), 2, new DashStyle([3, 2], 0)), new Point(x, y), new Point(x + 12, y));
                else ctx.DrawEllipse(B(s.Colour), null, new Point(x + 4, y), 4, 4);
                var ft = new FormattedText(s.Name, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, ui, 11, Tokens.Brush("MutedB"));
                ctx.DrawText(ft, new Point(x + (s.Kind == "dash" ? 17 : 12), y - ft.Height / 2));
                x += ft.Width + 34;
            }
        }
    }
}
