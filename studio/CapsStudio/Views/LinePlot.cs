using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>Minimal line plot for analysis curves: axes, ticks, a dashed reference line, the curve. By default the axes
/// start at zero with a reference at y = 1 (g(r)); with AutoRange both axes follow the data.</summary>
public sealed class LinePlot : Control
{
    private (double X, double Y)[] _data = [];
    public string XLabel { get; set; } = "";
    public string YLabel { get; set; } = "";
    public bool AutoRange { get; set; }
    public double? RefY { get; set; } = 1.0;

    public void SetData((double X, double Y)[] data)
    {
        _data = data;
        _overlay = [];
        InvalidateVisual();
    }

    /// <summary>Data as points with a line through them (a fit or a smoothed curve).</summary>
    public void SetData((double X, double Y)[] points, (double X, double Y)[] line)
    {
        _data = points;
        _overlay = line;
        Markers = true;
        InvalidateVisual();
    }
    private (double X, double Y)[] _overlay = [];
    public bool Markers { get; set; }

    // design tokens of the current theme (plot in lib.py: grid Bg3, axis line, labels dim, series amber, points cyan)
    private static IBrush Grid => Tokens.Brush("Bg3B");
    private static IBrush Axis => Tokens.Brush("LineB");
    private static IBrush Label => Tokens.Brush("DimB");
    private static IPen Curve => new Pen(Tokens.Brush("AccB"), 2, lineJoin: PenLineJoin.Round);
    private static IBrush Dot => Tokens.Brush("SelB");
    private static IPen RefPen => new Pen(Tokens.Brush("DimB"), 1, new DashStyle([4, 4], 0));

    public override void Render(DrawingContext ctx)
    {
        var b = Bounds;
        const double L = 44, R = 12, T = 18, B = 24;
        var w = b.Width - L - R;
        var h = b.Height - T - B;
        if (w < 20 || h < 20) return;
        var tf = new Typeface(Tokens.Mono);
        void Text(string s, double x, double y, bool right = false, bool centre = false)
        {
            var ft = new FormattedText(s, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, tf, 10, Label);
            ctx.DrawText(ft, new Point(right ? x - ft.Width : centre ? x - ft.Width / 2 : x, y - ft.Height / 2));
        }

        if (_data.Length < 2)
        {
            Text("no data", L + w / 2, T + h / 2, centre: true);
            return;
        }
        double xmin = 0, xmax, ymin = 0, ymax;
        if (AutoRange)
        {
            var all = _data.Concat(_overlay).ToArray();
            xmin = all.Min(p => p.X);
            xmax = all.Max(p => p.X);
            if (xmax - xmin < 1e-9) xmax = xmin + 1;
            ymin = all.Min(p => p.Y);
            ymax = all.Max(p => p.Y);
            if (RefY is double r) { ymin = Math.Min(ymin, r); ymax = Math.Max(ymax, r); }
            var pad = Math.Max(1e-9, (ymax - ymin) * 0.08);
            ymin -= pad;
            ymax += pad;
        }
        else
        {
            xmax = _data[^1].X + (_data[1].X - _data[0].X) / 2;
            ymax = Math.Max(1.2, Math.Ceiling(_data.Max(p => p.Y) * 2) / 2);
        }
        double X(double v) => L + (v - xmin) / (xmax - xmin) * w;
        double Y(double v) => T + h - (v - ymin) / (ymax - ymin) * h;
        string Fmt(double v) => Math.Abs(v) >= 1000 ? v.ToString("0", CultureInfo.InvariantCulture) : v.ToString("0.##", CultureInfo.InvariantCulture);

        foreach (var t in new[] { ymin, (ymin + ymax) / 2, ymax })
        {
            ctx.DrawLine(new Pen(Grid, 1), new Point(L, Y(t)), new Point(L + w, Y(t)));
            Text(Fmt(t), L - 6, Y(t), right: true);
        }
        if (AutoRange)
        {
            var xf = xmax - xmin >= 20 ? "0" : xmax - xmin >= 2 ? "0.#" : "0.###";
            Text(xmin.ToString(xf, CultureInfo.InvariantCulture), L, T + h + 12);
            Text(((xmin + xmax) / 2).ToString(xf, CultureInfo.InvariantCulture), X((xmin + xmax) / 2), T + h + 12, centre: true);
            Text(xmax.ToString(xf, CultureInfo.InvariantCulture), L + w, T + h + 12, right: true);
        }
        else
        {
            for (var t = 0.0; t <= xmax + 1e-9; t += xmax > 8 ? 2 : 1)
                if (X(t) < L + w - 40) Text(t.ToString("0", CultureInfo.InvariantCulture), X(t), T + h + 12, centre: true);
        }
        ctx.DrawLine(new Pen(Axis, 1), new Point(L, T + h), new Point(L + w, T + h));
        if (RefY is double rv && rv >= ymin && rv <= ymax) ctx.DrawLine(RefPen, new Point(L, Y(rv)), new Point(L + w, Y(rv)));

        StreamGeometry Line((double X, double Y)[] d)
        {
            var geo = new StreamGeometry();
            using var g = geo.Open();
            g.BeginFigure(new Point(X(d[0].X), Y(d[0].Y)), false);
            foreach (var p in d.Skip(1)) g.LineTo(new Point(X(p.X), Y(p.Y)));
            g.EndFigure(false);
            return geo;
        }
        if (Markers)
        {
            foreach (var p in _data) ctx.DrawEllipse(Dot, null, new Point(X(p.X), Y(p.Y)), 2.2, 2.2);
            if (_overlay.Length > 1) ctx.DrawGeometry(null, Curve, Line(_overlay));
        }
        else
        {
            ctx.DrawGeometry(null, Curve, Line(_data));
        }
        Text(YLabel, L, T - 9);
        if (AutoRange) Text("x: " + XLabel, L + w, T - 9, right: true);
        else Text(XLabel, L + w, T + h + 12 + 0, right: true);
    }
}
