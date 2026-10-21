using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;
using Avalonia.Styling;

namespace CapsStudio.Views;

/// <summary>How one series of a SciFigure is drawn.</summary>
public enum SciKind { Line, Markers, LineMarkers, Dashed, Dotted, Steps }

/// <summary>One series of a figure. Colour is an index into the figure palette (−1: the ink colour). Err: ± per point.</summary>
public sealed record SciSeries(string Name, double[] X, double[] Y)
{
    public SciKind Kind { get; init; } = SciKind.Line;
    public double[]? Err { get; init; }
    public int Colour { get; init; }
    public bool Open { get; init; }
    public bool Square { get; init; }
    public bool Right { get; init; }
    public bool Legend { get; init; } = true;
    public double Width { get; init; } = 1.6;
    /// <summary>Fill the area under the curve (a histogram drawn as steps), lightly.</summary>
    public bool Fill { get; init; }
}

/// <summary>A vertical marker: a dashed line at X with a label at its top (Tg, a peak, a fit boundary).</summary>
public sealed record SciMark(double X, string Label);

/// <summary>A scientific figure: axis titles with units (a small markup: *italic*, _{sub}, ^{sup}), linear or log axes,
/// series with markers and error bars, a reference line, a horizontal band (mean ± standard error), a shaded x window (a
/// fit range), vertical marks and notes printed in the legend.</summary>
public sealed record SciFigure
{
    public string XLabel { get; init; } = "";
    public string YLabel { get; init; } = "";
    public string Y2Label { get; init; } = "";
    public bool XLog { get; init; }
    public bool YLog { get; init; }
    public SciSeries[] Series { get; init; } = [];
    public double? RefY { get; init; }
    public (double Lo, double Hi, string Label)? YBand { get; init; }
    public (double From, double To, string Label)? XBand { get; init; }
    public SciMark[] Marks { get; init; } = [];
    public string[] Notes { get; init; } = [];
    public double? XMin { get; init; }
    public double? XMax { get; init; }
    public double? YMin { get; init; }
    public double? YMax { get; init; }
    /// <summary>The y axis starts at zero when the data come near it (distributions, g(r), MSD).</summary>
    public bool YZero { get; init; }
    public string Empty { get; init; } = "no data";
}

/// <summary>Draws a SciFigure the way a journal figure looks: a closed frame, inward major and minor ticks on all four
/// sides, round tick values, axis titles as "quantity / unit" (the y title turned), a legend in the emptiest corner,
/// markers with error bars, in the colours of the current theme. Compact draws the small version on a result card.</summary>
public sealed class SciPlot : Control
{
    public static readonly StyledProperty<SciFigure?> FigureProperty = AvaloniaProperty.Register<SciPlot, SciFigure?>(nameof(Figure));
    public static readonly StyledProperty<bool> CompactProperty = AvaloniaProperty.Register<SciPlot, bool>(nameof(Compact));
    public SciFigure? Figure { get => GetValue(FigureProperty); set => SetValue(FigureProperty, value); }
    public bool Compact { get => GetValue(CompactProperty); set => SetValue(CompactProperty, value); }
    /// <summary>A click inside the frame, with the x it lands on in data units.</summary>
    public event Action<double>? Clicked;

    static SciPlot() => AffectsRender<SciPlot>(FigureProperty, CompactProperty);

    private Rect _frame;
    private Func<double, double>? _xinv;

    protected override void OnPointerPressed(Avalonia.Input.PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        var p = e.GetPosition(this);
        if (_xinv != null && _frame.Contains(p)) Clicked?.Invoke(_xinv(p.X));
    }

    // ------------------------------------------------------------------ colours
    private static bool Light => Application.Current?.ActualThemeVariant == ThemeVariant.Light;
    // Okabe & Ito's colour-blind safe set, lightened on the dark theme so every colour reads on the background
    private static readonly string[] DarkPalette = ["#5AB4E5", "#F0A030", "#3CC49A", "#E07CB4", "#E8DA5A", "#F07050"];
    private static readonly string[] LightPalette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#B8860B", "#E69F00"];
    private static IBrush Ink => Tokens.Brush("TextB");
    private static IBrush Soft => Tokens.Brush("MutedB");
    internal static IBrush Colour(int k) => k < 0 ? Ink : new SolidColorBrush(Color.Parse((Light ? LightPalette : DarkPalette)[k % DarkPalette.Length]));
    private static Color ColourC(int k) => k < 0 ? ((Ink as ISolidColorBrush)?.Color ?? Colors.White) : Color.Parse((Light ? LightPalette : DarkPalette)[k % DarkPalette.Length]);

    // ------------------------------------------------------------------ axes
    private sealed class Axis
    {
        public double Lo, Hi;
        public bool Log;
        public double Step;            // major step (decades on a log axis)
        public List<double> Major = [], Minor = [];
        public double Map(double v, double a, double b) => Log ? a + (Math.Log10(v) - Lo) / (Hi - Lo) * (b - a) : a + (v - Lo) / (Hi - Lo) * (b - a);
        public double Unmap(double p, double a, double b) { var t = Lo + (p - a) / (b - a) * (Hi - Lo); return Log ? Math.Pow(10, t) : t; }
        public string Label(double v)
        {
            if (Log) return Pow10(v);
            if (Math.Abs(v) < Step * 1e-6) return "0";
            var dec = 0;
            while (dec < 6 && Math.Abs(Step * Math.Pow(10, dec) - Math.Round(Step * Math.Pow(10, dec))) > 1e-6 * Math.Pow(10, dec)) dec++;
            var a = Math.Abs(v);
            if (a >= 1e5 || (a < 1e-3 && a > 0)) return v.ToString("0.#E0", CultureInfo.InvariantCulture).Replace("-", "−");
            return v.ToString(dec == 0 ? "0" : "0." + new string('0', dec), CultureInfo.InvariantCulture).Replace("-", "−");
        }
    }

    private static string Pow10(double e)
    {
        var n = (int)Math.Round(e);
        if (n == 0) return "1";
        if (n == 1) return "10";
        const string sup = "⁰¹²³⁴⁵⁶⁷⁸⁹";
        var s = Math.Abs(n).ToString(CultureInfo.InvariantCulture);
        return "10" + (n < 0 ? "⁻" : "") + string.Concat(s.Select(c => sup[c - '0']));
    }

    private static double NiceStep(double span, int target)
    {
        var raw = span / Math.Max(1, target);
        var mag = Math.Pow(10, Math.Floor(Math.Log10(Math.Max(1e-300, raw))));
        var f = raw / mag;
        return (f < 1.5 ? 1 : f < 2.25 ? 2 : f < 3.5 ? 2.5 : f < 7.5 ? 5 : 10) * mag;
    }

    /// <summary>A linear axis over [lo, hi] with about `target` major ticks; the ends snap to the minor step when that wastes
    /// little room, else stay at a small padding past the data.</summary>
    private static Axis Linear(double lo, double hi, int target, bool pad, double? fixLo, double? fixHi)
    {
        if (!(hi > lo)) { var c = double.IsFinite(lo) ? lo : 0; var d = Math.Max(Math.Abs(c) * 0.05, 1e-6); lo = c - d; hi = c + d; }
        var span = hi - lo;
        if (pad) { lo -= span * 0.04; hi += span * 0.04; }
        var step = NiceStep(hi - lo, target);
        var minor = step / (Math.Abs(step / Math.Pow(10, Math.Floor(Math.Log10(step))) - 2) < 1e-9 ? 4 : 5);
        var sLo = Math.Floor(lo / minor + 1e-9) * minor;
        var sHi = Math.Ceiling(hi / minor - 1e-9) * minor;
        var ax = new Axis { Lo = fixLo ?? sLo, Hi = fixHi ?? sHi, Step = step };
        for (var t = Math.Ceiling(ax.Lo / step - 1e-9) * step; t <= ax.Hi + step * 1e-9; t += step) ax.Major.Add(Math.Abs(t) < step * 1e-9 ? 0 : t);
        for (var t = Math.Ceiling(ax.Lo / minor - 1e-9) * minor; t <= ax.Hi + minor * 1e-9; t += minor) ax.Minor.Add(t);
        return ax;
    }

    private static Axis Logarithmic(double lo, double hi, double? fixLo, double? fixHi)
    {
        lo = Math.Max(lo, 1e-300);
        var a = Math.Floor(Math.Log10(fixLo ?? lo) + 1e-9);
        var b = Math.Ceiling(Math.Log10(fixHi ?? hi) - 1e-9);
        if (b <= a) b = a + 1;
        var ax = new Axis { Log = true, Lo = fixLo is double fl ? Math.Log10(fl) : a, Hi = fixHi is double fh ? Math.Log10(fh) : b, Step = 1 };
        var every = Math.Max(1, (int)Math.Ceiling((b - a) / 6));   // at most about six decade labels
        for (var e = a; e <= b + 1e-9; e += 1)
        {
            if (((int)(e - a)) % every == 0) ax.Major.Add(e);
            for (var m = 2; m <= 9; ++m) { var v = e + Math.Log10(m); if (v <= ax.Hi + 1e-9 && v >= ax.Lo - 1e-9) ax.Minor.Add(v); }
            if (every > 1 && ((int)(e - a)) % every != 0) ax.Minor.Add(e);
        }
        return ax;
    }

    // ------------------------------------------------------------------ text with a small markup
    private sealed record Run(string Text, bool Italic, int Shift);   // Shift: −1 subscript, +1 superscript

    private static List<Run> Parse(string s)
    {
        var runs = new List<Run>();
        var italic = false;
        var buf = new System.Text.StringBuilder();
        void Flush(int shift = 0) { if (buf.Length > 0) { runs.Add(new Run(buf.ToString(), italic, shift)); buf.Clear(); } }
        for (var i = 0; i < s.Length; ++i)
        {
            var c = s[i];
            if (c == '*') { Flush(); italic = !italic; continue; }
            if ((c == '_' || c == '^') && i + 1 < s.Length && s[i + 1] == '{')
            {
                Flush();
                var end = s.IndexOf('}', i + 2);
                if (end < 0) end = s.Length;
                var inner = s.Substring(i + 2, end - i - 2);
                var it = inner.StartsWith('*') && inner.EndsWith('*') && inner.Length > 2;
                runs.Add(new Run(it ? inner[1..^1] : inner, it, c == '_' ? -1 : 1));
                i = end;
                continue;
            }
            buf.Append(c);
        }
        Flush();
        return runs;
    }

    private sealed record Laid(FormattedText Ft, double Dx, double Dy);

    /// <summary>The markup laid out from (0, 0) (top left of the line box); returns the pieces, the width and the height.</summary>
    private static (List<Laid> Parts, double W, double H) Layout(string s, double size, IBrush brush, FontWeight weight = FontWeight.Normal)
    {
        var parts = new List<Laid>();
        double x = 0, h = 0;
        foreach (var r in Parse(s))
        {
            var sz = r.Shift == 0 ? size : size * 0.72;
            var tf = new Typeface(Tokens.Sans, r.Italic ? FontStyle.Italic : FontStyle.Normal, weight);
            var ft = new FormattedText(r.Text, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, tf, sz, brush);
            var dy = r.Shift == 0 ? 0 : r.Shift < 0 ? size * 0.42 : -size * 0.12;
            parts.Add(new Laid(ft, x, dy));
            x += ft.WidthIncludingTrailingWhitespace + (r.Italic ? size * 0.06 : 0);
            h = Math.Max(h, size * 1.25);
        }
        return (parts, x, h);
    }

    private static void DrawRich(DrawingContext ctx, string s, double size, IBrush brush, Point at, double alignX = 0, double alignY = 0, FontWeight weight = FontWeight.Normal)
    {
        var (parts, w, h) = Layout(s, size, brush, weight);
        var o = new Point(at.X - w * alignX, at.Y - h * alignY);
        foreach (var p in parts) ctx.DrawText(p.Ft, new Point(o.X + p.Dx, o.Y + p.Dy));
    }

    private static double RichWidth(string s, double size) => Layout(s, size, Brushes.White).W;

    // ------------------------------------------------------------------ render
    public override void Render(DrawingContext ctx)
    {
        var fig = Figure;
        var b = Bounds;
        var compact = Compact;
        double fTick = compact ? 9.5 : 11.5, fLabel = compact ? 10.5 : 13, fLegend = compact ? 9 : 11;
        var tickLen = compact ? 4.0 : 6.0;
        var series = fig?.Series.Where(s => s.X.Length > 0).ToArray() ?? [];
        bool Ok(double x, double y, bool right) => double.IsFinite(x) && double.IsFinite(y) && (!fig!.XLog || x > 0) && (right || !fig.YLog || y > 0);
        var left = series.Where(s => !s.Right).SelectMany(s => s.X.Zip(s.Y, (x, y) => (x, y, e: 0.0)).Select((p, i) => (p.x, p.y, e: s.Err is { } er && i < er.Length && double.IsFinite(er[i]) ? er[i] : 0)))
            .Where(p => Ok(p.x, p.y, false)).ToArray();
        if (fig == null || left.Length < 2)
        {
            var ft = new FormattedText(fig?.Empty ?? "no data", CultureInfo.InvariantCulture, FlowDirection.LeftToRight, new Typeface(Tokens.Sans), fTick, Soft);
            ctx.DrawText(ft, new Point((b.Width - ft.Width) / 2, (b.Height - ft.Height) / 2));
            _xinv = null;
            return;
        }
        var rightPts = series.Where(s => s.Right).SelectMany(s => s.X.Zip(s.Y)).Where(p => Ok(p.First, p.Second, true)).ToArray();
        var hasRight = rightPts.Length >= 2;

        // ranges
        double xlo = left.Min(p => p.x), xhi = left.Max(p => p.x);
        if (hasRight) { xlo = Math.Min(xlo, rightPts.Min(p => p.First)); xhi = Math.Max(xhi, rightPts.Max(p => p.First)); }
        double ylo = left.Min(p => p.y - p.e), yhi = left.Max(p => p.y + p.e);
        if (fig.RefY is double ry && !fig.YLog) { ylo = Math.Min(ylo, ry); yhi = Math.Max(yhi, ry); }
        if (fig.YBand is { } yb) { ylo = Math.Min(ylo, yb.Lo); yhi = Math.Max(yhi, yb.Hi); }
        if (fig.YZero && !fig.YLog && ylo >= 0 && ylo < 0.3 * (yhi - ylo)) ylo = 0;
        var nX = compact ? 4 : Math.Clamp((int)(b.Width / 110), 4, 8);
        var nY = compact ? 4 : Math.Clamp((int)(b.Height / 60), 3, 7);
        var xa = fig.XLog ? Logarithmic(xlo, xhi, fig.XMin, fig.XMax) : Linear(xlo, xhi, nX, false, fig.XMin, fig.XMax);
        var ya = fig.YLog ? Logarithmic(ylo, yhi, fig.YMin, fig.YMax) : Linear(ylo, yhi, nY, !(fig.YZero && ylo == 0), fig.YMin ?? (fig.YZero && ylo == 0 ? 0 : null), fig.YMax);
        Axis? y2 = hasRight ? Linear(Math.Min(0, rightPts.Min(p => p.Second)), rightPts.Max(p => p.Second), nY, false, null, null) : null;

        // margins from the tick labels and titles
        var yTickW = ya.Major.Select(v => RichWidth(ya.Label(v), fTick)).DefaultIfEmpty(0).Max();
        var y2TickW = y2?.Major.Select(v => RichWidth(y2.Label(v), fTick)).DefaultIfEmpty(0).Max() ?? 0;
        var gap = compact ? 3.0 : 5.0;
        var L = 4 + fLabel * 1.3 + gap + yTickW + gap;
        var R = hasRight ? 4 + fLabel * 1.3 + gap + y2TickW + gap : compact ? 8 : 14;
        var T = compact ? 6.0 : 10.0;
        var B = gap + fTick * 1.3 + gap + fLabel * 1.35 + 2;
        var w = b.Width - L - R;
        var h = b.Height - T - B;
        if (w < 40 || h < 30) return;
        _frame = new Rect(L, T, w, h);
        double PX(double v) => xa.Map(v, L, L + w);
        double PY(double v) => ya.Map(v, T + h, T);
        double PY2(double v) => y2!.Map(v, T + h, T);
        _xinv = px => xa.Unmap(px, L, L + w);

        var ink = Ink;
        var framePen = new Pen(ink, compact ? 0.9 : 1.1);
        var tickPen = new Pen(ink, compact ? 0.8 : 1.0);
        var bg = Tokens.Brush("Bg1B");

        // shaded x window and horizontal band under everything
        if (fig.XBand is { } xb)
        {
            var x0 = Math.Clamp(PX(xb.From), L, L + w);
            var x1 = Math.Clamp(PX(xb.To), L, L + w);
            if (x1 > x0) ctx.FillRectangle(new SolidColorBrush(ColourC(0), Light ? 0.10 : 0.13), new Rect(x0, T, x1 - x0, h));
        }
        if (fig.YBand is { } band)
        {
            var y0 = Math.Clamp(PY(band.Hi), T, T + h);
            var y1 = Math.Clamp(PY(band.Lo), T, T + h);
            ctx.FillRectangle(new SolidColorBrush(ColourC(1), Light ? 0.16 : 0.22), new Rect(L, y0, w, Math.Max(1, y1 - y0)));
            ctx.DrawLine(new Pen(Colour(1), 1.2), new Point(L, (y0 + y1) / 2), new Point(L + w, (y0 + y1) / 2));
        }
        if (fig.RefY is double rv && (!fig.YLog || rv > 0))
        {
            var py = PY(rv);
            if (py >= T && py <= T + h) ctx.DrawLine(new Pen(Soft, 0.9, new DashStyle([4, 3], 0)), new Point(L, py), new Point(L + w, py));
        }

        // data, clipped to the frame
        using (ctx.PushClip(new Rect(L, T, w, h)))
        {
            foreach (var s in series)
            {
                var brush = Colour(s.Colour);
                var pts = new List<Point>();
                var errs = new List<(Point Lo, Point Hi)>();
                for (var i = 0; i < Math.Min(s.X.Length, s.Y.Length); ++i)
                {
                    if (!Ok(s.X[i], s.Y[i], s.Right)) continue;
                    var px = PX(s.X[i]);
                    var py = s.Right ? PY2(s.Y[i]) : PY(s.Y[i]);
                    pts.Add(new Point(px, py));
                    if (s.Err is { } er && i < er.Length && er[i] > 0 && double.IsFinite(er[i]) && !s.Right)
                        errs.Add((new Point(px, PY(fig.YLog ? Math.Max(s.Y[i] - er[i], Math.Pow(10, ya.Lo)) : s.Y[i] - er[i])), new Point(px, PY(s.Y[i] + er[i]))));
                }
                if (pts.Count == 0) continue;
                var lw = compact ? Math.Max(1, s.Width * 0.75) : s.Width;
                IPen? pen = s.Kind switch
                {
                    SciKind.Line or SciKind.LineMarkers or SciKind.Steps => new Pen(brush, lw, lineJoin: PenLineJoin.Round),
                    SciKind.Dashed => new Pen(brush, lw, new DashStyle([5, 3], 0), lineCap: PenLineCap.Flat),
                    SciKind.Dotted => new Pen(brush, lw, new DashStyle([1, 2.5], 0), lineCap: PenLineCap.Round),
                    _ => null,
                };
                if (pen != null && pts.Count > 1)
                {
                    var geo = new StreamGeometry();
                    using (var g = geo.Open())
                    {
                        if (s.Kind == SciKind.Steps)
                        {
                            var half = pts.Count > 1 ? (pts[1].X - pts[0].X) / 2 : 0;
                            g.BeginFigure(new Point(pts[0].X - half, PY(Math.Max(ya.Log ? Math.Pow(10, ya.Lo) : 0, ya.Log ? Math.Pow(10, ya.Lo) : ya.Lo))), s.Fill);
                            foreach (var p in pts) { g.LineTo(new Point(p.X - half, p.Y)); g.LineTo(new Point(p.X + half, p.Y)); }
                            g.LineTo(new Point(pts[^1].X + half, PY(ya.Log ? Math.Pow(10, ya.Lo) : Math.Max(0, ya.Lo))));
                            g.EndFigure(s.Fill);
                        }
                        else
                        {
                            g.BeginFigure(pts[0], false);
                            foreach (var p in pts.Skip(1)) g.LineTo(p);
                            g.EndFigure(false);
                        }
                    }
                    if (s.Fill) ctx.DrawGeometry(new SolidColorBrush(ColourC(s.Colour), Light ? 0.18 : 0.22), null, geo);
                    ctx.DrawGeometry(null, pen, geo);
                }
                if (s.Kind is SciKind.Markers or SciKind.LineMarkers)
                {
                    var ep = new Pen(brush, compact ? 0.8 : 1.0);
                    var cap = compact ? 2.0 : 3.0;
                    foreach (var (lo, hi) in errs)
                    {
                        ctx.DrawLine(ep, lo, hi);
                        ctx.DrawLine(ep, new Point(lo.X - cap, lo.Y), new Point(lo.X + cap, lo.Y));
                        ctx.DrawLine(ep, new Point(hi.X - cap, hi.Y), new Point(hi.X + cap, hi.Y));
                    }
                    var r = compact ? 2.3 : pts.Count > 80 ? 2.2 : 3.4;
                    var stride = Math.Max(1, pts.Count / (compact ? 60 : 160));   // a dense curve keeps a readable marker spacing
                    for (var i = 0; i < pts.Count; i += stride) Marker(ctx, pts[i], r, brush, s.Open, s.Square, bg, compact);
                }
            }
        }
        // vertical marks
        foreach (var m in fig.Marks)
        {
            if (!double.IsFinite(m.X) || (fig.XLog && m.X <= 0)) continue;
            var px = PX(m.X);
            if (px < L || px > L + w) continue;
            ctx.DrawLine(new Pen(ink, compact ? 0.8 : 1.0, new DashStyle([3, 3], 0)), new Point(px, T), new Point(px, T + h));
            if (m.Label.Length > 0)
            {
                var lw = RichWidth(m.Label, fLegend);
                var right = px + 4 + lw < L + w - 2;
                DrawRich(ctx, m.Label, fLegend, ink, new Point(right ? px + 4 : px - 4, T + 3), right ? 0 : 1, 0);
            }
        }

        // frame and ticks (inward, on all four sides; the right side carries the second axis when there is one)
        ctx.DrawRectangle(null, framePen, new Rect(L, T, w, h));
        foreach (var v in xa.Minor) { var px = PX(xa.Log ? Math.Pow(10, v) : v); if (px < L - 0.5 || px > L + w + 0.5) continue; ctx.DrawLine(tickPen, new Point(px, T + h), new Point(px, T + h - tickLen * 0.55)); ctx.DrawLine(tickPen, new Point(px, T), new Point(px, T + tickLen * 0.55)); }
        foreach (var v in xa.Major)
        {
            var px = PX(xa.Log ? Math.Pow(10, v) : v);
            if (px < L - 0.5 || px > L + w + 0.5) continue;
            ctx.DrawLine(tickPen, new Point(px, T + h), new Point(px, T + h - tickLen));
            ctx.DrawLine(tickPen, new Point(px, T), new Point(px, T + tickLen));
            DrawRich(ctx, xa.Label(v), fTick, ink, new Point(px, T + h + gap), 0.5, 0);
        }
        foreach (var v in ya.Minor)
        {
            var py = PY(ya.Log ? Math.Pow(10, v) : v);
            if (py < T - 0.5 || py > T + h + 0.5) continue;
            ctx.DrawLine(tickPen, new Point(L, py), new Point(L + tickLen * 0.55, py));
            if (!hasRight) ctx.DrawLine(tickPen, new Point(L + w, py), new Point(L + w - tickLen * 0.55, py));
        }
        foreach (var v in ya.Major)
        {
            var py = PY(ya.Log ? Math.Pow(10, v) : v);
            if (py < T - 0.5 || py > T + h + 0.5) continue;
            ctx.DrawLine(tickPen, new Point(L, py), new Point(L + tickLen, py));
            if (!hasRight) ctx.DrawLine(tickPen, new Point(L + w, py), new Point(L + w - tickLen, py));
            DrawRich(ctx, ya.Label(v), fTick, ink, new Point(L - gap, py), 1, 0.5);
        }
        if (y2 != null)
        {
            foreach (var v in y2.Minor) { var py = PY2(v); if (py >= T - 0.5 && py <= T + h + 0.5) ctx.DrawLine(tickPen, new Point(L + w, py), new Point(L + w - tickLen * 0.55, py)); }
            foreach (var v in y2.Major)
            {
                var py = PY2(v);
                if (py < T - 0.5 || py > T + h + 0.5) continue;
                ctx.DrawLine(tickPen, new Point(L + w, py), new Point(L + w - tickLen, py));
                DrawRich(ctx, y2.Label(v), fTick, ink, new Point(L + w + gap, py), 0, 0.5);
            }
        }
        // axis titles
        DrawRich(ctx, fig.XLabel, fLabel, ink, new Point(L + w / 2, T + h + gap + fTick * 1.3 + gap), 0.5, 0);
        using (ctx.PushTransform(Matrix.CreateRotation(-Math.PI / 2) * Matrix.CreateTranslation(4, T + h / 2)))
            DrawRich(ctx, fig.YLabel, fLabel, ink, new Point(0, 0), 0.5, 0);
        if (y2 != null && fig.Y2Label.Length > 0)
            using (ctx.PushTransform(Matrix.CreateRotation(Math.PI / 2) * Matrix.CreateTranslation(b.Width - 4, T + h / 2)))
                DrawRich(ctx, fig.Y2Label, fLabel, ink, new Point(0, 0), 0.5, 0);

        // legend in the corner with the fewest data points under it
        var entries = series.Where(s => s.Legend && s.Name.Length > 0).ToArray();
        var notes = fig.Notes.Where(n => n.Length > 0).ToArray();
        if (fig.YBand is { Label.Length: > 0 } ybl) entries = [.. entries, new SciSeries(ybl.Label, [], []) { Colour = 1, Kind = SciKind.Line, Width = 6 }];
        if (fig.XBand is { Label.Length: > 0 } xbl) entries = [.. entries, new SciSeries(xbl.Label, [], []) { Colour = 0, Kind = SciKind.Line, Width = 6 }];
        if (entries.Length == 0 && notes.Length == 0) return;
        if (compact) return;   // a card's figure is too small for a legend: the card says what it shows
        var sample = compact ? 14.0 : 22.0;
        var rowH = fLegend * 1.45;
        var pad = compact ? 4.0 : 6.0;
        var lw2 = Math.Max(entries.Select(e => sample + 5 + RichWidth(e.Name, fLegend)).DefaultIfEmpty(0).Max(), notes.Select(n => RichWidth(n, fLegend)).DefaultIfEmpty(0).Max());
        var boxW = Math.Min(w - 8, lw2 + 2 * pad);
        var boxH = (entries.Length + notes.Length) * rowH + 2 * pad - (rowH - fLegend * 1.2);
        if (boxH > h - 8) return;
        var corners = new[] { new Point(L + w - boxW - 6, T + 6), new Point(L + 6, T + 6), new Point(L + w - boxW - 6, T + h - boxH - 6), new Point(L + 6, T + h - boxH - 6) };
        var allPts = series.SelectMany(s => s.X.Zip(s.Y, (x, y) => (x, y, s.Right))).Where(p => Ok(p.x, p.y, p.Right))
            .Select(p => new Point(PX(p.x), p.Right ? PY2(p.y) : PY(p.y))).ToArray();
        var markLines = fig.Marks.Where(m => double.IsFinite(m.X) && (!fig.XLog || m.X > 0)).Select(m => PX(m.X)).ToArray();
        int Cost(Point c)
        {
            var r = new Rect(c, new Size(boxW, boxH)).Inflate(3);
            return allPts.Count(r.Contains) + markLines.Count(x => x >= r.Left && x <= r.Right) * 25;
        }
        var at = corners.OrderBy(Cost).First();
        var box = new Rect(at, new Size(boxW, boxH));
        ctx.FillRectangle(new SolidColorBrush((bg as ISolidColorBrush)?.Color ?? Colors.Black, 0.88), box);
        ctx.DrawRectangle(null, new Pen(Tokens.Brush("LineB"), 0.8), box);
        var y = box.Top + pad + fLegend * 0.6;
        foreach (var e in entries)
        {
            var x0 = box.Left + pad;
            var brush = Colour(e.Colour);
            if (e.Width >= 5) ctx.FillRectangle(new SolidColorBrush(ColourC(e.Colour), Light ? 0.25 : 0.3), new Rect(x0, y - fLegend * 0.35, sample, fLegend * 0.7));
            else if (e.Kind is SciKind.Line or SciKind.LineMarkers or SciKind.Steps or SciKind.Dashed or SciKind.Dotted)
            {
                IPen pen = e.Kind == SciKind.Dashed ? new Pen(brush, Math.Min(2, e.Width), new DashStyle([5, 3], 0))
                    : e.Kind == SciKind.Dotted ? new Pen(brush, Math.Min(2, e.Width), new DashStyle([1, 2.5], 0), lineCap: PenLineCap.Round)
                    : new Pen(brush, Math.Min(2, e.Width));
                ctx.DrawLine(pen, new Point(x0, y), new Point(x0 + sample, y));
            }
            if (e.Kind is SciKind.Markers or SciKind.LineMarkers) Marker(ctx, new Point(x0 + sample / 2, y), compact ? 2.3 : 3.2, brush, e.Open, e.Square, bg, compact);
            DrawRich(ctx, e.Name, fLegend, ink, new Point(x0 + sample + 5, y), 0, 0.5);
            y += rowH;
        }
        foreach (var n in notes) { DrawRich(ctx, n, fLegend, ink, new Point(box.Left + pad, y), 0, 0.5); y += rowH; }
    }

    private static void Marker(DrawingContext ctx, Point p, double r, IBrush brush, bool open, bool square, IBrush bg, bool compact)
    {
        var pen = new Pen(brush, compact ? 0.9 : 1.2);
        if (square)
        {
            var rect = new Rect(p.X - r, p.Y - r, 2 * r, 2 * r);
            ctx.DrawRectangle(open ? bg : brush, pen, rect);
        }
        else ctx.DrawEllipse(open ? bg : brush, pen, p, r, r);
    }
}
