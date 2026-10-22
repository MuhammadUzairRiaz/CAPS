using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;
using Avalonia.Styling;

namespace CapsStudio.Views;

/// <summary>How one series of a SciFigure is drawn. Bars: one bar per point (grouped when several bar series share x).</summary>
public enum SciKind { Line, Markers, LineMarkers, Dashed, Dotted, Steps, Bars }

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

    /// <summary>A series from (x, y) pairs.</summary>
    public static SciSeries Of(string name, IReadOnlyList<(double X, double Y)> pts) => new(name, pts.Select(p => p.X).ToArray(), pts.Select(p => p.Y).ToArray());
}

/// <summary>A vertical marker: a dashed line at X with a label at its top (Tg, a peak, a fit boundary).</summary>
public sealed record SciMark(double X, string Label);
/// <summary>A horizontal guide: a dashed line at Y with a label at its right end (a target, a tolerance).</summary>
public sealed record SciGuide(double Y, string Label);
/// <summary>A labelled point (a critical point, the current state).</summary>
public sealed record SciPoint(double X, double Y, string Label, int Colour = -1);
/// <summary>A stretch of the x axis (a stage of a protocol), tinted when Shade, its label written along the top.</summary>
public sealed record SciStage(double From, double To, string Label, bool Shade);
/// <summary>A map over two axes: one cell per (x, y) on a regular grid, coloured by z.</summary>
public sealed record SciHeat(double[] X, double[] Y, double[] Z, string ZLabel);

/// <summary>A scientific figure: axis titles with units (a small markup: *italic*, _{sub}, ^{sup}), linear or log axes,
/// series with markers and error bars, a reference line, a horizontal band (mean ± standard error), a shaded x window (a
/// fit range), vertical marks, horizontal guides, stages along x, labelled points, categories along x, a heat map, and
/// notes printed in the legend.</summary>
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
    public SciGuide[] Guides { get; init; } = [];
    public SciPoint[] Points { get; init; } = [];
    public SciStage[] Stages { get; init; } = [];
    /// <summary>Names along x at 0, 1, 2, … (a bar chart over categories).</summary>
    public string[] Categories { get; init; } = [];
    public SciHeat? Heat { get; init; }
    public string[] Notes { get; init; } = [];
    public double? XMin { get; init; }
    public double? XMax { get; init; }
    public double? YMin { get; init; }
    public double? YMax { get; init; }
    public double? Y2Min { get; init; }
    /// <summary>The y axis starts at zero when the data come near it (distributions, g(r), MSD).</summary>
    public bool YZero { get; init; }
    /// <summary>The x axis counts whole things (blocks, cycles): its major step is never below 1.</summary>
    public bool XWhole { get; init; }
    /// <summary>The y axis counts whole things (links, bridges, beads): its major step is never below 1.</summary>
    public bool YWhole { get; init; }
    public string Empty { get; init; } = "no data";
}

/// <summary>Draws a SciFigure the way a journal figure looks: a closed frame, inward major and minor ticks on all four
/// sides, round tick values, axis titles as "quantity / unit" (the y title turned), a legend in the emptiest corner,
/// markers with error bars, in the colours of the current theme. Compact draws the small version on a result card.
/// A long series (a live monitor's thousands of rows) is drawn through its extremes per pixel column, so a redraw stays
/// cheap however long the run.</summary>
public sealed class SciPlot : Control
{
    public static readonly StyledProperty<SciFigure?> FigureProperty = AvaloniaProperty.Register<SciPlot, SciFigure?>(nameof(Figure));
    public static readonly StyledProperty<bool> CompactProperty = AvaloniaProperty.Register<SciPlot, bool>(nameof(Compact));
    public static readonly StyledProperty<string> EmptyTextProperty = AvaloniaProperty.Register<SciPlot, string>(nameof(EmptyText), "no data");
    public SciFigure? Figure { get => GetValue(FigureProperty); set => SetValue(FigureProperty, value); }
    public bool Compact { get => GetValue(CompactProperty); set => SetValue(CompactProperty, value); }
    /// <summary>What the plot says before it has a figure.</summary>
    public string EmptyText { get => GetValue(EmptyTextProperty); set => SetValue(EmptyTextProperty, value); }
    /// <summary>A click inside the frame, with the x it lands on in data units.</summary>
    public event Action<double>? Clicked;
    /// <summary>Allow a drag box; Brushed gets (x0, x1, y0, y1) in data units, NaN when cleared by a click.</summary>
    public bool Brushable { get; set; }
    public event Action<double, double, double, double>? Brushed;

    static SciPlot() => AffectsRender<SciPlot>(FigureProperty, CompactProperty, EmptyTextProperty);

    private Rect _frame;
    private Func<double, double>? _xinv, _yinv;
    private Point? _b0, _b1;

    protected override void OnPointerPressed(Avalonia.Input.PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        var p = e.GetPosition(this);
        if (_xinv == null || !_frame.Contains(p)) return;
        if (!Brushable) { Clicked?.Invoke(_xinv(p.X)); return; }
        _b0 = _b1 = p;
        e.Pointer.Capture(this);
        InvalidateVisual();
    }
    protected override void OnPointerMoved(Avalonia.Input.PointerEventArgs e)
    {
        base.OnPointerMoved(e);
        if (_b0 == null || !Brushable) return;
        var p = e.GetPosition(this);
        _b1 = new Point(Math.Clamp(p.X, _frame.Left, _frame.Right), Math.Clamp(p.Y, _frame.Top, _frame.Bottom));
        InvalidateVisual();
    }
    protected override void OnPointerReleased(Avalonia.Input.PointerReleasedEventArgs e)
    {
        base.OnPointerReleased(e);
        if (_b0 is not { } a || _b1 is not { } b || !Brushable || _xinv == null || _yinv == null) return;
        e.Pointer.Capture(null);
        if (Math.Abs(a.X - b.X) < 3 && Math.Abs(a.Y - b.Y) < 3) { _b0 = _b1 = null; InvalidateVisual(); Brushed?.Invoke(double.NaN, double.NaN, double.NaN, double.NaN); return; }
        double x0 = _xinv(a.X), x1 = _xinv(b.X), y0 = _yinv(a.Y), y1 = _yinv(b.Y);
        Brushed?.Invoke(Math.Min(x0, x1), Math.Max(x0, x1), Math.Min(y0, y1), Math.Max(y0, y1));
    }
    /// <summary>Clears the drawn brush box (new data).</summary>
    public void ClearBrush() { _b0 = _b1 = null; InvalidateVisual(); }

    // ------------------------------------------------------------------ colours
    private static bool Light => Application.Current?.ActualThemeVariant == ThemeVariant.Light;
    // Okabe & Ito's colour-blind safe set, lightened on the dark theme so every colour reads on the background
    private static readonly string[] DarkPalette = ["#5AB4E5", "#F0A030", "#3CC49A", "#E07CB4", "#E8DA5A", "#F07050"];
    private static readonly string[] LightPalette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#B8860B", "#E69F00"];
    private static IBrush Ink => Tokens.Brush("TextB");
    private static IBrush Soft => Tokens.Brush("MutedB");
    internal static IBrush Colour(int k) => k < 0 ? Ink : new SolidColorBrush(Color.Parse((Light ? LightPalette : DarkPalette)[k % DarkPalette.Length]));
    private static Color ColourC(int k) => k < 0 ? ((Ink as ISolidColorBrush)?.Color ?? Colors.White) : Color.Parse((Light ? LightPalette : DarkPalette)[k % DarkPalette.Length]);
    // viridis at five stops (van der Walt & Smith, matplotlib 2015)
    private static readonly Color[] Viridis = [Color.FromRgb(68, 1, 84), Color.FromRgb(59, 82, 139), Color.FromRgb(33, 145, 140), Color.FromRgb(94, 201, 98), Color.FromRgb(253, 231, 37)];
    private static Color Ramp(double t)
    {
        t = Math.Clamp(double.IsFinite(t) ? t : 0, 0, 1) * (Viridis.Length - 1);
        var k = Math.Min(Viridis.Length - 2, (int)t);
        var f = t - k;
        Color a = Viridis[k], c = Viridis[k + 1];
        return Color.FromRgb((byte)(a.R + (c.R - a.R) * f), (byte)(a.G + (c.G - a.G) * f), (byte)(a.B + (c.B - a.B) * f));
    }

    // ------------------------------------------------------------------ axes
    private sealed class Axis
    {
        public double Lo, Hi;
        public bool Log;
        public double Step;            // major step (decades on a log axis)
        public string[]? Names;        // category names at 0, 1, 2, …
        public List<double> Major = [], Minor = [];
        public double Map(double v, double a, double b) => Log ? a + (Math.Log10(v) - Lo) / (Hi - Lo) * (b - a) : a + (v - Lo) / (Hi - Lo) * (b - a);
        public double Unmap(double p, double a, double b) { var t = Lo + (p - a) / (b - a) * (Hi - Lo); return Log ? Math.Pow(10, t) : t; }
        /// <summary>A major tick's position in data units (log axes keep their ticks as exponents).</summary>
        public double Value(double t) => Log ? Math.Pow(10, t) : t;
        public string Label(double v)
        {
            if (Names != null) { var i = (int)Math.Round(v); return i >= 0 && i < Names.Length ? Names[i] : ""; }
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

    /// <summary>True when every value is a whole multiple of one step that is 1, 2, 2.5 or 5 times a power of ten: the
    /// tick values a reader expects (0, 0.5, 1.0 …; 0, 250, 500 …).</summary>
    public static bool RoundTicks(IReadOnlyList<double> v)
    {
        if (v.Count < 2) return v.Count == 1;
        var step = v[1] - v[0];
        if (!(step > 0)) return false;
        var mag = Math.Pow(10, Math.Floor(Math.Log10(step)));
        var f = step / mag;
        if (!new[] { 1, 2, 2.5, 5, 10 }.Any(n => Math.Abs(f - n) < 1e-6 * n)) return false;
        for (var i = 0; i < v.Count; ++i)
        {
            if (i > 0 && Math.Abs(v[i] - v[i - 1] - step) > 1e-6 * step) return false;
            var k = v[i] / step;
            if (Math.Abs(k - Math.Round(k)) > 1e-6) return false;
        }
        return true;
    }

    /// <summary>A linear axis over [lo, hi] with about `target` major ticks; the ends snap to the minor step when that wastes
    /// little room, else stay at a small padding past the data.</summary>
    private static Axis Linear(double lo, double hi, int target, double padLo, double padHi, double? fixLo, double? fixHi, double minStep = 0)
    {
        if (fixLo is double fl0 && fixHi is double fh0 && fh0 > fl0) { lo = fl0; hi = fh0; }
        // a flat series (all zero: a count not started yet) still gets a readable axis
        if (!(hi > lo)) { var c = double.IsFinite(lo) ? lo : 0; var d = c == 0 ? 1 : Math.Abs(c) * 0.05; lo = c - d; hi = c + d; }
        if (fixLo is double fl && fl < hi) lo = fl;
        if (fixHi is double fh && fh > lo) hi = fh;
        var span = hi - lo;
        if (lo != 0 && fixLo == null) lo -= span * padLo;   // an axis that starts at zero stays there
        if (fixHi == null) hi += span * padHi;
        var step = NiceStep(hi - lo, target);
        if (minStep > 0) step = step < minStep ? minStep : Math.Floor(step / minStep) * minStep;   // whole counts: 2 rather than 2.5
        var minor = step / (Math.Abs(step / Math.Pow(10, Math.Floor(Math.Log10(step))) - 2) < 1e-9 ? 4 : 5);
        var sLo = Math.Floor(lo / minor + 1e-9) * minor;
        var sHi = Math.Ceiling(hi / minor - 1e-9) * minor;
        var ax = new Axis { Lo = fixLo ?? sLo, Hi = fixHi ?? sHi, Step = step };
        if (!(ax.Hi > ax.Lo)) ax.Hi = ax.Lo + step;
        var k0 = Math.Ceiling(ax.Lo / step - 1e-9);
        for (var k = k0; k * step <= ax.Hi + step * 1e-9; k += 1) { var t = k * step; ax.Major.Add(Math.Abs(t) < step * 1e-9 ? 0 : t); }
        var m0 = Math.Ceiling(ax.Lo / minor - 1e-9);
        for (var k = m0; k * minor <= ax.Hi + minor * 1e-9; k += 1) ax.Minor.Add(k * minor);
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

    private static Axis Categorical(int n)
    {
        var ax = new Axis { Lo = -0.5, Hi = Math.Max(1, n) - 0.5, Step = 1 };
        for (var i = 0; i < n; ++i) ax.Major.Add(i);
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

    /// <summary>The text a title shows, without its markup (for checks and tooltips).</summary>
    public static string Plain(string s) => string.Concat(Parse(s).Select(r => r.Text));

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

    [ThreadStatic] private static bool _estimate;   // no font system (the headless self-test): widths from the character count
    private static double RichWidth(string s, double size) => _estimate ? 0.56 * size * Plain(s).Length : Layout(s, size, Brushes.White).W;

    // ------------------------------------------------------------------ layout: axes and margins for a size
    private sealed class Frame
    {
        public Axis X = null!, Y = null!;
        public Axis? Y2, Z;
        public double L, R, T, B, W, H;
        public double FTick, FLabel, FLegend, TickLen, Gap;
        public bool Compact;
        public double PX(double v) => X.Map(v, L, L + W);
        public double PY(double v) => Y.Map(v, T + H, T);
        public double PY2(double v) => Y2!.Map(v, T + H, T);
    }

    private static bool Ok(SciFigure fig, double x, double y, bool right) => double.IsFinite(x) && double.IsFinite(y) && (!fig.XLog || x > 0) && (right || !fig.YLog || y > 0);

    /// <summary>The axes and margins of a figure drawn at this size (null: nothing to draw).</summary>
    private static Frame? Compute(SciFigure fig, Size size, bool compact)
    {
        var small = !compact && size.Height < 175;
        var tiny = !compact && size.Height < 140;   // a monitor squeezed into a third of a column
        var f = new Frame
        {
            Compact = compact,
            FTick = compact || tiny ? 9.5 : small ? 10.5 : 11.5,
            FLabel = compact || tiny ? 10.5 : small ? 11.5 : 13,
            FLegend = compact || tiny ? 9 : small ? 10 : 11,
            TickLen = compact || tiny ? 4.0 : small ? 5.0 : 6.0,
            Gap = compact || tiny ? 3.0 : small ? 4.0 : 5.0,
        };
        if (fig.Heat is { } hm) return ComputeHeat(fig, hm, size, f);
        var series = fig.Series.Where(s => s.X.Length > 0).ToArray();
        var left = new List<(double x, double y, double e)>();
        foreach (var s in series.Where(s => !s.Right))
            for (var i = 0; i < Math.Min(s.X.Length, s.Y.Length); ++i)
                if (Ok(fig, s.X[i], s.Y[i], false)) left.Add((s.X[i], s.Y[i], s.Err is { } er && i < er.Length && double.IsFinite(er[i]) ? er[i] : 0));
        foreach (var p in fig.Points) if (Ok(fig, p.X, p.Y, false)) left.Add((p.X, p.Y, 0));
        var bars = series.Where(s => s.Kind == SciKind.Bars).ToArray();
        if (left.Count < 2 && !(bars.Length > 0 && left.Count >= 1) && fig.Categories.Length == 0) return null;
        if (left.Count == 0) left.Add((0, 0, 0));
        var rightPts = series.Where(s => s.Right).SelectMany(s => s.X.Zip(s.Y)).Where(p => Ok(fig, p.First, p.Second, true)).ToArray();
        var hasRight = rightPts.Length >= 2;

        double xlo = left.Min(p => p.x), xhi = left.Max(p => p.x);
        if (hasRight) { xlo = Math.Min(xlo, rightPts.Min(p => p.First)); xhi = Math.Max(xhi, rightPts.Max(p => p.First)); }
        foreach (var st in series.Where(s => s.Kind is SciKind.Steps or SciKind.Bars && s.X.Length > 1 && !fig.XLog))   // a histogram's outer bins drawn whole
        {
            var half = MinSpacing(st.X) / 2;
            xlo = Math.Min(xlo, st.X.Min() - half);
            xhi = Math.Max(xhi, st.X.Max() + half);
        }
        if (bars.Length > 0 && bars.All(s => s.X.Length == 1) && fig.Categories.Length == 0) { xlo -= 0.5; xhi += 0.5; }
        double ylo = left.Min(p => p.y - p.e), yhi = left.Max(p => p.y + p.e);
        if (fig.RefY is double ry && !fig.YLog)   // a reference line above the data keeps some room over it, off the frame
        {
            if (ry >= yhi) yhi = ry + 0.12 * (ry - Math.Min(ylo, ry));
            ylo = Math.Min(ylo, ry); yhi = Math.Max(yhi, ry);
        }
        if (fig.YBand is { } yb) { ylo = Math.Min(ylo, yb.Lo); yhi = Math.Max(yhi, yb.Hi); }
        foreach (var g in fig.Guides)
            if (double.IsFinite(g.Y) && (!fig.YLog || g.Y > 0))
            {
                if (!fig.YLog && g.Y >= yhi) yhi = g.Y + 0.08 * (g.Y - Math.Min(ylo, g.Y));   // a target above the data stays off the frame
                ylo = Math.Min(ylo, g.Y); yhi = Math.Max(yhi, g.Y);
            }
        if (bars.Length > 0 && !fig.YLog) { ylo = Math.Min(ylo, 0); yhi = Math.Max(yhi, 0); }
        if (fig.YZero && !fig.YLog && ylo >= 0 && ylo < 0.3 * (yhi - ylo)) ylo = 0;
        var nX = compact ? 4 : Math.Clamp((int)(size.Width / 110), 3, 8);
        var nY = compact ? 4 : Math.Clamp((int)(size.Height / 60), 3, 7);
        var xpad = series.Any(s => s.Kind is SciKind.Markers or SciKind.LineMarkers) || fig.Points.Length > 0 ? 0.03 : 0;
        f.X = fig.Categories.Length > 0 ? Categorical(fig.Categories.Length)
            : fig.XLog ? Logarithmic(xlo, xhi, fig.XMin, fig.XMax) : Linear(xlo, xhi, nX, xpad, xpad, fig.XMin, fig.XMax, fig.XWhole ? 1 : 0);
        if (fig.Categories.Length > 0) f.X.Names = fig.Categories;
        var yFromZero = fig.YMin == null && fig.YZero && ylo == 0 || bars.Length > 0 && ylo == 0;
        f.Y = fig.YLog ? Logarithmic(ylo, yhi, fig.YMin, fig.YMax) : Linear(ylo, yhi, nY, yFromZero ? 0 : 0.04, 0.04, fig.YMin ?? (yFromZero ? 0 : null), fig.YMax, fig.YWhole ? 1 : 0);
        f.Y2 = hasRight ? Linear(fig.Y2Min ?? rightPts.Min(p => p.Second), rightPts.Max(p => p.Second), nY, 0.04, 0.04, fig.Y2Min, null) : null;

        var yTickW = f.Y.Major.Select(v => RichWidth(f.Y.Label(v), f.FTick)).DefaultIfEmpty(0).Max();
        var y2TickW = f.Y2?.Major.Select(v => RichWidth(f.Y2.Label(v), f.FTick)).DefaultIfEmpty(0).Max() ?? 0;
        f.L = 4 + f.FLabel * 1.3 + f.Gap + yTickW + f.Gap;
        f.R = hasRight ? 4 + f.FLabel * 1.3 + f.Gap + y2TickW + f.Gap : compact ? 8 : 14;
        f.T = compact ? 6.0 : 10.0;
        f.B = f.Gap + f.FTick * 1.3 + f.Gap + f.FLabel * 1.35 + 2;
        f.W = size.Width - f.L - f.R;
        f.H = size.Height - f.T - f.B;
        return f.W < 40 || f.H < 30 ? null : f;
    }

    private static Frame? ComputeHeat(SciFigure fig, SciHeat hm, Size size, Frame f)
    {
        var xs = hm.X.Where(double.IsFinite).Distinct().OrderBy(v => v).ToArray();
        var ys = hm.Y.Where(double.IsFinite).Distinct().OrderBy(v => v).ToArray();
        if (xs.Length < 1 || ys.Length < 1) return null;
        double dx = xs.Length > 1 ? xs[1] - xs[0] : 1, dy = ys.Length > 1 ? ys[1] - ys[0] : 1;
        double x0 = xs[0] - dx / 2, x1 = xs[^1] + dx / 2, y0 = ys[0] - dy / 2, y1 = ys[^1] + dy / 2;
        var nX = Math.Clamp((int)(size.Width / 120), 3, 7);
        var nY = Math.Clamp((int)(size.Height / 60), 3, 7);
        f.X = Linear(x0, x1, nX, 0, 0, x0, x1);
        f.Y = Linear(y0, y1, nY, 0, 0, y0, y1);
        var fin = hm.Z.Where(double.IsFinite).ToArray();
        double zmin = fin.Length > 0 ? fin.Min() : 0, zmax = fin.Length > 0 ? fin.Max() : 1;
        if (zmin >= 0 && zmin < 0.3 * (zmax - zmin)) zmin = 0;
        if (zmax - zmin < 1e-12) zmax = zmin + 1;
        f.Z = Linear(zmin, zmax, nY, 0, 0, zmin, zmax);
        var yTickW = f.Y.Major.Select(v => RichWidth(f.Y.Label(v), f.FTick)).DefaultIfEmpty(0).Max();
        var zTickW = f.Z.Major.Select(v => RichWidth(f.Z.Label(v), f.FTick)).DefaultIfEmpty(0).Max();
        f.L = 4 + f.FLabel * 1.3 + f.Gap + yTickW + f.Gap;
        f.R = 14 + 11 + f.Gap + zTickW + f.Gap + f.FLabel * 1.3 + 4;
        f.T = 10;
        f.B = f.Gap + f.FTick * 1.3 + f.Gap + f.FLabel * 1.35 + 2;
        f.W = size.Width - f.L - f.R;
        f.H = size.Height - f.T - f.B;
        return f.W < 40 || f.H < 30 ? null : f;
    }

    private static double MinSpacing(double[] x)
    {
        var s = x.Where(double.IsFinite).Distinct().OrderBy(v => v).ToArray();
        var d = double.PositiveInfinity;
        for (var i = 1; i < s.Length; ++i) d = Math.Min(d, s[i] - s[i - 1]);
        return double.IsFinite(d) && d > 0 ? d : 1;
    }

    /// <summary>The major tick values of a figure drawn at this size, in data units (a log axis: its decades), and the
    /// tick labels as drawn. The self-test reads these to check every figure's ticks are round.</summary>
    public static (double[] X, double[] Y, string[] XText, string[] YText)? Ticks(SciFigure fig, double width, double height, bool compact = false)
    {
        Frame? f;
        _estimate = Application.Current == null;   // the headless self-test runs without a font system
        try { f = Compute(fig, new Size(width, height), compact); }
        finally { _estimate = false; }
        if (f == null) return null;
        return (f.X.Major.ToArray(), f.Y.Major.ToArray(), f.X.Major.Select(f.X.Label).ToArray(), f.Y.Major.Select(f.Y.Label).ToArray());
    }

    // ------------------------------------------------------------------ render
    public override void Render(DrawingContext ctx)
    {
        var fig = Figure;
        var b = Bounds;
        var f = fig == null ? null : Compute(fig, b.Size, Compact);
        if (fig == null || f == null)
        {
            var text = fig?.Empty ?? EmptyText;
            _xinv = _yinv = null;
            if (text.Length == 0 || b.Width < 20 || b.Height < 14) return;
            var ft = new FormattedText(text, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, new Typeface(Tokens.Sans), Compact ? 9.5 : 11.5, Soft);
            ctx.DrawText(ft, new Point((b.Width - ft.Width) / 2, (b.Height - ft.Height) / 2));
            return;
        }
        var compact = f.Compact;
        double L = f.L, T = f.T, w = f.W, h = f.H;
        var xa = f.X;
        var ya = f.Y;
        _frame = new Rect(L, T, w, h);
        _xinv = px => xa.Unmap(px, L, L + w);
        _yinv = py => ya.Unmap(py, T + h, T);
        var ink = Ink;
        var tickPen = new Pen(ink, compact ? 0.8 : 1.0);
        var bg = Tokens.Brush("Bg1B");
        if (fig.Heat is { } hm) { RenderHeat(ctx, fig, hm, f, tickPen); DrawAxes(ctx, fig, f, tickPen, false); DrawBrush(ctx); return; }

        var series = fig.Series.Where(s => s.X.Length > 0).ToArray();
        // shaded stages, x window and horizontal band under everything
        var keepOut = new List<Rect>();
        foreach (var st in fig.Stages)
        {
            var x0 = Math.Clamp(f.PX(st.From), L, L + w);
            var x1 = Math.Clamp(f.PX(st.To), L, L + w);
            if (x1 <= x0) continue;
            if (st.Shade) ctx.FillRectangle(new SolidColorBrush(ColourC(1), Light ? 0.08 : 0.10), new Rect(x0, T, x1 - x0, h));
            var prev = Array.IndexOf(fig.Stages, st) is var si and > 0 ? fig.Stages[si - 1] : null;
            if (x0 > L + 1 && prev != null && prev.Shade != st.Shade && x1 - x0 >= 12 && x0 - f.PX(prev.From) >= 12) ctx.DrawLine(new Pen(Soft, 0.8, new DashStyle([2, 3], 0)), new Point(x0, T), new Point(x0, T + h));
            if (st.Label.Length > 0)
            {
                var lw = RichWidth(st.Label, f.FLegend * 0.92);
                var lx = (x0 + x1 - lw) / 2;
                // a label that would crowd the one before it (two short stages side by side) is left out
                if (lw + 6 < x1 - x0 && !keepOut.Any(k => lx < k.Right + 8 && lx + lw > k.Left - 8))
                {
                    DrawRich(ctx, st.Label, f.FLegend * 0.92, Soft, new Point((x0 + x1) / 2, T + f.TickLen + 2), 0.5, 0);
                    keepOut.Add(new Rect((x0 + x1 - lw) / 2, T, lw, f.TickLen + 2 + f.FLegend * 1.2));
                }
            }
        }
        if (fig.XBand is { } xb)
        {
            var x0 = Math.Clamp(f.PX(xb.From), L, L + w);
            var x1 = Math.Clamp(f.PX(xb.To), L, L + w);
            if (x1 > x0) ctx.FillRectangle(new SolidColorBrush(ColourC(0), Light ? 0.10 : 0.13), new Rect(x0, T, x1 - x0, h));
        }
        if (fig.YBand is { } band)
        {
            var y0 = Math.Clamp(f.PY(band.Hi), T, T + h);
            var y1 = Math.Clamp(f.PY(band.Lo), T, T + h);
            ctx.FillRectangle(new SolidColorBrush(ColourC(1), Light ? 0.16 : 0.22), new Rect(L, y0, w, Math.Max(1, y1 - y0)));
            ctx.DrawLine(new Pen(Colour(1), 1.2), new Point(L, (y0 + y1) / 2), new Point(L + w, (y0 + y1) / 2));
        }
        if (fig.RefY is double rv && (!fig.YLog || rv > 0))
        {
            var py = f.PY(rv);
            if (py >= T && py <= T + h) ctx.DrawLine(new Pen(Soft, 0.9, new DashStyle([4, 3], 0)), new Point(L, py), new Point(L + w, py));
        }
        foreach (var g in fig.Guides)
        {
            if (!double.IsFinite(g.Y) || (fig.YLog && g.Y <= 0)) continue;
            var py = f.PY(g.Y);
            if (py < T - 0.5 || py > T + h + 0.5) continue;
            ctx.DrawLine(new Pen(ink, compact ? 0.8 : 1.0, new DashStyle([5, 3], 0)), new Point(L, py), new Point(L + w, py));
        }

        // data, clipped to the frame
        var drawn = new List<(SciSeries S, List<Point> Pts)>();
        var barSeries = series.Where(s => s.Kind == SciKind.Bars).ToList();
        var barW = barSeries.Count > 0 ? (fig.Categories.Length > 0 ? 1 : MinSpacing(barSeries.SelectMany(s => s.X).ToArray())) : 0;
        using (ctx.PushClip(new Rect(L, T, w, h)))
        {
            foreach (var s in series)
            {
                var brush = Colour(s.Colour);
                var pts = new List<Point>();
                var errs = new List<(Point Lo, Point Hi)>();
                if (s.Kind == SciKind.Bars)
                {
                    var k = barSeries.IndexOf(s);
                    var gw = 0.78 * barW / barSeries.Count;
                    var y0 = f.PY(ya.Log ? Math.Pow(10, ya.Lo) : Math.Clamp(0, ya.Lo, ya.Hi));
                    for (var i = 0; i < Math.Min(s.X.Length, s.Y.Length); ++i)
                    {
                        if (!Ok(fig, s.X[i], s.Y[i], false)) continue;
                        var xl = f.PX(s.X[i] - 0.39 * barW + k * gw) + 0.5;
                        var xr = f.PX(s.X[i] - 0.39 * barW + (k + 1) * gw) - 0.5;
                        var py = f.PY(s.Y[i]);
                        var r = new Rect(new Point(xl, Math.Min(py, y0)), new Point(Math.Max(xl + 1, xr), Math.Max(py, y0)));
                        ctx.FillRectangle(new SolidColorBrush(ColourC(s.Colour), Light ? 0.55 : 0.6), r);
                        ctx.DrawRectangle(null, new Pen(brush, compact ? 0.8 : 1.0), r);
                        pts.Add(new Point((xl + xr) / 2, py));
                        pts.Add(new Point((xl + xr) / 2, (py + y0) / 2));
                    }
                    drawn.Add((s with { Kind = SciKind.Markers }, pts));
                    continue;
                }
                var breaks = new HashSet<int>();   // a gap in the data (NaN) starts a new piece of line
                for (var i = 0; i < Math.Min(s.X.Length, s.Y.Length); ++i)
                {
                    if (!Ok(fig, s.X[i], s.Y[i], s.Right)) { if (pts.Count > 0) breaks.Add(pts.Count); continue; }
                    var px = f.PX(s.X[i]);
                    var py = s.Right ? f.PY2(s.Y[i]) : f.PY(s.Y[i]);
                    pts.Add(new Point(px, py));
                    if (s.Err is { } er && i < er.Length && er[i] > 0 && double.IsFinite(er[i]) && !s.Right)
                        errs.Add((new Point(px, f.PY(fig.YLog ? Math.Max(s.Y[i] - er[i], Math.Pow(10, ya.Lo)) : s.Y[i] - er[i])), new Point(px, f.PY(s.Y[i] + er[i]))));
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
                var line = s.Kind is SciKind.Line or SciKind.Dashed or SciKind.Dotted && breaks.Count == 0 ? Decimate(pts, w) : pts;
                drawn.Add((s, line));
                if (pen != null && line.Count > 1)
                {
                    var geo = new StreamGeometry();
                    using (var g = geo.Open())
                    {
                        if (s.Kind == SciKind.Steps)
                        {
                            var half = (pts[1].X - pts[0].X) / 2;
                            var baseY = f.PY(ya.Log ? Math.Pow(10, ya.Lo) : Math.Max(0, ya.Lo));
                            g.BeginFigure(new Point(pts[0].X - half, baseY), s.Fill);
                            foreach (var p in pts) { g.LineTo(new Point(p.X - half, p.Y)); g.LineTo(new Point(p.X + half, p.Y)); }
                            g.LineTo(new Point(pts[^1].X + half, baseY));
                            g.EndFigure(s.Fill);
                        }
                        else
                        {
                            g.BeginFigure(line[0], false);
                            for (var i = 1; i < line.Count; ++i)
                            {
                                if (breaks.Contains(i)) { g.EndFigure(false); g.BeginFigure(line[i], false); }
                                else g.LineTo(line[i]);
                            }
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
            foreach (var p in fig.Points)
            {
                if (!Ok(fig, p.X, p.Y, false)) continue;
                var pt = new Point(f.PX(p.X), f.PY(p.Y));
                Marker(ctx, pt, compact ? 3 : 4.2, Colour(p.Colour), false, false, bg, compact);
                if (p.Label.Length > 0)
                {
                    var lw = RichWidth(p.Label, f.FLegend);
                    var right = pt.X + 7 + lw < L + w - 2;
                    var up = pt.Y - f.FLegend * 1.4 > T + 2;
                    var at = new Point(right ? pt.X + 7 : pt.X - 7, up ? pt.Y - 4 : pt.Y + 4);
                    DrawRich(ctx, p.Label, f.FLegend, ink, at, right ? 0 : 1, up ? 1 : 0);
                    keepOut.Add(new Rect(right ? at.X : at.X - lw, up ? at.Y - f.FLegend * 1.25 : at.Y, lw, f.FLegend * 1.25));
                }
                drawn.Add((new SciSeries("", [], []) { Kind = SciKind.Markers }, [pt]));
            }
        }
        // vertical marks
        foreach (var m in fig.Marks)
        {
            if (!double.IsFinite(m.X) || (fig.XLog && m.X <= 0)) continue;
            var px = f.PX(m.X);
            if (px < L || px > L + w) continue;
            ctx.DrawLine(new Pen(ink, compact ? 0.8 : 1.0, new DashStyle([3, 3], 0)), new Point(px, T), new Point(px, T + h));
            if (m.Label.Length > 0)
            {
                var lw = RichWidth(m.Label, f.FLegend);
                var right = px + 4 + lw < L + w - 2;
                var top = T + f.TickLen + 1 + keepOut.Where(r => r.Top <= T + 1 && r.Left < (right ? px + 4 + lw : px) && r.Right > (right ? px : px - 4 - lw)).Select(r => r.Height).DefaultIfEmpty(0).Max();
                DrawRich(ctx, m.Label, f.FLegend, ink, new Point(right ? px + 4 : px - 4, top), right ? 0 : 1, 0);
                keepOut.Add(new Rect(right ? px + 4 : px - 4 - lw, top, lw, f.FLegend * 1.25));
            }
        }

        var hasRight = f.Y2 != null;
        DrawAxes(ctx, fig, f, tickPen, hasRight);
        DrawBrush(ctx);
        if (compact) return;   // a card's figure is too small for a legend: the card says what it shows

        // legend in the corner with the fewest data points under it
        var entries = series.Where(s => s.Legend && s.Name.Length > 0).ToArray();
        var notes = fig.Notes.Where(n => n.Length > 0).ToArray();
        if (fig.YBand is { Label.Length: > 0 } ybl) entries = [.. entries, new SciSeries(ybl.Label, [], []) { Colour = 1, Kind = SciKind.Line, Width = 6 }];
        if (fig.XBand is { Label.Length: > 0 } xbl) entries = [.. entries, new SciSeries(xbl.Label, [], []) { Colour = 0, Kind = SciKind.Line, Width = 6 }];
        foreach (var g in fig.Guides.Where(g => g.Label.Length > 0 && double.IsFinite(g.Y))) entries = [.. entries, new SciSeries(g.Label, [], []) { Colour = -1, Kind = SciKind.Dashed, Width = 1 }];   // a guide is named in the legend, not on its line
        if (entries.Length == 0 && notes.Length == 0) return;
        var sample = 22.0;
        var rowH = f.FLegend * 1.45;
        var pad = 6.0;
        var lw2 = Math.Max(entries.Select(e => sample + 5 + RichWidth(e.Name, f.FLegend)).DefaultIfEmpty(0).Max(), notes.Select(n => RichWidth(n, f.FLegend)).DefaultIfEmpty(0).Max());
        var boxW = Math.Min(w - 8, lw2 + 2 * pad);
        var boxH = (entries.Length + notes.Length) * rowH + 2 * pad - (rowH - f.FLegend * 1.2);
        if (boxH > h - 8) return;
        var inset = f.TickLen + 3;
        var corners = new[] { new Point(L + w - boxW - inset, T + inset), new Point(L + inset, T + inset), new Point(L + w - boxW - inset, T + h - boxH - inset), new Point(L + inset, T + h - boxH - inset) };
        // the points and the lines between them, every few pixels (a legend should not sit on a line either)
        var allPts = new List<Point>();
        foreach (var (s, pts) in drawn)
        {
            for (var i = 0; i < pts.Count; ++i)
            {
                allPts.Add(pts[i]);
                if (i == 0 || s.Kind == SciKind.Markers) continue;
                var d = pts[i] - pts[i - 1];
                var n = (int)Math.Min(40, Math.Sqrt(d.X * d.X + d.Y * d.Y) / 5);
                for (var k = 1; k < n; ++k) allPts.Add(pts[i - 1] + d * (k / (double)n));
            }
        }
        var markLines = fig.Marks.Where(m => double.IsFinite(m.X) && (!fig.XLog || m.X > 0)).Select(m => f.PX(m.X)).ToArray();
        int Cost(Point c)
        {
            var r = new Rect(c, new Size(boxW, boxH)).Inflate(3);
            return allPts.Count(r.Contains) + markLines.Count(x => x >= r.Left && x <= r.Right) * 25 + keepOut.Count(k => k.Intersects(r)) * 200;
        }
        var at0 = corners.OrderBy(Cost).First();
        var box = new Rect(at0, new Size(boxW, boxH));
        ctx.FillRectangle(new SolidColorBrush((bg as ISolidColorBrush)?.Color ?? Colors.Black, 0.88), box);
        ctx.DrawRectangle(null, new Pen(Tokens.Brush("LineB"), 0.8), box);
        var y = box.Top + pad + f.FLegend * 0.6;
        foreach (var e in entries)
        {
            var x0 = box.Left + pad;
            var brush = Colour(e.Colour);
            if (e.Width >= 5 || e.Kind == SciKind.Bars)
            {
                var r = new Rect(x0 + (e.Kind == SciKind.Bars ? 5 : 0), y - f.FLegend * 0.38, e.Kind == SciKind.Bars ? sample - 10 : sample, f.FLegend * 0.76);
                ctx.FillRectangle(new SolidColorBrush(ColourC(e.Colour), e.Kind == SciKind.Bars ? (Light ? 0.55 : 0.6) : Light ? 0.25 : 0.3), r);
                if (e.Kind == SciKind.Bars) ctx.DrawRectangle(null, new Pen(brush, 1), r);
            }
            else if (e.Kind is SciKind.Line or SciKind.LineMarkers or SciKind.Steps or SciKind.Dashed or SciKind.Dotted)
            {
                IPen pen = e.Kind == SciKind.Dashed ? new Pen(brush, Math.Min(2, e.Width), new DashStyle([5, 3], 0))
                    : e.Kind == SciKind.Dotted ? new Pen(brush, Math.Min(2, e.Width), new DashStyle([1, 2.5], 0), lineCap: PenLineCap.Round)
                    : new Pen(brush, Math.Min(2, e.Width));
                ctx.DrawLine(pen, new Point(x0, y), new Point(x0 + sample, y));
            }
            if (e.Kind is SciKind.Markers or SciKind.LineMarkers) Marker(ctx, new Point(x0 + sample / 2, y), 3.2, brush, e.Open, e.Square, bg, false);
            DrawRich(ctx, e.Name, f.FLegend, ink, new Point(x0 + sample + 5, y), 0, 0.5);
            y += rowH;
        }
        foreach (var n in notes) { DrawRich(ctx, n, f.FLegend, ink, new Point(box.Left + pad, y), 0, 0.5); y += rowH; }
    }

    /// <summary>A long polyline kept to its first, lowest, highest and last point in each pixel column: the drawn line is
    /// the same, the geometry stays a few thousand points.</summary>
    private static List<Point> Decimate(List<Point> pts, double width)
    {
        if (pts.Count < 4 * Math.Max(50, width)) return pts;
        for (var i = 1; i < pts.Count; ++i) if (pts[i].X < pts[i - 1].X) return pts;   // not along x: keep it whole
        var outp = new List<Point>((int)(4 * width) + 4);
        var i0 = 0;
        while (i0 < pts.Count)
        {
            var col = Math.Floor(pts[i0].X);
            int lo = i0, hi = i0, i1 = i0;
            while (i1 + 1 < pts.Count && Math.Floor(pts[i1 + 1].X) == col)
            {
                i1++;
                if (pts[i1].Y < pts[lo].Y) lo = i1;
                if (pts[i1].Y > pts[hi].Y) hi = i1;
            }
            outp.Add(pts[i0]);
            if (lo < hi) { if (lo != i0) outp.Add(pts[lo]); if (hi != i1) outp.Add(pts[hi]); }
            else { if (hi != i0 && hi != lo) outp.Add(pts[hi]); if (lo != i1 && lo != i0) outp.Add(pts[lo]); }
            if (i1 != i0) outp.Add(pts[i1]);
            i0 = i1 + 1;
        }
        return outp;
    }

    private void DrawBrush(DrawingContext ctx)
    {
        if (_b0 is not { } ba || _b1 is not { } bb) return;
        var r = new Rect(new Point(Math.Min(ba.X, bb.X), Math.Min(ba.Y, bb.Y)), new Point(Math.Max(ba.X, bb.X), Math.Max(ba.Y, bb.Y)));
        ctx.FillRectangle(new SolidColorBrush(ColourC(0), 0.14), r);
        ctx.DrawRectangle(null, new Pen(Colour(0), 1, new DashStyle([4, 3], 0)), r);
    }

    /// <summary>The frame, the inward ticks on all four sides (the right side carries the second axis when there is one), the
    /// tick labels and the axis titles.</summary>
    private static void DrawAxes(DrawingContext ctx, SciFigure fig, Frame f, IPen tickPen, bool hasRight)
    {
        double L = f.L, T = f.T, w = f.W, h = f.H, tickLen = f.TickLen, gap = f.Gap;
        var ink = Ink;
        var xa = f.X;
        var ya = f.Y;
        ctx.DrawRectangle(null, new Pen(ink, f.Compact ? 0.9 : 1.1), new Rect(L, T, w, h));
        foreach (var v in xa.Minor)
        {
            var px = f.PX(xa.Value(v));
            if (px < L - 0.5 || px > L + w + 0.5) continue;
            ctx.DrawLine(tickPen, new Point(px, T + h), new Point(px, T + h - tickLen * 0.55));
            ctx.DrawLine(tickPen, new Point(px, T), new Point(px, T + tickLen * 0.55));
        }
        // tick labels that would touch are thinned out (every second, third … label), never overlapped
        var xLabels = xa.Major.Select(v => (v, s: xa.Label(v))).ToArray();
        var widest = xLabels.Select(p => RichWidth(p.s, f.FTick)).DefaultIfEmpty(0).Max();
        var spacing = xLabels.Length > 1 ? Math.Abs(f.PX(xa.Value(xLabels[1].v)) - f.PX(xa.Value(xLabels[0].v))) : w;
        var every = Math.Max(1, (int)Math.Ceiling((widest + 6) / Math.Max(1, spacing)));
        for (var i = 0; i < xLabels.Length; ++i)
        {
            var px = f.PX(xa.Value(xLabels[i].v));
            if (px < L - 0.5 || px > L + w + 0.5) continue;
            ctx.DrawLine(tickPen, new Point(px, T + h), new Point(px, T + h - tickLen));
            ctx.DrawLine(tickPen, new Point(px, T), new Point(px, T + tickLen));
            if (i % every == 0) DrawRich(ctx, xLabels[i].s, f.FTick, ink, new Point(px, T + h + gap), 0.5, 0);
        }
        foreach (var v in ya.Minor)
        {
            var py = f.PY(ya.Value(v));
            if (py < T - 0.5 || py > T + h + 0.5) continue;
            ctx.DrawLine(tickPen, new Point(L, py), new Point(L + tickLen * 0.55, py));
            if (!hasRight && f.Z == null) ctx.DrawLine(tickPen, new Point(L + w, py), new Point(L + w - tickLen * 0.55, py));
        }
        foreach (var v in ya.Major)
        {
            var py = f.PY(ya.Value(v));
            if (py < T - 0.5 || py > T + h + 0.5) continue;
            ctx.DrawLine(tickPen, new Point(L, py), new Point(L + tickLen, py));
            if (!hasRight && f.Z == null) ctx.DrawLine(tickPen, new Point(L + w, py), new Point(L + w - tickLen, py));
            DrawRich(ctx, ya.Label(v), f.FTick, ink, new Point(L - gap, py), 1, 0.5);
        }
        if (f.Y2 is { } y2)
        {
            foreach (var v in y2.Minor) { var py = f.PY2(v); if (py >= T - 0.5 && py <= T + h + 0.5) ctx.DrawLine(tickPen, new Point(L + w, py), new Point(L + w - tickLen * 0.55, py)); }
            foreach (var v in y2.Major)
            {
                var py = f.PY2(v);
                if (py < T - 0.5 || py > T + h + 0.5) continue;
                ctx.DrawLine(tickPen, new Point(L + w, py), new Point(L + w - tickLen, py));
                DrawRich(ctx, y2.Label(v), f.FTick, ink, new Point(L + w + gap, py), 0, 0.5);
            }
        }
        DrawRich(ctx, fig.XLabel, f.FLabel, ink, new Point(L + w / 2, T + h + gap + f.FTick * 1.3 + gap), 0.5, 0);
        var ytw = RichWidth(fig.YLabel, f.FLabel);
        var yc = Math.Clamp(T + h / 2, Math.Min(2 + ytw / 2, T + h), Math.Max(T + h + f.B - 2 - ytw / 2, 2 + ytw / 2));
        using (ctx.PushTransform(Matrix.CreateRotation(-Math.PI / 2) * Matrix.CreateTranslation(4, yc)))
            DrawRich(ctx, fig.YLabel, f.FLabel, ink, new Point(0, 0), 0.5, 0);
        if (f.Y2 != null && fig.Y2Label.Length > 0)
            using (ctx.PushTransform(Matrix.CreateRotation(Math.PI / 2) * Matrix.CreateTranslation(L + w + f.R - 4, T + h / 2)))
                DrawRich(ctx, fig.Y2Label, f.FLabel, ink, new Point(0, 0), 0.5, 0);
    }

    /// <summary>The cells of a map and its colour bar (viridis) with its own ticks and title.</summary>
    private static void RenderHeat(DrawingContext ctx, SciFigure fig, SciHeat hm, Frame f, IPen tickPen)
    {
        var xs = hm.X.Where(double.IsFinite).Distinct().OrderBy(v => v).ToArray();
        var ys = hm.Y.Where(double.IsFinite).Distinct().OrderBy(v => v).ToArray();
        double dx = xs.Length > 1 ? xs[1] - xs[0] : 1, dy = ys.Length > 1 ? ys[1] - ys[0] : 1;
        var z = f.Z!;
        double zmin = z.Lo, zmax = z.Hi;
        using (ctx.PushClip(new Rect(f.L, f.T, f.W, f.H)))
            for (var k = 0; k < hm.Z.Length && k < hm.X.Length && k < hm.Y.Length; ++k)
            {
                var r = new Rect(new Point(f.PX(hm.X[k] - dx / 2), f.PY(hm.Y[k] + dy / 2)), new Point(f.PX(hm.X[k] + dx / 2) + 0.5, f.PY(hm.Y[k] - dy / 2) + 0.5));
                ctx.FillRectangle(new SolidColorBrush(Ramp((hm.Z[k] - zmin) / (zmax - zmin))), r);
            }
        var bx = f.L + f.W + 14;
        const double bw = 11;
        for (var k = 0; k < 64; ++k)
            ctx.FillRectangle(new SolidColorBrush(Ramp(k / 63.0)), new Rect(bx, f.T + f.H - (k + 1) * f.H / 64, bw, f.H / 64 + 0.5));
        ctx.DrawRectangle(null, new Pen(Ink, 1), new Rect(bx, f.T, bw, f.H));
        double PZ(double v) => z.Map(v, f.T + f.H, f.T);
        foreach (var v in z.Major)
        {
            var py = PZ(v);
            ctx.DrawLine(tickPen, new Point(bx + bw, py), new Point(bx + bw - 3, py));
            DrawRich(ctx, z.Label(v), f.FTick, Ink, new Point(bx + bw + f.Gap, py), 0, 0.5);
        }
        if (hm.ZLabel.Length > 0)
            using (ctx.PushTransform(Matrix.CreateRotation(Math.PI / 2) * Matrix.CreateTranslation(f.L + f.W + f.R - 4, f.T + f.H / 2)))
                DrawRich(ctx, hm.ZLabel, f.FLabel, Ink, new Point(0, 0), 0.5, 0);
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
