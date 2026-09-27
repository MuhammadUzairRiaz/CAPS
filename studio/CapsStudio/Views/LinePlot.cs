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
    /// <summary>What an empty plot says (default "no data").</summary>
    public string EmptyText { get; set; } = "no data";
    public bool AutoRange { get; set; }
    public double? RefY { get; set; } = 1.0;
    /// <summary>A dashed vertical line at this x (a trajectory player's current frame).</summary>
    public double? CursorX { get; set; }
    /// <summary>The curve's colour (default: the accent).</summary>
    public IBrush? LineBrush { get; set; }

    // ---- brushing (design/boards/HistScatter): drag a box on the plot to select what lies in it
    /// <summary>Allow a drag box; Brushed gets (x0, x1, y0, y1) in data units, NaN when cleared by a click.</summary>
    public bool Brushable { get; set; }
    public event Action<double, double, double, double>? Brushed;
    private Rect _plot;
    private double _xmin, _xmax, _ymin, _ymax;
    private Point? _b0, _b1;

    protected override void OnPointerPressed(Avalonia.Input.PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        if (!Brushable) return;
        var p = e.GetPosition(this);
        if (!_plot.Contains(p)) return;
        _b0 = _b1 = p;
        e.Pointer.Capture(this);
        InvalidateVisual();
    }
    protected override void OnPointerMoved(Avalonia.Input.PointerEventArgs e)
    {
        base.OnPointerMoved(e);
        if (_b0 == null || !Brushable) return;
        var p = e.GetPosition(this);
        _b1 = new Point(Math.Clamp(p.X, _plot.Left, _plot.Right), Math.Clamp(p.Y, _plot.Top, _plot.Bottom));
        InvalidateVisual();
    }
    protected override void OnPointerReleased(Avalonia.Input.PointerReleasedEventArgs e)
    {
        base.OnPointerReleased(e);
        if (_b0 is not { } a || _b1 is not { } b || !Brushable) return;
        e.Pointer.Capture(null);
        if (Math.Abs(a.X - b.X) < 3 && Math.Abs(a.Y - b.Y) < 3) { _b0 = _b1 = null; InvalidateVisual(); Brushed?.Invoke(double.NaN, double.NaN, double.NaN, double.NaN); return; }
        double Dx(double px) => _xmin + (px - _plot.Left) / _plot.Width * (_xmax - _xmin);
        double Dy(double py) => _ymin + (_plot.Bottom - py) / _plot.Height * (_ymax - _ymin);
        Brushed?.Invoke(Math.Min(Dx(a.X), Dx(b.X)), Math.Max(Dx(a.X), Dx(b.X)), Math.Min(Dy(a.Y), Dy(b.Y)), Math.Max(Dy(a.Y), Dy(b.Y)));
    }
    /// <summary>Clears the drawn brush box (new data).</summary>
    public void ClearBrush() { _b0 = _b1 = null; InvalidateVisual(); }

    public void SetData((double X, double Y)[] data)
    {
        _data = data;
        _overlay = [];
        _second = [];
        _fit = [];
        InvalidateVisual();
    }

    /// <summary>Two curves to compare: A solid (selection blue), B dashed (accent).</summary>
    public void SetCompare((double X, double Y)[] a, (double X, double Y)[] b)
    {
        _data = a.Length >= 2 ? a : b;
        _second = a.Length >= 2 ? b : [];
        _overlay = [];
        Markers = false;
        InvalidateVisual();
    }
    private (double X, double Y)[] _second = [];
    private IPen CurveA => new Pen(Tokens.Brush(AccentFirst ? "AccB" : "SelB"), 2, lineJoin: PenLineJoin.Round);
    private IPen CurveB => new Pen(Tokens.Brush(AccentFirst ? "SelB" : "AccB"), 2, new DashStyle([5, 3], 0), lineJoin: PenLineJoin.Round);
    /// <summary>SetCompare draws A solid in the accent and B dashed in selection blue (default: the other way round).</summary>
    public bool AccentFirst { get; set; }

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
    /// <summary>Error bars (± value) for the points drawn with Markers, one per point; they widen the y range.</summary>
    public double[]? Errors { get; set; }

    /// <summary>Draw the data as bars (a histogram: each point is a bin centre); the y axis starts at zero.</summary>
    public bool Bars { get; set; }

    /// <summary>A shaded x range (a fit window), drawn under the curves.</summary>
    public (double From, double To)? Band { get; set; }
    /// <summary>Stages along x (an equilibration protocol): from, to, a short label; shaded ones (constant pressure) tinted.</summary>
    public (double From, double To, string Label, bool Shade)[] Stages { get; set; } = [];
    private (double X, double Y)[] _fit = [];
    private (double X, double Y)[] _third = [];
    /// <summary>A third curve (measured data), drawn dotted in the text colour; it takes part in the range.</summary>
    public void SetThird((double X, double Y)[] data) { _third = data; InvalidateVisual(); }
    /// <summary>A curve (accent) with a dashed fit line (selection blue) over it.</summary>
    public void SetWithFit((double X, double Y)[] data, (double X, double Y)[] fit)
    {
        _data = data;
        _fit = fit;
        _overlay = [];
        _second = [];
        Markers = false;
        InvalidateVisual();
    }

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
            Text(EmptyText, L + w / 2, T + h / 2, centre: true);
            return;
        }
        double xmin = 0, xmax, ymin = 0, ymax;
        if (AutoRange)
        {
            var all = _data.Concat(_overlay).Concat(_second).Concat(_third).ToArray();   // a fit line may run past the data: it does not set the range
            if (Markers && Errors is { } er && er.Length == _data.Length)
                all = all.Concat(_data.Select((p, i) => (p.X, p.Y + (double.IsFinite(er[i]) ? er[i] : 0)))).Concat(_data.Select((p, i) => (p.X, p.Y - (double.IsFinite(er[i]) ? er[i] : 0)))).ToArray();
            xmin = all.Min(p => p.X);
            xmax = all.Max(p => p.X);
            if (xmax - xmin < 1e-9) xmax = xmin + 1;
            ymin = all.Min(p => p.Y);
            ymax = all.Max(p => p.Y);
            if (RefY is double r) { ymin = Math.Min(ymin, r); ymax = Math.Max(ymax, r); }
            if (Bars) { ymin = 0; var dx = _data.Length > 1 ? _data[1].X - _data[0].X : 1; xmin -= dx / 2; xmax += dx / 2; }
            var pad = Math.Max(Math.Max(1e-9, Math.Abs(ymax) * 1e-3), (ymax - ymin) * 0.08);   // a flat series still gets a readable axis
            if (!Bars) ymin -= pad;
            ymax += pad;
        }
        else
        {
            xmax = _data[^1].X + (_data[1].X - _data[0].X) / 2;
            ymax = Math.Max(1.2, Math.Ceiling(_data.Concat(_second).Max(p => p.Y) * 2) / 2);
        }
        double X(double v) => L + (v - xmin) / (xmax - xmin) * w;
        double Y(double v) => T + h - (v - ymin) / (ymax - ymin) * h;
        _plot = new Rect(L, T, w, h);
        _xmin = xmin; _xmax = xmax; _ymin = ymin; _ymax = ymax;
        if (_b0 is { } ba && _b1 is { } bb)   // the brush box
        {
            var r = new Rect(new Point(Math.Min(ba.X, bb.X), Math.Min(ba.Y, bb.Y)), new Point(Math.Max(ba.X, bb.X), Math.Max(ba.Y, bb.Y)));
            ctx.FillRectangle(new SolidColorBrush((Tokens.Brush("SelB") as ISolidColorBrush)?.Color ?? Colors.SteelBlue, 0.14), r);
            ctx.DrawRectangle(null, new Pen(Tokens.Brush("SelB"), 1, new DashStyle([4, 3], 0)), r);
        }
        // enough decimals that the three y ticks differ (a flat density still reads 0.3998 / 0.4 / 0.4002)
        var ydec = Math.Clamp((int)Math.Ceiling(-Math.Log10(Math.Max(1e-12, (ymax - ymin) / 2))) + 1, 0, 6);
        string Fmt(double v) => Math.Abs(v) >= 1000 && ydec <= 1 ? v.ToString("0", CultureInfo.InvariantCulture) : v.ToString(ydec <= 2 ? "0.##" : "0." + new string('#', ydec), CultureInfo.InvariantCulture);

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
        foreach (var st in Stages)   // NPT stages tinted, each stage's ensemble written where it fits
        {
            var x0 = Math.Clamp(X(st.From), L, L + w);
            var x1 = Math.Clamp(X(st.To), L, L + w);
            if (x1 <= x0) continue;
            if (st.Shade) ctx.FillRectangle(new SolidColorBrush((Tokens.Brush("AccB") as ISolidColorBrush)?.Color ?? Colors.Orange, 0.09), new Rect(x0, T, x1 - x0, h));
            ctx.DrawLine(new Pen(Grid, 1, new DashStyle([2, 3], 0)), new Point(x0, T), new Point(x0, T + h));
            var ft = new FormattedText(st.Label, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, tf, 9, Label);
            if (ft.Width + 4 < x1 - x0) ctx.DrawText(ft, new Point((x0 + x1 - ft.Width) / 2, T + 2));
        }
        if (Band is { } band)
        {
            var x0 = Math.Clamp(X(band.From), L, L + w);
            var x1 = Math.Clamp(X(band.To), L, L + w);
            if (x1 > x0) ctx.FillRectangle(new SolidColorBrush((Tokens.Brush("SelB") as ISolidColorBrush)?.Color ?? Colors.SteelBlue, 0.10), new Rect(x0, T, x1 - x0, h));
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
        if (Bars)
        {
            var dx = _data.Length > 1 ? _data[1].X - _data[0].X : 1;
            var fill = Tokens.Brush("SelB");
            foreach (var p in _data)
            {
                if (p.Y <= 0) continue;
                var x0 = X(p.X - dx / 2) + 1;
                var x1 = X(p.X + dx / 2) - 1;
                if (x1 > x0) ctx.FillRectangle(fill, new Rect(x0, Y(p.Y), x1 - x0, Y(0) - Y(p.Y)), 2);
            }
        }
        else if (Markers)
        {
            if (Errors is { } err && err.Length == _data.Length)   // ± one standard deviation, with caps
            {
                var ep = new Pen(Dot, 1);
                for (var i = 0; i < _data.Length; ++i)
                {
                    if (!(err[i] > 0)) continue;
                    var (px, py0, py1) = (X(_data[i].X), Y(_data[i].Y - err[i]), Y(_data[i].Y + err[i]));
                    ctx.DrawLine(ep, new Point(px, py0), new Point(px, py1));
                    ctx.DrawLine(ep, new Point(px - 3, py0), new Point(px + 3, py0));
                    ctx.DrawLine(ep, new Point(px - 3, py1), new Point(px + 3, py1));
                }
            }
            foreach (var p in _data) ctx.DrawEllipse(Dot, null, new Point(X(p.X), Y(p.Y)), 2.2, 2.2);
            if (_overlay.Length > 1) ctx.DrawGeometry(null, Curve, Line(_overlay));
        }
        else if (_second.Length > 1)
        {
            ctx.DrawGeometry(null, CurveA, Line(_data));
            ctx.DrawGeometry(null, CurveB, Line(_second));
        }
        else
        {
            ctx.DrawGeometry(null, LineBrush != null ? new Pen(LineBrush, 1.6, lineJoin: PenLineJoin.Round) : Curve, Line(_data));
        }
        if (_third.Length > 1)
        {
            var clipped = _third.Where(p => p.X >= xmin && p.X <= xmax).ToArray();
            if (clipped.Length > 1) ctx.DrawGeometry(null, new Pen(Tokens.Brush("MutedB"), 1.4, new DashStyle([1.5, 2.5], 0), lineCap: PenLineCap.Round), Line(clipped));
        }
        if (_fit.Length > 1)
        {
            var clipped = _fit.Where(p => p.Y >= ymin && p.Y <= ymax && p.X >= xmin && p.X <= xmax).ToArray();
            if (clipped.Length > 1) ctx.DrawGeometry(null, new Pen(Tokens.Brush("SelB"), 2, new DashStyle([6, 4], 0)), Line(clipped));
        }
        if (CursorX is double cx && cx >= xmin && cx <= xmax)
            ctx.DrawLine(new Pen(Tokens.Brush("TextB"), 1, new DashStyle([3, 3], 0)), new Point(X(cx), T), new Point(X(cx), T + h));
        Text(YLabel, L, T - 9);
        if (AutoRange) Text("x: " + XLabel, L + w, T - 9, right: true);
        else Text(XLabel, L + w, T + h + 12 + 0, right: true);
    }
}
