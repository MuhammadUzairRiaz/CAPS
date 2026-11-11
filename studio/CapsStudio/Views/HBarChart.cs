using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>Horizontal bars with a label on the left and the value at the bar's end, and one dashed guide
/// (design/boards/SolventScreen: χ per solvent, the χ = 0.5 line).</summary>
public sealed class HBarChart : Control
{
    private (string Label, double Value, string Colour, string Text)[] _bars = [];
    private double? _guide;
    private string _guideLabel = "";
    public double LabelWidth { get; set; } = 96;

    public void Set((string Label, double Value, string Colour, string Text)[] bars, double? guide = null, string guideLabel = "")
    {
        _bars = bars;
        _guide = guide;
        _guideLabel = guideLabel;
        InvalidateVisual();
    }

    public override void Render(DrawingContext ctx)
    {
        if (_bars.Length == 0) return;
        var b = Bounds;
        const double bottom = 22;
        var L = LabelWidth;
        var w = b.Width - L - 12;
        var rowH = (b.Height - bottom) / _bars.Length;
        if (w < 40 || rowH < 8) return;
        var max = Math.Max(_guide ?? 0, _bars.Max(x => x.Value)) * 1.05;
        double X(double v) => L + Math.Clamp(v / max, 0, 1) * w;
        var ui = new Typeface(Tokens.Sans);
        var mono = new Typeface(Tokens.Mono);
        var dim = Tokens.Brush("DimB");
        var muted = Tokens.Brush("MutedB");
        double? gx = _guide is double g ? X(g) : null;
        var values = new (FormattedText Text, double X1, double Y)[_bars.Length];
        for (var i = 0; i < _bars.Length; ++i)
        {
            var (label, value, colour, text) = _bars[i];
            var y = i * rowH;
            var bh = Math.Min(18, rowH * 0.62);
            var brush = colour.StartsWith('#') ? new SolidColorBrush(Color.Parse(colour)) : Tokens.Brush(colour);
            var lt = new FormattedText(label, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, ui, 11, muted);
            ctx.DrawText(lt, new Point(L - 8 - lt.Width, y + (rowH - lt.Height) / 2));
            var x1 = X(value);
            ctx.FillRectangle(brush, new Rect(L, y + (rowH - bh) / 2, Math.Max(2, x1 - L), bh), 3);
            values[i] = (new FormattedText(text, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, mono, 10.5, Tokens.Brush("TextB")), x1, y);
        }
        // The guide goes under the value labels, and a label that would sit on it moves clear of it.
        if (gx is double gl)
        {
            var pen = new Pen(dim, 1, new DashStyle([4, 4], 0));
            ctx.DrawLine(pen, new Point(gl, 0), new Point(gl, b.Height - bottom + 4));
            var gt = new FormattedText(_guideLabel, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, mono, 10, dim);
            ctx.DrawText(gt, new Point(gl - gt.Width / 2, b.Height - bottom + 6));
        }
        const double pad = 6, clear = 4;
        foreach (var (vt, x1, y) in values)
        {
            var right = L + w;
            var lx = x1 + pad + vt.Width > right ? x1 - pad - vt.Width : x1 + pad;
            if (gx is double gv && lx < gv + clear && lx + vt.Width > gv - clear)
            {
                if (gv + pad + vt.Width <= right) lx = gv + pad;                 // beyond the guide
                else if (x1 - pad - vt.Width >= L + 2) lx = x1 - pad - vt.Width; // inside the bar
            }
            ctx.DrawText(vt, new Point(lx, y + (rowH - vt.Height) / 2));
        }
    }
}
