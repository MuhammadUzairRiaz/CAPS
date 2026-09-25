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
            var vt = new FormattedText(text, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, mono, 10.5, Tokens.Brush("TextB"));
            var inside = x1 + 6 + vt.Width > L + w;
            ctx.DrawText(vt, new Point(inside ? x1 - 6 - vt.Width : x1 + 6, y + (rowH - vt.Height) / 2));
        }
        if (_guide is double g)
        {
            var pen = new Pen(dim, 1, new DashStyle([4, 4], 0));
            ctx.DrawLine(pen, new Point(X(g), 0), new Point(X(g), b.Height - bottom + 4));
            var gt = new FormattedText(_guideLabel, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, mono, 10, dim);
            ctx.DrawText(gt, new Point(X(g) - gt.Width / 2, b.Height - bottom + 6));
        }
    }
}
