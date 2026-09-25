using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>A chain's sequence as a strip of cells, one per unit, coloured A, B, C … (design chain colours).</summary>
public sealed class SequenceStrip : Control
{
    private static readonly string[] Colours = ["#F0A83C", "#6CC4D8", "#DE775D", "#9B7AD5", "#7DC884", "#D6AC5C", "#E9ECEF", "#2271DB"];
    private int[] _units = [];
    /// <summary>Colours for unit 0, 1, … (default: the chain colours); letters are drawn only with the default.</summary>
    public string[]? Palette { get; set; }

    public void SetUnits(int[] units)
    {
        _units = units;
        InvalidateVisual();
    }

    public override void Render(DrawingContext ctx)
    {
        if (_units.Length == 0) return;
        var w = Bounds.Width / _units.Length;
        var gap = w > 5 ? 1.5 : 0;
        for (var i = 0; i < _units.Length; i++)
        {
            var pal = Palette ?? Colours;
            var b = new SolidColorBrush(Color.Parse(pal[_units[i] % pal.Length]));
            ctx.FillRectangle(b, new Rect(i * w, 0, Math.Max(0.5, w - gap), Bounds.Height), 2);
        }
        if (w >= 12 && Palette == null)
        {
            var tf = new Typeface(Tokens.Mono);
            for (var i = 0; i < _units.Length; i++)
            {
                var ft = new FormattedText(((char)('A' + _units[i])).ToString(), CultureInfo.InvariantCulture, FlowDirection.LeftToRight, tf, 10, Brushes.Black);
                ctx.DrawText(ft, new Point(i * w + (w - ft.Width) / 2, (Bounds.Height - ft.Height) / 2));
            }
        }
    }
}
