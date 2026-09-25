using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>Atom labels over the 3D view (design/boards/Appearance): text beside each visible labelled atom, on a
/// translucent plate so it reads over any colour.</summary>
public sealed class LabelOverlay : Control
{
    private List<ViewLabel> _labels = new();
    private static readonly Typeface Face = new("avares://CapsStudio/Assets/Fonts#IBM Plex Mono");

    public LabelOverlay() { IsHitTestVisible = false; }

    public void SetLabels(List<ViewLabel> labels)
    {
        _labels = labels;
        InvalidateVisual();
    }

    private (double X, double Y, double R, string Label)? _lens;
    /// <summary>The all-atom lens outline (design/boards/LensView): a dashed circle and its label, or none.</summary>
    public void SetLens((double X, double Y, double R, string Label)? lens)
    {
        if (_lens == lens) return;
        _lens = lens;
        InvalidateVisual();
    }

    public override void Render(DrawingContext ctx)
    {
        if (_lens is { } lz)
        {
            var sel = Tokens.Brush("SelB");
            ctx.DrawEllipse(null, new Pen(sel, 1.5, new DashStyle([5, 4], 0)), new Point(lz.X, lz.Y), lz.R, lz.R);
            var ft = new FormattedText(lz.Label, System.Globalization.CultureInfo.InvariantCulture, FlowDirection.LeftToRight, Face, 10.5, sel);
            var x = lz.X + lz.R * 0.72 - 4;
            var y = lz.Y - lz.R * 0.72 - ft.Height - 4;
            ctx.FillRectangle(Tokens.Brush("Bg1B"), new Rect(x - 6, y - 2, ft.Width + 12, ft.Height + 4), 4);
            ctx.DrawRectangle(null, new Pen(sel, 1), new Rect(x - 6, y - 2, ft.Width + 12, ft.Height + 4), 4, 4);
            ctx.DrawText(ft, new Point(x, y));
        }
        if (_labels.Count == 0) return;
        var ink = Tokens.Brush("TextB");
        var plate = new SolidColorBrush(Color.FromArgb(0xB0, 0x16, 0x19, 0x1C));
        var light = Application.Current?.ActualThemeVariant == Avalonia.Styling.ThemeVariant.Light;
        if (light) plate = new SolidColorBrush(Color.FromArgb(0xC8, 0xFF, 0xFF, 0xFF));
        foreach (var l in _labels)
        {
            var ft = new FormattedText(l.Text, System.Globalization.CultureInfo.InvariantCulture, FlowDirection.LeftToRight, Face, 10.5, ink);
            var x = l.X + 5;
            var y = l.Y - ft.Height - 2;
            ctx.FillRectangle(plate, new Rect(x - 2, y, ft.Width + 4, ft.Height), 2);
            ctx.DrawText(ft, new Point(x, y));
        }
    }
}
