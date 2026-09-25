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

    public override void Render(DrawingContext ctx)
    {
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
