using System;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>The composer's page on screen: fitted to the control, with a ruler in inches, the selected panel outlined in the
/// accent, and a click selecting a panel.</summary>
public sealed class FigureCanvas : Control
{
    public FigureComposition? Composition { get; set; }
    public int Selected { get; set; }
    public event Action<int>? PanelClicked;
    private double _scale = 1;
    private Point _origin;

    public string ScaleText => $"{Math.Round(_scale / 0.75 * 100)} % · {Math.Round(_scale * 72)} px per inch";

    public override void Render(DrawingContext ctx)
    {
        var c = Composition;
        if (c == null) return;
        const double margin = 36;
        _scale = Math.Min((Bounds.Width - 2 * margin) / c.WidthPt, (Bounds.Height - 2 * margin - 18) / c.HeightPt);
        if (_scale <= 0) return;
        _origin = new Point((Bounds.Width - c.WidthPt * _scale) / 2, margin + 18 + (Bounds.Height - 2 * margin - 18 - c.HeightPt * _scale) / 2);
        // ruler
        var dim = Tokens.Brush("DimB");
        var tf = new Typeface(Tokens.Mono);
        for (var inch = 0; inch <= Math.Floor(c.WidthIn); ++inch)
        {
            var x = _origin.X + inch * 72 * _scale;
            ctx.DrawLine(new Pen(dim, 1), new Point(x, _origin.Y - 6), new Point(x, _origin.Y - 2));
            var ft = new FormattedText($"{inch} in", System.Globalization.CultureInfo.InvariantCulture, FlowDirection.LeftToRight, tf, 10, dim);
            ctx.DrawText(ft, new Point(x - ft.Width / 2, _origin.Y - 20));
        }
        using (ctx.PushTransform(Matrix.CreateTranslation(_origin.X, _origin.Y)))
        {
            ctx.DrawRectangle(null, new Pen(Tokens.Brush("LineB"), 1), new Rect(-0.5, -0.5, c.WidthPt * _scale + 1, c.HeightPt * _scale + 1));
            c.Draw(ctx, _scale);
            if (Selected >= 0 && Selected < c.Rows * c.Cols)
                ctx.DrawRectangle(null, new Pen(Tokens.Brush("AccB"), 1.5), c.PanelRect(Selected, _scale).Inflate(2), 3, 3);
        }
    }

    protected override void OnPointerPressed(PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        var c = Composition;
        if (c == null) return;
        var p = e.GetPosition(this) - _origin;
        for (var k = 0; k < c.Rows * c.Cols; ++k)
            if (c.PanelRect(k, _scale).Contains(new Point(p.X, p.Y))) { PanelClicked?.Invoke(k); return; }
    }
}
