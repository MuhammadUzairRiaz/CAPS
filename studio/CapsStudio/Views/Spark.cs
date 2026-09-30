using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>A small profile line (a layer's atoms along z): values 0 … 1 across the width, the area under it tinted.</summary>
public sealed class Spark : Control
{
    public static readonly StyledProperty<double[]?> ValuesProperty = AvaloniaProperty.Register<Spark, double[]?>(nameof(Values));
    public static readonly StyledProperty<IBrush?> StrokeProperty = AvaloniaProperty.Register<Spark, IBrush?>(nameof(Stroke));
    static Spark() { AffectsRender<Spark>(ValuesProperty, StrokeProperty); }
    public double[]? Values { get => GetValue(ValuesProperty); set => SetValue(ValuesProperty, value); }
    public IBrush? Stroke { get => GetValue(StrokeProperty); set => SetValue(StrokeProperty, value); }

    public override void Render(DrawingContext ctx)
    {
        var v = Values;
        if (v is not { Length: > 1 }) return;
        var w = Bounds.Width; var h = Bounds.Height;
        var brush = Stroke ?? Brushes.Gray;
        var pts = v.Select((y, k) => new Point(w * k / (v.Length - 1), h - 1 - Math.Clamp(y, 0, 1) * (h - 2))).ToList();
        var area = new StreamGeometry();
        using (var g = area.Open())
        {
            g.BeginFigure(new Point(0, h), true);
            foreach (var p in pts) g.LineTo(p);
            g.LineTo(new Point(w, h));
            g.EndFigure(true);
        }
        var c = (brush as ISolidColorBrush)?.Color ?? Colors.Gray;
        ctx.DrawGeometry(new SolidColorBrush(c, 0.18), null, area);
        var pen = new Pen(brush, 1.2, lineJoin: PenLineJoin.Round);
        for (var k = 1; k < pts.Count; ++k) ctx.DrawLine(pen, pts[k - 1], pts[k]);
    }
}
