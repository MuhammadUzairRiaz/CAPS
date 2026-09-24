using Avalonia;
using Avalonia.Controls;
using Avalonia.Controls.Documents;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>A CAPS icon (design system, 24 × 24 stroked, stroke 1.7) drawn from the geometry "Icon.{Kind}" in
/// Themes/Icons.axaml. The colour follows the inherited Foreground, so buttons colour their icons through styles.</summary>
public sealed class Icon : Control
{
    public static readonly StyledProperty<string> KindProperty = AvaloniaProperty.Register<Icon, string>(nameof(Kind), "hex");
    public static readonly StyledProperty<double> SizeProperty = AvaloniaProperty.Register<Icon, double>(nameof(Size), 16);
    public static readonly StyledProperty<double> StrokeWidthProperty = AvaloniaProperty.Register<Icon, double>(nameof(StrokeWidth), 1.7);
    public static readonly StyledProperty<IBrush?> ForegroundProperty = TextElement.ForegroundProperty.AddOwner<Icon>();

    static Icon()
    {
        AffectsRender<Icon>(KindProperty, SizeProperty, StrokeWidthProperty, ForegroundProperty);
        AffectsMeasure<Icon>(SizeProperty);
    }

    public string Kind { get => GetValue(KindProperty); set => SetValue(KindProperty, value); }
    public double Size { get => GetValue(SizeProperty); set => SetValue(SizeProperty, value); }
    public double StrokeWidth { get => GetValue(StrokeWidthProperty); set => SetValue(StrokeWidthProperty, value); }
    public IBrush? Foreground { get => GetValue(ForegroundProperty); set => SetValue(ForegroundProperty, value); }

    protected override Size MeasureOverride(Size availableSize) => new(Size, Size);

    public override void Render(DrawingContext ctx)
    {
        if (!this.TryFindResource("Icon." + Kind, ActualThemeVariant, out var res) || res is not Geometry g) return;
        var s = Size / 24.0;
        using (ctx.PushTransform(Matrix.CreateScale(s, s)))
            ctx.DrawGeometry(null, new Pen(Foreground ?? Brushes.Gray, StrokeWidth, lineCap: PenLineCap.Round, lineJoin: PenLineJoin.Round), g);
    }
}
