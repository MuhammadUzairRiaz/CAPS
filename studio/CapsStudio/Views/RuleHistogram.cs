using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Media;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>Bond rules (design/boards/BondRules): one element pair's distances as a histogram, the bonded side of the
/// cut-off in the accent, with a white handle to drag the cut-off; a pair never bonded is drawn dim.</summary>
public sealed class RuleHistogram : Control
{
    public static readonly StyledProperty<RulePair?> PairProperty = AvaloniaProperty.Register<RuleHistogram, RulePair?>(nameof(Pair));
    static RuleHistogram() { AffectsRender<RuleHistogram>(PairProperty); }
    public RulePair? Pair { get => GetValue(PairProperty); set => SetValue(PairProperty, value); }
    private bool _drag;

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs e)
    {
        base.OnPropertyChanged(e);
        if (e.Property == PairProperty)
        {
            if (e.OldValue is RulePair o) o.Changed -= InvalidateVisual;
            if (e.NewValue is RulePair n) n.Changed += InvalidateVisual;
        }
    }

    private IBrush Res(string key, IBrush fallback) => this.TryFindResource(key, ActualThemeVariant, out var r) && r is IBrush b ? b : fallback;
    private double X(RulePair p, double r) => (r - p.Lo) / Math.Max(1e-9, p.Hi - p.Lo) * Bounds.Width;
    private double R(RulePair p, double x) => p.Lo + Math.Clamp(x / Math.Max(1, Bounds.Width), 0, 1) * (p.Hi - p.Lo);

    public override void Render(DrawingContext ctx)
    {
        if (Pair is not { } p || p.Counts.Length == 0 || Bounds.Width <= 0) return;
        var w = Bounds.Width;
        var h = Bounds.Height;
        var line = Res("LineB", Brushes.DimGray);
        var acc = Res("AccB", Brushes.Orange);
        var text = Res("TextB", Brushes.White);
        var max = Math.Max(1, p.Counts.Max());
        var bw = w / p.Counts.Length;
        for (var k = 0; k < p.Counts.Length; ++k)
        {
            var c = p.Counts[k];
            if (c <= 0) continue;
            var y = Math.Max(1.5, Math.Sqrt(c / (double)max) * (h - 4));   // a square-root scale: the few far pairs still show
            var r = p.Lo + (k + 0.5) * p.Bin;
            var bonded = !p.Never && r <= p.Cutoff;
            ctx.FillRectangle(bonded ? acc : line, new Rect(k * bw + 0.5, h - y, Math.Max(1, bw - 1), y));
        }
        if (p.Never) return;
        var x = X(p, p.Cutoff);
        ctx.DrawLine(new Pen(text, 1.5), new Point(x, 0), new Point(x, h));
        ctx.FillRectangle(text, new Rect(x - 4, 0, 8, 10), 2);
    }

    protected override void OnPointerPressed(PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        if (Pair is not { Never: false } p) return;
        _drag = true;
        p.Cutoff = Math.Round(R(p, e.GetPosition(this).X), 2);
        e.Pointer.Capture(this);
        e.Handled = true;
    }
    protected override void OnPointerMoved(PointerEventArgs e)
    {
        base.OnPointerMoved(e);
        if (_drag && Pair is { } p) p.Cutoff = Math.Round(R(p, e.GetPosition(this).X), 2);
    }
    protected override void OnPointerReleased(PointerReleasedEventArgs e)
    {
        base.OnPointerReleased(e);
        if (!_drag) return;
        _drag = false;
        e.Pointer.Capture(null);
        Pair?.Committed?.Invoke();
    }
}
