using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Media;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>Brush to select (design/boards/BrushSelect): one atom property as a histogram — every atom in grey, the
/// selected ones over it in the selection colour, the brushed range as a band. Drag across the bars to brush a range;
/// a click without dragging clears it.</summary>
public sealed class BrushHistogram : Control
{
    public static readonly StyledProperty<BrushHist?> HistProperty = AvaloniaProperty.Register<BrushHistogram, BrushHist?>(nameof(Hist));
    static BrushHistogram() { AffectsRender<BrushHistogram>(HistProperty); }
    public BrushHist? Hist { get => GetValue(HistProperty); set => SetValue(HistProperty, value); }

    private double? _dragFrom;
    private double _dragTo;

    protected override void OnPropertyChanged(AvaloniaPropertyChangedEventArgs e)
    {
        base.OnPropertyChanged(e);
        if (e.Property == HistProperty)
        {
            if (e.OldValue is BrushHist o) o.Changed -= InvalidateVisual;
            if (e.NewValue is BrushHist n) n.Changed += InvalidateVisual;
        }
    }

    private IBrush Res(string key, IBrush fallback) => this.TryFindResource(key, ActualThemeVariant, out var r) && r is IBrush b ? b : fallback;

    public override void Render(DrawingContext ctx)
    {
        var h = Hist;
        var w = Bounds.Width;
        var ht = Bounds.Height;
        if (h == null || h.Bins.Length == 0 || w <= 0) return;
        var all = Res("LineB", Brushes.DimGray);
        var sel = Res("SelB", Brushes.DeepSkyBlue);
        var acc = Res("AccB", Brushes.Orange);
        var max = Math.Max(1, h.Bins.Max());
        var n = h.Bins.Length;
        var bw = w / n;
        // the brushed band (or the one being dragged)
        double? a = null, b = null;
        if (_dragFrom is { } f) { a = Math.Min(f, _dragTo); b = Math.Max(f, _dragTo); }
        else if (h.BrushFrom is { } lo && h.BrushTo is { } hi) { a = h.FractionOf(lo, false); b = h.FractionOf(hi, true); }
        if (a is { } x0 && b is { } x1)
        {
            var c = (acc as ISolidColorBrush)?.Color ?? Colors.Orange;
            ctx.FillRectangle(new SolidColorBrush(c, 0.16), new Rect(x0 * w, 0, Math.Max(2, (x1 - x0) * w), ht));
            var pen = new Pen(acc, 1, dashStyle: new DashStyle([3, 2], 0));
            ctx.DrawLine(pen, new Point(x0 * w, 0), new Point(x0 * w, ht));
            ctx.DrawLine(pen, new Point(x1 * w, 0), new Point(x1 * w, ht));
        }
        // a square-root scale, so a few atoms in a bin still show beside thousands
        double Hgt(int c) => c <= 0 ? 0 : Math.Max(1.5, Math.Sqrt(c / (double)max) * (ht - 2));
        for (var k = 0; k < n; ++k)
        {
            var x = k * bw + (h.Categorical ? bw * 0.15 : 0.5);
            var wd = h.Categorical ? bw * 0.7 : Math.Max(1, bw - 1);
            var y = Hgt(h.Bins[k]);
            if (y > 0) ctx.FillRectangle(all, new Rect(x, ht - y, wd, y));
            var ys = Hgt(h.Selected.Length > k ? h.Selected[k] : 0);
            if (ys > 0) ctx.FillRectangle(sel, new Rect(x, ht - ys, wd, ys));
        }
    }

    private double Frac(Point p) => Math.Clamp(p.X / Math.Max(1, Bounds.Width), 0, 1);

    protected override void OnPointerPressed(PointerPressedEventArgs e)
    {
        base.OnPointerPressed(e);
        if (Hist == null) return;
        _dragFrom = _dragTo = Frac(e.GetPosition(this));
        e.Pointer.Capture(this);
        e.Handled = true;
    }

    protected override void OnPointerMoved(PointerEventArgs e)
    {
        base.OnPointerMoved(e);
        if (_dragFrom == null) return;
        _dragTo = Frac(e.GetPosition(this));
        InvalidateVisual();
    }

    protected override void OnPointerReleased(PointerReleasedEventArgs e)
    {
        base.OnPointerReleased(e);
        if (_dragFrom is not { } f || Hist is not { } h) return;
        var to = Frac(e.GetPosition(this));
        _dragFrom = null;
        e.Pointer.Capture(null);
        if (Math.Abs(to - f) * Bounds.Width < 3 && !h.Categorical) h.SetBrush(null, null);       // a click: no brush
        else if (h.Categorical) h.BrushCategories(Math.Min(f, to), Math.Max(f, to));              // the bars touched
        else h.SetBrush(h.ValueAt(Math.Min(f, to)), h.ValueAt(Math.Max(f, to)));
        InvalidateVisual();
    }
}
