using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>The ring menu (design/boards/SelectionBar): wedges of actions around the pointer, the one the pointer points at
/// lit; the selection's count in the hub. Drawn, not templated: it opens under a held right button and follows the flick.</summary>
public sealed class RingMenu : Control
{
    public const double Radius = 86, Hub = 30, Outer = Radius + 26;
    private IReadOnlyList<RingItem> _items = [];
    private string _hub = "";
    private int _hot = -1;
    private Point _centre;

    public RingMenu() { IsHitTestVisible = false; IsVisible = false; }

    public IReadOnlyList<RingItem> Items => _items;
    public int Hot => _hot;
    public Point Centre => _centre;

    public void Open(Point centre, IReadOnlyList<RingItem> items, string hub)
    {
        _centre = centre; _items = items; _hub = hub; _hot = -1;
        IsVisible = true;
        InvalidateVisual();
    }

    public void Close() { IsVisible = false; _hot = -1; }

    /// <summary>The wedge the pointer points at (beyond the hub), else -1.</summary>
    public int Aim(Point p)
    {
        var d = p - _centre;
        var r = Math.Sqrt(d.X * d.X + d.Y * d.Y);
        var hot = -1;
        if (r > Hub && _items.Count > 0)
        {
            // wedge k is centred on -90° + 360°·k/n (the first at the top), clockwise on screen
            var a = Math.Atan2(d.Y, d.X) * 180 / Math.PI + 90;
            var n = _items.Count;
            hot = (int)Math.Floor(((a + 180.0 / n) % 360 + 360) % 360 / (360.0 / n));
        }
        if (hot != _hot) { _hot = hot; InvalidateVisual(); }
        return hot;
    }

    private IBrush B(string key, IBrush fallback) => this.TryFindResource(key, ActualThemeVariant, out var r) && r is IBrush b ? b : fallback;

    public override void Render(DrawingContext ctx)
    {
        if (_items.Count == 0) return;
        var bg = B("Bg2B", Brushes.DimGray);
        var line = B("LineB", Brushes.Gray);
        var text = B("TextB", Brushes.White);
        var acc = B("AccB", Brushes.Orange);
        var sel = B("SelB", Brushes.SkyBlue);
        var hub = B("Bg1B", Brushes.Black);
        var c = _centre;
        var bgc = (bg as ISolidColorBrush)?.Color ?? Colors.DimGray;
        ctx.DrawEllipse(new SolidColorBrush(bgc, 0.95), new Pen(line, 1), c, Outer, Outer);
        var n = _items.Count;
        var penLine = new Pen(line, 1);
        for (var k = 0; k < n; ++k)
        {
            var a0 = (-90 + 360.0 * k / n - 180.0 / n) * Math.PI / 180;
            var a1 = (-90 + 360.0 * k / n + 180.0 / n) * Math.PI / 180;
            if (k == _hot)
            {
                var g = new StreamGeometry();
                using (var s = g.Open())
                {
                    s.BeginFigure(c + new Point(Hub * Math.Cos(a0), Hub * Math.Sin(a0)), true);
                    s.LineTo(c + new Point(Outer * Math.Cos(a0), Outer * Math.Sin(a0)));
                    s.ArcTo(c + new Point(Outer * Math.Cos(a1), Outer * Math.Sin(a1)), new Size(Outer, Outer), 0, false, SweepDirection.Clockwise);
                    s.LineTo(c + new Point(Hub * Math.Cos(a1), Hub * Math.Sin(a1)));
                    s.ArcTo(c + new Point(Hub * Math.Cos(a0), Hub * Math.Sin(a0)), new Size(Hub, Hub), 0, false, SweepDirection.CounterClockwise);
                    s.EndFigure(true);
                }
                var ac = (acc as ISolidColorBrush)?.Color ?? Colors.Orange;
                ctx.DrawGeometry(new SolidColorBrush(ac, 0.2), null, g);
            }
            ctx.DrawLine(penLine, c + new Point(Hub * Math.Cos(a0), Hub * Math.Sin(a0)), c + new Point(Outer * Math.Cos(a0), Outer * Math.Sin(a0)));
        }
        for (var k = 0; k < n; ++k)
        {
            var a = (-90 + 360.0 * k / n) * Math.PI / 180;
            var p = c + new Point((Radius - 4) * Math.Cos(a), (Radius - 4) * Math.Sin(a));
            var col = k == _hot ? acc : text;
            if (this.TryFindResource("Icon." + _items[k].Icon, ActualThemeVariant, out var res) && res is Geometry geo)
                using (ctx.PushTransform(Matrix.CreateScale(16 / 24.0, 16 / 24.0) * Matrix.CreateTranslation(p.X - 8, p.Y - 16)))
                    ctx.DrawGeometry(null, new Pen(col, 1.7 * 24 / 16.0, lineCap: PenLineCap.Round, lineJoin: PenLineJoin.Round), geo);
            var ft = new FormattedText(_items[k].Word, System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface("Inter"), 10.5, col);
            ctx.DrawText(ft, new Point(p.X - ft.Width / 2, p.Y + 2));
        }
        ctx.DrawEllipse(hub, new Pen(sel, 1), c, Hub - 4, Hub - 4);
        if (_hub.Length > 0)
        {
            var mono = this.TryFindResource("Mono", ActualThemeVariant, out var m) && m is FontFamily ff ? ff : FontFamily.Default;
            var ft = new FormattedText(_hub, System.Globalization.CultureInfo.CurrentCulture, FlowDirection.LeftToRight, new Typeface(mono), 11, sel);
            ctx.DrawText(ft, new Point(c.X - ft.Width / 2, c.Y - ft.Height / 2));
        }
    }
}
