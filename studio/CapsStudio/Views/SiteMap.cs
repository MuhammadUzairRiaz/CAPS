using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Media;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>The top face of a sheet seen from above (2D sheets page): the cell repeated 3 × 3 (the middle one solid), its
/// sites marked by kind — fcc (above the third layer), hcp (above the second), top, bridge — and, under the pointer, what
/// a site is: the atom beneath, its in-plane distance r and the height h = √(d² − r²) the chosen termination takes there.</summary>
public sealed class SiteMap : Control
{
    private IReadOnlyList<DftSite> _sites = [];
    private (double Ax, double Ay, double Bx, double By) _cell;
    private DftSite? _hover;
    private readonly List<(Point P, DftSite S)> _drawn = new();

    public void Set(IReadOnlyList<DftSite> sites, (double Ax, double Ay, double Bx, double By) cell)
    {
        _sites = sites.Where(s => s.Face == "top").ToList();
        _cell = cell;
        InvalidateVisual();
    }

    private static IBrush KindBrush(string k) => Tokens.Brush(k switch { "fcc" => "SelB", "hcp" => "AccB", "top" => "MutedB", _ => "DimB" });

    public override void Render(DrawingContext ctx)
    {
        var b = Bounds;
        ctx.FillRectangle(Tokens.Brush("Bg0B"), new Rect(b.Size));
        _drawn.Clear();
        if (_sites.Count == 0 || _cell.Ax <= 0) return;
        // 3 × 3 cells, scaled to fit
        double minx = 1e9, maxx = -1e9, miny = 1e9, maxy = -1e9;
        for (var i = -1; i <= 1; ++i)
            for (var j = -1; j <= 1; ++j)
                foreach (var (fx, fy) in new[] { (0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0) })
                {
                    var x = (i + fx) * _cell.Ax + (j + fy) * _cell.Bx;
                    var y = (i + fx) * _cell.Ay + (j + fy) * _cell.By;
                    minx = Math.Min(minx, x); maxx = Math.Max(maxx, x); miny = Math.Min(miny, y); maxy = Math.Max(maxy, y);
                }
        var sc = Math.Min((b.Width - 24) / (maxx - minx), (b.Height - 40) / (maxy - miny));
        Point P(double x, double y) => new(12 + (x - minx) * sc, b.Height - 28 - (y - miny) * sc);
        var line = new Pen(Tokens.Brush("LineB"), 1);
        for (var i = -1; i <= 1; ++i)
            for (var j = -1; j <= 1; ++j)
            {
                var o = (X: i * _cell.Ax + j * _cell.Bx, Y: i * _cell.Ay + j * _cell.By);
                var g = new StreamGeometry();
                using (var c = g.Open())
                {
                    c.BeginFigure(P(o.X, o.Y), true);
                    c.LineTo(P(o.X + _cell.Ax, o.Y + _cell.Ay));
                    c.LineTo(P(o.X + _cell.Ax + _cell.Bx, o.Y + _cell.Ay + _cell.By));
                    c.LineTo(P(o.X + _cell.Bx, o.Y + _cell.By));
                    c.EndFigure(true);
                }
                ctx.DrawGeometry(i == 0 && j == 0 ? Tokens.Brush("Bg2B") : null, i == 0 && j == 0 ? new Pen(Tokens.Brush("MutedB"), 1.2) : line, g);
                foreach (var s in _sites)
                {
                    var p = P(s.X + o.X, s.Y + o.Y);
                    if (p.X < 0 || p.Y < 0 || p.X > b.Width || p.Y > b.Height - 20) continue;
                    var r = s.Kind == "top" ? 7.0 : s.Kind == "bridge" ? 3.0 : 5.0;
                    ctx.DrawEllipse(s.Kind == "top" ? null : KindBrush(s.Kind), new Pen(KindBrush(s.Kind), s == _hover ? 3 : 1.4), p, r, r);
                    _drawn.Add((p, s));
                }
            }
        // legend and the hovered site
        var face = new Typeface(Tokens.Sans);
        var x0 = 12.0;
        foreach (var k in new[] { "fcc", "hcp", "top", "bridge" })
        {
            ctx.DrawEllipse(k == "top" ? null : KindBrush(k), new Pen(KindBrush(k), 1.4), new Point(x0 + 5, b.Height - 10), 4, 4);
            var t = new FormattedText(k, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, face, 11, Tokens.Brush("MutedB"));
            ctx.DrawText(t, new Point(x0 + 13, b.Height - 18));
            x0 += 20 + t.Width + 10;
        }
        if (_hover is { } h)
        {
            var txt = $"{h.Kind} · beneath {(h.Beneath.Length > 0 ? h.Beneath : "—")} · r {h.R.ToString("0.000", CultureInfo.InvariantCulture)} Å" +
                      (double.IsFinite(h.H) && h.H > 0 ? $" · h {h.H.ToString("0.000", CultureInfo.InvariantCulture)} Å" : double.IsFinite(h.H) && h.H < 0 ? " · bond too short for this site" : "");
            var t = new FormattedText(txt, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, face, 11.5, Tokens.Brush("TextB"));
            ctx.FillRectangle(Tokens.Brush("Bg1B"), new Rect(8, 6, t.Width + 12, t.Height + 6), 4);
            ctx.DrawText(t, new Point(14, 9));
        }
    }

    protected override void OnPointerMoved(PointerEventArgs e)
    {
        base.OnPointerMoved(e);
        var p = e.GetPosition(this);
        DftSite? best = null;
        var bd = 10.0;
        foreach (var (q, s) in _drawn)
        {
            var d = Math.Sqrt((q.X - p.X) * (q.X - p.X) + (q.Y - p.Y) * (q.Y - p.Y));
            if (d < bd) { bd = d; best = s; }
        }
        if (best != _hover) { _hover = best; InvalidateVisual(); }
    }
    protected override void OnPointerExited(PointerEventArgs e) { base.OnPointerExited(e); if (_hover != null) { _hover = null; InvalidateVisual(); } }
}
