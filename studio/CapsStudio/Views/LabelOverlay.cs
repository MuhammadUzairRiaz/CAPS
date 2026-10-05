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

    private LabelStyle _style = new("IBM Plex Mono", 10.5, false, true, null, 0xFFF0A83C);
    public void SetLabels(List<ViewLabel> labels, LabelStyle? style = null)
    {
        _labels = labels;
        if (style != null) _style = style;
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

    private List<Point> _lasso = new();
    private (Point A, Point B)? _arrow;
    private List<MonitorMark> _monitors = new();
    /// <summary>The lasso being drawn (view points), or an empty list.</summary>
    public void SetLasso(List<Point> pts) { _lasso = pts; InvalidateVisual(); }
    /// <summary>The move tool's drag, from where it started to the cursor, or none.</summary>
    public void SetArrow((Point A, Point B)? a) { _arrow = a; InvalidateVisual(); }
    /// <summary>Pinned measurements: a dashed line and label for each distance, a label for angles and dihedrals.</summary>
    public void SetMonitors(List<MonitorMark> m) { _monitors = m; InvalidateVisual(); }
    /// <summary>Look: angle and torsion monitors drawn with their arcs.</summary>
    public bool Arcs { get; set; }

    /// <summary>An angle's arc at the vertex v between the directions to a and c (the short way), with short legs.</summary>
    private static void DrawArc(DrawingContext ctx, IPen pen, IPen legs, Point v, Point a, Point c, double r)
    {
        double a0 = Math.Atan2(a.Y - v.Y, a.X - v.X), a1 = Math.Atan2(c.Y - v.Y, c.X - v.X);
        var sweep = a1 - a0;
        while (sweep > Math.PI) sweep -= 2 * Math.PI;
        while (sweep < -Math.PI) sweep += 2 * Math.PI;
        ctx.DrawLine(legs, v, new Point(v.X + Math.Cos(a0) * r * 1.6, v.Y + Math.Sin(a0) * r * 1.6));
        ctx.DrawLine(legs, v, new Point(v.X + Math.Cos(a1) * r * 1.6, v.Y + Math.Sin(a1) * r * 1.6));
        var geo = new StreamGeometry();
        using (var g = geo.Open())
        {
            g.BeginFigure(new Point(v.X + Math.Cos(a0) * r, v.Y + Math.Sin(a0) * r), false);
            const int steps = 24;
            for (var k = 1; k <= steps; ++k)
            {
                var t = a0 + sweep * k / steps;
                g.LineTo(new Point(v.X + Math.Cos(t) * r, v.Y + Math.Sin(t) * r));
            }
            g.EndFigure(false);
        }
        ctx.DrawGeometry(null, pen, geo);
    }

    public override void Render(DrawingContext ctx)
    {
        var acc = Tokens.Brush("AccB");
        if (_lasso.Count > 1)
        {
            var geo = new StreamGeometry();
            using (var g = geo.Open())
            {
                g.BeginFigure(_lasso[0], true);
                foreach (var q in _lasso.Skip(1)) g.LineTo(q);
                g.EndFigure(true);
            }
            ctx.DrawGeometry(new SolidColorBrush(Color.FromArgb(0x24, 0xF0, 0xA8, 0x3C)), new Pen(acc, 1.4, new DashStyle([4, 3], 0)), geo);
        }
        if (_arrow is { } ar)
        {
            var pen = new Pen(acc, 2);
            ctx.DrawLine(pen, ar.A, ar.B);
            ctx.DrawEllipse(acc, null, ar.B, 4, 4);
            ctx.DrawEllipse(null, pen, ar.A, 5, 5);
        }
        foreach (var m in _monitors)
        {
            var at = m.Line ? new Point((m.X0 + m.X1) / 2, (m.Y0 + m.Y1) / 2) : new Point(m.X0, m.Y0);
            if (m.Line)
            {
                ctx.DrawLine(new Pen(acc, 1.4, new DashStyle([5, 4], 0)), new Point(m.X0, m.Y0), new Point(m.X1, m.Y1));
                ctx.DrawEllipse(acc, null, new Point(m.X0, m.Y0), 3, 3);
                ctx.DrawEllipse(acc, null, new Point(m.X1, m.Y1), 3, 3);
            }
            if (Arcs && m.Xs is { } xs && m.Ys is { } ys)
            {
                var arcPen = new Pen(acc, 1.6);
                var legPen = new Pen(acc, 1, new DashStyle([3, 3], 0));
                if (xs.Length == 3) DrawArc(ctx, arcPen, legPen, new Point(xs[1], ys[1]), new Point(xs[0], ys[0]), new Point(xs[2], ys[2]), 16);
                else if (xs.Length == 4)
                {
                    // a torsion: the central bond, and the arc between the two outer bonds seen down it
                    var b = new Point(xs[1], ys[1]); var c = new Point(xs[2], ys[2]);
                    ctx.DrawLine(legPen, b, c);
                    var mid = new Point((b.X + c.X) / 2, (b.Y + c.Y) / 2);
                    DrawArc(ctx, arcPen, legPen, mid, new Point(mid.X + xs[0] - b.X, mid.Y + ys[0] - b.Y), new Point(mid.X + xs[3] - c.X, mid.Y + ys[3] - c.Y), 18);
                }
            }
            var ft = new FormattedText(m.Text, System.Globalization.CultureInfo.InvariantCulture, FlowDirection.LeftToRight, Face, 11, acc);
            var r = new Rect(at.X + 8, at.Y - ft.Height / 2 - 3, ft.Width + 12, ft.Height + 6);
            ctx.FillRectangle(Tokens.Brush("Bg1B"), r, 4);
            ctx.DrawRectangle(null, new Pen(acc, 1), r, 4, 4);
            ctx.DrawText(ft, new Point(r.X + 6, r.Y + 3));
        }
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
        LabelDrawing.Draw(ctx, _labels, _style);
    }
}

/// <summary>Draws atom and bond labels in a style (the Studio view's overlay and the view window).</summary>
public static class LabelDrawing
{
    // the Studio's own fonts by their embedded resources, others by name (system fonts)
    private static FontFamily Family(string name) => name switch
    {
        "IBM Plex Mono" => new FontFamily("avares://CapsStudio/Assets/Fonts#IBM Plex Mono"),
        "IBM Plex Sans" => new FontFamily("avares://CapsStudio/Assets/Fonts#IBM Plex Sans"),
        "Inter" => new FontFamily("fonts:Inter#Inter"),
        _ => new FontFamily(name),
    };

    public static void Draw(DrawingContext ctx, IReadOnlyList<ViewLabel> labels, LabelStyle style)
    {
        if (labels.Count == 0) return;
        var ink = Tokens.Brush("TextB");
        var light = Application.Current?.ActualThemeVariant == Avalonia.Styling.ThemeVariant.Light;
        IBrush plate = light ? new SolidColorBrush(Color.FromArgb(0xC8, 0xFF, 0xFF, 0xFF)) : new SolidColorBrush(Color.FromArgb(0xB0, 0x16, 0x19, 0x1C));
        var face = new Typeface(Family(style.Font), FontStyle.Normal, style.Bold ? FontWeight.Bold : FontWeight.Normal);
        var brushes = new Dictionary<uint, IBrush>();
        IBrush Brush(uint? argb) => argb is { } c && c != 0 ? (brushes.TryGetValue(c, out var b) ? b : brushes[c] = new SolidColorBrush(Color.FromUInt32(c))) : ink;
        foreach (var l in labels)
        {
            var brush = l.Bond ? Brush(style.BondArgb) : Brush(l.Argb != 0 ? l.Argb : style.AtomArgb);
            var ft = new FormattedText(l.Text, System.Globalization.CultureInfo.InvariantCulture, FlowDirection.LeftToRight, face, style.Size, brush);
            // atoms: beside the atom; bonds: centred on the bond's midpoint
            var x = l.Bond ? l.X - ft.Width / 2 : l.X + 5;
            var y = l.Bond ? l.Y - ft.Height / 2 : l.Y - ft.Height - 2;
            if (style.Plate) ctx.FillRectangle(plate, new Rect(x - 2, y, ft.Width + 4, ft.Height), 2);
            ctx.DrawText(ft, new Point(x, y));
        }
    }
}
