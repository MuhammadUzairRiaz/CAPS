using System.Globalization;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>A crystal habit (a convex Wulff shape) drawn in orthographic projection: each facet shaded by its tilt to the
/// light and coloured by its form, the form's indices on its largest visible facet; drag to turn it.</summary>
public sealed class HabitView : Control
{
    private double[][] _v = [];
    private int[][] _f = [];
    private int[] _face = [];
    private string[] _labels = [];
    private double _yaw = 0.6, _pitch = 0.45;
    private Point? _drag;
    private static readonly Color[] Palette =
    [
        Color.Parse("#E69F00"), Color.Parse("#56B4E9"), Color.Parse("#009E73"), Color.Parse("#F0E442"), Color.Parse("#0072B2"),
        Color.Parse("#D55E00"), Color.Parse("#CC79A7"), Color.Parse("#999999"),
    ];

    /// <summary>The shape: vertices, facets (vertex indices), each facet's form index (−1: a cap with no form) and the forms' labels.</summary>
    public void Set(double[][] vertices, int[][] facets, int[] facetForm, string[] formLabels)
    {
        _v = vertices; _f = facets; _face = facetForm; _labels = formLabels;
        InvalidateVisual();
    }

    protected override void OnPointerPressed(PointerPressedEventArgs e) { _drag = e.GetPosition(this); e.Pointer.Capture(this); }
    protected override void OnPointerReleased(PointerReleasedEventArgs e) { _drag = null; e.Pointer.Capture(null); }
    protected override void OnPointerMoved(PointerEventArgs e)
    {
        if (_drag is not Point p) return;
        var q = e.GetPosition(this);
        _yaw += (q.X - p.X) * 0.01;
        _pitch = Math.Clamp(_pitch + (q.Y - p.Y) * 0.01, -1.5, 1.5);
        _drag = q;
        InvalidateVisual();
    }

    private (double X, double Y, double Z) Rot(double[] v)
    {
        double cy = Math.Cos(_yaw), sy = Math.Sin(_yaw), cp = Math.Cos(_pitch), sp = Math.Sin(_pitch);
        var x = cy * v[0] + sy * v[1];
        var y = -sy * v[0] + cy * v[1];
        var z = v[2];
        return (x, cp * y - sp * z, sp * y + cp * z);   // screen x, screen up, towards the viewer
    }

    public override void Render(DrawingContext ctx)
    {
        var b = Bounds;
        ctx.FillRectangle(Brushes.Transparent, new Rect(b.Size));   // the whole area takes the drag
        if (_v.Length == 0 || _f.Length == 0) return;
        var r = _v.Select(Rot).ToArray();
        var ext = r.Max(p => Math.Sqrt(p.X * p.X + p.Y * p.Y));
        if (ext <= 0) return;
        var s = 0.42 * Math.Min(b.Width, b.Height) / ext;
        Point P((double X, double Y, double Z) p) => new(b.Width / 2 + s * p.X, b.Height / 2 - s * p.Y);
        var edge = new Pen(Tokens.Brush("TextB"), 0.8);
        var font = new Typeface(Tokens.Mono);
        var order = Enumerable.Range(0, _f.Length).Select(k =>
        {
            var pts = _f[k].Select(i => r[i]).ToArray();
            // the facet's normal in view space (Newell)
            double nx = 0, ny = 0, nz = 0;
            for (var i = 0; i < pts.Length; ++i)
            {
                var a = pts[i]; var c = pts[(i + 1) % pts.Length];
                nx += (a.Y - c.Y) * (a.Z + c.Z); ny += (a.Z - c.Z) * (a.X + c.X); nz += (a.X - c.X) * (a.Y + c.Y);
            }
            var nn = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            return (k, pts, nx: nx / nn, ny: ny / nn, nz: nz / nn, depth: pts.Average(p => p.Z));
        }).Where(x => x.nz > 1e-6).OrderBy(x => x.depth).ToArray();
        var labelAt = new Dictionary<int, (double area, Point at)>();
        foreach (var (k, pts, nx, ny, nz, _) in order)
        {
            var form = k < _face.Length ? _face[k] : -1;
            var baseC = form >= 0 ? Palette[form % Palette.Length] : Color.Parse("#555555");
            var light = 0.45 + 0.55 * Math.Clamp(0.35 * nx + 0.45 * ny + 0.82 * nz, 0, 1);
            var c = Color.FromRgb((byte)(baseC.R * light), (byte)(baseC.G * light), (byte)(baseC.B * light));
            var g = new StreamGeometry();
            using (var gc = g.Open())
            {
                gc.BeginFigure(P(pts[0]), true);
                for (var i = 1; i < pts.Length; ++i) gc.LineTo(P(pts[i]));
                gc.EndFigure(true);
            }
            ctx.DrawGeometry(new SolidColorBrush(c), edge, g);
            if (form < 0) continue;
            var sp = pts.Select(P).ToArray();
            double area = 0;
            for (var i = 0; i < sp.Length; ++i) area += sp[i].X * sp[(i + 1) % sp.Length].Y - sp[(i + 1) % sp.Length].X * sp[i].Y;
            area = Math.Abs(area) / 2;
            var cen = new Point(sp.Average(p => p.X), sp.Average(p => p.Y));
            if (!labelAt.TryGetValue(form, out var best) || area > best.area) labelAt[form] = (area, cen);
        }
        foreach (var (form, (area, at)) in labelAt)
        {
            if (area < 1500 || form >= _labels.Length) continue;
            var t = new FormattedText(_labels[form], CultureInfo.InvariantCulture, FlowDirection.LeftToRight, font, 11, Brushes.Black);
            ctx.DrawText(t, new Point(at.X - t.Width / 2, at.Y - t.Height / 2));
        }
    }
}
