using System.Globalization;
using System.Text.Json.Nodes;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;

namespace CapsStudio.Views;

/// <summary>A reaction template's pattern drawn in 2D (design/boards/ReactionTemplate): element labels, map-number badges,
/// reacting atoms ringed; bonds formed by the reaction in green, bonds it breaks dashed red.</summary>
public sealed class TemplateDrawing : Control
{
    private readonly List<(int Map, string Symbol, double X, double Y, bool Reacting)> _atoms = new();
    private readonly List<(int A, int B)> _bonds = new();
    private HashSet<(int, int)> _formed = new(), _broken = new();

    public void Set(JsonNode? drawing, IEnumerable<(int, int)> formed, IEnumerable<(int, int)> broken)
    {
        _atoms.Clear();
        _bonds.Clear();
        _formed = formed.Select(Key).ToHashSet();
        _broken = broken.Select(Key).ToHashSet();
        if (drawing?["atoms"] is JsonArray atoms)
            foreach (var a in atoms)
                _atoms.Add((a!["map"]!.GetValue<int>(), a["symbol"]!.GetValue<string>(), a["x"]!.GetValue<double>(), a["y"]!.GetValue<double>(), a["reacting"]!.GetValue<bool>()));
        if (drawing?["bonds"] is JsonArray bonds)
            foreach (var b in bonds) _bonds.Add((b![0]!.GetValue<int>(), b[1]!.GetValue<int>()));
        InvalidateVisual();
    }

    private static (int, int) Key((int A, int B) e) => (Math.Min(e.A, e.B), Math.Max(e.A, e.B));

    private static Color ElementColour(string s) => s switch
    {
        "O" => Color.Parse("#E5534B"), "N" => Color.Parse("#6F9BF0"), "S" => Color.Parse("#E0B84C"), "Cl" => Color.Parse("#57B26A"),
        "F" => Color.Parse("#9ED36A"), "Si" => Color.Parse("#D6A45E"), _ => ((Tokens.Brush("TextB") as ISolidColorBrush)?.Color ?? Colors.White),
    };

    public override void Render(DrawingContext ctx)
    {
        var b = Bounds;
        if (_atoms.Count == 0 || b.Width < 40) return;
        double minx = _atoms.Min(a => a.X), maxx = _atoms.Max(a => a.X), miny = _atoms.Min(a => a.Y), maxy = _atoms.Max(a => a.Y);
        var span = Math.Max(Math.Max(maxx - minx, maxy - miny), 1.0);
        var scale = Math.Min((b.Width - 90) / Math.Max(maxx - minx, 1.0), (b.Height - 90) / Math.Max(maxy - miny, 1.0));
        scale = Math.Min(scale, 70);
        var cx = (minx + maxx) / 2;
        var cy = (miny + maxy) / 2;
        Point P(double x, double y) => new(b.Width / 2 + (x - cx) * scale, b.Height / 2 - (y - cy) * scale);
        var pos = _atoms.ToDictionary(a => a.Map, a => P(a.X, a.Y));
        var text = Tokens.Brush("TextB");
        foreach (var (a, c) in _bonds)
        {
            if (!pos.TryGetValue(a, out var p) || !pos.TryGetValue(c, out var q)) continue;
            var k = Key((a, c));
            var pen = _formed.Contains(k) ? new Pen(Tokens.Brush("OkB"), 2.6)
                : _broken.Contains(k) ? new Pen(Tokens.Brush("ErrB"), 2.2, new DashStyle([4, 3], 0)) : new Pen(text, 2.0);
            var d = q - p;
            var len = Math.Sqrt(d.X * d.X + d.Y * d.Y);
            if (len < 1) continue;
            var u = new Vector(d.X / len, d.Y / len);
            ctx.DrawLine(pen, p + u * 14, q - u * 14);
        }
        var face = new Typeface(Tokens.Sans, FontStyle.Normal, FontWeight.SemiBold);
        var mono = new Typeface(Tokens.Mono);
        foreach (var a in _atoms)
        {
            var p = pos[a.Map];
            if (a.Reacting) ctx.DrawEllipse(null, new Pen(Tokens.Brush("SelB"), 2), p, 16, 16);
            var ft = new FormattedText(a.Symbol, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, face, 17, new SolidColorBrush(ElementColour(a.Symbol)));
            ctx.DrawText(ft, new Point(p.X - ft.Width / 2, p.Y - ft.Height / 2));
            var badge = new Point(p.X + 16, p.Y - 15);
            ctx.DrawEllipse(Tokens.Brush("Bg3B"), new Pen(Tokens.Brush("AccB"), 1), badge, 9, 9);
            var n = new FormattedText(a.Map.ToString(CultureInfo.InvariantCulture), CultureInfo.InvariantCulture, FlowDirection.LeftToRight, mono, 10, Tokens.Brush("AccB"));
            ctx.DrawText(n, new Point(badge.X - n.Width / 2, badge.Y - n.Height / 2));
        }
    }
}
