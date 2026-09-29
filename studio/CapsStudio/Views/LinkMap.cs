using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>Which chains a React run linked, and where: each chain as a bar of its repeat units, each link a line between the
/// two units it joins, labelled with the molecule that bridges them (MAH #4) — the same links the live view highlights.</summary>
public sealed class LinkMap : Control
{
    private IReadOnlyList<MainViewModel.RxLink> _links = [];
    private int _chains, _units;
    // chain colours in the order chains are numbered
    private static readonly Color[] Palette =
    [
        Color.Parse("#E8A33D"), Color.Parse("#4FA3D8"), Color.Parse("#6CC08B"), Color.Parse("#D76B8E"), Color.Parse("#9B7AD5"),
        Color.Parse("#D8C24F"), Color.Parse("#5CC4C0"), Color.Parse("#E07A4F"), Color.Parse("#8FA0B8"), Color.Parse("#B8D86B"),
    ];

    public void SetLinks(IReadOnlyList<MainViewModel.RxLink> links, int chains, int units)
    {
        _links = links;
        _chains = chains;
        _units = Math.Max(units, links.Count == 0 ? 1 : (int)links.Max(l => Math.Max(l.UnitA, l.UnitB)));
        InvalidateVisual();
    }

    public override void Render(DrawingContext ctx)
    {
        var b = Bounds;
        var dim = Tokens.Brush("DimB");
        var tf = new Typeface(Tokens.Sans);
        void Text(string s, double x, double y, IBrush br, bool centre = false, double size = 10)
        {
            var ft = new FormattedText(s, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, tf, size, br);
            ctx.DrawText(ft, new Point(centre ? x - ft.Width / 2 : x, y - ft.Height / 2));
        }
        if (_chains <= 0)
        {
            Text("Crosslink to see which chains join, and where", b.Width / 2, b.Height / 2, dim, centre: true, size: 11);
            return;
        }
        const double left = 62, right = 16, top = 16, bottom = 22;
        var gap = _chains > 1 ? (b.Height - top - bottom) / (_chains - 1) : 0;
        double Y(int chain) => top + (Math.Clamp(chain, 1, _chains) - 1) * gap;
        double X(int unit) => left + (unit <= 0 ? 0.5 : (unit - 1) / Math.Max(1.0, _units - 1)) * (b.Width - left - right);
        IBrush ChainBrush(int c) => new SolidColorBrush(Palette[(Math.Max(1, c) - 1) % Palette.Length]);
        for (var c = 1; c <= _chains; ++c)
        {
            var y = Y(c);
            ctx.DrawLine(new Pen(ChainBrush(c), 5, lineCap: PenLineCap.Round), new Point(X(1), y), new Point(X(_units), y));
            Text($"chain {c}", 6, y, dim);
        }
        Text("unit 1", X(1), b.Height - 8, dim, centre: true);
        Text($"unit {_units}", X(_units), b.Height - 8, dim, centre: true);
        var acc = Tokens.Brush("AccB");
        var pen = new Pen(Tokens.Brush("TextB"), 1.6);
        foreach (var l in _links)
        {
            Point a = new(X(l.UnitA), Y(l.ChainA)), q = new(X(l.UnitB), Y(l.ChainB));
            if (l.ChainA == l.ChainB)
            {
                // two units of one chain (through a crosslinker bonded to the same chain twice would be a loop; shown as an arc)
                var g = new StreamGeometry();
                using (var gc = g.Open())
                {
                    gc.BeginFigure(a, false);
                    gc.QuadraticBezierTo(new Point((a.X + q.X) / 2, a.Y - 22), q);
                }
                ctx.DrawGeometry(null, pen, g);
            }
            else ctx.DrawLine(pen, a, q);
            ctx.DrawEllipse(acc, null, a, 4.5, 4.5);
            ctx.DrawEllipse(acc, null, q, 4.5, 4.5);
            var label = l.Via > 0 ? $"{l.ViaName} #{l.Via}" : "C–C";
            Text(label, (a.X + q.X) / 2 + 4, (a.Y + q.Y) / 2, Tokens.Brush("MutedB"), size: 9.5);
        }
    }
}
