using System.Globalization;
using System.Text;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Media;
using Avalonia.Media.Imaging;
using CapsStudio.ViewModels;

namespace CapsStudio.Views;

/// <summary>Figure overlays (design/boards/FigureBackground) drawn the same way on the preview tiles and in the exported
/// PNG; the SVG gets the same marks as vector elements.</summary>
public static class FigureDrawing
{
    private static readonly uint[] Viridis = [0x440154, 0x472D7B, 0x3B528B, 0x2C728E, 0x21918C, 0x28AE80, 0x5EC962, 0xADDC30, 0xFDE725];
    private static Color C(uint rgb) => Color.FromRgb((byte)(rgb >> 16), (byte)(rgb >> 8), (byte)rgb);
    private static FontFamily Sans => Tokens.Sans;
    private static FontFamily Mono => Tokens.Mono;

    private static FormattedText Text(string s, FontFamily f, FontWeight w, double size, IBrush b) =>
        new(s, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, new Typeface(f, FontStyle.Normal, w), size, b);

    /// <summary>Draws o in the image's own pixels (the caller scales the context to fit a smaller control).</summary>
    public static void Draw(DrawingContext ctx, FigureOverlay o)
    {
        var ink = new SolidColorBrush(C(o.Ink));
        var m = o.Margin;
        if (o.Title is { } title)
            ctx.DrawText(Text(title, Sans, FontWeight.SemiBold, o.TitlePx, ink), new Point(m, m * 0.8));
        if (o.BarPx > 0)
        {
            var y = o.Height - m;
            ctx.DrawLine(new Pen(ink, o.BarThickness, lineCap: PenLineCap.Flat), new Point(m, y), new Point(m + o.BarPx, y));
            var t = Text(o.BarLabel, Mono, FontWeight.Normal, o.LabelPx, ink);
            ctx.DrawText(t, new Point(m + o.BarPx / 2 - t.Width / 2, y - o.BarThickness - t.Height - o.LabelPx * 0.15));
        }
        if (o.Legend)
        {
            var w = o.Width * 0.2;
            var h = Math.Max(6, o.Width * 0.012);
            var x = o.Width - m - w;
            var y = o.Height - m - h - o.LabelPx * 1.4;
            var stops = new GradientStops();
            for (int k = 0; k < Viridis.Length; ++k) stops.Add(new GradientStop(C(Viridis[k]), k / (double)(Viridis.Length - 1)));
            ctx.DrawRectangle(new LinearGradientBrush { StartPoint = new RelativePoint(0, 0, RelativeUnit.Relative), EndPoint = new RelativePoint(1, 0, RelativeUnit.Relative), GradientStops = stops },
                              null, new Rect(x, y, w, h));
            var head = Text(o.LegendTitle, Sans, FontWeight.Normal, o.LabelPx, ink);
            ctx.DrawText(head, new Point(x, y - head.Height - o.LabelPx * 0.2));
            var lo = Text(o.LegendLo, Mono, FontWeight.Normal, o.LabelPx, ink);
            var hi = Text(o.LegendHi, Mono, FontWeight.Normal, o.LabelPx, ink);
            ctx.DrawText(lo, new Point(x, y + h + o.LabelPx * 0.2));
            ctx.DrawText(hi, new Point(x + w - hi.Width, y + h + o.LabelPx * 0.2));
        }
    }

    /// <summary>Render overlays (design/boards/RenderOverlays) in the image's own pixels: label box top left, vertical
    /// legend top right, scale bar bottom left, axis tripod bottom right.</summary>
    public static void DrawRender(DrawingContext ctx, RenderSpec o)
    {
        var u = o.Unit;
        var ink = new SolidColorBrush(C(o.Ink));
        var muted = new SolidColorBrush(o.Dark ? Color.Parse("#A5ABB1") : Color.Parse("#5A6168"));
        if (o.Lines.Length > 0)
        {
            var head = Text(o.Lines[0], Sans, FontWeight.SemiBold, 15 * u, ink);
            var sub = o.Lines.Length > 1 ? Text(o.Lines[1], Mono, FontWeight.Normal, 12 * u, muted) : null;
            var bw = Math.Max(head.Width, sub?.Width ?? 0) + 16 * u;
            var bh = head.Height + (sub?.Height ?? 0) + 12 * u;
            ctx.DrawRectangle(new SolidColorBrush(o.Dark ? Color.Parse("#0B0D0F") : Colors.White, 0.8), null, new Rect(10 * u, 12 * u, bw, bh), 6 * u, 6 * u);
            ctx.DrawText(head, new Point(18 * u, 17 * u));
            if (sub != null) ctx.DrawText(sub, new Point(18 * u, 17 * u + head.Height + 2 * u));
        }
        if (o.Legend)
        {
            var x = o.Width - 46 * u;
            var y = 24 * u;
            var stops = new GradientStops();
            for (int k = 0; k < Viridis.Length; ++k) stops.Add(new GradientStop(C(Viridis[Viridis.Length - 1 - k]), k / (double)(Viridis.Length - 1)));   // high at the top
            ctx.DrawRectangle(new LinearGradientBrush { StartPoint = new RelativePoint(0, 0, RelativeUnit.Relative), EndPoint = new RelativePoint(0, 1, RelativeUnit.Relative), GradientStops = stops },
                              new Pen(muted, Math.Max(1, u * 0.6)), new Rect(x, y, 14 * u, 160 * u));
            var hi = Text(o.LegendHi, Mono, FontWeight.Normal, 11 * u, ink);
            var lo = Text(o.LegendLo, Mono, FontWeight.Normal, 11 * u, ink);
            var name = Text(o.LegendName, Sans, FontWeight.Normal, 11 * u, muted);
            ctx.DrawText(hi, new Point(x - 6 * u - hi.Width, y));
            ctx.DrawText(lo, new Point(x - 6 * u - lo.Width, y + 160 * u - lo.Height));
            ctx.DrawText(name, new Point(x + 14 * u - name.Width, y + 166 * u));
        }
        if (o.BarPx > 0)
        {
            var y = o.Height - 24 * u;
            ctx.DrawLine(new Pen(ink, 3 * u, lineCap: PenLineCap.Flat), new Point(20 * u, y), new Point(20 * u + o.BarPx, y));
            var t = Text(o.BarLabel, Mono, FontWeight.Normal, 11 * u, ink);
            ctx.DrawText(t, new Point(20 * u + o.BarPx / 2 - t.Width / 2, y - 6 * u - t.Height));
        }
        if (o.Tripod)
        {
            // the renderer's rotation (yaw about y, then pitch about x) applied to the x, y and z unit vectors
            double cy = Math.Cos(o.Yaw), sy = Math.Sin(o.Yaw), cp = Math.Cos(o.Pitch), sp = Math.Sin(o.Pitch);
            var cx = o.Width - 40 * u;
            var cyy = o.Height - 40 * u;
            var len = 22 * u;
            var axes = new (string Name, double X, double Y, double Z, uint Col)[]
            {
                ("x", cy, sy * sp, -sy * cp, 0xE07A5F), ("y", 0, cp, sp, 0x7DC884), ("z", sy, -cy * sp, cy * cp, 0x5B8DEF),
            };
            foreach (var a in axes.OrderBy(a => a.Z))
            {
                var end = new Point(cx + a.X * len, cyy - a.Y * len);
                ctx.DrawLine(new Pen(new SolidColorBrush(C(a.Col)), 2 * u, lineCap: PenLineCap.Round), new Point(cx, cyy), end);
                var t = Text(a.Name, Mono, FontWeight.SemiBold, 10 * u, new SolidColorBrush(C(a.Col)));
                ctx.DrawText(t, new Point(end.X + (a.X >= 0 ? 2 * u : -2 * u - t.Width), end.Y - t.Height / 2 - (a.Y >= 0 ? 3 * u : -3 * u)));
            }
        }
    }

    /// <summary>A rendered image with the render overlays, written as PNG.</summary>
    public static void SaveRenderPng(byte[] rgba, int w, int h, RenderSpec o, string path)
    {
        using var img = MainViewModel.ToBitmap(rgba, w, h);
        using var rtb = new RenderTargetBitmap(new PixelSize(w, h), new Vector(96, 96));
        using (var ctx = rtb.CreateDrawingContext())
        {
            ctx.DrawImage(img, new Rect(0, 0, w, h));
            DrawRender(ctx, o);
        }
        rtb.Save(path);
    }

    /// <summary>The structure (straight-alpha RGBA) with the overlay, written as PNG with the dpi in a pHYs chunk.</summary>
    public static void SavePng(byte[] rgba, int w, int h, FigureOverlay o, double dpi, string path)
    {
        using var img = MainViewModel.ToBitmap(rgba, w, h);
        using (var rtb = new RenderTargetBitmap(new PixelSize(w, h), new Vector(96, 96)))
        {
            using (var ctx = rtb.CreateDrawingContext())
            {
                ctx.DrawImage(img, new Rect(0, 0, w, h));
                Draw(ctx, o);
            }
            rtb.Save(path);
        }
        if (dpi > 0) WriteDpi(path, dpi);
    }

    /// <summary>The overlay as SVG elements, placed before the closing tag of the structure's SVG.</summary>
    public static void AddToSvg(string path, FigureOverlay o)
    {
        var svg = File.ReadAllText(path);
        var inv = CultureInfo.InvariantCulture;
        string Hex(uint c) => "#" + c.ToString("X6", inv);
        var ink = Hex(o.Ink);
        var sb = new StringBuilder("<g id=\"caps-figure-overlay\">");
        var m = o.Margin;
        if (o.Title is { } title)
            sb.Append(string.Format(inv, "<text x=\"{0:F1}\" y=\"{1:F1}\" font-family=\"IBM Plex Sans, Helvetica, Arial, sans-serif\" font-weight=\"600\" font-size=\"{2:F1}\" fill=\"{3}\">{4}</text>",
                m, m * 0.8 + o.TitlePx, o.TitlePx, ink, System.Security.SecurityElement.Escape(title)));
        if (o.BarPx > 0)
        {
            var y = o.Height - m;
            sb.Append(string.Format(inv, "<line x1=\"{0:F1}\" y1=\"{1:F1}\" x2=\"{2:F1}\" y2=\"{1:F1}\" stroke=\"{3}\" stroke-width=\"{4:F1}\"/>", m, y, m + o.BarPx, ink, o.BarThickness));
            sb.Append(string.Format(inv, "<text x=\"{0:F1}\" y=\"{1:F1}\" text-anchor=\"middle\" font-family=\"IBM Plex Mono, Menlo, monospace\" font-size=\"{2:F1}\" fill=\"{3}\">{4}</text>",
                m + o.BarPx / 2, y - o.BarThickness - o.LabelPx * 0.4, o.LabelPx, ink, o.BarLabel));
        }
        if (o.Legend)
        {
            var w = o.Width * 0.2;
            var h = Math.Max(6, o.Width * 0.012);
            var x = o.Width - m - w;
            var y = o.Height - m - h - o.LabelPx * 1.4;
            sb.Append("<defs><linearGradient id=\"caps-legend\">");
            for (int k = 0; k < Viridis.Length; ++k) sb.Append(string.Format(inv, "<stop offset=\"{0:F3}\" stop-color=\"{1}\"/>", k / (double)(Viridis.Length - 1), Hex(Viridis[k])));
            sb.Append("</linearGradient></defs>");
            sb.Append(string.Format(inv, "<rect x=\"{0:F1}\" y=\"{1:F1}\" width=\"{2:F1}\" height=\"{3:F1}\" fill=\"url(#caps-legend)\"/>", x, y, w, h));
            sb.Append(string.Format(inv, "<text x=\"{0:F1}\" y=\"{1:F1}\" font-family=\"IBM Plex Sans, Helvetica, Arial, sans-serif\" font-size=\"{2:F1}\" fill=\"{3}\">{4}</text>", x, y - o.LabelPx * 0.4, o.LabelPx, ink, o.LegendTitle));
            sb.Append(string.Format(inv, "<text x=\"{0:F1}\" y=\"{1:F1}\" font-family=\"IBM Plex Mono, Menlo, monospace\" font-size=\"{2:F1}\" fill=\"{3}\">{4}</text>", x, y + h + o.LabelPx * 1.2, o.LabelPx, ink, o.LegendLo));
            sb.Append(string.Format(inv, "<text x=\"{0:F1}\" y=\"{1:F1}\" text-anchor=\"end\" font-family=\"IBM Plex Mono, Menlo, monospace\" font-size=\"{2:F1}\" fill=\"{3}\">{4}</text>", x + w, y + h + o.LabelPx * 1.2, o.LabelPx, ink, o.LegendHi));
        }
        sb.Append("</g>");
        var end = svg.LastIndexOf("</svg>", StringComparison.Ordinal);
        if (end < 0) return;
        File.WriteAllText(path, svg[..end] + sb + svg[end..]);
    }

    /// <summary>Inserts (or replaces) the pHYs chunk so the PNG carries its dpi (pixels per metre).</summary>
    public static void WriteDpi(string path, double dpi)
    {
        var png = File.ReadAllBytes(path);
        if (png.Length < 33 || png[12] != 'I' || png[13] != 'H') return;
        var ppm = (uint)Math.Round(dpi / 0.0254);
        var data = new byte[9];
        void Be(byte[] b, int at, uint v) { b[at] = (byte)(v >> 24); b[at + 1] = (byte)(v >> 16); b[at + 2] = (byte)(v >> 8); b[at + 3] = (byte)v; }
        Be(data, 0, ppm); Be(data, 4, ppm); data[8] = 1;
        var chunk = new byte[21];
        Be(chunk, 0, 9);
        Encoding.ASCII.GetBytes("pHYs").CopyTo(chunk, 4);
        data.CopyTo(chunk, 8);
        Be(chunk, 17, Crc(chunk.AsSpan(4, 13)));
        // drop an existing pHYs, then put ours right after IHDR (8-byte signature + 25-byte IHDR chunk)
        var outp = new List<byte>(png.Length + 21);
        outp.AddRange(png.AsSpan(0, 33).ToArray());
        outp.AddRange(chunk);
        int p = 33;
        while (p + 8 <= png.Length)
        {
            var len = (png[p] << 24) | (png[p + 1] << 16) | (png[p + 2] << 8) | png[p + 3];
            var type = Encoding.ASCII.GetString(png, p + 4, 4);
            var total = len + 12;
            if (p + total > png.Length) { outp.AddRange(png.AsSpan(p).ToArray()); break; }
            if (type != "pHYs") outp.AddRange(png.AsSpan(p, total).ToArray());
            p += total;
        }
        File.WriteAllBytes(path, outp.ToArray());
    }

    private static uint Crc(ReadOnlySpan<byte> b)
    {
        uint c = 0xFFFFFFFF;
        foreach (var x in b)
        {
            c ^= x;
            for (int k = 0; k < 8; ++k) c = (c & 1) != 0 ? 0xEDB88320 ^ (c >> 1) : c >> 1;
        }
        return c ^ 0xFFFFFFFF;
    }
}

/// <summary>The overlay of a preview tile, scaled from the tile's pixels to the control.</summary>
public sealed class FigureOverlayView : Control
{
    public static readonly StyledProperty<FigureOverlay?> OverlayProperty = AvaloniaProperty.Register<FigureOverlayView, FigureOverlay?>(nameof(Overlay));
    static FigureOverlayView() => AffectsRender<FigureOverlayView>(OverlayProperty);
    public FigureOverlay? Overlay { get => GetValue(OverlayProperty); set => SetValue(OverlayProperty, value); }

    public override void Render(DrawingContext ctx)
    {
        if (Overlay is not { } o || o.Width <= 0) return;
        var k = Math.Min(Bounds.Width / o.Width, Bounds.Height / o.Height);
        using (ctx.PushTransform(Matrix.CreateScale(k, k)))
            FigureDrawing.Draw(ctx, o);
    }
}

/// <summary>The render-frame guide over the 3D view: outside dimmed, the frame dashed, the overlays drawn inside.</summary>
public sealed class RenderGuideView : Control
{
    public MainViewModel? Vm { get; set; }

    public override void Render(DrawingContext ctx)
    {
        if (Vm?.RenderGuide(Bounds.Width, Bounds.Height) is not { } g) return;
        var dim = new SolidColorBrush(Colors.Black, 0.45);
        var (x, y, w, h) = (g.X, g.Y, g.W, g.H);
        ctx.FillRectangle(dim, new Rect(0, 0, Bounds.Width, Math.Max(0, y)));
        ctx.FillRectangle(dim, new Rect(0, y + h, Bounds.Width, Math.Max(0, Bounds.Height - y - h)));
        ctx.FillRectangle(dim, new Rect(0, y, Math.Max(0, x), h));
        ctx.FillRectangle(dim, new Rect(x + w, y, Math.Max(0, Bounds.Width - x - w), h));
        var acc = new SolidColorBrush(Color.Parse("#F5A524"));
        ctx.DrawRectangle(null, new Pen(acc, 1, new DashStyle([6, 4], 0)), new Rect(x, y, w, h));
        var cap = new FormattedText(Vm.RenderSizeText, CultureInfo.InvariantCulture, FlowDirection.LeftToRight, new Typeface(Tokens.Mono), 11, acc);
        ctx.DrawText(cap, new Point(x, Math.Max(2, y - cap.Height - 4)));
        var k = w / g.Spec.Width;
        using (ctx.PushClip(new Rect(x, y, w, h)))
        using (ctx.PushTransform(Matrix.CreateScale(k, k) * Matrix.CreateTranslation(x, y)))
            FigureDrawing.DrawRender(ctx, g.Spec);
    }
}

/// <summary>A checkerboard marking transparent pixels.</summary>
public sealed class CheckerView : Control
{
    public override void Render(DrawingContext ctx)
    {
        const double s = 10;
        ctx.FillRectangle(Brushes.White, new Rect(Bounds.Size));
        var grey = new SolidColorBrush(Color.Parse("#D9DCDF"));
        for (double y = 0; y < Bounds.Height; y += s)
            for (double x = ((int)(y / s) % 2) * s; x < Bounds.Width; x += 2 * s)
                ctx.FillRectangle(grey, new Rect(x, y, Math.Min(s, Bounds.Width - x), Math.Min(s, Bounds.Height - y)));
    }
}
