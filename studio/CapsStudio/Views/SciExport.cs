using System.Globalization;
using System.IO.Compression;
using System.Text;
using Avalonia;
using Avalonia.Controls;
using Avalonia.Input;
using Avalonia.Input.Platform;
using Avalonia.Media;
using Avalonia.Media.Imaging;
using Avalonia.Platform.Storage;

namespace CapsStudio.Views;

/// <summary>A figure written for print, the way caps.plotstyle writes the Python figures: one journal column (3.35 in)
/// or the page width (7 in), black on white, 8 pt Arial, inward ticks on four sides, Okabe and Ito's colours, as vector
/// SVG, vector PDF (the fonts embedded as subsets) or a PNG at 300 dpi. Every plot in the Studio has Save figure and
/// Copy figure on its context menu; the figure composer draws its plot panels the same way.</summary>
public static class SciExport
{
    public const double Single = 3.35, Double = 7.0;   // in
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    internal static string F(double v) => Math.Round(v, 3).ToString("0.###", Inv);

    /// <summary>The height for a width (in) from the shape on screen, kept between half and nine tenths of the width.</summary>
    public static double HeightFor(double widthIn, Size onScreen) =>
        widthIn * (onScreen.Width > 0 && onScreen.Height > 0 ? Math.Clamp(onScreen.Height / onScreen.Width, 0.5, 0.9) : 0.75);

    /// <summary>The figure as an SVG document widthIn × heightIn inches (coordinates in points).</summary>
    public static string ToSvg(SciFigure fig, double widthIn, double heightIn)
    {
        var s = new SvgSink(widthIn * 72, heightIn * 72);
        SciPlot.DrawTo(s, fig, new Size(widthIn * 72, heightIn * 72), SciStyle.Paper);
        return s.Document();
    }

    /// <summary>The figure as a one-page vector PDF widthIn × heightIn inches.</summary>
    public static byte[] ToPdf(SciFigure fig, double widthIn, double heightIn)
    {
        var p = new PdfSink(widthIn * 72, heightIn * 72);
        SciPlot.DrawTo(p, fig, new Size(widthIn * 72, heightIn * 72), SciStyle.Paper);
        return p.Document();
    }

    /// <summary>The figure drawn for print into an Avalonia context at scale (pixels per point), e.g. a PNG at dpi / 72.</summary>
    public static void Draw(DrawingContext ctx, SciFigure fig, Rect area, double scale)
    {
        using var t = ctx.PushTransform(Matrix.CreateScale(scale, scale) * Matrix.CreateTranslation(area.X, area.Y));
        SciPlot.DrawTo(new ScreenSink(ctx, SciStyle.Paper.Font), fig, new Size(area.Width / scale, area.Height / scale), SciStyle.Paper);
    }

    /// <summary>A PNG of the figure at dpi (white background, the dpi written into the file).</summary>
    public static void SavePng(SciFigure fig, string path, double widthIn, double heightIn, int dpi = 300)
    {
        var (w, h) = ((int)Math.Round(widthIn * dpi), (int)Math.Round(heightIn * dpi));
        using var rtb = new RenderTargetBitmap(new PixelSize(w, h), new Vector(96, 96));
        using (var ctx = rtb.CreateDrawingContext())
        {
            ctx.FillRectangle(Brushes.White, new Rect(0, 0, w, h));
            Draw(ctx, fig, new Rect(0, 0, w, h), dpi / 72.0);
        }
        rtb.Save(path);
        FigureFiles.SetPngDpi(path, dpi);
    }

    /// <summary>Writes the figure by the file's extension: .pdf, .svg or .png (300 dpi).</summary>
    public static string Save(SciFigure fig, string path, double widthIn, double heightIn)
    {
        var ext = Path.GetExtension(path).ToLowerInvariant();
        if (ext == ".svg") File.WriteAllText(path, ToSvg(fig, widthIn, heightIn));
        else if (ext == ".png") SavePng(fig, path, widthIn, heightIn);
        else File.WriteAllBytes(path, ToPdf(fig, widthIn, heightIn));
        return $"Wrote {Path.GetFileName(path)} · {F(widthIn)} × {F(heightIn)} in" + (ext == ".png" ? " at 300 dpi" : ext == ".svg" ? " · vector" : " · vector, fonts embedded");
    }

    /// <summary>The data behind the figure as CSV: one block per series (name, then x, y and ± when there are errors).</summary>
    public static string ToCsv(SciFigure fig)
    {
        var sb = new StringBuilder();
        string Q(string s) => "\"" + SciPlot.Plain(s).Replace("\"", "\"\"") + "\"";
        foreach (var s in fig.Series.Where(s => s.X.Length > 0))
        {
            sb.Append("# ").Append(s.Name.Length > 0 ? SciPlot.Plain(s.Name) : "series").Append('\n');
            sb.Append(Q(fig.XLabel)).Append(',').Append(Q(s.Right ? fig.Y2Label : fig.YLabel)).Append(s.Err != null ? ",\"±\"\n" : "\n");
            for (var i = 0; i < Math.Min(s.X.Length, s.Y.Length); ++i)
            {
                sb.Append(s.X[i].ToString("R", Inv)).Append(',').Append(s.Y[i].ToString("R", Inv));
                if (s.Err != null) sb.Append(',').Append(i < s.Err.Length ? s.Err[i].ToString("R", Inv) : "");
                sb.Append('\n');
            }
            sb.Append('\n');
        }
        return sb.ToString();
    }

    // ------------------------------------------------------------------ the context menu of every plot
    internal static void AttachMenu(SciPlot plot)
    {
        var menu = new ContextMenu();
        MenuItem Item(string header, Func<Task> run) { var m = new MenuItem { Header = header }; m.Click += async (_, _) => { try { await run(); } catch (Exception ex) { Status(plot, "Could not write the figure: " + ex.Message); } }; return m; }
        menu.Items.Add(Item("Save figure, one column (3.35 in)…", () => SaveAs(plot, Single)));
        menu.Items.Add(Item("Save figure, two columns (7 in)…", () => SaveAs(plot, Double)));
        menu.Items.Add(new Separator());
        menu.Items.Add(Item("Copy figure (SVG)", () => CopySvg(plot)));
        menu.Items.Add(Item("Copy data (CSV)", () => CopyCsv(plot)));
        menu.Opening += (_, e) => { if (plot.Figure == null) e.Cancel = true; };
        plot.ContextMenu = menu;
    }

    private static void Status(Control c, string text)
    {
        if (TopLevel.GetTopLevel(c)?.DataContext is ViewModels.MainViewModel vm) vm.Status = text;
    }

    private static async Task SaveAs(SciPlot plot, double widthIn)
    {
        if (plot.Figure is not { } fig || TopLevel.GetTopLevel(plot) is not { } top) return;
        var name = SciPlot.Plain(fig.YLabel).Split('/')[0].Trim();
        name = new string(name.Select(ch => char.IsLetterOrDigit(ch) ? ch : '_').ToArray()).Trim('_');
        var file = await top.StorageProvider.SaveFilePickerAsync(new FilePickerSaveOptions
        {
            Title = $"Save figure ({F(widthIn)} in wide)", SuggestedFileName = (name.Length > 0 ? name : "figure") + ".pdf", DefaultExtension = "pdf",
            FileTypeChoices =
            [
                new FilePickerFileType("PDF · vector") { Patterns = ["*.pdf"] },
                new FilePickerFileType("SVG · vector") { Patterns = ["*.svg"] },
                new FilePickerFileType("PNG · 300 dpi") { Patterns = ["*.png"] },
            ],
        });
        if (file?.TryGetLocalPath() is not { } path) return;
        Status(plot, Save(fig, path, widthIn, HeightFor(widthIn, plot.Bounds.Size)));
    }

    private static async Task CopySvg(SciPlot plot)
    {
        if (plot.Figure is not { } fig || TopLevel.GetTopLevel(plot)?.Clipboard is not { } cb) return;
        var svg = ToSvg(fig, Single, HeightFor(Single, plot.Bounds.Size));
        var data = new DataObject();
        data.Set(DataFormats.Text, svg);
        data.Set("image/svg+xml", Encoding.UTF8.GetBytes(svg));
        await cb.SetDataObjectAsync(data);
        Status(plot, "Copied the figure as SVG (one column, 3.35 in)");
    }

    private static async Task CopyCsv(SciPlot plot)
    {
        if (plot.Figure is not { } fig || TopLevel.GetTopLevel(plot)?.Clipboard is not { } cb) return;
        await cb.SetTextAsync(ToCsv(fig));
        Status(plot, "Copied the figure's data as CSV");
    }
}

/// <summary>Writes the drawing as SVG elements (points as user units, the page size in inches).</summary>
internal sealed class SvgSink(double widthPt, double heightPt) : ISciSink
{
    private readonly StringBuilder _b = new();
    private readonly StringBuilder _defs = new();
    private int _clips, _depth;
    private static string F(double v) => SciExport.F(v);
    private static string Hex(Color c) => $"#{c.R:X2}{c.G:X2}{c.B:X2}";
    private static string Esc(string s) => s.Replace("&", "&amp;").Replace("<", "&lt;").Replace(">", "&gt;").Replace("\"", "&quot;");
    private static string Dash(double[]? d, double w) => d == null ? "" : $" stroke-dasharray=\"{string.Join(" ", d.Select(v => F(v * w)))}\"";
    private void Indent() => _b.Append(' ', 2 * _depth);
    private static string Fill(Color? c, double a) => c is { } f ? $" fill=\"{Hex(f)}\"" + (a < 1 ? $" fill-opacity=\"{F(a)}\"" : "") : " fill=\"none\"";
    private static string Stroke(Color? c, double w) => c is { } s ? $" stroke=\"{Hex(s)}\" stroke-width=\"{F(w)}\"" : "";

    public void Line(Point a, Point b, Color c, double alpha, double w, double[]? dash = null, bool round = false)
    {
        Indent();
        _b.Append($"<line x1=\"{F(a.X)}\" y1=\"{F(a.Y)}\" x2=\"{F(b.X)}\" y2=\"{F(b.Y)}\"{Stroke(c, w)}{(alpha < 1 ? $" stroke-opacity=\"{F(alpha)}\"" : "")}{Dash(dash, w)}{(round ? " stroke-linecap=\"round\"" : "")}/>\n");
    }
    public void Rect(Rect r, Color? fill, double fillAlpha, Color? stroke, double w, double[]? dash = null)
    {
        Indent();
        _b.Append($"<rect x=\"{F(r.X)}\" y=\"{F(r.Y)}\" width=\"{F(r.Width)}\" height=\"{F(r.Height)}\"{Fill(fill, fillAlpha)}{Stroke(stroke, w)}{Dash(dash, w)}/>\n");
    }
    public void Path(IReadOnlyList<IReadOnlyList<Point>> pieces, bool closed, Color? fill, double fillAlpha, Color? stroke, double w, double[]? dash, bool round)
    {
        var d = new StringBuilder();
        foreach (var p in pieces.Where(p => p.Count > 1))
        {
            d.Append('M').Append(F(p[0].X)).Append(' ').Append(F(p[0].Y));
            for (var i = 1; i < p.Count; ++i) d.Append('L').Append(F(p[i].X)).Append(' ').Append(F(p[i].Y));
            if (closed) d.Append('Z');
        }
        if (d.Length == 0) return;
        Indent();
        _b.Append($"<path d=\"{d}\"{Fill(fill, fillAlpha)}{Stroke(stroke, w)}{Dash(dash, w)} stroke-linejoin=\"round\"{(round ? " stroke-linecap=\"round\"" : "")}/>\n");
    }
    public void Circle(Point c, double r, Color? fill, Color? stroke, double w)
    {
        Indent();
        _b.Append($"<circle cx=\"{F(c.X)}\" cy=\"{F(c.Y)}\" r=\"{F(r)}\"{Fill(fill, 1)}{Stroke(stroke, w)}/>\n");
    }
    public void Text(string s, Point topLeft, double size, bool italic, bool bold, Color c)
    {
        Indent();
        _b.Append($"<text x=\"{F(topLeft.X)}\" y=\"{F(topLeft.Y + Metrics(size).Ascent)}\" font-size=\"{F(size)}\"{(italic ? " font-style=\"italic\"" : "")}{(bold ? " font-weight=\"bold\"" : "")}{(c != Colors.Black ? $" fill=\"{Hex(c)}\"" : "")} xml:space=\"preserve\">{Esc(s)}</text>\n");
    }
    public double Measure(string s, double size, bool italic, bool bold) => PrintFonts.Available ? PrintFonts.Measure(s, size, italic, bold) : 0.56 * size * s.Length;
    public (double Height, double Ascent) Metrics(double size) => (PrintFonts.LineEm * size, PrintFonts.AscentEm * size);
    public void PushClip(Rect r)
    {
        var id = $"c{++_clips}";
        _defs.Append($"  <clipPath id=\"{id}\"><rect x=\"{F(r.X)}\" y=\"{F(r.Y)}\" width=\"{F(r.Width)}\" height=\"{F(r.Height)}\"/></clipPath>\n");
        Indent();
        _b.Append($"<g clip-path=\"url(#{id})\">\n");
        _depth++;
    }
    public void PushTransform(Matrix m)
    {
        Indent();
        _b.Append($"<g transform=\"matrix({F(m.M11)} {F(m.M12)} {F(m.M21)} {F(m.M22)} {F(m.M31)} {F(m.M32)})\">\n");
        _depth++;
    }
    public void Pop() { _depth--; Indent(); _b.Append("</g>\n"); }

    public string Document() =>
        $"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{F(widthPt / 72)}in\" height=\"{F(heightPt / 72)}in\" viewBox=\"0 0 {F(widthPt)} {F(heightPt)}\" " +
        $"font-family=\"{PrintFonts.SvgFamily}\">\n<defs>\n{_defs}</defs>\n<rect width=\"{F(widthPt)}\" height=\"{F(heightPt)}\" fill=\"#FFFFFF\"/>\n{_b}</svg>\n";

    /// <summary>The elements drawn so far, with their clip paths, for a page that holds several figures (the composer).</summary>
    internal (string Defs, string Body) Parts() => (_defs.ToString(), _b.ToString());
    internal void Raw(string s) => _b.Append(s);
}

/// <summary>Writes the drawing as one PDF page in points, y down from the top left (a flip on the page's matrix); text is
/// set in subsets of the print fonts embedded as CID TrueType fonts (Identity-H) with a ToUnicode map, so it stays text,
/// selectable and searchable; transparency through graphics states.</summary>
internal sealed class PdfSink : ISciSink
{
    private readonly double _w, _h;
    private readonly StringBuilder _c = new();
    private readonly Dictionary<TrueTypeFont, (string Res, HashSet<int> Gids, Dictionary<int, int> Uni)> _fonts = new();
    private readonly Dictionary<double, string> _alphas = new();
    private readonly List<(string Res, byte[] Rgb, int W, int H)> _images = new();
    private static string F(double v) => SciExport.F(v);
    private static string Rgb(Color c) => $"{F(c.R / 255.0)} {F(c.G / 255.0)} {F(c.B / 255.0)}";

    public PdfSink(double widthPt, double heightPt)
    {
        _w = widthPt; _h = heightPt;
        _c.Append($"1 0 0 -1 0 {F(heightPt)} cm\n1 1 1 rg 0 0 {F(widthPt)} {F(heightPt)} re f\n");
    }

    private string Alpha(double a)
    {
        a = Math.Round(Math.Clamp(a, 0, 1), 3);
        if (!_alphas.TryGetValue(a, out var r)) _alphas[a] = r = $"GS{_alphas.Count}";
        return r;
    }
    private void Stroke(Color c, double alpha, double w, double[]? dash, bool round)
    {
        _c.Append($"{Rgb(c)} RG {F(w)} w {(round ? 1 : 0)} J 1 j [{(dash == null ? "" : string.Join(" ", dash.Select(v => F(v * w))))}] 0 d");
        if (alpha < 1) _c.Append($" /{Alpha(alpha)} gs");
        _c.Append('\n');
    }

    public void Line(Point a, Point b, Color c, double alpha, double w, double[]? dash = null, bool round = false)
    {
        _c.Append("q ");
        Stroke(c, alpha, w, dash, round);
        _c.Append($"{F(a.X)} {F(a.Y)} m {F(b.X)} {F(b.Y)} l S Q\n");
    }
    public void Rect(Rect r, Color? fill, double fillAlpha, Color? stroke, double w, double[]? dash = null)
    {
        var rect = $"{F(r.X)} {F(r.Y)} {F(r.Width)} {F(r.Height)} re";
        if (fill is { } f) { _c.Append($"q {Rgb(f)} rg{(fillAlpha < 1 ? $" /{Alpha(fillAlpha)} gs" : "")} {rect} f Q\n"); }
        if (stroke is { } s) { _c.Append("q "); Stroke(s, 1, w, dash, false); _c.Append($"{rect} S Q\n"); }
    }
    public void Path(IReadOnlyList<IReadOnlyList<Point>> pieces, bool closed, Color? fill, double fillAlpha, Color? stroke, double w, double[]? dash, bool round)
    {
        var d = new StringBuilder();
        foreach (var p in pieces.Where(p => p.Count > 1))
        {
            d.Append($"{F(p[0].X)} {F(p[0].Y)} m\n");
            for (var i = 1; i < p.Count; ++i) d.Append($"{F(p[i].X)} {F(p[i].Y)} l\n");
            if (closed) d.Append("h\n");
        }
        if (d.Length == 0) return;
        if (fill is { } f) _c.Append($"q {Rgb(f)} rg{(fillAlpha < 1 ? $" /{Alpha(fillAlpha)} gs" : "")}\n{d}f Q\n");
        if (stroke is { } s) { _c.Append("q "); Stroke(s, 1, w, dash, round); _c.Append(d).Append("S Q\n"); }
    }
    public void Circle(Point c, double r, Color? fill, Color? stroke, double w)
    {
        // four Bézier quarters
        const double k = 0.5523;
        var d = $"{F(c.X + r)} {F(c.Y)} m {F(c.X + r)} {F(c.Y + k * r)} {F(c.X + k * r)} {F(c.Y + r)} {F(c.X)} {F(c.Y + r)} c " +
                $"{F(c.X - k * r)} {F(c.Y + r)} {F(c.X - r)} {F(c.Y + k * r)} {F(c.X - r)} {F(c.Y)} c " +
                $"{F(c.X - r)} {F(c.Y - k * r)} {F(c.X - k * r)} {F(c.Y - r)} {F(c.X)} {F(c.Y - r)} c " +
                $"{F(c.X + k * r)} {F(c.Y - r)} {F(c.X + r)} {F(c.Y - k * r)} {F(c.X + r)} {F(c.Y)} c h";
        _c.Append("q ");
        if (fill is { } f) _c.Append($"{Rgb(f)} rg ");
        if (stroke is { } s) Stroke(s, 1, w, null, false);
        _c.Append(d).Append(fill != null && stroke != null ? " B" : fill != null ? " f" : " S").Append(" Q\n");
    }
    public void Text(string s, Point topLeft, double size, bool italic, bool bold, Color c)
    {
        var x = topLeft.X;
        var y = topLeft.Y + Metrics(size).Ascent;
        foreach (var (font, synth, text) in PrintFonts.Runs(s, italic, bold))
        {
            if (!_fonts.TryGetValue(font, out var e)) _fonts[font] = e = ($"F{_fonts.Count + 1}", new HashSet<int>(), new Dictionary<int, int>());
            var hex = new StringBuilder();
            double adv = 0;
            foreach (var cp in PrintFonts.Codepoints(text))
            {
                var gid = font.Glyph(cp);
                e.Gids.Add(gid);
                if (gid != 0) e.Uni.TryAdd(gid, cp);
                hex.Append(gid.ToString("X4"));
                adv += font.Advance(gid) * size;
            }
            // the page is flipped (y down): the text matrix flips the glyphs back, and slants them for a synthetic italic
            _c.Append($"BT {Rgb(c)} rg /{e.Res} {F(size)} Tf 1 0 {(synth ? "0.21" : "0")} -1 {F(x)} {F(y)} Tm <{hex}> Tj ET\n");
            x += adv;
        }
    }
    public double Measure(string s, double size, bool italic, bool bold) => PrintFonts.Measure(s, size, italic, bold);
    public (double Height, double Ascent) Metrics(double size) => (PrintFonts.LineEm * size, PrintFonts.AscentEm * size);
    public void PushClip(Rect r) => _c.Append($"q {F(r.X)} {F(r.Y)} {F(r.Width)} {F(r.Height)} re W n\n");
    public void PushTransform(Matrix m) => _c.Append($"q {F(m.M11)} {F(m.M12)} {F(m.M21)} {F(m.M22)} {F(m.M31)} {F(m.M32)} cm\n");
    public void Pop() => _c.Append("Q\n");

    /// <summary>An RGBA image (top row first) drawn into r (the 3D panel of a composed page).</summary>
    internal void Image(byte[] rgba, int w, int h, Rect r)
    {
        var rgb = new byte[w * h * 3];
        for (var i = 0; i < w * h; ++i) { rgb[3 * i] = rgba[4 * i]; rgb[3 * i + 1] = rgba[4 * i + 1]; rgb[3 * i + 2] = rgba[4 * i + 2]; }
        var res = $"Im{_images.Count}";
        _images.Add((res, rgb, w, h));
        // image space is the unit square with its first row at the top in a y-up frame: flip it back on the y-down page
        _c.Append($"q {F(r.Width)} 0 0 {F(-r.Height)} {F(r.X)} {F(r.Y + r.Height)} cm /{res} Do Q\n");
    }

    private static byte[] Deflate(byte[] b)
    {
        using var ms = new MemoryStream();
        using (var z = new ZLibStream(ms, CompressionLevel.Optimal, true)) z.Write(b, 0, b.Length);
        return ms.ToArray();
    }

    public byte[] Document()
    {
        var objs = new List<byte[]>();
        int Add(string s) { objs.Add(Encoding.ASCII.GetBytes(s)); return objs.Count; }
        int AddStream(string dict, byte[] data) { var head = Encoding.ASCII.GetBytes($"{dict.TrimEnd('>', ' ')} /Length {data.Length} >>\nstream\n"); var tail = Encoding.ASCII.GetBytes("\nendstream"); objs.Add([.. head, .. data, .. tail]); return objs.Count; }
        var fontRefs = new StringBuilder();
        var tag = 0;
        foreach (var (font, (res, gids, uni)) in _fonts)
        {
            var sub = font.Subset(gids);
            var file = AddStream($"<< /Filter /FlateDecode /Length1 {sub.Length} >>", Deflate(sub));
            var name = new string(Enumerable.Range(0, 6).Select(k => (char)('A' + (tag * 7 + k * 3 + font.Name.Length) % 26)).ToArray()) + "+" + font.Name;
            tag++;
            var u = font.UnitsPerEm;
            var desc = Add($"<< /Type /FontDescriptor /FontName /{name} /Flags {(font.Italic ? 96 : 32)} /FontBBox [{font.BBox.X0 * 1000 / u} {font.BBox.Y0 * 1000 / u} {font.BBox.X1 * 1000 / u} {font.BBox.Y1 * 1000 / u}] " +
                           $"/ItalicAngle {(font.Italic ? -12 : 0)} /Ascent {F(font.Ascent * 1000)} /Descent {F(font.Descent * 1000)} /CapHeight {F(font.CapHeight * 1000)} /StemV 80 /FontFile2 {file} 0 R >>");
            var widths = string.Join(" ", gids.OrderBy(g => g).Select(g => $"{g} [{F(font.Advance(g) * 1000)}]"));
            var cid = Add($"<< /Type /Font /Subtype /CIDFontType2 /BaseFont /{name} /CIDSystemInfo << /Registry (Adobe) /Ordering (Identity) /Supplement 0 >> /FontDescriptor {desc} 0 R /CIDToGIDMap /Identity /DW 1000 /W [{widths}] >>");
            var cmap = new StringBuilder("/CIDInit /ProcSet findresource begin 12 dict begin begincmap /CIDSystemInfo << /Registry (Adobe) /Ordering (UCS) /Supplement 0 >> def /CMapName /Adobe-Identity-UCS def /CMapType 2 def\n1 begincodespacerange <0000> <FFFF> endcodespacerange\n");
            foreach (var chunk in uni.OrderBy(p => p.Key).Chunk(100))
            {
                cmap.Append($"{chunk.Length} beginbfchar\n");
                foreach (var (g, cp) in chunk) cmap.Append($"<{g:X4}> <{string.Concat(Encoding.BigEndianUnicode.GetBytes(char.ConvertFromUtf32(cp)).Select(b => b.ToString("X2")))}>\n");
                cmap.Append("endbfchar\n");
            }
            cmap.Append("endcmap CMapName currentdict /CMap defineresource pop end end\n");
            var tu = AddStream("<< >>", Encoding.ASCII.GetBytes(cmap.ToString()));
            var type0 = Add($"<< /Type /Font /Subtype /Type0 /BaseFont /{name} /Encoding /Identity-H /DescendantFonts [{cid} 0 R] /ToUnicode {tu} 0 R >>");
            fontRefs.Append($"/{res} {type0} 0 R ");
        }
        var gs = string.Concat(_alphas.Select(a => $"/{a.Value} << /ca {F(a.Key)} /CA {F(a.Key)} >> "));
        var imgs = new StringBuilder();
        foreach (var (res, rgb, w, h) in _images)
            imgs.Append($"/{res} {AddStream($"<< /Type /XObject /Subtype /Image /Width {w} /Height {h} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode >>", Deflate(rgb))} 0 R ");
        var content = AddStream("<< /Filter /FlateDecode >>", Deflate(Encoding.UTF8.GetBytes(_c.ToString())));
        var pagesId = objs.Count + 2;
        var page = Add($"<< /Type /Page /Parent {pagesId} 0 R /MediaBox [0 0 {F(_w)} {F(_h)}] /Resources << /Font << {fontRefs}>> /ExtGState << {gs}>> /XObject << {imgs}>> >> /Contents {content} 0 R >>");
        Add($"<< /Type /Pages /Kids [{page} 0 R] /Count 1 >>");
        var catalog = Add($"<< /Type /Catalog /Pages {pagesId} 0 R >>");
        var info = Add($"<< /Producer (CAPS Studio) /CreationDate (D:{DateTime.UtcNow:yyyyMMddHHmmss}Z) >>");
        using var o = new MemoryStream();
        void W(string s) { var b = Encoding.ASCII.GetBytes(s); o.Write(b, 0, b.Length); }
        W("%PDF-1.6\n");
        o.Write([(byte)'%', 0xE2, 0xE3, 0xCF, 0xD3, (byte)'\n']);   // binary content follows
        var offsets = new long[objs.Count + 1];
        for (var i = 0; i < objs.Count; ++i)
        {
            offsets[i + 1] = o.Position;
            W($"{i + 1} 0 obj\n");
            o.Write(objs[i]);
            W("\nendobj\n");
        }
        var xref = o.Position;
        W($"xref\n0 {objs.Count + 1}\n0000000000 65535 f \n");
        for (var i = 1; i <= objs.Count; ++i) W($"{offsets[i]:D10} 00000 n \n");
        W($"trailer\n<< /Size {objs.Count + 1} /Root {catalog} 0 R /Info {info} 0 R >>\nstartxref\n{xref}\n%%EOF\n");
        return o.ToArray();
    }
}
