using System.Text;

namespace CapsStudio.Views;

/// <summary>A TrueType font file read far enough to lay out text (character → glyph, advance widths, ascent) and to
/// embed a subset of it in a PDF (the glyphs used, ids kept, the rest emptied). Only fonts with TrueType outlines (glyf)
/// are read; a missing or unreadable file gives null.</summary>
internal sealed class TrueTypeFont
{
    public string Path { get; }
    public string Name { get; }
    public int UnitsPerEm { get; }
    public int NumGlyphs { get; }
    public double Ascent { get; }       // em
    public double Descent { get; }      // em, negative
    public double CapHeight { get; }    // em
    public (int X0, int Y0, int X1, int Y1) BBox { get; }
    public bool Italic { get; }
    private readonly byte[] _d;
    private readonly Dictionary<string, (int Off, int Len)> _tables = new();
    private readonly ushort[] _adv;
    private readonly Dictionary<int, int> _cmap = new();
    private readonly uint[] _loca;

    private TrueTypeFont(string path, byte[] d)
    {
        Path = path;
        _d = d;
        int n = U16(4);
        for (var i = 0; i < n; ++i)
        {
            var r = 12 + 16 * i;
            _tables[Encoding.ASCII.GetString(d, r, 4)] = ((int)U32(r + 8), (int)U32(r + 12));
        }
        foreach (var t in new[] { "head", "hhea", "maxp", "hmtx", "cmap", "loca", "glyf" })
            if (!_tables.ContainsKey(t)) throw new InvalidDataException("not a TrueType font with outlines: no " + t);
        var head = _tables["head"].Off;
        UnitsPerEm = U16(head + 18);
        BBox = (S16(head + 36), S16(head + 38), S16(head + 40), S16(head + 42));
        Italic = (U16(head + 44) & 2) != 0;
        var longLoca = S16(head + 50) == 1;
        var hhea = _tables["hhea"].Off;
        Ascent = S16(hhea + 4) / (double)UnitsPerEm;
        Descent = S16(hhea + 6) / (double)UnitsPerEm;
        int nh = U16(hhea + 34);
        NumGlyphs = U16(_tables["maxp"].Off + 4);
        _adv = new ushort[Math.Max(1, nh)];
        for (var i = 0; i < nh; ++i) _adv[i] = (ushort)U16(_tables["hmtx"].Off + 4 * i);
        _loca = new uint[NumGlyphs + 1];
        var lo = _tables["loca"].Off;
        for (var i = 0; i <= NumGlyphs; ++i) _loca[i] = longLoca ? U32(lo + 4 * i) : (uint)U16(lo + 2 * i) * 2;
        CapHeight = _tables.TryGetValue("OS/2", out var os2) && os2.Len >= 90 && U16(os2.Off) >= 2 ? S16(os2.Off + 88) / (double)UnitsPerEm : Ascent * 0.8;
        ReadCmap();
        Name = ReadName() ?? System.IO.Path.GetFileNameWithoutExtension(path);
    }

    private int U16(int o) => (_d[o] << 8) | _d[o + 1];
    private short S16(int o) => (short)U16(o);
    private uint U32(int o) => ((uint)_d[o] << 24) | ((uint)_d[o + 1] << 16) | ((uint)_d[o + 2] << 8) | _d[o + 3];

    private void ReadCmap()
    {
        var c = _tables["cmap"].Off;
        int n = U16(c + 2);
        int best = -1, rank = -1;
        for (var i = 0; i < n; ++i)
        {
            int pid = U16(c + 4 + 8 * i), eid = U16(c + 6 + 8 * i);
            var off = c + (int)U32(c + 8 + 8 * i);
            int fmt = U16(off);
            var r = pid == 3 && eid == 10 && fmt == 12 ? 4 : pid == 0 && fmt == 12 ? 3 : pid == 3 && eid == 1 && fmt == 4 ? 2 : pid == 0 && fmt == 4 ? 1 : -1;
            if (r > rank) { rank = r; best = off; }
        }
        if (best < 0) return;
        if (U16(best) == 12)
        {
            var groups = (int)U32(best + 12);
            for (var g = 0; g < groups; ++g)
            {
                var o = best + 16 + 12 * g;
                uint a = U32(o), b = U32(o + 4), gid = U32(o + 8);
                for (var cp = a; cp <= b && cp < 0x110000; ++cp) _cmap[(int)cp] = (int)(gid + cp - a);
            }
            return;
        }
        var segX2 = U16(best + 6);
        int ends = best + 14, starts = ends + segX2 + 2, deltas = starts + segX2, ranges = deltas + segX2;
        for (var s = 0; s < segX2 / 2; ++s)
        {
            int end = U16(ends + 2 * s), start = U16(starts + 2 * s), delta = S16(deltas + 2 * s), ro = U16(ranges + 2 * s);
            for (var cp = start; cp <= end && cp != 0xFFFF; ++cp)
            {
                int gid;
                if (ro == 0) gid = (cp + delta) & 0xFFFF;
                else
                {
                    var at = ranges + 2 * s + ro + 2 * (cp - start);
                    gid = U16(at);
                    if (gid != 0) gid = (gid + delta) & 0xFFFF;
                }
                if (gid != 0) _cmap[cp] = gid;
            }
        }
    }

    private string? ReadName()
    {
        if (!_tables.TryGetValue("name", out var t)) return null;
        int n = U16(t.Off + 2), strings = t.Off + U16(t.Off + 4);
        for (var i = 0; i < n; ++i)
        {
            var r = t.Off + 6 + 12 * i;
            int pid = U16(r), nameId = U16(r + 6), len = U16(r + 8), off = U16(r + 10);
            if (nameId != 6) continue;   // the PostScript name
            var raw = _d.AsSpan(strings + off, len).ToArray();
            var s = pid == 3 || pid == 0 ? Encoding.BigEndianUnicode.GetString(raw) : Encoding.ASCII.GetString(raw);
            s = new string(s.Where(ch => ch > ' ' && ch < 127 && ch != '/' && ch != '(' && ch != ')' && ch != '[' && ch != ']' && ch != '<' && ch != '>' && ch != '{' && ch != '}' && ch != '%').ToArray());
            if (s.Length > 0) return s;
        }
        return null;
    }

    public static TrueTypeFont? Load(string path)
    {
        try { return File.Exists(path) ? new TrueTypeFont(path, File.ReadAllBytes(path)) : null; }
        catch { return null; }
    }

    public bool Has(int cp) => _cmap.ContainsKey(cp);
    public int Glyph(int cp) => _cmap.TryGetValue(cp, out var g) ? g : 0;
    /// <summary>The advance of a glyph in em.</summary>
    public double Advance(int gid) => _adv[Math.Min(gid, _adv.Length - 1)] / (double)UnitsPerEm;

    /// <summary>The font with only these glyphs' outlines (and the parts of composite glyphs), ids unchanged: the tables a
    /// PDF reader needs (head, hhea, maxp, hmtx, loca, glyf and the hinting tables).</summary>
    public byte[] Subset(IEnumerable<int> used)
    {
        var keep = new HashSet<int> { 0 };
        var todo = new Stack<int>(used.Where(g => g > 0 && g < NumGlyphs));
        var glyf = _tables["glyf"].Off;
        while (todo.Count > 0)
        {
            var g = todo.Pop();
            if (!keep.Add(g)) continue;
            if (_loca[g + 1] <= _loca[g]) continue;
            var o = glyf + (int)_loca[g];
            if (S16(o) >= 0) continue;   // a simple glyph
            var p = o + 10;
            while (true)
            {
                int flags = U16(p), comp = U16(p + 2);
                todo.Push(comp);
                p += 4 + ((flags & 1) != 0 ? 4 : 2);
                p += (flags & 8) != 0 ? 2 : (flags & 0x40) != 0 ? 4 : (flags & 0x80) != 0 ? 8 : 0;
                if ((flags & 0x20) == 0) break;
            }
        }
        using var g2 = new MemoryStream();
        var loca = new uint[NumGlyphs + 1];
        for (var g = 0; g < NumGlyphs; ++g)
        {
            loca[g] = (uint)g2.Position;
            if (!keep.Contains(g) || _loca[g + 1] <= _loca[g]) continue;
            g2.Write(_d, glyf + (int)_loca[g], (int)(_loca[g + 1] - _loca[g]));
            while (g2.Position % 4 != 0) g2.WriteByte(0);
        }
        loca[NumGlyphs] = (uint)g2.Position;
        var locaBytes = new byte[4 * loca.Length];
        for (var i = 0; i < loca.Length; ++i) { locaBytes[4 * i] = (byte)(loca[i] >> 24); locaBytes[4 * i + 1] = (byte)(loca[i] >> 16); locaBytes[4 * i + 2] = (byte)(loca[i] >> 8); locaBytes[4 * i + 3] = (byte)loca[i]; }
        var head = _d.AsSpan(_tables["head"].Off, _tables["head"].Len).ToArray();
        head[50] = 0; head[51] = 1;                  // long offsets in loca
        head[8] = head[9] = head[10] = head[11] = 0;  // checksum adjustment: not used by PDF readers
        var tables = new SortedDictionary<string, byte[]>(StringComparer.Ordinal)
        {
            ["head"] = head, ["loca"] = locaBytes, ["glyf"] = g2.ToArray(),
        };
        foreach (var t in new[] { "hhea", "maxp", "hmtx", "cvt ", "fpgm", "prep" })
            if (_tables.TryGetValue(t, out var e)) tables[t] = _d.AsSpan(e.Off, e.Len).ToArray();
        using var o2 = new MemoryStream();
        var w = new BinaryWriter(o2);
        void BE16(int v) { w.Write((byte)(v >> 8)); w.Write((byte)v); }
        void BE32(uint v) { w.Write((byte)(v >> 24)); w.Write((byte)(v >> 16)); w.Write((byte)(v >> 8)); w.Write((byte)v); }
        var nt = tables.Count;
        var es = (int)Math.Floor(Math.Log2(nt));
        BE32(0x00010000); BE16(nt); BE16(16 << es); BE16(es); BE16(nt * 16 - (16 << es));
        var off = 12 + 16 * nt;
        foreach (var (tag, data) in tables)
        {
            uint sum = 0;
            for (var i = 0; i < data.Length; i += 4)
                sum += (uint)((i < data.Length ? data[i] << 24 : 0) | (i + 1 < data.Length ? data[i + 1] << 16 : 0) | (i + 2 < data.Length ? data[i + 2] << 8 : 0) | (i + 3 < data.Length ? data[i + 3] : 0));
            w.Write(Encoding.ASCII.GetBytes(tag)); BE32(sum); BE32((uint)off); BE32((uint)data.Length);
            off += (data.Length + 3) & ~3;
        }
        foreach (var data in tables.Values) { w.Write(data); for (var i = data.Length; i % 4 != 0; ++i) w.Write((byte)0); }
        w.Flush();
        return o2.ToArray();
    }
}

/// <summary>The fonts print figures are set in: Arial (or Liberation Sans, then the bundled IBM Plex Sans) with its
/// italic and bold, and fonts with wider coverage for characters Arial lacks (⟨ ⟩ …). Each character is taken from the
/// first font of the chain that has it.</summary>
internal static class PrintFonts
{
    private static readonly string Win = Environment.GetFolderPath(Environment.SpecialFolder.Fonts);
    private static string AppFont(string file) => System.IO.Path.Combine(AppContext.BaseDirectory, "Assets", "Fonts", file);

    private static readonly string[] RegularPaths =
    [
        "/System/Library/Fonts/Supplemental/Arial.ttf", "/Library/Fonts/Arial.ttf", System.IO.Path.Combine(Win, "arial.ttf"),
        "/usr/share/fonts/truetype/msttcorefonts/Arial.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf", "/usr/share/fonts/TTF/LiberationSans-Regular.ttf",
        AppFont("IBMPlexSans-Regular.ttf"),
    ];
    private static readonly string[] ItalicPaths =
    [
        "/System/Library/Fonts/Supplemental/Arial Italic.ttf", "/Library/Fonts/Arial Italic.ttf", System.IO.Path.Combine(Win, "ariali.ttf"),
        "/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf",
        "/usr/share/fonts/liberation-sans/LiberationSans-Italic.ttf", "/usr/share/fonts/TTF/LiberationSans-Italic.ttf",
    ];
    private static readonly string[] BoldPaths =
    [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf", "/Library/Fonts/Arial Bold.ttf", System.IO.Path.Combine(Win, "arialbd.ttf"),
        "/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf", "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf", "/usr/share/fonts/TTF/LiberationSans-Bold.ttf",
        AppFont("IBMPlexSans-Bold.ttf"),
    ];
    // wide coverage for what the first font lacks
    private static readonly string[] WidePaths =
    [
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf", "/Library/Fonts/Arial Unicode.ttf", "/System/Library/Fonts/Apple Symbols.ttf",
        System.IO.Path.Combine(Win, "seguisym.ttf"), System.IO.Path.Combine(Win, "cambria.ttc"),
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/dejavu/DejaVuSans.ttf", "/usr/share/fonts/TTF/DejaVuSans.ttf",
    ];

    private static readonly Lazy<TrueTypeFont?> Reg = new(() => First(RegularPaths));
    private static readonly Lazy<TrueTypeFont?> Ita = new(() => First(ItalicPaths));
    private static readonly Lazy<TrueTypeFont?> Bold = new(() => First(BoldPaths));
    private static readonly Lazy<List<TrueTypeFont>> Wide = new(() => WidePaths.Select(TrueTypeFont.Load).OfType<TrueTypeFont>().ToList());
    private static TrueTypeFont? First(IEnumerable<string> paths) => paths.Select(TrueTypeFont.Load).FirstOrDefault(f => f != null);

    /// <summary>The primary font for a style (null: no TrueType font found). Synthetic: an italic drawn slanted.</summary>
    public static (TrueTypeFont? Font, bool Synthetic) For(bool italic, bool bold)
    {
        if (bold && Bold.Value is { } b) return (b, false);
        if (italic && Ita.Value is { } i) return (i, false);
        return (Reg.Value, italic);
    }

    /// <summary>The family name an SVG asks for (the viewer substitutes per glyph like the PDF chain does).</summary>
    public static string SvgFamily => Reg.Value?.Path.Contains("Plex") == true ? "'IBM Plex Sans', Arial, Helvetica, sans-serif" : "Arial, Helvetica, 'Liberation Sans', sans-serif";

    /// <summary>The text split into pieces each set in one font: the style's font where it has the character, else the
    /// first wide font that has it, else the style's font (its missing-glyph box).</summary>
    public static List<(TrueTypeFont Font, bool Synthetic, string Text)> Runs(string s, bool italic, bool bold)
    {
        var (main, synth) = For(italic, bold);
        var outp = new List<(TrueTypeFont, bool, string)>();
        if (main == null) return outp;
        var sb = new StringBuilder();
        TrueTypeFont? cur = null;
        var curSynth = false;
        for (var i = 0; i < s.Length; ++i)
        {
            var cp = char.ConvertToUtf32(s, i);
            if (char.IsHighSurrogate(s[i])) ++i;
            var (f, sy) = main.Has(cp) || cp == ' ' ? (main, synth) : Wide.Value.FirstOrDefault(w => w.Has(cp)) is { } w2 ? (w2, italic) : (main, synth);
            if (cur != null && (f != cur || sy != curSynth)) { outp.Add((cur, curSynth, sb.ToString())); sb.Clear(); }
            cur = f; curSynth = sy;
            sb.Append(char.ConvertFromUtf32(cp));
        }
        if (cur != null && sb.Length > 0) outp.Add((cur, curSynth, sb.ToString()));
        return outp;
    }

    public static double Measure(string s, double size, bool italic, bool bold) =>
        Runs(s, italic, bold).Sum(r => Codepoints(r.Text).Sum(cp => r.Font.Advance(r.Font.Glyph(cp)))) * size;

    public static IEnumerable<int> Codepoints(string s)
    {
        for (var i = 0; i < s.Length; ++i)
        {
            yield return char.ConvertToUtf32(s, i);
            if (char.IsHighSurrogate(s[i])) ++i;
        }
    }

    public static bool Available => Reg.Value != null;
    public static double AscentEm => Reg.Value?.Ascent ?? 0.905;
    public static double LineEm => Reg.Value is { } f ? f.Ascent - f.Descent : 1.15;
}
