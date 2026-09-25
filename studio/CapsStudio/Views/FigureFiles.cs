using System;
using System.IO;
using System.IO.Compression;
using System.Text;

namespace CapsStudio.Views;

/// <summary>TIFF and PDF from a rendered page (RGBA, top row first). TIFF: baseline, uncompressed RGB with the
/// resolution tags. PDF: one page of the figure's size in points holding the page as a Flate-compressed image — a raster
/// at the chosen dpi, not vector (use SVG for vector plots).</summary>
public static class FigureFiles
{
    public static void WriteTiff(string path, byte[] rgba, int w, int h, int dpi)
    {
        using var f = new BinaryWriter(File.Create(path));
        var pixels = w * h * 3;
        const int entries = 12;
        var ifd = 8;
        var after = ifd + 2 + entries * 12 + 4;
        var bps = after;             // 3 × SHORT bits per sample
        var xres = bps + 6;          // RATIONAL
        var yres = xres + 8;
        var data = yres + 8;
        f.Write((byte)'I'); f.Write((byte)'I'); f.Write((ushort)42); f.Write(ifd);
        f.Write((ushort)entries);
        void Tag(ushort tag, ushort type, int count, int value) { f.Write(tag); f.Write(type); f.Write(count); f.Write(value); }
        Tag(256, 4, 1, w);                  // width
        Tag(257, 4, 1, h);                  // height
        Tag(258, 3, 3, bps);                // bits per sample → 8, 8, 8
        Tag(259, 3, 1, 1);                  // no compression
        Tag(262, 3, 1, 2);                  // RGB
        Tag(273, 4, 1, data);               // strip offset
        Tag(277, 3, 1, 3);                  // samples per pixel
        Tag(278, 4, 1, h);                  // rows per strip
        Tag(279, 4, 1, pixels);             // strip bytes
        Tag(282, 5, 1, xres);               // x resolution
        Tag(283, 5, 1, yres);               // y resolution
        Tag(296, 3, 1, 2);                  // inches
        f.Write(0);                         // no next IFD
        f.Write((ushort)8); f.Write((ushort)8); f.Write((ushort)8);
        f.Write(dpi); f.Write(1); f.Write(dpi); f.Write(1);
        var row = new byte[w * 3];
        for (var y = 0; y < h; ++y)
        {
            for (var x = 0; x < w; ++x) { var i = (y * w + x) * 4; row[x * 3] = rgba[i]; row[x * 3 + 1] = rgba[i + 1]; row[x * 3 + 2] = rgba[i + 2]; }
            f.Write(row);
        }
    }

    /// <summary>Adds (or replaces) the PNG's pHYs chunk so readers see the figure's dpi.</summary>
    public static void SetPngDpi(string path, int dpi)
    {
        var b = File.ReadAllBytes(path);
        if (b.Length < 33 || b[12] != 'I' || b[13] != 'H') return;
        var ppm = (uint)Math.Round(dpi / 0.0254);
        var chunk = new byte[21];
        void Be(int at, uint v) { chunk[at] = (byte)(v >> 24); chunk[at + 1] = (byte)(v >> 16); chunk[at + 2] = (byte)(v >> 8); chunk[at + 3] = (byte)v; }
        Be(0, 9);
        chunk[4] = (byte)'p'; chunk[5] = (byte)'H'; chunk[6] = (byte)'Y'; chunk[7] = (byte)'s';
        Be(8, ppm); Be(12, ppm); chunk[16] = 1;
        uint crc = 0xFFFFFFFF;
        for (var i = 4; i < 17; ++i)
        {
            crc ^= chunk[i];
            for (var k = 0; k < 8; ++k) crc = (crc & 1) != 0 ? 0xEDB88320 ^ (crc >> 1) : crc >> 1;
        }
        Be(17, crc ^ 0xFFFFFFFF);
        using var ms = new MemoryStream();
        ms.Write(b, 0, 33);   // signature + IHDR
        ms.Write(chunk);
        var pos = 33;
        while (pos + 8 <= b.Length)   // the rest, without an existing pHYs
        {
            var len = (b[pos] << 24) | (b[pos + 1] << 16) | (b[pos + 2] << 8) | b[pos + 3];
            var type = Encoding.ASCII.GetString(b, pos + 4, 4);
            if (type != "pHYs") ms.Write(b, pos, len + 12);
            pos += len + 12;
        }
        File.WriteAllBytes(path, ms.ToArray());
    }

    public static void WritePdf(string path, byte[] rgba, int w, int h, double widthPt, double heightPt)
    {
        var rgb = new byte[w * h * 3];
        for (var i = 0; i < w * h; ++i) { rgb[i * 3] = rgba[i * 4]; rgb[i * 3 + 1] = rgba[i * 4 + 1]; rgb[i * 3 + 2] = rgba[i * 4 + 2]; }
        byte[] packed;
        using (var ms = new MemoryStream())
        {
            using (var z = new ZLibStream(ms, CompressionLevel.Optimal, true)) z.Write(rgb, 0, rgb.Length);
            packed = ms.ToArray();
        }
        var inv = System.Globalization.CultureInfo.InvariantCulture;
        var content = Encoding.ASCII.GetBytes(string.Format(inv, "q {0:0.###} 0 0 {1:0.###} 0 0 cm /Im0 Do Q\n", widthPt, heightPt));
        using var f = File.Create(path);
        var offsets = new long[6];
        void W(string s) { var b = Encoding.ASCII.GetBytes(s); f.Write(b, 0, b.Length); }
        W("%PDF-1.4\n");
        offsets[1] = f.Position; W("1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj\n");
        offsets[2] = f.Position; W("2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj\n");
        offsets[3] = f.Position;
        W(string.Format(inv, "3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 {0:0.###} {1:0.###}] /Resources << /XObject << /Im0 5 0 R >> >> /Contents 4 0 R >> endobj\n", widthPt, heightPt));
        offsets[4] = f.Position; W($"4 0 obj << /Length {content.Length} >> stream\n"); f.Write(content); W("endstream endobj\n");
        offsets[5] = f.Position;
        W($"5 0 obj << /Type /XObject /Subtype /Image /Width {w} /Height {h} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length {packed.Length} >> stream\n");
        f.Write(packed); W("\nendstream endobj\n");
        var xref = f.Position;
        W("xref\n0 6\n0000000000 65535 f \n");
        for (var k = 1; k <= 5; ++k) W($"{offsets[k]:D10} 00000 n \n");
        W($"trailer << /Size 6 /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n");
    }
}
