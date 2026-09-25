using System.Collections.ObjectModel;
using System.Globalization;
using Avalonia;
using Avalonia.Media.Imaging;
using Avalonia.Platform;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A figure size: a journal column or poster panel (width in mm at a dpi), or a slide in pixels.</summary>
public sealed record FigurePreset(string Name, double Mm, int Dpi, int SlidePx)
{
    public string WidthText => Mm > 0 ? Mm.ToString("0.#", CultureInfo.InvariantCulture) + " mm" : "—";
    public string DpiText => Dpi > 0 ? Dpi.ToString(CultureInfo.InvariantCulture) : "—";
    public int Pixels => SlidePx > 0 ? SlidePx : (int)Math.Round(Mm / 25.4 * Dpi);
    public string PixelsText => $"{Pixels} × {(int)Math.Round(Pixels * 9 / 16.0)}";
}

/// <summary>Everything drawn over the structure: title, scale bar and colour legend, sized from the image width so a
/// preview tile and the exported file look the same.</summary>
public sealed record FigureOverlay(string? Title, double BarAngstrom, double BarPx, bool Legend, string LegendLo, string LegendHi, string LegendTitle,
                                  uint Ink, double Width, double Height)
{
    // ≈ 8 pt title and 6.5 pt labels on a 85 mm column (journal minimum sizes), in proportion on other widths
    public double Margin => Width * 0.03;
    public double TitlePx => Math.Max(9, Width * 0.03);
    public double LabelPx => Math.Max(8, Width * 0.024);
    public double BarThickness => Math.Max(2, Width * 0.005);
    public string BarLabel => BarAngstrom.ToString("0.##", CultureInfo.InvariantCulture) + " Å";
}

/// <summary>One background in the comparison row: its preview and the overlay drawn on it.</summary>
public sealed class FigureTile : System.ComponentModel.INotifyPropertyChanged
{
    public event System.ComponentModel.PropertyChangedEventHandler? PropertyChanged;
    private void Raise(string n) => PropertyChanged?.Invoke(this, new System.ComponentModel.PropertyChangedEventArgs(n));
    public int Background { get; init; }
    public string Name { get; init; } = "";
    public string Use { get; init; } = "";
    private Bitmap? _image;
    private FigureOverlay? _overlay;
    private bool _selected;
    private double _height = 247;
    public double Height { get => _height; set { _height = value; Raise(nameof(Height)); } }
    public Bitmap? Image { get => _image; set { _image = value; Raise(nameof(Image)); } }
    public FigureOverlay? Overlay { get => _overlay; set { _overlay = value; Raise(nameof(Overlay)); } }
    public bool Selected { get => _selected; set { _selected = value; Raise(nameof(Selected)); } }
    public bool Checker => Background == 2;
    public Avalonia.Media.IBrush Fill => Background switch
    {
        0 => new Avalonia.Media.SolidColorBrush(Avalonia.Media.Color.Parse("#0F1113")),
        1 => Avalonia.Media.Brushes.White,
        _ => Avalonia.Media.Brushes.Transparent,
    };
}

/// <summary>Export › Figure (design/boards/FigureBackground): the same frame on dark, white and transparent backgrounds,
/// journal and slide sizes, and a title, true scale bar and colour legend in ink that suits the background.</summary>
public sealed partial class MainViewModel
{
    public bool IsFigure => _module == 18;

    public static readonly FigurePreset[] FigurePresets =
    [
        new("Journal · single column", 85, 600, 0),
        new("Journal · double column", 170, 600, 0),
        new("Poster panel", 250, 300, 0),
        new("Slide 16:9", 0, 0, 1920),
    ];
    public static readonly string[] FigureFormats = ["PNG · RGBA", "SVG · vector"];
    public static readonly string[] FigureAspects = ["16:9", "4:3", "1:1", "As the view"];

    public ObservableCollection<FigureTile> FigureTiles { get; } =
    [
        new FigureTile { Background = 0, Name = "Dark", Use = "screen, talks" },
        new FigureTile { Background = 1, Name = "White", Use = "papers, print" },
        new FigureTile { Background = 2, Name = "Transparent", Use = "slides, posters, overlays" },
    ];

    private int _figFormat, _figAspect;
    private decimal _figWidthMm = 85, _figDpi = 600;
    private bool _figCell = true, _figBar = true, _figTitle = true, _figLegend;
    private string _figTitleText = "", _figCustom = "none";
    private double _viewAspect = 16 / 9.0;

    public void OpenFigure()
    {
        if (_doc == null) { Status = "Open or build a structure first, then export it as a figure"; return; }
        var p = FigurePresets[Math.Clamp(_settings.FigurePreset, 0, FigurePresets.Length - 1)];
        _figWidthMm = p.Mm > 0 ? (decimal)p.Mm : 325;
        _figDpi = p.Dpi > 0 ? p.Dpi : 150;
        if (_figTitleText.Length == 0) _figTitleText = FigureDefaultTitle();
        foreach (var t in FigureTiles) t.Selected = t.Background == FigBackground;
        RaiseFigure();
        SetModule(18);
        RefreshFigure();
    }

    private string FigureDefaultTitle()
    {
        if (_doc == null) return "";
        var s = _doc.Summary();
        return string.Format(CultureInfo.InvariantCulture, "{0} · {1:N0} atoms", Title.Replace(" (unsaved)", ""), s.Atoms);
    }

    public int FigBackground
    {
        get => _settings.FigureBackground;
        set
        {
            if (_settings.FigureBackground == value) return;
            _settings.FigureBackground = Math.Clamp(value, 0, 2);
            foreach (var t in FigureTiles) t.Selected = t.Background == value;
            Raise();
            Raise(nameof(FigBgDark)); Raise(nameof(FigBgWhite)); Raise(nameof(FigBgTransparent));
            Changed("Figure background");
        }
    }
    public bool FigBgDark { get => FigBackground == 0; set { if (value) FigBackground = 0; } }
    public bool FigBgWhite { get => FigBackground == 1; set { if (value) FigBackground = 1; } }
    public bool FigBgTransparent { get => FigBackground == 2; set { if (value) FigBackground = 2; } }

    /// <summary>A custom background colour (#RRGGBB) replaces the choice above; "none" keeps it.</summary>
    public string FigCustom { get => _figCustom; set { if (Set(ref _figCustom, value.Trim())) { RaiseFigure(); RefreshFigure(); } } }
    private uint? FigCustomRgb
    {
        get
        {
            var t = _figCustom.TrimStart('#');
            return t.Length == 6 && uint.TryParse(t, NumberStyles.HexNumber, CultureInfo.InvariantCulture, out var v) ? v : null;
        }
    }

    public int FigPreset
    {
        get => _settings.FigurePreset;
        set
        {
            var v = Math.Clamp(value, 0, FigurePresets.Length - 1);
            _settings.FigurePreset = v;
            var p = FigurePresets[v];
            if (p.Mm > 0) { _figWidthMm = (decimal)p.Mm; _figDpi = p.Dpi; }
            else { _figDpi = 150; _figWidthMm = Math.Round(1920m / 150m * 25.4m, 1); _figAspect = 0; Raise(nameof(FigAspect)); }
            Raise(); Changed("Figure size");
            RaiseFigure();
            RefreshFigure();
        }
    }
    public int FigFormat { get => _figFormat; set { if (Set(ref _figFormat, value)) RaiseFigure(); } }
    public int FigAspect { get => _figAspect; set { if (Set(ref _figAspect, value)) { RaiseFigure(); RefreshFigure(); } } }
    public decimal FigWidthMm { get => _figWidthMm; set { if (Set(ref _figWidthMm, Math.Clamp(value, 10, 2000))) RaiseFigure(); } }
    public decimal FigDpi { get => _figDpi; set { if (Set(ref _figDpi, Math.Clamp(Math.Round(value), 72, 2400))) RaiseFigure(); } }
    public bool FigCell { get => _figCell; set { if (Set(ref _figCell, value)) RefreshFigure(); } }
    public bool FigScaleBar { get => _figBar; set { if (Set(ref _figBar, value)) RefreshFigure(); } }
    public bool FigTitle { get => _figTitle; set { if (Set(ref _figTitle, value)) RefreshFigure(); } }
    public string FigTitleText { get => _figTitleText; set { if (Set(ref _figTitleText, value)) RefreshFigure(); } }
    public bool FigLegend { get => _figLegend; set { if (Set(ref _figLegend, value)) RefreshFigure(); } }
    public bool FigLegendAvailable => _colour == 3;
    public bool FigPerspective => _perspective;

    private double FigAspectRatio => _figAspect switch { 0 => 16 / 9.0, 1 => 4 / 3.0, 2 => 1.0, _ => _viewAspect };
    /// <summary>The exported image in pixels: width from mm × dpi (or the slide's 1920), height from the aspect.</summary>
    public (int W, int H) FigPixels
    {
        get
        {
            var w = FigurePresets[FigPreset].SlidePx > 0 && _figAspect == 0 ? 1920 : (int)Math.Round((double)_figWidthMm / 25.4 * (double)_figDpi);
            w = Math.Clamp(w, 64, 12000);
            return (w, Math.Clamp((int)Math.Round(w / FigAspectRatio), 64, 12000));
        }
    }
    public string FigPixelsText { get { var (w, h) = FigPixels; return $"{w} × {h}"; } }
    public string FigNote => _perspective ? "Perspective view: the scale bar is true only at the centre of the cell. Switch to orthographic for a figure." :
        FigFormat == 1 ? "SVG keeps atoms as vector shapes; transparent SVG has no background shape." :
        "PNG keeps an alpha channel on Transparent; the dpi is written into the file.";

    private void RaiseFigure()
    {
        Raise(nameof(FigWidthMm)); Raise(nameof(FigDpi)); Raise(nameof(FigPixelsText)); Raise(nameof(FigNote));
        Raise(nameof(FigLegendAvailable)); Raise(nameof(FigPerspective)); Raise(nameof(FigCustom)); Raise(nameof(FigTitleText));
    }

    public void SetViewAspect(double w, double h) { if (w > 10 && h > 10) _viewAspect = w / h; }

    /// <summary>Render options of the figure for one background at w × h pixels.</summary>
    public CapsRenderOpts FigureOptions(int background, int w, int h, int supersample)
    {
        var o = ExportOptions(w, h);
        o.Supersample = supersample;
        o.ShowCell = _figCell ? 1 : 0;
        o.Background = background;
        if (FigCustomRgb is { } rgb) { o.Background = 3; o.CustomRgb = rgb; }
        return o;
    }

    /// <summary>The overlay for an image of w × h pixels: the scale bar a round length near 12 % of the width.</summary>
    public FigureOverlay FigureOverlayFor(int background, int w, int h, double pxPerAngstrom)
    {
        uint ink = background == 0 ? 0xE9ECEFu : 0x141413u;
        if (FigCustomRgb is { } rgb)
        {
            var lum = 0.2126 * ((rgb >> 16) & 255) + 0.7152 * ((rgb >> 8) & 255) + 0.0722 * (rgb & 255);
            ink = lum < 128 ? 0xE9ECEFu : 0x141413u;
        }
        double barA = 0, barPx = 0;
        if (_figBar && pxPerAngstrom > 0)
        {
            var target = 0.12 * w / pxPerAngstrom;
            barA = new[] { 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000.0 }.LastOrDefault(x => x <= target);
            if (barA <= 0) barA = 1;
            barPx = barA * pxPerAngstrom;
        }
        return new FigureOverlay(_figTitle && _figTitleText.Length > 0 ? _figTitleText : null, barA, barPx, _figLegend && _colour == 3,
                                 LegendLo, LegendHi, "distance to molecule centre", ink, w, h);
    }

    private int _figGen;
    /// <summary>Re-renders the three preview tiles (off the UI thread), keeping the camera of the view.</summary>
    public void RefreshFigure()
    {
        if (!IsFigure || _doc == null) return;
        var gen = ++_figGen;
        var doc = _doc;
        const int tw = 440;
        var th = (int)Math.Round(tw / FigAspectRatio);
        var cam = Camera;
        var jobs = FigureTiles.Select(t => (Tile: t, Opt: FigureOptions(t.Background, tw, th, 2))).ToList();
        Task.Run(() =>
        {
            foreach (var (tile, opt) in jobs)
            {
                var rgba = new byte[tw * th * 4];
                double scale;
                try { doc.Render(cam, opt, rgba); scale = doc.ViewScale(cam, opt); }
                catch { return; }
                Avalonia.Threading.Dispatcher.UIThread.Post(() =>
                {
                    if (gen != _figGen || !IsFigure) return;   // the page was left: drop the stale preview
                    tile.Height = th;
                    tile.Image = ToBitmap(rgba, tw, th);
                    tile.Overlay = FigureOverlayFor(tile.Background, tw, th, scale);
                });
            }
        });
    }

    /// <summary>Renders the figure at full size and writes it (PNG with overlay and dpi, or SVG with the overlay as
    /// vector marks). compose draws the PNG overlay (Avalonia, in the view layer).</summary>
    public async Task<string> ExportFigure(string path, Action<byte[], int, int, FigureOverlay, double, string> composePng, Action<string, FigureOverlay> addSvg)
    {
        if (_doc == null) throw new InvalidOperationException("nothing open");
        var doc = _doc;
        var (w, h) = FigPixels;
        var bg = FigBackground;
        var cam = Camera;
        var opt = FigureOptions(bg, w, h, w * h > 6_000_000 ? 1 : 2);
        var dpi = FigurePresets[FigPreset].SlidePx > 0 && _figAspect == 0 ? 0 : (double)_figDpi;
        var svg = _figFormat == 1;
        var scale = doc.ViewScale(cam, opt);
        var overlay = FigureOverlayFor(bg, w, h, scale);
        Status = $"Exporting the figure · {w} × {h}";
        if (svg)
        {
            await Task.Run(() => doc.ExportSvg(cam, opt, path));
            addSvg(path, overlay);
        }
        else
        {
            var rgba = new byte[w * h * 4];
            await Task.Run(() => doc.Render(cam, opt, rgba));
            composePng(rgba, w, h, overlay, dpi, path);
        }
        var what = $"{w} × {h}" + (dpi > 0 ? $" · {dpi:0} dpi" : "") + $" · {FigureTiles[bg].Name.ToLowerInvariant()}" + (overlay.BarPx > 0 ? $" · scale bar {overlay.BarLabel}" : "");
        Status = $"Wrote {path} · {what}";
        return what;
    }

    internal static WriteableBitmap ToBitmap(byte[] rgba, int w, int h)
    {
        var bmp = new WriteableBitmap(new PixelSize(w, h), new Vector(96, 96), PixelFormat.Rgba8888, AlphaFormat.Unpremul);
        using var fb = bmp.Lock();
        for (int y = 0; y < h; ++y)
            System.Runtime.InteropServices.Marshal.Copy(rgba, y * w * 4, fb.Address + y * fb.RowBytes, w * 4);
        return bmp;
    }
}
