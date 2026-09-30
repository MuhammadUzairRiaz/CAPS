using System.Globalization;
using System.Text.Json;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>A size preset of the export dialog.</summary>
public sealed record ExportPreset(string Name, string Detail, int W, int H, int Dpi);

/// <summary>Export image / movie (design/boards/ExportDialog): the view re-rendered at full resolution with the same
/// camera and styles — size presets, PNG at 8 or 16 bits, background, supersampling, colour profile, labels and
/// measurements, the provenance manifest in the file's metadata — or a movie of the frames (or a turntable of one
/// structure) as an animated PNG, a PNG sequence, or MP4 when ffmpeg is installed.</summary>
public sealed partial class MainViewModel
{
    public static readonly ExportPreset[] ExportPresets =
    [
        new("Screen", "1920 × 1080", 1920, 1080, 0),
        new("4K", "3840 × 2160", 3840, 2160, 0),
        new("Journal figure", "3.5 in · 600 dpi", 2100, 1575, 600),
        new("Poster", "7680 × 4320", 7680, 4320, 0),
    ];
    public static readonly string[] ExportImageFormats = ["PNG · 8-bit", "PNG · 16-bit", "SVG · vector", "POV-Ray scene (.pov)", "glTF 3D · Blender, ParaView (.glb)", "OBJ 3D (.obj + .mtl)"];
    /// <summary>A 3D scene for another renderer or tool (POV-Ray, glTF, OBJ), not an image.</summary>
    public bool ExportIsScene => _expFormat >= 3;
    /// <summary>A raster image (PNG): the renderer section applies.</summary>
    public bool ExportIsRaster => _expTab == 0 && _expFormat <= 1;
    public static readonly string[] ExportBackgrounds = ["Dark", "White", "Transparent"];
    public static readonly string[] ExportSupersampling = ["1×", "4×", "9×", "16×"];
    public static readonly string[] ExportProfiles = ["sRGB", "Untagged"];
    public static string[] ExportMovieFormats => FfmpegPath != null
        ? ["Animated PNG", "PNG sequence", "MP4 · H.264 (ffmpeg)"]
        : ["Animated PNG", "PNG sequence"];
    public static readonly string[] ExportMovieSources = ["Trajectory frames", "Turntable (360° about the vertical)"];

    private bool _expOpen;
    public bool ExportDialogOpen { get => _expOpen; set { if (Set(ref _expOpen, value) && value) RaiseExport(); } }
    private int _expTab, _expPreset = 1, _expFormat = 1, _expBg = 1, _expSs = 1, _expProfile, _movFormat, _movSource, _movFps = 24, _movStep = 1, _movTurnFrames = 120;
    private decimal _expW = 3840, _expH = 2160;
    private bool _expLabels = true, _expProv = true, _expBusy;
    private string _expProgress = "", _expResult = "";

    public int ExportTab { get => _expTab; set { if (Set(ref _expTab, value)) { Raise(nameof(ExportIsImage)); Raise(nameof(ExportIsMovie)); RaiseExport(); } } }
    public bool ExportIsImage => _expTab == 0;
    public bool ExportIsMovie => _expTab == 1;
    public int ExportPreset
    {
        get => _expPreset;
        set
        {
            if (!Set(ref _expPreset, value) || value < 0 || value >= ExportPresets.Length) return;
            _expW = ExportPresets[value].W;
            _expH = ExportPresets[value].H;
            Raise(nameof(ExportWidth)); Raise(nameof(ExportHeight));
            RaiseExport();
        }
    }
    public decimal ExportWidth { get => _expW; set { if (Set(ref _expW, Math.Clamp(Math.Round(value), 16, 16384))) { MatchPreset(); RaiseExport(); } } }
    public decimal ExportHeight { get => _expH; set { if (Set(ref _expH, Math.Clamp(Math.Round(value), 16, 16384))) { MatchPreset(); RaiseExport(); } } }
    private void MatchPreset()
    {
        var k = Array.FindIndex(ExportPresets, p => p.W == (int)_expW && p.H == (int)_expH);
        if (k != _expPreset) { _expPreset = k; Raise(nameof(ExportPreset)); }
    }
    public int ExportFormat { get => _expFormat; set { if (Set(ref _expFormat, value)) RaiseExport(); } }
    public int ExportDlgBackground { get => _expBg; set { if (Set(ref _expBg, value)) RaiseExport(); } }
    public int ExportSupersample { get => _expSs; set => Set(ref _expSs, value); }
    // the renderer of an image export: as the view (the rasteriser) or CAPS's ray tracer (occlusion, shadows, depth of field)
    public static readonly string[] ExportEngines = ["As the view", "Ray traced"];
    public static readonly string[] ExportTraceQualities = ["Draft · 16 samples", "Good · 64 samples", "Final · 256 samples"];
    public static readonly string[] ExportDepthOfField = ["Off", "Soft", "Strong"];
    private int _expEngine, _expTraceQ = 1, _expDof;
    private bool _expShadows = true, _expOcclusion = true, _expTraceOutlines;
    public int ExportEngine { get => _expEngine; set { if (Set(ref _expEngine, Math.Clamp(value, 0, 1))) Raise(nameof(ExportIsTraced)); } }
    public bool ExportIsTraced => _expEngine == 1;
    public int ExportTraceQuality { get => _expTraceQ; set => Set(ref _expTraceQ, Math.Clamp(value, 0, 2)); }
    public int ExportDof { get => _expDof; set => Set(ref _expDof, Math.Clamp(value, 0, 2)); }
    public bool ExportShadows { get => _expShadows; set => Set(ref _expShadows, value); }
    public bool ExportOcclusion { get => _expOcclusion; set => Set(ref _expOcclusion, value); }
    public bool ExportTraceOutlines { get => _expTraceOutlines; set => Set(ref _expTraceOutlines, value); }
    public int ExportProfile { get => _expProfile; set => Set(ref _expProfile, value); }
    public bool ExportLabels { get => _expLabels; set { if (Set(ref _expLabels, value)) RaiseExport(); } }
    public bool ExportProvenance { get => _expProv; set => Set(ref _expProv, value); }
    public int MovieFormat { get => _movFormat; set { if (Set(ref _movFormat, value)) RaiseExport(); } }
    public int MovieSource { get => _movSource; set { if (Set(ref _movSource, value)) RaiseExport(); } }
    public decimal MovieFps { get => _movFps; set { if (Set(ref _movFps, (int)Math.Clamp(value, 1, 120))) RaiseExport(); } }
    public decimal MovieStep { get => _movStep; set { if (Set(ref _movStep, (int)Math.Clamp(value, 1, 100000))) RaiseExport(); } }
    public decimal MovieTurnFrames { get => _movTurnFrames; set { if (Set(ref _movTurnFrames, (int)Math.Clamp(value, 12, 1440))) RaiseExport(); } }
    public bool MovieIsTurntable => _movSource == 1;
    public bool ExportDlgBusy { get => _expBusy; private set { if (Set(ref _expBusy, value)) Raise(nameof(ExportDlgIdle)); } }
    public bool ExportDlgIdle => !_expBusy;
    public string ExportProgress { get => _expProgress; private set => Set(ref _expProgress, value); }
    public string ExportResult { get => _expResult; private set => Set(ref _expResult, value); }
    public event Action? ExportPreviewChanged;

    public string ExportSizeChip => $"{_expW:0} × {_expH:0} px";
    public string ExportAspectChip
    {
        get
        {
            static int Gcd(int a, int b) => b == 0 ? a : Gcd(b, a % b);
            var (w, h) = ((int)_expW, (int)_expH);
            var g = Gcd(w, h);
            var (a, b) = (w / g, h / g);
            return a <= 32 && b <= 32 ? $"{a} : {b}" : (w / (double)h).ToString("0.00", CultureInfo.InvariantCulture) + " : 1";
        }
    }
    public string ExportBgChip => "preview · " + ExportBackgrounds[_expBg].ToLowerInvariant() + " background";
    public string ExportButton => _expTab == 1
        ? _movFormat switch { 1 => "Export PNG sequence", 2 => "Export MP4", _ => "Export animated PNG" }
        : _expFormat switch { 2 => "Export SVG", 3 => "Export POV-Ray scene", 4 => "Export glTF", 5 => "Export OBJ", _ => "Export PNG" };
    public string MovieFramesText
    {
        get
        {
            if (_movSource == 1) return $"{_movTurnFrames} frames · {_movTurnFrames / (double)_movFps:0.0} s at {_movFps} fps";
            var n = Frames;
            var count = (n - 1) / _movStep + 1;
            return $"{count} of {n} frame{(n == 1 ? "" : "s")} · {count / (double)_movFps:0.0} s at {_movFps} fps" + (n == 1 ? " · one frame: use a turntable" : "");
        }
    }
    public bool ExportLabelsAvailable => _expTab == 0 && _expFormat != 2;

    private void RaiseExport()
    {
        foreach (var n in new[] { nameof(ExportSizeChip), nameof(ExportAspectChip), nameof(ExportBgChip), nameof(ExportButton), nameof(MovieFramesText),
                                  nameof(MovieIsTurntable), nameof(ExportLabelsAvailable), nameof(ExportIsRaster), nameof(ExportIsScene) }) Raise(n);
        ExportPreviewChanged?.Invoke();
    }

    public void OpenExportDialog(int tab = 0)
    {
        if (_doc == null) { Status = "Open a structure first"; return; }
        _expTab = tab;
        Raise(nameof(ExportTab)); Raise(nameof(ExportIsImage)); Raise(nameof(ExportIsMovie));
        ExportResult = "";
        ExportProgress = "";
        ExportDialogOpen = true;
        RaiseExport();
    }

    /// <summary>The options the core renders with, at w × h (the preview uses a smaller size with the same framing).</summary>
    public CapsRenderOpts ExportDialogOptions(int w, int h, bool final)
    {
        var o = ViewOptions(w, h, final ? _expSs + 1 : 2);
        if (final) o.LodNear = o.LodFar = 0;   // exports in full detail
        o.Background = _expBg;
        if (!(_expLabels && _expTab == 0)) o.Highlight0 = o.Highlight1 = o.Highlight2 = o.Highlight3 = -1;
        o.Focus = 0;
        return o;
    }

    /// <summary>The preview (the export's framing at a reduced size), RGBA.</summary>
    public (byte[] Rgba, int W, int H)? ExportPreview(int maxW, int maxH)
    {
        if (_doc == null) return null;
        var scale = Math.Min(maxW / (double)_expW, maxH / (double)_expH);
        var (w, h) = (Math.Max(16, (int)(_expW * (decimal)scale)), Math.Max(16, (int)(_expH * (decimal)scale)));
        var rgba = new byte[w * h * 4];
        _doc.Render(Camera, ExportDialogOptions(w, h, false), rgba);
        return (rgba, w, h);
    }

    /// <summary>The labels (the view's atom labels) and the measurement between picked atoms at w × h, for the overlay.</summary>
    public (List<ViewLabel> Labels, (double X1, double Y1, double X2, double Y2)? Line, string Measure) ExportMarks(int w, int h)
    {
        var labels = new List<ViewLabel>();
        if (_doc == null || !(_expLabels && _expTab == 0)) return (labels, null, "");
        var opt = ExportDialogOptions(w, h, true);
        labels = ViewLabels(Camera, opt, 1.0);
        (double, double, double, double)? line = null;
        if (_selection.Count >= 2)
        {
            var n = (int)_doc.Summary().Atoms;
            var p = _doc.ProjectAtoms(Camera, opt, n);
            int a = _selection[0], b = _selection[^1];
            if (a < n && b < n) line = (p[3 * a], p[3 * a + 1], p[3 * b], p[3 * b + 1]);
        }
        return (labels, line, MeasureText);
    }

    /// <summary>ffmpeg on the PATH (for MP4), or null.</summary>
    public static string? FfmpegPath
    {
        get
        {
            foreach (var dir in (Environment.GetEnvironmentVariable("PATH") ?? "").Split(Path.PathSeparator).Concat(["/opt/homebrew/bin", "/usr/local/bin"]))
            {
                var f = Path.Combine(dir, OperatingSystem.IsWindows() ? "ffmpeg.exe" : "ffmpeg");
                if (File.Exists(f)) return f;
            }
            return null;
        }
    }

    public string ExportDefaultName => Path.GetFileNameWithoutExtension(_doc?.Path ?? "structure") +
        (_expTab == 1 ? _movFormat switch { 1 => "_frames", 2 => ".mp4", _ => ".png" } : _expFormat switch { 2 => ".svg", 3 => ".pov", 4 => ".glb", 5 => ".obj", _ => ".png" });

    /// <summary>Writes the image: PNG through the core (bits, dpi, profile, manifest, overlay) or SVG.</summary>
    public async Task<string?> ExportDialogImage(string path, Func<int, int, byte[]?> overlay)
    {
        if (_doc == null || _expBusy) return null;
        var doc = _doc;
        var (w, h) = ((int)_expW, (int)_expH);
        var opt = ExportDialogOptions(w, h, true);
        var cam = Camera;
        ExportDlgBusy = true;
        ExportProgress = $"Rendering {w} × {h}…";
        try
        {
            if (_expFormat == 2)
            {
                await Task.Run(() => doc.ExportSvg(cam, opt, path));
            }
            else if (_expFormat >= 3)
            {
                var fmt = _expFormat switch { 3 => "pov", 5 => "obj", _ => "glb" };
                ExportProgress = "Writing the 3D scene…";
                var report = await Task.Run(() => doc.ExportScene(cam, opt, path, fmt));
                var r = System.Text.Json.Nodes.JsonNode.Parse(report);
                ExportResult = $"Wrote {Path.GetFileName(path)} · {(int?)(double?)r?["spheres"] ?? 0:N0} atoms · {(int?)(double?)r?["cylinders"] ?? 0:N0} bonds" +
                               (_expFormat == 3 ? $" · render with POV-Ray: povray +W{w} +H{h} +A {Path.GetFileName(path)}" : $" · {(int?)(double?)r?["triangles"] ?? 0:N0} triangles");
                Status = ExportResult;
                return path;
            }
            else
            {
                var dpi = _expPreset >= 0 ? ExportPresets[_expPreset].Dpi : 0;
                var options = JsonSerializer.Serialize(new
                {
                    bits = _expFormat == 1 ? 16 : 8, dpi, colour_profile = _expProfile == 0 ? "srgb" : "none", provenance = _expProv, source = doc.Path,
                    engine = _expEngine == 1 ? "raytrace" : "raster", samples = _expTraceQ switch { 0 => 16, 2 => 256, _ => 64 },
                    shadows = _expShadows, occlusion = _expOcclusion, outlines = _expTraceOutlines,
                    aperture_fraction = _expDof switch { 1 => 0.012, 2 => 0.035, _ => 0.0 },
                });
                if (_expEngine == 1) ExportProgress = $"Ray tracing {w} × {h} · {(_expTraceQ switch { 0 => 16, 2 => 256, _ => 64 })} samples a pixel on every core…";
                var layer = _expLabels ? overlay(w, h) : null;
                await Task.Run(() => doc.ExportImage(cam, opt, path, options, layer));
            }
            var what = $"{w} × {h}" + (_expFormat == 1 ? " · 16-bit" : _expFormat == 2 ? " · SVG" : "") + (_expEngine == 1 && _expFormat != 2 ? " · ray traced" : "") + (_expProv && _expFormat != 2 ? " · provenance embedded" : "");
            ExportResult = $"Wrote {Path.GetFileName(path)} · {what}";
            Status = ExportResult;
            Record($"doc.render({PyStr(path)}, width={w}, height={h}, background=\"{ExportBackgrounds[_expBg].ToLowerInvariant()}\")");
            return path;
        }
        catch (Exception e) { ExportResult = "Export failed: " + e.Message; return null; }
        finally { ExportDlgBusy = false; ExportProgress = ""; }
    }

    private CancellationTokenSource? _movieCancel;
    public void CancelExport() => _movieCancel?.Cancel();

    /// <summary>Writes the movie: an animated PNG or a PNG sequence through the core; MP4 encodes a sequence with ffmpeg.</summary>
    public async Task<string?> ExportDialogMovie(string path)
    {
        if (_doc == null || _expBusy) return null;
        var doc = _doc;
        var (w, h) = ((int)_expW, (int)_expH);
        var opt = ExportDialogOptions(w, h, true);
        var cam = Camera;
        var mp4 = _movFormat == 2;
        var seqDir = mp4 ? Path.Combine(Path.GetTempPath(), "caps-movie-" + Guid.NewGuid().ToString("N")[..8]) : path;
        var options = JsonSerializer.Serialize(new
        {
            format = _movFormat == 0 ? "apng" : "png_sequence", step = _movStep, fps = _movFps,
            turntable = _movSource == 1 ? 360 : 0, turntable_frames = _movTurnFrames, colour_profile = _expProfile == 0 ? "srgb" : "none",
            provenance = _expProv ? "1" : "0", source = doc.Path,
        });
        _movieCancel = new CancellationTokenSource();
        var token = _movieCancel.Token;
        ExportDlgBusy = true;
        try
        {
            var n = await Task.Run(() => doc.ExportMovie(cam, opt, mp4 ? seqDir : path, options, (done, total) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => ExportProgress = $"Frame {done} of {total}");
                return !token.IsCancellationRequested;
            }));
            if (mp4)
            {
                ExportProgress = "Encoding MP4 with ffmpeg…";
                var psi = new System.Diagnostics.ProcessStartInfo(FfmpegPath!) { UseShellExecute = false, RedirectStandardError = true };
                foreach (var a in new[] { "-y", "-framerate", _movFps.ToString(CultureInfo.InvariantCulture), "-i", Path.Combine(seqDir, "frame_%05d.png"),
                                          "-c:v", "libx264", "-preset", "medium", "-crf", "16", "-tune", "animation", "-pix_fmt", "yuv420p",
                                          "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-movflags", "+faststart", path }) psi.ArgumentList.Add(a);
                using var p = System.Diagnostics.Process.Start(psi)!;
                var err = await p.StandardError.ReadToEndAsync();
                await p.WaitForExitAsync();
                try { Directory.Delete(seqDir, true); } catch { }
                if (p.ExitCode != 0) throw new InvalidOperationException("ffmpeg: " + err.Trim().Split('\n').LastOrDefault());
            }
            ExportResult = $"Wrote {Path.GetFileName(path)} · {n} frames · {w} × {h}" + (token.IsCancellationRequested ? " (stopped)" : "");
            Status = ExportResult;
            return path;
        }
        catch (Exception e) { ExportResult = "Export failed: " + e.Message; return null; }
        finally { ExportDlgBusy = false; ExportProgress = ""; _movieCancel = null; }
    }
}
